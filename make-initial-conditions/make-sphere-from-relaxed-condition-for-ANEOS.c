/////////////////////////////////////////////////////////////
//このプログラムはFDPS-basalt-box-R=50km-N=50000-after-relaxation-original.binなどの
//一様密度の箱の粒子配置を球に切り取るものです.
//球の半径ほか必要な情報はmain関数の頭で指定できます。
//main関数の頭で与える温度を指定し, param.hで指定したANEOSテーブルに従って,
//圧力がなるべく0になるような密度と内部エネルギーを線形補間で求めて粒子に与えます.
//
//使い方///////////////////////
//(Makefileはsrc, my-srcの場所やコンパイラなど適切に書き換えてください)
//このプログラムとMakefileがあるディレクトリで
//$ make makespherefromrelaxedconditionANEOS
//$ ./make-sphere-from-relaxed-condition-for-ANEOS.out
/////////////////////////////////////////////////////////////

#include <particle_simulator.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "mathfunc.h"
#include "param.h"
#include "kernel.h"
#include "class.h"
#include "force.h"

#define D 3
#define X 0
#define Y 1
#define Z 2

static PS::F64 **temp_ane_table;
static PS::F64 **dens_ane_table;
static PS::F64 ***eng_ane_table;
static PS::F64 ***pres_ane_table;
static PS::F64 ***snds_ane_table;

//assignment of explicit flaws to particles
void flaw_assignment(PS::ParticleSystem<RealPtcl>& sph_system,double volume)
{
  long int i,j=0,i2;
  long int particle_number=sph_system.getNumberOfParticleLocal();
  long int** flaw_numbers;
  flaw_numbers=(long int**)malloc(particle_number*sizeof(long int*));//flaw numbers that are assigned to i-th particle
  int* flag_assigned;
  flag_assigned=(int*)malloc(particle_number*sizeof(int));//if any one flaw is assigned to i-th particle, i-th element of this array becomes 1
  int* number_of_flaw_max;
  number_of_flaw_max=(int*)malloc(particle_number*sizeof(int));//max number of flaws that are assigned to i-th particle
  int* number_of_flaw;
  number_of_flaw=(int*)malloc(particle_number*sizeof(int));//the number of flaws that are assigned to i-th particle
  int number_of_particle_assigned=0;//the number of particles that already have any one flaw
  int n;
  
  //plant the seed of random number
  srand(time(NULL));
  
  //initialize flaw conditions
  for(i=0;i<particle_number;i++){
    flag_assigned[i] = 0;
    number_of_flaw_max[i] = ((int)log((double)particle_number)+1);
    flaw_numbers[i] = (long int*)malloc(sizeof(long int)*number_of_flaw_max[i]);
    if(flaw_numbers[i]==NULL){
      printf("can not allocate flaw numbers array! \n");
      PS::Finalize();
      exit(1);
    }
    number_of_flaw[i] = 0;
    for(n=0;n<MAX_FLAW_NUMBER;n++){
      sph_system[i].eps_ij_act_flaw[n]=0.0;
    }
  }
  
  //flaw assignment
  while(number_of_particle_assigned < particle_number){
    i = (int)(((double)rand()/(double)RAND_MAX)*particle_number);
    flaw_numbers[i][number_of_flaw[i]] = j;
    number_of_flaw[i]++;
    if(flag_assigned[i]==0){
      number_of_particle_assigned++;
      flag_assigned[i]=1;
    }
    if(number_of_flaw[i] >= MAX_FLAW_NUMBER){
      printf("the number of flaws for one SPH particle exceeds defined max number! \n");
      PS::Finalize();
      exit(1);
    }
    if(number_of_flaw[i] >= number_of_flaw_max[i]){
      number_of_flaw_max[i]++;
      flaw_numbers[i] = (long int*)realloc(flaw_numbers[i],sizeof(long int)*number_of_flaw_max[i]);
      if(flaw_numbers[i]==NULL){
        printf("can not reallocate flaw numbers array! \n");
        PS::Finalize();
        exit(1);
      }
    }
    j++;
  }
  
  //flaw activation threshold strain assignment
  for(i=0;i<particle_number;i++){
    sph_system[i].ni_tot_flaw = number_of_flaw[i];
    for(i2=0;i2<number_of_flaw[i];i2++){
      sph_system[i].eps_ij_act_flaw[i2] = pow(flaw_numbers[i][i2]/(PARAM::K_WEIBULL[sph_system[i].property_tag]*volume),1.0/PARAM::M_WEIBULL[sph_system[i].property_tag]);
    }
  }
  
  //free flaw numbers array
  for(i=0;i<particle_number;i++){
    free(flaw_numbers[i]);
  }
  free(flaw_numbers);
  free(flag_assigned);
  free(number_of_flaw_max);
  free(number_of_flaw);

}

//read ANEOS table at initialization
void read_aneos_table(void)
{
  PS::S64 n,nd,nt;
  PS::S32 error=0;
  PS::S32 aneos_flag=0;
  char buf[256];
  PS::F64 temp1, temp2, temp3, temp4;
  
  FILE *fp_aneos_table;
  FILE *fp_temp_file;
  char temp_file_name[]="temp_aneos.txt";

  for(n=0 ; n<PARAM::N_MATERIAL ; n++){
    if(PARAM::EQUATION_OF_STATE[n]==5){
      aneos_flag = 1;
    }
  }
  if(aneos_flag == 1){
    //allocate arrays for storing aneos table////////////////////////////////////////////////
    temp_ane_table = (PS::F64**)malloc(sizeof(PS::F64*)*PARAM::N_MATERIAL);
    dens_ane_table = (PS::F64**)malloc(sizeof(PS::F64*)*PARAM::N_MATERIAL);
    eng_ane_table  = (PS::F64***)malloc(sizeof(PS::F64**)*PARAM::N_MATERIAL);
    pres_ane_table = (PS::F64***)malloc(sizeof(PS::F64**)*PARAM::N_MATERIAL);
    snds_ane_table = (PS::F64***)malloc(sizeof(PS::F64**)*PARAM::N_MATERIAL);
    if(temp_ane_table==NULL || dens_ane_table==NULL || eng_ane_table==NULL || pres_ane_table==NULL || snds_ane_table==NULL){
      error=1;
    }
    if(PS::Comm::getSum(error)>=1){
      if(PS::Comm::getRank()==0){
	std::cout << "cannot allocate aneos array for N_MATERIAL!" << std::endl;
      }
      PS::Finalize();
      exit(1);
    }

    for(n=0 ; n<PARAM::N_MATERIAL ; n++){
      if(PARAM::EQUATION_OF_STATE[n]==5){
	temp_ane_table[n] = (PS::F64*)malloc(sizeof(PS::F64)*PARAM::N_TEMP_GRID_ANE[n]);
	eng_ane_table[n]  = (PS::F64**)malloc(sizeof(PS::F64*)*PARAM::N_TEMP_GRID_ANE[n]);
	pres_ane_table[n] = (PS::F64**)malloc(sizeof(PS::F64*)*PARAM::N_TEMP_GRID_ANE[n]);
	snds_ane_table[n] = (PS::F64**)malloc(sizeof(PS::F64*)*PARAM::N_TEMP_GRID_ANE[n]);
	if(temp_ane_table[n]==NULL || eng_ane_table[n]==NULL || pres_ane_table[n]==NULL || snds_ane_table[n]==NULL){
	  error=1;
	}
	if(PS::Comm::getSum(error)>=1){
	  if(PS::Comm::getRank()==0){
	    std::cout << "cannot allocate aneos array for temperature!" << std::endl;
	  }
	  PS::Finalize();
	  exit(1);
	}

	dens_ane_table[n] = (PS::F64*)malloc(sizeof(PS::F64)*PARAM::N_DENS_GRID_ANE[n]);
	if(dens_ane_table[n]==NULL){
	  error=1;
	}
	if(PS::Comm::getSum(error)>=1){
	  if(PS::Comm::getRank()==0){
	    std::cout << "cannot allocate aneos array for density!" << std::endl;
	  }
	  PS::Finalize();
	  exit(1);
	}

	for(nt=0 ; nt<PARAM::N_TEMP_GRID_ANE[n] ; nt++){
	  eng_ane_table[n][nt]  = (PS::F64*)malloc(sizeof(PS::F64)*PARAM::N_DENS_GRID_ANE[n]);
	  pres_ane_table[n][nt] = (PS::F64*)malloc(sizeof(PS::F64)*PARAM::N_DENS_GRID_ANE[n]);
	  snds_ane_table[n][nt] = (PS::F64*)malloc(sizeof(PS::F64)*PARAM::N_DENS_GRID_ANE[n]);
	  if(eng_ane_table[n][nt]==NULL || pres_ane_table[n][nt]==NULL || snds_ane_table[n][nt]==NULL){
	    error=1;
	  }
	  if(PS::Comm::getSum(error)>=1){
	    if(PS::Comm::getRank()==0){
	      std::cout << "cannot allocate aneos array for density!" << std::endl;
	    }
	    PS::Finalize();
	    exit(1);
	  }
	}
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////

    //read aneos table//////////////////////////////////////////////////////////////////////
    for(n=0 ; n<PARAM::N_MATERIAL ; n++){
      //produce one-dimensional temporary aneos table//////////////////////
      if(PS::Comm::getRank()==0){
	if( (fp_aneos_table=fopen(PARAM::TABLE_NAME_ANE[n],"r"))==NULL ){
	  std::cout << "cannot open aneos table!" << std::endl;
	  PS::Abort(-1);
	}
	if( (fp_temp_file=fopen(temp_file_name,"w"))==NULL){
	  std::cout << "cannot open temporary file!" << std::endl;
	  PS::Abort(-1);
	}

	while(fgets(buf,sizeof(buf),fp_aneos_table)!=NULL){
	  if(sscanf(buf,"%le %le %le %le ",&temp1,&temp2,&temp3,&temp4)!=4){
	      std::cout << "aneos table format error!" << std::endl;
	      std::cout << "input aneos table should be 4 column table that does not include header information" << std::endl;
	      PS::Abort(-1);
	  }
	  fprintf(fp_temp_file,"%.12e \n",temp1);
	  fprintf(fp_temp_file,"%.12e \n",temp2);
	  fprintf(fp_temp_file,"%.12e \n",temp3);
	  fprintf(fp_temp_file,"%.12e \n",temp4);
	}
	fclose(fp_aneos_table);
	fclose(fp_temp_file);
      }
      PS::Comm::barrier();
      ////////////////////////////////////////////////////////////////////

      //read one-dimensional aneos table//////////////////////////////////
      if( (fp_temp_file=fopen(temp_file_name,"r"))==NULL){
	std::cout << "cannot open temporary file!" << std::endl;
	PS::Abort(-1);
      }

      //read temperature
      for(nt=0 ; nt<PARAM::N_TEMP_GRID_ANE[n] ; nt++){
	fgets(buf,sizeof(buf),fp_temp_file);
	sscanf(buf,"%le",&(temp_ane_table[n][nt]));
      }

      //read energy, pressure, and sound speed
      for(nt=0 ; nt<PARAM::N_TEMP_GRID_ANE[n] ; nt++){
	for(nd=0 ; nd<PARAM::N_DENS_GRID_ANE[n] ; nd++){
	  fgets(buf,sizeof(buf),fp_temp_file);
	  sscanf(buf,"%le",&(eng_ane_table[n][nt][nd]));
	  fgets(buf,sizeof(buf),fp_temp_file);
	  sscanf(buf,"%le",&(pres_ane_table[n][nt][nd]));
	  fgets(buf,sizeof(buf),fp_temp_file);
	  sscanf(buf,"%le",&(snds_ane_table[n][nt][nd]));
	}
      }

      //read density
      for(nd=0 ; nd<PARAM::N_DENS_GRID_ANE[n] ; nd++){
	fgets(buf,sizeof(buf),fp_temp_file);
	sscanf(buf,"%le",&(dens_ane_table[n][nd]));
      }
      
      fclose(fp_temp_file);
      ////////////////////////////////////////////////////////////////////
    }
    ////////////////////////////////////////////////////////////////////////////////////////
  }
}

//free ANEOS array at finalize
void free_aneos_array(void)
{
  PS::S64 n,nd,nt;
  PS::S32 aneos_flag=0;

  for(n=0 ; n<PARAM::N_MATERIAL ; n++){
    if(PARAM::EQUATION_OF_STATE[n]==5){
      aneos_flag = 1;
    }
  }
  if(aneos_flag == 1){
    for(n=0 ; n<PARAM::N_MATERIAL ; n++){
      if(PARAM::EQUATION_OF_STATE[n]==5){
	for(nt=0 ; nt<PARAM::N_TEMP_GRID_ANE[n] ; nt++){
	  free(eng_ane_table[n][nt]);
	  free(pres_ane_table[n][nt]);
	  free(snds_ane_table[n][nt]);
	}
	free(temp_ane_table[n]);
	free(eng_ane_table[n]);
	free(pres_ane_table[n]);
	free(snds_ane_table[n]);
	free(dens_ane_table[n]);
      } 
    }
    free(temp_ane_table);
    free(dens_ane_table);
    free(eng_ane_table);
    free(pres_ane_table);
    free(snds_ane_table);
  }
}

void calc_dens_and_eng_at_zero_pressure(double temp, int property_tag, double* densP, double* engP)
{
  int nt,nd;
  int nt_a,nt_b;
  int nd_a,nd_b;
  int flag_zero_pres=0;
  double* pres_at_temp;
  pres_at_temp = (double*)malloc(sizeof(double)*PARAM::N_DENS_GRID_ANE[property_tag]);
  if(pres_at_temp==NULL){
    printf("cannot allocate array for pres_at_temp! \n");
    exit(0);
  }
  double* eng_at_temp;
  eng_at_temp = (double*)malloc(sizeof(double)*PARAM::N_DENS_GRID_ANE[property_tag]);
  if(eng_at_temp==NULL){
    printf("cannot allocate array for eng_at_temp! \n");
    exit(0);
  }
  

  //find nt_a and nt_b that enclose temp
  if(temp < temp_ane_table[property_tag][0] || temp > temp_ane_table[property_tag][PARAM::N_TEMP_GRID_ANE[property_tag]-1] ){
    printf("temp:%.2e is out of ANEOS table! \n",temp);
    exit(0);
  }
  for(nt=0 ; nt<PARAM::N_TEMP_GRID_ANE[property_tag]-1 ; nt++){
    if(temp_ane_table[property_tag][nt] <= temp && temp < temp_ane_table[property_tag][nt+1]){
      nt_a = nt;
      nt_b = nt + 1;
    }
  }

  //calculate pressures and internal-energies at temp
  for(nd=0 ; nd<PARAM::N_DENS_GRID_ANE[property_tag] ; nd++){
    pres_at_temp[nd] = pres_ane_table[property_tag][nt_a][nd] + ( (pres_ane_table[property_tag][nt_b][nd] - pres_ane_table[property_tag][nt_a][nd]) / (temp_ane_table[property_tag][nt_b] - temp_ane_table[property_tag][nt_a]) ) * (temp - temp_ane_table[property_tag][nt_a]);
    eng_at_temp[nd] = eng_ane_table[property_tag][nt_a][nd] + ( (eng_ane_table[property_tag][nt_b][nd] - eng_ane_table[property_tag][nt_a][nd]) / (temp_ane_table[property_tag][nt_b] - temp_ane_table[property_tag][nt_a]) ) * (temp - temp_ane_table[property_tag][nt_a]);
  }

  //find nd_a and nd_b that enclose zero pressure
  for(nd=0 ; nd<PARAM::N_DENS_GRID_ANE[property_tag]-1 ; nd++){
    if(pres_at_temp[nd]<=0 && pres_at_temp[nd+1]>0){
      nd_a = nd;
      nd_b = nd + 1;
      flag_zero_pres = 1;
      break;
    }
  }
  if(flag_zero_pres==0){
    printf("cannot find zero-pressure density! \n");
    exit(0);
  }

  //calculate zero-pressure density and internal energy
  *densP = dens_ane_table[property_tag][nd_a] - ( (dens_ane_table[property_tag][nd_b] - dens_ane_table[property_tag][nd_a]) / (pres_at_temp[nd_b] - pres_at_temp[nd_a]) ) * pres_at_temp[nd_a];
  *engP = eng_at_temp[nd_a] - ( (eng_at_temp[nd_b] - eng_at_temp[nd_a]) / (pres_at_temp[nd_b] - pres_at_temp[nd_a]) ) * pres_at_temp[nd_a];

  free(pres_at_temp);
  free(eng_at_temp);
}

//calculate aneos pressure, sound speed, and temperature
void calc_p_Cs_and_temp_aneos(RealPtcl* par_i)
{
  PS::S64 nt,nd;
  PS::S64 nd_a,nd_b;
  PS::S64 nt_a,nt_b,nt_c,nt_d;
  PS::F64 dens_a,dens_b;
  PS::F64 eng_a,eng_b,eng_c,eng_d;
  PS::F64 pres_a,pres_b,pres_c,pres_d,pres_ac,pres_bd;
  PS::F64 snds_a,snds_b,snds_c,snds_d,snds_ac,snds_bd;
  PS::F64 temp_a,temp_b,temp_c,temp_d,temp_ac,temp_bd;

  PS::F64 dens_i,eng_i;
  if(PARAM::IS_POROSITY_MODEL[par_i->property_tag]){
    dens_i = par_i->alpha_por * par_i->dens;
  }
  else{
    dens_i = par_i->dens;
  }
  eng_i = par_i->eng;

  //determine density grids that enclose dens_i//////////////////////////////////////////////////////////////////////
  //if dens_i is out of ANEOS table, we use min/max value of density
  if(dens_i < dens_ane_table[par_i->property_tag][0]){
    printf("dens[%lld]=%.4e is out of ANEOS table! \n",par_i->id,dens_i);
    dens_i = dens_ane_table[par_i->property_tag][0];

    nd_a = 0;
    nd_b = 1;
  }
  else if(dens_i > dens_ane_table[par_i->property_tag][PARAM::N_DENS_GRID_ANE[par_i->property_tag]-1]){
    printf("dens[%lld]=%.4e is out of ANEOS table! \n",par_i->id,dens_i);
    dens_i = dens_ane_table[par_i->property_tag][PARAM::N_DENS_GRID_ANE[par_i->property_tag]-1];

    nd_a = PARAM::N_DENS_GRID_ANE[par_i->property_tag]-2;
    nd_b = PARAM::N_DENS_GRID_ANE[par_i->property_tag]-1;
  }
  else{
    for(nd=0 ; nd<PARAM::N_DENS_GRID_ANE[par_i->property_tag] ; nd++){
      if(dens_i >= dens_ane_table[par_i->property_tag][nd] && dens_i < dens_ane_table[par_i->property_tag][nd+1]){
	nd_a = nd;
	nd_b = nd+1;
      }
    }
  }
  dens_a = dens_ane_table[par_i->property_tag][nd_a];
  dens_b = dens_ane_table[par_i->property_tag][nd_b];
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //determine energy grids for nd_a (smaller density) that enclose eng_i/////////////////////////////////////////////
  //if eng_i is out of ANEOS table, we use min/max value of energy
  if(eng_i < eng_ane_table[par_i->property_tag][0][nd_a]){
    printf("eng[%lld]=%.4e is out of ANEOS table! \n",par_i->id,eng_i);
    eng_i = eng_ane_table[par_i->property_tag][0][nd_a];

    nt_a = 0;
    nt_c = 1;
  }
  else if(eng_i > eng_ane_table[par_i->property_tag][PARAM::N_TEMP_GRID_ANE[par_i->property_tag]-1][nd_a]){
    printf("eng[%lld]=%.4e is out of ANEOS table! \n",par_i->id,eng_i);
    eng_i = eng_ane_table[par_i->property_tag][PARAM::N_TEMP_GRID_ANE[par_i->property_tag]-1][nd_a];

    nt_a = PARAM::N_TEMP_GRID_ANE[par_i->property_tag]-2;
    nt_c = PARAM::N_TEMP_GRID_ANE[par_i->property_tag]-1;
  }
  else{
    for(nt=0 ; nt<PARAM::N_TEMP_GRID_ANE[par_i->property_tag] ; nt++){
      if(eng_i >= eng_ane_table[par_i->property_tag][nt][nd_a] && eng_i < eng_ane_table[par_i->property_tag][nt+1][nd_a]){
	nt_a = nt;
	nt_c = nt+1;
      }
    }
  }
  eng_a = eng_ane_table[par_i->property_tag][nt_a][nd_a];
  eng_c = eng_ane_table[par_i->property_tag][nt_c][nd_a];
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //determine energy grids for nd_b (larger density) that enclose eng_i/////////////////////////////////////////////
  //if eng_i is out of ANEOS table, we use min/max value of energy
  if(eng_i < eng_ane_table[par_i->property_tag][0][nd_b]){
    printf("eng[%lld]=%.4e is out of ANEOS table! \n",par_i->id,eng_i);
    eng_i = eng_ane_table[par_i->property_tag][0][nd_b];

    nt_b = 0;
    nt_d = 1;
  }
  else if(eng_i > eng_ane_table[par_i->property_tag][PARAM::N_TEMP_GRID_ANE[par_i->property_tag]-1][nd_b]){
    printf("eng[%lld]=%.4e is out of ANEOS table! \n",par_i->id,eng_i);
    eng_i = eng_ane_table[par_i->property_tag][PARAM::N_TEMP_GRID_ANE[par_i->property_tag]-1][nd_b];

    nt_b = PARAM::N_TEMP_GRID_ANE[par_i->property_tag]-2;
    nt_d = PARAM::N_TEMP_GRID_ANE[par_i->property_tag]-1;
  }
  else{
    for(nt=0 ; nt<PARAM::N_TEMP_GRID_ANE[par_i->property_tag] ; nt++){
      if(eng_i >= eng_ane_table[par_i->property_tag][nt][nd_b] && eng_i < eng_ane_table[par_i->property_tag][nt+1][nd_b]){
	nt_b = nt;
	nt_d = nt+1;
      }
    }
  }
  eng_b = eng_ane_table[par_i->property_tag][nt_b][nd_b];
  eng_d = eng_ane_table[par_i->property_tag][nt_d][nd_b];
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //calculate pressure/////////////////////////////////////////////////////////////////////////////////////////////
  pres_a = pres_ane_table[par_i->property_tag][nt_a][nd_a];
  pres_b = pres_ane_table[par_i->property_tag][nt_b][nd_b];
  pres_c = pres_ane_table[par_i->property_tag][nt_c][nd_a];
  pres_d = pres_ane_table[par_i->property_tag][nt_d][nd_b];

  pres_ac = pres_a + ( (pres_c - pres_a) / (eng_c - eng_a) ) * (eng_i - eng_a);
  pres_bd = pres_b + ( (pres_d - pres_b) / (eng_d - eng_b) ) * (eng_i - eng_b);

  par_i->pres = pres_ac + ( (pres_bd - pres_ac) / (dens_b - dens_a) ) * (dens_i - dens_a);
  if(PARAM::IS_POROSITY_MODEL[par_i->property_tag]){
    par_i->pres /= par_i->alpha_por;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //calculate sound speed//////////////////////////////////////////////////////////////////////////////////////////
  snds_a = snds_ane_table[par_i->property_tag][nt_a][nd_a];
  snds_b = snds_ane_table[par_i->property_tag][nt_b][nd_b];
  snds_c = snds_ane_table[par_i->property_tag][nt_c][nd_a];
  snds_d = snds_ane_table[par_i->property_tag][nt_d][nd_b];

  snds_ac = snds_a + ( (snds_c - snds_a) / (eng_c - eng_a) ) * (eng_i - eng_a);
  snds_bd = snds_b + ( (snds_d - snds_b) / (eng_d - eng_b) ) * (eng_i - eng_b);

  par_i->snds = snds_ac + ( (snds_bd - snds_ac) / (dens_b - dens_a) ) * (dens_i - dens_a);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //calculate temperature//////////////////////////////////////////////////////////////////////////////////////////
  temp_a = temp_ane_table[par_i->property_tag][nt_a];
  temp_b = temp_ane_table[par_i->property_tag][nt_b];
  temp_c = temp_ane_table[par_i->property_tag][nt_c];
  temp_d = temp_ane_table[par_i->property_tag][nt_d];

  temp_ac = temp_a + ( (temp_c - temp_a) / (eng_c - eng_a) ) * (eng_i - eng_a);
  temp_bd = temp_b + ( (temp_d - temp_b) / (eng_d - eng_b) ) * (eng_i - eng_b);

  par_i->temp = temp_ac + ( (temp_bd - temp_ac) / (dens_b - dens_a) ) * (dens_i - dens_a);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //we DO NOT calculate gamma for ANEOS table EoS,
  //so that it is dangerous to use Godunov SPH for the ANEOS EoS
  par_i->gamma = 1.0;

}


int main(int argc,char* argv[])
{
  //input parameters////////////////////////////////////////
  char input_filename[]="FDPS-basalt-box-R=50km-N=50000-after-relaxation-original-mks.bin"; //input particle information file
  char output_filename[]="FDPS-basalt-sphere-R=50km-N=50000-relaxed-basaltANEOS.bin"; //output particle information file
  //char output_filename[]="FDPS-basalt-sphere-Rt=50km-Nt=50000-q=1_50-relaxed-basaltANEOS.bin"; //output particle information file
  double time=0.0; //initial time
  double temp0=300.0; //temperature[K]
  double R=5.0e4; //radius of output sphere
  //double R=5.0e4*pow(1.0/50.0,1.0/3.0); //もし質量比(ex.)1/50のインパクタを作りたければこちらを使ってください
  int property_tag=0;
  /////////////////////////////////////////////////////////

  int i,ix,iy,iz,i2,d;
  double x,y,z;
  double M=0.0;
  double Volume=(4.0/3.0)*M_PI*R*R*R;
  double init_D=0.0;
  double dens,eng;
  
  PS::Initialize(argc,argv);
  PS::ParticleSystem<RealPtcl> sph_system_input,sph_system_output;
  sph_system_input.initialize();
  sph_system_output.initialize();
  FileHeader header_input,header_output;
  
  sph_system_input.readParticleBinary(input_filename, header_input);
  header_output.time = header_input.time;

  read_aneos_table();
  
  //calculate particle number for sphere
  header_output.Nbody=0;
  for(i=0;i<header_input.Nbody;i++){
    x=sph_system_input[i].pos.x;
    y=sph_system_input[i].pos.y;
    z=sph_system_input[i].pos.z;
    if(sqrt(x*x+y*y+z*z)<R) header_output.Nbody++;
  }
  sph_system_output.setNumberOfParticleLocal(header_output.Nbody);
  
  i2=0;
  //copy to real configuration
  for(i=0;i<header_input.Nbody;i++){
    x=sph_system_input[i].pos.x;
    y=sph_system_input[i].pos.y;
    z=sph_system_input[i].pos.z;
    if(sqrt(x*x+y*y+z*z)<R){
      sph_system_output[i2]=sph_system_input[i];
      sph_system_output[i2].id = i2;
      sph_system_output[i2].D_cbrt = init_D;
      sph_system_output[i2].property_tag = property_tag;
      sph_system_output[i2].Sab_rho = 0.0;

      sph_system_output[i2].temp = temp0;
      calc_dens_and_eng_at_zero_pressure(sph_system_output[i2].temp, sph_system_output[i2].property_tag, &dens, &eng);
      sph_system_output[i2].dens = dens;
      sph_system_output[i2].eng  = eng;
      sph_system_output[i2].mass = dens*pow(sph_system_output[i2].smth, 3.0);
      calc_p_Cs_and_temp_aneos(&sph_system_output[i2]);

      M += sph_system_output[i2].mass;
      i2++;
    }
  }

  //re-assign flaws
  flaw_assignment(sph_system_output,Volume);

  double km=1.0e3;
  FILE* fpout;
  char filenameout[]="test.txt";
  fpout=fopen(filenameout,"w");
  fprintf(fpout,"%f \n",0.0);
  fprintf(fpout,"%d \n",header_output.Nbody);
  for(i=0;i<header_output.Nbody;i++){
    fprintf(fpout,"%.2f %.2f %.2f %.2f %.2e %.6e %.6e %.6e %d \n",sph_system_output[i].pos.x/km,sph_system_output[i].pos.y/km,sph_system_output[i].pos.z/km,sph_system_output[i].smth/km,sph_system_output[i].mass,sph_system_output[i].dens,sph_system_output[i].eng,sph_system_output[i].pres,sph_system_output[i].id);
  }
  fclose(fpout);
  
  sph_system_output.writeParticleBinary(output_filename,header_output);
  printf("This program make basalt sphere : %s from \n",output_filename);
  printf("%s \n",input_filename);
  printf("total partice number: %d, total mass: %e \n",header_output.Nbody,M);

  free_aneos_array();
  PS::Finalize();
  return(0);
}
