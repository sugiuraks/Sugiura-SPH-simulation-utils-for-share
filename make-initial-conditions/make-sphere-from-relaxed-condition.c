/////////////////////////////////////////////////////////////
//このプログラムはFDPS-basalt-box-R=50km-N=50000-after-relaxation-original.binなどの
//一様密度の箱の粒子配置を球に切り取るものです.
//球の半径ほか必要な情報はmain関数の頭で指定できます。
//
//使い方///////////////////////
//(Makefileはsrc, my-srcの場所やコンパイラなど適切に書き換えてください)
//このプログラムとMakefileがあるディレクトリで
//$ make makespherefromrelaxedcondition
//$ ./make-sphere-from-relaxed-condition.out
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

int main(int argc,char* argv[])
{
  //input parameters////////////////////////////////////////
  char input_filename[]="FDPS-basalt-box-R=50km-N=50000-cgs-with-maxu.bin"; //input particle information file
  char output_filename[]="FDPS-basalt-sphere-R=50km-N=50000-cgs-with-maxu.bin"; //output particle information file
  double time=0.0; //initial time
  double rho0=2.7; //density
  double R=5.0e6; //radius of output sphere
  //double R=5.0e6*pow(1.0/16.0,1.0/3.0); //もし質量比(ex.)1/16のインパクタを作りたければこちらを使ってください
  /////////////////////////////////////////////////////////

  int i,ix,iy,iz,i2,d;
  double x,y,z;
  double M=0.0;
  double Volume=(4.0/3.0)*M_PI*R*R*R;
  double init_D=0.0;
  
  PS::Initialize(argc,argv);
  PS::ParticleSystem<RealPtcl> sph_system_input,sph_system_output;
  sph_system_input.initialize();
  sph_system_output.initialize();
  FileHeader header_input,header_output;
  
  sph_system_input.readParticleBinary(input_filename, header_input);
  header_output.time = header_input.time;
  
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
      sph_system_output[i2].dens = rho0;
      sph_system_output[i2].D_cbrt = init_D;
      sph_system_output[i2].property_tag = 0;

      M += sph_system_output[i2].mass;
      i2++;
    }
  }

  //re-assign flaws
  flaw_assignment(sph_system_output,Volume);

  double km=1.0e5;
  FILE* fpout;
  char filenameout[]="test.txt";
  fpout=fopen(filenameout,"w");
  fprintf(fpout,"%f \n",0.0);
  fprintf(fpout,"%d \n",header_output.Nbody);
  for(i=0;i<header_output.Nbody;i++){
    fprintf(fpout,"%.2f %.2f %.2f %.2f %.2e %f %d \n",sph_system_output[i].pos.x/km,sph_system_output[i].pos.y/km,sph_system_output[i].pos.z/km,sph_system_output[i].smth/km,sph_system_output[i].mass,sph_system_output[i].dens,sph_system_output[i].id);
  }
  fclose(fpout);
  
  sph_system_output.writeParticleBinary(output_filename,header_output);
  printf("This program make basalt sphere : %s from \n",output_filename);
  printf("%s \n",input_filename);
  printf("total partice number: %d, total mass: %e \n",header_output.Nbody,M);
  
  PS::Finalize();
  return(0);
}
