#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "SPH.h"
#include "param-mysorcecode.h"
#include "input.h"

//make initial condition file for use of FDPS
//meaning of mode...if mode is 0, make initial condition file for all particles, but if 1, make initial condition file for particles that have remaining flag of 1.
void make_initial_condition_file(int particle_number,int dim,double time,particle_t par[],char filename[],int remaining_flag[],int mode)
{
  FILE *fp;
  int i,d,d2;
  int max_flaw_number=get_max_flaw_number();
  int output_particle_number=0;
  int new_ID=0;
  if(mode==0){
    output_particle_number=particle_number;
  }
  if(mode==1){
    for(i=0;i<particle_number;i++){
      if(remaining_flag[i]==1){
	output_particle_number++;
      }
    }
  }
  
  if((fp=fopen(filename,"wb"))==NULL){
    printf("cannot open file! \n");
    exit(0);
  }
  else{
    fwrite(&time,sizeof(double),1,fp);
    fwrite(&output_particle_number,sizeof(int),1,fp);
    fwrite(&max_flaw_number,sizeof(int),1,fp);
    for(i=0;i<particle_number;i++){
      if(mode==0||remaining_flag[i]==1){
	fwrite(&par[i].r[X],sizeof(double),1,fp);
	fwrite(&par[i].r[Y],sizeof(double),1,fp);
	fwrite(&par[i].r[Z],sizeof(double),1,fp);
	fwrite(&par[i].v[X],sizeof(double),1,fp);
	fwrite(&par[i].v[Y],sizeof(double),1,fp);
	fwrite(&par[i].v[Z],sizeof(double),1,fp);
	fwrite(&par[i].mass,sizeof(double),1,fp);
	fwrite(&par[i].u,sizeof(double),1,fp);
	fwrite(&par[i].h,sizeof(double),1,fp);
	fwrite(&par[i].rho_eos,sizeof(double),1,fp);
	fwrite(&par[i].Sab_rho[X][X],sizeof(double),1,fp);
	fwrite(&par[i].Sab_rho[X][Y],sizeof(double),1,fp);
	fwrite(&par[i].Sab_rho[X][Z],sizeof(double),1,fp);
	fwrite(&par[i].Sab_rho[Y][X],sizeof(double),1,fp);
	fwrite(&par[i].Sab_rho[Y][Y],sizeof(double),1,fp);
	fwrite(&par[i].Sab_rho[Y][Z],sizeof(double),1,fp);
	fwrite(&par[i].Sab_rho[Z][X],sizeof(double),1,fp);
	fwrite(&par[i].Sab_rho[Z][Y],sizeof(double),1,fp);
	fwrite(&par[i].Sab_rho[Z][Z],sizeof(double),1,fp);
	fwrite(&par[i].property_tag,sizeof(int),1,fp);
	fwrite(&new_ID,sizeof(int),1,fp);
	new_ID++;
	fwrite(&par[i].D_cbrt,sizeof(double),1,fp);
	fwrite(&par[i].alpha_por,sizeof(double),1,fp);
	fwrite(&par[i].ni_tot_flaw,sizeof(long int),1,fp);
	fwrite(par[i].eps_ij_act_flaw,sizeof(double),max_flaw_number,fp);
      }
    }
  }
  fclose(fp);
}

//assignment of explicit flaws to particles
void flaw_assignment(particle_t par[],int particle_number,double volume)
{
  int i,j=0,i2;
  int* flaw_numbers[particle_number];//flaw numbers that are assigned to i-th particle
  int flag_assigned[particle_number];//if any one flaw is assigned to i-th particle, i-th element of this array becomes 1
  int number_of_flaw_max[particle_number];//max number of flaws that are assigned to i-th particle
  int number_of_flaw[particle_number];//the number of flaws that are assigned to i-th particle
  int number_of_particle_assigned=0;//the number of particles that already have any one flaw

  //plant the seed of random number
  srand(time(NULL));

  //initialize flaw conditions
  for(i=0;i<particle_number;i++){
    flag_assigned[i] = 0;
    number_of_flaw_max[i] = ((int)log((double)particle_number)+1);
    flaw_numbers[i] = (int*)malloc(sizeof(int)*number_of_flaw_max[i]);
    if(flaw_numbers[i]==NULL){
      printf("can not allocate flaw numbers array! \n");
      exit(FAILURE);
    }
    number_of_flaw[i] = 0;
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
    if(number_of_flaw[i] >= number_of_flaw_max[i]){
      number_of_flaw_max[i]++;
      flaw_numbers[i] = (int*)realloc(flaw_numbers[i],sizeof(int)*number_of_flaw_max[i]);
      if(flaw_numbers[i]==NULL){
	printf("can not reallocate flaw numbers array! \n");
	exit(FAILURE);
      }
    }
    j++;
  }

  //flaw activation threshold strain assignment
  for(i=0;i<particle_number;i++){
    par[i].ni_tot_flaw = number_of_flaw[i];
    par[i].eps_ij_act_flaw = (double*)malloc(sizeof(double)*number_of_flaw[i]);
    if(par[i].eps_ij_act_flaw==NULL){
      printf("can not allocate threshold strain array! \n");
      exit(FAILURE);
    }
    for(i2=0;i2<number_of_flaw[i];i2++){
      par[i].eps_ij_act_flaw[i2] = pow(flaw_numbers[i][i2]/(k_weibull[par[i].property_tag]*volume),1.0/m_weibull[par[i].property_tag]);
    }
  }

  //free flaw numbers array
  for(i=0;i<particle_number;i++){
    free(flaw_numbers[i]);
  }

}
