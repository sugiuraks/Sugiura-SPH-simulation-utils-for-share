#include "SPH.h"
#include "time-development.h"
#include "space.h"
#include "input.h"
#include "mirror.h"
#include "param-mysorcecode.h"
#include <stdio.h>
#include <stdlib.h>

static int particle_number;
static unsigned long long mem=0;
static int max_flaw_number;

//read particle number,dimension and particle information from given filename's binary file and allocate particle array and return pointer to particle array.
//this program is for reading FDPS file
particle_t* read_particle_information(char filename[])
{
  FILE *fp;
  int d,i,d2;
  particle_t *par;
  double t0;
  double Range[D];
  get_Range(Range);
  double Rmin[D];
  get_Rmin(Rmin);
  double Sab_rho_temp[D*D];
  double Csmooth_init=Csmooth;
  int is_free_boundary;
  
  if((fp=fopen(filename,"rb"))==NULL){
    printf("cannot open particle number file! \n");
    exit(FAILURE);
  }
  //read initial time
  fread(&t0,sizeof(double),1,fp);
  set_t_as_t0(t0);
  //read particle number
  fread(&particle_number,sizeof(int),1,fp);
  //read max_flaw_number
  fread(&max_flaw_number,sizeof(int),1,fp);
  
  //allocation of particle array
  if(get_is_mirror()){
    par = ( particle_t* )malloc( sizeof( particle_t ) * ( particle_number + MAX_MIRROR_NUM ) );
    count_memory((int)(sizeof(particle_t)*(particle_number+MAX_MIRROR_NUM)));
  }
  else{
    par = ( particle_t* )malloc( sizeof( particle_t ) * ( particle_number) );
    count_memory((int)(sizeof(particle_t)*(particle_number)));
  }
  if(par == NULL){
    printf("cannot allocate particle array! \n");
    exit(FAILURE);
  }
  
  //set particle information 
  for(i=0 ; i<particle_number ; i++){
    if(
       fread(&par[i].r[X],sizeof(double),1,fp)          != 1 ||
       fread(&par[i].r[Y],sizeof(double),1,fp)          != 1 ||
       fread(&par[i].r[Z],sizeof(double),1,fp)          != 1 ||
       fread(&par[i].v[X],sizeof(double),1,fp)          != 1 ||
       fread(&par[i].v[Y],sizeof(double),1,fp)          != 1 ||
       fread(&par[i].v[Z],sizeof(double),1,fp)          != 1 ||
       fread(&par[i].mass,sizeof(double),1,fp)          != 1 ||
       fread(&par[i].u,sizeof(double),1,fp)             != 1 ||
       fread(&par[i].h,sizeof(double),1,fp)             != 1 ||
       fread(&par[i].rho_eos,sizeof(double),1,fp)       != 1 ||
       fread(&par[i].Sab_rho[X][X],sizeof(double),1,fp) != 1 ||
       fread(&par[i].Sab_rho[X][Y],sizeof(double),1,fp) != 1 ||
       fread(&par[i].Sab_rho[X][Z],sizeof(double),1,fp) != 1 ||
       fread(&par[i].Sab_rho[Y][X],sizeof(double),1,fp) != 1 ||
       fread(&par[i].Sab_rho[Y][Y],sizeof(double),1,fp) != 1 ||
       fread(&par[i].Sab_rho[Y][Z],sizeof(double),1,fp) != 1 ||
       fread(&par[i].Sab_rho[Z][X],sizeof(double),1,fp) != 1 ||
       fread(&par[i].Sab_rho[Z][Y],sizeof(double),1,fp) != 1 ||
       fread(&par[i].Sab_rho[Z][Z],sizeof(double),1,fp) != 1 ||
       fread(&par[i].property_tag,sizeof(int),1,fp)     != 1 ||
       fread(&par[i].ID,sizeof(int),1,fp)               != 1 ||
       fread(&par[i].D_cbrt,sizeof(double),1,fp)        != 1 ||
       fread(&par[i].alpha_por,sizeof(double),1,fp)     != 1 ||
       fread(&par[i].ni_tot_flaw,sizeof(long int),1,fp)      != 1 
       ){
      printf("file form error! \n");
      exit(FAILURE);
    }
    
    //allocation of array for activation threshold strain
    par[i].eps_ij_act_flaw = (double*)malloc(sizeof(double)*max_flaw_number);
    count_memory((int)(sizeof(double)*par[i].ni_tot_flaw));
    if(par[i].eps_ij_act_flaw == NULL){
      printf("can not allocate array for activation threshold strain! \n");
      exit(FAILURE);
    }
    fread(par[i].eps_ij_act_flaw,sizeof(double),max_flaw_number,fp);

	if( fread(&par[i].u_max,sizeof(double),1,fp)        != 1 ){
	  printf("file form error! \n");
      exit(FAILURE);
	}
    
    par[i].inter_par = (particle_t**)malloc(sizeof(particle_t*)*INTER);
    count_memory((int)(sizeof(particle_t*)*INTER));
    par[i].inter_Nmax = INTER;
    par[i].hmax = par[i].h;
    par[i].Csmooth = Csmooth_init;
    par[i].low_density_flag = 0;
    par[i].monoto_flag = 0;
  }
  
  is_free_boundary=1;
  for(d=0 ; d<D ; d++){
    if(boundary_flags[d]!=0) is_free_boundary=0;
  }
  if(is_free_boundary){
    reset_space_data(par);
  }
  else{
    for(i=0 ; i<particle_number ; i++){
      for(d=0 ; d<D ; d++){
	if( par[i].r[d] < Rmin[d] || par[i].r[d] > Rmin[d] + Range[d] ){
	  printf("par %d is out of range! dir = %d r[%d]=%f\n",i,d,d,par[i].r[d]);
	  exit(0);
	}
      }
    }
  }
  
  fclose(fp);
  
  return(par);
}

//free_particle_array
void free_particle_array(particle_t par[])
{
  int i;
  int particle_number=get_particle_number();
  
  for(i=0;i<particle_number;i++){
    free(par[i].inter_par);
    count_memory((int)(-(sizeof(particle_t*)*par[i].inter_Nmax)));
    if(is_fracture_model[par[i].property_tag]&&(par[i].eps_ij_act_flaw!=NULL)){
      free(par[i].eps_ij_act_flaw);
      count_memory((int)(-(sizeof(double)*par[i].ni_tot_flaw)));
    }
  }
  
  free(par);
  if(get_is_mirror()){
    count_memory((int)(sizeof(particle_t)*(particle_number+MAX_MIRROR_NUM)));
  }
  else{
    count_memory((int)(sizeof(particle_t)*(particle_number)));
  }
}

//get particle number
int get_particle_number(void)
{
  return(particle_number);
}

//get max_flaw_number
int get_max_flaw_number(void)
{
  return(max_flaw_number);
}

//count memory.
void count_memory(int count)
{
  mem += (unsigned long long)count;
}

//return memory
unsigned long long get_memory(void)
{
  return(mem);
}

