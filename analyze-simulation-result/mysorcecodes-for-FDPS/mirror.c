#include "SPH.h"
#include "input.h"
#include "space.h"
#include "param-mysorcecode.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static int mirror_particle_number;
static double real_Range[D];
static double real_Rmin[D];
static double real_Rmax[D];
static double box_length[D][3]; //D direction's box size.If second argument is 1,this box is for real particles.If 0,this box locate in minas.If 2,in plus.
static int is_copied[MAX_MIRROR_NUM][D]; //if i particle is copied to D direction,this value is 1.
static int mirror_partner[MAX_MIRROR_NUM]; //element equals to mirror particle's number.This array have real particle' number that is copied to element number's mirror particle.

//read real particle space data and flags from param.h and set box length and examine where real particles are in real space.
void read_mirror_datas(particle_t par[])
{
  int d,i;
  int particle_number=get_particle_number();
  double Rmin[D];
  get_Rmin(Rmin);
  double Range[D];
  get_Range(Range);
  
  //read real particle's Range
  for(d=0 ; d<D ; d++){
    real_Range[d] = init_real_Range[d];
  }
  
  //read real particle's Rmin
  for(d=0 ; d<D ; d++){
    real_Rmin[d] = init_real_Rmin[d];
  }
  
  //real particle's space should be smaller than all space
  for(d=0 ; d<D ; d++){
    if(Rmin[d] > real_Rmin[d]){
      printf("Rmin[%d] is larger than real particle's \n",d);
      exit(FAILURE);
    }
    if(Rmin[d] + Range[d] < real_Rmin[d] + real_Range[d]){
      printf("Rmax[%d] is smaller than real particle's \n",d);
      exit(FAILURE);
    }
    
    //set real_Rmax and box length
    real_Rmax[d] = real_Rmin[d] + real_Range[d];
    box_length[d][0] = real_Rmin[d] - Rmin[d];
    box_length[d][1] = real_Range[d];
    box_length[d][2] = Rmin[d] + Range[d] - real_Rmax[d];
  }
  
  //examin that real particles are really in real space
  for(i=0 ; i<particle_number ; i++){
    for(d=0 ; d<D ; d++){
      if(par[i].r[d] < real_Rmin[d] || par[i].r[d] > real_Rmax[d]){
	printf("particle is out of real space! par %d dir %d r[dir]=%f\n",i,d,par[i].r[d]);
	exit(FAILURE);
      }
    }
  }
}

//get mirror particle number
int get_mirror_particle_number(void)
{
  return(mirror_particle_number);
}

//return 1 if this particle should be copied to a direction
int mirror_condition(particle_t par_i,int dir,int box_num)
{
  int ret;
  if(box_num == 1 || boundary_flags[dir] == 0){
    ret = 0;
  }
  if(box_num == 0){
    if(boundary_flags[dir] == 1){
      if(real_Rmin[dir] + box_length[dir][0] > par_i.r[dir]) ret = 1;
      else ret = 0;
    }
    if(boundary_flags[dir] == 2){
      if(real_Rmax[dir] - box_length[dir][0] < par_i.r[dir]) ret = 1;
      else ret = 0;
    }
  }
  if(box_num == 2){
    if(boundary_flags[dir] == 1){
      if(real_Rmax[dir] - box_length[dir][2] < par_i.r[dir]) ret = 1;
      else ret = 0;
    }
    if(boundary_flags[dir] == 2){
      if(real_Rmin[dir] + box_length[dir][2] > par_i.r[dir]) ret = 1;
      else ret = 0;
    }
  }
  return(ret);
}

//return delta of real particle's positon and mirror particle's position of the direction
double mirror_pos_delta(particle_t par_i,int dir,int box_num)
{
  double ret;
  if(box_num == 1 || boundary_flags[dir] == 0){
    ret = 0.0;
  }
  if(box_num == 0){
    if(boundary_flags[dir] == 1){
      ret = -2 * fabs(par_i.r[dir] - real_Rmin[dir]);
    }
    if(boundary_flags[dir] == 2){
      ret = -real_Range[dir];
    }
  }
  if(box_num == 2){
    if(boundary_flags[dir] == 1){
      ret = 2 * fabs(real_Rmax[dir] - par_i.r[dir]);
    }
    if(boundary_flags[dir] == 2){
      ret = real_Range[dir];
    }
  }
  
  return(ret);
}

//locate mirror particles
void locate_mirror_particle(particle_t par[])
{
  int i,d,Cx,Cy,Cz;
  int particle_number=get_particle_number();

  mirror_particle_number = 0;
  for(i=0 ; i<particle_number ; i++){
    
    for(Cx=0 ; Cx<=2 ; Cx++){
    if(mirror_condition(par[i],X,Cx) || Cx==1){
      if(D == 1){
	if(Cx != 1){
	  par[particle_number+mirror_particle_number].r[X] = par[i].r[X] + mirror_pos_delta(par[i],X,Cx);
	  
	  if(Cx != 1 && boundary_flags[X] == 1) par[particle_number+mirror_particle_number].v[X] = -par[i].v[X];
	  else                                  par[particle_number+mirror_particle_number].v[X] = par[i].v[X];
	  
	  par[particle_number+mirror_particle_number].mass = par[i].mass;
	  par[particle_number+mirror_particle_number].h = par[i].h;
	  par[particle_number+mirror_particle_number].ID = particle_number+mirror_particle_number;
	  
	  mirror_partner[mirror_particle_number] = i;
	  
	  is_copied[mirror_particle_number][X] = 1;
	  
	  mirror_particle_number++;
	  if(mirror_particle_number>=MAX_MIRROR_NUM){
	    printf("too much mirror particles! \n");
	    exit(FAILURE);
	  }
	}
	continue;
      }
      for(Cy=0 ; Cy<=2 ; Cy++){
      if(mirror_condition(par[i],Y,Cy) || Cy==1){
	if(D == 2){
	  if(Cx != 1 || Cy != 1){
	    par[particle_number+mirror_particle_number].r[X] = par[i].r[X] + mirror_pos_delta(par[i],X,Cx);
	    par[particle_number+mirror_particle_number].r[Y] = par[i].r[Y] + mirror_pos_delta(par[i],Y,Cy);
	    
	    if(Cx != 1 && boundary_flags[X] == 1) par[particle_number+mirror_particle_number].v[X] = -par[i].v[X];
	    else                                  par[particle_number+mirror_particle_number].v[X] = par[i].v[X];
	    if(Cy != 1 && boundary_flags[Y] == 1) par[particle_number+mirror_particle_number].v[Y] = -par[i].v[Y];
	    else                                  par[particle_number+mirror_particle_number].v[Y] = par[i].v[Y];
	    
	    par[particle_number+mirror_particle_number].mass = par[i].mass;
	    par[particle_number+mirror_particle_number].h = par[i].h;
	    par[particle_number+mirror_particle_number].ID = particle_number+mirror_particle_number;
	    
	    mirror_partner[mirror_particle_number] = i;
	    
	    if(Cx != 1) is_copied[mirror_particle_number][X] = 1;
	    else        is_copied[mirror_particle_number][X] = 0;
	    if(Cy != 1) is_copied[mirror_particle_number][Y] = 1;
	    else        is_copied[mirror_particle_number][Y] = 0;
	    
	    mirror_particle_number++;
	    if(mirror_particle_number>=MAX_MIRROR_NUM){
	      printf("too much mirror particles! \n");
	      exit(FAILURE);
	    }
	  }
	  continue;
	}
	for(Cz=0 ; Cz<=2 ; Cz++){
	if(mirror_condition(par[i],Z,Cz) || Cz==1){
	  if(Cx != 1 || Cy != 1 || Cz != 1){
	    par[particle_number+mirror_particle_number].r[X] = par[i].r[X] + mirror_pos_delta(par[i],X,Cx);
	    par[particle_number+mirror_particle_number].r[Y] = par[i].r[Y] + mirror_pos_delta(par[i],Y,Cy);
	    par[particle_number+mirror_particle_number].r[Z] = par[i].r[Z] + mirror_pos_delta(par[i],Z,Cz);
	    
	    if(Cx != 1 && boundary_flags[X] == 1) par[particle_number+mirror_particle_number].v[X] = -par[i].v[X];
	    else                                  par[particle_number+mirror_particle_number].v[X] = par[i].v[X];
	    if(Cy != 1 && boundary_flags[Y] == 1) par[particle_number+mirror_particle_number].v[Y] = -par[i].v[Y];
	    else                                  par[particle_number+mirror_particle_number].v[Y] = par[i].v[Y];
	    if(Cy != 1 && boundary_flags[Z] == 1) par[particle_number+mirror_particle_number].v[Z] = -par[i].v[Z];
	    else                                  par[particle_number+mirror_particle_number].v[Z] = par[i].v[Z];
	    
	    par[particle_number+mirror_particle_number].mass = par[i].mass;
	    par[particle_number+mirror_particle_number].h = par[i].h;
	    par[particle_number+mirror_particle_number].ID = particle_number+mirror_particle_number;
	    
	    mirror_partner[mirror_particle_number] = i;
	    if(Cx != 1) is_copied[mirror_particle_number][X] = 1;
	    else        is_copied[mirror_particle_number][X] = 0;
	    if(Cy != 1) is_copied[mirror_particle_number][Y] = 1;
	    else        is_copied[mirror_particle_number][Y] = 0;
	    if(Cz != 1) is_copied[mirror_particle_number][Z] = 1;
	    else        is_copied[mirror_particle_number][Z] = 0;
	    
	    mirror_particle_number++;
	    if(mirror_particle_number>=MAX_MIRROR_NUM){
	      printf("too much mirror particles! \n");
	      exit(FAILURE);
	    }
	  }
	}
	}
      }
      }
    }
    }
  }
}

//copy various variables that need to be copied.
void copy_variables(particle_t par[])
{
  int mil_i;
  int particle_number=get_particle_number();
  int d1,d2;
  int sign;
  
#ifdef _OPENMP
#pragma omp parallel for private(mil_i,d1,d2,sign)
#endif
  for(mil_i=0 ; mil_i<mirror_particle_number ; mil_i++){
    par[particle_number+mil_i].rho = par[mirror_partner[mil_i]].rho;
    par[particle_number+mil_i].p = par[mirror_partner[mil_i]].p;
    par[particle_number+mil_i].Cs = par[mirror_partner[mil_i]].Cs;
    par[particle_number+mil_i].h = par[mirror_partner[mil_i]].h;
    par[particle_number+mil_i].gamma = par[mirror_partner[mil_i]].gamma;
    par[particle_number+mirror_particle_number].property_tag = par[mirror_partner[mil_i]].property_tag;
    par[particle_number+mil_i].monoto_flag = par[mirror_partner[mil_i]].monoto_flag;
    par[particle_number+mil_i].damage = par[mirror_partner[mil_i]].damage;
    par[particle_number+mil_i].low_density_flag = par[mirror_partner[mil_i]].low_density_flag;
    for( d1=0 ; d1<D ; d1++ ){
      for( d2=0 ; d2<D ; d2++ ){
	sign=1;
	if(is_copied[mil_i][d1] && boundary_flags[d1]==1) sign*=-1;
	if(is_copied[mil_i][d2] && boundary_flags[d2]==1) sign*=-1;
	par[particle_number+mil_i].Sab[d1][d2] = sign * par[mirror_partner[mil_i]].Sab[d1][d2];
      }
    }
  }
}
    
//copy gradients
void copy_gradients(particle_t par[])
{
  int mil_i,d,d2;
  int particle_number=get_particle_number();
  int sign;
  
#ifdef _OPENMP
#pragma omp parallel for private(mil_i,d,d2)
#endif
  for(mil_i=0 ; mil_i<mirror_particle_number ; mil_i++){
    for(d=0 ; d<D ; d++){
      if(is_copied[mil_i][d] && boundary_flags[d]==1){
	par[particle_number+mil_i].nablaV[d] = -par[mirror_partner[mil_i]].nablaV[d];
	par[particle_number+mil_i].nablap[d] = -par[mirror_partner[mil_i]].nablap[d];
	par[particle_number+mil_i].nablar[d] = -par[mirror_partner[mil_i]].nablar[d];
      }
      else{
	par[particle_number+mil_i].nablaV[d] = par[mirror_partner[mil_i]].nablaV[d];
	par[particle_number+mil_i].nablap[d] = par[mirror_partner[mil_i]].nablap[d];
	par[particle_number+mil_i].nablar[d] = par[mirror_partner[mil_i]].nablar[d];
      }
      for(d2=0 ; d2<D ; d2++){
	sign=1;
	if(is_copied[mil_i][d] && boundary_flags[d]==1) sign*=-1;
	if(is_copied[mil_i][d2] && boundary_flags[d2]==1) sign*=-1;
	par[particle_number+mil_i].dvdx[d][d2] = par[mirror_partner[mil_i]].dvdx[d][d2];
      }
    }
  }
}

//copy d2Vdx2
void copy_d2Vdx2(particle_t par[])
{
  int mil_i,d,d2;
  int particle_number=get_particle_number();
  int sign;
  
#ifdef _OPENMP
#pragma omp parallel for private(mil_i,d,d2,sign)
#endif
  for(mil_i=0 ; mil_i<mirror_particle_number ; mil_i++){
    for(d=0 ; d<D ; d++){
      for(d2=0 ; d2<D ; d2++){
	sign=1;
	if(is_copied[mil_i][d] && boundary_flags[d]==1) sign*=-1;
	if(is_copied[mil_i][d2] && boundary_flags[d2]==1) sign*=-1;
	par[particle_number+mil_i].d2Vdx2[d][d2] = sign*par[mirror_partner[mil_i]].d2Vdx2[d][d2];
      }
    }
  }
}

//if periodic or wall voundary and particle go out of real space for a direction,particle relocate for satisfying periodic or wall condition.
void relocate_particle_for_periodic_and_wall(particle_t par[])
{
  int i,d;
  int particle_number=get_particle_number();
  
#ifdef _OPENMP
#pragma omp parallel for private(i,d)
#endif
  for(i=0 ; i<particle_number ; i++){
    for(d=0 ; d<D ; d++){
      //for periodic voundary condition
      if(par[i].r[d] > real_Rmax[d] && boundary_flags[d] == 2){
	par[i].r[d] -= real_Range[d];
      }
      if(par[i].r[d] < real_Rmin[d] && boundary_flags[d] == 2){
	par[i].r[d] += real_Range[d];
      }
      //for wall voundaty condition
      if(par[i].r[d] > real_Rmax[d] && boundary_flags[d] == 1){
	par[i].r[d] = 2*real_Rmax[d] - par[i].r[d];
	par[i].v[d] *= -1.0;
      }
      if(par[i].r[d] < real_Rmin[d] && boundary_flags[d] ==1){
	par[i].r[d] = 2*real_Rmin[d] - par[i].r[d];
	par[i].v[d] *= -1.0;
      }
    }
  }
}
