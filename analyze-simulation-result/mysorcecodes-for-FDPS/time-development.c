#include "SPH.h"
#include "input.h"
#include "mirror.h"
#include "param-mysorcecode.h"
#include "fracture-and-porosity-model.h"
#include "time-development.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static double dt;
static double t;
static int time_step_sign=1;
static int leapfrog_init=0;

//get dt
double get_dt(void)
{
  return(dt);
}

//set dt
void set_dt(particle_t par[])
{
  int i,d;
  int particle_number=get_particle_number();
  double nabla_dot_v;
  double temp;
  double g_i2;
  
  nabla_dot_v = 0.0;
  for(d=0 ; d<D ; d++){
    nabla_dot_v += par[0].dvdx[d][d];
  }
  dt = Ccfl * ( par[0].h / ( par[0].Cs + par[0].h * fabs(nabla_dot_v) ) );
  for(i=0 ; i<particle_number ; i++){
    nabla_dot_v = 0.0;
    for(d=0 ; d<D ; d++){
      nabla_dot_v += par[i].dvdx[d][d];
    }
    if( dt > ( temp = Ccfl * ( par[i].h / ( par[i].Cs + par[i].h * fabs(nabla_dot_v) ) ) ) ) dt = temp;
    if(is_selfgravity){
      g_i2 = 0.0;
      for(d=0 ; d<D ; d++){
	g_i2 += par[i].g[d] * par[i].g[d];
      }
      if( dt > ( temp = Ccfl * sqrt( par[i].h / sqrt(g_i2) ) ) ) dt = temp;
    }
  }
}

//calculate time development. This function only select time development method.
void time_development(particle_t par[])
{
  switch(time_development_method){
  case 0 : 
    time_development_for_Euler(par);
    break;
  case 1 :
    time_development_for_second_order_rungekutta(par);
    break;
  case 2 :
    if(leapfrog_init==0){
      initial_kick_for_leapfrog(par);
      leapfrog_init=1;
    }
    else{
      time_development_for_leapfrog(par);
    }
    break;
  }
}

//calculate time development of all particles and check internal energy not to be negative.
void time_development_for_Euler(particle_t par[])
{
  int i,d,d2;
  int particle_number=get_particle_number();
  int error=0;
  double a[D];

  calc_effective_strain(par);
  
#ifdef _OPENMP
#pragma omp parallel for private(i,d,a,d2)
#endif
  for(i=0 ; i<particle_number ; i++){

    if(is_fracture_model[par[i].property_tag]) time_development_of_D_i(&par[i],dt);

    for(d=0 ; d<D ; d++){
      if(is_selfgravity){
	a[d] = par[i].a[d] + par[i].g[d];
      }
      else{
	a[d] = par[i].a[d];
      }
      par[i].r[d] += par[i].v[d] * dt + 0.5 * a[d] * dt * dt;
      par[i].v[d] += a[d] * dt;
    }
    par[i].u += par[i].dudt * dt;
    for(d=0 ; d<D ; d++){
      for(d2=0 ; d2<D ; d2++){
	par[i].Sab_rho[d][d2] += par[i].dSab_rho_dt[d][d2] * dt;
      }
    }
    if(is_density_devel){
      par[i].rho_eos += par[i].drhodt * dt;
    }
    if(is_porosity_model[par[i].property_tag]){
      par[i].alpha_por += par[i].dalpha_dt * dt;
      if(par[i].alpha_por < 1.0) par[i].alpha_por = 1.0;
    }

    //copy velocity to old array for energy conservation//////
    for(d=0 ; d<D ; d++){
      par[i].v_old[d] = par[i].v[d];
    }
    //////////////////////////////////////////////////////////
    
    //if internal energy is negative
    if(par[i].u < ACC){
      //printf("par %d's u is negative! u=%f t=%f\n",i,par[i].u,t);
      par[i].u = ACC;
    }
    //if density becomes lower than critical value
    if(par[i].rho_eos < rho_min_crit && is_density_devel!=0){
      par[i].rho_eos = rho_min_crit;
      par[i].low_density_flag=1;
    }
    else{
      par[i].low_density_flag=0;
    }

    //calculate Sab from Sab/rho and rho
    /*for(d=0 ; d<D ; d++){
      for(d2=0 ; d2<D ; d2++){
	par[i].Sab[d][d2] = par[i].Sab_rho[d][d2] * par[i].rho;
      }
      }*/

    par[i].p_old = par[i].p;//hold pressure of previous time step for time development of distension parameter
  }
  
  //apply wall or periodic voundary condition for particle position
  relocate_particle_for_periodic_and_wall(par);
  
  t += dt;
}

//calculate time development of all particles for scond order rungekutta method
void time_development_for_second_order_rungekutta(particle_t par[])
{
  int i,d,d2;
  int particle_number=get_particle_number();
  int error=0;
  double a[D];

  calc_effective_strain(par);
  
#ifdef _OPENMP
#pragma omp parallel for private(i,d,a,d2)
#endif
  for(i=0 ; i<particle_number ; i++){

    if(is_fracture_model[par[i].property_tag]) time_development_of_D_i_for_second_order_rungekutta(&par[i],dt,time_step_sign);

    for(d=0 ; d<D ; d++){
      if(is_selfgravity){
	a[d] = par[i].a[d] + par[i].g[d];
      }
      else{
	a[d] = par[i].a[d];
      }
    }
    if(time_step_sign==1){
      for( d=0 ; d<D ; d++ ){
	par[i].r_old[d] = par[i].r[d];
	par[i].v_old[d] = par[i].v[d];
	par[i].r[d] += ( par[i].v[d] + 0.5 * a[d] * 0.5 * dt ) * 0.5 * dt;
	par[i].v[d] += a[d] * 0.5 * dt;
	for(d2=0 ; d2<D ; d2++){
	  par[i].Sab_rho_old[d][d2] = par[i].Sab_rho[d][d2];
	  par[i].Sab_rho[d][d2] += par[i].dSab_rho_dt[d][d2] * 0.5 * dt;
	}
      }
      par[i].u_old = par[i].u;
      par[i].u += par[i].dudt * 0.5 * dt;
      if(is_density_devel){
	par[i].rho_old = par[i].rho_eos;
	par[i].rho_eos += par[i].drhodt * 0.5 * dt;
      }
      if(is_porosity_model[par[i].property_tag]){
	par[i].alpha_por_old = par[i].alpha_por;
	par[i].alpha_por += par[i].dalpha_dt * 0.5 * dt;
	if(par[i].alpha_por < 1.0) par[i].alpha_por = 1.0;
	par[i].dalpha_dt_old = par[i].dalpha_dt;
      }
    }
    else{
      for( d=0 ; d<D ; d++ ){
	par[i].r[d] = par[i].r_old[d] + ( par[i].v_old[d] + 0.5 * dt * a[d] ) * dt;
	par[i].v[d] = par[i].v_old[d] + a[d] * dt;
	for( d2=0 ; d2<D ; d2++ ){
	  par[i].Sab_rho[d][d2] = par[i].Sab_rho_old[d][d2] + par[i].dSab_rho_dt[d][d2] * dt;
	}
      }
      par[i].u = par[i].u_old + par[i].dudt * dt;
      if(is_density_devel){
	par[i].rho_eos = par[i].rho_old + par[i].drhodt * dt;
      }
      if(is_porosity_model[par[i].property_tag]){
	par[i].alpha_por = par[i].alpha_por_old + par[i].dalpha_dt_old * dt;
	if(par[i].alpha_por < 1.0) par[i].alpha_por = 1.0;
      }

      //copy velocity to old array for energy conservation//////
      for(d=0 ; d<D ; d++){
	par[i].v_old[d] = par[i].v[d];
      }
      //////////////////////////////////////////////////////////
    }

    //if internal energy is negative
    if(par[i].u < ACC){
      //printf("par %d's u is negative! u=%f t=%f\n",i,par[i].u,t);
      par[i].u = ACC;
    }
    //if density becomes lower than critical value
    if(par[i].rho_eos < rho_min_crit && is_density_devel!=0){
      par[i].rho_eos = rho_min_crit;
      par[i].low_density_flag=1;
    }
    else{
      par[i].low_density_flag=0;
    }

    //calculate Sab from Sab/rho and rho
    /*for(d=0 ; d<D ; d++){
      for(d2=0 ; d2<D ; d2++){
	par[i].Sab[d][d2] = par[i].Sab_rho[d][d2] * par[i].rho;
      }
      }*/

    par[i].p_old = par[i].p;//hold pressure of previous time step for time development of distension parameter
  }
  
  //apply wall or periodic voundary condition for particle position
  relocate_particle_for_periodic_and_wall(par);
  
  if(time_step_sign==-1){
    t += dt;
  }
  //Change the odd or even time step
  time_step_sign *= -1;
}

//initial kick for leapfrogmethod. this function should be called before time roop begin.
void initial_kick_for_leapfrog(particle_t par[])
{
  int i,d,d2;
  int particle_number=get_particle_number();
  int error=0;
  double a[D];

  calc_effective_strain(par);
  
#ifdef _OPENMP
#pragma omp parallel for private(i,d,a,d2)
#endif
  for(i=0 ; i<particle_number ; i++){

    for(d=0 ; d<D ; d++){
      if(is_selfgravity){
	a[d] = par[i].a[d] + par[i].g[d];
      }
      else{
	a[d] = par[i].a[d];
      }
    }
    //store the information of initial timestep (v,u,rho,Sab/rho)
    for( d=0 ; d<D ; d++ ){
      par[i].v_old[d] = par[i].v[d];
      for(d2=0 ; d2<D ; d2++){
	par[i].Sab_rho_old[d][d2] = par[i].Sab_rho[d][d2];
      }
    }
    par[i].u_old = par[i].u;
    if(is_density_devel){
      par[i].rho_old = par[i].rho_eos;
    }
    if(is_fracture_model[par[i].property_tag]){
      par[i].D_cbrt_old = par[i].D_cbrt;
    }
    if(is_porosity_model[par[i].property_tag]){
      par[i].alpha_por_old = par[i].alpha_por;
    }

    //initial kick
    for( d=0 ; d<D ; d++ ){
      par[i].r[d] += ( par[i].v[d] + 0.5 * a[d] * dt ) * dt;
      par[i].v[d] = par[i].v_old[d] + a[d] * dt;
      for( d2=0 ; d2<D ; d2++ ){
	par[i].Sab_rho[d][d2] = par[i].Sab_rho_old[d][d2] + par[i].dSab_rho_dt[d][d2] * dt;
      }
    }
    par[i].u = par[i].u_old + par[i].dudt * dt;
    if(is_density_devel){
      par[i].rho_eos = par[i].rho_old + par[i].drhodt * dt;
    }
    if(is_porosity_model[par[i].property_tag]){
      par[i].alpha_por = par[i].alpha_por_old + par[i].dalpha_dt * dt;
      if(par[i].alpha_por<1.0) par[i].alpha_por=1.0;
    }
    if(is_fracture_model[par[i].property_tag]) time_development_of_D_i_for_leapfrog(&par[i],dt);

    //store the time derivatives of initial timestep (a,dudt,drhodt,dSab_rho_dt)
    for( d=0 ; d<D ; d++){
      par[i].a_old[d] = par[i].a[d];
      par[i].g_old[d] = par[i].g[d];
      for( d2=0 ; d2<D ; d2++){
	par[i].dSab_rho_dt_old[d][d2] = par[i].dSab_rho_dt[d][d2];
      }
    }
    par[i].dudt_old = par[i].dudt;
    if(is_density_devel){
      par[i].drhodt_old = par[i].drhodt;
    }
    if(is_fracture_model[par[i].property_tag]) par[i].dD_cbrt_dt_old = 0.0;
    if(is_porosity_model[par[i].property_tag]){
      par[i].dalpha_dt_old = par[i].dalpha_dt;
    }

    //if internal energy is negative
    if(par[i].u < ACC){
      //printf("par %d's u is negative! u=%f t=%f\n",i,par[i].u,t);
      par[i].u = ACC;
    }
    
    //calculate Sab from Sab/rho and rho
    /*for(d=0 ; d<D ; d++){
      for(d2=0 ; d2<D ; d2++){
	par[i].Sab[d][d2] = par[i].Sab_rho[d][d2] * par[i].rho;
      }
      }*/

    par[i].p_old = par[i].p;//hold pressure of previous time step for time development of distension parameter
  }
  //apply wall or periodic voundary condition for particle position
  relocate_particle_for_periodic_and_wall(par);
  
}

//time development for leapfrog method
void time_development_for_leapfrog(particle_t par[])
{
  int i,d,d2;
  int particle_number=get_particle_number();
  int error=0;
  double a[D],a_old[D];
  double dudt_old,dudt_new;
  double t_center_vel[D];

  calc_effective_strain(par);
  
#ifdef _OPENMP
#pragma omp parallel for private(i,d,a,a_old,d2,dudt_old,dudt_new,t_center_vel)
#endif
  for(i=0 ; i<particle_number ; i++){

    if(is_fracture_model[par[i].property_tag]) time_development_of_D_i_for_leapfrog(&par[i],dt);

    for(d=0 ; d<D ; d++){
      if(is_selfgravity){
	a[d] = par[i].a[d] + par[i].g[d];
	a_old[d] = par[i].a_old[d] + par[i].g_old[d];
      }
      else{
	a[d] = par[i].a[d];
	a_old[d] = par[i].a_old[d];
      }
    }

    //calculate time centered velocity for total energy conservation
    for( d=0 ; d<D ; d++){
      t_center_vel[d] = par[i].v_old[d] + 0.5 * 0.5 * ( a_old[d] + a[d] ) * dt;
    }

    //time development
    for( d=0 ; d<D ; d++ ){
      par[i].v_old[d] += 0.5 * ( a_old[d] + a[d] ) * dt;
      for(d2=0 ; d2<D ; d2++){
	par[i].Sab_rho_old[d][d2] += 0.5 * ( par[i].dSab_rho_dt_old[d][d2] + par[i].dSab_rho_dt[d][d2] ) * dt;
      }
    }
    par[i].u_old += 0.5 * ( par[i].dudt_old + par[i].dudt ) * dt;
    if(is_density_devel){
      par[i].rho_old += 0.5 * ( par[i].drhodt_old + par[i].drhodt ) * dt;
    }
    if(is_porosity_model[par[i].property_tag]){
      par[i].alpha_por_old += 0.5 * (par[i].dalpha_dt_old + par[i].dalpha_dt ) * dt;
      if(par[i].alpha_por_old < 1.0) par[i].alpha_por_old=1.0;
    }

    //calculate variables for calculation of time derivatives
    for( d=0 ; d<D ; d++ ){
      par[i].r[d] += ( par[i].v_old[d] + 0.5 * a[d] * dt ) * dt;
      par[i].v[d] = par[i].v_old[d] + a[d] * dt;
      for( d2=0 ; d2<D ; d2++ ){
	par[i].Sab_rho[d][d2] = par[i].Sab_rho_old[d][d2] + par[i].dSab_rho_dt[d][d2] * dt;
      }
    }
    par[i].u = par[i].u_old + par[i].dudt * dt;
    if(is_density_devel){
      par[i].rho_eos = par[i].rho_old + par[i].drhodt * dt;
    }
    if(is_porosity_model[par[i].property_tag]){
      par[i].alpha_por = par[i].alpha_por/*_old*/ + par[i].dalpha_dt * dt;
      if(par[i].alpha_por < 1.0) par[i].alpha_por=1.0;
    }

    //store the variables of previous timestep (a,dudt,drhodt,dSab_rho_dt)
    for( d=0 ; d<D ; d++){
      par[i].a_old[d] = par[i].a[d];
      par[i].g_old[d] = par[i].g[d];
      for( d2=0 ; d2<D ; d2++){
	par[i].dSab_rho_dt_old[d][d2] = par[i].dSab_rho_dt[d][d2];
      }
    }
    par[i].dudt_old = par[i].dudt;
    if(is_density_devel){
      par[i].drhodt_old = par[i].drhodt;
    }
    if(is_porosity_model[par[i].property_tag]){
      par[i].dalpha_dt_old = par[i].dalpha_dt;
    }

    //if internal energy is negative
    if(par[i].u < ACC){
      //printf("par %d's u is negative! u=%f t=%f\n",i,par[i].u,t);
      par[i].u = ACC;
    }
    //if density becomes lower than critical value
    if(par[i].rho_eos < rho_min_crit && is_density_devel!=0){
      par[i].rho_eos = rho_min_crit;
      par[i].low_density_flag=1;
    }
    else{
      par[i].low_density_flag=0;
    }
    
    //calculate Sab from Sab/rho and rho
    /*for(d=0 ; d<D ; d++){
      for(d2=0 ; d2<D ; d2++){
	par[i].Sab[d][d2] = par[i].Sab_rho[d][d2] * par[i].rho;
      }
      }*/

    par[i].p_old = par[i].p;//hold pressure of previous time step for time development of distension parameter
  }

  //apply wall or periodic voundary condition for particle position
  relocate_particle_for_periodic_and_wall(par);
  
  t += dt;
}

//set t as t0
void set_t_as_t0(double t0)
{
  t = t0;
}

//get t
double get_t(void)
{
  return(t);
}

//set dt by argument of this function
void set_dt_as_const(double dt_const)
{
  dt = dt_const;
}
