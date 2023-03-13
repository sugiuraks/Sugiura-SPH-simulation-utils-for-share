#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include "SPH.h"
#include "param-mysorcecode.h"
#include "input.h"
#include "calc-pressure-force.h"
#include "time-development.h"
#include "make-initial-condition-file.h"
#include "tillotson.h"

//calculate effective strain for all particles. Before this function, we should call calc_plastic_Sab if we adopt plastic model.
void calc_effective_strain(particle_t par[])
{
  int i,d,d2;
  int particle_number=get_particle_number();

  double sigma_eff_i[D][D];
  double a,b,c,q,p;
  double _Complex alpha_plus,alpha_minus;
  double sigma_1,sigma_2,sigma_3;
  double sigma_max;

#ifdef _OPENMP
#pragma omp parallel for private(i,d,d2,sigma_eff_i,a,b,c,p,q,alpha_plus,alpha_minus,sigma_1,sigma_2,sigma_3,sigma_max)
#endif
  for(i=0;i<particle_number;i++){
    for(d=0 ; d<D ; d++){
      for(d2=0 ; d2<D ; d2++){
	sigma_eff_i[d][d2] = - ( (d == d2) ? 1.0 : 0.0 ) * par[i].p + par[i].Sab[d][d2];
      }
    }

#if D == 3
    a = - (sigma_eff_i[X][X] + sigma_eff_i[Y][Y] + sigma_eff_i[Z][Z]);
    b = sigma_eff_i[X][X]*sigma_eff_i[Y][Y] + sigma_eff_i[Y][Y]*sigma_eff_i[Z][Z] + sigma_eff_i[Z][Z]*sigma_eff_i[X][X] - sigma_eff_i[X][Y]*sigma_eff_i[X][Y] - sigma_eff_i[Y][Z]*sigma_eff_i[Y][Z] - sigma_eff_i[Z][X]*sigma_eff_i[Z][X];
    c = - (sigma_eff_i[X][X]*sigma_eff_i[Y][Y]*sigma_eff_i[Z][Z] + 2.0*sigma_eff_i[X][Y]*sigma_eff_i[Y][Z]*sigma_eff_i[Z][X] - sigma_eff_i[X][X]*sigma_eff_i[Y][Z]*sigma_eff_i[Y][Z] - sigma_eff_i[Y][Y]*sigma_eff_i[Z][X]*sigma_eff_i[Z][X] - sigma_eff_i[Z][Z]*sigma_eff_i[X][Y]*sigma_eff_i[X][Y]);
    q = ( 27.0*c + 2.0*a*a*a - 9.0*a*b ) / 54.0;
    p = ( 3.0*b - a*a ) / 9.0;
    alpha_plus = cpow( -q + csqrt( (double _Complex)(q*q + p*p*p) ) , 1.0/3.0 );
    alpha_minus = cpow( -q - csqrt( (double _Complex)(q*q + p*p*p) ) , 1.0/3.0 );
    sigma_1 = creal( alpha_plus + alpha_minus - (a/3.0) );
    sigma_2 = creal( 0.5 * ( -1.0 + I*sqrt(3.0) ) * alpha_plus + 0.5 * ( -1.0 - I*sqrt(3.0) ) * alpha_minus - (a/3.0) );
    sigma_3 = creal( 0.5 * ( -1.0 - I*sqrt(3.0) ) * alpha_plus + 0.5 * ( -1.0 + I*sqrt(3.0) ) * alpha_minus - (a/3.0) );
    sigma_max = ( sigma_1 > sigma_2 ? sigma_1 : sigma_2 );
    sigma_max = ( sigma_max > sigma_3 ? sigma_max : sigma_3 );
#endif 
#if D == 2
    sigma_max = 0.5 * ( sigma_eff_i[X][X] + sigma_eff_i[Y][Y] ) + sqrt( ( 0.5 * ( sigma_eff_i[X][X] - sigma_eff_i[Y][Y]) ) * ( 0.5 * ( sigma_eff_i[X][X] - sigma_eff_i[Y][Y]) ) + sigma_eff_i[X][Y]*sigma_eff_i[X][Y] );
#endif
#if D == 1
    sigma_max = sigma_eff_i[X][X];
#endif

    par[i].eps_i_eff = sigma_max / ( K_bulk[par[i].property_tag] + (4.0/3.0) * mu_shear[par[i].property_tag] );

  }

}

//time developmen of D for i-th particle
void time_development_of_D_i(particle_t* par_iP,double dt)
{
  double D_i_max,Cg,Rs;
  int ni=0;
  int j;
  double delta_D=0.01;
  double alpha_0=alpha_por_0[par_iP->property_tag];
  
  for(j=0 ; j<par_iP->ni_tot_flaw ; j++){
    if( par_iP->eps_ij_act_flaw[j] < par_iP->eps_i_eff ) ni++;
  }
  if(par_iP->ni_tot_flaw == 0){
    D_i_max = 0.0;
  }
  else{
    D_i_max = pow( ni / par_iP->ni_tot_flaw, 1.0/3.0 );
  }
  Cg = 0.4 * sqrt( par_iP->Cs*par_iP->Cs + 4.0*mu_shear[par_iP->property_tag]/(3.0*par_iP->rho) );
  Rs = 3.0 * par_iP->h;
  
  if(par_iP->damage >= D_i_max) par_iP->dD_cbrt_dt = 0.0;
  else par_iP->dD_cbrt_dt = ( Cg / Rs ) *ni;
  //add damage increasing due to pore compaction
  if(is_porosity_model[par_iP->property_tag]){
    par_iP->dD_cbrt_dt += - (1.0/3.0) * ( pow( 1.0 - (par_iP->alpha_por-1.0)/(alpha_0-1.0) + delta_D , -(2.0/3.0) ) / ( pow(1.0+delta_D,1.0/3.0) - pow(delta_D,1.0/3.0) ) ) * ( 1.0 / (alpha_0-1.0) ) * par_iP->dalpha_dt;
  }

  par_iP->D_cbrt += par_iP->dD_cbrt_dt * dt;

  par_iP->damage = pow(par_iP->D_cbrt,3.0);
  if(par_iP->damage >= 1.0){
    par_iP->damage = 1.0;
    par_iP->D_cbrt = 1.0;
  }
}

//time developmen of D for i-th particle
void time_development_of_D_i_for_second_order_rungekutta(particle_t* par_iP,double dt,int time_step_sign)
{
  double D_i_max,Cg,Rs;
  int ni=0;
  int j;
  double delta_D=0.01;
  double alpha_0=alpha_por_0[par_iP->property_tag];

  for(j=0 ; j<par_iP->ni_tot_flaw ; j++){
    if( par_iP->eps_ij_act_flaw[j] < par_iP->eps_i_eff ) ni++;
  }
  if(par_iP->ni_tot_flaw == 0){
    D_i_max = 0.0;
  }
  else{
    D_i_max = pow( ni / par_iP->ni_tot_flaw, 1.0/3.0 );
  }
  Cg = 0.4 * sqrt( par_iP->Cs*par_iP->Cs + 4.0*mu_shear[par_iP->property_tag]/(3.0*par_iP->rho) );
  Rs = 3.0 * par_iP->h;

  if(time_step_sign==1){
    par_iP->D_cbrt_old = par_iP->D_cbrt;

    if(par_iP->damage >= D_i_max) par_iP->dD_cbrt_dt = 0.0;
    else par_iP->dD_cbrt_dt = ( Cg / Rs ) * ni;
    //add damage increasing due to pore compaction
    if(is_porosity_model[par_iP->property_tag]){
      par_iP->dD_cbrt_dt += - (1.0/3.0) * ( pow( 1.0 - (par_iP->alpha_por-1.0)/(alpha_0-1.0) + delta_D , -(2.0/3.0) ) / ( pow(1.0+delta_D,1.0/3.0) - pow(delta_D,1.0/3.0) ) ) * ( 1.0 / (alpha_0-1.0) ) * par_iP->dalpha_dt;
    }
    par_iP->dD_cbrt_dt_old = par_iP->dD_cbrt_dt;

    par_iP->D_cbrt += par_iP->dD_cbrt_dt * 0.5 * dt;

    par_iP->damage = pow(par_iP->D_cbrt,3.0);
    if(par_iP->damage >= 1.0) par_iP->damage = 1.0;
  }
  else{
    if(par_iP->damage >= D_i_max) par_iP->dD_cbrt_dt = 0.0;
    else par_iP->dD_cbrt_dt = ( Cg / Rs ) * ni;
    //add damage increasing due to pore compaction
    if(is_porosity_model[par_iP->property_tag]){
      par_iP->dD_cbrt_dt += - (1.0/3.0) * ( pow( 1.0 - (par_iP->alpha_por_old-1.0)/(alpha_0-1.0) + delta_D , -(2.0/3.0) ) / ( pow(1.0+delta_D,1.0/3.0) - pow(delta_D,1.0/3.0) ) ) * ( 1.0 / (alpha_0-1.0) ) * par_iP->dalpha_dt_old;
    }
    
    par_iP->D_cbrt = par_iP->D_cbrt_old + par_iP->dD_cbrt_dt_old * dt;

    par_iP->damage = pow(par_iP->D_cbrt,3.0);
    if(par_iP->damage >= 1.0){
      par_iP->damage = 1.0;
      par_iP->D_cbrt = 1.0;
    }
  }
}

//time developmen of D for i-th particle
void time_development_of_D_i_for_leapfrog(particle_t* par_iP,double dt)
{
  double D_i_max,Cg,Rs;
  int ni=0;
  int j;
  double delta_D=0.01;
  double alpha_0=alpha_por_0[par_iP->property_tag];
  
  for(j=0 ; j<par_iP->ni_tot_flaw ; j++){
    if( par_iP->eps_ij_act_flaw[j] < par_iP->eps_i_eff ) ni++;
  }
  if(par_iP->ni_tot_flaw == 0){
    D_i_max = 0.0;
  }
  else{
    D_i_max = pow( ni / par_iP->ni_tot_flaw, 1.0/3.0 );
  }
  Cg = 0.4 * sqrt( par_iP->Cs*par_iP->Cs + 4.0*mu_shear[par_iP->property_tag]/(3.0*par_iP->rho) );
  Rs = 3.0 * par_iP->h;
  
  if(par_iP->damage >= D_i_max) par_iP->dD_cbrt_dt = 0.0;
  else par_iP->dD_cbrt_dt = ( Cg / Rs ) * ni;
  //add damage increasing due to pore compaction
  if(is_porosity_model[par_iP->property_tag]){
    par_iP->dD_cbrt_dt += - (1.0/3.0) * ( pow( 1.0 - (par_iP->alpha_por-1.0)/(alpha_0-1.0) + delta_D , -(2.0/3.0) ) / ( pow(1.0+delta_D,1.0/3.0) - pow(delta_D,1.0/3.0) ) ) * ( 1.0 / (alpha_0-1.0) ) * par_iP->dalpha_dt;
  }
  
  par_iP->D_cbrt_old += 0.5 * ( par_iP->dD_cbrt_dt_old + par_iP->dD_cbrt_dt ) * dt;

  par_iP->D_cbrt = par_iP->D_cbrt/*_old*/ + par_iP->dD_cbrt_dt * dt;

  par_iP->dD_cbrt_dt_old = par_iP->dD_cbrt_dt;

  par_iP->damage = pow(par_iP->D_cbrt,3.0);
  if(par_iP->damage >= 1.0){
    par_iP->damage = 1.0;
    par_iP->D_cbrt = 1.0;
  }

}

//calculate time derivative of distension parameter and modify time derivative of deviatoric stress tensor for porosity model.
void calc_dalpha_dt_and_modify_dSab_dt_for_porosity(particle_t par[])
{
  double dalpha_dp;
  double delp_delrho,delp_delu;
  double f,dalpha_drho;
  int i;
  int particle_number=get_particle_number();
  int d1,d2;
  int count=0;
  int error=0;

#ifdef _OPENMP
#pragma omp parallel for private(i,d1,d2,dalpha_dp,delp_delrho,delp_delu,f,dalpha_drho) shared(error) 
#endif
  for(i=0;i<particle_number;i++){
    if(is_porosity_model[par[i].property_tag]){
      //if pressure is smaller than elastic limit or unloading, time derivative of distension should be 0.
      if( par[i].p < Pe_por[par[i].property_tag] || (par[i].p - par[i].p_old) < 0.0 ){
	dalpha_dp = 0.0;
	par[i].dalpha_dt = 0.0;
	f = 1.0;
      }
      else{
	dalpha_dp = -2.0 * ( alpha_por_0[par[i].property_tag] - 1.0 ) * ( Ps_por[par[i].property_tag] - par[i].p ) / pow(Ps_por[par[i].property_tag]-Pe_por[par[i].property_tag],2.0);
	//calculate delp_delrho and delp_delu
	switch(equation_of_state[par[i].property_tag]){
	case 0 : calc_delp_delrho_and_delp_delu_ideal_gas(par[i].alpha_por*par[i].rho,par[i].u,par[i].property_tag,&delp_delrho,&delp_delu); break;
	case 1 : calc_delp_delrho_and_delp_delu_elastic_eos(par[i].alpha_por*par[i].rho,par[i].u,par[i].property_tag,&delp_delrho,&delp_delu); break;
	case 2 : calc_delp_delrho_and_delp_delu_stiffened_gas_EoS(par[i].alpha_por*par[i].rho,par[i].u,par[i].property_tag,&delp_delrho,&delp_delu); break;
	case 3 : calc_delp_delrho_and_delp_delu_tillotson(par[i].u,par[i].alpha_por*par[i].rho,par[i].property_tag,&delp_delrho,&delp_delu); break;
	case 4 : calc_delp_delrho_and_delp_delu_Mie_Gruneisen_EoS(par[i].alpha_por*par[i].rho,par[i].u,par[i].property_tag,&delp_delrho,&delp_delu); break;
	}
	
	par[i].dalpha_dt = ( ( par[i].dudt*delp_delu + par[i].alpha_por*par[i].drhodt*delp_delrho )/( par[i].alpha_por + dalpha_dp*(par[i].p - par[i].rho*delp_delrho) ) )*dalpha_dp;
	if(par[i].dalpha_dt>0.0){
	  par[i].dalpha_dt=0.0;
	  count++;
	}
	dalpha_drho = ( ( (par[i].p/(par[i].rho*par[i].rho))*delp_delu + par[i].alpha_por*delp_delrho )/( par[i].alpha_por + dalpha_dp*(par[i].p - par[i].rho*delp_delrho) ) )*dalpha_dp;
	f = 1.0 + dalpha_drho * ( par[i].rho / par[i].alpha_por );
      }
      
      //modify time derivative of deviatoric stress tensor for porosity model
      for(d1=0 ; d1<D ; d1++){
	for(d2=0 ; d2<D ; d2++){
	  par[i].dSab_rho_dt[d1][d2] = ( f / par[i].alpha_por ) * par[i].dSab_rho_dt[d1][d2] - ( 1.0 / par[i].alpha_por ) * ( par[i].Sab[d1][d2] / par[i].rho ) * par[i].dalpha_dt;
	  if( isnan(par[i].dSab_rho_dt[d1][d2]) ){
	    printf("par %d's dSabdt[%d][%d] is not a number! \n",i,d1,d2);
	    error=1;
	  }
	}
      }
    }
  }
  
  if(error==1){
    printf("there are not a numbers\n");
    if(!DEBUG){
      char filename[]="debug-particleinfo.bin";
      //make_initial_condition_file(particle_number,D,get_t(),par,filename);
    }
    exit(FAILURE);
  }

}
