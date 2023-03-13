#include "SPH.h"
#include "space.h"
#include "input.h"
#include "mirror.h"
#include "time-development.h"
#include "tillotson.h"
#include "calc-self-gravity.h"
#include "calc-pressure-force.h"
#include "param-mysorcecode.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <limits.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "make-initial-condition-file.h"

extern int time_step;

//fanctions for calculating hydro variables////////////////////////////////////////////////////////////////////////////////////////

//set h regulator as max of h
/*void set_h_regulator(particle_t par[])
{
int i;
int particle_number=get_particle_number();

h_regulator=par[0].h;
for(i=0;i<particle_number;i++){
if( h_regulator < par[i].h ) h_regulator = par[i].h;
}

h_regulator *= hmax_constant;
}*/

 //calculate smoothinglength 
void calc_h(particle_t *par)
{
  int i,d,n;
  int particle_number=get_particle_number();
  double delta2;
  double rhotmp;
  particle_t *par_j;
  
#ifdef _OPENMP
#pragma omp parallel for private(i,d,n,rhotmp,delta2,par_j)
#endif
  for(i=0;i<particle_number;i++){
    //set more smooth h and J
    par[i].h *= par[i].Csmooth;
    
    //calc density
    rhotmp=0.0;
    for(n=0 ; n<par[i].inter_N ; n++){
      par_j = par[i].inter_par[n];
      delta2=0.0;
      for(d=0;d<D;d++){
	delta2 += (par[i].r[d]-par_j->r[d]) * (par[i].r[d]-par_j->r[d]);
      }
      rhotmp += par_j->mass*pow(par[i].h*sqrt(PI),-D)*exp(-delta2/(par[i].h*par[i].h));
    }
    //set h 
    par[i].h = eta*pow(par[i].mass/rhotmp,1.0/D);
    //h regulation
    if(par[i].h > h_regulator){
      par[i].h = h_regulator;
    }
  }
}


//judge where particle's h overlap to the node or not.
int is_overlap(double par_pos[],double node_center[],double box_size[],double h_overlap,double par_i_Csmooth)
{
  double dif[D];
  int d;
  
  for( d=0 ; d<D ; d++){
    dif[d] = par_pos[d] - node_center[d];
    if( fabs(dif[d]) > 1.001 * ( par_i_Csmooth * Ncell * h_overlap + 0.5 * box_size[d] ) ) return(0);
  }
  
  return(1);
}

//calculate hmax recursively
void calc_hmax_recursive(particle_t *par_i,tree_node_t *n_node,double hmax_of_per_node[])
{
  int d,n;
  double h_overlap;
  if(n_node->level<2) h_overlap=hmax_of_per_node[n_node->level];
  else h_overlap=hmax_of_per_node[n_node->level-2];
  
  //if difference bitween i particle and this node is small
  if( is_overlap(par_i->r,n_node->center,n_node->box_size,h_overlap,par_i->Csmooth) ){
    //if this node level is bottom level
    if(n_node->is_leaf){
      //calc hmax
      if(par_i->hmax < n_node->hmax) par_i->hmax=n_node->hmax;
    }
    else{
      for(n=0;n<N_CHILD;n++){
	//if this node is not reaf node,call same fanction recursively to valid child.
	if(n_node->child[n]!=NULL){
	  calc_hmax_recursive(par_i,n_node->child[n],hmax_of_per_node);
	}
      }
    }
  }
}

//calculate hmax 
void calc_hmax(particle_t par[],tree_node_t *root_node)
{
  int i,d,n;
  int particle_number=get_particle_number();
  double hmax_of_per_node[MAX_TREE_LEVEL]; //hmax of node that have particle i.Array number equals to that node level.
  tree_node_t *n_node;
  
  //particle roop
#ifdef _OPENMP
#pragma omp parallel for private(i,d,n_node,hmax_of_per_node)
#endif
  for(i=0;i<particle_number;i++){
    //set hmax of node that have particle i to hmax_of_per_node array
    n_node=par[i].tree_parent;
    do{
      hmax_of_per_node[n_node->level]=n_node->hmax;
      n_node=n_node->parent;
    }while(n_node);
    
    //calculate hmax recursively
    par[i].hmax=par[i].h;
    calc_hmax_recursive(&par[i],root_node,hmax_of_per_node);
  }
}

//calculate density or time derivative of density
void calc_density_or_drhodt(particle_t par[])
{
  int i,d,n;
  int particle_number=get_particle_number();
  particle_t *par_j;
  double dif[D],eij[D];
  double delta2,hij,delta;
  double dif_v[D];
  double dWdx[D];
  int porosity_flag=0;
  
  //if density is calculated from summation
  if(is_density_devel==0||is_density_devel==1){
#ifdef _OPENMP
#pragma omp parallel for private(i,n,d,delta2,par_j)
#endif
    for( i=0 ; i<particle_number ; i++ ){
      par[i].rho=0.0;
      for( n=0 ; n<par[i].inter_N ; n++ ){
	par_j = par[i].inter_par[n];
	delta2=0.0;
	for(d=0;d<D;d++){
	  delta2 += (par[i].r[d]-par_j->r[d])*(par[i].r[d]-par_j->r[d]);
	}
	par[i].rho += par_j->mass*pow(par[i].h*sqrt(PI),-D)*exp(-delta2/(par[i].h*par[i].h));
      }
      if(is_density_devel==0){
	par[i].rho_eos = par[i].rho;
      }
      if(is_porosity_model[par[i].property_tag]){
	porosity_flag=1;
      }
    }
  }
  //if density is time developed
  if(is_density_devel==1||is_density_devel==2||porosity_flag==1){
#ifdef _OPENMP
#pragma omp parallel for private(i,n,d,delta2,par_j,dif,eij,delta,hij,dif_v,dWdx)
#endif
    for( i=0 ; i<particle_number ; i++ ){
      par[i].drhodt = 0.0;
      for( n=0 ; n<par[i].inter_N ; n++ ){
	par_j = par[i].inter_par[n];
	delta2 = 0.0;
	for( d=0 ; d<D ; d++ ){
	  dif[d] = par[i].r[d] - par_j->r[d];
	  delta2 += dif[d] * dif[d];
	  dif_v[d] = par[i].v[d] - par_j->v[d];
	}
	
	nablaW_gauss(delta2,dif,par[i].h,dWdx);
	
	for(d=0;d<D;d++){
	  par[i].drhodt += /*( par[i].rho / par_j->rho ) */ par_j->mass * dif_v[d] * dWdx[d];
	}
      }
      if(is_density_devel==2){
	par[i].rho = par[i].rho_eos;
      }
    }
  }
}

//calculate pressure, sound speed and gamma (for Riemann solver of ideal gas). This function only select appropriate EoS for each particle.
void calc_pressures(particle_t par[])
{
  int i;
  int particle_number=get_particle_number();

#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
  for(i=0;i<particle_number;i++){
    switch(equation_of_state[par[i].property_tag]){
    case 0 : calc_p_and_Cs_ideal_gas(&par[i]); break;
    case 1 : calc_p_elastic_eos(&par[i]); break;
    case 2 : calc_p_Cs_and_gamma_stiffened_gas_EoS(&par[i]); break;
    case 3 : calc_p_Cs_and_gamma_tillotson(&par[i]); break;
    case 4 : calc_p_Cs_and_gamma_Mie_Gruneisen_EoS(&par[i]); break;
    }
    if(par[i].p < 0.0) par[i].p *= ( 1.0 - par[i].damage);//In this simulation, we throughly use reduced pressure. (2016/4/13)
  }
}

//calculate pressure,sound speed using ideal gas eos.
void calc_p_and_Cs_ideal_gas(particle_t* par_i)
{
  par_i->p = ( gamma_ideal[par_i->property_tag] - 1.0 ) * par_i->rho_eos * par_i->u;
  par_i->gamma = gamma_ideal[par_i->property_tag];
  par_i->Cs = sqrt( gamma_ideal[par_i->property_tag] * par_i->p  / par_i->rho_eos );
}

//calculate delp_delrho and delp_delu for ideal gas
void calc_delp_delrho_and_delp_delu_ideal_gas(double rho,double u,int p_tag,double* delp_delrho_P,double* delp_delu_P)
{
  *delp_delrho_P = ( gamma_ideal[p_tag] - 1.0 ) * u;
  *delp_delu_P = ( gamma_ideal[p_tag] - 1.0 ) * rho;
}

//calculate pressure,sound speed and specific heat ratio from internal energy and density.
void calc_p_Cs_and_gamma_tillotson(particle_t* par_i)
{
  int sign;
  double rho_eos_i;
  if(is_porosity_model[par_i->property_tag]){
    rho_eos_i = par_i->alpha_por * par_i->rho_eos;
  }
  else{
    rho_eos_i = par_i->rho_eos;
  }
  
  par_i->p = calc_tillotson_pressure(par_i->u,rho_eos_i,par_i->property_tag);
  if(is_porosity_model[par_i->property_tag]){
    par_i->p /= par_i->alpha_por;
  }
  /*par_i->gamma = calc_tillotson_gamma(par_i->u,par_i->rho_eos,par_i->p,par_i->property_tag);
    par_i->Cs = sqrt( fabs( par_i->gamma * par_i->p ) / par_i->rho_eos );*/
  par_i->Cs = sqrt( calc_dp_drho_tillotson(par_i->u,rho_eos_i,par_i->p,par_i->property_tag) );
  if( isnan(par_i->Cs ) ) par_i->Cs=sqrt(A_til[par_i->property_tag]/rho0_til[par_i->property_tag]);
  if( signbit(par_i->p) ) sign = -1;
  else                    sign =  1;
  par_i->gamma = ( rho_eos_i / ( par_i->p + sign*ACC ) ) * par_i->Cs * par_i->Cs;
  if( fabs( par_i->gamma ) > 100.0 ){
    if( signbit(par_i->gamma) ) par_i->gamma = -100.0;
    else                        par_i->gamma =  100.0;
  }
}

//calculate pressure, Cs and gamma of stiffened gas equation of state
void calc_p_Cs_and_gamma_stiffened_gas_EoS(particle_t* par_i)
{
  double C0 = C0_stiffened[par_i->property_tag];
  double gamma0 = gamma0_stiffened[par_i->property_tag];
  double rho0 = rho0_stiffened[par_i->property_tag];
  int sign;
  double rho_eos_i;
  if(is_porosity_model[par_i->property_tag]){
    rho_eos_i = par_i->alpha_por * par_i->rho_eos;
  }
  else{
    rho_eos_i = par_i->rho_eos;
  }

  par_i->p = C0 * C0 * ( rho_eos_i - rho0 ) + ( gamma0 - 1.0 ) * rho_eos_i * par_i->u;
  
  if(is_porosity_model[par_i->property_tag]){
    par_i->p /= par_i->alpha_por;
  }
  par_i->Cs = sqrt( C0 * C0 + ( gamma0 - 1.0 ) * ( par_i->u + par_i->p / rho_eos_i ) );
  if( isnan(par_i->Cs) ){
    par_i->Cs = C0;
  }
  if( signbit(par_i->p) ) sign = -1;
  else                    sign =  1;
  par_i->gamma = ( rho_eos_i / par_i->p ) * par_i->Cs;
  if( fabs( par_i->gamma ) > 100.0 ){
    if( signbit(par_i->gamma) ) par_i->gamma = -100.0;
    else                        par_i->gamma =  100.0;
  }
}

//calculate delp_delrho and delp_delu for stiffened gas EoS
void calc_delp_delrho_and_delp_delu_stiffened_gas_EoS(double rho,double u,int p_tag,double* delp_delrho_P,double* delp_delu_P)
{
  double C0 = C0_stiffened[p_tag];
  double gamma0 = gamma0_stiffened[p_tag];
  double rho0 = rho0_stiffened[p_tag];

  *delp_delrho_P = C0*C0 + ( gamma0 - 1.0 )*u;
  *delp_delu_P = ( gamma0 - 1.0 )*rho;
}

//calculate pressure, Cs and gamma of Mie Gruneisen EoS
void calc_p_Cs_and_gamma_Mie_Gruneisen_EoS(particle_t* par_i)
{
  double rho0=rho0_MG[par_i->property_tag];
  double C0=C0_MG[par_i->property_tag];
  double S=S_MG[par_i->property_tag];
  double Gamma=Gamma_MG[par_i->property_tag];
  double eta;
  double rho_eos_i;
  if(is_porosity_model[par_i->property_tag]){
    rho_eos_i = par_i->alpha_por * par_i->rho_eos;
  }
  else{
    rho_eos_i = par_i->rho_eos;
  }
  
  eta = 1.0 - rho0 / rho_eos_i;
  
  par_i->p = ( rho0 * C0 * C0 * eta / pow(1.0-S*eta,2.0) ) * ( 1.0 - 0.5 * Gamma * eta ) + Gamma * rho0 * par_i->u;
  if(is_porosity_model[par_i->property_tag]){
    par_i->p /= par_i->alpha_por;
  }
  par_i->Cs = ( C0 > sqrt(fabs(Gamma*par_i->p/rho_eos_i)) ? C0 : sqrt(fabs(Gamma*par_i->p/rho_eos_i)) );
  par_i->gamma = Gamma;
}

//calculate delp_delrho and delp_delu for Mie Gruneisen EoS
void calc_delp_delrho_and_delp_delu_Mie_Gruneisen_EoS(double rho,double u,int p_tag,double* delp_delrho_P,double* delp_delu_P)
{
  double rho0=rho0_MG[p_tag];
  double C0=C0_MG[p_tag];
  double S=S_MG[p_tag];
  double Gamma=Gamma_MG[p_tag];
  double eta=1.0-(rho0/rho);

  *delp_delrho_P = pow(rho0*C0/((1.0-S*eta)*rho),2.0)*( (1.0+(2.0*eta*S/(1.0-S*eta)))*(1.0-Gamma*eta/2.0) - Gamma*eta/2.0 );
  *delp_delu_P = Gamma*rho0;
}

//calculate pressure using elastic equation of state (P=Cs^2(rho-rho0eos)).
void calc_p_elastic_eos(particle_t* par_i)
{
  double rho_eos_i;
  if(is_porosity_model[par_i->property_tag]){
    rho_eos_i = par_i->alpha_por * par_i->rho_eos;
  }
  else{
    rho_eos_i = par_i->rho_eos;
  }

  par_i->p = Cs_elastic[par_i->property_tag]*Cs_elastic[par_i->property_tag]*(rho_eos_i-rho0_elastic[par_i->property_tag]);
  if(is_porosity_model[par_i->property_tag]){
    par_i->p /= par_i->alpha_por;
  }
  par_i->Cs = Cs_elastic[par_i->property_tag];
  par_i->gamma = 0.0;
}

//calculate delp_delrho and delp_delu for elastic EoS
void calc_delp_delrho_and_delp_delu_elastic_eos(double rho,double u,int p_tag,double* delp_delrho_P,double* delp_delu_P)
{
  *delp_delrho_P = Cs_elastic[p_tag]*Cs_elastic[p_tag];
  *delp_delu_P = 0.0;
}

//calculate nabla kernel fanction and set dWdx array.It is gaussian kernel varsion.
void nablaW_gauss(double delta2,double dif[],double h,double dWdx[])
{
  int d;
  double temp;
  double h_inv2;
  h_inv2 = 1.0 / ( h * h );
  
  temp = ( -2.0 * pow( h * sqrt(PI) , -D ) * h_inv2 ) * exp( -delta2 * h_inv2 );
  for(d=0 ; d<D ; d++){
    dWdx[d] = temp * dif[d];
  }
}

//calculate gradients of density,pressure,1/rho and velosity.
void calc_gradients(particle_t par[])
{
  int i,d,d2,n;
  int particle_number=get_particle_number();
  double dWdx[D];
  double rhoi_inv2,rhoj_inv;
  double delta2,dif[D];
  particle_t *par_j;
  
#ifdef _OPENMP
#pragma omp parallel for private(i,d,d2,n,dWdx,rhoj_inv,rhoi_inv2,delta2,dif,par_j)
#endif
  for(i=0 ; i<particle_number ; i++){
    rhoi_inv2 = 1.0 /( par[i].rho * par[i].rho );
    //initialize gradients
    for(d=0 ; d<D ; d++){
      par[i].nablaV[d] = 0.0;
      par[i].nablar[d] = 0.0;
      par[i].nablap[d] = 0.0;
      for(d2=0 ; d2<D ; d2++){
	par[i].dvdx[d][d2] = 0.0;
      }
    }
    //calculate gradients
    for(n=0 ; n<par[i].inter_N ; n++){
      par_j = par[i].inter_par[n];
      delta2 = 0.0;
      for(d=0 ; d<D ; d++){
	dif[d] = par[i].r[d] - par_j->r[d];
	delta2 += dif[d] * dif[d];
      }
      if( sqrt(delta2) < Ncell * par[i].h ){
	nablaW_gauss(delta2,dif,par[i].h,dWdx);
	rhoj_inv = 1.0 / par_j->rho;
	
	for(d=0 ; d<D ; d++){
	  par[i].nablaV[d] += -rhoi_inv2 * par_j->mass * dWdx[d];
	  //calculate the gradient of density by relative density for free surface (2015/10/23)
	  par[i].nablar[d] += par_j->mass * ( par_j->rho - par[i].rho ) * rhoj_inv * dWdx[d]; 
	  par[i].nablap[d] += par_j->mass * ( par_j->p - par[i].p ) * rhoj_inv * dWdx[d];
	  for(d2=0 ; d2<D ; d2++){
	    par[i].dvdx[d][d2] += par_j->mass * ( par_j->v[d] - par[i].v[d] ) * rhoj_inv * dWdx[d2];
	  }
	}
      }
    }
  }
}

//calculate d2Vdx2 in easy way
void calc_d2Vdx2_in_easy_way(particle_t par[])
{
  int i,d1,d2,n;
  int particle_number=get_particle_number();
  double delta2,dif[D],dWdx[D];
  particle_t par_j;

  for(i=0;i<particle_number;i++){
    for(d1=0;d1<D;d1++){
      for(d2=0;d2<D;d2++){
	par[i].d2Vdx2[d1][d2]=0.0;
      }
    }
    for(n=0;n<par[i].inter_N;n++){
      if(&par[i]!=par[i].inter_par[n]){
	par_j=*par[i].inter_par[n];
	delta2=0.0;
	for(d2=0;d2<D;d2++){
	  dif[d2]=par[i].r[d2]-par_j.r[d2];
	  delta2+=dif[d2]*dif[d2];
	}
	nablaW_gauss(delta2,dif,par[i].h,dWdx);
	
	for(d1=0;d1<D;d1++){
	  for(d2=0;d2<D;d2++){
	    par[i].d2Vdx2[d1][d2] += par_j.nablaV[d1]*(par_j.mass/par_j.rho)*dWdx[d2];
	  }
	}
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






//fanctions for calculate accerelation of pressure force////////////////////////////////////////////////////////////////////////////////

void calc_sigma_for_other_structure_stand(particle_t par_i,particle_t par_j,double *temp)
{
  int d,d2;
  double Sssi,Sssj;
  double eij[D],delta2=0.0;
  double cos_theta,sin_theta;
  double cos_phi,sin_phi;

  for(d=0;d<D;d++){
    delta2 += (par_i.r[d]-par_j.r[d])*(par_i.r[d]-par_j.r[d]);
  }
  for(d=0;d<D;d++){
    eij[d] = (par_i.r[d]-par_j.r[d])/sqrt(delta2);
  }
#if D == 1
  Sssi = par_i.Sab[X][X];
  Sssj = par_j.Sab[X][X];
#endif
#if D == 2
  cos_theta = eij[X];
  sin_theta = eij[Y];
  Sssi = cos_theta*cos_theta*par_i.Sab[X][X] + sin_theta*sin_theta*par_i.Sab[Y][Y] + 2*sin_theta*cos_theta*par_i.Sab[X][Y];
  Sssj = cos_theta*cos_theta*par_j.Sab[X][X] + sin_theta*sin_theta*par_j.Sab[Y][Y] + 2*sin_theta*cos_theta*par_j.Sab[X][Y];
#endif
#if D == 3
  cos_theta = eij[Z];
  sin_theta = sqrt( 1 - eij[Z] * eij[Z] );
  if( fabs(eij[X]) < ACC && fabs(eij[Y]) < ACC ){
    cos_phi = 1.0;
    sin_phi = 0.0;
  }
  else{
    cos_phi = eij[X] / sqrt( eij[X]*eij[X] + eij[Y]*eij[Y] );
    sin_phi = eij[Y] / sqrt( eij[X]*eij[X] + eij[Y]*eij[Y] );
  }
  Sssi = sin_theta*sin_theta*cos_phi*cos_phi*par_i.Sab[X][X] + sin_theta*sin_theta*sin_phi*sin_phi*par_i.Sab[Y][Y] + cos_theta*cos_theta*par_i.Sab[Z][Z] + 2*sin_theta*sin_theta*cos_phi*sin_phi*par_i.Sab[X][Y] + 2*sin_theta*cos_theta*cos_phi*par_i.Sab[X][Z] + 2*sin_theta*cos_theta*sin_phi*par_i.Sab[Y][Z];
  Sssj = sin_theta*sin_theta*cos_phi*cos_phi*par_j.Sab[X][X] + sin_theta*sin_theta*sin_phi*sin_phi*par_j.Sab[Y][Y] + cos_theta*cos_theta*par_j.Sab[Z][Z] + 2*sin_theta*sin_theta*cos_phi*sin_phi*par_j.Sab[X][Y] + 2*sin_theta*cos_theta*cos_phi*par_j.Sab[X][Z] + 2*sin_theta*cos_theta*sin_phi*par_j.Sab[Y][Z];
#endif
  *temp += par_j.mass * ( Sssi / ( par_i.rho * par_i.rho ) + Sssj / ( par_j.rho * par_j.rho ) );
  
  if(*temp>0.0) *temp=0.0;

}

//calculate aij using standard SPH
void calc_aij_standard(particle_t *par_iP,particle_t par_j,double *commonP)
{
  int d,d2;
  double hij = 0.5*( par_iP->h + par_j.h );
  double hij_inv2 = 1.0 / ( hij * hij );
  double Csij = 0.5*( par_iP->Cs + par_j.Cs );
  double muij,PIij;
  double temp;
  double dif[D];
  double dif_v[D];
  double delta2=0.0;
  double dif_dot_difv=0.0;
  double dWdx[D];
  
  //calculate artificial viscosity
  for(d=0 ; d<D ; d++){
    dif[d] = par_iP->r[d] - par_j.r[d];
    dif_v[d] = par_iP->v[d] - par_j.v[d];
    delta2 += dif[d] * dif[d];
    dif_dot_difv += dif[d] * dif_v[d];
  }
  //if two particle go away each other,artificial viscosity is nonactive.
  if( dif_dot_difv > ACC ){
    PIij = 0.0;
  }
  else{
    muij = ( hij * dif_dot_difv ) / ( delta2 + 0.01 * hij * hij );
    PIij = ( -alpha_vis * Csij * muij + beta_vis * muij * muij ) * 2.0 / ( par_iP->rho + par_j.rho ) ;
  }
  
  //calculate aij
  nablaW_gauss( delta2 , dif , hij , dWdx );
  temp = - par_j.mass * ( par_iP->p / ( par_iP->rho * par_iP->rho ) + par_j.p / ( par_j.rho * par_j.rho ) + PIij );
  for(d=0 ; d<D ; d++){
    par_iP->a[d] += temp * dWdx[d];
    //add deviatric part
    for(d2=0 ; d2<D ; d2++){
      par_iP->a[d] += par_j.mass * ( par_iP->Sab[d][d2] / ( par_iP->rho * par_iP->rho ) + par_j.Sab[d][d2] / ( par_j.rho * par_j.rho ) ) * dWdx[d2];
    }
  }
  
  *commonP = temp;
}

//calculate dudtij using standard SPH
double calc_dudtij_standard(particle_t par_i,particle_t par_j,double next_v[],double common)
{	
  int d,d2;
  double delta2=0.0;
  double dif[D];
  double difv_dot_dWdx=0.0;
  double dWdx[D];
  double dev_part=0.0;
  
  //calculate dudtij
  common *= -0.5; //times raito of common factor of a and dudt
  for(d=0 ; d<D ; d++){
    dif[d] = par_i.r[d] - par_j.r[d];
    delta2 += dif[d] * dif[d];
  }
  nablaW_gauss( delta2 , dif , 0.5 * ( par_i.h + par_j.h ) , dWdx );
  for(d=0 ; d<D ; d++){
    difv_dot_dWdx += ( /*next_v[d]*/par_i.v[d] - par_j.v[d] ) * dWdx[d];
  }
  
  //calculate deviatric part
  for(d=0 ; d<D ; d++){
    for(d2=0 ; d2<D ; d2++){
      dev_part += ( ( par_i.Sab[d][d2] * par_j.mass ) / ( par_i.rho * par_i.rho * 2.0 ) ) * ( ( par_j.v[d] - par_i.v[d] ) * dWdx[d2] + ( par_j.v[d2] - par_i.v[d2] ) * dWdx[d] );
    }
  }
    
  return( common * difv_dot_dWdx + dev_part);
}

//calculate dSabdt standard
void calc_dSab_rho_dt_standard(particle_t *par_iP, particle_t par_j)
{
  int d,d2,d3;
  double depsij_rho[D][D];
  double Rij_rho[D][D];
  double depsggij_rho=0.0; //diagonal part of dot epsilon / density
  double SiRijgg_rho[D][D]; //S^{alpha gamma}_{i}*R^{beta gamma}_{ij} / density
  for(d=0 ; d<D ; d++){
    for(d2=0 ; d2<D ; d2++){
      SiRijgg_rho[d][d2] = 0.0;
    }
  }
  double dif[D];
  double delta2=0.0;
  double dWdx[D];
  for(d=0 ; d<D ; d++){
    for(d2=0 ; d2<D ; d2++){
      SiRijgg_rho[d][d2] = 0.0;
    }
  }

  for(d=0; d<D ; d++){
    dif[d] = par_iP->r[d] - par_j.r[d];
    delta2 += dif[d] * dif[d];
  }
  nablaW_gauss( delta2, dif, 0.5 * ( par_iP->h + par_j.h ), dWdx );

  for(d=0 ; d<D ; d++){
    for(d2=0 ; d2<=d ; d2++){
      depsij_rho[d][d2] = depsij_rho[d2][d] = ( par_j.mass / ( 2.0 * par_iP->rho * par_iP->rho ) ) * ( ( par_j.v[d] - par_iP->v[d] ) * dWdx[d2] + ( par_j.v[d2] - par_iP->v[d2] ) * dWdx[d] );
      Rij_rho[d][d2] = ( par_j.mass / ( 2.0 * par_iP->rho * par_iP->rho ) ) * ( ( par_j.v[d] - par_iP->v[d] ) * dWdx[d2] - ( par_j.v[d2] - par_iP->v[d2] ) * dWdx[d] );
      Rij_rho[d2][d] = -Rij_rho[d][d2];
    }
  }
  for(d=0 ; d<D ; d++){
    depsggij_rho += depsij_rho[d][d];
  }
  for(d=0 ; d<D ; d++){
    for(d2=0 ; d2<D ; d2++){
      for(d3=0 ; d3<D ; d3++){
	SiRijgg_rho[d][d2] += par_iP->Sab[d][d3] * Rij_rho[d2][d3];
      }
    }
  }

  //calculate time change rate of deviatric stress tensor / density
  for(d=0 ; d<D ; d++){
    for(d2=0 ; d2<D ; d2++){
      par_iP->dSab_rho_dt[d][d2] += 2.0 * mu_shear[par_iP->property_tag] * (depsij_rho[d][d2] - (1.0/3.0) * ((d==d2)?1.0:0.0) * depsggij_rho) + SiRijgg_rho[d][d2] + SiRijgg_rho[d2][d] + par_iP->Sab[d][d2] * depsggij_rho;
    }
  }

}

//calclulate Vij2 and s star using linear interpolation
void calc_Vij2_using_linear(particle_t par_i,particle_t par_j,inter_t *inter_N,double *ssP,double delta,double eij[])
{
  int d;
  double Cij,Dij;
  double Vi = 1.0 / par_i.rho;
  double Vj = 1.0 / par_j.rho;
  double hi2 = par_i.h * par_i.h;
  double hj2 = par_j.h * par_j.h;

  Cij = ( Vi - Vj ) / delta;
  Dij = 0.5 * ( Vi + Vj );
  inter_N->Vij2_hi = 0.25 * hi2 * Cij * Cij + Dij * Dij;
  inter_N->Vij2_hj = 0.25 * hj2 * Cij * Cij + Dij * Dij;
  *ssP = 0.5 * ( hi2 * Cij * Dij / ( 2 * inter_N->Vij2_hi ) + hj2 * Cij * Dij / ( 2 * inter_N->Vij2_hj ) );
}

//calclulate Vij2 and s star using cubic spline interpolation
void calc_Vij2_using_cubic_spline(particle_t par_i,particle_t par_j,inter_t *inter_N,double *ssP,double delta,double eij[])
{
  int d;
  double dVi=0.0;
  double dVj=0.0;
  double Aij,Bij,Cij,Dij;
  double Vi = 1.0 / par_i.rho;
  double Vj = 1.0 / par_j.rho;
  double hi2 = par_i.h * par_i.h;
  double hj2 = par_j.h * par_j.h;
  
  for(d=0 ; d<D ; d++){
    dVi += par_i.nablaV[d] * eij[d];
    dVj += par_j.nablaV[d] * eij[d];
  }
  
  Aij = -2 * ( Vi - Vj ) * pow(delta,-3) + ( dVi + dVj ) * pow(delta,-2);
  Bij = 0.5 * ( dVi - dVj ) / delta;
  Cij = 1.5 *( Vi - Vj ) / delta - 0.25 * ( dVi + dVj );
  Dij = 0.5 * ( Vi + Vj ) - 0.125 * ( dVi - dVj ) * delta;
  inter_N->Vij2_hi = ( ( 0.234375 * hi2 * Aij * Aij + 0.1875 * ( 2 * Aij * Cij + Bij * Bij ) ) * hi2 + 0.25 * ( 2 * Bij * Dij + Cij * Cij ) ) * hi2 + Dij * Dij;
  inter_N->Vij2_hj = ( ( 0.234375 * hj2 * Aij * Aij + 0.1875 * ( 2 * Aij * Cij + Bij * Bij ) ) * hj2 + 0.25 * ( 2 * Bij * Dij + Cij * Cij ) ) * hj2 + Dij * Dij;
  *ssP = 0.5 * ( ( ( 0.46875 * hi2 * Aij * Bij + 0.375 * ( Aij * Dij + Bij * Cij ) ) * hi2 + 0.5 * Cij * Dij ) * hi2 / inter_N->Vij2_hi + ( ( 0.46875 * hj2 * Aij * Bij + 0.375 * ( Aij * Dij + Bij * Cij ) ) * hj2 + 0.5 * Cij * Dij ) * hj2 / inter_N->Vij2_hj );
}

//calculate Vij using quintic spline interpolation.But now it is 1 dimension version.
void calc_Vij2_using_quintic_spline(particle_t par_i,particle_t par_j,inter_t *inter_N,double *ssP,double delta,double eij[])
{
  int d,d1,d2;
  double dVi=0.0;
  double dVj=0.0;
  double dV2i=0.0;
  double dV2j=0.0;
  double Aij,Bij,Cij,Dij,Eij,Fij;
  double Vi = 1.0 / par_i.rho;
  double Vj = 1.0 / par_j.rho;
  double hi2 = par_i.h * par_i.h;
  double hj2 = par_j.h * par_j.h;
  double ssi,ssj;

  for(d=0 ; d<D ; d++){
    dVi += par_i.nablaV[d] * eij[d];
    dVj += par_j.nablaV[d] * eij[d];
  }
  for(d1=0 ; d1<D ; d1++){
    for(d2=0 ; d2<D ; d2++){
      dV2i += par_i.d2Vdx2[d1][d2]*eij[d1]*eij[d2];
      dV2j += par_j.d2Vdx2[d1][d2]*eij[d1]*eij[d2];
    }
  }
  
  //quintic spline interpolation
  Aij = 6.0 * ( Vi - Vj ) * pow(delta,-5.0) - 3.0 * ( dVi + dVj ) * pow(delta,-4.0) + 0.5 * ( dV2i - dV2j ) * pow(delta,-3.0);
  Bij = -0.5 * ( dVi - dVj ) * pow(delta,-3.0) + 0.25 * ( dV2i + dV2j ) * pow(delta,-2.0);
  Cij = -5.0 * ( Vi - Vj ) * pow(delta,-3.0) + 2.5 * ( dVi + dVj ) * pow(delta,-2.0) - 0.25 * ( dV2i - dV2j ) / delta;
  Dij = 0.75 * ( dVi - dVj ) / delta - 0.125 * ( dV2i + dV2j );
  Eij = 1.875 * ( Vi - Vj ) / delta - 0.4375 * ( dVi + dVj ) + 0.03125 * ( dV2i - dV2j ) * delta;
  Fij = 0.5 * ( Vi + Vj ) - 0.15625 * ( dVi - dVj ) * delta + 0.015625 * ( dV2i + dV2j ) * delta * delta;

  inter_N->Vij2_hi = ( ( ( ( 0.92285156 * Aij * Aij * hi2 + 0.41015625 * ( 2.0 * Aij * Cij + Bij * Bij ) ) * hi2 + 0.234375 * ( 2.0 * Aij * Eij + 2.0 * Bij * Dij + Cij * Cij ) ) * hi2 + 0.1875 * ( 2.0 * Bij * Fij + 2.0 * Cij * Eij + Dij * Dij ) ) * hi2 + 0.25 * ( 2.0 * Dij * Fij + Eij * Eij ) ) * hi2 + Fij * Fij;
  inter_N->Vij2_hj = ( ( ( ( 0.92285156 * Aij * Aij * hj2 + 0.41015625 * ( 2.0 * Aij * Cij + Bij * Bij ) ) * hj2 + 0.234375 * ( 2.0 * Aij * Eij + 2.0 * Bij * Dij + Cij * Cij ) ) * hj2 + 0.1875 * ( 2.0 * Bij * Fij + 2.0 * Cij * Eij + Dij * Dij ) ) * hj2 + 0.25 * ( 2.0 * Dij * Fij + Eij * Eij ) ) * hj2 + Fij * Fij;

  ssi = ( ( ( ( ( 1.84570313 * Aij * Bij * hi2 + 0.8203125 * ( Aij * Dij + Bij * Cij ) ) * hi2 + 0.46875 * ( Aij * Fij + Bij * Eij + Cij * Dij ) ) * hi2 + 0.375 * ( Cij * Fij + Dij * Eij ) ) * hi2 + 0.5 * Eij * Fij ) * hi2 ) / inter_N->Vij2_hi;
  ssj = ( ( ( ( ( 1.84570313 * Aij * Bij * hj2 + 0.8203125 * ( Aij * Dij + Bij * Cij ) ) * hj2 + 0.46875 * ( Aij * Fij + Bij * Eij + Cij * Dij ) ) * hj2 + 0.375 * ( Cij * Fij + Dij * Eij ) ) * hj2 + 0.5 * Eij * Fij ) * hj2 ) / inter_N->Vij2_hj;
  *ssP = 0.5 * ( ssi + ssj );
  
}


//solve riemann problem for ideal gas eos
void solve_riemann_problem_for_ideal_gas_eos(double p1,double rho1,double v1,double gamma1,double p2,double rho2,double v2,double gamma2,double* p_riemannP,double* v_riemannP)
{
  double ppre,p;
  double W1,W2;
  int i;
  double gamma = ( gamma1 < gamma2 ? gamma1 : gamma2 );
  double alpha= ( 2.0 * gamma ) / ( gamma - 1.0 );
  int roop=4;
  p = ( p1 + p2 ) / 2.0;
  for(i=0 ; i<roop ; i++){
    ppre = p;
    if( p < p1 - ACC ) W1 = sqrt( p1 * rho1 * gamma ) * ( ( gamma - 1.0 ) / ( 2.0 * gamma ) ) * ( 1.0 - p / p1 ) / ( 1.0 - pow(p/p1 , 1.0/alpha) );
    else W1 = sqrt(p1 * rho1 * gamma ) * sqrt( 0.5 * ( gamma + 1.0 ) * ( ( p - p1 ) / p1 ) + 1.0 );
    if( p < p2 - ACC ) W2 = sqrt( p2 * rho2 * gamma ) * ( ( gamma - 1.0 ) / ( 2.0 * gamma ) ) * ( 1.0 - p / p2 ) / ( 1.0 - pow(p/p2 , 1.0/alpha) );
    else W2 = sqrt(p2 * rho2 * gamma ) * sqrt( 0.5 * ( gamma + 1.0 ) * ( ( p - p2 ) / p2 ) + 1.0 );
    p = ( ( p2 / W2 + p1 / W1 ) + v2 - v1 ) / ( 1.0 / W2 + 1.0 / W1 );
    if( fabs( p - ppre ) < ACC ) break;
  }
  *p_riemannP = p;
  *v_riemannP = ( ( W1 * v1 + W2 * v2 ) + p2 - p1 ) / ( W1 + W2 );
  
  if( isnan(p)  || isnan(*v_riemannP) ){
    *p_riemannP = 0.5 * ( p1 + p2 );
    *v_riemannP = 0.5 * ( v1 + v2 );
  }
}

//solve riemann problem for elastic EoS of P=Cs^2(rho-rho0eos)
void solve_riemann_problem_for_elastic_eos(double p1,double rho1,double v1,double Csi,double p2,double rho2,double v2,double Csj,double* p_riemannP,double* v_riemannP)
{
  double ppre,p;
  double W1,W2;
  int i;	
  double Csij = 0.5*(Csi+Csj);
  double rho0eosij = 0.5*( ( rho1 - p1/(Csij*Csij) ) + ( rho2 - p2/(Csij*Csij) ) );
  
  p = ( p1 + p2 ) / 2.0;
  for(i=0 ; i<4 ; i++){
    ppre = p;
    if( p < p1 - ACC ) W1 = fabs( ( (p-p1) / Csij ) * ( 1.0 / log( Csij*Csij*rho1 / (p+Csij*Csij*rho0eosij) ) ) );
    else W1 = Csij * sqrt( rho1 * ( rho1 + (p-p1)/(Csij*Csij) ) );
    if( p < p2 - ACC ) W2 = fabs( ( (p-p2) / Csij ) * ( 1.0 / log( Csij*Csij*rho2 / (p+Csij*Csij*rho0eosij) ) ) );
    else W2 = Csij * sqrt( rho2 * ( rho2 + (p-p2)/(Csij*Csij) ) );
    p = ( ( p2 / W2 + p1 / W1 ) + v2 - v1 ) / ( 1.0 / W2 + 1.0 / W1 );
  }
  *p_riemannP = p;
  *v_riemannP = ( ( W1 * v1 + W2 * v2 ) + p2 - p1 ) / ( W1 + W2 );
  
  //if riemann solver become not a number in elastic eos, this situation maybe occurs in rarefaction-like place. So we do not need to solve riemann solver
  if( isnan(p)  || isnan(*v_riemannP) ){
    *p_riemannP = 0.5 * ( p1 + p2 );
    *v_riemannP = 0.5 * ( v1 + v2 );
  }
}

//strengthen monotonicity constraint for riemann solver
void set_monotonicity_constraint(particle_t par[])
{
  int i,n,d;
  particle_t par_j;
  int particle_number=get_particle_number();
  double tent_hi;
  double drdsi,drdsj,dpdsi,dpdsj;
  double eij[D];
  double delta2,dif[D],delta;

  for(i=0;i<particle_number;i++){
    par[i].monoto_flag=0;
  }

  for(i=0;i<particle_number;i++){
    tent_hi=pow(par[i].mass/par[i].rho,1.0/D);
    for(n=0;n<par[i].inter_N;n++){
      if(&par[i]!=par[i].inter_par[n]){
	par_j = *(par[i].inter_par[n]);
	
	delta2=0.0;
	drdsi=drdsj=dpdsi=dpdsj=0.0;
	for(d=0;d<D;d++){
	  dif[d]=par[i].r[d]-par_j.r[d];
	  delta2+=dif[d]*dif[d];
	}
	delta=sqrt(delta2);
	for(d=0;d<D;d++){
	  eij[d]=dif[d]/delta;
	}
	if( delta < search_region_monoto*tent_hi ){
	  for(d=0 ; d<D ; d++){
	    drdsi += par[i].nablar[d] * eij[d];
	    drdsj += par_j.nablar[d] * eij[d];
	    dpdsi += par[i].nablap[d] * eij[d];
	    dpdsj += par_j.nablap[d] * eij[d];
	  }
	  if( (drdsi*drdsj < ACC) || (dpdsi*dpdsj < ACC) ){
	    par[i].monoto_flag=1;
	  }
	}
      }
    }
  }
}

//decide riemann solver criterion. If we use riemann solver for elastic equation of state, return 1. If we use that for ideal gas equation of state, return 0.
int riemann_solver_criterion(particle_t par_i,particle_t par_j)
{
  if(equation_of_state[par_i.property_tag]==1||equation_of_state[par_j.property_tag]==1){
    //for elastic equation of state
    return(1);
  }
  else if(equation_of_state[par_i.property_tag]==0||equation_of_state[par_j.property_tag]==0){
    //for ideal gas equation of state
    return(0);
  }
  else{
    //for other equation of state
    if( par_i.p<0.0 || par_j.p<0.0 || 0.5*(u_crit_Riemann[par_i.property_tag]+u_crit_Riemann[par_j.property_tag]) > 0.5*(par_i.u+par_j.u) ) return(1);
    else return(0);
  }
}

//calculate riemann solver
void calc_riemann_solver(particle_t par_i,particle_t par_j,inter_t *inter_N,double ss,double eij[],double delta,double dt)
{
  int d,d2;
  double drdsi=0.0,drdsj=0.0,dpdsi=0.0,dpdsj=0.0,dvdsi=0.0,dvdsj=0.0,vi=0.0,vj=0.0,gi=0.0,gj=0.0;
  double rhoR,rhoL,pR,pL,vR,vL,gammaR,gammaL;
  
  for(d=0 ; d<D ; d++){
    vi += par_i.v[d] * eij[d];
    vj += par_j.v[d] * eij[d];
    if(is_selfgravity){
      gi += par_i.g[d] * eij[d];
      gj += par_j.g[d] * eij[d];
    }
  }
  //if space second order or not in the shock wave,calc s direction's gradients 
  if( (space_order==2) && (3*(vj-vi) < (par_i.Cs<par_j.Cs)?par_i.Cs:par_j.Cs) && par_i.monoto_flag==0 && par_j.monoto_flag==0 ){
    for(d=0 ; d<D ; d++){
      drdsi += par_i.nablar[d] * eij[d];
      drdsj += par_j.nablar[d] * eij[d];
      dpdsi += par_i.nablap[d] * eij[d];
      dpdsj += par_j.nablap[d] * eij[d];
      for(d2=0 ; d2<D ; d2++){
	dvdsi += par_i.dvdx[d][d2] * eij[d] * eij[d2];
	dvdsj += par_j.dvdx[d][d2] * eij[d] * eij[d2];
      }
    }
    //monotonisity constraint
    if( dvdsi * dvdsj < 0.0 ) dvdsi = dvdsj = 0.0;
    if( drdsi * drdsj < 0.0 ) drdsi = drdsj = 0.0;
    if( dpdsi * dpdsj < 0.0 ) dpdsi = dpdsj = 0.0;
  }
  
  //set initial condition of riemann problem
  rhoR = par_i.rho + drdsi * ( ss + par_i.Cs * 0.5 * dt -0.5 * delta );
  pR = par_i.p + dpdsi * ( ss + par_i.Cs * 0.5 * dt -0.5 * delta );
  vR = vi + dvdsi * ( ss + par_i.Cs * 0.5 * dt -0.5 * delta );
  gammaR = par_i.gamma;
  rhoL = par_j.rho + drdsj * ( ss - par_j.Cs * 0.5 * dt + 0.5 * delta );
  pL = par_j.p + dpdsj * ( ss - par_j.Cs * 0.5 * dt + 0.5 * delta );
  vL = vj + dvdsj * ( ss - par_j.Cs * 0.5 * dt + 0.5 * delta );
  gammaL = par_j.gamma;
  if(rhoR < 0.0) rhoR = par_i.rho;
  if(rhoL < 0.0) rhoL = par_j.rho;
  if(is_selfgravity){
    pR += - 0.5 * par_i.Cs * dt * rhoR * gi;
    pL += 0.5 * par_j.Cs * dt * rhoL * gj;
  }
  
  if( riemann_solver_criterion(par_i,par_j) ){
    solve_riemann_problem_for_elastic_eos(pR,rhoR,vR,par_i.Cs,pL,rhoL,vL,par_j.Cs,&inter_N->p_riemann,&inter_N->v_riemann);
  }
  else{
    solve_riemann_problem_for_ideal_gas_eos(pR,rhoR,vR,par_i.gamma,pL,rhoL,vL,par_j.gamma,&inter_N->p_riemann,&inter_N->v_riemann);
  }
  //maybe riemann solver of velosity is not good in simulation of collision for non-ideal gas equation of state
  if(equation_of_state[par_i.property_tag]!=0||equation_of_state[par_j.property_tag]!=0){
    inter_N->v_riemann=0.5*(vi+vj);
  }
  //in the case of second-oreder riemann solver, this can unphysically makes resultant pressure of riemann solver negative at shock surface, and negative internal energy. (2015/11/25)
  if( par_i.p>0.0 && par_j.p>0.0 && inter_N->p_riemann<0.0) inter_N->p_riemann = 0.5*(par_i.p+par_j.p);
  
}

void calc_sigma_for_other_structure(double eij[],inter_t *inter_N)
{
  int d,d2;
  double Sssijast;
  double cos_theta,sin_theta;
  double cos_phi,sin_phi;
#if D == 1
  Sssijast = inter_N->Sijs[X][X];
#endif
#if D == 2
  cos_theta = eij[X];
  sin_theta = eij[Y];
  Sssijast = cos_theta*cos_theta*inter_N->Sijs[X][X] + sin_theta*sin_theta*inter_N->Sijs[Y][Y] + 2*sin_theta*cos_theta*inter_N->Sijs[X][Y];
#endif
#if D == 3
  cos_theta = eij[Z];
  sin_theta = sqrt( 1 - eij[Z] * eij[Z] );
  if( fabs(eij[X]) < ACC && fabs(eij[Y]) < ACC ){
    cos_phi = 1.0;
    sin_phi = 0.0;
  }
  else{
    cos_phi = eij[X] / sqrt( eij[X]*eij[X] + eij[Y]*eij[Y] );
    sin_phi = eij[Y] / sqrt( eij[X]*eij[X] + eij[Y]*eij[Y] );
  }
  Sssijast = sin_theta*sin_theta*cos_phi*cos_phi*inter_N->Sijs[X][X] + sin_theta*sin_theta*sin_phi*sin_phi*inter_N->Sijs[Y][Y] + cos_theta*cos_theta*inter_N->Sijs[Z][Z] + 2*sin_theta*sin_theta*cos_phi*sin_phi*inter_N->Sijs[X][Y] + 2*sin_theta*cos_theta*cos_phi*inter_N->Sijs[X][Z] + 2*sin_theta*cos_theta*sin_phi*inter_N->Sijs[Y][Z];
#endif

  for(d=0;d<D;d++){
    for(d2=0;d2<D;d2++){
      inter_N->Sijs[d][d2] = 0.0;
    }
  }
  
  if(inter_N->p_riemann<0.0){
    inter_N->p_riemann = 0.0;
  }
  if(Sssijast>0.0){
    Sssijast = 0.0;
  }
  inter_N->p_riemann -= Sssijast;

}

//calculate aij using godnov SPH
void calc_aij_riemann(particle_t *par_iP,particle_t par_j,inter_t *inter_N,double dt)
{
  int d,d2;
  double eij[D];
  double delta2=0.0;
  double ss;
  double criterion=par_iP->p+par_j.p;
  double vij=0.0;
  double vij_vector[D];
  double rhoij=0.5*(par_iP->rho+par_j.rho);;
  double Vij2_crit=10.0*(1.0/(rhoij*rhoij));

  double aij;
  
  for(d=0 ; d<D ; d++){
    inter_N->dif[d] = par_iP->r[d] - par_j.r[d];
    delta2 += inter_N->dif[d] * inter_N->dif[d];
  }
  inter_N->delta = sqrt(delta2);
  
  for(d=0 ; d<D ; d++){
    eij[d] = inter_N->dif[d] / inter_N->delta;
  }

  /*//criterion that include deviatric stress tensor/////////////////////////////////////////////////////////
  double Sssi,Sssj;
  double cos_theta,sin_theta;
  double cos_phi,sin_phi;
#if D == 2
  cos_theta = eij[X];
  sin_theta = eij[Y];
  Sssi = cos_theta*cos_theta*par_iP->Sab[X][X] + sin_theta*sin_theta*par_iP->Sab[Y][Y] + 2*sin_theta*cos_theta*par_iP->Sab[X][Y];
  Sssj = cos_theta*cos_theta*par_j.Sab[X][X] + sin_theta*sin_theta*par_j.Sab[Y][Y] + 2*sin_theta*cos_theta*par_j.Sab[X][Y];
  criterion = ( par_iP->p - Sssi ) + ( par_j.p - Sssj);
#endif
#if D == 3
  cos_theta = eij[Z];
  sin_theta = sqrt( 1 - eij[Z] * eij[Z] );
  cos_phi = eij[X] / sqrt( eij[X]*eij[X] + eij[Y]*eij[Y] );
  sin_phi = eij[Y] / sqrt( eij[X]*eij[X] + eij[Y]*eij[Y] );
  Sssi = sin_theta*sin_theta*cos_phi*cos_phi*par_iP->Sab[X][X] + sin_theta*sin_theta*sin_phi*sin_phi*par_iP->Sab[Y][Y] + cos_theta*cos_theta*par_iP->Sab[Z][Z] + 2*sin_theta*sin_theta*cos_phi*sin_phi*par_iP->Sab[X][Y] + 2*sin_theta*cos_theta*cos_phi*par_iP->Sab[X][Z] + 2*sin_theta*cos_theta*sin_phi*par_iP->Sab[Y][Z];
  Sssj = sin_theta*sin_theta*cos_phi*cos_phi*par_j.Sab[X][X] + sin_theta*sin_theta*sin_phi*sin_phi*par_j.Sab[Y][Y] + cos_theta*cos_theta*par_j.Sab[Z][Z] + 2*sin_theta*sin_theta*cos_phi*sin_phi*par_j.Sab[X][Y] + 2*sin_theta*cos_theta*cos_phi*par_j.Sab[X][Z] + 2*sin_theta*cos_theta*sin_phi*par_j.Sab[Y][Z];
  criterion = ( par_iP->p - Sssi ) + ( par_j.p - Sssj);
#endif
  /////////////////////////////////////////////////////////////////////////////////////////////////////////*/
  
  //calculate Vij2 and s star
  //selecting interpolation type to suppress the tensile instability according to Sugiura and Inutsuka 2015. 
#if D == 1
  if(criterion>0.0){
    calc_Vij2_using_cubic_spline(*par_iP,par_j,inter_N,&ss,inter_N->delta,eij);
  }
  else{
    calc_Vij2_using_quintic_spline(*par_iP,par_j,inter_N,&ss,inter_N->delta,eij);
  }
#else
  if(criterion>0.0){
    calc_Vij2_using_linear(*par_iP,par_j,inter_N,&ss,inter_N->delta,eij);
  }
  else{
    calc_Vij2_using_cubic_spline(*par_iP,par_j,inter_N,&ss,inter_N->delta,eij);
  }
#endif
  //if Vij2 is too large because of interpolation, we limit Vij2.
  if( inter_N->Vij2_hi > Vij2_crit || inter_N->Vij2_hj > Vij2_crit ) calc_Vij2_using_linear(*par_iP,par_j,inter_N,&ss,inter_N->delta,eij);
    
  inter_N->ss = ss;
  //calclulate riemann solver
  calc_riemann_solver(*par_iP,par_j,inter_N,ss,eij,inter_N->delta,dt);
  //for low density particle, we do not use riemann solver (2015/12/21)
  if(par_iP->low_density_flag==1||par_j.low_density_flag==1){
    inter_N->p_riemann = 0.5 * ( par_iP->p + par_j.p );
  }
  
  nablaW_gauss( delta2 , inter_N->dif , sqrt(2.0)*par_iP->h , inter_N->dWdx_hi );
  nablaW_gauss( delta2 , inter_N->dif , sqrt(2.0)*par_j.h , inter_N->dWdx_hj );

  //test artificial viscosity//////////////////////////////////////////////////////////////////////////////////////
  /*double hij = 0.5*( par_iP->h + par_j.h );
  double Csij = 0.5*( par_iP->Cs + par_j.Cs );
  double muij,PIij;
  double dif_v[D];
  double dif_dot_difv=0.0;
  for(d=0 ; d<D ; d++){
    dif_v[d] = par_iP->v[d] - par_j.v[d];
    dif_dot_difv += inter_N->dif[d] * dif_v[d];
  }
  muij = ( hij * dif_dot_difv ) / ( delta2 + 0.1 * hij * hij );
  PIij = ( -alpha_vis * Csij * muij + beta_vis * muij * muij ) * 2.0 / ( par_iP->rho + par_j.rho );
  inter_N->p_riemann=0.5*(par_iP->p+par_j.p);*/
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //calculate Sij[D][D] star and vij[D] star
  for(d=0 ; d<D ; d++){
    vij_vector[d] = 0.5 * ( par_iP->v[d] + par_j.v[d] ) + ( par_iP->v[d] - par_j.v[d] ) * ( ss / inter_N->delta );
    vij += vij_vector[d] * eij[d];
  }
  for(d=0 ; d<D ; d++){
    for(d2=0 ; d2<=d ; d2++){
      inter_N->Sijs[d][d2] = inter_N->Sijs[d2][d] = 0.5 * ( par_iP->Sab[d][d2] + par_j.Sab[d][d2] ) + ( par_iP->Sab[d][d2] - par_j.Sab[d][d2] ) * ( ss / inter_N->delta );
    }
    inter_N->vijs[d] = inter_N->v_riemann * eij[d] + ( vij_vector[d] - vij * eij[d] );
  }
  
  //calculate aij with godnov method
  for(d=0 ; d<D ; d++){
    par_iP->a[d] += - par_j.mass * inter_N->p_riemann * ( inter_N->Vij2_hi * inter_N->dWdx_hi[d] + inter_N->Vij2_hj * inter_N->dWdx_hj[d] );
    //par_iP->a[d] += - par_j.mass *( inter_N->p_riemann * ( inter_N->Vij2_hi *inter_N-> dWdx_hi[d] + inter_N->Vij2_hj * inter_N->dWdx_hj[d] ) + PIij * 0.5 * (inter_N->dWdx_hi[d] + inter_N->dWdx_hj[d]) );
    //add deviatoric part
    for(d2=0 ; d2<D ; d2++){
      par_iP->a[d] += par_j.mass * inter_N->Sijs[d][d2] * ( inter_N->Vij2_hi * inter_N->dWdx_hi[d2] + inter_N->Vij2_hj * inter_N->dWdx_hj[d2] );
    }
  }
}

//calculate dudtij using godnov SPH
double calc_dudtij_riemann(particle_t par_i,particle_t par_j,double tcenter_v[],inter_t inter_n)
{
  int d,d2;
  double ret=0.0;
  
  for(d=0 ; d<D ; d++){
    for(d2=0 ; d2<=d ; d2++){
      ret += ((d==d2)?1.0:2.0) * 0.5 * par_j.mass * ( -inter_n.p_riemann * ((d==d2)?1.0:0.0) + inter_n.Sijs[d][d2] ) * ( inter_n.Vij2_hi * ( ( inter_n.vijs[d] - tcenter_v[d] ) * inter_n.dWdx_hi[d2] + ( inter_n.vijs[d2] - tcenter_v[d2] ) * inter_n.dWdx_hi[d] ) + inter_n.Vij2_hj * ( ( inter_n.vijs[d] - tcenter_v[d] ) * inter_n.dWdx_hj[d2] + ( inter_n.vijs[d2] - tcenter_v[d2] ) * inter_n.dWdx_hj[d] ) );
    }
  } 
  
  return(ret);
}

//calculate dSab_rho_dt using godunov SPH
void calc_dSab_rho_dt_riemann(particle_t *par_iP, particle_t par_j, inter_t inter_n)
{
  int d,d2,d3;
  double depsij_rho[D][D];
  double Rij_rho[D][D];
  double depsggij_rho=0.0; //diagonal part of dot epsilon/rho
  double SiRijgg_rho[D][D]; //S^{alpha gamma}_{i}*R^{beta gamma}_{ij}/rho
  for(d=0 ; d<D ; d++){
    for(d2=0 ; d2<D ; d2++){
      SiRijgg_rho[d][d2] = 0.0;
    }
  }

  for(d=0 ; d<D ; d++){
    for(d2=0 ; d2<=d ; d2++){
      depsij_rho[d][d2] = depsij_rho[d2][d] = 0.5 * par_j.mass * ( inter_n.Vij2_hi * ( ( inter_n.vijs[d] - par_iP->v[d] ) * inter_n.dWdx_hi[d2] + ( inter_n.vijs[d2] - par_iP->v[d2] ) * inter_n.dWdx_hi[d] ) + inter_n.Vij2_hj * ( ( inter_n.vijs[d] - par_iP->v[d] ) * inter_n.dWdx_hj[d2] + ( inter_n.vijs[d2] - par_iP->v[d2] ) * inter_n.dWdx_hj[d] ) );
      Rij_rho[d][d2] = 0.5 * par_j.mass * ( inter_n.Vij2_hi * ( ( inter_n.vijs[d] - par_iP->v[d] ) * inter_n.dWdx_hi[d2] - ( inter_n.vijs[d2] - par_iP->v[d2] ) * inter_n.dWdx_hi[d] ) + inter_n.Vij2_hj * ( ( inter_n.vijs[d] - par_iP->v[d] ) * inter_n.dWdx_hj[d2] - ( inter_n.vijs[d2] - par_iP->v[d2] ) * inter_n.dWdx_hj[d] ) );
      Rij_rho[d2][d] = -Rij_rho[d][d2];
    }
  }
  for(d=0 ; d<D ; d++){
    depsggij_rho += depsij_rho[d][d];
    for(d2=0 ; d2<D ; d2++){
      for(d3=0 ; d3<D ; d3++){
	SiRijgg_rho[d][d2] += inter_n.Sijs[d][d3] * Rij_rho[d2][d3];
      }
    }
  }

  //calculate time change of deviatoric stress tensor / density
  for(d=0 ; d<D ; d++){
    for(d2=0 ; d2<D ; d2++){
      par_iP->dSab_rho_dt[d][d2] += 2.0 * mu_shear[par_iP->property_tag] * ( depsij_rho[d][d2] - (1.0/3.0) * ((d==d2)?1.0:0.0) * depsggij_rho ) + SiRijgg_rho[d][d2] + SiRijgg_rho[d2][d] + inter_n.Sijs[d][d2] * depsggij_rho;
    }
  }
}

//calculate acceleration of pressure force, internal energy change rate and change rate of deviatric stress tensor.
void calc_acceleration_dudt_and_dSab_rho_dt(particle_t par[])
{
  int i,d,n,d2;
  int particle_number=get_particle_number();
  double next_v[D];
  double dt=get_dt();
  double delta2;
  inter_t *inter;
  int error=0;
  particle_t *par_j;
  double t_center_v_const;
  if(is_riemann) t_center_v_const=0.5;
  else t_center_v_const=1.0;
  int inter_N_max=0;
  for( i=0 ; i<particle_number ; i++ ){
    if( inter_N_max < par[i].inter_N ) inter_N_max = par[i].inter_N;
  }
  int this_thread;
#ifdef _OPENMP
  this_thread = omp_get_thread_num();
#else
  this_thread = 0;
#endif

  //calclate a, dudt and dSabdt for Godunov SPH
  if(is_riemann){
#ifdef _OPENMP
#pragma omp parallel private(i,d,d2,next_v,n,inter,delta2,par_j)
#endif
    {
      //allocation of inter array
      inter = (inter_t*)malloc(sizeof(inter_t)*inter_N_max);
      if(inter == NULL){
	printf("cannot allocate inter array! \n");
	exit(FAILURE);
      }
      count_memory((int)(sizeof(inter_t)*inter_N_max));

#ifdef _OPENMP
#pragma omp for
#endif
      for( i=0 ; i<particle_number ; i++ ){

	//initialize acceleration, dudt and dSabdt
	for(d=0 ; d<D ; d++){
	  par[i].a[d] = 0.0;
	  for(d2=0 ; d2<D ; d2++){
	    par[i].dSab_rho_dt[d][d2] = 0.0;
	  }
	}
	par[i].dudt = 0.0;

	//calculate acceleration of pressure force 
	for(n=0 ; n<par[i].inter_N ; n++){
	  par_j = par[i].inter_par[n];
	  if( &par[i] != par_j ){
	    calc_aij_riemann(&par[i],*par_j,&(inter[n]),dt);
	  }
	}

	//if a[i] is not a number
	for(d=0 ; d<D ; d++){
	  if( isnan(par[i].a[d]) ){
	    printf("par %d's a[%d] is not a number! \n",i,d);
	    error=1;
	  }
	}

	if(is_ienergy_needed){
	  //calculate next timestep's velosity of i particle for energy conservation
	  for(d=0 ; d<D ; d++){
	    next_v[d] = par[i].v_old[d] + t_center_v_const * par[i].a[d] * dt;
	  }
      
	  //calculate dudt 
	  for(n=0 ; n<par[i].inter_N ; n++){
	    par_j = par[i].inter_par[n];
	    if( &par[i] != par_j){
	      par[i].dudt += calc_dudtij_riemann(par[i],*par_j,next_v,inter[n]);
	    }
	  }
      
	  //if dudt is not a number
	  if( isnan(par[i].dudt) ){
	    printf("par %d's dudt is not a number! \n",i);
	    error=1;
	  }
	}

	if(is_cohesion){
	  //calculate dSab_rho/dt
	  for(n=0 ; n<par[i].inter_N ; n++){
	    par_j = par[i].inter_par[n];
	    if( &par[i] != par_j){
	      calc_dSab_rho_dt_riemann(&par[i],*par_j,inter[n]);
	    }
	  }
	  
	  //if dSabrhodt is not a number
	  for(d=0 ; d<D ; d++){
	    for(d2=0 ; d2<D ; d2++){
	      if( isnan(par[i].dSab_rho_dt[d][d2]) ){
		printf("par %d's dSab_rho_dt[%d][%d] is not a number! \n",i,d,d2);
		error=1;
	      }
	    }
	  }
	}
	else{
	  for(d=0; d<D ; d++){
	    for(d2=0 ; d2<D ; d2++){
	      par[i].Sab[d][d2]=0.0;
	      par[i].Sab_rho[d][d2]=0.0;
	    }
	  }
	}

      }
      free(inter);
      count_memory((int)(-(sizeof(inter_t)*inter_N_max)));
    }
  }
  //calculate a, dudt and dSab_rho_dt for standard SPH
  else{
#ifdef _OPENMP
#pragma omp parallel private(i,d,d2,next_v,n,inter,delta2,par_j)
#endif
    {
      //allocation of inter array
      inter = (inter_t*)malloc(sizeof(inter_t)*inter_N_max);
      if(inter == NULL){
	printf("cannot allocate inter array! \n");
	exit(FAILURE);
      }
      count_memory((int)(sizeof(inter_t)*inter_N_max));

#ifdef _OPENMP
#pragma omp for
#endif
      for( i=0 ; i<particle_number ; i++ ){

	//initialize acceleration, dudt and dSab_rho_dt
	for(d=0 ; d<D ; d++){
	  par[i].a[d] = 0.0;
	  for(d2=0 ; d2<D ; d2++){
	    par[i].dSab_rho_dt[d][d2] = 0.0;
	  }
	}
	par[i].dudt = 0.0;

	//calculate acceleration of pressure force 
	for(n=0 ; n<par[i].inter_N ; n++){
	  par_j = par[i].inter_par[n];
	  if( &par[i] != par_j ){
	    calc_aij_standard(&par[i],*par_j,&(inter[n].common));
	  }
	}

	//if a[i] is not a number
	for(d=0 ; d<D ; d++){
	  if( isnan(par[i].a[d]) ){
	    printf("par %d's a[%d] is not a number! \n",i,d);
	    error=1;
	  }
	}

	if(is_ienergy_needed){
	  //calculate next timestep's velosity of i particle for energy conservation
	  for(d=0 ; d<D ; d++){
	    next_v[d] = par[i].v_old[d] + t_center_v_const * par[i].a[d] * dt;
	  }
      
	  //calculate dudt 
	  for(n=0 ; n<par[i].inter_N ; n++){
	    par_j = par[i].inter_par[n];
	    if( &par[i] != par_j){
	      par[i].dudt += calc_dudtij_standard(par[i],*par_j,next_v,inter[n].common);
	    }
	  }
      
	  //if dudt is not a number
	  if( isnan(par[i].dudt) ){
	    printf("par %d's dudt is not a number! \n",i);
	    error=1;
	  }
	}

	if(is_cohesion){
	  //calculate dSab_rho/dt
	  for(n=0 ; n<par[i].inter_N ; n++){
	    par_j = par[i].inter_par[n];
	    if( &par[i] != par_j){
	      calc_dSab_rho_dt_standard(&par[i],*par_j);
	    }
	  }
	  
	  //if dSab_rho_dt is not a number
	  for(d=0 ; d<D ; d++){
	    for(d2=0 ; d2<D ; d2++){
	      if( isnan(par[i].dSab_rho_dt[d][d2]) ){
		printf("par %d's dSabdt[%d][%d] is not a number! \n",i,d,d2);
		error=1;
	      }
	    }
	  }
	}
	else{
	  for(d=0; d<D ; d++){
	    for(d2=0 ; d2<D ; d2++){
	      par[i].Sab[d][d2]=0.0;
	      par[i].Sab_rho[d][d2]=0.0;
	    }
	  }
	}

      }
      free(inter);
      count_memory((int)(-(sizeof(inter_t)*inter_N_max)));
    }
  }
	
  if(this_thread==0&&error==1){
    printf("there are not a numbers\n");
    if(!DEBUG){
      char filename[]="debug-particleinfo.bin";
      //make_initial_condition_file(particle_number,D,get_t(),par,filename);
    }
    exit(FAILURE);
  }
}

//calculate plastic Sab
void calc_plastic_Sab(particle_t par[])
{
  int i,d,d2;
  int particle_number=get_particle_number();
  double J2;//second invariand of deviatoric stress tensor
  double f;
  double Y_yield,Yi,Yd;
  double P_pos,u_yield;
  
#ifdef _OPENMP
#pragma omp parallel for private(i,d,d2,J2,Y_yield,Yi,Yd,P_pos,u_yield,f)
#endif
  for(i=0 ; i<particle_number ; i++){
    J2=0.0;
    for(d=0 ; d<D ; d++){
      for(d2=0 ; d2<D ; d2++){
	par[i].Sab[d][d2] = par[i].Sab_rho[d][d2] * par[i].rho;
	J2 += 0.5 * par[i].Sab[d][d2] * par[i].Sab[d][d2];
      }
    }
    switch(plastic_model[par[i].property_tag]){
    case 0 :
      f=1;
      if(is_fracture_model[par[i].property_tag]) f *= ( 1.0 - par[i].damage );
      break;
    case 1 :
      //the way to describe elastic-perfectly plastic material in Benz and Asphaug 1995
      f = ( (Y0_mises[par[i].property_tag]*Y0_mises[par[i].property_tag]/(3.0*J2))<1.0 ) ? (Y0_mises[par[i].property_tag]*Y0_mises[par[i].property_tag]/(3.0*J2)) : 1.0 ;
      if(is_fracture_model[par[i].property_tag]) f *= ( 1.0 - par[i].damage );
      break;
    case 2 :
      //the way to describe elastic-perfectly plastic material in Libersky and Petschek 1990
      f = ( (sqrt(Y0_mises[par[i].property_tag]*Y0_mises[par[i].property_tag]/(3.0*J2)))<1.0 ) ? (sqrt(Y0_mises[par[i].property_tag]*Y0_mises[par[i].property_tag]/(3.0*J2))) : 1.0 ;
      if(is_fracture_model[par[i].property_tag]) f *= ( 1.0 - par[i].damage );
      break;
    case 3 :
      P_pos = ( ( par[i].p > 0 ) ? par[i].p : 0 );
      u_yield = ( ( par[i].u < u_melt[par[i].property_tag] ) ? par[i].u : u_melt[par[i].property_tag] ); 
      Yi = ( Y0_cohesion[par[i].property_tag] + ( mu_i_fric[par[i].property_tag] * P_pos ) / ( 1.0 + mu_i_fric[par[i].property_tag] * P_pos / ( Y0_mises[par[i].property_tag] - Y0_cohesion[par[i].property_tag] ) ) ) * ( 1.0 - u_yield / u_melt[par[i].property_tag] );
      Yd = mu_d_fric[par[i].property_tag] * P_pos;
      if( Yd > Yi ) Yd = Yi;
      if(!is_friction_model[par[i].property_tag]) Yd = 0.0;
      Y_yield = ( 1.0 - par[i].damage ) * Yi + par[i].damage * Yd;
      if( Y_yield > Yi || !is_fracture_model[par[i].property_tag] ) Y_yield = Yi;
      f = ( Y_yield/sqrt(J2) < 1.0 ) ? (Y_yield/sqrt(J2)) : 1.0 ;
      break;
    }
    if( isnan(f) ){
      f=0.0;
    }
    par[i].f_test=f;
    for(d=0 ; d<D ; d++){
      for(d2=0 ; d2<D ; d2++){
	par[i].Sab[d][d2] = ( par[i].Sab_rho[d][d2] * par[i].rho ) * f;
      }
    }
  }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//fanctions for neighbor searching//////////////////////////////////////////////////////////////////////////////////////////////////

void neighbor_searching_recursively(particle_t *par_i,tree_node_t *n_node)
{
  particle_t *p_list;
  double delta2;
  int d,n;
  
  //if n_node is leaf node
  if(n_node->is_leaf){
    for(p_list = n_node->first_particle ; p_list ; p_list = p_list->tree_next){
      delta2=0.0;
      for(d=0 ; d<D ; d++){
	delta2 += ( par_i->r[d] - p_list->r[d] ) * ( par_i->r[d] - p_list->r[d] );
      }
      if( sqrt(delta2) < 1.0001 * Ncell * par_i->Csmooth * ( (par_i->h > p_list->h)?par_i->h:p_list->h ) ){
	par_i->inter_par[par_i->inter_N] = p_list;
	par_i->inter_N++;
	//if the number of inter particles for i particle is larger than max number,increase nmax number and realloc inter array
	if( par_i->inter_N >= par_i->inter_Nmax){
	  par_i->inter_Nmax += INTER;
	  par_i->inter_par = (particle_t**)realloc(par_i->inter_par,sizeof(particle_t*)*par_i->inter_Nmax);
	  if(par_i->inter_par==NULL){
	    printf("cannot realloc inter particle array! \n");
	    exit(FAILURE);
	  }
	  count_memory((int)(sizeof(particle_t*)*INTER));
	}
      }
    }
  }
  //if n_node is not leaf node
  else{
    //if particle i's hmax is overlapped to n_node,let n_node's children seach neighbor particles recursively.
    if(is_overlap(par_i->r,n_node->center,n_node->box_size,par_i->hmax,par_i->Csmooth)){
      for(n=0 ; n<N_CHILD ; n++){
	if(n_node->child[n]){
	  neighbor_searching_recursively(par_i,n_node->child[n]);
	}
      }
    }
  }
}

void neighbor_searching(particle_t par[],tree_node_t *root_node)
{
  int i;
  int particle_number=get_particle_number();
  
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
  for(i=0;i<particle_number;i++){
    par[i].inter_N=0;
    neighbor_searching_recursively(&par[i],root_node);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//calculate density at any point that have position r and smoothing length h
double calc_rho_at_r(double r[],double h,tree_node_t *root_node)
{
  particle_t cell_par;
  particle_t *par_j;
  int d,n;
  double delta2;
  double rho=0.0;

  for( d=0 ; d<D ; d++){
    cell_par.r[d] = r[d];
  }
  cell_par.h = cell_par.hmax = h;
  cell_par.inter_N = 0;
  cell_par.inter_Nmax = INTER;
  cell_par.inter_par = (particle_t**)malloc(INTER*sizeof(particle_t*));
  if(cell_par.inter_par==NULL){
    printf("cannot allocate cell particle neighbor! \n");
    return(0.0);
  }
  count_memory((int)(INTER*sizeof(particle_t*)));

  neighbor_searching_recursively(&cell_par,root_node);

  for( n=0 ; n<cell_par.inter_N ; n++ ){
    par_j = cell_par.inter_par[n];
    delta2=0.0;
    for( d=0 ; d<D ; d++ ){
      delta2 += ( cell_par.r[d] - par_j->r[d] ) * ( cell_par.r[d] - par_j->r[d] );
    }
    rho += par_j->mass * pow( 1.0 / ( h * sqrt(PI) ), D ) * exp( -delta2 / ( h * h ) );
  }

  free(cell_par.inter_par);
  count_memory((int)(-(INTER*sizeof(particle_t*))));

  return(rho);
}

//set Csmooth for 2D and cubic spline case in negative pressure to suppress the tensile instability.
void set_Csmooth_for_2D_cubic(particle_t par[])
{
  int particle_number=get_particle_number();
  int i;
  double Csmooth_crit;
  double A=0.925887;
  double B=2.37512;
  double C=-0.89341;
  double val;

#ifdef _OPENMP
#pragma omp parallel for private(i,Csmooth_crit,val)
#endif
  for( i=0 ; i<particle_number ; i++){
    if(par[i].p>0.0){
      par[i].Csmooth = 1.0;
    }
    else{
      val = fabs( ( par[i].p / par[i].rho ) / ( par[i].Cs * par[i].Cs ) );
      Csmooth_crit = A * log( B * ( val - C ) );
      par[i].Csmooth = ( Csmooth_crit + 0.3 > 1.0 ) ? Csmooth_crit + 0.3 : 1.0 ;
    }
  }

}

