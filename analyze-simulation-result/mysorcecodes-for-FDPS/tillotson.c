#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SPH.h"
#include "param-mysorcecode.h"

//calculate tillotson pressure for compression or cold expansion state
double calc_ps(double E,double rho,int p_tag)
{
  double eta=rho/rho0_til[p_tag];
  double mu=eta-1.0;
  double ps;

  ps = ( a_til[p_tag] + b_til[p_tag]/((E/(E0_til[p_tag]*eta*eta))+1.0) )*rho*E+A_til[p_tag]*mu+B_til[p_tag]*mu*mu;

  return(ps);
}

//calculate tillotson pressure for hot expansion state
double calc_pg(double E,double rho,int p_tag)
{
  double eta=rho/rho0_til[p_tag];
  double mu=eta-1.0;
  double pg;

  pg = a_til[p_tag]*rho*E + ( b_til[p_tag]*rho*E/((E/(E0_til[p_tag]*eta*eta))+1.0) + A_til[p_tag]*mu*exp(-beta_til[p_tag]*(rho0_til[p_tag]/rho-1.0)) )*exp(-alpha_til[p_tag]*pow(rho0_til[p_tag]/rho-1.0,2.0));

  return(pg);
}

//calculate tillotson EOS's dp/drho for compression or cold expansion state
double calc_dps_drho(double E,double rho,double dE_drho,int p_tag)
{
  double eta=rho/rho0_til[p_tag];
  double mu=eta-1.0;
  double dps_drho;

  dps_drho = ( -b_til[p_tag]*(dE_drho/(E0_til[p_tag]*eta*eta)-2.0*E/(rho0_til[p_tag]*E0_til[p_tag]*pow(eta,3.0))) / (pow(E/(E0_til[p_tag]*eta*eta)+1.0,2.0)) )*rho*E + (a_til[p_tag]+b_til[p_tag]/(E/(E0_til[p_tag]*eta*eta)+1.0))*(E+rho*dE_drho) + A_til[p_tag]/rho0_til[p_tag] + (2.0*B_til[p_tag]/rho0_til[p_tag])*(rho/rho0_til[p_tag]-1.0);

  return(dps_drho);
}

//calculate tillotson EOS's dp/drho for hot expansion state
double calc_dpg_drho(double E,double rho,double dE_drho,int p_tag)
{
  double eta=rho/rho0_til[p_tag];
  double mu=eta-1.0;
  double dpg_drho;

  dpg_drho = a_til[p_tag]*(E+rho*dE_drho) + ( b_til[p_tag]*(E+rho*dE_drho)/(E/(E0_til[p_tag]*eta*eta)+1.0) - b_til[p_tag]*rho*E*(dE_drho/(E0_til[p_tag]*eta*eta)-2*E/(E0_til[p_tag]*rho0_til[p_tag]*pow(eta,3.0)))/pow(E/(E0_til[p_tag]*eta*eta)+1.0,2.0) + (A_til[p_tag]/rho0_til[p_tag])*exp(-beta_til[p_tag]*(rho0_til[p_tag]/rho-1.0)) + A_til[p_tag]*mu*beta_til[p_tag]*(rho0_til[p_tag]/(rho*rho))*exp(-beta_til[p_tag]*(rho0_til[p_tag]/rho-1.0)) )*exp(-alpha_til[p_tag]*pow(rho0_til[p_tag]/rho-1.0,2.0)) + ( b_til[p_tag]*rho*E/(E/(E0_til[p_tag]*eta*eta)+1.0) + A_til[p_tag]*mu*exp(-beta_til[p_tag]*(rho0_til[p_tag]/rho-1.0)) )*2.0*(rho0_til[p_tag]/rho-1.0)*alpha_til[p_tag]*(rho0_til[p_tag]/(rho*rho))*exp(-alpha_til[p_tag]*pow(rho0_til[p_tag]/rho-1.0,2.0));

  return(dpg_drho);
}

//calculate tillotson pressure
double calc_tillotson_pressure(double E,double rho,int p_tag)
{
  double p;

  //for compuression or cold expansion state
  if( rho > rho0_til[p_tag] || E < Eiv_til[p_tag] ){
    p = calc_ps(E,rho,p_tag);
  }
  //for hot expansion state
  else if( E > Ecv_til[p_tag] ){
    p = calc_pg(E,rho,p_tag);
  }
  //for intermediate state
  else{
    p = ( ( E - Eiv_til[p_tag] ) * calc_pg(E,rho,p_tag) + ( Ecv_til[p_tag] - E ) * calc_ps(E,rho,p_tag) )/( Ecv_til[p_tag] - Eiv_til[p_tag] );
  }

  return(p);
}

//calculate tillotson gamma
double calc_tillotson_gamma(double E,double rho,double p,int p_tag)
{
  double pg,ps;
  double dE_drho=p/(rho*rho);
  double dp_drho;
  double gamma;

  //for compression or cold expansion state
  if( rho > rho0_til[p_tag] || E < Eiv_til[p_tag] ){
    dp_drho = calc_dps_drho(E,rho,dE_drho,p_tag);
  }
  //for hot expansion state
  else if( E > Ecv_til[p_tag] ){
    dp_drho = calc_dpg_drho(E,rho,dE_drho,p_tag);
  }
  //for intermediate state
  else{
    pg = calc_pg(E,rho,p_tag);
    ps = calc_ps(E,rho,p_tag);
    dp_drho = ( ( E - Eiv_til[p_tag] ) * calc_dpg_drho(E,rho,dE_drho,p_tag) + ( Ecv_til[p_tag] - E ) * calc_dps_drho(E,rho,dE_drho,p_tag) + ( pg - ps ) * dE_drho )/( Ecv_til[p_tag] - Eiv_til[p_tag] );
  }

  gamma = ( rho / p ) * dp_drho;
  
  return(gamma);
}

//calculate dp/drho for tillotson EoS (=square of sound speed)
double calc_dp_drho_tillotson(double E,double rho,double p, int p_tag)
{
  double pg,ps;
  double dE_drho=p/(rho*rho);
  double dp_drho;
  double gamma;

  //for compression or cold expansion state
  if( rho > rho0_til[p_tag] || E < Eiv_til[p_tag] ){
    dp_drho = calc_dps_drho(E,rho,dE_drho,p_tag);
  }
  //for hot expansion state
  else if( E > Ecv_til[p_tag] ){
    dp_drho = calc_dpg_drho(E,rho,dE_drho,p_tag);
  }
  //for intermediate state
  else{
    pg = calc_pg(E,rho,p_tag);
    ps = calc_ps(E,rho,p_tag);
    dp_drho = ( ( E - Eiv_til[p_tag] ) * calc_dpg_drho(E,rho,dE_drho,p_tag) + ( Ecv_til[p_tag] - E ) * calc_dps_drho(E,rho,dE_drho,p_tag) + ( pg - ps ) * dE_drho )/( Ecv_til[p_tag] - Eiv_til[p_tag] );
  }

  return(dp_drho);
}

//calculate delp_delrho and delp_delu for compression or cold expansion state
void calc_delps_delrho_and_delps_delu(double E,double rho,double* delp_delrho_P,double* delp_delu_P,int p_tag)
{
  double rho0=rho0_til[p_tag];
  double A=A_til[p_tag];
  double B=B_til[p_tag];
  double E0=E0_til[p_tag];
  double Eiv=Eiv_til[p_tag];
  double Ecv=Ecv_til[p_tag];
  double a=a_til[p_tag];
  double b=b_til[p_tag];
  double alpha=alpha_til[p_tag];
  double beta=beta_til[p_tag];

  double eta = rho / rho0;
  double mu = eta - 1.0;

  *delp_delrho_P = ( (2.0*b/pow((E/(E0*eta*eta))+1.0,2.0))*(E*rho0*rho0/(E0*pow(rho,3.0))) )*rho*E + ( a + b/((E/(E0*eta*eta))+1.0) )*E + (A/rho0) + (2.0*B/rho0)*( (rho/rho0) - 1.0 );
  *delp_delu_P = ( (-b/pow((E/(E0*eta*eta))+1.0,2.0))*(1.0/(E0*eta*eta)) )*rho*E + ( a + b/((E/(E0*eta*eta))+1.0) )*rho;

}

//calculate delp_delrho and delp_delu for hot expansion state
void calc_delpg_delrho_and_delpg_delu(double E,double rho,double* delp_delrho_P,double* delp_delu_P,int p_tag)
{
  double rho0=rho0_til[p_tag];
  double A=A_til[p_tag];
  double B=B_til[p_tag];
  double E0=E0_til[p_tag];
  double Eiv=Eiv_til[p_tag];
  double Ecv=Ecv_til[p_tag];
  double a=a_til[p_tag];
  double b=b_til[p_tag];
  double alpha=alpha_til[p_tag];
  double beta=beta_til[p_tag];

  double eta = rho / rho0;
  double mu = eta - 1.0;
  
  *delp_delrho_P = a*E + ( b*E/((E/(E0*eta*eta))+1.0) + (2.0*b*rho*E/pow((E/(E0*eta*eta))+1.0,2.0))*(E*rho0*rho0/(E0*pow(rho,3.0))) + (A/rho0)*exp(-beta*((rho0/rho)-1.0)) + A*mu*((rho0*beta)/(rho*rho))*exp(-beta*((rho0/rho)-1.0)) )*exp(-alpha*pow((rho0/rho)-1.0,2.0)) + ( b*rho*E/((E/(E0*eta*eta))+1.0) + A*mu*exp(-beta*((rho0/rho)-1.0)) )*2.0*((rho0/rho)-1.0)*alpha*(rho0/(rho*rho))*exp(-alpha*pow((rho0/rho)-1.0,2.0));
  *delp_delu_P = a*rho + ( b*rho/((E/(E0*eta*eta))+1.0) - (b*rho*E/pow((E/(E0*eta*eta))+1.0,2.0))*(1.0/(E0*eta*eta)) )*exp(-alpha*pow((rho0/rho)-1.0,2.0));

}

//calculate delp_delrho and delp_delu for tillotson EoS
void calc_delp_delrho_and_delp_delu_tillotson(double E,double rho,int p_tag,double* delp_delrho_P,double* delp_delu_P)
{
  double pg,ps;
  double delps_delrho,delps_delu;
  double delpg_delrho,delpg_delu;

  //for compression or cold expansion state
  if( rho > rho0_til[p_tag] || E < Eiv_til[p_tag] ){
    calc_delps_delrho_and_delps_delu(E,rho,delp_delrho_P,delp_delu_P,p_tag);
  }
  //for hot expansion state
  else if( E > Ecv_til[p_tag] ){
    calc_delpg_delrho_and_delpg_delu(E,rho,delp_delrho_P,delp_delu_P,p_tag);
  }
  //for intermediate state
  else{
    pg = calc_pg(E,rho,p_tag);
    ps = calc_ps(E,rho,p_tag);
    calc_delps_delrho_and_delps_delu(E,rho,&delps_delrho,&delps_delu,p_tag);
    calc_delpg_delrho_and_delpg_delu(E,rho,&delpg_delrho,&delpg_delu,p_tag);
    *delp_delrho_P = ( ( E - Eiv_til[p_tag] ) * delps_delrho + ( Ecv_til[p_tag] - E ) * delpg_delrho )/( Ecv_til[p_tag] - Eiv_til[p_tag] );
    *delp_delu_P = ( ( E - Eiv_til[p_tag] ) * delps_delu + ( Ecv_til[p_tag] - E ) * delpg_delu + ( pg - ps ) )/( Ecv_til[p_tag] - Eiv_til[p_tag] );
  }

}

double calc_tillotson_E_from_p_and_rho(double rho,double p,int p_tag)
{
  double alpha2,beta2,gamma2,x,y;
  double E,El,Er,Emid,ptmp;
  double psclit,pgclit,acc;
  double eta = rho / rho0_til[p_tag];
  double mu = eta - 1.0;

  psclit = calc_ps(rho,Eiv_til[p_tag],p_tag);
  pgclit = calc_pg(rho,Ecv_til[p_tag],p_tag);

  //for compression or cold expansion state
  if(rho > rho0_til[p_tag] || p < psclit){
    alpha2 = a_til[p_tag] * rho / ( E0_til[p_tag] * eta * eta );
    beta2 = ( a_til[p_tag] + b_til[p_tag] ) * rho + ( A_til[p_tag] * mu + B_til[p_tag] * mu * mu - p ) / ( E0_til[p_tag] * eta * eta );
    gamma2 = A_til[p_tag] * mu + B_til[p_tag] * mu * mu - p;
  
    E = ( -beta2 + sqrt( beta2 * beta2 - 4.0 * alpha2 * gamma2 ) ) / ( 2.0 * alpha2 );
  }
  //for hot expansion state
  else if( p > pgclit ){
    x = exp( -beta_til[p_tag] * ( rho0_til[p_tag] / rho - 1.0 ) );
    y = exp( -alpha_til[p_tag] * pow( rho0_til[p_tag] / rho - 1.0 , 2.0 ) );
    alpha2 = a_til[p_tag] * rho / ( E0_til[p_tag] * eta * eta );
    beta2 = ( a_til[p_tag] + b_til[p_tag] * y ) * rho + ( A_til[p_tag] * mu * x * y - p ) / ( E0_til[p_tag] * eta * eta );
    gamma2 = A_til[p_tag] * mu * x * y - p;

    E = ( -beta2 + sqrt( beta2 * beta2 - 4.0 * alpha2 * gamma2 ) ) / (2.0 * alpha2 );
  }
  //for intermediate state
  else{
    acc = fabs( 1.0E-4 * 0.5 * ( psclit + pgclit ) );
    El = Eiv_til[p_tag];
    Er = Ecv_til[p_tag];
    Emid = 0.5 * ( Eiv_til[p_tag] + Ecv_til[p_tag] );
    ptmp = calc_tillotson_pressure(Emid,rho,p_tag);
    do{
      if( ptmp > p ){
	Er = Emid;
      }
      else{
	El = Emid;
      }
      Emid = 0.5 * ( Er + El );
      ptmp = calc_tillotson_pressure(Emid,rho,p_tag);
    }while( fabs( p - ptmp ) > acc );

    E = Emid;
  }

  if( E != E ){
    printf("E is not a number! \n");
    exit(FAILURE);
  }
  if(E < 0.0 ){
    printf("E is negative! \n");
    exit(FAILURE);
  }

  return(E);
}


