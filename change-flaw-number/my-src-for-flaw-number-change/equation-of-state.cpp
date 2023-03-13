#include "header.h"

//calculate pressure,sound speed using ideal gas eos.
void calc_p_and_Cs_ideal_gas(RealPtcl* par_i)
{
  PS::F64 dens_i;
  if(PARAM::IS_POROSITY_MODEL[par_i->property_tag]){
    dens_i = par_i->alpha_por * par_i->dens;
  }
  else{
    dens_i = par_i->dens;
  }

  par_i->pres = ( PARAM::GAMMA_IDEAL[par_i->property_tag] - 1.0 ) * dens_i * par_i->eng;
  if(PARAM::IS_POROSITY_MODEL[par_i->property_tag]){
    par_i->pres /= par_i->alpha_por;
  }
  par_i->gamma = PARAM::GAMMA_IDEAL[par_i->property_tag];
  par_i->snds = sqrt( PARAM::GAMMA_IDEAL[par_i->property_tag] * par_i->pres  / dens_i );
}

//calculate delp_delrho and delp_delu for ideal gas
void calc_delp_delrho_and_delp_delu_ideal_gas(PS::F64 rho,PS::F64 u,PS::S32 p_tag,PS::F64* delp_delrho_P,PS::F64* delp_delu_P)
{
  *delp_delrho_P = ( PARAM::GAMMA_IDEAL[p_tag] - 1.0 ) * u;
  *delp_delu_P = ( PARAM::GAMMA_IDEAL[p_tag] - 1.0 ) * rho;
}

//calculate pressure using elastic equation of state (P=Cs^2(rho-rho0eos)).
void calc_p_elastic_eos(RealPtcl* par_i)
{
  PS::F64 dens_i;
  if(PARAM::IS_POROSITY_MODEL[par_i->property_tag]){
    dens_i = par_i->alpha_por * par_i->dens;
  }
  else{
    dens_i = par_i->dens;
  }

  par_i->pres = PARAM::CS_ELASTIC[par_i->property_tag]*PARAM::CS_ELASTIC[par_i->property_tag]*(dens_i-PARAM::RHO0_ELASTIC[par_i->property_tag]);
  if(PARAM::IS_POROSITY_MODEL[par_i->property_tag]){
    par_i->pres /= par_i->alpha_por;
  }
  par_i->snds = PARAM::CS_ELASTIC[par_i->property_tag];
  par_i->gamma = 0.0;
}

//calculate delp_delrho and delp_delu for elastic EoS
void calc_delp_delrho_and_delp_delu_elastic_eos(PS::F64 rho,PS::F64 u,PS::S32 p_tag,PS::F64* delp_delrho_P,PS::F64* delp_delu_P)
{
  *delp_delrho_P = PARAM::CS_ELASTIC[p_tag]*PARAM::CS_ELASTIC[p_tag];
  *delp_delu_P = 0.0;
}

//calculate pressure, Cs and gamma of stiffened gas equation of state
void calc_p_Cs_and_gamma_stiffened_gas_EoS(RealPtcl* par_i)
{
  PS::F64 C0 = PARAM::C0_STIFFENED[par_i->property_tag];
  PS::F64 gamma0 = PARAM::GAMMA0_STIFFENED[par_i->property_tag];
  PS::F64 rho0 = PARAM::RHO0_STIFFENED[par_i->property_tag];
  PS::S32 sign;

  PS::F64 dens_i;
  if(PARAM::IS_POROSITY_MODEL[par_i->property_tag]){
    dens_i = par_i->alpha_por * par_i->dens;
  }
  else{
    dens_i = par_i->dens;
  }

  par_i->pres = C0 * C0 * ( dens_i - rho0 ) + ( gamma0 - 1.0 ) * dens_i * par_i->eng;
  if(PARAM::IS_POROSITY_MODEL[par_i->property_tag]){
    par_i->pres /= par_i->alpha_por;
  }
  par_i->snds = sqrt( C0 * C0 + ( gamma0 - 1.0 ) * ( par_i->eng + par_i->pres / dens_i ) );
  if( std::isnan(par_i->snds) ){
    par_i->snds = C0;
  }
  if( math::sign(par_i->pres) ) sign = -1;
  else                       sign =  1;
  par_i->gamma = ( dens_i / par_i->pres ) * par_i->snds * par_i->snds;
  if( fabs( par_i->gamma ) > 100.0 || std::isnan( par_i->gamma ) || std::isinf( par_i->gamma ) ){
    par_i->gamma = sign * 100;
  }
}

//calculate delp_delrho and delp_delu for stiffened gas EoS
void calc_delp_delrho_and_delp_delu_stiffened_gas_EoS(PS::F64 rho,PS::F64 u,PS::S32 p_tag,PS::F64* delp_delrho_P,PS::F64* delp_delu_P)
{
  PS::F64 C0 = PARAM::C0_STIFFENED[p_tag];
  PS::F64 gamma0 = PARAM::GAMMA0_STIFFENED[p_tag];
  PS::F64 rho0 = PARAM::RHO0_STIFFENED[p_tag];

  *delp_delrho_P = C0*C0 + ( gamma0 - 1.0 )*u;
  *delp_delu_P = ( gamma0 - 1.0 )*rho;
}

//calculate delp_delrho and delp_delu for Mie Gruneisen EoS
void calc_delp_delrho_and_delp_delu_Mie_Gruneisen_EoS(PS::F64 rho,PS::F64 u,PS::S32 p_tag,PS::F64* delp_delrho_P,PS::F64* delp_delu_P)
{
  PS::F64 rho0=PARAM::RHO0_MG[p_tag];
  PS::F64 C0=PARAM::C0_MG[p_tag];
  PS::F64 S=PARAM::S_MG[p_tag];
  PS::F64 Gamma=PARAM::GAMMA_MG[p_tag];
  PS::F64 eta=1.0-(rho0/rho);

  *delp_delrho_P = pow(rho0*C0/((1.0-S*eta)*rho),2.0)*( (1.0+(2.0*eta*S/(1.0-S*eta)))*(1.0-Gamma*eta/2.0) - Gamma*eta/2.0 );
  *delp_delu_P = Gamma*rho0;
}

//calculate pressure, Cs and gamma of Mie Gruneisen EoS
void calc_p_Cs_and_gamma_Mie_Gruneisen_EoS(RealPtcl* par_i)
{
  PS::F64 rho0=PARAM::RHO0_MG[par_i->property_tag];
  PS::F64 C0=PARAM::C0_MG[par_i->property_tag];
  PS::F64 S=PARAM::S_MG[par_i->property_tag];
  PS::F64 Gamma=PARAM::GAMMA_MG[par_i->property_tag];
  PS::F64 eta,delp_delrho,delp_delu;
  PS::S32 sign;

  PS::F64 dens_i;
  if(PARAM::IS_POROSITY_MODEL[par_i->property_tag]){
    dens_i = par_i->alpha_por * par_i->dens;
  }
  else{
    dens_i = par_i->dens;
  }
  
  eta = 1.0 - rho0 / dens_i;
  calc_delp_delrho_and_delp_delu_Mie_Gruneisen_EoS(dens_i,par_i->eng,par_i->property_tag,&delp_delrho,&delp_delu);
  
  par_i->pres = ( rho0 * C0 * C0 * eta / pow(1.0-S*eta,2.0) ) * ( 1.0 - 0.5 * Gamma * eta ) + Gamma * rho0 * par_i->eng;
  if(PARAM::IS_POROSITY_MODEL[par_i->property_tag]){
    par_i->pres /= par_i->alpha_por;
  }
  par_i->snds = sqrt( delp_delrho + (par_i->pres/(dens_i*dens_i))*delp_delu);
  if( std::isnan(par_i->snds) ){
    par_i->snds = C0;
  }
  if( math::sign(par_i->pres) ) sign = -1;
  else                       sign =  1;
  par_i->gamma = ( dens_i / par_i->pres ) * par_i->snds * par_i->snds;
  if( fabs( par_i->gamma ) > 100.0 || std::isnan( par_i->gamma ) || std::isinf( par_i->gamma ) ){
    par_i->gamma = sign * 100;
  }
}

//calculate tillotson pressure for compression or cold expansion state
PS::F64 calc_ps_tillotson(PS::F64 E,PS::F64 rho,PS::S32 p_tag)
{
  PS::F64 rho0_til=PARAM::RHO0_TIL[p_tag];
  PS::F64 A_til=PARAM::A_TIL[p_tag];
  PS::F64 B_til=PARAM::B_TIL[p_tag];
  PS::F64 E0_til=PARAM::E0_TIL[p_tag];
  PS::F64 Eiv_til=PARAM::EIV_TIL[p_tag];
  PS::F64 Ecv_til=PARAM::ECV_TIL[p_tag];
  PS::F64 a_til=PARAM::a_TIL[p_tag];
  PS::F64 b_til=PARAM::b_TIL[p_tag];
  PS::F64 alpha_til=PARAM::ALPHA_TIL[p_tag];
  PS::F64 beta_til=PARAM::BETA_TIL[p_tag];

  PS::F64 eta=rho/rho0_til;
  PS::F64 mu=eta-1.0;
  PS::F64 ps;

  ps = ( a_til + b_til/((E/(E0_til*eta*eta))+1.0) )*rho*E+A_til*mu+B_til*mu*mu;

  return(ps);
}

//calculate tillotson pressure for hot expansion state
PS::F64 calc_pg_tillotson(PS::F64 E,PS::F64 rho,PS::S32 p_tag)
{
  PS::F64 rho0_til=PARAM::RHO0_TIL[p_tag];
  PS::F64 A_til=PARAM::A_TIL[p_tag];
  PS::F64 B_til=PARAM::B_TIL[p_tag];
  PS::F64 E0_til=PARAM::E0_TIL[p_tag];
  PS::F64 Eiv_til=PARAM::EIV_TIL[p_tag];
  PS::F64 Ecv_til=PARAM::ECV_TIL[p_tag];
  PS::F64 a_til=PARAM::a_TIL[p_tag];
  PS::F64 b_til=PARAM::b_TIL[p_tag];
  PS::F64 alpha_til=PARAM::ALPHA_TIL[p_tag];
  PS::F64 beta_til=PARAM::BETA_TIL[p_tag];

  PS::F64 eta=rho/rho0_til;
  PS::F64 mu=eta-1.0;
  PS::F64 pg;

  pg = a_til*rho*E + ( b_til*rho*E/((E/(E0_til*eta*eta))+1.0) + A_til*mu*exp(-beta_til*(rho0_til/rho-1.0)) )*exp(-alpha_til*pow(rho0_til/rho-1.0,2.0));

  return(pg);
}

//calculate tillotson EOS's dp/drho for compression or cold expansion state
PS::F64 calc_dps_drho_tillotson(PS::F64 E,PS::F64 rho,PS::F64 dE_drho,PS::S32 p_tag)
{
  PS::F64 rho0_til=PARAM::RHO0_TIL[p_tag];
  PS::F64 A_til=PARAM::A_TIL[p_tag];
  PS::F64 B_til=PARAM::B_TIL[p_tag];
  PS::F64 E0_til=PARAM::E0_TIL[p_tag];
  PS::F64 Eiv_til=PARAM::EIV_TIL[p_tag];
  PS::F64 Ecv_til=PARAM::ECV_TIL[p_tag];
  PS::F64 a_til=PARAM::a_TIL[p_tag];
  PS::F64 b_til=PARAM::b_TIL[p_tag];
  PS::F64 alpha_til=PARAM::ALPHA_TIL[p_tag];
  PS::F64 beta_til=PARAM::BETA_TIL[p_tag];

  PS::F64 eta=rho/rho0_til;
  PS::F64 mu=eta-1.0;
  PS::F64 dps_drho;

  dps_drho = ( -b_til*(dE_drho/(E0_til*eta*eta)-2.0*E/(rho0_til*E0_til*pow(eta,3.0))) / (pow(E/(E0_til*eta*eta)+1.0,2.0)) )*rho*E + (a_til+b_til/(E/(E0_til*eta*eta)+1.0))*(E+rho*dE_drho) + A_til/rho0_til + (2.0*B_til/rho0_til)*(rho/rho0_til-1.0);

  return(dps_drho);
}

//calculate tillotson EOS's dp/drho for hot expansion state
PS::F64 calc_dpg_drho_tillotson(PS::F64 E,PS::F64 rho,PS::F64 dE_drho,PS::S32 p_tag)
{
  PS::F64 rho0_til=PARAM::RHO0_TIL[p_tag];
  PS::F64 A_til=PARAM::A_TIL[p_tag];
  PS::F64 B_til=PARAM::B_TIL[p_tag];
  PS::F64 E0_til=PARAM::E0_TIL[p_tag];
  PS::F64 Eiv_til=PARAM::EIV_TIL[p_tag];
  PS::F64 Ecv_til=PARAM::ECV_TIL[p_tag];
  PS::F64 a_til=PARAM::a_TIL[p_tag];
  PS::F64 b_til=PARAM::b_TIL[p_tag];
  PS::F64 alpha_til=PARAM::ALPHA_TIL[p_tag];
  PS::F64 beta_til=PARAM::BETA_TIL[p_tag];

  PS::F64 eta=rho/rho0_til;
  PS::F64 mu=eta-1.0;
  PS::F64 dpg_drho;

  dpg_drho = a_til*(E+rho*dE_drho) + ( b_til*(E+rho*dE_drho)/(E/(E0_til*eta*eta)+1.0) - b_til*rho*E*(dE_drho/(E0_til*eta*eta)-2*E/(E0_til*rho0_til*pow(eta,3.0)))/pow(E/(E0_til*eta*eta)+1.0,2.0) + (A_til/rho0_til)*exp(-beta_til*(rho0_til/rho-1.0)) + A_til*mu*beta_til*(rho0_til/(rho*rho))*exp(-beta_til*(rho0_til/rho-1.0)) )*exp(-alpha_til*pow(rho0_til/rho-1.0,2.0)) + ( b_til*rho*E/(E/(E0_til*eta*eta)+1.0) + A_til*mu*exp(-beta_til*(rho0_til/rho-1.0)) )*2.0*(rho0_til/rho-1.0)*alpha_til*(rho0_til/(rho*rho))*exp(-alpha_til*pow(rho0_til/rho-1.0,2.0));

  return(dpg_drho);
}

//calculate tillotson gamma and sound speed
void calc_tillotson_gamma_and_Cs(PS::F64 E,PS::F64 rho,PS::F64 p,PS::S32 p_tag,PS::F64* gamma_P,PS::F64* Cs_P)
{
  PS::F64 pg,ps;
  PS::F64 dE_drho=p/(rho*rho);
  PS::F64 dp_drho;
  PS::F64 gamma;
  PS::F64 Eiv_til=PARAM::EIV_TIL[p_tag];
  PS::F64 Ecv_til=PARAM::ECV_TIL[p_tag];
  PS::F64 rho0_til=PARAM::RHO0_TIL[p_tag];
  PS::F64 A_til=PARAM::A_TIL[p_tag];
  PS::S32 sign;

  //for compression or cold expansion state
  if( rho > rho0_til || E < Eiv_til ){
    dp_drho = calc_dps_drho_tillotson(E,rho,dE_drho,p_tag);
  }
  //for hot expansion state
  else if( E > Ecv_til ){
    dp_drho = calc_dpg_drho_tillotson(E,rho,dE_drho,p_tag);
  }
  //for intermediate state
  else{
    pg = calc_pg_tillotson(E,rho,p_tag);
    ps = calc_ps_tillotson(E,rho,p_tag);
    dp_drho = ( ( E - Eiv_til ) * calc_dpg_drho_tillotson(E,rho,dE_drho,p_tag) + ( Ecv_til - E ) * calc_dps_drho_tillotson(E,rho,dE_drho,p_tag) + ( pg - ps ) * dE_drho )/( Ecv_til - Eiv_til );
  }

  *Cs_P = sqrt( dp_drho );
  if( std::isnan(*Cs_P) ){
    *Cs_P = sqrt( A_til / rho0_til );
  }
  if( math::sign(p) ) sign = -1;
  else             sign =  1;
  *gamma_P = ( rho / p ) * (*Cs_P) * (*Cs_P);
  if( fabs( *gamma_P ) > 100.0 || std::isnan( *gamma_P ) || std::isinf( *gamma_P ) ){
    *gamma_P = sign * 100;
  }
}

//calculate tillotson pressure, gamma and sound speed
void calc_p_Cs_and_gamma_tillotson(RealPtcl* par_i)
{
  PS::F64 p;
  PS::F64 Eiv_til=PARAM::EIV_TIL[par_i->property_tag];
  PS::F64 Ecv_til=PARAM::ECV_TIL[par_i->property_tag];
  PS::F64 rho0_til=PARAM::RHO0_TIL[par_i->property_tag];
  PS::F64 gamma_temp,Cs_temp;

  PS::F64 dens_i;
  if(PARAM::IS_POROSITY_MODEL[par_i->property_tag]){
    dens_i = par_i->alpha_por * par_i->dens;
  }
  else{
    dens_i = par_i->dens;
  }

  //for compuression or cold expansion state
  if( dens_i > rho0_til || par_i->eng < Eiv_til ){
    p = calc_ps_tillotson(par_i->eng,dens_i,par_i->property_tag);
  }
  //for hot expansion state
  else if( par_i->eng > Ecv_til ){
    p = calc_pg_tillotson(par_i->eng,dens_i,par_i->property_tag);
  }
  //for intermediate state
  else{
    p = ( ( par_i->eng - Eiv_til ) * calc_pg_tillotson(par_i->eng,dens_i,par_i->property_tag) + ( Ecv_til - par_i->eng ) * calc_ps_tillotson(par_i->eng,dens_i,par_i->property_tag) )/( Ecv_til - Eiv_til );
  }
  par_i->pres = p;
  if(PARAM::IS_POROSITY_MODEL[par_i->property_tag]){
    par_i->pres /= par_i->alpha_por;
  }

  calc_tillotson_gamma_and_Cs(par_i->eng,dens_i,par_i->pres,par_i->property_tag,&gamma_temp,&Cs_temp);
  par_i->snds = Cs_temp;
  par_i->gamma = gamma_temp;
}

//calculate delp_delrho and delp_delu for compression or cold expansion state
void calc_delps_delrho_and_delps_delu(PS::F64 E,PS::F64 rho,PS::F64* delp_delrho_P,PS::F64* delp_delu_P,PS::S32 p_tag)
{
  PS::F64 rho0=PARAM::RHO0_TIL[p_tag];
  PS::F64 A=PARAM::A_TIL[p_tag];
  PS::F64 B=PARAM::B_TIL[p_tag];
  PS::F64 E0=PARAM::E0_TIL[p_tag];
  PS::F64 Eiv=PARAM::EIV_TIL[p_tag];
  PS::F64 Ecv=PARAM::ECV_TIL[p_tag];
  PS::F64 a=PARAM::a_TIL[p_tag];
  PS::F64 b=PARAM::b_TIL[p_tag];
  PS::F64 alpha=PARAM::ALPHA_TIL[p_tag];
  PS::F64 beta=PARAM::BETA_TIL[p_tag];

  PS::F64 eta = rho / rho0;
  PS::F64 mu = eta - 1.0;

  *delp_delrho_P = ( (2.0*b/pow((E/(E0*eta*eta))+1.0,2.0))*(E*rho0*rho0/(E0*pow(rho,3.0))) )*rho*E + ( a + b/((E/(E0*eta*eta))+1.0) )*E + (A/rho0) + (2.0*B/rho0)*( (rho/rho0) - 1.0 );
  *delp_delu_P = ( (-b/pow((E/(E0*eta*eta))+1.0,2.0))*(1.0/(E0*eta*eta)) )*rho*E + ( a + b/((E/(E0*eta*eta))+1.0) )*rho;

}

//calculate delp_delrho and delp_delu for hot expansion state
void calc_delpg_delrho_and_delpg_delu(PS::F64 E,PS::F64 rho,PS::F64* delp_delrho_P,PS::F64* delp_delu_P,PS::S32 p_tag)
{
  PS::F64 rho0=PARAM::RHO0_TIL[p_tag];
  PS::F64 A=PARAM::A_TIL[p_tag];
  PS::F64 B=PARAM::B_TIL[p_tag];
  PS::F64 E0=PARAM::E0_TIL[p_tag];
  PS::F64 Eiv=PARAM::EIV_TIL[p_tag];
  PS::F64 Ecv=PARAM::ECV_TIL[p_tag];
  PS::F64 a=PARAM::a_TIL[p_tag];
  PS::F64 b=PARAM::b_TIL[p_tag];
  PS::F64 alpha=PARAM::ALPHA_TIL[p_tag];
  PS::F64 beta=PARAM::BETA_TIL[p_tag];  

  PS::F64 eta = rho / rho0;
  PS::F64 mu = eta - 1.0;
  
  *delp_delrho_P = a*E + ( b*E/((E/(E0*eta*eta))+1.0) + (2.0*b*rho*E/pow((E/(E0*eta*eta))+1.0,2.0))*(E*rho0*rho0/(E0*pow(rho,3.0))) + (A/rho0)*exp(-beta*((rho0/rho)-1.0)) + A*mu*((rho0*beta)/(rho*rho))*exp(-beta*((rho0/rho)-1.0)) )*exp(-alpha*pow((rho0/rho)-1.0,2.0)) + ( b*rho*E/((E/(E0*eta*eta))+1.0) + A*mu*exp(-beta*((rho0/rho)-1.0)) )*2.0*((rho0/rho)-1.0)*alpha*(rho0/(rho*rho))*exp(-alpha*pow((rho0/rho)-1.0,2.0));
  *delp_delu_P = a*rho + ( b*rho/((E/(E0*eta*eta))+1.0) - (b*rho*E/pow((E/(E0*eta*eta))+1.0,2.0))*(1.0/(E0*eta*eta)) )*exp(-alpha*pow((rho0/rho)-1.0,2.0));

}

//calculate delp_delrho and delp_delu for tillotson EoS
void calc_delp_delrho_and_delp_delu_tillotson(PS::F64 E,PS::F64 rho,PS::S32 p_tag,PS::F64* delp_delrho_P,PS::F64* delp_delu_P)
{
  PS::F64 pg,ps;
  PS::F64 delps_delrho,delps_delu;
  PS::F64 delpg_delrho,delpg_delu;

  //for compression or cold expansion state
  if( rho > PARAM::RHO0_TIL[p_tag] || E < PARAM::EIV_TIL[p_tag] ){
    calc_delps_delrho_and_delps_delu(E,rho,delp_delrho_P,delp_delu_P,p_tag);
  }
  //for hot expansion state
  else if( E > PARAM::ECV_TIL[p_tag] ){
    calc_delpg_delrho_and_delpg_delu(E,rho,delp_delrho_P,delp_delu_P,p_tag);
  }
  //for intermediate state
  else{
    pg = calc_pg_tillotson(E,rho,p_tag);
    ps = calc_ps_tillotson(E,rho,p_tag);
    calc_delps_delrho_and_delps_delu(E,rho,&delps_delrho,&delps_delu,p_tag);
    calc_delpg_delrho_and_delpg_delu(E,rho,&delpg_delrho,&delpg_delu,p_tag);
    *delp_delrho_P = ( ( E - PARAM::EIV_TIL[p_tag] ) * delps_delrho + ( PARAM::ECV_TIL[p_tag] - E ) * delpg_delrho )/( PARAM::ECV_TIL[p_tag] - PARAM::EIV_TIL[p_tag] );
    *delp_delu_P = ( ( E - PARAM::EIV_TIL[p_tag] ) * delps_delu + ( PARAM::ECV_TIL[p_tag] - E ) * delpg_delu + ( pg - ps ) )/( PARAM::ECV_TIL[p_tag] - PARAM::EIV_TIL[p_tag] );
  }

}

//calculate pressure, gamma and sound speed. 
//we select appropriate EoS for each particle depending on param.h
void calc_pressures(PS::ParticleSystem<RealPtcl>& sph_system)
{
  PS::S64 i;
  PS::S64 particle_number_local = sph_system.getNumberOfParticleLocal();

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for private(i) shared(particle_number_local)
#endif	
  for(i=0 ; i<particle_number_local ; i++){
    switch(PARAM::EQUATION_OF_STATE[sph_system[i].property_tag]){
    case 0 : calc_p_and_Cs_ideal_gas(&(sph_system[i])); break;
    case 1 : calc_p_elastic_eos(&(sph_system[i])); break;
    case 2 : calc_p_Cs_and_gamma_stiffened_gas_EoS(&(sph_system[i])); break;
    case 3 : calc_p_Cs_and_gamma_tillotson(&(sph_system[i])); break;
    case 4 : calc_p_Cs_and_gamma_Mie_Gruneisen_EoS(&(sph_system[i])); break;
    }
    if(sph_system[i].pres < 0.0) sph_system[i].pres *= (1.0 - sph_system[i].damage); //fracture model
  }
}
