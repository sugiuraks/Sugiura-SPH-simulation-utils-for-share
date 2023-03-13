#include "header.h"

static PS::F64 **temp_ane_table;
static PS::F64 **dens_ane_table;
static PS::F64 ***eng_ane_table;
static PS::F64 ***pres_ane_table;
static PS::F64 ***snds_ane_table;

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
    case 5 : calc_p_Cs_and_temp_aneos(&(sph_system[i])); break;
    }
    if(sph_system[i].pres < 0.0) sph_system[i].pres *= (1.0 - sph_system[i].damage); //fracture model
  }
}
