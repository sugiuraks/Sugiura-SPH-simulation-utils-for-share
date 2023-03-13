#ifndef TILLOTSON
#define TILLOTSON

//fanctions which is need to be called in another fanction/////////////////////

//calculate tillotson pressure
//argument1...internal energy  argument2...density  argument3...this particle's property tag
double calc_tillotson_pressure(double E,double rho,int p_tag);

//calculate tillotson gamma
//argument1...internal energy  argument2...density  argument3...pressure  argument4...this particle's property tag
double calc_tillotson_gamma(double E,double rho,double p,int p_tag);

//calculate dp/drho for tillotson EoS (=square of sound speed)
//argument1...internal energy  argument2...density  argument3...pressure  argument4...this particle's property tag
double calc_dp_drho_tillotson(double E,double rho,double p,int p_tag);

//calculate E from p and rho using tillotson EoS
//argument1...density rho  argument2...pressure p  argument3...this particle's property tag
double calc_tillotson_E_from_p_and_rho(double rho,double p,int p_tag);

//calculate delp_delrho and delp_delu for tillotson EoS
//argument1...internal energy, argument2...density, argument3...propertytag of i-th particle, argument4...pointer to delp_delrho, argument5...pointer to delp_delu
void calc_delp_delrho_and_delp_delu_tillotson(double E,double rho,int p_tag,double* delp_delrho_P,double* delp_delu_P);

//////////////////////////////////////////////////////////////////////////////


//fanctions which is not need to be called in another fanction////////////////

//calculate tillotson pressure for compression or cold expansion state
//argument1...internal energy  argument2...density  argument3...this particle's property tag
double calc_ps(double E,double rho,int p_tag);

//calculate tillotson pressure for hot expansion state
//argument1...internal energy  argument2...density  argument3...this particle's property tag
double calc_pg(double E,double rho,int p_tag);

//calculate tillotson EOS's dp/drho for compression or cold expansion state
//argument1...internal energy  argument2...density  argument3...dE/drho  argument4...this particle's property tag
double calc_dps_drho(double E,double rho,double dE_drho,int p_tag);

//calculate tillotson EOS's dp/drho for hot expansion state
//argument1...internal energy  argument2...density  argument3...dE/drho  argument4...this particle's property tag
double calc_dpg_drho(double E,double rho,double dE_drho,int p_tag);

//calculate delp_delrho and delp_delu for compression or cold expansion state
//argument1...internal energy, argument2...density, argument3...pointer to delp_delrho, argument4...pointer to delp_delu, argument5...propertytag of i-th particle
void calc_delps_delrho_and_delps_delu(double E,double rho,double* delp_delrho_P,double* delp_delu_P,int p_tag);

//calculate delp_delrho and delp_delu for hot expansion state
//argument1...internal energy, argument2...density, argument3...pointer to delp_delrho, argument4...pointer to delp_delu, argument5...propertytag of i-th particle
void calc_delpg_delrho_and_delpg_delu(double E,double rho,double* delp_delrho_P,double* delp_delu_P,int p_tag);

/////////////////////////////////////////////////////////////////////////////
#endif
