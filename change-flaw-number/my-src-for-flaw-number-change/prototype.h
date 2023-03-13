#pragma once

//equation-of-state.cpp
void calc_p_and_Cs_ideal_gas(RealPtcl* par_i);
void calc_delp_delrho_and_delp_delu_ideal_gas(PS::F64 rho,PS::F64 u,PS::S32 p_tag,PS::F64* delp_delrho_P,PS::F64* delp_delu_P);
void calc_p_elastic_eos(RealPtcl* par_i);
void calc_delp_delrho_and_delp_delu_elastic_eos(PS::F64 rho,PS::F64 u,PS::S32 p_tag,PS::F64* delp_delrho_P,PS::F64* delp_delu_P);
void calc_p_Cs_and_gamma_stiffened_gas_EoS(RealPtcl* par_i);
void calc_delp_delrho_and_delp_delu_stiffened_gas_EoS(PS::F64 rho,PS::F64 u,PS::S32 p_tag,PS::F64* delp_delrho_P,PS::F64* delp_delu_P);
void calc_delp_delrho_and_delp_delu_Mie_Gruneisen_EoS(PS::F64 rho,PS::F64 u,PS::S32 p_tag,PS::F64* delp_delrho_P,PS::F64* delp_delu_P);
void calc_p_Cs_and_gamma_Mie_Gruneisen_EoS(RealPtcl* par_i);
PS::F64 calc_ps_tillotson(PS::F64 E,PS::F64 rho,PS::S32 p_tag);
PS::F64 calc_pg_tillotson(PS::F64 E,PS::F64 rho,PS::S32 p_tag);
PS::F64 calc_dps_drho_tillotson(PS::F64 E,PS::F64 rho,PS::F64 dE_drho,PS::S32 p_tag);
PS::F64 calc_dpg_drho_tillotson(PS::F64 E,PS::F64 rho,PS::F64 dE_drho,PS::S32 p_tag);
void calc_tillotson_gamma_and_Cs(PS::F64 E,PS::F64 rho,PS::F64 p,PS::S32 p_tag,PS::F64* gamma_P,PS::F64* Cs_P);
void calc_p_Cs_and_gamma_tillotson(RealPtcl* par_i);
void calc_delps_delrho_and_delps_delu(PS::F64 E,PS::F64 rho,PS::F64* delp_delrho_P,PS::F64* delp_delu_P,PS::S32 p_tag);
void calc_delpg_delrho_and_delpg_delu(PS::F64 E,PS::F64 rho,PS::F64* delp_delrho_P,PS::F64* delp_delu_P,PS::S32 p_tag);
void calc_delp_delrho_and_delp_delu_tillotson(PS::F64 E,PS::F64 rho,PS::S32 p_tag,PS::F64* delp_delrho_P,PS::F64* delp_delu_P);
void calc_pressures(PS::ParticleSystem<RealPtcl>& sph_system);

//init.cpp
void read_particle_information(PS::ParticleSystem<RealPtcl>& sph_system, char particlebin_file[]);
void set_domain(PS::ParticleSystem<RealPtcl>& sph_system, PS::DomainInfo& dinfo);
void output_parameters(char particlebin_file[]);
void My_SPH_initialize(PS::ParticleSystem<RealPtcl>& sph_system, PS::DomainInfo& dinfo, PS::TreeForForceShort<RESULT::Dens, EPI::Dens, EPJ::Dens>::Gather_or_Scatter& dens_tree, PS::TreeForForceShort<RESULT::Gradients, EPI::Gradients, EPJ::Gradients>::Gather_or_Scatter& gradients_tree, PS::TreeForForceShort<RESULT::Elastic, EPI::Elastic, EPJ::Elastic>::Symmetry_or_Scatter& hydr_tree, PS::TreeForForceLong <RESULT::Grav , EPI::Grav , EPJ::Grav >::Monopole& grav_tree, char particlebin_file[], int argc, char* argv[]);

//integral.cpp
PS::F64 getTimeStep(void);
void setTimeStepGlobal(const PS::ParticleSystem<RealPtcl>& sph_system);
void setTime(PS::F64 t);
PS::F64 getTime(void);
void TimeDevelopment_for_Euler(PS::ParticleSystem<RealPtcl>& sph_system);
void TimeDevelopment_for_rungekutta(PS::ParticleSystem<RealPtcl>& sph_system);
void initial_kick_for_leapfrog(PS::ParticleSystem<RealPtcl>& sph_system);
void TimeDevelopment_for_leapfrog(PS::ParticleSystem<RealPtcl>& sph_system);
void time_development(PS::ParticleSystem<RealPtcl>& sph_system);
void check_not_a_number(PS::ParticleSystem<RealPtcl>& sph_system);
void one_time_development_roop(PS::ParticleSystem<RealPtcl>& sph_system, PS::DomainInfo& dinfo, PS::TreeForForceShort<RESULT::Dens, EPI::Dens, EPJ::Dens>::Gather_or_Scatter& dens_tree, PS::TreeForForceShort<RESULT::Gradients, EPI::Gradients, EPJ::Gradients>::Gather_or_Scatter& gradients_tree, PS::TreeForForceShort<RESULT::Elastic, EPI::Elastic, EPJ::Elastic>::Symmetry_or_Scatter& hydr_tree, PS::TreeForForceLong <RESULT::Grav , EPI::Grav , EPJ::Grav >::Monopole& grav_tree);

//Tensor_Product.cpp
PS::matrix CalcTensorProduct(const PS::F64vec arg1, const PS::F64vec arg2);
PS::F64 CalcTensorToScalar(const PS::matrix arg1, const PS::matrix arg2);

//fracture-and-porosity-model.cpp
void calc_effective_strain(PS::ParticleSystem<RealPtcl>& sph_system);
void time_development_of_D_and_alpha(PS::ParticleSystem<RealPtcl>& sph_system);
void calc_dalpha_dt_and_modify_dSab_dt_for_porosity(PS::ParticleSystem<RealPtcl>& sph_system);
