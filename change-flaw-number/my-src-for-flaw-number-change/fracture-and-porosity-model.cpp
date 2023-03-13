#include "header.h"

//calculate effective strain for all particles. 
void calc_effective_strain(PS::ParticleSystem<RealPtcl>& sph_system)
{
  PS::S64 i,d,d2;
  PS::S64 particle_number = sph_system.getNumberOfParticleLocal();

  PS::matrix sigma_eff_i;
  PS::matrix deltaab = PS::matrix(1, 0, 0,
	  0, 1, 0,
	  0, 0, 1 );//Kronecker delta
  PS::F64 a,b,c,q,p;
  double _Complex alpha_plus,alpha_minus;
  PS::F64 sigma_1,sigma_2,sigma_3;
  PS::F64 sigma_max;

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for private(i,d,d2,sigma_eff_i,a,b,c,p,q,alpha_plus,alpha_minus,sigma_1,sigma_2,sigma_3,sigma_max)
#endif
  for(i=0;i<particle_number;i++){
    sigma_eff_i = - deltaab * sph_system[i].pres + sph_system[i].Sab;
    
    a = - (sigma_eff_i.xx + sigma_eff_i.yy + sigma_eff_i.zz);
    b = sigma_eff_i.xx*sigma_eff_i.yy + sigma_eff_i.yy*sigma_eff_i.zz + sigma_eff_i.zz*sigma_eff_i.xx - sigma_eff_i.xy*sigma_eff_i.xy - sigma_eff_i.yz*sigma_eff_i.yz - sigma_eff_i.zx*sigma_eff_i.zx;
    c = - (sigma_eff_i.xx*sigma_eff_i.yy*sigma_eff_i.zz + 2.0*sigma_eff_i.xy*sigma_eff_i.yz*sigma_eff_i.zx - sigma_eff_i.xx*sigma_eff_i.yz*sigma_eff_i.yz - sigma_eff_i.yy*sigma_eff_i.zx*sigma_eff_i.zx - sigma_eff_i.zz*sigma_eff_i.xy*sigma_eff_i.xy);
    q = ( 27.0*c + 2.0*a*a*a - 9.0*a*b ) / 54.0;
    p = ( 3.0*b - a*a ) / 9.0;
    alpha_plus = cpow( -q + csqrt( (double _Complex)(q*q + p*p*p) ) , 1.0/3.0 );
    alpha_minus = cpow( -q - csqrt( (double _Complex)(q*q + p*p*p) ) , 1.0/3.0 );
    sigma_1 = creal( alpha_plus + alpha_minus - (a/3.0) );
    sigma_2 = creal( 0.5 * ( -1.0 + I*sqrt(3.0) ) * alpha_plus + 0.5 * ( -1.0 - I*sqrt(3.0) ) * alpha_minus - (a/3.0) );
    sigma_3 = creal( 0.5 * ( -1.0 - I*sqrt(3.0) ) * alpha_plus + 0.5 * ( -1.0 + I*sqrt(3.0) ) * alpha_minus - (a/3.0) );
    sigma_max = ( sigma_1 > sigma_2 ? sigma_1 : sigma_2 );
    sigma_max = ( sigma_max > sigma_3 ? sigma_max : sigma_3 );

    sph_system[i].eps_i_eff = sigma_max / ( PARAM::K_BULK[sph_system[i].property_tag] + (4.0/3.0) * PARAM::MU_SHEAR[sph_system[i].property_tag] );
  }
  
}

//time development of damage parameter
void time_development_of_D_and_alpha(PS::ParticleSystem<RealPtcl>& sph_system)
{
  PS::F64 D_i_max,Cg,Rs;
  PS::S64 ni;
  PS::S64 i,j;
  PS::F64 delta_D=0.01;
  PS::S64 particle_number=sph_system.getNumberOfParticleLocal();
  PS::F64 TimeStep=getTimeStep();

  calc_effective_strain(sph_system);
  
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for private(i,j,ni,D_i_max,Cg,Rs)
#endif
  for(i=0;i<particle_number;i++){
    if(PARAM::IS_FRACTURE_MODEL[sph_system[i].property_tag]){
      //calculate the number of activated flaws for this particle
      ni=0;
      for(j=0 ; j<sph_system[i].ni_tot_flaw ; j++){
        if( sph_system[i].eps_ij_act_flaw[j] < sph_system[i].eps_i_eff ) ni++;
      }
      //calculate Dmax of this particle
      if(sph_system[i].ni_tot_flaw == 0){
        D_i_max = 0.0;
      }
      else{
        D_i_max = pow( ni / sph_system[i].ni_tot_flaw, 1.0/3.0 );
      }
      
      Cg = 0.4 * sqrt( sph_system[i].snds*sph_system[i].snds + 4.0*PARAM::MU_SHEAR[sph_system[i].property_tag]/(3.0*sph_system[i].dens) );
      Rs = 3.0 * sph_system[i].smth;
      
      if(sph_system[i].damage >= D_i_max) sph_system[i].dD_cbrt_dt = 0.0;
      else sph_system[i].dD_cbrt_dt = ( Cg / Rs ) *ni;
      //add damage increasing due to pore compaction
      if(PARAM::IS_POROSITY_MODEL[sph_system[i].property_tag]){
        sph_system[i].dD_cbrt_dt += - (1.0/3.0) * ( pow( 1.0 - (sph_system[i].alpha_por-1.0)/(PARAM::ALPHA_POR_0[sph_system[i].property_tag]-1.0) + delta_D , -(2.0/3.0) ) / ( pow(1.0+delta_D,1.0/3.0) - pow(delta_D,1.0/3.0) ) ) * ( 1.0 / (PARAM::ALPHA_POR_0[sph_system[i].property_tag]-1.0) ) * sph_system[i].dalpha_dt;
      }
      
      sph_system[i].D_cbrt += sph_system[i].dD_cbrt_dt * TimeStep;
      
      sph_system[i].damage = pow(sph_system[i].D_cbrt,3.0);
      if(sph_system[i].damage >= 1.0){
        sph_system[i].damage = 1.0;
        sph_system[i].D_cbrt = 1.0;
      }
    }
    else{
      sph_system[i].damage = sph_system[i].D_cbrt = 0.0;
    }
    if(PARAM::IS_POROSITY_MODEL[sph_system[i].property_tag]){
      //time development of distension parameter
      sph_system[i].alpha_por += sph_system[i].dalpha_dt * TimeStep; 
      if(sph_system[i].alpha_por < 1.0) sph_system[i].alpha_por=1.0;
    }
  }
}

//calculate time derivative of distension parameter and modify time derivative of deviatoric stress tensor for porosity model.
void calc_dalpha_dt_and_modify_dSab_dt_for_porosity(PS::ParticleSystem<RealPtcl>& sph_system)
{
  PS::F64 dalpha_dp;
  PS::F64 delp_delrho,delp_delu;
  PS::F64 f,dalpha_drho;
  PS::S64 i;
  PS::S64 particle_number=sph_system.getNumberOfParticleLocal();
  PS::S64 count=0;
  PS::S32 error=0;

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for private(i,dalpha_dp,delp_delrho,delp_delu,f,dalpha_drho) shared(error) reduction(+:count)
#endif
  for(i=0;i<particle_number;i++){
    if(PARAM::IS_POROSITY_MODEL[sph_system[i].property_tag]){
      //if pressure is smaller than elastic limit or unloading, time derivative of distension should be 0.
      if( sph_system[i].pres < PARAM::PE_POR[sph_system[i].property_tag] || (sph_system[i].pres - sph_system[i].pres_old) < 0.0 ){
        dalpha_dp = 0.0;
        sph_system[i].dalpha_dt = 0.0;
        f = 1.0;
      }
      else if( sph_system[i].pres > PARAM::PS_POR[sph_system[i].property_tag] ){
        dalpha_dp = 0.0;
        sph_system[i].dalpha_dt = 0.0;
        f = 1.0;
        sph_system[i].alpha_por = 1.0;
      }
      else{
        dalpha_dp = -2.0 * ( PARAM::ALPHA_POR_0[sph_system[i].property_tag] - 1.0 ) * ( PARAM::PS_POR[sph_system[i].property_tag] - sph_system[i].pres ) / pow(PARAM::PS_POR[sph_system[i].property_tag]-PARAM::PE_POR[sph_system[i].property_tag],2.0);
        //calculate delp_delrho and delp_delu
        switch(PARAM::EQUATION_OF_STATE[sph_system[i].property_tag]){
        case 0 : calc_delp_delrho_and_delp_delu_ideal_gas(sph_system[i].alpha_por*sph_system[i].dens,sph_system[i].eng,sph_system[i].property_tag,&delp_delrho,&delp_delu); break;
        case 1 : calc_delp_delrho_and_delp_delu_elastic_eos(sph_system[i].alpha_por*sph_system[i].dens,sph_system[i].eng,sph_system[i].property_tag,&delp_delrho,&delp_delu); break;
        case 2 : calc_delp_delrho_and_delp_delu_stiffened_gas_EoS(sph_system[i].alpha_por*sph_system[i].dens,sph_system[i].eng,sph_system[i].property_tag,&delp_delrho,&delp_delu); break;
        case 3 : calc_delp_delrho_and_delp_delu_tillotson(sph_system[i].eng,sph_system[i].alpha_por*sph_system[i].dens,sph_system[i].property_tag,&delp_delrho,&delp_delu); break;
        case 4 : calc_delp_delrho_and_delp_delu_Mie_Gruneisen_EoS(sph_system[i].alpha_por*sph_system[i].dens,sph_system[i].eng,sph_system[i].property_tag,&delp_delrho,&delp_delu); break;
        }
        
        sph_system[i].dalpha_dt = ( ( sph_system[i].eng_dot*delp_delu + sph_system[i].alpha_por*sph_system[i].drho_dt*delp_delrho )/( sph_system[i].alpha_por + dalpha_dp*(sph_system[i].pres - sph_system[i].dens*delp_delrho) ) )*dalpha_dp;
        if(sph_system[i].dalpha_dt>0.0){
          sph_system[i].dalpha_dt=0.0;
          count++;
        }
        dalpha_drho = ( ( (sph_system[i].pres/(sph_system[i].dens*sph_system[i].dens))*delp_delu + sph_system[i].alpha_por*delp_delrho )/( sph_system[i].alpha_por + dalpha_dp*(sph_system[i].pres - sph_system[i].dens*delp_delrho) ) )*dalpha_dp;
        f = 1.0 + dalpha_drho * ( sph_system[i].dens / sph_system[i].alpha_por );
      }
      
      //modify time derivative of deviatoric stress tensor for porosity model
      sph_system[i].dSab_rho_dt = ( f / sph_system[i].alpha_por ) * sph_system[i].dSab_rho_dt - ( 1.0 / sph_system[i].alpha_por ) * ( sph_system[i].Sab / sph_system[i].dens ) * sph_system[i].dalpha_dt;
      if( std::isnan(f) ){
        printf("par %d's f is not a number! \n",sph_system[i].id);
        error=1;
      }
    }
  }
  
  if(PS::Comm::getSum(error)>=1){
    if(PS::Comm::getRank()==0){
      printf("there are not a numbers\n");
    }
    PS::Finalize();
    exit(1);
  }

}
