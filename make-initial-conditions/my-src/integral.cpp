#include "header.h"

static PS::F64 TimeStep;
static PS::F64 Time;
static PS::S32 time_step_sign=1;
static PS::S32 leapfrog_init=0;

PS::F64 getTimeStep(void)
{
  return(TimeStep);
}

void setTimeStepGlobal(const PS::ParticleSystem<RealPtcl>& sph_system)
{
  PS::F64 dt = 1.0e+30;//set VERY LARGE VALUE
  PS::F64 dt_temp;
  PS::S64 particle_number_local = sph_system.getNumberOfParticleLocal();
  PS::S64 i;
  PS::F64 Cspi;

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for private(i, dt_temp, Cspi) shared(particle_number_local) reduction(min:dt)
#endif
  for(i = 0 ; i < particle_number_local ; ++ i){
    //in the case of fluid
    if(PARAM::PLASTIC_MODEL[sph_system[i].property_tag]==4||PARAM::CONSIDER_PERP_WAVE_FOR_DT==0){
      Cspi = sph_system[i].snds;
    }
    //for elastic case
    else{
      Cspi = sqrt( sph_system[i].snds*sph_system[i].snds + (4.0*PARAM::MU_SHEAR[sph_system[i].property_tag])/(3.0*sph_system[i].dens) );
    }
    dt_temp = PARAM::C_CFL * ( sph_system[i].smth / ( Cspi + sph_system[i].smth * fabs(sph_system[i].gradvel.getTrace()) ) );
    dt = std::min(dt, dt_temp);
    if(PARAM::IS_SELF_GRAVITY){
      dt_temp = PARAM::C_CFL * sqrt( sph_system[i].smth / sqrt(sph_system[i].grav.x*sph_system[i].grav.x + sph_system[i].grav.y*sph_system[i].grav.y + sph_system[i].grav.z*sph_system[i].grav.z) );
      dt = std::min(dt, dt_temp);
    }
  }
  TimeStep = PS::Comm::getMinValue(dt);
}

void setTime(PS::F64 t)
{
  Time = t;
}

PS::F64 getTime(void)
{
  return(Time);
}

//calculate time development of all particles for Euler method
void TimeDevelopment_for_Euler(PS::ParticleSystem<RealPtcl>& sph_system)
{
  PS::S64 i;
  PS::S64 particle_number_local = sph_system.getNumberOfParticleLocal();
  PS::F64vec acc_tot;

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for private(i, acc_tot) shared(particle_number_local)
#endif
  for(i = 0 ; i < particle_number_local ; ++ i){ 
    acc_tot = sph_system[i].acc + sph_system[i].grav;
    
    sph_system[i].pos += ( sph_system[i].vel + 0.5 * acc_tot * TimeStep ) * TimeStep;
    sph_system[i].vel += acc_tot * TimeStep;
    sph_system[i].eng += sph_system[i].eng_dot * TimeStep;
    sph_system[i].Sab_rho += sph_system[i].dSab_rho_dt * TimeStep;
    if(PARAM::DENSITY_DEVELOPMENT_METHOD==1||PARAM::DENSITY_DEVELOPMENT_METHOD==2){
      sph_system[i].dens += sph_system[i].drho_dt * TimeStep;
      if(sph_system[i].dens < PARAM::MINIMUM_DENSITY) sph_system[i].dens = PARAM::MINIMUM_DENSITY;
      if(PARAM::IS_VARIABLE_H==1){
        sph_system[i].smth = PARAM::SMTH * pow( sph_system[i].mass / sph_system[i].dens , 1.0/PARAM::Dim );
      }
    }
    
    //copy velocity to old array regarding energy conservation/////////////
    sph_system[i].vel_old = sph_system[i].vel;
    //////////////////////////////////////////////////////////////////////
    
    if(sph_system[i].eng < PARAM::ACC) sph_system[i].eng = PARAM::ACC;
    
    //save old pressure value for porosity model/////////////////////////
    sph_system[i].pres_old = sph_system[i].pres;
    /////////////////////////////////////////////////////////////////////
  }
  Time += TimeStep;
}

//calculate time development of all particles for second order rungekutta method
void TimeDevelopment_for_rungekutta(PS::ParticleSystem<RealPtcl>& sph_system)
{
  PS::S64 i;
  PS::S64 particle_number_local = sph_system.getNumberOfParticleLocal();
  PS::F64vec acc_tot;

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for private(i, acc_tot) shared(particle_number_local)
#endif
  for(i = 0 ; i < particle_number_local ; ++ i){
    acc_tot = sph_system[i].acc + sph_system[i].grav;

    if( time_step_sign == 1 ){
      sph_system[i].pos_old = sph_system[i].pos;
      sph_system[i].vel_old = sph_system[i].vel;
      sph_system[i].eng_old = sph_system[i].eng;
      sph_system[i].Sab_rho_old = sph_system[i].Sab_rho;
      
      sph_system[i].pos += ( sph_system[i].vel + 0.5 * acc_tot * 0.5 * TimeStep ) * 0.5 * TimeStep;
      sph_system[i].vel += acc_tot * 0.5 * TimeStep;
      sph_system[i].eng += sph_system[i].eng_dot * 0.5 * TimeStep;
      sph_system[i].Sab_rho += sph_system[i].dSab_rho_dt * 0.5 * TimeStep;

      if(sph_system[i].eng < PARAM::ACC) sph_system[i].eng = PARAM::ACC;

      if(PARAM::DENSITY_DEVELOPMENT_METHOD==1||PARAM::DENSITY_DEVELOPMENT_METHOD==2){
        sph_system[i].dens_old = sph_system[i].dens;
        sph_system[i].dens += sph_system[i].drho_dt * 0.5 * TimeStep;
        if(sph_system[i].dens < PARAM::MINIMUM_DENSITY) sph_system[i].dens = PARAM::MINIMUM_DENSITY;
        if(PARAM::IS_VARIABLE_H==1){
          sph_system[i].smth = PARAM::SMTH * pow( sph_system[i].mass / sph_system[i].dens , 1.0/PARAM::Dim );
        }
      }
      
      //save old pressure value for porosity model/////////////////////////
      sph_system[i].pres_old = sph_system[i].pres;
      /////////////////////////////////////////////////////////////////////
    }
    else{
      sph_system[i].pos = sph_system[i].pos_old + ( sph_system[i].vel_old + 0.5 * TimeStep * acc_tot ) * TimeStep;
      sph_system[i].vel = sph_system[i].vel_old + acc_tot * TimeStep;
      sph_system[i].eng = sph_system[i].eng_old + sph_system[i].eng_dot * TimeStep;
      sph_system[i].Sab_rho = sph_system[i].Sab_rho_old + sph_system[i].dSab_rho_dt * TimeStep;

      if(PARAM::DENSITY_DEVELOPMENT_METHOD==1||PARAM::DENSITY_DEVELOPMENT_METHOD==2){
        sph_system[i].dens = sph_system[i].dens_old + sph_system[i].drho_dt * TimeStep;
        if(sph_system[i].dens < PARAM::MINIMUM_DENSITY) sph_system[i].dens = PARAM::MINIMUM_DENSITY;
        if(PARAM::IS_VARIABLE_H==1){
          sph_system[i].smth = PARAM::SMTH * pow( sph_system[i].mass / sph_system[i].dens , 1.0/PARAM::Dim );
        }
      }
      
      //copy velocity to old array regarding energy conservation/////////////
      sph_system[i].vel_old = sph_system[i].vel;
      //////////////////////////////////////////////////////////////////////
    }
    if(sph_system[i].eng < PARAM::ACC) sph_system[i].eng = PARAM::ACC;
  }
  if(time_step_sign==-1) Time += TimeStep;
  time_step_sign *= -1;
}

//initial kick for leapfrog method
void initial_kick_for_leapfrog(PS::ParticleSystem<RealPtcl>& sph_system)
{
  PS::S64 i;
  PS::S64 particle_number_local = sph_system.getNumberOfParticleLocal();
  PS::F64vec acc_tot;

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for private(i, acc_tot) shared(particle_number_local)
#endif
  for(i = 0 ; i < particle_number_local ; ++ i){ 
    acc_tot = sph_system[i].acc + sph_system[i].grav;

    sph_system[i].pos_old = sph_system[i].pos;
    sph_system[i].vel_old = sph_system[i].vel;
    sph_system[i].eng_old = sph_system[i].eng;
    sph_system[i].Sab_rho_old = sph_system[i].Sab_rho;

    sph_system[i].pos += ( sph_system[i].vel + 0.5 * acc_tot * TimeStep ) * TimeStep;
    sph_system[i].vel += acc_tot * TimeStep;
    sph_system[i].eng += sph_system[i].eng_dot * TimeStep;
    sph_system[i].Sab_rho += sph_system[i].dSab_rho_dt * TimeStep;

    sph_system[i].acc_old = sph_system[i].acc;
    sph_system[i].eng_dot_old = sph_system[i].eng_dot;
    sph_system[i].grav_old = sph_system[i].grav;
    sph_system[i].dSab_rho_dt_old = sph_system[i].dSab_rho_dt;

    if(PARAM::DENSITY_DEVELOPMENT_METHOD==1||PARAM::DENSITY_DEVELOPMENT_METHOD==2){
      sph_system[i].dens_old = sph_system[i].dens;
      sph_system[i].dens += sph_system[i].drho_dt * TimeStep;
      if(sph_system[i].dens < PARAM::MINIMUM_DENSITY) sph_system[i].dens = PARAM::MINIMUM_DENSITY;
      if(PARAM::IS_VARIABLE_H==1){
        sph_system[i].smth = PARAM::SMTH * pow( sph_system[i].mass / sph_system[i].dens , 1.0/PARAM::Dim );
      }
      sph_system[i].drho_dt_old = sph_system[i].drho_dt;
    }
    
    if(sph_system[i].eng < PARAM::ACC) sph_system[i].eng = PARAM::ACC;

    //save old pressure value for porosity model/////////////////////////
    sph_system[i].pres_old = sph_system[i].pres;
    /////////////////////////////////////////////////////////////////////
  }
}

//calculate time development of all particles for leapfrog method
void TimeDevelopment_for_leapfrog(PS::ParticleSystem<RealPtcl>& sph_system)
{
  PS::S64 i;
  PS::S64 particle_number_local = sph_system.getNumberOfParticleLocal();
  PS::F64vec acc_tot, acc_tot_old;

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for private(i, acc_tot, acc_tot_old) shared(particle_number_local)
#endif
  for(i = 0 ; i < particle_number_local ; ++ i){
    acc_tot = sph_system[i].acc + sph_system[i].grav;
    acc_tot_old = sph_system[i].acc_old + sph_system[i].grav_old;
 
    sph_system[i].vel_old += 0.5 * ( acc_tot_old + acc_tot ) * TimeStep;
    sph_system[i].eng_old += 0.5 * ( sph_system[i].eng_dot_old + sph_system[i].eng_dot ) * TimeStep;
    sph_system[i].Sab_rho_old += 0.5 * ( sph_system[i].dSab_rho_dt_old + sph_system[i].dSab_rho_dt ) * TimeStep;

    sph_system[i].pos += ( sph_system[i].vel_old + 0.5 * acc_tot * TimeStep ) * TimeStep;
    sph_system[i].vel = sph_system[i].vel_old + acc_tot * TimeStep;
    sph_system[i].eng = sph_system[i].eng_old + sph_system[i].eng_dot * TimeStep;
    sph_system[i].Sab_rho = sph_system[i].Sab_rho_old + sph_system[i].dSab_rho_dt * TimeStep;

    sph_system[i].acc_old = sph_system[i].acc;
    sph_system[i].eng_dot_old = sph_system[i].eng_dot;
    sph_system[i].grav_old = sph_system[i].grav;
    sph_system[i].dSab_rho_dt_old = sph_system[i].dSab_rho_dt;

    if(PARAM::DENSITY_DEVELOPMENT_METHOD==1||PARAM::DENSITY_DEVELOPMENT_METHOD==2){
      sph_system[i].dens_old += 0.5 * ( sph_system[i].drho_dt_old + sph_system[i].drho_dt ) * TimeStep;
      if(sph_system[i].dens_old < PARAM::MINIMUM_DENSITY) sph_system[i].dens_old = PARAM::MINIMUM_DENSITY;
      sph_system[i].dens = sph_system[i].dens_old + sph_system[i].drho_dt * TimeStep;
      if(sph_system[i].dens < PARAM::MINIMUM_DENSITY) sph_system[i].dens = PARAM::MINIMUM_DENSITY;
      if(PARAM::IS_VARIABLE_H==1){
        sph_system[i].smth = PARAM::SMTH * pow( sph_system[i].mass / sph_system[i].dens , 1.0/PARAM::Dim );
      }
      sph_system[i].drho_dt_old = sph_system[i].drho_dt;
    }
    
    if(sph_system[i].eng < PARAM::ACC) sph_system[i].eng = PARAM::ACC;

    //save old pressure value for porosity model/////////////////////////
    sph_system[i].pres_old = sph_system[i].pres;
    /////////////////////////////////////////////////////////////////////
  }
  Time += TimeStep;
}

//calculate time development. This function only select time development method.
void time_development(PS::ParticleSystem<RealPtcl>& sph_system)
{
  switch(PARAM::TIME_DEVELOPMENT_METHOD){
  case 0 : 
    TimeDevelopment_for_Euler(sph_system);
    break;
  case 1 :
    TimeDevelopment_for_rungekutta(sph_system);
    break;
  case 2 :
    if(leapfrog_init==0){
      initial_kick_for_leapfrog(sph_system);
      leapfrog_init=1;
    }
    else{
      TimeDevelopment_for_leapfrog(sph_system);
    }
    break;
  }
}

//check not a number
void check_not_a_number(PS::ParticleSystem<RealPtcl>& sph_system)
{
  PS::S64 i;
  PS::S64 particle_number_local=sph_system.getNumberOfParticleLocal();
  PS::S32 error=0;

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for private(i) shared(particle_number_local,error)
#endif
  for(i=0 ; i<particle_number_local ; i++){
    if( std::isnan(sph_system[i].acc.x) || std::isnan(sph_system[i].acc.y) || std::isnan(sph_system[i].acc.z) ){
      //printf("%d-th particle's acceleration is not a number! \n",sph_system[i].id);
      error=1;
    }
    if( std::isnan(sph_system[i].eng_dot) ){
      //printf("%d-th particle's dudt is not a number! \n",sph_system[i].id);
      error=1;
    }
  }

  if(PS::Comm::getSum(error)>=1){
    if(PS::Comm::getRank()==0){
      printf("There are not a numbers! \n");
    }
    PS::Finalize();
    exit(1);
  }
}

//calculate Sab from Sab/density and density
//In this fanction, plastic model and friction model are also applied
void calc_Sab(PS::ParticleSystem<RealPtcl>& sph_system)
{
  PS::S64 i;
  PS::S64 particle_number_local = sph_system.getNumberOfParticleLocal();
  PS::F64 J2; //second invariant of deviatoric stress tensor
  PS::F64 f;
  PS::F64 Y_yield,Yi,Yd;
  PS::F64 P_pos,u_yield;
  PS::F64 Y0_mises,Y0_cohesion,mu_i_fric,mu_d_fric,u_melt;

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for private(i,J2,f,Y_yield,Yi,Yd,P_pos,u_yield,Y0_mises,Y0_cohesion,mu_i_fric,mu_d_fric,u_melt) shared(particle_number_local)
#endif
  for(i=0 ; i<particle_number_local ; i++){
    sph_system[i].Sab = sph_system[i].dens * sph_system[i].Sab_rho;
    J2 = sph_system[i].Sab.getJ2();
    Y0_mises=PARAM::Y0_MISES[sph_system[i].property_tag];
    Y0_cohesion=PARAM::Y0_COHESION[sph_system[i].property_tag];
    mu_i_fric=PARAM::MU_I_FRIC[sph_system[i].property_tag];
    mu_d_fric=PARAM::MU_D_FRIC[sph_system[i].property_tag];
    u_melt=PARAM::U_MELT[sph_system[i].property_tag];

    switch(PARAM::PLASTIC_MODEL[sph_system[i].property_tag]){
    case 0 :
      //elastic only
      f=1;
      if(PARAM::IS_FRACTURE_MODEL[sph_system[i].property_tag]) f *= ( 1.0 - sph_system[i].damage );
      break;
    case 1 :
      //the way to describe elastic-perfectly plastic material in Benz and Asphaug 1995
      f = ( (Y0_mises*Y0_mises/(3.0*J2))<1.0 ) ? (Y0_mises*Y0_mises/(3.0*J2)) : 1.0 ;
      if(PARAM::IS_FRACTURE_MODEL[sph_system[i].property_tag]) f *= ( 1.0 - sph_system[i].damage );
      break;
    case 2 :
      //the way to describe elastic-perfectly plastic material in Libersky and Petschek 1990
      f = ( (sqrt(Y0_mises*Y0_mises/(3.0*J2)))<1.0 ) ? (sqrt(Y0_mises*Y0_mises/(3.0*J2))) : 1.0 ;
      if(PARAM::IS_FRACTURE_MODEL[sph_system[i].property_tag]) f *= ( 1.0 - sph_system[i].damage );
      break;
    case 3 :
      //pressure dependent failure model with friction in Jutzi 2015
      P_pos = ( ( sph_system[i].pres > 0 ) ? sph_system[i].pres : 0 );
      u_yield = ( ( sph_system[i].eng < u_melt ) ? sph_system[i].eng : u_melt ); 
      Yi = ( Y0_cohesion + ( mu_i_fric * P_pos ) / ( 1.0 + mu_i_fric * P_pos / ( Y0_mises - Y0_cohesion ) ) ) * ( 1.0 - u_yield / u_melt );
      Yd = mu_d_fric * P_pos;
      if( Yd > Yi ) Yd = Yi;
      if(!PARAM::IS_FRICTION_MODEL[sph_system[i].property_tag]) Yd = 0.0;
      Y_yield = ( 1.0 - sph_system[i].damage ) * Yi + sph_system[i].damage * Yd;
      if( Y_yield > Yi || !PARAM::IS_FRACTURE_MODEL[sph_system[i].property_tag] ) Y_yield = Yi;
      f = ( Y_yield/sqrt(J2) < 1.0 ) ? (Y_yield/sqrt(J2)) : 1.0 ;
      break;
    case 4 :
      //if we want to consider fluid case
      f=0.0;
      break;
    }
    if( std::isnan(f) ){
      f=0.0;
    }
    sph_system[i].f_test = f;
    sph_system[i].Sab_rho *= f;
    if(PARAM::TIME_DEVELOPMENT_METHOD==2){
      sph_system[i].Sab_rho_old *= f;
    }
    sph_system[i].Sab = sph_system[i].Sab_rho * sph_system[i].dens;
  }
}

//remember every particle's max internal energy
void set_eng_max(PS::ParticleSystem<RealPtcl>& sph_system)
{
  PS::S64 i;
  PS::S64 particle_number_local = sph_system.getNumberOfParticleLocal();

#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for private(i) shared(particle_number_local)
#endif
  for(i=0 ; i<particle_number_local ; i++){
	if(sph_system[i].eng > sph_system[i].eng_max){
	  sph_system[i].eng_max = sph_system[i].eng;
	}
  }
}

//one time development roop
void one_time_development_roop(PS::ParticleSystem<RealPtcl>& sph_system, PS::DomainInfo& dinfo, PS::TreeForForceShort<RESULT::Dens, EPI::Dens, EPJ::Dens>::Gather_or_Scatter& dens_tree, PS::TreeForForceShort<RESULT::Gradients, EPI::Gradients, EPJ::Gradients>::Gather_or_Scatter& gradients_tree, PS::TreeForForceShort<RESULT::Elastic, EPI::Elastic, EPJ::Elastic>::Symmetry_or_Scatter& elast_tree, PS::TreeForForceLong <RESULT::Grav , EPI::Grav , EPJ::Grav >::Monopole& grav_tree)
{
  PS::S64 i;
  PS::S64 particle_number_local=sph_system.getNumberOfParticleLocal();

  //if we use periodic boundary condition, particles that are out of region should be re-located to inside of region regarding periodic boundary.
  if(PARAM::IS_PERIODIC_BOUNDARY.x==1||PARAM::IS_PERIODIC_BOUNDARY.y==1||PARAM::IS_PERIODIC_BOUNDARY.z==1){
    sph_system.adjustPositionIntoRootDomain(dinfo);
  }
  
  dinfo.decomposeDomainAll(sph_system);//decompose domain
  
  sph_system.exchangeParticle(dinfo);//exchange particle
  
  if(PARAM::DENSITY_DEVELOPMENT_METHOD==0){
    dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo); //calculate density
  }
  calc_pressures(sph_system);//calculate pressure, sound speed and gamma
  calc_Sab(sph_system);//calculate Sab from density and Sab/density
  gradients_tree.calcForceAllAndWriteBack(CalcGradients(), sph_system, dinfo);//calculate gradients
  
  setTimeStepGlobal(sph_system);//set time step
  if(PARAM::IS_GODUNOV_SPH){
    elast_tree.calcForceAllAndWriteBack(CalcElasticForceGodunov(), sph_system, dinfo);//calculate force and dudt for Godunov SPH method
  }
  else{
    elast_tree.calcForceAllAndWriteBack(CalcElasticForceStand(), sph_system, dinfo);//calculate force and dudt for Standard SPH method
  }
  if(PARAM::IS_SELF_GRAVITY){
    grav_tree.calcForceAllAndWriteBack(CalcGravityForceCloseParticle(), CalcGravityForceSuperParticle(), sph_system, dinfo);
  }
  else{
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL	    
#pragma omp parallel for private(i) shared(particle_number_local)
#endif
    for(i=0 ; i<particle_number_local ; i++){
      sph_system[i].grav = 0.0;
    }
  }
  check_not_a_number(sph_system);//check there is no not a number made in calc force function
  
  calc_dalpha_dt_and_modify_dSab_dt_for_porosity(sph_system);//calculate time derivative of distension parameter and modify time derivative of deviatoric stress tensor for porosity model.
  
  time_development(sph_system);//time development

  set_eng_max(sph_system);//remember every particle's max internal energy
  
  if(PARAM::TIME_DEVELOPMENT_METHOD!=1||time_step_sign==-1){
    time_development_of_D_and_alpha(sph_system);
  }
  
}


