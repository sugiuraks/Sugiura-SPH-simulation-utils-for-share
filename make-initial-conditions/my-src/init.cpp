#include "header.h"

void read_particle_information(PS::ParticleSystem<RealPtcl>& sph_system, char particlebin_file[])
{
  FileHeader header;
  PS::S64 i,particle_number_local;
  PS::S32 error=0;

  //read time, particle number and all information of particles
  sph_system.readParticleBinary(particlebin_file, header);
  //scatter initial time
  PS::Comm::broadcast(&(header.time),1,0);
  //scatter max_flaw_number of this inputfile
  PS::Comm::broadcast(&(header.max_flaw_number),1,0);
  //set simulation time as input time
  setTime(header.time);
  //exmine property tag of each particle is correct
  particle_number_local = sph_system.getNumberOfParticleLocal();
  for(i=0 ; i<particle_number_local ; i++){
    if(sph_system[i].property_tag<0 || sph_system[i].property_tag>=PARAM::N_MATERIAL){
      error=1;
    }
  }
  if(PS::Comm::getSum(error)>=1){
    if(PS::Comm::getRank()==0){
      std::cout << "some particles have irregal property tag! " << std::endl;
    }
    PS::Finalize();
    exit(1);
  }
  //examine all particles' ni_tot_flaw are really smaller than MAX_FLAW_NUMBER
  for(i=0 ; i<particle_number_local ; i++){
    if(sph_system[i].ni_tot_flaw > MAX_FLAW_NUMBER){
      error=1;
    }
  }
  if(PS::Comm::getSum(error)>=1){
    if(PS::Comm::getRank()==0){
      std::cout << "some particles have ni_tot_flaw that is larger than MAX_FLAW_NUMBER! " << std::endl;
    }
    PS::Finalize();
    exit(1);
  }
  //examine max_flaw_number of this input file is really the same as defined MAX_FLAW_NUMBER
  if(header.max_flaw_number!=MAX_FLAW_NUMBER){
    if(PS::Comm::getRank()==0){
      std::cout << "max_flaw_number of this input file is different from defined MAX_FLAW_NUMBER! " << std::endl;
      std::cout << "max_flaw_number of this input file:" << header.max_flaw_number << "defined MAX_FLAW_NUMBER:" << MAX_FLAW_NUMBER << std::endl;
    }
    PS::Finalize();
    exit(1);
  }
}

void set_domain(PS::ParticleSystem<RealPtcl>& sph_system, PS::DomainInfo& dinfo)
{
  PS::S64 i;
  PS::S64 particle_number_local = sph_system.getNumberOfParticleLocal();
  PS::S32 error=0;

  //initialize domain
  dinfo.initialize();
  //set periodic boundary condition as in param.h
  if(PARAM::IS_PERIODIC_BOUNDARY.x==1){
    if(PARAM::IS_PERIODIC_BOUNDARY.y==1){
      if(PARAM::IS_PERIODIC_BOUNDARY.z==1){
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
      }
      else{
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XY);
      }
    }
    else{
      if(PARAM::IS_PERIODIC_BOUNDARY.z==1){
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XZ);
      }
      else{
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_X);
      }
    }
  }
  else{
    if(PARAM::IS_PERIODIC_BOUNDARY.y==1){
      if(PARAM::IS_PERIODIC_BOUNDARY.z==1){
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_YZ);
      }
      else{
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_Y);
      }
    }
    else{
      if(PARAM::IS_PERIODIC_BOUNDARY.z==1){
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_Z);
      }
      else{
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
      }
    }
  }
  //If we use periodic boundary condition, boundaries are set as in param.h
  if(PARAM::IS_PERIODIC_BOUNDARY.x==1||PARAM::IS_PERIODIC_BOUNDARY.y==1||PARAM::IS_PERIODIC_BOUNDARY.z==1){
    dinfo.setPosRootDomain(PARAM::LOW_BOUNDARY, PARAM::HIGH_BOUNDARY);
    //examine all particles are really located in boundaries
    for(i=0 ; i<particle_number_local ; i++){
      if( sph_system[i].pos.x < PARAM::LOW_BOUNDARY.x || sph_system[i].pos.x > PARAM::HIGH_BOUNDARY.x ||
          sph_system[i].pos.y < PARAM::LOW_BOUNDARY.y || sph_system[i].pos.y > PARAM::HIGH_BOUNDARY.y ||
          sph_system[i].pos.z < PARAM::LOW_BOUNDARY.z || sph_system[i].pos.z > PARAM::HIGH_BOUNDARY.z 
          ){
        error=1;
        //printf("%d %f %f %f \n",sph_system[i].id,sph_system[i].pos.x,sph_system[i].pos.y,sph_system[i].pos.z);
      }
    }
    if(PS::Comm::getSum(error)>=1){
      if(PS::Comm::getRank()==0){
        std::cout << "some particles are out of range! " << std::endl;
      }
      PS::Finalize();
      exit(1);
    }
  }
}

//output all parameters for check
//parameters will be written in "allparameters-and-conditions.txt"
void output_parameters(char particlebin_file[])
{
  FILE* fp;
  char filename[]="allparameters-and-conditions.txt";
  time_t ct;
  struct tm *now;
  int d,n;
  int error=0;

  if(PS::Comm::getRank()==0){
    if((fp=fopen(filename,"w"))==NULL){
      printf("cannot open parameters and conditions file! \n");
      error=1;
    }
    else{
      //output start time
      ct = time(NULL);
      now = localtime(&ct);
      fprintf(fp,"simulation date : %d/%d/%d %2d:%2d:%2d \n",(now->tm_year)+1900,(now->tm_mon)+1,(now->tm_mday),(now->tm_hour),(now->tm_min),(now->tm_sec));
      //output name of inputfile
      fprintf(fp,"inputfile : %s \n",particlebin_file);
      //output spatial dimension
      fprintf(fp,"number of spatial dimension : %d \n",PARAM::Dim);
      fprintf(fp,"\n");
      
      //parallelization information
      fprintf(fp,"------parallelization informations------\n");
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
      fprintf(fp,"MPI parallelization: # of processes = %d \n",PS::Comm::getNumberOfProc());
#endif
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
      fprintf(fp,"OpenMP parallelization: # of threads / process = %d \n",PS::Comm::getNumberOfThread());
#endif
#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL
#ifndef PARTICLE_SIMULATOR_THREAD_PARALLEL
      fprintf(fp,"not parallelized \n");
#endif
#endif
      fprintf(fp,"\n");
      
      //output space information
      fprintf(fp,"------space informations------\n");
      fprintf(fp,"boundary conditions : ");
      fprintf(fp,"X = ");
      if(PARAM::IS_PERIODIC_BOUNDARY.x == 0){
        fprintf(fp,"free, ");
      }
      else if(PARAM::IS_PERIODIC_BOUNDARY.x == 1){
        fprintf(fp,"periodic, ");
      }
      else{
        std::cout << "invalid boundary flag for x direction! " << std::endl;
        error=1;
      }
      fprintf(fp,"Y = ");
      if(PARAM::IS_PERIODIC_BOUNDARY.y == 0){
        fprintf(fp,"free, ");
      }
      else if(PARAM::IS_PERIODIC_BOUNDARY.y == 1){
        fprintf(fp,"periodic, ");
      }
      else{
        std::cout << "invalid boundary flag for y direction! " << std::endl;
        error=1;
      }
      fprintf(fp,"Z = ");
      if(PARAM::IS_PERIODIC_BOUNDARY.z == 0){
        fprintf(fp,"free");
      }
      else if(PARAM::IS_PERIODIC_BOUNDARY.z == 1){
        fprintf(fp,"periodic");
      }
      else{
        std::cout << "invalid boundary flag for x direction! " << std::endl;
        error=1;
      }
      fprintf(fp,"\n");
      if(PARAM::IS_PERIODIC_BOUNDARY.x!=0||PARAM::IS_PERIODIC_BOUNDARY.y!=0||PARAM::IS_PERIODIC_BOUNDARY.z!=0){
        fprintf(fp,"range for particles : X [%e:%e] Y [%e:%e] Z [%e:%e] \n",PARAM::LOW_BOUNDARY.x,PARAM::HIGH_BOUNDARY.x,PARAM::LOW_BOUNDARY.y,PARAM::HIGH_BOUNDARY.y,PARAM::LOW_BOUNDARY.z,PARAM::HIGH_BOUNDARY.z);
      }
      fprintf(fp,"\n");
      
      //output simulation conditions
      fprintf(fp,"------simulation conditions------\n");
      if(PARAM::IS_GODUNOV_SPH==1){
        fprintf(fp,"Godunov SPH method \n");
        if(PARAM::SPACE_ORDER!=1&&PARAM::SPACE_ORDER!=2){
          std::cout << "invalid space order for riemann solver! " << std::endl;
          error=1;
        }
        fprintf(fp,"space order for riemann solver : %d \n",PARAM::SPACE_ORDER);
      }
      else if(PARAM::IS_GODUNOV_SPH==0){
        fprintf(fp,"Standard SPH method \nparameters for artificial viscosity : alpha = %.2f, beta = %.2f \n",PARAM::ALPHA_VIS,PARAM::BETA_VIS);
        if(PARAM::IS_CORRECTED_VELOCITY_GRADIENT==1){
          fprintf(fp,"correct velocity gradient for drho/dt and dSab/rho/dt: Yes \n");
        }
        else if(PARAM::IS_CORRECTED_VELOCITY_GRADIENT==0){
          fprintf(fp,"correct velocity gradient for drho/dt and dSab/rho/dt: No \n");
        }
        else{
          std::cout << "invalid switch IS_CORRECTED_VELOCITY_GRADIENT! " << std::endl;
          error=1;
        }
      }
      else{
        std::cout << "invalid switch IS_GODUNOV_SPH! " << std::endl;
        error=1;
      }
      
      if(PARAM::IS_VARIABLE_H==1){
        fprintf(fp,"variable smoothing length, Csmooth = %.2f, eta = %.2f \n",PARAM::Csmooth,PARAM::SMTH);
      }
      else if(PARAM::IS_VARIABLE_H==0){
        fprintf(fp,"constant smoothing length \n");
      }
      else{
        std::cout << "invalid switch IS_VARIABLE_H! " << std::endl;
        error=1;
      }

      if(PARAM::IS_GODUNOV_SPH==1&&PARAM::KERNEL_FUNCTION!=0){
        std::cout << "Godunov SPH can be used only with gaussian kernel! " << std::endl;
        error=1;
      }
      
      if(PARAM::KERNEL_FUNCTION==0){
        fprintf(fp,"kernel function: gaussian, search radius = %.2f\n",PARAM::SEARCH_RADIUS_H);
      }
      else if(PARAM::KERNEL_FUNCTION==1){
        fprintf(fp,"kernel function: cubic spline, search radius = 2.00\n");
      }
      else{
        std::cout << "invalid switch KERNEL_FUNCTION! " << std::endl;
        error=1;
      }
      
      if(PARAM::IS_SELF_GRAVITY!=0&&PARAM::IS_SELF_GRAVITY!=1){
        std::cout << "invalid flag IS_SELF_GRAVITY! " << std::endl;
        error=1;
      }
      else{
        fprintf(fp,"calculate self gravity : ");
        if(PARAM::IS_SELF_GRAVITY==0){
          fprintf(fp,"No \n");
        }
        else{
          fprintf(fp,"Yes \nG_CONST=%e \n",PARAM::G_CONST);
        }
      }
      fprintf(fp,"calculate deviatoric stress tensor : Yes \n");
      
      fprintf(fp,"time development method : ");
      switch(PARAM::TIME_DEVELOPMENT_METHOD){
      case 0 : fprintf(fp,"Euler method \n"); break;
      case 1 : fprintf(fp,"second order runge-kutta method \n"); break;
      case 2 : fprintf(fp,"second order leapfrog method \n"); break;
      default : printf("invalid time development method flag! \n"); error=1; break;
      }

      fprintf(fp,"density development method : ");
      switch(PARAM::DENSITY_DEVELOPMENT_METHOD){
      case 0 : fprintf(fp,"summation \n"); break;
      case 1 : fprintf(fp,"follow time evolution of density by EoC, minimum density = %f \n",PARAM::MINIMUM_DENSITY); break;
      case 2 : fprintf(fp,"follow time evolution of density by EoC (drho_dt = sumj mj * (vi - vj) *gradW), minimum density = %f \n",PARAM::MINIMUM_DENSITY); break;
      default : printf("invalid density development method flag! \n"); error=1; break;
      }
      
      fprintf(fp,"Ccfl = %.2f \n",PARAM::C_CFL);
      if(PARAM::CONSIDER_PERP_WAVE_FOR_DT!=0&&PARAM::CONSIDER_PERP_WAVE_FOR_DT!=1){
        std::cout << "invalid flag CONSIDER_PERP_WAVE_FOR_DT! " << std::endl;
        error=1;
      }
      else{
        fprintf(fp,"consider perpendicular wave for dt : ");
        if(PARAM::CONSIDER_PERP_WAVE_FOR_DT==0){
          fprintf(fp,"No \n");
        }
        else{
          fprintf(fp,"Yes \n\n");
        }
      }
      fprintf(fp,"maximum smoothing length = %e \n",PARAM::MAXIMUM_H);
      fprintf(fp,"\n");
      
      //output material dependent parameters
      fprintf(fp,"------material dependent parameters------\n");
      fprintf(fp,"the number of material type : %d \n",PARAM::N_MATERIAL);
      for( n=0 ; n<PARAM::N_MATERIAL ; n++){
        fprintf(fp,"# material type %d \n",n);
        fprintf(fp,"Equation of State : ");
        switch(PARAM::EQUATION_OF_STATE[n]){
        case 0 : fprintf(fp,"ideal gas, gamma = %f \n",PARAM::GAMMA_IDEAL[n]); break;
        case 1 : fprintf(fp,"simple elastic, Cs = %e, rho0=%f \n",PARAM::CS_ELASTIC[n],PARAM::RHO0_ELASTIC[n]); break;
        case 2 : fprintf(fp,"stiffened gas, C0 = %e, gamma0 = %f, rho0 = %f \n",PARAM::C0_STIFFENED[n],PARAM::GAMMA0_STIFFENED[n],PARAM::RHO0_STIFFENED[n]); break;
        case 3 : fprintf(fp,"tillotson, rho0 = %f, A = %e, B = %e, E0 = %e, Eiv = %e, Ecv = %e, a = %.2f, b = %.2f, alpha = %.2f, beta = %.2f \n",PARAM::RHO0_TIL[n],PARAM::A_TIL[n],PARAM::B_TIL[n],PARAM::E0_TIL[n],PARAM::EIV_TIL[n],PARAM::ECV_TIL[n],PARAM::a_TIL[n],PARAM::b_TIL[n],PARAM::ALPHA_TIL[n],PARAM::BETA_TIL[n]); break;
        case 4 : fprintf(fp,"Mie Gruneisen, rho0 = %f, C0 = %e, S = %e, Gamma = %f \n",PARAM::RHO0_MG[n],PARAM::C0_MG[n],PARAM::S_MG[n],PARAM::GAMMA_MG[n]); break;
	case 5 : fprintf(fp,"ANEOS table, table name = %s, temperature grid number = %lld, density grid number = %lld \n",PARAM::TABLE_NAME_ANE[n],PARAM::N_TEMP_GRID_ANE[n],PARAM::N_DENS_GRID_ANE[n]); break;
        default : printf("invalid EoS flag for material number %d! \n",n); error=1; break;
        }
        if(PARAM::EQUATION_OF_STATE[n]!=0&&PARAM::EQUATION_OF_STATE[n]!=1&&PARAM::IS_GODUNOV_SPH==1) fprintf(fp,"critical internal energy for Riemann solver = %e \n",PARAM::U_CRIT_RIEMANN[n]);
        if(PARAM::PLASTIC_MODEL[n]!=4){
          fprintf(fp,"shear modulus = %e \n",PARAM::MU_SHEAR[n]);
          switch(PARAM::PLASTIC_MODEL[n]){
          case 0 : fprintf(fp,"elastic only and do not consider about plasticity \n"); break;
          case 1 : fprintf(fp,"use plastic model of Benz and Asphaug 1995, yielding criterion Y0_mises = %e \n",PARAM::Y0_MISES[n]); break;
          case 2 : fprintf(fp,"use plastic model of Libersky and Petchek 1990, yielding criterion Y0_mises = %e \n",PARAM::Y0_MISES[n]); break;
          case 3 : fprintf(fp,"use plastic model including friction, fracture and pressure dependent yielding stress in Jutzi 2015 \nY0_mises = %e, u_melt = %e, mu_i_fric = %.2f, mu_d_fric = %.2f, Y0_cohesion = %e \n",PARAM::Y0_MISES[n],PARAM::U_MELT[n],PARAM::MU_I_FRIC[n],PARAM::MU_D_FRIC[n],PARAM::Y0_COHESION[n]); 
            if(PARAM::IS_FRICTION_MODEL[n]!=0&&PARAM::IS_FRICTION_MODEL[n]!=1){
              printf("invalid friction model flag! \n"); 
              error=1; 
            }
            if(PARAM::IS_FRICTION_MODEL[n]==1) fprintf(fp,"friction is considered \n"); 
            break;
          default : printf("invalid plastic model flag for material number %d! \n",n); error=1; break;
          }
        }
        else{
          fprintf(fp,"this calculation is for fluid and do not consider Sab \n");
        }
        if(PARAM::IS_FRACTURE_MODEL[n]==1){
          fprintf(fp,"fracture model is considered \nweibull k parameter = %e, weibull m parameter = %e, bulk modulus = %e \n",PARAM::K_WEIBULL[n],PARAM::M_WEIBULL[n],PARAM::K_BULK[n]);
        }
        else if(PARAM::IS_FRACTURE_MODEL[n]!=0){
          printf("invalid fracture model flag for material number %d! \n",n);
          error=1;
        }
        if(PARAM::IS_POROSITY_MODEL[n]==1){
          fprintf(fp,"porosity model is considered \nalpha0 = %f, Pe = %e, Ps = %e \n",PARAM::ALPHA_POR_0[n],PARAM::PE_POR[n],PARAM::PS_POR[n]);
        }
        else if(PARAM::IS_POROSITY_MODEL[n]!=0){
          printf("invalid porosity model flag for material number %d! \n",n);
          error=1;
        }
        fprintf(fp,"\n");
      }
      fclose(fp);
    }
  }
  PS::Comm::broadcast(&error,1,0);
  if(error==1){
    PS::Finalize(); 
    exit(1);
  }
}

void My_SPH_initialize(PS::ParticleSystem<RealPtcl>& sph_system, PS::DomainInfo& dinfo, PS::TreeForForceShort<RESULT::Dens, EPI::Dens, EPJ::Dens>::Gather_or_Scatter& dens_tree, PS::TreeForForceShort<RESULT::Gradients, EPI::Gradients, EPJ::Gradients>::Gather_or_Scatter& gradients_tree, PS::TreeForForceShort<RESULT::Elastic, EPI::Elastic, EPJ::Elastic>::Symmetry_or_Scatter& elast_tree, PS::TreeForForceLong <RESULT::Grav , EPI::Grav , EPJ::Grav >::Monopole& grav_tree, char particlebin_file[], int argc, char* argv[])
{
  PS::Initialize(argc, argv);

  if(argc!=3){
    if(PS::Comm::getRank()==0){
      std::cout << "Please input the name of this program, the name of particle informatin binary file, and the name of output directory! " << std::endl;
    }
    PS::Finalize();
    exit(1);
  }

  sph_system.initialize();
  
  read_particle_information(sph_system, particlebin_file);
  
  set_domain(sph_system, dinfo);

  dens_tree.initialize(3 * sph_system.getNumberOfParticleLocal());
  gradients_tree.initialize(3 * sph_system.getNumberOfParticleLocal());
  elast_tree.initialize(3 * sph_system.getNumberOfParticleLocal());
  if(PARAM::IS_SELF_GRAVITY){
    grav_tree.initialize(3 * sph_system.getNumberOfParticleLocal());
  }

  output_parameters(particlebin_file);
  read_aneos_table();

}
