#include "SPH.h"
#include "input.h"
#include "calc-pressure-force.h"
#include "calc-self-gravity.h"
#include "time-development.h"
#include "space.h"
#include "param-mysorcecode.h"
#include "mirror.h"
#include "fracture-and-porosity-model.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static int is_mirror;

//initialization. Set boundary flag, read space data, read particle information and read mirror data.
particle_t* initialize(char particlebin_file[],int argc)
{
  int d;
  int is_free_boundary=1;
  particle_t* par;

  if(argc!=2){
    printf("please input name of this program, name of particle information binary file! \n");
    exit(FAILURE);
  }

  for(d=0 ; d<D ; d++){
    if(boundary_flags[d]!=0) is_free_boundary=0;
  }
  if(!is_free_boundary){
    is_mirror=1;
  }
  else{
    is_mirror=0;
  }

  //read space data from param.h
  read_space_data();
  
  //read particle informations
  par = read_particle_information(particlebin_file);

  //if the boundary is not free boundary, read mirror datas from param.h
  if(!is_free_boundary){
    read_mirror_datas(par);
  }

  output_parameters(particlebin_file);
  
  return(par);
}

//get is_mirror
int get_is_mirror(void)
{
  return(is_mirror);
}

//one time development roop
void one_time_development_roop(particle_t par[])
{
  tree_node_t* root_node;

  if(is_mirror) locate_mirror_particle(par); //If we use wall or periodic boundary condition, we locate mirror particles
  else          reset_space_data(par);       //If not, we reset whole space to just include all particles.
  
  root_node = make_tree(par); //make tree for neighbor searching or calculating self gravity.
  neighbor_searching(par,root_node); //neighbor searching

  if(is_variable_h) calc_h(par); //If we adopt variable smoothing length, we calculate smoothing length for each particle.
  if(is_variable_h || is_selfgravity) calc_tree_variables(root_node); //If we adopt variable smoothing length or want to use self-gravity, we should calculate variables of tree nodes.
  if(is_variable_h) calc_hmax(par,root_node); //If we adopt variable smoothing length, we should calculate max h in neighbor of each particle.

  calc_density_or_drhodt(par); //calculate density or time derivative of density.
  calc_pressures(par); //calculate pressure, sound speed and gamma with decided equation of state.
  calc_plastic_Sab(par); //calculate deviatoric stress tensor with plastic model. If we select plastic_model[] = 0 (=elastic model), then we do not do anything.
  if(use_strong_monoto) set_monotonicity_constraint(par); //If we use strong monotonicity constraint for second-order Riemann solver.
  if(is_mirror) copy_variables(par); //If we use wall or periodic boundary condition, we copy various variables of real particles to mirror particles.

  calc_gradients(par); //calculate gradient of specific volume (for cubic spline interpolation) and density, pressure, velocity (for second-order Riemann solver).
  if(is_mirror) copy_gradients(par); //If we use wall or periodic boundary condition, we copy gradients of variables.

  if(D==1){//if spatial dimension is one, we should use quintic spline interpolation for negative pressure region.
    calc_d2Vdx2_in_easy_way(par); //calculate second derivative of specific volume for quintic spline interpolation
    if(is_mirror) copy_d2Vdx2(par); //If we use wall or periodic boundary condition, we copy second derivative of specific volume.
  }

  set_dt(par); //set time stepping using Courant condition.
  if(is_selfgravity) calc_self_grav_acceleration(par,root_node); //If we want to calculate self-gravity, do so.

  calc_acceleration_dudt_and_dSab_rho_dt(par); //calculate acceleration and time derivative of internal energy by pressure force, and time derivative of Sab/rho.
  calc_dalpha_dt_and_modify_dSab_dt_for_porosity(par);//calculate time derivative of distension parameter and modify time derivative of deviatoric stress tensor for porosity model.
  
  time_development(par); //time development
  
  break_tree(root_node); //break tree

}

char dimension_c(int d)
{
  if(d == X) return('X');
  if(d == Y) return('Y');
  if(d == Z) return('Z');

  return('C');
}

char* yes_or_no(int ans)
{
  static char yes_no[5];
  if(ans){
    yes_no[0]='Y';
    yes_no[1]='e';
    yes_no[2]='s';
    yes_no[3]='\0';
  }
  else{
    yes_no[0]='N';
    yes_no[1]='o';
    yes_no[2]='\0';
  }

  return(yes_no);
}

//examine any flag is really 1 or 0. If not, print error message and exit.
void examine_flag(int flag,const char message[])
{
  if(flag!=0&&flag!=1){
    printf("invalid %s! \n",message);
    exit(FAILURE);
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
  if((fp=fopen(filename,"w"))==NULL){
    printf("cannot open parameters and conditions file! \n");
    return;
  }

  //output start time
  ct = time(NULL);
  now = localtime(&ct);
  fprintf(fp,"simulation date : %d/%d/%d %2d:%2d:%2d \n",(now->tm_year)+1900,(now->tm_mon)+1,(now->tm_mday),(now->tm_hour),(now->tm_min),(now->tm_sec));
  //output name of inputfile
  fprintf(fp,"inputfile : %s \n",particlebin_file);
  //output spatial dimension
  fprintf(fp,"number of spatial dimension : %d \n",D);
  fprintf(fp,"\n");

  //output space information
  fprintf(fp,"------space informations------\n");
  fprintf(fp,"boundary conditions : ");
  for( d=0 ; d<D ; d++){
    fprintf(fp,"%c = ",dimension_c(d));
    switch (boundary_flags[d]){
    case 0 : fprintf(fp,"free "); break;
    case 1 : fprintf(fp,"wall "); break;
    case 2 : fprintf(fp,"periodic "); break;
    default : printf("invalid boundary flag! \n"); exit(FAILURE);
    }
  }
  fprintf(fp,"\n");
  fprintf(fp,"range for real particles : ");
  for( d=0 ; d<D ; d++){
    fprintf(fp,"%c [%.4e:%.4e] ",dimension_c(d),init_real_Rmin[d],init_real_Rmin[d]+init_real_Range[d]);
  }
  fprintf(fp,"\n");
  if(is_mirror){
    fprintf(fp,"range for mirror particles : ");
    for( d=0 ; d<D ; d++){
      fprintf(fp,"%c [%.4e:%.4e] ",dimension_c(d),init_Rmin[d],init_Rmin[d]+init_Range[d]);
    }
    fprintf(fp,"\n");
  }
  fprintf(fp,"\n");

  //output simulation conditions
  fprintf(fp,"------simulation conditions------\n");
  if(is_riemann==1){
    fprintf(fp,"Godunov SPH method \n");
    if(space_order!=1&&space_order!=2){
      printf("invalid space order for riemann solver! \n");
      exit(0);
    }
    fprintf(fp,"space order for riemann solver : %d \n",space_order);
    examine_flag(use_strong_monoto,"strong monotonicity constraint");
    fprintf(fp,"strong monotonicity constraint : %s \n",yes_or_no(use_strong_monoto));
  }
  else if(is_riemann==0){
    fprintf(fp,"Standard SPH method \nparameters for artificial viscosity : alpha = %.2f, beta = %.2f \n",alpha_vis,beta_vis);
  }
  else{
    printf("invalid swich is_riemann! \n");
    exit(FAILURE);
  }
  if(is_variable_h==1){
    fprintf(fp,"variable smoothing length, Csmooth = %.2f, eta = %.2f \n",Csmooth,eta);
  }
  else if(is_variable_h==0){
    fprintf(fp,"constant smoothing length \n");
  }
  else{
    printf("invalid swich is_variable_h! \n");
    exit(FAILURE);
  }
  fprintf(fp,"neighbor searching range : %d h \n",Ncell);
  examine_flag(is_ienergy_needed,"internal energy flag");
  fprintf(fp,"calculate internal energy : %s \n",yes_or_no(is_ienergy_needed));
  examine_flag(is_cohesion,"deviatoric stress tensor flag");
  fprintf(fp,"calculate deviatoric stress tensor : %s \n",yes_or_no(is_cohesion));
  switch(is_density_devel){
  case 0 : fprintf(fp,"calculate both rho and rho_eos with summation \n"); break;
  case 1 : fprintf(fp,"calculate rho with summation and rho_eos with time development \n"); break;
  case 2 : fprintf(fp,"calculate both rho and rho_eos with time development \n"); break;
  default : printf("invalid density development flag! \n"); exit(FAILURE);
  }
  examine_flag(is_selfgravity,"self gravity flag");
  fprintf(fp,"calculate self gravity : %s \n",yes_or_no(is_selfgravity));
  if(is_selfgravity) fprintf(fp,"theta = %.2f, g_const = %e\n",theta,g_const);
  fprintf(fp,"time development method : ");
  switch(time_development_method){
  case 0 : fprintf(fp,"Euler method \n"); break;
  case 1 : fprintf(fp,"second order runge-kutta method \n"); break;
  case 2 : fprintf(fp,"second order leapfrog method \n"); break;
  default : printf("invalid time development method flag! \n"); exit(FAILURE);
  }
  fprintf(fp,"Ccfl = %.2f \n",Ccfl);
  fprintf(fp,"minimum density = %f \n",rho_min_crit);
  fprintf(fp,"maximum smoothing length = %e \n",h_regulator);
  fprintf(fp,"\n");
  
  //output material dependent parameters
  fprintf(fp,"------material dependent parameters------\n");
  fprintf(fp,"the number of material type : %d \n",N_material);
  for( n=0 ; n<N_material ; n++){
    fprintf(fp,"# material type %d \n",n);
    if(is_cohesion){
      fprintf(fp,"mu_shear = %e \n",mu_shear[n]);
      switch(plastic_model[n]){
      case 0 : fprintf(fp,"elastic only and do not consider about plasticity \n"); break;
      case 1 : fprintf(fp,"use plastic model of Benz and Asphaug 1995, yielding criterion Y0_mises = %e \n",Y0_mises[n]); break;
      case 2 : fprintf(fp,"use plastic model of Libersky and Petchek 1990, yielding criterion Y0_mises = %e \n",Y0_mises[n]); break;
      case 3 : fprintf(fp,"use plastic model including friction, fracture and pressure dependent yielding stress in Jutzi 2015 \nY0_mises = %e, u_melt = %e, mu_i_fric = %.2f, mu_d_fric = %.2f, Y0_cohesion = %e \n",Y0_mises[n],u_melt[n],mu_i_fric[n],mu_d_fric[n],Y0_cohesion[n]); examine_flag(is_friction_model[n],"friction model flag"); if(is_friction_model[n]) fprintf(fp,"friction is considered \n"); break;
      default : printf("invalid plastic model flag for material number %d! \n",n); exit(FAILURE);
      }
    }
    examine_flag(is_fracture_model[n],"fracture model flag");
    if(is_fracture_model[n]==1) fprintf(fp,"fracture model is considered \nweibull k parameter = %e, weibull m parameter = %e, bulk modulus = %e \n",k_weibull[n],m_weibull[n],K_bulk[n]);
    examine_flag(is_porosity_model[n],"porosity model flag");
    if(is_porosity_model[n]==1) fprintf(fp,"porosity model is considered \nalpha0 = %f, Pe = %e, Ps = %e \n",alpha_por_0[n],Pe_por[n],Ps_por[n]);
    fprintf(fp,"Equation of State : ");
    switch(equation_of_state[n]){
    case 0 : fprintf(fp,"ideal gas, gamma = %f \n",gamma_ideal[n]); break;
    case 1 : fprintf(fp,"simple elastic, Cs = %e, rho0=%f \n",Cs_elastic[n],rho0_elastic[n]); break;
    case 2 : fprintf(fp,"stiffened gas, C0 = %e, gamma0 = %f, rho0 = %f \n",C0_stiffened[n],gamma0_stiffened[n],rho0_stiffened[n]); break;
    case 3 : fprintf(fp,"tillotson, rho0 = %f, A = %e, B = %e, E0 = %e, Eiv = %e, Ecv = %e, a = %.2f, b = %.2f, alpha = %.2f, beta = %.2f \n",rho0_til[n],A_til[n],B_til[n],E0_til[n],Eiv_til[n],Ecv_til[n],a_til[n],b_til[n],alpha_til[n],beta_til[n]); break;
    case 4 : fprintf(fp,"Mie Gruneisen, rho0 = %f, C0 = %e, S = %e, Gamma = %f \n",rho0_MG[n],C0_MG[n],S_MG[n],Gamma_MG[n]);
    default : printf("invalid EoS flag for material number %d! \n",n); exit(FAILURE);
    }
    if(equation_of_state[n]!=0&&equation_of_state[n]!=1&&is_riemann==1) fprintf(fp,"critical internal energy for Riemann solver = %e \n",u_crit_Riemann[n]);
  }

  fclose(fp);
}
