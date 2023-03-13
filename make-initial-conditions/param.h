#pragma once

#define MAX_FLAW_NUMBER 40 //max number of flaws for one SPH particle

namespace PARAM{
  const PS::S32 Dim = 3; //number of dimension
  const PS::F64 ACC = 1.0e-12; //infinitesimal constant

  const PS::F64 G_CONST = 6.67259e-8; //gravitational constant
  const PS::S32 IS_SELF_GRAVITY = 1; //If this value is 1, we calculate self-gravity

  const PS::F64 SMTH = 1.0; //eta
  const PS::F64 C_CFL = 0.5; //CFL number
  const PS::F64 Csmooth = 1.0; //Csmooth for calculation of smoothing length
  //#define VARIABLE_SMOOTHING_LENGTH //If this macro is defined, we use variable smoothing length
#ifdef VARIABLE_SMOOTHING_LENGTH
  const PS::S32 IS_VARIABLE_H = 1; 
#else
  const PS::S32 IS_VARIABLE_H = 0; 
#endif
  const PS::F64 MAXIMUM_H = 5.0e6; //maximum smoothing length
  const PS::F64 MINIMUM_DENSITY = 0.2; //minimum density in the case of time-evolved density
  const PS::F64 SEARCH_RADIUS_H = 2.0; //This value defines how far we search SPH particles. The search radius becomes this value times smoothing length h.
  const PS::S32 IS_GODUNOV_SPH = 0; //If this value is 1, we use Godunov SPH method. Else we use Standard SPH method.
  const PS::S32 SPACE_ORDER = 2; //space order for riemann solver
  const PS::S32 TIME_DEVELOPMENT_METHOD = 2; //Select time development method. If 0, we use simple Euler method. If 1, we use second order runge-kutta method. If 2 we use second order leapfrog method.
  const PS::S32 DENSITY_DEVELOPMENT_METHOD = 1;//Select density development method. If 0, we calculate density by summation. If 1, we follow the time evolution of density by EoC. If 2, we follow the time evolution of density by EoC (drho/dt = sumj mj * (vi - vj) * gradW).
const PS::S32 IS_CORRECTED_VELOCITY_GRADIENT=1;//If this value is 1, we corrected velocity gradients for calculating drho/dt and dSab/rho/dt. Now,for dSab/rho/dt, this is only valid for SSPH.
  const PS::S32 CONSIDER_PERP_WAVE_FOR_DT=0;//If this value is 1, we use the speed of perpendicular wave for elastic dynamics (Csp = sqrt( Cs*Cs + 4.0*mu/3.0*rho )) to determine dt.
  const PS::S32 KERNEL_FUNCTION=1;//If this value is 0, we use gaussian kernel. If 1, we use cubic spline kernel. Note that Godunov SPH can be used only with gaussian kernel.

  const PS::S32vec IS_PERIODIC_BOUNDARY = PS::S32vec(0, 0, 0); //Decide which we use periodic boundary condition or not for each direction
  const PS::F64vec LOW_BOUNDARY = PS::F64vec(-1.0e7, -1.0e7, -1.0e7); //Lowest boundary for each direction
  const PS::F64vec HIGH_BOUNDARY = PS::F64vec(1.0e7, 1.0e7, 1.0e7); //Highest boundary for each direction

  const PS::S32 N_MATERIAL = 1; //the number of type of materials

  const PS::S32 PLASTIC_MODEL[N_MATERIAL]={3}; //this array decides what plastic model should be adopted to each material type.
  //0...elastic only, 1...plastic model of Benz and Asphaug 1995, 2...plastic model of Libersky and Petchek 1990, 3...plastic model including friction, fracture and pressure dependent yielding stress (Jutzi 2015), 4...fluid (Sab=0)
  const PS::S32 IS_FRACTURE_MODEL[N_MATERIAL]={1}; //If this simulation uses fracture model, this value should be 1.
  const PS::S32 IS_FRICTION_MODEL[N_MATERIAL]={1}; //If this simulation uses friction model, this value should be 1.
  const PS::S32 IS_POROSITY_MODEL[N_MATERIAL]={0}; //If this simulation uses porosity model, this value should be 1.

  const PS::F64 MU_SHEAR[N_MATERIAL] = {2.27e11};//shear modulus

  //material parameters for plastic models
  const PS::F64 Y0_MISES[N_MATERIAL]={3.5e10}; //yield stress for Mises criterion
  const PS::F64 Y0_COHESION[N_MATERIAL]={1.0e8}; //cohesion (yield stress for zero pressure)
  const PS::F64 MU_I_FRIC[N_MATERIAL]={1.0}; //coefficient of internal friction
  const PS::F64 MU_D_FRIC[N_MATERIAL]={0.839}; //coefficient of friction
  const PS::F64 U_MELT[N_MATERIAL]={3.4e10}; //melting specific internal energy

  //material parameters for fracture model
  const PS::F64 K_WEIBULL[N_MATERIAL]={4.0e29}; //Weibull k parameter
  const PS::F64 M_WEIBULL[N_MATERIAL]={9.0}; //Weibull m parameter
  const PS::F64 K_BULK[N_MATERIAL]={2.67e11}; //bulk modulus

  //material parameters for porosity model
  const PS::F64 ALPHA_POR_0[N_MATERIAL]={1.0}; //initial distension
  const PS::F64 PS_POR[N_MATERIAL]={2.0e8}; //pressure for complete compaction
  const PS::F64 PE_POR[N_MATERIAL]={1.0e6}; //threshold pressure between elastic and plastic regime

  const PS::S32 EQUATION_OF_STATE[N_MATERIAL] = {3}; //This array decide what equation of state is applied for each material.
  //0...ideal gas, 1...simple elastic, 2...stiffened gas 3...tillotson, 4...Mie Gruneisen, 5...ANEOS table
  const PS::F64 U_CRIT_RIEMANN[N_MATERIAL] = {3.04e6}; //threshold internal energy for Riemann solver between ideal and elastic eos

  //parameters for ideal gas equation of state///////////////////////
  const PS::F64 GAMMA_IDEAL[N_MATERIAL] = {1.666666666666}; //specific heat ratio
  //////////////////////////////////////////////////////////////////

  //parameters for simple elastic equation of state/////////////////
  const PS::F64 CS_ELASTIC[N_MATERIAL] = {5.0e4}; //sound speed for elastic equation of state
  const PS::F64 RHO0_ELASTIC[N_MATERIAL] = {1.5}; //reference density
  //////////////////////////////////////////////////////////////////

  //parameters for stiffened gas equation of state//////////////////
  const PS::F64 C0_STIFFENED[N_MATERIAL] = {1.0}; //sound speed for elastic part
  const PS::F64 GAMMA0_STIFFENED[N_MATERIAL] = {1.4}; //gruneisen parameter
  const PS::F64 RHO0_STIFFENED[N_MATERIAL] = {0.1}; //reference density of elastic part
  //////////////////////////////////////////////////////////////////

  //parameters for tillotson equation of state//////////////////////
  const PS::F64 RHO0_TIL[N_MATERIAL] = {2.7}; 
  const PS::F64 A_TIL[N_MATERIAL] = {2.67e11};
  const PS::F64 B_TIL[N_MATERIAL] = {2.67e11};
  const PS::F64 E0_TIL[N_MATERIAL] = {4.87e12};   
  const PS::F64 EIV_TIL[N_MATERIAL] = {4.72e10};
  const PS::F64 ECV_TIL[N_MATERIAL] = {1.82e11};  
  const PS::F64 a_TIL[N_MATERIAL] = {0.5};
  const PS::F64 b_TIL[N_MATERIAL] = {1.5};
  const PS::F64 ALPHA_TIL[N_MATERIAL] = {5.0};
  const PS::F64 BETA_TIL[N_MATERIAL] = {5.0};
  //////////////////////////////////////////////////////////////////

  //parameters for Mie Gruneisen equation of state//////////////////
  const PS::F64 RHO0_MG[N_MATERIAL] = {0.1}; 
  const PS::F64 C0_MG[N_MATERIAL] = {1.0};
  const PS::F64 S_MG[N_MATERIAL] = {1.5};
  const PS::F64 GAMMA_MG[N_MATERIAL] = {1.4};
  //////////////////////////////////////////////////////////////////

  //parameters for ANEOS table equation of state////////////////////
  const char TABLE_NAME_ANE[N_MATERIAL][256] = {"basalt_.aneos"};
  const PS::S64 N_TEMP_GRID_ANE[N_MATERIAL] = {60};
  const PS::S64 N_DENS_GRID_ANE[N_MATERIAL] = {120};
  //////////////////////////////////////////////////////////////////

  //constants for artificial viscosity for standard SPH method/////
  const PS::F64 ALPHA_VIS = 1.0;
  const PS::F64 BETA_VIS = 2.0;
  /////////////////////////////////////////////////////////////////
};
