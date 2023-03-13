#ifndef PARAM
#define PARAM

//parameters!////////////////////////////////////////////////////////////////////////////////////////////////////////////
static const double g_const=6.67259E-8;//gravitational constant

static const int boundary_flags[D]={0,0,0};//each direction's flags for periodic or wall boundary.If 0, free boundary.If 1, wall boundary.If 2, periodic boundary.
static const double init_Range[D]={2.0e7,2.0e7,2.0e7};//initial range for whole space
static const double init_Rmin[D]={-1.0e7,-1.0e7,-1.0e7};//initial Rmin for whole space
static const double init_real_Range[D]={2.0e7,2.0e7,2.0e7};//initial range for real particles
static const double init_real_Rmin[D]={-1.0e7,-1.0e7,-1.0e7};//initial Rmin for real particles

static const int space_order=2;//space order for Riemann solver
static const int use_strong_monoto=0;//If we want to use strong monotonicity constraint, this value should be 1.
static const double Ccfl=0.5;//CFL number
static const int roop_number=200000;//roop number for time roop
static const int Ncell=3;//this number determines the ratio between neighbor searching renge and the smoothing length
static const int is_variable_h=0;//If we want to use variable smoothing length, this value should be 1. If 0, we use constant smoothing length.
static const double Csmooth=1.0;//Csmooth
static const double theta=0.5;//criterion angle for calculation of gravity using tree 
static const double eta=1.0;//ratio between the smoothing length and average particle spacing
static const int is_selfgravity=1;//if this number is 1, we calculate self gravity
static const int is_riemann=0;//if this number is 1, we use Godunov SPH method. else we use standard SPH method
static const int is_density_devel=2;//If this is 0, both rho and rho_eos are calculated by summaiton. If 1, only rho_eos time develops. If 2, both rho and rho_eos time develop.
static const int is_ienergy_needed=1;//If time change of internal energy is needed, this value become 1. If not needed, become 0.
static const int is_cohesion=1;//If this simulation calculate the time evolution of deviatoric stress tensor, this calue should be 1.

static const int time_development_method=2;//Select time development method. If 0, we use simple Euler method. If 1, we use second order runge-kutta method. If 2 we use second order leapfrog method.

#define N_material 1 //the number of material type. If this simulation includes only one material, this should be 1. 
//For parameters which is in form of array with N_material (e.g. mu_shear[N_material]), parameter of n-th element means parameter for n-th material.

static const int is_fracture_model[N_material]={1};//If this simulation uses fracture model, this value becomes 1.
static const int is_friction_model[N_material]={1};//If this simulation includes friction between granular material, this value becomes 1.
static const int plastic_model[N_material]={3};//This array decides what plastic model should be adopted for each material type. 
//0...elastic only, 1...plastic model of Benz and Asphaug 1995, 2...plastic model of Libersky and Petchek 1990, 3...plastic model including friction, fracture and pressure dependent yielding stress (Jutzi 2015)
static const int is_porosity_model[N_material]={0};//If this simulation uses porosity model, this value should be 1.

//mateial parameters for elastic-plastic model (plastic model is Benz and Asphaug 1995 or Libersky and Petchek 1990)
static const double mu_shear[N_material]={2.27e11};//shear modulus
static const double Y0_mises[N_material]={3.50e10};//yield stress for Mises criterion

//material parameters for plastic model including friction, fracture and pressure dependent yielding stress (Jutzi 2015)
static const double u_melt[N_material]={3.4E10};//melting specific internal energy
static const double mu_i_fric[N_material]={1.0};//coefficient of internal friction
static const double mu_d_fric[N_material]={0.849};//coefficient of friction
static const double Y0_cohesion[N_material]={1.0E9};//cohision

//material parameters for fracture model
static const double k_weibull[N_material]={4.0E29};//Weibull k parameter
static const double m_weibull[N_material]={9.0};//Weibull m parameter
static const double K_bulk[N_material]={2.67e11};//bulk modulus

//material parameters for porosity model
static const double alpha_por_0[N_material]={1.0};//initial distension
static const double Ps_por[N_material]={1.0e3};//pressure for complete compaction
static const double Pe_por[N_material]={1.0e5};//threshold pressure between elastic and plastic regime


static const int equation_of_state[N_material]={3};//This array decide what equation of state is applied for each material.
//0...ideal gas, 1...simple elastic, 2...stiffened gas 3...tillotson, 4...Mie Gruneisen

//parameters for ideal gas equation of state//
static const double gamma_ideal[N_material]={1.4};//specific heat ratio for ideal gas equation of state
//////////////////////////////////////////////

//parameters for elastic equation of state////
static const double Cs_elastic[N_material]={1.0};//sound speed for elastic equation of state
static const double rho0_elastic[N_material]={0.1};//reference density 
//////////////////////////////////////////////

//parameters for stiffened gas equation of state//
static const double C0_stiffened[N_material]={541.67};//sound speed of elastic part 
static const double gamma0_stiffened[N_material]={1.4};//grunuisen paramter
static const double rho0_stiffened[N_material]={0.9152};//reference density of elastic part
//////////////////////////////////////////////////

//tillotson parameters(cgs units [Benz and Asphaug 1999])//////////////////////////////////////////////
static const double rho0_til[N_material]={2.7};
static const double A_til[N_material]={2.67E11};
static const double B_til[N_material]={2.67E11};
static const double E0_til[N_material]={4.87E12};
static const double Eiv_til[N_material]={4.72E10};
static const double Ecv_til[N_material]={1.82E11};
static const double a_til[N_material]={0.5};
static const double b_til[N_material]={1.5};
static const double alpha_til[N_material]={5.0};
static const double beta_til[N_material]={5.0};
///////////////////////////////////////////////////////////////////////////////////////////////////////

//parameters for Mie Gruneisen EoE////////////////////////////////////////////////////////////////////
static const double rho0_MG[N_material]={0.1};
static const double C0_MG[N_material]={1.0};
static const double S_MG[N_material]={1.5};
static const double Gamma_MG[N_material]={1.4};
////////////////////////////////////////////////////////////////////////////////////////////////////

//critical internal energy for Riemann solver
static const double u_crit_Riemann[N_material]={3.04E10};

//constants for artificial viscosity//////////////
static const double alpha_vis=1.0;
static const double beta_vis=2.0;
//////////////////////////////////////////////////

static const double search_region_monoto=3.0;//search region for modified monotonicity constraint
static const double rho_min_crit=0.1;//mininum density
static const double h_regulator=5.0e6;//maximum smoothing length

#endif
