#ifndef INCLUDE_CALC_PRESSURE_FORCE
#define INCLUDE_CALC_PRESSURE_FORCE

//struct for calculation of acceleration and dudt 
typedef struct{
	double p_riemann;
	double v_riemann;
	double Vij2_hi;
	double Vij2_hj;
	double common;
	int num;
	double Sijs[D][D];
	double vijs[D];
	double delta;
	double dif[D];
	double dWdx_hi[D];
	double dWdx_hj[D];
	double ss;
} inter_t;

//fanctions that need to be called mainly

//fanctions for neighbor searching////////////////////////////////////////////////////////////

//seach neighbor particles
//argument 1....particle array  argument 2....pointer to root node
void neighbor_searching(particle_t par[],tree_node_t *root_node);

//fanctions for calculating hydro variables///////////////////////////////////////////////////

//set h ragulator as max of h
//argument...particle array
//void set_h_regulator(particle_t par[]);

//calculate smoothing length and J
//argument1...particle array  
void calc_h(particle_t par[]);

//calculate hmax and J correspond to hmax
//argument1...particle array  argument2...pointer to tree's root node.
void calc_hmax(particle_t par[],tree_node_t *root_node);

//calculate density or time derivative of density
//argument1...particle array  
void calc_density_or_drhodt(particle_t par[]);

//calculate pressure, sound speed and gamma. This function only select appropriate EoS for each particle.
//argument...particle array
void calc_pressures(particle_t par[]);

//calculate gradients of density,pressure,1/rho and velosity
//argument1...particle array 
void calc_gradients(particle_t par[]);

//calculate d2Vdx2 in easy way
void calc_d2Vdx2_in_easy_way(particle_t par[]);

//calculate density at any point that have position r and smoothing length h
double calc_rho_at_r(double r[],double h,tree_node_t *root_node);


//fanctions for calculate pressure force and internal enegy change rate/////////////////////////////

//calculate acceleration of pressure force, internal enrgy change rate and change rate of deviatric stress tensor. If is_riemann flag is 1,calculate with godnov method.if 0,calculate with standard SPH.
//argument1...particlearray  
void calc_acceleration_dudt_and_dSab_rho_dt(particle_t par[]);

//calculate plastic Sab 
//argument...particle array
void calc_plastic_Sab(particle_t par[]);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




//fanctions that NOT need to be called in main fanction////////////////////////////////////////////////////////////////////////////////////////////////////

//fanctions for neighbor seaching///////////////////////////////////////////////////////////////////

//seach neighbor particles recursively
//argument 1...pointer to i particle  argument 2....pointer to a tree node
void neighbor_searching_recursively(particle_t *par_i,tree_node_t *n_node);


//fanctions for calculating hydro variables/////////////////////////////////////////////////////////

//calculate pressure and sound speed using ideal gas equation of state
//argument...pointer to i-th particle
void calc_p_and_Cs_ideal_gas(particle_t* par_i);

//calculate delp_delrho and delp_delu for ideal gas
//argument1...density, argument2...internal energy, argument3...property tag of i-th particle, argument4...pointer to delp_delrho, argument5...pointer to delp_delu
void calc_delp_delrho_and_delp_delu_ideal_gas(double rho,double u,int p_tag,double* delp_delrho_P,double* delp_delu_P);

//calculate pressure,sound speed and specific heat ratio from internal energy and density for tillotson equation of state.
//argument1...pointer to i-th particle
void calc_p_Cs_and_gamma_tillotson(particle_t* par_i);

//calculate pressure, Cs and gamma of stiffened gas equation of state
//argument1...pointer to i-th particle
void calc_p_Cs_and_gamma_stiffened_gas_EoS(particle_t* par_i);

//calculate delp_delrho and delp_delu for stiffened gas EoS
//argument1...density, argument2...internal energy, argument3...property tag of i-th particle, argument4...pointer to delp_delrho, argument5...pointer to delp_delu
void calc_delp_delrho_and_delp_delu_stiffened_gas_EoS(double rho,double u,int p_tag,double* delp_delrho_P,double* delp_delu_P);

//calculate pressure using elastic equation of state (P=Cs^2(rho-rho0eos)).
//argument1...pointer to i-th particle
void calc_p_elastic_eos(particle_t* par_i);

//calculate delp_delrho and delp_delu for elastic EoS
//argument1...density, argument2...internal energy, argument3...property tag of i-th particle, argument4...pointer to delp_delrho, argument5...pointer to delp_delu
void calc_delp_delrho_and_delp_delu_elastic_eos(double rho,double u,int p_tag,double* delp_delrho_P,double* delp_delu_P);

//calculate pressure, Cs and gamma of Mie Gruneisen EoS
//argument1...pointer to i-th particle
void calc_p_Cs_and_gamma_Mie_Gruneisen_EoS(particle_t* par_i);

//calculate delp_delrho and delp_delu for Mie Gruneisen EoS
//argument1...density, argument2...internal energy, argument3...property tag of i-th particle, argument4...pointer to delp_delrho, argument5...pointer to delp_delu
void calc_delp_delrho_and_delp_delu_Mie_Gruneisen_EoS(double rho,double u,int p_tag,double* delp_delrho_P,double* delp_delu_P);

//judge where particle's h overlap to n_node or not.
//return velue...if orverlap 1,if not 0. argument1...particle position  argument2...n_node's center  argument3...n_node's box size  argument4...h for overlap
int is_overlap(double par_pos[],double node_center[],double box_size[],double h_overlap,double par_i_Csmooth);

//calculate hmax recursively
//argument1...pointer to particle i  argument2...pointer to tree's now_node  argument3...array of hmax of tree nodes that have particle i
void calc_hmax_recursive(particle_t *par_i,tree_node_t *n_node,double hmax_of_per_node[]);

//calculate nabla kernel fanction and set dWdx array.It is gaussian kernel varsion.
//argument1...square of distance between particle i and j  argument2...distance array between i and j of each direction  argument3...smoothing length  argument4...array that store nabla kernel
void nablaW_gauss(double delta2,double dif[],double h,double dWdx[]);


//fanctions for calculating acceleration and dudt///////////////////////////////////////////////////

//calculate aij using standard SPH
//argument1...pointer to particle i  argument2...particle j's struct   argument3...pointer to common part value of aij and dudtij
void calc_aij_standard(particle_t *par_iP,particle_t par_j,double *commonP);

//calculate dudtij using standard SPH
//return value dudtij  argument1...particle i's struct   argument2...particle j's struct  argument3...array of next time step's velosity of i particle  argument4...common part value of aij and dudtij
double calc_dudtij_standard(particle_t par_i,particle_t par_j,double next_v[],double common);

//calculate dSabdt standard
//argument 1...pointer to particle i. argument2...particle j's struct
void calc_dSab_rho_dt_standard(particle_t *par_iP, particle_t par_j);

//calculate Vij2 and s star using linear interpolation
//argument1...particle i's struct  argument2...particle j's struct  argument3...pointer to interaction struct of i and j particle  argument4...pointer to s star  argument5...diferrence between i and j  argument6...unit vector eij
void calc_Vij2_using_linear(particle_t par_i,particle_t par_j,inter_t *inter_N,double *ssP,double delta,double eij[]);

//calculate Vij2 and s star using cubic spline interpolation
//argument1...particle i's struct  argument2...particle j's struct  argument3...pointer to interaction struct of i and j particle  argument4...pointer to s star  argument5...diferrence between i and j  argument6...unit vector eij
void calc_Vij2_using_cubic_spline(particle_t par_i,particle_t par_j,inter_t *inter_N,double *ssP,double delta,double eij[]);

//calculate Vij2 and s star using quintic spline interpolation
//argument1...particle i's struct  argument2...particle j's struct  argument3...pointer to interaction struct of i and j particle  argument4...pointer to s star  argument5...diferrence between i and j  argument6...unit vector eij
void calc_Vij2_using_quintic_spline(particle_t par_i,particle_t par_j,inter_t *inter_N,double *ssP,double delta,double eij[]);

//solve riemann problem
//argument1~8...initial condition of riemann problem  argument9...pressure of riemann solver  argument10...velosity of riemann solver argument11...key that become 1 if invalid riemann solver
void solve_riemann_problem_for_ideal_gas_eos(double p1,double rho1,double v1,double gamma1,double p2,double rho2,double v2,double gamma2,double* p_riemannP,double* v_riemannP);

//solve riemann problem for elastic EoS of P=Cs^2(rho-rho0eos)
void solve_riemann_problem_for_elastic_eos(double p1,double rho1,double v1,double Csi,double p2,double rho2,double v2,double Cs2,double* p_riemannP,double* v_riemannP);

//strengthen monotonicity constraint for riemann solver
void set_monotonicity_constraint(particle_t par[]);

//calculate riemann solver
//argument1...particle i's struct  argument2...particle j's struct  argument3...pointer to interaction struct of i and j particle  argument4...s star  argument5...unit vector eij  argument6...difference between i and j  argument7...delta t
void calc_riemann_solver(particle_t par_i,particle_t par_j,inter_t *inter_N,double ss,double eij[],double delta,double dt);

//calculate aij using godnov method
//argument1...pointer to particle i  argument2...particle j's struct  argument3...pointer to interaction struct of i and j particle  argument4...delta t
void calc_aij_riemann(particle_t *par_iP,particle_t par_j,inter_t *inter_N,double dt);

//calculate dudtij using godnov method
//return value...dudtij  argument1...particle i's struct  argument2...particle j's struct  argument3...array of time centered velosity of i particle  argument4...interaction struct of i and j particle
double calc_dudtij_riemann(particle_t par_i,particle_t par_j,double tcenter_v[],inter_t inter_n);

//calculate dSabdt using godunov SPH
void calc_dSab_rho_dt_riemann(particle_t *par_iP, particle_t par_j, inter_t inter_n);

//set Csmooth for 2D and cubic spline case in negative pressure to suppress the tensile instability.
void set_Csmooth_for_2D_cubic(particle_t par[]);



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
