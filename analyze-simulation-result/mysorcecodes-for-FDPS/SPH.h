#ifndef INCLUDE_SPH
#define INCLUDE_SPH

#define PI 3.141592653589793238 //the ratio of the circumference of a circle to its diameter (‰~Žü—¦)

#define X 0 //coodinate constant for each X Y Z
#define Y 1
#define Z 2 

#define D 3 //dimension

#if D == 3
#define N_CHILD 8 //number of tree node's children
#elif D == 2
#define N_CHILD 4
#elif D == 1
#define N_CHILD 2
#endif

#define MAX_TREE_LEVEL 40 //max of tree level

#define INTER 150 //one unit of inter particle number

#define MAX_MIRROR_NUM 0 //max of mirror particle number

#define FAILURE 1
#define SUCCESS 0

#define ACC 1.0E-6 //infinitesimal constant

#define EMPTY -1 //if hydro box is empty or a particle's next is vacant,this macro is settled.

#define NODE_PAR_MAX 5 //max particle number that one node can have

#define DEBUG 0 //if this value is 1,it is debug mode.

//definition of particle's structure
typedef struct part_tag{
  int ID;
  double mass; //mass
  double r[D]; //position
  double r_old[D]; //position of previous time step
  double v[D]; //velocity
  double v_old[D]; //velocity of previous time step
  double rho; //density
  double rho_eos;//density for equation of state
  double rho_old; //density of previous time step
  double p; //pressure
  double p_old; //pressure of previous time step
  double u; //internal energy
  double u_old; //internal energy of previous time step
  double u_max; //maximum internal energy that this particle experienced
  double gamma; //specific heat ratio
  double a[D]; //acceleration
  double a_old[D]; //acceleration of previous time step
  double dudt; //internal energy change rate
  double dudt_old; //internal energy change reta of previous time step
  double nablaV[D]; //gradient of 1/rho
  double nablar[D]; //gradient of density
  double nablap[D]; //gradient of pressure
  double dvdx[D][D]; //gradient of each direction's velosity
  double h; //smoothing length
  double hmax; //max of smoothing length in the vicinity
  double Cs; //sound speed of this particle's position
  double g[D]; //gravitational acceleration
  double g_old[D]; //gravitational acceleration of previous time step
  double grav_pot;//gravitational potential
  struct part_tag **inter_par; //pointer array to inter particle
  int inter_N; //number of inter particle
  int inter_Nmax; //max number of inter particle;
  struct tree_node_tag *tree_parent; //pointer to reef tree node that this particle join to
  struct part_tag *tree_next; //pointer to next particle that join to same tree node
  double d2Vdx2[D][D]; //second differential of specific volume in one dimension
  int monoto_flag; //the flag for monotonicity condition. If the sign of the differencial of i particle is different from neighbor particles, this value become 1. 
  double Sab[D][D]; //deviatric stress tensor
  double Sab_rho[D][D]; //deviatoric stress tensor / density
  double Sab_rho_old[D][D]; //deviatoric stress tensor / density of previous time step
  double dSab_rho_dt[D][D]; //change rate of deviatric stress tensor / density
  double dSab_rho_dt_old[D][D]; //change rate of deviatric stress tensor / density of previous time step
  double drhodt; //change rate of density
  double drhodt_old; //change rate of density of previous time step
  double Csmooth; //Csmooth for this particle
  double damage;//damage parameter D
  double D_cbrt;//cubic root of damage parameter D
  double D_cbrt_old;//cubic root of damage parameter D of previous time step
  double dD_cbrt_dt;//time change rate of cubic root of D
  double dD_cbrt_dt_old;//time change rate of cubic root of D of previous time step
  double* eps_ij_act_flaw;//flaw activation threshold strain
  long int ni_tot_flaw;//the number of flaws that are assined to this particle
  double eps_i_eff;//effective strain of this particle
  int property_tag;//property tag. If 0, this particle is for sphere. If 1, this particle is for plate. If 2, this particle is for bottom of plate.
  int low_density_flag;//if density of this particle is lower than minimum density, this flag becomes 1, and do not use riemann solver for this particle.
  double alpha_por;//distension parameter
  double alpha_por_old;//distension patemeter of previous time step
  double dalpha_dt;//time derivative of distension patameter
  double dalpha_dt_old;//time derivative of distension parameter of previous time step
  double f_test;
} particle_t;

//definition of tree node's structure
typedef struct tree_node_tag{
  double CoM[D]; //center of mass of all particles this tree node contain
  double center[D]; //center position of this tree node
  double box_size[D]; //side length of this node's box
  double mass; //all mass of particles this tree node contain
  double hmax; //max smoothing length of all particles this tree node contain
  struct tree_node_tag *parent; //pointer to parent node of this node
  struct tree_node_tag *child[N_CHILD]; //pointer to children of this node
  struct part_tag *first_particle; //pointer to first particle that join to this node
  int is_leaf; //if this value is 1,this node is reaf node.if 0,not read node
  int level; //this tree node's node level;
  int C[3]; //each value show this node's address 
  int num_par; //number of particle this node have
} tree_node_t;

typedef struct{
  double r[D];
  double rho;
  double p;
  double h;
} cell_t;

//initialization. Set boundary flag, read space data, read particle information, set h regulator and read mirror data.
particle_t* initialize(char particlebin_file[],int argc);

//get is_mirror
int get_is_mirror(void);

//one time development roop
void one_time_development_roop(particle_t par[]);

//output all parameters for check
void output_parameters(char particlebin_file[]);

#endif
