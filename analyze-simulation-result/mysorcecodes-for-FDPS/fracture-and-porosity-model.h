#ifndef FRACTURE_MODEL
#define FRACTURE_MODEL

//calculate effective strain for all particles
//argument...particle array
void calc_effective_strain(particle_t par[]);

//time developmen of D for i-th particle
//1st argument...pointer to particle i struct, 2nd argument...time stepping
void time_development_of_D_i(particle_t* par_iP,double dt);

//time developmen of D for i-th particle for second order rungekutta method
//1st argument...pointer to particle i struct, 2nd argument...time stepping, 3rd argument...time step sign (if even time step, this value is 1. else -1)
void time_development_of_D_i_for_second_order_rungekutta(particle_t* par_iP,double dt,int time_step_sign);

//time developmen of D for i-th particle for leapfrog method
//1st argument...pointer to particle i struct, 2nd argument...time stepping
void time_development_of_D_i_for_leapfrog(particle_t* par_iP,double dt);

//calculate time derivative of distension parameter and modify time derivative of deviatoric stress tensor for porosity model.
//argument...particle array
void calc_dalpha_dt_and_modify_dSab_dt_for_porosity(particle_t par[]);

#endif
