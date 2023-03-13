#ifndef INCLUDE_TIME_DEVELOPMENT
#define INCLUDE_TIME_DEVELOPMENT

//get dt
double get_dt(void);

//set dt using particle information
void set_dt(particle_t par[]);

//calculate time development. This function only select time development method.
void time_development(particle_t par[]);

//calculate time development of all particles and check internal energy not to be negative.
void time_development_for_Euler(particle_t par[]);

//calculate time development of all particles for scond order rungekutta method
void time_development_for_second_order_rungekutta(particle_t par[]);

//initial kick for leapfrog method. this function should be called before time roop begin.
void initial_kick_for_leapfrog(particle_t par[]);

//time development for leapfrog method
void time_development_for_leapfrog(particle_t par[]);

//set t as t0
void set_t_as_t0(double t);

//get t
double get_t(void);

//set dt by argument of this function
void set_dt_as_const(double dt_const);

#endif
