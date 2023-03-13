#ifndef MAKE_INITIAL
#define MAKE_INITIAL

void make_initial_condition_file(int particle_number,int dim,double time,particle_t par[],char filename[],int remaining_flag[],int mode);

//assignment of explicit flaws to particles
//argument...particle array
void flaw_assignment(particle_t par[],int particle_number,double volume);

#endif
