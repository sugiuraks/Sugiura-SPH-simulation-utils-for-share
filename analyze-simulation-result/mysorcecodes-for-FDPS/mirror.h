#ifndef INCLUDE_MIRROR
#define INCLUDE_MIRROR

//fanctions that need to be called in main fanction//////////////////////////////////////////////////////////////////////////////////////////////

//read real particle space data and flags from "for-mirror.txt" and set box length and examine where real particles are in real space.
//argument...particle array
void read_mirror_datas(particle_t par[]);

//get mirror particle number
//return value...mirror particle number
int get_mirror_particle_number(void);

//locate mirror particles
//argument...particle array
void locate_mirror_particle(particle_t par[]);

//copy various variables that need to be copied
//argument...particle array
void copy_variables(particle_t par[]);

//copy gradients
//argument...particle array
void copy_gradients(particle_t par[]);

//copy d2Vdx2
void copy_d2Vdx2(particle_t par[]);

//if periodic or wall voundary and particle go out of real space for a direction,particle relocate for satisfying periodic or wall consition.
//argument...particle array
void relocate_particle_for_periodic_and_wall(particle_t par[]);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//fanctions that NOT need to be called in main fanction///////////////////////////////////////////////////////////////////////////////////////////

//return 1 if this particle should be copied to a direction.
//return value...flag 1or0  argument1...particle i's struct  argument2...direction (each X Y Z)  argument3...box number
int mirror_condition(particle_t par_i,int dir,int box_num);

//return delta of real particle's position and mirror particle's position of the direction
//return value...delta argument1...particle i's struct  argument2...direction (each X Y Z)  argument3...box number
double mirror_pos_delta(particle_t par_i,int dir,int box_num);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
