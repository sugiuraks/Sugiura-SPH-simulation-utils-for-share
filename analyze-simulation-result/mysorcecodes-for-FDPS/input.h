#ifndef INCLUDE_INPUT
#define INCLUDE_INPUT

//read particle number,dimension and particle information from given filename's binary file and allocate particle array and return pointer to particle array.
//return value...pointer to particle array  argument...filename of particle information binary file
particle_t* read_particle_information(char filename[]);

//return particle number
int get_particle_number(void);

//get max_flaw_number
int get_max_flaw_number(void);

//free particle array
void free_particle_array(particle_t par[]);

//count memory. argument is added to memory.
void count_memory(int count);

//return memory
unsigned long long get_memory(void);

#endif
