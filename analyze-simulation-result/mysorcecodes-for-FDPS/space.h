#ifndef INCLUDE_SPACE
#define INCLUDE_SPACE

//read space data from 'space.txt'
void read_space_data(void);

//reset space data
//argument...particleinfo array
void reset_space_data(particle_t par[]);

//get Range data
//argument...Range array
void get_Range(double *RangeP);

//get Rmin data
//argument...Rmin array
void get_Rmin(double *RminP);

#endif
