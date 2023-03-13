#include "SPH.h"
#include "input.h"
#include "param-mysorcecode.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static double Range[D]; //space range 
static double Rmin[D]; //space minimam

//read space data from param.h
void read_space_data(void)
{
	int d; //now dimension

	for(d=0;d<D;d++){
	  Range[d] = init_Range[d];
	}
	
	for(d=0;d<D;d++){
	  Rmin[d] = init_Rmin[d];
	}
}

//reset space data 
void reset_space_data(particle_t par[])
{
	double center[D]; //center of space
	double delta2,deltamax=0.0;
	double temp;
	int d,i;
	int particle_number=get_particle_number();

	//calculating center!!
	for(d=0;d<D;d++){
		temp=0.0;
		#ifdef _OPENMP
		#pragma omp parallel for private(i) reduction(+:temp)
		#endif
		for(i=0;i<particle_number;i++){
			temp+=par[i].r[d];
		}
		center[d]=temp/particle_number;
	}

	//calculating deltamax!!
	for(i=0;i<particle_number;i++){
		delta2=0.0;
		for(d=0;d<D;d++){
			delta2+=(center[d]-par[i].r[d])*(center[d]-par[i].r[d]);
		}
		if(deltamax<sqrt(delta2)) deltamax=sqrt(delta2);
	}

	//reset space data!!
	for(d=0;d<D;d++){
		Range[d]=2.0*deltamax;
		Rmin[d]=center[d]-deltamax;
	}
}

//get Range data
void get_Range(double *RangeP)
{
	int d;
	for(d=0;d<D;d++){
		RangeP[d]=Range[d];
	}
}

//get Rmin data
void get_Rmin(double *RminP)
{
	int d;
	for(d=0;d<D;d++){
		RminP[d]=Rmin[d];
	}
}
