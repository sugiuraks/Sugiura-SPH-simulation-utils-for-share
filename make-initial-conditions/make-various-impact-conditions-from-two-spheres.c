//////////////////////////////////////////////////////////
//make-sphere-from-relaxed-condition.cなどで作った2つの球を衝突させる初期条件を,
//衝突速度と衝突角度を一定の間隔で増やしながらたくさん作るプログラムです.
//(コードをうまく改造すれば他の条件などもたくさん作れるようになると思います)
//衝突速度と衝突角度の増やし方はfor文の頭で指定します.
//他, main関数の頭の情報は適切に与える必要があります.

//使い方/////////////////
//(Makefileはsrc, my-srcの場所やコンパイラなど適切に書き換えてください)
//このプログラムとMakefileがあるディレクトリで
//$ make makevariousinitials
//$ ./make-various-initials.out
/////////////////////////////////////////////////////////
#include <particle_simulator.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "mathfunc.h"
#include "param.h"
#include "kernel.h"
#include "class.h"
#include "force.h"

#define D 3
#define X 0
#define Y 1
#define Z 2

PS::matrix set_rotation_matrix_along_z(double theta)
{
  PS::matrix M_rotate = PS::matrix(0, 0, 0,
	   0, 0, 0,
	   0, 0, 1 );
  M_rotate.xx =  cos(theta);
  M_rotate.xy = -sin(theta);
  M_rotate.yx =  sin(theta);
  M_rotate.yy =  cos(theta);

  return(M_rotate);
}

PS::matrix set_rotation_matrix_along_x(double theta)
{
  PS::matrix M_rotate = PS::matrix(1, 0, 0,
	   0, 0, 0,
	   0, 0, 0 );
  M_rotate.yy =  cos(theta);
  M_rotate.yz = -sin(theta);
  M_rotate.zy =  sin(theta);
  M_rotate.zz =  cos(theta);

  return(M_rotate);
}

void calc_initial_position_and_velocity(double init_r1[],double init_v1[],double init_r2[],double init_v2[],double vimp,double R1,double M1,double R2,double M2,double theta_degree)
{
  double r1[D],r2[D];
  double v1[D],v2[D];
  double a1[D],a2[D];
  double delta,dr[D];

  double theta=theta_degree*(M_PI/180);
  double dt=-1.0;

  int d,n=0;

  r1[X]=R1*cos(theta);
  r1[Y]=R1*sin(theta);
  r1[Z]=0.0;
  r2[X]=-R2*cos(theta);
  r2[Y]=-R2*sin(theta);
  r2[Z]=0.0;
  v1[X]=-(M2/(M1+M2))*vimp;
  v1[Y]=0.0;
  v1[Z]=0.0;
  v2[X]=(M1/(M1+M2))*vimp;
  v2[Y]=0.0;
  v2[Z]=0.0;
  
  do{
    delta=0.0;
    for(d=0;d<D;d++){
      dr[d]=r1[d]-r2[d];
      delta+=dr[d]*dr[d];
    }
    delta=sqrt(delta);

    for(d=0;d<D;d++){
      a1[d] = -PARAM::G_CONST*M2*dr[d]/pow(delta,3.0);
      a2[d] = PARAM::G_CONST*M1*dr[d]/pow(delta,3.0);
    }

    for(d=0;d<D;d++){
      r1[d] += v1[d]*dt + 0.5*a1[d]*dt*dt;
      r2[d] += v2[d]*dt + 0.5*a2[d]*dt*dt;

      v1[d] += a1[d]*dt;
      v2[d] += a2[d]*dt;
    }

    //printf("%e %e %e %e %e %e\n",r1[X],r1[Y],r1[Z],r2[X],r2[Y],r2[Z]);
  }while(delta<2.0*(R1+R2));

  for(d=0;d<D;d++){
    init_r1[d]=r1[d];
    init_v1[d]=v1[d];
    init_r2[d]=r2[d];
    init_v2[d]=v2[d];
  }

}

int main(int argc,char* argv[])
{
  //input parameters////////////////////////////////////////////////////////
  char input_filename_target[]="FDPS-basalt-sphere-R=50km-N=100000-relaxed.bin"; //input filename for a target
  char input_filename_impactor[]="FDPS-basalt-sphere-R=50km-N=100000-relaxed.bin"; //input filename for an impactor
  double R1=50.0e5; //radius of a target [cm]
  double M1=(4.0/3.0)*pi*2.7*R1*R1*R1; //mass of a target [g]
  double R2=50.0e5; //radius of an impactor [cm]
  double M2=1.0*M1; //mass of an impactor [g]
  double time_init=0.0; //time for initial timestep [s]

  //in addition to these parameters, you should specify the output filename in the following for loops
  //////////////////////////////////////////////////////////////////////////
  
  int i,ix,iy,iz,i2,d,count;
  double x,y,z;
  double vesc=sqrt(2.0*PARAM::G_CONST*(M1+M2)/(R1+R2));
  double init_r1[D],init_v1[D],init_r2[D],init_v2[D];
  double rotate_theta;
  PS::matrix M_rotate_x = PS::matrix(1, 0, 0,
	     0, 1, 0,
	     0, 0, 1 );
  PS::matrix M_rotate_z = PS::matrix(1, 0, 0,
	     0, 1, 0,
	     0, 0, 1 );
  double realMtarget=0.0,realMimpactor=0.0;
  int realNtarget=0,realNimpactor=0;
  char filename[256];
  double vimp;
  int vimp_int;
  double theta_degree;
  int theta_degree_int;
  int flag=0;
  
  PS::Initialize(argc,argv);
  PS::ParticleSystem<RealPtcl> sph_system_target,sph_system_impactor,sph_system;
  sph_system_target.initialize();
  sph_system_impactor.initialize();
  sph_system.initialize();
  FileHeader header_target,header_impactor,header;
  
  sph_system_target.readParticleBinary(input_filename_target, header_target);
  sph_system_impactor.readParticleBinary(input_filename_impactor, header_impactor);
  header.time = time_init;
  header.Nbody = header_target.Nbody + header_impactor.Nbody;
  sph_system.setNumberOfParticleLocal(header.Nbody);

  printf("This program make initial condition file from \n");
  printf("%s and %s \n\n",input_filename_impactor,input_filename_target);

  //you can specify range and step of velocity and angle grids here/////
  for(theta_degree_int=10;theta_degree_int<=50;theta_degree_int+=10){
    for(vimp_int=2000;vimp_int<=15000;vimp_int+=1000){
	  //output filename//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      sprintf(filename,"FDPS-bsc-Rt=50km-q=1-theta=%d-v=%dms-Ntot=100000.bin",theta_degree_int,vimp_int);
	  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      vimp = vimp_int*100.0;
      theta_degree = theta_degree_int;
  
      calc_initial_position_and_velocity(init_r1,init_v1,init_r2,init_v2,vimp,R1,M1,R2,M2,theta_degree);
      
      //copy to real configuration for target
      for(i=0;i<header_target.Nbody;i++){
        sph_system[i]=sph_system_target[i];
        sph_system[i].pos.x += init_r1[X];
        sph_system[i].pos.y += init_r1[Y];
        sph_system[i].pos.z += init_r1[Z];
        sph_system[i].vel.x += init_v1[X];
        sph_system[i].vel.y += init_v1[Y];
        sph_system[i].vel.z += init_v1[Z];
        sph_system[i].id = i;
        if(flag==0){
          realMtarget+=sph_system[i].mass;
          realNtarget++;
        }
      }
	  
      //copy to real configutation for impactor
      for(i=0;i<header_impactor.Nbody;i++){
        i2=i+header_target.Nbody;
        sph_system[i2]=sph_system_impactor[i];
        sph_system[i2].pos.x += init_r2[X];
        sph_system[i2].pos.y += init_r2[Y];
        sph_system[i2].pos.z += init_r2[Z];
        sph_system[i2].vel.x += init_v2[X];
        sph_system[i2].vel.y += init_v2[Y];
        sph_system[i2].vel.z += init_v2[Z];
        sph_system[i2].id = i2;
        if(flag==0){
          realMimpactor+=sph_system[i].mass;
          realNimpactor++;
        }
      }
      flag=1;
      
      sph_system.writeParticleBinary(filename,header);
      printf("initial condition file: %s \n",filename);
      printf("condition: angle = %.2f[degree] speed = %.2e[cm/s] \n",theta_degree,vimp);
    }
  }

  printf("\n--------------------------------------------------------\n");
  printf("total partice number: %d \n",header.Nbody);
  printf("real target mass: %e, real impactor mass: %e\n",realMtarget,realMimpactor);
  printf("real target N: %d, real impactor N: %d\n",realNtarget,realNimpactor);
  printf("--------------------------------------------------------\n");
  
  PS::Finalize();
  return(0);
}
