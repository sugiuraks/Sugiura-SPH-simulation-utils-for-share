/////////////////////////////////////////////////////////////
//このプログラムは, input_filenameの粒子配置を
//lengthRatioの長さだけ適切に拡大縮小することで,
//サイズの違う初期条件の粒子配置の箱を準備するプログラムです。
//
//使い方///////////////////////
//(Makefileはsrc, my-srcの場所やコンパイラなど適切に書き換えてください)
//このプログラムとMakefileがあるディレクトリで
//$ make makedifferentscale
//$ ./make-different-scale-from-relaxed-condition.out
//////////////////////////////////////////////////////////////

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

int main(int argc,char* argv[])
{
  //input parameters////////////////////////////////////////////////////////////////
  char input_filename[]="FDPS-basalt-box-R=50km-N=50000-after-relaxation-original-mks.bin"; //input particle information file
  char output_filename[]="FDPS-basalt-box-R=25km-N=50000-after-relaxation.bin"; //output particle information file
  double time=0.0; //initial time
  double rho0=2.7; //density
  double lengthRatio=0.5; //scaling for length
  //////////////////////////////////////////////////////////////////////////////////

  int i,ix,iy,iz,i2,d;
  double x,y,z;
  double massRatio=pow(lengthRatio,3.0);
  
  PS::Initialize(argc,argv);
  PS::ParticleSystem<RealPtcl> sph_system_input,sph_system_output;
  sph_system_input.initialize();
  sph_system_output.initialize();
  FileHeader header_input,header_output;
  
  sph_system_input.readParticleBinary(input_filename, header_input);
  header_output.time = header_input.time;
  
  //calculate particle number for output file
  header_output.Nbody=header_input.Nbody;
  sph_system_output.setNumberOfParticleLocal(header_output.Nbody);
  
  i2=0;
  //copy to real configuration
  for(i=0;i<header_output.Nbody;i++){
    sph_system_output[i2]=sph_system_input[i];
    sph_system_output[i2].pos.x = sph_system_input[i].pos.x*lengthRatio;
    sph_system_output[i2].pos.y = sph_system_input[i].pos.y*lengthRatio;
    sph_system_output[i2].pos.z = sph_system_input[i].pos.z*lengthRatio;
    sph_system_output[i2].mass  = sph_system_input[i].mass*massRatio;
    sph_system_output[i2].smth  = sph_system_input[i].smth*lengthRatio;
    sph_system_output[i2].id = i2;
    sph_system_output[i2].dens = rho0;
    i2++;
  }
  
  double km=1.0e5;
  FILE* fpout;
  char filenameout[]="FDPS-basalt-sphere-collision-init.txt";
  fpout=fopen(filenameout,"w");
  //fprintf(fpout,"%f \n",0.0);
  //fprintf(fpout,"%d \n",header_output.Nbody);
  for(i=0;i<header_output.Nbody;i++){
    fprintf(fpout,"%.2e %.2e %.2e %.2e %.2f %f %d %.2e\n",sph_system_output[i].pos.x,sph_system_output[i].pos.y,sph_system_output[i].pos.z,sph_system_output[i].smth,sph_system_output[i].damage,sph_system_output[i].dens,sph_system_output[i].id,sph_system_output[i].mass);
    if(sph_system_output[i].property_tag!=0) printf("i=%d p=%d \n",i,sph_system_output[i].property_tag);
  }
  fclose(fpout);
  
  sph_system_output.writeParticleBinary(output_filename,header_output);
  printf("This program make different scale : %s from \n",output_filename);
  printf("%s \n",input_filename);
  printf("lengthRatio : %f \n",lengthRatio);
  
  PS::Finalize();
  return(0);
}
