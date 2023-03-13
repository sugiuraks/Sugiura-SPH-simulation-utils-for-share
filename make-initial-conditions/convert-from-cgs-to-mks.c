////////////////////////////////////////////////////////
//cgs単位系の初期条件ファイルである input_filename を
//mks単位系の初期条件ファイルである output_filename に変換するプログラムです
//
//使い方///////////////////////
//(Makefileはsrc, my-srcの場所やコンパイラなど適切に書き換えてください)
//このプログラムとMakefileがあるディレクトリで
//$ make convertfromcgstomks
//$ ./convert-from-cgs-to-mks.out
////////////////////////////////////////////////////////

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
  double time=0.0;
  int i,ix,iy,iz,i2,d;
  double x,y,z;
  char input_filename[]="FDPS-basalt-box-R=50km-N=50000-after-relaxation-original-cgs.bin";
  char output_filename[]="FDPS-basalt-box-R=50km-N=50000-after-relaxation-original-mks.bin";
  
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
    sph_system_output[i2].pos.x = sph_system_input[i].pos.x*0.01;
    sph_system_output[i2].pos.y = sph_system_input[i].pos.y*0.01;
    sph_system_output[i2].pos.z = sph_system_input[i].pos.z*0.01;
    sph_system_output[i2].vel.x = sph_system_input[i].vel.x*0.01;
    sph_system_output[i2].vel.y = sph_system_input[i].vel.y*0.01;
    sph_system_output[i2].vel.z = sph_system_input[i].vel.z*0.01;
    sph_system_output[i2].mass  = sph_system_input[i].mass*0.001;
    sph_system_output[i2].eng   = sph_system_input[i].eng*1.0e-4;
    sph_system_output[i2].smth  = sph_system_input[i].smth*0.01;
    sph_system_output[i2].dens  = sph_system_input[i].dens*1000.0;
    sph_system_output[i2].Sab_rho = sph_system_input[i].Sab_rho*1.0e-4;
    
    sph_system_output[i2].id = i2;
    i2++;
  }
  
  sph_system_output.writeParticleBinary(output_filename,header_output);
  printf("This program converts cgs file: %s \n",input_filename);
  printf("to mks file: %s \n",output_filename);
  
  PS::Finalize();
  return(0);
}
