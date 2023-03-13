/////////////////////////////////////////////////////////////
//FDPS-basalt-box-R=50km-N=50000-after-relaxation-original.binは,
//乱雑な粒子配置だが一様密度(rho=2.7g/cc)になるような粒子配置で,
//1辺が100kmの箱の中に粒子を詰めてあり,
//半径50kmの球状にこの粒子配置を切り取ることで粒子数が約5万になります.
//このプログラムは, この粒子配置を適切に拡大縮小することで,
//解像度の異なる箱を作るものです.
//main関数の頭で適宜情報を入力しますが, 重要なのは
//Ninput:  input_filenameの粒子配置を半径50kmに切り取った時の粒子数
//Noutput: output_filenameの粒子配置を半径50kmに切り取った時の粒子数
//
//使い方///////////////////////
//(Makefileはsrc, my-srcの場所やコンパイラなど適切に書き換えてください)
//このプログラムとMakefileがあるディレクトリで
//$ make makedifferentresolution
//$ ./make-different-resolution-from-relaxed-condition.out
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
  //input parameters//////////////////////////
  char input_filename[]="FDPS-basalt-box-R=50km-N=50000-after-relaxation-original-cgs.bin"; //input particle information file
  char output_filename[]="FDPS-basalt-box-R=50km-N=100000-after-relaxation.bin"; //output particle information file
  double Ninput=5.0e4; //resolution for a sphere of input file
  double Noutput=10.0e4; //desired resolution for a sphere of output file
  double time=0.0; //initial time
  double rho0=2.7; //density
  double boxSize=100.0e5; //size of the box of input_filename
  ///////////////////////////////////////////
  
  int i,ix,iy,iz,i2,d;
  double x,y,z;
  double lengthRatio=pow(Ninput/Noutput,1.0/3.0); //scaling for length
  int Ncell=(int)ceil(1.0/lengthRatio); //Number of cells of initial box to fulfill 100km box of output resolution
  double iCellX,iCellY,iCellZ;
  int counterX,counterY,counterZ;
  double initI = 0.5 - ((double)Ncell/2.0);
  double boxLength = lengthRatio*boxSize;
  
  PS::Initialize(argc,argv);
  PS::ParticleSystem<RealPtcl> sph_system_input,sph_system_output;
  sph_system_input.initialize();
  sph_system_output.initialize();
  FileHeader header_input,header_output;
  
  sph_system_input.readParticleBinary(input_filename, header_input);
  header_output.time = header_input.time;
  
  //calculate particle number for output file
  header_output.Nbody=header_input.Nbody*Ncell*Ncell*Ncell;
  sph_system_output.setNumberOfParticleLocal(header_output.Nbody);
  
  i2=0;
  //copy to real configuration
  for(iCellX=initI,counterX=0 ; counterX<Ncell ; iCellX+=1.0,counterX++){
    for(iCellY=initI,counterY=0 ; counterY<Ncell ; iCellY+=1.0,counterY++){
      for(iCellZ=initI,counterZ=0 ; counterZ<Ncell ; iCellZ+=1.0,counterZ++){
        for(i=0;i<header_input.Nbody;i++){
          sph_system_output[i2]=sph_system_input[i];
          sph_system_output[i2].pos.x = sph_system_input[i].pos.x*lengthRatio + iCellX*boxLength;
          sph_system_output[i2].pos.y = sph_system_input[i].pos.y*lengthRatio + iCellY*boxLength;
          sph_system_output[i2].pos.z = sph_system_input[i].pos.z*lengthRatio + iCellZ*boxLength;
          sph_system_output[i2].mass  = sph_system_input[i].mass*lengthRatio*lengthRatio*lengthRatio;
          sph_system_output[i2].smth  = sph_system_input[i].smth*lengthRatio;
          sph_system_output[i2].id = i2;
          sph_system_output[i2].dens = rho0;
          i2++;
        }
      }
    }
  }
  
  double km=1.0e5;
  FILE* fpout;
  char filenameout[]="test.txt";
  fpout=fopen(filenameout,"w");
  fprintf(fpout,"%f \n",0.0);
  fprintf(fpout,"%d \n",header_output.Nbody);
  for(i=0;i<header_output.Nbody;i++){
    fprintf(fpout,"%.2f %.2f %.2f %.2f %.2e %f %d\n",sph_system_output[i].pos.x/km,sph_system_output[i].pos.y/km,sph_system_output[i].pos.z/km,sph_system_output[i].smth/km,sph_system_output[i].mass,sph_system_output[i].dens,sph_system_output[i].id);
  }
  fclose(fpout);
  
  sph_system_output.writeParticleBinary(output_filename,header_output);
  printf("This program make different resolution : %s from \n",output_filename);
  printf("%s \n",input_filename);
  printf("N for original sphere: %d, for output file: %d \n",(int)Ninput,(int)Noutput);
  
  PS::Finalize();
  return(0);
}
