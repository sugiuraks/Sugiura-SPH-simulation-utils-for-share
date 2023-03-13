///////////////////////////////////////////////////////////////////////////////////////////////////////
//箱の粒子配置の初期条件ファイルの MAX_FLAW_NUMBER を変更するプログラムです
//input_filename(MAX_FLAW_NUMBER変更前)とoutput_filename(MAX_FLAW_NUMBER変更後)はmain関数の頭で指定しています
//箱の粒子配置であれば, input_filenameにはcgs単位系に変更後, 解像度変更後, サイズ変更後のものを使用しても大丈夫です

//使い方///////////////////////////////////
//my-src-for-flaw-number-change/ 以下にある param.h の頭に書いてある
//MAX_FLAW_NUMBER に変更前の MAX_FLAW_NUMBER を (私が渡した初期条件ファイルでは通常40にしてあるはずです)
//MAX_FLAW_NUMBER_FOR_OUTPUT に 変更後の MAX_FLAW_NUMBER を指定してください (1000万粒子なら50程度で行けるはずです)
//(Makefileはsrc, my-src-for-flaw-number-changeの場所やコンパイラなど適切に書き換えてください)
//このプログラムとMakefile, MAX_FLAW_NUMBER変更前の初期条件ファイルがあるディレクトリで
//$ make changeflawnumber
//$ ./change-flaw-number.out
//MAX_FLAW_NUMBERが変更されたあとの初期条件ファイル (output_filename) が生成されます
//////////////////////////////////////////

//output_filenameの初期条件ファイルは箱の粒子配置の初期条件ファイルなので
//いつも通り球に切り取ってターゲットとインパクタを作り, 初期速度と角度を与えて衝突条件を作ってください
//その際, 対応する param.h の頭に記入してある MAX_FLAW_NUMBER に変更後の値を入れてください
///////////////////////////////////////////////////////////////////////////////////////////////////////

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
  int i;
  //input parameters///////////////////////////////////////////////////////////////////////
  char input_filename[]="FDPS-basalt-box-R=50km-N=50000-after-relaxation.bin";
  char output_filename[]="FDPS-basalt-box-R=50km-N=50000-Nflaw=50-after-relaxation.bin";
  /////////////////////////////////////////////////////////////////////////////////////////
  
  PS::Initialize(argc,argv);
  PS::ParticleSystem<RealPtcl> sph_system_input,sph_system_output;
  sph_system_input.initialize();
  sph_system_output.initialize();
  FileHeader header_input,header_output;
  
  sph_system_input.readParticleBinary(input_filename, header_input);
  header_output = header_input;

  //copy
  sph_system_output.setNumberOfParticleLocal(header_output.Nbody);
  for(i=0;i<header_input.Nbody;i++){
    sph_system_output[i]=sph_system_input[i];
  }

  sph_system_output.writeParticleBinary(output_filename,header_output);
  printf("This program changes flaw number from %d to %d \n",MAX_FLAW_NUMBER,MAX_FLAW_NUMBER_OUTPUT);
  
  PS::Finalize();
  return(0);
}
