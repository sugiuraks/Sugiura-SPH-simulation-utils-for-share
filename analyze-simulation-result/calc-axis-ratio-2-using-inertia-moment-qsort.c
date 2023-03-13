/////////////////////////////////////////////////////////////////////////////
//衝突計算結果の粒子情報が含まれるbinaryファイルを読み込んで,
//n番目に大きい天体を同定し, さらにその軸比などを計算するプログラムです.
//その後にその天体に重力的に束縛されているSPH粒子も同定してくれ(るはず)です.
//main関数一番上の nth_largest に解析したい天体が何番目に大きいものか指定してください.
//ただし, c言語の配列は0から始まるため, 最大天体は nth_largest=0 です.
//352行目からの部分でn番目に大きい天体を構成する粒子を辿るfor文を書いてありますので, ここを書き換えれば例えばn番目に大きい天体全体のエネルギーや運動量などを求めることができます.
//
//使い方///////////////////
//このプログラム と make-calc-axis-ratio-2-using-inertia-moment.sh と mysorcecodes-for-FDPS があるディレクトリで
//$ ./make-calc-axis-ratio-2-using-inertia-moment-qsort.sh
//$ ./calc-axis-ratio-2-using-inertia-moment-qsort 第二引数:解析したいbinaryファイルの名前
/////////////////////////////////////////////////////////////////////////////

#include "SPH.h"
#include "space.h"
#include "calc-pressure-force.h"
#include "calc-self-gravity.h"
#include "input.h"
#include "mirror.h"
#include "time-development.h"
#include "param-mysorcecode.h"
#include "tillotson.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

typedef struct clump_tag{
  int num_part;//number of particles that belong to this clump
  double Mass;//mass of all particles that belong to this clump
  double CoM[D];//center of mass of this clump
  int first_particle_ID;//remember the ID of first particle belonging to this clump for link list
  double mean_v[D];
} clump_t;

int compare_mass(const void *a, const void *b)
{
  if( ((clump_t *)a)->Mass > ((clump_t *)b)->Mass ) return(-1);
  else if( ((clump_t *)a)->Mass < ((clump_t *)b)->Mass ) return(1);
  else return(0);
}

//make clumps with close particles. Alone particles also make one clump.
clump_t* make_close_clumps(particle_t par[],int master_clump_num[],int* num_clumpP,int next_ID[])
{
  //important parameter///////////////////////////
  double nearby_criteria=1.5; //粒子間距離がこの値xスムージング長なら同じクランプに属するとする
  ////////////////////////////////////////////////
  
  int i,d,n;
  int particle_number=get_particle_number();
  clump_t* clumps=NULL;
  int smallest_num_of_clump;
  double delta,hij;
  particle_t* par_jP;
  int eaterID,killedID,lastID;
  
  for(i=0;i<particle_number;i++){
    //search the smallest number of clump among surrounding particles
    smallest_num_of_clump=particle_number;
    for(n=0;n<par[i].inter_N;n++){
      par_jP = par[i].inter_par[n];
      hij = 0.5*(par[i].h+par_jP->h);
      delta=0.0;
      for(d=0;d<D;d++){
        delta += (par[i].r[d]-par_jP->r[d])*(par[i].r[d]-par_jP->r[d]);
      }
      delta = sqrt(delta);
      if( delta<nearby_criteria*hij ){
        if(smallest_num_of_clump>master_clump_num[par_jP->ID]) smallest_num_of_clump=master_clump_num[par_jP->ID];
      }
    }
    
    //if there is no existing clumps nearby i-th particle
    if(smallest_num_of_clump==particle_number){
      //make new clump 
      clumps = (clump_t*)realloc(clumps,((*num_clumpP)+1)*sizeof(clump_t));
      clumps[(*num_clumpP)].num_part = 0;
      clumps[(*num_clumpP)].Mass = 0.0;
      clumps[(*num_clumpP)].first_particle_ID=-1;
      for(d=0;d<D;d++){
        clumps[(*num_clumpP)].CoM[d] = 0.0;
        clumps[(*num_clumpP)].mean_v[d] = 0.0;
      }
      eaterID = (*num_clumpP);  
     
      (*num_clumpP)++;
    }
    else{
      eaterID = smallest_num_of_clump;
    }

    //add surrounding particles to clump with the number of eaterID
    for(n=0;n<par[i].inter_N;n++){
      par_jP = par[i].inter_par[n];
      hij = 0.5*(par[i].h+par_jP->h);
      delta=0.0;
      for(d=0;d<D;d++){
        delta += (par[i].r[d]-par_jP->r[d])*(par[i].r[d]-par_jP->r[d]);
      }
      delta = sqrt(delta);
      if( delta<nearby_criteria*hij ){
        //if j-th particle is not belonging to any clump, add this particle.
        if(master_clump_num[par_jP->ID]==particle_number){
          master_clump_num[par_jP->ID] = eaterID;
          for(d=0;d<D;d++){
            clumps[eaterID].CoM[d] =  ( clumps[eaterID].CoM[d]*clumps[eaterID].Mass + par_jP->r[d]*par_jP->mass) / ( clumps[eaterID].Mass + par_jP->mass);
            clumps[eaterID].mean_v[d] = ( clumps[eaterID].mean_v[d]*clumps[eaterID].num_part + par_jP->v[d] ) / ( clumps[eaterID].num_part + 1 );
          }
          clumps[eaterID].num_part++;
          clumps[eaterID].Mass += par_jP->mass;
          next_ID[par_jP->ID]=clumps[eaterID].first_particle_ID;
          clumps[eaterID].first_particle_ID=par_jP->ID;
        }
        //if j-th particle is already belonging to any other clump, add all particles of that clump
        else if(master_clump_num[par_jP->ID]!=eaterID){
          killedID = master_clump_num[par_jP->ID];
          for(d=0;d<D;d++){
            clumps[eaterID].CoM[d] =  ( clumps[eaterID].CoM[d]*clumps[eaterID].Mass + clumps[killedID].CoM[d]*clumps[killedID].Mass) / ( clumps[eaterID].Mass + clumps[killedID].Mass);
            clumps[eaterID].mean_v[d] = ( clumps[eaterID].mean_v[d]*clumps[eaterID].num_part + clumps[killedID].mean_v[d]*clumps[killedID].num_part ) / ( clumps[eaterID].num_part + clumps[killedID].num_part );
          }
          clumps[eaterID].num_part += clumps[killedID].num_part;
          clumps[eaterID].Mass += clumps[killedID].Mass;
          //search last ID of clumps[killedID], and replace the master of all particles belonging to clumps[killedID] to eaterID
          lastID = clumps[killedID].first_particle_ID;
          while(next_ID[lastID]!=-1){
            master_clump_num[lastID] = eaterID;
            lastID = next_ID[lastID];
	  }
          master_clump_num[lastID] = eaterID;
          //add particles of clumps[killedID] to link list of clumps[eaterID]
          next_ID[lastID] = clumps[eaterID].first_particle_ID;
          clumps[eaterID].first_particle_ID = clumps[killedID].first_particle_ID; 
          
          //discard clump[killedID]
          clumps[killedID].num_part=0;
          clumps[killedID].Mass=0.0;
          for(d=0;d<D;d++){
            clumps[killedID].CoM[d]=0.0;
            clumps[killedID].mean_v[d]=0.0;
          }
          clumps[killedID].first_particle_ID=-1;
        }
      }
    }

  }
  return(clumps);
}

void make_gravitational_clump_iteratively(int largest_clump_ID,clump_t clumps[],int num_clump,int next_ID[])
{
  int n,d,i;
  int flag=0;
  double delta;
  double rel_v2;
  double mu;
  double Kinetic_E;
  double Grav_E;
  int last_ID;

  for(n=0;n<num_clump;n++){
    if(n!=largest_clump_ID && clumps[n].num_part!=0){
      delta=rel_v2=0.0;
      for(d=0;d<D;d++){
        delta += (clumps[largest_clump_ID].CoM[d]-clumps[n].CoM[d])*(clumps[largest_clump_ID].CoM[d]-clumps[n].CoM[d]);
        rel_v2 += (clumps[largest_clump_ID].mean_v[d]-clumps[n].mean_v[d])*(clumps[largest_clump_ID].mean_v[d]-clumps[n].mean_v[d]);
      }
      delta = sqrt(delta);
      mu = (clumps[largest_clump_ID].Mass*clumps[n].Mass) / (clumps[largest_clump_ID].Mass+clumps[n].Mass);

      Kinetic_E = 0.5 * mu * rel_v2;
      Grav_E = - g_const * clumps[largest_clump_ID].Mass * clumps[n].Mass / delta;

      //If clump n is gravitationally bounded by largest clump, clump n and largest clump marge.
      if(Kinetic_E+Grav_E<0.0){
        flag=1;
        for(d=0;d<D;d++){
          clumps[largest_clump_ID].CoM[d] =  ( clumps[largest_clump_ID].CoM[d]*clumps[largest_clump_ID].Mass + clumps[n].CoM[d]*clumps[n].Mass) / ( clumps[largest_clump_ID].Mass + clumps[n].Mass);
          clumps[largest_clump_ID].mean_v[d] = ( clumps[largest_clump_ID].mean_v[d]*clumps[largest_clump_ID].num_part + clumps[n].mean_v[d]*clumps[n].num_part ) / ( clumps[largest_clump_ID].num_part + clumps[n].num_part );
        }
        clumps[largest_clump_ID].num_part += clumps[n].num_part;
        clumps[largest_clump_ID].Mass += clumps[n].Mass;
        //search last ID of clumps[n]
        i = clumps[n].first_particle_ID;
        while(next_ID[i]!=-1){
          i = next_ID[i];
        }
        last_ID = i;
        //add particles of clumps[n] to link list of clumps[largest_clump_ID]
        next_ID[last_ID] = clumps[largest_clump_ID].first_particle_ID;
        clumps[largest_clump_ID].first_particle_ID = clumps[n].first_particle_ID;
        
        //discard clump n
        clumps[n].num_part=0;
        clumps[n].Mass=0.0;
        for(d=0;d<D;d++){
          clumps[n].CoM[d]=0.0;
          clumps[n].mean_v[d]=0.0;
        }
        clumps[n].first_particle_ID=-1;
      }
    }
  }
  
  //If one or more clumps marge to largest clump, the same fanction is called iteratively.
  if(flag==1){
    make_gravitational_clump_iteratively(largest_clump_ID,clumps,num_clump,next_ID);
  }

}

void calc_axis_ratio_using_inertia_moment(double** position,int particle_number,double Mt,double* b_aP,double* c_aP)
{
  double aAxis,bAxis,cAxis;
  int i,d,d2;
  double mass=Mt/particle_number;
  double r0[D];
  double Ixx=0.0,Ixy=0.0,Ixz=0.0,Iyx=0.0,Iyy=0.0,Iyz=0.0,Izx=0.0,Izy=0.0,Izz=0.0;
  double Ia,Ib,Ic;
  double I1,I2,I3;
  double km=1.0e5;
  double a,b,c,q,p;
  double _Complex alpha_plus,alpha_minus;
  double delta2;

  //calculate center of mass of clump
  r0[X]=r0[Y]=r0[Z]=0.0;
  for(i=0;i<particle_number;i++){
    r0[X] += position[i][X]/particle_number;
    r0[Y] += position[i][Y]/particle_number;
    r0[Z] += position[i][Z]/particle_number;
  }

  //calculate inertia moment tensor
  for(i=0;i<particle_number;i++){
    delta2 = (position[i][X]-r0[X])*(position[i][X]-r0[X]) + (position[i][Y]-r0[Y])*(position[i][Y]-r0[Y]) + (position[i][Z]-r0[Z])*(position[i][Z]-r0[Z]);
    Ixx += ( delta2 - (position[i][X]-r0[X])*(position[i][X]-r0[X]) )*mass;
    Ixy +=        - (position[i][X]-r0[X])*(position[i][Y]-r0[Y])*mass;
    Ixz +=        - (position[i][X]-r0[X])*(position[i][Z]-r0[Z])*mass;
    Iyx +=        - (position[i][Y]-r0[Y])*(position[i][X]-r0[X])*mass;
    Iyy += ( delta2 - (position[i][Y]-r0[Y])*(position[i][Y]-r0[Y]) )*mass;
    Iyz +=        - (position[i][Y]-r0[Y])*(position[i][Z]-r0[Z])*mass;
    Izx +=        - (position[i][Z]-r0[Z])*(position[i][X]-r0[X])*mass;
    Izy +=        - (position[i][Z]-r0[Z])*(position[i][Y]-r0[Y])*mass;
    Izz += ( delta2 - (position[i][Z]-r0[Z])*(position[i][Z]-r0[Z]) )*mass;
  }

  //calculate principal moment of inertia
  a = - (Ixx + Iyy + Izz);
  b = Ixx*Iyy + Iyy*Izz + Izz*Ixx - Ixy*Ixy - Iyz*Iyz - Izx*Izx;
  c = - (Ixx*Iyy*Izz + 2.0*Ixy*Iyz*Izx - Ixx*Iyz*Iyz - Iyy*Izx*Izx - Izz*Ixy*Ixy);
  q = ( 27.0*c + 2.0*a*a*a - 9.0*a*b ) / 54.0;
  p = ( 3.0*b - a*a ) / 9.0;
  alpha_plus = cpow( -q + csqrt( (double _Complex)(q*q + p*p*p) ) , 1.0/3.0 );
  alpha_minus = cpow( -q - csqrt( (double _Complex)(q*q + p*p*p) ) , 1.0/3.0 );
  Ia = creal( alpha_plus + alpha_minus - (a/3.0) );
  Ib = creal( 0.5 * ( -1.0 + I*sqrt(3.0) ) * alpha_plus + 0.5 * ( -1.0 - I*sqrt(3.0) ) * alpha_minus - (a/3.0) );
  Ic = creal( 0.5 * ( -1.0 - I*sqrt(3.0) ) * alpha_plus + 0.5 * ( -1.0 + I*sqrt(3.0) ) * alpha_minus - (a/3.0) );

  I1 = fmax( fmax(Ia,Ib), Ic);
  I3 = fmin( fmin(Ia,Ib), Ic);
  I2 = Ia + Ib + Ic - I1 - I3;

  //calculate axis length using ellipsoid approximation
  aAxis = sqrt( 10.0*( I1+I2-I3)/Mt );
  bAxis = sqrt( 10.0*( I1-I2+I3)/Mt );
  cAxis = sqrt( 10.0*(-I1+I2+I3)/Mt );

  *(b_aP) = bAxis/aAxis;
  *(c_aP) = cAxis/aAxis;

  printf("a=%ekm, b=%ekm, c=%ekm\n",aAxis/km,bAxis/km,cAxis/km);

}

int main(int argc,char *argv[])
{
  //define mass of target//////
  double Rt=50.0e5;
  double rho0=2.7;
  double Mt=(4.0/3.0)*M_PI*rho0*Rt*Rt*Rt;
  /////////////////////////////

  //nth largest clump that you want to analyze
  int nth_largest=0;
  ////////////////////////////////////////////
  
  particle_t* par = initialize(argv[1],argc);
  int i,n,l,d;
  int particle_number=get_particle_number();
  int num_clump=0; //number of all clumps
  clump_t* clumps; //pointer to clump allay
  int* master_clump_num=(int*)malloc(sizeof(int)*particle_number);//The number of clump this particle belongs to
  int* next_ID=(int*)malloc(sizeof(int)*particle_number);//next particle ID 
  char name_largest[]="largest.txt";
  FILE *fp_largest;
  for(i=0;i<particle_number;i++){
    master_clump_num[i]=particle_number;
    next_ID[i]=-1;
    //re-give the ID for each particle
    par[i].ID=i;
  }
  tree_node_t* root_node;
  double b_a,c_a;

  //identify clumps with friends-of-friends algorithm///////////////////////////////////////
  root_node = make_tree(par);
  neighbor_searching(par,root_node);
  clumps = make_close_clumps(par,master_clump_num,&num_clump,next_ID);
  //////////////////////////////////////////////////////////////////////////////////////////
  
  //sort clumps with respect to those mass//////////////////////////////////////////////////
  qsort(clumps,num_clump,sizeof(clump_t),compare_mass);
  //////////////////////////////////////////////////////////////////////////////////////////

  //print position of particles belonging to nth largest clump//////////////////////////////////
  fp_largest=fopen(name_largest,"w");
  for(i=clumps[nth_largest].first_particle_ID;i!=-1;i=next_ID[i]){
    fprintf(fp_largest,"%.3e %.3e %.3e %d\n",par[i].r[X],par[i].r[Y],par[i].r[Z],par[i].ID);
  }
  fclose(fp_largest);
  ///////////////////////////////////////////////////////////////////////////////////////////

  //calculate and print axis ratio of the largest remnant////////////////////////////////////
  int particle_number_nth_largest=clumps[nth_largest].num_part;
  double **position;
  if((position = (double**)malloc(sizeof(double*)*particle_number_nth_largest))==NULL){
    printf("cannot allocate position array! \n");
    exit(0);
  }
  for(i=0;i<particle_number_nth_largest;i++){
    if((position[i] = (double*)malloc(sizeof(double)*D))==NULL){
      printf("cannot allocate position array! \n");
      exit(0);
    }
  }
  n=0;
  for(i=clumps[nth_largest].first_particle_ID;i!=-1;i=next_ID[i]){
    position[n][X]=par[i].r[X];
    position[n][Y]=par[i].r[Y];
    position[n][Z]=par[i].r[Z];
    n++;
  }
  calc_axis_ratio_using_inertia_moment(position,particle_number_nth_largest,Mt,&b_a,&c_a);
  
  printf("%d th largest: b/a=%f c/a=%f mass=%e fMlr=%f \n",nth_largest,b_a,c_a,clumps[nth_largest].Mass,clumps[nth_largest].Mass/Mt);
  /////////////////////////////////////////////////////////////////////////////////////////////

  //テンプレート, n番目に大きい集積天体を構成するSPH粒子をたどる部分です////////////////////////////////////////////
  //ここでは例としてn番目に大きい集積天体の質量を計算してみます(clumps[nth_largest].Massにすでに保存されているが)/
  double MRemnant=0;
  for(i=clumps[nth_largest].first_particle_ID;i!=-1;i=next_ID[i]){
    MRemnant += par[i].mass;
  }
  printf("MRemnant = %.4e clumps[nth_largest].Mass = %.4e \n",MRemnant,clumps[nth_largest].Mass);
  /////////////////////////////////////////////////////////////////////////////////////////////////
  
  //calculate the mass fraction of particles that is gravitationally bounded by nth largest clump
  make_gravitational_clump_iteratively(nth_largest,clumps,num_clump,next_ID);

  //テンプレート, 最大集積天体に重力的に束縛されているSPH粒子(最大集積天体そのものを含む)をたどる部分です///////
  //ここでも例として最大集積天体に束縛されているSPH粒子の質量を計算してみます///////////////////////////////
  double MRemnantGrav=0;
  for(i=clumps[nth_largest].first_particle_ID;i!=-1;i=next_ID[i]){
    MRemnantGrav += par[i].mass;
  }
  printf("MRemnantGrav = %.4e clumps[nth_largest].Mass = %.4e \n",MRemnantGrav,clumps[nth_largest].Mass);
  /////////////////////////////////////////////////////////////////////////////////////////////////
  
  free_particle_array(par);
  break_tree(root_node);
  free(clumps);
  for(i=0;i<particle_number_nth_largest;i++){
    free(position[i]);
  }
  free(position);
  
  return(0);
}
