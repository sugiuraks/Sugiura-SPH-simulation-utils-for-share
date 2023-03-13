#include "SPH.h"
#include "space.h"
#include "input.h"
#include "mirror.h"
#include "param-mysorcecode.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static tree_node_t* root_node_s;

//fanctions for tree construction and deconstuction///////////////////////////////////////////////////////////////////////////

//return X Y Z 's value correspond to index. index -> (X,Y,Z) is about
//0->(0,0,0) 1->(1,0,0) 2->(0,1,0) 3->(1,1,0) 4->(0,0,1) 5->(1,0,1) 6->(0,1,1) 7->(1,1,1)
int convert_index_to_address(int index,int dir)
{
  if(dir==X) return(index%2);
  if(dir==Y) return((index/2)%2);
  if(dir==Z) return(index/4);
  
  return(0);
}

//return index from relation between tree node center and particle position
int convert_pos_to_index(double node_center[],double p_position[])
{
  int i[3],d;
  for(d=0;d<D;d++){
    if(node_center[d]>p_position[d]) i[d]=0;
    else i[d]=1;
  }
  if(D==1) i[Y]=i[Z]=0;
  if(D==2) i[Z]=0;
  
  return(4*i[Z]+2*i[Y]+i[X]);
}

//make tree node and initialize.if you want to make root node,please give NULL pointer
tree_node_t* make_new_node(tree_node_t* parent,int child_index)
{
  tree_node_t* new_node;
  int d,n;
  
  //allocation of new tree node
  new_node=(tree_node_t*)malloc(sizeof(tree_node_t));
  count_memory((int)sizeof(tree_node_t));
  if(new_node==NULL){
    return NULL;
  }
  
  //initialization of new tree node
  if(parent==NULL){
    //if new node is root node
    double Range[D];
    get_Range(Range);
    double Rmin[D];
    get_Rmin(Rmin);
    
    for(d=0;d<D;d++){
      new_node->box_size[d]=Range[d];
      new_node->center[d]=Rmin[d]+0.5*Range[d];
    }
    new_node->parent=NULL;
    new_node->level=0;
    new_node->C[X]=new_node->C[Y]=new_node->C[Z]=0;
  }
  else{
    //if new node is not root node
    for(d=0;d<D;d++){
      new_node->box_size[d] = 0.5 * parent->box_size[d];
      new_node->center[d]=parent->center[d]+(convert_index_to_address(child_index,d)-0.5)*new_node->box_size[d];
    }
    new_node->parent=parent;
    new_node->level=parent->level+1;
    for(d=0;d<3;d++){
      new_node->C[d]=2*parent->C[d]+convert_index_to_address(child_index,d);
    }
  }
  //for all new node
  for(d=0;d<D;d++){
    new_node->CoM[d]=0.0;
  }
  new_node->mass=0.0;
  new_node->hmax=0.0;
  for(n=0;n<N_CHILD;n++){
    new_node->child[n]=NULL;
  }
  new_node->first_particle=NULL;
  new_node->is_leaf=1;
  new_node->num_par=0;
  
  //return new node pointer
  return(new_node);
}

//insert particle to tree node and make list and so on
void insert_particle_to_node(particle_t* par_i,tree_node_t* n_node)
{
  particle_t* p_temp;
  particle_t* p_list;
  particle_t* p_next_temp;
  int n_index;
  
  n_node->num_par++;
  //if n_node is not full or level of n_node is max
  if( n_node->num_par<NODE_PAR_MAX || n_node->level==MAX_TREE_LEVEL){
    //insert particle to list
    p_temp=n_node->first_particle;
    n_node->first_particle=par_i;
    par_i->tree_next=p_temp;
    
    par_i->tree_parent=n_node;
  }
  else{
    //if n_node becomes full first
    if(n_node->is_leaf){
      //distribute list particle to child node and allocate child
      for(p_list=n_node->first_particle;p_list!=NULL;p_list=p_next_temp){
	p_next_temp=p_list->tree_next;
	n_index=convert_pos_to_index(n_node->center,p_list->r);
	//if n_node's child is not allocated
	if(n_node->child[n_index]==NULL){
	  n_node->child[n_index]=make_new_node(n_node,n_index);
	  //allocation failure
	  if(n_node->child[n_index]==NULL){
	    printf("cannot allocate child node! \n");
	    exit(FAILURE);
	  }
	}
	//call same fanction recursively to list particle
	insert_particle_to_node(p_list,n_node->child[n_index]);
      }
      n_node->is_leaf=0;
    }
    n_index=convert_pos_to_index(n_node->center,par_i->r);
    //if n_node's child is not allocated
    if(n_node->child[n_index]==NULL){
      n_node->child[n_index]=make_new_node(n_node,n_index);
      //allocation failure
      if(n_node->child[n_index]==NULL){
	printf("cannot allocate child node! \n");
	exit(FAILURE);
      }
    }
    //call same fanction recursively
    insert_particle_to_node(par_i,n_node->child[n_index]);
  }
}

//make tree from particle info and return pointer to root node
tree_node_t* make_tree(particle_t par[])
{
  int i,n,d;
  int particle_number=get_particle_number();
  if(get_is_mirror()){
    particle_number += get_mirror_particle_number();
  }
  tree_node_t* root_node;
  
  //make root node
  root_node=make_new_node(NULL,0);
  //allocation failure
  if(root_node==NULL){
    printf("cannot allocate root node! \n");
    exit(FAILURE);
  }
  
  //particle roop
  for(i=0;i<particle_number;i++){
    insert_particle_to_node(&par[i],root_node);
  }
  
  //check node level is really smaller than MAX_TREE_LEVEL
  for(i=0;i<particle_number;i++){
    if( par[i].tree_parent->level > MAX_TREE_LEVEL ){
      printf("tree level is too large! \n");
      exit(FAILURE);
    }
  }
  
  root_node_s=root_node;
  return(root_node);
}

//break tree nodes
void break_tree(tree_node_t *n_node)
{
  int n;
  
  //free all child node
  for(n=0;n<N_CHILD;n++){
    if(n_node->child[n]!=NULL){
      break_tree(n_node->child[n]);
    }
  }
  //free itself
  free(n_node);
  count_memory(-(int)sizeof(tree_node_t));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//fanctions for calculating tree variables and gravitational force///////////////////////////////////////////////////////

//calclate all tree variables
void calc_tree_variables(tree_node_t *n_node)
{
  particle_t *p_list;
  int d,n;
  
  //if n_node is leaf node,calculate variables from all particles this node have.
  if(n_node->is_leaf){
    for(p_list=n_node->first_particle ; p_list ; p_list=p_list->tree_next){
      n_node->mass += p_list->mass;
      for(d=0 ; d<D ; d++){
	n_node->CoM[d] += p_list->mass * p_list->r[d];
      }
      if(n_node->hmax < p_list->h) n_node->hmax = p_list->h;
    }
    for(d=0 ; d<D ; d++){
      n_node->CoM[d] /= n_node->mass;
    }
  }
  //if n_node is not reaf node,let all n_node's child calculate variables and calculate n_node's variables.
  else{
    for(n=0 ; n<N_CHILD ; n++){
      if(n_node->child[n]){
	calc_tree_variables(n_node->child[n]);
      }
    }
    for(n=0 ; n<N_CHILD ; n++){
      if(n_node->child[n]){
	n_node->mass += (n_node->child[n])->mass;
	for(d=0 ; d<D ; d++){
	  n_node->CoM[d] += ((n_node->child[n])->CoM[d]) * ((n_node->child[n])->mass);
	}
	if(n_node->hmax < (n_node->child[n])->hmax) n_node->hmax = (n_node->child[n])->hmax;
      }
    }
    for(d=0 ; d<D ; d++){
      n_node->CoM[d] /= n_node->mass;
    }
  }
}

//calculate self gravitational force recursively using tree method. 
void calc_tree_grav(particle_t *par_i,tree_node_t *n_node)
{
  double delta2=0.0;
  double dif[3];
  double temp;
  double Xtilda,hmaxi,x;
  double box_size=0.0;
  int d,n;
  particle_t *p_list;
  for(d=0 ; d<D ; d++){
    dif[d] = n_node->CoM[d] - par_i->r[d];
    delta2 += dif[d] * dif[d];
    if( box_size < n_node->box_size[d] ) box_size = n_node->box_size[d];
  }
  
  //if difference between n_node and i particle is large
  if((box_size/sqrt(delta2)) < theta){
    temp = g_const * n_node->mass * pow(delta2,-1.5);
    for(d=0 ; d<D ; d++){
      par_i->g[d] += dif[d] * temp;
    }
    par_i->grav_pot += -g_const * par_i->mass * n_node->mass * pow(delta2,-0.5);
  }
  //if difference is small and n_node is not leaf node
  else if(!n_node->is_leaf){
    for(n=0 ; n<N_CHILD ; n++){
      //if n_node->child[n] is not NULL,call same fanction recursively.
      if(n_node->child[n]){
	calc_tree_grav(par_i,n_node->child[n]);
      }
    }
  }
  //if n_node is leaf node
  else{
    for(p_list = n_node->first_particle ; p_list ; p_list = p_list->tree_next){
      delta2 = 0.0;
      for(d=0 ; d<D ; d++){
	dif[d] = p_list->r[d] - par_i->r[d];
	delta2 += dif[d] * dif[d];
      }
      //if 2 particles are far away
      if( sqrt(delta2) > ((par_i->h>p_list->h)?par_i->h:p_list->h) ){
	temp = g_const * p_list->mass * pow(delta2,-1.5);
	for(d=0 ; d<D ; d++){
	  par_i->g[d] += dif[d] * temp;
	}
      }
      //if 2 particles are close,calculate gravitational force considering particle's expanse.
      else{
	//should not calculate gravitational force from same particle
	if(par_i!=p_list){
	  //calculate 1.0 / hmax,hmax is larger of par_i->h or p_list->h
	  hmaxi=1.0 / ((par_i->h > p_list->h) ? par_i->h : p_list->h);
	  //x is ratio of hmin and hmax
	  x = ((par_i->h < p_list->h) ? par_i->h : p_list->h) * hmaxi;
	  //calculate fitting formula 
	  Xtilda = sqrt(delta2) * hmaxi * (1.0 / (1.0 + x * (0.007697 + x * (0.53799 - x * 0.125256))));
	  temp = g_const * p_list->mass * pow(delta2,-1.5) *( erf(Xtilda) - (2.0 / sqrt(PI)) * Xtilda * exp(- Xtilda * Xtilda));
	  for(d=0 ; d<D ; d++){
	    par_i->g[d] += dif[d] * temp;
	  }
	}
      }
      //when we calculate gravitational potential, we simply assume point particle
      if(par_i!=p_list){
	par_i->grav_pot += -g_const * par_i->mass * p_list->mass * pow(delta2,-0.5);
      }
    }
  }
}

//calculate all particle's self gravitational acceleration
void calc_self_grav_acceleration(particle_t par[],tree_node_t *root_node)
{
  int i,d;
  int particle_number=get_particle_number();
  
#ifdef _OPENMP
#pragma omp parallel for private(i,d)
#endif
  for(i=0;i<particle_number;i++){
    //initialize gravilational force and potential
    for(d=0;d<D;d++){
      par[i].g[d]=0.0;
    }
    par[i].grav_pot=0.0;
    //calculate particle i's gravitational force recursively using tree method
    calc_tree_grav(&par[i],root_node);
  }
}

//get root_node
tree_node_t* get_root_node(void)
{
  return(root_node_s);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
