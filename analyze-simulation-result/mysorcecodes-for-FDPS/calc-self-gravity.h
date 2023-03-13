#ifndef INCLUDE_CALC_SELF_GRAVITY
#define INCLUDE_CALC_SELF_GRAVITY

//fanctions that need to be called in main fanction///////////////////////////////////////////////////////////////////////////////////

//about tree construction and deconstruction///////////////////////////

//make all tree nodes that have particle.
//return value...pointer to root node  argument...particle array
tree_node_t* make_tree(particle_t par[]);

//break all tree nodes.
//argument...pointer to tree node.If call it in main fanction,argument is pointer to root node.
void break_tree(tree_node_t *n_node);

//about calculation gravitational force or valiables//////////////////

//calculate tree node's variables.
//argument...pointer to tree node.If call it in main fanction,argument is pointer to root node.
void calc_tree_variables(tree_node_t *n_node);

//calculate all particle's self gravitational acceleration.
//argument 1...particle array  argument 2...pointer to root node.
void calc_self_grav_acceleration(particle_t par[],tree_node_t *root_node);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//fanctions that NOT need to be called in main fanction//////////////////////////////////////////////////////////////////////////////

//return X Y Z 's value correspond to index. index -> (X,Y,Z) is about
//0->(0,0,0) 1->(1,0,0) 2->(0,1,0) 3->(1,1,0) 4->(0,0,1) 5->(1,0,1) 6->(0,1,1) 7->(1,1,1)
//return value...'dir' value correspond to index  argument 1...index  argument 2...direction,each of X Y or Z
int convert_index_to_address(int index,int dir);

//return index from relation between tree node center and particle position
//return value...index  argument 1...a node center position array  argument 2...a particle position array
int convert_pos_to_index(double node_center[],double p_position[]);

//make tree node and initialize.if you want to make root node,please give NULL pointer
//return value...pointer to new node.if allocation failed,return NULL.  argument 1...pointer to new node' parent node.  argument 2...new node's index
tree_node_t* make_new_node(tree_node_t* parent,int child_index);

//insert particle to tree node and make list and so on
//argument 1...pointer to a particle.  argument 2...pointer to a tree node.
void insert_particle_to_node(particle_t* par_i,tree_node_t* n_node);

//calculate self gravitational force recursively using tree method. 
//argument 1...pointer to a particle.  argument 2...pointer to a tree node.
void calc_tree_grav(particle_t *par_i,tree_node_t *n_node);

//get root_node
tree_node_t* get_root_node(void);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
