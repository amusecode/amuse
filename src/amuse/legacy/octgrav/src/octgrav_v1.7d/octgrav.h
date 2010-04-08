#ifndef _OCTGRAV_H_
#define _OCTGRAV_H_

#include<iostream>
#include<vector>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>

using namespace std;

#include <sys/time.h>
#include <builtin_types.h>

struct double4 
{
  double x, y, z, w;
};

#include "n_per_cell.h"


#define INT_AS_FLOAT(x) (*((float*)&(x)))
inline int n_norm(int n, int j) {
  n = ((n-1)/j) * j + j;
  if (n == 0) n = j;
  return n;
}

#include "octgravdefs.h"
#include "oct_tree_base.h"
#include "tree_manip.h"

double get_time();

struct octgrav_data {
  
  /***************/
  /* bodies data */
  /***************/

  float4 *bodies_pos, *bodies_grav;

  /*************/
  /* tree data */
  /*************/

  int4 *children;  

  float4 *node_com, *node_pos, *node_Qu, *node_Qd;
  int    *n_in_node, *node_bodies_offset;
  float4 *Oct1, *Oct2;
  float2 *Oct3;
  
//   float4 *leaf_com, *leaf_pos;
//   int    *n_in_leaf, *leaf_bodies_offset;

  float4 *cell_com, *cell_pos;
  int    *n_in_cell, *cell_bodies_offset;
  
//   int *leaves, *leaf_pos_offset;

  /********************/
  /* interaction list */
  /********************/

  int *interaction_node, *interaction_node_offset, *interaction_node_len;
  int *interaction_leaf, *interaction_leaf_offset, *interaction_leaf_len;
};

extern "C"
{
  void initCUDA();
  void allocateCUDAarray(void** pos, int n);
  void deleteCUDAarray(void* pos);
  void copyArrayToDevice(void* device, const void* host, int n);
  void copyArrayFromDevice(void* host, const void* device, int n);
  double host_evaluate_gravity(float  inv_opening_angle,
			       float  softening_squared,
			       
			       int    n_bodies,
			       float4 *bodies_pos,
			       float4 *bodies_grav,

			       int    n_children,
			       int4   *children,

			       int    n_nodes,
			       float4 root_pos,
			       float4 root_com,
			       float4 *node_pos,
			       float4 *node_com,
			       float4 *node_Qu,
			       float4 *node_Qd,
			       float4 *Oct1,
			       float4 *Oct2,
			       float2 *Oct3,
			       int    *n_in_node,
			       int    *node_bodies_offset,
			       
			       int    n_cells,
			       float4 *cell_pos,
			       float4 *cell_com,
			       int    *n_in_cell,
			       int    *leaf_bodies_offset);
}

class octgrav {
protected:
  bool rebuild_tree, old_tree, first_run;

  oct_tree<N_PER_CELL> tree;
  octgrav_data hst, dev;
  
  int n_bodies;
  vector<float4> bodies_pos;
  vector<int>    bodies_id, bodies_id_orig;
  
  float softening, opening_angle;

  /*************************/

  int n_crit;
  vector<int4> children_list;
  vector<oct_tree<N_PER_CELL>*> node_list, leaf_list, cell_list;
//   vector<int> sorted_leaves;

  /*************************/

  float4 load_data(vector<float4>&);
  void reorder_bodies();

  /*************************/

  double CUDA_evaluate_gravity();

  /*************************/

  void allocate_bodies_list();
  void allocate_host_memory();
  void allocate_cuda_memory();
  void free_host_memory();
  void free_cuda_memory();

  /*************************/

  void prepare_data_for_device();
  void copy_data_to_device();
  void copy_data_from_device();
  
public:

  void evaluate_gravity(vector<float4>&, vector<float4>&);
  

  void  set_opening_angle(float theta) {opening_angle = theta;}
  float get_opening_angle() {return opening_angle;}

  void  set_softening(float eps) {softening = eps;}
  float get_softening() {return softening;}
  
  
  octgrav(float theta = 1.0,
	  float eps   = 0.25) {
    rebuild_tree = true;
    old_tree = false;
    opening_angle = theta;
    softening     = eps;
    first_run = true;

    cuda_interaction_list_len = 0;
    cuda_interaction_node_len = 0;
    cuda_interaction_leaf_len = 0;
    cuda_interaction_node_list = 0;
    cuda_interaction_leaf_list = 0;
    cuda_n_node = 0;
    cuda_n_leaf = 0;
    n_alloc = 0;
  }
	  
  ~octgrav() {

    if (cuda_interaction_list_len > 0) deleteCUDAarray((void*)dev_interaction_list_len);
    if (cuda_interaction_node_len > 0) deleteCUDAarray((void*)dev_interaction_node_len);
    if (cuda_interaction_leaf_len > 0) deleteCUDAarray((void*)dev_interaction_leaf_len);
    if (cuda_interaction_node_list > 0) deleteCUDAarray((void*)dev_interaction_node_list);
    if (cuda_interaction_leaf_list > 0) deleteCUDAarray((void*)dev_interaction_leaf_list);
    if (cuda_n_node > 0) deleteCUDAarray((void*)dev_n_node);
    if (cuda_n_leaf > 0) deleteCUDAarray((void*)dev_n_leaf);
    
  };
  
};

#endif 
