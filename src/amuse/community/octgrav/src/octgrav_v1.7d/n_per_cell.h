#ifndef N_PER_CELL

#define N_PER_CELL 8
#define NCRIT   128
#define NTHREADS 256
#define NBLOCKS  128

extern "C" {
  extern int  cuda_interaction_list_len;
  extern int  cuda_interaction_node_len;
  extern int  cuda_interaction_leaf_len;
  extern int  cuda_interaction_node_list;
  extern int  cuda_interaction_leaf_list;
  extern int  cuda_n_node;
  extern int  cuda_n_leaf;
  extern int  n_alloc;

  extern int3 *dev_interaction_list_len;
  extern int2 *dev_interaction_node_len;
  extern int2 *dev_interaction_leaf_len;
  extern int  *dev_interaction_node_list;
  extern int  *dev_interaction_leaf_list;
  extern int  *dev_n_node;
  extern int  *dev_n_leaf;
}

#endif
