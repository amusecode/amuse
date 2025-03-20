#include "octgrav.h"

void octgrav::allocate_host_memory() {
//   hst.bodies_pos  = (float4*)malloc(n_bodies * sizeof(float4));
//   hst.bodies_grav = (float4*)malloc(n_bodies * sizeof(float4));

  hst.children = (int4*)malloc(children_list.size() * sizeof(int4));
  
  hst.node_Qu       = (float4*)malloc(node_list.size() * sizeof(float4));
  hst.node_Qd       = (float4*)malloc(node_list.size() * sizeof(float4));

  hst.Oct1       = (float4*)malloc(node_list.size() * sizeof(float4));
  hst.Oct2       = (float4*)malloc(node_list.size() * sizeof(float4));
  hst.Oct3       = (float2*)malloc(node_list.size() * sizeof(float2));

  hst.node_com      = (float4*)malloc(node_list.size() * sizeof(float4));
  hst.node_pos      = (float4*)malloc(node_list.size() * sizeof(float4));
  hst.n_in_node     = (int*   )malloc(node_list.size() * sizeof(int));
  hst.node_bodies_offset = (int*)malloc(node_list.size() * sizeof(int));

//   hst.leaf_com  = (float4*)malloc(leaf_list.size() * sizeof(float4));
//   hst.leaf_pos  = (float4*)malloc(leaf_list.size() * sizeof(float4));
//   hst.n_in_leaf = (int*   )malloc(leaf_list.size() * sizeof(int));
//   hst.leaf_bodies_offset = (int*)malloc(leaf_list.size() * sizeof(int));
  
  hst.cell_com  = (float4*)malloc(cell_list.size() * sizeof(float4));
  hst.cell_pos  = (float4*)malloc(cell_list.size() * sizeof(float4));
  hst.n_in_cell = (int*   )malloc(cell_list.size() * sizeof(int));
  hst.cell_bodies_offset = (int*)malloc(cell_list.size() * sizeof(int));
}

void octgrav::free_host_memory() {
  free(hst.bodies_pos);
  free(hst.bodies_grav);

  free(hst.children);

  free(hst.node_Qu);
  free(hst.node_Qd);

  free(hst.Oct1);
  free(hst.Oct2);
  free(hst.Oct3);

  free(hst.node_com);
  free(hst.node_pos);
  free(hst.n_in_node);
  free(hst.node_bodies_offset);
  
//   free(hst.leaf_com);
//   free(hst.leaf_pos);
//   free(hst.n_in_leaf);
//   free(hst.leaf_bodies_offset);

  free(hst.cell_com);
  free(hst.cell_pos);
  free(hst.n_in_cell);
  free(hst.cell_bodies_offset);
}

void octgrav::allocate_cuda_memory() {
  allocateCUDAarray((void**)&dev.bodies_pos,  n_bodies * sizeof(float4));
  allocateCUDAarray((void**)&dev.bodies_grav, n_norm(n_bodies, 256) * sizeof(float4));

  allocateCUDAarray((void**)&dev.children, children_list.size() * sizeof(int4));

  allocateCUDAarray((void**)&dev.node_Qu,      node_list.size() * sizeof(float4));
  allocateCUDAarray((void**)&dev.node_Qd,      node_list.size() * sizeof(float4));


  allocateCUDAarray((void**)&dev.Oct1,      node_list.size() * sizeof(float4));
  allocateCUDAarray((void**)&dev.Oct2,      node_list.size() * sizeof(float4));
  allocateCUDAarray((void**)&dev.Oct3,      node_list.size() * sizeof(float2));

  allocateCUDAarray((void**)&dev.node_pos,      node_list.size() * sizeof(float4));
  allocateCUDAarray((void**)&dev.node_com,      node_list.size() * sizeof(float4));
  allocateCUDAarray((void**)&dev.n_in_node,     node_list.size() * sizeof(int));
  allocateCUDAarray((void**)&dev.node_bodies_offset, node_list.size() * sizeof(int));

//   allocateCUDAarray((void**)&dev.leaf_pos,  leaf_list.size() * sizeof(float4));
//   allocateCUDAarray((void**)&dev.leaf_com,  leaf_list.size() * sizeof(float4));
//   allocateCUDAarray((void**)&dev.n_in_leaf, leaf_list.size() * sizeof(int));
//   allocateCUDAarray((void**)&dev.leaf_bodies_offset, leaf_list.size() * sizeof(int));

  allocateCUDAarray((void**)&dev.cell_pos,  cell_list.size() * sizeof(float4));
  allocateCUDAarray((void**)&dev.cell_com,  cell_list.size() * sizeof(float4));
  allocateCUDAarray((void**)&dev.n_in_cell, cell_list.size() * sizeof(int));
  allocateCUDAarray((void**)&dev.cell_bodies_offset, cell_list.size() * sizeof(int));
}

void octgrav::free_cuda_memory() {

   deleteCUDAarray((void*)dev.bodies_pos);
   deleteCUDAarray((void*)dev.bodies_grav);

   deleteCUDAarray((void*)dev.children);

   deleteCUDAarray((void*)dev.node_Qu);
   deleteCUDAarray((void*)dev.node_Qd);

   deleteCUDAarray((void*)dev.Oct1);
   deleteCUDAarray((void*)dev.Oct2);
   deleteCUDAarray((void*)dev.Oct3);

   deleteCUDAarray((void*)dev.node_pos);
   deleteCUDAarray((void*)dev.node_com);
   deleteCUDAarray((void*)dev.n_in_node);
   deleteCUDAarray((void*)dev.node_bodies_offset);

//    deleteCUDAarray((void*)dev.leaf_pos);
//    deleteCUDAarray((void*)dev.leaf_com);
//    deleteCUDAarray((void*)dev.n_in_leaf);
//    deleteCUDAarray((void*)dev.leaf_bodies_offset);

   deleteCUDAarray((void*)dev.cell_pos);
   deleteCUDAarray((void*)dev.cell_com);
   deleteCUDAarray((void*)dev.n_in_cell);
   deleteCUDAarray((void*)dev.cell_bodies_offset);

}
