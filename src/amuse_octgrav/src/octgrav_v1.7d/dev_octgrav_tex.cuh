#ifndef _DEV_OCTGRAV_TEX_CUH_
#define _DEV_OCTGRAV_TEX_CUH_

/* bodies array */

texture<float4, 1, cudaReadModeElementType> bodies_pos_tex;

/* for tree-walk */

texture<int4,   1, cudaReadModeElementType> children_tex;
texture<float4, 1, cudaReadModeElementType> cell_pos_tex;
texture<float4, 1, cudaReadModeElementType> cell_com_tex;
texture<int,    1, cudaReadModeElementType> n_in_cell_tex;

/* for evaluating interaction list */

texture<float4, 1, cudaReadModeElementType> node_pos_tex;
texture<float4, 1, cudaReadModeElementType> node_com_tex;
texture<float4, 1, cudaReadModeElementType> node_Qu_tex;
texture<float4, 1, cudaReadModeElementType> node_Qd_tex;

texture<float4, 1, cudaReadModeElementType> Oct1_tex;
texture<float4, 1, cudaReadModeElementType> Oct2_tex;
texture<float2, 1, cudaReadModeElementType> Oct3_tex;

texture<int,    1, cudaReadModeElementType> n_in_node_tex;
texture<int,    1, cudaReadModeElementType> node_bodies_offset_tex;
texture<int,    1, cudaReadModeElementType> cell_bodies_offset_tex;

texture<int,    1, cudaReadModeElementType> interaction_node_tex;
texture<int,    1, cudaReadModeElementType> interaction_leaf_tex;

__constant__ int    n_cells;
__constant__ int    n_bodies;
__constant__ int    n_nodes;
__constant__ float  inv_opening_angle;
__constant__ float  softening_squared;
__constant__ float4 root_pos;
__constant__ float4 root_com;

#endif
