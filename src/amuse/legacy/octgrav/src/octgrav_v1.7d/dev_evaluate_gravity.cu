#ifndef _DEV_EVALUATE_GRAVITY_CU
#define _DEV_EVALUATE_GRAVITY_CU

#ifndef __DEVICE_EMULATION__
#define _DEVICE_CODE_
#endif

#include "octgravdefs.h"
#include "dev_octgrav_tex.cuh"

#ifndef LMEM_STACK_SIZE
#define LMEM_STACK_SIZE 256
#endif

#define LEAF_BIT (1 << (24))

__device__ bool open_node(float4 cell_com,
			  float4 cell_pos,
			  float4 node_pos,
			  float4 node_com) {
  
  float3 dr = {fabs(node_com.x - cell_pos.x) - cell_pos.w,
	       fabs(node_com.y - cell_pos.y) - cell_pos.w,
	       fabs(node_com.z - cell_pos.z) - cell_pos.w};
  dr.x += fabs(dr.x); dr.x *= 0.5f;
  dr.y += fabs(dr.y); dr.y *= 0.5f;
  dr.z += fabs(dr.z); dr.z *= 0.5f;

  float ds = sqrtf(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
  return ( 2.0f*node_pos.w*inv_opening_angle > ds - cell_com.w);
  

}

/**********************************************
 *   compute length of the interaction list   *
 **********************************************/

template<int octant>
__device__ int4 interact_len(int node, int node_old,
			     float4 cell_com,
			     float4 cell_pos,
			     int    *ids_stack,
			     int4 stack) {
  /* if empty, exit */
  if (node == 0) return stack;
  
  /* check if the leaf or node has to be opened */
  float4 node_pos = tex1Dfetch(node_pos_tex, (node_old << (2)) + octant);
  float4 node_com = tex1Dfetch(node_com_tex, (node_old << (2)) + octant);
  stack.w += 8;

  if (open_node(cell_com, cell_pos, node_pos, node_com)) {
    if ((node & LEAF_BIT) == 0) {                        /* if node, */
      ids_stack[stack.x] = node;                         /* store it in stack */
      stack.x++;
      stack.w += 1;
    } else {
      stack.z++;                            /* otherwise account for this leaf */
    }
  } else {
    stack.y++;                              /* account for the node */
  }

  return stack;
}

__device__ int3 compute_interaction_list_len(float4 cell_com,
					     float4 cell_pos) {
  int     ids_stack[LMEM_STACK_SIZE];
  
  int    node     = 0;

  int4 stack = {0,0,0,0};
  ids_stack[stack.x] = node;
  stack.w += 1;
  
  stack.x++;
  
  while(stack.x > 0) {
    
    /* read node id & pos */
    stack.x--;
    node     = ids_stack[stack.x];
    stack.w += 1;    /* 1 for id & 4 for pos */
    
    int4 up = tex1Dfetch(children_tex, node + 0);
    int4 dn = tex1Dfetch(children_tex, node + 1);
    stack.w += 8;

#define INTERACT_LEN(oct, child) \
    {stack = interact_len<oct>(child, \
                               node, \
                               cell_com, \
			       cell_pos, \
			       ids_stack, \
			       stack);}
    INTERACT_LEN(0, up.x);
    INTERACT_LEN(1, up.y);
    INTERACT_LEN(2, up.z);
    INTERACT_LEN(3, up.w);
    INTERACT_LEN(4, dn.x);
    INTERACT_LEN(5, dn.y);
    INTERACT_LEN(6, dn.z);
    INTERACT_LEN(7, dn.w);
  }
  
  /*
   *  number of nodes,
   *  number of leaves,
   *  number of reads from + writes to memory.
   */
  
  return make_int3(stack.y, stack.z, stack.w);
}

__global__ void dev_compute_interaction_list_len(int3  *interaction_list_len) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (index >= n_cells)
    index = threadIdx.x;

  float4 cell_com = tex1Dfetch(cell_com_tex, index);
  float4 cell_pos = tex1Dfetch(cell_pos_tex, index);
  cell_com.w = sqrtf((cell_com.x - cell_pos.x)*(cell_com.x - cell_pos.x)+
		     (cell_com.y - cell_pos.y)*(cell_com.y - cell_pos.y)+
		     (cell_com.z - cell_pos.z)*(cell_com.z - cell_pos.z));
  
  interaction_list_len[index] = 
    compute_interaction_list_len(cell_com, cell_pos);
}

/****************************
 *  build interaction list  *
 ****************************/
template<int octant>
__device__ int3 interact_bld(int node, int node_old,
			     float4 cell_com,
			     float4 cell_pos,
			     int    *ids_stack,
			     int *interaction_node_list,
			     int *interaction_leaf_list,
			     int3 stack) {
  if (node == 0) return stack;
  
  float4 node_pos = tex1Dfetch(node_pos_tex, (node_old << (2)) + octant);
  float4 node_com = tex1Dfetch(node_com_tex, (node_old << (2)) + octant);

  if (open_node(cell_com, cell_pos, node_pos, node_com)) {
    if ((node & LEAF_BIT) == 0) {           /* if node, */
      ids_stack[stack.x] = node;            /* store it in stack */
      stack.x++;
    } else {				       
      interaction_leaf_list[stack.z++] = (node_old << (2)) + octant;
    }
  } else {
    interaction_node_list[stack.y++] = (node_old << (2)) + octant;
  }
  
  return stack;
}

__device__ void build_interaction_list(float4 cell_com,
				       float4 cell_pos,
				       int *interaction_node_list,
				       int *interaction_leaf_list) {
  int    ids_stack[LMEM_STACK_SIZE];
  
  int    node     = 0;
  int3 stack = {0, 0, 0};
  ids_stack[stack.x] = node;
  stack.x++;
  
  while(stack.x > 0) {
    
    /* read node id */
    stack.x--;
    node     = ids_stack[stack.x];
    
    int4 up = tex1Dfetch(children_tex, node + 0);
    int4 dn = tex1Dfetch(children_tex, node + 1);

#define INTERACT_BUILD(oct, child) \
    {stack = interact_bld<oct>(child, \
			       node,  \
                               cell_com, \
			       cell_pos, \
			       ids_stack, \
			       interaction_node_list, \
			       interaction_leaf_list, \
			       stack);}
    INTERACT_BUILD(0, up.x);
    INTERACT_BUILD(1, up.y);
    INTERACT_BUILD(2, up.z);
    INTERACT_BUILD(3, up.w);
    INTERACT_BUILD(4, dn.x);
    INTERACT_BUILD(5, dn.y);
    INTERACT_BUILD(6, dn.z);
    INTERACT_BUILD(7, dn.w);
  }
  
}

__global__ void dev_build_interaction_list(int cell_offset,
					   int *interaction_node_list,
					   int2 *interaction_node_offset,
					   int *interaction_leaf_list,
					   int2 *interaction_leaf_offset) {
  int index = cell_offset + (blockIdx.x * blockDim.x + threadIdx.x);
  

  if (index < n_cells) {
    float4 cell_com = tex1Dfetch(cell_com_tex, index);
    float4 cell_pos = tex1Dfetch(cell_pos_tex, index);
    cell_com.w = sqrtf((cell_com.x - cell_pos.x)*(cell_com.x - cell_pos.x)+
		       (cell_com.y - cell_pos.y)*(cell_com.y - cell_pos.y)+
		       (cell_com.z - cell_pos.z)*(cell_com.z - cell_pos.z));
      
      build_interaction_list(cell_com, cell_pos,
			     &interaction_node_list[interaction_node_offset[index].x],
			     &interaction_leaf_list[interaction_leaf_offset[index].x]);
  }
}

/**************************************************
***************************************************
***                                             ***
**   evaluate gravity via the interaction list  ***
***                                             ***
***************************************************
***************************************************/

/***************************/
/*  body-body interaction  */
/***************************/

__device__ float4 body_body_interaction(float4 grav, float4 body_i, float4 body_j) {
  float3 dr;
  dr.x = body_i.x - body_j.x;
  dr.y = body_i.y - body_j.y;
  dr.z = body_i.z - body_j.z;
  
  float ds2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
  float inv_ds = rsqrtf(ds2 + softening_squared) * (ds2 != 0.0f);
    
  float inv_s3 = body_j.w * inv_ds*inv_ds*inv_ds;
  grav.x -= inv_s3 * dr.x;
  grav.y -= inv_s3 * dr.y;
  grav.z -= inv_s3 * dr.z;
  grav.w -= body_j.w * inv_ds;
  return grav;

}


/***************************/
/*  body-node Octupole interaction  */
/***************************/

__device__ float4 body_node_Octupole(float4 grav,
				     float4 body_i, 
				     float4 com,
				     float4 Oct1, 
				     float4 Oct2,
				     float2 Oct3) {
  float3 dr;
  dr.x = body_i.x - com.x;    // 1 FLOP
  dr.y = body_i.y - com.y;    // 1 FLOP
  dr.z = body_i.z - com.z;    // 1 FLOP
  float ds2     = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;   // 5 FLOP
  float inv_ds  = rsqrt(ds2 + softening_squared) * (ds2 != 0.0f);      // 3 FLOP
  float inv_ds2 = inv_ds*inv_ds;                       // 1 FLOP
  float inv_ds3 = inv_ds *inv_ds2;                     // 1 FLOP
  float inv_ds5 = inv_ds3*inv_ds2;                     // 1 FLOP
  float inv_ds7 = 0.5f*inv_ds5*inv_ds2;                 // 2 FLOP

  float SijRj1 = Oct1.x*dr.x + Oct1.y*dr.y + Oct1.z*dr.z;   // 5 FLOP
  float SijRj2 = Oct2.x*dr.x + Oct2.y*dr.y + Oct2.z*dr.z;   // 5 FLOP
  float SijRj3 = Oct1.w*dr.x + Oct2.w*dr.y + Oct3.x*dr.z;   // 5 FLOP
  float SijRjRi_sq = SijRj1 * dr.x*dr.x + SijRj2 * dr.y*dr.y + SijRj3 * dr.z*dr.z;  // 8 FLOP


  /******************/
  /***  POTENTIAL ***/
  /******************/
 
  float pot = inv_ds7 * (SijRjRi_sq + Oct3.y*dr.x*dr.y*dr.z);  // 5 FLOP
  grav.w -= pot;   // 1 FLOP
  
  /*********************/
  /***  ACCELERATION ***/
  /*********************/

        /*** part 1 ***/

  float3 grav0 = {0.0f,0.0f,0.0f};
  grav0.x -= 7.0f*inv_ds2 * dr.x * pot;    //  4 FLOP
  grav0.y -= 7.0f*inv_ds2 * dr.y * pot;    //  4 FLOP
  grav0.z -= 7.0f*inv_ds2 * dr.z * pot;    //  4 FLOP

        /*** part 2 ***/

        /* S11*dx^2 + S21*dy^2 + S31*dz^2 */
        /* S12*dx^2 + S22*dy^2 + S32*dz^2 */
        /* S13*dx^2 + S23*dy^2 + S33*dz^2 */
  grav0.x += inv_ds7 * (2.0f*SijRj1*dr.x + dr.x*dr.x*Oct1.x + dr.y*dr.y*Oct2.x + dr.z*dr.z*Oct1.w); // 13 FLOP
  grav0.y += inv_ds7 * (2.0f*SijRj2*dr.y + dr.x*dr.x*Oct1.y + dr.y*dr.y*Oct2.y + dr.z*dr.z*Oct2.w); // 13 FLOP
  grav0.y += inv_ds7 * (2.0f*SijRj3*dr.z + dr.x*dr.x*Oct1.z + dr.y*dr.y*Oct2.z + dr.z*dr.z*Oct3.x); // 13 FLOP
  
        /*** part 2 ***/
  
  grav0.x += inv_ds7*Oct3.y * dr.y*dr.z;    // 4 FLOP
  grav0.y += inv_ds7*Oct3.y * dr.z*dr.x;    // 4 FLOP
  grav0.z += inv_ds7*Oct3.y * dr.x*dr.y;    // 4 FLOP

  grav.x += grav0.x;
  grav.y += grav0.y;
  grav.z += grav0.z;

  // TOTAL 108 FLOP

  return grav;
}

__device__ float4 evaluate_body_node_Octupole(float4 acc,
					      float4 body_pos,
					      int  &n_inter,
					      int2 list_len) {
  extern __shared__ float4 shared_com[];
  float4 *shared_Oct1 = &shared_com[blockDim.x];
  float4 *shared_Oct2 = &shared_Oct1[blockDim.x];
  float2 *shared_Oct3 = (float2*)&shared_Oct2[blockDim.x];

  n_inter = 0;
  for (int i = list_len.x; i < list_len.x + list_len.y; i += blockDim.x) {
    
    int node = tex1Dfetch(interaction_node_tex, i + threadIdx.x);
    if ( (node < 0) || (node >= n_nodes) ) node = 0;
    shared_com[threadIdx.x] = tex1Dfetch(node_com_tex, node);
    shared_Oct1[threadIdx.x] = tex1Dfetch(Oct1_tex, node);
    shared_Oct2[threadIdx.x] = tex1Dfetch(Oct2_tex, node);
    shared_Oct3[threadIdx.x] = tex1Dfetch(Oct3_tex, node);

    if (i + threadIdx.x >= list_len.x + list_len.y) {
      float4 null4 = {0,0,0,0};
      float2 null2 = {0,0};
      shared_Oct1[threadIdx.x] = null4;	
      shared_Oct2[threadIdx.x] = null4;
      shared_Oct3[threadIdx.x] = null2;
    }
    
    __syncthreads();
    
    /* check for body-node interaction */
    for (int j = 0; j < blockDim.x; j++) {
      n_inter++;
      acc = body_node_Octupole(acc, body_pos, shared_com[j],
				  shared_Oct1[j], shared_Oct2[j], shared_Oct3[j]);
    }
    __syncthreads();

  }

  return acc;
}

/***************************/
/*  body-node interaction  */
/***************************/

__device__ float4 body_node_interaction(float4 grav,
					float4 body_i, 
					float4 com, 
					float4 Qu, float4 Qd) {
  float3 dr;
  dr.x = body_i.x - com.x;    // 1 FLOP
  dr.y = body_i.y - com.y;    // 1 FLOP
  dr.z = body_i.z - com.z;    // 1 FLOP
  float ds2     = (((dr.x*dr.x) + dr.y*dr.y) + dr.z*dr.z);   // 5 FLOP
  float inv_ds  = rsqrt(ds2 + softening_squared) * (ds2 != 0.0f) ;      // 3 FLOP
  float inv_ds2 = inv_ds*inv_ds;                       // 1 FLOP
  float inv_ds3 = inv_ds *inv_ds2;                     // 1 FLOP
  float inv_ds5 = inv_ds3*inv_ds2;                     // 1 FLOP

  /************/
  /* potential */
  /************/

  grav.w -= com.w * inv_ds;    // 2 FLOP
  float Qy0 = inv_ds5 * (Qd.x*dr.x + Qu.x*dr.y + Qu.y*dr.z);  // 6 FLOP
  float Qy1 = inv_ds5 * (Qu.x*dr.x + Qd.y*dr.y + Qu.z*dr.z);  // 6 FLOP
  float Qy2 = inv_ds5 * (Qu.y*dr.x + Qu.z*dr.y + Qd.z*dr.z);  // 6 FLOP
  float yQy = Qy0 * dr.x + Qy1 * dr.y + Qy2 * dr.z;           // 5 FLOP

  grav.w -= 0.5f * yQy;    // 2 FLOP

  /* acceleartion */
  yQy = com.w * inv_ds3 + inv_ds2*2.5f * yQy;     // 4 FLOP

  grav.x += Qy0 - yQy * dr.x;  // 3 FLOPS
  grav.y += Qy1 - yQy * dr.y;  // 3 FLOPS
  grav.z += Qy2 - yQy * dr.z;  // 3 FLOPS


  // TOTAL 54 FLOP

  return grav;
}

__device__ float4 evaluate_body_node(float4 acc,
				     float4 body_pos,
				     int  &n_inter,
				     int2 list_len,
				     int n_in_cell) {
  extern __shared__ float4 shared_com[];
  float4 *shared_Qu = &shared_com[blockDim.x];
  float4 *shared_Qd = &shared_Qu [blockDim.x];
  
  n_inter = 0;
  
  int i_thread     = threadIdx.x/n_in_cell;
  int n_threads    = blockDim.x/n_in_cell;
  int n_per_thread = blockDim.x/n_threads;
  int j0 =  i_thread    * n_per_thread;
  int j1 = (i_thread+1) * n_per_thread;
  if (i_thread + 1 == n_threads) 
    j1 = blockDim.x;

  for (int i = list_len.x; i < list_len.x + list_len.y; i += blockDim.x) {
    
    int node = tex1Dfetch(interaction_node_tex, i + threadIdx.x);
    if ( (node < 0) || (node >= n_nodes) ) node = 0;
    shared_com[threadIdx.x] = tex1Dfetch(node_com_tex, node);
    shared_Qu[threadIdx.x]  = tex1Dfetch(node_Qu_tex, node);
    shared_Qd[threadIdx.x]  = tex1Dfetch(node_Qd_tex, node);

    if (i + threadIdx.x >= list_len.x + list_len.y) {
      float4 null4 = {0.0f,0.0f,0.0f,0.0f};
      shared_com[threadIdx.x] = null4;	
      shared_Qu [threadIdx.x] = null4;
      shared_Qd [threadIdx.x] = null4;
    }
    
    __syncthreads();
    
    /* check for body-node interaction */
    for (int j = j0; j < j1; j++) {
      n_inter++;
      acc = body_node_interaction(acc, body_pos, 
				  shared_com[j],
				  shared_Qu[j], shared_Qd[j]);
    }
    __syncthreads();
    
  }

  /*** now combine accelarations ****/
  int *n_inter_sh = (int*)&shared_com[blockDim.x + 1];
  shared_com[threadIdx.x] = acc;
  n_inter_sh[threadIdx.x] = n_inter;
  __syncthreads();

  if (threadIdx.x < n_in_cell) {
    for (int i = n_in_cell + threadIdx.x; i < n_in_cell*n_threads; i += n_in_cell) {
      float4 acc1 = shared_com[i];

      acc.x += acc1.x;
      acc.y += acc1.y;
      acc.z += acc1.z;
      acc.w += acc1.w;
    }

    for (int i = n_in_cell + threadIdx.x; i < blockDim.x; i += n_in_cell) {
      n_inter += n_inter_sh[i];
    }
  }
  __syncthreads();
  

  return acc;
}

__device__ float4 evaluate_body_leaf(float4 acc,
				     float4 body_pos,
				     int &n_inter,
				     int2 list_len) {
  extern __shared__ int shared_offset[];
  int *shared_len    = (int*)&shared_offset[blockDim.x];
  float4 *shared_pos = (float4*)&shared_len[blockDim.x];
  
  n_inter = 0;
  int tile   = 0;
  for (int i = list_len.x; i < list_len.x + list_len.y; i += blockDim.x, tile++) {
    int node_id                = tex1Dfetch(interaction_leaf_tex, i + threadIdx.x);
    shared_len   [threadIdx.x] = tex1Dfetch(n_in_node_tex,          node_id);
    shared_offset[threadIdx.x] = tex1Dfetch(node_bodies_offset_tex, node_id);
    __syncthreads();

    int j = min(blockDim.x, list_len.y - tile*blockDim.x);
    while (j-- > 0) {
      int len = shared_len[j];
      __syncthreads();

      shared_pos[threadIdx.x] = tex1Dfetch(bodies_pos_tex, shared_offset[j] + threadIdx.x);
      __syncthreads();

      while(len-- > 0) {
	n_inter++;
	acc = body_body_interaction(acc, body_pos, shared_pos[len]);
      }
      __syncthreads();
    }

    __syncthreads();
    
  }

  return acc;
}


__global__ void dev_evaluate_gravity_node(int cell_offset,
					  float4 *grav_acc,
					  int    *n_interactions,
					  int2   *interaction_node_len) {
  
  int cellId = cell_offset + blockIdx.x;
  bool write_flag = true;
  if (cellId >= n_cells) {
    cellId = blockIdx.x;
    write_flag = false;
  }
      
  int index     = tex1Dfetch(cell_bodies_offset_tex, cellId);
  int n_in_cell = tex1Dfetch(n_in_cell_tex, cellId);
  
  float4 body_pos = tex1Dfetch(bodies_pos_tex, index + threadIdx.x%n_in_cell);

  float4 acc = {0,0,0,0};

  int n_inter;
#ifdef OCTUPOLE
  acc = evaluate_body_node_Octupole(acc, body_pos, n_inter,
 				    interaction_node_len[cellId]);
#endif

#ifdef QUADRUPOLE
  acc = evaluate_body_node(acc, body_pos, n_inter,
			   interaction_node_len[cellId],
			   n_in_cell);
#endif
  
  if (threadIdx.x < n_in_cell) {
    if (write_flag) {
      grav_acc[index + threadIdx.x] = acc;
//       fprintf(stderr, "cellId= %d  index= %d  n_in_cell= %d\n",
// 	      cellId, index + threadIdx.x, n_in_cell);
//       fprintf(stderr, " acc= [%f %f %f %f]\n", acc.x, acc.y, acc.z, acc.w);
    }
     n_interactions[index + threadIdx.x] = n_inter;
  }
}

__global__ void dev_evaluate_gravity_leaf(int cell_offset,
					  float4 *grav_acc,
					  int    *n_interactions,
					  int2   *interaction_leaf_len) {
  

  int cellId = cell_offset + blockIdx.x;
  bool write_flag = true;
  if (cellId >= n_cells) {
    cellId = blockIdx.x;
    write_flag = false;
  }
      
  int index     = tex1Dfetch(cell_bodies_offset_tex, cellId);
  int n_in_cell = tex1Dfetch(n_in_cell_tex, cellId);

  float4 body_pos = tex1Dfetch(bodies_pos_tex, index + threadIdx.x%n_in_cell);
  
  float4 acc = grav_acc[index + threadIdx.x%n_in_cell];
  
  int n_inter;
  acc = evaluate_body_leaf(acc, body_pos, n_inter,
			   interaction_leaf_len[cellId]);

  if (threadIdx.x < n_in_cell) {
    if (write_flag)
      grav_acc[index + threadIdx.x] = acc;
    n_interactions[index + threadIdx.x] = n_inter;
  }
}

#endif 
