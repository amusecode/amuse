#include "support_kernels.cu"
#include "dev_shared_traverse_functions.cu"

#include <stdio.h>


texture<float4, 1, cudaReadModeElementType> texNodeSize;
texture<float4, 1, cudaReadModeElementType> texNodeCenter;
texture<float4, 1, cudaReadModeElementType> texMultipole;
texture<float4, 1, cudaReadModeElementType> texBody;

__device__ int ngb_cnt(float3 pos_i,
                       float h_i,
                       float3 pos_j)
{
  
  float3 dr = {pos_i.x - pos_j.x,
               pos_i.y - pos_j.y,
               pos_i.z - pos_j.z};
  float ds2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
  
  return (ds2 <= h_i*h_i);
}

// __device__ bool split_node_ngb(float4 node_pos,
//                                float4 group_pos)
// {
//   float s = node_pos.w + group_pos.w;
//   
//   float3 dr = {fabs(group_pos.x - node_pos.x),
//                fabs(group_pos.y - node_pos.y),
//                fabs(group_pos.z - node_pos.z)};
//   
//   return ((dr.x < s) && (dr.y < s) && (dr.z < s));
// }

 __device__ bool split_node_ngb(float4 nodeCenter, float4 nodeSize, float4 groupCenter, float4 groupSize)
{
  //Compute the distance between the group and the cell
  float3 dr = {fabs(groupCenter.x - nodeCenter.x) - (groupSize.x + nodeSize.x),
               fabs(groupCenter.y - nodeCenter.y) - (groupSize.y + nodeSize.y),
               fabs(groupCenter.z - nodeCenter.z) - (groupSize.z + nodeSize.z)};

  dr.x += fabs(dr.x); dr.x *= 0.5f;
  dr.y += fabs(dr.y); dr.y *= 0.5f;
  dr.z += fabs(dr.z); dr.z *= 0.5f;

  //Distance squared, no need to do sqrt since opening criteria has been squared
  float ds2    = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

  return (ds2 <= groupCenter.w);
}



#define TEXTURES
#define OLDPREFIX
#define DOGRAV


template<int DIM2, int SHIFT>
__device__ uint get_n_ngb(int DIM2x, int DIM2y,
				      int tid, int tx, int ty,
				      int body_i, float4 pos_i,
				      real4 group_pos,
                                      uint2 node_begend,
                                      real4 *body_pos,
				      int *shmem,
                                      int *lmem,
				      int &ngb,
                                      int &direCount,
                                      volatile float4 *boxSizeInfo,
                                      float4 groupSize,
				      volatile float4 *boxCenterInfo,
                                      float &ds2_min,
                                      float h_group) {

  ngb = -1;
  
  int n_ngb = 0;
  
// float ds2_min = 1.0e10f;
  ds2_min = 1.0e10f;


  /*********** set necessary thread constants **********/

  const int DIMx = 1  << DIM2x;
  const int DIMy = 1  << DIM2y;
  const int DIM  = 1  << DIM2;
  const int offs = ty << DIM2x;

  /*********** shared memory distribution **********/

                                                  //  begin,    end,   size
                                                  // -----------------------
  int *approx = (int*)&shmem [     0];            //  0*DIM,  2*DIM,  2*DIM
  int *direct = (int*)&approx[ 2*DIM];            //  2*DIM,  3*DIM,  1*DIM
  int *nodes  = (int*)&direct[   DIM];            //  3*DIM, 13*DIM, 10*DIM
  int *prefix = (int*)&nodes [10*DIM];            // 13*DIM, 15*DIM,  2*DIM

  int    *body_list = (int*   )&nodes    [  DIM]; //  4*DIM,  8*DIM,  4*DIM
  float  *sh_mass   = (float* )&body_list[4*DIM]; //  8*DIM,  9*DIM,  1*DIM
  float3 *sh_pos    = (float3*)&sh_mass  [  DIM]; //  9*DIM, 12*DIM   3*DIM

  int    *sh_jid    = (int*  )&sh_pos[DIM];
  
  //Reduction at the end
  float  *sh_ds2    = (float*)&shmem[DIM];
  int    *sh_ngb    = (int*  )&sh_ds2[DIM];


  /*********** stack **********/

  int *nstack = lmem;

  /*********** begin tree-walk **********/

  int n_direct = 0;


  for (int root_node = node_begend.x; root_node < node_begend.y; root_node += DIM) {
    int n_nodes0 = min(node_begend.y - root_node, DIM);
    int n_stack0 = 0;
    int n_stack_pre = 0;

    { nstack[ACCS<SHIFT>(n_stack0)] = root_node + tid;   n_stack0++; }

    /*********** walk each level **********/
    while (n_nodes0 > 0) {


      int n_nodes1 = 0;
      int n_offset = 0;

      int n_stack1 = n_stack0;
      int c_stack0 = n_stack_pre;

      /*********** walk a level **********/
      while(c_stack0 < n_stack0) {

	/***
	**** --> fetch the list of nodes rom LMEM
	***/
	bool use_node = tid <  n_nodes0;
 	{ prefix[tid] = nstack[ACCS<SHIFT>(c_stack0)];   c_stack0++; }
	__syncthreads();
	int node  = prefix[min(tid, n_nodes0 - 1)];

        if(n_nodes0 > 0){       //Work around pre 4.1 compiler bug
          n_nodes0 -= DIM;
        }

	/***
	**** --> process each of the nodes in the list in parallel
	***/

        #ifndef TEXTURES
          float4 nodeSize = get_float4(boxSizeInfo[node]);                   //Fetch the size of the box. Size.w = child info
          float4 node_pos = get_float4(boxCenterInfo[node]);                 //Fetch the center of the box. center.w = opening info
        #else
          float4 nodeSize =  tex1Dfetch(texNodeSize, node);
          float4 node_pos =  tex1Dfetch(texNodeCenter, node);
        #endif

        int node_data = __float_as_int(nodeSize.w);

        //TODO Fix This
        group_pos.w = h_group;        //The looking radus

        //Check if a cell has to be opened
        bool   split   = split_node_ngb(node_pos, nodeSize, group_pos, groupSize);          //Check if node should be split
        bool leaf       = node_pos.w <= 0;  //Small AND equal incase of a 1 particle cell       //Check if it is a leaf
//         split = true;


	uint mask    = BTEST((split && !leaf) && use_node);               // mask = #FFFFFFFF if use_node+split+not_a_leaf==true, otherwise zero
        int child    =    node_data & 0x0FFFFFFF;                         //Index to the first child of the node
        int nchild   = (((node_data & 0xF0000000) >> 28)) & mask;         //The number of children this node has


  	/***
	**** --> calculate prefix
	***/

	int *prefix0 = &prefix[  0];
	int *prefix1 = &prefix[DIM];

        #ifdef OLDPREFIX
          int n_total = calc_prefix<DIM2>(prefix, tid,  nchild);
          prefix[tid] += n_offset - nchild;
          __syncthreads();
        #else
          inclusive_scan_block<ADDOP<int>, int>(prefix, nchild, tid);        // inclusive scan to compute memory offset of each child
          int n_total = prefix[blockDim.x - 1];                              // fetch total number of children, i.e. offset of the last child -1
          __syncthreads();                                                   // thread barrier to make sure that warps completed their jobs
          prefix[tid] += n_offset - nchild;                                  // convert inclusive into exclusive scan for referencing purpose
          __syncthreads();                                                   // thread barrier
        #endif

        for (int i = n_offset; i < n_offset + n_total; i += DIM)         //nullify part of the array that will be filled with children
            nodes[tid + i] = 0;                                          //but do not touch those parts which has already been filled
        __syncthreads();                                                 //Thread barrier to make sure all warps finished writing data

        bool flag = (split && !leaf) && use_node;                        //Flag = use_node + split + not_a_leaf;Use only non_leaf nodes that are to be split
        if (flag) nodes[prefix[tid]] = child;                            //Thread with the node that is about to be split
        __syncthreads();                                                 //writes the first child in the array of nodes

        /*** in the following 8 lines, we calculate indexes of all the children that have to be walked from the index of the first child***/
        if (flag && nodes[prefix[tid] + 1] == 0) nodes[prefix[tid] + 1] = child + 1; __syncthreads();
	if (flag && nodes[prefix[tid] + 2] == 0) nodes[prefix[tid] + 2] = child + 2; __syncthreads();
	if (flag && nodes[prefix[tid] + 3] == 0) nodes[prefix[tid] + 3] = child + 3; __syncthreads();
	if (flag && nodes[prefix[tid] + 4] == 0) nodes[prefix[tid] + 4] = child + 4; __syncthreads();
	if (flag && nodes[prefix[tid] + 5] == 0) nodes[prefix[tid] + 5] = child + 5; __syncthreads();
	if (flag && nodes[prefix[tid] + 6] == 0) nodes[prefix[tid] + 6] = child + 6; __syncthreads();
	if (flag && nodes[prefix[tid] + 7] == 0) nodes[prefix[tid] + 7] = child + 7; __syncthreads();

	n_offset += n_total;    //Increase the offset in the array by the number of newly added nodes


	/***
	**** --> save list of nodes to LMEM
	***/

        /*** if half of shared memory or more is filled with the the nodes, dump these into slowmem stack ***/
	while(n_offset >= DIM) {
	  n_offset -= DIM;
	  const int offs1 = ACCS<SHIFT>(n_stack1);
	  nstack[offs1]   = nodes[n_offset + tid];   n_stack1++;
	  n_nodes1       += DIM;

          if((n_stack1 - c_stack0) >= (LMEM_STACK_SIZE << SHIFT))
          {
            //We overwrote our current stack
	    direCount = -1; return 0;	 
          }
	}

	__syncthreads();


#if 1
	/***********************************/
	/******       DIRECT          ******/
	/***********************************/

        int *sh_body = &approx[DIM];

        flag         = split && leaf && use_node;                                //flag = split + leaf + use_node
        int  jbody   = node_data & BODYMASK;                                     //the first body in the leaf
        int  nbody   = (((node_data & INVBMASK) >> LEAFBIT)+1) & BTEST(flag);    //number of bodies in the leaf masked with the flag

        body_list[tid] = direct[tid];                                            //copy list of bodies from previous pass to body_list
        sh_body  [tid] = jbody;                                                  //store the leafs first body id into shared memory

        // step 1
        #ifdef OLDPREFIX
          calc_prefix<DIM2>(prefix0, tid, flag);
        #else
          inclusive_scan_block<ADDOP<int>, int>(prefix0, (int)flag, tid);       // inclusive scan on flags to construct array
        #endif

        if (flag) prefix1[prefix0[tid] - 1] = tid;                             //with tidś whose leaves have to be opened
          __syncthreads();                                                      //thread barrier, make sure all warps completed the job

        // step 2
        #ifdef OLDPREFIX
          int n_bodies  = calc_prefix<DIM2>(prefix0, tid, nbody);
        #else
          inclusive_scan_block<ADDOP<int>, int>(prefix0, nbody, tid);        // inclusive scan to compute memory offset for each body
          int n_bodies = prefix0[blockDim.x - 1];                            //Total number of bides extract from the leaves
          __syncthreads();                                                   // thread barrier to make sure that warps completed their jobs
        #endif

	direct [tid]  = prefix0[tid];                                       //Store a copy of inclusive scan in direct
	prefix0[tid] -= nbody;                                              //convert inclusive int oexclusive scan
	prefix0[tid] += 1;                                                  //add unity, since later prefix0[tid] == 0 used to check barrier

	int nl_pre = 0;                                                     //Number of leaves that have already been processed

        #define NJMAX (DIM*4)
	while (n_bodies > 0) {
	  int nb    = min(n_bodies, NJMAX - n_direct);                    //Make sure number of bides to be extracted does not exceed
                                                                          //the amount of allocated shared memory

          // step 0                                                      //nullify part of the body_list that will be filled with bodies
	  for (int i = n_direct; i < n_direct + nb; i += DIM){           //from the leaves that are being processed
            body_list[i + tid] = 0;
          }
	  __syncthreads();

          //step 1:
	  if (flag && (direct[tid] <= nb) && (prefix0[tid] > 0))        //make sure that the thread indeed carries a leaf
	    body_list[n_direct + prefix0[tid] - 1] = 1;                 //whose bodies will be extracted
	  __syncthreads();

          //step 2:
         #ifdef OLDPREFIX
            int nl = calc_prefix<DIM2>(nb, &body_list[n_direct], tid);
          #else
            int nl = inclusive_scan_array<ADDOP<int>, int>              // inclusive scan to compute number of leaves to process
                            (&body_list[n_direct], nb, tid);            // to make sure that there is enough shared memory for bodies
          #endif
	  nb = direct[prefix1[nl_pre + nl - 1]];                        // number of bodies stored in these leaves

	  // step 3:
	  for (int i = n_direct; i < n_direct + nb; i += DIM) {          //segmented fill of the body_list
	    int j = prefix1[nl_pre + body_list[i + tid] - 1];            // compute the first body in shared j-body array
	    body_list[i + tid] = (i + tid - n_direct) -                 //add to the index of the first j-body in a child
                                 (prefix0[j] - 1) + sh_body[j];         //the index of the first child in body_list array
	  }
	  __syncthreads();


         /**************************************************
          *  example of what is accomplished in steps 0-4   *
          *       ---------------------------               *
          * step 0: body_list = 000000000000000000000       *
          * step 1: body_list = 100010001000000100100       *
          * step 2: body_list = 111122223333333444555       *
          * step 3: body_list = 012301230123456012012       *
          *         assuming that sh_body[j] = 0            *
         ***************************************************/

	  n_bodies     -= nb;                                   //subtract from n_bodies number of bodies that have been extracted
   	  nl_pre       += nl;                                   //increase the number of leaves that where processed
	  direct [tid] -= nb;                                   //subtract the number of extracted bodies in this pass
 	  prefix0[tid] = max(prefix0[tid] - nb, 0);             //same here, but do not let the number be negative (GT200 bug!?)
	  n_direct     += nb;                                  //increase the number of bodies to be procssed

	  while(n_direct >= DIM) {
	    n_direct -= DIM;


	    float4 posj  = body_pos[body_list[n_direct + tid]];
//             float4 posj  = tex1Dfetch(texBody, body_list[n_direct + tid]);
	    sh_pos [tid] = (float3){posj.x, posj.y, posj.z};
	    sh_jid [tid] = body_list[n_direct + tid];

            __syncthreads();
#pragma unroll
	    for (int j = 0; j < DIMx; j++)
            {
              //TODO should we check on selfGrav?
              int selfGrav = (body_i != sh_jid[offs + j]);
// 	      if (body_i != sh_jid[offs + j]) //If statement replaced by multiplication
              {
                direCount++;
#ifdef DONGBCOUNT
                n_ngb += ngb_count(pos_i, h_i, sh_pos[offs + j]);
#endif
              }
	    }//End for j < DIMx
	    __syncthreads();
	  }// end while n_direct >= DIM
	}// end while n_bodies > 0
	direct[tid] = body_list[tid];
	__syncthreads();
#endif
      } //end lvl


      n_nodes1 += n_offset;
      if (n_offset > 0)
      { 
        nstack[ACCS<SHIFT>(n_stack1)] = nodes[tid];   n_stack1++; 
        if((n_stack1 - c_stack0) >= (LMEM_STACK_SIZE << SHIFT))
        {
          //We overwrote our current stack
	  direCount = -1; return 0;	 
	}
      }
      __syncthreads();


      /***
      **** --> copy nodes1 to nodes0: done by reassigning the pointers
      ***/
      n_nodes0    = n_nodes1;

      n_stack_pre = n_stack0;
      n_stack0    = n_stack1;

    }//end while   levels
  }//end for

  if(n_direct > 0)
  {
    if (tid < n_direct) {

      float4 posj = body_pos[direct[tid]];

  //     float4 posj  = tex1Dfetch(texBody, direct[tid]);
      sh_pos [tid] = (float3){posj.x, posj.y, posj.z};
      sh_jid [tid] = direct[tid];
    } else {
      sh_jid [tid] = -1;
      sh_pos [tid] = (float3){1.0e10f, 1.0e10f, 1.0e10f};
    }

    __syncthreads();
  #pragma unroll
    for (int j = 0; j < DIMx; j++) {
      if ((sh_jid[offs + j] >= 0)) {
        int selfGrav = (body_i != sh_jid[offs + j]);
        direCount++;
  #ifdef DONGBCOUNT
          n_ngb += ngb_count(pos_i, h_i, sh_pos[offs + j]);          
  #endif
      }
    }
    __syncthreads();
  }

  /***
  **** --> reduce data between threads
  ***/

  //Sum the interaction counters
  //and the number of neighbours
  sh_ds2[tid] = direCount;
  sh_ngb[tid] = n_ngb;
  __syncthreads();


  if (ty == 0) {
  #pragma unroll
    for (int i = 1; i < DIMy; i++){
        int idx = (i << DIM2x) + tx;
        direCount += sh_ds2[idx];
        n_ngb     += sh_ngb[idx];
    }
  }
  __syncthreads();

  return n_ngb;
}

 
extern "C" __global__ void
__launch_bounds__(NTHREAD)
dev_get_n_ngb(const int n_active_groups,
              int    n_bodies,
              uint2 node_begend,                        
              real4  *body_pos,                        
              real4  *group_body_pos, 
              int    *ngb_out,
              int    *active_inout,
              int2   *interactions,
              uint2  *group_list,
              float4  *boxSizeInfo,                        
              float4  *boxCenterInfo,                        
              int     *MEM_BUF) {


  const int blockDim2 = NTHREAD2;
   __shared__ int shmem[15*(1 << blockDim2)];
//    __shared__ int shmem[24*(1 << blockDim2)]; is possible on FERMI
//    int             lmem[LMEM_STACK_SIZE];



  /*********** check if this block is linked to a leaf **********/

  int bid = gridDim.x * blockIdx.y + blockIdx.x;

  while(true)
  {

    if(threadIdx.x == 0)
    {
      bid         = atomicAdd(&active_inout[n_bodies], 1);
      shmem[0]    = bid;
    }
  __syncthreads();

  bid   = shmem[0];

  if (bid >= n_active_groups) return;


  int tid = threadIdx.y * blockDim.x + threadIdx.x;

  int grpOffset = 0;

//   volatile int *lmem = &MEM_BUF[blockIdx.x*LMEM_STACK_SIZE*blockDim.x + threadIdx.x*LMEM_STACK_SIZE];
//   int *lmem = &MEM_BUF[blockIdx.x*LMEM_STACK_SIZE*blockDim.x + threadIdx.x*LMEM_STACK_SIZE];
   int *lmem = &MEM_BUF[blockIdx.x*LMEM_STACK_SIZE*blockDim.x];


  /*********** set necessary thread constants **********/

    uint2 grpInfo =  group_list[bid];
    uint body_i   =  grpInfo.x;
    uint nb_i     = (grpInfo.y - grpInfo.x) + 1;

    int DIM2x = 0;
    while (((nb_i - 1) >> DIM2x) > 0) DIM2x++;

    DIM2x     = max(DIM2x,4);
    int DIM2y = blockDim2 - DIM2x;

    int tx = tid & ((1 << DIM2x) - 1);
    int ty = tid >> DIM2x;

    body_i += tx%nb_i;


    //float4 pos_i = tex1Dfetch(bodies_pos_ref, body_i);   // texture read: 4 floats


    float4 pos_i = group_body_pos[body_i];

    real4 group_pos;
    real4 curGroupSize;

    computeGroupProps(group_pos, curGroupSize, pos_i, shmem);

    int  ngb_i;
    uint n_ngb;
    float ds2;
    float h_group = 0;

    int direCount = 0;

    n_ngb = get_n_ngb<blockDim2, 0>( DIM2x, DIM2y, tid, tx, ty,
                                            body_i, pos_i, group_pos,
                                            node_begend, body_pos,
                                            shmem, lmem, ngb_i, direCount, boxSizeInfo,
                                            curGroupSize, boxCenterInfo,ds2, 
                                            h_group);
    if(direCount < 0)
    {

      //Try to get access to the big stack, only one block per time is allowed
      if(threadIdx.x == 0)
      {
        int res = atomicExch(&active_inout[n_bodies+1], 1); //If the old value (res) is 0 we can go otherwise sleep
        int waitCounter  = 0;
        while(res != 0)
        {
            //Sleep
            for(int i=0; i < (1024); i++)
            {
                    waitCounter += 1;
            }
            //Test again
            shmem[0] = waitCounter;
            res = atomicExch(&active_inout[n_bodies+1], 1); 
        }
      }

      __syncthreads();

      lmem = &MEM_BUF[gridDim.x*LMEM_STACK_SIZE*blockDim.x];    //Use the extra large buffer
      direCount = 0;
      n_ngb = get_n_ngb<blockDim2, 8>(DIM2x, DIM2y, tid, tx, ty,
                                                body_i, pos_i, group_pos,
                                                node_begend, body_pos,
                                                shmem, lmem, ngb_i, direCount, boxSizeInfo,
                                                curGroupSize, boxCenterInfo,ds2,
                                                h_group);

      lmem = &MEM_BUF[blockIdx.x*LMEM_STACK_SIZE*blockDim.x]; //Back to normal location

      if(threadIdx.x == 0)
      {
              atomicExch(&active_inout[n_bodies+1], 0); //Release the lock
      }
    }//end if apprCount < 0

    if (tid < nb_i) {
        ngb_out     [body_i] = ngb_i;
        active_inout[body_i] = 1;
        interactions[body_i].x = n_ngb;
        interactions[body_i].y = direCount ;
      }
  }     //end while
}


