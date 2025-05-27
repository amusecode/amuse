#include "support_kernels.cu"
#include "dev_shared_traverse_functions.cu"

#include <stdio.h>

//Some settings

#define TEXTURES
#define OLDPREFIX
#define DOGRAV



texture<float4, 1, cudaReadModeElementType> texNodeSize;
texture<float4, 1, cudaReadModeElementType> texNodeCenter;
texture<float4, 1, cudaReadModeElementType> texMultipole;
texture<float4, 1, cudaReadModeElementType> texBody;

#ifdef INDSOFT
__device__ float4 get_D04(float ds2, float epsP, float epsQ, int selfGrav = 1) {

//   float eps     = fmaxf(epsP, epsQ);
  float epseff  = epsP + epsQ;
  float ids, ids2, ids3;

  if(ds2 >= (epseff*epseff))
  {
     ids  = rsqrtf(ds2);
     //if(isnan(ids)) ids = 0; not needed if we use non-zero softening
     ids  = ids*selfGrav;       //Prevent selfGravity, instead of using if-statement
     ids3 = ids*ids*ids;
  }
  else
  {
    //these two lines are faster than a real sqrt
    float dist = ds2*selfGrav*rsqrtf(ds2); //Gives NaN is ds is 0
    if(isnan(dist)) dist = 0.0f;
    //float    dist = sqrtf(ds2); //Slower than the two lines above
    //assert(!isnan(dist));

    float epseffi = 1.f/epseff;
    float rhinv   = dist*epseffi;

    if(rhinv < 0.5f)
    {
      ids3  = 4.f/3.f + (rhinv*rhinv)*(4.f*rhinv-4.8f);
      ids   = 1.4f - (rhinv*rhinv)*(8.f/3.f+(rhinv*rhinv)*(3.2f*rhinv-4.8f));
    }
    else
    {
      ids3 = 8.f/3.f-6.f*rhinv+4.8f*(rhinv*rhinv)-4.f/3.f*(rhinv*rhinv*rhinv)-1.f/120.f/(rhinv*rhinv*rhinv);
      ids  = 1.6f-1/(30.f*rhinv)-(rhinv*rhinv)*(16.f/3.f+rhinv*(-8.f+rhinv*(4.8f-rhinv*16.f/15.f)));
    }//end if rhin < 0.5

    ids  = ids*2.f*epseffi;
    ids3 =  ids3*8.f*(epseffi*epseffi*epseffi);

    //Self gravity prevention (and NaN prevention)
    ids  *= selfGrav;
    ids3 *= selfGrav;

  } //end dist >= epseff
  
  ids2 = ids*ids;
  float ids5 = ids3*ids2;
  float ids7 = ids5*ids2;
  return (float4){ids, -ids3, +3.0f*ids5, -15.0f*ids7};
} 

#else

__device__ float4 get_D04(float ds2, int selfGrav = 1) {
#if 1
  float ids  = rsqrtf(ds2);  //Does not work with zero-softening
  //   if(isnan(ids)) ids = 0;               //This does work with zero-softening, few percent performance drop
  //float ids  = (1.0f / sqrtf(ds2)) * selfGrav; Slower in Pre CUDA4.1
  ids *= selfGrav;
#else
  const float ids = (ds2 > 0.0f) ? rsqrtf(ds2) : 0.0f;
#endif
  const float ids2 = ids*ids;
  float ids3 = ids *ids2;
  float ids5 = ids3*ids2;
  float ids7 = ids5*ids2;
  return (float4){ids, -ids3, +3.0f*ids5, -15.0f*ids7};
}  // 9 flops

#endif

#ifdef INDSOFT
__device__ float4 add_acc(float4 acc,  float4 pos,
                          float massj, float3 posj, float epsQ,
                          float &ds2, float eps2P, int selfGrav) {

#else

__device__ float4 add_acc(float4 acc,  float4 pos,
			  float massj, float3 posj,
			  float &ds2, float eps2, int selfGrav) {
#endif

  float3 dr = {pos.x - posj.x,
	       pos.y - posj.y,
	       pos.z - posj.z};

  #ifdef INDSOFT
    ds2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
    float4 D04 = get_D04(ds2, eps2P, epsQ, selfGrav);
  #else
    ds2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z + eps2;
    float4 D04 = get_D04(ds2, selfGrav);
  #endif

  float  D0  = D04.x*massj;
  float  D1  = D04.y*massj;

  acc.w -= D0;
  acc.x += D1*dr.x;
  acc.y += D1*dr.y;
  acc.z += D1*dr.z;

  return acc;
}

#ifdef INDSOFT
__device__ float4 add_acc(float4 acc, float4 pos,
                          float mass, float3 com,
                          float3 Q0,  float3 Q1, float epsNode, float eps2) {
#else

__device__ float4 add_acc(float4 acc, float4 pos,
			  float mass, float3 com,
			  float3 Q0,  float3 Q1, float eps2) {
#endif
  float3 dr = {pos.x - com.x,
	       pos.y - com.y,
	       pos.z - com.z};

  #ifdef INDSOFT
    float  ds2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
    float4 D04 = get_D04(ds2, epsNode, eps2);
  #else
    float  ds2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z + eps2;
    float4 D04 = get_D04(ds2);
  #endif

  float  D0  = D04.x*mass;
  float  D1  = D04.y*mass;
  float  D2  = D04.z*mass;
  float  D3  = D04.w*mass;

  float oct_q11 = Q0.x;
  float oct_q22 = Q0.y;
  float oct_q33 = Q0.z;
  float oct_q12 = Q1.x;
  float oct_q13 = Q1.y;
  float oct_q23 = Q1.z;

  float Qii = oct_q11 + oct_q22 + oct_q33;
  float QijRiRj =
          (oct_q11*dr.x*dr.x + oct_q22*dr.y*dr.y + oct_q33*dr.z*dr.z) +
    2.0f*(oct_q12*dr.y*dr.x + oct_q13*dr.z*dr.x + oct_q23*dr.y*dr.z);
//2.0f was
  //volatile float QijRiRj_1 = oct_q12*dr.y*dr.x + oct_q13*dr.z*dr.x + oct_q23*dr.y*dr.z;
 //QijRiRj = dr.x*dr.x;// + oct_q22*dr.y*dr.y + oct_q33*dr.z*dr.z;  //oct_q11;//*

  acc.w        -= D0 + 0.5f*D1*Qii + 0.5f*D2*QijRiRj;
  float C01a    = D1 + 0.5f*D2*Qii + 0.5f*D3*QijRiRj;
  acc.x         += C01a*dr.x + D2*(oct_q11*dr.x + oct_q12*dr.y + oct_q13*dr.z);
  acc.y         += C01a*dr.y + D2*(oct_q12*dr.x + oct_q22*dr.y + oct_q23*dr.z);
  acc.z         += C01a*dr.z + D2*(oct_q13*dr.x + oct_q23*dr.y + oct_q33*dr.z);

  return acc;
}


template<int DIM2, int SHIFT>
__device__ float4 approximate_gravity(int DIM2x, int DIM2y,
				      int tid, int tx, int ty,
				      int body_i, float4 pos_i,
				      real4 group_pos,
                                      float eps2,
                                      uint2 node_begend,
                                      real4 *multipole_data,
                                      real4 *body_pos,
				      int *shmem,
                                      int *lmem,
				      int &ngb,
                                      int &apprCount, int &direCount,
                                      volatile float4 *boxSizeInfo,
                                      float4 groupSize,
				      volatile float4 *boxCenterInfo,
                                      float group_eps,
                                      real4 *body_vel,
                                      float &ds2_min) {

  float4 acc_i = {0.0f, 0.0f, 0.0f, 0.0f};
  ngb = -1;
  
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

  float  *node_mon0 = (float* )&nodes    [DIM];   //  4*DIM,  5*DIM,  1*DIM
  float3 *node_mon1 = (float3*)&node_mon0[DIM];   //  5*DIM,  8*DIM,  3*DIM
  float3 *node_oct0 = (float3*)&node_mon1[DIM];   //  8*DIM, 11*DIM,  3*DIM
  float3 *node_oct1 = (float3*)&node_oct0[DIM];   // 11*DIM, 14*DIM,  3*DIM

  int    *body_list = (int*   )&nodes    [  DIM]; //  4*DIM,  8*DIM,  4*DIM
  float  *sh_mass   = (float* )&body_list[4*DIM]; //  8*DIM,  9*DIM,  1*DIM
  float3 *sh_pos    = (float3*)&sh_mass  [  DIM]; //  9*DIM, 12*DIM   3*DIM

  float  *sh_pot = sh_mass;
  float3 *sh_acc = sh_pos;

  int    *sh_jid    = (int*  )&sh_pos[DIM];
  float  *sh_ds2 = (float*)&sh_acc[DIM];
  int    *sh_ngb = (int*  )&sh_ds2[DIM];

  #ifdef INDSOFT
    //This works with shmem of dimx15
    float *sh_eps     = (float* )&sh_jid    [DIM];     //JB Partially overwrites the prefix part
    float  *node_eps  = (float*) &approx[ 2*DIM];   //JB 11*DIM, 14*DIM,  3*DIM

  #endif


  /*********** stack **********/

  int *nstack = lmem;

  /*********** begin tree-walk **********/

  int n_approx = 0;
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

        #ifdef INDSOFT
          //Very inefficient this but for testing I have to live with it...
          float node_eps_val = multipole_data[node*3 + 1].w;
        #endif

	//Check if a cell has to be opened
	#ifdef IMPBH
	  //Improved barnes hut method
          #ifndef TEXTURES
            float4 nodeCOM = multipole_data[node*3];
          #else
            float4 nodeCOM = tex1Dfetch(texMultipole,node*3);
          #endif

	  nodeCOM.w      = node_pos.w;
          #ifdef INDSOFT
            bool   split   = split_node_grav_impbh(nodeCOM, group_pos, groupSize, group_eps, node_eps_val);
          #else
            bool   split   = split_node_grav_impbh(nodeCOM, group_pos, groupSize);
          #endif

	#else
          //Minimum distance method
          #ifdef INDSOFT
            //bool   split   = split_node_grav_md(node_pos, nodeSize, group_pos, groupSize, group_eps, node_eps_val);  //Check if node should be split
            bool   split   = split_node_grav_md_rcut(node_pos, nodeSize, group_pos, groupSize, group_eps, node_eps_val, rcut_sq);  // modified by M.I.
          #else
            //bool   split   = split_node_grav_md(node_pos, nodeSize, group_pos, groupSize);  //Check if node should be split
            bool   split   = split_node_grav_md(node_pos, nodeSize, group_pos, groupSize, rcut_sq);  // modified by M.I.
          #endif
	#endif

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
	  nstack[offs1] = nodes[n_offset + tid];   n_stack1++;
	  n_nodes1 += DIM;

          if((n_stack1 - c_stack0) >= (LMEM_STACK_SIZE << SHIFT))
          {
            //We overwrote our current stack
	    apprCount = -1; return acc_i;	 
          }
	}

	__syncthreads();



	/******************************/
	/******************************/
	/*****     EVALUATION     *****/
	/******************************/
	/******************************/
#if 1
	/***********************************/
	/******       APPROX          ******/
	/***********************************/

        #ifdef OLDPREFIX
          n_total = calc_prefix<DIM2>(prefix, tid,  1 - (split || !use_node));
        #else
          inclusive_scan_block<ADDOP<int>, int>(prefix, 1 - (split || !use_node), tid);
          n_total = prefix[blockDim.x - 1];
        #endif


// 	n_total = calc_prefix<DIM2>(prefix, tid,  !split && use_node);         // for some unkown reason this does not work right on the GPU
	if (!split && use_node) approx[n_approx + prefix[tid] - 1] = node;
 	__syncthreads();
	n_approx += n_total;

 	while (n_approx >= DIM) {
 	  n_approx -= DIM;
	  int address      = (approx[n_approx + tid] << 1) + approx[n_approx + tid];
          #ifndef TEXTURES
            float4 monopole  = multipole_data[address    ];
            float4 octopole0 = multipole_data[address + 1];
            float4 octopole1 = multipole_data[address + 2];
          #else
            float4 monopole  = tex1Dfetch(texMultipole, address);
            float4 octopole0 = tex1Dfetch(texMultipole, address + 1);
            float4 octopole1 = tex1Dfetch(texMultipole, address + 2);
          #endif

 	  node_mon0[tid] = monopole.w;
 	  node_mon1[tid] = (float3){monopole.x,  monopole.y,  monopole.z};
 	  node_oct0[tid] = (float3){octopole0.x, octopole0.y, octopole0.z};
 	  node_oct1[tid] = (float3){octopole1.x, octopole1.y, octopole1.z};

          #ifdef INDSOFT
            float temp          = node_eps[tid]; //Backup value in the shmem into register
            node_eps[tid]       = octopole0.w;
          #endif

	  __syncthreads();
#pragma unroll
	  for (int i = 0; i < DIMx; i++)
          {
            apprCount++;
#ifdef DOGRAV
            #ifdef INDSOFT
              acc_i = add_acc(acc_i, pos_i,
                              node_mon0[offs + i], node_mon1[offs + i],
                              node_oct0[offs + i], node_oct1[offs + i], node_eps[offs + i], eps2);

            #else
              acc_i = add_acc(acc_i, pos_i,
                              node_mon0[offs + i], node_mon1[offs + i],
                              node_oct0[offs + i], node_oct1[offs + i], eps2);
            #endif

#endif
          }
 	  __syncthreads();
          #ifdef INDSOFT
             node_eps[tid] = temp; //Restore original value in shmem
             __syncthreads();
          #endif
 	}
	__syncthreads();
#endif

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
	    sh_mass[tid] = posj.w;
	    sh_pos [tid] = (float3){posj.x, posj.y, posj.z};
	    sh_jid [tid] = body_list[n_direct + tid];

            #ifdef INDSOFT
              float   temp = sh_eps [tid];        //Store the value from Shmem into a register
              sh_eps [tid] = body_vel[body_list[n_direct + tid]].w;  //Load the softening
            #endif
	    __syncthreads();
#pragma unroll
	    for (int j = 0; j < DIMx; j++)
            {
              int selfGrav = (body_i != sh_jid[offs + j]);
// 	      if (body_i != sh_jid[offs + j]) //If statement replaced by multiplication
              {
		float ds2 =  1.0e10f;
                direCount++;
#ifdef DOGRAV
                #ifdef INDSOFT
                  acc_i = add_acc(acc_i, pos_i, sh_mass[offs + j], sh_pos[offs + j], sh_eps[offs + j], ds2, eps2, selfGrav);
                #else
                  acc_i = add_acc(acc_i, pos_i, sh_mass[offs + j], sh_pos[offs + j], ds2, eps2, selfGrav);
                #endif
#endif
                ds2 += (!selfGrav)*1.0e10f;
		if (ds2 < ds2_min) {
		  ngb     = sh_jid[offs + j];
		  ds2_min = ds2;
		}
	      }
	    }
            #ifdef INDSOFT
              //Restore the shmem value after all threads have used it
              __syncthreads();
              sh_eps[tid] = temp;
            #endif
	    __syncthreads();
	  }

	}
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
	  apprCount = -1; return acc_i;	 
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


  if(n_approx > 0)
  {
    #ifdef INDSOFT
      float temp = node_eps[tid];
    #endif

    if (tid < n_approx) {
      int address      = (approx[tid] << 1) + approx[tid];
      #ifndef TEXTURES
        float4 monopole  = multipole_data[address    ];
        float4 octopole0 = multipole_data[address + 1];
        float4 octopole1 = multipole_data[address + 2];
      #else
        float4 monopole  = tex1Dfetch(texMultipole, address);
        float4 octopole0 = tex1Dfetch(texMultipole, address + 1);
        float4 octopole1 = tex1Dfetch(texMultipole, address + 2);
      #endif

      node_mon0[tid] = monopole.w;
      node_mon1[tid] = (float3){monopole.x,  monopole.y,  monopole.z};
      node_oct0[tid] = (float3){octopole0.x, octopole0.y, octopole0.z};
      node_oct1[tid] = (float3){octopole1.x, octopole1.y, octopole1.z};

      #ifdef INDSOFT
        node_eps[tid]  = octopole0.w;
      #endif

    } else {
      //Set non-active memory locations to zero
      node_mon0[tid] = 0.0f;
      node_mon1[tid] = (float3){0.0f, 0.0f, 0.0f};
      node_oct0[tid] = (float3){0.0f, 0.0f, 0.0f};
      node_oct1[tid] = (float3){0.0f, 0.0f, 0.0f};

      #ifdef INDSOFT
        node_eps[tid]  = 0.01f;
      #endif
    }
    __syncthreads();
  #pragma unroll
    for (int i = 0; i < DIMx; i++)
    {
      apprCount++;
  #ifdef DOGRAV
      #ifdef INDSOFT
        acc_i = add_acc(acc_i, pos_i,
                        node_mon0[offs + i], node_mon1[offs + i],
                        node_oct0[offs + i], node_oct1[offs + i], node_eps[offs + i], eps2);

      #else
        acc_i = add_acc(acc_i, pos_i,
                        node_mon0[offs + i], node_mon1[offs + i],
                        node_oct0[offs + i], node_oct1[offs + i], eps2);
      #endif
  #endif
    }

    __syncthreads();
    #ifdef INDSOFT
      node_eps[tid]  = temp;
      __syncthreads();
    #endif
  } //if n_approx > 0

  if(n_direct > 0)
  {
    if (tid < n_direct) {

      float4 posj = body_pos[direct[tid]];

  //     float4 posj  = tex1Dfetch(texBody, direct[tid]);
      sh_mass[tid] = posj.w;
      sh_pos [tid] = (float3){posj.x, posj.y, posj.z};
      sh_jid [tid] = direct[tid];

      #ifdef INDSOFT
        sh_eps [tid] = body_vel[direct[tid]].w;
      #endif
    } else {
      sh_mass[tid] = 0.0f;
      sh_jid [tid] = -1;

      #ifdef INDSOFT
        sh_eps [tid] = 0.01f;
      #endif
    }

    __syncthreads();
  #pragma unroll
    for (int j = 0; j < DIMx; j++) {
      if ((sh_jid[offs + j] >= 0)) {
        int selfGrav = (body_i != sh_jid[offs + j]);
        float ds2 =  1.0e10f;
        direCount++;
  #ifdef DOGRAV
        #ifdef INDSOFT
          acc_i = add_acc(acc_i, pos_i, sh_mass[offs + j], sh_pos[offs + j], sh_eps[offs + j], ds2, eps2, selfGrav);
        #else
          acc_i = add_acc(acc_i, pos_i, sh_mass[offs + j], sh_pos[offs + j], ds2, eps2, selfGrav);
        #endif
  #endif
        ds2 += (!selfGrav)*1.0e10f;
        if (ds2 < ds2_min) {
          ngb     = sh_jid[offs + j];
          ds2_min = ds2;
        }
      }
    }
    __syncthreads();
  }

  /***
  **** --> reduce data between threads
  ***/
  sh_pot[tid] = acc_i.w;
  sh_acc[tid] = (float3){acc_i.x, acc_i.y, acc_i.z};
  sh_ds2[tid] = ds2_min;
  sh_ngb[tid] = ngb;
  __syncthreads();

  if (ty == 0) {
#pragma unroll
    for (int i = 1; i < DIMy; i++) {
      int idx = (i << DIM2x) + tx;
      acc_i.w += sh_pot[idx];
      acc_i.x += sh_acc[idx].x;
      acc_i.y += sh_acc[idx].y;
      acc_i.z += sh_acc[idx].z;
      if (sh_ds2[idx] < ds2_min) {
	ds2_min = sh_ds2[idx];
	ngb     = sh_ngb[idx];
      }
    }
  }
  __syncthreads();


  //Sum the interaction counters
  sh_ds2[tid] = direCount;
  sh_ngb[tid] = apprCount;

  __syncthreads();


  if (ty == 0) {
  #pragma unroll
    for (int i = 1; i < DIMy; i++){
        int idx = (i << DIM2x) + tx;
        direCount  += sh_ds2[idx];
        apprCount  += sh_ngb[idx];
    }
  }
  __syncthreads();

  return acc_i;
}

 
extern "C" __global__ void
__launch_bounds__(NTHREAD)
dev_approximate_gravity(const int n_active_groups,
                        int    n_bodies,
                        float eps2,
                        uint2 node_begend,                        
                        real4  *body_pos,                        
                        float4 *acc_out,
                        real4  *group_body_pos, 
                        real   *ds2_out,
                        int    *ngb_out,
                        int    *active_inout,
                        int2   *interactions,
                        uint2  *group_list,
                        real4  *multipole_data,
                        float4  *boxSizeInfo,                        
                        float4  *boxCenterInfo,                        
                        int     *MEM_BUF) {
//                                                    int     grpOffset){


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

//   volatile int *lmem = &MEM_BUF[blockIdx.x*LMEM_STACK_SIZE*blockDim.x + threadIdx.x*LMEM_STACK_SIZE];
//   int *lmem = &MEM_BUF[blockIdx.x*LMEM_STACK_SIZE*blockDim.x + threadIdx.x*LMEM_STACK_SIZE];
   int *lmem = &MEM_BUF[blockIdx.x*LMEM_STACK_SIZE*blockDim.x];


  /*********** set necessary thread constants **********/

//   real4 curGroupSize    = groupSizeInfo[active_groups[bid + grpOffset]];
//   int   groupData       = __float_as_int(curGroupSize.w);
//   uint body_i           =   groupData & CRITMASK;
//   uint nb_i             = ((groupData & INVCMASK) >> CRITBIT) + 1;
// 
//   real4 group_pos       = groupCenterInfo[active_groups[bid + grpOffset]];

//   if(tid == 0)
//   printf("[%f %f %f %f ] \n [%f %f %f %f ] %d %d \n",
//           curGroupSize.x, curGroupSize.y, curGroupSize.z, curGroupSize.w,
//           group_pos.x, group_pos.y, group_pos.z, group_pos.w, body_i, nb_i);
    uint2 grpInfo = group_list[bid];
    uint body_i = grpInfo.x;
    uint nb_i   = (grpInfo.y - grpInfo.x) + 1;

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
//     float group_eps = computeGroupSoftening(body_vel, body_i, shmem);
    float group_eps = 0;        //This is disabled for the moment

    int ngb_i;
    float ds2;

    float4 acc_i = {0.0f, 0.0f, 0.0f, 0.0f};

    int apprCount = 0;
    int direCount = 0;

    float4* body_vel = body_pos;

    acc_i = approximate_gravity<blockDim2, 0>( DIM2x, DIM2y, tid, tx, ty,
                                            body_i, pos_i, group_pos,
                                            eps2, node_begend,
                                            multipole_data, body_pos,
                                            shmem, lmem, ngb_i, apprCount, direCount, boxSizeInfo, curGroupSize, boxCenterInfo,
                                          group_eps, body_vel, ds2);
    if(apprCount < 0)
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
      apprCount = direCount = 0;
      acc_i = approximate_gravity<blockDim2, 8>( DIM2x, DIM2y, tid, tx, ty,
                                              body_i, pos_i, group_pos,
                                              eps2, node_begend,
                                              multipole_data, body_pos,
                                              shmem, lmem, ngb_i, apprCount, direCount, boxSizeInfo, curGroupSize, boxCenterInfo,
                                            group_eps, body_vel, ds2);

      lmem = &MEM_BUF[blockIdx.x*LMEM_STACK_SIZE*blockDim.x]; //Back to normal location

      if(threadIdx.x == 0)
      {
              atomicExch(&active_inout[n_bodies+1], 0); //Release the lock
      }
    }//end if apprCount < 0


    // changed by M.I.
    //acc_i.w = ds2;

    if (tid < nb_i) {
        acc_out     [body_i] = acc_i;
        ngb_out     [body_i] = ngb_i;
        ds2_out     [body_i] = ds2;
        active_inout[body_i] = 1;
        interactions[body_i].x = apprCount;
        interactions[body_i].y = direCount ;

      }
  }     //end while
}


