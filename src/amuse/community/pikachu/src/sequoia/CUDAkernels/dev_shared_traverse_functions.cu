#include "node_specs.h"

/*

This file contains device functions that are used in the various methods that
traverse the tree. 

*/


/* 

Utility functions

*/

//Computes the storage index of the tree-traverse stack
template<int SHIFT>
__forceinline__ __device__ int ACCS(const int i)
{
  return (i & ((LMEM_STACK_SIZE << SHIFT) - 1))*blockDim.x + threadIdx.x;
}


#define BTEST(x) (-(int)(x))

__forceinline__ __device__ float4 get_float4(float4 const volatile &v)
{
  return make_float4(v.x, v.y, v.z, v.w);
}


/*

End utility functions

*/

/*

Compute the properties of the current group on the fly 

*/

__forceinline__ __device__ void computeGroupProps(real4 &group_pos, 
                                             real4 &curGroupSize,
                                             real4 pos_i,
                                             int* shmem)
{
    const int tid = threadIdx.x;
    volatile float3 *sh_rmin = (float3*)&shmem [ 0];                                                                                                                               
    volatile float3 *sh_rmax = (float3*)&sh_rmin[NTHREAD];                                                                                                                         
                                                                                                                                                                                  
    float3 r_min = (float3){+1e10f, +1e10f, +1e10f};                                                                                                                               
    float3 r_max = (float3){-1e10f, -1e10f, -1e10f};                                                                                                                               
                                                                                                                                                                                  
                                                                                                                                                                                  
    //Set the shared memory with the data                                                                                                                                          
  //  if (tid >= nb_i)                                                                                                                                                             
    {                                                                                                                                                                              
  //    sh_rmin[tid].x = r_min.x; sh_rmin[tid].y = r_min.y; sh_rmin[tid].z = r_min.z;                                                                                              
  //    sh_rmax[tid].x = r_max.x; sh_rmax[tid].y = r_max.y; sh_rmax[tid].z = r_max.z;                                                                                              
    }                                                                                                                                                                              
  //  else                                                                                                                                                                         
    {                                                                                                                                                                              
      sh_rmin[tid].x = r_min.x = pos_i.x; sh_rmin[tid].y = r_min.y = pos_i.y; sh_rmin[tid].z = r_min.z = pos_i.z;                                                                  
      sh_rmax[tid].x = r_max.x = pos_i.x; sh_rmax[tid].y = r_max.y = pos_i.y; sh_rmax[tid].z = r_max.z = pos_i.z;                                                                  
    }                                                                                                                                                                              
                                                                                                                                                                                  
    __syncthreads();                        

    if(blockDim.x >= 512) if (tid < 256) {sh_MinMax(tid, tid + 256, &r_min, &r_max, sh_rmin, sh_rmax); } __syncthreads();
    if(blockDim.x >= 256) if (tid < 128) {sh_MinMax(tid, tid + 128, &r_min, &r_max, sh_rmin, sh_rmax); } __syncthreads();
    if(blockDim.x >= 128) if (tid < 64)  {sh_MinMax(tid, tid + 64,  &r_min, &r_max, sh_rmin, sh_rmax); } __syncthreads();                                                         
    if(blockDim.x >= 64)  if (tid < 32)  {sh_MinMax(tid, tid + 32,  &r_min, &r_max, sh_rmin, sh_rmax); }                                                                             
    if(blockDim.x >= 32)  if (tid < 16)  {sh_MinMax(tid, tid + 16,  &r_min, &r_max, sh_rmin, sh_rmax); }                                                                             
                                                                                                                                                                                  
    if(tid < 8)                                                                                                                                                                    
    {                                                                                                                                                                              
          sh_MinMax(tid, tid +  8, &r_min, &r_max, sh_rmin, sh_rmax);                                                                                                              
          sh_MinMax(tid, tid +  4, &r_min, &r_max, sh_rmin, sh_rmax);                                                                                                              
          sh_MinMax(tid, tid +  2, &r_min, &r_max, sh_rmin, sh_rmax);                                                                                                              
          sh_MinMax(tid, tid +  1, &r_min, &r_max, sh_rmin, sh_rmax);                                                                                                              
    }                                                                                                                                                                              
                                                                                                                                                                                  
    __syncthreads();               

    r_min.x = sh_rmin[0].x;                                                                                                                                                        
    r_min.y = sh_rmin[0].y;                                                                                                                                                        
    r_min.z = sh_rmin[0].z;                                                                                                                                                        
    r_max.x = sh_rmax[0].x;                                                                                                                                                        
    r_max.y = sh_rmax[0].y;                                                                                                                                                        
    r_max.z = sh_rmax[0].z;                                                                                                                                                        
                                                                                                                                                                                  
    //Compute the group center and size                                                                                                                                          
    group_pos.x = 0.5*(r_min.x + r_max.x);                                                                                                                                       
    group_pos.y = 0.5*(r_min.y + r_max.y);                                                                                                                                       
    group_pos.z = 0.5*(r_min.z + r_max.z);                                                                                                                                       
                                                                                                                                                                                  
    float3 grpSize = (float3){fmaxf(fabs(group_pos.x-r_min.x), fabs(group_pos.x-r_max.x)),                                                                                       
                              fmaxf(fabs(group_pos.y-r_min.y), fabs(group_pos.y-r_max.y)),                                                                                       
                              fmaxf(fabs(group_pos.z-r_min.z), fabs(group_pos.z-r_max.z))};                                                                                      
                                                                                                                                                                                  
                                                                                                                                                                                  
    //Store the box size and opening criteria 
    curGroupSize.x = grpSize.x; 
    curGroupSize.y = grpSize.y;                                                                                                                                                  
    curGroupSize.z = grpSize.z;                                                                                                                                                  
    float l = max(grpSize.x, max(grpSize.y, grpSize.z));                                                                                                                         

    group_pos.w = l;          
}

/*

Compute the softening of the group

*/

__forceinline__ __device__ float computeGroupSoftening(real4 *body_vel, 
                                                      int body_i,                                                     
                                                      int* shmem)
{
    float group_eps;
    #ifdef INDSOFT
      eps2 = body_vel[body_i].w;
      group_eps = eps2;

      volatile float *reduc = (float*) &shmem[0];
      reduc[threadIdx.x] = eps2;

      //Find the maximum softening value for the particles in this group
      __syncthreads();
      // do reduction in shared mem
      if(blockDim.x >= 512) if (tid < 256) {reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 256]);} __syncthreads();
      if(blockDim.x >= 256) if (tid < 128) {reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 128]);} __syncthreads();
      if(blockDim.x >= 128) if (tid < 64)  {reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 64]);} __syncthreads();
      if(blockDim.x >= 64) if (tid < 32) { reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 32]);}
      if(blockDim.x >= 32) if (tid < 16) { reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 16]);}

      if(tid < 8)
      {
        reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 8]);
        reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 4]);
        reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 2]);
        reduc[threadIdx.x] = group_eps = fmaxf(group_eps, reduc[threadIdx.x + 1]);
      }
      __syncthreads();

      group_eps = reduc[0];
    #else
      group_eps  = 0;
    #endif
  
  return group_eps;
}


/*

Prefix sum functions

*/

template<class T>
 struct ADDOP {
  __device__ static inline T identity()           {return (T)(0);}
  __device__ static inline T apply(T a, T b)      {return (T)(a + b);};
  __device__ static inline T unapply(T a, T b)    {return (T)(a - b);};
  __device__ static inline T mask(bool flag, T b) {return (T)(-(int)(flag) & b);};
};

template<class OP, class T>
// __device__ T inclusive_scan_warp(volatile T *ptr, T mysum,  const unsigned int idx = threadIdx.x) {
__device__ __forceinline__ T inclusive_scan_warp(volatile T *ptr, T mysum,  const unsigned int idx ) {
  const unsigned int lane = idx & 31;

  if (lane >=  1) ptr[idx] = mysum = OP::apply(ptr[idx -  1], mysum);
  if (lane >=  2) ptr[idx] = mysum = OP::apply(ptr[idx -  2], mysum);
  if (lane >=  4) ptr[idx] = mysum = OP::apply(ptr[idx -  4], mysum);
  if (lane >=  8) ptr[idx] = mysum = OP::apply(ptr[idx -  8], mysum);
  if (lane >= 16) ptr[idx] = mysum = OP::apply(ptr[idx - 16], mysum);

  return ptr[idx];
}


__device__ __forceinline__ int inclusive_scan_warp(volatile int *ptr, int mysum, const unsigned int idx) {

  const unsigned int lane = idx & 31;

  if (lane >=  1) ptr[idx] = mysum = ptr[idx -  1]   + mysum;
  if (lane >=  2) ptr[idx] = mysum = ptr[idx -  2]   + mysum;
  if (lane >=  4) ptr[idx] = mysum = ptr[idx -  4]   + mysum;
  if (lane >=  8) ptr[idx] = mysum = ptr[idx -  8]   + mysum;
  if (lane >= 16) ptr[idx] = mysum = ptr[idx -  16]  + mysum;

  return ptr[idx];
}


template<class OP, class T>
__device__ __inline__ T inclusive_scan_block(volatile T *ptr, const T v0, const unsigned int idx) {
  const unsigned int lane   = idx & 31;
  const unsigned int warpid = idx >> 5;

  // step 0: Write the valume from the thread to the memory
  ptr[idx] = v0;
  T mysum = v0;
  __syncthreads();

  // step 1: Intra-warp scan in each warp
//   T val = inclusive_scan_warp<OP, T>(ptr, mysum, idx);
  T val = inclusive_scan_warp(ptr, mysum, idx);
  __syncthreads();

  // step 2: Collect per-warp particle results
  if (lane == 31) ptr[warpid] =  ptr[idx];
  __syncthreads();

  mysum =  ptr[idx];

  // step 3: Use 1st warp to scan per-warp results
  if (warpid == 0) inclusive_scan_warp<OP, T>(ptr,mysum, idx);
  __syncthreads();

  // step 4: Accumulate results from Steps 1 and 3;
  if (warpid > 0) val = OP::apply(ptr[warpid - 1], val);
  __syncthreads();

  // Step 5: Write and return the final result
  ptr[idx] = val;
  __syncthreads();

  return val; //ptr[blockDim.x - 1];
}



template<class OP, class T>
// __device__ T inclusive_scan_block(volatile T *ptr, const unsigned int idx = threadIdx.x) {
__device__ T inclusive_scan_block(volatile T *ptr, const unsigned int idx) {
  const unsigned int lane   = idx & 31;
  const unsigned int warpid = idx >> 5;

   T mysum = ptr[idx];
   __syncthreads();

  // step 1: Intra-warp scan in each warp
  T val = inclusive_scan_warp<OP, T>(ptr, mysum, idx);
  __syncthreads();

  // step 2: Collect per-warp particle results
  if (lane == 31) ptr[warpid] = ptr[idx];
  __syncthreads();

  mysum = ptr[idx];

  // step 3: Use 1st warp to scan per-warp results
  if (warpid == 0) inclusive_scan_warp<OP, T>(ptr,mysum, idx);
  __syncthreads();

  // step 4: Accumulate results from Steps 1 and 3;
  if (warpid > 0) val = OP::apply(ptr[warpid - 1], val);
  __syncthreads();

  // Step 5: Write and return the final result
  ptr[idx] = val;
  __syncthreads();

  return val; //ptr[blockDim.x - 1];
}


template<class OP, class T>
// __device__ T inclusive_scan_array(volatile T *ptr_global, const int N, const unsigned int idx = threadIdx.x) {
__device__ T inclusive_scan_array(volatile T *ptr_global, const int N, const unsigned int idx) {


  T y = OP::identity();
  volatile T *ptr = ptr_global;

  for (int p = 0; p < N; p += blockDim.x) {
    ptr = &ptr_global[p];
    inclusive_scan_block<OP, T>(ptr, idx);
    ptr[idx] = OP::apply(ptr[idx], y);
    __syncthreads();

    y = ptr[blockDim.x - 1];
    __syncthreads();
  }

  return y;

}


/* 

Opening criteria functions

*/


//1) Minimum distance opening criteria

#ifdef INDSOFT
  __device__ bool split_node_grav_md(float4 nodeCenter, float4 nodeSize, float4 groupCenter, float4 groupSize,
                                     float group_eps, float node_eps)
#else
  __device__ bool split_node_grav_md(float4 nodeCenter, float4 nodeSize, float4 groupCenter, float4 groupSize)
#endif
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

  #ifdef INDSOFT
    //Naar idee van Inti nu minder overbodige openingen
    if(ds2 <= ((group_eps + node_eps ) * (group_eps + node_eps) )) return true;
  #endif

  return (ds2 <= fabs(nodeCenter.w));
}


 // modified by M.I.
 //2) with rsearch
#ifdef INDSOFT
  __device__ bool split_node_grav_md_rsearch(float4 nodeCenter, 
					     float4 nodeSize, 
					     float4 groupCenter, 
					     float4 groupSize,
					     float group_eps, 
					     float node_eps, 
					     float rsearch_sq)
#else
  __device__ bool split_node_grav_md_rsearch(float4 nodeCenter, 
					     float4 nodeSize, 
					     float4 groupCenter, 
					     float4 groupSize, 
					     float rsearch_sq)
#endif
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

  #ifdef INDSOFT
    //Naar idee van Inti nu minder overbodige openingen
    if(ds2 <= ((group_eps + node_eps ) * (group_eps + node_eps) )) return true;
  #endif

    return ( (ds2 <= fabs(nodeCenter.w)) || ds2 < rsearch_sq);
}


// 3) Improved Barnes Hut criterium
#ifdef INDSOFT
__device__ bool split_node_grav_impbh(float4 nodeCOM, float4 groupCenter, float4 groupSize,
                                      float group_eps, float node_eps)
#else
__device__ bool split_node_grav_impbh(float4 nodeCOM, float4 groupCenter, float4 groupSize)
#endif
{
  //Compute the distance between the group and the cell
  float3 dr = {fabs(groupCenter.x - nodeCOM.x) - (groupSize.x),
               fabs(groupCenter.y - nodeCOM.y) - (groupSize.y),
               fabs(groupCenter.z - nodeCOM.z) - (groupSize.z)};

  dr.x += fabs(dr.x); dr.x *= 0.5f;
  dr.y += fabs(dr.y); dr.y *= 0.5f;
  dr.z += fabs(dr.z); dr.z *= 0.5f;

  //Distance squared, no need to do sqrt since opening criteria has been squared
  float ds2    = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

  #ifdef INDSOFT
    //Extra test
    if(ds2 <= ((group_eps + node_eps ) * (group_eps + node_eps) )) return true;
  #endif

  return (ds2 <= fabs(nodeCOM.w));
}

// by M.I.
// 4) Improved Barnes Hut criterium with rsearch
#ifdef INDSOFT
__device__ bool split_node_grav_impbh_rsearch(float4 nodeCOM, 
					      float4 nodeCenter, 
					      float4 nodeSize, 
					      float4 groupCenter, 
					      float4 groupSize,
					      float group_eps, 
					      float node_eps,
					      float rsearch_sq)
#else
__device__ bool split_node_grav_impbh_rsearch(float4 nodeCOM,
					      float4 nodeCenter, 
					      float4 nodeSize,  
					      float4 groupCenter, 
					      float4 groupSize,
					      float rsearch_sq)
#endif
{
  //Compute the distance between the group and the cell
  float3 dr_impbh = {fabs(groupCenter.x - nodeCOM.x) - (groupSize.x),
		     fabs(groupCenter.y - nodeCOM.y) - (groupSize.y),
		     fabs(groupCenter.z - nodeCOM.z) - (groupSize.z)};

  dr_impbh.x += fabs(dr_impbh.x); dr_impbh.x *= 0.5f;
  dr_impbh.y += fabs(dr_impbh.y); dr_impbh.y *= 0.5f;
  dr_impbh.z += fabs(dr_impbh.z); dr_impbh.z *= 0.5f;

  //Distance squared, no need to do sqrt since opening criteria has been squared
  float ds2_impbh = dr_impbh.x*dr_impbh.x + dr_impbh.y*dr_impbh.y + dr_impbh.z*dr_impbh.z;

  //Compute the distance between the group and the cell
  float3 dr_md = {fabs(groupCenter.x - nodeCenter.x) - (groupSize.x + nodeSize.x),
		  fabs(groupCenter.y - nodeCenter.y) - (groupSize.y + nodeSize.y),
		  fabs(groupCenter.z - nodeCenter.z) - (groupSize.z + nodeSize.z)};

  dr_md.x += fabs(dr_md.x); dr_md.x *= 0.5f;
  dr_md.y += fabs(dr_md.y); dr_md.y *= 0.5f;
  dr_md.z += fabs(dr_md.z); dr_md.z *= 0.5f;

  //Distance squared, no need to do sqrt since opening criteria has been squared
  float ds2_md    = dr_md.x*dr_md.x + dr_md.y*dr_md.y + dr_md.z*dr_md.z;

  #ifdef INDSOFT
    //Extra test
    if(ds2 <= ((group_eps + node_eps ) * (group_eps + node_eps) )) return true;
  #endif
    return ( (ds2_impbh <= fabs(nodeCOM.w)) || ds2_md < rsearch_sq);
}

