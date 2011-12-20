// //#include "/home/jbedorf/papers/GBPZ2010/codes/jb/build_tree/CUDA/support_kernels.cu"
#include "support_kernels.cu"

#include <stdio.h>

//////////////////////////////
//////////////////////////////
//////////////////////////////
#define LEVEL_MIN 3

extern "C" __global__ void boundaryReduction(const int n_particles,
                                            real4      *positions,
                                            float3     *output_min,
                                            float3     *output_max)
{
  const uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  const uint tid = threadIdx.x;
  //const uint idx = bid * blockDim.x + tid;

  volatile __shared__ float3 shmem[512];
  float3 r_min = (float3){+1e10f, +1e10f, +1e10f};
  float3 r_max = (float3){-1e10f, -1e10f, -1e10f};

  volatile float3 *sh_rmin = (float3*)&shmem [ 0];
  volatile float3 *sh_rmax = (float3*)&shmem[256];
  sh_rmin[tid].x = r_min.x; sh_rmin[tid].y = r_min.y; sh_rmin[tid].z = r_min.z;
  sh_rmax[tid].x = r_max.x; sh_rmax[tid].y = r_max.y; sh_rmax[tid].z = r_max.z;

  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  const int blockSize   = blockDim.x;
//   unsigned int tid      = threadIdx.x;
  unsigned int i        = blockIdx.x*(blockSize*2) + threadIdx.x;
  unsigned int gridSize = blockSize*2*gridDim.x;

  real4 pos;
  // we reduce multiple elements per thread.  The number is determined by the
  // number of active thread blocks (via gridSize).  More blocks will result
  // in a larger gridSize and therefore fewer elements per thread
  //based on reduce6 example
  while (i < n_particles) {
    if (i             < n_particles)
    {
      pos = positions[i];
      r_min.x = fminf(pos.x, r_min.x);
      r_min.y = fminf(pos.y, r_min.y);
      r_min.z = fminf(pos.z, r_min.z);
      r_max.x = fmaxf(pos.x, r_max.x);
      r_max.y = fmaxf(pos.y, r_max.y);
      r_max.z = fmaxf(pos.z, r_max.z);
    }
    if (i + blockSize < n_particles)
    {
      pos = positions[i + blockSize];
      r_min.x = fminf(pos.x, r_min.x);
      r_min.y = fminf(pos.y, r_min.y);
      r_min.z = fminf(pos.z, r_min.z);
      r_max.x = fmaxf(pos.x, r_max.x);
      r_max.y = fmaxf(pos.y, r_max.y);
      r_max.z = fmaxf(pos.z, r_max.z);
    }
    i += gridSize;
  }

  sh_rmin[tid].x = r_min.x; sh_rmin[tid].y = r_min.y; sh_rmin[tid].z = r_min.z;
  sh_rmax[tid].x = r_max.x; sh_rmax[tid].y = r_max.y; sh_rmax[tid].z = r_max.z;

  __syncthreads();
  // do reduction in shared mem  
  if(blockDim.x >= 512) if (tid < 256) {sh_MinMax(tid, tid + 256, &r_min, &r_max, sh_rmin, sh_rmax);} __syncthreads();
  if(blockDim.x >= 256) if (tid < 128) {sh_MinMax(tid, tid + 128, &r_min, &r_max, sh_rmin, sh_rmax);} __syncthreads();
  if(blockDim.x >= 128) if (tid < 64)  {sh_MinMax(tid, tid + 64,  &r_min, &r_max, sh_rmin, sh_rmax);} __syncthreads();

  if (tid < 32) 
  {
    sh_MinMax(tid, tid + 32, &r_min, &r_max, sh_rmin,sh_rmax);
    sh_MinMax(tid, tid + 16, &r_min, &r_max, sh_rmin,sh_rmax);
    sh_MinMax(tid, tid +  8, &r_min, &r_max, sh_rmin,sh_rmax);
    sh_MinMax(tid, tid +  4, &r_min, &r_max, sh_rmin,sh_rmax);
    sh_MinMax(tid, tid +  2, &r_min, &r_max, sh_rmin,sh_rmax);
    sh_MinMax(tid, tid +  1, &r_min, &r_max, sh_rmin,sh_rmax);
  }

  // write result for this block to global mem
  if (tid == 0)
  {
    //Compiler doesnt allow: volatile float3 = float3
    output_min[bid].x = sh_rmin[0].x; output_min[bid].y = sh_rmin[0].y; output_min[bid].z = sh_rmin[0].z;
    output_max[bid].x = sh_rmax[0].x; output_max[bid].y = sh_rmax[0].y; output_max[bid].z = sh_rmax[0].z;
  }

}


//Get the domain size, by taking into account the group size
extern "C" __global__ void boundaryReductionGroups(const int n_groups,
                                                   real4      *positions,
                                                   real4      *sizes,
                                                   float3     *output_min,
                                                   float3     *output_max)
{
  const uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  const uint tid = threadIdx.x;
  //const uint idx = bid * blockDim.x + tid;

  volatile __shared__ float3 shmem[512];
  float3 r_min = (float3){+1e10f, +1e10f, +1e10f};
  float3 r_max = (float3){-1e10f, -1e10f, -1e10f};

  volatile float3 *sh_rmin = (float3*)&shmem [ 0];
  volatile float3 *sh_rmax = (float3*)&shmem[256];
  sh_rmin[tid].x = r_min.x; sh_rmin[tid].y = r_min.y; sh_rmin[tid].z = r_min.z;
  sh_rmax[tid].x = r_max.x; sh_rmax[tid].y = r_max.y; sh_rmax[tid].z = r_max.z;

  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  const int blockSize   = blockDim.x;
//   unsigned int tid      = threadIdx.x;
  unsigned int i        = blockIdx.x*(blockSize*2) + threadIdx.x;
  unsigned int gridSize = blockSize*2*gridDim.x;

  real4 pos;
  real4 size;
  // we reduce multiple elements per thread.  The number is determined by the
  // number of active thread blocks (via gridSize).  More blocks will result
  // in a larger gridSize and therefore fewer elements per thread
  //based on reduce6 example
  while (i < n_groups) {
    if (i             < n_groups)
    {
      pos = positions[i];
      size = sizes[i];
      r_min.x = fminf(pos.x-size.x, r_min.x);
      r_min.y = fminf(pos.y-size.y, r_min.y);
      r_min.z = fminf(pos.z-size.z, r_min.z);
      r_max.x = fmaxf(pos.x+size.x, r_max.x);
      r_max.y = fmaxf(pos.y+size.y, r_max.y);
      r_max.z = fmaxf(pos.z+size.z, r_max.z);
    }
    if (i + blockSize < n_groups)
    {
      pos = positions[i + blockSize];
      size = sizes[i + blockSize];
      r_min.x = fminf(pos.x-size.x, r_min.x);
      r_min.y = fminf(pos.y-size.y, r_min.y);
      r_min.z = fminf(pos.z-size.z, r_min.z);
      r_max.x = fmaxf(pos.x+size.x, r_max.x);
      r_max.y = fmaxf(pos.y+size.y, r_max.y);
      r_max.z = fmaxf(pos.z+size.z, r_max.z);
    }
    i += gridSize;
  }

  sh_rmin[tid].x = r_min.x; sh_rmin[tid].y = r_min.y; sh_rmin[tid].z = r_min.z;
  sh_rmax[tid].x = r_max.x; sh_rmax[tid].y = r_max.y; sh_rmax[tid].z = r_max.z;

  __syncthreads();
  // do reduction in shared mem  
  if(blockDim.x >= 512) if (tid < 256) {sh_MinMax(tid, tid + 256, &r_min, &r_max, sh_rmin, sh_rmax);} __syncthreads();
  if(blockDim.x >= 256) if (tid < 128) {sh_MinMax(tid, tid + 128, &r_min, &r_max, sh_rmin, sh_rmax);} __syncthreads();
  if(blockDim.x >= 128) if (tid < 64)  {sh_MinMax(tid, tid + 64,  &r_min, &r_max, sh_rmin, sh_rmax);} __syncthreads();

  if (tid < 32) 
  {
    sh_MinMax(tid, tid + 32, &r_min, &r_max, sh_rmin,sh_rmax);
    sh_MinMax(tid, tid + 16, &r_min, &r_max, sh_rmin,sh_rmax);
    sh_MinMax(tid, tid +  8, &r_min, &r_max, sh_rmin,sh_rmax);
    sh_MinMax(tid, tid +  4, &r_min, &r_max, sh_rmin,sh_rmax);
    sh_MinMax(tid, tid +  2, &r_min, &r_max, sh_rmin,sh_rmax);
    sh_MinMax(tid, tid +  1, &r_min, &r_max, sh_rmin,sh_rmax);
  }

  // write result for this block to global mem
  if (tid == 0)
  {
    //Compiler doesnt allow: volatile float3 = float3
    output_min[bid].x = sh_rmin[0].x; output_min[bid].y = sh_rmin[0].y; output_min[bid].z = sh_rmin[0].z;
    output_max[bid].x = sh_rmax[0].x; output_max[bid].y = sh_rmax[0].y; output_max[bid].z = sh_rmax[0].z;
  }

}

//#define EXACT_KEY


#if 0
extern "C" __global__ void cl_build_key_list(uint2  *body_key,
                                            real4  *body_pos,
                                            int   n_bodies,
                                            real4  corner) {
  
  uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  uint tid = threadIdx.x;
  uint id  = bid * blockDim.x + tid;
  
  if (id > n_bodies) return;

  real4 pos = body_pos[id];
  int4 crd;
  
  real domain_fac = corner.w;
  
  #ifndef EXACT_KEY
    crd.x = (int)roundf(__fdividef((pos.x - corner.x), domain_fac));
    crd.y = (int)roundf(__fdividef((pos.y - corner.y) , domain_fac));
    crd.z = (int)roundf(__fdividef((pos.z - corner.z) , domain_fac));
  #else            
    crd.x = (int)((pos.x - corner.x) / domain_fac);
    crd.y = (int)((pos.y - corner.y) / domain_fac);
    crd.z = (int)((pos.z - corner.z) / domain_fac);
  #endif

//   crd.x = (int)((pos.x - corner.x) / domain_fac + 0.5);
//   crd.y = (int)((pos.y - corner.y) / domain_fac + 0.5);
//   crd.z = (int)((pos.z - corner.z) / domain_fac + 0.5);

   uint2 key = get_key(crd);


  if (id == n_bodies) key = (uint2){0xFFFFFFFF, 0xFFFFFFFF};

  body_key[id] = key;

}

#endif

extern "C" __global__ void cl_build_key_list(uint4  *body_key,
                                            real4  *body_pos,
                                            int   n_bodies,
                                            real4  corner) {
  
  uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  uint tid = threadIdx.x;
  uint id  = bid * blockDim.x + tid;
  
  if (id > n_bodies) return;

  real4 pos = body_pos[id];
  int4 crd;
  
  real domain_fac = corner.w;
  
  #ifndef EXACT_KEY
    crd.x = (int)roundf(__fdividef((pos.x - corner.x), domain_fac));
    crd.y = (int)roundf(__fdividef((pos.y - corner.y) , domain_fac));
    crd.z = (int)roundf(__fdividef((pos.z - corner.z) , domain_fac));
  #else            
    crd.x = (int)((pos.x - corner.x) / domain_fac);
    crd.y = (int)((pos.y - corner.y) / domain_fac);
    crd.z = (int)((pos.z - corner.z) / domain_fac);
  #endif

//   crd.x = (int)((pos.x - corner.x) / domain_fac + 0.5);
//   crd.y = (int)((pos.y - corner.y) / domain_fac + 0.5);
//   crd.z = (int)((pos.z - corner.z) / domain_fac + 0.5);

   uint4 key = get_key(crd);


  if (id == n_bodies) key = (uint4){0xFFFFFFFF, 0xFFFFFFFF, 0, 0};

  body_key[id] = key;

}

                
#if 1
extern "C" __global__ void build_phkey_list(uint4  *body_key,
                                            real4  *body_pos,
                                            int    n_bodies,
                                            real4  corner,
                                            uint   *reorder) {
 
  uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  uint tid = threadIdx.x;
  uint id  = bid * blockDim.x + tid;
  
  if (id > n_bodies) return;

  real4 pos = body_pos[id];
//   real4 pos = body_pos[reorder[id]];
  int4 crd;
  
  real domain_fac = corner.w;
  
  //Get the integer position, will be used for the key calculation
  #if 1
    crd.x = (int)roundf(__fdividef((pos.x - corner.x), domain_fac));
    crd.y = (int)roundf(__fdividef((pos.y - corner.y) , domain_fac));
    crd.z = (int)roundf(__fdividef((pos.z - corner.z) , domain_fac));
  #else
            
    crd.x = (int)((pos.x - corner.x) / domain_fac);
    crd.y = (int)((pos.y - corner.y) / domain_fac);
    crd.z = (int)((pos.z - corner.z) / domain_fac);
  #endif


  uint4 key_new = get_key(crd);

  if (id == n_bodies) key_new = (uint4){0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF};

  body_key[id] = key_new;

}
#else

extern "C" __global__ void build_phkey_list(uint2  *body_key,
                                            real4  *body_pos,
                                            int   n_bodies,
                                            real4  corner) {
 
  uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  uint tid = threadIdx.x;
  uint id  = bid * blockDim.x + tid;
  
  if (id > n_bodies) return;

  real4 pos = body_pos[id];
  int4 crd;
  
  real domain_fac = corner.w;
  
  //Get the integer position, will be used for the key calculation
  #ifndef EXACT_KEY
    crd.x = (int)roundf(__fdividef((pos.x - corner.x), domain_fac));
    crd.y = (int)roundf(__fdividef((pos.y - corner.y) , domain_fac));
    crd.z = (int)roundf(__fdividef((pos.z - corner.z) , domain_fac));
  #else            
    crd.x = (int)((pos.x - corner.x) / domain_fac);
    crd.y = (int)((pos.y - corner.y) / domain_fac);
    crd.z = (int)((pos.z - corner.z) / domain_fac);
  #endif

  
  const int bits = 18;
  int i,xi, yi, zi;
  int mask;
  long key;
    
  //0= 000, 1=001, 2=011, 3=010, 4=110, 5=111, 6=101, 7=100
  //000=0=0, 001=1=1, 011=3=2, 010=2=3, 110=6=4, 111=7=5, 101=5=6, 100=4=7
  const int C[8] = {0, 1, 7, 6, 3, 2, 4, 5};
    
  int temp;
    
  mask = 1 << (bits - 1);
  key  = 0;
    
    
    for(i = 0; i < bits; i++, mask >>= 1)
      {
        xi = (crd.x & mask) ? 1 : 0;
        yi = (crd.y & mask) ? 1 : 0;
        zi = (crd.z & mask) ? 1 : 0;        
               
        if(xi == 0 && yi == 0 && zi == 0)
        {
          temp = crd.z; crd.z = crd.y; crd.y = temp;
        }
        else  if(xi == 0 && yi == 0 && zi == 1)
        {
          temp = crd.x; crd.x = crd.y; crd.y = temp;
        }
        else  if(xi == 1 && yi == 0 && zi == 1)
        {
          temp = crd.x; crd.x = crd.y; crd.y = temp;
        }
        else  if(xi == 1 && yi == 0 && zi == 0)
        {
          crd.x = (crd.x) ^ (-1);
          crd.z = (crd.z) ^ (-1);
        }
        else  if(xi == 1 && yi == 1 && zi == 0)
        {
         crd.x = (crd.x) ^ (-1);
         crd.z = (crd.z) ^ (-1);
        }
        else  if(xi == 1 && yi == 1 && zi == 1)
        {
         temp = (crd.x) ^ (-1);         
         crd.x = (crd.y) ^ (-1);
         crd.y = temp;
        }
        else  if(xi == 0 && yi == 1 && zi == 1)
        {
         temp = (crd.x) ^ (-1);         
         crd.x = (crd.y) ^ (-1);
         crd.y = temp;
        }
        else
        {
         temp = (crd.z) ^ (-1);         
         crd.z = (crd.y) ^ (-1);
         crd.y = temp;          
        }
        
        int index = (xi << 2) + (yi << 1) + zi;
        key = (key << 3) + C[index];
      }

  uint2 key_new;
  key_new.x = key & 0xFFFFFFFF;
  key_new.y = (key >> 32) & 0xFFFFFFFF;

  if (id == n_bodies) key_new = (uint2){0xFFFFFFFF, 0xFFFFFFFF};

  body_key[id] = key_new;

}


#endif         
  

extern "C" __global__ void cl_build_valid_list(int n_bodies,
                                               int level,
                                               uint4  *body_key,
                                               uint *valid_list){
//                                                uint2 *test_key_data) {
  
  const uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  const uint tid = threadIdx.x;
  const uint id  = bid * blockDim.x + tid;
  const uint4 key_F = {0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF};
  const uint4 key_B = {0xFFFFFFF1, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF}; //A border, valid0 will become 1
  const uint4 key_I = {0xFFFFFFF2, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF}; //Ignore
  const uint4 key_E = {0xFFFFFFF3, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF}; //End
  const uint4 key_A = {0xFFFFFFF4, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF}; //Start and End
//   const uint2 key_TEST = {0x0, 0x0}; //Start and End

//TODO clean this if we dont use it
  
  if (id >= n_bodies) return;   // >=   since the last particle is extra boudnary particle
  
  uint4 mask = get_mask(level);
  mask.x = mask.x | ((uint)1 << 30) | ((uint)1 << 31);

  uint4 key_m;

  uint4 key_c    = body_key[id];


  uint4 key_p;
  if (id == 0)
  {
    key_m = key_F;
  }
  else
  {
    key_m = body_key[id-1];
  }

  if((id+1) <  n_bodies) //The last particle gets a different key to compare with
  {
    key_p = body_key[id+1];
  }
  else
    key_p = (uint4){0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF};


  int valid0 = 0;
  int valid1 = 0;

  if (cmp_uint4(key_c, key_A) == 0) {
    valid0 = 1; //Set a border
    valid1 = 1; //Set a border
  }
  else if (cmp_uint4(key_c, key_B) == 0) {
    valid0 = 1; //Set a border
  }
  else if (cmp_uint4(key_c, key_E) == 0) {
    valid1 = 1; //Set a border
  }
  else if (cmp_uint4(key_c, key_I) == 0) {
    //Do nothing
  }
  else if (cmp_uint4(key_c, key_F) != 0) {
    key_c.x = key_c.x & mask.x;
    key_c.y = key_c.y & mask.y;
    key_c.z = key_c.z & mask.z;

    key_p.x = key_p.x & mask.x;
    key_p.y = key_p.y & mask.y;
    key_p.z = key_p.z & mask.z;

    key_m.x = key_m.x & mask.x;
    key_m.y = key_m.y & mask.y;
    key_m.z = key_m.z & mask.z;

    valid0 = abs(cmp_uint4(key_c, key_m));
    valid1 = abs(cmp_uint4(key_c, key_p));
  }

   valid_list[id*2]   = id | ((valid0) << 31);
   valid_list[id*2+1] = id | ((valid1) << 31);

}


//////////////////////////////
//////////////////////////////
//////////////////////////////


extern "C" __global__ void cl_build_nodes(uint level,
                             uint  compact_list_len,
                             uint offset,
                             uint *compact_list,
//                              uint *compact_list_end,
                             uint4 *bodies_key,
                             uint4 *node_key,
                             uint  *n_children,
                             uint2 *node_bodies){
//                              uint  *testValidList) {

  uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  uint tid = threadIdx.x;
  uint id  = bid * blockDim.x + tid;

  if (id >= compact_list_len) return;

  uint  bi   = compact_list[id*2];
  uint  bj   = compact_list[id*2+1] + 1;
  

  uint4 key  = bodies_key[bi];
  uint4 mask = get_mask(level);
  key = (uint4){key.x & mask.x, key.y & mask.y, key.z & mask.z, 0}; 


  node_bodies[offset+id] = (uint2){bi | (level << BITLEVELS), bj};
  node_key   [offset+id] = key;
  n_children [offset+id] = 0;
  
  if ((int)level > (int)(LEVEL_MIN - 1)) 
    if (bj - bi <= NLEAF)                            //Leaf can only have NLEAF particles, if its more there will be a split
      for (int i = bi; i < bj; i++)
        bodies_key[i] = (uint4){0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF}; //sets the key to FF to indicate the body is used

}

//////////////////////////////
//////////////////////////////
//////////////////////////////


extern "C" __global__ void cl_link_tree(int n_nodes,
                            uint *n_children,
                            uint2 *node_bodies,
                            real4 *bodies_pos,
                            real4 corner,
                            uint2 *level_list,           //TODO could make this constant if it proves usefull
//                             uint* parent_id_list,
                            uint* valid_list,
                            uint4 *node_keys,
                            uint4 *bodies_key,
                            int   maxLevel) {

  const uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  const uint tid = threadIdx.x;
  uint id  = bid * blockDim.x + tid;
  
  if (id >= n_nodes) return;

  uint2 bij  = node_bodies[id];
  uint level = (bij.x &  LEVELMASK) >> BITLEVELS;
  uint bi    =  bij.x & ILEVELMASK;
  uint bj    =  bij.y;

  real4 pos  = bodies_pos[bi];
  int4 crd;
  real domain_fac = corner.w;

  #ifndef EXACT_KEY
    crd.x = (int)roundf(__fdividef((pos.x - corner.x), domain_fac));
    crd.y = (int)roundf(__fdividef((pos.y - corner.y) , domain_fac));
    crd.z = (int)roundf(__fdividef((pos.z - corner.z) , domain_fac));
  #else            
    crd.x = (int)((pos.x - corner.x) / domain_fac);
    crd.y = (int)((pos.y - corner.y) / domain_fac);
    crd.z = (int)((pos.z - corner.z) / domain_fac);
  #endif


  uint4 key = get_key(crd);


  /********* accumulate children *****/
  
  uint4 mask = get_mask(level - 1);
  key = (uint4){key.x & mask.x, key.y & mask.y,  key.z & mask.z, 0}; 

  uint2 cij;

  
  if(id > 0)
    cij = level_list[level-1];

  int ci;
  //Jeroen, modified this since we dont use textures in find_key,
  //the function will fail because out of bound memory access when id==0
  if(id > 0)
    ci = find_key(key, cij, node_keys);
  else
    ci = 0;

  //ci now points to the node that is the parent, was used in previous group method
//   parent_id_list[id] = ci;

  mask = get_imask(mask);
  key = (uint4) {key.x | mask.x, key.y | mask.y, key.z | mask.z, 0 };
  if (id > 0)   
    atomicAdd(&n_children[ci], (1 << 28));

  key = get_key(crd);
  mask = get_mask(level);
  key = (uint4) {key.x & mask.x, key.y & mask.y, key.z & mask.z, 0}; 

  /********* store the 1st child *****/

  cij = level_list[level+1];
  int cj = -1;

  cj = find_key(key, cij, node_keys);

  atomicOr(&n_children[id], cj); //Atomic since multiple threads can work on this

  uint valid =  id | (uint)(0 << 31); 

  
  if ((int)level > (int)(LEVEL_MIN - 1)) 
    if ((bj - bi) <= NLEAF)    
      valid = id | (uint)(1 << 31);   //Distinguish leaves and nodes

 valid_list[id] = valid; 
}

//Determines which level of node starts at which offset
extern "C" __global__ void build_level_list(const int n_nodes,
                                            const int n_leafs,
                                            uint *leafsIdxs,
                                            uint2 *node_bodies,                                      
                                            uint* valid_list)
{

  const uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  const uint tid = threadIdx.x;
  const uint id  = bid * blockDim.x + tid;
  
  if (id >= n_nodes-n_leafs) return;

  const int nodeID = leafsIdxs[id+n_leafs];   //Get the idx into the node_bodies array

  int level_c, level_m, level_p;


  uint2 bij   = node_bodies[leafsIdxs[id+n_leafs]];    //current non-leaf
  level_c     = (bij.x &  LEVELMASK) >> BITLEVELS;

  if((id+1) < (n_nodes-n_leafs))        //The last node gets a default lvl
  {
    bij         = node_bodies[leafsIdxs[id+1+n_leafs]]; //next non-leaf
    level_p     = (bij.x &  LEVELMASK) >> BITLEVELS;
  }
  else
    level_p     = MAXLEVELS+5;  //Last is always an end

  //Compare level with the node before and node after
  if(nodeID == 0)
  {
    level_m = -1;    
  }
  else
  {
    bij         = node_bodies[ leafsIdxs[id-1+n_leafs]]; //Get info of previous non-leaf node
    level_m     =  (bij.x &  LEVELMASK) >> BITLEVELS;   
  }

  int valid0 = 0;
  int valid1 = 0;

  valid0 = (level_c != level_m) << 31 | (id+n_leafs);
  valid1 = (level_c != level_p) << 31 | (id+n_leafs);

  valid_list[id*2]   = valid0;
  valid_list[id*2+1] = valid1;

} //end build_level_list


//Finds nodes/leafs that will become groups
//After executions valid_list contains the 
//valid nodes/leafs that form groups
extern "C" __global__ void build_group_list(int n_nodes,
                                            uint* parent_id_list,
                                            uint2 *node_bodies,
                                            uint* valid_list)
{

  uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  uint tid = threadIdx.x;
  uint id  = bid * blockDim.x + tid;
  
  if (id >= n_nodes) return;

  uint2 bij          = node_bodies[id];
  int ownChildren    =  bij.y - (bij.x & ILEVELMASK);
  

  bij  = node_bodies[parent_id_list[id]];  
  int parentChildren    =  bij.y - (bij.x & ILEVELMASK);


  //group if nchild <= NCRIT AND parent_nchild > NCRIT
  //if((ownChildren <= NCRIT) && (parentChildren > NCRIT))
  if((ownChildren <= NCRIT) && (parentChildren > NCRIT))
    valid_list[id] = id | (uint)(1 << 31);  //Group
  else
    valid_list[id] = id | (0 << 31);  //Not a group
}

//Finds nodes/leafs that will become groups
//After executions valid_list contains the 
//valid nodes/leafs that form groups
extern "C" __global__ void build_group_list2(int    n_particles,
                                             uint  *validList,
                                             real4  *bodies_pos,
                                             const float DIST)
{
  uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  uint tid = threadIdx.x;
  uint idx = bid * blockDim.x + tid;

  //TODO use shared mem ffor the positions
//since we use them multiple times?


  //Note that we do not include the final particle
  //Since there is no reason to check it
  if (idx >= n_particles) return;

  //Get the current 
  float4 curPos, nexPos, prevPos;

  curPos  =  bodies_pos[idx];

  //Have to check the first and last to prevent out of bound access
  if(idx+1 == n_particles)
    nexPos  =  curPos;
  else
    nexPos = bodies_pos[idx+1];

  if(idx == 0)
    prevPos = curPos;
  else
    prevPos =  bodies_pos[idx-1];

  //Compute geometrical distance
  float dsPlus = ((curPos.x-nexPos.x)*(curPos.x-nexPos.x)) + 
                 ((curPos.y-nexPos.y)*(curPos.y-nexPos.y)) + 
                 ((curPos.z-nexPos.z)*(curPos.z-nexPos.z));

  float dsMin = ((curPos.x-prevPos.x)*(curPos.x-prevPos.x)) + 
                ((curPos.y-prevPos.y)*(curPos.y-prevPos.y)) + 
                ((curPos.z-prevPos.z)*(curPos.z-prevPos.z));

  //Multiples of the preferred group size are _always_ valid
  int validStart = ((idx     % NCRIT) == 0);
  int validEnd   = (((idx+1) % NCRIT) == 0);

//   const int DIST = 1;
//   const float DIST = 44;

  //The extra possible split(s) if the distance between two particles is too large
  if(dsPlus > DIST) validEnd     = 1;
  if(dsMin  > DIST) validStart   = 1;
  
  //Last particle is always the end, n_particles dont have to be a multiple of NCRIT
  //so this is required
  if(idx+1 == n_particles) validEnd = 1;

  //Set valid
  validList[2*idx + 0] = (idx)   | (uint)(validStart << 31);
  validList[2*idx + 1] = (idx+1) | (uint)(validEnd   << 31);    
}

   
extern "C" __global__ void store_group_list(int    n_particles,
                                            int n_groups,
                                            uint  *validList,
                                            uint  *body2group_list,
                                            uint2 *group_list)
{
  uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  uint tid = threadIdx.x;
//   uint idx = bid * blockDim.x + tid;
  
  if(bid >= n_groups) return;

  int start = validList[2*bid];
  int end   = validList[2*bid+1];

  if((start + tid) < end)
  {
    body2group_list[start + tid] = bid;
  }

  if(tid == 0)
  {
     group_list[bid] = (uint2){start,end};
  }
}

extern "C" __global__ void expandLeafList(int n_leafs,
                                          uint  *leaf2NodeIdx,
                                          uint2 *node_bodies,
                                          uint  *leafPart2Body)
{
  uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  uint tid = threadIdx.x;
  uint idx = bid * blockDim.x + tid;


  if(bid >= n_leafs) return;

  uint2 bij  =  node_bodies[leaf2NodeIdx[bid]];
  uint bi    =  bij.x & ILEVELMASK;
  uint bj    =  bij.y;

  //Write the particle id at the correct location, only if we are 
  //below the end particle id
  if(bi+tid < bj)
  {
    leafPart2Body[idx] = idx;
  }
}
    

//Assign a grp id to each particle of that grp to 
//create particle -> group relation using the
//group -> particle relation
extern "C" __global__ void  build_body2group_list(const int n_groups,
                                  uint *group_list,
                                  uint2 *node_bodies,
                                  uint  *body2group_list)
{
  const int bid = gridDim.x   * blockIdx.y + blockIdx.x;
  const int tid = threadIdx.x;

  if (bid >= n_groups) return;

  const int nodeID = group_list[bid];

  uint2 bij          = node_bodies[nodeID];

  const uint firstChild = bij.x & ILEVELMASK;
  const uint nChildren  = bij.y - (bij.x & ILEVELMASK);

  int idx = firstChild+tid;

  //Save the group id for this particle
  if (tid < nChildren) 
    body2group_list[idx] = bid;  
}


