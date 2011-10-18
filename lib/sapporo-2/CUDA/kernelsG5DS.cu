/*

Sapporo 2 device kernels
GRAPE5

Version 1.0
CUDA DoubleSingle kernels

*/

// #include <stdio.h>

#include "../include/defines.h"

typedef float2 DS;  // double single;

struct DS4 {
  DS x, y, z, w;
};
struct DS2 {
  DS x, y;
};

__device__ DS to_DS(double a) {
  DS b;
  b.x = (float)a;
  b.y = (float)(a - b.x);
  return b;
}

__device__ double to_double(DS a) {
  double b;
  b = (double)((double)a.x + (double)a.y);
  return b;
}


// This function computes c = a + b.
__device__ DS dsadd(DS a, DS b) {
  // Compute dsa + dsb using Knuth's trick.
  float t1 = a.x + b.x;
  float e = t1 - a.x;
  float t2 = ((b.x - e) + (a.x - (t1 - e))) + a.y + b.y;
  
  // The result is t1 + t2, after normalization.
  DS c;
  c.x = e = t1 + t2;
  c.y = t2 - (e - t1);
  return c;
} // dsadd

// This function computes c = a + b.
__device__ DS dsadd(DS a, float b) {
  // Compute dsa + dsb using Knuth's trick.
  float t1 = a.x + b;
  float e = t1 - a.x;
  float t2 = ((b - e) + (a.x - (t1 - e))) + a.y;
  
  // The result is t1 + t2, after normalization.
  DS c;
  c.x = e = t1 + t2;
  c.y = t2 - (e - t1);
  return c;
} // dsadd


__device__ void body_body_interaction(float4 &acc_i,                                       
                                      DS4     pos_i,                                       
                                      DS4     pos_j, 
                                      float &EPS2) {

  float3 dr = {(pos_j.x.x - pos_i.x.x) + (pos_j.x.y - pos_i.x.y),
               (pos_j.y.x - pos_i.y.x) + (pos_j.y.y - pos_i.y.y),
               (pos_j.z.x - pos_i.z.x) + (pos_j.z.y - pos_i.z.y)};   // 3x3 = 9 FLOP

  float ds2 = ((dr.x*dr.x + (dr.y*dr.y)) + dr.z*dr.z);

  //EPS is in GRAPE5 always non-zero, if it is zero well then behaviour is undefined
  float inv_ds  = rsqrt(ds2 + EPS2);

/*if((ds2 + EPS2) == 0)
  inv_ds = 0;
*/

  float mass    = pos_j.w.x;
  float inv_ds2 = inv_ds*inv_ds;                         // 1 FLOP
  float inv_ds3 = mass * inv_ds*inv_ds2;                 // 2 FLOP
  
  // 3*4 + 3 = 15 FLOP
  acc_i.x = ((inv_ds3 * dr.x) + acc_i.x);
  acc_i.y = ((inv_ds3 * dr.y) + acc_i.y);
  acc_i.z = ((inv_ds3 * dr.z) + acc_i.z);
  
  acc_i.w = (mass * inv_ds  + acc_i.w);
}


/*
 *  blockDim.x = ni
 *  gridDim.x  = 16, 32, 64, 128, etc. 
 */ 

//TODO should make this depending on if we use Fermi or GT80/GT200
// #define ajc(i, j) (i + __mul24(blockDim.x,j))
#define ajc(i, j) (i + blockDim.x*j)
extern "C" __global__ void dev_evaluate_gravity(
                                     int        nj_total, 
                                     int        nj,
                                     int        offset,
                                     double4    *pos_j,                                      
                                     double4    *pos_i,
                                     double4    *acc_i,                                      
                                     double     EPS2_d) {

  //Divide by 32 since that is the size of DS4
//   __shared__ DS4 shared_pos[(256*(sizeof(DS4) + sizeof(float4))) / 32];

//TODO for some reason fixed shared memory does not work with starlab?
  __shared__ char shared_mem[NTHREADS*(sizeof(DS4))];
  DS4 *shared_pos = (DS4*)&shared_mem[0];
//   extern __shared__ DS4 shared_pos[];
  
  float EPS2 = (float)EPS2_d;

  DS4 pos;
  pos.x = to_DS(pos_i[threadIdx.x].x); pos.y = to_DS(pos_i[threadIdx.x].y);
  pos.z = to_DS(pos_i[threadIdx.x].z); pos.w = to_DS(pos_i[threadIdx.x].w);

  #define LARGEnum 1e10f
  float4 acc = {0.0f, 0.0f, 0.0f, 0.0f};
  
  int i    = blockIdx.x * (nj*blockDim.y) + nj*threadIdx.y;
  int tile = 0;

  while (i <  blockIdx.x * (nj*blockDim.y) + nj*threadIdx.y + nj) {
    
    if (i + threadIdx.x < nj_total) {
      shared_pos[ajc(threadIdx.x, threadIdx.y)].x = to_DS(pos_j[i + threadIdx.x].x);
      shared_pos[ajc(threadIdx.x, threadIdx.y)].y = to_DS(pos_j[i + threadIdx.x].y);
      shared_pos[ajc(threadIdx.x, threadIdx.y)].z = to_DS(pos_j[i + threadIdx.x].z);
      shared_pos[ajc(threadIdx.x, threadIdx.y)].w = to_DS(pos_j[i + threadIdx.x].w);
    } else {
      shared_pos[ajc(threadIdx.x, threadIdx.y)].x = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx.x, threadIdx.y)].y = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx.x, threadIdx.y)].z = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx.x, threadIdx.y)].w = (float2){0.0f,  -1.0f}; 
    }
    __syncthreads();
    
    int j  = min(nj - tile*blockDim.x, blockDim.x);
    int j1 = (j/16)*16;

    #pragma unroll 16
      for (int k = 0; k < j1; k++) {
          body_body_interaction(acc, pos, shared_pos[ajc(k, threadIdx.y)], EPS2);
        }
        
        for (int k = j1; k < j; k++) {
          body_body_interaction(acc, pos, shared_pos[ajc(k, threadIdx.y)], EPS2);
        }
  
    __syncthreads();

    i += blockDim.x;
    tile++;
  } //end while

  
  float4 *shared_acc = (float4*)&shared_pos[0];  
  acc.w = -acc.w;
  
  shared_acc[ajc(threadIdx.x, threadIdx.y)] = acc;
  __syncthreads();

  if (threadIdx.y == 0) {

    for (int i = 1; i < blockDim.y; i++) {
      float4 acc1 = shared_acc[ajc(threadIdx.x, i)];
      acc.x += acc1.x;
      acc.y += acc1.y;
      acc.z += acc1.z;
      acc.w += acc1.w;
    }
  }
  __syncthreads();
  
  if (threadIdx.y == 0) {
    //Convert results to double and write  
    acc_i[blockIdx.x * blockDim.x + threadIdx.x] = (double4){acc.x, acc.y, acc.z, acc.w};    
  }

}



/*
 *  blockDim.x = #of block in previous kernel
 *  gridDim.x  = ni
 */ 
extern "C" __global__ void dev_reduce_forces(
                                              double4 *acc_i) {
//  extern __shared__ float4 shared_acc[];
  __shared__ float4     shared_acc[NBLOCKS];
  
  int index = threadIdx.x * gridDim.x + blockIdx.x;

  //Convert the data to floats
  shared_acc[threadIdx.x] = (float4){acc_i[index].x, acc_i[index].y, acc_i[index].z, acc_i[index].w};
         
  __syncthreads();

  if (threadIdx.x == 0) {
    float4 acc0 = shared_acc[0];

    for (int i = 1; i < blockDim.x; i++) {
      acc0.x += shared_acc[i].x;
      acc0.y += shared_acc[i].y;
      acc0.z += shared_acc[i].z;
      acc0.w += shared_acc[i].w;
    }
    //Store the results
    acc_i[blockIdx.x] = (double4){acc0.x, acc0.y, acc0.z, acc0.w};
  }
  __syncthreads();  
}


/*
 * Function that moves the (changed) j-particles
 * to the correct address location.
*/
extern "C" __global__ void dev_copy_particles(int nj, int nj_max,                                                                                        
                                              double4   *pos_j, 
                                              double4   *pos_j_temp,
                                              int       *address_j) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  //Copy the changed particles
  if (index < nj)
  {   
     pos_j[address_j[index]] = pos_j_temp[index];
  }
}

/*

Function to predict the particles
DS version

*/
extern "C" __global__ void dev_predictor(int nj) {
//   int index = blockIdx.x * blockDim.x + threadIdx.x;
  
  //NOt used in GRAPE5
}

