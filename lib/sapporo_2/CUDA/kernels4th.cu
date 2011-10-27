/*

Sapporo 2 device kernels

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


__device__ void body_body_interaction(float &ds_min,
                                      int &n_ngb,
                                      int *ngb_list,
                                      float4 &acc_i, 
                                      float4 &jrk_i,
                                      DS4     pos_i, 
                                      float4  vel_i,
                                      DS4     pos_j, 
                                      float4  vel_j,
                                      float &EPS2) {

  float3 dr = {(pos_j.x.x - pos_i.x.x) + (pos_j.x.y - pos_i.x.y),
               (pos_j.y.x - pos_i.y.x) + (pos_j.y.y - pos_i.y.y),
               (pos_j.z.x - pos_i.z.x) + (pos_j.z.y - pos_i.z.y)};   // 3x3 = 9 FLOP

  float ds2 = ((dr.x*dr.x + (dr.y*dr.y)) + dr.z*dr.z);

  bool ngb = true;

  if (ngb) {
    if (ds2 <= pos_i.w.x) {
      if (n_ngb < NGB_PB) {
        if(__float_as_int(pos_i.w.y) != __float_as_int(pos_j.w.y))      //Jeroen, is dit echt nodig?
          ngb_list[n_ngb++] = __float_as_int(pos_j.w.y);
      }
    }

    if (ds2 < ds_min*(__float_as_int(pos_i.w.y) != __float_as_int(pos_j.w.y))) {
      ds_min  = ds2;
      ngb_list[NGB_PB] = __float_as_int(pos_j.w.y);
//       ngb_list[NGB_PB] = (pos_j.w.y);
    }    
  }


  float inv_ds  = rsqrt(ds2 + EPS2) * (__float_as_int(pos_i.w.y) != __float_as_int(pos_j.w.y));
//   float inv_ds  = rsqrt(ds2 + EPS2) * (pos_i.w.y != pos_j.w.y);3

//   if(ds2 == 0)
//   {
// //     inv_ds = 0;
//      printf("ds 0, result: %d\t%d\n", __float_as_int(pos_i.w.y),__float_as_int(pos_j.w.y));
//   }

//TODO make sure the above trick still works on Fermi devices 
//and especially for double precision calculations
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

  float3 dv;    // 3 FLOP
  dv.x = vel_j.x - vel_i.x;
  dv.y = vel_j.y - vel_i.y;
  dv.z = vel_j.z - vel_i.z;

  float drdv = -3.0f * (inv_ds3*inv_ds2) * (((dr.x*dv.x) + dr.y*dv.y) + dr.z*dv.z);

//   jrk_i.x += 1 * (mass > 0);


  jrk_i.x = (jrk_i.x + inv_ds3 * dv.x) + drdv * dr.x;  
  jrk_i.y = (jrk_i.y + inv_ds3 * dv.y) + drdv * dr.y;
  jrk_i.z = (jrk_i.z + inv_ds3 * dv.z) + drdv * dr.z;

  // TOTAL 50 FLOP (or 60 FLOP if compared against GRAPE6)  
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
                                      double     EPS2_d,
                                      double4    *vel_j,
                                      int        *id_j,                                     
                                      double4    *vel_i,                                     
                                      double4    *jrk_i,
                                      int        *id_i,
                                      int        *ngb_list) {

  //Divide by 32 since that is the size of DS4
//   __shared__ DS4 shared_pos[(256*(sizeof(DS4) + sizeof(float4))) / 32];
//__shared__ double4 shared_pos[(sizeof(double4) + sizeof(double4) + sizeof(int)*2 + sizeof(double))];
 //__shared__ DS4 shared_pos[sizeof(DS4) + sizeof(float4)];
// __shared__ int shared_mem[256*(sizeof(double4) + sizeof(double4) + sizeof(int)*2 + sizeof(double))];

  
//  DS4 *shared_pos = (DS4*)&shared_mem[0];

//TODO for some reason fixed shared memory does not work with starlab?
   __shared__ char shared_mem[NTHREADS*(sizeof(DS4) + sizeof(float4))];
   DS4* shared_pos = (DS4*)&shared_mem[0];


//    extern __shared__ DS4 shared_pos[];
  float4 *shared_vel = (float4*)&shared_pos[blockDim.x*blockDim.y];

  int local_ngb_list[NGB_PB + 1];
  int n_ngb = 0;

  float EPS2 = (float)EPS2_d;

  DS4 pos;
  pos.x = to_DS(pos_i[threadIdx.x].x); pos.y = to_DS(pos_i[threadIdx.x].y);
  pos.z = to_DS(pos_i[threadIdx.x].z); pos.w = to_DS(pos_i[threadIdx.x].w);

  //Combine the particle id into the w part of the position
  pos.w.y = __int_as_float(id_i[threadIdx.x]);

  float4 vel = (float4){vel_i[threadIdx.x].x, vel_i[threadIdx.x].y, vel_i[threadIdx.x].z, vel_i[threadIdx.x].w};

  #define LARGEnum 1e10f
  float ds_min = LARGEnum;
  
  float4 acc = {0.0f, 0.0f, 0.0f, 0.0f};
  float4 jrk = {0.0f, 0.0f, 0.0f, 0.0f};

  
  int i    = blockIdx.x * (nj*blockDim.y) + nj*threadIdx.y;
  int tile = 0;

  while (i <  blockIdx.x * (nj*blockDim.y) + nj*threadIdx.y + nj) {
    
    if (i + threadIdx.x < nj_total) {
      shared_pos[ajc(threadIdx.x, threadIdx.y)].x = to_DS(pos_j[i + threadIdx.x].x);
      shared_pos[ajc(threadIdx.x, threadIdx.y)].y = to_DS(pos_j[i + threadIdx.x].y);
      shared_pos[ajc(threadIdx.x, threadIdx.y)].z = to_DS(pos_j[i + threadIdx.x].z);
      shared_pos[ajc(threadIdx.x, threadIdx.y)].w = to_DS(pos_j[i + threadIdx.x].w);
      //Combine the particle id into the w part of the position
      shared_pos[ajc(threadIdx.x, threadIdx.y)].w.y = __int_as_float (id_j[i + threadIdx.x]); 

      shared_vel[ajc(threadIdx.x, threadIdx.y)] = 
                (float4){vel_j[i + threadIdx.x].x, vel_j[i + threadIdx.x].y,
                         vel_j[i + threadIdx.x].z, vel_j[i + threadIdx.x].w};

    } else {
      shared_pos[ajc(threadIdx.x, threadIdx.y)].x = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx.x, threadIdx.y)].y = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx.x, threadIdx.y)].z = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx.x, threadIdx.y)].w = (float2){0.0f,  -1.0f}; 
      shared_vel[ajc(threadIdx.x, threadIdx.y)]   = (float4){0.0f, 0.0f, 0.0f, 0.0f};
    }
    __syncthreads();
    
    int j  = min(nj - tile*blockDim.x, blockDim.x);
    int j1 = (j/16)*16;

    #pragma unroll 16
      for (int k = 0; k < j1; k++) {
          body_body_interaction(ds_min, n_ngb, local_ngb_list,
                                    acc, jrk, pos, vel,
                                    shared_pos[ajc(k, threadIdx.y)], shared_vel[ajc(k, threadIdx.y)], EPS2);
        }
        
        for (int k = j1; k < j; k++) {
          body_body_interaction(ds_min, n_ngb, local_ngb_list,
                                    acc, jrk, pos, vel,
                                    shared_pos[ajc(k, threadIdx.y)], shared_vel[ajc(k, threadIdx.y)], EPS2);
        }
  

    __syncthreads();

    i += blockDim.x;
    tile++;
  } //end while

  
  float4 *shared_acc = (float4*)&shared_pos[0];
  float4 *shared_jrk = (float4*)&shared_acc[blockDim.x*blockDim.y];
  int    *shared_ngb = (int*   )&shared_jrk[blockDim.x*blockDim.y];
  int    *shared_ofs = (int*   )&shared_ngb[blockDim.x*blockDim.y];
  float  *shared_ds  = (float* )&shared_ofs[blockDim.x*blockDim.y];
  acc.w = -acc.w;
  jrk.w = __int_as_float(local_ngb_list[NGB_PB]);
//  jrk.w = local_ngb_list[NGB_PB];

  shared_acc[ajc(threadIdx.x, threadIdx.y)] = acc;
  shared_jrk[ajc(threadIdx.x, threadIdx.y)] = jrk;
  shared_ngb[ajc(threadIdx.x, threadIdx.y)] = n_ngb;
  shared_ofs[ajc(threadIdx.x, threadIdx.y)] = 0;
  shared_ds [ajc(threadIdx.x, threadIdx.y)] = ds_min;
  __syncthreads();

  if (threadIdx.y == 0) {

    for (int i = 1; i < blockDim.y; i++) {
      float4 acc1 = shared_acc[ajc(threadIdx.x, i)];
      float4 jrk1 = shared_jrk[ajc(threadIdx.x, i)];
      float  ds1  = shared_ds [ajc(threadIdx.x, i)];
      
      acc.x += acc1.x;
      acc.y += acc1.y;
      acc.z += acc1.z;
      acc.w += acc1.w;
      
      jrk.x += jrk1.x;
      jrk.y += jrk1.y;
      jrk.z += jrk1.z;
      
      if (ds1  < ds_min) {
        jrk.w   = jrk1.w;
        ds_min  = ds1;
      }

      shared_ofs[ajc(threadIdx.x, i)] = min(n_ngb + 1, NGB_PB);
      n_ngb += shared_ngb[ajc(threadIdx.x, i)];
    }
    n_ngb  = min(n_ngb, NGB_PB);
  }
  __syncthreads();
  
  if (threadIdx.y == 0) {
    //Convert results to double and write
    vel_i[offset  + blockIdx.x * blockDim.x + threadIdx.x].w = ds_min;
    acc_i[blockIdx.x * blockDim.x + threadIdx.x] = (double4){acc.x, acc.y, acc.z, acc.w};
    jrk_i[blockIdx.x * blockDim.x + threadIdx.x] = (double4){jrk.x, jrk.y, jrk.z, jrk.w};
  }



  offset  = threadIdx.x * NBLOCKS*NGB_PB + blockIdx.x * NGB_PB;
  offset += shared_ofs[ajc(threadIdx.x, threadIdx.y)];
  
  if (threadIdx.y == 0)
    ngb_list[offset++] = n_ngb;
  
  n_ngb = shared_ngb[ajc(threadIdx.x, threadIdx.y)];
  for (int i = 0; i < n_ngb; i++) 
    ngb_list[offset + i] = local_ngb_list[i];

}


extern "C" __global__ void dev_evaluate_gravity_new(                                     
                                     int        nj_total, 
                                     int        njPerBlock,
                                     int        offset,
                                     int        readOffset,
                                     double4    *pos_j, 
                                     double4    *vel_j,
                                     int        *id_j,
                                     double4    *pos_i,
                                     double4    *vel_i,
                                     double4    *acc_i, 
                                     double4    *jrk_i,
                                     int        *id_i,
                                     int        *ngb_list,
                                     double     EPS2_d,
                                     int        ni)
{
  extern __shared__ DS4 shared_pos[];
  float4 *shared_vel = (float4*)&shared_pos[blockDim.x*blockDim.y];

  //The i-particle for which we perform the force calculation
  int iReadLocation = threadIdx.x % ni;
  
  //The number of loops we have to perform 
  int iThread           = threadIdx.x   / ni;
  int nThreads          = blockDim.x    / ni;
  int nPerThread        = blockDim.x    / nThreads;
  int j0                = iThread       * nPerThread; //Start location of the loop
  int j1                = (iThread+1)   * nPerThread; //End location of the loop
  if(iThread+1 >= nThreads)
    j1 = blockDim.x;    //End location of the last iThread block is different!

//   j1 = min(blockIdx.x * (njPerBlock*blockDim.y) + njPerBlock*threadIdx.y + njPerBlock, j1);


  int local_ngb_list[NGB_PB + 1];
  int n_ngb = 0;

  float EPS2 = (float)EPS2_d;

  DS4 pos;
  pos.x = to_DS(pos_i[iReadLocation].x); pos.y = to_DS(pos_i[iReadLocation].y);
  pos.z = to_DS(pos_i[iReadLocation].z); pos.w = to_DS(pos_i[iReadLocation].w);

  //Combine the particle id into the w part of the position
  pos.w.y = __int_as_float(id_i[iReadLocation]);

  float4 vel = (float4){vel_i[iReadLocation].x, vel_i[iReadLocation].y, vel_i[iReadLocation].z, vel_i[iReadLocation].w};

  #define LARGEnum 1e10f
  float ds_min = LARGEnum;
  
  float4 acc = {0.0f, 0.0f, 0.0f, 0.0f};
  float4 jrk = {0.0f, 0.0f, 0.0f, 0.0f};


//   int interAct = 0;
//   printf("[%d,%d] %d\t%d\t%d\t%d \n", threadIdx.x, blockIdx.x, iThread, nThreads, j0, j1);

//   if(threadIdx.x == 0 && blockIdx.x == 0)    printf("njPerBlock: %d \n", njPerBlock);

  //Now to calculate which data to use, use the sapporov1 method first
//   int i = blockIdx.x * (njPerBlock*blockDim.y) + njPerBlock*threadIdx.y;
    int i = blockIdx.x * njPerBlock;
//   int tile = 0;

//    if(threadIdx.x != 0 && blockIdx.x == 0 && ni == 95)
//    printf("[%d\t%d, %d], njPerBlock: %d i: %d test: %d\tnj_total: %d\tj: %d\tj1: %d \t %d\n", threadIdx.x, threadIdx.y, blockIdx.x, njPerBlock, i,
//                 blockIdx.x * (njPerBlock*blockDim.y) + njPerBlock*threadIdx.y + njPerBlock, nj_total, j0, j1, ajc(threadIdx.x, threadIdx.y));

//    if(threadIdx.x != 0 && blockIdx.x == 0 && ni == 95)
//         printf("[%d,%d], njPerBlock: %d i: %d test: %d\n", threadIdx.x, blockIdx.x, njPerBlock, i,
//                blockIdx.x * (njPerBlock*blockDim.y) + njPerBlock*threadIdx.y + njPerBlock);
  while (i <  blockIdx.x * njPerBlock + njPerBlock) { 
    
    
    if (i + threadIdx.x < nj_total) {
      shared_pos[threadIdx.x].x = to_DS(pos_j[readOffset + i + threadIdx.x].x);
      shared_pos[threadIdx.x].y = to_DS(pos_j[readOffset + i + threadIdx.x].y);
      shared_pos[threadIdx.x].z = to_DS(pos_j[readOffset + i + threadIdx.x].z);
      shared_pos[threadIdx.x].w = to_DS(pos_j[readOffset + i + threadIdx.x].w);
      //Combine the particle id into the w part of the position
      shared_pos[threadIdx.x].w.y = __int_as_float (id_j[readOffset + i + threadIdx.x]); 

      shared_vel[threadIdx.x] = 
                (float4){vel_j[readOffset + i + threadIdx.x].x, vel_j[readOffset + i + threadIdx.x].y,
                         vel_j[readOffset + i + threadIdx.x].z, vel_j[readOffset + i + threadIdx.x].w};

    } else {


//       if(blockIdx.x == 0)    printf("Read a null value using thread: %d \n", threadIdx.x);

      shared_pos[ajc(threadIdx.x, threadIdx.y)].x = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx.x, threadIdx.y)].y = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx.x, threadIdx.y)].z = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx.x, threadIdx.y)].w = (float2){0.0f,  -1.0f}; 
      shared_vel[ajc(threadIdx.x, threadIdx.y)]   = (float4){0.0f, 0.0f, 0.0f, 0.0f};
    }


  //j1 cant be larger than the number of items that have been read into 
  //shared memory j1-j0
  //Todo figure out how to change it so it only loops the correct values
  j1 = min((blockIdx.x * njPerBlock + njPerBlock)-i, j1);
  __syncthreads();  


//   if(blockIdx.x == 1)
//   if(threadIdx.x == 0)
//    if(threadIdx.x != 0 && blockIdx.x == 0 && ni == 95)
//    printf("[%d\t%d, %d], njPerBlock: %d i: %d test: %d\tnj_total: %d\tj: %d\tj1: %d\n", threadIdx.x, threadIdx.y, blockIdx.x, njPerBlock, i,
//                 blockIdx.x * (njPerBlock*blockDim.y) + njPerBlock*threadIdx.y + njPerBlock, nj_total, j0, j1);

  



    #pragma unroll 16
    for (int k = j0; k < j1; k++) {
              body_body_interaction(ds_min, n_ngb, local_ngb_list,
                                    acc, jrk, pos, vel,
                                    shared_pos[k], shared_vel[k], EPS2);
        }


    __syncthreads();

    i += blockDim.x;
//     tile++;
  }

  
  //Reduce the data 
  float4 *shared_acc = (float4*)&shared_pos[0];
  float4 *shared_jrk = (float4*)&shared_acc[blockDim.x*blockDim.y];
  int    *shared_ngb = (int*   )&shared_jrk[blockDim.x*blockDim.y];
  int    *shared_ofs = (int*   )&shared_ngb[blockDim.x*blockDim.y];
  float  *shared_ds  = (float* )&shared_ofs[blockDim.x*blockDim.y];
  acc.w = -acc.w;
  jrk.w = __int_as_float(local_ngb_list[NGB_PB]);
//  jrk.w = local_ngb_list[NGB_PB];

  shared_acc[threadIdx.x] = acc;
  shared_jrk[threadIdx.x] = jrk;
  shared_ngb[threadIdx.x] = n_ngb;
  shared_ofs[threadIdx.x] = 0;
  shared_ds [threadIdx.x] = ds_min;
  __syncthreads();


  for (int i = ni + threadIdx.x; i < ni*nThreads; i += ni) {

//       float4 acc1 = shared_acc[ajc(threadIdx.x, i)];
//       float4 jrk1 = shared_jrk[ajc(threadIdx.x, i)];
//       float  ds1  = shared_ds [ajc(threadIdx.x, i)];
      float4 acc1 = shared_acc[i];
      float4 jrk1 = shared_jrk[i];
      float  ds1  = shared_ds [i];
      
      acc.x += acc1.x;
      acc.y += acc1.y;
      acc.z += acc1.z;
      acc.w += acc1.w;
      
      jrk.x += jrk1.x;
      jrk.y += jrk1.y;
      jrk.z += jrk1.z;
      
      if (ds1  < ds_min) {
        jrk.w   = jrk1.w;
        ds_min  = ds1;
      }

//       shared_ofs[ajc(threadIdx.x, i)] = min(n_ngb + 1, NGB_PB);
//       n_ngb += shared_ngb[ajc(threadIdx.x, i)];
      shared_ofs[i] = min(n_ngb + 1, NGB_PB);
      n_ngb += shared_ngb[i];
    
  }//end reduction

 __syncthreads();

  
  if (threadIdx.x < ni) {
    //Convert results to double and write
    vel_i[offset  + blockIdx.x * blockDim.x + threadIdx.x].w = ds_min;
//     acc_i[blockIdx.x * blockDim.x + threadIdx.x] = (double4){acc.x, acc.y, acc.z, acc.w};
    acc_i[blockIdx.x * ni + threadIdx.x] = (double4){acc.x, acc.y, acc.z, acc.w};
    jrk_i[blockIdx.x * ni + threadIdx.x] = (double4){jrk.x, jrk.y, jrk.z, jrk.w};
  }
 

  offset  = threadIdx.x * NBLOCKS*NGB_PB + blockIdx.x * NGB_PB;
  offset += shared_ofs[threadIdx.x];
  
  if (threadIdx.x < ni)
    ngb_list[offset++] = n_ngb;
  
  n_ngb = shared_ngb[threadIdx.x];
  for (int i = 0; i < n_ngb; i++) 
    ngb_list[offset + i] = local_ngb_list[i];

}









__device__ void body_body_interaction_nongb(float &ds_min,                                     
                                      float4 &acc_i, 
                                      float4 &jrk_i,
                                      DS4     pos_i, 
                                      float4  vel_i,
                                      DS4     pos_j, 
                                      float4  vel_j,
                                      float &EPS2) {

  float3 dr = {(pos_j.x.x - pos_i.x.x) + (pos_j.x.y - pos_i.x.y),
               (pos_j.y.x - pos_i.y.x) + (pos_j.y.y - pos_i.y.y),
               (pos_j.z.x - pos_i.z.x) + (pos_j.z.y - pos_i.z.y)};   // 3x3 = 9 FLOP

  float ds2 = ((dr.x*dr.x + (dr.y*dr.y)) + dr.z*dr.z);

  float inv_ds  = rsqrt(ds2 + EPS2) * (__float_as_int(pos_i.w.y) != __float_as_int(pos_j.w.y));
//   float inv_ds  = rsqrt(ds2 + EPS2) * (pos_i.w.y != pos_j.w.y);3

//TODO make sure the above trick still works on Fermi devices 
//and especially for double precision calculations
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

  float3 dv;    // 3 FLOP
  dv.x = vel_j.x - vel_i.x;
  dv.y = vel_j.y - vel_i.y;
  dv.z = vel_j.z - vel_i.z;

  float drdv = -3.0f * (inv_ds3*inv_ds2) * (((dr.x*dv.x) + dr.y*dv.y) + dr.z*dv.z);


  jrk_i.x = (jrk_i.x + inv_ds3 * dv.x) + drdv * dr.x;
  jrk_i.y = (jrk_i.y + inv_ds3 * dv.y) + drdv * dr.y;
  jrk_i.z = (jrk_i.z + inv_ds3 * dv.z) + drdv * dr.z;

  // TOTAL 50 FLOP (or 60 FLOP if compared against GRAPE6)  
}



//TODO should make this depending on if we use Fermi or GT80/GT200
//#define ajc(i, j) (i + __mul24(blockDim.x,j))
// #define ajc(i, j) (i + __mul24(blockDim.x,j))
#define ajc(i, j) (i + blockDim.x*j)
extern "C" __global__ void dev_evaluate_gravity_nongb(
                                     int        nj_total, 
                                     int        nj,
                                     int        offset,
                                     int        readOffset,
                                     double4    *pos_j, 
                                     double4    *vel_j,
                                     int        *id_j,
                                     double4    *pos_i,
                                     double4    *vel_i,
                                     double4    *acc_i, 
                                     double4    *jrk_i,
                                     int        *id_i,
                                     int        *ngb_list,
                                     double     EPS2_d) {
  extern __shared__ DS4 shared_pos[];
  float4 *shared_vel = (float4*)&shared_pos[blockDim.x*blockDim.y];

  float EPS2 = (float)EPS2_d;

  DS4 pos;
  pos.x = to_DS(pos_i[threadIdx.x].x); pos.y = to_DS(pos_i[threadIdx.x].y);
  pos.z = to_DS(pos_i[threadIdx.x].z); pos.w = to_DS(pos_i[threadIdx.x].w);

  //Combine the particle id into the w part of the position
  pos.w.y = __int_as_float(id_i[threadIdx.x]);

  float4 vel = (float4){vel_i[threadIdx.x].x, vel_i[threadIdx.x].y, vel_i[threadIdx.x].z, vel_i[threadIdx.x].w};

  #define LARGEnum 1e10f
  float ds_min = LARGEnum;
  
  float4 acc = {0.0f, 0.0f, 0.0f, 0.0f};
  float4 jrk = {0.0f, 0.0f, 0.0f, 0.0f};

  int i = blockIdx.x * (nj*blockDim.y) + nj*threadIdx.y;
  int tile = 0;


  while (i <  blockIdx.x * (nj*blockDim.y) + nj*threadIdx.y + nj) { 
    
    
    if (i + threadIdx.x < nj_total) {
      shared_pos[ajc(threadIdx.x, threadIdx.y)].x = to_DS(pos_j[readOffset + i + threadIdx.x].x);
      shared_pos[ajc(threadIdx.x, threadIdx.y)].y = to_DS(pos_j[readOffset + i + threadIdx.x].y);
      shared_pos[ajc(threadIdx.x, threadIdx.y)].z = to_DS(pos_j[readOffset + i + threadIdx.x].z);
      shared_pos[ajc(threadIdx.x, threadIdx.y)].w = to_DS(pos_j[readOffset + i + threadIdx.x].w);
      //Combine the particle id into the w part of the position
      shared_pos[ajc(threadIdx.x, threadIdx.y)].w.y = __int_as_float (id_j[readOffset + i + threadIdx.x]); 

      shared_vel[ajc(threadIdx.x, threadIdx.y)] = 
                (float4){vel_j[readOffset + i + threadIdx.x].x, vel_j[readOffset + i + threadIdx.x].y,
                         vel_j[readOffset + i + threadIdx.x].z, vel_j[readOffset + i + threadIdx.x].w};

    } else {
      shared_pos[ajc(threadIdx.x, threadIdx.y)].x = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx.x, threadIdx.y)].y = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx.x, threadIdx.y)].z = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx.x, threadIdx.y)].w = (float2){0.0f,  -1.0f}; 
      shared_vel[ajc(threadIdx.x, threadIdx.y)]   = (float4){0.0f, 0.0f, 0.0f, 0.0f};
    }
    __syncthreads();
    
    int j  = min(nj - tile*blockDim.x, blockDim.x);
    int j1 = (j/16)*16;

    #pragma unroll 16
       for (int k = 0; k < j1; k++) {
          body_body_interaction_nongb(ds_min, acc, jrk, pos, vel,
                                    shared_pos[ajc(k, threadIdx.y)], shared_vel[ajc(k, threadIdx.y)], EPS2);
        }
        
        for (int k = j1; k < j; k++) {
          body_body_interaction_nongb(ds_min, acc, jrk, pos, vel,
                                    shared_pos[ajc(k, threadIdx.y)], shared_vel[ajc(k, threadIdx.y)], EPS2);
        }
  

    __syncthreads();

    i += blockDim.x;
    tile++;
  }

  
  float4 *shared_acc = (float4*)&shared_pos[0];
  float4 *shared_jrk = (float4*)&shared_acc[blockDim.x*blockDim.y];
  int    *shared_ngb = (int*   )&shared_jrk[blockDim.x*blockDim.y];
  int    *shared_ofs = (int*   )&shared_ngb[blockDim.x*blockDim.y];
  float  *shared_ds  = (float* )&shared_ofs[blockDim.x*blockDim.y];
  acc.w = -acc.w;
//   jrk.w = __int_as_float(local_ngb_list[NGB_PB]);
//  jrk.w = local_ngb_list[NGB_PB];

  shared_acc[ajc(threadIdx.x, threadIdx.y)] = acc;
  shared_jrk[ajc(threadIdx.x, threadIdx.y)] = jrk;
//   shared_ngb[ajc(threadIdx.x, threadIdx.y)] = n_ngb;
  shared_ofs[ajc(threadIdx.x, threadIdx.y)] = 0;
  shared_ds [ajc(threadIdx.x, threadIdx.y)] = ds_min;
  __syncthreads();

  if (threadIdx.y == 0) {

    for (int i = 1; i < blockDim.y; i++) {
      float4 acc1 = shared_acc[ajc(threadIdx.x, i)];
      float4 jrk1 = shared_jrk[ajc(threadIdx.x, i)];
      float  ds1  = shared_ds [ajc(threadIdx.x, i)];
      
      acc.x += acc1.x;
      acc.y += acc1.y;
      acc.z += acc1.z;
      acc.w += acc1.w;
      
      jrk.x += jrk1.x;
      jrk.y += jrk1.y;
      jrk.z += jrk1.z;
      
      if (ds1  < ds_min) {
        jrk.w   = jrk1.w;
        ds_min  = ds1;
      }    
    }
  }
  __syncthreads();
  
  if (threadIdx.y == 0) {
    //Convert results to double and write
    vel_i[offset  + blockIdx.x * blockDim.x + threadIdx.x].w = ds_min;
    acc_i[blockIdx.x * blockDim.x + threadIdx.x] = (double4){acc.x, acc.y, acc.z, acc.w};
    jrk_i[blockIdx.x * blockDim.x + threadIdx.x] = (double4){jrk.x, jrk.y, jrk.z, jrk.w};
  }
}


extern "C" __global__ void dev_evaluate_gravity_new_nongb(                                     
                                     int        nj_total, 
                                     int        njPerBlock,
                                     int        offset,
                                     int        readOffset,
                                     double4    *pos_j, 
                                     double4    *vel_j,
                                     int        *id_j,
                                     double4    *pos_i,
                                     double4    *vel_i,
                                     double4    *acc_i, 
                                     double4    *jrk_i,
                                     int        *id_i,
                                     int        *ngb_list,
                                     double     EPS2_d,
                                     int        ni)
{
  extern __shared__ DS4 shared_pos[];
  float4 *shared_vel = (float4*)&shared_pos[blockDim.x*blockDim.y];

  //The i-particle for which we perform the force calculation
  int iReadLocation = threadIdx.x % ni;
  
  //The number of loops we have to perform 
  int iThread           = threadIdx.x   / ni;
  int nThreads          = blockDim.x    / ni;
  int nPerThread        = blockDim.x    / nThreads;
  int j0                = iThread       * nPerThread; //Start location of the loop
  int j1                = (iThread+1)   * nPerThread; //End location of the loop
  if(iThread+1 >= nThreads)
    j1 = blockDim.x;    //End location of the last iThread block is different!

  float EPS2 = (float)EPS2_d;

  DS4 pos;
  pos.x = to_DS(pos_i[iReadLocation].x); pos.y = to_DS(pos_i[iReadLocation].y);
  pos.z = to_DS(pos_i[iReadLocation].z); pos.w = to_DS(pos_i[iReadLocation].w);

  //Combine the particle id into the w part of the position
  pos.w.y = __int_as_float(id_i[iReadLocation]);

  float4 vel = (float4){vel_i[iReadLocation].x, vel_i[iReadLocation].y, vel_i[iReadLocation].z, vel_i[iReadLocation].w};

  #define LARGEnum 1e10f
  float ds_min = LARGEnum;
  
  float4 acc = {0.0f, 0.0f, 0.0f, 0.0f};
  float4 jrk = {0.0f, 0.0f, 0.0f, 0.0f};

  //Now to calculate which data to use, use the sapporov1 method first
//   int i = blockIdx.x * (njPerBlock*blockDim.y) + njPerBlock*threadIdx.y;
  int i = blockIdx.x * njPerBlock;

  while (i <  blockIdx.x * njPerBlock + njPerBlock) { 
    
    
    if (i + threadIdx.x < nj_total) {
      shared_pos[threadIdx.x].x = to_DS(pos_j[readOffset + i + threadIdx.x].x);
      shared_pos[threadIdx.x].y = to_DS(pos_j[readOffset + i + threadIdx.x].y);
      shared_pos[threadIdx.x].z = to_DS(pos_j[readOffset + i + threadIdx.x].z);
      shared_pos[threadIdx.x].w = to_DS(pos_j[readOffset + i + threadIdx.x].w);
      //Combine the particle id into the w part of the position
      shared_pos[threadIdx.x].w.y = __int_as_float (id_j[readOffset + i + threadIdx.x]); 

      shared_vel[threadIdx.x] = 
                (float4){vel_j[readOffset + i + threadIdx.x].x, vel_j[readOffset + i + threadIdx.x].y,
                         vel_j[readOffset + i + threadIdx.x].z, vel_j[readOffset + i + threadIdx.x].w};

    } else {
      shared_pos[ajc(threadIdx.x, threadIdx.y)].x = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx.x, threadIdx.y)].y = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx.x, threadIdx.y)].z = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx.x, threadIdx.y)].w = (float2){0.0f,  -1.0f}; 
      shared_vel[ajc(threadIdx.x, threadIdx.y)]   = (float4){0.0f, 0.0f, 0.0f, 0.0f};
    }


  //j1 cant be larger than the number of items that have been read into 
  //shared memory j1-j0
  //Todo figure out how to change it so it only loops the correct values
  j1 = min((blockIdx.x * njPerBlock + njPerBlock)-i, j1);
  __syncthreads();  


    #pragma unroll 16
    for (int k = j0; k < j1; k++) {
              body_body_interaction_nongb(ds_min, acc, jrk, pos, vel,
                                    shared_pos[k], shared_vel[k], EPS2);
        }

    __syncthreads();

    i += blockDim.x;
  }
  
  //Reduce the data 
  float4 *shared_acc = (float4*)&shared_pos[0];
  float4 *shared_jrk = (float4*)&shared_acc[blockDim.x*blockDim.y];
  int    *shared_ngb = (int*   )&shared_jrk[blockDim.x*blockDim.y];
  int    *shared_ofs = (int*   )&shared_ngb[blockDim.x*blockDim.y];
  float  *shared_ds  = (float* )&shared_ofs[blockDim.x*blockDim.y];
  acc.w = -acc.w;
//   jrk.w = __int_as_float(local_ngb_list[NGB_PB]);
//  jrk.w = local_ngb_list[NGB_PB];

  shared_acc[threadIdx.x] = acc;
  shared_jrk[threadIdx.x] = jrk;
//   shared_ngb[threadIdx.x] = n_ngb;
  shared_ofs[threadIdx.x] = 0;
  shared_ds [threadIdx.x] = ds_min;
  __syncthreads();


  for (int i = ni + threadIdx.x; i < ni*nThreads; i += ni) {

//       float4 acc1 = shared_acc[ajc(threadIdx.x, i)];
//       float4 jrk1 = shared_jrk[ajc(threadIdx.x, i)];
//       float  ds1  = shared_ds [ajc(threadIdx.x, i)];
      float4 acc1 = shared_acc[i];
      float4 jrk1 = shared_jrk[i];
      float  ds1  = shared_ds [i];
      
      acc.x += acc1.x;
      acc.y += acc1.y;
      acc.z += acc1.z;
      acc.w += acc1.w;
      
      jrk.x += jrk1.x;
      jrk.y += jrk1.y;
      jrk.z += jrk1.z;
      
      if (ds1  < ds_min) {
        jrk.w   = jrk1.w;
        ds_min  = ds1;
      }
    
  }//end reduction

 __syncthreads();

  
  if (threadIdx.x < ni) {
    //Convert results to double and write
    vel_i[offset  + blockIdx.x * blockDim.x + threadIdx.x].w = ds_min;
//     acc_i[blockIdx.x * blockDim.x + threadIdx.x] = (double4){acc.x, acc.y, acc.z, acc.w};
    acc_i[blockIdx.x * ni + threadIdx.x] = (double4){acc.x, acc.y, acc.z, acc.w};
    jrk_i[blockIdx.x * ni + threadIdx.x] = (double4){jrk.x, jrk.y, jrk.z, jrk.w};
  }

}


/*
 *  blockDim.x = #of block in previous kernel
 *  gridDim.x  = ni
 */ 
extern "C" __global__ void dev_reduce_forces(double4 *acc_i, 
                                             double4 *jrk_i,
                                             double  *ds_i,
                                             double4 *vel_i,
                                             int     offset_ds,
                                             int     offset,
                                             int     *ngb_list) {
//  extern __shared__ float4 shared_acc[];
//   __shared__ char shared_mem[NBLOCKS*(2*sizeof(float4) + 3*sizeof(int))];
//   float4* shared_acc = (float4*)&shared_mem;

  __shared__ float4     shared_acc[NBLOCKS];
  __shared__ float4     shared_jrk[NBLOCKS];
  __shared__ int        shared_ngb[NBLOCKS];
  __shared__ int        shared_ofs[NBLOCKS];
  __shared__ float      shared_ds[NBLOCKS];
  
//   float4 *shared_jrk = (float4*)&shared_acc[blockDim.x];
//   int    *shared_ngb = (int*   )&shared_jrk[blockDim.x];
//   int    *shared_ofs = (int*   )&shared_ngb[blockDim.x];
//   float  *shared_ds  = (float* )&shared_ofs[blockDim.x];
  
  int index = threadIdx.x * gridDim.x + blockIdx.x;

//   shared_acc[threadIdx.x] = acc_i[index];
//   shared_jrk[threadIdx.x] = jrk_i[index];
//   shared_ds [threadIdx.x] = vel_i[offset_ds + index].w;

  //Convert the data to floats
  shared_acc[threadIdx.x] = (float4){acc_i[index].x, acc_i[index].y, acc_i[index].z, acc_i[index].w};
  shared_jrk[threadIdx.x] = (float4){jrk_i[index].x, jrk_i[index].y, jrk_i[index].z, jrk_i[index].w};
  shared_ds [threadIdx.x] = (float)vel_i[offset_ds + index].w;  //TODO JB dont we miss the value at vel_i[0 + x] this way?


  int ngb_index = threadIdx.x * NGB_PB + blockIdx.x * NGB_PB*NBLOCKS;
  shared_ngb[threadIdx.x] = ngb_list[ngb_index];
  shared_ofs[threadIdx.x] = 0;
         
  __syncthreads();


  int n_ngb = shared_ngb[threadIdx.x];
  if (threadIdx.x == 0) {
    float4 acc0 = shared_acc[0];
    float4 jrk0 = shared_jrk[0];
    float  ds0 = shared_ds [0];

    for (int i = 1; i < blockDim.x; i++) {
      acc0.x += shared_acc[i].x;
      acc0.y += shared_acc[i].y;
      acc0.z += shared_acc[i].z;
      acc0.w += shared_acc[i].w;

      jrk0.x += shared_jrk[i].x;
      jrk0.y += shared_jrk[i].y;
      jrk0.z += shared_jrk[i].z;

      if (shared_ds[i] < ds0) {
        ds0    = shared_ds[i];
        jrk0.w = shared_jrk[i].w;
      }

      shared_ofs[i] = min(n_ngb + 1, NGB_PP);
      n_ngb += shared_ngb[i];

    }
    n_ngb = min(n_ngb, NGB_PP);

    jrk0.w = (int)__float_as_int(jrk0.w);

    //Store the results
    acc_i[blockIdx.x] = (double4){acc0.x, acc0.y, acc0.z, acc0.w};
    jrk_i[blockIdx.x] = (double4){jrk0.x, jrk0.y, jrk0.z, jrk0.w};;
    ds_i [blockIdx.x] = ds0;
  }
  __syncthreads();

  offset += blockIdx.x * NGB_PP + shared_ofs[threadIdx.x];
  int offset_end;
  if (threadIdx.x == 0) {
    shared_ofs[0] = offset + NGB_PP;
    ngb_list[offset++] = n_ngb;
  }
  __syncthreads();
  
  offset_end = shared_ofs[0];
  
  n_ngb = shared_ngb[threadIdx.x];

  for (int i = 0; i < n_ngb; i++)
    if (offset + i < offset_end){
      ngb_list[offset + i] = ngb_list[ngb_index + 1 + i];
    }
  
}


/*
 * Function that moves the (changed) j-particles
 * to the correct address location.
*/

extern "C" __global__ void dev_copy_particles(int nj, int nj_max,
                                              double4   *pos_j, 
                                              double4   *pos_j_temp,
                                              int       *address_j,
                                              double2   *t_j,
                                              double4   *Ppos_j, 
                                              double4   *Pvel_j,                                              
                                              double4   *vel_j,
                                              double4   *acc_j,
                                              double4   *jrk_j,
                                              int       *id_j,
                                              double2   *t_j_temp,                                              
                                              double4   *vel_j_temp,
                                              double4   *acc_j_temp,
                                              double4   *jrk_j_temp,
                                              int       *id_j_temp) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  //Copy the changed particles
  if (index < nj)
  {
    t_j  [address_j[index]] = t_j_temp[index];

    Ppos_j[address_j[index]] = pos_j_temp[index];
     pos_j[address_j[index]] = pos_j_temp[index];

    Pvel_j[address_j[index]] = vel_j_temp[index];
     vel_j[address_j[index]] = vel_j_temp[ index];

    acc_j[address_j[index]]  = acc_j_temp[index];
    jrk_j[address_j[index]]  = jrk_j_temp[index];

    id_j[address_j[index]]   = id_j_temp[index];
  }
}

/*

Function to predict the particles
DS version

*/
extern "C" __global__ void dev_predictor(int nj,
                              double  t_i_d,
                              double2 *t_j,
                              double4 *Ppos_j,
                              double4 *Pvel_j,
                              double4 *pos_j, 
                              double4 *vel_j,
                              double4 *acc_j,
                              double4 *jrk_j) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (index < nj) {

    //Convert the doubles to DS
    DS2 t;
    t.x = to_DS(t_j[index].x);
    t.y = to_DS(t_j[index].y);

    DS t_i;
    t_i = to_DS(t_i_d);

    DS4 pos;
    pos.x = to_DS(pos_j[index].x); pos.y = to_DS(pos_j[index].y);
    pos.z = to_DS(pos_j[index].z); pos.w = to_DS(pos_j[index].w);

    float4 vel = (float4){vel_j[index].x, vel_j[index].y, vel_j[index].z, vel_j[index].w};
    float4 acc = (float4){acc_j[index].x, acc_j[index].y, acc_j[index].z, acc_j[index].w};
    float4 jrk = (float4){jrk_j[index].x, jrk_j[index].y, jrk_j[index].z, jrk_j[index].w};
  
    float dt = (t_i.x - t.x.x) + (t_i.y - t.x.y);
    float dt2 = dt*dt/2.0f;
    float dt3 = dt2*dt/3.0f;
    
    pos.x  = dsadd(pos.x, vel.x * dt + acc.x * dt2 + jrk.x * dt3);
    pos.y  = dsadd(pos.y, vel.y * dt + acc.y * dt2 + jrk.y * dt3);
    pos.z  = dsadd(pos.z, vel.z * dt + acc.z * dt2 + jrk.z * dt3);

    vel.x += acc.x * dt + jrk.x * dt2;
    vel.y += acc.y * dt + jrk.y * dt2;
    vel.z += acc.z * dt + jrk.z * dt2;


    Ppos_j[index].x = to_double(pos.x); Ppos_j[index].y = to_double(pos.y);
    Ppos_j[index].z = to_double(pos.z); Ppos_j[index].w = to_double(pos.w);            

    Pvel_j[index] = (double4){vel.x, vel.y, vel.z, vel.w};
  }
  __syncthreads();
}

