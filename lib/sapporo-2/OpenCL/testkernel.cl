/*

Sapporo 2 device kernels

Version 1.0
CUDA DoubleSingle kernels

*/

#pragma OPENCL EXTENSION cl_khr_fp64: enable

#define __syncthreads() barrier(CLK_LOCAL_MEM_FENCE)                                                                                                        
#define blockIdx_x  get_group_id(0)                                                                                                                         
#define blockIdx_y  get_group_id(1)                                                                                                                         
#define threadIdx_x get_local_id(0)                                                                                                                         
#define threadIdx_y get_local_id(1)                                                                                                                         
#define gridDim_x   get_num_groups(0)                                                                                                                       
#define gridDim_y   get_num_groups(1)                                                                                                                       
#define blockDim_x  get_local_size(0)                                                                                                                       
#define blockDim_y  get_local_size(1)  


// #define NGB_PB 256
// #define NGB_PP 256
// 
// #define NBLOCKS 30

#include "include/defines.h"

typedef float2 DS;  // double single;

typedef struct DS4 {
  DS x, y, z, w;
} DS4;
typedef struct DS2 {
  DS x, y;
} DS2;


__inline DS to_DS(double a) {
  DS b;
  b.x = (float)a;
  b.y = (float)(a - b.x);
  return b;
}

__inline double to_double(DS a) {
  double b;
  b = (double)((double)a.x + (double)a.y);
  return b;
}


// This function computes c = a + b.
__inline DS dsaddds(DS a, DS b) {
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
__inline DS dsadd(DS a, float b) {
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


__inline void body_body_interaction(float *ds_min,
                                    int   *n_ngb,
                                    __private int *ngb_list,
                                    float4 *acc_i, 
                                    float4 *jrk_i,
                                    DS4     pos_i, 
                                    float4  vel_i,
                                    DS4     pos_j, 
                                    float4  vel_j,
                                    float  EPS2) {

  float4 dr = {(pos_j.x.x - pos_i.x.x) + (pos_j.x.y - pos_i.x.y),
               (pos_j.y.x - pos_i.y.x) + (pos_j.y.y - pos_i.y.y),
               (pos_j.z.x - pos_i.z.x) + (pos_j.z.y - pos_i.z.y), 0.0f};   // 3x3 = 9 FLOP

  float ds2 = ((dr.x*dr.x + (dr.y*dr.y)) + dr.z*dr.z);

  bool ngb = true;

  if (ngb) {
    if (ds2 <= pos_i.w.x) {
      if ((*n_ngb) < NGB_PB) {
        if(__float_as_int(pos_i.w.y) != __float_as_int(pos_j.w.y))      //Jeroen, is dit echt nodig?
          ngb_list[(*n_ngb)++] = __float_as_int(pos_j.w.y);
      }
    }

    if (ds2 < (*ds_min)*(__float_as_int(pos_i.w.y) != __float_as_int(pos_j.w.y))) {
      (*ds_min)  = ds2;
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
if((ds2 + EPS2) == 0)
  inv_ds = 0;


  float mass    = pos_j.w.x;
  float inv_ds2 = inv_ds*inv_ds;                         // 1 FLOP
  float inv_ds3 = mass * inv_ds*inv_ds2;                 // 2 FLOP
  
  // 3*4 + 3 = 15 FLOP
  (*acc_i).x = ((inv_ds3 * dr.x) +  (*acc_i).x);
  (*acc_i).y = ((inv_ds3 * dr.y) +  (*acc_i).y);
  (*acc_i).z = ((inv_ds3 * dr.z) +  (*acc_i).z);
  
  (*acc_i).w = (mass * inv_ds  + (*acc_i).w);

  float4 dv;    // 3 FLOP
  dv.x = vel_j.x - vel_i.x;
  dv.y = vel_j.y - vel_i.y;
  dv.z = vel_j.z - vel_i.z;

  float drdv = -3.0f * (inv_ds3*inv_ds2) * (((dr.x*dv.x) + dr.y*dv.y) + dr.z*dv.z);


  (*jrk_i).x = ((*jrk_i).x + inv_ds3 * dv.x) + drdv * dr.x;
  (*jrk_i).y = ((*jrk_i).y + inv_ds3 * dv.y) + drdv * dr.y;
  (*jrk_i).z = ((*jrk_i).z + inv_ds3 * dv.z) + drdv * dr.z;

  // TOTAL 50 FLOP (or 60 FLOP if compared against GRAPE6)  
}

/*
 *  blockDim.x = ni
 *  gridDim.x  = 16, 32, 64, 128, etc. 
 */ 

//TODO should make this depending on if we use Fermi or GT80/GT200
//#define ajc(i, j) (i + __mul24(blockDim.x,j))
// #define ajc(i, j) (i + __mul24(blockDim.x,j))
#define ajc(i, j) (i + blockDim_x*j)
__kernel void dev_evaluate_gravity(
                                     int        nj_total, 
                                     int        nj,
                                     int        offset,
                                     int        readOffset,
                                     __global double4    *pos_j, 
                                     __global double4    *vel_j,
                                     __global int        *id_j,
                                     __global double4    *pos_i,
                                     __global double4    *vel_i,
                                     __global double4    *acc_i, 
                                     __global double4    *jrk_i,
                                     __global int        *id_i,
                                     __global int        *ngb_list,
                                     double     EPS2_d,
                                     __local  DS4   *shared_pos) {
//   extern __shared__ DS4 shared_pos[];
  __local float4 *shared_vel = (__local float4*)&shared_pos[blockDim_x*blockDim_y];

//   int local_ngb_list[NGB_PB + 1];
  __private int local_ngb_list[NGB_PB + 1];
  int n_ngb = 0;

  float EPS2 = (float)EPS2_d;

  DS4 pos;
  pos.x = to_DS(pos_i[threadIdx_x].x); pos.y = to_DS(pos_i[threadIdx_x].y);
  pos.z = to_DS(pos_i[threadIdx_x].z); pos.w = to_DS(pos_i[threadIdx_x].w);

  //Combine the particle id into the w part of the position
  pos.w.y = __int_as_float(id_i[threadIdx_x]);

  float4 vel = (float4){vel_i[threadIdx_x].x, vel_i[threadIdx_x].y, vel_i[threadIdx_x].z, vel_i[threadIdx_x].w};

  #define LARGEnum 1e10f
  float ds_min = LARGEnum;
  
  float4 acc = {0.0f, 0.0f, 0.0f, 0.0f};
  float4 jrk = {0.0f, 0.0f, 0.0f, 0.0f};

  int i = blockIdx_x * (nj*blockDim_y) + nj*threadIdx_y;
  int tile = 0;

  int count = 0;

  while (i <  blockIdx_x * (nj*blockDim_y) + nj*threadIdx_y + nj) { 
    

    if (i + threadIdx_x < nj_total) {
      shared_pos[ajc(threadIdx_x, threadIdx_y)].x = to_DS(pos_j[readOffset + i + threadIdx_x].x);
      shared_pos[ajc(threadIdx_x, threadIdx_y)].y = to_DS(pos_j[readOffset + i + threadIdx_x].y);
      shared_pos[ajc(threadIdx_x, threadIdx_y)].z = to_DS(pos_j[readOffset + i + threadIdx_x].z);
      shared_pos[ajc(threadIdx_x, threadIdx_y)].w = to_DS(pos_j[readOffset + i + threadIdx_x].w);
      //Combine the particle id into the w part of the position
      shared_pos[ajc(threadIdx_x, threadIdx_y)].w.y = __int_as_float (id_j[readOffset + i + threadIdx_x]); 

      shared_vel[ajc(threadIdx_x, threadIdx_y)] = 
                (float4){vel_j[readOffset + i + threadIdx_x].x, vel_j[readOffset + i + threadIdx_x].y,
                         vel_j[readOffset + i + threadIdx_x].z, vel_j[readOffset + i + threadIdx_x].w};

    } else {
      shared_pos[ajc(threadIdx_x, threadIdx_y)].x = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx_x, threadIdx_y)].y = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx_x, threadIdx_y)].z = (float2){LARGEnum, 0.0f};
      shared_pos[ajc(threadIdx_x, threadIdx_y)].w = (float2){0.0f,  -1.0f}; 
      shared_vel[ajc(threadIdx_x, threadIdx_y)]   = (float4){0.0f, 0.0f, 0.0f, 0.0f};
    }
    __syncthreads();

    
    int j  = min(nj - tile*blockDim_x, blockDim_x);
    int j1 = (j/16)*16;

    #pragma unroll 16
       for (int k = 0; k < j1; k++) {
          body_body_interaction(&ds_min, &n_ngb, local_ngb_list,
                                    &acc, &jrk, pos, vel,
                                    shared_pos[ajc(k, threadIdx_y)], shared_vel[ajc(k, threadIdx_y)], EPS2);
        }
        
        for (int k = j1; k < j; k++) {
          body_body_interaction(&ds_min, &n_ngb, local_ngb_list,
                                    &acc, &jrk, pos, vel,
                                    shared_pos[ajc(k, threadIdx_y)], shared_vel[ajc(k, threadIdx_y)], EPS2);
        }
  

    __syncthreads();

    i += blockDim_x;
    tile++;
  } //end while


  __local float4 *shared_acc = (__local float4*)&shared_pos[0];
  __local float4 *shared_jrk = (__local float4*)&shared_acc[blockDim_x*blockDim_y];
  __local int    *shared_ngb = (__local int*   )&shared_jrk[blockDim_x*blockDim_y];
  __local int    *shared_ofs = (__local int*   )&shared_ngb[blockDim_x*blockDim_y];
  __local float  *shared_ds  = (__local float* )&shared_ofs[blockDim_x*blockDim_y];
  acc.w = -acc.w;
  jrk.w = __int_as_float(local_ngb_list[NGB_PB]);
//  jrk.w = local_ngb_list[NGB_PB];

  shared_acc[ajc(threadIdx_x, threadIdx_y)] = acc;
  shared_jrk[ajc(threadIdx_x, threadIdx_y)] = jrk;
  shared_ngb[ajc(threadIdx_x, threadIdx_y)] = n_ngb;
  shared_ofs[ajc(threadIdx_x, threadIdx_y)] = 0;
  shared_ds [ajc(threadIdx_x, threadIdx_y)] = ds_min;
  __syncthreads();

  if (threadIdx_y == 0) {

    for (int i = 1; i < blockDim_y; i++) {
      float4 acc1 = shared_acc[ajc(threadIdx_x, i)];
      float4 jrk1 = shared_jrk[ajc(threadIdx_x, i)];
      float  ds1  = shared_ds [ajc(threadIdx_x, i)];
      
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

      shared_ofs[ajc(threadIdx_x, i)] = min(n_ngb + 1, NGB_PB);
      n_ngb += shared_ngb[ajc(threadIdx_x, i)];
    }
    n_ngb  = min(n_ngb, NGB_PB);
  }
  __syncthreads();

  if (threadIdx_y == 0) {
    //Convert results to double and write
    vel_i[offset  + blockIdx_x * blockDim_x + threadIdx_x].w = ds_min;
    acc_i[blockIdx_x * blockDim_x + threadIdx_x] = (double4){acc.x, acc.y, acc.z, acc.w};
    jrk_i[blockIdx_x * blockDim_x + threadIdx_x] = (double4){jrk.x, jrk.y, jrk.z, jrk.w};
  }


  offset  = threadIdx_x * NBLOCKS*NGB_PB + blockIdx_x * NGB_PB;
  offset += shared_ofs[ajc(threadIdx_x, threadIdx_y)];
  
  if (threadIdx_y == 0)
    ngb_list[offset++] = n_ngb;
  
  n_ngb = shared_ngb[ajc(threadIdx_x, threadIdx_y)];
  for (int i = 0; i < n_ngb; i++) 
    ngb_list[offset + i] = local_ngb_list[i];


}




/*
 *  blockDim.x = #of block in previous kernel
 *  gridDim.x  = ni
 */ 
__kernel void dev_reduce_forces(__global double4 *acc_i, 
                                __global double4 *jrk_i,
                                __global double  *ds_i,
                                __global double4 *vel_i,
                                         int     offset_ds,
                                         int     offset,
                               __global  int     *ngb_list,
                               __local  float4   *shared_acc ) {
  
//    extern __shared__ float4 shared_acc[];
 __local  float4 *shared_jrk = (__local float4*)&shared_acc[blockDim_x];
 __local  int    *shared_ngb = (__local int*   )&shared_jrk[blockDim_x];
  __local int    *shared_ofs = (__local int*   )&shared_ngb[blockDim_x];
  __local float  *shared_ds  = (__local float* )&shared_ofs[blockDim_x];
  
  int index = threadIdx_x * gridDim_x + blockIdx_x;

//   shared_acc[threadIdx.x] = acc_i[index];
//   shared_jrk[threadIdx.x] = jrk_i[index];
//   shared_ds [threadIdx.x] = vel_i[offset_ds + index].w;

  //Convert the data to floats
  shared_acc[threadIdx_x] = (float4){acc_i[index].x, acc_i[index].y, acc_i[index].z, acc_i[index].w};
  shared_jrk[threadIdx_x] = (float4){jrk_i[index].x, jrk_i[index].y, jrk_i[index].z, jrk_i[index].w};
  shared_ds [threadIdx_x] = (float)vel_i[offset_ds + index].w;  //TODO JB dont we miss the value at vel_i[0 + x] this way?


  int ngb_index = threadIdx_x * NGB_PB + blockIdx_x * NGB_PB*NBLOCKS;
  shared_ngb[threadIdx_x] = ngb_list[ngb_index];
  shared_ofs[threadIdx_x] = 0;
         
  __syncthreads();


  int n_ngb = shared_ngb[threadIdx_x];
  if (threadIdx_x == 0) {
    float4 acc0 = shared_acc[0];
    float4 jrk0 = shared_jrk[0];
    float  ds0 = shared_ds [0];

    for (int i = 1; i < blockDim_x; i++) {
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
    acc_i[blockIdx_x] = (double4){acc0.x, acc0.y, acc0.z, acc0.w};
    jrk_i[blockIdx_x] = (double4){jrk0.x, jrk0.y, jrk0.z, jrk0.w};;
    ds_i [blockIdx_x] = ds0;
  }
  __syncthreads();


  offset += blockIdx_x * NGB_PP + shared_ofs[threadIdx_x];
  int offset_end;
  if (threadIdx_x == 0) {
    shared_ofs[0] = offset + NGB_PP;
    ngb_list[offset++] = n_ngb;
  }
  __syncthreads();
  
  offset_end = shared_ofs[0];
  
  n_ngb = shared_ngb[threadIdx_x];
  for (int i = 0; i < n_ngb; i++)
    if (offset + i < offset_end)
      ngb_list[offset + i] = ngb_list[ngb_index + 1 + i];
  
}


/*
 * Function that moves the (changed) j-particles
 * to the correct address location.
*/
__kernel void dev_copy_particles(int nj, int nj_max,
                                 __global             int       *address_j,
                                 __global             double2   *t_j,
                                 __global             double4   *Ppos_j, 
                                 __global             double4   *Pvel_j,
                                 __global             double4   *pos_j, 
                                 __global             double4   *vel_j,
                                 __global             double4   *acc_j,
                                 __global             double4   *jrk_j,
                                 __global             int       *id_j,
                                 __global             double2   *t_j_temp,
                                 __global             double4   *pos_j_temp,
                                 __global             double4   *vel_j_temp,
                                 __global             double4   *acc_j_temp,
                                 __global             double4   *jrk_j_temp,
                                  __global            int       *id_j_temp) {
  int index = blockIdx_x * blockDim_x + threadIdx_x;
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
__kernel void dev_predictor(int nj,
                              double  t_i_d,
                            __global  double2 *t_j,
                            __global  double4 *Ppos_j,
                            __global  double4 *Pvel_j,
                            __global  double4 *pos_j, 
                            __global  double4 *vel_j,
                            __global  double4 *acc_j,
                            __global  double4 *jrk_j) {
  int index = blockIdx_x * blockDim_x + threadIdx_x;
  
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

