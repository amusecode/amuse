#ifndef _DEV_EVALUATE_GRAVITY_CU_
#define _DEV_EVALUATE_GRAVITY_CU_

__constant__ float EPS2;

typedef float2 DS;  // double single;

struct DS4 {
  DS x, y, z, w;
};
struct DS2 {
  DS x, y;
};


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



template<bool ngb>
__device__ void body_body_interaction(float &ds_min,
                                      int &n_ngb, int *ngb_list,
                                      float4 &acc_i, float4 &jrk_i,
                                      DS4     pos_i, float4  vel_i,
                                      DS4     pos_j, float4  vel_j) {

  float3 dr = {(pos_j.x.x - pos_i.x.x) + (pos_j.x.y - pos_i.x.y),
                (pos_j.y.x - pos_i.y.x) + (pos_j.y.y - pos_i.y.y),
                (pos_j.z.x - pos_i.z.x) + (pos_j.z.y - pos_i.z.y)};   // 3x3 = 9 FLOP

  
  float ds2 = ((dr.x*dr.x + (dr.y*dr.y)) + dr.z*dr.z);

  if (ngb) {

    if (ds2 <= pos_i.w.x) {
      if (n_ngb < NGB_PB) {
        if(__float_as_int(pos_i.w.y) != __float_as_int(pos_j.w.y))      //Jeroen, prevent self on neighbour list
          ngb_list[n_ngb++] = __float_as_int(pos_j.w.y);
      }
    }

    if (ds2 < ds_min*(__float_as_int(pos_i.w.y) != __float_as_int(pos_j.w.y))) {
      ds_min  = ds2;
      ngb_list[NGB_PB] = __float_as_int(pos_j.w.y);
    }
    
  }
  
  float inv_ds  = rsqrt(ds2 + EPS2) * (__float_as_int(pos_i.w.y) != __float_as_int(pos_j.w.y));
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

/*
 *  blockDim.x = ni
 *  gridDim.x  = 16, 32, 64, 128, etc. 
 */ 

#define ajc(i, j) (i + __mul24(blockDim.x,j))
template<bool ngb>
__global__ void dev_evaluate_gravity(int nj_total, int nj,
                                     int offset,
                                     DS4    *pos_j, 
                                     float4 *vel_j,
                                     DS4    *pos_i,
                                     float4 *vel_i,
                                     float4 *acc_i, 
                                     float4 *jrk_i,
                                     int *ngb_list) {
  extern __shared__ DS4 shared_pos[];
  float4 *shared_vel = (float4*)&shared_pos[blockDim.x*blockDim.y];

  int local_ngb_list[NGB_PB + 1];
  int n_ngb = 0;

  DS4    pos = pos_i[threadIdx.x];
  float4 vel = vel_i[threadIdx.x];

#define LARGEnum 1e10f
  float ds_min = LARGEnum;
  
  float4 acc = {0.0f, 0.0f, 0.0f, 0.0f};
  float4 jrk = {0.0f, 0.0f, 0.0f, 0.0f};

  int i = blockIdx.x * (nj*blockDim.y) + nj*threadIdx.y;
  int tile = 0;
  while (i <  blockIdx.x * (nj*blockDim.y) + nj*threadIdx.y + nj) { 
    
    
    if (i + threadIdx.x < nj_total) {
      shared_pos[ajc(threadIdx.x, threadIdx.y)] = pos_j[i + threadIdx.x];
      shared_vel[ajc(threadIdx.x, threadIdx.y)] = vel_j[i + threadIdx.x];
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

// #pragma unroll 16
    for (int k = 0; k < j1; k++) {
      body_body_interaction<ngb>(ds_min, n_ngb, local_ngb_list,
                                 acc, jrk, pos, vel,
                                 shared_pos[ajc(k, threadIdx.y)], shared_vel[ajc(k, threadIdx.y)]);
    }
    
    for (int k = j1; k < j; k++) {
      body_body_interaction<ngb>(ds_min, n_ngb, local_ngb_list,
                                 acc, jrk, pos, vel,
                                 shared_pos[ajc(k, threadIdx.y)], shared_vel[ajc(k, threadIdx.y)]);
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
  jrk.w = __int_as_float(local_ngb_list[NGB_PB]);

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
    vel_i[offset  + blockIdx.x * blockDim.x + threadIdx.x].w = ds_min;
    acc_i[blockIdx.x * blockDim.x + threadIdx.x]             = acc;
    jrk_i[blockIdx.x * blockDim.x + threadIdx.x]             = jrk;
  }
  
  offset  = threadIdx.x * NBLOCKS*NGB_PB + blockIdx.x * NGB_PB;
  offset += shared_ofs[ajc(threadIdx.x, threadIdx.y)];
  
  if (threadIdx.y == 0)
    ngb_list[offset++] = n_ngb;
  
  n_ngb = shared_ngb[ajc(threadIdx.x, threadIdx.y)];
  for (int i = 0; i < n_ngb; i++) 
    ngb_list[offset + i] = local_ngb_list[i];

}


/*
 *  blockDim.x = #of block in previous kernel
 *  gridDim.x  = ni
 */ 
__global__ void dev_reduce_forces(float4 *acc_i, 
                                  float4 *jrk_i,
                                  float  *ds_i,
                                  float4 *vel_i,
                                  int offset_ds,
                                  int offset,
                                  int *ngb_list) {
  
  extern __shared__ float4 shared_acc[];
  float4 *shared_jrk = (float4*)&shared_acc[blockDim.x];
  int    *shared_ngb = (int*   )&shared_jrk[blockDim.x];
  int    *shared_ofs = (int*   )&shared_ngb[blockDim.x];
  float  *shared_ds  = (float* )&shared_ofs[blockDim.x];
  
  int index = threadIdx.x * gridDim.x + blockIdx.x;
  shared_acc[threadIdx.x] = acc_i[index];
  shared_jrk[threadIdx.x] = jrk_i[index];
  shared_ds [threadIdx.x] = vel_i[offset_ds + index].w;

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

    acc_i[blockIdx.x] = acc0;
    jrk_i[blockIdx.x] = jrk0;
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
    if (offset + i < offset_end)
      ngb_list[offset + i] = ngb_list[ngb_index + 1 + i];
    
}

__global__ void dev_copy_particles(int nj, int nj_max,
                                   int *address_j,
                                   DS2 *t_j,
                                   DS4    *Ppos_j, 
                                   float4 *Pvel_j,
                                   DS4    *pos_j, 
                                   float4 *vel_j,
                                   float4 *acc_j,
                                   float4 *jrk_j) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index < nj) {
    t_j  [address_j[index]] = t_j  [nj_max + index];

    DS4    pos = pos_j[nj_max + index];
    float4 vel = vel_j[nj_max + index];

    Ppos_j[address_j[index]] = pos;
     pos_j[address_j[index]] = pos;

    Pvel_j[address_j[index]] = vel;
     vel_j[address_j[index]] = vel;

    acc_j[address_j[index]] = acc_j[nj_max + index];
    jrk_j[address_j[index]] = jrk_j[nj_max + index];
  }
  __syncthreads();
};
                                   

__global__ void dev_predictor(int nj,
                              DS  t_i,
                              DS2    *t_j,
                              DS4    *Ppos_j,
                              float4 *Pvel_j,
                              DS4    *pos_j, 
                              float4 *vel_j,
                              float4 *acc_j,
                              float4 *jrk_j) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (index < nj) {

    DS2    t   = t_j  [index];
    DS4    pos = pos_j[index];
    float4 vel = vel_j[index];
    float4 acc = acc_j[index];
    float4 jrk = jrk_j[index];

    float dt = (t_i.x - t.x.x) + (t_i.y - t.x.y);
    float dt2 = dt*dt/2.0f;
    float dt3 = dt2*dt/3.0f;
    
    pos.x  = dsadd(pos.x, vel.x * dt + acc.x * dt2 + jrk.x * dt3);
    pos.y  = dsadd(pos.y, vel.y * dt + acc.y * dt2 + jrk.y * dt3);
    pos.z  = dsadd(pos.z, vel.z * dt + acc.z * dt2 + jrk.z * dt3);

    vel.x += acc.x * dt + jrk.x * dt2;
    vel.y += acc.y * dt + jrk.y * dt2;
    vel.z += acc.z * dt + jrk.z * dt2;
   
    Ppos_j[index] = pos;
    Pvel_j[index] = vel;
  }
  __syncthreads();
}


#endif
