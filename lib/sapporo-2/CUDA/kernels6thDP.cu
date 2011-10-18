/*

Sapporo 2 device kernels

Version 1.0
CUDA Double Precision 6th order hermite 

*/

//#include <stdio.h>

#include "../include/defines.h"

#if 1
__device__ void body_body_interaction(double  &ds_min,
                                      int     &n_ngb,
                                      int     *ngb_list,
                                      double4 &accNew_i, 
                                      double4 &jrkNew_i,
                                      double4 &snpNew_i,
                                      double4 pos_i, 
                                      double4 vel_i,
                                      double4 acc_i,                                        
                                      double4 pos_j, 
                                      double4 vel_j,
                                      double3 acc_j,
                                      int pjID,
                                      int piID,
                                      double &EPS2, double &accTest) {

  double3 dr = {pos_j.x - pos_i.x, pos_j.y - pos_i.y, pos_j.z - pos_i.z}; //3FLOP

  double ds2 = ((dr.x*dr.x + (dr.y*dr.y)) + dr.z*dr.z) + EPS2; //5 FLOP

  bool ngb = true;

  if (ngb) {
    //Neighbour list
    if (ds2 <= pos_i.w) {
      if (n_ngb < NGB_PB) {
        if(pjID != piID) //Prevent adding self as neighbour
          ngb_list[n_ngb++] =  pjID;
      }
    }

    //Nearest neighbour
    if (ds2 < ds_min*((pjID != piID))) {
      ds_min  = ds2;
      ngb_list[NGB_PB] = pjID;
    }    
  }

//   double inv_ds  = rsqrt(ds2 + EPS2) * (pjID != piID);
//      double inv_ds  = rsqrt(ds2) * (pjID != piID);
  double inv_ds  = (1.0 / sqrt(ds2)) * (pjID != piID);

  //TODO make sure the above trick still works on Fermi devices 
  //and especially for double precision calculations. Otherwise we can use:
//   if((ds2 + EPS2) == 0)
//   if((ds2) == 0)
//     inv_ds = 0;


  double mass    = pos_j.w;
  double inv_ds2 = inv_ds*inv_ds;                         // 1 FLOP
  double inv_ds3 = mass * inv_ds*inv_ds2;                 // 2 FLOP
  
  // 3*4 + 3 = 15 FLOP
  accTest = (inv_ds3 * dr.x);
  accNew_i.x = ((inv_ds3 * dr.x) + accNew_i.x);
  accNew_i.y = ((inv_ds3 * dr.y) + accNew_i.y);
  accNew_i.z = ((inv_ds3 * dr.z) + accNew_i.z);
  
  accNew_i.w = (mass * inv_ds  + accNew_i.w); //Potential

  double3 dv  = {vel_j.x - vel_i.x, vel_j.y - vel_i.y, vel_j.z - vel_i.z}; //3 FLOP
  double3 da  = {acc_j.x - acc_i.x, acc_j.y - acc_i.y, acc_j.z - acc_i.z}; //3 FLOP
  double  v2  = dv.x*dv.x + dv.y*dv.y + dv.z*dv.z;
  double  ra  = (dr.x*da.x) + (dr.y*da.y) + (dr.z*da.z);

//   if(threadIdx.x == 0 && blockIdx.x == 0){
//   printf("posx: %g %g %g posj %g %g %g dx: %g %g %g eps: %g  ds2: %g inv_ds: %g inv_ds3: %g\n", 
//           pos_i.x, pos_i.y, pos_i.z, 
//           pos_j.x, pos_j.y, pos_j.z, 
//           dr.x, dr.y, dr.z, 
//           EPS2, ds2, inv_ds, inv_ds3);
// 
//   printf("accix: %f accj: %f accdax: %f  dvx: %f  v2: %f \n", acc_i.x, acc_j.x, da.x, dv.x, v2);
// }


  double alpha = (((dr.x*dv.x) + dr.y*dv.y) + dr.z*dv.z) * inv_ds2;
  double beta  = (v2 + ra) * inv_ds2 + alpha * alpha;

  //Jerk
  alpha       *= -3.0;
  double3     jerk;
  jerk.x     = (inv_ds3 * dv.x) + alpha * (inv_ds3 * dr.x);
  jerk.y     = (inv_ds3 * dv.y) + alpha * (inv_ds3 * dr.y);
  jerk.z     = (inv_ds3 * dv.z) + alpha * (inv_ds3 * dr.z);

  //Snap
  alpha       *= 2.0;
  beta        *= -3.0;
  snpNew_i.x = snpNew_i.x + (inv_ds3 * da.x) + alpha * jerk.x + beta * (inv_ds3 * dr.x);
  snpNew_i.y = snpNew_i.y + (inv_ds3 * da.y) + alpha * jerk.y + beta * (inv_ds3 * dr.y);
  snpNew_i.z = snpNew_i.z + (inv_ds3 * da.z) + alpha * jerk.z + beta * (inv_ds3 * dr.z);

  //Had to reuse jerk for snap so only add to total now
  jrkNew_i.x += jerk.x;
  jrkNew_i.y += jerk.y;
  jrkNew_i.z += jerk.z;


//    if(threadIdx.x == 0 && blockIdx.x == 0){
// //     printf("ax: %f jx: %f sx: %f  \n", (inv_ds3 * dr.x), jerk.x, ( inv_ds3 * da.x) + alpha * jerk.x + beta * (inv_ds3 * dr.x));
//     printf("snpNew_i.x %f \t %f \n", snpNew_i.x,(inv_ds3 * da.x) + alpha * jerk.x + beta * (inv_ds3 * dr.x) );
// 
//  }
//    if(threadIdx.x == 0 && blockIdx.x == 0){
//     printf("acc.x %f \t new: %f \t %f %f\n", (inv_ds3 * dr.x), accNew_i.x, pos_j.x, pos_i.x);
// 
//  }


  // TOTAL 50 FLOP (or 60 FLOP if compared against GRAPE6)  
}
#else



__device__ void body_body_interaction(double  &ds_min,
                                      int     &n_ngb,
                                      int     *ngb_list,
                                      double4 &accNew_i, 
                                      double4 &jrkNew_i,
                                      double4 &snpNew_i,
                                      double4 pos_i, 
                                      double4 vel_i,
                                      double4 acc_i,                                        
                                      double4 pos_j, 
                                      double4 vel_j,
                                      double4 acc_j,
                                      int pjID,
                                      int piID,
                                      double &EPS2, double &accTest) {
/*
  double3 dr = {pos_j.x - pos_i.x, pos_j.y - pos_i.y, pos_j.z - pos_i.z}; //3FLOP

  double ds2 = ((dr.x*dr.x + (dr.y*dr.y)) + dr.z*dr.z) + EPS2; //5 FLOP*/

 

  double dx = pos_j.x - pos_i.x;
  double dy = pos_j.y - pos_i.y;
  double dz = pos_j.z - pos_i.z;

  double dvx = vel_j.x - vel_i.x;
  double dvy = vel_j.y - vel_i.y;
  double dvz = vel_j.z - vel_i.z;

  double dax = vel_j.x - vel_i.x;
  double day = acc_j.x - acc_i.x;
  double daz = acc_j.x - acc_i.x;

  double r2 = EPS2 + dx*dx + dy*dy + dz*dz;

 bool ngb = true;

  if (ngb) {
    //Neighbour list
    if (r2 <= pos_i.w) {
      if (n_ngb < NGB_PB) {
        if(pjID != piID) //Prevent adding self as neighbour
          ngb_list[n_ngb++] = pjID;
      }
    }

    //Nearest neighbour
    if (r2 < ds_min*(pjID != piID)) {
      ds_min  = r2;
      ngb_list[NGB_PB] = pjID;
    }    
  }



  double rv =  dx*dvx +  dy*dvy +  dz*dvz;
  double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
  double ra =  dx*dax +  dy*day +  dz*daz;

  double rinv2 = 1.0 / r2;
  double rinv = sqrt(rinv2);

  if(r2 == 0) rinv = 0;


  double alpha = rv * rinv2;
  double beta  = (v2 + ra) * rinv2 + alpha * alpha;
  rinv *= pos_j.w; //posjw is mass
  
  accNew_i.w += rinv; // potential

  double rinv3 = rinv * rinv2;

  double ax = rinv3 * dx;
  double ay = rinv3 * dy;
  double az = rinv3 * dz;
  alpha *= -3.0;
  double jx = rinv3 * dvx + alpha * ax;
  double jy = rinv3 * dvy + alpha * ay;
  double jz = rinv3 * dvz + alpha * az;
  alpha *= 2.0;
  beta  *= -3.0;
//   double sx = __dmul_rn(rinv3, dax) + __dmul_rn(alpha, jx) + __dmul_rn(beta, ax);
//   double sy = __dmul_rn(rinv3, day) + __dmul_rn(alpha, jy) + __dmul_rn(beta, ay);
//   double sz = __dmul_rn(rinv3, daz) + __dmul_rn(alpha, jz) + __dmul_rn(beta, az);
  double sx = rinv3*dax + (alpha*jx) + (beta*ax);
  double sy = rinv3*day + (alpha*jy) + (beta*ay);
  double sz = rinv3*daz + (alpha*jz) + (beta*az);

  accNew_i.x += ax;
  accNew_i.y += ay;
  accNew_i.z += az;

  jrkNew_i.x += jx;
  jrkNew_i.y += jy;
  jrkNew_i.z += jz;

  snpNew_i.x += sx;
  snpNew_i.y += sy;
  snpNew_i.z += sz;


//   double rinv1 = rsqrtf(r2);
// 
//  if(r2 == 0) rinv1 = 0.0;
// 
// 
//   double rinv2 = rinv1 * rinv1;
//   double alpha = (drdv)*rinv2;
//   double beta = (dvdv + drda)*rinv2 + alpha*alpha;
//   rinv1 *= pos_j.w;
//   double rinv3 = rinv1 * rinv2;
// 
//   accNew_i.w += rinv1; // potential
// 
//   // float pot = rinv1;
//   double ax = rinv3*dx;
//   double ay = rinv3*dy;
//   double az = rinv3*dz;
//   double jx = rinv3*dvx + (-3.0*alpha)*ax;
//   double jy = rinv3*dvy + (-3.0*alpha)*ay;
//   double jz = rinv3*dvz + (-3.0*alpha)*az;
//   double sx = rinv3*dax + (-6.0*alpha)*jx + (-3.0*beta)*ax;
//   double sy = rinv3*day + (-6.0*alpha)*jy + (-3.0*beta)*ay;
//   double sz = rinv3*daz + (-6.0*alpha)*jz + (-3.0*beta)*az;

//   accNew_i.x += ax;
//   accNew_i.y += ay;
//   accNew_i.z += az;
// 
//   jrkNew_i.x += jx;
//   jrkNew_i.y += jy;
//   jrkNew_i.z += jz;
// 
//   snpNew_i.x += sx;
//   snpNew_i.y += sy;
//   snpNew_i.z += sz;

// acc 203 37.000000: -0.1866005300987792 -0.2235977706495966 1.58119      phi: -0.257288  jrk:116.635 96.8278 -224.822    snp:4030.343359079666 2884.44 -17115.6  r2: 0.000244325
// acc 351 765.000000: 0.5314722333495085 1.089017857875183 0.777513       phi: -0.257942  jrk:-92.2404 162.941 249.098    snp:-1820.064827138607 -11322 -11174.2  r2: 0.000244345


//    if(threadIdx.x == 0 && blockIdx.x == 0){
// //     printf("ax: %f jx: %f sx: %f  \n", (inv_ds3 * dr.x), jerk.x, ( inv_ds3 * da.x) + alpha * jerk.x + beta * (inv_ds3 * dr.x));
//     printf("snpNew_i.x %f \t %f \n", snpNew_i.x,(inv_ds3 * da.x) + alpha * jerk.x + beta * (inv_ds3 * dr.x) );
// 
//  }
//    if(threadIdx.x == 0 && blockIdx.x == 0){
//     printf("acc.x %f \t new: %f \t %f %f\n", (inv_ds3 * dr.x), accNew_i.x, pos_j.x, pos_i.x);
// 
//  }


  // TOTAL 50 FLOP (or 60 FLOP if compared against GRAPE6)  
}
#endif

/*
 *  blockDim.x = ni
 *  gridDim.x  = 16, 32, 64, 128, etc. 
 */ 

//Kernel for same softening value for all particles

//TODO should make this depending on if we use Fermi or GT80/GT200
//#define ajc(i, j) (i + __mul24(blockDim.x,j))
#define ajc(i, j) (i + blockDim.x*j)


extern "C" __global__ void dev_evaluate_gravity2(
                                     int        nj_total, 
                                     int        nj,
                                     int        offset,
                                     double4    *pos_j, 
                                     double4    *vel_j,
                                     int        *id_j,
                                     double4    *pos_i,
                                     double4    *vel_i,
                                     double4    *acc_i, 
                                     double4    *jrk_i,
                                     int        *id_i,
                                     int        *ngb_list,
                                     double     EPS2,
                                     double4    *acc_j,       
                                     double4    *snp_i) {

if(threadIdx.x == 0)
  int x =1;
// printf("threadIdx.x: %d \n", blockIdx.x);

}

extern "C" __global__ void dev_evaluate_gravity(
                                     int        nj_total, 
                                     int        nj,
                                     int        offset,
                                     double4    *pos_j, 
                                     double4    *pos_i,
                                     double4    *acc_i,
                                     double     EPS2,
                                     double4    *vel_j,
                                     int        *id_j,
                                     double4    *vel_i,
                                     double4    *jrk_i,
                                     int        *id_i,
                                     int        *ngb_list,
                                     double4    *acc_j,       
                                     double4    *snp_i) {
//   extern __shared__ double4 shared_pos[];
//    __shared__ char shared_mem[NTHREADS*(sizeof(double4) + sizeof(double4) + sizeof(double4) + sizeof(int)*2 + sizeof(double))];
//   __shared__ char shared_mem[NTHREADS*(sizeof(double4) + sizeof(double4) + sizeof(double4))];
  __shared__ char shared_mem[NTHREADS*(sizeof(double4) + sizeof(double4) + sizeof(double3) + sizeof(int))];
  double4 *shared_pos = (double4*)&shared_mem[0];
  double4 *shared_vel = (double4*)&shared_pos[blockDim.x*blockDim.y];
  double3 *shared_acc = (double3*)&shared_vel[blockDim.x*blockDim.y];
  int     *shared_id  = (int*)&shared_acc[blockDim.x*blockDim.y];


  int local_ngb_list[NGB_PB + 1];
  int n_ngb = 0;

  //Read the i-particle properties for the particle that belongs to
  //this thread
  double4 pos = pos_i[threadIdx.x];
  double4 vel = vel_i[threadIdx.x];
  double4 acc = acc_i[threadIdx.x];

  acc.w = __int_as_float(id_i[threadIdx.x]);

  //Set the softening for the i-particle
  EPS2 = vel.w;

  //Combine the particle id into the w part of the position
  int particleID = id_i[threadIdx.x];

  #define LARGEnum 1e10f
  double ds_min = LARGEnum;
  
  double4 accNew = {0.0f, 0.0f, 0.0f, 0.0f};
  double4 jrkNew = {0.0f, 0.0f, 0.0f, 0.0f};
  double4 snpNew = {0.0f, 0.0f, 0.0f, 0.0f};

  int count = 0;
  int i = blockIdx.x * (nj*blockDim.y) + nj*threadIdx.y;
  int tile = 0;
  while (i <  blockIdx.x * (nj*blockDim.y) + nj*threadIdx.y + nj) { 
    if (i + threadIdx.x < nj_total) {
      //Read j-particles into shared memory
      shared_pos[ajc(threadIdx.x, threadIdx.y)] = pos_j[i + threadIdx.x];
      shared_vel[ajc(threadIdx.x, threadIdx.y)] = vel_j[i + threadIdx.x];
      shared_acc[ajc(threadIdx.x, threadIdx.y)] = (double3){acc_j[i + threadIdx.x].x, acc_j[i + threadIdx.x].y, acc_j[i + threadIdx.x].z};
      shared_id [ajc(threadIdx.x, threadIdx.y)] = id_j [i + threadIdx.x]; 

    } else {
      shared_pos[ajc(threadIdx.x, threadIdx.y)] = (double4){LARGEnum,LARGEnum,LARGEnum,0};

      shared_id[ajc(threadIdx.x, threadIdx.y)]    = -1; 
      shared_vel[ajc(threadIdx.x, threadIdx.y)]   = (double4){0.0, 0.0, 0.0, 0.0};
      shared_acc[ajc(threadIdx.x, threadIdx.y)]   = (double3){0.0, 0.0, 0.0};
    }
    __syncthreads();
    
    int j  = min(nj - tile*blockDim.x, blockDim.x);
    int j1 = (j/16)*16;

    double accTest = 0;

    #pragma unroll 16
    for (int k = 0; k < j1; k++) {
          body_body_interaction(ds_min, n_ngb, local_ngb_list,
                                accNew, jrkNew, snpNew, pos, vel, acc,
                                shared_pos[ajc(k, threadIdx.y)], shared_vel[ajc(k, threadIdx.y)], 
                                shared_acc[ajc(k, threadIdx.y)], 
                                shared_id[ajc(k, threadIdx.y)], particleID, EPS2, accTest);
        }
        
        for (int k = j1; k < j; k++) {
          body_body_interaction(ds_min, n_ngb, local_ngb_list,
                                accNew, jrkNew, snpNew, pos, vel, acc,
                                shared_pos[ajc(k, threadIdx.y)], shared_vel[ajc(k, threadIdx.y)], 
                                shared_acc[ajc(k, threadIdx.y)], 
                                shared_id[ajc(k, threadIdx.y)], particleID, EPS2, accTest);
        }
  
    __syncthreads();

    i += blockDim.x;
    tile++;
  }

  //Combine seperate results if more than one thread is used
  //per particle
  double4 *shared_acc2 = (double4*)&shared_pos[0];
  double4 *shared_snp = (double4*)&shared_acc2[blockDim.x*blockDim.y];
  
  acc.w = -acc.w;
  shared_acc2[ajc(threadIdx.x, threadIdx.y)] = accNew;
  shared_snp[ajc(threadIdx.x, threadIdx.y)]  = snpNew;
  __syncthreads();

  if (threadIdx.y == 0) {
    for (int i = 1; i < blockDim.y; i++) {
      double4 acc1 = shared_acc2[ajc(threadIdx.x, i)];
      double4 snp1 = shared_snp[ajc(threadIdx.x, i)];
     
      accNew.x += acc1.x;
      accNew.y += acc1.y;
      accNew.z += acc1.z;
      accNew.w += acc1.w;
      
      snpNew.x += snp1.x;
      snpNew.y += snp1.y;
      snpNew.z += snp1.z;
    }
  }
  __syncthreads();


  double4 *shared_jrk = (double4*)&shared_pos[0];
  int    *shared_ngb  = (int*   )&shared_jrk[blockDim.x*blockDim.y];
  int    *shared_ofs  = (int*   )&shared_ngb[blockDim.x*blockDim.y];
  double  *shared_ds  = (double* )&shared_ofs[blockDim.x*blockDim.y];
  
  jrkNew.w = __int_as_float(local_ngb_list[NGB_PB]);
  shared_jrk[ajc(threadIdx.x, threadIdx.y)] = jrkNew;
  shared_ngb[ajc(threadIdx.x, threadIdx.y)] = n_ngb;
  shared_ofs[ajc(threadIdx.x, threadIdx.y)] = 0;
  shared_ds [ajc(threadIdx.x, threadIdx.y)] = ds_min;
  __syncthreads();

  if (threadIdx.y == 0) {

    for (int i = 1; i < blockDim.y; i++) {
      double4 jrk1 = shared_jrk[ajc(threadIdx.x, i)];
      double  ds1  = shared_ds [ajc(threadIdx.x, i)];

      jrkNew.x += jrk1.x;
      jrkNew.y += jrk1.y;
      jrkNew.z += jrk1.z;


      if (ds1  < ds_min) {
        jrkNew.w   = jrk1.w;
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
    acc_i[blockIdx.x * blockDim.x + threadIdx.x] = accNew;
    jrk_i[blockIdx.x * blockDim.x + threadIdx.x] = jrkNew;
    snp_i[blockIdx.x * blockDim.x + threadIdx.x] = snpNew;
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

Double precision version
 */ 
extern "C" __global__ void dev_reduce_forces(double4 *acc_i, 
                                             double4 *jrk_i,
                                             double  *ds_i,
                                             double4 *vel_i,
                                             int     offset_ds,
                                             int     offset,
                                             int     *ngb_list,
                                             double4 *snp_i) {
  //NBLOCKS*(3*sizeof(double4) + 2*sizeof(int) + sizeof(double));   
//   extern __shared__ double4 shared_acc[];
 /* double4 *shared_jrk = (double4*)&shared_acc[blockDim.x];
  double4 *shared_snp = (double4*)&shared_jrk[blockDim.x];
  int    *shared_ngb = (int*   )&shared_snp[blockDim.x];
  int    *shared_ofs = (int*   )&shared_ngb[blockDim.x];
  double *shared_ds  = (double* )&shared_ofs[blockDim.x];*/

  __shared__ double4 shared_acc[NBLOCKS];
  __shared__ double4 shared_jrk[NBLOCKS];
  __shared__ double4 shared_snp[NBLOCKS];
  __shared__ int     shared_ngb[NBLOCKS];
  __shared__ int     shared_ofs[NBLOCKS];
  __shared__ double  shared_ds[NBLOCKS];

  int index = threadIdx.x * gridDim.x + blockIdx.x;

  //Convert the data to floats
  shared_acc[threadIdx.x] = acc_i[index];
  shared_jrk[threadIdx.x] = jrk_i[index];
  shared_snp[threadIdx.x] = snp_i[index];
  shared_ds [threadIdx.x] = (double)vel_i[offset_ds + index].w;  //TODO JB dont we miss the value at vel_i[0 + x] this way?


  int ngb_index = threadIdx.x * NGB_PB + blockIdx.x * NGB_PB*NBLOCKS;
  shared_ngb[threadIdx.x] = ngb_list[ngb_index];
  shared_ofs[threadIdx.x] = 0;
         
  __syncthreads();

  int n_ngb = shared_ngb[threadIdx.x];
  if (threadIdx.x == 0) {
    double4 acc0 = shared_acc[0];
    double4 jrk0 = shared_jrk[0];
    double4 snp0 = shared_snp[0];
    double  ds0 = shared_ds [0];

    for (int i = 1; i < blockDim.x; i++) {
      acc0.x += shared_acc[i].x;
      acc0.y += shared_acc[i].y;
      acc0.z += shared_acc[i].z;
      acc0.w += shared_acc[i].w;

      jrk0.x += shared_jrk[i].x;
      jrk0.y += shared_jrk[i].y;
      jrk0.z += shared_jrk[i].z;

      snp0.x += shared_snp[i].x;
      snp0.y += shared_snp[i].y;
      snp0.z += shared_snp[i].z;


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
    acc_i[blockIdx.x] = acc0;
    jrk_i[blockIdx.x] = jrk0;
    snp_i[blockIdx.x] = snp0;
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


/*
 * Function that moves the (changed) j-particles
 * to the correct address location.
*/
extern "C" __global__ void dev_copy_particles(int nj, 
                                              int nj_max,
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
                                              int       *id_j_temp,
                                              double4   *Pacc_j,
                                              double4   *snp_j,
                                              double4   *crk_j,
                                              double4   *snp_j_temp,
                                              double4   *crk_j_temp) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  //Copy the changed particles
  if (index < nj)
  {
    t_j  [address_j[index]] = t_j_temp[index];

    Ppos_j[address_j[index]] = pos_j_temp[index];
     pos_j[address_j[index]] = pos_j_temp[index];

    Pvel_j[address_j[index]] = vel_j_temp[index];
     vel_j[address_j[index]] = vel_j_temp[ index];

    Pacc_j[address_j[index]] = acc_j_temp[index];
     acc_j[address_j[index]] = acc_j_temp[index];

    jrk_j[address_j[index]]  = jrk_j_temp[index];
    snp_j[address_j[index]]  = snp_j_temp[index];
    crk_j[address_j[index]]  = crk_j_temp[index];

    id_j[address_j[index]]   = id_j_temp[index];
  }
}

/*

Function to predict the particles
Double Precision version

6th order hermite
*/
extern "C" __global__ void dev_predictor(int nj,
                                        double  t_i_d,
                                        double2 *t_j,
                                        double4 *Ppos_j,
                                        double4 *Pvel_j,
                                        double4 *pos_j, 
                                        double4 *vel_j,
                                        double4 *acc_j,
                                        double4 *jrk_j,
                                        double4 *Pacc_j,
                                        double4 *snp_j,
                                        double4 *crk_j){
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (index < nj) {

    double dt = t_i_d  - t_j[index].x;
    double dt2 = (1./2.)*dt;
    double dt3 = (1./3.)*dt;
    double dt4 = (1./4.)*dt;
    double dt5 = (1./5.)*dt;

    double4  pos         = pos_j[index];
    double4  vel         = vel_j[index];
    double4  acc         = acc_j[index];
    double4  jrk         = jrk_j[index];
    double4  snp         = snp_j[index];
    double4  crk         = crk_j[index];

    //Positions
    pos.x += dt  * (vel.x +  dt2 * (acc.x + dt3 * (jrk.x + 
             dt4 * (snp.x +  dt5 * (crk.x)))));
    pos.y += dt  * (vel.y +  dt2 * (acc.y + dt3 * (jrk.y + 
             dt4 * (snp.y +  dt5 * (crk.y)))));
    pos.z += dt  * (vel.z +  dt2 * (acc.z + dt3 * (jrk.z + 
             dt4 * (snp.z +  dt5 * (crk.z)))));
    Ppos_j[index] = pos;

    //Velocities
    vel.x += dt * (acc.x + dt2 * (jrk.x + 
             dt3 * (snp.x +  dt4 * (crk.x))));
    vel.y += dt * (acc.y + dt2 * (jrk.y + 
             dt3 * (snp.y +  dt4 * (crk.y))));
    vel.z += dt * (acc.z + dt2 * (jrk.z + 
             dt3 * (snp.z +  dt4 * (crk.z))));
    Pvel_j[index] = vel;


    //Accelerations
    acc.x += dt * (jrk.x + dt2 * (snp.x +  dt3 * (crk.x)));
    acc.y += dt * (jrk.y + dt2 * (snp.y +  dt3 * (crk.y)));
    acc.z += dt * (jrk.z + dt2 * (snp.z +  dt3 * (crk.z)));
    Pacc_j[index] = acc;

  }
}

extern "C" __global__ void dev_no_predictor(int nj,
                                        double  t_i_d,
                                        double2 *t_j,
                                        double4 *Ppos_j,
                                        double4 *Pvel_j,
                                        double4 *pos_j, 
                                        double4 *vel_j,
                                        double4 *acc_j,
                                        double4 *jrk_j,
                                        double4 *Pacc_j,
                                        double4 *snp_j,
                                        double4 *crk_j){
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (index < nj) {

    double4  pos         = pos_j[index];
    double4  vel         = vel_j[index];
    double4  acc         = acc_j[index];
    //Positions

    Ppos_j[index] = pos;
    Pvel_j[index] = vel;
    Pacc_j[index] = acc;
  }
}

