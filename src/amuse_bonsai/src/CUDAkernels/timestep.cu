#ifdef __DEVICE_EMULATION__
  #define EMUSYNC __syncthreads();
#else
  #define EMUSYNC
#endif

#include "support_kernels.cu"

//Reduce function to get the minimum timestep
__device__ void get_TnextD(const int n_bodies,
                           double2 *time,
                           double *tnext, volatile double *sdata) {
  //float2 time : x is time begin, y is time end

  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  const int blockSize   = blockDim.x;
  unsigned int tid      = threadIdx.x;
  unsigned int i        = blockIdx.x*(blockSize*2) + threadIdx.x;
  unsigned int gridSize = blockSize*2*gridDim.x;
  sdata[tid] = 1.0e10f;
  double tmin = 1.0e10f;

  // we reduce multiple elements per thread.  The number is determined by the
  // number of active thread blocks (via gridSize).  More blocks will result
  // in a larger gridSize and therefore fewer elements per thread
  while (i < n_bodies) {
    if (i             < n_bodies) tmin = fmin(tmin, time[i            ].y);
    if (i + blockSize < n_bodies) tmin = fmin(tmin, time[i + blockSize].y);

    i += gridSize;
  }

  sdata[tid] = tmin;
  __syncthreads();

  // do reduction in shared mem
  if (blockSize >= 512) { if (tid < 256) { sdata[tid] = tmin = fmin(tmin, sdata[tid + 256]); } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdata[tid] = tmin = fmin(tmin, sdata[tid + 128]); } __syncthreads(); }
  if (blockSize >= 128) { if (tid <  64) { sdata[tid] = tmin = fmin(tmin, sdata[tid +  64]); } __syncthreads(); }
#ifndef __DEVICE_EMULATION__
  if (tid < 32)
#endif
    {
      if (blockSize >=  64) { sdata[tid] = tmin = fmin(tmin, sdata[tid + 32]); EMUSYNC; }
      if (blockSize >=  32) { sdata[tid] = tmin = fmin(tmin, sdata[tid + 16]); EMUSYNC; }
      if (blockSize >=  16) { sdata[tid] = tmin = fmin(tmin, sdata[tid +  8]); EMUSYNC; }
      if (blockSize >=   8) { sdata[tid] = tmin = fmin(tmin, sdata[tid +  4]); EMUSYNC; }
      if (blockSize >=   4) { sdata[tid] = tmin = fmin(tmin, sdata[tid +  2]); EMUSYNC; }
      if (blockSize >=   2) { sdata[tid] = tmin = fmin(tmin, sdata[tid +  1]); EMUSYNC; }
  }

  // write result for this block to global mem
  if (tid == 0) tnext[blockIdx.x] = sdata[0];
}

extern "C" __global__ void get_Tnext(const int n_bodies,
                                     double2 *time,
                                     double *tnext) {
  extern __shared__ double sdata[];
  get_TnextD(n_bodies, time, tnext, sdata);
}


//Reduce function to get the number of active particles
__device__ void get_nactiveD(const int n_bodies,
                                       uint *valid,
                                       uint *tnact, volatile int *sdataInt) {
  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  const int blockSize   = blockDim.x;
  unsigned int tid      = threadIdx.x;
  unsigned int i        = blockIdx.x*(blockSize*2) + threadIdx.x;
  unsigned int gridSize = blockSize*2*gridDim.x;
  sdataInt[tid] = 0;
  int sum       = 0;

  // we reduce multiple elements per thread.  The number is determined by the
  // number of active thread blocks (via gridSize).  More blocks will result
  // in a larger gridSize and therefore fewer elements per thread
  while (i < n_bodies) {
    if (i             < n_bodies) sum = sum + valid[i            ];
    if (i + blockSize < n_bodies) sum = sum + valid[i + blockSize];

    i += gridSize;
  }
  sdataInt[tid] = sum;
  __syncthreads();

  // do reduction in shared mem
  if (blockSize >= 512) { if (tid < 256) { sdataInt[tid] = sum = sum + sdataInt[tid + 256]; } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdataInt[tid] = sum = sum + sdataInt[tid + 128]; } __syncthreads(); }
  if (blockSize >= 128) { if (tid <  64) { sdataInt[tid] = sum = sum + sdataInt[tid +  64]; } __syncthreads(); }


#ifndef __DEVICE_EMULATION__
  if (tid < 32)
#endif
    {
      if (blockSize >=  64) { sdataInt[tid] = sum = sum + sdataInt[tid + 32]; EMUSYNC; }
      if (blockSize >=  32) { sdataInt[tid] = sum = sum + sdataInt[tid + 16]; EMUSYNC; }
      if (blockSize >=  16) { sdataInt[tid] = sum = sum + sdataInt[tid +  8]; EMUSYNC; }
      if (blockSize >=   8) { sdataInt[tid] = sum = sum + sdataInt[tid +  4]; EMUSYNC; }
      if (blockSize >=   4) { sdataInt[tid] = sum = sum + sdataInt[tid +  2]; EMUSYNC; }
      if (blockSize >=   2) { sdataInt[tid] = sum = sum + sdataInt[tid +  1]; EMUSYNC; }
  }

  // write result for this block to global mem
  if (tid == 0) tnact[blockIdx.x] = sdataInt[0];
}

//Reduce function to get the number of active particles
extern "C" __global__ void get_nactive(const int n_bodies,
                                       uint *valid,
                                       uint *tnact) {
  extern __shared__ int sdataInt[];
  get_nactiveD(n_bodies, valid, tnact, sdataInt);
}

#if 0

extern "C" __global__ void predict_particles(const int n_bodies,
                                             float tc,
                                             float tp,
                                             real4 *pos,
                                             real4 *vel,
                                             real4 *acc,
                                             float2 *time,
                                             uint  *body2grouplist,
                                             uint  *valid_list){
  const uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  const uint tid = threadIdx.x;
  const uint idx = bid * blockDim.x + tid;


  if (idx >= n_bodies) return;

  float4 p = pos [idx];
  float4 v = vel [idx];
  float4 a = acc [idx];
  float tb = time[idx].x;
  float te = time[idx].y;

  float dt_cb  = tc - tb;
  float dt_pb  = tp - tb;

  v.x -= a.x*dt_pb;
  v.y -= a.y*dt_pb;
  v.z -= a.z*dt_pb;

  p.x -= (v.x*dt_pb + a.x*dt_pb*dt_pb*0.5f);
  p.y -= (v.y*dt_pb + a.y*dt_pb*dt_pb*0.5f);
  p.z -= (v.z*dt_pb + a.z*dt_pb*dt_pb*0.5f);

  p.x += (v.x*dt_cb + a.x*dt_cb*dt_cb*0.5f);
  p.y += (v.y*dt_cb + a.y*dt_cb*dt_cb*0.5f);
  p.z += (v.z*dt_cb + a.z*dt_cb*dt_cb*0.5f);

  v.x += a.x*dt_cb;
  v.y += a.y*dt_cb;
  v.z += a.z*dt_cb;

  pos[idx] = p;
  vel[idx] = v;

  //Set the group to active if the time current = time end of
  //this particle. Can be that multiple particles write to the
  //same location but the net result is the same
  int grpID = body2grouplist[idx];
  if(tc == te)
  {
    valid_list[grpID] = grpID | (1 << 31);
  }
}
#endif

#if 1
extern "C" __global__ void predict_particles(const int n_bodies,
                                             double tc,
                                             double tp,
                                             real4 *pos,
                                             real4 *vel,
                                             real4 *acc,
                                             double2 *time,
                                             uint  *body2grouplist,
                                             uint  *valid_list,
                                             real4 *pPos,
                                             real4 *pVel){
  const uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  const uint tid = threadIdx.x;
  const uint idx = bid * blockDim.x + tid;


  if (idx >= n_bodies) return;

  float4 p = pos [idx];
  float4 v = vel [idx];
  float4 a = acc [idx];
  double tb = time[idx].x;
  double te = time[idx].y;

  float dt_cb  = (float)(tc - tb);
//   float dt_pb  = tp - tb;

  p.x += v.x*dt_cb + a.x*dt_cb*dt_cb*0.5f;
  p.y += v.y*dt_cb + a.y*dt_cb*dt_cb*0.5f;
  p.z += v.z*dt_cb + a.z*dt_cb*dt_cb*0.5f;

  v.x += a.x*dt_cb;
  v.y += a.y*dt_cb;
  v.z += a.z*dt_cb;

  pPos[idx] = p;
  pVel[idx] = v;


  //Set the group to active if the time current = time end of
  //this particle. Can be that multiple particles write to the
  //same location but the net result is the same
  int grpID = body2grouplist[idx];
  if(tc == te)
  {
    valid_list[grpID] = grpID | (1 << 31);
  }
}
#endif


extern "C"  __global__ void setActiveGroups(const int n_bodies,
                                            double tc,
                                            double2 *time,
                                            uint  *body2grouplist,
                                            uint  *valid_list){
  const uint bid = blockIdx.y * gridDim.x + blockIdx.x;
  const uint tid = threadIdx.x;
  const uint idx = bid * blockDim.x + tid;

  if (idx >= n_bodies) return;

  double te = time[idx].y;

  //Set the group to active if the time current = time end of
  //this particle. Can be that multiple particles write to the
  //same location but the net result is the same
  int grpID = body2grouplist[idx];


  //Test not only the article with current time, but any particle
  //with time diff less then 1/16k
/*  if(te-tc <= (1./16384))
  {
    valid_list[grpID] = grpID | (1 << 31);
  }
*/


  if(te <= tc)
  {
    valid_list[grpID] = grpID | (1 << 31);
  }

}


#ifdef _AMUSE_STOPPING_CONDITIONS_
extern "C" __global__ void correct_particles(const int n_bodies,
                                             double tc,
                                             double2 *time,
                                             uint   *active_list,
                                             real4 *vel,
                                             real4 *acc0,
                                             real4 *acc1,
                                             real4 *pos,
                                             real4 *pPos,
                                             real4 *pVel,
                                             int   *ngb,
                                             int   *pairDetection) {
  const int bid =  blockIdx.y *  gridDim.x +  blockIdx.x;
  const int tid =  threadIdx.y * blockDim.x + threadIdx.x;
  const int dim =  blockDim.x * blockDim.y;

  int idx = bid * dim + tid;
  if (idx >= n_bodies) return;

  //Check if particle is set to active during approx grav
  if (active_list[idx] != 1) return;


  float4 v  = vel [idx];
  float4 a0 = acc0[idx];
  float4 a1 = acc1[idx];
  double  tb = time[idx].x;
//   float  dt = time[idx].y;

  float dt_cb = (float)(tc - tb);

  //Store the predicted position as the one to use
  pos[idx] = pPos[idx];

  //Correct the position
  v = pVel[idx];

  dt_cb *= 0.5f;
  v.x += (a1.x - a0.x)*dt_cb;
  v.y += (a1.y - a0.y)*dt_cb;
  v.z += (a1.z - a0.z)*dt_cb;

  //Store the corrected velocity, accelaration and the new time step info
  vel [idx] = v;
  acc0[idx] = a1;

  //Code specific to stopping conditions
  int j = ngb[idx];
  
#if 1 
  if(j >= 0)    //Only check if we have a valid nearby neighbour
  {
    float4 posi = pPos[idx];
    float4 posj = pPos[j];
    float  radj = vel[j].w; //Particle radius is stored in w component of velocity
    float  radi = v.w;

    //Compute distance and compare to summed radius
    float ds2 = ((posi.x-posj.x)*(posi.x-posj.x)) +
                ((posi.y-posj.y)*(posi.y-posj.y)) +
                ((posi.z-posj.z)*(posi.z-posj.z));

    float rsum = radi + radj;
    if (ds2 <= rsum*rsum)
    {
      float4 veli = pVel[idx];
      float4 velj = pVel[j];    
  
      //Compute distance and compare to summed radius
      float r =         ((posi.x-posj.x)*(posi.x-posj.x)) +
                        ((posi.y-posj.y)*(posi.y-posj.y)) +
                        ((posi.z-posj.z)*(posi.z-posj.z));
      float v =         ((veli.x-velj.x)*(veli.x-velj.x)) +
                        ((veli.y-velj.y)*(veli.y-velj.y)) +
                        ((veli.z-velj.z)*(veli.z-velj.z));
      float vr =        ((posi.x-posj.x)*(veli.x-velj.x)) +
                        ((posi.y-posj.y)*(veli.y-velj.y)) +
                        ((posi.z-posj.z)*(veli.z-velj.z));
      //TODO remove these expensive operations instead just 
      //do vr*vr and EPS*EPS
      r = sqrt(r);
      v = sqrt(v);

      #define EPS 0.001   // see couple/multiples.py
//       if (abs(vr) < EPS*r*v)
      if(1) //JB: 9 sept 13 . Disabled untill we figure out why tests fail
      {
        //Collision detected, store the indices of the involved particles
        //Note that this will create double items in the final list
        //if j is nearest neighbour of i and i nearest neighbour of j
        pairDetection[2*idx+0] = idx | (1 << 31);
        pairDetection[2*idx+1] = j   | (1 << 31);

      //Another option is to store it like this, but this destroys the
      //info about pairs
      }
    }//if ds2 <=
  }//if j >= 0
#endif

}




#else
extern "C" __global__ void correct_particles(const int n_bodies,
                                             double tc,
                                             double2 *time,
                                             uint   *active_list,
                                             real4 *vel,
                                             real4 *acc0,
                                             real4 *acc1,
                                             real4 *pos,
                                             real4 *pPos,
                                             real4 *pVel) {
  const int bid =  blockIdx.y *  gridDim.x +  blockIdx.x;
  const int tid =  threadIdx.y * blockDim.x + threadIdx.x;
  const int dim =  blockDim.x * blockDim.y;

  int idx = bid * dim + tid;
  if (idx >= n_bodies) return;

  //Check if particle is set to active during approx grav
  if (active_list[idx] != 1) return;


  float4 v  = vel [idx];
  float4 a0 = acc0[idx];
  float4 a1 = acc1[idx];
  float  tb = time[idx].x;
//   float  dt = time[idx].y;

  float dt_cb = (float)(tc - tb);

  //Store the predicted position as the one to use
  pos[idx] = pPos[idx];

  //Correct the position
  v = pVel[idx];

  dt_cb *= 0.5f;
  v.x += (a1.x - a0.x)*dt_cb;
  v.y += (a1.y - a0.y)*dt_cb;
  v.z += (a1.z - a0.z)*dt_cb;

  //Store the corrected velocity, accelaration and the new time step info
  vel [idx] = v;
  acc0[idx] = a1;

//   time[idx] = (float2){tc, tc + dt};

}
#endif







#if 0
extern "C" __global__ void correct_particles(const int n_bodies,
                                             float tc,
                                             float2 *time,
                                             uint   *active_list,
                                             real4 *vel,
                                             real4 *acc0,
                                             real4 *acc1) {
  const int bid =  blockIdx.y *  gridDim.x +  blockIdx.x;
  const int tid =  threadIdx.y * blockDim.x + threadIdx.x;
  const int dim =  blockDim.x * blockDim.y;

  int idx = bid * dim + tid;
  if (idx >= n_bodies) return;

  //Check if particle is set to active during approx grav
  if (active_list[idx] != 1) return;


  float4 v  = vel [idx];
  float4 a0 = acc0[idx];
  float4 a1 = acc1[idx];
  float  tb = time[idx].x;

  float dt_cb = tc - tb;

  v.x -= a0.x * dt_cb;
  v.y -= a0.y * dt_cb;
  v.z -= a0.z * dt_cb;

  dt_cb *= 0.5f;
  v.x += (a0.x + a1.x)*dt_cb;
  v.y += (a0.y + a1.y)*dt_cb;
  v.z += (a0.z + a1.z)*dt_cb;

  vel [idx] = v;
  acc0[idx] = a1;
}
#endif

extern "C"  __global__ void compute_dt(const int n_bodies,
                                       double    tc,
                                       float    eta,
                                       int      dt_limit,
                                       float    eps2,
                                       double2   *time,
                                       real4    *vel,
                                       int      *ngb,
                                       real4    *bodies_pos,
                                       real4    *bodies_acc,
                                       uint     *active_list,
                                       float    timeStep){
  const int bid =  blockIdx.y *  gridDim.x +  blockIdx.x;
  const int tid =  threadIdx.y * blockDim.x + threadIdx.x;
  const int dim =  blockDim.x * blockDim.y;

  int idx = bid * dim + tid;
  if (idx >= n_bodies) return;

  //Check if particle is set to active during approx grav
  if (active_list[idx] != 1) return;

  int j = ngb[idx];

  float4 ri, rj;
  float4 vi, vj;
  float4 ai, aj;

  float ds2, mi, mj;
  ri = bodies_pos[idx];
  mi = ri.w;
  vi = vel[idx];
  ai = bodies_acc[idx];
  int j1, j2;

  if (j >= 0) {
    rj = bodies_pos[j];
    float3 dr = {ri.x - rj.x,
                 ri.y - rj.y,
                 ri.z - rj.z};
    ds2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
  } else  {
    j1 = max(0, idx - 1);
    rj = bodies_pos[j1];
    mj = rj.w;
    float3 dr = {ri.x - rj.x,
                 ri.y - rj.y,
                 ri.z - rj.z};
    if (idx != j1)  ds2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
    else            ds2 = 1.0e10f;

    j2 = min(n_bodies-1, idx + 1);
    rj = bodies_pos[j2];
    dr = (float3){ri.x - rj.x,
                  ri.y - rj.y,
                  ri.z - rj.z};
    if (idx != j2) {
      if (dr.x*dr.x + dr.y*dr.y + dr.z*dr.z < ds2) {
        ds2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
        j = j2;
	mj = rj.w;
      } else {
        j = j1;
      };
    } else {
      j = j1;
    }
  }

  //Add softening to the distance between the chosen particles
  ds2 += eps2;


  vj = vel[j];
  aj = bodies_acc[j];

  const float3 vda = make_float3(ai.x - aj.x,
                           ai.y - aj.y,
                           ai.z - aj.z);
  const float3 vdv = make_float3(vi.x - vj.x,
                           vi.y - vj.y,
                           vi.z - vj.z);

  const float vs2 = vdv.x*vdv.x + vdv.y*vdv.y  + vdv.z*vdv.z;

  //Compute the minimum crossing time
  const float mct = (ds2*ds2) / (vs2*vs2);

  //Free fall time
  float da2 = vda.x*vda.x + vda.y*vda.y + vda.z*vda.z;
  float mij = mi + mj; //Sum masses
  da2      *= (mij*mij);

  const float fft = (ds2 / da2);

  //Time step is minimum of the free fall time and minimum crossing time
  float dt_est = sqrt(sqrt(min(mct, fft)));

  //Make it a power of 2
  float dt_param = eta; //eta
 // float dt_param = 1.0; //eta
  float dt = dt_est*dt_param;

  int power = -(int)__log2f(dt) + 1;
  power     = max(power, dt_limit);


  int count = 0;
  dt = 1.0f/(1 << power);
  while(fmodf(tc, dt) != 0.0f)
  {
	  dt *= 0.5f;      // could be slow!
	  count++;
	  if(count > 30)
	  {
		  dt = timeStep;
		  break;
	  }
  }

  //if(dt < 1./16384) dt = 1./16384;
  //if(dt < 1./1048576) dt = 1./1048576;


  time[idx].x = tc;

#ifdef ADAPTIVE_TIMESTEP

  //Prevent a time-step smaller than specified through the interface
  if(dt < timeStep)
    dt = timeStep;


  time[idx].y = tc + (double)dt;
#else
  time[idx].y = tc + timeStep;
#endif
//  if(idx % 1000 == 0)
//    time[idx].y = tc + 1./2048 ;
//  else
//    time[idx].y = tc + timeStep;




#if 0

  ds2 = ds2*__powf(10.0f, 0.666667f) + eps2;
//   ds2 += eps2;
  vj = vel[j];
  aj = bodies_acc[j];

  float3 vda = {ai.x - aj.x,
                ai.y - aj.y,
                ai.z - aj.z};
  float3 vdv = {vi.x - vj.x,
                vi.y - vj.y,
                vi.z - vj.z};
  float da = sqrtf(vda.x*vda.x + vda.y*vda.y + vda.z*vda.z);
  float dv = sqrtf(vdv.x*vdv.x + vdv.y*vdv.y + vdv.z*vdv.z);
  float ds = sqrtf(ds2);

  float dt = eta * dv/da*(sqrt(2*da*ds/(dv*dv) + 1) - 1);

  int power = -(int)__log2f(dt) + 1;
  power     = max(power, dt_limit);

  dt = 1.0f/(1 << power);
  while(fmodf(tc, dt) != 0.0f) dt *= 0.5f;      // could be slow!

//  dt = 0.015625;
  dt = 1.0f/(1 << 8);
  dt = 1.0f/(1 << 6);
  dt = 1.0f/(1 << 7);
  dt = timeStep;
  time[idx].x = tc;
  //time[idx].y = tc + dt;
  time[idx].y = tc + dt;
#endif
}


//Reduce function to get the energy of the system in single precision
__device__ void compute_energyD(const int n_bodies,
                                            real4 *pos,
                                            real4 *vel,
                                            real4 *acc,
                                            float2 *energy, volatile float *shDataKin) {

  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  const int blockSize   = blockDim.x;
  unsigned int tid      = threadIdx.x;
  unsigned int i        = blockIdx.x*(blockSize*2) + threadIdx.x;
  unsigned int gridSize = blockSize*2*gridDim.x;

  volatile float *shDataPot = (float*)&shDataKin [blockSize];
  float eKin, ePot;
  shDataKin[tid] = eKin = 0;   //Stores Ekin
  shDataPot[tid] = ePot = 0;   //Stores Epot

  real4 temp;
  // we reduce multiple elements per thread.  The number is determined by the
  // number of active thread blocks (via gridSize).  More blocks will result
  // in a larger gridSize and therefore fewer elements per thread
  while (i < n_bodies) {
    if (i             < n_bodies)
    {
      //Ekin
      temp  = vel[i];
      eKin += pos[i].w*0.5*(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);

      //Epot
      ePot += pos[i].w*0.5*acc[i].w;
    }

    if (i + blockSize < n_bodies)
    {
      temp = vel[i + blockSize];
      eKin += pos[i + blockSize].w*0.5*(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);

      //Epot
      ePot += pos[i + blockSize].w*0.5*acc[i + blockSize].w;
    }

    i += gridSize;
  }
  shDataKin[tid] = eKin;
  shDataPot[tid] = ePot;

  __syncthreads();

  // do reduction in shared mem
  if (blockSize >= 512) { if (tid < 256) {
    shDataPot[tid] = ePot = ePot + shDataPot[tid + 256];
    shDataKin[tid] = eKin = eKin + shDataKin[tid + 256];   } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) {
    shDataPot[tid] = ePot = ePot + shDataPot[tid + 128];
    shDataKin[tid] = eKin = eKin + shDataKin[tid + 128];   } __syncthreads(); }
  if (blockSize >= 128) { if (tid <  64) {
    shDataPot[tid] = ePot = ePot + shDataPot[tid + 64];
    shDataKin[tid] = eKin = eKin + shDataKin[tid + 64];   } __syncthreads(); }


#ifndef __DEVICE_EMULATION__
  if (tid < 32)
#endif
    {
      if (blockSize >=  64) {shDataKin[tid] = eKin = eKin + shDataKin[tid + 32]; shDataPot[tid] = ePot = ePot + shDataPot[tid + 32];  EMUSYNC; }
      if (blockSize >=  32) {shDataKin[tid] = eKin = eKin + shDataKin[tid + 16]; shDataPot[tid] = ePot = ePot + shDataPot[tid + 16];  EMUSYNC; }
      if (blockSize >=  16) {shDataKin[tid] = eKin = eKin + shDataKin[tid +  8]; shDataPot[tid] = ePot = ePot + shDataPot[tid +  8];  EMUSYNC; }
      if (blockSize >=   8) {shDataKin[tid] = eKin = eKin + shDataKin[tid +  4]; shDataPot[tid] = ePot = ePot + shDataPot[tid +  4];  EMUSYNC; }
      if (blockSize >=   4) {shDataKin[tid] = eKin = eKin + shDataKin[tid +  2]; shDataPot[tid] = ePot = ePot + shDataPot[tid +  2];  EMUSYNC; }
      if (blockSize >=   2) {shDataKin[tid] = eKin = eKin + shDataKin[tid +  1]; shDataPot[tid] = ePot = ePot + shDataPot[tid +  1];  EMUSYNC; }
  }

  // write result for this block to global mem
  if (tid == 0) energy[blockIdx.x] = (float2){shDataKin[0], shDataPot[0] };
}

extern "C" __global__ void compute_energy(const int n_bodies,
                                            real4 *pos,
                                            real4 *vel,
                                            real4 *acc,
                                            float2 *energy) {

  extern __shared__ float shDataKin[];
  compute_energyD(n_bodies, pos, vel, acc, energy,shDataKin);
}




//Reduce function to get the energy of the system in double precision
__device__ void compute_energy_doubleD(const int n_bodies,
                                            real4 *pos,
                                            real4 *vel,
                                            real4 *acc,
                                            double2 *energy, volatile double *shDDataKin) {

  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  const int blockSize   = blockDim.x;
  unsigned int tid      = threadIdx.x;
  unsigned int i        = blockIdx.x*(blockSize*2) + threadIdx.x;
  unsigned int gridSize = blockSize*2*gridDim.x;

  volatile double *shDDataPot = (double*)&shDDataKin [blockSize];
  double eKin, ePot;
  shDDataKin[tid] = eKin = 0;   //Stores Ekin
  shDDataPot[tid] = ePot = 0;   //Stores Epot

  real4 temp;
  // we reduce multiple elements per thread.  The number is determined by the
  // number of active thread blocks (via gridSize).  More blocks will result
  // in a larger gridSize and therefore fewer elements per thread
  while (i < n_bodies) {
    if (i             < n_bodies)
    {
      //Ekin
      temp  = vel[i];
      eKin += pos[i].w*0.5*(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);

      //Epot
      ePot += pos[i].w*0.5*acc[i].w;
    }

    if (i + blockSize < n_bodies)
    {
      temp = vel[i + blockSize];
      eKin += pos[i + blockSize].w*0.5*(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);

      //Epot
      ePot += pos[i + blockSize].w*0.5*acc[i + blockSize].w;
    }

    i += gridSize;
  }
  shDDataKin[tid] = eKin;
  shDDataPot[tid] = ePot;

  __syncthreads();

  // do reduction in shared mem
  if (blockSize >= 512) { if (tid < 256) {
    shDDataPot[tid] = ePot = ePot + shDDataPot[tid + 256];
    shDDataKin[tid] = eKin = eKin + shDDataKin[tid + 256];   } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) {
    shDDataPot[tid] = ePot = ePot + shDDataPot[tid + 128];
    shDDataKin[tid] = eKin = eKin + shDDataKin[tid + 128];   } __syncthreads(); }
  if (blockSize >= 128) { if (tid <  64) {
    shDDataPot[tid] = ePot = ePot + shDDataPot[tid + 64];
    shDDataKin[tid] = eKin = eKin + shDDataKin[tid + 64];   } __syncthreads(); }


#ifndef __DEVICE_EMULATION__
  if (tid < 32)
#endif
    {
      if (blockSize >=  64) {shDDataKin[tid] = eKin = eKin + shDDataKin[tid + 32]; shDDataPot[tid] = ePot = ePot + shDDataPot[tid + 32];  EMUSYNC; }
      if (blockSize >=  32) {shDDataKin[tid] = eKin = eKin + shDDataKin[tid + 16]; shDDataPot[tid] = ePot = ePot + shDDataPot[tid + 16];  EMUSYNC; }
      if (blockSize >=  16) {shDDataKin[tid] = eKin = eKin + shDDataKin[tid +  8]; shDDataPot[tid] = ePot = ePot + shDDataPot[tid +  8];  EMUSYNC; }
      if (blockSize >=   8) {shDDataKin[tid] = eKin = eKin + shDDataKin[tid +  4]; shDDataPot[tid] = ePot = ePot + shDDataPot[tid +  4];  EMUSYNC; }
      if (blockSize >=   4) {shDDataKin[tid] = eKin = eKin + shDDataKin[tid +  2]; shDDataPot[tid] = ePot = ePot + shDDataPot[tid +  2];  EMUSYNC; }
      if (blockSize >=   2) {shDDataKin[tid] = eKin = eKin + shDDataKin[tid +  1]; shDDataPot[tid] = ePot = ePot + shDDataPot[tid +  1];  EMUSYNC; }
  }


  // write result for this block to global mem
  if (tid == 0) energy[blockIdx.x] = (double2){shDDataKin[0], shDDataPot[0] };

}


//Reduce function to get the energy of the system
extern "C" __global__ void compute_energy_double(const int n_bodies,
                                            real4 *pos,
                                            real4 *vel,
                                            real4 *acc,
                                            double2 *energy) {
  extern __shared__ double shDDataKin[];
  compute_energy_doubleD(n_bodies, pos, vel, acc, energy, shDDataKin);
}

extern "C"  __global__  void distanceCheck(const int n_bodies,
                              real4 *pos,
                              int   *ids,
                              real4 *out,
                              const int numberOfBH,
                              real4 *vel)
{
  const int bid =  blockIdx.y *  gridDim.x +  blockIdx.x;
  const int tid =  threadIdx.y * blockDim.x + threadIdx.x;
  const int dim =  blockDim.x * blockDim.y;

  int idx = bid * dim + tid;
  if (idx >= n_bodies) return;

  int partID = ids[idx];

  if(partID < numberOfBH)
  {
    real4 curPos = pos[idx];
    //curPos.w     = partID;
    out[partID*2+0]  = curPos;
    out[partID*2+1]  = vel[idx];
  }

}

