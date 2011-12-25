#define __out 

__device__ __constant__ float EPS2;
__device__ __constant__ float DT_TICK;

struct ds64
{
  union
  {
    float2 val;
    double dbl;
  };
  __device__ ds64() {}
  __device__ ds64(float x) : val(make_float2(x, x)) {}
  __device__ ds64 operator+=(const float x) 
  {
    const float vx = val.x + x;
    const float vy = val.y - ((vx - val.x) - x);
    val = make_float2(vx, vy);
    return *this;
  }
  __device__ double to_double() const { return (double)val.x + (double)val.y; }
};

template<class REAL>
struct cuvec3
{
  REAL x, y, z;
  __host__ __device__ cuvec3() {}
  __host__ __device__ cuvec3(const REAL v) : x(v), y(v), z(v) {}
  __host__ __device__ cuvec3(const REAL _x, const REAL _y, const REAL _z) : x(_x), y(_y), z(_z) {}
  
  __host__ __device__ cuvec3 operator=(const cuvec3<float> v) {x = v.x; y = v.y; z = v.z; return *this;};
  __host__ __device__ cuvec3 operator=(const cuvec3<double > v) {x = v.x; y = v.y; z = v.z; return *this;};
  

  __host__ __device__ REAL   operator*(const cuvec3<REAL> v) const {return       (x*v.x + y*v.y + z*v.z);}
  __host__ __device__ cuvec3 operator*(const        REAL  v) const {return cuvec3(x*v,   y*v,   z*v);}
//  __host__ __device__ cuvec3 operator+(const cuvec3<REAL> v) const {return cuvec3(x+v.x, y+v.y, z+v.z);}
  __host__ __device__ cuvec3 operator-(const cuvec3<REAL> v) const {return cuvec3(x-v.x, y-v.y, z-v.z);}
  __host__ __device__ cuvec3 operator%(const cuvec3<REAL> v) const {return cuvec3(x*v.y - y*v.x, y*v.z-z*v.y, z*v.x - x*v.z);}
  __host__ __device__ cuvec3 operator-() const {return cuvec3(-x, -y, -z);}
  
  __host__ __device__ cuvec3 operator+(const cuvec3<float> v) const {return cuvec3(x+v.x, y+v.y, z+v.z);}
  __host__ __device__ cuvec3 operator+(const cuvec3<double > v) const {return cuvec3(x+v.x, y+v.y, z+v.z);}


	__host__ __device__ cuvec3 operator += (const cuvec3<REAL> v)
  {
		*this = *this + v;
		return *this;
	}

	__host__ __device__ cuvec3 operator -= (const cuvec3<REAL> v)
  {
		*this = *this - v;
		return *this;
	}
	__host__ __device__ cuvec3 operator *= (const REAL s)
  {
		*this = *this * s;
		return *this;
	}
  __host__ __device__ friend cuvec3 operator * (const REAL s ,const cuvec3<REAL> v)
  {
    return v*s;
  }


  __host__ __device__ REAL norm2() const {return (*this)*(*this);};
};

typedef cuvec3<double > dcuvec3;
typedef cuvec3<float> fcuvec3;

__device__ float sqr(const float x)
{
  return x*x;
}

/****************************/
/****************************/
/****************************/

template<class T>
struct ADDOP 
{
  __device__ static inline T identity()           {return (T)(0);}
  __device__ static inline T apply(T a, T b)      {return (T)(a + b);};
  __device__ static inline T unapply(T a, T b)    {return (T)(a - b);};
  __device__ static inline T mask(bool flag, T b) {return (T)(-(int)(flag) & b);};
};


template<class OP, class T>
__device__ __forceinline__ T inclusive_scan_warp(volatile T *ptr, T mysum,  const unsigned int idx ) 
{
  const unsigned int lane = idx & 31;

  if (lane >=  1) ptr[idx] = mysum = OP::apply(ptr[idx -  1], mysum);
  if (lane >=  2) ptr[idx] = mysum = OP::apply(ptr[idx -  2], mysum);
  if (lane >=  4) ptr[idx] = mysum = OP::apply(ptr[idx -  4], mysum);
  if (lane >=  8) ptr[idx] = mysum = OP::apply(ptr[idx -  8], mysum);
  if (lane >= 16) ptr[idx] = mysum = OP::apply(ptr[idx - 16], mysum);

  return ptr[idx];
}

template<class OP, class T>
__device__ T inclusive_scan_block(volatile T *ptr, const unsigned int idx) 
{
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
__device__ T inclusive_scan_array(volatile T *ptr_global, const int N, const unsigned int idx) 
{
  T y = OP::identity();
  volatile T *ptr = ptr_global;

  for (int p = 0; p < N; p += blockDim.x) 
  {
    ptr = &ptr_global[p];
    inclusive_scan_block<OP, T>(ptr, idx);
    ptr[idx] = OP::apply(ptr[idx], y);
    __syncthreads();

    y = ptr[blockDim.x - 1];
    __syncthreads();
  }

  return y;
}

/****************************/
/****************************/
/****************************/

struct  dev_particle 
{
  dcuvec3 pos;              // 6
  fcuvec3 vel;              // 9
  fcuvec3 acc;              // 12
  fcuvec3 jrk;              // 15
  float mass;               // 16
  float h2;                 // 17
  unsigned int time;        // 18
  int id;                   // 19
  int iPad;                 // 20
  int iPadX[12];

  __host__ __device__ dev_particle() {}
  __host__ dev_particle(const regf4::Particle&);
};
#define PTCL_LEN (sizeof(dev_particle) / sizeof(float4))

struct dev_predictor
{
	fcuvec3 pos;        // 3
	fcuvec3 vel;        // 6
  union 
  {
    float   mass;       // 7
    float dt;
  };
  float   h2;         // 8
}; 
#define PRED_LEN (sizeof(dev_predictor) / sizeof(float4))

struct dev_force
{
  ds64 accx, accy, accz;  // 6
  fcuvec3 jrk;        // 9
  float h2;           // 10
  int   nngb;         // 11
  int   iPad;         // 12
  __device__ dev_force() : accx(0.0f), accy(0.0f), accz(0.0f), jrk(0.0f), nngb(0) {}
}; 

/********************************/
/********************************/
/********************************/

__global__ void dev_predict_ptcl(
    const int ni,
    const unsigned int tsys,
    const dev_particle   *ptcl_in,
    __out dev_predictor  *pred_out,
    __out float          *dt_out)
{
  const int id = blockIdx.x*blockDim.x + threadIdx.x;
  const int addr = id < ni ? id : ni-1;

  const dev_particle ip = ptcl_in[addr];

  dev_predictor ipred;

  const float dt  = DT_TICK*(tsys - ip.time);
  const float dt2 = dt*(1.0f/2.0f);
  const float dt3 = dt*(1.0f/3.0f);

  ipred.pos  = ip.pos + dt*(ip.vel + dt2*(ip.acc + dt3*ip.jrk));
  ipred.vel  = ip.vel + dt*(ip.acc + dt2* ip.jrk);
  ipred.mass = ip.mass;
  ipred.h2   = ip.h2;

  if (id < ni)
  {
    pred_out[addr] = ipred;
    dt_out  [addr] = dt;
  }
}

/********************************/
/********************************/
/********************************/

template<int NGB_PER_BLOCK>
__forceinline__ __device__ dev_force dev_regfij(
    const unsigned int jidx,
    const dev_predictor pi,
    const dev_predictor pj,
    __out	dev_force     fi,
    __out unsigned int *ngb_list)
{
  const fcuvec3 dr = pj.pos - pi.pos;
  const fcuvec3 dv = pj.vel - pi.vel;

  const float r2 = dr*dr;
  const float r2p = fminf(r2, (dr + pi.dt*dv).norm2());
  if (r2p < (pi.h2 + pj.h2)*0.5f)
  {
    if (pj.mass > 0.0f)
    {
      ngb_list[fi.nngb & (NGB_PER_BLOCK-1)] = jidx;
      fi.nngb += (r2 > 0.0f);
    }
  }
  else
  {
    const float rv    = dr*dv;
    const float rinv1 = rsqrt(r2 + EPS2);
    const float rinv2 = rinv1*rinv1;
    const float rinv3 = pj.mass*(rinv1*rinv2);
    const float alpha = rv*rinv2;

    const fcuvec3 Aij = rinv3*dr;
    const fcuvec3 Jij = rinv3*dv - Aij*(3.0f*alpha);
    fi.accx += Aij.x;
    fi.accy += Aij.y;
    fi.accz += Aij.z;
    fi.jrk  += Jij;
  }

  return fi;
}


/********************************/

template<int NTHREAD, int NJBLOCK, int NJBLOCK2, int NGB_PER_BLOCK>
__global__ void dev_regf(
    const int ni,
    const int nj_per_block,
    const int *active_list,
    const dev_predictor *pred_in,
    const float *dt_in,
    __out dev_force     *force_out,
    __out unsigned int  *ngb_out)
{
  __shared__ dev_predictor jpshared[NTHREAD];

  // compute iblock & jblock offset
  const int iblock =  blockIdx.x*NTHREAD;
  const int jblock =  blockIdx.y;
  const int tid    = threadIdx.x;

  // read i-particle into registers

  const int idx  = iblock + tid;
  const int addr = active_list[idx < ni ? idx : ni - 1];
  dev_predictor ipred = pred_in[addr];

  ipred.dt = dt_in[addr];

    // initialize i-particle's force

    dev_force iforce;

  // obtain beginning & end of j particles for this block

  const int jbeg = jblock*nj_per_block;
  const int jend = jbeg + nj_per_block;

  unsigned int *ingb_ptr = ngb_out + NGB_PER_BLOCK*(jblock + NJBLOCK*idx);


  for (int j = jbeg; j < jend; j += NTHREAD)
  {
#if 0
    jpshared[tid] = pred_in[j + tid];
#else
    float4 *src = (float4*)&pred_in[j];
    float4 *dst = (float4*)jpshared;
#pragma unroll
    for (int it = 0; it < PRED_LEN; it++)
    {
      dst[tid] = src[tid];
      dst += NTHREAD;
      src += NTHREAD;
    }
#endif
    __syncthreads();

    if (idx < ni)
    {
#pragma unroll 8
      for (int jj = 0; jj < NTHREAD; jj++)
        iforce = dev_regfij<NGB_PER_BLOCK>(j+jj, ipred, jpshared[jj], iforce, ingb_ptr);
    }
    __syncthreads();
  }

  if (idx < ni)
  {
    iforce.h2 = ipred.h2;
    force_out[jblock + idx*NJBLOCK2] = iforce;
  }
}



/********************************/
/********************************/
/********************************/
  template<class OP, class T, int NTHREAD>
__device__ T reduce_block(volatile T *ptr, T mySum, const unsigned int tid)
{
  ptr[tid] = mySum;
  __syncthreads();

  if (NTHREAD >= 512) { if (tid < 256) { ptr[tid] = mySum = OP::apply(mySum, ptr[tid+256]); } __syncthreads(); }
  if (NTHREAD >= 256) { if (tid < 128) { ptr[tid] = mySum = OP::apply(mySum, ptr[tid+128]); } __syncthreads(); }
  if (NTHREAD >= 128) { if (tid <  64) { ptr[tid] = mySum = OP::apply(mySum, ptr[tid+ 64]); } __syncthreads(); }

  if (tid < 32)
  {
    if (NTHREAD >= 64) ptr[tid] = mySum = OP::apply(mySum, ptr[tid+32]);
    if (NTHREAD >= 32) ptr[tid] = mySum = OP::apply(mySum, ptr[tid+16]);
    if (NTHREAD >= 16) ptr[tid] = mySum = OP::apply(mySum, ptr[tid+ 8]);
    if (NTHREAD >=  8) ptr[tid] = mySum = OP::apply(mySum, ptr[tid+ 4]);
    if (NTHREAD >=  4) ptr[tid] = mySum = OP::apply(mySum, ptr[tid+ 2]);
    if (NTHREAD >=  2) ptr[tid] = mySum = OP::apply(mySum, ptr[tid+ 1]);
  }
  __syncthreads();

  return ptr[0];
}


// here each particle is assigned to a single block...
// for 60 active blocks in dev_regf, and 64 threads the max efficiency is 60/64...
template<int NTHREAD, int NJBLOCK, int NJBLOCK2>
__global__ void dev_reduce_regf(
    const dev_force *force_in,
    __out int2      *ngb_offset,
    __out dev_force *force_out)
{
  // we use parallel prefix sum to obtain reduce forces 

  const int idx =  blockIdx.x;   // body id
  const int tid = threadIdx.x;   // block id

  __shared__ float shdata[2*NTHREAD];
  double *shdbl = (double*)shdata;

  dev_force iforce;
  if (tid < NJBLOCK)
    iforce = force_in[tid + idx*NJBLOCK2];

  iforce.accx.dbl = reduce_block<ADDOP<double>, double, NTHREAD>(shdbl, iforce.accx.to_double(), tid);
  iforce.accy.dbl = reduce_block<ADDOP<double>, double, NTHREAD>(shdbl, iforce.accy.to_double(), tid);
  iforce.accz.dbl = reduce_block<ADDOP<double>, double, NTHREAD>(shdbl, iforce.accz.to_double(), tid);

  iforce.jrk.x = reduce_block<ADDOP<float>, float, NTHREAD>(shdata, iforce.jrk.x, tid);
  iforce.jrk.y = reduce_block<ADDOP<float>, float, NTHREAD>(shdata, iforce.jrk.y, tid);
  iforce.jrk.z = reduce_block<ADDOP<float>, float, NTHREAD>(shdata, iforce.jrk.z, tid);

  int *shint = (int*)shdata;
  shint[tid] = iforce.nngb;
  inclusive_scan_block<ADDOP<int>, int>(shint, tid);
  const int nngb = shint[NTHREAD-1];

  /* #ngb in a block, memory offset */
#if 0
  if (idx == 0)
  {
    for (int t = 0; t < NTHREAD; t++)
    {
      __syncthreads();
      if (t == tid)
        printf(" nnbb= %d  offset= %d  addr= %d tid= %d NJBLOCK= %d\n",
            iforce.nngb, shint[tid] - iforce.nngb,
            idx + tid*NJBLOCK, tid, NJBLOCK);
    }

  }
#endif

  if (tid < NJBLOCK)
    ngb_offset[tid + idx*NJBLOCK] = (int2){iforce.nngb, shint[tid] - iforce.nngb};

  if (tid == 0)
  {
    iforce.nngb = nngb;
    force_out[idx] = iforce;
  }
}

/********************************/

template<int NTHREAD, int NJBLOCK, int NGB_PER_BLOCK, int NGB_MAX>
__global__ void dev_reduce_ngb(
    const int2         *ngb_offset,
    const unsigned int *ngb_in,
    __out unsigned int *ngb_out
    )
{
  const int idx =  blockIdx.x;   // body id
  const int tid = threadIdx.x;   

  for (int i = 0; i < NJBLOCK; i++)
  {
    const int2 ingb   = ngb_offset[i + idx*NJBLOCK];
    const int  nngb   = ingb.x;
    const int  offset = ingb.y;

    if (tid < nngb)
    {
#if 0
      if (idx == 0)
      {
        for (int t = 0; t < NTHREAD; t++)
        {
          __syncthreads();
          if (tid == t)
            printf("block= %d tid= %d: addr= %d   offset= %d  nngb= %d  newx= %d\n", i, tid, idx*NGB_MAX+offset+tid, offset, nngb ,offset+nngb);
        }
      }
#endif
      const int offset_tot = min(offset+tid, NGB_MAX-1);
      ngb_out[idx*NGB_MAX+offset_tot] = ngb_in[NGB_PER_BLOCK*(i + NJBLOCK*idx)+tid];
    }
  }
}

/********************************/

__global__ void dev_move_particles(
    const int nj,
    const int *addr_in,
    const dev_particle *ptcl_in,
    __out dev_particle *ptcl_out)
{
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= nj) return;

  const int addr = addr_in[idx];
  ptcl_out[addr] = ptcl_in[idx];
}

/********************************/

struct gpot_struct
{
  dcuvec3 pos;
  float mass;
};

template<int BLOCKSIZE>
__global__ void dev_compute_potential(
    const int ni,
    const int nj,
    const dev_particle *ptcl_in,
    __out float        *gpot_out)
{
  const int idx = blockDim.x*blockIdx.x + threadIdx.x;
  const int addr = idx < ni ? idx : ni - 1;

  const int tid = threadIdx.x;

  __shared__ gpot_struct shmem[BLOCKSIZE];

  ds64 gpot(0.0f);

  const dcuvec3 ipos = ptcl_in[addr].pos;

  for  (int j = 0; j < nj; j += BLOCKSIZE)
  {
    dev_particle pj = ptcl_in[j+tid];
    shmem[tid].pos  = pj.pos;
    shmem[tid].mass = pj.mass;
    __syncthreads();

#pragma unroll
    for (int jj = 0; jj < BLOCKSIZE; jj++)
    {
      const dcuvec3 jpos  = shmem[jj].pos;
      const float   jmass = shmem[jj].mass;
      const fcuvec3  dr  = fcuvec3(jpos.x - ipos.x, jpos.y - ipos.y, jpos.z - ipos.z);
      const float  r2  = dr*dr;
      const float rinv = (r2 > 0.0f) ? rsqrt(r2 + EPS2) : 0.0f;
      gpot += jmass * rinv;
    }
    __syncthreads();
  }

  if (idx < ni)
    gpot_out[idx] = -gpot.to_double();
}



