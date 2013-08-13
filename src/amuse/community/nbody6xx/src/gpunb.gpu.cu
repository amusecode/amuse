// #include <iostream>
#include <stdio.h>
// #include <vector>
#include <cmath>
#include <cassert>
#ifdef CUDA_5
#include <helper_cuda.h>
#define CUDA_SAFE_CALL checkCudaErrors
#else
#include <cutil.h>
#endif
#include "cuda_pointer.h"

#define NTHREAD 64 // 64, 96, 128 or 192
#define NJBLOCK 28 // 8800GTS/512 has 16
#define NIBLOCK 16 // 16 or 32 
#define NIMAX (NTHREAD * NIBLOCK) // 1024

#define NXREDUCE 32 // must be >NJBLOCK
#define NYREDUCE 8

#define NB_PER_BLOCK 256 // NNB per block
#define NB_BUF_SIZE (1<<20)

#define NAN_CHECK(val) assert((val) == (val));

typedef unsigned short uint16;

// template <class T>
// struct myvector{
// 	int num;
// 	T *val;
// 	myvector(){
// 		num = 0;
// 		val = NULL;
// 	}
// 	~myvector(){
// 		delete [] val;
// 	}
// 	void clear(){
// 		num = 0;
// 	}
// 	void reserve(size_t count){
// 		val = new T[count];
// 	}
// 	void free(){
// 		delete [] val;
// 	}
// 	void push_back(const T &t){
// 		val[num++] = t;
// 	}
// 	size_t size(){
// 		return num;
// 	}
// 	T &operator[](int i){
// 		return val[i];
// 	}
// };

#define PROFILE
#ifdef PROFILE
#include <sys/time.h>
static double get_wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.e-6 * tv.tv_usec;
}
#else
static double get_wtime(){
	return 0.0;
}
#endif

static double time_send, time_grav, time_out, time_nb;
static long long numInter;

struct Jparticle{
  float3 pos;
  float  mass;
  float3 vel;
  float  pad;
  Jparticle() {}
  Jparticle(double mj, double xj[3], double vj[3]){
    pos.x = xj[0];
    pos.y = xj[1];
    pos.z = xj[2];
    mass  = mj;
    vel.x = vj[0];
    vel.y = vj[1];
    vel.z = vj[2];

    NAN_CHECK(xj[0]);
    NAN_CHECK(xj[1]);
    NAN_CHECK(xj[2]);
    NAN_CHECK(mj);
    NAN_CHECK(vj[0]);
    NAN_CHECK(vj[1]);
    NAN_CHECK(vj[2]);
  }
};
struct Iparticle{
  float3 pos;
  float  h2;
  float3 vel;
  float  dtr;
  Iparticle() {}
  Iparticle(double h2i, double dtri,double xi[3], double vi[3]){
    pos.x = xi[0];
    pos.y = xi[1];
    pos.z = xi[2];
    h2    = h2i;
    vel.x = vi[0];
    vel.y = vi[1];
    vel.z = vi[2];
    dtr   = dtri;

    NAN_CHECK(xi[0]);
    NAN_CHECK(xi[1]);
    NAN_CHECK(xi[2]);
    NAN_CHECK(h2i);
    NAN_CHECK(vi[0]);
    NAN_CHECK(vi[1]);
    NAN_CHECK(vi[2]);
  }
};
struct Force{
	float3 acc;
	float  pot;
	float3 jrk;
	int    nnb;          //  8 words
  //	unsigned short  neib[NB_PER_BLOCK]; // 24 words
	// __device__  Force(){
	// 	acc.x = acc.y = acc.z = 0.f;
	// 	jrk.x = jrk.y = jrk.z = 0.f;
	// 	pot = 0.f;
	// 	nnb = 0;
	// }
  	__device__  void clear(){
		acc.x = acc.y = acc.z = 0.f;
		jrk.x = jrk.y = jrk.z = 0.f;
		pot = 0.f;
		nnb = 0;
	}
  	__device__ void operator+=(const Force &rhs){
		acc.x += rhs.acc.x;
		acc.y += rhs.acc.y;
		acc.z += rhs.acc.z;
#ifdef POTENTIAL
		pot   += rhs.pot;
#endif        
		jrk.x += rhs.jrk.x;
		jrk.y += rhs.jrk.y;
		jrk.z += rhs.jrk.z;
		if(nnb>=0 && rhs.nnb>=0){
			nnb += rhs.nnb;
		}else{
			nnb = -1;
		}
	}
#if __CUDA_ARCH__ >= 300
	__device__ void reduce_with(const int mask){
		acc.x += __shfl_xor(acc.x, mask);
		acc.y += __shfl_xor(acc.y, mask);
		acc.z += __shfl_xor(acc.z, mask);
#ifdef POTENTIAL        
		pot   += __shfl_xor(pot  , mask);
#endif        
		jrk.x += __shfl_xor(jrk.x, mask);
		jrk.y += __shfl_xor(jrk.y, mask);
		jrk.z += __shfl_xor(jrk.z, mask);
		int ntmp = __shfl_xor(nnb, mask);
		if(nnb>=0 && ntmp>=0){
			nnb += ntmp;
		}else{
			nnb = -1;
		}
	}
#endif
};

// __device__ float rsqrtfNR(float x){
// 	float y = rsqrtf(x);
// 	return (-0.5f * y) * (x*y*y - 3.0f);
// }

__device__ void h4_kernel(
		const int j,
		const Iparticle &ip, 
		const Jparticle &jp, 
		Force &fo,
        uint16 nblist[]){
	float dx = jp.pos.x - ip.pos.x;
	float dy = jp.pos.y - ip.pos.y;
	float dz = jp.pos.z - ip.pos.z;
	float dvx = jp.vel.x - ip.vel.x;
	float dvy = jp.vel.y - ip.vel.y;
	float dvz = jp.vel.z - ip.vel.z;

	float r2 = dx*dx + dy*dy + dz*dz;
    //Add Velocity criterion============================//
	float dxp = dx + ip.dtr * dvx;
	float dyp = dy + ip.dtr * dvy;
	float dzp = dz + ip.dtr * dvz;
	float r2p = dxp*dxp + dyp*dyp + dzp*dzp;
    //==================================================//
    
	float rv = dx*dvx + dy*dvy + dz*dvz;
	float rinv1 = rsqrtf(r2);
	if(min(r2,r2p) < ip.h2){
		// fo.neib[fo.nnb++ % NB_PER_BLOCK] = j;
		nblist[fo.nnb & (NB_PER_BLOCK-1)] = (uint16)j;
		fo.nnb++;
		rinv1 = 0.f;
	}
	float rinv2 = rinv1 * rinv1;
	float mrinv1 = jp.mass * rinv1;
	float mrinv3 = mrinv1 * rinv2;
	rv *= -3.f * rinv2;
	
#ifdef POTENTIAL
	fo.pot += mrinv1;
#endif
	fo.acc.x += mrinv3 * dx;
	fo.acc.y += mrinv3 * dy;
	fo.acc.z += mrinv3 * dz;
	// fo.acc.z += 1.0;
	fo.jrk.x += mrinv3 * (dvx + rv * dx);
	fo.jrk.y += mrinv3 * (dvy + rv * dy);
	fo.jrk.z += mrinv3 * (dvz + rv * dz);
}
__global__ void h4_gravity(
                           const int nbody,
                           const Iparticle ipbuf[],
                           const Jparticle jpbuf[],
                           Force fobuf[][NJBLOCK],
                           uint16 nbbuf[][NJBLOCK][NB_PER_BLOCK]){
  int ibid = blockIdx.x;
  int jbid = blockIdx.y;
  int tid = threadIdx.x;
  int iaddr = tid + NTHREAD * ibid;
  int jstart = (nbody * (jbid  )) / NJBLOCK;
  int jend   = (nbody * (jbid+1)) / NJBLOCK;

  Iparticle ip = ipbuf[iaddr];
  Force fo;
  fo.clear();
  uint16 *nblist = nbbuf[iaddr][jbid];
#if __CUDA_ARCH__ >= 300 // just some trial
	for(int j=jstart; j<jend; j+=32){
		__shared__ Jparticle jpshare[32];
		__syncthreads();
		float4 *src = (float4 *)&jpbuf[j];
		float4 *dst = (float4 *)jpshare;
		dst[tid] = src[tid];
		__syncthreads();
		if(jend-j < 32){
#pragma unroll 4
			for(int jj=0; jj<jend-j; jj++){
				const Jparticle jp = jpshare[jj];
				// const Jparticle jp( (float4 *)jpshare + 2*jj);
				h4_kernel(j-jstart+jj, ip, jp, fo, nblist);
			}
		}else{
#pragma unroll 8
			for(int jj=0; jj<32; jj++){
				const Jparticle jp = jpshare[jj];
				// const Jparticle jp( (float4 *)jpshare + 2*jj);
				h4_kernel(j-jstart+jj, ip, jp, fo, nblist);
			}
		}
	}
#else
  for(int j=jstart; j<jend; j+=NTHREAD){
    __shared__ Jparticle jpshare[NTHREAD];
    __syncthreads();

    float4 *src = (float4 *)&jpbuf[j];
    float4 *dst = (float4 *)jpshare;
    dst[        tid] = src[        tid];
    dst[NTHREAD+tid] = src[NTHREAD+tid];

    __syncthreads();

    if(jend-j < NTHREAD){
#pragma unroll 4
      for(int jj=0; jj<jend-j; jj++){
        Jparticle jp = jpshare[jj];
        h4_kernel(j-jstart+jj, ip, jp, fo, nblist);
      }
    }else{
#pragma unroll 8
      for(int jj=0; jj<NTHREAD; jj++){
        Jparticle jp = jpshare[jj];
        h4_kernel(j-jstart+jj, ip, jp, fo, nblist);
      }
    }
  }
#endif
  if(fo.nnb > NB_PER_BLOCK) fo.nnb = -1;
  fobuf[iaddr][jbid] = fo;
}

#if __CUDA_ARCH__ >= 300
__device__ void warp_reduce_int(int inp, int *out){
	inp += __shfl_xor(inp, 1);
	inp += __shfl_xor(inp, 2);
	inp += __shfl_xor(inp, 4);
	inp += __shfl_xor(inp, 8);
# if NXREDUCE==32
	inp += __shfl_xor(inp, 16);
# endif
	*out = inp;
}

__device__ void warp_reduce_float8(float4 inp1, float4 inp2, float *out){
	const int tid = threadIdx.x;
	float4 tmp4L = (4&tid) ? inp2 : inp1;
	float4 tmp4R = (4&tid) ? inp1 : inp2;
	tmp4L.x += __shfl_xor(tmp4R.x, 4);
	tmp4L.y += __shfl_xor(tmp4R.y, 4);
	tmp4L.z += __shfl_xor(tmp4R.z, 4);
	tmp4L.w += __shfl_xor(tmp4R.w, 4);
	float4 tmp4;
	tmp4.x = (2&tid) ? tmp4L.z : tmp4L.x;
	tmp4.y = (2&tid) ? tmp4L.w : tmp4L.y;
	tmp4.z = (2&tid) ? tmp4L.x : tmp4L.z;
	tmp4.w = (2&tid) ? tmp4L.y : tmp4L.w;
	tmp4.x += __shfl_xor(tmp4.z, 2);
	tmp4.y += __shfl_xor(tmp4.w, 2);
	float2 tmp2;
	tmp2.x = (1&tid) ? tmp4.y : tmp4.x;
	tmp2.y = (1&tid) ? tmp4.x : tmp4.y;
	tmp2.x += __shfl_xor(tmp2.y, 1);

	tmp2.x += __shfl_xor(tmp2.x, 8);
# if NXREDUCE==32
	tmp2.x += __shfl_xor(tmp2.x, 16);
# endif
	if(tid < 8){
		out[tid] = tmp2.x;
	}
}
#endif

__global__ void force_reduce_kernel(
		const int ni,
		const Force fpart[][NJBLOCK],
        Force ftot []){
	const int xid = threadIdx.x;
	const int yid = threadIdx.y;
	const int bid = blockIdx.x;
	const int iaddr = yid + blockDim.y * bid;

#if __CUDA_ARCH__ >= 300
	Force f;
	if(xid < NJBLOCK){
		f = fpart[iaddr][xid];
	}else{
		f.clear();
	}
	if(iaddr < ni){
		const float4 tmp1 = make_float4(f.acc.x, f.acc.y, f.acc.z, f.pot);
		const float4 tmp2 = make_float4(f.jrk.x, f.jrk.y, f.jrk.z, 0.0f);
		const int    itmp = f.nnb;
		float *dst  = (float *)(ftot + iaddr);
		int   *idst = (int *)(dst + 7);
		warp_reduce_float8(tmp1, tmp2, dst);
		warp_reduce_int(itmp, idst);
	}
#else
	__shared__ Force fshare[NYREDUCE][NXREDUCE];
	if(xid < NJBLOCK){
		fshare[yid][xid] = fpart[iaddr][xid];
	}else{
		fshare[yid][xid].clear();
	}
	Force *fs = fshare[yid];
#if NXREDUCE==32
	if(xid < 16) fs[xid] += fs[xid + 16];
#endif
	if(xid < 8) fs[xid] += fs[xid + 8];
	if(xid < 4) fs[xid] += fs[xid + 4];
	if(xid < 2) fs[xid] += fs[xid + 2];
	if(xid < 1) fs[xid] += fs[xid + 1];
	
	if(iaddr < ni){
		ftot[iaddr] = fs[0];
	}
#endif
}

__global__ void gather_nb_kernel(
		const int    ni,
		const int    nj,
		const Force  fpart[][NJBLOCK],
		const Force  ftot [],
		const int    nboff[],
		const uint16 nbpart[][NJBLOCK][NB_PER_BLOCK],
        int  nblist[])
{
	const int xid = threadIdx.x;
	const int yid = threadIdx.y;
	const int bid = blockIdx.x;
	const int iaddr = yid + blockDim.y * bid;
	if(iaddr >= ni) return;
	if(ftot[iaddr].nnb < 0) return;

	const int mynnb = (xid < NJBLOCK) ? fpart[iaddr][xid].nnb
	                                  : 0;

	// now performe prefix sum
#if __CUDA_ARCH__ >= 300
	int ix = mynnb;
#pragma unroll
	for(int ioff=1; ioff<NXREDUCE; ioff*=2){
		int iy = __shfl_up(ix, ioff);
		if(xid>=ioff) ix += iy;
	}
	int iz = __shfl_up(ix, 1);
	const int off = (xid == 0) ? 0 : iz;
#else
	__shared__ int ishare[NYREDUCE][NXREDUCE];
	ishare[yid][xid] = mynnb;
	volatile int *ish = ishare[yid];
	if(xid>=1)  ish[xid] += ish[xid-1];
	if(xid>=2)  ish[xid] += ish[xid-2];
	if(xid>=4)  ish[xid] += ish[xid-4];
	if(xid>=8)  ish[xid] += ish[xid-8];
#if NXREDUCE==32
	if(xid>=16)  ish[xid] += ish[xid-16];
#endif
	const int off = (xid == 0) ? 0 
	                           : ish[xid-1];
#endif
	int *nbdst = nblist + nboff[iaddr] + off;

	const int jstart = (nj * xid) / NJBLOCK;
	if(xid < NJBLOCK){
		for(int k=0; k<mynnb; k++){
			const int nbid = jstart + int(nbpart[iaddr][xid][k]);
			// const int nbid = iaddr * 1000 + k;
			nbdst[k] = nbid;
		}
	}
}

static cudaPointer <Jparticle> jpbuf;
static cudaPointer <Iparticle> ipbuf;
static cudaPointer <Force[NJBLOCK]> fopart;
static cudaPointer <Force> fobuf;
static cudaPointer <uint16[NJBLOCK][NB_PER_BLOCK]>nbpart;
static cudaPointer <int> nblist;
static cudaPointer <int> nboff;

//static myvector<int> nblist;
static int nbody, nbodymax;
static int devid, numGPU;
static bool is_open = false;
static bool devinit = false;
// static int *nblist;

void GPUNB_devinit(int irank){
  if(devinit) return;
  
  cudaGetDeviceCount(&numGPU);
  assert(numGPU > 0);
  char *gpu_list = getenv("GPU_LIST");
  if(gpu_list)
  {
    numGPU = 0;
    char *p = strtok(gpu_list, " ");
    if (p) {
      devid = atoi(p);
      numGPU++;
    }
    assert(numGPU > 0);
  }else{
    devid=irank%numGPU;
  }
  cudaSetDevice(devid);

#ifdef PROFILE
  fprintf(stderr, "***********************\n");
  fprintf(stderr, "Initializing NBODY6/GPU library\n");
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, devid);
  fprintf(stderr, "#GPU %d; device: %d %s\n", numGPU, devid, prop.name);
  fprintf(stderr, "***********************\n");
#endif
  devinit = true;
}

void GPUNB_open(int nbmax,int irank){
	time_send = time_grav = time_nb = time_out = 0.0;
	numInter = 0;
    
    //select GPU========================================//
    GPUNB_devinit(irank);
    
    if(is_open){
      fprintf(stderr, "gpunb: it is already open\n");
      return;
	}
	is_open = true;

    //==================================================//
    // CUT_DEVICE_INIT();
	// size_t jpsize = nbmax * sizeof(Jparticle);
	// size_t ipsize = NIMAX * sizeof(Iparticle);
	// size_t fosize = NIBLOCK * NJBLOCK * NTHREAD * sizeof(Force);
	// cudaMallocHost((void **)&jp_host, jpsize);
	// jpsize += NTHREAD * sizeof(Jparticle);
	// cudaMalloc    ((void **)&jp_dev , jpsize);
	// cudaMallocHost((void **)&ip_host, ipsize);
	// cudaMalloc    ((void **)&ip_dev , ipsize);
	// cudaMallocHost((void **)&fo_host, fosize);
	// cudaMalloc    ((void **)&fo_dev , fosize);
	jpbuf.allocate(nbmax + NTHREAD);
	ipbuf.allocate(NIMAX);
    fopart.allocate(NIMAX);
	fobuf.allocate(NIMAX);
    nbpart.allocate(NIMAX);
    nblist.allocate(NB_BUF_SIZE);
    nboff.allocate(NIMAX+1);
	nbodymax = nbmax;

    //    nblist.reserve(nbmax);
#ifdef PROFILE
	fprintf(stderr, "RANK: %d ******************\n",irank);
	fprintf(stderr, "Opened NBODY6/GPU library\n");
	fprintf(stderr, "nbmax = %d\n", nbmax);
	fprintf(stderr, "***********************\n");
#endif
}

void GPUNB_close(){
  if(!is_open){
		fprintf(stderr, "gpunb: it is already close\n");
		return;
	}
	is_open = false;
    
	// cudaFreeHost(jp_host);
	// cudaFree    (jp_dev);
	// cudaFreeHost(ip_host);
	// cudaFree    (ip_dev);
	// cudaFreeHost(fo_host);
	// cudaFree    (fo_dev);
	jpbuf.free();
	ipbuf.free();
    fopart.free();
	fobuf.free();
    nbpart.free();
    nblist.free();
    nboff.free();
	nbodymax = 0;

#ifdef PROFILE
	fprintf(stderr, "Closed NBODY6/GPU library\n");
	fprintf(stderr, "%d*********************\n",devid);
	fprintf(stderr, "time send : %f sec\n", time_send);
	fprintf(stderr, "time grav : %f sec\n", time_grav);
    fprintf(stderr, "time nb   : %f sec\n", time_nb);
    fprintf(stderr, "time out  : %f sec\n", time_out);
	fprintf(stderr, "%f Gflops (gravity part only)\n", 60.e-9 * numInter / time_grav);
	fprintf(stderr, "***********************\n");
#endif
}

void GPUNB_send(
		int nj,
		double mj[],
		double xj[][3],
		double vj[][3]){
	time_send -= get_wtime();
	nbody = nj;
	assert(nbody <= nbodymax);
    //    time_send -= get_wtime();
	for(int j=0; j<nj; j++){
		// jp_host[j] = Jparticle(mj[j], xj[j], vj[j]);
		jpbuf[j] = Jparticle(mj[j], xj[j], vj[j]);
	}
	// size_t jpsize = nj * sizeof(Jparticle);
	// cudaMemcpy(jp_dev, jp_host, jpsize, cudaMemcpyHostToDevice);
	jpbuf.htod(nj);
	time_send += get_wtime();
}

void GPUNB_regf(
                int ni,
                double h2[],
                double dtr[],
                double xi[][3],
                double vi[][3],
                double acc[][3],
                double jrk[][3],
                double pot[],
                int lmax,
                int nnbmax,
                int *listbase){
  assert(is_open);

  time_grav -= get_wtime();
  numInter += ni * nbody;
  assert(0 < ni && ni <= NIMAX);
   
/*        printf(" ni lm %d %d %d \t %e %e %e\n",ni, lmax, nnbmax, h2[0], xi[0][0], vi[0][0]);*/
  for(int i=0; i<ni; i++){
    // ip_host[i] = Iparticle(h2[i], xi[i], vi[i]);
    ipbuf[i] = Iparticle(h2[i],dtr[i], xi[i], vi[i]);
  }
  // set i-particles
  // size_t ipsize = ni * sizeof(Iparticle);
  // cudaMemcpy(ip_dev, ip_host, ipsize, cudaMemcpyHostToDevice);
  ipbuf.htod(ni);
  
  // gravity kernel
  int niblock = 1 + (ni-1) / NTHREAD;
  dim3 grid(niblock, NJBLOCK, 1);
  dim3 threads(NTHREAD, 1, 1);
  // h4_gravity <<< grid, threads >>> 
  //	(nbody, ip_dev, jp_dev, fo_dev);
  h4_gravity <<< grid, threads >>> 
    (nbody, ipbuf, jpbuf, fopart, nbpart);

  const int ni8 = 1 + (ni-1) / NYREDUCE;
  dim3 rgrid   (ni8, 1, 1);
  dim3 rthreads(NXREDUCE, NYREDUCE, 1);
  force_reduce_kernel <<< rgrid, rthreads >>>
    (ni, fopart, fobuf);

  // recieve force
  // size_t fosize = ni * NJBLOCK * sizeof(Force);
  // cudaMemcpy(fo_host, fo_dev, fosize, cudaMemcpyDeviceToHost);
  fobuf.dtoh(ni);

  double wt = get_wtime();
  time_grav += wt;
  time_nb -= wt;

  // now make prefix sum
  int nbsum = 0;
  for(int i=0; i<ni; i++){
    nboff[i] = nbsum;
    const int nnb = fobuf[i].nnb;
    if(nnb >= 0) nbsum += nnb;
  }
  assert(nbsum <= NB_BUF_SIZE);
  nboff.htod(ni);
  
  gather_nb_kernel <<< rgrid, rthreads>>>
    (ni, nbody, fopart, fobuf, nboff, nbpart, nblist);
  nblist.dtoh(nbsum);

  wt = get_wtime();
  time_nb += get_wtime();
  time_out -= get_wtime();
  
  // out data
  for(int i=0; i<ni; i++){
    Force &fo = fobuf[i];
    acc[i][0] = fo.acc.x;
    acc[i][1] = fo.acc.y;
    acc[i][2] = fo.acc.z;
    jrk[i][0] = fo.jrk.x;
    jrk[i][1] = fo.jrk.y;
    jrk[i][2] = fo.jrk.z;
    //    fprintf(stderr, "%f %f %f %f %f %f\n", acc[i][0], acc[i][1], acc[i][2], jrk[i][0], jrk[i][1], jrk[i][2]);
        //        exit(0);
#ifdef POTENTIAL
    pot[i] = fo.pot;
#endif
    int *nnbp = listbase + lmax * i;
    int *nblistp = nnbp + 1;
    if(fo.nnb >=0 && fo.nnb <= nnbmax){
      *nnbp = fo.nnb;
      //      fprintf(stderr, "nnb %d\n", fo.nnb);
      const int off = nboff[i];
      for(int k=0; k<fo.nnb; k++) nblistp[k]=nblist[off + k];
    }
    else *nnbp = fo.nnb ? -abs(fo.nnb) : -1;
  }
  time_out += get_wtime();
}

extern "C" {
  void gpunb_devinit_ (int *irank){
    GPUNB_devinit(*irank);
  }
  void gpunb_open_(int *nbmax, int *irank){
    GPUNB_open(*nbmax, *irank);
  }
  void gpunb_close_(){
    GPUNB_close();
  }
	void gpunb_send_(
			int *nj,
			double mj[],
			double xj[][3],
			double vj[][3]){
		GPUNB_send(*nj, mj, xj, vj);
	}
	void gpunb_regf_(
			int *ni,
			double h2[],
			double dtr[],
			double xi[][3],
			double vi[][3],
			double acc[][3],
			double jrk[][3],
			double pot[],
			int *lmax,
			int *nnbmax,
			int *list){ // list[][lmax]
      GPUNB_regf(*ni, h2, dtr, xi, vi, acc, jrk, pot, *lmax, *nnbmax, list);
	}
}
