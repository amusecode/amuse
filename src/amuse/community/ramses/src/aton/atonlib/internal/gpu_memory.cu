// GPU memory functions.

#include <stdio.h>
#include <stdlib.h>

#include "aton_fortran.h"
#include "gpu.h"

#define NCELLS3 (NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2)
#define NBUFF (max(NCELLX, max(NCELLY, NCELLZ)))


// GPU global array definitions.
double *cuegy, *cuegy_new;
double *cuflx, *cuflx_new;
double *cusrc0;
int *cusrc0pos;
double *cutemperature;
double *cuxion;
double *cudensity;
double *cu_photon_source;
double *cu_boundary_values;


extern "C" void aton_get_grid_size_(int *ncellx, int *ncelly, int *ncellz, int *nbound) {
  *ncellx = NCELLX;
  *ncelly = NCELLY;
  *ncellz = NCELLZ;
  *nbound = NBOUND;
}

namespace aton {
  int get_boundary_buffer_size() {
    return 4*NBUFF*NBUFF;
  }
}

extern "C" void aton_gpu_malloc_(int* source_count) {
  cudaMalloc((void**)&cuegy, ((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double));
  cudaMalloc((void**)&cuflx, ((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2)*3)*sizeof(double));

  cudaMalloc((void**)&cuegy_new, ((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double)); 
  cudaMalloc((void**)&cuflx_new, ((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2)*3)*sizeof(double));
 
  cudaMalloc((void**)&cusrc0, (*source_count)*sizeof(double));
  cudaMalloc((void**)&cusrc0pos, 3*(*source_count)*sizeof(int));

  cudaMalloc((void**)&cuxion, ((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double));
  cudaMalloc((void**)&cudensity, ((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double));
  cudaMalloc((void**)&cutemperature, ((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double));

  cudaMalloc((void**)&cu_boundary_values, 6*4*NBUFF*NBUFF * sizeof(double));

  cudaMalloc((void**)&cu_photon_source, ((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double));

  cudaThreadSynchronize();

  fprintf(stderr, "Allocation error?: %s\n", cudaGetErrorString(cudaGetLastError()));
}

// Packs a single boundary face.
//
// i1: axis to keep fixed (0, 1 or 2)
// i2, i3: axes to loop over (0, 1, or 2). i1, i2, and i3 should be distinct.
// fixed: fixed coordinate along axis i1.
//
__global__ void cuPackBoundary(int i1, int i2, int i3,
			       int fixed,
			       const double* cuegy, const double* cuflx,
			       double *cu_boundary_values) {
  int bound_idx = threadIdx.x + NBUFF*blockIdx.x;

  int pos[3];
  pos[i1] = fixed;
  pos[i2] = threadIdx.x;
  pos[i3] = blockIdx.x;
  int idx = ((pos[0]+NBOUND) +
	     (pos[1]+NBOUND) * (NCELLX + 2*NBOUND) +
	     (pos[2]+NBOUND) * (NCELLX + 2*NBOUND) * (NCELLY + 2*NBOUND));

  cu_boundary_values[bound_idx + 0*NBUFF*NBUFF] = cuegy[idx];
  cu_boundary_values[bound_idx + 1*NBUFF*NBUFF] = cuflx[idx + 0*NCELLS3];
  cu_boundary_values[bound_idx + 2*NBUFF*NBUFF] = cuflx[idx + 1*NCELLS3];
  cu_boundary_values[bound_idx + 3*NBUFF*NBUFF] = cuflx[idx + 2*NCELLS3];
}

// Reverse of cuPackBoundary.
__global__ void cuUnpackBoundary(int i1, int i2, int i3,
			       int fixed,
			       double* cuegy, double* cuflx,
			       const double *cu_boundary_values) {
  int bound_idx = threadIdx.x + NBUFF*blockIdx.x;

  int pos[3];
  pos[i1] = fixed;
  pos[i2] = threadIdx.x;
  pos[i3] = blockIdx.x;
  int idx = ((pos[0]+NBOUND) +
	     (pos[1]+NBOUND) * (NCELLX + 2*NBOUND) +
	     (pos[2]+NBOUND) * (NCELLX + 2*NBOUND) * (NCELLY + 2*NBOUND));

  cuegy[idx] = cu_boundary_values[bound_idx + 0*NBUFF*NBUFF];
  cuflx[idx + 0*NCELLS3] = cu_boundary_values[bound_idx + 1*NBUFF*NBUFF];
  cuflx[idx + 1*NCELLS3] = cu_boundary_values[bound_idx + 2*NBUFF*NBUFF];
  cuflx[idx + 2*NCELLS3] = cu_boundary_values[bound_idx + 3*NBUFF*NBUFF];
}

extern "C" void aton_gpu_to_cpu_boundary_(double *boundary_values) {
  dim3 xdim(NCELLX);
  dim3 ydim(NCELLY);
  dim3 zdim(NCELLZ);

  int stride = 4*NBUFF*NBUFF;

  cuPackBoundary<<<zdim,ydim>>>(0, 1, 2, 0, cuegy, cuflx,
				&cu_boundary_values[0*stride]);
  cuPackBoundary<<<zdim,ydim>>>(0, 1, 2, NCELLX-1, cuegy, cuflx,
				&cu_boundary_values[1*stride]);
  cuPackBoundary<<<zdim,xdim>>>(1, 0, 2, 0, cuegy, cuflx,
				&cu_boundary_values[2*stride]);
  cuPackBoundary<<<zdim,xdim>>>(1, 0, 2, NCELLY-1, cuegy, cuflx,
				&cu_boundary_values[3*stride]);
  cuPackBoundary<<<ydim,xdim>>>(2, 0, 1, 0, cuegy, cuflx,
				&cu_boundary_values[4*stride]);
  cuPackBoundary<<<ydim,xdim>>>(2, 0, 1, NCELLZ-1, cuegy, cuflx,
				&cu_boundary_values[5*stride]);

  cudaMemcpy(boundary_values,
	     cu_boundary_values,
	     6*4*NBUFF*NBUFF*sizeof(double),
	     cudaMemcpyDeviceToHost);
}

extern "C" void aton_cpu_to_gpu_boundary_(const double *boundary_values) {
  cudaMemcpy(cu_boundary_values,
	     boundary_values,
	     6*4*NBUFF*NBUFF*sizeof(double),
	     cudaMemcpyHostToDevice);

  dim3 xdim(NCELLX);
  dim3 ydim(NCELLY);
  dim3 zdim(NCELLZ);

  int stride = 4*NBUFF*NBUFF;

  cuUnpackBoundary<<<zdim,ydim>>>(0, 1, 2, -1, cuegy, cuflx,
				  &cu_boundary_values[0*stride]);
  cuUnpackBoundary<<<zdim,ydim>>>(0, 1, 2, NCELLX, cuegy, cuflx,
				  &cu_boundary_values[1*stride]);
  cuUnpackBoundary<<<zdim,xdim>>>(1, 0, 2, -1, cuegy, cuflx,
				  &cu_boundary_values[2*stride]);
  cuUnpackBoundary<<<zdim,xdim>>>(1, 0, 2, NCELLY, cuegy, cuflx,
				  &cu_boundary_values[3*stride]);
  cuUnpackBoundary<<<ydim,xdim>>>(2, 0, 1, -1, cuegy, cuflx,
				  &cu_boundary_values[4*stride]);
  cuUnpackBoundary<<<ydim,xdim>>>(2, 0, 1, NCELLZ, cuegy, cuflx,
				  &cu_boundary_values[5*stride]);
}

namespace aton {
  void cpu_to_gpu_boundary_values(const double *values) {
    aton_cpu_to_gpu_boundary_(values);
  }
  
  void gpu_to_cpu_boundary_values(double *values) {
    aton_gpu_to_cpu_boundary_(values);
  }
}

extern "C" void aton_gpu_to_cpu_full_(double *e, double *f, double *xion, double *temp, double *dens, double *src, int *srcpos, int *nsrc) {
  cudaMemcpy(e,cuegy,((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double),cudaMemcpyDeviceToHost);
  cudaMemcpy(f,cuflx,((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double)*3,cudaMemcpyDeviceToHost);
  cudaMemcpy(xion,cuxion,((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double),cudaMemcpyDeviceToHost);
  cudaMemcpy(temp,cutemperature,((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double),cudaMemcpyDeviceToHost);
  cudaMemcpy(dens,cudensity,((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double),cudaMemcpyDeviceToHost);
  cudaMemcpy(src,cusrc0,(*nsrc)*sizeof(double),cudaMemcpyDeviceToHost); 
  cudaMemcpy(srcpos,cusrc0pos,3*(*nsrc)*sizeof(int),cudaMemcpyDeviceToHost);
}

extern "C" void aton_cpu_to_gpu_full_(double *e, double *f, double *xion, double *temp, double *dens, double *src, int *srcpos, int *nsrc, double *photon_source) {
  cudaMemcpy(cuegy,e,((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(cuflx,f,((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double)*3,cudaMemcpyHostToDevice);
  cudaMemcpy(cuxion,xion,((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(cutemperature,temp,((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(cudensity,dens,((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(cusrc0,src,(*nsrc)*sizeof(double),cudaMemcpyHostToDevice);
  cudaMemcpy(cusrc0pos,srcpos,3*(*nsrc)*sizeof(int),cudaMemcpyHostToDevice);

  cudaMemcpy(cuegy_new,cuegy,((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double),cudaMemcpyDeviceToDevice);
  cudaMemcpy(cuflx_new,cuflx,((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double)*3,cudaMemcpyDeviceToDevice);

  cudaMemcpy(cu_photon_source,photon_source,((NCELLX+NBOUND2)*(NCELLY+NBOUND2)*(NCELLZ+NBOUND2))*sizeof(double),cudaMemcpyHostToDevice);

  // Clear the boundary value buffer.
  cudaMemset(cu_boundary_values,
	     0,
	     6*4*NBUFF*NBUFF*sizeof(double));
}

extern "C" void aton_debug_dump_(double *e, double *f, double *x, double *temp, double *src0, int *src0pos, double *dens, int *nsource, double *photon_source, double *time, int *isnap)
{
  char fname[256];
  FILE *fp;

  sprintf(fname,"runs/out.%05d",*isnap);

  int nc=NCELLX+NBOUND2;
  double aexp=1.;

  fp=fopen(fname,"wb");
  fwrite(&nc,sizeof(int),1,fp);
  fwrite(nsource,sizeof(int),1,fp);
  fwrite(time,sizeof(double),1,fp);
  fwrite(e,sizeof(double),NCELLS3,fp);
  fwrite(f,sizeof(double),3*NCELLS3,fp);
  fwrite(x,sizeof(double),NCELLS3,fp);
  fwrite(temp,sizeof(double),NCELLS3,fp);
  fwrite(src0,sizeof(double),*nsource,fp);
  fwrite(src0pos,sizeof(int),3*(*nsource),fp);
  fwrite(dens,sizeof(double),NCELLS3,fp);
  fwrite(&aexp,sizeof(double),1,fp);
  fclose(fp);

  sprintf(fname,"runs/test.%05d",*isnap);
  printf("writing %s\n", fname);

  FILE* ft = fopen(fname, "w");
  fprintf(ft, "# i j k  xion dens energy temperature photon_source\n");
  for (int i=0; i<NCELLX; ++i) {
    for (int j=0; j<NCELLY; ++j) {
      for (int k=0; k<NCELLZ; ++k) {
	const int index =
	  (i+NBOUND) +
	  (j+NBOUND)*(NCELLX+2*NBOUND) +
	  (k+NBOUND)*(NCELLX+2*NBOUND)*(NCELLY+2*NBOUND);
	fprintf(ft, "%d %d %d %e %e %e %e %e\n",
		i, j, k,
		x[index], dens[index], e[index], temp[index],
		photon_source[index]);
      }
    }
  }
  fclose(ft);

  printf("Dump %s done\n", fname);
}
