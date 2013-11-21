// Written by: Pablo Bena Llambay / Mark Gieles
// Calculates specific potential of a cluster (N1 = 0) or a cluster pair (0 < N1 < N)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cuda.h>
#include "mameclot.h"

#define BLOCKSIZE 256
	
//Computing the Potential on the Device
__global__ void compute_potential_gpu(float *m, 
	   float *x, float *y, float *z, float *phi, int N, int N1) {
  int i,j; 
  float rijx, rijy, rijz;
  float xi, yi, zi;
  float potential;	
  i = threadIdx.x + blockIdx.x*blockDim.x;
  if (i < (N1 == 0 ? N : N1)) 
  {
    xi = x[i];
    yi = y[i];
    zi = z[i];

    for (j = (N1 == 0 ? 0 : N1); j < N; j++) 
    {
      rijx = xi - x[j];
      rijy = yi - y[j];
      rijz = zi - z[j];
     
      if (i!=j)  
         potential -= m[j]/sqrt(rijx*rijx + rijy*rijy + rijz*rijz);
    }
   phi[i] = potential;
   }   
}


extern "C" void calculate_potential(float *m, float *x, float *y, float *z,
        float *phi, int N, int N1) 
{
  float *m_d,*x_d,*y_d,*z_d,*phi_d; // Device variables!

  //Allocating memory on the Device
  cudaMalloc(&m_d  , sizeof(float)*N); 
  cudaMalloc(&x_d  , sizeof(float)*N); 
  cudaMalloc(&y_d  , sizeof(float)*N);
  cudaMalloc(&z_d  , sizeof(float)*N); 
  cudaMalloc(&phi_d, sizeof(float)*N);

  cudaMemcpy(m_d,m    , sizeof(float)*N, cudaMemcpyHostToDevice); // Host -> Device
  cudaMemcpy(x_d,x    , sizeof(float)*N, cudaMemcpyHostToDevice); // Host -> Device
  cudaMemcpy(y_d,y    , sizeof(float)*N, cudaMemcpyHostToDevice); // Host -> Device
  cudaMemcpy(z_d,z    , sizeof(float)*N, cudaMemcpyHostToDevice); // Host -> Device
  cudaMemcpy(phi_d,phi, sizeof(float)*N, cudaMemcpyHostToDevice); // Host -> Device

  compute_potential_gpu <<<((N+BLOCKSIZE-1))/BLOCKSIZE,BLOCKSIZE >>>(m_d,x_d, y_d, z_d, phi_d,N,N1);
  cudaMemcpy(phi,phi_d, sizeof(float)*N, cudaMemcpyDeviceToHost); // Host -> Device
    
  //Freeing memory
  cudaFree(m_d);
  cudaFree(x_d);
  cudaFree(y_d);
  cudaFree(z_d);
  cudaFree(phi_d);
  
  return;
}

