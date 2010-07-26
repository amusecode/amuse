#include <cutil.h>
#include <multithreading.h>
#include <stdio.h>
#include "sapporo_defs.h"

#include "dev_evaluate_gravity.cu"


double get_time();

inline int n_norm(int n, int j) {
  n = ((n-1)/j) * j + j;
  if (n == 0) n = j;
  return n;
}

#include "sapporo_multi.h"

extern "C"
{

  int get_device_count() {
    int gpuCount = 0;
    CUDA_SAFE_CALL(cudaGetDeviceCount(&gpuCount));
    CUT_CHECK_ERROR("Failed to get CUDA device count\n");
    return gpuCount;
  }
  


//#define _DEBUG_PRINT_

  cudaError_t host_evaluate_gravity(sapporo_multi_struct gpu) {
    
    double t0 = get_time();
    
    int ofs = gpu.offset;
    
    
    DS t_i = (DS){gpu.t_i_x, gpu.t_i_y};
    if (gpu.nj_modified > 0) {
    
      int nj_scaled = n_norm(gpu.nj_modified, NTHREADS);
      dim3 threads(NTHREADS, 1, 1);
      dim3 grid(nj_scaled/NTHREADS, 1, 1);
      if (nj_scaled < threads.x) {
        threads.x = nj_scaled;
        grid.x    = 1;
      };
#ifdef _DEBUG_PRINT_
      double t1 = get_time();
#endif
      dev_copy_particles<<<grid, threads>>>(gpu.nj_modified,
                                            gpu.nj_max,
                                            gpu.address_j,
                                            gpu.t_j,
                                            gpu.Ppos_j,
                                            gpu.Pvel_j,
                                            gpu.pos_j,
                                            gpu.vel_j,
                                            gpu.acc_j,
                                            gpu.jrk_j);
#ifdef _DEBUG_PRINT_
      fprintf(stderr, " dev_copy_particles:    %lf sec\n", get_time() - t1);
#endif
    }

    if (gpu.predict) {
      
      int nj_scaled = n_norm(gpu.nj, NTHREADS);
      dim3 threads_p(NTHREADS, 1, 1);
      dim3 grid_p(nj_scaled/NTHREADS, 1, 1);
      
      double t1 = get_time();
      dev_predictor<<<grid_p, threads_p>>>(gpu.nj, 
                                           t_i, 
                                           gpu.t_j + ofs, 
                                           gpu.Ppos_j + ofs, 
                                           gpu.Pvel_j + ofs,
                                           gpu.pos_j  + ofs,
                                           gpu.vel_j  + ofs, 
                                           gpu.acc_j  + ofs, 
                                           gpu.jrk_j  + ofs);

#ifdef _DEBUG_PRINT_
      fprintf(stderr, "  dev_predict:          %lf sec\n", get_time() - t1);
#endif

    };
      
    int p = gpu.ni;
    int q = min(NTHREADS/gpu.ni, 32);
//     q = 1;
    dim3 threads(p, q, 1);
    dim3 grid(NBLOCKS, 1, 1);

    int shared_mem_size = p*q*(sizeof(DS4) + sizeof(float4));
    int nj_scaled = n_norm(gpu.nj, q*NBLOCKS);
    
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(
        "EPS2", 
        &gpu.EPS2, 
        sizeof(float), 
        0, 
        cudaMemcpyHostToDevice)
    );
    
    double t1 = get_time();
    if (gpu.ngb)
      dev_evaluate_gravity<true><<<grid, threads, shared_mem_size>>>(gpu.nj, 
                                                                     nj_scaled/(NBLOCKS*q),
                                                                     NTHREADS,
                                                                     gpu.Ppos_j+ ofs,
                                                                     gpu.Pvel_j + ofs,
                                                                     gpu.pos_i,  gpu.vel_i,
                                                                     gpu.acc_i,  gpu.jrk_i,
                                                                     gpu.ngb_list);
    else
      dev_evaluate_gravity<false><<<grid, threads, shared_mem_size>>>(gpu.nj, 
                                                                      nj_scaled/(NBLOCKS*q),
                                                                      NTHREADS,
                                                                      gpu.Ppos_j + ofs, 
                                                                      gpu.Pvel_j + ofs,
                                                                      gpu.pos_i,  gpu.vel_i,
                                                                      gpu.acc_i,  gpu.jrk_i,
                                                                      gpu.ngb_list);

#ifdef _DEBUG_PRINT_
     fprintf(stderr, "  dev_evaluate_gravity: %lf sec\n", get_time() - t1);
#endif
     
    dim3 threads_r(NBLOCKS, 1, 1);
    dim3 grid_r(gpu.ni, 1, 1);
    int shared_mem_size_r= NBLOCKS*(2*sizeof(float4) + 3*sizeof(int)); 
    
    t1 = get_time();
    dev_reduce_forces<<<grid_r, threads_r, shared_mem_size_r>>>(gpu.acc_i, 
                                                                gpu.jrk_i, 
                                                                gpu.ds_i,
                                                                gpu.vel_i,
                                                                NTHREADS,
                                                                NGB_PB*NBLOCKS*NTHREADS,
                                                                gpu.ngb_list);

#ifdef _DEBUG_PRINT_
    fprintf(stderr, "  dev_reduce_forces:    %lf sec\n", get_time() - t1);
#endif
    
    return cudaSuccess;

 }


}
