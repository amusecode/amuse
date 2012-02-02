#include "sapporo.h"

void sapporo::free_cuda_memory(int ignore) {
    dev_struct &dev = device;
    
    CUDA_SAFE_CALL(cudaFree( (void*)dev.Ppos_j));
    CUDA_SAFE_CALL(cudaFree( (void*)dev.Pvel_j));
    
    CUDA_SAFE_CALL(cudaFree( (void*)dev.address_j));
    
    CUDA_SAFE_CALL(cudaFree( (void*)dev.t_j));
    CUDA_SAFE_CALL(cudaFree( (void*)dev.pos_j));
    CUDA_SAFE_CALL(cudaFree( (void*)dev.vel_j));
    CUDA_SAFE_CALL(cudaFree( (void*)dev.acc_j));
    CUDA_SAFE_CALL(cudaFree( (void*)dev.jrk_j));
    
    CUDA_SAFE_CALL(cudaFree( (void*)dev.pos_i));
    CUDA_SAFE_CALL(cudaFree( (void*)dev.vel_i));
    CUDA_SAFE_CALL(cudaFree( (void*)dev.acc_i));
    CUDA_SAFE_CALL(cudaFree( (void*)dev.jrk_i));
    CUDA_SAFE_CALL(cudaFree( (void*)dev.ds_i));
    
    CUDA_SAFE_CALL(cudaFree( (void*)dev.ngb_list_i));
    
    
    CUDA_SAFE_CALL(cudaThreadExit());
    CUT_CHECK_ERROR("Failedn");
}

void sapporo::allocate_cuda_memory(int ignore) {
    dev_struct &dev = device;
    
    int memory_needed = nj_max * sizeof(DS4) * (1 + 2);
    memory_needed += nj_max * sizeof(float4) * (1 + 2 + 2 + 2);
    memory_needed += nj_max * sizeof(int) * (1);
    memory_needed += nj_max * sizeof(DS2) * (2);
    memory_needed += n_pipes * sizeof(DS4);
    memory_needed += n_pipes * (1 + NBLOCKS + NBLOCKS + NBLOCKS) * sizeof(float4);
    memory_needed += n_pipes * sizeof(float);
    memory_needed += (n_pipes*(NGB_PP + 1) + n_pipes*NBLOCKS*(NGB_PP+1)) * sizeof(int);
    cerr << "bytes of memory needed :" << memory_needed << endl;
    
    
    CUDA_SAFE_CALL(cudaMalloc( (void**)((void*)&dev.Ppos_j), nj_max * sizeof(DS4)));
    CUDA_SAFE_CALL(cudaMalloc( (void**)((void*)&dev.Pvel_j), nj_max * sizeof(float4)));
    
    CUDA_SAFE_CALL(cudaMalloc( (void**)((void*)&dev.address_j), nj_max * sizeof(int)));
    
    CUDA_SAFE_CALL(cudaMalloc( (void**)((void*)&dev.t_j),   2*nj_max * sizeof(DS2)));
    CUDA_SAFE_CALL(cudaMalloc( (void**)((void*)&dev.pos_j), 2*nj_max * sizeof(DS4)));
    CUDA_SAFE_CALL(cudaMalloc( (void**)((void*)&dev.vel_j), 2*nj_max * sizeof(float4)));
    CUDA_SAFE_CALL(cudaMalloc( (void**)((void*)&dev.acc_j), 2*nj_max * sizeof(float4)));
    CUDA_SAFE_CALL(cudaMalloc( (void**)((void*)&dev.jrk_j), 2*nj_max * sizeof(float4)));
    
    CUDA_SAFE_CALL(cudaMalloc( (void**)((void*)&dev.pos_i), n_pipes                 * sizeof(DS4)));
    CUDA_SAFE_CALL(cudaMalloc( (void**)((void*)&dev.vel_i), n_pipes * (1 + NBLOCKS) * sizeof(float4)));
    CUDA_SAFE_CALL(cudaMalloc( (void**)((void*)&dev.acc_i), n_pipes *      NBLOCKS  * sizeof(float4)));
    CUDA_SAFE_CALL(cudaMalloc( (void**)((void*)&dev.jrk_i), n_pipes *      NBLOCKS  * sizeof(float4)));
    CUDA_SAFE_CALL(cudaMalloc( (void**)((void*)&dev.ds_i),  n_pipes * sizeof(float)));
    
    CUDA_SAFE_CALL(cudaMalloc( (void**)((void*)&dev.ngb_list_i), (n_pipes*(NGB_PP + 1) + n_pipes*NBLOCKS*(NGB_PP+1)) * sizeof(int)));
}

void sapporo::send_j_particles_to_device(int ignore) {
    dev_struct &dev = device;
    
    int nj = address_j.size();
    
    CUDA_SAFE_CALL(cudaMemcpy( dev.address_j,      &address_j[0], nj * sizeof(int),    cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy( &dev.t_j  [nj_max], &t_j[0],       nj * sizeof(DS2),    cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy( &dev.pos_j[nj_max], &pos_j[0],     nj * sizeof(DS4),    cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy( &dev.vel_j[nj_max], &vel_j[0],     nj * sizeof(float4), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy( &dev.acc_j[nj_max], &acc_j[0],     nj * sizeof(float4), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy( &dev.jrk_j[nj_max], &jrk_j[0],     nj * sizeof(float4), cudaMemcpyHostToDevice));
    
}

void sapporo::send_i_particles_to_device(int ignore, int ni) {
    dev_struct &dev = device;
    
    CUDA_SAFE_CALL(cudaMemcpy( dev.pos_i, &pos_i[0], ni * sizeof(DS4),    cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy( dev.vel_i, &vel_i[0], ni * sizeof(float4), cudaMemcpyHostToDevice));
}

void sapporo::fetch_data_from_device(int ignore, int ni) {
    dev_struct &dev = device;
    
    CUDA_SAFE_CALL(cudaMemcpy( &acc_i[0], dev.acc_i, ni * sizeof(float4), cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy( &jrk_i[0], dev.jrk_i, ni * sizeof(float4), cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy( &ds_i[0],  dev.ds_i,  ni * sizeof(float),  cudaMemcpyDeviceToHost));
}

int sapporo::fetch_ngb_list_from_device(int ignore) {
    dev_struct &dev = device;
    
    int ni = dev.ni;
    
    CUDA_SAFE_CALL(cudaMemcpy( 
               &ngb_list_i[n_pipes*NGB_PP*0], 
               &dev.ngb_list_i[NTHREADS*NGB_PB*NBLOCKS], 
               ni*NGB_PP*sizeof(int), 
               cudaMemcpyDeviceToHost
               ));
    return ni;
}


double sapporo::evaluate_gravity(int ni, int nj) {
#ifdef NGB
    bool ngb = true;
#else
    bool ngb = false;
#endif
    
    
    sapporo_multi_struct sapporo_multi_data[1];
    
    
    dev_struct &dev = device;
    sapporo_multi_struct &gpu = sapporo_multi_data[0];
    
    int nj_i   = nj;
    int offset = 0;
    
    dev.ni  = ni;
    dev.nj  = nj_i;
    
    gpu.device = device_id;
    
    gpu.EPS2 = EPS2;
    
    gpu.ngb = ngb;
    gpu.nj  = dev.nj;
    gpu.ni  = dev.ni;
    gpu.nj_total = nj;
    
    gpu.nj_max = nj_max;
    gpu.nj_modified = nj_modified;
    gpu.predict     = predict;
    
    gpu.t_i_x = t_i.x;
    gpu.t_i_y = t_i.y;
    
    gpu.address_j = dev.address_j;
    gpu.offset = offset;
    
    gpu.t_j = dev.t_j;
    gpu.Ppos_j = dev.Ppos_j;
    gpu.Pvel_j = dev.Pvel_j;
    
    gpu.pos_j = dev.pos_j;
    gpu.vel_j = dev.vel_j;
    gpu.acc_j = dev.acc_j;
    gpu.jrk_j = dev.jrk_j;
    
    gpu.pos_i = dev.pos_i;
    gpu.vel_i = dev.vel_i;
    gpu.acc_i = dev.acc_i;
    gpu.jrk_i = dev.jrk_i;
    gpu.ds_i  = dev.ds_i;
    
    gpu.ngb_list = dev.ngb_list_i;
    host_evaluate_gravity(sapporo_multi_data[0]);
    nj_modified = 0;
    predict     = false;
    return 0;
}
