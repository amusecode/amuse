#include "cuda_fastkick.h"


#define HANDLE_CUDA_ERROR(cuda_error) { \
    cudaError error = cuda_error; \
    if (error != cudaSuccess) { \
        cerr << "Cuda error in file '" << __FILE__ << "' in line " << __LINE__ << ": "; \
        cerr << cudaGetErrorString(error) << endl; \
        return -1; \
    } \
}

#define THREADS_PER_BLOCK 512

__global__ void dev_get_gravity_at_point(
        float eps2, float *eps, float *xh, float *yh, float *zh, float *xt, float *yt, float *zt, 
        float *ax, float *ay, float *az, int n, 
        float *field_m, float *fxh, float *fyh, float *fzh, float *fxt, float *fyt, float *fzt, int n_field) {
    float dx, dy, dz, r2, tmp, dr2, eps2_total;
    for (int tid=threadIdx.x + blockIdx.x*blockDim.x; tid < n; tid += blockDim.x*gridDim.x){
        eps2_total = eps2 + eps[tid]*eps[tid];
        ax[tid] = 0;
        ay[tid] = 0;
        az[tid] = 0;
        for (int i=0; i < n_field; i++){
            dx = (fxh[i] - xh[tid]) + (fxt[i] - xt[tid]);
            dy = (fyh[i] - yh[tid]) + (fyt[i] - yt[tid]);
            dz = (fzh[i] - zh[tid]) + (fzt[i] - zt[tid]);
            dr2 = dx*dx + dy*dy + dz*dz;
            if (dr2 > 0) {
                r2 = eps2_total + dr2;
                tmp = field_m[i] / (r2 * sqrt(r2));
                ax[tid] += tmp * dx;
                ay[tid] += tmp * dy;
                az[tid] += tmp * dz;
            }
        }
    }
}

__global__ void dev_get_potential_at_point(
        float eps2, float *eps, float *xh, float *yh, float *zh, float *xt, float *yt, float *zt, 
        float *phi, int n, 
        float *field_m, float *fxh, float *fyh, float *fzh, float *fxt, float *fyt, float *fzt, int n_field) {
    float dx, dy, dz, r, dr2, eps2_total;
    for (int tid=threadIdx.x + blockIdx.x*blockDim.x; tid < n; tid += blockDim.x*gridDim.x){
        eps2_total = eps2 + eps[tid]*eps[tid];
        phi[tid] = 0;
        for (int i=0; i < n_field; i++){
            dx = (fxh[i] - xh[tid]) + (fxt[i] - xt[tid]);
            dy = (fyh[i] - yh[tid]) + (fyt[i] - yt[tid]);
            dz = (fzh[i] - zh[tid]) + (fzt[i] - zt[tid]);
            dr2 = dx*dx + dy*dy + dz*dz;
            if (dr2 > 0) {
                r = sqrt(eps2_total + dr2);
                phi[tid] -= field_m[i] / r;
            }
        }
    }
}

__global__ void dev_get_potential_energy(
        float *partial_results, float eps2, float *field_m, 
        float *fxh, float *fyh, float *fzh, float *fxt, float *fyt, float *fzt, int n_field) {
    extern __shared__ float thread_results[];
    unsigned int i, j;
    float dx, dy, dz, r, dr2, potential_energy = 0;
    for (j=threadIdx.x + blockIdx.x*blockDim.x; j < n_field; j += blockDim.x*gridDim.x){
        for (i=0; i<j; i++){
            dx = (fxh[i] - fxh[j]) + (fxt[i] - fxt[j]);
            dy = (fyh[i] - fyh[j]) + (fyt[i] - fyt[j]);
            dz = (fzh[i] - fzh[j]) + (fzt[i] - fzt[j]);
            dr2 = dx*dx + dy*dy + dz*dz;
            r = sqrt(eps2 + dr2);
            potential_energy -= field_m[i]*field_m[j] / r;
        }
    }
    
    // Reduce results from all threads within this block
    thread_results[threadIdx.x] = potential_energy;
    __syncthreads();
    for (i = blockDim.x/2; i>0; i>>=1) {
        if (threadIdx.x < i) {
            thread_results[threadIdx.x] += thread_results[threadIdx.x + i];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        partial_results[blockIdx.x] = thread_results[0];
    }
}

// field particle data
int number_of_field_particles;
float *dev_m;
float *dev_x_head, *dev_y_head, *dev_z_head;
float *dev_x_tail, *dev_y_tail, *dev_z_tail;
int particle_data_allocated_on_device;

// data for points in get_..._at_point()
float *dev_eps;
float *dev_px_head, *dev_py_head, *dev_pz_head;
float *dev_px_tail, *dev_py_tail, *dev_pz_tail;
float *dev_ax, *dev_ay, *dev_az;
float *dev_phi;


int cuda_initialize_code() {
    particle_data_allocated_on_device = 0;
    return 0;
}

int cuda_commit_particles(vector<double>& m, vector<double>& x, vector<double>& y, vector<double>& z) {
    cerr << "cuda_commit_particles" << endl;
    int N = x.size();
    float *x_head = new float[N];
    float *y_head = new float[N];
    float *z_head = new float[N];
    float *x_tail = new float[N];
    float *y_tail = new float[N];
    float *z_tail = new float[N];
    float *m_float = new float[N];
    for (int i=0; i<N; i++) {
        x_head[i] = (float) x[i];
        y_head[i] = (float) y[i];
        z_head[i] = (float) z[i];
        x_tail[i] = (float) (x[i] - x_head[i]);
        y_tail[i] = (float) (y[i] - y_head[i]);
        z_tail[i] = (float) (z[i] - z_head[i]);
        m_float[i] = (float) m[i];
    }
    
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_x_head, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_y_head, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_z_head, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_x_tail, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_y_tail, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_z_tail, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_m, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_x_head, x_head, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_y_head, y_head, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_z_head, z_head, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_x_tail, x_tail, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_y_tail, y_tail, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_z_tail, z_tail, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_m, m_float, N*sizeof(float), cudaMemcpyHostToDevice));
    delete[] x_head;
    delete[] y_head;
    delete[] z_head;
    delete[] x_tail;
    delete[] y_tail;
    delete[] z_tail;
    delete[] m_float;
    number_of_field_particles = N;
    particle_data_allocated_on_device = 1;
    cerr << "cuda_commit_particles finished successfully" << endl;
    return 0;
}


int cuda_cleanup_code() {
    if (particle_data_allocated_on_device == 1) {
        HANDLE_CUDA_ERROR(cudaFree(dev_x_head));
        HANDLE_CUDA_ERROR(cudaFree(dev_y_head));
        HANDLE_CUDA_ERROR(cudaFree(dev_z_head));
        HANDLE_CUDA_ERROR(cudaFree(dev_x_tail));
        HANDLE_CUDA_ERROR(cudaFree(dev_y_tail));
        HANDLE_CUDA_ERROR(cudaFree(dev_z_tail));
        HANDLE_CUDA_ERROR(cudaFree(dev_m));
        particle_data_allocated_on_device = 0;
    }
    return 0;
}

int cuda_get_gravity_at_point(double eps2, double *eps, double *x, double *y, double *z, double *ax, double *ay, double *az, int N) {
    cerr << "cuda_get_gravity_at_point" << endl;
    float *x_head = new float[N];
    float *y_head = new float[N];
    float *z_head = new float[N];
    float *x_tail = new float[N];
    float *y_tail = new float[N];
    float *z_tail = new float[N];
    float *eps_float = new float[N];
    for (int i=0; i<N; i++) {
        x_head[i] = (float) x[i];
        y_head[i] = (float) y[i];
        z_head[i] = (float) z[i];
        x_tail[i] = (float) (x[i] - x_head[i]);
        y_tail[i] = (float) (y[i] - y_head[i]);
        z_tail[i] = (float) (z[i] - z_head[i]);
        eps_float[i] = (float) eps[i];
    }
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_eps, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_px_head, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_py_head, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_pz_head, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_px_tail, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_py_tail, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_pz_tail, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_ax, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_ay, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_az, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_eps, eps_float, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_px_head, x_head, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_py_head, y_head, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_pz_head, z_head, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_px_tail, x_tail, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_py_tail, y_tail, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_pz_tail, z_tail, N*sizeof(float), cudaMemcpyHostToDevice));
    
    int blocks_per_grid = (N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    if (blocks_per_grid > 32) blocks_per_grid = 32;
    dev_get_gravity_at_point<<<blocks_per_grid, THREADS_PER_BLOCK>>>(eps2, dev_eps,
        dev_px_head, dev_py_head, dev_pz_head, dev_px_tail, dev_py_tail, dev_pz_tail,
        dev_ax, dev_ay, dev_az, N,
        dev_m, dev_x_head, dev_y_head, dev_z_head, dev_x_tail, dev_y_tail, dev_z_tail, number_of_field_particles);
    
    HANDLE_CUDA_ERROR(cudaMemcpy(x_head, dev_ax, N*sizeof(float), cudaMemcpyDeviceToHost));
    HANDLE_CUDA_ERROR(cudaMemcpy(y_head, dev_ay, N*sizeof(float), cudaMemcpyDeviceToHost));
    HANDLE_CUDA_ERROR(cudaMemcpy(z_head, dev_az, N*sizeof(float), cudaMemcpyDeviceToHost));
    for (int i=0; i<N; i++) {
        ax[i] = x_head[i];
        ay[i] = y_head[i];
        az[i] = z_head[i];
    }
    HANDLE_CUDA_ERROR(cudaFree(dev_eps));
    HANDLE_CUDA_ERROR(cudaFree(dev_px_head));
    HANDLE_CUDA_ERROR(cudaFree(dev_py_head));
    HANDLE_CUDA_ERROR(cudaFree(dev_pz_head));
    HANDLE_CUDA_ERROR(cudaFree(dev_px_tail));
    HANDLE_CUDA_ERROR(cudaFree(dev_py_tail));
    HANDLE_CUDA_ERROR(cudaFree(dev_pz_tail));
    HANDLE_CUDA_ERROR(cudaFree(dev_ax));
    HANDLE_CUDA_ERROR(cudaFree(dev_ay));
    HANDLE_CUDA_ERROR(cudaFree(dev_az));
    delete[] x_head;
    delete[] y_head;
    delete[] z_head;
    delete[] x_tail;
    delete[] y_tail;
    delete[] z_tail;
    delete[] eps_float;
    cerr << "cuda_get_gravity_at_point finished successfully" << endl;
    return 0;
}

int cuda_get_potential_at_point(double eps2, double *eps, double *x, double *y, double *z, double *phi, int N) {
    cerr << "cuda_get_potential_at_point" << endl;
    float *x_head = new float[N];
    float *y_head = new float[N];
    float *z_head = new float[N];
    float *x_tail = new float[N];
    float *y_tail = new float[N];
    float *z_tail = new float[N];
    float *eps_float = new float[N];
    for (int i=0; i<N; i++) {
        x_head[i] = (float) x[i];
        y_head[i] = (float) y[i];
        z_head[i] = (float) z[i];
        x_tail[i] = (float) (x[i] - x_head[i]);
        y_tail[i] = (float) (y[i] - y_head[i]);
        z_tail[i] = (float) (z[i] - z_head[i]);
        eps_float[i] = (float) eps[i];
    }
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_eps, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_px_head, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_py_head, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_pz_head, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_px_tail, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_py_tail, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_pz_tail, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_phi, N*sizeof(float)));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_eps, eps_float, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_px_head, x_head, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_py_head, y_head, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_pz_head, z_head, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_px_tail, x_tail, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_py_tail, y_tail, N*sizeof(float), cudaMemcpyHostToDevice));
    HANDLE_CUDA_ERROR(cudaMemcpy(dev_pz_tail, z_tail, N*sizeof(float), cudaMemcpyHostToDevice));
    
    int blocks_per_grid = (N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    if (blocks_per_grid > 32) blocks_per_grid = 32;
    dev_get_potential_at_point<<<blocks_per_grid, THREADS_PER_BLOCK>>>(eps2, dev_eps,
        dev_px_head, dev_py_head, dev_pz_head, dev_px_tail, dev_py_tail, dev_pz_tail,
        dev_phi, N,
        dev_m, dev_x_head, dev_y_head, dev_z_head, dev_x_tail, dev_y_tail, dev_z_tail, number_of_field_particles);
    
    HANDLE_CUDA_ERROR(cudaMemcpy(x_head, dev_phi, N*sizeof(float), cudaMemcpyDeviceToHost));
    for (int i=0; i<N; i++) {
        phi[i] = x_head[i];
    }
    HANDLE_CUDA_ERROR(cudaFree(dev_eps));
    HANDLE_CUDA_ERROR(cudaFree(dev_px_head));
    HANDLE_CUDA_ERROR(cudaFree(dev_py_head));
    HANDLE_CUDA_ERROR(cudaFree(dev_pz_head));
    HANDLE_CUDA_ERROR(cudaFree(dev_px_tail));
    HANDLE_CUDA_ERROR(cudaFree(dev_py_tail));
    HANDLE_CUDA_ERROR(cudaFree(dev_pz_tail));
    HANDLE_CUDA_ERROR(cudaFree(dev_phi));
    delete[] x_head;
    delete[] y_head;
    delete[] z_head;
    delete[] x_tail;
    delete[] y_tail;
    delete[] z_tail;
    delete[] eps_float;
    cerr << "cuda_get_potential_at_point finished successfully" << endl;
    return 0;
}
int cuda_get_potential_energy(double eps2, double *potential_energy) {
    cerr << "cuda_get_potential_energy" << endl;
    int blocks_per_grid = (number_of_field_particles + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    if (blocks_per_grid > 32) blocks_per_grid = 32;
    float *partial_results = new float[blocks_per_grid];
    float *dev_partial_results;
    HANDLE_CUDA_ERROR(cudaMalloc((void **) &dev_partial_results, blocks_per_grid*sizeof(float)));
    int shared_memory = 2 * THREADS_PER_BLOCK * sizeof(float);
    
    dev_get_potential_energy<<<blocks_per_grid, THREADS_PER_BLOCK, shared_memory>>>(
        dev_partial_results, eps2, 
        dev_m, dev_x_head, dev_y_head, dev_z_head, dev_x_tail, dev_y_tail, dev_z_tail, 
        number_of_field_particles);
    
    HANDLE_CUDA_ERROR(cudaMemcpy(partial_results, dev_partial_results, blocks_per_grid*sizeof(float), cudaMemcpyDeviceToHost));
    HANDLE_CUDA_ERROR(cudaFree(dev_partial_results));
    
    *potential_energy = 0;
    for (int i=0; i < blocks_per_grid; i++) {
        *potential_energy += partial_results[i];
    }
    delete[] partial_results;
    cerr << "cuda_get_potential_energy finished successfully" << endl;
    return 0;
}
