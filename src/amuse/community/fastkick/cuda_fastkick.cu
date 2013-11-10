#include "cuda_fastkick.h"


#define HANDLE_CUDA_ERROR(cuda_error) { \
    cudaError error = cuda_error; \
    if (error != cudaSuccess) { \
        cerr << "Cuda error in file '" << __FILE__ << "' in line " << __LINE__ << ": "; \
        cerr << cudaGetErrorString(error) << endl; \
        return -1; \
    } \
}

#define THREADS_PER_BLOCK 128

__global__ void dev_get_gravity_at_point(
        float eps2, float *eps, float *xh, float *yh, float *zh, float *xt, float *yt, float *zt, 
        float *ax, float *ay, float *az, int n, 
        float *field_m, float *fxh, float *fyh, float *fzh, float *fxt, float *fyt, float *fzt, int n_field) {
    float dx, dy, dz, r2, tmp;
    for (int tid=threadIdx.x + blockIdx.x*blockDim.x; tid < n; tid += blockDim.x*gridDim.x){
        ax[tid] = 0;
        ay[tid] = 0;
        az[tid] = 0;
        for (int i=0; i < n_field; i++){
            dx = (fxh[i] - xh[tid]) + (fxt[i] - xt[tid]);
            dy = (fyh[i] - yh[tid]) + (fyt[i] - yt[tid]);
            dz = (fzh[i] - zh[tid]) + (fzt[i] - zt[tid]);
            r2 = eps2 + eps[tid]*eps[tid] + dx*dx + dy*dy + dz*dz;
            tmp = field_m[i] / (r2 * sqrt(r2));
            ax[tid] += tmp * dx;
            ay[tid] += tmp * dy;
            az[tid] += tmp * dz;
        }
    }
}

__global__ void dev_get_potential_at_point(
        float eps2, float *eps, float *xh, float *yh, float *zh, float *xt, float *yt, float *zt, 
        float *phi, int n, 
        float *field_m, float *fxh, float *fyh, float *fzh, float *fxt, float *fyt, float *fzt, int n_field) {
    float dx, dy, dz, r;
    for (int tid=threadIdx.x + blockIdx.x*blockDim.x; tid < n; tid += blockDim.x*gridDim.x){
        phi[tid] = 0;
        for (int i=0; i < n_field; i++){
            dx = (fxh[i] - xh[tid]) + (fxt[i] - xt[tid]);
            dy = (fyh[i] - yh[tid]) + (fyt[i] - yt[tid]);
            dz = (fzh[i] - zh[tid]) + (fzt[i] - zt[tid]);
            r = sqrt(eps2 + eps[tid]*eps[tid] + dx*dx + dy*dy + dz*dz);
            phi[tid] -= field_m[i] / r;
        }
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
    int N = x.size();
    float x_head[N], y_head[N], z_head[N];
    float x_tail[N], y_tail[N], z_tail[N];
    float m_float[N];
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
    number_of_field_particles = N;
    particle_data_allocated_on_device = 1;
    return 0;
}

int cuda_recommit_particles(vector<double>& m, vector<double>& x, vector<double>& y, vector<double>& z) {
    HANDLE_CUDA_ERROR(cudaFree(dev_x_head));
    HANDLE_CUDA_ERROR(cudaFree(dev_y_head));
    HANDLE_CUDA_ERROR(cudaFree(dev_z_head));
    HANDLE_CUDA_ERROR(cudaFree(dev_x_tail));
    HANDLE_CUDA_ERROR(cudaFree(dev_y_tail));
    HANDLE_CUDA_ERROR(cudaFree(dev_z_tail));
    HANDLE_CUDA_ERROR(cudaFree(dev_m));
    return cuda_commit_particles(m, x, y, z);
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
    float x_head[N], y_head[N], z_head[N];
    float x_tail[N], y_tail[N], z_tail[N];
    float eps_float[N], resultx[N], resulty[N], resultz[N];
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
    
    HANDLE_CUDA_ERROR(cudaMemcpy(resultx, dev_ax, N*sizeof(float), cudaMemcpyDeviceToHost));
    HANDLE_CUDA_ERROR(cudaMemcpy(resulty, dev_ay, N*sizeof(float), cudaMemcpyDeviceToHost));
    HANDLE_CUDA_ERROR(cudaMemcpy(resultz, dev_az, N*sizeof(float), cudaMemcpyDeviceToHost));
    for (int i=0; i<N; i++) {
        ax[i] = resultx[i];
        ay[i] = resulty[i];
        az[i] = resultz[i];
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
    return 0;
}

int cuda_get_potential_at_point(double eps2, double *eps, double *x, double *y, double *z, double *phi, int N) {
    float x_head[N], y_head[N], z_head[N];
    float x_tail[N], y_tail[N], z_tail[N];
    float eps_float[N], result[N];
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
    
    HANDLE_CUDA_ERROR(cudaMemcpy(result, dev_phi, N*sizeof(float), cudaMemcpyDeviceToHost));
    for (int i=0; i<N; i++) {
        phi[i] = result[i];
    }
    HANDLE_CUDA_ERROR(cudaFree(dev_eps));
    HANDLE_CUDA_ERROR(cudaFree(dev_px_head));
    HANDLE_CUDA_ERROR(cudaFree(dev_py_head));
    HANDLE_CUDA_ERROR(cudaFree(dev_pz_head));
    HANDLE_CUDA_ERROR(cudaFree(dev_px_tail));
    HANDLE_CUDA_ERROR(cudaFree(dev_py_tail));
    HANDLE_CUDA_ERROR(cudaFree(dev_pz_tail));
    HANDLE_CUDA_ERROR(cudaFree(dev_phi));
    return 0;
}
