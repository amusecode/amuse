#ifndef KERNEL_H
#define KERNEL_H
__global__ void Reconstruct ( int *nex,
                              unsigned long nextsize,
                              double4 *A,
                              double4 *B,
										unsigned int offset);

__global__ void initvectors(double4 *acc3, float4 *apred);

__global__ void Predictor (const double TIME,
                           double4 *p_pred,
                           float4  *v_pred,
                           float4  *a_pred,
                           double4 *p_corr,
                           double4 *v_corr,
                           double  *loc_time,
                           double4 *acc,
									double4 *acc1,
									double4 *acc2,
                           double4 *acc3,
									int istart,
									int* nvec,
									int ppgpus,
									unsigned int N);

__global__ void reduce(double4 *ac, unsigned int bf_real, unsigned int dimension);
__global__ void reposition (double4 *ac, double4 *af, unsigned int offset, unsigned long nextsize);

__device__ void AddPlummerEnergy(double4 dr, double *pot, double *b, double *M, double mul);

__global__ void energy(double4 *posCD, double4 *velCD, double *E, unsigned int N, unsigned int istart, unsigned int ppG, double plummer_core, double plummer_mass);

__global__ void potential_energy(double4 *posD, double4 *velCD, double *E, unsigned int N, unsigned int istart, unsigned int ppG, double plummer_core, double plummer_mass);

__device__ void AddPlummer(double4 dr, float4 dv, float4 da, double4 *acc, double4 *jrk, double4 *snp, double *b, double *M, double mul);

__global__  void evaluation(unsigned int N,
									 const double4 *const globalX,
                            const float4 *const globalV,
                            const float4 *const globalA,
									 double4 *aC,
									 double4 *aC1,
									 double4 *aC2,
                            unsigned int istart,
                            int ppGpus,
                            int Bfactor,
                            int dimension,
                            int *next,
                            double *loc_time,
                            double TIME,
									 double pl_core,
									 double pl_mass);
#endif
