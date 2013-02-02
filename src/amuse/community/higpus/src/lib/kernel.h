#ifndef KERNEL_H
#define KERNEL_H

__global__ void sum_partial(double4 *a, double4 *b, unsigned int nextsize);
__global__ void update_local_time(int *next, double *local_time, double GTIME);

__global__ void Corrector_gpu(double GTIME,
                              double *local_time,
                              double *step,
                              int *next,
                              unsigned long nextsize,
                              double4 *pos_CH,
                              double4 *vel_CH,
                              double4 *a_tot_D,
                              double4 *a1_tot_D,
                              double4 *a2_tot_D,
                              double4 *a_H0,
                              double4 *a3_H,
                              double ETA6,
                              double ETA4,
                              double DTMAX,
                              double DTMIN,
                              unsigned int N);

__global__ void Reconstruct(int *nex,
                            unsigned long nextsize,
                            double4 *pc,
                            double4 *vc,
                            double4 *a3,
                            double4 *a,
                            double4 *a1,
                            double4 *a2,
                            double4 *pva3,
                            double4 *aaa);

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

__global__ void reduce(double4 *ac, double4 *ac1, double4 *ac2, unsigned int bf_real, unsigned int dimension);
__global__ void reposition (double4 *ac, double4 *ac1, double4 *ac2, double4 *af, unsigned long nextsize);

__device__ void AddPlummerEnergy(double4 dr, double *pot, double *b, double *M, double mul);


__global__ void energy(double4 *posD, double4 *velD, double *K, double *P, unsigned int N, double EPS, unsigned int istart, unsigned int ppG, double plummer_core, double plummer_mass, double rscale, double mscale);

__device__ void AddPlummer(double4 dr, float4 dv, float4 da, double4 *acc, double4 *jrk, double4 *snp, double *b, double *M, double mul);
__device__ void AddGalacticEnergy(double4 dr, double *pot, double Mscale, double Rscale, double mul);
__device__ void AddGalacticBulge(double4 dr, float4 dv, float4 da, double4 *acc, double4 *jrk, double4 *snp, float M, float b, double mul);
__device__ void AddGalacticHalo(double4 dr, float4 dv, float4 da, double4 *acc, double4 *jrk, double4 *snp, float M3, float a3, double mul);
__device__ void AddGalacticDisk(double4 dr, float4 dv, float4 da, double4 *acc, double4 *jrk, double4 *snp, float M2, float a2, float b2, double mul);
__device__ void AddAllenSantillan(double4 dr, float4 dv, float4 da, float Mscale, float Rscale, double4 *acc, double4 *jrk, double4 *snp, double mul);

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
									 double EPS,
									 double pl_core,
									 double pl_mass,
									 double rscale,
									 double mscale);
#endif
