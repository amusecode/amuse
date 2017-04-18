#pragma once
#include "common.hpp"
#include "mathaux.hpp"
namespace etics {
    namespace scf {
        void InitializeCache(int N);
        void UpdateN(int N);
        __global__ void LoadParticlesToCache(Particle *P, int N);
        __global__ void CalculatePhi0l(int l);
        __global__ void CalculateCoefficientsPartial(int n, int l, Complex *PartialSum);
        void CalculateCoefficients(int n, int l, Complex *A_h);
        void CalculateCoefficients(Complex *A_h);
        template<int Mode> __device__ void CalculateGravityTemplate(int i, Complex *A, vec3 *F, Real *Potential);
        __global__ void CalculateGravityFromCoefficients(Real *Potential, vec3 *F);
        void SendCoeffsToGPU(Complex *A_h);
        void CalculateGravity(Particle *P, int N, Real *Potential, vec3 *F);
        void Init(int N, int k3gs_new, int k3bs_new, int k4gs_new, int k4bs_new);
        void GuessLaunchConfiguration(int N, int *k3gs_new, int *k3bs_new, int *k4gs_new, int *k4bs_new);
    }

    struct CacheStruct {
        int N;
        Real *xi;
        Real *Phi0l;
        Real *Wprev1;
        Real *Wprev2;
        Real *costheta;
        Real *sintheta_I;
        Complex *Exponent;
        Real *mass;
    };
}
