/**
 * @file
 * @author  Yohai Meiron <ymeiron@pku.edu.cn>
 * @brief   Functions to calculate gravitational force using the SCF method.
 */
#include "common.hpp"
#include "mathaux.hpp"
#include "scf.hpp"
#include "ic.hpp"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <mpi.h>

namespace etics {
    namespace scf {
        __constant__ Real RadCoeff[(NMAX+1)*(LMAX+1)];             /*!< Stores coefficients related to the G */
        __constant__ Real AngCoeff[(LMAX+1)*(LMAX+2)/2];           /*!< used blab bla222 */
        __constant__ Complex A[(NMAX+1)*(LMAX+1)*(LMAX+2)/2];
        __constant__ CacheStruct Cache;
                     Complex *PartialSum;

        Real RadCoeff_h[(NMAX+1)*(LMAX+1)];        /*!< Stores coefficients related to the G */
        Real AngCoeff_h[(LMAX+1)*(LMAX+2)/2];      /*!< used blab bla222 */
        Complex A_h[(NMAX+1)*(LMAX+1)*(LMAX+2)/2];
        CacheStruct Cache_h;
        Complex *PartialSum_h;

        int k3gs, k3bs, k4gs, k4bs;
    }
}

void etics::scf::InitializeCache(int N) { // not sure why it's a separate function, the instructions can be in etics::scf::Init()
    Cache_h.N = N;
    cudaMalloc((void**)&Cache_h.xi,         N * sizeof(Real));
    cudaMalloc((void**)&Cache_h.Phi0l,      N * sizeof(Real));
    cudaMalloc((void**)&Cache_h.Wprev1,     N * sizeof(Real));
    cudaMalloc((void**)&Cache_h.Wprev2,     N * sizeof(Real));
    cudaMalloc((void**)&Cache_h.costheta,   N * sizeof(Real));
    cudaMalloc((void**)&Cache_h.sintheta_I, N * sizeof(Real));
    cudaMalloc((void**)&Cache_h.Exponent,   N * sizeof(Complex));
    cudaMalloc((void**)&Cache_h.mass,       N * sizeof(Real));
}

void etics::scf::UpdateN(int N) {
    Cache_h.N = N;
    cudaMemcpyToSymbol(Cache, &Cache_h, sizeof(CacheStruct));
}

__global__ void etics::scf::LoadParticlesToCache(Particle *P, int N) { // formerly "Kernel1"
    int i = threadIdx.x + blockIdx.x *  blockDim.x;
    while (i < N) {
        vec3 Pos = P[i].pos;
        Real r = sqrt(Pos.x*Pos.x + Pos.y*Pos.y + Pos.z*Pos.z);
        Real xi = (r-1)/(r+1);
        Real costheta = Pos.z/r;
        Real sintheta_I = rsqrt(1-costheta*costheta);

        Cache.xi[i] = xi;
        Cache.Phi0l[i] = 0.5 * (1 - xi);
        Cache.costheta[i] = costheta;
        Cache.sintheta_I[i] = sintheta_I;

        Real Normal_I = rsqrt(Pos.x*Pos.x + Pos.y*Pos.y);
        Complex Exponent = make_Complex(Pos.x*Normal_I, -Pos.y*Normal_I);
        Cache.Exponent[i] = Exponent;
        Cache.mass[i] = P[i].m;

        i += blockDim.x * gridDim.x;
    }
}

__global__ void etics::scf::CalculatePhi0l(int l) { // formerly "Kernel2"
    int i = threadIdx.x + blockIdx.x *  blockDim.x;
    while (i < Cache.N) {
        Real xi = Cache.xi[i];
        Cache.Phi0l[i] *= 0.25*(1-xi*xi);

        i += blockDim.x * gridDim.x;
    }
}

__global__ void etics::scf::CalculateCoefficientsPartial(int n, int l, Complex *PartialSum) { // formerly "Kernel3"
    extern __shared__ Complex ReductionCache[]; // size is determined in kernel launch
    int tid = threadIdx.x;
    for (int m = 0; m <= l; m++) ReductionCache[m*blockDim.x+tid] = make_Complex(0, 0);
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < Cache.N) {
        Real xi = Cache.xi[i];
        Real Wnl;
        if (n == 0)      Wnl = 1;
        else if (n == 1) {Wnl = (4*l+3)*xi; Cache.Wprev2[i] = Wnl;}
        else if (n == 2) {Wnl = -(2*l+1.5)+( 8*l*(l+2) +7.5)*xi*xi; Cache.Wprev1[i] = Wnl;}
        else {
            Real Wprev1 = Cache.Wprev1[i];
            Wnl = (xi*(2*n+4*l+1)*Wprev1 - (n+4*l+1)*Cache.Wprev2[i])/(Real)n;
            if (n < NMAX) { // Writing is expensive, avoid if possible.
                Cache.Wprev2[i] = Wprev1;
                Cache.Wprev1[i] = Wnl;
            }
        }
        Real RadialPart = - Cache.mass[i] * SQRT_4_PI * Cache.Phi0l[i] * Wnl * RadCoeff[(LMAX+1)*n+l];
        Real costheta = Cache.costheta[i];
        Real Plm = Pl(l, costheta);
        ReductionCache[tid] = Complex_add(ReductionCache[tid], make_Complex(RadialPart * Plm * AngCoeff[(l+1)*l/2],0));
        if (l == 0) {i += blockDim.x * gridDim.x; continue;}

        //////////////////////////////// ugly fix
        if ((costheta < -0.999) || (costheta > +0.999)) {
            i += blockDim.x * gridDim.x;
            continue;
        }
        //////////////////////////////// ugly fix

        Real Plm_prev1 = Plm;
        Real sintheta_I = Cache.sintheta_I[i];
        Plm = (costheta*Plm - Pl(l-1, costheta))*l*sintheta_I;
        Complex Exponent = Cache.Exponent[i];
        Real tmp0 = RadialPart * Plm * AngCoeff[(l+1)*l/2+1];
        ReductionCache[blockDim.x+tid] = Complex_add(ReductionCache[blockDim.x+tid], make_Complex(tmp0 * Exponent.x, tmp0 * Exponent.y));

        if (l == 1) {i += blockDim.x * gridDim.x; continue;}

        Complex TorodialPart = Exponent;
        for (int m = 2; m <= l; m++) { // make sure no redundancy at the end of the loop
            Real Plm_prev2 = Plm_prev1;
            Plm_prev1 = Plm;
            Plm = - 2*(m-1)*costheta*sintheta_I*Plm_prev1 - (l+m-1)*(l-m+2)*Plm_prev2;
            TorodialPart = Complex_mul(TorodialPart, Exponent);
            tmp0 = RadialPart * Plm * AngCoeff[(l+1)*l/2+m];
            ReductionCache[m*blockDim.x+tid] = Complex_add(ReductionCache[m*blockDim.x+tid], make_Complex(tmp0 * TorodialPart.x, tmp0 * TorodialPart.y));
        }
        i += blockDim.x * gridDim.x;
    }
    __syncthreads();
    for (int m = 0; m <= l; m++) {
        i = blockDim.x/2;
        while (i != 0) {
            if (tid < i)
                ReductionCache[m*blockDim.x+tid] = Complex_add(ReductionCache[m*blockDim.x+tid], ReductionCache[m*blockDim.x+tid+i]);
            __syncthreads();
            i /= 2;
        }
        if (tid == 0)
            PartialSum[blockIdx.x*(l+1) + m] = ReductionCache[m*blockDim.x];
    }
}

void etics::scf::CalculateCoefficients(int n, int l, Complex *A_h) {
    int BaseAddress = n*(LMAX+1)*(LMAX+2)/2 + l*(l+1)/2;
    CalculateCoefficientsPartial<<<k3gs,k3bs,k3bs*sizeof(Complex)*(LMAX+1)>>>(n, l, PartialSum);
    cudaMemcpy(PartialSum_h, PartialSum, k3gs*(l+1)*sizeof(Complex), cudaMemcpyDeviceToHost);
    for (int m = 0; m <= l; m++)
        for (int Block=0; Block<k3gs; Block++)
            A_h[BaseAddress + m] = Complex_add(A_h[BaseAddress + m], PartialSum_h[Block*(l+1) + m]);
}

void etics::scf::CalculateCoefficients(Complex *A_h) {
    memset(A_h, 0, (NMAX+1)*(LMAX+1)*(LMAX+2)/2 * sizeof(Complex));
    for (int l = 0; l <= LMAX; l++) {
        if (l > 0) CalculatePhi0l<<<128,128>>>(l); // wouldn't it make sense just putting it after the n-loop finishes? Probably not becasue then we need to skip at the last iter
        for (int n = 0; n <= NMAX; n++) {
            CalculateCoefficients(n, l, A_h);
        }
    }
}

template<int Mode>
__device__ void etics::scf::CalculateGravityTemplate(int i, Complex *A, vec3 *F, Real *Potential) {
// it gets A as parameter because it can be either on host or device
#warning !!! This cannot really be a host function because it needs device cahce, angular coefficients which are on device!!
// 0 = both force and potential, 1 = only force, 2 = only pot
#define A(n,l,m) A[n*(LMAX+1)*(LMAX+2)/2 + l*(l+1)/2 + m]
    Real dPhiLeft;
    Real dPhiLeftMul;
    Real dPhiRight;
    Real dPhiRightAdd;
    Real dPhi;
    Real RadialPart2;
    Real PlmDerivTheta;

    Real Pot = 0;
    Real Fr = 0, Ftheta = 0, Fphi = 0;
    Real xi = Cache.xi[i];
    Real OneOverXiPlusOne = 1/(1+xi);
    Real r_I = (1-xi)*OneOverXiPlusOne;
    Real r = 1/r_I; // It's quite likely we can do everything without r.
    Real costheta = Cache.costheta[i];
    Real sintheta_I = rsqrt(1-costheta*costheta); // faster than using cachei // You sure??? in K3 it's the opposite

    Complex ExponentTmp[LMAX];
    Complex Exponent = Complex_conj(Cache.Exponent[i]);
    ExponentTmp[0] = Exponent;
    for (int m = 1; m < LMAX; m++) ExponentTmp[m] = Complex_mul(ExponentTmp[m-1],Exponent);
    if (Mode != 2) {
        Real xi2 = xi*xi;
        Real xi3 = xi2*xi;
        dPhiLeft = -0.25*OneOverXiPlusOne;
        dPhiLeftMul = 0.25*(1-xi2);
        dPhiRight = xi3 - xi2 - xi + 1;
        dPhiRightAdd = 2*(xi3 - 2*xi2 + xi);
    }
    Real Phi0l = 1/(1+r);
    Real tmp1 = Phi0l*Phi0l*r;
    for (int l = 0; l <= LMAX; l++) {
        if (Mode != 2) {
            if (l > 0) {
                dPhiLeft  *= dPhiLeftMul;
                dPhiRight += dPhiRightAdd;
            }
        }
        if (Mode != 2) dPhi = dPhiLeft * dPhiRight;
        for (int n = 0; n <= NMAX; n++) {
            Real Wnl, Wprev1, Wprev2;
            if (n == 0)      Wnl = 1;
            else if (n == 1) {Wnl = (4*l+3)*xi; Wprev2 = Wnl;}
            else if (n == 2) {Wnl = -(2*l+1.5)+( 8*l*(l+2) +7.5)*xi*xi; Wprev1 = Wnl;}
            else {
                Wnl = (xi*(2*n+4*l+1)*Wprev1 - (n+4*l+1)*Wprev2)/(Real)n;
                Wprev2 = Wprev1;
                Wprev1 = Wnl;
            }

            Real Wderiv = 0;
            if (n == 1) {Wderiv = 4*l + 3;}
            else if (n > 1) {
                Wderiv = (-n*xi*Wnl + (n+4*l+2)*Wprev2)/(1-xi*xi);
            } // From an unknown reason it's faster to have this Block separate from the previous one.

            Real RadialPart  = - SQRT_4_PI * Phi0l * Wnl;
            if (Mode != 2) RadialPart2 = SQRT_4_PI * (dPhi*Wnl + Phi0l*Wderiv*2/pow(1+r,2));

            Real Plm = Pl(l, costheta);
            Real tmp2 = Complex_real(A(n,l,0)) * AngCoeff[(l+1)*l/2] * Plm;
            if (Mode != 1) Pot += RadialPart  * tmp2;
            if (Mode != 2) Fr  += RadialPart2 * tmp2;

            if (l == 0) continue;

            //////////////////////////////// ugly fix
            if ((costheta < -0.999) || (costheta > +0.999)) {
                continue;
            }
            //////////////////////////////// ugly fix

            // The Block below is l>=1, m=0.
            if (Mode != 2) {
                PlmDerivTheta = (costheta*Plm - Pl(l-1, costheta))*l*sintheta_I; //TODO check if storing Pl(l-1) somewhere makes it faster
                Ftheta += - PlmDerivTheta * AngCoeff[(l+1)*l/2] * Complex_real(A(n,l,0)) * RadialPart * r_I;
            }

            // The Block below is l>=1, m=1.
            if (Mode == 2) PlmDerivTheta = (costheta*Plm - Pl(l-1, costheta))*l*sintheta_I; //TODO see above regarding storing Pl(l-1)
            Real Plm_prev1 = Plm;
            Plm = PlmDerivTheta; // PlmDerivTheta equals Plm for m=1.
            if (Mode != 2) PlmDerivTheta = - Plm*costheta*sintheta_I - l*(l+1)*Plm_prev1;
            tmp2 = 2 * AngCoeff[(l+1)*l/2+1];
            Complex tmp3 = Complex_mul(ExponentTmp[0], A(n,l,1));
            Complex tmp4 = make_Complex(tmp2 * tmp3.x, tmp2 * tmp3.y);
            Complex tmp5 = make_Complex(Plm *  tmp4.x, Plm *  tmp4.y);
            Complex tmp6 = make_Complex(RadialPart * tmp5.x, RadialPart * tmp5.y);
            if (Mode != 1) Pot += Complex_real(tmp6);
            if (Mode != 2) {
                Fr +=       RadialPart2 * Complex_real(tmp5);
                Fphi +=     Complex_imag(tmp6) * sintheta_I * r_I;
                Ftheta += - RadialPart * PlmDerivTheta *  Complex_real(tmp4) * r_I;
            }

            if (l == 1) continue;

            for (int m = 2; m <= l; m++) {
                Real Plm_prev2 = Plm_prev1;
                Plm_prev1 = Plm;
                Plm = - 2*(m-1)*costheta*sintheta_I*Plm_prev1 - (l+m-1)*(l-m+2)*Plm_prev2;
                tmp2 = 2 * AngCoeff[(l+1)*l/2+m];
                tmp3 = Complex_mul(ExponentTmp[m-1], A(n,l,m));
                tmp4 = make_Complex(tmp2 * tmp3.x, tmp2 * tmp3.y);
                tmp5 = make_Complex(Plm *  tmp4.x, Plm *  tmp4.y);
                tmp6 = make_Complex(RadialPart * tmp5.x, RadialPart * tmp5.y);
                if (Mode != 1) Pot  += Complex_real(tmp6);
                if (Mode != 2) {
                    PlmDerivTheta = - m*Plm*costheta*sintheta_I - (l+m)*(l-m+1)*Plm_prev1;
                    Fr +=       RadialPart2 * Complex_real(tmp5);
                    Fphi +=     m * Complex_imag(tmp6) * sintheta_I * r_I;
                    Ftheta += - RadialPart * PlmDerivTheta *  Complex_real(tmp4) * r_I;
                }
            }
        }
        Phi0l *= tmp1;
    }

    if (Mode != 2) {
        Real sintheta = 1/sintheta_I;
        Real tanphi = Exponent.y/Exponent.x;
        Real cosphi = ((Exponent.x >= 0)?(+1):(-1)) * rsqrt(1+tanphi*tanphi); // no simpler way to get sign bit?
        Real sinphi = tanphi*cosphi;
        *F = vec3(sintheta*cosphi*Fr + costheta*cosphi*Ftheta - sinphi*Fphi, sintheta*sinphi*Fr + costheta*sinphi*Ftheta + cosphi*Fphi,   costheta*Fr - sintheta*Ftheta);
    }
    if (Mode != 1) *Potential = Pot;
#undef A
}

__global__ void etics::scf::CalculateGravityFromCoefficients(Real *Potential, vec3 *F) { // formerly "Kernel4"
#define A(n,l,m) A[n*(LMAX+1)*(LMAX+2)/2 + l*(l+1)/2 + m]
#ifdef A_ON_SHARED_MEMORY
    __shared__ Complex Buffer[(NMAX+1)*(LMAX+1)*(LMAX+2)/2];
    if (threadIdx.x < warpSize) {
        for(int i = threadIdx.x; i  < (NMAX+1)*(LMAX+1)*(LMAX+2)/2; i += warpSize) {
            Buffer[i] = A[i];
        }
    }
    __syncthreads();
    #define A(n,l,m) Buffer[n*(LMAX+1)*(LMAX+2)/2 + l*(l+1)/2 + m]
#endif

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < Cache.N) {
        CalculateGravityTemplate<0>(i, A, &F[i], &Potential[i]);
#warning if we have A_ON_SHARED_MEMORY the above won't work
        i += blockDim.x * gridDim.x;
    }
}
#undef A

void etics::scf::SendCoeffsToGPU(Complex *A_h) {
        cudaMemcpyToSymbol(A, A_h, (NMAX+1)*(LMAX+1)*(LMAX+2)/2 * sizeof(Complex));
}

void etics::scf::CalculateGravity(Particle *P, int N, Real *Potential, vec3 *F) {
    LoadParticlesToCache<<<128,128>>>(P, N);
    CalculateCoefficients(A_h);
    Complex ATotal[(NMAX+1)*(LMAX+1)*(LMAX+2)/2];
    MPI_Allreduce(&A_h, &ATotal, (NMAX+1)*(LMAX+1)*(LMAX+2)/2*2, MPI_ETICS_REAL, MPI_SUM, MPI_COMM_WORLD);
    std::copy ( ATotal, ATotal+(NMAX+1)*(LMAX+1)*(LMAX+2)/2, A_h);
#warning not really need this copy, just calculate the coefficients in some local array, then sum it into a global array (A_h or somthing) and copy it to GPUs
//     cudaMemcpyToSymbol(A, A_h, (NMAX+1)*(LMAX+1)*(LMAX+2)/2 * sizeof(Complex));
    SendCoeffsToGPU(A_h);
    CalculateGravityFromCoefficients<<<k4gs,k4bs>>>(Potential, F);
}

namespace etics {
        namespace scf {
        int blockSizeToDynamicSMemSize(int BlockSize) { // Should be a lambda function
            return (LMAX+1)*sizeof(Complex)*BlockSize;
        }
    }
}

void etics::scf::Init(int N, int k3gs_new, int k3bs_new, int k4gs_new, int k4bs_new) {
    if ((k3gs_new==0) || (k3bs_new==0)) {
        cerr << "Warning: launch configuration for CalculateCoefficientsPartial(...) is unspecified; performance can be improved by optimizing it for this device." << endl;
        int blockSize;
        int minGridSize;
        int gridSize;
        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize, CalculateCoefficientsPartial, blockSizeToDynamicSMemSize, 128);
        cerr << "Warning: setting blockSizeLimit=128 for cudaOccupancyMaxPotentialBlockSizeVariableSMem." << endl;
        gridSize = minGridSize;
        cerr << "Using the following launch configuration: <<<" << gridSize << "," << blockSize << ">>>" << endl;
        k3gs = gridSize;
        k3bs = blockSize;
    } else {
        k3gs = k3gs_new;
        k3bs = k3bs_new;
    }

    if ((k4gs_new==0) || (k4bs_new==0)) {
        cerr << "Warning: launch configuration for CalculateGravityFromCoefficients is unspecified; performance can be improved by optimizing it for this device." << endl;
        int blockSize;
        int minGridSize;
        int gridSize;
        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, CalculateGravityFromCoefficients, 0, N);
        gridSize = (N + blockSize - 1) / blockSize;
        cerr << "Using the following launch configuration: <<<" << gridSize << "," << blockSize << ">>>" << endl;
        k4gs = gridSize;
        k4bs = blockSize;
    } else {
        k4gs = k4gs_new;
        k4bs = k4bs_new;
    }

    RadialCoefficients(RadCoeff_h);
    cudaMemcpyToSymbol(RadCoeff, RadCoeff_h, (NMAX+1)*(LMAX+1) * sizeof(Real));
    AngularCoefficients(AngCoeff_h);
    cudaMemcpyToSymbol(AngCoeff, AngCoeff_h, (LMAX+1)*(LMAX+2)/2 * sizeof(Real));
    InitializeCache(N);
    cudaMemcpyToSymbol(Cache, &Cache_h, sizeof(CacheStruct));
    PartialSum_h = (Complex*)malloc(k3gs*(LMAX+1)*sizeof(Complex)); // why not use "new"?
    cudaMalloc((void**)&PartialSum, k3gs*(LMAX+1)*sizeof(Complex));
}
