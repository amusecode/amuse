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
#include <unistd.h> // for gethostname

namespace etics {
    namespace scf {
        extern __constant__ CacheStruct Cache;
        extern              Complex *PartialSum;
        extern Complex A_h[(NMAX+1)*(LMAX+1)*(LMAX+2)/2];
        extern Complex *PartialSum_h;
        extern int k3gs, k3bs, k4gs, k4bs;

        int blockSizeToDynamicSMemSize(int BlockSize);
        void TestK3(Particle *ParticleList, int N, int numberoftries, double *Average, double *StandardDeviation, bool *Success);
        void TestK4(Particle *ParticleList, int N, int numberoftries, double *Average, double *StandardDeviation, bool *Success);
        void OptimizeLaunchConfiguration(int N);
    }
}

double *Potential;
vec3 *F;
vec3 FirstParticleForce;
const double A000 = 9*(1-0.75*log(3)); // ~1.584, the theoretical A000 coefficient for a Hernquist sphere with a=1/3, which is what we get after noramalizing our initial conditions to Henon units.
const double ExpectedForceX = -1e-5;

void etics::scf::GuessLaunchConfiguration(int N, int *k3gs_new, int *k3bs_new, int *k4gs_new, int *k4bs_new) {
    int blockSize;
    int minGridSize;
    int gridSize;
    cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize, CalculateCoefficientsPartial, blockSizeToDynamicSMemSize, 128);
    cerr << "Warning: setting blockSizeLimit=128 for cudaOccupancyMaxPotentialBlockSizeVariableSMem." << endl;
    gridSize = minGridSize;
    *k3gs_new = gridSize;
    *k3bs_new = blockSize;

    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, CalculateGravityFromCoefficients, 0, N);
    gridSize = (N + blockSize - 1) / blockSize;
    *k4gs_new = gridSize;
    *k4bs_new = blockSize;
}


void etics::scf::TestK3(Particle *ParticleList, int N, int numberoftries, double *Average, double *StandardDeviation, bool *Success) {
    double Average_tmp=0, StandardDeviation_tmp = 0;
    for (int k=0; k<numberoftries; k++) {
        LoadParticlesToCache<<<128,128>>>(ParticleList, N); // need to clear the cache
        DeviceTimer Timer;
        Timer.Start();

        CalculateCoefficients(A_h);

        Timer.Stop();
        double Milliseconds = Timer.Difference()*1000;
        Average_tmp += Milliseconds;
        StandardDeviation_tmp += Milliseconds*Milliseconds;
    }
    Average_tmp /= numberoftries;
    StandardDeviation_tmp = sqrt(StandardDeviation_tmp/numberoftries - Average_tmp*Average_tmp);
    *Average = Average_tmp;
    *StandardDeviation = StandardDeviation_tmp;
    double A000 = 9*(1-0.75*log(3)); // ~1.584, the theoretical A000 coefficient for a Hernquist sphere with a=1/3, which is what we get after noramalizing our initial conditions to Henon units.
    *Success = (0.8 < A_h[0].x/A000) && (A_h[0].x/A000 < 1.2); // very rough success criterion.
}


void etics::scf::TestK4(Particle *ParticleList, int N, int numberoftries, double *Average, double *StandardDeviation, bool *Success) {
    // Need to make sure A_h is loaded to GPU and the first particle is at (187.79, 187.79, 0); also, global arrays Potential and F should be allocated on GPU
    double Average_tmp=0, StandardDeviation_tmp = 0;
    SendCoeffsToGPU(A_h);
    cudaMemset(F, 0, sizeof(vec3));

    for (int k=0; k<numberoftries; k++) {
        LoadParticlesToCache<<<128,128>>>(ParticleList, N); // need to clear the cache
        DeviceTimer Timer;
        Timer.Start();

        CalculateGravityFromCoefficients<<<k4gs,k4bs>>>(Potential, F);

        Timer.Stop();
        double Milliseconds = Timer.Difference()*1000;
        Average_tmp += Milliseconds;
        StandardDeviation_tmp += Milliseconds*Milliseconds;
    }
    Average_tmp /= numberoftries;
    StandardDeviation_tmp = sqrt(StandardDeviation_tmp/numberoftries - Average_tmp*Average_tmp);
    *Average = Average_tmp;
    *StandardDeviation = StandardDeviation_tmp;
    cudaMemcpy(&FirstParticleForce, F, sizeof(vec3), cudaMemcpyDeviceToHost);
    *Success = (0.8 < FirstParticleForce.x/ExpectedForceX) && (FirstParticleForce.x/ExpectedForceX < 1.2); // very rough success criterion.
}

void etics::scf::OptimizeLaunchConfiguration(int N) {
    cout << "We are going to try to optimize the launch configuration for the main ETICS (SCF) kernels by a brute force search." << endl << endl;

    cudaDeviceProp DeviceProperties;
    cudaGetDeviceProperties(&DeviceProperties, 0); // should be DevID!!!
    char HostName[256];
    gethostname(HostName, 256);
    cout << "Probing device GPU" << 0 << " on host " << HostName << ": " << DeviceProperties.name << endl << endl;

    const int ComputeCapability = ((DeviceProperties.major << 4) + DeviceProperties.minor);
    int CorePerSM = 0;
    switch (ComputeCapability) { // We count FP32 cores here... not exactly what we want.
        case 0x10 : CorePerSM = 8;   break; // Tesla   Generation (SM 1.0) G80   class
        case 0x11 : CorePerSM = 8;   break; // Tesla   Generation (SM 1.1) G8x   class
        case 0x12 : CorePerSM = 8;   break; // Tesla   Generation (SM 1.2) G9x   class
        case 0x13 : CorePerSM = 8;   break; // Tesla   Generation (SM 1.3) GT200 class
        case 0x20 : CorePerSM = 32;  break; // Fermi   Generation (SM 2.0) GF100 class
        case 0x21 : CorePerSM = 48;  break; // Fermi   Generation (SM 2.1) GF10x class
        case 0x30 : CorePerSM = 192; break; // Kepler  Generation (SM 3.0) GK10x class
        case 0x32 : CorePerSM = 192; break; // Kepler  Generation (SM 3.2) GK10x class
        case 0x35 : CorePerSM = 192; break; // Kepler  Generation (SM 3.5) GK11x class
        case 0x37 : CorePerSM = 192; break; // Kepler  Generation (SM 3.7) GK21x class
        case 0x50 : CorePerSM = 128; break; // Maxwell Generation (SM 5.0) GM10x class
        case 0x52 : CorePerSM = 128; break; // Maxwell Generation (SM 5.2) GM20x class
        case 0x60 : CorePerSM = 64; break;  // Pascal  Generation (SM 6.0) GP10x class
    }

    int Cores = CorePerSM * DeviceProperties.multiProcessorCount ;
    if (Cores == 0)  {
        Cores = 3584;
        cout << "Could not count cores! Your GPU is possibly too new. We'll take " << Cores << " but it doesn't really matter." << endl << endl;
    }

    const int RealMaxGS = Cores/2;
    cout << "We are only going to examine grids smaller than " << RealMaxGS << " blocks, for absolutely no good reason." << endl << endl;

    const int numberoftries = 10;

    const int WarpSize = DeviceProperties.warpSize, MinBS=WarpSize;
    const int ShmemSizePerBlock = sizeof(Complex)*(LMAX+1); // for Kernel3
    int MaxBS = (int)(DeviceProperties.sharedMemPerBlock / (WarpSize*ShmemSizePerBlock))*WarpSize;
    MaxBS = (MaxBS>DeviceProperties.maxThreadsPerBlock)?(DeviceProperties.maxThreadsPerBlock):(MaxBS);

    Particle *ParticleList, *ParticleList_h;
    cout << "Generating initial conditions (this may take a while)." << endl << endl;
    etics::ic::hernquist(N, 0, &ParticleList_h);
    ParticleList_h[0].pos = vec3(187.79445239392416256476, 187.79445239392416256476, 0);
    cudaMalloc((void**)&ParticleList, N * sizeof(Particle));
    cudaMemcpy(ParticleList, ParticleList_h, N * sizeof(Particle), cudaMemcpyHostToDevice);
    cudaMalloc((void**)&Potential, N * sizeof(Real));
    cudaMalloc((void**)&F, N * sizeof(vec3));

    int k3gs_tmp, k3bs_tmp, k4gs_tmp, k4bs_tmp;
    GuessLaunchConfiguration(N, &k3gs_tmp, &k3bs_tmp, &k4gs_tmp, &k4bs_tmp);
    printf("Recommended launch configuration for Kernel3 (CalculateCoefficientsPartial): <<<%d,%d>>>\n", k3gs_tmp, k3bs_tmp);
    Init(N, k3gs_tmp, k3bs_tmp, k4gs_tmp, k4bs_tmp);
    double Average=0, StandardDeviation=0;
    cout << "Testing..." << endl;
    bool Success;
    TestK3(ParticleList, N, numberoftries, &Average, &StandardDeviation, &Success);

    printf("Executed in %.2f ms +/- %.2f\n\n", Average, StandardDeviation);
    double k3_normal_time = Average;


    printf("Recommended launch configuration for Kernel4 (CalculateGravityFromCoefficients): <<<%d,%d>>>\n", k4gs_tmp, k4bs_tmp);
    cout << "Testing..." << endl;
    TestK4(ParticleList, N, numberoftries, &Average, &StandardDeviation, &Success);
    double k4_normal_time = Average;

    printf("Executed in %.2f ms +/- %.2f\n\n", Average, StandardDeviation);


    free(PartialSum_h);
    cudaFree(PartialSum);
    PartialSum_h = (Complex*)malloc(RealMaxGS*(LMAX+1)*sizeof(Complex)); // why not use "new"?
    cudaMalloc((void**)&PartialSum, RealMaxGS*(LMAX+1)*sizeof(Complex));

    int TotalTests = RealMaxGS * (MaxBS/32);
    double Average_arr[TotalTests], StandardDeviation_arr[TotalTests];
    int BlockSize_arr[TotalTests], GridSize_arr[TotalTests];

    cout << "Optimizing K3 (block size is a power of two due to summation algorithm)" << endl;
    int i = 0;
    Success = true;
    for (k3bs = MinBS; k3bs <= MaxBS; k3bs *= 2) {
        if (!Success) break;
        int MinGS = (Cores/k3bs>0)?(Cores/k3bs):1; // honestly can't remember why
        int MaxGS = (N+k3bs-1)/k3bs;
        MaxGS = (MaxGS>RealMaxGS)?(RealMaxGS):(MaxGS);
        for (k3gs = MinGS; k3gs <= MaxGS; k3gs++) {
            TestK3(ParticleList, N, numberoftries, &Average, &StandardDeviation, &Success);
            if (!Success) break;
            printf("<<<%04d,%03d>>> %7.2f %7.2f %.3e\n", k3gs, k3bs, Average, StandardDeviation, A_h[0].x/A000);
            fflush(stdout);
            BlockSize_arr[i]=k3bs; GridSize_arr[i]=k3gs;
            Average_arr[i] = Average;
            StandardDeviation_arr[i] = StandardDeviation;
            i++;
        }
    }

    int MaxIndex = i;
    int IndexOfMininum = 0;
    for (int i = 0; i < MaxIndex; i++) if (Average_arr[i] < Average_arr[IndexOfMininum]) IndexOfMininum = i;
    int k3gs_opt=GridSize_arr[IndexOfMininum], k3bs_opt=BlockSize_arr[IndexOfMininum];
    double k3_opt_time = Average_arr[IndexOfMininum];
    printf("Fastest configuration is: <<<%04d,%03d>>>\n", GridSize_arr[IndexOfMininum], BlockSize_arr[IndexOfMininum]);
    cout << "Other options:" << endl;
    for (int i = 0; i < MaxIndex; i++) {
        if ((Average_arr[i]-StandardDeviation_arr[i] < Average_arr[IndexOfMininum]+StandardDeviation_arr[IndexOfMininum]) && (i!=IndexOfMininum)) {
            printf("    <<<%04d,%03d>>>\n", GridSize_arr[i], BlockSize_arr[i]);
        }
    }

    cout << "Optimizing K4" << endl;
    Init(N, k3gs_tmp, k3bs_tmp, k4gs_tmp, k4bs_tmp);
    TestK3(ParticleList, N, 1, &Average, &StandardDeviation, &Success);

    i = 0;
    Success = true;
    for (k4bs = MinBS; k4bs <= DeviceProperties.maxThreadsPerBlock; k4bs += WarpSize) {
        if (!Success) break;
        int MinGS = (Cores/k4bs>0)?(Cores/k4bs):1; // honestly can't remember why
        int MaxGS = (N+k4bs-1)/k4bs;
        MaxGS = (MaxGS>RealMaxGS)?(RealMaxGS):(MaxGS);
        for (k4gs = MinGS; k4gs <= MaxGS; k4gs++) {
            TestK4(ParticleList, N, numberoftries, &Average, &StandardDeviation, &Success);
            if (!Success) break;
            printf("<<<%04d,%03d>>> %7.2f %7.2f %.3e\n", k4gs, k4bs, Average, StandardDeviation, FirstParticleForce.x/ExpectedForceX);
            fflush(stdout);
            BlockSize_arr[i]=k4bs; GridSize_arr[i]=k4gs;
            Average_arr[i] = Average;
            StandardDeviation_arr[i] = StandardDeviation;
            i++;
        }
    }

    MaxIndex = i;
    IndexOfMininum = 0;
    for (int i = 0; i < MaxIndex; i++) if (Average_arr[i] < Average_arr[IndexOfMininum]) IndexOfMininum = i;
    int k4gs_opt=GridSize_arr[IndexOfMininum], k4bs_opt=BlockSize_arr[IndexOfMininum];
    double k4_opt_time = Average_arr[IndexOfMininum];
    printf("Fastest configuration is: <<<%04d,%03d>>>\n", GridSize_arr[IndexOfMininum], BlockSize_arr[IndexOfMininum]);
    cout << "Other options:" << endl;
    for (int i = 0; i < MaxIndex; i++) {
        if ((Average_arr[i]-StandardDeviation_arr[i] < Average_arr[IndexOfMininum]+StandardDeviation_arr[IndexOfMininum]) && (i!=IndexOfMininum)) {
            printf("    <<<%04d,%03d>>>\n", GridSize_arr[i], BlockSize_arr[i]);
        }
    }

    printf("===================================  SUMMARY  ==================================\n");
    printf("Parameters: LMAX=%d, NMAX=%d, N=%d (GPU%d=\"%s\" on %s)\n", LMAX, NMAX, N, 0, DeviceProperties.name, HostName);
    printf("Recommended launch configuration for K3: <<<%d,%d>>>; execution time: %.2f ms.\n", k3gs_tmp, k3bs_tmp, k3_normal_time);
    printf("Optimal launch configuration for K3: <<<%d,%d>>>; execution time: %.2f ms.\n", k3gs_opt, k3bs_opt, k3_opt_time);
    printf("Recommended launch configuration for K4: <<<%d,%d>>>; execution time: %.2f ms.\n", k4gs_tmp, k4bs_tmp, k4_normal_time);
    printf("Optimal launch configuration for K4: <<<%d,%d>>>; execution time: %.2f ms.\n", k4gs_opt, k4bs_opt, k4_opt_time);
    printf("================================================================================\n");

    cudaFree(PartialSum);
    cudaFree(ParticleList);
    cudaFree(Potential);
    cudaFree(F);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Please specify number of particles (> 10000)." << endl;
        return 1;
    }
    int N = atoi(argv[1]);
    if (N < 10000) {
        cout << "Please use more than 10000 particles." << endl;
        return 1;
    }
    etics::scf::OptimizeLaunchConfiguration(N);
    return 0;
}