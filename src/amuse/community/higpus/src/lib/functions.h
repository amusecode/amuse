#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "types.h"
#include <vector>
#include <string>

#include <mpi.h>

using namespace std;

HostError getEnergyLine_HiGPUslog(string *lastLine, string path);
HostError AcquireEnergy(double *E, string path);
HostError CheckHDDMemory(bool *cleanstop __attribute__((unused)), string path);
HostError NextParticles(unsigned int N, unsigned int ompthreads, unsigned int* counter, unsigned int* vetint, double ATIME, double *local_time, double* step, int* next, unsigned long* nextsize);

HostError CudaInit(unsigned int *M, int NGPU, int rank, string gpu_name, const bool setdev, vector<unsigned int> &dev, string path);

HostError Corrector(double *GTIME, double *ATIME, double *local_time, double *step, int *next, unsigned long nextsize, double4 *pos_CH, double4 *vel_CH, double4 *a_H1, double4 *a_H0, double4 *a3_H, double4 *p_v_a3, double ETA6, double ETA4, double DTMAX, double DTMIN, unsigned int N);

HostError append_file(string argument, string path);

HostError open_file(string argument, string path);
HostError print_info(double plc __attribute__((unused)), double plm __attribute__((unused)), double rs __attribute__((unused)), double ms __attribute__((unused)), string path);
HostError adjust_param_ifwarmstart(bool *CDM, bool *CDV, unsigned *FMAX, string warm_start_file, string path, string *FINP, double *GTIME_WARM, double DTPRINT);
HostError __MPIbcast_otherparams(bool *warm_start, bool *VIR, bool *setdev, double *plummer_core, double *plummer_mass, double *rscale, double *mscale, double *ratio, vector<unsigned> &dev);

extern "C"
HostError Hermite6th(const double TTIME, double* GTIME, double* ATIME, double* local_time, double* step, const unsigned int N, const unsigned int M, double4* pos_PH, float4* vel_PH, double4* pos_CH, double4* vel_CH, double4* a_H0, const unsigned int MAXDIM, unsigned int NGPU, unsigned int TPB, int rank, int size, unsigned int BFMAX, double ETA6, double ETA4, double DTMAX, double DTMIN, double EPS, double DTPRINT, unsigned int FMAX, const bool warm, double GTW, unsigned int GPUMINTHREADS, double plummer_core, double plummer_mass, double rscale, double mscale, vector<unsigned int> devices, bool *cleanstop, string path);

extern "C"
HostError InitBlocks(double4 *pos_PH, float4 *vel_PH, unsigned int TPB, unsigned int N, unsigned int M, unsigned int BFMAX, double ETA4, double DTMIN, double DTMAX, unsigned int NGPU, double EPS, unsigned int *MAXDIM, unsigned int *GPUMINTHREADS, string gpu_name, int rank, int nodes, double4* pos_CH, double4* vel_CH, double4* a_H0, double* step, double* local_time, double* ACTUAL_TIME, const bool vir, const double ratio, const bool warm, const bool setdev, vector<unsigned int>& devices, double plummer_core, double plummer_mass, double rscale, double mscale, string path);

HostError Calculate_Energy(double4 **pos_CD, double4 **vel_CD, unsigned int N, double EPS, unsigned int TPB, unsigned int NGPU, int rank, unsigned int ppG, double *kin, double *pot, double plummer_core, double plummer_mass, double rscale, double mscale, vector<unsigned int> devices);

HostError Calculate_potential_Energy(double4 *pos_CH, double4 *vel_CH, unsigned int N, double EPS, unsigned int TPB, unsigned int NGPU, int rank, unsigned int *devices, unsigned int ppG, double *Energy, double plummer_core, double plummer_mass);

HostError isDivisible(unsigned int *N, unsigned int *M, int size, unsigned int NGPU, unsigned int TPB, unsigned int *BFMAX);

HostError ReduceAll(unsigned int cpy_size, unsigned int N, unsigned int NGPU, unsigned long nextsize, double4* a_H, double4 *a_H0, double *mpi_red_aux, double *mpi_red, int *);

HostError DetermineSteps(double stp, unsigned int N, unsigned int M, double4* a_H0, double4* a_H1, double ETA4, double DTMAX, double *step, double *ACTUAL_TIME, double *local_time);
HostError Max_dimension(unsigned int TPB, unsigned int BFMAX, unsigned int N, unsigned int *dim_max, unsigned int GPUMINTHREADS);

HostError check_argv(int ac, char *av[], string *param, bool *warm_start, string *warm_start_file, bool *vir, double *ratio, bool *setdev, vector<unsigned int>& dev, double *plummer_core, double *plummer_mass, double *rscale, double *mscale, string *argument, bool *cleanstop);

HostError CPU_memcheck(const char *file, const int line, string path);
HostError CheckBlocks(double *step, unsigned int N, string path);

HostError cpu_read_external(const string file_name, double4 *pos, float4 *vel, double4 *, const unsigned int N, const unsigned int M, const bool CDM, const bool CDV, const bool);


HostError cpu_read_params(

      const string file_to_read,
      unsigned int *N,
      unsigned int *gpus,
      unsigned int *threads,
      double       *totaltime,
      double       *dtmax,
      double       *dtmin,
      double       *epsilon,
      double       *eta6,
      double       *eta4,
      double       *dtprint,
      unsigned int *filemax,
      bool         *cdm,
      bool         *cdv,
      string       *file,
      string       *gpu_name,
		string       path);

HostError __MPIstart(int *rank, int *size, MPI_Datatype *mpi_float4, MPI_Datatype *mpi_double4);
HostError __MPIbcastParam(unsigned int *N, unsigned int *M, unsigned int *NGPU, unsigned int *TPB, double *TTIME,  double *DTMAX, double *DTMIN, double *EPS, double *ETA6, double *ETA4, double *DTPRINT, unsigned int *BMAX, string *gpu_name, int rank, int size);

HostError random_check_Bcast( int size, int rank, const char* Format ... );
HostError __MPIbcastPosVel(double4 *pos, float4 *vel, double4*, unsigned int N, int rank, int size);
#endif
