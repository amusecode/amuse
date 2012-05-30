#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <types.h>
#include <vector>
#include <string>

#include <mpi.h>

using namespace std;

extern "C"
HostError InitBlocks(double4 *pos, float4 *vel, unsigned int TPB, unsigned int N, unsigned int M, unsigned int BFMAX, double ETA4, unsigned int NGPU, unsigned int *MAXDIM, double DTMAX, double DTMIN, unsigned int *GPUMINTHREADS, unsigned int *devices, string gpu_name, int rank, int nodes, double4* pos_CH, double4* vel_CH, double4* a_H0, double* step, double *local_time, double* ACTUAL_TIME, double plummer_core, double plummer_mass, string path);

HostError CudaInit(unsigned int *M, int NGPU, int rank, unsigned int *devices, string gpu_name, string path);

HostError Calculate_Energy(double4 **pos_CD, double4 **vel_CD, float4 **vel_PD, unsigned int N, unsigned int TPB, unsigned int NGPU, int rank, unsigned int *devices, unsigned int ppG, double *Energy, double plummer_core, double plummer_mass);

HostError Calculate_potential_Energy(double4 *pos_CH, float4 *vel_PH, unsigned int N, unsigned int TPB, unsigned int NGPU, int rank, unsigned int *devices, unsigned int ppG, double *Energy, double plummer_core, double plummer_mass);

HostError isDivisible(unsigned int *N, unsigned int *M, int size, unsigned int NGPU, unsigned int TPB, unsigned int *BFMAX);

HostError ReduceAll(unsigned int cpy_size, unsigned int N, unsigned int NGPU, unsigned long nextsize, double4* a_H, double4 *a_H0, double *mpi_red_aux, double *mpi_red, int *);

HostError DetermineSteps(double stp, unsigned int N, unsigned int M, double4* a_H0, double4* a_H1, double ETA4, double DTMAX, double *step, double *ACTUAL_TIME, double *local_time);
HostError Max_dimension(unsigned int TPB, unsigned int BFMAX, unsigned int N, unsigned int *dim_max, unsigned int GPUMINTHREADS);
HostError check_argv(int ac, char *av[], std::string *param, bool *, std::string *, double *plummer_core, double *plummer_mass);
HostError CPU_memcheck(const char *file, const int line, string path);
HostError CheckBlocks(double *step, unsigned int N, string path);

extern "C"
HostError Hermite6th(const double TTIME, double* GTIME, double* ATIME, double* local_time, double* step, const unsigned int N, const unsigned int M, double4* pos_PH, float4* vel_PH, double4* pos_CH, double4* vel_CH, double4* a_H0, unsigned int MAXDIM, unsigned int NGPU, unsigned int *devices, unsigned int TPB, int rank, int size, unsigned int BFMAX, double ETA6, double ETA4, double DTMAX, double DTMIN, double DTPRINT, unsigned int FMAX, const bool warm, double GTIME_WARM, unsigned int GPUMINTHREADS, double plummer_core, double plummer_mass, string path);

HostError cpu_read_external(const string file_name, double4 *pos, float4 *vel, double4 *, const unsigned int N, const unsigned int M, const bool CDM, const bool CDV, const bool);

HostError cpu_read_params(
      const std::string file_to_read,
      unsigned int *N,
      unsigned int *gpus,
      unsigned int *threads,
      double       *totaltime,
      double       *dtmax,
      double       *dtmin,
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
HostError __MPIbcastParam(unsigned int *N, unsigned int *NGPU, unsigned int *TPB, double *TTIME,  double *DTMAX, double *DTMIN, double *ETA6, double *ETA4, double *DTPRINT, unsigned int *BMAX, string *gpu_name, int rank, int size);

HostError random_check_Bcast( int size, int rank, const char* Format ... );
HostError __MPIbcastPosVel(double4 *pos, float4 *vel, double4*, unsigned int N, int rank, int size);
#endif
