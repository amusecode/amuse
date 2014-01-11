#ifndef MPI_INTERFACE_H
#define MPI_INTERFACE_H

#ifndef NOMPI
#include"mpi.h"
#else

#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

typedef int MPI_Comm;
#define MPI_COMM_WORLD 1

inline double MPI_Wtime() {
#ifdef WIN32
    FILETIME filetime;
    GetSystemTimeAsFileTime(&filetime);
    unsigned long long longtime = filetime.dwHighDateTime;
    longtime <<=32;
    longtime |= filetime.dwLowDateTime;
    longtime /=10;
    longtime -= 11644473600000000ULL;
    return longtime / 1000000.0;
#else
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec;
#endif
}

#endif

//#include"/data/iwasawa/work/amuse/prerequisites/include/mpi.h"
#include<iostream>
#include<fstream>
#include<cmath>
#include<ctime>
//#include<unistd.h>
#include<typeinfo>
#include<cstring>
#include"Particle.h"

void allocate_mpi_buffer(const int &N);

void divide_proc(int Narray[], 
		 int Ndisp[], 
		 const int &Nprt,
		 const int &Nproc);

void divide_proc(const int &Nprt,
		 int &Njp_loc,
		 int &first_address);

void reduce_force(Particle prt[],
		  int address[],
		  const int &Ni0,
		  const int &Ni,
		  const MPI_Comm &icomm);

template<class T>  inline void bcast(T *x,
				     const int &nwords,
				     const MPI_Comm &comm){

#ifndef NOMPI
  MPI_Bcast(x, nwords*sizeof(T), MPI_BYTE, 0, comm);
#endif
}

void sum_double(double &x, 
		double &xsum,
		const int &nwords,
		const MPI_Comm &comm);

void mpi_abort();

#ifndef NOMPI
inline int get_mpi_rank(){
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  return myrank;
}
#else
inline int get_mpi_rank(){
  return 0;
}
#endif

#endif //MPI_INTERFACE_H
