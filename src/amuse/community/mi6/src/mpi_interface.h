#ifndef MPI_INTERFACE_H
#define MPI_INTERFACE_H

#include<iostream>
#include<fstream>
#include<cmath>
#include<ctime>
//#include<unistd.h>
#include<typeinfo>
#include<cstring>
#include"mpi.h"
//#include"/data/iwasawa/work/amuse/prerequisites/include/mpi.h"
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
  MPI_Bcast(x, nwords*sizeof(T), MPI_BYTE, 0, comm);
}

void sum_double(double &x, 
		double &xsum,
		const int &nwords,
		const MPI_Comm &comm);

void mpi_abort();

inline int get_mpi_rank(){
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  return myrank;
}

#endif //MPI_INTERFACE_H
