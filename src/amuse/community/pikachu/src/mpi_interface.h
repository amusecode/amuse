#ifndef MPI_INTERFACE_H
#define MPI_INTERFACE_H

#include<iostream>
#include<fstream>
#include<cmath>
#include<ctime>
#include<unistd.h>
#include<typeinfo>
#include<cstring>
#include"mpi.h"
#include"Vector3.h"
#include"const.h"

inline void make_type_vec3(MPI_Datatype* mpi_vec3){

  int len_blk[3];
  len_blk[0] = len_blk[1] = len_blk[2] = 1;

  MPI_Datatype type[3];
  type[0] =  type[1] = type[2] = MPI_DOUBLE;

  MPI_Aint disp[3];
  disp[0] = 0;
  Vector3 vtmp;
  MPI_Aint add_s, add;
  MPI_Address(&(vtmp[0]), &add_s);
  MPI_Address(&(vtmp[1]), &add);
  disp[1] = add - add_s;
  MPI_Address(&(vtmp[2]), &add);
  disp[2] = add - add_s;

  MPI_Type_struct(3, len_blk, disp, type, mpi_vec3);
  MPI_Type_commit(mpi_vec3);

}

inline int mpi_get_rank(){
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  return myrank;
}

inline int mpi_get_size(){
  int Nproc = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &Nproc);
  return Nproc;
}

inline void mpi_initialize(int argc, 
			   char *argv[], 
			   int *myrank,
			   int *Nproc){
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, Nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, myrank);
  MPI_Barrier(MPI_COMM_WORLD);
  std::cerr << "MPI Initialize: myrank = " << *myrank
	    << "                Nproc = " << *Nproc <<std::endl;

}

inline int mpi_max_int(const int *localval){
  int localval_tmp = *localval;
  int glovalval = 0;
  MPI_Allreduce(&localval_tmp, &glovalval, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
  return glovalval;
}


inline double mpi_max_double(const double *localval){
  double localval_tmp = *localval;
  double glovalval = 0;
  MPI_Allreduce(&localval_tmp, &glovalval, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); 
  return glovalval;
}

inline double mpi_min_double(const double *localval){
  double localval_tmp = *localval;
  double glovalval = 0;
  MPI_Allreduce(&localval_tmp, &glovalval, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); 
  return glovalval;
}

inline int mpi_sum_int_np(int x){
  int xsum;
  int xsent = x;
  MPI_Allreduce(&xsent, &xsum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
  return xsum;
}

inline double mpi_sum_double_np(const double &x){
  double xsum;
  double xsent = x;
  MPI_Allreduce(&xsent, &xsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
  return xsum;
}

template<class T> inline void mpi_bcast_T(T *x, 
					const int &nwords, 
					const int &source=0){
  MPI_Bcast(x, nwords*sizeof(T), MPI_BYTE, source, MPI_COMM_WORLD);
}

template<class T> inline void mpi_gatherv_T(T xsend[], 
					    const int &sendcnt, 
					    T xrecv[],
					    int recvcnt[], 
					    int displacements[]){ 

  static int recv_cnt_buf[NPROC_MAX];
  static int disp_buf[NPROC_MAX+1];
  int t_size = sizeof(T);
  int Nproc = mpi_get_size();
  for(int i=0; i<Nproc; i++){
    recv_cnt_buf[i] = recvcnt[i] * t_size;
    disp_buf[i] = displacements[i] * t_size;
  }
  disp_buf[Nproc] = displacements[Nproc] * t_size;
  MPI_Gatherv(xsend, sendcnt*t_size, MPI_BYTE, xrecv, recv_cnt_buf, disp_buf, MPI_BYTE, 0, MPI_COMM_WORLD);
}


template<class T> inline void mpi_allgatherv_T(T val_send[],
					       const int &_cnt_send,
					       T val_recv[],
					       int &cnt_glb,
					       const int &_Nproc){
  int cnt_send = _cnt_send;
  int Nproc = _Nproc;
  static int len_array[NPROC_MAX];
  static int len_disp[NPROC_MAX+1];
  MPI_Allgather(&cnt_send, 1, MPI_INT, len_array, 1, MPI_INT, MPI_COMM_WORLD);
  len_disp[0] = 0;
  for(int i=0; i<Nproc; i++){
    len_disp[i+1] = len_disp[i] + len_array[i];
  }
  cnt_glb = len_disp[Nproc];
  int t_size=sizeof(T);
  for(int i=0; i<Nproc; i++){
    len_array[i] *= t_size;
    len_disp[i+1] *= t_size;
  }
  MPI_Allgatherv(val_send, cnt_send*t_size, MPI_BYTE, val_recv, len_array, len_disp, MPI_BYTE, MPI_COMM_WORLD);
}


template<class T> inline void mpi_allgatherv_T_2(T val_send[],
						 const int &_cnt_send,
						 T val_recv[],
						 int &cnt_glb,
						 const int &_Nproc,
						 int len_disp[]){
  int cnt_send = _cnt_send;
  int Nproc = _Nproc;
  static int len_array[NPROC_MAX];

  MPI_Allgather(&cnt_send, 1, MPI_INT, len_array, 1, MPI_INT, MPI_COMM_WORLD);
  len_disp[0] = 0;
  for(int i=0; i<Nproc; i++){
    len_disp[i+1] = len_disp[i] + len_array[i];
  }
  cnt_glb = len_disp[Nproc];
  static int len_disp_byte[NPROC_MAX+1];
  int t_size=sizeof(T);
  len_disp_byte[0] = 0;
  for(int i=0; i<Nproc; i++){
    len_array[i] *= t_size;
    len_disp_byte[i+1] = len_disp[i+1]*t_size;
  }
  MPI_Allgatherv(val_send, cnt_send*t_size, MPI_BYTE, val_recv, len_array, len_disp_byte, MPI_BYTE, MPI_COMM_WORLD);
}


inline double mpi_wtime(){
    //MPI_Barrier(MPI_COMM_WORLD);
    return MPI_Wtime();
}

#endif
