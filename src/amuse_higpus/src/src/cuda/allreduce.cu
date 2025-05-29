#include <mpi.h>
#include <my_errors.h>
#include <functions.h>

using namespace std;

HostError ReduceAll(unsigned int cpy_size, unsigned int N, unsigned int NGPU, unsigned long nextsize, double4* a_H, double4 *a_H0, double *mpi_red_aux, double *mpi_red, int *next){

	  for(unsigned int i = 0; i < cpy_size; i++){
      for(unsigned int dev = 1; dev < NGPU; dev++){
         unsigned int p = i + dev*cpy_size;
         a_H[i].x += a_H[p].x;
         a_H[i].y += a_H[p].y;
         a_H[i].z += a_H[p].z;
      }
   }

      for(unsigned int i = 0; i < cpy_size; i++){
         mpi_red_aux[i] = a_H[i].x;
         mpi_red[i] = 0.0;
      }

   MPISafeCall(MPI_Allreduce(mpi_red_aux, mpi_red, cpy_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));

   for(unsigned int k = 0; k < 3; k++){
      unsigned int p = k*nextsize;
      for(unsigned int i = 0; i < nextsize; i++){
         int who = next[i] + k*N;
         a_H0[who].x = a_H[i+p].x = mpi_red[i+p];
      }
   }
      for(unsigned int i = 0; i < cpy_size; i++){
         mpi_red_aux[i] = a_H[i].y;
         mpi_red[i] = 0.0;
      }

   MPISafeCall(MPI_Allreduce(mpi_red_aux, mpi_red, cpy_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
     
	for(unsigned int k = 0; k < 3; k++){
      unsigned int p = k*nextsize;
      for(unsigned int i = 0; i < nextsize; i++){
         int who = next[i] + k*N;
         a_H0[who].y = a_H[i+p].y = mpi_red[i+p];
      }
   }

      for(unsigned int i = 0; i < cpy_size; i++){
         mpi_red_aux[i] = a_H[i].z;
         mpi_red[i] = 0.0;
      }

   MPISafeCall(MPI_Allreduce(mpi_red_aux, mpi_red, cpy_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));

   for(unsigned int k = 0; k < 3; k++){
      unsigned int p = k*nextsize;
      for(unsigned int i = 0; i < nextsize; i++){
         int who = next[i] + k*N;
         a_H0[who].z = a_H[i+p].z = mpi_red[i+p];
      }
   }


   return HNoError;

}

