#include <functions.h>
#include <my_errors.h>
#include <kernel.h>

HostError Calculate_Energy(double4 **pos_CD, double4 **vel_CD, unsigned int N, double EPS, unsigned int THREADS, unsigned int NGPU, int rank, unsigned int ppG, double *kinetic, double *potential, double plummer_core, double plummer_mass, double rscale, double mscale, vector<unsigned int> devices){

   double *KH = new double [N*NGPU];
   double *PH = new double [N*NGPU];
   double *Kmpi = new double [N];
   double *Pmpi = new double [N];

   double *K [NGPU];
   double *P [NGPU];

	for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      DeviceSafeCall(cudaMalloc(( void**)&K[i], N*sizeof(double)));
		DeviceSafeCall(cudaMalloc(( void**)&P[i], N*sizeof(double)));
   }

   *kinetic = 0.0;
   *potential = 0.0;

   for(unsigned int i = 0; i < N; i++)
      Kmpi[i] = Pmpi[i] = 0.0;


   int BLOCKS = N/THREADS;
	int SHARED = THREADS*sizeof(double4);

	for(unsigned int i = 0; i < NGPU; i++)
		DeviceSafeCall(cudaThreadSynchronize());

	for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      int istart = ppG*(i+rank*NGPU);
      energy<<<BLOCKS, THREADS, SHARED>>>(pos_CD[i], vel_CD[i], K[i], P[i], N, EPS, istart, ppG, plummer_core, plummer_mass, rscale, mscale);
   }

	for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      DeviceSafeCall(cudaMemcpy(&KH[i*N], K[i], N*sizeof(double), cudaMemcpyDeviceToHost));
		DeviceSafeCall(cudaMemcpy(&PH[i*N], P[i], N*sizeof(double), cudaMemcpyDeviceToHost));
   }

   for(unsigned int i = 0; i < N; i++){
      for(unsigned int dev = 1; dev < NGPU; dev++){
         unsigned int p = i + dev*N;
         KH[i] += KH[p];
         PH[i] += PH[p];
      }
   }

   MPISafeCall(MPI_Allreduce(KH, Kmpi, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
   MPISafeCall(MPI_Allreduce(PH, Pmpi, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));

   for(unsigned int i = 0; i < N; i++){
      *kinetic += Kmpi[i];
      *potential += Pmpi[i];
   }

   delete [] Kmpi;
   delete [] KH;
   delete [] Pmpi;
   delete [] PH;

	for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      DeviceSafeCall(cudaFree(P[i]));
		DeviceSafeCall(cudaFree(K[i]));
   }


	   return HNoError;

}


