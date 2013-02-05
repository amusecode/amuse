#include <functions.h>
#include <cmath>
#include <my_errors.h>
#include <fstream>
#include <kernel.h>
#include <utilis.h>
#
HostError CheckBlocks(double *step, unsigned int M, string path){

   int local_rank;
	string temp;
	char *output_name;
   
	MPISafeCall(MPI_Comm_rank(MPI_COMM_WORLD, &local_rank));

   if(local_rank == 0){
      ofstream blocks;
      temp = path + "Blocks.dat";
      output_name = to_char(temp);
      blocks.open(output_name);		
      
		int *bl = new int [35];
      double conto;

      for(int i = 0; i < 30; i++)
         bl[i]=0;

      for(unsigned int i = 0; i < M; i++){
         conto = -1.00001*log(step[i])/log(2.0);
         if(conto>30)
            return HInvalidBlocks;
         bl[(int)conto]++;
      }

      for(int i = 0; i < 30; i++)
         blocks<<-i<<"  "<<bl[i]<<endl;

      blocks.close();
      delete [] bl;
   }

   return HNoError;

}


HostError DetermineSteps(double stp, unsigned int N, unsigned int M, double4* a_H0, double4* a_H1, double ETA4, double DTMIN, double DTMAX, double *step, double *ACTUAL_TIME, double *local_time){

   for(unsigned int who = 0; who < M; who++){
      double overh3 = 1.0/(stp*stp*stp);
      double overh2 = 1.0/(stp*stp);
      unsigned int who1 = who + N;

    double a2dx = overh2 * (-6.0 * (a_H0[who].x - a_H1[who].x) -
             stp * (4.0 * a_H1[who1].x + 2.0 * a_H0[who1].x));
    double a2dy = overh2 * (-6.0 * (a_H0[who].y - a_H1[who].y) -
             stp * (4.0 * a_H1[who1].y + 2.0 * a_H0[who1].y));
    double a2dz = overh2 * (-6.0 * (a_H0[who].z - a_H1[who].z) -
             stp * (4.0 * a_H1[who1].z + 2.0 * a_H0[who1].z));

    double a3dx = overh3 * (12.0 * (a_H0[who].x - a_H1[who].x) +
             6.0 * stp * (a_H1[who1].x + a_H0[who1].x));
    double a3dy = overh3 * (12.0 * (a_H0[who].y - a_H1[who].y) +
             6.0 * stp * (a_H1[who1].y + a_H0[who1].y));
    double a3dz = overh3 * (12.0 * (a_H0[who].z - a_H1[who].z) +
             6.0 * stp * (a_H1[who1].z + a_H0[who1].z));

    double a2dotsmod =  a2dx*a2dx + a2dy*a2dy + a2dz*a2dz;

    double a3dotsmod = a3dx*a3dx + a3dy*a3dy + a3dz*a3dz;

    double amod = a_H1[who].x * a_H1[who].x +
      a_H1[who].y * a_H1[who].y +
      a_H1[who].z * a_H1[who].z ;

    double adotmod = a_H1[who1].x * a_H1[who1].x +
      a_H1[who1].y * a_H1[who1].y +
      a_H1[who1].z * a_H1[who1].z ;

    double dt = sqrt(ETA4*(sqrt(amod*a2dotsmod) + adotmod) / (sqrt(adotmod*a3dotsmod) + a2dotsmod));

    if(dt>DTMAX)
       dt = DTMAX;

    if(dt<DTMIN)
       dt = DTMIN;

	 int exponent = log(dt)/log(2.0) - 1;

    step[who] = pow(2.0,exponent);

    *ACTUAL_TIME = min(step[who],*ACTUAL_TIME);
    local_time[who] = 0.0;
   }
 
   for(unsigned int who = M; who < N; who++){
      step[who] = 1.0e+10;
      local_time[who] = 1.0e+10;
   }

   return HNoError;
}

extern "C"
HostError InitBlocks(double4 *pos_PH, float4 *vel_PH, unsigned int TPB, unsigned int N, unsigned int M, unsigned int BFMAX, double ETA4, double DTMIN, double DTMAX, unsigned int NGPU, double EPS, unsigned int *MAXDIM, unsigned int *GPUMINTHREADS, string gpu_name, int rank, int nodes, double4* pos_CH, double4* vel_CH, double4* a_H0, double* step, double* local_time, double* ACTUAL_TIME, const bool vir, const double ratio, const bool warm, const bool setdev, vector<unsigned int>& devices, double plummer_core, double plummer_mass, double rscale, double mscale, string path){

   HostSafeCall(CudaInit(GPUMINTHREADS, NGPU, rank, gpu_name, setdev, devices, path));
   
	HostSafeCall(Max_dimension(TPB, BFMAX, N, MAXDIM, *GPUMINTHREADS));
cout<<N<<"  "<<*MAXDIM<<endl;
    double4 **a_D    = new double4* [NGPU];
   double4 **a1_D   = new double4* [NGPU];
   double4 **a2_D   = new double4* [NGPU];

   double4 **a_tot_D = new double4* [NGPU];
   double4 **a1_tot_D = new double4* [NGPU];
   double4 **a2_tot_D = new double4* [NGPU];

   double4 **pos_PD = new double4* [NGPU];
   double  **loc_D  = new double*  [NGPU];
   float4  **vel_PD = new float4*  [NGPU];
   float4  **acc_PD = new float4*  [NGPU];
   double4 **pos_CD = new double4* [NGPU];
   double4 **vel_CD = new double4* [NGPU];
   double4 **a3_D    = new double4* [NGPU];

   double4 **a_temp_Dev = new double4* [NGPU];

   double4 *a_H = new double4 [NGPU*3*N];
   double4 *a_H1 = new double4 [3*N];
 
   string temp;
	char *output_name;

	int **next_D = new int* [NGPU];

   for(unsigned int i = 0; i < N; i++)
      pos_CH[i] = pos_PH[i];

   int * next = new int [N];
   for(unsigned int i = 0; i < N; i++)
      next[i] = i;
   unsigned long nextsize = N;

   double *mpi_red_aux = new double [3*N];
   double *mpi_red = new double [3*N];
   double stp = 0.001;

   HostSafeCall(CPU_memcheck(__FILE__, __LINE__, path));


   unsigned long BL = ceil((double)nextsize/TPB);
   int dim = TPB*BL;
   unsigned int threads = TPB;
   unsigned int bfmax = BFMAX;

   unsigned int Bfactor = 1;

   while(threads*BL < *GPUMINTHREADS && threads > 32){
      threads /= 2;
      bfmax *= 2;
      BL = ceil((double)nextsize/threads);
   }
   dim = threads * BL;


   while(threads*BL < *GPUMINTHREADS && Bfactor < bfmax){
      BL *= 2;
      Bfactor *= 2;
   }

   unsigned int malloc_size  = (*MAXDIM)*sizeof(double4); //it contains a, adot, a2dots sequentially
   unsigned int malloc_db4   = nextsize*sizeof(double4);
   unsigned int malloc_fl4   = nextsize*sizeof(float4);
   unsigned int malloc_ui    = nextsize*sizeof(unsigned int);
   unsigned int malloc_db    = nextsize*sizeof(double);
   unsigned int malloc_db4_N = N*sizeof(double4);

   for(unsigned int i = 0; i < NGPU; i++){
 		DeviceSafeCall(cudaSetDevice(devices[i]));

      DeviceSafeCall(cudaMalloc((void **)&a_D[i],    malloc_size));
      DeviceSafeCall(cudaMalloc((void **)&a1_D[i],    malloc_size));
      DeviceSafeCall(cudaMalloc((void **)&a2_D[i],    malloc_size));
      DeviceSafeCall(cudaMalloc((void **)&a_tot_D[i], malloc_db4_N));
      DeviceSafeCall(cudaMalloc((void **)&a1_tot_D[i], malloc_db4_N));
      DeviceSafeCall(cudaMalloc((void **)&a2_tot_D[i], malloc_db4_N));

      DeviceSafeCall(cudaMalloc((void **)&pos_PD[i], malloc_db4));
      DeviceSafeCall(cudaMalloc((void **)&pos_CD[i], malloc_db4));
      DeviceSafeCall(cudaMalloc((void **)&vel_CD[i], malloc_db4));
      DeviceSafeCall(cudaMalloc((void **)&vel_PD[i], malloc_fl4));
      DeviceSafeCall(cudaMalloc((void **)&acc_PD[i], malloc_fl4));
      DeviceSafeCall(cudaMalloc((void **)&next_D[i], malloc_ui));
      DeviceSafeCall(cudaMalloc((void **)&loc_D[i],  malloc_db));
      DeviceSafeCall(cudaMalloc((void **)&a3_D[i],   malloc_db4_N));
      DeviceSafeCall(cudaMalloc((void **)&a_temp_Dev[i], 3*malloc_db4_N));

      DeviceSafeCall(cudaMemcpy( pos_PD[i], pos_PH, malloc_db4, cudaMemcpyHostToDevice ));
      DeviceSafeCall(cudaMemcpy( vel_PD[i], vel_PH, malloc_fl4, cudaMemcpyHostToDevice ));
      DeviceSafeCall(cudaMemcpy( pos_CD[i], pos_CH, malloc_db4, cudaMemcpyHostToDevice ));
      DeviceSafeCall(cudaMemcpy( vel_CD[i], vel_CH, malloc_db4, cudaMemcpyHostToDevice ));
  }

   int ppG = N/(NGPU*nodes);

   if(vir && (!warm)){
      double kk, pp, virial_old;
      HostSafeCall(Calculate_Energy(pos_CD, vel_CD, N, EPS, TPB, NGPU, rank, ppG, &kk, &pp, plummer_core, plummer_mass, rscale, mscale, devices));
      virial_old = 2. * kk / fabs(pp);

      ofstream hlog;
		temp = path + "HiGPUslog.dat";
		output_name = to_char(temp);
		hlog.open(output_name, ios::app);
      hlog<<"==============================================="<<endl;
      hlog<<" Old virial ratio    : "<<virial_old<<endl;
      hlog.close();

      double scale = sqrt(ratio / virial_old);
      for(unsigned int i = 0; i < N; i++){
          vel_CH[i].x *= scale;
          vel_CH[i].y *= scale;
          vel_CH[i].z *= scale;
          vel_PH[i].x = vel_CH[i].x;
          vel_PH[i].y = vel_CH[i].y;
          vel_PH[i].z = vel_CH[i].z;
      }
   }

	for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
		DeviceSafeCall(cudaMemcpy( vel_PD[i], vel_PH, malloc_fl4, cudaMemcpyHostToDevice ));
		DeviceSafeCall(cudaMemcpy( vel_CD[i], vel_CH, malloc_db4, cudaMemcpyHostToDevice ));
	}


	for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      int BLT = ceil((double)N/threads);
      initvectors<<<BLT, threads>>>(a3_D[i], acc_PD[i]);
   }

   unsigned int dim2 = ceil((double)nextsize/TPB)*TPB;

   for(unsigned int i = nextsize; i < dim2; i++)
      next[i] = -1;

	for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      DeviceSafeCall(cudaMemcpy( next_D[i], next, malloc_ui, cudaMemcpyHostToDevice ));
   }

	int SHR = threads * (sizeof(double4) + 2 * sizeof(float4));

	for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      int istart = ppG*(i+rank*NGPU);
      evaluation<<< BL, threads, SHR >>> ( N, pos_PD[i], vel_PD[i], acc_PD[i],  a_D[i], a1_D[i], a2_D[i],
                                       istart, ppG, Bfactor, dim, next_D[i], loc_D[i], 0.0, EPS, plummer_core, plummer_mass, rscale, mscale);
   }

	for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
		update_local_time<<<dim/threads, threads>>>(next_D[i], loc_D[i], 0.0);
	}

   for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      DeviceSafeCall(cudaThreadSynchronize());
      DeviceCheckErrors();
   }


   int bl = BL;
   int bf = Bfactor;
	SHR = threads * sizeof(double4);

	while(bf != 1){
      bl>>=1;
		bf>>=1;
      for(unsigned int i = 0; i < NGPU; i++){
         DeviceSafeCall(cudaSetDevice(devices[i]));
         reduce<<< 3*bl, threads, SHR>>>(a_D[i], a1_D[i], a2_D[i], bf, dim);
      }

      for(unsigned int i = 0; i < NGPU; i++){
         DeviceSafeCall(cudaSetDevice(devices[i]));
         DeviceSafeCall(cudaThreadSynchronize());
         DeviceCheckErrors();
      }
   }

	for(unsigned int i = 0; i < NGPU; i++){
         DeviceSafeCall(cudaSetDevice(devices[i]));
         DeviceSafeCall(cudaThreadSynchronize());
	}
	
	for(unsigned int i = 0; i < NGPU; i++){
		 DeviceSafeCall(cudaSetDevice(devices[i]));
        reposition<<<bl, threads>>>(a_D[i], a1_D[i], a2_D[i], a_temp_Dev[i], nextsize);
      }



   unsigned int cpy_size = 3*nextsize;

	for(unsigned int i = 0; i < NGPU; i++){
         DeviceSafeCall(cudaSetDevice(devices[i]));
         DeviceSafeCall(cudaThreadSynchronize());
   }

	for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      DeviceCheckErrors();
      DeviceSafeCall(cudaMemcpy(&a_H[i*cpy_size], a_temp_Dev[i], cpy_size*sizeof(double4), cudaMemcpyDeviceToHost));
   }


   HostSafeCall(ReduceAll(cpy_size, N, NGPU, nextsize, a_H, a_H0, mpi_red_aux, mpi_red, next));

   for(unsigned int i = 2*N; i < 3*N; i++) //it puts the snap to 0.0
      a_H0[i].x = a_H0[i].y = a_H0[i].z = a_H0[i].w = 0.0;

	for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      DeviceSafeCall(cudaMemcpy(a_tot_D[i], a_H0, N*sizeof(double4), cudaMemcpyHostToDevice));
      DeviceSafeCall(cudaMemcpy(a1_tot_D[i], &a_H0[N], N*sizeof(double4), cudaMemcpyHostToDevice));
      DeviceSafeCall(cudaMemcpy(a2_tot_D[i], &a_H0[2*N], N*sizeof(double4), cudaMemcpyHostToDevice));
	}

	for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      int BLT = ppG/TPB + ceil((double)nextsize/TPB);
      int istart = ppG*(i+rank*NGPU);
      Predictor <<<BLT, TPB>>> (stp, pos_PD[i], vel_PD[i], acc_PD[i], pos_CD[i], vel_CD[i], loc_D[i], a_tot_D[i], a1_tot_D[i], a2_tot_D[i], a3_D[i], istart, next_D[i], ppG, N);
    }
    for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      DeviceCheckErrors();
   }


	 SHR = threads * (sizeof(double4) + 2 * sizeof(float4));


   for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      int istart = ppG*(i+rank*NGPU);
      evaluation<<< BL, threads, SHR >>> ( N, pos_PD[i], vel_PD[i], acc_PD[i],  a_D[i], a1_D[i], a2_D[i],
                                       istart, ppG, Bfactor, dim, next_D[i], loc_D[i], 0.0, EPS, plummer_core, plummer_mass, rscale, mscale);
   }

   for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      DeviceSafeCall(cudaThreadSynchronize());
      DeviceCheckErrors();
   }


   bl = BL;
   bf = Bfactor;
   SHR = threads * sizeof(double4);

   while(bf != 1){
      bl>>=1;
		bf>>=1;
      for(unsigned int i = 0; i < NGPU; i++){
         DeviceSafeCall(cudaSetDevice(devices[i]));
         reduce<<< 3*bl, threads, SHR>>>(a_D[i], a1_D[i], a2_D[i], bf, dim);
      }

      for(unsigned int i = 0; i < NGPU; i++){
         DeviceSafeCall(cudaSetDevice(devices[i]));
         DeviceSafeCall(cudaThreadSynchronize());
         DeviceCheckErrors();
      }
   }

   for(unsigned int i = 0; i < NGPU; i++){
         DeviceSafeCall(cudaSetDevice(devices[i]));
         DeviceSafeCall(cudaThreadSynchronize());
   }

   for(unsigned int i = 0; i < NGPU; i++){
       DeviceSafeCall(cudaSetDevice(devices[i]));
        reposition<<<bl, threads>>>(a_D[i], a1_D[i], a2_D[i], a_temp_Dev[i], nextsize);
      }

   cpy_size = 3*nextsize;

	   for(unsigned int i = 0; i < NGPU; i++){
         DeviceSafeCall(cudaSetDevice(devices[i]));
         DeviceSafeCall(cudaThreadSynchronize());
   }

   for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      DeviceCheckErrors();
      DeviceSafeCall(cudaMemcpy(&a_H[i*cpy_size], a_temp_Dev[i], cpy_size*sizeof(double4), cudaMemcpyDeviceToHost));
   }


   HostSafeCall(ReduceAll(cpy_size, N, NGPU, nextsize,  a_H, a_H1, mpi_red_aux, mpi_red, next));
   HostSafeCall(DetermineSteps(stp, N, M, a_H0, a_H1, ETA4, DTMIN, DTMAX, step, ACTUAL_TIME, local_time));
   HostSafeCall(CheckBlocks(step, M, path));

	for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      DeviceSafeCall(cudaFree(a_D[i]));
      DeviceSafeCall(cudaFree(a1_D[i]));
      DeviceSafeCall(cudaFree(a2_D[i]));

      DeviceSafeCall(cudaFree(a_tot_D[i]));
      DeviceSafeCall(cudaFree(a1_tot_D[i]));
      DeviceSafeCall(cudaFree(a2_tot_D[i]));

      DeviceSafeCall(cudaFree(pos_PD[i]));
      DeviceSafeCall(cudaFree(pos_CD[i]));
      DeviceSafeCall(cudaFree(vel_CD[i]));
      DeviceSafeCall(cudaFree(vel_PD[i]));
      DeviceSafeCall(cudaFree(acc_PD[i]));
      DeviceSafeCall(cudaFree(next_D[i]));
      DeviceSafeCall(cudaFree(loc_D[i]));
      DeviceSafeCall(cudaFree(a3_D[i]));
   }
   
	delete [] a_H;
   delete [] a_H1;
   delete [] mpi_red_aux;
   delete [] mpi_red;

   delete [] next;

   HostSafeCall(CPU_memcheck(__FILE__, __LINE__,path));

   return HNoError;

}
