#include <iomanip>
#include <omp.h>

#include <types.h>
#include <my_errors.h>
#include <kernel.h>
#include <functions.h>
#include <utilis.h>

using namespace std;

extern "C"
HostError Hermite6th(const double TTIME, double* GTIME, double* ATIME, double* local_time, double* step, const unsigned int N, const unsigned int M, double4* pos_PH, float4* vel_PH, double4* pos_CH, double4* vel_CH, double4* a_H0, const unsigned int MAXDIM, unsigned int NGPU, unsigned int TPB, int rank, int size, unsigned int BFMAX, double ETA6, double ETA4, double DTMAX, double DTMIN, double EPS, double DTPRINT, unsigned int FMAX, const bool warm, double GTW, unsigned int GPUMINTHREADS, double plummer_core, double plummer_mass, double rscale, double mscale, vector<unsigned int> devices, bool *cleanstop, string path){


	unsigned int ompthreads = 1; // N must be integer multiple of this number
   omp_set_num_threads( ompthreads );
	unsigned int* vetint = new unsigned int [N];
	unsigned int *counter = new unsigned int [ompthreads];
	int *next = new int [N];
	unsigned long nextsize = N;
	double NEXTOUT = DTPRINT;

	double4 **pos_PD = new double4* [NGPU];
	float4  **vel_PD = new float4*  [NGPU];
	float4  **acc_PD = new float4*  [NGPU];
	double4 **pos_CD = new double4* [NGPU];
	double4 **vel_CD = new double4* [NGPU];
	int **next_D = new int* [NGPU];
	double  **loc_D  = new double*  [NGPU];

#ifdef APE
	double **step_D  = new double* [NGPU];
#endif
#ifdef GPUCORR
	double **step_D  = new double* [NGPU];
#endif

	double4 **a_D    = new double4* [NGPU];
   double4 **a1_D   = new double4* [NGPU];
   double4 **a2_D   = new double4* [NGPU];

   double4 **a_tot_D = new double4* [NGPU];
   double4 **a1_tot_D = new double4* [NGPU];
   double4 **a2_tot_D = new double4* [NGPU];


	double4 **a3_D    = new double4* [NGPU];
	double4 **p_v_a3_Dev = new double4* [NGPU];
	double4 **a_temp_Dev = new double4* [NGPU];


	unsigned int malloc_size  = MAXDIM*sizeof(double4); //it contains a, adot, a2dots sequentially
   unsigned int malloc_db4   = nextsize*sizeof(double4);
   unsigned int malloc_fl4   = nextsize*sizeof(float4);
   unsigned int malloc_ui    = nextsize*sizeof(unsigned int);
   unsigned int malloc_db    = nextsize*sizeof(double);
   unsigned int malloc_db4_N = N*sizeof(double4);

   for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
#ifdef GPUCORR
		DeviceSafeCall(cudaHostAlloc((void**)&step_D[i], malloc_db, cudaHostAllocMapped));
#endif

		DeviceSafeCall(cudaMalloc((void **)&a_D[i],    malloc_size));
      DeviceSafeCall(cudaMalloc((void **)&a1_D[i],    malloc_size));
      DeviceSafeCall(cudaMalloc((void **)&a2_D[i],    malloc_size));
      DeviceSafeCall(cudaMalloc((void **)&a_tot_D[i], malloc_db4_N));
      DeviceSafeCall(cudaMalloc((void **)&a1_tot_D[i], malloc_db4_N));
      DeviceSafeCall(cudaMalloc((void **)&a2_tot_D[i], malloc_db4_N));

      DeviceSafeCall(cudaMalloc((void **)&pos_PD[i],  malloc_db4));
      DeviceSafeCall(cudaMalloc((void **)&pos_CD[i],  malloc_db4));
      DeviceSafeCall(cudaMalloc((void **)&vel_CD[i],  malloc_db4));
      DeviceSafeCall(cudaMalloc((void **)&vel_PD[i],  malloc_fl4));
      DeviceSafeCall(cudaMalloc((void **)&acc_PD[i],  malloc_fl4));
      DeviceSafeCall(cudaMalloc((void **)&next_D[i],  malloc_ui));
      DeviceSafeCall(cudaMalloc((void **)&loc_D[i],   malloc_db));
#ifdef APE
		DeviceSafeCall(cudaMalloc((void **)&step_D[i],   malloc_db));
#endif
      DeviceSafeCall(cudaMalloc((void **)&a3_D[i],    malloc_db4_N));
		DeviceSafeCall(cudaMalloc((void **)&p_v_a3_Dev[i], 3*malloc_db4_N));
		DeviceSafeCall(cudaMalloc((void **)&a_temp_Dev[i], 3*malloc_db4_N));

      DeviceSafeCall(cudaMemcpy( pos_PD[i], pos_PH,     malloc_db4,    cudaMemcpyHostToDevice ));
      DeviceSafeCall(cudaMemcpy( vel_PD[i], vel_PH,     malloc_fl4,    cudaMemcpyHostToDevice ));
      DeviceSafeCall(cudaMemcpy( pos_CD[i], pos_CH,     malloc_db4,    cudaMemcpyHostToDevice ));
      DeviceSafeCall(cudaMemcpy( vel_CD[i], vel_CH,     malloc_db4,    cudaMemcpyHostToDevice ));

		DeviceSafeCall(cudaMemcpy(a_tot_D[i], a_H0,        malloc_db4_N,  cudaMemcpyHostToDevice));
		DeviceSafeCall(cudaMemcpy(a1_tot_D[i], &a_H0[N],   malloc_db4_N,  cudaMemcpyHostToDevice));
		DeviceSafeCall(cudaMemcpy(a2_tot_D[i], &a_H0[2*N], malloc_db4_N,  cudaMemcpyHostToDevice));

		DeviceSafeCall(cudaMemcpy(  loc_D[i], local_time, malloc_db,     cudaMemcpyHostToDevice));

#ifdef APE
		DeviceSafeCall(cudaMemcpy(  step_D[i], step, malloc_db,     cudaMemcpyHostToDevice));
#endif

#ifdef APE
		for(unsigned int j = 0; j < NGPU; j++){
			if(j != i)
				DeviceSafeCall(cudaDeviceEnablePeerAccess(devices[j], 0));
		}
#endif

   }

#ifdef GPUCORR
	for(unsigned int i = 0; i < NGPU; i++){
		DeviceSafeCall(cudaSetDevice(devices[i]));
		for(unsigned int j = 0; j < N; j++)
			step_D[i][j] = step[j];
	}
#endif

	for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      int BL = ceil((double)N/TPB);
      initvectors<<<BL, TPB>>>(a3_D[i], acc_PD[i]);
   }

	unsigned int ppG = N/(NGPU*size);
	unsigned int Bfactor;
	unsigned long BLOCKS;
	unsigned int DIMENSION, THREADS, bfmax, SHARED;
	double *mpi_red_aux = new double [3*N];
   double *mpi_red = new double [3*N];
	double4* a_H1 = new double4 [3*N];
	double4* p_v_a3 = new double4 [3*N];
	double4* a3_H = new double4 [N];
	double4 *a_H = new double4 [NGPU*3*N];

	for(unsigned int i = 0; i < 3*N; i++)
		a_H1[i] = a_H0[i];

	double E0, kk0, pp0;
	int out_index = 1;
	ofstream stream;
	string temp;
	char *output_name;

	if(warm){
	   while(NEXTOUT <= GTW) NEXTOUT += DTPRINT;
      while(NEXTOUT <= *GTIME) NEXTOUT += DTPRINT;
		if(rank == 0)
	      HostSafeCall(AcquireEnergy(&E0, path));
	}
	else{
	   HostSafeCall(Calculate_Energy(pos_CD, vel_CD, N, EPS, TPB, NGPU, rank, ppG, &kk0, &pp0, plummer_core, plummer_mass, rscale, mscale, devices));
		E0 = kk0 + pp0;
      temp = path + "energy.dat";
      output_name = to_char(temp);
      stream.open(output_name, ios::app);
      stream<<scientific<<setprecision(16);
      stream<<0.0<<"  "<<0.0<<"  "<<kk0<<"  "<<pp0<<"  "<<2.*kk0/fabs(pp0)<<endl;
      stream.close();
	}

	if(rank == 0){
		temp = path + "HiGPUslog.dat";
		output_name = to_char(temp);
		stream.open(output_name, ios::app);
		stream<<"==============================================="<<endl;
		stream<<scientific<<setprecision(16);
		stream<<"#Initial Total Energy : #"<<E0<<"#"<<endl;
		stream.close();

		string str = to_string(FMAX) + ".dat";
		stream.open(to_char(str), ios::out);
		for(unsigned int i = 0; i < M; i++)
			stream<<pos_CH[i].x<<"  "<<pos_CH[i].y<<"  "<<pos_CH[i].z<<"  "<<vel_CH[i].x<<"  "<<vel_CH[i].y<<"  "<<vel_CH[i].z<<"  "<<pos_CH[i].w<<endl;
		stream.close();

	}

	MPISafeCall(MPI_Barrier(MPI_COMM_WORLD));


	double start = 0.0;
	double end = 0.0;
	double start_program = 0.0;
	double end_program = 0.0;
	HiGPUsTimes *Times;
	Times = new HiGPUsTimes [N+1];
	for(unsigned int i = 0; i <= N; i++){
      Times[i].next_time = 0.0;
      Times[i].cpynext_time = 0.0;
      Times[i].predictor_time = 0.0;
      Times[i].evaluation_time = 0.0;
      Times[i].reduce_time = 0.0;
      Times[i].reposition_time = 0.0;
      Times[i].memcpy2_time = 0.0;
      Times[i].mpireduce_time = 0.0;
      Times[i].corrector_time = 0.0;
      Times[i].reconstruct_time = 0.0;
      Times[i].energy_time = 0.0;
		Times[i].rtime = 0.0;
		Times[i].thr = 0.0;
		Times[i].totthr = 0.0;
		Times[i].bfac = 0.0;
   }


//	HostSafeCall(GPU_memcheck(NGPU, devices, __FILE__, __LINE__));
	HostSafeCall(CPU_memcheck(__FILE__, __LINE__, path));


	struct timeval tv;
	gettimeofday(&tv, NULL);
	int sec = tv.tv_sec;
	int microsec = tv.tv_usec;
	start_program = sec + microsec * 0.000001;

	if(rank == 0)
		HostSafeCall(CheckHDDMemory(cleanstop, path));

	MPISafeCall(MPI_Barrier(MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast(cleanstop, 1, MPI::BOOL, 0, MPI_COMM_WORLD));

	if(*cleanstop)
		return HNoError;

#ifdef EXPANDING
	double plum0 = plummer_core;
#endif
	
	do{
		if(*GTIME >= TTIME){

			//ClearAll();
			gettimeofday(&tv, NULL);
			int sec = tv.tv_sec;
			int microsec = tv.tv_usec;
			end_program = sec + microsec * 0.000001;

			delete [] p_v_a3;
         delete [] a3_H;
         delete [] a_H;
         delete [] a_H1;
         delete [] mpi_red_aux;
         delete [] mpi_red;
         delete [] next;
         delete [] vetint;
         delete [] counter;
         delete [] Times;

			if(rank == 0){
			   temp = path + "HiGPUslog.dat";
            output_name = to_char(temp);
            stream.open(output_name, ios::app);
				stream<<scientific<<setprecision(6);
				stream<<" \n Total integration time : "<<end_program-start_program<<" seconds "<<endl;
				stream.close();
			}

			return HNoError;
		}

		get_times(&start);

#ifdef GPUCORR
			HostSafeCall(NextParticles(N, ompthreads, counter, vetint, *ATIME, local_time, step_D[0], next, &nextsize));
#else
        HostSafeCall(NextParticles(N, ompthreads, counter, vetint, *ATIME, local_time, step, next, &nextsize));
#endif 
         *GTIME = *ATIME;
#ifdef EXPANDING
			plummer_core = plum0*exp(GTIME);
#endif

			unsigned int dim2 = ceil((double)nextsize/TPB)*TPB;
			for(unsigned int i = nextsize; i < dim2; i++)
				next[i] = -1;
		get_times(&end);
		set_times(end-start, &(Times[nextsize].next_time));

		get_times(&start);
			for(unsigned int i = 0; i < NGPU; i++){
				DeviceSafeCall(cudaSetDevice(devices[i]));
				DeviceSafeCall(cudaMemcpy(next_D[i], next,  dim2 * sizeof( int ), cudaMemcpyHostToDevice));
			}
		get_times(&end);
		set_times(end-start, &(Times[nextsize].cpynext_time));

		get_times(&start);
			for(unsigned int i = 0; i < NGPU; i++){
		      DeviceSafeCall(cudaSetDevice(devices[i]));
			   int BL = ppG/TPB + ceil((double)nextsize/TPB);
				int istart = ppG*(i+rank*NGPU);
	
				Predictor <<<BL, TPB>>> (*GTIME, pos_PD[i], vel_PD[i], acc_PD[i], pos_CD[i], vel_CD[i], 
												loc_D[i], a_tot_D[i], a1_tot_D[i], a2_tot_D[i], a3_D[i], istart, 
												next_D[i], ppG, N);
			}
		get_times(&end);
		set_times(end-start, &(Times[nextsize].predictor_time));

		
 	   THREADS = TPB;
  		bfmax = BFMAX;
		*ATIME = 1.0e+10;

		Bfactor = 1;

		BLOCKS = ceil((double)nextsize/THREADS);

		while(THREADS*BLOCKS < GPUMINTHREADS && THREADS > 32){
			THREADS /= 2;
			bfmax *= 2;
			BLOCKS = ceil((double)nextsize/THREADS);
		}

		DIMENSION = THREADS*BLOCKS;

		while(THREADS*BLOCKS < GPUMINTHREADS && Bfactor < bfmax){
         BLOCKS *= 2;
         Bfactor *= 2;
      }

		SHARED = THREADS * (sizeof(double4) + 2 * sizeof(float4));

      for(unsigned int i = 0; i < NGPU; i++){
         DeviceSafeCall(cudaSetDevice(devices[i]));
         DeviceSafeCall(cudaThreadSynchronize());
		}

		set_times((double)Bfactor, &(Times[nextsize].bfac));
      set_times((double)THREADS, &(Times[nextsize].thr));
      set_times((double)BLOCKS*THREADS, &(Times[nextsize].totthr));


		get_times(&start);
			for(unsigned int i = 0; i < NGPU; i++){
		      DeviceSafeCall(cudaSetDevice(devices[i]));
	         DeviceCheckErrors();
		      int istart = ppG*(i+rank*NGPU);
		      evaluation<<< BLOCKS, THREADS, SHARED >>> ( N, pos_PD[i], vel_PD[i], acc_PD[i],  a_D[i], a1_D[i], a2_D[i],
	                                       istart, ppG, Bfactor, DIMENSION, next_D[i], loc_D[i], *GTIME, EPS, plummer_core, plummer_mass, rscale, mscale);
		   }
		get_times(&end);
		set_times(end-start, &(Times[nextsize].evaluation_time));
		
		int bl = BLOCKS;
	   int bf = Bfactor;
	   SHARED = THREADS * sizeof(double4);


		for(unsigned int i = 0; i < NGPU; i++){
		  DeviceSafeCall(cudaSetDevice(devices[i]));
		  DeviceSafeCall(cudaThreadSynchronize());
		  DeviceCheckErrors();
	   }


		get_times(&start);
			while(bf != 1){
		      bl>>=1;
				bf>>=1;
		      for(unsigned int i = 0; i < NGPU; i++){
		         DeviceSafeCall(cudaSetDevice(devices[i]));
		         reduce<<< 3*bl, THREADS, SHARED>>>(a_D[i], a1_D[i], a2_D[i], bf, DIMENSION);
		      }
		
		      for(unsigned int i = 0; i < NGPU; i++){
		         DeviceSafeCall(cudaSetDevice(devices[i]));
		         DeviceSafeCall(cudaThreadSynchronize());
		         DeviceCheckErrors();
		      }
		   }
		get_times(&end);
		set_times(end-start, &(Times[nextsize].reduce_time));

		get_times(&start);
			for(unsigned int i = 0; i < NGPU; i++){
       DeviceSafeCall(cudaSetDevice(devices[i]));
        reposition<<<DIMENSION/THREADS, THREADS>>>(a_D[i], a1_D[i], a2_D[i], a_temp_Dev[i], nextsize);
      }
		get_times(&end);
		set_times(end-start, &(Times[nextsize].reposition_time));

		unsigned int cpy_size = 3*nextsize;

#ifdef APE
		for(unsigned int i = 1; i < NGPU; i++){
			int SHRD = THREADS*sizeof(double4);
			DeviceSafeCall(cudaSetDevice(devices[0]));
			DeviceSafeCall(cudaMemcpyPeer(p_v_a3_Dev[0], devices[0], a_temp_Dev[i], devices[i], cpy_size*sizeof(double4)));
			sum_partial<<<3*DIMENSION/THREADS, THREADS, SHRD>>>(p_v_a3_Dev[0], a_temp_Dev[0], 3*nextsize);
			DeviceSafeCall(cudaThreadSynchronize());
		}



		// QUI VA AGGIUNTA LA FUNZIONE CHE RIDUCE A_TEMP_DEV SU TUTTE LE GPU 0 E LO DAI POI A TUTTE LE GPU DI TUTTI I NODI : not yet implemented



#else
		get_times(&start);
	   for(unsigned int i = 0; i < NGPU; i++){
	      DeviceSafeCall(cudaSetDevice(devices[i]));
	      DeviceCheckErrors();
			DeviceSafeCall(cudaMemcpy(&a_H[i*cpy_size], a_temp_Dev[i], cpy_size*sizeof(double4), cudaMemcpyDeviceToHost));
	   }
		get_times(&end);
		set_times(end-start, &(Times[nextsize].memcpy2_time));

		get_times(&start);
	   HostSafeCall(ReduceAll(cpy_size, N, NGPU, nextsize, a_H, a_H1, mpi_red_aux, mpi_red, next));
		get_times(&end);
		set_times(end-start, &(Times[nextsize].mpireduce_time));
#endif

#ifdef GPUCORR
		get_times(&start);
		for(unsigned int i = 0; i < NGPU; i++){
			DeviceSafeCall(cudaSetDevice(devices[i]));
			DeviceSafeCall(cudaMemcpy(a_temp_Dev[i], a_H, 3*nextsize*sizeof( double4 ), cudaMemcpyHostToDevice ));
		}

		for(unsigned int i = 0; i < NGPU; i++){
         DeviceSafeCall(cudaSetDevice(devices[i]));
         Corrector_gpu<<<DIMENSION/THREADS, THREADS>>>(*GTIME, loc_D[i], step_D[i], next_D[i], nextsize, pos_CD[i], vel_CD[i], 
					a_tot_D[i], a1_tot_D[i], a2_tot_D[i], a_temp_Dev[i], a3_D[i], ETA6, ETA4, DTMAX, DTMIN, N);// chiama direttamete corrector
      }

		DeviceSafeCall(cudaSetDevice(devices[0]));
		DeviceSafeCall(cudaThreadSynchronize());

		for(unsigned int i = 0; i < nextsize; i++){
         int who = next[i];
			local_time[who] = *GTIME;
         *ATIME = min (local_time[who] + step_D[0][who], *ATIME);
      }

		get_times(&end);
      set_times(end-start, &(Times[nextsize].corrector_time));

#elif APE

	for(unsigned int i = 0; i < NGPU; i++){
         DeviceSafeCall(cudaSetDevice(devices[i]));
         Corrector_gpu<<<DIMENSION/THREADS, THREADS>>>(*GTIME, loc_D[i], step_D[i], next_D[i], nextsize, pos_CD[i], vel_CD[i],
               a_tot_D[i], a1_tot_D[i], a2_tot_D[i], a_temp_Dev[i], a3_D[i], ETA6, ETA4, DTMAX, DTMIN, N);// chiama direttamete corrector
      }

      DeviceSafeCall(cudaSetDevice(devices[0]));
      DeviceSafeCall(cudaThreadSynchronize());
      DeviceSafeCall(cudaMemcpy(step, step_D[0], N*sizeof(double), cudaMemcpyDeviceToHost));

      for(unsigned int i = 0; i < nextsize; i++){
         int who = next[i];
         local_time[who] = *GTIME;
         *ATIME = min (local_time[who] + step[who], *ATIME);
      }

#else
		// corrector su cpu
		for(unsigned int i = 0; i < NGPU; i++){
			DeviceSafeCall(cudaSetDevice(devices[i]));
			update_local_time<<<DIMENSION/THREADS, THREADS>>>(next_D[i], loc_D[i], *GTIME);
		}
		get_times(&start);
		HostSafeCall(Corrector(GTIME, ATIME, local_time, step, next, nextsize, pos_CH, vel_CH, a_H0,
         a_H1, a3_H, p_v_a3, ETA6, ETA4, DTMAX, DTMIN, N));
		get_times(&end);
		set_times(end-start, &(Times[nextsize].corrector_time));
#endif

#ifndef APE
#ifndef GPUCORR
		get_times(&start);
		for(unsigned int i = 0; i < NGPU; i++){
				DeviceSafeCall(cudaSetDevice(devices[i]));
				DeviceSafeCall(cudaMemcpy(p_v_a3_Dev[i], p_v_a3, 3*nextsize*sizeof( double4 ), cudaMemcpyHostToDevice ));
				DeviceSafeCall(cudaMemcpy(a_temp_Dev[i], a_H,    3*nextsize*sizeof( double4 ), cudaMemcpyHostToDevice ));
		}
		get_times(&end);
		set_times(end-start, &(Times[nextsize].rtime));

		for(unsigned int i = 0; i < NGPU; i++){
			DeviceSafeCall(cudaSetDevice(devices[i]));
			int BB = 6*DIMENSION/THREADS;
				Reconstruct<<< BB, THREADS >>>(next_D[i], nextsize, pos_CD[i], vel_CD[i], a3_D[i], a_tot_D[i], a1_tot_D[i], a2_tot_D[i], p_v_a3_Dev[i], a_temp_Dev[i]);
		}
		get_times(&end);
		set_times(end-start, &(Times[nextsize].reconstruct_time));
#endif
#endif
		if((*GTIME+GTW) >= NEXTOUT ){

			DeviceSafeCall(cudaSetDevice(devices[0]));
			DeviceSafeCall(cudaMemcpy( pos_CH, pos_CD[0],     malloc_db4,    cudaMemcpyDeviceToHost ));
			DeviceSafeCall(cudaMemcpy( vel_CH, vel_CD[0],     malloc_db4,    cudaMemcpyDeviceToHost ));
			CheckBlocks(step, M, path);
			double kk,pp;
			get_times(&start);
				HostSafeCall(Calculate_Energy(pos_CD, vel_CD, N, EPS, TPB, NGPU, rank, ppG, &kk, &pp, plummer_core, plummer_mass, rscale, mscale, devices));
			get_times(&end);
			set_times(end-start, &(Times[nextsize].energy_time));

			if(rank == 0){
				HostSafeCall(CheckHDDMemory(cleanstop, path));

#ifdef CHECK_TIMES
            string ffff = to_string(out_index + FMAX);
            ffff = path + "times_"+ffff+".dat";
            
            stream.open(to_char(ffff), ios::out);

            stream<<scientific<<setprecision(6);
            stream<<"N  "<<"    NEXT  "<<"      CPY_NEXT"<<"        PRED "<<"       EVAL "<<"           REDU "<<"        REPOS "<<"       CPY_ACC "<<"        MPI "<<"         CORR   "<<"     CPY_REC "<<"        RECON "<<"     THREADS "<<"    TOTTHREAD "<<"      BFACT "<<endl;
            for(unsigned int i = 1; i <= N; i++){
               if(Times[i].next_time != 0.0)
                  stream<<i<<"  "<<
                     Times[i].next_time<<"  "<<
                     Times[i].cpynext_time<<"  "<<
                     Times[i].predictor_time<<"  "<<
                     Times[i].evaluation_time<<"  "<<
                     Times[i].reduce_time<<"  "<<
                     Times[i].reposition_time<<"  "<<
                     Times[i].memcpy2_time<<"  "<<
                     Times[i].mpireduce_time<<"  "<<
                     Times[i].corrector_time<<"  "<<
                     Times[i].rtime<<"  "<<
                     Times[i].reconstruct_time<<"  "<<
                     Times[i].thr<<"  "<<
                     Times[i].totthr<<"  "<<
                     Times[i].bfac<<endl;

					Times[i].next_time = 0.0;

            }
            stream.close();
#endif

				double E = kk + pp;
				temp = path + "energy.dat";
            output_name = to_char(temp);
            stream.open(output_name, ios::app);
            stream<<scientific<<setprecision(16);
            stream<<*GTIME+GTW<<"  "<<fabs((E-E0)/E0)<<"  "<<kk<<"  "<<pp<<"  "<<2.*kk/fabs(pp)<<endl;
            stream.close();

				string file_name = path + to_string(out_index + FMAX);
				file_name += ".dat";
				stream.open(to_char(file_name), ios::out);
				stream<<scientific<<setprecision(16);
				for(unsigned int i = 0; i < M; i++)
					stream<<pos_CH[i].x<<"  "<<pos_CH[i].y<<"  "<<pos_CH[i].z<<"  "<<vel_CH[i].x<<"  "<<vel_CH[i].y<<"  "<<vel_CH[i].z<<"  "<<pos_CH[i].w<<endl;
				stream.close();
				out_index++;
			}

			MPISafeCall(MPI_Barrier(MPI_COMM_WORLD));
		   MPISafeCall(MPI_Bcast(cleanstop, 1, MPI::BOOL, 0, MPI_COMM_WORLD));

			if(*cleanstop)
				return HNoError;

			NEXTOUT+=DTPRINT;
		}

	}while(1);

}
