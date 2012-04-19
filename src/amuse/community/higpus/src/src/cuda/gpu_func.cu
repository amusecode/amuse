#include <iomanip>
#include <omp.h>

#include <types.h>
#include <my_errors.h>
#include <kernel.h>
#include <functions.h>
#include <utilis.h>

using namespace std;

HostError NextParticles(unsigned int N, unsigned int ompthreads, unsigned int* counter, unsigned int* vetint, double* ATIME, double *local_time, double* step, int* next, unsigned long* nextsize){

	*nextsize = 0;

#pragma omp parallel
	{
		unsigned int cpu_thread_id = omp_get_thread_num();
		unsigned int istart = cpu_thread_id*N/ompthreads;
		counter[cpu_thread_id] = 0;

#pragma omp for
		for(unsigned int i = 0; i < N; i++){ 
			if( (local_time[i] + step[i]) == *ATIME){
				vetint[counter[cpu_thread_id] + istart] = i;
				counter[cpu_thread_id]++;
			} 
		}

#pragma omp barrier 
		unsigned int from = 0;
      unsigned int to = 0;

	   for(unsigned int i = cpu_thread_id; i > 0; i--)
			from += counter[i-1];

	   to = from + counter[cpu_thread_id];

	   for(unsigned int i = from; i < to; i++)
			next[i] = vetint[istart + i - from];
	}

	for(unsigned int i = 0; i < ompthreads; i++)
		*nextsize += counter[i];

	return HNoError;
}

HostError Corrector(double *GTIME, double *ATIME, double *local_time, double *step, int *next, unsigned long nextsize, double4 *pos_CH, double4 *vel_CH, double4 *a_H1, double4 *a_H0, double4 *a3_H, double4 *p_v_a3, double ETA6, double ETA4, double DTMAX, double DTMIN, unsigned int N){

	for(unsigned int i = 0; i < nextsize; i++){

		double dt;
		int who = next[i];
		int who1 = who + N;
		int who2 = who1 + N;

		double amod = a_H0[who].x * a_H0[who].x +
           			  a_H0[who].y * a_H0[who].y +
           			  a_H0[who].z * a_H0[who].z ;

		double adotmod = a_H0[who1].x * a_H0[who1].x +
							  a_H0[who1].y * a_H0[who1].y +
							  a_H0[who1].z * a_H0[who1].z ;

		double a2dotsmod =  a_H0[who2].x * a_H0[who2].x +
								  a_H0[who2].y * a_H0[who2].y +
								  a_H0[who2].z * a_H0[who2].z ;

	   double h = *GTIME-local_time[who];
		local_time[who] = *GTIME;

		double h1 = 0.5*h;
		double h2 = h*h*0.1;
		double h3 = h*h*h / 120.0;


		pos_CH[who].x = pos_CH[who].x +
			h1 * vel_CH[who].x -
			h2 * (a_H0[who].x - a_H1[who].x) +
         h3 * (a_H1[who1].x + a_H0[who1].x);
      pos_CH[who].y = pos_CH[who].y +
         h1 * vel_CH[who].y -
         h2 * (a_H0[who].y - a_H1[who].y) +
         h3 * (a_H1[who1].y + a_H0[who1].y);
      pos_CH[who].z = pos_CH[who].z +
         h1 * vel_CH[who].z -
         h2 * (a_H0[who].z - a_H1[who].z) +
         h3 * (a_H1[who1].z + a_H0[who1].z);
      vel_CH[who].x = vel_CH[who].x +
         h1 * (a_H0[who].x + a_H1[who].x) -
         h2 * (a_H0[who1].x - a_H1[who1].x) +
         h3 * (a_H0[who2].x + a_H1[who2].x);
	   vel_CH[who].y = vel_CH[who].y +
         h1 * (a_H0[who].y  + a_H1[who].y) -
         h2 * (a_H0[who1].y - a_H1[who1].y) +
         h3 * (a_H0[who2].y + a_H1[who2].y);
	   vel_CH[who].z = vel_CH[who].z +
         h1 * (a_H0[who].z  + a_H1[who].z) -
         h2 * (a_H0[who1].z - a_H1[who1].z) +
         h3 * (a_H0[who2].z + a_H1[who2].z);

		pos_CH[who].x += h1*vel_CH[who].x;
      pos_CH[who].y += h1*vel_CH[who].y;
      pos_CH[who].z += h1*vel_CH[who].z;

      p_v_a3[i].x = pos_CH[who].x;
	   p_v_a3[i].y = pos_CH[who].y;
	   p_v_a3[i].z = pos_CH[who].z;

      p_v_a3[i+nextsize].x = vel_CH[who].x;
	   p_v_a3[i+nextsize].y = vel_CH[who].y;
	   p_v_a3[i+nextsize].z = vel_CH[who].z;

      h2 = h1*h1;
      h3 = 1.0/(h2*h1);
      double h4 = h3/h1;
      double h5 = h4/h1;

      h3 *= 0.75;
      h4 *= 1.5;
      h5 *= 7.5;

      double Amin = a_H0[who].x - a_H1[who].x;
      double Jmin = h1 * (a_H0[who1].x - a_H1[who1].x);
      double Jplu = h1 * (a_H0[who1].x + a_H1[who1].x);
      double Smin = h2 * (a_H0[who2].x - a_H1[who2].x);
      double Splu = h2 * (a_H0[who2].x + a_H1[who2].x);

      a3_H[who].x = h3*(-5.0*Amin + 5.0*Jplu - Smin);
      double a4halfx = h4*(-Jmin + Splu);
      double a5halfx = h5*(3.0*Amin - 3.0*Jplu + Smin);
      a3_H[who].x += h1*a4halfx + h2/2.0*a5halfx;
      a4halfx += h1*a5halfx;

	   Amin = a_H0[who].y - a_H1[who].y;
      Jmin = h1 * (a_H0[who1].y - a_H1[who1].y);
      Jplu = h1 * (a_H0[who1].y + a_H1[who1].y);
      Smin = h2 * (a_H0[who2].y - a_H1[who2].y);
      Splu = h2 * (a_H0[who2].y + a_H1[who2].y);

      a3_H[who].y = h3*(-5.0*Amin + 5.0*Jplu - Smin);
      double a4halfy = h4*(-Jmin + Splu);
      double a5halfy = h5*(3.0*Amin - 3.0*Jplu + Smin);

      a3_H[who].y += h1*a4halfy + h2/2.0*a5halfy;
      a4halfy += h1*a5halfy;

      Amin = a_H0[who].z - a_H1[who].z;
      Jmin = h1 * (a_H0[who1].z - a_H1[who1].z);
      Jplu = h1 * (a_H0[who1].z + a_H1[who1].z);
      Smin = h2 * (a_H0[who2].z - a_H1[who2].z);
      Splu = h2 * (a_H0[who2].z + a_H1[who2].z);

      a3_H[who].z = h3*(-5.0*Amin + 5.0*Jplu - Smin);
      double a4halfz = h4*(-Jmin + Splu);
      double a5halfz = h5*(3.0*Amin - 3.0*Jplu + Smin);

      a3_H[who].z += h1*a4halfz + h2/2.0*a5halfz;
      a4halfz += h1*a5halfz;

      double a3dotsmod = sqrt(a3_H[who].x*a3_H[who].x + a3_H[who].y*a3_H[who].y + a3_H[who].z*a3_H[who].z);

      double a4mod = sqrt(a4halfx*a4halfx + a4halfy*a4halfy + a4halfz*a4halfz);
      double a5mod = sqrt(a5halfx*a5halfx + a5halfy*a5halfy + a5halfz*a5halfz);

      double    dt3 = (sqrt(amod*a2dotsmod) + adotmod) / (a5mod*a3dotsmod + a4mod*a4mod);
      dt3 = ETA6 * pow(dt3,1.0/6.0);

	   dt = dt3;

     double rest = *GTIME / (2.0 * step[who]);
     rest = (double)((int)(rest)) - rest;

	  if(dt<2.0e-5)
		  dt = max(sqrt(ETA4 * sqrt(amod) / sqrt(a2dotsmod)), dt3);

	  if(dt > 2.0*step[who] && rest == 0.0 && 2.0*step[who] <= DTMAX)
		  step[who] *= 2.0;
	  else if (dt < step[who] && 0.5*step[who] >= DTMIN)
		  step[who] *= 0.5;

	  p_v_a3[i+2*nextsize].x = a3_H[who].x;
	  p_v_a3[i+2*nextsize].y = a3_H[who].y;
	  p_v_a3[i+2*nextsize].z = a3_H[who].z;

     *ATIME = min (local_time[who] + step[who], *ATIME);

	  a_H1[who].x = a_H0[who].x;
     a_H1[who].y = a_H0[who].y;
     a_H1[who].z = a_H0[who].z;
     a_H1[who1].x = a_H0[who1].x;
     a_H1[who1].y = a_H0[who1].y;
     a_H1[who1].z = a_H0[who1].z;
     a_H1[who2].x = a_H0[who2].x;
     a_H1[who2].y = a_H0[who2].y;
     a_H1[who2].z = a_H0[who2].z;
	}

	return HNoError;
}


extern "C"
HostError GPU_memcheck(unsigned int NGPU, unsigned int *devices, const char *file, const int line){

	size_t free, total;
   double Mb_conv = 1024.0*1024.0;

   int local_rank;
   MPISafeCall(MPI_Comm_rank(MPI_COMM_WORLD, &local_rank));

   if(local_rank == 0){

      ofstream meminfo;
      meminfo.open("gpu_memory.dat", ios::app);
      if(!meminfo)
         return HNoFile;

   meminfo<<std::fixed<<std::setprecision(1);
   for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      DeviceSafeCall(cudaMemGetInfo( &free, &total ));
      meminfo<<" GPU # "<<i<<" : "<<" Used GPU memory @ file/line :"<<file<<"/"<<line<<" : "<<((double)total - (double)free)/((double)total) * 100.0<<" % "<<std::endl;
//    std::cout<<"            Free GPU memory : "<<(double)free/Mb_conv<<" MB "<<std::endl;
//    std::cout<<"           Total GPU memory : "<<(double)total/Mb_conv<<" MB "<<std::endl;

      if((double)free/Mb_conv < 20.0)
         return HNoMemory;
   }
   meminfo.close();
   }

   return HNoError;
}



HostError getEnergyLine_H6Blog(string *lastLine)
{

    ifstream data;
    data.open("HiGPUslog.dat");

	 if(!data)
		 return HNoFile;

	 char *name;

	getline(data, *lastLine);
	name = to_char(*lastLine);
		while(name[0]!='#'){
			getline(data, *lastLine);
			name = to_char(*lastLine);
		}

    data.close();

    return HNoError;
}



HostError AcquireEnergy(double *E){

	string last_line;
	HostSafeCall(getEnergyLine_H6Blog(&last_line));

	size_t found  = last_line.find('#',1);
	size_t found2 = last_line.find('#',found+1);

	if (found==string::npos)
		return HNotFound;

	string energy = last_line.substr(int(found+1), int(found2)-int(found+1));
	*E = to_double(energy);

	return HNoError;
}


extern "C"
HostError Hermite6th(const double TTIME, double* GTIME, double* ATIME, double* local_time, double* step, const unsigned int N, const unsigned int M, double4* pos_PH, float4* vel_PH, double4* pos_CH, double4* vel_CH, double4* a_H0, const unsigned int MAXDIM, unsigned int NGPU, unsigned int *devices, unsigned int TPB, int rank, int size, unsigned int BFMAX, double ETA6, double ETA4, double DTMAX, double DTMIN, double EPS, double DTPRINT, unsigned int FMAX, const bool warm, double GTW, unsigned int GPUMINTHREADS, double plummer_core, double plummer_mass){

	unsigned int ompthreads = 4; // N must be integer multiple of this number
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
	}


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

	double E0;
	int out_index = 1;

	if(warm){
		while(NEXTOUT <= GTW) NEXTOUT += DTPRINT;
		while(NEXTOUT <= *GTIME) NEXTOUT += DTPRINT;
		if(rank == 0)
	      HostSafeCall(AcquireEnergy(&E0));
	}
	else{
	   HostSafeCall(Calculate_Energy(pos_CD, vel_CD, N, EPS, TPB, NGPU, rank, devices, ppG, &E0, plummer_core, plummer_mass));
	}

	ofstream stream;

	if(rank == 0){
		stream.open("HiGPUslog.dat", ios::app);
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
	H6BTimes *Times;
	Times = new H6BTimes [N+1];
	for(unsigned int i = 0; i <= N; i++){
		Times[i].next_time = 0.0;
		Times[i].memcpy1_time = 0.0;
		Times[i].predictor_time = 0.0;
		Times[i].evaluation_time = 0.0;
		Times[i].reduce_time = 0.0;
		Times[i].reposition_time = 0.0;
		Times[i].memcpy2_time = 0.0;
		Times[i].mpireduce_time = 0.0;
		Times[i].corrector_time = 0.0;
		Times[i].reconstruct_time = 0.0;
		Times[i].energy_time = 0.0;
	}

	HostSafeCall(GPU_memcheck(NGPU, devices, __FILE__, __LINE__));
	HostSafeCall(CPU_memcheck(__FILE__, __LINE__));

	struct timeval tv;
	gettimeofday(&tv, NULL);
	int sec = tv.tv_sec;
	int microsec = tv.tv_usec;
	start_program = sec + microsec * 0.000001;
	
	do{
		
		if(*GTIME >= TTIME){
			
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
				DeviceSafeCall(cudaFree(p_v_a3_Dev[i]));
				DeviceSafeCall(cudaFree(a_temp_Dev[i]));
			}

			delete [] p_v_a3;
			delete [] a3_H;
			delete [] a_H;
        	delete [] a_H1;
        	delete [] mpi_red_aux;
        	delete [] mpi_red;
        	delete [] next;
			delete [] vetint;
			delete [] counter;
			
        	HostSafeCall(GPU_memcheck(NGPU, devices, __FILE__, __LINE__));
			HostSafeCall(CPU_memcheck(__FILE__, __LINE__));
			
			gettimeofday(&tv, NULL);
			int sec = tv.tv_sec;
			int microsec = tv.tv_usec;
			end_program = sec + microsec * 0.000001;

#ifdef CHECK_TIMES
			if(rank == 0){
				stream.open("times.dat", ios::app);
				stream<<scientific<<setprecision(3);
				stream<<" \n Total integration time : "<<end_program-start_program<<" seconds "<<endl;
				stream.close();
				stream.open("HiGPUslog.dat", ios::app);
				stream<<scientific<<setprecision(3);
				stream<<" \n Total integration time : "<<end_program-start_program<<" seconds "<<endl;
				stream.close();
			}
#endif

			if(rank == 0){
         	ofstream out_times;
            string str;
            str = "times_N_" + to_string(N) + "_nodes_" + to_string(size) + ".dat";
            char* file;
            file = to_char(str);
            out_times.open(file, ios::out);
            out_times<<"Executed with N = "<<N<<" , using "<<size<<" nodes "<<endl;
            out_times<<"global time of execution T = "<<scientific<<setprecision(4)<<end_program-start_program<<endl;
				out_times.close();
			}

			return HNoError;
		}

		get_times(&start);	
		HostSafeCall(NextParticles(N, ompthreads, counter, vetint, ATIME, local_time, step, next, &nextsize));

		*GTIME = *ATIME;
		
		unsigned int dim2 = ceil((double)nextsize/TPB)*TPB;
		for(unsigned int i = nextsize; i < dim2; i++)
			next[i] = -1;

		get_times(&end);
		set_times(end-start, &(Times[nextsize].next_time));

		get_times(&start);
		for(unsigned int i = 0; i < NGPU; i++){
	      DeviceSafeCall(cudaSetDevice(devices[i]));
		   int BL = ppG/TPB + ceil((double)nextsize/TPB);
			int istart = ppG*(i+rank*NGPU);

			DeviceSafeCall(cudaMemcpy(next_D[i], next,  dim2 * sizeof( int ), cudaMemcpyHostToDevice));

			Predictor <<<BL, TPB>>> (*GTIME, pos_PD[i], vel_PD[i], acc_PD[i], pos_CD[i], vel_CD[i], loc_D[i], a_tot_D[i], a1_tot_D[i], a2_tot_D[i], a3_D[i], istart, next_D[i], ppG, N);
		}

		get_times(&end);
		set_times(end-start, &(Times[nextsize].predictor_time));

 	   THREADS = TPB;
  		bfmax = BFMAX;
		*ATIME = 1.0e+10;

		Bfactor = 1;
		BLOCKS = ceil((double)nextsize/THREADS);

		while(nextsize*BLOCKS < GPUMINTHREADS && THREADS > 32){
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
 
		get_times(&start);
		for(unsigned int i = 0; i < NGPU; i++){
	      DeviceSafeCall(cudaSetDevice(devices[i]));
         DeviceCheckErrors();
	      int istart = ppG*(i+rank*NGPU);
	      evaluation<<< BLOCKS, THREADS, SHARED >>> ( N, pos_PD[i], vel_PD[i], acc_PD[i],  a_D[i], a1_D[i], a2_D[i],
                                       istart, ppG, Bfactor, DIMENSION, next_D[i], loc_D[i], *GTIME, EPS, plummer_core, plummer_mass);
	   }

		get_times(&end);
		set_times(end-start, &(Times[nextsize].evaluation_time));
		
		int bl = BLOCKS;
	   int bf = Bfactor;
	   SHARED = THREADS * sizeof(double4);


		get_times(&start);
	   while(bf != 1){
	      bl /= 2;
			for(unsigned int i = 0; i < NGPU; i++){
            DeviceSafeCall(cudaSetDevice(devices[i]));
            DeviceSafeCall(cudaThreadSynchronize());
            DeviceCheckErrors();
         }

	      for(unsigned int i = 0; i < NGPU; i++){
	         DeviceSafeCall(cudaSetDevice(devices[i]));
				reduce<<< bl, THREADS, SHARED>>>(a_D[i],  bf, DIMENSION);
	         reduce<<< bl, THREADS, SHARED>>>(a1_D[i], bf, DIMENSION);
   	      reduce<<< bl, THREADS, SHARED>>>(a2_D[i], bf, DIMENSION);
			}
	      bf /= 2;
	   }
		get_times(&end);
		set_times(end-start, &(Times[nextsize].reduce_time));

	
      for(unsigned int i = 0; i < NGPU; i++){
         DeviceSafeCall(cudaSetDevice(devices[i]));
         DeviceSafeCall(cudaThreadSynchronize());
      }	

		get_times(&start);
	   for(unsigned int i = 0; i < NGPU; i++){
	      DeviceSafeCall(cudaSetDevice(devices[i]));
			reposition<<<DIMENSION/THREADS, THREADS>>>(a_D[i], a_temp_Dev[i], 0, nextsize);
         reposition<<<DIMENSION/THREADS, THREADS>>>(a1_D[i], a_temp_Dev[i], 1, nextsize);
         reposition<<<DIMENSION/THREADS, THREADS>>>(a2_D[i], a_temp_Dev[i], 2, nextsize);
	   }
		get_times(&end);
		set_times(end-start, &(Times[nextsize].reposition_time));

	   unsigned int cpy_size = 3*nextsize;
	
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

		get_times(&start);
		HostSafeCall(Corrector(GTIME, ATIME, local_time, step, next, nextsize, pos_CH, vel_CH, a_H0,
         a_H1, a3_H, p_v_a3, ETA6, ETA4, DTMAX, DTMIN, N));
		get_times(&end);
		set_times(end-start, &(Times[nextsize].corrector_time));

		get_times(&start);
		for(unsigned int i = 0; i < NGPU; i++){
			DeviceSafeCall(cudaSetDevice(devices[i]));
			DeviceSafeCall(cudaMemcpy(p_v_a3_Dev[i], p_v_a3, 3*nextsize*sizeof( double4 ), cudaMemcpyHostToDevice ));
			DeviceSafeCall(cudaMemcpy(a_temp_Dev[i], a_H,    3*nextsize*sizeof( double4 ), cudaMemcpyHostToDevice ));

			int BB = DIMENSION/THREADS;

			Reconstruct<<< BB, THREADS >>>(next_D[i], nextsize, pos_CD[i], p_v_a3_Dev[i], 0);
			Reconstruct<<< BB, THREADS >>>(next_D[i], nextsize, vel_CD[i], p_v_a3_Dev[i], 1);
			Reconstruct<<< BB, THREADS >>>(next_D[i], nextsize, a3_D[i], p_v_a3_Dev[i], 2);
			Reconstruct<<< BB, THREADS >>>(next_D[i], nextsize, a_tot_D[i], a_temp_Dev[i], 0);
			Reconstruct<<< BB, THREADS >>>(next_D[i], nextsize, a1_tot_D[i], a_temp_Dev[i], 1);
			Reconstruct<<< BB, THREADS >>>(next_D[i], nextsize, a2_tot_D[i], a_temp_Dev[i], 2);
		}
		get_times(&end);
		set_times(end-start, &(Times[nextsize].reconstruct_time));

		if((*GTIME+GTW) >= NEXTOUT ){

			CheckBlocks(step, M);
			double E;
			get_times(&start);
			HostSafeCall(Calculate_Energy(pos_CD, vel_CD, N, EPS, TPB, NGPU, rank, devices, ppG, &E, plummer_core, plummer_mass));
			get_times(&end);
			set_times(end-start, &(Times[nextsize].energy_time));

			if(rank == 0){
				stream.open("times.dat", ios::out);
				stream<<scientific<<setprecision(3);
				stream<<"  N  "<<"  NEXT  "<<"	  PRED "<<"       EVAL "<<"      REDU "<<"     REPOS "<<"    CPY_ACC "<<"    MPI "<<"        CORR   "<<"    RECON "<<endl;  
				for(unsigned int i = 1; i <= M; i++){
					if(Times[i].next_time != 0.0)
						stream<<i<<"  "<<
							Times[i].next_time<<"  "<<
							Times[i].memcpy1_time<<"  "<<
							Times[i].predictor_time<<"  "<<
							Times[i].evaluation_time<<"  "<<
							Times[i].reduce_time<<"  "<<
							Times[i].reposition_time<<"  "<<
							Times[i].memcpy2_time<<"  "<<
							Times[i].mpireduce_time<<"  "<<
							Times[i].corrector_time<<"  "<<
							Times[i].memcpy3_time<<"  "<<
							Times[i].reconstruct_time<<"  "<<endl;
				}
				stream.close();

				stream.open("energy.dat", ios::app);
				stream<<scientific<<setprecision(5);
				stream<<*GTIME+GTW<<"  "<<(E-E0)/E0<<endl;
				stream.close();

				string file_name = to_string(out_index + FMAX);
				file_name += ".dat";
				stream.open(to_char(file_name), ios::out);
				stream<<scientific<<setprecision(16);
				for(unsigned int i = 0; i < M; i++)
					stream<<pos_CH[i].x<<"  "<<pos_CH[i].y<<"  "<<pos_CH[i].z<<"  "<<vel_CH[i].x<<"  "<<vel_CH[i].y<<"  "<<vel_CH[i].z<<"  "<<pos_CH[i].w<<endl;
				stream.close();
				out_index++;
			}

			NEXTOUT+=DTPRINT;
		}

	}while(1);

}

HostError Calculate_Energy(double4 **pos_CD, double4 **vel_CD, unsigned int N, double EPS, unsigned int THREADS, unsigned int NGPU, int rank, unsigned int *devices, unsigned int ppG, double *Energy, double plummer_core, double plummer_mass){

	double *E [NGPU];
	double *EH = new double [N*NGPU];
	double *Empi = new double [N];

	*Energy = 0.0;
        
	for(unsigned int i = 0; i < N; i++)
		Empi[i] = 0.0;

	for(unsigned int i = 0; i < NGPU; i++){
		DeviceSafeCall(cudaSetDevice(devices[i]));
		DeviceSafeCall(cudaMalloc(( void**)&E[i], N*sizeof(double)));
	}

	int BLOCKS = N/THREADS;
	int SHARED = THREADS*sizeof(double4);

	for(unsigned int i = 0; i < NGPU; i++){
		DeviceSafeCall(cudaSetDevice(devices[i]));
		DeviceSafeCall(cudaThreadSynchronize());
		int istart = ppG*(i+rank*NGPU);
		energy<<<BLOCKS, THREADS, SHARED>>>(pos_CD[i], vel_CD[i], E[i], N, EPS, istart, ppG, plummer_core, plummer_mass);
	}

	for(unsigned int i = 0; i < NGPU; i++){
		DeviceSafeCall(cudaSetDevice(devices[i]));
		DeviceSafeCall(cudaMemcpy(&EH[i*N], E[i], N*sizeof(double), cudaMemcpyDeviceToHost));
	}

	for(unsigned int i = 0; i < N; i++){
      for(unsigned int dev = 1; dev < NGPU; dev++){
         unsigned int p = i + dev*N;
      	EH[i] += EH[p];
		}
   }

	MPISafeCall(MPI_Allreduce(EH, Empi, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));

	for(unsigned int i = 0; i < N; i++)
		*Energy += Empi[i];

	delete [] Empi;
	delete [] EH;

	for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
		DeviceSafeCall(cudaFree(E[i]));
		
	}

	return HNoError;
}

HostError Calculate_potential_Energy(double4 *pos_CH, unsigned int N, double EPS, unsigned int THREADS, unsigned int NGPU, int rank, unsigned int *devices, unsigned int ppG, double *Energy, double plummer_core, double plummer_mass){

	double *E [NGPU];
	double4 **pos_CD = new double4* [NGPU];
	double *EH = new double [N*NGPU];
	double *Empi = new double [N];

	*Energy = 0.0;

	for(unsigned int i = 0; i < N; i++)
		Empi[i] = 0.0;

	for(unsigned int i = 0; i < NGPU; i++){
		DeviceSafeCall(cudaSetDevice(devices[i]));
		DeviceSafeCall(cudaMalloc(( void**)&E[i], N*sizeof(double)));
		DeviceSafeCall(cudaMalloc((void **)&pos_CD[i],  N*sizeof(double4)));	
      DeviceSafeCall(cudaMemcpy( pos_CD[i], pos_CH,     N*sizeof(double4),    cudaMemcpyHostToDevice ));
	}

	int BLOCKS = N/THREADS;
	int SHARED = THREADS*sizeof(double4);

	for(unsigned int i = 0; i < NGPU; i++){
		DeviceSafeCall(cudaSetDevice(devices[i]));
		DeviceSafeCall(cudaThreadSynchronize());
		int istart = ppG*(i+rank*NGPU);
		potential_energy<<<BLOCKS, THREADS, SHARED>>>(pos_CD[i], E[i], N, EPS, istart, ppG, plummer_core, plummer_mass);
	}

	for(unsigned int i = 0; i < NGPU; i++){
		DeviceSafeCall(cudaSetDevice(devices[i]));
		DeviceSafeCall(cudaMemcpy(&EH[i*N], E[i], N*sizeof(double), cudaMemcpyDeviceToHost));
	}

	for(unsigned int i = 0; i < N; i++){
		for(unsigned int dev = 1; dev < NGPU; dev++){
			unsigned int p = i + dev*N;
			EH[i] += EH[p];
		}
	}

	MPISafeCall(MPI_Allreduce(EH, Empi, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));

	for(unsigned int i = 0; i < N; i++)
		*Energy += Empi[i];

	delete [] Empi;
	delete [] EH;

	for(unsigned int i = 0; i < NGPU; i++){
		DeviceSafeCall(cudaSetDevice(devices[i]));
		DeviceSafeCall(cudaFree(E[i]));
		DeviceSafeCall(cudaFree(pos_CD[i]));
	}

	return HNoError;
}


HostError CudaInit(unsigned int *M, int NGPU, int rank, unsigned int *devices, string gpu_name){
	cudaDeviceProp *properties;
	int count;
	DeviceSafeCall(cudaGetDeviceCount(&count));
	properties = new cudaDeviceProp [count];
	ofstream hlog;
        
	if(rank == 0)
		hlog.open("HiGPUslog.dat", ios::app);

	if(count < NGPU || count <= 0)
		return HNoGpus;

	for(int i = 0; i < count; i++){
		DeviceSafeCall(cudaGetDeviceProperties(&properties[i], i));
		if(rank == 0)
			hlog<<" Available : "<<properties[i].name<<" as device : "<<i<<endl;
	}

	if(rank == 0)
		hlog<<"============================================="<<endl;
        
	int k = 0;
	for(int i = 0; i < count; i++){
		if(to_string(properties[i].name) != gpu_name)
			continue;
		else{
			devices[k] = i;
			k++;
		}
	}
	if(k<NGPU)
		return HNoGpus;

	if(rank==0) {
		for(int i = 0; i < NGPU; i++)
			hlog<<" Using : "<<properties[devices[i]].name<<" (device "<<devices[i]<<")"<<endl;
	}
	if(rank == 0){
		if(properties[devices[0]].major == 2)
			*M = properties[devices[0]].multiProcessorCount * 1536;
		else if(properties[devices[0]].major == 1){
			if(properties[devices[0]].minor == 3)
				*M = properties[devices[0]].multiProcessorCount * 1024;
			else
				return HNoDouble;
		}
	//	cout<<" Maximum number of parallel threads on the gpu : "<<*M<<endl;
	}

	MPISafeCall(MPI_Bcast(M, 1, MPI_INT, 0, MPI_COMM_WORLD));

	delete [] properties;

	if(rank == 0)
		hlog.close();

	return HNoError;
}

HostError DetermineSteps(double stp, unsigned int N, unsigned int M, double4* a_H0, double4* a_H1, double ETA4, double DTMAX, double DTMIN, double *step, double *ACTUAL_TIME, double *local_time){

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

HostError CheckBlocks(double *step, unsigned int M){

	int local_rank;

	MPISafeCall(MPI_Comm_rank(MPI_COMM_WORLD, &local_rank));
   
   for(unsigned int i = 0; i < M; i++)
	if(local_rank == 0){
      ofstream blocks;
		blocks.open("Blocks.dat");
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

extern "C"
HostError InitBlocks(double4 *pos_PH, float4 *vel_PH, unsigned int TPB, unsigned int N, unsigned int M, unsigned int BFMAX, double ETA4, unsigned int NGPU, double EPS, unsigned int *MAXDIM, double DTMAX, double DTMIN, unsigned int *GPUMINTHREADS, unsigned int *devices, string gpu_name, int rank, int nodes, double4* pos_CH, double4* vel_CH, double4* a_H0, double* step, double* local_time, double* ACTUAL_TIME, double plummer_core, double plummer_mass){
	
	HostSafeCall(CudaInit(GPUMINTHREADS, NGPU, rank, devices, gpu_name));
	
	HostSafeCall(Max_dimension(TPB, BFMAX, N, MAXDIM, *GPUMINTHREADS));
 	
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

	for(unsigned int i = 0; i < N; i++)
		pos_CH[i] = pos_PH[i];

	int * next = new int [N];
	for(unsigned int i = 0; i < N; i++)
		next[i] = i;
   
	unsigned long nextsize = N;
	int **next_D = new int* [NGPU];

	double *mpi_red_aux = new double [3*N];
	double *mpi_red = new double [3*N];
	double stp = 0.001;

	HostSafeCall(GPU_memcheck(NGPU, devices, __FILE__, __LINE__));
	HostSafeCall(CPU_memcheck(__FILE__, __LINE__));
   
	unsigned long BL = ceil((double)nextsize/TPB);
	int dim = TPB*BL;
	unsigned int threads = TPB;
	unsigned int bfmax = BFMAX;

	unsigned int Bfactor = 1;
	
	while(nextsize*BL < *GPUMINTHREADS && threads > 32){
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

	for(unsigned int i = 0; i < NGPU; i++){
		DeviceSafeCall(cudaSetDevice(devices[i]));
		int BLT = ceil((double)N/threads);
		initvectors<<<BLT, threads>>>(a3_D[i], acc_PD[i]);
	}
	
	HostSafeCall(GPU_memcheck(NGPU, devices, __FILE__, __LINE__));
	HostSafeCall(CPU_memcheck(__FILE__, __LINE__));

	unsigned int dim2 = ceil((double)nextsize/TPB)*TPB;

	for(unsigned int i = nextsize; i < dim2; i++)
		next[i] = -1;

	for(unsigned int i = 0; i < NGPU; i++){
		DeviceSafeCall(cudaSetDevice(devices[i]));
		DeviceSafeCall(cudaMemcpy( next_D[i], next, malloc_ui, cudaMemcpyHostToDevice ));
	}

	int ppG = N/(NGPU*nodes);
	int SHR = threads * (sizeof(double4) + 2 * sizeof(float4));
	
	for(unsigned int i = 0; i < NGPU; i++){
		DeviceSafeCall(cudaSetDevice(devices[i]));

		int istart = ppG*(i+rank*NGPU);
		evaluation<<< BL, threads, SHR >>> ( N, pos_PD[i], vel_PD[i], acc_PD[i],  a_D[i], a1_D[i], a2_D[i], 
													istart, ppG, Bfactor, dim, next_D[i], loc_D[i], 0.0, EPS, plummer_core, plummer_mass);
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
		bl /= 2;
		for(unsigned int i = 0; i < NGPU; i++){
			DeviceSafeCall(cudaSetDevice(devices[i]));
			reduce<<< bl, threads, SHR>>>(a_D[i],  bf, dim);
			reduce<<< bl, threads, SHR>>>(a1_D[i], bf, dim);
			reduce<<< bl, threads, SHR>>>(a2_D[i], bf, dim);
		}
		for(unsigned int i = 0; i < NGPU; i++){
			DeviceSafeCall(cudaSetDevice(devices[i]));
			DeviceSafeCall(cudaThreadSynchronize());
			DeviceCheckErrors();
		}
		bf /= 2;
	}

	for(unsigned int i = 0; i < NGPU; i++){
		DeviceSafeCall(cudaSetDevice(devices[i]));
		DeviceSafeCall(cudaThreadSynchronize());
		reposition<<<dim/threads, threads>>>(a_D[i], a_temp_Dev[i], 0, nextsize);
		reposition<<<dim/threads, threads>>>(a1_D[i], a_temp_Dev[i], 1, nextsize);
		reposition<<<dim/threads, threads>>>(a2_D[i], a_temp_Dev[i], 2, nextsize);   
	}
	
	unsigned int cpy_size = 3*nextsize;

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
                                       istart, ppG, Bfactor, dim, next_D[i], loc_D[i], 0.0, EPS, plummer_core, plummer_mass);
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
      bl /= 2;
      for(unsigned int i = 0; i < NGPU; i++){
         DeviceSafeCall(cudaSetDevice(devices[i]));
			reduce<<< bl, threads, SHR>>>(a_D[i],  bf, dim);
         reduce<<< bl, threads, SHR>>>(a1_D[i], bf, dim);
         reduce<<< bl, threads, SHR>>>(a2_D[i], bf, dim);
		}
      for(unsigned int i = 0; i < NGPU; i++){
         DeviceSafeCall(cudaSetDevice(devices[i]));
         DeviceSafeCall(cudaThreadSynchronize());
         DeviceCheckErrors();
      }
      bf /= 2;
   }
   
	for(unsigned int i = 0; i < NGPU; i++){
      DeviceSafeCall(cudaSetDevice(devices[i]));
      DeviceSafeCall(cudaThreadSynchronize());
		reposition<<<dim/threads, threads>>>(a_D[i], a_temp_Dev[i], 0, nextsize);
      reposition<<<dim/threads, threads>>>(a1_D[i], a_temp_Dev[i], 1, nextsize);
      reposition<<<dim/threads, threads>>>(a2_D[i], a_temp_Dev[i], 2, nextsize);   
	}

   cpy_size = 3*nextsize;

	for(unsigned int i = 0; i < NGPU; i++){
		DeviceSafeCall(cudaSetDevice(devices[i]));
      DeviceCheckErrors();
      DeviceSafeCall(cudaMemcpy(&a_H[i*cpy_size], a_temp_Dev[i], cpy_size*sizeof(double4), cudaMemcpyDeviceToHost));
   }
	
	HostSafeCall(ReduceAll(cpy_size, N, NGPU, nextsize,  a_H, a_H1, mpi_red_aux, mpi_red, next));
   
	HostSafeCall(DetermineSteps(stp, N, M, a_H0, a_H1, ETA4, DTMAX, DTMIN, step, ACTUAL_TIME, local_time));
   
	HostSafeCall(CheckBlocks(step, M));
	
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
	
	HostSafeCall(GPU_memcheck(NGPU, devices, __FILE__, __LINE__));
   HostSafeCall(CPU_memcheck(__FILE__, __LINE__));
  
	return HNoError;
}
