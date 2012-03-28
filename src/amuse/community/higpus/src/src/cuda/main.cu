#include <iostream>
#include <iomanip>
#include <string>
#include <mpi.h>


#include <my_errors.h>
#include <functions.h>
#include <utilis.h>

using namespace std;

int main(int argc, char *argv[]){

	MPISafeCall(MPI_Init(&argc, &argv));

	int rank, size;
	HostSafeCall(__MPIstart(&rank, &size, &mpi_float4, &mpi_double4));

	unsigned int N, NGPU, TPB, FMAX, BFMAX, MAXDIM, GPUMINTHREADS;
	unsigned int *devices;
	bool CDM, CDV, warm_start = 0;
	double EPS, ETA6, ETA4, TTIME, DTMAX, DTMIN, DTPRINT;
	string FINP, GPUNAME, warm_start_file;
	double plummer_core = 0.0;
	double plummer_mass = 0.0;

	if(rank == 0){

	string param;
	HostSafeCall(check_argv(argc, argv, &param, &warm_start, &warm_start_file, &plummer_core, &plummer_mass));

	if(!warm_start){
		ofstream generic;
		generic.open("gpu_memory.dat", ios::out);
		generic.close();
		generic.open("cpu_memory.dat", ios::out);
		generic.close();
		generic.open("H6Blog.dat", ios::out);
		generic.close();
		generic.open("energy.dat", ios::out);
		generic<<scientific<<setprecision(5);
		double value = 0.0;
		generic<<value<<"  "<<value<<endl;
		generic.close();
	}
	else{
		ofstream generic;
      generic.open("gpu_memory.dat", ios::app);
		generic<<" RESTART AT THIS POINT ***************************************************************"<<endl;
      generic.close();
      generic.open("cpu_memory.dat", ios::app);
		generic<<" RESTART AT THIS POINT ***************************************************************"<<endl;
      generic.close();
      generic.open("H6Blog.dat", ios::app);
		generic<<" RESTART AT THIS POINT ***************************************************************"<<endl;
      generic.close();
      generic.open("energy.dat", ios::app);
		generic<<"\n \n"<<endl;
      generic.close();
	}

	HostSafeCall(cpu_read_params(param, &N, &NGPU, &TPB, &TTIME, &DTMAX, &DTMIN, &EPS, &ETA6, &ETA4, &DTPRINT, &FMAX, &CDM, &CDV, &FINP, &GPUNAME));

	HostSafeCall(isDivisible(N, size, NGPU, TPB, &BFMAX));

#ifdef CHECK_ERRORS
	cout<<" Check errors option ENABLED "<<endl;
#else
	cout<<" Check errors option DISABLED "<<endl;
#endif

#ifdef CHECK_TIMES
	ofstream generic;
	generic.open("times.dat", ios::out);
	generic.close();
	cout<<" Check times option ENABLED "<<endl;
#else
	cout<<" Check times option DISABLED "<<endl;
#endif

#ifdef PLUMMER
	cout<<" Plummer Potential option ENABLED "<<endl;
	ofstream hlog;
	hlog.open("H6Blog.dat", ios::app);
	hlog<<" Plummer (core) : "<<plummer_core<<endl;
	hlog<<" Plummer (mass) : "<<plummer_mass<<endl;
	hlog<<" =============================================== "<<endl;
	hlog.close();
#else
	cout<<" Plummer Potential option DISABLED "<<endl;
#endif

	}
	
	MPISafeCall(MPI_Barrier(MPI_COMM_WORLD));

	HostSafeCall(__MPIbcastParam(&N, &NGPU, &TPB, &TTIME, &DTMAX, &DTMIN, &EPS, &ETA6, &ETA4, &DTPRINT, &BFMAX, &GPUNAME, rank, size)); // i parametri delle stampe sono SOLO di rank 0

	double4 *pos_PH = new double4 [N];
	float4  *vel_PH = new  float4 [N];
	double4 *pos_CH = new double4 [N];
	double4 *vel_CH = new double4 [N];
	double4 *a_H0 = new double4 [3*N];
	double  *step = new double [N];
	double  *local_time = new double [N];
	double ATIME = 1.0e+10;
	double GTIME = 0.0;
	double GTIME_WARM = 0.0;

	devices = new unsigned int [NGPU];

	if(rank == 0){
		if(warm_start){
			FINP = warm_start_file;
			CDM = 0;
			CDV = 0;
			unsigned int FMAX_old = FMAX;
			string loc = to_string(FMAX);
			loc = warm_start_file.substr(0, loc.length());
			FMAX = to_uint(loc);
			cout<<" FMAX NEW : "<<FMAX<<endl;
			GTIME_WARM = (FMAX-FMAX_old)*DTPRINT;
			cout<<" GTIME WARM : "<<GTIME_WARM<<endl;
		}
		HostSafeCall(cpu_read_external(FINP, pos_PH, vel_PH, vel_CH, N, CDM, CDV, warm_start));
	}

	MPISafeCall(MPI_Barrier(MPI_COMM_WORLD));
	MPISafeCall(MPI_Bcast(&warm_start, 1, MPI::BOOL, 0, MPI_COMM_WORLD));
	MPISafeCall(MPI_Bcast(&plummer_core, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
	MPISafeCall(MPI_Bcast(&plummer_mass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD));

   HostSafeCall(__MPIbcastPosVel(pos_PH, vel_PH, vel_CH, N, rank, size));
	
	HostSafeCall(InitBlocks(pos_PH, vel_PH, TPB, N, BFMAX, ETA4, NGPU, EPS, &MAXDIM, DTMAX, &GPUMINTHREADS, devices, GPUNAME, rank, size, pos_CH, vel_CH, a_H0, step, local_time, &ATIME, plummer_core, plummer_mass));
	
	HostSafeCall(Hermite6th(TTIME, &GTIME, &ATIME, local_time, step, N, pos_PH, vel_PH, pos_CH, vel_CH, a_H0, MAXDIM, NGPU, devices, TPB, rank, size, BFMAX, ETA6, ETA4, DTMAX, DTMIN, EPS, DTPRINT, FMAX, warm_start, GTIME_WARM, GPUMINTHREADS, plummer_core, plummer_mass));

	MPI_Finalize();

	return 0;
}
