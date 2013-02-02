#include <iostream>
#include <iomanip>
#include <string>
#include <mpi.h>
#include <my_errors.h>
#include <functions.h>
#include <utilis.h>

using namespace std;

int main(int argc, char *argv[]){

	// MPI initialization ////////////////////////////////////////////////////////
	int rank, size;
	MPISafeCall(MPI_Init(&argc, &argv));
	HostSafeCall(__MPIstart(&rank, &size, &mpi_float4, &mpi_double4));
	//////////////////////////////////////////////////////////////////////////////

	// General Parameters : declaration //////////////////////////////////////////////////////
	unsigned int			N, M, NGPU, TPB, FMAX, BFMAX, MAXDIM, GPUMINTHREADS;
	bool						CDM, CDV, VIR = 0, warm_start = 0, setdev = 0, cleanstop = 0;
	double					EPS, ETA6, ETA4, TTIME, DTMAX, DTMIN, DTPRINT, plummer_core = 0.0, plummer_mass = 0.0, ratio = 0.0, mscale = 0.0, rscale = 0.0;
	string					FINP, GPUNAME, warm_start_file, argument = " ", path = "./";
	vector<unsigned int> dev; // it will contain the devices that will be used for the simulations
   /////////////////////////////////////////////////////////////////////////////////////////

	if(rank == 0){
		string param;
		// read command line options (in argv[][])
		HostSafeCall(check_argv(argc, argv, &param, &warm_start, &warm_start_file, &VIR, &ratio, &setdev, dev, &plummer_core, &plummer_mass, &mscale, &rscale, &argument, &cleanstop));

		if(!cleanstop){
			// open output files
		   if(!warm_start)
	         HostSafeCall(open_file(argument, path));
	      else
	         HostSafeCall(append_file(argument, path));
	
	      // read simulation parameters from file (input_param.txt)
	      HostSafeCall(cpu_read_params(param, &N, &NGPU, &TPB, &TTIME, &DTMAX, &DTMIN, &EPS,
	                                    &ETA6, &ETA4, &DTPRINT, &FMAX, &CDM, &CDV, &FINP, &GPUNAME, path));
	      // check input correctness
	      HostSafeCall(isDivisible(&N, &M, size, NGPU, TPB, &BFMAX));
	
	      // print some info
	      HostSafeCall(print_info(plummer_core, plummer_mass, rscale, mscale, path));
		}
   }

	MPISafeCall(MPI_Barrier(MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast(&cleanstop, 1, MPI::BOOL, 0, MPI_COMM_WORLD));

	if(cleanstop){
		MPI_Finalize();
		return 0;
	}
	
	MPISafeCall(MPI_Barrier(MPI_COMM_WORLD));
	// Broadcast parameters among MPI processes
	HostSafeCall(__MPIbcastParam(&N, &M, &NGPU, &TPB, &TTIME, &DTMAX, &DTMIN, &EPS, &ETA6, &ETA4, &DTPRINT, &BFMAX, &GPUNAME, rank, size));

	// Maindeclarations //////////////////////////////////////////////////////////////
	double4 *pos_PH		= new double4 [N];
	float4  *vel_PH		= new  float4 [N];
	double4 *pos_CH		= new double4 [N];
	double4 *vel_CH		= new double4 [N];
	double4 *a_H0			= new double4 [3*N];
	double  *step			= new double [N];
	double  *local_time  = new double [N];
	double ATIME = 1.0e+10, GTIME = 0.0, GTIME_WARM = 0.0;
	///////////////////////////////////////////////////////////////////////////////////

   if(rank == 0){
      if(warm_start)
         HostSafeCall(adjust_param_ifwarmstart(&CDM, &CDV, &FMAX, warm_start_file, path, &FINP, &GTIME_WARM, DTPRINT));
      // Read positions, velocities and masses from external file : initial conditions
      HostSafeCall(cpu_read_external(FINP, pos_PH, vel_PH, vel_CH, N, M, CDM, CDV, warm_start));
   }

	// Broadcast (again) new parameters and data among MPI processes
   HostSafeCall(__MPIbcast_otherparams(&warm_start, &VIR, &setdev, &plummer_core, &plummer_mass, &rscale, &mscale, &ratio, dev));
   HostSafeCall(__MPIbcastPosVel(pos_PH, vel_PH, vel_CH, N, rank, size));

	// Initialize Blocks for Hermite Block Time Steps integration
   HostSafeCall(InitBlocks(pos_PH, vel_PH, TPB, N, M, BFMAX, ETA4, DTMIN, DTMAX, NGPU, EPS, &MAXDIM, &GPUMINTHREADS, GPUNAME, rank, size, pos_CH, vel_CH, a_H0, step, local_time, &ATIME, VIR, ratio, warm_start, setdev, dev, plummer_core, plummer_mass, rscale, mscale, path));

   // Start integration
   HostSafeCall(Hermite6th(TTIME, &GTIME, &ATIME, local_time, step, N, M, pos_PH, vel_PH, pos_CH, vel_CH, a_H0, MAXDIM, NGPU, TPB, rank, size, BFMAX, ETA6, ETA4, DTMAX, DTMIN, EPS, DTPRINT, FMAX, warm_start, GTIME_WARM, GPUMINTHREADS, plummer_core, plummer_mass, rscale, mscale, dev, &cleanstop, path));

	MPISafeCall(MPI_Barrier(MPI_COMM_WORLD));
   MPISafeCall(MPI_Bcast(&cleanstop, 1, MPI::BOOL, 0, MPI_COMM_WORLD));

   if(cleanstop){
      MPI_Finalize();
      return 0;
   }


	
	MPI_Finalize();

	return 0;
}
