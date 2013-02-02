#ifndef MY_ERRORS_H
#define MY_ERRORS_H

#include "types.h"
#include <mpi.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <cstdlib>

#define HostSafeCall(err)      __mySafeCall( err, __FILE__, __LINE__)
#define MPISafeCall(err)         __myMPISafeCall(err, __FILE__, __LINE__)
#define DeviceSafeCall(err)    __myCudaSafeCall( err, __FILE__, __LINE__ )
#define DeviceCheckErrors()    __myCudaLastError(__FILE__, __LINE__)

inline const std::string __getErrorString(HostError err){

	std::string error_string = "Undefined string error";

	switch(err)
	{
		case(HNoError):
			error_string = " No error detected in function ;";
			break;
		case(HNoFile):
			error_string = " Cannot read the specified file ; ";
			break;	
		case(HInvalidArgv):
			error_string = " Invalid option passed to *argv[] : lunch with option -h to have an help ;";
			break;
		case(HNoLines):
			error_string = " Not enough lines in the specified file ;";
			break;
		case(HInvalidFormat):
			error_string = " Unrecognized format passed to function random_check_Bcast ;";
			break;
		case(HBcastFailed):
			error_string = " Broadcast failed ;";
			break;
		case(HInvalidN):
			error_string = " Invalid number of particles : it must be divisible by : \n 1) Number of computational nodes\n 2) Number of GPUs per node\n 3) Number of threads per block\n";
			break;
		case(HNoMemory):
			error_string = " Not enough memory on one of the available GPUs \n";
			break;
		case(HInvalidMeminfo):
			error_string = " Invalid reading of file /proc/meminfo, modify the function CPU_memcheck() \n";
			break;
		case(HNoGpus):
			error_string = " Not enough (or already tested) gpus available on the computational node \n";
			break;
		case(HNoDouble):
			error_string = " The gpu choosen does not support Double-precision floating-point numbers \n";
			break;
		case(HNotFound):
			error_string = " The name of the gpu in the input file must be terminated with the # character \n";
			break;
		case(HInvalidBlocks):
			error_string = " Invalid particle time step \n";
			break;
		case(HNoDefPlummer):
			error_string = " You did not define PLUMMER but you are passing option -p at lunch time : Please compile with -DPLUMMER if you want to add a plummer potential to the stellar environment, otherwise do not pass the option -p at lunch time \n";
			break;
		case(HNoDefPlummer2):
			error_string = " You did not pass the -p argument at lunch time but you compiled the code with the option -DPLUMMER. This implies you want to add a plummer potential to the stellar environment. Please add the -p option or compile without -DPLUMMER \n";
			case(HNoDefGalaxy):
         error_string = " You did not define GALAXY but you are passing option -gal at lunch time : Please compile with -DGALAXY if you want to add a Milky Way potential to the stellar environment, otherwise do not pass the option -gal at lunch time \n";
         break;
      case(HNoDefGalaxy2):
         error_string = " You did not pass the -gal argument at lunch time but you compiled the code with the option -DGALAXY. This implies you want to add a Milky Way potential to the stellar environment. Please add the -gal option or compile without -DGALAXY \n";
			break;
		case(HNoSpace):
			error_string = " There is not enough space in your Hard Disk \n";
			break;
		case(HNoNumber):
			error_string = " The number of GPUs specified with the option -d must be equal to that specified in the file 'input_param.txt' \n ";
			break;
	}

	return error_string;
}

inline void __myCudaSafeCall(cudaError err __attribute__((unused)), const char *file __attribute__((unused)), const int line __attribute__((unused))){

#ifdef CHECK_ERRORS

	do
	{
		if ( cudaSuccess != err ){
			std::cout<<" cudaSafeCall() FAILED at :"<<file<<"::"<<line<<":: "<<cudaGetErrorString(err)<<std::endl;
			exit(1);
		}
	}while(0);

#endif

	return;

}

inline void __myMPISafeCall(int err __attribute__((unused)), const char *file __attribute__((unused)), const int line __attribute__((unused))){

#ifdef CHECK_ERRORS

	char mpi_err_string[10000];
	int len;

	MPI_Error_string(err, mpi_err_string, &len);

	do
	{
		if (err)

		{
			std::cout<<" MPISafeCall() FAILED at : "<<file<<"::"<<line<<"\n"<<mpi_err_string<<std::endl;
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			exit(1);
		}
	}while(0);
#endif
	return;

}

inline void __myCudaLastError(const char *file __attribute__((unused)), const int line __attribute__((unused))){

#ifdef CHECK_ERRORS

	do{ 
		cudaError_t err = cudaGetLastError();
      
		if( cudaSuccess != err ){
			std::cout<<" cudaCheckError() FAILED at :"<<file<<"::"<<line<<":: "<<cudaGetErrorString(err)<<std::endl;
         exit( -1 );
      }
    } while ( 0 );

#endif
	return;
}

inline void __mySafeCall(HostError err __attribute__((unused)), const char *file __attribute__((unused)), const int line __attribute__((unused))){

#ifdef CHECK_ERRORS

    do
    {
        if ( err )
        {
			  std::cout<<" HostSafeCall() FAILED at : "<<file<<"::"<<line<<"\n"<<__getErrorString(err)<<std::endl;
           exit(1);
        }
		
    } while ( 0 );
#endif

    return;
}

#endif
