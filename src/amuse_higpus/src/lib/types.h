#ifndef TYPES_H
#define TYPES_H

#include <mpi.h>

typedef enum { HNoError, HNoFile, HInvalidArgv, HNoLines, HInvalidFormat, HBcastFailed, HInvalidN, HNoMemory, HInvalidMeminfo, HNoGpus, HNoDouble, HNotFound, HInvalidBlocks, HNoDefPlummer, HNoDefPlummer2, HNoSpace, HNoNumber, HNoDefGalaxy, HNoDefGalaxy2} HostError; //new type

extern MPI_Datatype mpi_double4; // variable
extern MPI_Datatype mpi_float4; // variable

struct HiGPUsTimes {
	double next_time;
	double rtime;
	double bfac;
	double thr;
	double totthr;
	double memcpy1_time;
	double cpynext_time;
	double predictor_time;
	double evaluation_time;
	double reduce_time;
	double reposition_time;
	double memcpy2_time;
	double mpireduce_time;
	double corrector_time;
	double memcpy3_time;
	double reconstruct_time;
	double energy_time;
};


#endif
