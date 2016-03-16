#ifndef _AMUSE_MPI_H_
#define _AMUSE_MPI_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <mpi.h>

int set_comm_world(MPI_Comm comm);
int get_comm_world(MPI_Comm * comm);

#ifdef __cplusplus
}
#endif

#endif
