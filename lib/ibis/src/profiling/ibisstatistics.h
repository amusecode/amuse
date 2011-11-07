#ifndef IBISSTATISTICS_H
#define IBISSTATISTICS_H

#include "mpi.h"

#define STATS_BARRIER   0
#define STATS_SEND      1
#define STATS_RECV      2
#define STATS_ISEND     3
#define STATS_IRECV     4
#define STATS_BCAST     5
#define STATS_SCATTER   6
#define STATS_GATHER    7
#define STATS_ALLGATHER 8
#define STATS_ALLTOALL  9
#define STATS_REDUCE    10
#define STATS_ALLREDUCE 11
#define STATS_SCAN      12
#define STATS_TOTAL     13

#define STATS_NAME_BARRIER   "barrier"
#define STATS_NAME_SEND      "send"
#define STATS_NAME_RECV      "receive"
#define STATS_NAME_ISEND     "isend"
#define STATS_NAME_IRECV     "ireceive"
#define STATS_NAME_BCAST     "bcast"
#define STATS_NAME_SCATTER   "scatter"
#define STATS_NAME_GATHER    "gather"
#define STATS_NAME_ALLGATHER "allgather"
#define STATS_NAME_ALLTOALL  "alltoall"
#define STATS_NAME_REDUCE    "reduce"
#define STATS_NAME_ALLREDUCE "allreduce"
#define STATS_NAME_SCAN      "scan"

#define STATS_MAX_COMM 32

#define STATS_OK        0
#define STATS_ERROR     1
#define STATS_NOT_FOUND 2

typedef struct s {
   MPI_Comm comm;
   int rank;
   int size;
   unsigned long counters[STATS_TOTAL];
} stats;

int create_communicator_statistics(MPI_Comm comm, int rank, int size);
int inc_communicator_statistics(MPI_Comm comm, int field);
int print_communicator_statistics(MPI_Comm comm);
int print_all_communicator_statistics();

#endif

