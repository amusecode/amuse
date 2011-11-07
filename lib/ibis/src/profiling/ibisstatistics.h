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

#define MAGIC_NUMBER 449682

void ibis_statistics_init(int rank, int size);
void ibis_statistics_add(int type, int src_dst_root, int length);
void ibis_statistics_finalize();

#endif

