#include "mpi.h"
#include "ibisstatistics.h"
#include <stdio.h>

/*****************************************************************************/
/*                        Initialization / Finalization                      */
/*****************************************************************************/

int MPI_Init(int *argc, char **argv[]) {
	int result, rank, size;

	result = PMPI_Init(argc, argv);

	PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
	PMPI_Comm_size(MPI_COMM_WORLD, &size);

	ibis_statistics_init(rank, size);

	return result;
}

int MPI_Init_thread(int *argc, char ***argv, int required, int *provided) {
	int result, rank, size;

	result = PMPI_Init_thread(argc, argv, required, provided);

	PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
	PMPI_Comm_size(MPI_COMM_WORLD, &size);

	ibis_statistics_init(rank, size);

	return result;
}

int MPI_Finalize(void) {
	ibis_statistics_finalize();
	return PMPI_Finalize();
}

int MPI_Abort(MPI_Comm comm, int errorcode) {
	ibis_statistics_finalize();
	return PMPI_Abort(comm, errorcode);
}

/*****************************************************************************/
/*                             Send / Receive                                */
/*****************************************************************************/

int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
		MPI_Comm comm, MPI_Request *request) {
	ibis_statistics_add(STATS_ISEND, dest, count);
	return PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
}

int MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
		MPI_Comm comm) {
	ibis_statistics_add(STATS_SEND, dest, count);
	return PMPI_Send(buf, count, datatype, dest, tag, comm);
}

int MPI_Ssend(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
		MPI_Comm comm) {
	ibis_statistics_add(STATS_SEND, dest, count);
	return PMPI_Ssend(buf, count, datatype, dest, tag, comm);
}

int MPI_Rsend(void* buf, int count, MPI_Datatype datatype, int dest, int tag,
		MPI_Comm comm) {
	ibis_statistics_add(STATS_SEND, dest, count);
	return PMPI_Rsend(buf, count, datatype, dest, tag, comm);
}

int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
		MPI_Comm comm, MPI_Request *request) {
	ibis_statistics_add(STATS_IRECV, source, count);
	return PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
}

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
		MPI_Comm comm, MPI_Status *status) {
	ibis_statistics_add(STATS_RECV, source, count);
	return PMPI_Recv(buf, count, datatype, source, tag, comm, status);
}

/*****************************************************************************/
/*                             Waits / Polling                               */
/*****************************************************************************/

int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status) {
	return PMPI_Test(request, flag, status);
}

int MPI_Wait(MPI_Request *request, MPI_Status *status) {
	return PMPI_Wait(request, status);
}

int MPI_Waitany(int count, MPI_Request *array_of_requests, int *index,
		MPI_Status *status) {
	return PMPI_Waitany(count, array_of_requests, index, status);
}

int MPI_Waitall(int count, MPI_Request *array_of_requests,
		MPI_Status *array_of_statuses) {
	return PMPI_Waitall(count, array_of_requests, array_of_statuses);
}

/*****************************************************************************/
/*                              Collectives                                  */
/*****************************************************************************/

int MPI_Barrier(MPI_Comm comm) {
	ibis_statistics_add(STATS_BARRIER, 0, 0);

	return PMPI_Barrier(comm);
}

int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype, int root,
		MPI_Comm comm) {
	ibis_statistics_add(STATS_BCAST, root, count);

	return PMPI_Bcast(buffer, count, datatype, root, comm);
}

int MPI_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
		void* recvbuf, int recvcount, MPI_Datatype recvtype, int root,
		MPI_Comm comm) {
	int typesize;
	PMPI_Type_size(sendtype, &typesize);

	ibis_statistics_add(STATS_GATHER, root, sendcount * typesize);

	return PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
			recvtype, root, comm);
}

int MPI_Gatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		void *recvbuf, int *recvcounts, int *displs, MPI_Datatype recvtype,
		int root, MPI_Comm comm) {
	int typesize;
	PMPI_Type_size(sendtype, &typesize);
	ibis_statistics_add(STATS_GATHER, root, sendcount * typesize);

	return PMPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts,
			displs, recvtype, root, comm);
}

int MPI_Allgather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
		void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
	int typesize;
	PMPI_Type_size(sendtype, &typesize);
	ibis_statistics_add(STATS_ALLGATHER, 0, sendcount * typesize);

	return PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
			recvtype, comm);
}

int MPI_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		void *recvbuf, int *recvcounts, int *displs, MPI_Datatype recvtype,
		MPI_Comm comm) {
	int typesize;
	PMPI_Type_size(sendtype, &typesize);
	ibis_statistics_add(STATS_ALLGATHER, 0, sendcount * typesize);

	return PMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts,
			displs, recvtype, comm);
}

int MPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype,
		void* recvbuf, int recvcount, MPI_Datatype recvtype, int root,
		MPI_Comm comm) {
	int typesize;
	PMPI_Type_size(sendtype, &typesize);
	ibis_statistics_add(STATS_SCATTER, root, sendcount * typesize);

	return PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount,
			recvtype, root, comm);
}

int MPI_Scatterv(void* sendbuf, int *sendcounts, int *displs,
		MPI_Datatype sendtype, void* recvbuf, int recvcount,
		MPI_Datatype recvtype, int root, MPI_Comm comm) {
	int size, i, typesize;
	int sendcount = 0;

	PMPI_Type_size(sendtype, &typesize);
	PMPI_Comm_size(MPI_COMM_WORLD, &size);

	for (i = 0; i < size; i++) {
		sendcount += sendcounts[i];
	}

	ibis_statistics_add(STATS_SCATTER, root, sendcount * typesize);

	return PMPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf,
			recvcount, recvtype, root, comm);
}

int MPI_Reduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
		MPI_Op op, int root, MPI_Comm comm) {
	ibis_statistics_add(STATS_REDUCE, root, count);

	return PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
}

int MPI_Allreduce(void* sendbuf, void* recvbuf, int count,
		MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
	ibis_statistics_add(STATS_ALLREDUCE, 0, count);

	return PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}

int MPI_Scan(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
		MPI_Op op, MPI_Comm comm) {
	ibis_statistics_add(STATS_SCAN, 0, count);

	return PMPI_Scan(sendbuf, recvbuf, count, datatype, op, comm);
}

int MPI_Alltoall(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
	int typesize;
	PMPI_Type_size(sendtype, &typesize);
	ibis_statistics_add(STATS_ALLTOALL, 0, sendcount * typesize);

	return PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount,
			recvtype, comm);
}

int MPI_Alltoallv(void *sendbuf, int *sendcounts, int *sdispls,
		MPI_Datatype sendtype, void *recvbuf, int *recvcounts, int *rdispls,
		MPI_Datatype recvtype, MPI_Comm comm) {
	int size, i, typesize;
	int sendcount = 0;

	PMPI_Type_size(sendtype, &typesize);
	PMPI_Comm_size(MPI_COMM_WORLD, &size);

	for (i = 0; i < size; i++) {
		sendcount += sendcounts[i];
	}

	ibis_statistics_add(STATS_ALLTOALL, 0, sendcount * typesize);

	return PMPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf,
			recvcounts, rdispls, recvtype, comm);
}

/*****************************************************************************/
/*                         Communicators and Groups                          */
/*****************************************************************************/

int MPI_Comm_free(MPI_Comm *comm) {
	return PMPI_Comm_free(comm);
}

int MPI_Comm_size(MPI_Comm comm, int *size) {
	return PMPI_Comm_size(comm, size);
}

int MPI_Comm_rank(MPI_Comm comm, int *rank) {
	return PMPI_Comm_rank(comm, rank);
}

int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm) {
	return PMPI_Comm_dup(comm, newcomm);
}

int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm) {
	return PMPI_Comm_create(comm, group, newcomm);
}

int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm) {
	return PMPI_Comm_split(comm, color, key, newcomm);
}

int MPI_Comm_group(MPI_Comm comm, MPI_Group *group) {
	return PMPI_Comm_group(comm, group);
}

int MPI_Group_incl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup) {
	return PMPI_Group_incl(group, n, ranks, newgroup);
}

int MPI_Group_range_incl(MPI_Group group, int n, int ranges[][3],
		MPI_Group *newgroup) {
	return PMPI_Group_range_incl(group, n, ranges, newgroup);
}

int MPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup) {
	return PMPI_Group_union(group1, group2, newgroup);
}

int MPI_Group_intersection(MPI_Group group1, MPI_Group group2,
		MPI_Group *newgroup) {
	return PMPI_Group_intersection(group1, group2, newgroup);
}

int MPI_Group_difference(MPI_Group group1, MPI_Group group2,
		MPI_Group *newgroup) {
	return PMPI_Group_difference(group1, group2, newgroup);
}

int MPI_Group_free(MPI_Group *group) {
	return PMPI_Group_free(group);
}

int MPI_Group_translate_ranks(MPI_Group group1, int n, int *ranks1,
		MPI_Group group2, int *ranks2) {
	return PMPI_Group_translate_ranks(group1, n, ranks1, group2, ranks2);
}

/*****************************************************************************/
/*                                Utilities                                  */
/*****************************************************************************/

double MPI_Wtime(void) {
	return PMPI_Wtime();
}

int MPI_Error_string(int errorcode, char *string, int *resultlen) {
	return PMPI_Error_string(errorcode, string, resultlen);
}

int MPI_Get_processor_name(char *name, int *resultlen) {
	return PMPI_Get_processor_name(name, resultlen);
}

int MPI_Initialized(int *flag) {
	return PMPI_Initialized(flag);
}

int MPI_Pack(void *inbuf, int incount, MPI_Datatype datatype, void *outbuf,
		int outsize, int *position, MPI_Comm comm) {
	return PMPI_Pack(inbuf, incount, datatype, outbuf, outsize, position, comm);
}

int MPI_Unpack(void *inbuf, int insize, int *position, void *outbuf,
		int outcount, MPI_Datatype datatype, MPI_Comm comm) {
	return PMPI_Unpack(inbuf, insize, position, outbuf, outcount, datatype,
			comm);
}
