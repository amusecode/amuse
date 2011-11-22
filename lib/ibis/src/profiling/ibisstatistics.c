#include "mpi.h"
#include "ibisstatistics.h"
#include <stdio.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <string.h>
#include <stdlib.h>

static int64_t *statistics;

static int socketfd = -1;

static int rank = 0;

static int size = 0;

static int add_count = 0;

void ibis_statistics_send_buffer(void *buffer, int length) {
	int total_written = 0;
	int bytes_written;

	if (socketfd <= 0) {
		fprintf(stderr, "could not write data, socket not valid\n");
		return;
	}

	while (total_written < length) {
		bytes_written = write(socketfd, ((char *) buffer) + total_written,
				length - total_written);

		if (bytes_written == -1) {
			fprintf(stderr, "could not write data\n");
			socketfd = -1;
			return;
		}

		total_written = total_written + bytes_written;
	}
}

/* send rank of process, total size, and sent count per process */
void ibis_statistics_send_statistics() {
	int i;

	ibis_statistics_send_buffer(&rank, sizeof(int32_t));
	ibis_statistics_send_buffer(&size, sizeof(int32_t));
	ibis_statistics_send_buffer(statistics, size * sizeof(int64_t));

	for (i = 0; i < size; i++) {
		/*fprintf(stderr, "rank %2d send %10ld bytes to rank %2d\n", rank, statistics[i], i);*/
		statistics[i] = 0;
	}
}

void ibis_statistics_socket_connect() {
	int portno;
	struct sockaddr_in serv_addr;
	struct hostent *server;
	char *portenv;

	socketfd = socket(AF_INET, SOCK_STREAM, 0);

	if (socketfd <= 0) {
		fprintf(stderr, "statistics: cannot open socket for MPI monitoring\n");
		return;
	}

	portenv = getenv("IBIS_MPI_COLLECTOR_PORT");

	if (portenv == NULL) {
		fprintf(stderr, "statistics: mpi collector port unknown\n");
		socketfd = -1;
		return;
	}

	portno = atoi(portenv);

	if (portno <= 0) {
		fprintf(stderr, "statistics: mpi collector port unknown\n");
		socketfd = -1;
		return;
	}

	/*fprintf(stderr, "statistics: connecting to socket on port %d\n", portno);*/

	server = gethostbyname("localhost");

	bzero((char *) &serv_addr, sizeof(serv_addr));
	serv_addr.sin_family = AF_INET;
	bcopy((char *) server->h_addr, (char *) &serv_addr.sin_addr.s_addr,
			server->h_length);
	serv_addr.sin_port = htons(portno);
	if (connect(socketfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr))
			< 0) {
		/* set fd to negative number to denote it does not work */
		fprintf(stderr, "statistics: cannot connect socket\n");
		socketfd = -1;
		return;
	}

	int32_t magic = MAGIC_NUMBER;

	ibis_statistics_send_buffer(&magic, sizeof(int32_t));
}

void ibis_statistics_init(int myrank, int mysize) {
	int i;
	rank = myrank;
	size = mysize;

	fprintf(stderr, "ibis_statistics_init(%d, %d)\n", rank, size);

	if (size < 0) {
		fprintf(stderr, "statistics: negative size pool!\n");
		socketfd = -1;
		return;
	}

	statistics = malloc(size * sizeof(int64_t));

	for (i = 0; i < size; i++) {
		statistics[i] = 0;
	}

	ibis_statistics_socket_connect();

}

void ibis_statistics_add_all_others(int length) {
	int i;

	for (i = 0; i < size; i++) {
		if (i != rank) {
			statistics[i] += length;
		}
	}
}

void ibis_statistics_add(int type, int src_dst_root, int length) {
	/*	fprintf(stderr, "ibis_statistics_add(%d, %d, %d)\n", type, src_dst_root, length);*/

	switch (type) {
	case STATS_BARRIER:
		/* guestimate, mostly looks nice in a visualization */
		ibis_statistics_add_all_others(4);
		break;
	case STATS_SEND:
	case STATS_ISEND:
		statistics[src_dst_root] += length;
		break;
	case STATS_BCAST:
	case STATS_SCATTER:
		if (rank == src_dst_root) {
			ibis_statistics_add_all_others(length);
		}
		break;
	case STATS_GATHER:
	case STATS_REDUCE:
		if (rank != src_dst_root) {
			statistics[src_dst_root] += length;
		}
		break;
	case STATS_ALLGATHER:
	case STATS_ALLTOALL:
	case STATS_ALLREDUCE:
	case STATS_SCAN:
		ibis_statistics_add_all_others(length);
		break;
	case STATS_RECV:
	case STATS_IRECV:
		/* we only track sending bytes, not receiving them */
		break;
	default:
		fprintf(stderr, "ibis_statistics_add(%d, %d, %d): don't know type %d\n",
				type, src_dst_root, length, type);
		break;
	}

	add_count++;

	/* send statistics once in a while */
	if (add_count > 100) {
		ibis_statistics_send_statistics();
		add_count = 0;
	}
}

void ibis_statistics_finalize() {
	int i;
	fprintf(stderr, "ibis_statistics_finalize()\n");

	/* send one last time */
	ibis_statistics_send_statistics();
	if (socketfd > 0) {
		close(socketfd);
	}
}

