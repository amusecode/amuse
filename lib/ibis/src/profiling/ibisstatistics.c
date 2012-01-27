#include "mpi.h"
#include "ibisstatistics.h"
#include <stdio.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

static int64_t *statistics;

static int socketfd = -1;

static int rank = 0;

static int size = 0;

static int add_count = 0;

//get port from file, specified by environment variable
//only used at rank 0
void ibis_statistics_get_ports(int ports[]) {
	char *portfile;
	FILE *portfilefp;
	int i;

	portfile = getenv("OMPI_IBIS_PROFILING_PORT_FILE");

	if (portfile == NULL) {
		fprintf(stderr, "statistics: mpi collector portfile not specified\n");
		return;
	}

	portfilefp = fopen(portfile, "r");

	if (portfilefp == NULL) {
		fprintf(stderr, "statistics: cannot open portfile %s\n", portfile);
		return;
	}

	for (i = 0; i < size; i++) {
		if (fscanf(portfilefp, "%d\n", &ports[i]) != 1) {
			fprintf(stderr,
					"statistics: cannot read port number from file %s\n",
					portfile);
			fclose(portfilefp);
			return;
		}
	}

	fclose(portfilefp);
}

void ibis_statistics_send_buffer(void *buffer, int length) {
	int total_written = 0;
	int bytes_written;

	if (socketfd <= 0) {
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

void ibis_statistics_socket_connect(int port) {
	struct sockaddr_in serv_addr;
	struct hostent *server;

	if (port <= 0) {
		fprintf(stderr, "statistics: ibis monitor port unknown\n");
		socketfd = -1;
		return;
	}

	socketfd = socket(AF_INET, SOCK_STREAM, 0);

	if (socketfd <= 0) {
		fprintf(stderr, "statistics: cannot open socket for monitoring\n");
		return;
	}

	/*fprintf(stderr, "statistics: connecting to socket on port %d\n", portno);*/

	server = gethostbyname("localhost");

	bzero((char *) &serv_addr, sizeof(serv_addr));
	serv_addr.sin_family = AF_INET;
	bcopy((char *) server->h_addr, (char *) &serv_addr.sin_addr.s_addr,
			server->h_length);
	serv_addr.sin_port = htons(port);
	if (connect(socketfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr))
			< 0) {
		/* set fd to negative number to denote it does not work */
		fprintf(stderr, "statistics: cannot connect socket\n");
		close(socketfd);
		socketfd = -1;
		return;
	}

	int32_t magic = MAGIC_NUMBER;

	ibis_statistics_send_buffer(&magic, sizeof(int32_t));
}

void ibis_statistics_init() {
	int port;
	int * ports;

	//set global rank and size
	PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
	PMPI_Comm_size(MPI_COMM_WORLD, &size);

	fprintf(stderr, "ibis_statistics_init(%d, %d)\n", rank, size);

	if (size < 0) {
		fprintf(stderr, "statistics: negative size pool!\n");
		socketfd = -1;
		return;
	}

	//allocate array for ports, get ports from file (if root)
	ports = calloc(size, sizeof(int));
	if (rank == 0) {
		ibis_statistics_get_ports(ports);
	}

	//broacast to all nodes, get local port from list, free list
	PMPI_Bcast(ports, size, MPI_INT, 0, MPI_COMM_WORLD);
	port = ports[rank];
	free(ports);

	//connect to monitoring socket, set socketfd to -1 if failure
	ibis_statistics_socket_connect(port);

	//inialize statistics to 0's
	statistics = calloc(size, sizeof(int64_t));
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

	/* send one last time */
	ibis_statistics_send_statistics();
	if (socketfd > 0) {
		close(socketfd);
	}
}

