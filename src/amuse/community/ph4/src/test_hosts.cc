#include "stdinc.h"

// Standalone parallel function to return the local rank of the
// calling process and the number of sister processes on the same node
// (as defined by the HOST or HOSTNAME environment variables).

void get_local_rank(MPI_Intracomm comm,
		    int& mylocalrank, int& mylocalsize,
		    bool verbose = false)
{
    MPI_Barrier(comm);

    int rank;
    MPI_Comm_rank(comm, &rank);
    int size;
    MPI_Comm_size(comm, &size);

    // Get my hostname.

    const int NHOSTNAME = 1024;
    char hostname[NHOSTNAME];
    char *h = getenv("HOST");
    if (!h) h = getenv("HOSTNAME");
    if (h)
	strcpy(hostname, h);
    else
	strcpy(hostname, "unknown host");

    if (rank > 0) {

	// Send my hostname to process 0.

	int nh = strlen(hostname);
	MPI_Send( &nh, 1, MPI_INT, 0, 990, comm);
	MPI_Send( hostname, strlen(hostname), MPI_CHAR, 0, 991, comm);

	// Get my local size and rank from process 0.

	MPI_Recv(&mylocalrank, 1, MPI_INT, 0, 992, comm);
	MPI_Recv(&mylocalsize, 1, MPI_INT, 0, 993, comm);

    } else {

	// Collect all hostnames.

	int nh = strlen(hostname);
	char hostnames[size][NHOSTNAME];
	for (int i = 0; i < nh; i++)
	    hostnames[0][i] = hostname[i];
	hostnames[0][nh] = '\0';

	for (int ip = 1; ip < size; ip++) {
	    MPI_Recv( &nh, 1, MPI_INT, ip, 990, comm);
	    MPI_Recv( hostnames[ip], nh, MPI_CHAR, ip, 991, comm);
	    hostnames[ip][nh] = '\0';	    
	}

	// Count the total number of processes sharing a host with
	// each process, and determine the "local rank" of each.

	int localsize[size], locallist[size][size];
	for (int ip = 0; ip < size; ip++) localsize[ip] = 0;

	// A full n^2 loop here ensures that the local lists will be
	// ordered.

	for (int ip = 0; ip < size; ip++)
	    for (int jp = 0; jp < size; jp++)
		if (!strcmp(hostnames[ip], hostnames[jp]))
		    locallist[ip][localsize[ip]++] = jp;

	// Determine local ranks.

	int localrank[size];
	for (int ip = 0; ip < size; ip++)
	    for (int i = 0; i < localsize[ip]; i++)
		if (locallist[ip][i] == ip) localrank[ip] = i;

	if (verbose) {
	    cout << "hosts:" << endl;
	    for (int ip = 0; ip < size; ip++) {
		cout << "    process " << ip << " on " << hostnames[ip]
		     << ", local np = " << localsize[ip] << ":";
		for (int i = 0; i < localsize[ip]; i++)
		    cout << " " << locallist[ip][i];
		cout << " (" << localrank[ip] << ")" << endl;
	    }
	}

	mylocalrank = localrank[0];
	mylocalsize = localsize[0];

	// Distribute localsize and localrank to all processes.

	for (int ip = 1; ip < size; ip++) {
	    MPI_Send( localrank+ip, 1, MPI_INT, ip, 992, comm);
	    MPI_Send( localsize+ip, 1, MPI_INT, ip, 993, comm);
	}
    }

    MPI_Barrier(comm);
}

void run_test(bool verbose = false)
{
    MPI_Intracomm comm = MPI_COMM_WORLD;
    int mylocalrank, mylocalsize;
    get_local_rank(comm, mylocalrank, mylocalsize, verbose);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, & rank)
    PRC(rank); PRL(mylocalrank);
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    run_test();
    MPI_Finalize();
}
