#include "stdinc.h"

// Standalone parallel function to return the local rank of the
// calling process and the number of sister processes on the same node
// (as defined by the HOST or HOSTNAME environment variables).

void get_local_rank(MPI_Intracomm comm,
		    int& mylocalrank, int& mylocalsize,
		    bool verbose = false)
{
    MPI_Barrier(comm);

    int rank = MPI_Get_rank(comm);
    int size = MPI_Get_size(comm);

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
	MPI_Send(comm, &nh, 1, MPI_INT, 0, 990);
	MPI_Send(comm, hostname, strlen(hostname), MPI_CHAR, 0, 991);

	// Get my local size and rank from process 0.

	MPI_Recv(comm, &mylocalrank, 1, MPI_INT, 0, 992);
	MPI_Recv(comm, &mylocalsize, 1, MPI_INT, 0, 993);

    } else {

	// Collect all hostnames.

	int nh = strlen(hostname);
	char hostnames[size][NHOSTNAME];
	for (int i = 0; i < nh; i++)
	    hostnames[0][i] = hostname[i];
	hostnames[0][nh] = '\0';

	for (int ip = 1; ip < size; ip++) {
	    MPI_Recv(comm, &nh, 1, MPI_INT, ip, 990);
	    MPI_Recv(comm, hostnames[ip], nh, MPI_CHAR, ip, 991);
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
	    MPI_Send(comm, localrank+ip, 1, MPI_INT, ip, 992);
	    MPI_Send(comm, localsize+ip, 1, MPI_INT, ip, 993);
	}
    }

    MPI_Barrier(comm);
}

void run_test(bool verbose = false)
{
    MPI_Intracomm comm = MPI_COMM_WORLD;
    int mylocalrank, mylocalsize;
    get_local_rank(comm, mylocalrank, mylocalsize, verbose);
    PRC(MPI_COMM_WORLD.Get_rank()); PRL(mylocalrank);
}

int main(int argc, char *argv[])
{
    MPI_Init(argc, argv);
    run_test();
    MPI_Finalize();
}
