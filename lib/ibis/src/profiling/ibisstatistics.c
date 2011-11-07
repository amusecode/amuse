#include "mpi.h"
#include "ibisstatistics.h"
#include <stdio.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <string.h>
#include <stdlib.h>


static char *statistic_names[STATS_TOTAL] = { 
   STATS_NAME_BARRIER, 
   STATS_NAME_SEND,
   STATS_NAME_RECV,
   STATS_NAME_ISEND,
   STATS_NAME_IRECV,
   STATS_NAME_BCAST,
   STATS_NAME_SCATTER,
   STATS_NAME_GATHER,
   STATS_NAME_ALLGATHER,
   STATS_NAME_ALLTOALL,
   STATS_NAME_REDUCE,
   STATS_NAME_ALLREDUCE,
   STATS_NAME_SCAN };

static stats statistics[STATS_MAX_COMM];

static int nextStat = 0;

static int socketfd = 0;

void statistics_socket_connect() {
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

    fprintf(stderr, "statistics: connecting to socket on port %d\n", portno);

    server = gethostbyname("localhost");

    bzero((char *) &serv_addr, sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    bcopy((char *) server->h_addr, (char *) &serv_addr.sin_addr.s_addr, server->h_length);
    serv_addr.sin_port = htons(portno);
    if (connect(socketfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) {
        /* set fd to negative number to denote it does not work */
        fprintf(stderr, "cannot connect socket\n");
        socketfd = -1;
        return;
    }

    write(socketfd, (char *)&portno, 4);
}


int create_communicator_statistics(MPI_Comm comm, int rank, int size)
{
   int i;


   fprintf(stderr,"creating statistics for %d\n", comm);

   statistics_socket_connect();

   if (nextStat == STATS_MAX_COMM) {
      return STATS_ERROR;
   }

   statistics[nextStat].comm = comm;
   statistics[nextStat].rank = rank;
   statistics[nextStat].size = size;

   for (i=0;i<STATS_TOTAL;i++) {
      statistics[nextStat].counters[i] = 0L;
   }

   nextStat++;

   return STATS_OK;
}

int inc_communicator_statistics(MPI_Comm comm, int field)
{
   int i;

/*   fprintf(stderr,"gathering statistics for %d\n", comm);*/

   if (field < 0 || field >= STATS_TOTAL) { 
      return STATS_ERROR;
   }

   for (i=0;i<nextStat;i++) {
      if (statistics[i].comm == comm) {
         statistics[i].counters[field]++;
         return STATS_OK;
      }
   }

   return STATS_NOT_FOUND;
}

static void print_comm_stats(int index) 
{
   char bufC[MPI_MAX_OBJECT_NAME+1];
   int lenC, i;

   MPI_Comm_get_name(statistics[index].comm, bufC, &lenC);

   fprintf(stderr, "MPIBIS: COMM STATS %s %d %d ", bufC, statistics[index].size, statistics[index].rank);

   for (i=0;i<STATS_TOTAL;i++) {
      fprintf(stderr, "/ %s %lu", statistic_names[i], statistics[index].counters[i]);
   }

   fprintf(stderr,"\n");
}

int print_communicator_statistics(MPI_Comm comm)
{
   int i;

   for (i=0;i<nextStat;i++) {
      if (statistics[i].comm == comm) {
         print_comm_stats(i);
         return STATS_OK;
      }
   }

   return STATS_NOT_FOUND;
}

int print_all_communicator_statistics()
{
   int i;
 
   for (i=0;i<nextStat;i++) {
      print_comm_stats(i);
   }

   return STATS_OK;
}

