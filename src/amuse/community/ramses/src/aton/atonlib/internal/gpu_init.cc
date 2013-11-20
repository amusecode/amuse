#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <cuda_runtime_api.h>

#define HOSTNAME_BUFFER_SIZE 50

extern "C" int aton_init_gpu_(const int* allow_gpu_overload) {
  fprintf(stderr, "init_cuda_device\n");

  int mpi_size, mpi_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  int gpu_device_count;
  cudaGetDeviceCount(&gpu_device_count);

  char my_hostname[HOSTNAME_BUFFER_SIZE] = {0};
  if (gethostname(my_hostname, HOSTNAME_BUFFER_SIZE) != 0) {
    perror("initdevice error:");
    return 0;
  }

  char* all_hostnames = (char*) malloc(mpi_size * HOSTNAME_BUFFER_SIZE);
  
  MPI_Allgather(my_hostname, HOSTNAME_BUFFER_SIZE, MPI_CHAR,
		all_hostnames, HOSTNAME_BUFFER_SIZE, MPI_CHAR,
		MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    fprintf(stderr, "mpi_size = %d\n", mpi_size);
    int i;
    for (i = 0; i < mpi_size; ++i) {
      fprintf(stderr, "%d: %s\n", i, &all_hostnames[i*HOSTNAME_BUFFER_SIZE]);
    }
  }

  // Count the number of hostnames equal to our own. Also find the position of
  // our process amongs the other processes running on this host.
  int count = 0;
  int my_index = 0;
  int i;
  for (i = 0; i < mpi_size; ++i) {
    if (strncmp(my_hostname,
		&all_hostnames[i*HOSTNAME_BUFFER_SIZE],
		HOSTNAME_BUFFER_SIZE) == 0) {
      if (i == mpi_rank) {
	my_index = count;
      }
      count++;
    }
  }
  free(all_hostnames);
  all_hostnames = NULL;

  fprintf(stderr, "%d mpi processes on hostname '%s'.\n", count, my_hostname);

  if (count > gpu_device_count) {
    fprintf(stderr, "initdevice.c warning: %d mpi processes on hostname '%s' but "
	   "only %d gpu devices.\n", count, my_hostname, gpu_device_count);
    if (*allow_gpu_overload == 0) {
      fprintf(stderr, "initdevice.c error: This warning is an error because "
	     "allow_gpu_overload=false.\n");
      return 0;
    }
  }

  int device_num = my_index % gpu_device_count;
  cudaSetDevice(device_num);

  return 1;
}
