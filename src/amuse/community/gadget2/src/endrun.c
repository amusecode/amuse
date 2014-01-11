#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef NOMPI
#include <mpi.h>
#endif
#include <signal.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"


/*! \file endrun.c
 *  \brief Termination of simulation
 *
 *  This file contains routines for termination of the simulation.
 */

/*!  This function aborts the simulations. If a single processors wants an
 *   immediate termination, the function needs to be called with ierr>0. A
 *   bunch of MPI-error messages may also appear in this case.  For ierr=0,
 *   MPI is gracefully cleaned up, but this requires that all processors
 *   call endrun().
 */

// adapted for AMUSE:
void endrun(int ierr)
{
  if(ierr)
    {
      printf("task %d: endrun called with an error level of %d\n\n\n", ThisTask, ierr);
      fflush(stdout);
    }
  exit(0);
}

void original_endrun(int ierr)
{
  if(ierr)
    {
      printf("task %d: endrun called with an error level of %d\n\n\n", ThisTask, ierr);
      fflush(stdout);
#ifdef DEBUG
      terminate_processes();
      raise(SIGABRT);
      sleep(60);
#else
#ifndef NOMPI
      MPI_Abort(MPI_COMM_WORLD, ierr);
#endif
#endif
      exit(0);
    }

#ifndef NOMPI
  MPI_Finalize();
#endif
  exit(0);
}
