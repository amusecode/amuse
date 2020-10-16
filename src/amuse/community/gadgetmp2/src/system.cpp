/* ################################################################################## */
/* ###                                                                            ### */
/* ###                                 Gadgetmp2                                  ### */
/* ###                                                                            ### */
/* ###   Original: Gadget2 in the version used in Amuse                           ### */
/* ###   Author: Gadget2 and Amuse contributors                                   ### */
/* ###                                                                            ### */
/* ###   Modified: July 2020                                                      ### */
/* ###   Author: Thomas Schano                                                    ### */
/* ###                                                                            ### */
/* ###   Changes are intended to enable precise calculations in                   ### */
/* ###   non periodic small domain simulations in which comoving parts            ### */
/* ###   are simulated in std precision                                           ### */
/* ###                                                                            ### */
/* ################################################################################## */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#ifdef _WIN32
#else
#include <sys/resource.h>
#endif
#include <unistd.h>
#include <signal.h>
#include <gsl/gsl_rng.h>
#ifndef NOMPI
#include <mpi.h>
#endif

//#include "allvars.hpp"
#include "proto.hpp"


/*! \file system.c
 *  \brief contains miscellaneous routines, e.g. elapsed time measurements
 */


/*! This routine returns a random number taken from a table of random numbers,
 *  which is refilled every timestep.  This method is used to allow random
 *  number application to particles independent of the number of processors
 *  used, and independent of the particular order the particles have. In order
 *  to work properly, the particle IDs should be set properly to unique
 *  integer values.
 */
my_float gadgetmp2::get_random_number(int id)
{
  return RndTable[(id % RNDTABLE)];
}


/*! This routine fills the random number table.
 */
void gadgetmp2::set_random_numbers(void)
{
  int i;

  for(i = 0; i < RNDTABLE; i++)
    RndTable[i] = gsl_rng_uniform(random_generator);
}


/*! returns the number of cpu-ticks in seconds that have elapsed, or the
 *  wall-clock time obtained with MPI_Wtime().
 */
double gadgetmp2::second(void)
{
#ifdef WALLCLOCK
#ifndef NOMPI
  return MPI_Wtime();
#else
  return ( clock()) / CLOCKS_PER_SEC;
#endif
#else
  return (clock()) / CLOCKS_PER_SEC;
#endif

  /* note: on AIX and presumably many other 32bit systems,
   * clock() has only a resolution of 10ms=0.01sec
   */
}


/*! returns the time difference between two measurements obtained with
 *  second(). The routine takes care of the possible overflow of the tick
 *  counter on 32bit systems, but depending on the system, this may not always
 *  work properly. Similarly, in some MPI implementations, the MPI_Wtime()
 *  function may also overflow, in which case a negative time difference would
 *  be returned. The routine returns instead a time difference equal to 0.
 */
double gadgetmp2::timediff(double t0, double t1)
{
  double dt;

  dt = t1 - t0;

  if(dt < 0)	/* overflow has occured (for systems with 32bit tick counter) */
    {
#ifdef WALLCLOCK
      dt = 0;
#else
      dt = t1 + pow(2, 32) / CLOCKS_PER_SEC - t0;
#endif
    }

  return dt;
}


/*! returns the maximum of two my_float
 */
my_float gadgetmp2::dmax(my_float x, my_float y)
{
  if(x > y)
    return x;
  else
    return y;
}

/*! returns the minimum of two my_float
 */
my_float gadgetmp2::dmin(my_float x, my_float y)
{
  if(x < y)
    return x;
  else
    return y;
}

/*! returns the maximum of two integers
 */
int gadgetmp2::imax(int x, int y)
{
  if(x > y)
    return x;
  else
    return y;
}

/*! returns the minimum of two integers
 */
int gadgetmp2::imin(int x, int y)
{
  if(x < y)
    return x;
  else
    return y;
}
