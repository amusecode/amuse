#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef NOMPI
#include <mpi.h>
#endif

#include "allvars.h"
#include "proto.h"


/*! \file longrange.c
 *  \brief driver routines for computation of long-range gravitational PM force
 */

#ifdef PMGRID

/*! Calls initializiation routines of periodic or/and non-periodic FFT
 *  routines.
 */
void long_range_init(void)
{
#ifdef PERIODIC
  pm_init_periodic();
#ifdef PLACEHIGHRESREGION
  pm_init_nonperiodic();
#endif
#else
  pm_init_nonperiodic();
#endif
}


/*! This function calls subroutines that determine the spatial region covered
 *  by the PM mesh.
 */
void long_range_init_regionsize(void)
{
#ifdef PERIODIC
#ifdef PLACEHIGHRESREGION
  if(RestartFlag != 1)
    pm_init_regionsize();
  pm_setup_nonperiodic_kernel();
#endif
#else
  if(RestartFlag != 1)
    pm_init_regionsize();
  pm_setup_nonperiodic_kernel();
#endif
}


/*! This function is a driver routine for the long-range PM force
 *  computation. It calls periodic and/or non-periodic FFT routines as needed
 *  for the present simulation set-up.
 */
void long_range_force(void)
{
  int i;

#ifndef PERIODIC
  int j;
  double fac;
#endif


  for(i = 0; i < NumPart; i++)
    P[i].GravPM[0] = P[i].GravPM[1] = P[i].GravPM[2] = 0;

#ifdef NOGRAVITY
  return;
#endif


#ifdef PERIODIC
  pmforce_periodic();
#ifdef PLACEHIGHRESREGION
  i = pmforce_nonperiodic(1);
  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmforce_nonperiodic(1);	/* try again */
    }
  if(i == 1)
    endrun(68686);
#endif
#else
  i = pmforce_nonperiodic(0);
  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmforce_nonperiodic(0);	/* try again */
    }
  if(i == 1)
    endrun(68687);
#ifdef PLACEHIGHRESREGION
  i = pmforce_nonperiodic(1);
  if(i == 1)			/* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();

      /* try again */

      for(i = 0; i < NumPart; i++)
	P[i].GravPM[0] = P[i].GravPM[1] = P[i].GravPM[2] = 0;

      i = pmforce_nonperiodic(0) + pmforce_nonperiodic(1);
    }
  if(i != 0)
    endrun(68688);
#endif
#endif


#ifndef PERIODIC
  if(All.ComovingIntegrationOn)
    {
      fac = 0.5 * All.Hubble * All.Hubble * All.Omega0;

      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].GravPM[j] += fac * P[i].Pos[j];
    }


  /* Finally, the following factor allows a computation of cosmological simulation
     with vacuum energy in physical coordinates */

  if(All.ComovingIntegrationOn == 0)
    {
      fac = All.OmegaLambda * All.Hubble * All.Hubble;

      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].GravPM[j] += fac * P[i].Pos[j];
    }
#endif

}


#endif
