/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/gravity/longrange.c
 * \date        05/2018
 * \brief       Driver routines for computation of long-range gravitational
 *              PM force
 * \details     contains functions:
 *                void long_range_init(void)
 *                void long_range_init_regionsize(void)
 *                void long_range_force(void)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 06.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef PMGRID
/*! \brief Driver routine to call initialization of periodic or/and
 *         non-periodic FFT routines.
 *
 *  \return void
 */
void long_range_init(void)
{
#ifndef GRAVITY_NOT_PERIODIC
  pm_init_periodic();
#ifdef TWODIMS
  pm2d_init_periodic();
#endif /* #ifdef TWODIMS */
#ifdef PLACEHIGHRESREGION
  pm_init_nonperiodic();
#endif /* #ifdef PLACEHIGHRESREGION */
#else  /* #ifndef GRAVITY_NOT_PERIODIC */
  pm_init_nonperiodic();
#endif /* #ifndef GRAVITY_NOT_PERIODIC #else */
}

/*! \brief Driver routine to determine the extend of the non-
 *         periodic or high resolution region.
 *
 *  The initialization is done by pm_init_regionsize(). Afterwards
 *  the convolution kernels are computed by pm_setup_nonperiodic_kernel().
 *
 *  \return void
 */
void long_range_init_regionsize(void)
{
#ifndef GRAVITY_NOT_PERIODIC
#ifdef PLACEHIGHRESREGION
  if(RestartFlag != 1)
    pm_init_regionsize();
  pm_setup_nonperiodic_kernel();
#endif /* #ifdef PLACEHIGHRESREGION */

#else  /* #ifndef GRAVITY_NOT_PERIODIC */
  if(RestartFlag != 1)
    pm_init_regionsize();
  pm_setup_nonperiodic_kernel();
#endif /* #ifndef GRAVITY_NOT_PERIODIC #else */
}

/*! \brief This function computes the long-range PM force for all particles.
 *
 *  In case of a periodic grid the force is calculated by pmforce_periodic()
 *  otherwise by pmforce_nonperiodic(). If a high resolution region is
 *  specified for the PM force, pmforce_nonperiodic() calculates that force in
 *  both cases.
 *
 *  \return void
 */
void long_range_force(void)
{
  int i;

  TIMER_START(CPU_PM_GRAVITY);

#ifdef GRAVITY_NOT_PERIODIC
  int j;
  double fac;
#endif /* #ifdef GRAVITY_NOT_PERIODIC */

  for(i = 0; i < NumPart; i++)
    {
      P[i].GravPM[0] = P[i].GravPM[1] = P[i].GravPM[2] = 0;
#ifdef EVALPOTENTIAL
      P[i].PM_Potential = 0;
#endif /* #ifdef EVALPOTENTIAL */
    }

#ifndef SELFGRAVITY
  return;
#endif /* #ifndef SELFGRAVITY */

#ifndef GRAVITY_NOT_PERIODIC

#ifdef TWODIMS
  pm2d_force_periodic(0);
#else  /* #ifdef TWODIMS */
  pmforce_periodic(0, NULL);
#endif /* #ifdef TWODIMS #else */

#ifdef PLACEHIGHRESREGION
  i = pmforce_nonperiodic(1);

  if(i == 1) /* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmforce_nonperiodic(1); /* try again */
    }
  if(i == 1)
    terminate("despite we tried to increase the region, we still don't fit all particles in it");
#endif /* #ifdef PLACEHIGHRESREGION */

#else /* #ifndef GRAVITY_NOT_PERIODIC */
  i = pmforce_nonperiodic(0);

  if(i == 1) /* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();
      i = pmforce_nonperiodic(0); /* try again */
    }
  if(i == 1)
    terminate("despite we tried to increase the region, somehow we still don't fit all particles in it");
#ifdef PLACEHIGHRESREGION
  i = pmforce_nonperiodic(1);

  if(i == 1) /* this is returned if a particle lied outside allowed range */
    {
      pm_init_regionsize();
      pm_setup_nonperiodic_kernel();

      /* try again */

      for(i = 0; i < NumPart; i++)
        P[i].GravPM[0] = P[i].GravPM[1] = P[i].GravPM[2] = 0;

      i = pmforce_nonperiodic(0) + pmforce_nonperiodic(1);
    }
  if(i != 0)
    terminate("despite we tried to increase the region, somehow we still don't fit all particles in it");
#endif /* #ifdef PLACEHIGHRESREGION */
#endif /* #ifndef GRAVITY_NOT_PERIODIC #else */

#ifdef GRAVITY_NOT_PERIODIC
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
#endif /* #ifdef GRAVITY_NOT_PERIODIC */

  TIMER_STOP(CPU_PM_GRAVITY);

  find_long_range_step_constraint();
}
#endif /* #ifdef PMGRID */
