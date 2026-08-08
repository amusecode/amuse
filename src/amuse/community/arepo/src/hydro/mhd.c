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
 * \file        src/mhd.c
 * \date        05/2018
 * \brief       Source terms for MHD implementation needed for cosmological
 *              MHD equations as well as Powell source terms.
 * \details     contains functions:
 *                void do_mhd_source_terms_first_half(void)
 *                void do_mhd_source_terms_second_half(void)
 *                void do_mhd_source_terms(void)
 *                void do_mhd_powell_source_terms(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 04.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef MHD

static void do_mhd_source_terms(void);

/*! \brief First half of the MHD source terms.
 *
 *  Before hydrodynamics timestep.
 *
 *  \return void
 */
void do_mhd_source_terms_first_half(void)
{
  do_mhd_source_terms();
  update_primitive_variables();
}

/*! \brief Second half of the MHD source terms.
 *
 *  After hydrodynamics timestep.
 *
 *  \return void
 */
void do_mhd_source_terms_second_half(void)
{
  do_mhd_source_terms();
  update_primitive_variables();
}

/*! \brief Adds source terms of MHD equations in expanding spacetime (i.e.
 *         in cosmological simulations) to energy.
 *
 *  \return void
 */
void do_mhd_source_terms(void)
{
  TIMER_START(CPU_MHD);

  if(All.ComovingIntegrationOn)
    {
      double atime    = All.Time;
      double hubble_a = hubble_function(atime);

      int idx, i;
      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;

          double dt_cell = 0.5 * (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval /
                           hubble_a; /* half the timestep of the cell */
          SphP[i].Energy += dt_cell * 0.5 * (SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]) *
                            SphP[i].Volume * atime * hubble_a;
        }
    }

  TIMER_STOP(CPU_MHD);
}

#endif /* #ifdef MHD */
