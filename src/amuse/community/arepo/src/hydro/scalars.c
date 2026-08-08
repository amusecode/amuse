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
 * \file        src/scalars.c
 * \date        05/2018
 * \brief       Routines to initialize passive scalars which are advected with
 *              the fluid.
 * \details     contains functions:
 *                void init_scalars()
 *                int scalar_init(MyFloat * addr, MyFloat * addr_mass, int
 *                  type)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 06.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../mesh/voronoi/voronoi.h"

#ifdef MAXSCALARS
int N_Scalar = 0;
struct scalar_elements scalar_elements[MAXSCALARS];
struct scalar_index ScalarIndex;
#endif /* #ifdef MAXSCALARS */

/*! \brief Main routine to initialize passive scalar quantities.
 *
 *  \return void
 */
void init_scalars()
{
#ifdef MAXSCALARS

#if defined(REFINEMENT_HIGH_RES_GAS)
  ScalarIndex.HighResMass = scalar_init(&SphP[0].HighResDensity, &SphP[0].HighResMass, SCALAR_TYPE_PASSIVE);
  if(ScalarIndex.HighResMass == -1)
    terminate("ScalarIndex.HighResMass initialized incorrectly\n");
#endif /* #if defined(REFINEMENT_HIGH_RES_GAS) */

#ifdef PASSIVE_SCALARS
  for(int i = 0; i < PASSIVE_SCALARS; i++)
    {
      scalar_init(&SphP[0].PScalars[i], &SphP[0].PConservedScalars[i], SCALAR_TYPE_PASSIVE);
    }
#endif /* #ifdef PASSIVE_SCALARS */

  mpi_printf("INIT: %d/%d Scalars used.\n", N_Scalar, MAXSCALARS);
#endif /* MAXSCALARS */
}

/*! \brief Initialize a specific scalar property.
 *
 *  \param[in] addr Pointer to (primitive) scalar in SphP[0] struct.
 *  \param[in] addr_mass Pointer to conserved scalar quantity in SphP[0].
 *  \param[in] type Type of scalar (e.g. SCALAR_TYPE_PASSIVE for passive
 *             scalar)
 *
 *  \return Number of scalars - 1
 */
int scalar_init(MyFloat *addr, MyFloat *addr_mass, int type)
{
#ifdef MAXSCALARS
  if(N_Scalar == MAXSCALARS)
    {
      mpi_printf("Failed to register scalar, maximum of %d already reached\n", MAXSCALARS);
      terminate("MAXSCALARS reached");
    }

  /* save type and relative address */
  scalar_elements[N_Scalar].type        = type;
  scalar_elements[N_Scalar].offset      = ((char *)addr) - ((char *)&SphP[0]);
  scalar_elements[N_Scalar].offset_mass = ((char *)addr_mass) - ((char *)&SphP[0]);

  N_Scalar++;

  return N_Scalar - 1;
  /* note: gradients are initialized in init_gradients */
#else  /* #ifdef MAXSCALARS */
  return -1;
#endif /* #ifdef MAXSCALARS #else */
}
