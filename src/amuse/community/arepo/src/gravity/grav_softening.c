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
 * \file        src/gravity/gravtree.c
 * \date        05/2018
 * \brief       Routines for setting the gravitational softening lengths.
 * \details     contains functions:
 *                void set_softenings(void)
 *                int get_softeningtype_for_hydro_cell(int i)
 *                double get_default_softening_of_particletype(int type)
 *                int get_softening_type_from_mass(double mass)
 *                double get_desired_softening_from_mass(double mass)
 *                void init_individual_softenings(void)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 06.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

/*! \brief Sets the (comoving) softening length of all particle
 *         types in the table All.SofteningTable[...].
 *
 *  A check is performed that the physical softening length is bounded by the
 *  Softening-MaxPhys values.
 *
 *  \return void
 */
void set_softenings(void)
{
  int i;

  if(All.ComovingIntegrationOn)
    {
      for(i = 0; i < NSOFTTYPES; i++)
        if(All.SofteningComoving[i] * All.Time > All.SofteningMaxPhys[i])
          All.SofteningTable[i] = All.SofteningMaxPhys[i] / All.Time;
        else
          All.SofteningTable[i] = All.SofteningComoving[i];
    }
  else
    {
      for(i = 0; i < NSOFTTYPES; i++)
        All.SofteningTable[i] = All.SofteningComoving[i];
    }

#ifdef ADAPTIVE_HYDRO_SOFTENING
  for(i = 0; i < NSOFTTYPES_HYDRO; i++)
    All.SofteningTable[i + NSOFTTYPES] = All.MinimumComovingHydroSoftening * pow(All.AdaptiveHydroSofteningSpacing, i);

  if(All.AdaptiveHydroSofteningSpacing < 1)
    terminate("All.AdaptiveHydroSofteningSpacing < 1");

#ifdef MULTIPLE_NODE_SOFTENING
  /* we check that type=0 has its own slot 0 in the softening types, so that only gas masses are stored there */
  if(All.SofteningTypeOfPartType[0] != 0)
    terminate("All.SofteningTypeOfPartType[0] != 0");

  for(i = 1; i < NTYPES; i++)
    if(All.SofteningTypeOfPartType[i] == All.SofteningTypeOfPartType[0])
      terminate("i=%d: All.SofteningTypeOfPartType[i] == All.SofteningTypeOfPartType[0]", i);
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */

  for(i = 0; i < NSOFTTYPES + NSOFTTYPES_HYDRO; i++)
    All.ForceSoftening[i] = 2.8 * All.SofteningTable[i];

  All.ForceSoftening[NSOFTTYPES + NSOFTTYPES_HYDRO] = 0; /* important - this entry is actually used */
}

#ifdef ADAPTIVE_HYDRO_SOFTENING
/*! \brief Finds the index of the softening table for a given cell depending
 *         on its radius.
 *
 *  \param[in] i Index of cell in SphP array.
 *
 *  \return Index of corresponding softening in softening lookup-table.
 */
int get_softeningtype_for_hydro_cell(int i)
{
  double soft = All.GasSoftFactor * get_cell_radius(i);

  if(soft <= All.ForceSoftening[NSOFTTYPES])
    return NSOFTTYPES;

  int k = 0.5 + log(soft / All.ForceSoftening[NSOFTTYPES]) / log(All.AdaptiveHydroSofteningSpacing);
  if(k >= NSOFTTYPES_HYDRO)
    k = NSOFTTYPES_HYDRO - 1;

  return NSOFTTYPES + k;
}
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */

/*! \brief Returns the default softening length for particle type 'type'.
 *
 * \param[in] type Type of the local particle.
 *
 * \return The softening length of particle with type 'type'.
 */
double get_default_softening_of_particletype(int type) { return All.SofteningTable[All.SofteningTypeOfPartType[type]]; }

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
/*! \brief Determines the softening type from the mass of a particle.
 *
 *  \param[in] mass Mass of the particle.
 *
 *  \return Index in gravitational softening table.
 */
int get_softening_type_from_mass(double mass)
{
  int i, min_type = -1;
  double eps     = get_desired_softening_from_mass(mass);
  double min_dln = MAX_FLOAT_NUMBER;

#if defined(MULTIPLE_NODE_SOFTENING) && defined(ADAPTIVE_HYDRO_SOFTENING)
  i = 1;
#else  /* #if defined(MULTIPLE_NODE_SOFTENING) && defined(ADAPTIVE_HYDRO_SOFTENING) */
  i = 0;
#endif /* #if defined(MULTIPLE_NODE_SOFTENING) && defined(ADAPTIVE_HYDRO_SOFTENING) #else */

  for(; i < NSOFTTYPES; i++)
    {
      if(All.ForceSoftening[i] > 0)
        {
          double dln = fabs(log(eps) - log(All.ForceSoftening[i]));

          if(dln < min_dln)
            {
              min_dln  = dln;
              min_type = i;
            }
        }
    }
  if(min_type < 0)
    terminate("min_type < 0  mass=%g  eps=%g   All.AvgType1Mass=%g  All.ForceSoftening[1]=%g", mass, eps, All.AvgType1Mass,
              All.ForceSoftening[1]);

  return min_type;
}

/*! \brief Returns the softening length of softening type 1
 *  particles depending on the particle mass.
 *
 *  \param[in] mass Particle mass.
 *
 *  \return Softening length for a softening type 1 particle of mass 'mass'.
 */
double get_desired_softening_from_mass(double mass)
{
  if(mass <= All.AvgType1Mass)
    return 2.8 * All.SofteningComoving[1];
  else
    return 2.8 * All.SofteningComoving[1] * pow(mass / All.AvgType1Mass, 1.0 / 3);
}

/*! \brief Initializes the mass dependent softening calculation for Type 1
 *         particles.
 *
 *  The average mass of Type 1 particles is calculated.
 *
 *  \return void
 */
void init_individual_softenings(void)
{
  int i, ndm;
  double mass, masstot;
  long long ndmtot;

  for(i = 0, ndm = 0, mass = 0; i < NumPart; i++)
    if(P[i].Type == 1)
      {
        ndm++;
        mass += P[i].Mass;
      }
  sumup_large_ints(1, &ndm, &ndmtot);
  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  All.AvgType1Mass = masstot / ndmtot;

  mpi_printf("INIT: AvgType1Mass = %g\n", All.AvgType1Mass);

  for(i = 0; i < NumPart; i++)
    {
      if(((1 << P[i].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
        P[i].SofteningType = get_softening_type_from_mass(P[i].Mass);
    }
}
#endif /* #ifdef INDIVIDUAL_GRAVITY_SOFTENING */
