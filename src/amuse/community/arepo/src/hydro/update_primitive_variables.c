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
 * \file        src/update_primitive_variables.c
 * \date        05/2018
 * \brief       Routines to recover the primitive hydrodynamical variables from
 *              the conserved ones.
 * \details     contains functions:
 *                void update_primitive_variables(void)
 *                void set_pressure_of_cell(int i)
 *                void set_pressure_of_cell_internal(struct particle_data
 *                  *localP, struct sph_particle_data *localSphP, int i)
 *                void do_validity_checks(struct particle_data *localP, struct
 *                  sph_particle_data *localSphP, int i, struct pv_update_data
 *                  *pvd)
 *                void update_primitive_variables_single(struct particle_data
 *                  *localP, struct sph_particle_data *localSphP, int i,
 *                  struct pv_update_data *pvd)
 *                void update_internal_energy(struct particle_data *localP,
 *                  struct sph_particle_data *localSphP, int i, struct
 *                  pv_update_data *pvd)
 *                double get_sound_speed(int p)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 11.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_linalg.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/*! \brief Main routine to update the primitive hydrodynamics variables from
 *         the conserved ones.
 *
 *  Note that the primitive variables are inconsistent with the (new)
 *  conserved variables after the hydro integration up to the point this
 *  function is called.
 *
 *  \return void
 */
void update_primitive_variables(void)
{
  TIMER_START(CPU_CELL_UPDATES);

  struct pv_update_data pvd;
  int idx, i;

  if(All.ComovingIntegrationOn)
    {
      pvd.atime    = All.Time;
      pvd.hubble_a = hubble_function(All.Time);
      pvd.a3inv    = 1 / (All.Time * All.Time * All.Time);
    }
  else
    pvd.atime = pvd.hubble_a = pvd.a3inv = 1.0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      do_validity_checks(P, SphP, i, &pvd);

      update_primitive_variables_single(P, SphP, i, &pvd);

      update_internal_energy(P, SphP, i, &pvd);

      set_pressure_of_cell_internal(P, SphP, i); /* calculate the pressure from Density and Utherm (and composition) */

      SphP[i].OldMass = P[i].Mass;

      SphP[i].TimeLastPrimUpdate = All.Time;
    }

  TIMER_STOP(CPU_CELL_UPDATES);
}

/*! \brief Wrapper function to calculate pressure of a cell from its internal
 *         energy.
 *
 *  \param[in] i Index of cell in P and SphP arrays.
 *
 *  \return void
 */
void set_pressure_of_cell(int i) { set_pressure_of_cell_internal(P, SphP, i); }

/*! \brief Function to calculate pressure from other hydrodynamics quantities.
 *
 *  How this is done depends on the adiabatic index and potentially on sub-
 *  resolution physics. Note that this is just the thermal pressure (i.e. not
 *  including magnetic fields).
 *
 *  \param[in] localP Pointer to particle data array.
 *  \param[in,out] localSphP Pointer to cell data array.
 *  \param[in] i Index in localP and localSphP arrays.
 *
 *  \return void
 */
void set_pressure_of_cell_internal(struct particle_data *localP, struct sph_particle_data *localSphP, int i)
{
#ifdef ISOTHERM_EQS
  localSphP[i].Pressure = localSphP[i].Density * All.IsoSoundSpeed * All.IsoSoundSpeed;
#else  /* #ifdef ISOTHERM_EQS */

  if(localSphP[i].Utherm >= 0)
    localSphP[i].Pressure = GAMMA_MINUS1 * localSphP[i].Density * localSphP[i].Utherm;
  else
    localSphP[i].Pressure = 0;
#endif /* #ifdef ISOTHERM_EQS */

#ifdef ENFORCE_JEANS_STABILITY_OF_CELLS
#if defined(USE_SFR)
  if(get_starformation_rate(i) == 0)
#endif /* #if defined(USE_SFR) */
    {
#ifdef ADAPTIVE_HYDRO_SOFTENING
      double cell_soft = All.ForceSoftening[localP[i].SofteningType];
#else  /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
    double cell_soft = All.GasSoftFactor * get_cell_radius(i);
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING #else */

      localSphP[i].Pressure =
          dmax(localSphP[i].Pressure, GAMMA_MINUS1 * localSphP[i].Density * 2 * All.G * localP[i].Mass / (All.cf_atime * cell_soft));
    }
#endif /* #ifdef ENFORCE_JEANS_STABILITY_OF_CELLS */
}

/*! \brief Validity checks for a gas cell.
 *
 *  So far, only a positive mass constraint implemented. Terminates if not
 *  successful.
 *
 *  \param[in] localP Pointer to particle data array
 *  \param[in,out] localSphP Pointer to cell data array
 *  \param[in] i Index in localP and localSphP arrays
 *  \param[in] pvd (unused)
 *
 *  \return void
 */
void do_validity_checks(struct particle_data *localP, struct sph_particle_data *localSphP, int i, struct pv_update_data *pvd)
{
  if(localP[i].Mass < 0)
    {
      printf("very bad...i=%d ID=%d mass=%g oldMass=%g utherm=%g pos=%g|%g|%g\n", i, (int)localP[i].ID, localP[i].Mass,
             localSphP[i].OldMass, localSphP[i].Utherm, localP[i].Pos[0], localP[i].Pos[1], localP[i].Pos[2]);

      terminate("stop");
    }
}

/*! \brief Updates primitive variables in a specified cell.
 *
 *  \param[in] localP Pointer to particle data array.
 *  \param[in,out] localSphP Pointer to cell data array.
 *  \param[in] i Index of cell in localP and localSphP arrays.
 *  \param[in] pvd additional data that is needed for update (e.g. cosmological
 *             factors).
 *
 *  \return void
 */
void update_primitive_variables_single(struct particle_data *localP, struct sph_particle_data *localSphP, int i,
                                       struct pv_update_data *pvd)
{
  localSphP[i].Density = localP[i].Mass / localSphP[i].Volume;

  if(localP[i].Mass > 0)
    {
      localP[i].Vel[0] = localSphP[i].Momentum[0] / localP[i].Mass;
      localP[i].Vel[1] = localSphP[i].Momentum[1] / localP[i].Mass;
      localP[i].Vel[2] = localSphP[i].Momentum[2] / localP[i].Mass;

#ifdef MAXSCALARS
      for(int k = 0; k < N_Scalar; k++)
        {
          *(MyFloat *)(((char *)(&localSphP[i])) + scalar_elements[k].offset) =
              *(MyFloat *)(((char *)(&localSphP[i])) + scalar_elements[k].offset_mass) / localP[i].Mass;
        }
#endif /* #ifdef MAXSCALARS */

#ifdef MHD
      localSphP[i].B[0] = localSphP[i].BConserved[0] / localSphP[i].Volume;
      localSphP[i].B[1] = localSphP[i].BConserved[1] / localSphP[i].Volume;
      localSphP[i].B[2] = localSphP[i].BConserved[2] / localSphP[i].Volume;
#endif /* #ifdef MHD */
    }
  else /* P[i].Mass <= 0 */
    {
      localP[i].Vel[0] = 0;
      localP[i].Vel[1] = 0;
      localP[i].Vel[2] = 0;

#ifdef MAXSCALARS
      for(int k = 0; k < N_Scalar; k++)
        *(MyFloat *)(((char *)(&localSphP[i])) + scalar_elements[k].offset) = 0;
#endif /* #ifdef MAXSCALARS */
    }
}

/*! \brief Updates the internal energy field in a specified cell
 *
 *  \param[in] localP Pointer to particle data array
 *  \param[in,out] localSphP Pointer to cell data array
 *  \param[in] i Index of cell in localP and localSphP arrays
 *  \param[in] pvd additional data that is needed for update (e.g. cosmological
 *             factors)
 *
 *  \return void
 */
void update_internal_energy(struct particle_data *localP, struct sph_particle_data *localSphP, int i, struct pv_update_data *pvd)
{
#ifndef ISOTHERM_EQS
  double ulimit;

  if(localP[i].Mass > 0)
    {
#ifdef MESHRELAX
      localSphP[i].Utherm = localSphP[i].Energy / localP[i].Mass;
#else  /* #ifdef MESHRELAX */
      localSphP[i].Utherm =
          (localSphP[i].Energy / localP[i].Mass -
           0.5 * (localP[i].Vel[0] * localP[i].Vel[0] + localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2])) /
          (pvd->atime * pvd->atime);
#endif /* #ifdef MESHRELAX #else */

#ifdef MHD
      localSphP[i].Utherm -=
          0.5 *
          (localSphP[i].B[0] * localSphP[i].B[0] + localSphP[i].B[1] * localSphP[i].B[1] + localSphP[i].B[2] * localSphP[i].B[2]) /
          localSphP[i].Density / pvd->atime;
#endif /* #ifdef MHD */

      ulimit = All.MinEgySpec;

      if(localSphP[i].Utherm < ulimit)
        {
          EgyInjection -= localSphP[i].Energy;

          localSphP[i].Utherm = ulimit;

#ifdef MESHRELAX
          localSphP[i].Energy = localP[i].Mass * localSphP[i].Utherm;
#else  /* #ifdef MESHRELAX */
          localSphP[i].Energy =
              pvd->atime * pvd->atime * localP[i].Mass * localSphP[i].Utherm +
              0.5 * localP[i].Mass *
                  (localP[i].Vel[0] * localP[i].Vel[0] + localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2]);
#endif /* #ifdef MESHRELAX */

#ifdef MHD
          localSphP[i].Energy +=
              0.5 *
              (localSphP[i].B[0] * localSphP[i].B[0] + localSphP[i].B[1] * localSphP[i].B[1] + localSphP[i].B[2] * localSphP[i].B[2]) *
              localSphP[i].Volume * pvd->atime;
#endif /* #ifdef MHD */

          EgyInjection += localSphP[i].Energy;
        }
    }
  else
    localSphP[i].Utherm = 0;

  if(localSphP[i].Density < All.LimitUBelowThisDensity && localSphP[i].Utherm > All.LimitUBelowCertainDensityToThisValue)
    {
      localSphP[i].Utherm = All.LimitUBelowCertainDensityToThisValue;
      localSphP[i].Energy =
          pvd->atime * pvd->atime * localP[i].Mass * localSphP[i].Utherm +
          0.5 * localP[i].Mass *
              (localP[i].Vel[0] * localP[i].Vel[0] + localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2]);
#ifdef MHD
      localSphP[i].Energy +=
          0.5 *
          (localSphP[i].B[0] * localSphP[i].B[0] + localSphP[i].B[1] * localSphP[i].B[1] + localSphP[i].B[2] * localSphP[i].B[2]) *
          localSphP[i].Volume * pvd->atime;
#endif /* #ifdef MHD */
    }

  if(localSphP[i].Utherm < 0)
    {
      printf("negative utherm %g\n", localSphP[i].Utherm);
      terminate("stop");
    }

#endif /* #ifndef ISOTHERM_EQS */
}

/*! \brief Calculates the sound speed of a specified cell
 *
 *  Depends on equation of state and potential sub-resolution physics.
 *
 *  \param[in] p Index of gas cell in P and SphP arrays
 *
 *  \return Sound speed
 */
double get_sound_speed(int p)
{
  double csnd;

#ifdef ISOTHERM_EQS
  csnd = All.IsoSoundSpeed;
#else  /* #ifdef ISOTHERM_EQS */

  double gamma;
  gamma = GAMMA;

  if(SphP[p].Density > 0)
    csnd = sqrt(gamma * SphP[p].Pressure / SphP[p].Density);
  else
    csnd = 0;
#endif /* #ifdef ISOTHERM_EQS #else */

#ifdef MHD
  /* for MHD, this is an upper bound to the signal velocity
     to do it more precisely, the magnet field in normal direction to the
     interfaces has to be taken into account */
  double Bsqr = SphP[p].B[0] * SphP[p].B[0] + SphP[p].B[1] * SphP[p].B[1] + SphP[p].B[2] * SphP[p].B[2];
  if(All.ComovingIntegrationOn)
    Bsqr /= All.Time;
  csnd = sqrt(csnd * csnd + Bsqr / SphP[p].Density);
#endif /* #ifdef MHD */

  return csnd;
}
