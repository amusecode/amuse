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
 * \file        src/time_integration/do_gravity_hydro.c
 * \date        05/2018
 * \brief       Contains the two half step kick operators.
 * \details     This file contains the functions applying the gravitational
 *              acceleration to the particles (both gas and gravity only).
 *              The functions
 *              find_gravity_timesteps_and_do_gravity_step_first_half and
 *              do_gravity_step_second_half are directly called in the main
 *              time-evolution loop in run.c.
 *              contains functions:
 *                static inline void kick_particle(int i, double dt_gravkick,
 *                  MySingle * Grav)
 *                void find_gravity_timesteps_and_do_gravity_step_first_half(
 *                  void)
 *                void do_gravity_step_second_half(void)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 04.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../mesh/voronoi/voronoi.h"

/*! \brief Applies gravity kick to particles.
 *
 *  Apply change of velocity due to gravitational acceleration.
 *  For hydrodynamic cells, both velocity and momentum are updated.
 *
 *  \param[in] i Index of particle in P and SphP arrays.
 *  \param[in] dt_gravkick Timestep of gravity kick operation.
 *  \param[in] Grav Gravitational acceleration of particle.
 *
 *  \return void
 */
static inline void kick_particle(int i, double dt_gravkick, MySingle* Grav)
{
  int j;
  double dvel[3];
  if(P[i].Type == 0)
    {
      SphP[i].Energy -= 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
      for(j = 0; j < 3; j++) /* do the kick for gas cells */
        {
          dvel[j] = Grav[j] * dt_gravkick;
          P[i].Vel[j] += dvel[j];
          SphP[i].Momentum[j] += P[i].Mass * dvel[j];
        }
      SphP[i].Energy += 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
    }
  else
    {
      for(j = 0; j < 3; j++) /* do the kick, only collisionless particles */
        P[i].Vel[j] += Grav[j] * dt_gravkick;
    }
}

/*! \brief Performs the first half step kick operator.
 *
 *  This function applies a half step kick similar to
 *  do_gravity_step_second_half(). If we are on a PM step the kick due to
 *  the particle mesh's long range gravity is applied first. Afterwards the
 *  short range kick due to the tree force is added.
 *  In both cases the momentum and energy for gas cells is updated.
 *
 *  \return void
 */
void find_gravity_timesteps_and_do_gravity_step_first_half(void)
{
#if(defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)) && !defined(MESHRELAX)

  TIMER_START(CPU_DRIFTS);

  int idx, i;
  integertime ti_step, tstart, tend;
  double dt_gravkick;

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current) /* need to do long-range kick */
    {
      ti_step = get_timestep_pm();

      All.PM_Ti_begstep = All.PM_Ti_endstep;
      All.PM_Ti_endstep = All.PM_Ti_begstep + ti_step;

      tstart = All.PM_Ti_begstep;
      tend   = tstart + ti_step / 2;

      if(All.ComovingIntegrationOn)
        dt_gravkick = get_gravkick_factor(tstart, tend);
      else
        dt_gravkick = (tend - tstart) * All.Timebase_interval;

      for(i = 0; i < NumPart; i++)
        kick_particle(i, dt_gravkick, P[i].GravPM);
    }
#endif /* #ifdef PMGRID */

#ifdef HIERARCHICAL_GRAVITY
  /* First, move all active particles to the highest allowed timestep for this synchronization time.
   * They will then cascade down to smaller timesteps as needed.
   */

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      int i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;
      int bin    = All.HighestSynchronizedTimeBin;
      int binold = P[i].TimeBinGrav;

      timebin_move_particle(&TimeBinsGravity, i, binold, bin);
      P[i].TimeBinGrav = bin;
    }

  long long Previous_GlobalNActiveGravity = TimeBinsGravity.GlobalNActiveParticles;

  double dt_gravsum = 0;

  int bin_highest_occupied = 0;
  int timebin;
  /* go over all timebins */

  for(timebin = All.HighestSynchronizedTimeBin; timebin >= 0; timebin--)
    {
      TimeBinsGravity.NActiveParticles = 0;
      timebin_add_particles_of_timebin_to_list_of_active_particles(&TimeBinsGravity, timebin);
      sumup_large_ints(1, &TimeBinsGravity.NActiveParticles, &TimeBinsGravity.GlobalNActiveParticles);

      if(TimeBinsGravity.GlobalNActiveParticles == 0) /* we are done at this point */
        break;

      /* calculate gravity for all active particles */
      if(TimeBinsGravity.GlobalNActiveParticles != Previous_GlobalNActiveGravity)
        {
          TIMER_STOP(CPU_DRIFTS);

          compute_grav_accelerations(timebin, FLAG_PARTIAL_TREE);

          TIMER_START(CPU_DRIFTS);
        }

      int nfine = 0;
      for(int i = 0; i < TimeBinsGravity.NActiveParticles; i++)
        {
          int target = TimeBinsGravity.ActiveParticleList[i];
          int binold = P[target].TimeBinGrav;

          if(test_if_grav_timestep_is_too_large(target, binold))
            nfine++;
        }

      long long nfine_tot;
      sumup_large_ints(1, &nfine, &nfine_tot);

      int push_down_flag = 0;
      if(nfine_tot > 0.33 * TimeBinsGravity.GlobalNActiveParticles)
        push_down_flag = 1;

      for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          int i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0)
            continue;
          int binold = P[i].TimeBinGrav;

          if(push_down_flag || test_if_grav_timestep_is_too_large(i, binold))
            {
              int bin = binold - 1;
              if(bin == 0)
                {
                  print_particle_info(i);
                  terminate("timestep too small");
                }

              timebin_move_particle(&TimeBinsGravity, i, binold, bin);
              P[i].TimeBinGrav = bin;
            }
          else if(binold > bin_highest_occupied)
            bin_highest_occupied = binold;
        }

      if(All.HighestOccupiedTimeBin == 0)
        {
          MPI_Allreduce(&bin_highest_occupied, &All.HighestOccupiedTimeBin, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

          if(All.HighestOccupiedTimeBin > 0)
            {
              mpi_printf("KICKS: Special Start-up Fix: All.HighestOccupiedGravTimeBin=%d\n", All.HighestOccupiedTimeBin);

              for(i = 0; i < GRAVCOSTLEVELS; i++)
                {
                  if(All.LevelToTimeBin[i] == 0)
                    All.LevelToTimeBin[i] = All.HighestOccupiedTimeBin;
                }
            }
        }

      if(TimeBinsGravity.GlobalNActiveParticles)
        {
          ti_step = timebin ? (((integertime)1) << timebin) : 0;
          tstart  = All.Ti_begstep[timebin]; /* beginning of step */
          tend    = tstart + ti_step / 2;    /* midpoint of step */

          if(All.ComovingIntegrationOn)
            dt_gravkick = get_gravkick_factor(tstart, tend);
          else
            dt_gravkick = (tend - tstart) * All.Timebase_interval;

          if(timebin < All.HighestSynchronizedTimeBin)
            {
              ti_step = (timebin + 1) ? (((integertime)1) << (timebin + 1)) : 0;

              tstart = All.Ti_begstep[timebin + 1]; /* beginning of step */
              tend   = tstart + ti_step / 2;        /* midpoint of step */

              if(All.ComovingIntegrationOn)
                dt_gravkick -= get_gravkick_factor(tstart, tend);
              else
                dt_gravkick -= (tend - tstart) * All.Timebase_interval;
            }

          dt_gravsum += dt_gravkick;

          mpi_printf("KICKS: 1st gravity for hierarchical timebin=%d:  %lld particles\n", timebin,
                     TimeBinsGravity.GlobalNActiveParticles);

          for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
            {
              int i = TimeBinsGravity.ActiveParticleList[idx];
              if(i < 0)
                continue;

              kick_particle(i, dt_gravkick, P[i].GravAccel);
            }
          Previous_GlobalNActiveGravity = TimeBinsGravity.GlobalNActiveParticles;
        }
    }

  /* reconstruct list of active particles because it is used for other things too (i.e. wind particles) */
  timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, All.HighestActiveTimeBin);
  sumup_large_ints(1, &TimeBinsGravity.NActiveParticles, &TimeBinsGravity.GlobalNActiveParticles);
#else /* #ifdef HIERARCHICAL_GRAVITY */

#ifdef FORCE_EQUAL_TIMESTEPS
  // gravity timebin is already set, and not anymore 0 as All.HighestActiveTimeBin, but all particles should receive a first half kick
  // in the 0-th timestep
  if(All.NumCurrentTiStep == 0)
    timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, TIMEBINS);
  else
#endif /* #ifdef FORCE_EQUAL_TIMESTEPS */
    timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, All.HighestActiveTimeBin);
  sumup_large_ints(1, &TimeBinsGravity.NActiveParticles, &TimeBinsGravity.GlobalNActiveParticles);

  mpi_printf("KICKS: 1st gravity for highest active timebin=%d:  particles %lld\n", All.HighestActiveTimeBin,
             TimeBinsGravity.GlobalNActiveParticles);

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifndef FORCE_EQUAL_TIMESTEPS
      int binold = P[i].TimeBinGrav;
      int bin    = -1;

      ti_step = get_timestep_gravity(i);
      timebins_get_bin_and_do_validity_checks(ti_step, &bin, P[i].TimeBinGrav);

      if(P[i].Type == 0)
        {
          int bin_hydro = P[i].TimeBinHydro;
          if(bin_hydro < bin)
            bin = bin_hydro;
        }

      ti_step = bin ? (((integertime)1) << bin) : 0;

      timebin_move_particle(&TimeBinsGravity, i, binold, bin);
      P[i].TimeBinGrav = bin;
#else  /* #ifndef FORCE_EQUAL_TIMESTEPS */
      int bin = P[i].TimeBinGrav;
      ti_step = bin ? (((integertime)1) << bin) : 0;
#endif /* #ifndef FORCE_EQUAL_TIMESTEPS #else */

      tstart = All.Ti_begstep[bin];  /* beginning of step */
      tend   = tstart + ti_step / 2; /* midpoint of step */

      if(All.ComovingIntegrationOn)
        dt_gravkick = get_gravkick_factor(tstart, tend);
      else
        dt_gravkick = (tend - tstart) * All.Timebase_interval;

      kick_particle(i, dt_gravkick, P[i].GravAccel);
    }
#endif /* #ifdef HIERARCHICAL_GRAVITY #else */

  TIMER_STOP(CPU_DRIFTS);
#endif
}

/*! \brief Performs the second half step kick operator.
 *
 * This function applies a half step kick similar to
 * do_gravity_step_first_half(). First the short range kick due to the tree
 * force is added. If we are on a PM step the kick due to the particle mesh's
 * long range gravity is applied too. In both cases the momentum and energy
 * for gas cells is updated.
 */
void do_gravity_step_second_half(void)
{
#if(defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)) && !defined(MESHRELAX)
  TIMER_START(CPU_DRIFTS);
  int idx;
  char fullmark[8];

  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    sprintf(fullmark, "(*)");
  else
    fullmark[0] = 0;

  if(ThisTask == 0)
    fprintf(FdTimings, "\nStep%s: %d, t: %g, dt: %g, highest active timebin: %d  (lowest active: %d, highest occupied: %d)\n",
            fullmark, All.NumCurrentTiStep, All.Time, All.TimeStep, All.HighestActiveTimeBin, All.LowestActiveTimeBin,
            All.HighestOccupiedTimeBin);

  double dt_gravkick;
#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current) /* need to do long-range kick */
    {
      TIMER_STOP(CPU_DRIFTS);
      long_range_force();
      TIMER_START(CPU_DRIFTS);
    }
#endif /* #ifdef PMGRID */
#ifdef HIERARCHICAL_GRAVITY
  /* go over all timebins, in inverse sequence so that we end up getting the cumulative force at the end */
  for(int timebin = 0; timebin <= All.HighestActiveTimeBin; timebin++)
    {
      if(TimeBinSynchronized[timebin])
        {
          /* need to make all timebins below the current one active */
          timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, timebin);
          sumup_large_ints(1, &TimeBinsGravity.NActiveParticles, &TimeBinsGravity.GlobalNActiveParticles);

          if(TimeBinsGravity.GlobalNActiveParticles)
            {
              TIMER_STOP(CPU_DRIFTS);

              compute_grav_accelerations(timebin, (timebin == All.HighestActiveTimeBin) ? FLAG_FULL_TREE : FLAG_PARTIAL_TREE);

              TIMER_START(CPU_DRIFTS);

              mpi_printf("KICKS: 2nd gravity for hierarchical timebin=%d:  particles %lld\n", timebin,
                         TimeBinsGravity.GlobalNActiveParticles);

              integertime ti_step = timebin ? (((integertime)1) << timebin) : 0;

              integertime tend = All.Ti_begstep[timebin]; /* end of step (Note: All.Ti_begstep[] has already been advanced for the next
                                                             step at this point)   */
              integertime tstart = tend - ti_step / 2;    /* midpoint of step */

              if(All.ComovingIntegrationOn)
                dt_gravkick = get_gravkick_factor(tstart, tend);
              else
                dt_gravkick = (tend - tstart) * All.Timebase_interval;

              if(timebin < All.HighestActiveTimeBin)
                {
                  ti_step = (timebin + 1) ? (((integertime)1) << (timebin + 1)) : 0;

                  tend = All.Ti_begstep[timebin + 1]; /* end of step (Note: All.Ti_begstep[] has already been advanced for the next
                                                         step at this point)   */
                  tstart = tend - ti_step / 2;        /* midpoint of step */

                  if(All.ComovingIntegrationOn)
                    dt_gravkick -= get_gravkick_factor(tstart, tend);
                  else
                    dt_gravkick -= (tend - tstart) * All.Timebase_interval;
                }

              for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
                {
                  int i = TimeBinsGravity.ActiveParticleList[idx];
                  if(i < 0)
                    continue;

                  kick_particle(i, dt_gravkick, P[i].GravAccel);

                  if(P[i].Type == 0)
                    {
                      if(All.HighestOccupiedTimeBin == timebin)
                        for(int j = 0; j < 3; j++)
                          SphP[i].FullGravAccel[j] = P[i].GravAccel[j];
                    }
                }
            }
        }
    }

#else  /* #ifdef HIERARCHICAL_GRAVITY */
  timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, All.HighestActiveTimeBin);
  sumup_large_ints(1, &TimeBinsGravity.NActiveParticles, &TimeBinsGravity.GlobalNActiveParticles);

  if(TimeBinsGravity.GlobalNActiveParticles)
    {
      TIMER_STOP(CPU_DRIFTS);

      /* calculate gravity for all active particles */
      compute_grav_accelerations(All.HighestActiveTimeBin, FLAG_FULL_TREE);

      TIMER_START(CPU_DRIFTS);

      mpi_printf("KICKS: 2nd gravity for highest active timebin=%d:  particles %lld\n", All.HighestActiveTimeBin,
                 TimeBinsGravity.GlobalNActiveParticles);

      for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          int i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0)
            continue;

          integertime ti_step = P[i].TimeBinGrav ? (((integertime)1) << P[i].TimeBinGrav) : 0;
          integertime tend    = All.Ti_begstep[P[i].TimeBinGrav];
          integertime tstart  = tend - ti_step / 2; /* midpoint of step */

          if(All.ComovingIntegrationOn)
            dt_gravkick = get_gravkick_factor(tstart, tend);
          else
            dt_gravkick = (tend - tstart) * All.Timebase_interval;

          kick_particle(i, dt_gravkick, P[i].GravAccel);
        }
    }
#endif /* #ifdef HIERARCHICAL_GRAVITY #else */

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current) /* need to do long-range kick */
    {
      integertime ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;
      integertime tstart  = All.PM_Ti_begstep + ti_step / 2;
      integertime tend    = tstart + ti_step / 2;

      if(All.ComovingIntegrationOn)
        dt_gravkick = get_gravkick_factor(tstart, tend);
      else
        dt_gravkick = (tend - tstart) * All.Timebase_interval;

      for(int i = 0; i < NumPart; i++)
        kick_particle(i, dt_gravkick, P[i].GravPM);
    }
#endif /* #ifdef PMGRID */

  TIMER_STOP(CPU_DRIFTS);
#endif /* #if (defined(SELFGRAVITY) || defined(EXTERNALGRAVITY)|| defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)) && !defined(MESHRELAX) */
}
