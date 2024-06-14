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
 * \file        src/time_integration/predict.c
 * \date        05/2018
 * \brief       Routines to find the next sync point, manage the list
 *              of active timebins/active particles and to drift particles.
 * \details     contains functions:
 *                void reconstruct_timebins(void)
 *                void find_next_sync_point(void)
 *                void mark_active_timebins(void)
 *                void drift_all_particles(void)
 *                void drift_particle(int i, integertime time1)
 *                static int int_compare(const void *a, const void *b)
 *                void make_list_of_active_particles(void)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 08.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/*! \brief This function (re)builds the time bin lists.
 *
 *  It counts the number of particles in each timebin and updates the
 *  linked lists containing the particles of each time bin. Afterwards the
 *  linked list of active particles is updated by
 *  make_list_of_active_particles().
 *
 *  The linked lists for each timebin are stored in 'FirstInTimeBin',
 *  'LastInTimeBin', 'PrevInTimeBin' and 'NextInTimeBin'. The counters
 *  of particles per timebin are 'TimeBinCount' and 'TimeBinCountSph'.
 *
 *  \return void
 */
void reconstruct_timebins(void)
{
  TIMER_START(CPU_TIMELINE);

  int i, bin;

  for(bin = 0; bin < TIMEBINS; bin++)
    {
      TimeBinsHydro.TimeBinCount[bin]   = 0;
      TimeBinsHydro.FirstInTimeBin[bin] = -1;
      TimeBinsHydro.LastInTimeBin[bin]  = -1;

      TimeBinsGravity.TimeBinCount[bin]   = 0;
      TimeBinsGravity.FirstInTimeBin[bin] = -1;
      TimeBinsGravity.LastInTimeBin[bin]  = -1;

#ifdef USE_SFR
      TimeBinSfr[bin] = 0;
#endif
    }

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].ID == 0 && P[i].Mass == 0)
        continue;

      if(P[i].Type != 0)
        continue;

      bin = P[i].TimeBinHydro;

      if(TimeBinsHydro.TimeBinCount[bin] > 0)
        {
          TimeBinsHydro.PrevInTimeBin[i]                                = TimeBinsHydro.LastInTimeBin[bin];
          TimeBinsHydro.NextInTimeBin[i]                                = -1;
          TimeBinsHydro.NextInTimeBin[TimeBinsHydro.LastInTimeBin[bin]] = i;
          TimeBinsHydro.LastInTimeBin[bin]                              = i;
        }
      else
        {
          TimeBinsHydro.FirstInTimeBin[bin] = TimeBinsHydro.LastInTimeBin[bin] = i;
          TimeBinsHydro.PrevInTimeBin[i] = TimeBinsHydro.NextInTimeBin[i] = -1;
        }
      TimeBinsHydro.TimeBinCount[bin]++;

#ifdef USE_SFR
      TimeBinSfr[bin] += SphP[i].Sfr;
#endif
    }

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].ID == 0 && P[i].Mass == 0)
        continue;

      bin = P[i].TimeBinGrav;

      if(TimeBinsGravity.TimeBinCount[bin] > 0)
        {
          TimeBinsGravity.PrevInTimeBin[i]                                  = TimeBinsGravity.LastInTimeBin[bin];
          TimeBinsGravity.NextInTimeBin[i]                                  = -1;
          TimeBinsGravity.NextInTimeBin[TimeBinsGravity.LastInTimeBin[bin]] = i;
          TimeBinsGravity.LastInTimeBin[bin]                                = i;
        }
      else
        {
          TimeBinsGravity.FirstInTimeBin[bin] = TimeBinsGravity.LastInTimeBin[bin] = i;
          TimeBinsGravity.PrevInTimeBin[i] = TimeBinsGravity.NextInTimeBin[i] = -1;
        }
      TimeBinsGravity.TimeBinCount[bin]++;
    }

  make_list_of_active_particles();

  TIMER_STOP(CPU_TIMELINE);
}

/*! \brief This function finds the next synchronization point of the system.
 *         (i.e. the earliest point of time any of the particles needs a force
 *         computation).
 *
 *  \return void
 */
void find_next_sync_point(void)
{
  int n;
  integertime ti_next_kick, ti_next_kick_global, ti_next_for_bin, dt_bin;
  double timeold;

  TIMER_START(CPU_DRIFTS);

  timeold = All.Time;

  All.NumCurrentTiStep++;

  /* find the next kick time */
  ti_next_kick = TIMEBASE;

  for(n = 0; n < TIMEBINS; n++)
    {
      int active = TimeBinsHydro.TimeBinCount[n];

#if(defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)) && !defined(MESHRELAX)
      active += TimeBinsGravity.TimeBinCount[n];
#endif /* #if (defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)) && !defined(MESHRELAX) \
        */
      if(active)
        {
          if(n > 0)
            {
              dt_bin          = (((integertime)1) << n);
              ti_next_for_bin = (All.Ti_Current / dt_bin) * dt_bin + dt_bin; /* next kick time for this timebin */
            }
          else
            {
              dt_bin          = 0;
              ti_next_for_bin = All.Ti_Current;
            }

          if(ti_next_for_bin < ti_next_kick)
            ti_next_kick = ti_next_for_bin;
        }
    }

#ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME
  minimum_large_ints(1, &ti_next_kick, &ti_next_kick_global);
#else  /* #ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME */
  MPI_Allreduce(&ti_next_kick, &ti_next_kick_global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif /* #ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME #else */

  All.Previous_Ti_Current = All.Ti_Current;
  All.Ti_Current          = ti_next_kick_global;

  if(All.ComovingIntegrationOn)
    All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
  else
    All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;

  set_cosmo_factors_for_current_time();

  All.TimeStep = All.Time - timeold;

  mark_active_timebins();

  TIMER_STOP(CPU_DRIFTS);
}

/*! \brief Sets active timebins for current time-step in global variables.
 *
 *  \return void
 */
void mark_active_timebins(void)
{
  int n;
  int lowest_active_bin = TIMEBINS, highest_active_bin = 0;
  int lowest_occupied_bin = TIMEBINS, highest_occupied_bin = 0;
  int lowest_occupied_gravity_bin = TIMEBINS, highest_occupied_gravity_bin = 0;
  int highest_synchronized_bin = 0;
  int nsynchronized_gravity = 0, nsynchronized_hydro = 0;
  integertime dt_bin;

  /* mark the bins that will be synchronized/active */

  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinsGravity.TimeBinCount[n])
        {
          if(highest_occupied_gravity_bin < n)
            highest_occupied_gravity_bin = n;

          if(lowest_occupied_gravity_bin > n)
            lowest_occupied_gravity_bin = n;
        }

      int active = TimeBinsHydro.TimeBinCount[n] + TimeBinsGravity.TimeBinCount[n];

      if(active)
        {
          if(highest_occupied_bin < n)
            highest_occupied_bin = n;

          if(lowest_occupied_bin > n)
            lowest_occupied_bin = n;
        }

      dt_bin = (((integertime)1) << n);

      if((All.Ti_Current % dt_bin) == 0)
        {
          TimeBinSynchronized[n] = 1;
          All.Ti_begstep[n]      = All.Ti_Current;

          nsynchronized_gravity += TimeBinsGravity.TimeBinCount[n];
          nsynchronized_hydro += TimeBinsHydro.TimeBinCount[n];

          if(highest_synchronized_bin < n)
            highest_synchronized_bin = n;

          if(active)
            {
              if(highest_active_bin < n)
                highest_active_bin = n;

              if(lowest_active_bin > n)
                lowest_active_bin = n;
            }
        }
      else
        TimeBinSynchronized[n] = 0;
    }

  int lowest_in[3], lowest_out[3];
  lowest_in[0] = lowest_occupied_bin;
  lowest_in[1] = lowest_occupied_gravity_bin;
  lowest_in[2] = lowest_active_bin;
  MPI_Allreduce(lowest_in, lowest_out, 3, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  All.LowestOccupiedTimeBin     = lowest_out[0];
  All.LowestOccupiedGravTimeBin = lowest_out[1];
  All.LowestActiveTimeBin       = lowest_out[2];

  int highest_in[4], highest_out[4];
  highest_in[0] = highest_occupied_bin;
  highest_in[1] = highest_occupied_gravity_bin;
  highest_in[2] = highest_active_bin;
  highest_in[3] = highest_synchronized_bin;
  MPI_Allreduce(highest_in, highest_out, 4, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  All.HighestOccupiedTimeBin     = highest_out[0];
  All.HighestOccupiedGravTimeBin = highest_out[1];
  All.HighestActiveTimeBin       = highest_out[2];
  All.HighestSynchronizedTimeBin = highest_out[3];

  /* note: the lowest synchronized bin is always 1 */

  int input_ints[2 + 2 * TIMEBINS];
  long long output_longs[2 + 2 * TIMEBINS];

  input_ints[0] = nsynchronized_hydro;
  input_ints[1] = nsynchronized_gravity;
  memcpy(input_ints + 2, TimeBinsGravity.TimeBinCount, TIMEBINS * sizeof(int));
  memcpy(input_ints + 2 + TIMEBINS, TimeBinsHydro.TimeBinCount, TIMEBINS * sizeof(int));

  sumup_large_ints(2 + 2 * TIMEBINS, input_ints, output_longs);

  All.GlobalNSynchronizedHydro   = output_longs[0];
  All.GlobalNSynchronizedGravity = output_longs[1];
  long long *tot_count_grav      = output_longs + 2;
  long long *tot_count_sph       = output_longs + 2 + TIMEBINS;

  long long tot_grav = 0, tot_sph = 0;

  for(n = 0; n < TIMEBINS; n++)
    {
      tot_grav += tot_count_grav[n];
      tot_sph += tot_count_sph[n];

      if(n > 0)
        {
          tot_count_grav[n] += tot_count_grav[n - 1];
          tot_count_sph[n] += tot_count_sph[n - 1];
        }
    }

  All.SmallestTimeBinWithDomainDecomposition = All.HighestOccupiedTimeBin;

  for(n = All.HighestOccupiedTimeBin; n >= All.LowestOccupiedTimeBin; n--)
    {
      if(tot_count_grav[n] > All.ActivePartFracForNewDomainDecomp * tot_grav ||
         tot_count_sph[n] > All.ActivePartFracForNewDomainDecomp * tot_sph)
        All.SmallestTimeBinWithDomainDecomposition = n;
    }
}

/*! \brief Applies drift operation to all particles to current time.
 *
 *  \return void
 */
void drift_all_particles(void)
{
  int i;

  TIMER_START(CPU_DRIFTS);

  for(i = 0; i < NumPart; i++)
    drift_particle(i, All.Ti_Current);

  TIMER_STOP(CPU_DRIFTS);
}

/*! \brief This function drifts drifts a particle i to time1.
 *
 * \param[in] i Particle/cell index.
 * \param[in] time1 Time to which particles get drifted.
 *
 * \return void
 */
void drift_particle(int i, integertime time1)
{
  int j;

  if(i < 0)
    terminate("i=%d  NumPart=%d", i, NumPart);

  integertime time0 = P[i].Ti_Current;

  if(time1 == time0)
    return;

  if(time1 < time0)
    terminate("no prediction into past allowed: time0=%lld time1=%lld\n", (long long)time0, (long long)time1);

  double dt_drift;

  if(All.ComovingIntegrationOn)
    dt_drift = get_drift_factor(time0, time1);
  else
    dt_drift = (time1 - time0) * All.Timebase_interval;

  if(P[i].Type == 0)
    {
      for(j = 0; j < 3; j++)
        {
          P[i].Pos[j] += SphP[i].VelVertex[j] * dt_drift;
        }
    }
  else
    {
#ifndef MESHRELAX
      for(j = 0; j < 3; j++)
        P[i].Pos[j] += P[i].Vel[j] * dt_drift;

#if defined(REFLECTIVE_X)
      if(P[i].Pos[0] < 0 || P[i].Pos[0] > boxSize_X)
        {
          P[i].Pos[0] = 2 * (P[i].Pos[0] > boxSize_X ? 1 : 0) * boxSize_X - P[i].Pos[0];
          P[i].Vel[0] *= -1;
        }
#endif /* #if defined(REFLECTIVE_X) */
#if defined(REFLECTIVE_Y)
      if(P[i].Pos[1] < 0 || P[i].Pos[1] > boxSize_Y)
        {
          P[i].Pos[1] = 2 * (P[i].Pos[1] > boxSize_Y ? 1 : 0) * boxSize_Y - P[i].Pos[1];
          P[i].Vel[1] *= -1;
        }
#endif /* #if defined(REFLECTIVE_Y) */
#if defined(REFLECTIVE_Z)
      if(P[i].Pos[2] < 0 || P[i].Pos[2] > boxSize_Z)
        {
          P[i].Pos[2] = 2 * (P[i].Pos[2] > boxSize_Z ? 1 : 0) * boxSize_Z - P[i].Pos[2];
          P[i].Vel[2] *= -1;
        }
#endif /* #if defined(REFLECTIVE_Z) */

#endif /* #ifndef MESHRELAX */
    }

  P[i].Ti_Current = time1;
}

/*! \brief Comparison function for two integer values.
 *
 *  \param[in] a First value.
 *  \param[in] b Second value.
 *
 *  \return (-1,0,1); -1 if a < b
 */
static int int_compare(const void *a, const void *b)
{
  if(*((int *)a) < *((int *)b))
    return -1;

  if(*((int *)a) > *((int *)b))
    return +1;

  return 0;
}

/*! \brief This function builds the linear list of active particles.
 *
 *  The list is stored in the array ActiveParticleList of the TimeBinData
 *  structs.
 *
 *  \return void
 */
void make_list_of_active_particles(void)
{
  TIMER_START(CPU_DRIFTS);

  int i, n;
  /* make a link list with the particles in the active time bins */
  TimeBinsHydro.NActiveParticles = 0;

  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinSynchronized[n])
        {
          for(i = TimeBinsHydro.FirstInTimeBin[n]; i >= 0; i = TimeBinsHydro.NextInTimeBin[i])
            if((P[i].Type == 0) && !((P[i].ID == 0) && (P[i].Mass == 0)))
              {
                if(P[i].Ti_Current != All.Ti_Current)
                  drift_particle(i, All.Ti_Current);

                TimeBinsHydro.ActiveParticleList[TimeBinsHydro.NActiveParticles] = i;
                TimeBinsHydro.NActiveParticles++;
              }
        }
    }

  TimeBinsGravity.NActiveParticles = 0;

  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinSynchronized[n])
        {
          for(i = TimeBinsGravity.FirstInTimeBin[n]; i >= 0; i = TimeBinsGravity.NextInTimeBin[i])
            {
              if(!((P[i].ID == 0) && (P[i].Mass == 0)))
                {
                  if(P[i].Ti_Current != All.Ti_Current)
                    drift_particle(i, All.Ti_Current);

                  TimeBinsGravity.ActiveParticleList[TimeBinsGravity.NActiveParticles] = i;
                  TimeBinsGravity.NActiveParticles++;
                }
            }
        }
    }

  /* sort both lists for better memory efficiency */
  mysort(TimeBinsHydro.ActiveParticleList, TimeBinsHydro.NActiveParticles, sizeof(int), int_compare);
  mysort(TimeBinsGravity.ActiveParticleList, TimeBinsGravity.NActiveParticles, sizeof(int), int_compare);

  int in[6];
  long long out[6];

  n     = 2;
  in[0] = TimeBinsGravity.NActiveParticles;
  in[1] = TimeBinsHydro.NActiveParticles;

  sumup_large_ints(n, in, out);

  TimeBinsGravity.GlobalNActiveParticles = out[0];
  TimeBinsHydro.GlobalNActiveParticles   = out[1];

  TIMER_STOP(CPU_DRIFTS);
}
