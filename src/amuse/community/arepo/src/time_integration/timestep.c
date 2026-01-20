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
 * \file        src/time_integration/timestep.c
 * \date        05/2018
 * \brief       Routines for 'kicking' particles in
 *              momentum space and assigning new timesteps.
 * \details     contains functions:
 *                void set_cosmo_factors_for_current_time(void)
 *                void find_timesteps_without_gravity(void)
 *                void update_timesteps_from_gravity(void)
 *                integertime get_timestep_pm(void)
 *                integertime get_timestep_gravity(int p)
 *                integertime get_timestep_hydro(int p)
 *                void validate_timestep(double dt, integertime ti_step, int p)
 *                int test_if_grav_timestep_is_too_large(int p, int bin)
 *                void find_long_range_step_constraint(void)
 *                int get_timestep_bin(integertime ti_step)
 *                double get_time_difference_in_Gyr(double a0, double a1)
 *                void timebins_init(struct TimeBinData *tbData, const char
 *                  *name, int *MaxPart)
 *                void timebins_allocate(struct TimeBinData *tbData)
 *                void timebins_reallocate(struct TimeBinData *tbData)
 *                void timebins_get_bin_and_do_validity_checks(integertime
 *                  ti_step, int *bin_new, int bin_old)
 *                void timebin_move_particle(struct TimeBinData *tbData, int p,
 *                  int timeBin_old, int timeBin_new)
 *                void timebin_remove_particle(struct TimeBinData *tbData,
 *                  int idx, int bin)
 *                void timebin_add_particle(struct TimeBinData *tbData, int
 *                  i_new, int i_old, int timeBin, int
 *                  addToListOfActiveParticles)
 *                void timebin_cleanup_list_of_active_particles(struct
 *                  TimeBinData *tbData)
 *                void timebin_move_sfr(int p, int timeBin_old, int
 *                  timeBin_new)
 *                void timebin_make_list_of_active_particles_up_to_timebin(
 *                  struct TimeBinData *tbData, int timebin)
 *                void timebin_add_particles_of_timebin_to_list_of_active_
 *                  particles(struct TimeBinData *tbData, int timebin)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 11.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/*! \brief Sets various cosmological factors for the current simulation time.
 *
 *  \return void
 */
void set_cosmo_factors_for_current_time(void)
{
  if(All.ComovingIntegrationOn)
    {
      All.cf_atime    = All.Time;
      All.cf_a2inv    = 1 / (All.Time * All.Time);
      All.cf_a3inv    = 1 / (All.Time * All.Time * All.Time);
      All.cf_afac1    = pow(All.Time, 3 * GAMMA_MINUS1);
      All.cf_afac2    = 1 / pow(All.Time, 3 * GAMMA - 2);
      All.cf_afac3    = pow(All.Time, 3 * (1 - GAMMA) / 2.0);
      All.cf_hubble_a = All.cf_H = All.cf_Hrate = hubble_function(All.Time);
      All.cf_time_hubble_a                      = All.Time * All.cf_hubble_a;
      All.cf_redshift                           = 1 / All.Time - 1;
    }
  else
    {
      All.cf_atime         = 1;
      All.cf_a2inv         = 1;
      All.cf_a3inv         = 1;
      All.cf_afac1         = 1;
      All.cf_afac2         = 1;
      All.cf_afac3         = 1;
      All.cf_hubble_a      = 1;
      All.cf_H             = All.Hubble;
      All.cf_time_hubble_a = 1;
      All.cf_Hrate         = 0;
      All.cf_redshift      = 0;
    }
}

/*! \brief Finds hydrodynamic timesteps for all particles.
 *
 *  Validates the timestep and moves particles to appropriate timebin/ linked
 *  list of particles.
 *
 *  \return void
 */
void find_timesteps_without_gravity(void)
{
#ifdef TREE_BASED_TIMESTEPS
  tree_based_timesteps();
#endif /* #ifdef TREE_BASED_TIMESTEPS */

  TIMER_START(CPU_TIMELINE);

  int idx, i, bin, binold;
  integertime ti_step;

#ifdef FORCE_EQUAL_TIMESTEPS
  integertime globTimeStep = TIMEBASE;

#ifdef PMGRID
  globTimeStep = get_timestep_pm();
#endif /* #ifdef PMGRID */

#if(defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)) && !defined(MESHRELAX)
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      ti_step = get_timestep_gravity(i);
      if(ti_step < globTimeStep)
        globTimeStep = ti_step;
    }
#endif /* #if (defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)) && !defined(MESHRELAX) \
        */

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      ti_step = get_timestep_hydro(i);
      if(ti_step < globTimeStep)
        globTimeStep = ti_step;
    }

#ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME
  minimum_large_ints(1, &globTimeStep, &All.GlobalTimeStep);
#else  /* #ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME */
  MPI_Allreduce(&globTimeStep, &All.GlobalTimeStep, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif /* #ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME #else */

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      timebins_get_bin_and_do_validity_checks(All.GlobalTimeStep, &bin, P[i].TimeBinGrav);
      binold = P[i].TimeBinGrav;
      timebin_move_particle(&TimeBinsGravity, i, binold, bin);
      P[i].TimeBinGrav = bin;
    }

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      timebins_get_bin_and_do_validity_checks(All.GlobalTimeStep, &bin, P[i].TimeBinHydro);
      binold = P[i].TimeBinHydro;
      timebin_move_particle(&TimeBinsHydro, i, binold, bin);
      P[i].TimeBinHydro = bin;
    }

#else  /* #ifdef FORCE_EQUAL_TIMESTEPS */
  /* Calculate and assign hydro timesteps */

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];

      if(i < 0)
        continue;

      ti_step = get_timestep_hydro(i);

      binold = P[i].TimeBinHydro;

      timebins_get_bin_and_do_validity_checks(ti_step, &bin, binold);

      timebin_move_particle(&TimeBinsHydro, i, binold, bin);

      P[i].TimeBinHydro = bin;
    }
#endif /* #ifdef FORCE_EQUAL_TIMESTEPS #else */

  TIMER_STOP(CPU_TIMELINE);
}

/*! \brief Moves particles to lower timestep bin if required by gravity
 *         timestep criterion.
 *
 *  \return void
 */
void update_timesteps_from_gravity(void)
{
#ifdef FORCE_EQUAL_TIMESTEPS
  return; /* don't need to do this */
#endif    /* #ifdef FORCE_EQUAL_TIMESTEPS */

#if !((defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE))) || defined(MESHRELAX)
  return;
#endif /* #if !((defined(SELFGRAVITY) || defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE))) || defined(MESHRELAX) \
        */

  TIMER_START(CPU_TIMELINE);

  int idx, i, binold;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].TimeBinGrav < P[i].TimeBinHydro)
        {
          binold = P[i].TimeBinHydro;
          timebin_move_particle(&TimeBinsHydro, i, binold, P[i].TimeBinGrav);
          P[i].TimeBinHydro = P[i].TimeBinGrav;
        }
    }

  TIMER_STOP(CPU_TIMELINE);
}

#ifdef PMGRID
/*! \brief Returns particle-mesh timestep as an integer-time variable.
 *
 *  \return Integer timestep of particle-mesh algorithm.
 */
integertime get_timestep_pm(void)
{
  integertime ti_step = TIMEBASE;
  while(ti_step > (All.DtDisplacement / All.Timebase_interval))
    ti_step >>= 1;

  if(ti_step > (All.PM_Ti_endstep - All.PM_Ti_begstep)) /* PM-timestep wants to increase */
    {
      int bin    = get_timestep_bin(ti_step);
      int binold = get_timestep_bin(All.PM_Ti_endstep - All.PM_Ti_begstep);

      while(TimeBinSynchronized[bin] == 0 && bin > binold) /* make sure the new step is synchronized */
        bin--;

      ti_step = bin ? (((integertime)1) << bin) : 0;
    }

  if(All.Ti_Current == TIMEBASE) /* we here finish the last timestep. */
    ti_step = 0;

  return ti_step;
}
#endif /* #ifdef PMGRID */

/*! \brief Returns gravity timestep as an integer-time variable.
 *
 *  \param[in] p Index of particle in P array.
 *
 *  \return Integer timestep limited due to gravitational acceleration.
 */
integertime get_timestep_gravity(int p)
{
  double dt;
  integertime ti_step;

  double ax, ay, az, ac;
  {
    /* calculate total acceleration */
    ax = All.cf_a2inv * P[p].GravAccel[0];
    ay = All.cf_a2inv * P[p].GravAccel[1];
    az = All.cf_a2inv * P[p].GravAccel[2];

#if defined(PMGRID) && !defined(NO_PMFORCE_IN_SHORT_RANGE_TIMESTEP)
    ax += All.cf_a2inv * P[p].GravPM[0];
    ay += All.cf_a2inv * P[p].GravPM[1];
    az += All.cf_a2inv * P[p].GravPM[2];
#endif /* #if defined(PMGRID) && !defined(NO_PMFORCE_IN_SHORT_RANGE_TIMESTEP) */

    ac = sqrt(ax * ax + ay * ay + az * az); /* this is now the physical acceleration */

    if(ac == 0)
      ac = 1.0e-30;

    switch(All.TypeOfTimestepCriterion)
      {
        case 0:
          /* only type 0 implemented at the moment -> remove type ? */
          dt = sqrt(2 * All.ErrTolIntAccuracy * All.cf_atime * All.ForceSoftening[P[p].SofteningType] / 2.8 / ac);
          break;
        default:
          terminate("Undefined timestep criterion");
          break;
      }

#ifdef EXTERNALGRAVITY
    double dt_ext = sqrt(All.ErrTolIntAccuracy / P[p].dGravAccel);
    if(dt_ext < dt)
      dt = dt_ext;
#endif
  }

  dt *= All.cf_hubble_a;

  if(P[p].Mass == 0 && P[p].ID == 0)
    dt = All.MaxSizeTimestep; /* this particle has been swallowed or eliminated */

  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;

  if(dt < All.MinSizeTimestep)
    {
#ifdef NOSTOP_WHEN_BELOW_MINTIMESTEP
      dt = All.MinSizeTimestep;
#else  /* #ifdef NOSTOP_WHEN_BELOW_MINTIMESTEP */
      print_particle_info(p);
      terminate("Timestep dt=%g below All.MinSizeTimestep=%g", dt, All.MinSizeTimestep);
#endif /* #ifdef NOSTOP_WHEN_BELOW_MINTIMESTEP #else */
    }

#ifdef PMGRID
  if(dt >= All.DtDisplacement)
    dt = All.DtDisplacement;
#endif /* #ifdef PMGRID */

  ti_step = (integertime)(dt / All.Timebase_interval);

  validate_timestep(dt, ti_step, p);

  return ti_step;
}

/*! \brief Returns hydrodynamics timestep as an integer-time variable.
 *
 *  \param[in] p Index of particle in P and SphP array.
 *
 *  \return Integer timestep limited due to CFL condition.
 */
integertime get_timestep_hydro(int p)
{
  double dt = 0, dt_courant = 0;
  integertime ti_step;

  assert(P[p].Type == 0);

  double csnd = get_sound_speed(p);

#if defined(VORONOI_STATIC_MESH)
  csnd += sqrt(P[p].Vel[0] * P[p].Vel[0] + P[p].Vel[1] * P[p].Vel[1] + P[p].Vel[2] * P[p].Vel[2]) / All.cf_atime;
#endif /* #if defined(VORONOI_STATIC_MESH) */

  double rad = get_cell_radius(p);

  if(csnd <= 0)
    csnd = 1.0e-30;

  dt_courant = rad / csnd;

#ifdef TREE_BASED_TIMESTEPS
  if(dt_courant > SphP[p].CurrentMaxTiStep)
    dt_courant = SphP[p].CurrentMaxTiStep;
#endif /* #ifdef TREE_BASED_TIMESTEPS */

  dt_courant *= All.CourantFac;

  if(All.ComovingIntegrationOn)
    dt_courant *= All.Time;

  dt = dt_courant;

#if defined(USE_SFR)

  if(P[p].Type == 0) /* to protect using a particle that has been turned into a star */
    {
      double sfr = get_starformation_rate(p);

      double dt_sfr = 0.1 * P[p].Mass / (sfr / ((All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR)));
      if(dt_sfr < dt)
        dt = dt_sfr;
    }
#endif /* #if defined(USE_SFR) */

#ifdef MHD_POWELL_LIMIT_TIMESTEP
  double b         = sqrt(SphP[p].B[0] * SphP[p].B[0] + SphP[p].B[1] * SphP[p].B[1] + SphP[p].B[2] * SphP[p].B[2]);
  double bmin      = sqrt(2 * 0.01 * SphP[p].Utherm * SphP[p].Density * All.cf_atime);
  double v         = sqrt(P[p].Vel[0] * P[p].Vel[0] + P[p].Vel[1] * P[p].Vel[1] + P[p].Vel[2] * P[p].Vel[2]) / All.cf_atime;
  double dt_powell = 0.5 * (b + bmin) / (fabs(SphP[p].DivB / All.cf_atime * v));

  if(dt_powell < dt)
    dt = dt_powell;
#endif /* #ifdef MHD_POWELL_LIMIT_TIMESTEP */

  /* convert the physical timestep to dloga if needed. Note: If comoving integration has not been selected,
     All.cf_hubble_a=1.
   */

  dt *= All.cf_hubble_a;

  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;

#ifdef TIMESTEP_OUTPUT_LIMIT
  if(dt >= All.TimestepOutputLimit)
    dt = All.TimestepOutputLimit;
#endif /* #ifdef TIMESTEP_OUTPUT_LIMIT */

  if(dt < All.MinSizeTimestep)
    {
#ifdef NOSTOP_WHEN_BELOW_MINTIMESTEP
      dt = All.MinSizeTimestep;
#else  /* #ifdef NOSTOP_WHEN_BELOW_MINTIMESTEP */
      print_particle_info(p);
      terminate("Timestep dt=%g below All.MinSizeTimestep=%g", dt, All.MinSizeTimestep);
#endif /* #ifdef NOSTOP_WHEN_BELOW_MINTIMESTEP #else */
    }

#ifdef PMGRID
  if(dt >= All.DtDisplacement)
    dt = All.DtDisplacement;
#endif /* #ifdef PMGRID */

  ti_step = (integertime)(dt / All.Timebase_interval);

  validate_timestep(dt, ti_step, p);

  return ti_step;
}

/*! \brief Checks if timestep is a valid one.
 *
 *  Terminates the simulation with error message otherwise.
 *
 *  \return void
 */
void validate_timestep(double dt, integertime ti_step, int p)
{
  if(!(ti_step > 0 && ti_step < TIMEBASE))
    {
      printf(
          "\nError: An invalid timestep was assigned on the integer timeline!\n"
          "We better stop.\n"
          "Task=%d Part-ID=%lld type=%d",
          ThisTask, (long long)P[p].ID, P[p].Type);

      printf("tibase=%g dt=%g ti_step=%d, xyz=(%g|%g|%g) vel=(%g|%g|%g) tree=(%g|%g|%g) mass=%g\n\n", All.Timebase_interval, dt,
             (int)ti_step, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], P[p].Vel[0], P[p].Vel[1], P[p].Vel[2], P[p].GravAccel[0],
             P[p].GravAccel[1], P[p].GravAccel[2], P[p].Mass);

      print_particle_info(p);
      myflush(stdout);
      terminate("integer timestep outside of allowed range");
    }

  if(ti_step == 1)
    {
      printf("Time-step of integer size 1 found for particle i=%d, pos=(%g|%g|%g), ID=%lld, dt=%g\n", p, P[p].Pos[0], P[p].Pos[1],
             P[p].Pos[2], (long long)P[p].ID, dt);
      print_particle_info(p);
    }
}

/*! \brief Checks if timestep according to its present timebin is too large
 *         compared to the requirements from gravity and hydrodynamics
 *
 *  I.e. does the cell need to be moved to a finer timebin?
 *
 *  \param[in] p Index of particle/cell.
 *  \param[in] bin Timebin to compare to.
 *
 *  \return 0: not too large; 1: too large.
 */
int test_if_grav_timestep_is_too_large(int p, int bin)
{
  integertime ti_step_bin = bin ? (((integertime)1) << bin) : 0;

  integertime ti_step = get_timestep_gravity(p);

  if(P[p].Type == 0)
    {
      if((P[p].ID != 0) && (P[p].Mass != 0))
        {
          int bin_hydro             = P[p].TimeBinHydro;
          integertime ti_step_hydro = bin_hydro ? (((integertime)1) << bin_hydro) : 0;
          if(ti_step_hydro < ti_step)
            ti_step = ti_step_hydro;
        }
    }

  if(ti_step < ti_step_bin)
    return 1;
  else
    return 0;
}

#ifdef PMGRID
/*! \brief Sets the global timestep for the long-range force calculation.
 *
 *  Evaluates timestep constraints due to long range force acceleration of all
 *  simulation particles and finds its global minimum.
 *
 *  \return void
 */
void find_long_range_step_constraint(void)
{
  int p;
  double ax, ay, az, ac;
  double dt, dtmin = MAX_DOUBLE_NUMBER;

  for(p = 0; p < NumPart; p++)
    {
      if(P[p].Type == 0)
        continue;

#ifdef PM_TIMESTEP_BASED_ON_TYPES
      if(((1 << P[p].Type) & (PM_TIMESTEP_BASED_ON_TYPES)))
#endif /* #ifdef PM_TIMESTEP_BASED_ON_TYPES */
        {
          /* calculate acceleration */
          ax = All.cf_a2inv * P[p].GravPM[0];
          ay = All.cf_a2inv * P[p].GravPM[1];
          az = All.cf_a2inv * P[p].GravPM[2];

          ac = sqrt(ax * ax + ay * ay + az * az); /* this is now the physical acceleration */

          if(ac < MIN_FLOAT_NUMBER)
            ac = MIN_FLOAT_NUMBER;

          dt = sqrt(2.0 * All.ErrTolIntAccuracy * All.cf_atime * All.ForceSoftening[P[p].SofteningType] / (2.8 * ac));

          dt *= All.cf_hubble_a;

          if(dt < dtmin)
            dtmin = dt;
        }
    }

  dtmin *= 2.0; /* move it one timebin higher to prevent being too conservative */

  MPI_Allreduce(&dtmin, &All.DtDisplacement, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  mpi_printf("TIMESTEPS: displacement time constraint: %g  (%g)\n", All.DtDisplacement, All.MaxSizeTimestep);

  if(All.DtDisplacement > All.MaxSizeTimestep)
    All.DtDisplacement = All.MaxSizeTimestep;

  if(All.DtDisplacement < All.MinSizeTimestep)
    All.DtDisplacement = All.MinSizeTimestep;
}
#endif /* #ifdef PMGRID */

/*! \brief Converts an integer time to a time bin.
 *
 *  \param[in] ti_step Timestep as integertime variable.
 *
 *  \return Associated time-bin.
 */
int get_timestep_bin(integertime ti_step)
{
  int bin = -1;

  if(ti_step == 0)
    return 0;

  if(ti_step == 1)
    terminate("time-step of integer size 1 not allowed\n");

  while(ti_step)
    {
      bin++;
      ti_step >>= 1;
    }

  return bin;
}

/*! \brief Calculates time difference in Gyr between two time integration unit
 *         values.
 *
 *  If simulation non-cosmological, a0 and a1 are proper time in code units,
 *  for cosmological simulation a0 and a1 are scalefactors.
 *
 *  \param[in] a0 First time or scalefactor.
 *  \param[in] a1 Second time or scalefactor.
 *
 *  \return Time difference in Gyr.
 */
double get_time_difference_in_Gyr(double a0, double a1)
{
  double result, time_diff = 0, t0, t1, factor1, factor2, term1, term2;

  if(All.ComovingIntegrationOn)
    {
      if(All.OmegaLambda + All.Omega0 != 1)
        printf("only implemented for flat cosmology so far.");

      factor1 = 2.0 / (3.0 * sqrt(All.OmegaLambda));

      term1   = sqrt(All.OmegaLambda / All.Omega0) * pow(a0, 1.5);
      term2   = sqrt(1 + All.OmegaLambda / All.Omega0 * pow(a0, 3));
      factor2 = log(term1 + term2);

      t0 = factor1 * factor2;

      term1   = sqrt(All.OmegaLambda / All.Omega0) * pow(a1, 1.5);
      term2   = sqrt(1 + All.OmegaLambda / All.Omega0 * pow(a1, 3));
      factor2 = log(term1 + term2);

      t1 = factor1 * factor2;

      result = t1 - t0;

      time_diff = result / (HUBBLE * All.HubbleParam); /* now in seconds */
      time_diff /= SEC_PER_MEGAYEAR * 1000;            /* now in gigayears */
    }
  else
    {
      time_diff = (a1 - a0) * All.UnitTime_in_s / All.HubbleParam; /* now in seconds */
      time_diff /= SEC_PER_MEGAYEAR * 1000;                        /* now in gigayears */
    }

  return time_diff;
}

/*! \brief Initializes time bin data.
 *
 *  Does not allocate anything!
 *
 *  \param[out] tbData Time bin data to be initialized.
 *  \param[in] name Name stored in time bin data.
 *  \param[in] MaxPart Maximum number of particles in time bin data.
 *
 *  \return void
 */
void timebins_init(struct TimeBinData *tbData, const char *name, int *MaxPart)
{
  int i;
  tbData->NActiveParticles   = 0;
  tbData->ActiveParticleList = 0;

  for(i = 0; i < TIMEBINS; i++)
    {
      tbData->FirstInTimeBin[i] = -1;
      tbData->LastInTimeBin[i]  = -1;
    }

  tbData->NextInTimeBin = 0;
  tbData->PrevInTimeBin = 0;

  strncpy(tbData->Name, name, 99);
  tbData->Name[99] = 0;
  tbData->MaxPart  = MaxPart;
}

/*! \brief Allocates linked lists in time bin data.
 *
 *  With tbData->MaxPart elements.
 *
 *  \param[in, out] tbData Pointer to time bin data to be allocated.
 *
 *  \return void
 */
void timebins_allocate(struct TimeBinData *tbData)
{
  char Identifier[200];
  Identifier[199] = 0;

  snprintf(Identifier, 199, "NextActiveParticle%s", tbData->Name);
  tbData->ActiveParticleList = (int *)mymalloc_movable(&tbData->ActiveParticleList, Identifier, *(tbData->MaxPart) * sizeof(int));

  snprintf(Identifier, 199, "NextInTimeBin%s", tbData->Name);
  tbData->NextInTimeBin = (int *)mymalloc_movable(&tbData->NextInTimeBin, Identifier, *(tbData->MaxPart) * sizeof(int));

  snprintf(Identifier, 199, "PrevInTimeBin%s", tbData->Name);
  tbData->PrevInTimeBin = (int *)mymalloc_movable(&tbData->PrevInTimeBin, Identifier, *(tbData->MaxPart) * sizeof(int));
}

/*! \brief Re-allocates linked lists in time bin data.
 *
 *  With tbData->MaxPart elements.
 *
 *  \param[out] tbData Pointer to time bin data to be re-allocated.
 *
 *  \return void
 */
void timebins_reallocate(struct TimeBinData *tbData)
{
  tbData->ActiveParticleList = (int *)myrealloc_movable(tbData->ActiveParticleList, *(tbData->MaxPart) * sizeof(int));
  tbData->NextInTimeBin      = (int *)myrealloc_movable(tbData->NextInTimeBin, *(tbData->MaxPart) * sizeof(int));
  tbData->PrevInTimeBin      = (int *)myrealloc_movable(tbData->PrevInTimeBin, *(tbData->MaxPart) * sizeof(int));
}

/*! \brief Gets timebin and checks if bin is valid.
 *
 *  Checks for example if old bin is synchronized with the bin it should be
 *  moved to.
 *
 *  \param[in] ti_step Timestep in integertime.
 *  \param[out] bin_new New time bin.
 *  \param[in] bin_old Old time bin.
 *
 *  \return void
 */
void timebins_get_bin_and_do_validity_checks(integertime ti_step, int *bin_new, int bin_old)
{
  /* make it a power 2 subdivision */
  integertime ti_min = TIMEBASE;
  while(ti_min > ti_step)
    ti_min >>= 1;
  ti_step = ti_min;

  /* get timestep bin */
  int bin = -1;

  if(ti_step == 0)
    bin = 0;

  if(ti_step == 1)
    terminate("time-step of integer size 1 not allowed\n");

  while(ti_step)
    {
      bin++;
      ti_step >>= 1;
    }

  if(bin > bin_old) /* timestep wants to increase */
    {
      while(TimeBinSynchronized[bin] == 0 && bin > bin_old) /* make sure the new step is synchronized */
        bin--;

      ti_step = bin ? (((integertime)1) << bin) : 0;
    }

  if(All.Ti_Current >= TIMEBASE) /* we here finish the last timestep. */
    {
      ti_step = 0;
      bin     = 0;
    }

  if((TIMEBASE - All.Ti_Current) < ti_step) /* check that we don't run beyond the end */
    {
      terminate("we are beyond the end of the timeline"); /* should not happen */
    }

  *bin_new = bin;
}

/*! \brief Move particle from one time bin to another.
 *
 *  \param[in, out] tbData Time bin data structure to operate on.
 *  \param[in] p Index of the particle to be moved.
 *  \param[in] timeBin_old Old time bin of particle to be moved.
 *  \param[in] timeBin_new New time bin of particle to be moved.
 *
 *  \return void
 */
void timebin_move_particle(struct TimeBinData *tbData, int p, int timeBin_old, int timeBin_new)
{
  if(timeBin_old == timeBin_new)
    return;

  tbData->TimeBinCount[timeBin_old]--;

  int prev = tbData->PrevInTimeBin[p];
  int next = tbData->NextInTimeBin[p];

  if(tbData->FirstInTimeBin[timeBin_old] == p)
    tbData->FirstInTimeBin[timeBin_old] = next;
  if(tbData->LastInTimeBin[timeBin_old] == p)
    tbData->LastInTimeBin[timeBin_old] = prev;
  if(prev >= 0)
    tbData->NextInTimeBin[prev] = next;
  if(next >= 0)
    tbData->PrevInTimeBin[next] = prev;

  if(tbData->TimeBinCount[timeBin_new] > 0)
    {
      tbData->PrevInTimeBin[p]                                  = tbData->LastInTimeBin[timeBin_new];
      tbData->NextInTimeBin[tbData->LastInTimeBin[timeBin_new]] = p;
      tbData->NextInTimeBin[p]                                  = -1;
      tbData->LastInTimeBin[timeBin_new]                        = p;
    }
  else
    {
      tbData->FirstInTimeBin[timeBin_new] = tbData->LastInTimeBin[timeBin_new] = p;
      tbData->PrevInTimeBin[p] = tbData->NextInTimeBin[p] = -1;
    }

  tbData->TimeBinCount[timeBin_new]++;

#ifdef USE_SFR
  if((P[p].Type == 0) && (tbData == &TimeBinsHydro))
    timebin_move_sfr(p, timeBin_old, timeBin_new);
#endif /* #ifdef USE_SFR */
}

/*! \brief Removes a particle from time bin structure.
 *
 *  Can only be done with active particles.
 *
 *  \param[in, out] tbData Time bin structure to be operated on.
 *  \param[in] idx Index of particle in ActiveParticleList.
 *  \param[in] bin Timebin in which particle is currently. If left -1, function
 *             will determine bin by itself.
 *
 *  \return void
 */
void timebin_remove_particle(struct TimeBinData *tbData, int idx, int bin)
{
  int p                           = tbData->ActiveParticleList[idx];
  tbData->ActiveParticleList[idx] = -1;

  if(bin == -1)
    {
      if(tbData == &TimeBinsGravity)
        bin = P[p].TimeBinGrav;
      else
        bin = P[p].TimeBinHydro;
    }

  tbData->TimeBinCount[bin]--;

  if(p >= 0)
    {
      int prev = tbData->PrevInTimeBin[p];
      int next = tbData->NextInTimeBin[p];

      if(prev >= 0)
        tbData->NextInTimeBin[prev] = next;
      if(next >= 0)
        tbData->PrevInTimeBin[next] = prev;

      if(tbData->FirstInTimeBin[bin] == p)
        tbData->FirstInTimeBin[bin] = next;
      if(tbData->LastInTimeBin[bin] == p)
        tbData->LastInTimeBin[bin] = prev;
    }
}

/* \brief Inserts a particle into the timebin struct behind another already
 *        existing particle.
 *
 *  \param[in, out] tbData Time bin structure to be operated on.
 *  \param[in] i_new New index in linked lists of time bin data.
 *  \param[in] i_old old index in linked lists of time bin data.
 *  \param[in] timeBin Time bin to which it should be added.
 *  \param[in] addToListOfActiveParticles Flag if particle should be added as
 *             an active particle.
 *
 *  \return void
 */
void timebin_add_particle(struct TimeBinData *tbData, int i_new, int i_old, int timeBin, int addToListOfActiveParticles)
{
  tbData->TimeBinCount[timeBin]++;

  if(i_old < 0)
    {
      /* if we don't have an existing particle to add if after, let's take the last one in this timebin */
      i_old = tbData->LastInTimeBin[timeBin];

      if(i_old < 0)
        {
          /* the timebin is empty at the moment, so just add the new particle */
          tbData->FirstInTimeBin[timeBin] = i_new;
          tbData->LastInTimeBin[timeBin]  = i_new;
          tbData->NextInTimeBin[i_new]    = -1;
          tbData->PrevInTimeBin[i_new]    = -1;
        }
    }

  if(i_old >= 0)
    {
      /* otherwise we added it already */
      tbData->PrevInTimeBin[i_new] = i_old;
      tbData->NextInTimeBin[i_new] = tbData->NextInTimeBin[i_old];
      if(tbData->NextInTimeBin[i_old] >= 0)
        tbData->PrevInTimeBin[tbData->NextInTimeBin[i_old]] = i_new;
      tbData->NextInTimeBin[i_old] = i_new;
      if(tbData->LastInTimeBin[timeBin] == i_old)
        tbData->LastInTimeBin[timeBin] = i_new;
    }

  if(addToListOfActiveParticles)
    {
      tbData->ActiveParticleList[tbData->NActiveParticles] = i_new;
      tbData->NActiveParticles++;
    }
}

/*! \brief Removes active particles that have ID and Mass 0, i.e. that were
 *         flagged as deleted from time bin data structure.
 *
 *  \param[in, out] tbData Time bin data structure to be operated on.
 *
 *  \return void
 */
void timebin_cleanup_list_of_active_particles(struct TimeBinData *tbData)
{
  int idx, i;
  for(idx = 0; idx < tbData->NActiveParticles; idx++)
    {
      i = tbData->ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].ID == 0 && P[i].Mass == 0)
        timebin_remove_particle(tbData, idx, -1);
    }
}

#ifdef USE_SFR
/*! \brief Updates TimeBinSfr when a gas cell changes timebin.
 *
 *  \param[in] p Index of cell in SphP array.
 *  \param[in] timeBin_old Old time bin.
 *  \param[in] timeBin_new New time bin.
 *
 *  \return void
 */
void timebin_move_sfr(int p, int timeBin_old, int timeBin_new)
{
  TimeBinSfr[timeBin_old] -= SphP[p].Sfr;
  TimeBinSfr[timeBin_new] += SphP[p].Sfr;
}
#endif /* #ifdef USE_SFR */

/*! \brief Crates list of active particles up to a specified timebin.
 *
 *  \param[in, out] tbData Time bin data to be operated on.
 *  \param[in] timebin Up to which timebin should particles be included.
 *
 *  \return void
 */
void timebin_make_list_of_active_particles_up_to_timebin(struct TimeBinData *tbData, int timebin)
{
  int tbin;
  tbData->NActiveParticles = 0;
  for(tbin = timebin; tbin >= 0; tbin--)
    timebin_add_particles_of_timebin_to_list_of_active_particles(tbData, tbin);
}

/*! \brief Add particles of a specific timebin to active particle list.
 *
 *  \param[in, out] tbData Time bin data to be operated on.
 *  \param[in] timebin Time bin which should be included.
 *
 *  \return void
 */
void timebin_add_particles_of_timebin_to_list_of_active_particles(struct TimeBinData *tbData, int timebin)
{
  int i;
  for(i = tbData->FirstInTimeBin[timebin]; i >= 0; i = tbData->NextInTimeBin[i])
    if(!(P[i].ID == 0 && P[i].Mass == 0))
      {
        tbData->ActiveParticleList[tbData->NActiveParticles] = i;
        tbData->NActiveParticles++;
      }
}
