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
 * \file        src/star_formation/starformation.c
 * \date        05/2018
 * \brief       Generic creation routines for star particles.
 * \details     Star formation rates are calculated in sfr_eEOS for the
 *              multiphase model.
 *              contains functions:
 *                void sfr_init()
 *                void sfr_create_star_particles(void)
 *                void convert_cell_into_star(int i, double birthtime)
 *                void spawn_star_from_cell(int igas, double birthtime, int
 *                  istar, MyDouble mass_of_star)
 *                void make_star(int idx, int i, double prob, MyDouble
 *                  mass_of_star, double *sum_mass_stars)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 07.06.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../gravity/forcetree.h"

#ifdef USE_SFR

static int stars_spawned;           /*!< local number of star particles spawned in the time step */
static int tot_stars_spawned;       /*!< global number of star paricles spawned in the time step */
static int stars_converted;         /*!< local number of gas cells converted into stars in the time step */
static int tot_stars_converted;     /*!< global number of gas cells converted into stars in the time step */
static int altogether_spawned;      /*!< local number of star+wind particles spawned in the time step */
static int tot_altogether_spawned;  /*!< global number of star+wind particles spawned in the time step */
static double cum_mass_stars = 0.0; /*!< cumulative mass of stars created in the time step (global value) */

static int sfr_init_called = 0;

/*! \brief Initialization routine.
 *
 *  \return void
 */
void sfr_init()
{
  if(sfr_init_called)
    return;

  sfr_init_called = 1;

  init_clouds();
}

/*! \brief This routine creates star particles according to their
 *         respective rates.
 *
 *  This function loops over all the active gas cells. If in a given cell the
 *  SFR is greater than zero, the probability of forming a star is computed
 *  and the corresponding particle is created stichastically according to the
 *  model in Springel & Hernquist (2003, MNRAS). It also saves information
 *  about the formed stellar mass and the star formation rate in the file
 *  FdSfr.
 *
 *  \return void
 */
void sfr_create_star_particles(void)
{
  TIMER_START(CPU_COOLINGSFR);

  int idx, i, bin;
  double dt, dtime;
  MyDouble mass_of_star;
  double sum_sm, total_sm, rate, sum_mass_stars, total_sum_mass_stars;
  double p = 0, pall = 0, prob, p_decide;
  double rate_in_msunperyear;
  double sfrrate, totsfrrate;

  stars_spawned = stars_converted = 0;
  sum_sm = sum_mass_stars = 0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i >= 0)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue; /* skip cells that have been swallowed or eliminated */

#ifdef SFR_KEEP_CELLS
          if(P[i].Mass < 0.3 * All.TargetGasMass)
            continue;
#endif /* #ifdef SFR_KEEP_CELLS */

          dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

          /*  the actual time-step */

          dtime = All.cf_atime * dt / All.cf_time_hubble_a;

          mass_of_star = 0;
          prob         = 0;
          p            = 0;
          pall         = 0;

          if(SphP[i].Sfr > 0)
            {
              p    = SphP[i].Sfr / ((All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR)) * dtime / P[i].Mass;
              pall = p;
              sum_sm += P[i].Mass * (1 - exp(-p));

#if defined(REFINEMENT_SPLIT_CELLS) && defined(REFINEMENT_MERGE_CELLS)

              if(P[i].Mass < 2.0 * All.TargetGasMass)
#ifdef SFR_KEEP_CELLS
                mass_of_star = 0.9 * P[i].Mass;
#else  /* #ifdef SFR_KEEP_CELLS */
                mass_of_star = P[i].Mass;
#endif /* #ifdef SFR_KEEP_CELLS */
              else
                mass_of_star = All.TargetGasMass;

#ifdef REFINEMENT_HIGH_RES_GAS
              if(SphP[i].HighResMass < HIGHRESMASSFAC * P[i].Mass)
                {
                  /* this cell does not appear to be in the high-res region.
                     If we form a star, then it is given the mass of the cell,
                     and later we give the star the SofteningType=3 particle to give it large softening */
#ifdef SFR_KEEP_CELLS
                  mass_of_star = 0.9 * P[i].Mass;
#else  /* #ifdef SFR_KEEP_CELLS */
                  mass_of_star = P[i].Mass;
#endif /* #ifdef SFR_KEEP_CELLS #else */
                }

#endif /* #ifdef REFINEMENT_HIGH_RES_GAS */

#else  /* #if defined(REFINEMENT_SPLIT_CELLS) && defined(REFINEMENT_MERGE_CELLS) */
              mass_of_star = P[i].Mass;
#endif /* #if defined(REFINEMENT_SPLIT_CELLS) && defined(REFINEMENT_MERGE_CELLS) #else */

#ifdef SFR_KEEP_CELLS
              if(P[i].Mass < 0.5 * All.TargetGasMass)
                continue; /* do not make stars from cells that should be derefined */
#endif                    /* #ifdef SFR_KEEP_CELLS */

              prob = P[i].Mass / mass_of_star * (1 - exp(-pall));
            }

          if(prob == 0)
            continue;

          if(prob < 0)
            terminate("prob < 0");

          if(prob > 1)
            {
              printf(
                  "SFR: Warning, need to make a heavier star than desired. Task=%d prob=%g P[i].Mass=%g mass_of_star=%g "
                  "mass_of_star_new=%g p=%g pall=%g\n",
                  ThisTask, prob, P[i].Mass, mass_of_star, P[i].Mass * (1 - exp(-pall)), p, pall);
              mass_of_star = P[i].Mass * (1 - exp(-pall));
              prob         = 1.0;
            }

          /* decide what process to consider (currently available: make a star or kick to wind) */
          p_decide = get_random_number();

          if(p_decide < p / pall) /* ok, it is decided to consider star formation */
            make_star(idx, i, prob, mass_of_star, &sum_mass_stars);
        }
    } /* end of main loop over active gas particles */

  int in[4], out[4], cnt = 2;
  in[0] = stars_spawned;
  in[1] = stars_converted;

  MPI_Allreduce(in, out, cnt, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  tot_stars_spawned   = out[0];
  tot_stars_converted = out[1];

  if(tot_stars_spawned > 0 || tot_stars_converted > 0)
    mpi_printf("SFR: spawned %d stars, converted %d gas particles into stars\n", tot_stars_spawned, tot_stars_converted);

  tot_altogether_spawned = tot_stars_spawned;
  altogether_spawned     = stars_spawned;

  if(tot_altogether_spawned)
    {
      /* need to assign new unique IDs to the spawned stars */

      int *list;

      if(All.MaxID == 0) /* MaxID not calculated yet */
        calculate_maxid();

      list = mymalloc("list", NTask * sizeof(int));

      MPI_Allgather(&altogether_spawned, 1, MPI_INT, list, 1, MPI_INT, MPI_COMM_WORLD);

      MyIDType newid = All.MaxID + 1;

      for(i = 0; i < ThisTask; i++)
        newid += list[i];

      myfree(list);

      for(i = 0; i < altogether_spawned; i++)
        {
          P[NumPart + i].ID = newid;

          newid++;
        }

      All.MaxID += tot_altogether_spawned;
    }

  /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
  if(tot_stars_spawned > 0 || tot_stars_converted > 0)
    {
      All.TotNumPart += tot_stars_spawned;
      All.TotNumGas -= tot_stars_converted;
      NumPart += stars_spawned;
    }

  for(bin = 0, sfrrate = 0; bin < TIMEBINS; bin++)
    if(TimeBinsHydro.TimeBinCount[bin])
      sfrrate += TimeBinSfr[bin];

  double din[3] = {sfrrate, sum_sm, sum_mass_stars}, dout[3];

  MPI_Reduce(din, dout, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      totsfrrate           = dout[0];
      total_sm             = dout[1];
      total_sum_mass_stars = dout[2];

      if(All.TimeStep > 0)
        rate = total_sm / (All.TimeStep / All.cf_time_hubble_a);
      else
        rate = 0;

      /* compute the cumulative mass of stars */
      cum_mass_stars += total_sum_mass_stars;

      /* convert to solar masses per yr */
      rate_in_msunperyear = rate * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

      fprintf(FdSfr, "%14e %14e %14e %14e %14e %14e\n", All.Time, total_sm, totsfrrate, rate_in_msunperyear, total_sum_mass_stars,
              cum_mass_stars);
      myflush(FdSfr);
    }

  TIMER_STOP(CPU_COOLINGSFR);
}

/*! \brief Convert a cell into a star.
 *
 *  This function converts an active star-forming gas cell into a star.
 *  The particle information of the gas cell is copied to the
 *  location star and the fields necessary for the creation of the star
 *  particle are initialized.
 *
 *  \param[in] i Index of the gas cell to be converted.
 *  \param[in] birthtime Time of birth (in code units) of the stellar particle.
 *
 *  \return void
 */
void convert_cell_into_star(int i, double birthtime)
{
  P[i].Type          = 4;
  P[i].SofteningType = All.SofteningTypeOfPartType[P[i].Type];

#if defined(REFINEMENT_HIGH_RES_GAS)
  if(SphP[i].HighResMass < HIGHRESMASSFAC * P[i].Mass)
    {
      /* this cell does not appear to be in the high-res region.
         We give the star the SofteningType=3 particle to give it large softening */
      P[i].SofteningType = All.SofteningTypeOfPartType[3];
    }
#endif /* #if defined(REFINEMENT_HIGH_RES_GAS) */

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  if(((1 << P[i].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
    P[i].SofteningType = get_softening_type_from_mass(P[i].Mass);
#endif /* #ifdef INDIVIDUAL_GRAVITY_SOFTENING */

  TimeBinSfr[P[i].TimeBinHydro] -= SphP[i].Sfr;

  voronoi_remove_connection(i);

  return;
}

/*! \brief Spawn a star particle from a gas cell.
 *
 *  This function spawns a star particle from an active star-forming
 *  cell. The particle information of the gas cell is copied to the
 *  location istar and the fields necessary for the creation of the star
 *  particle are initialized. The conserved variables of the gas cell
 *  are then updated according to the mass ratio between the two components
 *  to ensure conservation.
 *
 *  \param[in] igas Index of the gas cell from which the star is spawned.
 *  \param[in] birthtime Time of birth (in code units) of the stellar particle.
 *  \param[in] istar Index of the spawned stellar particle.
 *  \param[in] mass_of_star The mass of the spawned stellar particle.
 *
 *  \return void
 */
void spawn_star_from_cell(int igas, double birthtime, int istar, MyDouble mass_of_star)
{
  P[istar]               = P[igas];
  P[istar].Type          = 4;
  P[istar].SofteningType = All.SofteningTypeOfPartType[P[istar].Type];
  P[istar].Mass          = mass_of_star;

#if defined(REFINEMENT_HIGH_RES_GAS)
  if(SphP[igas].HighResMass < HIGHRESMASSFAC * P[igas].Mass)
    {
      /* this cell does not appear to be in the high-res region.
         We give the star the SofteningType=3 particle to give it large softening */
      P[istar].SofteningType = All.SofteningTypeOfPartType[3];
    }
#endif /* #if defined(REFINEMENT_HIGH_RES_GAS) */

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  if(((1 << P[istar].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
    P[istar].SofteningType = get_softening_type_from_mass(P[istar].Mass);
#endif /* #ifdef INDIVIDUAL_GRAVITY_SOFTENING */

  timebin_add_particle(&TimeBinsGravity, istar, igas, P[istar].TimeBinGrav, TimeBinSynchronized[P[istar].TimeBinGrav]);

  /* now change the conserved quantities in the cell in proportion */
  double fac = (P[igas].Mass - P[istar].Mass) / P[igas].Mass;

#ifdef MHD
  double Emag = 0.5 * (SphP[igas].B[0] * SphP[igas].B[0] + SphP[igas].B[1] * SphP[igas].B[1] + SphP[igas].B[2] * SphP[igas].B[2]) *
                SphP[igas].Volume * All.cf_atime;
  SphP[igas].Energy -= Emag;
#endif /* #ifdef MHD */

  P[igas].Mass *= fac;
  SphP[igas].Energy *= fac;
  SphP[igas].Momentum[0] *= fac;
  SphP[igas].Momentum[1] *= fac;
  SphP[igas].Momentum[2] *= fac;

#ifdef MHD
  SphP[igas].Energy += Emag;
#endif /* #ifdef MHD */

#ifdef MAXSCALARS
  for(int s = 0; s < N_Scalar; s++) /* Note, the changes in MATERIALS, HIGHRESGASMASS, etc., are treated as part of the Scalars */
    *(MyFloat *)(((char *)(&SphP[igas])) + scalar_elements[s].offset_mass) *= fac;
#endif /* #ifdef MAXSCALARS */

  return;
}

/*! \brief Make a star particle from a gas cell.
 *
 *  Given a gas cell where star formation is active and the probability
 *  of forming a star, this function selectes either to convert the gas
 *  cell into a star particle or to spawn a star depending on the
 *  target mass for the star.
 *
 *  \param[in] idx Index of the gas cell in the hydro list of active cells.
 *  \param[in] i Index of the gas cell.
 *  \param[in] prob Probability of making a star.
 *  \param[in] mass_of_star Desired mass of the star particle.
 *  \param[in, out] sum_mass_stars Holds the mass of all the stars created at the
 *             current time-step (for the local task)
 *
 *  \return void
 */
void make_star(int idx, int i, double prob, MyDouble mass_of_star, double *sum_mass_stars)
{
  if(mass_of_star > P[i].Mass)
    terminate("mass_of_star > P[i].Mass");

  if(get_random_number() < prob)
    {
      if(mass_of_star == P[i].Mass)
        {
          /* here we turn the gas particle itself into a star particle */
          Stars_converted++;
          stars_converted++;

          *sum_mass_stars += P[i].Mass;

          convert_cell_into_star(i, All.Time);
          timebin_remove_particle(&TimeBinsHydro, idx, P[i].TimeBinHydro);
        }
      else
        {
          /* in this case we spawn a new star particle, only reducing the mass in the cell by mass_of_star */
          altogether_spawned = stars_spawned;
          if(NumPart + altogether_spawned >= All.MaxPart)
            terminate("NumPart=%d spwawn %d particles no space left (All.MaxPart=%d)\n", NumPart, altogether_spawned, All.MaxPart);

          int j = NumPart + altogether_spawned; /* index of new star */

          spawn_star_from_cell(i, All.Time, j, mass_of_star);

          *sum_mass_stars += mass_of_star;
          stars_spawned++;
        }
    }
}

#endif /* #ifdef USE_SFR */
