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
 * \file        src/time_integration/timestep.h
 * \date        05/2018
 * \brief       Header for timestep criteria.
 * \details
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 29.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifndef TIMESTEP_H
#define TIMESTEP_H

#include "../main/allvars.h"

#ifdef ENLARGE_DYNAMIC_RANGE_IN_TIME
typedef long long integertime;
#define TIMEBINS 60
#define TIMEBASE                                                                                           \
  (((long long)1) << TIMEBINS) /* The simulated timespan is mapped onto the integer interval [0,TIMESPAN], \
                                *  where TIMESPAN needs to be a power of 2. */
#else                          /* #ifdef   ENLARGE_DYNAMIC_RANGE_IN_TIME */
typedef int integertime;
#define TIMEBINS 29
#define TIMEBASE (1 << TIMEBINS)
#endif /* #ifdef   ENLARGE_DYNAMIC_RANGE_IN_TIME #else */

/*! \brief Linked list for particles in specific timebin.
 */
struct TimeBinData
{
  int NActiveParticles;
  long long GlobalNActiveParticles;
  int *ActiveParticleList;
  int TimeBinCount[TIMEBINS];

  int FirstInTimeBin[TIMEBINS];
  int LastInTimeBin[TIMEBINS];
  int *NextInTimeBin;
  int *PrevInTimeBin;
  char Name[100];
  int *MaxPart;
};

void find_timesteps_without_gravity(void);
void update_timesteps_from_gravity(void);
integertime get_timestep_gravity(int p);
integertime get_timestep_hydro(int p);
integertime get_timestep_pm(void);
int test_if_grav_timestep_is_too_large(int p, int bin);
void validate_timestep(double dt, integertime ti_step, int p);
int get_timestep_bin(integertime ti_step);
double get_time_difference_in_Gyr(double a0, double a1);

/* TimeBinData stuff */
void timebins_init(struct TimeBinData *tbData, const char *name, int *MaxPart);
void timebins_allocate(struct TimeBinData *tbData);
void timebins_reallocate(struct TimeBinData *tbData);
void timebins_get_bin_and_do_validity_checks(integertime ti_step, int *bin_new, int bin_old);
void timebin_move_particle(struct TimeBinData *tbData, int p, int timeBin_old, int timeBin_new);
void timebin_add_particle(struct TimeBinData *tbData, int i_new, int i_old, int timeBin, int addToListOfActiveParticles);
void timebin_remove_particle(struct TimeBinData *tbData, int idx, int bin);
void timebin_cleanup_list_of_active_particles(struct TimeBinData *tbData);
void timebin_move_sfr(int p, int timeBin_old, int timeBin_new);
void timebin_make_list_of_active_particles_up_to_timebin(struct TimeBinData *tbData, int timebin);
void timebin_add_particles_of_timebin_to_list_of_active_particles(struct TimeBinData *tbData, int timebin);

#endif /* TIMESTEP */
