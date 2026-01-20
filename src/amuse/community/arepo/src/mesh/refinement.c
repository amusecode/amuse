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
 * \file        src/mesh/refinement.c
 * \date        05/2018
 * \brief       Driver routines that handle refinement and de-refinement.
 * \details     contains functions:
 *                void do_derefinements_and_refinements()
 *                void refinement_prepare()
 *                void refinement_cleanup()
 *                void move_collisionless_particle(int new_i, int old_i)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 06.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include "../main/allvars.h"

#ifdef REFINEMENT
#include "../main/proto.h"

#if defined(REFINEMENT_MERGE_CELLS) && defined(REFINEMENT_SPLIT_CELLS)
char *FlagDoNotRefine;
#endif /* #if defined (REFINEMENT_MERGE_CELLS) && defined (REFINEMENT_SPLIT_CELLS) */

static void refinement_prepare();
static void refinement_cleanup();

/*! \brief Main routine to trigger refinement and de-refinements.
 *
 *  Called in main run loop (run.c).
 *
 *  \return void
 */
void do_derefinements_and_refinements()
{
  refinement_prepare();

#ifdef REFINEMENT_MERGE_CELLS
  do_derefinements();
#endif /* #ifdef REFINEMENT_MERGE_CELLS */

#ifdef REFINEMENT_SPLIT_CELLS
  do_refinements();
#endif /* #ifdef REFINEMENT_SPLIT_CELLS */

  refinement_cleanup();
}

/*! \brief Prepares for refinement.
 *
 *  Determines quantities needed by refinement routine;
 *  Allocates additional arrays.
 *
 *  \return void
 */
void refinement_prepare()
{
  TIMER_START(CPU_REFINE);

#ifdef REFINEMENT_VOLUME_LIMIT
  int idx, i;
#endif /* #ifdef REFINEMENT_VOLUME_LIMIT */

#if defined(REFINEMENT_MERGE_CELLS) && defined(REFINEMENT_SPLIT_CELLS)
  FlagDoNotRefine = mymalloc_movable(&FlagDoNotRefine, "FlagDoNotRefine", NumGas * sizeof(char));
#endif /* #if defined (REFINEMENT_MERGE_CELLS) && defined (REFINEMENT_SPLIT_CELLS) */

#ifdef REFINEMENT_VOLUME_LIMIT
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      SphP[i].MinNgbVolume = MAX_REAL_NUMBER;

      int q = SphP[i].first_connection;
      while(q >= 0)
        {
          int dp       = DC[q].dp_index;
          int particle = Mesh.DP[dp].index;

          if(particle < 0)
            {
              if(q == SphP[i].last_connection)
                break;

              q = DC[q].next;
              continue;
            }

          if(particle >= NumGas && Mesh.DP[dp].task == ThisTask)
            particle -= NumGas;

          double Volume;
          if(DC[q].task == ThisTask)
            Volume = SphP[particle].Volume;
          else
            {
#ifndef OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT
              Volume = PrimExch[particle].Volume;
#else  /* #ifndef OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT */
              Volume = RefExch[particle].Volume;
#endif /* #ifndef OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT #else */
            }

          if(Volume < SphP[i].MinNgbVolume)
            SphP[i].MinNgbVolume = Volume;

          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }
    }
#endif /* #ifdef REFINEMENT_VOLUME_LIMIT */

  TIMER_STOP(CPU_REFINE);
}

/*! \brief Cleans up after refinement.
 *
 *  Frees memory allocated by refinement_prepare().
 *
 *  \return void
 */
void refinement_cleanup()
{
#if defined(REFINEMENT_MERGE_CELLS) && defined(REFINEMENT_SPLIT_CELLS)
  myfree(FlagDoNotRefine);
#endif /* #if defined (REFINEMENT_MERGE_CELLS) && defined (REFINEMENT_SPLIT_CELLS) */
}

/*! \brief Moves collisionless particle from index old_i to new_i.
 *
 *  Needed if new cell is introduced, as cells have to be at the beginning of
 *  the P array and all other particles have to be located after the last
 *  gas cell. This routine moves not only data in P and SphP, but also updates
 *  the time-bin data consistently.
 *
 *  \param[in] new_i New index of particle in P.
 *  \param[in] old_i Previous index of particle in P.
 *
 *  \return void
 */
void move_collisionless_particle(int new_i, int old_i)
{
  int prev, next, bin;
  struct TimeBinData *tbData;

  P[new_i] = P[old_i];

  if(P[old_i].Mass == 0 && P[old_i].ID == 0)
    return;

  if(P[old_i].Mass == 0 && P[old_i].Type == 4)
    return;

  tbData = &TimeBinsGravity;
  bin    = P[old_i].TimeBinGrav;

  if(TimeBinSynchronized[bin])
    {
      /* particle is active, need to add it to the list of active particles again
         we assume here, that the new particle at the old index in this list is also active! */
      tbData->ActiveParticleList[tbData->NActiveParticles] = new_i;
      tbData->NActiveParticles++;
    }

  /* now move it in the link list of its timebin
     we only need to change the gravity timebin here */

  tbData->NextInTimeBin[new_i] = tbData->NextInTimeBin[old_i];
  tbData->PrevInTimeBin[new_i] = tbData->PrevInTimeBin[old_i];

  prev = tbData->PrevInTimeBin[old_i];
  next = tbData->NextInTimeBin[old_i];

  if(prev >= 0)
    tbData->NextInTimeBin[prev] = new_i;
  else
    {
      if(tbData->FirstInTimeBin[bin] != old_i)
        terminate("strange");
      tbData->FirstInTimeBin[bin] = new_i;
    }

  if(next >= 0)
    tbData->PrevInTimeBin[next] = new_i;
  else
    {
      if(tbData->LastInTimeBin[bin] != old_i)
        terminate("strange");
      tbData->LastInTimeBin[bin] = new_i;
    }
}

#endif /* REFINEMENT */
