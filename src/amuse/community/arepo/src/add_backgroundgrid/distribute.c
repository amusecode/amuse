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
 * \file        src/add_backgroundgrid/distribute.c
 * \date        05/2018
 * \brief       Distributes the cell properties in an SPH kernel weighted
 *              fashion to neighboring cells.
 * \details     contains functions:
 *                static void particle2in(data_in * in, int i, int firstnode)
 *                static void out2particle(data_out * out, int i, int mode)
 *                static void kernel_local(void)
 *                static void kernel_imported(void)
 *                void distribute_particles(void)
 *                int find_cells_evaluate(int target, int mode, int thread_id)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 11.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <mpi.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "add_bggrid.h"

#ifdef ADDBACKGROUNDGRID

static int find_cells_evaluate(int target, int mode, int thread_id);

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;
  MyFloat Weight;
  MyFloat Mass;
  MyFloat InternalEnergy;
  MyFloat Momentum[3];
#ifdef MHD
  MyFloat B[3];
#endif /* #ifdef MHD */
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/*! \brief Routine that fills the relevant particle/cell data into the input
 *         structure defined above. Needed by generic_comm_helpers2.
 *
 *  \param[out] in Data structure to fill.
 *  \param[in] i Index of particle in P and SphP arrays.
 *  \param[in] firstnode First note of communication.
 *
 *  \return void
 */
static void particle2in(data_in *in, int i, int firstnode)
{
  in->Pos[0] = P[i].Pos[0];
  in->Pos[1] = P[i].Pos[1];
  in->Pos[2] = P[i].Pos[2];

  in->Hsml = SphP[i].Hsml;

  in->Weight         = SphP[i].Weight;
  in->Mass           = P[i].Mass;
  in->InternalEnergy = SphP[i].Utherm * P[i].Mass;

  int k;
  for(k = 0; k < 3; k++)
    in->Momentum[k] = P[i].Vel[k] * P[i].Mass;

#ifdef MHD
  for(k = 0; k < 3; k++)
    in->B[k] = SphP[i].B[k];
#endif /* #ifdef MHD */

  in->Firstnode = firstnode;
}

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{
  char nothing;
} data_out;

static data_out *DataResult, *DataOut;

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *  particle and cell data arrays (P, SphP,...)
 *  \param[in] i Index of particle in P and SphP arrays
 *  \param[in] mode Mode of function: local particles or information that was
 *  communicated from other tasks and has to be added locally?
 *
 *  \return void
 */
static void out2particle(data_out *out, int i, int mode) { return; }

#include "../utils/generic_comm_helpers2.h"

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int idx;
  {
    int j, threadid = get_thread_num();
    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

        idx = NextParticle++;

        if(idx >= TimeBinsGravity.NActiveParticles)
          break;

        int i = TimeBinsGravity.ActiveParticleList[idx];
        if(i < 0)
          continue;

        find_cells_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

/*! \brief Routine that defines what to do with imported particles.
 *
 *  Calls the *_evaluate function in MODE_IMPORTED_PARTICLES.
 *
 *  \return void
 */
static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
  {
    int threadid = get_thread_num();

    while(1)
      {
        i = cnt++;

        if(i >= Nimport)
          break;

        find_cells_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/*! \brief Main function to distribute hydro quantities over a kernel average.
 *
 *  \return void
 */
void distribute_particles(void)
{
  mpi_printf("ADD BACKGROUND GRID: distributing the fluid quantities\n");

  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsGravity.NActiveParticles, kernel_local, kernel_imported);

#ifdef MHD
  /* now divide the B field in each cell by the weight (sum of the wk's,
     which we stored in SphP.divB */
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].ID >= IDNew)
        {
          int j;
          if(SphP[i].DivB > 0)
            for(j = 0; j < 3; j++)
              SphP[i].B[j] /= SphP[i].DivB;
        }
    }
#endif /* #ifdef MHD */

  mpi_printf("ADD BACKGROUND GRID: done\n");
}

/*! \brief Distributes imported properties on neighbouring cells.
 *
 *  \param[in] target Index of particle/cell.
 *  \param[in] mode Flag if it operates on local or imported data.
 *  \param[in] threadid ID of thread.
 *
 *  \return 0
 */
int find_cells_evaluate(int target, int mode, int thread_id)
{
  int j, n, numnodes, *firstnode;
  double h, h2, hinv, hinv3;
  MyDouble dx, dy, dz, r;
  MyDouble *pos;
  double xtmp, ytmp, ztmp;

  data_in local, *target_data;
  data_out out;
  out.nothing = 0;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos  = target_data->Pos;
  h    = target_data->Hsml;
  h2   = h * h;
  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else  /* #ifndef  TWODIMS */
  hinv3 = hinv * hinv / boxSize_Z;
#endif /* #ifndef  TWODIMS #else */

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  double wsum = 0;

  for(n = 0; n < nfound; n++)
    {
      j = Thread[thread_id].Ngblist[n];

      if(P[j].ID < IDNew)
        continue;

      dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
      dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
      dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

      double r2 = dx * dx + dy * dy + dz * dz;

      if(r2 < h2)
        {
          r = sqrt(r2);

          double u = r * hinv;
          double wk;
          if(u < 0.5)
            wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
          else
            wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

          double weight = SphP[j].Volume * wk / target_data->Weight;

          wsum += weight;

          P[j].Mass += target_data->Mass * weight;
          SphP[j].Energy += target_data->InternalEnergy * weight;

          int k;
          for(k = 0; k < 3; k++)
            SphP[j].Momentum[k] += target_data->Momentum[k] * weight;

#ifdef MHD
          for(k = 0; k < 3; k++)
            SphP[j].B[k] += target_data->B[k] * weight;
          SphP[j].DivB += wk;
#endif /* #ifdef MHD */
        }
    }

  if(wsum > 1.01)
    {
      printf("wsum=%g, Weight=%g, target=%d\n", wsum, target_data->Weight, target);
      terminate("bla");
    }

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif /* #ifdef ADDBACKGROUNDGRID */
