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
 * \file        src/ngbtree/ngbtree_search.c
 * \date        05/2018
 * \brief       This file contains a search routine on the neighbor tree.
 * \details     contains functions:
 *                static void particle2in(data_in * in, int i, int firstnode)
 *                static void out2particle(data_out * out, int i, int mode)
 *                static void kernel_local(void)
 *                static void kernel_imported(void)
 *                void find_nearest_meshpoint_global(mesh_search_data *
 *                  searchdata_input, int nn, int hsmlguess, int verbose)
 *                int ngbsearch_primary_cell_evaluate(int target, int mode,
 *                  int threadid)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 21.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/* temporary particle arrays */
static MyDouble *ngbsearch_nearest_dist;
static MyDouble *ngbsearch_hsml;
static mesh_search_data *searchdata;

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble pos[3];   /* tracer particle position */
  MyDouble hsml;     /* current search radius */
  MyDouble distance; /* nearest neighbor distance */

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
  in->pos[0] = searchdata[i].Pos[0];
  in->pos[1] = searchdata[i].Pos[1];
  in->pos[2] = searchdata[i].Pos[2];

  in->hsml     = ngbsearch_hsml[i];
  in->distance = ngbsearch_nearest_dist[i];

  in->Firstnode = firstnode;
}

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Distance; /* distance to closest cell on task */
  int Task;
  int Index;
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
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      if(out->Index >= 0)
        {
          ngbsearch_nearest_dist[i] = out->Distance;
          searchdata[i].Task        = out->Task;
          searchdata[i].u.Index     = out->Index;
        }
    }
  else /* combine */
    {
      /* closer cell on other task? */
      if(out->Distance < ngbsearch_nearest_dist[i])
        {
          ngbsearch_nearest_dist[i] = out->Distance;
          searchdata[i].Task        = out->Task;
          searchdata[i].u.Index     = out->Index;
        }
    }
}

#include "../utils/generic_comm_helpers2.h"

static int ngbsearch_primary_cell_evaluate(int target, int mode, int threadid);
static int n;

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int i;

  /* do local particles */
  {
    int j, threadid = get_thread_num();

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

        i = NextParticle++;

        if(i >= n)
          break;

        if(searchdata[i].Task == -1)
          ngbsearch_primary_cell_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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
  int i, cnt = 0;
  {
    int threadid = get_thread_num();

    while(1)
      {
        i = cnt++;

        if(i >= Nimport)
          break;

        ngbsearch_primary_cell_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/*! \brief Searches the cells at the positions in searchdata.
 *
 *  This function searches the cells which are at the positions specified in
 *  searchdata. The Pos field must be set. After the search is performed the
 *  Task and Index field contain the task/index of the cell at position Pos.
 *  If hsmlguess=1 initial search radius is read from Index/Hsml union in
 *  searchdata.
 *
 *  \param[in] searchdata_input Contains the search positions, after function
 *             call the fields Task and Index are set.
 *  \param[in] nn Number of items in searchdata.
 *  \param[in] hsmlguess Guess for initial search radius;
 *             1: from searchdata; else from MeanVolume of cells.
 *  \param[in] verbose More output.
 *
 *  \return void
 */
void find_nearest_meshpoint_global(mesh_search_data *searchdata_input, int nn, int hsmlguess, int verbose)
{
  int i;
  n                      = nn;
  ngbsearch_nearest_dist = mymalloc("ngbsearch_nearest_dist", n * sizeof(MyDouble));
  ngbsearch_hsml         = mymalloc("ngbsearch_hsml", n * sizeof(MyDouble));
  searchdata             = searchdata_input;

  for(i = 0; i < n; i++)
    {
      ngbsearch_nearest_dist[i] = MAX_REAL_NUMBER;

      if(hsmlguess)
        ngbsearch_hsml[i] = searchdata[i].u.hsmlguess;
      else
        ngbsearch_hsml[i] = 1e-6 * pow(All.MeanVolume, 1.0 / 3);

      searchdata[i].Task = -1;  // None found yet
    }

  generic_set_MaxNexport();

  int ntot, iter = 0;

  /* we will repeat the whole thing for those points where we did not find a nearest neighbor */
  do
    {
      generic_comm_pattern(n, kernel_local, kernel_imported);

      int npleft = 0;

      /* do final operations on results */
      for(i = 0; i < n; i++)
        {
          if(searchdata[i].Task == -1)
            {
              npleft++;
              ngbsearch_hsml[i] *= 2.0;

              if(iter >= MAXITER - 10)
                {
                  printf("i=%d task=%d hsml=%g nearest dist=%g pos=(%g|%g|%g)\n", i, ThisTask, ngbsearch_hsml[i],
                         ngbsearch_nearest_dist[i], searchdata[i].Pos[0], searchdata[i].Pos[1], searchdata[i].Pos[2]);
                  myflush(stdout);
                }
              if(iter > MAXITER)
                terminate("NGBSEARCH: iter > MAXITER");
            }
        }

      /* sum up the left overs */
      MPI_Allreduce(&npleft, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if(ntot > 0) /* ok, we need to repeat for a few particles */
        {
          iter++;
          if(iter > 0 && ThisTask == 0 && verbose)
            {
              printf("NGBSEARCH: iteration %d: need to repeat for %d points.\n", iter, ntot);
              myflush(stdout);
            }

          if(iter > MAXITER)
            terminate("NGBSEARCH: failed to converge in tracer particles\n");
        }
    }
  while(ntot > 0);

  myfree(ngbsearch_hsml);
  myfree(ngbsearch_nearest_dist);
}

/*! \brief Performs the neighbor search.
 *
 *  \param[in] target the index of the particle to process(mode 0: in
 *             searchdata, mode 1: in NgbSearchDataGet/Result).
 *  \param[in] mode either 0 (handle local particles) or 1 (handle particles
 *             sent to us).
 *  \param[in] treadid Id of thread.
 *
 *  \return 0
 */
int ngbsearch_primary_cell_evaluate(int target, int mode, int threadid)
{
  int j, n;
  int numnodes, *firstnode;
  MyDouble h, distmax;
  MyDouble dx, dy, dz, r;
  MyDouble *pos;
  data_in local, *target_data;
  data_out out;

  int index = -1;

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

  pos     = target_data->pos;
  h       = target_data->hsml;
  distmax = target_data->distance;

  int numngb = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < numngb; n++)
    {
      j = Thread[threadid].Ngblist[n];

      dx = pos[0] - P[j].Pos[0];
      dy = pos[1] - P[j].Pos[1];
      dz = pos[2] - P[j].Pos[2];

      if(dx > boxHalf_X)
        dx -= boxSize_X;
      if(dx < -boxHalf_X)
        dx += boxSize_X;
      if(dy > boxHalf_Y)
        dy -= boxSize_Y;
      if(dy < -boxHalf_Y)
        dy += boxSize_Y;
      if(dz > boxHalf_Z)
        dz -= boxSize_Z;
      if(dz < -boxHalf_Z)
        dz += boxSize_Z;

      r = sqrt(dx * dx + dy * dy + dz * dz);
      if(r < distmax && r < h && P[j].ID != 0 && P[j].Mass > 0)
        {
          distmax = r;
          index   = j;
        }
    }

  out.Distance = distmax;
  out.Task     = ThisTask;
  out.Index    = index;

  if(index < 0)
    {
      out.Distance = MAX_REAL_NUMBER;
      out.Task     = -1;
      out.Index    = -1;
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}
