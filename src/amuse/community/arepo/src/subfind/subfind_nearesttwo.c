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
 * \file        src/subfind/subfind_nearesttwo.c
 * \date        05/2018
 * \brief       Neighbor finding of particles in group.
 * \details     contains functions:
 *                static void particle2in(data_in * in, int i, int firstnode)
 *                static void out2particle(data_out * out, int i, int mode)
 *                static void kernel_local(void)
 *                static void kernel_imported(void)
 *                void subfind_find_nearesttwo(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 14.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef SUBFIND
#include "subfind.h"

static int subfind_nearesttwo_evaluate(int target, int mode, int threadid);

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
  MyIDType ID;
  MyFloat Hsml;
  MyFloat Density;
  MyFloat Dist[2];
  int Count;
  long long Index[2];

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
  int k;

#ifdef CELL_CENTER_GRAVITY
  if(P[i].Type == 0)
    {
      in->Pos[0] = PS[i].Center[0];
      in->Pos[1] = PS[i].Center[1];
      in->Pos[2] = PS[i].Center[2];
    }
  else
#endif /* #ifdef CELL_CENTER_GRAVITY */
    {
      in->Pos[0] = P[i].Pos[0];
      in->Pos[1] = P[i].Pos[1];
      in->Pos[2] = P[i].Pos[2];
    }

  in->Hsml    = PS[i].Hsml;
  in->ID      = P[i].ID;
  in->Density = PS[i].Density;
  in->Count   = NgbLoc[i].count;
  for(k = 0; k < NgbLoc[i].count; k++)
    {
      in->Dist[k]  = R2Loc[i].dist[k];
      in->Index[k] = NgbLoc[i].index[k];
    }
  in->Firstnode = firstnode;
}

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{
  MyFloat Dist[2];
  long long Index[2];
  int Count;
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
      int k;

      NgbLoc[i].count = out->Count;

      for(k = 0; k < out->Count; k++)
        {
          R2Loc[i].dist[k]   = out->Dist[k];
          NgbLoc[i].index[k] = out->Index[k];
        }
    }
  else /* combine */
    {
      int k, l;

      for(k = 0; k < out->Count; k++)
        {
          if(NgbLoc[i].count >= 1)
            if(NgbLoc[i].index[0] == out->Index[k])
              continue;

          if(NgbLoc[i].count == 2)
            if(NgbLoc[i].index[1] == out->Index[k])
              continue;

          if(NgbLoc[i].count < 2)
            {
              l = NgbLoc[i].count;
              NgbLoc[i].count++;
            }
          else
            {
              if(R2Loc[i].dist[0] > R2Loc[i].dist[1])
                l = 0;
              else
                l = 1;

              if(out->Dist[k] >= R2Loc[i].dist[l])
                continue;
            }

          R2Loc[i].dist[l]   = out->Dist[k];
          NgbLoc[i].index[l] = out->Index[k];

          if(NgbLoc[i].count == 2)
            if(NgbLoc[i].index[0] == NgbLoc[i].index[1])
              terminate("this is not supposed to happen");
        }
    }
}

#define USE_SUBCOMM_COMMUNICATOR
#include "../utils/generic_comm_helpers2.h"

static double *Dist2list;
static int *Ngblist;

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int i;
  {
    int j, threadid = get_thread_num();

    for(j = 0; j < SubNTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

        i = NextParticle++;

        if(i >= NumPartGroup)
          break;

        subfind_nearesttwo_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        subfind_nearesttwo_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/*! \brief Neighbour finding for each particle in group.
 *
 *  \return void
 */
void subfind_find_nearesttwo(void)
{
  if(SubThisTask == 0)
    printf("SUBFIND-COLLECTIVE, root-task=%d: Start finding nearest two.\n", ThisTask);

  /* allocate buffers to arrange communication */

  Ngblist   = (int *)mymalloc("Ngblist", NumPartGroup * sizeof(int));
  Dist2list = (double *)mymalloc("Dist2list", NumPartGroup * sizeof(double));

  generic_set_MaxNexport();

  for(int i = 0; i < NumPartGroup; i++)
    NgbLoc[i].count = 0;

  generic_comm_pattern(NumPartGroup, kernel_local, kernel_imported);

  myfree(Dist2list);
  myfree(Ngblist);

  if(SubThisTask == 0)
    printf("SUBFIND-COLLECTIVE, root-task=%d: Done with nearest two.\n", ThisTask);
}

/*! \brief Neighbor finding routine on local particles.
 *
 *  \param[in] target Index of particle/cell.
 *  \param[in] mode Flag if it operates on local or imported data.
 *  \param[in] threadid ID of thread.
 *
 * \return 0
 */
static int subfind_nearesttwo_evaluate(int target, int mode, int threadid)
{
  int j, k, n, no, count;
  MyIDType ID;
  long long index[2];
  double dist[2];
  int numngb, numnodes, *firstnode;
  double hsml;
  double density;
  MyDouble *pos;
  struct NODE *current;
  double dx, dy, dz, disthsml, r2;
  MyDouble xtmp, ytmp, ztmp;

  data_in local, *in;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  ID      = in->ID;
  density = in->Density;
  pos     = in->Pos;
  hsml    = in->Hsml;
  count   = in->Count;
  for(k = 0; k < count; k++)
    {
      dist[k]  = in->Dist[k];
      index[k] = in->Index[k];
    }

  if(count == 2)
    if(index[0] == index[1])
      {
        terminate("task=%d/%d target=%d mode=%d  index_0=%lld  index_1=%lld\n", SubThisTask, ThisTask, target, mode, index[0],
                  index[1]);
      }

  numngb = 0;
  count  = 0;

  hsml *= 1.00001; /* prevents that the most distant neighbour on the edge of the search region may not be found.
                    * (needed for consistency with serial algorithm)
                    */

  for(k = 0; k < numnodes; k++)
    {
      if(mode == MODE_LOCAL_PARTICLES)
        {
          no = SubTree_MaxPart; /* root node */
        }
      else
        {
          no = firstnode[k];
          no = SubNodes[no].u.d.nextnode; /* open it */
        }
      while(no >= 0)
        {
          if(no < SubTree_MaxPart) /* single particle */
            {
              int p = no;
              no    = SubNextnode[no];

              disthsml = hsml;
              dx       = FOF_NEAREST_LONG_X(SubTree_Pos_list[3 * p + 0] - pos[0]);
              if(dx > disthsml)
                continue;
              dy = FOF_NEAREST_LONG_Y(SubTree_Pos_list[3 * p + 1] - pos[1]);
              if(dy > disthsml)
                continue;
              dz = FOF_NEAREST_LONG_Z(SubTree_Pos_list[3 * p + 2] - pos[2]);
              if(dz > disthsml)
                continue;
              if((r2 = (dx * dx + dy * dy + dz * dz)) > disthsml * disthsml)
                continue;

              Dist2list[numngb] = r2;
              Ngblist[numngb++] = p;
            }
          else if(no < SubTree_MaxPart + SubTree_MaxNodes) /* internal node */
            {
              if(mode == 1)
                {
                  if(no < SubTree_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the
                                                           branch */
                    {
                      break;
                    }
                }

              current = &SubNodes[no];

              no = current->u.d.sibling; /* in case the node can be discarded */

              disthsml = hsml + 0.5 * current->len;

              dx = FOF_NEAREST_LONG_X(current->center[0] - pos[0]);
              if(dx > disthsml)
                continue;
              dy = FOF_NEAREST_LONG_Y(current->center[1] - pos[1]);
              if(dy > disthsml)
                continue;
              dz = FOF_NEAREST_LONG_Z(current->center[2] - pos[2]);
              if(dz > disthsml)
                continue;
              /* now test against the minimal sphere enclosing everything */
              disthsml += FACT1 * current->len;
              if(dx * dx + dy * dy + dz * dz > disthsml * disthsml)
                continue;

              no = current->u.d.nextnode; /* ok, we need to open the node */
            }
          else if(no >= SubTree_ImportedNodeOffset) /* point from imported nodelist */
            {
              terminate("do not expect imported points here");
            }
          else /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES");

              if(target >= 0) /* note: if no target is given, export will not occur */
                subfind_treefind_collective_export_node_threads(no, target, threadid);

              no = SubNextnode[no - SubTree_MaxNodes];
            }
        }
    }

  for(n = 0; n < numngb; n++)
    {
      j  = Ngblist[n];
      r2 = Dist2list[n];

      if(P[j].ID != ID) /* exclude the self-particle */
        {
          if(PS[j].Density > density) /* we only look at neighbours that are denser */
            {
              if(count < 2)
                {
                  dist[count]  = r2;
                  index[count] = (((long long)SubThisTask) << 32) + j;
                  count++;
                }
              else
                {
                  if(dist[0] > dist[1])
                    k = 0;
                  else
                    k = 1;

                  if(r2 < dist[k])
                    {
                      dist[k]  = r2;
                      index[k] = (((long long)SubThisTask) << 32) + j;
                    }
                }
            }
        }
    }

  out.Count = count;
  for(k = 0; k < count; k++)
    {
      out.Dist[k]  = dist[k];
      out.Index[k] = index[k];
    }

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif /* #ifdef SUBFIND */
