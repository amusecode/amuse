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
 * \file        src/fof/fof_nearest.c
 * \date        05/2018
 * \brief       Routine to find nearest primary link type particle to link
 *              secondary link type to FoF groups.
 * \details     contains functions:
 *                static void particle2in(data_in * in, int i, int firstnode)
 *                static void out2particle(data_out * out, int i, int mode)
 *                static void kernel_local(void)
 *                static void kernel_imported(void)
 *                double fof_find_nearest_dmparticle(MyIDType * vMinID, int
 *                  *vHead, int *vLen, int *vNext, int *vTail, int *vMinIDTask)
 *                static int fof_find_nearest_dmparticle_evaluate(int target,
 *                  int mode, int threadid)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 24.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_math.h>
#include <inttypes.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"
#include "../subfind/subfind.h"
#include "fof.h"

#ifdef FOF

static MyFloat *fof_nearest_distance;
static MyFloat *fof_nearest_hsml;

static MyIDType *MinID;
static int *Head, *Len, *Next, *Tail, *MinIDTask;

static int fof_find_nearest_dmparticle_evaluate(int target, int mode, int threadid);

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;

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
  in->Hsml   = fof_nearest_hsml[i];

  in->Firstnode = firstnode;
}

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{
  MyFloat Distance;
  MyIDType MinID;
  int MinIDTask;
#if defined(SUBFIND)
  MyFloat DM_Hsml;
#endif /* #if defined(SUBFIND) */
} data_out;

static data_out *DataResult, *DataOut;

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *             particle and cell data arrays (PS)
 *  \param[in] i Index of particle in P and SphP arrays
 *  \param[in] mode Mode of function: local particles or information that was
 *             communicated from other tasks and has to be added locally?
 *
 *  \return void
 */
static void out2particle(data_out *out, int i, int mode)
{
  if(out->Distance < fof_nearest_distance[i])
    {
      fof_nearest_distance[i] = out->Distance;
      MinID[i]                = out->MinID;
      MinIDTask[i]            = out->MinIDTask;
#if defined(SUBFIND)
      PS[i].Hsml = out->DM_Hsml;
#endif /* #if defined(SUBFIND) */
    }
}

#include "../utils/generic_comm_helpers2.h"

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

        if(i >= NumPart)
          break;

        if((1 << P[i].Type) & (FOF_SECONDARY_LINK_TYPES))
          {
            if(fof_nearest_distance[i] > 1.0e29) /* we haven't found any neighbor yet */
              {
                fof_find_nearest_dmparticle_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
              }
          }
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

        fof_find_nearest_dmparticle_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/*! \brief Finds nearest dark matter particle for secondary link types
 *
 *  \param[out] vMinID Pointer to MinID array.
 *  \param[in] vHead Pointer to Head array.
 *  \param[in] vLen Pointer to Len array.
 *  \param[in] vNext Pointer to Next array.
 *  \param[in] vTail Pointer to Tail array.
 *  \param[out] vMinIDTask Pointer to MinIDTask array.
 *
 *  \return Time spent in this function.
 */
double fof_find_nearest_dmparticle(MyIDType *vMinID, int *vHead, int *vLen, int *vNext, int *vTail, int *vMinIDTask)
{
  MinID     = vMinID;
  Head      = vHead;
  Len       = vLen;
  Next      = vNext;
  Tail      = vTail;
  MinIDTask = vMinIDTask;

  int i, n, npleft, iter;
  long long ntot;
  double tstart = second();

  mpi_printf("FOF: Start finding nearest dm-particle (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  fof_nearest_distance = (MyFloat *)mymalloc("fof_nearest_distance", sizeof(MyFloat) * NumPart);
  fof_nearest_hsml     = (MyFloat *)mymalloc("fof_nearest_hsml", sizeof(MyFloat) * NumPart);

  for(n = 0; n < NumPart; n++)
    {
      if((1 << P[n].Type) & (FOF_SECONDARY_LINK_TYPES))
        {
          fof_nearest_distance[n] = 1.0e30;
          if(P[n].Type == 0)
#ifdef USE_AREPO_FOF_WITH_GADGET_FIX
            fof_nearest_hsml[n] = SphP[n].Hsml;
#else  /* #ifdef USE_AREPO_FOF_WITH_GADGET_FIX */
            fof_nearest_hsml[n] = get_cell_radius(n);
#endif /* #ifdef USE_AREPO_FOF_WITH_GADGET_FIX #else */
          else
            fof_nearest_hsml[n] = 0.1 * LinkL;
        }
    }

  generic_set_MaxNexport();

  iter = 0;
  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      double t0 = second();

      generic_comm_pattern(NumPart, kernel_local, kernel_imported);

      /* do final operations on results */
      for(i = 0, npleft = 0; i < NumPart; i++)
        {
          if((1 << P[i].Type) & (FOF_SECONDARY_LINK_TYPES))
            {
              if(fof_nearest_distance[i] > 1.0e29)
                {
                  if(fof_nearest_hsml[i] < 4 * LinkL) /* we only search out to a maximum distance */
                    {
                      /* need to redo this particle */
                      npleft++;
                      fof_nearest_hsml[i] *= 2.0;
                      if(iter >= MAXITER - 10)
                        {
                          printf("FOF: i=%d task=%d ID=%d P[i].Type=%d Hsml=%g LinkL=%g nearest=%g pos=(%g|%g|%g)\n", i, ThisTask,
                                 (int)P[i].ID, P[i].Type, fof_nearest_hsml[i], LinkL, fof_nearest_distance[i], P[i].Pos[0],
                                 P[i].Pos[1], P[i].Pos[2]);
                          myflush(stdout);
                        }
                    }
                  else
                    {
                      fof_nearest_distance[i] = 0; /* we do not continue to search for this particle */
                    }
                }
            }
        }

      sumup_large_ints(1, &npleft, &ntot);

      double t1 = second();
      if(ntot > 0)
        {
          iter++;
          if(iter > 0)
            mpi_printf("FOF: fof-nearest iteration %d: need to repeat for %lld particles. (took = %g sec)\n", iter, ntot,
                       timediff(t0, t1));

          if(iter > MAXITER)
            terminate("FOF: failed to converge in fof-nearest\n");
        }
    }
  while(ntot > 0);

  myfree(fof_nearest_hsml);
  myfree(fof_nearest_distance);

  mpi_printf("FOF: done finding nearest dm-particle\n");

  double tend = second();
  return timediff(tstart, tend);
}

/*! \brief Evaluate function to finding nearest dark matter particle for
 *         secondary link types.
 *
 *  \param[in] target Index of particle/cell.
 *  \param[in] mode Flag if it operates on local or imported data.
 *  \param[in] threadid ID of thread.
 *
 *  \return 0
 */
static int fof_find_nearest_dmparticle_evaluate(int target, int mode, int threadid)
{
  int k, no, index, numnodes, *firstnode;
  double h, r2max, dist;
  double dx, dy, dz, r2;
  MyDouble *pos;
  data_in local, *target_data;
  data_out out;

  double xtmp, ytmp, ztmp;

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

  pos = target_data->Pos;
  h   = target_data->Hsml;

  index = -1;
  r2max = 1.0e30;

  /* Now start the actual tree-walk computation for this particle */

  for(k = 0; k < numnodes; k++)
    {
      if(mode == MODE_LOCAL_PARTICLES)
        {
          no = Tree_MaxPart; /* root node */
        }
      else
        {
          no = firstnode[k];
          no = Nodes[no].u.d.nextnode; /* open it */
        }

      while(no >= 0)
        {
          if(no < Tree_MaxPart) /* single particle */
            {
              int p = no;
              no    = Nextnode[no];

              if(!((1 << P[p].Type) & (FOF_SECONDARY_LINK_TARGET_TYPES)))
                continue;

              dist = h;
              dx   = FOF_NEAREST_LONG_X(Tree_Pos_list[3 * p + 0] - pos[0]);
              if(dx > dist)
                continue;
              dy = FOF_NEAREST_LONG_Y(Tree_Pos_list[3 * p + 1] - pos[1]);
              if(dy > dist)
                continue;
              dz = FOF_NEAREST_LONG_Z(Tree_Pos_list[3 * p + 2] - pos[2]);
              if(dz > dist)
                continue;

              r2 = dx * dx + dy * dy + dz * dz;
              if(r2 < r2max && r2 < h * h)
                {
                  index = p;
                  r2max = r2;
                }
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes) /* internal node */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no <
                     Tree_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the branch */
                    break;
                }

              struct NODE *current = &Nodes[no];

              no = current->u.d.sibling; /* in case the node can be discarded */

              dist = h + 0.5 * current->len;
              dx   = FOF_NEAREST_LONG_X(current->center[0] - pos[0]);
              if(dx > dist)
                continue;
              dy = FOF_NEAREST_LONG_Y(current->center[1] - pos[1]);
              if(dy > dist)
                continue;
              dz = FOF_NEAREST_LONG_Z(current->center[2] - pos[2]);
              if(dz > dist)
                continue;

              /* now test against the minimal sphere enclosing everything */
              dist += FACT1 * current->len;
              if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;

              no = current->u.d.nextnode; /* ok, we need to open the node */
            }
          else if(no >= Tree_ImportedNodeOffset) /* point from imported nodelist */
            {
              terminate("do not expect imported points here");
            }
          else /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES");

              if(target >= 0)
                tree_treefind_export_node_threads(no, target, threadid);

              no = Nextnode[no - Tree_MaxNodes];
            }
        }
    }

  if(index >= 0)
    {
      out.Distance  = sqrt(r2max);
      out.MinID     = MinID[Head[index]];
      out.MinIDTask = MinIDTask[Head[index]];
#if defined(SUBFIND)
      out.DM_Hsml = PS[index].Hsml;
#endif /* #if defined(SUBFIND) */
    }
  else
    {
      out.Distance  = 2.0e30;
      out.MinID     = 0;
      out.MinIDTask = -1;
#if defined(SUBFIND)
      out.DM_Hsml = 0;
#endif /* #if defined(SUBFIND) */
    }

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif /* #ifdef FOF */
