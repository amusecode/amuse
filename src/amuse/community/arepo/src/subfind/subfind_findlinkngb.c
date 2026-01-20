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
 * \file        src/subfind/subfind_findlinkngb.c
 * \date        05/2018
 * \brief       Algorithm to find smoothing lengths of particles to get a
 *              desried number of neighbours.
 * \details     contains functions:
 *                static void particle2in(data_in * in, int i, int firstnode)
 *                static void out2particle(data_out * out, int i, int mode)
 *                static void kernel_local(void)
 *                static void kernel_imported(void)
 *                void subfind_find_linkngb(void)
 *                static int subfind_linkngb_evaluate(int target, int mode,
 *                  int threadid)
 *                int subfind_treefind_collective_export_node_threads(int no,
 *                  int i, int thread_id)
 *                static int subfind_ngb_compare_dist(const void *a, const
 *                  void *b)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 15.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef SUBFIND
#include "subfind.h"

static int subfind_ngb_compare_dist(const void *a, const void *b);
static int subfind_linkngb_evaluate(int target, int mode, int threadid);

static int *DM_NumNgb;
static double *Dist2list;
static int *Ngblist;
static MyFloat *Left, *Right;
static char *Todo;

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
  MyFloat DM_Hsml;

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

  in->DM_Hsml = PS[i].Hsml;

  in->Firstnode = firstnode;
}

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{
  int Ngb;
} data_out;

static data_out *DataResult, *DataOut;

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *  particle and cell data arrays.
 *  \param[in] i Index of particle.
 *  \param[in] mode Mode of function: local particles or information that was
 *  communicated from other tasks and has to be added locally?
 *
 *  \return void
 */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      DM_NumNgb[i] = out->Ngb;
    }
  else /* combine */
    {
      DM_NumNgb[i] += out->Ngb;
    }
}

#define USE_SUBCOMM_COMMUNICATOR
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

        if(Todo[i])
          subfind_linkngb_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        subfind_linkngb_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/*! \brief Iteratvie search for particle smoothing length to enclose a given
 *         number of neighbours.
 *
 *  \return void
 */
void subfind_find_linkngb(void)
{
  long long ntot;
  int i, npleft, iter = 0;
  double t0, t1;

  if(SubThisTask == 0)
    printf("SUBFIND-COLLECTIVE, root-task=%d: Start find_linkngb. (%d particles on root-task)\n", ThisTask, NumPartGroup);

  /* allocate buffers to arrange communication */

  Ngblist   = (int *)mymalloc("Ngblist", NumPartGroup * sizeof(int));
  Dist2list = (double *)mymalloc("Dist2list", NumPartGroup * sizeof(double));

  generic_set_MaxNexport();

  Left      = (MyFloat *)mymalloc("Left", sizeof(MyFloat) * NumPartGroup);
  Right     = (MyFloat *)mymalloc("Right", sizeof(MyFloat) * NumPartGroup);
  Todo      = (char *)mymalloc("Todo", sizeof(char) * NumPartGroup);
  DM_NumNgb = (int *)mymalloc_movable(&DM_NumNgb, "DM_NumNgb", sizeof(int) * NumPartGroup);

  for(i = 0; i < NumPartGroup; i++)
    {
      Left[i] = Right[i] = 0;
      Todo[i]            = 1;
    }

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      t0 = second();

      generic_comm_pattern(NumPartGroup, kernel_local, kernel_imported);

      /* do final operations on results */
      for(i = 0, npleft = 0; i < NumPartGroup; i++)
        {
          /* now check whether we had enough neighbours */
          if(Todo[i])
            {
              if(DM_NumNgb[i] != All.DesLinkNgb && ((Right[i] - Left[i]) > 1.0e-6 * Left[i] || Left[i] == 0 || Right[i] == 0))
                {
                  /* need to redo this particle */
                  npleft++;

                  if(DM_NumNgb[i] < All.DesLinkNgb)
                    Left[i] = dmax(PS[i].Hsml, Left[i]);
                  else
                    {
                      if(Right[i] != 0)
                        {
                          if(PS[i].Hsml < Right[i])
                            Right[i] = PS[i].Hsml;
                        }
                      else
                        Right[i] = PS[i].Hsml;
                    }

                  if(iter >= MAXITER - 10)
                    {
                      printf("i=%d task=%d ID=%d DM_Hsml=%g Left=%g Right=%g Right-Left=%g\n   pos=(%g|%g|%g)\n", i, ThisTask,
                             (int)P[i].ID, PS[i].Hsml, Left[i], Right[i], (double)(Right[i] - Left[i]), P[i].Pos[0], P[i].Pos[1],
                             P[i].Pos[2]);
                      fflush(stdout);
                    }

                  if(Right[i] > 0 && Left[i] > 0)
                    PS[i].Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
                  else
                    {
                      if(Right[i] == 0 && Left[i] == 0)
                        terminate("can't occur");

                      if(Right[i] == 0 && Left[i] > 0)
                        PS[i].Hsml *= 1.26;

                      if(Right[i] > 0 && Left[i] == 0)
                        PS[i].Hsml /= 1.26;
                    }
                }
              else
                Todo[i] = 0;
            }
        }

      sumup_large_ints_comm(1, &npleft, &ntot, SubComm);

      t1 = second();

      if(ntot > 0)
        {
          iter++;

          if(iter > 0 && SubThisTask == 0)
            {
              printf("SUBFIND-COLLECTIVE, root-task=%d: find linkngb iteration %d, need to repeat for %lld particles. (took %g sec)\n",
                     ThisTask, iter, ntot, timediff(t0, t1));
              fflush(stdout);
            }

          if(iter > MAXITER)
            terminate("failed to converge in neighbour iteration in density()\n");
        }
    }
  while(ntot > 0);

  myfree(DM_NumNgb);
  myfree(Todo);
  myfree(Right);
  myfree(Left);

  myfree(Dist2list);
  myfree(Ngblist);

  if(SubThisTask == 0)
    printf("SUBFIND-COLLECTIVE, root-task=%d: Done with find_linkngb\n", ThisTask);
}

/*! \brief Evaluate function for the neighbor search algorithm.
 *
 *  \param[in] target Index of particle of interest.
 *  \param[in] mode Local or imported particles?
 *  \param[in] treadid ID of thread.
 *
 *  \return 0
 */
static int subfind_linkngb_evaluate(int target, int mode, int threadid)
{
  int no, numnodes, *firstnode, numngb;
  double hsml;
  MyDouble *pos;
  int i, k, p, exported = 0;
  struct NODE *current;
  double dx, dy, dz, dist, r2;
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

  pos  = in->Pos;
  hsml = in->DM_Hsml;

  numngb = 0;

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
              p  = no;
              no = SubNextnode[no];

              dist = hsml;
              dx   = FOF_NEAREST_LONG_X(SubTree_Pos_list[3 * p + 0] - pos[0]);
              if(dx > dist)
                continue;
              dy = FOF_NEAREST_LONG_Y(SubTree_Pos_list[3 * p + 1] - pos[1]);
              if(dy > dist)
                continue;
              dz = FOF_NEAREST_LONG_Z(SubTree_Pos_list[3 * p + 2] - pos[2]);
              if(dz > dist)
                continue;
              if((r2 = (dx * dx + dy * dy + dz * dz)) > dist * dist)
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
                    break;
                }

              current = &SubNodes[no];

              no = current->u.d.sibling; /* in case the node can be discarded */

              dist = hsml + 0.5 * current->len;
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
          else
            { /* pseudo particle */
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES");

              if(target >= 0) /* if no target is given, export will not occur */
                {
                  exported = 1;

                  if(mode == MODE_LOCAL_PARTICLES)
                    subfind_treefind_collective_export_node_threads(no, target, threadid);
                }

              no = SubNextnode[no - SubTree_MaxNodes];
            }
        }
    }

  if(mode == MODE_LOCAL_PARTICLES) /* local particle */
    if(exported == 0)              /* completely local */
      if(numngb >= All.DesLinkNgb)
        {
          R2list = (r2type *)mymalloc("R2list", sizeof(r2type) * numngb);
          for(i = 0; i < numngb; i++)
            {
              R2list[i].index = Ngblist[i];
              R2list[i].r2    = Dist2list[i];
            }

          qsort(R2list, numngb, sizeof(r2type), subfind_ngb_compare_dist);

          PS[target].Hsml = sqrt(R2list[All.DesLinkNgb - 1].r2);
          numngb          = All.DesLinkNgb;

          for(i = 0; i < numngb; i++)
            {
              Ngblist[i]   = R2list[i].index;
              Dist2list[i] = R2list[i].r2;
            }

          myfree(R2list);
        }

  out.Ngb = numngb;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

/*! \brief Prepares node export.
 *
 *  \param[in] no Index of node.
 *  \param[in] i Index of particle.
 *  \param[in] thread_id Export thread.
 *
 *  \return 0
 */
int subfind_treefind_collective_export_node_threads(int no, int i, int thread_id)
{
  /* The task indicated by the pseudoparticle node */
  int task = SubDomainTask[no - (SubTree_MaxPart + SubTree_MaxNodes)];

  if(Thread[thread_id].Exportflag[task] != i)
    {
      Thread[thread_id].Exportflag[task]     = i;
      int nexp                               = Thread[thread_id].Nexport++;
      Thread[thread_id].PartList[nexp].Task  = task;
      Thread[thread_id].PartList[nexp].Index = i;
      Thread[thread_id].ExportSpace -= Thread[thread_id].ItemSize;
    }

  int nexp                      = Thread[thread_id].NexportNodes++;
  nexp                          = -1 - nexp;
  struct datanodelist *nodelist = (struct datanodelist *)(((char *)Thread[thread_id].PartList) + Thread[thread_id].InitialSpace);
  nodelist[nexp].Task           = task;
  nodelist[nexp].Index          = i;
  nodelist[nexp].Node           = SubDomainNodeIndex[no - (SubTree_MaxPart + SubTree_MaxNodes)];
  Thread[thread_id].ExportSpace -= sizeof(struct datanodelist) + sizeof(int);
  return 0;
}

/*! \brief Comparison function for r2type objects.
 *
 *  Compares element r2.
 *
 *  \param[in] a First object.
 *  \param[in] b Second object.
 *
 *  \return (-1,0,1) -1 if a < b.
 */
static int subfind_ngb_compare_dist(const void *a, const void *b)
{
  if(((r2type *)a)->r2 < (((r2type *)b)->r2))
    return -1;

  if(((r2type *)a)->r2 > (((r2type *)b)->r2))
    return +1;

  return 0;
}

#endif /* #ifdef SUBFIND */
