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
 * \file        src/fof/fof_findgroups.c
 * \date        05/2018
 * \brief       Routine to identify friend of friends groups.
 * \details     contains functions:
 *                static void particle2in(data_in * in, int i, int firstnode)
 *                static void out2particle(data_out * out, int i, int mode)
 *                static void kernel_local(void)
 *                static void kernel_imported(void)
 *                double fof_find_groups(MyIDType * vMinID, int *vHead,
 *                  int *vLen, int *vNext, int *vTail, int *vMinIDTask)
 *                static int fof_find_dmparticles_evaluate(int target,
 *                  int mode, int threadid)
 *                static int fof_treefind_fof_primary(MyDouble searchcenter[3],
 *                  MyFloat hsml, int target, int numnodes, int *firstnode,
 *                  int mode, int threadid)
 *                void fof_check_for_full_nodes_recursive(int no)
 *                int fof_return_a_particle_in_cell_recursive(int no)
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

static int fof_find_dmparticles_evaluate(int target, int mode, int threadid);
static int fof_treefind_fof_primary(MyDouble searchcenter[3], MyFloat hsml, int target, int numnodes, int *firstnode, int mode,
                                    int threadid);

static int *Tree_Head;

static MyIDType *MinID;
static int *Head, *Len, *Next, *Tail, *MinIDTask;

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];

  MyIDType MinID;
  int MinIDTask;

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/*! \brief Routine that fills the relevant particle/cell data into the input
 *         structure defined above. Needed by generic_comm_helpers2.
 *
 *  \param[out] in Data structure to fill.
 *  \param[in] i Index of particle in P array.
 *  \param[in] firstnode First note of communication.
 *
 *  \return void
 */
static void particle2in(data_in *in, int i, int firstnode)
{
  in->Pos[0] = P[i].Pos[0];
  in->Pos[1] = P[i].Pos[1];
  in->Pos[2] = P[i].Pos[2];

  in->MinID     = MinID[Head[i]];
  in->MinIDTask = MinIDTask[Head[i]];

  in->Firstnode = firstnode;
}

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{
  char link_count_flag;
} data_out;

static data_out *DataResult, *DataOut;

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *             particle and cell data arrays (P, SphP,...)
 *  \param[in] i Index of particle in P and SphP arrays
 *  \param[in] mode Mode of function: local particles or information that was
 *             communicated from other tasks and has to be added locally?
 *
 *  \return void
 */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      terminate("here not used");
    }
  else /* combine */
    {
      if(out->link_count_flag)
        Flags[i].Marked = 1;
    }
}

#include "../utils/generic_comm_helpers2.h"

static int link_across;
static int nprocessed;

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

        if(((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES)))
          {
            if(Flags[i].Nonlocal && Flags[i].Changed)
              {
                fof_find_dmparticles_evaluate(i, MODE_LOCAL_PARTICLES, threadid);

                nprocessed++;
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

        link_across += fof_find_dmparticles_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/*! \brief Links particles to groups.
 *
 *  \param[in, out] vMinID Pointer to MinID array.
 *  \param[in, out] vHead Pointer to Head array.
 *  \param[in, out] vLen Pointer to Len array.
 *  \param[in, out] vNext Pointer to Next array.
 *  \param[in, out] vTail Pointer to Tail array.
 *  \param[in, out] vMinIDTask Pointer to MinIDTask array.
 *
 *  \return Time spent in this function.
 */
double fof_find_groups(MyIDType *vMinID, int *vHead, int *vLen, int *vNext, int *vTail, int *vMinIDTask)
{
  MinID     = vMinID;
  Head      = vHead;
  Len       = vLen;
  Next      = vNext;
  Tail      = vTail;
  MinIDTask = vMinIDTask;

  int i, npart, marked;
  long long totmarked, totnpart;
  long long link_across_tot, ntot;
  double t0, t1, tstart, tend;

  tstart = second();

  mpi_printf("FOF: Start linking particles (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  /* allocate a flag field that is used to mark nodes that are fully inside the linking length */
  flag_node_inside_linkinglength = (unsigned char *)mymalloc("flag_node_inside_linkinglength", Tree_MaxNodes * sizeof(unsigned char));
  memset(flag_node_inside_linkinglength, 0, Tree_MaxNodes * sizeof(unsigned char));
  flag_node_inside_linkinglength -= Tree_MaxPart;

  Flags = (struct bit_flags *)mymalloc("Flags", NumPart * sizeof(struct bit_flags));

  generic_set_MaxNexport();

  Tree_Head = mymalloc("Tree_Head", Tree_NumNodes * sizeof(int));
  Tree_Head -= Tree_MaxPart;

  /* allocate buffers to arrange communication */
  generic_alloc_partlist_nodelist_ngblist_threadbufs();

  t0 = second();

  /* first, link only among local particles */
  for(i = 0, marked = 0, npart = 0; i < NumPart; i++)
    {
      if(((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES)))
        {
          fof_find_dmparticles_evaluate(i, MODE_LOCAL_NO_EXPORT, 0);

          npart++;

          if(Flags[i].Nonlocal)
            marked++;
        }
    }

  sumup_large_ints(1, &marked, &totmarked);
  sumup_large_ints(1, &npart, &totnpart);
  t1 = second();
  mpi_printf("FOF: links on local processor done (took %g sec).\nFOF: Marked=%lld out of the %lld primaries which are linked\n",
             timediff(t0, t1), totmarked, totnpart);

  generic_free_partlist_nodelist_ngblist_threadbufs();

  t0 = second();
  fof_check_for_full_nodes_recursive(Tree_MaxPart);
  t1 = second();
  mpi_printf("FOF: fully linked nodes determined (took %g sec).\n", timediff(t0, t1));
  mpi_printf("FOF: begin linking across processors (presently allocated=%g MB) \n", AllocatedBytes / (1024.0 * 1024.0));

  for(i = 0; i < NumPart; i++)
    Flags[i].Marked = 1;

  do
    {
      t0 = second();

      for(i = 0; i < NumPart; i++)
        {
          Flags[i].Changed      = Flags[i].Marked;
          Flags[i].Marked       = 0;
          Flags[i].MinIDChanged = 0;
        }

      NextParticle = 0; /* begin with this index */

      link_across = 0;
      nprocessed  = 0;

      generic_comm_pattern(NumPart, kernel_local, kernel_imported);

      sumup_large_ints(1, &link_across, &link_across_tot);
      sumup_large_ints(1, &nprocessed, &ntot);

      t1 = second();

      mpi_printf("FOF: have done %15lld cross links (processed %14lld, took %g sec)\n", link_across_tot, ntot, timediff(t0, t1));

      /* let's check out which particles have changed their MinID */
      for(i = 0; i < NumPart; i++)
        if(Flags[i].Nonlocal)
          {
            if(Flags[Head[i]].MinIDChanged)
              Flags[i].Marked = 1;
          }
    }
  while(link_across_tot > 0);

  Tree_Head += Tree_MaxPart;
  myfree(Tree_Head);
  myfree(Flags);
  /* free flag */
  myfree(flag_node_inside_linkinglength + Tree_MaxPart);

  mpi_printf("FOF: Local groups found.\n");

  tend = second();
  return timediff(tstart, tend);
}

/*! \brief Links dark matter particles.
 *
 *  \param[in] target Index of particle/cell.
 *  \param[in] mode Flag if it operates on local or imported data.
 *  \param[in] threadid ID of thread.
 *
 *  \return Number of links.
 */
static int fof_find_dmparticles_evaluate(int target, int mode, int threadid)
{
  int j, n, links, p, s, ss, numnodes, *firstnode;
  int numngb;
  MyDouble *pos;
  data_in local, *target_data;

  links = 0;

  if(mode == MODE_LOCAL_NO_EXPORT || mode == MODE_LOCAL_PARTICLES)
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

  numngb = fof_treefind_fof_primary(pos, LinkL, target, numnodes, firstnode, mode, threadid);

  if(mode == MODE_LOCAL_PARTICLES || mode == MODE_LOCAL_NO_EXPORT)
    for(n = 0; n < numngb; n++)
      {
        j = Thread[threadid].Ngblist[n];

        if(Head[target] != Head[j]) /* only if not yet linked */
          {
            if(Len[Head[target]] > Len[Head[j]]) /* p group is longer */
              {
                p = target;
                s = j;
              }
            else
              {
                p = j;
                s = target;
              }
            Next[Tail[Head[p]]] = Head[s];

            Tail[Head[p]] = Tail[Head[s]];

            Len[Head[p]] += Len[Head[s]];

            if(MinID[Head[s]] < MinID[Head[p]])
              {
                MinID[Head[p]]     = MinID[Head[s]];
                MinIDTask[Head[p]] = MinIDTask[Head[s]];
              }

            ss = Head[s];
            do
              Head[ss] = Head[p];
            while((ss = Next[ss]) >= 0);
          }
      }

  if(mode == MODE_IMPORTED_PARTICLES)
    {
      if(numngb > 0)
        DataResult[target].link_count_flag = 1;
      else
        DataResult[target].link_count_flag = 0;
    }

  links += numngb;

  return links;
}

/*! \brief Finds the neighbors among the primary link types which are within a
 *         certain distance.
 *
 *  \param[in] searchcenter Position of search center.
 *  \param[in] hsml Search radius.
 *  \param[in] target Index of partcle.
 *  \param[in] numnodes Number of nodes.
 *  \param[in] fistnode First node.
 *  \param[in] mode
 *             -1: only local particles should be found and no export occurs;
 *              0: export occurs, but local particles are ignored;
 *              1: particles are found for an imported point.
 *  \param[in] threadid ID of thread.
 *
 *  \return Number of particles found.
 */
static int fof_treefind_fof_primary(MyDouble searchcenter[3], MyFloat hsml, int target, int numnodes, int *firstnode, int mode,
                                    int threadid)
{
  int k, numngb, no, p, nexport_flag = 0;
  MyDouble dx, dy, dz, dist, r2;

#define FACT2 0.866025403785 /* sqrt(3)/2 */
#define FACT3 (2.0 * FACT2)  /* sqrt(3)   */

  MyDouble xtmp, ytmp, ztmp;

  numngb = 0;

  for(k = 0; k < numnodes; k++)
    {
      if(mode == MODE_LOCAL_PARTICLES || mode == MODE_LOCAL_NO_EXPORT)
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
              p  = no;
              no = Nextnode[no];

              if(!((1 << P[p].Type) & (FOF_PRIMARY_LINK_TYPES)))
                continue;

              if(mode == MODE_LOCAL_PARTICLES)
                continue;

              dist = hsml;
              dx   = FOF_NEAREST_LONG_X(Tree_Pos_list[3 * p + 0] - searchcenter[0]);
              if(dx > dist)
                continue;
              dy = FOF_NEAREST_LONG_Y(Tree_Pos_list[3 * p + 1] - searchcenter[1]);
              if(dy > dist)
                continue;
              dz = FOF_NEAREST_LONG_Z(Tree_Pos_list[3 * p + 2] - searchcenter[2]);
              if(dz > dist)
                continue;
              if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;

              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(MinID[Head[p]] > DataGet[target].MinID)
                    {
                      MinID[Head[p]]              = DataGet[target].MinID;
                      MinIDTask[Head[p]]          = DataGet[target].MinIDTask;
                      Flags[Head[p]].MinIDChanged = 1;
                      numngb++;
                    }
                }
              else
                {
                  /* this will only be done for MODE_LOCAL_NO_EXPORT */
                  Thread[threadid].Ngblist[numngb++] = p;
                }
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes) /* internal node */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no <
                     Tree_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the branch */
                    break;

                  if(Tree_Head[no] >= 0)
                    if(MinID[Tree_Head[no]] <= DataGet[target].MinID)
                      {
                        no = Nodes[no].u.d.sibling; /* the node can be discarded */
                        continue;
                      }
                }

              struct NODE *current = &Nodes[no];
              int nocur            = no;
              no                   = current->u.d.sibling; /* in case the node can be discarded */

              if(mode == MODE_LOCAL_PARTICLES)
                {
                  if(nocur >= Tree_FirstNonTopLevelNode)
                    {
                      /* we have a node with only local particles, hence we can skip it for mode == 0 */
                      continue;
                    }
                }

              dist = hsml + 0.5 * current->len;
              dx   = FOF_NEAREST_LONG_X(current->center[0] - searchcenter[0]);
              if(dx > dist)
                continue;
              dy = FOF_NEAREST_LONG_Y(current->center[1] - searchcenter[1]);
              if(dy > dist)
                continue;
              dz = FOF_NEAREST_LONG_Z(current->center[2] - searchcenter[2]);
              if(dz > dist)
                continue;

              /* now test against the minimal sphere enclosing everything */
              dist += FACT1 * current->len;
              r2 = dx * dx + dy * dy + dz * dz;
              if(r2 > dist * dist)
                continue;

              if(mode != MODE_LOCAL_PARTICLES)
                {
                  /* test whether the node is contained within the sphere */
                  dist = hsml - FACT2 * current->len;
                  if(dist > 0)
                    if(r2 < dist * dist && hsml > FACT3 * current->len)
                      {
                        if(flag_node_inside_linkinglength[nocur] & (1 << BITFLAG_INSIDE_LINKINGLENGTH)) /* already flagged */
                          {
                            /* sufficient to return only one particle inside this cell */
                            p = fof_return_a_particle_in_cell_recursive(nocur);

                            if(p >= 0)
                              {
                                if(mode == MODE_IMPORTED_PARTICLES)
                                  {
                                    if(MinID[Head[p]] > DataGet[target].MinID)
                                      {
                                        MinID[Head[p]]              = DataGet[target].MinID;
                                        MinIDTask[Head[p]]          = DataGet[target].MinIDTask;
                                        Flags[Head[p]].MinIDChanged = 1;
                                        numngb++;
                                      }
                                  }
                                else
                                  Thread[threadid].Ngblist[numngb++] = p;
                              }

                            continue;
                          }
                        else
                          {
                            /* flag it now */
                            flag_node_inside_linkinglength[nocur] |= (1 << BITFLAG_INSIDE_LINKINGLENGTH);
                          }
                      }
                }

              no = current->u.d.nextnode; /* ok, we need to open the node */
            }
          else if(no >= Tree_ImportedNodeOffset) /* point from imported nodelist */
            {
              terminate("do not expect imported points here");
            }
          else
            {
              if(mode == MODE_LOCAL_PARTICLES)
                {
                  if(target >= 0)
                    tree_treefind_export_node_threads(no, target, threadid);
                }
              else if(mode == MODE_LOCAL_NO_EXPORT)
                {
                  nexport_flag = 1;
                }
              else if(mode == MODE_IMPORTED_PARTICLES)
                terminate("stop no=%d Tree_MaxPart=%d Tree_MaxNodes=%d", no, Tree_MaxPart, Tree_MaxNodes);

              no = Nextnode[no - Tree_MaxNodes];
              continue;
            }
        }
    }

  if(mode == MODE_LOCAL_NO_EXPORT)
    {
      if(nexport_flag == 0)
        Flags[target].Nonlocal = 0;
      else
        Flags[target].Nonlocal = 1;
    }

  return numngb;
}

/*! \brief Walks a tree recursively and sets Tree_Head of node.
 *
 *  \param[in] no Index of node we are in.
 *
 *  \return void
 */
void fof_check_for_full_nodes_recursive(int no)
{
  if(no >= Tree_MaxPart && no < Tree_MaxPart + Tree_MaxNodes) /* internal node */
    {
      int head = -1; /* no particle yet */

      int p = Nodes[no].u.d.nextnode;

      while(p != Nodes[no].u.d.sibling)
        {
          if(p < Tree_MaxPart) /* a particle */
            {
              if((1 << P[p].Type) & (FOF_PRIMARY_LINK_TYPES))
                {
                  if(head == -1)
                    head = Head[p];
                  else if(head >= 0)
                    {
                      if(head != Head[p])
                        head = -2;
                    }
                }

              p = Nextnode[p];
            }
          else if(p < Tree_MaxPart + Tree_MaxNodes) /* an internal node  */
            {
              fof_check_for_full_nodes_recursive(p);

              if(head == -1)
                head = Tree_Head[p];
              else if(head >= 0)
                {
                  if(head != Tree_Head[p])
                    head = -2;
                }

              p = Nodes[p].u.d.sibling;
            }
          else /* a pseudo particle */
            p = Nextnode[p - Tree_MaxNodes];
        }

      Tree_Head[no] = head;
    }
}

/*! \brief Finds a particle in node.
 *
 *  \param[in] no Index of node.
 *
 *  \return Particle index; -1 if no particle was found.
 */
int fof_return_a_particle_in_cell_recursive(int no)
{
  if(no >= Tree_MaxPart && no < Tree_MaxPart + Tree_MaxNodes) /* internal node */
    {
      int p = Nodes[no].u.d.nextnode;

      while(p != Nodes[no].u.d.sibling)
        {
          if(p < Tree_MaxPart) /* a particle */
            {
              if((1 << P[p].Type) & (FOF_PRIMARY_LINK_TYPES))
                {
                  return p;
                }

              p = Nextnode[p];
            }
          else if(p < Tree_MaxPart + Tree_MaxNodes) /* an internal node  */
            {
              int ret = fof_return_a_particle_in_cell_recursive(p);

              if(ret >= 0)
                return ret;

              p = Nodes[p].u.d.sibling;
            }
          else /* a pseudo particle */
            p = Nextnode[p - Tree_MaxNodes];
        }
    }

  return -1;
}

#endif /* #ifdef FOF */
