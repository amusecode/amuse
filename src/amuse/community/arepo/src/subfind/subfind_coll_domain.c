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
 * \file        src/subfind/subfind_coll_domain.c
 * \date        05/2018
 * \brief       Domain decomposition for collective subfind algorithm.
 * \details     contains functions:
 *                static int mydata_cmp(struct mydata *lhs, struct mydata *rhs)
 *                void subfind_coll_domain_decomposition(void)
 *                void subfind_coll_findExtent(void)
 *                int subfind_coll_domain_determineTopTree(void)
 *                void subfind_domain_do_local_refine(int n, int *list)
 *                void subfind_coll_domain_walktoptree(int no)
 *                void subfind_coll_domain_combine_topleaves_to_domains(int ncpu, int ndomain)
 *                void subfind_coll_domain_allocate(void)
 *                void subfind_coll_domain_free(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 15.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef SUBFIND
#include "../domain/bsd_tree.h"
#include "../domain/domain.h"
#include "subfind.h"

/*! \brief Define structure of my tree nodes.
 */
struct mydata
{
  double workload;
  int topnode_index;

  RB_ENTRY(mydata) linkage; /* this creates the linkage pointers needed by the RB tree, using symbolic name 'linkage' */
};

/*! \brief Comparison function of mydata objects (i.e. tree elements).
 *
 *  Compares the elements (most important first):
 *   workload, topnode_index.
 *
 *  \param[in] lhs First object to compare.
 *  \param[in] rhs Second object to compare.
 *
 *  \return (-1,0,1) -1 if lhs.workload > rhs.workload or lhs.topnode_index <
 *          rhs.topnode_index.
 */
static int mydata_cmp(struct mydata *lhs, struct mydata *rhs)
{
  if(lhs->workload > rhs->workload)
    return -1;
  else if(lhs->workload < rhs->workload)
    return 1;
  else if(lhs->topnode_index < rhs->topnode_index)
    return -1;
  else if(lhs->topnode_index > rhs->topnode_index)
    return 1;

  return 0;
}

/* the following macro declares 'struct mytree', which is the header element
 * needed as handle for a tree
 */
RB_HEAD(mytree, mydata);

static struct mydata *nload;
static struct mytree queue_load;

/* the following macros declare appropriate function prototypes and functions
 * needed for this type of tree
 */
RB_PROTOTYPE_STATIC(mytree, mydata, linkage, mydata_cmp);
RB_GENERATE_STATIC(mytree, mydata, linkage, mydata_cmp);

/*! \brief Performs domain decomposition for subfind collective.
 *
 *  \return void
 */
void subfind_coll_domain_decomposition(void)
{
  int i;
  int col_grouplen, col_partcount;

  subfind_coll_domain_allocate();
  subfind_coll_findExtent();

  Key             = (peanokey *)mymalloc_movable(&Key, "Key", (sizeof(peanokey) * NumPart));
  Sub_LocTopNodes = (struct local_topnode_data *)mymalloc_movable(&Sub_LocTopNodes, "Sub_LocTopNodes",
                                                                  (MaxTopNodes * sizeof(struct local_topnode_data)));

  MPI_Allreduce(&NumPartGroup, &col_grouplen, 1, MPI_INT, MPI_SUM, SubComm);
  MPI_Allreduce(&NumPart, &col_partcount, 1, MPI_INT, MPI_SUM, SubComm);

  fac_work = 0.5 / col_grouplen;
  fac_load = 0.5 / col_partcount;

  subfind_coll_domain_determineTopTree();

  /* find the split of the top-level tree */
  subfind_coll_domain_combine_topleaves_to_domains(SubNTask, SubNTopleaves);

  /* determine the particles that need to be exported, and to which CPU they need to be sent */
  for(i = 0; i < NumPart; i++)
    {
      if(PS[i].GrNr == GrNr)
        {
          int no = 0;
          while(Sub_LocTopNodes[no].Daughter >= 0)
            no = Sub_LocTopNodes[no].Daughter + (Key[i] - Sub_LocTopNodes[no].StartKey) / (Sub_LocTopNodes[no].Size >> 3);

          no = Sub_LocTopNodes[no].Leaf;

          int task = SubDomainTask[no];

          PS[i].TargetTask = task;
        }
      else
        PS[i].TargetTask = SubThisTask;

      PS[i].TargetIndex = 0; /* unimportant here */
    }

  fof_subfind_exchange(SubComm);

  /* note that the domain decomposition leads to an invalid values of NumPartGroup. This will however be redetermined in the main
   * routine of the collective subfind, after the domain decomposition has been done.
   */

  /* copy what we need for the topnodes */
  for(i = 0; i < SubNTopnodes; i++)
    {
      SubTopNodes[i].StartKey = Sub_LocTopNodes[i].StartKey;
      SubTopNodes[i].Size     = Sub_LocTopNodes[i].Size;
      SubTopNodes[i].Daughter = Sub_LocTopNodes[i].Daughter;
      SubTopNodes[i].Leaf     = Sub_LocTopNodes[i].Leaf;

      int j;
      int bits   = my_ffsll(SubTopNodes[i].Size);
      int blocks = (bits - 1) / 3 - 1;

      for(j = 0; j < 8; j++)
        {
          peano1D xb, yb, zb;
          peano_hilbert_key_inverse(SubTopNodes[i].StartKey + j * (SubTopNodes[i].Size >> 3), BITS_PER_DIMENSION, &xb, &yb, &zb);
          xb >>= blocks;
          yb >>= blocks;
          zb >>= blocks;
          int idx = (xb & 1) | ((yb & 1) << 1) | ((zb & 1) << 2);
          if(idx < 0 || idx > 7)
            terminate("j=%d  idx=%d", j, idx);

          SubTopNodes[i].MortonToPeanoSubnode[idx] = j;
        }
    }

  myfree(Sub_LocTopNodes);
  myfree(Key);

  SubTopNodes   = (struct topnode_data *)myrealloc_movable(SubTopNodes, SubNTopnodes * sizeof(struct topnode_data));
  SubDomainTask = (int *)myrealloc_movable(SubDomainTask, SubNTopleaves * sizeof(int));
}

/*! \brief Determines extent of local data and writes it to global variables.
 *
 *  \return void
 */
void subfind_coll_findExtent(void)
{
  int i, j;
  double len, xmin[3], xmax[3], xmin_glob[3], xmax_glob[3];

  /* determine extension */
  for(i = 0; i < 3; i++)
    {
      xmin[i] = MAX_REAL_NUMBER;
      xmax[i] = -MAX_REAL_NUMBER;
    }

  for(i = 0; i < NumPart; i++)
    {
      if(PS[i].GrNr == GrNr)
        {
          for(j = 0; j < 3; j++)
            {
#ifdef CELL_CENTER_GRAVITY
              if(P[i].Type == 0)
                {
                  if(xmin[j] > PS[i].Center[j])
                    xmin[j] = PS[i].Center[j];

                  if(xmax[j] < PS[i].Center[j])
                    xmax[j] = PS[i].Center[j];
                }
              else
#endif /* #ifdef CELL_CENTER_GRAVITY */
                {
                  if(xmin[j] > P[i].Pos[j])
                    xmin[j] = P[i].Pos[j];

                  if(xmax[j] < P[i].Pos[j])
                    xmax[j] = P[i].Pos[j];
                }
            }
        }
    }

  MPI_Allreduce(xmin, xmin_glob, 3, MPI_DOUBLE, MPI_MIN, SubComm);
  MPI_Allreduce(xmax, xmax_glob, 3, MPI_DOUBLE, MPI_MAX, SubComm);

  len = 0;
  for(j = 0; j < 3; j++)
    if(xmax_glob[j] - xmin_glob[j] > len)
      len = xmax_glob[j] - xmin_glob[j];

  len *= 1.001;

  SubDomainLen        = len;
  SubDomainInverseLen = 1.0 / SubDomainLen;
  SubDomainFac        = 1.0 / len * (((peanokey)1) << (BITS_PER_DIMENSION));
  SubDomainBigFac     = (SubDomainLen / (((long long)1) << 52));

  for(j = 0; j < 3; j++)
    {
      SubDomainCenter[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]);
      SubDomainCorner[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]) - 0.5 * len;
    }
}

/*! \brief Determines extent of the subfind top-tree.
 *
 *  \return void
 */
int subfind_coll_domain_determineTopTree(void)
{
  int i, count;

  mp = (struct domain_peano_hilbert_data *)mymalloc("mp", sizeof(struct domain_peano_hilbert_data) * NumPartGroup);

  for(i = 0, count = 0; i < NumPart; i++)
    {
      if(PS[i].GrNr == GrNr)
        {
          peano1D xb, yb, zb;

#ifdef CELL_CENTER_GRAVITY
          if(P[i].Type == 0)
            {
              xb = domain_double_to_int(((PS[i].Center[0] - SubDomainCorner[0]) * SubDomainInverseLen) + 1.0);
              yb = domain_double_to_int(((PS[i].Center[1] - SubDomainCorner[1]) * SubDomainInverseLen) + 1.0);
              zb = domain_double_to_int(((PS[i].Center[2] - SubDomainCorner[2]) * SubDomainInverseLen) + 1.0);
            }
          else
#endif /* #ifdef CELL_CENTER_GRAVITY */
            {
              xb = domain_double_to_int(((P[i].Pos[0] - SubDomainCorner[0]) * SubDomainInverseLen) + 1.0);
              yb = domain_double_to_int(((P[i].Pos[1] - SubDomainCorner[1]) * SubDomainInverseLen) + 1.0);
              zb = domain_double_to_int(((P[i].Pos[2] - SubDomainCorner[2]) * SubDomainInverseLen) + 1.0);
            }

          mp[count].key = Key[i] = peano_hilbert_key(xb, yb, zb, BITS_PER_DIMENSION);
          mp[count].index        = i;
          count++;
        }
    }

  if(count != NumPartGroup)
    terminate("cost != NumPartGroup");

  mysort_domain(mp, count, sizeof(struct domain_peano_hilbert_data));

  SubNTopnodes                = 1;
  SubNTopleaves               = 1;
  Sub_LocTopNodes[0].Daughter = -1;
  Sub_LocTopNodes[0].Parent   = -1;
  Sub_LocTopNodes[0].Size     = PEANOCELLS;
  Sub_LocTopNodes[0].StartKey = 0;
  Sub_LocTopNodes[0].PIndex   = 0;
  Sub_LocTopNodes[0].Cost     = NumPartGroup;
  Sub_LocTopNodes[0].Count    = NumPartGroup;

  int limitNTopNodes = 2 * imax(1 + (NTask / 7 + 1) * 8, All.TopNodeFactor * SubNTask);

  if(limitNTopNodes > MaxTopNodes)
    terminate("limitNTopNodes > MaxTopNodes");

  RB_INIT(&queue_load);
  nload     = mymalloc("nload", limitNTopNodes * sizeof(struct mydata));
  int *list = mymalloc("list", limitNTopNodes * sizeof(int));

  double limit = 1.0 / (All.TopNodeFactor * SubNTask);

  /* insert the root node */
  nload[0].workload      = 1.0;
  nload[0].topnode_index = 0;
  RB_INSERT(mytree, &queue_load, &nload[0]);

  int iter = 0;

  do
    {
      count = 0;

      double first_workload = 0;

      for(struct mydata *nfirst = RB_MIN(mytree, &queue_load); nfirst != NULL; nfirst = RB_NEXT(mytree, &queue_load, nfirst))
        {
          if(Sub_LocTopNodes[nfirst->topnode_index].Size >= 8)
            {
              first_workload = nfirst->workload;
              break;
            }
        }

      for(struct mydata *np = RB_MIN(mytree, &queue_load); np != NULL; np = RB_NEXT(mytree, &queue_load, np))
        {
          if(np->workload < 0.125 * first_workload)
            break;

          if(SubNTopnodes + 8 * (count + 1) >= limitNTopNodes)
            break;

          if(np->workload > limit || (SubNTopleaves < SubNTask && count == 0))
            {
              if(Sub_LocTopNodes[np->topnode_index].Size >= 8)
                {
                  list[count] = np->topnode_index;
                  count++;
                }
            }
        }

      if(count > 0)
        {
          subfind_domain_do_local_refine(count, list);
          iter++;
        }
    }
  while(count > 0);

  myfree(list);
  myfree(nload);
  myfree(mp);

  /* count toplevel leaves */

  /* count the number of top leaves */
  SubNTopleaves = 0;
  subfind_coll_domain_walktoptree(0);

  if(SubNTopleaves < SubNTask)
    terminate("SubNTopleaves = %d < SubNTask = %d", SubNTopleaves, SubNTask);

  return 0;
}

/*! \brief Refines top-tree locally.
 *
 *  \param[in] n Number of new nodes.
 *  \param[in] list Array with indices of new nodes.
 *
 *  \return void
 */
void subfind_domain_do_local_refine(int n, int *list)
{
  double *worktotlist = mymalloc("worktotlist", 8 * n * sizeof(double));
  double *worklist    = mymalloc("worklist", 8 * n * sizeof(double));

  /* create the new nodes */
  for(int k = 0; k < n; k++)
    {
      int i = list[k];

      Sub_LocTopNodes[i].Daughter = SubNTopnodes;
      SubNTopnodes += 8;
      SubNTopleaves += 7;

      for(int j = 0; j < 8; j++)
        {
          int sub = Sub_LocTopNodes[i].Daughter + j;

          Sub_LocTopNodes[sub].Daughter = -1;
          Sub_LocTopNodes[sub].Parent   = i;
          Sub_LocTopNodes[sub].Size     = (Sub_LocTopNodes[i].Size >> 3);
          Sub_LocTopNodes[sub].StartKey = Sub_LocTopNodes[i].StartKey + j * Sub_LocTopNodes[sub].Size;
          Sub_LocTopNodes[sub].PIndex   = Sub_LocTopNodes[i].PIndex;
          Sub_LocTopNodes[sub].Cost     = 0;
          Sub_LocTopNodes[sub].Count    = 0;
        }

      int sub = Sub_LocTopNodes[i].Daughter;

      for(int p = Sub_LocTopNodes[i].PIndex, j = 0; p < Sub_LocTopNodes[i].PIndex + Sub_LocTopNodes[i].Count; p++)
        {
          if(PS[mp[p].index].GrNr != GrNr)
            terminate("Houston, we have a problem.");

          if(j < 7)
            while(mp[p].key >= Sub_LocTopNodes[sub + 1].StartKey)
              {
                j++;
                sub++;
                Sub_LocTopNodes[sub].PIndex = p;
                if(j >= 7)
                  break;
              }

          Sub_LocTopNodes[sub].Count++;
          Sub_LocTopNodes[sub].Cost++;
        }

      for(int j = 0; j < 8; j++)
        {
          sub                 = Sub_LocTopNodes[i].Daughter + j;
          worklist[k * 8 + j] = fac_work * Sub_LocTopNodes[sub].Cost + fac_load * Sub_LocTopNodes[sub].Count;
        }
    }

  MPI_Allreduce(worklist, worktotlist, 8 * n, MPI_DOUBLE, MPI_SUM, SubComm);

  for(int k = 0; k < n; k++)
    {
      int i = list[k];
      RB_REMOVE(mytree, &queue_load, &nload[i]);
    }

  for(int k = 0, l = 0; k < n; k++)
    {
      int i = list[k];

      for(int j = 0; j < 8; j++, l++)
        {
          int sub = Sub_LocTopNodes[i].Daughter + j;

          /* insert the  node */
          nload[sub].workload      = worktotlist[l];
          nload[sub].topnode_index = sub;
          RB_INSERT(mytree, &queue_load, &nload[sub]);
        }
    }

  myfree(worklist);
  myfree(worktotlist);
}

/*! \brief Walk the top tree and set reference to leaf node.
 *
 *  \param[in] no Node index.
 *
 *  \return void
 */
void subfind_coll_domain_walktoptree(int no)
{
  int i;

  if(Sub_LocTopNodes[no].Daughter == -1)
    {
      Sub_LocTopNodes[no].Leaf = SubNTopleaves;
      SubNTopleaves++;
    }
  else
    {
      for(i = 0; i < 8; i++)
        subfind_coll_domain_walktoptree(Sub_LocTopNodes[no].Daughter + i);
    }
}

/*! \brief Uses the cumulative cost function (which weights work-load and
 *         memory-load equally) to subdivide the list of top-level leave
 *         nodes into pieces that are (approximately) equal in size.
 *
 *  \param[in] ncpu Number of tasks.
 *  \param[in] ndomain Number of domains.
 *
 *  \return void
 */
void subfind_coll_domain_combine_topleaves_to_domains(int ncpu, int ndomain)
{
  int i, j, start, end, n, no;
  double work, workavg, work_before, workavg_before, workhalfnode;
  float *domainWork, *local_domainWork;
  int *domainCount, *local_domainCount;

  /* sum the costs for each top leave */

  domainWork  = (float *)mymalloc("local_domainWork", SubNTopleaves * sizeof(float));
  domainCount = (int *)mymalloc("local_domainCount", SubNTopleaves * sizeof(int));

  local_domainWork  = (float *)mymalloc("local_domainWork", SubNTopleaves * sizeof(float));
  local_domainCount = (int *)mymalloc("local_domainCount", SubNTopleaves * sizeof(int));

  for(i = 0; i < SubNTopleaves; i++)
    {
      local_domainWork[i]  = 0;
      local_domainCount[i] = 0;
    }

  /* find for each particle its top-leave, and then add the associated cost with it */
  for(n = 0; n < NumPart; n++)
    {
      if(PS[n].GrNr == GrNr)
        {
          no = 0;
          while(Sub_LocTopNodes[no].Daughter >= 0)
            no = Sub_LocTopNodes[no].Daughter + (Key[n] - Sub_LocTopNodes[no].StartKey) / (Sub_LocTopNodes[no].Size >> 3);

          no = Sub_LocTopNodes[no].Leaf;

          local_domainCount[no] += 1;
          local_domainWork[no] += 1;
        }
    }

  MPI_Allreduce(local_domainWork, domainWork, SubNTopleaves, MPI_FLOAT, MPI_SUM, SubComm);
  MPI_Allreduce(local_domainCount, domainCount, SubNTopleaves, MPI_INT, MPI_SUM, SubComm);

  myfree(local_domainCount);
  myfree(local_domainWork);

  /* now combine the top leaves to form the individual domains */

  workhalfnode = 0.5 / ndomain;
  workavg      = 1.0 / ncpu;
  work_before = workavg_before = 0;

  start = 0;

  for(i = 0; i < ncpu; i++)
    {
      work = 0;
      end  = start;

      work += fac_work * domainWork[end] + fac_load * domainCount[end];

      while((work + work_before + (end + 1 < ndomain ? fac_work * domainWork[end + 1] + fac_load * domainCount[end + 1] : 0) <
             workavg + workavg_before + workhalfnode) ||
            (i == ncpu - 1 && end < ndomain - 1))
        {
          if((ndomain - end) > (ncpu - i))
            end++;
          else
            break;

          work += fac_work * domainWork[end] + fac_load * domainCount[end];
        }

      for(j = start; j <= end; j++)
        SubDomainTask[j] = i;

      work_before += work;
      workavg_before += workavg;
      start = end + 1;
    }

  myfree(domainCount);
  myfree(domainWork);
}

/*! \brief Allocates all the stuff that will be required for the
 *         tree-construction/walk later on.
 *
 *  \return void
 */
void subfind_coll_domain_allocate(void)
{
  MaxTopNodes = (int)(All.TopNodeAllocFactor * All.MaxPart + 1);

  if(SubDomainTask)
    terminate("subfind collective domain storage already allocated");

  SubTopNodes   = (struct topnode_data *)mymalloc_movable(&SubTopNodes, "SubTopNodes", (MaxTopNodes * sizeof(struct topnode_data)));
  SubDomainTask = (int *)mymalloc_movable(&SubDomainTask, "SubDomainTask", (MaxTopNodes * sizeof(int)));
}

/*! \brief Free memory used for subfind collective domain decomposition.
 *
 *  \return void
 */
void subfind_coll_domain_free(void)
{
  if(!SubDomainTask)
    terminate("subfind collective domain storage not allocated");

  myfree(SubDomainTask);
  myfree(SubTopNodes);

  SubDomainTask = NULL;
  SubTopNodes   = NULL;
}

#endif /* #ifdef SUBFIND */
