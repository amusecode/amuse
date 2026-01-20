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
 * \file        src/domain_toplevel.c
 * \date        05/2018
 * \brief       Top level tree construction and walk routines used for the
 *              domain decomposition.
 * \details     Uses BSD macros.
 *              contains functions:
 *                static int mydata_cmp(struct mydata *lhs, struct mydata *rhs)
 *                int domain_determineTopTree(void)
 *                void domain_do_local_refine(int n, int *list)
 *                void domain_walktoptree(int no)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 17.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../mesh/voronoi/voronoi.h"
#include "bsd_tree.h"
#include "domain.h"

/*! \brief Structure of tree nodes.
 */
struct mydata
{
  double workload;
  int topnode_index;

  RB_ENTRY(mydata) linkage; /* this creates the linkage pointers needed by the RB tree, using symbolic name 'linkage' */
};

/*! \brief Comparison function of tree elements.
 *
 *  Compares elements workload and topnode_index.
 *
 *  \param[in] lhs pointer to left hand side top level tree node.
 *  \param[in] rhs pointer to right hand side top level tree node.
 *
 *  \return -1: left is larger or lower topnode index, 1 opposite, 0 equal.
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

static double *list_cost, *list_sphcost;

/*! \brief Construct top-level tree.
 *
 *  This function constructs the global top-level tree node that is used
 *  for the domain decomposition. This is done by considering the string of
 *  Peano-Hilbert keys for all particles, which is recursively chopped off
 *  in pieces of eight segments until each segment holds at most a certain
 *  number of particles.
 *
 *  \return 0
 */
int domain_determineTopTree(void)
{
  double t0 = second();
  int count = 0, message_printed = 0;

  mp           = (struct domain_peano_hilbert_data *)mymalloc_movable(&mp, "mp", sizeof(struct domain_peano_hilbert_data) * NumPart);
  list_cost    = mymalloc_movable(&list_cost, "list_cost", sizeof(double) * NumPart);
  list_sphcost = mymalloc_movable(&list_sphcost, "listsph_cost", sizeof(double) * NumPart);

  for(int i = 0; i < NumPart; i++)
    {
      peano1D xb = domain_double_to_int(((P[i].Pos[0] - DomainCorner[0]) * DomainInverseLen) + 1.0);
      peano1D yb = domain_double_to_int(((P[i].Pos[1] - DomainCorner[1]) * DomainInverseLen) + 1.0);
      peano1D zb = domain_double_to_int(((P[i].Pos[2] - DomainCorner[2]) * DomainInverseLen) + 1.0);

      mp[count].key = Key[i] = peano_hilbert_key(xb, yb, zb, BITS_PER_DIMENSION);
      mp[count].index        = i;
      count++;

      list_cost[i]    = domain_grav_tot_costfactor(i);
      list_sphcost[i] = domain_hydro_tot_costfactor(i);
    }

  /* sort according to key (local particles!) */
  mysort_domain(mp, count, sizeof(struct domain_peano_hilbert_data));

  NTopnodes            = 1;
  NTopleaves           = 1;
  topNodes[0].Daughter = -1;
  topNodes[0].Parent   = -1;
  topNodes[0].Size     = PEANOCELLS;
  topNodes[0].StartKey = 0;
  topNodes[0].PIndex   = 0;
  topNodes[0].Count    = count;
  topNodes[0].Cost     = gravcost;
  topNodes[0].SphCost  = sphcost;

  int limitNTopNodes = 2 * imax(1 + (NTask / 7 + 1) * 8, All.TopNodeFactor * All.MultipleDomains * NTask);

#ifdef ADDBACKGROUNDGRID
  limitNTopNodes = imax(limitNTopNodes, 2 * All.GridSize * All.GridSize * All.GridSize);
#endif /* #ifdef ADDBACKGROUNDGRID */

  while(limitNTopNodes > MaxTopNodes)
    {
      mpi_printf("DOMAIN: Increasing TopNodeAllocFactor=%g  ", All.TopNodeAllocFactor);
      All.TopNodeAllocFactor *= 1.3;
      mpi_printf("new value=%g\n", All.TopNodeAllocFactor);
      if(All.TopNodeAllocFactor > 1000)
        terminate("something seems to be going seriously wrong here. Stopping.\n");

      MaxTopNodes = (int)(All.TopNodeAllocFactor * All.MaxPart + 1);

      topNodes        = (struct local_topnode_data *)myrealloc_movable(topNodes, (MaxTopNodes * sizeof(struct local_topnode_data)));
      TopNodes        = (struct topnode_data *)myrealloc_movable(TopNodes, (MaxTopNodes * sizeof(struct topnode_data)));
      DomainTask      = (int *)myrealloc_movable(DomainTask, (MaxTopNodes * sizeof(int)));
      DomainLeaveNode = (struct domain_cost_data *)myrealloc_movable(DomainLeaveNode, (MaxTopNodes * sizeof(struct domain_cost_data)));
    }

  RB_INIT(&queue_load);
  nload     = mymalloc("nload", limitNTopNodes * sizeof(struct mydata));
  int *list = mymalloc("list", limitNTopNodes * sizeof(int));

#ifdef ADDBACKGROUNDGRID
  peanokey MaxTopleaveSize = (PEANOCELLS / (All.GridSize * All.GridSize * All.GridSize));
#else  /* #ifdef ADDBACKGROUNDGRID */
  double limit = 1.0 / (All.TopNodeFactor * All.MultipleDomains * NTask);
#endif /* #ifdef ADDBACKGROUNDGRID #else */

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
          if(topNodes[nfirst->topnode_index].Size >= 8)
            {
              first_workload = nfirst->workload;
              break;
            }
        }

      for(struct mydata *np = RB_MIN(mytree, &queue_load); np != NULL; np = RB_NEXT(mytree, &queue_load, np))
        {
#ifndef ADDBACKGROUNDGRID
          if(np->workload < 0.125 * first_workload)
            break;

          if(NTopnodes + 8 * (count + 1) >= limitNTopNodes)
            break;
#endif /* #ifndef ADDBACKGROUNDGRID */

#ifdef ADDBACKGROUNDGRID
          if(topNodes[np->topnode_index].Size > MaxTopleaveSize)
#else  /* #ifdef ADDBACKGROUNDGRID */
          if(np->workload > limit || (NTopleaves < All.MultipleDomains * NTask && count == 0))
#endif /* #ifdef ADDBACKGROUNDGRID #else */
            {
              if(topNodes[np->topnode_index].Size < 8)
                {
                  if(message_printed == 0)
                    {
                      mpi_printf("DOMAIN: Note: we would like to refine top-tree, but PEANOGRID is not fine enough\n");
#ifndef OVERRIDE_PEANOGRID_WARNING
                      terminate(
                          "Consider setting BITS_PER_DIMENSION up to a value of 42 to get a fine enough PEANOGRID, or force a "
                          "continuation by activating OVERRIDE_PEANOGRID_WARNING");
#endif /* #ifndef OVERRIDE_PEANOGRID_WARNING */
                      message_printed = 1;
                    }
                }
              else
                {
                  list[count] = np->topnode_index;
                  count++;
                }
            }
        }

      if(count > 0)
        {
          domain_do_local_refine(count, list);
          iter++;
        }
    }
  while(count > 0);

  myfree(list);
  myfree(nload);
  myfree(list_sphcost);
  myfree(list_cost);
  myfree(mp);

  /* count the number of top leaves */
  NTopleaves = 0;
  domain_walktoptree(0);

  double t1 = second();
  mpi_printf("DOMAIN: NTopleaves=%d, determination of top-level tree involved %d iterations and took %g sec\n", NTopleaves, iter,
             timediff(t0, t1));

  t0 = second();

  domain_sumCost();

  t1 = second();
  mpi_printf("DOMAIN: cost summation for top-level tree took %g sec\n", timediff(t0, t1));

  return 0;
}

/*! \brief Refine top-level tree locally.
 *
 *  Requires arrays list_cost and list_sphcost, mp.
 *
 *  \param[in] n Number of nodes that should be refined.
 *  \param[in] list List of node indices that should be refined.
 *
 *  \return void
 */
void domain_do_local_refine(int n, int *list)
{
  double *worktotlist = mymalloc("worktotlist", 8 * n * sizeof(double));
  double *worklist    = mymalloc("worklist", 8 * n * sizeof(double));

  double non_zero = 0, non_zero_tot;

  /* create the new nodes */
  for(int k = 0; k < n; k++)
    {
      int i                = list[k];
      topNodes[i].Daughter = NTopnodes;
      NTopnodes += 8;
      NTopleaves += 7;

      for(int j = 0; j < 8; j++)
        {
          int sub = topNodes[i].Daughter + j;

          topNodes[sub].Daughter = -1;
          topNodes[sub].Parent   = i;
          topNodes[sub].Size     = (topNodes[i].Size >> 3);
          topNodes[sub].StartKey = topNodes[i].StartKey + j * topNodes[sub].Size;
          topNodes[sub].PIndex   = topNodes[i].PIndex;
          topNodes[sub].Count    = 0;
          topNodes[sub].Cost     = 0;
          topNodes[sub].SphCost  = 0;
        }

      int sub = topNodes[i].Daughter;

      for(int p = topNodes[i].PIndex, j = 0; p < topNodes[i].PIndex + topNodes[i].Count; p++)
        {
          if(j < 7)
            while(mp[p].key >= topNodes[sub + 1].StartKey)
              {
                j++;
                sub++;
                topNodes[sub].PIndex = p;
                if(j >= 7)
                  break;
              }

          topNodes[sub].Cost += list_cost[mp[p].index];
          topNodes[sub].SphCost += list_sphcost[mp[p].index];
          topNodes[sub].Count++;
        }

      for(int j = 0; j < 8; j++)
        {
          int sub             = topNodes[i].Daughter + j;
          worklist[k * 8 + j] = fac_work * topNodes[sub].Cost + fac_worksph * topNodes[sub].SphCost + fac_load * topNodes[sub].Count;

          if(worklist[k * 8 + j] != 0)
            non_zero++;
        }
    }

  MPI_Allreduce(&non_zero, &non_zero_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(non_zero_tot > 0.05 * (NTask * 8 * n))
    MPI_Allreduce(worklist, worktotlist, 8 * n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  else
    allreduce_sparse_double_sum(worklist, worktotlist, 8 * n);

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
          int sub = topNodes[i].Daughter + j;

          /* insert the  node */
          nload[sub].workload      = worktotlist[l];
          nload[sub].topnode_index = sub;
          RB_INSERT(mytree, &queue_load, &nload[sub]);
        }
    }

  myfree(worklist);
  myfree(worktotlist);
}

/*! \brief Walks top level tree recursively.
 *
 *  This function walks the global top tree in order to establish the
 *  number of leaves it has, and for assigning the leaf numbers along the
 *  Peano-Hilbert Curve. These leaves are later combined to domain pieces,
 *  which are distributed to different processors.
 *
 *  \param[in] no Present node.
 *
 *  \return void
 */
void domain_walktoptree(int no)
{
  if(topNodes[no].Daughter == -1)
    {
      topNodes[no].Leaf = NTopleaves;
      NTopleaves++;
    }
  else
    {
      for(int i = 0; i < 8; i++)
        domain_walktoptree(topNodes[no].Daughter + i);
    }
}
