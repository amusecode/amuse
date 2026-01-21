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
 * \file        src/gravity/forcetree.c
 * \date        05/2018
 * \brief       Gravitational tree build.
 * \details     This file contains the construction of the tree used for
 *              calculating the gravitational force. The type tree implemented
 *              is a geometrical oct-tree, starting from a cube encompassing
 *              all particles. This cube is automatically found in the domain
 *              decomposition, which also splits up the global "top-level"
 *              tree along node boundaries, moving the particles of different
 *              parts of the tree to separate processors. In this version of
 *              the code, the tree construction may be repeated every timestep
 *              without a renewed domain decomposition. If particles are on
 *              the "wrong" processor because a new domain decomposition has
 *              not been carried out, they are sent as temporary points to the
 *              right insertion processor according to the layout of the
 *              top-level nodes. In addition, the mapping of the top-level
 *              nodes to processors may be readjusted in order to improve
 *              work-load balance for the current time step.
 *              contains functions:
 *                int construct_forcetree(int mode, int
 *                  optimized_domain_mapping, int insert_only_primary,
 *                  int timebin)
 *                int force_treebuild(int npart, int optimized_domain_mapping,
 *                  int insert_only_primary, int timebin)
 *                int force_treebuild_construct(int npart, int
 *                  optimized_domain_mapping, int insert_only_primary,
 *                  int timebin)
 *                int force_treebuild_insert_single_point(int i, unsigned
 *                  long long *intpos, int th, unsigned char levels)
 *                void force_assign_cost_values(void)
 *                int force_create_empty_nodes(int no, int topnode, int bits,
 *                  int x, int y, int z)
 *                void force_insert_pseudo_particles(void)
 *                void force_update_node_recursive(int no, int sib, int father,
 *                  int *last)
 *                void force_exchange_topleafdata(void)
 *                void force_treeupdate_toplevel(int no, int topnode, int bits,
 *                  int x, int y, int z)
 *                void force_treeallocate(int maxpart, int maxindex)
 *                void force_treefree(void)
 *                void dump_particles(void)
 *                int force_add_empty_nodes(void)
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
#include <time.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

static int *th_list;
static unsigned char *level_list;
int NTreeInsert;

#ifdef FOF
#ifndef FOF_SECONDARY_LINK_TARGET_TYPES
#define FOF_SECONDARY_LINK_TARGET_TYPES FOF_PRIMARY_LINK_TYPES
#endif /* #ifndef FOF_SECONDARY_LINK_TARGET_TYPES */
#endif /* #ifdef FOF */

#ifdef HIERARCHICAL_GRAVITY
#define INDEX(idx) (TimeBinsGravity.ActiveParticleList[idx])
#else /* #ifdef HIERARCHICAL_GRAVITY */
#define INDEX(idx) (idx)
#endif /* #ifdef HIERARCHICAL_GRAVITY #else */

/*! \brief Triggers forcetree construction until successful.
 *
 *  Allocates memory and constructs forcetree until successful;
 *  currently, there are two valid modes: forcetree only for gas or for all
 *  particles.
 *
 *  \param[in] mode Mode: all particles or just gas cells.
 *  \param[in] optimized_domain_mapping Handed over to force_treebuild.
 *  \param[in] insert_only_primary Handed over to force_treebuild.
 *  \param[in] timebin Handed over to force_treebuild.
 *
 *  \return Number of nodes in tree.
 */
int construct_forcetree(int mode, int optimized_domain_mapping, int insert_only_primary, int timebin)
{
  int npart, Tree_NumNodes = 0;

  do
    {
      /* Note: force_treebuild will call force_treefree if it is about to return a negative value!
       * Therefore, this has to be allocated within the loop! The only exception is when
       * insert_only_primary == 2, in which case the code assumes that the forcetree is already
       * allocated (this happens only in fof.c). In this case, force_treeallocate is not called
       * during the first loop.
       */
      if(insert_only_primary != 2 || Tree_NumNodes < 0)
        force_treeallocate(NumPart, All.MaxPart); /* reallocate force tree structure */

      /* prepare variables for force_treebuild call */
      switch(mode)
        {
          case 0: /* all particles */
            {
              npart = NumPart;
              break;
            }
          case 1: /* only gas particles */
            {
              npart = NumGas;
              break;
            }
          default:
            {
              mpi_terminate("FORCETREE: construct_forcetree: invalid mode!\n");
            }
        }

      Tree_NumNodes = force_treebuild(npart, optimized_domain_mapping, insert_only_primary, timebin);
    }
  while(Tree_NumNodes < 0);

  return Tree_NumNodes;
}

/*! \brief Constructs the gravitational oct-tree and handles errors.
 *
 *  \param[in] npart Number of particles on local task.
 *  \param[in] optimized_domain_mapping Specifies if mapping of the top-level
 *             nodes to processors may be readjusted.
 *  \param[in] insert_only_primary If this is set, only particles of the types
 *             set in FOF_PRIMARY_LINK_TYPES are inserted.
 *  \param[in] timebin Current timebin; needed for HIERARCHICAL_GRAVITY.
 *
 *  \return number of local+top nodes of the constructed tree.
 */
int force_treebuild(int npart, int optimized_domain_mapping, int insert_only_primary, int timebin)
{
  int i, flag;

#ifdef HIERARCHICAL_GRAVITY
  NTreeInsert              = TimeBinsGravity.NActiveParticles;
  optimized_domain_mapping = 0;
#else  /* #ifdef HIERARCHICAL_GRAVITY */
  NTreeInsert = npart;
#endif /* #ifdef HIERARCHICAL_GRAVITY #else */

  TIMER_START(CPU_TREEBUILD);

  long long loc_insert = NTreeInsert, tot_insert;
  MPI_Reduce(&loc_insert, &tot_insert, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  mpi_printf("FORCETREE: Tree construction.  (inserting %lld points)\n", tot_insert);

  TIMER_STOPSTART(CPU_TREEBUILD, CPU_TREEBUILD_INSERT);

  int flag_single = force_treebuild_construct(npart, optimized_domain_mapping, insert_only_primary, timebin);

  TIMER_STOPSTART(CPU_TREEBUILD_INSERT, CPU_TREEBUILD);

  MPI_Allreduce(&flag_single, &flag, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if(flag < 0)
    {
      /* tree construction was not successful and needs to be repeated */
      if(flag_single != -2)
        {
          myfree(Tree_Points);
        }

      force_treefree();

      if(flag == -3)
        {
          /* we need to do an extra domain decomposition to recover from an out-of-box condition for a particle,
             which can happen if GRAVITY_NOT_PERIODIC is used */
          ngb_treefree();
          domain_free();

          domain_Decomposition();

          ngb_treeallocate();
          ngb_treebuild(NumGas);
        }
      else
        {
          All.TreeAllocFactor *= 1.15;
          mpi_printf("FORCETREE: Increasing TreeAllocFactor, new value=%g\n", All.TreeAllocFactor);

          if(All.TreeAllocFactor > MAX_TREE_ALLOC_FACTOR)
            {
              char buf[500];
              sprintf(buf,
                      "task %d: looks like a serious problem in tree construction, stopping with particle dump.  Tree_NumNodes=%d "
                      "Tree_MaxNodes=%d  Tree_NumPartImported=%d NumPart=%d\n",
                      ThisTask, Tree_NumNodes, Tree_MaxNodes, Tree_NumPartImported, NumPart);
              dump_particles();
              terminate(buf);
            }
        }

      TIMER_STOP(CPU_TREEBUILD); /* stop timer before returning */
      return -1; /* stop right here with error code to invoke a new call of this function, possibly with changed values for npart */
    }            /* if(flag < 0) */

  Nextnode = (int *)mymalloc_movable(&Nextnode, "Nextnode", (Tree_MaxPart + NTopleaves + Tree_NumPartImported) * sizeof(int));
  Father   = (int *)mymalloc_movable(&Father, "Father", (Tree_MaxPart + Tree_NumPartImported) * sizeof(int));

  for(i = 0; i < Tree_MaxPart + Tree_NumPartImported; i++)
    Father[i] = -1;

  TIMER_STOPSTART(CPU_TREEBUILD, CPU_TREEBUILD_BRANCHES);

  /* insert the pseudo particles that represent the mass distribution of other domains */
  force_insert_pseudo_particles();

  /* now compute the multipole moments recursively */
  int last = -1;

  force_update_node_recursive(Tree_MaxPart, -1, -1, &last);

  if(last >= Tree_MaxPart)
    {
      if(last >= Tree_MaxPart + Tree_MaxNodes) /* a pseudo-particle or imported particle */
        Nextnode[last - Tree_MaxNodes] = -1;
      else
        Nodes[last].u.d.nextnode = -1;
    }
  else
    Nextnode[last] = -1;

  TIMER_STOPSTART(CPU_TREEBUILD_BRANCHES, CPU_TREEBUILD_TOPLEVEL);

  force_exchange_topleafdata();

  Tree_NextFreeNode = Tree_MaxPart + 1;
  force_treeupdate_toplevel(Tree_MaxPart, 0, 1, 0, 0, 0);

  TIMER_STOPSTART(CPU_TREEBUILD_TOPLEVEL, CPU_LOGS);

#ifdef HIERARCHICAL_GRAVITY
  if(timebin == All.HighestOccupiedGravTimeBin)
#endif /* #ifdef HIERARCHICAL_GRAVITY */
    {
      double locdata[2] = {Tree_NumPartImported, Tree_NumNodes}, sumdata[2];
      MPI_Reduce(locdata, sumdata, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      double tot_imported = sumdata[0];
      double tot_numnodes = sumdata[1];

      mpi_printf(
          "FORCETREE: Tree construction done.  <avg imported/local ratio>=%g <numnodes>=%g NTopnodes=%d NTopleaves=%d "
          "tree-build-scalability=%g\n",
          tot_imported / (All.TotNumPart + 1.0e-60), tot_numnodes / NTask, NTopnodes, NTopleaves,
          ((double)((tot_numnodes - NTask * ((double)NTopnodes)) + NTopnodes)) / (tot_numnodes + 1.0e-60));
    }
#ifdef HIERARCHICAL_GRAVITY
  else
    mpi_printf("FORCETREE: Tree construction done.\n");
#endif /* #ifdef HIERARCHICAL_GRAVITY */

  TIMER_STOP(CPU_LOGS);

  return Tree_NumNodes;
}

/*! \brief Constructs the gravitational oct-tree.
 *
 *  The index convention for accessing tree nodes is the following:
 *  node index
 *  [0...            Tree_MaxPart-1]  references single particles,
 *  the indices
 *  [Tree_MaxPart... Tree_MaxPart+Tree_MaxNodes-1]  references tree nodes
 *  [Tree_MaxPart+Tree_MaxNodes...  Tree_MaxPart+Tree_MaxNodes+NTopleaves-1]
 *     references "pseudo particles", i.e. mark branches on foreign CPUs
 *  [Tree_MaxPart+Tree_MaxNodes+NTopleaves...
 *     Tree_MaxPart+Tree_MaxNodes+NTopleaves+Tree_NumPartImported-1]
 *     references imported points.
 *
 *  the pointer `Nodes' is shifted such that Nodes[Tree_MaxPart] gives the
 *  first tree node (i.e. the root node).
 *
 *  \param[in] npart Number of particles on local task.
 *  \param[in] optimized_domain_mapping Specifies if mapping of the top-level
 *             nodes to processors may be readjusted.
 *  \param[in] insert_only_primary If this is set, only particles of the types
 *             set in FOF_PRIMARY_LINK_TYPES are inserted.
 *  \param[in] timebin (unused).
 *
 *  \return if successful returns the number of local+top nodes of the
 *             constructed tree;
 *          -1 if the number of allocated tree nodes is too small;
 *          -2 if the number of allocated tree nodes is even too small to fit
 *             the top nodes;
 *          -3 if a particle out of domain box condition was encountered.
 */
int force_treebuild_construct(int npart, int optimized_domain_mapping, int insert_only_primary, int timebin)
{
  int idx, i, j, no, flag = 0;
  int ngrp, recvTask, count_ListNoData, *no_place = NULL;
  unsigned long long *intposp;
  MyDouble *posp;

#ifdef DISABLE_OPTIMIZE_DOMAIN_MAPPING
  optimized_domain_mapping = 0;
#endif /* #ifdef DISABLE_OPTIMIZE_DOMAIN_MAPPING */

#if !defined(GRAVITY_NOT_PERIODIC)
  double boxsize[3];
  boxsize[0] = boxSize_X;
  boxsize[1] = boxSize_Y;
  boxsize[2] = boxSize_Z;
#endif /* #if !defined(GRAVITY_NOT_PERIODIC) */

  /* create an empty root node  */
  Tree_NextFreeNode   = Tree_MaxPart;              /* index of first free node */
  struct NODE *nfreep = &Nodes[Tree_NextFreeNode]; /* select first node        */

  for(j = 0; j < 8; j++)
    nfreep->u.suns[j] = -1;

  nfreep->len = DomainLen;
  for(j = 0; j < 3; j++)
    nfreep->center[j] = DomainCenter[j];

  Tree_NumNodes = 1;
  Tree_NextFreeNode++;

  /* create a set of empty nodes corresponding to the top-level domain
   * grid. We need to generate these nodes first to make sure that we have a
   * complete top-level tree which allows the easy insertion of the
   * pseudo-particles at the right place
   */
  if(force_create_empty_nodes(Tree_MaxPart, 0, 1, 0, 0, 0) < 0)
    return -2;

  Tree_FirstNonTopLevelNode = Tree_NextFreeNode;

  /* if a high-resolution region in a global tree is used, we need to generate
   * an additional set of empty nodes to make sure that we have a complete
   * top-level tree for the high-resolution inset
   */

  /* we first do a dummy allocation here that we'll resize later if needed, in which case the following arrays will have to be moved
   * once. */
  int guess_nimported = 1.2 * NumPart;

  Tree_Points =
      (struct treepoint_data *)mymalloc_movable(&Tree_Points, "Tree_Points", guess_nimported * sizeof(struct treepoint_data));

  th_list    = (int *)mymalloc_movable(&th_list, "th_list", NumPart * sizeof(int));
  level_list = (unsigned char *)mymalloc_movable(&level_list, "level_list", NumPart * sizeof(unsigned char));
  Tree_IntPos_list =
      (unsigned long long *)mymalloc_movable(&Tree_IntPos_list, "Tree_IntPos_list", 3 * NumPart * sizeof(unsigned long long));

  if(NumPart < NTreeInsert)
    {
      terminate("ERROR: NumPart %d, NTreeInsert %d! This should not happen!", NumPart, NTreeInsert);
    }

  /* first check whether particles are still in domain box */
  for(idx = 0; idx < NTreeInsert; idx++)
    {
      i = INDEX(idx);
      if(i < 0)
        continue;

      if(P[i].Ti_Current != All.Ti_Current)
        drift_particle(i, All.Ti_Current);

      posp = &Tree_Pos_list[i * 3];

      for(j = 0; j < 3; j++, posp++)
        {
#ifdef CELL_CENTER_GRAVITY
          if(P[i].Type == 0)
            *posp = SphP[i].Center[j];
          else
#endif /* #ifdef CELL_CENTER_GRAVITY */
            *posp = P[i].Pos[j];

#if !defined(GRAVITY_NOT_PERIODIC)
          if(*posp < 0)
            *posp += boxsize[j];
          if(*posp >= boxsize[j])
            *posp -= boxsize[j];
#endif /* #if !defined(GRAVITY_NOT_PERIODIC) */
          if(*posp < DomainCorner[j] || *posp >= DomainCorner[j] + DomainLen)
            {
              flag = 1;
              break;
            }
        }
    }

#if defined(GRAVITY_NOT_PERIODIC)
  int flag_sum;
  MPI_Allreduce(&flag, &flag_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(flag_sum)
    {
      mpi_printf(
          "FORCETREE: Particle out of domain box condition was triggered. Need to do an (unplanned) new domain decomposition.\n");
      myfree(Tree_IntPos_list);
      myfree(level_list);
      myfree(th_list);
      return -3;
    }
#else  /* #if defined(GRAVITY_NOT_PERIODIC) */
  if(flag)
    {
      char buf[1000];
      sprintf(buf, "i=%d ID=%lld type=%d moved out of box. Pos[j=%d]=%g DomainCorner[%d]=%g DomainLen=%g", i, (long long)P[i].ID,
              P[i].Type, j, P[i].Pos[j], j, DomainCorner[j], DomainLen);
      terminate(buf);
    }
#endif /* #if defined(GRAVITY_NOT_PERIODIC) #else */

#if defined(EVALPOTENTIAL) && defined(PMGRID) && !defined(GRAVITY_NOT_PERIODIC)
  double mass_highres = 0, mass_lowres = 0;
  for(int idx = 0; idx < NTreeInsert; idx++)
    {
      int i = INDEX(idx);
      if(i < 0)
        continue;

#ifdef PLACEHIGHRESREGION
      if(pmforce_is_particle_high_res(P[i].Type, &Tree_Pos_list[3 * i]))
        mass_highres += P[i].Mass;
      else
#endif /* #ifdef PLACEHIGHRESREGION */
        mass_lowres += P[i].Mass;
    }
  double mass_pmregions[2] = {mass_lowres, mass_highres};
  MPI_Allreduce(mass_pmregions, All.MassPMregions, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif /* #if defined(EVALPOTENTIAL) && defined(PMGRID) && !defined(GRAVITY_NOT_PERIODIC) */

  /* now we determine for each point the insertion top-level node, and the task on which this lies */
  if(optimized_domain_mapping)
    {
      TaskCost    = mymalloc("TaskCost", NTask * sizeof(double));
      TaskCount   = mymalloc("TaskCount", NTask * sizeof(int));
      DomainCost  = mymalloc("DomainCost", NTopleaves * sizeof(double));
      DomainCount = mymalloc("DomainCount", NTopleaves * sizeof(int));
      ListNoData  = mymalloc("ListNoData", NTopleaves * sizeof(struct no_list_data));
      no_place    = mymalloc("no_place", NTopleaves * sizeof(int));

      memset(no_place, -1, NTopleaves * sizeof(int));

      for(j = 0; j < NTopleaves; j++)
        DomainCost[j] = 0;
      for(j = 0; j < NTopleaves; j++)
        DomainCount[j] = 0;
      for(j = 0; j < NTask; j++)
        TaskCost[j] = 0;

      for(j = 0; j < NTask; j++)
        Send_count[j] = 0;

      count_ListNoData = 0;
    }

  for(idx = 0; idx < NTreeInsert; idx++)
    {
      i = INDEX(idx);
      if(i < 0)
        continue;

      posp = &Tree_Pos_list[i * 3];

      unsigned long long xxb  = force_double_to_int(((*posp++ - DomainCorner[0]) * DomainInverseLen) + 1.0);
      unsigned long long yyb  = force_double_to_int(((*posp++ - DomainCorner[1]) * DomainInverseLen) + 1.0);
      unsigned long long zzb  = force_double_to_int(((*posp++ - DomainCorner[2]) * DomainInverseLen) + 1.0);
      unsigned long long mask = ((unsigned long long)1) << (52 - 1);
      unsigned char shiftx    = (52 - 1);
      unsigned char shifty    = (52 - 2);
      unsigned char shiftz    = (52 - 3);
      unsigned char levels    = 0;

      intposp    = &Tree_IntPos_list[i * 3];
      *intposp++ = xxb;
      *intposp++ = yyb;
      *intposp++ = zzb;

      no = 0;
      while(TopNodes[no].Daughter >= 0) /* walk down top tree to find correct leaf */
        {
          unsigned char subnode = (((unsigned char)((xxb & mask) >> (shiftx--))) | ((unsigned char)((yyb & mask) >> (shifty--))) |
                                   ((unsigned char)((zzb & mask) >> (shiftz--))));

          mask >>= 1;
          levels++;

          no = TopNodes[no].Daughter + TopNodes[no].MortonToPeanoSubnode[subnode];
        }

      no = TopNodes[no].Leaf;

      th_list[i]    = no;
      level_list[i] = levels;

      if(optimized_domain_mapping)
        {
          /* find costs for all top leaves */

          int bin = All.HighestActiveTimeBin;
          double cost;

          if(domain_bintolevel[bin] >= 0)
            cost = MIN_FLOAT_NUMBER + P[i].GravCost[domain_bintolevel[bin]] * domain_grav_weight[bin];
          else
            {
              if(domain_refbin[bin] >= 0)
                cost = MIN_FLOAT_NUMBER + P[i].GravCost[domain_bintolevel[domain_refbin[bin]]] * domain_grav_weight[bin];
              else
                cost = 1.0;
            }

          int task = DomainTask[no];
          TaskCost[task] += cost;

          if(task == ThisTask)
            {
              DomainCost[no] += cost;
              DomainCount[no]++;
            }
          else
            {
              int p = no_place[no];
              if(p >= 0)
                {
                  ListNoData[p].domainCost += cost;
                  ListNoData[p].domainCount++;
                }
              else
                {
                  Send_count[task]++;
                  p                         = count_ListNoData++;
                  no_place[no]              = p;
                  ListNoData[p].task        = task;
                  ListNoData[p].no          = no;
                  ListNoData[p].domainCost  = cost;
                  ListNoData[p].domainCount = 1;
                }
            }
        }
    }

  if(optimized_domain_mapping)
    {
      /* if necessary, re-adjust the mapping of the top-level nodes to the processors */

      if(All.Ti_Current > 0)
        {
          double current_balance, impact;
          current_balance = force_get_current_balance(&impact);

          mpi_printf("FORCETREE: current balance=  %g | %g\n", current_balance, impact);

          if(All.HighestActiveTimeBin <
             All.SmallestTimeBinWithDomainDecomposition) /* only do this for steps which did not do a domain decomposition */
            {
              if(impact > MAX_IMPACT_BEFORE_OPTIMIZATION)
                {
                  force_get_global_cost_for_leavenodes(count_ListNoData);
                  force_optimize_domain_mapping();
                }
              else
                {
                  mpi_printf(
                      "FORCETREE: we're not trying to optimize further because overall imbalance impact is only %g (threshold is "
                      "%g)\n",
                      impact, MAX_IMPACT_BEFORE_OPTIMIZATION);
                  memcpy(DomainNewTask, DomainTask, NTopleaves * sizeof(int));
                }
            }
          else
            {
              mpi_printf("FORCETREE: we're not trying to optimize futher because we just did a domain decomposition\n");
              memcpy(DomainNewTask, DomainTask, NTopleaves * sizeof(int));
            }
        }
      else
        memcpy(DomainNewTask, DomainTask, NTopleaves * sizeof(int));
    }
  else
    memcpy(DomainNewTask, DomainTask, NTopleaves * sizeof(int));

  if(optimized_domain_mapping)
    {
      myfree(no_place);
      myfree(ListNoData);
      myfree(DomainCount);
      myfree(DomainCost);
      myfree(TaskCount);
      myfree(TaskCost);
    }

  for(j = 0; j < NTask; j++)
    {
      Force_Send_count[j] = 0;
    }

  for(idx = 0; idx < NTreeInsert; idx++) /* make list of insertion top leaf and task for all particles */
    {
      i = INDEX(idx);
      if(i < 0)
        continue;

      no         = th_list[i];
      th_list[i] = DomainNodeIndex[no];

      int task = DomainNewTask[no];

      Tree_Task_list[i] = task;

      if(task != ThisTask)
        {
          Force_Send_count[task]++;
        }
    }

  MPI_Alltoall(Force_Send_count, 1, MPI_INT, Force_Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, Tree_NumPartImported = 0, Tree_NumPartExported = 0, Force_Recv_offset[0] = 0, Force_Send_offset[0] = 0; j < NTask; j++)
    {
      Tree_NumPartImported += Force_Recv_count[j];
      Tree_NumPartExported += Force_Send_count[j];
      if(j > 0)
        {
          Force_Send_offset[j] = Force_Send_offset[j - 1] + Force_Send_count[j - 1];
          Force_Recv_offset[j] = Force_Recv_offset[j - 1] + Force_Recv_count[j - 1];
        }
    }

  if(Tree_NumPartImported > guess_nimported)
    {
      printf("ThisTask=%d: Tree_NumPartImported=%d  NumPart=%d\n", ThisTask, Tree_NumPartImported, NumPart);
      Tree_Points = (struct treepoint_data *)myrealloc_movable(Tree_Points, Tree_NumPartImported * sizeof(struct treepoint_data));
    }

  if(Tree_NumPartImported > 0.25 * NumPart)
    {
      Tree_MaxNodes = (int)(All.TreeAllocFactor * (NumPart + Tree_NumPartImported)) + NTopnodes;

      Nodes += Tree_MaxPart;
      Nodes = (struct NODE *)myrealloc_movable(Nodes, (Tree_MaxNodes + 1) * sizeof(struct NODE));
      Nodes -= Tree_MaxPart;

#ifdef MULTIPLE_NODE_SOFTENING
      ExtNodes += Tree_MaxPart;
      ExtNodes = (struct ExtNODE *)myrealloc_movable(ExtNodes, (Tree_MaxNodes + 1) * sizeof(struct ExtNODE));
      ExtNodes -= Tree_MaxPart;
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
    }

  struct treepoint_data *export_Tree_Points =
      (struct treepoint_data *)mymalloc("export_Tree_Points", Tree_NumPartExported * sizeof(struct treepoint_data));

  for(j = 0; j < NTask; j++)
    {
      Force_Send_count[j] = 0;
    }

  for(idx = 0; idx < NTreeInsert; idx++) /* prepare particle data to be copied to other tasks */
    {
      i = INDEX(idx);
      if(i < 0)
        continue;

      int task = Tree_Task_list[i];

      if(task != ThisTask)
        {
          int n = Force_Send_offset[task] + Force_Send_count[task]++;

          /* this point has to go to another task */
          export_Tree_Points[n].Pos[0]        = Tree_Pos_list[3 * i + 0];
          export_Tree_Points[n].Pos[1]        = Tree_Pos_list[3 * i + 1];
          export_Tree_Points[n].Pos[2]        = Tree_Pos_list[3 * i + 2];
          export_Tree_Points[n].IntPos[0]     = Tree_IntPos_list[3 * i + 0];
          export_Tree_Points[n].IntPos[1]     = Tree_IntPos_list[3 * i + 1];
          export_Tree_Points[n].IntPos[2]     = Tree_IntPos_list[3 * i + 2];
          export_Tree_Points[n].Mass          = P[i].Mass;
          export_Tree_Points[n].OldAcc        = P[i].OldAcc;
          export_Tree_Points[n].SofteningType = P[i].SofteningType;
          export_Tree_Points[n].index         = i;
          export_Tree_Points[n].Type          = P[i].Type;
          export_Tree_Points[n].th            = th_list[i];
          export_Tree_Points[n].level         = level_list[i];
#ifndef HIERARCHICAL_GRAVITY
          if(TimeBinSynchronized[P[i].TimeBinGrav])
            export_Tree_Points[n].ActiveFlag = 1;
          else
            export_Tree_Points[n].ActiveFlag = 0;
#endif /* #ifndef HIERARCHICAL_GRAVITY */
        }
    }

  /* exchange  data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Force_Send_count[recvTask] > 0 || Force_Recv_count[recvTask] > 0)
          MPI_Sendrecv(&export_Tree_Points[Force_Send_offset[recvTask]], Force_Send_count[recvTask] * sizeof(struct treepoint_data),
                       MPI_BYTE, recvTask, TAG_DENS_A, &Tree_Points[Force_Recv_offset[recvTask]],
                       Force_Recv_count[recvTask] * sizeof(struct treepoint_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
    }

  myfree(export_Tree_Points);

  Tree_ImportedNodeOffset = Tree_MaxPart + Tree_MaxNodes + NTopleaves;

  int full_flag = 0;

  /* now we insert all particles */
  for(idx = 0; idx < NTreeInsert; idx++)
    {
      i = INDEX(idx);
      if(i < 0)
        continue;

#ifdef NO_GAS_SELFGRAVITY
      if(P[i].Type == 0)
        continue;
#endif /* #ifdef NO_GAS_SELFGRAVITY */
#ifdef NO_SELFGRAVITY_TYPE
      if(P[i].Type == NO_SELFGRAVITY_TYPE)
        continue;
#endif /* #ifdef NO_SELFGRAVITY_TYPE */
#if defined(FOF) || defined(SUBFIND)
      if(insert_only_primary == 1)
        {
          if(!((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES)))
            continue;
        }
      else if(insert_only_primary == 2)
        {
          if(!((1 << P[i].Type) & (FOF_SECONDARY_LINK_TARGET_TYPES)))
            continue;
        }
#endif /* #if defined(FOF) || defined(SUBFIND) */
      if(Tree_Task_list[i] == ThisTask)
        {
          if(force_treebuild_insert_single_point(i, &Tree_IntPos_list[3 * i], th_list[i], level_list[i]) < 0)
            {
              full_flag = 1;
              break;
            }
        }
    }

  if(full_flag == 0) /* only continue if previous step was successful */
    {
      for(i = 0; i < Tree_NumPartImported; i++)
        {
#ifdef NO_GAS_SELFGRAVITY
          if(Tree_Points[i].Type == 0)
            continue;
#endif /* #ifdef NO_GAS_SELFGRAVITY */
#ifdef NO_SELFGRAVITY_TYPE
          if(Tree_Points[i].Type == NO_SELFGRAVITY_TYPE)
            continue;
#endif /* #ifdef NO_SELFGRAVITY_TYPE */
#if defined(FOF) || defined(SUBFIND)
          if(insert_only_primary == 1)
            {
              if(!((1 << Tree_Points[i].Type) & (FOF_PRIMARY_LINK_TYPES)))
                continue;
            }
          else if(insert_only_primary == 2)
            {
              if(!((1 << Tree_Points[i].Type) & (FOF_SECONDARY_LINK_TARGET_TYPES)))
                continue;
            }
#endif /* #if defined(FOF) || defined(SUBFIND) */
          if(force_treebuild_insert_single_point(i + Tree_ImportedNodeOffset, Tree_Points[i].IntPos, Tree_Points[i].th,
                                                 Tree_Points[i].level) < 0)
            {
              full_flag = 1;
              break;
            }
        }
    }

  myfree_movable(Tree_IntPos_list);
  myfree_movable(level_list);
  myfree_movable(th_list);

  if(full_flag)
    return -1;

#ifdef ADDBACKGROUNDGRID
  if(force_add_empty_nodes())
    return -1;
#endif /* #ifdef ADDBACKGROUNDGRID */

  return Tree_NumNodes;
}

/*! \brief Inserts a single particle into the gravitational tree.
 *
 *  \param[in] i Index of particle.
 *  \param[in] intpos Integer representation of particle position.
 *  \param[in] th Target node.
 *  \param[in] levels Level of target node.
 *
 *  \return 0 if successful;
 *          -1 if too few nodes have been allocated in the Nodes array
 */
int force_treebuild_insert_single_point(int i, unsigned long long *intpos, int th, unsigned char levels)
{
  int j, parent = -1;
  unsigned char subnode       = 0;
  unsigned long long xxb      = intpos[0];
  unsigned long long yyb      = intpos[1];
  unsigned long long zzb      = intpos[2];
  unsigned long long mask     = ((unsigned long long)1) << ((52 - 1) - levels);
  unsigned char shiftx        = (52 - 1) - levels;
  unsigned char shifty        = (52 - 2) - levels;
  unsigned char shiftz        = (52 - 3) - levels;
  signed long long centermask = (0xFFF0000000000000llu);
  unsigned long long *intppos;
  centermask >>= levels;

  while(1)
    {
      if(th >= Tree_MaxPart && th < Tree_ImportedNodeOffset) /* we are dealing with an internal node */
        {
          subnode = (((unsigned char)((xxb & mask) >> (shiftx--))) | ((unsigned char)((yyb & mask) >> (shifty--))) |
                     ((unsigned char)((zzb & mask) >> (shiftz--))));

          centermask >>= 1;
          mask >>= 1;
          levels++;

          if(levels > MAX_TREE_LEVEL)
            {
              /* seems like we're dealing with particles at identical (or extremely close)
               * locations. Shift subnode index to allow tree construction. Note: Multipole moments
               * of tree are still correct, but one should MAX_TREE_LEVEL large enough to have
               *      DomainLen/2^MAX_TREE_LEEL  < gravitational softening length
               */
              for(j = 0; j < 8; j++)
                {
                  if(Nodes[th].u.suns[subnode] < 0)
                    break;

                  subnode++;
                  if(subnode >= 8)
                    subnode = 7;
                }
            }

          int nn = Nodes[th].u.suns[subnode];

          if(nn >= 0) /* ok, something is in the daughter slot already, need to continue */
            {
              parent = th;
              th     = nn;
            }
          else
            {
              /* here we have found an empty slot where we can attach
               * the new particle as a leaf.
               */
              Nodes[th].u.suns[subnode] = i;
              break; /* done for this particle */
            }
        }
      else
        {
          /* We try to insert into a leaf with a single particle.  Need
           * to generate a new internal node at this point.
           */
          Nodes[parent].u.suns[subnode] = Tree_NextFreeNode;
          struct NODE *nfreep           = &Nodes[Tree_NextFreeNode];

          double len = ((double)(mask << 1)) * DomainBigFac;
          double cx  = ((double)((xxb & centermask) | mask)) * DomainBigFac + DomainCorner[0];
          double cy  = ((double)((yyb & centermask) | mask)) * DomainBigFac + DomainCorner[1];
          double cz  = ((double)((zzb & centermask) | mask)) * DomainBigFac + DomainCorner[2];

          nfreep->len       = len;
          nfreep->center[0] = cx;
          nfreep->center[1] = cy;
          nfreep->center[2] = cz;

          for(j = 0; j < 8; j++)
            nfreep->u.suns[j] = -1;

          if(th >= Tree_ImportedNodeOffset)
            intppos = Tree_Points[th - Tree_ImportedNodeOffset].IntPos;
          else
            intppos = &Tree_IntPos_list[3 * th];

          subnode = (((unsigned char)((intppos[0] & mask) >> shiftx)) | ((unsigned char)((intppos[1] & mask) >> shifty)) |
                     ((unsigned char)((intppos[2] & mask) >> shiftz)));

          nfreep->u.suns[subnode] = th;

          th = Tree_NextFreeNode; /* resume trying to insert the new particle the newly created internal node */
          Tree_NumNodes++;
          Tree_NextFreeNode++;

          if(Tree_NumNodes >= Tree_MaxNodes)
            {
              return -1;
            }
        }
    }

  return 0;
}

/*! \brief Distributes the gravity costs of each node among the particles it
 *         contains.
 *
 *  \return void
 */
void force_assign_cost_values(void)
{
  int idx, i, ngrp, recvTask;

  if(TakeLevel >= 0)
    {
      int thread;

      /* consolidate the cost measurements done by the different threads */
      for(thread = 1; thread < NUM_THREADS; thread++)
        for(i = 0; i < NumPart; i++)
          Thread[0].P_CostCount[i] += Thread[thread].P_CostCount[i];

      for(thread = 1; thread < NUM_THREADS; thread++)
        for(i = 0; i < Tree_NumNodes; i++)
          Thread[0].Node_CostCount[i + Tree_MaxPart] += Thread[thread].Node_CostCount[i + Tree_MaxPart];

      for(thread = 1; thread < NUM_THREADS; thread++)
        for(i = 0; i < Tree_NumPartImported; i++)
          Thread[0].TreePoints_CostCount[i] += Thread[thread].TreePoints_CostCount[i];

#ifdef VERBOSE
      /* calculate some check sums to validate the total cost assignment */
      double sumbefore = 0, sumbeforetot;
      for(i = 0; i < NumPart; i++)
        sumbefore += P[i].GravCost[TakeLevel];
      MPI_Allreduce(&sumbefore, &sumbeforetot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      double nodecost = 0, nodecosttot;
      for(i = 0; i < Tree_NumNodes; i++)
        nodecost += Thread[0].Node_CostCount[i + Tree_MaxPart];
      MPI_Allreduce(&nodecost, &nodecosttot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      double importedcost = 0, importedcosttot;
      for(i = 0; i < Tree_NumPartImported; i++)
        importedcost += Thread[0].TreePoints_CostCount[i];
      MPI_Allreduce(&importedcost, &importedcosttot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      double partcost = 0, partcosttot;
      for(idx = 0; idx < NTreeInsert; idx++)
        {
          i = INDEX(idx);
          if(i < 0)
            continue;

          {
            int no = Father[i];

            if(no >= 0)
              partcost += Thread[0].P_CostCount[i];
          }
        }
      MPI_Allreduce(&partcost, &partcosttot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif /* #ifdef VERBOSE */

      double *loc_cost  = mymalloc("loc_cost", NTopnodes * sizeof(double));
      double *glob_cost = mymalloc("glob_cost", NTopnodes * sizeof(double));

      for(i = 0; i < NTopnodes; i++)
        loc_cost[i] = Thread[0].Node_CostCount[i + Tree_MaxPart];

      MPI_Allreduce(loc_cost, glob_cost, NTopnodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for(i = 0; i < NTopnodes; i++)
        Thread[0].Node_CostCount[i + Tree_MaxPart] = glob_cost[i];

      myfree(glob_cost);
      myfree(loc_cost);

      for(i = 0; i < NumPart; i++)
        P[i].GravCost[TakeLevel] = 0;

      /* distribute costs of parent nodes to particles */
      for(idx = 0; idx < NTreeInsert; idx++)
        {
          i = INDEX(idx);
          if(i < 0)
            continue;

          {
            double sum = Thread[0].P_CostCount[i];

            int no = Father[i];

            while(no >= 0)
              {
                if(Nodes[no].u.d.mass > 0)
                  sum += Thread[0].Node_CostCount[no] * (P[i].Mass / Nodes[no].u.d.mass);

                no = Nodes[no].u.d.father;
              }

            P[i].GravCost[TakeLevel] = sum;
          }
        }

      /* Now, if we moved points to other CPUs, we need to collect these cost values */
      struct gravcost_data
      {
        float GravCost;
        int index;
      } * gdata_export, *gdata_import;

      gdata_export = mymalloc("grav_data_export", Tree_NumPartExported * sizeof(struct gravcost_data));
      gdata_import = mymalloc("grav_data_import", Tree_NumPartImported * sizeof(struct gravcost_data));

      for(i = 0; i < Tree_NumPartImported; i++)
        {
          double sum = Thread[0].TreePoints_CostCount[i];

          int no = Father[i + Tree_MaxPart];

          while(no >= 0)
            {
              if(Nodes[no].u.d.mass > 0)
                sum += Thread[0].Node_CostCount[no] * Tree_Points[i].Mass / Nodes[no].u.d.mass;

              no = Nodes[no].u.d.father;
            }

          gdata_import[i].GravCost = sum;
          gdata_import[i].index    = Tree_Points[i].index;
        }

      /* exchange  data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;

          if(recvTask < NTask)
            {
              if(Force_Send_count[recvTask] > 0 || Force_Recv_count[recvTask] > 0)
                {
                  MPI_Sendrecv(&gdata_import[Force_Recv_offset[recvTask]], Force_Recv_count[recvTask] * sizeof(struct gravcost_data),
                               MPI_BYTE, recvTask, TAG_DENS_A, &gdata_export[Force_Send_offset[recvTask]],
                               Force_Send_count[recvTask] * sizeof(struct gravcost_data), MPI_BYTE, recvTask, TAG_DENS_A,
                               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

      for(i = 0; i < Tree_NumPartExported; i++)
        P[gdata_export[i].index].GravCost[TakeLevel] = gdata_export[i].GravCost;

      myfree(gdata_import);
      myfree(gdata_export);

#ifdef VERBOSE
      double sum = 0, sumtot;
      for(i = 0; i < NumPart; i++)
        sum += P[i].GravCost[TakeLevel];
      MPI_Allreduce(&sum, &sumtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      mpi_printf(
          "FORCETREE: Cost assignment for TakeLevel=%d, highest active-TimeBin=%d   yields cost=%g|%g (before %g)  nodecosttot=%g  "
          "partcosttot=%g importedcosttot=%g\n",
          TakeLevel, All.HighestActiveTimeBin, sumtot, nodecosttot + partcosttot + importedcosttot, sumbeforetot, nodecosttot,
          partcosttot, importedcosttot);
#else  /* #ifdef VERBOSE */
      mpi_printf("FORCETREE: Cost assignment for TakeLevel=%d, highest active-TimeBin=%d\n", TakeLevel, All.HighestActiveTimeBin);
#endif /* #ifdef VERBOSE #else */
    }
}

/*! \brief Recursively creates a set of empty tree nodes which
 *         corresponds to the top-level tree for the domain grid.
 *
 *  This is done to ensure that this top-level tree is always "complete" so
 *  that we can easily associate the pseudo-particles of other CPUs with
 *  tree-nodes at a given level in the tree, even when the particle population
 *  is so sparse that some of these nodes are actually empty.
 *
 * \param[in] no Parent node for which daughter nodes shall be created.
 * \param[in] topnode Index of the parent node in the 'TopNodes' array.
 * \param[in] bits 2^bits is the number of nodes per dimension at the level of
 *            the daughter nodes.
 * \param[in] x Position of the parent node in the x direction, falls in the
 *            range [0,2^(bits-1) - 1].
 * \param[in] y Position of the parent node in the y direction, falls in the
 *            range [0,2^(bits-1) - 1].
 * \param[in] z Position of the parent node in the z direction, falls in the
 *            range [0,2^(bits-1) - 1].
 *
 * \return 0 if successful;
 *         -1 if number of allocated tree nodes is too small to fit the newly
 *         created nodes.
 */
int force_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z)
{
  if(TopNodes[topnode].Daughter >= 0)
    {
      for(int i = 0; i < 2; i++) /* loop over daughter nodes */
        for(int j = 0; j < 2; j++)
          for(int k = 0; k < 2; k++)
            {
              if(Tree_NumNodes >= Tree_MaxNodes)
                {
                  if(All.TreeAllocFactor > MAX_TREE_ALLOC_FACTOR)
                    {
                      char buf[500];
                      sprintf(buf, "task %d: looks like a serious problem (NTopnodes=%d), stopping with particle dump.\n", ThisTask,
                              NTopnodes);
                      dump_particles();
                      terminate(buf);
                    }
                  return -1;
                }

              int sub = 7 & peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);

              int count = i + 2 * j + 4 * k;

              Nodes[no].u.suns[count] = Tree_NextFreeNode;

              double lenhalf                     = 0.25 * Nodes[no].len;
              Nodes[Tree_NextFreeNode].len       = 0.5 * Nodes[no].len;
              Nodes[Tree_NextFreeNode].center[0] = Nodes[no].center[0] + (2 * i - 1) * lenhalf;
              Nodes[Tree_NextFreeNode].center[1] = Nodes[no].center[1] + (2 * j - 1) * lenhalf;
              Nodes[Tree_NextFreeNode].center[2] = Nodes[no].center[2] + (2 * k - 1) * lenhalf;

              for(int n = 0; n < 8; n++)
                Nodes[Tree_NextFreeNode].u.suns[n] = -1;

              if(TopNodes[TopNodes[topnode].Daughter + sub].Daughter == -1)
                DomainNodeIndex[TopNodes[TopNodes[topnode].Daughter + sub].Leaf] = Tree_NextFreeNode;

              Tree_NextFreeNode++;
              Tree_NumNodes++;

              if(force_create_empty_nodes(Tree_NextFreeNode - 1, TopNodes[topnode].Daughter + sub, bits + 1, 2 * x + i, 2 * y + j,
                                          2 * z + k) < 0)
                return -1; /* create granddaughter nodes for current daughter node */
            }
    }

  return 0;
}

/*! \brief Inserts pseudo particles.
 *
 *  This function inserts pseudo-particles which will represent the mass
 *  distribution of the other CPUs. Initially, the mass of the
 *  pseudo-particles is set to zero, and their coordinate is set to the
 *  center of the domain-cell they correspond to. These quantities will be
 *  updated later on.
 *
 *  \return void
 */
void force_insert_pseudo_particles(void)
{
  for(int i = 0; i < NTopleaves; i++)
    {
      int index = DomainNodeIndex[i];

      if(DomainNewTask[i] != ThisTask)
        Nodes[index].u.suns[0] = Tree_MaxPart + Tree_MaxNodes + i;
    }
}

/*! \brief Determines multipole moments.
 *
 *  This routine determines the multipole moments for a given internal node
 *  and all its subnodes using a recursive computation.  The result is
 *  stored in the Nodes[] structure in the sequence of this tree-walk.
 *
 *  \param[in] no Node for which the moments shall be found.
 *  \param[in] sib Sibling of node no.
 *  \param[in] father Father node of node no.
 *  \param[in, out] last Last node for which this function was called, or -1
 *                  when called for root node.
 *
 *  \return void
 */
void force_update_node_recursive(int no, int sib, int father, int *last)
{
  int j, jj, p, pp, nextsib, suns[8];
  double s[3], mass;
  unsigned char maxsofttype;
#ifdef MULTIPLE_NODE_SOFTENING
  double mass_per_type[NSOFTTYPES];
#ifdef ADAPTIVE_HYDRO_SOFTENING
  unsigned char maxhydrosofttype;
  unsigned char minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

  if(no >= Tree_MaxPart && no < Tree_MaxPart + Tree_MaxNodes) /* internal node */
    {
      for(j = 0; j < 8; j++)
        suns[j] = Nodes[no].u.suns[j]; /* this "backup" is necessary because the nextnode entry will
                                          overwrite one element (union!) */
      if(*last >= 0)
        {
          if(*last >= Tree_MaxPart)
            {
              if(*last >= Tree_MaxPart + Tree_MaxNodes)
                Nextnode[*last - Tree_MaxNodes] = no; /* a pseudo-particle or imported point */
              else
                Nodes[*last].u.d.nextnode = no;
            }
          else
            Nextnode[*last] = no;
        }

      *last = no;

      mass        = 0;
      s[0]        = 0;
      s[1]        = 0;
      s[2]        = 0;
      maxsofttype = NSOFTTYPES + NSOFTTYPES_HYDRO;

#ifdef MULTIPLE_NODE_SOFTENING
      for(j = 0; j < NSOFTTYPES; j++)
        mass_per_type[j] = 0;

#ifdef ADAPTIVE_HYDRO_SOFTENING
      maxhydrosofttype = NSOFTTYPES;
      minhydrosofttype = NSOFTTYPES + NSOFTTYPES_HYDRO - 1;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

      for(j = 0; j < 8; j++)
        {
          if((p = suns[j]) >= 0)
            {
              /* check if we have a sibling on the same level */
              for(jj = j + 1; jj < 8; jj++)
                if((pp = suns[jj]) >= 0)
                  break;

              if(jj < 8) /* yes, we do */
                nextsib = pp;
              else
                nextsib = sib;

              force_update_node_recursive(p, nextsib, no, last);

              if(p < Tree_MaxPart) /* a particle */
                {
                  MyDouble *pos = &Tree_Pos_list[3 * p];

                  mass += P[p].Mass;
                  s[0] += P[p].Mass * pos[0];
                  s[1] += P[p].Mass * pos[1];
                  s[2] += P[p].Mass * pos[2];

                  if(All.ForceSoftening[maxsofttype] < All.ForceSoftening[P[p].SofteningType])
                    maxsofttype = P[p].SofteningType;

#ifdef MULTIPLE_NODE_SOFTENING
#ifdef ADAPTIVE_HYDRO_SOFTENING
                  mass_per_type[P[p].Type == 0 ? 0 : P[p].SofteningType] += P[p].Mass;

                  if(P[p].Type == 0)
                    {
                      if(maxhydrosofttype < P[p].SofteningType)
                        maxhydrosofttype = P[p].SofteningType;
                      if(minhydrosofttype > P[p].SofteningType)
                        minhydrosofttype = P[p].SofteningType;
                    }
#else  /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
                  mass_per_type[P[p].SofteningType] += P[p].Mass;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING #else */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
                }
              else if(p < Tree_MaxPart + Tree_MaxNodes) /* an internal node  */
                {
                  mass += Nodes[p].u.d.mass;
                  s[0] += Nodes[p].u.d.mass * Nodes[p].u.d.s[0];
                  s[1] += Nodes[p].u.d.mass * Nodes[p].u.d.s[1];
                  s[2] += Nodes[p].u.d.mass * Nodes[p].u.d.s[2];

                  if(All.ForceSoftening[maxsofttype] < All.ForceSoftening[Nodes[p].u.d.maxsofttype])
                    maxsofttype = Nodes[p].u.d.maxsofttype;

#ifdef MULTIPLE_NODE_SOFTENING
                  int k;
                  for(k = 0; k < NSOFTTYPES; k++)
                    mass_per_type[k] += ExtNodes[p].mass_per_type[k];

#ifdef ADAPTIVE_HYDRO_SOFTENING
                  if(maxhydrosofttype < Nodes[p].u.d.maxhydrosofttype)
                    maxhydrosofttype = Nodes[p].u.d.maxhydrosofttype;
                  if(minhydrosofttype > Nodes[p].u.d.minhydrosofttype)
                    minhydrosofttype = Nodes[p].u.d.minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
                }
              else if(p < Tree_MaxPart + Tree_MaxNodes + NTopleaves) /* a pseudo particle */
                {
                  /* nothing to be done here because the mass of the
                   *  pseudo-particle is still zero. This will be changed
                   * later.
                   */
                }
              else
                { /* an imported point */
                  int n = p - (Tree_MaxPart + Tree_MaxNodes + NTopleaves);

                  if(n >= Tree_NumPartImported)
                    terminate("n >= Tree_NumPartImported");

                  mass += Tree_Points[n].Mass;
                  s[0] += Tree_Points[n].Mass * Tree_Points[n].Pos[0];
                  s[1] += Tree_Points[n].Mass * Tree_Points[n].Pos[1];
                  s[2] += Tree_Points[n].Mass * Tree_Points[n].Pos[2];

                  /* Might not need the following routine */
                  if(All.ForceSoftening[maxsofttype] < All.ForceSoftening[Tree_Points[n].SofteningType])
                    maxsofttype = Tree_Points[n].SofteningType;

#ifdef MULTIPLE_NODE_SOFTENING
#ifdef ADAPTIVE_HYDRO_SOFTENING
                  mass_per_type[Tree_Points[n].Type == 0 ? 0 : Tree_Points[n].SofteningType] += Tree_Points[n].Mass;

                  if(Tree_Points[n].Type == 0)
                    {
                      if(maxhydrosofttype < Tree_Points[n].SofteningType)
                        maxhydrosofttype = Tree_Points[n].SofteningType;
                      if(minhydrosofttype > Tree_Points[n].SofteningType)
                        minhydrosofttype = Tree_Points[n].SofteningType;
                    }
#else  /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
                  mass_per_type[Tree_Points[n].SofteningType] += Tree_Points[n].Mass;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING #else */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
                }
            }
        }

      if(mass)
        {
          s[0] /= mass;
          s[1] /= mass;
          s[2] /= mass;
        }
      else
        {
          s[0] = Nodes[no].center[0];
          s[1] = Nodes[no].center[1];
          s[2] = Nodes[no].center[2];
        }

      Nodes[no].u.d.mass        = mass;
      Nodes[no].u.d.s[0]        = s[0];
      Nodes[no].u.d.s[1]        = s[1];
      Nodes[no].u.d.s[2]        = s[2];
      Nodes[no].u.d.maxsofttype = maxsofttype;
#ifdef MULTIPLE_NODE_SOFTENING
      int k;
      for(k = 0; k < NSOFTTYPES; k++)
        ExtNodes[no].mass_per_type[k] = mass_per_type[k];

#ifdef ADAPTIVE_HYDRO_SOFTENING
      Nodes[no].u.d.maxhydrosofttype = maxhydrosofttype;
      Nodes[no].u.d.minhydrosofttype = minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

      Nodes[no].u.d.sibling = sib;
      Nodes[no].u.d.father  = father;
    }
  else /* single particle or pseudo particle */
    {
      if(*last >= 0)
        {
          if(*last >= Tree_MaxPart)
            {
              if(*last >= Tree_MaxPart + Tree_MaxNodes)
                Nextnode[*last - Tree_MaxNodes] = no; /* a pseudo-particle or an imported point */
              else
                Nodes[*last].u.d.nextnode = no;
            }
          else
            Nextnode[*last] = no;
        }

      *last = no;

      if(no < Tree_MaxPart) /* only set it for single particles... */
        Father[no] = father;
      if(no >= Tree_MaxPart + Tree_MaxNodes + NTopleaves) /* ...or for imported points */
        Father[no - Tree_MaxNodes - NTopleaves] = father;
    }
}

/*! \brief Communicates the values of the multipole moments of the
 *         top-level tree-nodes of the domain grid.
 *
 *  This data can then be used to update the pseudo-particles on each CPU
 *  accordingly.
 *
 *  \return void
 */
void force_exchange_topleafdata(void)
{
  struct DomainNODE
  {
    MyDouble s[3];
    MyDouble mass;
#ifdef MULTIPLE_NODE_SOFTENING
    MyDouble mass_per_type[NSOFTTYPES];
#ifdef ADAPTIVE_HYDRO_SOFTENING
    unsigned char maxhydrosofttype;
    unsigned char minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
    unsigned char maxsofttype;
#if defined(SUBFIND) && defined(SUBFIND_EXTENDED_PROPERTIES)
    int NodeGrNr;
#endif /* #if defined(SUBFIND) && defined(SUBFIND_EXTENDED_PROPERTIES) */
  };

  struct DomainNODE *DomainMoment = (struct DomainNODE *)mymalloc("DomainMoment", NTopleaves * sizeof(struct DomainNODE));

  /* share the pseudo-particle data accross CPUs */
  int *recvcounts = (int *)mymalloc("recvcounts", sizeof(int) * NTask);
  int *recvoffset = (int *)mymalloc("recvoffset", sizeof(int) * NTask);
  int *bytecounts = (int *)mymalloc("bytecounts", sizeof(int) * NTask);
  int *byteoffset = (int *)mymalloc("byteoffset", sizeof(int) * NTask);

  for(int task = 0; task < NTask; task++)
    recvcounts[task] = 0;

  for(int n = 0; n < NTopleaves; n++)
    recvcounts[DomainNewTask[n]]++;

  for(int task = 0; task < NTask; task++)
    bytecounts[task] = recvcounts[task] * sizeof(struct DomainNODE);

  recvoffset[0] = 0, byteoffset[0] = 0;
  for(int task = 1; task < NTask; task++)
    {
      recvoffset[task] = recvoffset[task - 1] + recvcounts[task - 1];
      byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];
    }

  struct DomainNODE *loc_DomainMoment =
      (struct DomainNODE *)mymalloc("loc_DomainMoment", recvcounts[ThisTask] * sizeof(struct DomainNODE));

  int idx = 0;
  for(int n = 0; n < NTopleaves; n++)
    {
      if(DomainNewTask[n] == ThisTask)
        {
          int no = DomainNodeIndex[n];

          /* read out the multipole moments from the local base cells */
          loc_DomainMoment[idx].s[0]        = Nodes[no].u.d.s[0];
          loc_DomainMoment[idx].s[1]        = Nodes[no].u.d.s[1];
          loc_DomainMoment[idx].s[2]        = Nodes[no].u.d.s[2];
          loc_DomainMoment[idx].mass        = Nodes[no].u.d.mass;
          loc_DomainMoment[idx].maxsofttype = Nodes[no].u.d.maxsofttype;

#ifdef MULTIPLE_NODE_SOFTENING
          for(int k = 0; k < NSOFTTYPES; k++)
            loc_DomainMoment[idx].mass_per_type[k] = ExtNodes[no].mass_per_type[k];

#ifdef ADAPTIVE_HYDRO_SOFTENING
          loc_DomainMoment[idx].maxhydrosofttype = Nodes[no].u.d.maxhydrosofttype;
          loc_DomainMoment[idx].minhydrosofttype = Nodes[no].u.d.minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
          idx++;
        }
    }

  MPI_Allgatherv(loc_DomainMoment, bytecounts[ThisTask], MPI_BYTE, DomainMoment, bytecounts, byteoffset, MPI_BYTE, MPI_COMM_WORLD);

  for(int task = 0; task < NTask; task++)
    recvcounts[task] = 0;

  for(int n = 0; n < NTopleaves; n++)
    {
      int task = DomainNewTask[n];
      if(task != ThisTask)
        {
          int no  = DomainNodeIndex[n];
          int idx = recvoffset[task] + recvcounts[task]++;

          Nodes[no].u.d.s[0]        = DomainMoment[idx].s[0];
          Nodes[no].u.d.s[1]        = DomainMoment[idx].s[1];
          Nodes[no].u.d.s[2]        = DomainMoment[idx].s[2];
          Nodes[no].u.d.mass        = DomainMoment[idx].mass;
          Nodes[no].u.d.maxsofttype = DomainMoment[idx].maxsofttype;

#ifdef MULTIPLE_NODE_SOFTENING
          for(int k = 0; k < NSOFTTYPES; k++)
            ExtNodes[no].mass_per_type[k] = DomainMoment[idx].mass_per_type[k];
#ifdef ADAPTIVE_HYDRO_SOFTENING
          Nodes[no].u.d.maxhydrosofttype = DomainMoment[idx].maxhydrosofttype;
          Nodes[no].u.d.minhydrosofttype = DomainMoment[idx].minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
        }
    }

  myfree(loc_DomainMoment);
  myfree(byteoffset);
  myfree(bytecounts);
  myfree(recvoffset);
  myfree(recvcounts);
  myfree(DomainMoment);
}

/*! \brief Updates the top-level tree after the multipole moments of the
 *         pseudo-particles have been updated.
 *
 *  \param[in] no Node to be updated.
 *  \param[in] topnode Index of the node no in the 'TopNodes' array.
 *  \param[in] bits 2^bits is the number of nodes per dimension at the level of
 *             the daughter nodes of node no.
 *  \param[in] x Position of the node no in the x direction, falls in the
 *             range [0,2^(bits-1) - 1].
 *  \param[in] y Position of the node no in the y direction, falls in the
 *             range [0,2^(bits-1) - 1].
 *  \param[in] z Position of the node no in the z direction, falls in the
 *             range [0,2^(bits-1) - 1].
 *
 *  \return void
 */
void force_treeupdate_toplevel(int no, int topnode, int bits, int x, int y, int z)
{
  double s[3], mass;
  unsigned char maxsofttype;
#ifdef MULTIPLE_NODE_SOFTENING
  double mass_per_type[NSOFTTYPES];
#ifdef ADAPTIVE_HYDRO_SOFTENING
  unsigned char maxhydrosofttype;
  unsigned char minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

  if(TopNodes[topnode].Daughter >= 0)
    {
      for(int i = 0; i < 2; i++)
        for(int j = 0; j < 2; j++)
          for(int k = 0; k < 2; k++)
            {
              int sub = 7 & peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);

              Tree_NextFreeNode++;
              force_treeupdate_toplevel(Tree_NextFreeNode - 1, TopNodes[topnode].Daughter + sub, bits + 1, 2 * x + i, 2 * y + j,
                                        2 * z + k);
            }

      mass        = 0;
      s[0]        = 0;
      s[1]        = 0;
      s[2]        = 0;
      maxsofttype = NSOFTTYPES + NSOFTTYPES_HYDRO;
#ifdef MULTIPLE_NODE_SOFTENING
      for(int j = 0; j < NSOFTTYPES; j++)
        mass_per_type[j] = 0;

#ifdef ADAPTIVE_HYDRO_SOFTENING
      maxhydrosofttype = NSOFTTYPES;
      minhydrosofttype = NSOFTTYPES + NSOFTTYPES_HYDRO - 1;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */

      int p = Nodes[no].u.d.nextnode;

      for(int j = 0; j < 8; j++) /* since we are dealing with top-level nodes, we know that there are 8 consecutive daughter nodes */
        {
          if(p >= Tree_MaxPart && p < Tree_MaxPart + Tree_MaxNodes) /* internal node */
            {
              mass += Nodes[p].u.d.mass;
              s[0] += Nodes[p].u.d.mass * Nodes[p].u.d.s[0];
              s[1] += Nodes[p].u.d.mass * Nodes[p].u.d.s[1];
              s[2] += Nodes[p].u.d.mass * Nodes[p].u.d.s[2];

              if(All.ForceSoftening[maxsofttype] < All.ForceSoftening[Nodes[p].u.d.maxsofttype])
                maxsofttype = Nodes[p].u.d.maxsofttype;
#ifdef MULTIPLE_NODE_SOFTENING
              for(int k = 0; k < NSOFTTYPES; k++)
                mass_per_type[k] += ExtNodes[p].mass_per_type[k];

#ifdef ADAPTIVE_HYDRO_SOFTENING
              if(maxhydrosofttype < Nodes[p].u.d.maxhydrosofttype)
                maxhydrosofttype = Nodes[p].u.d.maxhydrosofttype;
              if(minhydrosofttype > Nodes[p].u.d.minhydrosofttype)
                minhydrosofttype = Nodes[p].u.d.minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
            }
          else
            terminate("may not happen");

          p = Nodes[p].u.d.sibling;
        }

      if(mass)
        {
          s[0] /= mass;
          s[1] /= mass;
          s[2] /= mass;
        }
      else
        {
          s[0] = Nodes[no].center[0];
          s[1] = Nodes[no].center[1];
          s[2] = Nodes[no].center[2];
        }

      Nodes[no].u.d.s[0]        = s[0];
      Nodes[no].u.d.s[1]        = s[1];
      Nodes[no].u.d.s[2]        = s[2];
      Nodes[no].u.d.mass        = mass;
      Nodes[no].u.d.maxsofttype = maxsofttype;
#ifdef MULTIPLE_NODE_SOFTENING
      for(int k = 0; k < NSOFTTYPES; k++)
        ExtNodes[no].mass_per_type[k] = mass_per_type[k];
#ifdef ADAPTIVE_HYDRO_SOFTENING
      Nodes[no].u.d.maxhydrosofttype = maxhydrosofttype;
      Nodes[no].u.d.minhydrosofttype = minhydrosofttype;
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
    }
}

/*! \brief Allocates the memory used for storage of the tree nodes.
 *
 *  Usually, the number of required nodes is of order 0.7*maxpart, but if this
 *  is insufficient, the code will try to allocated more space.
 *
 *  \param[in] maxpart Number of particles on the current task.
 *  \param[in] maxindex The Nodes pointer will be shifted such that the index
 *             of the first element is maxindex.
 *
 *  \return void
 */
void force_treeallocate(int maxpart, int maxindex)
{
  if(Nodes)
    terminate("already allocated");

  Tree_MaxPart  = maxindex;
  Tree_MaxNodes = (int)(All.TreeAllocFactor * maxpart) + NTopnodes;

  DomainNewTask   = (int *)mymalloc_movable(&DomainNewTask, "DomainNewTask", NTopleaves * sizeof(int));
  DomainNodeIndex = (int *)mymalloc_movable(&DomainNodeIndex, "DomainNodeIndex", NTopleaves * sizeof(int));
  Tree_Task_list  = (int *)mymalloc_movable(&Tree_Task_list, "Tree_Task_list", maxpart * sizeof(int));
  Tree_Pos_list   = (MyDouble *)mymalloc_movable(&Tree_Pos_list, "Tree_Pos_list", 3 * maxpart * sizeof(MyDouble));

  Nodes = (struct NODE *)mymalloc_movable(&Nodes, "Nodes", (Tree_MaxNodes + 1) * sizeof(struct NODE));
  Nodes -= Tree_MaxPart;
#ifdef MULTIPLE_NODE_SOFTENING
  ExtNodes = (struct ExtNODE *)mymalloc_movable(&ExtNodes, "ExtNodes", (Tree_MaxNodes + 1) * sizeof(struct ExtNODE));
  ExtNodes -= Tree_MaxPart;
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
}

/*! \brief Frees the memory allocated for the tree.
 *
 *  I.e. it frees the space allocated by the function force_treeallocate().
 *
 *  \return void
 */
void force_treefree(void)
{
  if(Nodes)
    {
#ifdef MULTIPLE_NODE_SOFTENING
      myfree(ExtNodes + Tree_MaxPart);
      ExtNodes = NULL;
#endif /* #ifdef MULTIPLE_NODE_SOFTENING */
      myfree(Nodes + Tree_MaxPart);
      myfree(Tree_Pos_list);
      myfree(Tree_Task_list);
      myfree(DomainNodeIndex);
      myfree(DomainNewTask);

      Nodes           = NULL;
      DomainNodeIndex = NULL;
      DomainNewTask   = NULL;
      Tree_Task_list  = NULL;
      Nextnode        = NULL;
      Father          = NULL;
    }
  else
    terminate("trying to free the tree even though it's not allocated");
}

/*! \brief Dump particle data into file.
 *
 *  This function dumps some of the basic particle data to a file. In case
 *  the tree construction fails, it is called just before the run
 *  terminates with an error message. Examination of the generated file may
 *  then give clues to what caused the problem.
 *
 *  \return void
 */
void dump_particles(void)
{
  char buffer[200];
  sprintf(buffer, "particles%d.dat", ThisTask);
  FILE *fd = fopen(buffer, "w");
  my_fwrite(&NumPart, 1, sizeof(int), fd);
  for(int i = 0; i < NumPart; i++)
    my_fwrite(&P[i].Pos[0], 3, sizeof(MyDouble), fd);
  for(int i = 0; i < NumPart; i++)
    my_fwrite(&P[i].Vel[0], 3, sizeof(MyFloat), fd);
  for(int i = 0; i < NumPart; i++)
    my_fwrite(&P[i].ID, 1, sizeof(int), fd);
  fclose(fd);
}

#ifdef ADDBACKGROUNDGRID
/*! \brief Add additional empty nodes.
 *
 *  Called during tree construction if ADDBACKGROUNDGRID is active.
 *
 *  \return 0: default; 1: number of nodes > max number of nodes.
 */
int force_add_empty_nodes(void)
{
  int nempty = 0;
  int no, j, subnode;

  for(no = Tree_MaxPart; no < Tree_MaxPart + Tree_NumNodes; no++)
    {
      int count = 0;

      for(subnode = 0; subnode < 8; subnode++)
        if(Nodes[no].u.suns[subnode] == -1)
          count++;

      if(count < 8)
        {
          for(subnode = 0, count = 0; subnode < 8; subnode++)
            if(Nodes[no].u.suns[subnode] == -1)
              {
                Nodes[no].u.suns[subnode] = Tree_NextFreeNode;
                struct NODE *nfreep       = &Nodes[Tree_NextFreeNode];

                nfreep->len    = 0.5 * Nodes[no].len;
                double lenhalf = 0.25 * Nodes[no].len;

                if(subnode & 1)
                  nfreep->center[0] = Nodes[no].center[0] + lenhalf;
                else
                  nfreep->center[0] = Nodes[no].center[0] - lenhalf;

                if(subnode & 2)
                  nfreep->center[1] = Nodes[no].center[1] + lenhalf;
                else
                  nfreep->center[1] = Nodes[no].center[1] - lenhalf;

                if(subnode & 4)
                  nfreep->center[2] = Nodes[no].center[2] + lenhalf;
                else
                  nfreep->center[2] = Nodes[no].center[2] - lenhalf;

                for(j = 0; j < 8; j++)
                  nfreep->u.suns[j] = -1;

                Tree_NumNodes++;
                Tree_NextFreeNode++;

                if(Tree_NumNodes >= Tree_MaxNodes)
                  {
                    if(All.TreeAllocFactor > 5.0)
                      {
                        char buf[500];
                        sprintf(
                            buf,
                            "task %d: looks like a serious problem, stopping with particle dump. Tree_NumNodes=%d Tree_MaxNodes=%d\n",
                            ThisTask, Tree_NumNodes, Tree_MaxNodes);
                        dump_particles();
                        terminate(buf);
                      }
                    return 1;
                  }
                nempty++;
              }
        }
    }

  printf("FORCETREE: Task %d has added %d empty nodes\n", ThisTask, nempty);
  return 0;
}
#endif /* #ifdef ADDBACKGROUNDGRID */
