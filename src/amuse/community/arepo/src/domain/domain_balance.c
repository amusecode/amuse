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
 * \file        src/domain/domain_balance.c
 * \date        05/2018
 * \brief       Load-balancing algorithms.
 * \details     Algorithms to estimate cost of different particles and cells
 *              and to balance the workload and memory usage equally over the
 *              mpi tasks.
 *              contains functions:
 *                double domain_grav_tot_costfactor(int i)
 *                double domain_hydro_tot_costfactor(int i)
 *                void domain_init_sum_cost(void)
 *                void domain_sumCost(void)
 *                void domain_combine_topleaves_to_domains(int ncpu, int
 *                  ndomain)
 *                int domain_sort_task(const void *a, const void *b)
 *                int domain_sort_load(const void *a, const void *b)
 *                static int mydata_cmp(struct mydata *lhs, struct mydata *rhs)
 *                void domain_combine_multipledomains(void)
 *                void domain_optimize_domain_to_task_mapping(void)
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

#include "../domain/bsd_tree.h"
#include "../domain/domain.h"
#include "../mesh/voronoi/voronoi.h"

/* do some preparation work for use of red-black ordered binary tree based on BSD macros */

/*! \brief Defines structure of mytree nodes.
 */
struct mydata
{
  double pri;
  int target;
  RB_ENTRY(mydata) linkage; /* this creates the linkage pointers needed by the RB tree, using symbolic name 'linkage' */
};

/* prototype of comparison function of tree elements */
static int mydata_cmp(struct mydata *lhs, struct mydata *rhs);

/* the following macro declares 'struct mytree', which is the header element needed as handle for a tree */
RB_HEAD(mytree, mydata);

/* the following macros declare appropriate function prototypes and functions needed for this type of tree */
RB_PROTOTYPE_STATIC(mytree, mydata, linkage, mydata_cmp);
RB_GENERATE_STATIC(mytree, mydata, linkage, mydata_cmp);

/*! \brief Computes gravity cost.
 *
 *  All timebins in which the particle appears are summed, and the relative
 *  frequency with which this timebin is executed is taken into account.
 *
 *  \param[in] i Index of cell in P and SphP array.
 *
 *  \return cost-factor.
 */
double domain_grav_tot_costfactor(int i)
{
  double w = MIN_FLOAT_NUMBER;

#ifdef SELFGRAVITY
  for(int bin = All.LowestOccupiedTimeBin; bin <= All.HighestActiveTimeBin; bin++)
    {
      if(domain_to_be_balanced[bin])
        {
#ifdef HIERARCHICAL_GRAVITY
          if(bin >= P[i].TimeBinGrav)
#endif /* #ifdef HIERARCHICAL_GRAVITY */
            {
              if(domain_bintolevel[bin] >= 0)
                w += domain_grav_weight[bin] * P[i].GravCost[domain_bintolevel[bin]];
              else
                {
                  if(domain_refbin[bin] >= 0)
                    w += domain_grav_weight[bin] * P[i].GravCost[domain_bintolevel[domain_refbin[bin]]];
                  else
                    w += domain_grav_weight[bin];
                }
            }
        }
    }
#endif /* #ifdef SELFGRAVITY */

  return w;
}

/*! \brief Computes hydro cost.
 *
 *  If a cell is active on a certain timebin, it is assigned a cost of "1".
 *  All active timebins are summed, and the frequency with which each timebin
 *  is executed is taken into account.
 *
 *  \param[in] i Index of cell in P and SphP array.
 *
 *  \return cost-factor.
 */
double domain_hydro_tot_costfactor(int i)
{
  double w = 0;

  if(P[i].Type == 0)
    for(int bin = P[i].TimeBinHydro; bin <= All.HighestOccupiedTimeBin; bin++)
      if(domain_to_be_balanced[bin])
        w += domain_hydro_weight[bin];

  return w;
}

/*! \brief Prepares cost measurement.
 *
 *  This function prepares the measurement of the total cost on each domain.
 *  In particular, we determine how the timebins are mapped to the explicit
 *  measurements of the gravity cost stored in the P.GravCost[] array (which
 *  in general will only be available for a subset of all timebins). For the
 *  unmatched timebins, a closest bin is selected that is the most similar in
 *  terms of particle number on the bin. Finally, the routine also determines
 *  how often each timebin is executed in one cycle associated with the
 *  highest occupied timebin.
 *
 *  \return void
 */
void domain_init_sum_cost(void)
{
  long long tot_count[TIMEBINS], tot_count_sph[TIMEBINS];

  sumup_large_ints(TIMEBINS, TimeBinsGravity.TimeBinCount, tot_count);
  sumup_large_ints(TIMEBINS, TimeBinsHydro.TimeBinCount, tot_count_sph);

  for(int i = 0; i < TIMEBINS; i++)
    {
      domain_bintolevel[i] = -1;
      domain_refbin[i]     = -1;
    }

  for(int j = 0; j < GRAVCOSTLEVELS; j++) /* bins that have known levels at this point */
    if(All.LevelToTimeBin[j] >= 0)
      domain_bintolevel[All.LevelToTimeBin[j]] = j;

  for(int i = 0; i < TIMEBINS; i++)
    if(tot_count[i] > 0 && domain_bintolevel[i] < 0) /* need to find a reference bin for this one */
      {
        double mindiff = MAX_REAL_NUMBER;
        int ref_bin    = -1;
        for(int j = 0; j < TIMEBINS; j++)
          if(domain_bintolevel[j] >= 0 && tot_count[j] > 0)
            {
              if(mindiff > llabs(tot_count[i] - tot_count[j]))
                {
                  mindiff = llabs(tot_count[i] - tot_count[j]);
                  ref_bin = j;
                }
            }

        if(ref_bin >= 0)
          domain_refbin[i] = ref_bin;
      }

  for(int i = 0; i < TIMEBINS; i++)
    {
      domain_to_be_balanced[i] = 0;
      domain_grav_weight[i]    = 1;
      domain_hydro_weight[i]   = 1;
    }

#ifdef HIERARCHICAL_GRAVITY

  domain_to_be_balanced[All.HighestActiveTimeBin] = 1;
  domain_grav_weight[All.HighestActiveTimeBin]    = 1;
  domain_hydro_weight[All.HighestActiveTimeBin]   = 1;

  for(int j = All.HighestActiveTimeBin - 1; j >= All.LowestOccupiedTimeBin; j--)
    {
      if(tot_count[j] > 0 || tot_count_sph[j] > 0)
        domain_to_be_balanced[j] = 1;

      domain_grav_weight[j] += 2;
    }

  for(int i = All.SmallestTimeBinWithDomainDecomposition - 1, weight = 1; i >= All.LowestOccupiedTimeBin; i--, weight *= 2)
    {
      if(tot_count[i] > 0)
        {
          domain_grav_weight[i] = weight;

          for(int j = i - 1; j >= All.LowestOccupiedTimeBin; j--)
            domain_grav_weight[j] += 2 * weight;
        }

      if(tot_count_sph[i] > 0)
        domain_hydro_weight[i] = weight;
    }

#else /* #ifdef HIERARCHICAL_GRAVITY */

  domain_to_be_balanced[All.HighestActiveTimeBin] = 1;
  domain_grav_weight[All.HighestActiveTimeBin]    = 1;
  domain_hydro_weight[All.HighestActiveTimeBin]   = 1;

  for(int i = All.SmallestTimeBinWithDomainDecomposition - 1, weight = 1; i >= All.LowestOccupiedTimeBin; i--, weight *= 2)
    {
      if(tot_count[i] > 0 || tot_count_sph[i] > 0)
        domain_to_be_balanced[i] = 1;

      if(tot_count[i] > 0)
        domain_grav_weight[i] = weight;

      if(tot_count_sph[i] > 0)
        domain_hydro_weight[i] = weight;
    }

#endif /* #ifdef HIERARCHICAL_GRAVITY #else */
}

/*! \brief Determine cost and load
 *
 *  This function determines the cost and load associated with each top-level
 *  leaf node of the tree. These leave nodes can be distributed among the
 *  processors in order to reach a good work-load and memory-load balance.
 *
 *  \return void
 */
void domain_sumCost(void)
{
  int i, j, n, no, nexport = 0, nimport = 0, ngrp, task, loc_first_no;

  struct domain_cost_data *loc_DomainLeaveNode, *listCost, *export_node_data, *import_node_data;

  int *blocksize = mymalloc("blocksize", sizeof(int) * NTask);
  int blk        = NTopleaves / NTask;
  int rmd        = NTopleaves - blk * NTask; /* remainder */
  int pivot_no   = rmd * (blk + 1);

  for(task = 0, loc_first_no = 0; task < NTask; task++)
    {
      if(task < rmd)
        blocksize[task] = blk + 1;
      else
        blocksize[task] = blk;

      if(task < ThisTask)
        loc_first_no += blocksize[task];
    }

  loc_DomainLeaveNode = mymalloc("loc_DomainLeaveNode", blocksize[ThisTask] * sizeof(struct domain_cost_data));
  memset(loc_DomainLeaveNode, 0, blocksize[ThisTask] * sizeof(struct domain_cost_data));

  listCost = mymalloc("listCost", NTopleaves * sizeof(struct domain_cost_data));

  int *no_place = mymalloc("no_place", NTopleaves * sizeof(int));
  memset(no_place, -1, NTopleaves * sizeof(int));

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  /* find for each particle its top-leave, and then add the associated cost with it */
  for(n = 0; n < NumPart; n++)
    {
#ifdef ADDBACKGROUNDGRID
      if(P[n].Type != 0)
        continue;
#endif /* #ifdef ADDBACKGROUNDGRID */
      no = 0;

      peanokey mask = ((peanokey)7) << (3 * (BITS_PER_DIMENSION - 1));
      int shift     = 3 * (BITS_PER_DIMENSION - 1);

      while(topNodes[no].Daughter >= 0)
        {
          no = topNodes[no].Daughter + (int)((Key[n] & mask) >> shift);
          mask >>= 3;
          shift -= 3;
        }

      no = topNodes[no].Leaf;

      int p = no_place[no];
      if(p < 0)
        {
          p            = nexport++;
          no_place[no] = p;

          memset(&listCost[p], 0, sizeof(struct domain_cost_data));
          listCost[p].no = no;

          if(no < pivot_no)
            task = no / (blk + 1);
          else
            task = rmd + (no - pivot_no) / blk; /* note: if blk=0, then this case can not occur, since then always no < pivot_no */

          if(task < 0 || task > NTask)
            terminate("task < 0 || task > NTask");

          Send_count[task]++;
        }

      listCost[p].Count += 1;
      listCost[p].Work += domain_grav_tot_costfactor(n);
      listCost[p].WorkSph += domain_hydro_tot_costfactor(n);

      if(P[n].Type == 0)
        listCost[p].CountSph += 1;
    }

  myfree(no_place);

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];
      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  export_node_data = mymalloc("export_node_data", nexport * sizeof(struct domain_cost_data));
  import_node_data = mymalloc("import_node_data", nimport * sizeof(struct domain_cost_data));

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < nexport; i++)
    {
      if(listCost[i].no < pivot_no)
        task = listCost[i].no / (blk + 1);
      else
        task = rmd +
               (listCost[i].no - pivot_no) / blk; /* note: if blk=0, then this case can not occur, since then always no < pivot_no */

      int ind               = Send_offset[task] + Send_count[task]++;
      export_node_data[ind] = listCost[i];
    }

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++) /* note: here we also have a transfer from each task to itself (for ngrp=0) */
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&export_node_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct domain_cost_data), MPI_BYTE,
                       recvTask, TAG_DENS_B, &import_node_data[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(struct domain_cost_data), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
    }

  for(i = 0; i < nimport; i++)
    {
      int j = import_node_data[i].no - loc_first_no;

      if(j < 0 || j >= blocksize[ThisTask])
        terminate("j=%d < 0 || j>= blocksize[ThisTask]=%d   loc_first_no=%d  import_node_data[i].no=%d  i=%d  nimport=%d", j,
                  blocksize[ThisTask], loc_first_no, import_node_data[i].no, i, nimport);

      loc_DomainLeaveNode[j].Count += import_node_data[i].Count;
      loc_DomainLeaveNode[j].Work += import_node_data[i].Work;
      loc_DomainLeaveNode[j].CountSph += import_node_data[i].CountSph;
      loc_DomainLeaveNode[j].WorkSph += import_node_data[i].WorkSph;
    }

  myfree(import_node_data);
  myfree(export_node_data);

  /* now share the cost data across all processors */
  int *bytecounts = (int *)mymalloc("bytecounts", sizeof(int) * NTask);
  int *byteoffset = (int *)mymalloc("byteoffset", sizeof(int) * NTask);

  for(task = 0; task < NTask; task++)
    bytecounts[task] = blocksize[task] * sizeof(struct domain_cost_data);

  for(task = 1, byteoffset[0] = 0; task < NTask; task++)
    byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];

  MPI_Allgatherv(loc_DomainLeaveNode, bytecounts[ThisTask], MPI_BYTE, DomainLeaveNode, bytecounts, byteoffset, MPI_BYTE,
                 MPI_COMM_WORLD);

  myfree(byteoffset);
  myfree(bytecounts);
  myfree(listCost);
  myfree(loc_DomainLeaveNode);
  myfree(blocksize);
}

/*! \brief Uses cost function to combine top-level nodes to domains.
 *
 *  This function uses the cumulative cost function (which weights work-load
 *  and memory-load equally) to subdivide the list of top-level leave nodes
 *  into pieces that are (approximately) equal in size.
 *
 *  \param[in] ncpu Number of chunks/damains.
 *  \param[in] ndomain Number of topleaves.
 *
 *  \return void
 */
void domain_combine_topleaves_to_domains(int ncpu, int ndomain)
{
  double t0 = second();

  double max_work     = 0;
  double workhalfnode = 0.5 / ndomain;
  double workavg      = 1.0 / ncpu;
  double work_before = 0, workavg_before = 0;
  int start = 0;

  int nabove_grav = 0, nabove_sph = 0;
  double todistribute_grav = 0.0;
  double todistribute_sph  = 0.0;
  double weightsum_grav    = 0.0;
  double weightsum_sph     = 0.0;

  for(int i = 0; i < ndomain; i++)
    {
      if(fac_work * DomainLeaveNode[i].Work > normsum_work / ncpu)
        {
          nabove_grav++;
          todistribute_grav += DomainLeaveNode[i].Work - normsum_work / ncpu / fac_work;
        }
      else
        weightsum_grav += DomainLeaveNode[i].Count;

      if(fac_worksph * DomainLeaveNode[i].WorkSph > normsum_worksph / ncpu)
        {
          nabove_sph++;
          todistribute_sph += DomainLeaveNode[i].WorkSph - normsum_worksph / ncpu / fac_worksph;
        }
      else
        weightsum_sph += DomainLeaveNode[i].Count;
    }

  struct leafnode_data
  {
    double workgrav;
    double worksph;
  };

  struct leafnode_data *leaf = (struct leafnode_data *)mymalloc("leaf", ndomain * sizeof(struct leafnode_data));

  for(int i = 0; i < ndomain; i++)
    {
      leaf[i].workgrav = DomainLeaveNode[i].Work;
      leaf[i].worksph  = DomainLeaveNode[i].WorkSph;

      if(fac_work > 0 && weightsum_grav > 0)
        {
          if(fac_work * DomainLeaveNode[i].Work > normsum_work / ncpu)
            leaf[i].workgrav = normsum_work / ncpu / fac_work;
          else
            leaf[i].workgrav += (DomainLeaveNode[i].Count / weightsum_grav) * todistribute_grav;
        }

      if(fac_worksph > 0 && weightsum_sph > 0)
        {
          if(fac_worksph * DomainLeaveNode[i].WorkSph > normsum_worksph / ncpu)
            leaf[i].worksph = normsum_worksph / ncpu / fac_worksph;
          else
            leaf[i].worksph += (DomainLeaveNode[i].Count / weightsum_sph) * todistribute_sph;
        }
    }

  for(int i = 0; i < ncpu; i++)
    {
      double work = 0;
      int end     = start;

      work += fac_work * leaf[end].workgrav + fac_load * DomainLeaveNode[end].Count + fac_worksph * leaf[end].worksph;

      while((work + work_before +
                 (end + 1 < ndomain ? fac_work * leaf[end + 1].workgrav + fac_load * DomainLeaveNode[end + 1].Count +
                                          fac_worksph * leaf[end + 1].worksph
                                    : 0) <
             workavg + workavg_before + workhalfnode) ||
            (i == ncpu - 1 && end < ndomain - 1))
        {
          if((ndomain - end) > (ncpu - i))
            end++;
          else
            break;

          work += fac_work * leaf[end].workgrav + fac_load * DomainLeaveNode[end].Count + fac_worksph * leaf[end].worksph;
        }

      DomainStartList[i] = start;
      DomainEndList[i]   = end;

      work_before += work;
      workavg_before += workavg;
      start = end + 1;

      if(max_work < work)
        max_work = work;
    }

  myfree(leaf);

  double t1 = second();
  mpi_printf("DOMAIN: balance reached among multiple-domains=%g, average leave-nodes per domain=%g  (took %g sec)\n",
             max_work / workavg, ((double)ndomain) / ncpu, timediff(t0, t1));
}

/*! \brief Structure containing data for segments.
 */
static struct domain_segments_data
{
  int task, start, end;
  double bin_GravCost[TIMEBINS];
  double bin_HydroCost[TIMEBINS];
  double work;
  double load;
  double worksph;
  double normalized_load;
} * domainAssign;

/*! \brief Structure containing data for task list.
 */
struct tasklist_data
{
  double bin_GravCost[TIMEBINS];
  double bin_HydroCost[TIMEBINS];
  double work;
  double load;
  double worksph;
  int count;
} * tasklist;

/*! \brief Comparison function for domain_segments_data structure.
 *
 *  Compares field task.
 *
 *  \param a Pointer to fist object.
 *  \param b Pointer to second object.
 *
 *  \return (-1,0,1); -1 if a < b.
 */
int domain_sort_task(const void *a, const void *b)
{
  if(((struct domain_segments_data *)a)->task < (((struct domain_segments_data *)b)->task))
    return -1;

  if(((struct domain_segments_data *)a)->task > (((struct domain_segments_data *)b)->task))
    return +1;

  return 0;
}

/*! \brief Comparison functions for domain_segmens_data structures.
 *
 *  Compares field normalized_load.
 *
 *  \param a Pointer to fist object.
 *  \param b Pointer to second object.
 *
 *  \return (-1,0,1) -1 if a>b.
 */
int domain_sort_load(const void *a, const void *b)
{
  if(((struct domain_segments_data *)a)->normalized_load > (((struct domain_segments_data *)b)->normalized_load))
    return -1;

  if(((struct domain_segments_data *)a)->normalized_load < (((struct domain_segments_data *)b)->normalized_load))
    return +1;

  return 0;
}

/*! \brief Comparison function for objects of type mydata.
 *
 *  Compares elements pri and target.
 *
 *  \param lhs Pointer to fist object.
 *  \param rhs Pointer to second object.
 *
 *  \return (-1,0,1); -1 if lhs < rhs.
 */
static int mydata_cmp(struct mydata *lhs, struct mydata *rhs)
{
  if(lhs->pri < rhs->pri)
    return -1;
  else if(lhs->pri > rhs->pri)
    return 1;
  else if(lhs->target < rhs->target)
    return -1;
  else if(lhs->target > rhs->target)
    return 1;

  return 0;
}

/*! \brief Assigns the domain pieces to individual MPI tasks with the goal to
 *         balance the work-load on different timebins.
 *
 *  The algorithm used works as follows:
 *  The domains are assigned to the CPUs in sequence of decreasing "effective
 *  load", which is a simple combined measure of relative total gravity, hydro
 *  and memory load. For each assignment, a number of possible target CPUs are
 *  evaluated, and the assignment leading to the lowest total runtime is
 *  adopted. The set of target CPUs that is tested in each step is the one
 *  that consists of the CPUs that currently have the lowest load in the set
 *  of primary tasks that are examined.
 *
 *  \return void
 */
void domain_combine_multipledomains(void)
{
  double t0 = second();

  int ndomains = All.MultipleDomains * NTask;

  domainAssign = (struct domain_segments_data *)mymalloc("domainAssign", ndomains * sizeof(struct domain_segments_data));

  tasklist = mymalloc("tasklist", NTask * sizeof(struct tasklist_data));

  for(int ta = 0; ta < NTask; ta++)
    {
      tasklist[ta].load    = 0;
      tasklist[ta].work    = 0;
      tasklist[ta].worksph = 0;
      tasklist[ta].count   = 0;

      for(int i = 0; i < TIMEBINS; i++)
        {
          tasklist[ta].bin_GravCost[i]  = 0;
          tasklist[ta].bin_HydroCost[i] = 0;
        }
    }

  for(int n = 0; n < ndomains; n++)
    for(int i = DomainStartList[n]; i <= DomainEndList[n]; i++)
      DomainTask[i] = n;

  /* we first determine the grav-cost and hydro-cost separately for each
   * timebin of all the domain-pieces that are available for a
   * mapping to individual MPI tasks
   */

  struct cost_data
  {
    double GravCost;
    double HydroCost;
  } * loc_bin_Cost, *glob_bin_Cost;

  loc_bin_Cost  = mymalloc_clear("loc_bin_Cost", sizeof(struct cost_data) * ndomains * TIMEBINS);
  glob_bin_Cost = mymalloc_clear("glob_bin_Cost", sizeof(struct cost_data) * ndomains * TIMEBINS);

  for(int i = 0; i < NumPart; i++)
    {
#ifdef ADDBACKGROUNDGRID
      if(P[i].Type != 0)
        continue;
#endif /* #ifdef ADDBACKGROUNDGRID */
      int no = 0;

      peanokey mask = ((peanokey)7) << (3 * (BITS_PER_DIMENSION - 1));
      int shift     = 3 * (BITS_PER_DIMENSION - 1);

      while(topNodes[no].Daughter >= 0)
        {
          no = topNodes[no].Daughter + (int)((Key[i] & mask) >> shift);
          mask >>= 3;
          shift -= 3;
        }

      no = topNodes[no].Leaf;

      int n = DomainTask[no];

#ifdef SELFGRAVITY
      for(int bin = All.LowestOccupiedTimeBin; bin <= All.HighestActiveTimeBin; bin++)
        {
          if(domain_to_be_balanced[bin])
            {
#ifdef HIERARCHICAL_GRAVITY
              if(bin >= P[i].TimeBinGrav)
#endif /* #ifdef HIERARCHICAL_GRAVITY */
                {
                  if(domain_bintolevel[bin] >= 0)
                    loc_bin_Cost[bin * ndomains + n].GravCost +=
                        MIN_FLOAT_NUMBER + domain_grav_weight[bin] * P[i].GravCost[domain_bintolevel[bin]];
                  else
                    {
                      if(domain_refbin[bin] >= 0)
                        loc_bin_Cost[bin * ndomains + n].GravCost +=
                            MIN_FLOAT_NUMBER + domain_grav_weight[bin] * P[i].GravCost[domain_bintolevel[domain_refbin[bin]]];
                      else
                        loc_bin_Cost[bin * ndomains + n].GravCost += domain_grav_weight[bin];
                    }
                }
            }
        }
#endif /* #ifdef SELFGRAVITY */

      if(P[i].Type == 0)
        {
          for(int bin = P[i].TimeBinHydro; bin <= All.HighestActiveTimeBin; bin++)
            if(domain_to_be_balanced[bin])
              loc_bin_Cost[bin * ndomains + n].HydroCost += domain_hydro_weight[bin];
        }
    }

  allreduce_sparse_double_sum((double *)(loc_bin_Cost + All.LowestOccupiedTimeBin * ndomains),
                              (double *)(glob_bin_Cost + All.LowestOccupiedTimeBin * ndomains),
                              2 * ndomains * (All.HighestOccupiedTimeBin - All.LowestOccupiedTimeBin + 1));

  /* now assign this cost to the domainAssign-structure, which keeps track of the different pieces */
  double tot_work    = 0;
  double tot_load    = 0;
  double tot_worksph = 0;

  for(int n = 0; n < ndomains; n++)
    {
      domainAssign[n].start   = DomainStartList[n];
      domainAssign[n].end     = DomainEndList[n];
      domainAssign[n].work    = 0;
      domainAssign[n].load    = 0;
      domainAssign[n].worksph = 0;

      for(int i = 0; i < TIMEBINS; i++)
        {
          domainAssign[n].bin_GravCost[i]  = glob_bin_Cost[i * ndomains + n].GravCost;
          domainAssign[n].bin_HydroCost[i] = glob_bin_Cost[i * ndomains + n].HydroCost;
        }

      for(int i = DomainStartList[n]; i <= DomainEndList[n]; i++)
        {
          domainAssign[n].work += DomainLeaveNode[i].Work;
          domainAssign[n].load += DomainLeaveNode[i].Count;
          domainAssign[n].worksph += DomainLeaveNode[i].WorkSph;
        }

      tot_work += domainAssign[n].work;
      tot_load += domainAssign[n].load;
      tot_worksph += domainAssign[n].worksph;
    }

  for(int n = 0; n < ndomains; n++)
    {
      domainAssign[n].normalized_load = domainAssign[n].work / (tot_work + MIN_FLOAT_NUMBER) +
                                        domainAssign[n].worksph / (tot_worksph + MIN_FLOAT_NUMBER) +
                                        domainAssign[n].load / ((double)tot_load + MIN_FLOAT_NUMBER);
    }

  myfree(glob_bin_Cost);
  myfree(loc_bin_Cost);

  /* sort the pieces according to their normalized work-load, with the most heavily loaded coming first */
  mysort(domainAssign, ndomains, sizeof(struct domain_segments_data), domain_sort_load);

  /* initialize a structure that stores the maximum gravity and hydro cost load for each timebin */
  double max_GravCost[TIMEBINS], max_HydroCost[TIMEBINS];
  for(int i = 0; i < TIMEBINS; i++)
    {
      max_GravCost[i]  = 0;
      max_HydroCost[i] = 0;
    }

  double max_load = 0;

  /* create priority trees, one for the cost of each occupied timebin,
   * one for the hydro cost of each occupied timebin */
  struct mytree queue_gravcost[TIMEBINS];
  struct mytree queue_hydrocost[TIMEBINS];
  struct mytree queue_load;
  struct mydata *ngrav[TIMEBINS];
  struct mydata *nhydro[TIMEBINS];
  struct mydata *nload;

  for(int bin = All.LowestOccupiedTimeBin; bin <= All.HighestOccupiedTimeBin; bin++)
    {
      if(domain_to_be_balanced[bin])
        {
          RB_INIT(&queue_gravcost[bin]);
          ngrav[bin] = mymalloc("ngrav[bin]", NTask * sizeof(struct mydata));

          RB_INIT(&queue_hydrocost[bin]);
          nhydro[bin] = mymalloc("nhydro[bin]", NTask * sizeof(struct mydata));
        }
    }

  RB_INIT(&queue_load);
  nload = mymalloc("nload", NTask * sizeof(struct mydata));
  for(int i = 0; i < NTask; i++)
    {
      nload[i].pri    = 0;
      nload[i].target = i;
      RB_INSERT(mytree, &queue_load, &nload[i]);
    }

  /* fill in all the tasks into each queue. The priority will be the current cost of the bin, the tag 'val' is used to label the task
   */
  for(int bin = All.LowestOccupiedTimeBin; bin <= All.HighestOccupiedTimeBin; bin++)
    {
      if(!domain_to_be_balanced[bin])
        continue;

      for(int i = 0; i < NTask; i++)
        {
          ngrav[bin][i].pri    = 0;
          ngrav[bin][i].target = i;
          RB_INSERT(mytree, &queue_gravcost[bin], &ngrav[bin][i]);

          nhydro[bin][i].pri    = 0;
          nhydro[bin][i].target = i;
          RB_INSERT(mytree, &queue_hydrocost[bin], &nhydro[bin][i]);
        }
    }

  int n_lowest = MAX_FIRST_ELEMENTS_CONSIDERED;
  if(n_lowest > NTask)
    n_lowest = NTask;

  int rep, *candidates = mymalloc("candidates", n_lowest * sizeof(int));
  struct mydata *np;

  /* now assign each of the domains to a CPU, trying to minimize the overall runtime */
  for(int n = 0; n < ndomains; n++)
    {
      double best_runtime = MAX_FLOAT_NUMBER;
      int best_target     = -1;

      for(int bin = All.LowestOccupiedTimeBin; bin <= All.HighestOccupiedTimeBin; bin++)
        {
          if(!domain_to_be_balanced[bin])
            continue;

          int target;

          for(int set = 0; set < 2; set++)
            {
              if(set == 0)
                {
#ifndef SELFGRAVITY
                  continue;
#endif /* #ifndef SELFGRAVITY */
                  /* look up the n_lowest smallest elements from the tree */
                  for(np = RB_MIN(mytree, &queue_gravcost[bin]), rep = 0; np != NULL && rep < n_lowest;
                      np = RB_NEXT(mytree, &queue_gravcost[bin], np), rep++)
                    candidates[rep] = np->target;
                }
              else
                {
                  for(np = RB_MIN(mytree, &queue_hydrocost[bin]), rep = 0; np != NULL && rep < n_lowest;
                      np = RB_NEXT(mytree, &queue_hydrocost[bin], np), rep++)
                    candidates[rep] = np->target;
                }

              for(rep = 0; rep < n_lowest; rep++)
                {
                  target = candidates[rep];

                  double runtime = 0;

                  for(int i = 0; i < TIMEBINS; i++)
                    {
                      double sum = domainAssign[n].bin_GravCost[i] + tasklist[target].bin_GravCost[i];
                      if(sum < max_GravCost[i])
                        sum = max_GravCost[i];

                      runtime += sum / (totgravcost + MIN_FLOAT_NUMBER);
                    }

                  for(int i = 0; i < TIMEBINS; i++)
                    {
                      double sum = domainAssign[n].bin_HydroCost[i] + tasklist[target].bin_HydroCost[i];
                      if(sum < max_HydroCost[i])
                        sum = max_HydroCost[i];

                      runtime += sum / (totsphcost + MIN_FLOAT_NUMBER);
                    }

                  double load = domainAssign[n].load + tasklist[target].load;
                  if(load < max_load)
                    load = max_load;

                  runtime += ((double)load) / totpartcount;

                  if(runtime < best_runtime || best_target < 0)
                    {
                      best_runtime = runtime;
                      best_target  = target;
                    }
                }
            }
        }

      /* now check also the load queue */
      for(np = RB_MIN(mytree, &queue_load), rep = 0; np != NULL && rep < n_lowest; np = RB_NEXT(mytree, &queue_load, np), rep++)
        candidates[rep] = np->target;

      int target;

      for(rep = 0; rep < n_lowest; rep++)
        {
          target = candidates[rep];

          double runtime = 0;

          for(int i = 0; i < TIMEBINS; i++)
            {
              double sum = domainAssign[n].bin_GravCost[i] + tasklist[target].bin_GravCost[i];
              if(sum < max_GravCost[i])
                sum = max_GravCost[i];

              runtime += sum / (totgravcost + 1.0e-60);
            }

          for(int i = 0; i < TIMEBINS; i++)
            {
              double sum = domainAssign[n].bin_HydroCost[i] + tasklist[target].bin_HydroCost[i];
              if(sum < max_HydroCost[i])
                sum = max_HydroCost[i];

              runtime += sum / (totsphcost + 1.0e-60);
            }

          double load = domainAssign[n].load + tasklist[target].load;
          if(load < max_load)
            load = max_load;

          runtime += ((double)load) / totpartcount;

          if(runtime < best_runtime || best_target < 0)
            {
              best_runtime = runtime;
              best_target  = target;
            }
        }

      if(best_target < 0)
        terminate("best_target < 0");

      target = best_target;

      domainAssign[n].task = target;
      tasklist[target].work += domainAssign[n].work;
      tasklist[target].load += domainAssign[n].load;
      tasklist[target].worksph += domainAssign[n].worksph;
      tasklist[target].count++;

      /* now update the elements in the sorted trees */

      RB_REMOVE(mytree, &queue_load, &nload[target]);
      nload[target].pri = tasklist[target].load;
      RB_INSERT(mytree, &queue_load, &nload[target]);

      if(max_load < tasklist[target].load)
        max_load = tasklist[target].load;

      for(int bin = All.LowestOccupiedTimeBin; bin <= All.HighestOccupiedTimeBin; bin++)
        {
          if(domain_to_be_balanced[bin])
            {
              tasklist[target].bin_GravCost[bin] += domainAssign[n].bin_GravCost[bin];
              tasklist[target].bin_HydroCost[bin] += domainAssign[n].bin_HydroCost[bin];

              double eps_grav = 1.0e-9 * (domainAssign[n].load / totpartcount) *
                                totgravcost; /* these will be added in order to break degeneracies in the sort-order in case the
                                                grav/hydro cost in certain cells is zero */
              double eps_hydro = 1.0e-9 * (domainAssign[n].load / totpartcount) * totsphcost;

              RB_REMOVE(mytree, &queue_gravcost[bin], &ngrav[bin][target]);
              ngrav[bin][target].pri = ngrav[bin][target].pri + domainAssign[n].bin_GravCost[bin] + eps_grav;
              RB_INSERT(mytree, &queue_gravcost[bin], &ngrav[bin][target]);

              RB_REMOVE(mytree, &queue_hydrocost[bin], &nhydro[bin][target]);
              nhydro[bin][target].pri = nhydro[bin][target].pri + domainAssign[n].bin_HydroCost[bin] + eps_hydro;
              RB_INSERT(mytree, &queue_hydrocost[bin], &nhydro[bin][target]);

              if(max_GravCost[bin] < tasklist[target].bin_GravCost[bin])
                max_GravCost[bin] = tasklist[target].bin_GravCost[bin];

              if(max_HydroCost[bin] < tasklist[target].bin_HydroCost[bin])
                max_HydroCost[bin] = tasklist[target].bin_HydroCost[bin];
            }
        }
    }

  myfree(candidates);

  /* free the elements for the RB tree again */
  myfree(nload);
  for(int bin = All.HighestOccupiedTimeBin; bin >= All.LowestOccupiedTimeBin; bin--)
    {
      if(domain_to_be_balanced[bin])
        {
          myfree(nhydro[bin]);
          myfree(ngrav[bin]);
        }
    }

  mysort(domainAssign, ndomains, sizeof(struct domain_segments_data), domain_sort_task);

  for(int n = 0; n < ndomains; n++)
    {
      DomainStartList[n] = domainAssign[n].start;
      DomainEndList[n]   = domainAssign[n].end;

      for(int i = DomainStartList[n]; i <= DomainEndList[n]; i++)
        DomainTask[i] = domainAssign[n].task;
    }

  myfree(tasklist);
  myfree(domainAssign);

  double t1 = second();
  mpi_printf("DOMAIN: combining multiple-domains took %g sec\n", timediff(t0, t1));
}

/*! \brief Assign domains to tasks to minimize communication.
 *
 *  This function determines a permutation of the new assignment of domains to
 *  CPUs such that the number of particles that has to be moved given the
 *  current distribution of particles is minimized.
 *
 *  \return void
 */
void domain_optimize_domain_to_task_mapping(void)
{
  double t0 = second();

  int *count_per_task = mymalloc_clear("count_per_task", NTask * sizeof(int));

  /* count how many we want to send to each task */
  for(int i = 0; i < NumPart; i++)
    {
      int no = 0;

      while(topNodes[no].Daughter >= 0)
        no = topNodes[no].Daughter + (Key[i] - topNodes[no].StartKey) / (topNodes[no].Size >> 3);

      no = topNodes[no].Leaf;

      int task = DomainTask[no];
      count_per_task[task]++;
    }

  /* find the task that holds most of our particles (we really would like to be this task) */

  int maxcount = count_per_task[0], maxtask = 0;
  for(int i = 1; i < NTask; i++)
    if(count_per_task[i] > maxcount)
      {
        maxcount = count_per_task[i];
        maxtask  = i;
      }

  struct domain_count_data loc_count;
  struct domain_count_data *domain_count = mymalloc("domain_count", NTask * sizeof(struct domain_count_data));

  loc_count.task       = maxtask;
  loc_count.count      = maxcount;
  loc_count.origintask = ThisTask;

  MPI_Allgather(&loc_count, sizeof(struct domain_count_data), MPI_BYTE, domain_count, sizeof(struct domain_count_data), MPI_BYTE,
                MPI_COMM_WORLD);

  qsort(domain_count, NTask, sizeof(struct domain_count_data), domain_compare_count);

  /* this array will hold a permutation of all tasks constructed such that
     particle exchange should be minimized */

  int *new_task = mymalloc("new_task", NTask * sizeof(int));

  /* this array will now flag tasks that have been assigned */
  for(int i = 0; i < NTask; i++)
    {
      count_per_task[i] = 0;
      new_task[i]       = -1;
    }

  for(int i = 0; i < NTask; i++)
    {
      int task   = domain_count[i].task;
      int origin = domain_count[i].origintask;

      if(new_task[task] == -1 && count_per_task[origin] == 0)
        {
          count_per_task[origin] = 1; /* taken */
          new_task[task]         = origin;
        }
    }

  /* now we have to fill up still unassigned ones in case there were collisions */
  for(int i = 0, j = 0; i < NTask; i++)
    {
      if(new_task[i] == -1)
        {
          while(count_per_task[j])
            j++;

          new_task[i]       = j;
          count_per_task[j] = 1;
        }
    }

  int *copy_DomainStartList = mymalloc("copy_DomainStartList", All.MultipleDomains * NTask * sizeof(int));
  int *copy_DomainEndList   = mymalloc("copy_DomainEndList", All.MultipleDomains * NTask * sizeof(int));

  memcpy(copy_DomainStartList, DomainStartList, All.MultipleDomains * NTask * sizeof(int));
  memcpy(copy_DomainEndList, DomainEndList, All.MultipleDomains * NTask * sizeof(int));

  /* apply permutation to DomainTask assignment */

  for(int i = 0; i < NTask; i++)
    for(int m = 0; m < All.MultipleDomains; m++)
      {
        DomainStartList[new_task[i] * All.MultipleDomains + m] = copy_DomainStartList[i * All.MultipleDomains + m];

        DomainEndList[new_task[i] * All.MultipleDomains + m] = copy_DomainEndList[i * All.MultipleDomains + m];
      }

  myfree(copy_DomainEndList);
  myfree(copy_DomainStartList);

  for(int i = 0; i < NTopleaves; i++)
    DomainTask[i] = new_task[DomainTask[i]];

  myfree(new_task);
  myfree(domain_count);
  myfree(count_per_task);

  double t1 = second();
  mpi_printf("DOMAIN: task reshuffling took %g sec\n", timediff(t0, t1));
}
