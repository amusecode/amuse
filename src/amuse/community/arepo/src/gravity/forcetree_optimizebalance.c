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
 * \file        src/gravity/forcetree_optimizebalance.c
 * \date        05/2018
 * \brief       Does some preparation work for use of red-black ordered binary
 *              tree based on BSD macros.
 * \details     contains functions:
 *                int force_sort_load(const void *a, const void *b)
 *                double force_get_current_balance(double *impact)
 *                void force_get_global_cost_for_leavenodes(int nexport)
 *                static int mydata_cmp(struct mydata *lhs, struct mydata *rhs)
 *                void force_optimize_domain_mapping(void)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 20.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/bsd_tree.h"
#include "../domain/domain.h"

/* \brief Structure of my tree nodes.
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

/*! \brief Data structure that describes force-segment.
 */
static struct force_segments_data
{
  int start, end, task;
  double work, cost, count, normalized_load;
} * force_domainAssign;

/*! \brief Comparison function for force_segments_data.
 *
 *  Sorting kernel.
 *
 *  \param[in] a First object.
 *  \param[in] b Second object.
 *
 *  \return (-1,0,1), -1 if a->normalized_load > b->normalized_load.
 */
int force_sort_load(const void *a, const void *b)
{
  if(((struct force_segments_data *)a)->normalized_load > (((struct force_segments_data *)b)->normalized_load))
    return -1;

  if(((struct force_segments_data *)a)->normalized_load < (((struct force_segments_data *)b)->normalized_load))
    return +1;

  return 0;
}

static double oldmax, oldsum;

/*! \brief Calculates current balance.
 *
 *  \param[out] impact Impact factor of imbalance (1 if optimally balanced).
 *
 *  \return Domain balance = max(cost) / average(cost).
 */
double force_get_current_balance(double *impact)
{
#ifndef NO_MPI_IN_PLACE
  MPI_Allreduce(MPI_IN_PLACE, TaskCost, NTask, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else  /* #ifndef NO_MPI_IN_PLACE */
  double *inTaskCost = mymalloc("inTaskCost", NTask * sizeof(double));
  ;
  memcpy(inTaskCost, TaskCost, NTask * sizeof(double));
  MPI_Allreduce(inTaskCost, TaskCost, NTask, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  myfree(inTaskCost);
#endif /* #ifndef NO_MPI_IN_PLACE #else */

  int i;
  for(i = 0, oldmax = oldsum = 0; i < NTask; i++)
    {
      oldsum += TaskCost[i];
      if(oldmax < TaskCost[i])
        oldmax = TaskCost[i];
    }

  *impact = 1.0 + domain_grav_weight[All.HighestActiveTimeBin] * (oldmax - oldsum / NTask) / All.TotGravCost;

  return oldmax / (oldsum / NTask);
}

/*! \brief Gather cost data of all leaf-nodes and communicate result.
 *
 *  \param[in] nexport Number of exported nodes.
 *
 *  \return void
 */
void force_get_global_cost_for_leavenodes(int nexport)
{
  int i, j, n, nimport, idx, task, ngrp;

  struct node_data
  {
    double domainCost;
    int domainCount;
    int no;
  } * export_node_data, *import_node_data;

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

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  export_node_data = mymalloc("export_node_data", nexport * sizeof(struct node_data));
  import_node_data = mymalloc("import_node_data", nimport * sizeof(struct node_data));

  for(i = 0; i < nexport; i++)
    {
      int task = ListNoData[i].task;
      int ind  = Send_offset[task] + Send_count[task]++;

      export_node_data[ind].domainCost  = ListNoData[i].domainCost;
      export_node_data[ind].domainCount = ListNoData[i].domainCount;
      export_node_data[ind].no          = ListNoData[i].no;
    }

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&export_node_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct node_data), MPI_BYTE, recvTask,
                       TAG_DENS_B, &import_node_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct node_data), MPI_BYTE,
                       recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  for(i = 0; i < nimport; i++)
    {
      int no = import_node_data[i].no;
      DomainCost[no] += import_node_data[i].domainCost;
      DomainCount[no] += import_node_data[i].domainCount;
    }

  myfree(import_node_data);
  myfree(export_node_data);

  /* now share the cost data across all processors */
  struct DomainNODE
  {
    double domainCost;
    int domainCount;
  } * DomainMoment, *loc_DomainMoment;

  DomainMoment = (struct DomainNODE *)mymalloc("DomainMoment", NTopleaves * sizeof(struct DomainNODE));

  /* share the cost data accross CPUs */
  int *recvcounts = (int *)mymalloc("recvcounts", sizeof(int) * NTask);
  int *recvoffset = (int *)mymalloc("recvoffset", sizeof(int) * NTask);
  int *bytecounts = (int *)mymalloc("bytecounts", sizeof(int) * NTask);
  int *byteoffset = (int *)mymalloc("byteoffset", sizeof(int) * NTask);

  for(task = 0; task < NTask; task++)
    recvcounts[task] = 0;

  for(n = 0; n < NTopleaves; n++)
    recvcounts[DomainTask[n]]++;

  for(task = 0; task < NTask; task++)
    bytecounts[task] = recvcounts[task] * sizeof(struct DomainNODE);

  for(task = 1, recvoffset[0] = 0, byteoffset[0] = 0; task < NTask; task++)
    {
      recvoffset[task] = recvoffset[task - 1] + recvcounts[task - 1];
      byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];
    }

  loc_DomainMoment = (struct DomainNODE *)mymalloc("loc_DomainMoment", recvcounts[ThisTask] * sizeof(struct DomainNODE));

  for(n = 0, idx = 0; n < NTopleaves; n++)
    {
      if(DomainTask[n] == ThisTask)
        {
          loc_DomainMoment[idx].domainCost  = DomainCost[n];
          loc_DomainMoment[idx].domainCount = DomainCount[n];
          idx++;
        }
    }

  MPI_Allgatherv(loc_DomainMoment, bytecounts[ThisTask], MPI_BYTE, DomainMoment, bytecounts, byteoffset, MPI_BYTE, MPI_COMM_WORLD);

  for(task = 0; task < NTask; task++)
    recvcounts[task] = 0;

  for(n = 0; n < NTopleaves; n++)
    {
      task = DomainTask[n];
      if(task != ThisTask)
        {
          idx = recvoffset[task] + recvcounts[task]++;

          DomainCost[n]  = DomainMoment[idx].domainCost;
          DomainCount[n] = DomainMoment[idx].domainCount;
        }
    }

  myfree(loc_DomainMoment);
  myfree(byteoffset);
  myfree(bytecounts);
  myfree(recvoffset);
  myfree(recvcounts);
  myfree(DomainMoment);
}

/*! \brief Comparison function of tree elements.
 *
 *  Compares
 *    - pri and if this is equal
 *    - target
 *
 *  \param[in] lhs First mydata object.
 *  \param[in] rhs Second mydata object.
 *
 *  \return (-1,0,1) -1 if lhs < rhs.
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

/*! \brief Optimization algorithm for the workload balance.
 *
 *  \return void
 */
void force_optimize_domain_mapping(void)
{
  int i, j;

  double fac_cost  = 0.5 / oldsum;
  double fac_count = 0.5 / All.TotNumPart;

  int ncpu              = NTask * All.MultipleDomains;
  int ndomain           = NTopleaves;
  double workavg        = 1.0 / ncpu;
  double workhalfnode   = 0.5 / NTopleaves;
  double work_before    = 0;
  double workavg_before = 0;

  int start = 0;

  force_domainAssign = mymalloc("force_domainAssign", ncpu * sizeof(struct force_segments_data));

  for(i = 0; i < ncpu; i++)
    {
      double work = 0, cost = 0, count = 0;
      int end = start;

      cost += fac_cost * DomainCost[end];
      count += fac_count * DomainCount[end];
      work += fac_cost * DomainCost[end] + fac_count * DomainCount[end];

      while((work + work_before + (end + 1 < NTopleaves ? fac_cost * DomainCost[end + 1] + fac_count * DomainCount[end + 1] : 0) <
             workavg + workavg_before + workhalfnode) ||
            (i == ncpu - 1 && end < ndomain - 1))
        {
          if((ndomain - end) > (ncpu - i))
            end++;
          else
            break;

          cost += fac_cost * DomainCost[end];
          count += fac_count * DomainCount[end];
          work += fac_cost * DomainCost[end] + fac_count * DomainCount[end];
        }

      force_domainAssign[i].start = start;
      force_domainAssign[i].end   = end;
      force_domainAssign[i].work  = work;
      force_domainAssign[i].cost  = cost;
      force_domainAssign[i].count = count;

      force_domainAssign[i].normalized_load = cost + count; /* note: they are already multiplied by fac_cost/fac_count */

      work_before += work;
      workavg_before += workavg;
      start = end + 1;
    }

  qsort(force_domainAssign, ncpu, sizeof(struct force_segments_data), force_sort_load);

  /* create three priority trees, one for the cost load, one for the particle count, and one for the combined cost */
  struct mytree queues[3]; /* 0=cost, 1=count, 2=combi */

  struct mydata *ncost  = mymalloc("ncost", NTask * sizeof(struct mydata));
  struct mydata *ncount = mymalloc("ncount", NTask * sizeof(struct mydata));
  struct mydata *ncombi = mymalloc("ncombi", NTask * sizeof(struct mydata));

  RB_INIT(&queues[0]);
  RB_INIT(&queues[1]);
  RB_INIT(&queues[2]);

  /* fill in all the tasks into the trees. The priority will be the current cost/count, the tag 'val' is used to label the task */
  for(i = 0; i < NTask; i++)
    {
      ncost[i].pri    = 0;
      ncost[i].target = i;
      RB_INSERT(mytree, &queues[0], &ncost[i]);

      ncount[i].pri    = 0;
      ncount[i].target = i;
      RB_INSERT(mytree, &queues[1], &ncount[i]);

      ncombi[i].pri    = 0;
      ncombi[i].target = i;
      RB_INSERT(mytree, &queues[2], &ncombi[i]);
    }

  double max_load = 0;
  double max_cost = 0;

  int n_lowest = MAX_FIRST_ELEMENTS_CONSIDERED;
  if(n_lowest > NTask)
    n_lowest = NTask;

  int rep, *candidates = mymalloc("candidates", n_lowest * sizeof(int));
  struct mydata *np;

  for(i = 0; i < ncpu; i++)
    {
      /* pick the least work-loaded target from the queue, and the least particle-loaded, and then decide which choice
         gives the smallest load overall */
      double cost, load;
      double bestwork = 1.0e30;
      int q, target = -1;

      for(q = 0; q < 3; q++)
        {
          /* look up the n_lowest smallest elements from the tree */
          for(np = RB_MIN(mytree, &queues[q]), rep = 0; np != NULL && rep < n_lowest; np = RB_NEXT(mytree, &queues[q], np), rep++)
            candidates[rep] = np->target;

          for(rep = 0; rep < n_lowest; rep++)
            {
              int t = candidates[rep];

              cost = ncost[t].pri + force_domainAssign[i].cost;
              load = ncount[t].pri + force_domainAssign[i].count;
              if(cost < max_cost)
                cost = max_cost;
              if(load < max_load)
                load = max_load;
              double w = cost + load;
              if(w < bestwork)
                {
                  bestwork = w;
                  target   = t;
                }
            }
        }

      force_domainAssign[i].task = target;

      cost = ncost[target].pri + force_domainAssign[i].cost;
      load = ncount[target].pri + force_domainAssign[i].count;

      RB_REMOVE(mytree, &queues[0], &ncost[target]);
      ncost[target].pri = cost;
      RB_INSERT(mytree, &queues[0], &ncost[target]);

      RB_REMOVE(mytree, &queues[1], &ncount[target]);
      ncount[target].pri = load;
      RB_INSERT(mytree, &queues[1], &ncount[target]);

      RB_REMOVE(mytree, &queues[2], &ncombi[target]);
      ncombi[target].pri = cost + load;
      RB_INSERT(mytree, &queues[2], &ncombi[target]);

      if(max_cost < cost)
        max_cost = cost;

      if(max_load < load)
        max_load = load;
    }

  myfree(candidates);

  /* free tree nodes again */
  myfree(ncombi);
  myfree(ncount);
  myfree(ncost);

  for(i = 0; i < ncpu; i++)
    for(j = force_domainAssign[i].start; j <= force_domainAssign[i].end; j++)
      DomainNewTask[j] = force_domainAssign[i].task;

  myfree(force_domainAssign);

  for(i = 0; i < NTask; i++)
    {
      TaskCost[i]  = 0;
      TaskCount[i] = 0;
    }

  for(i = 0; i < NTopleaves; i++)
    {
      TaskCost[DomainNewTask[i]] += DomainCost[i];
      TaskCount[DomainNewTask[i]] += DomainCount[i];
    }

  double max, sum, maxload, sumload;
  for(i = 0, max = sum = 0, maxload = sumload = 0; i < NTask; i++)
    {
      sum += TaskCost[i];
      if(max < TaskCost[i])
        max = TaskCost[i];
      sumload += TaskCount[i];
      if(maxload < TaskCount[i])
        maxload = TaskCount[i];
    }

  mpi_printf("FORCETREE: Active-TimeBin=%d  [unoptimized work-balance=%g]  new work-balance=%g, new load-balance=%g\n",
             All.HighestActiveTimeBin, oldmax / (oldsum / NTask), max / (sum / NTask), maxload / (sumload / NTask));

  if((max / (sum / NTask) > oldmax / (oldsum / NTask)) || (maxload > All.MaxPart))
    {
      mpi_printf(
          "FORCETREE: The work-load is either worse than before or the memory-balance is not viable. We keep the old distribution.\n");
      memcpy(DomainNewTask, DomainTask, NTopleaves * sizeof(int));
    }
}
