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
 * \file        src/subfind/subfind.c
 * \date        05/2018
 * \brief       Main routines of the subfind sub-halo finder.
 * \details     contains functions:
 *                double subfind_get_particle_balance(void)
 *                void subfind(int num)
 *                void subfind_reorder_according_to_submp(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 11.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"
#include "../fof/fof.h"

#ifdef SUBFIND
#include "subfind.h"

/*! \brief Gets a measure of the particle load balance.
 *
 *  \return Maximum number of particle at one core divided by its average.
 */
double subfind_get_particle_balance(void)
{
  int maxpart;
  long long sum;
  MPI_Allreduce(&NumPart, &maxpart, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  sumup_large_ints(1, &NumPart, &sum);
  return maxpart / (((double)sum) / NTask);
}

/*! \brief Main subfind algorithm.
 *
 *  \param[in] num Index of this snapshot output.
 *
 *  \return void
 */
void subfind(int num)
{
  double t0, t1, tstart, tend, cputime;
  int i, gr, nlocid, offset;

  TIMER_START(CPU_SUBFIND);

  tstart = second();

  mpi_printf("\nSUBFIND: We now execute a parallel version of SUBFIND.\n");

  /* let's determine the local dark matter densities */

  TIMER_STOP(CPU_SUBFIND);
  construct_forcetree(0, 0, 1, All.HighestOccupiedTimeBin); /* build forcetree with all particles */
  TIMER_START(CPU_SUBFIND);

  cputime = subfind_density(FIND_SMOOTHING_LENGTHS);
  mpi_printf("SUBFIND: iteration to correct primary neighbor count took %g sec\n", cputime);

  /* free the tree storage again */
  myfree(Father);
  myfree(Nextnode);
  myfree(Tree_Points);
  force_treefree();

  TIMER_STOP(CPU_SUBFIND);
  construct_forcetree(0, 0, 0, All.HighestOccupiedTimeBin); /* build forcetree with all particles */
  TIMER_START(CPU_SUBFIND);

  cputime = subfind_density(FIND_TOTAL_DENSITIES);
  mpi_printf("SUBFIND: density() took %g sec\n", cputime);

  /* free the tree storage again */
  myfree(Father);
  myfree(Nextnode);
  myfree(Tree_Points);
  force_treefree();

  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 0)
      {
#ifdef CELL_CENTER_GRAVITY
        for(int j = 0; j < 3; j++)
          PS[i].Center[j] = SphP[i].Center[j];
#endif /* #ifdef CELL_CENTER_GRAVITY */
        PS[i].Utherm = SphP[i].Utherm;
      }
    else
      PS[i].Utherm = 0;

  SubTreeAllocFactor = All.TreeAllocFactor;

  /* Count, how many groups are above this limit, and how many processors we need for them */
  int ncount = 0, nprocs = 0;
  int seriallen = 0;
  long long sum_seriallen;

  double GroupSize = 0.6;

  do
    {
      ncount    = 0;
      nprocs    = 0;
      seriallen = 0;

      /* Let's set a fiducial size for the maximum group size before we select the collective subfind algorithm */
      MaxSerialGroupLen = (int)(GroupSize * All.TotNumPart / NTask);

      for(i = 0; i < Ngroups; i++)
        if(Group[i].Len > MaxSerialGroupLen)
          {
            ncount++;
            nprocs += ((Group[i].Len - 1) / MaxSerialGroupLen) + 1;
          }
        else
          seriallen += Group[i].Len;

      MPI_Allreduce(&ncount, &Ncollective, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&nprocs, &NprocsCollective, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      sumup_large_ints(1, &seriallen, &sum_seriallen);

      GroupSize += 0.05;
    }
  while(NprocsCollective > 0 && NprocsCollective >= NTask - 1);

  if(GroupSize > 0.65)
    {
      mpi_printf("Increased GroupSize to %g.\n", GroupSize);
    }

  MPI_Allreduce(&ncount, &Ncollective, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nprocs, &NprocsCollective, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  sumup_large_ints(1, &seriallen, &sum_seriallen);

  mpi_printf("SUBFIND: Number of FOF halos treated with collective SubFind code = %d\n", Ncollective);
  mpi_printf("SUBFIND: Number of processors used in different partitions for the collective SubFind code = %d\n", NprocsCollective);
  mpi_printf("SUBFIND: (The adopted size-limit for the collective algorithm was %d particles.)\n", MaxSerialGroupLen);
  mpi_printf("SUBFIND: The other %d FOF halos are treated in parallel with serial code\n", TotNgroups - Ncollective);

  /* set up a global table that informs about the processor assignment of the groups that are treated collectively */
  ProcAssign                             = mymalloc_movable(&ProcAssign, "ProcAssign", Ncollective * sizeof(struct proc_assign_data));
  struct proc_assign_data *locProcAssign = mymalloc("locProcAssign", ncount * sizeof(struct proc_assign_data));

  for(i = 0, ncount = 0; i < Ngroups; i++)
    if(Group[i].Len > MaxSerialGroupLen)
      {
        locProcAssign[ncount].GrNr = Group[i].GrNr;
        locProcAssign[ncount].Len  = Group[i].Len;
        ncount++;
      }

  /* gather the information on the collective groups accross all CPUs */
  int *recvcounts = (int *)mymalloc("recvcounts", sizeof(int) * NTask);
  int *bytecounts = (int *)mymalloc("bytecounts", sizeof(int) * NTask);
  int *byteoffset = (int *)mymalloc("byteoffset", sizeof(int) * NTask);

  MPI_Allgather(&ncount, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);

  int task;
  for(task = 0; task < NTask; task++)
    bytecounts[task] = recvcounts[task] * sizeof(struct proc_assign_data);

  for(task = 1, byteoffset[0] = 0; task < NTask; task++)
    byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];

  MPI_Allgatherv(locProcAssign, bytecounts[ThisTask], MPI_BYTE, ProcAssign, bytecounts, byteoffset, MPI_BYTE, MPI_COMM_WORLD);

  myfree(byteoffset);
  myfree(bytecounts);
  myfree(recvcounts);
  myfree(locProcAssign);

  /* make sure, the table is sorted in ascending group-number order */
  qsort(ProcAssign, Ncollective, sizeof(struct proc_assign_data), subfind_compare_procassign_GrNr);

  /* assign the processor sets for the collective groups and set disjoint color-flag to later split the processors into different
   * communicators */
  for(i = 0, nprocs = 0, CommSplitColor = Ncollective; i < Ncollective; i++)
    {
      ProcAssign[i].FirstTask = nprocs;
      ProcAssign[i].NTask     = ((ProcAssign[i].Len - 1) / MaxSerialGroupLen) + 1;
      nprocs += ProcAssign[i].NTask;

      if(ThisTask >= ProcAssign[i].FirstTask && ThisTask < (ProcAssign[i].FirstTask + ProcAssign[i].NTask))
        CommSplitColor = i;
    }

  /* Now assign a target task for the group. For collective groups, the target task is the master in the CPU set, whereas
   * the serial ones are distributed in a round-robin fashion to the remaining CPUs
   */
  for(i = 0; i < Ngroups; i++)
    {
      if(Group[i].Len > MaxSerialGroupLen) /* we have a collective group */
        {
          if(Group[i].GrNr >= Ncollective || Group[i].GrNr < 0)
            terminate("odd");
          Group[i].TargetTask = ProcAssign[Group[i].GrNr].FirstTask;
        }
      else
        Group[i].TargetTask = ((Group[i].GrNr - Ncollective) % (NTask - NprocsCollective)) + NprocsCollective;
    }

  /* distribute the groups */
  subfind_distribute_groups();
  qsort(Group, Ngroups, sizeof(struct group_properties), fof_compare_Group_GrNr);

  /* assign target CPUs for the particles in groups */
  /* the particles not in groups will be distributed such that a uniform particle load results */
  t0                  = second();
  int *count_loc_task = mymalloc_clear("count_loc_task", NTask * sizeof(int));
  int *count_task     = mymalloc("count_task", NTask * sizeof(int));
  int *count_free     = mymalloc("count_free", NTask * sizeof(int));
  int count_loc_free  = 0;

  for(i = 0; i < NumPart; i++)
    {
      if(PS[i].GrNr < TotNgroups) /* particle is in a group */
        {
          if(PS[i].GrNr < Ncollective) /* we are in a collective group */
            PS[i].TargetTask = ProcAssign[PS[i].GrNr].FirstTask + (i % ProcAssign[PS[i].GrNr].NTask);
          else
            PS[i].TargetTask = ((PS[i].GrNr - Ncollective) % (NTask - NprocsCollective)) + NprocsCollective;

          count_loc_task[PS[i].TargetTask]++;
        }
      else
        count_loc_free++;

      PS[i].TargetIndex = 0; /* unimportant here */
    }

  MPI_Allgather(&count_loc_free, 1, MPI_INT, count_free, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allreduce(count_loc_task, count_task, NTask, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  long long sum = 0;
  for(i = 0; i < NTask; i++)
    sum += count_task[i] + count_free[i];

  int maxload = (sum + NTask - 1) / NTask;
  for(i = 0; i < NTask; i++)
    {
      count_task[i] = maxload - count_task[i]; /* this is the amount that can fit on this task */
      if(count_task[i] < 0)
        count_task[i] = 0;
    }

  int current_task = 0;

  for(i = 0; i < ThisTask; i++)
    {
      while(count_free[i] > 0 && current_task < NTask)
        {
          if(count_free[i] < count_task[current_task])
            {
              count_task[current_task] -= count_free[i];
              count_free[i] = 0;
            }
          else
            {
              count_free[i] -= count_task[current_task];
              count_task[current_task] = 0;
              current_task++;
            }
        }
    }

  for(i = 0; i < NumPart; i++)
    {
      if(PS[i].GrNr >=
         TotNgroups) /* particle not in a group. Can in principle stay but we move it such that a good load balance is obtained. */
        {
          while(count_task[current_task] == 0 && current_task < NTask - 1)
            current_task++;

          PS[i].TargetTask = current_task; /* particle not in any group, move it here so that uniform load is achieved */
          count_task[current_task]--;
        }
    }

  myfree(count_free);
  myfree(count_task);
  myfree(count_loc_task);

#ifdef SUBFIND_EXTENDED_PROPERTIES
  int ngroups_cat = 42;     // dummy. not used for any calculation but fct needs to receive a value and we want to keep fct universal.
#endif                      /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
  int nsubgroups_cat = 42;  // dummy. not used for any calculation but fct needs to receive a value and we want to keep fct universal.

  double balance = subfind_get_particle_balance();
  mpi_printf("SUBFIND: particle balance=%g\n", balance);

  /* distribute particles such that groups are completely on the CPU(s) that do the corresponding group(s) */
  fof_subfind_exchange(MPI_COMM_WORLD);
  t1 = second();
  mpi_printf("SUBFIND: subfind_exchange() took %g sec\n", timediff(t0, t1));

  balance = subfind_get_particle_balance();
  mpi_printf("SUBFIND: particle balance for processing=%g\n", balance);

  /* lets estimate the maximum number of substructures we need to store on the local CPU */
  if(ThisTask < NprocsCollective)
    {
      MaxNsubgroups = (ProcAssign[CommSplitColor].Len / ProcAssign[CommSplitColor].NTask) / All.DesLinkNgb;
    }
  else
    {
      for(i = 0, nlocid = 0; i < Ngroups; i++)
        nlocid += Group[i].Len;

      MaxNsubgroups = nlocid / All.DesLinkNgb; /* should be a quite conservative upper limit */
    }

  Nsubgroups = 0;
  SubGroup = (struct subgroup_properties *)mymalloc_movable(&SubGroup, "SubGroup", MaxNsubgroups * sizeof(struct subgroup_properties));

  /* we can now split the communicator to give each collectively treated group its own processor set */
  MPI_Comm_split(MPI_COMM_WORLD, CommSplitColor, ThisTask, &SubComm);
  MPI_Comm_size(SubComm, &SubNTask);
  MPI_Comm_rank(SubComm, &SubThisTask);
  SubTagOffset = TagOffset;

  /* here the execution paths for collective groups and serial groups branch. The collective CPUs work in small sets that each
   * deal with one large group. The serial CPUs each deal with several halos by themselves
   */
  if(CommSplitColor < Ncollective) /* we are one of the CPUs that does a collective group */
    {
      /* we now apply a collective version of subfind to the group split across the processors belonging to communicator SubComm
       * The relevant group is the one stored in Group[0] on SubThisTask==0.
       */
      subfind_process_group_collectively(nsubgroups_cat);
    }
  else
    {
      /* now let us sort according to GrNr and Density. This step will temporarily break the association with SphP[] and other arrays!
       */
      submp = (struct submp_data *)mymalloc("submp", sizeof(struct submp_data) * NumPart);
      for(i = 0; i < NumPart; i++)
        {
          PS[i].SubNr         = TotNgroups + 1; /* set a default that is larger than reasonable group number */
          PS[i].OldIndex      = i;
          submp[i].index      = i;
          submp[i].GrNr       = PS[i].GrNr;
          submp[i].DM_Density = PS[i].Density;
        }
      qsort(submp, NumPart, sizeof(struct submp_data), subfind_compare_submp_GrNr_DM_Density);
      subfind_reorder_according_to_submp();
      myfree(submp);

      /* now we have the particles in each group consecutively */
      if(SubThisTask == 0)
        printf(
            "SUBFIND-SERIAL: Start to do %d small groups (cumulative length %lld) with serial subfind algorithm on %d processors "
            "(root-node=%d)\n",
            TotNgroups - Ncollective, sum_seriallen, SubNTask, ThisTask);

      /* we now apply a serial version of subfind to the local groups */
      t0 = second();
      for(gr = 0, offset = 0; gr < Ngroups; gr++)
        {
          if(((Group[gr].GrNr - Ncollective) % (NTask - NprocsCollective)) + NprocsCollective == ThisTask)
            offset = subfind_process_group_serial(gr, offset, nsubgroups_cat);
          else
            terminate("how come that we have this group number?");
        }

      MPI_Barrier(SubComm);
      t1 = second();
      if(SubThisTask == 0)
        printf("SUBFIND-SERIAL: processing of serial groups took %g sec\n", timediff(t0, t1));

      /* undo local rearrangement that made groups consecutive. After that, the association of SphP[] will be correct again */
      submp = (struct submp_data *)mymalloc("submp", sizeof(struct submp_data) * NumPart);
      for(i = 0; i < NumPart; i++)
        {
          submp[i].index    = i;
          submp[i].OldIndex = PS[i].OldIndex;
        }
      qsort(submp, NumPart, sizeof(struct submp_data), subfind_compare_submp_OldIndex);
      subfind_reorder_according_to_submp();
      myfree(submp);
    }

  /* free the communicator */
  MPI_Comm_free(&SubComm);

  /* make common allocation on all tasks */
  int max_load, max_loadsph, load;

  /* for resize */
  load = All.MaxPart;
  MPI_Allreduce(&load, &max_load, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  load = All.MaxPartSph;
  MPI_Allreduce(&load, &max_loadsph, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  /* do resize */
  All.MaxPart = max_load;
  reallocate_memory_maxpart();
  PS = (struct subfind_data *)myrealloc_movable(PS, All.MaxPart * sizeof(struct subfind_data));

  All.MaxPartSph = max_loadsph;
  reallocate_memory_maxpartsph();

  /* distribute particles back to original CPU */
  t0 = second();
  for(i = 0; i < NumPart; i++)
    {
      PS[i].TargetTask  = PS[i].OriginTask;
      PS[i].TargetIndex = PS[i].OriginIndex;
    }

  fof_subfind_exchange(MPI_COMM_WORLD);
  t1 = second();
  if(ThisTask == 0)
    printf("SUBFIND: subfind_exchange() (for return to original CPU)  took %g sec\n", timediff(t0, t1));

  TIMER_STOP(CPU_SUBFIND);
  construct_forcetree(0, 0, 0, All.HighestOccupiedTimeBin); /* build forcetree with all particles */
  TIMER_START(CPU_SUBFIND);

  /* compute spherical overdensities for FOF groups */
  cputime = subfind_overdensity();
  mpi_printf("SUBFIND: determining spherical overdensity masses took %g sec\n", cputime);

  myfree(Father);
  myfree(Nextnode);
  myfree(Tree_Points);
  force_treefree();

#ifdef SUBFIND_EXTENDED_PROPERTIES
  subfind_add_grp_props_calc_fof_angular_momentum(num, ngroups_cat);
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */

  MPI_Allreduce(&Nsubgroups, &TotNsubgroups, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  /* sort the groups according to group/subgroup-number */
  t0 = second();
  parallel_sort(Group, Ngroups, sizeof(struct group_properties), fof_compare_Group_GrNr);
  parallel_sort(SubGroup, Nsubgroups, sizeof(struct subgroup_properties), subfind_compare_SubGroup_GrNr_SubNr);
  t1 = second();
  mpi_printf("SUBFIND: assembled and ordered groups and subgroups (took %g sec)\n", timediff(t0, t1));

  /* determine largest subgroup and total particle/cell count in substructures */
  int lenmax, glob_lenmax, totlen;
  long long totsublength;
  for(i = 0, totlen = 0, lenmax = 0; i < Nsubgroups; i++)
    {
      totlen += SubGroup[i].Len;

      if(SubGroup[i].Len > lenmax)
        lenmax = SubGroup[i].Len;
    }
  sumup_large_ints(1, &totlen, &totsublength);
  MPI_Reduce(&lenmax, &glob_lenmax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

  /* set binding energy of fuzz to zero, was overwritten with Hsml before; needed for proper snapshot sorting of fuzz */
  for(i = 0; i < NumPart; i++)
    if(PS[i].SubNr == TotNgroups + 1)
      PS[i].BindingEnergy = 0;

  TIMER_STOP(CPU_SUBFIND);
  TIMER_START(CPU_SNAPSHOT);

  /* now final output of catalogue */
  subfind_save_final(num);

  TIMER_STOP(CPU_SNAPSHOT);
  TIMER_START(CPU_SUBFIND);

  tend = second();

  if(ThisTask == 0)
    {
      printf("SUBFIND: Finished with SUBFIND.  (total time=%g sec)\n", timediff(tstart, tend));
      printf("SUBFIND: Total number of subhalos with at least %d particles: %d\n", All.DesLinkNgb, TotNsubgroups);
      if(TotNsubgroups > 0)
        {
          printf("SUBFIND: Largest subhalo has %d particles/cells.\n", glob_lenmax);
          printf("SUBFIND: Total number of particles/cells in subhalos: %lld\n", totsublength);
        }
    }

  myfree_movable(SubGroup);
  myfree_movable(ProcAssign);

  TIMER_STOP(CPU_SUBFIND);
}

/*! \brief Reorders particles in P and SphP array.
 *
 *  Reordering given by the submp array.
 *
 *  \return void
 */
void subfind_reorder_according_to_submp(void)
{
  int i;
  struct particle_data Psave, Psource;
  struct subfind_data PSsave, PSsource;
  int idsource, idsave, dest;
  int *Id;

  Id = (int *)mymalloc("Id", sizeof(int) * (NumPart));

  for(i = 0; i < NumPart; i++)
    Id[submp[i].index] = i;

  for(i = 0; i < NumPart; i++)
    {
      if(Id[i] != i)
        {
          Psource  = P[i];
          PSsource = PS[i];
          idsource = Id[i];

          dest = Id[i];

          do
            {
              Psave  = P[dest];
              PSsave = PS[dest];
              idsave = Id[dest];

              P[dest]  = Psource;
              PS[dest] = PSsource;
              Id[dest] = idsource;

              if(dest == i)
                break;

              Psource  = Psave;
              PSsource = PSsave;
              idsource = idsave;

              dest = idsource;
            }
          while(1);
        }
    }

  myfree(Id);
}

#endif /* #ifdef SUBFIND */
