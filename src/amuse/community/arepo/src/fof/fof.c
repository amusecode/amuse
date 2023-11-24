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
 * \file        src/fof/fof.c
 * \date        05/2018
 * \brief       Parallel friend of friends (FoF) group finder.
 * \details     contains functions:
 *                void fof_fof(int num)
 *                void fof_prepare_output_order(void)
 *                double fof_get_comoving_linking_length(void)
 *                void fof_compile_catalogue(void)
 *                void fof_assign_group_numbers(void)
 *                void fof_compute_group_properties(int gr, int start, int len)
 *                void fof_exchange_group_data(void)
 *                void fof_finish_group_properties(void)
 *                double fof_periodic(double x)
 *                double fof_periodic_wrap(double x)
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

static MyIDType *MinID;
static int *Head, *Len, *Next, *Tail, *MinIDTask;

/*! \brief Main routine to execute the friend of friends group finder.
 *
 *  If called with num == -1 as argument, only FOF is carried out and no group
 *  catalogs are saved to disk. If num >= 0, the code will store the
 *  group/subgroup catalogs, and bring the particles into output order.
 *  In this case, the calling routine (which is normally savepositions()) will
 *  need to free PS[] and bring the particles back into the original order,
 *  as well as reestablished the mesh.
 *
 *  \param[in] num Index of output; if negative, no output written.
 *
 *  \return void
 */
void fof_fof(int num)
{
  int i, start, lenloc, largestgroup;
  double t0, t1, cputime;

  TIMER_START(CPU_FOF);

  mpi_printf("FOF: Begin to compute FoF group catalogue...  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  if(num >= 0 && RestartFlag != 3 && RestartFlag != 6)
    {
      /* let's discard an existing mesh - we do this here to reduce the peak memory usage, even at the price of
       * having to recreate it later */
      free_mesh();
    }

  if(RestartFlag != 6)
    {
      ngb_treefree();

      domain_free();
    }

  domain_Decomposition();

  ngb_treeallocate();
  ngb_treebuild(NumGas);

  /* check */
  for(i = 0; i < NumPart; i++)
    if((P[i].Mass == 0 && P[i].ID == 0) || (P[i].Type == 4 && P[i].Mass == 0))
      terminate("this should not happen");

  /* this structure will hold auxiliary information for each particle, needed only during group finding */
  PS = (struct subfind_data *)mymalloc_movable(&PS, "PS", All.MaxPart * sizeof(struct subfind_data));

  memset(PS, 0, NumPart * sizeof(struct subfind_data));

  /* First, we save the original location of the particles, in order to be able to revert to this layout later on */
  for(i = 0; i < NumPart; i++)
    {
      PS[i].OriginTask  = ThisTask;
      PS[i].OriginIndex = i;
    }

  fof_OldMaxPart    = All.MaxPart;
  fof_OldMaxPartSph = All.MaxPartSph;

  LinkL = fof_get_comoving_linking_length();

  mpi_printf("FOF: Comoving linking length: %g    (presently allocated=%g MB)\n", LinkL, AllocatedBytes / (1024.0 * 1024.0));

  MinID     = (MyIDType *)mymalloc("MinID", NumPart * sizeof(MyIDType));
  MinIDTask = (int *)mymalloc("MinIDTask", NumPart * sizeof(int));

  Head = (int *)mymalloc("Head", NumPart * sizeof(int));
  Len  = (int *)mymalloc("Len", NumPart * sizeof(int));
  Next = (int *)mymalloc("Next", NumPart * sizeof(int));
  Tail = (int *)mymalloc("Tail", NumPart * sizeof(int));

#ifdef HIERARCHICAL_GRAVITY
  timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, All.HighestOccupiedTimeBin);
#endif /* #ifdef HIERARCHICAL_GRAVITY */

  construct_forcetree(0, 0, 1, All.HighestOccupiedTimeBin); /* build tree for all particles */

#if defined(SUBFIND)
  subfind_density_hsml_guess();
#endif /* #if defined(SUBFIND) */

  /* initialize link-lists */
  for(i = 0; i < NumPart; i++)
    {
      Head[i] = Tail[i] = i;
      Len[i]            = 1;
      Next[i]           = -1;
      MinID[i]          = P[i].ID;
      MinIDTask[i]      = ThisTask;
    }

  /* call routine to find primary groups */
  cputime = fof_find_groups(MinID, Head, Len, Next, Tail, MinIDTask);
  mpi_printf("FOF: group finding took = %g sec\n", cputime);

#ifdef FOF_SECONDARY_LINK_TARGET_TYPES
  myfree(Father);
  myfree(Nextnode);
  myfree(Tree_Points);

  /* now rebuild the tree with all the types selected as secondary link targets */
  construct_forcetree(0, 0, 2, All.HighestOccupiedTimeBin);
#endif /* #ifdef FOF_SECONDARY_LINK_TARGET_TYPES */

#ifdef HIERARCHICAL_GRAVITY
  timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, All.HighestActiveTimeBin);
#endif /* #ifdef HIERARCHICAL_GRAVITY */

  /* call routine to attach secondary particles/cells to primary groups */
  cputime = fof_find_nearest_dmparticle(MinID, Head, Len, Next, Tail, MinIDTask);

  mpi_printf("FOF: attaching gas and star particles to nearest dm particles took = %g sec\n", cputime);

  myfree(Father);
  myfree(Nextnode);
  myfree(Tree_Points);
  force_treefree();

  myfree(Tail);
  myfree(Next);
  myfree(Len);

  t0 = second();

  FOF_PList = (struct fof_particle_list *)mymalloc_movable(&FOF_PList, "FOF_PList", NumPart * sizeof(struct fof_particle_list));

  for(i = 0; i < NumPart; i++)
    {
      FOF_PList[i].MinID     = MinID[Head[i]];
      FOF_PList[i].MinIDTask = MinIDTask[Head[i]];
      FOF_PList[i].Pindex    = i;
    }

  myfree_movable(Head);
  myfree_movable(MinIDTask);
  myfree_movable(MinID);

  FOF_GList = (struct fof_group_list *)mymalloc_movable(&FOF_GList, "FOF_GList", sizeof(struct fof_group_list) * NumPart);

  fof_compile_catalogue();

  t1 = second();
  mpi_printf("FOF: compiling local group data and catalogue took = %g sec\n", timediff(t0, t1));

  MPI_Allreduce(&Ngroups, &TotNgroups, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  sumup_large_ints(1, &Nids, &TotNids);

  if(TotNgroups > 0)
    {
      int largestloc = 0;

      for(i = 0; i < NgroupsExt; i++)
        if(FOF_GList[i].LocCount + FOF_GList[i].ExtCount > largestloc)
          largestloc = FOF_GList[i].LocCount + FOF_GList[i].ExtCount;
      MPI_Allreduce(&largestloc, &largestgroup, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
  else
    largestgroup = 0;

  mpi_printf("FOF: Total number of FOF groups with at least %d particles: %d\n", FOF_GROUP_MIN_LEN, TotNgroups);
  mpi_printf("FOF: Largest FOF group has %d particles.\n", largestgroup);
  mpi_printf("FOF: Total number of particles in FOF groups: %lld\n", TotNids);

  t0 = second();

  MaxNgroups = 2 * imax(NgroupsExt, TotNgroups / NTask + 1);

  Group = (struct group_properties *)mymalloc_movable(&Group, "Group", sizeof(struct group_properties) * MaxNgroups);

  mpi_printf("FOF: group properties are now allocated.. (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  for(i = 0, start = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID < FOF_GList[i].MinID)
        {
          start++;
          if(start > NumPart)
            terminate("start > NumPart");
        }

      if(FOF_PList[start].MinID != FOF_GList[i].MinID)
        terminate("ID mismatch");

      for(lenloc = 0; start + lenloc < NumPart;)
        if(FOF_PList[start + lenloc].MinID == FOF_GList[i].MinID)
          lenloc++;
        else
          break;

      Group[i].MinID     = FOF_GList[i].MinID;
      Group[i].MinIDTask = FOF_GList[i].MinIDTask;

      fof_compute_group_properties(i, start, lenloc);

      start += lenloc;
    }

  fof_exchange_group_data();

  fof_finish_group_properties();

  t1 = second();
  mpi_printf("FOF: computation of group properties took = %g sec\n", timediff(t0, t1));

  fof_assign_group_numbers();

  mpi_printf("FOF: Finished computing FoF groups.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  myfree_movable(FOF_GList);
  myfree_movable(FOF_PList);

#ifdef SUBFIND
  if(num >= 0)
    {
      TIMER_STOP(CPU_FOF);

      subfind(num);

      TIMER_START(CPU_FOF);
    }
#else  /* #ifdef SUBFIND */
  Nsubgroups    = 0;
  TotNsubgroups = 0;
  if(num >= 0)
    {
      TIMER_STOP(CPU_FOF);
      TIMER_START(CPU_SNAPSHOT);

      fof_save_groups(num);

      TIMER_STOP(CPU_SNAPSHOT);
      TIMER_START(CPU_FOF);
    }
#endif /* #ifdef SUBFIND #else */

  myfree_movable(Group);

  mpi_printf("FOF: All FOF related work finished.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

#ifndef FOF_STOREIDS
  if(num >= 0)
    {
      TIMER_STOP(CPU_FOF);
      TIMER_START(CPU_SNAPSHOT);

      /* now distribute the particles into output order */
      t0 = second();
      fof_prepare_output_order();
      fof_subfind_exchange(
          MPI_COMM_WORLD); /* distribute particles such that FOF groups will appear in coherent way in snapshot files */
      t1 = second();
      mpi_printf("FOF: preparing output order of particles took %g sec\n", timediff(t0, t1));

      TIMER_STOP(CPU_SNAPSHOT);
      TIMER_START(CPU_FOF);
    }
  else
    myfree(PS);
#else  /* #ifndef FOF_STOREIDS */
  myfree(PS);
#endif /* #ifndef FOF_STOREIDS #else */

  TIMER_STOP(CPU_FOF);
}

/*! \brief Sorts groups by the desired output order.
 *
 *  \return void
 */
void fof_prepare_output_order(void)
{
  int i, off, ntype[NTYPES];

  struct data_aux_sort *aux_sort = (struct data_aux_sort *)mymalloc("aux_sort", sizeof(struct data_aux_sort) * NumPart);

  for(i = 0; i < NTYPES; i++)
    ntype[i] = 0;

  for(i = 0; i < NumPart; i++)
    {
      aux_sort[i].OriginTask  = ThisTask;
      aux_sort[i].OriginIndex = i;
      aux_sort[i].GrNr        = PS[i].GrNr;
#ifdef SUBFIND
      aux_sort[i].SubNr            = PS[i].SubNr;
      aux_sort[i].DM_BindingEnergy = PS[i].BindingEnergy;
#endif /* #ifdef SUBFIND */
      aux_sort[i].Type = P[i].Type;
      aux_sort[i].ID   = P[i].ID;
#if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT)
      aux_sort[i].FileOrder = P[i].FileOrder;
#endif /* #if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) */

      ntype[P[i].Type]++;
    }

  qsort(aux_sort, NumPart, sizeof(struct data_aux_sort), fof_compare_aux_sort_Type);

  if(RestartFlag == 18)
    {
#if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT)
      for(i = 0, off = 0; i < NTYPES; off += ntype[i], i++)
        parallel_sort(aux_sort + off, ntype[i], sizeof(struct data_aux_sort), fof_compare_aux_sort_FileOrder);
#endif /* #if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) */
    }
  else
    {
      for(i = 0, off = 0; i < NTYPES; off += ntype[i], i++)
        parallel_sort(aux_sort + off, ntype[i], sizeof(struct data_aux_sort), fof_compare_aux_sort_GrNr);
    }

  for(i = 0; i < NumPart; i++)
    {
      aux_sort[i].TargetTask  = ThisTask;
      aux_sort[i].TargetIndex = i;
    }

  /* now bring back into starting order */
  parallel_sort(aux_sort, NumPart, sizeof(struct data_aux_sort), fof_compare_aux_sort_OriginTask_OriginIndex);

  for(i = 0; i < NumPart; i++)
    {
      PS[i].TargetTask  = aux_sort[i].TargetTask;
      PS[i].TargetIndex = aux_sort[i].TargetIndex;
    }

  myfree(aux_sort);
}

/*! \brief Calculate linking length based on mean particle separation.
 *
 *  \return Linking length.
 */
double fof_get_comoving_linking_length(void)
{
  int i, ndm;
  long long ndmtot;
  double mass, masstot, rhodm;

  for(i = 0, ndm = 0, mass = 0; i < NumPart; i++)
    if(((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES)))
      {
        ndm++;
        mass += P[i].Mass;
      }
  sumup_large_ints(1, &ndm, &ndmtot);
  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  rhodm = (All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  return FOF_LINKLENGTH * pow(masstot / ndmtot / rhodm, 1.0 / 3);
}

/*! \brief Compiles the group catalogue.
 *
 *  Combines results from all tasks.
 *
 *  \return void
 */
void fof_compile_catalogue(void)
{
  int i, j, start, nimport, ngrp, recvTask;
  struct fof_group_list *get_FOF_GList;

  /* sort according to MinID */
  mysort(FOF_PList, NumPart, sizeof(struct fof_particle_list), fof_compare_FOF_PList_MinID);

  for(i = 0; i < NumPart; i++)
    {
      FOF_GList[i].MinID     = FOF_PList[i].MinID;
      FOF_GList[i].MinIDTask = FOF_PList[i].MinIDTask;
      if(FOF_GList[i].MinIDTask == ThisTask)
        {
          FOF_GList[i].LocCount = 1;
          FOF_GList[i].ExtCount = 0;
        }
      else
        {
          FOF_GList[i].LocCount = 0;
          FOF_GList[i].ExtCount = 1;
        }
    }

  /* eliminate duplicates in FOF_GList with respect to MinID */

  if(NumPart)
    NgroupsExt = 1;
  else
    NgroupsExt = 0;

  for(i = 1, start = 0; i < NumPart; i++)
    {
      if(FOF_GList[i].MinID == FOF_GList[start].MinID)
        {
          FOF_GList[start].LocCount += FOF_GList[i].LocCount;
          FOF_GList[start].ExtCount += FOF_GList[i].ExtCount;
        }
      else
        {
          start            = NgroupsExt;
          FOF_GList[start] = FOF_GList[i];
          NgroupsExt++;
        }
    }

  /* sort the remaining ones according to task */
  mysort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinIDTask);

  /* count how many we have of each task */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;
  for(i = 0; i < NgroupsExt; i++)
    Send_count[FOF_GList[i].MinIDTask]++;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      if(j == ThisTask) /* we will not exchange the ones that are local */
        Recv_count[j] = 0;
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  get_FOF_GList = (struct fof_group_list *)mymalloc("get_FOF_GList", nimport * sizeof(struct fof_group_list));

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the group info */
              MPI_Sendrecv(&FOF_GList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct fof_group_list), MPI_BYTE, recvTask,
                           TAG_DENS_A, &get_FOF_GList[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct fof_group_list),
                           MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  for(i = 0; i < nimport; i++)
    get_FOF_GList[i].MinIDTask = i;

  /* sort the groups according to MinID */
  mysort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinID);
  mysort(get_FOF_GList, nimport, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinID);

  /* merge the imported ones with the local ones */
  for(i = 0, start = 0; i < nimport; i++)
    {
      while(FOF_GList[start].MinID < get_FOF_GList[i].MinID)
        {
          start++;
          if(start >= NgroupsExt)
            terminate("start >= NgroupsExt");
        }

      if(get_FOF_GList[i].LocCount != 0)
        terminate("start >= NgroupsExt");

      if(FOF_GList[start].MinIDTask != ThisTask)
        terminate("FOF_GList[start].MinIDTask != ThisTask");

      if(FOF_GList[start].MinID != get_FOF_GList[i].MinID)
        terminate(
            "FOF_GList[start].MinID != get_FOF_GList[i].MinID start=%d i=%d FOF_GList[start].MinID=%llu get_FOF_GList[i].MinID=%llu\n",
            start, i, (long long)FOF_GList[start].MinID, (long long)get_FOF_GList[i].MinID);

      FOF_GList[start].ExtCount += get_FOF_GList[i].ExtCount;
    }

  /* copy the size information back into the list, to inform the others */
  for(i = 0, start = 0; i < nimport; i++)
    {
      while(FOF_GList[start].MinID < get_FOF_GList[i].MinID)
        {
          start++;
          if(start >= NgroupsExt)
            terminate("start >= NgroupsExt");
        }

      get_FOF_GList[i].ExtCount = FOF_GList[start].ExtCount;
      get_FOF_GList[i].LocCount = FOF_GList[start].LocCount;
    }

  /* sort the imported/exported list according to MinIDTask */
  mysort(get_FOF_GList, nimport, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinIDTask);
  mysort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinIDTask);

  for(i = 0; i < nimport; i++)
    get_FOF_GList[i].MinIDTask = ThisTask;

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the group info */
              MPI_Sendrecv(&get_FOF_GList[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct fof_group_list), MPI_BYTE,
                           recvTask, TAG_DENS_A, &FOF_GList[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct fof_group_list), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(get_FOF_GList);

  /* eliminate all groups that are too small, and count local groups */
  for(i = 0, Ngroups = 0, Nids = 0; i < NgroupsExt; i++)
    {
      if(FOF_GList[i].LocCount + FOF_GList[i].ExtCount < FOF_GROUP_MIN_LEN)
        {
          FOF_GList[i] = FOF_GList[NgroupsExt - 1];
          NgroupsExt--;
          i--;
        }
      else
        {
          if(FOF_GList[i].MinIDTask == ThisTask)
            {
              Ngroups++;
              Nids += FOF_GList[i].LocCount + FOF_GList[i].ExtCount;
            }
        }
    }

  /* sort the group list according to MinID */
  mysort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinID);
}

/*! \brief Assigns each group a global group number.
 *
 *  \return void
 */
void fof_assign_group_numbers(void)
{
  int i, j, ngr, start, lenloc;
  long long totNids;
  double t0, t1;

  mpi_printf("FOF: start assigning group numbers\n");

  t0 = second();

  /* assign group numbers (at this point, both Group and FOF_GList are sorted by MinID) */
  for(i = 0; i < NgroupsExt; i++)
    {
      FOF_GList[i].LocCount += FOF_GList[i].ExtCount; /* total length */
      FOF_GList[i].ExtCount = ThisTask;               /* original task */
    }

  parallel_sort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_LocCountTaskDiffMinID);

  for(i = 0, ngr = 0; i < NgroupsExt; i++)
    {
      if(FOF_GList[i].ExtCount == FOF_GList[i].MinIDTask)
        ngr++;

      FOF_GList[i].GrNr = ngr - 1;
    }

  MPI_Allgather(&ngr, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);

  /* count how many groups there are on earlier CPUs */
  long long ngr_sum;
  for(j = 0, ngr_sum = 0; j < ThisTask; j++)
    ngr_sum += Send_count[j];

  for(i = 0; i < NgroupsExt; i++)
    FOF_GList[i].GrNr += ngr_sum;

  sumup_large_ints(1, &ngr, &ngr_sum);
  if(ngr_sum != TotNgroups)
    {
      printf("ngr_sum=%d\n", (int)ngr_sum);
      terminate("inconsistency");
    }

  /* bring the group list back into the original order */
  parallel_sort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_ExtCountMinID);

  /* Assign the group numbers to the group properties array */
  for(i = 0, start = 0; i < Ngroups; i++)
    {
      while(FOF_GList[start].MinID < Group[i].MinID)
        {
          start++;
          if(start >= NgroupsExt)
            terminate("start >= NgroupsExt");
        }
      Group[i].GrNr = FOF_GList[start].GrNr;
    }

  /* sort the groups according to group-number */
  parallel_sort(Group, Ngroups, sizeof(struct group_properties), fof_compare_Group_GrNr);

  for(i = 0; i < NumPart; i++)
    PS[i].GrNr = TotNgroups + 1; /* this marks all particles that are not in any group */

  for(i = 0, start = 0, Nids = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID < FOF_GList[i].MinID)
        {
          start++;
          if(start > NumPart)
            terminate("start > NumPart");
        }

      if(FOF_PList[start].MinID != FOF_GList[i].MinID)
        terminate("FOF_PList[start=%d].MinID=%lld != FOF_GList[i=%d].MinID=%lld", start, (long long)FOF_PList[start].MinID, i,
                  (long long)FOF_GList[i].MinID);

      for(lenloc = 0; start + lenloc < NumPart;)
        if(FOF_PList[start + lenloc].MinID == FOF_GList[i].MinID)
          {
            PS[FOF_PList[start + lenloc].Pindex].GrNr = FOF_GList[i].GrNr;
            Nids++;
            lenloc++;
          }
        else
          break;

      start += lenloc;
    }

  sumup_large_ints(1, &Nids, &totNids);

  if(totNids != TotNids)
    {
      char buf[1000];
      sprintf(buf, "Task=%d Nids=%d totNids=%d TotNids=%d\n", ThisTask, Nids, (int)totNids, (int)TotNids);
      terminate(buf);
    }

  t1 = second();

  mpi_printf("FOF: Assigning of group numbers took = %g sec\n", timediff(t0, t1));
}

/*! \brief Computes all kind of properties of groups.
 *
 *  Not complete after calling this. There is still the function
 *  fof_finish_group_properties, which finalizes the calculation
 *  (with normalization, averages, unit conversions and other operations).
 *
 *  \param[in] gr Index in Group array.
 *  \param[in] start Start index in FOF_PList.
 *  \param[in] len Number of particles in this group.
 *
 *  \return void
 */
void fof_compute_group_properties(int gr, int start, int len)
{
  int j, k, index, type, start_index = FOF_PList[start].Pindex;
  double xyz[3];

  Group[gr].Len  = 0;
  double gr_Mass = 0;
#ifdef USE_SFR
  double gr_Sfr = 0;
#endif /* #ifdef USE_SFR */

  double gr_CM[3], gr_Vel[3];
  for(k = 0; k < 3; k++)
    {
      gr_CM[k]              = 0;
      gr_Vel[k]             = 0;
      Group[gr].FirstPos[k] = P[start_index].Pos[k];
    }

  double gr_MassType[NTYPES];
  for(k = 0; k < NTYPES; k++)
    {
      Group[gr].LenType[k] = 0;
      gr_MassType[k]       = 0;
    }

  // calculate
  for(k = 0; k < len; k++)
    {
      index = FOF_PList[start + k].Pindex;

      Group[gr].Len++;
      gr_Mass += P[index].Mass;
      type = P[index].Type;

      Group[gr].LenType[type]++;

      gr_MassType[type] += P[index].Mass;

#ifdef USE_SFR
      if(P[index].Type == 0)
        gr_Sfr += SphP[index].Sfr;
#endif /* #ifdef USE_SFR */

      for(j = 0; j < 3; j++)
        {
          xyz[j] = P[index].Pos[j];
          xyz[j] = fof_periodic(xyz[j] - P[start_index].Pos[j]);
          gr_CM[j] += P[index].Mass * xyz[j];
          gr_Vel[j] += P[index].Mass * P[index].Vel[j];
        }
    }

  // put values into group struct
  Group[gr].Mass = gr_Mass;
#ifdef USE_SFR
  Group[gr].Sfr = gr_Sfr;
#endif /* #ifdef USE_SFR */

  for(k = 0; k < 3; k++)
    {
      Group[gr].CM[k]  = gr_CM[k];
      Group[gr].Vel[k] = gr_Vel[k];
    }

  for(k = 0; k < NTYPES; k++)
    Group[gr].MassType[k] = gr_MassType[k];
}

/*! \brief Global exchange of identified groups to their appropriate task.
 *
 *  \return void
 */
void fof_exchange_group_data(void)
{
  struct group_properties *get_Group;
  int i, j, ngrp, recvTask, nimport, start;
  double xyz[3];

  /* sort the groups according to task */
  mysort(Group, NgroupsExt, sizeof(struct group_properties), fof_compare_Group_MinIDTask);

  /* count how many we have of each task */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;
  for(i = 0; i < NgroupsExt; i++)
    Send_count[FOF_GList[i].MinIDTask]++;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      if(j == ThisTask) /* we will not exchange the ones that are local */
        Recv_count[j] = 0;
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  get_Group = (struct group_properties *)mymalloc("get_Group", sizeof(struct group_properties) * nimport);

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the group data */
              MPI_Sendrecv(&Group[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct group_properties), MPI_BYTE, recvTask,
                           TAG_DENS_A, &get_Group[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct group_properties),
                           MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  /* sort the groups again according to MinID */
  mysort(Group, NgroupsExt, sizeof(struct group_properties), fof_compare_Group_MinID);
  mysort(get_Group, nimport, sizeof(struct group_properties), fof_compare_Group_MinID);

  /* now add in the partial imported group data to the main ones */
  for(i = 0, start = 0; i < nimport; i++)
    {
      while(Group[start].MinID < get_Group[i].MinID)
        {
          start++;
          if(start >= NgroupsExt)
            terminate("start >= NgroupsExt");
        }

      Group[start].Len += get_Group[i].Len;
      Group[start].Mass += get_Group[i].Mass;

      for(j = 0; j < NTYPES; j++)
        {
          Group[start].LenType[j] += get_Group[i].LenType[j];
          Group[start].MassType[j] += get_Group[i].MassType[j];
        }

#ifdef USE_SFR
      Group[start].Sfr += get_Group[i].Sfr;
#endif /* #ifdef USE_SFR */

      for(j = 0; j < 3; j++)
        {
          xyz[j] = get_Group[i].CM[j] / get_Group[i].Mass;
          xyz[j] = fof_periodic(xyz[j] + get_Group[i].FirstPos[j] - Group[start].FirstPos[j]);
          Group[start].CM[j] += get_Group[i].Mass * xyz[j];
          Group[start].Vel[j] += get_Group[i].Vel[j];
        }
    }

  myfree(get_Group);
}

/*! \brief Finalizes group property calculation.
 *
 *  Called after a loop over all particles of a group is already completed.
 *
 *  \return void
 */
void fof_finish_group_properties(void)
{
  double cm[3];
  int i, j, ngr;

  for(i = 0; i < NgroupsExt; i++)
    {
      if(Group[i].MinIDTask == ThisTask)
        {
          for(j = 0; j < 3; j++)
            {
              Group[i].Vel[j] /= Group[i].Mass;
              cm[j]          = Group[i].CM[j] / Group[i].Mass;
              cm[j]          = fof_periodic_wrap(cm[j] + Group[i].FirstPos[j]);
              Group[i].CM[j] = cm[j];
            }
        }
    }

  /* eliminate the non-local groups */
  for(i = 0, ngr = NgroupsExt; i < ngr; i++)
    {
      if(Group[i].MinIDTask != ThisTask)
        {
          Group[i] = Group[ngr - 1];
          i--;
          ngr--;
        }
    }

  if(ngr != Ngroups)
    terminate("ngr != Ngroups");

  mysort(Group, Ngroups, sizeof(struct group_properties), fof_compare_Group_MinID);
}

/*! \brief Do periodic wrap for coordinate.
 *
 *  Note that his works only for cubic box.
 *
 *  \param[in] x Coordinate.
 *
 *  \return coordinate within [-0.5*BoxSize,0.5*BoxSize).
 */
double fof_periodic(double x)
{
#ifndef GRAVITY_NOT_PERIODIC
  if(x >= 0.5 * All.BoxSize)
    x -= All.BoxSize;
  if(x < -0.5 * All.BoxSize)
    x += All.BoxSize;
#endif /* #ifndef GRAVITY_NOT_PERIODIC */
  return x;
}

/*! \brief Do periodic wrap for coordinate.
 *
 *  Note that his works only for cubic box.
 *
 *  \param[in] x Coordinate.
 *
 *  \return coordinate within [0,BoxSize).
 */
double fof_periodic_wrap(double x)
{
#ifndef GRAVITY_NOT_PERIODIC
  while(x >= All.BoxSize)
    x -= All.BoxSize;
  while(x < 0)
    x += All.BoxSize;
#endif /* #ifndef GRAVITY_NOT_PERIODIC */
  return x;
}

#endif /* of FOF */
