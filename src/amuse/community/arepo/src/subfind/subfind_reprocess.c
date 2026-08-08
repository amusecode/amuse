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
 * \file        src/subfind/subfind_fof_reprocess.c
 * \date        05/2018
 * \brief       Routines to calculate additional group properties.
 * \details     contains functions:
 *                void subfind_add_grp_props_calc_fof_angular_momentum(int num,
 *                  int ngroups_cat)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 14.05.2018 Prepared file for public release -- Rainer Weinberger
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
#include "../fof/fof.h"
#include "subfind.h"

#ifdef SUBFIND_EXTENDED_PROPERTIES
/*! \brief Angular Momentum calculation for groups.
 *
 *  \param[in] num Index of snapshot.
 *  \param[in] ngroups_cat Number of groups in group file.
 *
 *  \return void
 */
void subfind_add_grp_props_calc_fof_angular_momentum(int num, int ngroups_cat)
{
  mpi_printf("FOF: Begin Angular Momentum Calculation for FOF Groups.\n");

  /* assign target CPUs for the particles in groups */
  /* the particles not in groups will be distributed such that a uniform particle load results */
  double t0           = second();
  int *count_loc_task = mymalloc_clear("count_loc_task", NTask * sizeof(int));
  int *count_task     = mymalloc("count_task", NTask * sizeof(int));
  int *count_free     = mymalloc("count_free", NTask * sizeof(int));
  int count_loc_free  = 0;

  for(int i = 0; i < NumPart; i++)
    {
      if(PS[i].GrNr < 0)
        terminate("PS[i].GrNr=%d", PS[i].GrNr);

      if(PS[i].GrNr < TotNgroups) /* particle is in a group */
        {
          if(PS[i].GrNr < Ncollective) /* we are in a collective group */
            PS[i].TargetTask = ProcAssign[PS[i].GrNr].FirstTask + (i % ProcAssign[PS[i].GrNr].NTask);
          else
            PS[i].TargetTask = ((PS[i].GrNr - Ncollective) % (NTask - NprocsCollective)) + NprocsCollective;

          if(PS[i].TargetTask < 0 || PS[i].TargetTask >= NTask)
            terminate("PS[i].TargetTask=%d PS[i].GrNr=%d", PS[i].TargetTask, PS[i].GrNr);

          count_loc_task[PS[i].TargetTask]++;
        }
      else
        count_loc_free++;

      PS[i].TargetIndex = 0; /* unimportant here */
    }

  MPI_Allgather(&count_loc_free, 1, MPI_INT, count_free, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allreduce(count_loc_task, count_task, NTask, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  long long sum = 0;
  for(int i = 0; i < NTask; i++)
    sum += count_task[i] + count_free[i];

  int maxload = (sum + NTask - 1) / NTask;
  for(int i = 0; i < NTask; i++)
    {
      count_task[i] = maxload - count_task[i]; /* this is the amount that can fit on this task */
      if(count_task[i] < 0)
        count_task[i] = 0;
    }

  int current_task = 0;

  for(int i = 0; i < ThisTask; i++)
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

  for(int i = 0; i < NumPart; i++)
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

  double balance = subfind_get_particle_balance();
  mpi_printf("SUBFIND: particle balance=%g\n", balance);

  /* distribute particles such that groups are completely on the CPU(s) that do the corresponding group(s) */
  fof_subfind_exchange(MPI_COMM_WORLD);
  double t1 = second();
  mpi_printf("SUBFIND: subfind_exchange() took %g sec\n", timediff(t0, t1));

  balance = subfind_get_particle_balance();
  mpi_printf("SUBFIND: particle balance for AM processing=%g\n", balance);

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
      subfind_fof_calc_am_collective(num, ngroups_cat);
    }
  else
    {
      /* now let us sort according to GrNr and Density. This step will temporarily break the association with SphP[] and other arrays!
       */
      submp = (struct submp_data *)mymalloc("submp", sizeof(struct submp_data) * NumPart);
      for(int i = 0; i < NumPart; i++)
        {
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
        printf("SUBFIND-SERIAL: Start to do AM for %d small groups with serial subfind algorithm on %d processors (root-node=%d)\n",
               TotNgroups - Ncollective, SubNTask, ThisTask);

      /* we now apply a serial version of subfind to the local groups */

      t0 = second();
      for(int gr = 0, offset = 0; gr < Ngroups; gr++)
        {
          if(((Group[gr].GrNr - Ncollective) % (NTask - NprocsCollective)) + NprocsCollective == ThisTask)
            offset = subfind_fof_calc_am_serial(gr, offset, num, ngroups_cat);
          else
            terminate("how come that we have this group number?");
        }

      MPI_Barrier(SubComm);
      t1 = second();
      if(SubThisTask == 0)
        printf("SUBFIND-SERIAL: processing AM of serial groups took %g sec\n", timediff(t0, t1));

      /* undo local rearrangement that made groups consecutive. After that, the association of SphP[] will be correct again */
      submp = (struct submp_data *)mymalloc("submp", sizeof(struct submp_data) * NumPart);
      for(int i = 0; i < NumPart; i++)
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

  /* distribute particles back to original CPU */
  t0 = second();
  for(int i = 0; i < NumPart; i++)
    {
      PS[i].TargetTask  = PS[i].OriginTask;
      PS[i].TargetIndex = PS[i].OriginIndex;
    }

  fof_subfind_exchange(MPI_COMM_WORLD);
  t1 = second();
  if(ThisTask == 0)
    printf("SUBFIND: subfind_exchange() (for return to original CPU after AM)  took %g sec\n", timediff(t0, t1));

  mpi_printf("FOF: Angular Momentum Calculation for FOF Groups finished successfully.\n");
}
#endif /* #ifdef SUBFIND_EXTENDED_PROPERTIES */
