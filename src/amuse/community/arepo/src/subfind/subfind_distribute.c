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
 * \file        src/subfind/subfind_distribute.c
 * \date        05/2018
 * \brief       Moves grops and particles across MPI tasks form their
 *              simulation ordering to a subfind ordering.
 * \details     contains functions:
 *                void subfind_distribute_groups(void)
 *                void subfind_distribute_particles(MPI_Comm Communicator)
 *                void subfind_reorder_P(int *Id, int Nstart, int N)
 *                void subfind_reorder_PS(int *Id, int Nstart, int N)
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

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../fof/fof.h"
#include "subfind.h"

#ifdef SUBFIND
static struct group_properties *send_Group;

/*! \brief Distributes groups equally on MPI tasks.
 *
 *  \return void
 */
void subfind_distribute_groups(void)
{
  int i, nexport, nimport, target, ngrp, recvTask;

  /* count how many we have of each task */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < Ngroups; i++)
    {
      target = Group[i].TargetTask;

      if(target < 0 || target >= NTask)
        terminate("target < 0 || target >= NTask");

      if(target != ThisTask)
        Send_count[target]++;
    }

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, nexport = 0, nimport = 0, Recv_offset[0] = Send_offset[0] = 0; i < NTask; i++)
    {
      nimport += Recv_count[i];
      nexport += Send_count[i];

      if(i > 0)
        {
          Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
          Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
        }
    }

  send_Group = (struct group_properties *)mymalloc_movable(&send_Group, "send_Group", nexport * sizeof(struct group_properties));

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < Ngroups; i++)
    {
      target = Group[i].TargetTask;

      if(target != ThisTask)
        {
          send_Group[Send_offset[target] + Send_count[target]] = Group[i];
          Send_count[target]++;

          Group[i] = Group[Ngroups - 1];
          Ngroups--;
          i--;
        }
    }

  if(Ngroups + nimport > MaxNgroups)
    {
#ifdef VERBOSE
      printf("SUBFIND: Task=%d: (Ngroups=%d) + (nimport=%d) > (MaxNgroups=%d). Will increase MaxNgroups.\n", ThisTask, Ngroups,
             nimport, MaxNgroups);
#endif /* #ifdef VERBOSE */
      MaxNgroups = Ngroups + nimport;
      Group      = (struct group_properties *)myrealloc_movable(Group, sizeof(struct group_properties) * MaxNgroups);
    }

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the group info */
              MPI_Sendrecv(&send_Group[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct group_properties), MPI_BYTE,
                           recvTask, TAG_DENS_A, &Group[Ngroups + Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct group_properties), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  Ngroups += nimport;

  myfree_movable(send_Group);
}

static struct particle_data *partBuf;
static struct subfind_data *subBuf;

/* \brief Distributes particles on MPI tasks.
 *
 *  This function redistributes the particles in P[] and PS[] according to what
 *  is stored in PS[].TargetTask, and PS[].TargetIndex. NOTE: The associated
 *  SphP[] is not moved, i.e. the association is broken until the particles are
 *  moved back into the original order!
 *
 *  \param[in] Communicator MPI communicator.
 *
 *  \return void
 */
void subfind_distribute_particles(MPI_Comm Communicator)
{
  int nimport, nexport;
  int i, j, n, ngrp, target;
  int max_load, load;
  int CommThisTask, CommNTask;

  MPI_Comm_size(Communicator, &CommNTask);
  MPI_Comm_rank(Communicator, &CommThisTask);

  for(n = 0; n < CommNTask; n++)
    Send_count[n] = 0;

  for(n = 0; n < NumPart; n++)
    {
      target = PS[n].TargetTask;

      if(target != CommThisTask)
        {
          if(target < 0 || target >= CommNTask)
            terminate("n=%d targettask=%d", n, target);

          Send_count[target]++;
        }
    }

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);

  for(j = 0, nimport = 0, nexport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < CommNTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* for resize */
  load = (NumPart + nimport - nexport);
  MPI_Allreduce(&load, &max_load, 1, MPI_INT, MPI_MAX, Communicator);

  partBuf = (struct particle_data *)mymalloc_movable(&partBuf, "partBuf", nexport * sizeof(struct particle_data));
  subBuf  = (struct subfind_data *)mymalloc_movable(&subBuf, "subBuf", nexport * sizeof(struct subfind_data));

  for(i = 0; i < CommNTask; i++)
    Send_count[i] = 0;

  for(n = 0; n < NumPart; n++)
    {
      target = PS[n].TargetTask;

      if(target != CommThisTask)
        {
          partBuf[Send_offset[target] + Send_count[target]] = P[n];
          subBuf[Send_offset[target] + Send_count[target]]  = PS[n];

          P[n]  = P[NumPart - 1];
          PS[n] = PS[NumPart - 1];

          Send_count[target]++;
          NumPart--;
          n--;
        }
    }

  /* do resize */
  if(max_load > (1.0 - ALLOC_TOLERANCE) * All.MaxPart)
    {
      All.MaxPart = max_load / (1.0 - 2 * ALLOC_TOLERANCE);
      reallocate_memory_maxpart();
      PS = (struct subfind_data *)myrealloc_movable(PS, All.MaxPart * sizeof(struct subfind_data));
    }

  for(i = 0; i < CommNTask; i++)
    Recv_offset[i] += NumPart;

#ifndef NO_ISEND_IRECV_IN_DOMAIN

  MPI_Request *requests = (MPI_Request *)mymalloc("requests", 8 * CommNTask * sizeof(MPI_Request));
  int n_requests        = 0;

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = CommThisTask ^ ngrp;

      if(target < CommNTask)
        {
          if(Recv_count[target] > 0)
            {
              MPI_Irecv(P + Recv_offset[target], Recv_count[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA,
                        Communicator, &requests[n_requests++]);
              MPI_Irecv(PS + Recv_offset[target], Recv_count[target] * sizeof(struct subfind_data), MPI_BYTE, target, TAG_KEY,
                        Communicator, &requests[n_requests++]);
            }
        }
    }

  MPI_Barrier(Communicator); /* not really necessary, but this will guarantee that all receives are
                                posted before the sends, which helps the stability of MPI on
                                bluegene, and perhaps some mpich1-clusters */

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = CommThisTask ^ ngrp;

      if(target < CommNTask)
        {
          if(Send_count[target] > 0)
            {
              MPI_Isend(partBuf + Send_offset[target], Send_count[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA,
                        Communicator, &requests[n_requests++]);
              MPI_Isend(subBuf + Send_offset[target], Send_count[target] * sizeof(struct subfind_data), MPI_BYTE, target, TAG_KEY,
                        Communicator, &requests[n_requests++]);
            }
        }
    }

  MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);

#else  /* #ifndef NO_ISEND_IRECV_IN_DOMAIN */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = CommThisTask ^ ngrp;

      if(target < CommNTask)
        {
          if(Send_count[target] > 0 || Recv_count[target] > 0)
            {
              MPI_Sendrecv(partBuf + Send_offset[target], Send_count[target] * sizeof(struct particle_data), MPI_BYTE, target,
                           TAG_PDATA, P + Recv_offset[target], Recv_count[target] * sizeof(struct particle_data), MPI_BYTE, target,
                           TAG_PDATA, Communicator, MPI_STATUS_IGNORE);

              MPI_Sendrecv(subBuf + Send_offset[target], Send_count[target] * sizeof(struct subfind_data), MPI_BYTE, target, TAG_KEY,
                           PS + Recv_offset[target], Recv_count[target] * sizeof(struct subfind_data), MPI_BYTE, target, TAG_KEY,
                           Communicator, MPI_STATUS_IGNORE);
            }
        }
    }
#endif /* #ifndef NO_ISEND_IRECV_IN_DOMAIN #else */

  NumPart += nimport;
  myfree_movable(subBuf);
  myfree_movable(partBuf);

  /* finally, let's also address the desired local order according to PS[].TargetIndex */

  struct fof_local_sort_data *mp;
  int *Id;

  mp = (struct fof_local_sort_data *)mymalloc("mp", sizeof(struct fof_local_sort_data) * (NumPart));
  Id = (int *)mymalloc("Id", sizeof(int) * (NumPart));

  for(i = 0; i < NumPart; i++)
    {
      mp[i].index       = i;
      mp[i].targetindex = PS[i].TargetIndex;
    }

  qsort(mp, NumPart, sizeof(struct fof_local_sort_data), fof_compare_local_sort_data_targetindex);

  for(i = 0; i < NumPart; i++)
    Id[mp[i].index] = i;

  subfind_reorder_P(Id, 0, NumPart);

  for(i = 0; i < NumPart; i++)
    Id[mp[i].index] = i;

  subfind_reorder_PS(Id, 0, NumPart);

  myfree(Id);
  myfree(mp);
}

/*! \brief Reorders elements in the P array.
 *
 * \param[in] Id Array containing ordering.
 * \param[in] Nstart Start index (in Id and P).
 * \param[in] N Final element index + 1.
 *
 *  \return void
 */
void subfind_reorder_P(int *Id, int Nstart, int N)
{
  int i;
  struct particle_data Psave, Psource;
  int idsource, idsave, dest;

  for(i = Nstart; i < N; i++)
    {
      if(Id[i] != i)
        {
          Psource  = P[i];
          idsource = Id[i];

          dest = Id[i];

          do
            {
              Psave  = P[dest];
              idsave = Id[dest];

              P[dest]  = Psource;
              Id[dest] = idsource;

              if(dest == i)
                break;

              Psource  = Psave;
              idsource = idsave;

              dest = idsource;
            }
          while(1);
        }
    }
}

/*! \brief Reorders elements in the PS array.
 *
 * \param[in] Id Array containing ordering.
 * \param[in] Nstart Start index (in Id and P).
 * \param[in] N Final element index + 1.
 *
 *  \return void
 */
void subfind_reorder_PS(int *Id, int Nstart, int N)
{
  int i;
  struct subfind_data PSsave, PSsource;
  int idsource, idsave, dest;

  for(i = Nstart; i < N; i++)
    {
      if(Id[i] != i)
        {
          PSsource = PS[i];

          idsource = Id[i];
          dest     = Id[i];

          do
            {
              PSsave = PS[dest];
              idsave = Id[dest];

              PS[dest] = PSsource;
              Id[dest] = idsource;

              if(dest == i)
                break;

              PSsource = PSsave;
              idsource = idsave;

              dest = idsource;
            }
          while(1);
        }
    }
}

#endif /* #ifdef SUBFIND */
