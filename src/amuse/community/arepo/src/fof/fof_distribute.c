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
 * \file        src/fof/fof_distribute.c
 * \date        05/2018
 * \brief       Communication and reordering routines for FoF.
 * \details     contains functions:
 *                void fof_subfind_exchange(MPI_Comm Communicator)
 *                void fof_reorder_PS(int *Id, int Nstart, int N)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 24.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"
#include "../subfind/subfind.h"
#include "fof.h"

#ifdef FOF

/*! \brief Redistributes the particles according to what is stored in
 *         PS[].TargetTask, and PS[].TargetIndex.
 *
 *  \param[in] Communicator MPI communicator.
 *
 *  \return void
 */
void fof_subfind_exchange(MPI_Comm Communicator)
{
  int nimport, nexport;
  int i, j, n, type, ngrp, target;
  int max_load, max_loadsph, load;
  struct particle_data *partBuf;
  struct subfind_data *subBuf;
  struct sph_particle_data *sphBuf;

  int CommThisTask, CommNTask;

  MPI_Comm_size(Communicator, &CommNTask);
  MPI_Comm_rank(Communicator, &CommThisTask);

  int old_AllMaxPart    = All.MaxPart;
  int old_AllMaxPartSph = All.MaxPartSph;

  for(type = 0; type < NTYPES; type++)
    {
      size_t ExportSpace = 0.5 * (FreeBytes); /* we will try to grab at most half of the still available memory  */
      size_t PartSpace   = sizeof(struct particle_data) + sizeof(struct subfind_data) + sizeof(struct sph_particle_data);
      if(PartSpace > ExportSpace)
        terminate("seems like we have insufficient storage, PartSpace=%lld ExportSpace=%lld", (long long)PartSpace,
                  (long long)ExportSpace);

      int glob_flag = 0;

      do
        {
          for(n = 0; n < CommNTask; n++)
            {
              Send_count[n] = 0;
            }

          ptrdiff_t AvailableSpace = ExportSpace; /* this must be a type that can become negative */

          for(n = 0; n < NumPart; n++)
            {
              if(AvailableSpace < 0)
                break;

              if(P[n].Type == type && PS[n].TargetTask != CommThisTask)
                {
                  target = PS[n].TargetTask;

                  if(target < 0 || target >= CommNTask)
                    terminate("n=%d targettask=%d", n, target);

                  AvailableSpace -= PartSpace;

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

          if(type == 0)
            {
              load = (NumGas + nimport - nexport);
              MPI_Allreduce(&load, &max_loadsph, 1, MPI_INT, MPI_MAX, Communicator);
            }

          partBuf = (struct particle_data *)mymalloc_movable(&partBuf, "partBuf", nexport * sizeof(struct particle_data));
          subBuf  = (struct subfind_data *)mymalloc_movable(&subBuf, "subBuf", nexport * sizeof(struct subfind_data));
          if(type == 0)
            sphBuf = (struct sph_particle_data *)mymalloc_movable(&sphBuf, "sphBuf", nexport * sizeof(struct sph_particle_data));

          for(i = 0; i < CommNTask; i++)
            {
              Send_count[i] = 0;
            }

          AvailableSpace = ExportSpace; /* this must be allowed to become negative */

          int nstay         = 0;
          int delta_numpart = 0;
          int delta_numgas  = 0;

          for(n = 0; n < NumPart; n++)
            {
              if(AvailableSpace < 0)
                break;

              if(P[n].Type == type && PS[n].TargetTask != CommThisTask)
                {
                  target = PS[n].TargetTask;

                  AvailableSpace -= PartSpace;

                  partBuf[Send_offset[target] + Send_count[target]] = P[n];
                  subBuf[Send_offset[target] + Send_count[target]]  = PS[n];

                  if(P[n].Type == 0)
                    {
                      sphBuf[Send_offset[target] + Send_count[target]] = SphP[n];
                      delta_numgas++;
                    }

                  Send_count[target]++;
                  delta_numpart++;
                }
              else
                {
                  if(nstay != n)
                    {
                      /* now move P[n] to P[nstay] */

                      P[nstay]  = P[n];
                      PS[nstay] = PS[n];

                      if(P[nstay].Type == 0)
                        SphP[nstay] = SphP[n];
                    }

                  nstay++;
                }
            }

          if(delta_numgas > 0)
            if(delta_numpart != delta_numgas)
              terminate("delta_numpart=%d != delta_numgas=%d", delta_numpart, delta_numgas);

          /* now close gap (if present) */
          memmove(P + nstay, P + nstay + delta_numpart, (NumPart - (nstay + delta_numpart)) * sizeof(struct particle_data));
          memmove(PS + nstay, PS + nstay + delta_numpart, (NumPart - (nstay + delta_numpart)) * sizeof(struct subfind_data));

          if(delta_numgas > 0)
            if(NumGas - (nstay + delta_numgas) > 0)
              memmove(SphP + nstay, SphP + nstay + delta_numpart,
                      (NumGas - (nstay + delta_numgas)) * sizeof(struct sph_particle_data));

          NumPart -= delta_numpart;
          NumGas -= delta_numgas;

          /* do resize, but only increase arrays!! (otherwise data in ActiveParticleList etc. gets lost */
          if(max_load > (1.0 - ALLOC_TOLERANCE) * All.MaxPart)
            {
              All.MaxPart = max_load / (1.0 - 2 * ALLOC_TOLERANCE);
              reallocate_memory_maxpart();
              PS = (struct subfind_data *)myrealloc_movable(PS, All.MaxPart * sizeof(struct subfind_data));
            }

          if(type == 0)
            {
              if(max_loadsph > (1.0 - ALLOC_TOLERANCE) * All.MaxPartSph)
                {
                  All.MaxPartSph = max_loadsph / (1.0 - 2 * ALLOC_TOLERANCE);
                  reallocate_memory_maxpartsph();
                }
            }

          /* create a gap behind the existing gas particles where we will insert the incoming particles */
          memmove(P + NumGas + nimport, P + NumGas, (NumPart - NumGas) * sizeof(struct particle_data));
          memmove(PS + NumGas + nimport, PS + NumGas, (NumPart - NumGas) * sizeof(struct subfind_data));

          for(i = 0; i < CommNTask; i++)
            Recv_offset[i] += NumGas;

          for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
            {
              target = CommThisTask ^ ngrp;

              if(target < CommNTask)
                {
                  if(Send_count[target] > 0 || Recv_count[target] > 0)
                    {
                      MPI_Sendrecv(partBuf + Send_offset[target], Send_count[target] * sizeof(struct particle_data), MPI_BYTE, target,
                                   TAG_PDATA, P + Recv_offset[target], Recv_count[target] * sizeof(struct particle_data), MPI_BYTE,
                                   target, TAG_PDATA, Communicator, MPI_STATUS_IGNORE);

                      MPI_Sendrecv(subBuf + Send_offset[target], Send_count[target] * sizeof(struct subfind_data), MPI_BYTE, target,
                                   TAG_KEY, PS + Recv_offset[target], Recv_count[target] * sizeof(struct subfind_data), MPI_BYTE,
                                   target, TAG_KEY, Communicator, MPI_STATUS_IGNORE);

                      if(type == 0)
                        MPI_Sendrecv(sphBuf + Send_offset[target], Send_count[target] * sizeof(struct sph_particle_data), MPI_BYTE,
                                     target, TAG_SPHDATA, SphP + Recv_offset[target],
                                     Recv_count[target] * sizeof(struct sph_particle_data), MPI_BYTE, target, TAG_SPHDATA,
                                     Communicator, MPI_STATUS_IGNORE);
                    }
                }
            }

          if(type == 0)
            NumGas += nimport;

          NumPart += nimport;

          if(type == 0)
            myfree_movable(sphBuf);

          myfree_movable(subBuf);
          myfree_movable(partBuf);

          int loc_flag = 0;
          if(AvailableSpace < 0)
            loc_flag = 1;

          MPI_Allreduce(&loc_flag, &glob_flag, 1, MPI_INT, MPI_SUM, Communicator);
          if(glob_flag > 0 && CommThisTask == 0)
            {
              printf(
                  "FOF-DISTRIBUTE: Need to cycle in particle exchange due to memory shortage. type=%d glob_flag=%d ThisTask=%d "
                  "CommThisTask=%d   PartSpace=%lld  ExportSpace=%lld\n",
                  type, glob_flag, ThisTask, CommThisTask, (long long)PartSpace, (long long)ExportSpace);
              fflush(stdout);
            }
        }
      while(glob_flag);
    }

  /* if there was a temporary memory shortage during the exchange, we may had to increase the maximum allocations. Go back to smaller
   * values again if possible */

  load = NumPart;
  MPI_Allreduce(&load, &max_load, 1, MPI_INT, MPI_MAX, Communicator);
  max_load = max_load / (1.0 - 2 * ALLOC_TOLERANCE);
  if(max_load < old_AllMaxPart)
    max_load = old_AllMaxPart;
  if(max_load != All.MaxPart)
    {
      All.MaxPart = max_load;
      reallocate_memory_maxpart();
      PS = (struct subfind_data *)myrealloc_movable(PS, All.MaxPart * sizeof(struct subfind_data));
    }

  load = NumGas;
  MPI_Allreduce(&load, &max_loadsph, 1, MPI_INT, MPI_MAX, Communicator);
  max_loadsph = max_loadsph / (1.0 - 2 * ALLOC_TOLERANCE);
  if(max_loadsph < old_AllMaxPartSph)
    max_loadsph = old_AllMaxPartSph;
  if(max_loadsph != All.MaxPartSph)
    {
      All.MaxPartSph = max_loadsph;
      reallocate_memory_maxpartsph();
    }

  /* finally, let's also address the desired local order according to PS[].TargetIndex */

  struct fof_local_sort_data *mp;
  int *Id;

  if(NumGas)
    {
      mp = (struct fof_local_sort_data *)mymalloc("mp", sizeof(struct fof_local_sort_data) * NumGas);
      Id = (int *)mymalloc("Id", sizeof(int) * NumGas);

      for(i = 0; i < NumGas; i++)
        {
          mp[i].index       = i;
          mp[i].targetindex = PS[i].TargetIndex;
        }

      qsort(mp, NumGas, sizeof(struct fof_local_sort_data), fof_compare_local_sort_data_targetindex);

      for(i = 0; i < NumGas; i++)
        Id[mp[i].index] = i;

      reorder_gas(Id);

      for(i = 0; i < NumGas; i++)
        Id[mp[i].index] = i;

      fof_reorder_PS(Id, 0, NumGas);

      myfree(Id);
      myfree(mp);
    }

  if(NumPart - NumGas > 0)
    {
      mp = (struct fof_local_sort_data *)mymalloc("mp", sizeof(struct fof_local_sort_data) * (NumPart - NumGas));
      mp -= NumGas;

      Id = (int *)mymalloc("Id", sizeof(int) * (NumPart - NumGas));
      Id -= NumGas;

      for(i = NumGas; i < NumPart; i++)
        {
          mp[i].index       = i;
          mp[i].targetindex = PS[i].TargetIndex;
        }

      qsort(mp + NumGas, NumPart - NumGas, sizeof(struct fof_local_sort_data), fof_compare_local_sort_data_targetindex);

      for(i = NumGas; i < NumPart; i++)
        Id[mp[i].index] = i;

      reorder_particles(Id);

      for(i = NumGas; i < NumPart; i++)
        Id[mp[i].index] = i;

      fof_reorder_PS(Id, NumGas, NumPart);

      Id += NumGas;
      myfree(Id);
      mp += NumGas;
      myfree(mp);
    }
}

/*! \brief Reorders the elements in the PS array according to the indices given
 *         in the ID array.
 *
 *  \param[in, out] ID Array that specifies new index of element in PS array;
 *                  i.e. PS[i] -> PS[ ID[i] ].
 *  \param[in] Nstart Starting index in ID and PS arrays.
 *  \param[in] N Final element +1 in ID and PS arrays.
 *
 *  \return void
 */
void fof_reorder_PS(int *Id, int Nstart, int N)
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

#endif /* #ifdef FOF */
