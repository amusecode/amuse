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
 * \file        src/domain_exchange.c
 * \date        05/2018
 * \brief       Algorithms for exchanging particle data and associated
 *              rearrangements.
 * \details     This includes changing the size of the P and SphP arrays as
 *              well as the particle exchange routine itself.
 *              contains functions:
 *                void domain_resize_storage(int count_get, int count_get_sph,
 *                  int option_flag)
 *                void domain_exchange(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 05.05.2018 Prepared file for public release -- Rainer Weinberger
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
#include "domain.h"

/*! \brief Changes memory allocation if necessary for particle and cell data.
 *
 *  If the memory usage due to a net import or export of particles changes
 *  above a certain tolerance, the P and SphP structures need to be
 *  reallocated.
 *
 *  \param[in] count get How many particles are imported?
 *  \param[in] count_get_sph How many cells are imported?
 *  \param[in] option_flag Options for reallocating peanokey or ngbtree.
 *
 *  \return void
 */
void domain_resize_storage(int count_get, int count_get_sph, int option_flag)
{
  int load        = NumPart + count_get;
  int sphload     = NumGas + count_get_sph;
  int loc_data[2] = {load, sphload}, res[2];

  MPI_Allreduce(loc_data, res, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  int max_load    = res[0];
  int max_sphload = res[1];

  if(max_load > (1.0 - ALLOC_TOLERANCE) * All.MaxPart || max_load < (1.0 - 3 * ALLOC_TOLERANCE) * All.MaxPart)
    {
      All.MaxPart = max_load / (1.0 - 2 * ALLOC_TOLERANCE);
      reallocate_memory_maxpart();

      if(option_flag == 1)
        Key = (peanokey *)myrealloc_movable(Key, sizeof(peanokey) * All.MaxPart);
    }

  if(max_sphload >= (1.0 - ALLOC_TOLERANCE) * All.MaxPartSph || max_sphload < (1.0 - 3 * ALLOC_TOLERANCE) * All.MaxPartSph)
    {
      All.MaxPartSph = max_sphload / (1.0 - 2 * ALLOC_TOLERANCE);
      if(option_flag == 2)
        {
          if(All.MaxPartSph > Ngb_MaxPart)
            ngb_treemodifylength(All.MaxPartSph - Ngb_MaxPart);
        }
      reallocate_memory_maxpartsph();
    }
}

/*! \brief Exchanges particles and cells according to new domain decomposition.
 *
 *  Communicates particles and cells to their new task. P and SphP arrays are
 *  changed in size accordingly.
 *
 *  \return void
 */
void domain_exchange(void)
{
  double t0 = second();

  int count_togo = 0, count_togo_sph = 0, count_get = 0, count_get_sph = 0;
  int *count, *count_sph, *offset, *offset_sph;
  int *count_recv, *count_recv_sph, *offset_recv, *offset_recv_sph;
  int i, n, no, target;
  struct particle_data *partBuf;
  struct sph_particle_data *sphBuf;

  peanokey *keyBuf;

  long long sumtogo = 0;

  for(i = 0; i < NTask; i++)
    sumtogo += toGo[i];

  sumup_longs(1, &sumtogo, &sumtogo);

  count           = (int *)mymalloc_movable(&count, "count", NTask * sizeof(int));
  count_sph       = (int *)mymalloc_movable(&count_sph, "count_sph", NTask * sizeof(int));
  offset          = (int *)mymalloc_movable(&offset, "offset", NTask * sizeof(int));
  offset_sph      = (int *)mymalloc_movable(&offset_sph, "offset_sph", NTask * sizeof(int));
  count_recv      = (int *)mymalloc_movable(&count_recv, "count_recv", NTask * sizeof(int));
  count_recv_sph  = (int *)mymalloc_movable(&count_recv_sph, "count_recv_sph", NTask * sizeof(int));
  offset_recv     = (int *)mymalloc_movable(&offset_recv, "offset_recv", NTask * sizeof(int));
  offset_recv_sph = (int *)mymalloc_movable(&offset_recv_sph, "offset_recv_sph", NTask * sizeof(int));

  int prec_offset;
  int *decrease;

  decrease = (int *)mymalloc_movable(&decrease, "decrease", NTask * sizeof(int));

  for(i = 1, offset_sph[0] = 0, decrease[0] = 0; i < NTask; i++)
    {
      offset_sph[i] = offset_sph[i - 1] + toGoSph[i - 1];
      decrease[i]   = toGoSph[i - 1];
    }

  prec_offset = offset_sph[NTask - 1] + toGoSph[NTask - 1];

  offset[0] = prec_offset;
  for(i = 1; i < NTask; i++)
    offset[i] = offset[i - 1] + (toGo[i - 1] - decrease[i]);

  myfree(decrease);

  for(i = 0; i < NTask; i++)
    {
      count_togo += toGo[i];
      count_togo_sph += toGoSph[i];
      count_get += toGet[i];
      count_get_sph += toGetSph[i];
    }

  partBuf = (struct particle_data *)mymalloc_movable(&partBuf, "partBuf", count_togo * sizeof(struct particle_data));
  sphBuf  = (struct sph_particle_data *)mymalloc_movable(&sphBuf, "sphBuf", count_togo_sph * sizeof(struct sph_particle_data));

  keyBuf = (peanokey *)mymalloc_movable(&keyBuf, "keyBuf", count_togo * sizeof(peanokey));

  for(i = 0; i < NTask; i++)
    {
      count[i] = count_sph[i] = 0;
    }

  for(n = 0; n < NumPart; n++)
    {
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

      target = DomainTask[no];

      if(target != ThisTask)
        {
          /* copy this particle into the exchange buffer */
          if(P[n].Type == 0)
            {
              partBuf[offset_sph[target] + count_sph[target]] = P[n];
              keyBuf[offset_sph[target] + count_sph[target]]  = Key[n];
              sphBuf[offset_sph[target] + count_sph[target]]  = SphP[n];
              count_sph[target]++;
            }
          else
            {
              partBuf[offset[target] + count[target]] = P[n];
              keyBuf[offset[target] + count[target]]  = Key[n];
              count[target]++;
            }

          if(P[n].Type == 0)
            {
              P[n]          = P[NumGas - 1];
              P[NumGas - 1] = P[NumPart - 1];

              Key[n]          = Key[NumGas - 1];
              Key[NumGas - 1] = Key[NumPart - 1];

              SphP[n] = SphP[NumGas - 1];

              NumGas--;
            }
          else
            {
              P[n]   = P[NumPart - 1];
              Key[n] = Key[NumPart - 1];
            }

          NumPart--;
          n--;

        } /* target != ThisTask */
    }     /* n < NumPart */

  /**** now resize the storage for the P[] and SphP[] arrays if needed ****/
  domain_resize_storage(count_get, count_get_sph, 1);

  /*****  space has been created, now can do the actual exchange *****/
  int count_totget = count_get_sph;

  if(count_totget)
    {
      memmove(P + NumGas + count_totget, P + NumGas, (NumPart - NumGas) * sizeof(struct particle_data));
      memmove(Key + NumGas + count_totget, Key + NumGas, (NumPart - NumGas) * sizeof(peanokey));
    }

  for(i = 0; i < NTask; i++)
    {
      count_recv_sph[i] = toGetSph[i];
      count_recv[i]     = toGet[i] - toGetSph[i];
    }

  int prec_count;
  for(i = 1, offset_recv_sph[0] = NumGas; i < NTask; i++)
    offset_recv_sph[i] = offset_recv_sph[i - 1] + count_recv_sph[i - 1];
  prec_count = NumGas + count_get_sph;

  offset_recv[0] = NumPart - NumGas + prec_count;

  for(i = 1; i < NTask; i++)
    offset_recv[i] = offset_recv[i - 1] + count_recv[i - 1];

#ifndef USE_MPIALLTOALLV_IN_DOMAINDECOMP

  int ngrp;
#ifdef NO_ISEND_IRECV_IN_DOMAIN /* synchronous communication */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
        {
          if(count_sph[target] > 0 || count_recv_sph[target] > 0)
            {
              MPI_Sendrecv(partBuf + offset_sph[target], count_sph[target] * sizeof(struct particle_data), MPI_BYTE, target,
                           TAG_PDATA_SPH, P + offset_recv_sph[target], count_recv_sph[target] * sizeof(struct particle_data), MPI_BYTE,
                           target, TAG_PDATA_SPH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              MPI_Sendrecv(sphBuf + offset_sph[target], count_sph[target] * sizeof(struct sph_particle_data), MPI_BYTE, target,
                           TAG_SPHDATA, SphP + offset_recv_sph[target], count_recv_sph[target] * sizeof(struct sph_particle_data),
                           MPI_BYTE, target, TAG_SPHDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              MPI_Sendrecv(keyBuf + offset_sph[target], count_sph[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_SPH,
                           Key + offset_recv_sph[target], count_recv_sph[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_SPH,
                           MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

          if(count[target] > 0 || count_recv[target] > 0)
            {
              MPI_Sendrecv(partBuf + offset[target], count[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA,
                           P + offset_recv[target], count_recv[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA,
                           MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              MPI_Sendrecv(keyBuf + offset[target], count[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY,
                           Key + offset_recv[target], count_recv[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

#else  /* #ifdef NO_ISEND_IRECV_IN_DOMAIN */
  /* asynchronous communication */

  MPI_Request *requests = (MPI_Request *)mymalloc_movable(&requests, "requests", 30 * NTask * sizeof(MPI_Request));
  int n_requests        = 0;

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
        {
          if(count_recv_sph[target] > 0)
            {
              MPI_Irecv(P + offset_recv_sph[target], count_recv_sph[target] * sizeof(struct particle_data), MPI_BYTE, target,
                        TAG_PDATA_SPH, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Irecv(SphP + offset_recv_sph[target], count_recv_sph[target] * sizeof(struct sph_particle_data), MPI_BYTE, target,
                        TAG_SPHDATA, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Irecv(Key + offset_recv_sph[target], count_recv_sph[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_SPH,
                        MPI_COMM_WORLD, &requests[n_requests++]);
            }

          if(count_recv[target] > 0)
            {
              MPI_Irecv(P + offset_recv[target], count_recv[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA,
                        MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Irecv(Key + offset_recv[target], count_recv[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY, MPI_COMM_WORLD,
                        &requests[n_requests++]);
            }
        }
    }

  MPI_Barrier(MPI_COMM_WORLD); /* not really necessary, but this will guarantee that all receives are
                                  posted before the sends, which helps the stability of MPI on
                                  bluegene, and perhaps some mpich1-clusters */

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
        {
          if(count_sph[target] > 0)
            {
              MPI_Isend(partBuf + offset_sph[target], count_sph[target] * sizeof(struct particle_data), MPI_BYTE, target,
                        TAG_PDATA_SPH, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Isend(sphBuf + offset_sph[target], count_sph[target] * sizeof(struct sph_particle_data), MPI_BYTE, target,
                        TAG_SPHDATA, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Isend(keyBuf + offset_sph[target], count_sph[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_SPH,
                        MPI_COMM_WORLD, &requests[n_requests++]);
            }

          if(count[target] > 0)
            {
              MPI_Isend(partBuf + offset[target], count[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA,
                        MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Isend(keyBuf + offset[target], count[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY, MPI_COMM_WORLD,
                        &requests[n_requests++]);
            }
        }
    }

  MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);
#endif /* #ifdef NO_ISEND_IRECV_IN_DOMAIN #else */

#else /* #ifndef USE_MPIALLTOALLV_IN_DOMAINDECOMP */
  /* begins block of myMPI_Alltoallv communications */

  myMPI_Alltoallv(partBuf, count_sph, offset_sph, P, count_recv_sph, offset_recv_sph, sizeof(struct particle_data), 0, MPI_COMM_WORLD);

  myMPI_Alltoallv(sphBuf, count_sph, offset_sph, SphP, count_recv_sph, offset_recv_sph, sizeof(struct sph_particle_data), 0,
                  MPI_COMM_WORLD);

  myMPI_Alltoallv(keyBuf, count_sph, offset_sph, Key, count_recv_sph, offset_recv_sph, sizeof(peanokey), 0, MPI_COMM_WORLD);

  myMPI_Alltoallv(partBuf, count, offset, P, count_recv, offset_recv, sizeof(struct particle_data), 0, MPI_COMM_WORLD);

  myMPI_Alltoallv(keyBuf, count, offset, Key, count_recv, offset_recv, sizeof(peanokey), 0, MPI_COMM_WORLD);

#endif /* #ifndef USE_MPIALLTOALLV_IN_DOMAINDECOMP #else */
       /* close block of myMPI_Alltoallv communications */

  NumPart += count_get;
  NumGas += count_get_sph;

  myfree(keyBuf);
  myfree(sphBuf);
  myfree(partBuf);
  myfree(offset_recv_sph);
  myfree(offset_recv);
  myfree(count_recv_sph);
  myfree(count_recv);
  myfree(offset_sph);
  myfree(offset);
  myfree(count_sph);
  myfree(count);

  double t1 = second();
  mpi_printf("DOMAIN: exchange of %lld particles done. (took %g sec)\n", sumtogo, timediff(t0, t1));
}
