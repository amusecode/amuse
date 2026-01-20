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
 * \file        src/gravdirect.c
 * \date        05/2018
 * \brief       Main driver routines for gravitational (short-range) force
 *              computation through direct summation
 * \details     Note that this is not the same thing as
 *              EXACT_GRAVITY_FOR_PARTICLE_TYPE!
 *              ALLOW_DIRECT_SUMMATION does direct summation for performance
 *              reasons if there is only a small number of interactions to be
 *              calculated and the overhead of a tree-construction would be
 *              more expensive than the direct summation calculation, while
 *              EXACT_GRAVITY_FOR_PARTICLE_TYPE always enforces a direct
 *              summation for all particle pairs of a given type.
 *              contains functions:
 *                void gravity_direct(int timebin)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 06.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

#ifdef ALLOW_DIRECT_SUMMATION
static int Nimport;

/*! \brief Computes the gravitational forces for all active particles through
 *         direct summation.
 *
 *  \param[in] timebin (unused)
 *
 *  \return void
 */
void gravity_direct(int timebin)
{
  int i, j, k, idx;

  TIMER_START(CPU_TREEDIRECT);

  if(TimeBinsGravity.GlobalNActiveParticles <= 1)
    {
      if(TimeBinsGravity.NActiveParticles > 0)
        {
          i = TimeBinsGravity.ActiveParticleList[0];
          if(i >= 0)
            {
              for(k = 0; k < 3; k++)
                P[i].GravAccel[k] = 0;

#ifdef EVALPOTENTIAL
              P[i].Potential = 0;
#endif /* #ifdef EVALPOTENTIAL */
            }
        }

      mpi_printf("Found only %d particles to do direct summation -> SKIPPING IT\n", TimeBinsGravity.GlobalNActiveParticles);
      TIMER_STOP(CPU_TREEDIRECT);
      return;
    }

  mpi_printf("GRAVDIRECT: direct summation.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  double tstart = second();

  DirectDataIn = (struct directdata *)mymalloc("DirectDataIn", TimeBinsGravity.NActiveParticles * sizeof(struct directdata));

  Nforces = 0;

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifdef CELL_CENTER_GRAVITY
      if(P[i].Type == 0)
        {
          for(k = 0; k < 3; k++)
            DirectDataIn[Nforces].Pos[k] = SphP[i].Center[k];
        }
      else
#endif /* #ifdef CELL_CENTER_GRAVITY */
        {
          for(k = 0; k < 3; k++)
            DirectDataIn[Nforces].Pos[k] = P[i].Pos[k];
        }

      DirectDataIn[Nforces].Mass = P[i].Mass;

      DirectDataIn[Nforces].Type          = P[i].Type;
      DirectDataIn[Nforces].SofteningType = P[i].SofteningType;

      Nforces++;
    }

  MPI_Allgather(&Nforces, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, Nimport = 0, Recv_offset[0] = 0; j < NTask; j++)
    {
      Nimport += Recv_count[j];

      if(j > 0)
        Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
    }

  DirectDataAll = (struct directdata *)mymalloc("DirectDataAll", Nimport * sizeof(struct directdata));

  for(j = 0; j < NTask; j++)
    {
      Send_count[j]  = Recv_count[j] * sizeof(struct directdata);
      Send_offset[j] = Recv_offset[j] * sizeof(struct directdata);
    }

  MPI_Allgatherv(DirectDataIn, Nforces * sizeof(struct directdata), MPI_BYTE, DirectDataAll, Send_count, Send_offset, MPI_BYTE,
                 MPI_COMM_WORLD);

  /* subdivide the work evenly */
  int first, count;
  subdivide_evenly(Nimport, NTask, ThisTask, &first, &count);

  DirectAccOut = (struct accdata *)mymalloc("DirectDataOut", count * sizeof(struct accdata));

  /* now calculate the forces */
  for(i = 0; i < count; i++)
    force_evaluate_direct(i + first, i, Nimport);

  /* now send the forces to the right places */

  DirectAccIn = (struct accdata *)mymalloc("DirectDataIn", Nforces * sizeof(struct accdata));

  MPI_Request *requests = (MPI_Request *)mymalloc_movable(&requests, "requests", 2 * NTask * sizeof(MPI_Request));
  int n_requests        = 0;

  int recvTask = 0;
  int sendTask = 0;
  int send_first, send_count;
  subdivide_evenly(Nimport, NTask, sendTask, &send_first, &send_count);

  while(recvTask < NTask && sendTask < NTask) /* go through both lists */
    {
      while(send_first + send_count < Recv_offset[recvTask])
        {
          if(sendTask >= NTask - 1)
            terminate("sendTask >= NTask  recvTask=%d sendTask=%d", recvTask, sendTask);

          sendTask++;
          subdivide_evenly(Nimport, NTask, sendTask, &send_first, &send_count);
        }

      while(Recv_offset[recvTask] + Recv_count[recvTask] < send_first)
        {
          if(recvTask >= NTask - 1)
            terminate("recvTask >= NTask  recvTask=%d sendTask=%d", recvTask, sendTask);

          recvTask++;
        }

      int start = imax(Recv_offset[recvTask], send_first);
      int next  = imin(Recv_offset[recvTask] + Recv_count[recvTask], send_first + send_count);

      if(next - start >= 1)
        {
          if(ThisTask == sendTask)
            MPI_Isend(DirectAccOut + start - send_first, (next - start) * sizeof(struct accdata), MPI_BYTE, recvTask, TAG_PDATA_SPH,
                      MPI_COMM_WORLD, &requests[n_requests++]);

          if(ThisTask == recvTask)
            MPI_Irecv(DirectAccIn + start - Recv_offset[recvTask], (next - start) * sizeof(struct accdata), MPI_BYTE, sendTask,
                      TAG_PDATA_SPH, MPI_COMM_WORLD, &requests[n_requests++]);
        }

      if(next == Recv_offset[recvTask] + Recv_count[recvTask])
        recvTask++;
      else
        {
          sendTask++;
          if(sendTask >= NTask)
            break;

          subdivide_evenly(Nimport, NTask, sendTask, &send_first, &send_count);
        }
    }

  MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);

  Nforces = 0;

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      for(k = 0; k < 3; k++)
        P[i].GravAccel[k] = DirectAccIn[Nforces].Acc[k];

#ifdef EVALPOTENTIAL
      P[i].Potential = DirectAccIn[Nforces].Potential;
#endif /* #ifdef EVALPOTENTIAL */
      Nforces++;
    }

  myfree(DirectAccIn);
  myfree(DirectAccOut);
  myfree(DirectDataAll);
  myfree(DirectDataIn);

  mpi_printf("GRAVDIRECT: force is done.\n");

  All.TotNumOfForces += TimeBinsGravity.GlobalNActiveParticles;

  double tend = second();

  double timedirect, sumt;
  timedirect = tend - tstart;

  MPI_Reduce(&timedirect, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdTimings, "Nf=%9lld   active part/task: avg=%g   total-Nf=%lld\n", TimeBinsGravity.GlobalNActiveParticles,
              ((double)TimeBinsGravity.GlobalNActiveParticles) / NTask, All.TotNumOfForces);
      fprintf(FdTimings, "  (direct) part/sec:  %g   ia/sec: %g\n", TimeBinsGravity.GlobalNActiveParticles / (sumt + 1.0e-20),
              TimeBinsGravity.GlobalNActiveParticles / (sumt + 1.0e-20) * TimeBinsGravity.GlobalNActiveParticles);
      myflush(FdTimings);
    }

  TIMER_STOP(CPU_TREEDIRECT);
}

#endif /* #ifdef ALLOW_DIRECT_SUMMATION */
