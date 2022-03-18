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
 * \file        src/mpi_utils/myIBarrier.c
 * \date        05/2018
 * \brief       Home-made MPI_Ibarrier routine.
 * \details     Non-blocking version of MPI_Barrier; Once reaching this point,
 *              a process notifies this to other tasks.
 *              contains functions:
 *                void myIBarrier(MPI_Comm comm, struct sMyIBarrier *barrier)
 *                void myIBarrierTest(struct sMyIBarrier *barrier, int *flag,
 *                  MPI_Status * unused)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 04.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifdef MYIBARRIER

#include <strings.h>

#include "myIBarrier.h"

/*! \brief Non-blocking MPI barrier; Notifies other tasks once it is called.
 *
 *  \param[in] comm MPI communicator.
 *  \param[in, out] Object containing information about the barrier.
 *
 *  \return void
 */
void myIBarrier(MPI_Comm comm, struct sMyIBarrier *barrier)
{
  barrier->comm = comm;
  MPI_Comm_rank(comm, &barrier->rank);
  MPI_Comm_size(comm, &barrier->nTasks);

  barrier->nLevels   = fls(barrier->rank - 1);
  barrier->LevelDone = mymalloc("myIBarrier", barrier->nLevels);
  memset(barrier->LevelDone, 0, barrier->nLevels);

  /* find messages we would expect from nonexisting tasks */
  for(level = 0; level < barrier->nLevels; level++)
    if((barrier->rank & (1 << level) == 0) && (barrier->rank + (1 << level) >= barrier->nTasks))
      barrier->LevelDone[level] = 1;

  /* find out if we have to send or wait */
  int level = 0;
  while(level < barrier->nLevels)
    {
      if(barrier->rank & (1 << level))
        {
          /* we need to send our result */
          int target = barrier->rank - (1 << level);
          int level  = barrier->nLevels;
          MPI_Isend(&level, 1, MPI_INT, target, MPI_TAG_IBARRIER, barrier->comm);
          break;
        }
      else
        {
          /* check if there is something to recieve in which case we have to wait, otherwise go down one level */
          if(barrier->rank + (1 << level) < barrier->nTasks)
            {
              barrier->levelDone[level] = 1;
              break;
            }
          else
            level++;
        }
    }
}

/*! \brief Test function for myIBarrier.
 *
 *  \param[in] barrier Object containing information about the barrier.
 *  \param[out] flag Was test successful?
 *  \param[in] unused Unused MPI_Status.
 *
 *  \return void
 */
void myIBarrierTest(struct sMyIBarrier *barrier, int *flag, MPI_Status *unused)
{
  flag = 0;

  int rflag;
  MPI_Status status;

  MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_IBARRIER, barrier->comm, &rflag, &status);

  if(rflag)
    {
      int source = status.MPI_SOURCE;

      int level;
      MPI_Recv(&level, 1, MPI_INT, source, MPI_TAG_IBARRIER, barrier->comm, MPI_STATUS_IGNORE);

      if(source > barrier->rank)
        {
          /* we got another result, so lets check if we can send out further */
          while((level < barrier->nLevels) && barrier->LevelDone[level])
            level++;

          if(level == barrier->nLevels)
            {
              if(barrier->rank != 0)
                terminate("fail");
              /* ok, the barrier resolved, tell everyone */

              for(level = 0; level < barrier->nLevels; level++)
                {
                  if(barrier->rank & (1 << level) == 0)
                    {
                      int target = barrier->rank + (1 << level);
                      if(target < barrier->nTasks)
                        MPI_Isend(&level, 1, MPI_INT, target, MPI_TAG_IBARRIER, barrier->comm);
                    }
                  else
                    break;
                }

              flag = 1;
            }
          else
            {
              if(barrier->rank & (1 << level))
                {
                  /* we need to send our result */
                  int target = barrier->rank - (1 << level);
                  int level  = barrier->nLevels;
                  MPI_Isend(&level, 1, MPI_INT, target, MPI_TAG_IBARRIER, barrier->comm);
                }
              else
                {
                  barrier->LevelDone[level] = 1;
                }
            }
        }
      else
        {
          for(; level < barrier->nLevels; level++)
            {
              if(barrier->rank & (1 << level) == 0)
                {
                  int target = barrier->rank + (1 << level);
                  if(target < barrier->nTasks)
                    MPI_Isend(&level, 1, MPI_INT, target, MPI_TAG_IBARRIER, barrier->comm);
                }
              else
                break;
            }

          flag = 1;
        }
    }
}

#endif /* #ifdef MYIBARRIER */
