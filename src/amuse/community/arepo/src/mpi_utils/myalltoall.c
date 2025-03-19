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
 * \file        src/mpi_utils/myalltoall.c
 * \date        05/2018
 * \brief       Specialized all-to-all MPI communication functions.
 * \details     contains functions:
 *                void myMPI_Alltoallv(void *sendb, size_t * sendcounts,
 *                  size_t * sdispls, void *recvb, size_t * recvcounts,
 *                  size_t * rdispls, int len, int big_flag, MPI_Comm comm)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 24.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/*! \brief A wrapper around MPI_Alltoallv that can deal with data in
 *         individual sends that are very big.
 *
 *  \param[in] sendb Starting address of send buffer.
 *  \param[in] sendcounts Integer array equal to the group size specifying the
 *             number of elements to send to each processor.
 *  \param[in] sdispls Integer array (of length group size). Entry j specifies
 *             the displacement (relative to sendbuf) from which to take the
 *             outgoing data destined for process j.
 *  \param[out] recvb Starting address of receive buffer.
 *  \param[in] recvcounts Integer array equal to the group size specifying the
 *             maximum number of elements that can be received from each
 *             processor.
 *  \param[in] rdispls Integer array (of length group size). Entry i specifies
 *             the displacement (relative to recvbuf at which to place the
 *             incoming data from process i.
 *  \param[in] len Size of single element in send array.
 *  \param[in] big_flag Flag if cummunication of large data. If not, the normal
 *             MPI_Alltoallv function is used.
 *  \param[in] comm MPI communicator.
 *
 *  \return void
 */
void myMPI_Alltoallv(void *sendb, size_t *sendcounts, size_t *sdispls, void *recvb, size_t *recvcounts, size_t *rdispls, int len,
                     int big_flag, MPI_Comm comm)
{
  char *sendbuf = (char *)sendb;
  char *recvbuf = (char *)recvb;

  if(big_flag == 0)
    {
      int ntask;
      MPI_Comm_size(comm, &ntask);

      int *scount = (int *)mymalloc("scount", ntask * sizeof(int));
      int *rcount = (int *)mymalloc("rcount", ntask * sizeof(int));
      int *soff   = (int *)mymalloc("soff", ntask * sizeof(int));
      int *roff   = (int *)mymalloc("roff", ntask * sizeof(int));

      for(int i = 0; i < ntask; i++)
        {
          scount[i] = sendcounts[i] * len;
          rcount[i] = recvcounts[i] * len;
          soff[i]   = sdispls[i] * len;
          roff[i]   = rdispls[i] * len;
        }

      MPI_Alltoallv(sendbuf, scount, soff, MPI_BYTE, recvbuf, rcount, roff, MPI_BYTE, comm);

      myfree(roff);
      myfree(soff);
      myfree(rcount);
      myfree(scount);
    }
  else
    {
      /* here we definitely have some large messages. We default to the
       * pair-wise protocoll, which should be most robust anyway.
       */

      int ntask, thistask;
      MPI_Comm_size(comm, &ntask);
      MPI_Comm_rank(comm, &thistask);

      for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
        {
          int target = thistask ^ ngrp;

          if(target < ntask)
            {
              if(sendcounts[target] > 0 || recvcounts[target] > 0)
                myMPI_Sendrecv(sendbuf + sdispls[target] * len, sendcounts[target] * len, MPI_BYTE, target, TAG_PDATA + ngrp,
                               recvbuf + rdispls[target] * len, recvcounts[target] * len, MPI_BYTE, target, TAG_PDATA + ngrp, comm,
                               MPI_STATUS_IGNORE);
            }
        }
    }
}
