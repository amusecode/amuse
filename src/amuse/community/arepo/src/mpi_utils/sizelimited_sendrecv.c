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
 * \file        src/mpi_utils/sizelimited_sendrecv.c
 * \date        05/2018
 * \brief       MPI_Sendrecv operations split into chunks of maximum size.
 * \details     If the number of elements in the MPI_Sendrecv is larger than
 *              count_limit, the function will split up the communication into
 *              multiple chunks communicated by the usual MPI_Sendrecv routine.
 *              contains functions:
 *                int myMPI_Sendrecv(void *sendb, size_t sendcount,
 *                  MPI_Datatype sendtype, int dest, int sendtag, void *recvb,
 *                  size_t recvcount, MPI_Datatype recvtype, int source,
 *                  int recvtag, MPI_Comm comm, MPI_Status * status)
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

/*! \brief Self-made sendrecv function with limiter to the number of elements
 *         that can be sent in one go.
 *
 *  If the total message is longer, multiple MPI_Sendrecv calls are executed
 *  until the entire message has been communicated.
 *
 *  \param[in] sendb Initial address of send buffer.
 *  \param[in] sendcount Number of elements in send buffer.
 *  \param[in] sendtype Type of elements in send buffer (handle).
 *  \param[in] dest Rank of destination.
 *  \param[in] sendtag Send tag.
 *  \param[out] recvb Initial address of receive buffer.
 *  \param[in] recvcount Number of elements in receive buffer.
 *  \param[in] recvtype Type of elements in receive buffer (handle).
 *  \param[in] source Rank of source.
 *  \param[in] recvtag Receive tag.
 *  \param[in] comm MPI communicator.
 *  \param[out] status Status, referring to receive operation.
 *
 *  \return 0
 */
int myMPI_Sendrecv(void *sendb, size_t sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvb, size_t recvcount,
                   MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status)
{
  int iter      = 0, size_sendtype, size_recvtype, send_now, recv_now;
  char *sendbuf = (char *)sendb;
  char *recvbuf = (char *)recvb;

  if(dest != source)
    terminate("dest != source");

  MPI_Type_size(sendtype, &size_sendtype);
  MPI_Type_size(recvtype, &size_recvtype);

  if(dest == ThisTask)
    {
      memcpy(recvbuf, sendbuf, recvcount * size_recvtype);
      return 0;
    }

  size_t count_limit = MPI_MESSAGE_SIZELIMIT_IN_BYTES / size_sendtype;

  while(sendcount > 0 || recvcount > 0)
    {
      if(sendcount > count_limit)
        {
          send_now = count_limit;
          iter++;
        }
      else
        send_now = sendcount;

      if(recvcount > count_limit)
        recv_now = count_limit;
      else
        recv_now = recvcount;

      MPI_Sendrecv(sendbuf, send_now, sendtype, dest, sendtag, recvbuf, recv_now, recvtype, source, recvtag, comm, status);

      sendcount -= send_now;
      recvcount -= recv_now;

      sendbuf += send_now * size_sendtype;
      recvbuf += recv_now * size_recvtype;
    }

  return 0;
}
