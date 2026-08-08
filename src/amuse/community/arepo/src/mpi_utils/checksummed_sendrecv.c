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
 * \file        src/mpi_utils/checksummed_sendrecv.c
 * \date        05/2018
 * \brief       MPI send-receive communication with checksum to verify
 *              communication.
 * \details     contains functions:
 *                int MPI_Check_Sendrecv(void *sendbuf, int sendcount,
 *                  MPI_Datatype sendtype, int dest, int sendtag,
 *                  void *recvbufreal, int recvcount, MPI_Datatype recvtype,
 *                  int source, int recvtag, MPI_Comm comm,
 *                  MPI_Status * status)
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

#ifdef MPISENDRECV_CHECKSUM

#undef MPI_Sendrecv

/*! \brief MPI_Sendrecv with built-in check if message arrived properly.
 *
 *  \param[in] sendbuf Initial address of send buffer.
 *  \param[in] sendcount Number of elements in send buffer.
 *  \param[in] sendtype Type of elements in send buffer.
 *  \param[in] dest Rank of destination.
 *  \param[in] sendtag Send tag.
 *  \param[out] recvbufreal Initial adress of receive buffer.
 *  \param[in] recvcount Number of elements in receive buffer.
 *  \param[in] recvtype Type of elements in receive buffer .
 *  \param[in] source Rank of source.
 *  \param[in] recvtag Receive tag.
 *  \param[in] comm Communicator
 *  \param[out] status Status object; this refers to receive operation.
 *
 *  \return 0
 */
int MPI_Check_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbufreal, int recvcount,
                       MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status)
{
  int checksumtag = 1000, errtag = 2000;
  int i, iter                    = 0, err_flag, err_flag_imported, size_sendtype, size_recvtype;
  long long sendCheckSum, recvCheckSum, importedCheckSum;
  unsigned char *p, *buf, *recvbuf;
  char msg[500];

  if(dest != source)
    terminate("destination task different from source task");

  MPI_Type_size(sendtype, &size_sendtype);
  MPI_Type_size(recvtype, &size_recvtype);

  if(dest == ThisTask)
    {
      memcpy(recvbufreal, sendbuf, recvcount * size_recvtype);
      return 0;
    }

  if(!(buf = mymalloc(recvcount * size_recvtype + 1024)))
    terminate("not enough memory to allocate the buffer buf");

  for(i = 0, p = buf; i < recvcount * size_recvtype + 1024; i++)
    *p++ = 255;

  recvbuf = buf + 512;

  MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status);

  for(i = 0, p = buf; i < 512; i++, p++)
    {
      if(*p != 255)
        {
          sprintf(msg, "MPI-ERROR: Task=%d/%s: Recv occured before recv buffer. message-size=%d from %d, i=%d c=%d\n", ThisTask,
                  getenv("HOST"), recvcount, dest, i, *p);
          terminate(msg);
        }
    }

  for(i = 0, p = recvbuf + recvcount * size_recvtype; i < 512; i++, p++)
    {
      if(*p != 255)
        {
          sprintf(msg, "MPI-ERROR: Task=%d/%s: Recv occured after recv buffer. message-size=%d from %d, i=%d c=%d\n", ThisTask,
                  getenv("HOST"), recvcount, dest, i, *p);
          terminate(msg);
        }
    }

  for(i = 0, p = sendbuf, sendCheckSum = 0; i < sendcount * size_sendtype; i++, p++)
    sendCheckSum += *p;

  importedCheckSum = 0;

  if(dest > ThisTask)
    {
      if(sendcount > 0)
        MPI_Ssend(&sendCheckSum, sizeof(sendCheckSum), MPI_BYTE, dest, checksumtag, MPI_COMM_WORLD);
      if(recvcount > 0)
        MPI_Recv(&importedCheckSum, sizeof(importedCheckSum), MPI_BYTE, dest, checksumtag, MPI_COMM_WORLD, status);
    }
  else
    {
      if(recvcount > 0)
        MPI_Recv(&importedCheckSum, sizeof(importedCheckSum), MPI_BYTE, dest, checksumtag, MPI_COMM_WORLD, status);
      if(sendcount > 0)
        MPI_Ssend(&sendCheckSum, sizeof(sendCheckSum), MPI_BYTE, dest, checksumtag, MPI_COMM_WORLD);
    }

  checksumtag++;

  for(i = 0, p = recvbuf, recvCheckSum = 0; i < recvcount * size_recvtype; i++, p++)
    recvCheckSum += *p;

  err_flag = err_flag_imported = 0;

  if(recvCheckSum != importedCheckSum)
    {
      printf(
          "MPI-ERROR: Receive error on task=%d/%s from task=%d, message size=%d, sendcount=%d checksums= %d %d  %d %d. Try to fix "
          "it...\n",
          ThisTask, getenv("HOST"), source, recvcount, sendcount, (int)(recvCheckSum >> 32), (int)recvCheckSum,
          (int)(importedCheckSum >> 32), (int)importedCheckSum);
      myflush(stdout);

      err_flag = 1;
    }

  if(dest > ThisTask)
    {
      MPI_Ssend(&err_flag, 1, MPI_INT, dest, errtag, MPI_COMM_WORLD);
      MPI_Recv(&err_flag_imported, 1, MPI_INT, dest, errtag, MPI_COMM_WORLD, status);
    }
  else
    {
      MPI_Recv(&err_flag_imported, 1, MPI_INT, dest, errtag, MPI_COMM_WORLD, status);
      MPI_Ssend(&err_flag, 1, MPI_INT, dest, errtag, MPI_COMM_WORLD);
    }
  errtag++;

  if(err_flag > 0 || err_flag_imported > 0)
    {
      printf("Task=%d is on %s, wants to send %d and has checksum=%d %d of send data\n", ThisTask, getenv("HOST"), sendcount,
             (int)(sendCheckSum >> 32), (int)sendCheckSum);
      myflush(stdout);

      do
        {
          sendtag++;
          recvtag++;

          for(i = 0, p = recvbuf; i < recvcount * size_recvtype; i++, p++)
            *p = 0;

          if((iter & 1) == 0)
            {
              if(dest > ThisTask)
                {
                  if(sendcount > 0)
                    MPI_Ssend(sendbuf, sendcount, sendtype, dest, sendtag, MPI_COMM_WORLD);
                  if(recvcount > 0)
                    MPI_Recv(recvbuf, recvcount, recvtype, dest, recvtag, MPI_COMM_WORLD, status);
                }
              else
                {
                  if(recvcount > 0)
                    MPI_Recv(recvbuf, recvcount, recvtype, dest, recvtag, MPI_COMM_WORLD, status);
                  if(sendcount > 0)
                    MPI_Ssend(sendbuf, sendcount, sendtype, dest, sendtag, MPI_COMM_WORLD);
                }
            }
          else
            {
              if(iter > 5)
                {
                  printf("we're trying to send each byte now on task=%d (iter=%d)\n", ThisTask, iter);
                  myflush(stdout);
                  if(dest > ThisTask)
                    {
                      for(i = 0, p = sendbuf; i < sendcount * size_sendtype; i++, p++)
                        MPI_Ssend(p, 1, MPI_BYTE, dest, i, MPI_COMM_WORLD);
                      for(i = 0, p = recvbuf; i < recvcount * size_recvtype; i++, p++)
                        MPI_Recv(p, 1, MPI_BYTE, dest, i, MPI_COMM_WORLD, status);
                    }
                  else
                    {
                      for(i = 0, p = recvbuf; i < recvcount * size_recvtype; i++, p++)
                        MPI_Recv(p, 1, MPI_BYTE, dest, i, MPI_COMM_WORLD, status);
                      for(i = 0, p = sendbuf; i < sendcount * size_sendtype; i++, p++)
                        MPI_Ssend(p, 1, MPI_BYTE, dest, i, MPI_COMM_WORLD);
                    }
                }
              else
                {
                  MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm,
                               status);
                }
            }

          importedCheckSum = 0;

          for(i = 0, p = sendbuf, sendCheckSum = 0; i < sendcount * size_sendtype; i++, p++)
            sendCheckSum += *p;

          printf("Task=%d gas send_checksum=%d %d\n", ThisTask, (int)(sendCheckSum >> 32), (int)sendCheckSum);
          myflush(stdout);

          if(dest > ThisTask)
            {
              if(sendcount > 0)
                MPI_Ssend(&sendCheckSum, sizeof(sendCheckSum), MPI_BYTE, dest, checksumtag, MPI_COMM_WORLD);
              if(recvcount > 0)
                MPI_Recv(&importedCheckSum, sizeof(importedCheckSum), MPI_BYTE, dest, checksumtag, MPI_COMM_WORLD, status);
            }
          else
            {
              if(recvcount > 0)
                MPI_Recv(&importedCheckSum, sizeof(importedCheckSum), MPI_BYTE, dest, checksumtag, MPI_COMM_WORLD, status);
              if(sendcount > 0)
                MPI_Ssend(&sendCheckSum, sizeof(sendCheckSum), MPI_BYTE, dest, checksumtag, MPI_COMM_WORLD);
            }

          for(i = 0, p = recvbuf, recvCheckSum = 0; i < recvcount; i++, p++)
            recvCheckSum += *p;

          err_flag = err_flag_imported = 0;

          if(recvCheckSum != importedCheckSum)
            {
              printf(
                  "MPI-ERROR: Again (iter=%d) a receive error on task=%d/%s from task=%d, message size=%d, checksums= %d %d  %d %d. "
                  "Try to fix it...\n",
                  iter, ThisTask, getenv("HOST"), source, recvcount, (int)(recvCheckSum >> 32), (int)recvCheckSum,
                  (int)(importedCheckSum >> 32), (int)importedCheckSum);
              myflush(stdout);
              err_flag = 1;
            }

          if(dest > ThisTask)
            {
              MPI_Ssend(&err_flag, 1, MPI_INT, dest, errtag, MPI_COMM_WORLD);
              MPI_Recv(&err_flag_imported, 1, MPI_INT, dest, errtag, MPI_COMM_WORLD, status);
            }
          else
            {
              MPI_Recv(&err_flag_imported, 1, MPI_INT, dest, errtag, MPI_COMM_WORLD, status);
              MPI_Ssend(&err_flag, 1, MPI_INT, dest, errtag, MPI_COMM_WORLD);
            }

          if(err_flag == 0 && err_flag_imported == 0)
            break;

          errtag++;
          checksumtag++;
          iter++;
        }
      while(iter < 10);

      if(iter >= 10)
        {
          char buf[1000];
          int length;
          FILE *fd;

          sprintf(buf, "send_data_%d.dat", ThisTask);
          fd     = fopen(buf, "w");
          length = sendcount * size_sendtype;
          fwrite(&length, 1, sizeof(int), fd);
          fwrite(sendbuf, sendcount, size_sendtype, fd);
          fclose(fd);

          sprintf(buf, "recv_data_%d.dat", ThisTask);
          fd     = fopen(buf, "w");
          length = recvcount * size_recvtype;
          fwrite(&length, 1, sizeof(int), fd);
          fwrite(recvbuf, recvcount, size_recvtype, fd);
          fclose(fd);

          sprintf(msg, "MPI-ERROR: Even 10 trials proved to be insufficient on task=%d/%s. Stopping\n", ThisTask, getenv("HOST"));
          terminate(msg);
        }
    }

  memcpy(recvbufreal, recvbuf, recvcount * size_recvtype);

  myfree(buf);

  return 0;
}

#endif /* #ifdef MPISENDRECV_CHECKSUM */
