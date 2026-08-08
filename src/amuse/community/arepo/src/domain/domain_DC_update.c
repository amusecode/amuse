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
 * \file        src/domain_DC_update.c
 * \date        05/2018
 * \brief       Algorithms for voronoi dynamic update
 * \details     contains functions:
 *                void domain_mark_in_trans_table(int i, int task)
 *                void domain_exchange_and_update_DC(void)
 *                int domain_compare_connection_ID(const void *a,
 *                  const void *b)
 *                int domain_compare_local_trans_data_ID(const void *a,
 *                  const void *b)
 *                int domain_compare_recv_trans_data_ID(const void *a,
 *                  const void *b)
 *                int domain_compare_recv_trans_data_oldtask(const void *a,
 *                  const void *b)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 17.05.2018 Prepared file for public release -- Rainer Weinberger
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

struct trans_data *trans_table;
int N_trans;

/*! \brief Data structure for local auxiliary translation table.
 */
static struct local_aux_trans_data
{
  MyIDType ID;
  int new_index;
} * local_trans_data;

/*! \brief Data structure for communicating the translation table.
 */
static struct aux_trans_data
{
  MyIDType ID;
  int old_task;
  int old_index;
  int new_index;
} * send_trans_data, *recv_trans_data;

/*! \brief Data structure for transcribing data.
 */
static struct aux_transscribe_data
{
  int old_index;
  int new_task;
  int new_index;
  int image_flags;
} * send_transscribe_data, *recv_transscribe_data;

/*! \brief Fill translation table.
 *
 *  Mark where cells are moved to and mark in DC accordingly to make sure
 *  they get communicated to the same task.
 *
 *  \param[in] i Index in P and SphP arrays.
 *  \param[in] task Task to which particle i is exported.
 *
 *  \return void
 */
void domain_mark_in_trans_table(int i, int task)
{
  if(Largest_Nvc > 0)
    {
      if(i < NumGas)
        {
          trans_table[i].ID       = P[i].ID;
          trans_table[i].new_task = task;

          int q = SphP[i].first_connection;

          while(q >= 0)
            {
              int qq = DC[q].next;
              if(q == qq)
                terminate("preventing getting stuck in a loop due to q == DC[q].next : i=%d q=%d last_connection=%d", i, q,
                          SphP[i].last_connection);

              if((P[i].Mass == 0 && P[i].ID == 0) || P[i].Type != 0) /* this cell has been deleted or turned into a star */
                DC[q].next = -1;
              else
                DC[q].next = task; /* we will temporarily use the next variable to store the new task */

              if(q == SphP[i].last_connection)
                break;

              q = qq;
            }
        }
      else if(i < N_trans)
        trans_table[i].new_task = -1; /* this one has been removed by rerrange_particle_sequence() */
    }
}

/*! \brief Communicates connections.
 *
 *  This algorithms communicates Delauny connections and updates them on the
 *  new task.
 *
 *  \return void
 */
void domain_exchange_and_update_DC(void)
{
  double t0 = second();

#if !defined(GRAVITY_NOT_PERIODIC) && !defined(DO_NOT_RANDOMIZE_DOMAINCENTER) && defined(SELFGRAVITY)
  /* remove all image flags, after our box movement stunt they are all incorrect anyway */
  for(int i = 0; i < MaxNvc; i++)
    {
      DC[i].image_flags = 1;
    }
#endif /* #if !defined(GRAVITY_NOT_PERIODIC) && !defined(DO_NOT_RANDOMIZE_DOMAINCENTER) && defined(SELFGRAVITY) */

  /* first, we need to complete the translation table */
  for(int j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(int i = 0; i < N_trans; i++)
    if(trans_table[i].new_task >= 0)
      Send_count[trans_table[i].new_task]++;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  int nimport = 0, nexport = 0;
  Recv_offset[0] = Send_offset[0] = 0;

  for(int j = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  send_trans_data = mymalloc("send_trans_data", nexport * sizeof(struct aux_trans_data));
  recv_trans_data = mymalloc("recv_trans_data", nimport * sizeof(struct aux_trans_data));

  for(int j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(int i = 0; i < N_trans; i++)
    {
      int task = trans_table[i].new_task;
      if(task >= 0)
        {
          send_trans_data[Send_offset[task] + Send_count[task]].ID        = trans_table[i].ID;
          send_trans_data[Send_offset[task] + Send_count[task]].old_index = i;
          send_trans_data[Send_offset[task] + Send_count[task]].old_task  = ThisTask;
          Send_count[task]++;
        }
    }

  /* exchange the data */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&send_trans_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct aux_trans_data), MPI_BYTE,
                       recvTask, TAG_DENS_B, &recv_trans_data[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(struct aux_trans_data), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
    }

  /* let's now sort the incoming list according to ID */
  mysort(recv_trans_data, nimport, sizeof(struct aux_trans_data), domain_compare_recv_trans_data_ID);

  /* make an auxiliary list for the local particles that we will also sort according to ID */
  local_trans_data = mymalloc("local_trans_data", NumGas * sizeof(struct local_aux_trans_data));
  for(int i = 0; i < NumGas; i++)
    {
      local_trans_data[i].ID        = P[i].ID;
      local_trans_data[i].new_index = i;
    }
  mysort(local_trans_data, NumGas, sizeof(struct local_aux_trans_data), domain_compare_local_trans_data_ID);

  int i, j;
  /* now we go through and put in the new index for matching IDs */
  for(i = 0, j = 0; i < nimport && j < NumGas;)
    {
      if(recv_trans_data[i].ID < local_trans_data[j].ID)
        {
          recv_trans_data[i].new_index = -1; /* this particle has been eliminated */
          i++;
        }
      else if(recv_trans_data[i].ID > local_trans_data[j].ID)
        j++;
      else
        {
          recv_trans_data[i].new_index = local_trans_data[j].new_index;
          i++;
          j++;
        }
    }

  for(; i < nimport; i++)
    recv_trans_data[i].new_index = -1; /* this particle has been eliminated */

  myfree(local_trans_data);

  /* now order the received data by sending task, so that we can return it */
  mysort(recv_trans_data, nimport, sizeof(struct aux_trans_data), domain_compare_recv_trans_data_oldtask);

  /* return the data */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&recv_trans_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct aux_trans_data), MPI_BYTE,
                       recvTask, TAG_DENS_B, &send_trans_data[Send_offset[recvTask]],
                       Send_count[recvTask] * sizeof(struct aux_trans_data), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
    }

  /* now let's fill in the new_index entry into the translation table */
  for(int i = 0; i < nexport; i++)
    trans_table[send_trans_data[i].old_index].new_index = send_trans_data[i].new_index;

  myfree(recv_trans_data);
  myfree(send_trans_data);

  /* it's now time to transcribe the task and index fields in the DC list */
  for(int j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(int i = 0; i < MaxNvc; i++)
    {
      int task = DC[i].task;
      if(task >= 0)
        {
          if(task >= NTask)
            terminate("i=%d Nvc=%d MaxNvc=%d task=%d\n", i, Nvc, MaxNvc, task);

          Send_count[task]++;
        }
    }

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  nimport = nexport = 0;
  Recv_offset[0] = Send_offset[0] = 0;

  for(int j = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  send_transscribe_data = mymalloc("send_transscribe_data", nexport * sizeof(struct aux_transscribe_data));
  recv_transscribe_data = mymalloc("recv_transscribe_data", nimport * sizeof(struct aux_transscribe_data));

  for(int j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(int i = 0; i < MaxNvc; i++)
    {
      int task = DC[i].task;
      if(task >= 0)
        {
          send_transscribe_data[Send_offset[task] + Send_count[task]].old_index   = DC[i].index;
          send_transscribe_data[Send_offset[task] + Send_count[task]].image_flags = DC[i].image_flags;
          Send_count[task]++;
        }
    }

  /* exchange the data */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&send_transscribe_data[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct aux_transscribe_data),
                       MPI_BYTE, recvTask, TAG_DENS_B, &recv_transscribe_data[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(struct aux_transscribe_data), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
    }

  for(int i = 0; i < nimport; i++)
    {
      if(recv_transscribe_data[i].old_index >= N_trans)
        terminate("recv_transscribe_data[i].old_index >= N_trans");

      if(recv_transscribe_data[i].old_index < 0)
        terminate("recv_transscribe_data[i].old_index < 0");

      int old_index = recv_transscribe_data[i].old_index;

      recv_transscribe_data[i].new_task  = trans_table[old_index].new_task;
      recv_transscribe_data[i].new_index = trans_table[old_index].new_index;

#if !defined(GRAVITY_NOT_PERIODIC) && !defined(DO_NOT_RANDOMIZE_DOMAINCENTER) && defined(SELFGRAVITY)
      // Nothing to do here
#else  /* #if !defined(GRAVITY_NOT_PERIODIC) && !defined(DO_NOT_RANDOMIZE_DOMAINCENTER) && defined(SELFGRAVITY) */
      if(recv_transscribe_data[i].new_task >= 0)
        {
          if(trans_table[old_index].wrapped)
            {
              int bitflags = ffs(recv_transscribe_data[i].image_flags) - 1;
              int zbits    = (bitflags / 9);
              int ybits    = (bitflags - zbits * 9) / 3;
              int xbits    = bitflags - zbits * 9 - ybits * 3;

              if(trans_table[old_index].wrapped & 1)
                {
                  if(xbits == 1)
                    xbits = 0;
                  else if(xbits == 0)
                    xbits = 2;
                  else /* xbits == 2 */
                    terminate("b");
                }
              else if(trans_table[old_index].wrapped & 2)
                {
                  if(xbits == 1)
                    {
                      terminate("a");
                    }
                  else if(xbits == 0)
                    xbits = 1;
                  else /* xbits == 2 */
                    xbits = 0;
                }

              if(trans_table[old_index].wrapped & 4)
                {
                  if(ybits == 1)
                    ybits = 0;
                  else if(ybits == 0)
                    ybits = 2;
                  else
                    {
                      terminate("b");
                    }
                }
              else if(trans_table[old_index].wrapped & 8)
                {
                  if(ybits == 1)
                    {
                      terminate("a");
                    }
                  else if(ybits == 0)
                    ybits = 1;
                  else
                    ybits = 0;
                }

              if(trans_table[old_index].wrapped & 16)
                {
                  if(zbits == 1)
                    zbits = 0;
                  else if(zbits == 0)
                    zbits = 2;
                  else
                    {
                      terminate("b");
                    }
                }
              else if(trans_table[old_index].wrapped & 32)
                {
                  if(zbits == 1)
                    {
                      terminate("a");
                    }
                  else if(zbits == 0)
                    zbits = 1;
                  else
                    zbits = 0;
                }

              recv_transscribe_data[i].image_flags = (1 << (zbits * 9 + ybits * 3 + xbits));
            }
        }
#endif /* #if !defined(GRAVITY_NOT_PERIODIC) && !defined(DO_NOT_RANDOMIZE_DOMAINCENTER) && defined(SELFGRAVITY) #else */
    }

  /* now return the data */
  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&recv_transscribe_data[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct aux_transscribe_data),
                       MPI_BYTE, recvTask, TAG_DENS_B, &send_transscribe_data[Send_offset[recvTask]],
                       Send_count[recvTask] * sizeof(struct aux_transscribe_data), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
    }

  for(int j = 0; j < NTask; j++)
    Send_count[j] = 0;

  /* copy the results over to the DC structure */
  for(int i = 0; i < MaxNvc; i++)
    {
      int task = DC[i].task;
      if(task >= 0)
        {
          DC[i].task        = send_transscribe_data[Send_offset[task] + Send_count[task]].new_task;
          DC[i].index       = send_transscribe_data[Send_offset[task] + Send_count[task]].new_index;
          DC[i].image_flags = send_transscribe_data[Send_offset[task] + Send_count[task]].image_flags;
          Send_count[task]++;
        }
    }

  myfree(recv_transscribe_data);
  myfree(send_transscribe_data);

  /* now we can exchange the DC data. The task where each item should go is stored in 'next' at this point */
  for(int j = 0; j < NTask; j++)
    Send_count[j] = 0;

  /* count where they should go */
  for(int i = 0; i < MaxNvc; i++)
    {
      if(DC[i].task >= 0)
        {
          int task = DC[i].next;
          if(task >= 0)
            {
              if(task >= NTask)
                terminate("Thistask=%d  i=%d Nvc=%d MaxNvc=%d DC[i].task=%d DC[i].next=%d\n", ThisTask, i, Nvc, MaxNvc, DC[i].task,
                          DC[i].next);

              if(DC[i].index >= 0)
                Send_count[task]++;
            }
        }
    }

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  nimport = nexport = 0;
  Recv_offset[0] = Send_offset[0] = 0;

  for(int j = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* make sure that we have enough room to store the new DC list */
  while(nimport > MaxNvc)
    {
      int old_MaxNvc = MaxNvc;
      Mesh.Indi.AllocFacNvc *= ALLOC_INCREASE_FACTOR;
      MaxNvc = Mesh.Indi.AllocFacNvc;
#ifdef VERBOSE
      printf("Task=%d: increase memory allocation, MaxNvc=%d Indi.AllocFacNvc=%g\n", ThisTask, MaxNvc, Mesh.Indi.AllocFacNvc);
#endif /* #ifdef VERBOSE */
      DC = myrealloc_movable(DC, MaxNvc * sizeof(connection));
      for(int n = old_MaxNvc; n < MaxNvc; n++)
        DC[n].task = -1;
    }

  connection *tmpDC = mymalloc("tmpDC", nexport * sizeof(connection));

  for(int j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(int i = 0; i < MaxNvc; i++)
    {
      if(DC[i].task >= 0)
        {
          int task = DC[i].next;

          if(task >= 0 && DC[i].index >= 0)
            tmpDC[Send_offset[task] + Send_count[task]++] = DC[i];
        }
    }

  /* exchange the connection information */

  for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&tmpDC[Send_offset[recvTask]], Send_count[recvTask] * sizeof(connection), MPI_BYTE, recvTask, TAG_DENS_B,
                       &DC[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(connection), MPI_BYTE, recvTask, TAG_DENS_B,
                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  myfree(tmpDC);

  Nvc = nimport;

  /* mark the remaining ones as available */
  for(int i = Nvc; i < MaxNvc - 1; i++)
    {
      DC[i].next = i + 1;
      DC[i].task = -1;
    }
  DC[MaxNvc - 1].next = -1;
  DC[MaxNvc - 1].task = -1;

  if(Nvc < MaxNvc)
    FirstUnusedConnection = Nvc;
  else
    FirstUnusedConnection = -1;

  /* now we need to connect the information to the particles, this we do via the IDs */

  local_trans_data = mymalloc("local_trans_data", NumGas * sizeof(struct local_aux_trans_data));
  for(int i = 0; i < NumGas; i++)
    {
      local_trans_data[i].ID        = P[i].ID;
      local_trans_data[i].new_index = i; /* is here used as rank of the particle */
    }
  mysort(local_trans_data, NumGas, sizeof(struct local_aux_trans_data), domain_compare_local_trans_data_ID);

  mysort(DC, Nvc, sizeof(connection), domain_compare_connection_ID);

  int last = -1;
  for(i = 0, j = 0; i < NumGas && j < Nvc; i++)
    {
      int k = local_trans_data[i].new_index;

      if(P[k].ID < DC[j].ID)
        {
          /* this particle has no connection information (new cell) */
          SphP[k].first_connection = -1;
          SphP[k].last_connection  = -1;
        }
      else if(P[k].ID == DC[j].ID)
        {
          SphP[k].first_connection = j;

          while(j < Nvc)
            {
              SphP[k].last_connection = j;

              if(last >= 0)
                DC[last].next = j;

              last = j;
              j++;
              if(j >= Nvc)
                break;
              if(P[k].ID != DC[j].ID)
                break;
            }
        }
      else
        {
          terminate("strange");
        }
    }

  for(; i < NumGas; i++)
    {
      int k                    = local_trans_data[i].new_index;
      SphP[k].first_connection = -1;
      SphP[k].last_connection  = -1;
    }

  if(last >= 0)
    DC[last].next = -1;

  myfree(local_trans_data);

  double t1 = second();
  mpi_printf("DOMAIN: done with rearranging connection information (took %g sec)\n", timediff(t0, t1));
}

/*! \brief Compare which ID is larger.
 *
 *  For connection data.
 *
 *  \param[in] a Pointer to first object.
 *  \param[in] b Pointer to second object.
 *
 *  \return (-1,0,1) -1 if a->ID is smaller.
 */
int domain_compare_connection_ID(const void *a, const void *b)
{
  if(((connection *)a)->ID < (((connection *)b)->ID))
    return -1;

  if(((connection *)a)->ID > (((connection *)b)->ID))
    return +1;

  return 0;
}

/*! \brief Compare which ID is larger.
 *
 *  For local_aux_trans_data.
 *
 *  \param[in] a Pointer to first object.
 *  \param[in] b Pointer to second object.
 *
 *  \return (-1,0,1) -1 if a->ID is smaller.
 */
int domain_compare_local_trans_data_ID(const void *a, const void *b)
{
  if(((struct local_aux_trans_data *)a)->ID < (((struct local_aux_trans_data *)b)->ID))
    return -1;

  if(((struct local_aux_trans_data *)a)->ID > (((struct local_aux_trans_data *)b)->ID))
    return +1;

  return 0;
}

/*! \brief Compare which ID is larger.
 *
 *  For aux_trans_data.
 *
 *  \param[in] a Pointer to first object.
 *  \param[in] b Pointer to second object.
 *
 *  \return (-1,0,1) -1 if a->ID is smaller.
 */
int domain_compare_recv_trans_data_ID(const void *a, const void *b)
{
  if(((struct aux_trans_data *)a)->ID < (((struct aux_trans_data *)b)->ID))
    return -1;

  if(((struct aux_trans_data *)a)->ID > (((struct aux_trans_data *)b)->ID))
    return +1;

  return 0;
}

/*! \brief Compare which old_task is larger.
 *
 *  For aux_trans_data.
 *
 *  \param[in] a Pointer to first object.
 *  \param[in] b Pointer to second object.
 *
 *  \return (-1,0,1) -1 if a->old_task is smaller.
 */
int domain_compare_recv_trans_data_oldtask(const void *a, const void *b)
{
  if(((struct aux_trans_data *)a)->old_task < (((struct aux_trans_data *)b)->old_task))
    return -1;

  if(((struct aux_trans_data *)a)->old_task > (((struct aux_trans_data *)b)->old_task))
    return +1;

  return 0;
}
