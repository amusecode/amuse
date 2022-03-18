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
 * \file        src/mesh/voronoi/voronoi_exchange.c
 * \date        05/2018
 * \brief       Algorithms that handle communication of Voronoi mesh data
 *              between MPI tasks.
 * \details     contains functions:
 *                void mesh_setup_exchange(void)
 *                void exchange_primitive_variables(void)
 *                void exchange_primitive_variables_and_gradients(void)
 *                int compare_primexch(const void *a, const void *b)
 *                void voronoi_update_ghost_velvertex(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 22.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../main/allvars.h"
#include "../../main/proto.h"

#include "voronoi.h"

/*! \brief Auxiliary data structure for communication of primitive variables.
 *
 */
struct data_primexch_compare
{
  int rank, task, index;
} * SortPrimExch, *SortPrimExch2;

/*! \brief Prepares exchange of primitive variables.
 *
 *  \return void
 */
void mesh_setup_exchange(void)
{
  if(All.TotNumGas == 0)
    return;

  TIMER_START(CPU_MESH_EXCHANGE);

  int listp;
  struct indexexch
  {
    int task, index;
  } * tmpIndexExch, *IndexExch;
  int i, j, p, task, off, count;
  int ngrp, recvTask, place;

  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if(ListExports[listp].origin != ThisTask)
            {
              Mesh_Send_count[ListExports[listp].origin]++;
            }
          listp = ListExports[listp].nextexport;
        }
    }

  MPI_Alltoall(Mesh_Send_count, 1, MPI_INT, Mesh_Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, Mesh_nimport = 0, Mesh_nexport = 0, Mesh_Recv_offset[0] = 0, Mesh_Send_offset[0] = 0; j < NTask; j++)
    {
      Mesh_nimport += Mesh_Recv_count[j];
      Mesh_nexport += Mesh_Send_count[j];

      if(j > 0)
        {
          Mesh_Send_offset[j] = Mesh_Send_offset[j - 1] + Mesh_Send_count[j - 1];
          Mesh_Recv_offset[j] = Mesh_Recv_offset[j - 1] + Mesh_Recv_count[j - 1];
        }
    }

  IndexExch    = (struct indexexch *)mymalloc("IndexExch", Mesh_nimport * sizeof(struct indexexch));
  tmpIndexExch = (struct indexexch *)mymalloc("tmpIndexExch", Mesh_nexport * sizeof(struct indexexch));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if((task = ListExports[listp].origin) != ThisTask)
            {
              place = ListExports[listp].index;
              off   = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

              tmpIndexExch[off].task  = ThisTask;
              tmpIndexExch[off].index = place;
            }
          listp = ListExports[listp].nextexport;
        }
    }

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&tmpIndexExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] * sizeof(struct indexexch), MPI_BYTE,
                           recvTask, TAG_DENS_A, &IndexExch[Mesh_Recv_offset[recvTask]],
                           Mesh_Recv_count[recvTask] * sizeof(struct indexexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpIndexExch);

  /* now we need to associate the imported data with the points stored in the DP[] array */

  SortPrimExch = (struct data_primexch_compare *)mymalloc("SortPrimExch", Mesh_nimport * sizeof(struct data_primexch_compare));

  for(i = 0; i < Mesh_nimport; i++)
    {
      SortPrimExch[i].rank  = i;
      SortPrimExch[i].task  = IndexExch[i].task;
      SortPrimExch[i].index = IndexExch[i].index;
    }

  /* let sort the data according to task and index */
  mysort(SortPrimExch, Mesh_nimport, sizeof(struct data_primexch_compare), compare_primexch);

  SortPrimExch2 = (struct data_primexch_compare *)mymalloc("SortPrimExch2", Mesh.Ndp * sizeof(struct data_primexch_compare));

  for(i = 0, count = 0; i < Mesh.Ndp; i++)
    {
      if(Mesh.DP[i].task != ThisTask)
        {
          SortPrimExch2[count].rank  = i;
          SortPrimExch2[count].task  = Mesh.DP[i].task;
          SortPrimExch2[count].index = Mesh.DP[i].index;
          count++;
        }
    }

  /* let sort according to task and index */
  mysort(SortPrimExch2, count, sizeof(struct data_primexch_compare), compare_primexch);

  /* count can be larger than nimport because a foreigh particle can appear
     multiple times on the local domain, due to periodicity */

  for(i = 0, j = 0; i < count; i++)
    {
      if(SortPrimExch2[i].task != SortPrimExch[j].task || SortPrimExch2[i].index != SortPrimExch[j].index)
        j++;

      if(j >= Mesh_nimport)
        terminate("j >= Mesh_nimport");

      Mesh.DP[SortPrimExch2[i].rank].index =
          SortPrimExch[j].rank; /* note: this change is now permanent and available for next exchange */
    }

  myfree(SortPrimExch2);
  myfree(SortPrimExch);
  myfree(IndexExch);

  /* allocate structures needed to exchange the actual information for ghost cells */
  PrimExch = (struct primexch *)mymalloc_movable(&PrimExch, "PrimExch", Mesh_nimport * sizeof(struct primexch));
  GradExch = (struct grad_data *)mymalloc_movable(&GradExch, "GradExch", Mesh_nimport * sizeof(struct grad_data));

  TIMER_STOP(CPU_MESH_EXCHANGE);
}

/*! \brief Communicate primitive variables across MPI tasks.
 *
 *  This routine is called before gradient calculation, afterwards,
 *  exchange_primitive_variables_and_gradients is called.
 *
 *  \return void
 */
void exchange_primitive_variables(void)
{
  if(All.TotNumGas == 0)
    return;

  TIMER_START(CPU_MESH_EXCHANGE);

  int listp;
  struct primexch *tmpPrimExch;
  int i, j, p, task, off;
  int ngrp, recvTask, place;

  tmpPrimExch = (struct primexch *)mymalloc("tmpPrimExch", Mesh_nexport * sizeof(struct primexch));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if((task = ListExports[listp].origin) != ThisTask)
            {
              place = ListExports[listp].index;
              off   = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

              tmpPrimExch[off].Volume = SphP[place].Volume;

              tmpPrimExch[off].Density = SphP[place].Density;

              tmpPrimExch[off].Pressure = SphP[place].Pressure;

#ifdef MHD
              tmpPrimExch[off].B[0] = SphP[place].B[0];
              tmpPrimExch[off].B[1] = SphP[place].B[1];
              tmpPrimExch[off].B[2] = SphP[place].B[2];
#ifdef MHD_POWELL
              tmpPrimExch[off].DivB = SphP[place].DivB;
#endif /* #ifdef MHD_POWELL */
#endif /* #ifdef MHD */

              tmpPrimExch[off].OldMass      = SphP[place].OldMass;
              tmpPrimExch[off].SurfaceArea  = SphP[place].SurfaceArea;
              tmpPrimExch[off].ActiveArea   = SphP[place].ActiveArea;
              tmpPrimExch[off].TimeBinHydro = P[place].TimeBinHydro;

#ifdef MAXSCALARS
              for(j = 0; j < N_Scalar; j++)
                tmpPrimExch[off].Scalars[j] = *(MyFloat *)(((char *)(&SphP[place])) + scalar_elements[j].offset);
#endif /* #ifdef MAXSCALARS */

              tmpPrimExch[off].TimeLastPrimUpdate = SphP[place].TimeLastPrimUpdate;

              for(j = 0; j < 3; j++)
                {
                  tmpPrimExch[off].VelGas[j] = P[place].Vel[j];
                  tmpPrimExch[off].Center[j] = SphP[place].Center[j];
                }
              tmpPrimExch[off].Csnd = get_sound_speed(place);
            }
          listp = ListExports[listp].nextexport;
        }
    }

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&tmpPrimExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] * sizeof(struct primexch), MPI_BYTE,
                           recvTask, TAG_DENS_A, &PrimExch[Mesh_Recv_offset[recvTask]],
                           Mesh_Recv_count[recvTask] * sizeof(struct primexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpPrimExch);

  TIMER_STOP(CPU_MESH_EXCHANGE);
}

/*! \brief Communicate primitive variables and gradients across MPI tasks.
 *
 *  This routine is called after gradient calculation.
 *
 *  \return void
 */
void exchange_primitive_variables_and_gradients(void)
{
  if(All.TotNumGas == 0)
    return;

  TIMER_START(CPU_MESH_EXCHANGE);

  int listp;
  struct grad_data *tmpGradExch;
  struct primexch *tmpPrimExch;

  int i, j, p, task, off;
  int ngrp, recvTask, place;

  tmpPrimExch = (struct primexch *)mymalloc("tmpPrimExch", Mesh_nexport * sizeof(struct primexch));
  tmpGradExch = (struct grad_data *)mymalloc("tmpGradExch", Mesh_nexport * sizeof(struct grad_data));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      /* in case previous steps already lowered the Mass, update OldMass to yield together with metallicity vector conservative
       * estimate of metal mass of each species contained in cell */
      if(P[p].Mass < SphP[p].OldMass)
        SphP[p].OldMass = P[p].Mass;

      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if((task = ListExports[listp].origin) != ThisTask)
            {
              place = ListExports[listp].index;
              off   = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

              tmpPrimExch[off].Volume   = SphP[place].Volume;
              tmpPrimExch[off].Density  = SphP[place].Density;
              tmpPrimExch[off].Pressure = SphP[place].Pressure;

#ifdef MHD
              tmpPrimExch[off].B[0] = SphP[place].B[0];
              tmpPrimExch[off].B[1] = SphP[place].B[1];
              tmpPrimExch[off].B[2] = SphP[place].B[2];
#ifdef MHD_POWELL
              tmpPrimExch[off].DivB = SphP[place].DivB;
#endif /* #ifdef MHD_POWELL */
#endif /* #ifdef MHD */

              tmpPrimExch[off].OldMass     = SphP[place].OldMass;
              tmpPrimExch[off].SurfaceArea = SphP[place].SurfaceArea;
              tmpPrimExch[off].ActiveArea  = SphP[place].ActiveArea;

              tmpPrimExch[off].TimeBinHydro = P[place].TimeBinHydro;

#ifdef MAXSCALARS
              for(j = 0; j < N_Scalar; j++)
                tmpPrimExch[off].Scalars[j] = *(MyFloat *)(((char *)(&SphP[place])) + scalar_elements[j].offset);
#endif /* #ifdef MAXSCALARS */

              tmpPrimExch[off].TimeLastPrimUpdate = SphP[place].TimeLastPrimUpdate;

              for(j = 0; j < 3; j++)
                {
                  tmpPrimExch[off].VelGas[j]    = P[place].Vel[j];
                  tmpPrimExch[off].Center[j]    = SphP[place].Center[j];
                  tmpPrimExch[off].VelVertex[j] = SphP[place].VelVertex[j];
                }

              tmpGradExch[off] = SphP[place].Grad;

              tmpPrimExch[off].Csnd = get_sound_speed(place);
            }
          listp = ListExports[listp].nextexport;
        }
    }

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* exchange the data */
              MPI_Sendrecv(&tmpPrimExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] * sizeof(struct primexch), MPI_BYTE,
                           recvTask, TAG_DENS_A, &PrimExch[Mesh_Recv_offset[recvTask]],
                           Mesh_Recv_count[recvTask] * sizeof(struct primexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);

              MPI_Sendrecv(&tmpGradExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] * sizeof(struct grad_data), MPI_BYTE,
                           recvTask, TAG_HYDRO_A, &GradExch[Mesh_Recv_offset[recvTask]],
                           Mesh_Recv_count[recvTask] * sizeof(struct grad_data), MPI_BYTE, recvTask, TAG_HYDRO_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpGradExch);
  myfree(tmpPrimExch);

  TIMER_STOP(CPU_MESH_EXCHANGE);

  /* note: because the sequence is the same as before, we don't have to do the sorts again */
}

/*! \brief Compare two data primexch compare objects.
 *
 *  The following variables (most important first):
 *      task
 *      index
 *
 *  \param[in] a Pointer to first data primexch compare object.
 *  \param[in] b Pointer to second data primexch compare object.
 *
 *  \return (-1,0,1); -1 if a < b.
 */
int compare_primexch(const void *a, const void *b)
{
  if(((struct data_primexch_compare *)a)->task < ((struct data_primexch_compare *)b)->task)
    return -1;

  if(((struct data_primexch_compare *)a)->task > ((struct data_primexch_compare *)b)->task)
    return +1;

  if(((struct data_primexch_compare *)a)->index < ((struct data_primexch_compare *)b)->index)
    return -1;

  if(((struct data_primexch_compare *)a)->index > ((struct data_primexch_compare *)b)->index)
    return +1;

  return 0;
}

/*! \brief Communicates vertex velocity divergence data across MPI tasks.
 *
 *  \return 0
 */
#ifdef OUTPUT_VERTEX_VELOCITY_DIVERGENCE
void voronoi_update_ghost_velvertex(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  int listp;
  int i, j, p, task, off;
  int ngrp, recvTask, place;
  struct velvertex_data
  {
    MyFloat VelVertex[3];
  } * tmpVelVertexExch, *tmpVelVertexRecv;

  tmpVelVertexExch = (struct velvertex_data *)mymalloc("tmpVelVertexExch", Mesh_nexport * sizeof(struct velvertex_data));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(i = 0; i < NumGasInMesh; i++)
    {
      p = List_InMesh[i];

      listp = List_P[p].firstexport;
      while(listp >= 0)
        {
          if((task = ListExports[listp].origin) != ThisTask)
            {
              place = ListExports[listp].index;
              off   = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

              for(j = 0; j < 3; j++)
                {
                  tmpVelVertexExch[off].VelVertex[j] = SphP[place].VelVertex[j];
                }
            }
          listp = ListExports[listp].nextexport;
        }
    }

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              tmpVelVertexRecv =
                  (struct velvertex_data *)mymalloc("tmpVelVertexRecv", Mesh_Recv_count[recvTask] * sizeof(struct velvertex_data));

              /* get the values */
              MPI_Sendrecv(&tmpVelVertexExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask] * sizeof(struct velvertex_data),
                           MPI_BYTE, recvTask, TAG_DENS_A, tmpVelVertexRecv, Mesh_Recv_count[recvTask] * sizeof(struct velvertex_data),
                           MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              for(i = 0; i < Mesh_Recv_count[recvTask]; i++)
                {
                  for(j = 0; j < 3; j++)
                    {
                      PrimExch[Mesh_Recv_offset[recvTask] + i].VelVertex[j] = tmpVelVertexExch[i].VelVertex[j];
                    }
                }

              myfree(tmpVelVertexRecv);
            }
        }
    }

  myfree(tmpVelVertexExch);

  CPU_Step[CPU_SET_VERTEXVELS] += measure_time();
}
#endif /* #ifdef OUTPUT_VERTEX_VELOCITY_DIVERGENCE */
