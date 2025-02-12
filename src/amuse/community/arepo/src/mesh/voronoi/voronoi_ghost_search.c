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
 * \file        src/mesh/voronoi/voronoi_ghost_search.c
 * \date        05/2018
 * \brief       Algorithms to search for (ghost) cells from other domains.
 * \details     contains functions:
 *                static void particle2in(data_in * in, int i, int firstnode)
 *                static void out2particle(data_out * out, int i, int mode)
 *                static void kernel_local(void)
 *                static void kernel_imported(void)
 *                int voronoi_ghost_search(tessellation * TT)
 *                static void voronoi_pick_up_additional_DP_points(void)
 *                int voronoi_ghost_search_evaluate(tessellation * T,
 *                  int target, int mode, int q, int thread_id)
 *                int ngb_treefind_ghost_search(tessellation * T, MyDouble
 *                  searchcenter[3], MyDouble refpos[3], MyFloat hsml, MyFloat
 *                  maxdist, int target, int origin, int *startnode, int
 *                  bitflags, int mode, int *nexport, int *nsend_local)
 *                int ngb_treefind_ghost_search(tessellation * T, MyDouble
 *                  searchcenter[3], MyDouble refpos[3], MyFloat hsml, MyFloat
 *                  maxdist, int target, int origin, int mode, int thread_id,
 *                  int numnodes, int *firstnode)
 *                int count_undecided_tetras(tessellation * T)
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

#include "../../main/allvars.h"
#include "../../main/proto.h"

#include "voronoi.h"

#if !defined(ONEDIMS)

static void voronoi_pick_up_additional_DP_points(void);

static tessellation *T;

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
  MyDouble RefPos[3];
  MyFloat MaxDist;
  int Origin;

  int Firstnode;

#ifdef EXTENDED_GHOST_SEARCH
  unsigned char BitFlagList[NODELISTLENGTH];
#endif /* #ifdef EXTENDED_GHOST_SEARCH */
} data_in;

static data_in *DataGet, *DataIn;

/*! \brief Routine that fills the relevant particle/cell data into the input
 *         structure defined above. Needed by generic_comm_helpers2.
 *
 *  \param[out] in Data structure to fill.
 *  \param[in] i Index of particle in P and SphP arrays.
 *  \param[in] firstnode First note of communication.
 *
 *  \return void
 */
static void particle2in(data_in *in, int i, int firstnode)
{
  point *DP         = T->DP;
  tetra *DT         = T->DT;
  tetra_center *DTC = T->DTC;

  int k, q;

  for(k = 0, q = -1; k < (NUMDIMS + 1); k++)
    {
#ifndef DOUBLE_STENCIL
      if(DP[DT[i].p[k]].task == ThisTask)
        if(DP[DT[i].p[k]].index >= 0 && DP[DT[i].p[k]].index < NumGas)
          {
            if(TimeBinSynchronized[P[DP[DT[i].p[k]].index].TimeBinHydro])
              {
                q = DT[i].p[k];
                break;
              }
          }
#else  /* #ifndef DOUBLE_STENCIL */
      if(DP[DT[i].p[k]].flag_primary_triangle && DT[i].p[k] >= 0)
        {
          q = DT[i].p[k];
          break;
        }
#endif /* #ifndef DOUBLE_STENCIL #else */
    }

  if(q == -1)
    terminate("q=-1");

  in->Pos[0] = DTC[i].cx;
  in->Pos[1] = DTC[i].cy;
  in->Pos[2] = DTC[i].cz;

  in->RefPos[0] = DP[q].x;
  in->RefPos[1] = DP[q].y;
  in->RefPos[2] = DP[q].z;

  in->Origin = ThisTask;

  in->MaxDist = SphP[DP[q].index].Hsml;

  in->Firstnode = firstnode;
}

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{
  int Count; /* counts how many have been found */
} data_out;

static data_out *DataResult, *DataOut;

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *  particle and cell data arrays (P, SphP,...)
 *  \param[in] i Index of particle in P and SphP arrays
 *  \param[in] mode Mode of function: local particles or information that was
 *             communicated from other tasks and has to be added locally?
 *
 *  \return void
 */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES || mode == MODE_IMPORTED_PARTICLES)
    if(out->Count)
      T->DTF[i] -= (T->DTF[i] & 2);
}

#include "../../utils/generic_comm_helpers2.h"

#ifdef EXTENDED_GHOST_SEARCH
/*! Data structure for extended ghost search.
 */
static struct data_nodelist_special
{
  unsigned char BitFlagList[NODELISTLENGTH];
} * DataNodeListSpecial;
#endif /* #ifdef EXTENDED_GHOST_SEARCH */

static point *DP_Buffer;
static int MaxN_DP_Buffer, N_DP_Buffer;
static int NadditionalPoints;
static int *send_count_new;

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int i, j, q;

  /* do local particles and prepare export list */
  {
    int thread_id = get_thread_num();

    for(j = 0; j < NTask; j++)
      Thread[thread_id].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[thread_id].ExportSpace < MinSpace)
          break;

        i = NextParticle++;

        if(i >= T->Ndt)
          break;

        if((T->DTF[i] & 2) == 0) /* DT that is not flagged as tested ok */
          {
            T->DTF[i] |= 2; /* if we find a particle, need to clear this flag again! */

            point *DP = T->DP;
            tetra *DT = T->DT;

            if(DT[i].t[0] < 0) /* deleted ? */
              continue;

            if(DT[i].p[0] == DPinfinity || DT[i].p[1] == DPinfinity || DT[i].p[2] == DPinfinity)
              continue;

#ifndef TWODIMS
            if(DT[i].p[3] == DPinfinity)
              continue;
#endif /* #ifndef TWODIMS */

#ifndef DOUBLE_STENCIL
            for(j = 0, q = -1; j < (NUMDIMS + 1); j++)
              {
                if(DP[DT[i].p[j]].task == ThisTask)
                  if(DP[DT[i].p[j]].index >= 0 && DP[DT[i].p[j]].index < NumGas)
                    {
                      if(TimeBinSynchronized[P[DP[DT[i].p[j]].index].TimeBinHydro])
                        {
                          q = DT[i].p[j];
                          break;
                        }
                    }
              }

            if(j == (NUMDIMS + 1)) /* this triangle does not have a local point. No need to test it */
              continue;

            if(q == -1)
              terminate("q==-1");
#else  /* #ifndef DOUBLE_STENCIL */
            /* here comes the check for a double stencil */
            for(j = 0, q = -1; j < (NUMDIMS + 1); j++)
              {
                if(DP[DT[i].p[j]].flag_primary_triangle && DT[i].p[j] >= 0)
                  {
                    q = DT[i].p[j];
                    break;
                  }
              }

            if(j ==
               (NUMDIMS +
                1)) /* this triangle does not have a point which is not at least neighbor to a primary point. No need to test it */
              continue;

            if(q == -1)
              terminate("q==-1");
#endif /* #ifndef DOUBLE_STENCIL #else */
            voronoi_ghost_search_evaluate(T, i, MODE_LOCAL_PARTICLES, q, thread_id);
          }
      }
  }
}

/*! \brief Routine that defines what to do with imported particles.
 *
 *  Calls the *_evaluate function in MODE_IMPORTED_PARTICLES.
 *
 *  \return void
 */
static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, count = 0;
  {
    int threadid = get_thread_num();

    while(1)
      {
        i = count++;

        if(i >= Nimport)
          break;

        voronoi_ghost_search_evaluate(T, i, MODE_IMPORTED_PARTICLES, 0, threadid);
      }
  }
}

/*! \brief Main routine to perform ghost search.
 *
 *  \param[in, out] TT Pointer to tessellation.
 *
 *  \return Number of additional points.
 */
int voronoi_ghost_search(tessellation *TT)
{
  T = TT;
  int j, ndone, ndone_flag;

  NadditionalPoints = 0;

  /* allocate buffers to arrange communication */

  send_count_new = (int *)mymalloc_movable(&send_count_new, "send_count_new", NTask * sizeof(int));

  MaxN_DP_Buffer = T->Indi.AllocFacN_DP_Buffer;
  DP_Buffer      = (point *)mymalloc_movable(&DP_Buffer, "DP_Buffer", MaxN_DP_Buffer * sizeof(point));

#ifdef DOUBLE_STENCIL
  {
    point *DP = T->DP;
    tetra *DT = T->DT;
    int i;

    for(i = 0; i < T->Ndp; i++)
      DP[i].flag_primary_triangle = 0;

    for(i = 0; i < T->Ndt; i++)
      {
        for(j = 0; j < (NUMDIMS + 1); j++)
          {
            if(DP[DT[i].p[j]].task == ThisTask)
              if(DP[DT[i].p[j]].index >= 0 && DP[DT[i].p[j]].index < NumGas)
                if(TimeBinSynchronized[P[DP[DT[i].p[j]].index].TimeBinHydro])
                  break;
          }

        if(j != (NUMDIMS + 1)) /* this triangle does have a local point, so mark all its points */
          {
            for(j = 0; j < (NUMDIMS + 1); j++)
              DP[DT[i].p[j]].flag_primary_triangle = 1;
          }
      }
  }
#endif /* #ifdef DOUBLE_STENCIL */

  generic_set_MaxNexport();

  NextParticle = 0;

  do
    {
      for(j = 0; j < NTask; j++)
        send_count_new[j] = 0;

      N_DP_Buffer = 0;

      /* allocate buffers to arrange communication */
      generic_alloc_partlist_nodelist_ngblist_threadbufs();

      kernel_local();

      /* do all necessary bookkeeping and the data exchange */
      generic_exchange(kernel_imported);

      generic_free_partlist_nodelist_ngblist_threadbufs();

      voronoi_pick_up_additional_DP_points();

      if(NextParticle >= T->Ndt)
        ndone_flag = 1;
      else
        ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }
  while(ndone < NTask);

  myfree(DP_Buffer);
  myfree(send_count_new);

#ifdef EXTENDED_GHOST_SEARCH
  myfree(DataNodeListSpecial);
#endif /* #ifdef EXTENDED_GHOST_SEARCH */

  return NadditionalPoints;
}

/*! \brief Gets additional Delaunay points.
 *
 *  \return void
 */
static void voronoi_pick_up_additional_DP_points(void)
{
  int nimport;

  /* The data blocks stored in DP_Buffer is not ordered according to processor rank, but rather in a permutated way.
   * We need to take this into account in calculating the offsets to in the send buffer.
   */

  for(int ngrp = 0, ncnt = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        Send_count[ncnt++] = send_count_new[recvTask];
    }

  Recv_offset[0] = 0;
  for(int j = 1; j < NTask; j++)
    Recv_offset[j] = Recv_offset[j - 1] + Send_count[j - 1];

  for(int ngrp = 0, ncnt = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        Send_offset[recvTask] = Recv_offset[ncnt++];
    }

  memcpy(Send_count, send_count_new, NTask * sizeof(int));

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  Recv_offset[0] = 0;
  nimport        = Recv_count[0];

  for(int j = 1; j < NTask; j++)
    {
      nimport += Recv_count[j];
      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
    }

  while(nimport + T->Ndp > T->MaxNdp)
    {
      T->Indi.AllocFacNdp *= ALLOC_INCREASE_FACTOR;
      T->MaxNdp = T->Indi.AllocFacNdp;
#ifdef VERBOSE
      printf("Task=%d: increase memory allocation, MaxNdp=%d Indi.AllocFacNdp=%g\n", ThisTask, T->MaxNdp, T->Indi.AllocFacNdp);
#endif /* #ifdef VERBOSE */
      T->DP -= 5;
      T->DP = myrealloc_movable(T->DP, (T->MaxNdp + 5) * sizeof(point));
      T->DP += 5;

      if(nimport + T->Ndp > T->MaxNdp && NumGas == 0)
        terminate("nimport + Ndp > MaxNdp");
    }

  /* get the delaunay points */
  for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&DP_Buffer[Send_offset[recvTask]], Send_count[recvTask] * sizeof(point), MPI_BYTE, recvTask, TAG_DENS_B,
                           &T->DP[T->Ndp + Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(point), MPI_BYTE, recvTask,
                           TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  T->Ndp += nimport;
  NadditionalPoints += nimport;

  if(N_DP_Buffer > Largest_N_DP_Buffer)
    Largest_N_DP_Buffer = N_DP_Buffer;
}

/*! \brief Evaluate function for voronoi_ghost_search.
 *
 *  Called in both mode local particles and then in mode imported particles.
 *
 *  \param[] T Pointer to tessellation.
 *  \param[in] target index in DTC and DTF arrays.
 *  \param[in] mode Mode of call (local/imported).
 *  \param[in] q index in DP array.
 *  \param[in] thread_id Thread_id, needed for ngb_treefind_ghost_search.
 *
 *  \return 0
 */
int voronoi_ghost_search_evaluate(tessellation *T, int target, int mode, int q, int thread_id)
{
  int origin, numnodes, *firstnode;
  int numngb;
  double h, dx, dy, dz, maxdist;
  MyDouble pos[3], refpos[3];
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      pos[0]    = T->DTC[target].cx;
      pos[1]    = T->DTC[target].cy;
      pos[2]    = T->DTC[target].cz;
      refpos[0] = T->DP[q].x;
      refpos[1] = T->DP[q].y;
      refpos[2] = T->DP[q].z;
#ifndef DOUBLE_STENCIL
      maxdist = SphP[T->DP[q].index].Hsml;
#else  /* #ifndef DOUBLE_STENCIL */
      maxdist = T->DP[q].Hsml;
#endif /* #ifndef DOUBLE_STENCIL #else */
      origin = ThisTask;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      /* note: we do not use a pointer here to VoroDataGet[target].Pos, because VoroDataGet may be moved in a realloc operation */
      pos[0]    = DataGet[target].Pos[0];
      pos[1]    = DataGet[target].Pos[1];
      pos[2]    = DataGet[target].Pos[2];
      refpos[0] = DataGet[target].RefPos[0];
      refpos[1] = DataGet[target].RefPos[1];
      refpos[2] = DataGet[target].RefPos[2];
      maxdist   = DataGet[target].MaxDist;
      origin    = DataGet[target].Origin;

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  dx = refpos[0] - pos[0];
  dy = refpos[1] - pos[1];
  dz = refpos[2] - pos[2];

  h = 1.0001 * sqrt(dx * dx + dy * dy + dz * dz);

  if(mode == MODE_LOCAL_PARTICLES)
    if(maxdist < 2 * h)
      T->DTF[target] -=
          (T->DTF[target] &
           2); /* since we restrict the search radius, we are not guaranteed to search the full circumcircle of the triangle */

  numngb = ngb_treefind_ghost_search(T, pos, refpos, h, maxdist, target, origin, mode, thread_id, numnodes, firstnode);

  out.Count = numngb;

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#ifdef EXTENDED_GHOST_SEARCH /* this allowes for mirrored images in a full 3x3 grid in terms of the principal domain */
/*! \brief Tree-search algorithm for ghost cells in EXTENDED_GHOST_SEARCH mode.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] searchcenter[3] Postion of the search center.
 *  \param[in] refpos[3] Reference position.
 *  \param[in] hsml Search radius.
 *  \param[in] maxdist Maximum distance.
 *  \param[in] target Index in DTF array.
 *  \param[in] origin Original task.
 *  \param[in] startnode Startnode.
 *  \param[in] bitflags Bitflags for ghost search.
 *  \param[in] mode Mode.
 *  \param[in, out] nexport Number of exported particles.
 *  \param[out] nsend_local Array with number of particles to be sent.
 *
 *  \return Number of points found.
 */
int ngb_treefind_ghost_search(tessellation *T, MyDouble searchcenter[3], MyDouble refpos[3], MyFloat hsml, MyFloat maxdist, int target,
                              int origin, int *startnode, int bitflags, int mode, int *nexport, int *nsend_local)
{
  int i, numngb, no, p, task, nexport_save, ndp_save, nadditionalpoints_save;
  int image_flag;
  struct NgbNODE *current;
  MyDouble dx, dy, dz, hsml2, maxdist2;
  int listp;
  double dx_ref, dy_ref, dz_ref, mindistance, thisdistance;
  double min_x = 0, min_y = 0, min_z = 0;
  int min_p = 0, min_imageflag = 0;
  MyFloat search_min[3], search_max[3], newcenter[3], newrefpos[3];
  MyFloat refsearch_min[3], refsearch_max[3];

  nadditionalpoints_save = NadditionalPoints;
  ndp_save               = T->Ndp;
  nexport_save           = *nexport;

  numngb      = 0;
  mindistance = 1.0e70;

  int repx, repy, repz = 0;
  int repx_A, repy_A, repz_A;
  int repx_B, repy_B, repz_B;
  int xbits;
  int ybits;
  int zbits;
  int count;

  if(mode == 0)
    {
      repx_A = -1;
      repx_B = 1;
      repy_A = -1;
      repy_B = 1;
      repz_A = -1;
      repz_B = 1;
      xbits = ybits = zbits = 0;
    }
  else
    {
      zbits = (bitflags / 9);
      ybits = (bitflags - zbits * 9) / 3;
      xbits = bitflags - zbits * 9 - ybits * 3;

      if(xbits == 1)
        repx_A = repx_B = -1;
      else if(xbits == 2)
        repx_A = repx_B = 1;
      else
        repx_A = repx_B = 0;

      if(ybits == 1)
        repy_A = repy_B = -1;
      else if(ybits == 2)
        repy_A = repy_B = 1;
      else
        repy_A = repy_B = 0;

      if(zbits == 1)
        repz_A = repz_B = -1;
      else if(zbits == 2)
        repz_A = repz_B = 1;
      else
        repz_A = repz_B = 0;
    }

  hsml2    = hsml * hsml;
  maxdist2 = maxdist * maxdist;

  for(repx = repx_A; repx <= repx_B; repx++)
    for(repy = repy_A; repy <= repy_B; repy++)
#if !defined(TWODIMS)
      for(repz = repz_A; repz <= repz_B; repz++)
#endif /* #if !defined(TWODIMS) */
        {
          image_flag = 0; /* for each coordinate there are three possibilities.
                             We encodee them to basis three, i.e. x*3^0 + y*3^1 + z*3^2
                           */
          if(repx == 0)
            {
              newcenter[0] = searchcenter[0];
              newrefpos[0] = refpos[0];
            }
          else if(repx == -1)
            {
#ifndef REFLECTIVE_X
              newcenter[0] = searchcenter[0] - boxSize_X;
              newrefpos[0] = refpos[0] - boxSize_X;
#else  /* #ifndef REFLECTIVE_X */
            newcenter[0] = -searchcenter[0];
            newrefpos[0] = -refpos[0];
#endif /* #ifndef REFLECTIVE_X #else */
              image_flag += 1;
            }
          else /* repx == 1 */
            {
#ifndef REFLECTIVE_X
              newcenter[0] = searchcenter[0] + boxSize_X;
              newrefpos[0] = refpos[0] + boxSize_X;
#else  /* #ifndef REFLECTIVE_X */
            newcenter[0] = -searchcenter[0] + 2 * boxSize_X;
            newrefpos[0] = -refpos[0] + 2 * boxSize_X;
#endif /* #ifndef REFLECTIVE_X #else */
              image_flag += 2;
            }

          if(repy == 0)
            {
              newcenter[1] = searchcenter[1];
              newrefpos[1] = refpos[1];
            }
          else if(repy == -1)
            {
#ifndef REFLECTIVE_Y
              newcenter[1] = searchcenter[1] - boxSize_Y;
              newrefpos[1] = refpos[1] - boxSize_Y;
#else  /* #ifndef REFLECTIVE_Y */
            newcenter[1] = -searchcenter[1];
            newrefpos[1] = -refpos[1];
#endif /* #ifndef REFLECTIVE_Y #else */
              image_flag += 1 * 3;
            }
          else /*  repy == 1 */
            {
#ifndef REFLECTIVE_Y
              newcenter[1] = searchcenter[1] + boxSize_Y;
              newrefpos[1] = refpos[1] + boxSize_Y;
#else  /* #ifndef REFLECTIVE_Y */
            newcenter[1] = -searchcenter[1] + 2 * boxSize_Y;
            newrefpos[1] = -refpos[1] + 2 * boxSize_Y;
#endif /* #ifndef REFLECTIVE_Y #else */
              image_flag += 2 * 3;
            }

          if(repz == 0)
            {
              newcenter[2] = searchcenter[2];
              newrefpos[2] = refpos[2];
            }
#if !defined(TWODIMS)
          else if(repz == -1)
            {
#ifndef REFLECTIVE_Z
              newcenter[2] = searchcenter[2] - boxSize_Z;
              newrefpos[2] = refpos[2] - boxSize_Z;
#else  /* #ifndef REFLECTIVE_Z */
              newcenter[2] = -searchcenter[2];
              newrefpos[2] = -refpos[2];
#endif /* #ifndef REFLECTIVE_Z #else */
              image_flag += 1 * 9;
            }
          else /* repz == 1 */
            {
#ifndef REFLECTIVE_Z
              newcenter[2] = searchcenter[1] + boxSize_Z;
              newrefpos[2] = refpos[1] + boxSize_Z;
#else  /* #ifndef REFLECTIVE_Z */
              newcenter[2] = -searchcenter[2] + 2 * boxSize_Z;
              newrefpos[2] = -refpos[2] + 2 * boxSize_Z;
#endif /* #ifndef REFLECTIVE_Z #else */
              image_flag += 2 * 9;
            }
#endif /* #if !defined(TWODIMS) */

          for(i = 0; i < 3; i++)
            {
              search_min[i]    = newcenter[i] - hsml;
              search_max[i]    = newcenter[i] + hsml;
              refsearch_min[i] = newrefpos[i] - maxdist;
              refsearch_max[i] = newrefpos[i] + maxdist;
            }

          if(mode == 1)
            if(bitflags != image_flag)
              {
                printf("bitflags=%d image_flag=%d xbits=%d ybits=%d zbits=%d  \n", bitflags, image_flag, xbits, ybits, zbits);
                terminate("problem");
              }

          no    = *startnode;
          count = 0;

          while(no >= 0)
            {
              count++;
              if(no < Ngb_MaxPart) /* single particle */
                {
                  p  = no;
                  no = Ngb_Nextnode[no];

                  if(P[p].Type > 0)
                    continue;

                  if(P[p].Mass == 0 && P[p].ID == 0)
                    continue; /* skip cells that have been swallowed or dissolved */

                  dx = P[p].Pos[0] - newcenter[0];
                  dy = P[p].Pos[1] - newcenter[1];
                  dz = P[p].Pos[2] - newcenter[2];

                  if(dx * dx + dy * dy + dz * dz > hsml2)
                    continue;

                  dx_ref = P[p].Pos[0] - newrefpos[0];
                  dy_ref = P[p].Pos[1] - newrefpos[1];
                  dz_ref = P[p].Pos[2] - newrefpos[2];

                  if((thisdistance = dx_ref * dx_ref + dy_ref * dy_ref + dz_ref * dz_ref) > maxdist2)
                    continue;

                  /* now we need to check whether this particle has already been sent to
                     the requesting cpu for this particular image shift */

                  if(thisdistance >= mindistance)
                    continue;

                  if(Ngb_Marker[p] != Ngb_MarkerValue)
                    {
                      Ngb_Marker[p]           = Ngb_MarkerValue;
                      List_P[p].firstexport   = -1;
                      List_P[p].currentexport = -1;
                    }

                  if(List_P[p].firstexport >= 0)
                    {
                      if(ListExports[List_P[p].currentexport].origin != origin)
                        {
                          listp = List_P[p].firstexport;
                          while(listp >= 0)
                            {
                              if(ListExports[listp].origin == origin)
                                {
                                  List_P[p].currentexport = listp;
                                  break;
                                }

                              listp = ListExports[listp].nextexport;
                            }

                          if(listp >= 0)
                            if((ListExports[listp].image_bits & (1 << image_flag))) /* already in list */
                              continue;
                        }
                      else
                        {
                          if((ListExports[List_P[p].currentexport].image_bits & (1 << image_flag))) /* already in list */
                            continue;
                        }
                    }

                  /* here we have found a new closest particle that has not been inserted yet */

                  numngb        = 1;
                  mindistance   = thisdistance;
                  min_p         = p;
                  min_imageflag = image_flag;

                  /* determine the point coordinates in min_x, min_y, min_z */
                  if(repx == 0)
                    min_x = P[p].Pos[0];
                  else if(repx == -1)
                    {
#ifndef REFLECTIVE_X
                      min_x = P[p].Pos[0] + boxSize_X;
#else  /* #ifndef REFLECTIVE_X */
                    min_x = -P[p].Pos[0];
#endif /* #ifndef REFLECTIVE_X #else */
                    }
                  else if(repx == 1)
                    {
#ifndef REFLECTIVE_X
                      min_x = P[p].Pos[0] - boxSize_X;
#else  /* #ifndef REFLECTIVE_X */
                    min_x = -P[p].Pos[0] + 2 * boxSize_X;
#endif /* #ifndef REFLECTIVE_X #else */
                    }

                  if(repy == 0)
                    min_y = P[p].Pos[1];
                  else if(repy == -1)
                    {
#ifndef REFLECTIVE_Y
                      min_y = P[p].Pos[1] + boxSize_Y;
#else  /* #ifndef REFLECTIVE_Y */
                    min_y = -P[p].Pos[1];
#endif /* #ifndef REFLECTIVE_Y #else */
                    }
                  else if(repy == 1)
                    {
#ifndef REFLECTIVE_Y
                      min_y = P[p].Pos[1] - boxSize_Y;
#else  /* #ifndef REFLECTIVE_Y */
                    min_y = -P[p].Pos[1] + 2 * boxSize_Y;
#endif /* #ifndef REFLECTIVE_Y #else */
                    }

                  if(repz == 0)
                    min_z = P[p].Pos[2];
#if !defined(TWODIMS)
                  else if(repz == -1)
                    {
#ifndef REFLECTIVE_Z
                      min_z = P[p].Pos[2] + boxSize_Z;
#else  /* #ifndef REFLECTIVE_Z */
                      min_z = -P[p].Pos[2];
#endif /* #ifndef REFLECTIVE_Z #else */
                    }
                  else if(repz == 1)
                    {
#ifndef REFLECTIVE_Z
                      min_z = P[p].Pos[2] - boxSize_Z;
#else  /* #ifndef REFLECTIVE_Z */
                      min_z = -P[p].Pos[2] + 2 * boxSize_Z;
#endif /* #ifndef REFLECTIVE_Z #else */
                    }
#endif /* #if !defined(TWODIMS) */
                }
              else if(no < Ngb_MaxPart + Ngb_MaxNodes) /* internal node */
                {
                  if(mode == 1)
                    {
                      if(no < Ngb_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the
                                                           branch */
                        {
                          break;
                        }
                    }

                  current = &Ngb_Nodes[no];
                  no      = current->u.d.sibling; /* in case the node can be discarded */

                  if(search_min[0] > current->u.d.range_max[0])
                    continue;
                  if(search_max[0] < current->u.d.range_min[0])
                    continue;
                  if(refsearch_min[0] > current->u.d.range_max[0])
                    continue;
                  if(refsearch_max[0] < current->u.d.range_min[0])
                    continue;

                  if(search_min[1] > current->u.d.range_max[1])
                    continue;
                  if(search_max[1] < current->u.d.range_min[1])
                    continue;
                  if(refsearch_min[1] > current->u.d.range_max[1])
                    continue;
                  if(refsearch_max[1] < current->u.d.range_min[1])
                    continue;

                  if(search_min[2] > current->u.d.range_max[2])
                    continue;
                  if(search_max[2] < current->u.d.range_min[2])
                    continue;
                  if(refsearch_min[2] > current->u.d.range_max[2])
                    continue;
                  if(refsearch_max[2] < current->u.d.range_min[2])
                    continue;

                  no = current->u.d.nextnode; /* ok, we need to open the node */
                }
              else /* pseudo particle */
                {
                  if(mode == 1)
                    terminate("mode == 1");

                  if(target >= 0) /* if no target is given, export will not occur */
                    {
                      if(Exportflag[task = DomainTask[no - (Ngb_MaxPart + Ngb_MaxNodes)]] != target)
                        {
                          Exportflag[task]      = target;
                          Exportnodecount[task] = NODELISTLENGTH;
                        }

                      if(Exportnodecount[task] == NODELISTLENGTH)
                        {
                          if(*nexport >= All.BunchSize)
                            {
                              T->Ndp            = ndp_save;
                              NadditionalPoints = nadditionalpoints_save;
                              *nexport          = nexport_save;
                              if(nexport_save == 0)
                                terminate(
                                    "nexport_save == 0"); /* in this case, the buffer is too small to process even a single particle */
                              for(task = 0; task < NTask; task++)
                                nsend_local[task] = 0;
                              for(no = 0; no < nexport_save; no++)
                                nsend_local[DataIndexTable[no].Task]++;
                              return -1;
                            }
                          Exportnodecount[task]             = 0;
                          Exportindex[task]                 = *nexport;
                          DataIndexTable[*nexport].Task     = task;
                          DataIndexTable[*nexport].Index    = target;
                          DataIndexTable[*nexport].IndexGet = *nexport;
                          *nexport                          = *nexport + 1;
                          nsend_local[task]++;
                        }

                      DataNodeListSpecial[Exportindex[task]].BitFlagList[Exportnodecount[task]] = image_flag;
                      DataNodeListSpecial[Exportindex[task]].NodeList[Exportnodecount[task]++] =
                          Ngb_DomainNodeIndex[no - (Ngb_MaxPart + Ngb_MaxNodes)];

                      if(Exportnodecount[task] < NODELISTLENGTH)
                        DataNodeListSpecial[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
                    }

                  no = Ngb_Nextnode[no - Ngb_MaxNodes];
                  continue;
                }
            }
        }

  *startnode = -1;

  if(numngb)
    {
      p = min_p;

      image_flag = min_imageflag;

      if(Ngb_Marker[p] != Ngb_MarkerValue)
        {
          Ngb_Marker[p]           = Ngb_MarkerValue;
          List_P[p].firstexport   = -1;
          List_P[p].currentexport = -1;
        }

      if(List_P[p].firstexport >= 0)
        {
          if(ListExports[List_P[p].currentexport].origin != origin)
            {
              listp = List_P[p].firstexport;
              while(listp >= 0)
                {
                  if(ListExports[listp].origin == origin)
                    {
                      List_P[p].currentexport = listp;
                      break;
                    }

                  if(ListExports[listp].nextexport < 0)
                    {
                      if(Ninlist >= MaxNinlist)
                        {
                          T->Indi.AllocFacNinlist *= ALLOC_INCREASE_FACTOR;
                          MaxNinlist = T->Indi.AllocFacNinlist;
#ifdef VERBOSE
                          printf("Task=%d: increase memory allocation, MaxNinlist=%d Indi.AllocFacNinlist=%g\n", ThisTask, MaxNinlist,
                                 T->Indi.AllocFacNinlist);
#endif /* #ifdef VERBOSE */
                          ListExports = myrealloc_movable(ListExports, MaxNinlist * sizeof(struct list_export_data));

                          if(Ninlist >= MaxNinlist)
                            terminate("Ninlist >= MaxNinlist");
                        }

                      List_P[p].currentexport                         = Ninlist++;
                      ListExports[List_P[p].currentexport].image_bits = 0;
                      ListExports[List_P[p].currentexport].nextexport = -1;
                      ListExports[List_P[p].currentexport].origin     = origin;
                      ListExports[List_P[p].currentexport].index      = p;
                      ListExports[listp].nextexport                   = List_P[p].currentexport;
                      break;
                    }
                  listp = ListExports[listp].nextexport;
                }
            }
        }
      else
        {
          /* here we have a local particle that hasn't been made part of the mesh */

          if(Ninlist >= MaxNinlist)
            {
              T->Indi.AllocFacNinlist *= ALLOC_INCREASE_FACTOR;
              MaxNinlist = T->Indi.AllocFacNinlist;
#ifdef VERBOSE
              printf("Task=%d: increase memory allocation, MaxNinlist=%d Indi.AllocFacNinlist=%g\n", ThisTask, MaxNinlist,
                     T->Indi.AllocFacNinlist);
#endif /* #ifdef VERBOSE */
              ListExports = myrealloc_movable(ListExports, MaxNinlist * sizeof(struct list_export_data));

              if(Ninlist >= MaxNinlist)
                terminate("Ninlist >= MaxNinlist");
            }

          List_InMesh[NumGasInMesh++] = p;

          List_P[p].currentexport = List_P[p].firstexport = Ninlist++;
          ListExports[List_P[p].currentexport].image_bits = 0;
          ListExports[List_P[p].currentexport].nextexport = -1;
          ListExports[List_P[p].currentexport].origin     = origin;
          ListExports[List_P[p].currentexport].index      = p;
        }

      if((ListExports[List_P[p].currentexport].image_bits & (1 << image_flag)))
        terminate("this should not happen");

      ListExports[List_P[p].currentexport].image_bits |= (1 << image_flag);

      /* add the particle to the ones that need to be exported */

      if(origin == ThisTask)
        {
          if(mode == 1)
            terminate("mode==1: how can this be?");

          if(T->Ndp >= T->MaxNdp)
            {
              T->Indi.AllocFacNdp *= ALLOC_INCREASE_FACTOR;
              T->MaxNdp = T->Indi.AllocFacNdp;
#ifdef VERBOSE
              printf("Task=%d: increase memory allocation, MaxNdp=%d Indi.AllocFacNdp=%g\n", ThisTask, T->MaxNdp, T->Indi.AllocFacNdp);
#endif /* #ifdef VERBOSE */
              T->DP -= 5;
              T->DP = myrealloc_movable(T->DP, (T->MaxNdp + 5) * sizeof(point));
              T->DP += 5;

              if(T->Ndp >= T->MaxNdp)
                terminate("Ndp >= MaxNdp");
            }

          SphP[p].ActiveArea = 0;

          point *dp = &T->DP[T->Ndp];
          dp->x     = min_x;
          dp->y     = min_y;
          dp->z     = min_z;
          dp->task  = ThisTask;
          dp->ID    = P[p].ID;
          if(image_flag)
            dp->index = p + NumGas; /* this is a replicated/mirrored local point */
          else
            dp->index = p; /* this is actually a local point that wasn't made part of the mesh yet */
          dp->originalindex = p;
          dp->timebin       = P[p].TimeBinHydro;
          dp->image_flags   = (1 << image_flag);

#ifdef DOUBLE_STENCIL
          dp->Hsml             = SphP[p].Hsml;
          dp->first_connection = -1;
          dp->last_connection  = -1;
#endif /* #ifdef DOUBLE_STENCIL */
          T->Ndp++;
          NadditionalPoints++;
        }
      else
        {
          if(mode == 0)
            terminate("mode == 0: how can this be?");

          if(N_DP_Buffer >= MaxN_DP_Buffer)
            {
              T->Indi.AllocFacN_DP_Buffer *= ALLOC_INCREASE_FACTOR;
              MaxN_DP_Buffer = T->Indi.AllocFacN_DP_Buffer;
#ifdef VERBOSE
              printf("Task=%d: increase memory allocation, MaxN_DP_Buffer=%d Indi.AllocFacN_DP_Buffer=%g\n", ThisTask, MaxN_DP_Buffer,
                     T->Indi.AllocFacN_DP_Buffer);
#endif /* #ifdef VERBOSE */
              DP_Buffer = (point *)myrealloc_movable(DP_Buffer, MaxN_DP_Buffer * sizeof(point));

              if(N_DP_Buffer >= MaxN_DP_Buffer)
                terminate("(N_DP_Buffer >= MaxN_DP_Buffer");
            }

          SphP[p].ActiveArea = 0;

          DP_Buffer[N_DP_Buffer].x             = min_x;
          DP_Buffer[N_DP_Buffer].y             = min_y;
          DP_Buffer[N_DP_Buffer].z             = min_z;
          DP_Buffer[N_DP_Buffer].ID            = P[p].ID;
          DP_Buffer[N_DP_Buffer].task          = ThisTask;
          DP_Buffer[N_DP_Buffer].index         = p;
          DP_Buffer[N_DP_Buffer].originalindex = p;
          DP_Buffer[N_DP_Buffer].timebin       = P[p].TimeBinHydro;
          DP_Buffer[N_DP_Buffer].image_flags   = (1 << image_flag);
#ifdef DOUBLE_STENCIL
          DP_Buffer[N_DP_Buffer].Hsml             = SphP[p].Hsml;
          DP_Buffer[N_DP_Buffer].first_connection = -1;
          DP_Buffer[N_DP_Buffer].last_connection  = -1;
#endif /* #ifdef DOUBLE_STENCIL */
          send_count_new[origin]++;
          N_DP_Buffer++;
        }
    }

  return numngb;
}

#else /* #ifdef EXTENDED_GHOST_SEARCH */

/*! \brief Tree-search algorithm for ghost cells without EXTENDED_GHOST_SEARCH.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] searchcenter[3] Postion of the search center.
 *  \param[in] refpos[3] Reference position.
 *  \param[in] hsml Search radius.
 *  \param[in] maxdist Maximum distance.
 *  \param[in] target Index in DTF array.
 *  \param[in] origin Original task.
 *  \param[in] mode Mode (local/imported).
 *  \param[in] thread_id ID of this thread.
 *  \param[in] numnodes Number of nodes.
 *  \param[in] firstnode Index of first node.
 *
 *  \return Number of points found.
 */
int ngb_treefind_ghost_search(tessellation *T, MyDouble searchcenter[3], MyDouble refpos[3], MyFloat hsml, MyFloat maxdist, int target,
                              int origin, int mode, int thread_id, int numnodes, int *firstnode)
{
  int i, k, numngb, no, p;
  int image_flag = 0;
  struct NgbNODE *current;
  MyDouble x, y, z, dx, dy, dz;
  int listp;
  double dx_ref, dy_ref, dz_ref, mindistance, thisdistance, maxdistSquared, hsmlSquared;
  double min_x = 0, min_y = 0, min_z = 0;
  int min_p = 0, min_imageflag = 0;
  double offx, offy, offz;
  MyFloat search_min[3], search_max[3], search_max_Lsub[3], search_min_Ladd[3];
  MyFloat refsearch_min[3], refsearch_max[3], refsearch_max_Lsub[3], refsearch_min_Ladd[3];

  for(i = 0; i < 3; i++)
    {
      search_min[i] = searchcenter[i] - hsml;
      search_max[i] = searchcenter[i] + hsml;
      refsearch_min[i] = refpos[i] - maxdist;
      refsearch_max[i] = refpos[i] + maxdist;
    }

#if !defined(REFLECTIVE_X)
  search_max_Lsub[0] = search_max[0] - boxSize_X;
  search_min_Ladd[0] = search_min[0] + boxSize_X;
  refsearch_max_Lsub[0] = refsearch_max[0] - boxSize_X;
  refsearch_min_Ladd[0] = refsearch_min[0] + boxSize_X;
#else  /* #if !defined(REFLECTIVE_X) */
  search_max_Lsub[0]    = 2 * boxSize_X - search_max[0];
  search_min_Ladd[0]    = -search_min[0];
  refsearch_max_Lsub[0] = 2 * boxSize_X - refsearch_max[0];
  refsearch_min_Ladd[0] = -refsearch_min[0];
#endif /* #if !defined(REFLECTIVE_X) #else */

#if !defined(REFLECTIVE_Y)
  search_max_Lsub[1] = search_max[1] - boxSize_Y;
  search_min_Ladd[1] = search_min[1] + boxSize_Y;
  refsearch_max_Lsub[1] = refsearch_max[1] - boxSize_Y;
  refsearch_min_Ladd[1] = refsearch_min[1] + boxSize_Y;
#else  /* #if !defined(REFLECTIVE_Y) */
  search_max_Lsub[1]    = 2 * boxSize_Y - search_max[1];
  search_min_Ladd[1]    = -search_min[1];
  refsearch_max_Lsub[1] = 2 * boxSize_Y - refsearch_max[1];
  refsearch_min_Ladd[1] = -refsearch_min[1];
#endif /* #if !defined(REFLECTIVE_Y) #else */

#if !defined(REFLECTIVE_Z)
  search_max_Lsub[2] = search_max[2] - boxSize_Z;
  search_min_Ladd[2] = search_min[2] + boxSize_Z;
  refsearch_max_Lsub[2] = refsearch_max[2] - boxSize_Z;
  refsearch_min_Ladd[2] = refsearch_min[2] + boxSize_Z;
#else  /* #if !defined(REFLECTIVE_Z) */
  search_max_Lsub[2]    = 2 * boxSize_Z - search_max[2];
  search_min_Ladd[2]    = -search_min[2];
  refsearch_max_Lsub[2] = 2 * boxSize_Z - refsearch_max[2];
  refsearch_min_Ladd[2] = -refsearch_min[2];
#endif /* #if !defined(REFLECTIVE_Z) #else */

  numngb = 0;
  mindistance = 1.0e70;
  int count;

  count = 0;

  maxdistSquared = maxdist * maxdist;
  hsmlSquared = hsml * hsml;

  numngb = 0;

  for(k = 0; k < numnodes; k++)
    {
      if(mode == MODE_LOCAL_PARTICLES)
        {
          no = Ngb_MaxPart; /* root node */

#ifdef EXTENDED_GHOST_SEARCH
          bitflags = 0;
#endif /* #ifdef EXTENDED_GHOST_SEARCH */
        }
      else
        {
          no = firstnode[k];

#ifdef EXTENDED_GHOST_SEARCH
          bitflags = first_bitflag[k];
#endif /* #ifdef EXTENDED_GHOST_SEARCH */
          no = Ngb_Nodes[no].u.d.nextnode; /* open it */
        }

      while(no >= 0)
        {
          count++;
          if(no < Ngb_MaxPart) /* single particle */
            {
              p = no;
              no = Ngb_Nextnode[no];

              if(P[p].Type > 0)
                continue;

              if(P[p].Mass == 0 && P[p].ID == 0)
                continue; /* skip cells that have been swallowed or eliminated */

              if(P[p].Ti_Current != All.Ti_Current)
                {
                  drift_particle(p, All.Ti_Current);
                }

              offx = offy = offz = 0;

              image_flag = 0; /* for each coordinates there are three possibilities. We
                                 encode them to basis three, i.e. x*3^0 + y*3^1 + z*3^2 */

#if !defined(REFLECTIVE_X)
              if(P[p].Pos[0] - refpos[0] < -boxHalf_X)
                {
                  offx = boxSize_X;
                  image_flag += 1;
                }
              else if(P[p].Pos[0] - refpos[0] > boxHalf_X)
                {
                  offx = -boxSize_X;
                  image_flag += 2;
                }
#endif /* #if !defined(REFLECTIVE_X) */

#if !defined(REFLECTIVE_Y)
              if(P[p].Pos[1] - refpos[1] < -boxHalf_Y)
                {
                  offy = boxSize_Y;
                  image_flag += 1 * 3;
                }
              else if(P[p].Pos[1] - refpos[1] > boxHalf_Y)
                {
                  offy = -boxSize_Y;
                  image_flag += 2 * 3;
                }
#endif /* #if !defined(REFLECTIVE_Y) */

#if !defined(REFLECTIVE_Z) && !defined(TWODIMS)
              if(P[p].Pos[2] - refpos[2] < -boxHalf_Z)
                {
                  offz = boxSize_Z;
                  image_flag += 1 * 9;
                }
              else if(P[p].Pos[2] - refpos[2] > boxHalf_Z)
                {
                  offz = -boxSize_Z;
                  image_flag += 2 * 9;
                }
#endif /* #if !defined(REFLECTIVE_Z) && !defined(TWODIMS) */

              int image_flag_periodic_bnds = image_flag;

#if defined(REFLECTIVE_X)
              int repx;
              for(repx = -1; repx <= 1; repx++, offx = 0)
#endif /* #if defined(REFLECTIVE_X) */
                {
#if defined(REFLECTIVE_Y)
                  int repy;
                  for(repy = -1; repy <= 1; repy++, offy = 0)
#endif /* #if defined(REFLECTIVE_Y) */
                    {
#if defined(REFLECTIVE_Z) && !defined(TWODIMS)
                      int repz;
                      for(repz = -1; repz <= 1; repz++, offz = 0)
#endif /* #if defined(REFLECTIVE_Z) && !defined(TWODIMS) */
                        {
                          image_flag = image_flag_periodic_bnds;

                          x = P[p].Pos[0];
                          y = P[p].Pos[1];
                          z = P[p].Pos[2];

#if defined(REFLECTIVE_X)
                          if(repx == 1)
                            {
                              offx = 2 * boxSize_X;
                              image_flag += 2;
                            }
                          else if(repx == -1)
                            {
                              image_flag += 1;
                            }
                          if(repx != 0)
                            x = -x;
#endif /* #if defined(REFLECTIVE_X) */

#if defined(REFLECTIVE_Y)
                          if(repy == 1)
                            {
                              offy = 2 * boxSize_Y;
                              image_flag += 2 * 3;
                            }
                          else if(repy == -1)
                            {
                              image_flag += 1 * 3;
                            }
                          if(repy != 0)
                            y = -y;
#endif /* #if  defined(REFLECTIVE_Y) */

#if defined(REFLECTIVE_Z) && !defined(TWODIMS)
                          if(repz == 1)
                            {
                              offz = 2 * boxSize_Z;
                              image_flag += 2 * 9;
                            }
                          else if(repz == -1)
                            {
                              image_flag += 1 * 9;
                            }
                          if(repz != 0)
                            z = -z;
#endif /* #if  defined(REFLECTIVE_Z) && !defined(TWODIMS) */

                          x += offx;
                          y += offy;
                          z += offz;

                          dx_ref = x - refpos[0];
                          dy_ref = y - refpos[1];
                          dz_ref = z - refpos[2];

                          if((thisdistance = dx_ref * dx_ref + dy_ref * dy_ref + dz_ref * dz_ref) > maxdistSquared)
                            continue;

                          dx = x - searchcenter[0];
                          dy = y - searchcenter[1];
                          dz = z - searchcenter[2];

                          if(dx * dx + dy * dy + dz * dz > hsmlSquared)
                            continue;

                          /* now we need to check whether this particle has already been sent to
                             the requesting cpu for this particular image shift */

                          if(thisdistance >= mindistance)
                            continue;

                          if(Ngb_Marker[p] != Ngb_MarkerValue)
                            {
                              Ngb_Marker[p] = Ngb_MarkerValue;
                              List_P[p].firstexport = -1;
                              List_P[p].currentexport = -1;
                            }

                          if(List_P[p].firstexport >= 0)
                            {
                              if(ListExports[List_P[p].currentexport].origin != origin)
                                {
                                  listp = List_P[p].firstexport;
                                  while(listp >= 0)
                                    {
                                      if(ListExports[listp].origin == origin)
                                        {
                                          List_P[p].currentexport = listp;
                                          break;
                                        }

                                      listp = ListExports[listp].nextexport;
                                    }

                                  if(listp >= 0)
                                    if((ListExports[listp].image_bits & (1 << image_flag))) /* already in list */
                                      continue;
                                }
                              else
                                {
                                  if((ListExports[List_P[p].currentexport].image_bits & (1 << image_flag))) /* already in list */
                                    continue;
                                }
                            }

                          /* here we have found a new closest particle that has not been inserted yet */

                          numngb = 1;
                          mindistance = thisdistance;
                          min_p = p;
                          min_imageflag = image_flag;
                          min_x = x;
                          min_y = y;
                          min_z = z;

                          maxdistSquared = thisdistance;
                        }
                    }
                }
            }
          else if(no < Ngb_MaxPart + Ngb_MaxNodes) /* internal node */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no <
                     Ngb_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the branch */
                    break;
                }

              current = &Ngb_Nodes[no];
              no = current->u.d.sibling; /* in case the node can be discarded */

              if(current->Ti_Current != All.Ti_Current)
                {
                  drift_node(current, All.Ti_Current);
                }

#if !defined(REFLECTIVE_X)
              if(search_min[0] > current->u.d.range_max[0] && search_max_Lsub[0] < current->u.d.range_min[0])
                continue;
              if(search_min_Ladd[0] > current->u.d.range_max[0] && search_max[0] < current->u.d.range_min[0])
                continue;
#else  /* #if !defined(REFLECTIVE_X) */
              if(search_min[0] > current->u.d.range_max[0] && search_max_Lsub[0] > current->u.d.range_max[0])
                continue;
              if(search_min_Ladd[0] < current->u.d.range_min[0] && search_max[0] < current->u.d.range_min[0])
                continue;
#endif /* #if !defined(REFLECTIVE_X) #else */

#if !defined(REFLECTIVE_Y)
              if(search_min[1] > current->u.d.range_max[1] && search_max_Lsub[1] < current->u.d.range_min[1])
                continue;
              if(search_min_Ladd[1] > current->u.d.range_max[1] && search_max[1] < current->u.d.range_min[1])
                continue;
#else  /* #if !defined(REFLECTIVE_Y) */
              if(search_min[1] > current->u.d.range_max[1] && search_max_Lsub[1] > current->u.d.range_max[1])
                continue;
              if(search_min_Ladd[1] < current->u.d.range_min[1] && search_max[1] < current->u.d.range_min[1])
                continue;
#endif /* #if !defined(REFLECTIVE_Y) #else */

#if !defined(REFLECTIVE_Z)
              if(search_min[2] > current->u.d.range_max[2] && search_max_Lsub[2] < current->u.d.range_min[2])
                continue;
              if(search_min_Ladd[2] > current->u.d.range_max[2] && search_max[2] < current->u.d.range_min[2])
                continue;
#else  /* #if !defined(REFLECTIVE_Z) */
              if(search_min[2] > current->u.d.range_max[2] && search_max_Lsub[2] > current->u.d.range_max[2])
                continue;
              if(search_min_Ladd[2] < current->u.d.range_min[2] && search_max[2] < current->u.d.range_min[2])
                continue;
#endif /* #if !defined(REFLECTIVE_Z) #else */

                /* now deal with the search region of the reference point */

#if !defined(REFLECTIVE_X)
              if(refsearch_min[0] > current->u.d.range_max[0] && refsearch_max_Lsub[0] < current->u.d.range_min[0])
                continue;
              if(refsearch_min_Ladd[0] > current->u.d.range_max[0] && refsearch_max[0] < current->u.d.range_min[0])
                continue;
#else  /* #if !defined(REFLECTIVE_X) */
              if(refsearch_min[0] > current->u.d.range_max[0] && refsearch_max_Lsub[0] > current->u.d.range_max[0])
                continue;
              if(refsearch_min_Ladd[0] < current->u.d.range_min[0] && refsearch_max[0] < current->u.d.range_min[0])
                continue;
#endif /* #if !defined(REFLECTIVE_X) #else */

#if !defined(REFLECTIVE_Y)
              if(refsearch_min[1] > current->u.d.range_max[1] && refsearch_max_Lsub[1] < current->u.d.range_min[1])
                continue;
              if(refsearch_min_Ladd[1] > current->u.d.range_max[1] && refsearch_max[1] < current->u.d.range_min[1])
                continue;
#else  /* #if !defined(REFLECTIVE_Y) */
              if(refsearch_min[1] > current->u.d.range_max[1] && refsearch_max_Lsub[1] > current->u.d.range_max[1])
                continue;
              if(refsearch_min_Ladd[1] < current->u.d.range_min[1] && refsearch_max[1] < current->u.d.range_min[1])
                continue;
#endif /* #if !defined(REFLECTIVE_Y) #else */

#if !defined(REFLECTIVE_Z)
              if(refsearch_min[2] > current->u.d.range_max[2] && refsearch_max_Lsub[2] < current->u.d.range_min[2])
                continue;
              if(refsearch_min_Ladd[2] > current->u.d.range_max[2] && refsearch_max[2] < current->u.d.range_min[2])
                continue;
#else  /* #if !defined(REFLECTIVE_Z) */
              if(refsearch_min[2] > current->u.d.range_max[2] && refsearch_max_Lsub[2] > current->u.d.range_max[2])
                continue;
              if(refsearch_min_Ladd[2] < current->u.d.range_min[2] && refsearch_max[2] < current->u.d.range_min[2])
                continue;
#endif /* #if !defined(REFLECTIVE_Z) #else */

              no = current->u.d.nextnode; /* ok, we need to open the node */
            }
          else /* pseudo particle */
            {
              if(mode == 1)
                terminate("mode == 1");

              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES should not occur here");

              if(target >= 0) /* if no target is given, export will not occur */
                ngb_treefind_export_node_threads(no, target, thread_id, image_flag);

              no = Ngb_Nextnode[no - Ngb_MaxNodes];
              continue;
            }
        }
    }

  if(numngb)
    {
      p = min_p;

      image_flag = min_imageflag;

      if(Ngb_Marker[p] != Ngb_MarkerValue)
        {
          Ngb_Marker[p] = Ngb_MarkerValue;
          List_P[p].firstexport = -1;
          List_P[p].currentexport = -1;
        }

      if(List_P[p].firstexport >= 0)
        {
          if(ListExports[List_P[p].currentexport].origin != origin)
            {
              listp = List_P[p].firstexport;
              while(listp >= 0)
                {
                  if(ListExports[listp].origin == origin)
                    {
                      List_P[p].currentexport = listp;
                      break;
                    }

                  if(ListExports[listp].nextexport < 0)
                    {
                      if(Ninlist >= MaxNinlist)
                        {
                          T->Indi.AllocFacNinlist *= ALLOC_INCREASE_FACTOR;
                          MaxNinlist = T->Indi.AllocFacNinlist;
#ifdef VERBOSE
                          printf("Task=%d: increase memory allocation, MaxNinlist=%d Indi.AllocFacNinlist=%g\n", ThisTask, MaxNinlist,
                                 T->Indi.AllocFacNinlist);
#endif /* #ifdef VERBOSE */
                          ListExports = myrealloc_movable(ListExports, MaxNinlist * sizeof(struct list_export_data));

                          if(Ninlist >= MaxNinlist)
                            terminate("Ninlist >= MaxNinlist");
                        }

                      List_P[p].currentexport = Ninlist++;
                      ListExports[List_P[p].currentexport].image_bits = 0;
                      ListExports[List_P[p].currentexport].nextexport = -1;
                      ListExports[List_P[p].currentexport].origin = origin;
                      ListExports[List_P[p].currentexport].index = p;
                      ListExports[listp].nextexport = List_P[p].currentexport;
                      break;
                    }
                  listp = ListExports[listp].nextexport;
                }
            }
        }
      else
        {
          /* here we have a local particle that hasn't been made part of the mesh */

          if(Ninlist >= MaxNinlist)
            {
              T->Indi.AllocFacNinlist *= ALLOC_INCREASE_FACTOR;
              MaxNinlist = T->Indi.AllocFacNinlist;
#ifdef VERBOSE
              printf("Task=%d: increase memory allocation, MaxNinlist=%d Indi.AllocFacNinlist=%g\n", ThisTask, MaxNinlist,
                     T->Indi.AllocFacNinlist);
#endif /* #ifdef VERBOSE */
              ListExports = myrealloc_movable(ListExports, MaxNinlist * sizeof(struct list_export_data));

              if(Ninlist >= MaxNinlist)
                terminate("Ninlist >= MaxNinlist");
            }

          List_InMesh[NumGasInMesh++] = p;

          List_P[p].currentexport = List_P[p].firstexport = Ninlist++;
          ListExports[List_P[p].currentexport].image_bits = 0;
          ListExports[List_P[p].currentexport].nextexport = -1;
          ListExports[List_P[p].currentexport].origin = origin;
          ListExports[List_P[p].currentexport].index = p;
        }

      if((ListExports[List_P[p].currentexport].image_bits & (1 << image_flag)))
        terminate("this should not happen");

      ListExports[List_P[p].currentexport].image_bits |= (1 << image_flag);

      /* add the particle to the ones that need to be exported */

      if(P[p].Ti_Current != All.Ti_Current)
        terminate("surprise! we don't expect this here anymore");

      if(origin == ThisTask)
        {
          if(mode == 1)
            terminate("mode==1: how can this be?");

          if(T->Ndp >= T->MaxNdp)
            {
              T->Indi.AllocFacNdp *= ALLOC_INCREASE_FACTOR;
              T->MaxNdp = T->Indi.AllocFacNdp;
#ifdef VERBOSE
              printf("Task=%d: increase memory allocation, MaxNdp=%d Indi.AllocFacNdp=%g\n", ThisTask, T->MaxNdp, T->Indi.AllocFacNdp);
#endif /* #ifdef VERBOSE */
              T->DP -= 5;
              T->DP = myrealloc_movable(T->DP, (T->MaxNdp + 5) * sizeof(point));
              T->DP += 5;

              if(T->Ndp >= T->MaxNdp)
                terminate("Ndp >= MaxNdp");
            }

          SphP[p].ActiveArea = 0;

          point *dp = &T->DP[T->Ndp];
          dp->x = min_x;
          dp->y = min_y;
          dp->z = min_z;
          dp->task = ThisTask;
          dp->ID = P[p].ID;
          if(image_flag)
            dp->index = p + NumGas; /* this is a replicated/mirrored local point */
          else
            dp->index = p; /* this is actually a local point that wasn't made part of the mesh yet */
          dp->originalindex = p;
          dp->timebin = P[p].TimeBinHydro;
          dp->image_flags = (1 << image_flag);
#ifdef DOUBLE_STENCIL
          dp->Hsml = SphP[p].Hsml;
          dp->first_connection = -1;
          dp->last_connection = -1;
#endif /* #ifdef DOUBLE_STENCIL */
          T->Ndp++;
          NadditionalPoints++;
        }
      else
        {
          if(mode == 0)
            terminate("mode == 0: how can this be?");

          if(N_DP_Buffer >= MaxN_DP_Buffer)
            {
              T->Indi.AllocFacN_DP_Buffer *= ALLOC_INCREASE_FACTOR;
              MaxN_DP_Buffer = T->Indi.AllocFacN_DP_Buffer;
#ifdef VERBOSE
              printf("Task=%d: increase memory allocation, MaxN_DP_Buffer=%d Indi.AllocFacN_DP_Buffer=%g\n", ThisTask, MaxN_DP_Buffer,
                     T->Indi.AllocFacN_DP_Buffer);
#endif /* #ifdef VERBOSE */
              DP_Buffer = (point *)myrealloc_movable(DP_Buffer, MaxN_DP_Buffer * sizeof(point));

              if(N_DP_Buffer >= MaxN_DP_Buffer)
                terminate("(N_DP_Buffer >= MaxN_DP_Buffer");
            }

          SphP[p].ActiveArea = 0;

          DP_Buffer[N_DP_Buffer].x = min_x;
          DP_Buffer[N_DP_Buffer].y = min_y;
          DP_Buffer[N_DP_Buffer].z = min_z;
          DP_Buffer[N_DP_Buffer].ID = P[p].ID;
          DP_Buffer[N_DP_Buffer].task = ThisTask;
          DP_Buffer[N_DP_Buffer].index = p;
          DP_Buffer[N_DP_Buffer].originalindex = p;
          DP_Buffer[N_DP_Buffer].timebin = P[p].TimeBinHydro;
          DP_Buffer[N_DP_Buffer].image_flags = (1 << image_flag);
#ifdef DOUBLE_STENCIL
          DP_Buffer[N_DP_Buffer].Hsml = SphP[p].Hsml;
          DP_Buffer[N_DP_Buffer].first_connection = -1;
          DP_Buffer[N_DP_Buffer].last_connection = -1;
#endif /* #ifdef DOUBLE_STENCIL */
          send_count_new[origin]++;
          N_DP_Buffer++;
        }
    }

  return numngb;
}

#endif /* #ifdef EXTENDED_GHOST_SEARCH #else */

/*! \brief Counts up undecided tetrahedra.
 *
 *  \param[in] T Pointer to tessellation.
 *
 *  \return (Local) number of undecided tetrahedra.
 */
int count_undecided_tetras(tessellation *T)
{
  int i, count;

  for(i = 0, count = 0; i < T->Ndt; i++)
    if((T->DTF[i] & 2) == 0)
      count++;

  return count;
}

#endif /* #if !defined(ONEDIMS) */
