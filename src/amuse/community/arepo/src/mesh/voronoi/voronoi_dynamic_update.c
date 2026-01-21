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
 * \file        src/mesh/voronoi/voronoi_dynamic_update.c
 * \date        05/2018
 * \brief       Algorithms for Voronoi dynamic update.
 * \details     contains functions:
 *                int voronoi_get_connected_particles(tessellation * T)
 *                void voronoi_init_connectivity(tessellation * T)
 *                void voronoi_update_connectivity(tessellation * T)
 *                void voronoi_remove_connection(int i)
 *                int compare_foreign_connection(const void *a, const void *b)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 22.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../main/allvars.h"
#include "../../main/proto.h"

#include "voronoi.h"

int Nvc;    /* number of connections */
int MaxNvc; /* maximum number of connections */
int Largest_Nvc;
connection *DC; /* Connections */

/*! Data structure for non-local connection.
 */
struct foreign_connection
{
  int task;
  int origin;
  int index;
  int image_flags;
} * ForeignDC, *ImportedDC;

#define MASK_X_SHIFT_RIGHT 38347922
#define MASK_X_SHIFT_LEFT 76695844
#define MASK_Y_SHIFT_RIGHT 14708792
#define MASK_Y_SHIFT_LEFT 117670336
#define MASK_Z_SHIFT_RIGHT 261632
#define MASK_Z_SHIFT_LEFT 133955584
#define MASK ((1 << 27) - 1)

int FirstUnusedConnection;

/*! \brief Gets connected active cells from a mesh.
 *
 *  \param[in] T Pointer to tesselation.
 *
 *  \return Number of cells.
 */
int voronoi_get_connected_particles(tessellation *T)
{
  int idx, i, j, p, q, count = 0, duplicates, image_flags, listp, nexport, nimport, origin;
  int ngrp, recvTask;

  CPU_Step[CPU_MISC] += measure_time();

  /* first, let's add all the primary active points */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      p = TimeBinsHydro.ActiveParticleList[idx];
      if(p < 0)
        continue;

      if(P[p].Type == 0)
        {
          Ngb_Marker[p] = Ngb_MarkerValue;

          if(P[p].Mass == 0 && P[p].ID == 0) /* skip cells that have been swallowed or eliminated */
            {
              List_P[p].firstexport   = -1;
              List_P[p].currentexport = -1;
              continue;
            }

          if(Ninlist >= MaxNinlist)
            {
              T->Indi.AllocFacNinlist *= ALLOC_INCREASE_FACTOR;
              MaxNinlist = T->Indi.AllocFacNinlist;
#ifdef VERBOSE
              printf("VORONOI: Task=%d: increase memory allocation, MaxNinlist=%d Indi.AllocFacNinlist=%g\n", ThisTask, MaxNinlist,
                     T->Indi.AllocFacNinlist);
#endif /* #ifdef VERBOSE */
              ListExports = myrealloc_movable(ListExports, MaxNinlist * sizeof(struct list_export_data));

              if(Ninlist >= MaxNinlist)
                terminate("Ninlist >= MaxNinlist");
            }

          List_InMesh[NumGasInMesh++] = p;

          List_P[p].currentexport = List_P[p].firstexport = Ninlist++;
          ListExports[List_P[p].currentexport].image_bits = 1;
          ListExports[List_P[p].currentexport].nextexport = -1;
          ListExports[List_P[p].currentexport].origin     = ThisTask;
          ListExports[List_P[p].currentexport].index      = p;

          if(T->Ndp >= T->MaxNdp)
            {
              T->Indi.AllocFacNdp *= ALLOC_INCREASE_FACTOR;
              T->MaxNdp = T->Indi.AllocFacNdp;
#ifdef VERBOSE
              printf("VORONOI: Task=%d: increase memory allocation, MaxNdp=%d Indi.AllocFacNdp=%g\n", ThisTask, T->MaxNdp,
                     T->Indi.AllocFacNdp);
#endif /* #ifdef VERBOSE */
              T->DP -= 5;
              T->DP = myrealloc_movable(T->DP, (T->MaxNdp + 5) * sizeof(point));
              T->DP += 5;

              if(T->Ndp >= T->MaxNdp)
                terminate("Ndp >= MaxNdp");
            }

          SphP[p].ActiveArea = 0;

          point *dp = &T->DP[T->Ndp];

          dp->x             = P[p].Pos[0];
          dp->y             = P[p].Pos[1];
          dp->z             = P[p].Pos[2];
          dp->ID            = P[p].ID;
          dp->task          = ThisTask;
          dp->index         = p;
          dp->originalindex = -1;
          dp->timebin       = P[p].TimeBinHydro;
          dp->image_flags   = 1;
#ifdef DOUBLE_STENCIL
          dp->Hsml             = SphP[p].Hsml;
          dp->first_connection = -1;
          dp->last_connection  = -1;
#endif /* #ifdef DOUBLE_STENCIL */
          T->Ndp++;
          count++;
        }
    }

  /* now, we go through the connection list and see whether we have any additional points to add */
  int count_foreign = 0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      q = SphP[i].first_connection;

      while(q >= 0)
        {
          if(q < 0 || q >= MaxNvc)
            {
              char buf[1000];
              sprintf(buf, "strange connectivity q=%d Nvc=%d", q, MaxNvc);
              terminate(buf);
            }

          if(DC[q].task >= 0 && DC[q].task < NTask)
            {
              if(ThisTask == DC[q].task) /* this one is local */
                {
                  p = DC[q].index; /* particle index */

                  if(P[p].Type == 0)
                    {
                      if(!(P[p].Mass == 0 && P[p].ID == 0)) /* skip cells that have been swallowed or dissolved */
                        {
                          if(P[p].Ti_Current != All.Ti_Current)
                            {
                              drift_particle(p, All.Ti_Current);
                            }

                          if(p < 0 || p >= NumGas)
                            {
                              char buf[1000];
                              sprintf(buf, "strange p=%d (Ngas=%d) for q=%d Nvc=%d", p, NumGas, q, Nvc);
                              terminate(buf);
                            }

                          image_flags = (DC[q].image_flags & MASK);

                          if(Ngb_Marker[p] != Ngb_MarkerValue)
                            {
                              Ngb_Marker[p]           = Ngb_MarkerValue;
                              List_P[p].firstexport   = -1;
                              List_P[p].currentexport = -1;
                            }

                          listp = List_P[p].firstexport;

                          /* now we need to check whether this particle has already been made part of the list */
                          if(List_P[p].firstexport >= 0)
                            {
                              if(ListExports[List_P[p].currentexport].origin != ThisTask)
                                terminate("can't be");
                            }
                          else
                            {
                              /* this one apparently hasn't been added at all yet */
                              if(Ninlist >= MaxNinlist)
                                {
                                  T->Indi.AllocFacNinlist *= ALLOC_INCREASE_FACTOR;
                                  MaxNinlist = T->Indi.AllocFacNinlist;
#ifdef VERBOSE
                                  printf("Task=%d: increase memory allocation, MaxNinlist=%d Indi.AllocFacNinlist=%g\n", ThisTask,
                                         MaxNinlist, T->Indi.AllocFacNinlist);
#endif /* #ifdef VERBOSE */
                                  ListExports = myrealloc_movable(ListExports, MaxNinlist * sizeof(struct list_export_data));

                                  if(Ninlist >= MaxNinlist)
                                    terminate("Ninlist >= MaxNinlist");
                                }

                              List_InMesh[NumGasInMesh++] = p;

                              List_P[p].currentexport = List_P[p].firstexport = Ninlist++;
                              ListExports[List_P[p].currentexport].image_bits = 0;
                              ListExports[List_P[p].currentexport].nextexport = -1;
                              ListExports[List_P[p].currentexport].origin     = ThisTask;
                              ListExports[List_P[p].currentexport].index      = p;
                            }

                          if(!(ListExports[List_P[p].currentexport].image_bits & image_flags)) /* already in list */
                            {
                              ListExports[List_P[p].currentexport].image_bits |= image_flags;

                              if(T->Ndp >= T->MaxNdp)
                                {
                                  T->Indi.AllocFacNdp *= ALLOC_INCREASE_FACTOR;
                                  T->MaxNdp = T->Indi.AllocFacNdp;
#ifdef VERBOSE
                                  printf("Task=%d: increase memory allocation, MaxNdp=%d Indi.AllocFacNdp=%g\n", ThisTask, T->MaxNdp,
                                         T->Indi.AllocFacNdp);
#endif /* #ifdef VERBOSE */
                                  T->DP -= 5;
                                  T->DP = myrealloc_movable(T->DP, (T->MaxNdp + 5) * sizeof(point));
                                  T->DP += 5;

                                  if(T->Ndp >= T->MaxNdp)
                                    terminate("Ndp >= MaxNdp");
                                }

                              SphP[p].ActiveArea = 0;

                              MyDouble x = P[p].Pos[0];
                              MyDouble y = P[p].Pos[1];
                              MyDouble z = P[p].Pos[2];

                              /* for each coordinates there are three possibilities. They are encoded in image_flag to basis three,
                               * i.e. x*3^0 + y*3^1 + z*3^2 */

#ifndef REFLECTIVE_X
                              if((image_flags & MASK_X_SHIFT_RIGHT))
                                x += boxSize_X;
                              else if((image_flags & MASK_X_SHIFT_LEFT))
                                x -= boxSize_X;
#else  /* #ifndef REFLECTIVE_X */
                              if((image_flags & MASK_X_SHIFT_RIGHT))
                                x = -x;
                              else if((image_flags & MASK_X_SHIFT_LEFT))
                                x = 2 * boxSize_X - x;
#endif /* #ifndef REFLECTIVE_X #else */
#ifndef REFLECTIVE_Y
                              if((image_flags & MASK_Y_SHIFT_RIGHT))
                                y += boxSize_Y;
                              else if((image_flags & MASK_Y_SHIFT_LEFT))
                                y -= boxSize_Y;
#else  /* #ifndef REFLECTIVE_Y */
                              if((image_flags & MASK_Y_SHIFT_RIGHT))
                                y = -y;
                              else if((image_flags & MASK_Y_SHIFT_LEFT))
                                y = 2 * boxSize_Y - y;
#endif /* #ifndef REFLECTIVE_Y #else */
#ifndef REFLECTIVE_Z
                              if((image_flags & MASK_Z_SHIFT_RIGHT))
                                z += boxSize_Z;
                              else if((image_flags & MASK_Z_SHIFT_LEFT))
                                z -= boxSize_Z;
#else  /* #ifndef REFLECTIVE_Z */
                              if((image_flags & MASK_Z_SHIFT_RIGHT))
                                z = -z;
                              else if((image_flags & MASK_Z_SHIFT_LEFT))
                                z = 2 * boxSize_Z - z;
#endif /* #ifndef REFLECTIVE_Z #else */

                              point *dp = &T->DP[T->Ndp];

                              dp->x = x;
                              dp->y = y;
                              dp->z = z;

                              dp->task = ThisTask;
                              dp->ID   = P[p].ID;
                              if(image_flags != 1)
                                dp->index = p + NumGas; /* this is a replicated/mirrored local point */
                              else
                                dp->index = p; /* this is actually a local point that wasn't made part of the mesh yet */
                              dp->originalindex = p;
                              dp->timebin       = P[p].TimeBinHydro;

                              dp->image_flags = image_flags;
#ifdef DOUBLE_STENCIL
                              dp->Hsml             = SphP[p].Hsml;
                              dp->first_connection = -1;
                              dp->last_connection  = -1;
#endif /* #ifdef DOUBLE_STENCIL */
                              T->Ndp++;
                              count++;
                            }
                        }
                    }
                }
              else
                {
                  /* here we have a foreign neighbor that we want */
                  count_foreign++;
                }
            }

          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }
    }

  /* we now compile a list of the foreign neighbors we want in the mesh */

  ForeignDC = mymalloc_movable(&ForeignDC, "ForeignDC", count_foreign * sizeof(struct foreign_connection));

  int count_foreign_bak = count_foreign;

  count_foreign = 0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      q = SphP[i].first_connection;

      while(q >= 0)
        {
          if(DC[q].task >= 0 && DC[q].task < NTask)
            {
              if(ThisTask != DC[q].task) /* this one is not local */
                {
                  p = DC[q].index;

                  ForeignDC[count_foreign].task        = DC[q].task;
                  ForeignDC[count_foreign].origin      = ThisTask;
                  ForeignDC[count_foreign].index       = DC[q].index;
                  ForeignDC[count_foreign].image_flags = (DC[q].image_flags & MASK);

                  /* here we have a foreign neighbor that we want */
                  count_foreign++;
                }
            }

          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }
    }

  if(count_foreign_bak != count_foreign)
    terminate("bad");

  /* we sort this list by tasks, and then eliminate duplicates */
  mysort(ForeignDC, count_foreign, sizeof(struct foreign_connection), compare_foreign_connection);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0, j = -1, duplicates = 0; i < count_foreign; i++)
    {
      if(j >= 0)
        if(memcmp(&ForeignDC[i], &ForeignDC[j], sizeof(struct foreign_connection)) == 0)
          {
            duplicates++;
            continue;
          }

      j++;

      ForeignDC[j] = ForeignDC[i];
      Send_count[ForeignDC[j].task]++;
    }

  count_foreign -= duplicates;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, nexport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  if(nexport != count_foreign)
    {
      char buf[1000];
      sprintf(buf, "nexport=%d  count_foreign=%d\n", nexport, count_foreign);
      terminate(buf);
    }

  if(Send_count[ThisTask] != 0)
    terminate("bad");

  ImportedDC = mymalloc_movable(&ImportedDC, "ImportedDC", nimport * sizeof(struct foreign_connection));

  /* get the point requests */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              MPI_Sendrecv(&ForeignDC[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct foreign_connection), MPI_BYTE,
                           recvTask, TAG_DENS_B, &ImportedDC[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct foreign_connection), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  point *DP_Buffer = (point *)mymalloc_movable(&DP_Buffer, "DP_Buffer", nimport * sizeof(point));

  /* now we prepare the points */
  for(j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  for(i = 0; i < nimport; i++)
    {
      p           = ImportedDC[i].index;
      origin      = ImportedDC[i].origin;
      image_flags = ImportedDC[i].image_flags;

      /* it could happen that the requested point has been refined or was turned into a star, that's why
       * we not necessarily will find all the points requested.
       */
      if(P[p].Type != 0)
        continue;

      if(P[p].Mass == 0 && P[p].ID == 0)
        continue; /* skip cells that have been swallowed or dissolved */

      if(P[p].Ti_Current != All.Ti_Current)
        {
          drift_particle(p, All.Ti_Current);
        }

      /* mark the points in the export lists */

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

      ListExports[List_P[p].currentexport].image_bits |= image_flags;

      MyDouble x = P[p].Pos[0];
      MyDouble y = P[p].Pos[1];
      MyDouble z = P[p].Pos[2];

      /* for each coordinates there are three possibilities. They are encoded in image_flag to basis three, i.e. x*3^0 + y*3^1 + z*3^2
       */
#ifndef REFLECTIVE_X
      if((image_flags & MASK_X_SHIFT_RIGHT))
        x += boxSize_X;
      else if((image_flags & MASK_X_SHIFT_LEFT))
        x -= boxSize_X;
#else  /* #ifndef REFLECTIVE_X */
      if((image_flags & MASK_X_SHIFT_RIGHT))
        x = -x;
      else if((image_flags & MASK_X_SHIFT_LEFT))
        x = 2 * boxSize_X - x;
#endif /* #ifndef REFLECTIVE_X #else */

#ifndef REFLECTIVE_Y
      if((image_flags & MASK_Y_SHIFT_RIGHT))
        y += boxSize_Y;
      else if((image_flags & MASK_Y_SHIFT_LEFT))
        y -= boxSize_Y;
#else  /* #ifndef REFLECTIVE_Y */
      if((image_flags & MASK_Y_SHIFT_RIGHT))
        y = -y;
      else if((image_flags & MASK_Y_SHIFT_LEFT))
        y = 2 * boxSize_Y - y;
#endif /* #ifndef REFLECTIVE_Y #else */

#ifndef REFLECTIVE_Z
      if((image_flags & MASK_Z_SHIFT_RIGHT))
        z += boxSize_Z;
      else if((image_flags & MASK_Z_SHIFT_LEFT))
        z -= boxSize_Z;
#else  /* #ifndef REFLECTIVE_Z */
      if((image_flags & MASK_Z_SHIFT_RIGHT))
        z = -z;
      else if((image_flags & MASK_Z_SHIFT_LEFT))
        z = 2 * boxSize_Z - z;
#endif /* #ifndef REFLECTIVE_Z #else */

      int k = Recv_offset[origin] + Recv_count[origin]++;

      SphP[p].ActiveArea = 0;

      DP_Buffer[k].x             = x;
      DP_Buffer[k].y             = y;
      DP_Buffer[k].z             = z;
      DP_Buffer[k].ID            = P[p].ID;
      DP_Buffer[k].task          = ThisTask;
      DP_Buffer[k].index         = p;
      DP_Buffer[k].originalindex = p;
      DP_Buffer[k].timebin       = P[p].TimeBinHydro;

      DP_Buffer[k].image_flags = image_flags;
#ifdef DOUBLE_STENCIL
      DP_Buffer[k].Hsml             = SphP[p].Hsml;
      DP_Buffer[k].first_connection = -1;
      DP_Buffer[k].last_connection  = -1;
#endif /* #ifdef DOUBLE_STENCIL */
    }

  /* because we may have dropped some of the points because they were turned
   * into stars we need to redetermine the send-offsets and counts
   */

  MPI_Alltoall(Recv_count, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, nexport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          /* note: the Recv_offsets stay at this point */
        }
    }

  /* now get the additional Delaunay points from the other processors */

  while(nexport + T->Ndp > T->MaxNdp)
    {
      T->Indi.AllocFacNdp *= ALLOC_INCREASE_FACTOR;
      T->MaxNdp = T->Indi.AllocFacNdp;
#ifdef VERBOSE
      printf("Task=%d: increase memory allocation, MaxNdp=%d Indi.AllocFacNdp=%g\n", ThisTask, T->MaxNdp, T->Indi.AllocFacNdp);
#endif /* #ifdef VERBOSE */
      T->DP -= 5;
      T->DP = myrealloc_movable(T->DP, (T->MaxNdp + 5) * sizeof(point));
      T->DP += 5;
    }

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the Delaunay points */

              MPI_Sendrecv(&DP_Buffer[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(point), MPI_BYTE, recvTask, TAG_DENS_B,
                           &T->DP[T->Ndp + Send_offset[recvTask]], Send_count[recvTask] * sizeof(point), MPI_BYTE, recvTask,
                           TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  T->Ndp += nexport;
  count += nexport;

  myfree(DP_Buffer);
  myfree(ImportedDC);
  myfree(ForeignDC);

  mpi_printf("VORONOI: done with connected particles\n");

  CPU_Step[CPU_MESH_DYNAMIC] += measure_time();

  /* at this point, it might make sense to sort the Delaunay point again
   * according to Peano-Hilbert, in an extended region that allows for the
   * ghost regions
   */

  peano_hilbert_order_DP();

  CPU_Step[CPU_PEANO] += measure_time();

  return count;
}

/*! \brief Initialises connectivity.
 *
 *  \param[in] T Pointer to tessellation.
 *
 *  \return void
 */
void voronoi_init_connectivity(tessellation *T)
{
  int i;

  mpi_printf("VORONOI: init connectivity\n");

  MaxNvc = T->Indi.AllocFacNvc;
  DC     = mymalloc_movable(&DC, "DC", MaxNvc * sizeof(connection));

  Nvc = 0;

  /* we use a chaining list to keep track of unused entries in the list of connections */
  /* here we set it up to contain all available spaces */
  FirstUnusedConnection = 0;
  for(i = 0; i < MaxNvc - 1; i++)
    {
      DC[i].next = i + 1;
      DC[i].task = -1; /* mark that this is unused */
    }
  DC[MaxNvc - 1].next = -1;
  DC[MaxNvc - 1].task = -1;

  /* initially, all particle have empty connection lists */
  for(i = 0; i < NumGas; i++)
    SphP[i].first_connection = SphP[i].last_connection = -1;

  mpi_printf("VORONOI: done with init of connectivity\n");
}

/*! \brief Updates connectivity.
 *
 *  \param[in] T Pointer to tessellation.
 *
 *  \return void
 */
void voronoi_update_connectivity(tessellation *T)
{
  int idx, i, k, q, p_task, p_index, q_task, q_index, q_dp_index, q_image_flags;
  MyIDType p_ID;

  CPU_Step[CPU_MISC] += measure_time();

  /* let's clear the connection lists of active particles */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(i >= NumGas)
        terminate("i >= NumGas");

      q = SphP[i].first_connection;

      if(q >= 0) /* we have connections, let's add them to the free list */
        {
          while(q >= 0)
            {
              Nvc--;
              DC[q].task = -1; /* mark that this is unused */

              if(q == SphP[i].last_connection)
                break;

              q = DC[q].next;
            }

          /* we add the new free spots at the beginning of the free list */
          DC[SphP[i].last_connection].next = FirstUnusedConnection;
          FirstUnusedConnection            = SphP[i].first_connection;

          SphP[i].first_connection = -1;
          SphP[i].last_connection  = -1;
        }
    }

  for(i = 0; i < T->Nvf; i++)
    {
      for(k = 0; k < 2; k++)
        {
          point *DP = T->DP;
          face *VF  = T->VF;

          if(k == 0)
            {
              p_task        = DP[VF[i].p1].task;
              p_index       = DP[VF[i].p1].index;
              p_ID          = DP[VF[i].p1].ID;
              q_task        = DP[VF[i].p2].task;
              q_index       = DP[VF[i].p2].index;
              q_dp_index    = VF[i].p2;
              q_image_flags = (DP[VF[i].p2].image_flags & MASK);
            }
          else
            {
              p_task        = DP[VF[i].p2].task;
              p_index       = DP[VF[i].p2].index;
              p_ID          = DP[VF[i].p2].ID;
              q_task        = DP[VF[i].p1].task;
              q_index       = DP[VF[i].p1].index;
              q_dp_index    = VF[i].p1;
              q_image_flags = (DP[VF[i].p1].image_flags & MASK);
            }

          if(p_task == ThisTask && p_index >= 0 && p_index < NumGas)
            {
              if(TimeBinSynchronized[P[p_index].TimeBinHydro])
                {
                  if(P[p_index].Type != 0)
                    continue;

                  if(P[p_index].Mass == 0 && P[p_index].ID == 0)
                    continue; /* skip cells that have been swallowed or dissolved */

                  /* need to add the connection to the other point to this particle */

                  if(FirstUnusedConnection < 0 || Nvc == MaxNvc)
                    {
                      if(!(FirstUnusedConnection < 0 && Nvc == MaxNvc))
                        {
                          char buf[1000];
                          sprintf(buf, "strange: FirstUnusedConnection=%d Nvc=%d MaxNvc=%d\n", FirstUnusedConnection, Nvc, MaxNvc);
                          terminate(buf);
                        }

                      int n, old_MaxNvc = MaxNvc;
                      T->Indi.AllocFacNvc *= ALLOC_INCREASE_FACTOR;
                      MaxNvc = T->Indi.AllocFacNvc;
#ifdef VERBOSE
                      printf("Task=%d: increase memory allocation, MaxNvc=%d Indi.AllocFacNvc=%g\n", ThisTask, MaxNvc,
                             T->Indi.AllocFacNvc);
#endif /* #ifdef VERBOSE */
                      DC = myrealloc_movable(DC, MaxNvc * sizeof(connection));
                      DP = T->DP;
                      VF = T->VF;

                      FirstUnusedConnection = old_MaxNvc;
                      for(n = old_MaxNvc; n < MaxNvc - 1; n++)
                        {
                          DC[n].next = n + 1;
                          DC[n].task = -1;
                        }
                      DC[MaxNvc - 1].next = -1;
                      DC[MaxNvc - 1].task = -1;
                    }

                  if(SphP[p_index].last_connection >= 0)
                    {
                      DC[SphP[p_index].last_connection].next = FirstUnusedConnection;
                      SphP[p_index].last_connection          = FirstUnusedConnection;
                    }
                  else
                    {
                      SphP[p_index].last_connection  = FirstUnusedConnection;
                      SphP[p_index].first_connection = FirstUnusedConnection;
                    }

                  FirstUnusedConnection = DC[FirstUnusedConnection].next;
                  Nvc++;

                  DC[SphP[p_index].last_connection].task        = q_task;
                  DC[SphP[p_index].last_connection].image_flags = q_image_flags;
                  DC[SphP[p_index].last_connection].ID          = p_ID;

                  if(q_task == ThisTask && q_index >= NumGas)
                    DC[SphP[p_index].last_connection].index = q_index - NumGas;
                  else
                    DC[SphP[p_index].last_connection].index = q_index;

                  DC[SphP[p_index].last_connection].dp_index = q_dp_index;
#ifdef TETRA_INDEX_IN_FACE
                  DC[SphP[p_index].last_connection].dt_index = VF[i].dt_index;
#endif                                                            /* #ifdef TETRA_INDEX_IN_FACE */
                  DC[SphP[p_index].last_connection].vf_index = i; /* index to the corresponding face */

                  if(SphP[p_index].last_connection >= MaxNvc)
                    {
                      terminate("this is wrong");
                    }
                }
            }

#ifdef DOUBLE_STENCIL
          int index;
          if(k == 0)
            index = VF[i].p1;
          else
            index = VF[i].p2;

          if(!(p_task == ThisTask && p_index >= 0 && p_index < NumGas) && DP[index].flag_primary_triangle > 0 && index >= 0)
            {
              /* need to add the connection to the other point to this particle */

              if(FirstUnusedConnection < 0 || Nvc == MaxNvc)
                {
                  if(!(FirstUnusedConnection < 0 && Nvc == MaxNvc))
                    {
                      char buf[1000];
                      sprintf(buf, "strange: FirstUnusedConnection=%d Nvc=%d MaxNvc=%d\n", FirstUnusedConnection, Nvc, MaxNvc);
                      terminate(buf);
                    }

                  int n, old_MaxNvc = MaxNvc;
                  T->Indi.AllocFacNvc *= ALLOC_INCREASE_FACTOR;
                  MaxNvc = T->Indi.AllocFacNvc;
#ifdef VERBOSE
                  printf("Task=%d: increase memory allocation, MaxNvc=%d Indi.AllocFacNvc=%g\n", ThisTask, MaxNvc,
                         T->Indi.AllocFacNvc);
#endif /* #ifdef VERBOSE */
                  DC = myrealloc_movable(DC, MaxNvc * sizeof(connection));
                  DP = T->DP;
                  VF = T->VF;

                  FirstUnusedConnection = old_MaxNvc;
                  for(n = old_MaxNvc; n < MaxNvc - 1; n++)
                    {
                      DC[n].next = n + 1;
                      DC[n].task = -1;
                    }
                  DC[MaxNvc - 1].next = -1;
                  DC[MaxNvc - 1].task = -1;
                }

              if(DP[index].last_connection >= 0)
                {
                  DC[DP[index].last_connection].next = FirstUnusedConnection;
                  DP[index].last_connection          = FirstUnusedConnection;
                }
              else
                {
                  DP[index].last_connection  = FirstUnusedConnection;
                  DP[index].first_connection = FirstUnusedConnection;
                }

              FirstUnusedConnection = DC[FirstUnusedConnection].next;
              Nvc++;

              DC[DP[index].last_connection].task        = q_task;
              DC[DP[index].last_connection].image_flags = q_image_flags;
              DC[DP[index].last_connection].ID          = p_ID;

              if(q_task == ThisTask && q_index >= NumGas)
                DC[DP[index].last_connection].index = q_index - NumGas;
              else
                DC[DP[index].last_connection].index = q_index;

              DC[DP[index].last_connection].dp_index = q_dp_index;

              DC[DP[index].last_connection].vf_index = i; /* index to the corresponding face */

              if(DP[index].last_connection >= MaxNvc)
                {
                  terminate("this is wrong");
                }
            }
#endif /* #ifdef DOUBLE_STENCIL */
        }
    }

  mpi_printf("VORONOI: done with updating connectivity.\n");

  CPU_Step[CPU_MESH_DYNAMIC] += measure_time();
}

/*! \brief Remove connection from cell.
 *
 *  \param[in] i Index of cell.
 *
 *  \return void
 */
void voronoi_remove_connection(int i)
{
  int q;
  if((q = SphP[i].first_connection) >= 0) /* we have connections, let's add them to the free list */
    {
      while(q >= 0)
        {
          Nvc--;
          DC[q].task = -1; /* mark that this is unused */

          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }

      /* we add the new free spots at the beginning of the free list */
      DC[SphP[i].last_connection].next = FirstUnusedConnection;
      FirstUnusedConnection            = SphP[i].first_connection;

      SphP[i].first_connection = -1;
      SphP[i].last_connection  = -1;
    }
}

/*! \brief Compares two foreign connection objects.
 *
 *  Compares (highest priority first):
 *      task
 *      index
 *      image_flags
 *
 *  \param[in] a First object.
 *  \param[in] b Second object.
 *
 *  \return (-1,0,1); -1: a < b.
 */
int compare_foreign_connection(const void *a, const void *b)
{
  if(((struct foreign_connection *)a)->task < (((struct foreign_connection *)b)->task))
    return -1;

  if(((struct foreign_connection *)a)->task > (((struct foreign_connection *)b)->task))
    return +1;

  if(((struct foreign_connection *)a)->index < (((struct foreign_connection *)b)->index))
    return -1;

  if(((struct foreign_connection *)a)->index > (((struct foreign_connection *)b)->index))
    return +1;

  if(((struct foreign_connection *)a)->image_flags < (((struct foreign_connection *)b)->image_flags))
    return -1;

  if(((struct foreign_connection *)a)->image_flags > (((struct foreign_connection *)b)->image_flags))
    return +1;

  return 0;
}
