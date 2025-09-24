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
 * \file        src/mesh/voronoi/voronoi.c
 * \date        05/2018
 * \brief       Main file for Voronoi-mesh construction.
 * \details     contains functions:
 *                void create_mesh(void)
 *                int voronoi_get_local_particles(void)
 *                void free_mesh_structures_not_needed_for_derefinement_
 *                  refinement(void)
 *                void free_all_remaining_mesh_structures(void)
 *                void free_mesh(void)
 *                int compute_max_delaunay_radius(void)
 *                void compute_voronoi_faces_and_volumes(void)
 *                int area_list_data_compare(const void *a, const void *b)
 *                void apply_area_list(void)
 *                void derefine_refine_compute_volumes(double *vol)
 *                double nearest_x(double d)
 *                double nearest_y(double d)
 *                double nearest_z(double d)
 *                double get_cell_radius(int i)
 *                void dump_points(tessellation * T)
 *                int face_get_normals(tessellation * T, int i, struct
 *                  geometry *geom)
 *                double distance_to_border(int cell)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 21.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../main/allvars.h"
#include "../../main/proto.h"

#include "voronoi.h"

tessellation Mesh, DeRefMesh;

unsigned char *Edge_visited;
struct area_list_data *AreaList;
int Narea, MaxNarea;

int DPinfinity; /* marker for special infinity point */
double CentralOffsetX, CentralOffsetY, CentralOffsetZ, ConversionFac;

struct list_export_data *ListExports;
struct list_P_data *List_P;
int NumGasInMesh;
int *List_InMesh;

int CountInSphereTests, CountInSphereTestsExact;
int CountConvexEdgeTest, CountConvexEdgeTestExact;
int Ninlist, MaxNinlist;

int CountFlips, Count_1_to_3_Flips2d, Count_2_to_4_Flips2d;
int Count_1_to_4_Flips, Count_2_to_3_Flips, Count_3_to_2_Flips, Count_4_to_4_Flips;
int Count_EdgeSplits, Count_FaceSplits;
int Count_InTetra, Count_InTetraExact;
int Largest_N_DP_Buffer;

long long TotCountInSphereTests, TotCountInSphereTestsExact;
long long TotCountConvexEdgeTest, TotCountConvexEdgeTestExact;

long long TotCountFlips, TotCount_1_to_3_Flips2d, TotCount_2_to_4_Flips2d;
long long TotCount_1_to_4_Flips, TotCount_2_to_3_Flips, TotCount_3_to_2_Flips, TotCount_4_to_4_Flips;
long long TotCount_EdgeSplits, TotCount_FaceSplits;
long long TotCount_InTetra, TotCount_InTetraExact;

/*! \brief Creates the Voronoi mesh.
 *
 *  Routine which is called in run.
 *  If first creates a first, giant tetrahedron and than successively insert
 *  particles (first local, then ghost particles) compute their circumcircles
 *  and count the undecided tetrahedra. This procedure is repeated until all
 *  tetrahedra are decided. Then, the maximum Delauny radius is computed as
 *  well as the faces and volumes of the Voronoi-cells.
 *
 *  \return void
 */
void create_mesh(void)
{
#ifdef CREATE_FULL_MESH
  int k;

  short int *buTimeBin = mymalloc_movable(&buTimeBin, "buTimeBin", NumPart * sizeof(short int));
  static int buTimeBinActive[TIMEBINS];

  for(k = 0; k < NumPart; k++)
    {
      buTimeBin[k]      = P[k].TimeBinHydro;
      P[k].TimeBinHydro = 0;
    }

  for(k = 0; k < TIMEBINS; k++)
    {
      buTimeBinActive[k] = TimeBinSynchronized[k];

      TimeBinSynchronized[k] = 1;
    }

  reconstruct_timebins();
#endif /* #ifdef CREATE_FULL_MESH */

  int tlast;
  int idx, i, iter = 0, n, skip;
  double tstart, tend;
  long long ntot;

  if(All.TotNumGas == 0)
    return;

  TIMER_START(CPU_MESH);

  mpi_printf("VORONOI: create delaunay mesh\n");

  Ngb_MarkerValue++;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Ti_Current != All.Ti_Current)
        {
          terminate("surprise! we don't expect this here anymore");
          drift_particle(i, All.Ti_Current);
        }

      SphP[i].Hsml = 1.01 * SphP[i].MaxDelaunayRadius;
    }

  initialize_and_create_first_tetra(&Mesh);

  CountInSphereTests = CountInSphereTestsExact = 0;
  CountConvexEdgeTest = CountConvexEdgeTestExact = 0;
  CountFlips = Count_1_to_3_Flips2d = Count_2_to_4_Flips2d = 0;
  Count_1_to_4_Flips                                       = 0;
  Count_2_to_3_Flips                                       = 0;
  Count_3_to_2_Flips                                       = 0;
  Count_4_to_4_Flips                                       = 0;
  Count_EdgeSplits                                         = 0;
  Count_FaceSplits                                         = 0;
  Count_InTetra = Count_InTetraExact = 0;
  Largest_N_DP_Buffer                = 0;

  MaxNinlist  = Mesh.Indi.AllocFacNinlist;
  ListExports = mymalloc_movable(&ListExports, "ListExports", MaxNinlist * sizeof(struct list_export_data));

  NumGasInMesh = 0;
  List_InMesh  = mymalloc_movable(&List_InMesh, "List_InMesh", NumGas * sizeof(int));

  List_P = mymalloc_movable(&List_P, "List_P", NumGas * sizeof(struct list_P_data));

  Mesh.DTC = mymalloc_movable(&Mesh.DTC, "DTC", Mesh.MaxNdt * sizeof(tetra_center));
  Mesh.DTF = mymalloc_movable(&Mesh.DTF, "DTF", Mesh.MaxNdt * sizeof(char));
  for(i = 0; i < Mesh.Ndt; i++)
    Mesh.DTF[i] = 0;

  Ninlist = 0;

  tlast = 0;

  do
    {
      skip = Mesh.Ndp;

      TIMER_STOPSTART(CPU_MESH, CPU_MESH_FIND_DP);

      tstart = second();

      if(iter == 0)
        {
          MPI_Allreduce(&Nvc, &Largest_Nvc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

          if(Largest_Nvc > 0)
            n = voronoi_get_connected_particles(&Mesh);
          else
            n = voronoi_get_local_particles();
        }
      else
        {
          n = voronoi_ghost_search(&Mesh);
        }

      sumup_large_ints(1, &n, &ntot);

      tend = second();

      if(iter == 0)
        mpi_printf("VORONOI: iter=%d: %llu local points, points/sec/task = %g, took %g secs\n", iter, ntot,
                   ntot / (timediff(tstart, tend) + 1.0e-30) / NTask, timediff(tstart, tend));
      else
        {
          if(ntot)
            mpi_printf("VORONOI: iter=%d: %llu additional points, points/sec/task = %g, took %g secs\n", iter, ntot,
                       ntot / (timediff(tstart, tend) + 1.0e-30) / NTask, timediff(tstart, tend));
          else
            mpi_printf("VORONOI: iter=%d: %llu additional points, took %g secs\n", iter, ntot, timediff(tstart, tend));
        }

      TIMER_STOPSTART(CPU_MESH_FIND_DP, CPU_MESH_INSERT);

      for(i = 0; i < n; i++)
        {
#ifndef OPTIMIZE_MEMORY_USAGE
          set_integers_for_point(&Mesh, skip + i);
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE */
          tlast = insert_point(&Mesh, skip + i, tlast);
        }

      TIMER_STOPSTART(CPU_MESH_INSERT, CPU_MESH_CELLCHECK);

      compute_circumcircles(&Mesh);

      if(iter > 0)
        {
          n = count_undecided_tetras(&Mesh);

          sumup_large_ints(1, &n, &ntot);

          if(ntot)
            {
              mpi_printf("VORONOI: still undecided %llu tetrahedras\n", ntot);

#ifndef DOUBLE_STENCIL
              for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
                {
                  i = TimeBinsHydro.ActiveParticleList[idx];
                  if(i < 0)
                    continue;
                  SphP[i].Hsml *= HSML_INCREASE_FACTOR;
                }
#else  /* #ifndef DOUBLE_STENCIL */
              for(i = 0; i < Mesh.Ndp; i++)
                Mesh.DP[i].Hsml *= HSML_INCREASE_FACTOR;
#endif /* #ifndef DOUBLE_STENCIL #else */
            }
        }
      else
        {
          ntot = 1;
        }

      TIMER_STOPSTART(CPU_MESH_CELLCHECK, CPU_MESH);

      if(iter > MAX_VORONOI_ITERATIONS)
        terminate("too many iterations\n");

      iter++;
    }
  while(ntot > 0);

#if(REFLECTIVE_X == 2) || (REFLECTIVE_Y == 2) || (REFLECTIVE_Z == 2)
  for(i = 0; i < Mesh.Ndp; i++)
    {
#if(REFLECTIVE_X == 2)
      Mesh.DP[i].image_flags |= OUTFLOW_X;
#endif /* #if (REFLECTIVE_X == 2) */
#if(REFLECTIVE_Y == 2)
      Mesh.DP[i].image_flags |= OUTFLOW_Y;
#endif /* #if (REFLECTIVE_Y == 2) */
#if(REFLECTIVE_Z == 2)
      Mesh.DP[i].image_flags |= OUTFLOW_Z;
#endif /* #if (REFLECTIVE_Z == 2) */
    }
#endif /* #if (REFLECTIVE_X == 2) || (REFLECTIVE_Y == 2) || (REFLECTIVE_Z == 2) */

  compute_max_delaunay_radius();

  TIMER_STOPSTART(CPU_MESH, CPU_LOGS);

#ifdef VERBOSE
  long long TotNdp, TotNdt;

  int in[15];
  long long out[15];

  in[0] = Mesh.Ndp;
  in[1] = Mesh.Ndt;
  in[2] = CountInSphereTests;
  in[3] = CountInSphereTestsExact;
  in[4] = CountFlips;
  in[5] = Count_InTetra;
  in[6] = Count_InTetraExact;
#ifndef TWODIMS
  in[7]  = Count_1_to_4_Flips;
  in[8]  = Count_2_to_3_Flips;
  in[9]  = Count_3_to_2_Flips;
  in[10] = Count_4_to_4_Flips;
  in[11] = Count_FaceSplits;
  in[12] = Count_EdgeSplits;
  in[13] = CountConvexEdgeTest;
  in[14] = CountConvexEdgeTestExact;
  n      = 15;
#else  /* #ifndef TWODIMS */
  in[7]                   = Count_1_to_3_Flips2d;
  in[8]                   = Count_2_to_4_Flips2d;
  n                       = 9;
#endif /* #ifndef TWODIMS #else */

  sumup_large_ints(n, in, out);

  TotNdp                     = out[0];
  TotNdt                     = out[1];
  TotCountInSphereTests      = out[2];
  TotCountInSphereTestsExact = out[3];
  TotCountFlips              = out[4];
  TotCount_InTetra           = out[5];
  TotCount_InTetraExact      = out[6];
#ifndef TWODIMS
  TotCount_1_to_4_Flips       = out[7];
  TotCount_2_to_3_Flips       = out[8];
  TotCount_3_to_2_Flips       = out[9];
  TotCount_4_to_4_Flips       = out[10];
  TotCount_FaceSplits         = out[11];
  TotCount_EdgeSplits         = out[12];
  TotCountConvexEdgeTest      = out[13];
  TotCountConvexEdgeTestExact = out[14];
#else  /* #ifndef TWODIMS */
  TotCount_1_to_3_Flips2d = out[7];
  TotCount_2_to_4_Flips2d = out[8];
#endif /* #ifndef TWODIMS #else */

  if(ThisTask == 0)
    {
#ifndef TWODIMS
      printf(
          "VORONOI: Average D-Points=%llu  (NumGas=%llu)  D-Tetrahedra=%llu  InSphereTests=%llu  InSphereTestsExact=%llu  "
          "Flips=%llu\n",
          TotNdp / NTask, All.TotNumGas / NTask, TotNdt / NTask, TotCountInSphereTests / NTask, TotCountInSphereTestsExact / NTask,
          TotCountFlips / NTask);
      printf("VORONOI: 1_to_4_Flips=%llu  2_to_3_Flips=%llu  3_to_2_Flips=%llu  4_to_4_Flips=%llu  FaceSplits=%llu  EdgeSplits=%llu\n",
             TotCount_1_to_4_Flips / NTask, TotCount_2_to_3_Flips / NTask, TotCount_3_to_2_Flips / NTask,
             TotCount_4_to_4_Flips / NTask, TotCount_FaceSplits / NTask, TotCount_EdgeSplits / NTask);
      printf("VORONOI: InTetra=%llu  InTetraExact=%llu  ConvexEdgeTest=%llu  ConvexEdgeTestExact=%llu\n", TotCount_InTetra,
             TotCount_InTetraExact / NTask, TotCountConvexEdgeTest / NTask, TotCountConvexEdgeTestExact / NTask);
#else  /* #ifndef TWODIMS */
      printf(
          "VORONOI: Average D-Points=%llu  (NumGas=%llu)  D-Triangles=%llu  InCircleTests=%llu InCircleTestsExact=%llu  Flips=%llu\n",
          TotNdp / NTask, All.TotNumGas / NTask, TotNdt / NTask, TotCountInSphereTests / NTask, TotCountInSphereTestsExact / NTask,
          TotCountFlips / NTask);
      printf("VORONOI: 1_to_3_Flips=%llu  2_to_4_Flips=%llu  InTriangle=%llu  InTriangleExact=%llu\n", TotCount_1_to_3_Flips2d / NTask,
             TotCount_2_to_4_Flips2d / NTask, TotCount_InTetra / NTask, TotCount_InTetraExact / NTask);
#endif /* #ifndef TWODIMS #else */
      printf("VORONOI: Total D-Points: %llu Ratio=%g\n", TotNdp, ((double)TotNdp) / All.TotNumGas);
    }
#endif /* #ifdef VERBOSE */

  TIMER_STOPSTART(CPU_LOGS, CPU_MESH_GEOMETRY);

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      SphP[i].Volume      = 0;
      SphP[i].SurfaceArea = 0;
#if defined(REGULARIZE_MESH_FACE_ANGLE) || defined(OUTPUT_MESH_FACE_ANGLE)
      SphP[i].MaxFaceAngle = 0;
#endif /* #if defined(REGULARIZE_MESH_FACE_ANGLE) || defined(OUTPUT_MESH_FACE_ANGLE) */
#ifdef OUTPUT_SURFACE_AREA
      SphP[i].CountFaces = 0;
#endif /* #ifdef OUTPUT_SURFACE_AREA */
    }

  compute_voronoi_faces_and_volumes();

  double vol, voltot;

  vol = 0;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      vol += SphP[i].Volume;

#ifdef ADAPTIVE_HYDRO_SOFTENING
      P[i].SofteningType = get_softeningtype_for_hydro_cell(i);
#endif /* #ifdef ADAPTIVE_HYDRO_SOFTENING */
    }

  MPI_Reduce(&vol, &voltot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  mpi_printf("VORONOI: Total volume of active cells = %g\n", voltot);

  TIMER_STOP(CPU_MESH_GEOMETRY);

  voronoi_update_connectivity(&Mesh);

  myfree(Mesh.DTF);

  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin) /* only do this for full steps */
    {
      /* check whether we can reduce allocation factors */
      while(Mesh.Ndp < ALLOC_DECREASE_FACTOR * Mesh.Indi.AllocFacNdp && Mesh.Indi.AllocFacNdp > MIN_ALLOC_NUMBER)
        Mesh.Indi.AllocFacNdp /= ALLOC_INCREASE_FACTOR;

      while(Mesh.Ndt < ALLOC_DECREASE_FACTOR * Mesh.Indi.AllocFacNdt && Mesh.Indi.AllocFacNdt > MIN_ALLOC_NUMBER)
        Mesh.Indi.AllocFacNdt /= ALLOC_INCREASE_FACTOR;

      while(Mesh.Nvf < ALLOC_DECREASE_FACTOR * Mesh.Indi.AllocFacNvf && Mesh.Indi.AllocFacNvf > MIN_ALLOC_NUMBER)
        Mesh.Indi.AllocFacNvf /= ALLOC_INCREASE_FACTOR;

      while(Ninlist < ALLOC_DECREASE_FACTOR * Mesh.Indi.AllocFacNinlist && Mesh.Indi.AllocFacNinlist > MIN_ALLOC_NUMBER)
        Mesh.Indi.AllocFacNinlist /= ALLOC_INCREASE_FACTOR;

      while(Largest_N_DP_Buffer < ALLOC_DECREASE_FACTOR * Mesh.Indi.AllocFacN_DP_Buffer &&
            Mesh.Indi.AllocFacN_DP_Buffer > MIN_ALLOC_NUMBER)
        Mesh.Indi.AllocFacN_DP_Buffer /= ALLOC_INCREASE_FACTOR;
    }

#ifdef CREATE_FULL_MESH
  for(k = 0; k < TIMEBINS; k++)
    TimeBinSynchronized[k] = buTimeBinActive[k];

  for(k = 0; k < NumPart; k++)
    P[k].TimeBinHydro = buTimeBin[k];

  reconstruct_timebins();

  myfree_movable(buTimeBin);
#endif /* #if defined(CREATE_FULL_MESH) */
}

/*! \brief Routine that fetches local gas cells.
 *
 *  Runs through all active particles and inserts active gas cells into mesh
 *  structure. Increases length of Mesh.DP and ListExports arrays if needed.
 *
 *  \return Number of points.
 */
int voronoi_get_local_particles(void)
{
  int p, idx, count = 0;

  /* first, let's add all the primary active points */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      p = TimeBinsHydro.ActiveParticleList[idx];

      if(p < 0)
        continue;

      if(P[p].Type == 0)
        {
          Ngb_Marker[p] = Ngb_MarkerValue;

          if((P[p].Mass == 0) && (P[p].ID == 0)) /* skip cells that have been swallowed or eliminated */
            {
              List_P[p].firstexport   = -1;
              List_P[p].currentexport = -1;
              continue;
            }

          if(Ninlist >= MaxNinlist)
            {
              Mesh.Indi.AllocFacNinlist *= ALLOC_INCREASE_FACTOR;
              MaxNinlist = Mesh.Indi.AllocFacNinlist;
#ifdef VERBOSE
              printf("VORONOI: Task=%d: increase memory allocation, MaxNinlist=%d Indi.AllocFacNinlist=%g\n", ThisTask, MaxNinlist,
                     Mesh.Indi.AllocFacNinlist);
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

          if(Mesh.Ndp >= Mesh.MaxNdp)
            {
              Mesh.Indi.AllocFacNdp *= ALLOC_INCREASE_FACTOR;
              Mesh.MaxNdp = Mesh.Indi.AllocFacNdp;
#ifdef VERBOSE
              printf("VORONOI: Task=%d: increase memory allocation, MaxNdp=%d Indi.AllocFacNdp=%g\n", ThisTask, Mesh.MaxNdp,
                     Mesh.Indi.AllocFacNdp);
#endif /* #ifdef VERBOSE */
              Mesh.DP -= 5;
              Mesh.DP = myrealloc_movable(Mesh.DP, (Mesh.MaxNdp + 5) * sizeof(point));
              Mesh.DP += 5;

              if(Mesh.Ndp >= Mesh.MaxNdp)
                terminate("Ndp >= MaxNdp");
            }

          SphP[p].ActiveArea = 0;

          point *dp = &Mesh.DP[Mesh.Ndp];

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

          Mesh.Ndp++;
          count++;
        }
    }

  return count;
}

#ifdef REFINEMENT
struct refdata *RefExch;

/*! \brief Structures that are freed before refinement and derefinement step.
 *
 *  To Optimize the memory usage, this, in comubnation with
 *  free_all_remaining_mesh_structures() can be used instead of a free_mesh()
 *  after the refinement. This saves some memory.
 *
 *  \return void
 */
void free_mesh_structures_not_needed_for_derefinement_refinement(void)
{
  if(All.TotNumGas == 0)
    return;

  int i;

  myfree(GradExch);

  RefExch = (struct refdata *)mymalloc_movable(&RefExch, "RefExch", Mesh_nimport * sizeof(struct refdata));

  for(i = 0; i < Mesh_nimport; i++)
    {
#ifdef REFINEMENT_VOLUME_LIMIT
      RefExch[i].Volume = PrimExch[i].Volume;
#endif /* #ifdef REFINEMENT_VOLUME_LIMIT */
      RefExch[i].TimeBinHydro = PrimExch[i].TimeBinHydro;
    }

  myfree_movable(PrimExch);
}

/* \brief Structures that are freed after refinement and derefinement step.
 *
 *  To Optimize the memory usage, this, in comubnation with
 *  free_mesh_structures_not_needed_for_derefinement_refinement(void) can be
 *  used instead of a free_mesh() after the refinement. This saves some memory.
 *
 *  \return void
 */
void free_all_remaining_mesh_structures(void)
{
  if(All.TotNumGas == 0)
    return;

  myfree(RefExch);

  myfree(Mesh.DTC); /* here we can free the centers of the Delaunay triangles again */
  Mesh.DTC = NULL;
  myfree(List_P);
  myfree(List_InMesh);
  myfree(ListExports);
  myfree(Mesh.DT);
  myfree(Mesh.DP - 5);
  myfree(Mesh.VF);
}
#endif /* #ifdef REFINEMENT */

/*! \brief Frees arrays associated with Voronoi-mesh.
 *
 *  \return void
 */
void free_mesh(void)
{
  if(All.TotNumGas == 0)
    return;

#if defined(DOUBLE_STENCIL)
  mpi_printf("freeing double stencil connections...\n");
  int i;
  for(i = 0; i < Mesh.Ndp; i++)
    if(Mesh.DP[i].first_connection >= 0)
      {
        if(Mesh.DP[i].flag_primary_triangle == 0)
          terminate("Mesh.DP[i].flag_primary_triangle");

        int q = Mesh.DP[i].first_connection;

        if(q >= 0) /* we have connections, let's add them to the free list */
          {
            while(q >= 0)
              {
                Nvc--;
                DC[q].task = -1; /* mark that this is unused */

                if(q == Mesh.DP[i].last_connection)
                  break;

                q = DC[q].next;
              }

            /* we add the new free spots at the beginning of the free list */
            DC[Mesh.DP[i].last_connection].next = FirstUnusedConnection;
            FirstUnusedConnection               = Mesh.DP[i].first_connection;

            Mesh.DP[i].first_connection = -1;
            Mesh.DP[i].last_connection  = -1;
          }
      }
  mpi_printf("done with freeing double stencil connections.\n");
#endif /* #if defined(DOUBLE_STENCIL) */

  myfree_movable(GradExch);
  myfree_movable(PrimExch);

  myfree_movable(Mesh.DTC); /* here we can free the centers of the Delaunay triangles again */
  Mesh.DTC = NULL;
  myfree_movable(List_P);
  myfree_movable(List_InMesh);
  myfree_movable(ListExports);
  myfree_movable(Mesh.DT);
  myfree_movable(Mesh.DP - 5);
  myfree_movable(Mesh.VF);
}

/*! \brief Get the maximum Delaunay radius for all active cells.
 *
 *  Defined as the maximum distance between tetrahedron center and its
 *  neighboring points. Stores this radius in the respective field in the
 *  SphP structure.
 *
 *  \return 0 (unused).
 */
int compute_max_delaunay_radius(void)
{
  int idx, i, j, count = 0;
  point *p;
  double dx, dy, dz, r;

#ifdef ONEDIMS
  return 0;
#endif /* #ifdef ONEDIMS */

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      SphP[i].MaxDelaunayRadius = 0;
    }

  point *DP         = Mesh.DP;
  tetra *DT         = Mesh.DT;
  tetra_center *DTC = Mesh.DTC;

  for(i = 0; i < Mesh.Ndt; i++)
    {
      if(DT[i].t[0] < 0) /* deleted ? */
        continue;

      dx = DP[DT[i].p[0]].x - DTC[i].cx;
      dy = DP[DT[i].p[0]].y - DTC[i].cy;
      dz = DP[DT[i].p[0]].z - DTC[i].cz;

      r = 2 * sqrt(dx * dx + dy * dy + dz * dz);

      for(j = 0; j < (DIMS + 1); j++)
        {
          p = &DP[DT[i].p[j]];

          if(p->task == ThisTask && p->index < NumGas && p->index >= 0)
            if(TimeBinSynchronized[P[p->index].TimeBinHydro])
              if(r > SphP[p->index].MaxDelaunayRadius)
                SphP[p->index].MaxDelaunayRadius = r;
        }
    }

  return count;
}

#ifndef ONEDIMS
/*! \brief Computes interface areas volume of cells.
 *
 *  Loops over Delaunay tetrahedra to calculate interface area and volume
 *  contributions to the individual cells. Calculates as well the center of
 *  mass.
 *
 *  \return void
 */
void compute_voronoi_faces_and_volumes(void)
{
  int idx, i, bit, nr;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      SphP[i].Volume    = 0;
      SphP[i].Center[0] = 0;
      SphP[i].Center[1] = 0;
      SphP[i].Center[2] = 0;
#if defined(REFINEMENT_SPLIT_CELLS)
      SphP[i].MinimumEdgeDistance = MAX_FLOAT_NUMBER;
#endif /* #if defined(REFINEMENT_SPLIT_CELLS) */
    }

  Edge_visited = mymalloc_movable(&Edge_visited, "Edge_visited", Mesh.Ndt * sizeof(unsigned char));

  for(i = 0; i < Mesh.Ndt; i++)
    Edge_visited[i] = 0;

  MaxNarea = Mesh.Indi.AllocFacNflux;
  Narea    = 0;
  AreaList = mymalloc_movable(&AreaList, "AreaList", MaxNarea * sizeof(struct area_list_data));

  for(i = 0; i < Mesh.Ndt; i++)
    {
      if(Mesh.DT[i].t[0] < 0) /* deleted ? */
        continue;

      bit = 1;
      nr  = 0;

      while(Edge_visited[i] != EDGE_ALL)
        {
          if((Edge_visited[i] & bit) == 0)
            process_edge_faces_and_volumes(&Mesh, i, nr);

          bit <<= 1;
          nr++;
        }
    }

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(SphP[i].Volume)
        {
          SphP[i].Center[0] /= SphP[i].Volume;
          SphP[i].Center[1] /= SphP[i].Volume;
          SphP[i].Center[2] /= SphP[i].Volume;
        }
    }

  apply_area_list();
  myfree(AreaList);

  myfree(Edge_visited);
}

/*! \brief Compare task of two area_list_data structures.
 *
 *  \param[in] a Pointer to first area_list_data structure.
 *  \param[in] b Pointer to second area_list_data structure.
 *
 *  \return (-1,0,1), -1 if a.task<b.task.
 */
int area_list_data_compare(const void *a, const void *b)
{
  if(((struct area_list_data *)a)->task < (((struct area_list_data *)b)->task))
    return -1;

  if(((struct area_list_data *)a)->task > (((struct area_list_data *)b)->task))
    return +1;

  return 0;
}

/*! \brief Sorts all interface areas and adds them to respective mesh
 *         generating points (ActiveArea).
 *
 *  \return void
 */
void apply_area_list(void)
{
  int i, j, p, nimport, ngrp, recvTask;

  /* now exchange the area-list and apply where needed */

  mysort(AreaList, Narea, sizeof(struct area_list_data), area_list_data_compare);

  for(j = 0; j < NTask; j++)
    Send_count[j] = 0;

  for(i = 0; i < Narea; i++)
    Send_count[AreaList[i].task]++;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct area_list_data *AreaListGet = (struct area_list_data *)mymalloc("AreaListGet", nimport * sizeof(struct area_list_data));

  /* exchange particle data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&AreaList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct area_list_data), MPI_BYTE, recvTask,
                           TAG_DENS_A, &AreaListGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct area_list_data),
                           MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  /* apply the area */
  for(i = 0; i < nimport; i++)
    {
      p = AreaListGet[i].index;
      SphP[p].ActiveArea += AreaListGet[i].darea;
    }

  myfree(AreaListGet);
}

/*! \brief Calculates volumes of all cells that are created in refinement.
 *
 *  \param[out] vol Volumes of cells.
 *
 *  \return void
 */
void derefine_refine_compute_volumes(double *vol)
{
  int i, bit, nr;

  for(i = 0; i < DeRefMesh.Ndp; i++)
    vol[i] = 0;

  Edge_visited = mymalloc_movable(&Edge_visited, "Edge_visited", DeRefMesh.Ndt * sizeof(unsigned char));

  for(i = 0; i < DeRefMesh.Ndt; i++)
    Edge_visited[i] = 0;

  for(i = 0; i < DeRefMesh.Ndt; i++)
    {
      if(DeRefMesh.DT[i].t[0] < 0) /* deleted ? */
        continue;

      bit = 1;
      nr  = 0;

      while(Edge_visited[i] != EDGE_ALL)
        {
          if((Edge_visited[i] & bit) == 0)
            derefine_refine_process_edge(&DeRefMesh, vol, i, nr);

          bit <<= 1;
          nr++;
        }
    }

  myfree(Edge_visited);
}

#endif /* #ifndef ONEDIMS */

/*! \brief Nearest distance in x direction, accounting for periodicity.
 *
 *  \param[in] d Distance to be checked.
 *
 *  \return Nearest distance.
 */
double nearest_x(double d)
{
#if !defined(REFLECTIVE_X)
  if(d < -boxHalf_X)
    d += boxSize_X;
  if(d > boxHalf_X)
    d -= boxSize_X;
#endif /* #if !defined(REFLECTIVE_X) */
  return d;
}

/*! \brief Nearest distance in y direction, accounting for periodicity.
 *
 *  \param[in] d Distance to be checked.
 *
 *  \return Nearest distance.
 */
double nearest_y(double d)
{
#if !defined(REFLECTIVE_Y)
  if(d < -boxHalf_Y)
    d += boxSize_Y;
  if(d > boxHalf_Y)
    d -= boxSize_Y;
#endif /* #if !defined(REFLECTIVE_Y) */
  return d;
}

/* \brief Nearest distance in z direction, accounting for periodicity.
 *
 * \param[in] d Distance to be checked.
 *
 * \return Nearest distance.
 */
double nearest_z(double d)
{
#if !defined(REFLECTIVE_Z)
  if(d < -boxHalf_Z)
    d += boxSize_Z;
  if(d > boxHalf_Z)
    d -= boxSize_Z;
#endif /* #if !defined(REFLECTIVE_Z) */
  return d;
}

/*! \brief Gets "radius" of a cell.
 *
 *  Defined as the radius of a sphere with the same volume as the Voronoi cell.
 *
 *  \param[in] i Index of cell in P and SphP arrays.
 *
 *  \return radius of cell i.
 */
double get_cell_radius(int i)
{
  double cellrad;

#ifdef TWODIMS
  cellrad = sqrt(SphP[i].Volume / M_PI);
#else /* #ifdef TWODIMS */
#ifdef ONEDIMS
#ifdef ONEDIMS_SPHERICAL
  cellrad = 0.5 * (Mesh.VF[i + 1].cx - Mesh.VF[i].cx);
#else  /* #ifdef ONEDIMS_SPHERICAL */
  cellrad = 0.5 * SphP[i].Volume;
#endif /* #ifdef ONEDIMS_SPHERICAL #else */
#else  /* #ifdef ONEDIMS */
  cellrad = pow(SphP[i].Volume * 3.0 / (4.0 * M_PI), 1.0 / 3);
#endif /* #ifdef ONEDIMS #else */
#endif /* #ifdef TWODIMS */
  return cellrad;
}

/*! \brief Writes a file points_X.dat with Delaunay points.
 *
 *  Writes position as in DP structure.
 *
 *  \param[in] T tessellation for which Delaunay point positions should be
 *               written.
 *
 *  \return void
 */
void dump_points(tessellation *T)
{
  FILE *fd;
  int i;
  double xyz[3];
  char buf[1000];

  sprintf(buf, "points_%d.dat", ThisTask);
  fd = fopen(buf, "w");
  my_fwrite(&T->Ndp, sizeof(int), 1, fd);
  for(i = 0; i < T->Ndp; i++)
    {
      xyz[0] = T->DP[i].x;
      xyz[1] = T->DP[i].y;
      xyz[2] = T->DP[i].z;
      my_fwrite(xyz, sizeof(double), 3, fd);
    }
  fclose(fd);
}

/*! \brief Calculates the normals to given interfaces.
 *
 *  \param[in] T Pointer to tesslation data.
 *  \param[in] i Index of Voronoi-face in tesslation T.
 *  \param[out] geom Pointer to structure to which normal data is written.
 *
 *  \return 0 if success, -1 if interface can be ignored.
 */
int face_get_normals(tessellation *T, int i, struct geometry *geom)
{
  int li, ri;
  double surface, surface_l, surface_r;
  int present_left, present_right;
  double mm;

  face *VF  = T->VF;
  point *DP = T->DP;

  li = DP[VF[i].p1].index;
  ri = DP[VF[i].p2].index;

  if(li < 0 || ri < 0)
    return -1;

  if(li >= NumGas && DP[VF[i].p1].task == ThisTask)
    li -= NumGas;

  if(ri >= NumGas && DP[VF[i].p2].task == ThisTask)
    ri -= NumGas;

  if(DP[VF[i].p1].task == ThisTask)
    surface_l = SphP[li].SurfaceArea;
  else
    surface_l = PrimExch[li].SurfaceArea;

  if(DP[VF[i].p2].task == ThisTask)
    surface_r = SphP[ri].SurfaceArea;
  else
    surface_r = PrimExch[ri].SurfaceArea;

  if(surface_r > surface_l)
    surface = 1.0e-5 * surface_r;
  else
    surface = 1.0e-5 * surface_l;

  present_left = present_right = 0;

  /* if the area of this face is negligible compared to the surface
     of the larger cell, skip it */
  if(DP[VF[i].p1].task == ThisTask && DP[VF[i].p1].index < NumGas)
    if(TimeBinSynchronized[P[DP[VF[i].p1].index].TimeBinHydro])
      if(VF[i].area > surface)
        present_left = 1;

  if(DP[VF[i].p2].task == ThisTask && DP[VF[i].p2].index < NumGas)
    if(TimeBinSynchronized[P[DP[VF[i].p2].index].TimeBinHydro])
      if(VF[i].area > surface)
        present_right = 1;

  if(present_left == 0 && present_right == 0)
    {
#ifndef VORONOI_STATIC_MESH
      VF[i].area = 0;
#endif /* #ifndef VORONOI_STATIC_MESH */
      return -1;
    }

  /* center of face */
  geom->cx = VF[i].cx;
  geom->cy = VF[i].cy;
  geom->cz = VF[i].cz;

  /* normal vector pointing to "right" state */
  geom->nx = DP[VF[i].p2].x - DP[VF[i].p1].x;
  geom->ny = DP[VF[i].p2].y - DP[VF[i].p1].y;
  geom->nz = DP[VF[i].p2].z - DP[VF[i].p1].z;

  geom->nn = sqrt(geom->nx * geom->nx + geom->ny * geom->ny + geom->nz * geom->nz);
  geom->nx /= geom->nn;
  geom->ny /= geom->nn;
  geom->nz /= geom->nn;

  /* need an ortonormal basis */
  if(geom->nx != 0 || geom->ny != 0)
    {
      geom->mx = -geom->ny;
      geom->my = geom->nx;
      geom->mz = 0;
    }
  else
    {
      geom->mx = 1;
      geom->my = 0;
      geom->mz = 0;
    }

  mm = sqrt(geom->mx * geom->mx + geom->my * geom->my + geom->mz * geom->mz);
  geom->mx /= mm;
  geom->my /= mm;
  geom->mz /= mm;

  geom->px = geom->ny * geom->mz - geom->nz * geom->my;
  geom->py = geom->nz * geom->mx - geom->nx * geom->mz;
  geom->pz = geom->nx * geom->my - geom->ny * geom->mx;

  return 0;
}

/*! \brief Calculates distance of a cell to boundary of computational box.
 *
 *  \param[in] cell Index of cell in P and SphP structure.
 *
 *  \return Distance to border.
 */
double distance_to_border(int cell)
{
  double d1 = boxSize_X - P[cell].Pos[0];
  assert(d1 > 0);

  double d2 = P[cell].Pos[0];

  double min = fmin(d1, d2);

  d1 = boxSize_Y - P[cell].Pos[1];
  assert(d1 > 0);

  d2 = P[cell].Pos[1];

  double min2 = fmin(d1, d2);
  min         = fmin(min, min2);

  d1 = boxSize_Z - P[cell].Pos[2];
  assert(d1 > 0);

  d2   = P[cell].Pos[2];
  min2 = fmin(d1, d2);

  min = fmin(min, min2);

  return min;
}
