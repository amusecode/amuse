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
 * \file        src/mesh/voronoi/voronoi_3d.c
 * \date        05/2018
 * \brief       Routines to build a 3d Voronoi mesh.
 * \details     Note that some of these routines have the same name as the ones
 *              in voronoi_1d.c and voronoi_2d.c and just replace them in case
 *              neither the Config-option TWODIMS nor ONEDIMS is active.
 *              contains functions:
 *                void initialize_and_create_first_tetra(tessellation * T)
 *                void get_line_segments(int sphp_index, int dp_index, double
 *                  *segments, unsigned int *nof_elements, unsigned int
 *                  max_elements)
 *                void process_edge_faces_and_volumes(tessellation * T,
 *                  int tt, int nr)
 *                int derefine_refine_get_triangles(tessellation * T, int tt,
 *                  int nr, point * dtip, triangle * trilist, int ntri,
 *                  int max_n_tri)
 *                double get_tri_volume(int i, triangle * trilist)
 *                int derefine_add_point_and_split_tri(int q, triangle
 *                  * trilist, int ntri, int max_ntri, double vol)
 *                void derefine_refine_process_edge(tessellation * T,
 *                  double *vol, int tt, int nr)
 *                int insert_point(tessellation * T, int pp, int ttstart)
 *                int convex_edge_test(tessellation * T, int tt, int tip,
 *                  int *edgenr)
 *                void make_a_face_split(tessellation * T, int tt0,
 *                  int face_nr, int pp, int tt1, int tt2, int qq1, int qq2)
 *                void make_an_edge_split(tessellation * T, int tt0,
 *                  int edge_nr, int count, int pp, int *ttlist)
 *                void make_a_4_to_4_flip(tessellation * T, int tt,
 *                  int tip_index, int edge_nr)
 *                void make_a_1_to_4_flip(tessellation * T, int pp, int tt0,
 *                  int tt1, int tt2, int tt3)
 *                void make_a_3_to_2_flip(tessellation * T, int tt0, int tt1,
 *                  int tt2, int tip, int edge, int bottom)
 *                void make_a_2_to_3_flip(tessellation * T, int tt0, int tip,
 *                  int tt1, int bottom, int qq, int tt2)
 *                int get_tetra(tessellation * T, point * p, int *moves,
 *                  int ttstart, int *flag, int *edgeface_nr)
 *                int InTetra(tessellation * T, int tt, point * p,
 *                  int *edgeface_nr, int *nexttetra)
 *                void compute_circumcircles(tessellation * T)
 *                void calc_mpz_determinant(mpz_t det, mpz_t ax, mpz_t ay,
 *                  mpz_t az, mpz_t bx, mpz_t by, mpz_t bz, mpz_t cx,
 *                  mpz_t cy, mpz_t cz)
 *                void get_circumcircle_exact(tessellation * T, int tt,
 *                  double *x, double *y, double *z)
 *                void update_circumcircle(tessellation * T, int tt)
 *                int test_tetra_orientation(point * p0, point * p1,
 *                  point * p2, point * p3)
 *                double calculate_tetra_volume(point * p0, point * p1,
 *                  point * p2, point * p3)
 *                void add_row(double *m, int r1, int r2, double fac)
 *                int solve_linear_equations(double *m, double *res)
 *                void set_integers_for_pointer(point * p)
 *                int InSphere_Exact(point * p0, point * p1, point * p2,
 *                  point * p3, point * p)
 *                int InSphere_Quick(point * p0, point * p1, point * p2,
 *                  point * p3, point * p)
 *                int InSphere_Errorbound(point * p0, point * p1, point * p2,
 *                  point * p3, point * p)
 *                int Orient3d_Exact(point * p0, point * p1, point * p2,
 *                  point * p3)
 *                int Orient3d_Quick(point * p0, point * p1, point * p2,
 *                  point * p3)
 *                int Orient3d(point * p0, point * p1, point * p2, point * p3)
 *                int compare_face_sort(const void *a, const void *b)
 *                void get_voronoi_face_vertex_indices(tessellation * T)
 *                void get_voronoi_face_vertex_coordinates(tessellation * T)
 *                void sort_faces_by_ID(void)
 *                void write_voronoi_face_vertex_indices(tessellation * T,
 *                  char *fname1, char *fname2, int writeTask, int lastTask)
 *                void write_voronoi_face_vertex_coordinates(tessellation * T,
 *                  char *fname, int writeTask, int lastTask)
 *                void write_voronoi_mesh(tessellation * T, char *fname,
 *                  int writeTask, int lastTask)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 21.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gmp.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../main/allvars.h"
#include "../../main/proto.h"

#include "voronoi.h"

#if !defined(TWODIMS) && !defined(ONEDIMS) /* will only be compiled in 3D case */

#define INSIDE_EPS 1.0e-6
#define GAUSS_EPS 1.0e-8

const int access_triangles[4][3] = {{1, 3, 2}, {0, 2, 3}, {0, 3, 1}, {0, 1, 2}};

const int edge_start[6]     = {0, 0, 0, 1, 1, 2};
const int edge_end[6]       = {1, 2, 3, 2, 3, 3};
const int edge_opposite[6]  = {3, 1, 2, 3, 0, 1};
const int edge_nexttetra[6] = {2, 3, 1, 0, 2, 0};

/*! \brief Initializes 3d tessellation and create all-enclosing tetrahedron.
 *
 *  \param[out] T Pointer to tessellation structure which is set and its arrays
 *  are allocated in this routine.
 *
 *  \return void
 */
void initialize_and_create_first_tetra(tessellation *T)
{
  point *p;
  int i, n;

  T->MaxNdp = T->Indi.AllocFacNdp;
  T->MaxNdt = T->Indi.AllocFacNdt;
  T->MaxNvf = T->Indi.AllocFacNvf;

  T->Ndp = 0;
  T->Ndt = 0;
  T->Nvf = 0;

  T->VF = mymalloc_movable(&T->VF, "VF", T->MaxNvf * sizeof(face));

  T->DP = mymalloc_movable(&T->DP, "DP", (T->MaxNdp + 5) * sizeof(point));
  T->DP += 5;

  T->DT = mymalloc_movable(&T->DT, "DT", T->MaxNdt * sizeof(tetra));

  /* construct all encompassing huge tetrahedron */

  double box, tetra_incircle, tetra_sidelength, tetra_height, tetra_face_height;

  box = boxSize_X;
  if(box < boxSize_Y)
    box = boxSize_Y;
  if(box < boxSize_Z)
    box = boxSize_Z;

  tetra_incircle    = 1.5 * box;
  tetra_sidelength  = tetra_incircle * sqrt(24);
  tetra_height      = sqrt(2.0 / 3) * tetra_sidelength;
  tetra_face_height = sqrt(3.0) / 2.0 * tetra_sidelength;

  point *DP = T->DP;
  tetra *DT = T->DT;

  /* first, let's make the points */
  DP[-4].x = 0.5 * tetra_sidelength;
  DP[-4].y = -1.0 / 3 * tetra_face_height;
  DP[-4].z = -0.25 * tetra_height;

  DP[-3].x = 0;
  DP[-3].y = 2.0 / 3 * tetra_face_height;
  DP[-3].z = -0.25 * tetra_height;

  DP[-2].x = -0.5 * tetra_sidelength;
  DP[-2].y = -1.0 / 3 * tetra_face_height;
  DP[-2].z = -0.25 * tetra_height;

  DP[-1].x = 0;
  DP[-1].y = 0;
  DP[-1].z = 0.75 * tetra_height;

  for(i = -4; i <= -1; i++)
    {
      DP[i].x += 0.5 * box;
      DP[i].y += 0.5 * box;
      DP[i].z += 0.5 * box;
    }

  for(i = -4, p = &DP[-4]; i < 0; i++, p++)
    {
      p->index   = -1;
      p->task    = ThisTask;
      p->timebin = 0;
    }

  /* we also define a neutral element at infinity */
  DPinfinity = -5;

  DP[DPinfinity].x       = MAX_DOUBLE_NUMBER;
  DP[DPinfinity].y       = MAX_DOUBLE_NUMBER;
  DP[DPinfinity].z       = MAX_DOUBLE_NUMBER;
  DP[DPinfinity].index   = -1;
  DP[DPinfinity].task    = ThisTask;
  DP[DPinfinity].timebin = 0;

  /* now let's make the big tetrahedron */
  DT[0].p[0] = -4;
  DT[0].p[1] = -3;
  DT[0].p[2] = -2;
  DT[0].p[3] = -1;

  /* On the outer faces, we attach tetrahedra with the neutral element as tip.
   * This way we will be able to navigate nicely within the tesselation,
   * and all tetrahedra have defined neighbouring tetrahedra.
   */

  for(i = 0; i < 4; i++)
    {
      n = i + 1; /* tetra index */

      DT[0].t[i] = n;
      DT[0].s[i] = 3;

      DT[n].t[3] = 0;
      DT[n].s[3] = i;
      DT[n].p[3] = DPinfinity;
    }

  DT[1].p[0] = DT[0].p[1];
  DT[1].p[1] = DT[0].p[2];
  DT[1].p[2] = DT[0].p[3];

  DT[2].p[0] = DT[0].p[0];
  DT[2].p[1] = DT[0].p[3];
  DT[2].p[2] = DT[0].p[2];

  DT[3].p[0] = DT[0].p[0];
  DT[3].p[1] = DT[0].p[1];
  DT[3].p[2] = DT[0].p[3];

  DT[4].p[0] = DT[0].p[0];
  DT[4].p[1] = DT[0].p[2];
  DT[4].p[2] = DT[0].p[1];

  DT[1].t[0] = 2;
  DT[2].t[0] = 1;
  DT[1].s[0] = 0;
  DT[2].s[0] = 0;

  DT[1].t[1] = 3;
  DT[3].t[0] = 1;
  DT[1].s[1] = 0;
  DT[3].s[0] = 1;

  DT[1].t[2] = 4;
  DT[4].t[0] = 1;
  DT[1].s[2] = 0;
  DT[4].s[0] = 2;

  DT[2].t[2] = 3;
  DT[3].t[1] = 2;
  DT[2].s[2] = 1;
  DT[3].s[1] = 2;

  DT[2].t[1] = 4;
  DT[4].t[2] = 2;
  DT[2].s[1] = 2;
  DT[4].s[2] = 1;

  DT[3].t[2] = 4;
  DT[4].t[1] = 3;
  DT[3].s[2] = 1;
  DT[4].s[1] = 2;

  T->Ndt = 5; /* we'll start out with 5 tetras */

  CentralOffsetX = 0.5 * box - 0.5000001 * tetra_sidelength;
  CentralOffsetY = 0.5 * box - (1.0000001 / 3) * tetra_face_height;
  CentralOffsetZ = 0.5 * box - 0.25000001 * tetra_height;

  ConversionFac = 1.0 / (1.001 * tetra_sidelength);

#ifndef OPTIMIZE_MEMORY_USAGE
  for(i = -4; i < 0; i++)
    set_integers_for_point(T, i);
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE */
}

#ifdef TETRA_INDEX_IN_FACE
/*! \brief Gets the line segments of a Voronoi cell.
 *
 *  Warning: The correspondance sphp_index == dp_index holds only for a global
 *  timestep!
 *
 *  \param[in] sphp_index The index of the Voronoi cell.
 *  \param[in] dp_index The index of the corresponding Delaunay point.
 *  \param[out] segments The array in which the line segments are stored.
 *  \param[out] nof_elements The number of elements written in segments during
 *              this function call.
 *  \param[in] max_elements The maximum size of the segments array.
 *
 *  \return void
 */
void get_line_segments(int sphp_index, int dp_index, double *segments, unsigned int *nof_elements, unsigned int max_elements)
{
  // index for segments array
  unsigned int a = 0;

  int edge      = SphP[sphp_index].first_connection;
  int last_edge = SphP[sphp_index].last_connection;

  // loop over all interfaces of the cell
  while(1)
    {
      int dq_index = DC[edge].dp_index;

      // one of the tetrahedras around the Delaunay connection
      int tt   = DC[edge].dt_index;
      tetra *t = &Mesh.DT[tt];

      // find the local index of the edge
      int nr = 6;
      int e, dp_start_index, dp_end_index;

      for(e = 0; e < 6; e++)
        {
          dp_start_index = t->p[edge_start[e]];
          dp_end_index   = t->p[edge_end[e]];

          if((dp_start_index == dp_index && dp_end_index == dq_index) || (dp_start_index == dq_index && dp_end_index == dp_index))
            {
              nr = e;
              break;
            }
        }

      // ensure that the local edge index has been found
      assert(nr != 6);

      // already set: t,tt,nr
      int i, j, k, l, m, ii, jj, kk, ll, nn;
      tetra *prev, *next;
      tetra_center *prevc, *nextc;

      i = edge_start[nr];
      j = edge_end[nr];
      k = edge_opposite[nr];
      l = edge_nexttetra[nr];

      prev  = t;
      prevc = &Mesh.DTC[tt];

      do
        {
          nn    = prev->t[l];
          next  = &Mesh.DT[nn];
          nextc = &Mesh.DTC[nn];

          if(a > max_elements - 7)
            {
              terminate("termination in voronoi_3d.c get_line_segments: not enough memory!");
            }

          segments[a++] = prevc->cx;
          segments[a++] = prevc->cy;
          segments[a++] = prevc->cz;
          segments[a++] = nextc->cx;
          segments[a++] = nextc->cy;
          segments[a++] = nextc->cz;

          for(m = 0, ll = ii = jj = -1; m < 4; m++)
            {
              if(next->p[m] == prev->p[k])
                ll = m;
              if(next->p[m] == prev->p[i])
                ii = m;
              if(next->p[m] == prev->p[j])
                jj = m;
            }

          if(ll < 0 || ii < 0 || jj < 0)
            terminate("inconsistency");

          kk = 6 - (ll + ii + jj);

          prev  = next;
          prevc = nextc;

          i = ii;
          l = ll;
          j = jj;
          k = kk;
        }
      while(next != t);

      if(edge == last_edge)
        {
          break;
        }

      edge = DC[edge].next;

    }  // end of while loop

  *nof_elements = a;

  return;
}
#endif /* #ifdef TETRA_INDEX_IN_FACE */

/*! \brief Calculate cell volumes and face areas of mesh.
 *
 *  \param[in, out] T Pointer to tessellation.
 *  \param[in] tt Index in DT array.
 *  \param[in] nr Index in edges.
 *
 *  \return void
 */
void process_edge_faces_and_volumes(tessellation *T, int tt, int nr)
{
  int i, j, k, l, m, ii, jj, kk, ll, nn, count, nr_next, p1, p2;
  face *f;
  tetra *prev, *next;
  tetra_center *prevc, *nextc;
  double ax, ay, az;
  double bx, by, bz;
  double cx, cy, cz;
  double nx, ny, nz;
  double sx, sy, sz;
  double hhx, hhy, hhz;
  double darea, dvol, h;

  if(T->Nvf + 1 >= T->MaxNvf)
    {
      T->Indi.AllocFacNvf *= ALLOC_INCREASE_FACTOR;
      T->MaxNvf = T->Indi.AllocFacNvf;
#ifdef VERBOSE
      printf("Task=%d: increase memory allocation, MaxNvf=%d Indi.AllocFacNvf=%g\n", ThisTask, T->MaxNvf, T->Indi.AllocFacNvf);
#endif /* #ifdef VERBOSE */
      T->VF = myrealloc_movable(T->VF, T->MaxNvf * sizeof(face));

      if(T->Nvf + 1 >= T->MaxNvf)
        terminate("Nvf larger than MaxNvf");
    }

  tetra *DT         = T->DT;
  point *DP         = T->DP;
  face *VF          = T->VF;
  tetra_center *DTC = T->DTC;

  tetra *t = &DT[tt];

  i = edge_start[nr];
  j = edge_end[nr];
  k = edge_opposite[nr];
  l = edge_nexttetra[nr];

  Edge_visited[tt] |= (1 << nr);

  p1 = t->p[i];
  p2 = t->p[j];

  f = &VF[T->Nvf++];

  f->area = 0;
  f->p1   = p1;
  f->p2   = p2;

  f->cx = 0;
  f->cy = 0;
  f->cz = 0;

#ifdef TETRA_INDEX_IN_FACE
  f->dt_index = tt;
#endif /* #ifdef TETRA_INDEX_IN_FACE */

  hhx = 0.5 * (DP[p1].x - DP[p2].x);
  hhy = 0.5 * (DP[p1].y - DP[p2].y);
  hhz = 0.5 * (DP[p1].z - DP[p2].z);

  h = sqrt(hhx * hhx + hhy * hhy + hhz * hhz);

  cx = DTC[tt].cx;
  cy = DTC[tt].cy;
  cz = DTC[tt].cz;

  count = 0;

  prev  = t;
  prevc = &DTC[tt];
  do
    {
      nn    = prev->t[l];
      next  = &DT[nn];
      nextc = &DTC[nn];

      if(prev != t && next != t)
        {
          ax = prevc->cx - cx;
          ay = prevc->cy - cy;
          az = prevc->cz - cz;

          bx = nextc->cx - cx;
          by = nextc->cy - cy;
          bz = nextc->cz - cz;

          nx = ay * bz - az * by;
          ny = az * bx - ax * bz;
          nz = ax * by - ay * bx;

          sx = nextc->cx + prevc->cx + cx;
          sy = nextc->cy + prevc->cy + cy;
          sz = nextc->cz + prevc->cz + cz;

          darea = 0.5 * sqrt(nx * nx + ny * ny + nz * nz);
          f->area += darea;

          darea *= (1.0 / 3);

          f->cx += darea * sx;
          f->cy += darea * sy;
          f->cz += darea * sz;
        }

      for(m = 0, ll = ii = jj = -1; m < 4; m++)
        {
          if(next->p[m] == prev->p[k])
            ll = m;
          if(next->p[m] == prev->p[i])
            ii = m;
          if(next->p[m] == prev->p[j])
            jj = m;
        }

      if(ll < 0 || ii < 0 || jj < 0)
        terminate("inconsistency");

      kk = 6 - (ll + ii + jj);

      /* need to determine the edge number to be able to flag it */

      for(nr_next = 0; nr_next < 6; nr_next++)
        if((edge_start[nr_next] == ii && edge_end[nr_next] == jj) || (edge_start[nr_next] == jj && edge_end[nr_next] == ii))
          {
            if((Edge_visited[nn] & (1 << nr_next)) && next != t)
              terminate("inconsistency");

            Edge_visited[nn] |= (1 << nr_next);
            break;
          }

      prev  = next;
      prevc = nextc;
      i     = ii;
      l     = ll;
      j     = jj;
      k     = kk;

      count++;

      if(count > 1000)
        terminate("count is too large");
    }
  while(next != t);

  i = edge_start[nr];
  j = edge_end[nr];

  if(f->area)
    {
      f->cx /= f->area;
      f->cy /= f->area;
      f->cz /= f->area;
    }

#ifdef REFINEMENT_MERGE_CELLS
  f->t  = tt;
  f->nr = nr; /* delaunay tetra and edge number that generated this face */
#endif        /* #ifdef REFINEMENT_MERGE_CELLS */

  dvol = (1.0 / 3) * f->area * h;

#if defined(REGULARIZE_MESH_FACE_ANGLE) || defined(OUTPUT_MESH_FACE_ANGLE)
  double angle = sqrt(f->area / M_PI) / h;
#endif /* #if defined(REGULARIZE_MESH_FACE_ANGLE) || defined(OUTPUT_MESH_FACE_ANGLE) */

  if(DP[p1].task == ThisTask && DP[p1].index >= 0 && DP[p1].index < NumGas)
    {
      if(TimeBinSynchronized[P[DP[p1].index].TimeBinHydro])
        {
          SphP[DP[p1].index].Volume += dvol;
          SphP[DP[p1].index].SurfaceArea += f->area;

#if defined(REGULARIZE_MESH_FACE_ANGLE) || defined(OUTPUT_MESH_FACE_ANGLE)
          if(SphP[DP[p1].index].MaxFaceAngle < angle)
            SphP[DP[p1].index].MaxFaceAngle = angle;
#endif /* #if defined(REGULARIZE_MESH_FACE_ANGLE) || defined(OUTPUT_MESH_FACE_ANGLE) */

#ifdef OUTPUT_SURFACE_AREA
          if(f->area)
            SphP[DP[p1].index].CountFaces++;
#endif /* #ifdef OUTPUT_SURFACE_AREA */

#if defined(REFINEMENT_SPLIT_CELLS)
          if(SphP[DP[p1].index].MinimumEdgeDistance > h)
            SphP[DP[p1].index].MinimumEdgeDistance = h;
#endif /* #if defined(REFINEMENT_SPLIT_CELLS) */
          /* let's now compute the center-of-mass of the pyramid at the bottom top */
          sx = 0.75 * f->cx + 0.25 * DP[p1].x;
          sy = 0.75 * f->cy + 0.25 * DP[p1].y;
          sz = 0.75 * f->cz + 0.25 * DP[p1].z;

          SphP[DP[p1].index].Center[0] += dvol * sx;
          SphP[DP[p1].index].Center[1] += dvol * sy;
          SphP[DP[p1].index].Center[2] += dvol * sz;
        }
    }

  if(DP[p2].task == ThisTask && DP[p2].index >= 0 && DP[p2].index < NumGas)
    {
      if(TimeBinSynchronized[P[DP[p2].index].TimeBinHydro])
        {
          SphP[DP[p2].index].Volume += dvol;
          SphP[DP[p2].index].SurfaceArea += f->area;

#if defined(REGULARIZE_MESH_FACE_ANGLE) || defined(OUTPUT_MESH_FACE_ANGLE)
          if(SphP[DP[p2].index].MaxFaceAngle < angle)
            SphP[DP[p2].index].MaxFaceAngle = angle;
#endif /* #if defined(REGULARIZE_MESH_FACE_ANGLE) || defined(OUTPUT_MESH_FACE_ANGLE) */

#ifdef OUTPUT_SURFACE_AREA
          if(f->area)
            SphP[DP[p2].index].CountFaces++;
#endif /* #ifdef OUTPUT_SURFACE_AREA */
#if defined(REFINEMENT_SPLIT_CELLS)
          if(SphP[DP[p2].index].MinimumEdgeDistance > h)
            SphP[DP[p2].index].MinimumEdgeDistance = h;
#endif /* #if defined(REFINEMENT_SPLIT_CELLS) */
          /* let's now compute the center-of-mass of the pyramid on top */
          sx = 0.75 * f->cx + 0.25 * DP[p2].x;
          sy = 0.75 * f->cy + 0.25 * DP[p2].y;
          sz = 0.75 * f->cz + 0.25 * DP[p2].z;

          SphP[DP[p2].index].Center[0] += dvol * sx;
          SphP[DP[p2].index].Center[1] += dvol * sy;
          SphP[DP[p2].index].Center[2] += dvol * sz;
        }
    }

  int low_p, high_p;

  if(DP[p1].ID < DP[p2].ID)
    {
      low_p  = p1;
      high_p = p2;
    }
  else
    {
      low_p  = p2;
      high_p = p1;
    }

  int this_task_responsible_flag = 0;

  if(TimeBinSynchronized[DP[low_p].timebin]) /* the one with the lower ID is active */
    {
      /* we need to check whether the one with the lower ID is a local particle */
      if(DP[low_p].task == ThisTask && DP[low_p].index >= 0 && DP[low_p].index < NumGas)
        this_task_responsible_flag = 1;
    }
  else if(TimeBinSynchronized[DP[high_p].timebin]) /* only the side with the higher ID is active */
    {
      /* we need to check whether we hold the one with the higher ID, if yes, we'll do it */
      if(DP[high_p].task == ThisTask && DP[high_p].index >= 0 && DP[high_p].index < NumGas)
        this_task_responsible_flag = 1;
    }

  if(this_task_responsible_flag)
    {
      for(k = 0; k < 2; k++)
        {
          int p, q;

          if(k == 0)
            {
              q = p1;
              p = DP[q].index;
            }
          else
            {
              q = p2;
              p = DP[q].index;
            }

          if(DP[q].task == ThisTask)
            {
              if(DP[q].index >= NumGas) /* this is a local ghost point */
                p -= NumGas;

              SphP[p].ActiveArea += f->area;
            }
          else
            {
              /* here we have a foreign ghost point */
              if(DP[q].originalindex < 0)
                terminate("should not happen");

              if(Narea >= MaxNarea)
                {
                  T->Indi.AllocFacNflux *= ALLOC_INCREASE_FACTOR;
                  MaxNarea = T->Indi.AllocFacNflux;
                  AreaList = myrealloc_movable(AreaList, MaxNarea * sizeof(struct area_list_data));

                  if(Narea >= MaxNarea)
                    terminate("Narea >= MaxNarea");
                }

              AreaList[Narea].task  = DP[q].task;
              AreaList[Narea].index = DP[q].originalindex;
              AreaList[Narea].darea = f->area;
              Narea++;
            }
        }
    }
}

/*! \brief Gathers tetrahedron data as elements in array called 'trilist'.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] tt Index of tetrahedron in T->DT array.
 *  \param[in] nr Index in (global) edge arrays.
 *  \param[in] dtip Point representing tip of tetrahedron.
 *  \param[out] trilist List of triangles.
 *  \param[in] ntri Index in trilist which should be filled.
 *  \param[in] max_n_tri Maximum index in trilist.
 *
 *  \return New length of trilist data.
 */
int derefine_refine_get_triangles(tessellation *T, int tt, int nr, point *dtip, triangle *trilist, int ntri, int max_n_tri)
{
  tetra *DT         = T->DT;
  tetra_center *DTC = T->DTC;

  int i, j, k, l, m, ii, jj, kk, ll, nn, count;
  tetra *prev, *next;
  tetra_center *prevc, *nextc;
  double cx, cy, cz;

  tetra *t = &DT[tt];

  i = edge_start[nr];
  j = edge_end[nr];
  k = edge_opposite[nr];
  l = edge_nexttetra[nr];

  cx = DTC[tt].cx;
  cy = DTC[tt].cy;
  cz = DTC[tt].cz;

  count = 0;

  prev  = t;
  prevc = &DTC[tt];
  do
    {
      nn    = prev->t[l];
      next  = &DT[nn];
      nextc = &DTC[nn];

      if(prev != t && next != t)
        {
          if(ntri >= max_n_tri)
            terminate("ntri >= max_n_tri");

          trilist[ntri].p[0][0] = cx;
          trilist[ntri].p[0][1] = cy;
          trilist[ntri].p[0][2] = cz;

          trilist[ntri].p[1][0] = prevc->cx;
          trilist[ntri].p[1][1] = prevc->cy;
          trilist[ntri].p[1][2] = prevc->cz;

          trilist[ntri].p[2][0] = nextc->cx;
          trilist[ntri].p[2][1] = nextc->cy;
          trilist[ntri].p[2][2] = nextc->cz;

          trilist[ntri].p[3][0] = dtip->x;
          trilist[ntri].p[3][1] = dtip->y;
          trilist[ntri].p[3][2] = dtip->z;

          if(get_tri_volume(ntri, trilist) < 0)
            {
              /* swap two points to get proper orientation */
              trilist[ntri].p[3][0] = nextc->cx;
              trilist[ntri].p[3][1] = nextc->cy;
              trilist[ntri].p[3][2] = nextc->cz;

              trilist[ntri].p[2][0] = dtip->x;
              trilist[ntri].p[2][1] = dtip->y;
              trilist[ntri].p[2][2] = dtip->z;
            }

          ntri++;
        }

      for(m = 0, ll = ii = jj = -1; m < 4; m++)
        {
          if(next->p[m] == prev->p[k])
            ll = m;
          if(next->p[m] == prev->p[i])
            ii = m;
          if(next->p[m] == prev->p[j])
            jj = m;
        }

      if(ll < 0 || ii < 0 || jj < 0)
        terminate("inconsistency");

      kk = 6 - (ll + ii + jj);

      prev  = next;
      prevc = nextc;
      i     = ii;
      l     = ll;
      j     = jj;
      k     = kk;

      count++;

      if(count > 1000)
        terminate("count is too large");
    }
  while(next != t);

  return ntri;
}

/*! \brief Returns volume of a tetrahedron.
 *
 *  \param[in] i Index of tetrahedron in trilist.
 *  \param[in] trilist Array with tetrahedra.
 *
 *  \return Volume of tetrahedron.
 */
double get_tri_volume(int i, triangle *trilist)
{
  double nx, ny, nz;

  double *p0 = &trilist[i].p[0][0];
  double *p1 = &trilist[i].p[1][0];
  double *p2 = &trilist[i].p[2][0];
  double *p3 = &trilist[i].p[3][0];

  nx = (p1[1] - p0[1]) * (p2[2] - p0[2]) - (p1[2] - p0[2]) * (p2[1] - p0[1]);
  ny = (p1[2] - p0[2]) * (p2[0] - p0[0]) - (p1[0] - p0[0]) * (p2[2] - p0[2]);
  nz = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0]);

  return (nx * (p3[0] - p0[0]) + ny * (p3[1] - p0[1]) + nz * (p3[2] - p0[2])) / 6.0;
}

/*! \brief Add point and adjust tetrahedra accordingly.
 *
 *  \param[in] q Index of point in DP array.
 *  \param[in, out] trilist Array of tetrahedra.
 *  \param[in] ntri Number of elements in trilist before splitting.
 *  \param[in] max_ntri Maximum number of tetrahedron allowed.
 *  \param[in] vol Volume of tetrahedron to be split.
 *
 *  \return Updated number of triangles.
 */
int derefine_add_point_and_split_tri(int q, triangle *trilist, int ntri, int max_ntri, double vol)
{
#define MIN_VOL_FAC 1.0e-6
  double m[3], n[3], sc[4], *a;
  double cut[3][3], p[8][3], ed[3];
  int i, j, k, l, nnew, flag[4], count, oldq;
  double vvi, vlargest, vv[5];
  int ilargest, nadd;

  for(i = 0, nnew = ntri; i < ntri; i++)
    {
      if(q < 0 || q >= Mesh.Ndp)
        {
          char buf[1000];
          sprintf(buf, "q=%d\n", q);
          terminate(buf);
        }

      if(trilist[i].owner < 0 || trilist[i].owner >= Mesh.Ndp)
        {
          char buf[1000];
          sprintf(buf, "trilist[i].owner=%d\n", trilist[i].owner);
          terminate(buf);
        }

      /* midpoint */
      m[0] = 0.5 * (Mesh.DP[q].x + Mesh.DP[trilist[i].owner].x);
      m[1] = 0.5 * (Mesh.DP[q].y + Mesh.DP[trilist[i].owner].y);
      m[2] = 0.5 * (Mesh.DP[q].z + Mesh.DP[trilist[i].owner].z);

      n[0] = (Mesh.DP[q].x - Mesh.DP[trilist[i].owner].x);
      n[1] = (Mesh.DP[q].y - Mesh.DP[trilist[i].owner].y);
      n[2] = (Mesh.DP[q].z - Mesh.DP[trilist[i].owner].z);

      if(q == trilist[i].owner)
        terminate("q == trilist[i].owner");

      for(k = 0, count = 0; k < 4; k++) /* determine the side of each point */
        {
          a = &trilist[i].p[k][0];

          sc[k] = (a[0] - m[0]) * n[0] + (a[1] - m[1]) * n[1] + (a[2] - m[2]) * n[2];

          if(sc[k] > 0)
            {
              flag[k] = 1;
              count++;
            }
          else
            flag[k] = 0;
        }

      switch(count)
        {
          case 0: /* the whole tetra is on the side of current owner - nothing to be done */
            break;

          case 4:                 /* the whole tetra is on the side of new point */
            trilist[i].owner = q; /* change owner */
            break;

          case 1:
          case 3:

            /* we have one point on either side */
            /* for count=1 the tip of the tetra is cut off and assigned to the new point. */
            /* the rest is subdivided into three tetras */

            if(nnew + 3 > max_ntri)
              {
                terminate("nnew + 3 > max_ntri");
              }

            trilist[nnew]     = trilist[i];
            trilist[nnew + 1] = trilist[i];
            trilist[nnew + 2] = trilist[i];

            /* find the point index that is on the other side */
            for(k = 0; k < 4; k++)
              {
                if(flag[k] == 1 && count == 1)
                  break;
                if(flag[k] == 0 && count == 3)
                  break;
              }

            /* determine the cut-points on the corresponding edges */

            for(j = 0; j < 3; j++)
              {
                double *b = trilist[i].p[k];
                double *a = trilist[i].p[access_triangles[k][j]];

                for(l = 0; l < 3; l++)
                  ed[l] = a[l] - b[l];

                double prod = (ed[0] * n[0] + ed[1] * n[1] + ed[2] * n[2]);
                double t;

                if(prod)
                  t = -sc[k] / prod;
                else
                  t = 0.5;

                if(t < 0)
                  t = 0;
                if(t > 1)
                  t = 1;

                for(l = 0; l < 3; l++)
                  cut[j][l] = b[l] + t * ed[l];
              }

            /* modify the tetra that's assigned to the new point */
            for(j = 0; j < 3; j++)
              {
                double *a = trilist[i].p[access_triangles[k][j]];
                for(l = 0; l < 3; l++)
                  a[l] = cut[j][l];
              }

            oldq = trilist[i].owner;

            if(count == 1)
              trilist[i].owner = q;

            /* modify the three new tetras */

            for(l = 0; l < 3; l++)
              {
                trilist[nnew].p[k][l] = cut[0][l];

                trilist[nnew + 1].p[access_triangles[k][0]][l] = cut[0][l];
                trilist[nnew + 1].p[k][l]                      = cut[2][l];

                trilist[nnew + 2].p[access_triangles[k][0]][l] = cut[0][l];
                trilist[nnew + 2].p[access_triangles[k][2]][l] = cut[2][l];
                trilist[nnew + 2].p[k][l]                      = cut[1][l];
              }

            if(count == 1)
              {
                trilist[nnew].owner     = oldq;
                trilist[nnew + 1].owner = oldq;
                trilist[nnew + 2].owner = oldq;
              }
            else
              {
                trilist[nnew].owner     = q;
                trilist[nnew + 1].owner = q;
                trilist[nnew + 2].owner = q;
              }

            nadd = 3;

            vvi = fabs(get_tri_volume(i, trilist));
            for(l = 0; l < nadd; l++)
              vv[l] = fabs(get_tri_volume(nnew + l, trilist));

            /* determine largest */
            ilargest = i;
            vlargest = vvi;
            for(l = 0; l < nadd; l++)
              if(vv[l] > vlargest)
                {
                  vlargest = vv[l];
                  ilargest = nnew + l;
                }
            if(i != ilargest)
              {
                /* swap the largest to location i */
                triangle trisave  = trilist[i];
                trilist[i]        = trilist[ilargest];
                trilist[ilargest] = trisave;

                vv[ilargest - nnew] = vvi;
              }

            for(l = 0; l < nadd; l++)
              {
                if(vv[l] < MIN_VOL_FAC * vol)
                  {
                    vv[l]             = vv[nadd - 1];
                    trilist[nnew + l] = trilist[nnew + nadd - 1];
                    l--;
                    nadd--;
                  }
              }

            nnew += nadd;
            break;

          case 2:
            /* we have two points on either side */

            if(nnew + 5 > max_ntri)
              terminate("nnew + 5 > max_ntri");

            int kfirst, ksecond, jfirst, jsecond;

            if(flag[2] == 1 && flag[3] == 1)
              {
                kfirst  = 3;
                ksecond = 2;
                jfirst  = 0;
                jsecond = 1;
              }
            else if(flag[1] == 1 && flag[3] == 1)
              {
                kfirst  = 3;
                ksecond = 1;
                jfirst  = 2;
                jsecond = 0;
              }
            else if(flag[0] == 1 && flag[3] == 1)
              {
                kfirst  = 3;
                ksecond = 0;
                jfirst  = 1;
                jsecond = 2;
              }
            else if(flag[1] == 1 && flag[2] == 1)
              {
                kfirst  = 1;
                ksecond = 2;
                jfirst  = 3;
                jsecond = 0;
              }
            else if(flag[0] == 1 && flag[2] == 1)
              {
                kfirst  = 0;
                ksecond = 2;
                jfirst  = 1;
                jsecond = 3;
              }
            else if(flag[0] == 1 && flag[1] == 1)
              {
                kfirst  = 0;
                ksecond = 1;
                jfirst  = 3;
                jsecond = 2;
              }
            else
              terminate("can't be");

            int next = 0;

            for(l = 0; l < 3; l++)
              p[next][l] = trilist[i].p[kfirst][l];
            next++;

            /* determine cuts with the corresponding two edges */
            {
              double *b = trilist[i].p[kfirst];
              double *a = trilist[i].p[jfirst];

              for(l = 0; l < 3; l++)
                ed[l] = a[l] - b[l];

              double prod = (ed[0] * n[0] + ed[1] * n[1] + ed[2] * n[2]);
              double t;

              if(prod)
                t = -sc[kfirst] / prod;
              else
                t = 0.5;

              if(t < 0)
                t = 0;
              if(t > 1)
                t = 1;

              for(l = 0; l < 3; l++)
                p[next][l] = b[l] + t * ed[l];
              next++;

              for(l = 0; l < 3; l++)
                p[next][l] = a[l];
              next++;
            }

            {
              double *b = trilist[i].p[kfirst];
              double *a = trilist[i].p[jsecond];

              for(l = 0; l < 3; l++)
                ed[l] = a[l] - b[l];

              double prod = (ed[0] * n[0] + ed[1] * n[1] + ed[2] * n[2]);
              double t;

              if(prod)
                t = -sc[kfirst] / prod;
              else
                t = 0.5;

              if(t < 0)
                t = 0;
              if(t > 1)
                t = 1;

              for(l = 0; l < 3; l++)
                p[next][l] = b[l] + t * ed[l];
              next++;

              for(l = 0; l < 3; l++)
                p[next][l] = a[l];
              next++;
            }

            for(l = 0; l < 3; l++)
              p[next][l] = trilist[i].p[ksecond][l];
            next++;

            {
              double *b = trilist[i].p[ksecond];
              double *a = trilist[i].p[jfirst];

              for(l = 0; l < 3; l++)
                ed[l] = a[l] - b[l];

              double prod = (ed[0] * n[0] + ed[1] * n[1] + ed[2] * n[2]);
              double t;

              if(prod)
                t = -sc[ksecond] / prod;
              else
                t = 0.5;

              if(t < 0)
                t = 0;
              if(t > 1)
                t = 1;

              for(l = 0; l < 3; l++)
                p[next][l] = b[l] + t * ed[l];
              next++;
            }

            {
              double *b = trilist[i].p[ksecond];
              double *a = trilist[i].p[jsecond];

              for(l = 0; l < 3; l++)
                ed[l] = a[l] - b[l];

              double prod = (ed[0] * n[0] + ed[1] * n[1] + ed[2] * n[2]);
              double t;

              if(prod)
                t = -sc[ksecond] / prod;
              else
                t = 0.5;

              if(t < 0)
                t = 0;
              if(t > 1)
                t = 1;

              for(l = 0; l < 3; l++)
                p[next][l] = b[l] + t * ed[l];
              next++;
            }

            oldq = trilist[i].owner;

            /* now let's initialize the new triangles */
            for(l = 0; l < 3; l++)
              {
                /* first the ones that get to the new side */
                trilist[i].p[0][l] = p[0][l];
                trilist[i].p[1][l] = p[6][l];
                trilist[i].p[2][l] = p[5][l];
                trilist[i].p[3][l] = p[7][l];

                trilist[nnew].p[0][l] = p[1][l];
                trilist[nnew].p[1][l] = p[3][l];
                trilist[nnew].p[2][l] = p[7][l];
                trilist[nnew].p[3][l] = p[0][l];

                trilist[nnew + 1].p[0][l] = p[1][l];
                trilist[nnew + 1].p[1][l] = p[7][l];
                trilist[nnew + 1].p[2][l] = p[6][l];
                trilist[nnew + 1].p[3][l] = p[0][l];

                /* now the ones that are on the old side */
                trilist[nnew + 2].p[0][l] = p[1][l];
                trilist[nnew + 2].p[1][l] = p[2][l];
                trilist[nnew + 2].p[2][l] = p[6][l];
                trilist[nnew + 2].p[3][l] = p[4][l];

                trilist[nnew + 3].p[0][l] = p[3][l];
                trilist[nnew + 3].p[1][l] = p[1][l];
                trilist[nnew + 3].p[2][l] = p[6][l];
                trilist[nnew + 3].p[3][l] = p[4][l];

                trilist[nnew + 4].p[0][l] = p[3][l];
                trilist[nnew + 4].p[1][l] = p[6][l];
                trilist[nnew + 4].p[2][l] = p[7][l];
                trilist[nnew + 4].p[3][l] = p[4][l];
              }

            trilist[i].owner        = q;
            trilist[nnew].owner     = q;
            trilist[nnew + 1].owner = q;

            trilist[nnew + 2].owner = oldq;
            trilist[nnew + 3].owner = oldq;
            trilist[nnew + 4].owner = oldq;

            nadd = 5;

            vvi = fabs(get_tri_volume(i, trilist));
            for(l = 0; l < nadd; l++)
              vv[l] = fabs(get_tri_volume(nnew + l, trilist));

            /* determine largest */
            ilargest = i;
            vlargest = vvi;
            for(l = 0; l < nadd; l++)
              if(vv[l] > vlargest)
                {
                  vlargest = vv[l];
                  ilargest = nnew + l;
                }
            if(i != ilargest)
              {
                /* swap the largest to location i */
                triangle trisave  = trilist[i];
                trilist[i]        = trilist[ilargest];
                trilist[ilargest] = trisave;

                vv[ilargest - nnew] = vvi;
              }

            for(l = 0; l < nadd; l++)
              {
                if(vv[l] < MIN_VOL_FAC * vol)
                  {
                    vv[l]             = vv[nadd - 1];
                    trilist[nnew + l] = trilist[nnew + nadd - 1];
                    l--;
                    nadd--;
                  }
              }

            nnew += nadd;
            break;
        }
    }

  return nnew;
}

/*! \brief Processes edge for volume calculation.
 *
 *  Calculates the contribution of edge to volumes of neighboring
 *  Voronoi cells in vol array.
 *
 *  \param[in] T Pointer to tesselation.
 *  \param[in, out] volume of tetrahedra.
 *  \param[in] tt Index of triangle in DT array.
 *  \param[in] nr Index in edge array.
 *
 *  \return void
 */
void derefine_refine_process_edge(tessellation *T, double *vol, int tt, int nr)
{
  tetra *DT         = T->DT;
  point *DP         = T->DP;
  tetra_center *DTC = T->DTC;

  int i, j, k, l, m, ii, jj, kk, ll, nn, count, nr_next, p1, p2;
  tetra *prev, *next;
  tetra_center *prevc, *nextc;
  double ax, ay, az;
  double bx, by, bz;
  double cx, cy, cz;
  double nx, ny, nz;
  double hhx, hhy, hhz;
  double darea, dvol, h;

  tetra *t = &DT[tt];

  i = edge_start[nr];
  j = edge_end[nr];
  k = edge_opposite[nr];
  l = edge_nexttetra[nr];

  Edge_visited[tt] |= (1 << nr);

  p1 = t->p[i];
  p2 = t->p[j];

  double area = 0;

  cx = DTC[tt].cx;
  cy = DTC[tt].cy;
  cz = DTC[tt].cz;

  count = 0;

  prev  = t;
  prevc = &DTC[tt];
  do
    {
      nn    = prev->t[l];
      next  = &DT[nn];
      nextc = &DTC[nn];

      if(prev != t && next != t)
        {
          ax = prevc->cx - cx;
          ay = prevc->cy - cy;
          az = prevc->cz - cz;

          bx = nextc->cx - cx;
          by = nextc->cy - cy;
          bz = nextc->cz - cz;

          nx = ay * bz - az * by;
          ny = az * bx - ax * bz;
          nz = ax * by - ay * bx;

          darea = 0.5 * sqrt(nx * nx + ny * ny + nz * nz);
          area += darea;
        }

      for(m = 0, ll = ii = jj = -1; m < 4; m++)
        {
          if(next->p[m] == prev->p[k])
            ll = m;
          if(next->p[m] == prev->p[i])
            ii = m;
          if(next->p[m] == prev->p[j])
            jj = m;
        }

      if(ll < 0 || ii < 0 || jj < 0)
        terminate("inconsistency");

      kk = 6 - (ll + ii + jj);

      /* need to determine the edge number to be able to flag it */

      for(nr_next = 0; nr_next < 6; nr_next++)
        if((edge_start[nr_next] == ii && edge_end[nr_next] == jj) || (edge_start[nr_next] == jj && edge_end[nr_next] == ii))
          {
            if((Edge_visited[nn] & (1 << nr_next)) && next != t)
              terminate("inconsistency");

            Edge_visited[nn] |= (1 << nr_next);
            break;
          }

      prev  = next;
      prevc = nextc;
      i     = ii;
      l     = ll;
      j     = jj;
      k     = kk;

      count++;

      if(count > 1000)
        terminate("count is too large");
    }
  while(next != t);

  i = edge_start[nr];
  j = edge_end[nr];

  hhx = 0.5 * (DP[p1].x - DP[p2].x);
  hhy = 0.5 * (DP[p1].y - DP[p2].y);
  hhz = 0.5 * (DP[p1].z - DP[p2].z);

  h    = sqrt(hhx * hhx + hhy * hhy + hhz * hhz);
  dvol = (1.0 / 3) * area * h;

  if(p1 >= 0 && p1 < DeRefMesh.Ndp)
    vol[p1] += dvol;

  if(p2 >= 0 && p2 < DeRefMesh.Ndp)
    vol[p2] += dvol;
}

/*! \brief Insert a point into mesh.
 *
 *  Finds the tetrahedron that contains this point, splits the tetrahedron.
 *  After this, flip the edges if needed restore Delaunayhood (which is applied
 *  recursively) until a valid Delaunay mesh is restored.
 *
 *  \param[in, out] T Pointer to tessellation.
 *  \param[in] pp index of Delaunay point in DP array.
 *  \param[in] ttstart initial guess in which triangle it might be,
 *             index in DT array.
 *
 * \return index to tetra that (currently) contains the point pp.
 */
int insert_point(tessellation *T, int pp, int ttstart)
{
  int tt0, tt1, tt2, tt3, tt4, tetra_with_p, tt;
  int to_check[STACKSIZE_TETRA], freestack[STACKSIZE_TETRA];
  int n_faces_to_check = 0, nfree_on_stack = 0, moves;
  int tip_index, flag, edgeface_nr;
  int non_convex, convex_edge = 0, i, j;

  /* first, need to do a point location */
  tt0 = get_tetra(T, &T->DP[pp], &moves, ttstart, &flag, &edgeface_nr);

  tetra_with_p = tt0;

  if(flag == 1) /* that's the normal split of a tetrahedron into 4 */
    {
      if(n_faces_to_check >= STACKSIZE_TETRA - 4)
        terminate("stacksize exceeded");

      /* we now need to split this tetrahedron into four  */
      if(nfree_on_stack)
        tt1 = freestack[--nfree_on_stack];
      else
        tt1 = T->Ndt++;

      if(nfree_on_stack)
        tt2 = freestack[--nfree_on_stack];
      else
        tt2 = T->Ndt++;

      if(nfree_on_stack)
        tt3 = freestack[--nfree_on_stack];
      else
        tt3 = T->Ndt++;

      if(T->Ndt > T->MaxNdt)
        {
          T->Indi.AllocFacNdt *= ALLOC_INCREASE_FACTOR;
          T->MaxNdt = T->Indi.AllocFacNdt;
#ifdef VERBOSE
          printf("Task=%d: increase memory allocation, MaxNdt=%d Indi.AllocFacNdt=%g\n", ThisTask, T->MaxNdt, T->Indi.AllocFacNdt);
#endif /* #ifdef VERBOSE */
          T->DT  = myrealloc_movable(T->DT, T->MaxNdt * sizeof(tetra));
          T->DTC = myrealloc_movable(T->DTC, T->MaxNdt * sizeof(tetra_center));
          T->DTF = myrealloc_movable(T->DTF, T->MaxNdt * sizeof(char));

          if(T->Ndt > T->MaxNdt)
            terminate("Ndt > MaxNdt");
        }

      make_a_1_to_4_flip(T, pp, tt0, tt1, tt2, tt3);

      /* now we have a triangulation again - need to check whether there are
         facets that are not Delaunay */
      /* let's initialize a stack with the facets that we need to check */

      n_faces_to_check = 0;

      to_check[n_faces_to_check++] = tt0;
      to_check[n_faces_to_check++] = tt1;
      to_check[n_faces_to_check++] = tt2;
      to_check[n_faces_to_check++] = tt3;
      char *DTF                    = T->DTF;
      DTF[tt0]                     = 0;
      DTF[tt1]                     = 0;
      DTF[tt2]                     = 0;
      DTF[tt3]                     = 0;
    }

  if(flag == 2)
    {
      /* create four new tetra  */
      if(nfree_on_stack)
        tt1 = freestack[--nfree_on_stack];
      else
        tt1 = T->Ndt++;

      if(nfree_on_stack)
        tt2 = freestack[--nfree_on_stack];
      else
        tt2 = T->Ndt++;

      if(nfree_on_stack)
        tt3 = freestack[--nfree_on_stack];
      else
        tt3 = T->Ndt++;

      if(nfree_on_stack)
        tt4 = freestack[--nfree_on_stack];
      else
        tt4 = T->Ndt++;

      if(T->Ndt > T->MaxNdt)
        {
          T->Indi.AllocFacNdt *= ALLOC_INCREASE_FACTOR;
          T->MaxNdt = T->Indi.AllocFacNdt;
#ifdef VERBOSE
          printf("Task=%d: increase memory allocation, MaxNdt=%d Indi.AllocFacNdt=%g\n", ThisTask, T->MaxNdt, T->Indi.AllocFacNdt);
#endif /* #ifdef VERBOSE */
          T->DT  = myrealloc_movable(T->DT, T->MaxNdt * sizeof(tetra));
          T->DTC = myrealloc_movable(T->DTC, T->MaxNdt * sizeof(tetra_center));
          T->DTF = myrealloc_movable(T->DTF, T->MaxNdt * sizeof(char));

          if(T->Ndt > T->MaxNdt)
            terminate("Ndt > MaxNdt");
        }

      n_faces_to_check = 0;

      to_check[n_faces_to_check++] = tt0;
      to_check[n_faces_to_check++] = T->DT[tt0].t[edgeface_nr];
      to_check[n_faces_to_check++] = tt1;
      to_check[n_faces_to_check++] = tt2;
      to_check[n_faces_to_check++] = tt3;
      to_check[n_faces_to_check++] = tt4;

      char *DTF                      = T->DTF;
      DTF[tt0]                       = 0;
      DTF[T->DT[tt0].t[edgeface_nr]] = 0;
      DTF[tt1]                       = 0;
      DTF[tt2]                       = 0;
      DTF[tt3]                       = 0;
      DTF[tt4]                       = 0;

      make_a_face_split(T, tt0, edgeface_nr, pp, tt1, tt2, tt3, tt4);
    }

  if(flag == 3) /* here we need to split an edge */
    {
      int i, j, k, l, ii, jj, kk, ll, m, count;
      int prev, next;

      /* count how many triangles share the edge */
      i = edge_start[edgeface_nr];
      j = edge_end[edgeface_nr];
      k = edge_opposite[edgeface_nr];
      l = edge_nexttetra[edgeface_nr];

      count            = 0;
      n_faces_to_check = 0;

      prev = tt0;
      do
        {
          to_check[n_faces_to_check++] = prev;
          T->DTF[prev]                 = 0;

          tetra *DT = T->DT;
          next      = DT[prev].t[l];

          for(m = 0, ll = ii = jj = -1; m < 4; m++)
            {
              if(DT[next].p[m] == DT[prev].p[k])
                ll = m;
              if(DT[next].p[m] == DT[prev].p[i])
                ii = m;
              if(DT[next].p[m] == DT[prev].p[j])
                jj = m;
            }

          if(ll < 0 || ii < 0 || jj < 0)
            terminate("inconsistency");

          kk = 6 - (ll + ii + jj);

          prev = next;
          i    = ii;
          l    = ll;
          j    = jj;
          k    = kk;

          count++;

          if(count > 1000)
            terminate("count exceeded");
        }
      while(next != tt0);

      int *ttlist = mymalloc_movable(&ttlist, "ttlist", count * sizeof(int));

      for(i = 0; i < count; i++)
        {
          if(nfree_on_stack)
            ttlist[i] = freestack[--nfree_on_stack];
          else
            {
              ttlist[i] = T->Ndt++;

              if(T->Ndt > T->MaxNdt)
                {
                  T->Indi.AllocFacNdt *= ALLOC_INCREASE_FACTOR;
                  T->MaxNdt = T->Indi.AllocFacNdt;
#ifdef VERBOSE
                  printf("Task=%d: increase memory allocation, MaxNdt=%d Indi.AllocFacNdt=%g\n", ThisTask, T->MaxNdt,
                         T->Indi.AllocFacNdt);
#endif /* #ifdef VERBOSE */
                  T->DT  = myrealloc_movable(T->DT, T->MaxNdt * sizeof(tetra));
                  T->DTC = myrealloc_movable(T->DTC, T->MaxNdt * sizeof(tetra_center));
                  T->DTF = myrealloc_movable(T->DTF, T->MaxNdt * sizeof(char));

                  if(T->Ndt > T->MaxNdt)
                    terminate("Ndt > MaxNdt");
                }
            }

          to_check[n_faces_to_check++] = ttlist[i];
          T->DTF[ttlist[i]]            = 0;
        }

      make_an_edge_split(T, tt0, edgeface_nr, count, pp, ttlist);

      myfree(ttlist);
    }

  int iter = 0;

  while(n_faces_to_check)
    {
      iter++;
      if(iter > 200000)
        terminate("too many iterations");

      tt = to_check[--n_faces_to_check]; /* this is the current tetra to look at.
                                            The facet in question lies opposite to q */
      if(T->DT[tt].t[0] < 0)             /* deleted? */
        continue;

      for(tip_index = 0; tip_index < 4; tip_index++)
        if(T->DT[tt].p[tip_index] == pp)
          break;

      if(tip_index < 4) /* otherwise the facet has been removed in a 3-2 flip */
        {
          tetra *DT = T->DT;
          point *DP = T->DP;
          int qq    = DT[tt].t[tip_index];           /* tetrahedron that's opposite of ours and shares the facet */
          int ppp   = DT[qq].p[DT[tt].s[tip_index]]; /* point that's opposite of the facet in the other tetrahedron */

          int ret, ret_exact;

          ret = InSphere_Errorbound(&DP[DT[qq].p[0]], &DP[DT[qq].p[1]], &DP[DT[qq].p[2]], &DP[DT[qq].p[3]], &DP[pp]);
          CountInSphereTests++;

          if(ret != 0)
            ret_exact = ret;
          else
            {
              // let's decide with exact integer arithmetic
              ret_exact = InSphere_Exact(&DP[DT[qq].p[0]], &DP[DT[qq].p[1]], &DP[DT[qq].p[2]], &DP[DT[qq].p[3]], &DP[pp]);
              CountInSphereTestsExact++;
            }

          if(ret_exact > 0) /* facet is illegal, because point lies inside */
            {
              /* let's see whether the point lies in the triangle, or on a side, or opposite of one convex edge */

              non_convex = convex_edge_test(T, tt, tip_index, &convex_edge);

              if(non_convex == 0) /* we can make a 2-3 flip */
                {
                  int ww;

                  if(nfree_on_stack)
                    ww = freestack[--nfree_on_stack];
                  else
                    ww = T->Ndt++;

                  if(T->Ndt > T->MaxNdt)
                    {
                      T->Indi.AllocFacNdt *= ALLOC_INCREASE_FACTOR;
                      T->MaxNdt = T->Indi.AllocFacNdt;
#ifdef VERBOSE
                      printf("Task=%d: increase memory allocation, MaxNdt=%d Indi.AllocFacNdt=%g\n", ThisTask, T->MaxNdt,
                             T->Indi.AllocFacNdt);
#endif /* #ifdef VERBOSE */
                      T->DT  = myrealloc_movable(T->DT, T->MaxNdt * sizeof(tetra));
                      T->DTC = myrealloc_movable(T->DTC, T->MaxNdt * sizeof(tetra_center));
                      T->DTF = myrealloc_movable(T->DTF, T->MaxNdt * sizeof(char));

                      if(T->Ndt > T->MaxNdt)
                        terminate("Ndt > MaxNdt");
                    }

                  if(n_faces_to_check >= STACKSIZE_TETRA - 3)
                    terminate("stacksize exceeded");

                  make_a_2_to_3_flip(T, tt, tip_index, qq, T->DT[tt].s[tip_index], ppp, ww);

                  to_check[n_faces_to_check++] = tt;
                  to_check[n_faces_to_check++] = qq;
                  to_check[n_faces_to_check++] = ww;
                  T->DTF[tt]                   = 0;
                  T->DTF[qq]                   = 0;
                  T->DTF[ww]                   = 0;
                }
              else if(non_convex == 1) /* we might be able to make a 3-2 flip, or we deal with a convex edge on the outer hull */
                {
                  /* test whether the reflex edge is surrounded by exactly three tetrahedra */

                  i = convex_edge + 2;
                  if(i >= 3)
                    i -= 3;
                  i = access_triangles[tip_index][i];

                  for(j = 0; j < 4; j++)
                    if(DT[tt].p[i] == DT[qq].p[j])
                      break;

                  if(j >= 4)
                    {
                      terminate("not found");
                    }

                  if(DT[tt].t[i] == DT[qq].t[j]) /* this means there is exactly one tetrahedron between them, i.e. we have found the
                                                    third partner for the flip */
                    {
                      int ww;

                      ww = DT[tt].t[i];

                      make_a_3_to_2_flip(T, tt, qq, ww, tip_index, convex_edge, DT[tt].s[tip_index]);

                      DT[ww].t[0] = -1; /* mark as deleted */

                      if(nfree_on_stack < STACKSIZE_TETRA)
                        freestack[nfree_on_stack++] = ww;
                      else
                        terminate("stack full");

                      tetra_with_p = tt;
                      if(n_faces_to_check >= STACKSIZE_TETRA - 2)
                        terminate("stack too full");

                      to_check[n_faces_to_check++] = tt;
                      to_check[n_faces_to_check++] = qq;
                      T->DTF[tt]                   = 0;
                      T->DTF[qq]                   = 0;
                    }
                  else
                    {
                      if(DT[DT[tt].t[i]].p[DT[tt].s[i]] == DPinfinity && DT[DT[qq].t[j]].p[DT[qq].s[j]] == DPinfinity)
                        {
                          printf("convex edge between points=%d %d on outer hull found\n",
                                 (int)(DT[tt].p[access_triangles[tip_index][convex_edge]]),
                                 (int)(DT[tt].p[access_triangles[tip_index][convex_edge < 2 ? convex_edge + 1 : 0]]));

                          terminate("inconsistency"); /* this should not occur since we have embedded the points into a convex big
                                                         triangle */
                        }
                    }
                }
              else if(non_convex == 2) /* we might be able to make a 4-4 flip */
                {
                  i = convex_edge + 2;
                  if(i >= 3)
                    i -= 3;
                  i = access_triangles[tip_index][i]; /* this is the point opposite of edge (but not tip) */

                  tetra *DT = T->DT;
                  char *DTF = T->DTF;

                  for(j = 0; j < 4; j++)
                    if(DT[tt].p[i] == DT[qq].p[j])
                      break;

                  if(DT[DT[tt].t[i]].p[DT[tt].s[i]] == DT[DT[qq].t[j]].p[DT[qq].s[j]])
                    {
                      /* ok, so we really have 4 tetra. The opposite points match up */

                      to_check[n_faces_to_check++] = tt;
                      to_check[n_faces_to_check++] = qq;
                      to_check[n_faces_to_check++] = DT[tt].t[i];
                      to_check[n_faces_to_check++] = DT[qq].t[j];
                      DTF[tt]                      = 0;
                      DTF[qq]                      = 0;
                      DTF[DT[tt].t[i]]             = 0;
                      DTF[DT[qq].t[j]]             = 0;

                      make_a_4_to_4_flip(T, tt, tip_index, convex_edge);
                    }
                }
            }
          else
            tetra_with_p = tt;
        }
    }

  return tetra_with_p;
}

/*! \brief Tests edges and detects if a flip is needed.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] tt Index in DT array.
 *  \param[in] tip Index of forth point (tip of tetrahedron).
 *  \param[out] edgenr Index of edge.
 *
 *  \return (-1,0,1,2), depending on which flip is necessary.
 */
int convex_edge_test(tessellation *T, int tt, int tip, int *edgenr)
{
  tetra *DT = T->DT;
  point *DP = T->DP;
  tetra *t  = &DT[tt];
  int i0, i1, i2, i3;
  int vol, flag0, flag1, flag2;
  int count_zeros = 0;

  i0 = access_triangles[tip][0];
  i1 = access_triangles[tip][1];
  i2 = access_triangles[tip][2];
  i3 = tip;

  point *p0 = &DP[t->p[i0]];
  point *p1 = &DP[t->p[i1]];
  point *p2 = &DP[t->p[i2]];
  point *p3 = &DP[t->p[i3]];
  point *p4 = &DP[DT[t->t[i3]].p[t->s[i3]]];

  CountConvexEdgeTest++;

#ifndef OPTIMIZE_MEMORY_USAGE
  double ax = p1->xx - p0->xx;
  double ay = p1->yy - p0->yy;
  double az = p1->zz - p0->zz;

  double bx = p2->xx - p0->xx;
  double by = p2->yy - p0->yy;
  double bz = p2->zz - p0->zz;

  double cx = p3->xx - p0->xx;
  double cy = p3->yy - p0->yy;
  double cz = p3->zz - p0->zz;

  double qx = p4->xx - p0->xx;
  double qy = p4->yy - p0->yy;
  double qz = p4->zz - p0->zz;
#else  /* #ifndef OPTIMIZE_MEMORY_USAGE */
  double ax, ay, az, bx, by, bz, cx, cy, cz, qx, qy, qz;
  double pA_xyz[3], pB_xyz[3];
  IntegerMapType pA_ixyz[3], pB_ixyz[3];

  get_integers_for_point(p0, pA_ixyz, pA_xyz);

  get_integers_for_point(p1, pB_ixyz, pB_xyz);
  ax = pB_xyz[0] - pA_xyz[0];
  ay = pB_xyz[1] - pA_xyz[1];
  az = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p2, pB_ixyz, pB_xyz);
  bx = pB_xyz[0] - pA_xyz[0];
  by = pB_xyz[1] - pA_xyz[1];
  bz = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p3, pB_ixyz, pB_xyz);
  cx = pB_xyz[0] - pA_xyz[0];
  cy = pB_xyz[1] - pA_xyz[1];
  cz = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p4, pB_ixyz, pB_xyz);
  qx = pB_xyz[0] - pA_xyz[0];
  qy = pB_xyz[1] - pA_xyz[1];
  qz = pB_xyz[2] - pA_xyz[2];
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE */

  double mv_data[] = {ax, bx, cx, qx, ay, by, cy, qy, az, bz, cz, qz};
  double x[3];

  int status;

  status = solve_linear_equations(mv_data, x);

  /* x now contains the coordinates of the point p4 expanded in the basis (a,b,c) */
  /* the coordinates of point 3 in this basis are (0,0,1) */

  if(status >= 0)
    {
      if(fabs(1.0 - x[2]) < INSIDE_EPS)
        terminate("inconsistency");

      double u, v, w;

      w = 1.0 / (1.0 - x[2]);

      u = w * x[0];
      v = w * x[1];

      if(u > INSIDE_EPS && v > INSIDE_EPS && (1 - (u + v)) > INSIDE_EPS)
        {
          /* we have a point safely in the triangle: 2-3 flip should be fine */
          return 0;
        }

      if(u > INSIDE_EPS && v < -INSIDE_EPS && (1 - (u + v)) > INSIDE_EPS)
        {
          /* edge 0 is clearly reflect,  3-2 flip allowed around edge 0 */
          *edgenr = 0;
          return 1;
        }

      if(u > INSIDE_EPS && v > INSIDE_EPS && (1 - (u + v)) < -INSIDE_EPS)
        {
          // printf("3-2 flip allowed since edge 1 is reflex\n");
          *edgenr = 1;
          return 1;
        }

      if(u < -INSIDE_EPS && v > INSIDE_EPS && (1 - (u + v)) > INSIDE_EPS)
        {
          // printf("3-2 flip allowed since edge 2 is reflex\n");
          *edgenr = 2;
          return 1;
        }

      if(u < -INSIDE_EPS && v < -INSIDE_EPS && (1 - (u + v)) > INSIDE_EPS)
        return -1; /* two reflex edges */

      if(u < -INSIDE_EPS && v > INSIDE_EPS && (1 - (u + v)) < -INSIDE_EPS)
        return -1; /* two reflex edges */

      if(u > INSIDE_EPS && v < -INSIDE_EPS && (1 - (u + v)) < -INSIDE_EPS)
        return -1; /* two reflex edges */
    }

  CountConvexEdgeTestExact++;

  /* Now we need to test in more detail if we are on one of the edges */

  vol = Orient3d_Exact(p0, p1, p2, p3);

  if(vol <= 0)
    {
      printf("flat or negatively tetrahedron found (vol=%d)\n", vol);
      {
        printf("p0=%d  %g %g %g\n", (int)(p0 - DP), p0->x, p0->y, p0->z);
        printf("p1=%d  %g %g %g\n", (int)(p1 - DP), p1->x, p1->y, p1->z);
        printf("p2=%d  %g %g %g\n", (int)(p2 - DP), p2->x, p2->y, p2->z);
        printf("p3=%d  %g %g %g\n", (int)(p3 - DP), p3->x, p3->y, p3->z);
        dump_points(T);
        terminate("inconsistent tetrahedron");
      }
    }

  flag0 = Orient3d_Exact(p1, p3, p2, p4);
  flag1 = Orient3d_Exact(p0, p2, p3, p4);
  flag2 = Orient3d_Exact(p0, p3, p1, p4);

  if(flag0 == 0)
    count_zeros++;

  if(flag1 == 0)
    count_zeros++;

  if(flag2 == 0)
    count_zeros++;

  if(flag0 >= 0 && flag1 >= 0 && flag2 < 0)
    {
      //  printf("3-2 flip allowed since edge 0 is reflex\n");
      *edgenr = 0;
      return 1;
    }

  if(flag0 < 0 && flag1 >= 0 && flag2 >= 0)
    {
      // printf("3-2 flip allowed since edge 1 is reflex\n");
      *edgenr = 1;
      return 1;
    }

  if(flag0 >= 0 && flag1 < 0 && flag2 >= 0)
    {
      // printf("3-2 flip allowed since edge 2 is reflex\n");
      *edgenr = 2;
      return 1;
    }

  if(flag0 >= 0 && flag1 >= 0 && flag2 == 0)
    {
      // printf("4-4 flip around edge 0 may be possible\n");
      *edgenr = 0;
      return 2;
    }

  if(flag0 >= 0 && flag1 == 0 && flag2 >= 0)
    {
      // printf("4-4 flip around edge 2 may be possible\n");
      *edgenr = 2;
      return 2;
    }

  if(flag0 == 0 && flag1 >= 0 && flag2 >= 0)
    {
      // printf("4-4 flip around edge 1 may be possible\n");
      *edgenr = 1;
      return 2;
    }

  if(flag0 >= 0 && flag1 >= 0 && flag2 >= 0)
    {
      /* we seem to have a point in the triangle: 2-3 flip should be fine */
      return 0;
    }

  return -1;
}

/*! \brief Performs face split.
 *
 *  \param[in, out] T Pointer to tessellation.
 *  \param[in] tt0 First index in DT array.
 *  \param[in] face_nr Index of face.
 *  \param[in] pp Index of point.
 *  \param[in] tt1 Second index in DT array.
 *  \param[in] tt2 Third index in DT array.
 *  \param[in] qq1 Index in DT array.
 *  \param[in] qq2 Index in DT array.
 *
 *  \return void
 */
void make_a_face_split(tessellation *T, int tt0, int face_nr, int pp, int tt1, int tt2, int qq1, int qq2)
{
  tetra *DT = T->DT;
  tetra *t0 = &DT[tt0];
  tetra *t1 = &DT[tt1];
  tetra *t2 = &DT[tt2];
  int qq0   = t0->t[face_nr];
  tetra *q0 = &DT[qq0];
  tetra *q1 = &DT[qq1];
  tetra *q2 = &DT[qq2];

  int m, i0 = -1, i1 = -1, i2 = -1, i3 = -1, j0 = -1, j1 = -1, j2 = -1, j3 = -1;

  Count_FaceSplits++;
  CountFlips++;

  *t1 = *t0;
  *t2 = *t0;

  *q1 = *q0;
  *q2 = *q0;

  i3 = face_nr;
  j3 = t0->s[face_nr];

  switch(i3)
    {
      case 3:
        i0 = 0;
        i1 = 1;
        i2 = 2;
        break;
      case 2:
        i0 = 0;
        i1 = 3;
        i2 = 1;
        break;
      case 1:
        i0 = 0;
        i1 = 2;
        i2 = 3;
        break;
      case 0:
        i0 = 1;
        i1 = 3;
        i2 = 2;
        break;
    }

  for(m = 0; m < 4; m++)
    {
      if(q0->p[m] == t0->p[i0])
        j0 = m;
      if(q0->p[m] == t0->p[i1])
        j2 = m;
      if(q0->p[m] == t0->p[i2])
        j1 = m;
    }

  if(i0 < 0 || i1 < 0 || i2 < 0 || i3 < 0 || j0 < 0 || j1 < 0 || j2 < 0 || j3 < 0)
    terminate("inconsistency");

  t0->p[i2] = pp;
  t1->p[i0] = pp;
  t2->p[i1] = pp;

  q0->p[j1] = pp;
  q1->p[j0] = pp;
  q2->p[j2] = pp;

  t0->t[i0] = tt1;
  t1->t[i2] = tt0;
  t0->s[i0] = i2;
  t1->s[i2] = i0;

  t1->t[i1] = tt2;
  t2->t[i0] = tt1;
  t1->s[i1] = i0;
  t2->s[i0] = i1;

  t2->t[i2] = tt0;
  t0->t[i1] = tt2;
  t2->s[i2] = i1;
  t0->s[i1] = i2;

  q0->t[j0] = qq1;
  q1->t[j1] = qq0;
  q0->s[j0] = j1;
  q1->s[j1] = j0;

  q1->t[j2] = qq2;
  q2->t[j0] = qq1;
  q1->s[j2] = j0;
  q2->s[j0] = j2;

  q2->t[j1] = qq0;
  q0->t[j2] = qq2;
  q2->s[j1] = j2;
  q0->s[j2] = j1;

  t0->t[i3] = qq0;
  q0->t[j3] = tt0;
  t0->s[i3] = j3;
  q0->s[j3] = i3;

  t1->t[i3] = qq1;
  q1->t[j3] = tt1;
  t1->s[i3] = j3;
  q1->s[j3] = i3;

  t2->t[i3] = qq2;
  q2->t[j3] = tt2;
  t2->s[i3] = j3;
  q2->s[j3] = i3;

  DT[t0->t[i2]].t[t0->s[i2]] = tt0;
  DT[t1->t[i0]].t[t1->s[i0]] = tt1;
  DT[t2->t[i1]].t[t2->s[i1]] = tt2;

  DT[q0->t[j1]].t[q0->s[j1]] = qq0;
  DT[q1->t[j0]].t[q1->s[j0]] = qq1;
  DT[q2->t[j2]].t[q2->s[j2]] = qq2;
}

/*! \brief Performs edge split.
 *
 *  \param[in, out] T Pointer to tessellation
 *  \param[in] tt0 Index in DT array
 *  \param[in] edge_nr Index of edge
 *  \param[in] count Number of elements in lists.
 *  \param[in] pp Index to point.
 *  \param[in] ttlist List of indices in DT.
 */
void make_an_edge_split(tessellation *T, int tt0, int edge_nr, int count, int pp, int *ttlist)
{
  tetra *DT = T->DT;
  tetra *t0 = &DT[tt0];
  tetra *prev, *next;
  tetra **tlist, **t_orig_list;
  int *i_list, *j_list, *k_list, *l_list;
  int i, j, k, l, ii, jj, kk, ll, m, nr, nrm, nrp;

  Count_EdgeSplits++;
  CountFlips++;

  tlist       = mymalloc("tlist", count * sizeof(tetra *));
  t_orig_list = mymalloc("t_orig_list", count * sizeof(tetra *));
  i_list      = mymalloc("i_list", sizeof(int) * count);
  j_list      = mymalloc("j_list", sizeof(int) * count);
  k_list      = mymalloc("k_list", sizeof(int) * count);
  l_list      = mymalloc("l_list", sizeof(int) * count);

  for(i = 0; i < count; i++)
    tlist[i] = &DT[ttlist[i]];

  i = edge_start[edge_nr];
  j = edge_end[edge_nr];
  k = edge_opposite[edge_nr];
  l = edge_nexttetra[edge_nr];

  nr   = 0;
  prev = t0;
  do
    {
      t_orig_list[nr] = prev;
      i_list[nr]      = i;
      j_list[nr]      = j;
      k_list[nr]      = k;
      l_list[nr]      = l;

      next = &DT[prev->t[l]];

      for(m = 0, ll = ii = jj = -1; m < 4; m++)
        {
          if(next->p[m] == prev->p[k])
            ll = m;
          if(next->p[m] == prev->p[i])
            ii = m;
          if(next->p[m] == prev->p[j])
            jj = m;
        }

      if(ll < 0 || ii < 0 || jj < 0)
        terminate("inconsistency");

      kk = 6 - (ll + ii + jj);

      prev = next;
      i    = ii;
      l    = ll;
      j    = jj;
      k    = kk;

      nr++;
    }
  while(next != t0);

  for(nr = 0; nr < count; nr++)
    {
      *tlist[nr] = *t_orig_list[nr];

      t_orig_list[nr]->p[j_list[nr]] = pp;
      tlist[nr]->p[i_list[nr]]       = pp;

      t_orig_list[nr]->t[i_list[nr]] = tlist[nr] - DT;
      tlist[nr]->t[j_list[nr]]       = t_orig_list[nr] - DT;

      t_orig_list[nr]->s[i_list[nr]] = j_list[nr];
      tlist[nr]->s[j_list[nr]]       = i_list[nr];

      DT[tlist[nr]->t[i_list[nr]]].t[tlist[nr]->s[i_list[nr]]] = tlist[nr] - DT;

      nrp = nr + 1;
      if(nrp >= count)
        nrp -= count;

      nrm = nr - 1;
      if(nrm < 0)
        nrm += count;

      tlist[nr]->t[l_list[nr]] = tlist[nrp] - DT;
      tlist[nr]->s[l_list[nr]] = k_list[nrp];

      tlist[nr]->t[k_list[nr]] = tlist[nrm] - DT;
      tlist[nr]->s[k_list[nr]] = l_list[nrm];
    }

  myfree(l_list);
  myfree(k_list);
  myfree(j_list);
  myfree(i_list);

  myfree(t_orig_list);
  myfree(tlist);
}

/*! \brief Make a 4 to 4 flip.
 *
 *  See Springel (2010) for discussion on flips.
 *
 *  \param[in, out] T Pointer to tessellation.
 *  \param[in] tt Index in DT array.
 *  \param[in] tip_index Index of the point making up the tip of the
 *             tetrahedron.
 *  \param[in] edge_nr Index of edge.
 *
 *  \return void
 */
void make_a_4_to_4_flip(tessellation *T, int tt, int tip_index, int edge_nr)
{
  tetra *DT = T->DT;
  //  printf("4-to-4 flip\n");
  tetra *t = &DT[tt];
  int i0, i1, i2, j;
  int ww, qq, uu;
  tetra *w, *q, *u;
  tetra *t_top[4], *t_bottom[4];
  int s_top[4], s_bottom[4];
  int p[6];

  Count_4_to_4_Flips++;
  CountFlips++;

  uu = 0;
  u  = NULL;

  for(j = 0; j < 4; j++)
    {
      t_top[j]    = NULL;
      t_bottom[j] = NULL;
      s_top[j]    = -1;
      s_bottom[j] = -1;
    }

  i0 = access_triangles[tip_index][edge_nr];
  edge_nr += 1;
  if(edge_nr >= 3)
    edge_nr -= 3;
  i1 = access_triangles[tip_index][edge_nr];
  edge_nr += 1;
  if(edge_nr >= 3)
    edge_nr -= 3;
  i2 = access_triangles[tip_index][edge_nr];

  t_top[0] = &DT[t->t[i0]];
  s_top[0] = t->s[i0];

  t_top[1] = &DT[t->t[i1]];
  s_top[1] = t->s[i1];

  ww = t->t[i2];
  w  = &DT[ww];
  qq = t->t[tip_index];
  q  = &DT[qq];

  for(j = 0; j < 4; j++)
    {
      if(w->p[j] == t->p[i0])
        {
          t_top[3] = &DT[w->t[j]];
          s_top[3] = w->s[j];
        }

      if(w->p[j] == t->p[i1])
        {
          t_top[2] = &DT[w->t[j]];
          s_top[2] = w->s[j];
        }

      if(w->p[j] == t->p[tip_index])
        {
          uu = w->t[j];
          u  = &DT[uu];
        }
    }

  for(j = 0; j < 4; j++)
    {
      if(u->p[j] == t->p[i0])
        {
          t_bottom[3] = &DT[u->t[j]];
          s_bottom[3] = u->s[j];
        }

      if(u->p[j] == t->p[i1])
        {
          t_bottom[2] = &DT[u->t[j]];
          s_bottom[2] = u->s[j];
        }

      if(q->p[j] == t->p[i0])
        {
          t_bottom[0] = &DT[q->t[j]];
          s_bottom[0] = q->s[j];
        }

      if(q->p[j] == t->p[i1])
        {
          t_bottom[1] = &DT[q->t[j]];
          s_bottom[1] = q->s[j];
        }
    }

  p[0] = t->p[i1];
  p[1] = t->p[i2];
  p[2] = t->p[i0];
  p[3] = DT[t->t[i2]].p[t->s[i2]];
  p[4] = t->p[tip_index];
  p[5] = DT[t->t[tip_index]].p[t->s[tip_index]];

  for(j = 0; j < 4; j++)
    {
      if(t_top[j] == NULL || t_bottom[j] == NULL)
        {
          printf("bad!\n");
          terminate("inconsistency");
        }
    }

  for(j = 0; j < 4; j++)
    {
      if(t_top[j] == NULL || t_bottom[j] == NULL)
        {
          printf("bad!\n");
          terminate("inconsistency");
        }
    }

  t->p[0] = p[0];
  t->p[1] = p[1];
  t->p[2] = p[5];
  t->p[3] = p[4];

  q->p[0] = p[1];
  q->p[1] = p[2];
  q->p[2] = p[5];
  q->p[3] = p[4];

  u->p[0] = p[2];
  u->p[1] = p[3];
  u->p[2] = p[5];
  u->p[3] = p[4];

  w->p[0] = p[3];
  w->p[1] = p[0];
  w->p[2] = p[5];
  w->p[3] = p[4];

  t->t[0] = qq;
  q->t[1] = tt;
  t->s[0] = 1;
  q->s[1] = 0;

  q->t[0] = uu;
  u->t[1] = qq;
  q->s[0] = 1;
  u->s[1] = 0;

  u->t[0] = ww;
  w->t[1] = uu;
  u->s[0] = 1;
  w->s[1] = 0;

  w->t[0] = tt;
  t->t[1] = ww;
  w->s[0] = 1;
  t->s[1] = 0;

  t->t[2]                = t_top[0] - DT;
  t->s[2]                = s_top[0];
  DT[t->t[2]].t[t->s[2]] = tt;
  DT[t->t[2]].s[t->s[2]] = 2;

  t->t[3]                = t_bottom[0] - DT;
  t->s[3]                = s_bottom[0];
  DT[t->t[3]].t[t->s[3]] = tt;
  DT[t->t[3]].s[t->s[3]] = 3;

  q->t[2]                = t_top[1] - DT;
  q->s[2]                = s_top[1];
  DT[q->t[2]].t[q->s[2]] = qq;
  DT[q->t[2]].s[q->s[2]] = 2;

  q->t[3]                = t_bottom[1] - DT;
  q->s[3]                = s_bottom[1];
  DT[q->t[3]].t[q->s[3]] = qq;
  DT[q->t[3]].s[q->s[3]] = 3;

  u->t[2]                = t_top[2] - DT;
  u->s[2]                = s_top[2];
  DT[u->t[2]].t[u->s[2]] = uu;
  DT[u->t[2]].s[u->s[2]] = 2;

  u->t[3]                = t_bottom[2] - DT;
  u->s[3]                = s_bottom[2];
  DT[u->t[3]].t[u->s[3]] = uu;
  DT[u->t[3]].s[u->s[3]] = 3;

  w->t[2]                = t_top[3] - DT;
  w->s[2]                = s_top[3];
  DT[w->t[2]].t[w->s[2]] = ww;
  DT[w->t[2]].s[w->s[2]] = 2;

  w->t[3]                = t_bottom[3] - DT;
  w->s[3]                = s_bottom[3];
  DT[w->t[3]].t[w->s[3]] = ww;
  DT[w->t[3]].s[w->s[3]] = 3;
}

/*! \brief Make a 1 to 4 flip.
 *
 *  See Springel (2010) for discussion on flips.
 *
 *  \param[in, out] T Pointer to tessellation.
 *  \param[in] pp Index of new point.
 *  \param[in] tt0 Index or first point in DT array.
 *  \param[in] tt1 Index of second point in DT array.
 *  \param[in] tt2 Index of third point in DT array.
 *  \param[in] tt3 Index of forth point in DT array.
 *
 *  \return void
 */
void make_a_1_to_4_flip(tessellation *T, int pp, int tt0, int tt1, int tt2, int tt3)
{
  tetra *DT = T->DT;

  tetra *t0 = &DT[tt0];
  tetra *t1 = &DT[tt1];
  tetra *t2 = &DT[tt2];
  tetra *t3 = &DT[tt3];

  Count_1_to_4_Flips++;
  CountFlips++;

  *t1 = *t0;
  *t2 = *t0;
  *t3 = *t0;

  t0->p[0] = pp;
  t1->p[1] = pp;
  t2->p[2] = pp;
  t3->p[3] = pp;

  t0->t[1] = tt1;
  t1->t[0] = tt0;
  t0->s[1] = 0;
  t1->s[0] = 1;

  t1->t[2] = tt2;
  t2->t[1] = tt1;
  t1->s[2] = 1;
  t2->s[1] = 2;

  t2->t[0] = tt0;
  t0->t[2] = tt2;
  t2->s[0] = 2;
  t0->s[2] = 0;

  t0->t[3] = tt3;
  t3->t[0] = tt0;
  t0->s[3] = 0;
  t3->s[0] = 3;

  t1->t[3] = tt3;
  t3->t[1] = tt1;
  t1->s[3] = 1;
  t3->s[1] = 3;

  t2->t[3] = tt3;
  t3->t[2] = tt2;
  t2->s[3] = 2;
  t3->s[2] = 3;

  DT[t0->t[0]].t[t0->s[0]] = tt0;
  DT[t1->t[1]].t[t1->s[1]] = tt1;
  DT[t2->t[2]].t[t2->s[2]] = tt2;
  DT[t3->t[3]].t[t3->s[3]] = tt3;
}

/*! \brief Make a 3 to 2 flip.
 *
 *  See Springel (2010) for discussion on flips.
 *
 *  \param[in, out] T Pointer to tessellation.
 *  \param[in] pp Index of new point.
 *  \param[in] tt0 Index or first point in DT array.
 *  \param[in] tt1 Index of second point in DT array.
 *  \param[in] tt2 Index of third point in DT array.
 *  \param[in] tip Index of point making up tip of tetrahedron.
 *  \param[in] edge Index of edge.
 *  \param[in] bottom Tetrahedron on bottom.
 *
 *  \return void
 */
void make_a_3_to_2_flip(tessellation *T, int tt0, int tt1, int tt2, int tip, int edge, int bottom)
{
  tetra *DT = T->DT;
  tetra *t0 = &DT[tt0];
  tetra *t1 = &DT[tt1];
  tetra *t2 = &DT[tt2];

  int i, j, k, ii, jj, iii, jjj;
  tetra qbak, tbak, wbak;

  Count_3_to_2_Flips++;
  CountFlips++;

  tbak = *t0;
  qbak = *t1;
  wbak = *t2;

  i = edge;
  j = i + 1;
  k = i + 2;
  if(j >= 3)
    j -= 3;
  if(k >= 3)
    k -= 3;

  i = access_triangles[tip][i];
  j = access_triangles[tip][j];
  k = access_triangles[tip][k];

  for(ii = 0; ii < 4; ii++)
    if(tbak.p[i] == qbak.p[ii])
      break;

  for(iii = 0; iii < 4; iii++)
    if(tbak.p[i] == wbak.p[iii])
      break;

  for(jj = 0; jj < 4; jj++)
    if(tbak.p[j] == qbak.p[jj])
      break;

  for(jjj = 0; jjj < 4; jjj++)
    if(tbak.p[j] == wbak.p[jjj])
      break;

  t0->p[0] = qbak.p[bottom];
  t0->p[1] = tbak.p[k];
  t0->p[2] = tbak.p[i];
  t0->p[3] = tbak.p[tip];

  t1->p[0] = qbak.p[bottom];
  t1->p[1] = tbak.p[j];
  t1->p[2] = tbak.p[k];
  t1->p[3] = tbak.p[tip];

  t0->t[2] = tt1;
  t1->t[1] = tt0;
  t0->s[2] = 1;
  t1->s[1] = 2;

  t0->t[0]                 = tbak.t[j];
  t0->s[0]                 = tbak.s[j];
  DT[t0->t[0]].s[t0->s[0]] = 0;
  DT[t0->t[0]].t[t0->s[0]] = tt0;

  t0->t[3]                 = qbak.t[jj];
  t0->s[3]                 = qbak.s[jj];
  DT[t0->t[3]].s[t0->s[3]] = 3;
  DT[t0->t[3]].t[t0->s[3]] = tt0;

  t0->t[1]                 = wbak.t[jjj];
  t0->s[1]                 = wbak.s[jjj];
  DT[t0->t[1]].s[t0->s[1]] = 1;
  DT[t0->t[1]].t[t0->s[1]] = tt0;

  t1->t[0]                 = tbak.t[i];
  t1->s[0]                 = tbak.s[i];
  DT[t1->t[0]].s[t1->s[0]] = 0;
  DT[t1->t[0]].t[t1->s[0]] = tt1;

  t1->t[3]                 = qbak.t[ii];
  t1->s[3]                 = qbak.s[ii];
  DT[t1->t[3]].s[t1->s[3]] = 3;
  DT[t1->t[3]].t[t1->s[3]] = tt1;

  t1->t[2]                 = wbak.t[iii];
  t1->s[2]                 = wbak.s[iii];
  DT[t1->t[2]].s[t1->s[2]] = 2;
  DT[t1->t[2]].t[t1->s[2]] = tt1;

  CountFlips++;
}

/*! \brief Make a 2 to 3 flip
 *
 *  See Springel (2010) for discussion on flips.
 *
 *  \param[in, out] T Pointer to tessellation.
 *  \param[in] pp Index of new point.
 *  \param[in] tt0 Index or first point in DT array.
 *  \param[in] tip Index of point makting up tip of tetrahedron.
 *  \param[in] tt1 Index of second point in DT array.
 *  \param[in] bottom Tetrahedron on bottom.
 *  \param[in] qq Index of point.
 *  \param[in] tt2 Index of third point in DT array.
 *
 *  \return void
 */
void make_a_2_to_3_flip(tessellation *T, int tt0, int tip, int tt1, int bottom, int qq, int tt2)
{
  tetra *DT = T->DT;
  tetra *t0 = &DT[tt0];
  tetra *t1 = &DT[tt1];
  tetra *t2 = &DT[tt2];
  tetra qbak, tbak;
  int k;

  Count_2_to_3_Flips++;

  tbak = *t0;
  qbak = *t1; /* to save info */

  *t1 = *t0;
  *t2 = *t0;

  /* redefine points */
  t0->p[access_triangles[tip][0]] = qq;
  t1->p[access_triangles[tip][1]] = qq;
  t2->p[access_triangles[tip][2]] = qq;

  /* make neighbour connections */
  t0->t[access_triangles[tip][1]] = tt1;
  t1->t[access_triangles[tip][0]] = tt0;
  t0->s[access_triangles[tip][1]] = access_triangles[tip][0];
  t1->s[access_triangles[tip][0]] = access_triangles[tip][1];

  t0->t[access_triangles[tip][2]] = tt2;
  t2->t[access_triangles[tip][0]] = tt0;
  t0->s[access_triangles[tip][2]] = access_triangles[tip][0];
  t2->s[access_triangles[tip][0]] = access_triangles[tip][2];

  t1->t[access_triangles[tip][2]] = tt2;
  t2->t[access_triangles[tip][1]] = tt1;
  t1->s[access_triangles[tip][2]] = access_triangles[tip][1];
  t2->s[access_triangles[tip][1]] = access_triangles[tip][2];

  /* these are the ones on the top */
  DT[t0->t[access_triangles[tip][0]]].t[t0->s[access_triangles[tip][0]]] = tt0;
  DT[t1->t[access_triangles[tip][1]]].t[t1->s[access_triangles[tip][1]]] = tt1;
  DT[t2->t[access_triangles[tip][2]]].t[t2->s[access_triangles[tip][2]]] = tt2;

  /* now the one at the bottom */

  if(qbak.p[access_triangles[bottom][0]] == tbak.p[access_triangles[tip][0]])
    k = 0;
  else if(qbak.p[access_triangles[bottom][1]] == tbak.p[access_triangles[tip][0]])
    k = 1;
  else
    k = 2;

  t0->t[tip]                   = qbak.t[access_triangles[bottom][k]];
  t0->s[tip]                   = qbak.s[access_triangles[bottom][k]];
  DT[t0->t[tip]].t[t0->s[tip]] = tt0;
  DT[t0->t[tip]].s[t0->s[tip]] = tip;

  if(qbak.p[access_triangles[bottom][0]] == tbak.p[access_triangles[tip][1]])
    k = 0;
  else if(qbak.p[access_triangles[bottom][1]] == tbak.p[access_triangles[tip][1]])
    k = 1;
  else
    k = 2;

  t1->t[tip]                   = qbak.t[access_triangles[bottom][k]];
  t1->s[tip]                   = qbak.s[access_triangles[bottom][k]];
  DT[t1->t[tip]].t[t1->s[tip]] = tt1;
  DT[t1->t[tip]].s[t1->s[tip]] = tip;

  if(qbak.p[access_triangles[bottom][0]] == tbak.p[access_triangles[tip][2]])
    k = 0;
  else if(qbak.p[access_triangles[bottom][1]] == tbak.p[access_triangles[tip][2]])
    k = 1;
  else
    k = 2;

  t2->t[tip]                   = qbak.t[access_triangles[bottom][k]];
  t2->s[tip]                   = qbak.s[access_triangles[bottom][k]];
  DT[t2->t[tip]].t[t2->s[tip]] = tt2;
  DT[t2->t[tip]].s[t2->s[tip]] = tip;
}

static int ErrorFlag = 0;

/*! \brief Gets tetrahedron.
 *
 *  Returns the index of the tetrahedron containing the point DP[pp].
 *  The search is started from the tetrahedron DT[ttstart].
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] p Point.
 *  \param[out] moves The number of moves necessary to find tetrahedron.
 *  \param[out] flag The return value from InTetra, specifying whether
 *              the point is inside or on the edge/face.
 *  \param[out] edgeface_nr The edge/face number on the tetrahedron containing
 *              the point, in case flag is >1.
 *
 *  \return Index of tetrahedron.
 */
int get_tetra(tessellation *T, point *p, int *moves, int ttstart, int *flag, int *edgeface_nr)
{
  int ret, count_moves = 0;
  int tt, next_tetra;

  tt = ttstart;

#define MAX_COUNT_MOVES 1000000

  while((ret = InTetra(T, tt, p, edgeface_nr, &next_tetra)) == 0)
    {
      count_moves++;

      if(count_moves > MAX_COUNT_MOVES)
        {
          ErrorFlag = 1;

          if(count_moves > MAX_COUNT_MOVES + 10)
            terminate("too many moves");
        }

      tt = next_tetra;
    }

  *moves = count_moves;
  *flag  = ret;

  return tt;
}

/*! \brief Is point in tetrahedron?
 *
 *  Tests whether point DP[pp] lies in the tetrahedron DT[tt]. The
 *  return value is 0 if the point is outside, 1 if it's inside, 2 if
 *  it's on a face, and 3 if it's on an edge. If it's either of the
 *  last two, the edgeface_nr is set to the corresponding index of the
 *  edge or face. If the point is outside, nexttetra is set to the
 *  index of a neighboring tetrahedron in the direction of the
 *  point, otherwise it's unmodified.
 *
 *  \param[in] T Tesslation.
 *  \param[in] tt Index of tetrahedron in DT array.
 *  \param[in] p Point.
 *  \param[out] edgeface_nr The edge/face number on the tetrahedron containing
 *              the point, in case flag is >1.
 *  \param[out] nexttetra Index of tetrahedron.
 *
 *  \return Point in thetrahedron?
 *
 */
int InTetra(tessellation *T, int tt, point *p, int *edgeface_nr, int *nexttetra)
{
  tetra *DT = T->DT;
  point *DP = T->DP;
  tetra *t  = &DT[tt];

  point *p0 = &DP[t->p[0]];
  point *p1 = &DP[t->p[1]];
  point *p2 = &DP[t->p[2]];
  point *p3 = &DP[t->p[3]];

  // test if we are in an "infinity tetra", which are the ones that
  // bound the tesselated volume. Arepo terminates if this happens,
  // but for Sunrise this is a valid occurence so we'll return -1 to
  // indicate the point is outside the volume. XXX Actually it
  // shouldn't do this anymore because we now do box tests instead
  if(isInfinity(p0) || isInfinity(p1) || isInfinity(p2) || isInfinity(p3))
    {
#ifndef LONGIDS
      printf("task=%d: we are in a tetraeder with an infinity point. tetra=%d, coordinates of point=(%g|%g|%g) ID=%d\n", ThisTask, tt,
             p->x, p->y, p->z, p->ID);
#else  /* #ifndef LONGIDS */
      printf("task=%d: we are in a tetraeder with an infinity point. tetra=%d, coordinates of point=(%g|%g|%g) ID=%llu\n", ThisTask,
             tt, p->x, p->y, p->z, p->ID);
#endif /* #ifndef LONGIDS #else */
      terminate("invalid tetrahedron");
    }

  Count_InTetra++;

#ifndef OPTIMIZE_MEMORY_USAGE
  double ax = p1->xx - p0->xx;
  double ay = p1->yy - p0->yy;
  double az = p1->zz - p0->zz;

  double bx = p2->xx - p0->xx;
  double by = p2->yy - p0->yy;
  double bz = p2->zz - p0->zz;

  double cx = p3->xx - p0->xx;
  double cy = p3->yy - p0->yy;
  double cz = p3->zz - p0->zz;

  double qx = p->xx - p0->xx;
  double qy = p->yy - p0->yy;
  double qz = p->zz - p0->zz;
#else  /* #ifndef OPTIMIZE_MEMORY_USAGE */
  double ax, ay, az, bx, by, bz, cx, cy, cz, qx, qy, qz;
  double pA_xyz[3], pB_xyz[3];
  IntegerMapType pA_ixyz[3], pB_ixyz[3];

  get_integers_for_point(p0, pA_ixyz, pA_xyz);

  get_integers_for_point(p1, pB_ixyz, pB_xyz);
  ax = pB_xyz[0] - pA_xyz[0];
  ay = pB_xyz[1] - pA_xyz[1];
  az = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p2, pB_ixyz, pB_xyz);
  bx = pB_xyz[0] - pA_xyz[0];
  by = pB_xyz[1] - pA_xyz[1];
  bz = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p3, pB_ixyz, pB_xyz);
  cx = pB_xyz[0] - pA_xyz[0];
  cy = pB_xyz[1] - pA_xyz[1];
  cz = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p, pB_ixyz, pB_xyz);
  qx = pB_xyz[0] - pA_xyz[0];
  qy = pB_xyz[1] - pA_xyz[1];
  qz = pB_xyz[2] - pA_xyz[2];
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE #else */

  double mv_data[] = {ax, bx, cx, qx, ay, by, cy, qy, az, bz, cz, qz};
  double x[3];

  int ivol, flag3, flag2, flag1, flag0;
  int count_zeros = 0;

  int status;

  status = solve_linear_equations(mv_data, x);

  if(status < 0)
    {
      ivol = Orient3d_Exact(p0, p1, p2, p3);
      if(ivol <= 0)
        {
          printf("flat or negatively tetrahedron found (ivol=%d) tt=%d\n", ivol, tt);
          terminate("invalid tetrahedron");
        }
    }

  /* x now contains the coordinates of the point p expanded in the basis (a,b,c) */

  if(ErrorFlag)
    {
      ivol  = Orient3d_Exact(p0, p1, p2, p3);
      flag3 = Orient3d_Exact(p0, p1, p2, p);
      flag2 = Orient3d_Exact(p0, p3, p1, p);
      flag1 = Orient3d_Exact(p0, p2, p3, p);
      flag0 = Orient3d_Exact(p1, p3, p2, p);

      printf("\n\nTetra=%d\n", (int)(t - DT));
      printf("ivol=%d  flag0=%d %d %d %d\n", ivol, flag0, flag1, flag2, flag3);
      printf("xx = %g %g %g   1-sum=%g\n", x[0], x[1], x[2], 1 - (x[0] + x[1] + x[2]));
      printf("a= %g %g %g\n", ax, ay, az);
      printf("b= %g %g %g\n", bx, by, bz);
      printf("c= %g %g %g\n", cx, cy, cz);
      printf("q= %g %g %g\n", qx, qy, qz);
      printf("(axb)*c) = %g\n", (ay * bz - az * by) * cx + (az * bx - ax * bz) * cy + (ax * by - ay * bx) * cz);
      printf("next tetras=%d %d %d %d\n", t->t[0], t->t[1], t->t[2], t->t[3]);
    }

  if(status >= 0)
    {
      if(x[0] > INSIDE_EPS && x[1] > INSIDE_EPS && x[2] > INSIDE_EPS && (1 - (x[0] + x[1] + x[2])) > INSIDE_EPS)
        {
          /* looks like we are safely inside the tetrahedron */

          return 1; /* our point is really nicely inside the tetrahedron */
        }

      if(x[0] < -INSIDE_EPS || x[1] < -INSIDE_EPS || x[2] < -INSIDE_EPS || (1 - (x[0] + x[1] + x[2])) < -INSIDE_EPS)
        {
          /* looks like we are clearly outside the tetrahedron.
             Let's look for a good neighbouring tetrahedron to continue the search */

          /* note: in the (a,b,c) basis, the center-of-mass has coordinates (1/4, 1/4, 1/4) */

          double w, u, v;

          if(ErrorFlag)
            {
              w = 0.25 / (0.25 - x[2]);
              u = 0.25 + w * (x[0] - 0.25);
              v = 0.25 + w * (x[1] - 0.25);
              printf("[3] w=%g u=%g v=%g    fabs(x[2] - 0.25)=%g\n", w, u, v, fabs(x[2] - 0.25));

              w = 0.25 / (0.25 - x[1]);
              u = 0.25 + w * (x[0] - 0.25);
              v = 0.25 + w * (x[2] - 0.25);
              printf("[3] w=%g u=%g v=%g    fabs(x[1] - 0.25)=%g\n", w, u, v, fabs(x[1] - 0.25));

              w = 0.25 / (0.25 - x[0]);
              u = 0.25 + w * (x[1] - 0.25);
              v = 0.25 + w * (x[2] - 0.25);
              printf("[3] w=%g u=%g v=%g    fabs(x[0] - 0.25)=%g\n", w, u, v, fabs(x[0] - 0.25));
            }

          if(fabs(x[2] - 0.25) > INSIDE_EPS)
            {
              w = 0.25 / (0.25 - x[2]);
              if(w > 0)
                {
                  u = 0.25 + w * (x[0] - 0.25);
                  v = 0.25 + w * (x[1] - 0.25);
                  if(u > -INSIDE_EPS && v > -INSIDE_EPS && (1 - (u + v) > -INSIDE_EPS))
                    {
                      *nexttetra = t->t[3];
                      return 0;
                    }
                }
            }

          if(fabs(x[1] - 0.25) > INSIDE_EPS)
            {
              w = 0.25 / (0.25 - x[1]);
              if(w > 0)
                {
                  u = 0.25 + w * (x[0] - 0.25);
                  v = 0.25 + w * (x[2] - 0.25);
                  if(u > -INSIDE_EPS && v > -INSIDE_EPS && (1 - (u + v) > -INSIDE_EPS))
                    {
                      *nexttetra = t->t[2];
                      return 0;
                    }
                }
            }

          if(fabs(x[0] - 0.25) > INSIDE_EPS)
            {
              w = 0.25 / (0.25 - x[0]);
              if(w > 0)
                {
                  u = 0.25 + w * (x[1] - 0.25);
                  v = 0.25 + w * (x[2] - 0.25);
                  if(u > -INSIDE_EPS && v > -INSIDE_EPS && (1 - (u + v) > -INSIDE_EPS))
                    {
                      *nexttetra = t->t[1];
                      return 0;
                    }
                }
            }

          *nexttetra = t->t[0];
          return 0;
        }
    }

  /* here we need to decide whether we have a degenerate case, i.e.
     whether we think the point lies on a face or an edge of the tetrahedron */

  if(ErrorFlag)
    {
      printf("doing exact test for tetra=%d\n", (int)(t - DT));
    }

  Count_InTetraExact++;

  if((ivol = Orient3d_Exact(p0, p1, p2, p3)) <= 0)
    {
      printf("flat or negatively oriented tetrahedron found (vol=%d)\n", ivol);
      terminate("invalid tetrahedron");
    }

  flag3 = Orient3d_Exact(p0, p1, p2, p);
  flag2 = Orient3d_Exact(p0, p3, p1, p);
  flag1 = Orient3d_Exact(p0, p2, p3, p);
  flag0 = Orient3d_Exact(p1, p3, p2, p);

  if(flag0 == 0)
    count_zeros++;

  if(flag1 == 0)
    count_zeros++;

  if(flag2 == 0)
    count_zeros++;

  if(flag3 == 0)
    count_zeros++;

  if(count_zeros > 2)
    {
      printf("task=%d flags=%d %d %d %d  (axb)*c = %g\n", ThisTask, flag0, flag1, flag2, flag3,
             (ay * bz - az * by) * cx + (az * bx - ax * bz) * cy + (ax * by - ay * bx) * cz);

      printf(
          "task=%d pp0=%ld pp1=%ld pp2=%ld pp3=%ld p=%ld IDs=(%llu %llu %llu %llu %llu) pos_0=(%g|%g|%g) pos_1=(%g|%g|%g) "
          "pos_2=(%g|%g|%g) pos_3=(%g|%g|%g) pos=(%g|%g|%g)\n",
          ThisTask, p0 - DP, p1 - DP, p2 - DP, p3 - DP, p - DP, (long long)p0->ID, (long long)p1->ID, (long long)p2->ID,
          (long long)p3->ID, (long long)p->ID, p0->x, p0->y, p0->z, p1->x, p1->y, p1->z, p2->x, p2->y, p2->z, p3->x, p3->y, p3->z,
          p->x, p->y, p->z);

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
      printf("task=%d imageflags=(%d %d %d %d %d)\n", ThisTask, p0->image_flags, p1->image_flags, p2->image_flags, p3->image_flags,
             p->image_flags);
#endif /* #if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z) */
      terminate("strange zero count");
    }

  if(flag0 >= 0 && flag1 >= 0 && flag2 >= 0 && flag3 >= 0)
    {
      /* we have a point inside the tetra, but it may still be on one of the edges */

      if(count_zeros == 0)
        {
          /* ok, let's split the tetra in 4, we are apparently well enough inside */
          return 1;
        }

      if(count_zeros == 1) /* we lie on a face */
        {
          if(flag0 == 0)
            {
              *edgeface_nr = 0;
              return 2;
            }

          if(flag1 == 0)
            {
              *edgeface_nr = 1;
              return 2;
            }

          if(flag2 == 0)
            {
              *edgeface_nr = 2;
              return 2;
            }

          if(flag3 == 0)
            {
              *edgeface_nr = 3;
              return 2;
            }
        }

      if(count_zeros == 2) /* we lie on an edge */
        {
          if(flag0 == 0 && flag1 == 0)
            {
              *edgeface_nr = 5;
              return 3;
            }

          if(flag0 == 0 && flag2 == 0)
            {
              *edgeface_nr = 4;
              return 3;
            }

          if(flag0 == 0 && flag3 == 0)
            {
              *edgeface_nr = 3;
              return 3;
            }

          if(flag1 == 0 && flag2 == 0)
            {
              *edgeface_nr = 2;
              return 3;
            }

          if(flag1 == 0 && flag3 == 0)
            {
              *edgeface_nr = 1;
              return 3;
            }

          if(flag2 == 0 && flag3 == 0)
            {
              *edgeface_nr = 0;
              return 3;
            }
        }
    }

  /* we seem to be lying clearly outside the tetrahedron */
  /* Let's determine a suitable neighbour */

  /* if there is a single negative value, let's pick this side */

  if(flag0 < 0 && flag1 >= 0 && flag2 >= 0 && flag3 >= 0)
    {
      *nexttetra = t->t[0];
      return 0;
    }

  if(flag0 >= 0 && flag1 < 0 && flag2 >= 0 && flag3 >= 0)
    {
      *nexttetra = t->t[1];
      return 0;
    }

  if(flag0 >= 0 && flag1 >= 0 && flag2 < 0 && flag3 >= 0)
    {
      *nexttetra = t->t[2];
      return 0;
    }
  if(flag0 >= 0 && flag1 >= 0 && flag2 >= 0 && flag3 < 0)
    {
      *nexttetra = t->t[3];
      return 0;
    }

  /* there are at least two negative values. Let's pick a random one */

  int ind = -1;

  if(flag0 < 0)
    {
      if(ind < 0)
        ind = 0;
      else
        {
          if(get_random_number() < 0.5)
            ind = 0;
        }
    }

  if(flag1 < 0)
    {
      if(ind < 0)
        ind = 1;
      else
        {
          if(get_random_number() < 0.5)
            ind = 1;
        }
    }

  if(flag2 < 0)
    {
      if(ind < 0)
        ind = 2;
      else
        {
          if(get_random_number() < 0.5)
            ind = 2;
        }
    }

  if(flag3 < 0)
    {
      if(ind < 0)
        ind = 3;
      else
        {
          if(get_random_number() < 0.5)
            ind = 3;
        }
    }

  *nexttetra = t->t[ind];
  return 0;
}

/*! \brief Computes the circum-circle of all tetrahedra in mesh.
 *
 *  \param[in, out] T Pointer to tessellation.
 *
 *  \return void
 */
void compute_circumcircles(tessellation *T)
{
  tetra *DT = T->DT;
  char *DTF = T->DTF;
  int i;

  for(i = 0; i < T->Ndt; i++)
    {
      if(DTF[i] & 1)
        continue;
      DTF[i] |= 1;

      if(DT[i].t[0] < 0) /* deleted ? */
        continue;

      if(DT[i].p[0] == DPinfinity)
        continue;
      if(DT[i].p[1] == DPinfinity)
        continue;
      if(DT[i].p[2] == DPinfinity)
        continue;
      if(DT[i].p[3] == DPinfinity)
        continue;

      update_circumcircle(T, i);
    }
}

/*! \brief Determinant calculation with arbitrary precision arithmetics.
 *
 *  Auxiliary function for exact circum-circle calculation.
 *
 *  \return void
 */
void calc_mpz_determinant(mpz_t det, mpz_t ax, mpz_t ay, mpz_t az, mpz_t bx, mpz_t by, mpz_t bz, mpz_t cx, mpz_t cy, mpz_t cz)
{
  mpz_t bz_cy, by_cz, cz_ay, cy_az, az_by, ay_bz;

  mpz_init(bz_cy);
  mpz_mul(bz_cy, bz, cy);

  mpz_init(by_cz);
  mpz_mul(by_cz, by, cz);

  mpz_init(cz_ay);
  mpz_mul(cz_ay, cz, ay);

  mpz_init(cy_az);
  mpz_mul(cy_az, cy, az);

  mpz_init(az_by);
  mpz_mul(az_by, az, by);

  mpz_init(ay_bz);
  mpz_mul(ay_bz, ay, bz);

  mpz_t bzcy_bycz, czay_cyaz, azby_aybz;

  mpz_init(bzcy_bycz);
  mpz_init(czay_cyaz);
  mpz_init(azby_aybz);

  mpz_sub(bzcy_bycz, bz_cy, by_cz);
  mpz_sub(czay_cyaz, cz_ay, cy_az);
  mpz_sub(azby_aybz, az_by, ay_bz);

  mpz_t a, b, c, ab;

  mpz_init(a);
  mpz_init(b);
  mpz_init(c);

  mpz_mul(a, bzcy_bycz, ax);
  mpz_mul(b, czay_cyaz, bx);
  mpz_mul(c, azby_aybz, cx);

  mpz_init(ab);

  mpz_add(ab, a, b);
  mpz_add(det, ab, c);

  mpz_clear(ab);
  mpz_clear(c);
  mpz_clear(b);
  mpz_clear(a);
  mpz_clear(azby_aybz);
  mpz_clear(czay_cyaz);
  mpz_clear(bzcy_bycz);
  mpz_clear(ay_bz);
  mpz_clear(az_by);
  mpz_clear(cy_az);
  mpz_clear(cz_ay);
  mpz_clear(by_cz);
  mpz_clear(bz_cy);
}

/*! \brief Arbitrary precision calculation of circum-circle.
 *
 *  \param[in, out] T Pointer to tessellation.
 *  \param[in] tt Index in DT array.
 *  \param[out] x X coordinate of circum-circle center.
 *  \param[out] y Y coordinate of circum-circle center.
 *  \param[out] z Z coordinate of circum-circle center.
 *
 *  \return void
 */
void get_circumcircle_exact(tessellation *T, int tt, double *x, double *y, double *z)
{
  tetra *DT = T->DT;
  point *DP = T->DP;
  tetra *t  = &DT[tt];

  point *p0 = &DP[t->p[0]];
  point *p1 = &DP[t->p[1]];
  point *p2 = &DP[t->p[2]];
  point *p3 = &DP[t->p[3]];

  mpz_t det, detA, detB, detC;
  mpz_t qx, qy, qz;
  mpz_t a2, b2, c2, tmp, AA, BB, CC;
  mpz_t ax, ay, az, bx, by, bz, cx, cy, cz;

  mpz_init(det);
  mpz_init(detA);
  mpz_init(detB);
  mpz_init(detC);
  mpz_init(qx);
  mpz_init(qy);
  mpz_init(qz);

  mpz_init(a2);
  mpz_init(b2);
  mpz_init(c2);
  mpz_init(tmp);
  mpz_init(AA);
  mpz_init(BB);
  mpz_init(CC);

  mpz_init(ax);
  mpz_init(ay);
  mpz_init(az);
  mpz_init(bx);
  mpz_init(by);
  mpz_init(bz);
  mpz_init(cx);
  mpz_init(cy);
  mpz_init(cz);

#ifndef OPTIMIZE_MEMORY_USAGE
  MY_mpz_set_si(tmp, p1->ix);
  MY_mpz_sub_ui(ax, tmp, p0->ix);
  MY_mpz_set_si(tmp, p1->iy);
  MY_mpz_sub_ui(ay, tmp, p0->iy);
  MY_mpz_set_si(tmp, p1->iz);
  MY_mpz_sub_ui(az, tmp, p0->iz);

  MY_mpz_set_si(tmp, p2->ix);
  MY_mpz_sub_ui(bx, tmp, p0->ix);
  MY_mpz_set_si(tmp, p2->iy);
  MY_mpz_sub_ui(by, tmp, p0->iy);
  MY_mpz_set_si(tmp, p2->iz);
  MY_mpz_sub_ui(bz, tmp, p0->iz);

  MY_mpz_set_si(tmp, p3->ix);
  MY_mpz_sub_ui(cx, tmp, p0->ix);
  MY_mpz_set_si(tmp, p3->iy);
  MY_mpz_sub_ui(cy, tmp, p0->iy);
  MY_mpz_set_si(tmp, p3->iz);
  MY_mpz_sub_ui(cz, tmp, p0->iz);
#else  /* #ifndef OPTIMIZE_MEMORY_USAGE */
  IntegerMapType pA_ixyz[3], pB_ixyz[3];
  double pA_xyz[3], pB_xyz[3];

  get_integers_for_point(p0, pA_ixyz, pA_xyz);

  get_integers_for_point(p1, pB_ixyz, pB_xyz);
  MY_mpz_set_si(tmp, pB_ixyz[0]);
  MY_mpz_sub_ui(ax, tmp, pA_ixyz[0]);
  MY_mpz_set_si(tmp, pB_ixyz[1]);
  MY_mpz_sub_ui(ay, tmp, pA_ixyz[1]);
  MY_mpz_set_si(tmp, pB_ixyz[2]);
  MY_mpz_sub_ui(az, tmp, pA_ixyz[2]);

  get_integers_for_point(p2, pB_ixyz, pB_xyz);
  MY_mpz_set_si(tmp, pB_ixyz[0]);
  MY_mpz_sub_ui(bx, tmp, pA_ixyz[0]);
  MY_mpz_set_si(tmp, pB_ixyz[1]);
  MY_mpz_sub_ui(by, tmp, pA_ixyz[1]);
  MY_mpz_set_si(tmp, pB_ixyz[2]);
  MY_mpz_sub_ui(bz, tmp, pA_ixyz[2]);

  get_integers_for_point(p3, pB_ixyz, pB_xyz);
  MY_mpz_set_si(tmp, pB_ixyz[0]);
  MY_mpz_sub_ui(cx, tmp, pA_ixyz[0]);
  MY_mpz_set_si(tmp, pB_ixyz[1]);
  MY_mpz_sub_ui(cy, tmp, pA_ixyz[1]);
  MY_mpz_set_si(tmp, pB_ixyz[2]);
  MY_mpz_sub_ui(cz, tmp, pA_ixyz[2]);
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE #else */

  mpz_set(tmp, ax);
  mpz_mul(AA, tmp, ax);
  mpz_set(tmp, ay);
  mpz_mul(BB, tmp, ay);
  mpz_set(tmp, az);
  mpz_mul(CC, tmp, az);
  mpz_add(tmp, AA, BB);
  mpz_add(a2, tmp, CC);

  mpz_set(tmp, bx);
  mpz_mul(AA, tmp, bx);
  mpz_set(tmp, by);
  mpz_mul(BB, tmp, by);
  mpz_set(tmp, bz);
  mpz_mul(CC, tmp, bz);
  mpz_add(tmp, AA, BB);
  mpz_add(b2, tmp, CC);

  mpz_set(tmp, cx);
  mpz_mul(AA, tmp, cx);
  mpz_set(tmp, cy);
  mpz_mul(BB, tmp, cy);
  mpz_set(tmp, cz);
  mpz_mul(CC, tmp, cz);
  mpz_add(tmp, AA, BB);
  mpz_add(c2, tmp, CC);

  calc_mpz_determinant(det, ax, ay, az, bx, by, bz, cx, cy, cz);
  calc_mpz_determinant(detA, a2, ay, az, b2, by, bz, c2, cy, cz);
  calc_mpz_determinant(detB, ax, a2, az, bx, b2, bz, cx, c2, cz);
  calc_mpz_determinant(detC, ax, ay, a2, bx, by, b2, cx, cy, c2);

  mpz_cdiv_q(tmp, detA, det);
  mpz_tdiv_q_2exp(qx, tmp, 1);

  mpz_cdiv_q(tmp, detB, det);
  mpz_tdiv_q_2exp(qy, tmp, 1);

  mpz_cdiv_q(tmp, detC, det);
  mpz_tdiv_q_2exp(qz, tmp, 1);

#ifndef OPTIMIZE_MEMORY_USAGE
  MY_mpz_set_si(tmp, p0->ix);
  mpz_add(AA, qx, tmp);

  MY_mpz_set_si(tmp, p0->iy);
  mpz_add(BB, qy, tmp);

  MY_mpz_set_si(tmp, p0->iz);
  mpz_add(CC, qz, tmp);
#else  /* #ifndef OPTIMIZE_MEMORY_USAGE */
  MY_mpz_set_si(tmp, pA_ixyz[0]);
  mpz_add(AA, qx, tmp);

  MY_mpz_set_si(tmp, pA_ixyz[1]);
  mpz_add(BB, qy, tmp);

  MY_mpz_set_si(tmp, pA_ixyz[2]);
  mpz_add(CC, qz, tmp);
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE #else */
  double xx, yy, zz;

  xx = mpz_get_d(AA);
  yy = mpz_get_d(BB);
  zz = mpz_get_d(CC);

  xx /= (1LLu << USEDBITS);
  yy /= (1LLu << USEDBITS);
  zz /= (1LLu << USEDBITS);

  xx = xx / ConversionFac + CentralOffsetX;
  yy = yy / ConversionFac + CentralOffsetY;
  zz = zz / ConversionFac + CentralOffsetZ;

  *x = xx;
  *y = yy;
  *z = zz;

  mpz_clear(det);
  mpz_clear(detA);
  mpz_clear(detB);
  mpz_clear(detC);
  mpz_clear(qx);
  mpz_clear(qy);
  mpz_clear(qz);

  mpz_clear(a2);
  mpz_clear(b2);
  mpz_clear(c2);
  mpz_clear(tmp);
  mpz_clear(AA);
  mpz_clear(BB);
  mpz_clear(CC);

  mpz_clear(ax);
  mpz_clear(ay);
  mpz_clear(az);
  mpz_clear(bx);
  mpz_clear(by);
  mpz_clear(bz);
  mpz_clear(cx);
  mpz_clear(cy);
  mpz_clear(cz);
}

/*! \brief Computes the circum-circle of tetrahedron tt.
 *
 *  \param[in, out] T Pointer to tessellation.
 *  \param[in] tt Index of triangle in DT array.
 *
 *  \return void
 */
void update_circumcircle(tessellation *T, int tt)
{
  tetra *DT         = T->DT;
  tetra_center *DTC = T->DTC;
  point *DP         = T->DP;
  tetra *t          = &DT[tt];
  tetra_center *tc  = &DTC[tt];

  if(t->t[0] < 0) /* deleted ? */
    return;

  point *p0 = &DP[t->p[0]];
  point *p1 = &DP[t->p[1]];
  point *p2 = &DP[t->p[2]];
  point *p3 = &DP[t->p[3]];

  if(isInfinity(p0) || isInfinity(p1) || isInfinity(p2) || isInfinity(p3))
    return;

#ifndef OPTIMIZE_MEMORY_USAGE
  double ax = p1->xx - p0->xx;
  double ay = p1->yy - p0->yy;
  double az = p1->zz - p0->zz;

  double bx = p2->xx - p0->xx;
  double by = p2->yy - p0->yy;
  double bz = p2->zz - p0->zz;

  double cx = p3->xx - p0->xx;
  double cy = p3->yy - p0->yy;
  double cz = p3->zz - p0->zz;
#else  /* #ifndef OPTIMIZE_MEMORY_USAGE */
  double ax, ay, az, bx, by, bz, cx, cy, cz;
  double pA_xyz[3], pB_xyz[3];
  IntegerMapType pA_ixyz[3], pB_ixyz[3];

  get_integers_for_point(p0, pA_ixyz, pA_xyz);

  get_integers_for_point(p1, pB_ixyz, pB_xyz);
  ax = pB_xyz[0] - pA_xyz[0];
  ay = pB_xyz[1] - pA_xyz[1];
  az = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p2, pB_ixyz, pB_xyz);
  bx = pB_xyz[0] - pA_xyz[0];
  by = pB_xyz[1] - pA_xyz[1];
  bz = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p3, pB_ixyz, pB_xyz);
  cx = pB_xyz[0] - pA_xyz[0];
  cy = pB_xyz[1] - pA_xyz[1];
  cz = pB_xyz[2] - pA_xyz[2];
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE #else */

  double aa = 0.5 * (ax * ax + ay * ay + az * az);
  double bb = 0.5 * (bx * bx + by * by + bz * bz);
  double cc = 0.5 * (cx * cx + cy * cy + cz * cz);

  double mv_data[] = {ax, ay, az, aa, bx, by, bz, bb, cx, cy, cz, cc};
  double x[3];

  int status = solve_linear_equations(mv_data, x);

  if(status < 0)
    {
      if(Orient3d_Exact(p0, p1, p2, p3) != 1)
        {
          printf("p0 = %g %g %g\n", p0->x, p0->y, p0->z);
          printf("p1 = %g %g %g\n", p1->x, p1->y, p1->z);
          printf("p2 = %g %g %g\n", p2->x, p2->y, p2->z);
          printf("p3 = %g %g %g\n", p3->x, p3->y, p3->z);

          printf("Orient-Test=%d\n", Orient3d_Exact(p0, p1, p2, p3));
          printf("tetra-volume=%g  tetra=%d\n", calculate_tetra_volume(p0, p1, p2, p3), tt);

          return;
        }

      double xc, yc, zc;

      get_circumcircle_exact(T, tt, &xc, &yc, &zc);

      tc->cx = xc;
      tc->cy = yc;
      tc->cz = zc;
    }
  else
    {
#ifndef OPTIMIZE_MEMORY_USAGE
      x[0] += p0->xx;
      x[1] += p0->yy;
      x[2] += p0->zz;
#else  /* #ifndef OPTIMIZE_MEMORY_USAGE */
      x[0] += pA_xyz[0];
      x[1] += pA_xyz[1];
      x[2] += pA_xyz[2];
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE #else */

      tc->cx = (x[0] - 1.0) / ConversionFac + CentralOffsetX;
      tc->cy = (x[1] - 1.0) / ConversionFac + CentralOffsetY;
      tc->cz = (x[2] - 1.0) / ConversionFac + CentralOffsetZ;
    }
}

/*! \brief Returns the orientation of the tetrahedron.
 *
 *  \param[in] p0 Point spanning the tetrahedron.
 *  \param[in] p1 Point spanning the tetrahedron.
 *  \param[in] p2 Point spanning the tetrahedron.
 *  \param[in] p3 Point spanning the tetrahedron.
 *
 *  \return -1: negative orientation; +1 positive orientation.
 */
int test_tetra_orientation(point *p0, point *p1, point *p2, point *p3)
{
  double nx, ny, nz;

  if(isInfinity(p0) || isInfinity(p1) || isInfinity(p2) || isInfinity(p3))
    return +1;

#ifndef OPTIMIZE_MEMORY_USAGE
  nx = (p1->yy - p0->yy) * (p2->zz - p0->zz) - (p1->zz - p0->zz) * (p2->yy - p0->yy);
  ny = (p1->zz - p0->zz) * (p2->xx - p0->xx) - (p1->xx - p0->xx) * (p2->zz - p0->zz);
  nz = (p1->xx - p0->xx) * (p2->yy - p0->yy) - (p1->yy - p0->yy) * (p2->xx - p0->xx);
  if(nx * (p3->xx - p0->xx) + ny * (p3->yy - p0->yy) + nz * (p3->zz - p0->zz) >= 0)
    return +1;
  else
    return -1;
#else  /* #ifndef OPTIMIZE_MEMORY_USAGE */
  IntegerMapType p0_ixyz[3], p1_ixyz[3], p2_ixyz[3], p3_ixyz[3];
  double p0_xyz[3], p1_xyz[3], p2_xyz[3], p3_xyz[3];

  get_integers_for_point(p0, p0_ixyz, p0_xyz);
  get_integers_for_point(p1, p1_ixyz, p1_xyz);
  get_integers_for_point(p2, p2_ixyz, p2_xyz);
  get_integers_for_point(p3, p3_ixyz, p3_xyz);

  nx = (p1_xyz[1] - p0_xyz[1]) * (p2_xyz[2] - p0_xyz[2]) - (p1_xyz[2] - p0_xyz[2]) * (p2_xyz[1] - p0_xyz[1]);
  ny = (p1_xyz[2] - p0_xyz[2]) * (p2_xyz[0] - p0_xyz[0]) - (p1_xyz[0] - p0_xyz[0]) * (p2_xyz[2] - p0_xyz[2]);
  nz = (p1_xyz[0] - p0_xyz[0]) * (p2_xyz[1] - p0_xyz[1]) - (p1_xyz[1] - p0_xyz[1]) * (p2_xyz[0] - p0_xyz[0]);

  get_integers_for_point(p3, p3_ixyz, p3_xyz);

  if(nx * (p3_xyz[0] - p0_xyz[0]) + ny * (p3_xyz[1] - p0_xyz[1]) + nz * (p3_xyz[2] - p0_xyz[2]) >= 0)
    return +1;
  else
    return -1;
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE #else */
}

/*! \brief Calculate the volume of a tetrahedron.
 *
 *  \param[in] p0 Point spanning the tetrahedron.
 *  \param[in] p1 Point spanning the tetrahedron.
 *  \param[in] p2 Point spanning the tetrahedron.
 *  \param[in] p3 Point spanning the tetrahedron.
 *
 *  \return Volume of the tetrahedron.
 */
double calculate_tetra_volume(point *p0, point *p1, point *p2, point *p3)
{
  double nx, ny, nz;

  if(isInfinity(p0) || isInfinity(p1) || isInfinity(p2) || isInfinity(p3))
    return +1;

  nx = (p1->y - p0->y) * (p2->z - p0->z) - (p1->z - p0->z) * (p2->y - p0->y);
  ny = (p1->z - p0->z) * (p2->x - p0->x) - (p1->x - p0->x) * (p2->z - p0->z);
  nz = (p1->x - p0->x) * (p2->y - p0->y) - (p1->y - p0->y) * (p2->x - p0->x);

  return nx * (p3->x - p0->x) + ny * (p3->y - p0->y) + nz * (p3->z - p0->z);
}

/*! \brief Add row in matrix equation.
 *
 *  Auxiliary function for solve_linear_equations.
 *
 *  \param[in, out] m Matrix.
 *  \param[in] r1 Index of row to be modified.
 *  \param[in] r2 Index of row which is added to r1.
 *  \param[in] fac Factor by which row r2 is multiplied before adding to r1.
 *
 *  \return void
 */
void add_row(double *m, int r1, int r2, double fac)
{
  int i;

  for(i = 0; i < 4; i++)
    m[r1 * 4 + i] += fac * m[r2 * 4 + i];
}

/*! \brief Solve system of linear equations for 3d Voronoi construction.
 *
 *  \param[in, out] m Matrix.
 *  \param[out] res Result.
 *
 *  \return 0 if success, <0 else.
 */
int solve_linear_equations(double *m, double *res)
{
  int ix, iy, iz, itmp;

  if(fabs(m[4]) > fabs(m[0]))
    {
      ix = 1;
      iy = 0;
      iz = 2;
    }
  else
    {
      ix = 0;
      iy = 1;
      iz = 2;
    }

  if(fabs(m[8]) > fabs(m[ix * 4]))
    {
      ix = 2;
      iy = 0;
      iz = 1;
    }

  add_row(m, iy, ix, -m[iy * 4] / m[ix * 4]);
  add_row(m, iz, ix, -m[iz * 4] / m[ix * 4]);

  if(fabs(m[iz * 4 + 1]) > fabs(m[iy * 4 + 1]))
    {
      /* swap iy/iz */
      itmp = iy;
      iy   = iz;
      iz   = itmp;
    }

  if(fabs(m[iy * 4 + 1]) < GAUSS_EPS)
    return -1;

  add_row(m, iz, iy, -m[iz * 4 + 1] / m[iy * 4 + 1]);

  res[2] = m[iz * 4 + 3] / m[iz * 4 + 2];
  res[1] = (m[iy * 4 + 3] - res[2] * m[iy * 4 + 2]) / m[iy * 4 + 1];
  res[0] = (m[ix * 4 + 3] - res[2] * m[ix * 4 + 2] - res[1] * m[ix * 4 + 1]) / m[ix * 4];

  if(fabs(m[iz * 4 + 2]) < GAUSS_EPS)
    {
      return -1;
    }
  if(fabs(m[iy * 4 + 1]) < GAUSS_EPS)
    {
      return -2;
    }
  if(fabs(m[ix * 4]) < GAUSS_EPS)
    {
      return -3;
    }

  return 0;
}

/*! \brief Converts coordinates of point p to integer values.
 *
 *  \param[in, out] p Point.
 *
 *  \return void
 */
#ifndef OPTIMIZE_MEMORY_USAGE
void set_integers_for_pointer(point *p)
{
  p->xx = (p->x - CentralOffsetX) * ConversionFac + 1.0;
  p->yy = (p->y - CentralOffsetY) * ConversionFac + 1.0;
  p->zz = (p->z - CentralOffsetZ) * ConversionFac + 1.0;

  if(p->xx < 1.0 || p->xx >= 2.0 || p->yy < 1.0 || p->yy >= 2.0 || p->zz < 1.0 || p->zz >= 2.0)
    {
      printf("(%g, %g, %g) (%g, %g, %g)\n", p->x, p->y, p->z, p->xx, p->yy, p->zz);
      terminate("invalid coordinate range");
    }

  p->ix = double_to_voronoiint(p->xx);
  p->iy = double_to_voronoiint(p->yy);
  p->iz = double_to_voronoiint(p->zz);

  p->xx = mask_voronoi_int(p->xx);
  p->yy = mask_voronoi_int(p->yy);
  p->zz = mask_voronoi_int(p->zz);
}
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE */

/*! \brief Checks if point is within a sphere using arbitrary precision
 *         operations.
 *
 *  \param p0 Point 1 of tetrahedron.
 *  \param p1 Point 2 of tetrahedron.
 *  \param p2 Point 3 of tetrahedron.
 *  \param p3 Point 4 of tetrahedron.
 *  \param p Point to be checked if it is in cricumsphere.
 *
 *  \return (-1,1); -1 in sphere, 1 outside.
 */
int InSphere_Exact(point *p0, point *p1, point *p2, point *p3, point *p)
{
  IntegerMapType ax, bx, cx, dx;
  IntegerMapType ay, by, cy, dy;
  IntegerMapType az, bz, cz, dz;

  if(isInfinity(p0) || isInfinity(p1) || isInfinity(p2) || isInfinity(p3))
    return -1;

#ifndef OPTIMIZE_MEMORY_USAGE
  ax = p0->ix - p->ix;
  ay = p0->iy - p->iy;
  az = p0->iz - p->iz;

  bx = p1->ix - p->ix;
  by = p1->iy - p->iy;
  bz = p1->iz - p->iz;

  cx = p2->ix - p->ix;
  cy = p2->iy - p->iy;
  cz = p2->iz - p->iz;

  dx = p3->ix - p->ix;
  dy = p3->iy - p->iy;
  dz = p3->iz - p->iz;
#else  /* #ifndef OPTIMIZE_MEMORY_USAGE */
  double pA_xyz[3], pB_xyz[3];
  IntegerMapType pA_ixyz[3], pB_ixyz[3];

  get_integers_for_point(p, pA_ixyz, pA_xyz);

  get_integers_for_point(p0, pB_ixyz, pB_xyz);
  ax = pB_ixyz[0] - pA_ixyz[0];
  ay = pB_ixyz[1] - pA_ixyz[1];
  az = pB_ixyz[2] - pA_ixyz[2];

  get_integers_for_point(p1, pB_ixyz, pB_xyz);
  bx = pB_ixyz[0] - pA_ixyz[0];
  by = pB_ixyz[1] - pA_ixyz[1];
  bz = pB_ixyz[2] - pA_ixyz[2];

  get_integers_for_point(p2, pB_ixyz, pB_xyz);
  cx = pB_ixyz[0] - pA_ixyz[0];
  cy = pB_ixyz[1] - pA_ixyz[1];
  cz = pB_ixyz[2] - pA_ixyz[2];

  get_integers_for_point(p3, pB_ixyz, pB_xyz);
  dx = pB_ixyz[0] - pA_ixyz[0];
  dy = pB_ixyz[1] - pA_ixyz[1];
  dz = pB_ixyz[2] - pA_ixyz[2];
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE #else */

  mpz_t ab, bc, cd, da, ac, bd;

  mpz_init(ab);
  mpz_init(bc);
  mpz_init(cd);
  mpz_init(da);
  mpz_init(ac);
  mpz_init(bd);

  mpz_t tmp, AA, BB, CC;

  mpz_init(tmp);
  mpz_init(AA);
  mpz_init(BB);
  mpz_init(CC);

  MY_mpz_set_si(tmp, ax);
  MY_mpz_mul_si(AA, tmp, by);
  MY_mpz_set_si(tmp, bx);
  MY_mpz_mul_si(BB, tmp, ay);
  mpz_sub(ab, AA, BB);

  MY_mpz_set_si(tmp, bx);
  MY_mpz_mul_si(AA, tmp, cy);
  MY_mpz_set_si(tmp, cx);
  MY_mpz_mul_si(BB, tmp, by);
  mpz_sub(bc, AA, BB);

  MY_mpz_set_si(tmp, cx);
  MY_mpz_mul_si(AA, tmp, dy);
  MY_mpz_set_si(tmp, dx);
  MY_mpz_mul_si(BB, tmp, cy);
  mpz_sub(cd, AA, BB);

  MY_mpz_set_si(tmp, dx);
  MY_mpz_mul_si(AA, tmp, ay);
  MY_mpz_set_si(tmp, ax);
  MY_mpz_mul_si(BB, tmp, dy);
  mpz_sub(da, AA, BB);

  MY_mpz_set_si(tmp, ax);
  MY_mpz_mul_si(AA, tmp, cy);
  MY_mpz_set_si(tmp, cx);
  MY_mpz_mul_si(BB, tmp, ay);
  mpz_sub(ac, AA, BB);

  MY_mpz_set_si(tmp, bx);
  MY_mpz_mul_si(AA, tmp, dy);
  MY_mpz_set_si(tmp, dx);
  MY_mpz_mul_si(BB, tmp, by);
  mpz_sub(bd, AA, BB);

  mpz_t abc, bcd, cda, dab;

  mpz_init(abc);
  mpz_init(bcd);
  mpz_init(cda);
  mpz_init(dab);

  MY_mpz_mul_si(AA, bc, az);
  MY_mpz_mul_si(BB, ac, -bz);
  MY_mpz_mul_si(CC, ab, cz);
  mpz_add(tmp, AA, BB);
  mpz_add(abc, tmp, CC);

  MY_mpz_mul_si(AA, cd, bz);
  MY_mpz_mul_si(BB, bd, -cz);
  MY_mpz_mul_si(CC, bc, dz);
  mpz_add(tmp, AA, BB);
  mpz_add(bcd, tmp, CC);

  MY_mpz_mul_si(AA, da, cz);
  MY_mpz_mul_si(BB, ac, dz);
  MY_mpz_mul_si(CC, cd, az);
  mpz_add(tmp, AA, BB);
  mpz_add(cda, tmp, CC);

  MY_mpz_mul_si(AA, ab, dz);
  MY_mpz_mul_si(BB, bd, az);
  MY_mpz_mul_si(CC, da, bz);
  mpz_add(tmp, AA, BB);
  mpz_add(dab, tmp, CC);

  mpz_t a2, b2, c2, d2;

  mpz_init(a2);
  mpz_init(b2);
  mpz_init(c2);
  mpz_init(d2);

  MY_mpz_set_si(tmp, ax);
  MY_mpz_mul_si(AA, tmp, ax);
  MY_mpz_set_si(tmp, ay);
  MY_mpz_mul_si(BB, tmp, ay);
  MY_mpz_set_si(tmp, az);
  MY_mpz_mul_si(CC, tmp, az);
  mpz_add(tmp, AA, BB);
  mpz_add(a2, tmp, CC);

  MY_mpz_set_si(tmp, bx);
  MY_mpz_mul_si(AA, tmp, bx);
  MY_mpz_set_si(tmp, by);
  MY_mpz_mul_si(BB, tmp, by);
  MY_mpz_set_si(tmp, bz);
  MY_mpz_mul_si(CC, tmp, bz);
  mpz_add(tmp, AA, BB);
  mpz_add(b2, tmp, CC);

  MY_mpz_set_si(tmp, cx);
  MY_mpz_mul_si(AA, tmp, cx);
  MY_mpz_set_si(tmp, cy);
  MY_mpz_mul_si(BB, tmp, cy);
  MY_mpz_set_si(tmp, cz);
  MY_mpz_mul_si(CC, tmp, cz);
  mpz_add(tmp, AA, BB);
  mpz_add(c2, tmp, CC);

  MY_mpz_set_si(tmp, dx);
  MY_mpz_mul_si(AA, tmp, dx);
  MY_mpz_set_si(tmp, dy);
  MY_mpz_mul_si(BB, tmp, dy);
  MY_mpz_set_si(tmp, dz);
  MY_mpz_mul_si(CC, tmp, dz);
  mpz_add(tmp, AA, BB);
  mpz_add(d2, tmp, CC);

  /* now calculate final result */

  mpz_mul(AA, c2, dab);
  mpz_mul(BB, d2, abc);
  mpz_sub(tmp, AA, BB);

  mpz_mul(AA, a2, bcd);
  mpz_mul(BB, b2, cda);
  mpz_sub(CC, AA, BB);

  mpz_add(AA, tmp, CC);

  /* AA now contains the result */

  int sign = mpz_sgn(AA);

  mpz_clear(d2);
  mpz_clear(c2);
  mpz_clear(b2);
  mpz_clear(a2);
  mpz_clear(dab);
  mpz_clear(cda);
  mpz_clear(bcd);
  mpz_clear(abc);
  mpz_clear(CC);
  mpz_clear(BB);
  mpz_clear(AA);
  mpz_clear(tmp);
  mpz_clear(bd);
  mpz_clear(ac);
  mpz_clear(da);
  mpz_clear(cd);
  mpz_clear(bc);
  mpz_clear(ab);

  return sign;
}

/*! \brief Checks if point is within a sphere.
 *
 *  \param p0 Point 1 of tetrahedron.
 *  \param p1 Point 2 of tetrahedron.
 *  \param p2 Point 3 of tetrahedron.
 *  \param p3 Point 4 of tetrahedron.
 *  \param p Point to be checked if it is in cricumsphere.
 *
 *  \return (-1,0,1); -1: in sphere, 0: on surfrace, 1: outside.
 */
int InSphere_Quick(point *p0, point *p1, point *p2, point *p3, point *p)
{
  double ax, bx, cx, dx;
  double ay, by, cy, dy;
  double az, bz, cz, dz;
  double a2, b2, c2, d2;
  double ab, bc, cd, da, ac, bd;
  double abc, bcd, cda, dab;
  double x;

  if(isInfinity(p0) || isInfinity(p1) || isInfinity(p2) || isInfinity(p3))
    return -1;

#ifndef OPTIMIZE_MEMORY_USAGE
  ax = p0->xx - p->xx;
  ay = p0->yy - p->yy;
  az = p0->zz - p->zz;

  bx = p1->xx - p->xx;
  by = p1->yy - p->yy;
  bz = p1->zz - p->zz;

  cx = p2->xx - p->xx;
  cy = p2->yy - p->yy;
  cz = p2->zz - p->zz;

  dx = p3->xx - p->xx;
  dy = p3->yy - p->yy;
  dz = p3->zz - p->zz;
#else  /* #ifndef OPTIMIZE_MEMORY_USAGE */
  double pA_xyz[3], pB_xyz[3];
  IntegerMapType pA_ixyz[3], pB_ixyz[3];

  get_integers_for_point(p, pA_ixyz, pA_xyz);

  get_integers_for_point(p0, pB_ixyz, pB_xyz);
  ax = pB_xyz[0] - pA_xyz[0];
  ay = pB_xyz[1] - pA_xyz[1];
  az = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p1, pB_ixyz, pB_xyz);
  bx = pB_xyz[0] - pA_xyz[0];
  by = pB_xyz[1] - pA_xyz[1];
  bz = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p2, pB_ixyz, pB_xyz);
  cx = pB_xyz[0] - pA_xyz[0];
  cy = pB_xyz[1] - pA_xyz[1];
  cz = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p3, pB_ixyz, pB_xyz);
  dx = pB_xyz[0] - pA_xyz[0];
  dy = pB_xyz[1] - pA_xyz[1];
  dz = pB_xyz[2] - pA_xyz[2];
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE #else */

  ab = ax * by - bx * ay;
  bc = bx * cy - cx * by;
  cd = cx * dy - dx * cy;
  da = dx * ay - ax * dy;
  ac = ax * cy - cx * ay;
  bd = bx * dy - dx * by;

  abc = az * bc - bz * ac + cz * ab;
  bcd = bz * cd - cz * bd + dz * bc;
  cda = cz * da + dz * ac + az * cd;
  dab = dz * ab + az * bd + bz * da;

  a2 = ax * ax + ay * ay + az * az;
  b2 = bx * bx + by * by + bz * bz;
  c2 = cx * cx + cy * cy + cz * cz;
  d2 = dx * dx + dy * dy + dz * dz;

  x = ((c2 * dab - d2 * abc) + (a2 * bcd - b2 * cda));

  if(x < 0)
    return -1;
  if(x > 0)
    return +1;

  return 0;
}

/*! \brief Checks if point is within a sphere with some error margin.
 *
 *  \param p0 Point 1 of tetrahedron.
 *  \param p1 Point 2 of tetrahedron.
 *  \param p2 Point 3 of tetrahedron.
 *  \param p3 Point 4 of tetrahedron.
 *  \param p Point to be checked if it is in cricumsphere.
 *
 *  \return (-1,0,1); -1: in sphere, 0: on surfrace (within error margin),
 *                    +1: outside.
 */
int InSphere_Errorbound(point *p0, point *p1, point *p2, point *p3, point *p)
{
  double ax, bx, cx, dx;
  double ay, by, cy, dy;
  double az, bz, cz, dz;
  double a2, b2, c2, d2;
  double ab, bc, cd, da, ac, bd;
  double abc, bcd, cda, dab;
  double x;

  if(isInfinity(p0) || isInfinity(p1) || isInfinity(p2) || isInfinity(p3))
    return -1;

#ifndef OPTIMIZE_MEMORY_USAGE
  ax = p0->xx - p->xx;
  ay = p0->yy - p->yy;
  az = p0->zz - p->zz;

  bx = p1->xx - p->xx;
  by = p1->yy - p->yy;
  bz = p1->zz - p->zz;

  cx = p2->xx - p->xx;
  cy = p2->yy - p->yy;
  cz = p2->zz - p->zz;

  dx = p3->xx - p->xx;
  dy = p3->yy - p->yy;
  dz = p3->zz - p->zz;
#else  /* #ifndef OPTIMIZE_MEMORY_USAGE */
  double pA_xyz[3], pB_xyz[3];
  IntegerMapType pA_ixyz[3], pB_ixyz[3];

  get_integers_for_point(p, pA_ixyz, pA_xyz);

  get_integers_for_point(p0, pB_ixyz, pB_xyz);
  ax = pB_xyz[0] - pA_xyz[0];
  ay = pB_xyz[1] - pA_xyz[1];
  az = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p1, pB_ixyz, pB_xyz);
  bx = pB_xyz[0] - pA_xyz[0];
  by = pB_xyz[1] - pA_xyz[1];
  bz = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p2, pB_ixyz, pB_xyz);
  cx = pB_xyz[0] - pA_xyz[0];
  cy = pB_xyz[1] - pA_xyz[1];
  cz = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p3, pB_ixyz, pB_xyz);
  dx = pB_xyz[0] - pA_xyz[0];
  dy = pB_xyz[1] - pA_xyz[1];
  dz = pB_xyz[2] - pA_xyz[2];
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE #else */

  double axby = ax * by;
  double bxay = bx * ay;
  double bxcy = bx * cy;
  double cxby = cx * by;
  double cxdy = cx * dy;
  double dxcy = dx * cy;
  double dxay = dx * ay;
  double axdy = ax * dy;
  double axcy = ax * cy;
  double cxay = cx * ay;
  double bxdy = bx * dy;
  double dxby = dx * by;

  ab = axby - bxay;
  bc = bxcy - cxby;
  cd = cxdy - dxcy;
  da = dxay - axdy;
  ac = axcy - cxay;
  bd = bxdy - dxby;

  abc = az * bc - bz * ac + cz * ab;
  bcd = bz * cd - cz * bd + dz * bc;
  cda = cz * da + dz * ac + az * cd;
  dab = dz * ab + az * bd + bz * da;

  a2 = ax * ax + ay * ay + az * az;
  b2 = bx * bx + by * by + bz * bz;
  c2 = cx * cx + cy * cy + cz * cz;
  d2 = dx * dx + dy * dy + dz * dz;

  x = ((c2 * dab - d2 * abc) + (a2 * bcd - b2 * cda));

  /* calculate absolute maximum size */

  ab = fabs(axby) + fabs(bxay);
  bc = fabs(bxcy) + fabs(cxby);
  cd = fabs(cxdy) + fabs(dxcy);
  da = fabs(dxay) + fabs(axdy);
  ac = fabs(axcy) + fabs(cxay);
  bd = fabs(bxdy) + fabs(dxby);

  az = fabs(az);
  bz = fabs(bz);
  cz = fabs(cz);
  dz = fabs(dz);

  abc = az * bc + bz * ac + cz * ab;
  bcd = bz * cd + cz * bd + dz * bc;
  cda = cz * da + dz * ac + az * cd;
  dab = dz * ab + az * bd + bz * da;

  double sizelimit = ((c2 * dab + d2 * abc) + (a2 * bcd + b2 * cda));

  double errbound = 1.0e-14 * sizelimit;

  if(x < -errbound)
    return -1;
  else if(x > errbound)
    return +1;

  return 0;
}

/*! \brief Returns orientation of tetrahedron using arbitrary precision
 *         floating point operations.
 *
 *  \param[in] p0 First point of tetrahedron.
 *  \param[in] p1 Second point of tetrahedron.
 *  \param[in] p2 Third point of tetrahedron.
 *  \param[in] p3 Forth point of tetrahedron.
 *
 *  \return (-1,0,1) -1 if negatively oriented, 0 if degenerate and 1 if
 *                   positively oriented.
 */
int Orient3d_Exact(point *p0, point *p1, point *p2, point *p3)
{
  IntegerMapType ax, bx, cx;
  IntegerMapType ay, by, cy;
  IntegerMapType az, bz, cz;

#ifndef OPTIMIZE_MEMORY_USAGE
  ax = p0->ix - p3->ix;
  ay = p0->iy - p3->iy;
  az = p0->iz - p3->iz;

  bx = p1->ix - p3->ix;
  by = p1->iy - p3->iy;
  bz = p1->iz - p3->iz;

  cx = p2->ix - p3->ix;
  cy = p2->iy - p3->iy;
  cz = p2->iz - p3->iz;
#else  /* #ifndef OPTIMIZE_MEMORY_USAGE */
  double pA_xyz[3], pB_xyz[3];
  IntegerMapType pA_ixyz[3], pB_ixyz[3];

  get_integers_for_point(p3, pA_ixyz, pA_xyz);

  get_integers_for_point(p0, pB_ixyz, pB_xyz);
  ax = pB_ixyz[0] - pA_ixyz[0];
  ay = pB_ixyz[1] - pA_ixyz[1];
  az = pB_ixyz[2] - pA_ixyz[2];

  get_integers_for_point(p1, pB_ixyz, pB_xyz);
  bx = pB_ixyz[0] - pA_ixyz[0];
  by = pB_ixyz[1] - pA_ixyz[1];
  bz = pB_ixyz[2] - pA_ixyz[2];

  get_integers_for_point(p2, pB_ixyz, pB_xyz);
  cx = pB_ixyz[0] - pA_ixyz[0];
  cy = pB_ixyz[1] - pA_ixyz[1];
  cz = pB_ixyz[2] - pA_ixyz[2];
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE */

  mpz_t bz_cy, by_cz, cz_ay, cy_az, az_by, ay_bz;
  mpz_t bz2, by2, cz2, cy2, az2, ay2;

  mpz_init(bz_cy);
  mpz_init(bz2);
  MY_mpz_set_si(bz2, bz);
  MY_mpz_mul_si(bz_cy, bz2, cy);

  mpz_init(by_cz);
  mpz_init(by2);
  MY_mpz_set_si(by2, by);
  MY_mpz_mul_si(by_cz, by2, cz);

  mpz_init(cz_ay);
  mpz_init(cz2);
  MY_mpz_set_si(cz2, cz);
  MY_mpz_mul_si(cz_ay, cz2, ay);

  mpz_init(cy_az);
  mpz_init(cy2);
  MY_mpz_set_si(cy2, cy);
  MY_mpz_mul_si(cy_az, cy2, az);

  mpz_init(az_by);
  mpz_init(az2);
  MY_mpz_set_si(az2, az);
  MY_mpz_mul_si(az_by, az2, by);

  mpz_init(ay_bz);
  mpz_init(ay2);
  MY_mpz_set_si(ay2, ay);
  MY_mpz_mul_si(ay_bz, ay2, bz);

  mpz_t bzcy_bycz, czay_cyaz, azby_aybz;

  mpz_init(bzcy_bycz);
  mpz_init(czay_cyaz);
  mpz_init(azby_aybz);

  mpz_sub(bzcy_bycz, bz_cy, by_cz);
  mpz_sub(czay_cyaz, cz_ay, cy_az);
  mpz_sub(azby_aybz, az_by, ay_bz);

  mpz_t a, b, c, ab, res;

  mpz_init(a);
  mpz_init(b);
  mpz_init(c);

  MY_mpz_mul_si(a, bzcy_bycz, ax);
  MY_mpz_mul_si(b, czay_cyaz, bx);
  MY_mpz_mul_si(c, azby_aybz, cx);

  mpz_init(ab);
  mpz_init(res);

  mpz_add(ab, a, b);
  mpz_add(res, ab, c);

  int sign = mpz_sgn(res);

  mpz_clear(res);
  mpz_clear(ab);
  mpz_clear(c);
  mpz_clear(b);
  mpz_clear(a);
  mpz_clear(azby_aybz);
  mpz_clear(czay_cyaz);
  mpz_clear(bzcy_bycz);
  mpz_clear(ay2);
  mpz_clear(ay_bz);
  mpz_clear(az2);
  mpz_clear(az_by);
  mpz_clear(cy2);
  mpz_clear(cy_az);
  mpz_clear(cz2);
  mpz_clear(cz_ay);
  mpz_clear(by2);
  mpz_clear(by_cz);
  mpz_clear(bz2);
  mpz_clear(bz_cy);

  return sign;
}

/*! \brief Returns orientation of tetrahedron.
 *
 *  \param[in] p0 First point of tetrahedron.
 *  \param[in] p1 Second point of tetrahedron.
 *  \param[in] p2 Third point of tetrahedron.
 *  \param[in] p3 Forth point of tetrahedron.
 *
 *  \return (-1,0,1) -1 if negatively oriented, 0 if degenerate and 1 if
 *                   positively oriented.
 */
int Orient3d_Quick(point *p0, point *p1, point *p2, point *p3)
{
  double ax, bx, cx;
  double ay, by, cy;
  double az, bz, cz;

#ifndef OPTIMIZE_MEMORY_USAGE
  ax = p0->xx - p3->xx;
  ay = p0->yy - p3->yy;
  az = p0->zz - p3->zz;

  bx = p1->xx - p3->xx;
  by = p1->yy - p3->yy;
  bz = p1->zz - p3->zz;

  cx = p2->xx - p3->xx;
  cy = p2->yy - p3->yy;
  cz = p2->zz - p3->zz;
#else  /* #ifndef OPTIMIZE_MEMORY_USAGE */
  double pA_xyz[3], pB_xyz[3];
  IntegerMapType pA_ixyz[3], pB_ixyz[3];

  get_integers_for_point(p3, pA_ixyz, pA_xyz);

  get_integers_for_point(p0, pB_ixyz, pB_xyz);
  ax = pB_xyz[0] - pA_xyz[0];
  ay = pB_xyz[1] - pA_xyz[1];
  az = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p1, pB_ixyz, pB_xyz);
  bx = pB_xyz[0] - pA_xyz[0];
  by = pB_xyz[1] - pA_xyz[1];
  bz = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p2, pB_ixyz, pB_xyz);
  cx = pB_xyz[0] - pA_xyz[0];
  cy = pB_xyz[1] - pA_xyz[1];
  cz = pB_xyz[2] - pA_xyz[2];
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE #else */

  double x = (ax * (bz * cy - by * cz) + bx * (cz * ay - cy * az) + cx * (az * by - ay * bz));

  if(x < 0)
    return -1;
  else if(x > 0)
    return +1;

  return 0;
}

/* \brief Returns orientation of tetrahedron.
 *
 *  \param[in] p0 First point of tetrahedron.
 *  \param[in] p1 Second point of tetrahedron.
 *  \param[in] p2 Third point of tetrahedron.
 *  \param[in] p3 Forth point of tetrahedron.
 *
 *  \return (-1,0,1) the orientation of the 4 points as +/-1. If either of the
 *          points is an infinity point, return 0.
 */
int Orient3d(point *p0, point *p1, point *p2, point *p3)
{
  if(isInfinity(p0) || isInfinity(p1) || isInfinity(p2) || isInfinity(p3))
    return 0;

#ifndef OPTIMIZE_MEMORY_USAGE
  double ax = p0->xx - p3->xx;
  double ay = p0->yy - p3->yy;
  double az = p0->zz - p3->zz;

  double bx = p1->xx - p3->xx;
  double by = p1->yy - p3->yy;
  double bz = p1->zz - p3->zz;

  double cx = p2->xx - p3->xx;
  double cy = p2->yy - p3->yy;
  double cz = p2->zz - p3->zz;
#else  /* #ifndef OPTIMIZE_MEMORY_USAGE */
  double ax, ay, az, bx, by, bz, cx, cy, cz;
  double pA_xyz[3], pB_xyz[3];
  IntegerMapType pA_ixyz[3], pB_ixyz[3];

  get_integers_for_point(p3, pA_ixyz, pA_xyz);

  get_integers_for_point(p0, pB_ixyz, pB_xyz);
  ax = pB_xyz[0] - pA_xyz[0];
  ay = pB_xyz[1] - pA_xyz[1];
  az = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p1, pB_ixyz, pB_xyz);
  bx = pB_xyz[0] - pA_xyz[0];
  by = pB_xyz[1] - pA_xyz[1];
  bz = pB_xyz[2] - pA_xyz[2];

  get_integers_for_point(p2, pB_ixyz, pB_xyz);
  cx = pB_xyz[0] - pA_xyz[0];
  cy = pB_xyz[1] - pA_xyz[1];
  cz = pB_xyz[2] - pA_xyz[2];
#endif /* #ifndef OPTIMIZE_MEMORY_USAGE #else */

  double bzcy = bz * cy;
  double bycz = by * cz;
  double czay = cz * ay;
  double cyaz = cy * az;
  double azby = az * by;
  double aybz = ay * bz;

  double x = ax * (bzcy - bycz) + bx * (czay - cyaz) + cx * (azby - aybz);

  double sizelimit =
      fabs(ax) * (fabs(bzcy) + fabs(bycz)) + fabs(bx) * (fabs(czay) + fabs(cyaz)) + fabs(cx) * (fabs(azby) + fabs(aybz));

  double errbound = 1.0e-14 * sizelimit;

  if(x < -errbound)
    return -1;
  else if(x > errbound)
    return +1;

  return Orient3d_Exact(p0, p1, p2, p3);
}

/*! \brief Data structure for face sort
 */
struct data_face_sort /* for sorting faces */
{
  MyIDType ID;     /* ID of corresponding cell */
  float normal[3]; /* non-normalized normal vector */
  int start;       /* start index into vertex list */
  int len;         /* number of vertices */
};

static int *VertexEntries;       /* face index list */
static float *VertexCoordinates; /* Voronoi vertex coordinates (circumsphere centers of delaunay tetras) */
static float *FaceNormals;       /* normal vectors */
static int Nvertices;            /* number of Voronoi vertices */
static int Nnormals;             /* number of normals */
static int Nentries;             /* number of entries in Voronoi face vertex list (including IDs and face vertex count) */
static int Nsort;                /* number of ID sorted faces */
static int MaxEntries, MaxFaces; /* for allocation */
static struct data_face_sort *FaceSort;

/*! \brief  Face sorting kernel
 *
 *  Compares ID of data_face_sort types.
 *
 *  \param[in] a Fist element.
 *  \param[in] b Second element.
 *
 *  \return (-1,0,1), -1 if a->ID < b ->ID.
 */
int compare_face_sort(const void *a, const void *b)
{
  if(((struct data_face_sort *)a)->ID < ((struct data_face_sort *)b)->ID)
    return -1;

  if(((struct data_face_sort *)a)->ID > ((struct data_face_sort *)b)->ID)
    return +1;

  return 0;
}

/*! \brief Gathers faces in list.
 *
 *  \param[in] T Pointer to tessellation.
 *
 *  \return void
 */
void get_voronoi_face_vertex_indices(tessellation *T)
{
  int i, j, k, l, m, ii, jj, kk, ll, tetra_nr, edge_nr, next_tetra_nr, count, dp_1, dp_2;
  tetra *prev, *next;
  tetra *DT = T->DT;
  point *DP = T->DP;
  int bit, nr_next;

  /* loop over tetras */
  for(tetra_nr = 0; tetra_nr < Mesh.Ndt; tetra_nr++)
    {
      if(Mesh.DT[tetra_nr].t[0] < 0) /* skip deleted tetras */
        continue;

      /* edge flagging */
      bit     = 1;
      edge_nr = 0;

      /* loop over edges */
      while(Edge_visited[tetra_nr] != EDGE_ALL)
        {
          if((Edge_visited[tetra_nr] & bit) != 0)
            {
              bit <<= 1;
              edge_nr++;
              continue;
            }

          tetra *t = &DT[tetra_nr];

          /* edge-point relation */
          i = edge_start[edge_nr];
          j = edge_end[edge_nr];
          k = edge_opposite[edge_nr];
          l = edge_nexttetra[edge_nr];

          /* mark edge as visited */
          Edge_visited[tetra_nr] |= (1 << edge_nr);

          /* delaunay points on both side of face */
          dp_1 = t->p[i];
          dp_2 = t->p[j];

          /* skip large tetra */
          if(dp_1 < 0 || dp_2 < 0)
            {
              bit <<= 1;
              edge_nr++;
              continue;
            }

          /* skip ghost points (both local and foreign) */
          if((DP[dp_1].task != ThisTask || DP[dp_1].index < 0 || DP[dp_1].index >= NumGas) &&
             (DP[dp_2].task != ThisTask || DP[dp_2].index < 0 || DP[dp_2].index >= NumGas))
            {
              bit <<= 1;
              edge_nr++;
              continue;
            }

          /* count number of face vertices */
          count = 0;
          prev  = t;

          do
            {
              count++;
              next_tetra_nr = prev->t[l];
              next          = &DT[next_tetra_nr];

              for(m = 0, ll = ii = jj = -1; m < 4; m++)
                {
                  if(next->p[m] == prev->p[k])
                    ll = m;
                  if(next->p[m] == prev->p[i])
                    ii = m;
                  if(next->p[m] == prev->p[j])
                    jj = m;
                }

              if(ll < 0 || ii < 0 || jj < 0)
                terminate("inconsistency");

              kk = 6 - (ll + ii + jj);
              i  = ii;
              l  = ll;
              j  = jj;
              k  = kk;

              prev = next;
            }
          while(next != t);

          count++;

          /* get face normals (from both sides) */
          FaceNormals[Nnormals++] = (DP[dp_2].x - DP[dp_1].x);
          FaceNormals[Nnormals++] = (DP[dp_2].y - DP[dp_1].y);
          FaceNormals[Nnormals++] = (DP[dp_2].z - DP[dp_1].z);
          FaceNormals[Nnormals++] = (DP[dp_1].x - DP[dp_2].x);
          FaceNormals[Nnormals++] = (DP[dp_1].y - DP[dp_2].y);
          FaceNormals[Nnormals++] = (DP[dp_1].z - DP[dp_2].z);

          /* fill vertex entry list, first ID, count then tetra numbers */
          VertexEntries[Nentries++] = (int)DP[dp_1].ID;
          VertexEntries[Nentries++] = (int)DP[dp_2].ID;
          VertexEntries[Nentries++] = (int)count;
          VertexEntries[Nentries++] = (int)tetra_nr;

          /* get tetra indices of face vertices */
          count = 0;
          prev  = t;
          do
            {
              count++;
              next_tetra_nr = prev->t[l];
              next          = &DT[next_tetra_nr];

              VertexEntries[Nentries++] = (int)next_tetra_nr;

              for(m = 0, ll = ii = jj = -1; m < 4; m++)
                {
                  if(next->p[m] == prev->p[k])
                    ll = m;
                  if(next->p[m] == prev->p[i])
                    ii = m;
                  if(next->p[m] == prev->p[j])
                    jj = m;
                }

              if(ll < 0 || ii < 0 || jj < 0)
                terminate("inconsistency");

              kk = 6 - (ll + ii + jj);

              /* flag edge */
              for(nr_next = 0; nr_next < 6; nr_next++)
                if((edge_start[nr_next] == ii && edge_end[nr_next] == jj) || (edge_start[nr_next] == jj && edge_end[nr_next] == ii))
                  {
                    if((Edge_visited[next_tetra_nr] & (1 << nr_next)) && next != t)
                      terminate("inconsistency");

                    Edge_visited[next_tetra_nr] |= (1 << nr_next);
                    break;
                  }

              i = ii;
              l = ll;
              j = jj;
              k = kk;

              prev = next;

              if(Nentries > MaxEntries)
                terminate("Nentries > MaxEntries");

              if(Nnormals > MaxFaces)
                terminate("Nentries > MaxEntries");
            }
          while(next != t);

          bit <<= 1;
          edge_nr++;
        }
    }
}

/*! \brief Set Vertex coordinates in the respective array.
 *
 *  Copys the coordinates from the DTC array of the tessellation to a
 *  designated array VertexCoordinates.
 *
 *  \param[in] T Pointer to tessellation.
 *
 *  \return void
 */
void get_voronoi_face_vertex_coordinates(tessellation *T)
{
  int tetra_nr = 0;

  for(tetra_nr = 0; tetra_nr < T->Ndt; tetra_nr++)
    {
      VertexCoordinates[3 * Nvertices + 0] = T->DTC[tetra_nr].cx;
      VertexCoordinates[3 * Nvertices + 1] = T->DTC[tetra_nr].cy;
      VertexCoordinates[3 * Nvertices + 2] = T->DTC[tetra_nr].cz;
      Nvertices++;
    }
}

/*! \brief Function calls qsort for sorting faces by ID.
 *
 *  Uses compare_face_sort as comparison function. Requires array FaceSort.
 *
 *  \return void
 */
void sort_faces_by_ID(void)
{
  int i = 0, j = 0, k = 0;

  do
    {
      FaceSort[j].ID        = VertexEntries[i + 0];
      FaceSort[j].start     = i + 3;
      FaceSort[j].len       = VertexEntries[i + 2];
      FaceSort[j].normal[0] = FaceNormals[k++];
      FaceSort[j].normal[1] = FaceNormals[k++];
      FaceSort[j].normal[2] = FaceNormals[k++];
      j++;

      FaceSort[j].ID        = VertexEntries[i + 1];
      FaceSort[j].start     = i + 3;
      FaceSort[j].len       = VertexEntries[i + 2];
      FaceSort[j].normal[0] = FaceNormals[k++];
      FaceSort[j].normal[1] = FaceNormals[k++];
      FaceSort[j].normal[2] = FaceNormals[k++];
      j++;

      i += 3 + VertexEntries[i + 2];

      if(j > MaxFaces)
        terminate("j > MaxFaces");
    }
  while(i < Nentries);

  Nsort = j;

  /* sort faces by ID */
  qsort(FaceSort, Nsort, sizeof(struct data_face_sort), compare_face_sort);
}

/*! \brief Outputs Voronoi vertex indices to file.
 *
 *  Outputs the Voronoi vertex indices from task writeTask to lastTask in file
 *  fname.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] fname1 File name of file index data is written in.
 *  \param[in] fname2 File name of file face data is written in.
 *  \param[in] writeTask Task that gathers information and writes data.
 *  \param[in] lastTask Last task that is included in this dump.
 *
 *  \return void
 */
void write_voronoi_face_vertex_indices(tessellation *T, char *fname1, char *fname2, int writeTask, int lastTask)
{
  FILE *fd1, *fd2;
  MPI_Status status;
  int nVertices_tot, nEntries_tot, nNormals_tot;
  int nVertices_before, i, task, *tmp;
  int *Nvertices_list, *Nentries_list, *Nnormals_list, *Nsort_list;
  struct data_face_sort *tmp_sort;

  VertexEntries = mymalloc("VertexEntries", MaxEntries * sizeof(int));
  FaceNormals   = mymalloc("VertexEntries", MaxFaces * sizeof(int));

  /* get faces */
  get_voronoi_face_vertex_indices(T);

  FaceSort = (struct data_face_sort *)mymalloc("face_sort", sizeof(struct data_face_sort) * MaxFaces);

  /* sort faces */
  sort_faces_by_ID();

  Nentries = 0;
  for(i = 0; i < Nsort; i++)
    Nentries += FaceSort[i].len + 2;

  /* I/O */
  Nvertices_list = mymalloc("Nvertices_list", sizeof(int) * NTask);
  Nentries_list  = mymalloc("Nentries_list", sizeof(int) * NTask);
  Nsort_list     = mymalloc("Nsort_list", sizeof(int) * NTask);
  Nnormals_list  = mymalloc("Nnormals_list", sizeof(int) * NTask);

  if(ThisTask == writeTask)
    {
      nVertices_tot = Nvertices;
      nEntries_tot  = Nentries;
      nNormals_tot  = Nnormals;
      for(task = writeTask + 1; task <= lastTask; task++)
        {
          MPI_Recv(&Nvertices_list[task], 1, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
          MPI_Recv(&Nentries_list[task], 1, MPI_INT, task, TAG_LOCALN + 1, MPI_COMM_WORLD, &status);
          MPI_Recv(&Nsort_list[task], 1, MPI_INT, task, TAG_LOCALN + 2, MPI_COMM_WORLD, &status);
          MPI_Recv(&Nnormals_list[task], 1, MPI_INT, task, TAG_LOCALN + 3, MPI_COMM_WORLD, &status);
          MPI_Send(&nVertices_tot, 1, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
          nVertices_tot += Nvertices_list[task];
          nEntries_tot += Nentries_list[task];
          nNormals_tot += Nnormals_list[task];
        }
      if(!(fd1 = fopen(fname1, "w")))
        terminate("I/O error");

      if(!(fd2 = fopen(fname2, "w")))
        terminate("I/O error");

      my_fwrite(&nEntries_tot, sizeof(int), 1, fd1);
      my_fwrite(&nNormals_tot, sizeof(int), 1, fd2);
      for(i = 0; i < Nsort; i++)
        {
          my_fwrite(&FaceSort[i].ID, sizeof(int), 1, fd1);
          my_fwrite(&FaceSort[i].len, sizeof(int), 1, fd1);
          my_fwrite(&VertexEntries[FaceSort[i].start], sizeof(int) * FaceSort[i].len, 1, fd1);
          my_fwrite(FaceSort[i].normal, 3 * sizeof(float), 1, fd2);
        }

      for(task = writeTask + 1; task <= lastTask; task++)
        {
          tmp_sort = (struct data_face_sort *)mymalloc("tmp_sort", sizeof(struct data_face_sort) * Nsort_list[task]);
          tmp      = mymalloc("tmp", sizeof(int) * Nentries_list[task]);
          MPI_Recv(tmp, Nentries_list[task], MPI_INT, task, TAG_N + 1, MPI_COMM_WORLD, &status);
          MPI_Recv(tmp_sort, Nsort_list[task] * sizeof(struct data_face_sort), MPI_BYTE, task, TAG_N + 2, MPI_COMM_WORLD, &status);

          for(i = 0; i < Nsort_list[task]; i++)
            {
              my_fwrite(&tmp_sort[i].ID, sizeof(int), 1, fd1);
              my_fwrite(&tmp_sort[i].len, sizeof(int), 1, fd1);
              my_fwrite(&tmp[tmp_sort[i].start], sizeof(int) * tmp_sort[i].len, 1, fd1);
              my_fwrite(tmp_sort[i].normal, 3 * sizeof(float), 1, fd2);
            }
          myfree(tmp);
          myfree(tmp_sort);
        }
      fclose(fd2);
      fclose(fd1);
    }
  else
    {
      MPI_Send(&Nvertices, 1, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
      MPI_Send(&Nentries, 1, MPI_INT, writeTask, TAG_LOCALN + 1, MPI_COMM_WORLD);
      MPI_Send(&Nsort, 1, MPI_INT, writeTask, TAG_LOCALN + 2, MPI_COMM_WORLD);
      MPI_Send(&Nnormals, 1, MPI_INT, writeTask, TAG_LOCALN + 3, MPI_COMM_WORLD);
      MPI_Recv(&nVertices_before, 1, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
      for(i = 0; i < Nentries; i++)
        if(VertexEntries[i] >= 0)
          VertexEntries[i] += nVertices_before;
      MPI_Send(VertexEntries, Nentries, MPI_INT, writeTask, TAG_N + 1, MPI_COMM_WORLD);
      MPI_Send(FaceSort, Nsort * sizeof(struct data_face_sort), MPI_BYTE, writeTask, TAG_N + 2, MPI_COMM_WORLD);
    }

  myfree(Nnormals_list);
  myfree(Nsort_list);
  myfree(Nentries_list);
  myfree(Nvertices_list);
  myfree(FaceSort);
  myfree(FaceNormals);
  myfree(VertexEntries);
}

/*! \brief Outputs Voronoi vertex coordinates to file.
 *
 *  Outputs the Voronoi vertex coordinates from task write Task to lastTask in
 *  file fname.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] fname File name of file the data is written in.
 *  \param[in] writeTask Task that gathers information and writes data.
 *  \param[in] lastTask Last task that is included in this dump.
 *
 *  \return void
 */
void write_voronoi_face_vertex_coordinates(tessellation *T, char *fname, int writeTask, int lastTask)
{
  FILE *fd;
  MPI_Status status;
  int *Nvertices_list;
  int nVertices_tot, task;
  float *tmp;

  VertexCoordinates = mymalloc("VertexCoordinates", MaxEntries * 3 * sizeof(float));

  /* get coordinates */
  get_voronoi_face_vertex_coordinates(T);

  /* I/O */
  Nvertices_list = mymalloc("Nvertices_list", sizeof(int) * NTask);
  if(ThisTask == writeTask)
    {
      nVertices_tot = Nvertices;
      for(task = writeTask + 1; task <= lastTask; task++)
        {
          MPI_Recv(&Nvertices_list[task], 1, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
          nVertices_tot += Nvertices_list[task];
        }

      if(!(fd = fopen(fname, "w")))
        terminate("I/O error");

      my_fwrite(&nVertices_tot, sizeof(int), 1, fd);
      my_fwrite(VertexCoordinates, sizeof(float), 3 * Nvertices, fd);
      for(task = writeTask + 1; task <= lastTask; task++)
        {
          tmp = mymalloc("tmp", 3 * sizeof(float) * Nvertices_list[task]);
          MPI_Recv(tmp, 3 * Nvertices_list[task], MPI_FLOAT, task, TAG_N + 1, MPI_COMM_WORLD, &status);
          my_fwrite(tmp, sizeof(float), 3 * Nvertices_list[task], fd);
          myfree(tmp);
        }
      fclose(fd);
    }
  else
    {
      MPI_Send(&Nvertices, 1, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
      MPI_Send(VertexCoordinates, 3 * Nvertices, MPI_FLOAT, writeTask, TAG_N + 1, MPI_COMM_WORLD);
    }
  myfree(Nvertices_list);
  myfree(VertexCoordinates);
}

/*! \brief Outputs Voronoi mesh to file.
 *
 *  Outputs the Voronoi mesh data from task write Task to lastTask in file
 *  fname.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] fname File name of file the data is written in.
 *  \param[in] writeTask Task that gathers information and writes data.
 *  \param[in] lastTask Last task that is included in this dump.
 *
 *  \return void
 */
void write_voronoi_mesh(tessellation *T, char *fname, int writeTask, int lastTask)
{
  char buf1[255], buf2[255];

  MaxEntries = 1000 * NumGas;
  MaxFaces   = 100 * NumGas;

  /* coordinates */
  Nvertices = 0;
  sprintf(buf1, "%s_coordinates.dat", fname);
  write_voronoi_face_vertex_coordinates(T, buf1, writeTask, lastTask);

  /* indices */
  Edge_visited = mymalloc_movable(&Edge_visited, "Edge_visited", Mesh.Ndt * sizeof(unsigned char));
  int i;
  for(i = 0; i < Mesh.Ndt; i++)
    Edge_visited[i] = 0;

  Nentries = 0;
  Nnormals = 0;
  sprintf(buf1, "%s_indices.dat", fname);
  sprintf(buf2, "%s_normals.dat", fname);
  write_voronoi_face_vertex_indices(T, buf1, buf2, writeTask, lastTask);
  myfree(Edge_visited);
}

#endif /* #if !defined(TWODIMS) && !defined(ONEDIMS) */
