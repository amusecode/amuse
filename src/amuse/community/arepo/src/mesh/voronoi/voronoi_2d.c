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
 * \file        src/mesh/voronoi/voronoi_2d.c
 * \date        05/2018
 * \brief       Routines to build a 2d Voronoi mesh.
 * \details     Note that some of these routines have the same name as the ones
 *              in voronoi_1d.c and voronoi_3d.c and just replace them in case
 *              the Config-option TWODIMS is active. This is also the reason
 *              why some of these functions are empty but nonetheless have to
 *              exist in this file.
 *              contains functions:
 *                void initialize_and_create_first_tetra(tessellation * T)
 *                int insert_point(tessellation * T, int pp, int ttstart)
 *                void make_a_2_to_4_flip(tessellation * T, int pp, int tt0,
 *                  int tt1, int tt2, int tt3, int i0, int j0)
 *                void make_a_1_to_3_flip(tessellation * T, int pp, int tt0,
 *                  int tt1, int tt2)
 *                void check_edge_and_flip_if_needed(tessellation * T, int ip,
 *                  int it)
 *                int get_triangle(tessellation * T, int pp, int *moves, int
 *                  *degenerate_flag, int ttstart)
 *                static inline void add_row_2d(double *m, int r1, int r2,
 *                  double fac)
 *                int solve_linear_equations_2d(double *m, double *res)
 *                int FindTriangle(tessellation * T, int tt, int pp,
 *                  int *degnerate_flag, int *nexttetra)
 *                int InCircle_Quick(tessellation * T, int pp0, int pp1,
 *                  int pp2, int pp)
 *                int InCircle_Errorbound(tessellation * T, int pp0, int pp1,
 *                  int pp2, int pp)
 *                int InCircle_Exact(tessellation * T, int pp0, int pp1,
 *                  int pp2, int pp)
 *                double test_triangle_orientation(tessellation * T, int pp0,
 *                  int pp1, int pp2)
 *                int Orient2d_Quick(tessellation * T, int pp0, int pp1,
 *                  int pp2)
 *                int Orient2d_Exact(tessellation * T, int pp0, int pp1,
 *                  int pp2)
 *                void process_edge_faces_and_volumes(tessellation * T, int tt,
 *                  int nr)
 *                int derefine_refine_get_triangles(tessellation * T, int tt,
 *                  int nr, point * dtip, triangle * trilist, int ntri,
 *                  int max_n_tri)
 *                int derefine_add_point_and_split_tri(int q, triangle
 *                  * trilist, int ntri, int max_ntri, double vol)
 *                double get_tri_volume(int i, triangle * trilist)
 *                void derefine_refine_process_edge(tessellation * T, double
 *                  *vol, int tt, int nr)
 *                void compute_circumcircles(tessellation * T)
 *                void update_circumcircle(tessellation * T, int tt)
 *                void set_integers_for_pointer(point * p)
 *                void write_voronoi_mesh(tessellation * T, char *fname, int
 *                  writeTask, int lastTask)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 21.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../main/allvars.h"
#include "../../main/proto.h"

#include "voronoi.h"

#if defined(TWODIMS) && !defined(ONEDIMS) /* will only be compiled in 2D case */

#define INSIDE_EPS 1.0e-8
#define GAUSS_EPS 1.0e-8

/*! \brief Initializes 2d tessellation and create all-enclosing triangle.
 *
 *  \param[out] T Pointer to tessellation structure which is set and its arrays
 *              are allocated in this routine.
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
  T->Nvf = 0;
  T->Ndt = 0;

  T->VF = mymalloc_movable(&T->VF, "VF", T->MaxNvf * sizeof(face));

  T->DP = mymalloc_movable(&T->DP, "DP", (T->MaxNdp + 5) * sizeof(point));
  T->DP += 5;

  T->DT = mymalloc_movable(&T->DT, "DT", T->MaxNdt * sizeof(tetra));

  /* construct all encompassing huge triangle */
  double box, tetra_incircle, tetra_sidelength, tetra_height;

  box = boxSize_X;
  if(box < boxSize_Y)
    box = boxSize_Y;

  box *= 1.05;

  tetra_incircle = 2.001 * (1 + sqrt(3)) / 3.0 * box; /* to give room for ghost particles needed for periodic/reflective
                                                         boundary conditions, the incircle is twice as large, i.e.
                                                         [-0.5*box, 1.5*box,-0.5*box, 1.5*box] should be inside triangle */
  tetra_sidelength = tetra_incircle * sqrt(12);
  tetra_height     = sqrt(3.0) / 2 * tetra_sidelength;

  if(ThisTask == 0)
    printf("side-length of enclosing triangle=%g tetra_height=%g box=%g\n", tetra_sidelength, tetra_height, box);

  point *DP = T->DP;
  tetra *DT = T->DT;

  /* first, let's make the points */
  DP[-3].x = 0.5 * tetra_sidelength;
  DP[-3].y = -1.0 / 3 * tetra_height;
  DP[-3].z = 0;

  DP[-2].x = 0;
  DP[-2].y = 2.0 / 3 * tetra_height;
  DP[-2].z = 0;

  DP[-1].x = -0.5 * tetra_sidelength;
  DP[-1].y = -1.0 / 3 * tetra_height;
  DP[-1].z = 0;

  for(i = -3; i <= -1; i++)
    {
      DP[i].x += 0.5 * box;
      DP[i].y += 1.0 / 3 * tetra_height - 0.5 * box;
    }

  for(i = -3, p = &DP[-3]; i < 0; i++, p++)
    {
      p->index   = -1;
      p->task    = ThisTask;
      p->timebin = 0;
    }

  /* we also define a neutral element at infinity */
  DPinfinity = -4;

  DP[DPinfinity].x       = MAX_DOUBLE_NUMBER;
  DP[DPinfinity].y       = MAX_DOUBLE_NUMBER;
  DP[DPinfinity].z       = MAX_DOUBLE_NUMBER;
  DP[DPinfinity].index   = -1;
  DP[DPinfinity].task    = ThisTask;
  DP[DPinfinity].timebin = 0;

  /* now let's make the big triangle */
  DT[0].p[0] = -3;
  DT[0].p[1] = -2;
  DT[0].p[2] = -1;

  /* On the outer faces, we attach tetrahedra with the neutral element as tip.
   * This way we will be able to navigate nicely within the tesselation,
   * and all tetrahedra have defined neighbouring tetrahedra.
   */

  for(i = 0; i < 3; i++)
    {
      n = i + 1; /* tetra index */

      DT[0].t[i] = n;
      DT[0].s[i] = 2;

      DT[n].t[2] = 0;
      DT[n].s[2] = i;
      DT[n].p[2] = DPinfinity;
    }

  DT[1].p[0] = DT[0].p[2];
  DT[1].p[1] = DT[0].p[1];

  DT[2].p[0] = DT[0].p[0];
  DT[2].p[1] = DT[0].p[2];

  DT[3].p[0] = DT[0].p[1];
  DT[3].p[1] = DT[0].p[0];

  DT[1].t[0] = 3;
  DT[3].t[1] = 1;
  DT[1].s[0] = 1;
  DT[3].s[1] = 0;

  DT[1].t[1] = 2;
  DT[2].t[0] = 1;
  DT[1].s[1] = 0;
  DT[2].s[0] = 1;

  DT[2].t[1] = 3;
  DT[3].t[0] = 2;
  DT[2].s[1] = 0;
  DT[3].s[0] = 1;

  T->Ndt = 4; /* we'll start out with 4 triangles */

  CentralOffsetX = 0.5 * box - 0.5000001 * tetra_sidelength;
  CentralOffsetY = -0.5000001 * box;

  ConversionFac = 1.0 / (1.001 * tetra_sidelength);

  for(i = -3; i < 0; i++)
    set_integers_for_point(T, i);
}

/*! \brief Insert a point into mesh.
 *
 *  Finds the triangle that contains this point, splits the triangle (usually
 *  into three). After this, flip the edges if needed restore
 *  Delaunayhood (which is applied recursively) until a valid Delaunay mesh
 *  is restored.
 *
 *  \param[in, out] T Pointer to tessellation.
 *  \param[in] pp Index of Delaunay point in DP array.
 *  \param[in] ttstart Initial guess in which triangle it might be,
 *             index in DT array.
 *
 * \return Index of triangle containing point pp.
 */
int insert_point(tessellation *T, int pp, int ttstart)
{
  int tt0, tt1, tt2, tt3, ttetra_with_p;
  int moves, degenerate_flag;

  /* first, need to do a point location */
  tt0 = get_triangle(T, pp, &moves, &degenerate_flag, ttstart);

  ttetra_with_p = tt0;

  if(degenerate_flag == 1) /* that's the normal split of a triangle into 3 */
    {
      /* we now need to split this triangle into three  */
      tt1 = T->Ndt++;
      tt2 = T->Ndt++;

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

      T->DT[tt1] = T->DT[tt0];
      T->DT[tt2] = T->DT[tt0];

      make_a_1_to_3_flip(T, pp, tt0, tt1, tt2);

      T->DTF[tt0] = 0;
      T->DTF[tt1] = 0;
      T->DTF[tt2] = 0;

      check_edge_and_flip_if_needed(T, pp, tt0);
      check_edge_and_flip_if_needed(T, pp, tt1);
      check_edge_and_flip_if_needed(T, pp, tt2);
    }
  else
    {
      degenerate_flag -= 10;

      tt1 = T->DT[tt0].t[degenerate_flag];

      /* we now need to split this into two triangles */
      tt2 = T->Ndt++;
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

      T->DT[tt2] = T->DT[tt0];
      T->DT[tt3] = T->DT[tt1];

      make_a_2_to_4_flip(T, pp, tt0, tt1, tt2, tt3, degenerate_flag, T->DT[tt0].s[degenerate_flag]);

      T->DTF[tt0] = 0;
      T->DTF[tt1] = 0;
      T->DTF[tt2] = 0;
      T->DTF[tt3] = 0;

      check_edge_and_flip_if_needed(T, pp, tt0);
      check_edge_and_flip_if_needed(T, pp, tt1);
      check_edge_and_flip_if_needed(T, pp, tt2);
      check_edge_and_flip_if_needed(T, pp, tt3);
    }

  return ttetra_with_p;
}

/*! \brief Make a 2 to 4 flip needed if point is on edge of a Delaunay
 *         triangle.
 *
 *  If a new point is at the edge of a Delaunay triangle, both adjacent
 *  triangles need to be split into two. See Springel (2010) for a
 *  detailed discussion.
 *
 *  \param[in, out] T Pointer to tessellation.
 *  \param[in] pp Index of Delaunay point in DP array.
 *  \param[in] tt0 Index of point 0 in DT array.
 *  \param[in] tt1 Index of point 1 in DT array.
 *  \param[in] tt2 Index of point 2 in DT array.
 *  \param[in] tt3 Index of point 3 in DT array.
 *  \param[in] i0 Index (in DT->s) of point opposite to common face that needs
 *             to be involved in flip.
 *  \param[in] j0 Second Index (in DT->s) of point opposite to common face that
 *             needs to be involved in flip.
 *
 *  \return void
 */
void make_a_2_to_4_flip(tessellation *T, int pp, int tt0, int tt1, int tt2, int tt3, int i0, int j0)
{
  tetra *DT = T->DT;
  tetra *t0 = &DT[tt0];
  tetra *t1 = &DT[tt1];
  tetra *t2 = &DT[tt2];
  tetra *t3 = &DT[tt3];

  int i1, i2, j1, j2;

  CountFlips++;
  Count_2_to_4_Flips2d++;

  i1 = i0 + 1;
  i2 = i0 + 2;
  j1 = j0 + 1;
  j2 = j0 + 2;

  if(i1 > 2)
    i1 -= 3;
  if(i2 > 2)
    i2 -= 3;

  if(j1 > 2)
    j1 -= 3;
  if(j2 > 2)
    j2 -= 3;

  t0->p[i1] = pp;
  t1->p[j2] = pp;
  t2->p[i2] = pp;
  t3->p[j1] = pp;

  t0->t[i0] = tt1;
  t1->t[j0] = tt0;
  t0->s[i0] = j0;
  t1->s[j0] = i0;

  t1->t[j1] = tt3;
  t3->t[j2] = tt1;
  t1->s[j1] = j2;
  t3->s[j2] = j1;

  t2->t[i1] = tt0;
  t0->t[i2] = tt2;
  t2->s[i1] = i2;
  t0->s[i2] = i1;

  t2->t[i0] = tt3;
  t3->t[j0] = tt2;
  t2->s[i0] = j0;
  t3->s[j0] = i0;

  DT[t0->t[i1]].t[t0->s[i1]] = tt0;
  DT[t1->t[j2]].t[t1->s[j2]] = tt1;
  DT[t2->t[i2]].t[t2->s[i2]] = tt2;
  DT[t3->t[j1]].t[t3->s[j1]] = tt3;
}

/*! \brief Makes a 1 to 3 flip needed if point is in a Delaunay triangle.
 *
 *  If a new point is in a Delaunay triangle, this
 *  triangles need to be split into three.
 *
 *  \param[in, out] T Pointer to tessellation.
 *  \param[in] pp Index of Delaunay point in DP array.
 *  \param[in] tt0 Index of point 0 in DT array.
 *  \param[in] tt1 Index of point 1 in DT array.
 *  \param[in] tt2 Index of point 2 in DT array.
 *
 *  \return void
 */
void make_a_1_to_3_flip(tessellation *T, int pp, int tt0, int tt1, int tt2)
{
  tetra *DT = T->DT;
  tetra *t0 = &DT[tt0];
  tetra *t1 = &DT[tt1];
  tetra *t2 = &DT[tt2];

  CountFlips++;
  Count_1_to_3_Flips2d++;

  t0->p[0] = pp;
  t1->p[1] = pp;
  t2->p[2] = pp;

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

  DT[t0->t[0]].t[t0->s[0]] = tt0;
  DT[t1->t[1]].t[t1->s[1]] = tt1;
  DT[t2->t[2]].t[t2->s[2]] = tt2;
}

/*! \brief Flips trangle if needed.
 *
 *  See Springel (2010) for detailed discussion how mesh is constructed.
 *
 *  \param[in, out] T Pointer to tessellation.
 *  \param[in] ip Index to Delaunay point, DP array.
 *  \param[in] it Index to corner of triangle, DT array.
 *
 *  \return void
 */
void check_edge_and_flip_if_needed(tessellation *T, int ip, int it)
{
  tetra *DT = T->DT;
  char *DTF = T->DTF;

  tetra *t = &DT[it];

  int tt, pp, t0, t2;
  int pi, pi1, pi2;
  int ni, ni1, ni2;
  int st2, st0;

  if(t->p[0] == ip)
    pi = 0;
  else if(t->p[1] == ip)
    pi = 1;
  else
    pi = 2;

  /* get the point that lies accross the edge to obtain the quadriliteral */

  tt = t->t[pi];
  ni = t->s[pi];
  pp = DT[tt].p[ni];

  int ret, ret_exact;

  ret = InCircle_Errorbound(T, t->p[0], t->p[1], t->p[2], pp);
  CountInSphereTests++;

  if(ret != 0)
    ret_exact = ret;
  else
    {
      ret_exact = InCircle_Exact(T, t->p[0], t->p[1], t->p[2], pp);
      CountInSphereTestsExact++;
    }

  if(ret_exact > 0)
    {
      /* pp lies in the triangle, the edge is not Delaunay. Need to do a flip */

      CountFlips++;

      ni1 = ni + 1;
      if(ni1 > 2)
        ni1 -= 3;
      ni2 = ni + 2;
      if(ni2 > 2)
        ni2 -= 3;

      pi1 = pi + 1;
      if(pi1 > 2)
        pi1 -= 3;
      pi2 = pi + 2;
      if(pi2 > 2)
        pi2 -= 3;

      t0 = DT[tt].t[ni1];
      t2 = t->t[pi1];

      st0 = DT[tt].s[ni1];
      st2 = t->s[pi1];

      /* change the points of the triangles */
      t->p[pi2]     = pp;
      DT[tt].p[ni2] = ip;

      /* change the pointers to the neighbouring triangles, and fix
         the adjency relations */

      t->t[pi1]     = tt;
      DT[tt].t[ni1] = it;
      t->s[pi1]     = ni1;
      DT[tt].s[ni1] = pi1;

      t->t[pi]      = t0;
      DT[t0].t[st0] = it;
      t->s[pi]      = st0;
      DT[t0].s[st0] = pi;

      DT[tt].t[ni]  = t2;
      DT[t2].t[st2] = tt;
      DT[tt].s[ni]  = st2;
      DT[t2].s[st2] = ni;

      DTF[tt] = 0;
      DTF[it] = 0;

      /* now we need to test also the two sides opposite of p */
      check_edge_and_flip_if_needed(T, ip, it);
      check_edge_and_flip_if_needed(T, ip, tt);
    }
}

/*! \brief Finds triangle in which new Delaunay point is located.
 *
 *  Starts with a suggested triangle ttstart and checks if the point is
 *  contained in this triangle. If not, the procedure is repeated for the
 *  neighboring triangle.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] pp Index of Delaunay point in DP array.
 *  \param[out] moves Number of iterations to find the correct triangle.
 *  \param[out] degenerate_flag Flag if point lies on edge of a triangle.
 *  \param[in] ttstart Starting index for the search for the correct triangle.
 *
 *  \return Index of triangle in DT array.
 */
int get_triangle(tessellation *T, int pp, int *moves, int *degenerate_flag, int ttstart)
{
  int count_moves = 0;
  int ret;
  int tt, next_tetra;

  tt = ttstart;

#define MAX_COUNT_MOVES 1000000

  while((ret = FindTriangle(T, tt, pp, degenerate_flag, &next_tetra)) == 0)
    {
      /* we need to see in which of the three possible neighbouring triangles
         we should walk. We'll choose the one which lies along the face that
         is traversed by a line from the cm of the triangle to the point in
         question.
       */
      count_moves++;

      if(count_moves > MAX_COUNT_MOVES)
        {
          printf("ta=%d triangle=%d  xy=(%g|%g) ID=%d\n", ThisTask, (int)(tt), T->DP[pp].x, T->DP[pp].y, T->DP[pp].ID);
          if(count_moves > MAX_COUNT_MOVES + 10)
            terminate("too many moves, problem to find triangle");
        }

      tt = next_tetra;
    }

  *moves = count_moves;

  return tt;
}

/*! \brief Add row in matrix equation.
 *
 *  Auxiliary function for solve_linear_equations_2d.
 *
 *  \param[in, out] m Matrix.
 *  \param[in] r1 Index of row to be modified.
 *  \param[in] r2 Index of row which is added to r1.
 *  \param[in] fac Factor by which row r2 is multiplied before adding to r1.
 *
 *  \return void
 */
static inline void add_row_2d(double *m, int r1, int r2, double fac)
{
  int i;

  for(i = 0; i < 3; i++)
    m[r1 * 3 + i] += fac * m[r2 * 3 + i];
}

/*! \brief Solve system of linear equations for 2d Voronoi construction.
 *
 *  This is needed in get_triangle routine.
 *
 *  \param[in, out] m Matrix.
 *  \param[in, out] res Array for result.
 *
 *  \return 0 if success, -1 else.
 */
int solve_linear_equations_2d(double *m, double *res)
{
  int ix, iy;

  if(fabs(m[0]) > fabs(m[3]))
    {
      ix = 0;
      iy = 1;
    }
  else
    {
      ix = 1;
      iy = 0;
    }

  add_row_2d(m, iy, ix, -m[iy * 3] / m[ix * 3]);

  res[1] = m[iy * 3 + 2] / m[iy * 3 + 1];
  res[0] = (m[ix * 3 + 2] - res[1] * m[ix * 3 + 1]) / m[ix * 3];

  if(fabs(m[ix * 3]) < 1.0e-12)
    return -1;

  return 0;
}

/*! \brief Does point lie in triangle?
 *
 *  Tests whether point pp lies in the triangle, on an edge, or outside. In the
 *  latter case, a neighboring triangle is returned. First, a fast search is
 *  performed and if this yields that point might be on an edge, a (more
 *  expensive) exact determination is performed.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] tt Index of triangle in DT array.
 *  \param[in] pp Index of Delaunay point in DP array.
 *  \param[out] degenerate_flag Flag if point lies on edge of a triangle.
 *  \param[out] nexttetra Index of neighboring triangle in direction of point.
 *
 *  \return 1: point inside triangle; 0 outside; 10,11,12: on edge.
 */
int FindTriangle(tessellation *T, int tt, int pp, int *degnerate_flag, int *nexttetra)
{
  tetra *DT = T->DT;
  point *DP = T->DP;
  tetra *t  = &DT[tt];
  point *p  = &DP[pp];

  int pp0, pp1, pp2;
  point *p0, *p1, *p2;

  pp0 = t->p[0];
  pp1 = t->p[1];
  pp2 = t->p[2];

  p0 = &DP[pp0];
  p1 = &DP[pp1];
  p2 = &DP[pp2];

  if(pp0 == DPinfinity || pp1 == DPinfinity || pp2 == DPinfinity)
    {
      char buf[1000];
      sprintf(buf, "we are in a triangle with an infinity point. tetra=%d  p=(%g|%g)\n", (int)(tt), p->x, p->y);
      terminate(buf);
    }

  Count_InTetra++;

  double ax = p1->xx - p0->xx;
  double ay = p1->yy - p0->yy;

  double bx = p2->xx - p0->xx;
  double by = p2->yy - p0->yy;

  double qx = p->xx - p0->xx;
  double qy = p->yy - p0->yy;

  double mv_data[] = {ax, bx, qx, ay, by, qy};
  double x[2];

  int ivol, flag2, flag1, flag0;
  int count_zeros = 0;

  int status;

  status = solve_linear_equations_2d(mv_data, x);

  if(status < 0)
    {
      ivol = Orient2d_Exact(T, t->p[0], t->p[1], t->p[2]);
      if(ivol <= 0)
        {
          char buf[1000];
          sprintf(buf, "flat or negatively triangle found (ivol=%d)\n", ivol);
          terminate(buf);
        }
    }

  if(status >= 0)
    {
      if(x[0] > INSIDE_EPS && x[1] > INSIDE_EPS && (1 - (x[0] + x[1])) > INSIDE_EPS)
        {
          /* looks like we are safely inside the triangle */

          *degnerate_flag = 1;
          return 1;
        }

      if(x[0] < -INSIDE_EPS || x[1] < -INSIDE_EPS || (1 - (x[0] + x[1])) < -INSIDE_EPS)
        {
          /* looks like we are clearly outside the triangle.
             Let's look for a good neighbouring triangle to continue the search */

          /* note: in the (a,b) basis, the center-of-mass has coordinates (1/3, 1/3) */

          double w, u;

          if(fabs(x[1] - (1.0 / 3)) > INSIDE_EPS)
            {
              w = (1.0 / 3) / ((1.0 / 3) - x[1]);
              if(w > 0)
                {
                  u = (1.0 / 3) + w * (x[0] - (1.0 / 3));
                  if(u > -INSIDE_EPS && (1 - u) > -INSIDE_EPS)
                    {
                      *nexttetra = t->t[2];
                      return 0;
                    }
                }
            }

          if(fabs(x[0] - (1.0 / 3)) > INSIDE_EPS)
            {
              w = (1.0 / 3) / ((1.0 / 3) - x[0]);
              if(w > 0)
                {
                  u = (1.0 / 3) + w * (x[1] - (1.0 / 3));
                  if(u > -INSIDE_EPS && (1 - u) > -INSIDE_EPS)
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
     whether we think the point lies on an edge of the triangle */

  Count_InTetraExact++;

  ivol = Orient2d_Exact(T, t->p[0], t->p[1], t->p[2]);

  if(ivol <= 0)
    {
      char buf[1000];
      sprintf(buf, "flat or negatively oriented triangle found (ivol=%d)\n", ivol);
      terminate(buf);
    }

  flag0 = Orient2d_Exact(T, pp1, pp2, pp);
  flag1 = Orient2d_Exact(T, pp2, pp0, pp);
  flag2 = Orient2d_Exact(T, pp0, pp1, pp);

  if(flag0 == 0)
    count_zeros++;

  if(flag1 == 0)
    count_zeros++;

  if(flag2 == 0)
    count_zeros++;

  if(count_zeros >= 2)
    {
      printf("flags=%d %d %d\n", flag0, flag1, flag2);

      printf("points: %d %d %d %d\n", (int)(pp0), (int)(pp1), (int)(pp2), (int)(pp));
      printf("Ngas=%d\n", NumGas);
      printf("xyz, p=%d: (%g|%g)  index=%d task=%d ID=%d  flags\n", (int)(pp0), p0->x, p0->y, p0->index, p0->task,
             P[p0->index % NumGas].ID);
      printf("xyz, p=%d: (%g|%g)  index=%d task=%d ID=%d  flags\n", (int)(pp1), p1->x, p1->y, p1->index, p1->task,
             P[p1->index % NumGas].ID);
      printf("xyz, p=%d: (%g|%g)  index=%d task=%d ID=%d  flags\n", (int)(pp2), p2->x, p2->y, p2->index, p2->task,
             P[p2->index % NumGas].ID);
      printf("xyz, p=%d: (%g|%g)  index=%d task=%d ID=%d  flags\n", (int)(pp), p->x, p->y, p->index, p->task, P[p->index % NumGas].ID);
      terminate("too many zeros - (perhaps identical points inserted?)");
    }

  if(flag0 >= 0 && flag1 >= 0 && flag2 >= 0)
    {
      /* we have a point inside the triangle, but it may still be on one of the edges */

      if(count_zeros == 0)
        {
          /* ok, we are inside */
          *degnerate_flag = 1;
          return 1;
        }

      if(count_zeros == 1) /* we lie on a face */
        {
          if(flag2 == 0)
            {
              *degnerate_flag = 12;
              return 12; /* point lies on side A */
            }
          if(flag1 == 0)
            {
              *degnerate_flag = 11;
              return 11; /* point lies on side C */
            }

          if(flag0 == 0)
            {
              *degnerate_flag = 10;
              return 10; /* point lies on side B */
            }
        }
    }

  /* we are clearly outside, let's select the suitable neighbour */

  if(flag0 < 0 && flag1 >= 0 && flag2 >= 0)
    {
      *nexttetra = t->t[0];
      return 0;
    }

  if(flag0 >= 0 && flag1 < 0 && flag2 >= 0)
    {
      *nexttetra = t->t[1];
      return 0;
    }

  if(flag0 >= 0 && flag1 >= 0 && flag2 < 0)
    {
      *nexttetra = t->t[2];
      return 0;
    }

  /* there are apparently two negative values. Let's pick a random one */

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

  *nexttetra = t->t[ind];
  return 0;
}

/*! \brief Tests whether point pp lies in the circumcircle around triangle
 *        p0,p1,p2.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] pp0 Index in DP of first point in triangle.
 *  \param[in] pp1 Index in DP of second point in triangle.
 *  \param[in] pp2 Index in DP of third point in triangle.
 *  \param[in] pp Index in DP of point to be checked.
 *
 *  \return (-1,0,1); -1: in circle; 0 on circle, 1: outside circle.
 */
int InCircle_Quick(tessellation *T, int pp0, int pp1, int pp2, int pp)
{
  point *DP = T->DP;
  point *p0 = &DP[pp0];
  point *p1 = &DP[pp1];
  point *p2 = &DP[pp2];
  point *p  = &DP[pp];

  double ax, ay, bx, by, cx, cy;
  double ab, bc, ca, a2, b2, c2, x;

  if(pp0 == DPinfinity || pp1 == DPinfinity || pp2 == DPinfinity || pp == DPinfinity)
    return -1;

  ax = p0->xx - p->xx;
  ay = p0->yy - p->yy;
  bx = p1->xx - p->xx;
  by = p1->yy - p->yy;
  cx = p2->xx - p->xx;
  cy = p2->yy - p->yy;

  ab = ax * by - bx * ay;
  bc = bx * cy - cx * by;
  ca = cx * ay - ax * cy;

  a2 = ax * ax + ay * ay;
  b2 = bx * bx + by * by;
  c2 = cx * cx + cy * cy;

  x = a2 * bc + b2 * ca + c2 * ab;

  if(x < 0)
    return -1;
  if(x > 0)
    return +1;

  return 0;
}

/*! \brief Tests whether point pp lies in the circumcircle around triangle
 *        p0,p1,p2 with some error margin.
 *
 *  This error margin should be large enough to exclude that close cases are
 *  misclssified due to numerical round-off errors.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] pp0 Index in DP of first point in triangle.
 *  \param[in] pp1 Index in DP of second point in triangle.
 *  \param[in] pp2 Index in DP of third point in triangle.
 *  \param[in] pp Index in DP of point to be checked.
 *
 *  \return (-1,0,1); -1: in circle; 0 on circle (within tolerance),
 *          1: outside circle.
 */
int InCircle_Errorbound(tessellation *T, int pp0, int pp1, int pp2, int pp)
{
  point *DP = T->DP;
  point *p0 = &DP[pp0];
  point *p1 = &DP[pp1];
  point *p2 = &DP[pp2];
  point *p  = &DP[pp];

  if(pp0 == DPinfinity || pp1 == DPinfinity || pp2 == DPinfinity || pp == DPinfinity)
    return -1;

  double ax, ay, bx, by, cx, cy;
  double ab, bc, ca, a2, b2, c2, x;
  double axby, bxay, bxcy, cxby, cxay, axcy;

  ax = p0->xx - p->xx;
  ay = p0->yy - p->yy;
  bx = p1->xx - p->xx;
  by = p1->yy - p->yy;
  cx = p2->xx - p->xx;
  cy = p2->yy - p->yy;

  axby = ax * by;
  bxay = bx * ay;
  bxcy = bx * cy;
  cxby = cx * by;
  cxay = cx * ay;
  axcy = ax * cy;

  ca = cxay - axcy;
  ab = axby - bxay;
  bc = bxcy - cxby;

  a2 = ax * ax + ay * ay;
  b2 = bx * bx + by * by;
  c2 = cx * cx + cy * cy;

  x = a2 * bc + b2 * ca + c2 * ab;

  /* calculate absolute maximum size */

  double sizelimit = a2 * (fabs(bxcy) + fabs(cxby)) + b2 * (fabs(cxay) + fabs(axcy)) + c2 * (fabs(axby) + fabs(bxay));

  double errbound = 1.0e-14 * sizelimit;

  if(x < -errbound)
    return -1;
  else if(x > errbound)
    return +1;

  return 0;
}

/*! \brief Tests whether point pp lies in the circumcircle around triangle
 *  p0,p1,p2 using arbitrary precision operations.
 *
 *  This is the exact solution, but computationally very expensive, thus only
 *  called for the unclear cases.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] pp0 Index in DP of first point in triangle.
 *  \param[in] pp1 Index in DP of second point in triangle.
 *  \param[in] pp2 Index in DP of third point in triangle.
 *  \param[in] pp Index in DP of point to be checked.
 *
 *  \return (-1,0,1); -1: in circle; 0 on circle,
 *          1: outside circle.
 */
int InCircle_Exact(tessellation *T, int pp0, int pp1, int pp2, int pp)
{
  point *DP = T->DP;
  point *p0 = &DP[pp0];
  point *p1 = &DP[pp1];
  point *p2 = &DP[pp2];
  point *p  = &DP[pp];

  if(pp0 == DPinfinity || pp1 == DPinfinity || pp2 == DPinfinity || pp == DPinfinity)
    return -1;

  IntegerMapType ax, ay, bx, by, cx, cy;

  ax = p0->ix - p->ix;
  ay = p0->iy - p->iy;
  bx = p1->ix - p->ix;
  by = p1->iy - p->iy;
  cx = p2->ix - p->ix;
  cy = p2->iy - p->iy;

  mpz_t axby, bxay, bxcy, cxby, cxay, axcy, tmp;

  mpz_init(tmp);

  mpz_init(axby);
  MY_mpz_set_si(tmp, ax);
  MY_mpz_mul_si(axby, tmp, by);
  mpz_init(bxay);
  MY_mpz_set_si(tmp, bx);
  MY_mpz_mul_si(bxay, tmp, ay);
  mpz_init(bxcy);
  MY_mpz_set_si(tmp, bx);
  MY_mpz_mul_si(bxcy, tmp, cy);
  mpz_init(cxby);
  MY_mpz_set_si(tmp, cx);
  MY_mpz_mul_si(cxby, tmp, by);
  mpz_init(cxay);
  MY_mpz_set_si(tmp, cx);
  MY_mpz_mul_si(cxay, tmp, ay);
  mpz_init(axcy);
  MY_mpz_set_si(tmp, ax);
  MY_mpz_mul_si(axcy, tmp, cy);

  mpz_t ca, ab, bc;

  mpz_init(ca);
  mpz_init(ab);
  mpz_init(bc);

  mpz_sub(ca, cxay, axcy);
  mpz_sub(ab, axby, bxay);
  mpz_sub(bc, bxcy, cxby);

  mpz_t AA, BB, a2, b2, c2;

  mpz_init(AA);
  mpz_init(BB);
  mpz_init(a2);
  mpz_init(b2);
  mpz_init(c2);

  MY_mpz_set_si(tmp, ax);
  MY_mpz_mul_si(AA, tmp, ax);
  MY_mpz_set_si(tmp, ay);
  MY_mpz_mul_si(BB, tmp, ay);
  mpz_add(a2, AA, BB);

  MY_mpz_set_si(tmp, bx);
  MY_mpz_mul_si(AA, tmp, bx);
  MY_mpz_set_si(tmp, by);
  MY_mpz_mul_si(BB, tmp, by);
  mpz_add(b2, AA, BB);

  MY_mpz_set_si(tmp, cx);
  MY_mpz_mul_si(AA, tmp, cx);
  MY_mpz_set_si(tmp, cy);
  MY_mpz_mul_si(BB, tmp, cy);
  mpz_add(c2, AA, BB);

  /* now calculate the final result */

  mpz_mul(AA, a2, bc);
  mpz_mul(BB, b2, ca);
  mpz_add(tmp, AA, BB);
  mpz_mul(BB, c2, ab);
  mpz_add(AA, BB, tmp);

  int sign = mpz_sgn(AA);

  mpz_clear(c2);
  mpz_clear(b2);
  mpz_clear(a2);
  mpz_clear(BB);
  mpz_clear(AA);
  mpz_clear(bc);
  mpz_clear(ab);
  mpz_clear(ca);
  mpz_clear(axcy);
  mpz_clear(cxay);
  mpz_clear(cxby);
  mpz_clear(bxcy);
  mpz_clear(bxay);
  mpz_clear(axby);
  mpz_clear(tmp);

  return sign;
}

/*! \brief Returns the orientation of the triangle.
 *
 *  Defined as the determinant of the matrix of the position of the three edge
 *  points a, b and c:
 *  | ax, ay, 1 |
 *  | bx, by, 1 |
 *  | cx, cy, 1 |
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] pp0 Index in DP of first point in triangle.
 *  \param[in] pp1 Index in DP of second point in triangle.
 *  \param[in] pp2 Index in DP of third point in triangle.
 *
 *  \return Determinant of orientation matrix.
 */
double test_triangle_orientation(tessellation *T, int pp0, int pp1, int pp2)
{
  point *DP = T->DP;
  point *p0 = &DP[pp0];
  point *p1 = &DP[pp1];
  point *p2 = &DP[pp2];

  return (p1->x - p0->x) * (p2->y - p0->y) - (p1->y - p0->y) * (p2->x - p0->x);
}

/*! \brief Check if triangle is positively or negatively oriented.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] pp0 Index in DP of first point in triangle.
 *  \param[in] pp1 Index in DP of second point in triangle.
 *  \param[in] pp2 Index in DP of third point in triangle.
 *
 *  \return -1 if negatively, 0 if degenerate (in a line) and 1 if positively
 *          oriented.
 */
int Orient2d_Quick(tessellation *T, int pp0, int pp1, int pp2)
{
  point *DP = T->DP;
  point *p0 = &DP[pp0];
  point *p1 = &DP[pp1];
  point *p2 = &DP[pp2];

  double x;

  x = (p1->xx - p0->xx) * (p2->yy - p0->yy) - (p1->yy - p0->yy) * (p2->xx - p0->xx);

  if(x < 0)
    return -1;
  if(x > 0)
    return +1;
  return 0;
}

/*! \brief Check if triangle is positively or negatively oriented.
 *
 *  Uses arbitrary precision operations, which is computationally expensive but
 *  garantees the correct result.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] pp0 Index in DP of first point in triangle.
 *  \param[in] pp1 Index in DP of second point in triangle.
 *  \param[in] pp2 Index in DP of third point in triangle.
 *
 *  \return -1 if negatively, 0 if degenerate (in a line) and 1 if positively
 *          oriented.
 */
int Orient2d_Exact(tessellation *T, int pp0, int pp1, int pp2)
{
  point *DP = T->DP;
  point *p0 = &DP[pp0];
  point *p1 = &DP[pp1];
  point *p2 = &DP[pp2];

#if USEDBITS > 31
  IntegerMapType dx1, dy1, dx2, dy2;

  dx1 = (p1->ix - p0->ix);
  dy1 = (p1->iy - p0->iy);
  dx2 = (p2->ix - p0->ix);
  dy2 = (p2->iy - p0->iy);

  mpz_t dx1dy2, dx2dy1, tmp;

  mpz_init(tmp);
  mpz_init(dx1dy2);
  mpz_init(dx2dy1);

  MY_mpz_set_si(tmp, dx1);
  MY_mpz_mul_si(dx1dy2, tmp, dy2);

  MY_mpz_set_si(tmp, dx2);
  MY_mpz_mul_si(dx2dy1, tmp, dy1);

  mpz_sub(tmp, dx1dy2, dx2dy1);

  int sign = mpz_sgn(tmp);

  mpz_clear(dx2dy1);
  mpz_clear(dx1dy2);
  mpz_clear(tmp);

  return (sign);

#else  /* #if USEDBITS > 31 */
  signed long long dx1, dy1, dx2, dy2, x;

  dx1 = (p1->ix - p0->ix);
  dy1 = (p1->iy - p0->iy);
  dx2 = (p2->ix - p0->ix);
  dy2 = (p2->iy - p0->iy);

  x = dx1 * dy2 - dy1 * dx2;

  if(x < 0)
    return -1;
  if(x > 0)
    return +1;
  return 0;
#endif /* #if USEDBITS > 31 #else */
}

const int edge_start[3] = {1, 2, 0};
const int edge_end[3]   = {2, 0, 1};

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
  int i, j, qq, p1, p2, k;
  face *f;
  double nx, ny;
  double sx, sy;
  double hx, hy;
  double dvol, h;

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

  point *dpi = &DP[t->p[i]];
  point *dpj = &DP[t->p[j]];

  qq = t->t[nr];

  Edge_visited[tt] |= (1 << nr);
  Edge_visited[qq] |= (1 << (t->s[nr]));

  p1 = t->p[i];
  p2 = t->p[j];

  f = &VF[T->Nvf++];

  f->p1 = p1;
  f->p2 = p2;

  f->cx = 0.5 * (DTC[tt].cx + DTC[qq].cx);
  f->cy = 0.5 * (DTC[tt].cy + DTC[qq].cy);
  f->cz = 0;

#ifdef TETRA_INDEX_IN_FACE
  f->dt_index = tt;
#endif /* #ifdef TETRA_INDEX_IN_FACE */

#ifdef REFINEMENT_MERGE_CELLS
  f->t  = tt;
  f->nr = nr; /* delaunay tetra and edge number that generated this face */
#endif        /* #ifdef REFINEMENT_MERGE_CELLS */

  nx = DTC[tt].cx - DTC[qq].cx;
  ny = DTC[tt].cy - DTC[qq].cy;

  f->area = sqrt(nx * nx + ny * ny);

  hx = 0.5 * (dpi->x - dpj->x);
  hy = 0.5 * (dpi->y - dpj->y);

  h    = sqrt(hx * hx + hy * hy);
  dvol = 0.5 * f->area * h;

#if defined(REGULARIZE_MESH_FACE_ANGLE) || defined(OUTPUT_MESH_FACE_ANGLE)
  double angle = 0.5 * f->area / h;
#endif /* #if defined(REGULARIZE_MESH_FACE_ANGLE) || defined(OUTPUT_MESH_FACE_ANGLE) */

  if(dpi->task == ThisTask && dpi->index >= 0 && dpi->index < NumGas)
    {
      if(TimeBinSynchronized[P[dpi->index].TimeBinHydro])
        {
          SphP[dpi->index].Volume += dvol;
          SphP[dpi->index].SurfaceArea += f->area;

#if defined(REGULARIZE_MESH_FACE_ANGLE) || defined(OUTPUT_MESH_FACE_ANGLE)
          if(SphP[dpi->index].MaxFaceAngle < angle)
            SphP[dpi->index].MaxFaceAngle = angle;
#endif /* #if defined(REGULARIZE_MESH_FACE_ANGLE) || defined(OUTPUT_MESH_FACE_ANGLE) */

#ifdef OUTPUT_SURFACE_AREA
          if(f->area)
            SphP[dpi->index].CountFaces++;
#endif /* #ifdef OUTPUT_SURFACE_AREA */

#if defined(REFINEMENT_SPLIT_CELLS)
          if(SphP[dpi->index].MinimumEdgeDistance > h)
            SphP[dpi->index].MinimumEdgeDistance = h;
#endif /* #if defined(REFINEMENT_SPLIT_CELLS) */
          /* let's now compute the center-of-mass of the pyramid at the bottom top */
          sx = (2.0 / 3) * f->cx + (1.0 / 3) * dpi->x;
          sy = (2.0 / 3) * f->cy + (1.0 / 3) * dpi->y;

          SphP[dpi->index].Center[0] += dvol * sx;
          SphP[dpi->index].Center[1] += dvol * sy;
        }
    }

  if(dpj->task == ThisTask && dpj->index >= 0 && dpj->index < NumGas)
    {
      if(TimeBinSynchronized[P[dpj->index].TimeBinHydro])
        {
          SphP[dpj->index].Volume += dvol;
          SphP[dpj->index].SurfaceArea += f->area;

#if defined(REGULARIZE_MESH_FACE_ANGLE) || defined(OUTPUT_MESH_FACE_ANGLE)
          if(SphP[dpj->index].MaxFaceAngle < angle)
            SphP[dpj->index].MaxFaceAngle = angle;
#endif /* #if defined(REGULARIZE_MESH_FACE_ANGLE) || defined(OUTPUT_MESH_FACE_ANGLE) */

#ifdef OUTPUT_SURFACE_AREA
          if(f->area)
            SphP[dpj->index].CountFaces++;
#endif /* #ifdef OUTPUT_SURFACE_AREA */

#if defined(REFINEMENT_SPLIT_CELLS)
          if(SphP[dpj->index].MinimumEdgeDistance > h)
            SphP[dpj->index].MinimumEdgeDistance = h;
#endif /* #if defined(REFINEMENT_SPLIT_CELLS) */

          /* let's now compute the center-of-mass of the pyramid on top */
          sx = (2.0 / 3) * f->cx + (1.0 / 3) * dpj->x;
          sy = (2.0 / 3) * f->cy + (1.0 / 3) * dpj->y;

          SphP[dpj->index].Center[0] += dvol * sx;
          SphP[dpj->index].Center[1] += dvol * sy;
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

/*! \brief Copies triangle information from DTC array to trilist.
 *
 *  Performs an orientation check and swaps orientation if needed.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in] tt Index of triangle in DT array.
 *  \param[in] nr Index in DT[tt].t array (adjacent tetrahedrons).
 *  \param[in] dtip Pointer to point to be inserted.
 *  \param[out] trilist Array of triangles.
 *  \param[in] ntri Index in trilist array.
 *  \param[in] max_n_tri Maximum index in trilist array.
 *
 *  \return Next index in trilist array.
 */
int derefine_refine_get_triangles(tessellation *T, int tt, int nr, point *dtip, triangle *trilist, int ntri, int max_n_tri)
{
  tetra *DT         = T->DT;
  tetra_center *DTC = T->DTC;
  tetra *t          = &DT[tt];
  int qq            = t->t[nr];

  if(ntri >= max_n_tri)
    terminate("ntri >= max_n_tri");

  trilist[ntri].p[0][0] = DTC[tt].cx;
  trilist[ntri].p[0][1] = DTC[tt].cy;

  trilist[ntri].p[1][0] = DTC[qq].cx;
  trilist[ntri].p[1][1] = DTC[qq].cy;

  trilist[ntri].p[2][0] = dtip->x;
  trilist[ntri].p[2][1] = dtip->y;

  if(get_tri_volume(ntri, trilist) < 0)
    {
      /* swap two points to get proper orientation */
      trilist[ntri].p[1][0] = DTC[tt].cx;
      trilist[ntri].p[1][1] = DTC[tt].cy;

      trilist[ntri].p[0][0] = DTC[qq].cx;
      trilist[ntri].p[0][1] = DTC[qq].cy;
    }

  ntri++;

  return ntri;
}

/*! \brief Add point and adjust triangles accordingly.
 *
 *  \param[in] q Index of point in DP array.
 *  \param[in, out] trilist Array of triangles.
 *  \param[in] ntri Number of elements in trilist before splitting.
 *  \param[in] max_ntri Maximum number of triangles allowed.
 *  \param[in] vol (Unused)
 *
 *  \return Updated number of triangles.
 */
int derefine_add_point_and_split_tri(int q, triangle *trilist, int ntri, int max_ntri, double vol)
{
  double m[2], n[2], sc[3], *a;
  double cut[2][2], ed[2];
  int i, j, k, kk, l, nnew, flag[3], count, oldq;

  for(i = 0, nnew = ntri; i < ntri; i++)
    {
      if(trilist[i].owner < 0 || trilist[i].owner >= Mesh.Ndp)
        {
          char buf[1000];
          sprintf(buf, "i=%d trilist[i].owner=%d\n", i, trilist[i].owner);
          terminate(buf);
        }

      if(q < 0 || q >= Mesh.Ndp)
        {
          char buf[1000];
          sprintf(buf, "i=%d q=%d\n", i, q);
          terminate(buf);
        }

      /* midpoint */
      m[0] = 0.5 * (Mesh.DP[q].x + Mesh.DP[trilist[i].owner].x);
      m[1] = 0.5 * (Mesh.DP[q].y + Mesh.DP[trilist[i].owner].y);

      n[0] = (Mesh.DP[q].x - Mesh.DP[trilist[i].owner].x);
      n[1] = (Mesh.DP[q].y - Mesh.DP[trilist[i].owner].y);

      if(q == trilist[i].owner)
        terminate("q == trilist[i].owner");

      for(k = 0, count = 0; k < 3; k++) /* determine the side of each point */
        {
          a = &trilist[i].p[k][0];

          sc[k] = (a[0] - m[0]) * n[0] + (a[1] - m[1]) * n[1];

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

          case 3:                 /* the whole tetra is on the side of new point */
            trilist[i].owner = q; /* change owner */
            break;

          case 1:
          case 2:

            if(nnew + 2 > max_ntri)
              terminate("nnew + 2 > max_ntri");

            trilist[nnew]     = trilist[i];
            trilist[nnew + 1] = trilist[i];

            /* find the point index that is on the other side */
            for(k = 0; k < 3; k++)
              {
                if(flag[k] == 1 && count == 1)
                  break;
                if(flag[k] == 0 && count == 2)
                  break;
              }

            for(j = 0; j < 2; j++)
              {
                kk = k + j + 1;
                if(kk > 2)
                  kk -= 3;

                double *b = trilist[i].p[k];
                double *a = trilist[i].p[kk];

                for(l = 0; l < 2; l++)
                  ed[l] = a[l] - b[l];

                double prod = (ed[0] * n[0] + ed[1] * n[1]);
                double t;
                if(prod)
                  t = -sc[k] / prod;
                else
                  t = 0.5;

                if(t < 0)
                  t = 0;
                if(t > 1)
                  t = 1;

                for(l = 0; l < 2; l++)
                  cut[j][l] = b[l] + t * ed[l];
              }

            /* modify the tetra that's assigned to the new point */
            for(j = 0; j < 2; j++)
              {
                kk = k + j + 1;
                if(kk > 2)
                  kk -= 3;

                for(l = 0; l < 2; l++)
                  trilist[i].p[kk][l] = cut[j][l];
              }

            oldq = trilist[i].owner;

            if(count == 1)
              trilist[i].owner = q;

            /* modify the two new tetras */
            kk = k + 1;
            if(kk > 2)
              kk -= 3;

            for(l = 0; l < 2; l++)
              {
                trilist[nnew].p[k][l] = cut[0][l];

                trilist[nnew + 1].p[k][l]  = cut[1][l];
                trilist[nnew + 1].p[kk][l] = cut[0][l];
              }

            if(count == 1)
              {
                trilist[nnew].owner     = oldq;
                trilist[nnew + 1].owner = oldq;
              }
            else
              {
                trilist[nnew].owner     = q;
                trilist[nnew + 1].owner = q;
              }
            nnew += 2;
            break;
        }
    }

  return nnew;
}

/*! \brief Determines area of triangle (i.e. 2d Volume).
 *
 *  \param i Index in trilist array.
 *  \param trilist Array with triangles.
 *
 *  \return Area of triangle.
 */
double get_tri_volume(int i, triangle *trilist)
{
  double *p0 = &trilist[i].p[0][0];
  double *p1 = &trilist[i].p[1][0];
  double *p2 = &trilist[i].p[2][0];

  double nz = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0]);

  return 0.5 * nz;
}

/*! \brief Process edge for volume calculation.
 *
 *  Calculates the contribution of edge to volumes of neighboring
 *  Voronoi cells in vol array.
 *
 *  \param[in] T Pointer to tessellation.
 *  \param[in, out] vol Volume of tetrahedra.
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

  int i, j, qq, p1, p2;
  double nx, ny;
  double hx, hy;
  double dvol, h;

  tetra *t = &DT[tt];

  i = edge_start[nr];
  j = edge_end[nr];

  point *dpi = &DP[t->p[i]];
  point *dpj = &DP[t->p[j]];

  qq = t->t[nr];

  Edge_visited[tt] |= (1 << nr);
  Edge_visited[qq] |= (1 << (t->s[nr]));

  p1 = t->p[i];
  p2 = t->p[j];

  nx = DTC[tt].cx - DTC[qq].cx;
  ny = DTC[tt].cy - DTC[qq].cy;

  double area = sqrt(nx * nx + ny * ny);

  hx = 0.5 * (dpi->x - dpj->x);
  hy = 0.5 * (dpi->y - dpj->y);

  h    = sqrt(hx * hx + hy * hy);
  dvol = 0.5 * area * h;

  if(p1 >= 0 && p1 < DeRefMesh.Ndp)
    vol[p1] += dvol;

  if(p2 >= 0 && p2 < DeRefMesh.Ndp)
    vol[p2] += dvol;
}

/*! \brief Computes the circum-circle of all triangles in mesh.
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

      if(DT[i].p[0] == DPinfinity)
        continue;
      if(DT[i].p[1] == DPinfinity)
        continue;
      if(DT[i].p[2] == DPinfinity)
        continue;

      update_circumcircle(T, i);
    }
}

/*! \brief Computes the circum-circle of triangle tt.
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

  tetra *t = &DT[tt];
  point *p0, *p1, *p2;
  int pp0, pp1, pp2;

  pp0 = t->p[0];
  pp1 = t->p[1];
  pp2 = t->p[2];

  p0 = &DP[pp0];
  p1 = &DP[pp1];
  p2 = &DP[pp2];

  if(t->p[0] == DPinfinity)
    return;
  if(t->p[1] == DPinfinity)
    return;
  if(t->p[2] == DPinfinity)
    return;

  double ax = p1->xx - p0->xx;
  double ay = p1->yy - p0->yy;

  double bx = p2->xx - p0->xx;
  double by = p2->yy - p0->yy;

  double aa = 0.5 * (ax * ax + ay * ay);
  double bb = 0.5 * (bx * bx + by * by);

  double mv_data[] = {ax, ay, aa, bx, by, bb};
  double x[2];

  int status = solve_linear_equations_2d(mv_data, x);

  if(status < 0)
    {
      terminate("trouble in circum-circle calculation\n");
    }
  else
    {
      x[0] += p0->xx;
      x[1] += p0->yy;

      DTC[tt].cx = (x[0] - 1.0) / ConversionFac + CentralOffsetX;
      DTC[tt].cy = (x[1] - 1.0) / ConversionFac + CentralOffsetY;
      DTC[tt].cz = 0;
    }
}

/*! \brief Computes the integer coordinates from coordinates for a point.
 *
 *  \pararm[in, out] p Pointer to point.
 *
 *  \return void
 */
void set_integers_for_pointer(point *p)
{
  p->xx = (p->x - CentralOffsetX) * ConversionFac + 1.0;
  p->yy = (p->y - CentralOffsetY) * ConversionFac + 1.0;

  if(p->xx < 1.0 || p->xx >= 2.0 || p->yy < 1.0 || p->yy >= 2.0)
    {
      printf("(%g, %g) (%g, %g)\n", p->x, p->y, p->xx, p->yy);
      terminate("invalid coordinate range");
    }

  p->ix = double_to_voronoiint(p->xx);
  p->iy = double_to_voronoiint(p->yy);

  p->xx = mask_voronoi_int(p->xx);
  p->yy = mask_voronoi_int(p->yy);
}

/*! \brief Outputs Voronoi mesh to file.
 *
 *  Outputs the Voronoi mesh data from task write Task to lastTask in file
 *  fname.
 *
 *  \param[in] T Pointer to tesselation.
 *  \param[in] fname File name of file the data is written in.
 *  \param[in] writeTask Task that gathers information and writes data.
 *  \param[in] lastTask Last task that is included in this dump.
 *
 *  \return void
 */
void write_voronoi_mesh(tessellation *T, char *fname, int writeTask, int lastTask)
{
  CPU_Step[CPU_MISC] += measure_time();

  FILE *fd;
  char msg[1000];
  MPI_Status status;
  int i, j, k, MaxNel, Nel;
  int ngas_tot, nel_tot, ndt_tot, nel_before, ndt_before, task;
  int *EdgeList, *Nedges, *NedgesOffset, *whichtetra;
  int *ngas_list, *nel_list, *ndt_list, *tmp;
  float *xyz_edges;
  tetra *q, *qstart;

  tetra_center *DTC = T->DTC;
  tetra *DT         = T->DT;
  point *DP         = T->DP;

  MaxNel = 10 * NumGas; /* max edge list */
  Nel    = 0;           /* length of edge list */

  EdgeList     = mymalloc("EdgeList", MaxNel * sizeof(int));
  Nedges       = mymalloc("Nedges", NumGas * sizeof(int));
  NedgesOffset = mymalloc("NedgesOffset", NumGas * sizeof(int));
  whichtetra   = mymalloc("whichtetra", NumGas * sizeof(int));
  xyz_edges    = mymalloc("xyz_edges", T->Ndt * DIMS * sizeof(float));
  ngas_list    = mymalloc("ngas_list", sizeof(int) * NTask);
  nel_list     = mymalloc("nel_list", sizeof(int) * NTask);
  ndt_list     = mymalloc("ndt_list", sizeof(int) * NTask);

  for(i = 0; i < T->Ndt; i++)
    {
      xyz_edges[i * DIMS + 0] = DTC[i].cx;
      xyz_edges[i * DIMS + 1] = DTC[i].cy;
    }

  for(i = 0; i < NumGas; i++)
    {
      Nedges[i]     = 0;
      whichtetra[i] = -1;
    }

  for(i = 0; i < T->Ndt; i++)
    {
      for(j = 0; j < DIMS + 1; j++)
        if(DP[DT[i].p[j]].task == ThisTask && DP[DT[i].p[j]].index >= 0 && DP[DT[i].p[j]].index < NumGas)
          whichtetra[DP[DT[i].p[j]].index] = i;
    }

  for(i = 0; i < NumGas; i++)
    {
      if(whichtetra[i] < 0)
        continue;

      qstart = q = &DT[whichtetra[i]];

      do
        {
          Nedges[i]++;

          if(Nel >= MaxNel)
            terminate("Nel >= MaxNel");

          EdgeList[Nel++] = q - DT;

          for(j = 0; j < 3; j++)
            if(DP[q->p[j]].task == ThisTask && DP[q->p[j]].index == i)
              break;

          k = j + 1;
          if(k >= 3)
            k -= 3;

          q = &DT[q->t[k]];
        }
      while(q != qstart);
    }

  for(i = 1, NedgesOffset[0] = 0; i < NumGas; i++)
    NedgesOffset[i] = NedgesOffset[i - 1] + Nedges[i - 1];

  /* determine particle numbers and number of edges in file */

  if(ThisTask == writeTask)
    {
      ngas_tot = NumGas;
      nel_tot  = Nel;
      ndt_tot  = T->Ndt;

      for(task = writeTask + 1; task <= lastTask; task++)
        {
          MPI_Recv(&ngas_list[task], 1, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
          MPI_Recv(&nel_list[task], 1, MPI_INT, task, TAG_LOCALN + 1, MPI_COMM_WORLD, &status);
          MPI_Recv(&ndt_list[task], 1, MPI_INT, task, TAG_LOCALN + 2, MPI_COMM_WORLD, &status);

          MPI_Send(&nel_tot, 1, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
          MPI_Send(&ndt_tot, 1, MPI_INT, task, TAG_N + 1, MPI_COMM_WORLD);

          ngas_tot += ngas_list[task];
          nel_tot += nel_list[task];
          ndt_tot += ndt_list[task];
        }

      if(!(fd = fopen(fname, "w")))
        {
          sprintf(msg, "can't open file `%s' for writing snapshot.\n", fname);
          terminate(msg);
        }

      my_fwrite(&ngas_tot, sizeof(int), 1, fd);
      my_fwrite(&nel_tot, sizeof(int), 1, fd);
      my_fwrite(&ndt_tot, sizeof(int), 1, fd);

      my_fwrite(Nedges, sizeof(int), NumGas, fd);
      for(task = writeTask + 1; task <= lastTask; task++)
        {
          tmp = mymalloc("tmp", sizeof(int) * ngas_list[task]);
          MPI_Recv(tmp, ngas_list[task], MPI_INT, task, TAG_N + 2, MPI_COMM_WORLD, &status);
          my_fwrite(tmp, sizeof(int), ngas_list[task], fd);
          myfree(tmp);
        }

      my_fwrite(NedgesOffset, sizeof(int), NumGas, fd);
      for(task = writeTask + 1; task <= lastTask; task++)
        {
          tmp = mymalloc("tmp", sizeof(int) * ngas_list[task]);
          MPI_Recv(tmp, ngas_list[task], MPI_INT, task, TAG_N + 3, MPI_COMM_WORLD, &status);
          my_fwrite(tmp, sizeof(int), ngas_list[task], fd);
          myfree(tmp);
        }

      my_fwrite(EdgeList, sizeof(int), Nel, fd);
      for(task = writeTask + 1; task <= lastTask; task++)
        {
          tmp = mymalloc("tmp", sizeof(int) * nel_list[task]);
          MPI_Recv(tmp, nel_list[task], MPI_INT, task, TAG_N + 4, MPI_COMM_WORLD, &status);
          my_fwrite(tmp, sizeof(int), nel_list[task], fd);
          myfree(tmp);
        }

      my_fwrite(xyz_edges, sizeof(float), T->Ndt * DIMS, fd);
      for(task = writeTask + 1; task <= lastTask; task++)
        {
          tmp = mymalloc("tmp", sizeof(float) * DIMS * ndt_list[task]);
          MPI_Recv(tmp, sizeof(float) * DIMS * ndt_list[task], MPI_BYTE, task, TAG_N + 5, MPI_COMM_WORLD, &status);
          my_fwrite(tmp, sizeof(float), DIMS * ndt_list[task], fd);
          myfree(tmp);
        }

      fclose(fd);
    }
  else
    {
      MPI_Send(&NumGas, 1, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
      MPI_Send(&Nel, 1, MPI_INT, writeTask, TAG_LOCALN + 1, MPI_COMM_WORLD);
      MPI_Send(&T->Ndt, 1, MPI_INT, writeTask, TAG_LOCALN + 2, MPI_COMM_WORLD);

      MPI_Recv(&nel_before, 1, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
      MPI_Recv(&ndt_before, 1, MPI_INT, writeTask, TAG_N + 1, MPI_COMM_WORLD, &status);

      for(i = 0; i < NumGas; i++)
        NedgesOffset[i] += nel_before;
      for(i = 0; i < Nel; i++)
        EdgeList[i] += ndt_before;

      MPI_Send(Nedges, NumGas, MPI_INT, writeTask, TAG_N + 2, MPI_COMM_WORLD);
      MPI_Send(NedgesOffset, NumGas, MPI_INT, writeTask, TAG_N + 3, MPI_COMM_WORLD);
      MPI_Send(EdgeList, Nel, MPI_INT, writeTask, TAG_N + 4, MPI_COMM_WORLD);
      MPI_Send(xyz_edges, sizeof(float) * DIMS * T->Ndt, MPI_BYTE, writeTask, TAG_N + 5, MPI_COMM_WORLD);
    }

  myfree(ndt_list);
  myfree(nel_list);
  myfree(ngas_list);
  myfree(xyz_edges);
  myfree(whichtetra);
  myfree(NedgesOffset);
  myfree(Nedges);
  myfree(EdgeList);

  mpi_printf("wrote Voronoi mesh to file\n");

  CPU_Step[CPU_MAKEIMAGES] += measure_time();
}

#endif /* #if defined(TWODIMS) && !defined(ONEDIMS) */
