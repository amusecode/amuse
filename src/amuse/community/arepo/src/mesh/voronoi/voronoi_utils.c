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
 * \file        src/mesh/voronoi/voronoi_utils.c
 * \date        05/2018
 * \brief       Utilities for 3d Voronoi mesh
 * \details     contains functions:
 *                double cross_section_plane_cell(int sphp_index, int dp_index, double *center, double *n)
 *                void intersections_plane_cell(int sphp_index, int dp_index, double *center, double *n, double *polygon, unsigned int
 * *nof_polygon_elements) void intersection_plane_grid(double *center, double *n, const char *filename) static double
 * polygon_area(double *polygon, unsigned int nof_elements) static int qs_partition(double *A, int p, int r, double *B) static void
 * qs_sort(double *A, int p, int r, double *B) static double calc_phi(double x, double y) static void rotate_z(double *vec, const
 * double alpha)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 23.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include "../../main/allvars.h"
#include "../../main/proto.h"

#if !defined(TWODIMS) && !defined(ONEDIMS)

// helper functions for cross_section_plane_cell and intersections_plane_cell:
static int qs_partition(double *A, int p, int r, double *B);
static void qs_sort(double *A, int p, int r, double *B);
static double calc_phi(double x, double y);
static void rotate_z(double *vec, const double alpha);
static void rotate_y(double *vec, const double alpha);

#ifdef TETRA_INDEX_IN_FACE
static double polygon_area(double *polygon, unsigned int nof_elements);
static const unsigned int max_poly_elements = 128;

/*! \brief Calculates the cross section between a plane and a Voronoi cell(3D).
 *
 *  \param[in] sphp_index The hydro index of the cell.
 *  \param[in] dp_index The delaunay point index of the cell.
 *  \param[in] center A point in the plane.
 *  \param[in] n A vector starting at center and normal to the plane.
 *
 *  \return The cross section between the plane and the cell.
 */
double cross_section_plane_cell(int sphp_index, int dp_index, double *center, double *n)
{
  double polygon[max_poly_elements];
  unsigned int nof_elements = 0;

  intersections_plane_cell(sphp_index, dp_index, center, n, polygon, &nof_elements);

  // polygon has to contain at least 3 points
  if(nof_elements < 6)
    {
      return 0;
    }
  else
    {
      return polygon_area(polygon, nof_elements);
    }
}

/*! \brief Calculates the intersections between a plane and a cell.
 *
 *  \param[in] sphp_index The hydro index of the cell.
 *  \param[in] dp_index The Delaunay point index of the cell.
 *  \param[in] center A point in the plane.
 *  \param[in] n A vector starting at center and normal to the plane.
 *  \param[out] polygon Store the intersections (polygon) in this array.
 *  \param[out] nof_polygon_elements The number of stored elements in the
 *              polygon array.
 *
 *  \return void
 */
void intersections_plane_cell(int sphp_index, int dp_index, double *center, double *n, double *polygon,
                              unsigned int *nof_polygon_elements)
{
  // memory for the line segments
  unsigned int line_segments_max = 2000;
  double *ls                     = (double *)mymalloc("line_segments", line_segments_max * sizeof(double));

  // get the line segments
  unsigned int nof_elements = 0;
  get_line_segments(sphp_index, dp_index, ls, &nof_elements, line_segments_max);
  assert(nof_elements % 6 == 0);  // 6 doubles represent one line segment

  // start the calculation
  unsigned int i;
  double phi;

  if(n[0] == 0 && n[1] == 0)
    {
      phi = 0;
    }
  else
    {
      phi = calc_phi(n[0], n[1]);
    }

  double r = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  assert(r > 0);
  double theta = acos(n[2] / r);

  double lambda;  // z1 + lambda * (z2 - z1) = 0

  unsigned int max_phi_elms = max_poly_elements / 2;
  double phi_values[max_phi_elms];  // phi coordinates of the points of the polygon
  unsigned int p = 0;               // number of points of the polygon

  // balance point of the polygon
  double bal_p_x = 0;
  double bal_p_y = 0;

  for(i = 0; i < nof_elements; i += 6)
    {
      // transform line segment to the center frame
      ls[i] -= center[0];      // x1
      ls[i + 1] -= center[1];  // y1
      ls[i + 2] -= center[2];  // z1
      ls[i + 3] -= center[0];  // x2
      ls[i + 4] -= center[1];  // y2
      ls[i + 5] -= center[2];  // z2

      // rotate line segment such that the cross secting plane is in the x-y plane / the normal vector of the plane is on the z-axis
      rotate_z(&ls[i], -phi);
      rotate_y(&ls[i], -theta);

      rotate_z(&ls[i + 3], -phi);
      rotate_y(&ls[i + 3], -theta);

      if(ls[i + 2] == ls[i + 5])  // same z-coords
        {
          if(ls[i + 2] != 0)  // no intersection
            {
              lambda = -1;
            }
          else
            {
              lambda = 0;  // take first point as intersection
            }
        }
      else
        {
          lambda = ls[i + 2] / (ls[i + 2] - ls[i + 5]);
        }

      if(lambda >= 0 && lambda <= 1)  // line segment intersects plane
        {
          if(p == max_phi_elms)
            {
              terminate("termination in voronoi_utils.c: intersections_plane_cell: not enough memory!\n");
            }

          polygon[2 * p]     = ls[i] + lambda * (ls[i + 3] - ls[i]);          // x coordinate of the intersection
          polygon[2 * p + 1] = ls[i + 1] + lambda * (ls[i + 4] - ls[i + 1]);  // y coordinate of the intersection

          bal_p_x += polygon[2 * p];
          bal_p_y += polygon[2 * p + 1];

          p++;
        }
    }

  // free memory
  myfree(ls);

  // polygon has to contain at least 3 points
  if(p < 3)
    {
      return;
    }

  // switch frame to balance point of the polygon
  bal_p_x /= p;
  bal_p_y /= p;

  for(i = 0; i < p; i++)
    {
      polygon[2 * i] -= bal_p_x;
      polygon[2 * i + 1] -= bal_p_y;

      // calculate the phi values
      phi_values[i] = calc_phi(polygon[2 * i], polygon[2 * i + 1]);
    }

  // sort polygon
  qs_sort(phi_values, 0, p - 1, polygon);

  // close polygon
  polygon[2 * p]     = polygon[0];
  polygon[2 * p + 1] = polygon[1];
  phi_values[p]      = phi_values[0];
  p++;

  // transform back
  for(i = 0; i < p; i++)
    {
      polygon[2 * i] += bal_p_x;
      polygon[2 * i + 1] += bal_p_y;
    }

  *nof_polygon_elements = 2 * p;
}

/*! \brief Write out the intersections between a plane and the grid
 *         (for plotting).
 *
 *  Binary output:
 *  int: Number of elements in the first array.
 *  int: Number of elements in the second array.
 *  int[]: Array, which stores the number of intersections for each intersected
 *         cell.
 *         The j-th entry gives the number of elements in the intersections
 *         array which correspond to the j-th intersected cell.
 *  double[]: intersections array, all intersections are stored in the
 *            order x1,y1,x2,y2,x3,y3,...
 *
 *  The intersections are given in a coordinate system where n is the z-axis
 *  and which has its origin at center.
 *
 *  \param[in] center A point in the plane.
 *  \param[in] n A vector starting at center and normal to the plane.
 *  \param[in] filename Filename.
 *
 *  \return void
 */
void intersection_plane_grid(double *center, double *n, const char *filename)
{
  if(NTask != 1)
    {
      terminate("termination in voronoi_utils.c: intersection_plane_grid: not yet parallelized!\n");
    }

  double phi;

  if(n[0] == 0 && n[1] == 0)
    {
      phi = 0;
    }
  else
    {
      phi = calc_phi(n[0], n[1]);
    }

  double r = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  assert(r > 0);
  double theta = acos(n[2] / r);

  double xaxis[3] = {1, 0, 0};
  double yaxis[3] = {0, 1, 0};
  double zaxis[3] = {0, 0, 1};

  rotate_y(xaxis, theta);
  rotate_z(xaxis, phi);

  rotate_y(yaxis, theta);
  rotate_z(yaxis, phi);

  rotate_y(zaxis, theta);
  rotate_z(zaxis, phi);

  printf("normal vector: (%f, %f, %f)\n", n[0], n[1], n[2]);
  printf("Coordinate system of output data: \n");
  printf("center: (%f, %f, %f)\n", center[0], center[1], center[2]);
  printf("x-axis: (%f, %f, %f)\n", xaxis[0], xaxis[1], xaxis[2]);
  printf("y-axis: (%f, %f, %f)\n", yaxis[0], yaxis[1], yaxis[2]);
  printf("z-axis: (%f, %f, %f)\n", zaxis[0], zaxis[1], zaxis[2]);

  const int cells_max_elms = NumGas;
  int *nof_intersections   = (int *)mymalloc("number of intersections", cells_max_elms * sizeof(int));
  unsigned int l           = 0;

  const int polygons_max_elms = NumGas * 5;
  double *polygons            = (double *)mymalloc("polygons", polygons_max_elms * 5 * sizeof(int));
  unsigned int j              = 0;

  unsigned int nof_polygon_elements = 0;

  unsigned int k = 0;

  for(k = 0; k < NumGas; k++)
    {
      nof_polygon_elements = 0;
      intersections_plane_cell(k, k, center, n, &polygons[j], &nof_polygon_elements);

      if(nof_polygon_elements != 0)
        {
          nof_intersections[l] = (int)nof_polygon_elements;
          l++;

          j += nof_polygon_elements;

          if(j > polygons_max_elms - 100)
            {
              terminate("termination in voronoi_utils.c: intersection_plane_grid: not enough memory for the polygons!\n");
            }
        }
    }

  // binary output
  FILE *pFile;

  pFile = fopen(filename, "wb");

  fwrite(&l, sizeof(int), 1, pFile);  // number of intersected cells
  fwrite(&j, sizeof(int), 1, pFile);  // number of elements in polygons array
  fwrite(nof_intersections, sizeof(int), l, pFile);
  fwrite(polygons, sizeof(double), j, pFile);

  fclose(pFile);

  myfree(polygons);
  myfree(nof_intersections);
}

/*! \brief Calculate the area of a 2D polygon.
 *
 *  Formula (wikipedia):A = 0.5 * sum_i=0^{n-1}(x_i * y_{i+1} - x_{i+1} * y_i).
 *
 *  \param[in] polygon Array of points of the polygon: x1, y1, x2, y2, ...,
 *             has to be sorted counterclockwise and closed
 *             (x_n == x_0 && y_n == y_0).
 *  \param[in] nof_elements Number of elements in the array.
 *
 *  \return Area of polygon.
 */
static double polygon_area(double *polygon, unsigned int nof_elements)
{
  assert(nof_elements >= 8);

  double result = 0;

  unsigned int k;

  for(k = 0; k < nof_elements - 2; k += 2)
    {
      result += polygon[k] * polygon[k + 3] - polygon[k + 2] * polygon[k + 1];
    }

  result *= 0.5;

  assert(result >= 0);

  return result;
}

#endif /* #ifdef TETRA_INDEX_IN_FACE */

/*! \brief Quicksort partitioning function, helper for qs_sort.
 *
 *  \param[in, out] A array to be sorted, usually angle phi.
 *  \param[in] p Lower index for quicksort.
 *  \param[in] r Upper index for quicksort.
 *  \param[in, out] B Array that also changes ordering the same way as A.
 *
 *  \return Index for partitioning.
 */
static int qs_partition(double *A, int p, int r, double *B)
{
  double x = A[r];
  double tmp;
  double tmp2;
  int i = p - 1;
  int j;

  for(j = p; j < r; j++)
    {
      if(A[j] <= x)
        {
          // switch phi values ( i <-> j )
          i++;
          tmp  = A[i];
          A[i] = A[j];
          A[j] = tmp;

          // switch coordinates ( 2i, 2i+1 <-> 2j, 2j+1)
          tmp          = B[2 * i];
          tmp2         = B[2 * i + 1];
          B[2 * i]     = B[2 * j];
          B[2 * i + 1] = B[2 * j + 1];
          B[2 * j]     = tmp;
          B[2 * j + 1] = tmp2;
        }
    }

  // switch phi values
  tmp      = A[i + 1];
  A[i + 1] = A[r];
  A[r]     = tmp;

  // switch coordinates
  tmp  = B[(i + 1) * 2];
  tmp2 = B[(i + 1) * 2 + 1];

  B[(i + 1) * 2]     = B[2 * r];
  B[(i + 1) * 2 + 1] = B[2 * r + 1];

  B[2 * r]     = tmp;
  B[2 * r + 1] = tmp2;

  return i + 1;
}

/*! \brief Quick-sorts the points of the polygon with respect to phi.
 *
 *  \param[in, out] A array to be sorted, usually angle phi.
 *  \param[in] p lower index for quicksort.
 *  \param[in] r upper index for quicksort.
 *  \param[in, out] B array that also changes ordering the same way as A;
 *                  usually polygon.
 *
 *  \return void
 */
static void qs_sort(double *A, int p, int r, double *B)
{
  int q;

  if(p < r)
    {
      q = qs_partition(A, p, r, B);
      qs_sort(A, p, q - 1, B);
      qs_sort(A, q + 1, r, B);
    }
}

/*! \brief Calculates the phi coordinate of a point.
 *
 *  Calculates polar angle in a 2d coordinate system from Cartesian coordinate
 *  system.
 *
 *  \param[in] x X coordinate.
 *  \param[in] y Y coordinate.
 *
 *  \return Phi (polar angle).
 */
static double calc_phi(double x, double y)
{
  // if both arguments are zero an error occurs in atan2
  if((x == 0) && (y == 0))
    {
      fprintf(stderr, "ERROR in calc_phi: both arguments are zero\n");
      return 0;
    }

  double p = atan2(y, x);  // in [-pi,pi]

  if(p < 0)
    {
      return p + 2 * M_PI;
    }

  return p;
}

/*! \brief Rotate a vector around the z axis.
 *
 *  \param[in, out] vec Array to 3 dimensional vector to be rotated.
 *  \param[in] alpha Rotation angle.
 *
 *  \return void
 */
static void rotate_z(double *vec, const double alpha)
{
  double vx_tmp = vec[0];
  vec[0]        = cos(alpha) * vec[0] - sin(alpha) * vec[1];
  vec[1]        = sin(alpha) * vx_tmp + cos(alpha) * vec[1];
}

/*! \brief Rotate a vector around the y axis.
 *
 *  \param[in, out] vec Array to 3 dimensional vector to be rotated.
 *  \param[in] alpha Rotation angle.
 *
 *  \return void
 */
static void rotate_y(double *vec, const double alpha)
{
  double vx_tmp = vec[0];

  vec[0] = cos(alpha) * vec[0] + sin(alpha) * vec[2];
  vec[2] = -sin(alpha) * vx_tmp + cos(alpha) * vec[2];
}

#endif /* #if !defined(TWODIMS) && !defined(ONEDIMS) */
