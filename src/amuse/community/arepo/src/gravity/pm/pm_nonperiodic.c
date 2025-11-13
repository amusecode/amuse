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
 * \file        src/gravity/pm/pm_non_periodic.c
 * \date        05/2018
 * \brief       Code for non-periodic FFT to compute long-range PM force.
 * \details     contains functions:
 *                void pm_init_regionsize(void)
 *                void pm_init_nonperiodic(void)
 *                int pmforce_is_particle_high_res(int type, MyDouble * Pos)
 *                void pmforce_nonperiodic_zoom_optimized_prepare_density(int
 *                  grnr)
 *                void pmforce_nonperiodic_zoom_optimized_readout_forces_or_
 *                  potential(int grnr, int dim)
 *                void pmforce_nonperiodic_uniform_optimized_prepare_density(
 *                  int grnr)
 *                void pmforce_nonperiodic_uniform_optimized_readout_forces_or_
 *                  potential(int grnr, int dim)
 *                int pmforce_nonperiodic(int grnr)
 *                void pm_setup_nonperiodic_kernel(void)
 *                static int pm_periodic_compare_sortindex(const void *a,
 *                  const void *b)
 *                static void msort_pmperiodic_with_tmp(large_numpart_type * b,
 *                  size_t n, large_numpart_type * t)
 *                static void mysort_pmperiodic(void *b, size_t n, size_t s,
 *                  int (*cmp) (const void *, const void *))
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 15.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../main/allvars.h"
#include "../../main/proto.h"

#if defined(PMGRID) && (defined(PLACEHIGHRESREGION) || defined(GRAVITY_NOT_PERIODIC))

#if defined(LONG_X) || defined(LONG_Y) || defined(LONG_Z)
#error "LONG_X/Y/Z not supported for the non-periodic FFT gravity code"
#endif /* #if defined(LONG_X) || defined(LONG_Y) || defined (LONG_Z) */

#ifndef GRIDBOOST
#define GRIDBOOST 2
#endif /* #ifndef GRIDBOOST */

#define GRID (GRIDBOOST * PMGRID)
#define GRIDz (GRID / 2 + 1)
#define GRID2 (2 * GRIDz)

#if(GRID > 1024)
typedef long long large_array_offset; /* use a larger data type in this case so that we can always address all cells of the 3D grid
                                         with a single index */
#else                                 /* #if (GRID > 1024) */
typedef unsigned int large_array_offset;
#endif                                /* #if (GRID > 1024) #else */

#ifdef NUMPART_PER_TASK_LARGE
typedef long long large_numpart_type; /* if there is a risk that the local particle number times 8 overflows a 32-bit integer, this
                                         data type should be used */
#else                                 /* #ifdef NUMPART_PER_TASK_LARGE */
typedef int large_numpart_type;
#endif                                /* #ifdef NUMPART_PER_TASK_LARGE */

/* short-cut macros for accessing different 3D arrays */
#define FI(x, y, z) (((large_array_offset)GRID2) * (GRID * (x) + (y)) + (z))
#define FC(c, z) (((large_array_offset)GRID2) * ((c)-myplan.base_firstcol) + (z))
#define TI(x, y, z) (((large_array_offset)GRID) * ((x) + (y)*myplan.nslab_x) + (z))

static fft_plan myplan; /*!< In this structure, various bookkeeping variables for the distributed FFTs are stored */

/*! \var maxfftsize
 *  \brief maximum size of the local fft grid among all tasks
 */
static size_t maxfftsize;

/*! \var rhogrid
 *  \brief This array hold the local part of the density field and
 *  after the FFTs the local part of the potential
 *
 *  \var forcegrid
 *  \brief This array will contain the force field
 *
 *  \var workspace
 *  \brief Workspace array used during the FFTs
 */
static fft_real *rhogrid, *forcegrid, *workspace;

/*! \brief Array containing the FFT of 'rhogrid'
 *
 *  This pointer points to the same array as 'rhogrid',
 *  because in-place FFTs are used.
 */
static fft_complex *fft_of_rhogrid;

static fft_real *kernel[2];
static fft_complex *fft_of_kernel[2];

/*! \param Determine particle extent.
 *
 *  This function determines the particle extension of all particles, and for
 *  those types selected with PLACEHIGHRESREGION if this is used, and then
 *  determines the boundaries of the non-periodic FFT-mesh that can be placed
 *  on this region. Note that a sufficient buffer region at the rim of the
 *  occupied part of the mesh needs to be reserved in order to allow a correct
 *  finite differencing using a 4-point formula. In addition, to allow
 *  non-periodic boundaries, the actual FFT mesh used is twice as large in
 *  each dimension compared with GRID.
 *
 *  \return void
 */
void pm_init_regionsize(void)
{
  double meshinner[2], xmin[2][3], xmax[2][3];
  int i, j;

  /* find enclosing rectangle */

  for(j = 0; j < 3; j++)
    {
      xmin[0][j] = xmin[1][j] = 1.0e36;
      xmax[0][j] = xmax[1][j] = -1.0e36;
    }

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      {
        if(P[i].Pos[j] > xmax[0][j])
          xmax[0][j] = P[i].Pos[j];
        if(P[i].Pos[j] < xmin[0][j])
          xmin[0][j] = P[i].Pos[j];

#ifdef PLACEHIGHRESREGION
        if(((1 << P[i].Type) & (PLACEHIGHRESREGION)))
          {
            if(P[i].Pos[j] > xmax[1][j])
              xmax[1][j] = P[i].Pos[j];
            if(P[i].Pos[j] < xmin[1][j])
              xmin[1][j] = P[i].Pos[j];
          }
#endif /* #ifdef PLACEHIGHRESREGION */
      }

  MPI_Allreduce(xmin, All.Xmintot, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(xmax, All.Xmaxtot, 6, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  for(j = 0; j < 2; j++)
    {
      All.TotalMeshSize[j] = All.Xmaxtot[j][0] - All.Xmintot[j][0];
      All.TotalMeshSize[j] = dmax(All.TotalMeshSize[j], All.Xmaxtot[j][1] - All.Xmintot[j][1]);
      All.TotalMeshSize[j] = dmax(All.TotalMeshSize[j], All.Xmaxtot[j][2] - All.Xmintot[j][2]);
#ifdef ENLARGEREGION
      All.TotalMeshSize[j] *= ENLARGEREGION;
#endif /* #ifdef ENLARGEREGION */

      /* symmetrize the box onto the center */
      for(i = 0; i < 3; i++)
        {
          All.Xmintot[j][i] = (All.Xmintot[j][i] + All.Xmaxtot[j][i]) / 2 - All.TotalMeshSize[j] / 2;
          All.Xmaxtot[j][i] = All.Xmintot[j][i] + All.TotalMeshSize[j];
        }
    }

  /* this will produce enough room for zero-padding and buffer region to
     allow finite differencing of the potential  */

  for(j = 0; j < 2; j++)
    {
      meshinner[j] = All.TotalMeshSize[j];
      All.TotalMeshSize[j] *= 2.001 * (GRID) / ((double)(GRID - 2 - 8));
    }

  /* move lower left corner by two cells to allow finite differencing of the potential by a 4-point function */

  for(j = 0; j < 2; j++)
    for(i = 0; i < 3; i++)
      {
        All.Corner[j][i]      = All.Xmintot[j][i] - 2.0005 * All.TotalMeshSize[j] / GRID;
        All.UpperCorner[j][i] = All.Corner[j][i] + (GRID / 2 - 1) * (All.TotalMeshSize[j] / GRID);
      }

#ifdef PLACEHIGHRESREGION
  All.Asmth[1] = ASMTH * All.TotalMeshSize[1] / GRID;
  All.Rcut[1]  = RCUT * All.Asmth[1];
#endif /* #ifdef PLACEHIGHRESREGION */

#ifdef PLACEHIGHRESREGION
  if(2 * All.TotalMeshSize[1] / GRID < All.Rcut[0])
    {
      All.TotalMeshSize[1] = 2 * (meshinner[1] + 2 * All.Rcut[0]) * (GRID) / ((double)(GRID - 2));

      for(i = 0; i < 3; i++)
        {
          All.Corner[1][i]      = All.Xmintot[1][i] - 1.0001 * All.Rcut[0];
          All.UpperCorner[1][i] = All.Corner[1][i] + (GRID / 2 - 1) * (All.TotalMeshSize[1] / GRID);
        }

      if(2 * All.TotalMeshSize[1] / GRID > All.Rcut[0])
        {
          All.TotalMeshSize[1] = 2 * (meshinner[1] + 2 * All.Rcut[0]) * (GRID) / ((double)(GRID - 10));

          for(i = 0; i < 3; i++)
            {
              All.Corner[1][i]      = All.Xmintot[1][i] - 1.0001 * (All.Rcut[0] + 2 * All.TotalMeshSize[1] / GRID);
              All.UpperCorner[1][i] = All.Corner[1][i] + (GRID / 2 - 1) * (All.TotalMeshSize[1] / GRID);
            }
        }

      All.Asmth[1] = ASMTH * All.TotalMeshSize[1] / GRID;
      All.Rcut[1]  = RCUT * All.Asmth[1];

      mpi_printf("PM-NONPERIODIC: All.Asmth[0]=%g All.Asmth[1]=%g\n", All.Asmth[0], All.Asmth[1]);
    }
#endif /* #ifdef PLACEHIGHRESREGION */

#ifdef PLACEHIGHRESREGION
  mpi_printf(
      "PM-NONPERIODIC: Allowed region for isolated PM mesh (high-res): (%g|%g|%g)  -> (%g|%g|%g)   ext=%g  totmeshsize=%g  "
      "meshsize=%g\n\n",
      All.Xmintot[1][0], All.Xmintot[1][1], All.Xmintot[1][2], All.Xmaxtot[1][0], All.Xmaxtot[1][1], All.Xmaxtot[1][2], meshinner[1],
      All.TotalMeshSize[1], All.TotalMeshSize[1] / GRID);
#endif /* #ifdef PLACEHIGHRESREGION */
}

/*! \brief Initialization of the non-periodic PM routines.
 *
 *  The plan-files for FFTW are created. Finally, the routine to set-up the
 *  non-periodic Greens function is called.
 *
 *  \return void
 */
void pm_init_nonperiodic(void)
{
  /* Set up the FFTW-3 plan files. */
  int ndim[1] = {GRID}; /* dimension of the 1D transforms */

  /* temporarily allocate some arrays to make sure that out-of-place plans are created */
  rhogrid   = (fft_real *)mymalloc("rhogrid", GRID2 * sizeof(fft_real));
  forcegrid = (fft_real *)mymalloc("forcegrid", GRID2 * sizeof(fft_real));

#ifdef DOUBLEPRECISION_FFTW
  int alignflag = 0;
#else  /* #ifdef DOUBLEPRECISION_FFTW */
  /* for single precision, the start of our FFT columns is presently only guaranteed to be 8-byte aligned */
  int alignflag = FFTW_UNALIGNED;
#endif /* #ifdef DOUBLEPRECISION_FFTW #else */
#ifndef FFT_COLUMN_BASED
  int stride = GRIDz;
#else  /* #ifndef FFT_COLUMN_BASED */
  int stride    = 1;
#endif /* #ifndef FFT_COLUMN_BASED #else */

  myplan.forward_plan_zdir = FFTW(plan_many_dft_r2c)(1, ndim, 1, rhogrid, 0, 1, GRID2, (fft_complex *)forcegrid, 0, 1, GRIDz,
                                                     FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.forward_plan_xdir =
      FFTW(plan_many_dft)(1, ndim, 1, (fft_complex *)rhogrid, 0, stride, GRIDz * GRID, (fft_complex *)forcegrid, 0, stride,
                          GRIDz * GRID, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.forward_plan_ydir =
      FFTW(plan_many_dft)(1, ndim, 1, (fft_complex *)rhogrid, 0, stride, GRIDz * GRID, (fft_complex *)forcegrid, 0, stride,
                          GRIDz * GRID, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_zdir = FFTW(plan_many_dft_c2r)(1, ndim, 1, (fft_complex *)rhogrid, 0, 1, GRIDz, forcegrid, 0, 1, GRID2,
                                                      FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_xdir =
      FFTW(plan_many_dft)(1, ndim, 1, (fft_complex *)rhogrid, 0, stride, GRIDz * GRID, (fft_complex *)forcegrid, 0, stride,
                          GRIDz * GRID, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_ydir =
      FFTW(plan_many_dft)(1, ndim, 1, (fft_complex *)rhogrid, 0, stride, GRIDz * GRID, (fft_complex *)forcegrid, 0, stride,
                          GRIDz * GRID, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myfree(forcegrid);
  myfree(rhogrid);

#ifndef FFT_COLUMN_BASED

  my_slab_based_fft_init(&myplan, GRID, GRID, GRID);

  maxfftsize = myplan.largest_x_slab * GRID * ((size_t)GRID2);

#else /* #ifndef FFT_COLUMN_BASED */

  my_column_based_fft_init(&myplan, GRID, GRID, GRID);

  maxfftsize = myplan.max_datasize;

#endif /* #ifndef FFT_COLUMN_BASED #else */

  /* now allocate memory to hold the FFT fields */

  size_t bytes, bytes_tot = 0;

#if defined(GRAVITY_NOT_PERIODIC)
  kernel[0] = (fft_real *)mymalloc("kernel[0]", bytes = maxfftsize * sizeof(fft_real));
  bytes_tot += bytes;
  fft_of_kernel[0] = (fft_complex *)kernel[0];
#endif /* #if defined(GRAVITY_NOT_PERIODIC) */

#if defined(PLACEHIGHRESREGION)
  kernel[1] = (fft_real *)mymalloc("kernel[1]", bytes = maxfftsize * sizeof(fft_real));
  bytes_tot += bytes;
  fft_of_kernel[1] = (fft_complex *)kernel[1];
#endif /* #if defined(PLACEHIGHRESREGION) */

  mpi_printf("\nPM-NONPERIODIC: Allocated %g MByte for FFT kernel(s).\n\n", bytes_tot / (1024.0 * 1024.0));
}

#ifdef PLACEHIGHRESREGION
/*! \brief Is this a high res particle in high resolution region?
 *
 *  For cosmological zoom simulations.
 *
 *  \param[in] type Parcile type.
 *  \param[in] Pos Position of particle.
 *
 *  \return 0: not high res; 1: high res.
 */
int pmforce_is_particle_high_res(int type, MyDouble *Pos)
{
  int flag = 1;

  if((1 << type) & (PLACEHIGHRESREGION))
    return 1;

#if defined(PLACEHIGHRESREGION) && (FORCETEST_TESTFORCELAW == 1)
  double r2 = 0;
  for(int j = 0; j < 3; j++)
    r2 += pow(Pos[j] - 0.5 * (All.Xmintot[1][j] + All.Xmaxtot[1][j]), 2);

  if(sqrt(r2) > 0.5 * (All.Xmaxtot[1][0] - All.Xmintot[1][0]))
    return 0;
#else /* #if defined(PLACEHIGHRESREGION) && (FORCETEST_TESTFORCELAW == 1) */

  for(int j = 0; j < 3; j++)
    if(Pos[j] < All.Xmintot[1][j] || Pos[j] > All.Xmaxtot[1][j])
      {
        flag = 0; /* we are outside */
        break;
      }

#endif /* #if defined(PLACEHIGHRESREGION) && (FORCETEST_TESTFORCELAW == 1) #else */

  return flag;
}
#endif /* #ifdef PLACEHIGHRESREGION */

#ifdef PM_ZOOM_OPTIMIZED

static void mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *));
static int pm_periodic_compare_sortindex(const void *a, const void *b);

/*! \brief This structure links the particles to the mesh cells, to which they
 *         contribute their mass.
 *
 *  Each particle will have eight items of this structure in the 'part' array.
 *  For each of the eight mesh cells the CIC assignment will contribute,
 *  one item of this struct exists.
 */
static struct part_slab_data
{
  large_array_offset globalindex; /*!< index in the global density mesh */
  large_numpart_type partindex; /*!< contains the local particle index shifted by 2^3, the first three bits encode to which part of the
                                   CIC assignment this item belongs to */
  large_array_offset localindex; /*!< index to a local copy of the corresponding mesh cell of the global density array (used during
                                    local mass and force assignment) */
} * part;                        /*!< array of part_slab_data linking the local particles to their mesh cells */

static size_t *localfield_sendcount, *localfield_first, *localfield_offset, *localfield_recvcount;
static large_array_offset *localfield_globalindex, *import_globalindex;
static fft_real *localfield_data, *import_data;
static large_numpart_type num_on_grid;

/*! \brief Prepares density field for nonperiodic FFTs.
 *
 *  \param[in] grnr (0, 1) 0 if full mesh, 1 if highres grid.
 *
 *  \return void
 */
void pmforce_nonperiodic_zoom_optimized_prepare_density(int grnr)
{
  large_numpart_type i;
  int level, recvTask;
  MPI_Status status;

  double to_slab_fac = GRID / All.TotalMeshSize[grnr];

  part                               = (struct part_slab_data *)mymalloc("part", 8 * (NumPart * sizeof(struct part_slab_data)));
  large_numpart_type *part_sortindex = (large_numpart_type *)mymalloc("part_sortindex", 8 * (NumPart * sizeof(large_numpart_type)));

  int ngrid = 0;

  /* determine the cells each particle accesses */
  for(i = 0; i < NumPart; i++)
    {
      MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
      if(P[i].Type == 0)
        pos = SphP[i].Center;
      else
#endif /* #ifdef CELL_CENTER_GRAVITY */
        pos = P[i].Pos;

      if(pos[0] < All.Corner[grnr][0] || pos[0] >= All.UpperCorner[grnr][0])
        continue;
      if(pos[1] < All.Corner[grnr][1] || pos[1] >= All.UpperCorner[grnr][1])
        continue;
      if(pos[2] < All.Corner[grnr][2] || pos[2] >= All.UpperCorner[grnr][2])
        continue;

      int slab_x = (int)(to_slab_fac * (pos[0] - All.Corner[grnr][0]));
      int slab_y = (int)(to_slab_fac * (pos[1] - All.Corner[grnr][1]));
      int slab_z = (int)(to_slab_fac * (pos[2] - All.Corner[grnr][2]));
      int myngrid;

      {
        myngrid = ngrid;
        ngrid += 1;
      }

      large_numpart_type index_on_grid = ((large_numpart_type)myngrid) * 8;

      int xx, yy, zz;

      for(xx = 0; xx < 2; xx++)
        for(yy = 0; yy < 2; yy++)
          for(zz = 0; zz < 2; zz++)
            {
              int slab_xx = slab_x + xx;
              int slab_yy = slab_y + yy;
              int slab_zz = slab_z + zz;

              if(slab_xx >= GRID)
                slab_xx -= GRID;
              if(slab_yy >= GRID)
                slab_yy -= GRID;
              if(slab_zz >= GRID)
                slab_zz -= GRID;

              large_array_offset offset = FI(slab_xx, slab_yy, slab_zz);

              part[index_on_grid].partindex   = (i << 3) + (xx << 2) + (yy << 1) + zz;
              part[index_on_grid].globalindex = offset;
              part_sortindex[index_on_grid]   = index_on_grid;
              index_on_grid++;
            }
    }

  /* note: num_on_grid will be  8 times larger than the particle number, but num_field_points will generally be much smaller */
  num_on_grid = ((large_numpart_type)ngrid) * 8;

  /* bring the part-field into the order of the accessed cells. This allows the removal of duplicates */
  mysort_pmperiodic(part_sortindex, num_on_grid, sizeof(large_numpart_type), pm_periodic_compare_sortindex);

  large_array_offset num_field_points;

  if(num_on_grid > 0)
    num_field_points = 1;
  else
    num_field_points = 0;

  /* determine the number of unique field points */
  for(i = 1; i < num_on_grid; i++)
    {
      if(part[part_sortindex[i]].globalindex != part[part_sortindex[i - 1]].globalindex)
        num_field_points++;
    }

  /* allocate the local field */
  localfield_globalindex = (large_array_offset *)mymalloc_movable(&localfield_globalindex, "localfield_globalindex",
                                                                  num_field_points * sizeof(large_array_offset));
  localfield_data        = (fft_real *)mymalloc_movable(&localfield_data, "localfield_data", num_field_points * sizeof(fft_real));
  localfield_first       = (size_t *)mymalloc_movable(&localfield_first, "localfield_first", NTask * sizeof(size_t));
  localfield_sendcount   = (size_t *)mymalloc_movable(&localfield_sendcount, "localfield_sendcount", NTask * sizeof(size_t));
  localfield_offset      = (size_t *)mymalloc_movable(&localfield_offset, "localfield_offset", NTask * sizeof(size_t));
  localfield_recvcount   = (size_t *)mymalloc_movable(&localfield_recvcount, "localfield_recvcount", NTask * sizeof(size_t));

  for(i = 0; i < NTask; i++)
    {
      localfield_first[i]     = 0;
      localfield_sendcount[i] = 0;
    }

  /* establish the cross link between the part[ ]-array and the local list of
   * mesh points. Also, count on which CPU the needed field points are stored.
   */
  for(i = 0, num_field_points = 0; i < num_on_grid; i++)
    {
      if(i > 0)
        if(part[part_sortindex[i]].globalindex != part[part_sortindex[i - 1]].globalindex)
          num_field_points++;

      part[part_sortindex[i]].localindex = num_field_points;

      if(i > 0)
        if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
          continue;

      localfield_globalindex[num_field_points] = part[part_sortindex[i]].globalindex;

#ifndef FFT_COLUMN_BASED
      int slab = part[part_sortindex[i]].globalindex / (GRID * GRID2);
      int task = myplan.slab_to_task[slab];
#else  /* #ifndef FFT_COLUMN_BASED */
      int task, column = part[part_sortindex[i]].globalindex / (GRID2);

      if(column < myplan.pivotcol)
        task = column / myplan.avg;
      else
        task = (column - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;
#endif /* #ifndef FFT_COLUMN_BASED #else */

      if(localfield_sendcount[task] == 0)
        localfield_first[task] = num_field_points;

      localfield_sendcount[task]++;
    }
  num_field_points++;

  for(i = 1, localfield_offset[0] = 0; i < NTask; i++)
    localfield_offset[i] = localfield_offset[i - 1] + localfield_sendcount[i - 1];

  myfree_movable(part_sortindex);
  part_sortindex = NULL;

  /* now bin the local particle data onto the mesh list */
  for(i = 0; i < num_field_points; i++)
    localfield_data[i] = 0;

  for(i = 0; i < num_on_grid; i += 8)
    {
      int pindex = (part[i].partindex >> 3);

      MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
      if(P[pindex].Type == 0)
        pos = SphP[pindex].Center;
      else
#endif /* #ifdef CELL_CENTER_GRAVITY */
        pos = P[pindex].Pos;

      int slab_x = (int)(to_slab_fac * (pos[0] - All.Corner[grnr][0]));
      int slab_y = (int)(to_slab_fac * (pos[1] - All.Corner[grnr][1]));
      int slab_z = (int)(to_slab_fac * (pos[2] - All.Corner[grnr][2]));

      double dx = to_slab_fac * (pos[0] - All.Corner[grnr][0]) - slab_x;
      double dy = to_slab_fac * (pos[1] - All.Corner[grnr][1]) - slab_y;
      double dz = to_slab_fac * (pos[2] - All.Corner[grnr][2]) - slab_z;

      double weight = P[pindex].Mass;

      localfield_data[part[i + 0].localindex] += weight * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
      localfield_data[part[i + 1].localindex] += weight * (1.0 - dx) * (1.0 - dy) * dz;
      localfield_data[part[i + 2].localindex] += weight * (1.0 - dx) * dy * (1.0 - dz);
      localfield_data[part[i + 3].localindex] += weight * (1.0 - dx) * dy * dz;
      localfield_data[part[i + 4].localindex] += weight * (dx) * (1.0 - dy) * (1.0 - dz);
      localfield_data[part[i + 5].localindex] += weight * (dx) * (1.0 - dy) * dz;
      localfield_data[part[i + 6].localindex] += weight * (dx)*dy * (1.0 - dz);
      localfield_data[part[i + 7].localindex] += weight * (dx)*dy * dz;
    }

  rhogrid = (fft_real *)mymalloc("rhogrid", maxfftsize * sizeof(fft_real));

  /* clear local FFT-mesh density field */
  large_array_offset ii;
  for(ii = 0; ii < maxfftsize; ii++)
    rhogrid[ii] = 0;

  /* exchange data and add contributions to the local mesh-path */
  MPI_Alltoall(localfield_sendcount, sizeof(size_t), MPI_BYTE, localfield_recvcount, sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);

  for(level = 0; level < (1 << PTask); level++) /* note: for level=0, target is the same task */
    {
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
        {
          if(level > 0)
            {
              import_data = (fft_real *)mymalloc("import_data", localfield_recvcount[recvTask] * sizeof(fft_real));
              import_globalindex =
                  (large_array_offset *)mymalloc("import_globalindex", localfield_recvcount[recvTask] * sizeof(large_array_offset));

              if(localfield_sendcount[recvTask] > 0 || localfield_recvcount[recvTask] > 0)
                {
                  myMPI_Sendrecv(localfield_data + localfield_offset[recvTask], localfield_sendcount[recvTask] * sizeof(fft_real),
                                 MPI_BYTE, recvTask, TAG_NONPERIOD_A, import_data, localfield_recvcount[recvTask] * sizeof(fft_real),
                                 MPI_BYTE, recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

                  myMPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
                                 localfield_sendcount[recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask, TAG_NONPERIOD_B,
                                 import_globalindex, localfield_recvcount[recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask,
                                 TAG_NONPERIOD_B, MPI_COMM_WORLD, &status);
                }
            }
          else
            {
              import_data        = localfield_data + localfield_offset[ThisTask];
              import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
            }

          /* note: here every element in rhogrid is only accessed once, so there should be no race condition */
          for(i = 0; i < localfield_recvcount[recvTask]; i++)
            {
              /* determine offset in local FFT slab */
#ifndef FFT_COLUMN_BASED
              large_array_offset offset =
                  import_globalindex[i] - myplan.first_slab_x_of_task[ThisTask] * GRID * ((large_array_offset)GRID2);
#else  /* #ifndef FFT_COLUMN_BASED */
              large_array_offset offset = import_globalindex[i] - myplan.base_firstcol * ((large_array_offset)GRID2);
#endif /* #ifndef FFT_COLUMN_BASED #else */
              rhogrid[offset] += import_data[i];
            }

          if(level > 0)
            {
              myfree(import_globalindex);
              myfree(import_data);
            }
        }
    }
}

/*! \brief Reads out the force component corresponding to spatial dimension
 *         'dim'.
 *
 *  If dim is negative, potential values are read out and assigned to
 *  particles.
 *
 *  \param[in] grnr Number of grid (0: base, 1 high-res)
 *  \param[in] dim Dimension to be read out
 *             (<0: potential,>=0 force component).
 *
 *  \return void
 */
void pmforce_nonperiodic_zoom_optimized_readout_forces_or_potential(int grnr, int dim)
{
#ifdef EVALPOTENTIAL
  /* factor to get potential */
  double fac = All.G / pow(All.TotalMeshSize[grnr], 4) * pow(All.TotalMeshSize[grnr] / GRID, 3);
#endif /* #ifdef EVALPOTENTIAL */

  large_numpart_type i;
  int level, recvTask;
  MPI_Status status;

  fft_real *grid;

  if(dim < 0)
    grid = rhogrid;
  else
    grid = forcegrid;

  double to_slab_fac = GRID / All.TotalMeshSize[grnr];

  for(level = 0; level < (1 << PTask); level++) /* note: for level=0, target is the same task */
    {
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
        {
          if(level > 0)
            {
              import_data = (fft_real *)mymalloc("import_data", localfield_recvcount[recvTask] * sizeof(fft_real));
              import_globalindex =
                  (large_array_offset *)mymalloc("import_globalindex", localfield_recvcount[recvTask] * sizeof(large_array_offset));

              if(localfield_sendcount[recvTask] > 0 || localfield_recvcount[recvTask] > 0)
                {
                  myMPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
                                 localfield_sendcount[recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask, TAG_NONPERIOD_C,
                                 import_globalindex, localfield_recvcount[recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask,
                                 TAG_NONPERIOD_C, MPI_COMM_WORLD, &status);
                }
            }
          else
            {
              import_data        = localfield_data + localfield_offset[ThisTask];
              import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
            }

          for(i = 0; i < localfield_recvcount[recvTask]; i++)
            {
#ifndef FFT_COLUMN_BASED
              large_array_offset offset =
                  import_globalindex[i] - myplan.first_slab_x_of_task[ThisTask] * GRID * ((large_array_offset)GRID2);
#else  /* #ifndef FFT_COLUMN_BASED */
              large_array_offset offset = import_globalindex[i] - myplan.base_firstcol * ((large_array_offset)GRID2);
#endif /* #ifndef FFT_COLUMN_BASED #else */
              import_data[i] = grid[offset];
            }

          if(level > 0)
            {
              myMPI_Sendrecv(import_data, localfield_recvcount[recvTask] * sizeof(fft_real), MPI_BYTE, recvTask, TAG_NONPERIOD_A,
                             localfield_data + localfield_offset[recvTask], localfield_sendcount[recvTask] * sizeof(fft_real),
                             MPI_BYTE, recvTask, TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

              myfree(import_globalindex);
              myfree(import_data);
            }
        }
    }

  /* read out the force/potential values, which all have been assembled in localfield_data */

  int k, ngrid = (num_on_grid >> 3);

  for(k = 0; k < ngrid; k++)
    {
      large_numpart_type j = (((large_numpart_type)k) << 3);

      int i = (part[j].partindex >> 3);

      MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
      if(P[i].Type == 0)
        pos = SphP[i].Center;
      else
#endif /* #ifdef CELL_CENTER_GRAVITY */
        pos = P[i].Pos;

#ifdef PLACEHIGHRESREGION
      if(grnr == 1)
        if(!(pmforce_is_particle_high_res(P[i].Type, pos)))
          continue;
#endif /* #ifdef PLACEHIGHRESREGION */

      int slab_x = (int)(to_slab_fac * (pos[0] - All.Corner[grnr][0]));
      double dx  = to_slab_fac * (pos[0] - All.Corner[grnr][0]) - slab_x;

      int slab_y = (int)(to_slab_fac * (pos[1] - All.Corner[grnr][1]));
      double dy  = to_slab_fac * (pos[1] - All.Corner[grnr][1]) - slab_y;

      int slab_z = (int)(to_slab_fac * (pos[2] - All.Corner[grnr][2]));
      double dz  = to_slab_fac * (pos[2] - All.Corner[grnr][2]) - slab_z;

      double value = +localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                     localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz +
                     localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz) +
                     localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz +
                     localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz) +
                     localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz +
                     localfield_data[part[j + 6].localindex] * (dx)*dy * (1.0 - dz) +
                     localfield_data[part[j + 7].localindex] * (dx)*dy * dz;

      if(dim < 0)
        {
#ifdef EVALPOTENTIAL
          P[i].PM_Potential += value * fac;
#endif /* #ifdef EVALPOTENTIAL */
        }
      else
        P[i].GravPM[dim] += value;
    }
}

#else /* #ifdef PM_ZOOM_OPTIMIZED */
/* Here come the routines for a different communication algorithm that is better suited for a homogenuously loaded boxes.
 */

/*! \brief Particle buffer structure
 */
static struct partbuf
{
  MyFloat Mass;
  MyFloat Pos[3];
} * partin, *partout;

static size_t nimport, nexport;

static size_t *Sndpm_count, *Sndpm_offset;
static size_t *Rcvpm_count, *Rcvpm_offset;

/*! \brief Prepares density for pm calculation in algorithm optimized for
 *         uniform densities.
 *
 *  \param[in] grnr Number of grid (0: base grid, 1: high res grid).
 *
 *  \return void
 */
void pmforce_nonperiodic_uniform_optimized_prepare_density(int grnr)
{
  int i, j;

  double to_slab_fac = GRID / All.TotalMeshSize[grnr];

  /* We here enlarge NTask such that each thread gets his own cache line for send_count/send_offset.
   * This should hopefully prevent a performance penalty from 'false sharing' for these variables
   */
  int multiNtask = roundup_to_multiple_of_cacheline_size(NTask * sizeof(size_t)) / sizeof(size_t);

  Sndpm_count  = mymalloc("Sndpm_count", MaxThreads * multiNtask * sizeof(size_t));
  Sndpm_offset = mymalloc("Sndpm_offset", MaxThreads * multiNtask * sizeof(size_t));
  Rcvpm_count  = mymalloc("Rcvpm_count", NTask * sizeof(size_t));
  Rcvpm_offset = mymalloc("Rcvpm_offset", NTask * sizeof(size_t));

  /* determine the slabs/columns each particles accesses */
  {
    size_t *send_count = Sndpm_count + get_thread_num() * multiNtask;

    /* each threads needs to do theloop to clear its send_count[] array */
    for(j = 0; j < NTask; j++)
      send_count[j] = 0;

    for(i = 0; i < NumPart; i++)
      {
        MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
        if(P[i].Type == 0)
          pos = SphP[i].Center;
        else
#endif /* #ifdef CELL_CENTER_GRAVITY */
          pos = P[i].Pos;

        if(pos[0] < All.Corner[grnr][0] || pos[0] >= All.UpperCorner[grnr][0])
          continue;
        if(pos[1] < All.Corner[grnr][1] || pos[1] >= All.UpperCorner[grnr][1])
          continue;
        if(pos[2] < All.Corner[grnr][2] || pos[2] >= All.UpperCorner[grnr][2])
          continue;

        int slab_x  = (int)(to_slab_fac * (pos[0] - All.Corner[grnr][0]));
        int slab_xx = slab_x + 1;

#ifndef FFT_COLUMN_BASED
        int task0   = myplan.slab_to_task[slab_x];
        int task1   = myplan.slab_to_task[slab_xx];

        send_count[task0]++;
        if(task0 != task1)
          send_count[task1]++;
#else  /* #ifndef FFT_COLUMN_BASED */
        int slab_y  = (int)(to_slab_fac * (pos[1] - All.Corner[grnr][1]));
        int slab_yy = slab_y + 1;

        int column0 = slab_x * GRID + slab_y;
        int column1 = slab_x * GRID + slab_yy;
        int column2 = slab_xx * GRID + slab_y;
        int column3 = slab_xx * GRID + slab_yy;

        int task0, task1, task2, task3;

        if(column0 < myplan.pivotcol)
          task0 = column0 / myplan.avg;
        else
          task0 = (column0 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column1 < myplan.pivotcol)
          task1 = column1 / myplan.avg;
        else
          task1 = (column1 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column2 < myplan.pivotcol)
          task2 = column2 / myplan.avg;
        else
          task2 = (column2 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column3 < myplan.pivotcol)
          task3 = column3 / myplan.avg;
        else
          task3 = (column3 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        send_count[task0]++;
        if(task1 != task0)
          send_count[task1]++;
        if(task2 != task1 && task2 != task0)
          send_count[task2]++;
        if(task3 != task0 && task3 != task1 && task3 != task2)
          send_count[task3]++;
#endif /* #ifndef FFT_COLUMN_BASED #else */
      }
  }

  /* collect thread-specific offset table and collect the results from the other threads */
  for(i = 0, Sndpm_offset[0] = 0; i < NTask; i++)
    for(j = 0; j < MaxThreads; j++)
      {
        int ind_prev, ind = j * multiNtask + i;
        if(ind > 0)
          {
            if(j == 0)
              ind_prev = (MaxThreads - 1) * multiNtask + i - 1;
            else
              ind_prev = ind - multiNtask;

            Sndpm_offset[ind] = Sndpm_offset[ind_prev] + Sndpm_count[ind_prev];
          }
      }

  for(j = 1; j < MaxThreads; j++)
    for(i = 0; i < NTask; i++)
      Sndpm_count[i] += Sndpm_count[i + j * multiNtask];

  MPI_Alltoall(Sndpm_count, sizeof(size_t), MPI_BYTE, Rcvpm_count, sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, nexport = 0, Rcvpm_offset[0] = 0, Sndpm_offset[0] = 0; j < NTask; j++)
    {
      nexport += Sndpm_count[j];
      nimport += Rcvpm_count[j];

      if(j > 0)
        {
          Sndpm_offset[j] = Sndpm_offset[j - 1] + Sndpm_count[j - 1];
          Rcvpm_offset[j] = Rcvpm_offset[j - 1] + Rcvpm_count[j - 1];
        }
    }

  /* allocate import and export buffer */
  partin  = (struct partbuf *)mymalloc("partin", nimport * sizeof(struct partbuf));
  partout = (struct partbuf *)mymalloc("partout", nexport * sizeof(struct partbuf));

  {
    size_t *send_count  = Sndpm_count + get_thread_num() * multiNtask;
    size_t *send_offset = Sndpm_offset + get_thread_num() * multiNtask;

    for(j = 0; j < NTask; j++)
      send_count[j] = 0;

    /* fill export buffer */
    for(i = 0; i < NumPart; i++)
      {
        MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
        if(P[i].Type == 0)
          pos = SphP[i].Center;
        else
#endif /* #ifdef CELL_CENTER_GRAVITY */
          pos = P[i].Pos;

        if(pos[0] < All.Corner[grnr][0] || pos[0] >= All.UpperCorner[grnr][0])
          continue;
        if(pos[1] < All.Corner[grnr][1] || pos[1] >= All.UpperCorner[grnr][1])
          continue;
        if(pos[2] < All.Corner[grnr][2] || pos[2] >= All.UpperCorner[grnr][2])
          continue;

        int slab_x  = (int)(to_slab_fac * (pos[0] - All.Corner[grnr][0]));
        int slab_xx = slab_x + 1;

#ifndef FFT_COLUMN_BASED
        int task0   = myplan.slab_to_task[slab_x];
        int task1   = myplan.slab_to_task[slab_xx];

        size_t ind0        = send_offset[task0] + send_count[task0]++;
        partout[ind0].Mass = P[i].Mass;
        for(j = 0; j < 3; j++)
          partout[ind0].Pos[j] = pos[j];

        if(task0 != task1)
          {
            size_t ind1        = send_offset[task1] + send_count[task1]++;
            partout[ind1].Mass = P[i].Mass;
            for(j = 0; j < 3; j++)
              partout[ind1].Pos[j] = pos[j];
          }
#else  /* #ifndef FFT_COLUMN_BASED */
        int slab_y  = (int)(to_slab_fac * (pos[1] - All.Corner[grnr][1]));
        int slab_yy = slab_y + 1;

        int column0 = slab_x * GRID + slab_y;
        int column1 = slab_x * GRID + slab_yy;
        int column2 = slab_xx * GRID + slab_y;
        int column3 = slab_xx * GRID + slab_yy;

        int task0, task1, task2, task3;

        if(column0 < myplan.pivotcol)
          task0 = column0 / myplan.avg;
        else
          task0 = (column0 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column1 < myplan.pivotcol)
          task1 = column1 / myplan.avg;
        else
          task1 = (column1 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column2 < myplan.pivotcol)
          task2 = column2 / myplan.avg;
        else
          task2 = (column2 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column3 < myplan.pivotcol)
          task3 = column3 / myplan.avg;
        else
          task3 = (column3 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        size_t ind0        = send_offset[task0] + send_count[task0]++;
        partout[ind0].Mass = P[i].Mass;
        for(j = 0; j < 3; j++)
          partout[ind0].Pos[j] = pos[j];

        if(task1 != task0)
          {
            size_t ind1        = send_offset[task1] + send_count[task1]++;
            partout[ind1].Mass = P[i].Mass;
            for(j = 0; j < 3; j++)
              partout[ind1].Pos[j] = pos[j];
          }
        if(task2 != task1 && task2 != task0)
          {
            size_t ind2        = send_offset[task2] + send_count[task2]++;
            partout[ind2].Mass = P[i].Mass;
            for(j = 0; j < 3; j++)
              partout[ind2].Pos[j] = pos[j];
          }
        if(task3 != task0 && task3 != task1 && task3 != task2)
          {
            size_t ind3        = send_offset[task3] + send_count[task3]++;
            partout[ind3].Mass = P[i].Mass;
            for(j = 0; j < 3; j++)
              partout[ind3].Pos[j] = pos[j];
          }
#endif /* #ifndef FFT_COLUMN_BASED #else */
      }
  }

  /* collect the send_count[] results from the other threads */
  for(j = 1; j < MaxThreads; j++)
    for(i = 0; i < NTask; i++)
      Sndpm_count[i] += Sndpm_count[i + j * multiNtask];

  int flag_big = 0, flag_big_all;
  for(i = 0; i < NTask; i++)
    if(Sndpm_count[i] * sizeof(struct partbuf) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  /* exchange particle data */
  myMPI_Alltoallv(partout, Sndpm_count, Sndpm_offset, partin, Rcvpm_count, Rcvpm_offset, sizeof(struct partbuf), flag_big_all,
                  MPI_COMM_WORLD);

  myfree(partout);

  /* allocate density field */
  rhogrid = (fft_real *)mymalloc("rhogrid", maxfftsize * sizeof(fft_real));

  /* clear local FFT-mesh density field */
  large_array_offset ii;

  for(ii = 0; ii < maxfftsize; ii++)
    rhogrid[ii] = 0;

#ifndef FFT_COLUMN_BASED
  /* bin particle data onto mesh, in multi-threaded fashion */
  {
    int tid = get_thread_num();

    int first_y, count_y;
    subdivide_evenly(GRID, MaxThreads, tid, &first_y, &count_y);
    int last_y = first_y + count_y - 1;

    for(i = 0; i < nimport; i++)
      {
        int slab_y  = (int)(to_slab_fac * (partin[i].Pos[1] - All.Corner[grnr][1]));
        int slab_yy = slab_y + 1;
        double dy   = to_slab_fac * (partin[i].Pos[1] - All.Corner[grnr][1]) - slab_y;
        int flag_slab_y, flag_slab_yy;

        if(slab_y >= first_y && slab_y <= last_y)
          flag_slab_y = 1;
        else
          flag_slab_y = 0;

        if(slab_yy >= first_y && slab_yy <= last_y)
          flag_slab_yy = 1;
        else
          flag_slab_yy = 0;

        if(flag_slab_y || flag_slab_yy)
          {
            double mass = partin[i].Mass;

            int slab_x  = (int)(to_slab_fac * (partin[i].Pos[0] - All.Corner[grnr][0]));
            int slab_z  = (int)(to_slab_fac * (partin[i].Pos[2] - All.Corner[grnr][2]));
            int slab_xx = slab_x + 1;
            int slab_zz = slab_z + 1;

            double dx = to_slab_fac * (partin[i].Pos[0] - All.Corner[grnr][0]) - slab_x;
            double dz = to_slab_fac * (partin[i].Pos[2] - All.Corner[grnr][2]) - slab_z;

            int flag_slab_x, flag_slab_xx;

            if(myplan.slab_to_task[slab_x] == ThisTask)
              {
                slab_x -= myplan.first_slab_x_of_task[ThisTask];
                flag_slab_x = 1;
              }
            else
              flag_slab_x = 0;

            if(myplan.slab_to_task[slab_xx] == ThisTask)
              {
                slab_xx -= myplan.first_slab_x_of_task[ThisTask];
                flag_slab_xx = 1;
              }
            else
              flag_slab_xx = 0;

            if(flag_slab_x)
              {
                if(flag_slab_y)
                  {
                    rhogrid[FI(slab_x, slab_y, slab_z)] += (mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz));
                    rhogrid[FI(slab_x, slab_y, slab_zz)] += (mass * (1.0 - dx) * (1.0 - dy) * (dz));
                  }

                if(flag_slab_yy)
                  {
                    rhogrid[FI(slab_x, slab_yy, slab_z)] += (mass * (1.0 - dx) * (dy) * (1.0 - dz));
                    rhogrid[FI(slab_x, slab_yy, slab_zz)] += (mass * (1.0 - dx) * (dy) * (dz));
                  }
              }

            if(flag_slab_xx)
              {
                if(flag_slab_y)
                  {
                    rhogrid[FI(slab_xx, slab_y, slab_z)] += (mass * (dx) * (1.0 - dy) * (1.0 - dz));
                    rhogrid[FI(slab_xx, slab_y, slab_zz)] += (mass * (dx) * (1.0 - dy) * (dz));
                  }

                if(flag_slab_yy)
                  {
                    rhogrid[FI(slab_xx, slab_yy, slab_z)] += (mass * (dx) * (dy) * (1.0 - dz));
                    rhogrid[FI(slab_xx, slab_yy, slab_zz)] += (mass * (dx) * (dy) * (dz));
                  }
              }
          }
      }
  }

#else /* #ifndef FFT_COLUMN_BASED */

  struct data_cols
  {
    int col0, col1, col2, col3;
    double dx, dy;
  } * aux;

  aux = mymalloc("aux", nimport * sizeof(struct data_cols));

  for(i = 0; i < nimport; i++)
    {
      int slab_x = (int)(to_slab_fac * (partin[i].Pos[0] - All.Corner[grnr][0]));
      int slab_xx = slab_x + 1;

      int slab_y = (int)(to_slab_fac * (partin[i].Pos[1] - All.Corner[grnr][1]));
      int slab_yy = slab_y + 1;

      aux[i].dx = to_slab_fac * (partin[i].Pos[0] - All.Corner[grnr][0]) - slab_x;
      aux[i].dy = to_slab_fac * (partin[i].Pos[1] - All.Corner[grnr][1]) - slab_y;

      aux[i].col0 = slab_x * GRID + slab_y;
      aux[i].col1 = slab_x * GRID + slab_yy;
      aux[i].col2 = slab_xx * GRID + slab_y;
      aux[i].col3 = slab_xx * GRID + slab_yy;
    }

  {
    int tid = get_thread_num();

    int first_col, last_col, count_col;
    subdivide_evenly(myplan.base_ncol, MaxThreads, tid, &first_col, &count_col);
    last_col = first_col + count_col - 1;
    first_col += myplan.base_firstcol;
    last_col += myplan.base_firstcol;

    for(i = 0; i < nimport; i++)
      {
        int flag0, flag1, flag2, flag3;
        int col0 = aux[i].col0;
        int col1 = aux[i].col1;
        int col2 = aux[i].col2;
        int col3 = aux[i].col3;

        if(col0 >= first_col && col0 <= last_col)
          flag0 = 1;
        else
          flag0 = 0;

        if(col1 >= first_col && col1 <= last_col)
          flag1 = 1;
        else
          flag1 = 0;

        if(col2 >= first_col && col2 <= last_col)
          flag2 = 1;
        else
          flag2 = 0;

        if(col3 >= first_col && col3 <= last_col)
          flag3 = 1;
        else
          flag3 = 0;

        if(flag0 || flag1 || flag2 || flag3)
          {
            double mass = partin[i].Mass;

            double dx = aux[i].dx;
            double dy = aux[i].dy;

            int slab_z = (int)(to_slab_fac * (partin[i].Pos[2] - All.Corner[grnr][2]));
            int slab_zz = slab_z + 1;

            double dz = to_slab_fac * (partin[i].Pos[2] - All.Corner[grnr][2]) - slab_z;

            if(flag0)
              {
                rhogrid[FC(col0, slab_z)] += (mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz));
                rhogrid[FC(col0, slab_zz)] += (mass * (1.0 - dx) * (1.0 - dy) * (dz));
              }

            if(flag1)
              {
                rhogrid[FC(col1, slab_z)] += (mass * (1.0 - dx) * (dy) * (1.0 - dz));
                rhogrid[FC(col1, slab_zz)] += (mass * (1.0 - dx) * (dy) * (dz));
              }

            if(flag2)
              {
                rhogrid[FC(col2, slab_z)] += (mass * (dx) * (1.0 - dy) * (1.0 - dz));
                rhogrid[FC(col2, slab_zz)] += (mass * (dx) * (1.0 - dy) * (dz));
              }

            if(flag3)
              {
                rhogrid[FC(col3, slab_z)] += (mass * (dx) * (dy) * (1.0 - dz));
                rhogrid[FC(col3, slab_zz)] += (mass * (dx) * (dy) * (dz));
              }
          }
      }
  }

  myfree(aux);

#endif /* #ifndef FFT_COLUMN_BASED #else */
}

/*! \brief If dim<0, this function reads out the potential, otherwise
 *         Cartesian force components.
 *
 *  \param[in] grnr Grid number (0: base grid, 1: high res grid).
 *  \param[in] dim Dimension of component to be read out (< 0: potential).
 *
 *  \return void
 */
void pmforce_nonperiodic_uniform_optimized_readout_forces_or_potential(int grnr, int dim)
{
#ifdef EVALPOTENTIAL
  /* factor to get potential */
  double fac = All.G / pow(All.TotalMeshSize[grnr], 4) * pow(All.TotalMeshSize[grnr] / GRID, 3);
#endif /* #ifdef EVALPOTENTIAL */

  double to_slab_fac = GRID / All.TotalMeshSize[grnr];

  double *flistin  = (double *)mymalloc("flistin", nimport * sizeof(double));
  double *flistout = (double *)mymalloc("flistout", nexport * sizeof(double));

  fft_real *grid;

  if(dim < 0)
    grid = rhogrid;
  else
    grid = forcegrid;

  size_t i;
  for(i = 0; i < nimport; i++)
    {
      flistin[i] = 0;

      int slab_x = (int)(to_slab_fac * (partin[i].Pos[0] - All.Corner[grnr][0]));
      int slab_y = (int)(to_slab_fac * (partin[i].Pos[1] - All.Corner[grnr][1]));
      int slab_z = (int)(to_slab_fac * (partin[i].Pos[2] - All.Corner[grnr][2]));

      double dx = to_slab_fac * (partin[i].Pos[0] - All.Corner[grnr][0]) - slab_x;
      double dy = to_slab_fac * (partin[i].Pos[1] - All.Corner[grnr][1]) - slab_y;
      double dz = to_slab_fac * (partin[i].Pos[2] - All.Corner[grnr][2]) - slab_z;

      int slab_xx = slab_x + 1;
      int slab_yy = slab_y + 1;
      int slab_zz = slab_z + 1;

#ifndef FFT_COLUMN_BASED
      if(myplan.slab_to_task[slab_x] == ThisTask)
        {
          slab_x -= myplan.first_slab_x_of_task[ThisTask];

          flistin[i] += +grid[FI(slab_x, slab_y, slab_z)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FI(slab_x, slab_y, slab_zz)] * (1.0 - dx) * (1.0 - dy) * (dz) +
                        grid[FI(slab_x, slab_yy, slab_z)] * (1.0 - dx) * (dy) * (1.0 - dz) +
                        grid[FI(slab_x, slab_yy, slab_zz)] * (1.0 - dx) * (dy) * (dz);
        }

      if(myplan.slab_to_task[slab_xx] == ThisTask)
        {
          slab_xx -= myplan.first_slab_x_of_task[ThisTask];

          flistin[i] += +grid[FI(slab_xx, slab_y, slab_z)] * (dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FI(slab_xx, slab_y, slab_zz)] * (dx) * (1.0 - dy) * (dz) +
                        grid[FI(slab_xx, slab_yy, slab_z)] * (dx) * (dy) * (1.0 - dz) +
                        grid[FI(slab_xx, slab_yy, slab_zz)] * (dx) * (dy) * (dz);
        }
#else  /* #ifndef FFT_COLUMN_BASED */
      int column0 = slab_x * GRID + slab_y;
      int column1 = slab_x * GRID + slab_yy;
      int column2 = slab_xx * GRID + slab_y;
      int column3 = slab_xx * GRID + slab_yy;

      if(column0 >= myplan.base_firstcol && column0 <= myplan.base_lastcol)
        {
          flistin[i] += +grid[FC(column0, slab_z)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FC(column0, slab_zz)] * (1.0 - dx) * (1.0 - dy) * (dz);
        }
      if(column1 >= myplan.base_firstcol && column1 <= myplan.base_lastcol)
        {
          flistin[i] +=
              +grid[FC(column1, slab_z)] * (1.0 - dx) * (dy) * (1.0 - dz) + grid[FC(column1, slab_zz)] * (1.0 - dx) * (dy) * (dz);
        }

      if(column2 >= myplan.base_firstcol && column2 <= myplan.base_lastcol)
        {
          flistin[i] +=
              +grid[FC(column2, slab_z)] * (dx) * (1.0 - dy) * (1.0 - dz) + grid[FC(column2, slab_zz)] * (dx) * (1.0 - dy) * (dz);
        }

      if(column3 >= myplan.base_firstcol && column3 <= myplan.base_lastcol)
        {
          flistin[i] += +grid[FC(column3, slab_z)] * (dx) * (dy) * (1.0 - dz) + grid[FC(column3, slab_zz)] * (dx) * (dy) * (dz);
        }
#endif /* #ifndef FFT_COLUMN_BASED #else */
    }

  /* exchange the potential component data */
  int flag_big = 0, flag_big_all;
  for(i = 0; i < NTask; i++)
    if(Sndpm_count[i] * sizeof(double) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
      flag_big = 1;

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  /* exchange  data */
  myMPI_Alltoallv(flistin, Rcvpm_count, Rcvpm_offset, flistout, Sndpm_count, Sndpm_offset, sizeof(double), flag_big_all,
                  MPI_COMM_WORLD);

  /* now assign them to the correct particles */
  int multiNtask = roundup_to_multiple_of_cacheline_size(NTask * sizeof(size_t)) / sizeof(size_t);

  {
    size_t *send_count  = Sndpm_count + get_thread_num() * multiNtask;
    size_t *send_offset = Sndpm_offset + get_thread_num() * multiNtask;

    int j;
    for(j = 0; j < NTask; j++)
      send_count[j] = 0;

    int i;
    for(i = 0; i < NumPart; i++)
      {
        MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
        if(P[i].Type == 0)
          pos = SphP[i].Center;
        else
#endif /* #ifdef CELL_CENTER_GRAVITY */
          pos = P[i].Pos;

        if(pos[0] < All.Corner[grnr][0] || pos[0] >= All.UpperCorner[grnr][0])
          continue;
        if(pos[1] < All.Corner[grnr][1] || pos[1] >= All.UpperCorner[grnr][1])
          continue;
        if(pos[2] < All.Corner[grnr][2] || pos[2] >= All.UpperCorner[grnr][2])
          continue;

        int slab_x  = (int)(to_slab_fac * (pos[0] - All.Corner[grnr][0]));
        int slab_xx = slab_x + 1;

#ifndef FFT_COLUMN_BASED
        int task0   = myplan.slab_to_task[slab_x];
        int task1   = myplan.slab_to_task[slab_xx];

        double value = flistout[send_offset[task0] + send_count[task0]++];

        if(task0 != task1)
          value += flistout[send_offset[task1] + send_count[task1]++];
#else  /* #ifndef FFT_COLUMN_BASED */
        int slab_y = (int)(to_slab_fac * (pos[1] - All.Corner[grnr][1]));
        int slab_yy = slab_y + 1;

        int column0 = slab_x * GRID + slab_y;
        int column1 = slab_x * GRID + slab_yy;
        int column2 = slab_xx * GRID + slab_y;
        int column3 = slab_xx * GRID + slab_yy;

        int task0, task1, task2, task3;

        if(column0 < myplan.pivotcol)
          task0 = column0 / myplan.avg;
        else
          task0 = (column0 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column1 < myplan.pivotcol)
          task1 = column1 / myplan.avg;
        else
          task1 = (column1 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column2 < myplan.pivotcol)
          task2 = column2 / myplan.avg;
        else
          task2 = (column2 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        if(column3 < myplan.pivotcol)
          task3 = column3 / myplan.avg;
        else
          task3 = (column3 - myplan.pivotcol) / (myplan.avg - 1) + myplan.tasklastsection;

        double value = flistout[send_offset[task0] + send_count[task0]++];

        if(task1 != task0)
          value += flistout[send_offset[task1] + send_count[task1]++];

        if(task2 != task1 && task2 != task0)
          value += flistout[send_offset[task2] + send_count[task2]++];

        if(task3 != task0 && task3 != task1 && task3 != task2)
          value += flistout[send_offset[task3] + send_count[task3]++];
#endif /* #ifndef FFT_COLUMN_BASED */

#ifdef PLACEHIGHRESREGION
        if(grnr == 1)
          if(!(pmforce_is_particle_high_res(P[i].Type, pos)))
            continue;
#endif /* #ifdef PLACEHIGHRESREGION */

        if(dim < 0)
          {
#ifdef EVALPOTENTIAL
            P[i].PM_Potential += value * fac;
#endif /* #ifdef EVALPOTENTIAL */
          }
        else
          P[i].GravPM[dim] += value;
      }
  }

  int j;
  /* restore total Sndpm_count */
  for(j = 1; j < MaxThreads; j++)
    for(i = 0; i < NTask; i++)
      Sndpm_count[i] += Sndpm_count[i + j * multiNtask];

  myfree(flistout);
  myfree(flistin);
}
#endif /* #ifdef PM_ZOOM_OPTIMIZED #else */

/*! \brief Calculates the long-range non-periodic forces using the PM method.
 *
 *  The potential is Gaussian filtered with Asmth, given in mesh-cell units.
 *  The potential is finite differenced using a 4-point finite differencing
 *  formula to obtain the force fields, which are then interpolated to the
 *  particle positions. We carry out a CIC charge assignment, and compute the
 *  potenial by Fourier transform methods. The CIC kernel is deconvolved.
 *
 *  \param[in] grnr Grid number (0: base grid, 1 high res grid).
 *
 *  \return 0
 */
int pmforce_nonperiodic(int grnr)
{
  int i, j, flag, flagsum, dim;

  double tstart = second();

  mpi_printf("PM-NONPERIODIC: Starting non-periodic PM calculation (grid=%d)  presently allocated=%g MB).\n", grnr,
             AllocatedBytes / (1024.0 * 1024.0));

#ifndef NUMPART_PER_TASK_LARGE
  if((((long long)NumPart) << 3) >= (((long long)1) << 31))
    terminate("We are dealing with a too large particle number per MPI rank - enabling NUMPART_PER_TASK_LARGE might help.");
#endif /* #ifndef NUMPART_PER_TASK_LARGE */

  double fac = All.G / pow(All.TotalMeshSize[grnr], 4) * pow(All.TotalMeshSize[grnr] / GRID, 3); /* to get potential */
  fac *= 1 / (2 * All.TotalMeshSize[grnr] / GRID);                                               /* for finite differencing */

  /* first, check whether all particles lie in the allowed region */
  for(i = 0, flag = 0; i < NumPart; i++)
    {
      MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
      if(P[i].Type == 0)
        pos = SphP[i].Center;
      else
#endif /* #ifdef CELL_CENTER_GRAVITY */
        pos = P[i].Pos;

#ifdef PLACEHIGHRESREGION
      if(grnr == 0 || (grnr == 1 && pmforce_is_particle_high_res(P[i].Type, pos)))
#endif /* #ifdef PLACEHIGHRESREGION */
        {
          for(j = 0; j < 3; j++)
            {
              if(pos[j] < All.Xmintot[grnr][j] || pos[j] > All.Xmaxtot[grnr][j])
                {
                  if(flag == 0)
                    {
                      printf("Particle Id=%llu on task=%d with coordinates (%g|%g|%g) lies outside PM mesh.\n",
                             (unsigned long long)P[i].ID, ThisTask, pos[0], pos[1], pos[2]);
                      myflush(stdout);
                    }
                  flag++;
                  break;
                }
            }
        }
    }

  MPI_Allreduce(&flag, &flagsum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(flagsum > 0)
    {
      mpi_printf("PM-NONPERIODIC: In total %d particles were outside allowed range.\n", flagsum);
      return 1; /* error - need to return because particles were outside allowed range */
    }

#ifdef PM_ZOOM_OPTIMIZED
  pmforce_nonperiodic_zoom_optimized_prepare_density(grnr);
#else  /* #ifdef PM_ZOOM_OPTIMIZED */
  pmforce_nonperiodic_uniform_optimized_prepare_density(grnr);
#endif /* #ifdef PM_ZOOM_OPTIMIZED #else */

  /* allocate the memory to hold the FFT fields */
  forcegrid = (fft_real *)mymalloc("forcegrid", maxfftsize * sizeof(fft_real));

  workspace = forcegrid;

#ifndef FFT_COLUMN_BASED
  fft_of_rhogrid = (fft_complex *)&rhogrid[0];
#else  /* #ifndef FFT_COLUMN_BASED */
  fft_of_rhogrid = (fft_complex *)&workspace[0];
#endif /* #ifndef FFT_COLUMN_BASED #else */

  /* Do the FFT of the density field */
#ifndef FFT_COLUMN_BASED
  my_slab_based_fft(&myplan, &rhogrid[0], &workspace[0], 1);
#else  /* #ifndef FFT_COLUMN_BASED */
  my_column_based_fft(&myplan, rhogrid, workspace, 1); /* result is in workspace, not in rhogrid ! */
#endif /* #ifndef FFT_COLUMN_BASED #else */

  /* multiply with kernel in Fourier space */
  /* multiply with the Fourier transform of the Green's function (kernel) */
  /* multiply with Green's function in order to obtain the potential */

#ifdef FFT_COLUMN_BASED
  for(large_array_offset ip = 0; ip < myplan.second_transposed_ncells; ip++)
    {
#else  /* #ifdef FFT_COLUMN_BASED */
  for(int x = 0; x < GRID; x++)
    for(int y = myplan.slabstart_y; y < myplan.slabstart_y + myplan.nslab_y; y++)
      for(int z = 0; z < GRIDz; z++)
        {
#endif /* #ifdef FFT_COLUMN_BASED #else */

#ifndef FFT_COLUMN_BASED
      large_array_offset ip = ((large_array_offset)GRIDz) * (GRID * (y - myplan.slabstart_y) + x) + z;
#endif /* #ifndef FFT_COLUMN_BASED */

      double re = fft_of_rhogrid[ip][0] * fft_of_kernel[grnr][ip][0] - fft_of_rhogrid[ip][1] * fft_of_kernel[grnr][ip][1];
      double im = fft_of_rhogrid[ip][0] * fft_of_kernel[grnr][ip][1] + fft_of_rhogrid[ip][1] * fft_of_kernel[grnr][ip][0];

      fft_of_rhogrid[ip][0] = re;
      fft_of_rhogrid[ip][1] = im;
    }

    /* Do the inverse FFT to get the potential */

#ifndef FFT_COLUMN_BASED
  my_slab_based_fft(&myplan, rhogrid, workspace, -1);
#else  /* #ifndef FFT_COLUMN_BASED */
  my_column_based_fft(&myplan, workspace, rhogrid, -1);
#endif /* #ifndef FFT_COLUMN_BASED #else */

  /* Now rhogrid holds the potential */

#ifdef EVALPOTENTIAL
#ifdef PM_ZOOM_OPTIMIZED
  pmforce_nonperiodic_zoom_optimized_readout_forces_or_potential(grnr, -1);
#else  /* #ifdef PM_ZOOM_OPTIMIZED */
  pmforce_nonperiodic_uniform_optimized_readout_forces_or_potential(grnr, -1);
#endif /* #ifdef PM_ZOOM_OPTIMIZED #else */
#endif /* #ifdef EVALPOTENTIAL */

  /* get the force components by finite differencing of the potential for each dimension,
   * and send the results back to the right CPUs
   */
  for(dim = 2; dim >= 0; dim--) /* Calculate each component of the force. */
    {
      /* we do the x component last, because for differencing the potential in the x-direction, we need to construct the transpose */
#ifndef FFT_COLUMN_BASED
      if(dim == 0)
        my_slab_transposeA(&myplan, rhogrid, forcegrid); /* compute the transpose of the potential field for finite differencing */

      for(int y = 2; y < GRID / 2 - 2; y++)
        for(int x = 0; x < myplan.nslab_x; x++)
          if(x + myplan.slabstart_x >= 2 && x + myplan.slabstart_x < GRID / 2 - 2)
            for(int z = 2; z < GRID / 2 - 2; z++)
              {
                int yrr = y, yll = y, yr = y, yl = y;
                int zrr = z, zll = z, zr = z, zl = z;

                switch(dim)
                  {
                    case 0: /* note: for the x-direction, we difference the transposed direction (y) */
                    case 1:
                      yr  = y + 1;
                      yl  = y - 1;
                      yrr = y + 2;
                      yll = y - 2;

                      break;
                    case 2:
                      zr  = z + 1;
                      zl  = z - 1;
                      zrr = z + 2;
                      zll = z - 2;

                      break;
                  }

                if(dim == 0)
                  forcegrid[TI(x, y, z)] = fac * ((4.0 / 3) * (rhogrid[TI(x, yl, zl)] - rhogrid[TI(x, yr, zr)]) -
                                                  (1.0 / 6) * (rhogrid[TI(x, yll, zll)] - rhogrid[TI(x, yrr, zrr)]));
                else
                  forcegrid[FI(x, y, z)] = fac * ((4.0 / 3) * (rhogrid[FI(x, yl, zl)] - rhogrid[FI(x, yr, zr)]) -
                                                  (1.0 / 6) * (rhogrid[FI(x, yll, zll)] - rhogrid[FI(x, yrr, zrr)]));
              }

      if(dim == 0)
        my_slab_transposeB(&myplan, forcegrid, rhogrid); /* reverse the transpose from above */
#else                                                    /* #ifndef FFT_COLUMN_BASED */
      fft_real *scratch = NULL, *forcep, *potp;

      if(dim != 2)
        {
          scratch = mymalloc("scratch", myplan.fftsize * sizeof(fft_real)); /* need a third field as scratch space */
          memcpy(scratch, rhogrid, myplan.fftsize * sizeof(fft_real));

          if(dim == 1)
            my_fft_swap23(&myplan, scratch, forcegrid);
          else
            my_fft_swap13(&myplan, scratch, forcegrid);
        }

      int ncols;
      if(dim == 2)
        ncols = myplan.base_ncol;
      else if(dim == 1)
        ncols = myplan.ncol_XZ;
      else
        ncols = myplan.ncol_YZ;

      large_array_offset i;

      for(i = 0; i < ncols; i++)
        {
          if(dim != 2)
            {
              forcep = &scratch[GRID * i];
              potp   = &forcegrid[GRID * i];
            }
          else
            {
              forcep = &forcegrid[GRID2 * i];
              potp   = &rhogrid[GRID2 * i];
            }

          int z;
          for(z = 2; z < GRID / 2 - 2; z++)
            {
              int zr  = z + 1;
              int zl  = z - 1;
              int zrr = z + 2;
              int zll = z - 2;

              forcep[z] = fac * ((4.0 / 3) * (potp[zl] - potp[zr]) - (1.0 / 6) * (potp[zll] - potp[zrr]));
            }
        }

      if(dim != 2)
        {
          if(dim == 1)
            my_fft_swap23back(&myplan, scratch, forcegrid);
          else
            my_fft_swap13back(&myplan, scratch, forcegrid);

          myfree(scratch);
        }
#endif                                                   /* #ifndef FFT_COLUMN_BASED #else */

#ifdef PM_ZOOM_OPTIMIZED
      pmforce_nonperiodic_zoom_optimized_readout_forces_or_potential(grnr, dim);
#else  /* #ifdef PM_ZOOM_OPTIMIZED */
      pmforce_nonperiodic_uniform_optimized_readout_forces_or_potential(grnr, dim);
#endif /* #ifdef PM_ZOOM_OPTIMIZED #else */
    }

  /* free stuff */
  myfree(forcegrid);
  myfree(rhogrid);

#ifdef PM_ZOOM_OPTIMIZED
  myfree(localfield_recvcount);
  myfree(localfield_offset);
  myfree(localfield_sendcount);
  myfree(localfield_first);
  myfree(localfield_data);
  myfree(localfield_globalindex);
  myfree(part);
#else  /* #ifdef PM_ZOOM_OPTIMIZED */
  myfree(partin);
  myfree(Rcvpm_offset);
  myfree(Rcvpm_count);
  myfree(Sndpm_offset);
  myfree(Sndpm_count);
#endif /* #ifdef PM_ZOOM_OPTIMIZED #else */

  double tend = second();

  mpi_printf("PM-NONPERIODIC: done.  (took %g seconds)\n", timediff(tstart, tend));

  return 0;
}

/*! \brief Sets-up the Greens function for the non-periodic potential in real
 *         space, and then converts it to Fourier space by means of an FFT.
 *
 *  \return void
 */
void pm_setup_nonperiodic_kernel(void)
{
  int i, j, k, x, y, z;
  double xx, yy, zz, r, u, fac;

  mpi_printf("PM-NONPERIODIC: Setting up non-periodic PM kernel (GRID=%d)  presently allocated=%g MB).\n", (int)GRID,
             AllocatedBytes / (1024.0 * 1024.0));

  /* now set up kernel and its Fourier transform */

#if defined(GRAVITY_NOT_PERIODIC)
  for(i = 0; i < maxfftsize; i++) /* clear local field */
    kernel[0][i] = 0;

#ifndef FFT_COLUMN_BASED
  for(i = myplan.slabstart_x; i < (myplan.slabstart_x + myplan.nslab_x); i++)
    for(j = 0; j < GRID; j++)
      {
#else  /* #ifndef FFT_COLUMN_BASED */
  int c;
  for(c = myplan.base_firstcol; c < (myplan.base_firstcol + myplan.base_ncol); c++)
    {
      i = c / GRID;
      j = c % GRID;
#endif /* #ifndef FFT_COLUMN_BASED #else */
        for(k = 0; k < GRID; k++)
          {
            xx = ((double)i) / GRID;
            yy = ((double)j) / GRID;
            zz = ((double)k) / GRID;

            if(xx >= 0.5)
              xx -= 1.0;
            if(yy >= 0.5)
              yy -= 1.0;
            if(zz >= 0.5)
              zz -= 1.0;

            r = sqrt(xx * xx + yy * yy + zz * zz);

            u = 0.5 * r / (((double)ASMTH) / GRID);

            fac = 1 - erfc(u);

#ifndef FFT_COLUMN_BASED
            size_t ip = FI(i - myplan.slabstart_x, j, k);
#else  /* #ifndef FFT_COLUMN_BASED */
          size_t ip = FC(c, k);
#endif /* #ifndef FFT_COLUMN_BASED #else */
            if(r > 0)
              kernel[0][ip] = -fac / r;
            else
              kernel[0][ip] = -1 / (sqrt(M_PI) * (((double)ASMTH) / GRID));
          }
      }

  {
    fft_real *workspc = (fft_real *)mymalloc("workspc", maxfftsize * sizeof(fft_real));
    /* Do the FFT of the kernel */
#ifndef FFT_COLUMN_BASED
    my_slab_based_fft(&myplan, kernel[0], workspc, 1);
#else  /* #ifndef FFT_COLUMN_BASED */
    my_column_based_fft(&myplan, kernel[0], workspc, 1); /* result is in workspace, not in kernel */
    memcpy(kernel[0], workspc, maxfftsize * sizeof(fft_real));
#endif /* #ifndef FFT_COLUMN_BASED #else */
    myfree(workspc);
  }

#endif /* #if defined(GRAVITY_NOT_PERIODIC) */

#if defined(PLACEHIGHRESREGION)

  for(i = 0; i < maxfftsize; i++) /* clear local field */
    kernel[1][i] = 0;

#ifndef FFT_COLUMN_BASED
  for(i = myplan.slabstart_x; i < (myplan.slabstart_x + myplan.nslab_x); i++)
    for(j = 0; j < GRID; j++)
      {
#else  /* #ifndef FFT_COLUMN_BASED */
  int c;
  for(c = myplan.base_firstcol; c < (myplan.base_firstcol + myplan.base_ncol); c++)
    {
      i = c / GRID;
      j = c % GRID;
#endif /* #ifndef FFT_COLUMN_BASED #else */
        for(k = 0; k < GRID; k++)
          {
            xx = ((double)i) / GRID;
            yy = ((double)j) / GRID;
            zz = ((double)k) / GRID;

            if(xx >= 0.5)
              xx -= 1.0;
            if(yy >= 0.5)
              yy -= 1.0;
            if(zz >= 0.5)
              zz -= 1.0;

            r = sqrt(xx * xx + yy * yy + zz * zz);

            u = 0.5 * r / (((double)ASMTH) / GRID);

            fac = erfc(u * All.Asmth[1] / All.Asmth[0]) - erfc(u);

#ifndef FFT_COLUMN_BASED
            size_t ip = FI(i - myplan.slabstart_x, j, k);
#else  /* #ifndef FFT_COLUMN_BASED */
          size_t ip = FC(c, k);
#endif /* #ifndef FFT_COLUMN_BASED #else */

            if(r > 0)
              kernel[1][ip] = -fac / r;
            else
              {
                fac           = 1 - All.Asmth[1] / All.Asmth[0];
                kernel[1][ip] = -fac / (sqrt(M_PI) * (((double)ASMTH) / GRID));
              }
          }
      }

  {
    fft_real *workspc = (fft_real *)mymalloc("workspc", maxfftsize * sizeof(fft_real));
    /* Do the FFT of the kernel */
#ifndef FFT_COLUMN_BASED
    my_slab_based_fft(&myplan, kernel[1], workspc, 1);
#else  /* #ifndef FFT_COLUMN_BASED */
    my_column_based_fft(&myplan, kernel[1], workspc, 1); /* result is in workspace, not in kernel */
    memcpy(kernel[1], workspc, maxfftsize * sizeof(fft_real));
#endif /* #ifndef FFT_COLUMN_BASED #else */
    myfree(workspc);
  }

#endif /* #if defined(PLACEHIGHRESREGION) */

  /* deconvolve the Greens function twice with the CIC kernel */
#ifdef FFT_COLUMN_BASED

  large_array_offset ip, ipcell;

  for(ip = 0; ip < myplan.second_transposed_ncells; ip++)
    {
      ipcell = ip + myplan.transposed_firstcol * GRID;
      y      = ipcell / (GRID * GRIDz);
      int yr = ipcell % (GRID * GRIDz);
      z      = yr / GRID;
      x      = yr % GRID;
#else  /* #ifdef FFT_COLUMN_BASED */
  for(x = 0; x < GRID; x++)
    for(y = myplan.slabstart_y; y < myplan.slabstart_y + myplan.nslab_y; y++)
      for(z = 0; z < GRIDz; z++)
        {
#endif /* #ifdef FFT_COLUMN_BASED #else */

      double kx, ky, kz;

      if(x > GRID / 2)
        kx = x - GRID;
      else
        kx = x;
      if(y > GRID / 2)
        ky = y - GRID;
      else
        ky = y;
      if(z > GRID / 2)
        kz = z - GRID;
      else
        kz = z;

      double k2 = kx * kx + ky * ky + kz * kz;

      if(k2 > 0)
        {
          double fx = 1, fy = 1, fz = 1;

          if(kx != 0)
            {
              fx = (M_PI * kx) / GRID;
              fx = sin(fx) / fx;
            }
          if(ky != 0)
            {
              fy = (M_PI * ky) / GRID;
              fy = sin(fy) / fy;
            }
          if(kz != 0)
            {
              fz = (M_PI * kz) / GRID;
              fz = sin(fz) / fz;
            }

          double ff = 1 / (fx * fy * fz);
          ff        = ff * ff * ff * ff;

#ifndef FFT_COLUMN_BASED
          large_array_offset ip = ((large_array_offset)GRIDz) * (GRID * (y - myplan.slabstart_y) + x) + z;
#endif /* #ifndef FFT_COLUMN_BASED */
#if defined(GRAVITY_NOT_PERIODIC)
          fft_of_kernel[0][ip][0] *= ff;
          fft_of_kernel[0][ip][1] *= ff;
#endif /* #if defined(GRAVITY_NOT_PERIODIC) */
#if defined(PLACEHIGHRESREGION)
          fft_of_kernel[1][ip][0] *= ff;
          fft_of_kernel[1][ip][1] *= ff;
#endif /* #if defined(PLACEHIGHRESREGION) */
        }
    }

  /* end deconvolution */
}

#ifdef PM_ZOOM_OPTIMIZED

/*! \brief Sort function for 'part' array indices.
 *
 *  Sorts the indices into the 'part' array by the global index of the
 *  corresponding 'part_slab_data' struct.
 *
 *  \param[in] a index to be compared.
 *  \param[in] b index to be compared.
 *
 *  \return sort result
 */
static int pm_periodic_compare_sortindex(const void *a, const void *b)
{
  if(part[*(int *)a].globalindex < part[*(int *)b].globalindex)
    return -1;

  if(part[*(int *)a].globalindex > part[*(int *)b].globalindex)
    return +1;

  return 0;
}

/*! \brief Implements the sorting function for mysort_pmperiodic()
 *
 *  The index array is sorted using a merge sort algorithm.
 *
 *  \param[in, out] b Index array to sort.
 *  \param[in] n Number of elements to sort.
 *  \param[out] t Temporary buffer array.
 *
 *  \return void
 */
static void msort_pmperiodic_with_tmp(large_numpart_type *b, size_t n, large_numpart_type *t)
{
  large_numpart_type *tmp;
  large_numpart_type *b1, *b2;
  size_t n1, n2;

  if(n <= 1)
    return;

  n1 = n / 2;
  n2 = n - n1;
  b1 = b;
  b2 = b + n1;

  msort_pmperiodic_with_tmp(b1, n1, t);
  msort_pmperiodic_with_tmp(b2, n2, t);

  tmp = t;

  while(n1 > 0 && n2 > 0)
    {
      if(part[*b1].globalindex <= part[*b2].globalindex)
        {
          --n1;
          *tmp++ = *b1++;
        }
      else
        {
          --n2;
          *tmp++ = *b2++;
        }
    }

  if(n1 > 0)
    memcpy(tmp, b1, n1 * sizeof(large_numpart_type));

  memcpy(b, t, (n - n2) * sizeof(large_numpart_type));
}

/*! \brief Sorts the index array b of n entries using the sort kernel
 *         cmp.
 *
 *  The parameter s is set to sizeof(int). The index array b
 *  is sorted according to the globalindex field of the referenced item in the
 *  'part' array
 *
 *  \param[in, out] b The index array to sort.
 *  \param[in] n Number of entries in array b.
 *  \param[in] s Size of each entry (must be sizeof(int)).
 *  \param[in] cmp Comparison function.
 *
 *  \return void
 */
static void mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *))
{
  const size_t size = n * s;

  large_numpart_type *tmp = (large_numpart_type *)mymalloc("tmp", size);

  msort_pmperiodic_with_tmp((large_numpart_type *)b, n, tmp);

  myfree(tmp);
}
#endif /* #ifdef PM_ZOOM_OPTIMIZED */

#endif /* #if defined(PMGRID) && (defined(PLACEHIGHRESREGION) || defined(GRAVITY_NOT_PERIODIC)) */
