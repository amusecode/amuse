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
 * \file        src/gravity/pm/pm_periodic.c
 * \date        05/2018
 * \brief       Routines for periodic PM-force computation.
 * \details     These routines support two different strategies for doing the
 *              particle data exchange to assemble the density field and to
 *              read out the forces and potentials:
 *
 *              The default scheme sends the particle positions to the target
 *              slabs, and bins them there. This works usually well for
 *              homogeneously loaded boxes, but can be problematic for zoom-in
 *              runs. In the latter case, PM_ZOOM_OPTIMIZED can be activated,
 *              where the data is binned on the originating processor followed
 *              by assembly of the binned density field.
 *
 *              In addition, the routines can be either used with a slab-based
 *              FFT (as is traditionally done in FFTW), or with a column-based
 *              FFT. The latter requires more communication and is hence
 *              usually slower than the slab-based one. But if the number of
 *              MPI ranks exceeds the number of cells per dimension, then the
 *              column-based one can still scale and offers a balanced memory
 *              consumption, whereas this is not the case for the slab-based
 *              approach. To select the column-based FFT, the switch
 *              FFT_COLUMN_BASED can be activated.
 *
 *              The switches PM_ZOOM_OPTIMIZED and FFT_COLUMN_BASED may also
 *              be combined, such that there are 4 main modes of how the PM
 *              routines may operate.
 *
 *              It is also possible to use non-cubical boxes, by means of
 *              setting one or several of the LONG_X, LONG_Y, and LONG_Z
 *              options in the config file. The values need to be integers,
 *              and then BoxSize is stretched by that factor in the
 *              corresponding dimension.
 *
 *              Much of the code is multi-threaded, so there should be some
 *              speed-up if OpenMP is used with NUM_THREADS > 1, but the
 *              benefit may be limited because the data transfer steps (which
 *              weigh in quite heavily) are not accelerated by this.
 *
 *              If eight times the particle load per processor exceeds 2^31
 *              ~ 2 billion, one should activate NUMPART_PER_TASK_LARGE. The
 *              code will check this condition and terminate if this is
 *              violated, so there should hopefully be no severe risk to
 *              accidentally forget this.
 *
 *              contains functions:
 *                void pm_init_periodic(void)
 *                void pmforce_zoom_optimized_prepare_density(int mode, int
 *                  *typelist)
 *                void pmforce_zoom_optimized_readout_forces_or_potential(int
 *                  dim)
 *                static void pmforce_uniform_optimized_prepare_density(int
 *                  mode)
 *                static void pmforce_uniform_optimized_readout_forces_or_
 *                  potential(int dim)
 *                void pmforce_periodic(int mode, int *typelist)
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

#if defined(PMGRID)

#define GRIDX (PMGRID * STRETCHX * DBX + DBX_EXTRA)
#define GRIDY (PMGRID * STRETCHY * DBY + DBY_EXTRA)
#define GRIDZ (PMGRID * STRETCHZ * DBZ + DBZ_EXTRA)

#define GRIDz (GRIDZ / 2 + 1)
#define GRID2 (2 * GRIDz)

#if(GRIDX > 1024) || (GRIDY > 1024) || (GRIDZ > 1024)
typedef long long large_array_offset; /* use a larger data type in this case so that we can always address all cells of the 3D grid
                                         with a single index */
#else                                 /* #if (GRIDX > 1024) || (GRIDY > 1024) || (GRIDZ > 1024) */
typedef unsigned int large_array_offset;
#endif                                /* #if (GRIDX > 1024) || (GRIDY > 1024) || (GRIDZ > 1024) #else */

#ifdef NUMPART_PER_TASK_LARGE
typedef long long large_numpart_type; /* if there is a risk that the local particle number times 8 overflows a 32-bit integer, this
                                         data type should be used */
#else                                 /* #ifdef NUMPART_PER_TASK_LARGE */
typedef int large_numpart_type;
#endif                                /* #ifdef NUMPART_PER_TASK_LARGE #else */

/* short-cut macros for accessing different 3D arrays */
#define FI(x, y, z) (((large_array_offset)GRID2) * (GRIDY * (x) + (y)) + (z))
#define FC(c, z) (((large_array_offset)GRID2) * ((c)-myplan.base_firstcol) + (z))
#ifndef FFT_COLUMN_BASED
#define NI(x, y, z) (((large_array_offset)GRIDZ) * ((y) + (x)*myplan.nslab_y) + (z))
#endif /* #ifndef FFT_COLUMN_BASED */

/* variables for power spectrum estimation */
#ifndef BINS_PS
#define BINS_PS 2000 /* number of bins for power spectrum computation */
#endif               /* #ifndef BINS_PS */
#ifndef POWERSPEC_FOLDFAC
#define POWERSPEC_FOLDFAC 16. /* folding factor to obtain an estimate of the power spectrum on very small scales */
#endif                        /* #ifndef POWERSPEC_FOLDFAC */

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

/*! \brief Array containing the FFT of #rhogrid
 *
 *  This pointer points to the same array as #rhogrid,
 *  because in-place FFTs are used.
 */
static fft_complex *fft_of_rhogrid;

/* Variable for power spectrum calculation */
static double power_spec_totmass, power_spec_totmass2;
static long long power_spec_totnumpart;

/*! \brief This routine generates the FFT-plans to carry out the FFTs later on.
 *
 *  Some auxiliary variables for bookkeeping are also initialized.
 *
 *  \return void
 */
void pm_init_periodic(void)
{
#ifdef LONG_X
  if(LONG_X != (int)(LONG_X))
    terminate("LONG_X must be an integer if used with PMGRID");
#endif /* #ifdef LONG_X */

#ifdef LONG_Y
  if(LONG_Y != (int)(LONG_Y))
    terminate("LONG_Y must be an integer if used with PMGRID");
#endif /* #ifdef LONG_Y */

#ifdef LONG_Z
  if(LONG_Z != (int)(LONG_Z))
    terminate("LONG_Z must be an integer if used with PMGRID");
#endif /* #ifdef LONG_Z */

  All.Asmth[0] = ASMTH * All.BoxSize / PMGRID;
  All.Rcut[0]  = RCUT * All.Asmth[0];

  /* Set up the FFTW-3 plan files. */
  int ndimx[1] = {GRIDX}; /* dimension of the 1D transforms */
  int ndimy[1] = {GRIDY}; /* dimension of the 1D transforms */
  int ndimz[1] = {GRIDZ}; /* dimension of the 1D transforms */

  int max_GRID2 = 2 * (imax(imax(GRIDX, GRIDY), GRIDZ) / 2 + 1);

  /* temporarily allocate some arrays to make sure that out-of-place plans are created */
  rhogrid   = (fft_real *)mymalloc("rhogrid", max_GRID2 * sizeof(fft_real));
  forcegrid = (fft_real *)mymalloc("forcegrid", max_GRID2 * sizeof(fft_real));

#ifdef DOUBLEPRECISION_FFTW
  int alignflag = 0;
#else  /* #ifdef DOUBLEPRECISION_FFTW */
  /* for single precision, the start of our FFT columns is presently only guaranteed to be 8-byte aligned */
  int alignflag = FFTW_UNALIGNED;
#endif /* #ifdef DOUBLEPRECISION_FFTW #else */

  myplan.forward_plan_zdir = FFTW(plan_many_dft_r2c)(1, ndimz, 1, rhogrid, 0, 1, GRID2, (fft_complex *)forcegrid, 0, 1, GRIDz,
                                                     FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

#ifndef FFT_COLUMN_BASED
  int stride = GRIDz;
#else  /* #ifndef FFT_COLUMN_BASED */
  int stride    = 1;
#endif /* #ifndef FFT_COLUMN_BASED #else */

  myplan.forward_plan_ydir =
      FFTW(plan_many_dft)(1, ndimy, 1, (fft_complex *)rhogrid, 0, stride, GRIDz * GRIDY, (fft_complex *)forcegrid, 0, stride,
                          GRIDz * GRIDY, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.forward_plan_xdir =
      FFTW(plan_many_dft)(1, ndimx, 1, (fft_complex *)rhogrid, 0, stride, GRIDz * GRIDX, (fft_complex *)forcegrid, 0, stride,
                          GRIDz * GRIDX, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_xdir =
      FFTW(plan_many_dft)(1, ndimx, 1, (fft_complex *)rhogrid, 0, stride, GRIDz * GRIDX, (fft_complex *)forcegrid, 0, stride,
                          GRIDz * GRIDX, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_ydir =
      FFTW(plan_many_dft)(1, ndimy, 1, (fft_complex *)rhogrid, 0, stride, GRIDz * GRIDY, (fft_complex *)forcegrid, 0, stride,
                          GRIDz * GRIDY, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myplan.backward_plan_zdir = FFTW(plan_many_dft_c2r)(1, ndimz, 1, (fft_complex *)rhogrid, 0, 1, GRIDz, forcegrid, 0, 1, GRID2,
                                                      FFTW_ESTIMATE | FFTW_DESTROY_INPUT | alignflag);

  myfree(forcegrid);
  myfree(rhogrid);

#ifndef FFT_COLUMN_BASED

  my_slab_based_fft_init(&myplan, GRIDX, GRIDY, GRIDZ);

  maxfftsize = imax(myplan.largest_x_slab * GRIDY, myplan.largest_y_slab * GRIDX) * ((size_t)GRID2);

#else /* #ifndef FFT_COLUMN_BASED */

  my_column_based_fft_init(&myplan, GRIDX, GRIDY, GRIDZ);

  maxfftsize = myplan.max_datasize;

#endif /* #ifndef FFT_COLUMN_BASED #else */
}

/* Below, the two functions
 *
 *           pmforce_ ...... _prepare_density()
 * and
 *           pmforce_ ...... _readout_forces_or_potential(int dim)
 *
 * are defined in two different versions, one that works better for uniform
 * simulations, the other for zoom-in runs. Only one of the two sets is used,
 * depending on the setting of PM_ZOOM_OPTIMIZED.
 */
#ifdef PM_ZOOM_OPTIMIZED
static void mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *));
static int pm_periodic_compare_sortindex(const void *a, const void *b);

/*! \brief This structure links the particles to the mesh cells, to which they
 *         contribute their mass.
 *
 *  Each particle will have eight items of this structure in the #part array.
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

/*! \brief Prepares density field for PM calculation in zoom-optimized
 *         algorithm.
 *
 *  \param[in] mode Modes force calculation or power spectrum calculation.
 *  \param[in] typelist Which particles to include (only for power spectrum).
 *
 *  \return void
 */
void pmforce_zoom_optimized_prepare_density(int mode, int *typelist)
{
  large_numpart_type i;
  int level, recvTask;
  MPI_Status status;

  double to_slab_fac =
      PMGRID / All.BoxSize; /* note: This is the same as GRIDX / (All.BoxSize * LONG_X), and similarly for each dimension */

  if(mode == 2)
    to_slab_fac *= POWERSPEC_FOLDFAC;
  if(mode == 3)
    to_slab_fac *= POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC;

  part                               = (struct part_slab_data *)mymalloc("part", 8 * (NumPart * sizeof(struct part_slab_data)));
  large_numpart_type *part_sortindex = (large_numpart_type *)mymalloc("part_sortindex", 8 * (NumPart * sizeof(large_numpart_type)));

  /* determine the cells each particle accesses */
  for(i = 0; i < NumPart; i++)
    {
      MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
      MyDouble posw[3], xtmp, ytmp, ztmp;
      if(P[i].Type == 0)
        {
          posw[0] = WRAP_X(SphP[i].Center[0]);
          posw[1] = WRAP_Y(SphP[i].Center[1]);
          posw[2] = WRAP_Z(SphP[i].Center[2]);

          pos = posw;
        }
      else
#endif /* #ifdef CELL_CENTER_GRAVITY */
        pos = P[i].Pos;

      int slab_x = (int)(to_slab_fac * pos[0]);
      int slab_y = (int)(to_slab_fac * pos[1]);
      int slab_z = (int)(to_slab_fac * pos[2]);

      if(mode >= 2)
        {
          slab_x %= GRIDX;
          slab_y %= GRIDY;
          slab_z %= GRIDZ;
        }
      else
        {
          if(slab_x >= GRIDX)
            slab_x -= GRIDX;
          if(slab_y >= GRIDY)
            slab_y -= GRIDY;
          if(slab_z >= GRIDZ)
            slab_z -= GRIDZ;
        }

      large_numpart_type index_on_grid = ((large_numpart_type)i) << 3;

      for(int xx = 0; xx < 2; xx++)
        for(int yy = 0; yy < 2; yy++)
          for(int zz = 0; zz < 2; zz++)
            {
              int slab_xx = slab_x + xx;
              int slab_yy = slab_y + yy;
              int slab_zz = slab_z + zz;

              if(slab_xx >= GRIDX)
                slab_xx -= GRIDX;
              if(slab_yy >= GRIDY)
                slab_yy -= GRIDY;
              if(slab_zz >= GRIDZ)
                slab_zz -= GRIDZ;

              large_array_offset offset = FI(slab_xx, slab_yy, slab_zz);

              part[index_on_grid].partindex   = (i << 3) + (xx << 2) + (yy << 1) + zz;
              part[index_on_grid].globalindex = offset;
              part_sortindex[index_on_grid]   = index_on_grid;
              index_on_grid++;
            }
    }

  /* note: num_on_grid will be  8 times larger than the particle number, but num_field_points will generally be much smaller */

  large_array_offset num_field_points;
  large_numpart_type num_on_grid = ((large_numpart_type)NumPart) << 3;

  /* bring the part-field into the order of the accessed cells. This allows the removal of duplicates */
  mysort_pmperiodic(part_sortindex, num_on_grid, sizeof(large_numpart_type), pm_periodic_compare_sortindex);

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
      int slab = part[part_sortindex[i]].globalindex / (GRIDY * GRID2);
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
      MyDouble posw[3], xtmp, ytmp, ztmp;
      if(P[pindex].Type == 0)
        {
          posw[0] = WRAP_X(SphP[pindex].Center[0]);
          posw[1] = WRAP_Y(SphP[pindex].Center[1]);
          posw[2] = WRAP_Z(SphP[pindex].Center[2]);

          pos = posw;
        }
      else
#endif /* #ifdef CELL_CENTER_GRAVITY */
        pos = P[pindex].Pos;

      int slab_x = (int)(to_slab_fac * pos[0]);
      int slab_y = (int)(to_slab_fac * pos[1]);
      int slab_z = (int)(to_slab_fac * pos[2]);

      double dx = to_slab_fac * pos[0] - slab_x;
      double dy = to_slab_fac * pos[1] - slab_y;
      double dz = to_slab_fac * pos[2] - slab_z;

      double weight = P[pindex].Mass;

      if(mode) /* only for power spectrum calculation */
        if(typelist[P[pindex].Type] == 0)
          continue;

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
                  import_globalindex[i] - myplan.first_slab_x_of_task[ThisTask] * GRIDY * ((large_array_offset)GRID2);
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

/* \brief Function to read out the force component corresponding to spatial
 *        dimension 'dim'.
 *
 *  \param[in] dim Dimension to be read out; If dim is negative, potential
 *             values are read out and assigned to particles.
 *
 *  \return void
 */
void pmforce_zoom_optimized_readout_forces_or_potential(int dim)
{
#ifdef EVALPOTENTIAL
  double fac = 4 * M_PI * All.G / (pow(All.BoxSize, 3) * STRETCHX * STRETCHY * STRETCHZ); /* to get potential  */
#endif                                                                                    /* #ifdef EVALPOTENTIAL */

  large_numpart_type i;
  int level, recvTask;
  MPI_Status status;

  fft_real *grid;

  if(dim < 0)
    grid = rhogrid;
  else
    grid = forcegrid;

  double to_slab_fac = PMGRID / All.BoxSize;

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
                  import_globalindex[i] - myplan.first_slab_x_of_task[ThisTask] * GRIDY * ((large_array_offset)GRID2);
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

  /* read out the froce/potential values, which all have been assembled in localfield_data */
  for(i = 0; i < NumPart; i++)
    {
      large_numpart_type j = (i << 3);

      MyDouble *pos;

#ifdef CELL_CENTER_GRAVITY
      MyDouble posw[3], xtmp, ytmp, ztmp;
      if(P[i].Type == 0)
        {
          posw[0] = WRAP_X(SphP[i].Center[0]);
          posw[1] = WRAP_Y(SphP[i].Center[1]);
          posw[2] = WRAP_Z(SphP[i].Center[2]);

          pos = posw;
        }
      else
#endif /* #ifdef CELL_CENTER_GRAVITY */
        pos = P[i].Pos;

      int slab_x = (int)(to_slab_fac * pos[0]);
      double dx  = to_slab_fac * pos[0] - slab_x;

      int slab_y = (int)(to_slab_fac * pos[1]);
      double dy  = to_slab_fac * pos[1] - slab_y;

      int slab_z = (int)(to_slab_fac * pos[2]);
      double dz  = to_slab_fac * pos[2] - slab_z;

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

/*
 *  Here come the routines for a different communication algorithm that is
 *  better suited for a homogenuously loaded boxes.
 */

/*! \brief Structure for particle buffer.
 */
static struct partbuf
{
  MyFloat Mass;
  MyFloat Pos[3];
} * partin, *partout;

static size_t nimport, nexport;

static size_t *Sndpm_count, *Sndpm_offset;
static size_t *Rcvpm_count, *Rcvpm_offset;

/*! \brief Prepares density field for PM calculation in uniform box optimized
 *         algorithm.
 *
 *  \param[in] mode Modes force calculation.
 *
 *  \return void
 */
static void pmforce_uniform_optimized_prepare_density(int mode)
{
  int i, j;

  double to_slab_fac = PMGRID / All.BoxSize;

  if(mode == 2)
    to_slab_fac *= POWERSPEC_FOLDFAC;
  if(mode == 3)
    to_slab_fac *= POWERSPEC_FOLDFAC * POWERSPEC_FOLDFAC;

  /* We here enlarge NTask such that each thread gets his own cache line for send_count/send_offset.
   * This should hopefully prevent a performance penalty from 'false sharing' for these variables
   */
  int multiNtask = roundup_to_multiple_of_cacheline_size(NTask * sizeof(size_t)) / sizeof(size_t);

  Sndpm_count  = (size_t *)mymalloc("Sndpm_count", MaxThreads * multiNtask * sizeof(size_t));
  Sndpm_offset = (size_t *)mymalloc("Sndpm_offset", MaxThreads * multiNtask * sizeof(size_t));
  Rcvpm_count  = (size_t *)mymalloc("Rcvpm_count", NTask * sizeof(size_t));
  Rcvpm_offset = (size_t *)mymalloc("Rcvpm_offset", NTask * sizeof(size_t));

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
        MyDouble posw[3], xtmp, ytmp, ztmp;
        if(P[i].Type == 0)
          {
            posw[0] = WRAP_X(SphP[i].Center[0]);
            posw[1] = WRAP_Y(SphP[i].Center[1]);
            posw[2] = WRAP_Z(SphP[i].Center[2]);

            pos = posw;
          }
        else
#endif /* #ifdef CELL_CENTER_GRAVITY */
          pos = P[i].Pos;

        int slab_x  = (int)(to_slab_fac * pos[0]);
        int slab_xx = slab_x + 1;

        if(mode >= 2)
          {
            slab_x %= GRIDX;
            slab_xx %= GRIDX;
          }
        else
          {
            if(slab_x >= GRIDX)
              slab_x -= GRIDX;

            if(slab_xx >= GRIDX)
              slab_xx -= GRIDX;
          }

#ifndef FFT_COLUMN_BASED
        int task0 = myplan.slab_to_task[slab_x];
        int task1 = myplan.slab_to_task[slab_xx];

        send_count[task0]++;
        if(task0 != task1)
          send_count[task1]++;
#else  /* #ifndef FFT_COLUMN_BASED */
        int slab_y  = (int)(to_slab_fac * pos[1]);
        int slab_yy = slab_y + 1;

        if(mode >= 2)
          {
            slab_y %= GRIDY;
            slab_yy %= GRIDY;
          }
        else
          {
            if(slab_y >= GRIDY)
              slab_y -= GRIDY;

            if(slab_yy >= GRIDY)
              slab_yy -= GRIDY;
          }

        int column0 = slab_x * GRIDY + slab_y;
        int column1 = slab_x * GRIDY + slab_yy;
        int column2 = slab_xx * GRIDY + slab_y;
        int column3 = slab_xx * GRIDY + slab_yy;

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
        MyDouble posw[3], xtmp, ytmp, ztmp;
        if(P[i].Type == 0)
          {
            posw[0] = WRAP_X(SphP[i].Center[0]);
            posw[1] = WRAP_Y(SphP[i].Center[1]);
            posw[2] = WRAP_Z(SphP[i].Center[2]);

            pos = posw;
          }
        else
#endif /* #ifdef CELL_CENTER_GRAVITY */
          pos = P[i].Pos;

        int slab_x  = (int)(to_slab_fac * pos[0]);
        int slab_xx = slab_x + 1;

        if(mode >= 2)
          {
            slab_x %= GRIDX;
            slab_xx %= GRIDX;
          }
        else
          {
            if(slab_x >= GRIDX)
              slab_x -= GRIDX;

            if(slab_xx >= GRIDX)
              slab_xx -= GRIDX;
          }

#ifndef FFT_COLUMN_BASED
        int task0 = myplan.slab_to_task[slab_x];
        int task1 = myplan.slab_to_task[slab_xx];

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
        int slab_y  = (int)(to_slab_fac * pos[1]);
        int slab_yy = slab_y + 1;

        if(mode >= 2)
          {
            slab_y %= GRIDY;
            slab_yy %= GRIDY;
          }
        else
          {
            if(slab_y >= GRIDY)
              slab_y -= GRIDY;

            if(slab_yy >= GRIDY)
              slab_yy -= GRIDY;
          }

        int column0 = slab_x * GRIDY + slab_y;
        int column1 = slab_x * GRIDY + slab_yy;
        int column2 = slab_xx * GRIDY + slab_y;
        int column3 = slab_xx * GRIDY + slab_yy;

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
    subdivide_evenly(GRIDY, MaxThreads, tid, &first_y, &count_y);
    int last_y = first_y + count_y - 1;

    for(i = 0; i < nimport; i++)
      {
        int slab_y  = (int)(to_slab_fac * partin[i].Pos[1]);
        int slab_yy = slab_y + 1;
        double dy   = to_slab_fac * partin[i].Pos[1] - slab_y;

        if(mode >= 2)
          {
            slab_y %= GRIDY;
            slab_yy %= GRIDY;
          }
        else
          {
            if(slab_y >= GRIDY)
              slab_y -= GRIDY;

            if(slab_yy >= GRIDY)
              slab_yy -= GRIDY;
          }

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

            int slab_x  = (int)(to_slab_fac * partin[i].Pos[0]);
            int slab_z  = (int)(to_slab_fac * partin[i].Pos[2]);
            int slab_xx = slab_x + 1;
            int slab_zz = slab_z + 1;

            double dx = to_slab_fac * partin[i].Pos[0] - slab_x;
            double dz = to_slab_fac * partin[i].Pos[2] - slab_z;

            if(mode >= 2)
              {
                slab_x %= GRIDX;
                slab_z %= GRIDZ;
                slab_xx %= GRIDX;
                slab_zz %= GRIDZ;
              }
            else
              {
                if(slab_x >= GRIDX)
                  slab_x -= GRIDX;
                if(slab_z >= GRIDZ)
                  slab_z -= GRIDZ;

                if(slab_xx >= GRIDX)
                  slab_xx -= GRIDX;
                if(slab_zz >= GRIDZ)
                  slab_zz -= GRIDZ;
              }

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
      int slab_x = (int)(to_slab_fac * partin[i].Pos[0]);
      int slab_xx = slab_x + 1;

      int slab_y = (int)(to_slab_fac * partin[i].Pos[1]);
      int slab_yy = slab_y + 1;

      aux[i].dx = to_slab_fac * partin[i].Pos[0] - slab_x;
      aux[i].dy = to_slab_fac * partin[i].Pos[1] - slab_y;

      if(mode >= 2)
        {
          slab_x %= GRIDX;
          slab_xx %= GRIDX;
          slab_y %= GRIDY;
          slab_yy %= GRIDY;
        }
      else
        {
          if(slab_x >= GRIDX)
            slab_x -= GRIDX;
          if(slab_xx >= GRIDX)
            slab_xx -= GRIDX;

          if(slab_y >= GRIDY)
            slab_y -= GRIDY;
          if(slab_yy >= GRIDY)
            slab_yy -= GRIDY;
        }

      aux[i].col0 = slab_x * GRIDY + slab_y;
      aux[i].col1 = slab_x * GRIDY + slab_yy;
      aux[i].col2 = slab_xx * GRIDY + slab_y;
      aux[i].col3 = slab_xx * GRIDY + slab_yy;
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

            int slab_z = (int)(to_slab_fac * partin[i].Pos[2]);
            int slab_zz = slab_z + 1;

            double dz = to_slab_fac * partin[i].Pos[2] - slab_z;

            if(mode >= 2)
              {
                slab_z %= GRIDZ;
                slab_zz %= GRIDZ;
              }
            else
              {
                if(slab_z >= GRIDZ)
                  slab_z -= GRIDZ;

                if(slab_zz >= GRIDZ)
                  slab_zz -= GRIDZ;
              }

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

/* \brief Function to read out the force component corresponding to spatial
 *        dimension 'dim'.
 *
 *  \param[in] dim Dimension to be read out; If dim is negative, potential values
 *             are read out and assigned to  particles.
 *
 *  \return void
 */
static void pmforce_uniform_optimized_readout_forces_or_potential(int dim)
{
#ifdef EVALPOTENTIAL
  double fac = 4 * M_PI * All.G / (pow(All.BoxSize, 3) * STRETCHX * STRETCHY * STRETCHZ); /* to get potential  */
#endif /* #ifdef EVALPOTENTIAL */

  double to_slab_fac = PMGRID / All.BoxSize;

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

      int slab_x = (int)(to_slab_fac * partin[i].Pos[0]);
      int slab_y = (int)(to_slab_fac * partin[i].Pos[1]);
      int slab_z = (int)(to_slab_fac * partin[i].Pos[2]);

      double dx = to_slab_fac * partin[i].Pos[0] - slab_x;
      double dy = to_slab_fac * partin[i].Pos[1] - slab_y;
      double dz = to_slab_fac * partin[i].Pos[2] - slab_z;

      if(slab_x >= GRIDX)
        slab_x -= GRIDX;
      if(slab_y >= GRIDY)
        slab_y -= GRIDY;
      if(slab_z >= GRIDZ)
        slab_z -= GRIDZ;

      int slab_xx = slab_x + 1;
      int slab_yy = slab_y + 1;
      int slab_zz = slab_z + 1;

      if(slab_xx >= GRIDX)
        slab_xx -= GRIDX;
      if(slab_yy >= GRIDY)
        slab_yy -= GRIDY;
      if(slab_zz >= GRIDZ)
        slab_zz -= GRIDZ;

#ifndef FFT_COLUMN_BASED
      if(myplan.slab_to_task[slab_x] == ThisTask)
        {
          slab_x -= myplan.first_slab_x_of_task[ThisTask];

          flistin[i] += grid[FI(slab_x, slab_y, slab_z)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FI(slab_x, slab_y, slab_zz)] * (1.0 - dx) * (1.0 - dy) * (dz) +
                        grid[FI(slab_x, slab_yy, slab_z)] * (1.0 - dx) * (dy) * (1.0 - dz) +
                        grid[FI(slab_x, slab_yy, slab_zz)] * (1.0 - dx) * (dy) * (dz);
        }

      if(myplan.slab_to_task[slab_xx] == ThisTask)
        {
          slab_xx -= myplan.first_slab_x_of_task[ThisTask];

          flistin[i] += grid[FI(slab_xx, slab_y, slab_z)] * (dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FI(slab_xx, slab_y, slab_zz)] * (dx) * (1.0 - dy) * (dz) +
                        grid[FI(slab_xx, slab_yy, slab_z)] * (dx) * (dy) * (1.0 - dz) +
                        grid[FI(slab_xx, slab_yy, slab_zz)] * (dx) * (dy) * (dz);
        }
#else  /* #ifndef FFT_COLUMN_BASED */
      int column0 = slab_x * GRIDY + slab_y;
      int column1 = slab_x * GRIDY + slab_yy;
      int column2 = slab_xx * GRIDY + slab_y;
      int column3 = slab_xx * GRIDY + slab_yy;

      if(column0 >= myplan.base_firstcol && column0 <= myplan.base_lastcol)
        {
          flistin[i] += grid[FC(column0, slab_z)] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) +
                        grid[FC(column0, slab_zz)] * (1.0 - dx) * (1.0 - dy) * (dz);
        }
      if(column1 >= myplan.base_firstcol && column1 <= myplan.base_lastcol)
        {
          flistin[i] +=
              grid[FC(column1, slab_z)] * (1.0 - dx) * (dy) * (1.0 - dz) + grid[FC(column1, slab_zz)] * (1.0 - dx) * (dy) * (dz);
        }

      if(column2 >= myplan.base_firstcol && column2 <= myplan.base_lastcol)
        {
          flistin[i] +=
              grid[FC(column2, slab_z)] * (dx) * (1.0 - dy) * (1.0 - dz) + grid[FC(column2, slab_zz)] * (dx) * (1.0 - dy) * (dz);
        }

      if(column3 >= myplan.base_firstcol && column3 <= myplan.base_lastcol)
        {
          flistin[i] += grid[FC(column3, slab_z)] * (dx) * (dy) * (1.0 - dz) + grid[FC(column3, slab_zz)] * (dx) * (dy) * (dz);
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
        MyDouble posw[3], xtmp, ytmp, ztmp;
        if(P[i].Type == 0)
          {
            posw[0] = WRAP_X(SphP[i].Center[0]);
            posw[1] = WRAP_Y(SphP[i].Center[1]);
            posw[2] = WRAP_Z(SphP[i].Center[2]);

            pos = posw;
          }
        else
#endif /* #ifdef CELL_CENTER_GRAVITY */
          pos = P[i].Pos;

        int slab_x  = (int)(to_slab_fac * pos[0]);
        int slab_xx = slab_x + 1;

        if(slab_x >= GRIDX)
          slab_x -= GRIDX;

        if(slab_xx >= GRIDX)
          slab_xx -= GRIDX;

#ifndef FFT_COLUMN_BASED
        int task0 = myplan.slab_to_task[slab_x];
        int task1 = myplan.slab_to_task[slab_xx];

        double value = flistout[send_offset[task0] + send_count[task0]++];

        if(task0 != task1)
          value += flistout[send_offset[task1] + send_count[task1]++];
#else  /* #ifndef FFT_COLUMN_BASED */
        int slab_y = (int)(to_slab_fac * pos[1]);
        int slab_yy = slab_y + 1;

        if(slab_y >= GRIDY)
          slab_y -= GRIDY;

        if(slab_yy >= GRIDY)
          slab_yy -= GRIDY;

        int column0 = slab_x * GRIDY + slab_y;
        int column1 = slab_x * GRIDY + slab_yy;
        int column2 = slab_xx * GRIDY + slab_y;
        int column3 = slab_xx * GRIDY + slab_yy;

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

/*! \brief Calculates the long-range periodic force given the particle
 *         positions using the PM method.
 *
 *  The force is Gaussian filtered with Asmth, given in
 *  mesh-cell units. We carry out a CIC charge assignment, and compute the
 *  potential by fast Fourier transform methods. The potential is
 *  finite-differenced using a 4-point finite differencing formula, and the
 *  forces are interpolated tri-linearly to the particle positions. The CIC
 *  kernel is deconvolved.
 *
 *  \param[in] mode For mode=0, normal force calculation, mode=1, only density
 *             field construction for a power spectrum calculation. In the
 *             later case, typelist flags the particle types that should be
 *             included in the density field.
 *  \param[in] typelist Flags of particle types included in power spectrum
 *             calculation.
 *
 *  \return void
 */
void pmforce_periodic(int mode, int *typelist)
{
  int x, y, z, xx, yy, zz;

  double tstart = second();

  if(mode == 0)
    mpi_printf("PM-PERIODIC: Starting periodic PM calculation.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

#ifndef NUMPART_PER_TASK_LARGE
  if((((long long)NumPart) << 3) >= (((long long)1) << 31))
    terminate("We are dealing with a too large particle number per MPI rank - enabling NUMPART_PER_TASK_LARGE might help.");
#endif /* #ifndef NUMPART_PER_TASK_LARGE */

  double asmth2 = All.Asmth[0] * All.Asmth[0];
  double d      = All.BoxSize / PMGRID;
  double dhalf  = 0.5 * d;

  double fac = 4 * M_PI * All.G / (pow(All.BoxSize, 3) * STRETCHX * STRETCHY * STRETCHZ); /* to get potential  */

  fac *= 1 / (2 * d); /* for finite differencing */

#ifdef PM_ZOOM_OPTIMIZED
  pmforce_zoom_optimized_prepare_density(mode, typelist);
#else  /* #ifdef PM_ZOOM_OPTIMIZED */
  pmforce_uniform_optimized_prepare_density(mode);
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

  if(mode != 0)
    {
      /* used to measure powerspectrum */
    }
  else
    {
      /* multiply with Green's function in order to obtain the potential (or forces for spectral diffencing) */

      double kfacx = 2.0 * M_PI / (GRIDX * d);
      double kfacy = 2.0 * M_PI / (GRIDY * d);
      double kfacz = 2.0 * M_PI / (GRIDZ * d);

#ifdef FFT_COLUMN_BASED
      for(large_array_offset ip = 0; ip < myplan.second_transposed_ncells; ip++)
        {
          large_array_offset ipcell = ip + ((large_array_offset)myplan.second_transposed_firstcol) * GRIDX;
          y                         = ipcell / (GRIDX * GRIDz);
          int yr                    = ipcell % (GRIDX * GRIDz);
          z                         = yr / GRIDX;
          x                         = yr % GRIDX;
#else  /* #ifdef FFT_COLUMN_BASED */
      for(x = 0; x < GRIDX; x++)
        for(y = myplan.slabstart_y; y < myplan.slabstart_y + myplan.nslab_y; y++)
          for(z = 0; z < GRIDz; z++)
            {
#endif /* #ifdef FFT_COLUMN_BASED #else */
          if(x >= (GRIDX / 2))
            xx = x - GRIDX;
          else
            xx = x;
          if(y >= (GRIDY / 2))
            yy = y - GRIDY;
          else
            yy = y;
          if(z >= (GRIDZ / 2))
            zz = z - GRIDZ;
          else
            zz = z;

          double kx = kfacx * xx;
          double ky = kfacy * yy;
          double kz = kfacz * zz;

          double k2 = kx * kx + ky * ky + kz * kz;

          if(k2 > 0)
            {
              double smth = -exp(-k2 * asmth2) / k2;

              /* do deconvolution */

              double fx = 1, fy = 1, fz = 1;

              if(xx != 0)
                {
                  fx = kx * dhalf;
                  fx = sin(fx) / fx;
                }
              if(yy != 0)
                {
                  fy = ky * dhalf;
                  fy = sin(fy) / fy;
                }
              if(zz != 0)
                {
                  fz = kz * dhalf;
                  fz = sin(fz) / fz;
                }

              double ff     = 1 / (fx * fy * fz);
              double deconv = ff * ff * ff * ff;

              smth *= deconv; /* deconvolution */

#ifndef FFT_COLUMN_BASED
              large_array_offset ip = ((large_array_offset)GRIDz) * (GRIDX * (y - myplan.slabstart_y) + x) + z;
#endif /* #ifndef FFT_COLUMN_BASED */

              fft_of_rhogrid[ip][0] *= smth;
              fft_of_rhogrid[ip][1] *= smth;
            }
        }

#ifdef FFT_COLUMN_BASED
      if(myplan.second_transposed_firstcol == 0)
        fft_of_rhogrid[0][0] = fft_of_rhogrid[0][1] = 0.0;
#else  /* #ifdef FFT_COLUMN_BASED */
      if(myplan.slabstart_y == 0)
        fft_of_rhogrid[0][0] = fft_of_rhogrid[0][1] = 0.0;
#endif /* #ifdef FFT_COLUMN_BASED #else */

        /* Do the inverse FFT to get the potential/forces */

#ifndef FFT_COLUMN_BASED
      my_slab_based_fft(&myplan, &rhogrid[0], &workspace[0], -1);
#else  /* #ifndef FFT_COLUMN_BASED */
      my_column_based_fft(&myplan, workspace, rhogrid, -1);
#endif /* #ifndef FFT_COLUMN_BASED #else */

      /* Now rhogrid holds the potential/forces */

#ifdef EVALPOTENTIAL
#ifdef PM_ZOOM_OPTIMIZED
      pmforce_zoom_optimized_readout_forces_or_potential(-1);
#else  /* #ifdef PM_ZOOM_OPTIMIZED */
      pmforce_uniform_optimized_readout_forces_or_potential(-1);
#endif /* #ifdef PM_ZOOM_OPTIMIZED #else */
#endif /* #ifdef EVALPOTENTIAL */

      /* get the force components by finite differencing of the potential for each dimension,
       * and send the results back to the right CPUs
       */
      for(int dim = 2; dim >= 0; dim--) /* Calculate each component of the force. */
        {
          /* we do the x component last, because for differencing the potential in the x-direction, we need to construct the transpose
           */

#ifndef FFT_COLUMN_BASED
          if(dim == 0)
            {
              my_slab_transposeA(&myplan, rhogrid,
                                 forcegrid); /* compute the transpose of the potential field for finite differencing */
              /* note: for the x-direction, we difference the transposed field */

              for(x = 0; x < GRIDX; x++)
                for(y = 0; y < myplan.nslab_y; y++)
                  for(z = 0; z < GRIDZ; z++)
                    {
                      int xrr = x + 2, xll = x - 2, xr = x + 1, xl = x - 1;
                      if(xr >= GRIDX)
                        xr -= GRIDX;
                      if(xrr >= GRIDX)
                        xrr -= GRIDX;
                      if(xl < 0)
                        xl += GRIDX;
                      if(xll < 0)
                        xll += GRIDX;

                      forcegrid[NI(x, y, z)] = fac * ((4.0 / 3) * (rhogrid[NI(xl, y, z)] - rhogrid[NI(xr, y, z)]) -
                                                      (1.0 / 6) * (rhogrid[NI(xll, y, z)] - rhogrid[NI(xrr, y, z)]));
                    }

              my_slab_transposeB(&myplan, forcegrid, rhogrid); /* reverse the transpose from above */
            }
          else
            {
              for(y = 0; y < GRIDY; y++)
                for(x = 0; x < myplan.nslab_x; x++)
                  for(z = 0; z < GRIDZ; z++)
                    {
                      if(dim == 1)
                        {
                          int yr = y + 1, yl = y - 1, yrr = y + 2, yll = y - 2;
                          if(yr >= GRIDY)
                            yr -= GRIDY;
                          if(yrr >= GRIDY)
                            yrr -= GRIDY;
                          if(yl < 0)
                            yl += GRIDY;
                          if(yll < 0)
                            yll += GRIDY;

                          forcegrid[FI(x, y, z)] = fac * ((4.0 / 3) * (rhogrid[FI(x, yl, z)] - rhogrid[FI(x, yr, z)]) -
                                                          (1.0 / 6) * (rhogrid[FI(x, yll, z)] - rhogrid[FI(x, yrr, z)]));
                        }
                      else if(dim == 2)
                        {
                          int zr = z + 1, zl = z - 1, zrr = z + 2, zll = z - 2;
                          if(zr >= GRIDZ)
                            zr -= GRIDZ;
                          if(zrr >= GRIDZ)
                            zrr -= GRIDZ;
                          if(zl < 0)
                            zl += GRIDZ;
                          if(zll < 0)
                            zll += GRIDZ;

                          forcegrid[FI(x, y, z)] = fac * ((4.0 / 3) * (rhogrid[FI(x, y, zl)] - rhogrid[FI(x, y, zr)]) -
                                                          (1.0 / 6) * (rhogrid[FI(x, y, zll)] - rhogrid[FI(x, y, zrr)]));
                        }
                    }
            }

#else  /* #ifndef FFT_COLUMN_BASED */

          if(dim == 2)
            {
              for(large_array_offset i = 0; i < myplan.base_ncol; i++)
                {
                  fft_real *forcep = &forcegrid[GRID2 * i];
                  fft_real *potp   = &rhogrid[GRID2 * i];

                  for(int z = 0; z < GRIDZ; z++)
                    {
                      int zr  = z + 1;
                      int zl  = z - 1;
                      int zrr = z + 2;
                      int zll = z - 2;

                      if(zr >= GRIDZ)
                        zr -= GRIDZ;
                      if(zrr >= GRIDZ)
                        zrr -= GRIDZ;
                      if(zl < 0)
                        zl += GRIDZ;
                      if(zll < 0)
                        zll += GRIDZ;

                      forcep[z] = fac * ((4.0 / 3) * (potp[zl] - potp[zr]) - (1.0 / 6) * (potp[zll] - potp[zrr]));
                    }
                }
            }
          else if(dim == 1)
            {
              fft_real *scratch = mymalloc("scratch", myplan.fftsize * sizeof(fft_real)); /* need a third field as scratch space */
              memcpy(scratch, rhogrid, myplan.fftsize * sizeof(fft_real));

              my_fft_swap23(&myplan, scratch, forcegrid);

              for(large_array_offset i = 0; i < myplan.ncol_XZ; i++)
                {
                  fft_real *forcep = &scratch[GRIDY * i];
                  fft_real *potp   = &forcegrid[GRIDY * i];

                  for(int y = 0; y < GRIDY; y++)
                    {
                      int yr  = y + 1;
                      int yl  = y - 1;
                      int yrr = y + 2;
                      int yll = y - 2;

                      if(yr >= GRIDY)
                        yr -= GRIDY;
                      if(yrr >= GRIDY)
                        yrr -= GRIDY;
                      if(yl < 0)
                        yl += GRIDY;
                      if(yll < 0)
                        yll += GRIDY;

                      forcep[y] = fac * ((4.0 / 3) * (potp[yl] - potp[yr]) - (1.0 / 6) * (potp[yll] - potp[yrr]));
                    }
                }

              my_fft_swap23back(&myplan, scratch, forcegrid);
              myfree(scratch);
            }
          else if(dim == 0)
            {
              fft_real *scratch = mymalloc("scratch", myplan.fftsize * sizeof(fft_real)); /* need a third field as scratch space */
              memcpy(scratch, rhogrid, myplan.fftsize * sizeof(fft_real));

              my_fft_swap13(&myplan, scratch, forcegrid);

              for(large_array_offset i = 0; i < myplan.ncol_YZ; i++)
                {
                  fft_real *forcep = &scratch[GRIDX * i];
                  fft_real *potp   = &forcegrid[GRIDX * i];

                  for(int x = 0; x < GRIDX; x++)
                    {
                      int xr  = x + 1;
                      int xl  = x - 1;
                      int xrr = x + 2;
                      int xll = x - 2;

                      if(xr >= GRIDX)
                        xr -= GRIDX;
                      if(xrr >= GRIDX)
                        xrr -= GRIDX;
                      if(xl < 0)
                        xl += GRIDX;
                      if(xll < 0)
                        xll += GRIDX;

                      forcep[x] = fac * ((4.0 / 3) * (potp[xl] - potp[xr]) - (1.0 / 6) * (potp[xll] - potp[xrr]));
                    }
                }

              my_fft_swap13back(&myplan, scratch, forcegrid);
              myfree(scratch);
            }
#endif /* #ifndef FFT_COLUMN_BASED #else */

#ifdef PM_ZOOM_OPTIMIZED
          pmforce_zoom_optimized_readout_forces_or_potential(dim);
#else  /* #ifdef PM_ZOOM_OPTIMIZED */
          pmforce_uniform_optimized_readout_forces_or_potential(dim);
#endif /* #ifdef PM_ZOOM_OPTIMIZED #else */
        }
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
#endif /* #ifdef PM_ZOOM_OPTIMIZED */

  double tend = second();

  if(mode == 0)
    mpi_printf("PM-PERIODIC: done.  (took %g seconds)\n", timediff(tstart, tend));
}

#ifdef PM_ZOOM_OPTIMIZED

/*! \brief Sort function for 'part' array indices.
 *
 * Sorts the indices into the 'part' array by the global index of the
 * corresponding 'part_slab_data' struct.
 *
 * \param[in] a Index to be compared.
 * \param[in] b Index to be compared.
 *
 * \return sort result
 */
static int pm_periodic_compare_sortindex(const void *a, const void *b)
{
  if(part[*(int *)a].globalindex < part[*(int *)b].globalindex)
    return -1;

  if(part[*(int *)a].globalindex > part[*(int *)b].globalindex)
    return +1;

  return 0;
}

/*! \brief Implements the sorting function for mysort_pmperiodic().
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

/*! \brief Sort the index array b of n entries using the sort kernel
 *         cmp.
 *
 *  The parameter s is set to sizeof(int). The index array b is sorted
 *  according to the globalindex field of the referenced item in the 'part'
 *  array.
 *
 *  \param[in, out] b The index array to sort.
 *  \param[in] n Number of entries in array b.
 *  \param[in] s Size of each entry (must be sizeof(int)).
 *  \param[in] cmp Comparison function.
 */
static void mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *))
{
  const size_t size = n * s;

  large_numpart_type *tmp = (large_numpart_type *)mymalloc("tmp", size);

  msort_pmperiodic_with_tmp((large_numpart_type *)b, n, tmp);

  myfree(tmp);
}
#endif /* #ifdef PM_ZOOM_OPTIMIZED */

#endif /* #if defined(PMGRID) */
