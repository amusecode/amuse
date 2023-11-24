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
 * \file        src/gravity/pm/pm_periodic2d.c
 * \date        05/2018
 * \brief       Routines for periodic PM-force computation in 2d.
 * \details     contains functions:
 *                void pm2d_init_periodic(void)
 *                void pm2d_init_periodic_allocate(void)
 *                void pm2d_init_periodic_free(void)
 *                void pm2d_force_periodic(int mode)
 *                int pm2d_periodic_compare_sortindex(const void *a, const
 *                  void *b)
 *                static void pm2d_msort_pmperiodic_with_tmp(int *b, size_t n,
 *                  int *t)
 *                void pm2d_mysort_pmperiodic(void *b, size_t n, size_t s,
 *                  int (*cmp) (const void *, const void *))
 *                void pm2d_periodic_transposeA(fftw_real * field,
 *                  fftw_real * scratch)
 *                void pm2d_periodic_transposeB(fftw_real * field,
 *                  fftw_real * scratch)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 21.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef PMGRID
#ifndef GRAVITY_NOT_PERIODIC
#ifdef TWODIMS

#ifdef NOTYPEPREFIX_FFTW
#include <rfftw_mpi.h>
#else /* #ifdef NOTYPEPREFIX_FFTW */
#ifdef DOUBLEPRECISION_FFTW
#include <drfftw_mpi.h> /* double precision FFTW */
#else                   /* #ifdef DOUBLEPRECISION_FFTW */
#include <srfftw_mpi.h>
#endif /* #ifdef DOUBLEPRECISION_FFTW #else */
#endif /* #ifdef NOTYPEPREFIX_FFTW #else */

#include "../../main/allvars.h"
#include "../../main/proto.h"

#define PMGRID2 (2 * (PMGRID / 2 + 1))

#if(PMGRID > 1024)
typedef long long large_array_offset;
#else  /* #if (PMGRID > 1024) */
typedef unsigned int large_array_offset;
#endif /* #if (PMGRID > 1024) #else */

#define d_fftw_real fftw_real

static rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;

static int slab_to_task[PMGRID];
static int *slabs_x_per_task;
static int *first_slab_x_of_task;

static int slabstart_x, nslab_x, slabstart_y, nslab_y, smallest_slab;

static int fftsize, maxfftsize;

static fftw_real *rhogrid, *forcegrid, *workspace;
static d_fftw_real *d_rhogrid, *d_forcegrid, *d_workspace;

static fftw_complex *fft_of_rhogrid;

static MyFloat to_slab_fac;

void pm2d_periodic_transposeA(fftw_real *field, fftw_real *scratch);
void pm2d_periodic_transposeB(fftw_real *field, fftw_real *scratch);
int pm2d_periodic_compare_sortindex(const void *a, const void *b);

/*! \brief Data for fft slab.
 */
static struct part_slab_data
{
  large_array_offset globalindex;
  int partindex;
  int localindex;
} * part;

static int *part_sortindex;

/*! \brief This routines generates the FFTW-plans to carry out the parallel
 *         FFTs later on. Some auxiliary variables are also initialized.
 *
 *  \return void
 */
void pm2d_init_periodic(void)
{
  int i;
  int slab_to_task_local[PMGRID];

  All.Asmth[0] = ASMTH * All.BoxSize / PMGRID;
  All.Rcut[0]  = RCUT * All.Asmth[0];

  /* Set up the FFTW plan files. */

  fft_forward_plan = rfftw2d_mpi_create_plan(MPI_COMM_WORLD, PMGRID, PMGRID, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fft_inverse_plan = rfftw2d_mpi_create_plan(MPI_COMM_WORLD, PMGRID, PMGRID, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

  /* Workspace out the ranges on each processor. */

  rfftwnd_mpi_local_sizes(fft_forward_plan, &nslab_x, &slabstart_x, &nslab_y, &slabstart_y, &fftsize);

  for(i = 0; i < PMGRID; i++)
    slab_to_task_local[i] = 0;

  for(i = 0; i < nslab_x; i++)
    slab_to_task_local[slabstart_x + i] = ThisTask;

  MPI_Allreduce(slab_to_task_local, slab_to_task, PMGRID, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(&nslab_x, &smallest_slab, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  slabs_x_per_task = (int *)mymalloc("slabs_per_task", NTask * sizeof(int));
  MPI_Allgather(&nslab_x, 1, MPI_INT, slabs_x_per_task, 1, MPI_INT, MPI_COMM_WORLD);

  first_slab_x_of_task = (int *)mymalloc("first_slab_of_task", NTask * sizeof(int));
  MPI_Allgather(&slabstart_x, 1, MPI_INT, first_slab_x_of_task, 1, MPI_INT, MPI_COMM_WORLD);

  to_slab_fac = PMGRID / All.BoxSize;

  MPI_Allreduce(&fftsize, &maxfftsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  printf("maxfftsize=%d PMGRID=%d\n", maxfftsize, PMGRID);
}

/*! \brief Allocates memory for 2d PM algorithm.
 *
 *  This function allocates the memory neeed to compute the long-range PM
 *  force. Three fields are used, one to hold the density (and its FFT, and
 *  then the real-space potential), one to hold the force field obtained by
 *  finite differencing, and finally a workspace field, which is used both as
 *  workspace for the parallel FFT, and as buffer for the communication
 *  algorithm used in the force computation.
 *
 *  \return void
 */
void pm2d_init_periodic_allocate(void)
{
  double bytes_tot = 0;
  size_t bytes;

  /* allocate the memory to hold the FFT fields */

  rhogrid = (fftw_real *)mymalloc("rhogrid", bytes = maxfftsize * sizeof(d_fftw_real));
  bytes_tot += bytes;

  forcegrid = (fftw_real *)mymalloc("forcegrid", bytes = maxfftsize * sizeof(d_fftw_real));
  bytes_tot += bytes;

  part = (struct part_slab_data *)mymalloc("part", bytes = 4 * NumPart * sizeof(struct part_slab_data));
  bytes_tot += bytes;

  part_sortindex = (int *)mymalloc("part_sortindex", bytes = 4 * NumPart * sizeof(int));
  bytes_tot += bytes;

  if(ThisTask == 0)
    printf("Using %g MByte for periodic FFT computation. (presently allocated=%g MB)\n", bytes_tot / (1024.0 * 1024.0),
           AllocatedBytes / (1024.0 * 1024.0));

  workspace = forcegrid;

  fft_of_rhogrid = (fftw_complex *)&rhogrid[0];

  d_rhogrid   = (d_fftw_real *)rhogrid;
  d_forcegrid = (d_fftw_real *)forcegrid;
  d_workspace = (d_fftw_real *)workspace;
}

/*! \brief This routine frees the space allocated for the parallel FFT
 *         algorithm.
 *
 *  \return void
 */
void pm2d_init_periodic_free(void)
{
  /* allocate the memory to hold the FFT fields */
  myfree(part_sortindex);
  myfree(part);
  myfree(forcegrid);
  myfree(rhogrid);
}

/*! \brief Long range periodic 2d gravity.
 *
 *  Calculates the long-range periodic force given the particle positions
 *  using the PM method. The force is Gaussian filtered with Asmth, given in
 *  mesh-cell units. We carry out a CIC charge assignment, and compute the
 *  potenial by Fourier transform methods. The potential is finite differenced
 *  using a 4-point finite differencing formula, and the forces are
 *  interpolated tri-linearly to the particle positions. The CIC kernel is
 *  deconvolved. Note that the particle distribution is not in the slab
 *  decomposition that is used for the FFT. Instead, overlapping patches
 *  between local domains and FFT slabs are communicated as needed.
 *
 *  \param[in] mode 0: normal PM force; 1: calculate mesh correction vector.
 *
 *  \return void
 */
void pm2d_force_periodic(int mode)
{
  double k2, kx, ky, smth;
  double dx, dy, weight;
  double fx, fy, ff;
  double asmth2, fac, acc_dim;
  int i, j, N, slab, level, sendTask, recvTask, task;
  int x, y, yl, yr, yll, yrr, ip, dim;
  int slab_x, slab_y;
  int slab_xx, slab_yy;
  int num_on_grid, num_field_points, pindex, xx, yy;
  MPI_Status status;
  int *localfield_count, *localfield_first, *localfield_offset, *localfield_togo;
  large_array_offset offset, *localfield_globalindex, *import_globalindex;
  d_fftw_real *localfield_d_data, *import_d_data;
  fftw_real *localfield_data, *import_data;

  if(ThisTask == 0)
    {
      printf("Starting periodic PM-2d calculation.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));
      myflush(stdout);
    }

  asmth2 = (2 * M_PI) * All.Asmth[0] / All.BoxSize;
  asmth2 *= asmth2;

  fac = All.G / (M_PI * All.BoxSize);    /* to get potential */
  fac *= 1 / (2 * All.BoxSize / PMGRID); /* for finite differencing */

  if(mode == 1)
    {
      fac *= 1.0 / (All.G) * All.BoxSize;
    }
  else
    {
      fac *= All.BoxSize;
    }

  pm2d_init_periodic_allocate();

  if(mode == 0)
    N = NumPart;
  else
    N = NumGas;

  /* determine the cells each particles accesses */
  for(i = 0, num_on_grid = 0; i < N; i++)
    {
      slab_x = (int)(to_slab_fac * P[i].Pos[0]);
      slab_y = (int)(to_slab_fac * P[i].Pos[1]);

      if(slab_x >= PMGRID)
        slab_x = PMGRID - 1;
      if(slab_y >= PMGRID)
        slab_y = PMGRID - 1;

      for(xx = 0; xx < 2; xx++)
        for(yy = 0; yy < 2; yy++)
          {
            slab_xx = slab_x + xx;
            slab_yy = slab_y + yy;

            if(slab_xx >= PMGRID)
              slab_xx -= PMGRID;
            if(slab_yy >= PMGRID)
              slab_yy -= PMGRID;

            offset = (PMGRID2 * slab_xx + slab_yy);

            part[num_on_grid].partindex   = (i << 2) + (xx << 1) + yy;
            part[num_on_grid].globalindex = offset;
            part_sortindex[num_on_grid]   = num_on_grid;
            num_on_grid++;
          }
    }

  /* note: num_on_grid will be  4 times larger than the particle number,
     but num_field_points will generally be much smaller */

  /* bring the part-field into the order of the accessed cells. This allow the removal of duplicates */
  pm2d_mysort_pmperiodic(part_sortindex, num_on_grid, sizeof(int), pm2d_periodic_compare_sortindex);

  /* determine the number of unique field points */
  for(i = 0, num_field_points = 0; i < num_on_grid; i++)
    {
      if(i > 0)
        if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
          continue;

      num_field_points++;
    }

  /* allocate the local field */
  localfield_globalindex = (large_array_offset *)mymalloc("first_slab_of_task", num_field_points * sizeof(large_array_offset));
  localfield_d_data      = (d_fftw_real *)mymalloc("localfield_d_data", num_field_points * sizeof(d_fftw_real));
  localfield_data        = (fftw_real *)localfield_d_data;
  localfield_first       = (int *)mymalloc("localfield_d_data", NTask * sizeof(int));
  localfield_count       = (int *)mymalloc("localfield_count", NTask * sizeof(int));
  localfield_offset      = (int *)mymalloc("localfield_count", NTask * sizeof(int));
  localfield_togo        = (int *)mymalloc("localfield_togo", NTask * NTask * sizeof(int));

  for(i = 0; i < NTask; i++)
    {
      localfield_first[i] = 0;
      localfield_count[i] = 0;
    }

  /* establish the cross link between the part[] array and the local list of
     mesh points. Also, count on which CPU how many of the needed field points are stored */
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

      slab = part[part_sortindex[i]].globalindex / PMGRID2;
      task = slab_to_task[slab];
      if(localfield_count[task] == 0)
        localfield_first[task] = num_field_points;
      localfield_count[task]++;
    }
  num_field_points++;

  for(i = 1, localfield_offset[0] = 0; i < NTask; i++)
    localfield_offset[i] = localfield_offset[i - 1] + localfield_count[i - 1];

  /* now bin the local particle data onto the mesh list */

  for(i = 0; i < num_field_points; i++)
    localfield_d_data[i] = 0;

  for(i = 0; i < num_on_grid; i += 4)
    {
      pindex = (part[i].partindex >> 2);

      slab_x = (int)(to_slab_fac * P[pindex].Pos[0]);
      slab_y = (int)(to_slab_fac * P[pindex].Pos[1]);

      dx = to_slab_fac * P[pindex].Pos[0] - slab_x;
      dy = to_slab_fac * P[pindex].Pos[1] - slab_y;

      weight = P[pindex].Mass;

      localfield_d_data[part[i + 0].localindex] += weight * (1.0 - dx) * (1.0 - dy);
      localfield_d_data[part[i + 1].localindex] += weight * (1.0 - dx) * dy;
      localfield_d_data[part[i + 2].localindex] += weight * (dx) * (1.0 - dy);
      localfield_d_data[part[i + 3].localindex] += weight * (dx)*dy;
    }

  /* clear local FFT-mesh density field */
  for(i = 0; i < fftsize; i++)
    d_rhogrid[i] = 0;

  /* exchange data and add contributions to the local mesh-path */

  MPI_Allgather(localfield_count, NTask, MPI_INT, localfield_togo, NTask, MPI_INT, MPI_COMM_WORLD);

  for(level = 0; level < (1 << PTask); level++) /* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
        {
          if(level > 0)
            {
              import_d_data =
                  (d_fftw_real *)mymalloc("import_d_data", localfield_togo[recvTask * NTask + ThisTask] * sizeof(d_fftw_real));
              import_globalindex = (large_array_offset *)mymalloc(
                  "import_d_data", localfield_togo[recvTask * NTask + ThisTask] * sizeof(large_array_offset));

              if(localfield_togo[sendTask * NTask + recvTask] > 0 || localfield_togo[recvTask * NTask + sendTask] > 0)
                {
                  MPI_Sendrecv(localfield_d_data + localfield_offset[recvTask],
                               localfield_togo[sendTask * NTask + recvTask] * sizeof(d_fftw_real), MPI_BYTE, recvTask, TAG_NONPERIOD_A,
                               import_d_data, localfield_togo[recvTask * NTask + sendTask] * sizeof(d_fftw_real), MPI_BYTE, recvTask,
                               TAG_NONPERIOD_A, MPI_COMM_WORLD, &status);

                  MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
                               localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask,
                               TAG_NONPERIOD_B, import_globalindex,
                               localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset), MPI_BYTE, recvTask,
                               TAG_NONPERIOD_B, MPI_COMM_WORLD, &status);
                }
            }
          else
            {
              import_d_data      = localfield_d_data + localfield_offset[ThisTask];
              import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
            }

          for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
            {
              /* determine offset in local FFT slab */
              offset = import_globalindex[i] - first_slab_x_of_task[ThisTask] * PMGRID2;

              d_rhogrid[offset] += import_d_data[i];
            }

          if(level > 0)
            {
              myfree(import_globalindex);
              myfree(import_d_data);
            }
        }
    }

  /* Do the FFT of the density field */

  rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

  /* multiply with Green's function for the potential */

  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < PMGRID; x++)
      {
        if(x > PMGRID / 2)
          kx = x - PMGRID;
        else
          kx = x;
        if(y > PMGRID / 2)
          ky = y - PMGRID;
        else
          ky = y;

        k2 = kx * kx + ky * ky;

        if(k2 > 0)
          {
            smth = -exp(-k2 * asmth2) / k2;

            /* do deconvolution */

            fx = fy = 1;
            if(kx != 0)
              {
                fx = (M_PI * kx) / PMGRID;
                fx = sin(fx) / fx;
              }
            if(ky != 0)
              {
                fy = (M_PI * ky) / PMGRID;
                fy = sin(fy) / fy;
              }
            ff = 1 / (fx * fy);
            smth *= ff * ff * ff * ff;

            /* end deconvolution */

            ip = PMGRID * (y - slabstart_y) + x;
            fft_of_rhogrid[ip].re *= smth;
            fft_of_rhogrid[ip].im *= smth;
          }
      }

  if(slabstart_y == 0)
    fft_of_rhogrid[0].re = fft_of_rhogrid[0].im = 0.0;

  /* Do the inverse FFT to get the potential */

  rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

#ifdef EVALPOTENTIAL /* now read out the potential */
  if(mode == 0)
    {
      for(level = 0; level < (1 << PTask); level++) /* note: for level=0, target is the same task */
        {
          sendTask = ThisTask;
          recvTask = ThisTask ^ level;

          if(recvTask < NTask)
            {
              if(level > 0)
                {
                  import_data = (fftw_real *)mymalloc("import_data", localfield_togo[recvTask * NTask + ThisTask] * sizeof(fftw_real));
                  import_globalindex = (large_array_offset *)mymalloc(
                      "import_data", localfield_togo[recvTask * NTask + ThisTask] * sizeof(large_array_offset));

                  if(localfield_togo[sendTask * NTask + recvTask] > 0 || localfield_togo[recvTask * NTask + sendTask] > 0)
                    {
                      MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
                                   localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask,
                                   TAG_NONPERIOD_C, import_globalindex,
                                   localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset), MPI_BYTE, recvTask,
                                   TAG_NONPERIOD_C, MPI_COMM_WORLD, &status);
                    }
                }
              else
                {
                  import_data        = localfield_data + localfield_offset[ThisTask];
                  import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
                }

              for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
                {
                  offset         = import_globalindex[i] - first_slab_x_of_task[ThisTask] * ((large_array_offset)PMGRID2);
                  import_data[i] = rhogrid[offset];
                }

              if(level > 0)
                {
                  MPI_Sendrecv(import_data, localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real), MPI_BYTE, recvTask,
                               TAG_NONPERIOD_A, localfield_data + localfield_offset[recvTask],
                               localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real), MPI_BYTE, recvTask, TAG_NONPERIOD_A,
                               MPI_COMM_WORLD, &status);

                  myfree(import_globalindex);
                  myfree(import_data);
                }
            }
        }

      /* read out the potential values, which all have been assembled in localfield_data */

      double pot;

      for(i = 0, j = 0; i < N; i++)
        {
          while(j < num_on_grid && (part[j].partindex >> 2) != i)
            j++;

          slab_x = (int)(to_slab_fac * P[i].Pos[0]);
          dx     = to_slab_fac * P[i].Pos[0] - slab_x;

          slab_y = (int)(to_slab_fac * P[i].Pos[1]);
          dy     = to_slab_fac * P[i].Pos[1] - slab_y;

          pot = +localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) +
                localfield_data[part[j + 1].localindex] * (1.0 - dx) * dy + localfield_data[part[j + 2].localindex] * dx * (1.0 - dy) +
                localfield_data[part[j + 3].localindex] * dx * dy;

          P[i].PM_Potential += pot * fac * (2 * All.BoxSize / PMGRID);
          /* compensate the finite differencing factor */;
        }
    }
#endif /* #ifdef EVALPOTENTIAL */

  /* get the force components by finite differencing the potential for each dimension,
     and send back the results to the right CPUs */

  for(dim = 1; dim >= 0; dim--) /* Calculate each component of the force. */
    { /* we do the x component last, because for differencing the potential in the x-direction, we need to contruct the transpose */
      if(dim == 0)
        pm2d_periodic_transposeA(rhogrid, forcegrid); /* compute the transpose of the potential field */

      for(xx = slabstart_x; xx < (slabstart_x + nslab_x); xx++)
        for(y = 0; y < PMGRID; y++)
          {
            x = xx - slabstart_x;

            yrr = yll = yr = yl = y;

            yr  = y + 1;
            yl  = y - 1;
            yrr = y + 2;
            yll = y - 2;
            if(yr >= PMGRID)
              yr -= PMGRID;
            if(yrr >= PMGRID)
              yrr -= PMGRID;
            if(yl < 0)
              yl += PMGRID;
            if(yll < 0)
              yll += PMGRID;

            if(dim == 0)
              {
                forcegrid[x + y * nslab_x] = fac * ((4.0 / 3) * (rhogrid[(x + yl * nslab_x)] - rhogrid[(x + yr * nslab_x)]) -
                                                    (1.0 / 6) * (rhogrid[(x + yll * nslab_x)] - rhogrid[(x + yrr * nslab_x)]));
              }
            else
              {
                forcegrid[PMGRID2 * x + y] = fac * ((4.0 / 3) * (rhogrid[PMGRID2 * x + yl] - rhogrid[PMGRID2 * x + yr]) -
                                                    (1.0 / 6) * (rhogrid[PMGRID2 * x + yll] - rhogrid[PMGRID2 * x + yrr]));
              }
          }

      if(dim == 0)
        pm2d_periodic_transposeB(forcegrid, rhogrid); /* compute the transpose of the potential field */

      /* send the force components to the right processors */

      for(level = 0; level < (1 << PTask); level++) /* note: for level=0, target is the same task */
        {
          sendTask = ThisTask;
          recvTask = ThisTask ^ level;

          if(recvTask < NTask)
            {
              if(level > 0)
                {
                  import_data = (fftw_real *)mymalloc("import_data", localfield_togo[recvTask * NTask + ThisTask] * sizeof(fftw_real));
                  import_globalindex = (large_array_offset *)mymalloc(
                      "import_data", localfield_togo[recvTask * NTask + ThisTask] * sizeof(large_array_offset));

                  if(localfield_togo[sendTask * NTask + recvTask] > 0 || localfield_togo[recvTask * NTask + sendTask] > 0)
                    {
                      MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
                                   localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset), MPI_BYTE, recvTask,
                                   TAG_NONPERIOD_C, import_globalindex,
                                   localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset), MPI_BYTE, recvTask,
                                   TAG_NONPERIOD_C, MPI_COMM_WORLD, &status);
                    }
                }
              else
                {
                  import_data        = localfield_data + localfield_offset[ThisTask];
                  import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
                }

              for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
                {
                  /* determine offset in local FFT slab */
                  offset         = import_globalindex[i] - first_slab_x_of_task[ThisTask] * PMGRID2;
                  import_data[i] = forcegrid[offset];
                }

              if(level > 0)
                {
                  MPI_Sendrecv(import_data, localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real), MPI_BYTE, recvTask,
                               TAG_NONPERIOD_A, localfield_data + localfield_offset[recvTask],
                               localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real), MPI_BYTE, recvTask, TAG_NONPERIOD_A,
                               MPI_COMM_WORLD, &status);

                  myfree(import_globalindex);
                  myfree(import_data);
                }
            }
        }

      /* read out the forces, which all have been assembled in localfield_data */

      for(i = 0, j = 0; i < N; i++)
        {
          while(j < num_on_grid && (part[j].partindex >> 2) != i)
            j++;

          slab_x = (int)(to_slab_fac * P[i].Pos[0]);
          dx     = to_slab_fac * P[i].Pos[0] - slab_x;

          slab_y = (int)(to_slab_fac * P[i].Pos[1]);
          dy     = to_slab_fac * P[i].Pos[1] - slab_y;

          acc_dim = +localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) +
                    localfield_data[part[j + 1].localindex] * (1.0 - dx) * dy +
                    localfield_data[part[j + 2].localindex] * (dx) * (1.0 - dy) + localfield_data[part[j + 3].localindex] * (dx)*dy;

          P[i].GravPM[dim] += acc_dim;
        }
    }

  /* free locallist */
  myfree(localfield_togo);
  myfree(localfield_offset);
  myfree(localfield_count);
  myfree(localfield_first);
  myfree(localfield_d_data);
  myfree(localfield_globalindex);

  pm2d_init_periodic_free();

  mpi_printf("done PM-2d.\n");
}

/*! \brief Compares two objects of type part_slab_data.
 *
 *  According to element globalindex.
 *
 *  \param[in] a Index of first object in part array.
 *  \param[in] b Index of second object in part array.
 *
 *  \return (-1,0,1); -1 if part[a].globalindex < part[b].globalindex
 */
int pm2d_periodic_compare_sortindex(const void *a, const void *b)
{
  if(part[*(int *)a].globalindex < part[*(int *)b].globalindex)
    return -1;

  if(part[*(int *)a].globalindex > part[*(int *)b].globalindex)
    return +1;

  return 0;
}

/*! \brief Merge sort algorithm for 2d periodic particle mesh algorithm.
 *
 *  \param[in, out] b Array to be sorted.
 *  \param[in] n Size of array b.
 *  \param[in, out] t Temporary array.
 *
 *  \return void
 */
static void pm2d_msort_pmperiodic_with_tmp(int *b, size_t n, int *t)
{
  int *tmp;
  int *b1, *b2;
  size_t n1, n2;

  if(n <= 1)
    return;

  n1 = n / 2;
  n2 = n - n1;
  b1 = b;
  b2 = b + n1;

  pm2d_msort_pmperiodic_with_tmp(b1, n1, t);
  pm2d_msort_pmperiodic_with_tmp(b2, n2, t);

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
    memcpy(tmp, b1, n1 * sizeof(int));

  memcpy(b, t, (n - n2) * sizeof(int));
}

/*! \brief Wrapper for sorting algorithm in 2d periodic PM algorithm.
 *
 *  Uses pm2d_msort_pmperiodic_with_tmp.
 *
 *  \param[in, out] b Array to be sorted.
 *  \param[in] n Number of elements in array b.
 *  \param[in] s Size of individual element of b (for memory allocation).
 *  \param[in] cmp Compare function (unused).
 *
 *  \return void
 */
void pm2d_mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *))
{
  const size_t size = n * s;

  int *tmp = (int *)mymalloc("tmp", size);

  pm2d_msort_pmperiodic_with_tmp((int *)b, n, tmp);

  myfree(tmp);
}

/*! \brief Transpose operation for 2d fft.
 *
 *  Used for transposing rhogrid.
 *
 *  \param[in, out] field Field that needs to be transposed.
 *  \param[in, out] scratch Temporary data.
 *
 *  \return void
 */
void pm2d_periodic_transposeA(fftw_real *field, fftw_real *scratch)
{
  int x, y, task;

  for(task = 0; task < NTask; task++)
    for(x = 0; x < nslab_x; x++)
      for(y = first_slab_x_of_task[task]; y < first_slab_x_of_task[task] + slabs_x_per_task[task]; y++)
        {
          scratch[(first_slab_x_of_task[task] * nslab_x + x * slabs_x_per_task[task] + (y - first_slab_x_of_task[task]))] =
              field[PMGRID2 * x + y];
        }

#ifndef NO_ISEND_IRECV_IN_DOMAIN
  MPI_Request *requests;
  int nrequests = 0;

  requests = (MPI_Request *)mymalloc(2 * NTask * sizeof(MPI_Request));

  for(task = 0; task < NTask; task++)
    {
      MPI_Isend(scratch + first_slab_x_of_task[task] * nslab_x, nslab_x * slabs_x_per_task[task] * sizeof(fftw_real), MPI_BYTE, task,
                TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);

      MPI_Irecv(field + first_slab_x_of_task[task] * nslab_x, nslab_x * slabs_x_per_task[task] * sizeof(fftw_real), MPI_BYTE, task,
                TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);
    }

  MPI_Waitall(nrequests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);
#else  /* #ifndef NO_ISEND_IRECV_IN_DOMAIN */
  int ngrp;

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      task = ThisTask ^ ngrp;

      if(task < NTask)
        {
          MPI_Sendrecv(scratch + first_slab_x_of_task[task] * nslab_x, nslab_x * slabs_x_per_task[task] * sizeof(fftw_real), MPI_BYTE,
                       task, TAG_KEY, field + first_slab_x_of_task[task] * nslab_x,
                       nslab_x * slabs_x_per_task[task] * sizeof(fftw_real), MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
        }
    }
#endif /* #ifndef NO_ISEND_IRECV_IN_DOMAIN #else */
}

/*! \brief Transpose operation for 2d fft.
 *
 *  Used for forcegrid transpose.
 *
 *  \param[in, out] field Field that needs to be transposed.
 *  \param[in, out] scratch Temporary data.
 *
 *  \return void
 */
void pm2d_periodic_transposeB(fftw_real *field, fftw_real *scratch)
{
  int x, y, task;

#ifndef NO_ISEND_IRECV_IN_DOMAIN
  MPI_Request *requests;
  int nrequests = 0;

  requests = (MPI_Request *)mymalloc(2 * NTask * sizeof(MPI_Request));

  for(task = 0; task < NTask; task++)
    {
      MPI_Isend(field + first_slab_x_of_task[task] * nslab_x, nslab_x * slabs_x_per_task[task] * sizeof(fftw_real), MPI_BYTE, task,
                TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);

      MPI_Irecv(scratch + first_slab_x_of_task[task] * nslab_x, nslab_x * slabs_x_per_task[task] * sizeof(fftw_real), MPI_BYTE, task,
                TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);
    }

  MPI_Waitall(nrequests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);

#else  /* #ifndef NO_ISEND_IRECV_IN_DOMAIN */
  int ngrp;

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      task = ThisTask ^ ngrp;

      if(task < NTask)
        {
          MPI_Sendrecv(field + first_slab_x_of_task[task] * nslab_x, nslab_x * slabs_x_per_task[task] * sizeof(fftw_real), MPI_BYTE,
                       task, TAG_KEY, scratch + first_slab_x_of_task[task] * nslab_x,
                       nslab_x * slabs_x_per_task[task] * sizeof(fftw_real), MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD,
                       MPI_STATUS_IGNORE);
        }
    }
#endif /* #ifndef NO_ISEND_IRECV_IN_DOMAIN #else */

  for(task = 0; task < NTask; task++)
    for(x = 0; x < nslab_x; x++)
      for(y = first_slab_x_of_task[task]; y < first_slab_x_of_task[task] + slabs_x_per_task[task]; y++)
        {
          field[PMGRID2 * x + y] =
              scratch[(first_slab_x_of_task[task] * nslab_x + x * slabs_x_per_task[task] + (y - first_slab_x_of_task[task]))];
        }
}

#endif /* #ifdef TWODIMS */
#endif /* #ifndef GRAVITY_NOT_PERIODIC */
#endif /* #ifdef PMGRID */
