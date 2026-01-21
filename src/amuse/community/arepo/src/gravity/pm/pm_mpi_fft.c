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
 * \file        src/gravity/pm/pm_mpi_fft.c
 * \date        05/2018
 * \brief       Home-made parallel FFT transforms as needed by the code.
 * \details     We only use the one-dimensional FFTW3 routines, because the
 *              MPI versions of FFTW3 allocate memory for themselves during the
 *              transforms (which we want to strictly avoid), and because we
 *              want to allow transforms that are so big that more than 2GB
 *              may be transferred betweeen processors.
 *
 *              contains functions:
 *                void my_slab_based_fft_init(fft_plan * plan, int NgridX,
 *                  int NgridY, int NgridZ)
 *                void my_slab_transposeA(fft_plan * plan, fft_real * field,
 *                  fft_real * scratch)
 *                void my_slab_transposeB(fft_plan * plan, fft_real * field,
 *                  fft_real * scratch)
 *                static void my_slab_transpose(void *av, void *bv, int *sx,
 *                  int *firstx, int *sy, int *firsty, int nx, int ny, int nz,
 *                  int mode)
 *                void my_slab_based_fft(fft_plan * plan, void *data,
 *                  void *workspace, int forward)
 *                void my_slab_based_fft_c2c(fft_plan * plan, void *data,
 *                  void *workspace, int forward)
 *                void my_column_based_fft_init(fft_plan * plan, int NgridX,
 *                  int NgridY, int NgridZ)
 *                void my_column_based_fft_init_c2c(fft_plan * plan,
 *                  int NgridX, int NgridY, int NgridZ)
 *                void my_fft_swap23(fft_plan * plan, fft_real * data,
 *                  fft_real * out)
 *                void my_fft_swap23back(fft_plan * plan, fft_real * data,
 *                  fft_real * out)
 *                void my_fft_swap13(fft_plan * plan, fft_real * data,
 *                  fft_real * out)
 *                void my_fft_swap13back(fft_plan * plan, fft_real * data,
 *                  fft_real * out)
 *                void my_column_based_fft(fft_plan * plan, void *data,
 *                  void *workspace, int forward)
 *                void my_column_based_fft_c2c(fft_plan * plan, void *data,
 *                  void *workspace, int forward)#
 *                static void my_fft_column_remap(fft_complex * data,
 *                  int Ndims[3], int in_firstcol, int in_ncol,
 *                  fft_complex * out, int perm[3], int out_firstcol,
 *                  int out_ncol, size_t * offset_send, size_t * offset_recv,
 *                  size_t * count_send, size_t * count_recv,
 *                  size_t just_count_flag)
 *                static void my_fft_column_transpose(fft_real * data,
 *                  int Ndims[3], int in_firstcol, int in_ncol, fft_real * out,
 *                  int perm[3], int out_firstcol, int out_ncol,
 *                  size_t * offset_send, size_t * offset_recv,
 *                  size_t * count_send, size_t * count_recv,
 *                  size_t just_count_flag)
 *                static void my_fft_column_transpose_c(fft_complex * data,
 *                  int Ndims[3], int in_firstcol, int in_ncol,
 *                  fft_complex * out, int perm[3], int out_firstcol,
 *                  int out_ncol, size_t * offset_send, size_t * offset_recv,
 *                  size_t * count_send, size_t * count_recv,
 *                  size_t just_count_flag)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 26.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../main/allvars.h"
#include "../../main/proto.h"

#if defined(PMGRID)

#ifndef FFT_COLUMN_BASED
/*! \brief Initializes slab based FFT.
 *
 *  \param[out] plan FFT plan.
 *  \param[in] NgridX Number of grid points in X direction.
 *  \param[in] NgridY Number of grid points in Y direction.
 *  \param[in] NgridZ Number of grid points in Z direction.
 *
 *  \return void
 */
void my_slab_based_fft_init(fft_plan *plan, int NgridX, int NgridY, int NgridZ)
{
  subdivide_evenly(NgridX, NTask, ThisTask, &plan->slabstart_x, &plan->nslab_x);
  subdivide_evenly(NgridY, NTask, ThisTask, &plan->slabstart_y, &plan->nslab_y);

  plan->slab_to_task = (int *)mymalloc("slab_to_task", NgridX * sizeof(int));

  for(int task = 0; task < NTask; task++)
    {
      int start, n;

      subdivide_evenly(NgridX, NTask, task, &start, &n);

      for(int i = start; i < start + n; i++)
        plan->slab_to_task[i] = task;
    }

  MPI_Allreduce(&plan->nslab_x, &plan->largest_x_slab, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&plan->nslab_y, &plan->largest_y_slab, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  plan->slabs_x_per_task = (int *)mymalloc("slabs_x_per_task", NTask * sizeof(int));
  MPI_Allgather(&plan->nslab_x, 1, MPI_INT, plan->slabs_x_per_task, 1, MPI_INT, MPI_COMM_WORLD);

  plan->first_slab_x_of_task = (int *)mymalloc("first_slab_x_of_task", NTask * sizeof(int));
  MPI_Allgather(&plan->slabstart_x, 1, MPI_INT, plan->first_slab_x_of_task, 1, MPI_INT, MPI_COMM_WORLD);

  plan->slabs_y_per_task = (int *)mymalloc("slabs_y_per_task", NTask * sizeof(int));
  MPI_Allgather(&plan->nslab_y, 1, MPI_INT, plan->slabs_y_per_task, 1, MPI_INT, MPI_COMM_WORLD);

  plan->first_slab_y_of_task = (int *)mymalloc("first_slab_y_of_task", NTask * sizeof(int));
  MPI_Allgather(&plan->slabstart_y, 1, MPI_INT, plan->first_slab_y_of_task, 1, MPI_INT, MPI_COMM_WORLD);

  plan->NgridX = NgridX;
  plan->NgridY = NgridY;
  plan->NgridZ = NgridZ;

  int Ngridz = NgridZ / 2 + 1; /* dimension needed in complex space */

  plan->Ngridz = Ngridz;
  plan->Ngrid2 = 2 * Ngridz;
}

/*! \brief Transposes the array field.
 *
 *  The array field is transposed such that the data in x direction is local
 *  to only one task. This is done, so the force in x-direction can be
 *  obtained by finite differencing. However the array is not fully
 *  transposed, i.e. the x-direction is not the fastest running array index.
 *
 *  \param[in] plan FFT pan.
 *  \param[in, out] field The array to transpose.
 *  \param[out] scratch Scratch space used during communication (same size as
 *              field).
 *
 *  \return void
 */
void my_slab_transposeA(fft_plan *plan, fft_real *field, fft_real *scratch)
{
  int n, prod, task, flag_big = 0, flag_big_all = 0;

  prod = NTask * plan->nslab_x;

  for(n = 0; n < prod; n++)
    {
      int x    = n / NTask;
      int task = n % NTask;

      int y;

      for(y = plan->first_slab_y_of_task[task]; y < plan->first_slab_y_of_task[task] + plan->slabs_y_per_task[task]; y++)
        memcpy(scratch + ((size_t)plan->NgridZ) * (plan->first_slab_y_of_task[task] * plan->nslab_x +
                                                   x * plan->slabs_y_per_task[task] + (y - plan->first_slab_y_of_task[task])),
               field + ((size_t)plan->Ngrid2) * (plan->NgridY * x + y), plan->NgridZ * sizeof(fft_real));
    }

  size_t *scount = (size_t *)mymalloc("scount", NTask * sizeof(size_t));
  size_t *rcount = (size_t *)mymalloc("rcount", NTask * sizeof(size_t));
  size_t *soff   = (size_t *)mymalloc("soff", NTask * sizeof(size_t));
  size_t *roff   = (size_t *)mymalloc("roff", NTask * sizeof(size_t));

  for(task = 0; task < NTask; task++)
    {
      scount[task] = plan->nslab_x * plan->slabs_y_per_task[task] * (plan->NgridZ * sizeof(fft_real));
      rcount[task] = plan->nslab_y * plan->slabs_x_per_task[task] * (plan->NgridZ * sizeof(fft_real));

      soff[task] = plan->first_slab_y_of_task[task] * plan->nslab_x * (plan->NgridZ * sizeof(fft_real));
      roff[task] = plan->first_slab_x_of_task[task] * plan->nslab_y * (plan->NgridZ * sizeof(fft_real));

      if(scount[task] > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
        flag_big = 1;
    }

  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  myMPI_Alltoallv(scratch, scount, soff, field, rcount, roff, 1, flag_big_all, MPI_COMM_WORLD);

  myfree(roff);
  myfree(soff);
  myfree(rcount);
  myfree(scount);
}

/*! \brief Undo the transposition of the array field.
 *
 *  The transposition of the array field is undone such that the data in
 *  x direction is distributed among all tasks again. Thus the result of
 *  force computation in x-direction is sent back to the original task.
 *
 *  \param[in] plan FFT plan.
 *  \param[in, out] field The array to transpose.
 *  \param[out] scratch Scratch space used during communication (same size as
 *              field).
 *
 *  \return void
 */
void my_slab_transposeB(fft_plan *plan, fft_real *field, fft_real *scratch)
{
  int n, prod, task, flag_big = 0, flag_big_all = 0;

  size_t *scount = (size_t *)mymalloc("scount", NTask * sizeof(size_t));
  size_t *rcount = (size_t *)mymalloc("rcount", NTask * sizeof(size_t));
  size_t *soff   = (size_t *)mymalloc("soff", NTask * sizeof(size_t));
  size_t *roff   = (size_t *)mymalloc("roff", NTask * sizeof(size_t));

  for(task = 0; task < NTask; task++)
    {
      rcount[task] = plan->nslab_x * plan->slabs_y_per_task[task] * (plan->NgridZ * sizeof(fft_real));
      scount[task] = plan->nslab_y * plan->slabs_x_per_task[task] * (plan->NgridZ * sizeof(fft_real));

      roff[task] = plan->first_slab_y_of_task[task] * plan->nslab_x * (plan->NgridZ * sizeof(fft_real));
      soff[task] = plan->first_slab_x_of_task[task] * plan->nslab_y * (plan->NgridZ * sizeof(fft_real));

      if(scount[task] > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
        flag_big = 1;
    }

  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  myMPI_Alltoallv(field, scount, soff, scratch, rcount, roff, 1, flag_big_all, MPI_COMM_WORLD);

  myfree(roff);
  myfree(soff);
  myfree(rcount);
  myfree(scount);

  prod = NTask * plan->nslab_x;

  for(n = 0; n < prod; n++)
    {
      int x    = n / NTask;
      int task = n % NTask;

      int y;
      for(y = plan->first_slab_y_of_task[task]; y < plan->first_slab_y_of_task[task] + plan->slabs_y_per_task[task]; y++)
        memcpy(field + ((size_t)plan->Ngrid2) * (plan->NgridY * x + y),
               scratch + ((size_t)plan->NgridZ) * (plan->first_slab_y_of_task[task] * plan->nslab_x +
                                                   x * plan->slabs_y_per_task[task] + (y - plan->first_slab_y_of_task[task])),
               plan->NgridZ * sizeof(fft_real));
    }
}

/*  \brief Transpose a slab decomposed 3D field.
 *
 *  Given a slab-decomposed 3D field a[...] with total dimension
 *  [nx x ny x nz], whose first dimension is split across the processors, this
 *  routine outputs in b[] the transpose where then the second dimension is
 *  split across the processors. sx[] gives for each MPI task how many slabs
 *  it has, and firstx[] is the first slab for a given task. Likewise,
 *  sy[]/firsty[] gives the same thing for the transposed order. Note, the
 *  contents of the array a[] will be destroyed by the routine.
 *
 *  An element (x,y,z) is accessed in a[] with index
 *  [([x - firstx] * ny + y) * nz + z] and in b[] as
 *  [((y - firsty) * nx + x) * nz + z]
 *
 *  \param[in, out] av Pointer to array a.
 *  \param[in, out] bv Pointer to array b.
 *  \param[in] sx Array storing number of slabs in each task.
 *  \param[in] fristx Array with first slab in each task.
 *  \param[in] sy Array storing number of transposed slabs in each task.
 *  \param[in] firsty Array storing first transposed slab in each task.
 *  \param[in] nx Number of elements in x direction.
 *  \param[in] ny Number of elements in y direction.
 *  \param[in] nz Number of elements in z direction.
 *  \param[in] mode If mode = 1, the reverse operation is carried out.
 *
 *  \return void
 */
static void my_slab_transpose(void *av, void *bv, int *sx, int *firstx, int *sy, int *firsty, int nx, int ny, int nz, int mode)
{
  char *a = (char *)av;
  char *b = (char *)bv;

  size_t *scount = (size_t *)mymalloc("scount", NTask * sizeof(size_t));
  size_t *rcount = (size_t *)mymalloc("rcount", NTask * sizeof(size_t));
  size_t *soff   = (size_t *)mymalloc("soff", NTask * sizeof(size_t));
  size_t *roff   = (size_t *)mymalloc("roff", NTask * sizeof(size_t));
  int i, n, prod, flag_big = 0, flag_big_all = 0;

  for(i = 0; i < NTask; i++)
    {
      scount[i] = sy[i] * sx[ThisTask] * ((size_t)nz);
      rcount[i] = sy[ThisTask] * sx[i] * ((size_t)nz);
      soff[i]   = firsty[i] * sx[ThisTask] * ((size_t)nz);
      roff[i]   = sy[ThisTask] * firstx[i] * ((size_t)nz);

      if(scount[i] * sizeof(fft_complex) > MPI_MESSAGE_SIZELIMIT_IN_BYTES)
        flag_big = 1;
    }

  /* produce a flag if any of the send sizes is above our transfer limit, in this case we will
   * transfer the data in chunks.
   */
  MPI_Allreduce(&flag_big, &flag_big_all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if(mode == 0)
    {
      /* first pack the data into contiguous blocks */
      prod = NTask * sx[ThisTask];
      for(n = 0; n < prod; n++)
        {
          int k = n / NTask;
          int i = n % NTask;
          int j;

          for(j = 0; j < sy[i]; j++)
            memcpy(b + (k * sy[i] + j + firsty[i] * sx[ThisTask]) * (nz * sizeof(fft_complex)),
                   a + (k * ny + (firsty[i] + j)) * (nz * sizeof(fft_complex)), nz * sizeof(fft_complex));
        }

      /* tranfer the data */
      myMPI_Alltoallv(b, scount, soff, a, rcount, roff, sizeof(fft_complex), flag_big_all, MPI_COMM_WORLD);

      /* unpack the data into the right order */
      prod = NTask * sy[ThisTask];
      for(n = 0; n < prod; n++)
        {
          int j = n / NTask;
          int i = n % NTask;
          int k;

          for(k = 0; k < sx[i]; k++)
            memcpy(b + (j * nx + k + firstx[i]) * (nz * sizeof(fft_complex)),
                   a + ((k + firstx[i]) * sy[ThisTask] + j) * (nz * sizeof(fft_complex)), nz * sizeof(fft_complex));
        }
    }
  else
    {
      /* first pack the data into contiguous blocks */
      prod = NTask * sy[ThisTask];
      for(n = 0; n < prod; n++)
        {
          int j = n / NTask;
          int i = n % NTask;
          int k;

          for(k = 0; k < sx[i]; k++)
            memcpy(b + ((k + firstx[i]) * sy[ThisTask] + j) * (nz * sizeof(fft_complex)),
                   a + (j * nx + k + firstx[i]) * (nz * sizeof(fft_complex)), nz * sizeof(fft_complex));
        }

      /* tranfer the data */
      myMPI_Alltoallv(b, rcount, roff, a, scount, soff, sizeof(fft_complex), flag_big_all, MPI_COMM_WORLD);

      /* unpack the data into the right order */
      prod = NTask * sx[ThisTask];
      for(n = 0; n < prod; n++)
        {
          int k = n / NTask;
          int i = n % NTask;
          int j;

          for(j = 0; j < sy[i]; j++)
            memcpy(b + (k * ny + (firsty[i] + j)) * (nz * sizeof(fft_complex)),
                   a + (k * sy[i] + j + firsty[i] * sx[ThisTask]) * (nz * sizeof(fft_complex)), nz * sizeof(fft_complex));
        }
    }
  /* now the result is in b[] */

  myfree(roff);
  myfree(soff);
  myfree(rcount);
  myfree(scount);
}

/*! \brief Performs a slab-based Fast Fourier transformation.
 *
 *  \param[in] plan FFT plan.
 *  \param[in, out] data Array to be Fourier transformed.
 *  \param[out] workspace Workspace to temporary operate in.
 *  \param[in] forward Forward (1) or backward (-1) Fourier transformaiton?
 *
 *  \return void
 */
void my_slab_based_fft(fft_plan *plan, void *data, void *workspace, int forward)
{
  int n, prod;
  int slabsx = plan->slabs_x_per_task[ThisTask];
  int slabsy = plan->slabs_y_per_task[ThisTask];

  int ngridx  = plan->NgridX;
  int ngridy  = plan->NgridY;
  int ngridz  = plan->Ngridz;
  int ngridz2 = 2 * ngridz;

  size_t ngridx_long  = ngridx;
  size_t ngridy_long  = ngridy;
  size_t ngridz_long  = ngridz;
  size_t ngridz2_long = ngridz2;

  fft_real *data_real       = (fft_real *)data;
  fft_complex *data_complex = (fft_complex *)data, *workspace_complex = (fft_complex *)workspace;

  if(forward == 1)
    {
      /* do the z-direction FFT, real to complex */
      prod = slabsx * ngridy;
      for(n = 0; n < prod; n++)
        {
          FFTW(execute_dft_r2c)(plan->forward_plan_zdir, data_real + n * ngridz2_long, workspace_complex + n * ngridz_long);
        }

      /* do the y-direction FFT, complex to complex */
      prod = slabsx * ngridz;
      for(n = 0; n < prod; n++)
        {
          int i = n / ngridz;
          int j = n % ngridz;

          FFTW(execute_dft)
          (plan->forward_plan_ydir, workspace_complex + i * ngridz * ngridy_long + j, data_complex + i * ngridz * ngridy_long + j);
        }

      /* now our data resides in data_complex[] */

      /* do the transpose */
      my_slab_transpose(data_complex, workspace_complex, plan->slabs_x_per_task, plan->first_slab_x_of_task, plan->slabs_y_per_task,
                        plan->first_slab_y_of_task, ngridx, ngridy, ngridz, 0);

      /* now the data is in workspace_complex[] */

      /* finally, do the transform along the x-direction (we are in transposed order, x and y have interchanged */
      prod = slabsy * ngridz;
      for(n = 0; n < prod; n++)
        {
          int i = n / ngridz;
          int j = n % ngridz;

          FFTW(execute_dft)
          (plan->forward_plan_xdir, workspace_complex + i * ngridz * ngridx_long + j, data_complex + i * ngridz * ngridx_long + j);
        }

      /* now the result is in data_complex[] */
    }
  else
    {
      prod = slabsy * ngridz;

      for(n = 0; n < prod; n++)
        {
          int i = n / ngridz;
          int j = n % ngridz;

          FFTW(execute_dft)
          (plan->backward_plan_xdir, data_complex + i * ngridz * ngridx_long + j, workspace_complex + i * ngridz * ngridx_long + j);
        }

      my_slab_transpose(workspace_complex, data_complex, plan->slabs_x_per_task, plan->first_slab_x_of_task, plan->slabs_y_per_task,
                        plan->first_slab_y_of_task, ngridx, ngridy, ngridz, 1);

      prod = slabsx * ngridz;

      for(n = 0; n < prod; n++)
        {
          int i = n / ngridz;
          int j = n % ngridz;

          FFTW(execute_dft)
          (plan->backward_plan_ydir, data_complex + i * ngridz * ngridy_long + j, workspace_complex + i * ngridz * ngridy_long + j);
        }

      prod = slabsx * ngridy;

      for(n = 0; n < prod; n++)
        {
          FFTW(execute_dft_c2r)(plan->backward_plan_zdir, workspace_complex + n * ngridz_long, data_real + n * ngridz2_long);
        }

      /* now the result is in data[] */
    }
}

/*! \brief Performs a slab-based complex to complex Fast Fourier
 *         transformation.
 *
 *  \param[in] plan FFT plan.
 *  \param[in, out] data Array to be Fourier transformed.
 *  \param[out] workspace Workspace to temporary operate in.
 *  \param[in] forward Forward (1) or backward (-1) Fourier transformaiton?
 *
 *  \return void
 */
void my_slab_based_fft_c2c(fft_plan *plan, void *data, void *workspace, int forward)
{
  int n, prod;
  int slabsx = plan->slabs_x_per_task[ThisTask];
  int slabsy = plan->slabs_y_per_task[ThisTask];

  int ngridx = plan->NgridX;
  int ngridy = plan->NgridY;
  int ngridz = plan->NgridZ;

  size_t ngridx_long = ngridx;
  size_t ngridy_long = ngridy;
  size_t ngridz_long = ngridz;

  fft_complex *data_start   = (fft_complex *)data;
  fft_complex *data_complex = (fft_complex *)data, *workspace_complex = (fft_complex *)workspace;

  if(forward == 1)
    {
      /* do the z-direction FFT, complex to complex */
      prod = slabsx * ngridy;
      for(n = 0; n < prod; n++)
        {
          FFTW(execute_dft)(plan->forward_plan_zdir, data_start + n * ngridz, workspace_complex + n * ngridz);
        }

      /* do the y-direction FFT, complex to complex */
      prod = slabsx * ngridz;
      for(n = 0; n < prod; n++)
        {
          int i = n / ngridz;
          int j = n % ngridz;

          FFTW(execute_dft)
          (plan->forward_plan_ydir, workspace_complex + i * ngridz * ngridy_long + j, data_complex + i * ngridz * ngridy_long + j);
        }

      /* now our data resides in data_complex[] */

      /* do the transpose */
      my_slab_transpose(data_complex, workspace_complex, plan->slabs_x_per_task, plan->first_slab_x_of_task, plan->slabs_y_per_task,
                        plan->first_slab_y_of_task, ngridx, ngridy, ngridz, 0);

      /* now the data is in workspace_complex[] */

      /* finally, do the transform along the x-direction (we are in transposed order, x and y have interchanged */
      prod = slabsy * ngridz;
      for(n = 0; n < prod; n++)
        {
          int i = n / ngridz;
          int j = n % ngridz;

          FFTW(execute_dft)
          (plan->forward_plan_xdir, workspace_complex + i * ngridz * ngridx_long + j, data_complex + i * ngridz * ngridx_long + j);
        }

      /* now the result is in data_complex[] */
    }
  else
    {
      prod = slabsy * ngridz;

      for(n = 0; n < prod; n++)
        {
          int i = n / ngridz;
          int j = n % ngridz;

          FFTW(execute_dft)
          (plan->backward_plan_xdir, data_complex + i * ngridz * ngridx_long + j, workspace_complex + i * ngridz * ngridx_long + j);
        }

      my_slab_transpose(workspace_complex, data_complex, plan->slabs_x_per_task, plan->first_slab_x_of_task, plan->slabs_y_per_task,
                        plan->first_slab_y_of_task, ngridx, ngridy, ngridz, 1);

      prod = slabsx * ngridz;

      for(n = 0; n < prod; n++)
        {
          int i = n / ngridz;
          int j = n % ngridz;

          FFTW(execute_dft)
          (plan->backward_plan_ydir, data_complex + i * ngridz * ngridy_long + j, workspace_complex + i * ngridz * ngridy_long + j);
        }

      prod = slabsx * ngridy;

      for(n = 0; n < prod; n++)
        {
          FFTW(execute_dft)(plan->backward_plan_zdir, workspace_complex + n * ngridz, data_start + n * ngridz);
        }

      /* now the result is in data[] */
    }
}

#else /* #ifndef FFT_COLUMN_BASED */

static void my_fft_column_remap(fft_complex *data, int Ndims[3], int in_firstcol, int in_ncol, fft_complex *out, int perm[3],
                                int out_firstcol, int out_ncol, size_t *offset_send, size_t *offset_recv, size_t *count_send,
                                size_t *count_recv, size_t just_count_flag);

static void my_fft_column_transpose(fft_real *data, int Ndims[3], /* global dimensions of data cube */
                                    int in_firstcol, int in_ncol, /* first column and number of columns */
                                    fft_real *out, int perm[3], int out_firstcol, int out_ncol, size_t *offset_send,
                                    size_t *offset_recv, size_t *count_send, size_t *count_recv, size_t just_count_flag);

static void my_fft_column_transpose_c(fft_complex *data, int Ndims[3], /* global dimensions of data cube */
                                      int in_firstcol, int in_ncol,    /* first column and number of columns */
                                      fft_complex *out, int perm[3], int out_firstcol, int out_ncol, size_t *offset_send,
                                      size_t *offset_recv, size_t *count_send, size_t *count_recv, size_t just_count_flag);

/*! \brief Initializes column based FFT.
 *
 *  \param[out] plan FFT plan.
 *  \param[in] NgridX Number of grid points in X direction.
 *  \param[in] NgridY Number of grid points in Y direction.
 *  \param[in] NgridZ Number of grid points in Z direction.
 *
 *  \return void
 */
void my_column_based_fft_init(fft_plan *plan, int NgridX, int NgridY, int NgridZ)
{
  plan->NgridX = NgridX;
  plan->NgridY = NgridY;
  plan->NgridZ = NgridZ;

  int Ngridz = NgridZ / 2 + 1;

  plan->Ngridz = Ngridz;
  plan->Ngrid2 = 2 * Ngridz;

  int columns, avg, exc, tasklastsection, pivotcol;

  columns         = NgridX * NgridY;
  avg             = (columns - 1) / NTask + 1;
  exc             = NTask * avg - columns;
  tasklastsection = NTask - exc;
  pivotcol        = tasklastsection * avg;

  plan->pivotcol        = pivotcol;
  plan->avg             = avg;
  plan->tasklastsection = tasklastsection;

  if(ThisTask < tasklastsection)
    {
      plan->base_firstcol = ThisTask * avg;
      plan->base_ncol     = avg;
    }
  else
    {
      plan->base_firstcol = ThisTask * avg - (ThisTask - tasklastsection);
      plan->base_ncol     = avg - 1;
    }

  plan->base_lastcol = plan->base_firstcol + plan->base_ncol - 1;

  subdivide_evenly(NgridX * Ngridz, NTask, ThisTask, &plan->transposed_firstcol, &plan->transposed_ncol);

  subdivide_evenly(NgridY * Ngridz, NTask, ThisTask, &plan->second_transposed_firstcol, &plan->second_transposed_ncol);

  subdivide_evenly(plan->NgridX * plan->Ngrid2, NTask, ThisTask, &plan->firstcol_XZ, &plan->ncol_XZ);

  subdivide_evenly(plan->NgridY * plan->Ngrid2, NTask, ThisTask, &plan->firstcol_YZ, &plan->ncol_YZ);

  plan->second_transposed_ncells = ((size_t)plan->NgridX) * plan->second_transposed_ncol;

  plan->max_datasize = ((size_t)plan->Ngrid2) * plan->base_ncol;
  plan->max_datasize = smax(plan->max_datasize, 2 * ((size_t)plan->NgridY) * plan->transposed_ncol);
  plan->max_datasize = smax(plan->max_datasize, 2 * ((size_t)plan->NgridX) * plan->second_transposed_ncol);
  plan->max_datasize = smax(plan->max_datasize, ((size_t)plan->ncol_XZ) * plan->NgridY);
  plan->max_datasize = smax(plan->max_datasize, ((size_t)plan->ncol_YZ) * plan->NgridX);

  plan->fftsize = plan->max_datasize;

  plan->offsets_send_A      = mymalloc_clear("offsets_send_A", NTask * sizeof(size_t));
  plan->offsets_recv_A      = mymalloc_clear("offsets_recv_A", NTask * sizeof(size_t));
  plan->offsets_send_B      = mymalloc_clear("offsets_send_B", NTask * sizeof(size_t));
  plan->offsets_recv_B      = mymalloc_clear("offsets_recv_B", NTask * sizeof(size_t));
  plan->offsets_send_C      = mymalloc_clear("offsets_send_C", NTask * sizeof(size_t));
  plan->offsets_recv_C      = mymalloc_clear("offsets_recv_C", NTask * sizeof(size_t));
  plan->offsets_send_D      = mymalloc_clear("offsets_send_D", NTask * sizeof(size_t));
  plan->offsets_recv_D      = mymalloc_clear("offsets_recv_D", NTask * sizeof(size_t));
  plan->offsets_send_13     = mymalloc_clear("offsets_send_13", NTask * sizeof(size_t));
  plan->offsets_recv_13     = mymalloc_clear("offsets_recv_13", NTask * sizeof(size_t));
  plan->offsets_send_23     = mymalloc_clear("offsets_send_23", NTask * sizeof(size_t));
  plan->offsets_recv_23     = mymalloc_clear("offsets_recv_23", NTask * sizeof(size_t));
  plan->offsets_send_13back = mymalloc_clear("offsets_send_13back", NTask * sizeof(size_t));
  plan->offsets_recv_13back = mymalloc_clear("offsets_recv_13back", NTask * sizeof(size_t));
  plan->offsets_send_23back = mymalloc_clear("offsets_send_23back", NTask * sizeof(size_t));
  plan->offsets_recv_23back = mymalloc_clear("offsets_recv_23back", NTask * sizeof(size_t));

  plan->count_send_A      = mymalloc_clear("count_send_A", NTask * sizeof(size_t));
  plan->count_recv_A      = mymalloc_clear("count_recv_A", NTask * sizeof(size_t));
  plan->count_send_B      = mymalloc_clear("count_send_B", NTask * sizeof(size_t));
  plan->count_recv_B      = mymalloc_clear("count_recv_B", NTask * sizeof(size_t));
  plan->count_send_C      = mymalloc_clear("count_send_C", NTask * sizeof(size_t));
  plan->count_recv_C      = mymalloc_clear("count_recv_C", NTask * sizeof(size_t));
  plan->count_send_D      = mymalloc_clear("count_send_D", NTask * sizeof(size_t));
  plan->count_recv_D      = mymalloc_clear("count_recv_D", NTask * sizeof(size_t));
  plan->count_send_13     = mymalloc_clear("count_send_13", NTask * sizeof(size_t));
  plan->count_recv_13     = mymalloc_clear("count_recv_13", NTask * sizeof(size_t));
  plan->count_send_23     = mymalloc_clear("count_send_23", NTask * sizeof(size_t));
  plan->count_recv_23     = mymalloc_clear("count_recv_23", NTask * sizeof(size_t));
  plan->count_send_13back = mymalloc_clear("count_send_13back", NTask * sizeof(size_t));
  plan->count_recv_13back = mymalloc_clear("count_recv_13back", NTask * sizeof(size_t));
  plan->count_send_23back = mymalloc_clear("count_send_23back", NTask * sizeof(size_t));
  plan->count_recv_23back = mymalloc_clear("count_recv_23back", NTask * sizeof(size_t));

  int dimA[3]  = {plan->NgridX, plan->NgridY, plan->Ngridz};
  int permA[3] = {0, 2, 1};

  my_fft_column_remap(NULL, dimA, plan->base_firstcol, plan->base_ncol, NULL, permA, plan->transposed_firstcol, plan->transposed_ncol,
                      plan->offsets_send_A, plan->offsets_recv_A, plan->count_send_A, plan->count_recv_A, 1);

  int dimB[3]  = {plan->NgridX, plan->Ngridz, plan->NgridY};
  int permB[3] = {2, 1, 0};

  my_fft_column_remap(NULL, dimB, plan->transposed_firstcol, plan->transposed_ncol, NULL, permB, plan->second_transposed_firstcol,
                      plan->second_transposed_ncol, plan->offsets_send_B, plan->offsets_recv_B, plan->count_send_B, plan->count_recv_B,
                      1);

  int dimC[3]  = {plan->NgridY, plan->Ngridz, plan->NgridX};
  int permC[3] = {2, 1, 0};

  my_fft_column_remap(NULL, dimC, plan->second_transposed_firstcol, plan->second_transposed_ncol, NULL, permC,
                      plan->transposed_firstcol, plan->transposed_ncol, plan->offsets_send_C, plan->offsets_recv_C, plan->count_send_C,
                      plan->count_recv_C, 1);

  int dimD[3]  = {plan->NgridX, plan->Ngridz, plan->NgridY};
  int permD[3] = {0, 2, 1};

  my_fft_column_remap(NULL, dimD, plan->transposed_firstcol, plan->transposed_ncol, NULL, permD, plan->base_firstcol, plan->base_ncol,
                      plan->offsets_send_D, plan->offsets_recv_D, plan->count_send_D, plan->count_recv_D, 1);

  int dim23[3]  = {plan->NgridX, plan->NgridY, plan->Ngrid2};
  int perm23[3] = {0, 2, 1};

  my_fft_column_transpose(NULL, dim23, plan->base_firstcol, plan->base_ncol, NULL, perm23, plan->firstcol_XZ, plan->ncol_XZ,
                          plan->offsets_send_23, plan->offsets_recv_23, plan->count_send_23, plan->count_recv_23, 1);

  int dim23back[3]  = {plan->NgridX, plan->Ngrid2, plan->NgridY};
  int perm23back[3] = {0, 2, 1};

  my_fft_column_transpose(NULL, dim23back, plan->firstcol_XZ, plan->ncol_XZ, NULL, perm23back, plan->base_firstcol, plan->base_ncol,
                          plan->offsets_send_23back, plan->offsets_recv_23back, plan->count_send_23back, plan->count_recv_23back, 1);

  int dim13[3]  = {plan->NgridX, plan->NgridY, plan->Ngrid2};
  int perm13[3] = {2, 1, 0};

  my_fft_column_transpose(NULL, dim13, plan->base_firstcol, plan->base_ncol, NULL, perm13, plan->firstcol_YZ, plan->ncol_YZ,
                          plan->offsets_send_13, plan->offsets_recv_13, plan->count_send_13, plan->count_recv_13, 1);

  int dim13back[3]  = {plan->Ngrid2, plan->NgridY, plan->NgridX};
  int perm13back[3] = {2, 1, 0};

  my_fft_column_transpose(NULL, dim13back, plan->firstcol_YZ, plan->ncol_YZ, NULL, perm13back, plan->base_firstcol, plan->base_ncol,
                          plan->offsets_send_13back, plan->offsets_recv_13back, plan->count_send_13back, plan->count_recv_13back, 1);
}

/*! \brief Initializes complex to complex column based FFT.
 *
 *  \param[out] plan FFT plan.
 *  \param[in] NgridX Number of grid points in X direction.
 *  \param[in] NgridY Number of grid points in Y direction.
 *  \param[in] NgridZ Number of grid points in Z direction.
 *
 *  \return void
 */
void my_column_based_fft_init_c2c(fft_plan *plan, int NgridX, int NgridY, int NgridZ)
{
  plan->NgridX = NgridX;
  plan->NgridY = NgridY;
  plan->NgridZ = NgridZ;

  int columns, avg, exc, tasklastsection, pivotcol;

  columns         = NgridX * NgridY;
  avg             = (columns - 1) / NTask + 1;
  exc             = NTask * avg - columns;
  tasklastsection = NTask - exc;
  pivotcol        = tasklastsection * avg;

  plan->pivotcol        = pivotcol;
  plan->avg             = avg;
  plan->tasklastsection = tasklastsection;

  if(ThisTask < tasklastsection)
    {
      plan->base_firstcol = ThisTask * avg;
      plan->base_ncol     = avg;
    }
  else
    {
      plan->base_firstcol = ThisTask * avg - (ThisTask - tasklastsection);
      plan->base_ncol     = avg - 1;
    }

  plan->base_lastcol = plan->base_firstcol + plan->base_ncol - 1;

  subdivide_evenly(NgridX * NgridZ, NTask, ThisTask, &plan->transposed_firstcol, &plan->transposed_ncol);

  subdivide_evenly(NgridY * NgridZ, NTask, ThisTask, &plan->second_transposed_firstcol, &plan->second_transposed_ncol);

  subdivide_evenly(plan->NgridX * plan->NgridZ, NTask, ThisTask, &plan->firstcol_XZ, &plan->ncol_XZ);

  subdivide_evenly(plan->NgridY * plan->NgridZ, NTask, ThisTask, &plan->firstcol_YZ, &plan->ncol_YZ);

  plan->second_transposed_ncells = ((size_t)plan->NgridX) * plan->second_transposed_ncol;

  plan->max_datasize = 2 * ((size_t)plan->NgridZ) * plan->base_ncol;
  plan->max_datasize = smax(plan->max_datasize, 2 * ((size_t)plan->NgridY) * plan->transposed_ncol);
  plan->max_datasize = smax(plan->max_datasize, 2 * ((size_t)plan->NgridX) * plan->second_transposed_ncol);
  plan->max_datasize = smax(plan->max_datasize, ((size_t)plan->ncol_XZ) * plan->NgridY);
  plan->max_datasize = smax(plan->max_datasize, ((size_t)plan->ncol_YZ) * plan->NgridX);

  plan->fftsize = plan->max_datasize;

  plan->offsets_send_A      = mymalloc_clear("offsets_send_A", NTask * sizeof(size_t));
  plan->offsets_recv_A      = mymalloc_clear("offsets_recv_A", NTask * sizeof(size_t));
  plan->offsets_send_B      = mymalloc_clear("offsets_send_B", NTask * sizeof(size_t));
  plan->offsets_recv_B      = mymalloc_clear("offsets_recv_B", NTask * sizeof(size_t));
  plan->offsets_send_C      = mymalloc_clear("offsets_send_C", NTask * sizeof(size_t));
  plan->offsets_recv_C      = mymalloc_clear("offsets_recv_C", NTask * sizeof(size_t));
  plan->offsets_send_D      = mymalloc_clear("offsets_send_D", NTask * sizeof(size_t));
  plan->offsets_recv_D      = mymalloc_clear("offsets_recv_D", NTask * sizeof(size_t));
  plan->offsets_send_13     = mymalloc_clear("offsets_send_13", NTask * sizeof(size_t));
  plan->offsets_recv_13     = mymalloc_clear("offsets_recv_13", NTask * sizeof(size_t));
  plan->offsets_send_23     = mymalloc_clear("offsets_send_23", NTask * sizeof(size_t));
  plan->offsets_recv_23     = mymalloc_clear("offsets_recv_23", NTask * sizeof(size_t));
  plan->offsets_send_13back = mymalloc_clear("offsets_send_13back", NTask * sizeof(size_t));
  plan->offsets_recv_13back = mymalloc_clear("offsets_recv_13back", NTask * sizeof(size_t));
  plan->offsets_send_23back = mymalloc_clear("offsets_send_23back", NTask * sizeof(size_t));
  plan->offsets_recv_23back = mymalloc_clear("offsets_recv_23back", NTask * sizeof(size_t));

  plan->count_send_A      = mymalloc_clear("count_send_A", NTask * sizeof(size_t));
  plan->count_recv_A      = mymalloc_clear("count_recv_A", NTask * sizeof(size_t));
  plan->count_send_B      = mymalloc_clear("count_send_B", NTask * sizeof(size_t));
  plan->count_recv_B      = mymalloc_clear("count_recv_B", NTask * sizeof(size_t));
  plan->count_send_C      = mymalloc_clear("count_send_C", NTask * sizeof(size_t));
  plan->count_recv_C      = mymalloc_clear("count_recv_C", NTask * sizeof(size_t));
  plan->count_send_D      = mymalloc_clear("count_send_D", NTask * sizeof(size_t));
  plan->count_recv_D      = mymalloc_clear("count_recv_D", NTask * sizeof(size_t));
  plan->count_send_13     = mymalloc_clear("count_send_13", NTask * sizeof(size_t));
  plan->count_recv_13     = mymalloc_clear("count_recv_13", NTask * sizeof(size_t));
  plan->count_send_23     = mymalloc_clear("count_send_23", NTask * sizeof(size_t));
  plan->count_recv_23     = mymalloc_clear("count_recv_23", NTask * sizeof(size_t));
  plan->count_send_13back = mymalloc_clear("count_send_13back", NTask * sizeof(size_t));
  plan->count_recv_13back = mymalloc_clear("count_recv_13back", NTask * sizeof(size_t));
  plan->count_send_23back = mymalloc_clear("count_send_23back", NTask * sizeof(size_t));
  plan->count_recv_23back = mymalloc_clear("count_recv_23back", NTask * sizeof(size_t));

  int dimA[3]  = {plan->NgridX, plan->NgridY, plan->NgridZ};
  int permA[3] = {0, 2, 1};

  my_fft_column_remap(NULL, dimA, plan->base_firstcol, plan->base_ncol, NULL, permA, plan->transposed_firstcol, plan->transposed_ncol,
                      plan->offsets_send_A, plan->offsets_recv_A, plan->count_send_A, plan->count_recv_A, 1);

  int dimB[3]  = {plan->NgridX, plan->NgridZ, plan->NgridY};
  int permB[3] = {2, 1, 0};

  my_fft_column_remap(NULL, dimB, plan->transposed_firstcol, plan->transposed_ncol, NULL, permB, plan->second_transposed_firstcol,
                      plan->second_transposed_ncol, plan->offsets_send_B, plan->offsets_recv_B, plan->count_send_B, plan->count_recv_B,
                      1);

  int dimC[3]  = {plan->NgridY, plan->NgridZ, plan->NgridX};
  int permC[3] = {2, 1, 0};

  my_fft_column_remap(NULL, dimC, plan->second_transposed_firstcol, plan->second_transposed_ncol, NULL, permC,
                      plan->transposed_firstcol, plan->transposed_ncol, plan->offsets_send_C, plan->offsets_recv_C, plan->count_send_C,
                      plan->count_recv_C, 1);

  int dimD[3]  = {plan->NgridX, plan->NgridZ, plan->NgridY};
  int permD[3] = {0, 2, 1};

  my_fft_column_remap(NULL, dimD, plan->transposed_firstcol, plan->transposed_ncol, NULL, permD, plan->base_firstcol, plan->base_ncol,
                      plan->offsets_send_D, plan->offsets_recv_D, plan->count_send_D, plan->count_recv_D, 1);

  int dim23[3]  = {plan->NgridX, plan->NgridY, plan->NgridZ};
  int perm23[3] = {0, 2, 1};

  my_fft_column_transpose_c(NULL, dim23, plan->base_firstcol, plan->base_ncol, NULL, perm23, plan->firstcol_XZ, plan->ncol_XZ,
                            plan->offsets_send_23, plan->offsets_recv_23, plan->count_send_23, plan->count_recv_23, 1);

  int dim23back[3]  = {plan->NgridX, plan->NgridZ, plan->NgridY};
  int perm23back[3] = {0, 2, 1};

  my_fft_column_transpose_c(NULL, dim23back, plan->firstcol_XZ, plan->ncol_XZ, NULL, perm23back, plan->base_firstcol, plan->base_ncol,
                            plan->offsets_send_23back, plan->offsets_recv_23back, plan->count_send_23back, plan->count_recv_23back, 1);

  int dim13[3]  = {plan->NgridX, plan->NgridY, plan->NgridZ};
  int perm13[3] = {2, 1, 0};

  my_fft_column_transpose_c(NULL, dim13, plan->base_firstcol, plan->base_ncol, NULL, perm13, plan->firstcol_YZ, plan->ncol_YZ,
                            plan->offsets_send_13, plan->offsets_recv_13, plan->count_send_13, plan->count_recv_13, 1);

  int dim13back[3]  = {plan->NgridZ, plan->NgridY, plan->NgridX};
  int perm13back[3] = {2, 1, 0};

  my_fft_column_transpose_c(NULL, dim13back, plan->firstcol_YZ, plan->ncol_YZ, NULL, perm13back, plan->base_firstcol, plan->base_ncol,
                            plan->offsets_send_13back, plan->offsets_recv_13back, plan->count_send_13back, plan->count_recv_13back, 1);
}

/*! \brief YZ column transpose.
 *
 *  \param[in] plan FFT plan.
 *  \param[in] data Array with data to be swapped.
 *  \param[out] out Array with data output.
 *
 *  \return void
 */
void my_fft_swap23(fft_plan *plan, fft_real *data, fft_real *out)
{
  int dim23[3]  = {plan->NgridX, plan->NgridY, plan->Ngrid2};
  int perm23[3] = {0, 2, 1};

  my_fft_column_transpose(data, dim23, plan->base_firstcol, plan->base_ncol, out, perm23, plan->firstcol_XZ, plan->ncol_XZ,
                          plan->offsets_send_23, plan->offsets_recv_23, plan->count_send_23, plan->count_recv_23, 0);
}

/*! \brief Reverse YZ column transpose.
 *
 *  \param[in] plan FFT plan.
 *  \param[in] data Array with data to be swapped.
 *  \param[out] out Array with data output.
 *
 *  \return void
 */
void my_fft_swap23back(fft_plan *plan, fft_real *data, fft_real *out)
{
  int dim23back[3]  = {plan->NgridX, plan->Ngrid2, plan->NgridY};
  int perm23back[3] = {0, 2, 1};

  my_fft_column_transpose(data, dim23back, plan->firstcol_XZ, plan->ncol_XZ, out, perm23back, plan->base_firstcol, plan->base_ncol,
                          plan->offsets_send_23back, plan->offsets_recv_23back, plan->count_send_23back, plan->count_recv_23back, 0);
}

/*! \brief XZ column transpose.
 *
 *  \param[in] plan FFT plan.
 *  \param[in] data Array with data to be swapped.
 *  \param[out] out Array with data output.
 *
 *  \return void
 */
void my_fft_swap13(fft_plan *plan, fft_real *data, fft_real *out)
{
  int dim13[3]  = {plan->NgridX, plan->NgridY, plan->Ngrid2};
  int perm13[3] = {2, 1, 0};

  my_fft_column_transpose(data, dim13, plan->base_firstcol, plan->base_ncol, out, perm13, plan->firstcol_YZ, plan->ncol_YZ,
                          plan->offsets_send_13, plan->offsets_recv_13, plan->count_send_13, plan->count_recv_13, 0);
}

/*! \brief Reverse XZ column transpose.
 *
 *  \param[in] plan FFT plan.
 *  \param[in] data Array with data to be swapped.
 *  \param[out] out Array with data output.
 *
 *  \return void
 */
void my_fft_swap13back(fft_plan *plan, fft_real *data, fft_real *out)
{
  int dim13back[3]  = {plan->Ngrid2, plan->NgridY, plan->NgridX};
  int perm13back[3] = {2, 1, 0};

  my_fft_column_transpose(data, dim13back, plan->firstcol_YZ, plan->ncol_YZ, out, perm13back, plan->base_firstcol, plan->base_ncol,
                          plan->offsets_send_13back, plan->offsets_recv_13back, plan->count_send_13back, plan->count_recv_13back, 0);
}

/*! \brief Performs a column-based Fast Fourier transformation.
 *
 *  \param[in] plan FFT plan.
 *  \param[in, out] data Array to be Fourier transformed.
 *  \param[out] workspace Workspace to temporary operate in.
 *  \param[in] forward Forward (1) or backward (-1) Fourier transformaiton?
 *
 *  \return void
 */
void my_column_based_fft(fft_plan *plan, void *data, void *workspace, int forward)
{
  size_t n;
  fft_real *data_real = data, *workspace_real = workspace;
  fft_complex *data_complex = data, *workspace_complex = workspace;

  if(forward == 1)
    {
      /* do the z-direction FFT, real to complex */
      for(n = 0; n < plan->base_ncol; n++)
        FFTW(execute_dft_r2c)(plan->forward_plan_zdir, data_real + n * plan->Ngrid2, workspace_complex + n * plan->Ngridz);

      int dimA[3]  = {plan->NgridX, plan->NgridY, plan->Ngridz};
      int permA[3] = {0, 2, 1};

      my_fft_column_remap(workspace_complex, dimA, plan->base_firstcol, plan->base_ncol, data_complex, permA,
                          plan->transposed_firstcol, plan->transposed_ncol, plan->offsets_send_A, plan->offsets_recv_A,
                          plan->count_send_A, plan->count_recv_A, 0);

      /* do the y-direction FFT in 'data', complex to complex */
      for(n = 0; n < plan->transposed_ncol; n++)
        FFTW(execute_dft)(plan->forward_plan_ydir, data_complex + n * plan->NgridY, workspace_complex + n * plan->NgridY);

      int dimB[3]  = {plan->NgridX, plan->Ngridz, plan->NgridY};
      int permB[3] = {2, 1, 0};

      my_fft_column_remap(workspace_complex, dimB, plan->transposed_firstcol, plan->transposed_ncol, data_complex, permB,
                          plan->second_transposed_firstcol, plan->second_transposed_ncol, plan->offsets_send_B, plan->offsets_recv_B,
                          plan->count_send_B, plan->count_recv_B, 0);

      /* do the x-direction FFT in 'data', complex to complex */
      for(n = 0; n < plan->second_transposed_ncol; n++)
        FFTW(execute_dft)(plan->forward_plan_xdir, data_complex + n * plan->NgridX, workspace_complex + n * plan->NgridX);

      /* result is now in workspace */
    }
  else
    {
      /* do inverse FFT in 'data' */
      for(n = 0; n < plan->second_transposed_ncol; n++)
        FFTW(execute_dft)(plan->backward_plan_xdir, data_complex + n * plan->NgridX, workspace_complex + n * plan->NgridX);

      int dimC[3]  = {plan->NgridY, plan->Ngridz, plan->NgridX};
      int permC[3] = {2, 1, 0};

      my_fft_column_remap(workspace_complex, dimC, plan->second_transposed_firstcol, plan->second_transposed_ncol, data_complex, permC,
                          plan->transposed_firstcol, plan->transposed_ncol, plan->offsets_send_C, plan->offsets_recv_C,
                          plan->count_send_C, plan->count_recv_C, 0);

      /* do inverse FFT in 'data' */
      for(n = 0; n < plan->transposed_ncol; n++)
        FFTW(execute_dft)(plan->backward_plan_ydir, data_complex + n * plan->NgridY, workspace_complex + n * plan->NgridY);

      int dimD[3]  = {plan->NgridX, plan->Ngridz, plan->NgridY};
      int permD[3] = {0, 2, 1};

      my_fft_column_remap(workspace_complex, dimD, plan->transposed_firstcol, plan->transposed_ncol, data_complex, permD,
                          plan->base_firstcol, plan->base_ncol, plan->offsets_send_D, plan->offsets_recv_D, plan->count_send_D,
                          plan->count_recv_D, 0);

      /* do complex-to-real inverse transform on z-coordinates */
      for(n = 0; n < plan->base_ncol; n++)
        FFTW(execute_dft_c2r)(plan->backward_plan_zdir, data_complex + n * plan->Ngridz, workspace_real + n * plan->Ngrid2);
    }
}

/*! \brief Performs a slab-based complex to complex Fast Fourier
 *         transformation.
 *
 *  \param[in] plan FFT plan.
 *  \param[in, out] data Array to be Fourier transformed.
 *  \param[out] workspace Workspace to temporary operate in.
 *  \param[in] forward Forward (1) or backward (-1) Fourier transformaiton?
 *
 *  \return void
 */
void my_column_based_fft_c2c(fft_plan *plan, void *data, void *workspace, int forward)
{
  size_t n;
  fft_complex *data_complex = data, *workspace_complex = workspace;

  if(forward == 1)
    {
      /* do the z-direction FFT, complex to complex */
      for(n = 0; n < plan->base_ncol; n++)
        FFTW(execute_dft)(plan->forward_plan_zdir, data_complex + n * plan->NgridZ, workspace_complex + n * plan->NgridZ);

      int dimA[3]  = {plan->NgridX, plan->NgridY, plan->NgridZ};
      int permA[3] = {0, 2, 1};

      my_fft_column_remap(workspace_complex, dimA, plan->base_firstcol, plan->base_ncol, data_complex, permA,
                          plan->transposed_firstcol, plan->transposed_ncol, plan->offsets_send_A, plan->offsets_recv_A,
                          plan->count_send_A, plan->count_recv_A, 0);

      /* do the y-direction FFT in 'data', complex to complex */
      for(n = 0; n < plan->transposed_ncol; n++)
        FFTW(execute_dft)(plan->forward_plan_ydir, data_complex + n * plan->NgridY, workspace_complex + n * plan->NgridY);

      int dimB[3]  = {plan->NgridX, plan->NgridZ, plan->NgridY};
      int permB[3] = {2, 1, 0};

      my_fft_column_remap(workspace_complex, dimB, plan->transposed_firstcol, plan->transposed_ncol, data_complex, permB,
                          plan->second_transposed_firstcol, plan->second_transposed_ncol, plan->offsets_send_B, plan->offsets_recv_B,
                          plan->count_send_B, plan->count_recv_B, 0);

      /* do the x-direction FFT in 'data', complex to complex */
      for(n = 0; n < plan->second_transposed_ncol; n++)
        FFTW(execute_dft)(plan->forward_plan_xdir, data_complex + n * plan->NgridX, workspace_complex + n * plan->NgridX);

      /* result is now in workspace */
    }
  else
    {
      /* do inverse FFT in 'data' */
      for(n = 0; n < plan->second_transposed_ncol; n++)
        FFTW(execute_dft)(plan->backward_plan_xdir, data_complex + n * plan->NgridX, workspace_complex + n * plan->NgridX);

      int dimC[3]  = {plan->NgridY, plan->NgridZ, plan->NgridX};
      int permC[3] = {2, 1, 0};

      my_fft_column_remap(workspace_complex, dimC, plan->second_transposed_firstcol, plan->second_transposed_ncol, data_complex, permC,
                          plan->transposed_firstcol, plan->transposed_ncol, plan->offsets_send_C, plan->offsets_recv_C,
                          plan->count_send_C, plan->count_recv_C, 0);

      /* do inverse FFT in 'data' */
      for(n = 0; n < plan->transposed_ncol; n++)
        FFTW(execute_dft)(plan->backward_plan_ydir, data_complex + n * plan->NgridY, workspace_complex + n * plan->NgridY);

      int dimD[3]  = {plan->NgridX, plan->NgridZ, plan->NgridY};
      int permD[3] = {0, 2, 1};

      my_fft_column_remap(workspace_complex, dimD, plan->transposed_firstcol, plan->transposed_ncol, data_complex, permD,
                          plan->base_firstcol, plan->base_ncol, plan->offsets_send_D, plan->offsets_recv_D, plan->count_send_D,
                          plan->count_recv_D, 0);

      /* do complex-to-complex inverse transform on z-coordinates */
      for(n = 0; n < plan->base_ncol; n++)
        FFTW(execute_dft)(plan->backward_plan_zdir, data_complex + n * plan->NgridZ, workspace_complex + n * plan->NgridZ);
    }
}

/*! \brief Remaps column-based FFT data.
 *
 *  \param[in] data Data to be transposed.
 *  \param[in] Ndims Global number of dimensions of data cube.
 *  \param[in] in_firstcol First column.
 *  \param[in] in_ncol Number of columns.
 *  \param[out] out Data output.
 *  \param[in] perm Permutations in dimensions.
 *  \param[out] out_firstcol First column in output data.
 *  \param[out] out_ncol Number of columns in output data.
 *  \param[out] offset_send Offset in array for send operation to MPI tasks.
 *  \param[out] offset_recv Offset in array for receive operation from MPI
 *              tasks.
 *  \param[out] count_send Count how many elements have to be sent to each
 *              MPI task.
 *  \param[out] count_recv Count how many elements have to be received from
 *              each MPI task.
 *  \param[in] just_count_flag Do element counting for communication instead
 *             of data transfer.
 *
 *  \return void
 */
static void my_fft_column_remap(fft_complex *data, int Ndims[3], int in_firstcol, int in_ncol, fft_complex *out, int perm[3],
                                int out_firstcol, int out_ncol, size_t *offset_send, size_t *offset_recv, size_t *count_send,
                                size_t *count_recv, size_t just_count_flag)
{
  int j, target, origin, ngrp, recvTask, perm_rev[3], xyz[3], uvw[3];
  size_t nimport, nexport;

  /* determine the inverse permutation */
  for(j = 0; j < 3; j++)
    perm_rev[j] = perm[j];

  if(!(perm_rev[perm[0]] == 0 && perm_rev[perm[1]] == 1 && perm_rev[perm[2]] == 2)) /* not yet the inverse */
    {
      for(j = 0; j < 3; j++)
        perm_rev[j] = perm[perm[j]];

      if(!(perm_rev[perm[0]] == 0 && perm_rev[perm[1]] == 1 && perm_rev[perm[2]] == 2))
        terminate("bummer");
    }

  int in_colums          = Ndims[0] * Ndims[1];
  int in_avg             = (in_colums - 1) / NTask + 1;
  int in_exc             = NTask * in_avg - in_colums;
  int in_tasklastsection = NTask - in_exc;
  int in_pivotcol        = in_tasklastsection * in_avg;

  int out_colums          = Ndims[perm[0]] * Ndims[perm[1]];
  int out_avg             = (out_colums - 1) / NTask + 1;
  int out_exc             = NTask * out_avg - out_colums;
  int out_tasklastsection = NTask - out_exc;
  int out_pivotcol        = out_tasklastsection * out_avg;

  size_t i, ncells = ((size_t)in_ncol) * Ndims[2];

  xyz[0] = in_firstcol / Ndims[1];
  xyz[1] = in_firstcol % Ndims[1];
  xyz[2] = 0;

  memset(count_send, 0, NTask * sizeof(size_t));

  /* loop over all cells in input array and determine target processor */
  for(i = 0; i < ncells; i++)
    {
      /* determine target task */
      uvw[0] = xyz[perm[0]];
      uvw[1] = xyz[perm[1]];
      uvw[2] = xyz[perm[2]];

      int newcol = Ndims[perm[1]] * uvw[0] + uvw[1];
      if(newcol < out_pivotcol)
        target = newcol / out_avg;
      else
        target = (newcol - out_pivotcol) / (out_avg - 1) + out_tasklastsection;

      /* move data element to targettask */

      if(just_count_flag)
        count_send[target]++;
      else
        {
          size_t off  = offset_send[target] + count_send[target]++;
          out[off][0] = data[i][0];
          out[off][1] = data[i][1];
        }
      xyz[2]++;
      if(xyz[2] == Ndims[2])
        {
          xyz[2] = 0;
          xyz[1]++;
          if(xyz[1] == Ndims[1])
            {
              xyz[1] = 0;
              xyz[0]++;
            }
        }
    }

  if(just_count_flag)
    {
      MPI_Alltoall(count_send, sizeof(size_t), MPI_BYTE, count_recv, sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);

      for(j = 0, nimport = 0, nexport = 0, offset_send[0] = 0, offset_recv[0] = 0; j < NTask; j++)
        {
          nexport += count_send[j];
          nimport += count_recv[j];

          if(j > 0)
            {
              offset_send[j] = offset_send[j - 1] + count_send[j - 1];
              offset_recv[j] = offset_recv[j - 1] + count_recv[j - 1];
            }
        }

      if(nexport != ncells)
        terminate("nexport=%lld != ncells=%lld", (long long)nexport, (long long)ncells);
    }
  else
    {
      nimport = 0;

      /* exchange all the data */
      for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;

          if(recvTask < NTask)
            {
              if(count_send[recvTask] > 0 || count_recv[recvTask] > 0)
                myMPI_Sendrecv(&out[offset_send[recvTask]], count_send[recvTask] * sizeof(fft_complex), MPI_BYTE, recvTask, TAG_DENS_A,
                               &data[offset_recv[recvTask]], count_recv[recvTask] * sizeof(fft_complex), MPI_BYTE, recvTask,
                               TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              nimport += count_recv[recvTask];
            }
        }

      /* now loop over the new cell layout */
      /* find enclosing rectangle around columns in new plane */

      int first[3], last[3];

      first[0] = out_firstcol / Ndims[perm[1]];
      first[1] = out_firstcol % Ndims[perm[1]];
      first[2] = 0;

      last[0] = (out_firstcol + out_ncol - 1) / Ndims[perm[1]];
      last[1] = (out_firstcol + out_ncol - 1) % Ndims[perm[1]];
      last[2] = Ndims[perm[2]] - 1;

      if(first[1] + out_ncol >= Ndims[perm[1]])
        {
          first[1] = 0;
          last[1]  = Ndims[perm[1]] - 1;
        }

      /* now need to map this back to the old coordinates */

      int xyz_first[3], xyz_last[3];

      for(j = 0; j < 3; j++)
        {
          xyz_first[j] = first[perm_rev[j]];
          xyz_last[j]  = last[perm_rev[j]];
        }

      memset(count_recv, 0, NTask * sizeof(size_t));

      size_t count = 0;

      /* traverse an enclosing box around the new cell layout in the old order */
      for(xyz[0] = xyz_first[0]; xyz[0] <= xyz_last[0]; xyz[0]++)
        for(xyz[1] = xyz_first[1]; xyz[1] <= xyz_last[1]; xyz[1]++)
          for(xyz[2] = xyz_first[2]; xyz[2] <= xyz_last[2]; xyz[2]++)
            {
              /* check that the point is actually part of a column */
              uvw[0] = xyz[perm[0]];
              uvw[1] = xyz[perm[1]];
              uvw[2] = xyz[perm[2]];

              int col = uvw[0] * Ndims[perm[1]] + uvw[1];

              if(col >= out_firstcol && col < out_firstcol + out_ncol)
                {
                  /* determine origin task */
                  int newcol = Ndims[1] * xyz[0] + xyz[1];
                  if(newcol < in_pivotcol)
                    origin = newcol / in_avg;
                  else
                    origin = (newcol - in_pivotcol) / (in_avg - 1) + in_tasklastsection;

                  size_t index = ((size_t)Ndims[perm[2]]) * (col - out_firstcol) + uvw[2];

                  /* move data element from origin task */
                  size_t off    = offset_recv[origin] + count_recv[origin]++;
                  out[index][0] = data[off][0];
                  out[index][1] = data[off][1];

                  count++;
                }
            }

      if(count != nimport)
        {
          int fi = out_firstcol % Ndims[perm[1]];
          int la = (out_firstcol + out_ncol - 1) % Ndims[perm[1]];

          terminate("count=%lld nimport=%lld   ncol=%d fi=%d la=%d first=%d last=%d\n", (long long)count, (long long)nimport, out_ncol,
                    fi, la, first[1], last[1]);
        }
    }
}

/*! \brief Transposes column-based FFT data.
 *
 *  \param[in] data Data to be transposed.
 *  \param[in] Ndims Global number of dimensions of data cube.
 *  \param[in] in_firstcol First column.
 *  \param[in] in_ncol Number of columns.
 *  \param[out] out Data output.
 *  \param[in] perm Permutations in dimensions.
 *  \param[out] out_firstcol First column in output data.
 *  \param[out] out_ncol Number of columns in output data.
 *  \param[out] offset_send Offset in array for send operation to MPI tasks.
 *  \param[out] offset_recv Offset in array for receive operation from MPI
 *              tasks.
 *  \param[out] count_send Count how many elements have to be sent to each
 *              MPI task.
 *  \param[out] count_recv Count how many elements have to be received from
 *              each MPI task.
 *  \param[in] just_count_flag Do element counting for communication instead
 *             of data transfer.
 *
 *  \return void
 */
static void my_fft_column_transpose(fft_real *data, int Ndims[3], int in_firstcol, int in_ncol, fft_real *out, int perm[3],
                                    int out_firstcol, int out_ncol, size_t *offset_send, size_t *offset_recv, size_t *count_send,
                                    size_t *count_recv, size_t just_count_flag)
{
  int j, target, origin, ngrp, recvTask, perm_rev[3], xyz[3], uvw[3];
  size_t nimport, nexport;

  /* determine the inverse permutation */
  for(j = 0; j < 3; j++)
    perm_rev[j] = perm[j];

  if(!(perm_rev[perm[0]] == 0 && perm_rev[perm[1]] == 1 && perm_rev[perm[2]] == 2)) /* not yet the inverse */
    {
      for(j = 0; j < 3; j++)
        perm_rev[j] = perm[perm[j]];

      if(!(perm_rev[perm[0]] == 0 && perm_rev[perm[1]] == 1 && perm_rev[perm[2]] == 2))
        terminate("bummer");
    }

  int in_colums          = Ndims[0] * Ndims[1];
  int in_avg             = (in_colums - 1) / NTask + 1;
  int in_exc             = NTask * in_avg - in_colums;
  int in_tasklastsection = NTask - in_exc;
  int in_pivotcol        = in_tasklastsection * in_avg;

  int out_colums          = Ndims[perm[0]] * Ndims[perm[1]];
  int out_avg             = (out_colums - 1) / NTask + 1;
  int out_exc             = NTask * out_avg - out_colums;
  int out_tasklastsection = NTask - out_exc;
  int out_pivotcol        = out_tasklastsection * out_avg;

  size_t i, ncells = ((size_t)in_ncol) * Ndims[2];

  xyz[0] = in_firstcol / Ndims[1];
  xyz[1] = in_firstcol % Ndims[1];
  xyz[2] = 0;

  memset(count_send, 0, NTask * sizeof(size_t));

  /* loop over all cells in input array and determine target processor */
  for(i = 0; i < ncells; i++)
    {
      /* determine target task */
      uvw[0] = xyz[perm[0]];
      uvw[1] = xyz[perm[1]];
      uvw[2] = xyz[perm[2]];

      int newcol = Ndims[perm[1]] * uvw[0] + uvw[1];
      if(newcol < out_pivotcol)
        target = newcol / out_avg;
      else
        target = (newcol - out_pivotcol) / (out_avg - 1) + out_tasklastsection;

      /* move data element to targettask */

      if(just_count_flag)
        count_send[target]++;
      else
        {
          size_t off = offset_send[target] + count_send[target]++;
          out[off]   = data[i];
        }
      xyz[2]++;
      if(xyz[2] == Ndims[2])
        {
          xyz[2] = 0;
          xyz[1]++;
          if(xyz[1] == Ndims[1])
            {
              xyz[1] = 0;
              xyz[0]++;
            }
        }
    }

  if(just_count_flag)
    {
      MPI_Alltoall(count_send, sizeof(size_t), MPI_BYTE, count_recv, sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);

      for(j = 0, nimport = 0, nexport = 0, offset_send[0] = 0, offset_recv[0] = 0; j < NTask; j++)
        {
          nexport += count_send[j];
          nimport += count_recv[j];

          if(j > 0)
            {
              offset_send[j] = offset_send[j - 1] + count_send[j - 1];
              offset_recv[j] = offset_recv[j - 1] + count_recv[j - 1];
            }
        }

      if(nexport != ncells)
        terminate("nexport=%lld != ncells=%lld", (long long)nexport, (long long)ncells);
    }
  else
    {
      nimport = 0;

      /* exchange all the data */
      for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;

          if(recvTask < NTask)
            {
              if(count_send[recvTask] > 0 || count_recv[recvTask] > 0)
                myMPI_Sendrecv(&out[offset_send[recvTask]], count_send[recvTask] * sizeof(fft_real), MPI_BYTE, recvTask, TAG_DENS_A,
                               &data[offset_recv[recvTask]], count_recv[recvTask] * sizeof(fft_real), MPI_BYTE, recvTask, TAG_DENS_A,
                               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              nimport += count_recv[recvTask];
            }
        }

      /* now loop over the new cell layout */
      /* find enclosing rectangle around columns in new plane */

      int first[3], last[3];

      first[0] = out_firstcol / Ndims[perm[1]];
      first[1] = out_firstcol % Ndims[perm[1]];
      first[2] = 0;

      last[0] = (out_firstcol + out_ncol - 1) / Ndims[perm[1]];
      last[1] = (out_firstcol + out_ncol - 1) % Ndims[perm[1]];
      last[2] = Ndims[perm[2]] - 1;

      if(first[1] + out_ncol >= Ndims[perm[1]])
        {
          first[1] = 0;
          last[1]  = Ndims[perm[1]] - 1;
        }

      /* now need to map this back to the old coordinates */

      int xyz_first[3], xyz_last[3];

      for(j = 0; j < 3; j++)
        {
          xyz_first[j] = first[perm_rev[j]];
          xyz_last[j]  = last[perm_rev[j]];
        }

      memset(count_recv, 0, NTask * sizeof(size_t));

      size_t count = 0;

      /* traverse an enclosing box around the new cell layout in the old order */
      for(xyz[0] = xyz_first[0]; xyz[0] <= xyz_last[0]; xyz[0]++)
        for(xyz[1] = xyz_first[1]; xyz[1] <= xyz_last[1]; xyz[1]++)
          for(xyz[2] = xyz_first[2]; xyz[2] <= xyz_last[2]; xyz[2]++)
            {
              /* check that the point is actually part of a column */
              uvw[0] = xyz[perm[0]];
              uvw[1] = xyz[perm[1]];
              uvw[2] = xyz[perm[2]];

              int col = uvw[0] * Ndims[perm[1]] + uvw[1];

              if(col >= out_firstcol && col < out_firstcol + out_ncol)
                {
                  /* determine origin task */
                  int newcol = Ndims[1] * xyz[0] + xyz[1];
                  if(newcol < in_pivotcol)
                    origin = newcol / in_avg;
                  else
                    origin = (newcol - in_pivotcol) / (in_avg - 1) + in_tasklastsection;

                  size_t index = ((size_t)Ndims[perm[2]]) * (col - out_firstcol) + uvw[2];

                  /* move data element from origin task */
                  size_t off = offset_recv[origin] + count_recv[origin]++;
                  out[index] = data[off];

                  count++;
                }
            }

      if(count != nimport)
        {
          int fi = out_firstcol % Ndims[perm[1]];
          int la = (out_firstcol + out_ncol - 1) % Ndims[perm[1]];

          terminate("count=%lld nimport=%lld   ncol=%d fi=%d la=%d first=%d last=%d\n", (long long)count, (long long)nimport, out_ncol,
                    fi, la, first[1], last[1]);
        }
    }
}

/*! \brief Transposes column-based complex FFT data.
 *
 *  \param[in] data Data to be transposed.
 *  \param[in] Ndims Global number of dimensions of data cube.
 *  \param[in] in_firstcol First column.
 *  \param[in] in_ncol Number of columns.
 *  \param[out] out Data output.
 *  \param[in] perm Permutations in dimensions.
 *  \param[out] out_firstcol First column in output data.
 *  \param[out] out_ncol Number of columns in output data.
 *  \param[out] offset_send Offset in array for send operation to MPI tasks.
 *  \param[out] offset_recv Offset in array for receive operation from MPI
 *              tasks.
 *  \param[out] count_send Count how many elements have to be sent to each
 *              MPI task.
 *  \param[out] count_recv Count how many elements have to be received from
 *              each MPI task.
 *  \param[in] just_count_flag Do element counting for communication instead
 *             of data transfer.
 *
 *  \return void
 */
static void my_fft_column_transpose_c(fft_complex *data, int Ndims[3], int in_firstcol, int in_ncol, fft_complex *out, int perm[3],
                                      int out_firstcol, int out_ncol, size_t *offset_send, size_t *offset_recv, size_t *count_send,
                                      size_t *count_recv, size_t just_count_flag)
{
  int j, target, origin, ngrp, recvTask, perm_rev[3], xyz[3], uvw[3];
  size_t nimport, nexport;

  /* determine the inverse permutation */
  for(j = 0; j < 3; j++)
    perm_rev[j] = perm[j];

  if(!(perm_rev[perm[0]] == 0 && perm_rev[perm[1]] == 1 && perm_rev[perm[2]] == 2)) /* not yet the inverse */
    {
      for(j = 0; j < 3; j++)
        perm_rev[j] = perm[perm[j]];

      if(!(perm_rev[perm[0]] == 0 && perm_rev[perm[1]] == 1 && perm_rev[perm[2]] == 2))
        terminate("bummer");
    }

  int in_colums          = Ndims[0] * Ndims[1];
  int in_avg             = (in_colums - 1) / NTask + 1;
  int in_exc             = NTask * in_avg - in_colums;
  int in_tasklastsection = NTask - in_exc;
  int in_pivotcol        = in_tasklastsection * in_avg;

  int out_colums          = Ndims[perm[0]] * Ndims[perm[1]];
  int out_avg             = (out_colums - 1) / NTask + 1;
  int out_exc             = NTask * out_avg - out_colums;
  int out_tasklastsection = NTask - out_exc;
  int out_pivotcol        = out_tasklastsection * out_avg;

  size_t i, ncells = ((size_t)in_ncol) * Ndims[2];

  xyz[0] = in_firstcol / Ndims[1];
  xyz[1] = in_firstcol % Ndims[1];
  xyz[2] = 0;

  memset(count_send, 0, NTask * sizeof(size_t));

  /* loop over all cells in input array and determine target processor */
  for(i = 0; i < ncells; i++)
    {
      /* determine target task */
      uvw[0] = xyz[perm[0]];
      uvw[1] = xyz[perm[1]];
      uvw[2] = xyz[perm[2]];

      int newcol = Ndims[perm[1]] * uvw[0] + uvw[1];
      if(newcol < out_pivotcol)
        target = newcol / out_avg;
      else
        target = (newcol - out_pivotcol) / (out_avg - 1) + out_tasklastsection;

      /* move data element to targettask */

      if(just_count_flag)
        count_send[target]++;
      else
        {
          size_t off  = offset_send[target] + count_send[target]++;
          out[off][0] = data[i][0];
          out[off][1] = data[i][1];
        }
      xyz[2]++;
      if(xyz[2] == Ndims[2])
        {
          xyz[2] = 0;
          xyz[1]++;
          if(xyz[1] == Ndims[1])
            {
              xyz[1] = 0;
              xyz[0]++;
            }
        }
    }

  if(just_count_flag)
    {
      MPI_Alltoall(count_send, sizeof(size_t), MPI_BYTE, count_recv, sizeof(size_t), MPI_BYTE, MPI_COMM_WORLD);

      for(j = 0, nimport = 0, nexport = 0, offset_send[0] = 0, offset_recv[0] = 0; j < NTask; j++)
        {
          nexport += count_send[j];
          nimport += count_recv[j];

          if(j > 0)
            {
              offset_send[j] = offset_send[j - 1] + count_send[j - 1];
              offset_recv[j] = offset_recv[j - 1] + count_recv[j - 1];
            }
        }

      if(nexport != ncells)
        terminate("nexport=%lld != ncells=%lld", (long long)nexport, (long long)ncells);
    }
  else
    {
      nimport = 0;

      /* exchange all the data */
      for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;

          if(recvTask < NTask)
            {
              if(count_send[recvTask] > 0 || count_recv[recvTask] > 0)
                myMPI_Sendrecv(&out[offset_send[recvTask]], count_send[recvTask] * sizeof(fft_complex), MPI_BYTE, recvTask, TAG_DENS_A,
                               &data[offset_recv[recvTask]], count_recv[recvTask] * sizeof(fft_complex), MPI_BYTE, recvTask,
                               TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              nimport += count_recv[recvTask];
            }
        }

      /* now loop over the new cell layout */
      /* find enclosing rectangle around columns in new plane */

      int first[3], last[3];

      first[0] = out_firstcol / Ndims[perm[1]];
      first[1] = out_firstcol % Ndims[perm[1]];
      first[2] = 0;

      last[0] = (out_firstcol + out_ncol - 1) / Ndims[perm[1]];
      last[1] = (out_firstcol + out_ncol - 1) % Ndims[perm[1]];
      last[2] = Ndims[perm[2]] - 1;

      if(first[1] + out_ncol >= Ndims[perm[1]])
        {
          first[1] = 0;
          last[1]  = Ndims[perm[1]] - 1;
        }

      /* now need to map this back to the old coordinates */

      int xyz_first[3], xyz_last[3];

      for(j = 0; j < 3; j++)
        {
          xyz_first[j] = first[perm_rev[j]];
          xyz_last[j]  = last[perm_rev[j]];
        }

      memset(count_recv, 0, NTask * sizeof(size_t));

      size_t count = 0;

      /* traverse an enclosing box around the new cell layout in the old order */
      for(xyz[0] = xyz_first[0]; xyz[0] <= xyz_last[0]; xyz[0]++)
        for(xyz[1] = xyz_first[1]; xyz[1] <= xyz_last[1]; xyz[1]++)
          for(xyz[2] = xyz_first[2]; xyz[2] <= xyz_last[2]; xyz[2]++)
            {
              /* check that the point is actually part of a column */
              uvw[0] = xyz[perm[0]];
              uvw[1] = xyz[perm[1]];
              uvw[2] = xyz[perm[2]];

              int col = uvw[0] * Ndims[perm[1]] + uvw[1];

              if(col >= out_firstcol && col < out_firstcol + out_ncol)
                {
                  /* determine origin task */
                  int newcol = Ndims[1] * xyz[0] + xyz[1];
                  if(newcol < in_pivotcol)
                    origin = newcol / in_avg;
                  else
                    origin = (newcol - in_pivotcol) / (in_avg - 1) + in_tasklastsection;

                  size_t index = ((size_t)Ndims[perm[2]]) * (col - out_firstcol) + uvw[2];

                  /* move data element from origin task */
                  size_t off    = offset_recv[origin] + count_recv[origin]++;
                  out[index][0] = data[off][0];
                  out[index][1] = data[off][1];

                  count++;
                }
            }

      if(count != nimport)
        {
          int fi = out_firstcol % Ndims[perm[1]];
          int la = (out_firstcol + out_ncol - 1) % Ndims[perm[1]];

          terminate("count=%lld nimport=%lld   ncol=%d fi=%d la=%d first=%d last=%d\n", (long long)count, (long long)nimport, out_ncol,
                    fi, la, first[1], last[1]);
        }
    }
}

#endif /* #ifndef FFT_COLUMN_BASED #else */

#endif /* #if defined(PMGRID) */
