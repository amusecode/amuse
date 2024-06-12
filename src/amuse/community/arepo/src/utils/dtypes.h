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
 * \file        src/utils/dtypes.h
 * \date        05/2018
 * \brief       Definition of intrinsic datatypes.
 * \details
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 28.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifndef DTYPES_H
#define DTYPES_H

#ifndef FFTW
#define CONCAT(prefix, name) prefix##name
#ifdef DOUBLEPRECISION_FFTW
#define FFTW(x) CONCAT(fftw_, x)
#else /* #ifdef DOUBLEPRECISION_FFTW */
#define FFTW(x) CONCAT(fftwf_, x)
#endif /* #ifdef DOUBLEPRECISION_FFTW #else */
#endif /* #ifndef FFTW */

#ifndef LONGIDS
typedef unsigned int MyIDType;
#define MPI_MYIDTYPE MPI_UNSIGNED
#else /* #ifndef LONGIDS */
typedef unsigned long long MyIDType;
#define MPI_MYIDTYPE MPI_UNSIGNED_LONG_LONG
#endif /* #ifndef LONGIDS #else */

#ifndef DOUBLEPRECISION /* default is single-precision */
typedef float MySingle;
typedef float MyFloat;
typedef float MyDouble;
#define MPI_MYFLOAT MPI_FLOAT
#define MPI_MYDOUBLE MPI_FLOAT
#else                     /* #ifndef DOUBLEPRECISION */
#if(DOUBLEPRECISION == 2) /* mixed precision */
typedef float MySingle;
typedef float MyFloat;
typedef double MyDouble;
#define MPI_MYFLOAT MPI_FLOAT
#define MPI_MYDOUBLE MPI_DOUBLE
#else                     /* #if (DOUBLEPRECISION == 2) */
#if(DOUBLEPRECISION == 3) /* mixed precision, fewer single precision variables */
typedef float MySingle;
typedef double MyFloat;
typedef double MyDouble;
#define MPI_MYFLOAT MPI_FLOAT
#define MPI_MYDOUBLE MPI_DOUBLE
#else /* #if (DOUBLEPRECISION == 3) */
/* everything double-precision */
typedef double MySingle;
typedef double MyFloat;
typedef double MyDouble;
#define MPI_MYFLOAT MPI_DOUBLE
#define MPI_MYDOUBLE MPI_DOUBLE
#endif /* #if (DOUBLEPRECISION == 3) #else */
#endif /* #if (DOUBLEPRECISION == 2) #else */
#endif /* #ifndef DOUBLEPRECISION #else */

#ifdef OUTPUT_IN_DOUBLEPRECISION
typedef double MyOutputFloat;
#else  /* #ifdef OUTPUT_IN_DOUBLEPRECISION */
typedef float MyOutputFloat;
#endif /* #ifdef OUTPUT_IN_DOUBLEPRECISION #else */

#ifdef INPUT_IN_DOUBLEPRECISION
typedef double MyInputFloat;
#else  /* #ifdef INPUT_IN_DOUBLEPRECISION */
typedef float MyInputFloat;
#endif /* #ifdef INPUT_IN_DOUBLEPRECISION #else */

#ifndef NGB_TREE_DOUBLEPRECISION
typedef float MyNgbTreeFloat;
#define MAX_NGBRANGE_NUMBER MAX_FLOAT_NUMBER
#else /* #ifndef NGB_TREE_DOUBLEPRECISION */
typedef double MyNgbTreeFloat;
#define MAX_NGBRANGE_NUMBER MAX_DOUBLE_NUMBER
#endif /* #ifndef NGB_TREE_DOUBLEPRECISION #else */

#if defined(PMGRID)
#include <fftw3.h>

#ifdef DOUBLEPRECISION_FFTW
typedef double fft_real;
typedef fftw_complex fft_complex;
#else  /* #ifdef DOUBLEPRECISION_FFTW */
typedef float fft_real;
typedef fftwf_complex fft_complex;
#endif /* #ifdef DOUBLEPRECISION_FFTW #else */
typedef ptrdiff_t fft_ptrdiff_t;

typedef struct
{
  int NgridX, NgridY, NgridZ;
  int Ngridz, Ngrid2;

  FFTW(plan) forward_plan_zdir;
  FFTW(plan) forward_plan_xdir;
  FFTW(plan) forward_plan_ydir;

  FFTW(plan) backward_plan_zdir;
  FFTW(plan) backward_plan_ydir;
  FFTW(plan) backward_plan_xdir;

#ifndef FFT_COLUMN_BASED

  int *slab_to_task; /*!< Maps a slab index to the task responsible for the slab */
  int *slabs_x_per_task;
  int *first_slab_x_of_task; /*!< Array containing the index of the first slab of each task */
  int *slabs_y_per_task;     /*!< Array containing the number of slabs each task is responsible for */
  int *first_slab_y_of_task; /*!< Array containing the index of the first slab of each task */

  int nslab_x, slabstart_x, nslab_y, slabstart_y;
  int largest_x_slab; /*!< size of the largest slab in x direction */
  int largest_y_slab; /*!< size of the largest slab in y direction */

#else /* #ifndef FFT_COLUMN_BASED */

  size_t max_datasize;
  size_t fftsize;

  int base_firstcol, base_ncol, base_lastcol;
  int transposed_firstcol, transposed_ncol;
  int second_transposed_firstcol, second_transposed_ncol;
  size_t second_transposed_ncells;

  int firstcol_XZ, ncol_XZ;
  int firstcol_YZ, ncol_YZ;

  int pivotcol; /* to go from column number to task */
  int avg;
  int tasklastsection;

  size_t *offsets_send_A;
  size_t *offsets_recv_A;
  size_t *offsets_send_B;
  size_t *offsets_recv_B;
  size_t *offsets_send_C;
  size_t *offsets_recv_C;
  size_t *offsets_send_D;
  size_t *offsets_recv_D;
  size_t *offsets_send_13;
  size_t *offsets_recv_13;
  size_t *offsets_send_23;
  size_t *offsets_recv_23;
  size_t *offsets_send_13back;
  size_t *offsets_recv_13back;
  size_t *offsets_send_23back;
  size_t *offsets_recv_23back;

  size_t *count_send_A;
  size_t *count_recv_A;
  size_t *count_send_B;
  size_t *count_recv_B;
  size_t *count_send_C;
  size_t *count_recv_C;
  size_t *count_send_D;
  size_t *count_recv_D;
  size_t *count_send_13;
  size_t *count_recv_13;
  size_t *count_send_23;
  size_t *count_recv_23;
  size_t *count_send_13back;
  size_t *count_recv_13back;
  size_t *count_send_23back;
  size_t *count_recv_23back;

#endif /* #ifndef FFT_COLUMN_BASED */
} fft_plan;

#endif /* #if defined(PMGRID) */

#endif /* #ifndef DTYPES_H */
