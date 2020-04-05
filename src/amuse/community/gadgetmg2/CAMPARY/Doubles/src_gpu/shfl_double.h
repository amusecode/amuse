/*
 * shfl_double.h
 *
 * This file is part of CAMPARY Library
 *
 * Copyright (C) 2014 - 
 *
 * CAMPARY Library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * CAMPARY Library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MultiPrecGPU Library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, 
 * Boston, MA  02110-1301  USA
 */
 
/* Contributors: Valentina Popescu valentina.popescu@ens-lyon.fr
 *               Sylvain Collange sylvain.collange@inria.fr
 */

// CUDA device functions for shuffling doubles
__device__ inline int2 double_as_int2(double v) {
  union { double d; int2 i; } caster;
  caster.d = v;
  return caster.i;
}

__device__ inline double int2_as_double(int2 v) {
  union { double d; int2 i; } caster;
  caster.i = v;
  return caster.d;
}

__device__ inline double shfl(double var, int srclane, int width=warpSize) {
  int2 r, v = double_as_int2(var);
  r.x = __shfl(v.x, srclane, width);
  r.y = __shfl(v.y, srclane, width);
  return int2_as_double(r);
}

__device__ inline double shfl_up(double var, unsigned int delta, int width=warpSize) {
  int2 r, v = double_as_int2(var);
  r.x = __shfl_up(v.x, delta, width);
  r.y = __shfl_up(v.y, delta, width);
  return int2_as_double(r);
}

__device__ inline double shfl_down(double var, unsigned int delta, int width=warpSize) {
  int2 r, v = double_as_int2(var);
  r.x = __shfl_down(v.x, delta, width);
  r.y = __shfl_down(v.y, delta, width);
  return int2_as_double(r);
}

