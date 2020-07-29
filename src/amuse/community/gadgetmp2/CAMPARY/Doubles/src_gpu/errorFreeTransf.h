/*
 * errorFreeTransf.h
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
 *               Mioara Joldes joldes@laas.fr
 */

#ifndef _errorFreeTransf_h
#define _errorFreeTransf_h

/** CPU_FMA: FP_FAST_FMA --constant in math.h,
    if defined, it indicates that fma(S) generally executes about as fast as,
    or faster than a multiply and an add of double operands. **/ 

#include "FP_basicOp.h"

/* Computes fl(a+b) and err(a+b). Assumes |a| >= |b| */
__host__ __device__ static __forceinline__ double fast_two_sum(const double a, const double b, double &err){
  double s = FPadd_rn(a, b);
  double z = FPadd_rn(s, -a);
  err = FPadd_rn(b, -z);
  return s;
}

/* Computes fl(a-b) and err(a-b). Assumes |a| >= |b| */
__host__ __device__ static __forceinline__ double fast_two_diff(const double a, const double b, double &err){
  double s = FPadd_rn(a, -b);
  double z = FPadd_rn(a, -s);
  err = FPadd_rn(z, -b);
  return s;
}

/* Computes fl(a+b) and err(a+b). */
__host__ __device__ static __forceinline__ double two_sum(const double a, const double b, double &err){
  double s = FPadd_rn(a, b);
  double aa = FPadd_rn(s, -b);
  double bb = FPadd_rn(s, -aa);
  double da = FPadd_rn(a, -aa);
  double db = FPadd_rn(b, -bb);
  err = FPadd_rn(da, db);
  return s;
}

/* Computes fl(a-b) and err(a-b). */
__host__ __device__ static __forceinline__ double two_diff(const double a, const double b, double &err){
  double s = FPadd_rn(a, -b);
  double bb = FPadd_rn(s, -a);
  double aa = FPadd_rn(s, -bb);
  double da = FPadd_rn(a, -aa);
  double db = FPadd_rn(b, bb);
  err = FPadd_rn(da, -db);
  return s;
}

/* Computes fl(a*b) and err(a*b). */
__host__ __device__ static __forceinline__ double two_prod(const double a, const double b, double &err){
	const double p = FPmul_rn(a, b);
 	err = FPfma_rn(a, b, -p);
 	return p;
}

#endif
