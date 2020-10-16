/*
 * FP_basicOp.h
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

#ifndef _FPbasicOp_h
#define _FPBasicOp_h

/*****************************************************************************************/
/** Templated function for adding two FP numbers, even if they don't have the same type **/
/*****************************************************************************************/
__host__ __device__ static __forceinline__ float FPadd_rn(const float x, const float y){
  #if defined __CUDA_ARCH__
    return __fadd_rn(x, y);
  #else
    return x + y;
  #endif
}

/**********************************************************************************************/
/** Templated function for multiplying two FP numbers, even if they don't have the same type **/
/**********************************************************************************************/
__host__ __device__ static __forceinline__ float FPmul_rn(const float x, const float y){
  #if defined __CUDA_ARCH__
    return __fmul_rn(x, y);
  #else
    return x * y;
  #endif
}

/**********************************************************************************************/
/** Templated function for dividing two FP numbers, even if they don't have the same type **/
/**********************************************************************************************/
__host__ __device__ static __forceinline__ float FPdiv_rn(const float x, const float y){
  #if defined __CUDA_ARCH__
    return __fdiv_rn(x, y);
  #else
    return x / y;
  #endif
}

/*************************************************************************************************/
/** Templated function for the FMA instruction, even if the FP numbers don't have the same type **/
/*************************************************************************************************/
/*Dekker's algorithm
  @ requires xy == \round_float(\NearestEven,x*y) && \abs(x) <= 0x1.p995 && 
  @          \abs(y) <= 0x1.p995 && \abs(x*y) <=  0x1.p1021;
  @ ensures  ((x*y == 0 || 0x1.p-969 <= \abs(x*y)) ==> x*y == xy+\result); @*/
__host__ static __forceinline__ float fma_d_rn_cpu(const float x, const float y, float xy) {
  float C,px,qx,hx,py,qy,hy,lx,ly,fma;
  /*@ assert C == 0x1p27+1; */
	C = 0x1p13+1;

  px = x*C;
  qx = x-px;
  hx = px+qx;
  lx = x-hx;

  py = y*C;
  qy = y-py;
  hy = py+qy;
  ly = y-hy;

  fma = -x*y+hx*hy;
  fma += hx*ly;
  fma += hy*lx;
  fma += lx*ly;
  return fma;
}

__host__ __device__ static __forceinline__ float FPfma_rn(const float x, const float y, const float z){
  #if defined __CUDA_ARCH__
    return __fmaf_rn(x, y, z);
  #else
		#ifdef FP_FAST_FMA
    //#warning cpu has fma
			return fma(x, y, z);
		#else
    	return fma_d_rn_cpu(x, y, z);
		#endif
  #endif
}

#endif
