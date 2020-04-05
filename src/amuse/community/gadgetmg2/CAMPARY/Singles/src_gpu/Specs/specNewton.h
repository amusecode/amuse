/*
 * specNewton.h
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

#ifndef _specNewton_h
#define _specNewton_h

/****************************************************/
template<> __host__ __device__ void __forceinline__ invNewtonUp<1,1>(const float *x, float *res){ res[0] = FPdiv_rn( 1., x[0] ); }

/****************************************************/
template<> __host__ __device__ void __forceinline__ divNewton<1,1,1>(const float *x, const float *y, float *res){ res[0] = FPdiv_rn( x[0], y[0] ); }
template<> __host__ __device__ void __forceinline__ divNewton<2,1,2>(const float *x, const float *y, float *res){
// DWDivFP2 with relative error <= 3.5u^2
  float t, rh, rl;
  t = FPdiv_rn( x[0], y[0] );
  rh = two_prod( t, y[0], rl );
  rh = FPadd_rn( x[0], -rh );
  rl = FPadd_rn( x[1], -rl );
  rh = FPadd_rn( rh, rl );
  rh = FPdiv_rn( rh, y[0] );
  res[0] = fast_two_sum( t, rh, res[1] );
}
template<> __host__ __device__ void __forceinline__ divNewton<2,2,2>(const float *x, const float *y, float *res){
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
  // DWDivDW3 with relative error <= 9.8u^2
    float t, r[2], d[2];
    t = FPdiv_rn( 1, y[0] );
    r[0] = FPfma_rn( -y[0], t, 1 ); 
    r[1] = -FPmul_rn( y[1], t );
    r[0] = fast_two_sum( r[0], r[1], r[1] );
    certifiedMul<2,1,2>( r, &t, d );
    certifiedAdd<2,1,2>( d, &t, r );
    certifiedMul<2,2,2>( x, r, res );
  #else
  // DWDivDW2 with relative error <= 15u^2 + 56u^3
    float t, r[2];
    t = FPdiv_rn( x[0], y[0] );
    certifiedMul<2,1,2>( y, &t, r );
    r[0] = FPadd_rn( x[0], -r[0] );
    r[1] = FPadd_rn( x[1], -r[1] );
    r[0] = FPadd_rn( r[0], r[1] );
    r[0] = FPdiv_rn( r[0], y[0] );
    res[0] = fast_two_sum( t, r[0], res[1] );
  #endif
}

/****************************************************/
template<> __host__ __device__ void __forceinline__ invNewtonUp_fast<1,1>(const float *x, float *res){ res[0] = 1./x[0]; }

/****************************************************/
template<> __host__ __device__ void __forceinline__ divNewton_fast<1,1,1>(const float *x, const float *y, float *res){ res[0] = x[0] / y[0]; }
template<> __host__ __device__ void __forceinline__ divNewton_fast<2,1,2>(const float *x, const float *y, float *res){
// DWDivFP2 with relative error <= 3.5u^2
  float t, rh, rl;
  t = FPdiv_rn( x[0], y[0] );
  rh = two_prod( t, y[0], rl );
  rh = FPadd_rn( x[0], -rh );
  rl = FPadd_rn( x[1], -rl );
  rh = FPadd_rn( rh, rl );
  rh = FPdiv_rn( rh, y[0] );
  res[0] = fast_two_sum( t, rh, res[1] );
/* DWDivFP1, algo from QD lib, mathematically equivalent, but more FLOPS  
  float th, ph, pl, dh, dp, ds, d, dl, tl;
  th = FPdiv_rn( x[0], y[0] );
  ph = two_prod( th, y[0], pl );
  dh = two_sum( x[0], -ph, dp);
  ds = FPadd_rn( x[1], -pl );
  dl = FPadd_rn( dp, ds );
  d = FPadd_rn( dh, dl );
  tl = FPdiv_rn( d, y[0] );
  res[0] = fast_two_sum( th, tl, res[1] ); */
}
template<> __host__ __device__ void __forceinline__ divNewton_fast<2,2,2>(const float *x, const float *y, float *res){
// DWDivDW2 with relative error <= 15u^2 + 56u^3
  float t, r[2];
  t = FPdiv_rn( x[0], y[0] );
  certifiedMul<2,1,2>( y, &t, r );
  r[0] = FPadd_rn( x[0], -r[0] );
  r[1] = FPadd_rn( x[1], -r[1] );
  r[0] = FPadd_rn( r[0], r[1] );
  r[0] = FPdiv_rn( r[0], y[0] );
  res[0] = fast_two_sum( t, r[0], res[1] );
/* DWDivDW1, algo from QD lib, mathematically equivalent, but more FLOPS    
  float th, tl, r[2], ph, pl, dh, dl, d;
  th = FPdiv_rn( x[0], y[0] );
  certifiedMul<2,1,2>( y, &th, r );
  ph = two_sum( x[0], -r[0], pl );
  dh = FPadd_rn( pl, -r[1] );
  dl = FPadd_rn( dh, x[1] );
  d = FPadd_rn( ph, dl );
  tl = FPdiv_rn( d, y[0] );
  res[0] = fast_two_sum( th, tl, res[1] ); */
}

/****************************************************/
template<> __host__ __device__ void __forceinline__ sqrtInvNewton<1,1>(float const *x, float *res){ res[0] = 1/sqrt(x[0]); }

/****************************************************/
template<> __host__ __device__ void __forceinline__ sqrtNewton<1,1>(float const *x, float *res){ res[0] = sqrt(x[0]); }

/****************************************************/
template<> __host__ __device__ void __forceinline__ sqrtHeron<1,1>(float const *x, float *res){ res[0] = sqrt(x[0]); }

/****************************************************/
template<> __host__ __device__ void __forceinline__ sqrtInvNewton_fast<1,1>(float const *x, float *res){ res[0] = 1/sqrt(x[0]); }

/****************************************************/
template<> __host__ __device__ void __forceinline__ sqrtNewton_fast<1,1>(float const *x, float *res){ res[0] = sqrt(x[0]); }

/****************************************************/
template<> __host__ __device__ void __forceinline__ sqrtHeron_fast<1,1>(float const *x, float *res){ res[0] = sqrt(x[0]); }

#endif
