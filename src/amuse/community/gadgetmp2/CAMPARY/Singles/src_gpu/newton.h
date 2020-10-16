/*
 * newton.h
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

#ifndef _newton_h
#define _newton_h

/***************************************************************
 ----------------- DIVISION ALGORITHMS -------------------
 ***************************************************************/
/** Computes the reciprocal of sX floats and returns the result on an sR float ulp-nonoveralpping expansion.
    sX - size of x, sR - size of res.
    Uses a recursive Newton-Raphson iteration method for 1/y, with initial guess q_0; 
    The iteration is q_{n+1} = q_n(2-yq_n).
    Constraints: we assume that sR and sX >= 1. Upward version implemented **/
// Reciprocal with relative error <= 2^{-2^q(p-4)-2} / (1-2^{-p+1})
template<int sX, int sR, int pDone>
__host__ __device__ static void __forceinline__ invNewtonUp_aux(const float *x, float *res){
  // pDone is the precision we have already computed in res
  const int p  = (2*pDone>sR) ? sR : 2*pDone;
  const float two[] = {2.};
  float aux[p];

  const int pp = pDone;
	if(pp>((p>sX)?sX:p)) truncatedMul<pDone,(p>sX)?sX:p,p>(res, x, aux);
	else truncatedMul<(p>sX)?sX:p,pDone,p>(x, res, aux);
  for(int j=0; j<p; j++) aux[j] = -aux[j];

  certifiedAdd<p,1,p>(aux, two, aux);
	if(p>pp) truncatedMul<p,pDone,p>(aux, res, res);
	else truncatedMul<pDone,p,p>(res, aux, res);

  if(p<sR) invNewtonUp_aux<sX,sR,p>(x, res);
}
template<int sX, int sR>
__host__ __device__ static void __forceinline__ invNewtonUp(const float *x, float *res){
  // starts by computing an approximation of the reciprocal, 1/x_0
  res[0] = 1./x[0];
  if(sR > 1) invNewtonUp_aux<sX,sR,1>(x, res);
}

/** Computes the division of sX floats by sY floats and returns the result on an sR float ulp-nonoveralpping expansion.
    sX - size of x, sY - size of y, sR - size of res.
    Firstly computes the reciprocal of y and after that multiplies the result by x **/
// Division
template<int sX, int sY, int sR>
__host__ __device__ static void __forceinline__ divNewton(const float *x, const float *y, float *res){
  float aux[sR];
  invNewtonUp<sY,sR>(y, aux);
  const int ss = sX;
  if(ss > sR) truncatedMul<sX,sR,sR>(x, aux, res);
  else truncatedMul<sR,sX,sR>(aux, x, res);
}

/***************************************************************
 ----------------- DIVISION ALGORITHMS 2.0 -------------------
                 using the fast add and mult
 ***************************************************************/
/** Computes the reciprocal of sX floats and returns the result on an sR float ulp-nonoveralpping expansion.
    sX - size of x, sR - size of res.
    Uses a recursive Newton-Raphson iteration method for 1/y, with initial guess q_0; 
    The iteration is q_{n+1} = q_n(2-yq_n).
    Constraints: we assume that sR and sX >= 1. Upward version implemented **/
// Reciprocal
// quick and dirty method that doesn't guaranty the error bound and that the intermediate expansions are ulp-nonoverlapping
template<int sX, int sR, int pDone>
__host__ __device__ static void __forceinline__ invNewtonUp_fast_aux(const float *x, float *res){
  // pDone is the precision we have already computed in res
  const int p  = (2*pDone>sR) ? sR : 2*pDone;
  const float two[] = {2.};
  float aux[p];

  const int pp = pDone;
	if(pp>((p>sX)?sX:p))	baileyMul_fast<pDone,(p>sX)?sX:p,p>(res, x, aux);
	else baileyMul_fast<(p>sX)?sX:p,pDone,p>(x, res, aux);
  for(int j=0; j<p; j++) aux[j] = -aux[j];

	baileyAdd_fast<p,1,p>(aux, two, aux);

	if(p>pp) baileyMul_fast<p,pDone,p>(aux, res, res);
	else baileyMul_fast<pDone,p,p>(res, aux, res);

  if(p<sR) invNewtonUp_fast_aux<sX,sR,p>(x, res);
}
template<int sX, int sR>
__host__ __device__ static void __forceinline__ invNewtonUp_fast(const float *x, float *res){
  // starts by computing an approximation of the reciprocal, 1/x_0
  res[0] = 1./x[0];
  if(sR > 1) invNewtonUp_fast_aux<sX,sR,1>(x, res);
}

/** Computes the division of sX floats by sY floats and returns the result on an sR float ulp-nonoveralpping expansion.
    sX - size of x, sY - size of y, sR - size of res.
    Firstly computes the reciprocal of y and after that multiplies the result by x **/
// Division
// quick and dirty method that doesn't guaranty the error bound and that the intermediate expansions are ulp-nonoverlapping
template<int sX, int sY, int sR>
__host__ __device__ static void __forceinline__ divNewton_fast(const float *x, const float *y, float *res){
  float aux[sR];
  invNewtonUp_fast<sY,sR>(y, aux);
  const int ss = sX;
  if(ss > sR) baileyMul_fast<sX,sR,sR>(x, aux, res);
	else baileyMul_fast<sR,sX,sR>(aux, x, res);
}

/***************************************************************
 ----------------- SQUARE ROOT ALGORITHMS -------------------
 ***************************************************************/
/** Computes the reciprocal of the square root of sX floats and returns the result on an sR float ulp-nonoveralpping expansion.
    sX - size of x, sR - size of res.
    Uses a recursive Newton iteration method for 1/sqrt(y), with initial guess q_0;
    The 'division-free' iteration is q_{n+1} = q_n(3-x*q_n^2)/2
    Constraints: we assume that sR and sX >= 1 **/
// Reciprocal_SquareRoot with relative error <= 2^{-2^q(p-4)-1} / (1-2^{-p+1})
template<int sX, int sR, int pDone>
__host__ __device__ static void __forceinline__ sqrtInvNewton_aux(const float *x, float *res){
  const int p  = (2*pDone>sR) ? sR : 2*pDone;
  const float three[] = {3.};
  float aux[p];

  const int pp = pDone;
	if(((p>sX)?sX:p)>pp) truncatedMul<(p>sX)?sX:p,pDone,p>(x, res, aux);
	else truncatedMul<pDone,(p>sX)?sX:p,p>(res, x, aux);
  
  if(p>pp) truncatedMul<p,pDone,p>(aux, res, aux);
	else truncatedMul<pDone,p,p>(res, aux, aux);

  for(int j=0; j<p; j++) aux[j] = -aux[j];
	certifiedAdd<p,1,p>(aux, three, aux);
	if(p>pp) truncatedMul<p,pDone,p>(aux, res, res);
	else truncatedMul<pDone,p,p>(res, aux, res);
  for(int j=0; j<p; j++) res[j] *= 0.5;

  if(p<sR) sqrtInvNewton_aux<sX,sR,p>(x, res);
}
template<int sX, int sR>
__host__ __device__ static void __forceinline__ sqrtInvNewton(const float *x, float *res){
  if(x[0]<0.) for(int i=0; i<sR; i++) res[i] = nan("");
  else if(x[0]==0.) for(int i=0; i<sR; i++) res[i] = 0.;
  else{
    res[0] = 1./sqrt(x[0]);
    if(sR>1) sqrtInvNewton_aux<sX,sR,1>(x, res);
  }
}

/** Computes the square root of sX floats and returns the result on an sR float ulp-nonoveralpping expansion.
    sX - size of x, sR - size of res.
    Firstly computes the reciprocal of sqrt(x) and after that multiplies the result by x **/
// SquareRoot
template<int sX, int sR>
__host__ __device__ static void __forceinline__ sqrtNewton(const float *x, float *res){
	float aux[sR];
  sqrtInvNewton<sX,sR>(x, aux);
  const int ss = sX;
  if(ss > sR) truncatedMul<sX,sR,sR>(x, aux, res);
  else truncatedMul<sR,sX,sR>(aux, x, res);
}

/** Computes the square root of sX floats and returns the result on an sR float ulp-nonoveralpping expansion.
    Uses a recursive Newton iteration method for sqrt(y), with initial guess q_0;
    The 'Heron' iteration is q_{n+1} = (q_n+x/q_n)/2
    Constraints: we assume that sR and sX >= 1 **/
// SquareRoot_Heron with relative error <= 3 * 2^{-2^q(p-4)-2} / (1-2^{-p+1})
template<int sX, int sR, int pDone>
__host__ __device__ static void __forceinline__ sqrtHeron_aux(const float *x, float *res){
  const int p  = (2*pDone>sR) ? sR : 2*pDone;
  float aux[p];

  divNewton<(p>sX)?sX:p,pDone,p>(x, res, aux);
  certifiedAdd<(p<pDone)?p:pDone,p,p>(res, aux, res);

  for(int j=0; j<p; j++) res[j] *= 0.5;
  if(p<sR) sqrtHeron_aux<sX,sR,p>(x, res);
}
template<int sX, int sR>
__host__ __device__ static void __forceinline__ sqrtHeron(const float *x, float *res){
  if(x[0]<0.) for(int i=0; i<sR; i++) res[i] = nan("");
  else if(x[0]==0.) for(int i=0; i<sR; i++) res[i] = 0.;
  else{
    res[0] = sqrt(x[0]);
    if(sR>1) sqrtHeron_aux<sX,sR,1>(x, res);
  }
}

/***************************************************************
 ----------------- SQUARE ROOT ALGORITHMS 2.0 ------------------
 ***************************************************************/
/** Computes the reciprocal of the square root of sX floats and returns the result on an sR float ulp-nonoveralpping expansion.
    sX - size of x, sR - size of res.
    Uses a recursive Newton iteration method for 1/sqrt(y), with initial guess q_0;
    The 'division-free' iteration is q_{n+1} = q_n(3-x*q_n^2)/2
    Constraints: we assume that sR and sX >= 1 **/
// Reciprocal_SquareRoot
// quick and dirty method that doesn't guaranty the error bound and that the intermediate expansions are ulp-nonoverlapping
template<int sX, int sR, int pDone>
__host__ __device__ static void __forceinline__ sqrtInvNewton_fast_aux(const float *x, float *res){
  const int p  = (2*pDone>sR) ? sR : 2*pDone;
  const float three[] = {3.};
  float aux[p];

  const int pp = pDone;
	if(((p>sX)?sX:p)>pp) baileyMul_fast<(p>sX)?sX:p,pDone,p>(x, res, aux);
	else baileyMul_fast<pDone,(p>sX)?sX:p,p>(res, x, aux);

	if(p>pp) baileyMul_fast<p,pDone,p>(aux, res, aux);
	else baileyMul_fast<pDone,p,p>(res, aux, aux);
	
  for(int j=0; j<p; j++) aux[j] = -aux[j];

	baileyAdd_fast<p,1,p>(aux, three, aux);

	if(p>pp) baileyMul_fast<p,pDone,p>(aux, res, res);
	else baileyMul_fast<pDone,p,p>(res, aux, res);

  for(int j=0; j<p; j++) res[j] *= 0.5;

  if(p<sR) sqrtInvNewton_fast_aux<sX,sR,p>(x, res);
}
template<int sX, int sR>
__host__ __device__ static void __forceinline__ sqrtInvNewton_fast(const float *x, float *res){
  if(x[0]<0.) for(int i=0; i<sR; i++) res[i] = nan("");
  else if(x[0]==0.) for(int i=0; i<sR; i++) res[i] = 0.;
  else{
    res[0] = 1./sqrt(x[0]);
    if(sR>1) sqrtInvNewton_fast_aux<sX,sR,1>(x, res);
  }
}

/** Computes the square root of sX floats and returns the result on an sR float ulp-nonoveralpping expansion.
    sX - size of x, sR - size of res.
    Firstly computes the reciprocal of sqrt(x) and after that multiplies the result by x **/
// SquareRoot
// quick and dirty method that doesn't guaranty the error bound and that the intermediate expansions are ulp-nonoverlapping
template<int sX, int sR>
__host__ __device__ static void __forceinline__ sqrtNewton_fast(const float *x, float *res){
	float aux[sR];
  sqrtInvNewton_fast<sX,sR>(x, aux);
  const int ss = sX;
  if(ss > sR) baileyMul_fast<sX,sR,sR>(x, aux, res);
  else baileyMul_fast<sR,sX,sR>(aux, x, res);
}

/** Computes the square root of sX floats and returns the result on an sR float ulp-nonoveralpping expansion.
    Uses a recursive Newton iteration method for sqrt(y), with initial guess q_0;
    The 'Heron' iteration is q_{n+1} = (q_n+x/q_n)/2
    Constraints: we assume that sR and sX >= 1 **/
// SquareRoot_Heron
// quick and dirty method that doesn't guaranty the error bound and that the intermediate expansions are ulp-nonoverlapping
template<int sX, int sR, int pDone>
__host__ __device__ static void __forceinline__ sqrtHeron_fast_aux(const float *x, float *res){
  const int p  = (2*pDone>sR) ? sR : 2*pDone;
  float aux[p];

  divNewton_fast<(p>sX)?sX:p,pDone,p>(x, res, aux);
  
  const int pp = ((p<pDone)?p:pDone);
	if(pp>p)	baileyAdd_fast<(p<pDone)?p:pDone,p,p>(res, aux, aux);
	else baileyAdd_fast<p,(p<pDone)?p:pDone,p>(aux, res, aux);
	
  for(int j=0; j<p; j++) res[j] *= 0.5;
  if(p<sR) sqrtHeron_fast_aux<sX,sR,p>(x, res);
}
template<int sX, int sR>
__host__ __device__ static void __forceinline__ sqrtHeron_fast(const float *x, float *res){
  if(x[0]<0.) for(int i=0; i<sR; i++) res[i] = nan("");
  else if(x[0]==0.) for(int i=0; i<sR; i++) res[i] = 0.;
  else{
    res[0] = sqrt(x[0]);
    if(sR>1) sqrtHeron_fast_aux<sX,sR,1>(x, res);
  }
}

#endif
