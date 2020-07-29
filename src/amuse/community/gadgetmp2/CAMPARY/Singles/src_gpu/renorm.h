/*
 * renorm.h
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

#ifndef _renorm_h
#define _renorm_h

#ifndef UNROLL_RENORM
#define UNROLL_RENORM
#endif
/***************************************************************/
/**************** variants of VecSum algorithm *****************/
/***************************************************************/
//recursive case
template<int from, int to, int step=1>
struct unrollVecSum{
  __host__ __device__ static __forceinline__ void fast_VecSum(float *x, float s){
    s = fast_two_sum(x[from], s, x[from+1]);
    unrollVecSum<from+step,to,step>::fast_VecSum(x, s);
  }
  __host__ __device__ static __forceinline__ void fast_VecSum(const float *x, float s, float *r){
    s = fast_two_sum(x[from], s, r[from+1]);
    unrollVecSum<from+step,to,step>::fast_VecSum(x, s, r);
  }
  __host__ __device__ static __forceinline__ void VecSum(float *x, float s){
    s = two_sum(x[from], s, x[from+1]);
    unrollVecSum<from+step,to,step>::VecSum(x, s);
  }
  __host__ __device__ static __forceinline__ void VecSum(const float *x, float s, float *r){
    s = two_sum(x[from], s, r[from+1]);
    unrollVecSum<from+step,to,step>::VecSum(x, s, r);
  }
};
//terminal case
template<int from, int step>
struct unrollVecSum<from,from,step>{
  __host__ __device__ static __forceinline__ void fast_VecSum(float *x, float s){ x[0] = s; }
  __host__ __device__ static __forceinline__ void fast_VecSum(const float *x, float s, float *r){ r[0] = s; }
  __host__ __device__ static __forceinline__ void VecSum(float *x, float s){ x[0] = s; }
  __host__ __device__ static __forceinline__ void VecSum(const float *x, float s, float *r){ r[0] = s; }
};

template<int sX>
__host__ __device__ static __forceinline__ void fast_VecSum(float *x){
	float s = fast_two_sum(x[sX-2], x[sX-1], x[sX-1]);
	#ifdef UNROLL_RENORM
		unrollVecSum<sX-3,-1,-1>::fast_VecSum(x, s);
	#else
		for(int i=sX-3; i>0; i--) s = fast_two_sum(x[i], s, x[i+1]);
		x[0] = fast_two_sum(x[0], s, x[1]);
	#endif
}

template<int sX>
__host__ __device__ static __forceinline__ void fast_VecSum(const float *x, float *r){
	float s = fast_two_sum(x[sX-2], x[sX-1], r[sX-1]);
	#ifdef UNROLL_RENORM
		unrollVecSum<sX-3,-1,-1>::fast_VecSum(x, s, r);
	#else
		for(int i=sX-3; i>0; i--) s = fast_two_sum(x[i], s, r[i+1]);
  	r[0] = fast_two_sum(x[0], s, r[1]);
	#endif
}

template<int sX>
__host__ __device__ static __forceinline__ void VecSum(float *x){
  float s = fast_two_sum(x[sX-2], x[sX-1], x[sX-1]);
  #ifdef UNROLL_RENORM
		unrollVecSum<sX-3,-1,-1>::VecSum(x, s);
	#else
		for(int i=sX-3; i>0; i--) s = two_sum(x[i], s, x[i+1]);
		x[0] = two_sum(x[0], s, x[1]);
	#endif
}

template<int sX>
__host__ __device__ static __forceinline__ void VecSum(const float *x, float *r){
  float s = fast_two_sum(x[sX-2], x[sX-1], r[sX-1]);
	#ifdef UNROLL_RENORM
  	unrollVecSum<sX-3,-1,-1>::VecSum(x, s, r);
	#else
	  for(int i=sX-3; i>0; i--) s = two_sum(x[i], s, r[i+1]);
  	r[0] = two_sum(x[0], s, r[1]);
	#endif
}

/** VecSum for addition of a FP expansion with a FP number **/
template<int sX>
__host__ __device__ static __forceinline__ void VecSum_4Add1(const float *x, const float y, float *r){
  float s = two_sum(x[sX-1], y, r[sX]);
	#ifdef UNROLL_RENORM
  	unrollVecSum<sX-2,-1,-1>::VecSum(x, s, r);
	#else
	  for(int i=sX-2; i>0; i--) s = two_sum(x[i], s, r[i+1]);
  	r[0] = two_sum(x[0], s, r[1]);
	#endif
}

/***************************************************************/
/*********** variants of VecSumErrBranch algorithm *************/
/***************************************************************/
template<int sX, int sR>
__host__ __device__ static __forceinline__ void fast_VecSumErrBranch(float *x){
	int ptr = 0, i = 1;
	float e = x[0];

  while(i<sX && ptr<sR){
    x[ptr] = fast_two_sum(e, x[i], e); i++;
    if(e == 0.) e = x[ptr]; else ptr++;
  }
  if(ptr<sR && e!=0.){ x[ptr] = e; ptr++; }
	for(i=ptr; i<sR; i++) x[i] = 0.;
}

template<int sX, int sR>
__host__ __device__ static __forceinline__ void fast_VecSumErrBranch(const float *x, float *r){
  int ptr = 0, i = 1;
  float e = x[0];

  while(i<sX && ptr<sR){
    r[ptr] = fast_two_sum(e, x[i], e); i++;
    if(e == 0.) e = r[ptr]; else ptr++;
  }
  if(ptr<sR && e!=0.){ r[ptr] = e; ptr++; }
  for(i=ptr; i<sR; i++) r[i] = 0.;
}

/***************************************************************/
/************** variants of VecSumErr algorithm ****************/
/***************************************************************/
template<int sX>
__host__ __device__ static __forceinline__ void fast_VecSumErr(float *x){
	float e;
	x[0] = fast_two_sum(x[0], x[1], e);
	for(int i=2; i<sX-1; i++) x[i-1] = fast_two_sum(e, x[i], e);
	x[sX-2] = fast_two_sum(e, x[sX-1], x[sX-1]);
}

template<int sX>
__host__ __device__ static __forceinline__ void fast_VecSumErr(const float *x, float *r){
  float e;
  r[0] = fast_two_sum(x[0], x[1], e);
  for(int i=2; i<sX-1; i++) r[i-1] = fast_two_sum(e, x[i], e);
  r[sX-2] = fast_two_sum(e, x[sX-1], r[sX-1]);
}

/***************************************************************/
/********** variants of the renormalization algorithm **********/
/***************************************************************/

/** Algorithm for normalizing the array x. The result satisfies |x_i|<=ulp(x_{i-1}); ulp-nonoverlapping expansion.
    sX is the size of array x and sR is the required size for the normalized expansion
    Requirement: sR <= sX **/
// Renormalize
template<int sX, int sR>
__host__ __device__ static __forceinline__ void fast_renorm2L(float *x){
	float e;
	int ptr = 0, i = 2;

  fast_VecSum<sX>(x);

  if(x[1] == 0.) e = x[0]; 
  else { e = x[1]; ptr++;}
  while(ptr<sR && i<sX){
    x[ptr] = fast_two_sum(e, x[i], e); i++;
    if(e == 0.) e = x[ptr]; else ptr++;
  }
  if(ptr<sR && e!=0.){ x[ptr] = e; ptr++; }
  for(i=ptr; i<sX; i++) x[i] = 0.;
}
template<int sX, int sR>
__host__ __device__ static __forceinline__ void fast_renorm2L(const float *x, float *r){
  float e, f[sX];
  int ptr = 0, i = 2;

  fast_VecSum<sX>(x,f);

  if(f[1] == 0.) e = f[0];
  else { r[0] = f[0]; e = f[1]; ptr++; }
  while(ptr<sR && i<sX){
    r[ptr] = fast_two_sum(e, f[i], e); i++;
    if(e == 0.) e = r[ptr]; else ptr++;
  }
  if(ptr<sR && e!=0.){ r[ptr] = e; ptr++; }
  for(i=ptr; i<sR; i++) r[i] = 0.;
}

template<int sX, int sR>
__host__ __device__ static __forceinline__ void renorm2L(float *x){
  float e;
  int ptr = 0, i = 2;

  VecSum<sX>(x);

  if(x[1] == 0.) e = x[0];
  else { e = x[1]; ptr++;}
  while(ptr<sR && i<sX){
    x[ptr] = fast_two_sum(e, x[i], e); i++;
    if(e == 0.) e = x[ptr]; else ptr++;
  }
  if(ptr<sR && e!=0.){ x[ptr] = e; ptr++; }
  for(i=ptr; i<sX; i++) x[i] = 0.;
}
template<int sX, int sR>
__host__ __device__ static __forceinline__ void renorm2L(const float *x, float *r){
  float e, f[sX];
  int ptr = 0, i = 2;

  VecSum<sX>(x,f);

  if(f[1] == 0.) e = f[0];
  else { r[0] = f[0]; e = f[1]; ptr++; }
  while(ptr<sR && i<sX){
    r[ptr] = fast_two_sum(e, f[i], e); i++;
    if(e == 0.) e = r[ptr]; else ptr++;
  }
  if(ptr<sR && e!=0.){ r[ptr] = e; ptr++; }
  for(i=ptr; i<sR; i++) r[i] = 0.;
}

/** Addition of a FP expansion with a FP number **/
// Renormalize - special case for adding an expansion with a FP number
template<int sX, int sR>
__host__ __device__ static __forceinline__ void renorm2L_4Add1(const float *x, const float y, float *r){
  float e, f[sX+1];
  int ptr = 0, i = 2;

  VecSum_4Add1<sX>(x,y,f);

  if(f[1] == 0.) e = f[0];
  else { r[0] = f[0]; e = f[1]; ptr++; }
  while(ptr<sR && i<sX+1){
    r[ptr] = fast_two_sum(e, f[i], e); i++;
    if(e == 0.) e = r[ptr]; else ptr++;
  }
  if(ptr<sR && e!=0.){ r[ptr] = e; ptr++; }
  for(i=ptr; i<sR; i++) r[i] = 0.;
}

/** Algorithm for normalizing the array x. The result satisfies |x_i|<ulp(x_{i-1}); P-nonoverlapping expansion.
    sX is the size of array x and sR is the required size for the normalized expansion
    Requirement: sR <= sX **/
template<int sX, int sR>
__host__ __device__ static __forceinline__ void fast_renorm3L(float *x){
  float e;
	int ptr = 0, i = 2;

	fast_VecSum<sX>(x);

	if(x[1] == 0.) e = x[0];
  else { e = x[1]; ptr++; }
  while(ptr<((sR<sX)?sR+1:sR) && i<sX){
    x[ptr] = fast_two_sum(e, x[i], e); i++;
    if(e == 0.) e = x[ptr]; else ptr++;
  }
  if(i==sX && ptr<((sR<sX)?sR+1:sR) && e!=0.){ x[ptr] = e; ptr++; }
	for(i=ptr; i<((sR<sX)?sR+1:sR); i++) x[i] = 0.;

  for(int j=0; j<((sR<sX)?sR-1:sX-2); j++){
    x[j] = fast_two_sum(x[j], x[j+1], e);
    for(int i=j+1; i<((sR<sX)?sR-1:sX-2); i++) x[i] = fast_two_sum(e, x[i+1], e);
    x[(sR<sX)?sR-1:sX-2] = fast_two_sum(e, x[(sR<sX)?sR:sX-1], x[(sR<sX)?sR:sX-1]);
  }
  for(i=sR; i<sX; i++) x[i] = 0.;
}

template<int sX, int sR>
__host__ __device__ static __forceinline__ void fast_renorm3L(float *x, float *r){
  float e, f[sX];
  int ptr = 0, i = 2;

  fast_VecSum<sX>(x,f);

  if(f[1] == 0.) e = f[0]; else { e = f[1]; ptr++; }
  while(ptr<((sR<sX)?sR+1:sR) && i<sX){
    f[ptr] = fast_two_sum(e, f[i], e); i++;
    if(e == 0.) e = f[ptr]; else ptr++;
  }
  if(i==sX && ptr<((sR<sX)?sR+1:sR) && e!=0.){ f[ptr] = e; ptr++; }
  for(i=ptr; i<((sR<sX)?sR+1:sR); i++) f[i] = 0.;

  for(int j=0; j<((sR<sX)?sR-1:sX-2); j++){
    r[j] = fast_two_sum(f[j], f[j+1], e);
    for(int i=j+1; i<((sR<sX)?sR-1:sX-2); i++) f[i] = fast_two_sum(e, f[i+1], e);
    f[(sR<sX)?sR-1:sX-2] = fast_two_sum(e, f[(sR<sX)?sR:sX-1], f[(sR<sX)?sR:sX-1]);
  }
  r[(sR<sX)?sR-1:sX-2] = f[(sR<sX)?sR-1:sX-2];
  r[(sR<sX)?sR-1:sX-1] = f[(sR<sX)?sR-1:sX-1];

  for(i=sR; i<sX; i++) r[i] = 0.;
}

/** Algorithm for normalizing the array x that contains random numbers. 
    After the first level	the result satisfies |x_i|<uls(x_{i-1}); S-nonoverlapping expansion.
		In the end, the result satisfies |x_i|<=ulp(x_{i-1}); ulp-nonoverlapping expansion.
    sX is the size of array x and sR is the required size for the normalized expansion
    Requirement: sR <= sX **/
// Renormalize_random
template<int sX, int sR>
__host__ __device__ static __forceinline__ void renorm_rand2L(float *x){
	float pr;
	int i, j, ptr = 0;

  x[0] = two_sum(x[0], x[1], x[1]);
  for(i=2; i<sX; i++){
		pr = two_sum(x[i-1], x[i], x[i]);
    for(j=i-2; j>0; j--) pr = two_sum(x[j], pr, x[j+1]);
		x[0] = fast_two_sum(x[0], pr, x[1]);
	}

	i = 2;
  if(x[1] == 0.) pr = x[0];
  else { pr = x[1]; ptr++;}
  while(ptr<sR && i<sX){
    x[ptr] = fast_two_sum(pr, x[i], pr); i++;
    if(pr == 0.) pr = x[ptr]; else ptr++;
  }
  if(ptr<sR && pr!=0.){ x[ptr] = pr; ptr++; }
  for(i=ptr; i<sX; i++) x[i] = 0.;
}

template<int sX, int sR>
__host__ __device__ static __forceinline__ void renorm_rand2L(const float *x, float *r){
  float pr, f[sX];
  int i, j, ptr = 0;

  f[0] = two_sum(x[0], x[1], f[1]);
  for(i=2; i<sX; i++){
    pr = two_sum(f[i-1], x[i], f[i]);
    for(j=i-2; j>0; j--) pr = two_sum(f[j], pr, f[j+1]);
    f[0] = fast_two_sum(f[0], pr, f[1]);
  }

  i = 2;
  if(f[1] == 0.) pr = f[0];
  else { r[0] = f[0]; pr = f[1]; ptr++; }
  while(ptr<sR && i<sX){
    r[ptr] = fast_two_sum(pr, f[i], pr); i++;
    if(pr == 0.) pr = r[ptr]; else ptr++;
  }
  if(ptr<sR && pr!=0.){ r[ptr] = pr; ptr++; }
  for(i=ptr; i<sR; i++) r[i] = 0.;
}

#endif
