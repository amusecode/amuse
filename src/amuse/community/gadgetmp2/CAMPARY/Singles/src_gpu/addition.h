/*
 * addition.h
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

#ifndef _addition_h
#define _addition_h

/** Merges the arrays x of size K and y of size L, and puts the result in the array z
    All 3 arrays x, y and z are sorted in decreasing order of magnitude
		K - size of x, L - size of y **/
template <int K, int L>
__host__ __device__ static __forceinline__ void merge(float const *x, float const *y, float *z){
  int i=0, j=0;
  for(int n=0; n<K+L; n++){
    if(i==K || (j<L && fabs(y[j])>fabs(x[i]))){ z[n] = y[j]; j++; }
    else{ z[n] = x[i]; i++; }
  }
}

/** Adds x and y and returns the result in r, as an ulp-nonoverlapping expansion. 
		K - size of x, L - size of y, R - size of r **/
// Addition_accurate with relative error < 4.5 * 2^{-(prec-1)R}
template <int K, int L, int R>
__host__ __device__ static __forceinline__ void certifiedAdd(const float *x, const float *y, float *r){
  if(K == 1) renorm2L_4Add1<L,R>(y, x[0], r);
  else if(L == 1) renorm2L_4Add1<K,R>(x, y[0], r);
  else{
    float f[K+L];
    merge<K,L>( x, y, f );
    renorm2L<K+L,R>( f, r );
  }
}

/** Adds x and y and returns the result in r
	Constraints: K>=L, R<K+L (if R>=K+L the extra terms will be 0) **/
// Addition_quick with fast renormalization
// quick and dirty method that doesn't guaranty the result is ulp-nonoverlapping
template <int K, int L, int R>
__host__ __device__ static __forceinline__ void baileyAdd_fast(float const *x, float const *y, float *z){
  if(L == 1) renorm2L_4Add1<K,R>(x, y[0], z);
  else{
	  int i, n;
	  float e, f[R+1];
	
	  for(i=K; i<R; i++) f[i] = 0.0;

    const int nn = R;
	  if (nn<L) f[nn] = FPadd_rn(x[nn], y[nn]);
	  else if (nn<K) f[nn] = x[nn];
	  else f[nn] = 0.0;

	  //from n=L to K-1 we have to add only the elements from x
	  for (n=(R+1>K)?K-1:R-1; n>=L; n--) f[n] = x[n];

	  //until n=L-1 we need to compute the sum x_n + y_n	
	  for (n=(R+1>L)?L-1:R-1; n>=0; n--){
		  f[n] = two_sum(x[n], y[n], e);
		  for(i=n+1; i<R; i++) f[i] = two_sum(f[i], e, e);
		  f[R] = FPadd_rn(f[R], e);
	  }
	
	  fast_renorm2L<R+1,R>(f, z);
	}
}

/** Adds x and y and returns the result in r, as an ulp-nonoverlapping expansion. 
	Constraints: K>=L, R<K+L (if R>=K+L the extra terms will be 0) **/
// Addition_quick with relative error < 24 * 2^{-(prec-1)R}
template <int K, int L, int R>
__host__ __device__ static __forceinline__ void baileyAdd_renorm(float const *x, float const *y, float *z){
  if(L == 1) renorm2L_4Add1<K,R>(x, y[0], z);
  else{
	  int i, n;
	  float e, f[R+1];

	  for(i=K; i<R; i++) f[i] = 0.0;

    const int nn = R;
	  if (nn<L) f[nn] = FPadd_rn(x[nn], y[nn]);
	  else if (nn<K) f[nn] = x[nn];
	  else f[nn] = 0.0;

	  //from n=L to K-1 we have to add only the elements from x
	  for (n=(R+1>K)?K-1:R-1; n>=L; n--) f[n] = x[n];

	  //until n=L-1 we need to compute the sum x_n + y_n	
	  for (n=(R+1>L)?L-1:R-1; n>=0; n--){
		  f[n] = two_sum(x[n], y[n], e);
		  for(i=n+1; i<R; i++) f[i] = two_sum(f[i], e, e);
		  f[R] = FPadd_rn(f[R], e);
	  }
	
	  renorm_rand2L<R+1,R>(f, z);
	}
}

#endif
