/*
 * multiplication.h
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

#ifndef _multiplication_h
#define _multiplication_h

#ifndef UNROLL_MUL_MAIN
//#define UNROLL_MUL_MAIN
#endif

#ifndef UNROLL_MUL_ERROR
//#define UNROLL_MUL_ERROR
#endif

#define dbl_prec 24
#define binSize 18

//recursive case
template<int L, int R, int from, int to, int step=1>
struct unrollPartialProds{
	__host__ __device__ static void __forceinline__ unrollX_trunc(const int exp_start, const int LBN, const float *x, const int *exp_x, const float *y, const int *exp_y, float *B){
		unrollPartialProds< L, R, 0,         (R-from<L?R-from:L), step > :: unrollY      (exp_start, LBN, x[from], exp_x[from], y, exp_y, B);
		unrollPartialProds< L, R, from+step, to,                  step > :: unrollX_trunc(exp_start, LBN, x,       exp_x,       y, exp_y, B);
	}
	
	__host__ __device__ static void __forceinline__ unrollX_all(const int exp_start, const int LBN, const float *x, const int *exp_x, const float *y, const int *exp_y, float *B){
		unrollPartialProds< L, R, 0,          L, step > :: unrollY    (exp_start, LBN, x[from], exp_x[from], y, exp_y, B);
		unrollPartialProds< L, R, from+step, to, step > :: unrollX_all(exp_start, LBN, x,       exp_x,       y, exp_y, B);
	}

	__host__ __device__ static void __forceinline__ unrollY(const int exp_start, const int LBN, const float x, const int exp_x, const float *y, const int *exp_y, float *B){
		float p, e;
		int l  = exp_start - (exp_x+exp_y[from]);
		int sh = (int)(l/binSize);
		l  = l - sh*binSize;
		if(sh < LBN-1){
		  p = two_prod(x, y[from], e);
		  if(l < 30){ // binSize - 2*(dbl_prec-binSize-1) - 1){
		    B[sh] = fast_two_sum(B[sh], p, p);
		    B[sh+1] = FPadd_rn(B[sh+1], p);

		    B[sh+1] = fast_two_sum(B[sh+1], e, e);
		    B[sh+2] = FPadd_rn(B[sh+2], e);
		  }else if(l < 37){ // binSize - (dbl_prec-binSize-1) - 1){
		    B[sh] = fast_two_sum(B[sh], p, p);
		    B[sh+1] = FPadd_rn(B[sh+1], p);

		    B[sh+1] = fast_two_sum(B[sh+1], e, e);
		    B[sh+2] = fast_two_sum(B[sh+2], e, e);
		    B[sh+3] = FPadd_rn(B[sh+3], e);
		  }else{
		    B[sh] = fast_two_sum(B[sh], p, p);
		    B[sh+1] = fast_two_sum(B[sh+1], p, p);
		    B[sh+2] = FPadd_rn(B[sh+2], p);

		    B[sh+2] = fast_two_sum(B[sh+2], e, e);
		    B[sh+3] = FPadd_rn(B[sh+3], e);
		  }
		}
		unrollPartialProds< L, R, from+step, to, step > :: unrollY(exp_start, LBN, x, exp_x, y, exp_y, B);
	}
	
	__host__ __device__ static void __forceinline__ unroll_simple1(const int exp_start, const int LBN, const float *x, const int *exp_x, const float *y, const int *exp_y, float *B){
		float p;
		int sh = (int)((exp_start - (exp_x[from]+exp_y[R-from])) / binSize);
		if(sh < LBN){
			p = FPmul_rn(x[from], y[R-from]);
			B[sh] = fast_two_sum(B[sh], p, p);
		  B[sh+1] = fast_two_sum(B[sh+1], p, p);
		  B[sh+2] = FPadd_rn(B[sh+2], p);
		}
		unrollPartialProds< L, R, from+step, to, step > :: unroll_simple1(exp_start, LBN, x, exp_x, y, exp_y, B);
	}
	
	__host__ __device__ static void __forceinline__ unroll_simple2(const int exp_start, const int LBN, const float *x, const int *exp_x, const float *y, const int *exp_y, float *B){
		float p;
		int sh = (int)((exp_start - (exp_x[R-from]+exp_y[from])) / binSize);
		if(sh < LBN){
			p = FPmul_rn(x[R-from], y[from]);
			B[sh] = fast_two_sum(B[sh], p, p);
		  B[sh+1] = fast_two_sum(B[sh+1], p, p);
		  B[sh+2] = FPadd_rn(B[sh+2], p);
		}
		unrollPartialProds< L, R, from+step, to, step > :: unroll_simple2(exp_start, LBN, x, exp_x, y, exp_y, B);
	}
};

//terminal case
template<int L, int R, int from, int step>
struct unrollPartialProds<L, R, from, from, step>{
	__host__ __device__ static void __forceinline__ unrollX_trunc(const int exp_start, const int LBN, const float *x, const int *exp_x, const float *y, const int *exp_y, float *B){ }
	__host__ __device__ static void __forceinline__ unrollX_all(const int exp_start, const int LBN, const float *x, const int *exp_x, const float *y, const int *exp_y, float *B){ }
	__host__ __device__ static void __forceinline__ unrollY(const int exp_start, const int LBN, const float x, const int exp_x, const float *y, const int *exp_y, float *B){ }
	__host__ __device__ static void __forceinline__ unroll_simple1(const int exp_start, const int LBN, const float *x, const int *exp_x, const float *y, const int *exp_y, float *B){ }
	__host__ __device__ static void __forceinline__ unroll_simple2(const int exp_start, const int LBN, const float *x, const int *exp_x, const float *y, const int *exp_y, float *B){ }
};

/** Multiplies x and y and returns the result in r, as an ulp-nonoverlapping expansion. 
    K - size of x, L - size of y, R - size of r.
    Computes all partial products based on scalar products and then accumulates 
    them in a special designed structure that has a fixed-point flavour.
    float-precision = 53, bin size = 45; **/
// Multiplication_accurate with relative error <= 2^{-(p-1)R}( 1 + (r+1)2^{-p} )
template <int K, int L, int R>
__host__ __device__ static __forceinline__ void certifiedMul(const float *x, const float *y, float *r){
  int const LBN = R*dbl_prec/binSize + 2;
  float B[LBN+2], lim[LBN+2];
  int i;

  int exp_x[K], exp_y[L];
  for(i=0; i<K; i++) frexp(x[i], &exp_x[i]);
  for(i=0; i<L; i++) frexp(y[i], &exp_y[i]);

  float factor = ldexp(1.0, -binSize); /* 2^(-45) */
  int exp_start = exp_x[0] + exp_y[0];
  lim[0] = ldexp(1.5, exp_start - binSize + dbl_prec-1); B[0] = lim[0];
  for(i=1; i<LBN+2; i++){ lim[i] = FPmul_rn(lim[i-1], factor); B[i] = lim[i]; }

	#ifdef UNROLL_MUL_MAIN
		unrollPartialProds< L, R, 0, K, 1 > :: unrollX_all(exp_start, LBN, x, exp_x, y, exp_y, B);
	#else
		float p, e;
		int j, sh, l;
		for(i=0; i<K; i++){
		  for(j=0; j<L; j++){
		    l  = exp_start - (exp_x[i]+exp_y[j]);
		    sh = (int)(l/binSize);
		    l  = l - sh*binSize;
		    if(sh < LBN-1){
		      p = two_prod(x[i], y[j], e);
		      if(l < 30){ // binSize - 2*(dbl_prec-binSize-1) - 1){
		        B[sh] = fast_two_sum(B[sh], p, p);
		        B[sh+1] = FPadd_rn(B[sh+1], p);

		        B[sh+1] = fast_two_sum(B[sh+1], e, e);
		        B[sh+2] = FPadd_rn(B[sh+2], e);
		      }else if(l < 37){ // binSize - (dbl_prec-binSize-1) - 1){
		        B[sh] = fast_two_sum(B[sh], p, p);
		        B[sh+1] = FPadd_rn(B[sh+1], p);

		        B[sh+1] = fast_two_sum(B[sh+1], e, e);
		        B[sh+2] = fast_two_sum(B[sh+2], e, e);
		        B[sh+3] = FPadd_rn(B[sh+3], e);
		      }else{
		        B[sh] = fast_two_sum(B[sh], p, p);
		        B[sh+1] = fast_two_sum(B[sh+1], p, p);
		        B[sh+2] = FPadd_rn(B[sh+2], p);

		        B[sh+2] = fast_two_sum(B[sh+2], e, e);
		        B[sh+3] = FPadd_rn(B[sh+3], e);
		}}}}
	#endif

  /* unbias the B's */
  for (i=0; i<LBN; i++) B[i] = FPadd_rn(B[i], -lim[i]);
  fast_VecSumErrBranch<LBN,R>(B, r);
}

/** Multiplies x and y and returns the result in r, as an ulp-nonoverlapping expansion. 
    K - size of x, L - size of y, R - size of z. Constraints: K>=L.
    The algorithm computes the partial products in a paper-and-pencil fashion and then 
    accumulates them in a special designed structure that has a fixed-point flavour.
		float-precision = 53, bin size = 45;
    For operations using float-float, triple-float and quad-float we provide specialized 
    versions that use a generalization of the QD library's multiplication algorithm. **/
// Multiplication_accurate with relative error <= 2^{-(p-1)R}( 1 + (r+1)2^{-p} + 2^{-p+2}(K+L-R-2) )
template <int K, int L, int R>
__host__ __device__ static __forceinline__ void truncatedMul(const float *x, const float *y, float *r){
	int const LBN = R*dbl_prec/binSize + 2;
  float B[LBN+2], lim[LBN+2];
  int i;

	int exp_x[(R+1<K)?R+1:K], exp_y[(R+1<L)?R+1:L];
	for(i=0; i<((R+1<K)?R+1:K); i++) frexp(x[i], &exp_x[i]);
	for(i=0; i<((R+1<L)?R+1:L); i++) frexp(y[i], &exp_y[i]);

	float factor = ldexp(1.0, -binSize); // 2^(-45)
	int exp_start = exp_x[0] + exp_y[0];
	lim[0] = ldexp(1.5, exp_start - binSize + dbl_prec-1); B[0] = lim[0];
	for(i=1; i<LBN+2; i++){ lim[i] = FPmul_rn(lim[i-1], factor); B[i] = lim[i]; }

	#ifdef UNROLL_MUL_MAIN
	  unrollPartialProds< L, R, 0, (R<K?R:K), 1 > :: unrollX_trunc(exp_start, LBN, x, exp_x, y, exp_y, B);
	#else
		int j, l, sh;
		float p, e;
		for(i=0; i<(R<K?R:K); i++){
			for(j=0; j<(R-i<L?R-i:L); j++){
				l  = exp_start - (exp_x[i]+exp_y[j]);
				sh = (int)(l/binSize);
		  	l  = l - sh*binSize;
			  if(sh < LBN-1){
					p = two_prod(x[i], y[j], e);
			    if(l < 30){ // binSize - 2*(dbl_prec-binSize-1) - 1){
						B[sh] = fast_two_sum(B[sh], p, p);
						B[sh+1] = FPadd_rn(B[sh+1], p);

						B[sh+1] = fast_two_sum(B[sh+1], e, e);
						B[sh+2] = FPadd_rn(B[sh+2], e);
					}else if(l < 37){ // binSize - (dbl_prec-binSize-1) - 1){
						B[sh] = fast_two_sum(B[sh], p, p);
		        B[sh+1] = FPadd_rn(B[sh+1], p);

						B[sh+1] = fast_two_sum(B[sh+1], e, e);
						B[sh+2] = fast_two_sum(B[sh+2], e, e);
		        B[sh+3] = FPadd_rn(B[sh+3], e);
					}else{
						B[sh] = fast_two_sum(B[sh], p, p);
						B[sh+1] = fast_two_sum(B[sh+1], p, p);
		        B[sh+2] = FPadd_rn(B[sh+2], p);

						B[sh+2] = fast_two_sum(B[sh+2], e, e);
		        B[sh+3] = FPadd_rn(B[sh+3], e);
		}}}}
	#endif
	
	//computation of the error correction terms; using just simple multiplication
  if (R < L){
  	#ifdef UNROLL_MUL_ERROR
  		unrollPartialProds< L, R, 0, R+1, 1 > :: unroll_simple1(exp_start, LBN, x, exp_x, y, exp_y, B);
  	#else
  		float p;
  		int sh;
			for(i=0; i<=R; i++){
		  	sh = (int)((exp_start - (exp_x[i]+exp_y[R-i])) / binSize);
		    if(sh < LBN){
					p = FPmul_rn(x[i], y[R-i]);
					B[sh] = fast_two_sum(B[sh], p, p);
		      B[sh+1] = fast_two_sum(B[sh+1], p, p);
		      B[sh+2] = FPadd_rn(B[sh+2], p);
		  }}
    #endif
	}else if(R < K){
  	#ifdef UNROLL_MUL_ERROR
  		unrollPartialProds< L, R, 0, L, 1 > :: unroll_simple2(exp_start, LBN, x, exp_x, y, exp_y, B);
  	#else
  		float p;
  		int sh;
		  for(i=0; i<L; i++){
				sh = (int)((exp_start - (exp_x[R-i]+exp_y[i])) / binSize);
		    if(sh < LBN){
		      p = FPmul_rn(x[R-i], y[i]);
					B[sh] = fast_two_sum(B[sh], p, p);
		      B[sh+1] = fast_two_sum(B[sh+1], p, p);
		      B[sh+2] = FPadd_rn(B[sh+2], p);
		  }}
		#endif
  }else if(R < K+L-1){
  	#ifdef UNROLL_MUL_ERROR
			unrollPartialProds< L, R, R-L+1, K, 1 > :: unroll_simple1(exp_start, LBN, x, exp_x, y, exp_y, B);
  	#else
  		float p;
  		int sh;
		  for(i=R-L+1; i<K; i++){
				sh = (int)((exp_start - (exp_x[i]+exp_y[R-i])) / binSize);
		    if(sh < LBN){
		      p = FPmul_rn(x[i], y[R-i]);
					B[sh] = fast_two_sum(B[sh], p, p);
		      B[sh+1] = fast_two_sum(B[sh+1], p, p);
		      B[sh+2] = FPadd_rn(B[sh+2], p);
			}}
		#endif
	}

	/* unbias the B's */
	for (i=0; i<LBN; i++) B[i] = FPadd_rn(B[i], -lim[i]);
  fast_VecSumErrBranch<LBN,R>(B, r);
}

/** Multiplies x and y and returns the result in r
    K - size of x, L - size of y, R - size of z. Constraints: K>=L.
    Uses a generalization of Bailey's multiplication algorithm that simulates the paper-and-pencil method
		Constraints on the algorithm: K>=L; R<=2*K*L (if R>=2*K*L the extra terms will be 0) **/
// Multiplication_quick with fast renormalization
// quick and dirty method that doesn't guaranty the result is ulp-nonoverlapping
template <int K, int L, int R>
__host__ __device__ static __forceinline__ void baileyMul_fast(float const *x, float const *y, float *z){
  float p, e, f[R+1];
  int i, j, n = R;
  for(i=K+L-1; i<R; i++) f[i] = 0.0;

  //computes the last term of the result f[R]
  //computes the necessary products (if any) using just simple multiplication
  const int nn = R;
  if (nn<L){
    f[nn] = FPmul_rn(x[0], y[nn]);
    for(i=1; i<=nn; i++){
      #if defined __CUDA_ARCH__ || FP_FAST_FMA
        f[nn] = FPfma_rn(x[i], y[nn-i], f[nn]);
      #else
        f[nn] = FPadd_rn(f[nn], FPmul_rn(x[i], y[nn-i]));
      #endif
    }
  }else if (nn<K){
    f[nn] = FPmul_rn(x[nn], y[0]);
    for(j=1; j<L; j++){
      #if defined __CUDA_ARCH__ || FP_FAST_FMA
        f[nn] = FPfma_rn(x[nn-j], y[j], f[nn]);
      #else
        f[nn] = FPadd_rn(f[nn], FPmul_rn(x[nn-j], y[j]));
      #endif
    }
  }else if (nn<K+L-1){
    f[nn] = FPmul_rn(x[nn-L+1], y[L-1]);
    for(i=nn-L+2; i<K; i++){
      #if defined __CUDA_ARCH__ || FP_FAST_FMA
        f[nn] = FPfma_rn(x[i], y[nn-i], f[nn]);
      #else
        f[nn] = FPadd_rn(f[nn], FPmul_rn(x[i], y[nn-i]));
      #endif
    }
  }else
    f[nn] = 0.0;
  
  // computes the last R-K elements of the result
  // we will have (K+L-1 - n) products to compute and sum
  for (n=(R+1>K+L-1)?K+L-2:R-1; n>=K; n--){
    f[n] = two_prod(x[n-L+1], y[L-1], e);
    for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
    f[R] = FPadd_rn(f[R], e);

    for(i=n-L+2; i<K; i++){
      p = two_prod(x[i],y[n-i],e);
      for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
      f[R] = FPadd_rn(f[R], e);

      f[n] = two_sum(f[n], p, e);
      for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
      f[R] = FPadd_rn(f[R], e);
    }
  }
  // computes the elements with the same number of products inside
  // we will have L products to compute and sum
  for (n=(R+1>K)?K-1:R-1; n>=L; n--){
    f[n] = two_prod(x[n], y[0], e);
    for(i=n+1; i<R; i++) f[i] = two_sum(f[i], e, e);
    f[R] = FPadd_rn(f[R], e);

    for(j=1; j<L; j++){
      p = two_prod(x[n-j], y[j], e);
      for(i=n+1; i<R; i++) f[i] = two_sum(f[i], e, e);
      f[R] = FPadd_rn(f[R], e);

      f[n] = two_sum(f[n], p, e);
      for(i=n+1; i<R; i++) f[i] = two_sum(f[i], e, e);
      f[R] = FPadd_rn(f[R], e);
    }
  }
  // computes the first L floats of the result
  // we will have (n+1) prducts to compute and sum
  for (n=(R+1>L)?L-1:R-1; n>=0; n--){
    f[n] = two_prod(x[0], y[n], e);
    for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
    f[R] = FPadd_rn(f[R], e);
    for(i=1; i<=n; i++){
      p = two_prod(x[i], y[n-i], e);
      for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
      f[R] = FPadd_rn(f[R], e);

      f[n] = two_sum(f[n], p, e);
      for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
      f[R] = FPadd_rn(f[R], e);
    }
  }
  fast_renorm2L<R+1,R>(f, z);
}

/** Multiplies x and y and returns the result in r
    K - size of x, L - size of y, R - size of z. Constraints: K>=L.
    Uses a generalization of Bailey's multiplication algorithm that simulates the paper-and-pencil method
		Constraints on the algorithm: K>=L; R<=2*K*L (if R>=2*K*L the extra terms will be 0) **/
// Multiplication_quick with relative error <= 2^{-(p-1)(R+1)} (128/127 (K+L)- 129/254 R - 385/254 + 2^{p-1} + 2^{-p-r}(R^2+R)((R+1)!)^2 )
template <int K, int L, int R>
__host__ __device__ static __forceinline__ void baileyMul_renorm(float const *x, float const *y, float *z){
  float p, e, f[R+1];
  int i, j, n = R;
  for(i=K+L-1; i<R; i++) f[i] = 0.0;

  //computes the last term of the result f[R]
  //computes the necessary products (if any) using just simple multiplication
  const int nn = R;
  if (nn<L){
    f[nn] = FPmul_rn(x[0], y[nn]);
    for(i=1; i<=nn; i++){
      #if defined __CUDA_ARCH__ || FP_FAST_FMA
        f[nn] = FPfma_rn(x[i], y[nn-i], f[nn]);
      #else
        f[nn] = FPadd_rn(f[nn], FPmul_rn(x[i], y[nn-i]));
      #endif
    }
  }else if (nn<K){
    f[nn] = FPmul_rn(x[nn], y[0]);
    for(j=1; j<L; j++){
      #if defined __CUDA_ARCH__ || FP_FAST_FMA
        f[nn] = FPfma_rn(x[nn-j], y[j], f[nn]);
      #else
        f[nn] = FPadd_rn(f[nn], FPmul_rn(x[nn-j], y[j]));
      #endif
    }
  }else if (nn<K+L-1){
    f[nn] = FPmul_rn(x[nn-L+1], y[L-1]);
    for(i=nn-L+2; i<K; i++){
      #if defined __CUDA_ARCH__ || FP_FAST_FMA
        f[nn] = FPfma_rn(x[i], y[nn-i], f[nn]);
      #else
        f[nn] = FPadd_rn(f[nn], FPmul_rn(x[i], y[nn-i]));
      #endif
    }
  }else
    f[nn] = 0.0;

  // computes the last R-K elements of the result
  // we will have (K+L-1 - n) products to compute and sum
  for (n=(R+1>K+L-1)?K+L-2:R-1; n>=K; n--){
    f[n] = two_prod(x[n-L+1], y[L-1], e);
    for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
    f[R] = FPadd_rn(f[R], e);

    for(i=n-L+2; i<K; i++){
      p = two_prod(x[i],y[n-i],e);
      for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
      f[R] = FPadd_rn(f[R], e);

      f[n] = two_sum(f[n], p, e);
      for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
      f[R] = FPadd_rn(f[R], e);
    }
  }
  // computes the elements with the same number of products inside
  // we will have L products to compute and sum
  for (n=(R+1>K)?K-1:R-1; n>=L; n--){
    f[n] = two_prod(x[n], y[0], e);
    for(i=n+1; i<R; i++) f[i] = two_sum(f[i], e, e);
    f[R] = FPadd_rn(f[R], e);

    for(j=1; j<L; j++){
      p = two_prod(x[n-j], y[j], e);
      for(i=n+1; i<R; i++) f[i] = two_sum(f[i], e, e);
      f[R] = FPadd_rn(f[R], e);

      f[n] = two_sum(f[n], p, e);
      for(i=n+1; i<R; i++) f[i] = two_sum(f[i], e, e);
      f[R] = FPadd_rn(f[R], e);
    }
  }
  // computes the first L floats of the result
  // we will have (n+1) prducts to compute and sum
  for (n=(R+1>L)?L-1:R-1; n>=0; n--){
    f[n] = two_prod(x[0], y[n], e);
    for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
    f[R] = FPadd_rn(f[R], e);

    for(i=1; i<=n; i++){
      p = two_prod(x[i], y[n-i], e);
      for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
      f[R] = FPadd_rn(f[R], e);

      f[n] = two_sum(f[n], p, e);
      for(j=n+1; j<R; j++) f[j] = two_sum(f[j], e, e);
      f[R] = FPadd_rn(f[R], e);
    }
  }
  renorm_rand2L<R+1,R>(f, z);
}

#endif
