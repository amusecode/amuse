/*
 * specMultiplication.h
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

#ifndef _specMultiplication_h
#define _specMultiplication_h

/****************************************************/
template<> __host__ __device__ __forceinline__ void certifiedMul<1,1,1>(double const *x, double const *y, double *r){ r[0] = FPmul_rn(x[0], y[0]); }
template<> __host__ __device__ __forceinline__ void certifiedMul<1,1,2>(double const *x, double const *y, double *r){ r[0] = two_prod(x[0], y[0], r[1]); }
template<> __host__ __device__ __forceinline__ void certifiedMul<1,2,2>(double const *x, double const *y, double *r){
// DWTimesFP1 or DWTimesFP3 with relative error <= 2u^2
  double ch, cl1, cl2;
  ch = two_prod( x[0], y[0], cl1 );
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
    cl2 = FPfma_rn( x[0], y[1], cl1 );
  #else
    cl2 = FPmul_rn( x[0], y[1] );
    ch = fast_two_sum( ch, cl2, cl2 );
    cl2 = FPadd_rn( cl2, cl1 );
  #endif
  r[0] = fast_two_sum( ch, cl2, r[1] );
}
template<> 
__host__ __device__ __forceinline__ void certifiedMul<2,1,2>(double const *x, double const *y, double *r){
// DWTimesFP1 or DWTimesFP3 with relative error <= 2u^2
  double ch, cl1, cl2;
  ch = two_prod( x[0], y[0], cl1 );  
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
    cl2 = FPfma_rn( x[1], y[0], cl1 );
  #else
    cl2 = FPmul_rn( x[1], y[0] );
    ch = fast_two_sum( ch, cl2, cl2 );
    cl2 = FPadd_rn( cl2, cl1 );
  #endif
  r[0] = fast_two_sum( ch, cl2, r[1] );
}
template<> 
__host__ __device__ __forceinline__ void certifiedMul<2,2,2>(double const *x, double const *y, double *r){
// DWTimesDW1 with relative error <= 7u^2 or DWTimesDW3 with relative error <= 5u^2
  double ch, cl, tl1, tl2;
  ch = two_prod( x[0], y[0], cl );
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
    double tl0 = FPmul_rn( x[1], y[1] );
    tl1 = FPfma_rn( x[0], y[1], tl0 );
    tl2 = FPfma_rn( x[1], y[0], tl1 );
    cl = FPadd_rn( cl, tl2 );
  #else
    tl1 = FPmul_rn( x[0], y[1] );
    tl2 = FPmul_rn( x[1], y[0] );
    cl = FPadd_rn( cl, FPadd_rn( tl1, tl2 ) );
  #endif
  r[0] = fast_two_sum( ch, cl, r[1] );
}

/****************************************************/
template<> __host__ __device__ __forceinline__ void truncatedMul<1,1,1>(double const *x, double const *y, double *z){ z[0] = FPmul_rn(x[0], y[0]); }
template<> __host__ __device__ __forceinline__ void truncatedMul<1,1,2>(double const *x, double const *y, double *z){ z[0] = two_prod(x[0], y[0], z[1]); }
template<>
__host__ __device__ __forceinline__ void truncatedMul<2,1,2>(double const *x, double const *y, double *z){
// DWTimesFP2 with relative error <= 3u^2 or DWTimesFP3 with relative error <= 2u^2
  double ch, cl1, cl2;
  ch = two_prod( x[0], y[0], cl1 );
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
    cl2 = FPfma_rn( x[1], y[0], cl1 );
  #else
    cl2 = FPmul_rn( x[1], y[0] );
    cl2 = FPadd_rn( cl1, cl2 );
  #endif
  z[0] = fast_two_sum( ch, cl2, z[1] );
}
template<>
__host__ __device__ __forceinline__ void truncatedMul<2,2,2>(double const *x, double const *y, double *z){
// DWTimesDW1 with relative error <= 7u^2 or DWTimesDW2 with relative error <= 6u^2
  double ch, cl, tl1, tl2;
  ch = two_prod( x[0], y[0], cl );
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
    tl1 = FPmul_rn( x[0], y[1] );
    tl2 = FPfma_rn( x[1], y[0], tl1 );
    cl = FPadd_rn( cl, tl2 );
  #else
    tl1 = FPmul_rn( x[0], y[1] );
    tl2 = FPmul_rn( x[1], y[0] );
    cl = FPadd_rn( cl, FPadd_rn( tl1, tl2 ) );
  #endif
  z[0] = fast_two_sum( ch, cl, z[1] );
}
template<>
__host__ __device__ __forceinline__ void truncatedMul<3,1,3>(double const *x, double const *y, double *z){
  double e, f[4];
  f[2] = two_prod( x[2], y[0], f[3] );
  f[1] = two_prod( x[1], y[0], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
  f[0] = two_prod( x[0], y[0], e );
  f[1] = two_sum ( f[1], e,    e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
	renorm_rand2L<4,3>(f, z);
}
template<>
__host__ __device__ __forceinline__ void truncatedMul<3,3,3>(double const *x, double const *y, double *z){
  double p, e, f[4];
  f[3] = FPmul_rn( x[1], y[2] );
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
    f[3] = FPfma_rn( x[2], y[1], f[3] );
  #else
    f[3] = FPadd_rn( f[3], FPmul_rn( x[2], y[1] ) );
  #endif
  f[2] = two_prod( x[0], y[2], e );
  f[3] = FPadd_rn( f[3], e );
     p = two_prod( x[1], y[1], e );
  f[3] = FPadd_rn( f[3], e );
  f[2] = two_sum ( f[2], p,    e );
  f[3] = FPadd_rn( f[3], e );
     p = two_prod( x[2], y[0], e );
  f[3] = FPadd_rn( f[3], e );
  f[2] = two_sum ( f[2], p,    e );
  f[3] = FPadd_rn( f[3], e );
  f[1] = two_prod( x[0], y[1], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
     p = two_prod( x[1], y[0], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
  f[1] = two_sum ( f[1], p,    e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
  f[0] = two_prod( x[0], y[0], e );
  f[1] = two_sum ( f[1], e,    e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
  renorm_rand2L<4,3>(f, z);
}
template<>
__host__ __device__ __forceinline__ void truncatedMul<4,1,4>(double const *x, double const *y, double *z){
  double e, f[5];
  f[3] = two_prod( x[3], y[0], f[4] );
  f[2] = two_prod( x[2], y[0], e );
  f[3] = two_sum ( f[3], e,    e ); 
  f[4] = FPadd_rn( f[4], e );
  f[1] = two_prod( x[1], y[0], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  f[0] = two_prod( x[0], y[0], e );
  f[1] = two_sum ( f[1], e,    e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  renorm_rand2L<5,4>(f, z);
}
template<>
__host__ __device__ __forceinline__ void truncatedMul<4,4,4>(double const *x, double const *y, double *z){
  double p, e, f[5];
  f[4] = FPmul_rn( x[1], y[3] );
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
    f[4] = FPfma_rn( x[2], y[2], f[4] );
    f[4] = FPfma_rn( x[3], y[1], f[4] );
  #else
    f[4] = FPadd_rn( f[4], FPmul_rn( x[2], y[2] ) );
    f[4] = FPadd_rn( f[4], FPmul_rn( x[3], y[1] ) );
  #endif  
  f[3] = two_prod( x[0], y[3], e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[1], y[2], e );
  f[4] = FPadd_rn( f[4], e );
  f[3] = two_sum ( f[3], p,    e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[2], y[1], e );
  f[4] = FPadd_rn( f[4], e );
  f[3] = two_sum ( f[3], p,    e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[3], y[0], e );
  f[4] = FPadd_rn( f[4], e );
  f[3] = two_sum ( f[3], p,    e );
  f[4] = FPadd_rn( f[4], e );
  f[2] = two_prod( x[0], y[2], e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[1], y[1], e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  f[2] = two_sum ( f[2], p,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[2], y[0], e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  f[2] = two_sum ( f[2], p,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  f[1] = two_prod( x[0], y[1], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[1], y[0], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  f[1] = two_sum ( f[1], p,    e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  f[0] = two_prod( x[0], y[0], e );
  f[1] = two_sum ( f[1], e,    e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
	renorm_rand2L<5,4>(f, z);
}
template<>
__host__ __device__ __forceinline__ void truncatedMul<3,2,2>(double const *x, double const *y, double *z){
  double p, e, f[3];
  f[2] = FPmul_rn( x[1], y[1] );
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
    f[2] = FPfma_rn( x[2], y[0], f[2] );
  #else
    f[2] = FPadd_rn( f[2], FPmul_rn( x[2], y[0] ));
  #endif
  f[1] = two_prod( x[0], y[1], e );
  f[2] = FPadd_rn( f[2], e );
     p = two_prod( x[1], y[0], e );
  f[2] = FPadd_rn( f[2], e );
  f[1] = two_sum ( f[1], p,    e );
  f[2] = FPadd_rn( f[2], e );
  f[0] = two_prod( x[0], y[0], e );
  f[1] = two_sum ( f[1], e,    e );
  f[2] = FPadd_rn( f[2], e );
  renorm_rand2L<3,2>(f, z);
}
template<>
__host__ __device__ __forceinline__ void truncatedMul<3,2,3>(double const *x, double const *y, double *z){
  double p, e, f[4];
  f[3] = FPmul_rn( x[2], y[1] );
  
  f[2] = two_prod( x[2], y[0], e );
  f[3] = FPadd_rn( f[3], e );
     p = two_prod( x[1], y[1], e ); 
  f[3] = FPadd_rn( f[3], e );
  f[2] = two_sum ( f[2], p,    e );
  f[3] = FPadd_rn( f[3], e );
  
  f[1] = two_prod( x[0], y[1], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
     p = two_prod( x[1], y[0], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
  f[1] = two_sum ( f[1], p,    e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
  f[0] = two_prod( x[0], y[0], e );
  f[1] = two_sum ( f[1], e,    e );
  f[2] = two_sum ( f[2], e,   e );
  f[3] = FPadd_rn( f[3], e );
  renorm_rand2L<4,3>(f, z);
}
template<>
__host__ __device__ __forceinline__ void truncatedMul<4,2,2>(double const *x, double const *y, double *z){
  double p, e, f[3];
  f[2] = FPmul_rn( x[1], y[1] );
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
    f[2] = FPfma_rn( x[2], y[0], f[2] );
  #else
    f[2] = FPadd_rn( f[2], FPmul_rn( x[2], y[0] ));
  #endif
  f[1] = two_prod( x[0], y[1], e );
  f[2] = FPadd_rn( f[2], e );
     p = two_prod( x[1], y[0], e );
  f[2] = FPadd_rn( f[2], e );
  f[1] = two_sum ( f[1], p,    e );
  f[2] = FPadd_rn( f[2], e );
  f[0] = two_prod( x[0], y[0], e );
  f[1] = two_sum ( f[1], e,    e );
  f[2] = FPadd_rn( f[2], e );
  renorm_rand2L<3,2>(f, z);
}
template<>
__host__ __device__ __forceinline__ void truncatedMul<4,2,4>(double const *x, double const *y, double *z){
  double p, e, f[5];
  f[4] = FPmul_rn( x[3], y[1] );
  f[3] = two_prod( x[3], y[0], e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[2], y[1], e );
  f[4] = FPadd_rn( f[4], e );
  f[3] = two_sum ( f[3], p,    e );
  f[4] = FPadd_rn( f[4], e );
  f[2] = two_prod( x[2], y[0], e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[1], y[1], e );
  f[3] = two_sum ( f[3], e,    e ); 
  f[4] = FPadd_rn( f[4], e );
  f[2] = two_sum ( f[2], p,    e );
  f[3] = two_sum ( f[3], e,    e ); 
  f[4] = FPadd_rn( f[4], e );
  f[1] = two_prod( x[0], y[1], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[1], y[0], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  f[1] = two_sum ( f[1], p,    e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  f[0] = two_prod( x[0], y[0], e );
  f[1] = two_sum ( f[1], e,    e );
  f[2] = two_sum ( f[2], e,   e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  renorm_rand2L<5,4>(f, z);
}

/****************************************************/
template<> __host__ __device__ __forceinline__ void baileyMul_fast<1,1,1>(double const *x, double const *y, double *z){ z[0] = FPmul_rn(x[0], y[0]); }
template<> __host__ __device__ __forceinline__ void baileyMul_fast<1,1,2>(double const *x, double const *y, double *z){ z[0] = two_prod(x[0], y[0], z[1]); }
template<>
__host__ __device__ __forceinline__ void baileyMul_fast<2,1,2>(double const *x, double const *y, double *z){
// DWTimesFP2 with relative error <= 3u^2 or DWTimesFP3 with relative error <= 2u^2
  double ch, cl1, cl2;
  ch = two_prod( x[0], y[0], cl1 );
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
    cl2 = FPfma_rn( x[1], y[0], cl1 );
  #else
    cl2 = FPmul_rn( x[1], y[0] );
    cl2 = FPadd_rn( cl1, cl2 );
  #endif
  z[0] = fast_two_sum( ch, cl2, z[1] );
}
template<>
__host__ __device__ __forceinline__ void baileyMul_fast<2,2,2>(double const *x, double const *y, double *z){
// DWTimesDW1 with relative error <= 7u^2 or DWTimesDW2 with relative error <= 6u^2
  double ch, cl, tl1, tl2;
  ch = two_prod( x[0], y[0], cl );
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
    tl1 = FPmul_rn( x[0], y[1] );
    tl2 = FPfma_rn( x[1], y[0], tl1 );
    cl = FPadd_rn( cl, tl2 );
  #else
    tl1 = FPmul_rn( x[0], y[1] );
    tl2 = FPmul_rn( x[1], y[0] );
    cl = FPadd_rn( cl, FPadd_rn( tl1, tl2 ) );
  #endif
  z[0] = fast_two_sum( ch, cl, z[1] );
}
template<>
__host__ __device__ __forceinline__ void baileyMul_fast<3,1,3>(double const *x, double const *y, double *z){
  double e, f[4];
  f[2] = two_prod( x[2], y[0], f[3] );
  f[1] = two_prod( x[1], y[0], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
  f[0] = two_prod( x[0], y[0], e );
  f[1] = two_sum ( f[1], e,    e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
  fast_renorm2L<4,3>( f, z );
}
template<>
__host__ __device__ __forceinline__ void baileyMul_fast<3,3,3>(double const *x, double const *y, double *z){
  double p, e, f[4];
  f[3] = FPmul_rn( x[1], y[2] );
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
    f[3] = FPfma_rn( x[2], y[1], f[3] );
  #else
    f[3] = FPadd_rn( f[3], FPmul_rn( x[2], y[1] ) );
  #endif
  f[2] = two_prod( x[0], y[2], e );
  f[3] = FPadd_rn( f[3], e );
     p = two_prod( x[1], y[1], e );
  f[3] = FPadd_rn( f[3], e );
  f[2] = two_sum ( f[2], p,    e );
  f[3] = FPadd_rn( f[3], e );
     p = two_prod( x[2], y[0], e );
  f[3] = FPadd_rn( f[3], e );
  f[2] = two_sum ( f[2], p,    e );
  f[3] = FPadd_rn( f[3], e );
  f[1] = two_prod( x[0], y[1], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
     p = two_prod( x[1], y[0], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
  f[1] = two_sum ( f[1], p,    e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
  f[0] = two_prod( x[0], y[0], e );
  f[1] = two_sum ( f[1], e,    e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
  fast_renorm2L<4,3>( f, z );
}
template<>
__host__ __device__ __forceinline__ void baileyMul_fast<4,1,4>(double const *x, double const *y, double *z){
  double e, f[5];
  f[3] = two_prod( x[3], y[0], f[4] );
  f[2] = two_prod( x[2], y[0], e );
  f[3] = two_sum ( f[3], e,    e ); 
  f[4] = FPadd_rn( f[4], e );
  f[1] = two_prod( x[1], y[0], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  f[0] = two_prod( x[0], y[0], e );
  f[1] = two_sum ( f[1], e,    e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  fast_renorm2L<5,4>( f, z );
}
template<>
__host__ __device__ __forceinline__ void baileyMul_fast<4,4,4>(double const *x, double const *y, double *z){
  double p, e, f[5];
  f[4] = FPmul_rn( x[1], y[3] );
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
    f[4] = FPfma_rn( x[2], y[2], f[4] );
    f[4] = FPfma_rn( x[3], y[1], f[4] );
  #else
    f[4] = FPadd_rn( f[4], FPmul_rn( x[2], y[2] ) );
    f[4] = FPadd_rn( f[4], FPmul_rn( x[3], y[1] ) );
  #endif
  f[3] = two_prod( x[0], y[3], e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[1], y[2], e );
  f[4] = FPadd_rn( f[4], e );
  f[3] = two_sum ( f[3], p,    e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[2], y[1], e );
  f[4] = FPadd_rn( f[4], e );
  f[3] = two_sum ( f[3], p,    e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[3], y[0], e );
  f[4] = FPadd_rn( f[4], e );
  f[3] = two_sum ( f[3], p,    e );
  f[4] = FPadd_rn( f[4], e );
  f[2] = two_prod( x[0], y[2], e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[1], y[1], e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  f[2] = two_sum ( f[2], p,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[2], y[0], e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  f[2] = two_sum ( f[2], p,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  f[1] = two_prod( x[0], y[1], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[1], y[0], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  f[1] = two_sum ( f[1], p,    e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  f[0] = two_prod( x[0], y[0], e );
  f[1] = two_sum ( f[1], e,    e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  fast_renorm2L<5,4>( f, z );
}

template<>
__host__ __device__ __forceinline__ void baileyMul_fast<5,1,5>(double const *x, double const *y, double *z){
	double e, f[6];
	f[5] = 0.0;
	f[4] = two_prod(x[4], y[0], e);
	f[5] = FPadd_rn(f[5], e);
	f[3] = two_prod(x[3], y[0], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	f[2] = two_prod(x[2], y[0], e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	f[1] = two_prod(x[1], y[0], e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	f[0] = two_prod(x[0], y[0], e);
	f[1] = two_sum(f[1], e, e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
  fast_renorm2L<6,5>( f, z );
}

template<>
__host__ __device__ __forceinline__ void baileyMul_fast<5,5,5>(double const *x, double const *y, double *z){
	double p, e, f[6];
	f[5] = FPmul_rn(x[1], y[4]);
	#if defined __CUDA_ARCH__ || FP_FAST_FMA
    f[5] = FPfma_rn( x[2], y[3], f[5] );
    f[5] = FPfma_rn( x[3], y[2], f[5] );
    f[5] = FPfma_rn( x[4], y[1], f[5] );
  #else
    f[5] = FPadd_rn(f[5], FPmul_rn(x[2], y[3]));
  	f[5] = FPadd_rn(f[5], FPmul_rn(x[3], y[2]));
  	f[5] = FPadd_rn(f[5], FPmul_rn(x[4], y[1]));
  #endif
	f[4] = two_prod(x[0], y[4], e);
	f[5] = FPadd_rn(f[5], e);
	p = two_prod(x[1], y[3], e);
	f[5] = FPadd_rn(f[5], e);
	f[4] = two_sum(f[4], p, e);
	f[5] = FPadd_rn(f[5], e);
	p = two_prod(x[2], y[2], e);
	f[5] = FPadd_rn(f[5], e);
	f[4] = two_sum(f[4], p, e);
	f[5] = FPadd_rn(f[5], e);
	p = two_prod(x[3], y[1], e);
	f[5] = FPadd_rn(f[5], e);
	f[4] = two_sum(f[4], p, e);
	f[5] = FPadd_rn(f[5], e);
	p = two_prod(x[4], y[0], e);
	f[5] = FPadd_rn(f[5], e);
	f[4] = two_sum(f[4], p, e);
	f[5] = FPadd_rn(f[5], e);
	f[3] = two_prod(x[0], y[3], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	p = two_prod(x[1], y[2], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	f[3] = two_sum(f[3], p, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	p = two_prod(x[2], y[1], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	f[3] = two_sum(f[3], p, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	p = two_prod(x[3], y[0], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	f[3] = two_sum(f[3], p, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	f[2] = two_prod(x[0], y[2], e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	p = two_prod(x[1], y[1], e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	f[2] = two_sum(f[2], p, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	p = two_prod(x[2], y[0], e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	f[2] = two_sum(f[2], p, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	f[1] = two_prod(x[0], y[1], e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	p = two_prod(x[1], y[0], e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	f[1] = two_sum(f[1], p, e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
	f[0] = two_prod(x[0], y[0], e);
	f[1] = two_sum(f[1], e, e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = FPadd_rn(f[5], e);
  fast_renorm2L<6,5>( f, z );
}

template<>
__host__ __device__ __forceinline__ void baileyMul_fast<6,1,6>(double const *x, double const *y, double *z){
	double e, f[7];
	f[6] = 0.0;
	f[5] = two_prod(x[5], y[0], e);
	f[6] = FPadd_rn(f[6], e);
	f[4] = two_prod(x[4], y[0], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[3] = two_prod(x[3], y[0], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[2] = two_prod(x[2], y[0], e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[1] = two_prod(x[1], y[0], e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[0] = two_prod(x[0], y[0], e);
	f[1] = two_sum(f[1], e, e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	fast_renorm2L<7,6>( f, z );
}

template<>
__host__ __device__ __forceinline__ void baileyMul_fast<6,6,6>(double const *x, double const *y, double *z){
	double p, e, f[7];
	f[6] = FPmul_rn(x[1], y[5]);
	#if defined __CUDA_ARCH__ || FP_FAST_FMA
    f[6] = FPfma_rn( x[2], y[4], f[6] );
    f[6] = FPfma_rn( x[3], y[3], f[6] );
    f[6] = FPfma_rn( x[4], y[2], f[6] );
    f[6] = FPfma_rn( x[5], y[1], f[6] );
  #else
    f[6] = FPadd_rn(f[6], FPmul_rn(x[2], y[4]));
  	f[6] = FPadd_rn(f[6], FPmul_rn(x[3], y[3]));
  	f[6] = FPadd_rn(f[6], FPmul_rn(x[4], y[2]));
  	f[6] = FPadd_rn(f[6], FPmul_rn(x[5], y[1]));
  #endif
	f[5] = two_prod(x[0], y[5], e);
	f[6] = FPadd_rn(f[6], e);
	p = two_prod(x[1], y[4], e);
	f[6] = FPadd_rn(f[6], e);
	f[5] = two_sum(f[5], p, e);
	f[6] = FPadd_rn(f[6], e);
	p = two_prod(x[2], y[3], e);
	f[6] = FPadd_rn(f[6], e);
	f[5] = two_sum(f[5], p, e);
	f[6] = FPadd_rn(f[6], e);
	p = two_prod(x[3], y[2], e);
	f[6] = FPadd_rn(f[6], e);
	f[5] = two_sum(f[5], p, e);
	f[6] = FPadd_rn(f[6], e);
	p = two_prod(x[4], y[1], e);
	f[6] = FPadd_rn(f[6], e);
	f[5] = two_sum(f[5], p, e);
	f[6] = FPadd_rn(f[6], e);
	p = two_prod(x[5], y[0], e);
	f[6] = FPadd_rn(f[6], e);
	f[5] = two_sum(f[5], p, e);
	f[6] = FPadd_rn(f[6], e);
	f[4] = two_prod(x[0], y[4], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	p = two_prod(x[1], y[3], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[4] = two_sum(f[4], p, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	p = two_prod(x[2], y[2], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[4] = two_sum(f[4], p, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	p = two_prod(x[3], y[1], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[4] = two_sum(f[4], p, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	p = two_prod(x[4], y[0], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[4] = two_sum(f[4], p, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[3] = two_prod(x[0], y[3], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	p = two_prod(x[1], y[2], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[3] = two_sum(f[3], p, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	p = two_prod(x[2], y[1], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[3] = two_sum(f[3], p, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	p = two_prod(x[3], y[0], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[3] = two_sum(f[3], p, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[2] = two_prod(x[0], y[2], e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	p = two_prod(x[1], y[1], e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[2] = two_sum(f[2], p, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	p = two_prod(x[2], y[0], e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[2] = two_sum(f[2], p, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[1] = two_prod(x[0], y[1], e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	p = two_prod(x[1], y[0], e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[1] = two_sum(f[1], p, e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	f[0] = two_prod(x[0], y[0], e);
	f[1] = two_sum(f[1], e, e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = FPadd_rn(f[6], e);
	fast_renorm2L<7,6>( f, z );
}

template<>
__host__ __device__ __forceinline__ void baileyMul_fast<7,1,7>(double const *x, double const *y, double *z){
	double e, f[8];
	f[7] = 0.0;
	f[6] = two_prod(x[6], y[0], e);
	f[7] = FPadd_rn(f[7], e);
	f[5] = two_prod(x[5], y[0], e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[4] = two_prod(x[4], y[0], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[3] = two_prod(x[3], y[0], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[2] = two_prod(x[2], y[0], e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[1] = two_prod(x[1], y[0], e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[0] = two_prod(x[0], y[0], e);
	f[1] = two_sum(f[1], e, e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	fast_renorm2L<8,7>( f, z );
}

template<>
__host__ __device__ __forceinline__ void baileyMul_fast<7,7,7>(double const *x, double const *y, double *z){
	double p, e, f[8];
	f[7] = FPmul_rn(x[1], y[6]);
	#if defined __CUDA_ARCH__ || FP_FAST_FMA
    f[7] = FPfma_rn( x[2], y[5], f[7] );
    f[7] = FPfma_rn( x[3], y[4], f[7] );
    f[7] = FPfma_rn( x[4], y[3], f[7] );
    f[7] = FPfma_rn( x[5], y[2], f[7] );
    f[7] = FPfma_rn( x[6], y[1], f[7] );
  #else
    f[7] = FPadd_rn(f[7], FPmul_rn(x[2], y[5]));
	  f[7] = FPadd_rn(f[7], FPmul_rn(x[3], y[4]));
	  f[7] = FPadd_rn(f[7], FPmul_rn(x[4], y[3]));
	  f[7] = FPadd_rn(f[7], FPmul_rn(x[5], y[2]));
	  f[7] = FPadd_rn(f[7], FPmul_rn(x[6], y[1]));
  #endif
	f[6] = two_prod(x[0], y[6], e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[1], y[5], e);
	f[7] = FPadd_rn(f[7], e);
	f[6] = two_sum(f[6], p, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[2], y[4], e);
	f[7] = FPadd_rn(f[7], e);
	f[6] = two_sum(f[6], p, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[3], y[3], e);
	f[7] = FPadd_rn(f[7], e);
	f[6] = two_sum(f[6], p, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[4], y[2], e);
	f[7] = FPadd_rn(f[7], e);
	f[6] = two_sum(f[6], p, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[5], y[1], e);
	f[7] = FPadd_rn(f[7], e);
	f[6] = two_sum(f[6], p, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[6], y[0], e);
	f[7] = FPadd_rn(f[7], e);
	f[6] = two_sum(f[6], p, e);
	f[7] = FPadd_rn(f[7], e);
	f[5] = two_prod(x[0], y[5], e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[1], y[4], e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[5] = two_sum(f[5], p, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[2], y[3], e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[5] = two_sum(f[5], p, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[3], y[2], e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[5] = two_sum(f[5], p, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[4], y[1], e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[5] = two_sum(f[5], p, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[5], y[0], e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[5] = two_sum(f[5], p, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[4] = two_prod(x[0], y[4], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[1], y[3], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[4] = two_sum(f[4], p, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[2], y[2], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[4] = two_sum(f[4], p, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[3], y[1], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[4] = two_sum(f[4], p, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[4], y[0], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[4] = two_sum(f[4], p, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[3] = two_prod(x[0], y[3], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[1], y[2], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[3] = two_sum(f[3], p, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[2], y[1], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[3] = two_sum(f[3], p, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[3], y[0], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[3] = two_sum(f[3], p, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[2] = two_prod(x[0], y[2], e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[1], y[1], e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[2] = two_sum(f[2], p, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[2], y[0], e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[2] = two_sum(f[2], p, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[1] = two_prod(x[0], y[1], e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	p = two_prod(x[1], y[0], e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[1] = two_sum(f[1], p, e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	f[0] = two_prod(x[0], y[0], e);
	f[1] = two_sum(f[1], e, e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = FPadd_rn(f[7], e);
	fast_renorm2L<8,7>( f, z );
}

template<>
__host__ __device__ __forceinline__ void baileyMul_fast<8,1,8>(double const *x, double const *y, double *z){
	double e, f[9];
	f[8] = 0.0;
	f[7] = two_prod(x[7], y[0], e);
	f[8] = FPadd_rn(f[8], e);
	f[6] = two_prod(x[6], y[0], e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[5] = two_prod(x[5], y[0], e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[4] = two_prod(x[4], y[0], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[3] = two_prod(x[3], y[0], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[2] = two_prod(x[2], y[0], e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[1] = two_prod(x[1], y[0], e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[0] = two_prod(x[0], y[0], e);
	f[1] = two_sum(f[1], e, e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	fast_renorm2L<9,8>( f, z );
}

template<>
__host__ __device__ __forceinline__ void baileyMul_fast<8,8,8>(double const *x, double const *y, double *z){
	double p, e, f[9];
	f[8] = FPmul_rn(x[1], y[7]);
	#if defined __CUDA_ARCH__ || FP_FAST_FMA
    f[8] = FPfma_rn( x[2], y[6], f[8] );
    f[8] = FPfma_rn( x[3], y[5], f[8] );
    f[8] = FPfma_rn( x[4], y[4], f[8] );
    f[8] = FPfma_rn( x[5], y[3], f[8] );
    f[8] = FPfma_rn( x[6], y[2], f[8] );
    f[8] = FPfma_rn( x[7], y[1], f[8] );
  #else
    f[8] = FPadd_rn(f[8], FPmul_rn(x[2], y[6]));
	  f[8] = FPadd_rn(f[8], FPmul_rn(x[3], y[5]));
	  f[8] = FPadd_rn(f[8], FPmul_rn(x[4], y[4]));
	  f[8] = FPadd_rn(f[8], FPmul_rn(x[5], y[3]));
	  f[8] = FPadd_rn(f[8], FPmul_rn(x[6], y[2]));
	  f[8] = FPadd_rn(f[8], FPmul_rn(x[7], y[1]));
  #endif
	f[7] = two_prod(x[0], y[7], e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[1], y[6], e);
	f[8] = FPadd_rn(f[8], e);
	f[7] = two_sum(f[7], p, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[2], y[5], e);
	f[8] = FPadd_rn(f[8], e);
	f[7] = two_sum(f[7], p, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[3], y[4], e);
	f[8] = FPadd_rn(f[8], e);
	f[7] = two_sum(f[7], p, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[4], y[3], e);
	f[8] = FPadd_rn(f[8], e);
	f[7] = two_sum(f[7], p, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[5], y[2], e);
	f[8] = FPadd_rn(f[8], e);
	f[7] = two_sum(f[7], p, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[6], y[1], e);
	f[8] = FPadd_rn(f[8], e);
	f[7] = two_sum(f[7], p, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[7], y[0], e);
	f[8] = FPadd_rn(f[8], e);
	f[7] = two_sum(f[7], p, e);
	f[8] = FPadd_rn(f[8], e);
	f[6] = two_prod(x[0], y[6], e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[1], y[5], e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[6] = two_sum(f[6], p, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[2], y[4], e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[6] = two_sum(f[6], p, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[3], y[3], e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[6] = two_sum(f[6], p, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[4], y[2], e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[6] = two_sum(f[6], p, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[5], y[1], e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[6] = two_sum(f[6], p, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[6], y[0], e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[6] = two_sum(f[6], p, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[5] = two_prod(x[0], y[5], e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[1], y[4], e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[5] = two_sum(f[5], p, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[2], y[3], e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[5] = two_sum(f[5], p, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[3], y[2], e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[5] = two_sum(f[5], p, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[4], y[1], e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[5] = two_sum(f[5], p, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[5], y[0], e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[5] = two_sum(f[5], p, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[4] = two_prod(x[0], y[4], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[1], y[3], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[4] = two_sum(f[4], p, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[2], y[2], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[4] = two_sum(f[4], p, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[3], y[1], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[4] = two_sum(f[4], p, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[4], y[0], e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[4] = two_sum(f[4], p, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[3] = two_prod(x[0], y[3], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[1], y[2], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[3] = two_sum(f[3], p, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[2], y[1], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[3] = two_sum(f[3], p, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[3], y[0], e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[3] = two_sum(f[3], p, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[2] = two_prod(x[0], y[2], e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[1], y[1], e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[2] = two_sum(f[2], p, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[2], y[0], e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[2] = two_sum(f[2], p, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[1] = two_prod(x[0], y[1], e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	p = two_prod(x[1], y[0], e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[1] = two_sum(f[1], p, e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	f[0] = two_prod(x[0], y[0], e);
	f[1] = two_sum(f[1], e, e);
	f[2] = two_sum(f[2], e, e);
	f[3] = two_sum(f[3], e, e);
	f[4] = two_sum(f[4], e, e);
	f[5] = two_sum(f[5], e, e);
	f[6] = two_sum(f[6], e, e);
	f[7] = two_sum(f[7], e, e);
	f[8] = FPadd_rn(f[8], e);
	fast_renorm2L<9,8>( f, z );
}

template<>
__host__ __device__ __forceinline__ void baileyMul_fast<3,2,2>(double const *x, double const *y, double *z){
  double p, e, f[3];
  f[2] = FPmul_rn( x[1], y[1] );
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
    f[2] = FPfma_rn( x[2], y[0], f[2] );
  #else
    f[2] = FPadd_rn( f[2], FPmul_rn( x[2], y[0] ));
  #endif
  f[1] = two_prod( x[0], y[1], e );
  f[2] = FPadd_rn( f[2], e );
     p = two_prod( x[1], y[0], e );
  f[2] = FPadd_rn( f[2], e );
  f[1] = two_sum ( f[1], p,    e );
  f[2] = FPadd_rn( f[2], e );
  f[0] = two_prod( x[0], y[0], e );
  f[1] = two_sum ( f[1], e,    e );
  f[2] = FPadd_rn( f[2], e );
  fast_renorm2L<3,2>( f, z );
}
template<>
__host__ __device__ __forceinline__ void baileyMul_fast<3,2,3>(double const *x, double const *y, double *z){
  double p, e, f[4];
  f[3] = FPmul_rn( x[2], y[1] ); 
  
  f[2] = two_prod( x[2], y[0], e );
  f[3] = FPadd_rn( f[3], e );
     p = two_prod( x[1], y[1], e ); 
  f[3] = FPadd_rn( f[3], e );
  f[2] = two_sum ( f[2], p,    e );
  f[3] = FPadd_rn( f[3], e );
  
  f[1] = two_prod( x[0], y[1], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
     p = two_prod( x[1], y[0], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
  f[1] = two_sum ( f[1], p,    e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = FPadd_rn( f[3], e );
  f[0] = two_prod( x[0], y[0], e );
  f[1] = two_sum ( f[1], e,    e );
  f[2] = two_sum ( f[2], e,   e );
  f[3] = FPadd_rn( f[3], e );
  fast_renorm2L<4,3>( f, z );
}
template<>
__host__ __device__ __forceinline__ void baileyMul_fast<4,2,2>(double const *x, double const *y, double *z){
  double p, e, f[3];
  f[2] = FPmul_rn( x[1], y[1] );
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
    f[2] = FPfma_rn( x[2], y[0], f[2] );
  #else
    f[2] = FPadd_rn( f[2], FPmul_rn( x[2], y[0] ));
  #endif
  f[1] = two_prod( x[0], y[1], e );
  f[2] = FPadd_rn( f[2], e );
     p = two_prod( x[1], y[0], e );
  f[2] = FPadd_rn( f[2], e );
  f[1] = two_sum ( f[1], p,    e );
  f[2] = FPadd_rn( f[2], e );
  f[0] = two_prod( x[0], y[0], e );
  f[1] = two_sum ( f[1], e,    e );
  f[2] = FPadd_rn( f[2], e );
  fast_renorm2L<3,2>( f, z );
}
template<>
__host__ __device__ __forceinline__ void baileyMul_fast<4,2,4>(double const *x, double const *y, double *z){
  double p, e, f[5];
  f[4] = FPmul_rn( x[3], y[1] );
  f[3] = two_prod( x[3], y[0], e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[2], y[1], e );
  f[4] = FPadd_rn( f[4], e );
  f[3] = two_sum ( f[3], p,    e );
  f[4] = FPadd_rn( f[4], e );
  f[2] = two_prod( x[2], y[0], e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[1], y[1], e );
  f[3] = two_sum ( f[3], e,    e ); 
  f[4] = FPadd_rn( f[4], e );
  f[2] = two_sum ( f[2], p,    e );
  f[3] = two_sum ( f[3], e,    e ); 
  f[4] = FPadd_rn( f[4], e );
  f[1] = two_prod( x[0], y[1], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
     p = two_prod( x[1], y[0], e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  f[1] = two_sum ( f[1], p,    e );
  f[2] = two_sum ( f[2], e,    e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  f[0] = two_prod( x[0], y[0], e );
  f[1] = two_sum ( f[1], e,    e );
  f[2] = two_sum ( f[2], e,   e );
  f[3] = two_sum ( f[3], e,    e );
  f[4] = FPadd_rn( f[4], e );
  fast_renorm2L<5,4>( f, z );
}

/****************************************************/
template<> __host__ __device__ __forceinline__ void baileyMul_renorm<1,1,1>(double const *x, double const *y, double *z){ z[0] = FPmul_rn(x[0], y[0]); }
template<> __host__ __device__ __forceinline__ void baileyMul_renorm<1,1,2>(double const *x, double const *y, double *z){ z[0] = two_prod(x[0], y[0], z[1]); }
template<> 
__host__ __device__ __forceinline__ void baileyMul_renorm<2,1,2>(double const *x, double const *y, double *z){
// DWTimesFP1 or DWTimesFP3 with relative error <= 2u^2
  double ch, cl1, cl2;
  ch = two_prod( x[0], y[0], cl1 );
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
    cl2 = FPfma_rn( x[1], y[0], cl1 );
  #else
    cl2 = FPmul_rn( x[1], y[0] );
    ch = fast_two_sum( ch, cl2, cl2 );
    cl2 = FPadd_rn( cl2, cl1 );
  #endif
  z[0] = fast_two_sum( ch, cl2, z[1] );
}
template<> 
__host__ __device__ __forceinline__ void baileyMul_renorm<2,2,2>(double const *x, double const *y, double *z){
// DWTimesDW1 with relative error <= 7u^2 or DWTimesDW2 with relative error <= 6u^2
  double ch, cl, tl1, tl2;
  ch = two_prod( x[0], y[0], cl );
  #if defined __CUDA_ARCH__ || FP_FAST_FMA
    tl1 = FPmul_rn( x[0], y[1] );
    tl2 = FPfma_rn( x[1], y[0], tl1 );
    cl = FPadd_rn( cl, tl2 );
  #else
    tl1 = FPmul_rn( x[0], y[1] );
    tl2 = FPmul_rn( x[1], y[0] );
    cl = FPadd_rn( cl, FPadd_rn( tl1, tl2 ) );
  #endif
  z[0] = fast_two_sum( ch, cl, z[1] );
}

#endif

