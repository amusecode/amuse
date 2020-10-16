/*
 * helper_fct.h
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

#ifndef _helper_fct_h
#define _helper_fct_h
#if defined(__CUDACC__)
#include <curand_kernel.h>
#endif

/*#define ROUND_DOWN  FE_DOWNWARD
#define ROUND_UP    FE_UPWARD
#define ROUND_NEAR  FE_TONEAREST
#define ROUND_ZERO  FE_TOWARDZERO
#define setRound    fesetround
void setRoundDown() { setRound(ROUND_DOWN); }
void setRoundUp  () { setRound(ROUND_UP);   }
void setRoundNear() { setRound(ROUND_NEAR); }
void setRoundZero() { setRound(ROUND_ZERO); }

#define EPS 0x1p-53
#define ETA 0x0.0000000000001p-1022
#define DBL_MIN 2.2250738585072014e-308*/

/***************************************************************/
/******************** Print array functions ********************/
/***************************************************************/
template<int size>
__host__ __device__ static void printArray(const float *p){
	printf("Array size %d,\n ",size);
	for(int i=0; i<size; i++) printf("%e; ", p[i]); 
	printf("------\n");
}

template<int size>
__host__ __device__ static void printArrayBin(const float *p){
	printf("Array size %d,\n ",size);
	for(int i=0; i<size; i++) printf("%a; ", p[i]); 
	printf("------\n");
}

/***************************************************************/
/******************* Unified random function *******************/
/***************************************************************/
template<typename T>
__device__ static T unifRand_gpu(){ }

template<typename T>
__device__ static T unifRand_gpu(T a, T b){
  return (b-a) * unifRand_gpu<T>() + a;
}

/***************************************************/
template<>
__device__ float unifRand_gpu<float>(){
 	curandState s;
	unsigned int seed = (long int)&s;
	curand_init(seed,0,0,&s);
  return (float) curand_uniform_double(&s);
}
template<>
__device__ int unifRand_gpu<int>(){
  curandState s;
  unsigned int seed = (long int)&s;
  curand_init(seed,0,0,&s);
  return floor(curand_uniform(&s));
}
/***************************************************/
__host__ void seed(){
  srand(time(0));
}
template<typename T>
__host__ T unifRand_host(){
	return rand() / T(RAND_MAX);
}

template<typename T>
__host__ T unifRand_host(T a, T b){
	return (b-a) * unifRand_host<T>() + a;
}
template<>
__host__ int unifRand_host<int>(int a, int b){
	float aux = (float) rand() / (float)RAND_MAX;
	aux *= (b-a);
  return (int)a + aux;
}

template<typename T>
__host__ __device__ T unifRand(T a, T b){
  #if defined(__CUDA_ARCH__)
  	return unifRand_gpu<T>(a, b);
	#else
  	return unifRand_host<T>(a, b);
  #endif
}

/***************************************************************/
/**************** Generate expansion functions *****************/
/***************************************************************/
__host__ __device__ static int ulpExp(float x){
	int e;
	frexp(x, &e); // returns value in [0.5,1)
	return (e-1)-24+1;
}

template<int size>
__host__ __device__ static void genExpans_P(float *res){
  res[0] = unifRand<float>(-1.0, 1.0);
  int exp0 = 1;
  res[0] = ldexp(res[0], exp0);
  for(int i=1; i<size; i++){
  	res[i] = unifRand<float>(-1.0, 1.0);
  	exp0 = unifRand<int>(1,5);
  	res[i] = ldexp(res[i], ulpExp(res[i-1]) - exp0);
  }
}
template<> __host__ __device__ void genExpans_P<2>(float *res){
  res[0] = unifRand<float>(-1.0, 1.0);
  int exp0 = 1;
  res[0] = ldexp(res[0], exp0);
 	res[1] = unifRand<float>(-1.0, 1.0);
  exp0 = unifRand<int>(1,5);
  res[1] = ldexp(res[1], ulpExp(res[0]) - exp0);
  float s = res[0] + res[1];
  float z = s -res[0];
  res[0] = s;
  res[1] = res[1] -z;
//  res[0] = fast_two_sum(res[0], res[1], res[1]);
}

template<int size>
__host__ __device__ static void genExpans_ulp(float *res){
	res[0] = unifRand<float>(-1.0, 1.0);
  int exp0 = unifRand<int>(-2, 2);
  res[0] = ldexp(res[0], exp0);

	int eq;
  for(int i=1; i<size; i++){
		eq = unifRand<int>(0,100) % 2;
		if(eq == 0){
			res[i] = ldexp(1.0, ulpExp(res[i-1]));
		}else{
      res[i] = unifRand<float>(-1.0, 1.0);
      exp0 = unifRand<int>(1,5);
      res[i] = ldexp(res[i], ulpExp(res[i-1]) - exp0);
		}
  }
}
template<> __host__ __device__ void genExpans_ulp<2>(float *res){
  res[0] = unifRand<float>(-1.0, 1.0);
  int exp0 = 1;
  res[0] = ldexp(res[0], exp0);
 	res[1] = unifRand<float>(-1.0, 1.0);
  exp0 = unifRand<int>(1,5);
  res[1] = ldexp(res[1], ulpExp(res[0]) - exp0);
  float s = res[0] + res[1];
  float z = s -res[0];
  res[0] = s;
  res[1] = res[1] -z;
//  res[0] = fast_two_sum(res[0], res[1], res[1]);
}

template<int size>
__host__ __device__ static void genExpans_ulp(float *res, int order){
  res[0] = unifRand<float>(0.5, 1.0);
  res[0] = ldexp(res[0], order);

  int eq, exp0;
  for(int i=1; i<size; i++){
    eq = unifRand<int>(0,100) % 2;
    if(eq == 0){
      res[i] = ldexp(1.0, ulpExp(res[i-1]));
    }else{
      res[i] = unifRand<float>(-1.0, 1.0);
      exp0 = unifRand<int>(1,5);
      res[i] = ldexp(res[i], ulpExp(res[i-1]) - exp0);
    }
  }
}
template<> __host__ __device__ void genExpans_ulp<2>(float *res, int order){
  res[0] = unifRand<float>(0.5, 1.0);
  res[0] = ldexp(res[0], order);
  int eq, exp0;
  eq = unifRand<int>(0,100) % 2;
  if(eq == 0){
    res[1] = ldexp(1.0, ulpExp(res[0]));
  }else{
    res[1] = unifRand<float>(-1.0, 1.0);
    exp0 = unifRand<int>(1,5);
    res[1] = ldexp(res[1], ulpExp(res[0]) - exp0);
  }
  float s = res[0] + res[1];
  float z = s -res[0];
  res[0] = s;
  res[1] = res[1] -z;
//  res[0] = fast_two_sum(res[0], res[1], res[1]);
}

template<int size>
__host__ __device__ static void genExpans_ulp(float *res, float first){
  res[0] = first;

  int eq, exp0;
  for(int i=1; i<size; i++){
    eq = unifRand<int>(0,100) % 2;
    if(eq == 0){
      res[i] = ldexp(1.0, ulpExp(res[i-1]));
    }else{
      res[i] = unifRand<float>(-1.0, 1.0);
      exp0 = unifRand<int>(1,5);
      res[i] = ldexp(res[i], ulpExp(res[i-1]) - exp0);
    }
  }
}
template<> __host__ __device__ void genExpans_ulp<2>(float *res, float first){
  res[0] = first;
  int eq, exp0;
  eq = unifRand<int>(0,100) % 2;
  if(eq == 0){
    res[1] = ldexp(1.0, ulpExp(res[0]));
  }else{
    res[1] = unifRand<float>(-1.0, 1.0);
    exp0 = unifRand<int>(1,5);
    res[1] = ldexp(res[1], ulpExp(res[0]) - exp0);
  }
  float s = res[0] + res[1];
  float z = s -res[0];
  res[0] = s;
  res[1] = res[1] -z;
//  res[0] = fast_two_sum(res[0], res[1], res[1]);
}

template<int size>
__host__ __device__ static void genOverlappingExpans(float *res){
  res[0] = unifRand<float>(-1.0,1.0);
  int exp0 = 1;
  res[0] = ldexp(res[0],exp0);
  for(int i=1; i<size; i++){
  	res[i] = unifRand<float>(-1.0,1.0);
  	exp0 = unifRand<int>(8,15);
  	res[i] = ldexp(res[i],ulpExp(res[i-1]) + 53 - exp0);
  }
}

#endif
