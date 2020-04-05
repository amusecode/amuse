/*
 * gpu_mprec.h
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

#ifndef _gpu_mprec_h
#define _gpu_mprec_h

# define CUDA_SAFE_CALL_NO_SYNC( call) do {                              \
  cudaError err = call;                                                  \
  if( cudaSuccess != err) {                                              \
	  fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
	 	        __FILE__, __LINE__, cudaGetErrorString( err) );              \
		exit(EXIT_FAILURE);                                                  \
  }} while (0)

# define CUDA_SAFE_CALL( call) do {                                      \
  CUDA_SAFE_CALL_NO_SYNC(call);                                          \
  cudaError err = cudaThreadSynchronize();                               \
  if( cudaSuccess != err) {                                              \
    fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
           __FILE__, __LINE__, cudaGetErrorString( err) );               \
    exit(EXIT_FAILURE);                                                  \
  }} while (0)

#define cudaSafeCall(call) CUDA_SAFE_CALL(call)

class gpu_mprec{
private:
	double val;
public:
	//---------constructors--------------
	__device__ gpu_mprec(){}
	__device__ gpu_mprec(const double newVal):val(newVal){}

	//----------geters & setters--------------
	__device__ double getVal(){ return val; }
	__device__ void setVal(const double newVal){ val = newVal; }
	
	template <int prec>
  __device__ friend gpu_mprec loadExpans(multi_prec<prec> const &mp);
  template <int prec>
  __device__ friend void storeExpans(gpu_mprec gmp, multi_prec<prec> &mp);

	//----------friend functions--------------
  template <int L, int K, int R>
  __device__ friend gpu_mprec P_addExpans_safe(gpu_mprec x, gpu_mprec y);
  template <int R>
  __device__ friend gpu_mprec P_addExpans_quick(gpu_mprec x, gpu_mprec y);
  
  template <int K, int L, int R>  
  __device__ friend gpu_mprec P_mulExpans(gpu_mprec x, gpu_mprec y);
};

template <int prec>
__device__ gpu_mprec loadExpans(multi_prec<prec> const &mp){
	if (threadIdx.x < prec)
		return gpu_mprec(mp.getData()[threadIdx.x]);
	else
		return gpu_mprec(0.);
}

template <int prec>
__device__ void storeExpans(gpu_mprec gmp, multi_prec<prec> &mp){
	mp.setElement(gmp.getVal(), threadIdx.x);
}

/** Algorithm for adding FP expansions in parallel
		It requires L<=K and the x dimension of the GPU block size to be K (dim3 blockSize(K,1)) **/
// PAddition_safe with relative error <= 2^{-(prec-1)R} (2^{2R}(1+2^{-p}) + 1/(1-2^{-p+1}) )
template <int L, int K, int R>
__device__ gpu_mprec P_addExpans_safe(gpu_mprec x, gpu_mprec y){
	double a = x.getVal(), b = y.getVal();
	double err = 0, e = 0, s = 0, s_i;
	int i;

	if (threadIdx.x == 0){ s = a; err = b; }
	s = two_sum(s, err, e);
	
	for(i=1; i<L; i++){
		err = __shfl_up(e, 1, R);
		s_i = __shfl(a, i, R);
		if (threadIdx.x == 0) err = s_i;
		s = two_sum(s, err, e);

		err = __shfl_up(e, 1, R);
		s_i = __shfl(b, i, R);
		if (threadIdx.x == 0) err = s_i;
		s = two_sum(s, err, e);
	}

	for(i=L; i<K; i++){
		err = __shfl_up(e, 1, R);
		s_i = __shfl(b, i, R);
		if (threadIdx.x == 0) err = s_i;
		s = two_sum(s, err, e);
	}

	for(i=0; i<R-2; i++){
		err = __shfl_up(e, 1, R);
		if(threadIdx.x == 0) err = 0;
		s = two_sum(s, err, e);
	}
	err = __shfl_up(e, 1, R);
	if(threadIdx.x == 0) err = 0;
	s += err;

	return gpu_mprec(s);
}

// PAddition_quick
template <int R>
__device__ gpu_mprec P_addExpans_quick(gpu_mprec x, gpu_mprec y){
	double a = x.getVal(), b = y.getVal();
	double err, e = 0;

	a = two_sum(a, b, e);
	while(__any(e != 0.)){
		err = __shfl_up(e, 1, R);
		if (threadIdx.x == 0) err = 0;
		a = two_sum(a, err, e);
	}
	if(__any(e != 0.)){
		err = __shfl_up(e, 1, R);
		if (threadIdx.x == 0) err = 0;
		a = a + err;
	}
	
	return gpu_mprec(a);
}

/** Algorithm for multiplying FP expansions in parallel
		It requires L<=K and the x dimension of the GPU block size to be K (dim3 blockSize(K,1)) **/
// PMultiplication with relative error <= 2^{-(prec-1)R} (128/127 (K+L-1)- 129/254 R + 2^{-p-r+2}(R^2-R)(R!)^2 )
template <int K, int L, int R>
__device__ gpu_mprec P_mulExpans(gpu_mprec x, gpu_mprec y){
	double a = x.getVal(), b = y.getVal();
	double p, e_p, e_s, s = 0, r = 0;

	for(int l=0; l<L; l++){
		double b_l = __shfl(b, l, K);
		p = two_prod(a, b_l, e_p);
		s = two_sum(s, p, e_s);

		// save s0
		double tmp = __shfl(s, 0, K);
		if (threadIdx.x == l) r = tmp;
		
		s = __shfl_down(s, 1);
		if (threadIdx.x == K-1) s=0;
		
		// propagate the product's error
		while(__any(e_p != 0.)){
			s = two_sum(s, e_p, e_p); // we lose at thread K-1
			e_p = __shfl_up(e_p, 1, K);
			if(threadIdx.x == 0) e_p = 0;
		}

		//propagate the accumulation error
		while(__any(e_s != 0.)){
			s = two_sum(s, e_s, e_s); //we lose at th K-1
			e_s = __shfl_up(e_s, 1, K);
			if(threadIdx.x == 0) e_s = 0;
		}
	}
	return gpu_mprec(r);
}

#endif

