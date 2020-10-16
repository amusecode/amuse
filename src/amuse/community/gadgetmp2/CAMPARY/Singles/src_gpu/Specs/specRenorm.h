/*
 * specRenorm.h
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

#ifndef _specRenorm_h
#define _specRenorm_h

/****************************************************/
template<> __host__ __device__ __forceinline__ void fast_VecSum<1>(float *x){ }
template<> __host__ __device__ __forceinline__ void fast_VecSum<1>(const float *x, float *r){ r[0] = x[0]; }
template<> __host__ __device__ __forceinline__ void fast_VecSum<2>(float *x){ x[0] = fast_two_sum(x[0], x[1], x[1]); }
template<> __host__ __device__ __forceinline__ void fast_VecSum<2>(const float *x, float *r){ r[0] = fast_two_sum(x[0], x[1], r[1]); }

/****************************************************/
template<> __host__ __device__ __forceinline__ void VecSum<1>(float *x){ }
template<> __host__ __device__ __forceinline__ void VecSum<1>(const float *x, float *r){ r[0] = x[0]; }
template<> __host__ __device__ __forceinline__ void VecSum<2>(float *x){ x[0] = two_sum(x[0], x[1], x[1]); }
template<> __host__ __device__ __forceinline__ void VecSum<2>(const float *x, float *r){ r[0] = two_sum(x[0], x[1], r[1]); }

/****************************************************/
template<> __host__ __device__ __forceinline__ void fast_VecSumErrBranch<1,1>(float *x){ }
template<> __host__ __device__ __forceinline__ void fast_VecSumErrBranch<1,1>(const float *x, float *r){ r[0] = x[0]; }
template<> __host__ __device__ __forceinline__ void fast_VecSumErrBranch<2,1>(float *x){ x[0] = FPadd_rn(x[0], x[1]); }
template<> __host__ __device__ __forceinline__ void fast_VecSumErrBranch<2,1>(const float *x, float *r){ r[0] = FPadd_rn(x[0], x[1]); }
template<> __host__ __device__ __forceinline__ void fast_VecSumErrBranch<2,2>(float *x){ x[0] = fast_two_sum(x[0], x[1], x[1]); }
template<> __host__ __device__ __forceinline__ void fast_VecSumErrBranch<2,2>(const float *x, float *r){ r[0] = fast_two_sum(x[0], x[1], r[1]); }

/****************************************************/
template<> __host__ __device__ __forceinline__ void fast_VecSumErr<1>(float *x){ }
template<> __host__ __device__ __forceinline__ void fast_VecSumErr<1>(const float *x, float *r){ r[0] = x[0]; }
template<> __host__ __device__ __forceinline__ void fast_VecSumErr<2>(float *x){ x[0] = fast_two_sum(x[0], x[1], x[1]); }
template<> __host__ __device__ __forceinline__ void fast_VecSumErr<2>(const float *x, float *r){ r[0] = fast_two_sum(x[0], x[1], r[1]); }

/****************************************************/
template<> __host__ __device__ __forceinline__ void fast_renorm2L<1,1>(float *x){ }
template<> __host__ __device__ __forceinline__ void fast_renorm2L<1,1>(const float *x, float *r){ r[0] = x[0]; }
template<> __host__ __device__ __forceinline__ void fast_renorm2L<2,1>(float *x){ x[0] = FPadd_rn(x[0], x[1]); }
template<> __host__ __device__ __forceinline__ void fast_renorm2L<2,1>(const float *x, float *r){ r[0] = FPadd_rn(x[0], x[1]); }
template<> __host__ __device__ __forceinline__ void fast_renorm2L<2,2>(float *x){ x[0] = fast_two_sum(x[0], x[1], x[1]); }
template<> __host__ __device__ __forceinline__ void fast_renorm2L<2,2>(const float *x, float *r){ r[0] = fast_two_sum(x[0], x[1], r[1]); }
template<> __host__ __device__ __forceinline__ void fast_renorm2L<3,2>(float *x){
	int ptr=0;
	float pr = fast_two_sum(x[1], x[2], x[2]);
	x[0] = fast_two_sum(x[0], pr, x[1]);

	if(x[1] == 0.) pr = x[0];
	else { pr = x[1]; ptr++; }
	x[ptr] = fast_two_sum(pr, x[2], x[ptr+1]);

	for(int i=ptr+1; i<2; i++) x[i] = 0.;
}
template<> __host__ __device__ __forceinline__ void fast_renorm2L<3,2>(const float *x, float *r){
	float f[3];
	int ptr=0;
	float pr = fast_two_sum(x[1], x[2], f[2]);
	f[0] = fast_two_sum(x[0], pr, f[1]);

	if(f[1] == 0.) pr = f[0];
	else { r[0] = f[0]; pr = f[1]; ptr++; }

	for(int i=2; ptr<2 && i<3; i++){
		r[ptr] = fast_two_sum(pr, f[i], pr);
		if(pr == 0.) pr = r[ptr]; else ptr++;
	}

	if(ptr<2 && pr!=0.){ r[ptr] = pr; ptr++; }
	for(int i=ptr; i<2; i++) r[i] = 0.;
}
template<> __host__ __device__ __forceinline__ void fast_renorm2L<4,3>(float *x){
	int ptr=0;
	float pr = fast_two_sum(x[2], x[3], x[3]);
	pr = fast_two_sum(x[1], pr, x[2]);
	x[0] = fast_two_sum(x[0], pr, x[1]);

	if(x[1] == 0.) pr = x[0];
	else { pr = x[1]; ptr++; }
	x[ptr] = fast_two_sum(pr, x[2], pr);
	if(pr == 0.) pr = x[ptr]; else ptr++;
	x[ptr] = fast_two_sum(pr, x[3], x[ptr+1]);

	for(int i=ptr+1; i<3; i++) x[i] = 0.;
}

template<> __host__ __device__ __forceinline__ void fast_renorm2L<4,3>(const float *x, float *r){
	float f[4];
	int ptr=0;
	float pr = fast_two_sum(x[2], x[3], f[3]);
	pr = fast_two_sum(x[1], pr, f[2]);
	f[0] = fast_two_sum(x[0], pr, f[1]);

	if(f[1] == 0.) pr = f[0];
	else { r[0] = f[0]; pr = f[1]; ptr++; }
	r[ptr] = fast_two_sum(pr, f[2], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;

	for(int i=3; ptr<3 && i<4; i++){
		r[ptr] = fast_two_sum(pr, f[i], pr);
		if(pr == 0.) pr = r[ptr]; else ptr++;
	}

	if(ptr<3 && pr!=0.){ r[ptr] = pr; ptr++; }
	for(int i=ptr; i<3; i++) r[i] = 0.;
}
template<> __host__ __device__ __forceinline__ void fast_renorm2L<5,4>(float *x){
	int ptr=0;
	float pr = fast_two_sum(x[3], x[4], x[4]);
	pr = fast_two_sum(x[2], pr, x[3]);
	pr = fast_two_sum(x[1], pr, x[2]);
	x[0] = fast_two_sum(x[0], pr, x[1]);

	if(x[1] == 0.) pr = x[0];
	else { pr = x[1]; ptr++; }
	x[ptr] = fast_two_sum(pr, x[2], pr);
	if(pr == 0.) pr = x[ptr]; else ptr++;
	x[ptr] = fast_two_sum(pr, x[3], pr);
	if(pr == 0.) pr = x[ptr]; else ptr++;
	x[ptr] = fast_two_sum(pr, x[4], x[ptr+1]);

	for(int i=ptr+1; i<4; i++) x[i] = 0.;
}

template<> __host__ __device__ __forceinline__ void fast_renorm2L<5,4>(const float *x, float *r){
	float f[5];
	int ptr=0;
	float pr = fast_two_sum(x[3], x[4], f[4]);
	pr = fast_two_sum(x[2], pr, f[3]);
	pr = fast_two_sum(x[1], pr, f[2]);
	f[0] = fast_two_sum(x[0], pr, f[1]);

	if(f[1] == 0.) pr = f[0];
	else { r[0] = f[0]; pr = f[1]; ptr++; }
	r[ptr] = fast_two_sum(pr, f[2], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;
	r[ptr] = fast_two_sum(pr, f[3], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;

	for(int i=4; ptr<4 && i<5; i++){
		r[ptr] = fast_two_sum(pr, f[i], pr);
		if(pr == 0.) pr = r[ptr]; else ptr++;
	}

	if(ptr<4 && pr!=0.){ r[ptr] = pr; ptr++; }
	for(int i=ptr; i<4; i++) r[i] = 0.;
}

template<>
__host__ __device__ __forceinline__ void fast_renorm2L<6,5>(float *x){
	int ptr=0;
	float pr = fast_two_sum(x[4], x[5], x[5]);
	pr = fast_two_sum(x[3], pr, x[4]);
	pr = fast_two_sum(x[2], pr, x[3]);
	pr = fast_two_sum(x[1], pr, x[2]);
	x[0] = fast_two_sum(x[0], pr, x[1]);

	if(x[1] == 0.) pr = x[0];
	else { pr = x[1]; ptr++; }
	x[ptr] = fast_two_sum(pr, x[2], pr);
	if(pr == 0.) pr = x[ptr]; else ptr++;
	x[ptr] = fast_two_sum(pr, x[3], pr);
	if(pr == 0.) pr = x[ptr]; else ptr++;
	x[ptr] = fast_two_sum(pr, x[4], pr);
	if(pr == 0.) pr = x[ptr]; else ptr++;
	x[ptr] = fast_two_sum(pr, x[5], x[ptr+1]);

	for(int i=ptr+1; i<5; i++) x[i] = 0.;
}

template<>
__host__ __device__ __forceinline__ void fast_renorm2L<6,5>(const float *x, float *r){
	float f[6];
	int ptr=0;
	float pr = fast_two_sum(x[4], x[5], f[5]);
	pr = fast_two_sum(x[3], pr, f[4]);
	pr = fast_two_sum(x[2], pr, f[3]);
	pr = fast_two_sum(x[1], pr, f[2]);
	f[0] = fast_two_sum(x[0], pr, f[1]);

	if(f[1] == 0.) pr = f[0];
	else { r[0] = f[0]; pr = f[1]; ptr++; }
	r[ptr] = fast_two_sum(pr, f[2], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;
	r[ptr] = fast_two_sum(pr, f[3], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;
	r[ptr] = fast_two_sum(pr, f[4], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;

	for(int i=5; ptr<5 && i<6; i++){
		r[ptr] = fast_two_sum(pr, f[i], pr);
		if(pr == 0.) pr = r[ptr]; else ptr++;
	}

	if(ptr<5 && pr!=0.){ r[ptr] = pr; ptr++; }
	for(int i=ptr; i<5; i++) r[i] = 0.;
}

template<>
__host__ __device__ __forceinline__ void fast_renorm2L<7,6>(float *x){
	int ptr=0;
	float pr = fast_two_sum(x[5], x[6], x[6]);
	pr = fast_two_sum(x[4], pr, x[5]);
	pr = fast_two_sum(x[3], pr, x[4]);
	pr = fast_two_sum(x[2], pr, x[3]);
	pr = fast_two_sum(x[1], pr, x[2]);
	x[0] = fast_two_sum(x[0], pr, x[1]);

	if(x[1] == 0.) pr = x[0];
	else { pr = x[1]; ptr++; }
	x[ptr] = fast_two_sum(pr, x[2], pr);
	if(pr == 0.) pr = x[ptr]; else ptr++;
	x[ptr] = fast_two_sum(pr, x[3], pr);
	if(pr == 0.) pr = x[ptr]; else ptr++;
	x[ptr] = fast_two_sum(pr, x[4], pr);
	if(pr == 0.) pr = x[ptr]; else ptr++;
	x[ptr] = fast_two_sum(pr, x[5], pr);
	if(pr == 0.) pr = x[ptr]; else ptr++;
	x[ptr] = fast_two_sum(pr, x[6], x[ptr+1]);

	for(int i=ptr+1; i<6; i++) x[i] = 0.;
}

template<>
__host__ __device__ __forceinline__ void fast_renorm2L<7,6>(const float *x, float *r){
	float f[7];
	int ptr=0;
	float pr = fast_two_sum(x[5], x[6], f[6]);
	pr = fast_two_sum(x[4], pr, f[5]);
	pr = fast_two_sum(x[3], pr, f[4]);
	pr = fast_two_sum(x[2], pr, f[3]);
	pr = fast_two_sum(x[1], pr, f[2]);
	f[0] = fast_two_sum(x[0], pr, f[1]);

	if(f[1] == 0.) pr = f[0];
	else { r[0] = f[0]; pr = f[1]; ptr++; }
	r[ptr] = fast_two_sum(pr, f[2], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;
	r[ptr] = fast_two_sum(pr, f[3], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;
	r[ptr] = fast_two_sum(pr, f[4], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;
	r[ptr] = fast_two_sum(pr, f[5], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;

	for(int i=6; ptr<6 && i<7; i++){
		r[ptr] = fast_two_sum(pr, f[i], pr);
		if(pr == 0.) pr = r[ptr]; else ptr++;
	}

	if(ptr<6 && pr!=0.){ r[ptr] = pr; ptr++; }
	for(int i=ptr; i<6; i++) r[i] = 0.;
}

/****************************************************/
template<>
__host__ __device__ __forceinline__ void renorm2L_4Add1<2,2>(const float *x, const float y, float *r){
	float f[3];
	int ptr=0;
	float pr = two_sum(x[1], y, f[2]);
	f[0] = two_sum(x[0], pr, f[1]);

	if(f[1] == 0.) pr = f[0];
	else { r[0] = f[0]; pr = f[1]; ptr++; }

	for(int i=2; ptr<2 && i<3; i++){
		r[ptr] = fast_two_sum(pr, f[i], pr);
		if(pr == 0.) pr = r[ptr]; else ptr++;
	}

	if(ptr<2 && pr!=0.){ r[ptr] = pr; ptr++; }
	for(int i=ptr; i<2; i++) r[i] = 0.;
}

template<>
__host__ __device__ __forceinline__ void renorm2L_4Add1<3,3>(const float *x, const float y, float *r){
	float f[4];
	int ptr=0;
	float pr = two_sum(x[2], y, f[3]);
	pr = two_sum(x[1], pr, f[2]);
	f[0] = two_sum(x[0], pr, f[1]);

	if(f[1] == 0.) pr = f[0];
	else { r[0] = f[0]; pr = f[1]; ptr++; }
	r[ptr] = fast_two_sum(pr, f[2], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;

	for(int i=3; ptr<3 && i<4; i++){
		r[ptr] = fast_two_sum(pr, f[i], pr);
		if(pr == 0.) pr = r[ptr]; else ptr++;
	}

	if(ptr<3 && pr!=0.){ r[ptr] = pr; ptr++; }
	for(int i=ptr; i<3; i++) r[i] = 0.;
}

template<>
__host__ __device__ __forceinline__ void renorm2L_4Add1<4,4>(const float *x, const float y, float *r){
	float f[5];
	int ptr=0;
	float pr = two_sum(x[3], y, f[4]);
	pr = two_sum(x[2], pr, f[3]);
	pr = two_sum(x[1], pr, f[2]);
	f[0] = two_sum(x[0], pr, f[1]);

	if(f[1] == 0.) pr = f[0];
	else { r[0] = f[0]; pr = f[1]; ptr++; }
	r[ptr] = fast_two_sum(pr, f[2], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;
	r[ptr] = fast_two_sum(pr, f[3], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;

	for(int i=4; ptr<4 && i<5; i++){
		r[ptr] = fast_two_sum(pr, f[i], pr);
		if(pr == 0.) pr = r[ptr]; else ptr++;
	}

	if(ptr<4 && pr!=0.){ r[ptr] = pr; ptr++; }
	for(int i=ptr; i<4; i++) r[i] = 0.;
}

template<>
__host__ __device__ __forceinline__ void renorm2L_4Add1<5,5>(const float *x, const float y, float *r){
	float f[6];
	int ptr=0;
	float pr = two_sum(x[4], y, f[5]);
	pr = two_sum(x[3], pr, f[4]);
	pr = two_sum(x[2], pr, f[3]);
	pr = two_sum(x[1], pr, f[2]);
	f[0] = two_sum(x[0], pr, f[1]);

	if(f[1] == 0.) pr = f[0];
	else { r[0] = f[0]; pr = f[1]; ptr++; }
	r[ptr] = fast_two_sum(pr, f[2], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;
	r[ptr] = fast_two_sum(pr, f[3], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;
	r[ptr] = fast_two_sum(pr, f[4], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;

	for(int i=5; ptr<5 && i<6; i++){
		r[ptr] = fast_two_sum(pr, f[i], pr);
		if(pr == 0.) pr = r[ptr]; else ptr++;
	}

	if(ptr<5 && pr!=0.){ r[ptr] = pr; ptr++; }
	for(int i=ptr; i<5; i++) r[i] = 0.;
}

template<>
__host__ __device__ __forceinline__ void renorm2L_4Add1<6,6>(const float *x, const float y, float *r){
	float f[7];
	int ptr=0;
	float pr = two_sum(x[5], y, f[6]);
	pr = two_sum(x[4], pr, f[5]);
	pr = two_sum(x[3], pr, f[4]);
	pr = two_sum(x[2], pr, f[3]);
	pr = two_sum(x[1], pr, f[2]);
	f[0] = two_sum(x[0], pr, f[1]);

	if(f[1] == 0.) pr = f[0];
	else { r[0] = f[0]; pr = f[1]; ptr++; }
	r[ptr] = fast_two_sum(pr, f[2], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;
	r[ptr] = fast_two_sum(pr, f[3], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;
	r[ptr] = fast_two_sum(pr, f[4], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;
	r[ptr] = fast_two_sum(pr, f[5], pr);
	if(pr == 0.) pr = r[ptr]; else ptr++;

	for(int i=6; ptr<6 && i<7; i++){
		r[ptr] = fast_two_sum(pr, f[i], pr);
		if(pr == 0.) pr = r[ptr]; else ptr++;
	}

	if(ptr<6 && pr!=0.){ r[ptr] = pr; ptr++; }
	for(int i=ptr; i<6; i++) r[i] = 0.;
}

/****************************************************/
template<> __host__ __device__ __forceinline__ void renorm2L<1,1>(float *x){ }
template<> __host__ __device__ __forceinline__ void renorm2L<1,1>(const float *x, float *r){ r[0] = x[0]; }
template<> __host__ __device__ __forceinline__ void renorm2L<2,1>(float *x){ x[0] = FPadd_rn(x[0], x[1]); }
template<> __host__ __device__ __forceinline__ void renorm2L<2,1>(const float *x, float *r){ r[0] = FPadd_rn(x[0], x[1]); }
template<> __host__ __device__ __forceinline__ void renorm2L<2,2>(float *x){ x[0] = two_sum(x[0], x[1], x[1]); }
template<> __host__ __device__ __forceinline__ void renorm2L<2,2>(const float *x, float *r){ r[0] = two_sum(x[0], x[1], r[1]); }

/****************************************************/
template<> __host__ __device__ __forceinline__ void renorm_rand2L<1,1>(float *x){ }
template<> __host__ __device__ __forceinline__ void renorm_rand2L<1,1>(const float *x, float *r){ r[0] = x[0]; }
template<> __host__ __device__ __forceinline__ void renorm_rand2L<2,1>(float *x){ x[0] = FPadd_rn(x[0], x[1]); }
template<> __host__ __device__ __forceinline__ void renorm_rand2L<2,1>(const float *x, float *r){ r[0] = FPadd_rn(x[0], x[1]); }
template<> __host__ __device__ __forceinline__ void renorm_rand2L<2,2>(float *x){ x[0] = two_sum(x[0], x[1], x[1]); }
template<> __host__ __device__ __forceinline__ void renorm_rand2L<2,2>(const float *x, float *r){ r[0] = two_sum(x[0], x[1], r[1]); }

#endif

