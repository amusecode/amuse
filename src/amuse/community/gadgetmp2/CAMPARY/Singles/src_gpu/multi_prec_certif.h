/*
 * multi_prec_certif.h
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

#ifndef _multi_prec_h
#define _multi_prec_h

#include <stdio.h>
#include <math.h>
#include <fenv.h>
#include <cstdio>
#include <string.h>
#include <iostream>
#include <stdlib.h>
using namespace std;

#include "errorFreeTransf.h"
#include "helper_fct.h"
#include "renorm.h"
#include "Specs/specRenorm.h"
#include "addition.h"
#include "Specs/specAddition.h"
#include "multiplication.h"
#include "Specs/specMultiplication.h"
#include "newton.h"
#include "Specs/specNewton.h"

/**forward declarations**/
template <int prec> class multi_prec ;

/**template friends**/
template <int prec>
__host__ __device__ multi_prec<prec> abs(const multi_prec<prec> &mp);
template <int prec>
__host__ __device__ multi_prec<prec> sqrt(const multi_prec<prec> &mp);
template <int prec>
__host__ __device__ multi_prec<prec> invSqrt(const multi_prec<prec> &mp);

template <int prec>
__host__ __device__ void renorm(multi_prec<prec> &mp);
template <int prec>
__host__ __device__ void renorm_rand(multi_prec<prec> &mp);
template <int prec>
__host__ __device__ void renorm_2ndL(multi_prec<prec> &mp);

//********************* ADD *********************//
template <int pR, int p1, int p2>
__host__ __device__ void    certifAddExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 );
template <int pR, int p1>
__host__ __device__ void  certifAddExpans_d( multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val );
template <int pR, int p1, int p2>
__host__ __device__ void     truncAddExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 );
template <int pR, int p1>
__host__ __device__ void   truncAddExpans_d( multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val );
template <int pR, int p1, int p2>
__host__ __device__ void   QD_LikeAddExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 );
template <int pR, int p1>
__host__ __device__ void QD_LikeAddExpans_d( multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val );

//********************* MUL *********************//
template <int pR, int p1, int p2>
__host__ __device__ void    certifMulExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 );
template <int pR, int p1>
__host__ __device__ void  certifMulExpans_d( multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val );
template <int pR, int p1, int p2>
__host__ __device__ void     truncMulExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 );
template <int pR, int p1>
__host__ __device__ void   truncMulExpans_d( multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val );
template <int pR, int p1, int p2>
__host__ __device__ void   QD_LikeMulExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 );
template <int pR, int p1>
__host__ __device__ void QD_LikeMulExpans_d( multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val );

//********************* DIV *********************//
template <int pR, int p1>
__host__ __device__ void   invExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1 );
template <int pR>
__host__ __device__ void invExpans_d( multi_prec<pR> &res, const float val );
template <int pR, int p1, int p2>
__host__ __device__ void   divExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 );
template <int pR, int p1>
__host__ __device__ void divExpans_d( multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val );
template <int pR, int p1>
__host__ __device__ void divExpans_d( multi_prec<pR> &res, const float val, const multi_prec<p1> &mp1 );

//********************* SQRT *********************//
template <int pR, int p1>
__host__ __device__ void      invSqrtExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1 );
template <int pR>
__host__ __device__ void    invSqrtExpans_d( multi_prec<pR> &res, const float val );
template <int pR, int p1>
__host__ __device__ void   sqrtNewtonExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1 );
template <int pR>
__host__ __device__ void sqrtNewtonExpans_d( multi_prec<pR> &res, const float val );
template <int pR, int p1>
__host__ __device__ void    sqrtHeronExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1 );
template <int pR>
__host__ __device__ void  sqrtHeronExpans_d( multi_prec<pR> &res, const float val );

template <int prec>
class multi_prec{
private:
	float data[prec];
public:
	//---------constructors--------------
	__host__ __device__ multi_prec(){};
	__host__ __device__ multi_prec(const float x);
	__host__ __device__ multi_prec(const float *datap, const int precp);
	__host__ __device__ multi_prec(const char *s);
	template<int precP>	__host__ __device__ multi_prec(const multi_prec<precP> *mp);
	__host__ __device__ multi_prec(const multi_prec<prec> *mp);

	//----------geters & setters--------------
	__host__ __device__ int getPrec() const;
	__host__ __device__ const float* getData() const;
	__host__ __device__ void setData(const float *datap, const int precp);
	__host__ __device__ void setData(const float *datap);
	__host__ __device__ void setElement(const float datap, const int index);

 	/**pretty print**/
 	__host__ __device__ void prettyPrint(){
  	printf("Prec = %d\n", prec); 
  	for(int i=0; i<prec; i++) printf("   Data[%d] = %e\n",i,data[i]);
 	}
 	__host__ __device__ void prettyPrintBin(){
  	printf("Prec = %d\n", prec); 
  	for(int i=0; i<prec; i++)	printf("   Data[%d] = %a;\n",i,data[i]);
 	}
  __host__ __device__ void prettyPrintBin_UnevalSum(){
  	for(int i=0; i<prec-1; i++) printf("%a + ", data[i]);
  	printf("%a;", data[prec-1]);
  }
  
  __host__ __device__ char* prettyPrintBF(){
		size_t needed = snprintf(NULL, 0, "%e", data[0]);
		if (prec>1)
			for(int i=1; i<prec; i++)	needed += snprintf(NULL, 0, "%e ", data[i]);
    
		char *ch;
		ch =(char *) malloc((2*needed++)*sizeof(char));
     
		sprintf(ch, "%e ", data[0]);
		if (prec>1)
			for(int i=1; i<prec; i++) sprintf(&ch[strlen(ch)], "%e ", data[i]);

		sprintf(&ch[strlen(ch)], "%c", '\0');
		return ch;
 	}

	//----------operators--------------
	/** operator [] overloading  **/
	__host__ __device__ float operator [](const int i)const {return data[i];}

	/** Puts the value other in data[0] of the current multi_prec variable  **/
	__host__ __device__ multi_prec& operator =(const float other){
		setData(&other, 1);
		return *this;
	}
	/** Transfers the multi_prec parameter other in the current multi_prec variable **/
	__host__ __device__ multi_prec<prec>& operator =(const multi_prec<prec>& other){
		if(this != &other) setData(other.getData(), prec); /* check against self assingment */
		return *this;
	}
	template <int pS>
	__host__ __device__ multi_prec& operator =(const multi_prec<pS>& other){
	 	setData(other.getData(),pS);
	 	return *this;
	}
	__host__ __device__ multi_prec& operator =(const char *s){
		if (read(s, *this)){
			printf("(qd_real::operator=): INPUT ERROR.");
    	*this = 0.0;
  	}
  	return *this;
	}

	/** Equality overloading**/
	template <int pS> __host__ __device__ bool operator ==(const multi_prec<pS> &mp2) const;
	__host__ __device__ bool operator ==(const multi_prec<prec> &mp2) const;
	__host__ __device__ bool operator ==(const float &mp2) const;
	template <int pS> __host__ __device__ bool operator !=(const multi_prec<pS> &mp2) const;
	__host__ __device__ bool operator !=(const multi_prec<prec> &mp2) const;
	__host__ __device__ bool operator !=(const float &mp2) const;
	template <int pS> __host__ __device__ bool operator < (const multi_prec<pS> &mp2) const;
	__host__ __device__ bool operator < (const multi_prec<prec> &mp2) const;
	__host__ __device__ bool operator < (const float &mp2) const;
	template <int pS> __host__ __device__ bool operator <=(const multi_prec<pS> &mp2) const;
	__host__ __device__ bool operator <=(const multi_prec<prec> &mp2) const;
	__host__ __device__ bool operator <=(const float &mp2) const;
	template <int pS> __host__ __device__ bool operator > (const multi_prec<pS> &mp2) const;
	__host__ __device__ bool operator > (const multi_prec<prec> &mp2) const;
	__host__ __device__ bool operator > (const float &mp2) const;
	template <int pS> __host__ __device__ bool operator >=(const multi_prec<pS> &mp2) const;
	__host__ __device__ bool operator >=(const multi_prec<prec> &mp2) const;
	__host__ __device__ bool operator >=(const float &mp2) const;

	/**Random expansion generation **/
	__host__ __device__ void randomInit_P(){ genExpans_P<prec>(data);}
  __host__ __device__ void randomInit_ulp(){ genExpans_ulp<prec>(data);}
  __host__ __device__ void randomInit_ulp(int order){ genExpans_ulp<prec>(data, order);}
  __host__ __device__ void randomInit_ulp(float first){ genExpans_ulp<prec>(data, first);}

	/**Unary - overloading**/
	__host__ __device__ friend multi_prec operator -(const multi_prec<prec> &mp1){
		multi_prec<prec> res;
		for(int i=0; i<prec; i++) res.data[i] = -mp1.getData()[i];
		return res;
	}
/////////////////////////////////////////////////////////////////////////
	/** operator += overloading  **/
	template <int pS>
	__host__ __device__ multi_prec<prec>& operator +=(const multi_prec<pS> &mp){
		truncAddExpans<prec,prec,pS>( *this, *this, mp );
		return *this;
	}
	__host__ __device__ multi_prec<prec>& operator +=(const float val){
		truncAddExpans_d<prec,prec>( *this, *this, val );
		return *this;
	}

  /** operator -= overloading  **/
  template <int pS>
  __host__ __device__ multi_prec<prec>& operator -=(const multi_prec<pS> &mp){
    truncAddExpans<prec,prec,pS>( *this, *this, -mp );
    return *this;
  }
  __host__ __device__ multi_prec<prec>& operator -=(const float val){
    truncAddExpans_d<prec,prec>( *this, *this, -val );
    return *this;
  }

  /** operator *= overloading  **/
  template <int pS>
  __host__ __device__ multi_prec<prec> operator *=(const multi_prec<pS> &mp){
    truncMulExpans<prec,prec,pS>( *this, *this, mp );
    return *this;
  }
  __host__ __device__ multi_prec<prec>& operator *=(const float val){
    truncMulExpans_d<prec,prec>( *this, *this, val );
    return *this;
  }

	/** operator /= overloading  **/
	template <int pS>
	__host__ __device__ multi_prec<prec>& operator /=(const multi_prec<pS> &mp){
		divExpans<prec,prec,pS>( *this, *this, mp );
		return *this;
	}
	__host__ __device__ multi_prec<prec>& operator /=(const float val){
		divExpans_d<prec,prec>( *this, *this, val );
		return *this;
	}

	/**absolute value **/
	__host__ __device__ friend multi_prec<prec> abs<>(const multi_prec<prec> &mp);
	__host__ __device__ friend multi_prec<prec> sqrt<>(const multi_prec<prec> &mp);
	__host__ __device__ friend multi_prec<prec> invSqrt<>(const multi_prec<prec> &mp);
	
	template<int pS> __host__ __device__ friend int read(const char *s, multi_prec<pS> &mp);
	
	template <int pR, int p1, int p2>
	__host__ __device__ friend multi_prec<pR> max( const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 );
	template <int pR, int p1, int p2>
	__host__ __device__ friend multi_prec<pR> min( const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 );
	
	__host__ __device__ friend void renorm<>(multi_prec<prec> &mp);
	__host__ __device__ friend void renorm_rand<>(multi_prec<prec> &mp);
	__host__ __device__ friend void renorm_2ndL<>(multi_prec<prec> &mp);

	/** operator + overloading  **/
	template <int pS>
	__host__ __device__ friend multi_prec operator +(const multi_prec<prec> &mp1, const multi_prec<pS> &mp2){
		multi_prec<(prec>pS)?prec:pS> res(mp1.getData(), prec);
		return res += mp2;
	}
	__host__ __device__ friend multi_prec operator +(const multi_prec<prec> &mp1, const float val){
		multi_prec<prec> res(mp1.getData(), prec);
		return res += val;
	}
	__host__ __device__ friend multi_prec operator +(const float val, const multi_prec<prec> &mp1){
		multi_prec<prec> res(mp1.getData(), prec);
		return res += val;
	}

	/** operator - overloading  **/
	template <int pS>
	__host__ __device__ friend multi_prec operator -(const multi_prec<prec> &mp1, const multi_prec<pS> &mp2){
		multi_prec<(prec>pS)?prec:pS> res(mp1.getData(), prec);
		return res -= mp2;
	}
	__host__ __device__ friend multi_prec operator -(const multi_prec<prec> &mp1, const float val){
		multi_prec<prec> res(mp1.getData(), prec);
		return res -= val;
	}
	__host__ __device__ friend multi_prec operator -(const float val, const multi_prec<prec> &mp1){
		multi_prec<prec> res(val);
		return res -= mp1;
	}

	/** operator * overloading  **/
	template <int pS>
	__host__ __device__ friend multi_prec operator *(const multi_prec<prec> &mp1, const multi_prec<pS> &mp2){
		multi_prec<(prec>pS)?prec:pS> res(mp1.getData(), prec);
		return res *= mp2;
	}
	__host__ __device__ friend multi_prec operator *(const multi_prec<prec> &mp1, const float val){
		multi_prec<prec> res(mp1.getData(), prec);
		return res *= val;
	}
	__host__ __device__ friend multi_prec operator *(const float val, const multi_prec<prec> &mp1){
		multi_prec<prec> res(mp1.getData(), prec);
		return res *= val;
	}

	/** operator / overloading  **/
  template <int pS>
  __host__ __device__ friend multi_prec operator /(const multi_prec<prec> &mp1, const multi_prec<pS> &mp2){
    multi_prec<(prec>pS)?prec:pS> res(mp1.getData(), prec);
    return res /= mp2;
  }
  __host__ __device__ friend multi_prec operator /(const multi_prec<prec> &mp1, const float val){
    multi_prec<prec> res(mp1.getData(), prec);
    return res /= val;
  }
  __host__ __device__ friend multi_prec operator /(const float val, const multi_prec<prec> &mp1){
    multi_prec<prec> res(val);
    return res /= mp1;
  }

	//********************* ADD *********************//
	template <int pR, int p1, int p2>
	__host__ __device__ friend void    certifAddExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 );
	template <int pR, int p1>
	__host__ __device__ friend void  certifAddExpans_d( multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val );
	template <int pR, int p1, int p2>
  __host__ __device__ friend void     truncAddExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 );
  template <int pR, int p1>
  __host__ __device__ friend void   truncAddExpans_d( multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val );
	template <int pR, int p1, int p2>
	__host__ __device__ friend void   QD_LikeAddExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 );
	template <int pR, int p1>
	__host__ __device__ friend void QD_LikeAddExpans_d( multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val );

	//********************* MUL *********************//
	template <int pR, int p1, int p2>
	__host__ __device__ friend void    certifMulExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 );
	template <int pR, int p1>
	__host__ __device__ friend void  certifMulExpans_d( multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val );
  template <int pR, int p1, int p2>
  __host__ __device__ friend void     truncMulExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 );
  template <int pR, int p1>
  __host__ __device__ friend void   truncMulExpans_d( multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val );
  template <int pR, int p1, int p2>
	__host__ __device__ friend void   QD_LikeMulExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 );
  template <int pR, int p1>
	__host__ __device__ friend void QD_LikeMulExpans_d( multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val );

	//********************* DIV *********************//
	template <int pR, int p1>
	__host__ __device__ friend void   invExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1 );
	template <int pR>
	__host__ __device__ friend void invExpans_d( multi_prec<pR> &res, const float val );
	template <int pR, int p1, int p2>
	__host__ __device__ friend void   divExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 );
	template <int pR, int p1>
	__host__ __device__ friend void divExpans_d( multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val );
  template <int pR, int p1>
  __host__ __device__ friend void divExpans_d( multi_prec<pR> &res, const float val, const multi_prec<p1> &mp1 );

	//********************* SQRT *********************//
  template <int pR, int p1>
  __host__ __device__ friend void      invSqrtExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1 );
  template <int pR>
  __host__ __device__ friend void    invSqrtExpans_d( multi_prec<pR> &res, const float val );
	template <int pR, int p1>
	__host__ __device__ friend void   sqrtNewtonExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1 );
	template <int pR>
  __host__ __device__ friend void sqrtNewtonExpans_d( multi_prec<pR> &res, const float val );
	template <int pR, int p1>
  __host__ __device__ friend void    sqrtHeronExpans( multi_prec<pR> &res, const multi_prec<p1> &mp1 );
	template <int pR>
  __host__ __device__ friend void  sqrtHeronExpans_d( multi_prec<pR> &res, const float val );
};

/*definitions*/
//********************* Constructors *********************//
template <int prec>
__host__ __device__ multi_prec<prec>::multi_prec(const float x){
	data[0] = x;
  for(int i=1; i<prec; i++) data[i] = 0;
}
template <int prec>
__host__ __device__ multi_prec<prec>::multi_prec(const float *datap,const int precp){
  for(int i=0; i<min(prec, precp); i++)	data[i]=datap[i];
  for(int i=min(prec,precp); i<prec; i++) data[i]=0;
}
template <int prec>
__host__ __device__ multi_prec<prec>::multi_prec(const char *s){
  if (read(s, *this)) {
    printf("(qd_real::qd_real): INPUT ERROR.");
    *this = 0.0;
  }
}
template <int prec>
template<int precP>
__host__ __device__ multi_prec<prec>::multi_prec(const multi_prec<precP> *mp){
	for(int i=0; i<min(prec, precP); i++)	data[i] = mp->getData()[i];
  for(int i=min(prec,precP); i<prec; i++) data[i]=0;
}
template <int prec>
__host__ __device__ multi_prec<prec>::multi_prec(const multi_prec<prec> *mp){
	for(int i=0; i<prec; i++)	data[i] = mp->getData()[i];
}

//********************* Getters & Setters *********************//
template <int prec>
__host__ __device__ int multi_prec<prec>::getPrec() const{ return prec; }
template <int prec>
__host__ __device__ const float* multi_prec<prec>::getData()const{ return data; }
template <int prec>
__host__ __device__ void multi_prec<prec>::setData(const float *datap, const int precp){  
  for(int i=0; i<min(prec, precp); i++)	data[i]=datap[i];
  for(int i=min(prec, precp); i<prec; i++) data[i]=0;
}
template <int prec>
__host__ __device__ void multi_prec<prec>::setData(const float *datap){  
  for(int i=0; i<prec; i++)	data[i] = datap[i];
}
template <int prec>
__host__ __device__ void multi_prec<prec>::setElement(const float datap, const int index){
	data[index] = datap;
}

//********************* Addition *********************//
/** Computes mp1 + mp2 and puts the result in res (an ulp-nonoverlapping expansion). **/
/** merge + renormalize **/
template <int pR, int p1, int p2>
__host__ __device__ void certifAddExpans(multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2){
  certifiedAdd<p1,p2,pR>(mp1.data, mp2.data, res.data);
}
template <int pR, int p1>
__host__ __device__ void certifAddExpans_d(multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val){
  renorm2L_4Add1<p1,pR>(mp1.data, val, res.data);
}

template <int pR, int p1, int p2>
__host__ __device__ void truncAddExpans(multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2){
  certifiedAdd<(pR<p1)?pR:p1,(pR<p2)?pR:p2,pR>(mp1.data, mp2.data, res.data);
}
template <int pR, int p1>
__host__ __device__ void truncAddExpans_d(multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val){
  renorm2L_4Add1<p1,pR>(mp1.data, val, res.data);
}

/** Computes mp1 + mp2 and puts the result in res (an ulp-nonoverlapping expansion). **/
/** QD generalized addition + renormalize_random **/
template <int pR, int p1, int p2>
__host__ __device__ void QD_LikeAddExpans(multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2){
  const int p = ((p1>=p2)) ? p1 : p2;
	if(p==p1) baileyAdd_renorm<p1,p2,pR>(mp1.data, mp2.data, res.data);
  else baileyAdd_renorm<p2,p1,pR>(mp2.data, mp1.data, res.data);
}
template <int pR, int p1>
__host__ __device__ void QD_LikeAddExpans_d(multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val){
  renorm2L_4Add1<p1,pR>(mp1.data, val, res.data);
}

//********************* Multiplication *********************//
/** Computes mp1 * mp2 and puts the result in res (an ulp-nonoverlapping expansion). **/
/** patial products accumulated in bins **/
template <int pR, int p1, int p2>
__host__ __device__ void certifMulExpans(multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2){
	certifiedMul<p1,p2,pR>(mp1.data, mp2.data, res.data);
}
template <int pR, int p1>
__host__ __device__ void certifMulExpans_d(multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val){
	certifiedMul<p1,1,pR>(mp1.data, &val, res.data);
}

/** Computes mp1 * mp2 and puts the result in res (an ulp-nonoverlapping expansion). **/
/** truncated patial products accumulated in bins **/
template <int pR, int p1, int p2>
__host__ __device__ void truncMulExpans(multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2){
  const int p = ((p1>=p2)) ? p1 : p2;
  if(p==p1) truncatedMul<p1,p2,pR>(mp1.data, mp2.data, res.data);
  else truncatedMul<p2,p1,pR>(mp2.data, mp1.data, res.data);
}
template <int pR, int p1>
__host__ __device__ void truncMulExpans_d(multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val){
  truncatedMul<p1,1,pR>(mp1.data, &val, res.data);
}

/** Computes mp1 * mp2 and puts the result in res (an ulp-nonoverlapping expansion). **/
/** QD generalized multiplication + renormalize_random **/
template <int pR, int p1, int p2>
__host__ __device__ void QD_LikeMulExpans(multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2){
  const int p = ((p1>=p2)) ? p1 : p2;
	if(p==p1) baileyMul_renorm<p1,p2,pR>(mp1.data, mp2.data, res.data);
  else baileyMul_renorm<p2,p1,pR>(mp2.data, mp1.data, res.data);
}
template <int pR, int p1>
__host__ __device__ void QD_LikeMulExpans_d(multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val){
  baileyMul_renorm<p1,1,pR>(mp1.data, &val, res.data);
}

//********************* Division *********************//
/** Computes the reciprocal of mp1 and puts the result in res (an ulp-nonoverlapping expansion). **/
/** newton raphson based algorithms **/
template <int pR, int p1>
__host__ __device__ void invExpans(multi_prec<pR> &res, const multi_prec<p1> &mp1){
	invNewtonUp<p1,pR>(mp1.data, res.data);
}
template <int pR>
__host__ __device__ void invExpans_d(multi_prec<pR> &res, const float val){
	invNewtonUp<1,pR>(&val, res.data);
}

template <int pR, int p1, int p2>
__host__ __device__ void divExpans(multi_prec<pR> &res, const multi_prec<p1> &mp1, const multi_prec<p2> &mp2){
  divNewton<p1,p2,pR>(mp1.data, mp2.data, res.data);
}
template <int pR, int p1>
__host__ __device__ void divExpans_d(multi_prec<pR> &res, const multi_prec<p1> &mp1, const float val){
  divNewton<p1,1,pR>(mp1.data, &val, res.data);
}
template <int pR, int p1>
__host__ __device__ void divExpans_d(multi_prec<pR> &res, const float val, const multi_prec<p1> &mp1){
  divNewton<1,p1,pR>(&val, mp1.data, res.data);
}

//********************* Square Root *********************//
/** Computes the reciprocal of sqrt(mp1) and puts the result in res (an ulp-nonoverlapping expansion). **/
/** newton raphson based algorithms **/
template <int pR, int p1>
__host__ __device__ void invSqrtExpans(multi_prec<pR> &res, const multi_prec<p1> &mp1){
  sqrtInvNewton<p1,pR>(mp1.data, res.data);
}
template <int pR>
__host__ __device__ void invSqrtExpans_d(multi_prec<pR> &res, const float val){
  sqrtInvNewton<1,pR>(&val, res.data);
}

template <int pR, int p1>
__host__ __device__ void sqrtNewtonExpans(multi_prec<pR> &res, const multi_prec<p1> &mp1){
  sqrtNewton<p1,pR>(mp1.data, res.data);
}
template <int pR>
__host__ __device__ void sqrtNewtonExpans_d(multi_prec<pR> &res, const float val){
	sqrtNewton<1,pR>(&val, res.data);
}

/** Computes the reciprocal of sqrt(mp1) and puts the result in res (an ulp-nonoverlapping expansion). **/
/** newton raphson based -- Heron iteration **/
template <int pR, int p1>
__host__ __device__ void sqrtHeronExpans(multi_prec<pR> &res, const multi_prec<p1> &mp1){
  sqrtHeron<p1,pR>(mp1.data, res.data);
}
template <int pR>
__host__ __device__ void sqrtHeronExpans_d(multi_prec<pR> &res, const float val){
  sqrtHeron<1,pR>(&val, res.data);
}

/**-----------operators-------------**/
/** Checks equality between to multiprecs.
    If precs are different, the higher precision terms have to be zero **/
template <int p>
template<int pS>
__host__ __device__ bool multi_prec<p>::operator ==(const multi_prec<pS> &mp2) const{
  for(int i=0; i<min(p,pS); i++)
    if(data[i] != (mp2.getData())[i]) return false;
  if(p<pS)
    for(int i=p;i<pS;i++) if(0 != (mp2.getData())[i]) return false;  
  if(pS<p)
    for(int i=pS;i<p;i++) if(0 != data[i]) return false;  
  return true;
}
template <int p>
__host__ __device__ bool multi_prec<p>::operator ==(const multi_prec<p> &mp2) const{
	for(int i=0; i<p; i++)
    if(data[i] != (mp2.getData())[i]) return false;  
  return true;
}
template <int p>
__host__ __device__ bool multi_prec<p>::operator ==(const float &mp2) const{
  if (p==1 && data[0] == mp2) return true;
  else if (data[0] != mp2) return false;
  else if (p>1 && data[0] == mp2)
    for(int i=1; i<p; i++) if(data[i] != 0) return false;
  return true;
}

/** Checks inequality between to multiprecs.
    Works as NOT (operator ==) defined above **/
template <int p>
template <int pS>
__host__ __device__ bool multi_prec<p>::operator !=(const multi_prec<pS> &mp2) const{ return !((*this) == mp2); }
template <int p>
__host__ __device__ bool multi_prec<p>::operator !=(const multi_prec<p> &mp2) const{ return !((*this) == mp2); }
template <int p>
__host__ __device__ bool multi_prec<p>::operator !=(const float &mp2) const{ return !((*this) == mp2); }

/** Checks equality between to multiprecs.
    If precs are different, the higher precision terms have to be zero **/
template <int p>
template<int pS>
__host__ __device__ bool multi_prec<p>::operator <(const multi_prec<pS> &mp2) const{
  int i;
  for(i=0; i<min(p,pS); i++)
    if(data[i] > mp2.getData()[i]) return false;
    else if(data[i] < mp2.getData()[i]) return true;

  if(p < pS)
    for(; i<pS; i++){
      if(mp2.getData()[i] != 0.)    
        if ((mp2.getData())[i] < 0) return false;
        else return true;
    }
  else if(p > pS)
    for(; i<p; i++){
      if(data[i] < 0) return true;
      else return false;
    }
  return false;
}
template <int p>
__host__ __device__ bool multi_prec<p>::operator <(const multi_prec<p> &mp2) const{
  for(int i=0; i<p; i++)
    if (data[i] > mp2.getData()[i]) return false;
    else if (data[i] < mp2.getData()[i]) return true;
  return false;
}
template <int p>
__host__ __device__ bool multi_prec<p>::operator <(const float &mp2) const{
  if(data[0] < mp2) return true;
  else if (data[0] == mp2)
    if (data[1] < 0.) return true;
    else if (data[1] > 0.) return false;
         else
           for(int i=2; i<p; i++)
             if (data[i] != 0.) return false;
  return false;
}

/** Checks equality between to multiprecs.
    If precs are different, the higher precision terms have to be zero **/
template <int p>
template<int pS>
__host__ __device__ bool multi_prec<p>::operator <=(const multi_prec<pS> &mp2) const{
  return ((*this) < mp2) || ((*this) == mp2);
}
template <int p>
__host__ __device__ bool multi_prec<p>::operator <=(const multi_prec<p> &mp2) const{
  return ((*this) < mp2) || ((*this) == mp2);
}
template <int p>
__host__ __device__ bool multi_prec<p>::operator <=(const float &mp2) const{
  return ((*this) < mp2) || ((*this) == mp2);
}

/** Checks equality between to multiprecs.
    If precs are different, the higher precision terms have to be zero **/
template <int p>
template<int pS>
__host__ __device__ bool multi_prec<p>::operator >(const multi_prec<pS> &mp2) const{
  return !((*this) < mp2) && ((*this) != mp2);
}
template <int p>
__host__ __device__ bool multi_prec<p>::operator >(const multi_prec<p> &mp2) const{
  return !((*this) < mp2) && ((*this) != mp2);
}
template <int p>
__host__ __device__ bool multi_prec<p>::operator >(const float &mp2) const{
  return !((*this) < mp2) && ((*this) != mp2);
}

template <int p>
template<int pS>
__host__ __device__ bool multi_prec<p>::operator >=(const multi_prec<pS> &mp2) const{
  return ((*this) > mp2) || ((*this) == mp2);
}
template <int p>
__host__ __device__ bool multi_prec<p>::operator >=(const multi_prec<p> &mp2) const{
  return ((*this) > mp2) || ((*this) == mp2);
}
template <int p>
__host__ __device__ bool multi_prec<p>::operator >=(const float &mp2) const{
  return ((*this) > mp2) || ((*this) == mp2);
}

template <int pR, int p1, int p2>
__host__ __device__ multi_prec<pR> max( const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 ){
  if (mp1 < mp2) return mp2;
  else return mp1;
}
template <int pR, int p1, int p2>
__host__ __device__ multi_prec<pR> min( const multi_prec<p1> &mp1, const multi_prec<p2> &mp2 ){
  if (mp1 > mp2) return mp2;
  else return mp1;
}

/** Returns the absolute value of the multi_prec parameter mp
    Remark: if data[0] < 0 returns the negation of mp **/
template <int p>
__host__ __device__ multi_prec<p> abs(const multi_prec<p> &mp){
 	multi_prec<p> res;
 	res = mp;
	if((res.getData())[0] < 0) return -res;
	return res;
}

template <int prec>
multi_prec<prec> sqrt(const multi_prec<prec> &mp){
	multi_prec<prec> res;
	sqrtNewtonExpans(res, mp);
	return res;
}
template <int prec>
multi_prec<prec> invSqrt(const multi_prec<prec> &mp){
	multi_prec<prec> res;
	invSqrtExpans(res, mp);
	return res;
}

template <int p>
__host__ __device__ void renorm(multi_prec<p> &mp){ renorm2L<p, p>(mp.data); }
template <int p>
__host__ __device__ void renorm_rand(multi_prec<p> &mp){ renorm_rand2L<p, p>(mp.data); }
template <int p>
__host__ __device__ void renorm_2ndL(multi_prec<p> &mp){ fast_VecSumErrBranch<p, p>(mp.data); }

template<int pS>
__host__ __device__ int read(const char *s, multi_prec<pS> &mp) {
  const char *p = s;
  char ch;
  int sign = 0;
  int point = -1;  /* location of decimal point */
  int nd = 0;      /* number of digits read */
  int e = 0;       /* exponent. */
  bool done = false;
  multi_prec<pS> r = 0.0;  /* number being read */

  /* Skip any leading spaces */
  while (*p == ' ') p++;
  while (!done && (ch = *p) != '\0') {
    if (ch >= '0' && ch <= '9') { /* It's a digit */
      int d = ch - '0';
      r *= 10.0;
      r += static_cast<float>(d);
      nd++;
    } else { /* Non-digit */
      switch (ch) {
      case '.':
        if (point >= 0) return -1;   /* we've already encountered a decimal point. */
        point = nd;
        break;
      case '-':
      case '+':
        if (sign != 0 || nd > 0) return -1;  /* we've already encountered a sign, or if its
                                                not at first position. */
        sign = (ch == '-') ? -1 : 1;
        break;
      case 'E':
      case 'e':
        int nread;
        nread = std::sscanf(p+1, "%d", &e);
        done = true;
        if (nread != 1) return -1;  /* read of exponent failed. */
        break;
      case ' ':
        done = true;
        break;
      default:
        return -1;
      }
    }
    p++;
  }

  /* Adjust exponent to account for decimal point */
  if (point >= 0) e -= (nd - point);

  /* Multiply the the exponent */
  if (e != 0) r *= pow(10.0,e); //(qd_real(10.0) ^ e);

  mp = (sign < 0) ? -r : r;
  return 0;
}

#endif
