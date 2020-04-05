/*
 * genRenorm.cpp
 * This file is part of MultiPrecGPU Library
 *
 * Copyright (C) 2012 - 
 *
 * MultiPrecGPU Library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * MultiPrecGPU Library is distributed in the hope that it will be useful,
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
/*
 * Printing the renormalization algorithm in an unrolled form
 */

#include <stdio.h>
#include <stdlib.h>
#include <fenv.h>
#include <cstdio>
#include <iostream>
#include <fstream>
using namespace std;

static void print_fast_renorm2L_inPlace(int const sX, int const sR){
	std::ofstream renormFct;
  renormFct.open("../Specs/specRenorm.h",ios::app);

	renormFct << "template<>\n";
  renormFct << "__host__ __device__ __forceinline__ void fast_renorm2L<" << sX << "," << sR << ">(float *x){\n";

	renormFct << "\tint ptr=0;\n";

	renormFct << "\tfloat pr = fast_two_sum(x[" << sX-2 << "], x[" << sX-1 << "], x[" << sX-1 << "]);\n";
	for(int i=sX-3; i>0; i--)	
		renormFct << "\tpr = fast_two_sum(x[" << i << "], pr, x[" << i+1 << "]);\n";
	renormFct << "\tx[0] = fast_two_sum(x[0], pr, x[1]);\n\n";
	
	renormFct << "\tif(x[1] == 0.) pr = x[0];\n";
	renormFct << "\telse { pr = x[1]; ptr++; }\n";
	
	for(int i=2; i<sX-1; i++){
		renormFct << "\tx[ptr] = fast_two_sum(pr, x[" << i << "], pr);\n";
		renormFct << "\tif(pr == 0.) pr = x[ptr]; else ptr++;\n";
	}
	renormFct << "\tx[ptr] = fast_two_sum(pr, x[" << sX-1 << "], x[ptr+1]);\n\n";
	
	renormFct << "\tfor(int i=ptr+1; i<" << sR << "; i++) x[i] = 0.;\n";

  renormFct << "}\n\n";
	renormFct.close();
}

static void print_fast_renorm2L(int const sX, int const sR){
	std::ofstream renormFct;
  renormFct.open("../Specs/specRenorm.h",ios::app);

	renormFct << "template<>\n";
  renormFct << "__host__ __device__ __forceinline__ void fast_renorm2L<" << sX << "," << sR << ">(const float *x, float *r){\n";

	renormFct << "\tfloat f[" << sX << "];\n\tint ptr=0;\n";

	renormFct << "\tfloat pr = fast_two_sum(x[" << sX-2 << "], x[" << sX-1 << "], f[" << sX-1 << "]);\n";
	for(int i=sX-3; i>0; i--)	
		renormFct << "\tpr = fast_two_sum(x[" << i << "], pr, f[" << i+1 << "]);\n";
	renormFct << "\tf[0] = fast_two_sum(x[0], pr, f[1]);\n\n";
	
	renormFct << "\tif(f[1] == 0.) pr = f[0];\n";
	renormFct << "\telse { r[0] = f[0]; pr = f[1]; ptr++; }\n";
	
	for(int i=2; i<sR; i++){
		renormFct << "\tr[ptr] = fast_two_sum(pr, f[" << i << "], pr);\n";
		renormFct << "\tif(pr == 0.) pr = r[ptr]; else ptr++;\n";
	}
	
	renormFct << "\n\tfor(int i=" << sR << "; ptr<" << sR << " && i<" << sX << "; i++){\n";
	renormFct << "\t\tr[ptr] = fast_two_sum(pr, f[i], pr);\n";
	renormFct << "\t\tif(pr == 0.) pr = r[ptr]; else ptr++;\n";
	renormFct << "\t}\n\n";
	
	renormFct << "\tif(ptr<" << sR << " && pr!=0.){ r[ptr] = pr; ptr++; }\n";
	renormFct << "\tfor(int i=ptr; i<" << sR << "; i++) r[i] = 0.;\n";
	
  renormFct << "}\n\n";
	renormFct.close();
}

static void print_renorm2L_4Add1(int const sX, int const sR){
	std::ofstream renormFct;
  renormFct.open("../Specs/specRenorm.h",ios::app);

	renormFct << "template<>\n";
  renormFct << "__host__ __device__ __forceinline__ void renorm2L_4Add1<" << sX << "," << sR << ">(const float *x, const float y, float *r){\n";

	renormFct << "\tfloat f[" << sX+1 << "];\n\tint ptr=0;\n";

	renormFct << "\tfloat pr = two_sum(x[" << sX-1 << "], y, f[" << sX << "]);\n";
	for(int i=sX-2; i>0; i--)	
		renormFct << "\tpr = two_sum(x[" << i << "], pr, f[" << i+1 << "]);\n";
	renormFct << "\tf[0] = two_sum(x[0], pr, f[1]);\n\n";
	
	renormFct << "\tif(f[1] == 0.) pr = f[0];\n";
	renormFct << "\telse { r[0] = f[0]; pr = f[1]; ptr++; }\n";
	
	for(int i=2; i<sR; i++){
		renormFct << "\tr[ptr] = fast_two_sum(pr, f[" << i << "], pr);\n";
		renormFct << "\tif(pr == 0.) pr = r[ptr]; else ptr++;\n";
	}
	
	renormFct << "\n\tfor(int i=" << sR << "; ptr<" << sR << " && i<" << sX+1 << "; i++){\n";
	renormFct << "\t\tr[ptr] = fast_two_sum(pr, f[i], pr);\n";
	renormFct << "\t\tif(pr == 0.) pr = r[ptr]; else ptr++;\n";
	renormFct << "\t}\n\n";
	
	renormFct << "\tif(ptr<" << sR << " && pr!=0.){ r[ptr] = pr; ptr++; }\n";
	renormFct << "\tfor(int i=ptr; i<" << sR << "; i++) r[i] = 0.;\n";
	
  renormFct << "}\n\n";
	renormFct.close();
}

static void print_baileyAdd(int const K, int const L, int const R){
  if(L==1){
    std::ofstream addFct;
    addFct.open("../Specs/specAddition.h",ios::app);

	  addFct << "template<>\n";
    addFct << "__host__ __device__ __forceinline__ void baileyAdd_fast<" << K << "," << L << "," << R << ">(float const *x, float const *y, float *z){\n";
    addFct << "\trenorm2L_4Add1<" << K << "," << R << ">( x, y[0], z );\n";
  
    addFct << "}\n\n";
    addFct.close();
  } else {

	std::ofstream addFct;
  addFct.open("../Specs/specAddition.h",ios::app);

	addFct << "template<>\n";
  addFct << "__host__ __device__ __forceinline__ void baileyAdd_fast<" << K << "," << L << "," << R << ">(float const *x, float const *y, float *z){\n";
	addFct << "\tfloat e, f[" << R+1 << "];\n";

  int RR = R+1;
  int i, n;
  for(i=K; i<RR-1; i++)
    addFct << "\tf[" << i << "] = 0.0;\n";

  n = RR-1;
  if (n<L){
    addFct << "\tf[" << n << "] = FPadd_rn(x[" << n << "], y[" << n << "]);\n";
  }else if (n<K)
    addFct << "\tf[" << n << "] = x[" << n << "];\n";
  else
    addFct << "\tf[" << n << "] = 0.0;\n";

  for (n=(RR>K)?K-1:RR-2; n>=L; n--)
    addFct << "\tf[" << n << "] = x[" << n << "];\n";

  for (n=(RR>L)?L-1:RR-2; n>=0; n--){
    addFct << "\tf[" << n << "] = two_sum(x[" << n << "], y[" << n << "], e);\n";

    for(i=n+1; i<RR-1; i++)
      addFct << "\tf[" << i << "] = two_sum(f[" << i << "], e, e);\n";

    addFct << "\tf[" << RR-1 << "] = FPadd_rn(f[" << RR-1 << "], e);\n"; 
  }
  addFct << "\tfast_renorm2L<" << RR << "," << R << ">( f, z );\n";
  
  addFct << "}\n\n";

	addFct.close();}
}

static void print_baileyMul(int const K, int const L, int const R){
	std::ofstream mulFct;
  mulFct.open("../Specs/specMultiplication.h",ios::app);

	mulFct << "template<>\n";
  mulFct << "__host__ __device__ __forceinline__ void baileyMul_fast<" << K << "," << L << "," << R << ">(float const *x, float const *y, float *z){\n";
	if(L>1)	mulFct << "\tfloat p, e, f[" << R+1 << "];\n";
	else mulFct << "\tfloat e, f[" << R+1 << "];\n";

  int RR = R+1;
	int i, j, n;
	for(i=K+L-1; i<RR-1; i++)
		mulFct << "\tf[" << i << "] = 0.0;\n";

	n = RR-1;
	if (n<L){
		mulFct << "\tf[" << n << "] = FPmul_rn(x[0], y[" << n << "]);\n";
    mulFct << "\t#if defined __CUDA_ARCH__ || FP_FAST_FMA\n";
		for(i=1; i<=n; i++){
			mulFct << "\t\tf[" << n << "] = FPfma_rn(x[" << i << "], y[" << n-i << "], f[" << n << "]);\n";		
		}
		mulFct << "\t#else\n";
		for(i=1; i<=n; i++){
			mulFct << "\t\tf[" << n << "] = FPadd_rn(f[" << n << "], FPmul_rn(x[" << i << "], y[" << n-i << "]));\n";		
		}
		mulFct << "\t#endif\n";
	}else if (n<K){
		mulFct << "\tf[" << n << "] = FPmul_rn(x[" << n << "], y[0]);\n";
		mulFct << "\t#if defined __CUDA_ARCH__ || FP_FAST_FMA\n";
		for(j=1; j<L; j++){
			mulFct << "\t\tf[" << n << "] = FPfma_rn(x[" << n-j << "], y[" << j << "], f[" << n << "]);\n";		
		}
		mulFct << "\t#else\n";
		for(j=1; j<L; j++){
			mulFct << "\t\tf[" << n << "] = FPadd_rn(f[" << n << "], FPmul_rn(x[" << n-j << "], y[" << j << "]));\n";
		}
		mulFct << "\t#endif\n";
	}else if (n<K+L-1){
		mulFct << "\tf[" << n << "] = FPmul_rn(x[" << n-L+1 << "], y[" << L-1 << "]);\n";
		mulFct << "\t#if defined __CUDA_ARCH__ || FP_FAST_FMA\n";
		for(i=n-L+2; i<K; i++){
			mulFct << "\t\tf[" << n << "] = FPfma_rn(x[" << i << "], y[" << n-i << "], f[" << n << "]);\n";		
		}
		mulFct << "\t#else\n";
		for(i=n-L+2; i<K; i++){
			mulFct << "\t\tf[" << n << "] = FPadd_rn(f[" << n << "], FPmul_rn(x[" << i << "], y[" << n-i << "]));\n";	
		}
		mulFct << "\t#endif\n";
	}else
		mulFct << "\tf[" << n << "] = 0.0;\n";

	/***********************************************************/
	for (n=(RR>K+L-1)?K+L-2:RR-2; n>=K; n--){
		mulFct << "\tf[" << n << "] = two_prod(x[" << n-L+1 << "], y[" << L-1 << "], e);\n";

		for(j=n+1; j<RR-1; j++)
			mulFct << "\tf[" << j << "] = two_sum(f[" << j << "], e, e);\n";
		mulFct << "\tf[" << RR-1 << "] = FPadd_rn(f[" << RR-1 << "], e);\n";	

		for(i=n-L+2; i<K; i++){
			mulFct << "\tp = two_prod(x[" << i << "], y[" << n-i << "], e);\n";
			for(j=n+1; j<RR-1; j++)
				mulFct << "\tf[" << j << "] = two_sum(f[" << j << "], e, e);\n";
			mulFct << "\tf[" << RR-1	<< "] = FPadd_rn(f[" << RR-1 << "], e);\n";	

			mulFct << "\tf[" << n << "] = two_sum(f[" << n << "], p, e);\n";
			for(j=n+1; j<RR-1; j++)
				mulFct << "\tf[" << j << "] = two_sum(f[" << j << "], e, e);\n";
			mulFct << "\tf[" << RR-1 << "] = FPadd_rn(f[" << RR-1 << "], e);\n";	
		}
	}

	/***********************************************************/
	for (n=(RR>K)?K-1:RR-2; n>=L; n--){
		mulFct << "\tf[" << n << "] = two_prod(x[" << n << "], y[0], e);\n";

		for(i=n+1; i<RR-1; i++)
			mulFct << "\tf[" << i << "] = two_sum(f[" << i << "], e, e);\n";
		mulFct << "\tf[" << RR-1 << "] = FPadd_rn(f[" << RR-1 << "], e);\n";	

		for(j=1; j<L; j++){
			mulFct << "\tp = two_prod(x[" << n-j << "], y[" << j << "], e);\n";
			for(i=n+1; i<RR-1; i++)
				mulFct << "\tf[" << i << "] = two_sum(f[" << i << "], e, e);\n";
			mulFct << "\tf[" << RR-1 << "] = FPadd_rn(f[" << RR-1 << "], e);\n";	

			mulFct << "\tf[" << n << "] = two_sum(f[" << n << "], p, e);\n";
			for(i=n+1; i<RR-1; i++)
				mulFct << "\tf[" << i << "] = two_sum(f[" << i << "], e, e);\n";
			mulFct << "\tf[" << RR-1 << "] = FPadd_rn(f[" << RR-1 << "], e);\n";	
		}
	}

	/***********************************************************/
	for (n=(RR>L)?L-1:RR-2; n>=0; n--){
		mulFct << "\tf[" << n << "] = two_prod(x[0], y[" << n << "], e);\n";

		for(j=n+1; j<RR-1; j++)
			mulFct << "\tf[" << j << "] = two_sum(f[" << j << "], e, e);\n";
		mulFct << "\tf[" << RR-1 << "] = FPadd_rn(f[" << RR-1 << "], e);\n";	

		for(i=1; i<=n; i++){
			mulFct << "\tp = two_prod(x[" << i << "], y[" << n-i << "], e);\n";
			for(j=n+1; j<RR-1; j++)
				mulFct << "\tf[" << j << "] = two_sum(f[" << j << "], e, e);\n";
			mulFct << "\tf[" << RR-1 << "] = FPadd_rn(f[" << RR-1 << "], e);\n";	

			mulFct << "\tf[" << n << "] = two_sum(f[" << n << "], p, e);\n";
			for(j=n+1; j<RR-1; j++)
				mulFct << "\tf[" << j << "] = two_sum(f[" << j << "], e, e);\n";
			mulFct << "\tf[" << RR-1 << "] = FPadd_rn(f[" << RR-1 << "], e);\n";	
		}
	}
  mulFct << "\tfast_renorm2L<" << RR << "," << R << ">( f, z );\n";

	mulFct << "}\n\n";
	mulFct.close();
}

int main(int argc, char *argv[]) {
	if(argc<4){
		cout << "Usage: ./genRenorm K L R" << endl; 
   	return 0;
  }

	int K = atoi(argv[1]);
	int L = atoi(argv[2]);	
	int R = atoi(argv[3]);

  if (K==1) 
    print_renorm2L_4Add1(L, R);
  else if (L==1) 
    print_renorm2L_4Add1(K, R);
 	if (K>=L){
   	print_baileyAdd(K, L, R);
	  print_baileyMul(K, L, R);
  }else{
    print_baileyAdd(L, K, R);
	  print_baileyMul(L, K, R);
  }

	print_fast_renorm2L_inPlace(R+1, R);
	print_fast_renorm2L(R+1, R);

	return 0;
}

