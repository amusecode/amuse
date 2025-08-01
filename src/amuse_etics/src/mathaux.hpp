#pragma once
#include "common.hpp"
#include <cuComplex.h>
#ifdef ETICS_DOUBLE_PRECISION
    #define Complex      cuDoubleComplex
    #define make_Complex make_cuDoubleComplex
    #define Complex_add  cuCadd
    #define Complex_mul  cuCmul
    #define Complex_imag cuCimag
    #define Complex_real cuCreal
    #define Complex_conj cuConj
#else
    #define Complex      cuFloatComplex
    #define make_Complex make_cuFloatComplex
    #define Complex_add  cuCaddf
    #define Complex_mul  cuCmulf
    #define Complex_imag cuCimagf
    #define Complex_real cuCrealf
    #define Complex_conj cuConjf
#endif

#define SQRT_4_PI 3.5449077018110320545963349666822903655950989122447742564276


double FactorialSCF_tmpname(int x);
double GammaPlusHalf(int x);
void RadialCoefficients(Real *Result);
#define MAX_HARDCODED_PL 14
__host__ __device__ Real Pl(int l, Real x);
void AngularCoefficients(Real *Result);
