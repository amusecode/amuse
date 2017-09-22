#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdio.h>

#define CONFIG_USE_DOUBLE

typedef double REAL;
typedef long INT;

#define PI ((REAL)(3.141592653589793))

#ifdef CONFIG_USE_DOUBLE
    #define TOLERANCE ((REAL)(2.2737367544323205948e-13))     // 2^(-42)
#else
    #define TOLERANCE ((REAL)(1.52587890625e-5))              // (2^-16)
#endif
#define MAXITER 64
#define COMPARE(x, y) (((x) > (y)) - ((x) < (y)))
#define SIGN(x) COMPARE(x, 0)

#endif // __COMMON_H__
