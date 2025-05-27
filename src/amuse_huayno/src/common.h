#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdio.h>
#include <tgmath.h>

#define REAL DOUBLE

#define PI ( (REAL) 3.14159265358979323846L)

#define TOLERANCE  (sizeof(DOUBLE)<8? 1.e-6:1.e-14)

#define MAXITER 64

#endif // __COMMON_H__
