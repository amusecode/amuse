#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mameclot.h"

#define sqr(x) pow(x,2.0)

void calculate_potential(float *m, float *x, float *y, float *z, float *phi, int N, int N1);

void calculate_potential(float *m, float *x, float *y, float *z, float *phi, int N, int N1)
{
  // Calculates specific potential of a cluster (N1 = 0) or a cluster pair (0 < N1 < N)
  double d2;
  for (int i=0; i < (N1 == 0 ? N : N1); i++){
    for (int j = (N1 == 0 ? i + 1 : N1); j < N; j++){
      d2 = sqr(x[i]-x[j]) +  sqr(y[i]-y[j]) +  sqr(z[i]-z[j]);
      phi[i] -= m[j]/sqrt(d2);
      phi[j] -= m[i]/sqrt(d2);
    }   
  }
}


