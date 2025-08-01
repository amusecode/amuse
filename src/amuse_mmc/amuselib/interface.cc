#include "worker_code.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef __cplusplus
extern "C" {
#endif
int asl_add(double x, double y, double *z);
int asl_many_points_on_sphere(double *x_req, double *y_req, double *z_req, int N);
int asl_rnd_points_on_sphere(double *x_req, double *y_req, double *z_req, int N);
#ifdef __cplusplus
}
#endif

int add(double x, double y, double *z)
{
  asl_add(x, y, z);
  return 0;
}

int many_points_on_sphere(double *x, double *y, double *z, int len)
{
  asl_many_points_on_sphere(x, y, z, len);
  return 11;
}

int rnd_points_on_sphere(double *x, double *y, double *z, int len)
{
  fprintf(stdout, "start calculating points\n");
  fflush(stdout);
  asl_rnd_points_on_sphere(x, y, z, len);
  fprintf(stdout, "done calculating points\n");
  fflush(stdout);
  return 11;
}
