/*
  AMUSE interface to SEI
  ======================

  Hanno Rein, Scott Tremaine

*/
#include <stdio.h>
#include <math.h>
#ifdef __cplusplus
extern "C" {
#endif
#include "src/sei.h"
#ifdef __cplusplus
}
#endif
#include "interface.h"

int initialization() {
  time = 0.0;
  dt = 1e-3 * 2.0 * M_PI;
  return 0;
}

int new_particle(int *id, double x, double y, double z,
                 double vx, double vy, double vz) 
{
  *id = 0;

  p.x = x;
  p.y = y;
  p.z = z;
  p.vx = vx;
  p.vy = vy;
  p.vz = vz;

  return 0;
}
int delete_particle(int id) {
  id = 0;

  return 0;
}

int set_state(int id, double x, double y, double z, 
	      double vx, double vy, double vz) 
{
  p.x = x;
  p.y = y;
  p.z = z;
  p.vx = vx;
  p.vy = vy;
  p.vz = vz;
  return 0;
}

int get_state(int id, double *x, double *y, double *z, 
	      double *vx, double *vy, double *vz) 
{
  *x = p.x;
  *y = p.y;
  *z = p.z;
  *vx = p.vx;
  *vy = p.vy;
  *vz = p.vz;
  return 0;
}

int evolve(double time_end) 
{
  while (time<time_end) {
    operator_sei(dt, &p);
    time += dt;
    //printf("New time = %f\n", time);
  }

  return 0;
}

