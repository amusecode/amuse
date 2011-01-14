#include "octgrav.h"
#include <math.h>

extern int n_bodies;

extern double eps;   //softening
extern double theta;   //angle opening
extern double t_now;   //start time
extern double dtime; //time-step (hardcoded)
extern double t_end;   //end time

extern double dt_out;  //output interval
//int nstep = 0;
extern double E_init;

//void F4XYZ_ADDMUL(float4 &a, float4 &b, double c);

//void MASS_F4XYZ_ADDMUL(vector<float4> &a, vector<float4> &b, double c);

double calcEkin(vector<float4> &bodies_pos, vector<float4> &bodies_vel);

double calcEpot(vector<float4> &bodies_pos, vector<float4> &bodies_grav);

void write_output(vector<float4> &bodies_pos, vector<float4> &bodies_vel, vector<float4> &bodies_grav);

/* Performs one integration step using leapfrog scheme. */
void leapfrog(double dtime, vector<float4> &bodies_pos, vector<float4> &bodies_vel, vector<float4> &bodies_grav, octgrav &system);

