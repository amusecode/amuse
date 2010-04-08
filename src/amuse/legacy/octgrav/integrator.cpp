#include "integrator.h"
#include <math.h>

int n_bodies = 0;

double eps   = 0.1;   //softening
double theta = 0.8;   //angle opening
double t_now = 0.0;   //start time
double dtime = 0.01; //time-step (hardcoded)
double t_end = 0.1;   //end time

double dt_out = 1.0;  //output interval
//int nstep = 0;
double E_init = 0.0;
//

/* Vector operations */
void F4XYZ_ADDMUL(float4 &a, float4 &b, double c) { 
  a.x += b.x * c;
  a.y += b.y * c;
  a.z += b.z * c;
}

int nstep = 0;

void MASS_F4XYZ_ADDMUL(vector<float4> &a, vector<float4> &b, double c) { 
  for(int i=0;i<n_bodies;i++) { 
    F4XYZ_ADDMUL(a[i],b[i],c);  
  }
}

/* Arithmetic */
#define PWR2(a) ((a)*(a))

/* Calculate kinetic energy */
double calcEkin(vector<float4> &bodies_pos, vector<float4> &bodies_vel) {
  double E_kin = 0.0;

  cerr << "integrator.cpp: n_bodies = " << n_bodies << ", bodies_pos[0] = " << bodies_pos[0].w << endl;
  // Kinetic energy
  for (int i = 0; i < n_bodies ; i++) {
    E_kin += 0.5 * bodies_pos[i].w * (PWR2(bodies_vel[i].x)+PWR2(bodies_vel[i].y)+PWR2(bodies_vel[i].z));

  }
  return E_kin;
}

/* Calculate potential energy (using potential 'w' values from particles) */
double calcEpot(vector<float4> &bodies_pos, vector<float4> &bodies_grav) 
{
  double rx,ry,rz;
  double r2;
  double E_pot = 0.0;
  for (int i = 0; i < n_bodies ; i++) {
    E_pot += 0.5 * bodies_pos[i].w * bodies_grav[i].w;
  }
  return E_pot;
}

void write_output(vector<float4> &bodies_pos, vector<float4> &bodies_vel, vector<float4> &bodies_grav)
{
  /* write output */
  double E_kin = calcEkin(bodies_pos,bodies_vel);
  double E_pot = calcEpot(bodies_pos,bodies_grav);
  
  double E_tot = E_kin + E_pot;
  double E_err_abs = E_tot - E_init;
  double E_err_rel = E_err_abs / E_init;

  //fprintf(stdout, "t_now: %lg  E_kin: %lg E_pot: %lg E_init: %lg\n", t_now, E_kin,E_pot,E_init);
  //fprintf(stdout, "t_now: %lg, t_end: %lg, dtime: %lg\n", t_now, t_end, dtime);

  //fprintf(stderr, "Absolute energy error: %lg\n",E_err_abs);
  //fprintf(stderr, "Relative energy error: %lg\n",E_err_rel);
  fflush(stdout);
  /* end write */
}

/******************************************************************************/
/*** EXPERIMENTAL variable but shared time-step leapfrog with time-symmetry ***/
/******************************************************************************/

#if 0
/*void init_goodfrog(vector<float4> &bodies_pos, vector<float4> &bodies_vel,vector<float4> &bodies_grav)
{
  h = calc_h(vector<float4> &bodies_pos, vector<float4> &bodies_vel);
  //empty?
}*/

vector<float4> bodies_grav_old;
double old_h, h;

void backup_grav(vector<float4> bodies_grav)
{  bodies_grav_old = bodies_grav; }

void calc_h(vector<float4> &bodies_pos, vector<float4> &bodies_vel)
{
  //todo...
}

void goodfrog(float dtime, vector<float4> &bodies_pos, vector<float4> &bodies_vel, vector<float4> &bodies_grav)
{

  // drift
  MASS_F4XYZ_ADDMUL(bodies_pos, bodies_vel, dtime);
  MASS_F4XYZ_ADDMUL(bodies_pos, bodies_grav, 0.5 * dtime * dtime); /* advance r by 1 step    */

  fprintf(stderr,"starting copy.\n");
  backup_grav(bodies_grav);
  //copy(bodies_grav.begin(), bodies_grav.end(), bodies_grav_old.begin());
  fprintf(stderr,"completing copy.\n");

  /*1. Calculate forces */
  //debug stuff
  octgrav system;
  system.set_softening(eps);
  system.set_opening_angle(theta);
  //end of debug stuff

  system.evaluate_gravity(bodies_pos, bodies_grav);                /* update grav by 1 step   */

  //this one can be moved up
  MASS_F4XYZ_ADDMUL(bodies_vel, bodies_grav_old, 0.5 * dtime);     /* update v by 1/2 step    */

  MASS_F4XYZ_ADDMUL(bodies_vel, bodies_grav, 0.5 * dtime);         /* update v by 1/2 step    */

  t_now = t_now + dtime;                                           /* advance time             */
  write_output(bodies_pos,bodies_vel,bodies_grav);                 /* write diagnostics        */

  nstep++;

  old_h = h;
  calc_h(bodies_pos, bodies_vel);
}
#endif

/******************************************************************************/
/*** CLASSIC LEAPFROG STEPPING                                              ***/
/******************************************************************************/

/* Performs one integration step using leapfrog scheme. */
void leapfrog(double dtime, vector<float4> &bodies_pos, vector<float4> &bodies_vel, vector<float4> &bodies_grav, octgrav &system)
{
 
  MASS_F4XYZ_ADDMUL(bodies_vel, bodies_grav, 0.5 * dtime);         /* advance vel by 1/2 step  */

  MASS_F4XYZ_ADDMUL(bodies_pos, bodies_vel, 0.5 * dtime);          /* advance pos by 1/2 step  */

  // loop over all bodies (MASS prefix implies a looping function)
  //fprintf(stderr,"particle 0: xvel = %f xpos = %f\n", bodies_vel[0].x, bodies_pos[0].x);
  
  MASS_F4XYZ_ADDMUL(bodies_pos, bodies_vel, 0.5 * dtime);          /* advance pos by 1/2 step  */

  //fprintf(stderr,"xvel = %f, xpos = %f, dtime = %f\n", bodies_vel[0].x, bodies_pos[0].x, dtime);

  //debug stuff (until Evghenii fixes this)
  // octgrav system;
  // system.set_softening(eps);
  // system.set_opening_angle(theta);
  //end of debug stuff

  system.evaluate_gravity(bodies_pos, bodies_grav);                /* perform force calc.      */
  //  calc_force_host(bodies_pos,bodies_grav);

  MASS_F4XYZ_ADDMUL(bodies_vel, bodies_grav, 0.5 * dtime);         /* advance vel by 1/2 step  */
  t_now = t_now + dtime;                                           /* advance time             */
  //cello write_output(bodies_pos,bodies_vel,bodies_grav);                 /* write diagnostics        */

  nstep++;                                    /* count another time step  */
 
}


