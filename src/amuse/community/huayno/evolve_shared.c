/*
 * Reference integrators with single, global shared time step.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "evolve.h"

static FLOAT global_timestep(struct sys s);
static void dkd4(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);

void evolve_shared2(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  FLOAT dtsys;
  clevel++;
  if(etime <= stime ||  dt==0 || clevel>MAXLEVEL)
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel);
  if(calc_timestep) timestep(s,s);
  dtsys=global_timestep(s);
  if(dtsys < dt)
  {
    evolve_shared2(s,stime, stime+dt/2,dt/2,0);
    evolve_shared2(s,stime+dt/2, etime,dt/2,1);
  }
  else
  {
    deepsteps++;
    simtime+=dt;
    kdk(s,zerosys, stime, etime, dt);
  }
  clevel--;
}

void evolve_shared4(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep) {
  FLOAT dtsys;
  clevel++;
  if(etime <= stime ||  dt==0 || clevel>MAXLEVEL)
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel);
  if(calc_timestep) timestep(s,s);
  dtsys = global_timestep(s);
  if(dtsys < dt) {
    evolve_shared4(s,stime, stime+dt/2,dt/2,0);
    evolve_shared4(s,stime+dt/2, etime,dt/2,1);
  } else {
    deepsteps++;
    simtime+=dt;
    dkd4(s, stime, etime, dt);
  }
  clevel--;
}

static FLOAT global_timestep(struct sys s)
{
  UINT i;
  FLOAT mindt;
  mindt=HUGE_VAL;
  for(i=0;i<s.n;i++)
  {
    if(mindt>s.part[i].timestep) mindt=s.part[i].timestep;
  }
  return mindt;
}

static void dkd4(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  DOUBLE a1 = 0.0792036964311957;
  DOUBLE a2 = 0.353172906049774;
  DOUBLE a3 = -.0420650803577195;
  DOUBLE a4 = 1 - 2*(a1 + a2 + a3);

  DOUBLE b1 = 0.209515106613362;
  DOUBLE b2 = -.143851773179818;
  DOUBLE b3 = 0.5 - b1 - b2;

  DOUBLE dtime = 0;

  dtime += a1 * dt;
  drift(s, dtime, a1 * dt);
  dtime += b1 * dt;
  kick(s, s, b1 * dt);
  dtime += a2 * dt;
  drift(s, dtime, a2 * dt);
  dtime += b2 * dt;
  kick(s, s, b2 * dt);
  dtime += a3 * dt;
  drift(s, dtime, a3 * dt);
  dtime += b3 * dt;
  kick(s, s, b3 * dt);
  dtime += a4 * dt;
  drift(s, dtime, a4 * dt);
  dtime += b3 * dt;
  kick(s, s, b3 * dt);
  dtime += a3 * dt;
  drift(s, dtime, a3 * dt);
  dtime += b2 * dt;
  kick(s, s, b2 * dt);
  dtime += a2 * dt;
  drift(s, dtime, a2 * dt);
  dtime += b1 * dt;
  kick(s, s, b1 * dt);
  dtime += a1 * dt;
  drift(s, dtime, a1 * dt);
}
