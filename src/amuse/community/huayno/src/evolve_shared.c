/*
 * Reference integrators with single, global shared time step.
 */

#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#include "evolve.h"
#include "integrators_shared.h"

void evolve_shared2(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  FLOAT dtsys;
  if(etime == stime ||  dt==0 || clevel>=MAXLEVEL)
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel);
  if(calc_timestep) timestep(clevel,s,s,SIGN(dt));
  dtsys=global_timestep(s);
  if(dtsys < fabs(dt))
  {
    evolve_shared2(clevel+1,s,stime, stime+dt/2,dt/2,0);
    evolve_shared2(clevel+1,s,stime+dt/2, etime,dt/2,1);
  }
  else
  {
    diag->deepsteps++;
    diag->simtime+=dt;
    kdk(clevel,s,zerosys, stime, etime, dt);
  }
}

void evolve_shared4(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep) {
  FLOAT dtsys;
  if(etime == stime ||  dt==0 || clevel>=MAXLEVEL)
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel);
  if(calc_timestep) timestep(clevel,s,s,SIGN(dt));
  dtsys = global_timestep(s);
  if(dtsys < fabs(dt)) {
    evolve_shared4(clevel+1,s,stime, stime+dt/2,dt/2,0);
    evolve_shared4(clevel+1,s,stime+dt/2, etime,dt/2,1);
  } else {
    diag->deepsteps++;
    diag->simtime+=dt;
    dkd4(clevel,s, stime, etime, dt);
  }
}

void evolve_shared6(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep) {
  FLOAT dtsys;
  if(etime == stime ||  dt==0 || clevel>=MAXLEVEL)
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel);
  if(calc_timestep) timestep(clevel,s,s,SIGN(dt));
  dtsys = global_timestep(s);
  if(dtsys < fabs(dt)) {
    evolve_shared6(clevel+1,s,stime, stime+dt/2,dt/2,0);
    evolve_shared6(clevel+1,s,stime+dt/2, etime,dt/2,1);
  } else {
    diag->deepsteps++;
    diag->simtime+=dt;
    dkd6(clevel,s, stime, etime, dt);
  }
}

void evolve_shared8(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep) {
  FLOAT dtsys;
  if(etime == stime ||  dt==0 || clevel>=MAXLEVEL)
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel);
  if(calc_timestep) timestep(clevel,s,s,SIGN(dt));
  dtsys = global_timestep(s);
  if(dtsys < fabs(dt)) {
    evolve_shared8(clevel+1,s,stime, stime+dt/2,dt/2,0);
    evolve_shared8(clevel+1,s,stime+dt/2, etime,dt/2,1);
  } else {
    diag->deepsteps++;
    diag->simtime+=dt;
    dkd8(clevel,s, stime, etime, dt);
  }
}

void evolve_shared10(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep) {
  FLOAT dtsys;
  if(etime == stime ||  dt==0 || clevel>=MAXLEVEL)
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel);
  if(calc_timestep) timestep(clevel,s,s,SIGN(dt));
  dtsys = global_timestep(s);
  if(dtsys < fabs(dt)) {
    evolve_shared10(clevel+1,s,stime, stime+dt/2,dt/2,0);
    evolve_shared10(clevel+1,s,stime+dt/2, etime,dt/2,1);
  } else {
    diag->deepsteps++;
    diag->simtime+=dt;
    dkd10(clevel,s, stime, etime, dt);
  }
}
