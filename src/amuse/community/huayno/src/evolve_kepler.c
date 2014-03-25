/*
 * The Kepler solver evolves the two-body problem, using standard Huayno sys structures.
 * TODO add warning when using the kepler solver with softening?
 */

#include <stdio.h>
#include <stdlib.h>
#include "evolve.h"
#include "evolve_shared.h"
#include "universal_kepler_solver.h"

void evolve_kepler(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt) {
  CHECK_TIMESTEP(etime,stime,dt,clevel);
  if (s.n != 2) ENDRUN("two-body solver was called with sys.n=%u\n", s.n);
  // translate coordinates original frame to 2-body frame
  int k;
  DOUBLE dpos[3],dpos0[3],pos_cm[3];
  DOUBLE dvel[3],dvel0[3],vel_cm[3];
  DOUBLE m1 = s.part->mass;
  DOUBLE m2 = s.last->mass;
  DOUBLE mtot = s.part->mass + s.last->mass;
  DOUBLE f1 = m2 / mtot;
  DOUBLE f2 = m1 / mtot;
  if(mtot>0.) {
    for(k=0;k<3;k++) dpos0[k] = s.part->pos[k] - s.last->pos[k];
    for(k=0;k<3;k++) dvel0[k] = s.part->vel[k] - s.last->vel[k];
    for(k=0;k<3;k++) pos_cm[k] = (m1 * s.part->pos[k] + m2 * s.last->pos[k]) / mtot;    
    for(k=0;k<3;k++) vel_cm[k] = (m1 * s.part->vel[k] + m2 * s.last->vel[k]) / mtot;
    // evolve center of mass for dt
    for(k=0;k<3;k++) pos_cm[k] += vel_cm[k] * dt;
    // call kepler solver
    int err=universal_kepler_solver(dt,mtot,eps2,
                                    dpos0[0],dpos0[1],dpos0[2],
                                    dvel0[0],dvel0[1],dvel0[2],
                                    &dpos[0],&dpos[1],&dpos[2],
                                    &dvel[0],&dvel[1],&dvel[2]);
    if (err != 0) ENDRUN("kepler solver failure"); // failure of the kepler solver should be very rare now
    // translate coordinates from 2-body frame to original frame
    for(k=0;k<3;k++) s.part->pos[k] = pos_cm[k] + f1 * dpos[k];
    for(k=0;k<3;k++) s.part->vel[k] = vel_cm[k] + f1 * dvel[k];
    for(k=0;k<3;k++) s.last->pos[k] = pos_cm[k] - f2 * dpos[k];
    for(k=0;k<3;k++) s.last->vel[k] = vel_cm[k] - f2 * dvel[k];
  } else {
    for(k=0;k<3;k++) s.part->pos[k]+=s.part->vel[k]*dt;
    for(k=0;k<3;k++) s.last->pos[k]+=s.last->vel[k]*dt;
  } 
  diag->cecount[clevel]++;
}

// special solver to test kepler
void evolve_kepler_test(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt) {
  struct particle *central=s.part;
  struct particle p1,p2;
  struct sys s2;

  CHECK_TIMESTEP(etime,stime,dt,clevel);

  for(struct particle *i=s.part+1;i<=s.last;i++) if(central->mass<i->mass) central=i;

  for(struct particle *i=s.part;i<=s.last;i++)
  {
    if(i==central) continue;
    p1=*central;
    p2=*i;

    p1.mass=central->mass+i->mass;
    p2.mass=0;

    s2.n=2;
    s2.part=&p1;
    s2.last=&p2;
    evolve_kepler(clevel,s2,stime,etime,dt);
    p2.mass=i->mass;
    *i=p2;
  }
  for(int k=0;k<3;k++) central->pos[k]+=central->vel[k]*dt;
}
