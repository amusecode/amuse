/*
 * The Kepler solver evolves the two-body problem, using standard Huayno sys structures.
 * TODO add warning when using the kepler solver with softening?
 */

#include <stdio.h>
#include <stdlib.h>
#include "evolve.h"
#include "evolve_shared.h"
#include "universal_variable_kepler.h"

void evolve_kepler(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt) {
  if (etime == stime ||  dt==0 || clevel>=MAXLEVEL) ENDRUN("timestep too small\n");
  if (s.n != 2) ENDRUN("two-body solver was called with sys.n=%u\n", s.n);
  // translate coordinates original frame to 2-body frame
  DOUBLE dpos[3],dpos0[3],pos_cm[3];
  DOUBLE dvel[3],dvel0[3],vel_cm[3];
  DOUBLE m1 = s.part->mass;
  DOUBLE m2 = s.last->mass;
  DOUBLE mtot = s.part->mass + s.last->mass;
  DOUBLE f1 = m2 / mtot;
  DOUBLE f2 = m1 / mtot;
  for(int k=0;k<3;k++) dpos0[k] = s.part->pos[k] - s.last->pos[k];
  for(int k=0;k<3;k++) dvel0[k] = s.part->vel[k] - s.last->vel[k];  
  for(int k=0;k<3;k++) pos_cm[k] = (m1 * s.part->pos[k] + m2 * s.last->pos[k]) / mtot;    
  for(int k=0;k<3;k++) vel_cm[k] = (m1 * s.part->vel[k] + m2 * s.last->vel[k]) / mtot;
  // evolve center of mass for dt
  for(int k=0;k<3;k++) pos_cm[k] += vel_cm[k] * dt;
  // call kepler solver
  int err = universal_variable_kepler_solver(dt,mtot,dpos0,dvel0,dpos,dvel);
  if (err != 0) {
    // failsafe kepler solver
    LOG("Kepler solver failed\n");
    diag->cefail[clevel]++;
    evolve_shared4(clevel,s, stime, etime, dt, 1);
  } else {
    // translate coordinates from 2-body frame to original frame
    for(int k=0;k<3;k++) s.part->pos[k] = pos_cm[k] + f1 * dpos[k];
    for(int k=0;k<3;k++) s.part->vel[k] = vel_cm[k] + f1 * dvel[k];
    for(int k=0;k<3;k++) s.last->pos[k] = pos_cm[k] - f2 * dpos[k];
    for(int k=0;k<3;k++) s.last->vel[k] = vel_cm[k] - f2 * dvel[k];
    // update statistics
  }
  diag->cecount[clevel]++;
}
