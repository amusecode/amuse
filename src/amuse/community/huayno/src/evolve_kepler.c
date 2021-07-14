/*
 * The Kepler solver evolves the two-body problem, using standard Huayno sys structures.
 * TODO add warning when using the kepler solver with softening?
 */

#include <stdio.h>
#include <stdlib.h>
#include "evolve.h"
#include "evolve_shared.h"
#include "universal_kepler_solver.h"
#ifdef _OPENMP
#include <omp.h>
#endif

static void evolve_kepler_2(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
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

static void evolve_kepler_n(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  struct particle *ipart;
  DOUBLE dpos[3],dpos0[3],cmpos[3];
  DOUBLE dvel[3],dvel0[3];
  UINT err;

  CHECK_TIMESTEP(etime,stime,dt,clevel);
  if (s.n-s.nzero > 1) ENDRUN("kepler-n solver was called with too many massive particles sys.n=%u\n", s.n-s.nzero);
  if(s.n==s.nzero)
  {
    drift(clevel,s,etime, dt);
    return;
  }
  for(int k=0;k<3;k++) cmpos[k]= s.part->pos[k] + s.part->vel[k]*dt;

  err=0;
#pragma omp parallel for if((ULONG) s.n>omp_get_num_threads() && !omp_in_parallel()) default(none) \
 private(ipart, dpos,dvel,dpos0,dvel0) shared(clevel, dt,cmpos, s, eps2) reduction(+: err)
  for(UINT i=1;i<s.n;i++)
  {
    ipart=GETPART(s,i);
    for(int k=0;k<3;k++) dpos0[k] = s.part->pos[k] - ipart->pos[k];
    for(int k=0;k<3;k++) dvel0[k] = s.part->vel[k] - ipart->vel[k];
    err+=universal_kepler_solver(dt,s.part->mass,eps2,
                                      dpos0[0],dpos0[1],dpos0[2],
                                      dvel0[0],dvel0[1],dvel0[2],
                                      &dpos[0],&dpos[1],&dpos[2],
                                      &dvel[0],&dvel[1],&dvel[2]);


    for(int k=0;k<3;k++) ipart->pos[k] = cmpos[k] - dpos[k];
    for(int k=0;k<3;k++) ipart->vel[k] = s.part->vel[k] - dvel[k];
  } 
  if (err != 0) {
    ENDRUN("kepler solver failure"); // failure of the kepler solver should be very rare now
  }
  for(int k=0;k<3;k++) s.part->pos[k]=cmpos[k];
  diag->cecount[clevel]+=s.nzero;
}


// special solver to test kepler
static void evolve_kepler_test(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  struct particle *central;
  struct particle *ipart;
  struct particle p[2];
  struct sys s2;

  CHECK_TIMESTEP(etime,stime,dt,clevel);
  if (s.n <= 1) ENDRUN("kepler test solver was called with too few massive particles sys.n=%u\n this hsouldn't happen\n", s.n);

  central=s.part;
  for(UINT i=1;i<s.n;i++)
  {
    ipart=GETPART(s,i);
    if(central->mass<ipart->mass) central=ipart;
  }

  for(UINT i=0;i<s.n;i++)
  {
    ipart=GETPART(s,i);
    if(ipart==central) continue;
    p[0]=*central;
    p[1]=*ipart;

    p[0].mass=central->mass+ipart->mass;
    p[1].mass=0;

    s2.n=2;
    s2.part=&p[0];
    s2.last=&p[1];
    evolve_kepler_2(clevel,s2,stime,etime,dt);
    p[1].mass=ipart->mass;
    *ipart=p[1];
  }
  for(int k=0;k<3;k++) central->pos[k]+=central->vel[k]*dt;
}

void evolve_kepler(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  if(s.n==2) // 2 body 
  {
    evolve_kepler_2(clevel,s,stime,etime,dt);
    return;
  }
  if(s.n-s.nzero==1 && s.nzero>0) // 1 massive, n orbiters
  {
    evolve_kepler_n(clevel,s,stime,etime,dt);
    return;
  }
  if(s.n-s.nzero>1) // more than 1 massive particle, consider heaviest as central;  
  {
    evolve_kepler_test(clevel,s,stime,etime,dt);
    return;
  }
  drift(clevel,s,etime, dt); // 1 massive or only zero mass
}
