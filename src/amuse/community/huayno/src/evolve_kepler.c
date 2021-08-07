/*
 * The Kepler solver evolves the two-body problem, using standard Huayno sys structures.
 * TODO add warning when using the kepler solver with softening?
 */

#include <stdio.h>
#include <stdlib.h>
#include "evolve.h"
#include "universal_kepler_solver.h"
#ifdef _OPENMP
#include <omp.h>
#endif

static void evolve_kepler_2(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  struct particle *ipart,*jpart;
  CHECK_TIMESTEP(etime,stime,dt,clevel);
  if (s.n != 2) ENDRUN("two-body solver was called with sys.n=%u\n", s.n);
  // translate coordinates original frame to 2-body frame
  int k;
  ipart=GETPART(s,0);
  jpart=GETPART(s,1);
  DOUBLE dpos[3],dpos0[3],pos_cm[3];
  DOUBLE dvel[3],dvel0[3],vel_cm[3];
  DOUBLE m1 = ipart->mass;
  DOUBLE m2 = jpart->mass;
  DOUBLE mtot = ipart->mass + jpart->mass;
  DOUBLE f1 = m2 / mtot;
  DOUBLE f2 = m1 / mtot;
  if(mtot>0.) {
    for(k=0;k<3;k++) dpos0[k] = ipart->pos[k] - jpart->pos[k];
    for(k=0;k<3;k++) dvel0[k] = ipart->vel[k] - jpart->vel[k];
    for(k=0;k<3;k++) pos_cm[k] = (m1 * ipart->pos[k] + m2 * jpart->pos[k]) / mtot;    
    for(k=0;k<3;k++) vel_cm[k] = (m1 * ipart->vel[k] + m2 * jpart->vel[k]) / mtot;
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
    for(k=0;k<3;k++) ipart->pos[k] = pos_cm[k] + f1 * dpos[k];
    for(k=0;k<3;k++) ipart->vel[k] = vel_cm[k] + f1 * dvel[k];
    for(k=0;k<3;k++) jpart->pos[k] = pos_cm[k] - f2 * dpos[k];
    for(k=0;k<3;k++) jpart->vel[k] = vel_cm[k] - f2 * dvel[k];
#ifdef COMPENSATED_SUMMP
    for(k=0;k<3;k++) ipart->pos_e[k]=0.;
    for(k=0;k<3;k++) jpart->pos_e[k]=0.;
#endif
#ifdef COMPENSATED_SUMMV
    for(k=0;k<3;k++) ipart->vel_e[k]=0.;
    for(k=0;k<3;k++) jpart->vel_e[k]=0.;
#endif
    
  } else {
    for(k=0;k<3;k++) COMPSUMP(ipart->pos[k],ipart->pos_e[k],dt*ipart->vel[k]);
    for(k=0;k<3;k++) COMPSUMP(jpart->pos[k],jpart->pos_e[k],dt*jpart->vel[k]);   
  } 
  ipart->postime=etime;
  jpart->postime=etime;

  diag->cecount[clevel]++;
}

static void evolve_kepler_n(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  struct particle *ipart, *spart;
  DOUBLE dpos[3],dpos0[3],spos[3];
  DOUBLE dvel[3],dvel0[3];
  UINT err;

  CHECK_TIMESTEP(etime,stime,dt,clevel);
  if (s.n-s.nzero > 1) ENDRUN("kepler-n solver was called with too many massive particles sys.n=%u\n", s.n-s.nzero);
  if (s.n-s.nzero < 1) ENDRUN("kepler-n solver was called with too few massive particles sys.n=%u\n", s.n-s.nzero);
  if(s.n==s.nzero)
  {
    drift(clevel,s,etime, dt);
    return;
  }
  spart=GETPART(s,0);
  for(int k=0;k<3;k++) spos[k]= spart->pos[k]; // save initial pos
  for(int k=0;k<3;k++) COMPSUMP(spart->pos[k],spart->pos_e[k],dt*spart->vel[k]); //evolve central
  spart->postime=etime;

  err=0;
#pragma omp parallel for if((ULONG) s.n>omp_get_num_threads() && !omp_in_parallel()) default(none) \
 private(ipart, dpos,dvel,dpos0,dvel0) shared(etime,clevel, dt,spos, s, eps2, spart) reduction(|: err)
  for(UINT i=1;i<s.n;i++)
  {
    ipart=GETPART(s,i);
    for(int k=0;k<3;k++) dpos0[k] = spos[k] - ipart->pos[k];
    for(int k=0;k<3;k++) dvel0[k] = spart->vel[k] - ipart->vel[k];
    err|=universal_kepler_solver(dt,spart->mass,eps2,
                                      dpos0[0],dpos0[1],dpos0[2],
                                      dvel0[0],dvel0[1],dvel0[2],
                                      &dpos[0],&dpos[1],&dpos[2],
                                      &dvel[0],&dvel[1],&dvel[2]);


    for(int k=0;k<3;k++) ipart->pos[k] = spart->pos[k] - dpos[k];
    for(int k=0;k<3;k++) ipart->vel[k] = spart->vel[k] - dvel[k];
#ifdef COMPENSATED_SUMMP
    for(int k=0;k<3;k++) ipart->pos_e[k]=0.;
#endif
#ifdef COMPENSATED_SUMMV
    for(int k=0;k<3;k++) ipart->vel_e[k]=0.;
#endif
    ipart->postime=etime;
  } 
  if (err != 0) {
    ENDRUN("kepler solver failure"); // failure of the kepler solver should be very rare now
  }
  diag->cecount[clevel]+=s.nzero;
}

void evolve_kepler(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  if(s.n-s.nzero==2) // 2 body 
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
    ENDRUN("evolve_kepler called for a system with more than 1 massive particle");
    return;
  }
  drift(clevel,s,etime, dt); // 1 massive or only zero mass
}
