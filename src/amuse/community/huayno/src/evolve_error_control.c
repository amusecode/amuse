/*
 * integrator which combines multiple integrators at different timestep/order to get error estimates
 */


#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#include "evolve.h"
#include "evolve_shared.h"

#define error_control_parameter   1.e-12

static FLOAT error_function(struct sys s1, struct sys s2)
{
  FLOAT maxdiv=0.;
  UINT i;
  struct particle *part1, *part2;
  if(s1.n!=s2.n) ENDRUN("error_function length mismatch %d,%d\n",s1.n,s2.n);
  for(i=0;i<s1.n;i++)
  {
    part1=GETPART(s1, i);
    part2=GETPART(s2, i);
    maxdiv=fmax(maxdiv,fabs( (FLOAT) part1->pos[0] - (FLOAT) part2->pos[0])); 
    maxdiv=fmax(maxdiv,fabs( (FLOAT) part1->pos[1] - (FLOAT) part2->pos[1])); 
    maxdiv=fmax(maxdiv,fabs( (FLOAT) part1->pos[2] - (FLOAT) part2->pos[2])); 
    maxdiv=fmax(maxdiv,fabs( (FLOAT) part1->vel[0] - (FLOAT) part2->vel[0])); 
    maxdiv=fmax(maxdiv,fabs( (FLOAT) part1->vel[1] - (FLOAT) part2->vel[1])); 
    maxdiv=fmax(maxdiv,fabs( (FLOAT) part1->vel[2] - (FLOAT) part2->vel[2])); 
  }
  return maxdiv;
}

void evolve_error_control_sub(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  struct sys s1,s2;
  CHECK_TIMESTEP(etime,stime,dt,clevel);
  
  s1.n=s.n; 
  s1.nzero=s.nzero;
  s1.part=(struct particle*) malloc(s.n*sizeof(struct particle));
  if(!s1.part) ENDRUN("failed allocation of sys\n");
  if(s1.nzero>0) s1.zeropart=s1.part+s1.n-s1.nzero;
  for(UINT i=0;i<s.n;i++) *GETPART(s1,i)=*GETPART(s,i);

  s2.n=s.n; 
  s2.nzero=s.nzero;
  s2.part=(struct particle*) malloc(s.n*sizeof(struct particle));
  if(!s2.part) ENDRUN("failed allocation of sys\n");
  if(s2.nzero>0) s2.zeropart=s2.part+s2.n-s2.nzero;
  for(UINT i=0;i<s.n;i++) *GETPART(s2,i)=*GETPART(s,i);

  //~ evolve_constant8(clevel, s1, stime,etime,dt);
  evolve_constant10(clevel, s1, stime,etime,dt);

  evolve_constant10(clevel+1, s2, stime, stime+dt/2, dt/2);
  evolve_constant10(clevel+1, s2, stime+dt/2, etime, dt/2);

  if(error_function(s1,s2)>error_control_parameter)
  {
    //~ LOG("%d error control REJECT\n", clevel);
    evolve_error_control_sub(clevel+1,s2, stime, stime+dt/2, dt/2);
    evolve_error_control_sub(clevel+1,s2, stime+dt/2, etime, dt/2);
  }// else LOG("%d error control ACCEPT\n", clevel);

  for(UINT i=0;i<s.n;i++) *GETPART(s,i)=*GETPART(s2,i);
  free(s1.part);
  free(s2.part);  
}

void evolve_error_control(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, FLOAT dtsys)
{
  CHECK_TIMESTEP(etime,stime,dt,clevel);
  if(dtsys<0) 
  {
    timestep(clevel,s,s,SIGN(dt));
    dtsys=global_timestep(s);
  }
  if(dtsys < fabs(dt))
  {
    evolve_error_control(clevel+1,s,stime, stime+dt/2,dt/2, dtsys);
    evolve_error_control(clevel+1,s,stime+dt/2, etime,dt/2, -1.);
  }
  else
  {
    evolve_error_control_sub(clevel, s, stime,etime,dt);
  }
}
