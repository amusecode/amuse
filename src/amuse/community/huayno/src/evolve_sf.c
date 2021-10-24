#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#include "evolve.h"

static void drift_naive(int clevel,struct sys s, DOUBLE etime); /* drift/extrap sys to itime*/
static void split(FLOAT dt, struct sys s, struct sys *slow, struct sys *fast);

void evolve_split_pass(int clevel,struct sys sys1,struct sys sys2, 
                         DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  struct sys slow=zerosys,fast=zerosys;
  CHECK_TIMESTEP(etime,stime,dt,clevel);
  if(calc_timestep) timestep(clevel,sys1, join(sys1,sys2), SIGN(dt));
//  if(calc_timestep) timestep(clevel,sys1, sys1, SIGN(dt));
  split((FLOAT) dt, sys1, &slow, &fast);
  if(fast.n==0) 
  {
    diag->deepsteps++;
    diag->simtime+=dt;
  }  
  if(fast.n>0) evolve_split_pass(clevel+1,fast, join(slow,sys2), stime, stime+dt/2, dt/2,0);
  if(slow.n>0) kdk(clevel,slow,sys2, stime, etime, dt);
  if(fast.n>0) evolve_split_pass(clevel+1,fast, join(slow,sys2), stime+dt/2, etime, dt/2,1);
}

void evolve_split_naive(int clevel,struct sys sys1,struct sys sys2, 
                          DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  struct sys slow=zerosys,fast=zerosys;
  CHECK_TIMESTEP(etime,stime,dt,clevel);
  if(calc_timestep) timestep(clevel,sys1, join(sys1,sys2), SIGN(dt));
//  if(calc_timestep) timestep(clevel,sys1, sys1, SIGN(dt));
  split((FLOAT) dt, sys1, &slow, &fast);
  if(fast.n==0) 
  {
    diag->deepsteps++;
    diag->simtime+=dt;
  }  
  if(slow.n>0) kick(clevel,slow, join(sys1,sys2), dt/2);
  if(fast.n>0) evolve_split_naive(clevel+1,fast, join(slow,sys2), stime, stime+dt/2, dt/2,0);
  if(slow.n>0) drift_naive(clevel,join(slow,sys2),stime+dt/2);
  if(fast.n>0) evolve_split_naive(clevel+1,fast, join(slow,sys2), stime+dt/2, etime, dt/2,1);
  if(slow.n>0) drift_naive(clevel,join(slow,sys2),etime);
  if(slow.n>0) kick(clevel,slow, join(sys1,sys2), dt/2);
}

void evolve_split_bridge(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  struct sys slow=zerosys,fast=zerosys;
  CHECK_TIMESTEP(etime,stime,dt,clevel);
  if(calc_timestep) timestep(clevel,s,s, SIGN(dt));
  split((FLOAT) dt, s, &slow, &fast);
  if(fast.n==0) 
  {
    diag->deepsteps++;
    diag->simtime+=dt;
  }  
  if(slow.n>0 && fast.n>0) kick(clevel,slow, fast, dt/2);
  if(slow.n>0 && fast.n>0) kick(clevel,fast, slow, dt/2);
  if(fast.n>0) evolve_split_bridge(clevel+1,fast, stime, stime+dt/2, dt/2,0); /* note calc_timestep? */
  if(slow.n>0) kdk(clevel,slow,zerosys, stime, etime, dt);
  if(fast.n>0) evolve_split_bridge(clevel+1,fast, stime+dt/2, etime, dt/2,1);
  if(slow.n>0 && fast.n>0) kick(clevel,slow, fast, dt/2);
  if(slow.n>0 && fast.n>0) kick(clevel,fast, slow, dt/2);
}

void evolve_split_bridge_dkd(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  struct sys slow=zerosys,fast=zerosys;
  CHECK_TIMESTEP(etime,stime,dt,clevel);
  if(calc_timestep) timestep(clevel,s,s, SIGN(dt));
  split((FLOAT) dt, s, &slow, &fast);
  if(fast.n==0) 
  {
    diag->deepsteps++;
    diag->simtime+=dt;
  }  
  if(slow.n>0 && fast.n>0) kick(clevel,slow, fast, dt/2);
  if(slow.n>0 && fast.n>0) kick(clevel,fast, slow, dt/2);
  if(fast.n>0) evolve_split_bridge_dkd(clevel+1,fast, stime, stime+dt/2, dt/2,0); /* note calc_timestep? */
  if(slow.n>0) dkd(clevel,slow,zerosys, stime, etime, dt);
  if(fast.n>0) evolve_split_bridge_dkd(clevel+1,fast, stime+dt/2, etime, dt/2,1);
  if(slow.n>0 && fast.n>0) kick(clevel,slow, fast, dt/2);
  if(slow.n>0 && fast.n>0) kick(clevel,fast, slow, dt/2);
}

void evolve_split_hold(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  struct sys slow=zerosys,fast=zerosys;
  CHECK_TIMESTEP(etime,stime,dt,clevel);
  if(calc_timestep) timestep(clevel,s,s, SIGN(dt));
  split((FLOAT) dt, s, &slow, &fast);
  if(fast.n==0) 
  {
    diag->deepsteps++;
    diag->simtime+=dt;
  }  
  if(fast.n>0) evolve_split_hold(clevel+1,fast, stime, stime+dt/2, dt/2,0);
  if(slow.n>0) kdk(clevel,slow,fast,stime,etime,dt);
  if(fast.n>0) evolve_split_hold(clevel+1,fast, stime+dt/2, etime, dt/2,1);
}

void evolve_split_hold_dkd(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  struct sys slow=zerosys,fast=zerosys;
  CHECK_TIMESTEP(etime,stime,dt,clevel);
  if(calc_timestep) timestep(clevel,s,s, SIGN(dt));
  split((FLOAT) dt, s, &slow, &fast);
  if(fast.n==0) 
  {
    diag->deepsteps++;
    diag->simtime+=dt;
  }  
  if(fast.n>0) evolve_split_hold_dkd(clevel+1,fast, stime, stime+dt/2, dt/2,0);
  if(slow.n>0) dkd(clevel,slow,fast,stime,etime,dt);
  if(fast.n>0) evolve_split_hold_dkd(clevel+1,fast, stime+dt/2, etime, dt/2,1);
}

void evolve_split_pass_dkd(int clevel,struct sys sys1,struct sys sys2, 
                         DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  struct sys slow=zerosys,fast=zerosys;
  CHECK_TIMESTEP(etime,stime,dt,clevel);
  if(calc_timestep) timestep(clevel,sys1, join(sys1,sys2), SIGN(dt));
//  if(calc_timestep) timestep(clevel,sys1, sys1, SIGN(dt));
  split((FLOAT) dt, sys1, &slow, &fast);
  if(fast.n==0) 
  {
    diag->deepsteps++;
    diag->simtime+=dt;
  }  
  if(fast.n>0) evolve_split_pass_dkd(clevel+1,fast, join(slow,sys2), stime, stime+dt/2, dt/2,0);
  if(slow.n>0) dkd(clevel,slow,sys2, stime, etime, dt);
  if(fast.n>0) evolve_split_pass_dkd(clevel+1,fast, join(slow,sys2), stime+dt/2, etime, dt/2,1);
}

void evolve_split_ppass_dkd(int clevel,struct sys sys1,struct sys sys2, 
                         DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  struct sys slow=zerosys,fast=zerosys;
  CHECK_TIMESTEP(etime,stime,dt,clevel);
  if(calc_timestep) timestep(clevel,sys1, join(sys1,sys2), SIGN(dt));
//  if(calc_timestep) timestep(clevel,sys1, sys1, SIGN(dt));
  split((FLOAT) dt, sys1, &slow, &fast);
  if(fast.n==0) 
  {
    diag->deepsteps++;
    diag->simtime+=dt;
  }  
  if(fast.n>0) 
    evolve_split_ppass_dkd(clevel+1,fast, join(slow,sys2), stime, stime+dt/2, dt/2,0);
  else  
    drift(clevel,join(slow,sys2),etime,dt/2);
  if(slow.n>0) kick(clevel,slow,join(slow,sys2), dt);
  if(slow.n>0) kick(clevel,sys2,slow, dt);
  if(fast.n>0) 
    evolve_split_ppass_dkd(clevel+1,fast, join(slow,sys2), stime+dt/2, etime, dt/2,1);
  else  
    drift(clevel,join(slow,sys2),etime,dt/2);
}


static void drift_naive(int clevel,struct sys s, DOUBLE etime)
{
  UINT i;
  DOUBLE dt;
  struct particle *ipart;
  for(i=0;i<s.n;i++)
  {
    ipart=GETPART(s,i);
    dt=etime-ipart->postime;
    COMPSUMP(ipart->pos[0],ipart->pos_e[0],dt*ipart->vel[0])
    COMPSUMP(ipart->pos[1],ipart->pos_e[1],dt*ipart->vel[1])
    COMPSUMP(ipart->pos[2],ipart->pos_e[2],dt*ipart->vel[2])
    ipart->postime=etime;
  }
  diag->dstep[clevel]++;
  diag->dcount[clevel]+=s.n;
}

#define K1   ( (14-sqrt((DOUBLE) 19))/108 )
#define K2   ( (20-7*sqrt((DOUBLE) 19))/108 )
#define K3   ( (((DOUBLE) 1)/2)-(K1+K2) )
#define D1   ( ((DOUBLE) 2)/5 )
#define D2   ( -((DOUBLE) 1)/10 )
#define D3   ( 1-(2*D1+2*D2) )
#define N1   4
void evolve_sf_4m5(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  int i;
  struct sys slow=zerosys,fast=zerosys;
  CHECK_TIMESTEP(etime,stime,dt,clevel);
  if(calc_timestep) timestep(clevel,s,s, SIGN(dt));
  split((FLOAT) dt, s, &slow, &fast);
  if(fast.n==0) 
  {
    diag->deepsteps++;
    diag->simtime+=dt;
    diag->timetrack+=fabs(dt);
  }  
  if(slow.n>0) kick(clevel,slow,join(fast,slow), K1*dt);
  if(fast.n>0 && slow.n>0) kick(clevel,fast,slow, K1*dt);

  if(slow.n>0) drift(clevel,slow,stime+D1*dt,D1*dt);
  if(fast.n>0) 
    for(i=0;i<N1;i++)
    {
      evolve_sf_4m5(clevel+1,fast,stime,stime+D1*dt/N1,D1*dt/N1,i==0?0:1);
      stime+=D1*dt/N1; 
    }
  if(slow.n>0) kick(clevel,slow,join(fast,slow), K2*dt);
  if(fast.n>0 && slow.n>0) kick(clevel,fast,slow, K2*dt);

  if(slow.n>0) drift(clevel,slow,stime+D2*dt,D2*dt);
  if(fast.n>0) evolve_sf_4m5(clevel+1,fast,stime,stime+D2*dt,D2*dt,1);
  stime+=D2*dt;

  if(slow.n>0) kick(clevel,slow,join(fast,slow), K3*dt);
  if(fast.n>0 && slow.n>0) kick(clevel,fast,slow, K3*dt);

  if(slow.n>0) drift(clevel,slow,stime+D3*dt,D3*dt);
  if(fast.n>0) 
    for(i=0;i<N1;i++)
    {
      evolve_sf_4m5(clevel+1,fast,stime,stime+D3*dt/N1,D3*dt/N1,1);
      stime+=D3*dt/N1;
    }
  if(slow.n>0) kick(clevel,slow,join(fast,slow), K3*dt);
  if(fast.n>0 && slow.n>0) kick(clevel,fast,slow, K3*dt);

  if(slow.n>0) drift(clevel,slow,stime+D2*dt,D2*dt);
  if(fast.n>0) evolve_sf_4m5(clevel+1,fast,stime,stime+D2*dt,D2*dt,1);
  stime+=D2*dt;

  if(slow.n>0) kick(clevel,slow,join(fast,slow), K2*dt);
  if(fast.n>0 && slow.n>0) kick(clevel,fast,slow, K2*dt);

  if(slow.n>0) drift(clevel,slow,etime,D1*dt);
  if(fast.n>0) 
    for(i=0;i<N1;i++)
    {
      evolve_sf_4m5(clevel+1,fast,stime,i==N1-1?etime:stime+D1*dt/N1,D1*dt/N1,1);
      stime+=D1*dt/N1; 
    }
  if(slow.n>0) kick(clevel,slow,join(fast,slow), K1*dt);
  if(fast.n>0 && slow.n>0) kick(clevel,fast,slow, K1*dt);    
}
#undef K1
#undef K2
#undef K3
#undef D1
#undef D2
#undef D3
#undef N1

#define K1   ( (642+sqrt( (DOUBLE) 471 ))/3924 )
#define K2   (  121*(12- sqrt( (DOUBLE) 471 ) )/3924 )
#define K3   (  1-2*(K1+K2) )
#define D1   ( ((DOUBLE) 6)/11 )
#define D2   ( ((DOUBLE) 0.5)-D1 ) 
#define N1   12
void evolve_sf_4m4(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  int i;
  struct sys slow=zerosys,fast=zerosys;
  CHECK_TIMESTEP(etime,stime,dt,clevel);
  if(calc_timestep) timestep(clevel,s,s, SIGN(dt));
  split((FLOAT) dt, s, &slow, &fast);
  if(fast.n==0) 
  {
    diag->deepsteps++;
    diag->simtime+=dt;
    diag->timetrack+=fabs(dt);
  }  
  if(slow.n>0) kick(clevel,slow,join(fast,slow), K1*dt);
  if(fast.n>0 && slow.n>0) kick(clevel,fast,slow, K1*dt);

  if(slow.n>0) drift(clevel,slow,stime+D1*dt,D1*dt);
  if(fast.n>0) 
    for(i=0;i<N1;i++)
    {
      evolve_sf_4m5(clevel+1,fast,stime,stime+D1*dt/N1,D1*dt/N1,i==0?0:1);
      stime+=D1*dt/N1; 
    }
  if(slow.n>0) kick(clevel,slow,join(fast,slow), K2*dt);
  if(fast.n>0 && slow.n>0) kick(clevel,fast,slow, K2*dt);

  if(slow.n>0) drift(clevel,slow,stime+D2*dt,D2*dt);
  if(fast.n>0) evolve_sf_4m5(clevel+1,fast,stime,stime+D2*dt,D2*dt,1);
  stime+=D2*dt;

  if(slow.n>0) kick(clevel,slow,join(fast,slow), K3*dt);
  if(fast.n>0 && slow.n>0) kick(clevel,fast,slow, K3*dt);

  if(slow.n>0) drift(clevel,slow,stime+D2*dt,D2*dt);
  if(fast.n>0) evolve_sf_4m5(clevel+1,fast,stime,stime+D2*dt,D2*dt,1);
  stime+=D2*dt;

  if(slow.n>0) kick(clevel,slow,join(fast,slow), K2*dt);
  if(fast.n>0 && slow.n>0) kick(clevel,fast,slow, K2*dt);

  if(slow.n>0) drift(clevel,slow,etime,D1*dt);
  if(fast.n>0) 
    for(i=0;i<N1;i++)
    {
      evolve_sf_4m5(clevel+1,fast,stime,i==N1-1?etime:stime+D1*dt/N1,D1*dt/N1,1);
      stime+=D1*dt/N1; 
    }
  if(slow.n>0) kick(clevel,slow,join(fast,slow), K1*dt);
  if(fast.n>0 && slow.n>0) kick(clevel,fast,slow, K1*dt);    
}
#undef K1
#undef K2
#undef K3
#undef D1
#undef D2

// partitions a contiguous array bordered by left and right according to a pivot timestep dt 
static struct particle *partition(FLOAT dt, struct particle *left, struct particle *right)
{
  UINT i,n;
  dt=fabs(dt);
  n=right-left+1;
  while(1)
  {
    if(i>=n) ENDRUN("partition error");
    i++;
    while(left->timestep<dt && left<right) left++;
    while(right->timestep>=dt && left<right) right--;
    if(left<right) 
      {SWAP( *left, *right, struct particle);}
    else 
      break;  
  }
  if(left->timestep<dt) left++;
  return left;
}

static void split(FLOAT dt, struct sys s, struct sys *slow, struct sys *fast)
{
  struct particle *left, *right, *pivot;
  slow->n=0;  
  fast->n=0;  
  if(s.n-s.nzero>0)
  {
    left=s.part;
    right=LAST(s);
    pivot=partition(dt, left, right);
    slow->n=right-pivot+1;
    fast->n=(pivot-left);
  }
  slow->nzero=0;  
  fast->nzero=0;  
  if(s.nzero>0)
  {
    left=s.zeropart;
    right=LASTZERO(s);
    pivot=partition(dt, left, right);
    slow->nzero=right-pivot+1;
    fast->nzero=(pivot-left);  
    slow->n+=slow->nzero;
    fast->n+=fast->nzero;  
  }
  if(fast->n<=1)
  {
    *fast=zerosys;
    slow->n=s.n;
    slow->nzero=s.nzero;
  } 
  if(slow->n>0)
  {
    slow->part=s.part+fast->n-fast->nzero;
  }
  if(fast->n>0)
  {
    fast->part=s.part;
  }
  if(slow->nzero>0)
  {
    slow->zeropart=s.zeropart+fast->nzero;
  }
  if(fast->nzero > 0)
  {
    fast->zeropart=s.zeropart;
  }
  if(fast->n+slow->n !=s.n) ENDRUN( "split error 2");
  if(fast->nzero+slow->nzero !=s.nzero) ENDRUN( "split error 3");
}

