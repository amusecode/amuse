#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "evolve.h"

#ifdef EVOLVE_OPENCL
#include "evolve_cl.h"
#endif

FLOAT eps2;
FLOAT dt_param;
struct sys zerosys ={ 0, NULL,NULL};

/* diagnostics */
DOUBLE simtime;
/*static*/ int clevel;
/*static*/ unsigned long tcount[MAXLEVEL],kcount[MAXLEVEL],dcount[MAXLEVEL],deepsteps;
/*static*/ unsigned long tstep[MAXLEVEL],kstep[MAXLEVEL],dstep[MAXLEVEL];
#ifdef EVOLVE_OPENCL
/*static*/ unsigned long cpu_step,cl_step,cpu_count,cl_count;
#endif

static void split(FLOAT dt, struct sys s, struct sys *slow, struct sys *fast);
static struct sys join(struct sys s1,struct sys s2);
static void drift_naive(struct sys s, DOUBLE etime); /* drift/extrap sys to itime*/
static void potential(struct sys s1, struct sys s2);
/*static*/ void kdk(struct sys s1,struct sys s2, DOUBLE stime, DOUBLE etime, DOUBLE dt);
static void dkd(struct sys s1,struct sys s2, DOUBLE stime, DOUBLE etime, DOUBLE dt);
static void report(struct sys s,DOUBLE etime, int inttype);

FLOAT system_kinetic_energy(struct sys s)
{
 UINT i;
 DOUBLE e=0.;
 for(i=0;i<s.n;i++) e+=0.5*s.part[i].mass*(
                           s.part[i].vel[0]*s.part[i].vel[0]+
                           s.part[i].vel[1]*s.part[i].vel[1]+
                           s.part[i].vel[2]*s.part[i].vel[2]);
 return (FLOAT) e;
}
        
FLOAT system_potential_energy(struct sys s)
{
 UINT i;
 DOUBLE e=0.;
 for(i=0;i<s.n;i++) e+=s.part[i].mass*s.part[i].pot;
 return (FLOAT) e/2;
}

void evolve_split_pass(struct sys sys1,struct sys sys2, 
                         DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  struct sys slow=zerosys,fast=zerosys;
  clevel++;
  if(etime <= stime ||  dt==0 || clevel>=MAXLEVEL) ENDRUN("timestep too small");
  if(calc_timestep) timestep(sys1, join(sys1,sys2));
//  if(calc_timestep) timestep(sys1, sys1);
  split((FLOAT) dt, sys1, &slow, &fast);
  if(fast.n==0) 
  {
    deepsteps++;
    simtime+=dt;
  }  
  if(fast.n>0) evolve_split_pass(fast, join(slow,sys2), stime, stime+dt/2, dt/2,0);
  if(slow.n>0) kdk(slow,sys2, stime, etime, dt);
  if(fast.n>0) evolve_split_pass(fast, join(slow,sys2), stime+dt/2, etime, dt/2,1);
  clevel--;
}

void evolve_split_naive(struct sys sys1,struct sys sys2, 
                          DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  struct sys slow=zerosys,fast=zerosys;
  clevel++;
  if(etime <= stime ||  dt==0 || clevel>=MAXLEVEL) ENDRUN("timestep too small");
  if(calc_timestep) timestep(sys1, join(sys1,sys2));
//  if(calc_timestep) timestep(sys1, sys1);
  split((FLOAT) dt, sys1, &slow, &fast);
  if(fast.n==0) 
  {
    deepsteps++;
    simtime+=dt;
  }  
  if(slow.n>0) kick(slow, join(sys1,sys2), dt/2);
  if(fast.n>0) evolve_split_naive(fast, join(slow,sys2), stime, stime+dt/2, dt/2,0);
  if(slow.n>0) drift_naive(join(slow,sys2),stime+dt/2);
  if(fast.n>0) evolve_split_naive(fast, join(slow,sys2), stime+dt/2, etime, dt/2,1);
  if(slow.n>0) drift_naive(join(slow,sys2),etime);
  if(slow.n>0) kick(slow, join(sys1,sys2), dt/2);
  clevel--;
}

void evolve_split_bridge(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  struct sys slow=zerosys,fast=zerosys;
  clevel++;
  if(etime <= stime ||  dt==0 || clevel>=MAXLEVEL) ENDRUN("timestep too small");
  if(calc_timestep) timestep(s,s);
  split((FLOAT) dt, s, &slow, &fast);
  if(fast.n==0) 
  {
    deepsteps++;
    simtime+=dt;
  }  
  if(slow.n>0 && fast.n>0) kick(slow, fast, dt/2);
  if(slow.n>0 && fast.n>0) kick(fast, slow, dt/2);
  if(fast.n>0) evolve_split_bridge(fast, stime, stime+dt/2, dt/2,0); /* note calc_timestep? */
  if(slow.n>0) kdk(slow,zerosys, stime, etime, dt);
  if(fast.n>0) evolve_split_bridge(fast, stime+dt/2, etime, dt/2,1);
  if(slow.n>0 && fast.n>0) kick(slow, fast, dt/2);
  if(slow.n>0 && fast.n>0) kick(fast, slow, dt/2);
  clevel--;
}

void evolve_split_bridge_dkd(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  struct sys slow=zerosys,fast=zerosys;
  clevel++;
  if(etime <= stime ||  dt==0 || clevel>=MAXLEVEL) ENDRUN("timestep too small");
  if(calc_timestep) timestep(s,s);
  split((FLOAT) dt, s, &slow, &fast);
  if(fast.n==0) 
  {
    deepsteps++;
    simtime+=dt;
  }  
  if(slow.n>0 && fast.n>0) kick(slow, fast, dt/2);
  if(slow.n>0 && fast.n>0) kick(fast, slow, dt/2);
  if(fast.n>0) evolve_split_bridge_dkd(fast, stime, stime+dt/2, dt/2,0); /* note calc_timestep? */
  if(slow.n>0) dkd(slow,zerosys, stime, etime, dt);
  if(fast.n>0) evolve_split_bridge_dkd(fast, stime+dt/2, etime, dt/2,1);
  if(slow.n>0 && fast.n>0) kick(slow, fast, dt/2);
  if(slow.n>0 && fast.n>0) kick(fast, slow, dt/2);
  clevel--;
}

void evolve_split_hold(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  struct sys slow=zerosys,fast=zerosys;
  clevel++;
  if(etime <= stime ||  dt==0 || clevel>=MAXLEVEL) ENDRUN("timestep too small");
  if(calc_timestep) timestep(s,s);
  split((FLOAT) dt, s, &slow, &fast);
  if(fast.n==0) 
  {
    deepsteps++;
    simtime+=dt;
  }  
  if(fast.n>0) evolve_split_hold(fast, stime, stime+dt/2, dt/2,0);
  if(slow.n>0) kdk(slow,fast,stime,etime,dt);
  if(fast.n>0) evolve_split_hold(fast, stime+dt/2, etime, dt/2,1);
  clevel--;
}

void evolve_split_hold_dkd(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  struct sys slow=zerosys,fast=zerosys;
  clevel++;
  if(etime <= stime ||  dt==0 || clevel>=MAXLEVEL) ENDRUN("timestep too small");
  if(calc_timestep) timestep(s,s);
  split((FLOAT) dt, s, &slow, &fast);
  if(fast.n==0) 
  {
    deepsteps++;
    simtime+=dt;
  }  
  if(fast.n>0) evolve_split_hold_dkd(fast, stime, stime+dt/2, dt/2,0);
  if(slow.n>0) dkd(slow,fast,stime,etime,dt);
  if(fast.n>0) evolve_split_hold_dkd(fast, stime+dt/2, etime, dt/2,1);
  clevel--;
}

void evolve_split_pass_dkd(struct sys sys1,struct sys sys2, 
                         DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  struct sys slow=zerosys,fast=zerosys;
  clevel++;
  if(etime <= stime ||  dt==0 || clevel>=MAXLEVEL) ENDRUN("timestep too small");
  if(calc_timestep) timestep(sys1, join(sys1,sys2));
//  if(calc_timestep) timestep(sys1, sys1);
  split((FLOAT) dt, sys1, &slow, &fast);
  if(fast.n==0) 
  {
    deepsteps++;
    simtime+=dt;
  }  
  if(fast.n>0) evolve_split_pass_dkd(fast, join(slow,sys2), stime, stime+dt/2, dt/2,0);
  if(slow.n>0) dkd(slow,sys2, stime, etime, dt);
  if(fast.n>0) evolve_split_pass_dkd(fast, join(slow,sys2), stime+dt/2, etime, dt/2,1);
  clevel--;
}

void evolve_split_ppass_dkd(struct sys sys1,struct sys sys2, 
                         DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  struct sys slow=zerosys,fast=zerosys;
  clevel++;
  if(etime <= stime ||  dt==0 || clevel>=MAXLEVEL) ENDRUN("timestep too small");
  if(calc_timestep) timestep(sys1, join(sys1,sys2));
//  if(calc_timestep) timestep(sys1, sys1);
  split((FLOAT) dt, sys1, &slow, &fast);
  if(fast.n==0) 
  {
    deepsteps++;
    simtime+=dt;
  }  
  if(fast.n>0) 
    evolve_split_ppass_dkd(fast, join(slow,sys2), stime, stime+dt/2, dt/2,0);
  else  
    drift(join(slow,sys2),etime,dt/2);
  if(slow.n>0) kick(slow,join(slow,sys2), dt);
  if(slow.n>0) kick(sys2,slow, dt);
  if(fast.n>0) 
    evolve_split_ppass_dkd(fast, join(slow,sys2), stime+dt/2, etime, dt/2,1);
  else  
    drift(join(slow,sys2),etime,dt/2);
  clevel--;
}


static void drift_naive(struct sys s, DOUBLE etime)
{
  UINT i;
  DOUBLE dt;
  for(i=0;i<s.n;i++)
  {
    dt=etime-s.part[i].postime;
    s.part[i].pos[0]+=dt*s.part[i].vel[0];
    s.part[i].pos[1]+=dt*s.part[i].vel[1];
    s.part[i].pos[2]+=dt*s.part[i].vel[2];
    s.part[i].postime=etime;
    s.part[i].timestep=HUGE_VAL;
  }
  dstep[clevel]++;
  dcount[clevel]+=s.n;
}

void init_code()
{
#ifdef EVOLVE_OPENCL
 init_cl();
#endif
}

void stop_code()
{
#ifdef EVOLVE_OPENCL
  close_cl();
#endif
  evolve_ok_stop(); // safe to call even if ok was not used
}

void init_evolve(struct sys s)
{
  UINT i;
  for(i=0;i<s.n;i++)
  {
    s.part[i].postime=0.;
    s.part[i].pot=0.;
    s.part[i].level=0;
    s.part[i].timestep=HUGE_VAL;
  }
  potential(s,s);
  if (inttype == OK) {
    evolve_ok_init(s);
  }
}

void do_evolve(struct sys s, double dt, int inttype)
{
  UINT p;
  int i;
  if(dt<=0) return;
  for(p=0;p<s.n;p++) s.part[p].postime=0.;
  for(p=0;p<s.n;p++) s.part[p].level=-1;
  clevel=-1;
  deepsteps=0;
  simtime=0.;
  for(i=0;i<MAXLEVEL;i++)
  {
    tstep[i]=0;tcount[i]=0;
    kstep[i]=0;kcount[i]=0;
    dstep[i]=0;dcount[i]=0;
  }
  switch (inttype)
  {
    case SHARED2:
      evolve_shared2(s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case SHARED4:
      evolve_shared4(s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case PASS:
      evolve_split_pass(s, zerosys,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case HOLD:
      evolve_split_hold(s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case BRIDGE:
      evolve_split_bridge(s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case NAIVE:
      evolve_split_naive(s, zerosys,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case HOLD_DKD:
      evolve_split_hold_dkd(s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case PASS_DKD:
      evolve_split_pass_dkd(s, zerosys, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case PPASS_DKD:
      evolve_split_ppass_dkd(s, zerosys, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case BRIDGE_DKD:
      evolve_split_bridge_dkd(s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case CC:
      evolve_cc2(s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt);
      break;
    case CC_KEPLER:
      evolve_cc2_kepler(s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt);
      break;
    case OK:
      evolve_ok2(s, zeroforces, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case KEPLER:
      evolve_kepler(s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt);
      break;
    default:  
      ENDRUN("unknown integrator\n");
      break;
  } 
  for(p=0;p<s.n;p++) s.part[p].pot=0;
  potential(s,s);
  report(s,(DOUBLE) dt, inttype);
}

static void split(FLOAT dt, struct sys s, struct sys *slow, struct sys *fast)
{
  UINT i=0;
  struct particle *left, *right;
  left=s.part;
  right=s.last;
  while(1)
  {
    if(i>=s.n) ENDRUN("split error 1");
    i++;
    while(left->timestep<dt && left<right) left++;
    while(right->timestep>=dt && left<right) right--;
    if(left<right) 
      {SWAP( *left, *right, struct particle);}
    else 
      break;  
  }
  if(left->timestep<dt) left++;
  slow->n=s.last-left+1;
  fast->n=(left-s.part);  
  if(fast->n==1)
  {
    fast->n=0;
    slow->n=s.n;
  } 
  if(slow->n > 0)
  {
    slow->part=s.part+fast->n;
    slow->last=s.last;//slow->part+slow->n-1;
  }
  if(fast->n > 0)
  {
    fast->part=s.part;
    fast->last=s.part+fast->n-1;
  }
  if(fast->n+slow->n !=s.n) ENDRUN( "split error 2");
  for(i=0;i<s.n;i++) s.part[i].level=clevel;
}

static struct sys join(struct sys s1,struct sys s2)
{
  struct sys s=zerosys;
  if(s1.n == 0) return s2;
  if(s2.n == 0) return s1;  
  s.n=s1.n+s2.n;
  if(s1.part+s1.n == s2.part)
  {
    s.part=s1.part;
    s.last=s2.last;
  } else
  {
    if(s2.part+s2.n == s1.part)
    {
      s.part=s2.part;
      s.last=s1.last;
    } else
      ENDRUN("join error 1");
  }   
  if(s.last-s.part + 1 != s.n) ENDRUN("join error 2");
  return s;
}

/*static*/ void kdk(struct sys s1,struct sys s2, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  if(s2.n>0) kick(s2, s1, dt/2);
  kick(s1,join(s1,s2),dt/2);
  drift(s1,etime, dt);
  kick(s1,join(s1,s2),dt/2);
  if(s2.n>0) kick(s2, s1, dt/2);
}

static void dkd(struct sys s1,struct sys s2, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  drift(s1,stime+dt/2, dt/2);
  kick(s1,join(s1,s2),dt);
  if(s2.n>0) kick(s2, s1, dt);
  drift(s1,etime, dt/2);
}

/*static*/ void drift(struct sys s, DOUBLE etime, DOUBLE dt)
{
  UINT i;
  for(i=0;i<s.n;i++)
  {
    s.part[i].pos[0]+=dt*s.part[i].vel[0];
    s.part[i].pos[1]+=dt*s.part[i].vel[1];
    s.part[i].pos[2]+=dt*s.part[i].vel[2];
    s.part[i].postime=etime;
    s.part[i].timestep=HUGE_VAL;
  }
  dstep[clevel]++;
  dcount[clevel]+=s.n;
}

static void kick_cpu(struct sys s1, struct sys s2, DOUBLE dt)
{
  UINT i,j;
  FLOAT dx[3],dr3,dr2,dr,acci;
  FLOAT acc[3];

#pragma omp parallel for if((ULONG) s1.n*s2.n>MPWORKLIMIT) default(none) \
 private(i,j,dx,dr3,dr2,dr,acc,acci) \
 shared(dt,s1,s2,eps2)
  for(i=0;i<s1.n;i++)
  {
    acc[0]=0.;
    acc[1]=0.;
    acc[2]=0.;
    for(j=0;j<s2.n;j++)
    {
//      if(s1.part+i==s2.part+j) continue; 
      dx[0]=s1.part[i].pos[0]-s2.part[j].pos[0];
      dx[1]=s1.part[i].pos[1]-s2.part[j].pos[1];
      dx[2]=s1.part[i].pos[2]-s2.part[j].pos[2];
      dr2=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+eps2;
      if(dr2>0) 
      {
        dr=sqrt(dr2);
        dr3=dr*dr2;
        acci=s2.part[j].mass/dr3;

        acc[0]-=dx[0]*acci;
        acc[1]-=dx[1]*acci;
        acc[2]-=dx[2]*acci;  
      }
    }
    s1.part[i].vel[0]+=dt*acc[0];
    s1.part[i].vel[1]+=dt*acc[1];
    s1.part[i].vel[2]+=dt*acc[2];
  }
}

/*static*/ void kick(struct sys s1, struct sys s2, DOUBLE dt)
{
#ifdef EVOLVE_OPENCL
  if((ULONG) s1.n*s2.n>CLWORKLIMIT) 
  {
    kick_cl(s1,s2,dt);
    cl_step++;
    cl_count+=(ULONG) s1.n*s2.n;
  } else
  {
    kick_cpu(s1,s2,dt);
    cpu_step++;
    cpu_count+=(ULONG) s1.n*s2.n;
  }
#else
  kick_cpu(s1,s2,dt);
#endif  
  kstep[clevel]++;
  kcount[clevel]+=(ULONG) s1.n*s2.n;
}

static void potential_cpu(struct sys s1,struct sys s2)
{
  UINT i,j;
  FLOAT dx[3],dr2,dr;
  FLOAT pot;

#pragma omp parallel for if((ULONG) s1.n*s2.n>MPWORKLIMIT) default(none) \
 private(i,j,dx,dr2,dr,pot) \
 shared(s1,s2,eps2)
  for(i=0;i<s1.n;i++)
  {
    pot=0;
    for(j=0;j<s2.n;j++)
    {
      if(s1.part+i==s2.part+j) continue; 
      dx[0]=s1.part[i].pos[0]-s2.part[j].pos[0];
      dx[1]=s1.part[i].pos[1]-s2.part[j].pos[1];
      dx[2]=s1.part[i].pos[2]-s2.part[j].pos[2];
      dr2=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+eps2;
      if(dr2>0) 
      {
        dr=sqrt(dr2);
        pot-=s2.part[j].mass/dr;
      }
    }
    s1.part[i].pot+=pot;
  }
}

static void potential(struct sys s1, struct sys s2)
{
#ifdef EVOLVE_OPENCL
  if((ULONG) s1.n*s2.n>CLWORKLIMIT) 
  {
    potential_cl(s1,s2);
  } else
  {
    potential_cpu(s1,s2);
  }
#else
  potential_cpu(s1,s2);
#endif  
}

static void timestep_cpu(struct sys s1, struct sys s2)
{
  UINT i,j;
  FLOAT timestep;
  FLOAT dx[3],dr3,dr2,dr,dv[3],dv2,mu,vdotdr2,tau,dtau;
#pragma omp parallel for if((ULONG) s1.n*s2.n>MPWORKLIMIT) default(none) \
 private(i,j,dx,dr3,dr2,dr,dv,dv2,mu,vdotdr2,\
           tau,dtau,timestep) \
 shared(s1,s2,eps2,dt_param)
  for(i=0;i<s1.n;i++)
  {
    if(s1.part[i].timestep !=HUGE_VAL) ENDRUN("timestep??");
  
    timestep=HUGE_VAL;
    for(j=0;j<s2.n;j++)
    {
      if(s1.part+i==s2.part+j) continue; 
      dx[0]=s1.part[i].pos[0]-s2.part[j].pos[0];
      dx[1]=s1.part[i].pos[1]-s2.part[j].pos[1];
      dx[2]=s1.part[i].pos[2]-s2.part[j].pos[2];
      dr2=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+eps2;
      if(dr2>0) 
      {
        dr=sqrt(dr2);
        dr3=dr*dr2;
        dv[0]=s1.part[i].vel[0]-s2.part[j].vel[0];
        dv[1]=s1.part[i].vel[1]-s2.part[j].vel[1];
        dv[2]=s1.part[i].vel[2]-s2.part[j].vel[2];
        vdotdr2=(dv[0]*dx[0]+dv[1]*dx[1]+dv[2]*dx[2])/dr2;
        dv2=dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2];
        mu=s1.part[i].mass+s2.part[j].mass;

#ifdef RATIMESTEP
        tau=RARVRATIO*dt_param/M_SQRT2*sqrt(dr3/mu);
        dtau=3/2.*tau*vdotdr2;
        if(dtau>1.) dtau=1.;
        tau/=(1-dtau/2);
        if(tau < timestep) timestep=tau;
#endif
#ifdef RVTIMESTEP
        if(dv2>0)
        {
          tau=dt_param*dr/sqrt(dv2);
          dtau=tau*vdotdr2*(1+mu/(dv2*dr));
          if(dtau>1.) dtau=1.;
          tau/=(1-dtau/2);
          if(tau < timestep) timestep=tau;
        }  
#endif
      }
    }
    if(timestep<s1.part[i].timestep) s1.part[i].timestep=timestep;
  }
}

/*static*/ void timestep(struct sys s1, struct sys s2)
{
#ifdef EVOLVE_OPENCL
  if((ULONG) s1.n*s2.n>CLWORKLIMIT) 
  {
    timestep_cl(s1,s2);
  } else
  {
    timestep_cpu(s1,s2);
  }
#else
  timestep_cpu(s1,s2);
#endif  
  tstep[clevel]++;
  tcount[clevel]+=(ULONG) s1.n*s2.n;
}

static void report(struct sys s,DOUBLE etime, int inttype)
{
  int maxlevel=0,i;
  long int ktot=0,dtot=0, kstot=0,dstot=0,ttot=0,tstot=0;
  UINT n,p,err=0;
  n=s.n;
  printf("** report **\n");
  printf("interaction counts:\n");
  for(i=0;i<MAXLEVEL;i++)
  {
    printf(" %4i: %10li %18li, %10li %18li\n",i, kstep[i], kcount[i], dstep[i],dcount[i]);
    if(kcount[i]>0) maxlevel=i;
    ttot+=tcount[i];
    ktot+=kcount[i];
    dtot+=dcount[i];
    tstot+=tstep[i];
    kstot+=kstep[i];
    dstot+=dstep[i];
  }    
  printf("total: %18li %18li %18li\n",ktot,dtot,ttot);  
  if(inttype == PASS_DKD || inttype == HOLD_DKD || inttype == PPASS_DKD)
    printf("equiv: %18li %18li %18li\n",(long int) deepsteps*n*n,2*deepsteps*n,(long int) deepsteps*n*n);
  else
    printf("equiv: %18li %18li %18li\n",(long int) 2*deepsteps*n*n,deepsteps*n,(long int) deepsteps*n*n);  
  printf("ksteps: %18li, dsteps: %18li, tsteps: %18li\n", kstot,dstot,tstot);
  printf("steps: %18li, equiv: %18li, maxlevel: %i\n", 
    deepsteps,((long) 1)<<maxlevel,maxlevel); 

  for(p=0;p<s.n;p++)
  {
    if(s.part[p].postime != (DOUBLE) etime) err++;
  }
  printf("postime errors: %u \n",err);
  printf("target time, actual time: %12.8g %12.8g %12.8g\n", 
           (double) etime,(double) simtime,(double) ((DOUBLE) etime-simtime));
#ifdef EVOLVE_OPENCL
  printf("cpu step,count: %12li,%18li\n",cpu_step,cpu_count);
  printf("cl step,count:  %12li,%18li\n",cl_step,cl_count);
#endif
  fflush(stdout);
}
