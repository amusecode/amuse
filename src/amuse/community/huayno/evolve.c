#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "evolve.h"
#include "evolve_shared.h"
#include "evolve_sf.h"
#include "evolve_cc.h"
#include "evolve_kepler.h"
#include "evolve_ok.h"
#include "evolve_bs.h"

#ifdef EVOLVE_OPENCL
#include "evolve_cl.h"
#endif

FLOAT eps2;
FLOAT dt_param;
struct sys zerosys ={ 0, NULL,NULL};

/* diagnostics */
DOUBLE simtime;
DOUBLE timetrack;
int clevel;
unsigned long tcount[MAXLEVEL],kcount[MAXLEVEL],dcount[MAXLEVEL],deepsteps;
unsigned long tstep[MAXLEVEL],kstep[MAXLEVEL],dstep[MAXLEVEL];
#ifdef EVOLVE_OPENCL
unsigned long cpu_step,cl_step,cpu_count,cl_count;
#endif

static void potential(struct sys s1, struct sys s2);
static void report(struct sys s,DOUBLE etime, int inttype);

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880168872420969808L
#endif

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

void init_evolve(struct sys s,int inttype)
{
  UINT i;
  for(i=0;i<s.n;i++)
  {
    s.part[i].postime=0.;
    s.part[i].pot=0.;
    s.part[i].timestep=HUGE_VAL;
#ifdef COMPENSATED_SUMMP
    s.part[i].pos_e[0]=0.;s.part[i].pos_e[1]=0.;s.part[i].pos_e[2]=0.;
#endif
#ifdef COMPENSATED_SUMMV
    s.part[i].vel_e[0]=0.;s.part[i].vel_e[1]=0.;s.part[i].vel_e[2]=0.;
#endif
  }
  potential(s,s);
  
  evolve_ok_stop();
  if (inttype == OK) evolve_ok_init(s);
}

void evolve_constant(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  clevel++;
  if(etime == stime ||  dt==0 || clevel>=MAXLEVEL)
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel);
  deepsteps++;
  simtime+=dt;
  dkd(s,zerosys, stime, etime, dt);
  clevel--;
}


void do_evolve(struct sys s, double dt, int inttype)
{
  UINT p;
  int i;
  if(dt==0) return;
  for(p=0;p<s.n;p++) s.part[p].postime=0.;
  clevel=-1;
  deepsteps=0;
  simtime=0.;
  timetrack=0.;
  for(i=0;i<MAXLEVEL;i++)
  {
    tstep[i]=0;tcount[i]=0;
    kstep[i]=0;kcount[i]=0;
    dstep[i]=0;dcount[i]=0;
    cefail[i]=0;cecount[i]=0;
    bsstep[i]=0;jcount[i]=0;
  }
  switch (inttype)
  {
    case CONSTANT:
      evolve_constant(s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt);
      break;
    case SHARED2:
      evolve_shared2(s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case SHARED4:
      evolve_shared4(s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case SHARED6:
      evolve_shared6(s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case SHARED8:
      evolve_shared8(s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case SHARED10:
      evolve_shared10(s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case SHAREDBS:
      evolve_bs_adaptive(s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
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
    case FOURTH_M4:
      evolve_sf_4m4(s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case FOURTH_M5:
      evolve_sf_4m5(s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    default:  
      ENDRUN("unknown integrator\n");
      break;
  } 
  for(p=0;p<s.n;p++) s.part[p].pot=0;
  potential(s,s);
  report(s,(DOUBLE) dt, inttype);
}

#define COMPSUM(sum,err,delta) \
  { \
    DOUBLE a; \
    a=sum; \
    err=err+delta; \
    sum=a+err; \
    err=err+(a-sum); \
  }

#define COMPSUM1(sum,err,delta) \
  { \
    DOUBLE t,y; \
    y=(delta)-err; \
    t=sum+y; \
    err=(t-sum)-y; \
    sum=t; \
  }


void drift(struct sys s, DOUBLE etime, DOUBLE dt)
{
  UINT i;
  for(i=0;i<s.n;i++)
  {
#ifndef COMPENSATED_SUMMP
    s.part[i].pos[0]+=dt*s.part[i].vel[0];
    s.part[i].pos[1]+=dt*s.part[i].vel[1];
    s.part[i].pos[2]+=dt*s.part[i].vel[2];
#else
    COMPSUM(s.part[i].pos[0],s.part[i].pos_e[0],dt*s.part[i].vel[0])
    COMPSUM(s.part[i].pos[1],s.part[i].pos_e[1],dt*s.part[i].vel[1])
    COMPSUM(s.part[i].pos[2],s.part[i].pos_e[2],dt*s.part[i].vel[2])
#endif
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

#pragma omp parallel for if((ULONG) s1.n*s2.n>MPWORKLIMIT && !omp_in_parallel()) default(none) \
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
#ifndef COMPENSATED_SUMMV
    s1.part[i].vel[0]+=dt*acc[0];
    s1.part[i].vel[1]+=dt*acc[1];
    s1.part[i].vel[2]+=dt*acc[2];
#else
    COMPSUM(s1.part[i].vel[0],s1.part[i].vel_e[0],dt*acc[0]);
    COMPSUM(s1.part[i].vel[1],s1.part[i].vel_e[1],dt*acc[1]);
    COMPSUM(s1.part[i].vel[2],s1.part[i].vel_e[2],dt*acc[2]);
#endif
  }
}

void kick(struct sys s1, struct sys s2, DOUBLE dt)
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

#pragma omp parallel for if((ULONG) s1.n*s2.n>MPWORKLIMIT && !omp_in_parallel()) default(none) \
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

inline FLOAT timestep_ij(struct particle *i, struct particle *j,int dir) {
  FLOAT timestep;
  FLOAT dx[3],dr3,dr2,dr,dv[3],dv2,mu,vdotdr2,tau,dtau;
  timestep=HUGE_VAL;
  if(i==j) return timestep;
  dx[0]=i->pos[0] - j->pos[0];
  dx[1]=i->pos[1] - j->pos[1];
  dx[2]=i->pos[2] - j->pos[2];
  dr2=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+eps2;
  if(dr2>0) 
  {
    dr=sqrt(dr2);
    dr3=dr*dr2;
    dv[0]=i->vel[0] - j->vel[0];
    dv[1]=i->vel[1] - j->vel[1];
    dv[2]=i->vel[2] - j->vel[2];
    vdotdr2=(dv[0]*dx[0]+dv[1]*dx[1]+dv[2]*dx[2])/dr2;
    dv2=dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2];
    mu=i->mass + j->mass;
#ifdef RATIMESTEP
    tau=RARVRATIO*dt_param/M_SQRT2*sqrt(dr3/mu);
    dtau=3/2.*dir*tau*vdotdr2;
    if(dtau>1.) dtau=1.;
    tau/=(1-dtau/2);
    if(tau < timestep) timestep=tau;
#endif
#ifdef RVTIMESTEP
    if(dv2>0) 
    {
      tau=dt_param*dr/sqrt(dv2);
      dtau=dir*tau*vdotdr2*(1+mu/(dv2*dr));
      if(dtau>1.) dtau=1.;
      tau/=(1-dtau/2);
      if(tau < timestep) timestep=tau;
    }
#endif
  }
  if (timestep < 0) 
  {
    ENDRUN("negative timestep!\n");
  }
  return timestep;
}


static void timestep_cpu(struct sys s1, struct sys s2,int dir)
{
  UINT i,j;
  FLOAT timestep,tau;
#pragma omp parallel for if((ULONG) s1.n*s2.n>MPWORKLIMIT && !omp_in_parallel()) default(none) \
 private(i,j,tau,timestep) \
 shared(s1,s2,stdout,dir)
  for(i=0;i<s1.n;i++)
  {
    if(s1.part[i].timestep !=HUGE_VAL) ENDRUN("timestep??");
  
    timestep=HUGE_VAL;
    for(j=0;j<s2.n;j++)
    {
      tau=timestep_ij(s1.part+i,s2.part+j,dir);
      if(tau < timestep) timestep=tau;
    }
    if(timestep<s1.part[i].timestep) s1.part[i].timestep=timestep;
  }
}

void timestep(struct sys s1, struct sys s2,int dir)
{
#ifdef EVOLVE_OPENCL
  if((ULONG) s1.n*s2.n>CLWORKLIMIT) 
  {
    timestep_cl(s1,s2,dir);
  } else
  {
    timestep_cpu(s1,s2,dir);
  }
#else
  timestep_cpu(s1,s2,dir);
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
  printf("time track, ratio: %12.8g %12.8g\n", (double) timetrack,(double) (timetrack/simtime));

#ifdef EVOLVE_OPENCL
  printf("cpu step,count: %12li,%18li\n",cpu_step,cpu_count);
  printf("cl step,count:  %12li,%18li\n",cl_step,cl_count);
#endif
  if(inttype==SHAREDBS)
  {
    unsigned long totalbs,totalj;
    printf("bs counts:\n");
    for(i=0;i<MAXLEVEL;i++) 
    { 
      totalbs+=bsstep[i];
      totalj+=jcount[i];
      printf("%d: %18li %18li %f\n",i,bsstep[i],jcount[i],jcount[i]/(1.*bsstep[i]+1.e-20));
    }
    printf(" total, total j, mean j: %18li %18li %f",totalbs,totalj,totalj/(1.*totalbs));
  }
  fflush(stdout);
}

struct sys join(struct sys s1,struct sys s2)
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

FLOAT global_timestep(struct sys s)
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

void kdk(struct sys s1,struct sys s2, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  if(s2.n>0) kick(s2, s1, dt/2);
  kick(s1,join(s1,s2),dt/2);
  drift(s1,etime, dt);
  kick(s1,join(s1,s2),dt/2);
  if(s2.n>0) kick(s2, s1, dt/2);
}

void dkd(struct sys s1,struct sys s2, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  drift(s1,stime+dt/2, dt/2);
  kick(s1,join(s1,s2),dt);
  if(s2.n>0) kick(s2, s1, dt);
  drift(s1,etime, dt/2);
}
