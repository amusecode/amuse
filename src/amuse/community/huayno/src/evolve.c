#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "evolve.h"
#include "evolve_shared.h"
#include "evolve_shared_collisions.h"
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

struct diagnostics global_diag;
struct diagnostics *diag;

static void report(struct sys s,DOUBLE etime, int inttype);

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880168872420969808L
#endif

void move_system(struct sys s, DOUBLE dpos[3],DOUBLE dvel[3],int dir)
{
  for(UINT p=0;p<s.n;p++)
  {
    for(int i=0;i<3;i++)
    {
        COMPSUMP(s.part[p].pos[i],s.part[p].pos_e[i],dir*dpos[i])
        COMPSUMV(s.part[p].vel[i],s.part[p].vel_e[i],dir*dvel[i])
    }
  }  
}

void system_center_of_mass(struct sys s, DOUBLE *cmpos, DOUBLE *cmvel)
{
  DOUBLE mass=0.,pos[3]={0.,0.,0.},vel[3]={0.,0.,0.};
  for(UINT p=0;p<s.n;p++)
  {
    for(int i=0;i<3;i++)
    {
      pos[i]+=(DOUBLE) s.part[p].mass*s.part[p].pos[i];
      vel[i]+=(DOUBLE) s.part[p].mass*s.part[p].vel[i];
    }
    mass+=(DOUBLE) s.part[p].mass;
  }
  for(int i=0;i<3;i++)
  {
    cmpos[i]=pos[i]/mass;
    cmvel[i]=vel[i]/mass;
  }
}

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
  diag=&global_diag;
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

void evolve_constant(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  if(etime == stime ||  dt==0 || clevel>=MAXLEVEL)
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel);
  diag->deepsteps++;
  diag->simtime+=dt;
  dkd(clevel,s,zerosys, stime, etime, dt);
}

void zero_diagnostics(struct diagnostics* diag)
{
  diag->deepsteps=0;
  diag->simtime=0.;
  diag->timetrack=0.;
#ifdef EVOLVE_OPENCL
  diag->cpu_step=0;
  diag->cl_step=0;
  diag->cpu_count=0;
  diag->cl_count=0;
#endif
  for(int i=0;i<MAXLEVEL;i++)
  {
    diag->tstep[i]=0;diag->tcount[i]=0;
    diag->kstep[i]=0;diag->kcount[i]=0;
    diag->dstep[i]=0;diag->dcount[i]=0;
    diag->cefail[i]=0;diag->cecount[i]=0;
    diag->bsstep[i]=0;diag->jcount[i]=0;
    diag->ntasks[i]=0;diag->taskcount[i]=0;
  }  
  diag->taskdrift=0;
  diag->taskkick=0;
}

void sum_diagnostics(struct diagnostics* total,struct diagnostics* diag)
{
  int tasksum=0;
  unsigned long taskcountsum=0;
  total->simtime+=diag->simtime;
  total->timetrack+=diag->timetrack;
  total->deepsteps+=diag->deepsteps;
  for(int i=0;i<MAXLEVEL;i++)
  {
    total->tstep[i]+=diag->tstep[i];
    total->tcount[i]+=diag->tcount[i];
    total->kstep[i]+=diag->kstep[i];
    total->kcount[i]+=diag->kcount[i];
    total->dstep[i]+=diag->dstep[i];
    total->dcount[i]+=diag->dcount[i];
    total->cefail[i]+=diag->cefail[i];
    total->cecount[i]+=diag->cecount[i];
    total->bsstep[i]+=diag->bsstep[i];
    total->jcount[i]+=diag->jcount[i];
    total->ntasks[i]+=diag->ntasks[i];tasksum+=diag->ntasks[i]; 
    total->taskcount[i]+=diag->taskcount[i];taskcountsum+=diag->taskcount[i]; 
  }          
#ifdef EVOLVE_OPENCL
  total->cpu_step+=diag->cpu_step;
  total->cl_step+=diag->cl_step;
  total->cpu_count+=diag->cpu_count;
  total->cl_count+=diag->cl_count;
#endif
#ifdef _OPENMP
  printf("%d: %d %li %li %li\n",omp_get_thread_num(),tasksum,diag->taskdrift,
     diag->taskkick,taskcountsum);
#endif
}


void do_evolve(struct sys s, double dt, int inttype)
{
  UINT p;
  int i,clevel;
  if(dt==0) return;
  for(p=0;p<s.n;p++) s.part[p].postime=0.;
  clevel=0;
  zero_diagnostics(diag);
  switch (inttype)
  {
    case CONSTANT:
      evolve_constant(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt);
      break;
    case SHARED2:
      evolve_shared2(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case SHARED4:
      evolve_shared4(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case SHARED6:
      evolve_shared6(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case SHARED8:
      evolve_shared8(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case SHARED10:
      evolve_shared10(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case SHAREDBS:
      evolve_bs_adaptive(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case BS_CC_KEPLER:
      evolve_bs(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt);
      break;
    case PASS:
      evolve_split_pass(clevel,s, zerosys,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case HOLD:
      evolve_split_hold(clevel,s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case BRIDGE:
      evolve_split_bridge(clevel,s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case NAIVE:
      evolve_split_naive(clevel,s, zerosys,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case HOLD_DKD:
      evolve_split_hold_dkd(clevel,s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case PASS_DKD:
      evolve_split_pass_dkd(clevel,s, zerosys, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case PPASS_DKD:
      evolve_split_ppass_dkd(clevel,s, zerosys, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case BRIDGE_DKD:
      evolve_split_bridge_dkd(clevel,s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case CC:
    case CC_KEPLER:
    case CC_BS:
    case CC_BSA:
    case CCC:
    case CCC_KEPLER:
    case CCC_BS:
    case CCC_BSA:
#ifdef _OPENMP
#pragma omp parallel shared(global_diag,s,dt,clevel) copyin(dt_param) 
      {
        diag=(struct diagnostics *) malloc(sizeof( struct diagnostics));
        zero_diagnostics(diag);
#pragma omp single
#endif
        evolve_cc2(clevel,s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt, inttype, 1);
#ifdef _OPENMP
#pragma omp critical
        sum_diagnostics(&global_diag,diag);
      }
      diag=&global_diag;
#endif        
      break;
    case OK:
      evolve_ok2(clevel,s, zeroforces, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case KEPLER:
      evolve_kepler(clevel,s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt);
      break;
    case FOURTH_M4:
      evolve_sf_4m4(clevel,s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case FOURTH_M5:
      evolve_sf_4m5(clevel,s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt,1);
      break;
    case SHARED6_COLLISIONS:
      evolve_shared6_nonrecursive(s, (DOUBLE) dt);
      break;
    default:  
      ENDRUN("unknown integrator\n");
      break;
  } 
  for(p=0;p<s.n;p++) s.part[p].pot=0;
  potential(s,s);
  report(s,(DOUBLE) dt, inttype);
}

void drift(int clevel,struct sys s, DOUBLE etime, DOUBLE dt)
{
  UINT i;
  for(i=0;i<s.n;i++)
  {
    COMPSUMP(s.part[i].pos[0],s.part[i].pos_e[0],dt*s.part[i].vel[0])
    COMPSUMP(s.part[i].pos[1],s.part[i].pos_e[1],dt*s.part[i].vel[1])
    COMPSUMP(s.part[i].pos[2],s.part[i].pos_e[2],dt*s.part[i].vel[2])
    s.part[i].postime=etime;
  }
  diag->dstep[clevel]++;
  diag->dcount[clevel]+=s.n;
  diag->taskdrift+=s.n;
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
    COMPSUMV(s1.part[i].vel[0],s1.part[i].vel_e[0],dt*acc[0]);
    COMPSUMV(s1.part[i].vel[1],s1.part[i].vel_e[1],dt*acc[1]);
    COMPSUMV(s1.part[i].vel[2],s1.part[i].vel_e[2],dt*acc[2]);
  }
}

void kick(int clevel,struct sys s1, struct sys s2, DOUBLE dt)
{
#ifdef EVOLVE_OPENCL
  if((ULONG) s1.n*s2.n>CLWORKLIMIT) 
  {
#pragma omp critical
    kick_cl(s1,s2,dt);
    diag->cl_step++;
    diag->cl_count+=(ULONG) s1.n*s2.n;
  } else
  {
    kick_cpu(s1,s2,dt);
    diag->cpu_step++;
    diag->cpu_count+=(ULONG) s1.n*s2.n;
  }
#else
  kick_cpu(s1,s2,dt);
#endif  
  diag->kstep[clevel]++;
  diag->kcount[clevel]+=(ULONG) s1.n*s2.n;
  diag->taskkick+=(ULONG) s1.n*s2.n;
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

void potential(struct sys s1, struct sys s2)
{
#ifdef EVOLVE_OPENCL
  if((ULONG) s1.n*s2.n>CLWORKLIMIT) 
  {
#pragma omp critical
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
 private(i,j,tau,timestep) copyin(dt_param) \
 shared(s1,s2,stdout,dir)
  for(i=0;i<s1.n;i++)
  {  
    timestep=HUGE_VAL;
    for(j=0;j<s2.n;j++)
    {
      tau=timestep_ij(s1.part+i,s2.part+j,dir);
      if(tau < timestep) timestep=tau;
    }
//    if(timestep<s1.part[i].timestep) 
    s1.part[i].timestep=timestep;
  }
}

void timestep(int clevel,struct sys s1, struct sys s2,int dir)
{
#ifdef EVOLVE_OPENCL
  if((ULONG) s1.n*s2.n>CLWORKLIMIT) 
  {
#pragma omp critical
    timestep_cl(s1,s2,dir);
  } else
  {
    timestep_cpu(s1,s2,dir);
  }
#else
  timestep_cpu(s1,s2,dir);
#endif  
  diag->tstep[clevel]++;
  diag->tcount[clevel]+=(ULONG) s1.n*s2.n;
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
    printf(" %4i: %10li %18li, %10li %18li\n",i, diag->kstep[i], diag->kcount[i], diag->dstep[i],diag->dcount[i]);
    if(diag->kcount[i]>0) maxlevel=i;
    ttot+=diag->tcount[i];
    ktot+=diag->kcount[i];
    dtot+=diag->dcount[i];
    tstot+=diag->tstep[i];
    kstot+=diag->kstep[i];
    dstot+=diag->dstep[i];
  }    
  printf("total: %18li %18li %18li\n",ktot,dtot,ttot);  
  if(inttype == PASS_DKD || inttype == HOLD_DKD || inttype == PPASS_DKD)
    printf("equiv: %18li %18li %18li\n",(long int) diag->deepsteps*n*n,2*diag->deepsteps*n,(long int) diag->deepsteps*n*n);
  else
    printf("equiv: %18li %18li %18li\n",(long int) 2*diag->deepsteps*n*n,diag->deepsteps*n,(long int) diag->deepsteps*n*n);  
  printf("ksteps: %18li, dsteps: %18li, tsteps: %18li\n", kstot,dstot,tstot);
  printf("steps: %18li, equiv: %18li, maxlevel: %i\n", 
    diag->deepsteps,((long) 1)<<maxlevel,maxlevel); 

  for(p=0;p<s.n;p++)
  {
    if(s.part[p].postime != (DOUBLE) etime) err++;
  }
  printf("postime errors: %u \n",err);
  printf("target time, actual time: %12.8g %12.8g %12.8g\n", 
           (double) etime,(double) diag->simtime,(double) ((DOUBLE) etime-diag->simtime));
  printf("time track, ratio: %12.8g %12.8g\n", (double) diag->timetrack,(double) (diag->timetrack/diag->simtime));

#ifdef EVOLVE_OPENCL
  printf("cpu step,count: %12li,%18li\n",diag->cpu_step,diag->cpu_count);
  printf("cl step,count:  %12li,%18li\n",diag->cl_step,diag->cl_count);
#endif
  if(inttype==SHAREDBS || inttype==CC_BS || inttype==CCC_BS || inttype==CC_BSA || inttype==CCC_BSA)
  {
    unsigned long totalbs=0,totalj=0;
    printf("bs counts:\n");
    for(i=0;i<MAXLEVEL;i++) 
    { 
      totalbs+=diag->bsstep[i];
      totalj+=diag->jcount[i];
      printf("%d: %18li %18li %f\n",i,diag->bsstep[i],diag->jcount[i],diag->jcount[i]/(1.*diag->bsstep[i]+1.e-20));
    }
    printf(" total, total j, mean j: %18li %18li %f\n",totalbs,totalj,totalj/(1.*totalbs));
  }
  if(inttype==CC_KEPLER || inttype==CCC_KEPLER)
  {
    unsigned long totalcefail=0,totalcecount=0;
    printf("kepler solver counts:\n");
    for(i=0;i<MAXLEVEL;i++) 
    { 
      totalcefail+=diag->cefail[i];
      totalcecount+=diag->cecount[i];
      printf("%d: %18li %18li\n",i,diag->cefail[i],diag->cecount[i]);
    } 
    printf(" total, total j, mean j: %18li %18li\n",totalcefail,totalcecount);
  }
  
#ifdef _OPENMP
  {
    int totaltasks=0;
    printf("task counts:\n");
    for(i=0;i<MAXLEVEL;i++) 
    { 
      printf("%d: %18li %18li\n",i,diag->ntasks[i],diag->taskcount[i]);
      totaltasks+=diag->ntasks[i];
    } 
   printf("openmp tasks: %d\n",totaltasks);

  }

#endif  
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

void kdk(int clevel,struct sys s1,struct sys s2, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  if(s2.n>0) kick(clevel,s2, s1, dt/2);
  kick(clevel,s1,join(s1,s2),dt/2);
  drift(clevel,s1,etime, dt);
  kick(clevel,s1,join(s1,s2),dt/2);
  if(s2.n>0) kick(clevel,s2, s1, dt/2);
}

void dkd(int clevel,struct sys s1,struct sys s2, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  drift(clevel,s1,stime+dt/2, dt/2);
  kick(clevel,s1,join(s1,s2),dt);
  if(s2.n>0) kick(clevel,s2, s1, dt);
  drift(clevel,s1,etime, dt/2);
}
