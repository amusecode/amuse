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
#include "evolve_error_control.h"

#ifdef EVOLVE_OPENCL
#include "evolve_cl.h"
#endif

int verbosity=0;

struct sys debugsys;

FLOAT eps2;
FLOAT dt_param;
struct sys zerosys ={0,0,NULL,NULL};
int accel_zero_mass=1;
int opencl_device_type=0;

struct diagnostics global_diag;
struct diagnostics *diag;

static void report(struct sys s,DOUBLE etime, int inttype);

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880168872420969808L
#endif

void move_system(struct sys s, DOUBLE dpos[3],DOUBLE dvel[3],int dir)
{
  struct particle *ipart;
  for(UINT p=0;p<s.n;p++)
  {
    ipart=GETPART(s, p);
    for(int i=0;i<3;i++)
    {
        COMPSUMP(ipart->pos[i],ipart->pos_e[i],dir*dpos[i])
        COMPSUMV(ipart->vel[i],ipart->vel_e[i],dir*dvel[i])
    }
  }  
}

void system_center_of_mass(struct sys s, DOUBLE *cmpos, DOUBLE *cmvel)
{
  DOUBLE mass=0.,pos[3]={0.,0.,0.},vel[3]={0.,0.,0.};
  struct particle *ipart;
  for(UINT p=0;p<s.n;p++)
  {
    ipart=GETPART(s, p);
    for(int i=0;i<3;i++)
    {
      pos[i]+=(DOUBLE) ipart->mass*ipart->pos[i];
      vel[i]+=(DOUBLE) ipart->mass*ipart->vel[i];
    }
    mass+=(DOUBLE) ipart->mass;
  }
  for(int i=0;i<3;i++)
  {
    cmpos[i]=pos[i]/mass;
    cmvel[i]=vel[i]/mass;
  }
}

DOUBLE system_kinetic_energy(struct sys s)
{
 UINT i;
 DOUBLE e=0.;
 struct particle *ipart;
 for(i=0;i<s.n;i++)
 {
   ipart=GETPART(s, i);
   e+=0.5*ipart->mass*( ipart->vel[0]*ipart->vel[0]+
                        ipart->vel[1]*ipart->vel[1]+
                        ipart->vel[2]*ipart->vel[2] );
 }
 return e;
}
        
DOUBLE system_potential_energy(struct sys s)
{
 UINT i;
 DOUBLE e=0.;
 struct particle *ipart;
 for(i=0;i<s.n;i++)
 {
   ipart=GETPART(s, i);
   e+=ipart->mass*ipart->pot;
 }
 return e/2;
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
  struct particle *ipart;
  for(UINT i=0;i<s.n;i++)
  {
    ipart=GETPART(s,i);
    ipart->postime=0.;
    ipart->pot=0.;
#ifdef COMPENSATED_SUMMP
    ipart->pos_e[0]=0.;ipart->pos_e[1]=0.;ipart->pos_e[2]=0.;
#endif
#ifdef COMPENSATED_SUMMV
    ipart->vel_e[0]=0.;ipart->vel_e[1]=0.;ipart->vel_e[2]=0.;
#endif
  }
  if(accel_zero_mass) split_zeromass(&s); // because of potential calc
  potential(s,s);
  
  evolve_ok_stop();
  if (inttype == OK) evolve_ok_init(s);
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
  if(verbosity>0) printf("task %d: %d %li %li %li\n",omp_get_thread_num(),tasksum,diag->taskdrift, diag->taskkick,taskcountsum);
#endif
}


void do_evolve(struct sys s, double dt, int inttype)
{
  int i,clevel;
  struct particle *ipart;
  if(dt==0) return;
  for(UINT p=0;p<s.n;p++) GETPART(s,p)->postime=0.;
  clevel=0;
  if(accel_zero_mass) split_zeromass(&s);
  zero_diagnostics(diag);
  debugsys=s;
  switch (inttype)
  {
    case CONSTANT:
    case CONSTANT2:
      evolve_constant2(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt);
      break;
    case CONSTANT4:
      evolve_constant4(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt);
      break;
    case CONSTANT6:
      evolve_constant6(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt);
      break;
    case CONSTANT8:
      evolve_constant8(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt);
      break;
    case CONSTANT10:
      evolve_constant10(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt);
      break;
    case SHARED2:
      evolve_shared2(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt, -1.);
      break;
    case SHARED4:
      evolve_shared4(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt, -1.);
      break;
    case SHARED6:
      evolve_shared6(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt, -1.);
      break;
    case SHARED8:
      evolve_shared8(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt, -1.);
      break;
    case SHARED10:
      evolve_shared10(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt, -1.);
      break;
    case SHAREDBS:
      evolve_bs_adaptive(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt, -1.);
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
    case CC_SHARED10:
    case CCC_SHARED10:
#ifdef _OPENMP
#pragma omp parallel shared(global_diag,s,dt,clevel) copyin(dt_param) 
      {
        diag=(struct diagnostics *) malloc(sizeof( struct diagnostics));
        zero_diagnostics(diag);
#pragma omp master      
	      if(verbosity>0) printf("Total Threads # %d\n", omp_get_num_threads()); 
#pragma omp single
#endif
#ifdef CC2_SPLIT_SHORTCUTS
        evolve_cc2_shortcut(clevel,s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt, inttype, 1, -1.);
#else
        evolve_cc2(clevel,s,(DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt, inttype, 1);
#endif        
#ifdef _OPENMP
#pragma omp critical
        sum_diagnostics(&global_diag,diag);
        free(diag);
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
    case SHARED2_COLLISIONS:
      evolve_shared2_collision_detection(s, (DOUBLE) dt);
      break;
    case SHARED4_COLLISIONS:
      evolve_shared4_collision_detection(s, (DOUBLE) dt);
      break;
    case SHARED6_COLLISIONS:
      evolve_shared6_collision_detection(s, (DOUBLE) dt);
      break;
    case SHARED8_COLLISIONS:
      evolve_shared8_collision_detection(s, (DOUBLE) dt);
      break;
    case SHARED10_COLLISIONS:
      evolve_shared10_collision_detection(s, (DOUBLE) dt);
      break;
    case ERROR_CONTROL:
      evolve_error_control(clevel,s, (DOUBLE) 0.,(DOUBLE) dt,(DOUBLE) dt, -1.);
      break;
    default:  
      ENDRUN("unknown integrator %d\n", inttype);
      break;
  } 
  for(UINT p=0;p<s.n;p++) GETPART(s,p)->pot=0;
  potential(s,s);
  if(verbosity>0) report(s,(DOUBLE) dt, inttype);
}

void drift(int clevel,struct sys s, DOUBLE etime, DOUBLE dt)
{
  struct particle *ipart;
  for(UINT i=0;i<s.n;i++)
  {
    ipart=GETPART(s,i);
    COMPSUMP(ipart->pos[0],ipart->pos_e[0],dt*ipart->vel[0])
    COMPSUMP(ipart->pos[1],ipart->pos_e[1],dt*ipart->vel[1])
    COMPSUMP(ipart->pos[2],ipart->pos_e[2],dt*ipart->vel[2])
    ipart->postime=etime;
  }
  diag->dstep[clevel]++;
  diag->dcount[clevel]+=s.n;
  diag->taskdrift+=s.n;
}

static void kick_cpu(struct sys s1, struct sys s2, DOUBLE dt)
{
  FLOAT dx[3],dr3,dr2,dr,acci;
  FLOAT acc[3];
  struct particle *ipart, *jpart;

#pragma omp parallel for if((ULONG) s1.n*(s2.n-s2.nzero)>MPWORKLIMIT && !omp_in_parallel()) default(none) \
 private(dx,dr3,dr2,dr,acc,acci,ipart,jpart) \
 shared(dt,s1,s2,eps2)
  for(UINT i=0;i<s1.n;i++)
  {
    ipart=GETPART(s1,i);
    acc[0]=0.;
    acc[1]=0.;
    acc[2]=0.;
    for(UINT j=0;j<s2.n-s2.nzero;j++)
    {
      jpart=GETPART(s2,j);
//~ if(jpart->mass==0) continue;
//      if(ipart==jpart) continue; 
      dx[0]=ipart->pos[0]-jpart->pos[0];
      dx[1]=ipart->pos[1]-jpart->pos[1];
      dx[2]=ipart->pos[2]-jpart->pos[2];
      dr2=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+eps2;
      if(dr2>0) 
      {
        dr=sqrt(dr2);
        dr3=dr*dr2;
        acci=jpart->mass/dr3;

        acc[0]-=dx[0]*acci;
        acc[1]-=dx[1]*acci;
        acc[2]-=dx[2]*acci;  
      }
    }
    COMPSUMV(ipart->vel[0],ipart->vel_e[0],dt*acc[0]);
    COMPSUMV(ipart->vel[1],ipart->vel_e[1],dt*acc[1]);
    COMPSUMV(ipart->vel[2],ipart->vel_e[2],dt*acc[2]);
  }
}

void kick(int clevel,struct sys s1, struct sys s2, DOUBLE dt)
{
#ifdef EVOLVE_OPENCL
  if((ULONG) s1.n*(s2.n-s2.nzero)>CLWORKLIMIT) 
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
  FLOAT dx[3],dr2,dr;
  FLOAT pot;
  struct particle *ipart, *jpart;

#pragma omp parallel for if((ULONG) s1.n*(s2.n-s2.nzero)>MPWORKLIMIT && !omp_in_parallel()) default(none) \
 private(dx,dr2,dr,pot, ipart,jpart) \
 shared(s1,s2,eps2)
  for(UINT i=0;i<s1.n;i++)
  {
    pot=0;
    ipart=GETPART(s1,i);
    for(UINT j=0;j<s2.n-s2.nzero;j++)
    {
      jpart=GETPART(s2,j);
      if(ipart==jpart) continue; 
      dx[0]=ipart->pos[0]-jpart->pos[0];
      dx[1]=ipart->pos[1]-jpart->pos[1];
      dx[2]=ipart->pos[2]-jpart->pos[2];
      dr2=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+eps2;
      if(dr2>0) 
      {
        dr=sqrt(dr2);
        pot-=jpart->mass/dr;
      }
    }
    ipart->pot+=pot;
  }
}

void potential(struct sys s1, struct sys s2)
{
#ifdef EVOLVE_OPENCL
  if((ULONG) s1.n*(s2.n-s2.nzero)>CLWORKLIMIT) 
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
  mu=i->mass + j->mass;
  if(dr2>0 && mu>0) 
  {
    dr=sqrt(dr2);
    dr3=dr*dr2;
    dv[0]=i->vel[0] - j->vel[0];
    dv[1]=i->vel[1] - j->vel[1];
    dv[2]=i->vel[2] - j->vel[2];
    vdotdr2=(dv[0]*dx[0]+dv[1]*dx[1]+dv[2]*dx[2])/dr2;
    dv2=dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2];
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
  UINT i,j, jmax;
  FLOAT timestep,tau;
  struct particle *ipart;
#pragma omp parallel for if((ULONG) (s1.n*s2.n-s1.nzero*s2.nzero)>MPWORKLIMIT && !omp_in_parallel()) default(none) \
 private(i,j,tau,timestep, jmax, ipart) copyin(dt_param) \
 shared(s1,s2,stdout,dir)
  for(i=0;i<s1.n;i++)
  {  
    timestep=HUGE_VAL;
    ipart=GETPART(s1,i);
    jmax=s2.n;if(i>=s1.n-s1.nzero) jmax=s2.n-s2.nzero;
    for(j=0;j<jmax;j++)
    {
      tau=timestep_ij(ipart,GETPART(s2,j),dir);
      if(tau < timestep) timestep=tau;
    }
//    if(timestep<ipart->timestep) 
    ipart->timestep=timestep;
  }
}

void timestep(int clevel,struct sys s1, struct sys s2,int dir)
{
#ifdef EVOLVE_OPENCL
  if((ULONG) (s1.n*s2.n-s1.nzero*s2.nzero)>CLWORKLIMIT) 
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
  struct particle *ipart;
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

  for(p=0;p<s.n;p++) if(GETPART(s,p)->postime != (DOUBLE) etime) err++;
  printf("postime errors: %u \n",err);
  printf("target time, actual time: %12.8g %12.8g %12.8g\n", 
           (double) etime,(double) diag->simtime,(double) ((DOUBLE) etime-diag->simtime));
  printf("time track, ratio: %12.8g %12.8g\n", (double) diag->timetrack,
       (double) (diag->simtime!=0? (diag->timetrack/diag->simtime) :1));

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
  if(inttype==KEPLER || inttype==CC_KEPLER || inttype==CCC_KEPLER || inttype==CCC_BS || 
     inttype==CC_BS || inttype==CCC_BSA || inttype==CC_BSA || inttype==CC_SHARED10 || inttype==CC_SHARED10)
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

void join_array(UINT n1, struct particle *p1,
                UINT n2, struct particle *p2,
                UINT *n, struct particle **p)
{
  if(n1==0 && n2==0)
  {
    *n=0;*p=NULL;
  }
  if(n1!=0 && n2==0)
  {
    *n=n1;*p=p1;
  }
  if(n1==0 && n2!=0)
  {
    *n=n2;*p=p2;
  }
  if(n1!=0 && n2!=0)
  {
    *n=n1+n2;
    if(p1+n1==p2)
    {
      *p=p1;
    }
    else
    {
      if(p2+n2==p1)
      {
        *p=p2;
      } else ENDRUN("join_array error");
    }
  }
}                  

struct sys join(struct sys s1,struct sys s2)
{
  struct sys s=zerosys;
  if(s1.n==0) return s2;
  if(s2.n==0) return s1;
  join_array(s1.n-s1.nzero, s1.part, s2.n-s2.nzero, s2.part, &s.n, &s.part);
  join_array(s1.nzero, s1.zeropart, s2.nzero, s2.zeropart, &s.nzero, &s.zeropart);
  s.n=s.n+s.nzero;
  if(s.n-s.nzero>0 && LAST(s)-s.part + 1 != s.n-s.nzero) ENDRUN("join error 1");
  if(s.nzero>0 && LASTZERO(s)-s.zeropart + 1 != s.nzero) ENDRUN("join error 2");
  return s;
}

FLOAT global_timestep(struct sys s)
{
  FLOAT dt,mindt=HUGE_VAL;
  for(UINT i=0;i<s.n;i++)
  {
    dt=GETPART(s, i)->timestep;
    if(dt<mindt) mindt=dt;
  }
  return mindt;
}

FLOAT max_global_timestep(struct sys s)
{
  FLOAT dt,maxdt=0;
  for(UINT i=0;i<s.n;i++)
  {
    dt=GETPART(s, i)->timestep;
    if(dt>maxdt) maxdt=dt;
  }
  return maxdt;
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

void split_zeromass(struct sys *s)
{
  UINT i=0;
  struct particle *left, *right;
  if(s->n==0) return;
  if(s->part==NULL) ENDRUN("split_zeromass malformed input");
  if(s->n-s->nzero==0)
  {
    if(s->zeropart==NULL || s->part!=s->zeropart) ENDRUN("split_zeromass malformed input");
    if(LASTZERO(*s)-s->zeropart+1!=s->nzero) ENDRUN( "split_zeromass malformed input sys");
    return;
  }  
  if(s->nzero!=0 && LAST(*s)+1!=s->zeropart) 
    ENDRUN("split_zeromass can only work on fully contiguous sys");
  left=s->part;
  right=s->part+(s->n-1);
  while(1)
  {
    if(i>=s->n) ENDRUN("split_zeromass error 1");
    i++;
    while(left->mass!=0 && left<right) left++;
    while(right->mass==0 && left<right) right--;
    if(left<right) 
      {SWAP( *left, *right, struct particle);}
    else 
      break;  
  }
  if(left->mass!=0) left++;
  s->nzero=s->n-(left-s->part);
  if(s->nzero<0) ENDRUN("split_zeromass find negative number of part");
  if(s->nzero>0)
  {
    s->zeropart=left;
  }
  if((left-s->part)+s->nzero !=s->n) ENDRUN( "split_zeromass error 2");
  for(i=0;i<(s->n-s->nzero);i++) if(GETPART(*s,i)->mass==0) ENDRUN ("split_zromass error 3");
  for(i=s->n-s->nzero;i<s->n;i++) if(GETPART(*s,i)->mass!=0) ENDRUN ("split_zeromass error 4");
#ifdef CONSISTENCY_CHECKS
  verify_split_zeromass(*s);
#endif
}

void verify_split_zeromass(struct sys s)
{
  if(!accel_zero_mass) return;
  for(UINT i=0;i<s.n-s.nzero;i++) if(GETPART(s,i)->mass==0) ENDRUN("massless particle in main part\n") 
  for(UINT i=s.n-s.nzero;i<s.n;i++) if(GETPART(s,i)->mass!=0) ENDRUN("massive particle in massless part\n") 
}
