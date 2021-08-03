#include "evolve.h"
#ifdef EVOLVE_OPENCL

#include <tgmath.h>
#include <stdcl.h>
#include "evolve_cl.h"
#include "evolve_kern.clh"

#define cl_double4_x(f) (((cl_double*)&(f))[0])
#define cl_double4_y(f) (((cl_double*)&(f))[1])
#define cl_double4_z(f) (((cl_double*)&(f))[2])
#define cl_double4_w(f) (((cl_double*)&(f))[3])

static void *h;
static cl_kernel kick_krn;
static cl_kernel timestep_krn;
static cl_kernel potential_krn;


void init_cl()
{
  h = clsopen(CLCONTEXT,srcstr,CLLD_NOW);
  kick_krn = clsym(CLCONTEXT,h,"kick_kernel",CLLD_NOW);
  timestep_krn = clsym(CLCONTEXT,h,"timestep_kernel",CLLD_NOW);
  potential_krn = clsym(CLCONTEXT,h,"potential_kernel",CLLD_NOW);
}

void close_cl()
{
  clclose(CLCONTEXT,h);
}

void kick_cl(struct sys s1, struct sys s2, DOUBLE dt)
{
  int i, n2;
  int groupsize,nthread,blocksize;
  CLFLOAT cleps2=(CLFLOAT) eps2;
  struct particle *ipart;
  
  n2=s2.n-s2.nzero;
  if(s1.n==0 || n2==0) return;

  groupsize=NTHREAD;
  if(s1.n < NTHREAD) groupsize=s1.n;
  nthread=((s1.n-1)/groupsize+1)*groupsize;
  blocksize=BLOCKSIZE;
  if(n2<blocksize) blocksize=n2;
  clndrange_t ndr = clndrange_init1d(0,nthread,groupsize);

  CLFLOAT4* ipos = (CLFLOAT4*) clmalloc(CLCONTEXT,nthread*sizeof(CLFLOAT4),0);
  CLFLOAT4* acc = (CLFLOAT4*) clmalloc(CLCONTEXT,nthread*sizeof(CLFLOAT4),0);
  CLFLOAT4* jpos = (CLFLOAT4*) clmalloc(CLCONTEXT,n2*sizeof(CLFLOAT4),0);

  for(i=0;i<s1.n;i++)
  {
    ipart=GETPART(s1,i);
    ipos[i].s0=ipart->pos[0];
    ipos[i].s1=ipart->pos[1];
    ipos[i].s2=ipart->pos[2];
    ipos[i].s3=ipart->mass;
  }  
  for(i=s1.n;i<nthread;i++) ipos[i]=(CLFLOAT4) {{0.0,0.0,0.0,0.0}};  
  for(i=0;i<n2;i++)
  {
    ipart=GETPART(s2,i);
    jpos[i].s0=ipart->pos[0];
    jpos[i].s1=ipart->pos[1];
    jpos[i].s2=ipart->pos[2];
    jpos[i].s3=ipart->mass;
  }
  for(i=0;i<nthread;i++) acc[i]=(CLFLOAT4) {{0.0,0.0,0.0,0.0}};  
                        
  clarg_set(CLCONTEXT,kick_krn,0, n2);
  clarg_set(CLCONTEXT,kick_krn,1, blocksize);
  clarg_set(CLCONTEXT,kick_krn,2, cleps2);
  clarg_set_global(CLCONTEXT,kick_krn,3, ipos);
  clarg_set_global(CLCONTEXT,kick_krn,4, jpos);
  clarg_set_global(CLCONTEXT,kick_krn,5, acc);
  clarg_set_local(CLCONTEXT,kick_krn,6,blocksize*sizeof(CLFLOAT4));

  clmsync(CLCONTEXT,0,ipos,CL_MEM_DEVICE|CL_EVENT_NOWAIT);
  clmsync(CLCONTEXT,0,jpos,CL_MEM_DEVICE|CL_EVENT_NOWAIT);
  clfork(CLCONTEXT,0,kick_krn,&ndr,CL_EVENT_NOWAIT);
  clmsync(CLCONTEXT,0,acc,CL_MEM_HOST|CL_EVENT_NOWAIT);
  clwait(CLCONTEXT,0,CL_KERNEL_EVENT|CL_MEM_EVENT);

  for(i=0;i<s1.n;i++)
  {
    ipart=GETPART(s1,i);
    ipart->vel[0]+=dt*acc[i].s0;
    ipart->vel[1]+=dt*acc[i].s1;
    ipart->vel[2]+=dt*acc[i].s2;
  }

  clfree( ipos);
  clfree( jpos);
  clfree( acc);  
}

void _timestep_cl(struct sys s1, struct sys s2,int dir)
{
  int i,n2;
  int groupsize,nthread,blocksize;
  CLFLOAT cleps2=(CLFLOAT) eps2;
  CLFLOAT cldtparam=(CLFLOAT) dt_param;
  struct particle *ipart;
    
  n2=s2.n;  
  if(s1.n==0 || s2.n==0) return;

  groupsize=NTHREAD;
  if(s1.n < NTHREAD) groupsize=s1.n;
  nthread=((s1.n-1)/groupsize+1)*groupsize;
  blocksize=BLOCKSIZE/2; 
  if(s2.n<blocksize) blocksize=s2.n;
  clndrange_t ndr = clndrange_init1d(0,nthread,groupsize);

  CLFLOAT4* ipos = (CLFLOAT4*) clmalloc(CLCONTEXT,nthread*sizeof(CLFLOAT4),0);
  CLFLOAT4* ivel = (CLFLOAT4*) clmalloc(CLCONTEXT,nthread*sizeof(CLFLOAT4),0);
  CLFLOAT* timestep = (CLFLOAT*) clmalloc(CLCONTEXT,nthread*sizeof(CLFLOAT),0);
  CLFLOAT4* jpos = (CLFLOAT4*) clmalloc(CLCONTEXT,s2.n*sizeof(CLFLOAT4),0);
  CLFLOAT4* jvel = (CLFLOAT4*) clmalloc(CLCONTEXT,s2.n*sizeof(CLFLOAT4),0);

  for(i=0;i<s1.n;i++)
  {
    ipart=GETPART(s1,i);
    ipos[i].s0=ipart->pos[0]; ivel[i].s0=ipart->vel[0];
    ipos[i].s1=ipart->pos[1]; ivel[i].s1=ipart->vel[1];
    ipos[i].s2=ipart->pos[2]; ivel[i].s2=ipart->vel[2];
    ipos[i].s3=ipart->mass;   ivel[i].s3=0.;
  }  
  for(i=s1.n;i<nthread;i++) ipos[i]=(CLFLOAT4) {{0.0,0.0,0.0,0.0}};  
  for(i=s1.n;i<nthread;i++) ivel[i]=(CLFLOAT4) {{0.0,0.0,0.0,0.0}};  
  for(i=0;i<s2.n;i++)
  {
    ipart=GETPART(s2,i);
    jpos[i].s0= ipart->pos[0]; jvel[i].s0=ipart->vel[0];
    jpos[i].s1= ipart->pos[1]; jvel[i].s1=ipart->vel[1];
    jpos[i].s2= ipart->pos[2]; jvel[i].s2=ipart->vel[2];
    jpos[i].s3= ipart->mass;   jvel[i].s3=0.;
  }
  for(i=0;i<nthread;i++) timestep[i]=(CLFLOAT) 0.;  
                        
  clarg_set(CLCONTEXT,timestep_krn,0, n2);
  clarg_set(CLCONTEXT,timestep_krn,1, blocksize);
  clarg_set(CLCONTEXT,timestep_krn,2, cleps2);
  clarg_set(CLCONTEXT,timestep_krn,3, cldtparam);
  clarg_set_global(CLCONTEXT,timestep_krn,4, ipos);
  clarg_set_global(CLCONTEXT,timestep_krn,5, ivel);
  clarg_set_global(CLCONTEXT,timestep_krn,6, jpos);
  clarg_set_global(CLCONTEXT,timestep_krn,7, jvel);
  clarg_set_global(CLCONTEXT,timestep_krn,8, timestep);
  clarg_set_local(CLCONTEXT,timestep_krn,9,blocksize*sizeof(CLFLOAT4));
  clarg_set_local(CLCONTEXT,timestep_krn,10,blocksize*sizeof(CLFLOAT4));
  clarg_set(CLCONTEXT,timestep_krn,11, dir);
  
  clmsync(CLCONTEXT,0,ipos,CL_MEM_DEVICE|CL_EVENT_NOWAIT);
  clmsync(CLCONTEXT,0,ivel,CL_MEM_DEVICE|CL_EVENT_NOWAIT);
  clmsync(CLCONTEXT,0,jpos,CL_MEM_DEVICE|CL_EVENT_NOWAIT);
  clmsync(CLCONTEXT,0,jvel,CL_MEM_DEVICE|CL_EVENT_NOWAIT);
  clfork(CLCONTEXT,0,timestep_krn,&ndr,CL_EVENT_NOWAIT);
  clmsync(CLCONTEXT,0,timestep,CL_MEM_HOST|CL_EVENT_NOWAIT);
  clwait(CLCONTEXT,0,CL_KERNEL_EVENT|CL_MEM_EVENT);

  for(i=0;i<s1.n;i++) 
  {
    ipart=GETPART(s1,i);
    ipart->timestep=timestep[i]; //if(timestep[i]<ipart->timestep)
  }

  clfree( ipos);
  clfree( ivel);
  clfree( jpos);
  clfree( jvel);
  clfree( timestep);  
}

void timestep_cl(struct sys s1, struct sys s2,int dir)
{
  if(accel_zero_mass)
  {
    struct sys s1m=zerosys, s1ml=zerosys, s2m=zerosys;
    
    s1m.n=s1.n-s1.nzero;
    if(s1m.n>0) s1m.part=s1.part;
    
    s1ml.n=s1.nzero;
    s1ml.nzero=s1.nzero;
    if(s1ml.n>0) s1ml.part=s1.zeropart;
    if(s1ml.n>0) s1ml.zeropart=s1.zeropart;
    
    s2m.n=s2.n-s2.nzero;
    if(s2m.n>0) s2m.part=s2.part;
    _timestep_cl(s1m,  s2, dir);
    _timestep_cl(s1ml, s2m, dir);
  } else
  {
    _timestep_cl(s1,s2,dir);
  }
}


void potential_cl(struct sys s1, struct sys s2)
{
  int i, n2;
  int groupsize,nthread,blocksize;
  CLFLOAT cleps2=(CLFLOAT) eps2;
  struct particle *ipart;
  
  n2=s2.n-s2.nzero;

  if(s1.n==0 || n2==0) return;

  groupsize=NTHREAD;
  if(s1.n < NTHREAD) groupsize=s1.n;
  nthread=((s1.n-1)/groupsize+1)*groupsize;
  blocksize=BLOCKSIZE;
  if( n2<blocksize) blocksize=n2;
  clndrange_t ndr = clndrange_init1d(0,nthread,groupsize);

  CLFLOAT4* ipos = (CLFLOAT4*) clmalloc(CLCONTEXT,nthread*sizeof(CLFLOAT4),0);
  CLFLOAT* pot = (CLFLOAT*) clmalloc(CLCONTEXT,nthread*sizeof(CLFLOAT),0);
  CLFLOAT4* jpos = (CLFLOAT4*) clmalloc(CLCONTEXT,n2*sizeof(CLFLOAT4),0);

  for(i=0;i<s1.n;i++)
  {
    ipart=GETPART(s1,i);
    ipos[i].s0=ipart->pos[0];
    ipos[i].s1=ipart->pos[1];
    ipos[i].s2=ipart->pos[2];
    ipos[i].s3=ipart->mass;
  }  
  for(i=s1.n;i<nthread;i++) ipos[i]=(CLFLOAT4) {{0.0,0.0,0.0,0.0}};  
  for(i=0;i<n2;i++)
  {
    ipart=GETPART(s2,i);
    jpos[i].s0=ipart->pos[0];
    jpos[i].s1=ipart->pos[1];
    jpos[i].s2=ipart->pos[2];
    jpos[i].s3=ipart->mass;
  }
  for(i=0;i<nthread;i++) pot[i]=0.0;  
                        
  clarg_set(CLCONTEXT,potential_krn,0, n2);
  clarg_set(CLCONTEXT,potential_krn,1, blocksize);
  clarg_set(CLCONTEXT,potential_krn,2, cleps2);
  clarg_set_global(CLCONTEXT,potential_krn,3, ipos);
  clarg_set_global(CLCONTEXT,potential_krn,4, jpos);
  clarg_set_global(CLCONTEXT,potential_krn,5, pot);
  clarg_set_local(CLCONTEXT,potential_krn,6,blocksize*sizeof(CLFLOAT4));

  clmsync(CLCONTEXT,0,ipos,CL_MEM_DEVICE|CL_EVENT_NOWAIT);
  clmsync(CLCONTEXT,0,jpos,CL_MEM_DEVICE|CL_EVENT_NOWAIT);
  clfork(CLCONTEXT,0,potential_krn,&ndr,CL_EVENT_NOWAIT);
  clmsync(CLCONTEXT,0,pot,CL_MEM_HOST|CL_EVENT_NOWAIT);
  clwait(CLCONTEXT,0,CL_KERNEL_EVENT|CL_MEM_EVENT);

  for(i=0;i<s1.n;i++) GETPART(s1,i)->pot+=pot[i];

  clfree( ipos);
  clfree( jpos);
  clfree( pot);  
}

#endif
