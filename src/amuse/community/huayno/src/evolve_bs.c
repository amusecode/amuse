#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#include "evolve.h"
#include "evolve_cc.h"

#define FIXEDJ     10
#define JMAX       (FIXEDJ == 0 ? 16 : FIXEDJ)
#define BSTOL   1.e-6

int fixed_j=FIXEDJ;
DOUBLE bs_target_error=BSTOL;

struct bparticle
{
  DOUBLE pos[3];
  DOUBLE vel[3]; 
};

struct bsys
{
  UINT n; 
  struct bparticle *part;
};

int nsequence(int j)
{
  return 2*j; 
}

static int BulirschStoer(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int CCSUBSTEP);
static DOUBLE error_function(struct bsys s1, struct bsys s2);
static void aitkenneville(int j, int k, struct bsys* s, struct bsys s_jk, struct bsys s_j1k);
static void nkdk(int clevel,int n,struct sys s1,struct sys s2, DOUBLE stime, DOUBLE etime, DOUBLE dt);
static void ndkd(int clevel,int n,struct sys s1,struct sys s2, DOUBLE stime, DOUBLE etime, DOUBLE dt);
static void n_cc_kepler(int clevel,int n,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt);

void evolve_bs_adaptive(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep)
{
  FLOAT dtsys;
  int done=0;
  if(etime == stime ||  dt==0 || clevel>=MAXLEVEL)
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%d/%d\n", etime, stime, dt, clevel,MAXLEVEL);
  if(calc_timestep) timestep(clevel,s,s,SIGN(dt));
  dtsys=global_timestep(s);
  if(dtsys > fabs(dt))
  {
    done=BulirschStoer(clevel,s, stime, etime, dt, 0);
  }
  if(done==0)
  {      
    evolve_bs_adaptive(clevel+1,s,stime, stime+dt/2,dt/2,0);
    evolve_bs_adaptive(clevel+1,s,stime+dt/2, etime,dt/2,1);
  }
  else
  {
    diag->deepsteps++;
    diag->simtime+=dt;
  }
}

void evolve_bs(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  FLOAT dtsys;
  int done=0;
  if(etime == stime ||  dt==0 || clevel>=MAXLEVEL)
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%d/%d\n", etime, stime, dt, clevel,MAXLEVEL);
  done=BulirschStoer(clevel,s, stime, etime, dt, 1);
  if(done==0)
  {
    evolve_bs(clevel+1,s,stime, stime+dt/2,dt/2);
    evolve_bs(clevel+1,s,stime+dt/2, etime,dt/2);
  }
  else
  {
    diag->deepsteps++;
    diag->simtime+=dt;
  }
}


#define FREEBSYS_ARRAY(arr)  for(i=0; arr[i].part!=NULL && i<JMAX; i++) free(arr[i].part);
#define ZEROBSYS_ARRAY(arr)  for(i=0; i<JMAX; i++) arr[i].part=NULL;
#define ALLOCBSYS(bs,N) { \
                          if(bs.part==NULL) \
                          { \
                            bs.n=N; \
                            bs.part=(struct bparticle*) malloc(N*sizeof(struct bparticle)); \
                          } else \
                            if(bs.n != N) ENDRUN("bsys allocated but mismatch\n"); \
                          if(!bs.part) ENDRUN("bsys allocation error\n"); \
                        }

void sys_to_bsys(struct sys s, struct bsys bs)
{
  UINT i;
  if(s.n!=bs.n) ENDRUN("sys to bsys copy mismatch\n");
  if(!bs.part || !s.part) ENDRUN("sys to bsys unallocated error\n");
  for(i=0;i<s.n;i++)
  {
      bs.part[i].pos[0]=s.part[i].pos[0];
      bs.part[i].pos[1]=s.part[i].pos[1];
      bs.part[i].pos[2]=s.part[i].pos[2];
      bs.part[i].vel[0]=s.part[i].vel[0];
      bs.part[i].vel[1]=s.part[i].vel[1];
      bs.part[i].vel[2]=s.part[i].vel[2];   
  }  
}

void bsys_to_sys(struct bsys bs, struct sys s)
{
  UINT i;
  if(s.n!=bs.n) ENDRUN("bsys to sys copy mismatch %d,%d\n",s.n,bs.n);
  if(!bs.part || !s.part) ENDRUN("bsys to sys unallocated error\n");
  for(i=0;i<s.n;i++)
  {
      s.part[i].pos[0]=bs.part[i].pos[0];
      s.part[i].pos[1]=bs.part[i].pos[1];
      s.part[i].pos[2]=bs.part[i].pos[2];
      s.part[i].vel[0]=bs.part[i].vel[0];
      s.part[i].vel[1]=bs.part[i].vel[1];
      s.part[i].vel[2]=bs.part[i].vel[2];   
  }  
}

static int BulirschStoer(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int CCSUBSTEP)
{
  struct bsys bsys_array[JMAX];
  struct bsys bsys_array1[JMAX];
  struct bsys *jline;
  struct bsys *j1line;
  struct sys tmpsys;
  int j,k;
  UINT i;
  DOUBLE error;

  tmpsys.n=s.n;
  tmpsys.part=(struct particle*) malloc(s.n*sizeof(struct particle));
  if(!tmpsys.part) ENDRUN("failed allocation of tmpsys\n");

  ZEROBSYS_ARRAY(bsys_array)
  ZEROBSYS_ARRAY(bsys_array1)

  j=0;
  jline=bsys_array;
  j1line=bsys_array1;
  error=2*dt*bs_target_error;
  if(dt*bs_target_error <=0) ENDRUN("timestep or error too small\n");
  while(error > dt*bs_target_error)// || fixed_j!=0)
  {
    j=j+1;
    if(j>JMAX)
    {
      FREEBSYS_ARRAY(bsys_array);
      FREEBSYS_ARRAY(bsys_array1);
      free(tmpsys.part);
      return 0;
    }

    SWAP(jline,j1line, struct bsys*)
    for(i=0;i<s.n;i++) tmpsys.part[i]=s.part[i];
    if(CCSUBSTEP)
    {
      n_cc_kepler(clevel,nsequence(j),tmpsys,stime,etime,dt);
    } else
    {
      nkdk(clevel,nsequence(j),tmpsys,zerosys,stime,etime,dt);
    }
    ALLOCBSYS(jline[0], tmpsys.n)
    sys_to_bsys(tmpsys,jline[0]);
    for(k=1;k<j;k++) aitkenneville(j,k,jline+k,jline[k-1],j1line[k-1]) ;    
    if(j==1) continue;  
    error=error_function(jline[j-1],j1line[j-2]);
//    printf("err: %d %g\n", j, (double) error);
    if(j==fixed_j) break;
  }
  diag->bsstep[clevel]+=1;
  diag->jcount[clevel]+=j;
  bsys_to_sys(jline[j-1],tmpsys);
  for(i=0;i<s.n;i++) s.part[i]=tmpsys.part[i];  
#ifdef COMPENSATED_SUMMP
  for(i=0;i<s.n;i++) s.part[i].pos_e[0]=0.;
  for(i=0;i<s.n;i++) s.part[i].pos_e[1]=0.;
  for(i=0;i<s.n;i++) s.part[i].pos_e[2]=0.;
#endif
#ifdef COMPENSATED_SUMMV
  for(i=0;i<s.n;i++) s.part[i].vel_e[0]=0.;
  for(i=0;i<s.n;i++) s.part[i].vel_e[1]=0.;
  for(i=0;i<s.n;i++) s.part[i].vel_e[2]=0.;
#endif
  FREEBSYS_ARRAY(bsys_array);
  FREEBSYS_ARRAY(bsys_array1);
  free(tmpsys.part);
  return 1;
}

static DOUBLE error_function(struct bsys s1, struct bsys s2)
{
  DOUBLE maxdiv=0.;
  UINT i;
  if(s1.n!=s2.n) ENDRUN("error_function length mismatch %d,%d\n",s1.n,s2.n);
  for(i=0;i<s1.n;i++)
  {
    maxdiv=fmax(maxdiv,fabs(s1.part[i].pos[0]-s2.part[i].pos[0])); 
    maxdiv=fmax(maxdiv,fabs(s1.part[i].pos[1]-s2.part[i].pos[1])); 
    maxdiv=fmax(maxdiv,fabs(s1.part[i].pos[2]-s2.part[i].pos[2])); 
    maxdiv=fmax(maxdiv,fabs(s1.part[i].vel[0]-s2.part[i].vel[0])); 
    maxdiv=fmax(maxdiv,fabs(s1.part[i].vel[1]-s2.part[i].vel[1])); 
    maxdiv=fmax(maxdiv,fabs(s1.part[i].vel[2]-s2.part[i].vel[2])); 
  }
  return maxdiv;
}

void aitkenneville(int j, int k, struct bsys *s, struct bsys s_jk, struct bsys s_j1k)
{
  UINT i;
  DOUBLE fac;
  if(s_jk.n!=s_j1k.n) ENDRUN("aitken length mismatch\n");
  ALLOCBSYS( (*s), s_jk.n)
  fac=1./((nsequence(j)*nsequence(j))/((DOUBLE) nsequence(j-k)*nsequence(j-k))-1.);
  for(i=0;i<s->n;i++)
  {
    s->part[i].pos[0]=s_jk.part[i].pos[0]+fac*(s_jk.part[i].pos[0]-s_j1k.part[i].pos[0]);
    s->part[i].pos[1]=s_jk.part[i].pos[1]+fac*(s_jk.part[i].pos[1]-s_j1k.part[i].pos[1]);
    s->part[i].pos[2]=s_jk.part[i].pos[2]+fac*(s_jk.part[i].pos[2]-s_j1k.part[i].pos[2]);
    s->part[i].vel[0]=s_jk.part[i].vel[0]+fac*(s_jk.part[i].vel[0]-s_j1k.part[i].vel[0]);
    s->part[i].vel[1]=s_jk.part[i].vel[1]+fac*(s_jk.part[i].vel[1]-s_j1k.part[i].vel[1]);
    s->part[i].vel[2]=s_jk.part[i].vel[2]+fac*(s_jk.part[i].vel[2]-s_j1k.part[i].vel[2]);   
  }
}

static void nkdk(int clevel,int n,struct sys s1,struct sys s2, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  int i;
  if(s2.n>0) kick(clevel,s2, s1, dt/2/n);
  kick(clevel,s1,join(s1,s2),dt/2/n);
  for(i=0;i<n-1;i++)
  {  
    stime+=dt/n;
    drift(clevel,s1,stime, dt/n);
    kick(clevel,s1,join(s1,s2),dt/n);
    if(s2.n>0) kick(clevel,s2, s1, dt/n);
  }
  drift(clevel,s1,etime, dt/n);
  kick(clevel,s1,join(s1,s2),dt/2/n);
  if(s2.n>0) kick(clevel,s2, s1, dt/2/n);
}

static void ndkd(int clevel,int n,struct sys s1,struct sys s2, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  int i;
  stime+=dt/2/n;
  drift(clevel,s1,stime, dt/2/n);
  for(i=0;i<n-1;i++)
  {  
    kick(clevel,s1,join(s1,s2),dt/n);
    if(s2.n>0) kick(clevel,s2, s1, dt/n);
    stime+=dt/n;
    drift(clevel,s1,stime, dt/n);
  }
  kick(clevel,s1,join(s1,s2),dt/n);
  if(s2.n>0) kick(clevel,s2, s1, dt/n);
  drift(clevel,s1,etime, dt/2/n);
}

static void n_cc_kepler(int clevel,int n,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt)
{
  UINT id;
  struct sys tmpsys;
  FLOAT odt_param=dt_param;
  dt_param=dt_param/n;
  tmpsys.n=s.n;
  tmpsys.part=(struct particle*) malloc(s.n*sizeof(struct particle));
  tmpsys.last=tmpsys.part+s.n-1;
  for(UINT i=0;i<s.n;i++)
  {
    tmpsys.part[i]=s.part[i];
    tmpsys.part[i].id=i;
  }
  evolve_cc2(clevel,tmpsys, stime, etime, dt,CCC_KEPLER,1);
  for(UINT i=0;i<s.n;i++)
  {
    id=s.part[tmpsys.part[i].id].id;
    s.part[tmpsys.part[i].id]=tmpsys.part[i];
    s.part[tmpsys.part[i].id].id=id;
  }
  dt_param=odt_param;
  free(tmpsys.part);
}
