/*
 * The Connected Components Hamiltonian split uses a connected component search
 * on the time step graph of the system to find isolated subsystems with fast
 * interactions. These subsystems are then evolved at greater accuracy compared
 * to the rest system.
 * Equation numbers in comments refer to: J\"anes, Pelupessy, Portegies Zwart, A&A 2014 (doi:10.1051/0004-6361/201423831)
 */

#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <string.h>

#include "evolve.h"
#include "evolve_kepler.h"
#include "evolve_bs.h"
#include "evolve_shared.h"

#define BS_SUBSYS_SIZE   10
#define SHARED10_SUBSYS_SIZE   10
#define SHARED10_MIN_DT_RATIO   (1./128) // if dtmin/dtmax > this do shared10, else try to find further CC 
                                        // (ie see if there are hard binaries) 

struct ccsys
{
  struct sys s; 
  struct ccsys *next_cc; 
};

#define LOG_CC_SPLIT(C, R) \
{ \
  LOG("clevel = %d s.n = %d c.n = {", clevel, s.n); \
  for (struct ccsys *_ci = (C); _ci!=NULL; _ci = _ci->next_cc) printf(" %d ", _ci->s.n ); \
  printf("} r.n = %d\n", (R).n); \
};

#define PRINTSYS(s) \
{ \
  LOG("sys %d %d : {",s.n,s.nzero); \
  for (UINT i=0;i<s.n-s.nzero;i++) printf(" %d", GETPART(s, i)->id); printf(" | "); \
  for (UINT i=s.n-s.nzero;i<s.n;i++) printf(" %d", GETPART(s, i)->id); printf(" }\n"); \
};

#define PRINTOFFSETS(s) \
{ \
  LOG("sysoffsets %d %d : {",s.n,s.nzero); \
  for (UINT i=0;i<s.n-s.nzero;i++) printf(" %d", GETPART(s, i)-s.part); printf(" | "); \
  for (UINT i=s.n-s.nzero;i<s.n;i++) printf(" %d", GETPART(s, i)-s.part); printf(" }\n"); \
};


#define LOGSYS_ID(SYS) for (UINT i = 0; i < (SYS).n; i++) { printf("%u ", GETPART(SYS, i)->id); } printf("\n");
#define LOGSYSp_ID(SYS) LOGSYS_ID(*SYS);
#define LOGSYSC_ID(SYS) for (struct ccsys *_ci = &(SYS); _ci!=NULL; _ci = _ci->next_cc) \
 {printf("{"); for (UINT i = 0; i < _ci->s.n; i++) {printf("%u ", GETPART(_ci->s,i)->id); } printf("}\t");} printf("\n");

void split_cc(int clevel,struct sys s, struct ccsys **c, struct sys *r, DOUBLE dt) {
  /*
   * split_cc: run a connected component search on sys s with threshold dt,
   * creates a singly-linked list of connected components c and a rest system r
   * c or r is set to zerosys if no connected components/rest is found
   */
  int dir=SIGN(dt);
  dt=fabs(dt);
  diag->tstep[clevel]++; // not directly comparable to corresponding SF-split statistics
  struct ccsys **c_next;
  if(s.n<=1) ENDRUN("This does not look right...");
  c_next = c; 
  if(*c_next!=NULL) ENDRUN("should start with zero pointer");
  UINT processed = 0; // increase if something is added from the stack to the cc
  struct particle **active, *comp_next, *compzero_next; // current active particle, and next in the mass and massless part 
  UINT comp_size, compzero_size;
  struct particle *stack_next=NULL, *stackzero_next=NULL; // two pointers keeping track of stack in massive and massless part 
  UINT stack_size, stackzero_size; //  counter for total size of stack and zero 
  struct particle *rest_next=NULL, *restzero_next=NULL; // swap this to add to the rest-system
  // find connected components
  
  if(s.n-s.nzero>0) stack_next=s.part;
  if(s.n-s.nzero>0) rest_next=LAST(s);
  if(s.nzero>0) stackzero_next=s.zeropart;
  if(s.nzero>0) restzero_next=LASTZERO(s);
  comp_next=stack_next;
  compzero_next=stackzero_next;

  while (processed < s.n)
  {
    if(stack_next!=comp_next) ENDRUN("consistency error in split_cc\n")
    if(stackzero_next!=compzero_next) ENDRUN("consistency error in split_cc\n")
    //~ if(stack_next==rest_next && stackzero_next==restzero_next) ENDRUN("impossible")

    // startup stack 

    comp_size=0;
    compzero_size=0;
    if(stack_next!=NULL && stack_next<rest_next+1)
    {
      //~ LOG("stack_next init\n");
      stack_next++;
      stack_size=1;
    } 
    if(comp_next==stack_next && stackzero_next!=NULL  && stackzero_next<restzero_next+1)
    {
      //~ LOG("stackzero_next init\n");
      stackzero_next++;
      stack_size=1;
    } 
    if(stack_next==comp_next && stackzero_next==compzero_next) ENDRUN("impossible")
    
    // search for the next connected component
    while (stack_size > 0)
    {
      //~ LOG("stack_size %d\n", stack_size);
      active=NULL;
      if(stack_next!=NULL &&
         stack_next-comp_next>0) {active=&comp_next;}
      else
        if(stackzero_next!=NULL && 
           stackzero_next-compzero_next>0) {active=&compzero_next;}
  
      if(active==NULL) ENDRUN("no active, while stack still >0\n");

      // iterate over all unvisited elements
      if(stack_next!=NULL)
      {
        //~ LOG("check massive %d\n", rest_next-stack_next+1);
        for (struct particle *i = stack_next; i <= rest_next; i++)
        {
          diag->tcount[clevel]++;
          // if element is connected to the first element of the stack
          if ( ((DOUBLE) timestep_ij(*active, i,dir)) <= dt)
          {
            // add i to the end of the stack by swapping stack_next and i
            //~ LOG("stack_next add %d\n", i->id);
            //~ LOG("stack offsets: %d, %d\n", stack_next-s.part, i-s.part);
            SWAP( *stack_next , *i, struct particle );
            stack_next++;
            stack_size++;
          }
        }
      }
      // iterate over all unvisited elements, skip when active is zero mass
      if(stackzero_next!=NULL && active!=&compzero_next)
      {
        //~ LOG("check zero %d\n", restzero_next-stackzero_next+1);
        for (struct particle *i = stackzero_next; i <= restzero_next; i++)
        {
          diag->tcount[clevel]++;
          // if element is connected to the first element of the stack
          if ( ((DOUBLE) timestep_ij(*active, i,dir)) <= dt)
          {
            // add i to the end of the stack by swapping stack_next and i
            //~ LOG("stackzero_next add %d\n", i->id);
            //~ LOG("stack offsets: %d, %d\n", stackzero_next-s.part, i-s.part);
            SWAP( *stackzero_next , *i, struct particle );
            stackzero_next++;
            stack_size++;
          }
        }
      }
      // pop the stack
      (*active)++;
      if(active==&compzero_next) compzero_size++;
      comp_size++;
      stack_size--;
      //~ LOG("popped %d, %d, %d\n", stack_size, comp_size, compzero_size);

    }
    processed += comp_size;
    //~ LOG("comp finish %d, %d\n", comp_size, compzero_size);
    // new component is non-trivial: create a new sys
    if (comp_size > 1)
    {
      //~ LOG("split_cc: found component with size: %d %d\n", comp_size, compzero_size);
      //~ LOG("%d %d \n", comp_next-stack_next, compzero_next-stackzero_next);
      *c_next=(struct ccsys*) malloc( sizeof(struct ccsys) );
      struct sys *new=&((*c_next)->s); 
      *new=zerosys;
      new->n = comp_size;
      new->nzero = compzero_size;
      if(comp_size-compzero_size>0)
      {
        new->part = comp_next - (comp_size-compzero_size);
      }
      if(compzero_size>0)
      {
        new->zeropart = compzero_next - compzero_size;
      }
      if(new->part==NULL) new->part=new->zeropart;
      //~ PRINTSYS((*c_next));
      //~ PRINTOFFSETS((*c_next));
      (*c_next)->next_cc = NULL;
      c_next = &((*c_next)->next_cc);
    }
    else  // new component is trivial: add to rest, reset pointers
    {
      //~ LOG("split_cc: add to rest\n");
      //~ LOG("%d %d \n", comp_next-stack_next, compzero_next-stackzero_next);
      
      if(active==&comp_next)
      {
        comp_next--;
        //~ LOG("r1 check offsets: %d, %d\n", comp_next-s.part, rest_next-s.part);
        
        SWAP( *comp_next, *rest_next, struct particle );
        rest_next--;
        stack_next--;
      } else
      {
        compzero_next--;
        //~ LOG("r2 check offsets: %d, %d\n", compzero_next-s.part, restzero_next-s.part);
        SWAP( *compzero_next, *restzero_next, struct particle );
        restzero_next--;
        stackzero_next--; 
      }
    }
  }
  if(stack_next!=NULL && stack_next!=rest_next+1) ENDRUN("unexpected")
  if(stackzero_next!=NULL && stackzero_next!=restzero_next+1) ENDRUN("unexpected")

  // create the rest system 
  *r=zerosys;
  if(rest_next!=NULL) r->n = LAST(s) - rest_next;
  if(restzero_next!=NULL) r->nzero = LASTZERO(s) - restzero_next;
  r->n+=r->nzero;
  if (r->n-r->nzero > 0)
  {
    r->part = rest_next+1;
  }
  if (r->nzero > 0)
  {
    r->zeropart = restzero_next+1;
  }
  if(r->part==NULL) r->part=r->zeropart;

  if (processed != s.n)
  {
    ENDRUN("split_cc particle count mismatch: processed=%u s.n=%u r->n=%u\n", processed, s.n, r->n);
  }

  //~ LOG("exit with %d %d\n", s.n, s.nzero);

}



void split_cc_verify(int clevel,struct sys s, struct ccsys *c, struct sys r) {
  /*
   * split_cc_verify: explicit verification if connected components c and rest system r form a correct
   * connected components decomposition of the system.
   */
  //~ LOG("split_cc_verify ping s.n=%d r->n=%d\n", s.n, r->n);
  //LOG_CC_SPLIT(c, r);
  UINT pcount_check = 0;

  for (UINT i = 0; i < s.n; i++)
  {
    pcount_check = 0;
    UINT particle_found = 0;
    struct particle *p = GETPART(s, i);
    for (struct ccsys *cj = c; cj!=NULL; cj = cj->next_cc)
    {
      verify_split_zeromass(cj->s);
      pcount_check += cj->s.n;
      //~ //LOG("%d\n", pcount_check);
      // search for p in connected components
      for (UINT k = 0; k < cj->s.n; k++)
      {
        struct particle * pk = GETPART( cj->s,k);
        // is pk equal to p
        if (p->id == pk->id)
        {
          particle_found += 1;
          //~ LOG("split_cc_verify: found %d in a cc\n",i);
        }
      }
      if (cj->s.n-cj->s.nzero>0 &&  (  GETPART( cj->s, cj->s.n - cj->s.nzero - 1) != LAST(cj->s) ))
      {
        LOG("split_cc_verify: last pointer for c is not set correctly!\n");
        LOG_CC_SPLIT(c, r);
        ENDRUN("data structure corrupted\n");
      }
      if (cj->s.nzero>0 &&  (  GETPART( cj->s, cj->s.n-1) != LASTZERO(cj->s) ))
      {
        LOG("split_cc_verify: last pointer for c is not set correctly!\n");
        LOG_CC_SPLIT(c, r);
        ENDRUN("data structure corrupted\n");
      }
    }

    verify_split_zeromass(r);

    // search for p in rest
    for (UINT k = 0; k < r.n; k++)
    {
      struct particle * pk = GETPART( r, k);

      // is pk equal to p
      if (p->id == pk->id)
      {
        particle_found += 1;
        //~ LOG("found at r\n")
      }
    }

    if (particle_found != 1)
    {
      LOG("split_cc_verify: particle %d (%d) particle_found=%d\n", i, p->id, particle_found);
      LOG_CC_SPLIT(c, r);
      ENDRUN("data structure corrupted\n");
    }
  }

  if (pcount_check + r.n != s.n)
  {
    LOG("split_cc_verify: particle count mismatch (%d %d)\n", pcount_check + r.n, s.n);
    LOG_CC_SPLIT(c, r);
    ENDRUN("data structure corrupted\n");
    //ENDRUN("split_cc_verify: particle count mismatch\n");
  }
  else
  {
     //~ LOG("split_cc_verify pong\n");
  }
  //ENDRUN("Fin.\n");
}

void split_cc_verify_ts(int clevel,struct ccsys *c, struct sys r, DOUBLE dt)
{
  DOUBLE ts_ij;
  int dir=SIGN(dt);
  dt=fabs(dt);
  // verify C-C interactions
  for (struct ccsys *ci = c; ci!=NULL; ci = ci->next_cc)
  {
    for (UINT i = 0; i < ci->s.n; i++)
    {
      for (struct ccsys *cj = c; cj!=NULL; cj = cj->next_cc)
      {
        if (ci == cj)
        {
          continue;
        }
        for (UINT j = 0; j < cj->s.n; j++)
        {
          ts_ij = (DOUBLE) timestep_ij(GETPART( ci->s, i), GETPART( cj->s, j), dir);
          //LOG("comparing %d %d\n", GETPART( ci->s, i)-> id, GETPART( cj->s, j)->id);
          //LOG("%f %f \n", ts_ij, dt);
          if (dt > ts_ij)
          {
            ENDRUN("split_cc_verify_ts C-C timestep underflow\n");
          }
        }
      }
    }
  }

  // verify C-R interactions
  for (struct ccsys *ci = c; ci!=NULL; ci = ci->next_cc)
  {
    for (UINT i = 0; i < ci->s.n; i++)
    {
      for (UINT j = 0; j < r.n; j++)
      {
        ts_ij = (DOUBLE) timestep_ij( GETPART(ci->s, i), GETPART(r, j),dir);
        if (ts_ij < dt)
        {
          ENDRUN("split_cc_verify_ts C-R timestep underflow\n");
        }
      }
    }
  }

  // verify R interactions
  for (UINT i = 0; i < r.n; i++)
  {
    for (UINT j = 0; j < r.n; j++)
    {
      if (i == j) continue;
      ts_ij = (DOUBLE) timestep_ij( GETPART(r, i), GETPART(r,j),dir);
      if (ts_ij < dt)
      {
        ENDRUN("split_cc_verify_ts R-R timestep underflow\n");
      }
    }
  }
}

// TODO rename to cc_free_sys?
void free_sys(struct ccsys * s)
{
  if (s==NULL) return;
  if (s->next_cc != NULL)
  {
    free_sys(s->next_cc);
  }
  free(s);
}

#define TASKCONDITION    (nc > 1 && s.n>BS_SUBSYS_SIZE)
void evolve_cc2(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int inttype, int recenter)
{
  DOUBLE cmpos[3],cmvel[3];
  int recentersub=0;
  struct ccsys *c = NULL;
  struct sys r = zerosys;
  CHECK_TIMESTEP(etime,stime,dt,clevel);

  if ((s.n == 2 || s.n-s.nzero<=1 )&& 
     (inttype==CCC_KEPLER || inttype==CC_KEPLER || inttype==CCC_BS || 
      inttype==CC_BS || inttype==CCC_BSA || inttype==CC_BSA || inttype==CC_SHARED10 || inttype==CCC_SHARED10))
  //~ if (s.n == 2 && (inttype==CCC_KEPLER || inttype==CC_KEPLER))
  {
    evolve_kepler(clevel,s, stime, etime, dt);
    return;
  }

  if(recenter && (inttype==CCC || inttype==CCC_KEPLER || inttype==CCC_BS || inttype==CCC_BSA || inttype==CCC_SHARED10))
  {
    system_center_of_mass(s,cmpos,cmvel);
    move_system(s,cmpos,cmvel,-1);
    evolve_cc2(clevel, s, stime, etime, dt, inttype, 0);
    for(int i=0;i<3;i++) cmpos[i]+=cmvel[i]*dt;
    move_system(s,cmpos,cmvel,1);
    return;
  }

  if (s.n <= BS_SUBSYS_SIZE && (inttype==CCC_BS ||inttype==CC_BS))
  {
    evolve_bs(clevel,s, stime, etime, dt);
    return;
  }

  if (s.n <= BS_SUBSYS_SIZE && (inttype==CCC_BSA ||inttype==CC_BSA))
  {
    evolve_bs_adaptive(clevel,s, stime, etime, dt, -1.);
    return;
  }

  if (s.n <= SHARED10_SUBSYS_SIZE && (inttype==CCC_SHARED10 ||inttype==CC_SHARED10))
  {
    timestep(clevel,s,s,SIGN(dt));
    FLOAT dtmax=max_global_timestep(s);    
    FLOAT dtmin=global_timestep(s);    
    if(dtmin/dtmax>SHARED10_MIN_DT_RATIO) 
    {
      evolve_shared10(clevel,s, stime, etime, dt, -1.);
      return;
    }
  }



#ifdef CONSISTENCY_CHECKS
  if (clevel == 0)
  {
    printf("consistency_checks: %d %d \n", s.n, clevel);
  }
#endif

#ifdef CONSISTENCY_CHECKS
  // debug: make a copy of s to verify that the split has been done properly
  struct sys s_before=zerosys;
  s_before.n = s.n;
  s_before.nzero = s.nzero;
  s_before.part = (struct particle*) malloc(s.n*sizeof(struct particle));
  if(s_before.nzero>0) s_before.zeropart = s_before.part+(s_before.n-s_before.nzero);
  for(UINT i=0; i<s.n;i++) *GETPART(s_before, i)=*GETPART(s,i);
#endif


  /*
   split_cc() decomposes particles in H (eq 25) into:
   1) K non-trivial connected components C_1..C_K
   2) Rest set R
  */
  split_cc(clevel,s, &c, &r, dt);
  //if (s.n != c.n) LOG_CC_SPLIT(&c, &r); // print out non-trivial splits

#ifdef CONSISTENCY_CHECKS
/*
    if (s.n != r.n) {
    LOG("s: ");
    LOGSYS_ID(s_before);
    LOG("c: ");
    LOGSYSC_ID(*c);
    LOG("r: ");
    LOGSYS_ID(r);
  }
*/
  // verify the split
  split_cc_verify(clevel,s_before, c, r);
  split_cc_verify_ts(clevel, c, r, dt);
  free(s_before.part);
  if (clevel == 0) {
    printf("ok \n");
  }
#endif

  if (c==NULL) {
    diag->deepsteps++;
    diag->simtime+=dt;
  }

  // Independently integrate every C_i at reduced pivot time step h/2 (1st time)
  int nc=0; for (struct ccsys *ci = c; ci!=NULL; ci = ci->next_cc) nc++;
  
  if(nc>1 || r.n>0) recentersub=1;

  for (struct ccsys *ci = c; ci!=NULL; ci = ci->next_cc)
  {
#ifdef _OPENMP
    if( TASKCONDITION )
    {
      diag->ntasks[clevel]++;
      diag->taskcount[clevel]+=ci->s.n;
#pragma omp task firstprivate(clevel,ci,stime,dt,recentersub) untied
      {
        struct sys lsys=zerosys;
        lsys.n=ci->s.n;
        lsys.nzero=ci->s.nzero;
        struct particle* lpart=(struct particle*) malloc(lsys.n*sizeof(struct particle));
        lsys.part=lpart;
        if(lsys.nzero>0) lsys.zeropart=lsys.part+(lsys.n-lsys.nzero);
        
        for(UINT i=0;i<lsys.n;i++) *GETPART(lsys,i)=*GETPART(ci->s,i);
      
        evolve_cc2(clevel+1,lsys, stime, stime+dt/2, dt/2,inttype,recentersub);

        for(UINT i=0;i<lsys.n;i++) *GETPART(ci->s,i)=*GETPART(lsys,i);
      
        free(lpart);
      }
    } else
#endif
  {
      evolve_cc2(clevel+1,ci->s, stime, stime+dt/2, dt/2,inttype,recentersub);
  }
  }
#pragma omp taskwait

  // Apply drifts and kicks at current pivot time step (eq 30)
  if(r.n>0) drift(clevel,r, stime+dt/2, dt/2); // drift r, 1st time

  // kick ci <-> cj (eq 23)
  for (struct ccsys *ci = c; ci!=NULL; ci = ci->next_cc)
  {
    for (struct ccsys *cj = c; cj!=NULL; cj = cj->next_cc)
    {
      if (ci != cj)
      {
        kick(clevel,ci->s, cj->s, dt);
        //kick(*cj, *ci, dt);
      }
    }
  }

  // kick c <-> rest (eq 24)
  if(r.n>0) for (struct ccsys *ci = c; ci!=NULL; ci = ci->next_cc)
  {
    kick(clevel,r, ci->s, dt);
    kick(clevel,ci->s, r, dt);
  }

  if(r.n>0) kick(clevel,r, r, dt); // kick rest (V_RR)

  if(r.n>0) drift(clevel,r, etime, dt/2); // drift r, 2nd time

  // Independently integrate every C_i at reduced pivot time step h/2 (2nd time, eq 27)
  for (struct ccsys *ci = c; ci!=NULL; ci = ci->next_cc)
  {
#ifdef _OPENMP
    if (TASKCONDITION)
    {
      diag->ntasks[clevel]++;
      diag->taskcount[clevel]+=ci->s.n;
#pragma omp task firstprivate(clevel,ci,stime,etime,dt,recentersub) untied
      {
        struct sys lsys=zerosys;
        lsys.n=ci->s.n;
        lsys.nzero=ci->s.nzero;
        struct particle* lpart=(struct particle*) malloc(lsys.n*sizeof(struct particle));
        lsys.part=lpart;
        if(lsys.nzero>0) lsys.zeropart=lsys.part+(lsys.n-lsys.nzero);
        
        for(UINT i=0;i<lsys.n;i++) *GETPART(lsys,i)=*GETPART(ci->s,i);
        
        evolve_cc2(clevel+1,lsys, stime+dt/2, etime, dt/2,inttype,recentersub);
        
        for(UINT i=0;i<lsys.n;i++) *GETPART(ci->s,i)=*GETPART(lsys,i);
        
        free(lpart);
      }
    } else
#endif
    {
      evolve_cc2(clevel+1,ci->s, stime+dt/2, etime, dt/2,inttype,recentersub);
    }
  }
#pragma omp taskwait

  free_sys(c);

}

// not actually helpful I think; needs testing
void evolve_cc2_shortcut(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int inttype, int recenter, FLOAT dtsys)
{
  CHECK_TIMESTEP(etime,stime,dt,clevel);
  if(dtsys<0) 
  {
    timestep(clevel,s,s,SIGN(dt));
    dtsys=max_global_timestep(s);
  }
  if(dtsys < fabs(dt))
  {
    evolve_cc2_shortcut(clevel+1,s,stime, stime+dt/2,dt/2, inttype, recenter, dtsys);
    evolve_cc2_shortcut(clevel+1,s,stime+dt/2, etime,dt/2, inttype, recenter, -1.);
  }
  else
  {
    evolve_cc2(clevel,s, stime, etime, dt, inttype, recenter);
  }
}


#undef TASKCONDITION


