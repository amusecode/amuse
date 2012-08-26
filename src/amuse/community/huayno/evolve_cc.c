/*
 * The Connected Components Hamiltonian split uses a connected component search
 * on the time step graph of the system to find isolated subsystems with fast
 * interactions. These subsystems are then evolved at greater accuracy compared
 * to the rest system.
 */

#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "evolve.h"
#include "evolve_kepler.h"
#include "evolve_bs.h"

//#define CC_DEBUG // perform (time-consuming, but thorough) CC sanity checks

#define IS_ZEROSYS(SYS) (((SYS)->n == 0) && ((SYS)->part == NULL) && ((SYS)->last == NULL) && ((SYS)->next_cc == NULL))
#define IS_ZEROSYSs(SYS) (((SYS).n == 0) && ((SYS).part == NULL) && ((SYS).last == NULL) && ((SYS).next_cc == NULL))
#define LOG_CC_SPLIT(C, R) { \
	LOG("clevel = %d s.n = %d c.n = {", clevel, s.n); \
	for (struct sys *_ci = (C); !IS_ZEROSYS(_ci); _ci = _ci->next_cc) printf(" %d ", _ci->n ); \
	printf("} r.n = %d\n", (R)->n); \
};

#define LOGSYS_ID(SYS) for (UINT i = 0; i < (SYS).n; i++) { printf("%u ", (SYS).part[i].id); } printf("\n");
#define LOGSYSp_ID(SYS) for (UINT i = 0; i < (SYS)->n; i++) { printf("%u ", (SYS)->part[i].id); } printf("\n");
#define LOGSYSC_ID(SYS) for (struct sys *_ci = &(SYS); !IS_ZEROSYS(_ci); _ci = _ci->next_cc) {printf("{"); for (UINT i = 0; i < _ci->n; i++) {printf("%u ", _ci->part[i].id); } printf("}\t");} printf("\n");

void split_cc(int clevel,struct sys s, struct sys *c, struct sys *r, DOUBLE dt) {
  /*
   * split_cc: run a connected component search on sys s with threshold dt,
   * creates a singly-linked list of connected components c and a rest system r
   * c or r is set to zerosys if no connected components/rest is found
   */
  int dir=SIGN(dt); 
  dt=fabs(dt); 
  diag->tstep[clevel]++; // not directly comparable to corresponding SF-split statistics
	struct sys *c_next;
	c_next = c;
	*c_next = zerosys;
	UINT processed = 0; // increase if something is added from the stack to the cc
	UINT comp_next = 0; // increase to move a particle from stack to cc; points to the first element of the stack
	UINT comp_size = 0; // amount of particles added to the current cc
	UINT stack_next = 1; // swap this with s[i] to increase the stack
	UINT stack_size = 1; // first element of the stack is s[comp_next]
	                     // last element of the stack is s[comp_next + stack_size - 1]
	UINT rest_next = s.n - 1; // swap this to add to the rest-system
	// find connected components
	while (processed < s.n) {
		//LOG("split_cc: searching for connected components: %d / %d\n", processed, s.n);
		// search for the next connected component
		while (stack_size > 0) {
			// iterate over all unvisited elements
			for (UINT i = stack_next; i <= rest_next; i++) {
				// if element is connected to the first element of the stack
				DOUBLE timestep = (DOUBLE) timestep_ij(s.part+comp_next, s.part+i,dir);
				diag->tcount[clevel]++;
                                if ( timestep <= dt) {
					// add i to the end of the stack by swapping stack_next and i
					SWAP( s.part[ stack_next ], s.part[i], struct particle );
					stack_next++;
					stack_size++;
				}
			}
			// pop the stack; add to the connected component
			comp_size++;
			comp_next++;
			stack_size--;
		}
		processed += comp_size;
		// new component is non-trivial: create a new sys
		if (comp_size > 1) {
	    //LOG("split_cc: found component with size: %d\n", comp_size);
      // create new component c from u[0] to u[cc_visited - 1]
      // remove components from u (u.n, u.part)
		  c_next->n = comp_size;
			c_next->part = &( s.part[ comp_next - comp_size ]);
			c_next->last = &( s.part[ comp_next - 1 ]);
			c_next->next_cc = (struct sys*) malloc( sizeof(struct sys) );
			c_next = c_next->next_cc;
			*c_next = zerosys;
			comp_next = stack_next;
			comp_size = 0;
			stack_next = stack_next + 1;
			stack_size = 1;
		// new component is trivial: add to rest
		} else {
			//LOG("found trivial component; adding to rest\n");
			SWAP(s.part[ comp_next - 1 ], s.part[ rest_next ], struct particle );
			rest_next--;
			comp_next = comp_next - 1;
			comp_size = 0;
			stack_next = comp_next + 1;
			stack_size = 1;
		}
	}
	if (processed != s.n) {
	  ENDRUN("split_cc particle count mismatch: processed=%u s.n=%u r->n=%u\n", processed, s.n, r->n);
	}
	// create the rest system
	r->n = (s.n - 1) - rest_next;
	if (r->n > 0) {
		r->part = &( s.part[rest_next + 1] );
		r->last = s.last;
	} else {
		r->part = NULL;
		r->last = NULL;
	}
	//LOG("split_cc: rest system size: %d\n", r->n);
}

void split_cc_verify(int clevel,struct sys s, struct sys *c, struct sys *r) {
  /*
   * split_cc_verify: explicit verification if connected components c and rest system r form a correct
   * connected components decomposition of the system.
   */
  //LOG("split_cc_verify ping s.n=%d r->n=%d\n", s.n, r->n);
	//LOG_CC_SPLIT(c, r);
	UINT pcount_check = 0;
	for (UINT i = 0; i < s.n; i++) {
		pcount_check = 0;
		UINT particle_found = 0;
		struct particle *p = &( s.part[i] );
	    for (struct sys *cj = c; !IS_ZEROSYS(cj); cj = cj->next_cc) {
			pcount_check += cj->n;
			//LOG("%d\n", pcount_check);
			// search for p in connected components
			for (UINT k = 0; k < cj->n; k++) {
				struct particle * pk = &( cj->part[k] );
				// is pk equal to p
				if (p->id == pk->id) {
					particle_found += 1;
					//LOG("split_cc_verify: found in a cc\n");
				}
			}
			if (& ( cj->part[cj->n - 1] ) != cj->last) {
				LOG("split_cc_verify: last pointer for c is not set correctly!");
				LOG_CC_SPLIT(c, r);
				ENDRUN("data structure corrupted\n");
			}
	    }

		// search for p in rest
		for (UINT k = 0; k < r->n; k++) {
			struct particle * pk = &( r->part[k] );

			// is pk equal to p
			if (p->id == pk->id) {
				particle_found += 1;
			}
		}

		if (particle_found != 1) {
			LOG("split_cc_verify: particle %d particle_found=%d ", i, particle_found);
			LOG_CC_SPLIT(c, r);
			ENDRUN("data structure corrupted\n");
		}

	}

	//if (& ( r->part[r->n - 1] ) != r->last) {
	//	LOG("split_cc_verify: last pointer for r is not set correctly! %d %d",&( r->part[r->n - 1] ), r->last);
	//	LOG_CC_SPLIT(c, r);
	//	ENDRUN("data structure corrupted\n");
	//}

	if (pcount_check + r->n != s.n) {
		LOG("split_cc_verify: particle count mismatch (%d %d)\n", pcount_check + r->n, s.n);
		LOG_CC_SPLIT(c, r);
		ENDRUN("data structure corrupted\n");
		//ENDRUN("split_cc_verify: particle count mismatch\n");
	} else {
		//LOG("split_cc_verify pong\n");
	}

	//ENDRUN("Fin.\n");

}

void split_cc_verify_ts(int clevel,struct sys *c, struct sys *r, DOUBLE dt) {

	DOUBLE ts_ij;
        int dir=SIGN(dt);
        dt=fabs(dt);
	// verify C-C interactions
    for (struct sys *ci = c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
    	for (UINT i = 0; i < ci->n; i++) {
       		for (struct sys *cj = c; !IS_ZEROSYS(cj); cj = cj->next_cc) {
    			if (ci == cj) {
    				continue;
    			}
    	    	for (UINT j = 0; j < cj->n; j++) {
    	    		ts_ij = (DOUBLE) timestep_ij((*ci).part+i, (*cj).part+j, dir);
    	    		//LOG("comparing %d %d\n", ci->part[i].id, cj->part[j].id);
    	    		//LOG("%f %f \n", ts_ij, dt);
    	    		if (dt > ts_ij) {
    	    			ENDRUN("split_cc_verify_ts C-C timestep underflow\n");
    	    		}
        		}
    		}
    	}
    }

	// verify C-R interactions
    for (struct sys *ci = c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
    	for (UINT i = 0; i < ci->n; i++) {

    		for (UINT j = 0; j < r->n; j++) {
	    		ts_ij = (DOUBLE) timestep_ij( (*ci).part+ i, (*r).part+ j,dir);
    	   		if (ts_ij < dt) {
    	   			ENDRUN("split_cc_verify_ts C-R timestep underflow\n");
    	   		}
        	}
    	}
    }

    // verify R interactions
	for (UINT i = 0; i < r->n; i++) {
    	for (UINT j = 0; j < r->n; j++) {
    		if (i == j) continue;
    		ts_ij = (DOUBLE) timestep_ij( (*r).part+ i, (*r).part+j,dir);
    		if (ts_ij < dt) {
    			ENDRUN("split_cc_verify_ts R-R timestep underflow\n");
    		}
		}
	}

}

// TODO rename to cc_free_sys?
void free_sys(struct sys * s) {
	if (s==NULL) return;
	if (s->next_cc != NULL) {
		free_sys(s->next_cc);
	}
	free(s);
}

DOUBLE sys_forces_max_timestep(struct sys s,int dir) {
  DOUBLE ts = 0.0;
  DOUBLE ts_ij;
  for (UINT i = 0; i < s.n-1; i++) {
    for (UINT j = i+1; j < s.n; j++) {
        ts_ij = (DOUBLE) timestep_ij(s.part+ i, s.part+j,dir); // check symm.
        if (ts_ij >= ts) { ts = ts_ij; };
    }
  }
  return ts;
}

#ifdef COMPENSATED_SUMMP
#define COMPSUMP(sum,err,delta) \
  { \
    DOUBLE a; \
    a=sum; \
    err=err+delta; \
    sum=a+err; \
    err=err+(a-sum); \
  }
#else
#define COMPSUMP(sum,err,delta)  {sum+=delta;}
#endif

#ifdef COMPENSATED_SUMMV
#define COMPSUMV(sum,err,delta) \
  { \
    DOUBLE a; \
    a=sum; \
    err=err+delta; \
    sum=a+err; \
    err=err+(a-sum); \
  }
#else
#define COMPSUMV(sum,err,delta)  {sum+=delta;}
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

#define TASKCONDITION   (nc > 1)
void evolve_cc2(int clevel,struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt, int inttype, int recenter) {
  DOUBLE cmpos[3],cmvel[3];
  int recentersub=0;
	struct sys c = zerosys, r = zerosys;
  if(etime == stime ||  dt==0 || clevel>=MAXLEVEL)
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel);

  if (s.n == 2 && (inttype==CCC_KEPLER || inttype==CC_KEPLER)) {
    evolve_kepler(clevel,s, stime, etime, dt);
    return;
  }

  if (s.n <= 10 && (inttype==CCC_BS ||inttype==CC_BS)) {
    evolve_bs_adaptive(clevel,s, stime, etime, dt,1);
    return;
  }

  if(recenter && (inttype==CCC || inttype==CCC_KEPLER || inttype==CCC_BS))
  {
     system_center_of_mass(s,cmpos,cmvel);
     move_system(s,cmpos,cmvel,-1);
  }

// not actually helpful I think; needs testing
#ifdef CC2_SPLIT_SHORTCUTS
    int dir=SIGN(dt);
    DOUBLE initial_timestep = sys_forces_max_timestep(s, dir);
    if(fabs(dt) > initial_timestep)
    {
      DOUBLE dt_step = dt;
      while (fabs(dt_step) > initial_timestep)
      { 
        dt_step = dt_step / 2;
        clevel++;
      }
      LOG("CC2_SPLIT_SHORTCUTS clevel=%d dt/dt_step=%Le\n", clevel, dt / dt_step);
      for (DOUBLE dt_now = 0; dir*dt_now < dir*(dt-dt_step/2); dt_now += dt_step)
        evolve_cc2(clevel,s, dt_now, dt_now + dt_step, dt_step,inttype,0);
      return;
    }
#endif

#ifdef CC2_SPLIT_CONSISTENCY_CHECKS
	if (clevel == 0) {
		printf("consistency_checks: ", s.n, clevel);
	}
#endif

#ifdef CC2_SPLIT_CONSISTENCY_CHECKS
	// debug: make a copy of s to verify that the split has been done properly
	struct sys s_before_split;
	s_before_split.n = s.n;
	s_before_split.part = (struct particle*) malloc(s.n*sizeof(struct particle));
	s_before_split.last = &( s_before_split.part[s.n - 1] );
	s_before_split.next_cc = NULL;
	memcpy(s_before_split.part, s.part, s.n*sizeof(struct particle));
#endif

	split_cc(clevel,s, &c, &r, dt);
	//if (s.n != c.n) LOG_CC_SPLIT(&c, &r); // print out non-trivial splits

#ifdef CC2_SPLIT_CONSISTENCY_CHECKS
/*
  	if (s.n != r.n) {
		LOG("s: ");
		LOGSYS_ID(s_before_split);
		LOG("c: ");
		LOGSYSC_ID(c);
		LOG("r: ");
		LOGSYS_ID(r);
	}
*/
	// verify the split
	split_cc_verify(clevel,s_before_split, &c, &r);
	split_cc_verify_ts(clevel,&c, &r, dt);
	free(s_before_split.part);
	if (clevel == 0) {
		printf("ok ");
	}
#endif

	if (IS_ZEROSYSs(c)) {
		diag->deepsteps++;
		diag->simtime+=dt;
	}

	// evolve all fast components, 1st time

  int nc=0; for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) nc++;

  if(nc>1 || r.n>0) recentersub=1;

	for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
#ifdef _OPENMP
  if( TASKCONDITION )
  {
      diag->ntasks[clevel]++;
      diag->taskcount[clevel]+=ci->n;
#pragma omp task firstprivate(clevel,ci,stime,dt,recentersub) untied
    {
      struct sys lsys;
      lsys.n=ci->n;
      struct particle* lpart=(struct particle*) malloc(lsys.n*sizeof(struct particle));
      lsys.part=lpart;lsys.last=lpart+lsys.n-1;
      for(UINT i=0;i<lsys.n;i++) lsys.part[i]=ci->part[i];
      evolve_cc2(clevel+1,lsys, stime, stime+dt/2, dt/2,inttype,recentersub);
      for(UINT i=0;i<lsys.n;i++) ci->part[i]=lpart[i];
      free(lpart);
    }
  } else
#endif
  {
      evolve_cc2(clevel+1,*ci, stime, stime+dt/2, dt/2,inttype,recentersub);
  }
  }
#pragma omp taskwait  
  
	if(r.n>0) drift(clevel,r, stime+dt/2, dt/2); // drift r, 1st time

	// kick ci <-> cj
	for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
		for (struct sys *cj = &c; !IS_ZEROSYS(cj); cj = cj->next_cc) {
			if (ci != cj) {
				kick(clevel,*ci, *cj, dt);
				//kick(*cj, *ci, dt);
			}
		}
	}

	// kick c <-> rest
	if(r.n>0) for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
		kick(clevel,r, *ci, dt);
		kick(clevel,*ci, r, dt);
	}

	if(r.n>0) kick(clevel,r, r, dt); // kick rest

	if(r.n>0) drift(clevel,r, etime, dt/2); 	// drift r, 2nd time

	// evolve all fast components, 2nd time
	for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
#ifdef _OPENMP
  if( TASKCONDITION)
  {
      diag->ntasks[clevel]++;
      diag->taskcount[clevel]+=ci->n;
#pragma omp task firstprivate(clevel,ci,stime,etime,dt,recentersub) untied
    {
      struct sys lsys;
      lsys.n=ci->n;
      struct particle* lpart=(struct particle*) malloc(lsys.n*sizeof(struct particle));
      lsys.part=lpart;lsys.last=lpart+lsys.n-1;
      for(UINT i=0;i<lsys.n;i++) lsys.part[i]=ci->part[i];
      evolve_cc2(clevel+1,lsys, stime+dt/2, etime, dt/2,inttype,recentersub);
      for(UINT i=0;i<lsys.n;i++) ci->part[i]=lpart[i];
      free(lpart);
    }
  } else
#endif
  {
      evolve_cc2(clevel+1,*ci, stime, stime+dt/2, dt/2,inttype,recentersub);
  }
  }
#pragma omp taskwait  
  
  if(recenter && (inttype==CCC || inttype==CCC_KEPLER || inttype==CCC_BS))
  {
    for(int i=0;i<3;i++) cmpos[i]+=cmvel[i]*dt;
    move_system(s,cmpos,cmvel,1);
  }

	free_sys(c.next_cc);
}
#undef TASKCONDITION

