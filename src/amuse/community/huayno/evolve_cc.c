/*
 * The Connected Components Hamiltonian split uses a connected component search
 * on the time step graph of the system to find isolated subsystems with fast
 * interactions. These subsystems are then evolved at greater accuracy compared
 * to the rest system.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "evolve.h"

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

DOUBLE timestep_ij(struct sys r, UINT i, struct sys s, UINT j) {
  /*
   * timestep_ij: calculate timestep between particle i of sys r and particle j of sys s
   */
  FLOAT timestep;
	FLOAT dx[3],dr3,dr2,dr,dv[3],dv2,mu,vdotdr2,tau,dtau;
  timestep=HUGE_VAL;
  dx[0]=r.part[i].pos[0]-s.part[j].pos[0];
  dx[1]=r.part[i].pos[1]-s.part[j].pos[1];
  dx[2]=r.part[i].pos[2]-s.part[j].pos[2];
  dr2=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+eps2;
  if(dr2>0) {
    dr=sqrt(dr2);
    dr3=dr*dr2;
      dv[0]=r.part[i].vel[0]-s.part[j].vel[0];
      dv[1]=r.part[i].vel[1]-s.part[j].vel[1];
      dv[2]=r.part[i].vel[2]-s.part[j].vel[2];
      vdotdr2=(dv[0]*dx[0]+dv[1]*dx[1]+dv[2]*dx[2])/dr2;
      dv2=dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2];
      mu=r.part[i].mass+s.part[j].mass;
#ifdef RATIMESTEP
      tau=RARVRATIO*dt_param/M_SQRT2*sqrt(dr3/mu);
      dtau=3/2.*tau*vdotdr2;
      if(dtau>1.) dtau=1.;
      tau/=(1-dtau/2);
      if(tau < timestep) timestep=tau;
#endif
#ifdef RVTIMESTEP
      if(dv2>0) {
        tau=dt_param*dr/sqrt(dv2);
        dtau=tau*vdotdr2*(1+mu/(dv2*dr));
        if(dtau>1.) dtau=1.;
        tau/=(1-dtau/2);
        if(tau < timestep) timestep=tau;
      }
#endif
  }
  if (timestep < 0) {
    ENDRUN("negative timestep!\n");
  }
  tcount[clevel]++;
  return timestep;
}

void split_cc(struct sys s, struct sys *c, struct sys *r, DOUBLE dt) {
  /*
   * split_cc: run a connected component search on sys s with threshold dt,
   * creates a singly-linked list of connected components c and a rest system r
   * c or r is set to zerosys if no connected components/rest is found
   */
  tstep[clevel]++; // not directly comparable to corresponding SF-split statistics
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
				DOUBLE timestep = timestep_ij(s, comp_next, s, i);
				if (((dt > 0) && (timestep <= dt)) ||
				   ((dt < 0) && (timestep >= dt))) {
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

void split_cc_verify(struct sys s, struct sys *c, struct sys *r) {
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

void split_cc_verify_ts(struct sys *c, struct sys *r, DOUBLE dt) {

	DOUBLE ts_ij;
	// verify C-C interactions
    for (struct sys *ci = c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
    	for (UINT i = 0; i < ci->n; i++) {
       		for (struct sys *cj = c; !IS_ZEROSYS(cj); cj = cj->next_cc) {
    			if (ci == cj) {
    				continue;
    			}
    	    	for (UINT j = 0; j < cj->n; j++) {
    	    		ts_ij = timestep_ij(*ci, i, *cj, j);
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
	    		ts_ij = timestep_ij(*ci, i, *r, j);
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
    		ts_ij = timestep_ij(*r, i, *r, j);
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

DOUBLE sys_forces_max_timestep(struct sys s) {
  DOUBLE ts = 0.0;
  DOUBLE ts_ij;
  for (UINT i = 0; i < s.n; i++) {
    for (UINT j = 0; j < s.n; j++) {
      if (i != j) {
        ts_ij = timestep_ij(s, i, s, j);
        if (ts_ij >= ts) { ts = ts_ij; };
      }
    }
  }
  return ts;
}

void evolve_cc2(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt) {

	struct sys c = zerosys, r = zerosys;
	clevel++;
	if (etime <= stime ||  dt==0 || clevel>=MAXLEVEL) ENDRUN("timestep too small\n");

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

	split_cc(s, &c, &r, dt);
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
	split_cc_verify(s_before_split, &c, &r);
	split_cc_verify_ts(&c, &r, dt);
	free(s_before_split.part);
	if (clevel == 0) {
		printf("ok ");
	}
#endif

	if (IS_ZEROSYSs(c)) {
		deepsteps++;
		simtime+=dt;
	}

	// evolve all fast components, 1st time
	for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
		evolve_cc2(*ci, stime, stime+dt/2, dt/2);
	}

	drift(r, stime+dt/2, dt/2); // drift r, 1st time

	// kick ci <-> cj
	for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
		for (struct sys *cj = &c; !IS_ZEROSYS(cj); cj = cj->next_cc) {
			if (ci != cj) {
				kick(*ci, *cj, dt);
				//kick(*cj, *ci, dt);
			}
		}
	}

	// kick c <-> rest
	for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
		kick(r, *ci, dt);
		kick(*ci, r, dt);
	}

	kick(r, r, dt); // kick rest

	drift(r, etime, dt/2); 	// drift r, 2nd time

	// evolve all fast components, 2nd time
	for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
		evolve_cc2(*ci, stime+dt/2, etime, dt/2);
	}

	clevel--;
	free_sys(c.next_cc);
}


#ifdef SKIP
void evolve_cc2_twobody(struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt) {
  // TODO use kepler solver only if dt is larger than the orbital period
  if (s.n == 2) {
    //LOG("evolve: called \n");
    evolve_twobody(s, stime, etime, dt);
  } else {
    //LOG("evolve: regular!\n");
    struct sys c = zerosys, r = zerosys;
    clevel++;
    if (etime <= stime ||  dt==0 || clevel>=MAXLEVEL) ENDRUN("timestep too small\n");

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

    split_cc(s, &c, &r, dt);
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
    split_cc_verify(s_before_split, &c, &r);
    split_cc_verify_ts(&c, &r, dt);
    free(s_before_split.part);
    if (clevel == 0) {
      printf("ok ");
    }
#endif

#ifdef CC2_SPLIT_SHORTCUTS
    if (s.n == c.n) {
      DOUBLE initial_timestep = sys_forces_max_timestep(s);
      DOUBLE dt_step = dt;

      while (dt_step > initial_timestep) dt_step = dt_step / 2;

      LOG("CC2_SPLIT_SHORTCUTS clevel=%d dt/dt_step=%Le\n", clevel, dt / dt_step);
      for (DOUBLE dt_now = 0; dt_now < dt; dt_now += dt_step) {
        evolve_split_cc2(s, dt_now, dt_now + dt_step,(DOUBLE) dt_step);
      }

      clevel--;
      free_sys(c.next_cc);
      return;
    }
#endif

    if (IS_ZEROSYSs(c)) {
      deepsteps++;
      simtime+=dt;
    }

    // evolve all fast components, 1st time
    for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
      evolve_split_cc2_twobody(*ci, stime, stime+dt/2, dt/2);
    }

    drift(r, stime+dt/2, dt/2); // drift r, 1st time

    // kick ci <-> cj
    for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
      for (struct sys *cj = &c; !IS_ZEROSYS(cj); cj = cj->next_cc) {
        if (ci != cj) {
          kick(*ci, *cj, dt);
          //kick(*cj, *ci, dt);
        }
      }
    }

    // kick c <-> rest
    for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
      kick(r, *ci, dt);
      kick(*ci, r, dt);
    }

    kick(r, r, dt); // kick rest

    drift(r, etime, dt/2);  // drift r, 2nd time

    // evolve all fast components, 2nd time
    for (struct sys *ci = &c; !IS_ZEROSYS(ci); ci = ci->next_cc) {
      evolve_split_cc2_twobody(*ci, stime+dt/2, etime, dt/2);
    }

    clevel--;
    free_sys(c.next_cc);
  }
}
#endif
