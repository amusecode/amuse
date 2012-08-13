/*
 * The Optimal Kick Hamiltonian split evaluates every interaction in the system
 * at the optimal split time step, thereby using the least amount of force
 * evaluations to evolve the system.
 */

#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#include "evolve.h"
#include "evolve_ok.h"

struct forces zeroforces = {0, NULL, NULL};

#define IS_ZEROFORCES(F) (((F).n == 0) && ((F).forc == NULL) && ((F).last == NULL))

#define LOG_FORCES(F) { \
  for (UINT i = 0; i < (F).n; i++) { \
    printf("%u\t%u\t%f\n", (F).forc[i].parti->id, (F).forc[i].partj->id, (F).forc[i].timestep); \
  } \
};

static void ok_timestep_cpu(struct forces f, DOUBLE dt) {
  int dir=SIGN(dt);
  for (UINT i = 0; i < f.n; i++) 
  {
      //if (f.forc[i].timestep != HUGE_VAL) ENDRUN("timestep??");
      f.forc[i].timestep = timestep_ij(f.forc[i].parti, f.forc[i].partj,dir);
  }
  diag->tstep[clevel]++;
  diag->tcount[clevel] += f.n;
}


/*
 * split_ok_forces: split forces into smaller than dt, faster than dt
 */
static void ok_split(FLOAT dt, struct forces f, struct forces *slow, struct forces *fast) {
  //LOG("dt=%lf f.n=%u\n", dt, f.n);
  UINT i = 0;
  struct force *left, *right;
  left = f.forc;
  right = f.last;
  dt=fabs(dt);
  while (1) {
    if (i >= f.n) ENDRUN("forces split error 1\n");
    i++;
    while ((left->timestep < dt) && (left<right)) left++;
    while ((right->timestep >= dt) && (left<right)) right--;
    if (left < right) {
      SWAP(*left, *right, struct force);
    } else
      break;
  }
  if (left->timestep < dt) left++;
  slow->n = f.last - left + 1;
  fast->n = left - f.forc;
  if (fast->n == 1) {
    fast->n = 0;
    slow->n = f.n;
  }
  if (slow->n > 0) {
    slow->forc = f.forc + fast->n;
    slow->last = f.last;//slow->part+slow->n-1;
  }
  if (fast->n > 0) {
    fast->forc = f.forc;
    fast->last = f.forc + fast->n - 1;
  }
  if (fast->n + slow->n != f.n)
    ENDRUN("forces split error 2: fast->n=%u slow->n=%u f.n=%u\n", fast->n, slow->n, f.n);
  //for (i = 0; i < f.n; i++) f.forc[i].level = clevel;
}

struct forces ok_main_forces = {0, NULL, NULL};

void evolve_ok_init(struct sys s) {
  UINT n_forces = s.n * s.n - s.n;
  if (ok_main_forces.forc != NULL) ENDRUN("OK (re)allocation error");  
  ok_main_forces.forc = (struct force *) malloc(n_forces * sizeof(struct force));
  ok_main_forces.last = &(ok_main_forces.forc[n_forces - 1]);
  ok_main_forces.n = n_forces;

  // initialize pointers of the forces structure
  UINT k = 0;
  for (UINT i = 0; i < s.n; i++) {
    for (UINT j = 0; j < s.n; j++) {
      if (i != j) {
        ok_main_forces.forc[k].parti = &( s.part[i] );
        ok_main_forces.forc[k].partj = &( s.part[j] );
        k++;
      }
    }
  }
}

void evolve_ok_stop() {
  if (ok_main_forces.forc != NULL) {
    free(ok_main_forces.forc);
  }
}

static void ok_kick(struct forces f, DOUBLE dt) {
  FLOAT dx[3],dr3,dr2,dr,acci;
  FLOAT acc[3];

  for (UINT i = 0; i < f.n; i++) {
    acc[0] = 0.;
    acc[1] = 0.;
    acc[2] = 0.;

    dx[0] = f.forc[i].parti->pos[0] - f.forc[i].partj->pos[0];
    dx[1] = f.forc[i].parti->pos[1] - f.forc[i].partj->pos[1];
    dx[2] = f.forc[i].parti->pos[2] - f.forc[i].partj->pos[2];
    dr2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] + eps2;

    if (dr2 > 0) {
      dr = sqrt(dr2);
      dr3 = dr*dr2;
      acci = f.forc[i].partj->mass / dr3;

      f.forc[i].parti->vel[0] -= dt * dx[0] * acci;
      f.forc[i].parti->vel[1] -= dt * dx[1] * acci;
      f.forc[i].parti->vel[2] -= dt * dx[2] * acci;
    }
  }

  diag->kstep[clevel]++;
  diag->kcount[clevel] += f.n;
}

void evolve_ok2(struct sys s, struct forces f, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep) {
  if (IS_ZEROFORCES(f) && clevel == -1) { f = ok_main_forces; }
  clevel++;
  if ((etime == stime) || (dt == 0) || (clevel >= MAXLEVEL))
    ENDRUN("timestep too small: etime=%Le stime=%Le dt=%Le clevel=%u\n", etime, stime, dt, clevel);
  // all particles are drifted together
  if (f.n == 0) {
    diag->deepsteps++;
    diag->simtime += dt;
    drift(s, etime, dt);
    clevel--;
    return;
  }
  if (calc_timestep) ok_timestep_cpu(f, dt);
  struct forces slowf = zeroforces, fastf = zeroforces;
  ok_split((FLOAT) dt, f, &slowf, &fastf);
  evolve_ok2(s, fastf, stime, stime+dt/2, dt/2, 0);
  ok_kick(slowf, dt);
  evolve_ok2(s, fastf, stime+dt/2, etime, dt/2, 1);
  clevel--;
}
