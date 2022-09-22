/*
 * Reference integrators with single, global shared time step.
 */

#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "evolve.h"
#include "integrators_shared.h"
// AMUSE STOPPING CONDITIONS SUPPORT
#include <stopcond.h>

static void set_dt_levels(DOUBLE *dt_levels, DOUBLE dt){
    int i;
    dt_levels[0] = dt;
    for (i=1; i<MAXLEVEL; i++) {
        dt_levels[i] = dt_levels[i-1]/2;
    }
}

static DOUBLE time_from_steps(int *steps, DOUBLE *dt_levels) {
    DOUBLE time = 0.0L;
    int i;
    for (i=0; i<MAXLEVEL; i++) {
        time += dt_levels[i] * steps[i];
    }
    return time;
}

static int update_steps_and_get_next_level(int *steps, int current_level) {
    while (current_level > 0) {
        if (steps[current_level] == 0) {
            steps[current_level] = 1;
            break;
        } else {
            steps[current_level] = 0;
            current_level--;
        }
    }
    return current_level;
}

static void detect_collisions(struct sys s) {
    UINT i, j;
    FLOAT dx[3], dr2, radius_sum;
    struct particle *ipart, *jpart;
#pragma omp parallel for if((ULONG) s.n*s.n>MPWORKLIMIT && !omp_in_parallel()) default(none) \
 private(i,j,dx,dr2,radius_sum, ipart, jpart) \
 shared(s)
    for (i=0; i<s.n; i++) {
        ipart=GETPART(s,i);
        for (j=i+1; j<s.n; j++) {
            jpart=GETPART(s,j);
            dx[0] = ipart->pos[0] - jpart->pos[0];
            dx[1] = ipart->pos[1] - jpart->pos[1];
            dx[2] = ipart->pos[2] - jpart->pos[2];
            dr2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
            radius_sum = ipart->radius + jpart->radius;
            if (dr2 <= radius_sum*radius_sum) {
                #pragma omp critical
                {
                    int stopping_index = next_index_for_stopping_condition();
                    if (stopping_index >= 0) {
                        set_stopping_condition_info(stopping_index, COLLISION_DETECTION);
                        set_stopping_condition_particle_index(stopping_index, 0, ipart->id);
                        set_stopping_condition_particle_index(stopping_index, 1, jpart->id);
                    }
                }
            }
        }
    }
}

static void evolve_shared_collision_detection(struct sys s, DOUBLE dt, void (*dkd_func)(int, struct sys, DOUBLE, DOUBLE, DOUBLE)) {
    FLOAT dtsys;
    int next_level, current_level = 0;
    DOUBLE etime, stime;
    DOUBLE *dt_levels = (DOUBLE*) malloc (MAXLEVEL * sizeof(DOUBLE));
    int *step_at_level = (int*) calloc (MAXLEVEL, sizeof(int));
    int is_collision_detection_enabled;
    
    if (dt == 0.0L) {
        ENDRUN("timestep too small: dt=%Le\n", (long double) dt);
    }
    is_stopping_condition_enabled(COLLISION_DETECTION, &is_collision_detection_enabled);
    set_dt_levels(dt_levels, dt);
    do {
        timestep(current_level, s, s, SIGN(dt));
        dtsys = global_timestep(s);
        while (dtsys < fabs(dt_levels[current_level])) {
            current_level++;
            if (current_level >= MAXLEVEL) {
                stime = time_from_steps(step_at_level, dt_levels);
                ENDRUN("timestep too small: stime=%Le dt=%Le clevel=%u\n", 
                    (long double) stime, (long double) dt_levels[current_level], current_level);
            }
        }
        diag->deepsteps++;
        diag->simtime+=dt_levels[current_level];
        stime = time_from_steps(step_at_level, dt_levels);
        next_level = update_steps_and_get_next_level(step_at_level, current_level);
        etime = time_from_steps(step_at_level, dt_levels);
        dkd_func(current_level, s, stime, etime, dt_levels[current_level]);
        if (is_collision_detection_enabled) {
            detect_collisions(s);
            if (set_conditions & enabled_conditions) break;
        }
        current_level = next_level;
    } while (current_level > 0);
    free(dt_levels);
    free(step_at_level);
}

static void kdk_with_zerosys(int clevel, struct sys s, DOUBLE stime, DOUBLE etime, DOUBLE dt) {
    kdk(clevel, s, zerosys, stime, etime, dt);
}

void evolve_shared2_collision_detection(struct sys s, DOUBLE dt) {
    evolve_shared_collision_detection(s, dt, kdk_with_zerosys);
}
void evolve_shared4_collision_detection(struct sys s, DOUBLE dt) {
    evolve_shared_collision_detection(s, dt, dkd4);
}
void evolve_shared6_collision_detection(struct sys s, DOUBLE dt) {
    evolve_shared_collision_detection(s, dt, dkd6);
}
void evolve_shared8_collision_detection(struct sys s, DOUBLE dt) {
    evolve_shared_collision_detection(s, dt, dkd8);
}
void evolve_shared10_collision_detection(struct sys s, DOUBLE dt) {
    evolve_shared_collision_detection(s, dt, dkd10);
}

