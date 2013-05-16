/*
 * Reference integrators with single, global shared time step.
 */

#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>
#include "evolve.h"
#include "integrators_shared.h"


void set_dt_levels(DOUBLE *dt_levels, DOUBLE dt){
    int i;
    dt_levels[0] = dt;
    for (i=1; i<MAXLEVEL; i++) {
        dt_levels[i] = dt_levels[i-1]/2;
    }
}

DOUBLE time_from_steps(int *steps, DOUBLE *dt_levels) {
    DOUBLE time = 0.0L;
    int i;
    for (i=0; i<MAXLEVEL; i++) {
        time += dt_levels[i] * steps[i];
    }
    return time;
}

int update_steps_and_get_next_level(int *steps, int current_level) {
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

void evolve_shared6_nonrecursive(struct sys s, DOUBLE dt) {
    FLOAT dtsys;
    int next_level, current_level = 0;
    DOUBLE etime, stime;
    DOUBLE *dt_levels = (DOUBLE*) malloc (MAXLEVEL * sizeof(DOUBLE));
    int *step_at_level = (int*) calloc (MAXLEVEL, sizeof(int));
    
    printf("evolve_shared6_nonrecursive, dt: %Le \n", dt);
    if (dt == 0.0L) {
        ENDRUN("timestep too small: dt=%Le\n", dt);
    }
    set_dt_levels(dt_levels, dt);
    do {
        timestep(current_level, s, s, SIGN(dt));
        dtsys = global_timestep(s);
        while (dtsys < fabs(dt_levels[current_level])) {
            current_level++;
            if (current_level >= MAXLEVEL) {
                stime = time_from_steps(step_at_level, dt_levels);
                ENDRUN("timestep too small: stime=%Le dt=%Le clevel=%u\n", 
                    stime, dt_levels[current_level], current_level);
            }
        }
        diag->deepsteps++;
        diag->simtime+=dt_levels[current_level];
        stime = time_from_steps(step_at_level, dt_levels);
        next_level = update_steps_and_get_next_level(step_at_level, current_level);
        etime = time_from_steps(step_at_level, dt_levels);
        dkd6(current_level, s, stime, etime, dt_levels[current_level]);
        current_level = next_level;
    } while (current_level > 0);
    free(dt_levels);
    free(step_at_level);
    
}

