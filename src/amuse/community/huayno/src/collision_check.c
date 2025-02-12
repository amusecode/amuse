/*
 * Algorithm to detect collisions between particles
 */
#include "evolve.h"
#include "evolve_ok.h"
#include "collision_check.h"
// AMUSE STOPPING CONDITIONS SUPPORT
#include <stopcond.h>


void detect_collisions(struct sys s){
    UINT i, j;
    FLOAT dx[3], dr2, radius_sum;
    struct particle *ipart, *jpart;
    #pragma omp parallel for if((ULONG) s.n*s.n>MPWORKLIMIT && !omp_in_parallel()) default(none) \
    private(i,j,dx,dr2,radius_sum, ipart, jpart) \
    shared(s)
    for (i=0; i<s.n; i++) {
        ipart=GETPART(s,i);
        if (ipart->mass > 0){
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
}