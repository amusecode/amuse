#ifndef NOMPI
#include <mpi.h>
#endif

#include <stopcond.h>
#ifdef __cplusplus
extern "C" {
#endif
int initialize_code() {
    // AMUSE STOPPING CONDITIONS SUPPORT
    supported_conditions = COLLISION_DETECTION_BITMAP | PAIR_DETECTION_BITMAP | TIMEOUT_DETECTION_BITMAP | OUT_OF_BOX_DETECTION_BITMAP;
    // -----------------------
    return 0;
}


int fire_condition(int condition_to_set, int particle_index1, int particle_index2, int rank) {
    int my_rank;
    int error, stopping_index;

#ifndef NOMPI
    error = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
    error = 0;
    my_rank = rank;
#endif
    if (rank >= 0 && rank != my_rank) { return 0; }

    stopping_index = next_index_for_stopping_condition();

    error = set_stopping_condition_info(stopping_index, condition_to_set);
    if(particle_index1 > 0) {
        error = set_stopping_condition_particle_index(stopping_index, 0, particle_index1);
    }
    if(particle_index2 >  0) {
        error = set_stopping_condition_particle_index(stopping_index, 1, particle_index2);
    }
    return 0;
}
#ifdef __cplusplus
}
#endif

