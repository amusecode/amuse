#include "amuse_mpi.h"

static MPI_Comm local_world;
static int was_set = 0;

int set_comm_world(MPI_Comm comm) {
    local_world = comm;
    was_set = 1;
    return 0;
}

int get_comm_world(MPI_Comm * comm) {
    if(!was_set) {
        local_world = MPI_COMM_WORLD;
        was_set = 1;
    }
    *comm = local_world;
    return 0;
}
