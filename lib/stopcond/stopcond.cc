#include "stopcond.h"
#include <iostream>
#include <string.h>

int type_of_stopping_condition_set[MAX_NUMBER_OF_SIMULTANIOS_CONDITIONS_SET];
int index_of_particle_in_stopping_condition[MAX_NUMBER_OF_SIMULTANIOS_CONDITIONS_SET * MAX_NUMBER_OF_PARTICLES_PER_INDEX];

long enabled_conditions = 0;
long set_conditions = 0;
long number_of_stopping_conditions_set = 0;

double timeout_parameter = 4.0;

int enable_stopping_condition(int type)
{
    if(type > 32) {
        return -1;
    }
    enabled_conditions |= 1l << type;
    
    //std::cerr<< enabled_conditions << std::endl;
    //std::cerr<< type<< std::endl;
    return 0;
}

int get_number_of_stopping_conditions_set(int * result)
{
    *result = number_of_stopping_conditions_set;
    
    return 0;
}

int is_stopping_condition_set(int type, int * result)
{
    long mask =  1l << type;
    *result =  (mask & set_conditions) > 0;
    return 0;
}

int is_stopping_condition_enabled(int type, int * result)
{   
    long mask =  1l << type;
    *result = (mask & enabled_conditions) > 0;
    return 0;
}

int disable_stopping_condition(int type)
{
    long mask =  ~(1l << type);
    enabled_conditions &= mask;
    return 0;
}

int has_stopping_condition(int type, int * result)
{
    long mask =  1l << type;
    *result = (mask & supported_conditions) > 0;
    return 0;
}

int get_stopping_condition_info(int index, int * type)
{
    if(index >= number_of_stopping_conditions_set) {
        return -1;
    }
    
    *type = type_of_stopping_condition_set[index];

    return 0;
}

int get_stopping_condition_particle_index(
    int index, 
    int index_in_the_condition, 
    int * index_of_particle)
{
    if(index >= number_of_stopping_conditions_set) {
        return -1;
    }
    if(index_in_the_condition >= MAX_NUMBER_OF_PARTICLES_PER_INDEX) {
        return -1;
    }
    *index_of_particle = index_of_particle_in_stopping_condition[index * MAX_NUMBER_OF_PARTICLES_PER_INDEX + index_in_the_condition];

    return 0;
}


int reset_stopping_conditions() {
    number_of_stopping_conditions_set = 0;
    set_conditions = 0;
}

int next_index_for_stopping_condition() {
    number_of_stopping_conditions_set++;
    return number_of_stopping_conditions_set - 1;
}

int set_stopping_condition_info(int index, int type) {
    if(index >= number_of_stopping_conditions_set) {
        return -1;
    }
    long mask =  1l << type;
    set_conditions |= mask;
    
    type_of_stopping_condition_set[index] = type;
    return 0;

}

int set_stopping_condition_particle_index(int index, int index_in_the_condition, int index_of_the_particle) {
    if(index >= number_of_stopping_conditions_set) {
        return -1;
    }
    if(index_in_the_condition >= MAX_NUMBER_OF_PARTICLES_PER_INDEX -1 ) {
        return -1;
    }
    index_of_particle_in_stopping_condition[index * MAX_NUMBER_OF_PARTICLES_PER_INDEX + index_in_the_condition] = index_of_the_particle;
    index_of_particle_in_stopping_condition[index * MAX_NUMBER_OF_PARTICLES_PER_INDEX + index_in_the_condition + 1] = -1;
    return 0;
}


int set_stopping_condition_timeout_parameter(double value) {
    if(timeout_parameter < 0.0) {
        return -1;
    }
    timeout_parameter = value;
    return 0;
}

int get_stopping_condition_timeout_parameter(double * value){
    *value = timeout_parameter;
    return 0;
}


#ifdef MPILIB

#include <mpi.h>
int sc_mpi_rank;
int sc_mpi_size;

int mpi_setup_stopping_conditions() {
    int error;
    error = MPI_Comm_rank(MPI_COMM_WORLD, &sc_mpi_rank);
    if(error) {
        std::cerr << error << std::endl;
        return -1;
    }
    error = MPI_Comm_size(MPI_COMM_WORLD, &sc_mpi_size);
    if(error) {
        std::cerr << error << std::endl;
        return -1;
    }
    return 0;
}
int mpi_distribute_stopping_conditions() {
    if(sc_mpi_size <= 1) {
        return 0;
    }
    if(!enabled_conditions) {return 0;}
    
}

int local_type_of_stopping_condition_set[MAX_NUMBER_OF_SIMULTANIOS_CONDITIONS_SET];
int local_index_of_particle_in_stopping_condition[MAX_NUMBER_OF_SIMULTANIOS_CONDITIONS_SET*MAX_NUMBER_OF_PARTICLES_PER_INDEX];


int mpi_collect_stopping_conditions() {
    if(!enabled_conditions) {return 0;}
    
    int counts[sc_mpi_size];
    int displs[sc_mpi_size];
    memset(counts, 0, sizeof(counts));
    memset(displs, 0, sizeof(displs));
    long set = 0;
    MPI_Reduce(&set_conditions, &set, 1, MPI_LONG,  MPI_BOR, 0, MPI_COMM_WORLD);
    set_conditions = set;
    MPI_Gather(&number_of_stopping_conditions_set, 1, MPI_INTEGER, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(sc_mpi_rank == 0) {
        number_of_stopping_conditions_set = 0;
        for(int i = 0; i < sc_mpi_size; i++) {
            number_of_stopping_conditions_set += counts[i];
        } 
    }
    if(sc_mpi_rank == 0) {
        int x = 0;
        for(int i = 0; i < sc_mpi_size; i++) {
            displs[i] = x;
            x += counts[i];
        }
    }
    
    MPI_Gatherv(
        type_of_stopping_condition_set, number_of_stopping_conditions_set, MPI_INTEGER, 
        local_type_of_stopping_condition_set, counts, displs, MPI_INTEGER, 
        0, MPI_COMM_WORLD);
    
    if(sc_mpi_rank == 0) {
        int x = 0;
        for(int i = 0; i < sc_mpi_size; i++) {
            displs[i] = x;
            counts[i] *= MAX_NUMBER_OF_PARTICLES_PER_INDEX;
            x += counts[i];
        }
    }
    
    MPI_Gatherv(
        index_of_particle_in_stopping_condition, 
        number_of_stopping_conditions_set*MAX_NUMBER_OF_PARTICLES_PER_INDEX, 
        MPI_INTEGER, 
        local_index_of_particle_in_stopping_condition, 
        counts, displs, 
        MPI_INTEGER, 
        0, MPI_COMM_WORLD);
    
    if(sc_mpi_rank == 0) {
        memcpy(index_of_particle_in_stopping_condition, local_index_of_particle_in_stopping_condition, sizeof(index_of_particle_in_stopping_condition));
        memcpy(type_of_stopping_condition_set, local_type_of_stopping_condition_set, sizeof(type_of_stopping_condition_set));
    }
}

#else

int mpi_setup_stopping_conditions() {return 0;}
int mpi_collect_stopping_conditions() {return 0;}
int mpi_distribute_stopping_conditions() {return 0;}

#endif

