#include <stdio.h>
#include "stopcond.h"
#include <string.h>
#include <float.h>
#include <stdlib.h>

static int MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET = 0;
static int DELTA_MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET = 100;

static int * type_of_stopping_condition_set = 0;
static int * index_of_particle_in_stopping_condition = 0;


long enabled_conditions = 0;
long set_conditions = 0;
long supported_conditions = 0;
long number_of_stopping_conditions_set = 0;

double timeout_parameter = 4.0;
double out_of_box_parameter = 0.0;
long number_of_steps_parameter = 1;
double minimum_density_parameter = -1.0;
double maximum_density_parameter = DBL_MAX;
double minimum_internal_energy_parameter = -1.0;
double maximum_internal_energy_parameter = DBL_MAX;
double size_limit_parameter = 0.0;
int use_center_of_mass_parameter = 1;

static int sc_mpi_size;

int enable_stopping_condition(int type) {
    if(type > 32) {
        return -1;
    }
    enabled_conditions |= 1l << type;
    return 0;
}

int enable_stopping_condition_(int *type) {
    return enable_stopping_condition(*type);
}

int increase_memory() {
    int old_max = MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET;
    int i = 0;
    int j = 0;
    void * new_pointer = 0;
    MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET += DELTA_MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET;
    if(type_of_stopping_condition_set == 0) {
        new_pointer = malloc (MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET * sizeof(int));
    } else {
        new_pointer = realloc (type_of_stopping_condition_set, sizeof(int) *  MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET);
    }
    if (!new_pointer) {
        return -1;
    }
    type_of_stopping_condition_set =  (int *) new_pointer;
    
    if(index_of_particle_in_stopping_condition == 0) {
        new_pointer = malloc (MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET * MAX_NUMBER_OF_PARTICLES_PER_INDEX * sizeof(int));
    } else {
        new_pointer = realloc (index_of_particle_in_stopping_condition, sizeof(int) *  MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET * MAX_NUMBER_OF_PARTICLES_PER_INDEX);
    }
    if (!new_pointer) {
        return -1;
    }
    index_of_particle_in_stopping_condition =  (int *) new_pointer;
    
    
    
    
    for(i = old_max; i < MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET; i++) {
        type_of_stopping_condition_set[i] = 0;
        for(j = 0; j < MAX_NUMBER_OF_PARTICLES_PER_INDEX; j++) {
            index_of_particle_in_stopping_condition[i * MAX_NUMBER_OF_PARTICLES_PER_INDEX + j] = 0;
        }
    }
    return 0;
}

//fortran only
//but obsolete, use
//  set_support_for_condition()
//  is_stopping_condition_enabled()
//  is_stopping_condition_set()
//  is_any_condition_set()

int set_supported_conditions_(int *type) {
    supported_conditions = *type;
    return 0;
}

int set_support_for_condition(int type) {
    long mask =  1l << type;
    supported_conditions |= mask;
    return 0;
}

int set_support_for_condition_(int *type) {
    return set_support_for_condition(*type);
}

int get_enabled_conditions_() {
    return enabled_conditions;;
}

int get_set_conditions_() {
    return set_conditions;;
}
//------------


int is_condition_enabled() {
    return 0;
}

int get_number_of_stopping_conditions_set(int * result) {
    *result = number_of_stopping_conditions_set;
    return 0;
}

int get_number_of_stopping_conditions_set_(int *result) {
    return get_number_of_stopping_conditions_set(result);
}

int is_stopping_condition_set(int type, int * result) {
    long mask =  1l << type;
    *result =  (mask & set_conditions) > 0;
    return 0;
}

int is_stopping_condition_set_(int *type, int *result) {
    return is_stopping_condition_set(*type, result);
}

int is_stopping_condition_enabled(int type, int * result) {
    long mask =  1l << type;
    *result = (mask & enabled_conditions) > 0;
    return 0;
}

int is_stopping_condition_enabled_(int *type, int *result) {
    return is_stopping_condition_enabled(*type, result);
}

int is_any_condition_set() {
    if (set_conditions & enabled_conditions) {
	return 1;
    }
    else {
	return 0;
    }
}

int is_any_condition_set_() {
    return is_any_condition_set();
}

int disable_stopping_condition(int type) {
    long mask =  ~(1l << type);
    enabled_conditions &= mask;
    return 0;
}

int disable_stopping_condition_(int *type) {
    return disable_stopping_condition(*type);
}

int has_stopping_condition(int type, int * result) {
    long mask =  1l << type;
    *result = (mask & supported_conditions) > 0;
    return 0;
}

int has_stopping_condition_(int *type, int *result) {
    return has_stopping_condition(*type, result);
}

int get_stopping_condition_info(int index, int * type, int *number_of_particles) {
    int i, id;
    
    if (index >= number_of_stopping_conditions_set) {
        return -1;
    }
    
    *number_of_particles = 0;
    for (i = 0; i < MAX_NUMBER_OF_PARTICLES_PER_INDEX; i++){
        id = index_of_particle_in_stopping_condition[index * MAX_NUMBER_OF_PARTICLES_PER_INDEX + i];
        if (id != -1){
            // Move the particles to the front of the storage - just in case...
            index_of_particle_in_stopping_condition[index * MAX_NUMBER_OF_PARTICLES_PER_INDEX + i] = -1;
            index_of_particle_in_stopping_condition[index * MAX_NUMBER_OF_PARTICLES_PER_INDEX + *number_of_particles] = id;
            
            (*number_of_particles)++;
        }
    }
    
    *type = type_of_stopping_condition_set[index];
    return 0;
}

int get_stopping_condition_info_(int *index, int *type, int *number_of_particles){
    return get_stopping_condition_info(*index, type, number_of_particles);
}

int get_stopping_condition_particle_index(int index, int index_in_the_condition, int *index_of_particle) {
    if(index >= number_of_stopping_conditions_set) {
        return -1;
    }
    
    if(index_in_the_condition >= MAX_NUMBER_OF_PARTICLES_PER_INDEX) {
        return -1;
    }
    
    *index_of_particle = index_of_particle_in_stopping_condition[index * MAX_NUMBER_OF_PARTICLES_PER_INDEX + index_in_the_condition];
    return 0;
}

int get_stopping_condition_particle_index_(int *index, int *index_in_the_condition, int *index_of_particle) {
    return get_stopping_condition_particle_index(*index, *index_in_the_condition, index_of_particle);
}

int reset_stopping_conditions() {
    number_of_stopping_conditions_set = 0;
    set_conditions = 0;
    return 0;
}

int reset_stopping_conditions_() {
    return reset_stopping_conditions();
}

int next_index_for_stopping_condition() {
    int i;
    if (number_of_stopping_conditions_set == MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET) {
        if( increase_memory() ) {
            return -1;
        }
        
    }
    
    // Initialize the particle index storage for this stopping condition
    for (i = 0; i < MAX_NUMBER_OF_PARTICLES_PER_INDEX; i++) {
        index_of_particle_in_stopping_condition[number_of_stopping_conditions_set * MAX_NUMBER_OF_PARTICLES_PER_INDEX + i] = -1;
    }
    
    return number_of_stopping_conditions_set++; // Post increment - returns original value
}

int next_index_for_stopping_condition_() {
    return next_index_for_stopping_condition();
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

int set_stopping_condition_info_(int *index, int *type) {
    return set_stopping_condition_info(*index, *type);
}

int set_stopping_condition_particle_index(int index, int index_in_the_condition, int id_of_the_particle) {
    if (index >= number_of_stopping_conditions_set)
        return -1;
    
    if (index_in_the_condition >= MAX_NUMBER_OF_PARTICLES_PER_INDEX)
        return -1;
    
    index_of_particle_in_stopping_condition[index * MAX_NUMBER_OF_PARTICLES_PER_INDEX + index_in_the_condition] = id_of_the_particle;
    return 0;
}

int set_stopping_condition_particle_index_(int *index, int *index_in_the_condition, int *index_of_the_particle) {
    return set_stopping_condition_particle_index(*index, *index_in_the_condition, *index_of_the_particle);
}

int set_stopping_condition_timeout_parameter(double value) {
    if(value < 0.0) {
	return -1;
    }
    timeout_parameter = value;
    return 0;
}

int set_stopping_condition_timeout_parameter_(double *value) {
    return set_stopping_condition_timeout_parameter(*value);
}

int get_stopping_condition_timeout_parameter(double * value) {
    *value = timeout_parameter;
    return 0;
}

int get_stopping_condition_timeout_parameter_(double * value) {
    return get_stopping_condition_timeout_parameter(value);
}

int set_stopping_condition_number_of_steps_parameter(int value) {
    if (value<1) {
	return -1;
    }
    number_of_steps_parameter = value;
    return 0;
}

int set_stopping_condition_number_of_steps_parameter_(int *value) {
    return set_stopping_condition_number_of_steps_parameter(*value);
}


int get_stopping_condition_number_of_steps_parameter(int *value) {
    *value = number_of_steps_parameter;
    return 0;
}

int get_stopping_condition_number_of_steps_parameter_(int *value) {
    return get_stopping_condition_number_of_steps_parameter(value);
}

int set_stopping_condition_out_of_box_use_center_of_mass_parameter(int value) {
    use_center_of_mass_parameter = value;
    return 0;
}

int set_stopping_condition_out_of_box_use_center_of_mass_parameter_(int *value) {
    return set_stopping_condition_out_of_box_use_center_of_mass_parameter(*value);
}

int get_stopping_condition_out_of_box_use_center_of_mass_parameter(int *value) {
    *value = use_center_of_mass_parameter;
    return 0;
}

int get_stopping_condition_out_of_box_use_center_of_mass_parameter_(int *value) {
    return get_stopping_condition_out_of_box_use_center_of_mass_parameter(value);
}



int set_stopping_condition_out_of_box_parameter(double value) {
    out_of_box_parameter = value;
    return 0;
}

int set_stopping_condition_out_of_box_parameter_(double *value) {
    return set_stopping_condition_out_of_box_parameter(*value);
}

int get_stopping_condition_out_of_box_parameter(double *value) {
    *value = out_of_box_parameter;
    return 0;
}

int get_stopping_condition_out_of_box_parameter_(double *value) {
    return get_stopping_condition_out_of_box_parameter(value);
}




int set_stopping_condition_minimum_density_parameter(double value) {
    minimum_density_parameter = value;
    return 0;
}
int set_stopping_condition_minimum_density_parameter_(double *value) {
    return set_stopping_condition_minimum_density_parameter(*value);
}
int get_stopping_condition_minimum_density_parameter(double *value) {
    *value = minimum_density_parameter;
    return 0;
}
int get_stopping_condition_minimum_density_parameter_(double *value) {
    return get_stopping_condition_minimum_density_parameter(value);
}

int set_stopping_condition_maximum_density_parameter(double value) {
    if (value < 0.0) {
        maximum_density_parameter = DBL_MAX;
    } else {
        maximum_density_parameter = value;
    }
    return 0;
}
int set_stopping_condition_maximum_density_parameter_(double *value) {
    return set_stopping_condition_maximum_density_parameter(*value);
}
int get_stopping_condition_maximum_density_parameter(double *value) {
    *value = maximum_density_parameter;
    return 0;
}
int get_stopping_condition_maximum_density_parameter_(double *value) {
    return get_stopping_condition_maximum_density_parameter(value);
}

int set_stopping_condition_minimum_internal_energy_parameter(double value) {
    minimum_internal_energy_parameter = value;
    return 0;
}
int set_stopping_condition_minimum_internal_energy_parameter_(double *value) {
    return set_stopping_condition_minimum_internal_energy_parameter(*value);
}
int get_stopping_condition_minimum_internal_energy_parameter(double *value) {
    *value = minimum_internal_energy_parameter;
    return 0;
}
int get_stopping_condition_minimum_internal_energy_parameter_(double *value) {
    return get_stopping_condition_minimum_internal_energy_parameter(value);
}

int set_stopping_condition_maximum_internal_energy_parameter(double value) {
    if (value < 0.0) {
        maximum_internal_energy_parameter = DBL_MAX;
    } else {
        maximum_internal_energy_parameter = value;
    }
    return 0;
}
int set_stopping_condition_maximum_internal_energy_parameter_(double *value) {
    return set_stopping_condition_maximum_internal_energy_parameter(*value);
}
int get_stopping_condition_maximum_internal_energy_parameter(double *value) {
    *value = maximum_internal_energy_parameter;
    return 0;
}
int get_stopping_condition_maximum_internal_energy_parameter_(double *value) {
    return get_stopping_condition_maximum_internal_energy_parameter(value);
}




#if defined( MPILIB ) && !defined(NOMPI)

#include <mpi.h>
static int sc_mpi_rank;
static MPI_Comm world = MPI_COMM_WORLD;
static int is_world_set = 0;

int mpi_set_communicator(void * comm) {
    world = *(MPI_Comm *)comm;
    is_world_set = 1;
}

int mpi_setup_stopping_conditions() {
    if(!is_world_set) {
        world = MPI_COMM_WORLD;
        is_world_set = 1;
    }
    int error;
    error = MPI_Comm_rank(world, &sc_mpi_rank);
    if(error) {
        return -1;
    }
    error = MPI_Comm_size(world, &sc_mpi_size);
    if(error) {
        return -1;
    }
    return 0;
}



int mpi_distribute_stopping_conditions() {
    if(sc_mpi_size <= 1) {
        return 0;
    }
    if(!enabled_conditions) {
        return 0;
    }
    return 0;
}



static int * local_type_of_stopping_condition_set = 0;
static int * local_index_of_particle_in_stopping_condition = 0;


int mpi_collect_stopping_conditions() {
    if(!enabled_conditions) {return 0;}
    int i;
    int local_number_of_stopping_conditions_set;
    int counts[sc_mpi_size];
    int displs[sc_mpi_size];
    
    memset(counts, 0, sizeof(counts));
    memset(displs, 0, sizeof(displs));
    long set = 0;
    
    MPI_Allreduce(&set_conditions, &set, 1, MPI_INTEGER,  MPI_BOR, world);
    
    set_conditions = set;
    
    MPI_Gather(&number_of_stopping_conditions_set, 1, MPI_INTEGER, counts, 1, MPI_INT, 0, world);
    
    if(sc_mpi_rank == 0) {
        local_number_of_stopping_conditions_set = 0;
        for(i = 0; i < sc_mpi_size; i++) {
            local_number_of_stopping_conditions_set += counts[i];
        }
       
        if(local_number_of_stopping_conditions_set > MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET) {
            int n = ((local_number_of_stopping_conditions_set - MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET) / DELTA_MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET) + 1;
            int tmp = DELTA_MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET;
            DELTA_MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET = n * DELTA_MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET;
            increase_memory();
            DELTA_MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET = tmp;
        }
        local_type_of_stopping_condition_set = (int *) calloc (MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET , sizeof(int));;
        local_index_of_particle_in_stopping_condition = (int *) calloc (MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET , MAX_NUMBER_OF_PARTICLES_PER_INDEX * sizeof(int));
    }
    if(sc_mpi_rank == 0) {
        int x = 0;
        for(i = 0; i < sc_mpi_size; i++) {
            displs[i] = x;
            x += counts[i];
        }
    }

    MPI_Gatherv(
        type_of_stopping_condition_set, 
        number_of_stopping_conditions_set, 
        MPI_INTEGER,
        local_type_of_stopping_condition_set, 
        counts, 
        displs, 
        MPI_INTEGER,
        0,
        world
    );

    if(sc_mpi_rank == 0) {
        int x = 0;
        for(i = 0; i < sc_mpi_size; i++) {
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
        counts, 
        displs,
        MPI_INTEGER,
        0, 
        world
    );

    if(sc_mpi_rank == 0) {
        number_of_stopping_conditions_set = local_number_of_stopping_conditions_set;
        if(number_of_stopping_conditions_set > 0) {
            memcpy(index_of_particle_in_stopping_condition, local_index_of_particle_in_stopping_condition, MAX_NUMBER_OF_PARTICLES_PER_INDEX * MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET * sizeof(int));
            memcpy(type_of_stopping_condition_set, local_type_of_stopping_condition_set, MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET * sizeof(int));
        }
        free(local_type_of_stopping_condition_set);
        free(local_index_of_particle_in_stopping_condition);
        local_type_of_stopping_condition_set = 0;
        local_index_of_particle_in_stopping_condition = 0;
    }
}

#else

int mpi_set_communicator(void * comm) {return 0;}
int mpi_setup_stopping_conditions() {sc_mpi_size = 1; return 0;}
int mpi_collect_stopping_conditions() {return 0;}
int mpi_distribute_stopping_conditions() {return 0;}

#endif



int mpi_setup_stopping_conditions_() {
    return mpi_setup_stopping_conditions();
}

int mpi_distribute_stopping_conditions_() {
    return mpi_distribute_stopping_conditions();
}

int mpi_collect_stopping_conditions_() {
    return mpi_collect_stopping_conditions();
}
