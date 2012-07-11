module StoppingConditions

    implicit none
    
    integer :: MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET
    parameter(MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET = 256)
    
    integer :: MAX_NUMBER_OF_PARTICLES_PER_INDEX
    parameter(MAX_NUMBER_OF_PARTICLES_PER_INDEX = 256)
    
    double precision :: DBL_MAX
    parameter(DBL_MAX =  HUGE(0.0d0))
    
    integer :: type_of_stopping_condition_set(0:MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET)
    integer :: index_of_particle_in_stopping_condition(&
&      0:MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET,&
&      0:MAX_NUMBER_OF_PARTICLES_PER_INDEX)

    integer*8 :: enabled_conditions = 0
    integer*8 :: set_conditions = 0
    integer*8 :: supported_conditions = 0
    integer*8 :: number_of_stopping_conditions_set = 0

    double precision :: timeout_parameter = 4.0
    double precision :: out_of_box_parameter = 0.0
    integer*8 :: number_of_steps_parameter = 1
    double precision :: minimum_density_parameter = -1.0
    double precision :: maximum_density_parameter = DBL_MAX
    double precision :: minimum_internal_energy_parameter = -1.0
    double precision :: maximum_internal_energy_parameter = DBL_MAX
    integer :: sc_mpi_size;

contains
    
    
function enable_stopping_condition(type) 
    implicit none
    integer, intent(in) :: type
    integer :: enable_stopping_condition
    
    if(type > 32) then
        enable_stopping_condition = -1
        return
    end if
    enabled_conditions = IOR(enabled_conditions, ISHFT(1 , type))
    enable_stopping_condition = 0
end function

function set_support_for_condition(type) 
    implicit none
    integer, intent(in) :: type
    integer :: set_support_for_condition

    supported_conditions = IOR(supported_conditions, ISHFT(1 , type))
    set_support_for_condition = 0
end function

function get_enabled_conditions() 
    implicit none
    integer :: get_enabled_conditions
    get_enabled_conditions = enabled_conditions
end function

function get_set_conditions() 
    implicit none
    integer :: get_set_conditions
    get_set_conditions = set_conditions
end function

function get_number_of_stopping_conditions_set(result)
    implicit none
    integer, intent(out) :: result
    integer :: get_number_of_stopping_conditions_set
    
    result = number_of_stopping_conditions_set;
    get_number_of_stopping_conditions_set = 0
end function

function is_stopping_condition_set(type, result)
    implicit none
    integer, intent(in) :: type
    integer, intent(out) :: result
    integer :: is_stopping_condition_set
    
    if(IAND(ISHFT(1 , type) , set_conditions) > 0) then
        result = 1
    else
        result = 0
    end if
    is_stopping_condition_set = 0
end function

function is_stopping_condition_enabled(type,  result)
    implicit none
    integer, intent(in) :: type
    integer, intent(out) :: result
    integer :: is_stopping_condition_enabled
    
    if(IAND(ISHFT(1 , type) , enabled_conditions) > 0) then
        result = 1
    else
        result = 0
    end if
    is_stopping_condition_enabled = 0
end function


function is_any_condition_set()
    implicit none
    integer :: is_any_condition_set
    if (IAND(set_conditions ,enabled_conditions) > 0) then
	is_any_condition_set = 1
    else
	is_any_condition_set = 0
    end if
end function

function disable_stopping_condition(type)
    implicit none
    integer, intent(in) :: type
    integer :: disable_stopping_condition
    enabled_conditions = IAND(enabled_conditions, NOT(ISHFT(1 , type)));
    disable_stopping_condition = 0
end function

function has_stopping_condition(type, result)
    implicit none
    integer, intent(in) :: type
    integer, intent(out) :: result
    integer :: has_stopping_condition
    if(IAND(ISHFT(1 , type) , supported_conditions) > 0) then
        result = 1
    else
        result = 0
    end if
    has_stopping_condition = 0;
end function

function get_stopping_condition_info(index, type, number_of_particles)
    implicit none
    integer, intent(in) :: index
    integer, intent(out) :: type
    integer, intent(out) :: number_of_particles
    integer :: get_stopping_condition_info
    integer :: i, id;
    
    if (index >= number_of_stopping_conditions_set) then
        get_stopping_condition_info = -1
        return
    end if
    
    number_of_particles = 1    
    do i = 0, MAX_NUMBER_OF_PARTICLES_PER_INDEX
        id = index_of_particle_in_stopping_condition(index, i);
        if (id >= 0) then
            number_of_particles = number_of_particles + 1;
        else
            exit
        end if
    end do
    
    type = type_of_stopping_condition_set(index);
    get_stopping_condition_info = 0
    return
end function

function get_stopping_condition_particle_index(index, index_in_the_condition, index_of_particle) 
    implicit none
    integer, intent(in) :: index
    integer, intent(in) :: index_in_the_condition
    integer, intent(out) :: index_of_particle
    integer :: get_stopping_condition_particle_index
    integer :: i, id;
    
    if (index >= number_of_stopping_conditions_set) then
        get_stopping_condition_particle_index = -1
        return
    end if
    
    if (index_in_the_condition >= MAX_NUMBER_OF_PARTICLES_PER_INDEX) then
        get_stopping_condition_particle_index = -1
        return
    end if
    
    
    index_of_particle = index_of_particle_in_stopping_condition(index, index_in_the_condition + 1)
    get_stopping_condition_particle_index = 0
end function


function reset_stopping_conditions() 
    implicit none
    integer :: reset_stopping_conditions
    number_of_stopping_conditions_set = 0
    set_conditions = 0
    reset_stopping_conditions = 0
end function

function next_index_for_stopping_condition()
    implicit none
    integer :: next_index_for_stopping_condition, i
    
    if (number_of_stopping_conditions_set == MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET) then
        next_index_for_stopping_condition = -1
	return
    end if
    
    do i = 0 , MAX_NUMBER_OF_PARTICLES_PER_INDEX
        index_of_particle_in_stopping_condition(number_of_stopping_conditions_set, i) = -1
    end do
    
    next_index_for_stopping_condition = number_of_stopping_conditions_set
    number_of_stopping_conditions_set = number_of_stopping_conditions_set + 1
end function



function set_stopping_condition_info(index, type)
    implicit none
    integer, intent(in) :: index
    integer, intent(in) :: type
    integer :: set_stopping_condition_info
    
    if (index >= number_of_stopping_conditions_set) then
        set_stopping_condition_info = -1
	return
    end if
    
    set_conditions = IOR(set_conditions, ISHFT(1 , type))

    type_of_stopping_condition_set(index) = type;
    set_stopping_condition_info = 0
end function


function set_stopping_condition_particle_index(index, index_in_the_condition, id_of_the_particle)
    implicit none
    integer, intent(in) :: index
    integer, intent(in) :: index_in_the_condition
    integer, intent(in) :: id_of_the_particle
    integer :: set_stopping_condition_particle_index
    
    if (index >= number_of_stopping_conditions_set) then
        set_stopping_condition_particle_index = -1
	return
    end if
    
    if (index_in_the_condition >= MAX_NUMBER_OF_PARTICLES_PER_INDEX) then
        set_stopping_condition_particle_index = -1
	return
    end if
    
    index_of_particle_in_stopping_condition(index, index_in_the_condition) = id_of_the_particle
    set_stopping_condition_particle_index = 0
end function

#if 0
//------------


int is_condition_enabled() {
    return 0;
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

int mpi_setup_stopping_conditions() {
    int error;
    error = MPI_Comm_rank(MPI_COMM_WORLD, &sc_mpi_rank);
    if(error) {
	return -1;
    }
    error = MPI_Comm_size(MPI_COMM_WORLD, &sc_mpi_size);
    if(error) {
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

static int local_type_of_stopping_condition_set[MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET];
static int local_index_of_particle_in_stopping_condition[MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET*MAX_NUMBER_OF_PARTICLES_PER_INDEX];

int mpi_collect_stopping_conditions() {
    if(!enabled_conditions) {return 0;}
    int i;
    int local_number_of_stopping_conditions_set;
    int counts[sc_mpi_size];
    int displs[sc_mpi_size];
    memset(counts, 0, sizeof(counts));
    memset(displs, 0, sizeof(displs));
    long set = 0;
    
    MPI_Allreduce(&set_conditions, &set, 1, MPI_INTEGER,  MPI_BOR, MPI_COMM_WORLD);
    
    set_conditions = set;
    
    MPI_Gather(&number_of_stopping_conditions_set, 1, MPI_INTEGER, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if(sc_mpi_rank == 0) {
	    local_number_of_stopping_conditions_set = 0;
        for(i = 0; i < sc_mpi_size; i++) {
            local_number_of_stopping_conditions_set += counts[i];
        }
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
	    MPI_COMM_WORLD
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
        MPI_COMM_WORLD
    );

    if(sc_mpi_rank == 0) {
        number_of_stopping_conditions_set = local_number_of_stopping_conditions_set;
	    memcpy(index_of_particle_in_stopping_condition, local_index_of_particle_in_stopping_condition, sizeof(index_of_particle_in_stopping_condition));
	    memcpy(type_of_stopping_condition_set, local_type_of_stopping_condition_set, sizeof(type_of_stopping_condition_set));
    }
}

#else

int mpi_setup_stopping_conditions() {sc_mpi_size = 1; return 0;}
int mpi_collect_stopping_conditions() {return 0;}
int mpi_distribute_stopping_conditions() {return 0;}

#endif
#endif

end module
