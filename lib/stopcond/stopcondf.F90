module StoppingConditions


    implicit none
    
    integer :: MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET
    parameter(MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET = 256)
    
    integer :: MAX_NUMBER_OF_PARTICLES_PER_INDEX
    parameter(MAX_NUMBER_OF_PARTICLES_PER_INDEX = 256)
    
    double precision :: DBL_MAX
    parameter(DBL_MAX = 1.d308)
    
    integer, private, save :: type_of_stopping_condition_set(0:MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET)
    integer, private, save :: index_of_particle_in_stopping_condition(&
&      0:MAX_NUMBER_OF_PARTICLES_PER_INDEX,&
&      0:MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET)

    integer, private, save :: local_type_of_stopping_condition_set(0:MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET)
    integer, private, save :: local_index_of_particle_in_stopping_condition(&
&      0:MAX_NUMBER_OF_PARTICLES_PER_INDEX,&
&      0:MAX_NUMBER_OF_SIMULTANIOUS_CONDITIONS_SET)

    
    integer*4, private, save :: enabled_conditions = 0
    integer*4, private, save :: set_conditions = 0
    integer*4, private, save :: supported_conditions = 0
    integer*4, private, save :: number_of_stopping_conditions_set = 0

    double precision, private, save :: timeout_parameter = 4.0
    double precision, private, save :: out_of_box_parameter = 0.0
    integer*4, private, save :: number_of_steps_parameter = 1
    double precision, private, save :: minimum_density_parameter = -1.0
    double precision, private, save :: maximum_density_parameter = DBL_MAX
    double precision, private, save :: minimum_internal_energy_parameter = -1.0
    double precision, private, save :: maximum_internal_energy_parameter = DBL_MAX
    double precision, private, save :: size_limit_parameter = 0.0
    
    integer, private, save :: sc_mpi_size;
    integer, private, save :: sc_mpi_rank;
    
    INTEGER :: COLLISION_DETECTION_BITMAP
    INTEGER :: PAIR_DETECTION_BITMAP
    INTEGER :: ESCAPER_DETECTION_BITMAP
    INTEGER :: TIMEOUT_DETECTION_BITMAP
    INTEGER :: NUMBER_OF_STEPS_DETECTION_BITMAP
    INTEGER :: OUT_OF_BOX_DETECTION_BITMAP
    INTEGER :: DENSITY_LIMIT_DETECTION_BITMAP
    INTEGER :: INTERNAL_ENERGY_LIMIT_DETECTION_BITMAP
    INTEGER :: SIZE_LIMIT_DETECTION_BITMAP

    INTEGER :: COLLISION_DETECTION
    INTEGER :: PAIR_DETECTION
    INTEGER :: ESCAPER_DETECTION
    INTEGER :: TIMEOUT_DETECTION
    INTEGER :: NUMBER_OF_STEPS_DETECTION
    INTEGER :: OUT_OF_BOX_DETECTION
    INTEGER :: DENSITY_LIMIT_DETECTION
    INTEGER :: INTERNAL_ENERGY_LIMIT_DETECTION
    INTEGER :: SIZE_LIMIT_DETECTION

    PARAMETER (COLLISION_DETECTION_BITMAP=1)
    PARAMETER (PAIR_DETECTION_BITMAP=2)
    PARAMETER (ESCAPER_DETECTION_BITMAP=4)
    PARAMETER (TIMEOUT_DETECTION_BITMAP=8)
    PARAMETER (NUMBER_OF_STEPS_DETECTION_BITMAP=16)
    PARAMETER (OUT_OF_BOX_DETECTION_BITMAP=32)
    PARAMETER (DENSITY_LIMIT_DETECTION_BITMAP=64)
    PARAMETER (INTERNAL_ENERGY_LIMIT_DETECTION_BITMAP=128)
    PARAMETER (SIZE_LIMIT_DETECTION_BITMAP=256)

    PARAMETER (COLLISION_DETECTION=0)
    PARAMETER (PAIR_DETECTION=1)
    PARAMETER (ESCAPER_DETECTION=2)
    PARAMETER (TIMEOUT_DETECTION=3)
    PARAMETER (NUMBER_OF_STEPS_DETECTION=4)
    PARAMETER (OUT_OF_BOX_DETECTION=5)
    PARAMETER (DENSITY_LIMIT_DETECTION=6)
    PARAMETER (INTERNAL_ENERGY_LIMIT_DETECTION=7)
    PARAMETER (SIZE_LIMIT_DETECTION=8)


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
    
    number_of_particles = 0
    do i = 0, MAX_NUMBER_OF_PARTICLES_PER_INDEX
        id = index_of_particle_in_stopping_condition(i, index);
        if (id >= 0) then
        
            index_of_particle_in_stopping_condition(i, index) = -1;
            index_of_particle_in_stopping_condition(number_of_particles, index) = id;
            number_of_particles = number_of_particles + 1;
            
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
    
    
    index_of_particle = index_of_particle_in_stopping_condition(index_in_the_condition, index)
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
        index_of_particle_in_stopping_condition(i, number_of_stopping_conditions_set) = -1
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
    
    index_of_particle_in_stopping_condition(index_in_the_condition, index) = id_of_the_particle
    set_stopping_condition_particle_index = 0
end function



function set_stopping_condition_timeout_parameter(value) 
    double precision, intent(in) :: value
    integer :: set_stopping_condition_timeout_parameter
    if(value < 0.0) then
	set_stopping_condition_timeout_parameter = -1
	return
    end if
    timeout_parameter = value
    set_stopping_condition_timeout_parameter = 0
end function

function get_stopping_condition_timeout_parameter(value)
    integer :: get_stopping_condition_timeout_parameter
    double precision, intent(out) :: value
    value = timeout_parameter
    get_stopping_condition_timeout_parameter = 0
end function

function set_stopping_condition_number_of_steps_parameter(value) 
    integer, intent(in) :: value
    integer :: set_stopping_condition_number_of_steps_parameter
    if(value < 0.0) then
	set_stopping_condition_number_of_steps_parameter = -1
	return
    end if
    number_of_steps_parameter = value
    set_stopping_condition_number_of_steps_parameter = 0
end function

function get_stopping_condition_number_of_steps_parameter(value)
    integer :: get_stopping_condition_number_of_steps_parameter
    integer, intent(out) :: value
    value = number_of_steps_parameter
    get_stopping_condition_number_of_steps_parameter = 0
end function

function set_stopping_condition_out_of_box_parameter(value) 
    double precision, intent(in) :: value
    integer :: set_stopping_condition_out_of_box_parameter
    out_of_box_parameter = value
    set_stopping_condition_out_of_box_parameter = 0
end function

function get_stopping_condition_out_of_box_parameter(value) 
    integer :: get_stopping_condition_out_of_box_parameter
    double precision, intent(out) :: value
    value = out_of_box_parameter
    get_stopping_condition_out_of_box_parameter = 0
end function

function set_stopping_condition_minimum_density_parameter(value) 
    double precision, intent(in) :: value
    integer :: set_stopping_condition_minimum_density_parameter
    minimum_density_parameter = value
    set_stopping_condition_minimum_density_parameter = 0
end function

function get_stopping_condition_minimum_density_parameter(value) 
    integer :: get_stopping_condition_minimum_density_parameter
    double precision, intent(out) :: value
    value = minimum_density_parameter
    get_stopping_condition_minimum_density_parameter = 0
end function

function set_stopping_condition_maximum_density_parameter(value) 
    double precision, intent(in) :: value
    integer :: set_stopping_condition_maximum_density_parameter
    maximum_density_parameter = value
    
    if(value < 0.0) then
	maximum_density_parameter = DBL_MAX
    else
	maximum_density_parameter = value
    end if
    set_stopping_condition_maximum_density_parameter = 0
end function

function get_stopping_condition_maximum_density_parameter(value) 
    integer :: get_stopping_condition_maximum_density_parameter
    double precision, intent(out) :: value
    value = maximum_density_parameter
    get_stopping_condition_maximum_density_parameter = 0
end function



function set_stopping_condition_minimum_internal_energy_parameter(value) 
    double precision, intent(in) :: value
    integer :: set_stopping_condition_minimum_internal_energy_parameter
    minimum_internal_energy_parameter = value
    set_stopping_condition_minimum_internal_energy_parameter = 0
end function

function get_stopping_condition_minimum_internal_energy_parameter(value) 
    integer :: get_stopping_condition_minimum_internal_energy_parameter
    double precision, intent(out) :: value
    value = minimum_internal_energy_parameter
    get_stopping_condition_minimum_internal_energy_parameter = 0
end function

function set_stopping_condition_maximum_internal_energy_parameter(value) 
    double precision, intent(in) :: value
    integer :: set_stopping_condition_maximum_internal_energy_parameter
    maximum_internal_energy_parameter = value
    
    if(value < 0.0) then
	maximum_internal_energy_parameter = DBL_MAX
    else
	maximum_internal_energy_parameter = value
    end if
    set_stopping_condition_maximum_internal_energy_parameter = 0
end function

function get_stopping_condition_maximum_internal_energy_parameter(value) 
    integer :: get_stopping_condition_maximum_internal_energy_parameter
    double precision, intent(out) :: value
    value = maximum_internal_energy_parameter
    get_stopping_condition_maximum_internal_energy_parameter = 0
end function




function set_stopping_condition_size_limit_parameter(value) 
    double precision, intent(in) :: value
    integer :: set_stopping_condition_size_limit_parameter
    size_limit_parameter = value    
    set_stopping_condition_size_limit_parameter = 0
end function

function get_stopping_condition_size_limit_parameter(value) 
    integer :: get_stopping_condition_size_limit_parameter
    double precision, intent(out) :: value
    value = size_limit_parameter
    get_stopping_condition_size_limit_parameter = 0
end function



#if defined( MPILIB ) && !defined(NOMPI)


function mpi_setup_stopping_conditions()

#if defined( MPILIB ) && !defined(NOMPI)
    use mpi
#endif
    implicit none
    
    integer :: mpi_setup_stopping_conditions
    integer :: error
    
    call MPI_Comm_rank(MPI_COMM_WORLD, sc_mpi_rank, error)
    if (error .NE. 0) then
        mpi_setup_stopping_conditions = -1
        return
    end if
    call MPI_Comm_size(MPI_COMM_WORLD, sc_mpi_size, error)
    if (error .NE. 0) then
        mpi_setup_stopping_conditions = -1
        return
    end if
    mpi_setup_stopping_conditions = 0
end function

function mpi_collect_stopping_conditions()
#if defined( MPILIB ) && !defined(NOMPI)
    use mpi
#endif
    implicit none
    integer :: mpi_collect_stopping_conditions
    mpi_collect_stopping_conditions = 0
end function

function mpi_distribute_stopping_conditions()
#if defined( MPILIB ) && !defined(NOMPI)
    use mpi
#endif
    implicit none
    integer :: mpi_distribute_stopping_conditions
    integer :: i, error, offset
    integer :: local_number_of_stopping_conditions_set
    integer :: counts(sc_mpi_size)
    integer :: displs(sc_mpi_size)
    integer*4 :: set
    mpi_distribute_stopping_conditions = 0
    counts = 0
    displs = 0
    set = 0
    
    call MPI_Allreduce(set_conditions, set, 1, MPI_INTEGER,  MPI_BOR, MPI_COMM_WORLD, error)
    
    set_conditions = set
    
    call MPI_Gather(number_of_stopping_conditions_set, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, error)
    
    if (sc_mpi_rank .EQ. 0) then
        local_number_of_stopping_conditions_set = 0
    
        do i = 1, sc_mpi_size
            local_number_of_stopping_conditions_set = local_number_of_stopping_conditions_set + counts(i)
        end do
    end if
    
    if (sc_mpi_rank .EQ. 0) then
        offset = 0
        do i = 1, sc_mpi_size
            displs(i) = offset;
            offset = offset + counts(i)
        end do
    end if
    
    call MPI_Gatherv( &
        type_of_stopping_condition_set, &
        number_of_stopping_conditions_set, &
        MPI_INTEGER, &
        local_type_of_stopping_condition_set, &
        counts, &
        displs, &
        MPI_INTEGER, &
        0, &
        MPI_COMM_WORLD, &
        error&
    )
    
    if (sc_mpi_rank .EQ. 0) then
    
        offset = 0
        do i = 1, sc_mpi_size
            displs(i) = offset;
            counts(i) = counts(i) * (MAX_NUMBER_OF_PARTICLES_PER_INDEX+1);
            offset = offset + counts(i)
        end do
        
    end if


    call MPI_Gatherv( &
        index_of_particle_in_stopping_condition,&
        number_of_stopping_conditions_set*(MAX_NUMBER_OF_PARTICLES_PER_INDEX+1),&
        MPI_INTEGER,&
        local_index_of_particle_in_stopping_condition,&
        counts, &
        displs,&
        MPI_INTEGER,&
        0, &
        MPI_COMM_WORLD&
    );
    if (sc_mpi_rank .EQ. 0) then
    
        number_of_stopping_conditions_set = local_number_of_stopping_conditions_set;
        index_of_particle_in_stopping_condition = local_index_of_particle_in_stopping_condition
        type_of_stopping_condition_set = local_type_of_stopping_condition_set
    end if
    
end function

#else


function mpi_setup_stopping_conditions()
    integer :: mpi_setup_stopping_conditions
    
    sc_mpi_size = 1
    mpi_setup_stopping_conditions = 0
end function

function mpi_collect_stopping_conditions()
    integer :: mpi_collect_stopping_conditions
    mpi_collect_stopping_conditions = 0
end function

function mpi_distribute_stopping_conditions()
    integer :: mpi_distribute_stopping_conditions
    mpi_distribute_stopping_conditions = 0
    
    if(enabled_conditions .EQ. 0) then
            return
    end if
    
end function

#endif

end module
