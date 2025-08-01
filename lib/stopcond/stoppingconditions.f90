module StoppingConditions

    INTEGER :: COLLISION_DETECTION
    INTEGER :: PAIR_DETECTION
    INTEGER :: ESCAPER_DETECTION
    INTEGER :: TIMEOUT_DETECTION
    INTEGER :: NUMBER_OF_STEPS_DETECTION
    INTEGER :: OUT_OF_BOX_DETECTION
    INTEGER :: DENSITY_LIMIT_DETECTION
    INTEGER :: INTERNAL_ENERGY_LIMIT_DETECTION
    INTEGER :: INTERACTION_OVER_DETECTION
    INTEGER :: SUPERNOVA_DETECTION

    PARAMETER (COLLISION_DETECTION=0)
    PARAMETER (PAIR_DETECTION=1)
    PARAMETER (ESCAPER_DETECTION=2)
    PARAMETER (TIMEOUT_DETECTION=3)
    PARAMETER (NUMBER_OF_STEPS_DETECTION=4)
    PARAMETER (OUT_OF_BOX_DETECTION=5)
    PARAMETER (DENSITY_LIMIT_DETECTION=6)
    PARAMETER (INTERNAL_ENERGY_LIMIT_DETECTION=7)
    PARAMETER (INTERACTION_OVER_DETECTION=8)
    PARAMETER (SUPERNOVA_DETECTION=9)

    interface
! int enable_stopping_condition(int type);
        integer(c_int) function enable_stopping_condition(t) bind(c)
            use iso_c_binding
            integer (c_int), value :: t
        end function
! int set_support_for_condition(int type);
        integer(c_int) function set_support_for_condition(t) bind(c)
            use iso_c_binding
            integer (c_int), value :: t
        end function
! int get_number_of_stopping_conditions_set(int * result);
        integer(c_int) function get_number_of_stopping_conditions_set(r) bind(c)
            use iso_c_binding
            integer (c_int) :: r
        end function
! int is_stopping_condition_set(int type, int * result);
        integer(c_int) function is_stopping_condition_set(t,r) bind(c)
            use iso_c_binding
            integer (c_int), value :: t
            integer (c_int) :: r
        end function
! int is_stopping_condition_enabled(int type, int * result);
        integer(c_int) function is_stopping_condition_enabled(t,r) bind(c)
            use iso_c_binding
            integer (c_int), value :: t
            integer (c_int) :: r
        end function
! int disable_stopping_condition(int type);
        integer(c_int) function disable_stopping_condition(t) bind(c)
            use iso_c_binding
            integer (c_int), value :: t
        end function
! int has_stopping_condition(int type, int * result);
        integer(c_int) function has_stopping_condition(t,r) bind(c)
            use iso_c_binding
            integer (c_int), value :: t
            integer (c_int) :: r
        end function
! int get_stopping_condition_info(int index, int * index_of_the_condition, int *number_of_particles);
        integer(c_int) function get_stopping_condition_info(i, index_of_the_condition ,number_of_particles)  bind(c)
            use iso_c_binding
            integer (c_int), value :: i
            integer (c_int) ::  index_of_the_condition ,number_of_particles
        end function
! int get_stopping_condition_particle_index(int index, int index_in_the_condition, int * index_of_particle);
        integer(c_int) function get_stopping_condition_particle_index(i, index_of_the_condition , index_of_particle)  bind(c)
            use iso_c_binding
            integer (c_int), value :: i,index_of_the_condition 
            integer (c_int) :: index_of_particle
        end function
! int set_stopping_condition_timeout_parameter(double value);
        integer(c_int) function set_stopping_condition_timeout_parameter(v) bind(c)
            use iso_c_binding
            real (c_double), value :: v
        end function
! int get_stopping_condition_timeout_parameter(double * value);
        integer(c_int) function get_stopping_condition_timeout_parameter(v) bind(c)
            use iso_c_binding
            real (c_double) :: v
        end function
! int set_stopping_condition_number_of_steps_parameter(int value);
        integer(c_int) function set_stopping_condition_number_of_steps_parameter(v) bind(c)
            use iso_c_binding
            integer (c_int), value :: v
        end function
! int get_stopping_condition_number_of_steps_parameter(int *value);
        integer(c_int) function get_stopping_condition_number_of_steps_parameter(v) bind(c)
            use iso_c_binding
            integer (c_int) :: v
        end function
! int set_stopping_condition_out_of_box_parameter(double value);
        integer(c_int) function set_stopping_condition_out_of_box_parameter(v) bind(c)
            use iso_c_binding
            real (c_double), value :: v
        end function
! int get_stopping_condition_out_of_box_parameter(double *value);
        integer(c_int) function get_stopping_condition_out_of_box_parameter(v) bind(c)
            use iso_c_binding
            real (c_double) :: v
        end function
! int set_stopping_condition_minimum_density_parameter(double value);
        integer(c_int) function set_stopping_condition_minimum_density_parameter(v) bind(c)
            use iso_c_binding
            real (c_double), value :: v
        end function
! int get_stopping_condition_minimum_density_parameter(double *value);
        integer(c_int) function get_stopping_condition_minimum_density_parameter(v) bind(c)
            use iso_c_binding
            real (c_double) :: v
        end function
! int set_stopping_condition_minimum_internal_energy_parameter(double value);
        integer(c_int) function set_stopping_condition_minimum_internal_energy_parameter(v) bind(c)
            use iso_c_binding
            real (c_double), value :: v
        end function
! int get_stopping_condition_minimum_internal_energy_parameter(double *value);
        integer(c_int) function get_stopping_condition_minimum_internal_energy_parameter(v) bind(c)
            use iso_c_binding
            real (c_double) :: v
        end function
! int set_stopping_condition_maximum_density_parameter(double value);
        integer(c_int) function set_stopping_condition_maximum_density_parameter(v) bind(c)
            use iso_c_binding
            real (c_double), value :: v
        end function
! int get_stopping_condition_maximum_density_parameter(double *value);
        integer(c_int) function get_stopping_condition_maximum_density_parameter(v) bind(c)
            use iso_c_binding
            real (c_double) :: v
        end function
! int set_stopping_condition_maximum_internal_energy_parameter(double value);
        integer(c_int) function set_stopping_condition_maximum_internal_energy_parameter(v) bind(c)
            use iso_c_binding
            real (c_double), value :: v
        end function
! int get_stopping_condition_maximum_internal_energy_parameter(double *value);
        integer(c_int) function get_stopping_condition_maximum_internal_energy_parameter(v) bind(c)
            use iso_c_binding
            real (c_double) :: v
        end function
! int set_stopping_condition_out_of_box_use_center_of_mass_parameter(int value);
        integer(c_int) function set_stopping_condition_out_of_box_use_center_of_mass_parameter_(v) &
            bind(c, name='set_stopping_condition_out_of_box_use_center_of_mass_parameter')
            use iso_c_binding
            logical (c_bool), value :: v
        end function
! int get_stopping_condition_out_of_box_use_center_of_mass_parameter(int *value);
        integer(c_int) function get_stopping_condition_out_of_box_use_center_of_mass_parameter_(v) &
            bind(c, name='get_stopping_condition_out_of_box_use_center_of_mass_parameter')
            use iso_c_binding
            logical (c_bool) :: v
        end function

! int is_any_condition_set();
        integer(c_int) function is_any_condition_set() bind(c)
            use iso_c_binding
        end function

! int reset_stopping_conditions();
        integer(c_int) function reset_stopping_conditions() bind(c)
            use iso_c_binding
        end function

! int next_index_for_stopping_condition();
        integer(c_int) function next_index_for_stopping_condition() bind(c)
            use iso_c_binding
        end function

! int set_stopping_condition_info(int index, int type);
        integer(c_int) function set_stopping_condition_info(i,t) bind(c)
            use iso_c_binding
            integer (c_int), value :: i
            integer (c_int), value :: t
        end function

! int set_stopping_condition_particle_index(int index, int index_in_the_condition, int index_of_particle);
        integer(c_int) function set_stopping_condition_particle_index(i,index_in_the_condition, index_of_particle) bind(c)
            use iso_c_binding
            integer (c_int), value :: i,index_in_the_condition, index_of_particle
        end function
! int mpi_setup_stopping_conditions();
       integer(c_int) function mpi_setup_stopping_conditions() bind(c)
            use iso_c_binding
        end function
! int mpi_distribute_stopping_conditions();
       integer(c_int) function mpi_distribute_stopping_conditions() bind(c)
            use iso_c_binding
        end function
! int mpi_collect_stopping_conditions();
       integer(c_int) function mpi_collect_stopping_conditions() bind(c)
            use iso_c_binding
        end function

    end interface

contains 

    function set_stopping_condition_out_of_box_use_center_of_mass_parameter(v) result(ret)
        use iso_c_binding
        logical :: v
        logical (c_bool) :: cv
        integer :: ret
        cv=v
        ret=set_stopping_condition_out_of_box_use_center_of_mass_parameter_(cv)
    end function
    function get_stopping_condition_out_of_box_use_center_of_mass_parameter(v) result(ret)
        use iso_c_binding
        logical :: v
        logical(c_bool) :: cv
        ret=get_stopping_condition_out_of_box_use_center_of_mass_parameter_(cv)
        v=cv
    end function
    
    
    
end module
