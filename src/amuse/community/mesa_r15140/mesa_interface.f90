module mesa_interface
    use star_lib
    use star_def
    use const_def
    use const_lib
    use math_lib
    use chem_def
    use eos_lib
    use eos_def
    use kap_lib
    use run_star_support
    use run_star_extras, rse_extras_controls => extras_controls

    implicit none

    integer, parameter :: MESA_FAIL=-1, MESA_SUCESS=0
    logical :: use_gyre=.false.

    integer, parameter :: M_CENTER=0
    integer, parameter :: R_CENTER=1
    integer, parameter :: L_CENTER=2
    integer, parameter :: V_CENTER=3

    character(len=256),dimension(max_star_handles) :: restart_name = ''


    contains

! ***********************************************************************
! Routines for making a new model
! ***********************************************************************

    subroutine allocate_star(id, ierr)
        use ctrls_io, only: set_default_controls, store_controls
        use star_job_ctrls_io, only: set_default_star_job_controls, store_star_job_controls
        integer, intent(out) :: id
        integer, intent(out) :: ierr
        integer :: i
        type (star_info), pointer :: s

        call alloc_star(id, ierr)

        if (failed('alloc_star',ierr)) return

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        call set_default_star_job_controls()
        call store_star_job_controls(s, ierr)
        if (failed('store_star_job_controls',ierr)) return

        call set_default_controls()
        call store_controls(s, ierr)
        if (failed('store_controls',ierr)) return

        s% inlist_fname = ''

        call init_starting_star_data(s, ierr)

    end subroutine allocate_star

    subroutine load_inlist(id, inlist, ierr)
        use ctrls_io, only: read_controls_file
        use star_job_ctrls_io, only: read_star_job_file, check_star_job_controls
        integer, intent(in) :: id
        character (len=*), intent(in) :: inlist
        integer, intent(out) :: ierr
        type (star_info), pointer :: s

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        s% inlist_fname = inlist

        call read_star_job_file(s, s% inlist_fname, 1, ierr)
        if (failed('read_star_job_file',ierr)) return

        call check_star_job_controls(s, ierr)
        if (failed('check_star_job_controls',ierr)) return

        call read_controls_file(s, s% inlist_fname, 1, ierr)
        if (failed('read_controls_file',ierr)) return

    end subroutine load_inlist


    subroutine create_pre_main_sequence(id, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr

        type (star_info), pointer :: s

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        s% job% create_pre_main_sequence_model = .true.
        s% job% load_saved_model = .false.

    end subroutine create_pre_main_sequence


    subroutine load_saved_model(id, filename, ierr)
        integer, intent(in) :: id
        character (len=*), intent(in) :: filename
        integer, intent(out) :: ierr

        type (star_info), pointer :: s

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        s% job% load_saved_model = .true.
        s% job% saved_model_name = trim(filename)

        s% job% create_pre_main_sequence_model = .false.
        
    end subroutine load_saved_model


    subroutine load_zams_model(id, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr

        type (star_info), pointer :: s

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        s% job% load_saved_model = .false.
        s% job% create_pre_main_sequence_model = .false.
        
    end subroutine load_zams_model


    subroutine load_mesa_photo(id, filename, ierr)
        integer, intent(in) :: id
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ierr
        type (star_info), pointer :: s

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        s% photo_directory = '.'
        restart_name(id) = trim(filename)


    end subroutine load_mesa_photo

    subroutine save_mesa_photo(id, filename, ierr)
        integer, intent(in) :: id
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ierr

        call star_save_for_restart(id, filename, ierr)

    end subroutine save_mesa_photo


    subroutine save_mesa_model(id, filename, ierr)
        integer, intent(in) :: id
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ierr

        call star_write_model(id, filename, ierr)
        
    end subroutine save_mesa_model




    subroutine create_he_star(id, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr

        type (star_info), pointer :: s

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        s% job% create_pre_main_sequence_model = .true.

        s% initial_he3 = 0d0

        s% job% relax_initial_y = .true.
        s% job% new_Y = 1d0 - s% initial_z

        s% job% relax_initial_z = .true.
        s% job% new_Z = s% initial_z

    end subroutine create_he_star


    subroutine finish_init_star(id, init_callback, ierr)
        integer, intent(inout) :: id
        interface
            subroutine init_callback(x)
                integer, intent(in) :: x
            end subroutine init_callback
        end interface
        integer, intent(out) :: ierr
        logical :: restart
        type (star_info), pointer :: s
        logical, parameter :: dbg=.true.

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        call init_callback(id) ! Call here and in evolve_controls as there are options that 
                               ! need to be set before and after the other inlist options get set
        id_from_read_star_job = id

        call before_evolve_loop(.true., .true., restart, &
            null_binary_controls, extras_controls, &
            id_from_read_star_job, s% inlist_fname, restart_name(id), &
            dbg, 0, id, ierr)

        if (failed('before_evolve_loop',ierr)) return



        contains

        subroutine extras_controls(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            
           ! This is a local copy of extras_controls,
           ! Plese add your changes to the version in run_star_extras.f90 not this one

            s% doing_timing = .true.
            s% job% check_before_step_timing = 0
            s% job% check_step_loop_timing = 0
            s% job% check_after_step_timing = 0
            s% job% time0_initial = 0
            s% calculate_Brunt_N2 = .true.
    
            call init_callback(id)

            ! Set zbase if its not been set yet
            if(s%kap_rq% zbase == -1) s%kap_rq% zbase = s%initial_z            

            ! Override photo_directory otherwise we look in photos/ instead of ./
            ! This does nothing for 15140 but later versions can use this to set 
            ! the folder for the restart photo to the current directory.
            if(len_trim(restart_name(id)) > 0) then
                s%photo_directory = '.'
            end if

            call rse_extras_controls(id, ierr)
  
        
         end subroutine extras_controls

        
    end subroutine finish_init_star

    subroutine remove_star(id, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr 

        call free_star(id, ierr)
        if (failed('free_star',ierr)) return

    end subroutine remove_star

    subroutine max_num_stars(num)
        integer, intent(out) :: num

        num =  max_star_handles
 
    end subroutine max_num_stars


    integer function how_many_allocated_star_ids()
        integer :: id
        how_many_allocated_star_ids = 0
        if (have_initialized_star_handles) then
            do id = 1, max_star_handles
                if (star_handles(id)% in_use .eqv. .true.) &
                    how_many_allocated_star_ids = how_many_allocated_star_ids+1
            end do
        end if
    end function how_many_allocated_star_ids


! ***********************************************************************
! Routines for setting stellar parameters
! ***********************************************************************

    subroutine  set_init_options(id, &
                                mass, &
                                mesa_dir_in, &
                                mesa_data_dir_in, &
                                output_folder, &
                                inlist, &
                                metallicity,&
                                max_age_stop_condition, &
                                min_timestep_stop_condition, &
                                max_iter_stop_condition, &
                                gyre_in, &
                                temp_dir &
                                )
        integer, intent(in) :: id
        character(len=*), intent(in) :: mesa_dir_in, mesa_data_dir_in, output_folder, inlist, gyre_in, temp_dir
        real(dp), intent(in) :: mass, metallicity,&
                                max_age_stop_condition, &
                                min_timestep_stop_condition
        integer,intent(in) :: max_iter_stop_condition

        type (star_info), pointer :: s
        integer :: ierr

        ! Needs work as it breaks things
        !call chdir(output_folder)

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return   

        mesa_dir = mesa_dir_in

        call const_init(mesa_dir_in, ierr)
        if (failed('const_init',ierr)) return 

        mesa_data_dir = mesa_data_dir_in

        s% inlist_fname = inlist

        s% job% mesa_dir = mesa_dir_in
        s% initial_mass = mass
        s% initial_z = metallicity
        s% max_age = max_age_stop_condition
        s% min_timestep_limit = min_timestep_stop_condition
        s% max_model_number = max_iter_stop_condition

        ! Disable for now
        s% do_history_file = .false.
        s% write_profiles_flag = .false.
        s% history_interval = -1
        s% profile_interval = -1
        s% photo_interval = -1
        s% terminal_interval = 1
        s% write_header_frequency = 100

        mesa_temp_caches_dir = trim(temp_dir) 

        s% report_ierr = .false.

        if(len_trim(gyre_in) > 0) then
            call setup_gyre(gyre_in)
            use_gyre = .true.
        end if

    end subroutine set_init_options

    subroutine update_star_job(id, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        type (star_info), pointer :: s 

        ierr = MESA_SUCESS

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        call do_star_job_controls_after(id, s, .false., ierr)


    end subroutine update_star_job



    subroutine get_net_name(id, net_name, ierr)
        integer, intent(in) :: id
        character(len=*), intent(out) :: net_name
        integer, intent(out) :: ierr
        type (star_info), pointer :: s 

        ierr = MESA_SUCESS

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        net_name = s%net_name

    end subroutine get_net_name


    subroutine set_net_name(id, net_name, ierr)
        integer, intent(in) :: id
        character(len=*), intent(in) :: net_name
        integer, intent(out) :: ierr
        type (star_info), pointer :: s 

        ierr = MESA_SUCESS

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        call star_change_to_new_net(id, .true., net_name, ierr)

        if (failed('change_net',ierr)) return

    end subroutine set_net_name


    subroutine set_new_age(id, age, ierr)
        integer, intent(in) :: id
        real(dp) :: age ! in years
        integer, intent(out) :: ierr

        ierr = MESA_SUCESS

        call star_set_age(id, age, ierr)


    end subroutine set_new_age


    subroutine get_star_control_nml(id, name, val, ierr)
        integer, intent(in) :: id
        character(len=*) :: name, val
        integer, intent(out) :: ierr
        ierr = 0
        call star_get_control_namelist(id, name, val, ierr)
    end subroutine get_star_control_nml

    subroutine set_star_control_nml(id, name, val, ierr)
        integer, intent(in) :: id
        character(len=*) :: name, val
        integer, intent(out) :: ierr
        ierr = 0
        call star_set_control_namelist(id, name, val, ierr)
    end subroutine set_star_control_nml


    subroutine get_star_job_nml(id, name, val, ierr)
        integer, intent(in) :: id
        character(len=*) :: name, val
        integer, intent(out) :: ierr
        ierr = 0
        call star_get_star_job_namelist(id, name, val, ierr)
    end subroutine get_star_job_nml


    subroutine set_star_job_nml(id, name, val, ierr)
        integer, intent(in) :: id
        character(len=*) :: name, val
        integer, intent(out) :: ierr
        ierr = 0
        call star_set_star_job_namelist(id, name, val, ierr)
    end subroutine set_star_job_nml


    subroutine get_eos_nml(id, name, val, ierr)
        integer, intent(in) :: id
        character(len=*) :: name, val
        integer, intent(out) :: ierr
        ierr = 0
        call eos_get_control_namelist(id, name, val, ierr)
    end subroutine get_eos_nml


    subroutine set_eos_nml(id, name, val, ierr)
        integer, intent(in) :: id
        character(len=*) :: name, val
        integer, intent(out) :: ierr
        ierr = 0
        call eos_set_control_namelist(id, name, val, ierr)
    end subroutine set_eos_nml


    subroutine get_kap_nml(id, name, val, ierr)
        integer, intent(in) :: id
        character(len=*) :: name, val
        integer, intent(out) :: ierr
        ierr = 0
        call kap_get_control_namelist(id, name, val, ierr)
    end subroutine get_kap_nml

    subroutine set_kap_nml(id, name, val, ierr)
        integer, intent(in) :: id
        character(len=*) :: name, val
        integer, intent(out) :: ierr
        ierr = 0
        call kap_set_control_namelist(id, name, val, ierr)
    end subroutine set_kap_nml


! ***********************************************************************
! Routines for evovling a star
! ***********************************************************************

    subroutine do_evolve_one_step(id, result, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: result, ierr
        type (star_info), pointer :: s 

        integer :: model_number
        logical :: continue_evolve_loop, first_try
        logical,parameter :: dbg=.true.

        ierr = MESA_SUCESS

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        call before_step_loop(s% id, ierr)
        if (failed('before_step_loop',ierr)) return

        result = s% extras_start_step(id)  
        if (result /= keep_going) return      

        first_try = .true.
        
        step_loop: do ! may need to repeat this loop
        
           if (stop_is_requested(s)) then
              continue_evolve_loop = .false.
              result = terminate
              exit
           end if
        
           result = star_evolve_step(id, first_try)
           if (result == keep_going) result = star_check_model(id)
           if (result == keep_going) result = s% extras_check_model(id)
           if (result == keep_going) result = star_pick_next_timestep(id)            
           if (result == keep_going) exit step_loop
           
           model_number = get_model_number(id, ierr)
           if (failed('get_model_number',ierr)) return
                          
           if (result == retry .and. s% job% report_retries) then
              write(*,'(i6,3x,a,/)') model_number, &
                 'retry reason ' // trim(result_reason_str(s% result_reason))
           end if
           
           if (result == redo) then
              result = star_prepare_to_redo(id)
           end if
           if (result == retry) then
              result = star_prepare_to_retry(id)
           end if
           if (result == terminate) then
              continue_evolve_loop = .false.
              exit step_loop
           end if
           first_try = .false.
        end do step_loop
        s% doing_first_model_of_run = .false.

        call after_step_loop(s% id, s% inlist_fname, &
            dbg, result, ierr)
        if (failed('after_step_loop',ierr)) return
        
        call do_saves(id, ierr)
        if (failed('do_saves',ierr)) return

        call flush()

    end subroutine do_evolve_one_step


    subroutine evolve_until(id, delta_t, ierr)
        integer, intent(in) :: id
        real(dp), intent(in) :: delta_t
        integer, intent(out) :: ierr
        type (star_info), pointer :: s  
        integer :: result
        logical,parameter :: dbg=.false.
        real(dp) :: old_age

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        old_age = s% max_age 

        s% max_age = s% star_age + delta_t
        s% num_adjusted_dt_steps_before_max_age =  -1

        evolve_loop: do 
            call do_evolve_one_step(id, result, ierr)
            if (result /= keep_going) then
                if (s% result_reason == result_reason_normal) then
                    exit evolve_loop
                else
                    ierr = -1
                    exit evolve_loop
                end if
            end if
            s% doing_first_model_of_run = .false.
        end do evolve_loop

        ! In case you want to evolve further
        s% max_age = old_age
        if(s% dt_next==0) s% dt_next = s% dt

        s% doing_first_model_of_run = .false.

    end subroutine evolve_until




! ***********************************************************************
! Routines for modifying a star
! ***********************************************************************

    subroutine set_initial_mass(id, new_mass, ierr)
        integer, intent(in) :: id
        real(dp), intent(in) :: new_mass
        integer, intent(out) :: ierr
        type (star_info), pointer :: s

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        s% initial_mass = new_mass

    end subroutine set_initial_mass

    subroutine set_model_number(id, new_number, ierr)
        integer, intent(in) :: id
        integer, intent(in) :: new_number
        integer, intent(out) :: ierr
        type (star_info), pointer :: s

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        s% model_number = new_number

    end subroutine set_model_number
    subroutine set_timestep(id, dt_next, ierr)
        integer, intent(in) :: id
        real(dp), intent(in) :: dt_next
        integer, intent(out) :: ierr

        call set_dt_next(id, dt_next, ierr)
        if (failed('set_dt_next',ierr)) return

    end subroutine set_timestep

    subroutine set_new_mass(id, new_mass, ierr)
        integer, intent(in) :: id
        real(dp), intent(in) :: new_mass
        integer, intent(out) :: ierr
        type (star_info), pointer :: s

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        call star_relax_mass(id, new_mass, s% job% lg_max_abs_mdot, ierr)
        if (failed('star_relax_mass',ierr)) return

    end subroutine set_new_mass


    subroutine set_new_metallicity(id, new_z, ierr)
        integer, intent(in) :: id
        real(dp), intent(in) :: new_z
        integer, intent(out) :: ierr
        type (star_info), pointer :: s

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        call star_relax_z(id, new_z, s% relax_dlnZ, 0d0, 1d0, ierr)
        if (failed('star_relax_z',ierr)) return

    end subroutine set_new_metallicity


    subroutine relax_to_new_comp(id, xa, xqs, ierr)
        integer, intent(in) :: id
        real(dp), intent(in) :: xa(:,:), xqs(:)
        type (star_info), pointer :: s   
        integer, intent(out) :: ierr
        integer :: k
        real(dp) :: max_age
        ierr=0 
        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        if(s%species /= size(xa,dim=1)) then
            ierr = -50
            return 
        end if

        call star_relax_composition(id, s% job% num_steps_to_relax_composition, &
                                    size(xqs),size(xa,dim=1) ,xa, xqs, ierr )


    end subroutine relax_to_new_comp


    subroutine relax_to_new_entropy(id, xqs, temperature, rho, ierr)
        integer, intent(in) :: id
        real(dp), intent(in) :: xqs(:), temperature(:), rho(:)
        type (star_info), pointer :: s   
        integer, intent(out) :: ierr
        real(dp), allocatable, dimension(:) :: entropy
        real(dp), dimension(num_eos_basic_results) :: res,d_dlnd, d_dlnT, d_dabar, d_dzbar
        integer :: num_pts, k
        ierr=0
        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        s% doing_first_model_of_run = .false.
        num_pts = size(xqs)

        allocate(entropy(num_pts))

        do k = 1, num_pts
            ! get entropy
            call eosDT_get( &
                s% eos_handle, 1 - s% X(k) - s% Y(k), s% X(k), s% abar(k), s% zbar(k), &
                s% species, s% chem_id, s% net_iso, s% xa(:,k), &
                rho(k), log10(rho(k)), temperature(k), log10(temperature(k)), &
                res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)
            if (ierr /= 0) then
                write(*,*) "failed in eosDT_get"
                return
            end if
            entropy(k) = exp(res(i_lnS))
        end do

        s% min_timestep_limit = 1d-15

        call star_relax_entropy(id, s% job% max_steps_to_relax_entropy, num_pts, entropy, xqs, ierr)
        deallocate(entropy)

    end subroutine relax_to_new_entropy



    subroutine set_star_center_value(id, boundary, val, ierr)
        integer, intent(in) :: id,boundary
        real(dp),intent(in) :: val
        integer, intent(out) :: ierr
        type (star_info), pointer :: s

        ierr = MESA_SUCESS
        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        select case (boundary)
            case(M_CENTER)
                call star_relax_m_center(id,  val, s% job% dlgm_per_step, s% job% relax_M_center_dt,ierr)
            case(R_CENTER)
                call star_relax_r_center(id, val, s% job% dlgR_per_step, s% job% relax_R_center_dt,ierr)
            case(L_CENTER)
                call star_relax_l_center(id,  val, s% job% dlgL_per_step, s% job% relax_L_center_dt,ierr)
            case(V_CENTER)
                call star_relax_v_center(id,  val, s% job% dv_per_step, s% job% relax_v_center_dt,ierr)
            case default
                ierr = MESA_FAIL
        end select


    end subroutine set_star_center_value


    subroutine get_star_center_value(id, boundary, val, ierr)
        integer, intent(in) :: id, boundary
        real(dp),intent(out) :: val
        integer, intent(out) :: ierr
        type (star_info), pointer :: s

        ierr = MESA_SUCESS
        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        select case (boundary)
            case(M_CENTER)
                val = s% m_center
            case(R_CENTER)
                val = s% r_center
            case(L_CENTER)
                val = s% l_center
            case(V_CENTER)
                val = s% v_center
            case default
                ierr = MESA_FAIL
        end select


    end subroutine get_star_center_value


    subroutine change_species_one_zone(id, zone, species, value, ierr)
        integer, intent(in) :: id, species, zone
        real(dp),intent(in) :: value
        integer, intent(out) :: ierr
        type (star_info), pointer :: s
        integer :: net_id

        ierr = MESA_SUCESS
        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        call set_abundance_in_section(id, s% chem_id(species), value,zone, zone, ierr)
        if (failed('star_set_abundance',ierr)) return

    end subroutine change_species_one_zone




! ***********************************************************************
! Routines for acessing stellar properites
! ***********************************************************************

    subroutine get_history_value(id, name, val, ierr)
        integer, intent(in) :: id
        character(*),intent(in) :: name

        real(dp),intent(out) :: val
        integer, intent(out) :: ierr
        integer,dimension(1) :: specs
        type (star_info), pointer :: s

        ierr = MESA_FAIL

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return
        ! Check name exists
        call star_history_specs(s, 1, (/name/), specs, .false.)
        if(specs(1)==-1) then
            ierr = MESA_FAIL
            return
        end if

        val = star_get_history_output_by_id(id, name)
        if(val == -HUGE(val)) then
            ierr = MESA_FAIL
        else
            ierr = MESA_SUCESS
        end if

    end subroutine get_history_value


    subroutine get_profile_value_zone(id, name, k, val, ierr)
        integer, intent(in) :: id, k
        character(*),intent(in) :: name

        real(dp),intent(out) :: val
        integer, intent(out) :: ierr
        integer :: profile_id
        type (star_info), pointer :: s

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        if(k<=0 .or. k > s%nz)then
            ierr = MESA_FAIL
            return
        end if

        ! Check it exists
        call star_get_profile_id(s, name, profile_id)
        if(profile_id<0) then
            val = -HUGE(val)
            ierr = MESA_FAIL
            return
        end if

        val = star_get_profile_output_by_id(id, name, k)

        if(val == -HUGE(val)) then
            ierr = MESA_FAIL
        else
            ierr = MESA_SUCESS
        end if

    end subroutine get_profile_value_zone

    subroutine get_profile_values(id, name, val, ierr)
        integer, intent(in) :: id
        character(*),intent(in) :: name

        real(dp),dimension(:),intent(out) :: val
        integer, intent(out) :: ierr
        integer :: i, profile_id
        type(star_info), pointer :: s

        call star_ptr(id, s, ierr)
        if (ierr /= 0) then
            ierr = MESA_FAIL
            val = -1
            return
        end if

        if (size(val)/=s%nz) then
            ierr = MESA_FAIL
            val = -1
            return
        end if

        ! Check it exists
        call star_get_profile_id(s, name, profile_id)

        if(profile_id<0) then
            ierr = MESA_FAIL
            return
        end if

        ierr = MESA_SUCESS
        val = -1
        do i=1, s%nz
            val(i) = star_get_profile_output_by_id(id, name, i)

            if(val(i) == -HUGE(val)) then
                ierr = MESA_FAIL
                exit
            end if
        end do

    end subroutine get_profile_values

    subroutine get_next_timestep(id, dt_next, ierr)
        integer, intent(in) :: id
        real(dp), intent(out) :: dt_next
        integer, intent(out) :: ierr

        call get_dt_next(id, dt_next, ierr)

    end subroutine get_next_timestep

 
    integer function reverse_zone_id(id, zone, ierr)
        integer, intent(in) :: id, zone
        integer, intent(out) :: ierr
        type (star_info), pointer :: s

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        reverse_zone_id = s%nz - zone

        if(reverse_zone_id < 0 .or. reverse_zone_id > s%nz) then
            ierr = MESA_FAIL
            reverse_zone_id = -1
        end if

    end function reverse_zone_id


    subroutine get_species_name(id, net_id, name, ierr)
        use chem_def
        integer, intent(in) :: id, net_id
        character(*), intent(out) :: name
        integer, intent(out) :: ierr
        type (star_info), pointer :: s    

        name=''

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        if(net_id < 0 .or. net_id > s%species) then
            ierr = -2
            return
        end if

        name = chem_isos% name(s% chem_id(net_id))

    end subroutine get_species_name


    subroutine get_species_id(id, iso_name, net_id, ierr)
        use chem_lib, only: chem_get_iso_id
        integer, intent(in) :: id
        integer, intent(out) :: net_id
        character(len=*), intent(in) :: iso_name
        integer, intent(out) :: ierr
        type (star_info), pointer :: s
        integer :: k    

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        ierr = -2
        do k=1,s%species
            if(chem_isos% name(s% chem_id(k)) == trim(iso_name)) then
                ierr = 0
                net_id = k
                exit
            end if
        end do

    end subroutine get_species_id

    subroutine get_mass_number_species(id, net_id, mass, ierr)
        integer, intent(in) :: id
        integer, intent(in) :: net_id
        real(dp), intent(out) :: mass
        integer, intent(out) :: ierr
        type (star_info), pointer :: s   
        integer :: species_id 

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        if(net_id < 1 .or. net_id > s% species) then
            ierr =-2
            return
        end if

        mass = chem_isos% W(s% chem_id(net_id))

    end subroutine get_mass_number_species


    subroutine get_total_mass_species(id, net_id, mass, ierr)
        integer, intent(in) :: id
        integer, intent(in) :: net_id
        real(dp), intent(out) :: mass
        integer, intent(out) :: ierr
        type (star_info), pointer :: s   
        integer :: species_id 

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        if(net_id < 1 .or. net_id > s% species) then
            ierr =-2
            return
        end if

        mass = dot_product(s% xa(net_id,1:s% nz),s% dq(1:s% nz))/sum(s% dq(1:s% nz)) 
        mass = mass*s% xmstar/Msun

    end subroutine get_total_mass_species

    subroutine get_species_at_zone(id, net_id, zone, mass, ierr)
        integer, intent(in) :: id, zone
        integer, intent(in) :: net_id
        real(dp), intent(out) :: mass
        integer, intent(out) :: ierr
        type (star_info), pointer :: s   
        integer :: species_id 

        mass = -1

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        if(zone<0 .or. zone> s%nz) then
            ierr = -2
            return
        end if 

        if(net_id < 1 .or. net_id > s% species) then
            ierr =-3
            return
        end if

        mass  = s% xa(net_id, zone)

    end subroutine get_species_at_zone

    subroutine time_of_step(id, time, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: time
        integer, intent(out) :: ierr
        type (star_info), pointer :: s   

        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        time = s% job% after_step_timing 

    end subroutine time_of_step


    ! Routines for setting/getting MESA varaibles that are not available via history/profile output
    subroutine get_value(id, name, val, k, ierr)
        integer, intent(in) :: id, k
        character(len=*), intent(in) :: name
        real(dp), intent(out) :: val
        integer, intent(out) :: ierr
        type (star_info), pointer :: s   

        ierr = 0
        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        if(k<1 .or. k> s% nz) then
            ierr = -2
            return
        end if

        select case(name)
            case('extra_heat') 
                val = s% extra_heat(k)
            case default
                ierr = -3
        end select


    end subroutine get_value


    subroutine set_value(id, name, val, k, ierr)
        integer, intent(in) :: id, k
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: val
        integer, intent(out) :: ierr
        type (star_info), pointer :: s   

        ierr = 0
        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        if(k<1 .or. k> s% nz) then
            ierr = -2
            return
        end if

        select case(name)
            case('extra_heat') 
                s% extra_heat(k) = val
            case default
                ierr = -3
        end select


    end subroutine set_value



! ***********************************************************************
! Routines for GYRE
! ***********************************************************************

    subroutine setup_gyre(gyre_in)
        use gyre_lib
        use const_def
        character(len=*), intent(in) :: gyre_in

        ! Initialize GYRE

        call gyre_init(gyre_in)

        ! Set constants
    
        call gyre_set_constant('G_GRAVITY', standard_cgrav)
        call gyre_set_constant('C_LIGHT', clight)
        call gyre_set_constant('A_RADIATION', crad)
    
        call gyre_set_constant('M_SUN', Msun)
        call gyre_set_constant('R_SUN', Rsun)
        call gyre_set_constant('L_SUN', Lsun)
    
        call gyre_set_constant('GYRE_DIR', TRIM(mesa_dir)//'/gyre/gyre')


    end subroutine setup_gyre


    subroutine get_gyre_data(id, mode_l, &
         add_center_point, keep_surface_point, add_atmosphere, &
         fileout, ierr)
        use gyre_lib
        integer, intent(in) :: id, mode_l
        integer, intent(out) :: ierr
        logical, intent(in) :: add_center_point, keep_surface_point, add_atmosphere
        type (star_info), pointer :: s  
        real(dp), allocatable     :: global_data(:)
        real(dp), allocatable     :: point_data(:,:)
        integer                   :: ipar(1)
        real(dp)                  :: rpar(1)
        character(len=*),intent(in) :: fileout

        ierr = 0
        call star_ptr(id, s, ierr)
        if (failed('star_ptr',ierr)) return

        if(.not.use_gyre) then
            ierr = -99
            return
        end if

        call star_get_pulse_data(s%id, 'GYRE', &
            add_center_point, keep_surface_point, add_atmosphere, &
            global_data, point_data, ierr)

        if (ierr /= 0) return

        call gyre_set_model(global_data, point_data, 101)


        call gyre_get_modes(mode_l, process_mode_, ipar, rpar)

        ierr = ipar(1)

    contains

        subroutine process_mode_ (md, ipar, rpar, retcode)
    
         type(mode_t), intent(in) :: md
         integer, intent(inout)   :: ipar(:)
         real(dp), intent(inout)  :: rpar(:)
         integer, intent(out)     :: retcode
   
         integer               :: k, unit
         complex(dp)           :: cfreq
         real(dp)              :: freq, growth
         type(grid_t)          :: gr
          
         retcode= 0 
         ipar(1) = 0
         cfreq = md% freq('HZ')
         gr = md%grid()

         open(newunit=unit,file=fileout,action='write',status='old',position="append")

         write(unit,*) md%n_pg,md%n_p,md%n_g,md%n_k,REAL(cfreq),IMAG(cfreq)

         do k=1,md%n_k
            write(unit,'(I8,6(1X,E24.16))') k, gr%pt(k)%x, &
                                            md%xi_r(k), md%xi_h(k), md%dW_dx(k)
         end do
         flush(unit)
         close(unit)

        end subroutine process_mode_

    end subroutine get_gyre_data




end module mesa_interface