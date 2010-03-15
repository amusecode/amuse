      module amuse_support
         implicit none
         character (len=256) :: AMUSE_inlist_path
         character (len=256) :: AMUSE_mesa_data_dir
         character (len=256) :: AMUSE_zams_filename = 'zams_z2m2'
         ! (use the solar metallicity model from the MESA starting_models folder)
         double precision :: AMUSE_metallicity = 0.02d0
         double precision :: AMUSE_dmass = 0.1d0
         double precision :: AMUSE_mlo = -1.0d0
         double precision :: AMUSE_mhi = 1.7d0
         double precision :: AMUSE_max_age_stop_condition = 1.0d12
         double precision :: AMUSE_min_timestep_stop_condition = 1.0d-6
         integer :: AMUSE_max_iter_stop_condition = -1111
         double precision :: AMUSE_mixing_length_ratio = 2.0d0
         double precision :: AMUSE_semi_convection_efficiency = 0.0d0
         logical :: debugging = .true.
         contains
         logical function failed(str, ierr)
            character (len=*), intent(in) :: str
            integer, intent(in) :: ierr
            failed = (ierr /= 0)
            if (failed) write(*, *) trim(str) // ' ierr', ierr
         end function failed
         subroutine get_zams_filename(str, ierr)
            character (len=256), intent(out) :: str
            integer, intent(out) :: ierr
            character (len=256) :: metallicity_str
            integer :: metallicity_exp, metallicity_factor
            if (AMUSE_metallicity.eq.0.0d0) then
               str = trim(AMUSE_mesa_data_dir) // '/star_data/starting_models/zams_z0m0'
            elseif (AMUSE_metallicity.eq.0.02d0) then
               str = trim(AMUSE_mesa_data_dir) // '/star_data/starting_models/zams_z2m2'
            else
               metallicity_exp = floor(log10(AMUSE_metallicity))-1
               metallicity_factor = floor(0.5 + AMUSE_metallicity/(1.0d1**metallicity_exp))
               write(metallicity_str,'(I0, A, I0)') metallicity_factor, "m", &
                  -metallicity_exp
               str = trim(AMUSE_mesa_data_dir) // '/star_data/starting_models/zams_z' &
                  // trim(metallicity_str)
            endif
            ierr = 0
         end subroutine get_zams_filename
      end module amuse_support

! Set the paths to the inlist and the data directory
      integer function set_MESA_paths(AMUSE_inlist_path_in, AMUSE_mesa_data_dir_in)
         use amuse_support, only: AMUSE_inlist_path, AMUSE_mesa_data_dir
         implicit none
         character(*), intent(in) :: AMUSE_inlist_path_in, AMUSE_mesa_data_dir_in
         AMUSE_inlist_path = AMUSE_inlist_path_in
         AMUSE_mesa_data_dir = AMUSE_mesa_data_dir_in
         set_MESA_paths = 0
      end function set_MESA_paths
    
! Initialize the stellar evolution code
      integer function initialize_code()
         use amuse_support, only: failed, AMUSE_mesa_data_dir, AMUSE_inlist_path
         use run_star_support
         use ctrls_io, only: set_default_controls
         implicit none
         integer :: initialize_code
         integer :: ierr
         initialize_code = -1
         call set_default_controls
         call do_read_star_job(AMUSE_inlist_path, ierr)
         if (failed('do_read_star_job', ierr)) return
         ! Replace value of mesa_data_dir just read, with supplied path.
         mesa_data_dir = AMUSE_mesa_data_dir
         call star_init(mesa_data_dir, kappa_file_prefix, &
            net_reaction_filename, net_rate_list_fname, ppn_rate_numbers_fname, ierr)
         if (failed('star_init', ierr)) return
         profile_columns_file = trim(mesa_data_dir) // '/star_data/profile_columns.list'
         log_columns_file = trim(mesa_data_dir) // '/star_data/log_columns.list'
         initialize_code = 0
      end function initialize_code

! Create new ZAMS model for a different metallicity
   subroutine new_zams_model(ierr)
      use create_zams, only: AMUSE_do_create_zams
      use amuse_support
      use star_lib, only: alloc_star, star_setup
      use star_private_def, only: star_info, get_star_ptr
      use run_star_support, only: run_create_zams, zams_inlist
      implicit none
      integer :: ierr, AMUSE_id
      type (star_info), pointer :: s
      AMUSE_id = alloc_star(ierr)
      if (failed('alloc_star', ierr)) return
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) return
      call star_setup(AMUSE_id, AMUSE_inlist_path, ierr)
      if (failed('star_setup', ierr)) return      
      call get_zams_filename(AMUSE_zams_filename, ierr)
      if (failed('get_zams_filename', ierr)) return
      run_create_zams = .true. ! is this necessary?
      call AMUSE_do_create_zams(s, AMUSE_metallicity, AMUSE_zams_filename, &
         AMUSE_dmass, AMUSE_mlo, AMUSE_mhi, ierr)
      if (failed('AMUSE_do_create_zams', ierr)) return      
      ierr = 0
   end subroutine new_zams_model

! Create a new particle
   function new_particle(AMUSE_id, AMUSE_mass)
      use amuse_support
      use star_lib, only: alloc_star, star_setup, star_load_zams, show_log_header
      use star_private_def, only: star_info, get_star_ptr
      use run_star_support, only: setup_for_run_star, before_evolve
      implicit none
      integer, intent(out) :: AMUSE_id
      integer :: new_particle, ierr
      double precision, intent(in) :: AMUSE_mass
      type (star_info), pointer :: s
      new_particle = -1
      AMUSE_id = alloc_star(ierr)
      if (failed('alloc_star', ierr)) return
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) return
      call star_setup(AMUSE_id, AMUSE_inlist_path, ierr)
      if (failed('star_setup', ierr)) return
      ! Replace value of mass and metallicity just read, with supplied values.
      s% initial_mass = AMUSE_mass
      s% initial_z = AMUSE_metallicity
      s% zams_filename = trim(AMUSE_zams_filename) // '.data'
      s% max_age = AMUSE_max_age_stop_condition
      s% min_timestep_limit = AMUSE_min_timestep_stop_condition
      s% max_model_number = AMUSE_max_iter_stop_condition
      s% mixing_length_alpha = AMUSE_mixing_length_ratio
      s% alpha_semiconvection = AMUSE_semi_convection_efficiency
      if (debugging) then
         write (*,*) "Creating new particles with mass: ", s% initial_mass
         write (*,*) "Loading starting model from: ", s% zams_filename
      endif
      call star_load_zams(AMUSE_id, ierr)
      if (failed('star_load_zams', ierr)) return
      call setup_for_run_star(AMUSE_id, s, .false., ierr)
      if (failed('setup_for_run_star', ierr)) return
      call before_evolve(AMUSE_id, ierr)
      if (failed('before_evolve', ierr)) return
      call show_log_header(AMUSE_id, ierr)
      if (failed('show_log_header', ierr)) return
      new_particle = 0
   end function

! Remove a particle (doesn't do anything yet)
   function delete_star(AMUSE_id)
      implicit none
      integer, intent(in) :: AMUSE_id
      integer :: delete_star
      delete_star = 0
   end function

! Does nothing...
   function initialize_stars()
      implicit none
      integer :: initialize_stars
      initialize_stars = 0
   end function

! Get/setters for code parameters:

! Return the number of particles currently allocated in the code
   function get_number_of_particles(AMUSE_value)
      use ctrls_io
      use run_star_support
      implicit none
      integer :: get_number_of_particles
      integer, intent(out) :: AMUSE_value
      AMUSE_value = 1
      get_number_of_particles = -1
   end function

! Return the metallicity parameter
   integer function get_metallicity(AMUSE_value)
      use amuse_support, only: AMUSE_metallicity
      implicit none
      double precision, intent(out) :: AMUSE_value
      AMUSE_value = AMUSE_metallicity
      get_metallicity = 0
   end function

! Set the metallicity parameter
   integer function set_metallicity(AMUSE_value)
      use amuse_support, only: AMUSE_metallicity, &
         AMUSE_zams_filename, failed, get_zams_filename
      use utils_lib, only: alloc_iounit, free_iounit
      implicit none
      double precision, intent(in) :: AMUSE_value
      integer :: ierr, iounit
      character (len=256) :: file
      set_metallicity = -1
      AMUSE_metallicity = AMUSE_value
      call get_zams_filename(AMUSE_zams_filename, ierr)
      if (failed('get_zams_filename', ierr)) return
      ! Check if the ZAMS model file exists
      iounit = alloc_iounit(ierr)
      if (failed('alloc_iounit', ierr)) return
      file = trim(AMUSE_zams_filename) // '.data'
      open(iounit, file=trim(file), action='read', status='old', iostat=ierr)
      if (ierr == 0) then
         close(iounit)
         call free_iounit(iounit)
      else
         call free_iounit(iounit)
         call new_zams_model(ierr) ! Have to create a new model otherwise.
         if (failed('new_zams_model', ierr)) return
      endif
      set_metallicity = 0
   end function

! Return the current mass of the star
   function get_mass(AMUSE_id, AMUSE_value)
      use star_private_def, only: star_info, get_star_ptr
      use amuse_support, only: failed
      implicit none
      integer, intent(in) :: AMUSE_id
      double precision, intent(out) :: AMUSE_value
      integer :: get_mass, ierr
      type (star_info), pointer :: s
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) then
         AMUSE_value = -1.0
         get_mass = -1
      else
         AMUSE_value = s% star_mass
         get_mass = 0
      endif
   end function

! Return the current temperature of the star
   function get_temperature(AMUSE_id, AMUSE_value)
      use star_private_def, only: star_info, get_star_ptr
      use amuse_support, only: failed
      implicit none
      integer, intent(in) :: AMUSE_id
      double precision, intent(out) :: AMUSE_value
      integer :: get_temperature, ierr
      type (star_info), pointer :: s
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) then
         AMUSE_value = -1.0
         get_temperature = -1
      else
         AMUSE_value = s% Teff
         get_temperature = 0
      endif
   end function

! Return the current luminosity of the star
      function get_luminosity(AMUSE_id, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id
         double precision, intent(out) :: AMUSE_value
         integer :: get_luminosity, ierr
         type (star_info), pointer :: s
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            AMUSE_value = -1.0
            get_luminosity = -1
         else
            AMUSE_value = 10.0d0**s% log_surface_luminosity
            get_luminosity = 0
         endif
      end function

! Return the current age of the star
      function get_age(AMUSE_id, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id
         double precision, intent(out) :: AMUSE_value
         integer :: get_age, ierr
         type (star_info), pointer :: s
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            AMUSE_value = -1.0
            get_age = -1
         else
            AMUSE_value = s% star_age
            get_age = 0
         endif
      end function

! Return the next timestep for the star
      function get_time_step(AMUSE_id, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use star_lib, only: star_pick_next_timestep
         use const_def, only: secyer
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id
         double precision, intent(out) :: AMUSE_value
         integer :: get_time_step, ierr
         type (star_info), pointer :: s
         AMUSE_value = -1.0
         get_time_step = -1
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) return
!         write (*,*) s% star_age, s% dt_next/secyer, s% dt/secyer, s% time_step
         if (s% star_age > 0.001) then
            ierr = star_pick_next_timestep(AMUSE_id)
            if (failed('star_pick_next_timestep', ierr)) return
         endif
!         write (*,*) s% star_age, s% dt_next/secyer, s% dt/secyer, s% time_step
         AMUSE_value = s% dt_next/secyer
         get_time_step = 0
      end function

! Return the current radius of the star
      function get_radius(AMUSE_id, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id
         double precision, intent(out) :: AMUSE_value
         integer :: get_radius, ierr
         type (star_info), pointer :: s
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            AMUSE_value = -1.0
            get_radius = -1
         else
            AMUSE_value = 10.0d0**s% log_surface_radius
            get_radius = 0
         endif
      end function

! Return the current stellar type of the star
      function get_stellar_type(AMUSE_id, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         use do_one_utils, only: do_show_log_header, do_write_log_state
         use utils, only:eval_current_y, eval_current_z
         implicit none
         integer, intent(in) :: AMUSE_id
         integer, intent(out) :: AMUSE_value
         integer :: get_stellar_type, ierr
         double precision :: x_avg, y_avg, z_avg
         type (star_info), pointer :: s
         AMUSE_value = -99
         get_stellar_type = -1
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) return
         select case(s% phase_of_evolution)
            case(0,1,2)
               if (s% initial_mass < 0.75) then
                  AMUSE_value = 0
               else
                  AMUSE_value = 1
               endif
            case(3)
               AMUSE_value = 3
            case(4:)
               y_avg = eval_current_y(s, 1, s% nz, ierr)
               if (failed('eval_current_y', ierr)) return
               z_avg = eval_current_z(s, 1, s% nz, ierr)
               if (failed('eval_current_z', ierr)) return
               x_avg = max(0d0, min(1d0, 1 - (y_avg + z_avg)))
               if (x_avg > 1.0d-5) then
                  if (s% center_he3 + s% center_he4 > 1.0d-5) then
                     AMUSE_value = 4 ! Core He burning
                  else
                     if (y_avg < 0.75 * x_avg) then
                        AMUSE_value = 5 ! Early AGB (inert C/O core)
                     else
                        AMUSE_value = 6 ! Late (thermally pulsing) AGB (inert C/O core)
                     endif
                  endif
               else
                  if (s% center_he3 + s% center_he4 > 1.0d-5) then
                     AMUSE_value = 7 ! Helium MS star
                  else
                     AMUSE_value = 9 ! Helium giant
                  endif
               endif
            case default
               write(*,*) "Unable to determine the stellar type."
               write(*,*) "The following information might help:"
               call do_show_log_header
               call do_write_log_state(s)
               return
         end select
         get_stellar_type = 0
!      integer, parameter :: phase_starting = 0
!      integer, parameter :: phase_early_main_seq = 1
!      integer, parameter :: phase_mid_main_seq = 2
!      integer, parameter :: phase_wait_for_he = 3
!      integer, parameter :: phase_he_ignition_over = 4
!      integer, parameter :: phase_he_igniting = 5
!      integer, parameter :: phase_helium_burning = 6
!      integer, parameter :: phase_carbon_burning = 7
      end function

! Evolve the star for one step
   function evolve(AMUSE_id)
      use star_lib, only: star_finish_step
      use star_private_def, only: star_info, get_star_ptr
      use run_star_support
      use run_star, only: check_model
      use amuse_support, only: failed
      implicit none
      integer, intent(in) :: AMUSE_id
      integer :: evolve
      type (star_info), pointer :: s
      integer :: ierr, model_number, result, result_reason
      logical :: first_try
      evolve = -1
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) return
      if (auto_extend_net) then
         call extend_net(s, ierr)
         if (failed('extend_net', ierr)) return
      end if
      first_try = .true.
      model_number = get_model_number(AMUSE_id, ierr)
      if (failed('get_model_number', ierr)) return
      step_loop: do ! may need to repeat this loop for retry or backup
         result = star_evolve_step(AMUSE_id, first_try)
         if (result == keep_going) result = check_model(s, AMUSE_id, 0)
         if (result == keep_going) result = star_pick_next_timestep(AMUSE_id)
         if (result == keep_going) exit step_loop
         model_number = get_model_number(AMUSE_id, ierr)
         if (failed('get_model_number', ierr)) return
         result_reason = get_result_reason(AMUSE_id, ierr)
         if (result == retry) then
            if (failed('get_result_reason', ierr)) return
            if (report_retries) &
               write(*,'(i6,3x,a,/)') model_number, &
                  'retry reason ' // trim(result_reason_str(result_reason))
         else if (result == backup) then
            if (failed('get_result_reason', ierr)) return
            if (report_backups) &
               write(*,'(i6,3x,a,/)') model_number, &
                  'backup reason ' // trim(result_reason_str(result_reason))
         end if
         if (result == retry) result = star_prepare_for_retry(AMUSE_id)
         if (result == backup) result = star_do1_backup(AMUSE_id)
         if (result == terminate) then
            if (result_reason == result_reason_normal) then
               if (s% max_model_number > 0 .and. s% model_number >= &
                  s% max_model_number) evolve = -3 ! max iterations reached
               if (s% star_age >= s% max_age) evolve = -2 ! max_age reached
            end if
            return
         end if
         first_try = .false.
      end do step_loop
      if (result == keep_going) then
         result = star_finish_step(AMUSE_id, .false.)
         if (result /= keep_going) return
      else if (result == terminate) then
         if (result_reason == result_reason_normal) then
            result = star_finish_step(AMUSE_id, save_photo_when_terminate)
            if (result /= keep_going) return
         end if
         evolve = 0
         return
      end if
      if (s% model_number == save_model_number) then
         call star_write_model(AMUSE_id, save_model_filename, .true., ierr)
         if (failed('star_write_model', ierr)) return
         write(*, *) 'saved to ' // trim(save_model_filename)
      end if
      if (s% model_number == profile_model_number) then
         write(*, '(a, i7)') 'save profile for model number', s% model_number
         call save_profile(AMUSE_id, 3, ierr)
         if (failed('save_profile', ierr)) return
      end if
      if (internals_num >= 0) then
         write(*, '(a, i7)') 'write internals for model number', s% model_number
         call std_write_internals(AMUSE_id, internals_num)
         stop 'finished std_write_internals'
      end if
      evolve = 0
   end function

! Evolve the star until AMUSE_end_time
      subroutine evolve_to(AMUSE_id, AMUSE_end_time)
         use star_lib, only: star_finish_step
         use star_private_def, only: star_info, get_star_ptr
         use run_star_support
         use run_star, only: check_model
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id
         double precision, intent(in) :: AMUSE_end_time
         type (star_info), pointer :: s
         integer :: ierr, model_number, result, result_reason
         double precision :: old_max_age
         logical :: first_try, continue_evolve_loop
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) return
         old_max_age = s% max_age
         s% max_age = AMUSE_end_time
         continue_evolve_loop = .true.
         evolve_loop: do while(continue_evolve_loop) ! evolve one step per loop
            if (auto_extend_net) then
               call extend_net(s, ierr)
               if (failed('extend_net', ierr)) return
            end if
            first_try = .true.
            model_number = get_model_number(AMUSE_id, ierr)
            if (failed('get_model_number', ierr)) return
            step_loop: do ! may need to repeat this loop for retry or backup
               result = star_evolve_step(AMUSE_id, first_try)
               if (result == keep_going) result = check_model(s, AMUSE_id, 0)
               if (result == keep_going) result = star_pick_next_timestep(AMUSE_id)
               if (result == keep_going) exit step_loop
               model_number = get_model_number(AMUSE_id, ierr)
               if (failed('get_model_number', ierr)) return
               result_reason = get_result_reason(AMUSE_id, ierr)
               if (result == retry) then
                  if (failed('get_result_reason', ierr)) return
                  if (report_retries) &
                     write(*,'(i6,3x,a,/)') model_number, &
                        'retry reason ' // trim(result_reason_str(result_reason))
               else if (result == backup) then
                  if (failed('get_result_reason', ierr)) return
                  if (report_backups) &
                     write(*,'(i6,3x,a,/)') model_number, &
                        'backup reason ' // trim(result_reason_str(result_reason))
               end if
               if (result == retry) result = star_prepare_for_retry(AMUSE_id)
               if (result == backup) result = star_do1_backup(AMUSE_id)
               if (result == terminate) then
                  if (result_reason == result_reason_normal) then
                     write(*, '(a, i12)') 'save profile for model number ', s% model_number
                     call save_profile(AMUSE_id, 3, ierr)
                     if (failed('save_profile', ierr)) return
                  end if
                  continue_evolve_loop = .false.
                  exit step_loop
               end if
               first_try = .false.
            end do step_loop
            if (result == keep_going) then
               result = star_finish_step(AMUSE_id, .false.)
               if (result /= keep_going) exit evolve_loop
            else if (result == terminate) then
               if (result_reason == result_reason_normal) then
                  result = star_finish_step(AMUSE_id, save_photo_when_terminate)
               end if
               exit evolve_loop
            end if
            if (s% model_number == save_model_number) then
               call star_write_model(AMUSE_id, save_model_filename, .true., ierr)
               if (failed('star_write_model', ierr)) return
               write(*, *) 'saved to ' // trim(save_model_filename)
            end if
            if (s% model_number == profile_model_number) then
               write(*, '(a, i7)') 'save profile for model number', s% model_number
               call save_profile(AMUSE_id, 3, ierr)
               if (failed('save_profile', ierr)) return
            end if
            if (internals_num >= 0) then
               write(*, '(a, i7)') 'write internals for model number', s% model_number
               call std_write_internals(AMUSE_id, internals_num)
               stop 'finished std_write_internals'
            end if
         end do evolve_loop
         s% max_age = old_max_age
         return
      end
      
! Return the maximum age stop condition
      integer function get_max_age_stop_condition(AMUSE_value)
         use amuse_support, only: AMUSE_max_age_stop_condition
         implicit none
         double precision, intent(out) :: AMUSE_value
         AMUSE_value = AMUSE_max_age_stop_condition
         get_max_age_stop_condition = 0
      end function get_max_age_stop_condition

! Set the maximum age stop condition
      integer function set_max_age_stop_condition(AMUSE_value)
         use amuse_support, only: AMUSE_max_age_stop_condition
         implicit none
         double precision, intent(in) :: AMUSE_value
         AMUSE_max_age_stop_condition = AMUSE_value
         set_max_age_stop_condition = 0
      end function set_max_age_stop_condition

! Return the maximum age stop condition
      integer function get_max_iter_stop_condition(AMUSE_value)
         use amuse_support, only: AMUSE_max_iter_stop_condition
         implicit none
         integer, intent(out) :: AMUSE_value
         AMUSE_value = AMUSE_max_iter_stop_condition
         get_max_iter_stop_condition = 0
      end function get_max_iter_stop_condition

! Set the maximum age stop condition
      integer function set_max_iter_stop_condition(AMUSE_value)
         use amuse_support, only: AMUSE_max_iter_stop_condition
         implicit none
         integer, intent(in) :: AMUSE_value
         AMUSE_max_iter_stop_condition = AMUSE_value
         set_max_iter_stop_condition = 0
      end function set_max_iter_stop_condition

! Return the minimum timestep stop condition
      integer function get_min_timestep_stop_condition(AMUSE_value)
         use amuse_support, only: AMUSE_min_timestep_stop_condition
         implicit none
         double precision, intent(out) :: AMUSE_value
         AMUSE_value = AMUSE_min_timestep_stop_condition
         get_min_timestep_stop_condition = 0
      end function get_min_timestep_stop_condition

! Set the minimum timestep stop condition
      integer function set_min_timestep_stop_condition(AMUSE_value)
         use amuse_support, only: AMUSE_min_timestep_stop_condition
         implicit none
         double precision, intent(in) :: AMUSE_value
         AMUSE_min_timestep_stop_condition = AMUSE_value
         set_min_timestep_stop_condition = 0
      end function set_min_timestep_stop_condition

! Retrieve the current value of the mixing length ratio
      integer function get_mixing_length_ratio(AMUSE_value)
         use amuse_support, only: AMUSE_mixing_length_ratio
         implicit none
         double precision, intent(out) :: AMUSE_value
         AMUSE_value = AMUSE_mixing_length_ratio
         get_mixing_length_ratio = 0
      end function get_mixing_length_ratio

! Set the current value of the mixing length ratio
      integer function set_mixing_length_ratio(AMUSE_value)
         use amuse_support, only: AMUSE_mixing_length_ratio
         implicit none
         double precision, intent(in) :: AMUSE_value
         AMUSE_mixing_length_ratio = AMUSE_value
         set_mixing_length_ratio = 0
      end function set_mixing_length_ratio

! Retrieve the current value of the semi convection efficiency
      integer function get_semi_convection_efficiency(AMUSE_value)
         use amuse_support, only: AMUSE_semi_convection_efficiency
         implicit none
         double precision, intent(out) :: AMUSE_value
         AMUSE_value = AMUSE_semi_convection_efficiency
         get_semi_convection_efficiency = 0
      end function get_semi_convection_efficiency

! Set the current value of the semi convection efficiency
      integer function set_semi_convection_efficiency(AMUSE_value)
         use amuse_support, only: AMUSE_semi_convection_efficiency
         implicit none
         double precision, intent(in) :: AMUSE_value
         AMUSE_semi_convection_efficiency = AMUSE_value
         set_semi_convection_efficiency = 0
      end function set_semi_convection_efficiency

! Retrieve the maximum number of stars that can be allocated in the code
      integer function get_maximum_number_of_stars(AMUSE_value)
         use star_def, only: max_star_handles
         implicit none
         integer, intent(out) :: AMUSE_value
         AMUSE_value = max_star_handles
         get_maximum_number_of_stars = 0
      end function get_maximum_number_of_stars

