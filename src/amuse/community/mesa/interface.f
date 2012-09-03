      module amuse_support
         implicit none
         character (len=4096) :: AMUSE_inlist_path
         character (len=4096) :: AMUSE_mesa_data_dir
         character (len=4096) :: AMUSE_local_data_dir ! Used for output starting_models
         character (len=4096) :: AMUSE_zams_filename = 'zams_z2m2'
         ! (use the solar metallicity model from the MESA starting_models folder)
         double precision :: AMUSE_metallicity = 0.02d0
         double precision :: AMUSE_dmass = 0.1d0
         double precision :: AMUSE_mlo = -1.0d0
         double precision :: AMUSE_mhi = 2.0d0
         double precision :: AMUSE_max_age_stop_condition = 1.0d12
         double precision :: AMUSE_min_timestep_stop_condition = 1.0d-6
         integer :: AMUSE_max_iter_stop_condition = -1111
         double precision :: AMUSE_mixing_length_ratio = 2.0d0
         double precision :: AMUSE_semi_convection_efficiency = 0.0d0
         integer :: AMUSE_RGB_wind_scheme = 1
         integer :: AMUSE_AGB_wind_scheme = 1
         double precision :: AMUSE_reimers_wind_efficiency = 0.5d0
         double precision :: AMUSE_blocker_wind_efficiency = 0.1d0
         double precision :: AMUSE_de_jager_wind_efficiency = 0.8d0
         double precision :: AMUSE_dutch_wind_efficiency = 0.8d0
         
         double precision, allocatable :: target_times(:)
         integer :: number_of_particles ! Dead or alive...
         
         logical :: new_model_defined = .false.
         integer :: id_new_model
         
         logical :: debugging = .false.
         logical :: do_stabilize_new_stellar_model = .true.
         
         contains
         logical function failed(str, ierr)
            character (len=*), intent(in) :: str
            integer, intent(in) :: ierr
            failed = (ierr /= 0)
            if (failed) write(*, *) trim(str) // ' ierr', ierr
         end function failed
         logical function evolve_failed(str, ierr, return_var, errorcode)
            character (len=*), intent(in) :: str
            integer, intent(in) :: ierr, errorcode
            integer, intent(out) :: return_var
            evolve_failed = (ierr /= 0)
            if (evolve_failed) then
               write(*, *) trim(str) // ' ierr', ierr
               return_var = errorcode
            endif
         end function evolve_failed
         subroutine get_zams_filename(str, ierr)
            character (len=1024), intent(out) :: str
            integer, intent(out) :: ierr
            character (len=1024) :: metallicity_str
            integer :: metallicity_exp, metallicity_factor
            if (AMUSE_metallicity.eq.0.0d0) then
               str = trim(AMUSE_local_data_dir) // '/star_data/starting_models/zams_z0m0'
            elseif (AMUSE_metallicity.eq.0.02d0) then
               str = trim(AMUSE_mesa_data_dir) // '/star_data/starting_models/zams_z2m2'
            else
               metallicity_exp = floor(log10(AMUSE_metallicity))-1
               metallicity_factor = floor(0.5 + AMUSE_metallicity/(1.0d1**metallicity_exp))
               write(metallicity_str,'(I0, A, I0)') metallicity_factor, "m", &
                  -metallicity_exp
               str = trim(AMUSE_local_data_dir) // '/star_data/starting_models/zams_z' &
                  // trim(metallicity_str)
            endif
            ierr = 0
         end subroutine get_zams_filename
      end module amuse_support

! Set the paths to the inlist and the data directory
      integer function set_MESA_paths(AMUSE_inlist_path_in, &
            AMUSE_mesa_data_dir_in, AMUSE_local_data_dir_in)
         use amuse_support, only: AMUSE_inlist_path, &
            AMUSE_mesa_data_dir, AMUSE_local_data_dir
         implicit none
         character(*), intent(in) :: AMUSE_inlist_path_in, &
            AMUSE_mesa_data_dir_in, AMUSE_local_data_dir_in
         AMUSE_inlist_path = AMUSE_inlist_path_in
         AMUSE_mesa_data_dir = AMUSE_mesa_data_dir_in
         AMUSE_local_data_dir = AMUSE_local_data_dir_in
         set_MESA_paths = 0
      end function set_MESA_paths
    
! Initialize the stellar evolution code
      integer function initialize_code()
         use amuse_support, only: failed, AMUSE_mesa_data_dir, AMUSE_inlist_path
         use run_star_support
         use ctrls_io, only: set_default_controls
         implicit none
         integer :: ierr
         initialize_code = -1
         !call set_default_controls
         call do_read_star_job(AMUSE_inlist_path, ierr)
         if (failed('do_read_star_job', ierr)) return
         ! Replace value of mesa_data_dir just read, with supplied path.
         mesa_data_dir = AMUSE_mesa_data_dir
         call star_init(mesa_data_dir, kappa_file_prefix, &
            net_reaction_filename, rates_dir, ppn_rate_numbers_fname, ierr)
         if (failed('star_init', ierr)) return
         profile_columns_file = trim(mesa_data_dir) // '/star_data/profile_columns.list'
         log_columns_file = trim(mesa_data_dir) // '/star_data/log_columns.list'
         call flush()
         report_backups = .true.
         report_retries = .true.
         initialize_code = 0
      end function initialize_code

      
      integer function commit_parameters()
         commit_parameters = 0
      end function commit_parameters
      
      integer function recommit_parameters()
         recommit_parameters = 0
      end function recommit_parameters
      
      integer function cleanup_code()
         cleanup_code = 0
      end function cleanup_code
      
! Create new ZAMS model for a different metallicity
   subroutine new_zams_model(ierr)
      use create_zams, only: AMUSE_do_create_zams
      use amuse_support
      use run_star_support, only: run_create_zams, zams_inlist
      implicit none
      integer :: ierr
      call get_zams_filename(AMUSE_zams_filename, ierr)
      if (failed('get_zams_filename', ierr)) return
      run_create_zams = .true. ! is this necessary?
      call AMUSE_do_create_zams(AMUSE_metallicity, AMUSE_zams_filename, &
         AMUSE_inlist_path, &
         AMUSE_dmass, AMUSE_mlo, AMUSE_mhi, ierr)
      if (failed('AMUSE_do_create_zams', ierr)) return      
      call flush()
      ierr = 0
   end subroutine new_zams_model

! Create a new particle
   function new_particle(AMUSE_id, AMUSE_mass)
      use amuse_support
      use star_lib, only: alloc_star, star_setup, star_load_zams, &
         show_terminal_header, show_log_description
      use star_private_def, only: star_info, get_star_ptr
      use run_star_support, only: setup_for_run_star, before_evolve, &
         show_log_description_at_start
      implicit none
      integer, intent(out) :: AMUSE_id
      integer :: new_particle, ierr
      double precision, intent(in) :: AMUSE_mass
      type (star_info), pointer :: s
      new_particle = -1
      AMUSE_id = alloc_star(ierr)
      number_of_particles = AMUSE_id
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
      s% RGB_wind_scheme = AMUSE_RGB_wind_scheme
      s% AGB_wind_scheme = AMUSE_AGB_wind_scheme
      s% Reimers_wind_eta = AMUSE_reimers_wind_efficiency
      s% Blocker_wind_eta = AMUSE_blocker_wind_efficiency
      s% de_Jager_wind_eta = AMUSE_de_jager_wind_efficiency
      s% Dutch_wind_eta = AMUSE_dutch_wind_efficiency
      if (debugging) then
         write (*,*) "Creating new particles with mass: ", s% initial_mass
         write (*,*) "Loading starting model from: ", s% zams_filename
      endif
      if (show_log_description_at_start) then
         write(*,*)
         call show_log_description(AMUSE_id, ierr)
         if (failed('show_log_description', ierr)) return
      end if
      write(*,*)
      call star_load_zams(AMUSE_id, ierr)
      if (failed('star_load_zams', ierr)) return
      call setup_for_run_star(AMUSE_id, s, .false., ierr)
      if (failed('setup_for_run_star', ierr)) return
      call before_evolve(AMUSE_id, ierr)
      if (failed('before_evolve', ierr)) return
      call show_terminal_header(AMUSE_id, ierr)
      if (failed('show_terminal_header', ierr)) return
      call flush()
      new_particle = 0
   end function

! Create a new pre-main-sequence star
   function new_pre_ms_particle(AMUSE_id, AMUSE_mass)
      use amuse_support
      use star_lib, only: alloc_star, star_setup, star_create_pre_ms_model, &
         show_terminal_header, show_log_description
      use star_private_def, only: star_info, get_star_ptr
      use run_star_support, only: setup_for_run_star, before_evolve, &
         show_log_description_at_start
      implicit none
      integer, intent(out) :: AMUSE_id
      integer :: new_pre_ms_particle, ierr
      double precision, intent(in) :: AMUSE_mass
      type (star_info), pointer :: s
      new_pre_ms_particle = -1
      AMUSE_id = alloc_star(ierr)
      number_of_particles = AMUSE_id
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
      s% RGB_wind_scheme = AMUSE_RGB_wind_scheme
      s% AGB_wind_scheme = AMUSE_AGB_wind_scheme
      s% Reimers_wind_eta = AMUSE_reimers_wind_efficiency
      s% Blocker_wind_eta = AMUSE_blocker_wind_efficiency
      s% de_Jager_wind_eta = AMUSE_de_jager_wind_efficiency
      s% Dutch_wind_eta = AMUSE_dutch_wind_efficiency
      if (debugging) then
         write (*,*) "Creating new pre-main-sequence particles with mass: ", s% initial_mass
      endif
      if (show_log_description_at_start) then
         write(*,*)
         call show_log_description(AMUSE_id, ierr)
         if (failed('show_log_description', ierr)) return
      end if
      write(*,*)
      call star_create_pre_ms_model(AMUSE_id, 0.0d0, 0.0d0, 0.0d0, ierr)
      if (failed('create_pre_ms_model', ierr)) return
      call setup_for_run_star(AMUSE_id, s, .false., ierr)
      if (failed('setup_for_run_star', ierr)) return
      call before_evolve(AMUSE_id, ierr)
      if (failed('before_evolve', ierr)) return
      call show_terminal_header(AMUSE_id, ierr)
      if (failed('show_terminal_header', ierr)) return
      call flush()
      new_pre_ms_particle = 0
   end function

! Remove a particle (doesn't do anything yet)
   function delete_star(AMUSE_id)
      implicit none
      integer, intent(in) :: AMUSE_id
      integer :: delete_star
      delete_star = 0
   end function

   function commit_particles()
      use amuse_support, only: target_times, number_of_particles
      implicit none
      integer :: commit_particles
      allocate(target_times(number_of_particles))
      target_times = 0
      commit_particles = 0
   end function

   function recommit_particles()
      use amuse_support, only: target_times, number_of_particles
      implicit none
      integer :: recommit_particles
      double precision, allocatable :: temp(:)
      allocate(temp(size(target_times)))
      temp = target_times
      deallocate(target_times)
      allocate(target_times(number_of_particles))
      target_times = 0
      target_times(1:size(temp)) = temp
      deallocate(temp)
      recommit_particles = 0
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
      character (len=1024) :: file
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
! Set the current mass of the star
   function set_mass(AMUSE_id, AMUSE_value)
      use star_private_def, only: star_info, get_star_ptr
      use const_def, only: msol
      use amuse_support, only: failed
      implicit none
      integer, intent(in) :: AMUSE_id
      double precision, intent(in) :: AMUSE_value
      integer :: set_mass, ierr
      type (star_info), pointer :: s
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) then
         set_mass = -1
      else
         s% mstar = AMUSE_value * msol
         s% mstar_old = AMUSE_value * msol
         s% mstar_older = AMUSE_value * msol
         s% star_mass = AMUSE_value
         set_mass = 0
      endif
   end function

! Return the current core mass of the star, where hydrogen abundance is <= h1_boundary_limit
   function get_core_mass(AMUSE_id, AMUSE_value)
      use star_private_def, only: star_info, get_star_ptr
      use amuse_support, only: failed
      implicit none
      integer, intent(in) :: AMUSE_id
      double precision, intent(out) :: AMUSE_value
      integer :: get_core_mass, ierr
      type (star_info), pointer :: s
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) then
         AMUSE_value = -1.0
         get_core_mass = -1
      else
         AMUSE_value = s% h1_boundary_mass
         get_core_mass = 0
      endif
   end function

! Return the current mass loss rate of the star
   function get_mass_loss_rate(AMUSE_id, AMUSE_value)
      use star_private_def, only: star_info, get_star_ptr
      use amuse_support, only: failed
      implicit none
      integer, intent(in) :: AMUSE_id
      double precision, intent(out) :: AMUSE_value
      integer :: get_mass_loss_rate, ierr
      type (star_info), pointer :: s
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) then
         AMUSE_value = -1.0
         get_mass_loss_rate = -1
      else
         AMUSE_value = -s% mstar_dot
         get_mass_loss_rate = 0
      endif
   end function

! Return the current user-specified mass transfer rate of the star
   function get_manual_mass_transfer_rate(AMUSE_id, AMUSE_value)
      use star_private_def, only: star_info, get_star_ptr
      use amuse_support, only: failed
      implicit none
      integer, intent(in) :: AMUSE_id
      double precision, intent(out) :: AMUSE_value
      integer :: get_manual_mass_transfer_rate, ierr
      type (star_info), pointer :: s
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) then
         AMUSE_value = -1.0
         get_manual_mass_transfer_rate = -1
      else
         AMUSE_value = s% mass_change
         get_manual_mass_transfer_rate = 0
      endif
   end function

! Set a new user-specified mass transfer rate of the star
   function set_manual_mass_transfer_rate(AMUSE_id, AMUSE_value)
      use star_private_def, only: star_info, get_star_ptr
      use amuse_support, only: failed
      implicit none
      integer, intent(in) :: AMUSE_id
      double precision, intent(out) :: AMUSE_value
      integer :: set_manual_mass_transfer_rate, ierr
      type (star_info), pointer :: s
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) then
         set_manual_mass_transfer_rate = -1
      else
         s% mass_change = AMUSE_value
         set_manual_mass_transfer_rate = 0
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
         AMUSE_value = s% dt_next/secyer
         get_time_step = 0
      end function

! Set the next timestep for the star
      function set_time_step(AMUSE_id, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use const_def, only: secyer
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id
         double precision, intent(in) :: AMUSE_value
         integer :: set_time_step, ierr
         type (star_info), pointer :: s
         set_time_step = -1
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) return
         s% dt_next = AMUSE_value*secyer
         set_time_step = 0
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
         use do_one_utils, only: do_show_terminal_header, do_terminal_summary
         use star_utils, only:eval_current_y, eval_current_z
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
         ! if the star is not a stellar remnant...
         if (s% log_surface_radius > -1.0) then
            ! Use the MESA phase_of_evolution marker (explained below)
            select case(s% phase_of_evolution)
               case(0)
                  AMUSE_value = 17 ! Pre-main-sequence star
               case(1,2)
                  if (s% star_mass < 0.75) then
                     AMUSE_value = 0 ! Convective low mass star
                  else
                     AMUSE_value = 1 ! Main sequence star
                  endif
               case(3)
                  AMUSE_value = 3 ! Red giant branch
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
                  call do_show_terminal_header(s)
                  call do_terminal_summary(s)
                  return
            end select
         else ! stellar remnant
            if (s% star_mass < 1.44) then ! white dwarf
               ! Helium White Dwarf:
               if (s% center_he3 + s% center_he4 > 0.1) AMUSE_value = 10
               ! Carbon/Oxygen White Dwarf:
               if (s% center_c12 > 0.1) AMUSE_value = 11
               ! Oxygen/Neon White Dwarf:
               if (s% center_ne20 > 0.01) AMUSE_value = 12
                ! Else? Unknown kind of white dwarf... hopefully never reached.
               if (AMUSE_value == -99) AMUSE_value = -10
            else
               if (s% star_mass < 3.2) then
                  AMUSE_value = 13 ! Neutron Star
               else
                  AMUSE_value = 14 ! Black Hole
               endif
            endif
         endif
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

! Return the current number of zones/mesh-cells of the star
      integer function get_number_of_zones(AMUSE_id, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id
         integer, intent(out) :: AMUSE_value
         integer :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            AMUSE_value = -1
            get_number_of_zones = -1
         else
            AMUSE_value = s% nz
            get_number_of_zones = 0
         endif
      end function

! Return the number_of_backups_in_a_row of the star
      integer function get_number_of_backups_in_a_row(AMUSE_id, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id
         integer, intent(out) :: AMUSE_value
         integer :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            get_number_of_backups_in_a_row = -1
         else
            AMUSE_value = s% number_of_backups_in_a_row
            get_number_of_backups_in_a_row = 0
         endif
      end function

! Reset number_of_backups_in_a_row of the star
      integer function reset_number_of_backups_in_a_row(AMUSE_id)
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id
         integer :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            reset_number_of_backups_in_a_row = -1
         else
            s% number_of_backups_in_a_row = 0
            reset_number_of_backups_in_a_row = 0
         endif
      end function

! Return the mass fraction at the specified zone/mesh-cell of the star
      integer function get_mass_fraction_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id, AMUSE_zone
         double precision, intent(out) :: AMUSE_value
         integer :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            AMUSE_value = -1.0
            get_mass_fraction_at_zone = -1
         else
            if (AMUSE_zone >= s% nz .or. AMUSE_zone < 0) then
                AMUSE_value = -1.0
                get_mass_fraction_at_zone = -2
            else
               if (s% number_of_backups_in_a_row > s% max_backups_in_a_row ) then
                  AMUSE_value = s% dq_old(s% nz - AMUSE_zone)
               else
                  AMUSE_value = s% dq(s% nz - AMUSE_zone)
               endif
               get_mass_fraction_at_zone = 0
            endif
         endif
      end function
! Set the mass fraction at the specified zone/mesh-cell of the star
      integer function set_mass_fraction_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id, AMUSE_zone
         double precision, intent(in) :: AMUSE_value
         integer :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            set_mass_fraction_at_zone = -1
         else
            if (AMUSE_zone >= s% nz .or. AMUSE_zone < 0) then
               set_mass_fraction_at_zone = -2
            else
               s% dq(s% nz - AMUSE_zone) = AMUSE_value
               set_mass_fraction_at_zone = 0
            endif
         endif
      end function

! Return the temperature at the specified zone/mesh-cell of the star
      integer function get_temperature_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id, AMUSE_zone
         double precision, intent(out) :: AMUSE_value
         integer :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            AMUSE_value = -1.0
            get_temperature_at_zone = -1
         else
            if (AMUSE_zone >= s% nz .or. AMUSE_zone < 0) then
                AMUSE_value = -1.0
                get_temperature_at_zone = -2
            else
               if (s% number_of_backups_in_a_row > s% max_backups_in_a_row ) then
                  AMUSE_value = exp(s% xs_old(s% i_lnT, s% nz - AMUSE_zone))
               else
                  AMUSE_value = exp(s% xs(s% i_lnT, s% nz - AMUSE_zone))
               endif
               get_temperature_at_zone = 0
            endif
         endif
      end function
! Set the temperature at the specified zone/mesh-cell of the star
      integer function set_temperature_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id, AMUSE_zone
         double precision, intent(in) :: AMUSE_value
         integer :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            set_temperature_at_zone = -1
         else
            if (AMUSE_zone >= s% nz .or. AMUSE_zone < 0) then
               set_temperature_at_zone = -2
            else
               s% xs(s% i_lnT, s% nz - AMUSE_zone) = log(AMUSE_value)
               set_temperature_at_zone = 0
            endif
         endif
      end function

! Return the density at the specified zone/mesh-cell of the star
      integer function get_density_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id, AMUSE_zone
         double precision, intent(out) :: AMUSE_value
         integer :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            AMUSE_value = -1.0
            get_density_at_zone = -1
         else
            if (AMUSE_zone >= s% nz .or. AMUSE_zone < 0) then
               AMUSE_value = -1.0
               get_density_at_zone = -2
            else
               if (s% number_of_backups_in_a_row > s% max_backups_in_a_row ) then
                  AMUSE_value = exp(s% xs_old(s% i_lnd, s% nz - AMUSE_zone))
               else
                  AMUSE_value = exp(s% xs(s% i_lnd, s% nz - AMUSE_zone))
               endif
               get_density_at_zone = 0
            endif
         endif
      end function
! Set the density at the specified zone/mesh-cell of the star
      integer function set_density_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id, AMUSE_zone
         double precision, intent(in) :: AMUSE_value
         integer :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            set_density_at_zone = -1
         else
            if (AMUSE_zone >= s% nz .or. AMUSE_zone < 0) then
               set_density_at_zone = -2
            else
               s% xs(s% i_lnd, s% nz - AMUSE_zone) = log(AMUSE_value)
               set_density_at_zone = 0
            endif
         endif
      end function

! Return the radius at the specified zone/mesh-cell of the star
      integer function get_radius_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id, AMUSE_zone
         double precision, intent(out) :: AMUSE_value
         integer :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            AMUSE_value = -1.0
            get_radius_at_zone = -1
         else
            if (AMUSE_zone >= s% nz .or. AMUSE_zone < 0) then
                AMUSE_value = -1.0
                get_radius_at_zone = -2
            else
               if (s% number_of_backups_in_a_row > s% max_backups_in_a_row ) then
                  AMUSE_value = exp(s% xs_old(s% i_lnR, s% nz - AMUSE_zone))
               else
                  AMUSE_value = exp(s% xs(s% i_lnR, s% nz - AMUSE_zone))
               endif
               get_radius_at_zone = 0
            endif
         endif
      end function
! Set the radius at the specified zone/mesh-cell of the star
      integer function set_radius_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id, AMUSE_zone
         double precision, intent(in) :: AMUSE_value
         integer :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            set_radius_at_zone = -1
         else
            if (AMUSE_zone >= s% nz .or. AMUSE_zone < 0) then
               set_radius_at_zone = -2
            else
               s% xs(s% i_lnR, s% nz - AMUSE_zone) = log(AMUSE_value)
               set_radius_at_zone = 0
            endif
         endif
      end function

! Return the luminosity at the specified zone/mesh-cell of the star
      integer function get_luminosity_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id, AMUSE_zone
         double precision, intent(out) :: AMUSE_value
         integer :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            AMUSE_value = -1.0
            get_luminosity_at_zone = -1
         else
            if (AMUSE_zone >= s% nz .or. AMUSE_zone < 0) then
                AMUSE_value = -1.0
                get_luminosity_at_zone = -2
            else
               if (s% number_of_backups_in_a_row > s% max_backups_in_a_row ) then
                  AMUSE_value = s% xs_old(s% i_lum, s% nz - AMUSE_zone)
               else
                  AMUSE_value = s% xs(s% i_lum, s% nz - AMUSE_zone)
               endif
               get_luminosity_at_zone = 0
            endif
         endif
      end function
! Set the luminosity at the specified zone/mesh-cell of the star
      integer function set_luminosity_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id, AMUSE_zone
         double precision, intent(in) :: AMUSE_value
         integer :: ierr
         type (star_info), pointer :: s
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            set_luminosity_at_zone = -1
         else
            if (AMUSE_zone >= s% nz .or. AMUSE_zone < 0) then
               set_luminosity_at_zone = -2
            else
               s% xs(s% i_lum, s% nz - AMUSE_zone) = AMUSE_value
               set_luminosity_at_zone = 0
            endif
         endif
      end function

! Return the mean molecular weight per particle (ions + free electrons) at the specified zone/mesh-cell of the star
      integer function get_mu_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
!         use micro, only: do_eos_for_cell
         use chem_def, only: ih1, ihe3, ihe4
         
         use eos_lib, only: eosDT_get
         use eos_def, only: num_eos_basic_results, i_mu
         use const_def, only: ln10
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id, AMUSE_zone
         double precision, intent(out) :: AMUSE_value
         integer, pointer :: net_iso(:)
         integer :: ierr, k
         type (star_info), pointer :: s
         double precision, dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
         double precision :: z, xh, xhe, abar, zbar
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            AMUSE_value = -1.0
            get_mu_at_zone = -1
         else
            if (AMUSE_zone >= s% nz .or. AMUSE_zone < 0) then
                AMUSE_value = -1.0
                get_mu_at_zone = -2
            else
               k = s% nz - AMUSE_zone
               if (s% number_of_backups_in_a_row > s% max_backups_in_a_row ) then
                  call get_abar_zbar(s, k, abar, zbar)
                  net_iso => s% net_iso
                  xh = s% xa_old(net_iso(ih1),k)
                  xhe = s% xa_old(net_iso(ihe3),k) + s% xa_old(net_iso(ihe4),k)
                  z = max(0d0,1d0-(xh+xhe))
                  call eosDT_get( &
                     s% eos_handle, z, xh, abar, zbar, &
                     exp(s% xs_old(s% i_lnd, k)), s% xs_old(s% i_lnd, k)/ln10, &
                     exp(s% xs_old(s% i_lnT, k)), s% xs_old(s% i_lnT, k)/ln10, &
                     res, d_dlnd, d_dlnT, ierr)
                  if (failed('eosDT_get', ierr)) then
                     AMUSE_value = -1.0
                     get_mu_at_zone = -4
                     return
                  endif
                  s% mu(k) = res(i_mu)
               endif
               AMUSE_value = s% mu(k)
               get_mu_at_zone = 0
            endif
         endif
!      end function
         
         contains
         
         subroutine get_abar_zbar(s, k, abar, zbar)
!            use star_private_def, only: star_info
            use chem_lib, only: composition_info
            type (star_info), pointer :: s
            integer, intent(in) :: k
            double precision, intent(out) :: abar, zbar
            double precision :: z2bar, ye, xsum, dabar_dx(s% species), dzbar_dx(s% species)
            integer :: species
            species = s% species
            call composition_info(species, s% chem_id, s% xa_old(1:species,k), &
                abar, zbar, z2bar, ye, xsum, dabar_dx, dzbar_dx)
         end subroutine get_abar_zbar
         
      end function


! Return the total (gas + radiation) pressure at the specified zone/mesh-cell of the star
      integer function get_pressure_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
         use star_private_def, only: star_info, get_star_ptr
         use chem_def, only: ih1, ihe3, ihe4
         use eos_lib, only: eosDT_get
         use eos_def, only: num_eos_basic_results, i_lnPgas
         use const_def, only: ln10, crad
         use amuse_support, only: failed, debugging
         
         implicit none
         integer, intent(in) :: AMUSE_id, AMUSE_zone
         double precision, intent(out) :: AMUSE_value
         integer, pointer :: net_iso(:)
         integer :: ierr, k
         type (star_info), pointer :: s
         double precision, dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
         double precision :: z, xh, xhe, abar, zbar
         
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) then
            AMUSE_value = -1.0
            get_pressure_at_zone = -1
         else
            k = s% nz - AMUSE_zone
            if (s% number_of_backups_in_a_row > s% max_backups_in_a_row ) then
!               if (debugging) then
!                  write(*, *) 'Warning: pressure may not be up to date, since the last evolve failed.'
!               endif
                  call get_abar_zbar(s, k, abar, zbar)
                  net_iso => s% net_iso
                  xh = s% xa_old(net_iso(ih1),k)
                  xhe = s% xa_old(net_iso(ihe3),k) + s% xa_old(net_iso(ihe4),k)
                  z = max(0d0,1d0-(xh+xhe))
                  call eosDT_get( &
                     s% eos_handle, z, xh, abar, zbar, &
                     exp(s% xs_old(s% i_lnd, k)), s% xs_old(s% i_lnd, k)/ln10, &
                     exp(s% xs_old(s% i_lnT, k)), s% xs_old(s% i_lnT, k)/ln10, &
                     res, d_dlnd, d_dlnT, ierr)
                  if (failed('eosDT_get', ierr)) then
                     AMUSE_value = -1.0
                     get_pressure_at_zone = -4
                     return
                  endif
!                  s% mu(k) = res(i_mu)
                  s% lnPgas(k) = res(i_lnPgas)
                  s% Pgas(k) = exp(s% lnPgas(k))
                  
                  s% Prad(k) = crad * s% T(k)**4 / 3
                  s% P(k) = s% Prad(k) + s% Pgas(k)
!               get_pressure_at_zone = -1
!               return
            endif
            
            if (AMUSE_zone >= s% nz .or. AMUSE_zone < 0) then
                AMUSE_value = -1.0
                get_pressure_at_zone = -2
            else
               AMUSE_value = s% P(k)
               get_pressure_at_zone = 0
            endif
         endif
         
         contains
         
         subroutine get_abar_zbar(s, k, abar, zbar)
!            use star_private_def, only: star_info
            use chem_lib, only: composition_info
            type (star_info), pointer :: s
            integer, intent(in) :: k
            double precision, intent(out) :: abar, zbar
            double precision :: z2bar, ye, xsum, dabar_dx(s% species), dzbar_dx(s% species)
            integer :: species
            species = s% species
            call composition_info(species, s% chem_id, s% xa_old(1:species,k), &
                abar, zbar, z2bar, ye, xsum, dabar_dx, dzbar_dx)
         end subroutine get_abar_zbar
         
      end function

! Return the current number of chemical abundance variables per zone of the star
   integer function get_number_of_species(AMUSE_id, AMUSE_value)
      use star_private_def, only: star_info, get_star_ptr
      use amuse_support, only: failed
      implicit none
      integer, intent(in) :: AMUSE_id
      integer, intent(out) :: AMUSE_value
      integer :: ierr
      type (star_info), pointer :: s
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) then
         AMUSE_value = -1
         get_number_of_species = -1
      else
         AMUSE_value = s% nvar_chem
         get_number_of_species = 0
      endif
   end function

! Return the name of chemical abundance variable 'AMUSE_species' of the star
   integer function get_name_of_species(AMUSE_id, AMUSE_species, AMUSE_value)
      use star_private_def, only: star_info, get_star_ptr
      use amuse_support, only: failed
      use chem_def, only: num_chem_isos, chem_isos
      implicit none
      integer, intent(in) :: AMUSE_id, AMUSE_species
      character (len=6), intent(out) :: AMUSE_value
      integer :: ierr
      type (star_info), pointer :: s
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) then
         AMUSE_value = 'error'
         get_name_of_species = -1
      else if (AMUSE_species > s% nvar_chem .or. AMUSE_species < 1) then
         AMUSE_value = 'error'
         get_name_of_species = -3
      else
         AMUSE_value = chem_isos% name(s% chem_id(AMUSE_species))
         get_name_of_species = 0
      endif
!      do ierr=1,s% nvar
!         write(*,*) ierr, s% nameofvar(ierr)
!      end do
!      do ierr=1,num_chem_isos
!         write(*,*) ierr, s% net_iso(ierr), chem_isos% name(ierr)
!      end do
!      do ierr=1,s% nvar_chem
!         write(*,*) ierr, s% chem_id(ierr), chem_isos% name(s% chem_id(ierr))
!      end do
   end function

! Return the chem_ID of chemical abundance variable 'AMUSE_species' of the star
   integer function get_id_of_species(AMUSE_id, AMUSE_species, AMUSE_value)
      use star_private_def, only: star_info, get_star_ptr
      use amuse_support, only: failed
      implicit none
      integer, intent(in) :: AMUSE_id, AMUSE_species
      integer, intent(out) :: AMUSE_value
      integer :: ierr
      type (star_info), pointer :: s
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) then
         AMUSE_value = -1
         get_id_of_species = -1
      else if (AMUSE_species > s% nvar_chem .or. AMUSE_species < 1) then
         AMUSE_value = -1
         get_id_of_species = -3
      else
         AMUSE_value = s% chem_id(AMUSE_species)
         get_id_of_species = 0
      endif
   end function

! Return the mass number of chemical abundance variable 'AMUSE_species' of the star
   integer function get_mass_of_species(AMUSE_id, AMUSE_species, AMUSE_value)
      use star_private_def, only: star_info, get_star_ptr
      use amuse_support, only: failed
      use chem_def, only: chem_isos
      use const_def, only: mev_to_ergs, clight, amu
      implicit none
      integer, intent(in) :: AMUSE_id, AMUSE_species
      double precision, intent(out) :: AMUSE_value
      integer :: ierr
      type (star_info), pointer :: s
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) then
         AMUSE_value = -1.0
         get_mass_of_species = -1
      else if (AMUSE_species > s% nvar_chem .or. AMUSE_species < 1) then
         AMUSE_value = -1.0
         get_mass_of_species = -3
      else
         AMUSE_value = chem_isos% A(s% chem_id(AMUSE_species)) + &
            chem_isos% mass_excess(s% chem_id(AMUSE_species))*mev_to_ergs/(clight*clight*amu)
         get_mass_of_species = 0
      endif
   end function

! Return the mass fraction of species 'AMUSE_species' at the specified 
! zone/mesh-cell of the star
   integer function get_mass_fraction_of_species_at_zone(AMUSE_id, &
         AMUSE_species, AMUSE_zone, AMUSE_value)
      use star_private_def, only: star_info, get_star_ptr
      use amuse_support, only: failed
      implicit none
      integer, intent(in) :: AMUSE_id, AMUSE_zone, AMUSE_species
      double precision, intent(out) :: AMUSE_value
      integer :: ierr
      type (star_info), pointer :: s
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) then
         AMUSE_value = -1.0
         get_mass_fraction_of_species_at_zone = -1
      else if (AMUSE_zone >= s% nz .or. AMUSE_zone < 0) then
         AMUSE_value = -1.0
         get_mass_fraction_of_species_at_zone = -2
      else if (AMUSE_species > s% nvar_chem .or. AMUSE_species < 1) then
         AMUSE_value = -1.0
         get_mass_fraction_of_species_at_zone = -3
      else
         if (s% number_of_backups_in_a_row > s% max_backups_in_a_row ) then
            AMUSE_value = s% xa_old(AMUSE_species, s% nz - AMUSE_zone)
         else
            AMUSE_value = s% xa(AMUSE_species, s% nz - AMUSE_zone)
         endif
         get_mass_fraction_of_species_at_zone = 0
      endif
   end function
! Set the mass fraction of species 'AMUSE_species' at the specified 
! zone/mesh-cell of the star
   integer function set_mass_fraction_of_species_at_zone(AMUSE_id, &
         AMUSE_species, AMUSE_zone, AMUSE_value)
      use star_private_def, only: star_info, get_star_ptr
      use amuse_support, only: failed
      implicit none
      integer, intent(in) :: AMUSE_id, AMUSE_zone, AMUSE_species
      double precision, intent(in) :: AMUSE_value
      integer :: ierr
      type (star_info), pointer :: s
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) then
         set_mass_fraction_of_species_at_zone = -1
      else if (AMUSE_zone >= s% nz .or. AMUSE_zone < 0) then
         set_mass_fraction_of_species_at_zone = -2
      else if (AMUSE_species > s% nvar_chem .or. AMUSE_species < 1) then
         set_mass_fraction_of_species_at_zone = -3
      else
         s% xa(AMUSE_species, s% nz - AMUSE_zone) = AMUSE_value
         s% xa_pre_hydro(AMUSE_species, s% nz - AMUSE_zone) = AMUSE_value
         set_mass_fraction_of_species_at_zone = 0
      endif
   end function

! Erase memory of the star - xs_old(er), xa_old(er), q_old(er), etc.
! Useful after setting the stucture of the star, to prevent backup steps to undo changes
   integer function erase_memory(AMUSE_id)
      use star_private_def, only: star_info, get_star_ptr
      use amuse_support, only: failed
      implicit none
      integer, intent(in) :: AMUSE_id
      integer :: ierr
      type (star_info), pointer :: s
      
      erase_memory = -1
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) return
      if (s%generations > 1) then
         s% nz_old = s% nz
         call realloc2d_if_necessary(s% xa_old, s% species, s% nz, ierr)
         if (failed('realloc2d_if_necessary', ierr)) return
         s% xa_old(:,:) = s% xa(:,:)
         call realloc2d_if_necessary(s% xs_old, s% nvar, s% nz, ierr)
         if (failed('realloc2d_if_necessary', ierr)) return
         s% xs_old(:,:) = s% xs(:,:)
         call realloc1d_if_necessary(s% q_old, s% nz, ierr)
         if (failed('realloc1d_if_necessary', ierr)) return
         s% q_old(:) = s% q(:)
         call realloc1d_if_necessary(s% dq_old, s% nz, ierr)
         if (failed('realloc1d_if_necessary', ierr)) return
         s% dq_old(:) = s% dq(:)
         if (s%generations == 3) then
            s% nz_older = s% nz
            call realloc2d_if_necessary(s% xa_older, s% species, s% nz, ierr)
            if (failed('realloc2d_if_necessary', ierr)) return
            s% xa_older(:,:) = s% xa(:,:)
            call realloc2d_if_necessary(s% xs_older, s% nvar, s% nz, ierr)
            if (failed('realloc2d_if_necessary', ierr)) return
            s% xs_older(:,:) = s% xs(:,:)
            call realloc1d_if_necessary(s% q_older, s% nz, ierr)
            if (failed('realloc1d_if_necessary', ierr)) return
            s% q_older(:) = s% q(:)
            call realloc1d_if_necessary(s% dq_older, s% nz, ierr)
            if (failed('realloc1d_if_necessary', ierr)) return
            s% dq_older(:) = s% dq(:)
         end if
      end if
      erase_memory = 0
      
      contains
      
      subroutine realloc1d_if_necessary(ptr,new_size,ierr)
         double precision, pointer :: ptr(:)
         integer, intent(in) :: new_size
         integer, intent(out) :: ierr
         ierr = 0
         if (associated(ptr)) then
            if (size(ptr,1) == new_size) return
            deallocate(ptr)
         end if
         allocate(ptr(new_size),stat=ierr)
      end subroutine realloc1d_if_necessary
      
      subroutine realloc2d_if_necessary(ptr,ld,new_size,ierr)
         double precision, pointer :: ptr(:,:)
         integer, intent(in) :: ld, new_size
         integer, intent(out) :: ierr
         ierr = 0
         if (associated(ptr)) then
            if (size(ptr,1) == ld .and. size(ptr,2) == new_size) return
            deallocate(ptr)
         end if
         allocate(ptr(ld,new_size),stat=ierr)
      end subroutine realloc2d_if_necessary
  
   end function erase_memory

! Evolve the star for one step
   function evolve_one_step(AMUSE_id)
      use star_private_def, only: star_info, get_star_ptr
      use run_star_support
      use run_star, only: check_model
      use run_star_extras, only: &
         how_many_extra_profile_columns, data_for_extra_profile_columns, &
         how_many_extra_log_columns, data_for_extra_log_columns
      use amuse_support, only: evolve_failed
      use const_def, only: secyer
      implicit none
      integer, intent(in) :: AMUSE_id
      integer :: evolve_one_step
      type (star_info), pointer :: s
      integer :: ierr, model_number, result, result_reason
      logical :: first_try
      evolve_one_step = -1
      call get_star_ptr(AMUSE_id, s, ierr)
      if (evolve_failed('get_star_ptr', ierr, evolve_one_step, -1)) return
      if (auto_extend_net) then
         call extend_net(s, ierr)
         if (evolve_failed('extend_net', ierr, evolve_one_step, -2)) return
      end if
      first_try = .true.
      model_number = get_model_number(AMUSE_id, ierr)
      if (evolve_failed('get_model_number', ierr, evolve_one_step, -3)) return
      step_loop: do ! may need to repeat this loop for retry or backup
         result = star_evolve_step(AMUSE_id, first_try)
         if (result == keep_going) result = check_model(s, AMUSE_id, 0)
         if (result == keep_going) result = star_pick_next_timestep(AMUSE_id)
         if (result == keep_going) exit step_loop
         model_number = get_model_number(AMUSE_id, ierr)
         if (evolve_failed('get_model_number', ierr, evolve_one_step, -3)) return
         result_reason = get_result_reason(AMUSE_id, ierr)
         if (result == retry) then
            if (evolve_failed('get_result_reason', ierr, evolve_one_step, -4)) return
            if (report_retries) &
               write(*,'(i6,3x,a,/)') model_number, &
                  'retry reason ' // trim(result_reason_str(result_reason))
         else if (result == backup) then
            if (evolve_failed('get_result_reason', ierr, evolve_one_step, -4)) return
            if (report_backups) &
               write(*,'(i6,3x,a,/)') model_number, &
                  'backup reason ' // trim(result_reason_str(result_reason))
         end if
         if (result == retry) result = star_prepare_for_retry(AMUSE_id)
         if (result == backup) result = star_do1_backup(AMUSE_id)
         if (result == terminate) then
            evolve_one_step = -11 ! Unspecified stop condition reached, or:
            if (s% dt_next < s% min_timestep_limit) &
               evolve_one_step = -15 ! minimum timestep limit reached
            if (s% number_of_backups_in_a_row > s% max_backups_in_a_row ) &
               evolve_one_step = -14 ! max backups reached
            if (s% max_model_number > 0 .and. s% model_number >= s% max_model_number) &
               evolve_one_step = -13 ! max iterations reached
            if (s% star_age >= s% max_age) &
               evolve_one_step = -12 ! max_age reached
            return
         end if
         first_try = .false.
      end do step_loop
      result = star_finish_step(AMUSE_id, 0, .false., &
                     how_many_extra_profile_columns, data_for_extra_profile_columns, &
                     how_many_extra_log_columns, data_for_extra_log_columns, ierr)
      if (evolve_failed('star_finish_step', ierr, evolve_one_step, -5)) return
      if (s% model_number == save_model_number) then
         call star_write_model(AMUSE_id, save_model_filename, .true., ierr)
         if (evolve_failed('star_write_model', ierr, evolve_one_step, -6)) return
         write(*, *) 'saved to ' // trim(save_model_filename)
      end if
      evolve_one_step = 0
      call flush()
   end function

! Evolve the star until AMUSE_end_time
   integer function evolve_for(AMUSE_id, AMUSE_delta_t)
      use star_private_def, only: star_info, get_star_ptr
      use const_def, only: secyer
      use amuse_support, only: evolve_failed, target_times
      implicit none
      integer, intent(in) :: AMUSE_id
      double precision, intent(in) :: AMUSE_delta_t
      type (star_info), pointer :: s
      integer :: ierr
      integer :: evolve_one_step
      
      evolve_for = 0
      call get_star_ptr(AMUSE_id, s, ierr)
      if (evolve_failed('get_star_ptr', ierr, evolve_for, -1)) return
      
      target_times(AMUSE_id) = target_times(AMUSE_id) + AMUSE_delta_t * secyer
      
      evolve_loop: do while(evolve_for == 0 .and. &
            (s% time + s% min_timestep_limit < target_times(AMUSE_id))) ! evolve one step per loop
         evolve_for = evolve_one_step(AMUSE_id)
      end do evolve_loop
   end function evolve_for

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

! Return the wind (mass loss) scheme for RGB stars
      integer function get_RGB_wind_scheme(AMUSE_value)
         use amuse_support, only: AMUSE_RGB_wind_scheme
         implicit none
         integer, intent(out) :: AMUSE_value
         AMUSE_value = AMUSE_RGB_wind_scheme
         get_RGB_wind_scheme = 0
      end function get_RGB_wind_scheme

! Set the wind (mass loss) scheme for RGB stars
      integer function set_RGB_wind_scheme(AMUSE_value)
         use amuse_support, only: AMUSE_RGB_wind_scheme
         implicit none
         integer, intent(in) :: AMUSE_value
         AMUSE_RGB_wind_scheme = AMUSE_value
         set_RGB_wind_scheme = 0
      end function set_RGB_wind_scheme

! Return the wind (mass loss) scheme for AGB stars
      integer function get_AGB_wind_scheme(AMUSE_value)
         use amuse_support, only: AMUSE_AGB_wind_scheme
         implicit none
         integer, intent(out) :: AMUSE_value
         AMUSE_value = AMUSE_AGB_wind_scheme
         get_AGB_wind_scheme = 0
      end function get_AGB_wind_scheme

! Set the wind (mass loss) scheme for AGB stars
      integer function set_AGB_wind_scheme(AMUSE_value)
         use amuse_support, only: AMUSE_AGB_wind_scheme
         implicit none
         integer, intent(in) :: AMUSE_value
         AMUSE_AGB_wind_scheme = AMUSE_value
         set_AGB_wind_scheme = 0
      end function set_AGB_wind_scheme

! Retrieve the current value of the Reimers wind (mass loss) efficiency
      integer function get_reimers_wind_efficiency(AMUSE_value)
         use amuse_support, only: AMUSE_reimers_wind_efficiency
         implicit none
         double precision, intent(out) :: AMUSE_value
         AMUSE_value = AMUSE_reimers_wind_efficiency
         get_reimers_wind_efficiency = 0
      end function get_reimers_wind_efficiency

! Set the current value of the Reimers wind (mass loss) efficiency
      integer function set_reimers_wind_efficiency(AMUSE_value)
         use amuse_support, only: AMUSE_reimers_wind_efficiency
         implicit none
         double precision, intent(in) :: AMUSE_value
         AMUSE_reimers_wind_efficiency = AMUSE_value
         set_reimers_wind_efficiency = 0
      end function set_reimers_wind_efficiency

! Retrieve the current value of the Blocker wind (mass loss) efficiency
      integer function get_blocker_wind_efficiency(AMUSE_value)
         use amuse_support, only: AMUSE_blocker_wind_efficiency
         implicit none
         double precision, intent(out) :: AMUSE_value
         AMUSE_value = AMUSE_blocker_wind_efficiency
         get_blocker_wind_efficiency = 0
      end function get_blocker_wind_efficiency

! Set the current value of the Blocker wind (mass loss) efficiency
      integer function set_blocker_wind_efficiency(AMUSE_value)
         use amuse_support, only: AMUSE_blocker_wind_efficiency
         implicit none
         double precision, intent(in) :: AMUSE_value
         AMUSE_blocker_wind_efficiency = AMUSE_value
         set_blocker_wind_efficiency = 0
      end function set_blocker_wind_efficiency

! Retrieve the current value of the de Jager wind (mass loss) efficiency
      integer function get_de_jager_wind_efficiency(AMUSE_value)
         use amuse_support, only: AMUSE_de_jager_wind_efficiency
         implicit none
         double precision, intent(out) :: AMUSE_value
         AMUSE_value = AMUSE_de_jager_wind_efficiency
         get_de_jager_wind_efficiency = 0
      end function get_de_jager_wind_efficiency

! Set the current value of the de Jager wind (mass loss) efficiency
      integer function set_de_jager_wind_efficiency(AMUSE_value)
         use amuse_support, only: AMUSE_de_jager_wind_efficiency
         implicit none
         double precision, intent(in) :: AMUSE_value
         AMUSE_de_jager_wind_efficiency = AMUSE_value
         set_de_jager_wind_efficiency = 0
      end function set_de_jager_wind_efficiency

! Retrieve the current value of the Dutch wind (mass loss) efficiency
      integer function get_dutch_wind_efficiency(AMUSE_value)
         use amuse_support, only: AMUSE_dutch_wind_efficiency
         implicit none
         double precision, intent(out) :: AMUSE_value
         AMUSE_value = AMUSE_dutch_wind_efficiency
         get_dutch_wind_efficiency = 0
      end function get_dutch_wind_efficiency

! Set the current value of the Dutch wind (mass loss) efficiency
      integer function set_dutch_wind_efficiency(AMUSE_value)
         use amuse_support, only: AMUSE_dutch_wind_efficiency
         implicit none
         double precision, intent(in) :: AMUSE_value
         AMUSE_dutch_wind_efficiency = AMUSE_value
         set_dutch_wind_efficiency = 0
      end function set_dutch_wind_efficiency

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

! Retrieve the current value of the do_stabilize_new_stellar_model flag
      integer function get_stabilize_new_stellar_model_flag(AMUSE_value)
         use amuse_support, only: do_stabilize_new_stellar_model
         implicit none
         integer, intent(out) :: AMUSE_value
         if (do_stabilize_new_stellar_model) then
            AMUSE_value = 1
         else
            AMUSE_value = 0
         end if
         get_stabilize_new_stellar_model_flag = 0
      end function get_stabilize_new_stellar_model_flag

! Set the current value of the do_stabilize_new_stellar_model flag
      integer function set_stabilize_new_stellar_model_flag(AMUSE_value)
         use amuse_support, only: do_stabilize_new_stellar_model
         implicit none
         integer, intent(in) :: AMUSE_value
         if (AMUSE_value /= 0) then
            do_stabilize_new_stellar_model = .true.
         else
            do_stabilize_new_stellar_model = .false.
         end if
         set_stabilize_new_stellar_model_flag = 0
      end function set_stabilize_new_stellar_model_flag

! Retrieve the maximum number of stars that can be allocated in the code
      integer function get_maximum_number_of_stars(AMUSE_value)
         use star_def, only: max_star_handles
         implicit none
         integer, intent(out) :: AMUSE_value
         AMUSE_value = max_star_handles
         get_maximum_number_of_stars = 0
      end function get_maximum_number_of_stars

! Create a new particle from a user supplied model (non-ZAMS, e.g. merger product)
   integer function new_specified_stellar_model(d_mass, radius, rho, temperature, luminosity, &
         XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, n)
      use amuse_support
      use star_lib, only: alloc_star, star_setup, show_terminal_header
      use star_private_def, only: star_info, get_star_ptr
      use run_star_support, only: setup_for_run_star, before_evolve
      use read_model, only: set_zero_age_params, finish_load_model
      use alloc, only: set_var_info, set_q_flag, allocate_star_info_arrays
      use micro, only: init_mesa_micro
      use init_model, only: get_zams_model
      use chem_lib, only: get_nuclide_index
      use star_utils, only: set_qs, set_q_vars
      use do_one_utils, only: set_phase_of_evolution
      use evolve_support, only: yrs_for_init_timestep
      use const_def, only: secyer, Msun, Lsun
      
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: d_mass(n), radius(n), rho(n), &
         temperature(n), luminosity(n), XH(n), XHE(n), XC(n), XN(n), &
         XO(n), XNE(n), XMG(n), XSI(n), XFE(n)
      double precision :: x(n)
      integer :: ierr, k
      type (star_info), pointer :: s
      
      if (new_model_defined) then
         new_specified_stellar_model = -30
         return
      endif
      
      new_specified_stellar_model = -1
      id_new_model = alloc_star(ierr)
      if (failed('alloc_star', ierr)) return
      call get_star_ptr(id_new_model, s, ierr)
      if (failed('get_star_ptr', ierr)) return
      call star_setup(id_new_model, AMUSE_inlist_path, ierr)
      if (failed('star_setup', ierr)) return
      ! Replace value of mass and metallicity just read, with supplied values.
      s% initial_mass = sum(d_mass)
      s% initial_z = AMUSE_metallicity
      s% zams_filename = trim(AMUSE_zams_filename) // '.data'
      s% max_age = AMUSE_max_age_stop_condition
      s% min_timestep_limit = AMUSE_min_timestep_stop_condition
      s% max_model_number = AMUSE_max_iter_stop_condition
      s% mixing_length_alpha = AMUSE_mixing_length_ratio
      s% alpha_semiconvection = AMUSE_semi_convection_efficiency
      s% RGB_wind_scheme = AMUSE_RGB_wind_scheme
      s% AGB_wind_scheme = AMUSE_AGB_wind_scheme
      s% Reimers_wind_eta = AMUSE_reimers_wind_efficiency
      s% Blocker_wind_eta = AMUSE_blocker_wind_efficiency
      s% de_Jager_wind_eta = AMUSE_de_jager_wind_efficiency
      s% Dutch_wind_eta = AMUSE_dutch_wind_efficiency
      
      s% doing_first_model_of_run = .true.
      s% dt = 0
      s% dt_old = 0
      call set_zero_age_params(s)
            s% net_name = 'basic.net'
      s% species = 0
      s% v_flag = .false.
      s% q_flag = .false.
      s% mstar = s% initial_mass*Msun
      call set_var_info(s, ierr)
      call init_mesa_micro(s, ierr) ! uses s% net_name
      s% generations = 1
      
      if (n > s% max_allowed_nz) s% max_allowed_nz = n
      s% nz = n
      call allocate_star_info_arrays(s, ierr)
      if (failed('allocate_star_info_arrays', ierr)) return
      s% xs(s% i_lnd, :) = log(rho(:))
      s% xs(s% i_lnT, :) = log(temperature(:))
      s% xs(s% i_lnR, :) = log(radius(:))
      if (luminosity(1) <= 0) then
         ! No luminosities provided, make an educated guess
         do k = 1, s% nz - 3
            if (temperature(k) .gt. 1.0e7) exit
         end do
         if (debugging) write(*,*) "temperature(", k, ") = ", temperature(k)
         if (debugging) write(*,*) "radius(", k, ") = ", radius(k)
         x = radius / radius(k)
         if (debugging) write(*,*) "x(", k, ") = ", x(k), x(1), x(s% nz)
         s% xs(s% i_lum, :) = Lsun * s% initial_mass**3.5 * (1.0 - (1.0 + x) * exp(-x**2 - x))
      else
         s% xs(s% i_lum, :) = luminosity(:)
      endif
      s% dq(:) = d_mass(:) / s% initial_mass
      s% xa(s% net_iso(get_nuclide_index('h1')), :) = XH(:)
      s% xa(s% net_iso(get_nuclide_index('he3')), :) = 0.0d0
      s% xa(s% net_iso(get_nuclide_index('he4')), :) = XHE(:)
      s% xa(s% net_iso(get_nuclide_index('c12')), :) = XC(:)
      s% xa(s% net_iso(get_nuclide_index('n14')), :) = XN(:)
      s% xa(s% net_iso(get_nuclide_index('o16')), :) = XO(:)
      s% xa(s% net_iso(get_nuclide_index('ne20')), :) = XNE(:)
      s% xa(s% net_iso(get_nuclide_index('mg24')), :) = XMG(:) + XSI(:) + XFE(:) ! basic net for now...
      s% prev_Lmax = maxval(abs(s% xs(s% i_lum, 1:n)))
      call set_qs(s% nz, s% q, s% dq, ierr)
      if (failed('set_qs', ierr)) return
      if (s% q_flag) call set_q_vars(s)
      
      s% dt_next = yrs_for_init_timestep(s)*secyer
      s% dxs(:,:) = 0
      !
      s% extra_heat(:) = 0
      s% rate_factors(:) = 1
      call finish_load_model(s, ierr)
      call set_phase_of_evolution(s)
      if (s% q_flag) call set_q_flag(s% id, s% q_flag, ierr)
      
      call setup_for_run_star(id_new_model, s, .false., ierr)
      if (failed('setup_for_run_star', ierr)) return
      call before_evolve(id_new_model, ierr)
      if (failed('before_evolve', ierr)) return
      if (debugging) s% trace_evolve = .true.
      if (debugging) s% report_ierr = .true.
      call show_terminal_header(id_new_model, ierr)
      if (failed('show_terminal_header', ierr)) return
      call flush()
      new_model_defined = .true.
      new_specified_stellar_model = 0
   end function new_specified_stellar_model

   integer function new_stellar_model(d_mass, radius, rho, temperature, luminosity, &
         XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, n)
      use amuse_support
      use star_lib, only: alloc_star, star_setup, show_terminal_header
      use star_def, only: result_reason_str
      use star_private_def, only: star_info, get_star_ptr
      use run_star_support, only: setup_for_run_star, before_evolve
      use read_model, only: set_zero_age_params, finish_load_model
      use alloc, only: set_var_info, set_q_flag, allocate_star_info_arrays
      use micro, only: init_mesa_micro
      use init_model, only: get_zams_model
      use chem_lib, only: get_nuclide_index
      use star_utils, only: set_qs, set_q_vars
      use do_one_utils, only: set_phase_of_evolution
      use evolve_support, only: yrs_for_init_timestep
      use mesh_adjust, only: do_mesh_adjust
      use adjust_mesh, only: remesh
      use hydro_eqns, only: P_eqn_phot
      use hydro_vars, only: set_vars
      use const_def, only: secyer, Msun, Lsun
      use star_utils, only: set_xqs
      
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: d_mass(n), radius(n), rho(n), &
         temperature(n), luminosity(n), XH(n), XHE(n), XC(n), XN(n), &
         XO(n), XNE(n), XMG(n), XSI(n), XFE(n)
      double precision :: total_mass, original_timestep, f
      double precision :: original_timestep_limit, original_dxdt_nuc_f
      integer :: new_particle, ierr, tmp1_id_new_model, tmp2_id_new_model, &
         new_specified_stellar_model, finalize_stellar_model, match_mesh, &
         evolve_one_step, erase_memory, index_low, k1, k2
      logical :: do_T = .false.
      logical :: do_restore_timestep = .false.
      type (star_info), pointer :: s, s_tmp
      
      if (new_model_defined) then
         new_stellar_model = -30
         return
      endif
      
      ! *** Define a temporary star with the target 'new' structure: ***
      new_stellar_model = new_specified_stellar_model(d_mass, radius, rho, &
         temperature, luminosity, XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, n)
      if (failed('new_specified_stellar_model', new_stellar_model)) return
      
      if (do_stabilize_new_stellar_model) then
         new_stellar_model = -1
         if (debugging) write(*,*) 'tmp1_id_new_model', tmp1_id_new_model
         ierr = finalize_stellar_model(tmp1_id_new_model, 0.0d0)
         if (failed('finalize_stellar_model', ierr)) return
         call get_star_ptr(tmp1_id_new_model, s_tmp, ierr)
         if (failed('get_star_ptr', ierr)) return
         if (debugging) write(*,*) 'CHECK:', s_tmp% nz, n
         
         ! *** Now, first create a normal ZAMS star ***
         total_mass = sum(d_mass)
         ierr = new_particle(tmp2_id_new_model, total_mass)
         if (failed('new_particle', ierr)) return
         id_new_model = tmp2_id_new_model
         if (debugging) write(*,*) 'id_new_model', id_new_model
         call get_star_ptr(id_new_model, s, ierr)
         if (failed('get_star_ptr', ierr)) return
         
         ! *** Match the mesh containing the target structure to the mesh of the new particle ***
         ierr = match_mesh(tmp1_id_new_model, s% nz, s% dq)
         if (failed('match_mesh', ierr)) return
         
         ! *** Copy the relevant variables (chemical fractions only, or also hydro vars...)
         original_timestep_limit = s% min_timestep_limit
         s% min_timestep_limit = 1.0d-12
         original_dxdt_nuc_f = s% dxdt_nuc_factor
         s% dxdt_nuc_factor = 1.0d-99
         original_timestep = s% dt_next
         f = 1.0d-4
         if (debugging) then
            s% trace_evolve = .true.
            s% report_ierr = .true.
         end if
         do
            s% xa(:,:) = f * s_tmp% xa(:,:) + (1.0d0 - f) * s% xa(:,:)
            ierr = erase_memory(id_new_model)
            if (failed('erase_memory', ierr)) return
            s% dt_next = 10.0 * s% min_timestep_limit
            ierr = evolve_one_step(id_new_model)
            if (failed('evolve_one_step', ierr)) return
            if (debugging) write(*,*) 'f: ', f
            call check_remeshed(s% nz, s_tmp% nz, s% dq, s_tmp% dq, ierr)
            if (failed('check_remeshed', ierr)) then
               ierr = match_mesh(tmp1_id_new_model, s% nz, s% dq)
               if (failed('match_mesh', ierr)) return
               call check_remeshed(s% nz, s_tmp% nz, s% dq, s_tmp% dq, ierr)
               if (failed('check_remeshed 2', ierr)) return
            end if
            if (debugging) write(*,*) 'CHECK check_remeshed OK'
            if (debugging) write(*,*) 'Backups', s% number_of_backups_in_a_row
            if (s% number_of_backups_in_a_row > 15) exit
            if (f >= 1.0d0) exit
            if (f >= 0.1d0) then
               f = min(1.1d0 * f, 1.0d0)
            else
               f = 1.5d0 * f
            endif
         end do
         
         ! *** Give the model the opportunity to remesh ***
         s% mesh_delta_coeff = 0.5
         ierr = remesh(s, .true., .false., .false.) 
         if (failed('remesh', ierr)) return 
         ierr = erase_memory(id_new_model)
         if (failed('erase_memory', ierr)) return
         ierr = match_mesh(tmp1_id_new_model, s% nz, s% dq)
         if (failed('match_mesh', ierr)) return
         s% number_of_backups_in_a_row = 0
         s% mesh_delta_coeff = 1
         
         ! *** Optionally, also do hydro vars ***
         if (do_T) then
            f = 1.0d-8
            index_low = s% nz / 10 ! Do not meddle with the atmosphere!
            do
               s% xa(:,:) = s_tmp% xa(:,:)
               s% xs(s% i_lnT,index_low:) = f*s_tmp%xs(s_tmp%i_lnT,index_low:) + (1d0-f)*s%xs(s%i_lnT,index_low:)
               ierr = erase_memory(id_new_model)
               if (failed('erase_memory', ierr)) return
               s% dt_next = 10.0 * s% min_timestep_limit
               ierr = evolve_one_step(id_new_model)
               if (failed('evolve_one_step', ierr)) return
               if (debugging) write(*,*) 'f: ', f
               call check_remeshed(s% nz, s_tmp% nz, s% dq, s_tmp% dq, ierr)
               if (failed('check_remeshed', ierr)) return
               if (debugging) write(*,*) 'CHECK check_remeshed OK'
               if (debugging) write(*,*) 'Backups', s% number_of_backups_in_a_row
               if (f >= 1.0d0) exit
               f = min(1.5d0 * f, 1.0d0)
            end do
         end if
         
         ! *** Restore the original timestep ***
         if (debugging) write(*,*) 'timesteps', s% dt_old, s% dt, s% dt_next, original_timestep
         s% dt_next = 10.0 * s% min_timestep_limit
         s% dt = 10.0 * s% min_timestep_limit
         if (debugging) write(*,*) 'timesteps', s% dt_old, s% dt, s% dt_next, original_timestep
         ierr = evolve_one_step(id_new_model)
         if (debugging) write(*,*) ierr, s% result_reason, trim(result_reason_str(s% result_reason))
         if (debugging) write(*,*) 'timesteps', s% dt_old, s% dt, s% dt_next, original_timestep
         if (do_restore_timestep) then
            do k1 = 1, 10
               if (debugging) write(*,*) 'increasing timesteps', s% dt_old, s% dt, s% dt_next, original_timestep
               if (debugging) write(*,*) 'Backups', s% number_of_backups_in_a_row
               s% xa(:,:) = s_tmp% xa(:,:)
               if (do_T) s% xs(s% i_lnT,index_low:) = s_tmp% xs(s_tmp% i_lnT,index_low:)
               ierr = erase_memory(id_new_model)
               if (failed('erase_memory', ierr)) return
               do k2 = 1, 10
                  ierr = evolve_one_step(id_new_model)
                  if (debugging) write(*,*) ierr, s% result_reason, trim(result_reason_str(s% result_reason))
               end do
               if (s% number_of_backups_in_a_row > 0) exit
            end do
         end if
         
         call check_remeshed(s% nz, s_tmp% nz, s% dq, s_tmp% dq, ierr)
         if (failed('check_remeshed', ierr)) then
            ierr = match_mesh(tmp1_id_new_model, s% nz, s% dq)
            if (failed('match_mesh', ierr)) return
            call check_remeshed(s% nz, s_tmp% nz, s% dq, s_tmp% dq, ierr)
            if (failed('check_remeshed 2', ierr)) return
         end if
         if (do_T) s% xs(s% i_lnT,index_low:) = s_tmp% xs(s_tmp% i_lnT,index_low:)
         ierr = erase_memory(id_new_model)
         if (failed('erase_memory', ierr)) return
         
         s% dxdt_nuc_factor = original_dxdt_nuc_f
         
         if (s% dt_next > 10.0 * original_timestep_limit) then
            s% min_timestep_limit = original_timestep_limit
         else
            s% min_timestep_limit = s% dt_next / 10.0
         endif
         
         if (debugging) write(*,*) 'Backups:', s% number_of_backups_in_a_row
         s% number_of_backups_in_a_row = 0
         if (debugging) write(*,*) 'Backups reset:', s% number_of_backups_in_a_row
         s% trace_evolve = .false.
         s% report_ierr = .false.
         call flush()
         new_model_defined = .true.
         new_stellar_model = 0
      end if
      
      contains
      
      subroutine check_remeshed(nz, nz_orig, dq, dq_orig, ierr)
         implicit none
         integer, intent(in) :: nz, nz_orig
         double precision, intent(in) :: dq(nz), dq_orig(nz_orig)
         integer, intent(out) :: ierr
         integer :: i
         if (nz .ne. nz_orig) then
            ierr = -1
            return
         end if
         do i = 1, nz
            if (dq(i) .ne. dq_orig(i)) then
               ierr = -1
               return
            end if
         end do
         ierr = 0
      end subroutine check_remeshed
      
   end function new_stellar_model

   function finalize_stellar_model(star_id, age_tag)
      use amuse_support
      use evolve, only: set_age
      implicit none
      integer :: finalize_stellar_model, ierr
      integer, intent(out) :: star_id
      double precision, intent(in) :: age_tag
      
      if (.not. new_model_defined) then
         finalize_stellar_model = -35
         return
      endif
      
      finalize_stellar_model = -1
      star_id = id_new_model
      number_of_particles = star_id
      call set_age(id_new_model, age_tag, ierr)
      if (failed('set_age', ierr)) return
      call flush()
      
      new_model_defined = .false.
      finalize_stellar_model = 0
   end function

   ! matches/interpolates existing mesh based on supplied dq's
   integer function match_mesh(model_id, nz_target, dq_target)
      use amuse_support, only: failed, debugging
      use star_private_def, only: star_info, get_star_ptr
      use alloc, only: free_star_info_arrays, allocate_star_info_arrays
      use mesh_plan, only: do_mesh_plan
      use mesh_adjust, only: do_mesh_adjust
      use adjust_mesh_support, only: check_validity
      use hydro_vars, only: set_vars
      use rates_def, only: i_rate, ipp, icno, i3alf, iphoto
      use net_lib, only: clean_up_fractions
      use num_lib, only: safe_log10
      use utils_lib
      use star_utils, only: set_q_vars, report_xa_bad_nums, &
         std_dump_model_info_for_ndiff, set_qs, set_xqs
      use chem_def
      
      integer, intent(in) :: model_id, nz_target
      double precision, intent(inout) :: dq_target(nz_target)
      
      type (star_info), pointer :: s_tmp
      logical, parameter :: dbg_remesh = .true.
      logical, parameter :: skip_net = .false., check_for_bad_nums = .true.
      integer :: k, k2, ierr, species, nvar, nz, nz_new, nz_old, &
         unchanged, split, merged
      type (star_info), target :: prev_info
      type (star_info), pointer :: prv
      double precision, pointer, dimension(:) :: xq_old, xq_new, energy

      double precision, parameter :: max_sum_abs = 10d0
      double precision, parameter :: xsum_tol = 1d-2
      double precision, parameter :: h_cntr_limit = 0.5d0 ! for pre-MS decision
      double precision, parameter :: he_cntr_limit = 0.1d0 ! for RGB vs AGB decision
      
 3       format(a40,2i6,99(1pe26.16))
      
      call get_star_ptr(model_id, s_tmp, ierr)
      if (failed('get_star_ptr', ierr)) return
      if (debugging) write(*,*) 'enter match_mesh'
      ierr = 0
      match_mesh = -1
      
      species = s_tmp% species
      nz_old = s_tmp% nz
      nz = nz_old
      nz_new = nz_target
      
      call clean_up_fractions(1, nz, species, nz, s_tmp% xa, max_sum_abs, xsum_tol, ierr)
      if (failed('clean_up_fractions', ierr)) return
      
      nullify(xq_old, xq_new)
      allocate(energy(nz), stat=ierr)
      
      energy(1:nz) = exp(s_tmp% lnE(1:nz))
      
      s_tmp% mesh_call_number = s_tmp% mesh_call_number + 1
      
      ! save pointers to arrays that will need to be updated for new mesh
      prv => prev_info
      prv = s_tmp ! this makes copies of pointers and scalars
      
      if (associated(s_tmp% comes_from)) deallocate(s_tmp% comes_from)
      allocate(s_tmp% comes_from(nz_target), xq_old(nz), xq_new(nz_target), stat=ierr)
      if (failed('allocate', ierr)) return
      
      call check_validity(s_tmp, ierr)
      if (failed('check_validity', ierr)) return
      
      if (check_for_bad_nums) then
         if (has_bad_num(species*nz, s_tmp% xa)) then
            write(*,*) 'bad num in xa before calling mesh_plan: model_number', s_tmp% model_number
            call report_xa_bad_nums(s_tmp, ierr)
            stop 'remesh'
         end if
      end if
      
      call set_xqs(nz, xq_old, s_tmp% dq, ierr)
      if (failed('set_xqs xq_old', ierr)) return
      call set_xqs(nz_target, xq_new, dq_target, ierr)
      if (failed('set_xqs xq_new', ierr)) return
      
      ! Set comes_from
      !      ! xq_old(comes_from(k)+1) > xq_new(k) >= xq_old(comes_from(k)), if comes_from(k) < nz_old.
      s_tmp% comes_from(:) = 0
      k2 = 1
      s_tmp% comes_from(1) = k2
      do k = 2, nz_target
         do
            if (k2 == nz) exit
            if (xq_new(k) >= xq_old(k2+1)) then
               k2 = k2 + 1
            else
               exit
            end if
         end do
         s_tmp% comes_from(k) = k2
      end do
      nz = nz_new
      s_tmp% nz = nz
      nvar = s_tmp% nvar
      
      call allocate_star_info_arrays(s_tmp, ierr)
      if (failed('allocate_star_info_arrays', ierr)) return
      
      if (associated(s_tmp% cell_type)) deallocate(s_tmp% cell_type)
      allocate(s_tmp% cell_type(nz))
      call set_types_of_new_cells(s_tmp% cell_type)
      
      s_tmp% rate_factors(1:prv% num_reactions) = prv% rate_factors(1:prv% num_reactions)
      
      ! store new q and dq
      s_tmp% dq(:) = dq_target(:)
      call set_qs(nz, s_tmp% q, s_tmp% dq, ierr)
      if (failed('set_qs', ierr)) return
      
      ! testing -- check for q strictly decreasing
      do k = 2, nz
         if (xq_new(k) <= xq_new(k-1)) then
            write(*,3) 'bad xq_new before call do_mesh_adjust', &
               k, nz, xq_new(k), xq_new(k-1), dq_target(k-1), xq_new(k-1) + dq_target(k-1)
            stop 'adjust mesh'
         end if
      end do

      if (s_tmp% q_flag) call set_q_vars(s_tmp)
      
      if (dbg_remesh) write(*,*) 'call do_mesh_adjust'
      call do_mesh_adjust( &
         nz, nz_old, prv% xs, prv% xa, energy, prv% eta, prv% dq, xq_old, &
         s_tmp% xs, s_tmp% xa, s_tmp% dq, xq_new, s_tmp% species, s_tmp% chem_id, s_tmp% net_iso, s_tmp% eos_handle, &
         s_tmp% mesh_adjust_use_quadratic, s_tmp% mesh_adjust_get_T_from_E, &
         s_tmp% i_lnd, s_tmp% i_lnT, s_tmp% i_lnR, s_tmp% i_lum, s_tmp% i_vel, s_tmp% i_lndq, s_tmp% i_lnq, &
         s_tmp% q_flag, s_tmp% v_flag, &
         prv% mstar, s_tmp% comes_from, s_tmp% cell_type, ierr)
      if (failed('do_mesh_adjust', ierr)) return
      if (dbg_remesh) write(*,*) 'back from do_mesh_adjust'
      
      ! testing
      do k = 2, nz
         if (xq_new(k) <= xq_new(k-1)) then
            write(*,3) 'bad xq_new after call do_mesh_adjust', k, nz, xq_new(k), xq_new(k-1)
            stop 'adjust mesh'
         end if
      end do

      if (ierr /= 0 .and. s_tmp% report_ierr) then
         write(*,*) 'mesh_adjust problem'
         write(*,*) 'doing mesh_call_number', s_tmp% mesh_call_number
         write(*,*) 's_tmp% model_number', s_tmp% model_number
         write(*,*) 's_tmp% nz', s_tmp% nz
         write(*,*) 's_tmp% num_retries', s_tmp% num_retries
         write(*,*) 's_tmp% num_backups', s_tmp% num_backups
         write(*,*) 
      end if
      
      if (check_for_bad_nums) then
         if (has_bad_num(species*nz, s_tmp% xa)) then
            write(*,*) 'bad num in xa after calling mesh_adjust: model_number', s_tmp% model_number
            stop 'remesh'
         end if
      end if
      
      if (s_tmp% prev_cdc_tau > 0) then ! interpolate cdc
         call set_prev_cdc(ierr)
         if (ierr /= 0 .and. s_tmp% report_ierr) &
            write(*,*) 'mesh_adjust problem: ierr from set_prev_cdc'
      end if

      call free_star_info_arrays(prv)

      call dealloc
      match_mesh = 0

      contains
      
      subroutine set_prev_cdc(ierr)
         use interp_1d_def
         use interp_1d_lib
         integer, intent(out) :: ierr
         integer, parameter :: nwork = pm_work_size
         double precision, pointer :: work(:,:)
         ierr = 0
         allocate(work(nz_old, nwork), stat=ierr)
         if (ierr /= 0) return
         call interpolate_vector( &
            nz_old, prv% q, nz, s_tmp% q, prv% cdc, s_tmp% cdc_prev, interp_pm, nwork, work, ierr)
         deallocate(work)
      end subroutine set_prev_cdc

      
      subroutine set_types_of_new_cells(cell_type)
         use mesh_adjust, only: split_type, unchanged_type, merged_type
         integer, pointer :: cell_type(:)
         integer :: k, k_old, new_type
         
 2       format(a40,2i6,99(1pe26.16))
         
         unchanged=0; split=0; merged=0
      
         do k=1,nz_new
            k_old = s_tmp% comes_from(k)
            new_type = -111
            if (xq_new(k) < xq_old(k_old)) then
               write(*,*) 'xq_new(k) < xq_old(k_old)'
               write(*,2) 'xq_new(k)', k, xq_new(k)
               write(*,2) 'xq_old(k_old)', k_old, xq_old(k_old)
               write(*,*) 'adjust mesh set_types_of_new_cells'
               stop 1
            else if (xq_new(k) > xq_old(k_old)) then
               new_type = split_type
            else if (k_old == nz_old) then
               if (k == nz_new) then
                  new_type = unchanged_type
               else
                  new_type = split_type
               end if
            else if (k == nz_new) then
               new_type = split_type
            else ! k_old < nz_old .and. k < nz .and. xq_new(k) == xq_old(k_old)
               if (xq_new(k+1) == xq_old(k_old+1)) then
                  new_type = unchanged_type
               else if (xq_new(k+1) > xq_old(k_old+1)) then
                  new_type = merged_type
               else
                  new_type = split_type
               end if
            end if
            cell_type(k) = new_type
            select case (new_type)
               case (split_type)
                  split = split + 1
               case (unchanged_type)
                  unchanged = unchanged + 1
               case (merged_type)
                  merged = merged + 1
               case default
                  write(*,*) 'failed to set new_type in adjust mesh set_types_of_new_cells'
                  stop 'set_types_of_new_cells'
            end select
         end do
         
         if (unchanged + split + merged /= nz_new) then
            write(*,2) 'unchanged + split + merged', unchanged + split + merged
            write(*,2) 'nz_new', nz_new
            stop 'set_types_of_new_cells'
         end if
      
      end subroutine set_types_of_new_cells

      subroutine dealloc
         if (associated(xq_old)) deallocate(xq_old)
         if (associated(xq_new)) deallocate(xq_new)
         if (associated(energy)) deallocate(energy)
      end subroutine dealloc

   end function match_mesh
