module amuse_support
   implicit none
   character (len=4096) :: AMUSE_inlist_path
   character (len=4096) :: AMUSE_mesa_dir
   character (len=4096) :: AMUSE_local_data_dir ! Used for output starting_models
   character (len=4096) :: AMUSE_gyre_in_file 
   double precision :: AMUSE_mass
   double precision :: AMUSE_metallicity = 0.02d0
   double precision :: AMUSE_dmass = 0.1d0
   double precision :: AMUSE_mlo = -1.0d0
   double precision :: AMUSE_mhi = 2.0d0
   double precision :: AMUSE_max_age_stop_condition = 1.0d36
   double precision :: AMUSE_min_timestep_stop_condition = 1.0d-6
   integer :: AMUSE_max_iter_stop_condition = -1111

   double precision, allocatable :: target_times(:)
   integer :: number_of_particles ! Dead or alive...

   logical :: new_model_defined = .false.
   integer :: id_new_model

   logical :: debugging = .false.

end module amuse_support


module amuse_mesa
   use mesa_interface
   use amuse_support

   implicit none

   contains

! Set the paths to the inlist and the data directory
   integer function set_MESA_paths(AMUSE_inlist_path_in, &
         AMUSE_mesa_dir_in, AMUSE_local_data_dir_in, AMUSE_gyre_in_file_in)
      
      character(*), intent(in) :: AMUSE_inlist_path_in, &
         AMUSE_mesa_dir_in, AMUSE_local_data_dir_in,  AMUSE_gyre_in_file_in

      AMUSE_inlist_path = AMUSE_inlist_path_in
      AMUSE_mesa_dir = AMUSE_mesa_dir_in
      AMUSE_local_data_dir = AMUSE_local_data_dir_in
      AMUSE_gyre_in_file =  AMUSE_gyre_in_file_in

      set_MESA_paths = 0
   end function set_MESA_paths

! Initialize the stellar evolution code
   integer function initialize_code()
      
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

! Create a new particle
   integer function new_particle(AMUSE_id, AMUSE_value)
      integer, intent(out) :: AMUSE_id
      double precision, intent(in) :: AMUSE_value
      integer :: id, ierr
      
      new_particle = -1
      call allocate_star(id, ierr)
      if(ierr/=MESA_SUCESS) return

      AMUSE_mass = AMUSE_value

      AMUSE_id = id
      number_of_particles = AMUSE_id
      call flush()
      new_particle = 0
   end function new_particle

! load zams model 
   integer function new_zams_particle(AMUSE_id, AMUSE_value)
      integer, intent(out) :: AMUSE_id
      double precision, intent(in) :: AMUSE_value
      integer :: ierr

      new_zams_particle = -1

      new_zams_particle = new_particle(AMUSE_id, AMUSE_value)

      if(new_zams_particle/=0 ) return  

      call load_zams_model(AMUSE_id, ierr)
      if(ierr/=MESA_SUCESS) return

      call finish_init_star(AMUSE_id, set_amuse_options, ierr)
      if(ierr/=MESA_SUCESS) return

      new_zams_particle = 0

   end function new_zams_particle

! load an existing stellar model 
   integer function load_model(AMUSE_id, AMUSE_filename)
      integer, intent(out) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_filename
      integer :: ierr

      load_model = -1

      call allocate_star(AMUSE_id, ierr) 
      if(ierr/=MESA_SUCESS) return

      ierr = MESA_FAIL
      call load_saved_model(AMUSE_id, AMUSE_filename, ierr)
      if(ierr/=MESA_SUCESS) return

      call finish_init_star(AMUSE_id, set_amuse_options, ierr)
      if(ierr/=MESA_SUCESS) return

      load_model = 0

   end function load_model

! Create a new pre-main-sequence star
   integer function new_pre_ms_particle(AMUSE_id, AMUSE_value)
      integer, intent(out) :: AMUSE_id
      double precision :: AMUSE_value
      integer ::  ierr

      new_pre_ms_particle = -1

      new_pre_ms_particle = new_particle(AMUSE_id, AMUSE_value)

      if(new_pre_ms_particle/=0 ) return  

      ierr = MESA_FAIL

      call create_pre_main_sequence(AMUSE_id, ierr)
      if(ierr/=MESA_SUCESS) return

      call finish_init_star(AMUSE_id, set_amuse_options, ierr)
      if(ierr/=MESA_SUCESS) return

      new_pre_ms_particle = 0

   end function new_pre_ms_particle


! Create a new pure He star
   ! Must set initial_z before calling this
   integer function new_pure_he_particle(AMUSE_id, AMUSE_value)
      integer, intent(out) :: AMUSE_id
      double precision :: AMUSE_value
      integer ::  ierr

      new_pure_he_particle = -1

      new_pure_he_particle = new_particle(AMUSE_id, AMUSE_value)

      if(new_pure_he_particle/=0 ) return  

      ierr = MESA_FAIL

      call create_he_star(AMUSE_id, ierr)
      if(ierr/=MESA_SUCESS) return

      call finish_init_star(AMUSE_id, set_amuse_options, ierr)
      if(ierr/=MESA_SUCESS) return

      new_pure_he_particle = 0

   end function new_pure_he_particle

   subroutine set_amuse_options(AMUSE_id)
      ! Dont call this directly as variables will be reset during initlization
      ! Instead it must be used as a calback function in finish_init_star()
      integer, intent(in) :: AMUSE_id

      call set_init_options(AMUSE_id, &
                           AMUSE_mass, &
                           AMUSE_mesa_dir, &
                           AMUSE_local_data_dir, &
                           AMUSE_inlist_path, &
                           AMUSE_metallicity,&
                           AMUSE_max_age_stop_condition, &
                           AMUSE_min_timestep_stop_condition, &
                           AMUSE_max_iter_stop_condition, &
                           AMUSE_gyre_in_file &
                           )

   end subroutine set_amuse_options


! Remove a particle 
   integer function delete_star(AMUSE_id)
      integer, intent(in) :: AMUSE_id
      integer :: ierr
      delete_star = -1

      call remove_star(AMUSE_id, ierr)
      if(ierr/=MESA_SUCESS) return

      delete_star = 0

   end function delete_star

   integer function commit_particles()
      allocate(target_times(number_of_particles))
      target_times = 0
      commit_particles = 0
   end function commit_particles

   integer function recommit_particles()
      double precision, allocatable :: temp(:)
      allocate(temp(size(target_times)))
      temp = target_times
      deallocate(target_times)
      allocate(target_times(number_of_particles))
      target_times = 0
      target_times(1:size(temp)) = temp
      deallocate(temp)
      recommit_particles = 0
   end function recommit_particles

! Get/setters for code parameters:

! Return the number of particles currently allocated in the code
   integer function get_number_of_particles(AMUSE_value)
      integer, intent(out) :: AMUSE_value
      AMUSE_value = 1
      get_number_of_particles = how_many_allocated_star_ids()
   end function get_number_of_particles

! Return the metallicity parameter
   integer function get_metallicity(AMUSE_value)
      double precision, intent(out) :: AMUSE_value
      AMUSE_value = AMUSE_metallicity
      get_metallicity = 0
   end function get_metallicity

! Set the metallicity parameter
   integer function set_metallicity(AMUSE_value)
      double precision, intent(in) :: AMUSE_value
      integer :: ierr
      AMUSE_metallicity = AMUSE_value
      set_metallicity = 0
   end function set_metallicity

! Return the current mass of the star
   integer function get_mass(AMUSE_id, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      double precision, intent(out) :: AMUSE_value
      integer :: ierr

      get_mass = 0
      call get_history_value(AMUSE_id,'star_mass', AMUSE_value, ierr)

      if (ierr /= MESA_SUCESS) then
         AMUSE_value = -1.0
         get_mass = -1
      endif

   end function get_mass


! Return the current core mass of the star, where hydrogen abundance is <= h1_boundary_limit
   integer function get_core_mass(AMUSE_id, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      double precision, intent(out) :: AMUSE_value
      integer :: ierr

      get_core_mass = 0
      call get_history_value(AMUSE_id,'he_core_mass', AMUSE_value, ierr)

      if (ierr /= MESA_SUCESS) then
         AMUSE_value = -1.0
         get_core_mass = -1
      endif
   end function get_core_mass

! Return the current mass loss rate of the star
   integer function get_mass_loss_rate(AMUSE_id, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      double precision, intent(out) :: AMUSE_value
      integer :: ierr
      
      get_mass_loss_rate = 0
      call get_history_value(AMUSE_id,'star_mdot', AMUSE_value, ierr)

      if (ierr /= MESA_SUCESS) then
         AMUSE_value = -1.0
         get_mass_loss_rate = -1
      endif

   end function get_mass_loss_rate

! Return the current user-specified mass transfer rate of the star
   integer function get_manual_mass_transfer_rate(AMUSE_id, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      double precision, intent(out) :: AMUSE_value
      integer :: ierr

      get_manual_mass_transfer_rate = 0
      call get_control_opt_dble(AMUSE_id, 'mass_change', AMUSE_value, ierr)

      if (ierr /= MESA_SUCESS) then
         AMUSE_value = -1.0
         get_manual_mass_transfer_rate = -1
      endif
   end function get_manual_mass_transfer_rate

! Set a new user-specified mass transfer rate of the star
   integer function set_manual_mass_transfer_rate(AMUSE_id, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      double precision, intent(out) :: AMUSE_value
      integer :: ierr
      set_manual_mass_transfer_rate = 0
      
      call set_control_opt_dble(AMUSE_id, 'mass_change', AMUSE_value, ierr)

      if (ierr /= MESA_SUCESS) then
         AMUSE_value = -1.0
         set_manual_mass_transfer_rate = -1
      endif
   end function set_manual_mass_transfer_rate

! Return the current temperature of the star
   integer function get_temperature(AMUSE_id, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      double precision, intent(out) :: AMUSE_value
      integer ::  ierr

      get_temperature = 0
      call get_history_value(AMUSE_id,'effective_T', AMUSE_value, ierr)

      if (ierr /= MESA_SUCESS) then
         AMUSE_value = -1.0
         get_temperature = -1
      endif

   end function get_temperature

! Return the current luminosity of the star
   integer function get_luminosity(AMUSE_id, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      double precision, intent(out) :: AMUSE_value
      integer :: ierr
      get_luminosity = 0
      call get_history_value(AMUSE_id,'luminosity', AMUSE_value, ierr)

      if (ierr /= MESA_SUCESS) then
         AMUSE_value = -1.0
         get_luminosity = -1
      endif
   end function

! Return the current age of the star
   integer function get_age(AMUSE_id, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      double precision, intent(out) :: AMUSE_value
      integer :: ierr
      get_age = 0
      call get_history_value(AMUSE_id,'star_age', AMUSE_value, ierr)

      if (ierr /= MESA_SUCESS) then
         AMUSE_value = -1.0
         get_age = -1
      endif
   end function

! Return the next timestep for the star in years
   integer function get_time_step(AMUSE_id, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      double precision, intent(out) :: AMUSE_value
      integer :: ierr
      get_time_step = 0
      call get_next_timestep(AMUSE_id, AMUSE_value, ierr)
      if (ierr /= MESA_SUCESS) then
         AMUSE_value = -1.0
         get_time_step = -1
      endif
   end function

! Set the next timestep for the star in years
   integer function set_time_step(AMUSE_id, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      double precision, intent(in) :: AMUSE_value
      integer :: ierr
      set_time_step = 0
      call set_timestep(AMUSE_id, AMUSE_value, ierr)
      if (ierr /= MESA_SUCESS) then
         set_time_step = -1
      endif
   end function

! Given mesa controls option AMUSE_name set it to AMUSE_value Where AMUSE_value is a double
   integer function set_control_dble(AMUSE_id, AMUSE_name, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_name
      double precision, intent(in) :: AMUSE_value
      integer ::  ierr
      set_control_dble = 0
      call set_control_opt_dble(AMUSE_id, AMUSE_name, AMUSE_value, ierr)
      if (ierr /= MESA_SUCESS) then
         set_control_dble = -1
      endif

   end function set_control_dble

! Given mesa controls option AMUSE_name set it to AMUSE_value Where AMUSE_value is a integer
   integer function set_control_int(AMUSE_id, AMUSE_name, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_name
      integer, intent(in) :: AMUSE_value
      integer ::  ierr
      set_control_int = 0
      call set_control_opt_int(AMUSE_id, AMUSE_name, AMUSE_value, ierr)
      if (ierr /= MESA_SUCESS) then
         set_control_int = -1
      endif

   end function set_control_int

! Given mesa controls option AMUSE_name set it to AMUSE_value Where AMUSE_value is a logical flag
   integer function set_control_logical(AMUSE_id, AMUSE_name, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_name
      logical, intent(in) :: AMUSE_value
      integer ::  ierr
      set_control_logical = 0
      call set_control_opt_logical(AMUSE_id, AMUSE_name, AMUSE_value, ierr)
      if (ierr /= MESA_SUCESS) then
         set_control_logical = -1
      endif

   end function set_control_logical


   ! Setting options in star_job does not actual change anything, this function must be called
   ! to actually update the model
   integer function star_job_update(AMUSE_id)
      integer, intent(in) :: AMUSE_id
      integer ::  ierr
      star_job_update = 0
      call update_star_job(AMUSE_id, ierr)
      if (ierr /= MESA_SUCESS) then
         star_job_update = -1
      endif

   end function star_job_update


! Given mesa controls option AMUSE_name set it to AMUSE_value Where AMUSE_value is a string
   integer function set_control_str(AMUSE_id,AMUSE_name,  AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_name
      character(len=*), intent(in) :: AMUSE_value
      integer ::  ierr
      set_control_str = 0
      call set_control_opt_str(AMUSE_id, AMUSE_name, AMUSE_value, ierr)
      if (ierr /= MESA_SUCESS) then
         set_control_str = -1
      endif

   end function set_control_str

! Given mesa star_job option AMUSE_name set it to AMUSE_value Where AMUSE_value is a double
   integer function set_star_job_dble(AMUSE_id, AMUSE_name, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_name
      double precision, intent(in) :: AMUSE_value
      integer ::  ierr
      set_star_job_dble = 0
      call set_star_job_opt_dble(AMUSE_id, AMUSE_name, AMUSE_value, ierr)
      if (ierr /= MESA_SUCESS) then
         set_star_job_dble = -1
      endif

   end function set_star_job_dble

! Given mesa star_job option AMUSE_name set it to AMUSE_value Where AMUSE_value is a integer
   integer function set_star_job_int(AMUSE_id, AMUSE_name, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_name
      integer, intent(in) :: AMUSE_value
      integer ::  ierr
      set_star_job_int = 0
      call set_star_job_opt_int(AMUSE_id, AMUSE_name, AMUSE_value, ierr)
      if (ierr /= MESA_SUCESS) then
         set_star_job_int = -1
      endif

   end function set_star_job_int

! Given mesa star_job option AMUSE_name set it to AMUSE_value Where AMUSE_value is a logical flag
   integer function set_star_job_logical(AMUSE_id, AMUSE_name, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_name
      logical, intent(in) :: AMUSE_value
      integer ::  ierr
      set_star_job_logical = 0
      call set_star_job_opt_logical(AMUSE_id, AMUSE_name, AMUSE_value, ierr)
      if (ierr /= MESA_SUCESS) then
         set_star_job_logical = -1
      endif

   end function set_star_job_logical

! Given mesa star_job option AMUSE_name set it to AMUSE_value Where AMUSE_value is a string
   integer function set_star_job_str(AMUSE_id, AMUSE_name, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_name
      character(len=*), intent(in) :: AMUSE_value
      integer ::  ierr
      set_star_job_str = 0
      call set_star_job_opt_str(AMUSE_id, AMUSE_name, AMUSE_value, ierr)
      if (ierr /= MESA_SUCESS) then
         set_star_job_str = -1
      endif

   end function set_star_job_str 


! Given mesa controls option AMUSE_name put the value into AMUSE_value Where AMUSE_value is a double
   integer function get_control_dble(AMUSE_id, AMUSE_name, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_name
      double precision, intent(out) :: AMUSE_value
      integer ::  ierr
      get_control_dble = 0
      call get_control_opt_dble(AMUSE_id, AMUSE_name, AMUSE_value, ierr)
      if (ierr /= MESA_SUCESS) then
         get_control_dble = -1
      endif
   end function get_control_dble

! Given mesa controls option AMUSE_name put the value into AMUSE_value Where AMUSE_value is a integer
   integer function get_control_int(AMUSE_id, AMUSE_name, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_name
      integer, intent(out) :: AMUSE_value
      integer ::  ierr
      get_control_int = 0
      call get_control_opt_int(AMUSE_id, AMUSE_name, AMUSE_value, ierr)
      if (ierr /= MESA_SUCESS) then
         get_control_int = -1
      endif
   end function get_control_int

! Given mesa controls option AMUSE_name put the value into AMUSE_value Where AMUSE_value is a logical flag
   integer function get_control_logical(AMUSE_id, AMUSE_name, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_name
      logical, intent(out) :: AMUSE_value
      integer ::  ierr
      get_control_logical = 0
      call get_control_opt_logical(AMUSE_id, AMUSE_name, AMUSE_value, ierr)

      if (ierr /= MESA_SUCESS) then
         get_control_logical = -1
      endif
   end function get_control_logical

! Given mesa controls option AMUSE_name put the value into AMUSE_value Where AMUSE_value is a string
   integer function get_control_str(AMUSE_id,AMUSE_name,  AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_name
      character(len=*), intent(out) :: AMUSE_value
      integer ::  ierr
      get_control_str = 0
      call get_control_opt_str(AMUSE_id, AMUSE_name, AMUSE_value, ierr)
      if (ierr /= MESA_SUCESS) then
         get_control_str = -1
      endif

   end function get_control_str

! Given mesa star_job option AMUSE_name put the value into AMUSE_value Where AMUSE_value is a double
   integer function get_star_job_dble(AMUSE_id, AMUSE_name, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_name
      double precision, intent(out) :: AMUSE_value
      integer ::  ierr
      get_star_job_dble = 0
      call get_star_job_opt_dble(AMUSE_id, AMUSE_name, AMUSE_value, ierr)
      if (ierr /= MESA_SUCESS) then
         get_star_job_dble = -1
      endif

   end function get_star_job_dble

! Given mesa star_job option AMUSE_name put the value into AMUSE_value Where AMUSE_value is a integer
   integer function get_star_job_int(AMUSE_id, AMUSE_name, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_name
      integer, intent(out) :: AMUSE_value
      integer ::  ierr
      get_star_job_int = 0
      call get_star_job_opt_int(AMUSE_id, AMUSE_name, AMUSE_value, ierr)
      if (ierr /= MESA_SUCESS) then
         get_star_job_int = -1
      endif

   end function get_star_job_int

! Given mesa star_job option AMUSE_name put the value into AMUSE_value Where AMUSE_value is a logical flag
   integer function get_star_job_logical(AMUSE_id, AMUSE_name, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_name
      logical, intent(out) :: AMUSE_value
      integer ::  ierr
      get_star_job_logical = 0
      call get_star_job_opt_logical(AMUSE_id, AMUSE_name, AMUSE_value, ierr)
      if (ierr /= MESA_SUCESS) then
         get_star_job_logical = -1
      endif

   end function get_star_job_logical

! Given mesa star_job option AMUSE_name put the value into AMUSE_value Where AMUSE_value is a string
   integer function get_star_job_str(AMUSE_id, AMUSE_name, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_name
      character(len=*), intent(out) :: AMUSE_value
      integer ::  ierr
      get_star_job_str = 0
      call get_star_job_opt_str(AMUSE_id, AMUSE_name, AMUSE_value, ierr)
      if (ierr /= MESA_SUCESS) then
         get_star_job_str = -1
      endif

   end function get_star_job_str 

! Return the current radius of the star
   integer function get_radius(AMUSE_id, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      double precision, intent(out) :: AMUSE_value
      integer ::  ierr
      get_radius = 0
      call get_history_value(AMUSE_id,'radius', AMUSE_value, ierr)

      if (ierr /= MESA_SUCESS) then
         AMUSE_value = -1.0
         get_radius = -1
      endif
   end function

! Return the current stellar type of the star
   integer function get_stellar_type(AMUSE_id, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      integer, intent(out) :: AMUSE_value
      integer :: ierr
      double precision :: lgR, mass, ch1, che3, che4, cc12, cne20
      double precision :: lgL, lgLH, ah1, ahe3, ahe4, lgLHe
      AMUSE_value = -99
      get_stellar_type = -1

      call get_history_value(AMUSE_id, 'log_R', lgR, ierr)
      if(ierr/=MESA_SUCESS) return
      call get_history_value(AMUSE_id,'star_mass',mass, ierr)
      if(ierr/=MESA_SUCESS) return

      call get_history_value(AMUSE_id,'center h1',ch1, ierr)
      if(ierr/=MESA_SUCESS) return
      call get_history_value(AMUSE_id,'center he3',che3, ierr)
      if(ierr/=MESA_SUCESS) return
      call get_history_value(AMUSE_id,'center he4',che4, ierr)
      if(ierr/=MESA_SUCESS) return      
      call get_history_value(AMUSE_id,'center c12',cc12, ierr)
      if(ierr/=MESA_SUCESS) return     
      call get_history_value(AMUSE_id,'center ne20',cne20, ierr)
      if(ierr/=MESA_SUCESS) return  

      call get_history_value(AMUSE_id,'log_LH',lgLH, ierr)
      if(ierr/=MESA_SUCESS) return     
      call get_history_value(AMUSE_id,'log_LHe',lgLHe, ierr)
      if(ierr/=MESA_SUCESS) return     
      call get_history_value(AMUSE_id,'log_L',lgL, ierr)
      if(ierr/=MESA_SUCESS) return  

      call get_history_value(AMUSE_id,'average h1',ah1, ierr)
      if(ierr/=MESA_SUCESS) return
      call get_history_value(AMUSE_id,'average he3',ahe3, ierr)
      if(ierr/=MESA_SUCESS) return
      call get_history_value(AMUSE_id,'average he4',ahe4, ierr)
      if(ierr/=MESA_SUCESS) return



      if (lgR > -1.0) then
         if(ch1 > 1d-4) then ! MS
            if(lgLH-lgL < -1 ) then ! Pre-main-sequence star
               AMUSE_value = 17
            else if(mass < 0.75) then
               AMUSE_value = 0 ! Convective low mass star
            else
               AMUSE_value = 1 ! Main sequence star
            endif                  
         else 
            if(lgLHe - lgL < -3) then
               AMUSE_value = 3 ! Red giant branch
            else
               if(che4 > 1d-4) then 
                  AMUSE_value = 4 ! Core He burning
                  if(ah1 > 1d-5) then
                     if(ahe3+ahe4 < 0.75*ah1) then
                        AMUSE_value = 5 ! Early AGB (inert C/O core)
                     else
                        AMUSE_value = 6 ! Late (thermally pulsing) AGB (inert C/O core)
                     end if
                  end if
               else
                  if (che3 + che4 > 1.0d-5) then
                     AMUSE_value = 7 ! Helium MS star
                  else
                     AMUSE_value = 9 ! Helium giant
                  end if
               end if
            end if
         end if
      else ! stellar remnant?
         if (mass <= 0.0075) then ! Planet
            AMUSE_value = 18
         else if (mass > 0.0075 .and. mass < 0.075) then ! brown dwarf
            AMUSE_value = 19
         else if (mass > 0.075 .and. mass < 1.44) then ! white dwarf
            ! Helium White Dwarf:
            if (che3 + che4 > 0.1) AMUSE_value = 10
            ! Carbon/Oxygen White Dwarf:
            if (cc12 > 0.1) AMUSE_value = 11
            ! Oxygen/Neon White Dwarf:
            if (cne20 > 0.01) AMUSE_value = 12
               ! Else? Unknown kind of white dwarf... hopefully never reached.
            if (AMUSE_value == -99) AMUSE_value = -10
         else
            if (mass < 3.2) then
               AMUSE_value = 13 ! Neutron Star
            else
               AMUSE_value = 14 ! Black Hole
            endif
         endif
      endif
      get_stellar_type = 0
   end function

! Return the current number of zones/mesh-cells of the star
   integer function get_number_of_zones(AMUSE_id, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      integer, intent(out) :: AMUSE_value
      double precision :: val
      integer :: ierr
      get_number_of_zones = 0
      call get_history_value(AMUSE_id,'num_zones', val, ierr)

      if (ierr /= MESA_SUCESS) then
         AMUSE_value = -1.0
         get_number_of_zones = -1
      endif
      AMUSE_value = int(val)
   end function

! Return the profile at the specified zone/mesh-cell of the star
   integer function get_profile_at_zone(AMUSE_id, AMUSE_zone, AMUSE_name, AMUSE_value)
      integer, intent(in) :: AMUSE_id, AMUSE_zone
      character(len=*), intent(in) ::AMUSE_name
      double precision, intent(out) :: AMUSE_value
      integer :: ierr, zone
      get_profile_at_zone = 0
      zone = reverse_zone_id(AMUSE_id, AMUSE_zone, ierr)
      call get_profile_value_zone(AMUSE_id, AMUSE_name, zone, AMUSE_value, ierr)

      if (ierr /= MESA_SUCESS) then
         AMUSE_value = -1.0
         get_profile_at_zone = -1
      endif
   end function

! Return the history column of the star
   integer function get_history(AMUSE_id, AMUSE_name, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) ::AMUSE_name
      double precision, intent(out) :: AMUSE_value
      integer :: ierr

      get_history = 0

      call get_history_value(AMUSE_id, AMUSE_name, AMUSE_value, ierr)

      if (ierr /= MESA_SUCESS) then
         AMUSE_value = -1.0
         get_history = -1
      endif
   end function

   ! Return the current number of chemical abundance variables per zone of the star
   integer function get_number_of_species(AMUSE_id, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      integer, intent(out) :: AMUSE_value
      double precision :: val
      integer :: ierr

      get_number_of_species = 0
      call get_history_value(AMUSE_id,'species', val, ierr)

      if (ierr /= MESA_SUCESS) then
         AMUSE_value = -1
         get_number_of_species = -1
      endif    
      AMUSE_value = int(val)

   end function

   ! Return the mass number of chemical abundance variable 'AMUSE_species' of the star
   integer function get_mass_of_species(AMUSE_id, AMUSE_species, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      integer,intent(in) :: AMUSE_species
      double precision, intent(out) :: AMUSE_value
      integer :: ierr
      get_mass_of_species = 0
      
      call get_mass_number_species(AMUSE_id, AMUSE_species, AMUSE_value, ierr)

      if (ierr /= MESA_SUCESS) then
         AMUSE_value = -1
         get_mass_of_species = ierr
      endif   
   end function

   ! Return the name of chemical abundance given by 'AMUSE_species' of the star
   integer function get_name_of_species(AMUSE_id, AMUSE_species, AMUSE_value)
      integer, intent(in) :: AMUSE_id, AMUSE_species
      character(len=*), intent(out) :: AMUSE_value

      integer :: ierr
      
      get_name_of_species = 0

      call get_species_name(AMUSE_id, AMUSE_species, AMUSE_value, ierr)

      if (ierr /= MESA_SUCESS) then
         AMUSE_value = ''
         get_name_of_species = ierr
      endif   

   end function
   

   ! Return the chem_id of the species given by AMUSE_species 
   integer function get_id_of_species(AMUSE_id, AMUSE_species, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_species
      integer, intent(out) :: AMUSE_value
      integer :: ierr
      
      get_id_of_species = 0

      call get_species_id(AMUSE_id, AMUSE_species, AMUSE_value, ierr)

      if (ierr /= MESA_SUCESS) then
         AMUSE_value = -1
         get_id_of_species= ierr
      endif   

   end function

! Retrieve the name of the current nuclear network
   integer function get_nuclear_network(AMUSE_id, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(out) :: AMUSE_value
      integer :: ierr

      get_nuclear_network = 0

      call get_net_name(AMUSE_id, AMUSE_value, ierr)

      if (ierr /= MESA_SUCESS) then
         AMUSE_value = ''
         get_nuclear_network= ierr
      endif   

   end function get_nuclear_network

! Set the name of the current nuclear network
   integer function set_nuclear_network(AMUSE_id, AMUSE_value)
      integer, intent(in) :: AMUSE_id
      character(len=*), intent(in) :: AMUSE_value
      integer :: ierr

      set_nuclear_network = 0

      call set_net_name(AMUSE_id, AMUSE_value, ierr)

      if (ierr /= MESA_SUCESS) then
         set_nuclear_network = ierr
      endif   

   end function set_nuclear_network



   ! Evolve the star for one step (for calls from amuse)
   function evolve_one_step(AMUSE_id)
      integer, intent(in) :: AMUSE_id
      integer :: ierr, evolve_one_step, res
      double precision :: result

      evolve_one_step = -1

      call do_evolve_one_step(AMUSE_id, res, ierr)

      if (ierr /= MESA_SUCESS) return   

      result = get_age(AMUSE_id, target_times(AMUSE_id))

      evolve_one_step = res

   end function evolve_one_step

   ! Evolve the star until AMUSE_end_time
   integer function evolve_for(AMUSE_id, AMUSE_delta_t)
      integer, intent(in) :: AMUSE_id
      double precision, intent(in) :: AMUSE_delta_t
      integer :: ierr
      double precision :: result

      evolve_for = 0
      call evolve_until(AMUSE_id, AMUSE_delta_t, ierr)
      if (ierr /= MESA_SUCESS) then
         evolve_for = -1
         call flush()
         return
      endif        

      result = get_age(AMUSE_id, target_times(AMUSE_id))
   end function evolve_for

! Return the maximum age stop condition
   integer function get_max_age_stop_condition(AMUSE_value)
      double precision, intent(out) :: AMUSE_value
      AMUSE_value = AMUSE_max_age_stop_condition
      get_max_age_stop_condition = 0
   end function get_max_age_stop_condition

! Set the maximum age stop condition
   integer function set_max_age_stop_condition(AMUSE_value)
      double precision, intent(in) :: AMUSE_value
      AMUSE_max_age_stop_condition = AMUSE_value
      set_max_age_stop_condition = 0
   end function set_max_age_stop_condition

! Return the maximum age stop condition
   integer function get_max_iter_stop_condition(AMUSE_value)
      integer, intent(out) :: AMUSE_value
      AMUSE_value = AMUSE_max_iter_stop_condition
      get_max_iter_stop_condition = 0
   end function get_max_iter_stop_condition

! Set the maximum age stop condition
   integer function set_max_iter_stop_condition(AMUSE_value)
      integer, intent(in) :: AMUSE_value
      AMUSE_max_iter_stop_condition = AMUSE_value
      set_max_iter_stop_condition = 0
   end function set_max_iter_stop_condition

! Return the minimum timestep stop condition
   integer function get_min_timestep_stop_condition(AMUSE_value)
      double precision, intent(out) :: AMUSE_value
      AMUSE_value = AMUSE_min_timestep_stop_condition
      get_min_timestep_stop_condition = 0
   end function get_min_timestep_stop_condition

! Set the minimum timestep stop condition
   integer function set_min_timestep_stop_condition(AMUSE_value)
      double precision, intent(in) :: AMUSE_value
      AMUSE_min_timestep_stop_condition = AMUSE_value
      set_min_timestep_stop_condition = 0
   end function set_min_timestep_stop_condition

! Retrieve the maximum number of stars that can be allocated in the code
   integer function get_maximum_number_of_stars(AMUSE_value)
      integer, intent(out) :: AMUSE_value
      call max_num_stars(AMUSE_value)
      get_maximum_number_of_stars = 0
   end function get_maximum_number_of_stars

! Create a new particle from a user supplied model (non-ZAMS, e.g. merger product)
   ! This expects all arrays to be in MESA order (1==surface n==center)
      ! 1-q (where q is mass fraction at zone i)
   ! star_mass needs to be an array for the interface to work but we only need the total star mass in star_mass(1)
   
   ! Negative return value indicates an error while positive is the new id
   integer function new_stellar_model(star_mass, xq, rho, temperature, &
         XH1,XHe3,XHe4,XC12,XN14,XO16,XNe20,XMg24,XSi28,XS32, &
         XAr36,XCa40,XTi44,XCr48,XFe52,XFe54,XFe56,XCr56,XNi56,&
         prot,neut, n)
      integer, intent(in) :: n
      double precision,dimension(n), intent(in) :: star_mass, xq,&
         XH1,XHe3,XHe4,XC12,XN14,XO16,XNe20,XMg24, XSi28,XS32, &
         XAr36,XCa40,XTi44,XCr48,XFe52,XFe54,XFe56,XCr56,XNi56,&
         prot, neut,&
         rho, temperature
      integer :: id, ierr
      double precision,allocatable :: xa(:,:)

      new_stellar_model = -1
      new_model_defined = .false.

      ! Load ZAMS to start with
      new_stellar_model = new_zams_particle(id, star_mass(1))
      if(new_stellar_model/=0) then
         new_stellar_model = -10
         return
      end if
      ! Set network
      new_stellar_model = set_nuclear_network(id, 'approx21.net')
      if(new_stellar_model/=0) then
         new_stellar_model = -11
         return
      end if
      allocate(xa(21,n))

      ! Build composition array
      call set_comp('h1',XH1)
      if(new_stellar_model/=0) then
         new_stellar_model = -12
         return
      end if

      call set_comp('he3',XHe3)
      if(new_stellar_model/=0) then
         new_stellar_model = -13
         return
      end if

      call set_comp('he4',XHe4)
      if(new_stellar_model/=0) then
         new_stellar_model = -14
         return
      end if

      call set_comp('c12',XC12)
      if(new_stellar_model/=0) then
         new_stellar_model = -15
         return
      end if

      call set_comp('o16',XO16)
      if(new_stellar_model/=0) then
         new_stellar_model = -16
         return
      end if

      call set_comp('n14',XN14)
      if(new_stellar_model/=0) then
         new_stellar_model = -17
         return
      end if
      call set_comp('ne20',XNe20)
      if(new_stellar_model/=0) then
         new_stellar_model = -18
         return
      end if

      call set_comp('mg24',XMg24)
      if(new_stellar_model/=0) then
         new_stellar_model = -19
         return
      end if

      call set_comp('si28',XSi28)
      if(new_stellar_model/=0) then
         new_stellar_model = -20
         return
      end if

      call set_comp('s32',XS32)
      if(new_stellar_model/=0) then
         new_stellar_model = -21
         return
      end if

      call set_comp('ar36',XAr36)
      if(new_stellar_model/=0)then
         new_stellar_model = -22
         return
      end if

      call set_comp('ca40',XCa40)
      if(new_stellar_model/=0) then
         new_stellar_model = -23
         return
      end if

      call set_comp('ti44',Xti44)
      if(new_stellar_model/=0)then
         new_stellar_model = -24
         return
      end if

      call set_comp('cr48',Xcr48)
      if(new_stellar_model/=0) then
         new_stellar_model = -25
         return
      end if

      call set_comp('fe52',Xfe52)
      if(new_stellar_model/=0) then
         new_stellar_model = -26
         return
      end if

      call set_comp('fe54',Xfe54)
      if(new_stellar_model/=0) then
         new_stellar_model = -27
         return
      end if

      call set_comp('fe56',Xfe56)
      if(new_stellar_model/=0) then
         new_stellar_model = -28
         return
      end if

      call set_comp('cr56',Xcr56)
      if(new_stellar_model/=0) then
         new_stellar_model = -29
         return
      end if

      call set_comp('ni56',Xni56)
      if(new_stellar_model/=0) then
         new_stellar_model = -30
         return
      end if

      call set_comp('prot',prot)
      if(new_stellar_model/=0) then
         new_stellar_model = -31
         return
      end if

      call set_comp('neut',neut)
      if(new_stellar_model/=0) then
         new_stellar_model = -32
         return
      end if

      ! Relax composition
      call relax_to_new_comp(id, xa, xq, ierr)
      if(ierr/=MESA_SUCESS) then
         new_stellar_model = -33
         return
      end if

      ! Relax entropy profile
      call relax_to_new_entropy(id, xq, temperature, rho, ierr)
      if(ierr/=MESA_SUCESS) then
         new_stellar_model = -34
         return
      end if

      id_new_model = id
      new_model_defined = .true.


      contains 
      
         subroutine set_comp(iso,comp)
            character(len=*) :: iso
            double precision,dimension(n) :: comp
            integer :: net_id, ierr
            
            call get_species_id(id, iso, net_id, ierr)
            if(ierr/=MESA_SUCESS) then
               new_stellar_model = ierr
               return
            end if
            xa(net_id,:) = max(1d-50,comp)
            
         end subroutine set_comp

   end function new_stellar_model

   function finalize_stellar_model(star_id, age_tag)
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
      finalize_stellar_model =  set_age(id_new_model, age_tag)

      new_model_defined = .false.
      finalize_stellar_model = 0
   end function finalize_stellar_model

   function set_age(AMUSE_ID, AMUSE_value)
      implicit none
      integer :: set_age, ierr
      integer, intent(out) :: amuse_id
      double precision, intent(in) :: amuse_value

      ierr = MESA_SUCESS

      call set_new_age(amuse_id, amuse_value, ierr)

      set_age = ierr

   end function set_age


! Set the current mass of the star
   ! This changes an exisiting model, dont use this to set the initial mass of a star
   function set_mass(AMUSE_id, AMUSE_value)
      implicit none
      integer, intent(in) :: AMUSE_id
      double precision, intent(in) :: AMUSE_value
      integer :: set_mass, ierr

      set_mass = -1
      call set_new_mass(AMUSE_id, AMUSE_value, ierr)
      if(ierr/=MESA_SUCESS) return
      set_mass = 0
   end function

! Functions created to map to the se.py interface but they do nothing and allways return failure if called

   integer function set_temperature_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
      integer, intent(in) :: AMUSE_id, AMUSE_zone
      double precision, intent(in) :: AMUSE_value

      set_temperature_at_zone = -1

   end function set_temperature_at_zone

   integer function set_luminosity_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
      integer, intent(in) :: AMUSE_id, AMUSE_zone
      double precision, intent(in) :: AMUSE_value

      set_luminosity_at_zone = -1

   end function set_luminosity_at_zone


   integer function set_mass_fraction_of_species_at_zone(AMUSE_id, &
      AMUSE_species, AMUSE_zone, AMUSE_value)
      integer, intent(in) :: AMUSE_id, AMUSE_zone, AMUSE_species
      double precision, intent(in) :: AMUSE_value
      integer :: ierr

      set_mass_fraction_of_species_at_zone = -1

   end function set_mass_fraction_of_species_at_zone

   integer function set_density_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
      integer, intent(in) :: AMUSE_id, AMUSE_zone
      double precision, intent(in) :: AMUSE_value
      set_density_at_zone = -1
   end function set_density_at_zone


   ! Set the radius at the specified zone/mesh-cell of the star
   integer function set_radius_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
      integer, intent(in) :: AMUSE_id, AMUSE_zone
      double precision, intent(in) :: AMUSE_value
      set_radius_at_zone = -1
   end function set_radius_at_zone






end module AMUSE_mesa