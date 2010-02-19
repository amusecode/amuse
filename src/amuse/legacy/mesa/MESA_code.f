   module amuse_support
      implicit none
      character (len=256) :: AMUSE_inlist_path
      character (len=256) :: AMUSE_ZAMS_inlist
      double precision :: AMUSE_metallicity
      contains
      logical function failed(str, ierr)
         character (len=*), intent(in) :: str
         integer, intent(in) :: ierr
         failed = (ierr /= 0)
         if (failed) write(*, *) trim(str) // ' ierr', ierr
      end function failed
   end module amuse_support

! Initialize the stellar evolution code
   subroutine initialize(AMUSE_inlist_path_in, AMUSE_mesa_data_dir, &
         AMUSE_status)
      use amuse_support, only: AMUSE_inlist_path, AMUSE_metallicity, &
         AMUSE_ZAMS_inlist
      use run_star_support
      use ctrls_io, only: set_default_controls
      implicit none
      character(*), intent(in) :: AMUSE_inlist_path_in, AMUSE_mesa_data_dir
      integer :: ierr, AMUSE_status
      ierr = 0
      AMUSE_metallicity = 2.0d-2
      AMUSE_inlist_path = AMUSE_inlist_path_in
      AMUSE_ZAMS_inlist = AMUSE_inlist_path_in // '_zams'
      call set_default_controls
      call do_read_star_job(AMUSE_inlist_path, ierr)
      if (ierr /= 0) then
         write(*,*) "Error while reading inlist:"
         write(*,*) AMUSE_inlist_path
         AMUSE_status = -1
         return
      end if
      ! Replace value of mesa_data_dir just read, with supplied path.
      mesa_data_dir = AMUSE_mesa_data_dir
      call star_init(mesa_data_dir, kappa_file_prefix, &
         net_reaction_filename, net_rate_list_fname, ppn_rate_numbers_fname, ierr)
      if (ierr /= 0) then
         write(*,*) "Error while initializing (in star_init)."
         AMUSE_status = -1
         return
      end if
      AMUSE_status = 0
      return
   end

! Does nothing...
   function initialize_code()
      implicit none
      integer :: initialize_code
      initialize_code = 0
   end function

! Create new ZAMS model
   subroutine new_zams_model()
      use create_zams, only: do_create_zams
      use amuse_support, only: AMUSE_inlist_path, AMUSE_metallicity, &
         AMUSE_ZAMS_inlist, failed
      use star_lib, only: alloc_star, star_setup
      use star_private_def, only: star_info, get_star_ptr
      use run_star_support, only: log_columns_file, profile_columns_file, &
         run_create_zams
      implicit none
      integer :: ierr, AMUSE_id
      type (star_info), pointer :: s
      AMUSE_id = alloc_star(ierr)
      if (failed('alloc_star', ierr)) return
      call get_star_ptr(AMUSE_id, s, ierr)
      if (failed('get_star_ptr', ierr)) return
      call star_setup(AMUSE_id, AMUSE_inlist_path, ierr)
      if (failed('star_setup', ierr)) return
      run_create_zams = .true.
!      call do_create_zams(s, AMUSE_ZAMS_inlist, log_columns_file, &
!         profile_columns_file, ierr)
   end subroutine new_zams_model

! Create a new particle
   function new_particle(AMUSE_id, AMUSE_mass)
      use amuse_support, only: AMUSE_inlist_path, AMUSE_metallicity
      use star_lib, only: alloc_star, star_setup, star_load_zams, show_log_header
      use star_private_def, only: star_info, get_star_ptr
      use run_star_support, only: setup_for_run_star, before_evolve
      implicit none
      integer, intent(out) :: AMUSE_id
      integer :: new_particle, ierr
      double precision, intent(in) :: AMUSE_mass
      type (star_info), pointer :: s
      AMUSE_id = alloc_star(ierr)
      if (ierr /= 0) then
         write(*,*) "Error: could not allocate star."
         new_particle = -1
         return
      end if
      call get_star_ptr(AMUSE_id, s, ierr)
      if (ierr /= 0) then
         write(*,*) "Error: could not get star-pointer for star with ID:", AMUSE_id
         new_particle = -1
         return
      end if
      call star_setup(AMUSE_id, AMUSE_inlist_path, ierr)
      if (ierr /= 0) then
         write(*,*) "Error while setting up new star."
         new_particle = -1
         return
      end if
      ! Replace value of mass and metallicity just read, with supplied values.
      s% initial_mass = AMUSE_mass
      s% initial_z = AMUSE_metallicity
      call star_load_zams(AMUSE_id, ierr)
      if (ierr /= 0) then
         write(*,*) "Error: could not load ZAMS model."
         new_particle = -1
         return
      end if
      call setup_for_run_star(AMUSE_id, s, .false., ierr)
      if (ierr /= 0) then
         write(*,*) "Error: could not setup star."
         new_particle = -1
         return
      end if
      call before_evolve(AMUSE_id, ierr)
      if (ierr /= 0) then
         write(*,*) "Error in before_evolve."
         new_particle = -1
         return
      end if
      call show_log_header(AMUSE_id, ierr)
      if (ierr /= 0) then
         write(*,*) "Error in show_log_header."
         new_particle = -1
         return
      end if
      new_particle = 0
   end function

! Remove a particle
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
      get_number_of_particles = 0
   end function

! Return the metallicity parameter
   function get_metallicity(AMUSE_value)
      use ctrls_io
      use run_star_support
      implicit none
      integer :: get_metallicity
      double precision :: AMUSE_value
      AMUSE_value = initial_z
      get_metallicity = 0
   end function

! Set the metallicity parameter
   function set_metallicity(AMUSE_value)
      use ctrls_io
      use run_star_support
      implicit none
      integer :: set_metallicity
      double precision :: AMUSE_value
      initial_z = AMUSE_value
      set_metallicity = 0
   end function

! Return the current mass of the star
   function get_mass(AMUSE_id, AMUSE_value)
      use ctrls_io
      use run_star_support
      implicit none
      integer, intent(in) :: AMUSE_id
      integer :: get_mass
      double precision :: AMUSE_value
      AMUSE_value = initial_z
      get_mass = 0
   end function

! Return the current temperature of the star
   function get_temperature(AMUSE_id, AMUSE_value)
      use ctrls_io
      use run_star_support
      implicit none
      integer, intent(in) :: AMUSE_id
      integer :: get_temperature
      double precision :: AMUSE_value
      AMUSE_value = initial_z
      get_temperature = 0
   end function

! Return the current luminosity of the star
      function get_luminosity(AMUSE_id, AMUSE_value)
         use ctrls_io
         use run_star_support
         implicit none
         integer, intent(in) :: AMUSE_id
         integer :: get_luminosity
         double precision :: AMUSE_value
         AMUSE_value = initial_z
         get_luminosity = 0
      end function

! Return the current age of the star
      function get_age(AMUSE_id, AMUSE_value)
         use ctrls_io
         use run_star_support
         implicit none
         integer, intent(in) :: AMUSE_id
         integer :: get_age
         double precision :: AMUSE_value
         AMUSE_value = initial_z
         get_age = 0
      end function

! Return the current radius of the star
      function get_radius(AMUSE_id, AMUSE_value)
         use ctrls_io
         use run_star_support
         implicit none
         integer, intent(in) :: AMUSE_id
         integer :: get_radius
         double precision :: AMUSE_value
         AMUSE_value = initial_z
         get_radius = 0
      end function

! Return the current stellar type of the star
      function get_stellar_type(AMUSE_id, AMUSE_value)
         use ctrls_io
         use run_star_support
         implicit none
         integer, intent(in) :: AMUSE_id
         integer :: get_stellar_type
         double precision :: AMUSE_value
         AMUSE_value = initial_z
         get_stellar_type = 0
      end function

! Evolve the star until AMUSE_end_time
      subroutine evolve(AMUSE_id, AMUSE_end_time)
         use run_star_support
         use run_star, only: do_run_star
         use star_private_def, only: star_info, get_star_ptr
         use amuse_support, only: failed
         implicit none
         integer, intent(in) :: AMUSE_id
         double precision, intent(in) :: AMUSE_end_time
         type (star_info), pointer :: s
         integer :: ierr
         call get_star_ptr(AMUSE_id, s, ierr)
         if (failed('get_star_ptr', ierr)) return
         s% max_age = AMUSE_end_time
         call do_run_star

         return
      end
      
