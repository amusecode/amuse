! Library interface to TWIN: useful for putting TWIN into Muse

module twin_library
   use real_kind
   use mesh

   ! Datastructure for storing all the information for one star
   ! These need to be swapped into and out off the TWIN global variables
   ! We only care about single stars and TWIN mode binaries
   type, private :: twin_star_t
      sequence
      ! Flag to indicate wether this star still needs some initialisation
      logical :: virgin

      integer :: number_of_variables
      integer :: number_of_meshpoints

      ! Array of independent variables and increments since last timestep
      real(double), pointer :: h(:,:)
      real(double), pointer :: dh(:,:)

      real(double) :: zams_mass

      ! Stellar type, as in Hurley & al 2000
      integer :: stellar_type

      ! Iteration control
      integer :: startup_iter
      integer :: normal_iter

      ! Binary parameters; these need to be set because TWIN is a binary code
      real(double) :: bms, per, ecc, p

      ! Don't ask...
      integer :: jf

      ! Module current_model_properties
      real(double) :: ml, ql, xl, uc(21)
      integer :: jmod, jnn, jter, joc, jkh

      ! Module test_variables
      real(double) :: dt
      real(double) :: hspn(2), rlf(2), zet(2), xit(2), age, bm, mc(2), om, bper, sm, enc
      real(double) :: tc(2), tfr, t0, m0, mta, om0, omta, a0, ata, e0, eta, cdd
      real(double) :: bp, horb, ro, ra2, rs, secc, tn(2), wmh, wmhe, mh, mhe, mco
      real(double) :: vmg, be(2), lh, lhe, lc, lnu, lth, mcb(8), msb(6), rcb(8), tct(8)
      real(double) :: prev(81), pprev(81)
      integer :: jhold, jm2, jm1


      ! COMMON block SOLV
      real(double) :: er(nvar)

      ! COMMON block STORE
      real(double), pointer :: hpr(:, :)
      !real(double) :: ht(4, nvar)    ! Used in nucleosynthesis code
      !real(double) :: ms(9999)    ! Mass loss history of primary, non-TWIN mode
      !real(double) :: st(9999)    ! Mass loss history of primary, non-TWIN mode
   end type twin_star_t


   ! Inverse permutation: where each element in the global H is
   ! mapped to in the copy in the struct. This makes it possible to
   ! not store unused elements, reducing the memory requirements.
   integer, private :: inv_var_perm(nvar)
   integer, private :: var_perm(nvar)
   integer, private :: actual_number_of_variables

   ! Data structure to store a number of stars in memory.
   integer, private :: max_stars = -1        ! Maximum number of stars
   type(twin_star_t), private, allocatable, target :: star_list(:)
   integer, private :: current_star = 0      ! Currently selected star
   integer, private :: num_stars = 0         ! Number of allocated stars

   ! Some data read from init.dat that is not saved to COMMON blocks
   integer, private :: kh2, kr1, kr2, ksv, kt5, jch, iter

   ! Layout of the ZAMS library file
   real(double), private :: mlo, dm, mhi
   integer, private :: kdm

   ! Number of models to run before switching to "normal" number of
   ! iterations (from "startup" number of iterations)
   integer, parameter, private :: switch_iterations = 5

   ! Print verbose output to stdout yes or no.
   logical, private :: verbose = .false.

   ! Name of the init.dat input file, if not otherwise specified
   character(len=500), private :: init_dat_name = 'init.dat'

   ! Name of the init.run input file, if not otherwise specified
   character(len=500), private :: init_run_name = 'init.run'


contains


   ! set_init_dat_name:
   !  Change the name of the input file used for the numerical settings and
   !  physics options (default: ./init.dat)
   subroutine set_init_dat_name(new_init_dat_name)
      use real_kind
      use file_exists_module
      
      implicit none
      character(len=500), intent(in) :: new_init_dat_name;

      if (.not. file_exists(new_init_dat_name) ) then
         if (verbose) print *, "warning: file ",trim(new_init_dat_name)," for ", trim(init_dat_name), " does not exist!"
      end if
      init_dat_name = trim(new_init_dat_name)

   end subroutine set_init_dat_name


   ! set_init_run_name:
   !  Change the name of the input file used for the numerical settings and
   !  physics options (default: ./init.dat)
   subroutine set_init_run_name(new_init_run_name)
      use real_kind
      use file_exists_module
      
      implicit none
      character(len=500), intent(in) :: new_init_run_name;

      if (.not. file_exists(new_init_run_name) ) then
         if (verbose) print *, "warning: file ",trim(new_init_run_name)," for ", trim(init_run_name), " does not exist!"
      end if
      init_run_name = new_init_run_name;
   end subroutine set_init_run_name


   subroutine open_fortran_output_file(filename, file_number)
      use real_kind
      
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(in) :: file_number
      integer :: status;
      close(file_number)
      open(unit = file_number, action = "write", file = filename, iostat=status)
      return
   end subroutine open_fortran_output_file


   function read_atomic_data_file(filename)
      use real_kind
      
      implicit none
      character(len=*), intent(in) :: filename
      integer :: read_atomic_data_file
      close (26)
      open(unit = 26, action="read", file=filename, iostat=read_atomic_data_file)
      call load_atomic_data(26)
      close (26)
   end function read_atomic_data_file


   function read_reaction_rate_file(filename)
      use real_kind
      use constants
      
      implicit none
      character(len=*), intent(in) :: filename
      integer :: read_reaction_rate_file
      call initialise_constants
      close (42)
      open(unit = 42, action="read", file=filename, iostat=read_reaction_rate_file)
      call load_reaction_neutrino_rates(42)
      close (42)
   end function read_reaction_rate_file


   function read_opacity_file(filename)
      use real_kind
      
      implicit none
      character(len=*), intent(in) :: filename
      integer :: read_opacity_file
      close (20)
      open(unit = 20, action="read", file=filename, iostat=read_opacity_file)
      call load_opacity(20)
      close (20)
   end function read_opacity_file


   function read_co_opacity_file(filename)
      use real_kind
      
      implicit none
      character(len=*), intent(in) :: filename
      integer :: read_co_opacity_file
      close (41)
      open(unit = 41, action="read", file=filename, iostat=read_co_opacity_file)
      call load_opacity_co(41)
      close (41)
   end function read_co_opacity_file


   function read_zams_library_format(filename)
      use real_kind
      
      implicit none
      character(len=*), intent(in) :: filename
      integer :: read_zams_library_format
      close (19)
      open(unit = 19, action="read", file=filename, iostat=read_zams_library_format)
      read (19, *) mlo, dm, mhi, kdm
      close (19)
   end function read_zams_library_format


   subroutine close_zams_library
      close(16)
   end subroutine close_zams_library


   function open_zams_library(filename)
      use real_kind
      
      implicit none
      character(len=*), intent(in) :: filename
      integer :: open_zams_library
      close(16)
      open(unit = 16, action="read", file=filename, iostat=open_zams_library)
   end function open_zams_library


   subroutine set_number_of_startup_iterations(new_iter)
      use real_kind
      
      implicit none
      integer, intent(in) :: new_iter
      kr1 = new_iter
      iter = kr1
   end subroutine set_number_of_startup_iterations



   subroutine set_number_of_iterations(new_iter)
      use real_kind
      
      implicit none
      integer, intent(in) :: new_iter
      kr2 = new_iter
   end subroutine set_number_of_iterations



   function get_number_of_iterations()
      use real_kind
      
      implicit none
      integer :: get_number_of_iterations
      get_number_of_iterations = kr2
   end function get_number_of_iterations



   ! initialise_twin:
   !  General TWIN initialisation: load physics datafiles and ZAMS libraries.
   !  Input variables:
   !   evpath - the path to the stellar evolution code, location of the input data
   !            can be set to '.' if datafiles are in the current directory
   !   nstars - The total number of stars (and binaries) for which we want to
   !            allocate space
   !   zstr   - string describing the metallicity, without the leading 0 and decimal
   !            point. So '02' means '0.02'
   !  Returns value:
   !     0 on succes, non zero on failure
   function initialise_twin(evpath, nstars, zstr)
      use real_kind
      use settings
      use extra_elements
      use file_exists_module
      use constants
      use init_run
      use init_dat
      use current_model_properties
      
      implicit none
      integer :: initialise_twin
      character(len=*), intent(in) ::  evpath, zstr
      integer, intent(in) :: nstars

      integer :: i,ii
      
      integer :: ke1,ke2,ke3,kbc,kev,kfn,kl,kp_var(40),kp_eqn(40),kp_bc(40)
      
      logical :: status, override_run, override_dat
      integer, parameter :: n_inp_files = 15
      character(len=500) :: inputfilenames(n_inp_files), fname
      integer :: inputunits(n_inp_files)     = (/12, 24, 16, 18, 19, 20, 21, 26, 63, 22, 23, 41, 42, 122, 123/)
      integer :: input_required(n_inp_files) = (/ 0,  0,  1,  1,  1,  0,  1,  1, -1, -1, -1,  0,  1, 1, 1/)

      ke1            = id(1)
      ke2            = id(2)
      ke3            = id(3)
      kbc            = id(4)
      kev            = id(5)
      kfn            = id(6)
      kl             = id(7)
      kp_var(1:40)   = id(11:50)
      kp_eqn(1:40)   = id(51:90)
      kp_bc(1:40)    = id(91:130)

      max_stars = nstars
      if (verbose) then
         print *, 'twin initialisation.'
         print *, 'arguments: evpath =', evpath
         print *, '           nstars =', nstars
         print *, '           zstr =', zstr
      end if

      ! Setup input file names
      inputfilenames(1)=trim(evpath)//"/input/zahb"//trim(zstr)//".mod"
      inputfilenames(2)=trim(evpath)//"/input/zahb"//".dat"
      inputfilenames(3)=trim(evpath)//"/input/zams/zams"//trim(zstr)//".mod"
      inputfilenames(4)=trim(evpath)//"/input/zams/zams"//trim(zstr)//".out"
      inputfilenames(5)=trim(evpath)//"/input/zams/zams"//trim(zstr)//".mas"
      inputfilenames(6)=trim(evpath)//"/input/metals/z"//trim(zstr)//"/phys.z"//trim(zstr)
      inputfilenames(7)=trim(evpath)//"/input/lt2ubv.dat"
      inputfilenames(8)=trim(evpath)//"/input/nucdata.dat"
      !INPUTFILENAMES(9)=TRIM(EVPATH)//"/input/mutate.dat"
      inputfilenames(10)=trim(init_dat_name)
      inputfilenames(11)=trim(init_run_name)
      inputfilenames(12)=trim(evpath)//"/input/COtables/COtables_z"//trim(zstr)
      inputfilenames(13)=trim(evpath)//"/input/physinfo.dat"
      inputfilenames(14)=trim(evpath)//"/run/muse/init.dat" ! Sensible defaults
      inputfilenames(15)=trim(evpath)//"/run/muse/init.run" ! Sensible defaults

      ! We want to override the default init.run and init.dat settings if
      ! init.run and init.dat are present in the current directory.
      override_run = .false.
      override_dat = .false.

      ! Check if all input files exist and open them as needed
      do i=1, n_inp_files
         write (fname, '("fort.",i2)') inputunits(i)
         if (file_exists(inputfilenames(i)) .and. .not. file_exists(fname)) then
            if (verbose) print *, "opening ",trim(inputfilenames(i))," for ", trim(fname)
            open(unit = inputunits(i), action="read", file=inputfilenames(i))
         end if
         ! Check if files exist and act appropriately
         ! If the file is required (INPUT_REQUIRED>0), abort on error
         ! If the file is optional (INPUT_REQUIRED==0), give a warning
         ! If the file is probably unneeded (INPUT_REQUIRED<0), do nothing
         if (.not. (file_exists(inputfilenames(i)) .or. file_exists(fname))) then
            if (input_required(i) > 0) then
               write (0, *) 'required input file ', trim(inputfilenames(i)), ' (', trim(fname), ') not found'
               initialise_twin = -1
               return
            else if (input_required(i) == 0) then
               write (0, *) 'warning: input file ', trim(inputfilenames(i)), ' (', trim(fname), ') not found'
            end if
         end if
      end do

      if ( file_exists(inputfilenames(10)) ) override_dat = .true.
      if ( file_exists(inputfilenames(11)) ) override_run = .true.

      ! We need to have a fort.11 - no, we don't, since we're not using star12!
      !WRITE (11, *) 0
      !CLOSE (11)

      ! We need a scratch file for output we don't want or need
      if (.not. file_exists('fort.25'))&
           open (unit=25, action='write', status = 'scratch')

      ! Other redundant (for our purposes) output we do not want to save but cannot
      ! easily avoid: solver output to fort.1 and opacity table output to fort.10
      open (unit=10, action='write', status = 'scratch')
      !OPEN (UNIT=1, ACTION='WRITE', STATUS = 'SCRATCH')

      if (verbose) print *, 'files initialised, initialising equation of state'
      call setsup

      ! Read layout of the ZAMS library
      read (19, *) mlo, dm, mhi, kdm
      rewind (19)

      if (verbose) print *, 'initialised equation of state, reading run settings'
      call read_init_run(123)    ! Load defaults
      if (override_run) call read_init_run(23)
      ! Create short (`pruned') summary of ZAMS models from long file
      ! Are we using this? Not reall, I think...
      !      CALL PRUNER ( 18, 17, ISB )
      !      CLOSE (17)
      !      CLOSE (18)

      ! Read init.dat file
      !> \todo FIXME: this has an unwanted side effect of resetting options that are not set
      !! in the override file to their defaults, even if the default file we loaded
      !! has them set. In other words, we need an "override init.dat" function.
      !<
      if (verbose) print *, 'reading numerical and physical settings'
      status = read_init_dat(122, kh2, kr1, kr2, ksv, kt5, jch)   ! Defaults
      if (status .and. override_dat)&
           status = read_init_dat(22, kh2, kr1, kr2, ksv, kt5, jch)
      if (status .eqv. .false.) then
         initialise_twin = -1
         return
      end if
      if (verbose) print *, 'read settings'
      if (verbose) print *, 'using', kh2, 'meshpoints per star'

      !     Autodetect if we should solve for Mg24 or not by checking if the
      !     corresponding equation is in the list of equations
      use_mg24_eqn = .false.
      do ii = 1, 40
         if (kp_eqn(ii) == emg24 .or. kp_eqn(ii) == esumx) then
            use_mg24_eqn = .true.
            exit
         end if
         if (kp_eqn(ii) == 0) exit   ! Break loop if end of list found
      end do

      ! Convert some things to `cgs' units: 10**11 cm, 10**33 gm, 10**33 erg/s
      cmi = cmi/csy
      cmj = cmj*cmsn/csy
      cms = cms/csy
      cmt = cmt*1.0d-11

      ! Read opacity data and construct splines
      ! KOP (read from init.dat) sets which type of opacity tables
      if (verbose) print *, 'reading opacity tables'
      if (kop<=1) then
         if (verbose) print *, "using opal '92"
         ! Iglesias & Rogers (1992), as implemented by Pols & al. 1995
         call load_opacity(20)
      else
         if (verbose) print *, "using opal '96 with co enhancement"
         ! Iglesias & Rogers (1996), as implemented by Eldridge & Tout 2003
         call load_opacity_co(41)
      end if

      ! Abort if the requested number of meshpoints is larger than the
      ! size of the array
      if (kh2 > nm) then
         write (0, *) 'cannot rezone to ', kh2, 'meshpoints. maximum size is ', nm, 'meshpoints.'
         initialise_twin = -2
         return
      end if

      ! Initialise some more variables
      if ( isb.eq.1 ) ktw = 1
      jb = 1

      ! Setup inverse lookup table for number of variables
      actual_number_of_variables = ke1+ke2+kev
      if (verbose) print *,'will backup',actual_number_of_variables,'variables.'
      inv_var_perm(:) = 0
      do ii=1, actual_number_of_variables
         inv_var_perm(kp_var(ii)) = ii
         var_perm(ii) = kp_var(ii)
      end do

      ! Allocate memory for all the stars we're interested in
      if (allocated(star_list)) deallocate(star_list);
      allocate(star_list(1:max_stars))
      if (verbose) print *, 'allocated memory for ',max_stars, 'stars'

      ! We're not writing any output directly from the library functions
      ! Report success
      initialise_twin = 0;
   end function initialise_twin



   ! Close the temporary files when we're done with them
   ! Sometimes we may need to do this manually
   subroutine close_temporary_files
      use real_kind
      
      implicit none
      close (25)
      close (10)
   end subroutine close_temporary_files



   ! load_zams_star:
   !  load a ZAMS model from disk
   !  Input variables:
   !   m - the mass, in solar units
   !   t - the starting age, in years. Only used for reporting the star's age
   !  Return value:
   !   >0: The stars ID for identifying it in the array of models
   !   =0: No star allocated, out of memory
   function load_zams_star(mass, age)
      use real_kind
      use mesh
      use constants
      use settings
      use current_model_properties
      
      implicit none
      type(twin_star_t), pointer :: star
      integer :: load_zams_star
      real(double), intent(in) :: mass, age

      integer :: im1,kh1,kp1,jmod1,jb1,jn1,jf1,i
      real(double) :: tnuc
      real(double) :: hnuc(50, nm)

      !Common blocks:
      real(double) :: sm1, dty1, age1, per1, bms1, ecc1, p1, enc1

      
      if (verbose) print *, 'create star with mass', mass, 'and age tag', age

      ml = log10(mass)
      sm1 = mass
      im1 = (ml - mlo)/dm + 1.501d0
      call load_star_model(16,im1, h, dh, hnuc, sm1,dty1,age1,per1,bms1,ecc1,p1,enc1,kh1,kp1,jmod1,jb1,jn1,jf1)
      kh = kh1

      ! Special case: a ZAMS star can be loaded even if the library is not
      ! initialised, but then it will not be stored in the array.
      ! This is a convenience for applications that want to call this
      ! function without calling initialise_twin() first
      if (max_stars <= 0) then
         load_zams_star = 1
         return
      end if
      if (num_stars >= max_stars) then
         load_zams_star = 0
         return
      end if
      num_stars = num_stars + 1
      if (verbose) print *, 'setting local pointer to star', num_stars
      star => star_list(num_stars)

      ! Allocate memory to store the star
      star%number_of_variables = actual_number_of_variables
      star%number_of_meshpoints = kh
      allocate(star%h(star%number_of_variables, star%number_of_meshpoints))
      allocate(star%dh(star%number_of_variables, star%number_of_meshpoints))
      allocate(star%hpr(star%number_of_variables, star%number_of_meshpoints))

      ! Make sure we eliminate the redundancy in storing and loading the model
      ! that is in TWIN (main/beginn).
      ! Essentially, we push the initialised data onto the stack
      jnn = 0
      call swap_out_star(num_stars)

      ! Estimate nuclear timescale for this star
      tnuc = 1.0d10 * mass**(-2.8)
      !star%dt = CSY * DTY1
      !star%dt = 3.0e2*CSY
      star%dt = ct3*1.0d-4 * tnuc*csy
      if (mass > 1.0 .and. mass < 1.2) star%dt = 1.0d-2*star%dt
      if (verbose) print *, 'setting initial timestep for star',num_stars,'to',star%dt,'s'

      ! Initialise some more variables
      star%zams_mass = mass
      star%age = age
      if (p1>0) star%p = p1
      star%jhold = 0
      star%jmod = 0
      star%jf = jf1
      star%er(:) = 0.0d0

      star%startup_iter = kr1
      star%normal_iter = kr2

      ! Binary orbital parameters
      if (per1>0) star%per = per1
      if (bms1>0) star%bms = bms1
      if (ecc1>0) star%ecc = ecc1

      ! Determine whether I and phi are computed or not, for OUTPUT
      star%jf = 0
      do i = 11, 50
         if ( id(i)==12 .or. id(i)==14 ) star%jf = star%jf + 1
      end do
      if ( star%jf==1 ) star%jf = 0

      ! The star will still need some initialisation (initial timestep, say)
      star%virgin = .true.

      ! Save all settings
      call swap_in_star(num_stars)
      call swap_out_star(num_stars)

      load_zams_star = num_stars
   end function load_zams_star



   ! evolve_star:
   !  evolve a star for one timestep; essentially star12 without the loop
   !  TODO: ZAHB construction; maybe also for central carbon burning or white dwarfs?
   function twin_evolve()
      use real_kind
      use mesh
      use control
      use constants
      use settings
      use test_variables
      use current_model_properties
      
      implicit none
      integer :: twin_evolve,jo,it
      real(double) :: dty

      ! We need some kind of initialisation for the star, cf the calls
      ! to printb in star12
      ! Actually, we need to pull some data from printb anyway...

      jo = 0
      ! Solve for structure, mesh, and major composition variables
      joc = 1

      dty = dt/csy
      jter = 0
      if (verbose) print *, 'taking timestep ',dty,'yr for star', current_star

      call smart_solver ( iter, id, kt5, jo )
      if (verbose) print *, 'did', jter, 'out of', iter, 'iterations'

      ! Converged if JO == 0
      if ( jo /= 0 ) then
         if (verbose) print *, 'failed to converge on timestep'
         !        If no convergence, restart from 2 steps back, DT decreased substantially
         !        Abort if timestep below limit
         do while (jo /= 0)
            call backup ( dty, jo )
            if ( jo.eq.2 ) then
               twin_evolve = jo
               return
            end if
            call nextdt ( dty, jo, it )

            if (verbose) print *, 'timestep reduced to', dty,'yr'

            if ( jo.eq.3 ) then
               twin_evolve = jo
               return
            end if
            jnn = jnn + 1
            jo = 0
            call smart_solver ( iter, id, kt5, jo )
         end do
      end if
      if (jo == 0) then
         if (verbose) print *, 'converged on timestep'
         call timestep_heuristic ( jo, 1 )
         call update ( dty )
         call nextdt ( dty, jo, 22 )
         if ( jnn > switch_iterations ) iter = kr2
         jnn = jnn + 1
      end if

      ! Report results
      twin_evolve = 0;
   end function twin_evolve

   ! dump_twin_model:
   !  write the state of star id to a named file, for later retrieval.
   !  The file will be in the format of a TWIN input file
   subroutine dump_twin_model(id, filename)
      use real_kind
      
      implicit none
      integer, intent(in) :: id
      character(len=500), intent(in) :: filename

      call select_star(id)

      open (unit=40, file=filename, action='write')
      call output(200, 40, 0, 0)
      close (40);

      call flush_star()
   end subroutine dump_twin_model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! get_*
   !  Selection of getter functions for luminosity, mass, radius, time, next dt etc.
   !  These return retarded information from the *previous* timestep (cf Ross)
   !  This makes sure that we will never need to backtrack beyond this point in
   !  case the code needs to go back to a previous model.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! get_luminosity:
   !  Returns the star's luminosity, in solar units
   function get_luminosity(id)
      use real_kind
      use constants
      
      implicit none
      real(double) :: get_luminosity
      integer, intent(in) :: id

      if (id<1 .or. id>num_stars) then
         get_luminosity = -1.0
         return
      end if

      get_luminosity = star_list(id)%hpr(inv_var_perm(8), 1)/clsn
   end function get_luminosity

   ! get_mass:
   !  Returns the star's mass, in solar units
   function get_mass(id)
      use real_kind
      use constants
      
      implicit none
      real(double) :: get_mass
      integer, intent(in) :: id
      if (id<1 .or. id>num_stars) then
         get_mass = -1.0
         return
      end if

      get_mass = star_list(id)%hpr(inv_var_perm(4), 1)/cmsn
   end function get_mass

   ! get_radius:
   !  Returns the star's radius, in solar units
   function get_radius(id)
      use real_kind
      use constants
      
      implicit none
      real(double) :: get_radius
      integer, intent(in) :: id

      if (id<1 .or. id>num_stars) then
         get_radius = -1.0
         return
      end if

      get_radius = exp(star_list(id)%hpr(inv_var_perm(7), 1))/crsn
   end function get_radius

   ! get_temperature:
   !  Returns the star's effective temperature, in Kelvin
   function get_temperature(id)
      use real_kind
      use constants
      
      implicit none
      real(double) :: get_temperature
      integer, intent(in) :: id

      if (id<1 .or. id>num_stars) then
         get_temperature = -1.0
         return
      end if

      get_temperature = exp(star_list(id)%hpr(inv_var_perm(2), 1))
   end function get_temperature

   ! get_age:
   !  Returns the star's age, in years
   function get_age(id)
      use real_kind
      use constants
      
      implicit none
      real(double) :: get_age
      integer, intent(in) :: id

      if (id<1 .or. id>num_stars) then
         get_age = -1.0
         return
      end if

      get_age = star_list(id)%age
   end function get_age



   ! Return the stellar type of the specified star
   function get_stellar_type(id)
      use real_kind
      
      implicit none
      integer :: get_stellar_type
      integer, intent(in) :: id

      if (id<1 .or. id>num_stars) then
         get_stellar_type = -1
         return
      end if
      get_stellar_type = star_list(id)%stellar_type
   end function get_stellar_type



   subroutine timestep_heuristic ( jo, jstar )
      use real_kind
      
      implicit none
      integer, intent(out) :: jo
      integer, intent(in) :: jstar
      integer :: jcm
      integer, external :: find_stellar_type
      call compute_output_quantities ( jstar )
      call update_timestep_parameters ( jstar )
      call check_stop_conditions ( jstar, jo, jcm, 22 )
      star_list(current_star)%stellar_type = find_stellar_type()
   end subroutine timestep_heuristic



   subroutine set_star_iter_parameters( id, kt1, kt2, jnn )
      use real_kind
      
      implicit none
      integer, intent(in) :: id, kt1, kt2, jnn

      if (id<1 .or. id>num_stars) then
         return
      end if

      if (kt1>0) star_list(id)%startup_iter = kt1
      if (kt2>0) star_list(id)%normal_iter = kt2
      star_list(id)%jnn = jnn
   end subroutine set_star_iter_parameters



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Functions for swapping in and out a particular star: !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! swap_in_star:
   !  Swap in a particular star for evolution with the TWIN code
   !  Assumes all stars are evolved with the same set of equations and physical
   !  parameters. This is fine as long as these are not changed on the fly, as
   !  for the post He-flash generation code. This needs to be disabled from
   !  init.run, or we need to swap in and out the data from init.dat and init.run
   subroutine swap_in_star(star_id)
      use real_kind
      use mesh
      use constants
      use settings
      use solver_global_variables, only: er
      use test_variables
      use current_model_properties
      use binary_history, only: hpr
      
      implicit none
      integer, intent(in) :: star_id
      type(twin_star_t), pointer :: star
      integer :: n, nv, number_of_variables
      integer :: jo,it
      real(double) :: dty,bms,ecc,per,tm,toa
      
      
      ! Nothing to do if the star has already been swapped in
      if (current_star == star_id .and. .not. star_list(star_id)%virgin) return

      if (verbose) print *, 'swapping in star', star_id

      current_star = star_id
      star => star_list(star_id)

      ! Store structure and structure changes
      kh = star%number_of_meshpoints
      number_of_variables = star%number_of_variables
      do n=1, number_of_variables
         nv = var_perm(n)
         h(nv, 1:kh) = star%h(n,1 :kh)
         dh(nv, 1:kh) = star%dh(n, 1:kh)
      end do

      ! Module current_model_properties
      ml = star%ml
      ql = star%ql
      xl = star%xl
      uc(:) = star%uc(:)
      jmod = star%jmod
      jnn = star%jnn
      jter = star%jter
      joc = star%joc
      jkh = star%jkh

      ! Module test_variables
      dt = star%dt

      hspn = star%hspn
      rlf = star%rlf
      zet = star%zet
      xit = star%xit
      age = star%age
      bm = star%bm
      mc = star%mc
      om = star%om
      bper = star%bper
      sm = star%sm
      enc = star%enc

      tc = star%tc
      tfr = star%tfr
      t0 = star%t0
      m0 = star%m0
      mta = star%mta
      om0 = star%om0
      omta = star%omta
      a0 = star%a0
      ata = star%ata
      e0 = star%e0
      eta = star%eta
      cdd = star%cdd

      bp = star%bp
      horb = star%horb
      ro = star%ro
      ra2 = star%ra2
      rs = star%rs
      secc = star%secc
      tn = star%tn
      wmh = star%wmh
      wmhe = star%wmhe
      mh = star%mh
      mhe = star%mhe
      mco = star%mco

      vmg = star%vmg
      be = star%be
      lh = star%lh
      lhe = star%lhe
      lc = star%lc
      lnu = star%lnu
      lth = star%lth
      mcb = star%mcb
      msb = star%msb
      rcb = star%rcb
      tct = star%tct

      prev(:) = star%prev(:)
      pprev(:) = star%pprev(:)
      jhold = star%jhold
      jm2 = star%jm2
      jm1 = star%jm1

      ! COMMON block SOLV. We need to retrieve the typical size of the different
      ! variables. It's not a good idea to leak this information to other stars
      er(1:nvar) = star%er(1:nvar)

      ! COMMON block STORE. We don't need the information about the primary's
      ! mass loss rate
      do n=1, number_of_variables
         nv = var_perm(n)
         hpr(nv, 1:kh) = star%hpr(n,1 :kh)
      end do

      kr1 = star%startup_iter
      kr2 = star%normal_iter
      iter = star%normal_iter
      if (jnn < switch_iterations) iter = star%startup_iter

      ! Does this star need some initialisation?
      if (star%virgin) then
         if (verbose) print *, 'performing one time initialisation for star ', star_id
         star%virgin = .false.

         ! Some tasks from beginn
         jo = 0
         it = 22
         jnn = 0
         dty = dt/csy
         call nextdt ( dty, jo, it )
         ! Do we care about REMESH at this stage? Maybe not!
         bms = cmsn*star%bms
         ecc = star%ecc
         per = star%per
         tm = star%zams_mass * cmsn
         toa = cg1*tm*(bms - tm)*(cg2*per/bms)**c3rd*dsqrt(1.0d0 - ecc*ecc)
         call remesh ( kh2, jch, bms, tm, star%p, ecc, toa, 1, star%jf )

         hpr(:,1:kh) = h(:,1:kh)
         dh(:, 1:kh) = 0.0d0
         jhold = 2
         prev(2:81) = (/hspn, rlf, zet, xit, age, bm, mc, om, bper, sm, enc,  &
              tc, tfr, t0, m0, mta, om0, omta, a0, ata, e0, eta, cdd,  &
              bp, horb, ro, ra2, rs, secc, tn, wmh, wmhe, mh, mhe, mco,  &
              vmg, be, lh, lhe, lc, lnu, lth, mcb, msb, rcb, tct/)
         pprev(1:81) = prev(1:81)
         jm1 = jmod

         ! some tasks from printb, notably the timestep control
         ! For now, we just call timestep_heuristic, which means we get some unwanted output
         call timestep_heuristic ( jo, 1 )
         call nextdt ( dty, jo, it )
         jnn = 1

      end if
   end subroutine swap_in_star



   ! swap_out_star:
   !  Swap out a particular star for later evolution with the TWIN code
   !  Assumes all stars are evolved with the same set of equations and physical
   !  parameters. This is fine as long as these are not changed on the fly, as
   !  for the post He-flash generation code. This needs to be disabled from
   !  init.run, or we need to swap in and out the data from init.dat and init.run
   subroutine swap_out_star(star_id)
      use real_kind
      use mesh
      use solver_global_variables, only: er
      use test_variables
      use current_model_properties
      use binary_history, only: hpr
      
      implicit none
      integer, intent(in) :: star_id
      type(twin_star_t), pointer :: star
      integer :: n, nv, number_of_variables
      
      
      if (star_id < 1) return
      if (verbose) print *, 'swapping out star', star_id
      star => star_list(star_id)

      ! Retrieve structure and structure changes
      number_of_variables = star%number_of_variables
      do n=1, number_of_variables
         nv = var_perm(n)
         star%h(n, 1:kh) = h(nv,1 :kh)
         star%dh(n, 1:kh) = dh(nv, 1:kh)
      end do
      star%number_of_meshpoints = kh

      ! Module current_model_properties:
      star%ml = ml
      star%ql = ql
      star%xl = xl
      star%uc(:) = uc(:)
      star%jmod = jmod
      star%jnn = jnn
      star%jter = jter
      star%joc = joc
      star%jkh = jkh

      ! Module test_variables:
      star%dt  = dt

      star%hspn = hspn
      star%rlf = rlf
      star%zet = zet
      star%xit = xit
      star%age = age
      star%bm = bm
      star%mc = mc
      star%om = om
      star%bper = bper
      star%sm = sm
      star%enc = enc

      star%tc = tc
      star%tfr = tfr
      star%t0 = t0
      star%m0 = m0
      star%mta = mta
      star%om0 = om0
      star%omta = omta
      star%a0 = a0
      star%ata = ata
      star%e0 = e0
      star%eta = eta
      star%cdd = cdd

      star%bp = bp
      star%horb = horb
      star%ro = ro
      star%ra2 = ra2
      star%rs = rs
      star%secc = secc
      star%tn = tn
      star%wmh = wmh
      star%wmhe = wmhe
      star%mh = mh
      star%mhe = mhe
      star%mco = mco

      star%vmg = vmg
      star%be = be
      star%lh = lh
      star%lhe = lhe
      star%lc = lc
      star%lnu = lnu
      star%lth = lth
      star%mcb = mcb
      star%msb = msb
      star%rcb = rcb
      star%tct = tct


      star%prev(:) = prev(:)
      star%pprev(:) = pprev(:)
      star%jhold = jhold
      star%jm2 = jm2
      star%jm1 = jm1

      ! COMMON block SOLV. We need to retrieve the typical size of the different
      ! variables. It's not a good idea to leak this information to other stars
      star%er(1:nvar) = er(1:nvar)

      ! COMMON block STORE. We don't need the information about the primary's
      ! mass loss rate
      do n=1, number_of_variables
         nv = var_perm(n)
         star%hpr(n, 1:kh) = hpr(nv,1 :kh)
      end do

   end subroutine swap_out_star


   ! select_star:
   !  flushes the data from the current star back to the cache and selects the
   !  requested star as the new current star.
   subroutine select_star(id)
      use real_kind
      
      implicit none
      integer :: id
      call swap_out_star(current_star)
      call swap_in_star(id)
   end subroutine select_star


   ! flush_star:
   !  flushes the data from the currently selected star back to the cache
   subroutine flush_star()
      use real_kind
      
      implicit none
      call swap_out_star(current_star)
   end subroutine flush_star


end module twin_library


