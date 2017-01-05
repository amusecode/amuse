! New (updated) library interface, for linking TWIN to AMUSE.
! The main practical design change is that we no longer calculate one timestep "ahead" of where we claim to be.
! This is a useful practice if the top-layer can't handle "backups", but it introduces a lag in the response from
! the evolution code to any proposed changes. This is ok if the stellar evolution mainly talks in one direction to, say, a
! dynamics code, as in the work by Church et al, but it is not so useful in a tightly coupled framework like AMUSE.
! In general, this iteration is much simpler and closer to how TWIN works internally, which should make this (much) simpler in the
! long run at the expense of a (bit of) extra memory
!
! TODO:
!  initialise a new star from a file
!  write a model to a file
!  append a model to a file (are these really different?)
!  join/break binary and relevant query functions
!  Test whether timestepping is done correctly.
#include "assert.h"

module twinlib
   use real_kind
   use indices, only: nvar, nvstar

   type, private :: twin_star_t
      integer :: id        ! ID for this star
      integer :: pid, sid  ! ID for primary/secondary star (in a binary)
      integer :: bid       ! ID for binary that this star is a member of (-1 if none, id for the binary itself)

      ! Some properties
      logical :: nucleosynthesis ! Solve extended nucleosynthesis network for this model

      ! Flag to indicate wether this star still needs some initialisation
      logical :: virgin

      ! Flag that tells whether this star still exists (i.e. it has not been removed)
      logical :: exists

      integer :: number_of_variables
      integer :: number_of_meshpoints

      ! Array of independent variables and increments since last timestep
      real(double), pointer :: h(:,:)
      real(double), pointer :: dh(:,:)
      real(double), pointer :: hpr(:, :)
      real(double), pointer :: ht(:, :)
      real(double), pointer :: Hnucpr(:,:,:)
      real(double), pointer :: Hnuc(:,:,:)
      real(double), pointer :: DHnuc(:,:,:)
      real(double), pointer :: menc(:,:)

      real(double) :: maximum_mass
      real(double) :: zams_mass

      real(double) :: age

      ! Stellar type, as in Hurley & al 2000
      integer :: stellar_type

      ! Iteration control
      integer :: startup_iter
      integer :: normal_iter

      ! Binary parameters; these need to be set because TWIN is a binary code
      real(double) :: bms, per, ecc, p

      ! Timestep parameters
      real(double) :: rlf_prev(2)
      real(double) :: qcnv_prev(2)
      real(double) :: lnuc_prev(2)
      real(double) :: lhe_prev(2)
      real(double) :: lh_prev(2)

      ! Module current_model_properties
      real(double) :: uc(21)
      real(double) :: dt
      integer :: jmod, jnn

      integer :: eqns(130)

      ! Module test_variables
      ! Not all of these matter, but better safe than sorry
      real(double) :: mc(2)      ! Mass scale for the core (for mesh spacing function)
      real(double) :: hspn(2)    ! Spin angular momentum, star 1 and star 2
      real(double) :: rlf(2)     ! Roche lobe filling factor, ln R*/RL
      real(double) :: zet(2)     ! Wind mass loss from the system [1e33g/s]
      real(double) :: xit(2)     ! Mass transfer to companion, by RLOF and wind
      real(double) :: tn(2)      ! Nuclear timescale; set in printb [s]
      real(double) :: be(2)      ! Binding energy of the stellar envelope [erg/(1Mo)]
      real(double) :: be0(2)     ! Binding energy of the stellar envelope: gravity [erg/(1Mo)]
      real(double) :: be1(2)     ! Binding energy of the stellar envelope: internal energy [erg/(1Mo)]
      real(double) :: be2(2)     ! Binding energy of the stellar envelope: recombination energy [erg/(1Mo)]
      real(double) :: be3(2)     ! Binding energy of the stellar envelope: H2 association energy [erg/(1Mo)]
      real(double) :: spsi(2)    ! Scaled surface potential, in -G(M1+M2)/a
      real(double) :: bm         ! Total mass in the binary [1e33 g]
      real(double) :: om         ! Total mass of the secondary [1e33 g]
      real(double) :: bper       ! Binary period [days]
      real(double) :: enc        ! Artificial energy generation rate
      real(double) :: cdd        ! Timestep control parameter
      integer :: jhold           ! If jhold<3, then don't change the timestep
      real(double) :: prev(81), pprev(81)
      integer :: jm2, jm1


      ! Wind/mass accretion options
      real(double) :: cmi
      real(double) :: cmdot_wind
   end type twin_star_t

   ! Data structure to store a number of stars in memory.
   integer, private :: max_stars = -1        ! Maximum number of stars
   type(twin_star_t), private, allocatable, target :: star_list(:)
   integer, private :: current_star = 0      ! Currently selected star
   integer, private :: sx_updated_for_star = 0 ! Stellar structure data sx is up-to-date for this star

   logical, private :: initialised = .false.

   ! Print verbose output to stdout yes or no.
   logical, private :: verbose = .false.

   ! Number of models to run before switching to "normal" number of
   ! iterations (from "startup" number of iterations)
   integer, parameter, private :: switch_iterations = 5

   ! Static (local) init.dat options
   integer, private :: wanted_kh, ksv, kt5, jch

   ! Layout of the ZAMS library file
   real(double), private :: mlo, dm, mhi
   integer, private :: kdm

   ! Solver list/equations for single stars and binaries
   integer, private :: eqns_single(130)
   integer, private :: eqns_binary(130)

   ! Temporary storage of amuse parameters
   character(len = 1000), private :: amuse_ev_path
   integer, private :: amuse_nstars, amuse_nmesh, amuse_kion
   logical :: amuse_verbose
   real(double), private :: amuse_Z, amuse_csmc, amuse_calp, amuse_cos
   real(double), private :: amuse_cth, amuse_maxage, amuse_mindt
   real(double) :: amuse_entropy_accuracy, amuse_entropy_force

   ! List private subroutines that should not be called directly
   private initialise_stellar_parameters, allocate_star, swap_in, swap_out, select_star, make_zahb_model

contains

   ! initialise_twin:
   !  General TWIN initialisation: load physics datafiles and ZAMS libraries.
   !  Input variables:
   !   path   - the path to the stellar evolution code, location of the input data
   !            can be set to '.' if datafiles are in the current directory. Leave blank ('')
   !            to use the default (either from environment variable or from configure option)
   !   nstars - The total number of stars (and binaries) for which we want to
   !            allocate space
   !   Z      - Desired metallicity, in the range 0.0 <= Z <= 0.04
   !   verb   - verbose output (true/false); optional.
   !   nmesh  - maximum number of mesh points (optional)
   !  Returns value:
   !     0 on succes, non zero on failure:
   !    -1 initialised before (not critical)
   !    -2 bad metallicity requested (out-of-range)
   !    -3 bad number of stars requested (<1)
   !    -4 Cannot load init.dat settings file
   integer function initialise_twin(path, nstars, Z, verb, nmesh)
      use real_kind
      use settings
      use constants
      use current_model_properties
      use control
      use settings
      use filenames
      use install_location
      use opacity
      use distortion
      use allocate_arrays
      use mesh, only: ktw, isb, max_nm, id
      use solver_global_variables, only: solver_output_unit
      use polytrope, only: NMP
      implicit none
      character(len=*), intent(in) :: path
      integer, intent(in)          :: nstars
      real(double), intent(in)     :: Z
      logical, intent(in),optional :: verb
      integer, intent(in),optional :: nmesh
      character(len=80)            :: tmpstr
      integer                      :: i, j

      ! Assume success, we'll override when needed
      initialise_twin = 0

      if (present(verb)) verbose = verb

      ! Only initialise once
      if (initialised) then
         initialise_twin = -1
         return
      end if

      ! Test whether requested metallicity is valid
      if (Z < 0.0 .or. Z > 0.04) then
         initialise_twin = -2
         return
      end if

      ! Sensible number of stars?
      if (nstars < 1) then
         initialise_twin = -3
         return
      end if

      ! Initialise path
      if (len(trim(path)) > 0) then
         evpath = trim(path)
      else
         call get_environment_variable("evpath", evpath)
         if (len(trim(evpath)) == 0) evpath = twin_install_path;
      end if

      ! Set metallicity
      czs = Z
      ! Now convert the value of czs to a string
      write(tmpstr, '(f10.7)') czs
      i = 1
      j = len(trim(tmpstr))
      do while (tmpstr(i:i) /= '.')
         i = i+1
      end do
      i = i + 1
      do while (tmpstr(j:j) == '0')
         j = j-1
      end do
      j = max(j,i)
      zstr = tmpstr(i:j)

      if (present(nmesh)) max_nm = max(max_nm, nmesh)
      max_nm = max(max_nm, NMP)

      if (verbose) then
         print *, 'twin initialisation.'
         print *, 'arguments: evpath   = ', trim(evpath)
         print *, '           nstars   = ', nstars
         print *, '           zstr     = ', trim(zstr)
         print *, '           max mesh = ', max_nm
      end if

      ! Read constant data for the run (init.dat)
      ! The first one just sets up the equations for doing a binary run
      inputfilenames(9) = trim(evpath)//'/input/amuse/init_twin.dat'
      initialise_twin = read_initdat_settings(inputfilenames(9))
      if (initialise_twin /= 0) return
      eqns_binary = id

      inputfilenames(9) = trim(evpath)//'/input/amuse/init.dat'
      initialise_twin = read_initdat_settings(inputfilenames(9))
      if (initialise_twin /= 0) return
      eqns_single = id

      ! Allocate memory for global arrays
      call allocate_global_arrays(max_nm)
      assert(allocated(h))

      ! Decide what opacity tables are required and whether ZAMS libraries are required or optional
      ! ZAMS files are *not* required, but if they don't exist we will need to construct a starting model from scratch.
      input_required(3) = -1
      input_required(4) = -1
      if (kop == 4) then
         if (verbose) print *, '*** Warning: CO enhanced tables not recommended/tested under AMUSE'
         input_required(5) = -1
         input_required(11) = 1
      else
         input_required(5) = 0
         input_required(11) = -1
      end if

      ! Input file init.run is *not* required
      input_required(10) = -1

      ! Default names for input files
      if (verbose) print *, 'Set location of input files'
      call set_default_filenames

      ! Check whether all files that are required exist
      if (verbose) print *, 'Checking if input files exist'
      call assert_input_files_exist

      ! Read opacity data and construct splines
      ! KOP (read from init.dat) sets which type of opacity tables
      if (verbose) print *, 'Load opacity table for opacity option ', kop
      call load_opacity(kop, czs)

      ! Initialise remaining data tables
      call setsup
      call initialise_distortion_tables
      call load_nucleosynthesis_rates(inputunits(13), inputunits(14))

      ! Read format of ZAMS library
      if (have_zams_library()) then
         ! Read layout of the ZAMS library
         read (19, *) mlo, dm, mhi, kdm
         rewind (19)
         if (verbose) print *, 'ZAMS library format ', mlo, dm, mhi, kdm
      end if

      ! Initialise number of stars that will be solved concurrently (default: 1)
      ktw = 1
      isb = 1

      ! Redirect output from the solver to a terminal
      solver_output_unit = (/ 6, 6, 6 /)

      ! Reserve some extra space for stars.
      ! This is friendly, because we will need to allocate an extra "star" whenever a binary is created
      max_stars = 2*nstars

      allocate(star_list(1:max_stars))
      do i=1, max_stars
         star_list(i)%exists = .false.
         star_list(i)%nucleosynthesis = .false.
      end do
      if (verbose) print *, 'allocated memory for ',max_stars, 'stars+binaries'
   end function initialise_twin



   ! read_initdat_settings:
   !  read options from a specified init.dat file.
   !  NOTE: some settings are actually overridden below, meaning the configuration file is ignored.
   !  This is mostly for options that will (probably) not work (correctly) anyway.
   integer function read_initdat_settings(filename)
      use filenames
      use settings
      use control
      use init_dat
      implicit none
      character(len=*), intent(in) :: filename
      integer :: ioerror

      read_initdat_settings = 0

      if (verbose) print *, 'Reading settings from ', trim(filename)
      open(unit = inputunits(9), action="read", file=filename, iostat = ioerror)
      if (ioerror /= 0 .or. read_init_dat(inputunits(9), wanted_kh, ksv, kt5, jch) .eqv. .false.) then
         read_initdat_settings = -4
         return
      end if
      rewind (inputunits(9))

      ! Override some settings from init.dat
      ! TODO: some of these we do want to be able to change...
      cmi_mode = 2                        ! Interpret CMI as accretion rate in Msun/year
      cmdot_wind = 1.0d0                  ! Enable stellar winds
      store_changes = .true.              ! We do want to store (predicted) changes when asked to store a file
      use_quadratic_predictions = .false. ! Use linear prediction for variable updates
      jch = 3                             ! Always construct new mesh spacing, safer
   end function read_initdat_settings



   ! Returns true if a zams library file exists (so that we can load models from it)
   logical function have_zams_library()
      use file_exists_module
      use filenames
      implicit none

      have_zams_library = file_exists(inputfilenames(3)) .and. file_exists(inputfilenames(4))
   end function have_zams_library



   ! Allocate a star from the list
   integer function allocate_star()
      implicit none
      integer :: i

      allocate_star = 0
      do i = 1, max_stars
         if (.not. star_list(i)%exists) then
            star_list(i)%exists = .true.
            star_list(i)%id = i
            star_list(i)%pid = -1
            star_list(i)%sid = -1
            star_list(i)%bid = -1
            star_list(i)%nucleosynthesis = .false.
            allocate_star = i
            return
         end if
      end do
   end function allocate_star



   subroutine release_star(star_id)
      implicit none
      integer, intent(in) :: star_id
      type(twin_star_t), pointer :: star

      if (star_id_is_valid(star_id) < 0) return

      if (current_star == star_id) current_star = 0

      star => star_list(star_id)

      ! Free memory used by this star
      star%exists = .false.
      deallocate(star%h)
      deallocate(star%dh)
      deallocate(star%hpr)
      deallocate(star%menc)
      if (star%nucleosynthesis) then
         deallocate(star%ht)
      end if
      star%nucleosynthesis = .false.
   end subroutine
   integer function delete_star(star_id)
      implicit none
      integer, intent(in) :: star_id
      call release_star(star_id)
      delete_star = 0
   end function


   ! Set (initialise) basic stellar variables, like the equation list and timestep
   subroutine initialise_stellar_parameters(star_id)
      use real_kind
      use test_variables, only: dt, age, jhold
      use stopping_conditions, only: uc
      use current_model_properties
      use constants
      use settings
      use test_variables, only: enc
      implicit none
      integer, intent(in)        :: star_id
      type(twin_star_t), pointer :: star
      real(double)               :: tnuc, mass

      star => star_list(star_id)
      current_star = star_id

      mass = star%zams_mass

      star%eqns = eqns_single
      star%startup_iter = kr1
      star%normal_iter = kr2

      star%uc = (/ 1.00E-01, 0.00E+00, 1.00E+02, 0.00E+00, 3.00E+00, 5.30E+00, 1.20E+00, &
                   6.30E+00, 3.00E+02, 0.00E+00, 1.00E-06, 1.00E+06, 1.00E+03, 1.00E+03, &
                   0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00 /)
      star%uc(2) = amuse_maxage
      star%uc(12) = amuse_mindt
      uc = star%uc
      star%maximum_mass = -1.0
      jmod = 0
      jnn = 0
      star%jmod = jmod
      star%jnn = jnn

      rlf_prev = (/0.0, 0.0/)   ! Previous value of Roche-lobe filling factor
      qcnv_prev = (/0.0, 0.0/)  ! Previous value of mass fraction of convective envelope
      lnuc_prev = (/0.0, 0.0/)  ! Previous value of nuclear burning luminosity
      lhe_prev = (/0.0, 0.0/)   ! Previous value of He burning luminosity
      lh_prev = (/0.0, 0.0/)    ! Previous value of H burning luminosity

      star%exists = .true.

      ! Estimate nuclear timescale for this star
      tnuc = 1.0d10 * mass**(-2.8)
      star%dt = ct3*1.0d-3 * tnuc*csy
      if (mass > 1.0 .and. mass < 1.2) star%dt = 1.0d-2*star%dt
      if (verbose) print *, 'setting initial timestep for star',star_id,'to',star%dt/csy,'yr'
      dt = star%dt
      age = star%age
      enc = 0.0d0

      jhold = 2

      ! ----------------------------------------------------------
      ! Compute quantities
      ! ----------------------------------------------------------
      call compute_output_quantities ( 1 )

      ! -------------------------------------------------------------------
      ! Update quantities that have to be calculated explicitly
      ! -------------------------------------------------------------------
      call update_explicit_quantities( 1 )

      ! -------------------------------------------------------------------
      ! Update the control parameters for the next timestep
      ! -------------------------------------------------------------------
      call update_timestep_parameters( 1 )
   end subroutine initialise_stellar_parameters



   ! new_zams_star:
   !  Create a new ZAMS model for a star
   !  Output variables:
   !   star_id    - The stars ID for identifying it in the array of models
   !  Input variables:
   !   mass       - the initial mass, in solar units
   !   start_age  - (optional) the starting age, in years. Only used for reporting the star's age (default: 0.0)
   !   nmesh      - (optional) the number of gridpoints in the model (default: whatever was read in from init.dat)
   !   wrot       - (optional) the rotation rate for the model (default: no rotation)
   !  Return value:
   !    0: Success
   !   -1: No star allocated, requested mesh is too large
   !   -3: No star allocated, out of memory
   integer function new_zams_star(star_id, mass, start_age, nmesh, wrot)
      use real_kind
      use mesh, only: nm, h, hpr, dh, kh
      use mesh_enc, only: menc
      use nucleosynthesis, only: ht_nvar, hnuc
      use constants
      use settings
      use polytrope
      use indices
      use test_variables
      use stopping_conditions
      implicit none
      integer, intent(out)                :: star_id
      real(double), intent(in)            :: mass
      real(double), optional, intent(in)  :: start_age
      integer, optional, intent(in)       :: nmesh
      real(double), optional, intent(in)  :: wrot
      type(twin_star_t), pointer          :: star
      integer                             :: new_kh
      integer                             :: new_id
      real(double)                        :: w
      real(double)                        :: sm1,dty1,age1,per1,bms1,ecc1,p1,enc1, tm, oa
      integer                             :: kh1,kp1,jmod1,jb1,jn1,jf1, im1
      real(double)                        :: hn1(50, nm)
      star_id = -1
      new_zams_star = -1

      ! Test if we can actually load a model from disk.
      ! If we cannot, we construct a pre-mainsequence star (that is, a polytrope) instead and evolve it to ZAMS
      if (have_zams_library() .eqv. .false.) then
         ! Construct a pre-mainsequence star
         new_zams_star = new_prems_star(star_id, mass, start_age, nmesh, wrot)

         ! Evolve to ZAMS
         ! *TODO*
         if (verbose) print *, 'Evolve to ZAMS'
         return
      end if

      ! Set default value for number of gridpoints and rotation rate
      new_kh = wanted_kh
      if (present(nmesh)) new_kh = nmesh

      w = 0.0d0
      if (present(wrot)) w = wrot

      if (new_kh > NM .or. new_kh < 1) then
         if (verbose) print *, 'Cannot create model at ', new_kh, 'meshpoints. Maximum size is ', NM, 'meshpoints.'
         new_zams_star = -1
         return
      end if

      new_id = allocate_star()
      if (new_id == 0) then
         new_zams_star = -3
         return
      end if
      star => star_list(new_id)

      ! Allocate memory for this star
      star%number_of_variables = NVSTAR
      star%number_of_meshpoints = new_kh

      allocate(star%h(star%number_of_variables, star%number_of_meshpoints))
      allocate(star%dh(star%number_of_variables, star%number_of_meshpoints))
      allocate(star%hpr(star%number_of_variables, star%number_of_meshpoints))
      allocate(star%menc(2, star%number_of_meshpoints))
      star%menc = 0.0d0

      ! Nucleosynthesis?
      ! *TODO* Make sure we can actually set this per star
      if (star%nucleosynthesis) then
         allocate(star%ht(ht_nvar, star%number_of_meshpoints))
      end if

      star%cmdot_wind = 1.0d0   ! Enable stellar winds
      call select_star(new_id)

      star%zams_mass = mass
      if (present(start_age)) then
         star%age = start_age
      else
         star%age = 0.0d0
      endif

      ! Load model
      if (star%nucleosynthesis .and. verbose) print *, '*** Warning: ZAMS model+nucleosynthesis is not reliable.'
      if (verbose) print *, 'Load ZAMS model'
      im1 = (log10(mass) - mlo)/dm + 1.501d0
      call load_star_model(16,im1, h, dh, hn1, sm1,dty1,age1,per1,bms1,ecc1,p1,enc1,kh1,kp1,jmod1,jb1,jn1,jf1)

      ! Set desired options, this is a single star (by construction)
      kh = kh1
      tm = cmsn * mass
      bm = cmsn * bms1
      om = bm - tm
      p1 = 2.*cpi / (w * csday + 1.0e-16)
      bper = per1
      oa = cg1*tm*om*(cg2*bper/bm)**c3rd*sqrt(1.0d0 - ecc1*ecc1)
      ! Remesh to desired numer of mesh points
      call remesh ( new_kh, jch, bm, tm, p1, ecc1, oa, 1, 2 )

      hpr = h

      call initialise_stellar_parameters(new_id)

      call swap_out()

      star%uc = uc
      new_zams_star = 0
      star_id = new_id
   end function new_zams_star



   ! new_prems_star:
   !  Create a new pre-main sequence model for a star
   !  Output variables:
   !   star_id    - The stars ID for identifying it in the array of models
   !  Input variables:
   !   mass       - the initial mass, in solar units
   !   start_age  - (optional) the starting age, in years. Only used for reporting the star's age (default: 0.0)
   !   nmesh      - (optional) the number of gridpoints in the model (default: whatever was read in from init.dat)
   !   wrot       - (optional) the rotation rate for the model (default: no rotation)
   !  Return value:
   !    0: Success
   !   -1: No star allocated, requested mesh is too large
   !   -3: No star allocated, out of memory
   integer function new_prems_star(star_id, mass, start_age, nmesh, wrot)
      use real_kind
      use mesh, only: nm, h, hpr, dh, kh
      use mesh_enc, only: menc
      use nucleosynthesis, only: ht_nvar, hnuc
      use constants
      use settings
      use polytrope
      use test_variables
      use stopping_conditions
      implicit none
      integer, intent(out)                :: star_id
      real(double), intent(in)            :: mass
      real(double), optional, intent(in)  :: start_age
      integer, optional, intent(in)       :: nmesh
      real(double), optional, intent(in)  :: wrot
      type(twin_star_t), pointer          :: star
      integer                             :: new_kh
      integer                             :: new_id
      real(double)                        :: w
      real(double)                        :: sm1,dty1,age1,per1,bms1,ecc1,p1,enc1, tm, oa
      integer                             :: kh1,kp1,jmod1,jb1,jn1,jf1
      real(double)                        :: hn1(50, nm)

      star_id = -1
      new_prems_star = -3

      ! Set default value for number of gridpoints and rotation rate
      new_kh = wanted_kh
      if (present(nmesh)) kh1 = nmesh

      w = 0.0d0
      if (present(wrot)) w = wrot

      if (new_kh > NM .or. new_kh < 1) then
         if (verbose) print *, 'Cannot create model at ', new_kh, 'meshpoints. Maximum size is ', NM, 'meshpoints.'
         new_prems_star = -1
         return
      end if

      new_id = allocate_star()
      if (new_id == 0) then
         return
      end if
      star => star_list(new_id)

      ! Allocate memory for this star
      star%number_of_variables = NVSTAR
      star%number_of_meshpoints = new_kh
      allocate(star%h(star%number_of_variables, star%number_of_meshpoints))
      allocate(star%dh(star%number_of_variables, star%number_of_meshpoints))
      allocate(star%hpr(star%number_of_variables, star%number_of_meshpoints))
      allocate(star%menc(2, star%number_of_meshpoints))
      star%menc = 0.0d0

      ! Nucleosynthesis?
      ! *TODO* Make sure we can actually set this per star
      if (star%nucleosynthesis) then
         allocate(star%ht(ht_nvar, star%number_of_meshpoints))
      end if

      call select_star(new_id)

      star%zams_mass = mass
      if (present(start_age)) then
         star%age = start_age
      else
         star%age = 0.0d0
      endif

      ! Construct pre-main sequence model
      if (verbose) print *, 'Construct pre-main sequence model'
      call generate_starting_model(mass, h, dh, hn1, sm1,dty1,age1,per1,bms1,ecc1,p1,enc1,kh1,kp1,jmod1,jb1,jn1,jf1)

      ! Set desired options, this is a single star (by construction)
      kh = kh1
      tm = cmsn * mass
      bm = cmsn * bms1
      om = bm - tm
      p1 = 2.*cpi / (w * csday + 1.0e-16)
      if (w == 0.0d0) p1 = 1.0d6
      bper = per1
      oa = cg1*tm*om*(cg2*bper/bm)**c3rd*sqrt(1.0d0 - ecc1*ecc1)
      ecc1 = 0.0
      ! Remesh to desired numer of mesh points
      call remesh ( new_kh, jch, bm, tm, p1, ecc1, oa, 1, 2 )

      hpr = h

      call initialise_stellar_parameters(new_id)

      call swap_out()

      star%uc = uc
      new_prems_star = 0
      star_id = new_id
   end function new_prems_star



   ! new_star_from_file:
   !  Create a new pre-main sequence model for a star
   !  Output variables:
   !   star_id    - The stars ID for identifying it in the array of models
   !  Input variables:
   !   filename   - the name of the file the model is stored in
   !   start_age  - (optional) the starting age, in years (default: whatever the age of the saved model is)
   !   nmesh      - (optional) the number of gridpoints in the model (default: whatever was read in from init.dat)
   !   wrot       - (optional) the rotation rate for the model (default: whatever is in the file)
   !  Return value:
   !    0: Success
   !   -1: No star allocated, requested mesh is too large
   !   -2: No star allocated, file not found
   !   -3: No star allocated, out of memory
   integer function new_star_from_file(star_id, filename, start_age, nmesh, wrot)
      use real_kind
      use mesh, only: nm, h, hpr, dh, kh
      use nucleosynthesis, only: ht_nvar, hnuc
      use constants
      use settings
      use polytrope
      use indices
      use test_variables
      use filenames
      use stopping_conditions
      implicit none
      integer, intent(out)                :: star_id
      character(len=*), intent(in)        :: filename
      real(double), optional, intent(in)  :: start_age
      integer, optional, intent(in)       :: nmesh
      real(double), optional, intent(in)  :: wrot
      type(twin_star_t), pointer          :: star
      integer                             :: new_kh
      integer                             :: new_id
      integer                             :: ip1, ioerror
      real(double)                        :: sm1,dty1,age1,per1,bms1,ecc1,p1,enc1, tm, oa
      integer                             :: kh1,kp1,jmod1,jb1,jn1,jf1
      real(double)                        :: hn1(50, nm)
      star_id = -1

      ! Set default value for number of gridpoints and rotation rate
      new_kh = wanted_kh
      if (present(nmesh)) new_kh = nmesh

      if (new_kh > NM .or. new_kh < 1) then
         if (verbose) print *, 'Cannot create model at ', new_kh, 'meshpoints. Maximum size is ', NM, 'meshpoints.'
         new_star_from_file = -1
         return
      end if

      ip1 = get_free_file_unit()
      open(unit = ip1, action="read", file=filename, iostat=ioerror)
      if (ioerror /= 0) then
         if (verbose) print *, 'Cannot load file ', trim(filename), '.'
         new_star_from_file = -2
         return
      end if

      new_id = allocate_star()
      if (new_id == 0) then
         new_star_from_file = -3
         return
      end if
      star => star_list(new_id)

      ! Allocate memory for this star
      star%number_of_variables = NVSTAR
      star%number_of_meshpoints = new_kh
      allocate(star%h(star%number_of_variables, star%number_of_meshpoints))
      allocate(star%dh(star%number_of_variables, star%number_of_meshpoints))
      allocate(star%hpr(star%number_of_variables, star%number_of_meshpoints))
      allocate(star%menc(2, star%number_of_meshpoints))
      star%menc = 0.0d0

      ! Nucleosynthesis?
      ! *TODO* Make sure we can actually set this per star
      if (star%nucleosynthesis) then
         allocate(star%ht(ht_nvar, star%number_of_meshpoints))
      end if

      call select_star(new_id)

      ! Load model
      if (star%nucleosynthesis .and. verbose) print *, '*** Warning: ZAMS model+nucleosynthesis is not reliable.'
      call load_star_model(ip1,1, h, dh, hn1, sm1,dty1,age1,per1,bms1,ecc1,p1,enc1,kh1,kp1,jmod1,jb1,jn1,jf1)
      close(ip1)

      if (verbose) print *, 'Loaded model with mass ', sm1, ' and age ', age1
      star%zams_mass = sm1
      star%age       = age1

      ! Optionally override options
      if (present(start_age)) star%age = start_age
      if (present(wrot))      p1 = 2.*cpi / (wrot * csday + 1.0e-16)

      ! Set desired options, this is a single star (by construction)
      kh = kh1
      tm = cmsn * sm1
      bm = cmsn * bms1
      om = bm - tm
      bper = per1
      oa = cg1*tm*om*(cg2*bper/bm)**c3rd*sqrt(1.0d0 - ecc1*ecc1)
      ! Remesh to desired numer of mesh points
      call remesh ( new_kh, jch, bm, tm, p1, ecc1, oa, 1, 2 )

      hpr = h

      call initialise_stellar_parameters(new_id)

      call swap_out()

      star%uc = uc
      new_star_from_file = 0
      star_id = new_id
   end function new_star_from_file



!~! Create a new particle
!~   integer function new_particle(star_id, mass)
!~      implicit none
!~      integer, intent(out) :: star_id
!~      real(double), intent(in) :: mass
!~      real(double) :: start_age
!~      start_age = 0.0
!~      star_id = new_zams_star(mass, start_age)
!~      if (star_id .lt. 1) then
!~        new_particle = -1
!~      else
!~        new_particle = 0
!~      end if
!~   end function



   ! write_star_to_file:
   !  write the state of star id to a named file, for later retrieval.
   !  The file will be in the format of a TWIN input file
   integer function write_star_to_file(id, filename)
      use real_kind
      use filenames, only: get_free_file_unit

      implicit none
      integer, intent(in) :: id
      character(len=*), intent(in) :: filename
      integer :: ip

      call select_star(id)

      ip = get_free_file_unit()
      open (unit=ip, file=filename, action='write')
      call output(200, ip, 0, 4)
      close (ip);
      write_star_to_file = 0
   end function write_star_to_file




   logical function is_binary_system(star_id)
      implicit none
      integer, intent(in) :: star_id
      is_binary_system = star_list(star_id)%exists .and. star_list(star_id)%bid /= -1
   end function is_binary_system



   logical function is_single_star(star_id)
      implicit none
      integer, intent(in) :: star_id
      is_single_star = star_list(star_id)%exists .and. star_list(star_id)%bid == -1
   end function is_single_star




   ! Select the star with star_id as the "current" star
   subroutine swap_in(star_id)
      use real_kind
      use mesh, only: nm, kh, h, dh, hpr
      use mesh_enc, only: menc, ip_mesh_size
      use test_variables
      use current_model_properties, only: jmod, jnn, rlf_prev, qcnv_prev, lnuc_prev, lhe_prev, lh_prev
      use nucleosynthesis, only: nucleosynthesis_enabled
      use settings, only: cmi, cmdot_wind
      use solver_global_variables, only: es
      use indices
      implicit none
      integer, intent(in) :: star_id
      type(twin_star_t), pointer :: star, primary, secondary
      if (current_star == star_id) return

      if (star_id_is_valid(star_id) < 0) return
      star => star_list(star_id)

      sx_updated_for_star = 0

      assert(allocated(h))

      kh = star%number_of_meshpoints
      if (is_single_star(star_id)) then
         h(1:star%number_of_variables, 1:kh) = star%h(1:star%number_of_variables, 1:kh)
         dh(1:star%number_of_variables, 1:kh) = star%dh(1:star%number_of_variables, 1:kh)
         hpr(1:star%number_of_variables, 1:kh) = star%hpr(1:star%number_of_variables, 1:kh)
         menc(1:2, 1:kh) = star%menc(1:2, 1:kh)
         ip_mesh_size = kh
      else
         assert(star%pid > -1)
         assert(star%sid > -1)
         primary => star_list(star%pid)
         secondary => star_list(star%sid)

         ! Primary star
         h(1:primary%number_of_variables, 1:kh)   = primary%h(1:primary%number_of_variables, 1:kh)
         dh(1:primary%number_of_variables, 1:kh)  = primary%dh(1:primary%number_of_variables, 1:kh)
         hpr(1:primary%number_of_variables, 1:kh) = primary%hpr(1:primary%number_of_variables, 1:kh)
         menc(1, 1:kh) = primary%menc(1, 1:kh)

         ! Orbital elements
         h(INDEX_ORBIT_VAR_START+1:INDEX_ORBIT_VAR_START+star%number_of_variables, 1:kh)   = &
            star%h(1:star%number_of_variables, 1:kh)
         dh(INDEX_ORBIT_VAR_START+1:INDEX_ORBIT_VAR_START+star%number_of_variables, 1:kh)  = &
            star%dh(1:star%number_of_variables, 1:kh)
         hpr(INDEX_ORBIT_VAR_START+1:INDEX_ORBIT_VAR_START+star%number_of_variables, 1:kh) = &
            star%hpr(1:star%number_of_variables, 1:kh)

         ! Secondary star
         h(INDEX_SECONDARY_START+1:INDEX_SECONDARY_START+secondary%number_of_variables, 1:kh)   = &
            secondary%h(1:secondary%number_of_variables, 1:kh)
         dh(INDEX_SECONDARY_START+1:INDEX_SECONDARY_START+secondary%number_of_variables, 1:kh)  = &
            secondary%dh(1:secondary%number_of_variables, 1:kh)
         hpr(INDEX_SECONDARY_START+1:INDEX_SECONDARY_START+secondary%number_of_variables, 1:kh) = &
            secondary%hpr(1:secondary%number_of_variables, 1:kh)
         menc(2, 1:kh) = secondary%menc(1, 1:kh)
      end if
      age  = star%age
      dt   = star%dt
      jmod = star%jmod
      jnn  = star%jnn

      mc = star%mc
      hspn = star%hspn
      rlf = star%rlf
      zet = star%zet
      xit = star%xit
      tn = star%tn
      be = star%be
      be0 = star%be0
      be1 = star%be1
      be2 = star%be2
      be3 = star%be3
      spsi = star%spsi
      bm = star%bm
      om = star%om
      bper = star%bper
      enc = star%enc
      cdd = star%cdd
      jhold = star%jhold
      prev = star%prev
      pprev = star%pprev
      jm2 = star%jm2
      jm1  = star%jm1

      ! Correctly set binary mass eigenvalue (not otherwise preserved for single stars, but needed)
      h(VAR_BMASS, 1:kh) = bm

      rlf_prev = star%rlf_prev
      qcnv_prev = star%qcnv_prev
      lnuc_prev = star%lnuc_prev
      lhe_prev = star%lhe_prev
      lh_prev = star%lh_prev

      cmi = star%cmi
      cmdot_wind = star%cmdot_wind

      nucleosynthesis_enabled = star%nucleosynthesis

      ! TODO: we can do better than just setting this to 0 every time (for instance, we could just store it)
      es = 0

      current_star = star_id
   end subroutine swap_in



   ! Backup the properties of the current star
   subroutine swap_out()
      use real_kind
      use mesh, only: nm, kh, h, dh, hpr
      use mesh_enc, only: menc
      use test_variables
      use current_model_properties, only: jmod, jnn, rlf_prev, qcnv_prev, lnuc_prev, lhe_prev, lh_prev
      use settings, only: cmi, cmdot_wind
      use indices
      implicit none
      type(twin_star_t), pointer :: star, primary, secondary

      if (current_star == 0) return

      star => star_list(current_star)
      if (.not. star%exists) return

      star%number_of_meshpoints = kh
      if (is_single_star(current_star)) then
         star%h(1:star%number_of_variables, 1:kh) = h(1:star%number_of_variables, 1:kh)
         star%dh(1:star%number_of_variables, 1:kh) = dh(1:star%number_of_variables, 1:kh)
         star%hpr(1:star%number_of_variables, 1:kh) = hpr(1:star%number_of_variables, 1:kh)
         star%menc(1:2, 1:kh) = menc(1:2, 1:kh)
      else
         assert(star%pid > -1)
         assert(star%sid > -1)
         primary => star_list(star%pid)
         secondary => star_list(star%sid)

         ! Primary star
         primary%h(1:primary%number_of_variables, 1:kh) = h(1:primary%number_of_variables, 1:kh)
         primary%dh(1:primary%number_of_variables, 1:kh) = dh(1:primary%number_of_variables, 1:kh)
         primary%hpr(1:primary%number_of_variables, 1:kh) = hpr(1:primary%number_of_variables, 1:kh)
         primary%menc(1, 1:kh) = menc(1, 1:kh)

         ! Orbital elements
         star%h(1:star%number_of_variables, 1:kh) = &
            h(INDEX_ORBIT_VAR_START+1:INDEX_ORBIT_VAR_START+star%number_of_variables, 1:kh)
         star%dh(1:star%number_of_variables, 1:kh) = &
            dh(INDEX_ORBIT_VAR_START+1:INDEX_ORBIT_VAR_START+star%number_of_variables, 1:kh)
         star%hpr(1:star%number_of_variables, 1:kh) = &
            hpr(INDEX_ORBIT_VAR_START+1:INDEX_ORBIT_VAR_START+star%number_of_variables, 1:kh)

         ! Secondary star
         secondary%h(1:secondary%number_of_variables, 1:kh) = &
            h(INDEX_SECONDARY_START+1:INDEX_SECONDARY_START+secondary%number_of_variables, 1:kh)
         secondary%dh(1:secondary%number_of_variables, 1:kh) = &
            dh(INDEX_SECONDARY_START+1:INDEX_SECONDARY_START+secondary%number_of_variables, 1:kh)
         secondary%hpr(1:secondary%number_of_variables, 1:kh) = &
            hpr(INDEX_SECONDARY_START+1:INDEX_SECONDARY_START+secondary%number_of_variables, 1:kh)
         secondary%menc(1, 1:kh) = menc(2, 1:kh)
      end if
      star%age  = age
      star%dt   = dt
      star%jmod = jmod
      star%jnn  = jnn

      star%mc = mc
      star%hspn = hspn
      star%rlf = rlf
      star%zet = zet
      star%xit = xit
      star%tn = tn
      star%be = be
      star%be0 = be0
      star%be1 = be1
      star%be2 = be2
      star%be3 = be3
      star%spsi = spsi
      star%bm = bm
      star%om = om
      star%bper = bper
      star%enc = enc
      star%cdd = cdd
      star%jhold = jhold
      star%prev = prev
      star%pprev = pprev
      star%jm2 = jm2
      star%jm1  = jm1

      star%rlf_prev = rlf_prev
      star%qcnv_prev = qcnv_prev
      star%lnuc_prev = lnuc_prev
      star%lhe_prev = lhe_prev
      star%lh_prev = lh_prev

      star%cmi = cmi
      star%cmdot_wind = cmdot_wind
   end subroutine swap_out



   ! Swap out the current star, swap in the new star
   subroutine select_star(star_id)
      implicit none
      integer, intent(in) :: star_id

      call swap_out()
      call swap_in(star_id)
   end subroutine select_star



   subroutine make_zahb_model(star_id, jo, dty)
      use test_variables, only: dt, age, jhold, lhe, mhe
      use printb_global_variables, only: sdc
      use current_model_properties, only: jmod, jnn, joc, jb
      use stopping_conditions
      use constants, only: csy, cg, cmsn, cpi, csday
      use mesh, only: h, ktw, isb, kh
      use settings
      use indices
      use resolve_helium_flash
      implicit none
      integer, intent(in) :: star_id
      integer, intent(inout) :: jo
      real(double), intent(inout) :: dty
      real(double) :: tm, oa, ecc, bms, omega, vi, p1
      integer :: iter
      integer :: Jstar

      call backup ( dty, jo )
      call update_quantities_if_needed(star_id)

      if ( lhe > uc(5) .and. sdc > uc(6) .and. mhe == 0.0d0 ) then   ! Helium flash
         tm  = h(VAR_MASS,  1)
         oa  = h(VAR_HORB,  1)
         ecc = h(VAR_ECC,   1)
         bms = h(VAR_BMASS, 1)
         omega = h(VAR_OMEGA, 1)
         vi  = h(VAR_INERT, 1)
         p1 = 2.0*cpi/(h(VAR_OMEGA, 1) * csday)
         if (verbose) print *, 'He flash (construct new model)'
         call make_post_flash_model(jo)
         tm  = h(VAR_MASS,  1)
         call remesh ( kh, jch, bms, tm, p1, ecc, oa, 1, 0 )
         dty = 1.0d3
         ! Conserve angular momentum
         ! FIXME: this does not really work properly for differentially rotating stars
         h(VAR_OMEGA, 1:kh) = omega*vi / h(VAR_INERT, 1)
         if (verbose) print *, 'Start ZAHB'
         jo = 0
      end if
   end subroutine make_zahb_model



   ! evolve_one_timestep:
   !  Advance the evolution of a star one timestep.
   !  BEWARE: in the case of non-convergence, the code will do a backup and retry from the previous model
   !  with a reduced timestep. If this occurs, the age of the star after this call is lower than it was before.
   ! Returns value:
   !  returns 0 if converged ok, non-zero value (positive or negative) indicate convergence failure!
   integer function evolve_one_timestep(star_id, myverb)
      use test_variables, only: dt, age, jhold, lhe, mhe
      use printb_global_variables, only: sdc
      use current_model_properties, only: jmod, jnn, joc, jb
      use stopping_conditions
      use constants, only: csy, cg, cmsn, cpi, csday
      use mesh, only: h, ktw, isb, kh
      use settings
      use indices
      use resolve_helium_flash
      implicit none
      integer, intent(in) :: star_id
      logical, optional, intent(in) :: myverb
      type(twin_star_t), pointer :: star
      real(double) :: dty, tdyn
      real(double) :: tm, oa, ecc, bms, omega, vi, p1
      integer :: iter
      integer :: jo
      integer :: Jstar
      logical :: my_verbose

      if (star_id_is_valid(star_id) < 0) then
         evolve_one_timestep = -1
         return
      end if

      my_verbose = verbose
      if (present(myverb)) my_verbose = myverb

      evolve_one_timestep = 0
      if (star_id /= current_star) call select_star(star_id)
      assert(current_star > 0)

      star => star_list(star_id)

      ! Clear return code
      jo = 0
      ! Solve for structure, mesh, and major composition variables
      joc = 1

      ! Do star 1 in a binary
      Jstar = 1
      jb = 1
      isb = 1
      ktw = 1

      dty = dt/csy
      if (my_verbose) print *, 'taking timestep ',dty,'yr for star', current_star

      ! Determine number of allowed iterations before backup
      iter = star%normal_iter
      if (jnn < switch_iterations) iter = star%startup_iter

      ! Set maximum mass the star can reach through accretion, basically infinite
      if (star%maximum_mass < 0.0d0) then
         star%uc(13) = 2.d0 * star%H(VAR_MASS, 1) / cmsn
      else
         star%uc(13) = star%maximum_mass
      end if

      ! Set timestep control/minimum timestep
      tdyn = 1. / (csy*sqrt(cg * h(VAR_MASS, 1) / exp(3.*h(VAR_LNR, 1))))
      star%uc(12) = tdyn * csy
      uc = star%uc

      ! Don't show iteration output, unless we're printing verbose output.
      kt5 = iter
      if (my_verbose) kt5 = 0

      call smart_solver ( iter, star%eqns, kt5, jo )

      ! Converged if JO == 0
      if ( jo /= 0 ) then
         if (my_verbose) print *, 'failed to converge on timestep'
         ! If no convergence, restart from 2 steps back, DT decreased substantially
         do while (jo /= 0)
            call backup ( dty, jo )
            if ( jo == 2 ) then
               evolve_one_timestep = -1000-jo
               ! If not converged, figure out why.
               ! If we're at the He flash, build a ZAHB model.
               call make_zahb_model(star_id, jo, dty)
               tdyn = 1. / (csy*sqrt(cg * h(VAR_MASS, 1) / exp(3.*h(VAR_LNR, 1))))
               if (jo /= 0) return
            end if
            call nextdt ( dty, jo, 22 )

            if (my_verbose) print *, 'timestep reduced to', dty,'yr'

            ! If the timestep is (well) below the dynamical timescale for the star, abort
            if (dty < tdyn) then
               evolve_one_timestep = -1002
               return
            end if

            if ( jo == 3 ) then
               evolve_one_timestep = -1002
               return
            end if
            jnn = jnn + 1
            jo = 0
            call smart_solver ( iter, star%eqns, kt5, jo )
         end do
      end if
     
      if(jo /= 0) then
        evolve_one_timestep = -1000-jo
      else
         evolve_one_timestep = 0
      end if
      if (my_verbose) then
         call update_quantities_if_needed(star_id)
         call write_summary ( 1, jnn, 6 )
      end if

      ! If converged, update nucleosynthesis (if wanted)
      ! *TODO*

      ! If not converged, diagnose mode of failure.
      !
      ! If close to the helium flash, try to skip over the problem by "going around"
      ! This is a multi-step process:
      ! 1. Evolve a reference model (3 Msun) until He ignition
      ! 2. Strip the envelope until we are left with a low-mass He star model
      ! 3. Follow the steps of FGB2HB:
      !    Accrete material of the desired surface composition. Allow H to burn.
      !    Stop accretion when the star has reached the desired mass.
      !    Stop evolution when the star has reached the desired core mass.
      !    Fix the composition profile of the envelope, based on the last pre-flash model.
      ! *TODO*
      ! *TODO*: make a switch to enable or disable this behaviour
      ! *TODO*: we can simplify this process by constructing the core directly as a low-mass He star.
      !         Rewrite fgb2hb to work like this.
      !
      ! If stuck on the way to the white dwarf cooling track, eliminate the H shell.

      if (my_verbose) print *, 'converged on timestep'
      ! ----------------------------------------------------------
      ! Compute quantities
      ! ----------------------------------------------------------
      call compute_output_quantities ( Jstar )

      ! -------------------------------------------------------------------
      ! Update quantities that have to be calculated explicitly
      ! -------------------------------------------------------------------
      call update_explicit_quantities( Jstar )

      ! -------------------------------------------------------------------
      ! Update the control parameters for the next timestep
      ! -------------------------------------------------------------------
      call update_timestep_parameters( Jstar )

      call update ( dty )
      call nextdt ( dty, jo, 22 )
      jnn = jnn + 1

      if (dty < tdyn) dty = tdyn * 1.1
      dt = dty * csy

      if (my_verbose) print *, 'new timestep', dty,'yr'

      ! Synchronise the star-list with the evolution code
      call swap_out()
   end function evolve_one_timestep



   integer function evolve_until_model_time(star_id, time)
      use stopping_conditions, only: uc
      use test_variables, only: dt
      use constants
      implicit none
      integer, intent(in) :: star_id
      real(double), intent(in) :: time
      real(double) :: age, dty, dty_min, dty_max
      evolve_until_model_time = 0

      if (star_id /= current_star) call select_star(star_id)

      dty_min = uc(12)/csy
      do while (age_of(star_id) < time)
         ! Tweak the timestep when we are close to the target age so that we don't overshoot,
         ! but also don't suddenly change the timestep srastically
         age = age_of(star_id)
         dty = dt / csy
         dty_max = time - age
         if ( age+2*dty < time .and. age+3*dty >= time) then
            ! We expect three more timesteps, start constraining the timestep
            dty = 0.4*max(dty_max, dty_min )
         end if
         if ( age+dty < time .and. age+2*dty >= time) then
            ! We expect maybe two more timesteps, constrain
            dty = 0.6*max(dty_max, dty_min )
         end if
         if ( age+dty >= time .and. age<time) then
            ! This is our final timestep
            dty = min(dty_max, dty )
         end if

         ! Age reached
         if ( dty_max <= dty_min ) then
            evolve_until_model_time = 0
            return
         end if
         dt = dty * csy

         evolve_until_model_time = evolve_one_timestep(star_id)
         if (age_of(star_id) >= uc(2)) evolve_until_model_time = 5

         ! Abort in case of convergence failure
         if (evolve_until_model_time /= 0) return
      end do
   end function evolve_until_model_time



   integer function synchronise_stars(star_id1, star_id2)
      implicit none
      integer, intent(in) :: star_id1, star_id2

      synchronise_stars = 0
      if (age_of(star_id1) < age_of(star_id2)) then
         synchronise_stars = evolve_until_model_time(star_id1, age_of(star_id2))
      elseif (age_of(star_id2) < age_of(star_id1)) then
         synchronise_stars = evolve_until_model_time(star_id2, age_of(star_id1))
      endif
   end function synchronise_stars



   ! Join the two stars in a binary with the given orbital period (in days) and eccentricity
   integer function join_binary(id1, id2, period, ecc)
      use test_variables, only: dt, age
      use current_model_properties
      use constants
      use indices
      implicit none
      integer, intent(in) :: id1, id2
      real(double), intent(in) :: period, ecc
      type(twin_star_t), pointer :: binary, primary, secondary
      real(double) :: tm, om, bm, oa, bper
      integer :: new_id

      ! Make sure the two single stars are not already part of a binary system
      if (.not. (is_single_star(id1) .and. is_single_star(id2))) then
         join_binary = 0
         return
      end if

      ! Make sure arrays are current
      if (id1 == current_star .or. id2 == current_star) call swap_out()

      new_id = allocate_star()
      if (new_id == 0) then
         join_binary = 0
         return
      end if
      binary => star_list(new_id)

      if (star_list(id1)%zams_mass > star_list(id2)%zams_mass) then
         binary%pid = id1
         binary%sid = id2
      else
         binary%pid = id2
         binary%sid = id1
      end if

      primary => star_list(binary%pid)
      secondary => star_list(binary%sid)

      ! Calculate orbit variables from orbital elements
      tm = cmsn * mass_of(binary%pid)
      om = cmsn * mass_of(binary%sid)
      bm = tm + om
      bper = period
      oa = cg1*tm*om*(cg2*bper/bm)**c3rd*sqrt(1.0d0 - ecc*ecc)

      ! Allocate memory for this star
      binary%number_of_variables = NVBIN
      binary%number_of_meshpoints = primary%number_of_meshpoints
      allocate(binary%h(binary%number_of_variables, binary%number_of_meshpoints))
      allocate(binary%dh(binary%number_of_variables, binary%number_of_meshpoints))
      allocate(binary%hpr(binary%number_of_variables, binary%number_of_meshpoints))
      allocate(binary%menc(2, binary%number_of_meshpoints))
      binary%menc = 0.0d0

      binary%h(VAR_HORB, :) = oa
      binary%h(VAR_ECC, :) = ecc
      binary%h(VAR_XI, :) = 0.0d0
      binary%h(VAR_BMASS, :) = bm
      binary%h(VAR_PMASS, :) = tm

      binary%hpr(VAR_HORB, :) = oa
      binary%hpr(VAR_ECC, :) = ecc
      binary%hpr(VAR_XI, :) = 0.0d0
      binary%hpr(VAR_BMASS, :) = bm
      binary%hpr(VAR_PMASS, :) = tm

      binary%dh = 0.0d0

      current_star = new_id
      binary%eqns = eqns_binary
      binary%startup_iter = primary%startup_iter
      binary%normal_iter = max(primary%normal_iter, secondary%normal_iter)

      binary%uc = (/ 1.00E-01, 2.00E+12, 1.00E+02, 0.00E+00, 3.00E+00, 5.30E+00, 1.20E+00, &
                     6.30E+00, 3.00E+02, 0.00E+00, 1.00E-06, 1.00E+06, 1.00E+03, 1.00E+03, &
                     0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00 /)
      binary%maximum_mass = -1.0

      rlf_prev = (/0.0, 0.0/)   ! Previous value of Roche-lobe filling factor
      qcnv_prev = (/0.0, 0.0/)  ! Previous value of mass fraction of convective envelope
      lnuc_prev = (/0.0, 0.0/)  ! Previous value of nuclear burning luminosity
      lhe_prev = (/0.0, 0.0/)   ! Previous value of He burning luminosity
      lh_prev = (/0.0, 0.0/)    ! Previous value of H burning luminosity

      binary%exists = .true.

      ! Estimate nuclear timescale for this star
      binary%dt = min(primary%dt, secondary%dt)
      dt = binary%dt
      age = binary%age

      call swap_out()

      join_binary = new_id
   end function join_binary



   ! Get global stellar properties:
   ! age (in years)
   real(double) function age_of(star_id)
      implicit none
      integer, intent(in) :: star_id
      integer :: tmp
      tmp = get_age(star_id, age_of)
   end function age_of
   integer function get_age(star_id, age)
      use real_kind
      implicit none
      integer, intent(in) :: star_id
      real(double), intent(out) :: age
      get_age = -1

      age = -1.0

      if (star_id_is_valid(star_id) < 0) return
      age = star_list(star_id)%age
      get_age = 0
   end function get_age



   ! Get luminosity (in solar units)
   real(double) function luminosity_of(star_id)
      implicit none
      integer, intent(in) :: star_id
      integer :: tmp
      tmp = get_luminosity(star_id, luminosity_of)
   end function luminosity_of
   integer function get_luminosity(star_id, luminosity)
      use real_kind
      use constants
      use indices
      implicit none
      integer, intent(in) :: star_id
      real(double), intent(out) :: luminosity
      get_luminosity = -1

      luminosity = -1.0

      if (star_id_is_valid(star_id) < 0) return
      luminosity = star_list(star_id)%H(VAR_LUM, 1) / CLSN
      get_luminosity = 0
   end function get_luminosity




   ! Get mass (in solar units)
   real(double) function mass_of(star_id)
      implicit none
      integer, intent(in) :: star_id
      integer :: tmp
      tmp = get_mass(star_id, mass_of)
   end function mass_of
   integer function get_mass(star_id, mass)
      use indices
      use constants
      implicit none
      integer, intent(in) :: star_id
      real(double), intent(out) :: mass
      get_mass = -1
      mass = -1.0

      if (star_id_is_valid(star_id) < 0) return
      mass = star_list(star_id)%H(VAR_MASS, 1) / CMSN
      get_mass = 0
   end function get_mass



   ! Get radius (in solar units)
   real(double) function radius_of(star_id)
      implicit none
      integer, intent(in) :: star_id
      integer :: tmp
      tmp = get_radius(star_id, radius_of)
   end function radius_of
   integer function get_radius(star_id, radius)
      use indices
      use constants
      use settings, only: ct
      implicit none
      integer, intent(in) :: star_id
      real(double), intent(out) :: radius
      real(double) :: r
      get_radius = -1

      radius = -1.0

      if (star_id_is_valid(star_id) < 0) return
      r = sqrt(exp(2.*star_list(star_id)%H(VAR_LNR, 1)) - CT(8))
      radius = r / CRSN
      get_radius = 0
   end function get_radius



   ! Get effective temperature (in Kelvin)
   real(double) function temperature_of(star_id)
      implicit none
      integer, intent(in) :: star_id
      integer :: tmp
      tmp = get_temperature(star_id, temperature_of)
   end function temperature_of
   integer function get_temperature(star_id, temperature)
      use real_kind
      use constants
      use indices
      implicit none
      integer, intent(in) :: star_id
      real(double), intent(out) :: temperature
      get_temperature = -1

      temperature = -1.0

      if (star_id_is_valid(star_id) < 0) return
      temperature = exp(star_list(star_id)%h(VAR_LNT, 1))
      get_temperature = 0
   end function get_temperature



   ! Get timestep (in yr)
   real(double) function timestep_of(star_id)
      implicit none
      integer, intent(in) :: star_id
      integer :: tmp
      tmp = get_time_step(star_id, timestep_of)
   end function timestep_of
   integer function get_time_step(star_id, timestep)
      use real_kind
      use constants
      implicit none
      integer, intent(in) :: star_id
      real(double), intent(out) :: timestep
      get_time_step = -1

      timestep = -1.0

      if (star_id_is_valid(star_id) < 0) return
      timestep = star_list(star_id)%dt / csy
      get_time_step = 0
   end function get_time_step




   ! Set timestep (in yr)
   subroutine set_timestep(star_id, dty)
      use real_kind
      use constants
      use test_variables, only: dt
      implicit none
      integer, intent(in) :: star_id
      real(double), intent(in) :: dty

      if (star_id_is_valid(star_id) < 0) return

      call select_star(star_id)
      dt = dty * csy
      star_list(star_id)%dt = dt
   end subroutine set_timestep


   ! Set some evolution options
   ! Maximum mass the star can reach before accretion is turned off (in solar units).
   ! Can be used to increase the mass of the star to a particular point.
   ! Setting it to -1 allows the mass of the star to grow indefinitely (unless the code breaks first)
   subroutine set_maximum_mass_after_accretion(star_id, mmass)
      implicit none
      integer, intent(in) :: star_id
      real(double), intent(in) :: mmass

      if (star_id_is_valid(star_id) < 0) return

      star_list(star_id)%maximum_mass = mmass
   end subroutine set_maximum_mass_after_accretion



   ! Set the accretion rate for this star, in Msun/yr (negative for winds / mass loss)
   integer function set_manual_mass_transfer_rate(star_id, mdot)
      use constants, only: csy
      use settings, only: cmi
      implicit none
      integer, intent(in) :: star_id
      real(double), intent(in) :: mdot
      set_manual_mass_transfer_rate = -1
      if (star_id_is_valid(star_id) < 0) return

      call select_star(star_id)
      cmi = mdot / csy
      star_list(star_id)%cmi = cmi
      set_manual_mass_transfer_rate = 0
   end function set_manual_mass_transfer_rate
   integer function get_manual_mass_transfer_rate(star_id, mdot)
      use constants, only: csy
      implicit none
      integer, intent(in) :: star_id
      real(double), intent(out) :: mdot
      get_manual_mass_transfer_rate = -1
      if (star_id_is_valid(star_id) < 0) return
      mdot = star_list(star_id)%cmi * csy
      get_manual_mass_transfer_rate = 0
   end function get_manual_mass_transfer_rate

   ! Stellar wind switch: can be modulated between 0.0 (no wind) and 1.0 (full strength)
   integer function set_wind_multiplier(star_id, mdot_factor)
      use settings, only: cmdot_wind
      implicit none
      integer, intent(in) :: star_id
      real(double), intent(in) :: mdot_factor
      set_wind_multiplier = -1

      if (star_id_is_valid(star_id) < 0) return

      call select_star(star_id)
      cmdot_wind = mdot_factor
      star_list(star_id)%cmdot_wind = mdot_factor
      set_wind_multiplier = 0
   end function set_wind_multiplier
   integer function get_wind_multiplier(star_id, mdot_factor)
      implicit none
      integer, intent(in) :: star_id
      real(double), intent(out) :: mdot_factor
      get_wind_multiplier = -1
      if (star_id_is_valid(star_id) < 0) return
      mdot_factor = star_list(star_id)%cmdot_wind
      get_wind_multiplier = 0
   end function get_wind_multiplier



! get_AGB_wind_setting:
! Get the current setting for mass-loss (AGB)
   integer function get_AGB_wind_setting(value)
      use real_kind
      use massloss
      implicit none
      integer, intent(out) :: value
      real(double) :: tmpVW, tmpW
      tmpVW = multiplier_vasiliadis_wood
      tmpW = multiplier_wachter
      if ((tmpVW.eq.0.0) .and. (tmpW.eq.1.0)) then
        value = 1
        get_AGB_wind_setting = 0
      else if ((tmpVW.eq.1.0) .and. (tmpW.eq.0.0)) then
        value = 2
        get_AGB_wind_setting = 0
      else
        value = 0
        get_AGB_wind_setting = -1
      endif
   end function



! set_AGB_wind_setting:
! Set the current setting for mass-loss (AGB)
   integer function set_AGB_wind_setting(value)
      use massloss
      implicit none
      integer, intent(in) :: value
      if (value .eq. 1) then
        multiplier_vasiliadis_wood = 0.0
        multiplier_wachter = 1.0
        set_AGB_wind_setting = 0
      else if (value .eq. 2) then
        multiplier_vasiliadis_wood = 1.0
        multiplier_wachter = 0.0
        set_AGB_wind_setting = 0
      else
        set_AGB_wind_setting = -1
      endif
   end function



! get_RGB_wind_setting:
! Get the current setting for mass-loss (RGB)
   integer function get_RGB_wind_setting(value)
      use real_kind
      use massloss
      implicit none
      real(double), intent(out) :: value
      if (multiplier_schroeder .gt. 0.0) then
        value = multiplier_schroeder
      else
        value = -multiplier_reimers
      endif
      get_RGB_wind_setting = 0
   end function



! set_RGB_wind_setting:
! Set the current setting for mass-loss (RGB)
   integer function set_RGB_wind_setting(value)
      use real_kind
      use massloss
      implicit none
      real(double), intent(in) :: value
      if (value .ge. 0.0) then
        multiplier_schroeder = value
        multiplier_reimers = 0.0
      else
        multiplier_schroeder = 0.0
        multiplier_reimers = -value
      endif
      set_RGB_wind_setting = 0
   end function



! get_Ostar_wind_setting:
! Get the current setting for mass-loss (O/B stars)
   integer function get_Ostar_wind_setting(value)
      use real_kind
      use massloss
      implicit none
      real(double), intent(out) :: value
      value = multiplier_vink
      get_Ostar_wind_setting = 0
   end function

! set_Ostar_wind_setting:
! Set the current setting for mass-loss (O/B stars)
   integer function set_Ostar_wind_setting(value)
      use real_kind
      use massloss
      implicit none
      real(double), intent(in) :: value
      multiplier_vink = value
      set_Ostar_wind_setting = 0
   end function



   ! Return the stellar type of the specified star
   integer function stellar_type_of(star_id)
      implicit none
      integer, intent(in) :: star_id
      integer :: tmp
      tmp = get_stellar_type(star_id, stellar_type_of)
   end function stellar_type_of
   integer function get_stellar_type(star_id, stellar_type)
      use real_kind
      implicit none
      integer, intent(in) :: star_id
      integer, intent(out) :: stellar_type
      integer, external :: find_stellar_type
      get_stellar_type = -1
      stellar_type = -1

      if (star_id_is_valid(star_id) < 0) return
      call select_star(star_id)
      stellar_type = find_stellar_type()
      get_stellar_type = 0
   end function get_stellar_type




   integer function set_ev_path(new_ev_path)
      use file_exists_module
      implicit none
      character(len=*), intent(in) :: new_ev_path
      if (.not. file_exists(new_ev_path) ) then
         if (.not. file_exists(trim(new_ev_path)//'/input/amuse/init.dat') ) then
            if (verbose) print *, "Warning: file ",trim(new_ev_path)," for ", trim(amuse_ev_path), " does not exist!"
            set_ev_path = -1
            return
         end if
      end if
      amuse_ev_path = new_ev_path
      set_ev_path = 0
   end function

   ! Part of standard AMUSE interface, called when all parameters have been set
   integer function commit_parameters()
      use stopping_conditions
      use settings, only: csmc, calp, cos, cps, cth, kion
      implicit none
      commit_parameters = initialise_twin(amuse_ev_path, amuse_nstars, amuse_Z, &
         amuse_verbose, amuse_nmesh)
      csmc = amuse_csmc
      calp = amuse_calp
      cos = amuse_cos
      cps = amuse_cos
      cth = amuse_cth
      kion = amuse_kion
      uc(2) = amuse_maxage
      uc(12) = amuse_mindt
   end function commit_parameters
   integer function recommit_parameters()
      use stopping_conditions
      use settings, only: csmc, calp, cos, cps, cth, kion
      implicit none
      csmc = amuse_csmc
      calp = amuse_calp
      cos = amuse_cos
      cps = amuse_cos
      cth = amuse_cth
      kion = amuse_kion
      uc(2) = amuse_maxage
      uc(12) = amuse_mindt
      recommit_parameters = 0
   end function recommit_parameters
   integer function commit_particles()
      implicit none
      commit_particles = 0
   end function commit_particles
   integer function recommit_particles()
      implicit none
      recommit_particles = 0
   end function recommit_particles




! TODO: need to implement these:


! Return the maximum_number_of_stars parameter
   integer function get_maximum_number_of_stars(value)
      implicit none
      integer, intent(out) :: value
      value = amuse_nstars
      get_maximum_number_of_stars = 0
   end function
! Set the maximum_number_of_stars parameter
   integer function set_maximum_number_of_stars(value)
      implicit none
      integer, intent(in) :: value
      amuse_nstars = value
      set_maximum_number_of_stars = 0
   end function

! Retrieve the current value of the mixing length ratio
   integer function get_mixing_length_ratio(value)
      implicit none
      real(double), intent(out) :: value
      value = amuse_calp
      get_mixing_length_ratio = 0
   end function
! Set the current value of the mixing length ratio
   integer function set_mixing_length_ratio(value)
      implicit none
      real(double), intent(in) :: value
      amuse_calp = value
      set_mixing_length_ratio = 0
   end function

! Return the minimum timestep stop condition
   integer function get_min_timestep_stop_condition(value)
      implicit none
      real(double), intent(out) :: value
      value = amuse_mindt
      get_min_timestep_stop_condition = 0
   end function
! Set the minimum timestep stop condition
   integer function set_min_timestep_stop_condition(value)
      implicit none
      real(double), intent(in) :: value
      amuse_mindt = value
      set_min_timestep_stop_condition = 0
   end function

! Return the maximum age stop condition
   integer function get_max_age_stop_condition(value)
      implicit none
      real(double), intent(out) :: value
      value = amuse_maxage
      get_max_age_stop_condition = 0
   end function
! Set the maximum age stop condition
   integer function set_max_age_stop_condition(value)
      implicit none
      real(double), intent(in) :: value
      amuse_maxage = value
      set_max_age_stop_condition = 0
   end function

   integer function get_convection_efficiency(value)
      use settings, only: crd
      implicit none
      real(double), intent(out) :: value
      value = crd
      get_convection_efficiency = 0
   end function
   integer function set_convection_efficiency(value)
      use settings, only: crd
      implicit none
      real(double), intent(in) :: value
      crd = value
      set_convection_efficiency = 0
   end function

! Retrieve the current value of the efficiency of semi-convection
   integer function get_semi_convection_efficiency(value)
      implicit none
      real(double), intent(out) :: value
      value = amuse_csmc
      get_semi_convection_efficiency = 0
   end function
! Set the current value of the efficiency of semi-convection
   integer function set_semi_convection_efficiency(value)
      implicit none
      real(double), intent(in) :: value
      amuse_csmc = value
      set_semi_convection_efficiency = 0
   end function

   integer function get_thermohaline_efficiency(value)
      implicit none
      real(double), intent(out) :: value
      value = amuse_cth
      get_thermohaline_efficiency = 0
   end function
   integer function set_thermohaline_efficiency(value)
      implicit none
      real(double), intent(in) :: value
      amuse_cth = value
      set_thermohaline_efficiency = 0
   end function

! Retrieve the current number of elements used for ionization in the EoS
   integer function get_number_of_ionization_elements(value)
      implicit none
      integer, intent(out) :: value
      value = amuse_kion
      get_number_of_ionization_elements = 0
   end function
! Set the current number of elements used for ionization in the EoS
   integer function set_number_of_ionization_elements(value)
      implicit none
      integer, intent(in) :: value
      amuse_kion = value
      set_number_of_ionization_elements = 0
   end function

! Retrieve the current value of the convective overshoot parameter
   integer function get_convective_overshoot_parameter(value)
      implicit none
      real(double), intent(out) :: value
      value = amuse_cos
      get_convective_overshoot_parameter = 0
   end function
! Set the current value of the convective overshoot parameter
   integer function set_convective_overshoot_parameter(value)
      implicit none
      real(double), intent(in) :: value
      amuse_cos = value
      set_convective_overshoot_parameter = 0
   end function


   integer function get_metallicity(value)
      implicit none
      real(double), intent(out) :: value
      value = amuse_Z
      get_metallicity = 0
   end function
   integer function set_metallicity(value)
      implicit none
      real(double), intent(in) :: value
      if (initialised) then
         ! Can only be done once, on startup
         set_metallicity = -1
      else
         amuse_Z = value
         set_metallicity = 0
      end if
   end function

   integer function get_verbosity(value)
      implicit none
      integer, intent(out) :: value
      if (amuse_verbose) then
         value = 1
      else
         value = 0
      end if
      get_verbosity = 0
   end function get_verbosity
   integer function set_verbosity(value)
      implicit none
      integer, intent(in) :: value
      set_verbosity = 0
      if (value .eq. 0) then
         amuse_verbose = .false.
      else if (value .eq. 1) then
         amuse_verbose = .true.
      else
         set_verbosity = -1
      end if
   end function set_verbosity

   integer function get_import_model_entropy_accuracy(value)
      implicit none
      real(double), intent(out) :: value
      value = amuse_entropy_accuracy
      get_import_model_entropy_accuracy = 0
   end function
   integer function set_import_model_entropy_accuracy(value)
      implicit none
      real(double), intent(in) :: value
      amuse_entropy_accuracy = value
      set_import_model_entropy_accuracy = 0
   end function

   integer function get_import_model_entropy_force(value)
      implicit none
      real(double), intent(out) :: value
      value = amuse_entropy_force
      get_import_model_entropy_force = 0
   end function
   integer function set_import_model_entropy_force(value)
      implicit none
      real(double), intent(in) :: value
      amuse_entropy_force = value
      set_import_model_entropy_force = 0
   end function

   integer function initialize_code()
      implicit none
      amuse_ev_path = 'src/trunk'
      amuse_nstars = 10
      amuse_nmesh = 500
      amuse_verbose = .true.
      amuse_Z = 0.02d0
      amuse_kion = 5
      amuse_cth = 1.0d0
      amuse_cos = 0.12d0
      amuse_csmc = 0.04d0
      amuse_calp = 2.0d0
      amuse_mindt = 1.0d6
      amuse_maxage = 2.0d12
      amuse_entropy_accuracy = 1.0d-4
      amuse_entropy_force = 20.0d0
      initialize_code = 0
   end function

   integer function cleanup_code()
      implicit none
      cleanup_code = 0
   end function

   ! This code assumes that delta_t is in years
   integer function evolve_for(star_id, delta_t)
      implicit none
      integer, intent(in) :: star_id
      real(double), intent(in) :: delta_t
      evolve_for = evolve_until_model_time(star_id, age_of(star_id) + delta_t)
   end function
   integer function evolve_one_step(star_id)
      implicit none
      integer, intent(in) :: star_id
      evolve_one_step = evolve_one_timestep(star_id)
   end function

   integer function get_wind_mass_loss_rate(star_id, value)
      use constants
      implicit none
      integer, intent(in) :: star_id
      real(double), intent(out) :: value
      get_wind_mass_loss_rate = -1
      if (star_id_is_valid(star_id) < 0) return
      value = star_list(star_id)%zet(1)*csy/cmsn
      get_wind_mass_loss_rate = 0
   end function

   ! msun / yr
   function get_mass_transfer_rate(star_id, value)
      use constants
      implicit none
      integer :: get_mass_transfer_rate
      real(double) , intent(out) :: value
      integer, intent(in) :: star_id
      get_mass_transfer_rate = -1
      if (star_id_is_valid(star_id) < 0) return
      value = star_list(star_id)%xit(1)*csy/cmsn
      get_mass_transfer_rate = 0
   end function

! Returns the star's spin period, in days
   integer function get_spin(id, value)
      implicit none
      real(double), intent(out) :: value
      integer, intent(in) :: id
      get_spin = -1
      if (star_id_is_valid(id) < 0) return
      value = star_list(id)%p
      get_spin = 0
   end function

   ! Binding energy of the stellar envelope [erg/(1Mo)]
   integer function get_envelope_binding_energy(id, value)
      implicit none
      integer, intent(in)       :: id
      real(double), intent(out) :: value
      get_envelope_binding_energy = -1
      if (star_id_is_valid(id) < 0) return
      value = star_list(id)%be(1)
      get_envelope_binding_energy = 0
   end function get_envelope_binding_energy

   ! Binding energy of the stellar envelope: gravity [erg/(1Mo)]
   integer function get_envelope_gravitational_energy(id, value)
      implicit none
      integer, intent(in)       :: id
      real(double), intent(out) :: value
      get_envelope_gravitational_energy = -1
      if (star_id_is_valid(id) < 0) return
      value = star_list(id)%be(1)
      get_envelope_gravitational_energy = 0
   end function get_envelope_gravitational_energy

   ! Binding energy of the stellar envelope: internal energy [erg/(1Mo)]
   integer function get_envelope_internal_energy(id, value)
      implicit none
      integer, intent(in)       :: id
      real(double), intent(out) :: value
      get_envelope_internal_energy = -1
      if (star_id_is_valid(id) < 0) return
      value = star_list(id)%be(1)
      get_envelope_internal_energy = 0
   end function get_envelope_internal_energy

   ! Binding energy of the stellar envelope: recombination energy [erg/(1Mo)]
   integer function get_envelope_recombination_energy(id, value)
      implicit none
      integer, intent(in)       :: id
      real(double), intent(out) :: value
      get_envelope_recombination_energy = -1
      if (star_id_is_valid(id) < 0) return
      value = star_list(id)%be(1)
      get_envelope_recombination_energy = 0
   end function get_envelope_recombination_energy

   ! Binding energy of the stellar envelope: H2 association energy [erg/(1Mo)]
   integer function get_envelope_association_energy(id, value)
      implicit none
      integer, intent(in)       :: id
      real(double), intent(out) :: value
      get_envelope_association_energy = -1
      if (star_id_is_valid(id) < 0) return
      value = star_list(id)%be(1)
      get_envelope_association_energy = 0
   end function get_envelope_association_energy

   integer function get_number_of_particles(value)
      implicit none
      integer, intent(out) :: value
      integer :: i
      value = 0
      do i = 1, max_stars
         if (star_list(i)%exists) value = value+1
      end do
      get_number_of_particles = 0
   end function

   integer function star_id_is_valid(star_id)
      implicit none
      integer, intent(in) :: star_id
      star_id_is_valid = -1
      ! nb: it looks like we can combine the following two lines into one, but we cannot since fortran does not
      ! short-circuit boolean logic and all elements in the expression are evaluated. so we must only index the
      ! array *after* we have established that the index is in range.
      if (star_id<1 .or. star_id>max_stars) return
      if (.not. star_list(star_id)%exists)  return
      star_id_is_valid = 0
   end function star_id_is_valid

! Internal structure getters:

! Return the current number of zones/mesh-cells of the star
      integer function get_number_of_zones(star_id, value)
         implicit none
         integer, intent(in) :: star_id
         integer, intent(out) :: value
         if (star_id_is_valid(star_id) < 0) then
            value = 0
            get_number_of_zones = -21
            return
         end if
         value = star_list(star_id)%number_of_meshpoints
         get_number_of_zones = 0
      end function

! Update the stellar structure data sx for the specified star, if necessary
      subroutine update_quantities_if_needed(star_id)
         implicit none
         integer, intent(in) :: star_id
         if (sx_updated_for_star .ne. star_id) then
            call swap_in(star_id)
            call compute_output_quantities ( 1 )
            sx_updated_for_star = star_id
         end if
      end subroutine

! Return the requested quantity (index in SX array) at the specified zone/mesh-cell of the star
      integer function get_quantity_at_zone(star_id, zone, qty, value)
         use structure_variables
         implicit none
         integer, intent(in) :: star_id, zone, qty
         double precision, intent(out) :: value
         if (star_id_is_valid(star_id) < 0) then
            value = -1.0
            get_quantity_at_zone = -21
            return
         end if
         if (zone >= star_list(star_id)%number_of_meshpoints .or. zone < 0) then
            value = -1.0
            get_quantity_at_zone = -22
            return
         end if
         call update_quantities_if_needed(star_id)
         value = sx(qty, zone+2)
         get_quantity_at_zone = 0
      end function get_quantity_at_zone

! Return the temperature at the specified zone/mesh-cell of the star
      integer function get_temperature_at_zone(star_id, zone, value)
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_temperature_at_zone = get_quantity_at_zone(star_id, zone, 4, value)
      end function

! Return the density at the specified zone/mesh-cell of the star
      integer function get_density_at_zone(star_id, zone, value)
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_density_at_zone = get_quantity_at_zone(star_id, zone, 3, value)
      end function

! Return the density at the specified zone/mesh-cell of the star
      integer function get_pressure_at_zone(star_id, zone, value)
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_pressure_at_zone = get_quantity_at_zone(star_id, zone, 2, value)
      end function

! Return the entropy at the specified zone/mesh-cell of the star
      integer function get_entropy_at_zone(star_id, zone, value)
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_entropy_at_zone = get_quantity_at_zone(star_id, zone, 28, value)
      end function

! Return the internal energy at the specified zone/mesh-cell of the star
      integer function get_uint_at_zone(star_id, zone, value)
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_uint_at_zone = get_quantity_at_zone(star_id, zone, 27, value)
      end function

! Return the electron degeneracy parameter at the specified zone/mesh-cell of the star
      integer function get_psie_at_zone(star_id, zone, value)
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_psie_at_zone = get_quantity_at_zone(star_id, zone, 1, value)
      end function

! Return the opacity at the specified zone/mesh-cell of the star
      integer function get_opacity_at_zone(star_id, zone, value)
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_opacity_at_zone = get_quantity_at_zone(star_id, zone, 5, value)
      end function

! Return the adiabatic gamma (gamma1) at the specified zone/mesh-cell of the star
      integer function get_gamma_ad_at_zone(star_id, zone, value)
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_gamma_ad_at_zone = get_quantity_at_zone(star_id, zone, 62, value)
      end function

! Return the specific heat at constant pressure at the specified zone/mesh-cell of the star
      integer function get_scp_at_zone(star_id, zone, value)
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_scp_at_zone = get_quantity_at_zone(star_id, zone, 63, value)
      end function

! Return the adiabatic temperature gradient at the specified zone/mesh-cell of the star
      integer function get_nabla_ad_at_zone(star_id, zone, value)
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_nabla_ad_at_zone = get_quantity_at_zone(star_id, zone, 6, value)
      end function

! Return the radiative temperature gradient at the specified zone/mesh-cell of the star
      integer function get_nabla_rad_at_zone(star_id, zone, value)
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_nabla_rad_at_zone = get_quantity_at_zone(star_id, zone, 8, value) + get_quantity_at_zone(star_id, zone, 6, value)
      end function

! Return the actual temperature gradient at the specified zone/mesh-cell of the star
      integer function get_nabla_at_zone(star_id, zone, value)
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_nabla_at_zone = get_quantity_at_zone(star_id, zone, 7, value)
      end function

! Return the radius at the specified zone/mesh-cell of the star
      integer function get_radius_at_zone(star_id, zone, value)
         use constants
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_radius_at_zone = get_quantity_at_zone(star_id, zone, 17, value)
         value = value * CRSN * 1.0D11
      end function

! Return the luminosity at the specified zone/mesh-cell of the star
      integer function get_luminosity_at_zone(star_id, zone, value)
         use constants
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_luminosity_at_zone = get_quantity_at_zone(star_id, zone, 18, value)
         value = value * CLSN * 1.0D33
      end function

! Return the mass-coordinate at the specified zone/mesh-cell of the star
      integer function get_mass_at_zone(star_id, zone, value)
         use constants
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_mass_at_zone = get_quantity_at_zone(star_id, zone, 9, value)
         value = value * CMSN * 1.0D33
      end function

! Return the mass in the specified zone/mesh-cell of the star
      integer function get_mass_in_zone(star_id, zone, value)
         use constants
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_mass_in_zone = get_quantity_at_zone(star_id, zone, 22, value)
         value = value * CMSN * 1.0D33
      end function

! Return the thermal energy production rate at the specified zone/mesh-cell of the star
      integer function get_eth_at_zone(star_id, zone, value)
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_eth_at_zone = get_quantity_at_zone(star_id, zone, 19, value)
      end function

! Return the nuclear energy generation rate at the specified zone/mesh-cell of the star
      integer function get_enuc_at_zone(star_id, zone, value)
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_enuc_at_zone = get_quantity_at_zone(star_id, zone, 20, value)
      end function

! Return the neutrino energy loss rate at the specified zone/mesh-cell of the star
      integer function get_nuloss_at_zone(star_id, zone, value)
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_nuloss_at_zone = get_quantity_at_zone(star_id, zone, 21, value)
      end function

! Return the mean molecular weight per particle (ions + free electrons) at the specified zone/mesh-cell of the star
      integer function get_mu_at_zone(star_id, zone, value)
         implicit none
         integer, intent(in) :: star_id, zone
         double precision, intent(out) :: value
         get_mu_at_zone = get_quantity_at_zone(star_id, zone, 31, value)
      end function

! Return the current number of chemical abundance variables per zone of the star
   integer function get_number_of_species(star_id, value)
      implicit none
      integer, intent(in) :: star_id
      integer, intent(out) :: value
      value = 9
      get_number_of_species = 0
   end function

! Return the name of chemical abundance variable 'AMUSE_species' of the star
   integer function get_name_of_species(star_id, AMUSE_species, value)
      implicit none
      integer, intent(in) :: star_id, AMUSE_species
      character (len=6), intent(out) :: value
      get_name_of_species = 0
      select case (AMUSE_species)
         case (1)
            value = 'h1'
         case (2)
            value = 'he4'
         case (3)
            value = 'c12'
         case (4)
            value = 'n14'
         case (5)
            value = 'o16'
         case (6)
            value = 'ne20'
         case (7)
            value = 'mg24'
         case (8)
            value = 'si28'
         case (9)
            value = 'fe56'
         case default
            value = 'error'
            get_name_of_species = -23
      end select
   end function

! Return the mass fraction of species 'AMUSE_species' at the specified
! zone/mesh-cell of the star
   integer function get_mass_fraction_of_species_at_zone(star_id, &
         AMUSE_species, zone, value)
      use atomic_data
!~      use extra_elements
      use binary_history, only: hpr
      use indices
      use mesh
      implicit none
      integer, intent(in) :: star_id, zone, AMUSE_species
      double precision, intent(out) :: value
      real(double) :: xa(9), na(9)
      real(double) :: avm
      integer :: zone_index, i

      value = -1.0
      if (star_id_is_valid(star_id) < 0) then
         get_mass_fraction_of_species_at_zone = -21
         return
      end if
      if (zone >= star_list(star_id)% number_of_meshpoints .or. zone < 0) then
         get_mass_fraction_of_species_at_zone = -22
         return
      end if
      if (AMUSE_species > 9 .or. AMUSE_species < 1) then
         get_mass_fraction_of_species_at_zone = -23
         return
      end if
      call update_quantities_if_needed(star_id)
      zone_index = star_list(star_id)%number_of_meshpoints - zone
      xa(1) = h(VAR_H1, zone_index)
      xa(2) = h(VAR_HE4, zone_index)
      xa(3) = h(VAR_C12, zone_index)
      xa(4) = h(VAR_N14, zone_index)
      xa(5) = h(VAR_O16, zone_index)
      xa(6) = h(VAR_NE20, zone_index)
      xa(7) = h(VAR_MG24, zone_index)
      xa(8) = h(VAR_SI28, zone_index)
      xa(9) = h(VAR_FE56, zone_index)
      do i=1, 9
        na(i) = xa(i) * can(i)/cbn(i)
      end do
      avm = sum(na(1:9))
      value = na(AMUSE_species) / avm
      get_mass_fraction_of_species_at_zone = 0
   end function

! Return the internal structure of the star at a specific zone
   integer function get_stellar_model_element(zone, star_id, &
         d_mass, cumul_mass, radius, rho, pressure, entropy, temperature, &
         luminosity, molecular_weight, XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE)
      use mesh
      use structure_variables
      use atomic_data
!~      use extra_elements
      use binary_history, only: hpr
      use indices

      implicit none
      integer, intent(in) :: star_id, zone
      double precision, intent(out) :: d_mass, cumul_mass, radius, rho, &
         pressure, entropy, temperature, luminosity, molecular_weight, &
         XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE
      real(double) :: xa(9), na(9)
      real(double) :: avm
      integer :: zone_index, i

      if (star_id_is_valid(star_id) < 0) then
         get_stellar_model_element = -21
         return
      end if
      if (zone > star_list(star_id)% number_of_meshpoints .or. zone < 1) then
         get_stellar_model_element = -22
         return
      end if

      call update_quantities_if_needed(star_id)
      d_mass      = sx(22, 1 + zone)
      cumul_mass  = sx(9,  1 + zone)
      radius      = sx(17, 1 + zone)
      rho         = sx(3,  1 + zone)
      pressure    = sx(2,  1 + zone)

      entropy     = sx(28, 1 + zone)
      temperature = sx(4,  1 + zone)
      luminosity  = sx(18, 1 + zone)
      molecular_weight = sx(31, 1 + zone)

      ! Convert *all* abundances to mass fractions
      zone_index = star_list(star_id)%number_of_meshpoints - zone + 1
      xa(1) = h(VAR_H1, zone_index)
      xa(2) = h(VAR_HE4, zone_index)
      xa(3) = h(VAR_C12, zone_index)
      xa(4) = h(VAR_N14, zone_index)
      xa(5) = h(VAR_O16, zone_index)
      xa(6) = h(VAR_NE20, zone_index)
      xa(7) = h(VAR_MG24, zone_index)
      xa(8) = h(VAR_SI28, zone_index)
      xa(9) = h(VAR_FE56, zone_index)
      do i=1, 9
        na(i) = xa(i) * can(i)/cbn(i)
      end do
      avm = sum(na(1:9))

      XH = na(1) / avm
      XHE = na(2) / avm
      XC = na(3) / avm
      XN = na(4) / avm
      XO = na(5) / avm
      XNE = na(6) / avm
      XMG = na(7) / avm
      XSI = na(8) / avm
      XFE = na(9) / avm

      get_stellar_model_element = 0
   end function

   ! Inject energy into the star at the given zone.
   ! The value is in erg/g/s
   integer function set_extra_energy_at_zone(star_id, zone, value)
      use mesh_enc, only: menc
      implicit none
      integer, intent(in) :: star_id, zone
      double precision, intent(in) :: value
      integer :: zone_index
      if (star_id_is_valid(star_id) < 0) then
         set_extra_energy_at_zone = -21
         return
      end if
      if (zone >= star_list(star_id)%number_of_meshpoints .or. zone < 0) then
         set_extra_energy_at_zone = -22
         return
      end if
      set_extra_energy_at_zone = 0
      call select_star(star_id)
      zone_index = star_list(star_id)%number_of_meshpoints - zone
      menc(1, zone_index) = value
   end function set_extra_energy_at_zone

   integer function set_mass_fraction_of_species_at_zone(star_id, &
         AMUSE_species, zone, value)
      implicit none
      integer, intent(in) :: star_id, zone, AMUSE_species
      double precision, intent(in) :: value
      set_mass_fraction_of_species_at_zone = -4
   end function

   integer function set_density_at_zone(star_id, zone, value)
      implicit none
      integer, intent(in) :: star_id, zone
      double precision, intent(in) :: value
      set_density_at_zone = -4
   end function

   integer function set_temperature_at_zone(star_id, zone, value)
      implicit none
      integer, intent(in) :: star_id, zone
      double precision, intent(in) :: value
      set_temperature_at_zone = -4
   end function

   integer function set_radius_at_zone(star_id, zone, value)
      implicit none
      integer, intent(in) :: star_id, zone
      double precision, intent(in) :: value
      set_radius_at_zone = -4
   end function

   integer function set_mass(star_id, value)
      implicit none
      integer, intent(in) :: star_id
      double precision, intent(in) :: value
      set_mass = -4
   end function


end module twinlib
