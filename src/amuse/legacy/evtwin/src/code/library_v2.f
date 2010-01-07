! Library interface to TWIN: useful for putting TWIN into Muse
! TODO:
!  * Replace all IMPLICIT statements with IMPLICIT NONE
      module twin_library_v2
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
         double precision, pointer :: h(:,:)
         double precision, pointer :: dh(:,:)
         
         double precision :: zams_mass

         ! Stellar type, as in Hurley & al 2000
         integer :: stellar_type
         
         ! Iteration control
         integer :: startup_iter
         integer :: normal_iter
         
         ! Binary parameters; these need to be set because TWIN is a binary code
         double precision :: bms, per, ecc, p

         ! Don't ask...
         integer :: jf
         
         ! COMMON block QUERY
         double precision :: ml, ql, xl, uc(21)
         integer :: jmod, jnn, jter, joc, jkh
         
         ! COMMON block TVBLES
         double precision :: dt, zq(9), age, bm(71), pr(81), ppr(81)
         integer :: jhold, jm2, jm1
         
         ! COMMON block SOLV
         double precision :: er(nvar)

         ! COMMON block STORE
         double precision, pointer :: hpr(:, :)
         !double precision :: ht(4, nvar)    ! Used in nucleosynthesis code
         !double precision :: ms(9999)    ! Mass loss history of primary, non-TWIN mode
         !double precision :: st(9999)    ! Mass loss history of primary, non-TWIN mode
      end type
      
      ! Inverse permutation: where each element in the global H is
      ! mapped to in the copy in the struct. This makes it possible to
      ! not store unused elements, reducing the memory requirements.
      integer, private :: INV_VAR_PERM(NVAR)
      integer, private :: VAR_PERM(NVAR)
      integer, private :: actual_number_of_variables

      ! Data structure to store a number of stars in memory.
      integer, private :: max_stars = -1        ! Maximum number of stars
      type(twin_star_t), private, allocatable, target :: star_list(:)
      integer, private :: current_star = 0      ! Currently selected star
      integer, private :: num_stars = 0         ! Number of allocated stars

      ! Some data read from init.dat that is not saved to COMMON blocks
      integer, private :: KH2, KR1, KR2, KSV, KT5, JCH, iter
      
      ! Layout of the ZAMS library file
      double precision, private :: MLO, DM, MHI
      integer, private :: KDM 
      
      ! Number of models to run before switching to "normal" number of
      ! iterations (from "startup" number of iterations)
      integer, parameter, private :: switch_iterations = 5

      ! Print verbose output to stdout yes or no.
      logical, private :: verbose = .false.
      
      ! Name of the init.dat input file, if not otherwise specified
      character*500, private :: init_dat_name = 'init.dat'
      
      ! Name of the init.run input file, if not otherwise specified
      character*500, private :: init_run_name = 'init.run'
      
      
      character(len = 500), private :: ev_path = 'src'
      character*100, private :: metallicity_str = '02'
      integer, private :: maximum_number_of_stars = 10
      contains

! set_init_dat_name:
!  Change the name of the input file used for the numerical settings and
!  physics options (default: ./init.dat)
      function set_init_dat_name(new_init_dat_name)
      use file_exists_module
      implicit none
      integer :: set_init_dat_name
      character(len=*), intent(in) :: new_init_dat_name;

      IF (.NOT. FILE_EXISTS(new_init_dat_name) ) THEN
         IF (VERBOSE) PRINT *, "Warning: file ",TRIM(new_init_dat_name)," for ", TRIM(init_dat_name), " does not exist!"
         set_init_dat_name = -1
      ELSE
         set_init_dat_name = 0
      END IF
      init_dat_name = trim(new_init_dat_name)
      end function

! set_init_run_name:
!  Change the name of the input file used for the numerical settings and
!  physics options (default: ./init.dat)
      function set_init_run_name(new_init_run_name)
      use file_exists_module
      implicit none
      integer :: set_init_run_name
      character(len = *), intent(in) :: new_init_run_name;

      IF (.NOT. FILE_EXISTS(new_init_run_name) ) THEN
         IF (VERBOSE) PRINT *, "Warning: file ",TRIM(new_init_run_name)," for ", TRIM(init_run_name), " does not exist!"
         set_init_run_name = -1
      ELSE
         set_init_run_name = 0
      END IF
      init_run_name = trim(new_init_run_name);
      end function
      
      function set_ev_path(new_ev_path)
      use file_exists_module
      implicit none
      integer :: set_ev_path
      character(len=*), intent(in) :: new_ev_path;

      IF (.NOT. FILE_EXISTS(new_ev_path) ) THEN
         IF (VERBOSE) PRINT *, "Warning: file ",TRIM(new_ev_path)," for ", TRIM(ev_path), " does not exist!"
         set_ev_path = -1
      ELSE
         set_ev_path = 0
      END IF
      ev_path = new_ev_path;
      end function

      subroutine open_fortran_output_file(filename, file_number)
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(in) :: file_number
      integer :: status;
      close(file_number)
      open(unit = file_number, action = "write", file = filename, iostat=status)
      return
      end subroutine

      function read_atomic_data_file(filename)
      implicit none
      character(len=*), intent(in) :: filename
      integer :: read_atomic_data_file
      close (26)
      open(UNIT = 26, ACTION="READ", FILE=filename, IOSTAT=read_atomic_data_file)
      CALL LOAD_ATOMIC_DATA(26)
      close (26)
      end function

      function read_reaction_rate_file(filename)
      USE CONSTANTS
      implicit none
      character(len=*), intent(in) :: filename
      integer :: read_reaction_rate_file
      CALL INITIALISE_CONSTANTS
      close (42)
      open(UNIT = 42, ACTION="READ", FILE=filename, IOSTAT=read_reaction_rate_file)
      CALL LOAD_REACTION_AND_NEUTRINO_RATES(42)
      close (42)
      end function

      function read_opacity_file(filename)
      implicit none
      character(len=*), intent(in) :: filename
      integer :: read_opacity_file
      close (20)
      open(UNIT = 20, ACTION="READ", FILE=filename, IOSTAT=read_opacity_file)
      CALL LOAD_OPACITY(20)
      close (20)
      end function

      function read_co_opacity_file(filename)
      implicit none
      character(len=*), intent(in) :: filename
      integer :: read_co_opacity_file
      close (41)
      open(UNIT = 41, ACTION="READ", FILE=filename, IOSTAT=read_co_opacity_file)
      CALL LOAD_OPACITY_CO(41)
      close (41)
      end function

      function read_zams_library_format(filename)
      implicit none
      character(len=*), intent(in) :: filename
      integer :: read_zams_library_format
      close (19)
      open(UNIT = 19, ACTION="READ", FILE=filename, IOSTAT=read_zams_library_format)
      READ (19, *) MLO, DM, MHI, KDM 
      close (19)
      end function

      subroutine close_zams_library
      close(16)
      end subroutine

      function open_zams_library(filename)
      implicit none
      character(len=*), intent(in) :: filename
      integer :: open_zams_library
      close(16)
      open(UNIT = 16, ACTION="READ", FILE=filename, IOSTAT=open_zams_library)
      end function


      subroutine set_number_of_startup_iterations(new_iter)
      implicit none
      integer, intent(in) :: new_iter
      kr1 = new_iter
      iter = kr1
      end subroutine

      subroutine set_number_of_iterations(new_iter)
      implicit none
      integer, intent(in) :: new_iter
      kr2 = new_iter 
      end subroutine

      function get_number_of_iterations()
      implicit none
      integer :: get_number_of_iterations
      get_number_of_iterations = kr2
      end function


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
      USE SETTINGS
      USE EXTRA_ELEMENTS
      USE FILE_EXISTS_MODULE
      USE CONSTANTS
      IMPLICIT REAL*8 (A-H, L-Z)    ! Ugly hack, needs to be changed to none!
      integer :: initialise_twin
      CHARACTER(LEN=*), intent(in) ::  EVPATH, ZSTR
      integer, intent(in) :: nstars

      REAL*8 H1(NVAR,NM), H2(NVAR,NM)
      DOUBLE PRECISION :: SM0, DTY0, AGE0, PER0, BMS0, ECC0, P0, ENC0
      DOUBLE PRECISION :: SM, DTY, AGE, PER, BMS, ECC, P, ENC
      DOUBLE PRECISION :: SM1, DTY1, AGE1, PER1, BMS1, ECC1, P1, ENC1
      DOUBLE PRECISION :: SM2, DTY2, AGE2, PER2, BMS2, ECC2, P2, ENC2
      CHARACTER*500 :: STARTFILE
      INTEGER :: KH, KTW,
     &KE1,KE2,KE3,KBC,KEV,KFN,KL,JH(3),KP_VAR(40),KP_EQN(40),KP_BC(40),
     &KE1_2,KE2_2,KE3_2,KBC_2,KEV_2,KFN_2,KL_2,JH_2(3),KP_VAR_2(40),KP_EQN_2(40),KP_BC_2(40)
      COMMON H(NVAR,NM), DH(NVAR,NM), EP(3), KH, KTW,
     & KE1, KE2, KE3, KBC, KEV, KFN, KL, JH, KP_VAR, KP_EQN, KP_BC, 
     & KE1_2,KE2_2,KE3_2,KBC_2,KEV_2,KFN_2,KL_2,JH_2,KP_VAR_2,KP_EQN_2,KP_BC_2
      COMMON /TN1/ SM0, DTY0, AGE0, PER0, BMS0, ECC0, P0, ENC0
      COMMON /T0 / SM, DTY, AGE, PER, BMS, ECC, P, ENC
      COMMON /T1 / SM1, DTY1, AGE1, PER1, BMS1, ECC1, P1, ENC1
      COMMON /T2 / SM2, DTY2, AGE2, PER2, BMS2, ECC2, P2, ENC2
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /CINIT_RUN/ ISB, IP1, IM1, IP2, IM2, KPT, KP, ML1, DML, KML,
     & QL1, DQL, KQL, XL1, DXL, KXL, ROT, KR, EX, JMX, STARTFILE

      LOGICAL :: READ_INIT_DAT

      LOGICAL :: STATUS, OVERRIDE_RUN, OVERRIDE_DAT
      INTEGER, PARAMETER :: N_INP_FILES = 15
      CHARACTER*500 :: INPUTFILENAMES(N_INP_FILES), FNAME
      INTEGER :: INPUTUNITS(N_INP_FILES)     = (/12, 24, 16, 18, 19, 20, 21, 26, 63, 22, 23, 41, 42, 122, 123/)
      INTEGER :: INPUT_REQUIRED(N_INP_FILES) = (/ 0,  0,  1,  1,  1,  0,  1,  1, -1, -1, -1,  0,  1, 1, 1/)

      max_stars = nstars
      IF (VERBOSE) THEN
         PRINT *, 'TWIN initialisation.'
         PRINT *, 'Arguments: EVPATH =', EVPATH
         PRINT *, '           nstars =', nstars
         PRINT *, '           ZSTR =', ZSTR
      END IF

      ! Setup input file names
      INPUTFILENAMES(1)=TRIM(EVPATH)//"/input/zahb"//TRIM(ZSTR)//".mod"
      INPUTFILENAMES(2)=TRIM(EVPATH)//"/input/zahb"//".dat"
      INPUTFILENAMES(3)=TRIM(EVPATH)//"/input/zams/zams"//TRIM(ZSTR)//".mod"
      INPUTFILENAMES(4)=TRIM(EVPATH)//"/input/zams/zams"//TRIM(ZSTR)//".out"
      INPUTFILENAMES(5)=TRIM(EVPATH)//"/input/zams/zams"//TRIM(ZSTR)//".mas"
      INPUTFILENAMES(6)=TRIM(EVPATH)//"/input/metals/z"//TRIM(ZSTR)//"/phys.z"//TRIM(ZSTR)
      INPUTFILENAMES(7)=TRIM(EVPATH)//"/input/lt2ubv.dat"
      INPUTFILENAMES(8)=TRIM(EVPATH)//"/input/nucdata.dat"
      !INPUTFILENAMES(9)=TRIM(EVPATH)//"/input/mutate.dat"
      INPUTFILENAMES(10)=TRIM(init_dat_name)
      INPUTFILENAMES(11)=TRIM(init_run_name)
      INPUTFILENAMES(12)=TRIM(EVPATH)//"/input/COtables/COtables_z"//TRIM(ZSTR)
      INPUTFILENAMES(13)=TRIM(EVPATH)//"/input/physinfo.dat"
      INPUTFILENAMES(14)=TRIM(EVPATH)//"/run/muse/init.dat" ! Sensible defaults
      INPUTFILENAMES(15)=TRIM(EVPATH)//"/run/muse/init.run" ! Sensible defaults
      !INPUTFILENAMES(14)=TRIM(EVPATH)//"/xxx/init.dat" ! Sensible defaults
      !INPUTFILENAMES(15)=TRIM(EVPATH)//"/xxx/init.run" ! Sensible defaults

      
      ! We want to override the default init.run and init.dat settings if
      ! init.run and init.dat are present in the current directory.
      OVERRIDE_RUN = .FALSE.
      OVERRIDE_DAT = .FALSE.
      ! VERBOSE = 1
      ! Check if all input files exist and open them as needed
      DO I=1, N_INP_FILES
         WRITE (FNAME, '("fort.",I2)') INPUTUNITS(I)
         
         IF (FILE_EXISTS(INPUTFILENAMES(I)) .AND. .NOT. FILE_EXISTS(FNAME)) THEN
            IF (VERBOSE) PRINT *, "Opening ",TRIM(INPUTFILENAMES(I))," for ", TRIM(FNAME)
            OPEN(UNIT = INPUTUNITS(I), ACTION="READ", FILE=INPUTFILENAMES(I))
         END IF
         ! Check if files exist and act appropriately
         ! If the file is required (INPUT_REQUIRED>0), abort on error
         ! If the file is optional (INPUT_REQUIRED==0), give a warning
         ! If the file is probably unneeded (INPUT_REQUIRED<0), do nothing
         IF (.NOT. (FILE_EXISTS(INPUTFILENAMES(I)) .OR. FILE_EXISTS(FNAME))) THEN
            IF (INPUT_REQUIRED(I) > 0) THEN
               WRITE (0, *) 'Required input file ', TRIM(INPUTFILENAMES(I)), ' (', TRIM(FNAME), ') not found'
               initialise_twin = -1
               return
            ELSEIF (INPUT_REQUIRED(I) == 0) THEN
               WRITE (0, *) 'Warning: input file ', TRIM(INPUTFILENAMES(I)), ' (', TRIM(FNAME), ') not found'
            END IF
         END IF
      END DO

      IF ( FILE_EXISTS(INPUTFILENAMES(10)) ) OVERRIDE_DAT = .TRUE.
      IF ( FILE_EXISTS(INPUTFILENAMES(11)) ) OVERRIDE_RUN = .TRUE.
      
! We need to have a fort.11 - no, we don't, since we're not using star12!
      !WRITE (11, *) 0
      !CLOSE (11)

! We need a scratch file for output we don't want or need
      IF (.NOT. FILE_EXISTS('fort.25'))
     &   OPEN (UNIT=25, ACTION='WRITE', STATUS = 'SCRATCH')

! Other redundant (for our purposes) output we do not want to save but cannot
! easily avoid: solver output to fort.1 and opacity table output to fort.10
      OPEN (UNIT=10, ACTION='WRITE', STATUS = 'SCRATCH')
      !OPEN (UNIT=1, ACTION='WRITE', STATUS = 'SCRATCH')

      IF (VERBOSE) print *, 'Files initialised, initialising equation of state'
      CALL SETSUP
      
! Read layout of the ZAMS library
      READ (19, *) MLO, DM, MHI, KDM 
      REWIND (19)

      IF (VERBOSE) print *, 'Initialised equation of state, reading run settings'
      CALL READ_INIT_RUN(123)    ! Load defaults
      IF (OVERRIDE_RUN) CALL READ_INIT_RUN(23)
c Create short (`pruned') summary of ZAMS models from long file
! Are we using this? Not reall, I think...
!      CALL PRUNER ( 18, 17, ISB )
!      CLOSE (17)
!      CLOSE (18)

! Read init.dat file
! FIXME: this has an unwanted side effect of resetting options that are not set
! in the override file to their defaults, even if the default file we loaded 
! has them set. In other words, we need an "override init.dat" function.
      IF (VERBOSE) print *, 'Reading numerical and physical settings'
      STATUS = READ_INIT_DAT(122, KH2, KR1, KR2, KSV, KT5, JCH)   ! Defaults
      IF (STATUS .AND. OVERRIDE_DAT)
     &   STATUS = READ_INIT_DAT(22, KH2, KR1, KR2, KSV, KT5, JCH)
      IF (STATUS .EQV. .FALSE.) THEN
         initialise_twin = -1
         return
      END IF
      IF (VERBOSE) print *, 'Read settings'
      IF (VERBOSE) print *, 'Using', KH2, 'meshpoints per star'

!     Autodetect if we should solve for N14 or not by checking if the
!     corresponding equation is in the list of equations
      !USE_N14_EQN = .FALSE.
      !DO II = 1, 40
      !   IF (KP_EQN(II) == EN14) THEN
      !      USE_N14_EQN = .TRUE.
      !      EXIT
      !   ENDIF
      !   IF (KP_EQN(II) == 0) EXIT   ! Break loop if end of list found
      !ENDDO
!     Autodetect if we should solve for Mg24 or not by checking if the
!     corresponding equation is in the list of equations
      USE_MG24_EQN = .FALSE.
      DO II = 1, 40
         IF (KP_EQN(II) == EMG24 .OR. KP_EQN(II) == ESUMX) THEN
            USE_MG24_EQN = .TRUE.
            EXIT
         ENDIF
         IF (KP_EQN(II) == 0) EXIT   ! Break loop if end of list found
      ENDDO
      
C Convert some things to `cgs' units: 10**11 cm, 10**33 gm, 10**33 erg/s
      CMI = CMI/CSY
      !CMJ = CMJ*CMSN/CSY
      CMS = CMS/CSY
      CMT = CMT*1.0D-11

! Read opacity data and construct splines 
! KOP (read from init.dat) sets which type of opacity tables
      IF (VERBOSE) print *, 'Reading opacity tables'
      IF (KOP<=1) THEN
         IF (VERBOSE) print *, "Using OPAL '92"
         ! Iglesias & Rogers (1992), as implemented by Pols & al. 1995
         CALL LOAD_OPACITY(20)
      ELSE
         IF (VERBOSE) print *, "Using OPAL '96 with CO enhancement"
         ! Iglesias & Rogers (1996), as implemented by Eldridge & Tout 2003
         CALL LOAD_OPACITY_CO(41)
      ENDIF

      ! Abort if the requested number of meshpoints is larger than the
      ! size of the array
      IF (KH2 > NM) THEN
         WRITE (0, *) 'Cannot rezone to ', KH2, 'meshpoints. Maximum size is ', NM, 'meshpoints.'
         initialise_twin = -2
         RETURN
      END IF

! Initialise some more variables      
      IF ( ISB.EQ.1 ) KTW = 1
      JB = 1

! Setup inverse lookup table for number of variables
      actual_number_of_variables = KE1+KE2+KEV
      if (verbose) print *,'Will backup',actual_number_of_variables,'variables.'
      INV_VAR_PERM(:) = 0
      DO II=1, actual_number_of_variables
         INV_VAR_PERM(KP_VAR(II)) = II
         VAR_PERM(II) = KP_VAR(II)
      END DO

! Allocate memory for all the stars we're interested in
      if (allocated(star_list)) deallocate(star_list);
      allocate(star_list(1:max_stars))
      if (verbose) print *, 'Allocated memory for ',max_stars, 'stars'

! We're not writing any output directly from the library functions
      ! Report success
      initialise_twin = 0;
      end function
      
      
      
! Close the temporary files when we're done with them
! Sometimes we may need to do this manually
      subroutine close_temporary_files
      implicit none
      close (25)
      close (10)
      end subroutine



! load_zams_star:
!  load a ZAMS model from disk
!  Input variables:
!   m - the mass, in solar units
!   t - the starting age, in years. Only used for reporting the star's age
!  Return value:
!   >0: The stars ID for identifying it in the array of models
!   =0: No star allocated, out of memory
      function load_zams_star(mass, age)
      use constants
      use settings
      IMPLICIT REAL*8 (A-H, L-Z)    ! Ugly hack, needs to be changed to none!
      type(twin_star_t), pointer :: star
      integer :: load_zams_star
      double precision, intent(in) :: mass, age
      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, ID(130), IE(130) 
      COMMON /T1 / SM1, DTY1, AGE1, PER1, BMS1, ECC1, P1, ENC1
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      
      if (verbose) print *, 'Create star with mass', mass, 'and age tag', age
     
      ML = log10(mass)
      SM1 = mass
      IM1 = (ML - MLO)/DM + 1.501D0
      CALL LOAD_STAR_MODEL(16,IM1, H, SM1,DTY1,AGE1,PER1,BMS1,ECC1,P1,ENC1,KH1,JN1,JF1)
      KH = KH1

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
      if (verbose) print *, 'Setting local pointer to star', num_stars
      star => star_list(num_stars)

      ! Allocate memory to store the star
      star%number_of_variables = actual_number_of_variables
      star%number_of_meshpoints = KH
      allocate(star%h(star%number_of_variables, star%number_of_meshpoints))
      allocate(star%dh(star%number_of_variables, star%number_of_meshpoints))
      allocate(star%hpr(star%number_of_variables, star%number_of_meshpoints))

      ! Make sure we eliminate the redundancy in storing and loading the model
      ! that is in TWIN (main/beginn).
      ! Essentially, we push the initialised data onto the stack
      JNN = 0
      call swap_out_star(num_stars)
      
      ! Estimate nuclear timescale for this star
      TNUC = 1.0d10 * mass**(-2.8)
      !star%dt = CSY * DTY1
      !star%dt = 3.0e2*CSY
      star%dt = CT3*1.0d-4 * TNUC*CSY
      if (mass > 1.0 .and. mass < 1.2) star%dt = 1.0d-2*star%dt
      !star%dt = 100
      if (verbose) print *, 'Setting initial timestep for star',num_stars,'to',star%dt,'s'
      
      ! Initialise some more variables
      star%zams_mass = mass
      star%age = age
      if (P1>0) star%p = P1
      star%jhold = 0
      star%jmod = 1
      star%jf = JF1
      star%er(:) = 0.0D0
      
      star%startup_iter = KR1
      star%normal_iter = KR2

      ! Binary orbital parameters
      if (PER1>0) star%per = PER1
      if (BMS1>0) star%bms = BMS1
      if (ECC1>0) star%ecc = ECC1
      
c Determine whether I and phi are computed or not, for OUTPUT
      star%jf = 0
      DO I = 11, 50
         IF ( ID(I)==12 .OR. ID(I)==14 ) star%jf = star%jf + 1
      ENDDO
      IF ( star%jf==1 ) star%jf = 0

      ! The star will still need some initialisation (initial timestep, say)
      star%virgin = .true.
      
      ! Save all settings
      call swap_in_star(num_stars)
      call swap_out_star(num_stars)

      load_zams_star = num_stars
      end function

! evolve_star:
!  evolve a star for one timestep; essentially star12 without the loop
!  TODO: ZAHB construction; maybe also for central carbon burning or white dwarfs?
      function evolve(star_id)
      USE MESH
      USE FUDGE_CONTROL
      USE CONSTANTS
      USE SETTINGS
      IMPLICIT REAL*8 (A-H, L-Z)    ! Ugly hack, needs to be changed to none!
      integer :: evolve
      integer, intent(in) :: star_id
      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, ID(130), IE(130) 
      COMMON /QUERY / ML(3), UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /TVBLES/ DT, ZQ(81), PR(81), PPR(81), IHOLD, JM2(2)
      
      ! We need some kind of initialisation for the star, cf the calls
      ! to printb in star12
      ! Actually, we need to pull some data from printb anyway...
      
      call swap_in_star(star_id)
      
      JO = 0
! Solve for structure, mesh, and major composition variables
      JOC = 1

      DTY = DT/CSY
      JTER = 0
      if (verbose) print *, 'Taking timestep ',DTY,'yr for star', current_star
    
      
      CALL SMART_SOLVER ( ITER, ID, KT5, JO )
      if (verbose) print *, 'Did', jter, 'out of', iter, 'iterations'

      ! Converged if JO == 0
      IF ( JO /= 0 ) THEN
         if (verbose) print *, 'Failed to converge on timestep'
!        If no convergence, restart from 2 steps back, DT decreased substantially
!        Abort if timestep below limit
         do while (JO /= 0)
            CALL BACKUP ( DTY, JO )
            IF ( JO.EQ.2 ) THEN
               evolve = JO
               return
            endif
            CALL NEXTDT ( DTY, JO, IT )

            if (verbose) print *, 'Timestep reduced to', DTY,'yr'

            IF ( JO.EQ.3 ) THEN
               evolve = JO
               return
            endif
            JNN = JNN + 1
            JO = 0
            CALL SMART_SOLVER ( ITER, ID, KT5, JO )
         end do
      END IF
      
      evolve = JO
      IF (JO == 0) THEN
         if (verbose) print *, 'Converged on timestep'
         
         call timestep_heuristic ( JO, 1 )
         
         CALL UPDATE ( DTY )
         CALL NEXTDT ( DTY, JO, 22 )
         IF ( JNN > switch_iterations ) ITER = KR2
         JNN = JNN + 1
         
         call swap_out_star(star_id)
      END IF
      end function

! dump_twin_model:
!  write the state of star id to a named file, for later retrieval.
!  The file will be in the format of a TWIN input file
      subroutine dump_twin_model(id, filename)
      implicit none
      integer, intent(in) :: id
      character*500, intent(in) :: filename
      
      call select_star(id)

      open (unit=40, file=filename, action='write')
      call output(200, 40, 0, 0)
      close (40);
      
      call flush_star()
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get_*
!  Selection of getter functions for luminosity, mass, radius, time, next dt etc.
!  These return retarded information from the *previous* timestep (cf Ross)
!  This makes sure that we will never need to backtrack beyond this point in
!  case the code needs to go back to a previous model.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! get_luminosity:
!  Returns the star's luminosity, in solar units
      function get_luminosity(id, value)
      use constants
      implicit none
      integer :: get_luminosity
      double precision :: value
      integer, intent(in) :: id
      
      if (id<1 .or. id>num_stars) then
         get_luminosity = -1
         value = 0.0
         return
      end if
      
      get_luminosity = 0
      value = star_list(id)%hpr(INV_VAR_PERM(8), 1)/CLSN
      end function
      
! get_mass:
!  Returns the star's mass, in solar units
      function get_mass(id, mass)
      use constants
      implicit none
      integer :: get_mass
      integer, intent(in) :: id
      double precision, intent(out) :: mass
      if (id<1 .or. id>num_stars) then
         mass = 0.0
         get_mass = -1
         return
      end if
      
      mass = star_list(id)%hpr(INV_VAR_PERM(4), 1)/CMSN
      get_mass = 0
      end function
      
! get_radius:
!  Returns the star's radius, in solar units
      function get_radius(id, radius)
      use constants
      implicit none
      
      integer :: get_radius
      double precision, intent(out) :: radius
      integer, intent(in) :: id
      
      if (id<1 .or. id>num_stars) then
         radius = 0.0
         get_radius = -1
         return
      end if
      
      radius = exp(star_list(id)%hpr(INV_VAR_PERM(7), 1))/CRSN
      get_radius = 0
      end function
      
! get_temperature:
!  Returns the star's effective temperature, in Kelvin
      function get_temperature(id, value)
      use constants
      implicit none
      integer :: get_temperature
      double precision :: value
      integer, intent(in) :: id
      
      if (id<1 .or. id>num_stars) then
         get_temperature = -1
         value = 0
         return
      end if
      
      get_temperature = 0
      value = exp(star_list(id)%hpr(INV_VAR_PERM(2), 1))
      end function
      
! get_age:
!  Returns the star's age, in years
      function get_age(id, value)
      use constants
      implicit none
      integer :: get_age
      double precision :: value
      integer, intent(in) :: id
      
      if (id<1 .or. id>num_stars) then
         get_age = -1
         value = 0.0
         return
      end if
      
      value = star_list(id)%age
      get_age = 0
      end function
      


! Return the stellar type of the specified star
      function get_stellar_type(id, value)
      implicit none
      integer :: get_stellar_type
      integer :: value
      integer, intent(in) :: id
      
      if (id<1 .or. id>num_stars) then
         get_stellar_type = -1
         value = -1
         return
      end if
      value = star_list(id)%stellar_type
      get_stellar_type = 0
      end function
    
! Return the metallicity parameter
      function get_metallicity(value)
      implicit none
      integer :: get_metallicity
      double precision :: value
      character*503 :: metallicity_str_with_zero = '0.0'
      
      metallicity_str_with_zero(3:) = metallicity_str
      READ (metallicity_str_with_zero, '(F10.8)') value
      get_metallicity = 0
      end function
      
      function set_metallicity(value)
      implicit none
      integer :: set_metallicity, i
      integer :: must_remove_zero = 1
      double precision :: value
      character*503 :: metallicity_str_with_zero = '0.0'
      
      WRITE(metallicity_str_with_zero,'(F10.8)') value
      must_remove_zero  = 1
      DO i = 10, 3, -1
        if ((metallicity_str_with_zero(i:i) .EQ. '0') .AND. (must_remove_zero .EQ. 1)) then
            metallicity_str_with_zero(i:i) = ' '
         else
            must_remove_zero = 0
        end if
      END DO
      metallicity_str = metallicity_str_with_zero(3:)
      
      set_metallicity = 0
      end function

    ! Return the maximum_number_of_stars parameter
      function get_maximum_number_of_stars(value)
      implicit none
      integer :: get_maximum_number_of_stars
      integer :: value
      value = maximum_number_of_stars
      get_maximum_number_of_stars = 0
      end function
      
      function set_maximum_number_of_stars(value)
      implicit none
      integer :: set_maximum_number_of_stars
      integer :: value
      if (max_stars > 0) then
        set_maximum_number_of_stars = -1
        return
      end if
      maximum_number_of_stars = value
      
      set_maximum_number_of_stars = 0
      end function
      
      function get_number_of_particles(value)
      implicit none
      integer :: get_number_of_particles
      integer :: value
      value = num_stars
      get_number_of_particles = 0
      end function


      function initialize_stars()
      implicit none
      integer :: initialize_stars
      initialize_stars = 0
      end function
      
      function new_particle(star_id, mass)
      implicit none
      integer :: new_particle, star_id
      double precision :: mass, age
      age = 0.0
      star_id = load_zams_star(mass, age)
      if (star_id .EQ. 0) then
        new_particle = -1
      else
        new_particle = 0
      end if
      end function
      
      
      function initialize_code()
      implicit none
      integer :: initialize_code
      initialize_code = initialise_twin(ev_path, maximum_number_of_stars, metallicity_str)
      end function
      
      function delete_star(star_id)
      implicit none
      integer :: delete_star, star_id
      delete_star =  -1
      end function

      subroutine timestep_heuristic ( JO, JSTAR )
      implicit none
      integer, intent(out) :: JO
      integer, intent(in) :: JSTAR
      integer :: JCM
      integer, external :: find_stellar_type
      CALL COMPUTE_OUTPUT_QUANTITIES ( JSTAR )
      CALL UPDATE_TIMESTEP_PARAMETERS ( JSTAR )
      CALL CHECK_STOP_CONDITIONS ( JSTAR, JO, JCM )
      star_list(current_star)%stellar_type = find_stellar_type()
      end subroutine



      subroutine set_star_iter_parameters( ID, KT1, KT2, JNN )
      implicit none
      integer, intent(in) :: id, KT1, KT2, JNN
      
      if (id<1 .or. id>num_stars) then
         return
      end if
      
      if (KT1>0) star_list(id)%startup_iter = KT1
      if (KT2>0) star_list(id)%normal_iter = KT2
      star_list(id)%jnn = jnn
      end subroutine



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
      USE CONSTANTS
      USE SETTINGS
      IMPLICIT REAL*8 (A-H, L-Z)    ! Ugly hack, needs to be changed to none!
      integer, intent(in) :: star_id
      type(twin_star_t), pointer :: star
      COMMON H(NVAR,NM), DH(NVAR,NM), EP(3), KH, KTW, KW(260)
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /TVBLES/ DT, ZQ(9), AGE, BM(71), PR(81), PPR(81), JHOLD, JM2, JM1
      COMMON /STORE / HPR(NVAR,NM), HT(4,NM), MS(9999), ST(50026)
      COMMON /SOLV  / C(NM+1,NEQ,NVAR+1), S(NEQ,121), ER(NVAR), KEQ, KJ2, KJ5, 
     :         KJ6, KJ10, KJ12, KI4, KEE, IA(1), KE1, KE2, KE3, KBC, 
     :         KEV, KFN, KL, JH1, JH2, JH3, KD(120), KQ, KVB, KVC
      integer :: n, nv, nvar
      
      ! Nothing to do if the star has already been swapped in
      if (current_star == star_id .and. .not. star_list(star_id)%virgin) return
      
      if (verbose) print *, 'Swapping in star', star_id
   
      current_star = star_id
      star => star_list(star_id)
      
      ! Store structure and structure changes
      KH = star%number_of_meshpoints
      nvar = star%number_of_variables
      do n=1, nvar
         nv = VAR_PERM(n)
         H(nv, 1:KH) = star%h(n,1 :KH)
         DH(nv, 1:KH) = star%dh(n, 1:KH)
      end do

      ! COMMON block QUERY
      ML = star%ml
      QL = star%ql
      XL = star%xl
      UC(:) = star%uc(:)
      JMOD = star%jmod
      JNN = star%jnn
      JTER = star%jter
      JOC = star%joc
      JKH = star%jkh
      
      ! COMMON block TVBLES
      DT = star%dt
      ZQ(:) = star%zq(:)
      AGE = star%age
      BM(:) = star%bm(:)
      PR(:) = star%pr(:)
      PPR(:) = star%ppr(:)
      JHOLD = star%jhold
      JM2 = star%jm2
      JM1 = star%jm1
      
      ! COMMON block SOLV. We need to retrieve the typical size of the different
      ! variables. It's not a good idea to leak this information to other stars
      ER(:) = star%er(:)
      
      ! COMMON block STORE. We don't need the information about the primary's
      ! mass loss rate
      do n=1, nvar
         nv = VAR_PERM(n)
         HPR(nv, 1:KH) = star%hpr(n,1 :KH)
      end do

      KR1 = star%startup_iter
      KR2 = star%normal_iter
      ITER = star%normal_iter
      IF (JNN < switch_iterations) ITER = star%startup_iter
      
      ! Does this star need some initialisation?
      if (star%virgin) then
         if (verbose) print *, 'Performing one time initialisation for star ', star_id
         star%virgin = .false.
         
         ! Some tasks from beginn
         JO = 0
         IT = 22
         JNN = 0
         DTY = DT/CSY
         CALL NEXTDT ( DTY, JO, IT )
         ! Do we care about REMESH at this stage? Maybe not!
         BMS = CMSN*star%bms
         ECC = star%ecc
         PER = star%per
         TM = star%zams_mass * CMSN
         TOA = CG1*TM*(BMS - TM)*(CG2*PER/BMS)**C3RD*DSQRT(1.0D0 - ECC*ECC)
         CALL REMESH ( KH2, JCH, BMS, TM, star%p, ECC, TOA, 1, star%jf )

         HPR(:,1:KH) = H(:,1:KH) 
         DH(:, 1:KH) = 0.0D0
         JHOLD = 2
         PR(1:81) = ZQ(1:81)  ! May generate index out-of-bound warning
         PPR(1:81) = PR(1:81)
         JM1 = JMOD

         ! some tasks from printb, notably the timestep control
         ! For now, we just call timestep_heuristic, which means we get some unwanted output
         call timestep_heuristic ( JO, 1 )
         CALL NEXTDT ( DTY, JO, IT )
         JNN = 1

      end if
      end subroutine
      
! swap_out_star:
!  Swap out a particular star for later evolution with the TWIN code
!  Assumes all stars are evolved with the same set of equations and physical
!  parameters. This is fine as long as these are not changed on the fly, as
!  for the post He-flash generation code. This needs to be disabled from 
!  init.run, or we need to swap in and out the data from init.dat and init.run
      subroutine swap_out_star(star_id)
      IMPLICIT REAL*8 (A-H, L-Z)    ! Ugly hack, needs to be changed to none!
      integer, intent(in) :: star_id
      type(twin_star_t), pointer :: star
      COMMON H(NVAR,NM), DH(NVAR,NM), EP(3), KH, KTW, KW(260)
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /TVBLES/ DT, ZQ(9), AGE, BM(71), PR(81), PPR(81), JHOLD, JM2, JM1
      COMMON /STORE / HPR(NVAR,NM), HT(4,NM), MS(9999), ST(50026)
      COMMON /SOLV  / C(NM+1,NEQ,NVAR+1), S(NEQ,121), ER(NVAR), KEQ, KJ2, KJ5, 
     :         KJ6, KJ10, KJ12, KI4, KEE, IA(1), KE1, KE2, KE3, KBC, 
     :         KEV, KFN, KL, JH1, JH2, JH3, KD(120), KQ, KVB, KVC
      integer :: n, nv, nvar
 
      if (star_id < 1) return
      if (verbose) print *, 'Swapping out star', star_id
      star => star_list(star_id)
      
      ! Retrieve structure and structure changes
      nvar = star%number_of_variables
      do n=1, nvar
         nv = VAR_PERM(n)
         star%h(n, 1:KH) = H(nv,1 :KH)
         star%dh(n, 1:KH) = DH(nv, 1:KH)
      end do
      star%number_of_meshpoints = KH

      ! COMMON block QUERY
      star%ml = ML
      star%ql = QL
      star%xl = XL
      star%uc(:) = UC(:)
      star%jmod = JMOD
      star%jnn = JNN
      star%jter = JTER
      star%joc = JOC
      star%jkh = JKH
      
      ! COMMON block TVBLES
      star%dt  = DT
      star%zq(:) = ZQ(:)
      star%age = AGE
      star%bm(:) = BM(:)
      star%pr(:) = PR(:)
      star%ppr(:) = PPR(:)
      star%jhold = JHOLD
      star%jm2 = JM2
      star%jm1 = JM1
      
      ! COMMON block SOLV. We need to retrieve the typical size of the different
      ! variables. It's not a good idea to leak this information to other stars
      star%er(:) = ER(:)
      
      ! COMMON block STORE. We don't need the information about the primary's
      ! mass loss rate
      do n=1, nvar
         nv = VAR_PERM(n)
         star%hpr(n, 1:KH) = HPR(nv,1 :KH)
      end do

      end subroutine

! select_star:
!  flushes the data from the current star back to the cache and selects the
!  requested star as the new current star.
      subroutine select_star(id)
      implicit none
      integer :: id
      call swap_out_star(current_star)
      call swap_in_star(id)
      end subroutine

! flush_star:
!  flushes the data from the currently selected star back to the cache
      subroutine flush_star()
      implicit none
      call swap_out_star(current_star)
      end subroutine

! Helper functions, from aamain.f

! Initialise constants, read in tables of physical data that are unchanged
! during the course of one run      
      SUBROUTINE SETSUP
      USE SETTINGS
      USE CONSTANTS
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON /UBVDAT/ TGR(34), GGR(13), TAB(34,13,5)

! Initialise some constants that cannot be computed at compile time on some
! compilers
      CALL INITIALISE_CONSTANTS

! Read nuclear reaction rates and neutrino loss parameters
! Used to be read from the same file as the opacity data, which is stupid
!  since this is independent of metalicity
      CALL LOAD_REACTION_AND_NEUTRINO_RATES(42)
 
! Load opacity tables - needs to be done after reading init.dat, because we need
!  to know what tables to read.
!      CALL LOAD_OPACITY(20)

c Read Bol. Corr, U-B, B-V table. 
  991 FORMAT (3(10F10.5,/), 4F10.5,/, 10F10.5,/, 3F10.5,/, 442(5f8.3,/))
      READ (21,991) TGR, GGR, (((TAB(I,J,K), K=1,5), J=1,13), I=1,34)
      CLOSE (21)

c Read nuclear reaction (QRT) and neutrino (QNT) Q values, in MeV; constants 
c for electron screening (CZA, CZB, CZC, CZD, VZ); atomic parameters (CBN, KZN),
c with masses (CAN) consistent with Q-values; ionization potentials (CHI) and 
c statistical weights (COM); molecular hydrogen parameters (CH2)
      CALL LOAD_ATOMIC_DATA(26)

      RETURN
      END SUBROUTINE
      
      
      
      SUBROUTINE READ_INIT_RUN(IR)
      USE MESH
      USE MESH_ENC
      USE SETTINGS
      USE FILE_EXISTS_MODULE
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IR
      INTEGER :: S_KPT
      INTEGER :: IOERROR
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EP(3)
      DOUBLE PRECISION :: SM0, DTY0, AGE0, PER0, BMS0, ECC0, P0, ENC0
      DOUBLE PRECISION :: SM, DTY, AGE, PER, BMS, ECC, P, ENC
      DOUBLE PRECISION :: ML, QL, XL, UC(21)
      INTEGER :: JMOD, JB, JNN, JTER, JOC, JKH, KH, KTW, KW(260)
      INTEGER :: ISB, IP1, IM1, IP2, IM2, KPT, KP, KML, KQL, KXL, KR, JMX
      DOUBLE PRECISION :: ML1, DML, QL1, DQL, XL1, DXL, ROT, EX
      CHARACTER*500 :: STARTFILE, NAME, ZAMSFILE
      COMMON H, DH, EP, KH, KTW, KW
      COMMON /TN1/ SM0, DTY0, AGE0, PER0, BMS0, ECC0, P0, ENC0
      COMMON /T0 / SM, DTY, AGE, PER, BMS, ECC, P, ENC
      COMMON /QUERY / ML, QL, XL, UC, JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /SAVEINIT/ S_KPT
      COMMON /CINIT_RUN/ ISB, IP1, IM1, IP2, IM2, KPT, KP, ML1, DML, KML,
     & QL1, DQL, KQL, XL1, DXL, KXL, ROT, KR, EX, JMX, STARTFILE, ZAMSFILE
! EG: Namelist for namelist I/O of init.run (fort.23)      
      NAMELIST /INIT_RUN/ ISB, KTW, IP1, IM1, IP2, IM2, KPT, KP, ML1, DML, KML,
     & QL1, DQL, KQL, XL1, DXL, KXL, ROT, KR, EX, SM, DTY, AGE, PER, BMS, ECC, 
     & P, ENC, JMX, UC, START_MODEL_ONLY, STARTFILE, ZAMSFILE,
     & START_WITH_RIGID_ROTATION
c First try to use namelist I/O; if that fails use `normal' old-fashioned IO
      STARTFILE = ''
      READ (IR, NML=INIT_RUN, IOSTAT=IOERROR)
      REWIND (IR)
      IF (IOERROR /= 0) THEN
         READ (IR, *) ISB, KTW, IP1, IM1, IP2, IM2, KPT, KP, ML1, DML, KML, 
     :      QL1, DQL, KQL, XL1, DXL, KXL, ROT, KR, EX, SM, DTY, AGE, PER, BMS,
     :      ECC, P, ENC, JMX, UC
         REWIND (IR)
      END IF
c If a starting model was named in the input file, use that instead.
c No more fort.nnn symbolic links, yay!      
      IF ( LEN(TRIM(STARTFILE))>0 ) THEN
         IF (.NOT. FILE_EXISTS(STARTFILE)) THEN
            INQUIRE(UNIT=IR, NAME=NAME)
            WRITE (0, *) 'Start file "', TRIM(STARTFILE), '" not found. Check setting in ', TRIM(NAME), '.'
            STOP
         ELSE
            IP1 = 62
            OPEN(UNIT = IP1, ACTION="READ", FILE=STARTFILE)
         END IF
      END IF
! EG: store some information we may later need so we don't have to reread the file
      S_KPT = KPT
! EG: Do some magic for ISB == -1, -2 meaning we want to do a mutation rather than
!  a regular evolution run. A mutation run is a run where we start with a
!  normal star and end up with something weird that we want to evolve next.
!  Useful for getting merger remnants, for instance, to evolve
      IF (ISB == -1 .OR. ISB == -2) THEN
         ! Make sure all parts of the code that need to know know we're mutating
         MUTATE = .TRUE.
         ! Keep the type of mutation we want to do         
         MUTATION_MODE = ISB
         ! We really do want the code to calculate things for one star
         ISB = 1
         ! Read in target model global data
         TPROF_IN = IP1
         READ (TPROF_IN, *) SM0, DTY0, AGE0, PER0, BMS0, ECC0, P0, ENC0
         REWIND (TPROF_IN)
         ! To start our mutation run, we first need a suitable model from the
         !  main sequence library, so we have to setup some variables to the
         !  mass that we want.
         SM = SM0
         ML1 = LOG10(SM0)
         IP1 = 16
      END IF
      END SUBROUTINE



! Load model number IM from the file IP
      SUBROUTINE LOAD_STAR_MODEL(IP,IM, H,SM,DTY,AGE,PER,BMS,ECC,P,ENC,KH,JN,JF)
      USE MESH
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IP, IM
      INTEGER, INTENT(OUT) :: KH,JN,JF
      DOUBLE PRECISION, INTENT(OUT) :: SM,DTY,AGE,PER,BMS,ECC,P,ENC
      DOUBLE PRECISION, INTENT(OUT) :: H(NVAR,NM)
      INTEGER :: I, J, K, JM

      DO I = 1, IM
c Skip all models upto IM
         READ (IP, *) SM, DTY, AGE, PER, BMS, ECC, P, ENC, KH, K, JM, K, JN, JF
         DO K = 1, KH
            READ (IP, *) (H(J, K), J = 1, JN)
            END DO
      END DO
      REWIND (IP) 
      END SUBROUTINE

      end module
