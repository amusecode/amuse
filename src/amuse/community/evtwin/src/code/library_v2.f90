! Library interface to TWIN: useful for putting TWIN into Muse

module twin_library_v2
   use real_kind
   use mesh
   use extrapolate_dh

   ! Datastructure for storing all the information for one star
   ! These need to be swapped into and out off the TWIN global variables
   ! We only care about single stars and TWIN mode binaries
   type, private :: twin_star_t
      sequence
      ! Flag to indicate wether this star still needs some initialisation
      logical :: virgin
         
      ! Flag that tells whether this star still exists (i.e. it has not been removed)
      ! Flag that tells whether this star still exists (i.e. it has not been removed)
      logical :: star_exists

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
      
      type(extrapolate_dh_data) :: stored_extrapolate_dh_data
      
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
   integer, private :: sx_updated_for_star = 0 ! Stellar structure data sx is up-to-date for this star
   integer, private :: highest_star_index = 0  ! Highest index assigned to a star so far
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
      character(len=500), intent(in) :: new_init_dat_name;

      if (.not. file_exists(new_init_dat_name) ) then
         if (verbose) print *, "warning: file ",trim(new_init_dat_name)," for ", trim(init_dat_name), " does not exist!"
         set_init_dat_name = -1
         return
      end if
      init_dat_name = trim(new_init_dat_name)
      set_init_dat_name = 0
   end function set_init_dat_name


   ! set_init_run_name:
   !  Change the name of the input file used for the numerical settings and
   !  physics options (default: ./init.dat)
   function set_init_run_name(new_init_run_name)
      use file_exists_module
      
      implicit none
      integer :: set_init_run_name
      character(len=500), intent(in) :: new_init_run_name;

      if (.not. file_exists(new_init_run_name) ) then
         if (verbose) print *, "warning: file ",trim(new_init_run_name)," for ", trim(init_run_name), " does not exist!"
         set_init_run_name = -1
         return
      end if
      init_run_name = new_init_run_name;
      set_init_run_name = 0
   end function set_init_run_name

   function set_ev_path(new_ev_path)
      use file_exists_module
      
      implicit none
      integer :: set_ev_path
      character(len=*), intent(in) :: new_ev_path;

      if (.not. file_exists(new_ev_path) ) then
         if (.not. file_exists(trim(new_ev_path)//"/lt2ubv.dat") ) then
            if (verbose) print *, "Warning: file ",trim(new_ev_path)," for ", trim(ev_path), " does not exist!"
            set_ev_path = -1
            return
         end if
      end if
      ev_path = new_ev_path;
      set_ev_path = 0
   end function

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


    function open_input_file(inputfilename, inputunit, input_required)
        use file_exists_module
        implicit none
        character(len=*) :: inputfilename
        character(len=500) :: fname
        integer :: inputunit
        integer :: input_required
        integer :: open_input_file
        
        open_input_file = 0
        
        write (fname, '("fort.",i2)') inputunit
        if (file_exists(inputfilename) .and. .not. file_exists(fname)) then
            if (verbose) print *, "opening ",trim(inputfilename)," for ", trim(fname)
            open(unit = inputunit, action="read", file=inputfilename)
        end if
            ! Check if files exist and act appropriately
            ! If the file is required (INPUT_REQUIRED>0), abort on error
            ! If the file is optional (INPUT_REQUIRED==0), give a warning
            ! If the file is probably unneeded (INPUT_REQUIRED<0), do nothing
        if (.not. (file_exists(inputfilename)) .or. file_exists(fname)) then
            if (input_required > 0) then
               write (0, *) 'required input file ', trim(inputfilename), ' (', trim(fname), ') not found'
               open_input_file = -1
               return
            else if (input_required == 0) then
               write (0, *) 'warning: input file ', trim(inputfilename), ' (', trim(fname), ') not found'
            end if
        end if

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
      use real_kind
      use settings
      use extra_elements
      use constants
      use current_model_properties
      use control
      
      implicit none
      integer :: initialise_twin
      character(len=*), intent(in) ::  evpath, zstr
      integer, intent(in) :: nstars
      integer :: ii
      integer :: ke1,ke2,ke3,kbc,kev,kfn,kl,kp_var(40),kp_eqn(40),kp_bc(40)
      

     
      initialise_twin = open_input_file(trim(evpath)//"/zahb"//trim(zstr)//".mod", 12, 0)
      if(initialise_twin.LT.0) return

      initialise_twin = open_input_file(trim(evpath)//"/zahb"//".dat", 24, 0)
      if(initialise_twin.LT.0) return

      initialise_twin = open_input_file(trim(evpath)//"/zams/zams"//trim(zstr)//".mod", 16, 1)
      if(initialise_twin.LT.0) return

      initialise_twin = open_input_file(trim(evpath)//"/zams/zams"//trim(zstr)//".out", 18, 1)
      if(initialise_twin.LT.0) return

      initialise_twin = open_input_file(trim(evpath)//"/zams/zams"//trim(zstr)//".mas", 19, 1)
      if(initialise_twin.LT.0) return

      initialise_twin = open_input_file(trim(ev_path)//"/mutate.dat", 63, -1)
      if(initialise_twin.LT.0) return

      
      initialise_twin = open_input_file(trim(evpath)//"/metals/z"//trim(zstr)//"/phys.z"//trim(zstr), 20, 0)
      if(initialise_twin.LT.0) return

      initialise_twin = open_input_file(trim(evpath)//"/lt2ubv.dat", 21, 1)
      if(initialise_twin.LT.0) return

      initialise_twin = open_input_file(trim(evpath)//"/nucdata.dat", 26, 1)
      if(initialise_twin.LT.0) return

      initialise_twin = open_input_file(trim(evpath)//"/COtables/COtables_z"//trim(zstr), 41, 0)
      if(initialise_twin.LT.0) return

      initialise_twin = open_input_file(trim(evpath)//"/physinfo.dat", 42, 1)
      if(initialise_twin.LT.0) return
   
      max_stars = nstars
      if (verbose) then
         print *, 'twin initialisation.'
         print *, 'arguments: evpath =', evpath
         print *, '           nstars =', nstars
         print *, '           zstr =', zstr
      end if

      ! We want to override the default init.run and init.dat settings if
      ! init.run and init.dat are present in the current directory.
     
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
      
      if (verbose) print *, 'read settings'
      if (verbose) print *, 'using', kh2, 'meshpoints per star'
      
      if (verbose) print *, 'files initialised, initialising equation of state'
      call setsup

      ! Read layout of the ZAMS library
      read (19, *) mlo, dm, mhi, kdm
      rewind (19)

      
      ! Create short (`pruned') summary of ZAMS models from long file
      ! Are we using this? Not reall, I think...
      !      CALL PRUNER ( 18, 17, ISB )
      !      CLOSE (17)
      !      CLOSE (18)


      ! Autodetect if we should solve for Mg24 or not by checking if the
      ! corresponding equation is in the list of equations
      use_mg24_eqn = .false.
      do ii = 51, 100
         if (id(ii) == emg24 .or. id(ii) == esumx) then
            use_mg24_eqn = .true.
            exit
         end if
         if (id(ii) == 0) exit   ! Break loop if end of list found
      end do

       ! Detect whether rotation is treated as solid body rotation or
       ! whether we consider differential rotation. The switch is on whether
       ! or not the rotational period (var. 13) is listed as an eigenvalue.
       rigid_rotation = .true.
       do ii = 11, 10 + id(1)+id(2)
          if (id(ii) == 13) then  ! Rotational pertiod not an EV
             rigid_rotation = .false.
             exit
          end if
       end do

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
      ! Allocate data for nucleosynthesis calculations
      call allocate_nucleosynthesis_data(kh2)

      
      if(use_quadratic_predictions) then
        call initlse_parabola_storage_space(kh2, ke1+ke2+kev)
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
   function load_zams_star(mass, spin, age)
      use real_kind
      use mesh
      use constants
      use settings
      use current_model_properties
      use control
      use starting_values0
      
      implicit none
      type(twin_star_t), pointer :: star
      integer :: load_zams_star, new_id
      real(double), intent(in) :: mass, spin, age

      integer :: im1,kh1,kp1,jmod1,jb1,jn1,jf1
      real(double) :: tnuc
      real(double) :: hnuc(50, nm)

      !Common blocks:
      real(double) :: sm1, dty1, age1, per1, bms1, ecc1, p1, enc1
      real(double) :: tn, perc

      
      if (verbose) print *, 'create star with mass', mass, 'and age tag', age
      sm1 = 0
      ml = log10(mass)
      !sm1 = mass 
      !sm1 is not an input, so commented this line out
      im1 = (ml - mlo)/dm + 1.501d0
        
      ! input for load star model is 
      ! 1. the file input number (16 stands for zams star)
      ! 2. the log10 mass of the star
      ! mlo, dm are read from evpath/zams/zams/metallicity.mas
      
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
      
      ! Find an unused index for the new star
      if (num_stars < highest_star_index) then
        new_id = 1
        do while (star_list(new_id)%star_exists)
          new_id = new_id + 1
        end do
        if (new_id > max_stars) then
          print *, 'Error - If you get this message, blame Nathan.'
          load_zams_star = 0
          return
        end if          
      else
        highest_star_index = highest_star_index + 1
        new_id = highest_star_index
      end if
      num_stars = num_stars + 1
      if (verbose) print *, 'setting local pointer to star', new_id
      star => star_list(new_id)

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
      call swap_out_star(new_id)


      xl = 0.5
      tn = 1.0d10 * h(4, 1)/h(8, 1) * clsn/cmsn
      perc = sqrt( exp(3.0 * h(7, 1))/(8.157_dbl*h(4, 1)) * cmsn/crsn**3 )


      ! Estimate nuclear timescale for this star
      tnuc = 1.0d10 * mass**(-2.8)
      star%dt = (1.0d-4*tn) * csy
      star%per = perc*10.0_dbl**xl
      if (mass > 1.0 .and. mass < 1.2) star%dt = 1.0d-2*star%dt
      if (verbose) print *, 'setting initial timestep for star',new_id,'to',star%dt,'s'

      ! Initialise some more variables
      star%star_exists = .true.
      star%zams_mass = mass
      star%age = age
      if (p1>0) star%p = p1
! If the spin is supplied by the user, it will override the current value.
      if (spin>0.0) star%p = spin
      star%jhold = 0
      star%jmod = 0
      star%jf = jf1
      star%er(:) = 0.0d0
      star%enc = 0.0d0

      star%startup_iter = kr1
      star%normal_iter = kr2
      
      if (use_quadratic_predictions) then
        print *,'use_quadratic_predictions ',use_quadratic_predictions
        call push_extrapolate_data(star%stored_extrapolate_dh_data)
      end if

      ! Binary orbital parameters
      if (per1>0) star%per = per1
      if (bms1>0) star%bms = bms1
      if (ecc1>0) star%ecc = ecc1
      
      if (t0dty >= 0.0_dbl) star%dt= t0dty * csy
      if (t0per >= 0.0_dbl) star%per = t0per
      if (t0bms >= 0.0_dbl) star%bms = t0bms
      if (t0ecc >= 0.0_dbl) star%ecc = t0ecc
      if (t0p >= 0.0_dbl)   star%p   = t0p
      if (t0enc >= 0.0_dbl) star%enc = t0enc
      
      ! taken form begin.f90
      ! this code:if ( jip.eq.13 .or. jip.eq.14 ) dty = ct3*dty
      star%dt = ct3 * star%dt
      

      ! The star will still need some initialisation (initial timestep, say)
      star%virgin = .true.
      if (verbose) print *, 'setting initial timestep for star',new_id,'to',star%dt,'s'

      ! Save all settings
      call swap_in_star(new_id)
      call swap_out_star(new_id)
      if (verbose) print *, 'setting initial timestep for star',new_id,'to',star%dt,'s'

      load_zams_star = new_id
   end function load_zams_star



   ! evolve_star:
   !  evolve a star for one timestep; essentially star12 without the loop
   !  TODO: ZAHB construction; maybe also for central carbon burning or white dwarfs?
   function evolve(star_id)
      use constants
      implicit none
      integer :: evolve
      integer, intent(in) :: star_id
      
      ! Check whether the star exists, and that it hasn't been removed.
      if (star_id<1 .or. star_id>highest_star_index .or. .not. star_list(star_id)%star_exists) then
         if (verbose) print *, 'Error: no star to evolve. Star_id=', star_id, ', star_exists=', star_list(star_id)%star_exists
         evolve = -1
         return
      end if
      
      call swap_in_star(star_id)
      evolve = twin_evolve()
      if (evolve == 0) call swap_out_star(star_id)
   end function evolve
   
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
   function get_luminosity(id, value)
      use real_kind
      use constants
      
      implicit none
      real(double) :: value
      integer :: get_luminosity
      integer, intent(in) :: id

      if (id<1 .or. id>highest_star_index .or. .not. star_list(id)%star_exists) then
         get_luminosity = -1
         value = 0.0
         return
      end if
      get_luminosity = 0
      value = star_list(id)%hpr(inv_var_perm(8), 1)/clsn
   end function get_luminosity
   
   function get_wind_mass_loss_rate(id,value)
      use real_kind
      use constants
      
      implicit none
      integer :: get_wind_mass_loss_rate
      real(double) , intent(out) :: value
      integer, intent(in) :: id
      
      if (id<1 .or. id>highest_star_index .or. .not. star_list(id)%star_exists) then
         value = 0.0
         get_wind_mass_loss_rate = -1
         return
      end if
      
      value = star_list(id)%zet(1)*csy/cmsn
      
      get_wind_mass_loss_rate = 0
   end function
   
  
   ! msun / yr
   function get_mass_transfer_rate(id, value)
      use real_kind
      use constants
      
      implicit none
      integer :: get_mass_transfer_rate
      real(double) , intent(out) :: value
      integer, intent(in) :: id
      
      if (id<1 .or. id>highest_star_index .or. .not. star_list(id)%star_exists) then
         value = 0.0
         get_mass_transfer_rate = -1
         return
      end if
      
      value = star_list(id)%xit(1)*csy/cmsn
      
      get_mass_transfer_rate = 0
   end function
   
! get_mass:
!  Returns the star's mass, in solar units
   function get_mass(id, value)
      use real_kind
      use constants
      
      implicit none
      real(double) :: value
      integer :: get_mass
      integer, intent(in) :: id
      if (id<1 .or. id>highest_star_index .or. .not. star_list(id)%star_exists) then
         value = 0.0
         get_mass = -1
         return
      end if
      value = star_list(id)%hpr(inv_var_perm(4), 1)/cmsn
      get_mass = 0
   end function get_mass

! get_radius:
!  Returns the star's radius, in solar units
   function get_radius(id, value)
      use real_kind
      use constants
      
      implicit none
      real(double) :: value
      integer :: get_radius
      integer, intent(in) :: id

      if (id<1 .or. id>highest_star_index .or. .not. star_list(id)%star_exists) then
         value = 0.0
         get_radius = -1
         return
      end if
      value = exp(star_list(id)%hpr(inv_var_perm(7), 1))/crsn
      get_radius = 0
   end function get_radius

! get_temperature:
!  Returns the star's effective temperature, in Kelvin
   function get_temperature(id, value)
      use real_kind
      use constants
      
      implicit none
      real(double) :: value
      integer :: get_temperature
      integer, intent(in) :: id

      if (id<1 .or. id>highest_star_index .or. .not. star_list(id)%star_exists) then
         get_temperature = -1
         value = 0
         return
      end if
      get_temperature = 0
      value = exp(star_list(id)%hpr(inv_var_perm(2), 1))
   end function get_temperature

! get_age:
!  Returns the star's age, in years
   function get_age(id, value)
      use real_kind
      use constants
      
      implicit none
      real(double) :: value
      integer :: get_age
      integer, intent(in) :: id

      if (id<1 .or. id>highest_star_index .or. .not. star_list(id)%star_exists) then
         get_age = -1
         value = 0.0
         return
      end if
      value = star_list(id)%age
      get_age = 0
   end function get_age

! Return the stellar type of the specified star
   function get_stellar_type(id, value)
      use real_kind
      
      implicit none
      integer :: get_stellar_type, value
      integer, intent(in) :: id

      if (id<1 .or. id>highest_star_index .or. .not. star_list(id)%star_exists) then
         get_stellar_type = -1
         value = -1
         return
      end if
      value = star_list(id)%stellar_type
      get_stellar_type = 0
   end function get_stellar_type


! get_metallicity:
! Return the metallicity parameter
   function get_metallicity(value)
      implicit none
      integer :: get_metallicity
      real(double) :: value
      character*503 :: metallicity_str_with_zero = '0.0'
      metallicity_str_with_zero(3:) = metallicity_str
      read (metallicity_str_with_zero, '(F10.8)') value
      get_metallicity = 0
   end function
   
   
! set_metallicity:
! Set the metallicity parameter
   function set_metallicity(value)
      implicit none
      integer :: set_metallicity, i
      integer :: must_remove_zero = 1
      real(double) :: value
      character*503 :: metallicity_str_with_zero = '0.0'
      write (metallicity_str_with_zero,'(F10.8)') value
      must_remove_zero  = 1
      do i = 10, 3, -1
        if ((metallicity_str_with_zero(i:i) .eq. '0') .and. (must_remove_zero .eq. 1)) then
            metallicity_str_with_zero(i:i) = ' '
         else
            must_remove_zero = 0
        end if
      end do
      metallicity_str = metallicity_str_with_zero(3:)
      set_metallicity = 0
   end function

! get_maximum_number_of_stars:
! Return the maximum_number_of_stars parameter
   function get_maximum_number_of_stars(value)
      implicit none
      integer :: get_maximum_number_of_stars
      integer :: value
      value = maximum_number_of_stars
      get_maximum_number_of_stars = 0
   end function

! set_maximum_number_of_stars:
! Set the maximum_number_of_stars parameter
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

! get_number_of_particles:
! Return the number of stars created
   function get_number_of_particles(value)
      implicit none
      integer :: get_number_of_particles
      integer :: value
      value = num_stars
      get_number_of_particles = 0
   end function

! get_max_age_stop_condition:
! Return the maximum age stop condition
   function get_max_age_stop_condition(value)
      use real_kind
      use current_model_properties
      implicit none
      integer :: get_max_age_stop_condition
      real(double) :: value
      value = uc(2)
      get_max_age_stop_condition = 0
   end function

! set_max_age_stop_condition:
! Set the maximum age stop condition
   function set_max_age_stop_condition(value)
      use real_kind
      use current_model_properties
      implicit none
      integer :: set_max_age_stop_condition
      real(double) :: value
      uc(2) = value
      set_max_age_stop_condition = 0
   end function

! get_min_timestep_stop_condition:
! Return the minimum timestep stop condition
   function get_min_timestep_stop_condition(value)
      use real_kind
      use current_model_properties
      use constants
      implicit none
      integer :: get_min_timestep_stop_condition
      real(double) :: value
      value = uc(12)/csy ! Convert from seconds to yr
      get_min_timestep_stop_condition = 0
   end function

! set_min_timestep_stop_condition:
! Set the minimum timestep stop condition
   function set_min_timestep_stop_condition(value)
      use real_kind
      use current_model_properties
      use constants
      implicit none
      integer :: set_min_timestep_stop_condition
      real(double) :: value
      uc(12) = value*csy ! Convert from yr to seconds
      set_min_timestep_stop_condition = 0
   end function

! get_time_step:
! Return the current time step to be taken for the evolution of this star.
   function get_time_step(id, value)
      use real_kind
      use constants
      implicit none
      integer :: get_time_step
      integer, intent(in) :: id
      real(double) :: value
      if (id<1 .or. id>highest_star_index .or. .not. star_list(id)%star_exists) then
         get_time_step = -1
         value = -1.0
         return
      end if
      value = star_list(id)%dt/csy ! Convert from seconds to yr
      get_time_step = 0
   end function

! get_number_of_ionization_elements:
! Retrieve the current number of elements used for ionization in the EoS
   function get_number_of_ionization_elements(value)
      use settings
      implicit none
      integer :: get_number_of_ionization_elements, value
      value = kion
      get_number_of_ionization_elements = 0
   end function

! set_number_of_ionization_elements:
! Set the current number of elements used for ionization in the EoS
   function set_number_of_ionization_elements(value)
      use settings
      implicit none
      integer :: set_number_of_ionization_elements, value
      kion = value
      set_number_of_ionization_elements = 0
   end function

! get_convective_overshoot_parameter:
! Retrieve the current value of the convective overshoot parameter
   function get_convective_overshoot_parameter(value)
      use real_kind
      use settings
      implicit none
      integer :: get_convective_overshoot_parameter
      real(double) :: value
      value = cos
      get_convective_overshoot_parameter = 0
   end function

! set_convective_overshoot_parameter:
! Set the current value of the convective overshoot parameter
   function set_convective_overshoot_parameter(value)
      use real_kind
      use settings
      implicit none
      integer :: set_convective_overshoot_parameter
      real(double) :: value
      cos = value
      set_convective_overshoot_parameter = 0
   end function

! get_mixing_length_ratio:
! Retrieve the current value of the mixing length ratio
   function get_mixing_length_ratio(value)
      use real_kind
      use settings
      implicit none
      integer :: get_mixing_length_ratio
      real(double) :: value
      value = calp
      get_mixing_length_ratio = 0
   end function

! set_mixing_length_ratio:
! Set the current value of the mixing length ratio
   function set_mixing_length_ratio(value)
      use real_kind
      use settings
      implicit none
      integer :: set_mixing_length_ratio
      real(double) :: value
      calp = value
      set_mixing_length_ratio = 0
   end function

! get_semi_convection_efficiency:
! Retrieve the current value of the efficiency of semi-convection
   function get_semi_convection_efficiency(value)
      use real_kind
      use settings
      implicit none
      integer :: get_semi_convection_efficiency
      real(double) :: value
      value = csmc
      get_semi_convection_efficiency = 0
   end function

! set_semi_convection_efficiency:
! Set the current value of the efficiency of semi-convection
   function set_semi_convection_efficiency(value)
      use real_kind
      use settings
      implicit none
      integer :: set_semi_convection_efficiency
      real(double) :: value
      csmc = value
      set_semi_convection_efficiency = 0
   end function

! get_thermohaline_mixing_parameter:
! Retrieve the current value of the thermohaline mixing parameter
   function get_thermohaline_mixing_parameter(value)
      use real_kind
      use settings
      implicit none
      integer :: get_thermohaline_mixing_parameter
      real(double) :: value
      value = cth
      get_thermohaline_mixing_parameter = 0
   end function

! set_thermohaline_mixing_parameter:
! Set the current value of the thermohaline mixing parameter
   function set_thermohaline_mixing_parameter(value)
      use real_kind
      use settings
      implicit none
      integer :: set_thermohaline_mixing_parameter
      real(double) :: value
      cth = value
      set_thermohaline_mixing_parameter = 0
   end function

! get_AGB_wind_setting:
! Get the current setting for mass-loss (AGB)
   function get_AGB_wind_setting(value)
      use real_kind
      use massloss
      implicit none
      integer :: get_AGB_wind_setting
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
   function set_AGB_wind_setting(value)
      use massloss
      implicit none
      integer :: set_AGB_wind_setting
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
   function get_RGB_wind_setting(value)
      use real_kind
      use massloss
      implicit none
      integer :: get_RGB_wind_setting
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
   function set_RGB_wind_setting(value)
      use real_kind
      use massloss
      implicit none
      integer :: set_RGB_wind_setting
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

! Initialize the stars (does nothing)
   function initialize_stars()
      implicit none
      integer :: initialize_stars
      initialize_stars = 0
   end function

! Create a new particle (spin period not specified - will be set by code)
   function new_particle(star_id, mass)
      implicit none
      integer :: new_particle, star_id
      double precision, intent(in) :: mass
      double precision :: spin, age
      age = 0.0
      spin = -1.0
      star_id = load_zams_star(mass, spin, age)
      if (star_id .eq. 0) then
        new_particle = -1
      else
        new_particle = 0
      end if
   end function

! Create a new particle with specified spin period
   function new_spinning_particle(star_id, mass, spin)
      use real_kind
      implicit none
      integer :: new_spinning_particle, star_id
      real(double) :: mass, spin, age
      age = 0.0
      star_id = load_zams_star(mass, spin, age)
      if (star_id .eq. 0) then
        new_spinning_particle = -1
      else
        new_spinning_particle = 0
      end if
   end function

! get_spin:
! Returns the star's spin period, in days
   function get_spin(id, value)
      implicit none
      integer :: get_spin
      real(double), intent(out) :: value
      integer, intent(in) :: id      
      if (id<1 .or. id>highest_star_index .or. .not. star_list(id)%star_exists) then
         get_spin = -1
         value = 0.0
         return
      end if
      value = star_list(id)%p
      get_spin = 0
   end function


! Initialize the code
   function initialize_code()
      use real_kind
      use settings
      use file_exists_module
      use init_run
      use init_dat
      use control
      use constants
      use extra_elements
      use current_model_properties
      
      implicit none
      integer :: initialize_code
      logical :: status, override_run,  override_dat
      
      initialize_code = 0
      
      override_run = .false.
      override_dat = .false.
      
      initialize_code = open_input_file(trim(init_dat_name), 22, -1)
      if(initialize_code.LT.0) return

      initialize_code = open_input_file(trim(init_run_name), 23, -1)
      if(initialize_code.LT.0) return

      initialize_code = open_input_file(trim(ev_path)//"/init/init.dat", 122, 1)
      if(initialize_code.LT.0) return

      initialize_code = open_input_file(trim(ev_path)//"/init/init.run", 123, 1)
      if(initialize_code.LT.0) return
             

      if ( file_exists(trim(init_dat_name)) ) override_dat = .true.
      if ( file_exists(trim(init_run_name)) ) override_run = .true.

      ! We need to have a fort.11 - no, we don't, since we're not using star12!
      !WRITE (11, *) 0
      !CLOSE (11)

      ! We need a scratch file for output we don't want or need
      if (.not. file_exists('fort.25'))&
           open (unit=25, action='write', status = 'scratch')

      ! Other redundant (for our purposes) output we do not want to save but cannot
      ! easily avoid: solver output to fort.1 and opacity table output to fort.10
      open (unit=10, action='write', status = 'scratch')

      if (verbose) print *, 'initialised equation of state, reading run settings'
      call read_init_run(123)    ! Load defaults
      if (override_run) call read_init_run(23)
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
         initialize_code = -1
         return
      end if
      ! this might help...
      
      ! Convert some things to `cgs' units: 10**11 cm, 10**33 gm, 10**33 erg/s
      cmi = cmi/csy
      cmj = cmj*cmsn/csy
      cms = cms/csy
      cmt = cmt*1.0d-11
         
      
   end function
    
   function commit_parameters()
      implicit none
      integer :: commit_parameters
            
      commit_parameters = initialise_twin(ev_path, maximum_number_of_stars, metallicity_str)
    
    end function
   
   function recommit_parameters()
      implicit none
      integer :: recommit_parameters
      recommit_parameters = -2
   end function  
       
   function cleanup_code()
      implicit none
      integer :: cleanup_code
      cleanup_code = -2
   end function  
       
   function delete_star(star_id)
      implicit none
      integer :: delete_star, star_id
      type(twin_star_t), pointer :: star
      ! Check whether the star exists, and that it hasn't been removed.
      if (star_id<1 .or. star_id>highest_star_index) then
        delete_star = -1
        return
      end if
      star => star_list(star_id)
      if (.not. star%star_exists) then
        if (verbose) print *, 'Error: no star to delete (star was already removed).'
        delete_star = -1
        return
      end if
      ! Deallocate memory of the star
      deallocate(star%h)
      deallocate(star%dh)
      deallocate(star%hpr)
      ! Tell the bad news to the star 
      star%star_exists = .false.
      num_stars = num_stars - 1
      delete_star = 0
   end function


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

      if (id<1 .or. id>highest_star_index .or. .not. star_list(id)%star_exists) then
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
      use control
      
      implicit none
      integer, intent(in) :: star_id
      type(twin_star_t), pointer :: star
      integer :: n, nv, number_of_variables
      integer :: jo,it, i
      real(double) :: dty,bms,ecc,per,tm,toa
      
      sx_updated_for_star = 0
      
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

      if (use_quadratic_predictions) then
        call pop_extrapolate_data(star%stored_extrapolate_dh_data)
      end if
      
      ! COMMON block SOLV. We need to retrieve the typical size of the different
      ! variables. It's not a good idea to leak this information to other stars
      if (allocated(er)) er(1:nvar) = star%er(1:nvar)

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
         ! well the stand-alone code does...
         
         bms = cmsn*star%bms
         ecc = star%ecc
         per = star%per
         tm = star%zams_mass * cmsn
         toa = cg1*tm*(bms - tm)*(cg2*per/bms)**c3rd*dsqrt(1.0d0 - ecc*ecc)
         
         call remesh ( kh2, jch, bms, tm, star%p, ecc, toa, 1, star%jf )
        
       
         ! Determine whether I and phi are computed or not, for OUTPUT
         ! Note, do this afer remesh, jf is used in remesh and remesh
         ! expects that the jf comes from the load_star_model
         star%jf = 0
         do i = 11, 50
            if ( id(i)==12 .or. id(i)==14 ) star%jf = star%jf + 1
         end do
         if ( star%jf==1 ) star%jf = 0
         
         ! flush(6)
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
      use control
      
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


      if (use_quadratic_predictions) then
        call push_extrapolate_data(star%stored_extrapolate_dh_data)
      end if
      
      star%prev(:) = prev(:)
      star%pprev(:) = pprev(:)
      star%jhold = jhold
      star%jm2 = jm2
      star%jm1 = jm1

      ! COMMON block SOLV. We need to retrieve the typical size of the different
      ! variables. It's not a good idea to leak this information to other stars
      if (allocated(er)) then
         star%er(1:nvar) = er(1:nvar)
      else
         star%er(1:nvar) = 0.0d0
      end if

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


! Internal structure getters:

! Return the current number of zones/mesh-cells of the star
      integer function get_number_of_zones(AMUSE_id, AMUSE_value)
         implicit none
         integer, intent(in) :: AMUSE_id
         integer, intent(out) :: AMUSE_value
         if (AMUSE_id<1 .or. AMUSE_id>highest_star_index .or. .not. star_list(AMUSE_id)%star_exists) then
            AMUSE_value = 0
            get_number_of_zones = -21
            return
         end if
         AMUSE_value = star_list(AMUSE_id)% number_of_meshpoints
         get_number_of_zones = 0
      end function

! Update the stellar structure data sx for the specified star, if necessary
      subroutine update_quantities_if_needed(AMUSE_id)
         implicit none
         integer, intent(in) :: AMUSE_id
         if (sx_updated_for_star .ne. AMUSE_id) then
            call swap_in_star(AMUSE_id)
            call compute_output_quantities ( 1 )
            sx_updated_for_star = AMUSE_id
         end if
      end subroutine

! Return the temperature at the specified zone/mesh-cell of the star
      integer function get_temperature_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
         use structure_variables
         implicit none
         integer, intent(in) :: AMUSE_id, AMUSE_zone
         double precision, intent(out) :: AMUSE_value
         if (AMUSE_id<1 .or. AMUSE_id>highest_star_index .or. .not. star_list(AMUSE_id)%star_exists) then
            AMUSE_value = -1.0
            get_temperature_at_zone = -21
            return
         end if
         if (AMUSE_zone >= star_list(AMUSE_id)% number_of_meshpoints .or. AMUSE_zone < 0) then
            AMUSE_value = -1.0
            get_temperature_at_zone = -22
            return
         end if
         call update_quantities_if_needed(AMUSE_id)
         AMUSE_value = sx(4, AMUSE_zone+2)
!         AMUSE_value = exp(star_list(AMUSE_id)% hpr(inv_var_perm(2), star_list(AMUSE_id)% number_of_meshpoints-AMUSE_zone))
         get_temperature_at_zone = 0
      end function

! Return the density at the specified zone/mesh-cell of the star
      integer function get_density_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
         use structure_variables
         implicit none
         integer, intent(in) :: AMUSE_id, AMUSE_zone
         double precision, intent(out) :: AMUSE_value
         if (AMUSE_id<1 .or. AMUSE_id>highest_star_index .or. .not. star_list(AMUSE_id)%star_exists) then
            AMUSE_value = -1.0
            get_density_at_zone = -21
            return
         end if
         if (AMUSE_zone >= star_list(AMUSE_id)% number_of_meshpoints .or. AMUSE_zone < 0) then
            AMUSE_value = -1.0
            get_density_at_zone = -22
            return
         end if
         call update_quantities_if_needed(AMUSE_id)
         AMUSE_value = sx(3, AMUSE_zone+2)
         get_density_at_zone = 0
      end function

! Return the radius at the specified zone/mesh-cell of the star
      integer function get_radius_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
         use structure_variables
         use constants
         implicit none
         integer, intent(in) :: AMUSE_id, AMUSE_zone
         double precision, intent(out) :: AMUSE_value
         if (AMUSE_id<1 .or. AMUSE_id>highest_star_index .or. .not. star_list(AMUSE_id)%star_exists) then
            AMUSE_value = -1.0
            get_radius_at_zone = -21
            return
         end if
         if (AMUSE_zone >= star_list(AMUSE_id)% number_of_meshpoints .or. AMUSE_zone < 0) then
            AMUSE_value = -1.0
            get_radius_at_zone = -22
            return
         end if
         call update_quantities_if_needed(AMUSE_id)
         AMUSE_value = sx(17, AMUSE_zone+2) * CRSN * 1.0D11
         get_radius_at_zone = 0
      end function

! Return the mean molecular weight per particle (ions + free electrons) at the specified zone/mesh-cell of the star
      integer function get_mu_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
         use structure_variables
         implicit none
         integer, intent(in) :: AMUSE_id, AMUSE_zone
         double precision, intent(out) :: AMUSE_value
         if (AMUSE_id<1 .or. AMUSE_id>highest_star_index .or. .not. star_list(AMUSE_id)%star_exists) then
            AMUSE_value = -1.0
            get_mu_at_zone = -21
            return
         end if
         if (AMUSE_zone >= star_list(AMUSE_id)% number_of_meshpoints .or. AMUSE_zone < 0) then
            AMUSE_value = -1.0
            get_mu_at_zone = -22
            return
         end if
         call update_quantities_if_needed(AMUSE_id)
         AMUSE_value = sx(31, AMUSE_zone+2)
         get_mu_at_zone = 0
      end function

! Return the current number of chemical abundance variables per zone of the star
   integer function get_number_of_species(AMUSE_id, AMUSE_value)
      implicit none
      integer, intent(in) :: AMUSE_id
      integer, intent(out) :: AMUSE_value
      AMUSE_value = 9
      get_number_of_species = 0
   end function

! Return the name of chemical abundance variable 'AMUSE_species' of the star
   integer function get_name_of_species(AMUSE_id, AMUSE_species, AMUSE_value)
      implicit none
      integer, intent(in) :: AMUSE_id, AMUSE_species
      character (len=6), intent(out) :: AMUSE_value
      get_name_of_species = 0
      select case (AMUSE_species)
         case (1)
            AMUSE_value = 'h1'
         case (2)
            AMUSE_value = 'he4'
         case (3)
            AMUSE_value = 'c12'
         case (4)
            AMUSE_value = 'n14'
         case (5)
            AMUSE_value = 'o16'
         case (6)
            AMUSE_value = 'ne20'
         case (7)
            AMUSE_value = 'mg24'
         case (8)
            AMUSE_value = 'si28'
         case (9)
            AMUSE_value = 'fe56'
         case default
            AMUSE_value = 'error'
            get_name_of_species = -23
      end select
   end function

! Return the mass fraction of species 'AMUSE_species' at the specified 
! zone/mesh-cell of the star
   integer function get_mass_fraction_of_species_at_zone(AMUSE_id, &
         AMUSE_species, AMUSE_zone, AMUSE_value)
      use atomic_data
      use extra_elements
      use binary_history, only: hpr
      implicit none
      integer, intent(in) :: AMUSE_id, AMUSE_zone, AMUSE_species
      double precision, intent(out) :: AMUSE_value
      real(double) :: xa(9), na(9)
      real(double) :: avm
      integer :: zone_index, i
  
      AMUSE_value = -1.0
      if (AMUSE_id<1 .or. AMUSE_id>highest_star_index .or. .not. star_list(AMUSE_id)%star_exists) then
         get_mass_fraction_of_species_at_zone = -21
         return
      end if
      if (AMUSE_zone >= star_list(AMUSE_id)% number_of_meshpoints .or. AMUSE_zone < 0) then
         get_mass_fraction_of_species_at_zone = -22
         return
      end if
      if (AMUSE_species > 9 .or. AMUSE_species < 1) then
         get_mass_fraction_of_species_at_zone = -23
         return
      end if
      call update_quantities_if_needed(AMUSE_id)
      zone_index = star_list(AMUSE_id)% number_of_meshpoints - AMUSE_zone
      xa(1) = hpr(5, zone_index)
      xa(2) = hpr(9, zone_index)
      xa(3) = hpr(10, zone_index)
      xa(4) = hpr(16, zone_index)
      xa(5) = hpr(3, zone_index)
      xa(6) = hpr(11, zone_index)
      xa(8) = hpr(NSi28, zone_index)
      xa(9) = hpr(NFe56, zone_index)
      xa(7) = 1.0d0 - sum(xa(1:6)) - sum(xa(8:9))
      do i=1, 9
        na(i) = xa(i) * can(i)/cbn(i)
      end do
      avm = sum(na(1:9))
      AMUSE_value = na(AMUSE_species) / avm
      get_mass_fraction_of_species_at_zone = 0
   end function

   integer function set_mass_fraction_of_species_at_zone(AMUSE_id, &
         AMUSE_species, AMUSE_zone, AMUSE_value)
      implicit none
      integer, intent(in) :: AMUSE_id, AMUSE_zone, AMUSE_species
      double precision, intent(in) :: AMUSE_value
      set_mass_fraction_of_species_at_zone = -4
   end function
   
   integer function set_density_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
      implicit none
      integer, intent(in) :: AMUSE_id, AMUSE_zone
      double precision, intent(in) :: AMUSE_value
      set_density_at_zone = -4
   end function
   
   integer function set_temperature_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
      implicit none
      integer, intent(in) :: AMUSE_id, AMUSE_zone
      double precision, intent(in) :: AMUSE_value
      set_temperature_at_zone = -4
   end function
   
   integer function set_radius_at_zone(AMUSE_id, AMUSE_zone, AMUSE_value)
      implicit none
      integer, intent(in) :: AMUSE_id, AMUSE_zone
      double precision, intent(in) :: AMUSE_value
      set_radius_at_zone = -4
   end function
   
   integer function set_mass(AMUSE_id, AMUSE_value)
      implicit none
      integer, intent(in) :: AMUSE_id
      double precision, intent(in) :: AMUSE_value
      set_mass = -4
   end function


end module twin_library_v2


