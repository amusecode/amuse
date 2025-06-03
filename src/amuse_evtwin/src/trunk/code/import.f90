module import

   use real_kind
   
   ! variables for storing the user supplied model
   implicit none
   private
   
   real(double), pointer :: new_model(:,:)
   logical :: new_model_defined = .false.
   integer :: n_shells_new_model
   
   public :: import_stellar_merger
   public :: new_stellar_model
   public :: finalize_stellar_model
   

contains

! This module is designed to work along the MUSE library.
! It takes an array with pressure, density, mass coordinate and
! composition for each mesh point (to be extended with angular momentum)
! and constructs a new stellar model, in the vein of mkmergermod.
! This facility is forked off to its own module because (for the moment)
! it makes use of numerical recipes code that shouldn't taint the main
! library file.
! Order of variables stored at each meshpoint:
!  Mass coordinate [Msun], Radius [Rsun], log density [cgs],
!  log pressure [cgs], XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE
! Mespoints should be passed in from *surface* to *centre*
! The optional argument verb indicates whether verbose output to the console is wanted or not
! The optional argument force indicates how hard we try to match the entropy profile.
! Higher values give better agreement, but may be harder to evolve succesfully
! Likewise, the argument epsilon gives the least-square difference required for convergence (default: 1e-4)
   integer function import_stellar_merger(star_id, nmesh, numvar, model, age_tag, verb, force, epsilon)
      use real_kind
      use twinlib
      use control
      use mesh_enc
      use constants
      use settings
      use init_dat
      use atomic_data
      use test_variables
      use interpolate
      use binary_history, only: hpr
      use indices
      implicit none
      integer, intent(inout) :: star_id
      integer, intent(in) :: nmesh, numvar
      real(double), intent(in) :: model(numvar, nmesh)
      real(double), intent(in), optional :: age_tag
      real(double), intent(in), optional :: force
      real(double), intent(in), optional :: epsilon
      logical, intent(in), optional :: verb
      real(double) :: xa(9), na(9)
   
      integer :: n, i, ik, status
      real(double) :: mass, entropy_acc, vma, xc
      real(double) :: fac_min_backup, eps_backup, crd_backup, cth_backup, wanted_eps_backup, cdc_backup(1:5)
      integer :: kr2_backup
      logical :: verbose
   
      verbose = .false.
      if (present(verb)) verbose = verb
   
      entropy_acc = 1.0d-4
      if (present(epsilon)) entropy_acc = epsilon
   
      mass = model(1, 1)
   
      ! Invert equation of state
      if (verbose) print *, 'inverting eos...'
      do n=1, nmesh
         ! Convert mass fractions back to baryon number fractions
         na(1:9) = model(5:13, n)
         do i=1, 9
            xa(i) = na(i) * cbn(i)/can(i)
         end do
         vma = sum(xa(1:9))
         xa(1:9) = xa(1:9) / vma
         th(VAR_MASS, n) = model(1, n) * cmsn    ! Mass coordinate
         th(VAR_H1, n) = xa(1)                 ! Hydrogen abundance
         th(VAR_HE4, n) = xa(2)
         th(VAR_C12, n) = xa(3)
         th(VAR_N14, n) = xa(4)
         th(VAR_O16, n) = xa(5)
         th(VAR_NE20, n) = xa(6)
         th(VAR_MG24, n) = xa(7)
         th(VAR_SI28, n) = xa(8)
         th(VAR_FE56, n) = xa(9)
         call prtoft (model(4, n), model(3, n), th(VAR_LNF, n), th(VAR_LNT, n), xa)
      end do
      ip_mesh_size = nmesh
      xc = th(VAR_H1, nmesh)
   
      ! Construct interpolation tables
      if (verbose) print *, 'constructing interpolation tables'
      do n = 1, 16
         call iptable_init (nmesh, th(VAR_MASS,:),th(n,:),thb(n,:),thc(n,:),thd(n,:))
      end do
   
      ! Look for a local minumum in ln f
      ! This will probably correspond to an inversion in the entropy and/or
      ! density profile, indicating a loss of resolution in the (SPH) output. We
      ! discard the entropy profile from this point onward, since we cannot
      ! trust it.
      maximum_match_mass = TH(VAR_MASS, 1)
      do ik = ip_mesh_size-1, 1, -1
         if (TH(1, ik) > TH(1, ik+1) .and. TH(VAR_MASS, ik) > 0.5*TH(VAR_MASS,1)) then
            if (verbose) print *, "Found inversion at m = ", TH(VAR_MASS,ik+1) / CMSN
            maximum_match_mass = TH(VAR_MASS, ik+1)
            maximum_match_mass = TH(VAR_MASS, ik+1) - 0.5
            exit
         end if
      end do
      if (verbose) print *, 'Maximum match mass = ', maximum_match_mass / CMSN
   
      if (verbose) print *, 'loading zams star of mass', mass
      status = new_zams_star(star_id, mass, 0.0d0)
   
      ! Turn off all mass loss during model construction
      status = set_wind_multiplier(star_id, 0.0d0)
   
      ! Stage 1: evolve until correct core properties
      do while (H(VAR_H1, kh) > xc)
         status = evolve_one_timestep(star_id)
         if (status /= 0)  then
            import_stellar_merger = status
            return
         end if
      end do
   
      fac_min_backup = fac_min
      eps_backup = eps
      crd_backup = crd
      cth_backup = cth
      wanted_eps_backup = wanted_eps
      kr2_backup = kr2
      cdc_backup(1:5) = cdc(1:5)
   
      ! Stage 2: tweak the composition and entropy profiles
      adj_mea = .true.
      adj_comp = .true.
      usemenc = .false.
      curr_diffsqr = 1.0d3;
      best_diffsqr = curr_diffsqr
      mutant_h(1:INDEX_SECONDARY_START, 1:kh) = h(1:INDEX_SECONDARY_START, 1:kh)
      impose_entropy_factor = 0.0
      impose_composition_factor = 1.0d-4
      fac_min = 1.0d0
      eps = 1.0d-6
      crd = 0.0d0
      cth = 0.0d0
      wanted_eps = 1.0d0
      entropy_force = 20.d0
      if (present(force)) entropy_force = force
      cdc(1) = 0.01
      cdc(2) = 0.25
      cdc(3) = 1.00
      cdc(4) = 4.00
      cdc(5) = 1.00
      kr2 = 100
      call set_timestep(star_id, 1.0d3)
      status = 0
      i = 0
      do while (status == 0 .and. i < 2000)
         if (timestep_of(star_id) > 1.0d6) call set_timestep(star_id, 1.0d6)
         status = evolve_one_timestep(star_id, verbose)
         if ( best_diffsqr<entropy_acc .and. get_composition_mean_square() < 1.0d-4) status = 53
         i = i+1
      end do
   
      import_stellar_merger = -2
      if (status == 53) import_stellar_merger = status
   
      ! Cleanup
      h(1:INDEX_SECONDARY_START, 1:kh) = mutant_h(1:INDEX_SECONDARY_START, 1:kh)
      call set_timestep(star_id, 1.0d3)
   
      adj_mea = .false.
      adj_comp = .false.
      usemenc = .false.
      fac_min = fac_min_backup
      eps = eps_backup
      crd = crd_backup
      cth = cth_backup
      wanted_eps = wanted_eps_backup
      kr2 = kr2_backup
      cdc(1:5) = cdc_backup(1:5) 
   
      ! Turn on mass loss for subsequent evolution
      status = set_wind_multiplier(star_id, 1.0d0)
   
      if (present(age_tag)) age = age_tag
   
   end function import_stellar_merger
   
   
   ! Each shell has: Mass coordinate [Msun], Radius [Rsun], density [cgs],
   !  pressure [cgs], XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE
   ! Meshpoints should be passed in from *surface* to *centre*
   integer function new_stellar_model(mass, radius, rho, pressure, &
         XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, n)
      use real_kind
      use twinlib
      use mesh, only: max_nm
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: mass(n), radius(n), rho(n), pressure(n), &
         XH(n), XHE(n), XC(n), XN(n), XO(n), XNE(n), XMG(n), XSI(n), XFE(n)
      double precision :: logrho(n), logpressure(n)
!~      real(double), pointer :: new_model(:,:)
      
      if (n < 3) then
         new_stellar_model = -31   ! Too few shells in new empty model
         return
      else if (n > max_nm) then
         new_stellar_model = -32   ! Too many shells in new empty model
         return
      endif
      
      allocate(new_model(13, n))
      n_shells_new_model = n
      new_stellar_model = 0    ! A new empty model was defined
      
      new_model(1, 1:n) = mass(1:n)
      new_model(2, 1:n) = radius(1:n)
      logrho(1:n) = log(rho(1:n))
      logpressure(1:n) = log(pressure(1:n))
      new_model(3, 1:n) = logrho(1:n)
      new_model(4, 1:n) = logpressure(1:n)
      new_model(5, 1:n) = XH
      new_model(6, 1:n) = XHE
      new_model(7, 1:n) = XC
      new_model(8, 1:n) = XN
      new_model(9, 1:n) = XO
      new_model(10, 1:n) = XNE
      new_model(11, 1:n) = XMG
      new_model(12, 1:n) = XSI
      new_model(13, 1:n) = XFE
      new_model_defined = .true.
!~      new_stellar_model = import_stellar_merger(star_id, n, 13, new_model, 0.0d0, &
!~         amuse_verbose, amuse_entropy_force, amuse_entropy_accuracy)
!~      deallocate(new_model)
   end function


   integer function finalize_stellar_model(star_id, age_tag)
      use twinlib, only: amuse_verbose, amuse_entropy_force, amuse_entropy_accuracy
      implicit none
      integer, intent(out) :: star_id
      double precision, intent(in) :: age_tag
     
      if (.not. new_model_defined) then
         finalize_stellar_model = -35
         return
      endif
      finalize_stellar_model = import_stellar_merger(star_id, n_shells_new_model, 13, new_model, age_tag, &
         amuse_verbose, amuse_entropy_force, amuse_entropy_accuracy)
      call flush()
     
      if (star_id .eq. 0) then
         finalize_stellar_model = -36
      end if
      deallocate(new_model)
      new_model_defined = .false.
   end function

end module
