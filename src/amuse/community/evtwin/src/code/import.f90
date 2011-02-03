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
function import_stellar_merger(nmesh, numvar, model, age_tag)
   use real_kind
   use twin_library_v2
   use control
   use mesh_enc
   use constants
   use settings
   use init_dat
   use atomic_data
   use test_variables
   use extra_elements
   use interpolate
   use binary_history, only: hpr
   
   implicit none
   integer :: import_stellar_merger
   integer, intent(in) :: nmesh, numvar
   real(double), intent(in) :: model(numvar, nmesh)
   real(double), intent(in) :: age_tag
   real(double) :: xa(9), na(9)

   type(init_dat_settings) :: initdat
   integer :: kh2, kr1, kr2, ksv, kt5, jch
   integer :: n, i, star_id, iter, status
   real(double) :: mass, entropy_max, target_diffsqr, vma, tkh
   character(len=500) :: outputfilename, basename

   basename = "star"

   call push_init_dat(initdat, kh2, kr1, kr2, ksv, kt5, jch)
   kr2 = get_number_of_iterations()
   kx = 0
   ky = 0
   kz = 0
   !> \todo FIXME: there are more mixing and mass loss options that can be
   !! set in the init.dat file, these need to be stored/restored as
   !! well!
   !<
   cth = 0.0
   crd = 0.0
   cmr = 0.0
   cmj = 0.0
   cmv = 0.0
   cmk = 0.0
   cmnl = 0.0
   mass = model(1, 1)

   ! Invert equation of state
   print *, 'inverting eos...'
   do n=1, nmesh
      ! Convert mass fractions back baryon number fractions
      na(1:9) = model(5:13, n)
      do i=1, 9
         xa(i) = na(i) * cbn(i)/can(i)
      end do
      vma = sum(xa(1:9))
      xa(1:9) = xa(1:9) / vma
      th(4, n) = model(1, n) * cmsn    ! Mass coordinate
      th(5, n) = xa(1)                 ! Hydrogen abundance
      th(9, n) = xa(2)
      th(10, n) = xa(3)
      th(16, n) = xa(4)
      th(3, n) = xa(5)
      th(11, n) = xa(6)
      th(NMg24, n) = xa(7)
      th(NSi28, n) = xa(8)
      th(NFe56, n) = xa(9)
      call prtoft (model(4, n), model(3, n), th(1, n), th(2, n), xa)
   end do
   ip_mesh_size = nmesh

   ! Construct interpolation tables
   print *, 'constructing interpolation tables'
   do n = 1, 16
      call iptable_init (nmesh, th(4,:),th(n,:),thb(n,:),thc(n,:),thd(n,:))
   end do

   print *, 'loading zams star of mass', mass
   call flush_star
   star_id = load_zams_star(mass, -1.0d0, 0.0d0)
   call select_star(star_id)

   ! Stage 1: match composition
   adj_comp = .true.
   impose_composition_factor = 0.0d0
   do iter=1, 100
      print *, 'composition adjustment factor =', impose_composition_factor
      status = twin_evolve()
      call flush_star()
      if (status /= 0) then
         print *, '*** failed to converge on timestep', iter, 'with code', status
         stop
      end if

      if (impose_composition_factor>=1.0d0) exit
      if (impose_composition_factor<=1.0d-2) then
         impose_composition_factor = min(1.5d0*impose_composition_factor, 1.0d0);
      else
         impose_composition_factor = min(1.2d0*impose_composition_factor, 1.0d0);
      end if
      impose_composition_factor = max(impose_composition_factor, 1.0d-4);

      call flush(6)
   end do

   ! Store output
   call flush_star()
   outputfilename = trim(basename)//'.comp_mod'
   print *, 'writing output to ', trim(outputfilename)
   call dump_twin_model(star_id, outputfilename);

   ! Stage 2: adjust entropy profile, keep composition fixed
   usemenc = .true.
   impose_entropy_factor = 0.0d0
   entropy_max = 1.0d2
   curr_diffsqr = 1.0d3
   best_diffsqr = 1.0d3
   target_diffsqr = eps
   target_diffsqr = 1.0e-4
   call set_number_of_iterations(20)
   do iter=1, 100
      age = 0.0
      print *, 'entropy adjustment factor =', impose_entropy_factor
      ! Construct the next stellar model in the pseudo-evolution
      ! sequence. Make sure the timestep is always close to the
      ! thermal time scale for the best accuracy. Make sure the
      ! solution is stable at this timestep by iterating on it while
      ! only changing the timestep.
      status = twin_evolve()
      tkh = 1.0d22*cg*h(4, 1)*2/(exp(h(7, 1))*h(8, 1)*csy)
      do while (status == 0 .and. dt < 10.0*csy)
         !print *, 'Grow timestep'
         dt = 1.01*dt
         status = twin_evolve()
      end do
      do while (status == 0 .and. dt > tkh * csy)
         !print *, 'Shrink timestep'
         dt = dt*0.8
         status = twin_evolve()
      end do
      !print *, DT/CSY, TKH
      dt = tkh * csy
      call flush_star()
      if (status /= 0) then
         print *, '*** failed to converge on timestep', iter, 'with code', status
         if (impose_entropy_factor >= 1.0d0) exit;
         stop
      end if

      ! Check convergence, adjust entropy matching factor
      call check_conv_to_target_structure()
      if ( curr_diffsqr < best_diffsqr ) then
         best_diffsqr = curr_diffsqr
         best_mod = iter
         mutant_h(1:24, 1:kh) = h(1:24, 1:kh)
      end if
      !WRITE (6, '(1P, 3D16.9, I6)'), CURR_MAXDIFFSQR, CURR_DIFFSQR, BEST_DIFFSQR, BEST_MOD
      print *, 'converged to',best_diffsqr,'at',best_mod,'now', iter


      if ( ( impose_entropy_factor>=entropy_max .and. best_diffsqr>curr_diffsqr ) .or. best_diffsqr<target_diffsqr ) exit
      if (impose_entropy_factor < 1.0) then
         impose_entropy_factor = min(1.5d0*impose_entropy_factor, entropy_max);
      else if (impose_entropy_factor < 1.0d2) then
         impose_entropy_factor = min(1.25d0*impose_entropy_factor, entropy_max);
      else if (impose_entropy_factor < 1.0d3) then
         impose_entropy_factor = min(1.05d0*impose_entropy_factor, entropy_max);
      else
         impose_entropy_factor = min(1.01d0*impose_entropy_factor, entropy_max);
      end if
      impose_entropy_factor = max(impose_entropy_factor, 1.0d-8);

      call flush(6)
   end do

   call pop_init_dat(initdat, kh2, kr1, kr2, ksv, kt5, jch)
   h(:, 1:kh) = mutant_h(:, 1:kh)
   hpr(:, 1:kh) = mutant_h(:, 1:kh)
   adj_comp = .false.
   usemenc = .false.
   impose_entropy_factor = 0.0d0
   impose_composition_factor = 0.0d0
   ! Set timestep
   !DT = 1.0D3 * CSY
   dt = tkh * csy
   ! Set age: make sure this is not reset when we go back one model
   age = age_tag
   prev(10) = age
   pprev(10) = age
   call flush_star()
   call set_star_iter_parameters( star_id, 10, 20, 0 )

   ! Store output
   outputfilename = trim(basename)//'.pmutate'
   print *, 'writing output to ', trim(outputfilename)
   call dump_twin_model(star_id, outputfilename);

   import_stellar_merger = star_id
end function import_stellar_merger

