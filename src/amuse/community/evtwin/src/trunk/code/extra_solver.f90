!     Support functions for the solver.
!     There are a few types of these:
!
!      * Functions to assist solve() in finding a solution:
!           SCALE_EQN()
!      * Functions to clean up a solution found by solve():
!           CHECKS()
!      * Functions that wrap around solve() and that can ease it through
!        situations where convergence is difficult:
!        SMART_SOLVER(), FUDGED_SOLVER(), RELAX_TIMESTEP()

! ------------------------------------------------------------------------------
!  SCALE_EQN
!   Calculate the expected scaling for equations and (surface) boundary
!   conditions, based on the current values in the H() array and the
!   current timestep.
!   Calls FUNCS1() to compute quantities that enter the equations
! ------------------------------------------------------------------------------
!  Input:
!   H(:)          -  Current values of independent variables (taken from
!                    the nameless COMMON block)
!   DT            -  Current timestep, in s (taken from module test_variables)
!  Output:
!   EQN_SCALE(:)  -  Expected scaling for equations
!   SBC_SCALE(:)  -  Expected scaling for surface boundary conditions
! ------------------------------------------------------------------------------
subroutine scale_eqn ( eqn_scale, sbc_scale )
   use real_kind
   use mesh
   use test_variables
   use structure_functions
   use indices
   
   implicit none
   real(double), intent(out) :: eqn_scale(1:NVAR), sbc_scale(1:NVAR)

   ! COMMON block INF: pass input to FUNCS1, get return values from FUNCS1
   real(double) :: var(nvar), dvar(nvar), fn1(nfunc)
   real(double) ::  vpk, vrk, vtk, lk
   real(double) :: lq, mt, vmk, X1T, X16T, X4T
   real(double) :: X12T, X20T, vik, phik
   real(double) :: X14T
   real(double) :: X24T
   real(double) :: adam, samt

   ! Local variables:
   integer :: ik, i
   integer, parameter :: Jstar = 1  ! Only do single star for now

   ! Initialise with default settings: all equations have a natural value of 1
   eqn_scale(:) = 0.0d0
   sbc_scale(:) = 1.0d0

   ! FIXME: fix equation scaling for new equation/function indices.
   do ik=kh, 1, -1
      !        Copy current value of variables
      forall (i=1:nvar)
         dvar(i) = dh(idx_for_star(i, Jstar), ik)
         var(i) = h(idx_for_star(i, Jstar), ik)
      end forall
      var = var + dvar

      !        Compute function values
      call funcs1 ( ik, -1, var(:), dvar(:), fn1(:) )
      vpk  = fn1(fn_idx_for_star(fn_vpk, Jstar))
      vrk  = fn1(fn_idx_for_star(fn_vrk, Jstar))
      vtk  = fn1(fn_idx_for_star(fn_vtk, Jstar))
      lk   = fn1(fn_idx_for_star(fn_lk, Jstar))
      lq   = fn1(fn_idx_for_star(fn_lq, Jstar))
      mt   = fn1(fn_idx_for_star(fn_mt, Jstar))
      vmk  = fn1(fn_idx_for_star(fn_vmk, Jstar))
      x1t  = fn1(fn_idx_for_star(fn_x1t, Jstar))
      x16t = fn1(fn_idx_for_star(fn_x16t, Jstar))
      x4t  = fn1(fn_idx_for_star(fn_x4t, Jstar))
      x12t = fn1(fn_idx_for_star(fn_x12t, Jstar))
      x20t = fn1(fn_idx_for_star(fn_x20t, Jstar))
      vik  = fn1(fn_idx_for_star(fn_vik, Jstar))
      phik = fn1(fn_idx_for_star(fn_phik, Jstar))
      x14t = fn1(fn_idx_for_star(fn_x14t, Jstar))
      x24t = fn1(fn_idx_for_star(FN_X24, Jstar))
      adam = fn1(fn_idx_for_star(FN_ADAM, Jstar))
      samt = fn1(fn_idx_for_star(FN_SAMT, Jstar))

      !        Composition equations
      eqn_scale(EQN_H1) = max(eqn_scale(EQN_H1), abs(X1T))
      eqn_scale(EQN_HE4) = max(eqn_scale(EQN_HE4), abs(X4T))
      eqn_scale(EQN_C12) = max(eqn_scale(EQN_C12), abs(X12T))
      eqn_scale(EQN_N14) = max(eqn_scale(EQN_N14), abs(X14T))
      eqn_scale(EQN_O16) = max(eqn_scale(EQN_O16), abs(X16T))
      eqn_scale(EQN_NE20) = max(eqn_scale(EQN_NE20), abs(X20T))
      eqn_scale(EQN_MG24) = max(eqn_scale(EQN_MG24), abs(X24T))

      !        Structure equations
      eqn_scale(EQN_PRES) = max(eqn_scale(EQN_PRES), abs(vpk))
      eqn_scale(EQN_RADIUS) = max(eqn_scale(EQN_RADIUS), abs(vrk))
      eqn_scale(EQN_TEMP) = max(eqn_scale(EQN_TEMP), abs(vtk))
      eqn_scale(EQN_LUM) = max(eqn_scale(EQN_LUM), abs(lk - lq*mt))
      eqn_scale(EQN_MASS) = max(eqn_scale(EQN_MASS), abs(vmk))
      eqn_scale(EQN_INERT) = max(eqn_scale(EQN_INERT), abs(vik))
      eqn_scale(EQN_PHI) = max(eqn_scale(EQN_PHI), abs(phik))
      eqn_scale(EQN_OMEGA) = max(eqn_scale(EQN_PHI), abs(adam), abs(samt))
   end do
   !eqn_scale(EQN_H1) = 1.0d0
   !eqn_scale(EQN_HE4) = 1.0d0
   !eqn_scale(EQN_C12) = 1.0d0
   !eqn_scale(EQN_N14) = 1.0d0
   !eqn_scale(EQN_O16) = 1.0d0
   !eqn_scale(EQN_NE20) = 1.0d0
   !eqn_scale(EQN_MG24) = 1.0d0
   eqn_scale(EQN_SUMX) = max(eqn_scale(EQN_H1), eqn_scale(EQN_HE4), eqn_scale(EQN_C12),  &
         eqn_scale(EQN_O16), eqn_scale(EQN_NE20), eqn_scale(EQN_N14),  &
         eqn_scale(EQN_MG24))

   do i = 1, NVAR
      if (eqn_scale(i) == 0.0d0) eqn_scale(i) = 1.0d0
   end do
end subroutine scale_eqn





! ------------------------------------------------------------------------------
!  CHECKS
!   Set some very small, or negative, values to zero
!   In particular, this deals with composition variables
! ------------------------------------------------------------------------------
!  Input:
!     COMMON H(:,:), DH(:,:)  - A solution step from the solver
!  Output:
!     COMMON H(:,:), DH(:,:)  - small abundances set to zero, central boundary
!                               conditions set exactly to 0.
! ------------------------------------------------------------------------------
subroutine checks
   use real_kind
   use settings
   use indices
   use mesh
  
   implicit none
   integer, parameter :: ix(9) = (/VAR_H1, VAR_HE4, VAR_C12, VAR_N14, VAR_O16, VAR_NE20, VAR_MG24, VAR_SI28, VAR_FE56/)
   integer :: ik,i,ij,jstar

   ! Composition variables should not be negative and neither should the eccentricity.
   ! TODO: detect when composition variables "overflow" and "correct" for this?
   do ik = 1, kh
      do jstar = 1, isb
         do i = 1, 9
            ij = idx_for_star(ix(i), jstar)
            if ( h(ij, ik) + dh(ij, ik) < 0.0d0 ) then
               dh(ij, ik) = -h(ij, ik)
            end if
            if ( h(ij, ik) + dh(ij, ik) > 1.0d0 ) then
               dh(ij, ik) = 1.0d0-h(ij, ik)
            end if
         end do
      end do
      if ( h(VAR_ECC, ik) + dh(VAR_ECC, ik) < 0.0d0 ) dh(VAR_ECC, ik) = -h(VAR_ECC, ik)
   end do

   ! Sanity corrections on the rotation rate.
   ! Unfortunately this causes artifacts if rotation is supposed to be an eigenvalue, so
   ! disabled until I can think of the correct way to do this.
   !do ik = 2, kh-1
   !   do jstar = 1, isb
   !      ij = idx_for_star(VAR_OMEGA, jstar)
   !      if ( h(ij, ik) + dh(ij, ik) < 0.0 ) then
   !         dh(ij, ik) = -h(ij, ik) / 1.5d0
   !      end if
   !   end do
   !end do

   do ik = 1, kh-1
      do jstar = 1, isb
         ij = idx_for_star(VAR_MASS, jstar)
         if ( h(ij, ik) + dh(ij, ik) < 0.0d0 ) dh(ij, ik) = 0.0d0

         ij = idx_for_star(VAR_LNR, jstar)
         if ( h(ij, ik) + dh(ij, ik) < 0.5 * log(ct(8)) ) dh(ij, ik) = 0.0d0
      end do
   end do

   ! Surface boundary condition: primary mass (eigenvalue)
   !dh(VAR_PMASS, 1:kh) = h(VAR_MASS, 1) + dh(VAR_MASS, 1) - h(VAR_PMASS, 1)
   h(VAR_PMASS, 1:kh) = h(VAR_MASS, 1) + dh(VAR_MASS, 1)
   dh(VAR_PMASS, 1:kh) = 0.0d0

   ! Central boundary conditions
   h(VAR_MASS, kh)  = 0.0d0
   h(VAR_LUM, kh)   = 0.0d0
   h(VAR_LNR, kh)   = 0.5 * log(ct(8))
   dh(VAR_MASS, kh) = 0.0d0
   dh(VAR_LUM, kh)  = 0.0d0
   dh(VAR_LNR, kh)  = 0.0d0
end subroutine checks





! ------------------------------------------------------------------------------
!  SMART_SOLVER:
!   Wrapper around SOLVER() that will try to recover when SOLVER() does not
!   converge before reporting back a non-convergence.
!   The various options that can be tried are controlled through ALLOW_*
!   constants, many of which can be set through init.dat. Many of these
!   options take the form of multiplying "difficult" terms in the equations
!   by a constant (say f) that is initially set to 0 and then slowly brought up
!   to 1 when SOLVER() has converged for a particular value of f.
! ------------------------------------------------------------------------------
!  Input:
!   COMMON H(:,:), DH(:,:) - The previous solution and the estimated changes
!   ITER                   - The number of allowed iterations in SOLVE()
!   IG(:)                  - The list of independant variables and equations
!   KT5                    - Print solver information after this many iterations
!  Output:
!   JO                     - Return value; 0 means converged ok
!   COMMON DH(:,:)         - Corrections for the current timestep
! ------------------------------------------------------------------------------
subroutine smart_solver ( iter, ig, kt5, jo )
   use real_kind
   use mesh
   use mesh_enc
   use settings
   use control

   implicit none
   integer, intent(in) :: iter, kt5, ig(130)
   integer, intent(out) :: jo
   real(double) :: fh(NVAR,NM), dfh(NVAR,NM)
   real(double) :: backup, bckavmu_smooth
   integer :: kt6
   integer, parameter :: verbosity = 0

  ! Store H and DH for possible reuse
   fh(1:NVAR, 1:kh) = h(1:NVAR, 1:kh)
   dfh(1:NVAR, 1:kh) = dh(1:NVAR, 1:kh)
   call solver ( iter, ig, kt5, jo )

  ! Detect circumstances that would force a smaller timestep; try to work around
  !  that first.

  ! Try reducing DH0
   if (jo /= 0) then
      kt6 = 3*iter+kt5
      if (verbosity>0) write(1, *) 'attempting to reduce dh0'
      backup = dh0
      dh0 = dh0 / 16.0d0
     ! Restore H and DH from before the last unsuccesful call to SOLVER
      jo = 0
      h(1:NVAR, 1:kh) = fh(1:NVAR, 1:kh)
      dh(1:NVAR, 1:kh) = dfh(1:NVAR, 1:kh)
      call solver ( iter, ig, kt6, jo )
      if (jo == 0) then
         if (verbosity>0) write(1, *) 'succes! :)'
      else
         if (verbosity>0) write(1, *) 'failed :('
      end if
      dh0 = backup
   end if

   if (jo /= 0) then
      kt6 = 3*iter+kt5
      if (verbosity>0) write(1, *) 'attempting to reduce dh0 more'
      backup = dh0
      dh0 = dh0 / 128.0d0
     ! Restore H and DH from before the last unsuccesful call to SOLVER
      jo = 0
      h(1:NVAR, 1:kh) = fh(1:NVAR, 1:kh)
      dh(1:NVAR, 1:kh) = dfh(1:NVAR, 1:kh)
      call solver ( iter, ig, kt6, jo )
      if (jo == 0) then
         if (verbosity>0) write(1, *) 'succes! :)'
      else
         if (verbosity>0) write(1, *) 'failed :('
      end if
      dh0 = backup
   end if

  ! Try setting DH to 0; this shouldn't usually help, but sometimes reducing
  !  the timestep seems to help convergence more than it seems it should,
  !  and it is actually clearing DH that helps.
   if (jo /= 0) then
      kt6 = 3*iter+kt5
      if (verbosity>0) write(1, *) 'trying to clear dh'
     ! Restore H from before the last unsuccesful call to SOLVER
      jo = 0
      h(1:NVAR, 1:kh) = fh(1:NVAR, 1:kh)
      dh(1:NVAR, 1:kh) = 0.0d0
      call solver ( iter, ig, kt6, jo )
      if (jo == 0) then
         if (verbosity>0) write(1, *) 'succes! :)'
      else
         if (verbosity>0) write(1, *) 'failed :('
      end if
   end if

   if (allow_underrelaxation) then
     ! Try reducing the convective mixing
      if (jo /= 0 .and. verbosity>0) write (1, *) 'attempting to underrelax convective mixing'
      call fudged_solver ( iter, ig, kt5, jo, mixing_fudge, fh, dfh)

     ! Try smoothing out the thermal energy contribution to the luminosity equation
      if (jo /= 0 .and. verbosity>0) write (1, *) 'attempting to underrelax thermal energy release'
      call fudged_solver ( iter, ig, kt5, jo, luminosity_fudge, fh, dfh)
   end if

   if (allow_mdotrelaxation) then
      if (jo /= 0 .and. verbosity>0) write (1, *) 'attempting to underrelax cmi term'
      call fudged_solver ( iter, ig, kt5, jo, mdot_smooth, fh, dfh)
   end if

   if (allow_avmurelaxation .and. .not. use_previous_mu) then
      if (jo /= 0 .and. verbosity>0) write (1, *) 'attempting to underrelax mean molecular weight'
      call fudged_solver ( iter, ig, kt5, jo, avmu_smooth, fh, dfh)
   end if

  ! Various ways to improve the solution that we found, eg, things that were
  ! computed implicitly from values at the current timestep (for stbility
  ! reasons).
  ! Store previous values of some of these settings (to be restored when the
  ! solution has been polished; these options need to stay swtiched on/off
  ! during the solution polishing step)
   bckavmu_smooth = avmu_smooth

  ! If converged with a value of mu from the previous timestep, try with the
  ! value of mu from the current timestep now
   if (allow_avmurelaxation .and. use_previous_mu .and. jo==0) then
      if (verbosity>0) write (1, *) 'switching to current molecular weight'
      avmu_smooth = 1.0d0
      call solver ( iter, ig, kt5/2, jo )
   end if

  ! Restore previous values of some settings
   avmu_smooth = bckavmu_smooth
end subroutine smart_solver





! ------------------------------------------------------------------------------
!  FUDGED_SOLVER
!   The iterator for SMART_SOLVER(), changing a `fudge factor' from 0 to 1.
!   When set to 0, the "fudge" is fully active (normally this means a
!   problematic term in the equations is suppressed).
!   This is called with various factors by the top level solver to attempt
!   relaxation of the solution. It works because FORTRAN passes arguments by
!   reference, so changing the value in this function will change the global.
! ------------------------------------------------------------------------------
!  Input:
!   COMMON H(:,:), DH(:,:) - The previous solution and the estimated changes
!   ITER                   - The number of allowed iterations in SOLVE()
!   IG(:)                  - The list of independant variables and equations
!   KT5                    - Print solver information after this many iterations
!   FACTOR                 - The "fudge factor" that is switched on or off
!   FH(:,:)                - A copy of H(:,:), in case we need to restore it
!   DFH(:,:)               - A copy of DH(:,:), in case we need to restore it
!  Output:
!   JO                     - Return value; 0 means converged ok
!   COMMON DH(:,:)         - Corrections for the current timestep
! ------------------------------------------------------------------------------
subroutine fudged_solver ( iter, ig, kt5, jo, factor, fh, dfh )
   use real_kind
   use mesh
   use mesh_enc
   use settings
  
   implicit none
   integer, intent(in) :: iter, kt5, ig(130)
   integer, intent(inout) :: jo
   real(double), intent(inout)  :: fh(NVAR,NM), dfh(NVAR,NM)
   real(double), intent(inout) :: factor
   real(double), parameter :: start_factor = 1.0d-10
   integer, parameter :: verbosity = 0      ! Increase for more output
   real(double) :: lh(NVAR,NM), dlh(NVAR,NM)
   integer :: kt6
   integer :: iftry
   real(double) :: prev_factor, ofactor

   if (jo /= 0) then
      kt6 = iter+1      ! Supress all status output from SOLVER
      if (verbosity > 1) kt6 = kt5

     ! Restore H and DH from before the unsuccesful call to SOLVER
      h(1:NVAR, 1:kh) = fh(1:NVAR, 1:kh)
      dh(1:NVAR, 1:kh) = dfh(1:NVAR, 1:kh)

      ofactor = factor
      factor = 0.0
      jo = 0
      call solver ( iter, ig, kt6, jo )

      factor = start_factor
      iftry = 0
      lh(1:NVAR, 1:kh) = h(1:NVAR, 1:kh)
      dlh(1:NVAR, 1:kh) = dh(1:NVAR, 1:kh)
      prev_factor = factor
      if (verbosity>0) write(1, *) 'starting relaxation'
      do while (factor<1.0d0 .and. jo == 0)
         factor = min(1.0d0, max(sqrt(factor), 1.1d0*factor))
         if (verbosity>0) write(1, *) 'factor now at', factor

        ! Restore H from before the last call to SOLVER; keep DH
         h(1:NVAR, 1:kh) = fh(1:NVAR, 1:kh)
         call solver ( iter, ig, kt6, jo )
         if (jo /= 0) then
           ! Not converged, backtrack and try with a reduced factor
            iftry = iftry+1
            if (iftry < 5) jo = 0
            dh(1:NVAR, 1:kh) = dlh(1:NVAR, 1:kh)
            factor = prev_factor*factor
         else
            iftry=0
            dlh(1:NVAR, 1:kh) = dh(1:NVAR, 1:kh)
            prev_factor = factor
         end if
      end do
      if (jo == 0) then
         if (verbosity>0) write(1, *) 'succes! :)'
      else
         if (verbosity>0) write(1, *) 'failed :('
      end if
     ! Set the factor back to 1.0, which means to use normal physics
      factor = ofactor
   end if
end subroutine fudged_solver




! ------------------------------------------------------------------------------
!  RELAX_TIMESTEP
!   When SOLVER() fails to converge and we need to cut back the timestep,
!   try using a reduced timestep to get a new estimate for DH(:,:) and then
!   increase the timestep again to its current value.
!   If succesful, we may avoid having to reduce the timestep entirely,
!   which speeds up the calculation of the rest of the evolution.
! ------------------------------------------------------------------------------
!  Input:
!   COMMON H(:,:), DH(:,:) - The previous solution and the estimated changes
!   ITER                   - The number of allowed iterations in SOLVE()
!   IG(:)                  - The list of independant variables and equations
!   KT5                    - Print solver information after this many iterations
!   JO                     - Convergence flag, if 0 we don't do anything here
!   DTY                    - The current value of the timestep (in years)
!  Output:
!   JO                     - Return value; 0 means converged ok
!   COMMON DH(:,:)         - Corrections for the current timestep
! ------------------------------------------------------------------------------
subroutine relax_timestep ( iter, ig, kt5, jo, dty )
   use real_kind
   use mesh
   use mesh_enc
   use settings
   use control
   use extrapolate_dh
   use constants
   use test_variables
  
   implicit none
   integer, intent(in) :: iter, kt5, ig(130)
   integer, intent(inout) :: jo
   real(double), intent(in) :: dty
   integer, parameter :: verbosity = 0      ! Increase for more output
   real(double) :: fh(NVAR,NM), dfh(NVAR,NM)
   real(double) :: new_dty

   if (.not. recycle_timestep .or. jo == 0) return

  ! Store H and DH for possible reuse
   fh(1:NVAR, 1:kh) = h(1:NVAR, 1:kh)
   dfh(1:NVAR, 1:kh) = dh(1:NVAR, 1:kh)
   new_dty = ct3*dty

   if (verbosity > 0) write (1, *) 'cycling timestep'
   do while (.true.)
      jo = 0
      dt = new_dty*csy
      h(1:NVAR, 1:kh) = fh(1:NVAR, 1:kh)
      call smart_solver ( iter, ig, kt5, jo )
      if (jo /= 0) then
         dt = dty*csy
         h(1:NVAR, 1:kh) = fh(1:NVAR, 1:kh)
         dh(1:NVAR, 1:kh) = dfh(1:NVAR, 1:kh)
         return
      end if
      if (new_dty >= dty) exit
      new_dty = min(dty, 1.05*dty)
   end do
   if (verbosity > 0) write (1, *) 'success'
end subroutine relax_timestep


