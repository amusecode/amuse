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
   use extra_elements
   use test_variables
   use funcs1_interface
   
   implicit none
   real(double), intent(out) :: eqn_scale(1:NVAR), sbc_scale(1:NVAR)

   ! COMMON block INF: pass input to FUNCS1, get return values from FUNCS1
   real(double) :: var(nvar), dvar(nvar), fn1(nfunc)
   real(double) ::  bcp, bct, vp, vpk, vr, vrk, vt, vtk, vl, lk
   real(double) :: lq, mt, vm, vmk, sg, wt, x1, X1T, x16, X16T, x4, X4T, x12
   real(double) :: X12T, x20, X20T, bcm, vi, vik, vphi, phik, bcf, bcs, bcph
   real(double) :: x14, X14T, avmu, sgth, omega, omegat, sgam, si
   real(double) :: x24, X24T

   ! Local variables:
   integer :: ik, i
   integer, parameter :: Jstar = 1  ! Only do single star for now

   ! Initialise with default settings: all equations have a natural value of 1
   eqn_scale(:) = 0.0d0
   sbc_scale(:) = 1.0d0

  do ik=kh, 1, -1
      !        Copy current value of variables: single star
      dvar(1:16) = dh(1 + 24*(Jstar-1):16 + 24*(Jstar-1), ik)
      var(1:16) = h(1 + 24*(Jstar-1):16 + 24*(Jstar-1), ik) + dvar(1:16)
      dvar(41:40+nxfstar) = dh(41+nxfstar*(Jstar-1):40+nxfstar*(Jstar-1)+nxfstar, ik)
      var(41:40+nxfstar) = h(41+nxfstar*(Jstar-1):40+nxfstar*(Jstar-1)+nxfstar, ik) + dvar(41:40+nxfstar)
      
      !        Binary parameters
      dvar(17:24) = dh(17:24, ik)
      var(17:24) = h(17:24, ik) + dvar(17:24)
      
      !        Compute function values
      call funcs1 ( ik, -1, var(:), dvar(:), fn1(:) )
      bcp  = fn1(58*(Jstar-1) + 1)
      bct  = fn1(58*(Jstar-1) + 2)
      vp   = fn1(58*(Jstar-1) + 3)
      vpk  = fn1(58*(Jstar-1) + 4)
      vr   = fn1(58*(Jstar-1) + 5)
      vrk  = fn1(58*(Jstar-1) + 6)
      vt   = fn1(58*(Jstar-1) + 7)
      vtk  = fn1(58*(Jstar-1) + 8)
      vl   = fn1(58*(Jstar-1) + 9)
      lk   = fn1(58*(Jstar-1) + 10)
      lq   = fn1(58*(Jstar-1) + 11)
      mt   = fn1(58*(Jstar-1) + 12)
      vm   = fn1(58*(Jstar-1) + 13)
      vmk  = fn1(58*(Jstar-1) + 14)
      sg   = fn1(58*(Jstar-1) + 15)
      wt   = fn1(58*(Jstar-1) + 16)
      x1   = fn1(58*(Jstar-1) + 17)
      x1t  = fn1(58*(Jstar-1) + 18)
      x16  = fn1(58*(Jstar-1) + 19)
      x16t = fn1(58*(Jstar-1) + 20)
      x4   = fn1(58*(Jstar-1) + 21)
      x4t  = fn1(58*(Jstar-1) + 22)
      x12  = fn1(58*(Jstar-1) + 23)
      x12t = fn1(58*(Jstar-1) + 24)
      x20  = fn1(58*(Jstar-1) + 25)
      x20t = fn1(58*(Jstar-1) + 26)
      bcm  = fn1(58*(Jstar-1) + 27)
      vi   = fn1(58*(Jstar-1) + 28)
      vik  = fn1(58*(Jstar-1) + 29)
      vphi = fn1(58*(Jstar-1) + 30)
      phik = fn1(58*(Jstar-1) + 31)
      bcf  = fn1(58*(Jstar-1) + 32)
      bcs  = fn1(58*(Jstar-1) + 33)
      bcph = fn1(58*(Jstar-1) + 34)
      x14  = fn1(58*(Jstar-1) + 35)
      x14t = fn1(58*(Jstar-1) + 36)
      avmu = fn1(58*(Jstar-1) + 37)
      sgth = fn1(58*(Jstar-1) + 38)
      omega = fn1(58*(Jstar-1) + 39)
      omegat = fn1(58*(Jstar-1) + 40)
      sgam = fn1(58*(Jstar-1) + 41)
      si   = fn1(58*(Jstar-1) + 42)

      !> \todo FIXME: these are stored in the binary section
      x24  = fn1(49)
      x24t = fn1(50)
      !        Composition equations
      eqn_scale(1) = max(eqn_scale(1), dabs(X1T))
      eqn_scale(2) = max(eqn_scale(2), dabs(X16T))
      eqn_scale(3) = max(eqn_scale(3), dabs(X4T))
      eqn_scale(4) = max(eqn_scale(4), dabs(X12T))
      eqn_scale(5) = max(eqn_scale(5), dabs(X20T))
      eqn_scale(EN14) = max(eqn_scale(EN14), dabs(X14T))
      eqn_scale(EMG24) = max(eqn_scale(EMG24), dabs(X24T))

      !        Structure equations
      eqn_scale(6) = max(eqn_scale(6), dabs(vpk))
      eqn_scale(7) = max(eqn_scale(7), dabs(vrk))
      eqn_scale(8) = max(eqn_scale(8), dabs(vtk))
      eqn_scale(9) = max(eqn_scale(9), dabs(lk - lq*mt))
      eqn_scale(10) = max(eqn_scale(10), dabs(vmk))
      eqn_scale(11) = max(eqn_scale(11), dabs(vik))
      eqn_scale(12) = max(eqn_scale(12), dabs(phik))
   end do
   eqn_scale(1) = 1.0d0
   eqn_scale(2) = 1.0d0
   eqn_scale(3) = 1.0d0
   eqn_scale(4) = 1.0d0
   eqn_scale(5) = 1.0d0
   eqn_scale(EN14) = 1.0d0
   eqn_scale(EMG24) = 1.0d0
   eqn_scale(ESUMX) = max(eqn_scale(1), eqn_scale(2), eqn_scale(3),  &
         eqn_scale(4), eqn_scale(5), eqn_scale(EN14),  &
         eqn_scale(EMG24))

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
  use mesh
  
  implicit none
  integer :: ix(20),ik,i,ij

  data ix/3, 5, 9, 10, 11, 16, 18, 27, 29, 33, 34, 35, 40, 7*3/

  ! Composition variables - should not be negative!
  do ik = 1, kh
     do i = 1, 6
        ij = ix(i)
        if ( h(ij, ik) + dh(ij, ik) < 0.0 ) then
           dh(ij, ik) = -h(ij, ik)
        end if
     end do
  end do

  ! Other variables
  do ik = 1, kh
     do i = 7, 13
        ij = ix(i)
        if ( h(ij, ik) + dh(ij, ik) <= 1.0d-12 ) then
           dh(ij, ik) = -h(ij, ik)
        end if
     end do
  end do

  ! Central boundary conditions
  h(4, kh) = 0.0d0
  h(8, kh) = 0.0d0
  dh(4, kh) = 0.0d0
  dh(8, kh) = 0.0d0
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
  real(double) :: bckdh0, bckavmu_smooth
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
     bckdh0 = dh0
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
     dh0 = bckdh0
  end if

  if (jo /= 0) then
     kt6 = 3*iter+kt5
     if (verbosity>0) write(1, *) 'attempting to reduce dh0 more'
     bckdh0 = dh0
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
     dh0 = bckdh0
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

  if (allow_amadvrelaxation .and. jo /= 0) then
     if (jo /= 0 .and. verbosity>0) write (1, *) 'attempting to underrelax angular momentum transport'
     call fudged_solver ( iter, ig, kt5, jo, amadv_smooth, fh, dfh)
     if (jo == 0) then
        write (1, *) 'succes! [Angular momentum transport]'
     else
        write (1, *) 'failure! [Angular momentum transport]'
     end if
  end if

  if (allow_overrelaxation) then
     ! Try boosting the convective mixing (!) using some weird scheme that seems to
     !  work well in some cases.
     if (jo /= 0 .and. verbosity>0) write (1, *) 'attempting to overrelax convective mixing'
     call fudged_solver ( iter, ig, kt5, jo, mixing_boost, fh, dfh)
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
     kt6 = iter+30+kt5+1    ! Supress all status output from SOLVER
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
        if (factor < 0.8d0) then
           factor = sqrt(factor)
        else
           factor = min(1.0d0, 1.1d0*factor)
        end if
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


