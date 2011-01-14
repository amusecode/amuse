module precalc_eostate
   use real_kind
   use mesh, only: nm, nvar
   use eostate_types
   use atomic_data
   
   implicit none
   type(eostate), private, save, allocatable :: eos(:,:)
   type(abundance), private, save, allocatable :: abund(:,:)
   logical, private, save :: have_precomputed_values(0:NVAR)

   
   contains

   
   subroutine clear_cache
      use real_kind
      
      implicit none
      have_precomputed_values(:) = .false.
   end subroutine
   
   
   
   subroutine initialise_eostate_precalc
      use real_kind
      
      implicit none
      if (.not. allocated(eos)) allocate(eos(0:NVAR, NM))
      if (.not. allocated(abund)) allocate(abund(0:NVAR, NM))
   end subroutine initialise_eostate_precalc
   
   
   
   function get_precomputed_state(jk, ji, out_abund, out_eos)
      use real_kind
      
      implicit none
      integer, intent(in) :: jk, ji
      type(abundance), intent(out) :: out_abund
      type(eostate), intent(out) :: out_eos
      logical :: get_precomputed_state

      get_precomputed_state = .false.
      if (.not. allocated(eos) .or. ji < 0) return

      out_eos = eos(ji, jk)
      out_abund = abund(ji, jk)

      get_precomputed_state = have_precomputed_values(ji)
   end function
   
   
   subroutine precompute_eos
      use real_kind
      use mesh
      use settings
      use neutrinos
      use constants
      use extra_elements
      
      implicit none
      real(double) :: var(11)
      real(double) :: dvar(11)
      real(double) :: dx, vx, dvx
      integer, parameter :: var_perm(11) = (/ 1, 2, 5, 9, 10, 16, 3, 11, NMg24, NSi28, NFe56 /)
      integer :: jk, n

      have_precomputed_values(:) = .false.
      call initialise_eostate_precalc
      ! Precompute EoS values
      do jk = 1, kh
         ! Function values
         forall (n = 1:11)
            dvar(n) = dh(var_perm(n), jk)
            var(n) = h(var_perm(n), jk) + dvar(n)
         end forall
         ! Compute function values
         call statef (var(1),var(2),var(3:11), abund(0, jk),eos(0, jk))
         have_precomputed_values(0) = .true.
         ! Varying independent variables in turn to get derivatives
         do n = 1, 8
            vx = var(n)
            dvx = dvar(n)

            dx = sign(dh0*max(abs(vx), 1.0d0), dvx)
            var(n) = vx + dx
            dvar(n) = dvx + dx
            call statef (var(1),var(2),var(3:11),  &
                              abund(var_perm(n), jk),eos(var_perm(n), jk))
            have_precomputed_values(var_perm(n)) = .true.

            var(n) = vx
            dvar(n) = dvx
         end do
      end do

      ! Precompute nuclear reaction rates
      do jk = 1, kh
         do n = 0, 8
            call nucrat (eos(var_perm(n), jk)%at,  &
                           abund(var_perm(n), jk),eos(var_perm(n), jk))
         end do
      end do

   end subroutine

end module precalc_eostate

! ------------------------------------------------------------------------------
!  STATEL
!   Cahching of EoS and reaction rate calculations, for funcs1.
!   These values don't need to be recalculated if the EoS does not depend
!   on the variable that was varied.
!   TODO: make estate independent of opacity and neutrino rates by doing the
!   opacity-related calculations and the neutrino rates here, which makes
!   estate a "pure" EoS routine that doesn't require extra input.
!
!   CAUTION: values of AF, AT and /ABUND/ are NOT used to decide whether
!   the cached values can be used or not. The values last stored in the
!   cache will be returned based on JI, independant of AF and AT.
! ------------------------------------------------------------------------------
!  Input:
!     JK       - The meshpoint at which the value was set
!     JI       - The variable that was varied, see FUNCS1.
!                JI<=0 means that the EoS state is calculated and stored
!     AF       - ln f, degeneracy parameter
!     AT       - ln T, temperature
!     XA(:)    - Abundances, not used here put passed through to ESTATE
!     JSTAR    - The number of the star in the binary, 1 or 2.
!  Output:
!     EOS      - EoS variables and reaction rates for current parameters
!     ABUND    - Current number densities, average charge and ionisation state
! ------------------------------------------------------------------------------
subroutine statel ( jk, ji, af, at, xa, Jstar, abund, eos )
   use real_kind
   use eostate_types
   use precalc_eostate
   use neutrinos
   use constants
   
   implicit none
   integer, intent(in) :: jk, ji, Jstar
   real(double), intent(in) :: af, at
   real(double), intent(in) :: xa(9)
   type(eostate), intent(out) :: eos
   type(abundance), intent(out) :: abund
   
   integer :: i
   type(eostate), save :: eos_cache(2)
   type(abundance), save :: abundance_cache(2)
   
   ! IFN contains entries for the independent variables: a 2 is listed
   !  if a variation in the corresponding variable requires that the
   !  equation of state is reevaluated, a 1 means that this is not
   !  nescessary and the cached result can be used.
   ! JI=0,-1,-2 are `dummy variables'. 0 means that the EOS is calculated
   !  and stored, -1 and -2 indicate the request came from
   !  printb or remesh respectively.
   integer, parameter :: ifn(-2:44) = (/2, 2, 0,                     &
                     2, 2, 2, 1, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 2, &
                     1, 1, 1, 1, 1, 1, 1, 1,                         &
                     2, 2, 2, 1, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 2, &
                     2, 2, 2, 1/)
   ! What star a given variable belongs to, *1 or *2 (or neither/both)
   integer, parameter :: instar(-2:44) = (/ 0, 0, 0,                 &
                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                     0, 0, 0, 0, 0, 0, 0, 0,                         &
                     2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
                     1, 1, 1, 1/)

   ! Do we have this value pre-computed?
   if (get_precomputed_state(jk, ji, abund, eos)) return

   ! No, but maybe the EoS has not changed?
   ! NB: FORTRAN logical operators .and., .or. etc. don't short-circuit,
   ! apparently, so we need to use nested if's to short-circuit things be
   ! hand...
   if (ji <= 44) then
      if (ifn(ji) == 1) then     ! EoS does not depend on this variable
         if (get_precomputed_state(jk, 0, abund, eos)) return
      end if
   end if

   ! The state is not cached, the most likely reason is that we're not
   ! actually using the precompute cache. Evaluate the EoS as normal,
   ! using the local (non-thread safe!) cache as needed.

   !  Don't compute eos for *1 if *2 has changed and vice versa.
   !if (ji<1 .or. (Jstar==1.and.ji<=16).or.(Jstar==2.and.ji>=24) .or. ji>44) then
   if ( ji > 44 .or. instar(ji) == 0 .or. instar(ji) == Jstar ) then
      i = 2
      if (ji <= 44) i = ifn(ji)
      if ( i /= 1 ) then
         call statef ( af, at, xa, abund, eos )
         call nucrat ( at, abund, eos )
         if ( i > 0 ) return
         eos_cache(Jstar) = eos
         abundance_cache(Jstar) = abund
         return
      end if
   end if
   eos = eos_cache(Jstar)
   abund = abundance_cache(Jstar)
end subroutine statel

