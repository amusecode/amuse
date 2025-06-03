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
   use neutrinos
   use constants
   use indices
   use eostate_types
   use atomic_data
   use equation_of_state
   
   implicit none
   integer, intent(in) :: jk, ji, Jstar
   real(double), intent(in) :: af, at
   real(double), intent(in) :: xa(9)
   type(eostate), intent(out) :: eos
   type(abundance), intent(out) :: abund
   
   integer :: i, istar
   type(eostate), save :: eos_cache(2)
   type(abundance), save :: abundance_cache(2)
   logical :: calculate, store

   calculate = .true.
   store = .true.

   if (ji < 0) then
      ! The call came from printb or remesh - don't cache the results (no point)
      store = .false.
   elseif (ji > 0) then
      ! If the EoS depends on the changed variable, then calculate it
      istar = 1
      i = idx_primary(ji)
      if (ji /= i) istar = 2

      ! We only need to recalculate the EoS if:
      !  1. The changed variable is a composition or thermodynamic variable
      !  and
      !  2. The changed variable is for the current star,
      !  or 
      !  3. We are at the surface meshpoint (where the stars are linked by the surface boundary conditions)
      calculate = .false.
      if (is_abundance(i) .or. i == VAR_LNT .or. i == VAR_LNF) then
         if (istar == jstar .or. jk == 1) then
            calculate = .true.
            store = .false.
         end if
      end if
   end if


   if (calculate) then
      call statef ( af, at, xa, abund, eos )
      call nucrat ( at, abund, eos )
      if (store) then
         eos_cache(Jstar) = eos
         abundance_cache(Jstar) = abund
      end if
      return
   end if
   eos = eos_cache(Jstar)
   abund = abundance_cache(Jstar)
end subroutine statel

