module ionisation

! ------------------------------------------------------------------------------
! FIND_IONISATION_LEVELS:
!  Given elemental (or isotopic) abundances, the local temperature and the
!  (modified) electron chemical potential, calculate the abundance of the
!  different ionisation stages.
! ------------------------------------------------------------------------------
!  Input:
!     nelem       - The number of elements/isotopes to consider
!     Zmax        - The maximum charge of any element in the list
!     T           - The local temperature [K]
!     dV          - The electron chemical potential, corrected for pressure
!                   ionisation and Coulomb interactions [kT]
!     X           - Abundances of the different isotopes. Could be mass
!                   fractions, number fractions, or something else [arbitrary]
!     com(Zmax+1) - Degeneracy of the ground state for different numbers of electrons
!     kzn(nelem)  - The total charge of each atomic species [e]
!     amu(nelem)  - The atomic mass of each atomic species, used only for output [AMU]
!     chi(Zmax,Zmax) - Cumulative ionisation potentials for each species.
!                      The j-th ionisation level of species i is indexed as
!                      chi(j, kzn(i)) [eV]
!  Output:
!     ns          - The number of distinct ionisation stages calculated.
!                   Typically, this will be Zmax*(Zmax+3)/2, unless nelem > Zmax.
!     xai()       - Abundances for each ionisation stage. Units are the same
!                   as for X, and abundances of ionisation stages for each
!                   species i sum to X(i).
!                   FIXME: exactly why is this convenient?
!                   Only entries 1:ns will be filled, but the array should be
!                   big enough to hold Zmax*(Zmax+3)/2 entries.
!     kzi()       - Atomic number (kzn) for each ionisation state; to identify
!                   what species it belongs to.
!                   Only entries 1:ns will be filled, but the array should be
!                   big enough to hold Zmax*(Zmax+3)/2 entries.
!     aamu()      - Mass for this ionisation state. Passed through from amu().
!                   Only entries 1:ns will be filled, but the array should be
!                   big enough to hold Zmax*(Zmax+3)/2 entries.
!     zz()        - Ionisation stage (charge) this entry corresponds to:
!                    0 - neutral (no free electrons due to this isotope)
!                    1 - singly ionised (1 free electron due to this isotope)
!                    2 - double ionised (2 free electrons due to this isotope)
!                    etc.
!                   Only entries 1:ns will be filled, but the array should be
!                   big enough to hold Zmax*(Zmax+3)/2 entries.
! ------------------------------------------------------------------------------

contains

subroutine find_ionisation_levels(nelem,  Zmax, t, dv, x,   com, kzn, amu, chi,   ns, xai, kzi, aamu, zz) 
   use constants, only: cevb
   use real_kind
   implicit none
!  Input parameters
   integer,      intent(in)  :: nelem, Zmax
   real(double), intent(in)  :: x(nelem), t, dv
!  Constant data
   real(double), intent(in)  :: com(Zmax+1)
   integer, intent(in)       :: kzn(nelem)
   real(double), intent(in)  :: amu(nelem)
   real(double), intent(in)  :: chi(Zmax, Zmax)
!  Output
   integer,      intent(out) :: ns
   real(double), intent(out) :: xai(:)
   integer,      intent(out) :: kzi(:)
   integer,      intent(out) :: zz(:)
   real(double), intent(out) :: aamu(:)

   real(double)              :: va(zmax), ha(zmax), sum_ha, ti
   integer                   :: i, j, n, n_first, n_last

   ti = cevb/T

   n = 1
   kzi = 0
   xai = 0
   n_last = 0
   do i = 1, nelem
      ! Record the first ionisation state for this element
      n_first = n

      ! Neutral state
      sum_ha =1.0d0
      xai(n) = 1.0
      zz(n) = 0
      n = n + 1

      if (kzn(i) > 0) then

         ! Singly ionised state
         va(1) = -chi(1, kzn(i))*ti + dv
         ha(1) = fxp(va(1)) * com(kzn(i)) / com(kzn(i)+1)
         sum_ha = sum_ha + ha(1)

         xai(n) = ha(1)
         zz(n) = 1
         n = n + 1
         do j = 2, kzn(i)
            va(j) = -chi(j,kzn(i))*ti + j*dv
            ha(j) = ha(j-1)*fxp(va(j) - va(j-1)) * com(kzn(i)+1-j) / com(kzn(i)+2-j)
            sum_ha = sum_ha + ha(j)

            xai(n) = ha(j)
            zz(n) = j
            n = n + 1
         end do
      end if
      ! Record the last ionisation state for this element
      n_last = n-1

      ! Normalise all ionisation states so that they sum up to the elemental abundance X
      xai(n_first:n_last) = xai(n_first:n_last) * x(i) / sum_ha
      kzi(n_first:n_last) = kzn(i)
      aamu(n_first:n_last) = amu(i)
   end do
   ns = n_last - 1

contains

   ! FXP: Fudged exponential function, used to avoid too large or too small numbers
   ! in the ionisation state of the Saha equation. Usesprecomputed values of
   ! the limiting exponentials.
   function fxp(x)
      use real_kind

      implicit none
      real(double) :: fxp
      real(double), intent(in) :: x
      real(double), parameter :: fxp_low = -50.0d0
      real(double), parameter :: fxp_high = 50.0d0
      real(double), parameter :: low = 1.928749847963917820563444d-22
      real(double), parameter :: high = 5.184705528587072045056000d+21
      if (x>fxp_high) then
         fxp = high
         return
      else if (x>fxp_low) then
         fxp = exp(x)
         return
      end if
      fxp = low
   end function fxp

end subroutine

end module

