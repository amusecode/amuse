! ------------------------------------------------------------------------------
! DIFFUSION_COEFFICIENTS
! Compute atomic diffusion coefficients, according to Paquette et al. 1986
!
! This function only calculates the first order coefficients at the moment,
! but calculating the second order corrections (as well as the coefficient
! for thermodiffusion) is not difficult.
! ------------------------------------------------------------------------------
! Input:
!  * rho:      Density [g/cm^3]
!  * T:        Temperature
!  * abund:    Abundance struct
! Output:
!  * Ddiff(9): Diffusion coefficients [cm^2/s] for the 9 isotopes H, He, C,
!              n, O, Ne, Mg, Si and Fe.
! ------------------------------------------------------------------------------
subroutine diffusion_coefficients(rho, T, abund, Ddiff, Nref)
   use real_kind
   use constants
   use paquette_coefficients
   use atomic_data
   use eostate_types

   implicit none
  ! Input and output variables
   real(double), intent(in) :: rho, T
   type(abundance), intent(in) :: abund
   real(double), intent(out) :: Ddiff(9)
   integer, intent(in) :: Nref            ! Reference abundance (H)
  ! Local variables
  ! Lengthscales for the plasma
   real(double) :: lambda_I2, lambda_D2, lambda2
  ! Baryon density
   real(double) :: Nb
  ! Convenience, to precalculate as much as possible
   real(double) :: kT, gamma2
  ! Collision integrals Omega (eqn. (18) in Paquette et al)
   real(double) :: OMEGA1(3), OMEGA22
  ! Dimensionless collision integrals (eqn. (65) in Paquette et al)
   real(double) :: F1(3), F22
  ! Local variables
   real(double) :: psi_st, eps_st, gamma_st2, A, e_psi_st
   real(double) :: dpsi_n1, dpsi_n
   integer :: i, j, n
  ! Local copy of struct variables
   real(double) :: NA(9), NIO, NZZ, AVM, Ne

   na = abund%na
   nio = abund%nio
   nzz = abund%nzz
   avm = abund%avm
   na = abund%na
   Ne = abund%ne

   Nb = rho/(AVM*AMU)
   kT = BOLTZM*T

  !     Typical distance between ions (squared)
   lambda_I2 = (3.0/(CPI4*Nb*NIO))**(2.0D0*C3RD)
  !     Debye length (squared) 
   lambda_D2 = (kT/(CPI4 * ECHAR**2 * Nb * (NZZ+Ne)))
  !     Use max of lambda_i and lambda_D as typical distance
   lambda2 = MAX(lambda_D2, lambda_I2)
  !     Precompute some factors
   gamma2 = (4.0*kT/(ECHAR**2*KZN(Nref)))**2*lambda2
   do i=1, 9
      gamma_st2 = gamma2/KZN(i)**2
     !        Reduced mass of the two particles
      A = (CAN(i)*CAN(Nref)/(CAN(i)+CAN(Nref)))
     !        Units of Omega: length^2 * velocity
      eps_st = CPI4*lambda2/gamma_st2 * SQRT(2.0D0*kT/(CPI4 * AMU*A))
     !        Find interpolation interval for dimensionless collision integral
      e_psi_st = log(1.0D0+gamma_st2)
      psi_st = log(e_psi_st)
     !        Evaluate the collision integral, for repulsive coefficients
      if (psi_st < 3.0D0) then
        !           Use spline interpolation to evaluate the collision integrals
        !           Determine interval in the table
         n = MIN(50, 1+floor((psi_st +7d0)/0.2d0))
         dpsi_n1 = (-7.0D0 + 0.2D0*(n+1)) - psi_st
         dpsi_n = psi_st - (-7.0D0 + 0.2D0*n)
         do j=1,3
            F1(j) = dexp(DC(1,n,j)*dpsi_n1**3 + DC(2,n,j)*dpsi_n**3  &
                + DC(3,n,j)*dpsi_n1    + DC(4,n,j)*dpsi_n)
         end do
         F22 = dexp(DD(1,n)*dpsi_n1**3 + DD(2,n)*dpsi_n**3  &
             + DD(3,n)*dpsi_n1    + DD(4,n)*dpsi_n)
      else
         F1(1) = 1.00141D0*e_psi_st - 3.18209D0
         F1(2) = 0.99559D0*e_psi_st - 1.29553D0
         F1(3) = 1.99814D0*e_psi_st - 0.64413D0
         F22   = 1.99016D0*e_psi_st - 4.56958D0
      end if
      OMEGA1(:) = eps_st*F1(:)
      OMEGA22   = eps_st*F22

     !        Diffusion coefficient
      Ddiff(i) = 3.0*kT/(16.0*rho*(NA(Nref)+NA(i))*A*OMEGA1(1))
   end do
end subroutine diffusion_coefficients



! ------------------------------------------------------------------------------
! COMPUTE_FULL_DIFFUSION_COEFFS
! Compute atomic diffusion coefficients, according to Paquette et al. (1986)
! This is essentially the same as the above subroutine, but it computes
! extra terms needed for the Burgers formalism.
! The output are the resistance coefficients in Burgers equations K_ij,
! z_ij, z'_ij and z''_ij. It is also possible to output diffusion
! coefficients D_ij and thermal diffusion coeffcient A_th used in
! Cowling&Chapman formalism, note Ath(NN,i) is Ath_ei.
! ------------------------------------------------------------------------------
! Input:
! (rho         - Density [g/cm^3])
!  T           - Temperature [k]
!  NN          - Number of species (NN = electrons)
!  CAN(i)      - Mass numbers for species [amu]
!  KZN(i)      - Average charge of species [e]
!  NA(i)       - Number densities of species [cm^-3]
! Output:
!  Kdiff(i,j)  - Coefficient for Burgers equations
!  Zdiff(i,j)  - Coefficient for Burgers equations
!  Zdiff1(i,j) - Coefficient for Burgers equations
!  Zdiff2(i,j) - Coefficient for Burgers equations
! ------------------------------------------------------------------------------
!subroutine compute_full_diffusion_coeffs(rho, T, NN, CAN1, KZN1, NA1, Kdiff, Zdiff, Zdiff1, Zdiff2)  !Rho not used
subroutine compute_full_diffusion_coeffs(T, NN, CAN1, KZN1, NA1, Kdiff, Zdiff, Zdiff1, Zdiff2)
   use real_kind
   use constants
   use paquette_coefficients

   implicit none
  ! Input and output variables
   integer, intent(in) :: NN
   real(double), intent(in) :: T, CAN1(NN), KZN1(NN), NA1(NN)
   real(double), intent(out) :: Kdiff(NN,NN), Zdiff(NN,NN), Zdiff1(NN,NN), Zdiff2(NN,NN)

  ! Local variables
  ! Lengthscales for the plasma
   real(double) :: lambda_I2, lambda_D2, lambda2

  ! Ion density
   real(double) :: NI

  ! Convenience, to precalculate as much as possible
   real(double) :: kT, gamma2

  ! Collision integrals Omega (eqn. (18) in Paquette et al)
   real(double) :: OMEGA1(3), OMEGA22, OMEGA2(NN)

  ! Dimensionless collision integrals (eqn. (65) in Paquette et al)
   real(double) :: F1(3), F22

  ! Local variables
  !real(double) :: Ddiff(NN,NN), Ath(NN, NN)
   real(double) :: psi_st, eps_st, gamma_st2, A, e_psi_st
   real(double) :: AA, BB, CC, EE, Ps, Pt, Pst, Qs, Qt, Qst, Ss
   real(double) :: St, Ms, Mt, Xs, Xt
  !real(double) :: Xe, Me, Pe, Se, Qe, Qet, DELTA
   real(double) :: dpsi_n1, dpsi_n
   integer :: i, j, n, Nref
   real(double) :: NZZ, Ne
   integer :: kz(NN)

   Ne = NA1(NN)
   NI = 0.0D0
   NZZ = 0.0D0
   do i = 1, NN-1
      NI = NI + NA1(i)
      NZZ = NZZ + KZN1(i)**2 * NA1(i)
   end do

  ! Charge, bearing in mind that some species might be neutral.
  ! In this case, we simply treat them as singly ionised.
  ! TODO: we could include a fudge-factor to scale the singly ionised cross
  ! section.
   kz = kzn1
   where (kz == 0) kz = 1

   kT = BOLTZM*T
  ! Typical distance between ions (squared)
   lambda_I2 = (3.0D0/(CPI4*NI))**(2.0D0/3.0D0)
  ! Debye length (squared)
   lambda_D2 = (kT/(CPI4 * ECHAR**2 * (Ne + NZZ)))
  ! Use max of lambda_i and lambda_d as typical distance
   lambda2 = MAX(lambda_D2, lambda_I2)

  ! First evaluate the collision intergrals OMEGA22_SS, that appear in Ps (eq (9) of Paquette 1986)
   do i = 1, NN
      gamma2 = (4.0D0*kT/(ECHAR**2*KZ(i)))**2*lambda2
      gamma_st2 = gamma2/KZ(i)**2

     ! Reduced mass of the two particles [amu]
      A = 0.5d0 * CAN1(i)

     ! Units of Omega: length^2 * velocity
      eps_st = CPI4*lambda2/gamma_st2 * SQRT(2.0D0*kT/(CPI4 * AMU*A))

     ! Find interpolation interval for dimensionless collision integral
     ! Use boundary value if out of range (may never happen)
      e_psi_st = log(1.0D0+gamma_st2)
      psi_st = MAX(-7.0d0, log(e_psi_st))

     ! Evaluate the collision integral
      if (psi_st <= 3.0D0) then
        ! Use spline interpolation to evaluate the collision integrals
        ! Determine interval in the table
         n = MIN(50, 1+floor((psi_st +7d0)/0.2d0))
         dpsi_n1 = (-7.2D0 + 0.2D0*(n+1)) - psi_st
         dpsi_n = psi_st - (-7.2D0 + 0.2D0*n)

        ! repulsive potential (ion-ion or electron-electron)
         F22 = dexp(DD(1,n)*dpsi_n1**3 + DD(2,n)*dpsi_n**3  &
             + DD(3,n)*dpsi_n1    + DD(4,n)*dpsi_n)
      else if (psi_st > 3.0D0 .and. psi_st < 4.0d0) then
        ! repulsive potential (ion-ion or electron-electron)
         F22   = 1.99016D0*e_psi_st - 4.56958D0
      else
        ! repulsive and attractive coeffcients are the same in this range
         F22   = 1.99016D0*e_psi_st - 4.56958D0
      end if
      OMEGA2(i) = eps_st*F22
   end do    ! i

   do Nref = 1, NN
      gamma2 = (4.0D0*kT/(ECHAR**2*KZ(Nref)))**2 * lambda2
      do i = 1, NN
         gamma_st2 = gamma2/KZ(i)**2

        ! Reduced mass of the two particles:
         A = (CAN1(i)*CAN1(Nref)/(CAN1(i)+CAN1(Nref)))        ! in units of AMU

        ! Units of Omega: length^2 * velocity:
         eps_st = CPI4*lambda2/gamma_st2 * SQRT(2.0D0*kT/(CPI4 * AMU*A))

        ! Find interpolation interval for dimensionless collision integral:
         e_psi_st = dlog(1.0D0+gamma_st2)
         psi_st = dlog(e_psi_st)

        ! Evaluate the collision integral:
         if (psi_st <= 3.0D0) then

           ! If psi_st falls outside range of Paquette fit, then just take border value, will probably never happen:
            psi_st = MAX(-7.D0, psi_st)

           ! Use spline interpolation to evaluate the collision integrals
           ! Determine interval in the table
            n = MIN(50, 1+floor((psi_st +7d0)/0.2d0))
            dpsi_n1 = (-7.2D0 + 0.2D0*(n+1)) - psi_st
            dpsi_n = psi_st - (-7.2D0 + 0.2D0*n)
            if ((Nref==NN.and.Nref/=i) .or. (i==NN.and.Nref/=i)) then
              ! attractive potential (electron-ion)
               do j=1,3
                  F1(j) = dexp(DCAT(1,n,j)*dpsi_n1**3 + DCAT(2,n,j)*dpsi_n**3  &
                      + DCAT(3,n,j)*dpsi_n1    + DCAT(4,n,j)*dpsi_n)
               end do
               F22 = dexp(DDAT(1,n)*dpsi_n1**3 + DDAT(2,n)*dpsi_n**3  &
                   + DDAT(3,n)*dpsi_n1    + DDAT(4,n)*dpsi_n)
              ! repulsive potential (ion-ion or electron-electron)
            else
               do j=1,3
                  F1(j) = dexp(DC(1,n,j)*dpsi_n1**3 + DC(2,n,j)*dpsi_n**3  &
                      + DC(3,n,j)*dpsi_n1    + DC(4,n,j)*dpsi_n)
               end do
               F22 = dexp(DD(1,n)*dpsi_n1**3 + DD(2,n)*dpsi_n**3  &
                   + DD(3,n)*dpsi_n1    + DD(4,n)*dpsi_n)
            end if
         else if (psi_st > 3.0D0.and.psi_st < 4.0d0) then
            if ((Nref==NN.and.Nref/=i) .or. (i==NN.and.Nref/=i)) then
              ! attractive potential (electron-ion)
               F1(1) = 1.01101D0*e_psi_st - 3.19815D0
               F1(2) = 1.04230D0*e_psi_st - 1.89637D0
               F1(3) = 2.15672D0*e_psi_st - 2.81038D0
               F22   = 2.08699D0*e_psi_st - 5.81444D0
            else
              ! repulsive potential (ion-ion or electron-electron)
               F1(1) = 1.00141D0*e_psi_st - 3.18209D0
               F1(2) = 0.99559D0*e_psi_st - 1.29553D0
               F1(3) = 1.99814D0*e_psi_st - 0.64413D0
               F22   = 1.99016D0*e_psi_st - 4.56958D0
            end if
         else if (psi_st.ge.4.0D0) then
           ! repulsive and attractive coeffcients are the same in this range
            F1(1) = 1.00141D0*e_psi_st - 3.18209D0
            F1(2) = 0.99559D0*e_psi_st - 1.29553D0
            F1(3) = 1.99814D0*e_psi_st - 0.64413D0
            F22   = 1.99016D0*e_psi_st - 4.56958D0
         end if
         OMEGA1(:) = eps_st*F1(:)
         OMEGA22   = eps_st*F22                ! for particle species Nref & k

        ! Binary diffusion coefficient for Chapman and Cowling formalism (first approximation)
        !Ddiff(Nref, i) = 3.0D0*kT/(16.0D0*(NA1(Nref)+NA1(i))*(A*AMU)*OMEGA1(1))

         AA = OMEGA22/(5.D0*OMEGA1(1))
         BB = (5.D0*OMEGA1(2)-OMEGA1(3))/(5.D0*OMEGA1(1))
         CC = 2.D0*OMEGA1(2)/(5.D0*OMEGA1(1))-1.D0
         Xs = NA1(Nref)/(NA1(i)+NA1(Nref))                ! number concentration of reference species Nref
         Xt = NA1(i)/(NA1(i)+NA1(Nref))                   ! number concentration of species i
         Ms = CAN1(Nref)/(CAN1(Nref)+CAN1(i))
         Mt = CAN1(i)/(CAN1(Nref)+CAN1(i))
         Pst = 3.0D0*(Ms - Mt)**2 + 4.D0*Ms*Mt*AA
         EE = kT/(8.D0*Ms*Mt*OMEGA1(1))
         Ps = 8.D0*Ms*EE*OMEGA2(Nref)/(5.D0*kT)
         Pt = 8.D0*Mt*EE*OMEGA2(i)/(5.D0*kT)
         Ss = Ms*Ps-Mt*(3.D0*(Mt-Ms)+4.D0*Ms*AA)
         St = Mt*Pt-Ms*(3.D0*(Ms-Mt)+4.D0*Mt*AA)
         Qs = Ps*(6.D0*Mt*Mt + 5.D0*Ms*Ms - 4.D0*Ms*Ms*BB + 8.D0*Ms*Mt*AA)
         Qt = Pt*(6.D0*Ms*Ms + 5.D0*Mt*Mt - 4.D0*Mt*Mt*BB + 8.D0*Mt*Ms*AA)
         Qst = 3.D0*(Ms-Mt)**2*(5.D0-4.D0*BB) + 4.D0*Ms*Mt*AA*(11.D0-4.D0*BB) + 2.D0*Ps*Pt
        !DELTA = 5.0D0*CC*CC*(Ms*Ms*Ps*Xs*Xs + Mt*Mt*Pt*Xt*Xt + Pst*Xs*Xt)/(Xs*Xs*Qs + Xt*Xt*Qt + Xs*Xt*Qst)

        ! second approximation
        !Ddiff(Nref, i) = Ddiff(Nref, i)/(1.0D0 - DELTA)

        ! thermal coefficients
        !Ath(Nref,i) = 5.D0*CC*(Xs*Ss-Xt*St)/(Xs**2*Qs+Xt**2*Qt+Xs*Xt*Qst)
        ! resistance coefficients for Burgers equations
         Kdiff(Nref, i) = 16.D0/3.D0*NA1(Nref)*NA1(i)*(A*AMU)*OMEGA1(1)
         Zdiff(Nref, i) = -CC
         Zdiff1(Nref, i) = -2.D0*BB+2.5D0
         Zdiff2(Nref, i) = 5.D0*AA
      end do !i
   end do !k=Nref

end subroutine compute_full_diffusion_coeffs



! ------------------------------------------------------------------------------
! DIFFUSION_TERMS_BURGERS
!  Compute the diffusion coefficients by solving the Burgers equations,
!  using the subroutine of Thoul & al. (1994).
! ------------------------------------------------------------------------------
!  Input:
!   rho     - Density [g/cm^3]
!   T       - Temperature [k]
!   NN      - Number of particle species (not counting electrons)
!   n(i)    - Particle number densities [cm^-3]
!   Ne      - Electron number density [cm^-3]
!   Z(i)    - Average charge of particle species (partial ionisation)
!   A(i)    - Baryon number of particle species
!  Output (in cgs units):
!   AP(i)   - Pressure terms in the expression for the diffusion velocity
!   AT(i)   - Temperature terms in the expression for the diffusion velocity
!   AC(i,j) - composition terms in the expression for the diffusion velocity
! ------------------------------------------------------------------------------
subroutine diffusion_terms_burgers(rho,T,NN,n,Ne, Z,A, AP,AT,AC)
   use real_kind
   use constants

   implicit none
   integer, intent(in) :: NN
   real(double), intent(in) :: rho, T, n(NN), Ne, Z(NN), A(NN)
   real(double), intent(out) :: AP(NN), AT(NN), AC(NN, NN)
  !     Local variables
   real(double) :: C(NN+1)      ! Concentration
  !     Permutation of particle numbers, needed because the solver can only
  !     deal with species that have C(i) > 0
   integer :: Ndiff
   integer, allocatable :: NP(:)
   real(double), allocatable :: NX(:), AA(:), ZZ(:)
   real(double), allocatable :: AAP(:), AAT(:), AAC(:, :)
  !     Diffusion coefficients, for Burgers equations
   real(double), allocatable :: Kdiff(:, :), Zdiff(:, :)
   real(double), allocatable :: Zdiff1(:, :), Zdiff2(:, :)
  !     Miscellaneous
   real(double) :: u;
   integer :: i, j

   allocate(NP(NN+1))
   allocate(NX(NN+1))
   allocate(AA(NN+1))
   allocate(ZZ(NN+1))
   allocate(AAP(NN+1))
   allocate(AAT(NN+1))
   allocate(AAC(NN+1, NN+1))
   allocate(Kdiff(NN+1,NN+1))
   allocate(Zdiff(NN+1,NN+1))
   allocate(Zdiff1(NN+1,NN+1))
   allocate(Zdiff2(NN+1,NN+1))

  !     Initialise all return values to 0
   AP(:) = 0.0d0
   AT(:) = 0.0d0
   AC(:,:) = 0.0d0

  !     Eliminate elements that are not present (or have very small abundances)
   Ndiff = 0
   do i = 1, NN
      if (n(i) > 1.0d-10) then
         Ndiff = Ndiff+1
         NP(Ndiff) = i;
      end if
   end do

  !     We need to have at least two species to calculate diffusion
   if (Ndiff < 2) then
      deallocate(NP)
      deallocate(NX)
      deallocate(AA)
      deallocate(ZZ)
      deallocate(AAP)
      deallocate(AAT)
      deallocate(AAC)
      deallocate(Kdiff)
      deallocate(Zdiff)
      deallocate(Zdiff1)
      deallocate(Zdiff2)
      return
   end if

  !     We have one additional species: electrons
   Ndiff = Ndiff+1

  !     Copy data to temporary arrays of appropriate size
   do i = 1, Ndiff-1
      NX(i) = n(NP(i))
      AA(i) = A(NP(i))
      ZZ(i) = Z(NP(i))
   end do
  !     Add electrons
   NX(Ndiff) = Ne
   AA(Ndiff) = AME/AMU
   ZZ(Ndiff) = -1.0d0

  !     Concentrations, relative to electrons
   C(1:Ndiff) = NX(1:Ndiff) / Ne

  !     Compute atomic diffusion coefficients, Paquette et al. (1986)
   Kdiff(:,:) = 0.0d0
   Zdiff(:,:) = 0.0d0
   Zdiff1(:,:) = 0.0d0
   Zdiff2(:,:) = 0.0d0

  !call compute_full_diffusion_coeffs(rho, T, Ndiff, AA, ZZ, NX, Kdiff, Zdiff, Zdiff1, Zdiff2)  !Rho not used
   call compute_full_diffusion_coeffs(T, Ndiff, AA, ZZ, NX, Kdiff, Zdiff, Zdiff1, Zdiff2)

  !     Kdiff = kappa, equation (37) of Thoul et al. (1994)
   Kdiff(:,:) = Kdiff(:,:) * T**1.5 / (1.41d-25*NX(Ndiff)**2)

  !     Find pressure, temperature and composition terms for diffusion velocities
   call diffusion_screened(Ndiff, AA, ZZ, C, Kdiff, Zdiff, Zdiff1, Zdiff2,  &
       AAP(1:Ndiff), AAT(1:Ndiff), AAC(1:Ndiff,1:Ndiff))

  !     Fix units: solar radii, 1e7 k, 1e2 g/cm^2, 6e13 years
   u = (T/1.0D7)**2.5/(rho/1.0D2) * (CRSN*1.0D11)**2/(6.D13*CSY)
   AAP(1:Ndiff) = AAP(1:Ndiff) * u
   AAT(1:Ndiff) = AAT(1:Ndiff) * u
   AAC(1:Ndiff,1:Ndiff) = AAC(1:Ndiff,1:Ndiff) * u

  !     Unscramble the results, don't include electrons (we're not treating
  !     them separately anyway)
   do i = 1, Ndiff-1
      AP(NP(i)) = AAP(i)
      AT(NP(i)) = AAT(i)
      do j = 1, Ndiff-1
         AC(NP(i),NP(j)) = AAC(i, j)
      end do
   end do

   deallocate(NP)
   deallocate(NX)
   deallocate(AA)
   deallocate(ZZ)
   deallocate(AAP)
   deallocate(AAT)
   deallocate(AAC)
   deallocate(Kdiff)
   deallocate(Zdiff)
   deallocate(Zdiff1)
   deallocate(Zdiff2)
end subroutine diffusion_terms_burgers


!> ------------------------------------------------------------------------------
!! DIFFUSION_SCREENeD
!!  Calculate diffusion coefficients by solving Burgers' equations.
!!
!! This routine was originally written by Anne A. Thoul, at the Institute
!! for Advanced Study, Princeton, NJ 08540.
!! See Thoul et al., Ap.j. 421, p. 828 (1994)
!! The subroutines LU_BACK_SUBSTITUTION and LU_DECOMPOSITION are based on
!! Numerical Recipes.
!!
!! \todo FIXME: replace these with LAPACK functions, which are free. See the
!! equivalent subroutine in MESA to see how this should be done.
!!
!! It has been updated to solve the Burgers
!! equations with resistance coefficients from a screened Coulomb
!! potential (Paquette et al 1986) instead of a pure Coulomb potential,
!! by Hu&al 2009
!!
!! The system contains n equations with n unknowns.
!! The equations are: the M momentum equations,
!!                    the M energy equations,
!!                    two constraints: the current neutrality
!!                                     the zero fluid velocity.
!! The unknowns are: the M diffusion velocities,
!!                   the M heat fluxes,
!!                   the electric field E
!!                   the gravitational force g.
!! ------------------------------------------------------------------------------
!!  Input:
!!     M                    - Number of atomic species to consider (+electrons)
!!     A(1:M)               - Atomic mass numbers
!!     Z(1:M)               - Average charge for this species (from Saha equation)
!!     C(1:M)               - Concentrations (by number/electron)
!!     k(1:M, 1:M)          - Relative diffusion coefficients for species i against j
!!     Zdiff(1:M, 1:M)      - Relative diffusion coefficients for species i against j
!!     Zdiff1(1:M, 1:M)     - Relative diffusion coefficients for species i against j
!!     Zdiff2(1:M, 1:M)     - Relative diffusion coefficients for species i against j
!!     See Hu & al. eq (1) and (2) for the definition of these coefficients
!!  Output:
!!     AP(1:M)              - Coefficient for pressure term in diffusion equation
!!     AT(1:M)              - Coefficient for temperature term in diffusion equation
!!     AX(1:M,1:M)          - Coefficient for composition term in diffusion equation
!!     See (3)  of Thoul et al (1994). Take care of the units!
!! ------------------------------------------------------------------------------
!<
subroutine diffusion_screened(M, A, Z, C, k, Zdiff, Zdiff1, Zdiff2, AP, AT, AX)
  ! Note that kappa_st is obtained from the resistance coefficient Kdiff with eq (37) Thoul&al.
  ! The vector AP, AT, and array AX contains the results for the diffusion
  ! coefficients used in eq 3  of Thoul et al (1994). Take care of the units!
  ! Warning: this routine does not allow species with 0 mass fraction,
  ! so make sure you leave them out!

   use real_kind

   implicit none
   integer, intent(in) :: M
   real(double), intent(in) :: A(M), Z(M), C(M), k(M,M), Zdiff(M,M), Zdiff1(M,M), Zdiff2(M,M)
   real(double), intent(out) :: AP(M), AT(M), AX(M,M)
   integer :: n, i, j, L
   integer indx(2*m+2)
   real(double) :: CC, AC
   real(double), allocatable :: XX(:, :), Y(:, :), YY(:, :)
   real(double), allocatable :: ALPHA(:), NU(:), gamma(:, :), DELTA(:, :), GA(:)
   real(double) :: KO, D

  ! The vector C contains the concentrations
  ! CC is the total concentration: CC=sum(C_s)
  ! AC is proportional to the mass density: AC=sum(A_s C_s)
  ! The arrays XX, Y, YY and k are various parameters which appear in
  ! Burgers equations.
  ! The vectors and arrays ALPHA, NU, gamma, DELTA, and GA represent
  ! the "right- and left-hand-sides" of Burgers equations, and later
  ! the diffusion coefficients.

  ! Initialise parameters:
   KO = 2.D0
   n = 2*M + 2

   ! Allocate data structures
   allocate(XX(M, M))
   allocate(Y(M, M))
   allocate(YY(M, M))
   allocate(ALPHA(2*M+2))
   allocate(NU(2*M+2))
   allocate(gamma(2*M+2, 2*M+2))
   allocate(DELTA(2*M+2, 2*M+2))
   allocate(GA(2*M+2))

  ! Calculate CC and AC:
   cc = sum(C(:))
   ac = dot_product(A(:), C(:))

  ! Calculate the coefficients of the burgers equations
   do i = 1, m
      do j = 1, m
         xx(i,j) = a(j)/(a(i) + a(j))
         y(i,j) = a(i)/(a(i) + a(j))
         yy(i,j) = 3.0d0 * y(i,j) + zdiff1(i,j) * xx(i,j) * a(j)/a(i)
      end do
   end do

  ! Write the burgers equations and the two constraints as
  ! alpha_s dp + nu_s dT + sum_t(not 2 or M) gamma_st dC_t
  !                     = sum_t delta_st w_t
  ! Note that we do *not* eliminate the He gradient, as in the original
  ! subroutine.
   do i = 1, m
      alpha(i) = c(i)/cc
      nu(i) = 0.d0
      gamma(i,1:m) = 0.d0
      do j = 1, m-1
         gamma(i,j) = -c(j)/cc
         if (j == i) then
            gamma(i,j) = gamma(i,j) + 1.d0
         end if
         gamma(i,j) = gamma(i,j)*c(i)/cc
      end do

      gamma(i,m+1:n) = 0.d0
   end do

   alpha(m+1:n) = 0.d0
   do i = m+1, n-2
      nu(i) = 2.5d0 * c(i-m)/cc
   end do
   nu(n-1:n) = 0.d0
   gamma(m+1:n,1:n) = 0.d0
   delta(1:n,1:n)=0.d0

   do i = 1, m
      do j = 1, m
         if (j == i) then
            do l = 1, m
               if(l /= i) then
                  delta(i,j) = delta(i,j) - k(i,l)
               end if
            end do
         else
            delta(i,j) = k(i,j)
         end if
      end do

      do j = M+1, n-2
         if(j-M == i) then
            do L = 1, M
               if (L /= i) then
                  DELTA(i,j) = DELTA(i,j) + Zdiff(i,L)*XX(i,L)*k(i,L)
               end if
            end do
         else
            DELTA(i,j) = -Zdiff(i,j-M)*Y(i,j-M)*k(i,j-M)
         end if
      end do

      DELTA(i,n-1) = C(i)*Z(i)

      DELTA(i,n) = -C(i)*A(i)
   end do

   do i = M+1, n-2
      do j = 1, M
         if (j == i-M) then
            do L = 1, M
               if (L /= i-M) then        ! pure C potential
                  DELTA(i,j) = DELTA(i,j) + 2.5D0*Zdiff(i-M,L)*XX(i-M,L)*k(i-M,L)
               end if
            end do
         else                             ! pure C potential
            DELTA(i,j) = -2.5D0*Zdiff(i-M,j)*XX(i-M,j)*k(i-M,j)
         end if
      end do

      do j = M+1, n-2
         if (j-M == i-M) then
            do L = 1, M
               if (L /= i-M) then         ! pure C potential
                  DELTA(i,j) = DELTA(i,j) - Y(i-M,L)*k(i-M,L)*  &
                      (0.8D0*Zdiff2(i-M,L)*XX(i-M,L) + YY(i-M,L))
               end if
            end do                         ! pure C potential
            DELTA(i,j) = DELTA(i,j) - 0.4D0*Zdiff2(i-M,i-M)*k(i-M,i-M)
         else                             ! pure C potential
            DELTA(i,j) = (3.D0 + Zdiff1(i-M,j-M) - 0.8D0*Zdiff2(i-M,j-M))*k(i-M,j-M)*XX(i-M,j-M)*Y(i-M,j-M)
         end if
      end do

      DELTA(i,n-1) = 0.D0
      DELTA(i,n) = 0.D0
   end do

   DELTA(n-1,1:M) = C(1:M)*Z(1:M)
   DELTA(n-1,M+1:n) = 0.D0

   DELTA(n,1:m) = C(1:m)*DNINT(A(1:m)) ! Diffusion velocities are for baryon fractions and not mass fractions
   DELTA(n,M+1:n) = 0.D0

  ! Invert the system for each possible right-hand-side, i.e.,
  ! if alpha is the r.h.s., we obtain the coefficient A_p
  ! if nu    ---------------------------------------- A_T
  ! if gamma(i,j) ----------------------------------- A_Cj
  !
  ! If i=1, we obtain the hydrogen diffusion velocity
  ! If i=2, ------------- helium   ------------------
  ! If i=3,M-1, --------- heavy element -------------
  ! If i=M, ------------- electrons -----------------
  ! Is it necessary to keep the elements in this order? i think not (HH)
  ! For i=M,2M, we get the heat fluxes
  ! For i=n-1, we get the electric field
  ! For i=n, we get the gravitational force g

  ! DELTA is now LU decomposition of the original matrix
   call lu_decomposition      (DELTA, n, n, INDX, D)

  ! LU decomposition of DELTA is output to DELTA
   call lu_back_substitution  (DELTA, n, n, INDX, ALPHA)
   call lu_back_substitution  (DELTA, n, n, INDX, NU)
   do j = 1, n
      ga(1:n) = gamma(1:n,j)
      call lu_back_substitution(DELTA, n, n, INDX, GA)
      gamma(1:n,j) = ga(1:n)
   end do

  ! The results for the coefficients must be multiplied by p/K_0:
   AP(1:M) = ALPHA(1:M)*KO*AC*CC
   AT(1:M) = NU(1:M)*KO*AC*CC
   AX(1:M,1:M) = gamma(1:M,1:M)*KO*AC*CC

   deallocate(XX)
   deallocate(Y)
   deallocate(YY)
   deallocate(ALPHA)
   deallocate(NU)
   deallocate(gamma)
   deallocate(DELTA)
   deallocate(GA)

end subroutine diffusion_screened



subroutine lu_back_substitution(A,n,NP,INDX,B)
   use real_kind
   implicit none

  !     .. Scalar Arguments ..
   integer :: n,NP
  !     ..
  !     .. Array Arguments ..
   real(double) :: A(NP,NP),B(n)
   integer :: INDX(n)
  !     ..
  !     .. Local Scalars ..
   real(double) :: SUM
   integer i,II,j,LL
  !     ..
   II = 0

   do i = 1,n
      LL = INDX(i)
      SUM = B(LL)
      B(LL) = B(i)
      if (II /= 0) then
         do j = II,i - 1
            SUM = SUM - A(i,j)*B(j)
         end do

      else if (SUM /= 0.D0) then
         II = i
      end if

      B(i) = SUM
   end do

   do i = n,1,-1
      SUM = B(i)
      if (i < n) then
         do j = i + 1,n
            SUM = SUM - A(i,j)*B(j)
         end do
      end if

      B(i) = SUM/A(i,i)
   end do
end subroutine lu_back_substitution



subroutine lu_decomposition(A,n,NP,INDX,D)
   use real_kind
   implicit none

  !     .. parameters ..
   real(double), parameter :: TINY=1.0D-20
  !     ..
  !     .. Scalar Arguments ..
   real(double) :: D
   integer :: n,NP
  !     ..
  !     .. Array Arguments ..
   real(double) :: A(NP,NP)
   integer :: INDX(n)
  !     ..
  !     .. Local Scalars ..
   real(double) :: AAMAX,DUM,SUM
   integer :: i,IMAX,j,k
  !     ..
  !     .. Local Arrays ..
   real(double), allocatable :: VV(:)
  !     ..
  !     .. Intrinsic Functions ..
   intrinsic dabs
  !     ..

   allocate(VV(n))

   imax = 1
   D = 1.D0
   do i = 1,n
      AAMAX = 0.D0
      do j = 1,n
         if (dabs(A(i,j)) > AAMAX) AAMAX = dabs(A(i,j))
      end do
     !if (AAMAX == 0.D0) PAUSE 'Singular matrix.'
      if (AAMAX == 0.D0) write(0,'(A)') 'Singular matrix in lu_decomposition()'  !> \todo Is a reason to stop the code?
      VV(i) = 1./AAMAX
   end do

   do j = 1,n
      if (j > 1) then
         do i = 1,j - 1
            SUM = A(i,j)
            if (i > 1) then
               do k = 1,i - 1
                  SUM = SUM - A(i,k)*A(k,j)
               end do
               A(i,j) = SUM
            end if

         end do
      end if


      AAMAX = 0.D0
      do i = j,n
         SUM = A(i,j)
         if (j > 1) then
            do k = 1,j - 1
               SUM = SUM - A(i,k)*A(k,j)
            end do
            A(i,j) = SUM
         end if

         DUM = VV(i)*dabs(SUM)
         if (DUM.GE.AAMAX) then
            IMAX = i
            AAMAX = DUM
         end if

      end do

      if (j /= IMAX) then
         do k = 1,n
            DUM = A(IMAX,k)
            A(IMAX,k) = A(j,k)
            A(j,k) = DUM
         end do
         D = -D
         VV(IMAX) = VV(j)
      end if

      INDX(j) = IMAX
      if (j /= n) then
         if (A(j,j) == 0.D0) A(j,j) = TINY
         DUM = 1./A(j,j)
         do i = j + 1,n
            A(i,j) = A(i,j)*DUM
         end do
      end if

   end do

   if (A(n,n) == 0.D0) A(n,n) = TINY

   deallocate(vv)

end subroutine lu_decomposition

