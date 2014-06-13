module indices
   implicit none

   ! Symbolic constants for indices in the array H
   ! Stellar structure
   integer, parameter :: VAR_LNF    = 1         ! ln f, the degeneracy parameter
   integer, parameter :: VAR_LNT    = 2         ! ln T, the temperature
   integer, parameter :: VAR_O16    = 3         ! O16 abundance, by baryon fraction
   integer, parameter :: VAR_MASS   = 4         ! M, Mass coordinate
   integer, parameter :: VAR_H1     = 5         ! H1 abundance, by baryon fraction
   integer, parameter :: VAR_QK     = 6         ! dQ/dk, gradient of mesh spacing function (eigenvalue)
   integer, parameter :: VAR_LNR    = 7         ! ln r, radius coordinate (actually ln sqrt(r**2 + ct(8)) )
   integer, parameter :: VAR_LUM    = 8         ! luminosity
   integer, parameter :: VAR_HE4    = 9         ! He4 abundance, by baryon fraction
   integer, parameter :: VAR_C12    = 10        ! C12 abundance, by baryon fraction
   integer, parameter :: VAR_NE20   = 11        ! Ne20 abundanc, by baryon fraction
   integer, parameter :: VAR_INERT  = 12        ! Integrated moment of inertia
   integer, parameter :: VAR_OMEGA  = 13        ! Rotation rate
   integer, parameter :: VAR_PHI    = 14        ! Value of potential
   integer, parameter :: VAR_PHIS   = 15        ! Potential difference with L1 (eigenvalue)
   integer, parameter :: VAR_N14    = 16        ! N14 abundance, by baryon fraction

   integer, parameter :: VAR_MG24   = 17        ! Mg24 abundance, by baryon fraction
   integer, parameter :: VAR_SI28   = 18        ! Si28 abundance, by baryon fraction
   integer, parameter :: VAR_FE56   = 19        ! Fe56 abundance, by baryon fraction
   integer, parameter :: VAR_TAM    = 20        ! Integrated angular momentum

   ! Total number of variables per star
   integer, parameter :: NVSTAR = 20
   integer, parameter :: INDEX_ORBIT_VAR_START = NVSTAR

   ! Binary parameters
   integer, parameter :: VAR_HORB   = INDEX_ORBIT_VAR_START+1
   integer, parameter :: VAR_ECC    = INDEX_ORBIT_VAR_START+2
   integer, parameter :: VAR_XI     = INDEX_ORBIT_VAR_START+3
   integer, parameter :: VAR_BMASS  = INDEX_ORBIT_VAR_START+4
   integer, parameter :: VAR_PMASS  = INDEX_ORBIT_VAR_START+5

   integer, parameter :: NVBIN = 5
   integer, parameter :: INDEX_SECONDARY_START = NVSTAR + NVBIN

   ! Total number of variables
   ! This needs to be at least 50 to accomodate the nucleosynthesis data
   integer, parameter :: nvar = max(50, nvstar + nvbin + nvstar)

   ! Total number of equations
   integer, parameter :: nestar = nvstar
   integer, parameter :: NEBIN = NVBIN
   integer, parameter :: INDEX_ORBIT_EQN_START = NESTAR
   integer, parameter :: neq = nvar


   ! Symbolic constants for equation numbers
   ! Stellar structure
   integer, parameter :: EQN_H1     = 1      ! Reaction-diffusion equation for H
   integer, parameter :: EQN_HE4    = 2      ! Reaction-diffusion equation for He
   integer, parameter :: EQN_C12    = 3      ! Reaction-diffusion equation for C
   integer, parameter :: EQN_N14    = 4      ! Reaction-diffusion equation for N
   integer, parameter :: EQN_O16    = 5      ! Reaction-diffusion equation for O
   integer, parameter :: EQN_NE20   = 6      ! Reaction-diffusion equation for Ne
   integer, parameter :: EQN_MG24   = 7      ! Reaction-diffusion equation for Mg
   integer, parameter :: EQN_SI28   = 8      ! Reaction-diffusion equation for Si
   integer, parameter :: EQN_FE56   = 9      ! Reaction-diffusion equation for Fe
   integer, parameter :: EQN_SUMX   = 10     ! Total variation in all elements (should be 0)

   integer, parameter :: EQN_PRES   = 11     ! Pressure equation; hydrostatic equilibrium
   integer, parameter :: EQN_RADIUS = 12     ! Radius equation; conservation of mass
   integer, parameter :: EQN_TEMP   = 13     ! Temperature equation; energy transport
   integer, parameter :: EQN_LUM    = 14     ! Luminosity equation; energy generation
   integer, parameter :: EQN_MASS   = 15     ! Mass equation; mesh spacing
   integer, parameter :: EQN_INERT  = 16     ! Integral for moment of inertia
   integer, parameter :: EQN_PHI    = 17     ! Potential (gravitation + rotation)
   integer, parameter :: EQN_TAM    = 18     ! Integral for total angular momentum
   integer, parameter :: EQN_OMEGA  = 19     ! Angular momentum transport equation

   ! Binary parameters
   integer, parameter :: EQN_OAM    = INDEX_ORBIT_EQN_START + 1  ! Orbital angular momentum
   integer, parameter :: EQN_ECC    = INDEX_ORBIT_EQN_START + 2  ! Orbital eccentricity
   integer, parameter :: EQN_XI     = INDEX_ORBIT_EQN_START + 3  ! Mass flow in contact
   integer, parameter :: EQN_BMASS  = INDEX_ORBIT_EQN_START + 4  ! Binary mass
   integer, parameter :: EQN_PMASS  = INDEX_ORBIT_EQN_START + 5  ! Primary mass

   ! Symbolic constants for boundary conditions
   ! We reuse the equation numbers from the structure part, but the correspondence
   ! is arbitrary except that we don't reuse the composition equations.
   ! Distinguish between surface (SBC) and central (CBC) boundary conditions
   ! Stellar structure
   integer, parameter :: SBC_MDOT   = EQN_PRES     ! Mass loss
   integer, parameter :: SBC_PRES   = EQN_RADIUS   ! Pressure
   integer, parameter :: SBC_TEMP   = EQN_TEMP     ! Temperature
   integer, parameter :: SBC_PHI    = EQN_LUM      ! Potential
   integer, parameter :: SBC_SSPIN  = EQN_MASS     ! Spin, solid body
   integer, parameter :: SBC_PHIS   = EQN_INERT    ! Surface potential (eigenvalue)

   integer, parameter :: CBC_MASS   = EQN_PRES
   integer, parameter :: CBC_LUM    = EQN_RADIUS
   integer, parameter :: CBC_RADIUS = EQN_TEMP
   integer, parameter :: CBC_INERT  = EQN_LUM
   integer, parameter :: CBC_TAM    = EQN_TAM

   ! Binary parameters
   integer, parameter :: SBC_OAM    = EQN_OAM
   integer, parameter :: SBC_ECC    = EQN_ECC
   integer, parameter :: CBC_XI     = EQN_XI
   integer, parameter :: SBC_BMASS  = EQN_BMASS
   integer, parameter :: SBC_PMASS  = EQN_PMASS



   ! Order in which variables are stored in the TWIN input files
   integer, parameter :: var_input_order(24) = (/ VAR_LNF, VAR_LNT, VAR_O16,&
               VAR_MASS, VAR_H1, VAR_QK, VAR_LNR, VAR_LUM, VAR_HE4, VAR_C12,&
               VAR_NE20, VAR_INERT, VAR_OMEGA, VAR_PHI, VAR_PHIS, VAR_N14, &
               VAR_HORB, VAR_ECC, VAR_XI, VAR_BMASS, 21, 22, 23, 24 /)

   ! Legacy order of equations, single star and orbital elements
   integer, parameter :: eqn_list_order(24) = (/ EQN_H1, EQN_O16, EQN_HE4, &
               EQN_C12, EQN_NE20, EQN_PRES, EQN_RADIUS, EQN_TEMP, EQN_LUM, &
               EQN_MASS, EQN_INERT, EQN_PHI, EQN_N14, EQN_MG24, EQN_SUMX,  &
               16, 17, 18, EQN_XI, 20, 21, 22, 23, 24 /)

   ! Legacy order of central boundary conditions
   integer, parameter :: cbc_list_order(24) = (/ EQN_H1, EQN_O16, EQN_HE4, &
               EQN_C12, EQN_NE20, CBC_MASS, CBC_LUM, CBC_RADIUS, CBC_INERT,&
               10, 11, 12, EQN_N14, EQN_MG24, EQN_SUMX, 16, 17, 18, CBC_XI, 20, 21, 22, 23, 24 /)

   ! Legacy order of surface boundary conditions
   integer, parameter :: sbc_list_order(24) = (/ EQN_H1, EQN_O16, EQN_HE4, &
               EQN_C12, EQN_NE20, SBC_MDOT, SBC_PRES, SBC_TEMP, SBC_PHI,   &
               SBC_SSPIN, SBC_PHIS, 12, EQN_N14, EQN_MG24, EQN_SUMX, 16,   &
               SBC_OAM, SBC_ECC, 19, SBC_BMASS, 21, 22, 23, 24 /)
               


   ! Number of dependent variables (`functions')
   ! These follow the same pattern of *1/bin/2 as the independent variables

   ! Stellar structure
   integer, parameter :: FN_BCP = 1
   integer, parameter :: FN_BCT = 2
   integer, parameter :: FN_VP = 3
   integer, parameter :: FN_VPK = 4
   integer, parameter :: FN_VR = 5
   integer, parameter :: FN_VRK = 6
   integer, parameter :: FN_VT = 7
   integer, parameter :: FN_VTK = 8
   integer, parameter :: FN_VL = 9
   integer, parameter :: FN_LK = 10
   integer, parameter :: FN_LQ = 11
   integer, parameter :: FN_MT = 12
   integer, parameter :: FN_VM = 13
   integer, parameter :: FN_VMK = 14
   integer, parameter :: FN_SG = 15
   integer, parameter :: FN_X1 = 16
   integer, parameter :: FN_X1T = 17
   integer, parameter :: FN_X16 = 18
   integer, parameter :: FN_X16T = 19
   integer, parameter :: FN_X4 = 20
   integer, parameter :: FN_X4T = 21
   integer, parameter :: FN_X12 = 22
   integer, parameter :: FN_X12T = 23
   integer, parameter :: FN_X20 = 24
   integer, parameter :: FN_X20T = 25
   integer, parameter :: FN_BCM = 26
   integer, parameter :: FN_VI = 27
   integer, parameter :: FN_VIK = 28
   integer, parameter :: FN_PHI = 29
   integer, parameter :: FN_PHIK = 30
   integer, parameter :: FN_BCF = 31
   integer, parameter :: FN_BCS = 32
   integer, parameter :: FN_BCPH = 33
   integer, parameter :: FN_X14 = 34
   integer, parameter :: FN_X14T = 35
   integer, parameter :: FN_AVMU = 36
   integer, parameter :: FN_SGTH = 37
   integer, parameter :: FN_OMEGA = 38
   integer, parameter :: FN_OMEGAT = 39
   integer, parameter :: FN_SAM = 40
   integer, parameter :: FN_SAMT = 41
   integer, parameter :: FN_SI = 42

   ! Net accretion mass flux
   integer, parameter :: FN_MACC = 43  ! Net mass accretion rate, 1e33g/s

   ! Scaled Roche potential
   integer, parameter :: FN_PHI_ROCHE = 44   ! In units of -G(M1+M2)/a

   ! Enthalpy
   integer, parameter :: FN_ENT = 45         ! In erg

   ! Total angular momentum, similar to the moment of inertia
   integer, parameter :: FN_AM  = 46  ! Total angular momentum
   integer, parameter :: FN_AMK = 47  ! Total angular momentum gradient

   ! Angular momentum transport
   integer, parameter :: FN_ADAM = 48 ! Angular momentum advection flux
   integer, parameter :: FN_SGAM = 49 ! Diffusion coefficient for angular momentum

   ! Meridional circulation
   integer, parameter :: FN_FU2K = 50 ! Vertical component of the circulation x rho r**2
   integer, parameter :: FN_FV   = 51 ! Prefactor for the continuity equation, 4pi r (dk/dm)/6
   integer, parameter :: FN_Vmc2 = 52 ! Horizontal component of meridional circulation
   integer, parameter :: FN_Umc2 = 53 ! Vertical component of meridional circulation
   integer, parameter :: FN_rho  = 54 ! Density
   integer, parameter :: FN_P    = 55 ! Pressure
   integer, parameter :: FN_HP   = 56 ! Pressure scale height

   ! Gravitational settling
   ! To do this properly, we need 9 fluxes for the 9 main isotopes and
   ! 9x9 = 81 relative diffusion coefficients. For the nucleosynthesis
   ! this is even worse; we'd need a whopping 1640 diffusion terms....
   ! The only way around this is to export the composition gradients to
   ! FUNCS1/FUNCS2, maybe through the semi-implicit functions.
   integer, parameter :: FN_FH  = 65  ! Advective flux hydrogen
   integer, parameter :: FN_FHe = 66  ! Advective flux helium
   integer, parameter :: FN_FC  = 67  ! Advective flux carbon
   integer, parameter :: FN_FN  = 68  ! Advective flux nitrogen
   integer, parameter :: FN_FO  = 69  ! Advective flux ocygen
   integer, parameter :: FN_FNe = 70  ! Advective flux neon
   integer, parameter :: FN_FMg = 71  ! Advective flux magnesium
   integer, parameter :: FN_FSi = 72  ! Advective flux silicon
   integer, parameter :: FN_FFe = 73  ! Advective flux iron

   ! Additional isotopes (Mg, Si, Fe)
   integer, parameter :: FN_X24  = 74  ! Magnesium abundance
   integer, parameter :: FN_X24t = 75  ! Magnesium abundance derivative
   integer, parameter :: FN_X28  = 76  ! Magnesium abundance
   integer, parameter :: FN_X28t = 77  ! Magnesium abundance derivative
   integer, parameter :: FN_X56  = 78  ! Magnesium abundance
   integer, parameter :: FN_X56t = 79  ! Magnesium abundance derivative


   integer, parameter :: nfstar = 80
   integer, parameter :: index_orbit_fn_start = nfstar

   ! Binary orbit
   integer, parameter :: FN_BCA = INDEX_ORBIT_FN_START + 1
   integer, parameter :: FN_BCE = INDEX_ORBIT_FN_START + 2
   integer, parameter :: FN_XIM = INDEX_ORBIT_FN_START + 3
   integer, parameter :: FN_XIK = INDEX_ORBIT_FN_START + 4
   integer, parameter :: FN_DLRK = INDEX_ORBIT_FN_START + 5
   integer, parameter :: FN_BCMB = INDEX_ORBIT_FN_START + 6
   integer, parameter :: FN_PMASS = INDEX_ORBIT_FN_START + 7
   integer, parameter :: FN_PHIL1 = INDEX_ORBIT_FN_START + 8
   integer, parameter :: FN_PHIL2 = INDEX_ORBIT_FN_START + 9

   integer, parameter :: nfbin = 16

   integer, parameter :: nfunc = nfstar + nfbin + nfstar
   integer, parameter :: index_secondary_fn_start = nfstar + nfbin


   contains


   
   ! Return the index corresponding to variable i for star 1
   elemental function idx_primary(i)
      implicit none
      integer, intent(in) :: i
      integer :: idx_primary
      integer :: ii

      ii = i
      if (ii <= nvar .and. ii > nvstar+nvbin) ii = ii - (nvstar+nvbin)

      idx_primary = ii
   end function idx_primary



   elemental function is_abundance(i)
      implicit none
      integer, intent(in) :: i
      logical :: is_abundance
      integer :: ii

      ii = idx_primary(i);

      is_abundance = ii == VAR_H1   .or. &
                     ii == VAR_He4  .or. &
                     ii == VAR_C12  .or. &
                     ii == VAR_N14  .or. &
                     ii == VAR_O16  .or. &
                     ii == VAR_NE20 .or. &
                     ii == VAR_MG24 .or. &
                     ii == VAR_SI28 .or. &
                     ii == VAR_FE56
   end function is_abundance



   pure function is_binary_var(i)
      implicit none
      integer, intent(in) :: i
      logical :: is_binary_var

      is_binary_var = i == VAR_HORB .or. &
                      i == VAR_ECC  .or. &
                      i == VAR_XI   .or. &
                      i == VAR_BMASS
   end function is_binary_var



   elemental function idx_for_star(vari, Jstar)
      implicit none
      integer, intent(in) :: vari, Jstar
      integer :: idx_for_star
      integer :: n, i

      i = idx_primary(vari)
      n = nvstar + nvbin
      if (i > nvstar .and. i <= nvstar + nvbin) n = 0 ! Orbital elements, same for both stars

      idx_for_star = n * (Jstar-1) + i
   end function idx_for_star



   elemental function fn_idx_for_star(fni, Jstar)
      implicit none
      integer, intent(in) :: fni, Jstar
      integer :: fn_idx_for_star
      integer :: n, i

      i = fni
      n = index_secondary_fn_start
      if (i > index_orbit_fn_start .and. i <= index_secondary_fn_start) n = 0 ! Orbital elements, same for both stars

      fn_idx_for_star = n * (Jstar - 1) + i
   end function fn_idx_for_star

end module indices
