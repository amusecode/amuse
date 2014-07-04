! This module defines all run parameters and options
module settings
   use real_kind
   use mesh, only: npx
   implicit none

   ! Flag to detect whether Mg24 is actually being solved for
   logical :: use_mg24_eqn = .false.

   real(double) :: CZS
   real(double) :: CALPHA_OVER_FE = 0.0d0
   real(double) :: CH = -1.0
   real(double) :: CDC(11)
   real(double) :: CT1, CT2, CT3
   real(double) :: INITIAL_CT(20)
   real(double) :: CT(20)
   real(double) :: CC, CN, CO, CNE, CMG, CSI, CFE, CNUC

   real(double) :: CALP, CU, COS, CPS, CRD, CTH, CTHE
   real(double) :: CXB, CGR
   real(double) :: CEA, CET
   real(double) :: CMT, CMS, CMI, CML, CHL, CTF, CLT
   real(double) :: CPA, CBR, CSU, CSD, CDF, CGW, CSO, CMB

   real(double) :: mtr_limit = 1.0d1  ! Limit mass transfer rate [Msun/yr]

   real(double) :: CLIMIT = 1.0D-1    ! Limit changes in variables during iterations

   ! Desired accuracy. The solver will aim for an accuracy between
   !   EPS and WANTED_EPS
   ! No effect if WANTED_EPS >= EPS
   real(double) :: eps, del, dh0
   real(double) :: wanted_eps = 1.0D-8
   real(double) :: max_eps = 1.0D-4
   real(double) :: max_erk = 0.0d0

   real(double) :: fac_down, fac_up, fac_min
   integer :: kr1, kr2

   real(double) :: accept_res

   real(double) :: CPHOTONTIRE = 0.0  ! Switch to include photon tiring

   ! Fully mix inner and outer gridpoints.
   ! Set either to 0, or so that the inner 20% and out 7% of gridpoints are
   ! mixed.
   integer :: inner_points_mixed = 0
   integer :: outer_points_mixed = 0

   ! Should the equation-of-state include the effects of pair production?
   ! This is only important in the very late burning stages of very massive
   ! stars. Positrons are only calculated if their degeneracy parameter
   ! >= -15.0 - otherwise they are negligible anyway.
   logical :: EOS_INCLUDE_PAIRPRODUCTION = .false.

   ! Turn on smart mass loss routine, that picks an appropriate
   ! recipe depending on the stellar parameters. This is an
   ! alternative for the de Jager rate and replaces it when
   ! SMART_MASS_LOSS is switched on. Off by default for reasons of
   ! backward compatibility, but should be enabled in the input file.
   real(double) :: SMART_MASS_LOSS = 0.0

   ! Individual mass loss recipe switches.
   ! These also turn on recipes when SMART_MASS_LOSS is used,
   ! although that does store its own set of mass loss options (to
   ! keep it more modular).
   real(double) :: CMR       ! Switch for Reimers-like mass loss rate
   real(double) :: CMJ       ! Switch for de Jager mass loss rate
   real(double) :: CMV       ! Switch for Vink mass loss rate
   real(double) :: CMK       ! Switch for Kudritzki 2002 mass loss rate
   real(double) :: CMNL      ! Switch for Nugis&Lamers mass loss rate (WR stars)
   real(double) :: CMRR      ! Switch for Real Reimers mass loss rate
   real(double) :: CMVW      ! Switch for Vasiliadis&Wood (AGB) rate
   real(double) :: CMSC      ! Switch for Schroeder&Cuntz mass loss
   real(double) :: CMW       ! Switch for Wachter&al (AGB) mass loss
   real(double) :: CMAL      ! Switch for Achmad&Lamers (A supergiants)

   ! Catch-all switch to enable/disable stellar winds in one go.
   ! Easier to maintain than setting all values explicitly.
   real(double) :: CMDOT_WIND = 1.0d0

   ! CV-like mass loss prescription, Ritter 1988
   ! This allows for mass transfer to begin before the star fully fills its
   ! Roche lobe.
   real(double) :: cmt_cv

   ! Rotationally enhanced mass loss rates, two options: Heger & al,
   ! Maeder & Meynet. Set one of these!
   real(double) :: CMDOTROT_HLW ! Heger, Langer & Woosely
   real(double) :: CMDOTROT_MM  ! Maeder & Meynet

   ! Non-conservative mass transfer options (depending on stellar parameters)
   real(double) :: CMTEL     ! Eddington-limited accretion (0 or 1)
   real(double) :: CMTWL     ! Angular momentum limited accretion

   ! Scaling with metallicity applied to de Jager mass loss rate in funcs1
   real(double) :: ZSCALING_MDOT = 0.8

   real(double) :: ARTMIX = 0.0D0   ! Artificial mixing coefficient [cm^2/s]

   real(double) :: CCAC = 0.0D0     ! Switch for composition accretion
   real(double) :: CLAC = 0.0D0     ! Switch for accretion luminosity

   real(double) :: CGRS = 0.0D0     ! Switch for gravitational settling
   logical :: grs_burgers = .false. ! Switch to solve Burger's equations for gravitational settling
   real(double) :: CRLEV = 0.0d0    ! Switch for radiative levitation
   character :: op_data_path*(128)  ! Directory holding Opacity Project data files
   integer :: rlev_update           ! Use current (0) or previous (-1) values for accelerations

   real(double) :: CSMC = 0.04D0! Semi-convection efficiency, after Langer (1991)
   real(double) :: CSMCE = 0.0D0! Semi-convection energy transport on/off

   real(double) :: Chs = 1.0d0      ! Switch for horizontal shear
   real(double) :: Cvs = 1.0d0      ! Switch for vertical shear
   real(double) :: Cshear = 1.0d0   ! Switch for shear mixing of chemicals
   real(double) :: Cadam = 1.0d0    ! Switch for advection term in angular momentum equation

   integer ::KTH, KX, KY, KZ
   integer ::KT1, KT2, KT3, KT4
   integer :: KCL, KION, KAM, KOP, KCC, KNUC, KCN
   integer :: KSX(45)
   integer :: KN, KJN(40)

   ! Number of allowed iterations for the nucleosynthesis code
   integer :: KR_NUCSYN = 60

   ! Variables derived from Settings and never changed
   real(double) :: CLOGZ
   logical :: rigid_rotation = .true.     ! Rigid rotation or differential rotation?

   ! Simple Roche model or detailed Clairaut treatment of distorted stars?
   logical :: use_clairaut_distortion = .false.

   ! Detailed mass transfer calculation for contact binaries?
   logical :: use_contact_flow = .false.

   ! Should output files include the changes DH or not?
   logical :: store_changes = .false.

   ! Switches for the new "smooth" remesher.
   ! The first two are set in init.dat, the last one in init.run
   logical :: use_smooth_remesher = .false.
   logical :: relax_loaded_model = .true.
   logical :: start_with_rigid_rotation = .true.

   ! Construct a post-He flash model from scratch or not
   logical :: construct_zahb_model = .false.

   ! Column output for .mdl files
   integer :: mdl_columns
   integer :: mdl_column_variables(NPX)


   ! Unused, but the code relies on this being defined:
   real(double) :: CQ1, CQ(17)
end module settings

