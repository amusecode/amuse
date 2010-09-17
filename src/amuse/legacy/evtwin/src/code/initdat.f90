!> \file initdat.f90  Routines to read the settings file init.dat, transfer its variables and set default values


!> \brief Routines to read init.dat, transfer its variables and set default values
module init_dat
  
   use real_kind
   
   implicit none
   
   type init_dat_settings
     
     real(double) :: eps, del, dh0, wanted_eps
     
     real(double) :: cdc_ms, cdc_hec, cdc_hes, cdc_dblsh, cdc5
     
     integer :: ke1, ke2, ke3, kbc, kev, kfn, kl, jh(3)
     integer :: kp_var(40), kp_eqn(40), kp_bc(40)
     
     integer :: kh2, kr1, kr2, jch, kth, kx, ky, kz
     integer :: kcl, kion, kam, kop, kcc, knuc, kcn
     integer :: kt1, kt2, kt3, kt4, kt5, ksv
     
     integer :: kn, kjn(40)
     
     real(double) :: ct1, ct2, ct3
     real(double) :: ct(20)
     
     real(double) :: ch, cc, cn, co, cne, cmg, csi, cfe
     real(double) :: calp, cu, cos, cps, crd, cth, cgrs, ccac
     real(double) :: cdsi, cshi, cssi, cesc, cgsf, cfmu, cfc
     real(double) :: artmix
     
     real(double) :: cxb, cgr, cea, cet
     real(double) :: cmt, cms, cml, chl, ctf, clt
     real(double) :: cmi, cmr, cmj, cmv, cmk, cmnl
     real(double) :: cpa, cbr, csu, csd, cdf, cgw, cso, cmb
     real(double) :: cphotontire
     
     integer :: convection_scheme
     logical :: use_fudge_control, allow_extension, allow_underrelaxation
     logical :: allow_overrelaxation, allow_mdotrelaxation
     logical :: allow_egenrelaxation, use_previous_mu
     logical :: allow_avmurelaxation, use_quadratic_predictions
     
     real(double) :: off_centre_weight
     
     integer :: ke1_2, ke2_2, ke3_2, kbc_2, kev_2, kfn_2, kl_2, jh_2(3)
     integer :: kp_var_2(40), kp_eqn_2(40), kp_bc_2(40)
     
  end type init_dat_settings
  
contains
  
  !> \brief Function for reading init.dat (fort.22) in either `classical'
  !! numbers-only format or more modern `with names' format      
  function read_init_dat(it, kh2, kr1, kr2, ksv, kt5, jch)
    use real_kind
    use zams_nucleosynthesis_abundances
    use nucleosynthesis
    use mesh
    use constants
    use settings
    use control
    use extrapolate_dh
    use massloss
    use accretion_abundances
    
    implicit none
    logical :: read_init_dat
    integer, intent(in) :: it
    integer, intent(out) :: kh2, kr1, kr2, ksv, kt5, jch
    integer :: ioerror
    real(double) :: ch_opac
    integer :: ke1, ke2, ke3, kbc, kev, kfn, kl, jh(3),kp_var(40), kp_eqn(40), kp_bc(40),   &
         ke1_2, ke2_2, ke3_2, kbc_2, kev_2, kfn_2, kl_2, jh_2(3),kp_var_2(40), kp_eqn_2(40), kp_bc_2(40)
    real(double) :: cdc_ms, cdc_hec, cdc_hes, cdc_dblsh, cdc5,   &
         cdc_ems, cdc_hg, cdc_1dup, cdc_rlof, cdc_rlof_reduce,  &
         unused1, unused(17)
    
    logical :: accret_composition
    
    ! Namelist for the new init.dat I/O
    namelist /init_dat/ kh2, kr1, kr2, jch, kth, kx, ky, kz,   &
         kcl, kion, kam, kop, kcc, knuc, kcn,  &
         kt1, kt2, kt3, kt4, kt5, ksv,  &
         eps, wanted_eps, del, dh0, cdc_ms, cdc_hec, cdc_hes, cdc_dblsh, cdc5,   &
         cdc_ems, cdc_hg, cdc_1dup, cdc_rlof, cdc_rlof_reduce,  &
         ke1, ke2, ke3, kbc, kev, kfn, kl, jh, kp_var, kp_eqn, kp_bc,   &
         ke1_2, ke2_2, ke3_2, kbc_2, kev_2, kfn_2, kl_2, jh_2, kp_var_2, kp_eqn_2, kp_bc_2,  &
         ksx, kn, kjn,   &
         ct1, ct2, ct3, ct,   &
         ch, cc, cn, co, cne, cmg, csi, cfe,   &
         calp, cu, cos, cps, crd, cxb, cgr, cea, cet,  &
         cmt, cms, cmi, cmr, cmj, cmv, cmk, cmnl, cmrr, cmvw, cmsc, cmw,  &
         cmal, smart_mass_loss, cml, chl, ctf, clt, cmdotrot_hlw, cmdotrot_mm,  &
         cmtel, cmtwl,  &
         cpa, cbr, csu, csd, cdf, cgw, cso, cmb, cphotontire, unused1,  &
         cth, cthe, unused, use_fudge_control, allow_extension, allow_underrelaxation,  &
         allow_overrelaxation, convection_scheme, allow_egenrelaxation,  &
         use_previous_mu, off_centre_weight, allow_mdotrelaxation,  &
         use_smooth_remesher, relax_loaded_model, convection_ledoux,  &
         allow_avmurelaxation, store_changes, use_quadratic_predictions, use_linesearch, climit,  &
         cmi_mode, zscaling_mdot, cdsi, cshi, cssi,  &
         cesc, cgsf, cfmu, cfc, artmix, cgrs, grs_burgers, ccac, csmc, csmce, &
         crlev, op_data_path, rlev_update, kr_nucsyn
    
    ! Namelist for the accretion of matter
    namelist /accret/ x1ac, x4ac, x12ac, x14ac, x16ac, x20ac, x24ac, x28ac, x56ac, accret_composition
    
    ! Namelist for initial (ZAMS) nucleosynthesis abundances
    namelist /abund/ cxd, cxhe3, cxli7, cxbe7, cxb11, cxc12, cxc13,  &
         cxc14, cxn14, cxn15, cxo16, cxo18, cxo17, cxf19, cxne21, cxne20,  &
         cxne22, cxna22, cxna23, cxmg24, cxmg25, cxmg26, cxal26m, cxal27,  &
         cxal26g, cxsi28, cxsi30, cxsi29, cxp31, cxs33, cxs32, cxs34, cxfe57,  &
         cxfe60, cxfe56, cxfe58, cxfe59, cxco59, cxni61, cxni59, cxni58,  &
         cxni60, cxCa40
    
    ! Default values
    ch_opac = ch
    store_changes = .false.   ! Don't store changes DH for later input
    cphotontire = 0.0d0  ! Do not keep track of photon tiring (do not take kinetic energy of the wind from the luminosity)
    cmv = 0.0d0    ! Vink mass loss, disabled by default
    cmk = 0.0d0    ! Kudritzki 2002 mass loss, disabled by default
    cmnl = 0.0d0   ! Nugis&Lamers, for WR stars, disabled by default
    ch = -1.0d0    ! Use default value from opacity table 
    artmix = 0.0d0 ! No artificial mixing
    cgrs  = 0.0d0  ! No gravitational settling/diffusion
    grs_burgers = .false.  ! Default to trace element (Cowling) formalism for settling
    ccac = 0.0d0   ! Default (single star): don't accrete composition
    cdc(:) = 1.0d0 ! Default timestep control
    cdc_ems = 1.0d0! Timestep parameter at end of MS
    cdc_hg = 1.0d0 ! Timestep parameter in HG
    cdc_1dup = 1.0d0! Timestep parameter at 1DUP
    cdc_rlof = 1.0d0! Timestep parameter if star is close to RLOF
    cdc_rlof_reduce = 1.0d0! Timestep if star is close to RLOF but contracting
    csmc = 0.04d0  ! Efficiency of semi-convection
    csmce = 0.0d0  ! Energy transport by semi convection off
    cthe = 0.0d0   ! Energy transport by thermohaline mixing off
    convection_ledoux = 0.0d0  ! Default: Schwarzschild convection
    cmdotrot_hlw = 0.0d0 ! Enhanced mass loss from Heger, Langer & Woosely
    cmdotrot_mm = 0.0d0  ! Enhanced mass loss from Maeder & Meynet
    cmtel = 0.0d0  ! No Eddington-limited accretion
    cmtwl = 0.0d0  ! No angular momentum limited accretion
    ct(:) = 0.0d0  ! Make sure extra MSF terms are forced to 0
    crlev = 0.0d0  ! No radiative levitation
    rlev_update = -1! By default, use previous values for radiative levitation
    op_data_path = 'mono/' ! Default path for OP data files (semi-sensible)
    use_linesearch = .false.
    accret_composition = .false.
    kr_nucsyn = 60
    
    ! Assigned to cdc(5) below, which is never used(?), but can cause unexpected behaviour if undefined.
    ! It's still read from the old init.dat and needs to be retained (and defined) for now.
    cdc5 = 1.d0
    
    ! Improved mass loss rates, can be used in a general mass loss recipe
    ! instead of the current de Jager rate. These individual mass loss
    ! recipes can be switched on and off there.
    smart_mass_loss = 0.0
    cmrr = multiplier_reimers
    cmvw = multiplier_vasiliadis_wood
    cmsc = multiplier_schroeder
    cmw = multiplier_wachter
    cmal = multiplier_achmad
    
    if (ktw == 2) ccac = 1.0d0 ! Default (binary/TWIN): do accrete composition
    x1ac = -1.0
    x4ac = -1.0
    x12ac = -1.0
    x14ac = -1.0
    x16ac = -1.0
    x20ac = -1.0
    x24ac = -1.0
    x28ac = -1.0
    x56ac = -1.0
    cmi_mode = 1
    
    ! First, try to read the NAMELIST formatted init.dat
    read (it, nml=init_dat, iostat=ioerror)
    kfn = nfunc    ! Just copy all functions; common error to forget. 
    kfn_2 = nfunc  ! Just copy all functions; common error to forget. 
    if (ke2_2 /= 0) nucleosynthesis_enabled = .true.
    if (ch < 0.0) ch = ch_opac    ! Use CH from opacity table
    if (use_previous_mu) avmu_smooth = 0.0 ! Use mu from prev. timestep
    if (wanted_eps > eps) wanted_eps = eps
    
    ! Turn on mass loss recipes in smart mass loss routine, pass on the
    ! exponent for metallicity scaling that is to be used.
    multiplier_dejager = cmj
    multiplier_schroeder = cmsc
    multiplier_reimers = cmrr
    multiplier_vasiliadis_wood = cmvw
    multiplier_wachter = cmw
    multiplier_achmad = cmal
    multiplier_vink = cmv
    multiplier_kudritzki = cmk
    multiplier_nl = cmnl
    metallicity_scaling_exponent = zscaling_mdot
    
    ! Store a copy of the mesh spacing function coefficients, we want to
    ! dynamically change them during some evolutionary stages and we need to
    ! be able to restore the original settings.
    initial_ct(:) = ct(:)
    if (ioerror == 0) then
       read_init_dat = .true.
       
       ! Allocate memory for storage of previous models in parabola extrapolation
       if (use_quadratic_predictions) call initlse_parabola_storage_space(kh2, ke1+ke2+kev)
       
       ! Export some local names to their global counter parts
       cdc(1) = cdc_ms 
       cdc(2) = cdc_hec
       cdc(3) = cdc_hes
       cdc(4) = cdc_dblsh
       cdc(5) = cdc5
       cdc(6) = cdc_ems
       cdc(7) = cdc_hg
       cdc(8) = cdc_1dup
       cdc(9) = cdc_rlof
       cdc(10) = cdc_rlof_reduce
       
       ! Set some leftover junk that's never actually used
       unused1 = 0.0
       unused(:) = 0.0d0
       cq1 = unused1
       cq(:) = unused(:)
       ! Try to get accretion information from the same file
       ! This is a bit ugly, but do we want to use another fort.nnn file for
       ! this?
       read (it, nml=accret, iostat=ioerror)
       
       ! Finally, read the nucleosynthesis abundances
       read (it, nml=abund, iostat=ioerror)

       ! Copy read in variables and permutations to global data structures
       id(1) = ke1
       id(2) = ke2
       id(3) = ke3
       id(4) = kbc
       id(5) = kev
       id(6) = kfn
       id(7) = kl
       id(8:10) = jh(1:3)
       id(11:50) = kp_var
       id(51:90) = kp_eqn
       id(91:130) = kp_bc

       ie(1) = ke1_2
       ie(2) = ke2_2
       ie(3) = ke3_2
       ie(4) = kbc_2
       ie(5) = kev_2
       ie(6) = kfn_2
       ie(7) = kl_2
       ie(8:10) = jh_2(1:3)
       ie(11:50) = kp_var_2
       ie(51:90) = kp_eqn_2
       ie(91:130) = kp_bc_2

       return
    end if
    
    ! Fall back to the `old' bare numbers init.dat; close file to seek from start
    rewind (it)
    read  (it, 993, iostat=ioerror) kh2, kr1, kr2, jch, kth, kx, ky, kz,  &
         kcl, kion, kam, kop, kcc, knuc, kcn, kt1, kt2, kt3, kt4, kt5, ksv,   &
         eps, del, dh0, cdc_ms, cdc_hec, cdc_hes, cdc_dblsh, cdc5,  &
         ke1, ke2, ke3, kbc, kev, kfn, kl, jh, kp_var, kp_eqn, kp_bc,   &
         ke1_2, ke2_2, ke3_2, kbc_2, kev_2, kfn_2, kl_2, jh_2, kp_var_2, kp_eqn_2, kp_bc_2,  &
         ksx, kn, kjn, ct1, ct2, ct3, ct(1:10),   &
         cc, cn, co, cne, cmg, csi, cfe,   &
         calp, cu, cos, cps, crd, cxb, cgr, cea, cet,   &
         cmt, cms, cmi, cmr, cmj, cml, chl, ctf, clt,   &
         cpa, cbr, csu, csd, cdf, cgw, cso, cmb, unused1,  &
         cth, unused
    ! FORMAT specifier for the old init.dat format
993 format (8I5, /, 7I5, /, 6I5, /, 1P, 8D8.1, 0P, /, 2(10I4, /,  &
         6(20I3, /)), 3(15I3, /), I3, /, 2(20I3, /), 10F5.2, 1P, 3D8.1, /,  &
         0P, 7F6.3, /, 1P, 5(9D9.2, /), 0P)
    
    kfn = nfunc    ! Just copy all functions; common error to forget. 
    
    initial_ct(:) = ct(:)
    ! Export some local names to their global counter parts
    cq1 = unused1
    cq(:) = unused(:)
    cdc(1) = cdc_ms 
    cdc(2) = cdc_hec
    cdc(3) = cdc_hes
    cdc(4) = cdc_dblsh
    cdc(5) = cdc5
    wanted_eps = eps
    
    if (ioerror == 0) then
       if (use_quadratic_predictions) call initlse_parabola_storage_space(kh2, ke1+ke2+kev)
       read_init_dat = .true.
    else
       ! Error during read operation
       read_init_dat = .false.
    end if
    
    ! Rewind the file, in case we want to read it again later
    rewind (it)

    ! Copy read in variables and permutations to global data structures
    id(1) = ke1
    id(2) = ke2
    id(3) = ke3
    id(4) = kbc
    id(5) = kev
    id(6) = kfn
    id(7) = kl
    id(8:10) = jh(1:3)
    id(11:50) = kp_var
    id(51:90) = kp_eqn
    id(91:130) = kp_bc

    ie(1) = ke1_2
    ie(2) = ke2_2
    ie(3) = ke3_2
    ie(4) = kbc_2
    ie(5) = kev_2
    ie(6) = kfn_2
    ie(7) = kl_2
    ie(8:10) = jh_2(1:3)
    ie(11:50) = kp_var_2
    ie(51:90) = kp_eqn_2
    ie(91:130) = kp_bc_2
  end function read_init_dat
  
  
  !> \brief Store variables in struct init_dat
  subroutine push_init_dat(init_dat, KH2, KR1, KR2, KSV, KT5, JCH)
    use real_kind
    use mesh
    use constants
    use settings
    use control
    use extrapolate_dh
    
    implicit none
    type(init_dat_settings), intent(out) :: init_dat
    integer, intent(in) :: kh2, kr1, kr2, ksv, kt5, jch
    
    init_dat%eps = eps
    init_dat%del = del
    init_dat%dh0 = dh0
    init_dat%wanted_eps = wanted_eps
    init_dat%cdc_ms = cdc(1)
    init_dat%cdc_hec = cdc(2)
    init_dat%cdc_hes = cdc(3)
    init_dat%cdc_dblsh = cdc(4)
    init_dat%ke1 = id(1)
    init_dat%ke2 = id(2)
    init_dat%ke3 = id(3)
    init_dat%kbc = id(4)
    init_dat%kev = id(5)
    init_dat%kfn = id(6)
    init_dat%kl = id(7)
    init_dat%jh(1:3) = id(8:10)
    init_dat%kp_var(1:40) = id(11:50)
    init_dat%kp_eqn(1:40) = id(51:90)
    init_dat%kp_bc(1:40) = id(91:130)
    init_dat%kh2 = kh2
    init_dat%kr1 = kr1
    init_dat%kr2 = kr2
    init_dat%jch = jch
    init_dat%kth = kth
    init_dat%kx = kx
    init_dat%ky = ky
    init_dat%kz = kz
    init_dat%kcl = kcl
    init_dat%kion = kion
    init_dat%kam = kam
    init_dat%kop = kop
    init_dat%kcc = kcc
    init_dat%knuc = knuc
    init_dat%kcn = kcn
    init_dat%kt1 = kt1
    init_dat%kt2 = kt2
    init_dat%kt3 = kt3
    init_dat%kt4 = kt4
    init_dat%kt5 = kt5
    init_dat%ksv = ksv
    init_dat%kn = kn
    init_dat%kjn(1:40) = kjn(1:40)
    init_dat%ct1 = ct1
    init_dat%ct2 = ct2
    init_dat%ct3 = ct3
    init_dat%ct(1:20) = ct(1:20)
    init_dat%ch = ch
    init_dat%cc = cc
    init_dat%cn = cn
    init_dat%co = co
    init_dat%cne = cne
    init_dat%cmg = cmg
    init_dat%csi = csi
    init_dat%cfe = cfe
    init_dat%calp = calp
    init_dat%cu = cu
    init_dat%cos = cos
    init_dat%cps = cps
    init_dat%crd = crd
    init_dat%cth = cth
    init_dat%cgrs = cgrs
    init_dat%ccac = ccac
    init_dat%cdsi = cdsi
    init_dat%cshi = cshi
    init_dat%cssi = cssi
    init_dat%cesc = cesc
    init_dat%cgsf = cgsf
    init_dat%cfmu = cfmu
    init_dat%cfc = cfc
    init_dat%artmix = artmix
    init_dat%cxb = cxb
    init_dat%cgr = cgr
    init_dat%cea = cea
    init_dat%cet = cet
    init_dat%cmt = cmt
    init_dat%cms = cms
    init_dat%cml = cml
    init_dat%chl = chl
    init_dat%ctf = ctf
    init_dat%clt = clt
    init_dat%cmi = cmi
    init_dat%cmr = cmr
    init_dat%cmj = cmj
    init_dat%cmv = cmv
    init_dat%cmk = cmk
    init_dat%cmnl = cmnl
    init_dat%cpa = cpa
    init_dat%cbr = cbr
    init_dat%csu = csu
    init_dat%csd = csd
    init_dat%cdf = cdf
    init_dat%cgw = cgw
    init_dat%cso = cso
    init_dat%cmb = cmb
    init_dat%cphotontire = cphotontire
    init_dat%convection_scheme = convection_scheme
    init_dat%use_fudge_control = use_fudge_control
    init_dat%allow_extension = allow_extension
    init_dat%allow_underrelaxation = allow_underrelaxation
    init_dat%allow_overrelaxation = allow_overrelaxation
    init_dat%allow_mdotrelaxation = allow_mdotrelaxation
    init_dat%allow_egenrelaxation = allow_egenrelaxation
    init_dat%use_previous_mu = use_previous_mu
    init_dat%allow_avmurelaxation = allow_avmurelaxation
    init_dat%use_quadratic_predictions = use_quadratic_predictions
    init_dat%off_centre_weight = off_centre_weight
    init_dat%ke1_2 = ie(1)
    init_dat%ke2_2 = ie(2)
    init_dat%ke3_2 = ie(3)
    init_dat%kbc_2 = ie(4)
    init_dat%kev_2 = ie(5)
    init_dat%kfn_2 = ie(6)
    init_dat%kl_2 = ie(7)
    init_dat%jh_2(1:3) = ie(8:10)
    init_dat%kp_var_2(1:40) = ie(11:50)
    init_dat%kp_eqn_2(1:40) = ie(51:90)
    init_dat%kp_bc_2(1:40) = ie(91:130)
  end subroutine push_init_dat
  
  
  !> \brief Read variables from struct init_dat
  subroutine pop_init_dat(init_dat, kh2, kr1, kr2, ksv, kt5, jch)
    use real_kind
    use mesh
    use constants
    use settings
    use control
    use extrapolate_dh
    
    implicit none
    type(init_dat_settings), intent(in) :: init_dat
    integer, intent(out) :: kh2, kr1, kr2, ksv, kt5, jch
    real(double) :: cdc_ms, cdc_hec, cdc_hes, cdc_dblsh, cdc5
    
    eps = init_dat%eps
    del = init_dat%del
    dh0 = init_dat%dh0
    wanted_eps = init_dat%wanted_eps
    cdc_ms = init_dat%cdc_ms
    cdc_hec = init_dat%cdc_hec
    cdc_hes = init_dat%cdc_hes
    cdc_dblsh = init_dat%cdc_dblsh
    cdc5 = init_dat%cdc5
    id(1) = init_dat%ke1
    id(2) = init_dat%ke2
    id(3) = init_dat%ke3
    id(4) = init_dat%kbc
    id(5) = init_dat%kev
    id(6) = init_dat%kfn
    id(7) = init_dat%kl
    id(8:10) = init_dat%jh(1:3)
    id(11:50) = init_dat%kp_var(1:40)
    id(51:90) = init_dat%kp_eqn(1:40)
    id(91:130) = init_dat%kp_bc(1:40)
    kh2 = init_dat%kh2
    kr1 = init_dat%kr1
    kr2 = init_dat%kr2
    jch = init_dat%jch
    kth = init_dat%kth
    kx = init_dat%kx
    ky = init_dat%ky
    kz = init_dat%kz
    kcl = init_dat%kcl
    kion = init_dat%kion
    kam = init_dat%kam
    kop = init_dat%kop
    kcc = init_dat%kcc
    knuc = init_dat%knuc
    kcn = init_dat%kcn
    kt1 = init_dat%kt1
    kt2 = init_dat%kt2
    kt3 = init_dat%kt3
    kt4 = init_dat%kt4
    kt5 = init_dat%kt5
    ksv = init_dat%ksv
    kn = init_dat%kn
    kjn(1:40) = init_dat%kjn(1:40)
    ct1 = init_dat%ct1
    ct2 = init_dat%ct2
    ct3 = init_dat%ct3
    ct(1:20) = init_dat%ct(1:20)
    ch = init_dat%ch
    cc = init_dat%cc
    cn = init_dat%cn
    co = init_dat%co
    cne = init_dat%cne
    cmg = init_dat%cmg
    csi = init_dat%csi
    cfe = init_dat%cfe
    calp = init_dat%calp
    cu = init_dat%cu
    cos = init_dat%cos
    cps = init_dat%cps
    crd = init_dat%crd
    cth = init_dat%cth
    cgrs = init_dat%cgrs
    ccac = init_dat%ccac
    cdsi = init_dat%cdsi
    cshi = init_dat%cshi
    cssi = init_dat%cssi
    cesc = init_dat%cesc
    cgsf = init_dat%cgsf
    cfmu = init_dat%cfmu
    cfc = init_dat%cfc
    artmix = init_dat%artmix
    cxb = init_dat%cxb
    cgr = init_dat%cgr
    cea = init_dat%cea
    cet = init_dat%cet
    cmt = init_dat%cmt
    cms = init_dat%cms
    cml = init_dat%cml
    chl = init_dat%chl
    ctf = init_dat%ctf
    clt = init_dat%clt
    cmi = init_dat%cmi
    cmr = init_dat%cmr
    cmj = init_dat%cmj
    cmv = init_dat%cmv
    cmk = init_dat%cmk
    cmnl = init_dat%cmnl
    cpa = init_dat%cpa
    cbr = init_dat%cbr
    csu = init_dat%csu
    csd = init_dat%csd
    cdf = init_dat%cdf
    cgw = init_dat%cgw
    cso = init_dat%cso
    cmb = init_dat%cmb
    cphotontire = init_dat%cphotontire
    convection_scheme = init_dat%convection_scheme
    use_fudge_control = init_dat%use_fudge_control
    allow_extension = init_dat%allow_extension
    allow_underrelaxation = init_dat%allow_underrelaxation
    allow_overrelaxation = init_dat%allow_overrelaxation
    allow_mdotrelaxation = init_dat%allow_mdotrelaxation
    allow_egenrelaxation = init_dat%allow_egenrelaxation
    use_previous_mu = init_dat%use_previous_mu
    allow_avmurelaxation = init_dat%allow_avmurelaxation
    use_quadratic_predictions = init_dat%use_quadratic_predictions
    off_centre_weight = init_dat%off_centre_weight
    ie(1) = init_dat%ke1_2
    ie(2) = init_dat%ke2_2
    ie(3) = init_dat%ke3_2
    ie(4) = init_dat%kbc_2
    ie(5) = init_dat%kev_2
    ie(6) = init_dat%kfn_2
    ie(7) = init_dat%kl_2
    ie(8:10) = init_dat%jh_2(1:3)
    ie(11:50) = init_dat%kp_var_2(1:40)
    ie(51:90) = init_dat%kp_eqn_2(1:40)
    ie(91:130) = init_dat%kp_bc_2(1:40)
  end subroutine pop_init_dat
  
  
  !> \brief Set default values for init.dat variables
  subroutine load_basic_init_dat(kh2, kr1, kr2, ksv, kt5, jch)
    use real_kind
    use mesh
    use constants
    use settings
    use control
    use extrapolate_dh
    
    implicit none
    integer, intent(out) :: kh2, kr1, kr2, ksv, kt5, jch
    integer :: ke1, ke2, ke3, kbc, kev, kfn, kl, jh(3),kp_var(40), kp_eqn(40), kp_bc(40)
    !integer :: ke1_2, ke2_2, ke3_2, kbc_2, kev_2, kfn_2, kl_2, jh_2(3),kp_var_2(40), kp_eqn_2(40), kp_bc_2(40)
    real(double) :: cdc_ms, cdc_hec, cdc_hes, cdc_dblsh, cdc5
    
    ! Load some reasonable default values
    kh2   =  kh
    jch   =    2
    kl    =    1
    jh(1:3) =   (/ 0,0,0 /)
    
    kt1   =  100
    kt2   =    2    
    kt3   =    0    
    kt4   =    1 
    kt5   =    0
    ksv   = 1000
    
    kth   = 1 
    kx    = 1 
    ky    = 1 
    kz    = 1
    
    kcl   = 1
    kion  = 5
    kam   = 1
    !kop   = 4
    kcc   = 0
    knuc  = 0
    kcn   = 0
    
    kr1   = 20
    kr2   = 5
    
    eps   = 1.00e-006
    del   = 1.00e-002
    dh0   = 1.00e-007
    
    use_quadratic_predictions = .false.
    
    wanted_eps             = 1.00e-008
    use_fudge_control      = .true.
    allow_extension        = .false.
    allow_underrelaxation  = .false.
    allow_overrelaxation   = .false.
    allow_egenrelaxation   = .false.
    allow_mdotrelaxation   = .false.
    allow_avmurelaxation   = .false.
    use_previous_mu        = .true.
    off_centre_weight      = 1.0d16
    
    cdc_ms    = 0.01
    cdc_hec   = 0.25
    cdc_hes   = 1.00
    cdc_dblsh = 4.00
    cdc5      = 1.00
    
    ct1   =  0.8
    ct2   =  1.2
    ct3   =  0.3
    
    ke1   =  6
    ke2   =  5
    ke3   =  0
    kbc   =  3
    kev   =  1
    kfn = nfunc    ! Just copy all functions; common error to forget. 
    
    kp_var(1:12)  =  (/ 7, 8, 4, 5, 3, 9,10,11,16, 1, 2, 6 /)
    kp_eqn(1:11)  =  (/ 1, 2, 3, 4, 5,13, 6, 7, 8, 9,10    /)
    kp_bc(1:12)   =  (/ 1, 2, 3, 4, 5,13, 6, 7, 8, 6, 7, 8 /)
    
    cc    =  0.176e+000
    cn    =  5.200e-002
    co    =  0.502e+000
    cne   =  9.200e-002
    cmg   =  3.400e-002
    csi   =  7.200e-002
    cfe   =  7.200e-002
    
    ct(:) = 0.0d0
    ct(1:10) =  (/ 0.0e+00, 0.0e+00, 5.0e-02, 5.0e-02, 0.15e+00,  &
         2.0e-02, 0.45e+0, 1.0e-04, 1.0e+15, 2.0e+04 /)
    
    convection_scheme = 1
    calp  =  2.000e+000
    cu    =  0.100e+000
    cos   =  0.120e+000
    cps   =  0.120e+000
    crd   =  1.000e-002
    cxb   =  0.150e+000
    cgr   =  1.000e-003
    cth   =  1.000e+000
    cgrs  =  0.000e+000
    ccac  =  0.000e+000
    
    artmix = 0.0
    
    cdsi  =  0.000e+000
    cshi  =  0.000e+000
    cssi  =  0.000e+000
    cesc  =  0.000e+000
    cgsf  =  0.000e+000
    
    cea   =  1.0e+02
    cet   =  1.0e-08
    
    cmi_mode = 1
    zscaling_mdot = 0.8
    
    cmt   =  0.0e0
    cms   =  0.0e0
    cmi   =  0.0e0
    cmr   =  0.0e0
    cmj   =  0.0e0
    cmv   =  0.0e0
    cmk   =  0.0e0
    cmnl  =  0.0e0
    cml   =  0.0e0
    chl   =  0.0e0
    ctf   =  0.0e0
    clt   =  0.0e0
    
    cpa   =  0.0e0
    cbr   =  0.0e0
    csu   =  0.0e0
    csd   =  0.0e0
    cdf   =  1.0e-02
    cgw   =  1.0e0
    cso   =  1.0e0
    cmb   =  0.0e0
    
    kn      = 12
    kjn(1:12) =  (/ 1,  2,  3,  5,  6,  7, 25, 26, 27, 29, 30, 31 /)
    
    ! Export variables to common block
    cdc(1) = cdc_ms 
    cdc(2) = cdc_hec
    cdc(3) = cdc_hes
    cdc(4) = cdc_dblsh
    cdc(5) = cdc5
    
    ! Export values to global variables
    id(1) = ke1
    id(2) = ke2
    id(3) = ke3
    id(4) = kbc
    id(5) = kev
    id(6) = kfn
    id(7) = kl
    id(8:10) = jh(1:3)
    id(11:50) = kp_var(1:40)
    id(51:90) = kp_eqn(1:40)
    id(91:130) = kp_bc(1:40)

    ie = 0
  end subroutine load_basic_init_dat
  
  
end module init_dat

