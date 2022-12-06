module AmuseInterface
    use helpers, only: InitialGenecStar,GenecStar,genec_star
    use evol, only: kindreal,ldi,npondcouche

    type(genec_star) :: StarInGenec
    type(genec_star) :: BackupStarInGenec
    integer :: number_of_stars = 0
    public:: StarInGenec,BackupStarInGenec,number_of_stars
contains

function restore_star()
    implicit none
    integer:: restore_star
    restore_star = 0
end function

function initialize_code()
    !use WriteSaveClose, only: quitafterclosing
    use genec, only: initialise_genec
    use evol, only: input_dir
    use inputparam, only: libgenec
    implicit none
    integer:: initialize_code

    libgenec = .true.
    input_dir = "./src/GENEC/code"
    call initialise_genec()

    initialize_code = 0
end function

function read_genec_model(index_of_the_star, cardfilename)
    ! This should only be called if no star has been initialised yet!
    implicit none
    integer:: index_of_the_star
    character(256):: cardfilename
    integer:: read_genec_model
    read_genec_model = 0
end function

function cleanup_code()
    implicit none
    integer:: cleanup_code
    cleanup_code = 0
end function

! **** Parameters

function get_model_number(model_number)
    use genec, only: nwmd
    implicit none
    integer:: model_number
    integer:: get_model_number
    model_number = nwmd
    get_model_number = 0
end function

function set_model_number(model_number)
    use genec, only: nwmd
    implicit none
    integer:: model_number
    integer:: set_model_number
    nwmd = model_number
    set_model_number = 0
end function

function get_min_timestep_stop_condition(min_timestep_stop_condition)
    implicit none
    real(kindreal):: min_timestep_stop_condition
    integer:: get_min_timestep_stop_condition
    get_min_timestep_stop_condition = 0
end function

function get_par_n_snap(n_snap_out)
    use inputparam, only: n_snap
    implicit none
    integer:: n_snap_out
    integer:: get_par_n_snap
    n_snap_out = n_snap
    get_par_n_snap = 0
end function

function set_par_n_snap(n_snap_in)
    use inputparam, only: n_snap
    implicit none
    integer:: n_snap_in
    integer:: set_par_n_snap
    n_snap = n_snap_in
    set_par_n_snap= 0
end function

function get_par_ipoly(par_ipoly)
    use inputparam, only: ipoly
    implicit none
    integer:: par_ipoly
    integer:: get_par_ipoly
    par_ipoly = ipoly
    get_par_ipoly = 0
end function

function set_par_ipoly(par_ipoly)
    use inputparam, only: ipoly
    implicit none
    integer:: par_ipoly
    integer:: set_par_ipoly
    ipoly = par_ipoly
    set_par_ipoly = 0
end function

function get_par_nwseq(par_nwseq)
    use inputparam, only: nwseq
    implicit none
    integer:: par_nwseq
    integer:: get_par_nwseq
    par_nwseq = nwseq
    get_par_nwseq = 0
end function

function set_par_nwseq(par_nwseq)
    use inputparam, only: nwseq
    implicit none
    integer:: par_nwseq
    integer:: set_par_nwseq
    nwseq = par_nwseq
    set_par_nwseq = 0
end function

function get_par_modanf(par_modanf)
    use inputparam, only: modanf
    implicit none
    integer:: par_modanf
    integer:: get_par_modanf
    par_modanf = modanf
    get_par_modanf = 0
end function

function set_par_modanf(par_modanf)
    use inputparam, only: modanf
    implicit none
    integer:: par_modanf
    integer:: set_par_modanf
    modanf = par_modanf
    set_par_modanf = 0
end function

function get_par_nzmod(par_nzmod)
    use inputparam, only: nzmod
    implicit none
    integer:: par_nzmod
    integer:: get_par_nzmod
    par_nzmod = nzmod
    get_par_nzmod = 0
end function

function set_par_nzmod(par_nzmod)
    use inputparam, only: nzmod
    implicit none
    integer:: par_nzmod
    integer:: set_par_nzmod
    nzmod = par_nzmod
    set_par_nzmod = 0
end function

function get_par_irot(par_irot)
    use inputparam, only: irot
    implicit none
    integer:: par_irot
    integer:: get_par_irot
    par_irot = irot
    get_par_irot = 0
end function

function set_par_irot(par_irot)
    use inputparam, only: irot
    implicit none
    integer:: par_irot
    integer:: set_par_irot
    irot = par_irot
    set_par_irot = 0
end function

function get_par_isol(par_isol)
    use inputparam, only: isol
    implicit none
    integer:: par_isol
    integer:: get_par_isol
    par_isol = isol
    get_par_isol = 0
end function

function set_par_isol(par_isol)
    use inputparam, only: isol
    implicit none
    integer:: par_isol
    integer:: set_par_isol
    isol = par_isol
    set_par_isol = 0
end function

function get_par_imagn(par_imagn)
    use inputparam, only: imagn
    implicit none
    integer:: par_imagn
    integer:: get_par_imagn
    par_imagn = imagn
    get_par_imagn = 0
end function

function set_par_imagn(par_imagn)
    use inputparam, only: imagn
    implicit none
    integer:: par_imagn
    integer:: set_par_imagn
    imagn = par_imagn
    set_par_imagn = 0
end function

function get_par_ialflu(par_ialflu)
    use inputparam, only: ialflu
    implicit none
    integer:: par_ialflu
    integer:: get_par_ialflu
    par_ialflu = ialflu
    get_par_ialflu = 0
end function

function set_par_ialflu(par_ialflu)
    use inputparam, only: ialflu
    implicit none
    integer:: par_ialflu
    integer:: set_par_ialflu
    ialflu = par_ialflu
    set_par_ialflu = 0
end function

function get_par_ianiso(par_ianiso)
    use inputparam, only: ianiso
    implicit none
    integer:: par_ianiso
    integer:: get_par_ianiso
    par_ianiso = ianiso
    get_par_ianiso = 0
end function

function set_par_ianiso(par_ianiso)
    use inputparam, only: ianiso
    implicit none
    integer:: par_ianiso
    integer:: set_par_ianiso
    ianiso = par_ianiso
    set_par_ianiso = 0
end function

function get_par_ipop3(par_ipop3)
    use inputparam, only: ipop3
    implicit none
    integer:: par_ipop3
    integer:: get_par_ipop3
    par_ipop3 = ipop3
    get_par_ipop3 = 0
end function

function set_par_ipop3(par_ipop3)
    use inputparam, only: ipop3
    implicit none
    integer:: par_ipop3
    integer:: set_par_ipop3
    ipop3 = par_ipop3
    set_par_ipop3 = 0
end function

function get_par_ibasnet(par_ibasnet)
    use inputparam, only: ibasnet
    implicit none
    integer:: par_ibasnet
    integer:: get_par_ibasnet
    par_ibasnet = ibasnet
    get_par_ibasnet = 0
end function

function set_par_ibasnet(par_ibasnet)
    use inputparam, only: ibasnet
    implicit none
    integer:: par_ibasnet
    integer:: set_par_ibasnet
    ibasnet = par_ibasnet
    set_par_ibasnet = 0
end function

function get_par_phase(par_phase)
    use inputparam, only: phase
    implicit none
    integer:: par_phase
    integer:: get_par_phase
    par_phase = phase
    get_par_phase = 0
end function

function set_par_phase(par_phase)
    use inputparam, only: phase
    implicit none
    integer:: par_phase
    integer:: set_par_phase
    phase = par_phase
    set_par_phase = 0
end function

function get_par_iopac(par_iopac)
    use inputparam, only: iopac
    implicit none
    integer:: par_iopac
    integer:: get_par_iopac
    par_iopac = iopac
    get_par_iopac = 0
end function

function set_par_iopac(par_iopac)
    use inputparam, only: iopac
    implicit none
    integer:: par_iopac
    integer:: set_par_iopac
    iopac = par_iopac
    set_par_iopac = 0
end function

function get_par_ikappa(par_ikappa)
    use inputparam, only: ikappa
    implicit none
    integer:: par_ikappa
    integer:: get_par_ikappa
    par_ikappa = ikappa
    get_par_ikappa = 0
end function

function set_par_ikappa(par_ikappa)
    use inputparam, only: ikappa
    implicit none
    integer:: par_ikappa
    integer:: set_par_ikappa
    ikappa = par_ikappa
    set_par_ikappa = 0
end function

function get_par_idiff(par_idiff)
    use inputparam, only: idiff
    implicit none
    integer:: par_idiff
    integer:: get_par_idiff
    par_idiff = idiff
    get_par_idiff = 0
end function

function set_par_idiff(par_idiff)
    use inputparam, only: idiff
    implicit none
    integer:: par_idiff
    integer:: set_par_idiff
    idiff = par_idiff
    set_par_idiff = 0
end function

function get_par_iadvec(par_iadvec)
    use inputparam, only: iadvec
    implicit none
    integer:: par_iadvec
    integer:: get_par_iadvec
    par_iadvec = iadvec
    get_par_iadvec = 0
end function

function set_par_iadvec(par_iadvec)
    use inputparam, only: iadvec
    implicit none
    integer:: par_iadvec
    integer:: set_par_iadvec
    iadvec = par_iadvec
    set_par_iadvec = 0
end function

function get_par_istati(par_istati)
    use inputparam, only: istati
    implicit none
    integer:: par_istati
    integer:: get_par_istati
    par_istati = istati
    get_par_istati = 0
end function

function set_par_istati(par_istati)
    use inputparam, only: istati
    implicit none
    integer:: par_istati
    integer:: set_par_istati
    istati = par_istati
    set_par_istati = 0
end function

function get_par_icoeff(par_icoeff)
    use inputparam, only: icoeff
    implicit none
    integer:: par_icoeff
    integer:: get_par_icoeff
    par_icoeff = icoeff
    get_par_icoeff = 0
end function

function set_par_icoeff(par_icoeff)
    use inputparam, only: icoeff
    implicit none
    integer:: par_icoeff
    integer:: set_par_icoeff
    icoeff = par_icoeff
    set_par_icoeff = 0
end function

function get_par_igamma(par_igamma)
    use inputparam, only: igamma
    implicit none
    integer:: par_igamma
    integer:: get_par_igamma
    par_igamma = igamma
    get_par_igamma = 0
end function

function set_par_igamma(par_igamma)
    use inputparam, only: igamma
    implicit none
    integer:: par_igamma
    integer:: set_par_igamma
    igamma = par_igamma
    set_par_igamma = 0
end function

function get_par_idialo(par_idialo)
    use inputparam, only: idialo
    implicit none
    integer:: par_idialo
    integer:: get_par_idialo
    par_idialo = idialo
    get_par_idialo = 0
end function

function set_par_idialo(par_idialo)
    use inputparam, only: idialo
    implicit none
    integer:: par_idialo
    integer:: set_par_idialo
    idialo = par_idialo
    set_par_idialo = 0
end function

function get_par_idialu(par_idialu)
    use inputparam, only: idialu
    implicit none
    integer:: par_idialu
    integer:: get_par_idialu
    par_idialu = idialu
    get_par_idialu = 0
end function

function set_par_idialu(par_idialu)
    use inputparam, only: idialu
    implicit none
    integer:: par_idialu
    integer:: set_par_idialu
    idialu = par_idialu
    set_par_idialu = 0
end function

function get_par_imloss(par_imloss)
    use inputparam, only: imloss
    implicit none
    integer:: par_imloss
    integer:: get_par_imloss
    par_imloss = imloss
    get_par_imloss = 0
end function

function set_par_imloss(par_imloss)
    use inputparam, only: imloss
    implicit none
    integer:: par_imloss
    integer:: set_par_imloss
    imloss = par_imloss
    set_par_imloss = 0
end function

function get_par_ifitm(par_ifitm)
    use inputparam, only: ifitm
    implicit none
    integer:: par_ifitm
    integer:: get_par_ifitm
    par_ifitm = ifitm
    get_par_ifitm = 0
end function

function set_par_ifitm(par_ifitm)
    use inputparam, only: ifitm
    implicit none
    integer:: par_ifitm
    integer:: set_par_ifitm
    ifitm = par_ifitm
    set_par_ifitm = 0
end function

function get_par_nndr(par_nndr)
    use inputparam, only: nndr
    implicit none
    integer:: par_nndr
    integer:: get_par_nndr
    par_nndr = nndr
    get_par_nndr = 0
end function

function set_par_nndr(par_nndr)
    use inputparam, only: nndr
    implicit none
    integer:: par_nndr
    integer:: set_par_nndr
    nndr = par_nndr
    set_par_nndr = 0
end function

function get_par_iledou(par_iledou)
    use inputparam, only: iledou
    implicit none
    integer:: par_iledou
    integer:: get_par_iledou
    par_iledou = iledou
    get_par_iledou = 0
end function

function set_par_iledou(par_iledou)
    use inputparam, only: iledou
    implicit none
    integer:: par_iledou
    integer:: set_par_iledou
    iledou = par_iledou
    set_par_iledou = 0
end function

function get_par_idifcon(par_idifcon)
    use inputparam, only: idifcon
    implicit none
    integer:: par_idifcon
    integer:: get_par_idifcon
    par_idifcon = idifcon
    get_par_idifcon = 0
end function

function set_par_idifcon(par_idifcon)
    use inputparam, only: idifcon
    implicit none
    integer:: par_idifcon
    integer:: set_par_idifcon
    idifcon = par_idifcon
    set_par_idifcon = 0
end function

function get_par_my(par_my)
    use inputparam, only: my
    implicit none
    integer:: par_my
    integer:: get_par_my
    par_my = my
    get_par_my = 0
end function

function set_par_my(par_my)
    use inputparam, only: my
    implicit none
    integer:: par_my
    integer:: set_par_my
    my = par_my
    set_par_my = 0
end function

function get_par_iover(par_iover)
    use inputparam, only: iover
    implicit none
    integer:: par_iover
    integer:: get_par_iover
    par_iover = iover
    get_par_iover = 0
end function

function set_par_iover(par_iover)
    use inputparam, only: iover
    implicit none
    integer:: par_iover
    integer:: set_par_iover
    iover = par_iover
    set_par_iover = 0
end function

function get_par_iunder(par_iunder)
    use inputparam, only: iunder
    implicit none
    integer:: par_iunder
    integer:: get_par_iunder
    par_iunder = iunder
    get_par_iunder = 0
end function

function set_par_iunder(par_iunder)
    use inputparam, only: iunder
    implicit none
    integer:: par_iunder
    integer:: set_par_iunder
    iunder = par_iunder
    set_par_iunder = 0
end function

function get_par_nbchx(par_nbchx)
    use inputparam, only: nbchx
    implicit none
    integer:: par_nbchx
    integer:: get_par_nbchx
    par_nbchx = nbchx
    get_par_nbchx = 0
end function

function set_par_nbchx(par_nbchx)
    use inputparam, only: nbchx
    implicit none
    integer:: par_nbchx
    integer:: set_par_nbchx
    nbchx = par_nbchx
    set_par_nbchx = 0
end function

function get_par_nrband(par_nrband)
    use inputparam, only: nrband
    implicit none
    integer:: par_nrband
    integer:: get_par_nrband
    par_nrband = nrband
    get_par_nrband = 0
end function

function set_par_nrband(par_nrband)
    use inputparam, only: nrband
    implicit none
    integer:: par_nrband
    integer:: set_par_nrband
    nrband = par_nrband
    set_par_nrband = 0
end function

function get_par_islow(par_islow)
    use inputparam, only: islow
    implicit none
    integer:: par_islow
    integer:: get_par_islow
    par_islow = islow
    get_par_islow = 0
end function

function set_par_islow(par_islow)
    use inputparam, only: islow
    implicit none
    integer:: par_islow
    integer:: set_par_islow
    islow = par_islow
    set_par_islow = 0
end function

function get_par_icncst(par_icncst)
    use inputparam, only: icncst
    implicit none
    integer:: par_icncst
    integer:: get_par_icncst
    par_icncst = icncst
    get_par_icncst = 0
end function

function set_par_icncst(par_icncst)
    use inputparam, only: icncst
    implicit none
    integer:: par_icncst
    integer:: set_par_icncst
    icncst = par_icncst
    set_par_icncst = 0
end function

function get_par_tauH_fit(par_tauH_fit)
    use inputparam, only: tauH_fit
    implicit none
    integer:: par_tauH_fit
    integer:: get_par_tauH_fit
    par_tauH_fit = tauH_fit
    get_par_tauH_fit = 0
end function

function set_par_tauH_fit(par_tauH_fit)
    use inputparam, only: tauH_fit
    implicit none
    integer:: par_tauH_fit
    integer:: set_par_tauH_fit
    tauH_fit = par_tauH_fit
    set_par_tauH_fit = 0
end function

function get_par_iauto(par_iauto)
    use inputparam, only: iauto
    implicit none
    integer:: par_iauto
    integer:: get_par_iauto
    par_iauto = iauto
    get_par_iauto = 0
end function

function set_par_iauto(par_iauto)
    use inputparam, only: iauto
    implicit none
    integer:: par_iauto
    integer:: set_par_iauto
    iauto = par_iauto
    set_par_iauto = 0
end function

function get_par_iprn(par_iprn)
    use inputparam, only: iprn
    implicit none
    integer:: par_iprn
    integer:: get_par_iprn
    par_iprn = iprn
    get_par_iprn = 0
end function

function set_par_iprn(par_iprn)
    use inputparam, only: iprn
    implicit none
    integer:: par_iprn
    integer:: set_par_iprn
    iprn = par_iprn
    set_par_iprn = 0
end function

function get_par_iout(par_iout)
    use inputparam, only: iout
    implicit none
    integer:: par_iout
    integer:: get_par_iout
    par_iout = iout
    get_par_iout = 0
end function

function set_par_iout(par_iout)
    use inputparam, only: iout
    implicit none
    integer:: par_iout
    integer:: set_par_iout
    iout = par_iout
    set_par_iout = 0
end function

function get_par_itmin(par_itmin)
    use inputparam, only: itmin
    implicit none
    integer:: par_itmin
    integer:: get_par_itmin
    par_itmin = itmin
    get_par_itmin = 0
end function

function set_par_itmin(par_itmin)
    use inputparam, only: itmin
    implicit none
    integer:: par_itmin
    integer:: set_par_itmin
    itmin = par_itmin
    set_par_itmin = 0
end function

function get_par_idebug(par_idebug)
    use inputparam, only: idebug
    implicit none
    integer:: par_idebug
    integer:: get_par_idebug
    par_idebug = idebug
    get_par_idebug = 0
end function

function set_par_idebug(par_idebug)
    use inputparam, only: idebug
    implicit none
    integer:: par_idebug
    integer:: set_par_idebug
    idebug = par_idebug
    set_par_idebug = 0
end function

function get_par_itests(par_itests)
    use inputparam, only: itests
    implicit none
    integer:: par_itests
    integer:: get_par_itests
    par_itests = itests
    get_par_itests = 0
end function

function set_par_itests(par_itests)
    use inputparam, only: itests
    implicit none
    integer:: par_itests
    integer:: set_par_itests
    itests = par_itests
    set_par_itests = 0
end function

function get_par_var_rates(par_var_rates)
    use inputparam, only: var_rates
    implicit none
    logical:: par_var_rates
    integer:: get_par_var_rates
    par_var_rates = var_rates
    get_par_var_rates = 0
end function

function set_par_var_rates(par_var_rates)
    use inputparam, only: var_rates
    implicit none
    logical:: par_var_rates
    integer:: set_par_var_rates
    var_rates = par_var_rates
    set_par_var_rates = 0
end function

function get_par_bintide(par_bintide)
    use inputparam, only: bintide
    implicit none
    logical:: par_bintide
    integer:: get_par_bintide
    par_bintide = bintide
    get_par_bintide = 0
end function

function set_par_bintide(par_bintide)
    use inputparam, only: bintide
    implicit none
    logical:: par_bintide
    integer:: set_par_bintide
    bintide = par_bintide
    set_par_bintide = 0
end function

function get_par_const_per(par_const_per)
    use inputparam, only: const_per
    implicit none
    logical:: par_const_per
    integer:: get_par_const_per
    par_const_per = const_per
    get_par_const_per = 0
end function

function set_par_const_per(par_const_per)
    use inputparam, only: const_per
    implicit none
    logical:: par_const_per
    integer:: set_par_const_per
    const_per = par_const_per
    set_par_const_per = 0
end function

function get_par_Add_Flux(par_Add_Flux)
    use inputparam, only: Add_Flux
    implicit none
    logical:: par_Add_Flux
    integer:: get_par_Add_Flux
    par_Add_Flux = Add_Flux
    get_par_Add_Flux = 0
end function

function set_par_Add_Flux(par_Add_Flux)
    use inputparam, only: Add_Flux
    implicit none
    logical:: par_Add_Flux
    integer:: set_par_Add_Flux
    Add_Flux = par_Add_Flux
    set_par_Add_Flux = 0
end function

function get_par_diff_only(par_diff_only)
    use inputparam, only: diff_only
    implicit none
    logical:: par_diff_only
    integer:: get_par_diff_only
    par_diff_only = diff_only
    get_par_diff_only = 0
end function

function set_par_diff_only(par_diff_only)
    use inputparam, only: diff_only
    implicit none
    logical:: par_diff_only
    integer:: set_par_diff_only
    diff_only = par_diff_only
    set_par_diff_only = 0
end function

function get_par_RSG_Mdot(par_RSG_Mdot)
    use inputparam, only: RSG_Mdot
    implicit none
    integer:: par_RSG_Mdot
    integer:: get_par_RSG_Mdot
    par_RSG_Mdot = RSG_Mdot
    get_par_RSG_Mdot = 0
end function

function set_par_RSG_Mdot(par_RSG_Mdot)
    use inputparam, only: RSG_Mdot
    implicit none
    integer:: par_RSG_Mdot
    integer:: set_par_RSG_Mdot
    RSG_Mdot = par_RSG_Mdot
    set_par_RSG_Mdot = 0
end function

function get_par_display_plot(par_display_plot)
    use inputparam, only: display_plot
    implicit none
    logical:: par_display_plot
    integer:: get_par_display_plot
    par_display_plot = display_plot
    get_par_display_plot = 0
end function

function set_par_display_plot(par_display_plot)
    use inputparam, only: display_plot
    implicit none
    logical:: par_display_plot
    integer:: set_par_display_plot
    display_plot = par_display_plot
    set_par_display_plot = 0
end function

function get_par_xyfiles(par_xyfiles)
    use inputparam, only: xyfiles
    implicit none
    logical:: par_xyfiles
    integer:: get_par_xyfiles
    par_xyfiles = xyfiles
    get_par_xyfiles = 0
end function

function set_par_xyfiles(par_xyfiles)
    use inputparam, only: xyfiles
    implicit none
    logical:: par_xyfiles
    integer:: set_par_xyfiles
    xyfiles = par_xyfiles
    set_par_xyfiles = 0
end function

function get_par_verbose(par_verbose)
    use inputparam, only: verbose
    implicit none
    logical:: par_verbose
    integer:: get_par_verbose
    par_verbose = verbose
    get_par_verbose = 0
end function

function set_par_verbose(par_verbose)
    use inputparam, only: verbose
    implicit none
    logical:: par_verbose
    integer:: set_par_verbose
    verbose = par_verbose
    set_par_verbose = 0
end function

function get_par_stop_deg(par_stop_deg)
    use inputparam, only: stop_deg
    implicit none
    logical:: par_stop_deg
    integer:: get_par_stop_deg
    par_stop_deg = stop_deg
    get_par_stop_deg = 0
end function

function set_par_stop_deg(par_stop_deg)
    use inputparam, only: stop_deg
    implicit none
    logical:: par_stop_deg
    integer:: set_par_stop_deg
    stop_deg = par_stop_deg
    set_par_stop_deg = 0
end function

function get_par_index_poly(par_index_poly)
    use inputparam, only: index_poly
    implicit none
    real(kindreal):: par_index_poly
    integer:: get_par_index_poly
    par_index_poly = index_poly
    get_par_index_poly = 0
end function

function set_par_index_poly(par_index_poly)
    use inputparam, only: index_poly
    implicit none
    real(kindreal):: par_index_poly
    integer:: set_par_index_poly
    index_poly = par_index_poly
    set_par_index_poly = 0
end function

function get_par_binm2(par_binm2)
    use inputparam, only: binm2
    implicit none
    real(kindreal):: par_binm2
    integer:: get_par_binm2
    par_binm2 = binm2
    get_par_binm2 = 0
end function

function set_par_binm2(par_binm2)
    use inputparam, only: binm2
    implicit none
    real(kindreal):: par_binm2
    integer:: set_par_binm2
    binm2 = par_binm2
    set_par_binm2 = 0
end function

function get_par_periodini(par_periodini)
    use inputparam, only: periodini
    implicit none
    real(kindreal):: par_periodini
    integer:: get_par_periodini
    par_periodini = periodini
    get_par_periodini = 0
end function

function set_par_periodini(par_periodini)
    use inputparam, only: periodini
    implicit none
    real(kindreal):: par_periodini
    integer:: set_par_periodini
    periodini = par_periodini
    set_par_periodini = 0
end function

function get_par_zinit(par_zinit)
    use inputparam, only: zinit
    implicit none
    real(kindreal):: par_zinit
    integer:: get_par_zinit
    par_zinit = zinit
    get_par_zinit = 0
end function

function set_par_zinit(par_zinit)
    use inputparam, only: zinit
    implicit none
    real(kindreal):: par_zinit
    integer:: set_par_zinit
    zinit = par_zinit
    set_par_zinit = 0
end function

function get_par_zsol(par_zsol)
    use inputparam, only: zsol
    implicit none
    real(kindreal):: par_zsol
    integer:: get_par_zsol
    par_zsol = zsol
    get_par_zsol = 0
end function

function set_par_zsol(par_zsol)
    use inputparam, only: zsol
    implicit none
    real(kindreal):: par_zsol
    integer:: set_par_zsol
    zsol = par_zsol
    set_par_zsol = 0
end function

function get_par_z(par_z)
    use inputparam, only: z
    implicit none
    real(kindreal):: par_z
    integer:: get_par_z
    par_z = z
    get_par_z = 0
end function

function set_par_z(par_z)
    use inputparam, only: z
    implicit none
    real(kindreal):: par_z
    integer:: set_par_z
    z = par_z
    set_par_z = 0
end function

function get_par_fenerg(par_fenerg)
    use inputparam, only: fenerg
    implicit none
    real(kindreal):: par_fenerg
    integer:: get_par_fenerg
    par_fenerg = fenerg
    get_par_fenerg = 0
end function

function set_par_fenerg(par_fenerg)
    use inputparam, only: fenerg
    implicit none
    real(kindreal):: par_fenerg
    integer:: set_par_fenerg
    fenerg = par_fenerg
    set_par_fenerg = 0
end function

function get_par_richac(par_richac)
    use inputparam, only: richac
    implicit none
    real(kindreal):: par_richac
    integer:: get_par_richac
    par_richac = richac
    get_par_richac = 0
end function

function set_par_richac(par_richac)
    use inputparam, only: richac
    implicit none
    real(kindreal):: par_richac
    integer:: set_par_richac
    richac = par_richac
    set_par_richac = 0
end function

function get_par_frein(par_frein)
    use inputparam, only: frein
    implicit none
    real(kindreal):: par_frein
    integer:: get_par_frein
    par_frein = frein
    get_par_frein = 0
end function

function set_par_frein(par_frein)
    use inputparam, only: frein
    implicit none
    real(kindreal):: par_frein
    integer:: set_par_frein
    frein = par_frein
    set_par_frein = 0
end function

function get_par_K_Kawaler(par_K_Kawaler)
    use inputparam, only: K_Kawaler
    implicit none
    real(kindreal):: par_K_Kawaler
    integer:: get_par_K_Kawaler
    par_K_Kawaler = K_Kawaler
    get_par_K_Kawaler = 0
end function

function set_par_K_Kawaler(par_K_Kawaler)
    use inputparam, only: K_Kawaler
    implicit none
    real(kindreal):: par_K_Kawaler
    integer:: set_par_K_Kawaler
    K_Kawaler = par_K_Kawaler
    set_par_K_Kawaler = 0
end function

function get_par_Omega_saturation(par_Omega_saturation)
    use inputparam, only: Omega_saturation
    implicit none
    real(kindreal):: par_Omega_saturation
    integer:: get_par_Omega_saturation
    par_Omega_saturation = Omega_saturation
    get_par_Omega_saturation = 0
end function

function set_par_Omega_saturation(par_Omega_saturation)
    use inputparam, only: Omega_saturation
    implicit none
    real(kindreal):: par_Omega_saturation
    integer:: set_par_Omega_saturation
    Omega_saturation = par_Omega_saturation
    set_par_Omega_saturation = 0
end function

function get_par_rapcrilim(par_rapcrilim)
    use inputparam, only: rapcrilim
    implicit none
    real(kindreal):: par_rapcrilim
    integer:: get_par_rapcrilim
    par_rapcrilim = rapcrilim
    get_par_rapcrilim = 0
end function

function set_par_rapcrilim(par_rapcrilim)
    use inputparam, only: rapcrilim
    implicit none
    real(kindreal):: par_rapcrilim
    integer:: set_par_rapcrilim
    rapcrilim = par_rapcrilim
    set_par_rapcrilim = 0
end function

function get_par_vwant(par_vwant)
    use inputparam, only: vwant
    implicit none
    real(kindreal):: par_vwant
    integer:: get_par_vwant
    par_vwant = vwant
    get_par_vwant = 0
end function

function set_par_vwant(par_vwant)
    use inputparam, only: vwant
    implicit none
    real(kindreal):: par_vwant
    integer:: set_par_vwant
    vwant = par_vwant
    set_par_vwant = 0
end function

function get_par_xfom(par_xfom)
    use inputparam, only: xfom
    implicit none
    real(kindreal):: par_xfom
    integer:: get_par_xfom
    par_xfom = xfom
    get_par_xfom = 0
end function

function set_par_xfom(par_xfom)
    use inputparam, only: xfom
    implicit none
    real(kindreal):: par_xfom
    integer:: set_par_xfom
    xfom = par_xfom
    set_par_xfom = 0
end function

function get_par_omega(par_omega)
    use inputparam, only: omega
    implicit none
    real(kindreal):: par_omega
    integer:: get_par_omega
    par_omega = omega
    get_par_omega = 0
end function

function set_par_omega(par_omega)
    use inputparam, only: omega
    implicit none
    real(kindreal):: par_omega
    integer:: set_par_omega
    omega = par_omega
    set_par_omega = 0
end function

function get_par_xdial(par_xdial)
    use inputparam, only: xdial
    implicit none
    real(kindreal):: par_xdial
    integer:: get_par_xdial
    par_xdial = xdial
    get_par_xdial = 0
end function

function set_par_xdial(par_xdial)
    use inputparam, only: xdial
    implicit none
    real(kindreal):: par_xdial
    integer:: set_par_xdial
    xdial = par_xdial
    set_par_xdial = 0
end function

function get_par_B_initial(par_B_initial)
    use inputparam, only: B_initial
    implicit none
    real(kindreal):: par_B_initial
    integer:: get_par_B_initial
    par_B_initial = B_initial
    get_par_B_initial = 0
end function

function set_par_B_initial(par_B_initial)
    use inputparam, only: B_initial
    implicit none
    real(kindreal):: par_B_initial
    integer:: set_par_B_initial
    B_initial = par_B_initial
    set_par_B_initial = 0
end function

function get_par_add_diff(par_add_diff)
    use inputparam, only: add_diff
    implicit none
    real(kindreal):: par_add_diff
    integer:: get_par_add_diff
    par_add_diff = add_diff
    get_par_add_diff = 0
end function

function set_par_add_diff(par_add_diff)
    use inputparam, only: add_diff
    implicit none
    real(kindreal):: par_add_diff
    integer:: set_par_add_diff
    add_diff = par_add_diff
    set_par_add_diff = 0
end function

function get_par_fmlos(par_fmlos)
    use inputparam, only: fmlos
    implicit none
    real(kindreal):: par_fmlos
    integer:: get_par_fmlos
    par_fmlos = fmlos
    get_par_fmlos = 0
end function

function set_par_fmlos(par_fmlos)
    use inputparam, only: fmlos
    implicit none
    real(kindreal):: par_fmlos
    integer:: set_par_fmlos
    fmlos = par_fmlos
    set_par_fmlos = 0
end function

function get_par_fitm(par_fitm)
    use inputparam, only: fitm
    implicit none
    real(kindreal):: par_fitm
    integer:: get_par_fitm
    par_fitm = fitm
    get_par_fitm = 0
end function

function set_par_fitm(par_fitm)
    use inputparam, only: fitm
    implicit none
    real(kindreal):: par_fitm
    integer:: set_par_fitm
    fitm = par_fitm
    set_par_fitm = 0
end function

function get_par_fitmi(par_fitmi)
    use inputparam, only: fitmi
    implicit none
    real(kindreal):: par_fitmi
    integer:: get_par_fitmi
    par_fitmi = fitmi
    get_par_fitmi = 0
end function

function set_par_fitmi(par_fitmi)
    use inputparam, only: fitmi
    implicit none
    real(kindreal):: par_fitmi
    integer:: set_par_fitmi
    fitmi = par_fitmi
    set_par_fitmi = 0
end function

function get_par_fitmi_default(par_fitmi_default)
    use inputparam, only: fitmi_default
    implicit none
    real(kindreal):: par_fitmi_default
    integer:: get_par_fitmi_default
    par_fitmi_default = fitmi_default
    get_par_fitmi_default = 0
end function

function set_par_fitmi_default(par_fitmi_default)
    use inputparam, only: fitmi_default
    implicit none
    real(kindreal):: par_fitmi_default
    integer:: set_par_fitmi_default
    fitmi_default = par_fitmi_default
    set_par_fitmi_default = 0
end function

function get_par_deltal(par_deltal)
    use inputparam, only: deltal
    implicit none
    real(kindreal):: par_deltal
    integer:: get_par_deltal
    par_deltal = deltal
    get_par_deltal = 0
end function

function set_par_deltal(par_deltal)
    use inputparam, only: deltal
    implicit none
    real(kindreal):: par_deltal
    integer:: set_par_deltal
    deltal = par_deltal
    set_par_deltal = 0
end function

function get_par_deltat(par_deltat)
    use inputparam, only: deltat
    implicit none
    real(kindreal):: par_deltat
    integer:: get_par_deltat
    par_deltat = deltat
    get_par_deltat = 0
end function

function set_par_deltat(par_deltat)
    use inputparam, only: deltat
    implicit none
    real(kindreal):: par_deltat
    integer:: set_par_deltat
    deltat = par_deltat
    set_par_deltat = 0
end function

function get_par_elph(par_elph)
    use inputparam, only: elph
    implicit none
    real(kindreal):: par_elph
    integer:: get_par_elph
    par_elph = elph
    get_par_elph = 0
end function

function set_par_elph(par_elph)
    use inputparam, only: elph
    implicit none
    real(kindreal):: par_elph
    integer:: set_par_elph
    elph = par_elph
    set_par_elph = 0
end function

function get_par_dovhp(par_dovhp)
    use inputparam, only: dovhp
    implicit none
    real(kindreal):: par_dovhp
    integer:: get_par_dovhp
    par_dovhp = dovhp
    get_par_dovhp = 0
end function

function set_par_dovhp(par_dovhp)
    use inputparam, only: dovhp
    implicit none
    real(kindreal):: par_dovhp
    integer:: set_par_dovhp
    dovhp = par_dovhp
    set_par_dovhp = 0
end function

function get_par_dunder(par_dunder)
    use inputparam, only: dunder
    implicit none
    real(kindreal):: par_dunder
    integer:: get_par_dunder
    par_dunder = dunder
    get_par_dunder = 0
end function

function set_par_dunder(par_dunder)
    use inputparam, only: dunder
    implicit none
    real(kindreal):: par_dunder
    integer:: set_par_dunder
    dunder = par_dunder
    set_par_dunder = 0
end function

function get_par_gkorm(par_gkorm)
    use inputparam, only: gkorm
    implicit none
    real(kindreal):: par_gkorm
    integer:: get_par_gkorm
    par_gkorm = gkorm
    get_par_gkorm = 0
end function

function set_par_gkorm(par_gkorm)
    use inputparam, only: gkorm
    implicit none
    real(kindreal):: par_gkorm
    integer:: set_par_gkorm
    gkorm = par_gkorm
    set_par_gkorm = 0
end function

function get_par_alph(par_alph)
    use inputparam, only: alph
    implicit none
    real(kindreal):: par_alph
    integer:: get_par_alph
    par_alph = alph
    get_par_alph = 0
end function

function set_par_alph(par_alph)
    use inputparam, only: alph
    implicit none
    real(kindreal):: par_alph
    integer:: set_par_alph
    alph = par_alph
    set_par_alph = 0
end function

function get_par_agdr(par_agdr)
    use inputparam, only: agdr
    implicit none
    real(kindreal):: par_agdr
    integer:: get_par_agdr
    par_agdr = agdr
    get_par_agdr = 0
end function

function set_par_agdr(par_agdr)
    use inputparam, only: agdr
    implicit none
    real(kindreal):: par_agdr
    integer:: set_par_agdr
    agdr = par_agdr
    set_par_agdr = 0
end function

function get_par_faktor(par_faktor)
    use inputparam, only: faktor
    implicit none
    real(kindreal):: par_faktor
    integer:: get_par_faktor
    par_faktor = faktor
    get_par_faktor = 0
end function

function set_par_faktor(par_faktor)
    use inputparam, only: faktor
    implicit none
    real(kindreal):: par_faktor
    integer:: set_par_faktor
    faktor = par_faktor
    set_par_faktor = 0
end function

function get_par_dgrp(par_dgrp)
    use inputparam, only: dgrp
    implicit none
    real(kindreal):: par_dgrp
    integer:: get_par_dgrp
    par_dgrp = dgrp
    get_par_dgrp = 0
end function

function set_par_dgrp(par_dgrp)
    use inputparam, only: dgrp
    implicit none
    real(kindreal):: par_dgrp
    integer:: set_par_dgrp
    dgrp = par_dgrp
    set_par_dgrp = 0
end function

function get_par_dgrl(par_dgrl)
    use inputparam, only: dgrl
    implicit none
    real(kindreal):: par_dgrl
    integer:: get_par_dgrl
    par_dgrl = dgrl
    get_par_dgrl = 0
end function

function set_par_dgrl(par_dgrl)
    use inputparam, only: dgrl
    implicit none
    real(kindreal):: par_dgrl
    integer:: set_par_dgrl
    dgrl = par_dgrl
    set_par_dgrl = 0
end function

function get_par_dgry(par_dgry)
    use inputparam, only: dgry
    implicit none
    real(kindreal):: par_dgry
    integer:: get_par_dgry
    par_dgry = dgry
    get_par_dgry = 0
end function

function set_par_dgry(par_dgry)
    use inputparam, only: dgry
    implicit none
    real(kindreal):: par_dgry
    integer:: set_par_dgry
    dgry = par_dgry
    set_par_dgry = 0
end function

function get_par_dgrc(par_dgrc)
    use inputparam, only: dgrc
    implicit none
    real(kindreal):: par_dgrc
    integer:: get_par_dgrc
    par_dgrc = dgrc
    get_par_dgrc = 0
end function

function set_par_dgrc(par_dgrc)
    use inputparam, only: dgrc
    implicit none
    real(kindreal):: par_dgrc
    integer:: set_par_dgrc
    dgrc = par_dgrc
    set_par_dgrc = 0
end function

function get_par_dgro(par_dgro)
    use inputparam, only: dgro
    implicit none
    real(kindreal):: par_dgro
    integer:: get_par_dgro
    par_dgro = dgro
    get_par_dgro = 0
end function

function set_par_dgro(par_dgro)
    use inputparam, only: dgro
    implicit none
    real(kindreal):: par_dgro
    integer:: set_par_dgro
    dgro = par_dgro
    set_par_dgro = 0
end function

function get_par_dgr20(par_dgr20)
    use inputparam, only: dgr20
    implicit none
    real(kindreal):: par_dgr20
    integer:: get_par_dgr20
    par_dgr20 = dgr20
    get_par_dgr20 = 0
end function

function set_par_dgr20(par_dgr20)
    use inputparam, only: dgr20
    implicit none
    real(kindreal):: par_dgr20
    integer:: set_par_dgr20
    dgr20 = par_dgr20
    set_par_dgr20 = 0
end function

function get_par_xcn(par_xcn)
    use inputparam, only: xcn
    implicit none
    real(kindreal):: par_xcn
    integer:: get_par_xcn
    par_xcn = xcn
    get_par_xcn = 0
end function

function set_par_xcn(par_xcn)
    use inputparam, only: xcn
    implicit none
    real(kindreal):: par_xcn
    integer:: set_par_xcn
    xcn = par_xcn
    set_par_xcn = 0
end function

function get_par_starname(par_starname)
    use inputparam, only: starname
    implicit none
    character(256):: par_starname
    integer:: get_par_starname
    if (StarInGenec%initialised) then
        par_starname = StarInGenec%starname
    else
        par_starname = InitialGenecStar%starname
    endif
    get_par_starname = 0
end function

function set_par_starname(par_starname)
    use inputparam, only: starname
    implicit none
    character(256):: par_starname
    integer:: set_par_starname
    InitialGenecStar%starname = par_starname
    set_par_starname = 0
end function

! **** End Parameters

function commit_parameters()
    implicit none
    integer:: commit_parameters
    commit_parameters = 0
end function

function commit_particles()
    use makeini, only: make_initial_star
    use genec, only: initialise_star
    use WriteSaveClose, only: OpenAll
    use inputparam, only: nzmod
    implicit none
    integer:: commit_particles
    ! makeini will actually override some things from set_defaults now! FIXME
    call make_initial_star()
    nzmod = 1
    !write(*,*) 'makeini done'
    !call OpenAll()
    !write(*,*) 'OpenAll done'
    call initialise_star()
    !write(*,*) 'initialise_star done'
    commit_particles = 0
end function

function delete_star(index_of_the_star)
    implicit none
    integer:: index_of_the_star
    integer:: delete_star
    delete_star = -1  ! not supported
end function

function evolve_model(end_time)
    use timestep, only: alter
    implicit none
    real(kindreal):: end_time
    integer:: evolve_model
    !stopping_condition = ""
    type(genec_star) :: PreviousStarInGenec
    PreviousStarInGenec = StarInGenec
    do while (alter < end_time)
      !if (stopping_condition == "") then
        write(*,*) "Current time: ", alter, ", evolving to: ", end_time
        evolve_model = evolve_one_step(0)
    end do
    evolve_model = 0
end function

function evolve_for(index_of_the_star, time)
    ! get current time
    ! set max time to current time plus argument
    ! evolve
    use timestep, only: alter
    implicit none
    integer:: index_of_the_star
    real(kindreal):: time, end_time
    integer:: evolve_for
    end_time = alter+time
    do while (alter < end_time)
        write(*,*) "Current time: ", alter, ", evolving to: ", end_time
        evolve_for = evolve_one_step(index_of_the_star)
    end do

    evolve_for = 0
end function

function evolve_one_step(index_of_the_star)
    use timestep, only: alter
    use WriteSaveClose, only: OpenAll
    use genec, only: evolve, modell, finalise, veryFirst
    use inputparam,only: modanf,nwseq,nzmod,end_at_phase,end_at_model
    use genec, only: n_snap
    implicit none
    integer:: index_of_the_star
    integer:: evolve_one_step
    integer:: original_nzmod
    !original_nzmod = nzmod
    nzmod = 1
    n_snap = 0
    end_at_phase=4
    end_at_model=0

    !write(*,*) "Evolving one step, current time: ", alter
    !if (stopping_condition == "") then
      call evolve()
      !if (stopping_condition /= "") return
      call finalise()
      !call OpenAll()
      !call initialise_star() ! will set modell to 1
      !write(*,*) "Evolved one step, current time: ", alter
      !nzmod = original_nzmod
      !write(*,*) "*****modanf, nwseq, nzmod: ", modanf, nwseq, nzmod
    !else
    !  write(*,*) "stopped: ", stopping_condition
    !endif
    veryFirst = .false.
    evolve_one_step = 0
end function

function finalize_stellar_model()
    implicit none
    integer:: finalize_stellar_model
    finalize_stellar_model = 0
end function

function write_genec_model()
    !use inputparam, only: modanf
    !use WriteSaveClose, only: OpenAll
    !use genec, only: finalise
    !use helpers, only: initialise_star
    implicit none
    integer:: write_genec_model
    !call finalise()
    ! modanf = 0
    !call OpenAll()
    !call initialise_star()
    !call finalise()
    !call OpenAll()
    !call initialise_star()
    write_genec_model = -1
end function

function get_age(index_of_the_star, age)
    use timestep, only: alter
    implicit none
    integer:: index_of_the_star
    real(kindreal):: age
    integer:: get_age
    age = alter
    get_age = 0
end function

function set_age(index_of_the_star, age)
    use timestep, only: alter
    implicit none
    integer:: index_of_the_star
    real(kindreal):: age
    integer:: set_age
    alter = age
    set_age = 0
end function

function get_surface_velocity(index_of_the_star, surface_velocity)
    use genec, only: vequat
    implicit none
    integer:: index_of_the_star
    real(kindreal):: surface_velocity
    integer:: get_surface_velocity
    surface_velocity = vequat
    get_surface_velocity = 0
end function

function get_omegi_at_zone(index_of_the_star, zone, omegi_i)
    use rotmod, only: omegi
    use strucmod, only: m
    implicit none
    integer:: index_of_the_star
    integer:: get_omegi_at_zone
    real(kindreal):: omegi_i
    integer:: zone, i
    i = m - zone
    if (zone <= m) then
        omegi_i = omegi(i)
    end if
    get_omegi_at_zone = 0
end function

function get_density_at_zone(index_of_the_star, zone, rho_i)
    use strucmod, only: rho, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: rho_i
    integer:: get_density_at_zone
    i = m - zone
    if (zone <= m) then
        rho_i = exp(rho(i))
    end if
    get_density_at_zone = 0
end function

function set_density_at_zone(index_of_the_star, zone, rho_i)
    use strucmod, only: rho, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: rho_i
    integer:: set_density_at_zone
    i = m - zone
    if (zone <= m) then
        rho(i) = log(rho_i)
    end if
    set_density_at_zone = 0
end function

function get_luminosity(index_of_the_star, luminosity)
    !use strucmod, only: s, m
    use caramodele, only: gls
    implicit none
    integer:: index_of_the_star
    real(kindreal):: luminosity
    integer:: get_luminosity
    !luminosity = exp(s(m))  ! in cgs units, so erg/s?
    luminosity = gls
    get_luminosity = 0
end function

function set_luminosity(index_of_the_star, luminosity)
    !use strucmod, only: s, m
    use caramodele, only: gls
    implicit none
    integer:: index_of_the_star
    real(kindreal):: luminosity
    integer:: set_luminosity
    !luminosity = exp(s(m))  ! in cgs units, so erg/s?
    gls = luminosity
    set_luminosity = 0
end function

function get_luminosity_at_zone(index_of_the_star, zone, lum_i)
    use strucmod, only: s, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: lum_i
    integer:: get_luminosity_at_zone
    i = m - zone
    if (zone <= m) then
        !lum_i = exp(s(zone+1))
        lum_i = exp(s(i)) - 1
    end if
    get_luminosity_at_zone = 0
end function

function set_luminosity_at_zone(index_of_the_star, zone, lum_i)
    use strucmod, only: s, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: lum_i
    integer:: set_luminosity_at_zone
    i = m - zone
    if (zone <= m) then
        s(i) = log(lum_i + 1)
    end if
    set_luminosity_at_zone = 0
end function

function get_mass_fraction_at_zone(index_of_the_star, zone, dq_i)
    use strucmod, only: q, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: dq_i
    integer:: get_mass_fraction_at_zone
    i = m - zone
    if (i == 1) then
        dq_i = exp(q(i))
    else if (i <= m) then
        !dq_i = 1-exp(q(zone+1))
        dq_i = (exp(q(i)) - exp(q(i-1)))
    end if
    get_mass_fraction_at_zone = 0
end function

function set_mass_fraction_at_zone(index_of_the_star, zone, dq_i)
    use strucmod, only: q, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: dq_i
    integer:: set_mass_fraction_at_zone
    i = m - zone
    if (i == 1) then
        q(i) = log(dq_i)
    else if (i <= m) then
        !dq_i = 1-exp(q(zone+1))
        q(i) = log(exp(q(i-1)) + dq_i)  ! this won't do
    end if
    set_mass_fraction_at_zone = -1 ! This function is incomplete!
end function

function get_mass(index_of_the_star, mass)
    use caramodele, only: gms
    implicit none
    real(kindreal):: mass
    integer:: get_mass, index_of_the_star
    mass = gms
    get_mass = 0
end function

function set_mass(index_of_the_star, mass)
    implicit none
    integer:: set_mass, index_of_the_star
    real(kindreal):: mass
    InitialGenecStar%mstar = mass
    set_mass = 0
end function

function get_mass_of_species(index_of_the_star, species, species_mass)
    implicit none
    integer:: index_of_the_star
    integer:: species
    real(kindreal):: species_mass
    integer:: get_mass_of_species
    get_mass_of_species = 0
    select case(species)
    case(1)
        species_mass = 1.
    case(2)
        species_mass = 3.
    case(3)
        species_mass = 4.
    case(4)
        species_mass = 12.
    case(5)
        species_mass = 13.
    case(6)
        species_mass = 14.
    case(7)
        species_mass = 15.
    case(8)
        species_mass = 16.
    case(9)
        species_mass = 17.
    case(10)
        species_mass = 18.
    case(11)
        species_mass = 20.
    case(12)
        species_mass = 22.
    case(13)
        species_mass = 24.
    case(14)
        species_mass = 25.
    case(15)
        species_mass = 26.
    case(16)
        species_mass = 14.
    case(17)
        species_mass = 18.
    case(18)
        species_mass = 19.
    case(19)
        species_mass = 21.
    case(20)
        species_mass = 23.
    case(21)
        species_mass = 26.
    case(22)
        species_mass = 27.
    case(23)
        species_mass = 28.
    case(24)
        species_mass = 1. !neut(i)
    case(25)
        species_mass = 1. !prot(i)
!    case(26)
!        species_mass = bid(i)
!    case(27)
!        species_mass = bid1(i)
    case default
        species_mass = 0.
    end select
end function

function get_mass_fraction_of_species_at_zone(index_of_the_star, species, zone, Xj_i)
    use strucmod, only: m
    use abundmod, only: &
        x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24,&
        xmg25,xmg26,xal26,xal27,xsi28,xprot,xneut,xbid,xbid1
    implicit none
    integer:: index_of_the_star
    integer:: species, zone, i
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_species_at_zone
    i = m-zone
    if (zone <= m) then
      select case(species)
      case(1)
          Xj_i = x(i)
      case(2)
          Xj_i = y3(i)
      case(3)
          Xj_i = y(i)
      case(4)
          Xj_i = xc12(i)
      case(5)
          Xj_i = xc13(i)
      case(6)
          Xj_i = xn14(i)
      case(7)
          Xj_i = xn15(i)
      case(8)
          Xj_i = xo16(i)
      case(9)
          Xj_i = xo17(i)
      case(10)
          Xj_i = xo18(i)
      case(11)
          Xj_i = xne20(i)
      case(12)
          Xj_i = xne22(i)
      case(13)
          Xj_i = xmg24(i)
      case(14)
          Xj_i = xmg25(i)
      case(15)
          Xj_i = xmg26(i)
      case(16)
          Xj_i = xc14(i)
      case(17)
          Xj_i = xf18(i)
      case(18)
          Xj_i = xf19(i)
      case(19)
          Xj_i = xne21(i)
      case(20)
          Xj_i = xna23(i)
      case(21)
          Xj_i = xal26(i)
      case(22)
          Xj_i = xal27(i)
      case(23)
          Xj_i = xsi28(i)
      case(24)
          Xj_i = xneut(i)
      case(25)
          Xj_i = xprot(i)
      case(26)
          Xj_i = xbid(i)
      case(27)
          Xj_i = xbid1(i)
      case default
          Xj_i = 0
      end select
    end if

    get_mass_fraction_of_species_at_zone = 0
end function

function set_mass_fraction_of_species_at_zone(index_of_the_star, species, zone, Xj_i)
    use strucmod, only: m
    use abundmod, only: &
        x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24,&
        xmg25,xmg26,xal26,xal27,xsi28,xprot,xneut,xbid,xbid1
    implicit none
    integer:: index_of_the_star
    integer:: species, zone, i
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_species_at_zone
    i = m - zone
    if (zone <= m) then
        select case(species)
      case(1)
          x(i) = Xj_i
      case(2)
          y3(i) = Xj_i
      case(3)
          y(i) = Xj_i
      case(4)
          xc12(i) = Xj_i
      case(5)
          xc13(i) = Xj_i
      case(6)
          xn14(i) = Xj_i
      case(7)
          xn14(i) = Xj_i
      case(8)
          xo16(i) = Xj_i
      case(9)
          xo17(i) = Xj_i
      case(10)
          xo18(i) = Xj_i
      case(11)
          xne20(i) = Xj_i
      case(12)
          xne22(i) = Xj_i
      case(13)
          xmg24(i) = Xj_i
      case(14)
          xmg25(i) = Xj_i
      case(15)
          xmg26(i) = Xj_i
      case(16)
          xc14(i) = Xj_i
      case(17)
          xf18(i) = Xj_i
      case(18)
          xf19(i) = Xj_i
      case(19)
          xne21(i) = Xj_i
      case(20)
          xna23(i) = Xj_i
      case(21)
          xal26(i) = Xj_i
      case(22)
          xal27(i) = Xj_i
      case(23)
          xsi28(i) = Xj_i
      case(24)
          xneut(i) = Xj_i
      case(25)
          xprot(i) = Xj_i
      case(26)
          xbid(i) = Xj_i
      case(27)
          xbid1(i) = Xj_i
      !case default
      end select
    end if

    set_mass_fraction_of_species_at_zone = 0
end function

function get_metallicity(metallicity)
    implicit none
    real(kindreal):: metallicity
    integer:: get_metallicity
    metallicity = InitialGenecStar%zini
    get_metallicity = 0
end function

function set_metallicity(metallicity)
    implicit none
    real(kindreal):: metallicity
    integer:: set_metallicity
    InitialGenecStar%zini = metallicity
    set_metallicity = 0
end function

function get_mu_at_zone(index_of_the_star, zone, mu_i)
    use strucmod, only: m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: mu_i, X, Y3, Y
    integer:: err
    integer:: get_mu_at_zone
    i = (m - zone)
    if (zone <= m) then
        mu_i = 0.0d0
        err = get_mass_fraction_of_species_at_zone(index_of_the_star, 1, zone, X)
        err = get_mass_fraction_of_species_at_zone(index_of_the_star, 2, zone, Y3)
        err = get_mass_fraction_of_species_at_zone(index_of_the_star, 3, zone, Y)

        mu_i = (4.0 / (2.0 + 6 * X - (Y+Y3)))
    end if
    get_mu_at_zone = 0
end function

function get_name_of_species(index_of_the_star, species, species_name)
    implicit none
    integer:: index_of_the_star
    integer:: species
    character(len = 6):: species_name
    integer:: get_name_of_species

    character(len = 6), dimension(27):: species_names

    species_names(1)  = 'h'
    species_names(2)  = 'he3'
    species_names(3)  = 'he'
    species_names(4)  = 'c12'
    species_names(5)  = 'c13'
    species_names(6)  = 'n14'
    species_names(7)  = 'n15'
    species_names(8)  = 'o16'
    species_names(9)  = 'o17'
    species_names(10) = 'o18'
    species_names(11) = 'ne20'
    species_names(12) = 'ne22'
    species_names(13) = 'mg24'
    species_names(14) = 'mg25'
    species_names(15) = 'mg26'
    species_names(16) = 'c14'
    species_names(17) = 'f18'
    species_names(18) = 'f19'
    species_names(19) = 'ne21'
    species_names(20) = 'na23'
    species_names(21) = 'al26'
    species_names(22) = 'al27'
    species_names(23) = 'si28'
    species_names(24) = 'neut'
    species_names(25) = 'prot'
    species_names(26) = 'bid'
    species_names(27) = 'bid1'
    species_name = species_names(species)
    !x: H1
    !y: He4
    !y3: He3
    !xneut: neutron
    !xprot: proton
    !xc12: C12
    !xXYY: XYY - X element YY mass number
    !xXXYY: as above
    get_name_of_species = 0
end function

function get_mass_fraction_of_h_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_h_at_zone
    get_mass_fraction_of_h_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 1, zone, Xj_i)
end function


function get_mass_fraction_of_he3_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_he3_at_zone
    get_mass_fraction_of_he3_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 2, zone, Xj_i)
end function


function get_mass_fraction_of_he_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_he_at_zone
    get_mass_fraction_of_he_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 3, zone, Xj_i)
end function


function get_mass_fraction_of_c12_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_c12_at_zone
    get_mass_fraction_of_c12_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 4, zone, Xj_i)
end function


function get_mass_fraction_of_c13_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_c13_at_zone
    get_mass_fraction_of_c13_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 5, zone, Xj_i)
end function


function get_mass_fraction_of_n14_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_n14_at_zone
    get_mass_fraction_of_n14_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 6, zone, Xj_i)
end function


function get_mass_fraction_of_n15_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_n15_at_zone
    get_mass_fraction_of_n15_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 7, zone, Xj_i)
end function


function get_mass_fraction_of_o16_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_o16_at_zone
    get_mass_fraction_of_o16_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 8, zone, Xj_i)
end function


function get_mass_fraction_of_o17_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_o17_at_zone
    get_mass_fraction_of_o17_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 9, zone, Xj_i)
end function


function get_mass_fraction_of_o18_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_o18_at_zone
    get_mass_fraction_of_o18_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 10, zone, Xj_i)
end function


function get_mass_fraction_of_ne20_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_ne20_at_zone
    get_mass_fraction_of_ne20_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 11, zone, Xj_i)
end function


function get_mass_fraction_of_ne22_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_ne22_at_zone
    get_mass_fraction_of_ne22_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 12, zone, Xj_i)
end function


function get_mass_fraction_of_mg24_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_mg24_at_zone
    get_mass_fraction_of_mg24_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 13, zone, Xj_i)
end function


function get_mass_fraction_of_mg25_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_mg25_at_zone
    get_mass_fraction_of_mg25_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 14, zone, Xj_i)
end function


function get_mass_fraction_of_mg26_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_mg26_at_zone
    get_mass_fraction_of_mg26_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 15, zone, Xj_i)
end function


function get_mass_fraction_of_c14_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_c14_at_zone
    get_mass_fraction_of_c14_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 16, zone, Xj_i)
end function


function get_mass_fraction_of_f18_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_f18_at_zone
    get_mass_fraction_of_f18_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 17, zone, Xj_i)
end function


function get_mass_fraction_of_f19_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_f19_at_zone
    get_mass_fraction_of_f19_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 18, zone, Xj_i)
end function


function get_mass_fraction_of_ne21_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_ne21_at_zone
    get_mass_fraction_of_ne21_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 19, zone, Xj_i)
end function


function get_mass_fraction_of_na23_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_na23_at_zone
    get_mass_fraction_of_na23_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 20, zone, Xj_i)
end function


function get_mass_fraction_of_al26_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_al26_at_zone
    get_mass_fraction_of_al26_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 21, zone, Xj_i)
end function


function get_mass_fraction_of_al27_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_al27_at_zone
    get_mass_fraction_of_al27_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 22, zone, Xj_i)
end function


function get_mass_fraction_of_si28_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_si28_at_zone
    get_mass_fraction_of_si28_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 23, zone, Xj_i)
end function


function get_mass_fraction_of_neut_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_neut_at_zone
    get_mass_fraction_of_neut_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 24, zone, Xj_i)
end function


function get_mass_fraction_of_prot_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_prot_at_zone
    get_mass_fraction_of_prot_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 25, zone, Xj_i)
end function


function get_mass_fraction_of_bid_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_bid_at_zone
    get_mass_fraction_of_bid_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 26, zone, Xj_i)
end function


function get_mass_fraction_of_bid1_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_bid1_at_zone
    get_mass_fraction_of_bid1_at_zone = get_mass_fraction_of_species_at_zone(index_of_the_star, 27, zone, Xj_i)
end function


function set_mass_fraction_of_h_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_h_at_zone
    set_mass_fraction_of_h_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 1, zone, Xj_i)
end function


function set_mass_fraction_of_he3_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_he3_at_zone
    set_mass_fraction_of_he3_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 2, zone, Xj_i)
end function


function set_mass_fraction_of_he_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_he_at_zone
    set_mass_fraction_of_he_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 3, zone, Xj_i)
end function


function set_mass_fraction_of_c12_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_c12_at_zone
    set_mass_fraction_of_c12_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 4, zone, Xj_i)
end function


function set_mass_fraction_of_c13_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_c13_at_zone
    set_mass_fraction_of_c13_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 5, zone, Xj_i)
end function


function set_mass_fraction_of_n14_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_n14_at_zone
    set_mass_fraction_of_n14_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 6, zone, Xj_i)
end function


function set_mass_fraction_of_n15_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_n15_at_zone
    set_mass_fraction_of_n15_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 7, zone, Xj_i)
end function


function set_mass_fraction_of_o16_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_o16_at_zone
    set_mass_fraction_of_o16_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 8, zone, Xj_i)
end function


function set_mass_fraction_of_o17_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_o17_at_zone
    set_mass_fraction_of_o17_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 9, zone, Xj_i)
end function


function set_mass_fraction_of_o18_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_o18_at_zone
    set_mass_fraction_of_o18_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 10, zone, Xj_i)
end function


function set_mass_fraction_of_ne20_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_ne20_at_zone
    set_mass_fraction_of_ne20_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 11, zone, Xj_i)
end function


function set_mass_fraction_of_ne22_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_ne22_at_zone
    set_mass_fraction_of_ne22_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 12, zone, Xj_i)
end function


function set_mass_fraction_of_mg24_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_mg24_at_zone
    set_mass_fraction_of_mg24_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 13, zone, Xj_i)
end function


function set_mass_fraction_of_mg25_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_mg25_at_zone
    set_mass_fraction_of_mg25_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 14, zone, Xj_i)
end function


function set_mass_fraction_of_mg26_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_mg26_at_zone
    set_mass_fraction_of_mg26_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 15, zone, Xj_i)
end function


function set_mass_fraction_of_c14_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_c14_at_zone
    set_mass_fraction_of_c14_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 16, zone, Xj_i)
end function


function set_mass_fraction_of_f18_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_f18_at_zone
    set_mass_fraction_of_f18_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 17, zone, Xj_i)
end function


function set_mass_fraction_of_f19_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_f19_at_zone
    set_mass_fraction_of_f19_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 18, zone, Xj_i)
end function


function set_mass_fraction_of_ne21_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_ne21_at_zone
    set_mass_fraction_of_ne21_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 19, zone, Xj_i)
end function


function set_mass_fraction_of_na23_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_na23_at_zone
    set_mass_fraction_of_na23_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 20, zone, Xj_i)
end function


function set_mass_fraction_of_al26_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_al26_at_zone
    set_mass_fraction_of_al26_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 21, zone, Xj_i)
end function


function set_mass_fraction_of_al27_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_al27_at_zone
    set_mass_fraction_of_al27_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 22, zone, Xj_i)
end function


function set_mass_fraction_of_si28_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_si28_at_zone
    set_mass_fraction_of_si28_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 23, zone, Xj_i)
end function


function set_mass_fraction_of_neut_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_neut_at_zone
    set_mass_fraction_of_neut_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 24, zone, Xj_i)
end function


function set_mass_fraction_of_prot_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_prot_at_zone
    set_mass_fraction_of_prot_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 25, zone, Xj_i)
end function


function set_mass_fraction_of_bid_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_bid_at_zone
    set_mass_fraction_of_bid_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 26, zone, Xj_i)
end function


function set_mass_fraction_of_bid1_at_zone(index_of_the_star, zone, Xj_i)
    implicit none
    integer:: index_of_the_star, zone
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_bid1_at_zone
    set_mass_fraction_of_bid1_at_zone = set_mass_fraction_of_species_at_zone(index_of_the_star, 27, zone, Xj_i)
end function

function get_number_of_particles(n)
    implicit none
    integer:: n
    integer:: get_number_of_particles
    n = 1
    get_number_of_particles = 0
end function

function get_number_of_species(index_of_the_star, n_species)
    use inputparam, only: ialflu
    implicit none
    integer:: index_of_the_star
    integer:: n_species
    integer:: get_number_of_species
    !if (ialflu==1) then
        n_species = 27
    !else
    !    n_species = 15
    !end if
    get_number_of_species = 0
end function

function get_firstlast_species_number(first, last)
    use inputparam, only: ialflu
    implicit none
    !integer:: index_of_the_star
    integer:: first, last
    integer:: get_firstlast_species_number
    first = 1
    if (ialflu==1) then
        last = 27
    else
        last = 15
    end if
    get_firstlast_species_number = 0
end function

function get_number_of_zones(index_of_the_star, n_zones)
    use strucmod, only: m
    implicit none
    integer:: index_of_the_star
    integer:: n_zones
    integer:: get_number_of_zones
    n_zones = m
    get_number_of_zones = 0
end function

function set_number_of_zones(index_of_the_star, n_zones)
    use strucmod, only: m
    implicit none
    integer:: index_of_the_star
    integer:: n_zones
    integer:: set_number_of_zones
    m = n_zones
    set_number_of_zones = 0
end function

function get_firstlast_zone(first, last)
    use strucmod, only: m
    implicit none
    !integer:: index_of_the_star
    integer:: first, last
    integer:: get_firstlast_zone
    first = 0
    last = m-1
    get_firstlast_zone = 0
end function

function get_pressure_at_zone(index_of_the_star, zone, P_i)
    use strucmod, only: p, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: P_i
    integer:: get_pressure_at_zone
    if (zone <= m) then
        i = m - zone
        P_i = exp(p(i))
    end if
    get_pressure_at_zone = 0
end function

!function set_pressure_at_zone(index_of_the_star, zone, P_i)
!    use strucmod, only: p, m
!    implicit none
!    integer:: index_of_the_star
!    integer:: zone, i
!    real(kindreal):: P_i
!    integer:: set_pressure_at_zone
!    if (zone <= m) then
!        i = m - zone
!        p(i) = log(P_i)
!    end if
!    set_pressure_at_zone = 0
!end function

function get_radius(index_of_the_star, am_radius)
    use caramodele, only: radius
    !use strucmod, only: r, m
    implicit none
    integer:: index_of_the_star
    real(kindreal):: am_radius
    integer:: get_radius
    am_radius = 10**radius
    !radius = exp(r(1))  ! in cm
    get_radius = 0
end function

function set_radius(index_of_the_star, am_radius)
    use caramodele, only: radius
    !use strucmod, only: r, m
    implicit none
    integer:: index_of_the_star
    real(kindreal):: am_radius
    integer:: set_radius
    radius = log10(am_radius)
    set_radius = 0
end function

function get_radius_at_zone(index_of_the_star, zone, R_i)
    use strucmod, only: r, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: R_i
    integer:: get_radius_at_zone
    i = m - zone
    if (zone <= m) then
        R_i = exp(r(i))  ! in cm
    end if
    get_radius_at_zone = 0
end function

function set_radius_at_zone(index_of_the_star, zone, R_i)
    use strucmod, only: r, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: R_i
    integer:: set_radius_at_zone
    i = m - zone
    if (zone <= m) then
        r(i) = log(R_i)
    end if
    set_radius_at_zone = 0
end function

function get_stellar_type(index_of_the_star, stellar_type)
    implicit none
    integer:: index_of_the_star
    integer:: stellar_type
    integer:: get_stellar_type
    get_stellar_type = -1
end function

function get_temperature(index_of_the_star, temperature)
    use caramodele, only: teff
    implicit none
    integer:: index_of_the_star
    real(kindreal):: temperature
    integer:: get_temperature
    temperature = teff
    get_temperature = 0
end function

function get_temperature_at_zone(index_of_the_star, zone, T_i)
    use strucmod, only: t, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: T_i
    integer:: get_temperature_at_zone
    i = m - zone
    if (zone <= m) then
        T_i = exp(t(i))
    else
        get_temperature_at_zone = -2
        return
    end if
    get_temperature_at_zone = 0
end function

function get_time_step(index_of_the_star, time_step)
    use timestep, only: dzeitj
    implicit none
    integer:: index_of_the_star
    real(kindreal):: time_step
    integer:: get_time_step
    time_step = dzeitj
    get_time_step = 0
end function

function get_time(time)
    use timestep, only: alter
    implicit none
    real(kindreal):: time
    integer:: get_time
    time = alter
    get_time = 0
end function

function star_to_genec()
    ! copy values from StarInGenec to Genec
end function

function star_from_genec()
    ! copy values to StarInGenec from Genec
end function

function new_particle(index_of_the_star, mass, metallicity, am_starname)
    implicit none
    integer:: index_of_the_star, key
    real(kindreal):: mass, metallicity
    integer:: new_particle
    character(len=12):: am_starname
    number_of_stars = number_of_stars + 1
    InitialGenecStar%starname = am_starname
    InitialGenecStar%index_of_the_star = number_of_stars
    InitialGenecStar%mstar = mass
    InitialGenecStar%zini = metallicity
    InitialGenecStar%idefaut = 1
    
    new_particle = 0
end function

function recommit_parameters()
    implicit none
    integer:: recommit_parameters
    recommit_parameters = 0
end function

function recommit_particles()
    use genec, only: finalise, initialise_star
    use WriteSaveClose, only: OpenAll
    implicit none
    integer:: recommit_particles
    !call finalise()
    !call OpenAll()
    !call initialise_star()
    recommit_particles = 0
end function

function set_genec_path(path)
    use evol, only: input_dir
    implicit none
    character(len=200):: path
    integer:: set_genec_path

    input_dir = path
    set_genec_path = 0
end function

function set_starname(index_of_the_star)
    implicit none
    integer:: set_starname, index_of_the_star
    InitialGenecStar%starname = 'AmuseStar'
    set_starname = 0
end function

function set_temperature_at_zone(index_of_the_star, zone, T_i)
    use strucmod, only: t, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: T_i
    integer:: set_temperature_at_zone
    i = m - zone
    if (zone <= m) then
        t(i) = log(T_i)
    end if
    
    set_temperature_at_zone = 0
end function

end module AmuseInterface
