module AmuseInterface
    use storage, only:&
            InitialGenecStar,InitialNetwork,&
            GenecStar,genec_star
    use helpers, only:&
            copy_to_genec_star,copy_from_genec_star
    use evol, only: kindreal,ldi,npondcouche

    type(genec_star) :: BackupBackupGenecStar
    type(genec_star) :: BackupGenecStar
    integer :: number_of_stars = 0
    integer :: evolve_steps
    public:: BackupGenecStar,number_of_stars
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
    use io_definitions
    implicit none
    integer:: initialize_code

    libgenec = .true.
    !io_logs = 6
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
    implicit none
    integer:: model_number
    integer:: get_model_number
    model_number = GenecStar%nwmd
    get_model_number = 0
end function

function set_model_number(model_number)
    implicit none
    integer:: model_number
    integer:: set_model_number
    GenecStar%nwmd = model_number
    set_model_number = 0
end function

function get_min_timestep_stop_condition(min_timestep_stop_condition)
    implicit none
    real(kindreal):: min_timestep_stop_condition
    integer:: get_min_timestep_stop_condition
    get_min_timestep_stop_condition = 0
end function

function get_par_n_snap(n_snap_out)
    implicit none
    integer:: n_snap_out
    integer:: get_par_n_snap
    n_snap_out = GenecStar%n_snap
    get_par_n_snap = 0
end function

function set_par_n_snap(n_snap_in)
    implicit none
    integer:: n_snap_in
    integer:: set_par_n_snap
    GenecStar%n_snap = n_snap_in
    set_par_n_snap= 0
end function

function get_par_ipoly(par_ipoly)
    implicit none
    integer:: par_ipoly
    integer:: get_par_ipoly
    par_ipoly = InitialGenecStar%ipoly
    get_par_ipoly = 0
end function

function set_par_ipoly(par_ipoly)
    implicit none
    integer:: par_ipoly
    integer:: set_par_ipoly
    InitialGenecStar%ipoly = par_ipoly
    set_par_ipoly = 0
end function

function get_par_nwseq(par_nwseq)
    implicit none
    integer:: par_nwseq
    integer:: get_par_nwseq
    par_nwseq = GenecStar%nwseq
    get_par_nwseq = 0
end function

function set_par_nwseq(par_nwseq)
    implicit none
    integer:: par_nwseq
    integer:: set_par_nwseq
    GenecStar%nwseq = par_nwseq
    set_par_nwseq = 0
end function

function get_par_modanf(par_modanf)
    implicit none
    integer:: par_modanf
    integer:: get_par_modanf
    par_modanf = GenecStar%modanf
    get_par_modanf = 0
end function

function set_par_modanf(par_modanf)
    implicit none
    integer:: par_modanf
    integer:: set_par_modanf
    GenecStar%modanf = par_modanf
    set_par_modanf = 0
end function

function get_par_nzmod(par_nzmod)
    implicit none
    integer:: par_nzmod
    integer:: get_par_nzmod
    par_nzmod = GenecStar%nzmod
    get_par_nzmod = 0
end function

function set_par_nzmod(par_nzmod)
    implicit none
    integer:: par_nzmod
    integer:: set_par_nzmod
    GenecStar%nzmod = par_nzmod
    set_par_nzmod = 0
end function

function get_par_irot(par_irot)
    implicit none
    integer:: par_irot
    integer:: get_par_irot
    par_irot = GenecStar%irot
    get_par_irot = 0
end function

function set_par_irot(par_irot)
    implicit none
    integer:: par_irot
    integer:: set_par_irot
    GenecStar%irot = par_irot
    set_par_irot = 0
end function

function get_par_isol(par_isol)
    implicit none
    integer:: par_isol
    integer:: get_par_isol
    par_isol = GenecStar%isol
    get_par_isol = 0
end function

function set_par_isol(par_isol)
    implicit none
    integer:: par_isol
    integer:: set_par_isol
    GenecStar%isol = par_isol
    set_par_isol = 0
end function

function get_par_imagn(par_imagn)
    implicit none
    integer:: par_imagn
    integer:: get_par_imagn
    par_imagn = GenecStar%imagn
    get_par_imagn = 0
end function

function set_par_imagn(par_imagn)
    implicit none
    integer:: par_imagn
    integer:: set_par_imagn
    GenecStar%imagn = par_imagn
    set_par_imagn = 0
end function

function get_par_ialflu(par_ialflu)
    implicit none
    integer:: par_ialflu
    integer:: get_par_ialflu
    par_ialflu = GenecStar%ialflu
    get_par_ialflu = 0
end function

function set_par_ialflu(par_ialflu)
    implicit none
    integer:: par_ialflu
    integer:: set_par_ialflu
    GenecStar%ialflu = par_ialflu
    set_par_ialflu = 0
end function

function get_par_ianiso(par_ianiso)
    implicit none
    integer:: par_ianiso
    integer:: get_par_ianiso
    par_ianiso = GenecStar%ianiso
    get_par_ianiso = 0
end function

function set_par_ianiso(par_ianiso)
    implicit none
    integer:: par_ianiso
    integer:: set_par_ianiso
    GenecStar%ianiso = par_ianiso
    set_par_ianiso = 0
end function

function get_par_ipop3(par_ipop3)
    implicit none
    integer:: par_ipop3
    integer:: get_par_ipop3
    par_ipop3 = GenecStar%ipop3
    get_par_ipop3 = 0
end function

function set_par_ipop3(par_ipop3)
    implicit none
    integer:: par_ipop3
    integer:: set_par_ipop3
    GenecStar%ipop3 = par_ipop3
    set_par_ipop3 = 0
end function

function get_par_ibasnet(par_ibasnet)
    implicit none
    integer:: par_ibasnet
    integer:: get_par_ibasnet
    par_ibasnet = GenecStar%ibasnet
    get_par_ibasnet = 0
end function

function set_par_ibasnet(par_ibasnet)
    implicit none
    integer:: par_ibasnet
    integer:: set_par_ibasnet
    GenecStar%ibasnet = par_ibasnet
    set_par_ibasnet = 0
end function

function get_par_iopac(par_iopac)
    implicit none
    integer:: par_iopac
    integer:: get_par_iopac
    par_iopac = GenecStar%iopac
    get_par_iopac = 0
end function

function set_par_iopac(par_iopac)
    implicit none
    integer:: par_iopac
    integer:: set_par_iopac
    GenecStar%iopac = par_iopac
    set_par_iopac = 0
end function

function get_par_ikappa(par_ikappa)
    implicit none
    integer:: par_ikappa
    integer:: get_par_ikappa
    par_ikappa = GenecStar%ikappa
    get_par_ikappa = 0
end function

function set_par_ikappa(par_ikappa)
    implicit none
    integer:: par_ikappa
    integer:: set_par_ikappa
    GenecStar%ikappa = par_ikappa
    set_par_ikappa = 0
end function

function get_par_idiff(par_idiff)
    implicit none
    integer:: par_idiff
    integer:: get_par_idiff
    par_idiff = GenecStar%idiff
    get_par_idiff = 0
end function

function set_par_idiff(par_idiff)
    implicit none
    integer:: par_idiff
    integer:: set_par_idiff
    GenecStar%idiff = par_idiff
    set_par_idiff = 0
end function

function get_par_iadvec(par_iadvec)
    implicit none
    integer:: par_iadvec
    integer:: get_par_iadvec
    par_iadvec = GenecStar%iadvec
    get_par_iadvec = 0
end function

function set_par_iadvec(par_iadvec)
    implicit none
    integer:: par_iadvec
    integer:: set_par_iadvec
    GenecStar%iadvec = par_iadvec
    set_par_iadvec = 0
end function

function get_par_istati(par_istati)
    implicit none
    integer:: par_istati
    integer:: get_par_istati
    par_istati = GenecStar%istati
    get_par_istati = 0
end function

function set_par_istati(par_istati)
    implicit none
    integer:: par_istati
    integer:: set_par_istati
    GenecStar%istati = par_istati
    set_par_istati = 0
end function

function get_par_icoeff(par_icoeff)
    implicit none
    integer:: par_icoeff
    integer:: get_par_icoeff
    par_icoeff = GenecStar%icoeff
    get_par_icoeff = 0
end function

function set_par_icoeff(par_icoeff)
    implicit none
    integer:: par_icoeff
    integer:: set_par_icoeff
    GenecStar%icoeff = par_icoeff
    set_par_icoeff = 0
end function

function get_par_igamma(par_igamma)
    implicit none
    integer:: par_igamma
    integer:: get_par_igamma
    par_igamma = GenecStar%igamma
    get_par_igamma = 0
end function

function set_par_igamma(par_igamma)
    implicit none
    integer:: par_igamma
    integer:: set_par_igamma
    GenecStar%igamma = par_igamma
    set_par_igamma = 0
end function

function get_par_idialo(par_idialo)
    implicit none
    integer:: par_idialo
    integer:: get_par_idialo
    par_idialo = GenecStar%idialo
    get_par_idialo = 0
end function

function set_par_idialo(par_idialo)
    implicit none
    integer:: par_idialo
    integer:: set_par_idialo
    GenecStar%idialo = par_idialo
    set_par_idialo = 0
end function

function get_par_idialu(par_idialu)
    implicit none
    integer:: par_idialu
    integer:: get_par_idialu
    par_idialu = GenecStar%idialu
    get_par_idialu = 0
end function

function set_par_idialu(par_idialu)
    implicit none
    integer:: par_idialu
    integer:: set_par_idialu
    GenecStar%idialu = par_idialu
    set_par_idialu = 0
end function

function get_par_imloss(par_imloss)
    implicit none
    integer:: par_imloss
    integer:: get_par_imloss
    par_imloss = GenecStar%imloss
    get_par_imloss = 0
end function

function set_par_imloss(par_imloss)
    implicit none
    integer:: par_imloss
    integer:: set_par_imloss
    GenecStar%imloss = par_imloss
    set_par_imloss = 0
end function

function get_par_ifitm(par_ifitm)
    implicit none
    integer:: par_ifitm
    integer:: get_par_ifitm
    par_ifitm = GenecStar%ifitm
    get_par_ifitm = 0
end function

function set_par_ifitm(par_ifitm)
    implicit none
    integer:: par_ifitm
    integer:: set_par_ifitm
    GenecStar%ifitm = par_ifitm
    set_par_ifitm = 0
end function

function get_par_nndr(par_nndr)
    implicit none
    integer:: par_nndr
    integer:: get_par_nndr
    par_nndr = GenecStar%nndr
    get_par_nndr = 0
end function

function set_par_nndr(par_nndr)
    implicit none
    integer:: par_nndr
    integer:: set_par_nndr
    GenecStar%nndr = par_nndr
    set_par_nndr = 0
end function

function get_par_iledou(par_iledou)
    implicit none
    integer:: par_iledou
    integer:: get_par_iledou
    par_iledou = GenecStar%iledou
    get_par_iledou = 0
end function

function set_par_iledou(par_iledou)
    implicit none
    integer:: par_iledou
    integer:: set_par_iledou
    GenecStar%iledou = par_iledou
    set_par_iledou = 0
end function

function get_par_idifcon(par_idifcon)
    implicit none
    integer:: par_idifcon
    integer:: get_par_idifcon
    par_idifcon = GenecStar%idifcon
    get_par_idifcon = 0
end function

function set_par_idifcon(par_idifcon)
    implicit none
    integer:: par_idifcon
    integer:: set_par_idifcon
    GenecStar%idifcon = par_idifcon
    set_par_idifcon = 0
end function

function get_par_my(par_my)
    implicit none
    integer:: par_my
    integer:: get_par_my
    par_my = GenecStar%my
    get_par_my = 0
end function

function set_par_my(par_my)
    implicit none
    integer:: par_my
    integer:: set_par_my
    GenecStar%my = par_my
    set_par_my = 0
end function

function get_par_iover(par_iover)
    implicit none
    integer:: par_iover
    integer:: get_par_iover
    par_iover = GenecStar%iover
    get_par_iover = 0
end function

function set_par_iover(par_iover)
    implicit none
    integer:: par_iover
    integer:: set_par_iover
    GenecStar%iover = par_iover
    set_par_iover = 0
end function

function get_par_iunder(par_iunder)
    implicit none
    integer:: par_iunder
    integer:: get_par_iunder
    par_iunder = GenecStar%iunder
    get_par_iunder = 0
end function

function set_par_iunder(par_iunder)
    implicit none
    integer:: par_iunder
    integer:: set_par_iunder
    GenecStar%iunder = par_iunder
    set_par_iunder = 0
end function

function get_par_nbchx(par_nbchx)
    implicit none
    integer:: par_nbchx
    integer:: get_par_nbchx
    par_nbchx = GenecStar%nbchx
    get_par_nbchx = 0
end function

function set_par_nbchx(par_nbchx)
    implicit none
    integer:: par_nbchx
    integer:: set_par_nbchx
    GenecStar%nbchx = par_nbchx
    set_par_nbchx = 0
end function

function get_par_nrband(par_nrband)
    implicit none
    integer:: par_nrband
    integer:: get_par_nrband
    par_nrband = GenecStar%nrband
    get_par_nrband = 0
end function

function set_par_nrband(par_nrband)
    implicit none
    integer:: par_nrband
    integer:: set_par_nrband
    GenecStar%nrband = par_nrband
    set_par_nrband = 0
end function

function get_par_islow(par_islow)
    implicit none
    integer:: par_islow
    integer:: get_par_islow
    par_islow = GenecStar%islow
    get_par_islow = 0
end function

function set_par_islow(par_islow)
    implicit none
    integer:: par_islow
    integer:: set_par_islow
    GenecStar%islow = par_islow
    set_par_islow = 0
end function

function get_par_icncst(par_icncst)
    implicit none
    integer:: par_icncst
    integer:: get_par_icncst
    par_icncst = GenecStar%icncst
    get_par_icncst = 0
end function

function set_par_icncst(par_icncst)
    implicit none
    integer:: par_icncst
    integer:: set_par_icncst
    GenecStar%icncst = par_icncst
    set_par_icncst = 0
end function

function get_par_tauH_fit(par_tauH_fit)
    implicit none
    integer:: par_tauH_fit
    integer:: get_par_tauH_fit
    par_tauH_fit = GenecStar%tauH_fit
    get_par_tauH_fit = 0
end function

function set_par_tauH_fit(par_tauH_fit)
    implicit none
    integer:: par_tauH_fit
    integer:: set_par_tauH_fit
    GenecStar%tauH_fit = par_tauH_fit
    set_par_tauH_fit = 0
end function

function get_par_iauto(par_iauto)
    implicit none
    integer:: par_iauto
    integer:: get_par_iauto
    par_iauto = GenecStar%iauto
    get_par_iauto = 0
end function

function set_par_iauto(par_iauto)
    implicit none
    integer:: par_iauto
    integer:: set_par_iauto
    GenecStar%iauto = par_iauto
    set_par_iauto = 0
end function

function get_par_iprn(par_iprn)
    implicit none
    integer:: par_iprn
    integer:: get_par_iprn
    par_iprn = GenecStar%iprn
    get_par_iprn = 0
end function

function set_par_iprn(par_iprn)
    implicit none
    integer:: par_iprn
    integer:: set_par_iprn
    GenecStar%iprn = par_iprn
    set_par_iprn = 0
end function

function get_par_iout(par_iout)
    implicit none
    integer:: par_iout
    integer:: get_par_iout
    par_iout = GenecStar%iout
    get_par_iout = 0
end function

function set_par_iout(par_iout)
    implicit none
    integer:: par_iout
    integer:: set_par_iout
    GenecStar%iout = par_iout
    set_par_iout = 0
end function

function get_par_itmin(par_itmin)
    implicit none
    integer:: par_itmin
    integer:: get_par_itmin
    par_itmin = GenecStar%itmin
    get_par_itmin = 0
end function

function set_par_itmin(par_itmin)
    implicit none
    integer:: par_itmin
    integer:: set_par_itmin
    GenecStar%itmin = par_itmin
    set_par_itmin = 0
end function

function get_par_idebug(par_idebug)
    implicit none
    integer:: par_idebug
    integer:: get_par_idebug
    par_idebug = GenecStar%idebug
    get_par_idebug = 0
end function

function set_par_idebug(par_idebug)
    implicit none
    integer:: par_idebug
    integer:: set_par_idebug
    GenecStar%idebug = par_idebug
    set_par_idebug = 0
end function

function get_par_itests(par_itests)
    implicit none
    integer:: par_itests
    integer:: get_par_itests
    par_itests = GenecStar%itests
    get_par_itests = 0
end function

function set_par_itests(par_itests)
    implicit none
    integer:: par_itests
    integer:: set_par_itests
    GenecStar%itests = par_itests
    set_par_itests = 0
end function

function get_par_var_rates(par_var_rates)
    implicit none
    logical:: par_var_rates
    integer:: get_par_var_rates
    par_var_rates = GenecStar%var_rates
    get_par_var_rates = 0
end function

function set_par_var_rates(par_var_rates)
    implicit none
    logical:: par_var_rates
    integer:: set_par_var_rates
    GenecStar%var_rates = par_var_rates
    set_par_var_rates = 0
end function

function get_par_bintide(par_bintide)
    implicit none
    logical:: par_bintide
    integer:: get_par_bintide
    par_bintide = GenecStar%bintide
    get_par_bintide = 0
end function

function set_par_bintide(par_bintide)
    implicit none
    logical:: par_bintide
    integer:: set_par_bintide
    GenecStar%bintide = par_bintide
    set_par_bintide = 0
end function

function get_par_const_per(par_const_per)
    implicit none
    logical:: par_const_per
    integer:: get_par_const_per
    par_const_per = GenecStar%const_per
    get_par_const_per = 0
end function

function set_par_const_per(par_const_per)
    implicit none
    logical:: par_const_per
    integer:: set_par_const_per
    GenecStar%const_per = par_const_per
    set_par_const_per = 0
end function

function get_par_Add_Flux(par_Add_Flux)
    implicit none
    logical:: par_Add_Flux
    integer:: get_par_Add_Flux
    par_Add_Flux = GenecStar%Add_Flux
    get_par_Add_Flux = 0
end function

function set_par_Add_Flux(par_Add_Flux)
    implicit none
    logical:: par_Add_Flux
    integer:: set_par_Add_Flux
    GenecStar%Add_Flux = par_Add_Flux
    set_par_Add_Flux = 0
end function

function get_par_diff_only(par_diff_only)
    implicit none
    logical:: par_diff_only
    integer:: get_par_diff_only
    par_diff_only = GenecStar%diff_only
    get_par_diff_only = 0
end function

function set_par_diff_only(par_diff_only)
    implicit none
    logical:: par_diff_only
    integer:: set_par_diff_only
    GenecStar%diff_only = par_diff_only
    set_par_diff_only = 0
end function

function get_par_RSG_Mdot(par_RSG_Mdot)
    implicit none
    integer:: par_RSG_Mdot
    integer:: get_par_RSG_Mdot
    par_RSG_Mdot = GenecStar%RSG_Mdot
    get_par_RSG_Mdot = 0
end function

function set_par_RSG_Mdot(par_RSG_Mdot)
    implicit none
    integer:: par_RSG_Mdot
    integer:: set_par_RSG_Mdot
    GenecStar%RSG_Mdot = par_RSG_Mdot
    set_par_RSG_Mdot = 0
end function

function get_par_display_plot(par_display_plot)
    implicit none
    logical:: par_display_plot
    integer:: get_par_display_plot
    par_display_plot = GenecStar%display_plot
    get_par_display_plot = 0
end function

function set_par_display_plot(par_display_plot)
    implicit none
    logical:: par_display_plot
    integer:: set_par_display_plot
    GenecStar%display_plot = par_display_plot
    set_par_display_plot = 0
end function

function get_par_xyfiles(par_xyfiles)
    implicit none
    logical:: par_xyfiles
    integer:: get_par_xyfiles
    par_xyfiles = GenecStar%xyfiles
    get_par_xyfiles = 0
end function

function set_par_xyfiles(par_xyfiles)
    implicit none
    logical:: par_xyfiles
    integer:: set_par_xyfiles
    GenecStar%xyfiles = par_xyfiles
    set_par_xyfiles = 0
end function

function get_par_verbose(par_verbose)
    implicit none
    logical:: par_verbose
    integer:: get_par_verbose
    par_verbose = GenecStar%verbose
    get_par_verbose = 0
end function

function set_par_verbose(par_verbose)
    implicit none
    logical:: par_verbose
    integer:: set_par_verbose
    GenecStar%verbose = par_verbose
    set_par_verbose = 0
end function

function get_par_stop_deg(par_stop_deg)
    implicit none
    logical:: par_stop_deg
    integer:: get_par_stop_deg
    par_stop_deg = GenecStar%stop_deg
    get_par_stop_deg = 0
end function

function set_par_stop_deg(par_stop_deg)
    implicit none
    logical:: par_stop_deg
    integer:: set_par_stop_deg
    GenecStar%stop_deg = par_stop_deg
    set_par_stop_deg = 0
end function

function get_par_index_poly(par_index_poly)
    implicit none
    real(kindreal):: par_index_poly
    integer:: get_par_index_poly
    par_index_poly = InitialGenecStar%n
    get_par_index_poly = 0
end function

function set_par_index_poly(par_index_poly)
    implicit none
    real(kindreal):: par_index_poly
    integer:: set_par_index_poly
    InitialGenecStar%n = par_index_poly
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
    implicit none
    real(kindreal):: par_vwant
    integer:: get_par_vwant
    par_vwant = InitialGenecStar%vwant
    get_par_vwant = 0
end function

function set_par_vwant(par_vwant)
    implicit none
    real(kindreal):: par_vwant
    integer:: set_par_vwant
    InitialGenecStar%vwant = par_vwant
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
    if (GenecStar%initialised) then
        par_starname = GenecStar%starname
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

function finalize_stellar_model()
    implicit none
    integer:: finalize_stellar_model
    call copy_to_genec_star(GenecStar)
    finalize_stellar_model = 0
end function

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
    call copy_from_genec_star(GenecStar)

    !nzmod = 1
    !write(*,*) 'makeini done'
    !call OpenAll()
    !write(*,*) 'OpenAll done'
    call initialise_star()
    !call copy_network_to_star(InitialNetwork, GenecStar)
    call copy_to_genec_star(GenecStar)
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
    nzmod = 1
    n_snap = 0

    call evolve()
    call finalise()
    veryFirst = .false.
    evolve_one_step = 0
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
    implicit none
    integer:: index_of_the_star
    real(kindreal):: age
    integer:: get_age
    age = GenecStar%alter
    get_age = 0
end function

function set_age(index_of_the_star, age)
    implicit none
    integer:: index_of_the_star
    real(kindreal):: age
    integer:: set_age
    GenecStar%alter = age
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
    implicit none
    integer:: index_of_the_star
    integer:: get_omegi_at_zone
    real(kindreal):: omegi_i
    integer:: zone, i
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        omegi_i = GenecStar%omegi(i)
    end if
    get_omegi_at_zone = 0
end function

function get_density_at_zone(index_of_the_star, zone, rho_i)
    use strucmod, only: rho
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: rho_i
    integer:: get_density_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        rho_i = exp(rho(i))
    end if
    get_density_at_zone = 0
end function

function set_density_at_zone(index_of_the_star, zone, rho_i)
    use strucmod, only: rho
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: rho_i
    integer:: set_density_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        rho(i) = log(rho_i)
    end if
    set_density_at_zone = 0
end function

function get_luminosity(index_of_the_star, luminosity)
    implicit none
    integer:: index_of_the_star
    real(kindreal):: luminosity
    integer:: get_luminosity
    !luminosity = exp(s(m))  ! in cgs units, so erg/s?
    luminosity = GenecStar%gls
    get_luminosity = 0
end function

function set_luminosity(index_of_the_star, luminosity)
    implicit none
    integer:: index_of_the_star
    real(kindreal):: luminosity
    integer:: set_luminosity
    !luminosity = exp(s(m))  ! in cgs units, so erg/s?
    GenecStar%gls = luminosity
    set_luminosity = 0
end function

function get_luminosity_at_zone(index_of_the_star, zone, lum_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: lum_i
    integer:: get_luminosity_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        !lum_i = exp(s(zone+1))
        lum_i = exp(GenecStar%s(i)) - 1
    end if
    get_luminosity_at_zone = 0
end function

function set_luminosity_at_zone(index_of_the_star, zone, lum_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: lum_i
    integer:: set_luminosity_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        GenecStar%s(i) = log(lum_i + 1)
    end if
    set_luminosity_at_zone = 0
end function

function get_mass_fraction_at_zone(index_of_the_star, zone, dq_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: dq_i
    integer:: get_mass_fraction_at_zone
    i = GenecStar%m - zone
    if (i == 1) then
        dq_i = exp(GenecStar%q(i))
    else if (i <= GenecStar%m) then
        !dq_i = 1-exp(q(zone+1))
        dq_i = (exp(GenecStar%q(i)) - exp(GenecStar%q(i-1)))
    end if
    get_mass_fraction_at_zone = 0
end function

function set_mass_fraction_at_zone(index_of_the_star, zone, dq_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: dq_i
    integer:: set_mass_fraction_at_zone
    i = GenecStar%m - zone
    if (i == 1) then
        GenecStar%q(i) = log(dq_i)
    else if (i <= GenecStar%m) then
        !dq_i = 1-exp(q(zone+1))
        GenecStar%q(i) = log(exp(GenecStar%q(i-1)) + dq_i)  ! this won't do
    end if
    set_mass_fraction_at_zone = -1 ! This function is incomplete!
end function

function get_mass(index_of_the_star, mass)
    implicit none
    real(kindreal):: mass
    integer:: get_mass, index_of_the_star
    mass = GenecStar%gms
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
    implicit none
    integer:: index_of_the_star
    integer:: species, zone, i
    real(kindreal):: Xj_i
    integer:: get_mass_fraction_of_species_at_zone
    i = GenecStar%m-zone
    if (zone <= GenecStar%m) then
      select case(species)
      case(1)
          Xj_i = GenecStar%x(i)
      case(2)
          Xj_i = GenecStar%y3(i)
      case(3)
          Xj_i = GenecStar%y(i)
      case(4)
          Xj_i = GenecStar%xc12(i)
      case(5)
          Xj_i = GenecStar%xc13(i)
      case(6)
          Xj_i = GenecStar%xn14(i)
      case(7)
          Xj_i = GenecStar%xn15(i)
      case(8)
          Xj_i = GenecStar%xo16(i)
      case(9)
          Xj_i = GenecStar%xo17(i)
      case(10)
          Xj_i = GenecStar%xo18(i)
      case(11)
          Xj_i = GenecStar%xne20(i)
      case(12)
          Xj_i = GenecStar%xne22(i)
      case(13)
          Xj_i = GenecStar%xmg24(i)
      case(14)
          Xj_i = GenecStar%xmg25(i)
      case(15)
          Xj_i = GenecStar%xmg26(i)
      case(16)
          Xj_i = GenecStar%xc14(i)
      case(17)
          Xj_i = GenecStar%xf18(i)
      case(18)
          Xj_i = GenecStar%xf19(i)
      case(19)
          Xj_i = GenecStar%xne21(i)
      case(20)
          Xj_i = GenecStar%xna23(i)
      case(21)
          Xj_i = GenecStar%xal26(i)
      case(22)
          Xj_i = GenecStar%xal27(i)
      case(23)
          Xj_i = GenecStar%xsi28(i)
      case(24)
          Xj_i = GenecStar%xneut(i)
      case(25)
          Xj_i = GenecStar%xprot(i)
      case(26)
          Xj_i = GenecStar%xbid(i)
      case(27)
          Xj_i = GenecStar%xbid1(i)
      case default
          Xj_i = 0
      end select
    end if

    get_mass_fraction_of_species_at_zone = 0
end function

function set_mass_fraction_of_species_at_zone(index_of_the_star, species, zone, Xj_i)
    implicit none
    integer:: index_of_the_star
    integer:: species, zone, i
    real(kindreal):: Xj_i
    integer:: set_mass_fraction_of_species_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        select case(species)
      case(1)
          GenecStar%x(i) = Xj_i
      case(2)
          GenecStar%y3(i) = Xj_i
      case(3)
          GenecStar%y(i) = Xj_i
      case(4)
          GenecStar%xc12(i) = Xj_i
      case(5)
          GenecStar%xc13(i) = Xj_i
      case(6)
          GenecStar%xn14(i) = Xj_i
      case(7)
          GenecStar%xn14(i) = Xj_i
      case(8)
          GenecStar%xo16(i) = Xj_i
      case(9)
          GenecStar%xo17(i) = Xj_i
      case(10)
          GenecStar%xo18(i) = Xj_i
      case(11)
          GenecStar%xne20(i) = Xj_i
      case(12)
          GenecStar%xne22(i) = Xj_i
      case(13)
          GenecStar%xmg24(i) = Xj_i
      case(14)
          GenecStar%xmg25(i) = Xj_i
      case(15)
          GenecStar%xmg26(i) = Xj_i
      case(16)
          GenecStar%xc14(i) = Xj_i
      case(17)
          GenecStar%xf18(i) = Xj_i
      case(18)
          GenecStar%xf19(i) = Xj_i
      case(19)
          GenecStar%xne21(i) = Xj_i
      case(20)
          GenecStar%xna23(i) = Xj_i
      case(21)
          GenecStar%xal26(i) = Xj_i
      case(22)
          GenecStar%xal27(i) = Xj_i
      case(23)
          GenecStar%xsi28(i) = Xj_i
      case(24)
          GenecStar%xneut(i) = Xj_i
      case(25)
          GenecStar%xprot(i) = Xj_i
      case(26)
          GenecStar%xbid(i) = Xj_i
      case(27)
          GenecStar%xbid1(i) = Xj_i
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
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: mu_i, X, Y3, Y
    integer:: err
    integer:: get_mu_at_zone
    i = (GenecStar%m - zone)
    if (zone <= GenecStar%m) then
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
    implicit none
    integer:: index_of_the_star
    integer:: n_species
    integer:: get_number_of_species
    if (GenecStar%ialflu==1) then
       n_species = 27
    else
        n_species = 15
    end if
    get_number_of_species = 0
end function

function get_firstlast_species_number(first, last)
    implicit none
    !integer:: index_of_the_star
    integer:: first, last
    integer:: get_firstlast_species_number
    first = 1
    if (GenecStar%ialflu==1) then
        last = 27
    else
        last = 15
    end if
    get_firstlast_species_number = 0
end function

function get_number_of_zones(index_of_the_star, n_zones)
    implicit none
    integer:: index_of_the_star
    integer:: n_zones
    integer:: get_number_of_zones
    n_zones = GenecStar%m
    get_number_of_zones = 0
end function

function set_number_of_zones(index_of_the_star, n_zones)
    implicit none
    integer:: index_of_the_star
    integer:: n_zones
    integer:: set_number_of_zones
    GenecStar%m = n_zones
    set_number_of_zones = 0
end function

function get_firstlast_zone(first, last)
    implicit none
    !integer:: index_of_the_star
    integer:: first, last
    integer:: get_firstlast_zone
    first = 0
    last = GenecStar%m-1
    get_firstlast_zone = 0
end function

function get_pressure_at_zone(index_of_the_star, zone, P_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: P_i
    integer:: get_pressure_at_zone
    if (zone <= GenecStar%m) then
        i = GenecStar%m - zone
        P_i = exp(GenecStar%p(i))
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
    implicit none
    integer:: index_of_the_star
    real(kindreal):: am_radius
    integer:: get_radius
    am_radius = 10**GenecStar%radius
    !radius = exp(r(1))  ! in cm
    get_radius = 0
end function

function set_radius(index_of_the_star, am_radius)
    implicit none
    integer:: index_of_the_star
    real(kindreal):: am_radius
    integer:: set_radius
    GenecStar%radius = log10(am_radius)
    set_radius = 0
end function

function get_eps_at_zone(index_of_the_star, zone, eps)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: eps
    integer:: get_eps_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        eps = GenecStar%eps(i)
    end if
    get_eps_at_zone = 0
end function

function get_epsy_at_zone(index_of_the_star, zone, epsy)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: epsy
    integer:: get_epsy_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        epsy = GenecStar%epsy(i)
    end if
    get_epsy_at_zone = 0
end function

function get_eps_c_adv_at_zone(index_of_the_star, zone, eps_c_adv)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: eps_c_adv
    integer:: get_eps_c_adv_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        eps_c_adv = GenecStar%eps_c_adv(i)
    end if
    get_eps_c_adv_at_zone = 0
end function

function get_nabla_rad_at_zone(index_of_the_star, zone, nabla_rad)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: nabla_rad
    integer:: get_nabla_rad_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        nabla_rad = GenecStar%Nabla_rad(i)
    end if
    get_nabla_rad_at_zone = 0
end function

function get_nabla_ad_at_zone(index_of_the_star, zone, nabla_ad)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: nabla_ad
    integer:: get_nabla_ad_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        nabla_ad = GenecStar%Nabla_ad(i)
    end if
    get_nabla_ad_at_zone = 0
end function

function get_nabla_mu_at_zone(index_of_the_star, zone, nabla_mu)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: nabla_mu
    integer:: get_nabla_mu_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        nabla_mu = GenecStar%Nabla_mu(i)
    end if
    get_nabla_mu_at_zone = 0
end function

function get_radius_at_zone(index_of_the_star, zone, R_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: R_i
    integer:: get_radius_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        R_i = exp(GenecStar%r(i))  ! in cm
    end if
    get_radius_at_zone = 0
end function

function set_radius_at_zone(index_of_the_star, zone, R_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: R_i
    integer:: set_radius_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        GenecStar%r(i) = log(R_i)
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
    implicit none
    integer:: index_of_the_star
    real(kindreal):: temperature
    integer:: get_temperature
    temperature = GenecStar%teff
    get_temperature = 0
end function

function get_temperature_at_zone(index_of_the_star, zone, T_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: T_i
    integer:: get_temperature_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        T_i = exp(GenecStar%t(i))
    else
        get_temperature_at_zone = -2
        return
    end if
    get_temperature_at_zone = 0
end function

function get_time_step(index_of_the_star, time_step)
    implicit none
    integer:: index_of_the_star
    real(kindreal):: time_step
    integer:: get_time_step
    time_step = GenecStar%dzeitj
    get_time_step = 0
end function

function get_time(time)
    implicit none
    real(kindreal):: time
    integer:: get_time
    time = GenecStar%alter
    get_time = 0
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

function get_starname(index_of_the_star, am_starname)
    implicit none
    integer:: get_starname, index_of_the_star
    character(len=200):: am_starname
    am_starname = InitialGenecStar%starname
    get_starname = 0
end function

function set_starname(index_of_the_star, am_starname)
    implicit none
    integer:: set_starname, index_of_the_star
    character(len=200):: am_starname
    InitialGenecStar%starname = am_starname
    set_starname = 0
end function

function get_phase(index_of_the_star, phase)
    implicit none
    integer:: index_of_the_star, phase
    integer:: get_phase
    phase = GenecStar%phase
    get_phase = 0
end function

function set_phase(index_of_the_star, phase)
    implicit none
    integer:: index_of_the_star, phase
    integer:: set_phase
    GenecStar%phase = phase
    set_phase = 0
end function

function set_temperature_at_zone(index_of_the_star, zone, T_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: T_i
    integer:: set_temperature_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        GenecStar%t(i) = log(T_i)
    end if
    
    set_temperature_at_zone = 0
end function

end module AmuseInterface
