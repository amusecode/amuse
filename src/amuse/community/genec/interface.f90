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

integer function restore_star()
    implicit none
    restore_star = 0
end function

integer function initialize_code()
    !use WriteSaveClose, only: quitafterclosing
    use genec, only: initialise_genec
    use evol, only: input_dir
    use inputparam, only: libgenec
    use io_definitions
    implicit none

    !io_runfile = 6
    !io_logs = 6
    libgenec = .true.
    !io_logs = 6
    input_dir = "./src/GENEC/code"
    call initialise_genec()

    initialize_code = 0
end function

integer function read_genec_model(index_of_the_star, cardfilename)
    ! This should only be called if no star has been initialised yet!
    implicit none
    integer:: index_of_the_star
    character(256):: cardfilename
    read_genec_model = 0
end function

integer function cleanup_code()
    implicit none
    cleanup_code = 0
end function

! **** Parameters

integer function get_model_number(model_number)
    implicit none
    integer:: model_number
    model_number = GenecStar%nwmd
    get_model_number = 0
end function

integer function set_model_number(model_number)
    implicit none
    integer:: model_number
    GenecStar%nwmd = model_number
    set_model_number = 0
end function

integer function get_min_timestep_stop_condition(min_timestep_stop_condition)
    implicit none
    real(kindreal):: min_timestep_stop_condition
    get_min_timestep_stop_condition = 0
end function


!#########


integer function get_initialised(index_of_the_particle, initialised)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: initialised
    initialised = GenecStar%initialised
    get_initialised = 0
end function get_initialised

integer function set_initialised(index_of_the_particle, initialised)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: initialised
    GenecStar%initialised = initialised
    set_initialised = 0
end function set_initialised

integer function get_starname(index_of_the_particle, starname)
    implicit none
    integer, intent(in):: index_of_the_particle
    character(len=256), intent(out):: starname
    starname = GenecStar%starname
    get_starname = 0
end function get_starname

integer function set_starname(index_of_the_particle, starname)
    implicit none
    integer, intent(in):: index_of_the_particle
    character(len=256), intent(in):: starname
    GenecStar%starname = starname
    set_starname = 0
end function set_starname

integer function get_nwseq(index_of_the_particle, nwseq)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: nwseq
    nwseq = GenecStar%nwseq
    get_nwseq = 0
end function get_nwseq

integer function set_nwseq(index_of_the_particle, nwseq)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: nwseq
    GenecStar%nwseq = nwseq
    set_nwseq = 0
end function set_nwseq

integer function get_modanf(index_of_the_particle, modanf)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: modanf
    modanf = GenecStar%modanf
    get_modanf = 0
end function get_modanf

integer function set_modanf(index_of_the_particle, modanf)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: modanf
    GenecStar%modanf = modanf
    set_modanf = 0
end function set_modanf

integer function get_nzmod(index_of_the_particle, nzmod)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: nzmod
    nzmod = GenecStar%nzmod
    get_nzmod = 0
end function get_nzmod

integer function set_nzmod(index_of_the_particle, nzmod)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: nzmod
    GenecStar%nzmod = nzmod
    set_nzmod = 0
end function set_nzmod

integer function get_end_at_phase(index_of_the_particle, end_at_phase)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: end_at_phase
    end_at_phase = GenecStar%end_at_phase
    get_end_at_phase = 0
end function get_end_at_phase

integer function set_end_at_phase(index_of_the_particle, end_at_phase)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: end_at_phase
    GenecStar%end_at_phase = end_at_phase
    set_end_at_phase = 0
end function set_end_at_phase

integer function get_end_at_model(index_of_the_particle, end_at_model)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: end_at_model
    end_at_model = GenecStar%end_at_model
    get_end_at_model = 0
end function get_end_at_model

integer function set_end_at_model(index_of_the_particle, end_at_model)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: end_at_model
    GenecStar%end_at_model = end_at_model
    set_end_at_model = 0
end function set_end_at_model

integer function get_irot(index_of_the_particle, irot)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: irot
    irot = GenecStar%irot
    get_irot = 0
end function get_irot

integer function set_irot(index_of_the_particle, irot)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: irot
    GenecStar%irot = irot
    set_irot = 0
end function set_irot

integer function get_isol(index_of_the_particle, isol)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: isol
    isol = GenecStar%isol
    get_isol = 0
end function get_isol

integer function set_isol(index_of_the_particle, isol)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: isol
    GenecStar%isol = isol
    set_isol = 0
end function set_isol

integer function get_imagn(index_of_the_particle, imagn)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: imagn
    imagn = GenecStar%imagn
    get_imagn = 0
end function get_imagn

integer function set_imagn(index_of_the_particle, imagn)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: imagn
    GenecStar%imagn = imagn
    set_imagn = 0
end function set_imagn

integer function get_ialflu(index_of_the_particle, ialflu)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: ialflu
    ialflu = GenecStar%ialflu
    get_ialflu = 0
end function get_ialflu

integer function set_ialflu(index_of_the_particle, ialflu)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: ialflu
    GenecStar%ialflu = ialflu
    set_ialflu = 0
end function set_ialflu

integer function get_ianiso(index_of_the_particle, ianiso)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: ianiso
    ianiso = GenecStar%ianiso
    get_ianiso = 0
end function get_ianiso

integer function set_ianiso(index_of_the_particle, ianiso)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: ianiso
    GenecStar%ianiso = ianiso
    set_ianiso = 0
end function set_ianiso

integer function get_ipop3(index_of_the_particle, ipop3)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: ipop3
    ipop3 = GenecStar%ipop3
    get_ipop3 = 0
end function get_ipop3

integer function set_ipop3(index_of_the_particle, ipop3)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: ipop3
    GenecStar%ipop3 = ipop3
    set_ipop3 = 0
end function set_ipop3

integer function get_ibasnet(index_of_the_particle, ibasnet)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: ibasnet
    ibasnet = GenecStar%ibasnet
    get_ibasnet = 0
end function get_ibasnet

integer function set_ibasnet(index_of_the_particle, ibasnet)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: ibasnet
    GenecStar%ibasnet = ibasnet
    set_ibasnet = 0
end function set_ibasnet

integer function get_phase(index_of_the_particle, phase)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: phase
    phase = GenecStar%phase
    get_phase = 0
end function get_phase

integer function set_phase(index_of_the_particle, phase)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: phase
    GenecStar%phase = phase
    set_phase = 0
end function set_phase

integer function get_var_rates(index_of_the_particle, var_rates)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: var_rates
    var_rates = GenecStar%var_rates
    get_var_rates = 0
end function get_var_rates

integer function set_var_rates(index_of_the_particle, var_rates)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: var_rates
    GenecStar%var_rates = var_rates
    set_var_rates = 0
end function set_var_rates

integer function get_bintide(index_of_the_particle, bintide)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: bintide
    bintide = GenecStar%bintide
    get_bintide = 0
end function get_bintide

integer function set_bintide(index_of_the_particle, bintide)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: bintide
    GenecStar%bintide = bintide
    set_bintide = 0
end function set_bintide

integer function get_binm2(index_of_the_particle, binm2)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: binm2
    binm2 = GenecStar%binm2
    get_binm2 = 0
end function get_binm2

integer function set_binm2(index_of_the_particle, binm2)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: binm2
    GenecStar%binm2 = binm2
    set_binm2 = 0
end function set_binm2

integer function get_periodini(index_of_the_particle, periodini)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: periodini
    periodini = GenecStar%periodini
    get_periodini = 0
end function get_periodini

integer function set_periodini(index_of_the_particle, periodini)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: periodini
    GenecStar%periodini = periodini
    set_periodini = 0
end function set_periodini

integer function get_const_per(index_of_the_particle, const_per)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: const_per
    const_per = GenecStar%const_per
    get_const_per = 0
end function get_const_per

integer function set_const_per(index_of_the_particle, const_per)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: const_per
    GenecStar%const_per = const_per
    set_const_per = 0
end function set_const_per

integer function get_iprezams(index_of_the_particle, iprezams)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iprezams
    iprezams = GenecStar%iprezams
    get_iprezams = 0
end function get_iprezams

integer function set_iprezams(index_of_the_particle, iprezams)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iprezams
    GenecStar%iprezams = iprezams
    set_iprezams = 0
end function set_iprezams

integer function get_zinit(index_of_the_particle, zinit)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: zinit
    zinit = GenecStar%zinit
    get_zinit = 0
end function get_zinit

integer function set_zinit(index_of_the_particle, zinit)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: zinit
    GenecStar%zinit = zinit
    set_zinit = 0
end function set_zinit

integer function get_zsol(index_of_the_particle, zsol)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: zsol
    zsol = GenecStar%zsol
    get_zsol = 0
end function get_zsol

integer function set_zsol(index_of_the_particle, zsol)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: zsol
    GenecStar%zsol = zsol
    set_zsol = 0
end function set_zsol

integer function get_z(index_of_the_particle, z)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: z
    z = GenecStar%z
    get_z = 0
end function get_z

integer function set_z(index_of_the_particle, z)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: z
    GenecStar%z = z
    set_z = 0
end function set_z

integer function get_iopac(index_of_the_particle, iopac)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iopac
    iopac = GenecStar%iopac
    get_iopac = 0
end function get_iopac

integer function set_iopac(index_of_the_particle, iopac)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iopac
    GenecStar%iopac = iopac
    set_iopac = 0
end function set_iopac

integer function get_ikappa(index_of_the_particle, ikappa)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: ikappa
    ikappa = GenecStar%ikappa
    get_ikappa = 0
end function get_ikappa

integer function set_ikappa(index_of_the_particle, ikappa)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: ikappa
    GenecStar%ikappa = ikappa
    set_ikappa = 0
end function set_ikappa

integer function get_idiff(index_of_the_particle, idiff)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: idiff
    idiff = GenecStar%idiff
    get_idiff = 0
end function get_idiff

integer function set_idiff(index_of_the_particle, idiff)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: idiff
    GenecStar%idiff = idiff
    set_idiff = 0
end function set_idiff

integer function get_iadvec(index_of_the_particle, iadvec)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iadvec
    iadvec = GenecStar%iadvec
    get_iadvec = 0
end function get_iadvec

integer function set_iadvec(index_of_the_particle, iadvec)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iadvec
    GenecStar%iadvec = iadvec
    set_iadvec = 0
end function set_iadvec

integer function get_istati(index_of_the_particle, istati)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: istati
    istati = GenecStar%istati
    get_istati = 0
end function get_istati

integer function set_istati(index_of_the_particle, istati)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: istati
    GenecStar%istati = istati
    set_istati = 0
end function set_istati

integer function get_icoeff(index_of_the_particle, icoeff)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: icoeff
    icoeff = GenecStar%icoeff
    get_icoeff = 0
end function get_icoeff

integer function set_icoeff(index_of_the_particle, icoeff)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: icoeff
    GenecStar%icoeff = icoeff
    set_icoeff = 0
end function set_icoeff

integer function get_fenerg(index_of_the_particle, fenerg)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: fenerg
    fenerg = GenecStar%fenerg
    get_fenerg = 0
end function get_fenerg

integer function set_fenerg(index_of_the_particle, fenerg)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: fenerg
    GenecStar%fenerg = fenerg
    set_fenerg = 0
end function set_fenerg

integer function get_richac(index_of_the_particle, richac)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: richac
    richac = GenecStar%richac
    get_richac = 0
end function get_richac

integer function set_richac(index_of_the_particle, richac)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: richac
    GenecStar%richac = richac
    set_richac = 0
end function set_richac

integer function get_igamma(index_of_the_particle, igamma)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: igamma
    igamma = GenecStar%igamma
    get_igamma = 0
end function get_igamma

integer function set_igamma(index_of_the_particle, igamma)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: igamma
    GenecStar%igamma = igamma
    set_igamma = 0
end function set_igamma

integer function get_frein(index_of_the_particle, frein)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: frein
    frein = GenecStar%frein
    get_frein = 0
end function get_frein

integer function set_frein(index_of_the_particle, frein)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: frein
    GenecStar%frein = frein
    set_frein = 0
end function set_frein

integer function get_K_Kawaler(index_of_the_particle, K_Kawaler)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: K_Kawaler
    K_Kawaler = GenecStar%K_Kawaler
    get_K_Kawaler = 0
end function get_K_Kawaler

integer function set_K_Kawaler(index_of_the_particle, K_Kawaler)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: K_Kawaler
    GenecStar%K_Kawaler = K_Kawaler
    set_K_Kawaler = 0
end function set_K_Kawaler

integer function get_Omega_saturation(index_of_the_particle, Omega_saturation)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: Omega_saturation
    Omega_saturation = GenecStar%Omega_saturation
    get_Omega_saturation = 0
end function get_Omega_saturation

integer function set_Omega_saturation(index_of_the_particle, Omega_saturation)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: Omega_saturation
    GenecStar%Omega_saturation = Omega_saturation
    set_Omega_saturation = 0
end function set_Omega_saturation

integer function get_rapcrilim(index_of_the_particle, rapcrilim)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: rapcrilim
    rapcrilim = GenecStar%rapcrilim
    get_rapcrilim = 0
end function get_rapcrilim

integer function set_rapcrilim(index_of_the_particle, rapcrilim)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: rapcrilim
    GenecStar%rapcrilim = rapcrilim
    set_rapcrilim = 0
end function set_rapcrilim

integer function get_vwant(index_of_the_particle, vwant)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: vwant
    vwant = GenecStar%vwant
    get_vwant = 0
end function get_vwant

integer function set_vwant(index_of_the_particle, vwant)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: vwant
    GenecStar%vwant = vwant
    set_vwant = 0
end function set_vwant

integer function get_xfom(index_of_the_particle, xfom)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: xfom
    xfom = GenecStar%xfom
    get_xfom = 0
end function get_xfom

integer function set_xfom(index_of_the_particle, xfom)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: xfom
    GenecStar%xfom = xfom
    set_xfom = 0
end function set_xfom

integer function get_omega(index_of_the_particle, omega)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: omega
    omega = GenecStar%omega
    get_omega = 0
end function get_omega

integer function set_omega(index_of_the_particle, omega)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: omega
    GenecStar%omega = omega
    set_omega = 0
end function set_omega

integer function get_xdial(index_of_the_particle, xdial)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: xdial
    xdial = GenecStar%xdial
    get_xdial = 0
end function get_xdial

integer function set_xdial(index_of_the_particle, xdial)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: xdial
    GenecStar%xdial = xdial
    set_xdial = 0
end function set_xdial

integer function get_idialo(index_of_the_particle, idialo)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: idialo
    idialo = GenecStar%idialo
    get_idialo = 0
end function get_idialo

integer function set_idialo(index_of_the_particle, idialo)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: idialo
    GenecStar%idialo = idialo
    set_idialo = 0
end function set_idialo

integer function get_idialu(index_of_the_particle, idialu)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: idialu
    idialu = GenecStar%idialu
    get_idialu = 0
end function get_idialu

integer function set_idialu(index_of_the_particle, idialu)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: idialu
    GenecStar%idialu = idialu
    set_idialu = 0
end function set_idialu

integer function get_Add_Flux(index_of_the_particle, Add_Flux)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: Add_Flux
    Add_Flux = GenecStar%Add_Flux
    get_Add_Flux = 0
end function get_Add_Flux

integer function set_Add_Flux(index_of_the_particle, Add_Flux)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: Add_Flux
    GenecStar%Add_Flux = Add_Flux
    set_Add_Flux = 0
end function set_Add_Flux

integer function get_diff_only(index_of_the_particle, diff_only)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: diff_only
    diff_only = GenecStar%diff_only
    get_diff_only = 0
end function get_diff_only

integer function set_diff_only(index_of_the_particle, diff_only)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: diff_only
    GenecStar%diff_only = diff_only
    set_diff_only = 0
end function set_diff_only

integer function get_B_initial(index_of_the_particle, B_initial)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: B_initial
    B_initial = GenecStar%B_initial
    get_B_initial = 0
end function get_B_initial

integer function set_B_initial(index_of_the_particle, B_initial)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: B_initial
    GenecStar%B_initial = B_initial
    set_B_initial = 0
end function set_B_initial

integer function get_add_diff(index_of_the_particle, add_diff)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: add_diff
    add_diff = GenecStar%add_diff
    get_add_diff = 0
end function get_add_diff

integer function set_add_diff(index_of_the_particle, add_diff)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: add_diff
    GenecStar%add_diff = add_diff
    set_add_diff = 0
end function set_add_diff

integer function get_n_mag(index_of_the_particle, n_mag)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: n_mag
    n_mag = GenecStar%n_mag
    get_n_mag = 0
end function get_n_mag

integer function set_n_mag(index_of_the_particle, n_mag)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: n_mag
    GenecStar%n_mag = n_mag
    set_n_mag = 0
end function set_n_mag

integer function get_alpha_F(index_of_the_particle, alpha_F)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: alpha_F
    alpha_F = GenecStar%alpha_F
    get_alpha_F = 0
end function get_alpha_F

integer function set_alpha_F(index_of_the_particle, alpha_F)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: alpha_F
    GenecStar%alpha_F = alpha_F
    set_alpha_F = 0
end function set_alpha_F

integer function get_nsmooth(index_of_the_particle, nsmooth)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: nsmooth
    nsmooth = GenecStar%nsmooth
    get_nsmooth = 0
end function get_nsmooth

integer function set_nsmooth(index_of_the_particle, nsmooth)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: nsmooth
    GenecStar%nsmooth = nsmooth
    set_nsmooth = 0
end function set_nsmooth

integer function get_qminsmooth(index_of_the_particle, qminsmooth)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: qminsmooth
    qminsmooth = GenecStar%qminsmooth
    get_qminsmooth = 0
end function get_qminsmooth

integer function set_qminsmooth(index_of_the_particle, qminsmooth)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: qminsmooth
    GenecStar%qminsmooth = qminsmooth
    set_qminsmooth = 0
end function set_qminsmooth

integer function get_imloss(index_of_the_particle, imloss)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: imloss
    imloss = GenecStar%imloss
    get_imloss = 0
end function get_imloss

integer function set_imloss(index_of_the_particle, imloss)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: imloss
    GenecStar%imloss = imloss
    set_imloss = 0
end function set_imloss

integer function get_fmlos(index_of_the_particle, fmlos)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: fmlos
    fmlos = GenecStar%fmlos
    get_fmlos = 0
end function get_fmlos

integer function set_fmlos(index_of_the_particle, fmlos)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: fmlos
    GenecStar%fmlos = fmlos
    set_fmlos = 0
end function set_fmlos

integer function get_ifitm(index_of_the_particle, ifitm)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: ifitm
    ifitm = GenecStar%ifitm
    get_ifitm = 0
end function get_ifitm

integer function set_ifitm(index_of_the_particle, ifitm)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: ifitm
    GenecStar%ifitm = ifitm
    set_ifitm = 0
end function set_ifitm

integer function get_fitm(index_of_the_particle, fitm)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: fitm
    fitm = GenecStar%fitm
    get_fitm = 0
end function get_fitm

integer function set_fitm(index_of_the_particle, fitm)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: fitm
    GenecStar%fitm = fitm
    set_fitm = 0
end function set_fitm

integer function get_fitmi(index_of_the_particle, fitmi)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: fitmi
    fitmi = GenecStar%fitmi
    get_fitmi = 0
end function get_fitmi

integer function set_fitmi(index_of_the_particle, fitmi)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: fitmi
    GenecStar%fitmi = fitmi
    set_fitmi = 0
end function set_fitmi

integer function get_deltal(index_of_the_particle, deltal)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: deltal
    deltal = GenecStar%deltal
    get_deltal = 0
end function get_deltal

integer function set_deltal(index_of_the_particle, deltal)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: deltal
    GenecStar%deltal = deltal
    set_deltal = 0
end function set_deltal

integer function get_deltat(index_of_the_particle, deltat)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: deltat
    deltat = GenecStar%deltat
    get_deltat = 0
end function get_deltat

integer function set_deltat(index_of_the_particle, deltat)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: deltat
    GenecStar%deltat = deltat
    set_deltat = 0
end function set_deltat

integer function get_nndr(index_of_the_particle, nndr)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: nndr
    nndr = GenecStar%nndr
    get_nndr = 0
end function get_nndr

integer function set_nndr(index_of_the_particle, nndr)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: nndr
    GenecStar%nndr = nndr
    set_nndr = 0
end function set_nndr

integer function get_RSG_Mdot(index_of_the_particle, RSG_Mdot)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: RSG_Mdot
    RSG_Mdot = GenecStar%RSG_Mdot
    get_RSG_Mdot = 0
end function get_RSG_Mdot

integer function set_RSG_Mdot(index_of_the_particle, RSG_Mdot)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: RSG_Mdot
    GenecStar%RSG_Mdot = RSG_Mdot
    set_RSG_Mdot = 0
end function set_RSG_Mdot

integer function get_SupraEddMdot(index_of_the_particle, SupraEddMdot)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: SupraEddMdot
    SupraEddMdot = GenecStar%SupraEddMdot
    get_SupraEddMdot = 0
end function get_SupraEddMdot

integer function set_SupraEddMdot(index_of_the_particle, SupraEddMdot)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: SupraEddMdot
    GenecStar%SupraEddMdot = SupraEddMdot
    set_SupraEddMdot = 0
end function set_SupraEddMdot

integer function get_Be_mdotfrac(index_of_the_particle, Be_mdotfrac)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: Be_mdotfrac
    Be_mdotfrac = GenecStar%Be_mdotfrac
    get_Be_mdotfrac = 0
end function get_Be_mdotfrac

integer function set_Be_mdotfrac(index_of_the_particle, Be_mdotfrac)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: Be_mdotfrac
    GenecStar%Be_mdotfrac = Be_mdotfrac
    set_Be_mdotfrac = 0
end function set_Be_mdotfrac

integer function get_start_mdot(index_of_the_particle, start_mdot)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: start_mdot
    start_mdot = GenecStar%start_mdot
    get_start_mdot = 0
end function get_start_mdot

integer function set_start_mdot(index_of_the_particle, start_mdot)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: start_mdot
    GenecStar%start_mdot = start_mdot
    set_start_mdot = 0
end function set_start_mdot

integer function get_iledou(index_of_the_particle, iledou)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iledou
    iledou = GenecStar%iledou
    get_iledou = 0
end function get_iledou

integer function set_iledou(index_of_the_particle, iledou)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iledou
    GenecStar%iledou = iledou
    set_iledou = 0
end function set_iledou

integer function get_idifcon(index_of_the_particle, idifcon)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: idifcon
    idifcon = GenecStar%idifcon
    get_idifcon = 0
end function get_idifcon

integer function set_idifcon(index_of_the_particle, idifcon)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: idifcon
    GenecStar%idifcon = idifcon
    set_idifcon = 0
end function set_idifcon

integer function get_iover(index_of_the_particle, iover)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iover
    iover = GenecStar%iover
    get_iover = 0
end function get_iover

integer function set_iover(index_of_the_particle, iover)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iover
    GenecStar%iover = iover
    set_iover = 0
end function set_iover

integer function get_elph(index_of_the_particle, elph)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: elph
    elph = GenecStar%elph
    get_elph = 0
end function get_elph

integer function set_elph(index_of_the_particle, elph)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: elph
    GenecStar%elph = elph
    set_elph = 0
end function set_elph

integer function get_my(index_of_the_particle, my)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: my
    my = GenecStar%my
    get_my = 0
end function get_my

integer function set_my(index_of_the_particle, my)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: my
    GenecStar%my = my
    set_my = 0
end function set_my

integer function get_dovhp(index_of_the_particle, dovhp)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dovhp
    dovhp = GenecStar%dovhp
    get_dovhp = 0
end function get_dovhp

integer function set_dovhp(index_of_the_particle, dovhp)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dovhp
    GenecStar%dovhp = dovhp
    set_dovhp = 0
end function set_dovhp

integer function get_iunder(index_of_the_particle, iunder)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iunder
    iunder = GenecStar%iunder
    get_iunder = 0
end function get_iunder

integer function set_iunder(index_of_the_particle, iunder)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iunder
    GenecStar%iunder = iunder
    set_iunder = 0
end function set_iunder

integer function get_dunder(index_of_the_particle, dunder)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dunder
    dunder = GenecStar%dunder
    get_dunder = 0
end function get_dunder

integer function set_dunder(index_of_the_particle, dunder)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dunder
    GenecStar%dunder = dunder
    set_dunder = 0
end function set_dunder

integer function get_gkorm(index_of_the_particle, gkorm)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: gkorm
    gkorm = GenecStar%gkorm
    get_gkorm = 0
end function get_gkorm

integer function set_gkorm(index_of_the_particle, gkorm)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: gkorm
    GenecStar%gkorm = gkorm
    set_gkorm = 0
end function set_gkorm

integer function get_alph(index_of_the_particle, alph)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: alph
    alph = GenecStar%alph
    get_alph = 0
end function get_alph

integer function set_alph(index_of_the_particle, alph)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: alph
    GenecStar%alph = alph
    set_alph = 0
end function set_alph

integer function get_agdr(index_of_the_particle, agdr)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: agdr
    agdr = GenecStar%agdr
    get_agdr = 0
end function get_agdr

integer function set_agdr(index_of_the_particle, agdr)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: agdr
    GenecStar%agdr = agdr
    set_agdr = 0
end function set_agdr

integer function get_faktor(index_of_the_particle, faktor)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: faktor
    faktor = GenecStar%faktor
    get_faktor = 0
end function get_faktor

integer function set_faktor(index_of_the_particle, faktor)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: faktor
    GenecStar%faktor = faktor
    set_faktor = 0
end function set_faktor

integer function get_dgrp(index_of_the_particle, dgrp)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dgrp
    dgrp = GenecStar%dgrp
    get_dgrp = 0
end function get_dgrp

integer function set_dgrp(index_of_the_particle, dgrp)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dgrp
    GenecStar%dgrp = dgrp
    set_dgrp = 0
end function set_dgrp

integer function get_dgrl(index_of_the_particle, dgrl)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dgrl
    dgrl = GenecStar%dgrl
    get_dgrl = 0
end function get_dgrl

integer function set_dgrl(index_of_the_particle, dgrl)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dgrl
    GenecStar%dgrl = dgrl
    set_dgrl = 0
end function set_dgrl

integer function get_dgry(index_of_the_particle, dgry)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dgry
    dgry = GenecStar%dgry
    get_dgry = 0
end function get_dgry

integer function set_dgry(index_of_the_particle, dgry)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dgry
    GenecStar%dgry = dgry
    set_dgry = 0
end function set_dgry

integer function get_dgrc(index_of_the_particle, dgrc)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dgrc
    dgrc = GenecStar%dgrc
    get_dgrc = 0
end function get_dgrc

integer function set_dgrc(index_of_the_particle, dgrc)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dgrc
    GenecStar%dgrc = dgrc
    set_dgrc = 0
end function set_dgrc

integer function get_dgro(index_of_the_particle, dgro)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dgro
    dgro = GenecStar%dgro
    get_dgro = 0
end function get_dgro

integer function set_dgro(index_of_the_particle, dgro)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dgro
    GenecStar%dgro = dgro
    set_dgro = 0
end function set_dgro

integer function get_dgr20(index_of_the_particle, dgr20)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dgr20
    dgr20 = GenecStar%dgr20
    get_dgr20 = 0
end function get_dgr20

integer function set_dgr20(index_of_the_particle, dgr20)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dgr20
    GenecStar%dgr20 = dgr20
    set_dgr20 = 0
end function set_dgr20

integer function get_nbchx(index_of_the_particle, nbchx)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: nbchx
    nbchx = GenecStar%nbchx
    get_nbchx = 0
end function get_nbchx

integer function set_nbchx(index_of_the_particle, nbchx)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: nbchx
    GenecStar%nbchx = nbchx
    set_nbchx = 0
end function set_nbchx

integer function get_nrband(index_of_the_particle, nrband)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: nrband
    nrband = GenecStar%nrband
    get_nrband = 0
end function get_nrband

integer function set_nrband(index_of_the_particle, nrband)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: nrband
    GenecStar%nrband = nrband
    set_nrband = 0
end function set_nrband

integer function get_xcn(index_of_the_particle, xcn)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: xcn
    xcn = GenecStar%xcn
    get_xcn = 0
end function get_xcn

integer function set_xcn(index_of_the_particle, xcn)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: xcn
    GenecStar%xcn = xcn
    set_xcn = 0
end function set_xcn

integer function get_islow(index_of_the_particle, islow)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: islow
    islow = GenecStar%islow
    get_islow = 0
end function get_islow

integer function set_islow(index_of_the_particle, islow)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: islow
    GenecStar%islow = islow
    set_islow = 0
end function set_islow

integer function get_icncst(index_of_the_particle, icncst)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: icncst
    icncst = GenecStar%icncst
    get_icncst = 0
end function get_icncst

integer function set_icncst(index_of_the_particle, icncst)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: icncst
    GenecStar%icncst = icncst
    set_icncst = 0
end function set_icncst

integer function get_tauH_fit(index_of_the_particle, tauH_fit)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: tauH_fit
    tauH_fit = GenecStar%tauH_fit
    get_tauH_fit = 0
end function get_tauH_fit

integer function set_tauH_fit(index_of_the_particle, tauH_fit)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: tauH_fit
    GenecStar%tauH_fit = tauH_fit
    set_tauH_fit = 0
end function set_tauH_fit

integer function get_display_plot(index_of_the_particle, display_plot)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: display_plot
    display_plot = GenecStar%display_plot
    get_display_plot = 0
end function get_display_plot

integer function set_display_plot(index_of_the_particle, display_plot)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: display_plot
    GenecStar%display_plot = display_plot
    set_display_plot = 0
end function set_display_plot

integer function get_iauto(index_of_the_particle, iauto)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iauto
    iauto = GenecStar%iauto
    get_iauto = 0
end function get_iauto

integer function set_iauto(index_of_the_particle, iauto)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iauto
    GenecStar%iauto = iauto
    set_iauto = 0
end function set_iauto

integer function get_iprn(index_of_the_particle, iprn)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iprn
    iprn = GenecStar%iprn
    get_iprn = 0
end function get_iprn

integer function set_iprn(index_of_the_particle, iprn)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iprn
    GenecStar%iprn = iprn
    set_iprn = 0
end function set_iprn

integer function get_iout(index_of_the_particle, iout)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: iout
    iout = GenecStar%iout
    get_iout = 0
end function get_iout

integer function set_iout(index_of_the_particle, iout)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: iout
    GenecStar%iout = iout
    set_iout = 0
end function set_iout

integer function get_itmin(index_of_the_particle, itmin)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: itmin
    itmin = GenecStar%itmin
    get_itmin = 0
end function get_itmin

integer function set_itmin(index_of_the_particle, itmin)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: itmin
    GenecStar%itmin = itmin
    set_itmin = 0
end function set_itmin

integer function get_xyfiles(index_of_the_particle, xyfiles)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: xyfiles
    xyfiles = GenecStar%xyfiles
    get_xyfiles = 0
end function get_xyfiles

integer function set_xyfiles(index_of_the_particle, xyfiles)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: xyfiles
    GenecStar%xyfiles = xyfiles
    set_xyfiles = 0
end function set_xyfiles

integer function get_idebug(index_of_the_particle, idebug)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: idebug
    idebug = GenecStar%idebug
    get_idebug = 0
end function get_idebug

integer function set_idebug(index_of_the_particle, idebug)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: idebug
    GenecStar%idebug = idebug
    set_idebug = 0
end function set_idebug

integer function get_itests(index_of_the_particle, itests)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: itests
    itests = GenecStar%itests
    get_itests = 0
end function get_itests

integer function set_itests(index_of_the_particle, itests)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: itests
    GenecStar%itests = itests
    set_itests = 0
end function set_itests

integer function get_verbose(index_of_the_particle, verbose)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: verbose
    verbose = GenecStar%verbose
    get_verbose = 0
end function get_verbose

integer function set_verbose(index_of_the_particle, verbose)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: verbose
    GenecStar%verbose = verbose
    set_verbose = 0
end function set_verbose

integer function get_stop_deg(index_of_the_particle, stop_deg)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(out):: stop_deg
    stop_deg = GenecStar%stop_deg
    get_stop_deg = 0
end function get_stop_deg

integer function set_stop_deg(index_of_the_particle, stop_deg)
    implicit none
    integer, intent(in):: index_of_the_particle
    logical, intent(in):: stop_deg
    GenecStar%stop_deg = stop_deg
    set_stop_deg = 0
end function set_stop_deg

integer function get_n_snap(index_of_the_particle, n_snap)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: n_snap
    n_snap = GenecStar%n_snap
    get_n_snap = 0
end function get_n_snap

integer function set_n_snap(index_of_the_particle, n_snap)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: n_snap
    GenecStar%n_snap = n_snap
    set_n_snap = 0
end function set_n_snap

! **** End Parameters

! **** Begin Properties
integer function get_m(index_of_the_particle, m)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(out):: m
    m = GenecStar%m
    get_m = 0
end function get_m

integer function set_m(index_of_the_particle, m)
    implicit none
    integer, intent(in):: index_of_the_particle
    integer, intent(in):: m
    GenecStar%m = m
    set_m = 0
end function set_m

integer function get_gms(index_of_the_particle, gms)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: gms
    gms = GenecStar%gms
    get_gms = 0
end function get_gms

integer function set_gms(index_of_the_particle, gms)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: gms
    GenecStar%gms = gms
    set_gms = 0
end function set_gms

integer function get_alter(index_of_the_particle, alter)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: alter
    alter = GenecStar%alter
    get_alter = 0
end function get_alter

integer function set_alter(index_of_the_particle, alter)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: alter
    GenecStar%alter = alter
    set_alter = 0
end function set_alter

integer function get_gls(index_of_the_particle, gls)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: gls
    gls = GenecStar%gls
    get_gls = 0
end function get_gls

integer function set_gls(index_of_the_particle, gls)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: gls
    GenecStar%gls = gls
    set_gls = 0
end function set_gls

integer function get_teff(index_of_the_particle, teff)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: teff
    teff = GenecStar%teff
    get_teff = 0
end function get_teff

integer function set_teff(index_of_the_particle, teff)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: teff
    GenecStar%teff = teff
    set_teff = 0
end function set_teff

integer function get_glsv(index_of_the_particle, glsv)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: glsv
    glsv = GenecStar%glsv
    get_glsv = 0
end function get_glsv

integer function set_glsv(index_of_the_particle, glsv)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: glsv
    GenecStar%glsv = glsv
    set_glsv = 0
end function set_glsv

integer function get_teffv(index_of_the_particle, teffv)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: teffv
    teffv = GenecStar%teffv
    get_teffv = 0
end function get_teffv

integer function set_teffv(index_of_the_particle, teffv)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: teffv
    GenecStar%teffv = teffv
    set_teffv = 0
end function set_teffv

integer function get_dzeitj(index_of_the_particle, dzeitj)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dzeitj
    dzeitj = GenecStar%dzeitj
    get_dzeitj = 0
end function get_dzeitj

integer function set_dzeitj(index_of_the_particle, dzeitj)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dzeitj
    GenecStar%dzeitj = dzeitj
    set_dzeitj = 0
end function set_dzeitj

integer function get_dzeit(index_of_the_particle, dzeit)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dzeit
    dzeit = GenecStar%dzeit
    get_dzeit = 0
end function get_dzeit

integer function set_dzeit(index_of_the_particle, dzeit)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dzeit
    GenecStar%dzeit = dzeit
    set_dzeit = 0
end function set_dzeit

integer function get_dzeitv(index_of_the_particle, dzeitv)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dzeitv
    dzeitv = GenecStar%dzeitv
    get_dzeitv = 0
end function get_dzeitv

integer function set_dzeitv(index_of_the_particle, dzeitv)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dzeitv
    GenecStar%dzeitv = dzeitv
    set_dzeitv = 0
end function set_dzeitv

integer function get_xmini(index_of_the_particle, xmini)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: xmini
    xmini = GenecStar%xmini
    get_xmini = 0
end function get_xmini

integer function set_xmini(index_of_the_particle, xmini)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: xmini
    GenecStar%xmini = xmini
    set_xmini = 0
end function set_xmini

integer function get_summas(index_of_the_particle, summas)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: summas
    summas = GenecStar%summas
    get_summas = 0
end function get_summas

integer function set_summas(index_of_the_particle, summas)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: summas
    GenecStar%summas = summas
    set_summas = 0
end function set_summas

integer function get_ab(index_of_the_particle, ab)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: ab
    ab = GenecStar%ab
    get_ab = 0
end function get_ab

integer function set_ab(index_of_the_particle, ab)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: ab
    GenecStar%ab = ab
    set_ab = 0
end function set_ab

integer function get_dm_lost(index_of_the_particle, dm_lost)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(out):: dm_lost
    dm_lost = GenecStar%dm_lost
    get_dm_lost = 0
end function get_dm_lost

integer function set_dm_lost(index_of_the_particle, dm_lost)
    implicit none
    integer, intent(in):: index_of_the_particle
    real(kindreal), intent(in):: dm_lost
    GenecStar%dm_lost = dm_lost
    set_dm_lost = 0
end function set_dm_lost

! **** End Properties


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

function get_eps_ne_adv_at_zone(index_of_the_star, zone, eps_ne_adv)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: eps_ne_adv
    integer:: get_eps_ne_adv_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        eps_ne_adv = GenecStar%eps_ne_adv(i)
    end if
    get_eps_ne_adv_at_zone = 0
end function

function get_eps_o_adv_at_zone(index_of_the_star, zone, eps_o_adv)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: eps_o_adv
    integer:: get_eps_o_adv_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        eps_o_adv = GenecStar%eps_o_adv(i)
    end if
    get_eps_o_adv_at_zone = 0
end function

function get_eps_si_adv_at_zone(index_of_the_star, zone, eps_si_adv)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: eps_si_adv
    integer:: get_eps_si_adv_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        eps_si_adv = GenecStar%eps_si_adv(i)
    end if
    get_eps_si_adv_at_zone = 0
end function

function get_eps_grav_at_zone(index_of_the_star, zone, eps_grav)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: eps_grav
    integer:: get_eps_grav_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        eps_grav = GenecStar%eps_grav(i)
    end if
    get_eps_grav_at_zone = 0
end function

function get_eps_nu_at_zone(index_of_the_star, zone, eps_nu)
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    real(kindreal):: eps_nu
    integer:: get_eps_nu_at_zone
    i = GenecStar%m - zone
    if (zone <= GenecStar%m) then
        eps_nu = GenecStar%eps_nu(i)
    end if
    get_eps_nu_at_zone = 0
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

function new_stellar_model(&
      integer_of_the_star,&
      initialised, starname, nwseq, modanf, nzmod, end_at_phase, end_at_model, irot, isol, imagn, ialflu, ianiso, ipop3,&
      ibasnet, phase, var_rates, bintide, binm2, periodini, const_per, iprezams, zinit, zsol, z, iopac, ikappa, idiff, iadvec,&
      istati, icoeff, fenerg, richac, igamma, frein, K_Kawaler, Omega_saturation, rapcrilim, vwant, xfom, omega, xdial, idialo,&
      idialu, Add_Flux, diff_only, B_initial, add_diff, n_mag, alpha_F, nsmooth, qminsmooth, imloss, fmlos, ifitm, fitm, fitmi,&
      deltal, deltat, nndr, RSG_Mdot, SupraEddMdot, Be_mdotfrac, start_mdot, iledou, idifcon, iover, elph, my, dovhp, iunder,&
      dunder, gkorm, alph, agdr, faktor, dgrp, dgrl, dgry, dgrc, dgro, dgr20, nbchx, nrband, xcn, islow, icncst, tauH_fit,&
      display_plot, iauto, iprn, iout, itmin, xyfiles, idebug, itests, verbose, stop_deg, n_snap,&
      m,gms,alter,gls,teff,glsv,teffv,dzeitj,dzeit,dzeitv,xmini,summas,ab,dm_lost,&
      q,p,t,r,s,x,y3,y,xc12,xc13,xn14,xn15,xo16,xo17,xo18,xne20,xne22,xmg24,xmg25,xmg26,xf19,xne21,xna23,xal27,xsi28,&
      xc14,xf18,xal26,xneut,xprot,omegi,xbid,xbid1,vp,vt,vr,vs,vx,vy,vy3,vxc12,vxc13,vxn14,vxn15,vxo16,vxo17,vxo18,&
      vxne20,vxne22,vxmg24,vxmg25,vxmg26,vxf19,vxne21,vxna23,vxal27,vxsi28,vxc14,vxf18,vxal26g,vxneut,vxprot,vomegi,&
      vxbid,vxbid1&
      )
    implicit none
    integer:: integer_of_the_star
    logical, intent(in):: &
            initialised,var_rates,bintide,const_per,Add_Flux,diff_only,qminsmooth,SupraEddMdot,display_plot,xyfiles,verbose,&
            stop_deg
    integer, intent(in):: &
            nwseq,modanf,nzmod,end_at_phase,end_at_model,irot,isol,imagn,ialflu,ianiso,ipop3,ibasnet,phase,iprezams,iopac,ikappa,&
            idiff,iadvec,istati,icoeff,igamma,idialo,idialu,n_mag,nsmooth,imloss,ifitm,nndr,RSG_Mdot,iledou,idifcon,iover,my,&
            iunder,nbchx,nrband,islow,icncst,tauH_fit,iauto,iprn,iout,itmin,idebug,itests,n_snap
    real(kindreal), intent(in):: &
            binm2,periodini,zinit,zsol,z,fenerg,richac,frein,K_Kawaler,Omega_saturation,rapcrilim,vwant,xfom,omega,xdial,&
            B_initial,add_diff,alpha_F,fmlos,fitm,fitmi,deltal,deltat,Be_mdotfrac,start_mdot,elph,dovhp,dunder,gkorm,alph,&
            agdr,faktor,dgrp,dgrl,dgry,dgrc,dgro,dgr20,xcn
    character(len=200), intent(in):: &
            starname
    integer, intent(in) :: m
    real(kindreal), intent(in) :: &
            gms,alter,gls,teff,glsv,teffv,dzeitj,dzeit,dzeitv,xmini,summas,ab,&
            dm_lost
    !real(kindreal), dimension(ldi) :: &
    real(kindreal), dimension(m), intent(in):: &
            q,p,t,r,s,x,y3,y,xc12,xc13,xn14,xn15,xo16,xo17,xo18,xne20,xne22,xmg24,xmg25,xmg26,xf19,xne21,xna23,xal27,xsi28,&
            xc14,xf18,xal26,xneut,xprot,omegi,xbid,xbid1,vp,vt,vr,vs,vx,vy,vy3,vxc12,vxc13,vxn14,vxn15,vxo16,vxo17,vxo18,&
            vxne20,vxne22,vxmg24,vxmg25,vxmg26,vxf19,vxne21,vxna23,vxal27,vxsi28,vxc14,vxf18,vxal26g,vxneut,vxprot,vomegi,&
            vxbid,vxbid1
    integer:: new_stellar_model

    GenecStar%initialised      = initialised
    GenecStar%starname         = starname
    !GenecStar%nwmd             = nwmd
    GenecStar%nwseq            = nwseq            
    GenecStar%modanf           = modanf
    GenecStar%nzmod            = nzmod
    GenecStar%end_at_phase     = end_at_phase
    GenecStar%end_at_model     = end_at_model
    GenecStar%irot             = irot
    GenecStar%isol             = isol
    GenecStar%imagn            = imagn
    GenecStar%ialflu           = ialflu
    GenecStar%ianiso           = ianiso
    GenecStar%ipop3            = ipop3
    GenecStar%ibasnet          = ibasnet
    GenecStar%phase            = phase
    GenecStar%iprezams         = iprezams
    GenecStar%var_rates        = var_rates
    GenecStar%bintide          = bintide
    GenecStar%binm2            = binm2
    GenecStar%periodini        = periodini
    GenecStar%const_per        = const_per
    GenecStar%iopac            = iopac
    GenecStar%ikappa           = ikappa
    GenecStar%zinit            = zinit
    GenecStar%zsol             = zsol
    GenecStar%z                = z
    GenecStar%idiff            = idiff
    GenecStar%iadvec           = iadvec
    GenecStar%istati           = istati
    GenecStar%icoeff           = icoeff
    GenecStar%igamma           = igamma
    GenecStar%idialo           = idialo
    GenecStar%idialu           = idialu
    GenecStar%n_mag            = n_mag
    GenecStar%nsmooth          = nsmooth
    GenecStar%fenerg           = fenerg
    GenecStar%richac           = richac
    GenecStar%frein            = frein
    GenecStar%K_Kawaler        = K_Kawaler
    GenecStar%Omega_saturation = Omega_saturation
    GenecStar%rapcrilim        = rapcrilim
    GenecStar%vwant            = vwant
    GenecStar%xfom             = xfom
    GenecStar%omega            = omega
    GenecStar%xdial            = xdial
    GenecStar%B_initial        = B_initial
    GenecStar%add_diff         = add_diff
    GenecStar%alpha_F          = alpha_F
    GenecStar%Add_Flux         = Add_Flux
    GenecStar%diff_only        = diff_only
    GenecStar%qminsmooth       = qminsmooth
    GenecStar%imloss           = imloss
    GenecStar%ifitm            = ifitm
    GenecStar%nndr             = nndr
    GenecStar%RSG_Mdot         = RSG_Mdot
    GenecStar%fmlos            = fmlos
    GenecStar%fitm             = fitm
    GenecStar%fitmi            = fitmi
    GenecStar%deltal           = deltal
    GenecStar%deltat           = deltat
    GenecStar%Be_mdotfrac      = Be_mdotfrac
    GenecStar%start_mdot       = start_mdot
    GenecStar%SupraEddMdot     = SupraEddMdot
    GenecStar%iledou           = iledou
    GenecStar%idifcon          = idifcon
    GenecStar%my               = my
    GenecStar%iover            = iover
    GenecStar%iunder           = iunder
    GenecStar%elph             = elph
    GenecStar%dovhp            = dovhp
    GenecStar%dunder           = dunder
    GenecStar%nbchx            = nbchx
    GenecStar%nrband           = nrband
    GenecStar%gkorm            = gkorm
    GenecStar%alph             = alph
    GenecStar%agdr             = agdr
    GenecStar%faktor           = faktor
    GenecStar%dgrp             = dgrp
    GenecStar%dgrl             = dgrl
    GenecStar%dgry             = dgry
    GenecStar%dgrc             = dgrc
    GenecStar%dgro             = dgro
    GenecStar%dgr20            = dgr20
    GenecStar%islow            = islow
    GenecStar%icncst           = icncst
    GenecStar%tauH_fit         = tauH_fit
    GenecStar%xcn              = xcn
    GenecStar%iauto            = iauto
    GenecStar%iprn             = iprn
    GenecStar%iout             = iout
    GenecStar%itmin            = itmin
    GenecStar%idebug           = idebug
    GenecStar%itests           = itests
    GenecStar%n_snap           = n_snap
    GenecStar%display_plot     = display_plot
    GenecStar%xyfiles          = xyfiles
    GenecStar%verbose          = verbose
    GenecStar%stop_deg         = stop_deg

    GenecStar%m                = m
    GenecStar%gms              = gms
    GenecStar%alter            = alter
    GenecStar%gls              = gls
    GenecStar%teff             = teff
    GenecStar%glsv             = glsv
    GenecStar%teffv            = teffv
    GenecStar%dzeitj           = dzeitj
    GenecStar%dzeit            = dzeit
    GenecStar%dzeitv           = dzeitv
    GenecStar%summas           = summas
    GenecStar%xmini            = xmini
    GenecStar%ab               = ab
    GenecStar%dm_lost          = dm_lost
    GenecStar%q                = q
    GenecStar%p                = p
    GenecStar%t                = t
    GenecStar%r                = r
    GenecStar%s                = s
    GenecStar%x                = x
    GenecStar%y                = y
    GenecStar%xc12             = xc12
    GenecStar%vp               = vp
    GenecStar%vt               = vt
    GenecStar%vr               = vr
    GenecStar%vs               = vs
    GenecStar%xo16             = xo16
    GenecStar%vx               = vx
    GenecStar%vy               = vy
    GenecStar%vxc12            = vxc12
    GenecStar%vxo16            = vxo16
    GenecStar%y3               = y3
    GenecStar%xc13             = xc13
    GenecStar%xn14             = xn14
    GenecStar%xn15             = xn15
    GenecStar%xo17             = xo17
    GenecStar%xo18             = xo18
    GenecStar%vy3              = vy3
    GenecStar%vxc13            = vxc13
    GenecStar%vxn14            = vxn14
    GenecStar%vxn15            = vxn15
    GenecStar%vxo17            = vxo17
    GenecStar%vxo18            = vxo18
    GenecStar%xne20            = xne20
    GenecStar%xne22            = xne22
    GenecStar%xmg24            = xmg24
    GenecStar%xmg25            = xmg25
    GenecStar%xmg26            = xmg26
    GenecStar%vxne20           = vxne20
    GenecStar%vxne22           = vxne22
    GenecStar%vxmg24           = vxmg24
    GenecStar%vxmg25           = vxmg25
    GenecStar%vxmg26           = vxmg26
    GenecStar%omegi            = omegi
    GenecStar%vomegi           = vomegi
    GenecStar%xf19             = xf19
    GenecStar%xne21            = xne21
    GenecStar%xna23            = xna23
    GenecStar%xal26            = xal26
    GenecStar%xal27            = xal27
    GenecStar%xsi28            = xsi28
    GenecStar%vxf19            = vxf19
    GenecStar%vxne21           = vxne21
    GenecStar%vxna23           = vxna23
    GenecStar%vxal26g          = vxal26g
    GenecStar%vxal27           = vxal27
    GenecStar%vxsi28           = vxsi28
    GenecStar%xneut            = xneut
    GenecStar%xprot            = xprot
    GenecStar%xc14             = xc14
    GenecStar%xf18             = xf18
    GenecStar%xbid             = xbid
    GenecStar%xbid1            = xbid1
    GenecStar%vxneut           = vxneut
    GenecStar%vxprot           = vxprot
    GenecStar%vxc14            = vxc14
    GenecStar%vxf18            = vxf18
    GenecStar%vxbid            = vxbid
    GenecStar%vxbid1           = vxbid1

    new_stellar_model = 0
end function

function get_stellar_model(&
      integer_of_the_star,&
      initialised, starname, nwseq, modanf, nzmod, end_at_phase, end_at_model, irot, isol, imagn, ialflu, ianiso, ipop3,&
      ibasnet, phase, var_rates, bintide, binm2, periodini, const_per, iprezams, zinit, zsol, z, iopac, ikappa, idiff, iadvec,&
      istati, icoeff, fenerg, richac, igamma, frein, K_Kawaler, Omega_saturation, rapcrilim, vwant, xfom, omega, xdial, idialo,&
      idialu, Add_Flux, diff_only, B_initial, add_diff, n_mag, alpha_F, nsmooth, qminsmooth, imloss, fmlos, ifitm, fitm, fitmi,&
      deltal, deltat, nndr, RSG_Mdot, SupraEddMdot, Be_mdotfrac, start_mdot, iledou, idifcon, iover, elph, my, dovhp, iunder,&
      dunder, gkorm, alph, agdr, faktor, dgrp, dgrl, dgry, dgrc, dgro, dgr20, nbchx, nrband, xcn, islow, icncst, tauH_fit,&
      display_plot, iauto, iprn, iout, itmin, xyfiles, idebug, itests, verbose, stop_deg, n_snap,&
      m,gms,alter,gls,teff,glsv,teffv,dzeitj,dzeit,dzeitv,xmini,summas,ab,dm_lost,&
      q,p,t,r,s,x,y3,y,xc12,xc13,xn14,xn15,xo16,xo17,xo18,xne20,xne22,xmg24,xmg25,xmg26,xf19,xne21,xna23,xal27,xsi28,&
      xc14,xf18,xal26,xneut,xprot,omegi,xbid,xbid1,vp,vt,vr,vs,vx,vy,vy3,vxc12,vxc13,vxn14,vxn15,vxo16,vxo17,vxo18,&
      vxne20,vxne22,vxmg24,vxmg25,vxmg26,vxf19,vxne21,vxna23,vxal27,vxsi28,vxc14,vxf18,vxal26g,vxneut,vxprot,vomegi,&
      vxbid,vxbid1&
      )
    implicit none
    integer:: integer_of_the_star
    logical, intent(out):: &
            initialised,var_rates,bintide,const_per,Add_Flux,diff_only,qminsmooth,SupraEddMdot,display_plot,xyfiles,verbose,&
            stop_deg
    integer, intent(out):: &
            nwseq,modanf,nzmod,end_at_phase,end_at_model,irot,isol,imagn,ialflu,ianiso,ipop3,ibasnet,phase,iprezams,iopac,ikappa,&
            idiff,iadvec,istati,icoeff,igamma,idialo,idialu,n_mag,nsmooth,imloss,ifitm,nndr,RSG_Mdot,iledou,idifcon,iover,my,&
            iunder,nbchx,nrband,islow,icncst,tauH_fit,iauto,iprn,iout,itmin,idebug,itests,n_snap
    real(kindreal), intent(out):: &
            binm2,periodini,zinit,zsol,z,fenerg,richac,frein,K_Kawaler,Omega_saturation,rapcrilim,vwant,xfom,omega,xdial,&
            B_initial,add_diff,alpha_F,fmlos,fitm,fitmi,deltal,deltat,Be_mdotfrac,start_mdot,elph,dovhp,dunder,gkorm,alph,&
            agdr,faktor,dgrp,dgrl,dgry,dgrc,dgro,dgr20,xcn
    character(len=200), intent(out):: &
            starname
    integer, intent(out) :: m
    real(kindreal), intent(out) :: &
            gms,alter,gls,teff,glsv,teffv,dzeitj,dzeit,dzeitv,xmini,summas,ab,&
            dm_lost
    real(kindreal), dimension(GenecStar%m), intent(out):: &
            q,p,t,r,s,x,y3,y,xc12,xc13,xn14,xn15,xo16,xo17,xo18,xne20,xne22,xmg24,xmg25,xmg26,xf19,xne21,xna23,xal27,xsi28,&
            xc14,xf18,xal26,xneut,xprot,omegi,xbid,xbid1,vp,vt,vr,vs,vx,vy,vy3,vxc12,vxc13,vxn14,vxn15,vxo16,vxo17,vxo18,&
            vxne20,vxne22,vxmg24,vxmg25,vxmg26,vxf19,vxne21,vxna23,vxal27,vxsi28,vxc14,vxf18,vxal26g,vxneut,vxprot,vomegi,&
            vxbid,vxbid1
    integer:: get_stellar_model

    initialised      = GenecStar%initialised
    starname         = GenecStar%starname
    nwseq            = GenecStar%nwseq            
    modanf           = GenecStar%modanf
    nzmod            = GenecStar%nzmod
    end_at_phase     = GenecStar%end_at_phase
    end_at_model     = GenecStar%end_at_model
    irot             = GenecStar%irot
    isol             = GenecStar%isol
    imagn            = GenecStar%imagn
    ialflu           = GenecStar%ialflu
    ianiso           = GenecStar%ianiso
    ipop3            = GenecStar%ipop3
    ibasnet          = GenecStar%ibasnet
    phase            = GenecStar%phase
    iprezams         = GenecStar%iprezams
    var_rates        = GenecStar%var_rates
    bintide          = GenecStar%bintide
    binm2            = GenecStar%binm2
    periodini        = GenecStar%periodini
    const_per        = GenecStar%const_per
    iopac            = GenecStar%iopac
    ikappa           = GenecStar%ikappa
    zinit            = GenecStar%zinit
    zsol             = GenecStar%zsol
    z                = GenecStar%z
    idiff            = GenecStar%idiff
    iadvec           = GenecStar%iadvec
    istati           = GenecStar%istati
    icoeff           = GenecStar%icoeff
    igamma           = GenecStar%igamma
    idialo           = GenecStar%idialo
    idialu           = GenecStar%idialu
    n_mag            = GenecStar%n_mag
    nsmooth          = GenecStar%nsmooth
    fenerg           = GenecStar%fenerg
    richac           = GenecStar%richac
    frein            = GenecStar%frein
    K_Kawaler        = GenecStar%K_Kawaler
    Omega_saturation = GenecStar%Omega_saturation
    rapcrilim        = GenecStar%rapcrilim
    vwant            = GenecStar%vwant
    xfom             = GenecStar%xfom
    omega            = GenecStar%omega
    xdial            = GenecStar%xdial
    B_initial        = GenecStar%B_initial
    add_diff         = GenecStar%add_diff
    alpha_F          = GenecStar%alpha_F
    Add_Flux         = GenecStar%Add_Flux
    diff_only        = GenecStar%diff_only
    qminsmooth       = GenecStar%qminsmooth
    imloss           = GenecStar%imloss
    ifitm            = GenecStar%ifitm
    nndr             = GenecStar%nndr
    RSG_Mdot         = GenecStar%RSG_Mdot
    fmlos            = GenecStar%fmlos
    fitm             = GenecStar%fitm
    fitmi            = GenecStar%fitmi
    deltal           = GenecStar%deltal
    deltat           = GenecStar%deltat
    Be_mdotfrac      = GenecStar%Be_mdotfrac
    start_mdot       = GenecStar%start_mdot
    SupraEddMdot     = GenecStar%SupraEddMdot
    iledou           = GenecStar%iledou
    idifcon          = GenecStar%idifcon
    my               = GenecStar%my
    iover            = GenecStar%iover
    iunder           = GenecStar%iunder
    elph             = GenecStar%elph
    dovhp            = GenecStar%dovhp
    dunder           = GenecStar%dunder
    nbchx            = GenecStar%nbchx
    nrband           = GenecStar%nrband
    gkorm            = GenecStar%gkorm
    alph             = GenecStar%alph
    agdr             = GenecStar%agdr
    faktor           = GenecStar%faktor
    dgrp             = GenecStar%dgrp
    dgrl             = GenecStar%dgrl
    dgry             = GenecStar%dgry
    dgrc             = GenecStar%dgrc
    dgro             = GenecStar%dgro
    dgr20            = GenecStar%dgr20
    islow            = GenecStar%islow
    icncst           = GenecStar%icncst
    tauH_fit         = GenecStar%tauH_fit
    xcn              = GenecStar%xcn
    iauto            = GenecStar%iauto
    iprn             = GenecStar%iprn
    iout             = GenecStar%iout
    itmin            = GenecStar%itmin
    idebug           = GenecStar%idebug
    itests           = GenecStar%itests
    n_snap           = GenecStar%n_snap
    display_plot     = GenecStar%display_plot
    xyfiles          = GenecStar%xyfiles
    verbose          = GenecStar%verbose
    stop_deg         = GenecStar%stop_deg
                                 
    m                = GenecStar%m
    gms              = GenecStar%gms
    alter            = GenecStar%alter
    gls              = GenecStar%gls
    teff             = GenecStar%teff
    glsv             = GenecStar%glsv
    teffv            = GenecStar%teffv
    dzeitj           = GenecStar%dzeitj
    dzeit            = GenecStar%dzeit
    dzeitv           = GenecStar%dzeitv
    summas           = GenecStar%summas
    xmini            = GenecStar%xmini
    ab               = GenecStar%ab
    dm_lost          = GenecStar%dm_lost
    q                = GenecStar%q
    p                = GenecStar%p
    t                = GenecStar%t
    r                = GenecStar%r
    s                = GenecStar%s
    x                = GenecStar%x
    y                = GenecStar%y
    xc12             = GenecStar%xc12
    vp               = GenecStar%vp
    vt               = GenecStar%vt
    vr               = GenecStar%vr
    vs               = GenecStar%vs
    xo16             = GenecStar%xo16
    vx               = GenecStar%vx
    vy               = GenecStar%vy
    vxc12            = GenecStar%vxc12
    vxo16            = GenecStar%vxo16
    y3               = GenecStar%y3
    xc13             = GenecStar%xc13
    xn14             = GenecStar%xn14
    xn15             = GenecStar%xn15
    xo17             = GenecStar%xo17
    xo18             = GenecStar%xo18
    vy3              = GenecStar%vy3
    vxc13            = GenecStar%vxc13
    vxn14            = GenecStar%vxn14
    vxn15            = GenecStar%vxn15
    vxo17            = GenecStar%vxo17
    vxo18            = GenecStar%vxo18
    xne20            = GenecStar%xne20
    xne22            = GenecStar%xne22
    xmg24            = GenecStar%xmg24
    xmg25            = GenecStar%xmg25
    xmg26            = GenecStar%xmg26
    vxne20           = GenecStar%vxne20
    vxne22           = GenecStar%vxne22
    vxmg24           = GenecStar%vxmg24
    vxmg25           = GenecStar%vxmg25
    vxmg26           = GenecStar%vxmg26
    omegi            = GenecStar%omegi
    vomegi           = GenecStar%vomegi
    xf19             = GenecStar%xf19
    xne21            = GenecStar%xne21
    xna23            = GenecStar%xna23
    xal26            = GenecStar%xal26
    xal27            = GenecStar%xal27
    xsi28            = GenecStar%xsi28
    vxf19            = GenecStar%vxf19
    vxne21           = GenecStar%vxne21
    vxna23           = GenecStar%vxna23
    vxal26g          = GenecStar%vxal26g
    vxal27           = GenecStar%vxal27
    vxsi28           = GenecStar%vxsi28
    xneut            = GenecStar%xneut
    xprot            = GenecStar%xprot
    xc14             = GenecStar%xc14
    xf18             = GenecStar%xf18
    xbid             = GenecStar%xbid
    xbid1            = GenecStar%xbid1
    vxneut           = GenecStar%vxneut
    vxprot           = GenecStar%vxprot
    vxc14            = GenecStar%vxc14
    vxf18            = GenecStar%vxf18
    vxbid            = GenecStar%vxbid
    vxbid1           = GenecStar%vxbid1
    get_stellar_model = 0
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

!function get_starname(index_of_the_star, am_starname)
!    implicit none
!    integer:: get_starname, index_of_the_star
!    character(len=200):: am_starname
!    am_starname = InitialGenecStar%starname
!    get_starname = 0
!end function
!
!function set_starname(index_of_the_star, am_starname)
!    implicit none
!    integer:: set_starname, index_of_the_star
!    character(len=200):: am_starname
!    InitialGenecStar%starname = am_starname
!    set_starname = 0
!end function

!function get_phase(index_of_the_star, phase)
!    implicit none
!    integer:: index_of_the_star, phase
!    integer:: get_phase
!    phase = GenecStar%phase
!    get_phase = 0
!end function
!
!function set_phase(index_of_the_star, phase)
!    implicit none
!    integer:: index_of_the_star, phase
!    integer:: set_phase
!    GenecStar%phase = phase
!    set_phase = 0
!end function

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
