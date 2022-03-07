function initialize_code()
    use WriteSaveClose, only: quitafterclosing
    use inputparam, only: amuseinterface, writetofiles
    use genec, only: initialise_genec
    use amuse_helpers, only: set_defaults
    use amuse_helpers, only: mstar, zini, starname
    use evol, only: input_dir
    implicit none
    integer:: initialize_code

    amuseinterface = .true.
    writetofiles = .false.
    input_dir = "./src/GENEC/code"
    call initialise_genec()
    call set_defaults()
    starname = "AmuseDefaultStar"
    mstar = 7.0
    zini = 0.014
    quitafterclosing = .false.    

    initialize_code = 0
end function

!function read_genec_model()
!    ! This should only be called if no star has been initialised yet!
!    use genec, only: amuseinterface
!    use amuse_helpers, only: starname
!    use inputparam, only:CharacteristicsParams,PhysicsParams,CompositionParams,RotationParams,&
!        SurfaceParams,ConvectionParams,ConvergenceParams,TimeControle,VariousSettings
!    implicit none
!    character(*):: inputfilename
!    integer:: read_genec_model
!
!    amuseinterface = .false. ! This will make initialise_star read the b file
!    write(inputfilename, *) starname
!    write(inputfilename, *) ".input"
!    open(unit=4242, file=inputfilename)
!    ! * Parse the CharacteristicsParams namelist *
!    read(*,nml=CharacteristicsParams)
!    ! * Parse the PhysicsParams namelist *
!    read(*,nml=PhysicsParams)
!    ! * Parse the CompositionParams namelist *
!    read(*,nml=CompositionParams)
!    ! * Parse the RotationParams namelist *
!    read(*,nml=RotationParams)
!    ! * Parse the SurfaceParams namelist *
!    read(*,nml=SurfaceParams)
!    ! * Parse the ConvectionParams namelist *
!    read(*,nml=ConvectionParams)
!    ! * Parse the ConvergenceParams namelist *
!    read(*,nml=ConvergenceParams)
!    ! * Parse the TimeControle namelist *
!    read(*,nml=TimeControle)
!    ! * Parse the VariousSettings namelist *
!    read(*,nml=VariousSettings)
!
!    close(4242)
!    read_genec_model = 0
!end function

function cleanup_code()
    implicit none
    integer:: cleanup_code
    cleanup_code = 0
end function

function commit_parameters()
    implicit none
    integer:: commit_parameters
    commit_parameters = 0
end function

function commit_particles()
    use genec, only: amuseinterface
    use amuse_helpers, only: initialise_star, makeini
    use WriteSaveClose, only: OpenAll
    use inputparam, only: nzmod
    implicit none
    integer:: commit_particles
    call makeini()  ! this will actually override some things from set_defaults now! FIXME
    !write(*,*) 'makeini done'
    call OpenAll()
    !write(*,*) 'OpenAll done'
    call initialise_star()
    nzmod = 1
    !write(*,*) 'initialise_star done'
    amuseinterface = .true. ! If we just read a b file, disable this again for continuing
    commit_particles = 0
end function

function delete_star()
    implicit none
    integer:: index_of_the_star
    integer:: delete_star
    delete_star = -1  ! not supported
end function

function evolve_model(end_time)
    use timestep, only: alter
    use WriteSaveClose, only: OpenAll
    use genec, only: evolve, modell, finalise, initialise_star
    implicit none
    double precision:: end_time
    integer:: evolve_model
    do while (alter < end_time)
        write(*,*) "Current time: ", alter, ", evolving to: ", end_time
        !modell = 1
        call evolve()
        call finalise()
        call OpenAll()
        call initialise_star()
        write(*,*) "*****Modell: ", modell
    end do
    evolve_model = 0
end function

function evolve_for(index_of_the_star, time)
    ! get current time
    ! set max time to current time plus argument
    ! evolve
    use timestep, only: alter
    use WriteSaveClose, only: OpenAll
    use genec, only: evolve, modell, finalise, initialise_star
    implicit none
    integer:: index_of_the_star
    double precision:: time, end_time
    integer:: evolve_for
    !if (alter == 0.) then
    !    call initialise_star()
    !endif
    end_time = alter+time
    do while (alter < end_time)
        write(*,*) "Current time: ", alter, ", evolving to: ", end_time
        !modell = 1
        call evolve()
        call finalise()
        call OpenAll()
        call initialise_star()
        write(*,*) "*****Modell: ", modell
    end do

    evolve_for = 0
end function

function evolve_one_step(index_of_the_star)
    use timestep, only: alter
    use WriteSaveClose, only: OpenAll
    use genec, only: evolve, modell, finalise, initialise_star
    use inputparam,only: modanf,nwseq,nzmod
    implicit none
    integer:: index_of_the_star
    integer:: evolve_one_step
    integer:: original_nzmod
    original_nzmod = nzmod
    nzmod = 1
    write(*,*) "Evolving one step, current time: ", alter
    !modell = 1
    call evolve()
    call finalise()
    call OpenAll()
    call initialise_star() ! will set modell to 1
    write(*,*) "Evolved one step, current time: ", alter
    nzmod = original_nzmod
    write(*,*) "*****modanf, nwseq, nzmod: ", modanf, nwseq, nzmod
    evolve_one_step = 0
end function

function get_age(index_of_the_star, age)
    use timestep, only: alter
    implicit none
    integer:: index_of_the_star
    double precision:: age
    integer:: get_age
    age = alter
    get_age = 0
end function

function get_density_at_zone(index_of_the_star, zone, rho_i)
    use strucmod, only: rho, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    double precision:: rho_i
    integer:: get_density_at_zone
    i = m - zone
    if (zone <= m) then
        rho_i = exp(rho(i))
    end if
    get_density_at_zone = 0
end function

function get_luminosity(index_of_the_star, luminosity)
    !use strucmod, only: s, m
    use caramodele, only: gls
    implicit none
    integer:: index_of_the_star
    double precision:: luminosity
    integer:: get_luminosity
    !luminosity = exp(s(m))  ! in cgs units, so erg/s?
    luminosity = gls
    get_luminosity = 0
end function

function get_luminosity_at_zone(index_of_the_star, zone, lum_i)
    use strucmod, only: s, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    double precision:: lum_i
    integer:: get_luminosity_at_zone
    i = m - zone
    if (zone <= m) then
        !lum_i = exp(s(zone+1))
        lum_i = exp(s(i)) - 1
    end if
    get_luminosity_at_zone = 0
end function

function get_mass_fraction_at_zone(index_of_the_star, zone, dq_i)
    use strucmod, only: q, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    double precision:: dq_i
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

function get_mass(index_of_the_star, mass)
    use inputparam, only: starname
    use caramodele, only: gms
    implicit none
    double precision:: mass
    integer:: get_mass, index_of_the_star
    mass = gms
    get_mass = 0
end function

function get_mass_of_species(index_of_the_star, species, species_mass)
    implicit none
    integer:: index_of_the_star
    integer:: species
    double precision:: species_mass
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
    double precision:: Xj_i
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

function get_metallicity(metallicity)
    implicit none
    double precision:: metallicity
    integer:: get_metallicity
    get_metallicity = 0
end function

function get_mu_at_zone(index_of_the_star, zone, mu_i)
    use strucmod, only: m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    integer:: get_mass_fraction_of_species_at_zone
    double precision:: mu_i, X, Y3, Y
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

    species_names(1) = 'h'
    species_names(2) = 'he3'
    species_names(3) = 'he'
    species_names(4) = 'c12'
    species_names(5) = 'c13'
    species_names(6) = 'n14'
    species_names(7) = 'n15'
    species_names(8) = 'o16'
    species_names(9) = 'o17'
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

function get_number_of_particles()
    implicit none
    integer:: get_number_of_particles
    get_number_of_particles = 0
end function

function get_number_of_species(index_of_the_star, n_species)
    use inputparam, only: ialflu
    implicit none
    integer:: index_of_the_star
    integer:: n_species
    integer:: get_number_of_species
    if (ialflu==1) then
        n_species = 27
    else
        n_species = 15
    end if
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
    use inputparam, only: starname
    use strucmod, only: m
    implicit none
    integer:: index_of_the_star
    integer:: n_zones
    integer:: get_number_of_zones
    !if (char(index_of_the_star) /= starname) then
    !    get_number_of_zones = -1
    !    return
    !end if
    n_zones = m
    get_number_of_zones = 0
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
    double precision:: P_i
    integer:: get_pressure_at_zone
    if (zone <= m) then
        i = m - zone
        P_i = exp(p(i))
    end if
    get_pressure_at_zone = 0
end function

function get_radius(index_of_the_star, radius)
    use strucmod, only: r, m
    implicit none
    integer:: index_of_the_star
    double precision:: radius
    integer:: get_radius
    radius = exp(r(1))  ! in cm
    get_radius = 0
end function

function get_radius_at_zone(index_of_the_star, zone, R_i)
    use strucmod, only: r, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    double precision:: R_i
    integer:: get_radius_at_zone
    i = m - zone
    if (zone <= m) then
        R_i = exp(r(i))  ! in cm
    end if
    get_radius_at_zone = 0
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
    double precision:: temperature
    integer:: get_temperature
    temperature = teff
    get_temperature = 0
end function

function get_temperature_at_zone(index_of_the_star, zone, T_i)
    use strucmod, only: t, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    double precision:: T_i
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
    double precision:: time_step
    integer:: get_time_step
    time_step = dzeitj
    get_time_step = 0
end function

function get_time(time)
    use timestep, only: alter
    implicit none
    double precision:: time
    integer:: get_time
    time = alter
    get_time = 0
end function

function new_particle(index_of_the_star, mass, metallicity, am_starname)
    use amuse_helpers, only: starname, mstar, zini
    implicit none
    integer:: index_of_the_star, key
    double precision:: mass, metallicity
    integer:: new_particle
    character(len=12):: am_starname
    starname = am_starname !'AmuseStar'!write(starname, '(i10.10)') key
    index_of_the_star = 1
    mstar = mass
    zini = metallicity
    
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
    call finalise()
    call OpenAll()
    call initialise_star()
    recommit_particles = 0
end function

function set_density_at_zone(index_of_the_star, zone, rho_i)
    use strucmod, only: rho, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    double precision:: rho_i
    integer:: set_density_at_zone
    i = m - zone
    if (zone <= m) then
        rho(i) = log(rho_i)
    end if
    set_density_at_zone = 0
end function

function set_genec_path(path)
    use evol, only: input_dir
    implicit none
    character(len = :), allocatable:: path
    integer:: set_genec_path

    input_dir = path
    set_genec_path = 0
end function

function set_starname(index_of_the_star)
    use amuse_helpers, only:starname
    implicit none
    integer:: set_starname, index_of_the_star
    ! Only allow this at the initialisation time!
    starname = 'AmuseStar'
    write(starname, '(i10.10)') index_of_the_star
    set_starname = 0
end function

function set_mass(index_of_the_star, mass)
    use amuse_helpers, only:mstar  ! this is initial mass only?
    implicit none
    integer:: set_mass, index_of_the_star
    double precision:: mass
    mstar = mass
    set_mass = 0
end function

function set_mass_fraction_of_species_at_zone(index_of_the_star, species, zone, Xj_i)
    use strucmod, only: m
    use abundmod, only: &
        x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24,&
        xmg25,xmg26,xal26,xal27,xsi28,xprot,xneut,xbid,xbid1
    implicit none
    integer:: index_of_the_star
    integer:: species, zone, i
    double precision:: Xj_i
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

function set_metallicity(metallicity)
    use amuse_helpers, only: zini
    implicit none
    double precision:: metallicity
    integer:: set_metallicity
    zini = metallicity
    set_metallicity = 0
end function

function set_radius_at_zone(index_of_the_star, zone, R_i)
    use strucmod, only: r, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    double precision:: R_i
    integer:: set_radius_at_zone
    i = m - zone
    if (zone <= m) then
        r(i) = log(R_i)
    end if
    set_radius_at_zone = 0
end function

function set_temperature_at_zone(index_of_the_star, zone, T_i)
    use strucmod, only: t, m
    implicit none
    integer:: index_of_the_star
    integer:: zone, i
    double precision:: T_i
    integer:: set_temperature_at_zone
    i = m - zone
    if (zone <= m) then
        t(i) = log(T_i)
    end if
    
    set_temperature_at_zone = 0
end function
