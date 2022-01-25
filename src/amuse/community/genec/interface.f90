function initialize_code()
    use WriteSaveClose, only: quitafterclosing
    use genec, only: initialise_genec, amuseinterface
    use amuse_helpers, only: set_defaults
    use amuse_helpers, only: mstar, zini, starname
    use evol, only: input_dir
    implicit none
    integer:: initialize_code

    amuseinterface = .true.
    input_dir = "./src/GENEC/code"
    call initialise_genec()
    call set_defaults()
    starname = "AmuseDefaultStar"
    mstar = 7.0
    zini = 0.014
    quitafterclosing = .false.    

    initialize_code = 0
end function

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
    use amuse_helpers, only: initialise_star, makeini
    use WriteSaveClose, only: OpenAll
    implicit none
    integer:: commit_particles
    call makeini()
    !write(*,*) 'makeini done'
    call OpenAll()
    !write(*,*) 'OpenAll done'
    call initialise_star()
    !write(*,*) 'initialise_star done'
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
        modell = 1
        call evolve()
        call finalise()
        call OpenAll()
        call initialise_star()
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
        modell = 1
        call evolve()
        call finalise()
        call OpenAll()
        call initialise_star()
    end do

    evolve_for = 0
end function

function evolve_one_step(index_of_the_star)
    use timestep, only: alter
    use WriteSaveClose, only: OpenAll
    use genec, only: evolve, modell, finalise, initialise_star
    implicit none
    integer:: index_of_the_star
    integer:: evolve_one_step
    write(*,*) "Evolving one step, current time: ", alter
    modell = 1
    call evolve()
    call finalise()
    call OpenAll()
    call initialise_star()
    write(*,*) "Evolved one step, current time: ", alter
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
    integer:: zone
    double precision:: rho_i
    integer:: get_density_at_zone
    if (zone <= m) then
        rho_i = exp(rho(m-zone))
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
        lum_i = exp(s(m-zone))
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

function get_mass_fraction_of_species_at_zone(index_of_the_star, species, zone, Xj_i)
    use abundmod, only: &
        x,y3,y,xc12,xc13,xc14,xn14,xn15,xo16,xo17,xo18,xf18,xf19,xne20,xne21,xne22,xna23,xmg24,&
        xmg25,xmg26,xal26,xal27,xsi28,xprot,xneut,xbid,xbid1
    implicit none
    integer:: index_of_the_star
    integer:: species, zone
    double precision:: Xj_i
    integer:: get_mass_fraction_of_species_at_zone
    select case(species)
    case(1)
        Xj_i = x(zone)
    case(2)
        Xj_i = y3(zone)
    case(3)
        Xj_i = y(zone)
    case(4)
        Xj_i = xc12(zone)
    case(5)
        Xj_i = xc13(zone)
    case default
        Xj_i = 0
    end select

    get_mass_fraction_of_species_at_zone = 0
end function

function get_metallicity(metallicity)
    implicit none
    double precision:: metallicity
    integer:: get_metallicity
    get_metallicity = 0
end function

function get_mu_at_zone(index_of_the_star, zone, mu_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone
    double precision:: mu_i
    integer:: get_mu_at_zone
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
    use inputparam,only: ialflu
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
    integer:: zone
    double precision:: R_i
    integer:: get_radius_at_zone
    if (zone <= m) then
        !R_i = exp(r(zone+1))  ! in cm
        R_i = exp(r(m-zone))  ! in cm
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
    integer:: zone
    double precision:: T_i
    integer:: get_temperature_at_zone
    if (zone <= m) then
        T_i = exp(t(m-zone))
    else
        get_temperature_at_zone = -2
        return
    end if
    get_temperature_at_zone = 0
end function

function get_time_step(index_of_the_star, time_step)
    use timestep, only: dzeit
    implicit none
    integer:: index_of_the_star
    double precision:: time_step
    integer:: get_time_step
    time_step = dzeit
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

function new_particle(index_of_the_star, mass, metallicity, key)
    use amuse_helpers, only: starname, mstar, zini
    implicit none
    integer:: index_of_the_star, key
    double precision:: mass, metallicity
    integer:: new_particle
    starname = 'AmuseStar'!write(starname, '(i10.10)') key
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
    implicit none
    integer:: index_of_the_star
    integer:: zone
    double precision:: rho_i
    integer:: set_density_at_zone
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
    implicit none
    integer:: index_of_the_star
    integer:: species, zone
    double precision:: Xj_i
    integer:: set_mass_fraction_of_species_at_zone
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
    implicit none
    integer:: index_of_the_star
    integer:: zone
    double precision:: R_i
    integer:: set_radius_at_zone
    set_radius_at_zone = 0
end function

function set_temperature_at_zone(index_of_the_star, zone, T_i)
    implicit none
    integer:: index_of_the_star
    integer:: zone
    double precision:: T_i
    integer:: set_temperature_at_zone
    set_temperature_at_zone = 0
end function
