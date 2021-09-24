!MODULE PhantomInterface
!    
!CONTAINS

function initialize_code()
  use StoppingConditions
  implicit none
  integer :: initialize_code
  integer :: error
  double precision :: polyk
  call amuse_initialize_code()

  !error = set_support_for_condition(TIMEOUT_DETECTION)
  !error = set_support_for_condition(NUMBER_OF_STEPS_DETECTION)
  !error = set_support_for_condition(OUT_OF_BOX_DETECTION)
  error = set_support_for_condition(DENSITY_LIMIT_DETECTION)
  !error = set_support_for_condition(INTERNAL_ENERGY_LIMIT_DETECTION)
  initialize_code=0
end function

function cleanup_code()
  implicit none
  integer :: cleanup_code
  call amuse_cleanup_code()
  cleanup_code=0
end function

function commit_particles()
  implicit none
  integer :: commit_particles
  call amuse_commit_particles()
  commit_particles=0
end function

function get_time(time)
  implicit none
  double precision :: time
  integer :: get_time
  call amuse_get_time(time)
  get_time=0
end function

function get_mass(index_of_the_particle, mass)
  implicit none
  integer :: index_of_the_particle
  double precision :: mass
  integer :: get_mass
  call amuse_get_mass(index_of_the_particle, mass)
  get_mass=0
end function

function set_mass(index_of_the_particle, mass)
  implicit none
  integer :: index_of_the_particle
  double precision :: mass
  integer :: set_mass
  call amuse_set_mass(index_of_the_particle, mass)
  set_mass=0
end function

function set_smoothing_length(index_of_the_particle, h_smooth)
  implicit none
  integer :: index_of_the_particle
  double precision :: h_smooth
  integer :: set_smoothing_length
  call amuse_set_smoothing_length(index_of_the_particle, h_smooth)
  set_smoothing_length=0
end function

function set_internal_energy(index_of_the_particle, u)
  implicit none
  integer :: index_of_the_particle
  double precision :: u
  integer :: set_internal_energy
  call amuse_set_internal_energy(index_of_the_particle, u)
  set_internal_energy=0
end function

function set_h2ratio(index_of_the_particle, h2ratio)
  implicit none
  integer :: index_of_the_particle
  double precision :: h2ratio
  integer :: set_h2ratio
  call amuse_set_h2ratio(index_of_the_particle, h2ratio)
  set_h2ratio=0
end function

function set_hi_abundance(index_of_the_particle, hi_abundance)
  implicit none
  integer :: index_of_the_particle
  double precision :: hi_abundance
  integer :: set_hi_abundance
  call amuse_set_hi_abundance(index_of_the_particle, hi_abundance)
  set_hi_abundance=0
end function

function set_proton_abundance(index_of_the_particle, proton_abundance)
  implicit none
  integer :: index_of_the_particle
  double precision :: proton_abundance
  integer :: set_proton_abundance
  call amuse_set_proton_abundance(index_of_the_particle, proton_abundance)
  set_proton_abundance=0
end function

function set_electron_abundance(index_of_the_particle, electron_abundance)
  implicit none
  integer :: index_of_the_particle
  double precision :: electron_abundance
  integer :: set_electron_abundance
  call amuse_set_electron_abundance(index_of_the_particle, electron_abundance)
  set_electron_abundance=0
end function

function set_co_abundance(index_of_the_particle, co_abundance)
  implicit none
  integer :: index_of_the_particle
  double precision :: co_abundance
  integer :: set_co_abundance
  call amuse_set_co_abundance(index_of_the_particle, co_abundance)
  set_co_abundance=0
end function

function get_h2ratio(index_of_the_particle, h2ratio)
  implicit none
  integer :: index_of_the_particle
  double precision :: h2ratio
  integer :: get_h2ratio
  call amuse_get_h2ratio(index_of_the_particle, h2ratio)
  get_h2ratio=0
end function

function get_hi_abundance(index_of_the_particle, hi_abundance)
  implicit none
  integer :: index_of_the_particle
  double precision :: hi_abundance
  integer :: get_hi_abundance
  call amuse_get_hi_abundance(index_of_the_particle, hi_abundance)
  get_hi_abundance=0
end function

function get_proton_abundance(index_of_the_particle, proton_abundance)
  implicit none
  integer :: index_of_the_particle
  double precision :: proton_abundance
  integer :: get_proton_abundance
  call amuse_get_proton_abundance(index_of_the_particle, proton_abundance)
  get_proton_abundance=0
end function

function get_electron_abundance(index_of_the_particle, electron_abundance)
  implicit none
  integer :: index_of_the_particle
  double precision :: electron_abundance
  integer :: get_electron_abundance
  call amuse_get_electron_abundance(index_of_the_particle, electron_abundance)
  get_electron_abundance=0
end function

function get_co_abundance(index_of_the_particle, co_abundance)
  implicit none
  integer :: index_of_the_particle
  double precision :: co_abundance
  integer :: get_co_abundance
  call amuse_get_co_abundance(index_of_the_particle, co_abundance)
  get_co_abundance=0
end function

function get_pressure(index_of_the_particle, p)
  implicit none
  integer :: index_of_the_particle
  double precision :: p
  integer :: get_pressure
  call amuse_get_pressure(index_of_the_particle, p)
  get_pressure=0
end function

function get_density(index_of_the_particle, rho)
  implicit none
  integer :: index_of_the_particle
  double precision :: rho
  integer :: get_density
  call amuse_get_density(index_of_the_particle, rho)
  get_density=0
end function

function get_smoothing_length(index_of_the_particle, h_smooth)
  implicit none
  integer :: index_of_the_particle
  double precision :: h_smooth
  integer :: get_smoothing_length
  call amuse_get_smoothing_length(index_of_the_particle, h_smooth)
  get_smoothing_length=0
end function

function get_internal_energy(index_of_the_particle, u)
  implicit none
  integer :: index_of_the_particle
  double precision :: u
  integer :: get_internal_energy
  call amuse_get_internal_energy(index_of_the_particle, u)  
  get_internal_energy=0
end function

function get_index_of_first_particle(index_of_the_particle)
  implicit none
  integer :: index_of_the_particle
  integer :: get_index_of_first_particle
  get_index_of_first_particle=0
end function

function get_total_radius(radius)
  implicit none
  double precision :: radius
  integer :: get_total_radius
  get_total_radius=0
end function

function get_potential_at_point(eps, x, y, z, phi, npoints)
  implicit none
  integer :: npoints
  double precision :: eps, x, y, z, phi
  integer :: get_potential_at_point
  get_potential_at_point=0
end function

function get_total_mass(mass)
  implicit none
  double precision :: mass
  integer :: get_total_mass
  get_total_mass=0
end function

function evolve_model(tmax)
  use StoppingConditions
  implicit none
  double precision :: tmax
  integer :: evolve_model
  integer :: sc
  integer :: i, nmax
  integer :: is_density_limit_detection_enabled, stopping_index
  integer :: error
  double precision :: minimum_density_parameter, maximum_density_parameter, rho, radius

  error = reset_stopping_conditions()
  error = is_stopping_condition_enabled(&
      DENSITY_LIMIT_DETECTION, is_density_limit_detection_enabled)
  error = get_stopping_condition_minimum_density_parameter(minimum_density_parameter)
  error = get_stopping_condition_maximum_density_parameter(maximum_density_parameter)

  call amuse_evolve_model(tmax)
  if (is_density_limit_detection_enabled > 0) then
      call amuse_get_number_of_sph_particles(nmax)
      do i=1, nmax
          call amuse_get_radius(i, radius)
          if (radius > 0) then
              call amuse_get_density(i, rho)
              if (&
                  (rho > maximum_density_parameter) .or. &
                  (rho < minimum_density_parameter) &
                  ) then
                  stopping_index = next_index_for_stopping_condition()
                  if (stopping_index > 0) then
                      error = set_stopping_condition_info(stopping_index, DENSITY_LIMIT_DETECTION)
                      error = set_stopping_condition_particle_index(stopping_index, 0, i)
                  endif
              endif
          endif
      enddo
  endif
  evolve_model=0
end function

function set_state_sph(index_of_the_particle, mass, x, y, z, &
        vx, vy, vz, u, h_smooth)
  implicit none
  integer :: index_of_the_particle
  double precision :: mass, x, y, z, vx, vy, vz, u, h_smooth
  integer :: set_state_sph
  call amuse_set_state_gas(index_of_the_particle, mass, x, y, z, &
      vx, vy, vz, u, h_smooth)
  set_state_sph=0
end function

function set_eps2(epsilon_squared)
  implicit none
  double precision :: epsilon_squared
  integer :: set_eps2
  set_eps2=-1
end function

function set_state_star(index_of_the_particle, mass, x, y, z, &
        vx, vy, vz, tform, radius)
  implicit none
  integer :: index_of_the_particle
  double precision :: mass, x, y, z, vx, vy, vz, tform, radius
  integer :: set_state_star
  set_state_star=-1
end function

function get_begin_time(time)
  implicit none
  double precision :: time
  integer :: get_begin_time
  get_begin_time=-1
end function

function get_eps2(epsilon_squared)
  implicit none
  double precision :: epsilon_squared
  integer :: get_eps2
  get_eps2=-1
end function

function get_index_of_next_particle(index_of_the_particle,  &
    index_of_the_next_particle)
  implicit none
  integer :: index_of_the_particle, index_of_the_next_particle
  integer :: get_index_of_next_particle
  get_index_of_next_particle=-1
end function

function new_sph_particle(index_of_the_particle, mass, x, y, z, &
        vx, vy, vz, u, h_smooth)
  implicit none
  integer :: index_of_the_particle
  double precision :: mass, x, y, z, vx, vy, vz, u, h_smooth
  integer :: new_sph_particle
  call amuse_new_sph_particle(index_of_the_particle, mass, x, y, z, &
      vx, vy, vz, u, h_smooth)
  new_sph_particle=0
end function

function delete_particle(index_of_the_particle)
  implicit none
  integer :: index_of_the_particle
  integer :: delete_particle
  call amuse_delete_particle(index_of_the_particle)
  delete_particle=0
end function

function get_potential(index_of_the_particle, potential)
  implicit none
  integer :: index_of_the_particle
  double precision :: potential
  integer :: get_potential
  get_potential=0
end function

function synchronize_model()
  implicit none
  integer :: synchronize_model
  synchronize_model=0
end function

function set_state_sink(index_of_the_particle, mass, x, y, z, &
        vx, vy, vz, radius, h_smooth)
  implicit none
  integer :: index_of_the_particle
  double precision :: mass, x, y, z, vx, vy, vz, radius, h_smooth
  integer :: set_state_sink
  call amuse_set_state_sink(index_of_the_particle, mass, x, y, z, &
      vx, vy, vz, radius, h_smooth)
  set_state_sink=0
end function

function get_state_sink(index_of_the_particle, mass, x, y, z, &
        vx, vy, vz, radius, h_smooth)
  implicit none
  integer :: index_of_the_particle
  double precision :: mass, x, y, z, vx, vy, vz, radius, h_smooth
  integer :: get_state_sink
  call amuse_get_state_sink(index_of_the_particle, mass, x, y, z, &
      vx, vy, vz, radius, h_smooth)
  get_state_sink=0
end function

function set_state_dm(index_of_the_particle, mass, x, y, z, &
        vx, vy, vz)
  implicit none
  integer :: index_of_the_particle
  double precision :: mass, x, y, z, vx, vy, vz
  integer :: set_state_dm
  call amuse_set_state_dm(index_of_the_particle, mass, x, y, z, &
      vx, vy, vz)
  set_state_dm=0
end function

function get_state_dm(index_of_the_particle, mass, x, y, z, &
        vx, vy, vz)
  implicit none
  integer :: index_of_the_particle
  double precision :: mass, x, y, z, vx, vy, vz
  integer :: get_state_dm
  call amuse_get_state_dm(index_of_the_particle, mass, x, y, z, &
      vx, vy, vz)
  get_state_dm=0
end function

function get_time_step(time_step)
  implicit none
  double precision :: time_step
  integer :: get_time_step
  call amuse_get_time_step(time_step)
  get_time_step=0
end function

function set_time_step(time_step)
  implicit none
  double precision :: time_step
  integer :: set_time_step
  call amuse_set_time_step(time_step)
  set_time_step=0
end function

function recommit_particles()
  implicit none
  integer :: recommit_particles
  call amuse_recommit_particles()
  recommit_particles=0
end function

function get_kinetic_energy(kinetic_energy)
  implicit none
  double precision :: kinetic_energy
  integer :: get_kinetic_energy
  call amuse_get_kinetic_energy(kinetic_energy)
  get_kinetic_energy=0
end function

function get_thermal_energy(thermal_energy)
  implicit none
  double precision :: thermal_energy
  integer :: get_thermal_energy
  call amuse_get_thermal_energy(thermal_energy)
  get_thermal_energy=0
end function

function get_number_of_particles(n)
  implicit none
  integer :: n
  integer :: get_number_of_particles
  call amuse_get_number_of_particles(n)
  get_number_of_particles=0
end function

function set_acceleration(index_of_the_particle, ax, ay, az)
  implicit none
  integer :: index_of_the_particle
  double precision :: ax, ay, az
  integer :: set_acceleration
  set_acceleration=-1
end function

function get_center_of_mass_position(x, y, z)
  implicit none
  double precision :: x, y, z
  integer :: get_center_of_mass_position
  get_center_of_mass_position=-1
end function

function get_center_of_mass_velocity(vx, vy, vz)
  implicit none
  double precision :: vx, vy, vz
  integer :: get_center_of_mass_velocity
  get_center_of_mass_velocity=-1
end function

function get_radius(index_of_the_particle, radius)
  implicit none
  integer :: index_of_the_particle
  double precision :: radius
  integer :: get_radius
  call amuse_get_radius(index_of_the_particle, radius)
  get_radius=0
end function

function get_state_star(index_of_the_particle, mass, x, y, z, vx, vy, vz,  &
    tform, radius)
  implicit none
  integer :: index_of_the_particle
  double precision :: mass, x, y, z, vx, vy, vz, tform, radius
  integer :: get_state_star
  get_state_star=-1
end function

function set_begin_time(time)
  implicit none
  double precision :: time
  integer :: set_begin_time
  set_begin_time=-1
end function

function new_dm_particle(index_of_the_particle, mass, x, y, z, vx, vy, vz)
  implicit none
  integer :: index_of_the_particle
  double precision :: mass, x, y, z, vx, vy, vz, radius
  integer :: new_dm_particle
  call amuse_new_dm_particle(index_of_the_particle, mass, x, y, z, &
      vx, vy, vz, radius)
  new_dm_particle=0
end function

function new_sink_particle(index_of_the_particle, mass, x, y, z, vx, vy, vz, &
        radius, h_smooth)
  implicit none
  integer :: index_of_the_particle
  double precision :: mass, x, y, z, vx, vy, vz, radius, h_smooth
  integer :: new_sink_particle
  call amuse_new_sink_particle(index_of_the_particle, mass, x, y, z, &
      vx, vy, vz, radius, h_smooth)
  new_sink_particle=0
end function

function set_radius(index_of_the_particle, radius)
  implicit none
  integer :: index_of_the_particle
  double precision :: radius
  integer :: set_radius
  call amuse_set_radius(index_of_the_particle, radius)
  set_radius=0
end function

function recommit_parameters()
  implicit none
  integer :: recommit_parameters
  recommit_parameters=0
end function

function get_potential_energy(potential_energy)
  implicit none
  double precision :: potential_energy
  integer :: get_potential_energy
  call amuse_get_potential_energy(potential_energy)
  get_potential_energy=0
end function

function get_state_sph(index_of_the_particle, mass, x, y, z, vx, vy, vz, u,  &
    h_smooth)
  implicit none
  integer :: index_of_the_particle
  double precision :: mass, x, y, z, vx, vy, vz, u, h_smooth
  integer :: get_state_sph
  call amuse_get_state_gas(index_of_the_particle, mass, x, y, z, vx, vy, vz, u, h_smooth)
  get_state_sph=0
end function

function get_gravity_at_point(eps, x, y, z, ax, ay, az, npoints)
  implicit none
  integer :: npoints
  double precision :: eps, x, y, z, ax, ay, az
  integer :: get_gravity_at_point
  get_gravity_at_point=-1
end function

function get_velocity(index_of_the_particle, vx, vy, vz)
  implicit none
  integer :: index_of_the_particle
  double precision :: vx, vy, vz
  integer :: get_velocity
  call amuse_get_velocity(index_of_the_particle, vx, vy, vz)
  get_velocity=0
end function

function get_acceleration(index_of_the_particle, ax, ay, az)
  implicit none
  integer :: index_of_the_particle
  double precision :: ax, ay, az
  integer :: get_acceleration
  call amuse_get_acceleration(index_of_the_particle, ax, ay, az)
  get_acceleration=0
end function

function new_star_particle(index_of_the_particle, mass, x, y, z, vx, vy,  &
    vz, tform, radius)
  implicit none
  integer :: index_of_the_particle
  double precision :: mass, x, y, z, vx, vy, vz, tform, radius
  integer :: new_star_particle
  new_star_particle=-1
end function

function get_position(index_of_the_particle, x, y, z)
  implicit none
  integer :: index_of_the_particle, i
  double precision :: x, y, z
  integer :: get_position
  call amuse_get_position(index_of_the_particle, x, y, z)
  get_position=0
end function

function set_position(index_of_the_particle, x, y, z)
  implicit none
  integer :: index_of_the_particle
  double precision :: x, y, z
  integer :: set_position
  call amuse_set_position(index_of_the_particle, x, y, z)
  set_position=0
end function

function commit_parameters()
  implicit none
  integer :: commit_parameters
  commit_parameters=0
end function

function set_velocity(index_of_the_particle, vx, vy, vz)
    implicit none
    integer :: index_of_the_particle
    double precision :: vx, vy, vz
    integer :: set_velocity
    call amuse_set_velocity(index_of_the_particle, vx, vy, vz)
    set_velocity=0
end function

function set_c_courant(C_cour)
    implicit none
    double precision :: C_cour
    integer :: set_c_courant
    call amuse_set_c_courant(C_cour)
    set_c_courant=0
end function

function set_c_force(C_force)
    implicit none
    double precision :: C_force
    integer :: set_c_force
    call amuse_set_c_force(C_force)
    set_c_force=0
end function

function set_c_cool(C_cool)
    implicit none
    double precision :: C_cool
    integer :: set_c_cool
    call amuse_set_c_cool(C_cool)
    set_c_cool=0
end function

function set_tolv(tolv)
    implicit none
    double precision :: tolv
    integer :: set_tolv
    call amuse_set_tolv(tolv)
    set_tolv=0
end function

function set_hfact(hfact)
    implicit none
    double precision :: hfact
    integer :: set_hfact
    call amuse_set_hfact(hfact)
    set_hfact=0
end function

function set_tolh(tolh)
    implicit none
    double precision :: tolh
    integer :: set_tolh
    call amuse_set_tolh(tolh)
    set_tolh=0
end function

function set_tree_accuracy(tree_accuracy)
    implicit none
    double precision :: tree_accuracy
    integer :: set_tree_accuracy
    call amuse_set_tree_accuracy(tree_accuracy)
    set_tree_accuracy=0
end function

function set_alpha(alpha)
    implicit none
    double precision :: alpha
    integer :: set_alpha
    call amuse_set_alpha(alpha)
    set_alpha=0
end function

function set_alphamax(alphamax)
    implicit none
    double precision :: alphamax
    integer :: set_alphamax
    call amuse_set_alphamax(alphamax)
    set_alphamax=0
end function

function set_beta(beta)
    implicit none
    double precision :: beta
    integer :: set_beta
    call amuse_set_beta(beta)
    set_beta=0
end function

function set_avdecayconst(avdecayconst)
    implicit none
    double precision :: avdecayconst
    integer :: set_avdecayconst
    call amuse_set_avdecayconst(avdecayconst)
    set_avdecayconst=0
end function

function set_idamp(idamp)
    implicit none
    integer :: idamp
    integer :: set_idamp
    call amuse_set_idamp(idamp)
    set_idamp=0
end function

function set_ieos(ieos)
    implicit none
    integer :: ieos
    integer :: set_ieos
    call amuse_set_ieos(ieos)
    set_ieos=0
end function

function set_icooling(icooling)
    implicit none
    integer :: icooling
    integer :: set_icooling
    call amuse_set_icooling(icooling)
    set_icooling=0
end function

function set_polyk(polyk)
    implicit none
    double precision :: polyk
    integer :: set_polyk
    call amuse_set_polyk(polyk)
    set_polyk=0
end function

function set_mu(mu)
    implicit none
    double precision :: mu
    integer :: set_mu
    call amuse_set_mu(mu)
    set_mu=0
end function

function set_rhofinal(rhofinal)
    implicit none
    double precision :: rhofinal
    integer :: set_rhofinal
    call amuse_set_rhofinal(rhofinal)
    set_rhofinal=0
end function

function set_rho_crit(rho_crit)
    implicit none
    double precision :: rho_crit
    integer :: set_rho_crit
    call amuse_set_rho_crit(rho_crit)
    set_rho_crit=0
end function

function set_r_crit(r_crit)
    implicit none
    double precision :: r_crit
    integer :: set_r_crit
    call amuse_set_r_crit(r_crit)
    set_r_crit=0
end function

function set_h_acc(h_acc)
    implicit none
    double precision :: h_acc
    integer :: set_h_acc
    call amuse_set_h_acc(h_acc)
    set_h_acc=0
end function

function set_h_soft_sinkgas(h_soft_sinkgas)
    implicit none
    double precision :: h_soft_sinkgas
    integer :: set_h_soft_sinkgas
    call amuse_set_h_soft_sinkgas(h_soft_sinkgas)
    set_h_soft_sinkgas=0
end function

function set_h_soft_sinksink(h_soft_sinksink)
    implicit none
    double precision :: h_soft_sinksink
    integer :: set_h_soft_sinksink
    call amuse_set_h_soft_sinksink(h_soft_sinksink)
    set_h_soft_sinksink=0
end function

function set_f_acc(f_acc)
    implicit none
    double precision :: f_acc
    integer :: set_f_acc
    call amuse_set_f_acc(f_acc)
    set_f_acc=0
end function

function set_iexternalforce(iexternalforce)
    implicit none
    integer :: iexternalforce
    integer :: set_iexternalforce
    call amuse_set_iexternalforce(iexternalforce)
    set_iexternalforce=0
end function

function set_irealvisc(irealvisc)
    implicit none
    integer :: irealvisc
    integer :: set_irealvisc
    call amuse_set_irealvisc(irealvisc)
    set_irealvisc=0
end function

function set_shearparam(shearparam)
    implicit none
    double precision :: shearparam
    integer :: set_shearparam
    call amuse_set_shearparam(shearparam)
    set_shearparam=0
end function

function set_bulkvisc(bulkvisc)
    implicit none
    double precision :: bulkvisc
    integer :: set_bulkvisc
    call amuse_set_bulkvisc(bulkvisc)
    set_bulkvisc=0
end function

function set_gamma(gamma)
    implicit none
    double precision :: gamma
    integer :: set_gamma
    call amuse_set_gamma(gamma)
    set_gamma=0
end function

function get_c_courant(C_cour)
  implicit none
  double precision :: C_cour
  integer :: get_c_courant
  call amuse_get_c_courant(C_cour)
  get_c_courant=0
end function

function get_c_force(C_force)
  implicit none
  double precision :: C_force
  integer :: get_c_force
  call amuse_get_c_force(C_force)
  get_c_force=0
end function

function get_c_cool(C_cool)
  implicit none
  double precision :: C_cool
  integer :: get_c_cool
  call amuse_get_c_cool(C_cool)
  get_c_cool=0
end function

function get_tolv(tolv)
    implicit none
    double precision :: tolv
    integer :: get_tolv
    call amuse_get_tolv(tolv)
    get_tolv=0
end function

function get_hfact(hfact)
    implicit none
    double precision :: hfact
    integer :: get_hfact
    call amuse_get_hfact(hfact)
    get_hfact=0
end function

function get_tolh(tolh)
    implicit none
    double precision :: tolh
    integer :: get_tolh
    call amuse_get_tolh(tolh)
    get_tolh=0
end function

function get_tree_accuracy(tree_accuracy)
    implicit none
    double precision :: tree_accuracy
    integer :: get_tree_accuracy
    call amuse_get_tree_accuracy(tree_accuracy)
    get_tree_accuracy=0
end function

function get_alpha(alpha)
    implicit none
    double precision :: alpha
    integer :: get_alpha
    call amuse_get_alpha(alpha)
    get_alpha=0
end function

function get_alphamax(alphamax)
    implicit none
    double precision :: alphamax
    integer :: get_alphamax
    call amuse_get_alphamax(alphamax)
    get_alphamax=0
end function

function get_beta(beta)
    implicit none
    double precision :: beta
    integer :: get_beta
    call amuse_get_beta(beta)
    get_beta=0
end function

function get_avdecayconst(avdecayconst)
    implicit none
    double precision :: avdecayconst
    integer :: get_avdecayconst
    call amuse_get_avdecayconst(avdecayconst)
    get_avdecayconst=0
end function

function get_idamp(idamp)
    implicit none
    integer :: idamp
    integer :: get_idamp
    call amuse_get_idamp(idamp)
    get_idamp=0
end function

function get_ieos(ieos)
    implicit none
    integer :: ieos
    integer :: get_ieos
    call amuse_get_ieos(ieos)
    get_ieos=0
end function

function get_icooling(icooling)
    implicit none
    integer :: icooling
    integer :: get_icooling
    call amuse_get_icooling(icooling)
    get_icooling=0
end function

function get_polyk(polyk)
    implicit none
    double precision :: polyk
    integer :: get_polyk
    call amuse_get_polyk(polyk)
    get_polyk=0
end function

function get_mu(mu)
    implicit none
    double precision :: mu
    integer :: get_mu
    call amuse_get_mu(mu)
    get_mu=0
end function

function get_rhofinal(rhofinal)
    implicit none
    double precision :: rhofinal
    integer :: get_rhofinal
    call amuse_get_rhofinal(rhofinal)
    get_rhofinal=0
end function

function get_rho_crit(rho_crit)
    implicit none
    double precision :: rho_crit
    integer :: get_rho_crit
    call amuse_get_rho_crit(rho_crit)
    get_rho_crit=0
end function

function get_r_crit(r_crit)
    implicit none
    double precision :: r_crit
    integer :: get_r_crit
    call amuse_get_r_crit(r_crit)
    get_r_crit=0
end function

function get_h_acc(h_acc)
    implicit none
    double precision :: h_acc
    integer :: get_h_acc
    call amuse_get_h_acc(h_acc)
    get_h_acc=0
end function

function get_h_soft_sinkgas(h_soft_sinkgas)
    implicit none
    double precision :: h_soft_sinkgas
    integer :: get_h_soft_sinkgas
    call amuse_get_h_soft_sinkgas(h_soft_sinkgas)
    get_h_soft_sinkgas=0
end function

function get_h_soft_sinksink(h_soft_sinksink)
    implicit none
    double precision :: h_soft_sinksink
    integer :: get_h_soft_sinksink
    call amuse_get_h_soft_sinksink(h_soft_sinksink)
    get_h_soft_sinksink=0
end function

function get_f_acc(f_acc)
    implicit none
    double precision :: f_acc
    integer :: get_f_acc
    call amuse_get_f_acc(f_acc)
    get_f_acc=0
end function

function get_iexternalforce(iexternalforce)
    implicit none
    integer :: iexternalforce
    integer :: get_iexternalforce
    call amuse_get_iexternalforce(iexternalforce)
    get_iexternalforce=0
end function

function get_irealvisc(irealvisc)
    implicit none
    integer :: irealvisc
    integer :: get_irealvisc
    call amuse_get_irealvisc(irealvisc)
    get_irealvisc=0
end function

function get_shearparam(shearparam)
    implicit none
    double precision :: shearparam
    integer :: get_shearparam
    call amuse_get_shearparam(shearparam)
    get_shearparam=0
end function

function get_bulkvisc(bulkvisc)
    implicit none
    double precision :: bulkvisc
    integer :: get_bulkvisc
    call amuse_get_bulkvisc(bulkvisc)
    get_bulkvisc=0
end function

function get_gamma(gamma)
    implicit none
    double precision :: gamma
    integer :: get_gamma
    call amuse_get_gamma(gamma)
    get_gamma=0
end function

function get_unit_length(unit_length)
    implicit none
    double precision :: unit_length
    integer :: get_unit_length
    call amuse_get_unit_length(unit_length)
    get_unit_length=0
end function

function set_unit_length(unit_length)
    implicit none
    double precision :: unit_length
    integer :: set_unit_length
    call amuse_set_unit_length(unit_length)
    set_unit_length=0
end function

function get_unit_mass(unit_mass)
    implicit none
    double precision :: unit_mass
    integer :: get_unit_mass
    call amuse_get_unit_mass(unit_mass)
    get_unit_mass=0
end function

function set_unit_mass(unit_mass)
    implicit none
    double precision :: unit_mass
    integer :: set_unit_mass
    call amuse_set_unit_mass(unit_mass)
    set_unit_mass=0
end function

function get_unit_time(unit_time)
    implicit none
    double precision :: unit_time
    integer :: get_unit_time
    call amuse_get_unit_time(unit_time)
    get_unit_time=0
end function

function set_unit_time(unit_time)
    implicit none
    double precision :: unit_time
    integer :: set_unit_time
    call amuse_set_unit_time(unit_time)
    set_unit_time=0
end function

function get_constant_solarm(solarm)
    implicit none
    double precision :: solarm
    integer :: get_constant_solarm
    call amuse_get_constant_solarm(solarm)
    get_constant_solarm=0
end function

function get_constant_pc(pc)
    implicit none
    double precision :: pc
    integer :: get_constant_pc
    call amuse_get_constant_pc(pc)
    get_constant_pc=0
end function

function get_constant_planckh(planckh)
    implicit none
    double precision :: planckh
    integer :: get_constant_planckh
    call amuse_get_constant_planckh(planckh)
    get_constant_planckh=0
end function

!END MODULE
