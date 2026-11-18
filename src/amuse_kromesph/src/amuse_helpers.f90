module chem_mod
  use krome_main !use krome (mandatory)
  use krome_user !use utility (for krome_idx_* constants and others)
  use krome_commons
  use krome_constants
  use krome_constants_amuse
  implicit none

  type particle_type
    integer :: index_of_the_particle
    double precision :: rho
    double precision :: u
    double precision :: gamma
    double precision :: mu
    double precision :: ionrate
    double precision :: abundances(krome_nmols)
  end type

  type(particle_type), allocatable :: particles(:)

  double precision :: tcurrent  ! time unit = yr

  integer :: nparticle
  integer :: tot_id

  integer, parameter :: NMAX = 1000000

  logical :: particles_searcheable = .FALSE.
  integer, allocatable :: pid(:)

contains

  function chem_initialize()  result(ret)
    integer :: ret
    tcurrent = 0.d0
    nparticle = 0
    tot_id = 0
    if (.not.allocated(particles)) allocate(particles(NMAX))
    particles(:)%rho = 0.d0
    call krome_init()
    ret = 0
  end function chem_initialize

  function chem_end() result(ret)
     integer :: ret
     if (allocated(particles)) deallocate(particles)
     ret = 0
  end function chem_end

  function chem_commit_parameters() result(ret)
    integer :: ret
    ret = 0
  end function chem_commit_parameters

  function chem_commit_particles() result (ret)
    integer :: ret
    particles_searcheable = .FALSE.
    nparticle = clean_particles(particles)
    ret = 0
  end function chem_commit_particles

  function chem_model_time(outtime) result(ret)
    double precision, intent(out) :: outtime
    integer :: ret
    outtime = tcurrent
    ret = 0
  end function chem_model_time

  function evolve_chem_model(tend) result(ret)
    double precision, intent(in) :: tend
    double precision :: dt
    integer :: i, iret, ret
    ret = 0
    dt = tend - tcurrent
    if (dt .LE. 0) return
    do i=1, nparticle
      iret = evolve_1_particle(particles(i), dt)
      ret = min(iret, ret)
    enddo
    tcurrent = tend
  end function evolve_chem_model

  function evolve_1_particle(particle, dt) result(ret)
    type(particle_type), intent(inout) :: particle
    double precision, intent(in) :: dt
    double precision :: rho, u, gamma, mu, cr, T
    double precision :: n(krome_nmols)
    integer :: ret

    n = particle%abundances
    rho = particle%rho
    u = particle%u
    gamma = particle%gamma
    mu = particle%mu
    cr = particle%ionrate

    n(krome_idx_e) = krome_get_electrons(n(:))
    T = (u*(gamma - 1d0) * mu) / boltzmann_erg

    call krome(n, rho, T, dt)

    ! get amu -> g conversion from amuse
    particle%mu = krome_get_mu_x(n(:)) * amu_2_g
    particle%gamma = krome_get_gamma_x(n(:), T)
    particle%u = (boltzmann_erg*T) / ((particle%gamma - 1d0)*particle%mu)
    particle%abundances = n(:)
    ret = 0
  end function evolve_1_particle

  function get_particle_abundance(index_of_the_particle, abundance_index, abundance) result(ret)
    integer, intent(in) :: index_of_the_particle, abundance_index
    double precision, intent(out) :: abundance
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    if (abundance_index .LT. 1 .OR. abundance_index .GT. krome_nmols) then
      ret = -1
      return
    endif
    abundance = particles(index)%abundances(abundance_index)
    ret = 0
  end function get_particle_abundance

  function set_particle_abundance(index_of_the_particle, abundance_index, abundance) result(ret)
    integer, intent(in) :: index_of_the_particle, abundance_index
    double precision, intent(in) :: abundance
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    if (abundance_index .LT. 1 .OR. abundance_index .GT. krome_nmols) then
      ret = -1
      return
    endif
    particles(index)%abundances(abundance_index) = abundance
    ret = 0
  end function set_particle_abundance

  function set_particle_abundances(index_of_the_particle, abundances, N) result(ret)
    integer, intent(in) :: N
    integer, intent(in) :: index_of_the_particle(N)
    double precision, intent(in) :: abundances(N)
    integer :: index, ret
    index = find_particle(index_of_the_particle(1))
    if (index .LT. 0) then
      ret = index
      return
    endif
    if (N .NE. krome_nmols) then
      ret = -2
      return
    endif
    particles(index)%abundances = abundances
    ret = 0
  end function set_particle_abundances

  function set_particle_state(index_of_the_particle, rho, u, gamma, mu, ionrate) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(in) :: rho, u, gamma, mu, ionrate
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    particles(index)%rho = rho
    particles(index)%u = u
    particles(index)%gamma = gamma
    particles(index)%mu = mu
    particles(index)%ionrate = ionrate
    ret = 0
  end function set_particle_state

  function get_particle_state(index_of_the_particle, rho, u, gamma, mu, ionrate) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(out) :: rho, u, gamma, mu, ionrate
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    rho = particles(index)%rho
    u = particles(index)%u
    gamma = particles(index)%gamma
    mu = particles(index)%mu
    ionrate = particles(index)%ionrate
    ret = 0
  end function get_particle_state

  function get_particle_rho(index_of_the_particle, rho) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(out) :: rho
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    rho = particles(index)%rho
    ret = 0
  end function get_particle_rho

  function set_particle_rho(index_of_the_particle, rho) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(in) :: rho
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    particles(index)%rho = rho
    ret = 0
  end function set_particle_rho

  function get_particle_internal_energy(index_of_the_particle, u) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(out) :: u
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    u = particles(index)%u
    ret = 0
  end function get_particle_internal_energy

  function set_particle_internal_energy(index_of_the_particle, u) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(in) :: u
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    particles(index)%u = u
    ret = 0
  end function set_particle_internal_energy

  function get_particle_gamma(index_of_the_particle, gamma) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(out) :: gamma
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    gamma = particles(index)%gamma
    ret = 0
  end function get_particle_gamma

  function set_particle_gamma(index_of_the_particle, gamma) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(in) :: gamma
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    particles(index)%gamma = gamma
    ret = 0
  end function set_particle_gamma

  function get_particle_mu(index_of_the_particle, mu) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(out) :: mu
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    mu = particles(index)%mu
    ret = 0
  end function get_particle_mu

  function set_particle_mu(index_of_the_particle, mu) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(in) :: mu
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    particles(index)%mu = mu
    ret = 0
  end function set_particle_mu

  function get_particle_ionrate(index_of_the_particle, ionrate) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(out) :: ionrate
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    ionrate = particles(index)%ionrate
    ret = 0
  end function get_particle_ionrate

  function set_particle_ionrate(index_of_the_particle, ionrate) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(in) :: ionrate
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    particles(index)%ionrate = ionrate
    ret = 0
  end function set_particle_ionrate

  function add_particle(index_of_the_particle, rho, u, gamma, mu, ionrate) result(ret)
  integer, intent(out) :: index_of_the_particle
  double precision, intent(in) :: rho, u, gamma, mu, ionrate
  ! double precision :: x(krome_nmols)
  integer :: i, ret

  particles_searcheable = .FALSE.
  index_of_the_particle = new_id()
  i = nparticle + 1

  if (i .GT. NMAX) then
    ret = -1
    return
  endif
  particles(i)%index_of_the_particle = index_of_the_particle
  particles(i)%rho = rho
  particles(i)%u = u
  particles(i)%gamma = gamma
  particles(i)%mu = mu
  particles(i)%ionrate = ionrate
  particles(i)%abundances = 1.d-40

  ! if (density.GT.0) then
  !   particles(i)%abundances(KROME_idx_H)  = 1. !H
  !   particles(i)%abundances(KROME_idx_H2) = 1.d-6    !H2
  !   particles(i)%abundances(KROME_idx_Hj) = 1.d-4    !H+
  !   particles(i)%abundances(KROME_idx_He) = 0.0775d0 !He

  !   x=particles(i)%density*particles(i)%abundances

  !   call krome_scale_Z(x(:), 0.d0) ! scale to solar

  !   x(krome_idx_Cj) = x(krome_idx_C) !carbon is fully ionized
  !   x(krome_idx_C)  = 1d-40
  !   x(krome_idx_e) = krome_get_electrons(x(:))

  !   particles(i)%abundances=x/particles(i)%density
  ! endif

  nparticle = nparticle + 1
  ret = 0
  end function add_particle

  function remove_particle(index_of_the_particle) result(ret)
  integer, intent(in) :: index_of_the_particle
  integer :: i, ret
  i = find_particle(index_of_the_particle)
  if (i .LE. 0) then
    ret = i
    return
  endif
  if (particles(i)%rho .LT. 0) then
    ret = -4
    return
  endif
  particles(i)%rho = -1.
  ret = 0
  end function remove_particle

function clean_particles(par) result(np)
  type(particle_type), allocatable, intent(inout) :: par(:)
  type(particle_type) :: tmp
  integer :: left, right, np
  left = 1
  if (.NOT. allocated(par)) then
    np = 0
    return
  endif
  right = size(par)
  if (right .EQ. 0) then
    np = 0
    return
  endif
  do while(.TRUE.)
    do while(par(left)%rho .GT. 0 .AND. left .LT. right)
      left = left + 1
    enddo
    do while(par(right)%rho .LE. 0 .AND. left .LT. right)
      right = right - 1
    enddo
    if (left .LT. right) then
      tmp = par(left)
      par(left) = par(right)
      par(right) = tmp
    else
      exit
    endif
  enddo
  if (par(left)%rho .GT. 0) left = left + 1
  np = left - 1
end function clean_particles

function find_particle(index_of_the_particle) result(index)
  use hashMod, only : hash_type, initHash, find
  integer, intent(in) :: index_of_the_particle
  integer :: index
  type(hash_type), save :: hash
  integer, save :: nbod = 0

  if (.NOT. particles_searcheable) then
    nbod = nparticle
    if (allocated(pid)) deallocate(pid)
    allocate(pid(nbod))
    pid(1:nbod) = particles(1:nbod)%index_of_the_particle
    call initHash(nbod/2+1, nbod, pid, hash)
    particles_searcheable = .TRUE.
  endif
  index = find(index_of_the_particle, pid, hash)
  if (index .LE. 0) then
    index = -1
    return
  endif
  if (index .GT. nbod) then
    index = -2
    return
  endif
  if (pid(index) .NE. index_of_the_particle) then
    index = -3
    return
  endif
end function find_particle

subroutine extend_particles(buf, n)
  type(particle_type), allocatable, intent(inout) :: buf(:)
  type(particle_type), allocatable :: tmpbuf(:)
  integer, intent(in) :: n
  integer :: m

  m = 0
  if (allocated(buf)) then
    m = min(n, size(buf))
    allocate(tmpbuf(m))
    tmpbuf(1:m) = buf(1:m)
    deallocate(buf)
  endif

  allocate(buf(n))

  if (m .GT. 0 .and. allocated(tmpbuf)) then
    buf(1:m) = tmpbuf(1:m)
    deallocate(tmpbuf)
  endif

end subroutine extend_particles

function new_id()
  integer new_id
  tot_id = tot_id + 1
  new_id = tot_id
end function new_id

end module chem_mod
