module chem_mod
  use krome_main
  use krome_user ! use utility (for krome_idx_* constants and others)
  implicit none

  type particle_type
    integer :: index_of_the_particle
    double precision :: number_density
    double precision :: temperature
    double precision :: ionrate
    double precision :: abundances(krome_nmols)
  end type

  type(particle_type), allocatable :: particles(:)

  double precision :: tcurrent ! [yr]

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
    particles(:)%number_density = 0.d0
    call krome_init()
    ret=0
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
    double precision :: cr, T
    double precision :: n(krome_nmols)
    integer :: ret
    n = particle%abundances * particle%number_density
    T = particle%temperature
    cr = particle%ionrate
    call krome_set_user_crate(cr)
    call krome_set_user_Av(1.d0)
    call krome_set_user_Tdust(1.d1)
    call krome(n, T, dt)
    particle%temperature = T
    particle%abundances = n/particle%number_density
    ret = 0
  end function evolve_1_particle

  function get_particle_density(index_of_the_particle, number_density) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(out) :: number_density
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    number_density = particles(index)%number_density
    ret = 0
  end function get_particle_density

  function set_particle_density(index_of_the_particle, number_density) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(in) :: number_density
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    particles(index)%number_density = number_density
    ret = 0
  end function set_particle_density

  function get_particle_temperature(index_of_the_particle, temperature) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(out) :: temperature
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    temperature = particles(index)%temperature
    ret = 0
  end function get_particle_temperature

  function set_particle_temperature(index_of_the_particle, temperature) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(in) :: temperature
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    particles(index)%temperature = temperature
    ret = 0
  end function set_particle_temperature

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

  function get_particle_abundance(index_of_the_particle, species_index, abundance) result(ret)
    integer, intent(in) :: index_of_the_particle, species_index
    double precision, intent(out) :: abundance
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    if (species_index .LT. 1 .OR. species_index .GT. krome_nmols) then
      ret = -1
      return
    endif
    abundance = particles(index)%abundances(species_index)
    ret = 0
  end function get_particle_abundance

  function set_particle_abundance(index_of_the_particle, species_index, abundance) result(ret)
    integer, intent(in) :: index_of_the_particle, species_index
    double precision, intent(in) :: abundance
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    if (species_index .LT. 1 .OR. species_index .GT. krome_nmols) then
      ret = -1
      return
    endif
    particles(index)%abundances(species_index) = abundance
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

  function get_particle_state(index_of_the_particle, number_density, temperature, ionrate) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(out) :: number_density, temperature, ionrate
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    number_density = particles(index)%number_density
    temperature = particles(index)%temperature
    ionrate = particles(index)%ionrate
    ret = 0
  end function get_particle_state

  function set_particle_state(index_of_the_particle, number_density, temperature, ionrate) result(ret)
    integer, intent(in) :: index_of_the_particle
    double precision, intent(in) :: number_density, temperature, ionrate
    integer :: index, ret
    index = find_particle(index_of_the_particle)
    if (index .LT. 0) then
      ret = index
      return
    endif
    particles(index)%number_density = number_density
    particles(index)%temperature = temperature
    particles(index)%ionrate = ionrate
    ret = 0
  end function set_particle_state

  function add_particle(index_of_the_particle, number_density, temperature, ionrate) result(ret)
    integer, intent(out) :: index_of_the_particle
    double precision, intent(in) :: number_density, temperature, ionrate
    double precision :: x(krome_nmols)
    integer :: i, ret

    particles_searcheable = .FALSE.
    i = nparticle + 1
    if (i .GT. NMAX) then
      ret = -1
      return
    endif

    index_of_the_particle = new_id()
    particles(i)%index_of_the_particle = index_of_the_particle
    particles(i)%number_density = number_density
    particles(i)%temperature = temperature
    particles(i)%ionrate = ionrate
    particles(i)%abundances = 1.d-40

    if (number_density .GT. 0.d0) then
      particles(i)%abundances(KROME_idx_H)  = 1.d0     !H
      particles(i)%abundances(KROME_idx_H2) = 1.d-6    !H2
      particles(i)%abundances(KROME_idx_Hj) = 1.d-4    !H+
      particles(i)%abundances(KROME_idx_He) = 0.0775d0 !He

      x = particles(i)%number_density * particles(i)%abundances

      call krome_scale_Z(x(:), 0.d0) ! scale to solar

      x(krome_idx_Cj) = x(krome_idx_C) !carbon is fully ionized
      x(krome_idx_C)  = 1d-40
      x(krome_idx_e) = krome_get_electrons(x(:))

      particles(i)%abundances = x / particles(i)%number_density
    endif

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
    if (particles(i)%number_density .LT. 0.d0) then
      ret = -4
      return
    endif
    particles(i)%number_density = -1.d0
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
    do while(par(left)%number_density .GT. 0 .AND. left .LT. right)
      left = left + 1
    enddo
    do while(par(right)%number_density .LE. 0 .AND. left .LT. right)
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
  if(par(left)%number_density .GT. 0) left = left + 1
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
  type(particle_type), allocatable, intent (inout) :: buf(:)
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
