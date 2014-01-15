MODULE Mikkola

integer, DIMENSION(:), ALLOCATABLE :: particle_id
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: particle_m
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: particle_y
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: particle_x
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: particle_z
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: particle_vy
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: particle_vx
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: particle_vz
INTEGER, DIMENSION(:), ALLOCATABLE :: particle_child1
INTEGER, DIMENSION(:), ALLOCATABLE :: particle_child2
LOGICAL, DIMENSION(:), ALLOCATABLE :: particle_is_child
integer, DIMENSION(:), ALLOCATABLE :: particle_id_added
DOUBLE PRECISION :: current_time
DOUBLE PRECISION :: begin_time
integer :: maximum_number_of_particles
integer :: number_of_particles_allocated
DOUBLE PRECISION :: lightspeed, tolerance, timestep
integer :: number_of_particles_added
LOGICAL :: evolve_to_exact_time


CONTAINS

FUNCTION new_particle_id(index_of_the_particle)
  IMPLICIT NONE
  INTEGER :: index_of_the_particle, new_particle_id
  INTEGER :: i

  
  new_particle_id = -1
  index_of_the_particle = -1
  DO i = 1, maximum_number_of_particles
    IF (particle_id(i).EQ.-1) THEN
        particle_id(i) = i
        index_of_the_particle = i
        new_particle_id = 0
        number_of_particles_allocated = number_of_particles_allocated + 1
        RETURN
    END IF
  END DO
END FUNCTION
  
FUNCTION commit_parameters()
  IMPLICIT NONE
  INTEGER :: commit_parameters, i
  
  current_time = begin_time
  
  ALLOCATE(particle_id(maximum_number_of_particles))
  ALLOCATE(particle_m(maximum_number_of_particles))
  ALLOCATE(particle_x(maximum_number_of_particles))
  ALLOCATE(particle_y(maximum_number_of_particles))
  ALLOCATE(particle_z(maximum_number_of_particles))
  ALLOCATE(particle_vx(maximum_number_of_particles))
  ALLOCATE(particle_vy(maximum_number_of_particles))
  ALLOCATE(particle_vz(maximum_number_of_particles))
  ALLOCATE(particle_child1(maximum_number_of_particles))
  ALLOCATE(particle_child2(maximum_number_of_particles))
  ALLOCATE(particle_is_child(maximum_number_of_particles))
  
  ALLOCATE(particle_id_added(maximum_number_of_particles))

  DO i = 1, maximum_number_of_particles
    particle_id(i) = -1
  END DO
  
  particle_is_child = .FALSE.
  particle_child1 = 0
  particle_child2 = 0
  particle_m = 0
  
  number_of_particles_allocated = 0
  commit_parameters=0
END FUNCTION

FUNCTION new_particle(index_of_the_particle, m, x, y, z, vx, vy, vz, r)
  IMPLICIT NONE
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: m, r
  DOUBLE PRECISION :: x, y, z
  DOUBLE PRECISION :: vx, vy, vz
  
  INTEGER :: new_particle, i
  new_particle = new_particle_id(index_of_the_particle)
  
  
  IF (new_particle.EQ.-1) THEN
      new_particle = -1
  ELSE
      particle_m(index_of_the_particle) = m
      particle_x(index_of_the_particle) = x
      particle_y(index_of_the_particle) = y
      particle_z(index_of_the_particle) = z
      particle_vx(index_of_the_particle) = vx
      particle_vy(index_of_the_particle) = vy
      particle_vz(index_of_the_particle) = vz
      particle_child1(index_of_the_particle) = 0
      particle_child2(index_of_the_particle) = 0
      particle_is_child(index_of_the_particle) = .FALSE.
      new_particle = 0
  END IF
END FUNCTION

FUNCTION get_mass(index_of_the_particle, mass)
  IMPLICIT NONE
  INTEGER :: get_mass
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: mass
  IF ( (.NOT. ALLOCATED(particle_id)) .OR. (particle_id(index_of_the_particle).EQ.-1) ) THEN
      get_mass = -1
      RETURN
  ENDIF
  mass = particle_m(index_of_the_particle)
  get_mass=0
END FUNCTION

FUNCTION get_velocity(index_of_the_particle, vx, vy, vz)
  IMPLICIT NONE
  INTEGER :: get_velocity
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: vx, vy, vz
  IF ( (.NOT. ALLOCATED(particle_id)) .OR. (particle_id(index_of_the_particle).EQ.-1) ) THEN
      get_velocity = -1
      RETURN
  ENDIF
  vx = particle_vx(index_of_the_particle)
  vy = particle_vy(index_of_the_particle)
  vz = particle_vz(index_of_the_particle)
  get_velocity=0
END FUNCTION

FUNCTION set_velocity(index_of_the_particle, vx, vy, vz)
  IMPLICIT NONE
  INTEGER :: set_velocity
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: vx, vy, vz
  IF ( (.NOT. ALLOCATED(particle_id)) .OR. (particle_id(index_of_the_particle).EQ.-1) ) THEN
      set_velocity = -1
      RETURN
  ENDIF
  particle_vx(index_of_the_particle) = vx
  particle_vy(index_of_the_particle) = vy
  particle_vz(index_of_the_particle) = vz
  set_velocity=0
END FUNCTION

FUNCTION get_position(index_of_the_particle, x, y, z)
  IMPLICIT NONE
  INTEGER :: get_position
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: x, y, z
  IF ( (.NOT. ALLOCATED(particle_id)) .OR. (particle_id(index_of_the_particle).EQ.-1) ) THEN
      get_position = -1
      RETURN
  ENDIF
  x = particle_x(index_of_the_particle)
  y = particle_y(index_of_the_particle)
  z = particle_z(index_of_the_particle)
  get_position=0
END FUNCTION

FUNCTION set_position(index_of_the_particle, x, y, z)
  IMPLICIT NONE
  INTEGER :: set_position
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: x, y, z
  IF ( (.NOT. ALLOCATED(particle_id)) .OR. (particle_id(index_of_the_particle).EQ.-1) ) THEN
      set_position = -1
      RETURN
  ENDIF
  particle_x(index_of_the_particle) = x
  particle_y(index_of_the_particle) = y
  particle_z(index_of_the_particle) = z
  set_position=0
END FUNCTION

FUNCTION get_state(index_of_the_particle, mass, x, y, z, vx, vy,  &
    vz, radius)
  IMPLICIT NONE
  INTEGER :: get_state
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: mass, radius, x, y, z, vx, vy, vz
  IF ( (.NOT. ALLOCATED(particle_id)) .OR. (particle_id(index_of_the_particle).EQ.-1) ) THEN
      get_state = -1
      RETURN
  ENDIF
  radius = 0
  mass = particle_m(index_of_the_particle)
  x = particle_x(index_of_the_particle)
  y = particle_y(index_of_the_particle)
  z = particle_z(index_of_the_particle)
  vx = particle_vx(index_of_the_particle)
  vy = particle_vy(index_of_the_particle)
  vz = particle_vz(index_of_the_particle)
  get_state=0
END FUNCTION

FUNCTION set_mass(index_of_the_particle, mass)
  IMPLICIT NONE
  INTEGER :: set_mass
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: mass
  IF ( (.NOT. ALLOCATED(particle_id)) .OR. (particle_id(index_of_the_particle).EQ.-1) ) THEN
      set_mass = -1
      RETURN
  ENDIF
  particle_m(index_of_the_particle) = mass  
  set_mass=0
END FUNCTION

FUNCTION get_time(time)
  IMPLICIT NONE
  INTEGER :: get_time
  DOUBLE PRECISION :: time
  time = current_time
  get_time=0
END FUNCTION

FUNCTION set_begin_time(time)
  IMPLICIT NONE
  INTEGER :: set_begin_time
  DOUBLE PRECISION :: time
  begin_time = time
  set_begin_time=0
END FUNCTION

FUNCTION get_begin_time(time)
  IMPLICIT NONE
  INTEGER :: get_begin_time
  DOUBLE PRECISION :: time
  time = begin_time
  get_begin_time=0
END FUNCTION

FUNCTION evolve_model(end_time)
  IMPLICIT NONE
  INTEGER :: evolve_model
  DOUBLE PRECISION :: end_time
  DOUBLE PRECISION :: POS(3,maximum_number_of_particles*3)
  INTEGER :: mergers(3, maximum_number_of_particles*3)
  INTEGER :: nmergers
  DOUBLE PRECISION :: VEL(3,maximum_number_of_particles*3)
  DOUBLE PRECISION :: BODY(maximum_number_of_particles*3)
  INTEGER :: INDEX(maximum_number_of_particles*3), REVERSE_INDEX(maximum_number_of_particles*3)
  DOUBLE PRECISION :: IWRR, DELTAT, TEND, soft, cmet(3), tolerance
  DOUBLE PRECISION :: BHspin(3)
  INTEGER :: stepr, i, j, Mikkola_ARWV
  INTEGER :: Np, Nbh, Ixc, error
  INTEGER :: idparent, idchild1, idchild2, new_index
  
  j = 1
  Np = 0
  REVERSE_INDEX = 0
  DO i=1, maximum_number_of_particles
     IF (particle_id(i).NE.-1 .AND. (.NOT. particle_is_child(i))) THEN
         INDEX(j) = j
         REVERSE_INDEX(j) = i
         POS(1,j) = particle_x(i) 
         POS(2,j) = particle_y(i) 
         POS(3,j) = particle_z(i) 
         VEL(1,j) = particle_vx(i) 
         VEL(2,j) = particle_vy(i) 
         VEL(3,j) = particle_vz(i)
         BODY(j) = particle_m(i) 
         j = j + 1
         Np = Np + 1
     ENDIF
  ENDDO
  
  
  Nbh = Np
  IWRR = -0 !?
  DELTAT = MIN(timestep, end_time - current_time) ! Initial timestep, not used according to Mikkola
  if (DELTAT .LE. 0.0) then
     DELTAT = 0.001 ! if timestep or end_time invalid ensure a valid deltat
  end if  
  number_of_particles_added  = 0
  nmergers = 0
  
!  TMAX = 12560 ! Maximum integration time
  stepr = 0 ! Not used, should be maximum number of steps
  soft= 0.e-6 ! Softening parameter
  cmet= [1.e-0, 0.e-0, 0.e-0] !?
  if (evolve_to_exact_time) then
    Ixc=2 ! time output is exacte (2) or not exact (1, faster)
  else
    Ixc=1
  endif
  BHspin=[0.0, 0.0, 0.0] !spin of the first black hole (between 0 and 1) 
  
  evolve_model = Mikkola_ARWV(current_time, BODY, POS,VEL,INDEX, &
&                IWRR,Np,DELTAT,end_time,stepr,soft,cmet,  &
&                lightspeed,Ixc,Nbh,BHspin,tolerance, &
&                mergers, nmergers) 
  
  j = 1
! this print statement is needed for gfortran 4.8.x 
! otherwise nmergers seems to optimized out
  PRINT *, "nmergers:", nmergers 
  DO i=1, nmergers
    idparent = mergers(1,i)
    idchild1 = mergers(2,i)
    idchild2 = mergers(3,i)
    error = new_particle_id(new_index)
    if(error.LT.0) then
        evolve_model = -3
        return
    endif
    particle_m(new_index) = BODY(idparent)
    particle_child1(new_index) =  REVERSE_INDEX(idchild1)
    particle_child2(new_index) = REVERSE_INDEX(idchild2)
    REVERSE_INDEX(Np+i) = new_index
    particle_id_added(i) = new_index
  ENDDO
  
  number_of_particles_added = number_of_particles_added + nmergers
  DO i=1, (nmergers + Np)
      particle_x(REVERSE_INDEX(i)) = POS(1,i)
      particle_y(REVERSE_INDEX(i)) = POS(2,i)
      particle_z(REVERSE_INDEX(i)) = POS(3,i)
      particle_vx(REVERSE_INDEX(i)) = VEL(1,i)
      particle_vy(REVERSE_INDEX(i)) = VEL(2,i)
      particle_vz(REVERSE_INDEX(i)) = VEL(3,i)
  ENDDO
  
  
  DO i=1, nmergers
    idparent = mergers(1,i)
    idchild1 = mergers(2,i)
    idchild2 = mergers(3,i)
    
    particle_is_child(REVERSE_INDEX(idchild1)) = .TRUE.
    particle_is_child(REVERSE_INDEX(idchild2)) = .TRUE.
  ENDDO
  current_time = end_time
END FUNCTION

FUNCTION get_index_of_first_particle(index_of_the_particle)
  IMPLICIT NONE
  INTEGER :: get_index_of_first_particle
  INTEGER :: index_of_the_particle
  get_index_of_first_particle=0
END FUNCTION

FUNCTION get_total_radius(radius)
  IMPLICIT NONE
  INTEGER :: get_total_radius
  DOUBLE PRECISION :: radius
  get_total_radius=0
END FUNCTION

FUNCTION get_potential_at_point(eps, x, y, z, phi)
  IMPLICIT NONE
  INTEGER :: get_potential_at_point
  DOUBLE PRECISION :: eps, x, y, z, phi
  get_potential_at_point=0
END FUNCTION


FUNCTION get_total_mass(mass)
  IMPLICIT NONE
  INTEGER :: get_total_mass
  DOUBLE PRECISION :: mass
  get_total_mass=0
END FUNCTION

FUNCTION set_eps2(epsilon_squared)
  IMPLICIT NONE
  INTEGER :: set_eps2
  DOUBLE PRECISION :: epsilon_squared
  set_eps2=0
END FUNCTION

FUNCTION get_eps2(epsilon_squared)
  IMPLICIT NONE
  INTEGER :: get_eps2
  DOUBLE PRECISION :: epsilon_squared
  get_eps2=0
END FUNCTION

FUNCTION set_lightspeed(inputvalue)
  IMPLICIT NONE
  INTEGER :: set_lightspeed
  DOUBLE PRECISION, intent(in) :: inputvalue
  lightspeed = inputvalue
  set_lightspeed=0
END FUNCTION

FUNCTION get_lightspeed(outputvalue)
  IMPLICIT NONE
  INTEGER :: get_lightspeed
  DOUBLE PRECISION, intent(out) :: outputvalue
  outputvalue = lightspeed
  get_lightspeed=0
END FUNCTION

FUNCTION set_time_step(inputvalue)
  IMPLICIT NONE
  INTEGER :: set_time_step
  DOUBLE PRECISION, intent(in) :: inputvalue
  timestep = inputvalue
  set_time_step=0
END FUNCTION

FUNCTION get_time_step(outputvalue)
  IMPLICIT NONE
  INTEGER :: get_time_step
  DOUBLE PRECISION, intent(out) :: outputvalue
  outputvalue = timestep
  get_time_step=0
END FUNCTION

FUNCTION set_tolerance(inputvalue)
  IMPLICIT NONE
  INTEGER :: set_tolerance
  DOUBLE PRECISION, intent(in) :: inputvalue
  tolerance = inputvalue
  set_tolerance=0
END FUNCTION

FUNCTION get_tolerance(outputvalue)
  IMPLICIT NONE
  INTEGER :: get_tolerance
  DOUBLE PRECISION, intent(out) :: outputvalue
  outputvalue = tolerance
  get_tolerance=0
END FUNCTION

FUNCTION get_index_of_next_particle(index_of_the_particle,  &
    index_of_the_next_particle)
  IMPLICIT NONE
  INTEGER :: get_index_of_next_particle
  INTEGER :: index_of_the_particle, index_of_the_next_particle
  get_index_of_next_particle=0
END FUNCTION

FUNCTION delete_particle(index_of_the_particle)
  IMPLICIT NONE
  INTEGER :: delete_particle
  INTEGER :: index_of_the_particle
  IF ( (.NOT. ALLOCATED(particle_id)) .OR. (particle_id(index_of_the_particle).EQ.-1) ) THEN
      delete_particle = -1
      RETURN
  ENDIF
  particle_id(index_of_the_particle)  = -1
  number_of_particles_allocated = number_of_particles_allocated - 1
  delete_particle=0
END FUNCTION

FUNCTION get_potential(index_of_the_particle, potential)
  IMPLICIT NONE
  INTEGER ::get_potential
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: potential
  get_potential=0
END FUNCTION

FUNCTION synchronize_model()
  IMPLICIT NONE
  INTEGER ::synchronize_model
  synchronize_model=0
END FUNCTION

FUNCTION set_state(index_of_the_particle, mass, x, y, z, vx, vy,  &
    vz, radius)
  IMPLICIT NONE
  INTEGER :: set_state
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: mass, radius, x, y, z, vx, vy, vz
  IF ( (.NOT. ALLOCATED(particle_id)) .OR. (particle_id(index_of_the_particle).EQ.-1) ) THEN
      set_state = -1
      RETURN
  ENDIF
  particle_m(index_of_the_particle) = mass
  particle_x(index_of_the_particle) = x
  particle_y(index_of_the_particle) = y
  particle_z(index_of_the_particle) = z
  particle_vx(index_of_the_particle) = vx
  particle_vy(index_of_the_particle) = vy
  particle_vz(index_of_the_particle) = vz
  set_state=0
END FUNCTION


FUNCTION commit_particles()
  IMPLICIT NONE
  INTEGER :: commit_particles
  commit_particles=0
END FUNCTION

FUNCTION recommit_particles()
  IMPLICIT NONE
  INTEGER :: recommit_particles
  recommit_particles=0
END FUNCTION

FUNCTION get_kinetic_energy(kinetic_energy)
  IMPLICIT NONE
  INTEGER :: get_kinetic_energy
  DOUBLE PRECISION :: kinetic_energy
  REAL*8 get_tkin
  kinetic_energy = get_tkin()
  get_kinetic_energy= 0
END FUNCTION


FUNCTION get_radiated_gravitational_energy(output)
  IMPLICIT NONE
  INTEGER :: get_radiated_gravitational_energy
  DOUBLE PRECISION :: output
  REAL*8 get_energr
  output = get_energr()
  get_radiated_gravitational_energy= 0
END FUNCTION

FUNCTION get_total_energy(output)
  IMPLICIT NONE
  INTEGER :: get_total_energy
  DOUBLE PRECISION :: output
  REAL*8 get_energy
  output = get_energy()
  get_total_energy= 0
END FUNCTION

FUNCTION get_number_of_particles(number_of_particles)
  IMPLICIT NONE
  INTEGER :: get_number_of_particles
  INTEGER :: number_of_particles
  get_number_of_particles=number_of_particles_allocated
END FUNCTION

FUNCTION set_acceleration(index_of_the_particle, ax, ay, az)
  IMPLICIT NONE
  INTEGER :: set_acceleration
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: ax, ay, az
  set_acceleration=0
END FUNCTION

FUNCTION get_indices_of_colliding_particles(index_of_particle1,  &
    index_of_particle2)
  IMPLICIT NONE
  INTEGER :: get_indices_of_colliding_particles
  INTEGER :: index_of_particle1, index_of_particle2
  get_indices_of_colliding_particles=0
END FUNCTION

FUNCTION get_center_of_mass_position(x, y, z)
  IMPLICIT NONE
  INTEGER :: get_center_of_mass_position
  DOUBLE PRECISION :: x, y, z
  get_center_of_mass_position=0
END FUNCTION

FUNCTION get_center_of_mass_velocity(vx, vy, vz)
  IMPLICIT NONE
  INTEGER :: get_center_of_mass_velocity
  DOUBLE PRECISION :: vx, vy, vz
  get_center_of_mass_velocity=0
END FUNCTION

FUNCTION get_radius(index_of_the_particle, radius)
  IMPLICIT NONE
  INTEGER :: get_radius
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: radius
  get_radius=0
END FUNCTION

FUNCTION set_radius(index_of_the_particle, radius)
  IMPLICIT NONE
  INTEGER :: set_radius
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: radius
  set_radius=0
END FUNCTION

FUNCTION cleanup_code()
  IMPLICIT NONE
  INTEGER :: cleanup_code
  
  IF (ALLOCATED(particle_id)) THEN
      DEALLOCATE(particle_id)
      DEALLOCATE(particle_m)
      DEALLOCATE(particle_x)
      DEALLOCATE(particle_y)
      DEALLOCATE(particle_z)
      DEALLOCATE(particle_vx)
      DEALLOCATE(particle_vy)
      DEALLOCATE(particle_vz)
      DEALLOCATE(particle_child1)
      DEALLOCATE(particle_child2)
      DEALLOCATE(particle_is_child)
      
      DEALLOCATE(particle_id_added)
  END IF
  
  number_of_particles_allocated = 0
  cleanup_code=0
END FUNCTION






FUNCTION recommit_parameters()
  IMPLICIT NONE
  INTEGER :: recommit_parameters
  recommit_parameters=0
END FUNCTION

FUNCTION initialize_code()
  IMPLICIT NONE
  INTEGER :: initialize_code
  maximum_number_of_particles = 100
  initialize_code=0
  current_time = 0
  begin_time = 0
  lightspeed = 1
  tolerance = 1.e-13 ! accuracy parameter to which to integrate
  timestep = 1.0
  evolve_to_exact_time = .TRUE.
END FUNCTION

FUNCTION get_potential_energy(potential_energy)
  IMPLICIT NONE
  INTEGER :: get_potential_energy
  DOUBLE PRECISION :: potential_energy
  REAL*8 get_upot
  potential_energy = -get_upot() !mikkola store potential energy as a positive value, in amuse in needs to be negative
  get_potential_energy= 0
END FUNCTION

! eps is ignored, and eps2 = 0
FUNCTION get_gravity_at_point(eps, x, y, z, forcex, forcey, forcez)
  IMPLICIT NONE
  INTEGER :: i, get_gravity_at_point
  DOUBLE PRECISION :: eps, x, y, z, forcex, forcey, forcez
  DOUBLE PRECISION :: eps2, r, r2, r3, rx, ry, rz, F
  forcex = 0
  forcey = 0
  forcez = 0
  eps2 = 0
  DO i=1, maximum_number_of_particles
     IF (particle_id(i).NE.-1) THEN
         rx = particle_x(i) - x
         ry = particle_y(i) - y
         rz = particle_z(i) - z
         r2 = (rx*rx+ry*ry+rz*rz + eps2)
         r = sqrt(r2)
         r3 = r2*r
         F = particle_m(i)/r3
         forcex = forcex + F * rx
         forcey = forcey + F * ry
         forcez = forcez + F * rz
     ENDIF
  ENDDO
  get_gravity_at_point=0
END FUNCTION


FUNCTION get_number_of_particles_added(output)
  IMPLICIT NONE
  INTEGER :: get_number_of_particles_added
  INTEGER :: output
  
  output = number_of_particles_added
  get_number_of_particles_added=0
END FUNCTION


FUNCTION get_id_of_added_particle(index_in_list, output)
  IMPLICIT NONE
  INTEGER :: get_id_of_added_particle
  INTEGER :: index_in_list, output
  if(index_in_list >= number_of_particles_added) then
    output = -1
    get_id_of_added_particle = -1
    return
  end if
  output = particle_id_added(index_in_list + 1)
  
  get_id_of_added_particle=0
END FUNCTION

FUNCTION get_acceleration(index_of_the_particle, ax, ay, az)
  IMPLICIT NONE
  INTEGER :: get_acceleration
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: ax, ay, az
  get_acceleration=0
END FUNCTION



FUNCTION get_children_of_particle(index_of_the_particle, child1, child2)
  IMPLICIT NONE
  INTEGER :: get_children_of_particle
  INTEGER :: index_of_the_particle, child1, child2
  IF ( (.NOT. ALLOCATED(particle_id)) .OR. (particle_id(index_of_the_particle).EQ.-1) ) THEN
      get_children_of_particle = -1
      child1 = 0
      child2 = 0
      RETURN
  ENDIF
  child1 = particle_child1(index_of_the_particle)
  child2 = particle_child2(index_of_the_particle)
  get_children_of_particle=0
END FUNCTION


FUNCTION set_maximum_number_of_particles(inputvalue)
  IMPLICIT NONE
  INTEGER :: set_maximum_number_of_particles
  INTEGER, intent(in) :: inputvalue
  maximum_number_of_particles = inputvalue
  set_maximum_number_of_particles=0
END FUNCTION

FUNCTION get_maximum_number_of_particles(outputvalue)
  IMPLICIT NONE
  INTEGER :: get_maximum_number_of_particles
  INTEGER, intent(out) :: outputvalue
  outputvalue = maximum_number_of_particles
  get_maximum_number_of_particles=0
END FUNCTION


FUNCTION set_evolve_to_exact_time(inputvalue)
  IMPLICIT NONE
  INTEGER :: set_evolve_to_exact_time
  LOGICAL, intent(in) :: inputvalue
  evolve_to_exact_time = inputvalue
  set_evolve_to_exact_time=0
END FUNCTION

FUNCTION get_evolve_to_exact_time(outputvalue)
  IMPLICIT NONE
  INTEGER :: get_evolve_to_exact_time
  LOGICAL, intent(out) :: outputvalue
  outputvalue = evolve_to_exact_time
  get_evolve_to_exact_time=0
END FUNCTION

END MODULE 
