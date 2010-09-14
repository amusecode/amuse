MODULE NN

integer, DIMENSION(:), ALLOCATABLE :: identifiers
integer, DIMENSION(:,:), ALLOCATABLE :: nearest_neighbor
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: particle_y
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: particle_x
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: particle_z
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: distances
integer :: maximum_number_of_particles
integer :: number_of_particles_allocated

CONTAINS

FUNCTION get_nearest_neighbor(index_of_the_particle, index_of_the_neighbor,  &
    distance)
  IMPLICIT NONE
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: index_of_the_neighbor, distance
  
  INTEGER :: get_nearest_neighbor
  
  distance = distances(index_of_the_particle)
  index_of_the_neighbor = nearest_neighbor(index_of_the_particle, 0)
  
  get_nearest_neighbor=0
END FUNCTION

FUNCTION new_particle(index_of_the_particle, x, y, z)
  IMPLICIT NONE
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: x, y, z
  
  INTEGER :: new_particle, i
  
  index_of_the_particle = -1
  DO i = 1, maximum_number_of_particles
    IF (identifiers(i).EQ.-1) THEN
        identifiers(i) = i
        index_of_the_particle = i
        number_of_particles_allocated = number_of_particles_allocated + 1
        EXIT
    END IF
  END DO
  
  IF (index_of_the_particle.EQ.-1) THEN
      new_particle = -1
  ELSE
      particle_x(index_of_the_particle) = x
      particle_y(index_of_the_particle) = y
      particle_z(index_of_the_particle) = z
      new_particle = 0
  END IF
END FUNCTION

FUNCTION get_close_neighbors(index_of_the_particle,  &
    index_of_first_neighbor, index_of_second_neighbor,  &
    index_of_third_neighbor)
  IMPLICIT NONE
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: index_of_first_neighbor, index_of_second_neighbor
  DOUBLE PRECISION :: index_of_third_neighbor
  
  INTEGER :: get_close_neighbors
  
  index_of_first_neighbor = nearest_neighbor(index_of_the_particle, 1)
  index_of_second_neighbor = nearest_neighbor(index_of_the_particle, 2)
  index_of_third_neighbor = nearest_neighbor(index_of_the_particle, 3)
  
  get_close_neighbors=0
END FUNCTION

FUNCTION delete_particle(index_of_the_particle)
  IMPLICIT NONE
  INTEGER :: index_of_the_particle
  
  INTEGER :: delete_particle
  
  identifiers(index_of_the_particle) = -1
  number_of_particles_allocated = number_of_particles_allocated - 1
  
  delete_particle=0
END FUNCTION

FUNCTION set_state(index_of_the_particle, x, y, z)
  IMPLICIT NONE
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: x, y, z
  
  INTEGER :: set_state
  
  particle_x(index_of_the_particle) = x
  particle_y(index_of_the_particle) = y
  particle_z(index_of_the_particle) = z
  
  set_state=0
END FUNCTION

FUNCTION get_state(index_of_the_particle, x, y, z)
  IMPLICIT NONE
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: x, y, z
  
  INTEGER :: get_state
 
  x = particle_x(index_of_the_particle)
  y = particle_y(index_of_the_particle)
  z = particle_z(index_of_the_particle)
  
  get_state=0
END FUNCTION

FUNCTION get_number_of_particles(value)
  IMPLICIT NONE
  INTEGER :: value
  
  INTEGER :: get_number_of_particles
  value = number_of_particles_allocated
  get_number_of_particles=0
END FUNCTION

FUNCTION run()
  IMPLICIT NONE
  INTEGER :: run
  INTEGER, DIMENSION(:), ALLOCATABLE :: index_to_table
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: x
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: y
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: z
  INTEGER, DIMENSION(:), ALLOCATABLE :: nn1
  INTEGER, DIMENSION(:), ALLOCATABLE :: nn2
  INTEGER, DIMENSION(:), ALLOCATABLE :: nn3
  INTEGER :: i, j, find_nearest_neighbors
  
  ALLOCATE(x(number_of_particles_allocated))
  ALLOCATE(y(number_of_particles_allocated))
  ALLOCATE(z(number_of_particles_allocated))
  ALLOCATE(nn1(number_of_particles_allocated))
  ALLOCATE(nn2(number_of_particles_allocated))
  ALLOCATE(nn3(number_of_particles_allocated))
  ALLOCATE(index_to_table(number_of_particles_allocated))
  
  j = 0
  DO i = 1, maximum_number_of_particles
    IF (identifiers(i).NE.-1) THEN
        j = j + 1
        x(j) = particle_x(i)
        y(j) = particle_y(i)
        z(j) = particle_z(i)
        index_to_table(j) = i
    END IF
  END DO
  
  run = find_nearest_neighbors(number_of_particles_allocated, x, y, z, nn1, nn2, nn3)
  
  DO j = 1, number_of_particles_allocated
    i = index_to_table(j)
    nearest_neighbor(i, 1) = index_to_table(nn1(j))
    nearest_neighbor(i, 2) = index_to_table(nn2(j))
    nearest_neighbor(i, 3) = index_to_table(nn3(j))
  END DO
  
  DEALLOCATE(x)
  DEALLOCATE(y)
  DEALLOCATE(z)
  DEALLOCATE(nn1)
  DEALLOCATE(nn2)
  DEALLOCATE(nn3)
  DEALLOCATE(index_to_table)
END FUNCTION

FUNCTION set_maximum_number_of_particles(value)
  IMPLICIT NONE
  DOUBLE PRECISION :: value
  
  INTEGER :: set_maximum_number_of_particles
  maximum_number_of_particles = value
  set_maximum_number_of_particles=0
END FUNCTION

FUNCTION commit_parameters()
  IMPLICIT NONE
  
  INTEGER :: commit_parameters, i
  
  ALLOCATE(particle_x(maximum_number_of_particles))
  ALLOCATE(particle_y(maximum_number_of_particles))
  ALLOCATE(particle_z(maximum_number_of_particles))
  ALLOCATE(distances(maximum_number_of_particles))
  ALLOCATE(nearest_neighbor(maximum_number_of_particles, 3))
  ALLOCATE(identifiers(maximum_number_of_particles))
  
  DO i = 1, maximum_number_of_particles
    identifiers(i) = -1
  END DO
  
  number_of_particles_allocated = 0
  commit_parameters=0
END FUNCTION

FUNCTION get_maximum_number_of_particles(value)
  IMPLICIT NONE
  DOUBLE PRECISION :: value
  
  INTEGER :: get_maximum_number_of_particles
  value = maximum_number_of_particles
  get_maximum_number_of_particles=0
END FUNCTION

END MODULE
