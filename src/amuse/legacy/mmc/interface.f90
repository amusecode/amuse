MODULE MMC

CONTAINS

FUNCTION nonstandard_init()
  INTEGER :: nonstandard_init
  INTEGER :: init_sequence
  INTEGER :: res
  ! read initial parameters
  PRINT*,'calling input'
  res = init_sequence()
  nonstandard_init = res
  PRINT*,'init done'
END FUNCTION

FUNCTION set_mmc_data_directory(data_directory) 
  INTEGER :: set_mmc_data_directory
  CHARACTER(len=200) :: data_directory
  CALL amuse_set_mmc_data_directory(data_directory)
  set_mmc_data_directory = 0
END FUNCTION

FUNCTION get_time(time)
  IMPLICIT NONE
  INTEGER :: get_time
  INTEGER :: res
  INTEGER :: parameter_test
  DOUBLE PRECISION :: time

  res = parameter_test(time)
  get_time = 0
END FUNCTION

FUNCTION get_kinetic_energy(Ek)
  IMPLICIT NONE
  INTEGER :: res
  INTEGER :: total_kinetic_energy
  INTEGER :: get_kinetic_energy
  DOUBLE PRECISION :: Ek

  res = total_kinetic_energy(Ek)
  get_kinetic_energy = 0
END FUNCTION


FUNCTION new_particle( index_of_the_particle, x, y, z)
  IMPLICIT NONE
  INTEGER :: new_particle
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: x, y, z
  new_particle=0
END FUNCTION

FUNCTION delete_particle( index_of_the_particle)
  IMPLICIT NONE
  INTEGER :: delete_particle
  INTEGER :: index_of_the_particle
  delete_particle=0
END FUNCTION

FUNCTION set_state( index_of_the_particle, x, y, z)
  IMPLICIT NONE
  INTEGER :: set_state
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: x, y, z
  set_state=0
END FUNCTION

FUNCTION get_state( index_of_the_particle, x, y, z)
  IMPLICIT NONE
  INTEGER :: get_state
  INTEGER :: index_of_the_particle
  DOUBLE PRECISION :: x, y, z
  get_state=0
END FUNCTION

FUNCTION get_number_of_particles( VALUE)
  IMPLICIT NONE
  INTEGER :: get_number_of_particles
  INTEGER :: VALUE
  get_number_of_particles=0
END FUNCTION

FUNCTION internal__redirect_outputs( stdoutfile, stderrfile)
  IMPLICIT NONE
  INTEGER :: internal__redirect_outputs
  CHARACTER(LEN=*) :: stdoutfile, stderrfile
  internal__redirect_outputs=0
END FUNCTION

FUNCTION run( )
  IMPLICIT NONE
  INTEGER :: run
  run=0
END FUNCTION

FUNCTION commit_parameters( )
  IMPLICIT NONE
  INTEGER commit_parameters
  commit_parameters=0
END FUNCTION

END MODULE
