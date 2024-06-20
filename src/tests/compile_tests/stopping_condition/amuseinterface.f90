INCLUDE 'stoppingconditions.f90'

MODULE AmuseInterface
CONTAINS
      FUNCTION initialize_code()
          use StoppingConditions
          IMPLICIT NONE
          INTEGER :: initialize_code
          INTEGER :: return
          initialize_code = 0
          return = set_support_for_condition(COLLISION_DETECTION)
          return = set_support_for_condition(PAIR_DETECTION)
      END FUNCTION

      FUNCTION fire_condition(condition_to_set, particle_index1, particle_index2, rank)
          use StoppingConditions
          IMPLICIT NONE
          include "mpif.h"
          INTEGER :: fire_condition
          INTEGER :: my_rank
          INTEGER :: error, stopping_index
          INTEGER, intent(in) :: condition_to_set, particle_index1, particle_index2, rank
          fire_condition = 0
          call mpi_comm_rank(MPI_COMM_WORLD, my_rank, error)

          if (rank.GE.0 .AND. rank.NE.my_rank) then
            return
          end if
          stopping_index = next_index_for_stopping_condition()
          error = set_stopping_condition_info(stopping_index, condition_to_set)
          if (particle_index1 .GT. 0) then
            error = set_stopping_condition_particle_index(stopping_index, 0, particle_index1)
          end if
          if (particle_index2 .GT. 0) then
            error = set_stopping_condition_particle_index(stopping_index, 1, particle_index2)
          end if
      END FUNCTION
END MODULE

