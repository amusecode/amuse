      SUBROUTINE sendinpdata
C=========================================================================
C===
C===  task: broadcast necessary data to all nodes and scatter 
C===  ----- nbody data
C===
C=========================================================================
C===                                                   17/08/2006 by SH 
C=========================================================================
      IMPLICIT NONE
#include "general.inc"
#ifndef NOMPI
#include <mpif.h>
      INCLUDE 'mpi.inc'

C---  local vars
      INTEGER ierr

C--------------------------------------------------------------------
C---  distribute the data to all processes
C--------------------------------------------------------------------

C---  Broadcast all useful values to all processors... 
      CALL MPI_Bcast(eps,      1, MPI_DOUBLE, rootRank,
     $                                        MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(t_end,    1, MPI_DOUBLE, rootRank,
     $                                        MPI_COMM_WORLD,ierr)

      CALL MPI_Bcast(eta,      1, MPI_DOUBLE, rootRank,
     $                                        MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(eta_s,    1, MPI_DOUBLE, rootRank,
     $                                        MPI_COMM_WORLD,ierr)

      CALL MPI_Bcast(dt_max,   1, MPI_DOUBLE, rootRank,
     $                                        MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(dt_disk,  1, MPI_DOUBLE, rootRank,
     $                                        MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(dt_contr, 1, MPI_DOUBLE, rootRank,
     $                                        MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(dt_bh,    1, MPI_DOUBLE, rootRank,
     $                                        MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(dt_restart, 1, MPI_DOUBLE, rootRank,
     $                                        MPI_COMM_WORLD,ierr)

      CALL MPI_Bcast(nbh,      1, MPI_INT,    rootRank,
     $                                        MPI_COMM_WORLD,ierr)

      CALL MPI_Bcast(irestart, 1, MPI_INT,    rootRank,
     $                                        MPI_COMM_WORLD,ierr)


#ifdef DEBUG
       WRITE(STDERR,600) myRank
 600   FORMAT (' PE',I2.2,': input data broadcasted ')
#endif

#endif

      RETURN
      END
