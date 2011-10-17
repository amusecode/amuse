      SUBROUTINE sendbodydata
C=========================================================================
C===
C===  task: scatter the nbody data to all nodes, local body data is
C===  ----- stored in arrays ending in '_loc'
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

C---  Scatter the "local" vectors from "global" 
      
      CALL MPI_Scatter(ind, n_loc, MPI_INT,    ind_loc, n_loc,
     $                      MPI_INT,    rootRank, MPI_COMM_WORLD,ierr)
     
      CALL MPI_Scatter(m,   n_loc, MPI_DOUBLE, m_loc,   n_loc,
     $                      MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
      CALL MPI_Scatter(h,   n_loc, MPI_DOUBLE, h_loc,   n_loc,
     $                      MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
      CALL MPI_Scatter(t,   n_loc, MPI_DOUBLE, t_loc,   n_loc,
     $                      MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
      CALL MPI_Scatter(dt,  n_loc, MPI_DOUBLE, dt_loc,  n_loc,
     $                      MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
      CALL MPI_Scatter(x, 3*n_loc, MPI_DOUBLE, x_loc, 3*n_loc,
     $                      MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
      CALL MPI_Scatter(v, 3*n_loc, MPI_DOUBLE, v_loc, 3*n_loc,
     $                      MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)

C---  this data only needs scattering when restarting
      IF (irestart.EQ.1) THEN

         CALL MPI_Scatter(a, 3*n_loc, MPI_DOUBLE, a_loc, 3*n_loc,
     $                      MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
         CALL MPI_Scatter(adot, 3*n_loc, MPI_DOUBLE, adot_loc, 3*n_loc,
     $                      MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
         CALL MPI_Scatter(pot,  n_loc, MPI_DOUBLE, pot_loc,  n_loc,
     $                      MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
      ENDIF

#ifdef DEBUG
       WRITE(STDERR,600) myRank,n_loc,ind_loc(1),ind_loc(n_loc)
 600   FORMAT (' PE',I2.2,': ',I10,' local particles (first and last) '
     $      ,2I10)
#endif

#endif

      RETURN
      END
