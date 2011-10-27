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
      
      CALL MPI_Scatter(ind, n_loc, MPI_INTEGER,  ind_loc, n_loc,
     $                  MPI_INTEGER, rootRank, MPI_COMM_WORLD,ierr)
     
      CALL MPI_Scatter(m,   n_loc, MPI_DOUBLE_PRECISION, m_loc,   n_loc,
     $                      MPI_DOUBLE_PRECISION, rootRank, 
     $                      MPI_COMM_WORLD,ierr)
      CALL MPI_Scatter(h,   n_loc, MPI_DOUBLE_PRECISION, h_loc,   n_loc,
     $                      MPI_DOUBLE_PRECISION, rootRank, 
     $                      MPI_COMM_WORLD,ierr)
      CALL MPI_Scatter(t,   n_loc, MPI_DOUBLE_PRECISION, t_loc,   n_loc,
     $                      MPI_DOUBLE_PRECISION, rootRank, 
     $                      MPI_COMM_WORLD,ierr)
      CALL MPI_Scatter(dt,  n_loc, MPI_DOUBLE_PRECISION, dt_loc,  n_loc,
     $                      MPI_DOUBLE_PRECISION, rootRank, 
     $                      MPI_COMM_WORLD,ierr)
      CALL MPI_Scatter(x, 3*n_loc, MPI_DOUBLE_PRECISION, x_loc, 3*n_loc,
     $                      MPI_DOUBLE_PRECISION, rootRank, 
     $                      MPI_COMM_WORLD,ierr)
      CALL MPI_Scatter(v, 3*n_loc, MPI_DOUBLE_PRECISION, v_loc, 3*n_loc,
     $                      MPI_DOUBLE_PRECISION, rootRank, 
     $                      MPI_COMM_WORLD,ierr)

C---  this data only needs scattering when restarting
      IF (irestart.EQ.1) THEN

         CALL MPI_Scatter(a, 3*n_loc, MPI_DOUBLE_PRECISION, a_loc, 
     $                     3*n_loc, MPI_DOUBLE_PRECISION, rootRank,
     $                      MPI_COMM_WORLD,ierr)
         CALL MPI_Scatter(adot, 3*n_loc, MPI_DOUBLE_PRECISION, adot_loc,
     $                       3*n_loc, MPI_DOUBLE_PRECISION, rootRank, 
     $                      MPI_COMM_WORLD,ierr)
         CALL MPI_Scatter(pot,  n_loc, MPI_DOUBLE_PRECISION, pot_loc, 
     $                       n_loc, MPI_DOUBLE_PRECISION, rootRank,
     $                      MPI_COMM_WORLD,ierr)
      ENDIF

#ifdef DEBUG
       WRITE(STDERR,600) myRank,n_loc,ind_loc(1),ind_loc(n_loc)
 600   FORMAT (' PE',I2.2,': ',I10,' local particles (first and last) '
     $      ,2I10)
#endif

#endif

      RETURN
      END
