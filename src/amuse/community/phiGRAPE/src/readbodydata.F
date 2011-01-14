      SUBROUTINE readbodydata
      IMPLICIT NONE
#include "general.inc"
#ifndef NOMPI
#include <mpif.h>
      INCLUDE 'mpi.inc'
      INTEGER ierr
#endif

C---  local vars
      INTEGER i,j

C--------------------------------------------------------------------
C---  read the global data for particles to the rootRank
C--------------------------------------------------------------------         
      IF (myRank.EQ.rootRank) THEN 
         OPEN(11, FILE=inp_fname)
      
         READ(11,*) diskstep
         READ(11,*) N
         READ(11,*) time_cur
      
         DO i=1,N
            READ(11,*) ind(i), m(i), (x(j,i),j=1,3), (v(j,i),j=1,3)
            ind(i) = ind(i) + 1 ! shift index by one
         ENDDO

         CLOSE(11)

#ifndef SILENT
         WRITE(STDOUT,*) N,' particles have been read at t =',time_cur
         CALL flush(STDOUT)
#endif
      ENDIF

C---  need to broadcast N and time_cur
#ifndef NOMPI
      CALL MPI_Bcast(N,        1, MPI_INT,    rootRank,
     $                                        MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(time_cur, 1, MPI_DOUBLE, rootRank,
     $                                        MPI_COMM_WORLD,ierr)
#endif
      
      RETURN
      END
