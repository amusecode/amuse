      SUBROUTINE get_min_t
      IMPLICIT NONE
#include "general.inc"

      INTEGER          i
      DOUBLE PRECISION min_t_loc,tnext
#ifndef NOMPI
#include <mpif.h>
C---  some short cuts
      INTEGER   ierr
C      INTEGER   MPI_INT,MPI_DOUBLE
C      PARAMETER(MPI_INT=MPI_INTEGER,MPI_DOUBLE=MPI_DOUBLE_PRECISION)
#endif

C---  check for active particles
      IF (n_loc.EQ.0) RETURN    ! this would be bad

C---  get local minimum of next time step
      min_t_loc = t_loc(1)+dt_loc(1)
      DO i=2,n_loc
         tnext = t_loc(i) + dt_loc(i)
         IF (tnext.LT. min_t_loc) min_t_loc = tnext
      ENDDO

#ifndef NOMPI
C---  get global mininum of min_t_loc
      CALL MPI_Allreduce(min_t_loc, min_t, 1, MPI_DOUBLE, MPI_MIN,
     $                                        MPI_COMM_WORLD,ierr)
#endif
#ifdef NOMPI
      min_t = min_t_loc
#endif

#ifdef DEBUG
      IF (myRank.EQ.0) THEN
         WRITE(0,500) min_t
      ENDIF
 500  FORMAT('root: minum time step min_t =',1P,E15.6)
#endif

      RETURN
      END
