      FUNCTION cputime()
C=========================================================================
C===
C===  task: get some measure of time using either
C===  -----          a) MPI_Wtime() in case of parallel runs
C===                 b) etime()     in case of serial runs (NOMPI)
C===
C=========================================================================
C===                                                   17/08/2006 by SH 
C=========================================================================
      IMPLICIT NONE

      DOUBLE PRECISION cputime

#ifdef NOMPI
      REAL thelp(2),cpu

      CALL etime(thelp,cpu)
      cputime = DBLE(cpu)
#endif

#ifndef NOMPI
#include <mpif.h>
      cputime = MPI_Wtime()
#endif

      RETURN
      END
