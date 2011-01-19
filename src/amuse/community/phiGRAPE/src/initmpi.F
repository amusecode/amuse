      SUBROUTINE startmpi
#ifndef NOMPI
#include <mpif.h>
      INCLUDE 'mpi.inc'
      INTEGER ierr
        CALL mpi_init(ierr)
#endif
      END 

      SUBROUTINE initmpi
C=========================================================================
C===
C===  task:  initialise MPI and set some common MPI vars, set NOMPI at
C===  -----  compilation to turn off MPI
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
#endif

      myRank = 0                ! also needed for serial runs
      n_proc = 1

#ifndef NOMPI
C---  start MPI
C          --------

C---  get infos on this run
      CALL MPI_Comm_Rank(MPI_COMM_WORLD,myRank, ierr)
C          -------------
      CALL MPI_Comm_Size(MPI_COMM_WORLD,n_proc, ierr)
C          -------------
      CALL MPI_Get_Processor_Name(processor_name, name_proc, ierr)
C          ----------------------


#ifdef DEBUG
      IF (myRank.EQ.rootRank) THEN
         WRITE(processor_name,'(A20,A)') processor_name,' (root)'
         WRITE(STDERR,'("PE",I2.2,": running on ",A50)') myRank,
     $                                                   processor_name
      ELSE
         WRITE(STDERR,'("PE",I2.2,": running on ",A50)') myRank,
     $                                                   processor_name
      ENDIF
      CALL flush(STDERR)
#endif  

      CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
C          -----------
#endif


      RETURN
      END
