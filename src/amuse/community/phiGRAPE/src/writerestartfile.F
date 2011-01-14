      SUBROUTINE writerestartfile
C=========================================================================
C===
C===  task:   write all run data to file for later restart
C===  -----
C=========================================================================
      IMPLICIT NONE
#include "general.inc"
#ifndef NOMPI
#include <mpif.h>
      INCLUDE 'mpi.inc'
      INTEGER          ierr,ipe,idummy
      DOUBLE PRECISION xdummy
      INTEGER status(MPI_STATUS_SIZE)
#endif

C---  local vars
      INTEGER i,j,itag

      
      itag = 1


      IF (myRank.EQ.rootRank) THEN
         OPEN(urestf,FILE=restfile,FORM='UNFORMATTED')

c---  some header info
         WRITE(urestf) time_cur,N,diskstep

C---  main data
         WRITE(urestf) (ind(i),i=1,N)
         WRITE(urestf) (m(i),i=1,N)
         WRITE(urestf) ((x(j,i),j=1,3),i=1,N)
         WRITE(urestf) ((v(j,i),j=1,3),i=1,N)
         WRITE(urestf) ((a(j,i),j=1,3),i=1,N)
         WRITE(urestf) ((adot(j,i),j=1,3),i=1,N)
         WRITE(urestf) (pot(i),i=1,N)
         WRITE(urestf) (t(i),i=1,N)
         WRITE(urestf) (dt(i),i=1,N)
         WRITE(urestf) Timesteps,n_act_sum,n_act_sum2,cpu_total
         WRITE(urestf) t_disk,t_contr,t_bh,t_timing,t_restart
      ENDIF

C---  now need to store *_act_tmp_loc arrays for each PE so that
C---  local guessing for partial forces is correct after restart

C---  begin writing data from root 
      IF (myRank.EQ.rootRank) THEN
         WRITE(urestf) ((a_act_tmp_loc(j,i),j=1,3),i=1,N)
         WRITE(urestf) ((adot_act_tmp_loc(j,i),j=1,3),i=1,N)
         WRITE(urestf) (pot_act_tmp_loc(i),i=1,N)
      ENDIF

#ifndef NOMPI
      DO ipe=1,n_proc-1

C---  now send data from ipe to root 
C---  arrays have to be reset on root later by reading back from file
         IF (myRank.EQ.ipe) THEN 
            CALL MPI_SEND(a_act_tmp_loc, 3*N, MPI_DOUBLE, rootRank,
     $                                  itag, MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(adot_act_tmp_loc, 3*N, MPI_DOUBLE, rootRank,
     $                                  itag, MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(pot_act_tmp_loc, 3*N, MPI_DOUBLE, rootRank,
     $                                  itag, MPI_COMM_WORLD,ierr)
         ELSE IF (myRank.EQ.rootRank) THEN
            CALL MPI_RECV(a_act_tmp_loc, 3*N, MPI_DOUBLE, ipe,
     $                           itag, MPI_COMM_WORLD,status,ierr)
            CALL MPI_RECV(adot_act_tmp_loc, 3*N, MPI_DOUBLE, ipe,
     $                           itag, MPI_COMM_WORLD,status,ierr)
            CALL MPI_RECV(pot_act_tmp_loc, 3*N, MPI_DOUBLE, ipe,
     $                           itag, MPI_COMM_WORLD,status,ierr)
         ENDIF

C---  write data from ipe 
         IF (myRank.EQ.rootRank) THEN
            WRITE(urestf) ((a_act_tmp_loc(j,i),j=1,3),i=1,N)
            WRITE(urestf) ((adot_act_tmp_loc(j,i),j=1,3),i=1,N)
            WRITE(urestf) (pot_act_tmp_loc(i),i=1,N)
         ENDIF

      ENDDO

C---  now get correct data back on root
      IF (myRank.EQ.rootRank) THEN
         CLOSE(urestf)
         OPEN(urestf,FILE=restfile,FORM='UNFORMATTED') ! reopen to read data

C---  skipping the part not needed
         READ(urestf) xdummy
         READ(urestf) idummy
         DO i=1,8
            READ(urestf) xdummy
         ENDDO
         READ(urestf) idummy
         READ(urestf) xdummy


C---  reading the needed data
         READ(urestf) ((a_act_tmp_loc(j,i),j=1,3),i=1,N)
         READ(urestf) ((adot_act_tmp_loc(j,i),j=1,3),i=1,N)
         READ(urestf) (pot_act_tmp_loc(i),i=1,N)

      ENDIF

#endif
C---  closing restart-file
      IF (myRank.EQ.rootRank) CLOSE(urestf)


      RETURN
      END
