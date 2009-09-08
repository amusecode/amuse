      SUBROUTINE readrestartfile
C=========================================================================
C===
C===  task:   read all run data from restart file
C===  -----
C=========================================================================
      IMPLICIT NONE
#include "general.inc"
#ifndef NOMPI
#include <mpif.h>
      INCLUDE 'mpi.inc'
      INTEGER ierr,idummy,ipe
      DOUBLE PRECISION xdummy
      INTEGER status(MPI_STATUS_SIZE)
#endif

C---  local vars
      INTEGER i,j,itag
      
      itag = 1

      IF (myRank.EQ.rootRank) THEN
         OPEN(urestf,FILE=restfile,FORM='UNFORMATTED')

c---  some header info
         READ(urestf) time_cur,N,diskstep

C---  main data
         READ(urestf) (ind(i),i=1,N)
         READ(urestf) (m(i),i=1,N)
         READ(urestf) ((x(j,i),j=1,3),i=1,N)
         READ(urestf) ((v(j,i),j=1,3),i=1,N)
         READ(urestf) ((a(j,i),j=1,3),i=1,N)
         READ(urestf) ((adot(j,i),j=1,3),i=1,N)
         READ(urestf) (pot(i),i=1,N)
         READ(urestf) (t(i),i=1,N)
         READ(urestf) (dt(i),i=1,N)
         READ(urestf) Timesteps,n_act_sum,n_act_sum2,cpu_total
         READ(urestf) t_disk,t_contr,t_bh,t_timing,t_restart
      ENDIF

C---  need to broadcast N here so all nodes now how much data they'll recv.
#ifndef NOMPI
      CALL MPI_Bcast(N,        1, MPI_INT,    rootRank,
     $                                        MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(time_cur, 1, MPI_DOUBLE, rootRank,
     $                                        MPI_COMM_WORLD,ierr)
#endif


C---  now need to read *_act_tmp_loc arrays and send to each PE so that
C---  local guessing for partial forces is correct after restart

C---  begin by reading root data; will be read later because
C---  arrays will be used as send buffer for now
C---  (NOTE: in serial runs no re-reading is needed which is why
C---         here the actual data is read instead of skipping over it)
      IF (myRank.EQ.rootRank) THEN
         READ(urestf) ((a_act_tmp_loc(j,i),j=1,3),i=1,N)
         READ(urestf) ((adot_act_tmp_loc(j,i),j=1,3),i=1,N)
         READ(urestf) (pot_act_tmp_loc(i),i=1,N)
      ENDIF

#ifndef NOMPI
      DO ipe=1,n_proc-1

C---  read data for PE ipe 
         IF (myRank.EQ.rootRank) THEN
            READ(urestf) ((a_act_tmp_loc(j,i),j=1,3),i=1,N)
            READ(urestf) ((adot_act_tmp_loc(j,i),j=1,3),i=1,N)
            READ(urestf) (pot_act_tmp_loc(i),i=1,N)
         ENDIF

C---  now send data from root to ipe 
         IF (myRank.EQ.rootRank) THEN 
            CALL MPI_SEND(a_act_tmp_loc, 3*N, MPI_DOUBLE, ipe,
     $                                  itag, MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(adot_act_tmp_loc, 3*N, MPI_DOUBLE,ipe,
     $                                  itag, MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(pot_act_tmp_loc, 3*N, MPI_DOUBLE, ipe,
     $                                  itag, MPI_COMM_WORLD,ierr)
         ELSE IF (myRank.EQ.ipe) THEN
            CALL MPI_RECV(a_act_tmp_loc, 3*N, MPI_DOUBLE, rootRank,
     $                           itag, MPI_COMM_WORLD,status,ierr)
            CALL MPI_RECV(adot_act_tmp_loc, 3*N, MPI_DOUBLE, rootRank,
     $                           itag, MPI_COMM_WORLD,status,ierr)
            CALL MPI_RECV(pot_act_tmp_loc, 3*N, MPI_DOUBLE, rootRank,
     $                           itag, MPI_COMM_WORLD,status,ierr)
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

C---  some data needs broadcasting
      CALL MPI_Bcast(Timesteps, 1, MPI_INT,    rootRank,
     $                                        MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(t_disk,    1, MPI_DOUBLE,    rootRank,
     $                                        MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(t_contr,   1, MPI_DOUBLE,    rootRank,
     $                                        MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(t_bh,      1, MPI_DOUBLE,    rootRank,
     $                                        MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(t_timing,  1, MPI_DOUBLE,    rootRank,
     $                                        MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(t_restart, 1, MPI_DOUBLE,    rootRank,
     $                                        MPI_COMM_WORLD,ierr)


#endif
C---  closing restart-file
      IF (myRank.EQ.rootRank) CLOSE(urestf)


      RETURN
      END
