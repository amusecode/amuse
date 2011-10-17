      SUBROUTINE allgather_ap
C=========================================================================
C===
C===  task:  collect predicted positions and velocities of active  
C===  -----  particles from all nodes on root, then redistribute so that 
C===         all nodes know all active particles
C===
C=========================================================================
      IMPLICIT NONE
#include "general.inc"


#ifndef NOMPI
#include <mpif.h>
      INCLUDE 'mpi.inc'
      INTEGER i,ierr
#endif

C---  first get information about how many active particles are on each node         
      disp_act_v(1)  = 0           ! these two values are needed also in 
      disp3_act_v(1) = 0           ! the serial run (NOMPI)

#ifndef NOMPI
C---  collect all n_act_loc in one array n_act_v
      CALL MPI_Allgather(n_act_loc, 1,MPI_INT, n_act_v, 1, MPI_INT,
     $                   MPI_COMM_WORLD,ierr)
      
      
      DO i=1,n_proc
         n3_act_v(i) = 6*n_act_v(i)      ! 16 is the length of one particle
      ENDDO                               ! data set in the body_act array
    
C---  get a displacement vector so that all vectors can be connected seamless
C---  also get total number of active particles
      n_act          = n_act_v(1)
      DO i=2,n_proc

         disp_act_v(i)  = disp_act_v(i-1) + n_act_v(i-1)
         disp3_act_v(i) = disp3_act_v(i-1) + n3_act_v(i-1)
         n_act          = n_act + n_act_v(i)
         
      ENDDO
      
C---  allgather global indices of all active particles
C---  NOTE: Gatherv + Bcast is equivalent to Allgatherv but no extra 
C---        memory is needed 
      if(myRank .eq. rootRank) then
      
        CALL MPI_Gatherv(MPI_IN_PLACE, n_act_loc,MPI_INT, ind_act,
     $                 n_act_v,
     $                 disp_act_v, MPI_INT,rootRank,
     $                                           MPI_COMM_WORLD,ierr)
      else       
        CALL MPI_Gatherv(ind_act, n_act_loc, MPI_INT, ind_act, n_act_v,
     $                 disp_act_v, MPI_INT,rootRank,
     $                                           MPI_COMM_WORLD,ierr)
      endif
      
      
    
      CALL MPI_Bcast(ind_act, n_act, MPI_INT, rootRank,
     $                                           MPI_COMM_WORLD,ierr)

C---  allgather the predicted positions and velocities of active particles
      CALL MPI_Allgatherv(body_pred_loc, 6*n_act_loc,       MPI_DOUBLE, 
     $                    body_pred, n3_act_v, disp3_act_v, MPI_DOUBLE,
     $                                              MPI_COMM_WORLD,ierr)

      
    
     
      if(myRank .eq. rootRank) then
      CALL MPI_Gatherv(MPI_IN_PLACE, n_act_loc, MPI_DOUBLE, h_act, 
     $                 n_act_v, disp_act_v, MPI_DOUBLE,rootRank,
     $                                          MPI_COMM_WORLD,ierr)
      else
      CALL MPI_Gatherv(h_act, n_act_loc, MPI_DOUBLE, h_act, 
     $                 n_act_v, disp_act_v, MPI_DOUBLE,rootRank,
     $                           MPI_COMM_WORLD,ierr)
      endif
      
      CALL MPI_Bcast(h_act, n_act, MPI_DOUBLE, rootRank,
     $                                          MPI_COMM_WORLD,ierr)

     print *, "m_act: ", m_act(1:2)
      if(myRank .eq. rootRank) then
      CALL MPI_Gatherv(MPI_IN_PLACE, n_act_loc, MPI_DOUBLE, m_act, 
     $                 n_act_v, disp_act_v, MPI_DOUBLE,rootRank,
     $                                          MPI_COMM_WORLD,ierr)
      else
      CALL MPI_Gatherv(m_act, n_act_loc, MPI_DOUBLE, m_act, 
     $                 n_act_v, disp_act_v, MPI_DOUBLE,rootRank,
     $                                          MPI_COMM_WORLD,ierr)
      endif
      CALL MPI_Bcast(m_act, n_act, MPI_DOUBLE, rootRank,
     $                                          MPI_COMM_WORLD,ierr)

     print *, "m_act: ", m_act(1:2)
#else
      n_act = n_act_loc         ! set total number of active particles for serial runs (NOMPI)
#endif

#ifdef DEBUG
      
      WRITE(STDERR,500) myRank,n_act
 500  FORMAT('PE',I2.2,': received a total of ',I8,' active particles')
#endif


      RETURN
      END
