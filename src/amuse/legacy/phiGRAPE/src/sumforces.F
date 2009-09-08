
      SUBROUTINE sumforces(ifirst)
C=========================================================================
C===
C===  task:  using allreduce sum up all the partial local forces
C===  -----
C===
C=========================================================================
      IMPLICIT NONE

C---  input vars
      INTEGER ifirst,i,ip

#include "general.inc"
#ifndef NOMPI
#include <mpif.h>
      INCLUDE 'mpi.inc'

C---  local vars
      INTEGER ierr
#endif
#ifdef NOMPI
C---  equivalence here!!!
#endif

      IF (ifirst.EQ.1) THEN     ! initial force calculation

#ifndef NOMPI
C--------------------------------------------------------------------
C---  collect all forces on root and then scatter to all PE's
C--------------------------------------------------------------------
         CALL MPI_Reduce(body_act_tmp, body_act_new, 7*n_act,MPI_DOUBLE, 
     $                   MPI_SUM, rootRank, MPI_COMM_WORLD,ierr)
         IF (myRank.EQ.rootRank) THEN
            DO i=1,n_act        ! copy forces etc. to main arrays
               pot(i)    = body_act_new(ipnew,i)
               a(1,i)    = body_act_new(iaxnew,i)
               a(2,i)    = body_act_new(iaynew,i)
               a(3,i)    = body_act_new(iaznew,i)
               adot(1,i) = body_act_new(ijxnew,i)
               adot(2,i) = body_act_new(ijynew,i)
               adot(3,i) = body_act_new(ijznew,i)
            ENDDO
         ENDIF
          
         CALL MPI_Scatter(pot,    n_loc, MPI_DOUBLE, pot_loc,    n_loc,
     $                    MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
         CALL MPI_Scatter(a,    3*n_loc, MPI_DOUBLE, a_loc,    3*n_loc,
     $                    MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)
         CALL MPI_Scatter(adot, 3*n_loc, MPI_DOUBLE, adot_loc, 3*n_loc,
     $                    MPI_DOUBLE, rootRank, MPI_COMM_WORLD,ierr)

         CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
C             -----------
#endif
#ifdef NOMPI
         DO i=1,n_act           ! just copy forces etc., *_loc is set by EQUIVALENCE
            ip         = ind_act_loc(i)
            pot(ip)    = body_act_tmp(ipnew,i)
            a(1,ip)    = body_act_tmp(iaxnew,i)
            a(2,ip)    = body_act_tmp(iaynew,i)
            a(3,ip)    = body_act_tmp(iaznew,i)
            adot(1,ip) = body_act_tmp(ijxnew,i)
            adot(2,ip) = body_act_tmp(ijynew,i)
            adot(3,ip) = body_act_tmp(ijznew,i)
         ENDDO
#endif

      ELSE                      ! after force calculation in main loop

#ifndef NOMPI         
C---  Reduce the "global" vectors from "local" 
         CALL MPI_Allreduce(body_act_tmp, body_act_new, 7*n_act,
     $                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif


      ENDIF

      RETURN
      END
