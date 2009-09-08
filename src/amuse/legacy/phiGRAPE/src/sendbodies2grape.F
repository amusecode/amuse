      SUBROUTINE sendbodies2grape
C===  init the local GRAPE's
      IMPLICIT NONE
#include "general.inc"

C---  local vars
      INTEGER          j,k
      DOUBLE PRECISION a2by18(NDIM),a1by6(NDIM),aby2(NDIM)

C---  initial load the particles to the local GRAPE's 
      DO k=1,NDIM
         a2by18(k) = 0.d0
         a1by6(k)  = 0.d0 
         aby2(k)   = 0.d0
      ENDDO

      DO j=nbh+1,n_loc         

C---  in case of a restart acceleration and jerk are known
         IF (irestart.EQ.1) THEN
            DO k=1,NDIM
               a1by6(k)  = adot_loc(k,j)*over6 
               aby2(k)   = a_loc(k,j)*over2
            ENDDO
         ENDIF

         CALL g6_set_j_particle(clusterid, j-1-nbh, ind_loc(j), 
     $        t_loc(j), 
     $        dt_loc(j),m_loc(j),a2by18(1),a1by6(1), 
     $        aby2(1),v_loc(1,j), x_loc(1,j))
         
      ENDDO

#ifdef DEBUG
      WRITE(STDERR,600) myRank,n_loc
 600  FORMAT('PE',I2.2,':',I10,' particles loaded to local GRAPE')
      CALL flush(STDERR)
#endif

      RETURN
      END
