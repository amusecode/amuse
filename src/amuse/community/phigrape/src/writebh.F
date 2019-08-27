      SUBROUTINE writebh(time,ipred)
C=========================================================================
C===
C===  task: write detailed output for BHs at given time, ipred=1 
C===  ----- BH data is predicted if t_BH!=time and ipred=1
C===
C=========================================================================
      IMPLICIT NONE
#include "general.inc"

      INTEGER          ipred,i,k
      DOUBLE PRECISION time
      DOUBLE PRECISION dtbh(NBHMAX),xbh(NDIM,NBHMAX),vbh(NDIM,NBHMAX),
     $                 dt2half,dt3over6

C---  get current time difference of BHs
      IF (ipred.EQ.1) THEN 

         DO i=1,nbh
            dtbh(i) = time - t_loc(i)       ! can use t_loc() here because
                                            ! there are only a few BH particles
                                            ! that always should be local on root

C---  predict BH at current time
            dt2half  = 0.5d0 * dtbh(i)*dtbh(i)
            dt3over6 = dtbh(i)*dt2half / 3.d0

            DO k=1,3
               xbh(k,i) = x_loc(k,i) + v_loc(k,i)*dtbh(i)
     $                               + a_loc(k,i)*dt2half
     $                               + adot_loc(k,i)*dt3over6
               vbh(k,i) = v_loc(k,i) + a_loc(k,i)*dtbh(i)
     $                               + adot_loc(k,i)*dt2half

            ENDDO
         ENDDO
      ELSE

C---  just copy current position
         DO i=1,nbh
            DO k=1,3
               xbh(k,i) = x_loc(k,i)
               vbh(k,i) = v_loc(k,i)
            ENDDO
         ENDDO

      ENDIF

      DO i=1,nbh
         WRITE(ubh,504) ind_loc(i)-1,time,m_loc(i),
     $                 (xbh(k,i),k=1,3),(vbh(k,i),k=1,3)
      ENDDO
 504  FORMAT(I10.10,1P,E15.6,2X,E15.6,2X,3E14.6,2X,3E14.6)

      CALL flush(ubh)

      RETURN
      END
