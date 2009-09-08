      SUBROUTINE predictor(imode)
C=========================================================================
C===  predict active particles
C===
C===    imode = 0 true prediction
C===          = 1 no prediction, just copy current pos and vel
C=========================================================================
      IMPLICIT NONE
#include "general.inc"
      
C---  local vars
      INTEGER i,k,k3,imode
      DOUBLE PRECISION dt2half,dt3over6

C---  dont do anything if none of the local particles are active
      IF (n_act_loc.LT.1) RETURN

      IF (imode.NE.1) THEN

         DO i=1,n_act_loc
         
            dt2half  = 0.5d0 * dt_act(i) * dt_act(i)
            dt3over6 = dt_act(i) * dt2half / 3.d0

            DO k=1,3
               k3 = k + 3
               body_pred_loc(k,i)  = x_act(k,i) + v_act(k,i)*dt_act(i)  
     $                                          + a_act(k,i)*dt2half
     $                                          + adot_act(k,i)*dt3over6
               body_pred_loc(k3,i) = v_act(k,i) + a_act(k,i)*dt_act(i)
     $                                          + adot_act(k,i)*dt2half

            ENDDO

         ENDDO

      ELSE

         DO i=1,n_act_loc
            DO k=1,NDIM
               body_pred_loc(k,i)   = x_act(k,i) ! set x_act_new and v_act_new here
               body_pred_loc(k+3,i) = v_act(k,i) ! since predictor is not called
            ENDDO
         ENDDO

      ENDIF


      RETURN
      END
