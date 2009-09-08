      SUBROUTINE sync
C===  particles previouisly selected are sync'd by a full GRAPE integration
      IMPLICIT NONE
#include "general.inc"
      
      INTEGER          i,power
      DOUBLE PRECISION dtcmax,ttmp

      IF (n_act_loc.LE.0) RETURN ! no sync'ing needed

C--------------------------------------------------------------------
C---  predictor (for local active particles)
C--------------------------------------------------------------------
         CALL predictor(0)
C             ---------         

C--------------------------------------------------------------------
C---  gather active particles on all PE's
C--------------------------------------------------------------------
         CALL allgather_ap
C             ------------

C--------------------------------------------------------------------
C---  get local forces on all PE's
C---  (as this is not the first call set ifirst to 0 in this call)
C--------------------------------------------------------------------
         CALL gravity(time_cur,0)
C             -------

C--------------------------------------------------------------------
C---  sum up the partial forces
C--------------------------------------------------------------------
         CALL sumforces(0)      ! call with 0; just sum forces
C             ---------

C--------------------------------------------------------------------
C---  make the corrector step
C---  this also updates the time steps of active particles
C--------------------------------------------------------------------
         CALL corrector
C             ---------


C--------------------------------------------------------------------
C---  need to correct time step to a power of two value
C--------------------------------------------------------------------

C---  first get maximal allowed time step wrt to current time
         ttmp = time_cur
 100     CONTINUE
            power     = INT(LOG(ttmp)/LOG(2.d0))
            dtcmax    = 2.d0**power
            IF (dtcmax.GT.ttmp) dtcmax = dtcmax/2.d0
            ttmp = dmod(time_cur,dtcmax)
         IF (ttmp.GT.0.d0) GOTO 100

C---  now correct individual time step to powers of two or dtcmax
         DO i=1,n_act_loc

            power     = INT(LOG(dt_act(i))/LOG(2.d0)) - 1
            dt_act(i) = MIN(2.d0**power, dtcmax)

         ENDDO
 
C--------------------------------------------------------------------
C---  update local particles from local active particles
C--------------------------------------------------------------------
         CALL update_loc_p
C             ------------

C--------------------------------------------------------------------
C---  also update local particles on GRAPE from local active particles
C--------------------------------------------------------------------
         CALL update_grape(0)
C             ------------


      RETURN
      END
