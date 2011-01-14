      SUBROUTINE predict_potential(tnow)
C===  find particles not at current time to compute potential
C===  for predicted positions 
      IMPLICIT NONE
#include "general.inc"
      
      INTEGER          i,power
      DOUBLE PRECISION dtcmax,ttmp,tnow


C--------------------------------------------------------------------
C---  select particles not at current time
C--------------------------------------------------------------------
         CALL selectsync(tnow)
C             ----------         


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
         CALL gravity(tnow,0)
C             -------

C--------------------------------------------------------------------
C---  sum up the partial forces
C--------------------------------------------------------------------
         CALL sumforces(0)      ! call with 0; just sum forces
C             ---------



      RETURN
      END
