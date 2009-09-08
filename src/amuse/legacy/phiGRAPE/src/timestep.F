      SUBROUTINE timestep(imode)
C=========================================================================
C===
C===  task: compute initial time step 
C===  -----
C===
C=========================================================================
C===                                                   17/08/2006 by SH 
C=========================================================================
      IMPLICIT NONE
#include "general.inc"

C---  local vars
      INTEGER          i,ip,power,imode
      DOUBLE PRECISION a2_mod,adot2_mod,tmp,xhelp


      IF (imode.NE.1) THEN      ! get time step for all active particles

C---  loop over local particles and compute initial time step
         DO i=1,n_loc
            
            a2_mod    =   a_loc(1,i)*a_loc(1,i)
     $                  + a_loc(2,i)*a_loc(2,i)
     $                  + a_loc(3,i)*a_loc(3,i)
            adot2_mod =   adot_loc(1,i)*adot_loc(1,i)
     $                  + adot_loc(2,i)*adot_loc(2,i)
     $                  + adot_loc(3,i)*adot_loc(3,i)


C---  get time step 
            IF (adot2_mod.EQ.0.d0) THEN
               tmp = eta_s
            ELSE
               tmp = eta_s*SQRT(a2_mod/adot2_mod)
            ENDIF           

C---  set time step to power of two value
            power = INT(LOG(tmp)/log(2.d0)) - 1
            tmp   = 2.d0**power

            IF (tmp.GT.dt_max) tmp = dt_max ! limit maximum time step
            IF (tmp.LT.dt_min) THEN
               WRITE(0,500) myRank,i,tmp,dt_min,time_cur
 500           FORMAT('PE',I2.2,': Warning! time step for particle'
     $              ,I8,' smaller than dt_min',3E14.6)
               tmp = dt_min     
            ENDIF
            dt_loc(i) = tmp

         ENDDO


      ELSE                      ! get timestep for active particle 
                                ! (active = new)

C---  loop over active particles and compute initial time step
         DO i=1,n_act_loc

            ip        = ind_act_loc(i)

            a2_mod    =   a_loc(1,ip)*a_loc(1,ip)
     $                  + a_loc(2,ip)*a_loc(2,ip)
     $                  + a_loc(3,ip)*a_loc(3,ip)
            adot2_mod =   adot_loc(1,ip)*adot_loc(1,ip)
     $                  + adot_loc(2,ip)*adot_loc(2,ip)
     $                  + adot_loc(3,ip)*adot_loc(3,ip)


C---  get time step 
            IF (adot2_mod.EQ.0.d0) THEN
               tmp = eta_s
            ELSE
               tmp = eta_s*SQRT(a2_mod/adot2_mod)
            ENDIF           

C---  set time step to power of two value
            power = INT(LOG(tmp)/log(2.d0)) - 1
            tmp   = 2.d0**power

            IF (tmp.GT.dt_max) tmp = dt_max ! limit maximum time step
            IF (tmp.LT.dt_min) THEN
               WRITE(0,500) myRank,i,tmp,dt_min,time_cur
               tmp = dt_min     
            ENDIF

            xhelp = tmp
 100        CONTINUE            ! check if timestep is good for particle time
            tmp   = xhelp
            xhelp = 0.5d0*xhelp
            IF (tmp*DBLE(INT(t_loc(ip)/tmp)).NE.t_loc(ip)) GOTO 100

C---  setting time step
            dt_loc(ip) = tmp

         ENDDO
      ENDIF

      RETURN
      END
