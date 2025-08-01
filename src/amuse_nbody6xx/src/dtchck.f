***
      SUBROUTINE dtchck(time,dt,dtmin)
      IMPLICIT NONE
      REAL*8 time,dt,dtmin
*
* Find the largest block time-step, starting with dt, 
* that is commensurate with current time.
*
      dt = MAX(dt,dtmin)
 10   continue
      if(DMOD(time,dt).ne.0.d0)then
         dt = 0.5d0*dt
         if(dt.ge.dtmin) goto 10
         dt = dtmin
      endif
*
      RETURN
      END
***
