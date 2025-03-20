***
      SUBROUTINE deltat(kw,age,tm,tn,tscls,dt,dtr)
      implicit none
*
      INTEGER kw
      REAL*8 age,tm,tn,tscls(20)
      REAL*8 dt,dtr
      REAL*8 pts1,pts2,pts3
      COMMON /POINTS/ pts1,pts2,pts3
*
*     Base new time scale for changes in radius & mass on stellar type.
*
      if(kw.le.1)then
         dt = pts1*tm
         dtr = tm - age
      elseif(kw.eq.2)then
         dt = pts1*(tscls(1) - tm)
         dtr = tscls(1) - age
      elseif(kw.eq.3)then
         if(age.lt.tscls(6))then
            dt = pts2*(tscls(4) - age)
         else
            dt = pts2*(tscls(5) - age)
         endif
         dtr = MIN(tscls(2),tn) - age
      elseif(kw.eq.4)then
         dt = pts2*tscls(3)
         dtr = MIN(tn,tscls(2) + tscls(3)) - age
      elseif(kw.eq.5)then
         if(age.lt.tscls(9))then
            dt = pts3*(tscls(7) - age)
         else
            dt = pts3*(tscls(8) - age)
         endif
         dtr = MIN(tn,tscls(13)) - age
      elseif(kw.eq.6)then
         if(age.lt.tscls(12))then
            dt = pts3*(tscls(10) - age)
         else
            dt = pts3*(tscls(11) - age)
         endif
         dt = MIN(dt,0.005d0)
         dtr = tn - age
      elseif(kw.eq.7)then
         dt = pts1*tm
         dtr = tm - age
      elseif(kw.eq.8.or.kw.eq.9)then
         if(age.lt.tscls(6))then
            dt = pts2*(tscls(4) - age)
         else
            dt = pts2*(tscls(5) - age)
         endif
         dtr = tn - age
      else
*        dt = MAX(0.1d0,age*10.d0)
         dt = MAX(0.1d0,dt*10.d0)
         dt = MIN(dt,5.0d+02)
         dtr = dt
      endif
*
      RETURN
      END
***
