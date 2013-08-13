***
      SUBROUTINE MAGBRK(kw,m1,menv,r,ospin,djmb)
      implicit none
      integer kw
      real*8 m1,menv,r,ospin,djmb
*
* Calculate spin changes owing to magnetic braking.
*
      djmb = 0.d0
      if(kw.lt.10.and.m1.ge.0.35d0)then
         djmb = 5.83d-16*(menv/m1)*(r*ospin)**3
      endif
*
      RETURN
      END
***
