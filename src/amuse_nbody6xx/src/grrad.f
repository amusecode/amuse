***
      SUBROUTINE GRRAD(m1,m2,sep,ecc,jorb,djgr,de)
      implicit none
      real*8 m1,m2,sep,ecc,jorb,djgr,de,acrit
      real*8 ecc2,sqome2,sqome5,f1
*
* Include gravitational radiation for very close or massive systems.
*
      djgr = 0.d0
      de = 0.d0
      ACRIT = 10.0
      IF (m1 + m2.GT.10.0) ACRIT = 30.0
      if(sep.le.ACRIT)then
         djgr = 8.315d-10*m1*m2*(m1+m2)/(sep*sep*sep*sep)
         ecc2 = ecc*ecc
         sqome2 = SQRT(1.d0-ecc2)
         sqome5 = sqome2**5
         f1 = (19.d0/6.d0) + (121.d0/96.d0)*ecc2
         de = djgr*ecc*f1/sqome5
         djgr = djgr*jorb*(1.d0+0.875d0*ecc2)/sqome5
      endif
*
      RETURN
      END
***
