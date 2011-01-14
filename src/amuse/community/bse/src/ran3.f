***
      REAL FUNCTION ran3(IDUM)
*
* Random number generator from Numerical Recipes, Press et al. pg 272.
*
      IMPLICIT NONE
      INTEGER j,k,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      PARAMETER(im1=2147483563,im2=2147483399,ia1=40014,ia2=40692)
      PARAMETER(iq1=53668,iq2=52774,ir1=12211,ir2=3791,ntab=32)
      INTEGER idum
      INTEGER idum2,iy,ir(ntab)
      COMMON /RAND3/ idum2,iy,ir
      DATA idum2/123456789/, iy/0/, ir/ntab*0/
      REAL am
*
      am = 1.0/float(im1)
      imm1 = im1 - 1
      ndiv = 1 + imm1/ntab
*
      if(idum.le.0)then
         idum = MAX(-idum,1)
         idum2 = idum
         do 11 , j = ntab+8,1,-1
            k = idum/iq1
            idum = ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum = idum + im1
            if(j.le.ntab) ir(j) = idum
 11      continue
         iy = ir(1)
      endif
      k = idum/iq1
      idum = ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum = idum + im1
      k = idum2/iq2
      idum2 = ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2 = idum2 + im2
      j = 1 + iy/ndiv
      iy = ir(j) - idum2
      ir(j) = idum
      if(iy.lt.1) iy = iy + imm1
      ran3 = am*iy
*
      RETURN
      END
***
