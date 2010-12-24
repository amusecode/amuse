      function ran2(idum)
*
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
*
      real ran2,am,eps,rnmx
*
      parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     *    ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     *    ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.-eps)
*
      integer idum2,j,k,iv,iy
*
      common /randx/ idum2,iy,iv(ntab)
*
      if(idum.le.0)then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=ntab+8,1,-1
           k=idum/iq1
           idum=ia1*(idum-k*iq1)-k*ir1
           if(idum.lt.0)idum=idum+im1
           if(j.le.ntab)iv(j)=idum
 11     continue
        iy=iv(1)
      end if

      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
*
      if(idum.lt.0)idum=idum+im1
*
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
*
      if(idum2.lt.0)idum2=idum2+im2
*
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
*
      if(iy.lt.1)iy=iy+imm1
*
      ran2=min(am*iy,rnmx)
*
      return
      end
*
*
*

