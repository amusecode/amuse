cProgram to project a snap12.dat file (which must be so named)
      integer i,k,idd,ik1,ik2,ikb,ikind,nwhich,ibstra,iex
      real body,temp,sm1,sm2,slum,slum1,rad1,rad2,r,vr,vt,ap,
     &     teff,epoch1,epoch2,mv,mbv,mi,mv1,mbv1,mi1,mv2,
     &     mbv2,mi2,spin1,spin2,ecc,rp
               
      open (7,file=trim(datadir)//'/snap12.dat')
      open (8,file=trim(datadir)//'/snap12p.dat')
      idum = 2
 10   continue
*
      read(7,100,end=20) i,k,idd,body,temp,sm1,sm2,slum,slum1,rad1,rad2,
     &            r,vr,vt,ap,ik1,ik2,ikb,ikind,nwhich,teff,epoch1,
     &            epoch2,ibstra,mv,mbv,mi,mv1,mbv1,mi1,mv2,mbv2,mi2,
     &            spin1,spin2,ecc,iex
 100  format(1x,3i8,1p12e12.4,4i5,i8,1p3e12.4,i5,1p12e12.4,i6)
*                                   
      cosphi = ran2(idum)
      rp = r*sqrt(1.0-cosphi**2)
*
      write(8,100) i,k,idd,body,temp,sm1,sm2,slum,slum1,rad1,rad2,
     &            rp,vr,vt,ap,ik1,ik2,ikb,ikind,nwhich,teff,epoch1,
     &            epoch2,ibstra,mv,mbv,mi,mv1,mbv1,mi1,mv2,mbv2,mi2,
     &            spin1,spin2,ecc,iex
*                    
      go to 10
 20   continue
      close (7)
      close (8)
      stop
      end

      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]vrt43D04.
