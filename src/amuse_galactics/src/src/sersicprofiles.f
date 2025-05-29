      function sersicpot(rad)

      include 'commonblocks'

      real L1,L2
      external gammq,gammqm1,gammln

      if(rad.eq.0.) then
         sersicpot = v0bulge**2.
         return
      endif

      u = rad/Re
      un = u**(1./nnn)
      aaa = nnn*(3.-ppp)
      L1 = Rho0*Re*Re*nnn*butt**(nnn*(ppp-2.))*
     :     gammq(nnn*(2.-ppp),butt*un)*exp(gammln(nnn*(2.-ppp)))
      L2 = Rho0*(Re**3.)*nnn*(butt**(-aaa))*exp(gammln(aaa))
      if(aaa+1.gt.butt*un) then
         L2 = L2*gammqm1(aaa,butt*un)
      else
         L2 = L2*(1.-gammq(aaa,butt*un))
      endif
      sersicpot=4.*pi*(L2/rad + L1)

      return
      END

      subroutine setsersicparameters

      include 'commonblocks'
      common /MOUSE/ comnen
      EXTERNAL funcbutt
     
      comnen = nnn
      Re = abulge
      if (nnn.gt.0.50) then
         butt1=0.6*(2.0*nnn)
      elseif (nnn.le.0.50) then
         butt1=1.0E-4
      endif
      butt2=1.20*(2.0*nnn)
      buttacc=1.0E-4

      butt=RTBISS(funcbutt,butt1,butt2,buttacc)
c      ppp=1.0d0 - 0.6097d0/nnn + 0.05463d0/(nnn*nnn)

      Rho0 = (v0bulge**2.)/(4.*pi*Re*Re*nnn*butt**(nnn*(ppp-2.))*
     *     exp( gammln(nnn*(2.-ppp))))

      return
      end

      function sersicmass(rad)

      include 'commonblocks'

      u = rad/Re
      un = u**(1./nnn)
      aaa = nnn*(3.-ppp)
      if(aaa+1.gt.butt*un) then
         gm = gammqm1(aaa,butt*un)
      else
         gm = 1.-gammq(aaa,butt*un)
      endif
      sersicmass = 4.*pi*Rho0*(Re**3.)*nnn*(butt**(nnn*(ppp-3.)))*
     +     gm*exp(gammln(nnn*(3.-ppp)))
      
      return
      end

      function sersicdens(rad)

      include 'commonblocks'

      u = rad/Re
      un = u**(1./nnn)
      sersicdens = Rho0*(u**(-ppp))*exp(-butt*un)

      return
      end

      function sersicdensprime(rad)

      include 'commonblocks'

      u = rad/Re
      un = u**(1./nnn)
      sersicdensprime = -sersicdens(rad)/Re*(ppp*nnn+butt*un)/nnn/u

      return
      end

      function sersicdens2prime(rad)

      include 'commonblocks'

      u = rad/Re
      un = u**(1./nnn)
      sersicdens2prime = sersicdens(rad)/Re/Re*(
     :     (ppp*nnn)**2.+ppp*nnn*nnn+2.*ppp*butt*nnn*un
     :     +butt*un*(nnn-1)+(butt*un)**2.)/((nnn*u)**2.)

      return
      end

      function sersicforce(rad)

      include 'commonblocks'
      real L2

      u = rad/Re
      un = u**(1./nnn)
      aaa = nnn*(3.-ppp)

      L2 = Rho0*(Re**3.)*nnn*(butt**(-aaa))*exp(gammln(aaa))
      if(aaa+1.gt.butt*un) then
         L2 = L2*gammqm1(aaa,butt*un)
      else
         L2 = L2*(1.-gammq(aaa,butt*un))
      endif

      sersicforce = -4.*pi*(L2/rad/rad)

      return
      end
C     ------------------------------------------------------------------
       real function funcbutt(b)
       real b,abe,comnen
       external gammp
       common /MOUSE/ comnen
          abe=2.*comnen
          funcbutt=( gammp(abe,b)-0.5 )
       return
       end

       real FUNCTION rtbiss(funcbutt,x1,x2,xacc)
       INTEGER JMAX
       real x1,x2,xacc,lick
       EXTERNAL funcbutt
       common /MOUSE/ comnen
       PARAMETER (JMAX=100)
       INTEGER j
       real dx,f,fmid,xmid
       lick=comnen

c      write(6,*) 'b: ',2.*lick,x2

       fmid=gammp(2.*lick,x2)-0.5
       f=gammp(2.*lick,x1)-0.5

C      write(6,*) 'c:  F-s: ',fmid,f

       if(f*fmid.ge.0.) pause 'root must be bracketed in rtbiss'
       if(f.lt.0.)then
         rtbiss=x1
         dx=x2-x1
       else
         rtbiss=x2
         dx=x1-x2
       endif
       do 11 j=1,JMAX
         dx=dx*.5
         xmid=rtbiss+dx
         fmid=gammp(2.*lick,xmid)-0.5
         if(fmid.le.0.)rtbiss=xmid
         if(abs(dx).lt.xacc .or. fmid.eq.0.) return
11     continue
       pause 'too many bisections in rtbiss'
       END

      FUNCTION FACTLN (N)
C
      INTEGER  N
      REAL     A(100),FACTLN,GAMMLN
C
      DATA  A/100*-1.0/
C
      IF (N.LT.0) PAUSE 'negative factorial'
      IF (N.LE.99) THEN
        IF (A(N+1).LT.0.0) A(N+1)=GAMMLN(float(N)+1.0)
        FACTLN=A(N+1)
      ELSE
        FACTLN=GAMMLN(float(N)+1.0)
      ENDIF
      RETURN
      END
C
C ------------------------------------------------------------------------------
C
      FUNCTION FACTRL (N)
C
      INTEGER  N,NTOP,J
      REAL     A(33),FACTRL,GAMMLN
C
      DATA  NTOP,A(1)/0,1.0/
C
      IF (N.LT.0) THEN
        PAUSE 'negative factorial'
      ELSE IF (N.LE.NTOP) THEN
        FACTRL=A(N+1)
      ELSE IF (N.LE.32) THEN
        DO 11 J=NTOP+1,N
          A(J+1)=J*A(J)
11      CONTINUE
        NTOP=N
        FACTRL=A(N+1)
      ELSE
        FACTRL=EXP(GAMMLN(float(N)+1.0))
      ENDIF
      RETURN
      END
C
C ------------------------------------------------------------------------------
C
      FUNCTION GAMMLN (XX)
C
      INTEGER  J
      real GAMMLN
      real XX
      real   COF(6),STP,HALF,ONE,FPF,X,TMP,SER
C
      DATA  COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA  HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
C
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
C
C ------------------------------------------------------------------------------
C
      FUNCTION GAMMP (A,X)
C
      real  GAMMP,A,X,GAMMCF,GLN
C
      IF (X.LT.0.0.OR.A.LE.0.0) PAUSE 'neg argument of gammp'
      IF (X.LT.A+1.0) THEN
        CALL GSER (GAMMP,A,X,GLN)
      ELSE
        CALL GCF (GAMMCF,A,X,GLN)
        GAMMP=1.0-GAMMCF
      ENDIF
      RETURN
      END
C
C ------------------------------------------------------------------------------
C
      FUNCTION GAMMQ (A,X)
C
      real GAMMQ,X,A,GAMSER,GLN
C

      IF (X.LT.0.0.OR.A.LE.0.0) PAUSE 'neg argument of gammq'
      IF (X.LT.A+1.0) THEN
        CALL GSER (GAMSER,A,X,GLN)
        GAMMQ=1.0-GAMSER
      ELSE
        CALL GCF (GAMMQ,A,X,GLN)
      ENDIF
      RETURN
      END
C

      FUNCTION GAMMQm1 (A,X)
C
      real  GAMMQm1,X,A,GAMSER,GLN
C
      IF (X.LT.0.0.OR.A.LE.0.0) PAUSE 'neg argument of gammqm1'
      CALL GSER (GAMSER,A,X,GLN)
      GAMMQm1=GAMSER
      RETURN
      END
C

C ------------------------------------------------------------------------------
C
      SUBROUTINE GCF (GAMMCF,A,X,GLN)
C
      INTEGER    ITMAX,N
      real GAMMCF,A,X,GLN,EPS,GOLD,A0,A1,B0,B1,FAC,AN
      real ANA,ANF,G,GAMMLN
      PARAMETER  (ITMAX=200,EPS=3.E-8)
C

      if(x.gt.100.) then
         gammcf = 0.
         return
      endif
      
      GLN=GAMMLN(A)
      GOLD=0.0
      A0=1.0
      A1=X
      B0=0.0
      B1=1.0
      FAC=1.0
      DO 11 N=1,ITMAX
        AN=FLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF (A1.NE.0.0) THEN
          FAC=1.0/A1
          G=B1*FAC
          IF (ABS((G-GOLD)/G).LT.EPS) GOTO 1
          GOLD=G
        ENDIF
11    CONTINUE
      gammcf = 0.
      return
1     GAMMCF=EXP(-X+A*log(X)-GLN)*G
      
      RETURN
      END
C
C ------------------------------------------------------------------------------
C
      SUBROUTINE GSER (GAMSER,A,X,GLN)
C
      INTEGER    ITMAX,N
      real       GAMSER,A,X,GLN,AP,SUM,DEL,EPS,GAMMLN
      PARAMETER  (ITMAX=100,EPS=3.E-30)
C
      GLN=GAMMLN(A)
      IF (X.LE.0.0) THEN
        IF (X.LT.0.0) PAUSE
        GAMSER=0.0
        RETURN
      ENDIF
      AP=A
      SUM=1.0/A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.0
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF (ABS(DEL).LT.ABS(SUM)*EPS) GOTO 1
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
      stop
1     GAMSER=SUM*EXP(-X+A*log(X)-GLN)
      RETURN
      END
C
C ------------------------------------------------------------------------------
C
      FUNCTION PLGNDR (L,M,X)
C
      INTEGER  L,M,I,LL
      REAL     X,PMM,SOMX2,FACT,PLGNDR,PMMP1,PLL
C
      IF (M.LT.0.OR.M.GT.L.OR.ABS(X).GT.1.) PAUSE 'bad arguments'
      PMM=1.0
      IF (M.GT.0) THEN
        SOMX2=SQRT((1.0-X)*(1.0+X))
        FACT=1.0
        DO 11 I=1,M
          PMM=-PMM*FACT*SOMX2
          FACT=FACT+2.0
11      CONTINUE
      ENDIF
      IF (L.EQ.M) THEN
        PLGNDR=PMM
      ELSE
        PMMP1=X*(2*M+1)*PMM
        IF (L.EQ.M+1) THEN
          PLGNDR=PMMP1
        ELSE
          DO 12 LL=M+2,L
            PLL=(X*(2*LL-1)*PMMP1-(LL+M-1)*PMM)/(LL-M)
            PMM=PMMP1
            PMMP1=PLL
12        CONTINUE
          PLGNDR=PLL
        ENDIF
      ENDIF
      RETURN
      END

      subroutine findbrackets(func,psi,rmin,rmax)

      real rmin,rmax,psi,psirmin,psirmax

      rmin = 1.
      rmax = 1.

 77   psirmin = func(rmin)
      if(psirmin.le.psi) then
         rmin = rmin*0.1
         goto 77
      endif
 88   psirmax = func(rmax)
      if(psirmax.ge.psi) then
         rmax = rmax*10.
         goto 88
      endif
      return
      end


      subroutine RTBIS (func,psi,X1,X2,XACC,XBEST)
C
      INTEGER    JMAX,J
      PARAMETER  (JMAX=100)
      real X1,X2,XACC,FMID,F,DX,XMID,psi,XBEST
C
      external func

      fmid = func(x2)-psi
      f = func(x1)-psi
      if(f*fmid.ge.0.0) write(*,*) psi,x2,fmid,x1,f
      IF (F*FMID.GE.0.0) PAUSE 'Root must be bracketed for bisection.'      
      IF (F.LT.0.0) THEN
        XBEST=X1
        DX=X2-X1
      ELSE
        XBEST=X2
        DX=X1-X2
      ENDIF
      DO 11 J=1,JMAX
        DX=DX*0.5
        XMID=XBEST+DX
        fmid = func(xmid) - psi
        IF (FMID.LT.0.0) XBEST=XMID
        IF ((ABS(DX).LT.XACC).OR.(FMID.EQ.0.0)) RETURN
11    CONTINUE
      PAUSE 'too many bisections'
      END







