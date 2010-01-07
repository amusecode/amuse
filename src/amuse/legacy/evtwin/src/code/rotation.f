! Functions for computing corrections to the stellar sructure equations
!  after Endal & Sofia, 1976, ApJ 210:184.
! Based on code by Alexander Heger

c POTENTIAL: HAUPTROUTINE, BENOETIGT DIE VERTEILUNG VON DICHTE, MASSE
c            UND WINKELGESCHWINDIGKEIT SOWIE DIE RADIUSKOORDINATE
c            LIEFERT DIE FELDER FT UND FP
c RADAU: LIEFERT DEN WERT DER IN RADGL ABGELEGTEN RADAU-GLEICHUNG (ETA)
c        BENOETIGT DIE RADIUSKOORDINATE UND RHO/<RHO>
c NULLST: BERECHNET NACH VORGABE VON RPSI RNULL
c PSIINT: BERECHNET DAS ZUM POTENTIAL GEHOERIGE INTEGRAL
c GPSI: BERECHNET G ALS FUNKTION VON THETA

! compute correction factors ft and fp to the structure equations
! FIXME: think of a more sensible name than POTENTIAL
! Input variables:
!  N     Number of meshpoints used by the model (arrays may be larger)
!  XM    Mass coordinate [gram or 10^33 gram (both work)]
!  R     Radius coordinate [cm or 10^11 cm (both work)]
!  RHO   Density [gram/cm^3]
!  AW    Angular velocity [rad/sec]
! Output variables:
!  FTOUT Correction factor FT at each meshpoint
!  FPOUT Correction factor FP at each meshpoint
      SUBROUTINE POTENTIAL(N, XM, R, RHO, AW, FTOUT, FPOUT)
      USE CONSTANTS
      USE MESH
      IMPLICIT REAL*8 ( A-H, O-Z)
      INTEGER, PARAMETER :: INTMAX = 101
      DOUBLE PRECISION, INTENT(OUT) :: FTOUT(NM), FPOUT(NM)
      DOUBLE PRECISION, INTENT(IN) :: AW(NM), R(NM), RHO(NM), XM(NM)
      INTEGER, INTENT(IN) :: N
      EXTERNAL RADGL
      EXTERNAL RKQB
      EXTERNAL APOT

      COMMON/POTI/ETA(NM),RNULL(NM),GEFF(NM),SPGP(NM),SPGM(NM),
     &GP(INTMAX),GMH(INTMAX),FT(NM),FP(NM)
      COMMON/RADI/RAE(INTMAX),ROV,ETANU,WR,RPS,XMPS,A
      COMMON/MUS/GUR(NM)
      COMMON/WINKEL/THETA(INTMAX),COST(INTMAX),SINT(INTMAX),
     &     COST2(INTMAX),SINT2(INTMAX),UTHETA(INTMAX)
      COMMON /FXLIM/
     &     FTMIN,FPMIN,FTWARN,FPWARN
      COMMON/PATH/KMAX, KOUNT,DXSAV, XP(200), YP(50,200)
      DIMENSION VETA(1)
      
      LOGICAL, SAVE :: lfirst  = .TRUE.
 
      KMAX = 0
      DXSAV = 0.D0

      FT_MIN=1.0D0
      FP_MIN=1.0D0
      IFTLO=0
      IFTHI=0
      IFPLO=0
      IFPHI=0

      IF (lfirst) THEN
         lfirst = .FALSE.
         DTHETA=CPI / (2.D0 * (INTMAX-1))
         DTHETA2=0.5D0*DTHETA
         DO I=1,INTMAX
            THETA(i)=(i-1)*DTHETA
            COST(i)=COS(THETA(i))
            COST2(i)=COST(i)*COST(i)
            SINT(i)=SIN(THETA(i))
            SINT2(i)=SINT(i)*SINT(i)
            UTHETA(i)=0.5D0*(3.D0*COST2(i)-1.D0)
         ENDDO
      ENDIF
      
c MAIN LOOP
      ETA(1)=-1.00D-20
      TEGRA=0.D0
      RNULL(1)=0.D0

      RNORM=1.D0
      DO I=1,N-1,1
! For each shell, compute the distortion (in terms of R0) of each shell. The
! inner edge of the shell is at I, the outer edge at (I+1). The density is
! defined in the interior of each shell and needs to be interpolated when
! it is needed at the edge.
c Compute mean density at this point
         ROM1=XM(I+1)*0.75D0/(CPI*R(I+1)**3)
         ROV=RHO(I)/ROM1
         VETA(1)=ETA(I)
         RI=R(I)
         IF (I.EQ.1) RI=1.D-10
         RIP=R(I+1)
         RI=LOG(RI)
         RIP=LOG(RIP)
         KANZ=1
         DAT=1.D-04
         DOT=(RIP-RI)/8.D0
         DUT=(RIP-RI)/300.D0

         CALL ODEINT(VETA,KANZ,RI,RIP,DAT,DOT,DUT,NOK,NBAD,
     &               RADGL,RKQB,I)
         ETA(I+1)=VETA(1)

! Compute R0
! Parameters are passed through the COMMON block /RADI/

         XMPS=XM(I+1)
         RPS=R(I+1)
         WR=AW(I+1)
         ETANU=VETA(1)
         R1=2.00D0*R(I+1)
         R2=0.01D0*R(I+1)
         RNUV1=RTSAFE(APOT,R2,R1,1.0D-8,I)

         RNUV2=RNUV1
         RNULL(I+1) =  RNUV2

! Compute integral for the potential calculation. See Endal&Sofia for details

         TEGRA= PSIINT(RHO(I+1),RHO(I),XM(I+1),AW(I+1),ETA(I+1),
     &        RNULL(I), RNULL(I+1)) + TEGRA

! Calculate g and 1/g on a quarter circle.
         GUR(I+1)=0.0D0
         DO K=1,INTMAX,1
            GP(K)=GPSI(RNULL(I+1), TEGRA ,K )
            GMH(K)=1.D0 / GP(K)
            GUR(I+1)=GUR(I+1)+GP(K)*DTHETA
         END DO
         GUR(I+1)=GUR(I+1)/CPI*2.D0

! Find SPSI*<g> and SPSI*<1/g>. SPSI is the surface area of the equipotential
         RRS=RAE(1)**2*SINT(1)
         RRS1=RAE(INTMAX)**2*SINT(INTMAX)
         SPGP(I+1)=(GP (1)*RRS+GP (INTMAX)*RRS1)*DTHETA2
         SPGM(I+1)=(GMh(1)*RRS+GMh(INTMAX)*RRS1)*DTHETA2
         DO K=2,INTMAX-1
            RRS=RAE(K)*RAE(K)*SINT(k)
            SPGP(I+1)=SPGP(I+1)+GP (K)*RRS*DTHETA
            SPGM(I+1)=SPGM(I+1)+GMh(K)*RRS*DTHETA
         ENDDO
         SPGP(I+1)=SPGP(I+1)*CPI4
         SPGM(I+1)=SPGM(I+1)*CPI4

c  BERECHNEN VON FP UND FT

         FP(I+1) = CPI4 * RPS**4 / (CG*XMPS*SPGM(I+1))
         FT(I+1) = (CPI4 * RPS**2)**2 / (SPGP(I+1)*SPGM(I+1))

         IF (FT(I+1) < FTwarn) THEN
            IF (iftlo == 0) iftlo=i+1
            ifthi=i+1
            ft_min=MIN(ft_min,ft(i+1))
         ELSEIF (iftlo /= 0) THEN 
            iftlo=0
            ifthi=0
            ft_min=1.0D0
         ENDIF 

         IF (FP(I+1) < FPwarn) THEN
            IF (ifplo == 0) ifplo=i+1
            ifphi=i+1
            fp_min=MIN(fp_min,fp(i+1))
         ELSEIF (ifplo /= 0) THEN 
            ifplo=0
            ifphi=0
            fp_min=1.0D0
         ENDIF
 
         FT(I+1)=MAX(FTMIN, MIN(FT(I+1),1.D0))
         FP(I+1)=MAX(FPMIN, MIN(FP(I+1),1.D0))
      ENDDO

      !FT(2)=FT(3)
      !FP(2)=FP(3)
      FT(1)=FT(2)
      FP(1)=FP(2)
      !FT(N-2)=FT(N-3)
      !FT(N-1)=FT(N-2)
      FT(N)=FT(N-1)
      !FP(N-2)=FP(N-3)
      !FP(N-1)=FP(N-2)
      FP(N)=FP(N-1)

      FTOUT(1:N)=FT(1:N)
      FPOUT(1:N)=FP(1:N)
      RETURN
      END
c***********************************************************************
c     PROGRAMTEIL RADAU
c BERECHNET DEN WERT ETA DER RADAU GLEICHUNG, DIE IN RADGL ABGELEGT IST
c RHO/<RHO> MUSS IN COMMON/RADI/ROV ABGELEGT SEIN.
c DIE RUNGE-KUTTA METHODE ZUR LOESUNG DER DIFFERENTIALGLEICHUNG STAMMT
c AUS DEN NUMERICAL RECIPES
c***********************************************************************

      SUBROUTINE RADGL(X,Y,DYDX)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (Nmax=10, intmax=101)

      DIMENSION y(Nmax),dydx(Nmax)
      COMMON/RADI/RAE(intmax),ROV,etanu,WR,RPS,XMPS,A

	x=x
c         if (x.lt.1.d-30) x=1.d-30
c         IF(Y(1).GT.1.D30) Y(1)=1.D30
c         IF (Y(1).LT.-1.D30) Y(1)=-1.D30
c        WRITE(*,*)'RADAU Y=',Y(1)
c        WRITE(*,*)'X=',X
         dydx(1)=(6.D0-6.D0*ROV*(Y(1)+1.D0)-Y(1)*(Y(1)-1.D0))
c         write(44,*)x,y(1),dydx(1)
      RETURN
      END

c-----------------------------------------------------------------------

      SUBROUTINE ODEINT(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,DERIVS
     *,RKQC,nuts)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (MAXSTP=10000,NMAX=10,TWO=2.0,ZERO=0.0,TINY=1.E-30)
      EXTERNAL DERIVS
      EXTERNAL RKQC
      COMMON /PATH/ KMAX,KOUNT,DXSAV,XP(200),YP(50,200)
      common/errlim/nrerr 
      DIMENSION YSTART(NVAR),YSCAL(NMAX),Y(NMAX),DYDX(NMAX)
      X=X1
      H=DSIGN(H1,X2-X1)
      NOK=0
      NBAD=0
      KOUNT=0
      Y(1:NVAR)=YSTART(1:NVAR)
      XSAV=X-DXSAV*TWO
      DO NSTP=1,MAXSTP
        CALL DERIVS(X,Y,DYDX)
        DO I=1,NVAR
          YSCAL(I)=DABS(Y(I))+DABS(H*DYDX(I))+TINY
        ENDDO
        IF(KMAX.GT.0)THEN
          IF(DABS(X-XSAV).GT.DABS(DXSAV)) THEN
            IF(KOUNT.LT.KMAX-1)THEN
              KOUNT=KOUNT+1
              XP(KOUNT)=X
              YP(1:NVAR,KOUNT)=Y(1:NVAR)
              XSAV=X
            ENDIF
          ENDIF
        ENDIF
        IF((X+H-X2)*(X+H-X1).GT.ZERO) H=X2-X
        CALL RKQC(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
        IF(HDID.EQ.H)THEN
          NOK=NOK+1
        ELSE
          NBAD=NBAD+1
        ENDIF
        IF((X-X2)*(X2-X1).GE.ZERO)THEN
          YSTART(1:NVAR)=Y(1:NVAR)
          IF(KMAX.NE.0)THEN
            KOUNT=KOUNT+1
            XP(KOUNT)=X
            YP(1:NVAR,KOUNT)=Y(1:NVAR)
          ENDIF
          RETURN
        ENDIF
        IF(DABS(HNEXT).LT.HMIN) hnext=hmin
        H=HNEXT
      ENDDO
      RETURN
      END

      SUBROUTINE RKQB(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=10, ONE=1.D0,SAFETY=0.9D0,ERRCON=6.D-4)
      EXTERNAL DERIVS
      DIMENSION Y(N),DYDX(N),YSCAL(N),YTEMP(NMAX),YSAV(NMAX),DYSAV(NMAX)
      PGROW=-0.20D0
      PSHRNK=-0.25D0
      FCOR=2.D-1/3.D0
      XSAV=X
      YSAV(1:N)=Y(1:N)
      DYSAV(1:N)=DYDX(1:N)
      H=HTRY
1     HH=0.5D0*H
      CALL RK4(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)
      X=XSAV+HH
      CALL DERIVS(X,YTEMP,DYDX)
      CALL RK4(YTEMP,DYDX,N,X,HH,Y,DERIVS)
      X=XSAV+H
c      IF(X.EQ.XSAV)PRINT*, 'Stepsize not significant in RKQC.'
      CALL RK4(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
      ERRMAX=0.
      DO I=1,N
        YTEMP(I)=Y(I)-YTEMP(I)
        ERRMAX=MAX(ERRMAX,ABS(YTEMP(I)/YSCAL(I)))
      ENDDO
      ERRMAX=ERRMAX/EPS
      IF(ERRMAX.GT.ONE) THEN
        H=SAFETY*H*(ERRMAX**PSHRNK)
        GOTO 1
      ELSE
        HDID=H
        IF(ERRMAX.GT.ERRCON)THEN
          HNEXT=SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT=4.D0*H
        ENDIF
      ENDIF
      Y(1:N)=Y(1:N)+YTEMP(1:N)*FCOR
      RETURN
      END

      SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT,DERIVS)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=10)
      DIMENSION Y(N),DYDX(N),YOUT(N),YT(NMAX),DYT(NMAX),DYM(NMAX)
      HH=H*0.5D0
      H6=H/6.D0
      XH=X+HH
      YT(1:N)=Y(1:N)+HH*DYDX(1:N)
      CALL DERIVS(XH,YT,DYT)
      YT(1:N)=Y(1:N)+HH*DYT(1:N)
      CALL DERIVS(XH,YT,DYM)
      DO I=1,N
        YT(I)=Y(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
      END DO
      CALL DERIVS(X+H,YT,DYT)
      YOUT(1:N)=Y(1:N)+H6*(DYDX(1:N)+DYT(1:N)+2.D0*DYM(1:N))
      RETURN
      END

      SUBROUTINE bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs)
      implicit real*8(a-h,o-z)
      INTEGER nv,NMAX,KMAXX,IMAX
      REAL*8 eps,hdid,hnext,htry,x,dydx(nv),y(nv),yscal(nv),SAFE1,SAFE2,
     *REDMAX,REDMIN,TINY,SCALMX
      PARAMETER (NMAX=50,KMAXX=8,IMAX=KMAXX+1,SAFE1=.25,SAFE2=.7,
     *REDMAX=1.e-5,REDMIN=.7,TINY=1.e-30,SCALMX=.1)
cU    USES derivs,mmid,pzextr
      INTEGER i,iq,k,kk,km,kmax,kopt,nseq(IMAX)
      REAL*8 eps1,epsold,errmax,fact,h,red,scale,work,wrkmin,xest,xnew,
     *a(IMAX),alf(KMAXX,KMAXX),err(KMAXX),yerr(NMAX),ysav(NMAX),
     *yseq(NMAX)
      LOGICAL first,reduct
      SAVE a,alf,epsold,first,kmax,kopt,nseq,xnew
      EXTERNAL derivs
      DATA first/.true./,epsold/-1./
      DATA nseq /2,4,6,8,10,12,14,16,18/
      if(eps.ne.epsold)then
        hnext=-1.d29
        xnew=-1.d29
        eps1=SAFE1*eps
        a(1)=nseq(1)+1
        do k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
        end do
        do iq=2,KMAXX
          do k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/
     &       ((a(iq+1)-a(1)+1.0D0)*(2*k+1)))
          end do
        end do
        epsold=eps
        do kopt=2,KMAXX-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt)) exit;
        end do
        kmax=kopt
      endif
      h=htry
      ysav(1:nv)=y(1:nv)
      if(h.ne.hnext.or.x.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do 17 k=1,kmax
        xnew=x+h
        if(xnew.eq.x)STOP 'step size underflow in bsstep'
        call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
        xest=(h/nseq(k))**2
        call pzextr(k,xest,yseq,y,yerr,nv)
        if(k.ne.1)then
          errmax=TINY
          do i=1,nv
            errmax=max(errmax,abs(yerr(i)/yscal(i)))
          end do
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/SAFE1)**(1.0D0/(2*km+1))
        endif
        if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          if(errmax.lt.1.)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
            red=SAFE2/err(km)
            goto 3
          else if(k.eq.kopt)then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1.d0/err(km)
              goto 3
            endif
          else if(kopt.eq.kmax)then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*SAFE2/err(km)
              goto 3
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          endif
        endif
17    continue
3     red=min(red,REDMIN)
      red=max(red,REDMAX)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.d35
      do kk=1,km
        fact=max(err(kk),SCALMX)
        work=fact*a(kk+1)
        if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
        endif
      end do
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),SCALMX)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif
      return
      END

      SUBROUTINE mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)
      implicit real*8(a-h,o-z)
      INTEGER nstep,nvar,NMAX
      REAL*8 htot,xs,dydx(nvar),y(nvar),yout(nvar)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
      INTEGER i,n
      REAL*8 h,h2,swap,x,ym(NMAX),yn(NMAX)
      h=htot/nstep
      forall (i=1:nvar)
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
      end forall
      x=xs+h
      call derivs(x,yn,yout)
      h2=2.*h
      do n=2,nstep
        do i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
        end do 
        x=x+h
        call derivs(x,yn,yout)
      end do
      do i=1,nvar
        yout(i)=0.5d0*(ym(i)+yn(i)+h*yout(i))
      end do
      return
      END

      SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv)
      implicit real*8(a-h,o-z)
      INTEGER iest,nv,IMAX,NMAX
      REAL*8 xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13,NMAX=50)
      INTEGER j,k1
      REAL*8 delta,f1,f2,q,d(NMAX),qcol(NMAX,IMAX),x(IMAX)
      SAVE qcol,x
      x(iest)=xest
      do 11 j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
11    continue
      if(iest.eq.1) then
        do 12 j=1,nv
          qcol(j,1)=yest(j)
12      continue
      else
        do 13 j=1,nv
          d(j)=yest(j)
13      continue
        do 15 k1=1,iest-1
          delta=1./(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do 14 j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
14        continue
15      continue
        do 16 j=1,nv
          qcol(j,iest)=dy(j)
16      continue
      endif
      return
      END

      SUBROUTINE rzextr(iest,xest,yest,yz,dy,nv)
      implicit real*8(a-h,o-z)
      INTEGER iest,nv,IMAX,NMAX
      REAL*8 xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13,NMAX=50)
      INTEGER j,k
      REAL*8 b,b1,c,ddy,v,yy,d(NMAX,IMAX),fx(IMAX),x(IMAX)
      SAVE d,x
      x(iest)=xest
      if(iest.eq.1) then
        do 11 j=1,nv
          yz(j)=yest(j)
          d(j,1)=yest(j)
          dy(j)=yest(j)
11      continue
      else
        do 12 k=1,iest-1
          fx(k+1)=x(iest-k)/xest
12      continue
        do 14 j=1,nv
          yy=yest(j)
          v=d(j,1)
          c=yy
          d(j,1)=yy
          do 13 k=2,iest
            b1=fx(k)*v
            b=b1-c
            if(b.ne.0.) then
              b=(c-v)/b
              ddy=c*b
              c=b1*b
            else
              ddy=v
            endif
            if (k.ne.iest) v=d(j,k)
            d(j,k)=ddy
            yy=yy+ddy
13        continue
          dy(j)=ddy
          yz(j)=yy
14      continue
      endif
      return
      END

c***********************************************************************
c   PROGRAMMTEIL NULST
c   BERECHNET RNULL AUS DEM IN APOT ABGELEGTEN ZUSAMMENHANG.
c   DIE PARAMETER SIND IN COMMON/RADI/ ABGELEGT
c   DIE KOMBINIERTE NEWTON- UND BISEKTORMETHODE ZUR NULLSTELLENSUCHE
c   STAMMT AUS DEN NUMERICAL RECIPIES
c***********************************************************************

      SUBROUTINE APOT(RNU,RFK,DRFK)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (CG=6.6726D-8, intmax=101)
      COMMON/RADI/RAE(intmax),ROV,etanu,WR,RPS,XMPS,A
      LOGICAL
     &     lfirst
      DATA
     &     lfirst /.TRUE./

c     IF (lfirst) THEN
c        lfirst=.FALSE.
         c0=1.0D0/3.0D0
         c1=0.6D0*CG
         c2=2.0D0/35.0D0
         c3=3.0D0*c2
c     ENDIF

      WRNU2=WR*WR*RNU*RNU/(c1*XMPS*(2.D0+etanu))
      A=RNU*WRNU2
      A2=A*A
      A3=A2*A
      DRFK= MAX(1.0D-20,1.D0+0.6D0*A2-c2*A3)**c0
      RFK=RNU*DRFK
      RFK=RFK - RPS
      DADRN=WRNU2*A
      DRFK=DRFK+RNU/(DRFK*DRFK)*(1.2D0*DADRN-c3*A*DADRN)
      IF (RNU.GT.1.D30) RNU=1.D30
      IF (RFK.GT.1.D30) RFK=1.D30
      IF (DRFK.GT.1.D30) DRFK=1.D30
      RETURN
      END
***************************************************************************
      FUNCTION RTSAFE(FUNCD,X1,X2,RACC,IGRID)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (MAXIT=300)
      common/errlim/nrerr
      CALL FUNCD(X1,FL,DF)
      CALL FUNCD(X2,FH,DF)
      XACC=RACC*(ABS(x1)+ABS(x2))
      IF(FL.LT.0.D0)THEN
        XL=X1
        XH=X2
      ELSE
        XH=X1
        XL=X2
        SWAP=FL
        FL=FH
        FH=SWAP
      ENDIF
      RTSAFE=.5D0*(X1+X2)
      DXOLD=ABS(X2-X1)
      DX=DXOLD
      CALL FUNCD(RTSAFE,F,DF)
      DO J=1,MAXIT
         IF(((RTSAFE-XH)*DF-F)*((RTSAFE-XL)*DF-F).GE.0.D0
     &        .OR. ABS(2.D0*F).GT.ABS(DXOLD*DF) ) THEN
            DXOLD=DX
            DX=0.5D0*(XH-XL)
            RTSAFE=XL+DX
            IF(XL.EQ.RTSAFE)RETURN
         ELSE
            DXOLD=DX
            DX=F/DF
            TEMP=RTSAFE
            RTSAFE=RTSAFE-DX
            IF(TEMP.EQ.RTSAFE)RETURN
         ENDIF
         IF(ABS(DX).LT.XACC) RETURN
         CALL FUNCD(RTSAFE,F,DF)
         IF(F.LT.0.D0) THEN
            XL=RTSAFE
            FL=F
         ELSE
            XH=RTSAFE
            FH=F
         ENDIF
      ENDDO
      RETURN
      END

c***********************************************************************
c     FUNCTION PSIINT
c  SUMMIERT DAS IN DEM POTENTIAL ENTHALTENE INTEGRAL STUETZSTELLE FUER
c  STUETZSTELLE AUF
c***********************************************************************

      FUNCTION psiint( ROIP,ROIM, XMIP, AWIP, ETAIP,RNIM,RNIP)

      IMPLICIT REAL*8 (A-H,O-Z)

      RNIPV=RNIP
      RNIMV=RNIM
      PSIint=0.5D0*(ROIP+ROIM)*(RNIPV**8-RNIMV**8)/xmip
      PSIint=PSIint*(5.d0+ETAIP)/(2.D0+ETAIP)*0.125D0*AWIP*AWIP

      END

c***********************************************************************
c    FUNCTION GPSI
c    BERECHNET G ALS FUNKTION VON THETA
c    DAZU MUSS ZUNAECHST R(THETA) BESTIMMT WERDEN
c***********************************************************************

      FUNCTION GPSI(RNU, TEGRA, KR)
      USE CONSTANTS
      IMPLICIT REAL*8(A-H, O-Z)
      INTEGER, PARAMETER :: INTMAX = 101
      COMMON /RADI/ RAE(INTMAX),ROV,ETANU,WR,RPS,XMPS,A
      COMMON/winkel/theta(INTMAX),cost(INTMAX),sint(INTMAX),
     &     cost2(INTMAX),sint2(INTMAX),UTHETA(INTMAX)
      
      R= RNU * (1.D0 - A*UTHETA(KR))
      RI=1.0D0/R
      RI2=RI*RI
      WR2=WR*WR

      RAE(KR)=R
      DPDR=-CG*XMPS*RI2 +CPI4*RI2*RI2*UTHETA(KR)*TEGRA +WR2*R*SINT2(kr)
      DPDTHE= (CPI4*RI2*RI*TEGRA + (WR2*R*R))*COST(kr)*SINT(kr)
      GPSI=SQRT(DPDR*DPDR+(DPDTHE*DPDTHE*RI2))

      END
