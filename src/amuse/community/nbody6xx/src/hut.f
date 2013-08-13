      subroutine hut(es0,spin10,spin20,ecc,spin1,spin2,nsteps,dtau)

      implicit real*8 (A-H,O-Z)
      real*8 u(3),udot(3)

      u(1)=es0
      u(2)=spin10
      u(3)=spin20

      call deriv2(u,udot)

*     Include step reduction for large de/dt or primary spin rate.
      IT = 0
 1    IF (ABS(udot(1))*dtau.GT.0.01*MAX(es0,0.01d0)) THEN
         dtau = 0.5D0*dtau
         nsteps = 2*nsteps
         IT = IT + 1
         IF (IT.LT.10) GO TO 1
      END IF

    2 IF (ABS(udot(2))*dtau.GT.0.01*u(2).AND.u(2).GT.0.1) THEN
         dtau = 0.5D0*dtau
         nsteps = 2*nsteps
         IT = IT + 1
         IF (IT.LT.10) GO TO 2
      END IF

*     Treat large eccentricity carefully for rapid spin change.
      IF (es0.GT.0.9995.and.ABS(udot(2))*dtau.GT.
     &     0.01*u(2).AND.u(2).GT.0.1) THEN
         dtau = 0.5D0*dtau
         nsteps = 2*nsteps
         IT = IT + 1
      END IF

      IF (IT.GT.0) THEN
         WRITE (96,3)  nsteps, IT, u, udot, dtau
 3       FORMAT (' HUT REDUCE    # IT u ud dt ',I6,I3,F8.4,1P,7E9.1)
         CALL FLUSH(96)
      END IF

      do i=1,nsteps
*     Save spins in case eccentricity goes negative.
         usave2 = u(2)
         usave3 = u(3)
         call rk4b(dtau,u)
*     Note there are 4 calls to deriv2 when u(1) may go negative.
         IF (u(1).lt.0.002.or.u(1).gt.0.99999) then
*     Enforce circularization by adopting e=0.00199 and copying spins.
            IF (u(1).LT.0.0) WRITE (6,4) I, nsteps, u(1)
 4          FORMAT (' HUT SAFETY EXIT!   I nsteps u1 ',2I5,F8.4)
            u(1) = 0.00199
            u(2) = usave2
            u(3) = usave3
            GO TO 5
         END IF
      enddo

 5    ecc=u(1)
      spin1=u(2)
      spin2=u(3)

      end

      subroutine rk4b(dt,u)

*     Runge-Kutta integrator.
*     -----------------------

*     Author:  Rosemary Mardling (3/98).

      parameter (n=3)
      implicit real*8 (A-H,O-Z)
      real*8 u0(n),ut(n),du(n),u(n)
      real*8 a(4),b(4)

      a(1)=dt/2.0d0
      a(2)=a(1)
      a(3)=dt
      a(4)=0.0d0

      b(1)=dt/6.0d0
      b(2)=dt/3.0d0
      b(3)=b(2)
      b(4)=b(1)

      do i=1,n
         u0(i)=u(i)
         ut(i)=0.0d0
      enddo

      do j=1,4
         call deriv2(u,du)

         do i=1,n
            u(i)=u0(i)+a(j)*du(i)
            ut(i)=ut(i)+b(j)*du(i)
         enddo
      enddo

      do i=1,n
         u(i)=u0(i)+ut(i)
      enddo

      end
      
      subroutine deriv2(u,udot)

      implicit real*8 (A-H,M,O-Z)
      real*8 u(3),udot(3)
      common/spins/angmom0,rg2(2),m21,r21,semi0,C1,C2,C3,C4,semi
*     SAVE IC
*     DATA IC /0/

      e=u(1)
      spin1=u(2)
      spin2=u(3)
      semi1=semi

      e2=e**2
      e4=e2**2
      e6=e4*e2
      fac=1-e2

      f2=1+7.5*e2+5.625*e4+0.3125*e6
      f3=1+3.75*e2+1.875*e4+0.078125*e6
      f4=1+1.5*e2+0.125*e4
      f5=1+3*e2+0.375*e4

      semi=angmom0-rg2(1)*spin1-m21*r21**2*rg2(2)*spin2
      semi=(semi*(1+m21)/m21/semi0**2)**2/fac
      oa = 1.0/semi

      if(e.le.0.0.or.e.ge.1.0.or.oa.lt.0.0)then
         udot(1) = 0.0
         udot(2) = 0.0
         udot(3) = 0.0
      else
         udot(1)=-oa**8*(e/fac**6.5)*
     &        (C1*(f3-(11./18.)*fac**1.5*f4*spin1/oa**1.5)+
     &        C2*(f3-(11./18.)*fac**1.5*f4*spin2/oa**1.5))
         udot(2)=(oa/fac)**6*C3*(oa**1.5*f2-fac**1.5*f5*spin1)
         udot(3)=(oa/fac)**6*C4*(oa**1.5*f2-fac**1.5*f5*spin2)
      endif
      if (e.lt.0.002.and.semi.lt.0.0d0) semi=semi1

*     IC = IC + 1
*     IF (IC.EQ.1.OR.MOD(IC,10000).EQ.0) THEN
*     WRITE (6,1)  IC, e, spin1, spin2, udot
*     1       FORMAT (' HUT DERIV    # e s1 s2 udot ',I10,F8.4,1P,5E9.1)
*     END IF

      end
