      subroutine hut2(spin10,spin20,spin1,spin2,nsteps,dtau)
*     
*     Spin evolution of circular binary.
*     ----------------------------------

      implicit real*8 (A-H,O-Z)
      real*8 u(2),udot(2)
      COMMON/CFLAG/ IFLAG,IQ

      u(1)=spin10
      u(2)=spin20
      IFLAG = 0
      IQ = 0

      call deriv3(u,udot)

*     Include step reduction for large derivatives.
      IT = 0
 1    IF (ABS(udot(1))*dtau.GT.0.01*spin10) THEN
         dtau = 0.5D0*dtau
         nsteps = 2*nsteps
         IT = IT + 1
         IF (IT.LT.5) GO TO 1
      END IF

 5    IF (ABS(udot(2))*dtau.GT.0.01*spin20) THEN
         dtau = 0.5D0*dtau
         nsteps = 2*nsteps
         IT = IT + 1
         IF (IT.LT.5) GO TO 5
      END IF

*     do i=1,nsteps
*     call rk4c(dtau,u)
*     IF (IFLAG.GT.0) GO TO 10
*     enddo

      ITER = 0
 6    call rk4c(dtau,u)
      ITER = ITER + 1
      IF (IFLAG.EQ.0.AND.ITER.LT.NSTEPS) GO TO 6

      spin1=u(1)
      spin2=u(2)
      end

      subroutine rk4c(dt,u)

*     Runge-Kutta integrator.
*     -----------------------

*     Author:  Rosemary Mardling (3/98).

      parameter (n=2)
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
         call deriv3(u,du)

         do i=1,n
            u(i)=u0(i)+a(j)*du(i)
            ut(i)=ut(i)+b(j)*du(i)
         enddo
      enddo

      do i=1,n
         u(i)=u0(i)+ut(i)
      enddo

      end
      
      subroutine deriv3(u,udot)

      implicit real*8 (A-H,M,O-Z)
      real*8 u(2),udot(2)
      common/spins/angmom0,rg2(2),m21,r21,semi0,C1,C2,C30,C40,C5,semi
      common/radii/  R1,R2
      COMMON/CFLAG/ IFLAG,IQ
      SAVE IC,ID
      DATA IC,ID /0,0/

*     Assume e=0.
      spin1=u(1)
      spin2=u(2)

      C3=C30*R1**6
      C4=C40*R2**6

*     e2=e**2
*     e4=e2**2
*     e6=e4*e2
*     fac=1-e2

*     f2=1+7.5*e2+5.625*e4+0.3125*e6
*     f3=1+3.75*e2+1.875*e4+0.078125*e6
*     f4=1+1.5*e2+0.125*e4
*     f5=1+3*e2+0.375*e4
      

      semi=angmom0-rg2(1)*spin1-m21*r21**2*rg2(2)*spin2
*     semi=(semi*(1+m21)/m21/semi0**2)**2/fac
      semi=(semi*(1+m21)/m21/semi0**2)**2
      oa = 1.0/semi
      IF (IQ.EQ.0) THEN
         oa0 = oa
         IQ = 1
      END IF
      
*     udot(1)=(oa/fac)**6*C3*(oa**1.5*f2-fac**1.5*f5*spin1)
*     udot(2)=(oa/fac)**6*C4*(oa**1.5*f2-fac**1.5*f5*spin2)
      udot(1)=oa**6*C3*(oa**1.5-spin1) - C5*spin1
      udot(2)=oa**6*C4*(oa**1.5-spin2)

*     Quit on relative change exceeding 2 % (SJA 05/09).
      IF (ABS(oa - oa0).GT.0.02*oa0) IFLAG = 1

*     IC = IC + 1
*     IC = IC + 1
*     IF (IC.LT.10) THEN
*     WRITE (6,1)  IC, spin1, spin2, udot, C5
*     1       FORMAT (' HUT DERIV    # s1 s2 udot C5 ',I5,1P,6E9.1)
*     END IF

      end
