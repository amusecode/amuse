      subroutine rkint(dt,u)

*	Runge-Kutta integrator.
*       -----------------------

*       Author:  Rosemary Mardling (3/98).

      parameter (n=6)
      implicit real*8 (A-H,O-Z)
      real*8 u0(n),ut(n),du(n),u(n)
      real*8 a(4),b(4)

      iq = 0
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
         call deriv(u,du,iq)

         do i=1,n
            u(i)=u0(i)+a(j)*du(i)
            ut(i)=ut(i)+b(j)*du(i)
         enddo
      enddo

      if(iq.gt.0)then
         do i=1,n
            ut(i)=0.0d0
         enddo
      endif

      do i=1,n
         u(i)=u0(i)+ut(i)
      enddo

      end


      subroutine deriv(u,du,iq)

*	Differential equations for hierarchical binary.
*       ----------------------------------------------

*       Author:  Rosemary Mardling (3/98).

      implicit real*8 (A-H,M,O-Z)
      common/rksave/  coeff,HOhat(3),e,a,hh,mb
      common/tidal/  cq(2),ct(2),cgr,dedt
      real*8  ehat(3),hhat(3),qhat(3),edot(3),hdot(3),u(6),du(6)
      SAVE ITIME,IGR
      DATA ITIME,IGR /0,0/
	

*       Save initial basic elements.
      e0 = e
      a0 = a
      h0 = hh
 
*       Update eccentricity and angular momentum.
      e=sqrt(u(1)**2+u(2)**2+u(3)**2)
      if(e.ge.1.0)then
         iq = 1
         goto 10
      endif
      hh=sqrt(u(4)**2+u(5)**2+u(6)**2)
      a=hh**2/mb/(1.0-e**2)

*       Define unit vectors.
      do i=1,3
         ehat(i)=u(i)/e
         hhat(i)=u(i+3)/hh
      enddo

*	Calculate unit vector orthogonal to ehat and hhat (Peter's q).

      call cross(hhat,ehat,qhat)

*	Calculate components of Peter's Sij tensor (third body average).

      S11=-coeff*(3*(dot(HOhat,ehat))**2 - 1)
      S12=-3*coeff*dot(HOhat,ehat)*dot(HOhat,qhat)
      S13=-3*coeff*dot(HOhat,ehat)*dot(HOhat,hhat)
      S22=-coeff*(3*(dot(HOhat,qhat))**2 - 1)
      S23=-3*coeff*dot(HOhat,qhat)*dot(HOhat,hhat)

*	Calculate rate of change of evec and hvec.

      do i=1,3
         edot(i)=(e*a*hh/2/mb)*
     &           (-5*S12*ehat(i)+(4*S11-S22)*qhat(i)-S23*hhat(i))

         hdot(i)=(a**2/2)*((1-e**2)*S23*ehat(i)
     &          -(1+4*e**2)*S13*qhat(i)+5*e**2*S12*hhat(i))
      enddo

      EKD = edot(1)*ehat(1)+edot(2)*ehat(2)+edot(3)*ehat(3)
      TAV = e/sqrt(edot(1)**2+edot(2)**2+edot(3)**2)
*
*       Adopt eccentricity functions from Hut theory (A & A 99, 126, 1981).
      ge=30*e+45*e**3+3.75*e**5
      ge=0.5*ge
      he=(1.0+3.75*e**2+1.875*e**4+5.0/64.0*e**6)/(1.0-e**2)**6.5

*       Specify scaling factors for quadrupole and tidal interaction.
      slr=a*(1.0-e**2)
      zq=ge/(slr**5*a*sqrt(a))
      zz=he/a**8
      zt=9.0*e*zz*(ct(1) + ct(2))
      zq=zq*(cq(1) + cq(2))
*       Correct for wrong eccentricity dependence.
      zq=zq*(1.0 - e**2)
*       Note that angular momentum derivative is of same form as de/dt.
*     zh=hh*zz*(ct(1) + ct(2))
*
*       Form scaled expression for GR precession.
      zg = cgr/(slr*a*sqrt(a))

*     zq = 0.0
*     zg = 0.0
*     zt = 0.0
*     zh = 0.0
*     IF (e.GT.0.995) THEN
*     WRITE (6,9)  (edot(i),i=1,3),(zt*ehat(i),i=1,3)
*   9 FORMAT (' edot zt*ehat   ',1P,3E10.2,2X,3E10.2)
*     END IF
*
*       Include tidal and quadrupole terms for each star and also GR.
      do i=1,3
         edot(i) = edot(i) + zt*ehat(i) + (zq + zg)*qhat(i)
*        hdot(i) = hdot(i) + zh*hhat(i)
      enddo
*       Save scalar value of eccentricity derivative for decision-making.
      dedt = edot(1)*ehat(1)+edot(2)*ehat(2)+edot(3)*ehat(3)
*
*     IF (e.GT.0.995) THEN
*     EDF = edot(1)*ehat(1)+edot(2)*ehat(2)+edot(3)*ehat(3)
*     WRITE (6,3)  e, EKD, (EKD - EDF)/EKD, e/zt, e/zq, e/zg
*   3 FORMAT (' E EKD DED/ED TT TQ TGR  ',F9.5,1P,5E9.1)
*     END IF
*     IF (IGR.LT.10.AND.zg.GT.zq*(cq(1) + cq(2)).AND.e.GT.0.99) THEN
*     IF (e.GT.0.995) THEN
*     ZD = ABS(EKD)
*     ZD = MAX(ZD,zq*(cq(1)+cq(2)))
*     if (e.GT.0.05.AND.zg.GT.ZD) THEN
*         WRITE (6,4)  e, a, a*(1.0-e), zg, ZD, 6.28/zg
*   4     FORMAT ('  GR    E A QP ZGR QUAD P_gr ',
*    &                                F8.4,1P,5E10.2)
*         IGR = IGR + 1
*     END IF
*
*     TAU = e/sqrt(edot(1)**2+edot(2)**2+edot(3)**2)
*     WRITE (6,5)  e, TAV, TAU, edot, zt, zq
*   5 FORMAT (' DERIV    e TAV TAU ed zt zt  ',F9.5,2F9.3,1P,5E10.2)
      ITIME = ITIME + 1
      IF (e.GT.0.980.AND.ITIME.EQ.-4) THEN
      EDT = zt
      EDQ = zq
      TF = e/EDT
      TQ = e/EDQ
      EDTOT = edot(1)*ehat(1)+edot(2)*ehat(2)+edot(3)*ehat(3)
      IF (ITIME.EQ.4) WRITE (6,6)  e0,TF,TQ,EDT,EDQ,EDTOT,EKD
    6 FORMAT (' DERIV    e TF TQ EDT EDQ ED EK ',F9.5,1P,6E10.2)
      END IF
      IF (ITIME.GE.4) ITIME = 0

*       Copy derivatives for Runge-Kutta integrator.
      do i=1,3
         du(i)=edot(i)
         du(i+3)=hdot(i)
      enddo
 10   return
      end
