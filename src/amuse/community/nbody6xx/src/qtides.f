      SUBROUTINE QTIDES(I1,I2,IM,SEMI,ES0)
*
*
*       Quadrupole and tidal terms for hierarchical binary.
*       ---------------------------------------------------
*
*       Theory of Rosemary Mardling, Ap. J. XX, YYY, 1995.
*       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      common/tidal/  cq(2),ct(2),cgr,dedt
      REAL*8  WW(3),QQ(3),W(2),Q(2),AT0(2),M21,WG(2),QG(2),BOD(2),
     &        A(2),B(2),C(6),QD(2),k1,k2
      DATA  WW  /2.119,3.113,8.175/
      DATA  QQ  /0.4909,0.4219,0.2372/
      DATA  A  /6.306505,-7.297806/
      DATA  B  /32.17211,13.01598/
      DATA  C  /5.101417,24.71539,-9.627739,1.733964,
     &                            -2.314374,-4.127795/
      SAVE ITIME
      DATA ITIME /0/
*
*
*       Specify index J1 as biggest radius to be used with AT0(1).
      IF (RADIUS(I1).GE.RADIUS(I2)) THEN
          J1 = I1
          J2 = I2
          BOD(1) = CM(1,IM)
          BOD(2) = CM(2,IM)
      ELSE
          J1 = I2
          J2 = I1
          BOD(1) = CM(2,IM)
          BOD(2) = CM(1,IM)
      END IF
*
*       Define oscillation period (dimensionless time) and damping constants.
      DO 5 K = 1,2
          IF (K.EQ.1) THEN
              IK = J1
          ELSE
              IK = J2
          END IF
*       Specify polytropic index for each star (n = 3, 2 or 3/2).
          IF (KSTAR(IK).EQ.3.OR.KSTAR(IK).EQ.5) THEN
              BODI = BOD(K)
              CALL GIANT3(IK,BODI,WG,QG,ZN,QL)
              W(K) = WG(1)
              Q(K) = QG(1)
          ELSE
              QL = 1.0D+04
              IP = 3
              IF (KSTAR(IK).GE.3) IP = 2
              IF (KSTAR(IK).EQ.4.OR.KSTAR(IK).EQ.6) IP = 3
              IF (KSTAR(IK).EQ.0) IP = 1
              W(K) = WW(IP)
              Q(K) = QQ(IP)
          END IF
          TL = TWOPI*RADIUS(IK)*SQRT(RADIUS(IK)/BOD(K)/W(K))
          AT0(K) = 1.0/(QL*TL)
          QD(K) = QL
    5 CONTINUE
*
*       Form mass, radius & pericentre ratio.
      RP = SEMI*(1.0 - ES0)
      IF (RADIUS(I1).GE.RADIUS(I2)) THEN
          M21 = CM(2,IM)/CM(1,IM)
          R21 = RADIUS(I2)/RADIUS(I1)
	  RP1 = RP/RADIUS(I1)
      ELSE
	  M21 = CM(1,IM)/CM(2,IM)
          R21 = RADIUS(I1)/RADIUS(I2)
	  RP1 = RP/RADIUS(I2)
      END IF
*
*	Evaluate damping coefficient.
      RR = RP1*(1.0 + ES0)
      CONST = 2.0*(AT0(1)*(Q(1)/W(1))**2*(1.0 + M21)*M21 + 
     &             AT0(2)*(Q(2)/W(2))**2*((1.0 + M21)/M21**2)*R21**8)/
     &                                                         RR**8
*
*	Form rational function approximation to Hut solution.
      FF = (( A(2)*ES0 + A(1))*ES0 + 1.0 )/
     &     (( B(2)*ES0 + B(1))*ES0 + 1.0 )
      FF = MIN(FF,0.999D0)
*
*       Determine eccentricity corresponding to t_{circ} = t_{grow}.
*     Z = TG*CONST/TSTAR + FF
*     ECC1 = (-1.0 + C(1)*Z - SQRT(C(2)*Z**2 + C(3)*Z + C(4)))
*    &	                                      /(C(5) + C(6)*Z)
*
*       Evaluate circularization time (in units of 10**6 yrs).
      TC = TSTAR*(1.0 - FF)/CONST
*
*       Obtain (de/dt) due to tidal circularization.
      FE = 1.0 + 3.75*ES0**2 + 1.875*ES0**4 + (5.0/64.0)*ES0**6
      FE = (9.0*TWOPI/10.0)*ES0*(1.0 - ES0**2)**1.5*FE
      EDOT = -CONST*FE
*
*       Define mass ratio and apsidal motion constants.
      qin = BOD(2)/BOD(1)
      k1 = 0.2*TWOPI*Q(1)**2/W(1)
      k2 = 0.2*TWOPI*Q(2)**2/W(2)
*
*       Absorb constants to simplify derivatives (note: n^2=(1+qin)/a^3).
      CQ(1) = k1*qin*SQRT(BOD(1)*(1.0+qin))*RADIUS(J1)**5
      CQ(2) = k2*SQRT(BOD(1)*(1.0+qin))*RADIUS(J2)**5/qin
*       Note that BOD(1)*(1+qin) is used instead of BODY(I1) for clarity.
      CT(1) = -(1.0/TWOPI)*k1*qin*(1.0+qin)/(QD(1)*SQRT(W(1)))
      CT(2) = -(1.0/TWOPI)*k2*(1.0+qin)/qin**2/(QD(2)*SQRT(W(2)))
      CT(1) = CT(1)*SQRT(BOD(1))*RADIUS(J1)**6.5
      CT(2) = CT(2)*SQRT(BOD(2))*RADIUS(J2)**6.5
*
*       Form gravitational precession constant (Holman et al, Nature 1998).
      TGR = 3.4D+07*YRS*RAU/(BODY(I1)*SMU*SQRT(BODY(I1)))
*       Note: P_{gr} = 3.4E+07*(1-e^2)*P_{yr}*a_{AU}/m_{sun}.
      TGR = TGR/(1.0D+06*TSTAR)
*       Convert to angular velocity from N-body time units.
      CGR = TWOPI/TGR
*
      ITIME = ITIME + 1
      IF (ITIME.LE.1) THEN
      TF = TIME + TOFF + (1.0 - FF)/CONST
      WRITE (6,20)  CQ, CT, EDOT, TF, CGR
   20 FORMAT (' QTIDES    CQ CT EDT TF CGR  ',1P,7E10.2)
      WRITE (6,30)  k1,k2,QD,W
   30 FORMAT (' k1 k2 QD W   ',2F10.5,1P,4E10.2)
      TB = YRS*SEMI*SQRT(SEMI/BODY(I1))
      TGR = TGR*(1.0 - ES0**2)*SEMI**2.5
      TGR = 1.0D+06*TSTAR*TGR
      WRITE (6,40)  TB, TGR, SEMI*RAU, ES0
   40 FORMAT (' TB(yr) TGR A(AU) E  ',1P,3E10.2,0P,F8.4)
      CALL FLUSH(6)
      END IF
*
      RETURN
*
      END
