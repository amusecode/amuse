      SUBROUTINE HICIRC(RP,ES0,I1,I2,BODYI,TG,TC,ECC1,EDOT,W)
*
*
*       Eccentricity for given circularization time.
*       --------------------------------------------
*
*       Theory of Rosemary Mardling, Ap. J. XX, YYY, 1995.
*       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*
      INCLUDE 'common6.h'
      REAL*8  WW(3),QQ(3),W(2),Q(2),AT0(2),M21,WG(2),QG(2),
     &        BODYI(2),BOD(2)
      REAL*8  A(2),B(2),C(6)
      DATA  WW  /2.119,3.113,8.175/
      DATA  QQ  /0.4909,0.4219,0.2372/
      DATA  A  /6.306505,-7.297806/
      DATA  B  /32.17211,13.01598/
      DATA  C  /5.101417,24.71539,-9.627739,1.733964,
     &                            -2.314374,-4.127795/
*
*
*       Specify index J1 as biggest radius to be used with AT0(1).
      IF (RADIUS(I1).GE.RADIUS(I2)) THEN
          J1 = I1
          J2 = I2
          BOD(1) = BODYI(1)
          BOD(2) = BODYI(2)
      ELSE
          J1 = I2
          J2 = I1
          BOD(1) = BODYI(2)
          BOD(2) = BODYI(1)
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
          IF (KSTAR(IK).EQ.3.OR.KSTAR(IK).EQ.5.OR.
     &        KSTAR(IK).EQ.6.OR.KSTAR(IK).EQ.9) THEN
              BODI = BOD(K)
              CALL GIANT3(IK,BODI,WG,QG,ZN,QL)
              W(K) = WG(1)
              Q(K) = QG(1)
              CONTINUE
          ELSE
              QL = 1.0D+04
              IP = 3
              IF (KSTAR(IK).GE.3) IP = 2
              IF (KSTAR(IK).EQ.4.OR.KSTAR(IK).EQ.7) IP = 3
              IF (KSTAR(IK).EQ.8) IP = 3
              IF (KSTAR(IK).EQ.0) IP = 1
              W(K) = WW(IP)
              Q(K) = QQ(IP)
          END IF
          TL = TWOPI*RADIUS(IK)*SQRT(RADIUS(IK)/BOD(K)/W(K))
          AT0(K) = 1.0/(QL*TL)
    5 CONTINUE
*
*       Form mass, radius & pericentre ratio (note BODYI(2) is a ghost).
      IF (RADIUS(I1).GE.RADIUS(I2)) THEN
          M21 = BODYI(2)/BODYI(1)
          R21 = RADIUS(I2)/RADIUS(I1)
	  RP1 = RP/RADIUS(I1)
      ELSE
	  M21 = BODYI(1)/BODYI(2)
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
      Z = TG*CONST/TSTAR + FF
      ECC1 = (-1.0 + C(1)*Z - SQRT(C(2)*Z**2 + C(3)*Z + C(4)))
     &     /(C(5) + C(6)*Z)
*
*       Evaluate circularization time (in units of 10**6 yrs).
      TC = TSTAR*(1.0 - FF)/CONST
*
*       Obtain (de/dt) due to tidal circularization.
      FE = 1.0 + 3.75*ES0**2 + 1.875*ES0**4 + (5.0/64.0)*ES0**6
      FE = (9.0*TWOPI/10.0)*ES0*(1.0 - ES0**2)**1.5*FE
      EDOT = -CONST*FE
*
      RETURN
*
      END
