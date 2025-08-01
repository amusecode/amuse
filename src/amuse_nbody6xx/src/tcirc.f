      SUBROUTINE TCIRC(RP,ES0,I1,I2,ICIRC,TC)
*     
*     
*     Circularization time.
*     ---------------------
*     
*     Theory of Piet Hut, A & A 99, 126 (1981).
*     Developed by Rosemary Mardling (31/1/97).
*     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*     
      INCLUDE 'common6.h'
      REAL*8  WW(3),QQ(3),W(2),Q(2),AT0(2),M21,WG(2),QG(2),
     &     WSCALE(2),QSCALE(2),A(2),B(2),C(6)
      DATA  WW  /2.119,3.113,8.175/
      DATA  QQ  /0.4909,0.4219,0.2372/
      DATA  A  /6.306505,-7.297806/
      DATA  B  /32.17211,13.01598/
      DATA  C  /5.101417,24.71539,-9.627739,1.733964,
     &     -2.314374,-4.127795/
      DATA  ECCM,ECCM1 /0.002,0.002001/
      SAVE ITIME
      DATA ITIME /0/
*     
*     
*     Set large circularization time for merged binary.
      IF (RADIUS(I1).EQ.0.0D0.OR.RADIUS(I2).EQ.0.0D0) THEN
         TC = 1.0D+10
         GO TO 30
      END IF
*     
*     Specify index J1 as biggest radius to be used with AT0(1).
      IF (RADIUS(I1).GE.RADIUS(I2)) THEN
         J1 = I1
         J2 = I2
      ELSE
         J1 = I2
         J2 = I1
      END IF
*     
*     Include an appropriate numerical factor (Chris Tout 07/08).
      AR = RP/((1.0 - ES0)*RADIUS(J1))
      FF = 1.0
      IF (KSTAR(J1).EQ.1.AND.BODY(J1)*SMU.GT.1.25) THEN
         FF = 100.0
      ELSE IF (KSTAR(J1).EQ.4.OR.KSTAR(J1).EQ.7) THEN
         FF = 100.0
      ELSE IF (KSTAR(J1).GE.10) THEN
         FF = 1.0D+17/AR**2
      END IF
      Q12 = BODY(J1)/BODY(J2)
      TAUC = FF*2.0*Q12**2/(1.0 + Q12)*AR**8
*     Exit on TC > 1.0D+08 yr (Chris Tout, Cambody Book 2008).
      IF (TAUC.GT.1.0D+08) THEN
         TC = 1.0D+10
         GO TO 30
      END IF
*     
*     Define oscillation period (dimensionless time) and damping constants.
      ZN = 0.0
      DO 5 K = 1,2
         IF (K.EQ.1) THEN
            IK = J1
         ELSE
            IK = J2
         END IF
*     Specify polytropic index for each star (n = 3, 2 or 3/2).
         IF (KSTAR(IK).EQ.3.OR.KSTAR(IK).EQ.5.OR.
     &        KSTAR(IK).EQ.6.OR.KSTAR(IK).EQ.9) THEN
            IPAIR = KVEC(I1)
            CALL GIANT(IPAIR,IK,WG,QG,WSCALE,QSCALE,ZN,QL)
            W(K) = WG(1)
            Q(K) = QG(1)
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
         TL = TWOPI*RADIUS(IK)*SQRT(RADIUS(IK)/BODY(IK)/W(K))
         AT0(K) = 1.0/(QL*TL)
    5 CONTINUE
*     
*     Form mass, radius & pericentre ratio.
      IF (RADIUS(I1).GE.RADIUS(I2)) THEN
         M21 = BODY(I2)/BODY(I1)
         R21 = RADIUS(I2)/RADIUS(I1)
         RP1 = RP/RADIUS(I1)
      ELSE
         M21 = BODY(I1)/BODY(I2)
         R21 = RADIUS(I1)/RADIUS(I2)
         RP1 = RP/RADIUS(I2)
      END IF
*     
*     Evaluate damping coefficient.
      RR = RP1*(1.0 + ES0)
      CONST = 2.0*(AT0(1)*(Q(1)/W(1))**2*(1.0 + M21)*M21 + 
     &     AT0(2)*(Q(2)/W(2))**2*((1.0 + M21)/M21**2)*R21**8)/
     &     RR**8
*     
*     Adopt WD scaling for any NS to avoid numerical problem.
      IF (KSTAR(I1).EQ.13.OR.KSTAR(I2).EQ.13) THEN
         CONST = 1.0D-04*CONST
      END IF
*     
*     Form rational function approximation to Hut solution.
      FF = (( A(2)*ES0 + A(1))*ES0 + 1.0 )/
     &     (( B(2)*ES0 + B(1))*ES0 + 1.0 )
*     FF = MIN(FF,0.999D0)
*     
*     See whether we only want the modified eccentricity (routine BINPOP).
      IF (ICIRC.LE.0) GO TO 10
*     
      TIME0 = 0.0
*     Obtain the new eccentricity.
      Z = (TIME - TIME0)*CONST + FF
      ECC = (-1.0 + C(1)*Z - SQRT(C(2)*Z**2 + C(3)*Z + C(4)))
     &     /(C(5) + C(6)*Z)
*     
      ECC = MAX(ECC,ECCM)
      RP = RP*(1.0 + ES0)/(1.0 + ECC)
      ES0 = ECC
      GO TO 30
*     
*     Evaluate circularization time (in units of 10**6 yrs).
 10   TC = TSTAR*(1.0 - FF)/CONST
*     
*     Activate tidal indicator if TC < 2x10**9 yrs or hyperbolic orbit.
      IF (TC.LT.2000.0.OR.ES0.GT.1.0) THEN
         IP = KVEC(I1)
         ITIME = ITIME + 1
*     Note possibility of counter exceeding the limit.
         IF (ITIME.GT.2000000000) ITIME = 0
         IF (ICIRC.EQ.0.AND.KZ(27).EQ.2.AND.ITIME.LT.100) THEN
            SEMI = -0.5*BODY(N+IP)/H(IP)
            if(rank.eq.0)then
               WRITE (6,20)  I1, NCHAOS, ES0, RP1, M21, TC, SEMI, ZN
 20            FORMAT (' TCIRC:    I1 NCH E RP M21 TC A n ',
     &              2I5,F8.4,F8.1,F6.2,1P,2E10.2,0P,F5.1)
            end if
         END IF
         ICIRC = 1
*     Define Roche search indicator for circularized orbit (ECCM1 > 0.002).
         IF (ES0.LE.ECCM1.AND.KSTAR(N+IP).EQ.0) THEN
            KSTAR(N+IP) = 10
         END IF
      ELSE
*     Note ICIRC = -1 for some calls.
         ICIRC = 0
      END IF
*     
 30   RETURN
*     
      END
