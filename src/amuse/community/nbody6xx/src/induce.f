      SUBROUTINE INDUCE(IPAIR,EMAX,EMIN,ICIRC,TC2,ANGLE,TG,EDAV)
*
*
*       Induced eccentricity of hierarchical binary.
*       --------------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  A1(3),A2(3),XREL(3),VREL(3),EI(3),HI(3),HO(3),BHAT(3)
*
*
*       Define c.m. & KS indices and evaluate semi-major axis & eccentricity.
      I = N + IPAIR
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      SEMI = -0.5d0*BODY(I)/H(IPAIR)
      ECC2 = (1.d0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2
     &                                   /(BODY(I)*SEMI)
      ECC = SQRT(ECC2)
*
*       Ensure that KS components are resolved (perturbed not always done).
      CALL RESOLV(IPAIR,1)
*
*       Obtain inner & outer angular momentum by cyclic notation.
      RIJ2 = 0.d0
      VIJ2 = 0.d0
      RDOT = 0.d0
      A12 = 0.d0
      A22 = 0.d0
      A1A2 = 0.d0
      RI2 = 0.d0
      VI2 = 0.d0
      RVI = 0.d0
      DO 5 K = 1,3
          RIJ2 = RIJ2 + (X(K,I) - X(K,JCOMP))**2
          VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,JCOMP))**2
          RDOT = RDOT + (X(K,I) - X(K,JCOMP))*(XDOT(K,I) -XDOT(K,JCOMP))
          K1 = K + 1
          IF (K1.GT.3) K1 = 1
          K2 = K1 + 1
          IF (K2.GT.3) K2 = 1
          A1(K) = (X(K1,I1) - X(K1,I2))*(XDOT(K2,I1) - XDOT(K2,I2))
     &          - (X(K2,I1) - X(K2,I2))*(XDOT(K1,I1) - XDOT(K1,I2))
          A2(K) = (X(K1,JCOMP) - X(K1,I))*(XDOT(K2,JCOMP) - XDOT(K2,I))
     &          - (X(K2,JCOMP) - X(K2,I))*(XDOT(K1,JCOMP) - XDOT(K1,I))
          A12 = A12 + A1(K)**2
          A22 = A22 + A2(K)**2
          A1A2 = A1A2 + A1(K)*A2(K)
*       Form relative vectors and scalars for inner binary.
          XREL(K) = X(K,I1) - X(K,I2)
          VREL(K) = XDOT(K,I1) - XDOT(K,I2)
          RI2 = RI2 + XREL(K)**2
          VI2 = VI2 + VREL(K)**2
          RVI = RVI + XREL(K)*VREL(K)
    5 CONTINUE
*
*       Evaluate orbital parameters for outer orbit.
      RIJ = SQRT(RIJ2)
      ZMB = BODY(I) + BODY(JCOMP)
      SEMI1 = 2.d0/RIJ - VIJ2/ZMB
      SEMI1 = 1.d0/SEMI1
      ECC1 = SQRT((1.d0 - RIJ/SEMI1)**2 + RDOT**2/(SEMI1*ZMB))
*     PMIN1 = SEMI1*(1.d0 - ECC1)
*       Determine inclination (8 bins of 22.5 degrees).
      FAC = A1A2/SQRT(A12*A22)
      FAC = ACOS(FAC)
      ANGLE = FAC
***   ANGLE = FAC*360.d0/TWOPI
      IN = 1 + FAC*360.0/(TWOPI*22.5)
*
*       Exit on hyperbolic orbit or outer eccentricity near 1.
      IF (SEMI1.LT.0.0.OR.ECC1.GT.0.9999) THEN
          EMAX = 0.d0
          TC2 = 999.d0
          TG = 1.0d+04
          EDAV = 1.0d+04
          GO TO 30
      END IF
*
*       Construct the Runge-Lenz vector (Heggie & Rasio 1995, Eq.(5)).
      EI2 = 0.d0
      DO 10 K = 1,3
          EI(K) = (VI2*XREL(K) - RVI*VREL(K))/BODY(I) -
     &                                                 XREL(K)/SQRT(RI2)
          EI2 = EI2 + EI(K)**2
   10 CONTINUE
*
*       Define unit vectors for inner eccentricity and angular momenta.
      COSJ = 0.d0
      SJSG = 0.d0
      DO 15 K = 1,3
          EI(K) = EI(K)/SQRT(EI2)
          HI(K) = A1(K)/SQRT(A12)
          HO(K) = A2(K)/SQRT(A22)
          COSJ = COSJ + HI(K)*HO(K)
          SJSG = SJSG + EI(K)*HO(K)
   15 CONTINUE
*
*       Form unit vector BHAT and scalars AH & BH (Douglas Heggie, 10/9/96).
      AH = 0.d0
      BH = 0.d0
      DO 20 K = 1,3
          K1 = K + 1
          IF (K1.GT.3) K1 = 1
          K2 = K1 + 1
          IF (K2.GT.3) K2 = 1
          BHAT(K) = HI(K1)*EI(K2) - HI(K2)*EI(K1)
          AH = AH + EI(K)*HO(K)
          BH = BH + BHAT(K)*HO(K)
   20 CONTINUE
*
*       Evaluate the expressions A & Z.
      A = COSJ*SQRT(1.d0 - EI2)
      Z = (1.d0 - EI2)*(2.d0 - COSJ**2) + 5.d0*EI2*SJSG**2
*
*       Obtain maximum inner eccentricity (Douglas Heggie, Sept. 1995).
      Z2 = Z**2 + 25.d0 + 16.d0*A**4 - 10.d0*Z - 
     &     20.d0*A**2 - 8.d0*A**2*Z
      EMAX = ONE6*(Z + 1.d0 - 4.d0*A**2 + SQRT(Z2))
      EMAX = MAX(EMAX,0.0001d0)
      EMAX = SQRT(EMAX)
      ZI = FAC
      EM = SQRT(SIN(ZI)**2 + EI2*COS(ZI)**2)
*
*       Form minimum eccentricity (Douglas Heggie, Sept. 1996).
      AZ = A**2 + Z - 2.d0
      IF (AZ.GE.0.0) THEN
          AZ1 = 1.d0 + Z - 4.d0*A**2
          EMIN2 = ONE6*(AZ1 - SQRT(AZ1**2 - 12.d0*AZ))
      ELSE
          EMIN2 = 1.d0 - 0.5d0*(A**2 + Z)
      END IF
      EMIN2 = MAX(EMIN2,0.0001d0)
      EMIN = SQRT(EMIN2)
*
*       Set impact parameters and estimate eccentricity growth time-scale.
      PMIN = SEMI*(1.d0 - ECC)
      PMIN2 = SEMI*(1.d0 - EMAX)
      TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
      TK1 = TWOPI*SEMI1*SQRT(SEMI1/ZMB)
      TG = TK1**2*ZMB*(1.d0 - ECC1**2)**1.5/(BODY(JCOMP)*TK)
*
      CONST = PFAC(A,Z)
      CONST = CONST*4.d0/(1.5d0*TWOPI*SQRT(6.d0))
*       Convert growth time to units of 10**6 yrs.
      TG = CONST*TG*TSTAR
*
*       Form doubly averaged eccentricity derivative (Douglas Heggie 9/96).
      YFAC = 15.d0*BODY(JCOMP)/(4.d0*ZMB)*TWOPI*TK/TK1**2
      YFAC = YFAC*ECC*SQRT(1.d0 - ECC**2)/
     &                             (1.d0 - ECC1**2)**(1.5)
      EDAV = YFAC*AH*BH
*
*       Determine circularization time for current & smallest pericentre.
      TC2 = 2000.d0
      TC = 2000.d0
      IF (PMIN2.LT.10.0*MAX(RADIUS(I1),RADIUS(I2)).AND.
     &    KZ(27).EQ.2.AND.ECC.GT.0.002) THEN
          CALL TCIRC(PMIN,ECC,I1,I2,ICIRC,TC)
          IF (ICIRC.GE.0) ICIRC = -1
*       Note: only use PMIN call because ICIRC may become >0 after first!
          CALL TCIRC(PMIN2,EMAX,I1,I2,ICIRC,TC2)
      ELSE IF (ECC.LT.0.002) THEN
          TC = 0.0
          TC2 = 0.0
      END IF
*
*       Print diagnostics for high inclinations and TC2 < 10**7 yrs.
      IF (IN.GE.4.AND.IN.LE.5.AND.TC2.LT.10.0.AND.KSTAR(I).NE.-2.AND.
     &    EI2.GT.4.0D-06) THEN
          WRITE (44,25)  SQRT(EI2), EM, EMAX, KSTAR(I1), KSTAR(I2),
     &                   KSTAR(I), SEMI, PMIN2, IN, TG, TC, TC2, TPHYS
   25     FORMAT (' INDUCE:    E EMX EMAX K* SEMI PM2 IN TG TC TC2 TP ',
     &                         3F8.4,3I3,1P,2E9.1,0P,I3,F7.3,3F7.1)
          CALL FLUSH(44)
      END IF
*
*       Use first value for now.
      TC2 = TC
*
   30 RETURN
*
      END
