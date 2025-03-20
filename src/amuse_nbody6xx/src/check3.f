      SUBROUTINE CHECK3
*
*
*       Three-body configuration check.
*       -------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Consider possible unstable triple configurations.
      V2 = 0.5*ZMASS/RSCALE
      DO 50 IPAIR = 1,NPAIRS
          I = N + IPAIR
          IF (NAME(I).LT.0.OR.BODY(I).EQ.0.0D0) GO TO 50
          I1 = 2*IPAIR - 1
          NP = LIST(1,I1)
          IF (NP.EQ.0.OR.NP.GE.2) GO TO 50
          J = LIST(2,I1)
          ZR = BODY(J)/BODY(I)
          IF (NAME(J).LE.0.OR.ZR.GT.1.0) GO TO 50
          RIJ2 = 0.0
          VIJ2 = 0.0
          RD = 0.0
          DO 10 K = 1,3
              RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
              VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,J))**2
              RD = RD + (X(K,I)-X(K,J))*(XDOT(K,I)-XDOT(K,J))
   10     CONTINUE
          RIJ = SQRT(RIJ2)
          ZM = BODY(I) + BODY(J)
          SEMI = 2.0/RIJ - VIJ2/ZM
          SEMI = 1.0/SEMI
          SEMI0 = -0.5*BODY(I)/H(IPAIR)
          ECC2 = (1.0 - R(IPAIR)/SEMI0)**2 +
     &            TDOT2(IPAIR)**2/(SEMI0*BODY(I))
          ECC = SQRT(ECC2)
          E2 = (1.0 - RIJ/SEMI)**2 + RD**2/(SEMI*ZM)
          ECC1 = SQRT(E2)
*       Include small hyperbolic velocity to catch disrupting system.
          HI = -0.5*ZM/SEMI
          IF (HI.LT.0.01*V2) THEN
              WRITE (48,20)  TIME+TOFF, NAME(I), NAME(J), ZR, ECC1,
     &                       SEMI0, SEMI
   20         FORMAT (' CHECK3    T NM Z E A0 A ',
     &                            F9.1,2I7,F5.1,F7.3,1P,4E10.2)
              CALL FLUSH(48)
          END IF
   50 CONTINUE
*
      RETURN
*
      END
