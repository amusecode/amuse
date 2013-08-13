      SUBROUTINE TRANSK
*
*
*       Transformation to physical momenta & separations.
*       -------------------------------------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIJ(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
*
*
*       Obtain chain momenta & separations of configuration QK & PK.
      DO 1 I = 1,N-1
          KS1 = 4*(I - 1) + 1
          L1 = 3*(I - 1) + 1
          CALL KSPHYS(QK(KS1),PK(KS1),XC(L1),WC(L1))
          RIK = QK(KS1)**2 + QK(KS1+1)**2 + QK(KS1+2)**2 + QK(KS1+3)**2
          RINV(I) = 1.0/RIK
    1 CONTINUE
*
*       Obtain physical coordinates.
      DO 2 K = 1,3
          XI(K) = 0.0D0
    2 CONTINUE
      DO 5 I = 1,N-1
          L = 3*(I - 1)
          DO 4 K = 1,3
              XI(L+3+K) = XI(L+K) + XC(L+K)
    4     CONTINUE
    5 CONTINUE
*
*       Form non-singular inverse distances.
      LRI = N - 1
      DO 20 I = 1,N-2
          LI = 3*(I - 1)
          DO 15 J = I+2,N
              LJ = 3*(J - 1)
              RIJ2 = 0.0D0
              IF (J.GT.I + 2) THEN
                  DO 10 K = 1,3
                      RIJ2 = RIJ2 + (XI(LJ+K) - XI(LI+K))**2
   10             CONTINUE
              ELSE
                  DO 12 K = 1,3
                      RIJ2 = RIJ2 + (XC(LI+K) + XC(LI+K+3))**2
   12             CONTINUE
              END IF
              LRI = LRI + 1
              RINV(LRI) = 1.0/SQRT(RIJ2)
   15     CONTINUE
   20 CONTINUE
*
*	Transform to physical variables from chain quantities.
      L = 3*(N - 2)
      DO 25 K = 1,3
          PI(K) = -WC(K)
          PI(L+K+3) = WC(L+K)
   25 CONTINUE
*
      DO 30 I = 2,N-1
          L = 3*(I - 1)
          DO 28 K = 1,3
              PI(L+K) = WC(L+K-3) - WC(L+K)
   28     CONTINUE
   30 CONTINUE
*
      RETURN
*
      END
