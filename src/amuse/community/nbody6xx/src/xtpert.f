      SUBROUTINE XTPERT(ACC,CHTIME)
*
*
*       External perturbations on chain members.
*       ----------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  M,MASS,MC,MIJ,MKK
      PARAMETER  (NMX=10,NMX3=3*NMX,NMXm=NMX*(NMX-1)/2)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XJ(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      REAL*8  ACC(NMX3)
*
*
*       Predict current coordinates of perturbers & chain members.
      TIME0 = TIME
      ISUB = ISYS(5)
      TIME = CHTIME + T0S(ISUB)
      CALL XCPRED(1)
*
*       Initialize the external perturbations.
      NK = 3*NN
      DO 1 K = 1,NK
          ACC(K) = 0.0D0
    1 CONTINUE
*
*       Consider each perturber in turn.
      NPC = LISTC(1) + 1
      DO 20 L = 2,NPC
          J = LISTC(L)
*       Sum perturber contributions over each chain component.
          IK = -3
          KS = 0
          DO 10 I = 1,NN
              IK = IK + 3
              A1 = X(1,J) - XC(1,I)
              A2 = X(2,J) - XC(2,I)
              A3 = X(3,J) - XC(3,I)
              RIJ2 = A1*A1 + A2*A2 + A3*A3
              IF (J.LE.N) GO TO 5
              JPAIR = J - N
*       Check c.m. approximation (only resolve KS components once).
              IF (RIJ2.GT.CMSEP2*R(JPAIR)**2) GO TO 5
              IF (KS.EQ.0) THEN
                  CALL KSRES(JPAIR,J1,J2,RIJ2)
                  KS = 1
              END IF
              KDUM = J1
              K = KDUM
    4         A1 = X(1,K) - XC(1,I)
              A2 = X(2,K) - XC(2,I)
              A3 = X(3,K) - XC(3,I)
              RIJ2 = A1*A1 + A2*A2 + A3*A3
              A6 = BODY(K)/(RIJ2*SQRT(RIJ2))
              ACC(IK+1) = ACC(IK+1) + A1*A6
              ACC(IK+2) = ACC(IK+2) + A2*A6
              ACC(IK+3) = ACC(IK+3) + A3*A6
*       See whether the second component has been included.
              IF (K.EQ.KDUM) THEN
                  K = K + 1
                  GO TO 4
              END IF
              GO TO 10
*
    5         A6 = BODY(J)/(RIJ2*SQRT(RIJ2))
              ACC(IK+1) = ACC(IK+1) + A1*A6
              ACC(IK+2) = ACC(IK+2) + A2*A6
              ACC(IK+3) = ACC(IK+3) + A3*A6
   10     CONTINUE
   20 CONTINUE
*
*       Restore the global time.
      TIME = TIME0
*
      RETURN
      END
