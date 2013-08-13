      SUBROUTINE CHPERT(GAMX)
*
*
*       Final chain perturber selection.
*       --------------------------------
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
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      REAL*8  ACC(NMX3),FP(NMX3)
*
*
*       Skip on zero perturbers.
      GAMX = 0.0
      IF (LISTC(1).EQ.0) GO TO 40
*
*       Determine global values of XC first time (also called in REINIT).
      IF (TIMEC.EQ.0.0D0) THEN
          CALL XCPRED(0)
      END IF
*
*       Initialize the external perturbations.
      NK = 3*NN
      DO 1 K = 1,NK
          ACC(K) = 0.0D0
    1 CONTINUE
*
*       Copy provisional perturber list for evaluation loop.
      NPC = LISTC(1) + 1
      DO 2 L = 2,NPC
          JPERT(L) = LISTC(L)
    2 CONTINUE
*
*       Consider each provisional perturber in turn.
      NP = 1
      DO 20 L = 2,NPC
          J = JPERT(L)
          IK = -3
          ITIME = 0
*       Sum perturber contributions over each chain component.
          DO 10 I = 1,NN
              IK = IK + 3
              A1 = X(1,J) - XC(1,I)
              A2 = X(2,J) - XC(2,I)
              A3 = X(3,J) - XC(3,I)
              RIJ2 = A1*A1 + A2*A2 + A3*A3
              A6 = BODY(J)/(RIJ2*SQRT(RIJ2))
              FP(IK+1) = A1*A6
              FP(IK+2) = A2*A6
              FP(IK+3) = A3*A6
              ACC(IK+1) = ACC(IK+1) + FP(IK+1)
              ACC(IK+2) = ACC(IK+2) + FP(IK+2)
              ACC(IK+3) = ACC(IK+3) + FP(IK+3)
*
*       Use standard perturbation test (GMIN) for acceptance.
              IF (I.GT.1) THEN
                  DF2 = 0.0
                  DO 5 K = 1,3
                      DF = FP(IK+K) - FP(IK+K-3)
                      DF2 = DF2 + DF**2
    5             CONTINUE
                  GAM = SQRT(DF2)/((BODYC(I-1) + BODYC(I))*RINV(I-1)**2)
*       Add accepted perturber to LISTC first time only.
                  IF (GAM.GT.GMIN.AND.ITIME.EQ.0) THEN
                      ITIME = 1
                      NP = NP + 1
                      LISTC(NP) = J
                  END IF
              END IF
   10     CONTINUE
   20 CONTINUE
*
*       Save new membership.
      LISTC(1) = NP - 1
*
*       Evaluate the total relative perturbations and save maximum.
      DO 30 I = 1,NN-1
          LI = 3*(I - 1)
          DF2 = 0.0
          DO 25 K = 1,3
              DF = ACC(K+LI+3) - ACC(K+LI)
              DF2 = DF2 + DF**2
   25     CONTINUE
*       Note that rejected perturbers are included in each final value.
          GAM = SQRT(DF2)/((BODYC(I) + BODYC(I+1))*RINV(I)**2)
          GAMX = MAX(GAM,GAMX)
   30 CONTINUE
*
   40 RETURN
      END
