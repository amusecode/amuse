      SUBROUTINE HIRECT(IM)
*
*
*       Rectification of hierarchical binary.
*       -------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
*
*
*       Include diagnostic error check.
      RB = 0.0
      UPR2 = 0.0
      TD2 = 0.0
      DO 5 K = 1,4
          RB = RB + UM(K,IM)**2
          UPR2 = UPR2 + UMDOT(K,IM)**2
          TD2 = TD2 + 2.0*UM(K,IM)*UMDOT(K,IM)
    5 CONTINUE
*
      ZMB = CM(1,IM) + CM(2,IM)
      HI = (2.0*UPR2 - ZMB)/RB
      ERR = (HI - HM(IM))/HI
      ZMU = CM(1,IM)*CM(2,IM)/ZMB
      DB = ZMU*(HI - HM(IM))
      IF (ABS(DB).GT.1.0D-08) THEN
      SEMI = -0.5*ZMB/HM(IM)
      ECC2 = (1.0 - RB/SEMI)**2 + TD2**2/(ZMB*SEMI)
      ECC = SQRT(ECC2)
      RA = RB/SEMI
      WRITE (16,3) TTOT, NAMEG(IM), KSTARM(IM), ECC, RA, HM(IM), DB, ERR
    3 FORMAT (' HIRECT:   T NM K* E R/A H DB DH/H ',
     &                    F9.3,I6,I4,F8.4,F8.4,F9.2,1P,2E10.1)
      CALL FLUSH(16)
      END IF
*       Initialize iteration counter for difficult case (SJA 10/97).
      ITER = 0
*
*       Form square regularized velocity for the explicit binding energy.
   10 RB = 0.0D0
      UPR2 = 0.0
      DO 15 K = 1,4
          RB = RB + UM(K,IM)**2
          UPR2 = UPR2 + UMDOT(K,IM)**2
   15 CONTINUE
*
*       Form KS scaling factors from energy and angular momentum relation.
      A1 = 0.25D0*ZMB/UPR2
*       Solve for C1 from H = (2*U'*U'*C1**2 - M)/(U*U*C2**2) with C2 = 1/C1.
      A2 = A1**2 + 0.5D0*HM(IM)*RB/UPR2
*
*       Avoid negative round-off value on second try (NB! no change in CK).
      IF (ITER.EQ.2.AND.A2.LT.0.0) A2 = 0.0D0
*
*       Check for undefined case (circular orbit or eccentric anomaly = 90).
      IF (A2.GE.0.0D0) THEN
          IF (A1.LT.1.0) THEN
*       Choose square root sign from eccentric anomaly (e*cos(E) = 1 - R/a).
              C1 = SQRT(A1 + SQRT(A2))
          ELSE
              C1 = SQRT(A1 - SQRT(A2))
          END IF
          CK = 1.0
      ELSE
*       Adopt C1*C2 = CK for difficult case (Seppo's suggestion of 1991).
          C1 = 1.0
          CK = ZMB/SQRT(-8.0D0*HM(IM)*RB*UPR2)
          WRITE (6,20)  IM, KSTARM(IM), RB, HM(IM), UPR2, A2, CK-1.0
   20     FORMAT (' WARNING!    HIRECT    IM K* R H UPR2 A2 CK-1 ',
     &                                    2I4,1P,5E10.2)
          ITER = ITER + 1
      END IF
*
*       Specify KS coordinate scaling from angular momentum conservation.
      C2 = CK/C1
*
*       Transform KS variables to yield the prescribed elements.
      DO 25 K = 1,4
          UM(K,IM) = C2*UM(K,IM)
          UMDOT(K,IM) = C1*UMDOT(K,IM)
   25 CONTINUE
*
*       Improve solution by second iteration in case of CK procedure.
      ITER = ITER + 1
      IF (ITER.EQ.2) GO TO 10
*
      RETURN
*
      END
