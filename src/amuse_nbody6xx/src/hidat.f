      SUBROUTINE HIDAT
*
*
*       Hierarchical data bank.
*       -----------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      REAL*8  XX(3,3),VV(3,3),M1,M2,M3
*     REAL*4  EB(KMAX),ECM(KMAX)
      LOGICAL  FIRST
      SAVE  FIRST
      DATA  FIRST /.TRUE./
*
*
*       Write formatted data bank on unit 87.
      IF (FIRST.and.rank.eq.0) THEN
          OPEN (UNIT=87,STATUS='UNKNOWN',FORM='FORMATTED',FILE='HIDAT')
          FIRST = .FALSE.
          WRITE (87,1)
    1     FORMAT (/,'  NAM1  NAM2  NAM3  K*       M1   M2   M3    RI',
     &            '    EMAX   E0    E1      P0    P1')
      END IF
*
      IF (NMERGE.GT.0) THEN
*       Count higher-order systems.
          MULT = 0
          DO 2 I = N+1,NTOT
              IF (NAME(I).LT.-2*NZERO) MULT = MULT + 1
    2     CONTINUE
          if(rank.eq.0)
     &    WRITE (87,3)  NPAIRS, NRUN, N, NC, NMERGE, MULT, NEWHI, TTOT
    3     FORMAT (/,I6,I4,I6,3I4,I6,F9.1)
      END IF
*
*       Produce output for each merged KS pair (TRIPLE, QUAD, .. [[B,B],B]).
      IPAIR = 0
*     ISTAB = 0
      DO 30 JPAIR = 1,NPAIRS
          ICM = N + JPAIR
          IF (NAME(ICM).GT.0.AND.BODY(ICM).GT.0.0D0) GO TO 30
          J2 = 2*JPAIR
          J1 = J2 - 1
*       Determine merger & ghost index (delay ghost c.m.).
          CALL FINDJ(J1,J,IM)
          KCM = KSTARM(IM)
          IF (BODY(ICM).GT.0.0D0) THEN
              SEMI1 = -0.5*BODY(ICM)/H(JPAIR)
              ECC1 = (1.0 - R(JPAIR)/SEMI1)**2 +
     &                                 TDOT2(JPAIR)**2/(BODY(ICM)*SEMI1)
          END IF
*
*       Consider any multiple hierarchies first (c.m. mass > 0 also here).
          IF (NAME(ICM).LT.-2*NZERO) THEN
              IM = 0
*       Find the original merger and ghost index (J > N possible).
              DO 5 K = 1,NMERGE
                  IF (NAMEM(K).EQ.NAME(ICM) + 2*NZERO) IM = K
    5         CONTINUE
              IF (IM.EQ.0) GO TO 30
              J = 0
              DO 8 K = IFIRST,NTOT
                  IF (NAMEG(IM).EQ.NAME(K)) J = K
    8         CONTINUE
              IF (J.EQ.0) GO TO 30
              IPAIR = IPAIR + 1
*       Employ actual masses and two-body distance for energy & eccentricity.
              BODYCM = CM(1,IM) + CM(2,IM)
*             EB(IPAIR) = CM(1,IM)*CM(2,IM)*HM(IM)/BODYCM
              SEMI = -0.5*BODYCM/HM(IM)
              RJ = SQRT(XREL(1,IM)**2 + XREL(2,IM)**2 + XREL(3,IM)**2)
              TD2 = 0.0
              DO 10 K = 1,4
                  TD2 = TD2 + 2.0*UM(K,IM)*UMDOT(K,IM)
   10         CONTINUE
              ECC2 = (1.0 - RJ/SEMI)**2 + TD2**2/(BODYCM*SEMI)
              M1 = CM(1,IM)*SMU
              M2 = CM(2,IM)*SMU
              M3 = BODY(J2)*SMU
*       Determine EMAX and other relevant parameters.
              E0 = SQRT(ECC2)
              CALL HIMAX(J1,IM,E0,SEMI,EMAX,EMIN,ZI,TG,EDAV)
              NAM2 = NAME(J)
*       Check standard case of triple or quadruple (NAME(J2) <= N or > N).
          ELSE IF (BODY(ICM).GT.0.0D0) THEN
              IPAIR = IPAIR + 1
*       Employ actual masses and two-body distance for energy & eccentricity.
              BODYCM = CM(1,IM) + CM(2,IM)
*             EB(IPAIR) = CM(1,IM)*CM(2,IM)*HM(IM)/BODYCM
              SEMI = -0.5*BODYCM/HM(IM)
              RJ = SQRT(XREL(1,IM)**2 + XREL(2,IM)**2 + XREL(3,IM)**2)
              TD2 = 0.0
              DO 12 K = 1,4
                  TD2 = TD2 + 2.0*UM(K,IM)*UMDOT(K,IM)
   12         CONTINUE
              ECC2 = (1.0 - RJ/SEMI)**2 + TD2**2/(BODYCM*SEMI)
*       Include separate diagnostics for the hierarchy (inner comps J1 & J).
              M1 = CM(1,IM)*SMU
              M2 = CM(2,IM)*SMU
              M3 = BODY(J2)*SMU
*       Retain the next part just in case.
              E0 = SQRT(ECC2)
              E1 = SQRT(ECC1)
*       Determine EMAX and other relevant parameters (if needed).
              CALL HIMAX(J1,IM,E0,SEMI,EMAX,EMIN,ZI,TG,EDAV)
              PMIN = SEMI1*(1.0 - E1)
              IF (J.LT.0) J = J1
              RM = SEMI*(1.0 - E0)/MAX(RADIUS(J1),RADIUS(J),1.0D-20)
              RM = MIN(RM,99.9D0)
*       Obtain inclination between inner relative motion and outer orbit.
              DO 15 K = 1,3
                  XX(K,1) = XREL(K,IM)
                  XX(K,2) = 0.0
                  XX(K,3) = X(K,J2)
                  VV(K,1) = VREL(K,IM)
                  VV(K,2) = 0.0
                  VV(K,3) = XDOT(K,J2)
   15         CONTINUE
              CALL INCLIN(XX,VV,X(1,ICM),XDOT(1,ICM),ANGLE)
              PCR = stability(M1,M2,M3,E0,E1,ANGLE)*SEMI
              NAM2 = NAME(J)
*       Perform stability check.
              IF (PMIN*(1.0 - GAMMA(JPAIR)).LT.PCR.AND.E1.LT.0.96.AND.
     &            LIST(1,J1).GT.0) THEN
*                 ISTAB = JPAIR
                  if(rank.eq.0)
     &            WRITE (6,20)  NAME(J1), NAME(J), 180.0*ANGLE/3.14, E1,
     &                          PMIN, PCR, GAMMA(JPAIR), R0(JPAIR)
   20             FORMAT (' MERGE UNSTAB    NAM INC E1 PM PCR G R0 ',
     &                                      2I6,F7.1,F7.3,1P,4E10.2)
              END IF
          ELSE IF (BODY(J1).EQ.0.0D0) THEN
*       Treat case of ghost binary (JPAIR saved in CM(3,IM) & CM(4,IM).
              IM = 0
*       Find the original merger index (use ghost name for [B,B]).
              DO 22 K = 1,NMERGE
                  IF (NAMEG(K).EQ.NAME(ICM)) IM = K
   22         CONTINUE
              IF (IM.EQ.0) GO TO 30
              KCM = KSTAR(ICM)
              J = 0
*       Locate the first KS component (former c.m. hence subtract NZERO).
              DO 23 K = 1,IFIRST
                  IF (NAME(K) - NZERO.EQ.NAME(J1)) J = K
   23         CONTINUE
              IF (J.EQ.0) GO TO 30
*       Include the case of [[B,S],[B,S]] which requires more work.
              IF (BODY(J).EQ.0.0D0) THEN
                  JM = 0
                  DO 24 K = 1,NMERGE
                      IF (NAMEM(IM).EQ.NAMEG(K)) JM = K
   24             CONTINUE
                  IF (JM.EQ.0) GO TO 30
                  J = 0
*       Employ new ghost identification to find the corresponding c.m index.
                  DO 25 K = N+1,NTOT
                      IF (NAME(K).EQ.NAMEM(JM)) J = K
   25             CONTINUE
                  IF (J.EQ.0) GO TO 30
                  J = 2*(J - N)
*       Note that both the triple masses (M1,M2) are saved in CM(1->4,JM).
              END IF
              IPAIR = IPAIR + 1
              NAM2 = NAME(J2)
*       Set indices for the whole system (KS pair, c.m.; note NAME(J2) < 0).
              JP = KVEC(J)
              ICM = N + JP
              J2 = ICM
              J = ICM
*       Form outer semi-major axis and eccentricity.
              SEMI1 = -0.5*BODY(ICM)/H(JP)
              ECC1 = (1.0 - R(JP)/SEMI1)**2 +
     &                                   TDOT2(JP)**2/(BODY(ICM)*SEMI1)
*       Copy masses from second merger component of [B,B] system.
              BODYJ1 = CM(3,IM)
              BODYJ2 = CM(4,IM)
              BODYCM = BODYJ1 + BODYJ2
              BODYCM = MAX(BODYCM,1.0D-10)
              M1 = BODYJ1*SMU
              M2 = BODYJ2*SMU
              M3 = (BODY(ICM) - BODYCM)*SMU
*             EB(IPAIR) = BODYJ1*BODYJ2*H(JPAIR)/BODYCM
              SEMI = -0.5*BODYCM/H(JPAIR)
              ECC2 = (1.0 - R(JPAIR)/SEMI)**2 +
     &                                   TDOT2(JPAIR)**2/(BODYCM*SEMI)
              ANGLE = 0.0
              EMAX = 0.0
          END IF
*
*       Evaluate the potential energy of c.m.
          PHI = 0.0
          DO 26 K = IFIRST,NTOT
              IF (K.EQ.ICM) GO TO 26
              RIJ2 = (X(1,K) - X(1,ICM))**2 + (X(2,K) - X(2,ICM))**2
     &                                      + (X(3,K) - X(3,ICM))**2
              PHI = PHI + BODY(K)/SQRT(RIJ2)
   26     CONTINUE
*
*       Evaluate eccentricities and periods.
          E0 = SQRT(ECC2)
          E1 = SQRT(ECC1)
          P0 = DAYS*SEMI*SQRT(ABS(SEMI)/BODYCM)
          P1 = DAYS*SEMI1*SQRT(ABS(SEMI1)/BODY(ICM))
          P0 = MIN(P0,9999.0D0)
*       Obtain binding energy (per unit mass) of c.m. motion.
          VJ2 = XDOT(1,ICM)**2 + XDOT(2,ICM)**2 + XDOT(3,ICM)**2
          IF (BODY(ICM).EQ.0.0D0) VJ2 = 0.0
*         ECM(IPAIR) = 0.5*VJ2 - PHI
          RI = SQRT((X(1,ICM) - RDENS(1))**2 + (X(2,ICM) - RDENS(2))**2
     &                                       + (X(3,ICM) - RDENS(3))**2)
          if(rank.eq.0)
     &    WRITE (87,28)  NAME(J1), NAM2, NAME(J2), KSTAR(J1), KSTAR(J),
     &                   KCM, M1, M2, M3, RI, EMAX, E0, E1, P0, P1
   28     FORMAT (3I6,3I3,3F5.1,F7.2,F7.3,2F6.2,F8.1,1P,E9.1)
   30 CONTINUE
      CALL FLUSH(87)
*
*       Check merger termination of unstable system.
*     IF (ISTAB.GT.0) THEN
*         KSPAIR = ISTAB
*         IPHASE = 7
*         CALL RESET
*     END IF
*
      RETURN
*
      END
