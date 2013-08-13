      SUBROUTINE BINDAT
*
*
*       Binary data bank.
*       -----------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      REAL*4  EB(KMAX),ECC(KMAX),RCM(KMAX),ECM(KMAX),PB(KMAX),AS(30)
      REAL*8  XX(3,3),VV(3,3)
*
*       Form binding energy and central distance for each KS pair.
      ZMBIN = 0.0
      DO 10 JPAIR = 1,NPAIRS
          J2 = 2*JPAIR
          J1 = J2 - 1
          ICM = N + JPAIR
          ZMBIN = ZMBIN + BODY(ICM)
          BODYCM = BODY(ICM)
*       Determine merger & ghost index for negative c.m. name (skip ghost).
          IF (NAME(ICM).LT.0.AND.BODY(ICM).GT.0.0) THEN
              CALL FINDJ(J1,J,IMERGE)
*       Employ actual masses and two-body distance for energy & eccentricity.
              BODYCM = CM(1,IMERGE) + CM(2,IMERGE)
              EB(JPAIR) = CM(1,IMERGE)*CM(2,IMERGE)*HM(IMERGE)/BODYCM
              SEMI = -0.5*BODYCM/HM(IMERGE)
              RJ = SQRT(XREL(1,IMERGE)**2 + XREL(2,IMERGE)**2 +
     &                                      XREL(3,IMERGE)**2)
*       Assume that merged binary is near apo or peri (hence ignore TDOT2).
              ECC2 = (1.0 - RJ/SEMI)**2
*       Include separate diagnostics for the hierarchy (inner comps J1 & J).
              SEMI1 = -0.5*BODY(ICM)/H(JPAIR)
              ECC1 = (1.0 - R(JPAIR)/SEMI1)**2 +
     &                                 TDOT2(JPAIR)**2/(BODY(ICM)*SEMI1)
              E0 = SQRT(ECC2)
              E1 = SQRT(ECC1)
              IF (J.LT.0) J = J1
              RM = SEMI*(1.0 - E0)/MAX(RADIUS(J1),RADIUS(J),1.0D-20)
              RM = MIN(RM,99.9D0)
              P0 = DAYS*SEMI*SQRT(ABS(SEMI)/BODYCM)
              P1 = DAYS*SEMI1*SQRT(ABS(SEMI1)/BODY(ICM))
              DO 2 K = 1,3
                  XX(K,1) = XREL(K,IMERGE)
                  XX(K,2) = 0.0
                  XX(K,3) = X(K,J2)
                  VV(K,1) = VREL(K,IMERGE)
                  VV(K,2) = 0.0
                  VV(K,3) = XDOT(K,J2)
    2         CONTINUE
              CALL INCLIN(XX,VV,X(1,ICM),XDOT(1,ICM),ALPH)
              PCR = stability(CM(1,IMERGE),CM(2,IMERGE),BODY(ICM),E0,E1,
     &                                                       ALPH)*SEMI
              PM = SEMI1*(1.0 - E1)/PCR
              if(rank.eq.0)
     &        WRITE (84,3) TTOT, NAME(J1), NAME(J), KSTAR(J1), KSTAR(J),
     &                     KSTARM(IMERGE), E0, E1, PM, RM, P0, P1, SEMI1
    3         FORMAT (' BINDAT:    T NM K* E0 E1 PM/PC PM0/R* P0 P1 A1',
     &                             F8.1,2I5,3I4,2F7.3,2F6.1,1P,3E9.1)
              CALL FLUSH(84)
          ELSE IF (BODY(J1).GT.0.0D0) THEN
*       Form binding energy and eccentricity for standard case.
              EB(JPAIR) = BODY(J1)*BODY(J2)*H(JPAIR)/
     &                                             (BODY(J1) + BODY(J2))
              SEMI = -0.5*BODY(ICM)/H(JPAIR)
              ECC2 = (1.0 - R(JPAIR)/SEMI)**2 +
     &                                  TDOT2(JPAIR)**2/(BODY(ICM)*SEMI)
          ELSE
              IM = 0
*       Search merger table to identify corresponding index of c.m. name.
              DO 5 K = 1,NMERGE
                  IF (NAMEG(K).EQ.NAME(ICM)) THEN
                      IM = K
                  END IF
    5         CONTINUE
              BODYJ1 = CM(3,IM)
              BODYJ2 = CM(4,IM)
              BODYCM = BODYJ1 + BODYJ2
              BODYCM = MAX(BODYCM,1.0D-10)
              EB(JPAIR) = BODYJ1*BODYJ2*H(JPAIR)/BODYCM
              SEMI = -0.5*BODYCM/H(JPAIR)
              ECC2 = (1.0 - SEMI/R(JPAIR))**2
          END IF
          ECC(JPAIR) = SQRT(ECC2)
          EB(JPAIR) = MAX(EB(JPAIR),-9.99999)
          PB(JPAIR) = DAYS*SEMI*SQRT(ABS(SEMI)/BODYCM)
          PB(JPAIR) = MIN(PB(JPAIR),99999.9)
          IF (SEMI.LT.0.0) PB(JPAIR) = 0.0
*       Obtain binding energy (per unit mass) of c.m. motion.
          VJ2 = XDOT(1,ICM)**2 + XDOT(2,ICM)**2 + XDOT(3,ICM)**2
          IF (BODY(ICM).EQ.0.0D0) VJ2 = 0.0
          POTJ = 0.0
          DO 9 J = IFIRST,NTOT
              IF (J.EQ.ICM) GO TO 9
              RIJ2 = (X(1,ICM) - X(1,J))**2 + (X(2,ICM) - X(2,J))**2 +
     &                                        (X(3,ICM) - X(3,J))**2 
              POTJ = POTJ + BODY(J)/SQRT(RIJ2)
    9     CONTINUE
          ECM(JPAIR) = 0.5*VJ2 - POTJ
*       Check for external tidal field (note that HT includes mass).
          IF (KZ(14).GT.0) THEN
              CALL XTRNLV(ICM,ICM)
              ECM(JPAIR) = ECM(JPAIR) + HT/(BODY(ICM) + 1.0E-20)
          END IF
          RCM(JPAIR) = SQRT((X(1,ICM) - RDENS(1))**2 +
     &                      (X(2,ICM) - RDENS(2))**2 +
     &                      (X(3,ICM) - RDENS(3))**2)
          RCM(JPAIR) = MIN(RCM(JPAIR),99.9)
   10 CONTINUE
*
*       Copy relevant binary diagnostics to single precision.
      AS(1) = TIME + TOFF
      AS(2) = RSCALE
      AS(3) = RTIDE
      AS(4) = RC
      AS(5) = TPHYS
      AS(6) = -1.5*(TIDAL(1)*ZMASS**2)**0.3333
      AS(7) = 0.0
      DO 20 K = 1,10
          AS(K+7) = E(K)
   20 CONTINUE
      AS(18) = SBCOLL
      AS(19) = BBCOLL
      AS(20) = ZKIN
      AS(21) = POT
      AS(22) = EBIN0
      AS(23) = EBIN
      AS(24) = ESUB
      AS(25) = EMERGE
      AS(26) = BE(3)
      AS(27) = ZMASS
      AS(28) = ZMBIN
      AS(29) = CHCOLL
      AS(30) = ECOLL
*
*       Write formatted data bank on unit 9.
*
      if(rank.eq.0)then
      WRITE (9,30)  NPAIRS, MODEL, NRUN, N, NC, NMERGE, (AS(K),K=1,7)
   30 FORMAT (3I4,I6,2I4,2X,F7.1,2F7.2,F7.3,F8.1,2F9.4)
      WRITE (9,35)  (AS(K),K=8,17)
   35 FORMAT (10F11.6)
      WRITE (9,40)  (AS(K),K=18,30)
   40 FORMAT (13F10.5)
      end if
*
      DO 50 JPAIR = 1,NPAIRS
          J1 = 2*JPAIR - 1
          J2 = 2*JPAIR
          KCM = KSTAR(N+JPAIR)
          IF (NAME(N+JPAIR).LT.0) THEN
              KCM = -10
          END IF
          ICM = N + JPAIR
          SEMI = -0.5*(BODY(J1) + BODY(J2))/H(JPAIR)
          ZN = SQRT((BODY(J1) + BODY(J2))/SEMI**3)
          RP = U(1,JPAIR)**2 + U(2,JPAIR)**2 + U(3,JPAIR)**2 +
     &                                         U(4,JPAIR)**2
          RI = SQRT((X(1,ICM) - RDENS(1))**2 +
     &              (X(2,ICM) - RDENS(2))**2 +
     &              (X(3,ICM) - RDENS(3))**2)
          VI = SQRT(XDOT(1,ICM)**2 + XDOT(2,ICM)**2 + XDOT(3,ICM)**2)
          if(rank.eq.0)
     &    WRITE (9,45)  TTOT, EB(JPAIR), ECC(JPAIR), ECM(JPAIR),
     &                  BODY(J1)*ZMBAR, BODY(J2)*ZMBAR, PB(JPAIR),
     &                  SEMI, ZN, RP, STEP(J1), RI, VI,
     &                  NAME(J1), NAME(J2), NAME(N+JPAIR), 
     &                  KSTAR(J1), KSTAR(J2), KCM
   45     FORMAT (1P,13E13.5,3I8,3I4)
   50 CONTINUE
      CALL FLUSH(9)
*
      RETURN
*
      END
