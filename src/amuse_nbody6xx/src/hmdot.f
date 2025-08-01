      SUBROUTINE HMDOT(J,IMERGE,M1,KW,MC,DMS,RNEW,ITERM)
*
*
*       Mass loss from inner hierarchical binary.
*       -----------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      REAL*8  BODYI(2),W(2)
      REAL*8  M1,MC
      LOGICAL IQUIT
*
*
*       Define merger termination index and relevant KS/c.m. indices.
      ITERM = 0
      IPAIR = KSPAIR
      I1 = 2*IPAIR - 1
      ICM = N + IPAIR
*
*       Include quitting conditions to avoid large changes at end-point.
      IQUIT = .FALSE.
*     IF (RADIUS(J)*SU.GT.500.0) IQUIT = .TRUE.
*     IF (BODY0(J)*ZMBAR.GT.15.0.AND.KSTAR(J).GE.4) IQUIT = .TRUE.
*     IF (ABS(M1 - MC).LT.0.1*M1.OR.MC.GT.1.4) IQUIT = .TRUE.
*     IF (BODY0(J)*ZMBAR - M1.GT.0.5*M1) IQUIT = .TRUE.
*       Note that tidal dissipation needs full KS for rectification in MDOT.
*     IF (KSTARM(IMERGE).LT.0.OR.BODY(ICM).EQ.0.0D0) IQUIT = .TRUE.
      IF (BODY(ICM).EQ.0.0D0) IQUIT = .TRUE.
*
*       Quit for mis-identification, advanced evolution or double merger.
      IF (J.LE.0.OR.IQUIT) THEN
          ITERM = 1
          GO TO 50
      END IF
*
*       Update radius for #J and copy #I1 for MDOT & HCORR in case of ghost.
      RSTAR = RNEW
      IF (J.GT.I1) THEN
          RADIUS(J) = RNEW
          RSTAR = RADIUS(I1)
*       Return RNEW for copying back to RADIUS(I1) on ghost & DM/M < 0.01.
      END IF
*
*       Skip further modifications on zero mass loss of same type.
      IF (ABS(DMS/M1).EQ.0.0D0.AND.KW.EQ.KSTAR(J)) GO TO 50
*
*       Set two-body elements for outer orbit and inner semi-major axis.
      SEMI = -0.5*BODY(ICM)/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(ICM)*SEMI)
      ECC1 = SQRT(ECC2)
      ZMB0 = CM(1,IMERGE) + CM(2,IMERGE)
      SEMI0 = -0.5*ZMB0/HM(IMERGE)
*
*       Form modified mass ratio and semi-major axes from M*A = const.
      DM = DMS/ZMBAR
      ZMB = ZMB0 - DM
      SEMI1 = SEMI*(ZMB0 + BODY(2*IPAIR))/(ZMB + BODY(2*IPAIR))
      SEMI2 = SEMI0*ZMB0/ZMB
      PMIN = SEMI1*(1.0 - ECC1)
*
*       Evaluate old separation, square regularized velocity, t'' & ECC.
      RI = 0.0D0
      V20 = 0.0
      TD2 = 0.0
      DO 10 K = 1,4
          RI = RI + UM(K,IMERGE)**2
          V20 = V20 + UMDOT(K,IMERGE)**2
          TD2 = TD2 + 2.0*UM(K,IMERGE)*UMDOT(K,IMERGE)
   10 CONTINUE
      ECC2 = (1.0 - RI/SEMI0)**2 + TD2**2/(SEMI0*ZMB0)
      ECC = SQRT(ECC2)
*
*       Determine inclination (use #I1 as first KS component).
      CALL HIMAX(I1,IMERGE,ECC,SEMI0,EMAX,EMIN,ZI,TG,EDAV)
*
*       Evaluate the general stability function (Mardling 2008).
      IF (ECC1.LT.1.0) THEN
          EOUT = ECC1
*       Increase tolerance near sensitive stability boundary (RM 10/2008).
          IF (EOUT.GT.0.9) THEN
              DE = 0.5*(1.0 - EOUT)
              DE = MIN(DE,0.01D0)
*       Add extra amount 0.011 to avoid switching.
              IF (ECC1.GT.0.9) DE = DE + 0.011
              EOUT = EOUT - DE
              PMIN = SEMI1*(1.0 - EOUT)
          END IF
          NST = NSTAB(SEMI2,SEMI1,ECC,EOUT,ZI,CM(1,IMERGE),
     &                                     CM(2,IMERGE),BODY(2*IPAIR))
          IF (NST.EQ.0) THEN
              PCRIT = 0.99*PMIN
              PCR = stability(CM(1,IMERGE),CM(2,IMERGE),BODY(2*IPAIR),
     &                                          ECC,ECC1,ZI)*SEMI2
          ELSE
              PCRIT = 1.01*PMIN
          END IF
      ELSE
          PCRIT = stability(CM(1,IMERGE),CM(2,IMERGE),BODY(2*IPAIR),
     &                                          ECC,ECC1,ZI)*SEMI2
      END IF
*
*       Update pericentre distance on successful stability test or exit.
      IF (PMIN.GT.PCRIT) THEN
          R0(IPAIR) = PCRIT
      ELSE
          ITERM = 1
          GO TO 50
      END IF
*
*       Check for possible circularization (tidal or sequential).
      RP = SEMI2*(1.0 - ECC)
      IF (KZ(27).GT.0.AND.RP.LT.4.0*RADIUS(J).AND.ECC.GT.0.002.AND.
     &    KSTARM(IMERGE).GE.0) THEN
          TC = 1.0D+10
          IF (KZ(27).EQ.2) THEN
              I2 = 0
*       Identify the ghost by searching single bodies.
              DO 12 K = IFIRST,N
                  IF (BODY(K).EQ.0.0D0.AND.
     &                NAME(K).EQ.NAMEG(IMERGE)) THEN
                      I2 = K
                  END IF
   12         CONTINUE
              IF (I2.GT.0) THEN
                  DO 14 K = 1,2
                      BODYI(K) = CM(K,IMERGE)
   14             CONTINUE
                  TG = 1.0
*       Evaluate the circularization time.
                  CALL HICIRC(RP,ECC,I1,I2,BODYI,TG,TC,EC,ED,W)
                  IF (TC.GT.2000) GO TO 18
              END IF
          END IF
          if(rank.eq.0)
     &    WRITE (6,15)  NAME(J), KSTAR(J), KSTARM(IMERGE), ECC, DMS,
     &                  RP, RADIUS(J), TC
   15     FORMAT (' HMDOT TERM    NAM K* E DM RP R* TC ',
     &                            I6,2I4,F8.4,F7.3,1P,3E10.2)
          ITERM = 1
          GO TO 50
      ELSE IF (ECC.GT.0.002) THEN
          IF (RP.LT.2.0*RADIUS(J)) THEN
              if(rank.eq.0)
     &        WRITE (6,15)  NAME(J), KSTAR(J), KSTARM(IMERGE), ECC, DMS,
     &                      RP, RADIUS(J)
              ITERM = 1
              GO TO 50
          END IF
      END IF
*
*       Obtain energy change from M*A = const and H = -M/(2*A) (DCH 8/96).
   18 DH = DM/SEMI0*(1.0 - 0.5*DM/ZMB0)
      HM0 = HM(IMERGE)
      HM(IMERGE) = HM(IMERGE) + DH
*
*       Form KS coordinate & velocity scaling factors (general point is OK).
      SEMI2 = -0.5*ZMB/HM(IMERGE)
      C2 = SQRT(SEMI2/SEMI0)
      V2 = 0.5*(ZMB + HM(IMERGE)*RI*(SEMI2/SEMI0))
      C1 = SQRT(V2/V20)
*
*       Re-scale KS variables to new energy (H < 0: constant eccentricity).
      DO 20 K = 1,4
          UM(K,IMERGE) = C2*UM(K,IMERGE)
          UMDOT(K,IMERGE) = C1*UMDOT(K,IMERGE)
*       Note no need to update XREL & VREL for RESTART (c.m. error cancels).
   20 CONTINUE
*
*       Reduce mass of relevant component (ghost is original second member).
      ZMU0 = CM(1,IMERGE)*CM(2,IMERGE)/ZMB0
      KM = 1
      IF (J.GE.IFIRST) KM = 2
      CM(KM,IMERGE) = CM(KM,IMERGE) - DM
*
*       Include corrections to EMERGE & EMDOT (consistency; no net effect!).
      ZMU = CM(1,IMERGE)*CM(2,IMERGE)/ZMB
      DECORR = ZMU*HM(IMERGE) - ZMU0*HM0
      EMERGE = EMERGE + DECORR
      EMDOT = EMDOT - DECORR
*       Compensate for EMERGE contribution in DEGRAV (cf. EVENTS).
      EGRAV = EGRAV - DECORR
*
*       Correct outer orbit for mass loss (use #I1 in case #J is a ghost).
      CALL HCORR(I1,DM,RSTAR)
*
*       Print some diagnostics on non-zero mass loss.
      IF (rank.eq.0.and.DMS.GT.1.0D-03) THEN
          WRITE (6,30)  NAME(J), KSTAR(J), IMERGE, KM, M1, DMS, PMIN,
     &                  PCRIT, SEMI2, SEMI1, DECORR
   30     FORMAT (' HMDOT    NAM K* IM KM M1 DM PM PC A A1 DE ',
     &                       I6,3I4,F6.2,1P,5E10.2,0P,F10.6)
      END IF
*
   50 RETURN
*
      END
