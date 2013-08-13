      SUBROUTINE HMDOT2(J,IGHOST,IMERGE,M1,KW,MC,DMS,RNEW,ITERM)
*
*
*       Mass loss from outer hierarchical binary.
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
      ICM = N + IPAIR
      JPAIR = IGHOST - N
      I1 = 2*JPAIR - 1
      I2 = I1 + 1
*
*       Include quitting conditions to avoid large changes at end-point.
      IQUIT = .FALSE.
*     IF (RADIUS(J)*SU.GT.500.0) IQUIT = .TRUE.
*     IF (BODY0(J)*ZMBAR.GT.15.0.AND.KSTAR(J).GE.4) IQUIT = .TRUE.
*     IF (ABS(M1 - MC).LT.0.1*M1.OR.MC.GT.1.4) IQUIT = .TRUE.
*     IF (BODY0(J)*ZMBAR - M1.GT.0.5*M1) IQUIT = .TRUE.
*       Note that tidal dissipation needs full KS for rectification in MDOT.
*     IF (KSTAR(IGHOST).LT.0.OR.BODY(ICM).EQ.0.0D0) IQUIT = .TRUE.
      IF (BODY(ICM).EQ.0.0D0) IQUIT = .TRUE.
*
*       Quit for mis-identification or advanced evolution.
      IF (IGHOST.LE.N.OR.IQUIT) THEN
          IF (IGHOST.LE.N) THEN
              if(rank.eq.0)
     &        WRITE (6,5)  NAME(I1), IGHOST
    5         FORMAT (' DANGER!    HMDOT2    NAM IG ',2I5)
              STOP
          END IF
          ITERM = 1
          GO TO 50
      END IF
*
*       Update radius for #J and copy #I2 for MDOT & HCORR in case J = I2.
      RSTAR = RNEW
      IF (J.GT.I1) THEN
          RADIUS(J) = RNEW
          RNEW = RADIUS(I2)
*       Return RNEW for copying back to RADIUS(I2) on DM/M < 0.01.
      END IF
*
*       Skip further modifications on zero mass loss of same type.
      IF (ABS(DMS/M1).EQ.0.0D0.AND.KW.EQ.KSTAR(J)) GO TO 50
*
*       Set two-body elements for outer orbit and inner semi-major axis.
      SEMI = -0.5*BODY(ICM)/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(ICM)*SEMI)
      ECC1 = SQRT(ECC2)
      ZMB0 = CM(3,IMERGE) + CM(4,IMERGE)
*       Note old H(JPAIR) used with new inner SEMI0 (small effect).
      SEMI0 = -0.5*ZMB0/H(JPAIR)
*
*       Form modified mass ratio and semi-major axes from M*A = const.
      DM = DMS/ZMBAR
      ZMB = ZMB0 - DM
      SEMI1 = SEMI*(ZMB0 + BODY(2*IPAIR-1))/(ZMB + BODY(2*IPAIR-1))
      SEMI2 = SEMI0*ZMB0/ZMB
      PMIN = SEMI1*(1.0 - ECC1)
*
*       Evaluate old separation, square regularized velocity, t'' & ECC.
      RI = 0.0D0
      V20 = 0.0
      TD2 = 0.0
      DO 10 K = 1,4
          RI = RI + U(K,JPAIR)**2
          V20 = V20 + UDOT(K,JPAIR)**2
          TD2 = TD2 + 2.0*U(K,JPAIR)*UDOT(K,JPAIR)
   10 CONTINUE
      ECC2 = (1.0 - RI/SEMI0)**2 + TD2**2/(SEMI0*ZMB0)
      ECC = SQRT(ECC2)
*
*       Obtain stability parameters for the new configuration.
      IF (ECC1.LT.1.0) THEN
          NST = NSTAB(SEMI2,SEMI1,ECC,ECC1,ZI,CM(1,IMERGE),
     &                                     CM(2,IMERGE),BODY(2*IPAIR))
          IF (NST.EQ.0) THEN
              PCRIT = 0.99*PMIN
          ELSE
              PCRIT = 1.01*PMIN
          END IF
      ELSE
          PCRIT = stability(CM(3,IMERGE),CM(4,IMERGE),BODY(2*IPAIR-1),
     &                                          ECC,ECC1,0.0D0)*SEMI2
      END IF
*
*       Update pericentre distance on successful stability test or exit.
      IF (PMIN.GT.PCRIT.AND.H(JPAIR).LT.-ECLOSE) THEN
          R0(IPAIR) = PCRIT
      ELSE
          ITERM = 1
          GO TO 50
      END IF
*
*       Check for possible circularization (tidal or sequential).
      RP = SEMI2*(1.0 - ECC)
      IF (KZ(27).GT.0.AND.RP.LT.4.0*RADIUS(J).AND.ECC.GT.0.002) THEN
          TC = 1.0D+10
          IF (KZ(27).EQ.2) THEN
              J1 = 2*JPAIR - 1
              J2 = J1 + 1
              DO 14 K = 1,2
                  BODYI(K) = CM(K+2,IMERGE)
   14         CONTINUE
              TG = 1.0
*       Evaluate the circularization time.
              CALL HICIRC(RP,ECC,J1,J2,BODYI,TG,TC,EC,ED,W)
              IF (TC.GT.2000) GO TO 18
          END IF
          if(rank.eq.0)
     &    WRITE (6,15)  NAME(J), KSTAR(J), KSTARM(IMERGE), ECC, DMS,
     &                  RP, RADIUS(J), TC
   15     FORMAT (' HMDOT2 TERM    NM K* E DM RP R* TC ',
     &                             I6,2I4,F8.4,F7.3,1P,3E10.2)
          ITERM = 1
          GO TO 50
      END IF
*
*       Obtain energy change from M*A = const and H = -M/(2*A) (DCH 8/96).
   18 DH = DM/SEMI0*(1.0 - 0.5*DM/ZMB0)
      HI = H(JPAIR)
      H(JPAIR) = H(JPAIR) + DH
*
*       Form KS coordinate & velocity scaling factors (general point is OK).
      SEMI2 = -0.5*ZMB/H(JPAIR)
      C2 = SQRT(SEMI2/SEMI0)
      V2 = 0.5*(ZMB + H(JPAIR)*RI*(SEMI2/SEMI0))
      C1 = SQRT(V2/V20)
*
*       Re-scale KS variables to new energy (H < 0: constant eccentricity).
      R(JPAIR) = 0.0D0
      DO 20 K = 1,4
          U(K,JPAIR) = C2*U(K,JPAIR)
          UDOT(K,JPAIR) = C1*UDOT(K,JPAIR)
          R(JPAIR) = R(JPAIR) + U(K,JPAIR)**2
   20 CONTINUE
*
*       Reduce mass of relevant component (ghost is original second member).
      ZMU0 = CM(3,IMERGE)*CM(4,IMERGE)/ZMB0
      KM = 3
      IF (J.EQ.2*JPAIR) KM = 4
      CM(KM,IMERGE) = CM(KM,IMERGE) - DM
*
*       Include corrections to EMERGE & EMDOT (consistency; no net effect!).
      ZMU = CM(3,IMERGE)*CM(4,IMERGE)/ZMB
      DECORR = ZMU*H(JPAIR) - ZMU0*HI
      EMERGE = EMERGE + DECORR
      EMDOT = EMDOT - DECORR
      EGRAV = EGRAV - DECORR
*
*       Correct outer orbit for mass loss (use #I2 for consistency).
      CALL HCORR(I2,DM,RSTAR)
*
*       Print some diagnostics on non-zero mass loss.
      IF (rank.eq.0.and.DMS.GT.1.0D-03) THEN
          WRITE (6,30)  NAME(J), KSTAR(J), IMERGE, KM, M1, DMS, PMIN,
     &                  PCRIT, SEMI2, SEMI1, R(JPAIR), H(JPAIR), DECORR
   30     FORMAT (' HMDOT2    NAM K* IM KM M1 DM PM PC A A1 R H DE ',
     &                        I6,3I4,F6.2,1P,7E10.2,0P,F10.6)
      END IF
*
   50 RETURN
*
      END
