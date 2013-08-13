      SUBROUTINE ECCMOD(I,ITERM)
*
*
*       Eccentricity modulation of hierarchical binary.
*       -----------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      REAL*8  BODYI(2),W(2)
      DATA  ITIME,IDELAY /1,0/
      SAVE  ITIME,IDELAY
*
*
*       Determine merger & ghost index.
      ITERM = 0
      CALL FINDJ(I,IGHOST,IM)
*
*       Quit for tidal dissipation or circular orbit (TMDIS set by IMPACT).
      IF (KSTARM(IM).EQ.-2.OR.KSTARM(IM).GE.10) THEN
          GO TO 50
      END IF
*
*       Initialize delay indicator and pair index.
      IQ = 0
      IPAIR = I - N
*
*       Skip on hyperbolic outer orbit, double merger or circularized orbit.
      IF (H(IPAIR).GT.0.0.OR.NAME(I).LT.-2*NZERO.OR.
     &    KSTARM(IM).EQ.10) THEN
          TMDIS(IM) = TIME + 10.0/TSTAR
          GO TO 50
      END IF
*
*       Resolve coordinates & velocities (first time only).
      CALL RESOLV(IPAIR,1)
*
*       Form inner distance, square KS velocity and radial velocity term.
      RI = 0.0
      V20 = 0.0
      TD2 = 0.0
      DO 5 K = 1,4
          RI = RI + UM(K,IM)**2
          V20 = V20 + UMDOT(K,IM)**2
          TD2 = TD2 + 2.0*UM(K,IM)*UMDOT(K,IM)
    5 CONTINUE
*
*       Evaluate inner semi-major axis and eccentricity.
      ZMB = CM(1,IM) + CM(2,IM)
      SEMI = -0.5*ZMB/HM(IM)
      IF (SEMI.LE.0.0) THEN
          TMDIS(IM) = TIME + 1.0
          WRITE(3,*)' ECCMOD ERROR ',SEMI,ECC,H(IPAIR),HM(IM),ZMB
          GO TO 50
      END IF
      ECC2 = (1.0 - RI/SEMI)**2 + TD2**2/(ZMB*SEMI)
      ECC = SQRT(ECC2)
*
*       Obtain growth time and modify KS elements from de/dt & dh/dt.
      I1 = 2*IPAIR - 1
      CALL HIGROW(I1,IGHOST,IM,ECC,SEMI,EMAX,EMIN,TG,EDAV,ZI,IQ)
*
*       Check termination for new CHAOS & SPIRAL or collision.
      IF (IQ.LT.0) THEN
          ITERM = 1
          ITIME = 1
          GO TO 50
      END IF
*
*       Delay for long growth times or aged active SPIRAL.
      IF ((KSTARM(IM).GE.0.AND.TG.GT.20.0).OR.IQ.GT.0.OR.
     &    (KSTARM(IM).EQ.-2.AND.MAX(ECC,EMAX).LT.0.9).OR.
     &    (KSTARM(IM).LT.0.AND.ECC.LT.0.1)) THEN
          RM = MAX(RADIUS(I1),RADIUS(IGHOST),1.0D-20)
          IDELAY = IDELAY + 1
          IF (EMAX.GT.0.99.AND.MOD(IDELAY,10).EQ.0) THEN
              ALPH = 360.0*ZI/TWOPI
              WRITE (6,10)  NAME(I1), IQ, KSTARM(IM), LIST(1,I1), ECC,
     &                      EMAX, TG, SEMI*(1.0-ECC)/RM, ALPH
          END IF
   10     FORMAT (' ECCMOD DELAY    NAM IQ K* NP E EMAX TG QP/R* IN ',
     &                              I6,3I4,2F8.4,F8.3,2F7.1)
          DT = 10.0
          IF (LIST(1,I1).GT.0.AND.SEMI*(1.0-ECC).LT.5.0*RM) THEN
              DT = 1.0
          END IF
          TMDIS(IM) = TIME + MAX(TG,DT)/TSTAR
          ITIME = 1
          GO TO 50
      END IF
*
*       Estimate current t_{circ} and de/dt for relevant condition.
      PMIN = SEMI*(1.0 - ECC)
      RM = MAX(RADIUS(I1),RADIUS(IGHOST),1.0D-20)
      IF (PMIN.LT.50.0*RM) THEN
          BODYI(1) = CM(1,IM)
          BODYI(2) = CM(2,IM)
          CALL HICIRC(PMIN,ECC,I1,IGHOST,BODYI,TG,TC,EC,EDT,W)
      ELSE
          TC = 1.0D+10
          EDT = 1.0D-10
      END IF
*
*       Include diagnostics every 5000th time.
      IF (ITIME.EQ.1.OR.MOD(ITIME,5000).EQ.0) THEN
          A1 = -0.5*BODY(I)/H(IPAIR)
          E2 = (1.0 - R(IPAIR)/A1)**2 + TDOT2(IPAIR)**2/(A1*BODY(I))
          ECC1 = SQRT(E2)
          ZID = 360.0*ZI/TWOPI
          NP = LIST(1,I1)
          YC = R0(IPAIR)/SEMI
          WRITE (6,15)  NAME(I1), NP, TTOT, ECC, EMAX, ECC1, PMIN/RM,
     &                  YC, TG, TC, EDAV, ZID
   15     FORMAT (' ECCMOD    NM NP T E EX E1 QP/R* PC/A TG TC EDA IN ',
     &                        I6,I4,F11.4,3F8.4,2F7.1,1P,3E9.1,0P,F9.3)
          CALL FLUSH(3)
      END IF
*
      NEMOD = NEMOD + 1
      ITIME = ITIME + 1
*
*       Include termination for expected short circularization time.
      IF (KSTARM(IM).GE.0.AND.PMIN.LT.3.0*RM) THEN
          ZID = 360.0*ZI/TWOPI
          WRITE (6,40)  ITIME, EMAX, ECC, SEMI*(1.0 - ECC)/RM, ZID, TC
   40     FORMAT (' ECCMOD TERM    IT EX E QP/R IN TC ',
     &                             I5,2F9.5,F6.2,F8.2,F8.1)
          CALL FLUSH(3)
          ITERM = 1
          ITIME = 1
      END IF
*
   50 RETURN
*
      END
