      SUBROUTINE HIGROW(I,IG,IM,ECC,SEMI,EMAX,EMIN,TG,EDAV,ZI,IQ)
*
*
*       Induced change of hierarchical binary.
*       --------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),YREL(3,MMAX),ZREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &               BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &               RP(NTMAX),ES(NTMAX),ZM(2,NTMAX),IOSC(NTMAX),
     &               NAMEC(NTMAX)
      COMMON/SLOW0/  RANGE,ISLOW(10)
      common/tidal/  cq(2),ct(2),cgr,dedt
      common/rksave/  coeff,HOhat(3),e0,a0,hh,xmb
      REAL*8  A1(3),A2(3),XREL(3),VREL(3),EI(3),HI(3),HO(3),BHAT(3),
     &        UI(4),V(4),EVEC(3),XR0(3),VR0(3),EI0(3),WW(6),W(4)
      REAL*8  BODYI(2),WG(2)
      LOGICAL ICOLL
      DATA  WW  /2.119,3.113,8.175,3.742,4.953,9.413/
      SAVE  ICALL,ITIME,ITRY,NAMEI,DTPREV,TCHECK,ZFAC,PMIN1
      DATA  ICALL,ITRY,NAMEI,TCHECK  /0,0,0,0.0D0/
*
*
*       Copy KS variables to local scalars.
      DO 1 K = 1,4
          UI(K) = UM(K,IM)
          V(K) = UMDOT(K,IM)
    1 CONTINUE
*
*       Transform to physical variables and multiply by 4 (momentum formula).
      CALL KSPHYS(UI,V,XREL,VREL)
      DO 2 K = 1,3
          VREL(K) = 4.0*VREL(K)
    2 CONTINUE
*
*       Specify index for outer body (second KS component) & pair index.
      JCOMP = I + 1
      IPAIR = KVEC(I)
      DTSUM = 0.0
      ITIME = 0
      ICALL = 0
      ECC0 = ECC
      SEMI0 = SEMI
      a0 = SEMI
      HM0 = HM(IM)
      ZM3 = BODY(JCOMP)
      BODYI(1) = CM(1,IM)
      BODYI(2) = CM(2,IM)
*
*       Obtain inner & outer angular momentum by cyclic notation.
    4 A12 = 0.0
      A22 = 0.0
      A1A2 = 0.0
      RI2 = 0.0
      VI2 = 0.0
      RVI = 0.0
      DO 5 K = 1,3
          K1 = K + 1
          IF (K1.GT.3) K1 = 1
          K2 = K1 + 1
          IF (K2.GT.3) K2 = 1
          A1(K) = XREL(K1)*VREL(K2) - XREL(K2)*VREL(K1)
          A2(K) = (X(K1,JCOMP) - X(K1,I))*(XDOT(K2,JCOMP) - XDOT(K2,I))
     &          - (X(K2,JCOMP) - X(K2,I))*(XDOT(K1,JCOMP) - XDOT(K1,I))
          A12 = A12 + A1(K)**2
          A22 = A22 + A2(K)**2
          A1A2 = A1A2 + A1(K)*A2(K)
          RI2 = RI2 + XREL(K)**2
          VI2 = VI2 + VREL(K)**2
          RVI = RVI + XREL(K)*VREL(K)
    5 CONTINUE
*
*       Evaluate orbital parameters for outer orbit from KS elements.
      ZMB = BODY(I) + BODY(JCOMP)
      SEMI1 = -0.5*ZMB/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/SEMI1)**2 + TDOT2(IPAIR)**2/(SEMI1*ZMB)
      ECC1 = SQRT(ECC2)
*       Determine inclination in radians.
      FAC = A1A2/SQRT(A12*A22)
      IF(FAC.LT.-1.D0.OR.FAC.GT.1.D0)THEN
         IF(FAC.LT.-1.05D0.OR.FAC.GT.1.05D0)THEN
            WRITE (6,9)  IPAIR, I, JCOMP, ECC1, SEMI1*SU, FAC 
   9        FORMAT (' WARNING! HIGROW FAC: IP I J E A FAC ',
     &                                   3I6,F7.3,F8.1,F9.5)
         ENDIF
         FAC = MAX(FAC,-1.D0) 
         FAC = MIN(FAC,1.D0) 
      ENDIF
      ZI = ACOS(FAC)
      CC = (1.0 - ECC**2)*COS(ZI)**2
*
*       Construct the Runge-Lenz vector (Heggie & Rasio 1995, Eq.(5)).
      EI2 = 0.0
      DO 10 K = 1,3
          EI(K) = (VI2*XREL(K) - RVI*VREL(K))/BODY(I) -
     &                                                 XREL(K)/SQRT(RI2)
          EI2 = EI2 + EI(K)**2
          EVEC(K) = EI(K)
   10 CONTINUE
      EI2 = MIN(EI2,0.9999d0)
*
*       Define unit vectors for inner eccentricity and angular momenta.
      COSJ = 0.0
      SJSG = 0.0
      DO 15 K = 1,3
          EI(K) = EI(K)/SQRT(EI2)
          HI(K) = A1(K)/SQRT(A12)
          HO(K) = A2(K)/SQRT(A22)
          COSJ = COSJ + HI(K)*HO(K)
          SJSG = SJSG + EI(K)*HO(K)
   15 CONTINUE
*
*       Form unit vector BHAT and scalars AH & BH (Douglas Heggie, 10/9/96).
      AH = 0.0
      BH = 0.0
      DO 16 K = 1,3
          K1 = K + 1
          IF (K1.GT.3) K1 = 1
          K2 = K1 + 1
          IF (K2.GT.3) K2 = 1
          BHAT(K) = HI(K1)*EI(K2) - HI(K2)*EI(K1)
          AH = AH + EI(K)*HO(K)
          BH = BH + BHAT(K)*HO(K)
   16 CONTINUE
*
*       Evaluate the expressions A & Z.
      A = COSJ*SQRT(1.0 - EI2)
      Z = (1.0 - EI2)*(2.0 - COSJ**2) + 5.0*EI2*SJSG**2
*
*       Obtain maximum inner eccentricity (Douglas Heggie, Sept. 1995).
      Z2 = Z**2 + 25.0 + 16.0*A**4 - 10.0*Z - 20.0*A**2 - 8.0*A**2*Z
      EMAX = ONE6*(Z + 1.0 - 4.0*A**2 + SQRT(Z2))
      EMAX = MAX(EMAX,0.0001d0)
      EMAX = SQRT(EMAX)
*
*       Form minimum eccentricity (Douglas Heggie, Sept. 1996).
      AZ = A**2 + Z - 2.0
      IF (AZ.GE.0.0) THEN
          AZ1 = 1.0 + Z - 4.0*A**2
          EMIN2 = ONE6*(AZ1 - SQRT(AZ1**2 - 12.0*AZ))
      ELSE
          EMIN2 = 1.0 - 0.5*(A**2 + Z)
      END IF
      EMIN2 = MAX(EMIN2,0.0001D0)
      EMIN = SQRT(EMIN2)
*
*       Estimate eccentricity growth time-scale (N-body units).
      TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
      TK1 = TWOPI*ABS(SEMI1)*SQRT(ABS(SEMI1)/ZMB)
      TG = TK1**2*ZMB*(1.0 - ECC1**2)**1.5/(BODY(JCOMP)*TK)
*
*       Evaluate numerical precession factor (involves elliptic integral).
      CONST = PFAC(A,Z)
      CONST = CONST*4.0/(1.5*TWOPI*SQRT(6.0))
*
*       Convert growth time to units of 10**6 yrs.
      TG = CONST*TG*TSTAR
*
*       Form doubly averaged eccentricity derivative (Douglas Heggie 9/96).
      YFAC = 15.0*BODY(JCOMP)/(4.0*ZMB)*TWOPI*TK/TK1**2
      YFAC = YFAC*ECC*SQRT(1.0 - ECC**2)/(1.0 - ECC1**2)**(1.5)
      EDAV = YFAC*AH*BH
*
*       Skip small EMAX or large pericentre distance (TMDIS increased).
      QPERI = SEMI*(1.0 - ECC)
      PMIN = SEMI*(1.0 - EMAX)
      RM = MAX(RADIUS(I),RADIUS(IG),1.0D-20)
      IF (KSTARM(IM).GE.0.AND.EMAX.LT.0.9) THEN
          IQ = 1
          GO TO 40
      END IF
*
*       Delay Kozai cycle for EMAX > 0.9 and TC > 2000.
      IF (EMAX.GE.0.9) THEN
          CALL HICIRC(PMIN,EMAX,I,IG,BODYI,TG,TC,EC,EDT,W)
          IF (TC.GT.2000.OR.(ECC.LT.0.1.AND.EMAX.LT.0.95)) THEN
              IQ = 2
              GO TO 40
          END IF
      END IF
*
*       Check termination time for marginal stability criterion.
      IF (NAME(I).EQ.NAMEI.AND.(TIME+TOFF).GT.TCHECK) THEN
          WRITE (6,17)  NAME(I), ECC, EMAX, ALPH, ZFAC, ECC1, PMIN1,
     &                  PCRIT
   17     FORMAT (' ECCMOD UNSTAB    NM E EX IN YF E1 PM PC ',
     &                               I6,2F8.4,F7.1,F6.2,F7.3,1P,2E10.2)
          IQ = -4
          GO TO 40
      END IF
*
*       Ensure evaluation of EDT for large eccentricity and EDAV > 0.
      IF (ECC.GT.0.9.AND.EDAV.GT.0.0) THEN
          CALL HICIRC(QPERI,ECC,I,IG,BODYI,TG,TC,EC,EDT,W)
          IF (TC.LT.2000.0) THEN
              ICOLL = .TRUE.
          ELSE
              ICOLL = .FALSE.
              EDT = EDAV
          END IF
      ELSE
          ICOLL = .FALSE.
          EDT = EDAV
      END IF
*
*       Activate indicator on small QPERI to prevent collision before EMAX.
      IF (QPERI.LT.2.0*RM) ICOLL = .TRUE.
*
*       Check circularization time near turning point of normal binary.
      IF ((KSTARM(IM).GE.0.AND.EDAV.GT.0.0.AND.
     &   (DEDT.LT.0.0.OR.ECC.GT.ABS(EMAX-0.00005))).OR.ICOLL) THEN
          IF (.NOT.ICOLL) THEN
              CALL HICIRC(QPERI,ECC,I,IG,BODYI,TG,TC,EC,EDT,W)
          END IF
          IF (TC.LT.7.0D+04.AND.QPERI.LT.6.0*RM.AND.ITRY.LT.4.OR.
     &        TC.LT.1.0D+04.AND.QPERI.LT.5.0*RM.AND.ITRY.LT.8) THEN
              WRITE (6,18)  ECC, EMAX, EMIN, SEMI*SU, TG, TC, EDAV,
     &                      EDT, QPERI/RM
   18         FORMAT (' HICIRC TRY:    E EX EM A TG TC EDA EDT QP/R ',
     &                                 2F9.5,F7.3,F7.1,1P,5E9.1)
              ITRY = ITRY + 1
          END IF
          IF (ITRY.GT.8) ITRY = 0
*
*       Check termination for TC < 2000 to be consistent with routine TCIRC.
          IF (TC.GT.2000.0) GO TO 25
*
*       Evaluate chaos boundary parameters.
          EB = -0.5*CM(1,IM)*CM(2,IM)/a0
          CJ = CM(1,IM)*CM(2,IM)/BODY(I)*SQRT(BODY(I))
          ZJ = CJ*SQRT(QPERI*(1.0 + ECC))
          IC = NCHAOS + 1
          IDIS = 0
          ITRY = 0
*
*       Define oscillation period (dimensionless time).
          DO 20 K = 1,2
              IF (K.EQ.1) THEN
                  IK = I
              ELSE
                  IK = IG
              END IF
*       Specify polytropic index for each star (n = 3, 2 or 3/2).
              IF (KSTAR(IK).EQ.3.OR.KSTAR(IK).EQ.5) THEN
                  BODI = CM(K,IM)
                  CALL GIANT3(IK,BODI,WG,QG,ZN,QL)
                  W(K) = WG(1)
              ELSE
                  IP = 3
                  IF (KSTAR(IK).GE.3) IP = 2
                  IF (KSTAR(IK).EQ.4.OR.KSTAR(IK).EQ.6) IP = 3
                  IF (KSTAR(IK).EQ.0) IP = 1
                  W(K) = WW(IP)
              END IF
   20     CONTINUE
*
          CALL CHAOS0(QPERI,ECC,EB,ZJ,CM(1,IM),CM(2,IM),RADIUS(I),
     &                      RADIUS(IG),W,ECRIT(IC),AR(IC),BR(IC),IDIS)
*
*       Begin chaos/spiral stage if chaos boundary has been crossed.
          IF (IDIS.EQ.0) THEN
              WRITE (6,22)  TTOT, NAME(I), IC, ECC, EMAX, SEMI, QPERI,
     &                      EDAV, QPERI/RM
   22         FORMAT (' ECCMOD CHAOS    T NAM IC E EX A QP EDAV QP/R ',
     &                                  F9.2,I6,I4,2F9.5,1P,4E10.2)
              CALL FLUSH(6)
              IQ = -1
*       Activate chaos indicator for calling KSTIDE from RESET.
              KSTARM(IM) = -1
              GO TO 40
	  ELSE IF (IDIS.EQ.-1) THEN
              WRITE (6,23)  TTOT, NAME(I), IC, ECC, EMAX, SEMI, QPERI,
     &                      EDAV, QPERI/RM
   23         FORMAT (' ECCMOD SPIRAL    T NAM IC E EX A QP EDAV QP/R ',
     &                                   F9.2,I6,I4,2F9.5,1P,4E10.2)
              CALL FLUSH(6)
              IQ = -2
*       Note that KSTARM(IM) = -2 would be treated as existing spiral.
              KSTARM(IM) = -1
              GO TO 40
          ELSE IF (IDIS.EQ.1) THEN
*       Check collision condition to be sure.
              IF (QPERI.LT.RADIUS(I) + RADIUS(IG)) THEN
              WRITE (6,24)  TTOT, NAME(I), IC, ECC, SEMI, QPERI, EDAV,
     &                      QPERI/RM
   24         FORMAT (' ECCMOD DISRUPT    T NAM IC E A QP EDAV QP/R  ',
     &                                    F9.2,I6,I4,F8.4,1P,4E10.2)
                  IQ = -3
                  GO TO 40
              END IF
          END IF
      END IF
*
*       Form quantities for the outer orbit.
   25 RR = 0.0
      RV0 = 0.0
      V20 = 0.0
      DO 30 K = 1,3
          XR0(K) = X(K,JCOMP) - X(K,I)
          VR0(K) = XDOT(K,JCOMP) - XDOT(K,I)
          RV0 = RV0 + XR0(K)*VR0(K)
          V20 = V20 + VR0(K)**2
          RR = RR + XR0(K)**2
   30 CONTINUE
*
*       Construct the outer Runge-Lenz vector.
      DO 35 K = 1,3
          EI0(K) = (V20*XR0(K) - RV0*VR0(K))/ZMB - XR0(K)/SQRT(RR)
   35 CONTINUE
*
*       Specify total time interval (unperturbed or perturbed case).
      IF (LIST(1,I).EQ.0) THEN
          DT0 = TIME - T0(I)
*         DT = MIN(DT0,0.1*(1.0 - ECC)/(ABS(EDAV) + 1.0D-15))
      ELSE
*       Ensure that interval extends to next apocentre (including slow-down).
          IMOD = KSLOW(IPAIR)
          DT0 = 0.98*FLOAT(ISLOW(IMOD))*TK1
*         DT = MIN(DT0,0.1*(1.0 - ECC)*TG/TSTAR)
      END IF
*
*       Set partial interval and check remaining time.
      DT1 = MIN(DT0,0.1*TG*SQRT(1.0 - ECC**2)/TSTAR)
      IF (DTSUM + DT1.GT.DT0) DT1 = DT0 - DTSUM + 1.0D-12
*
*       Obtain quadrupole and tidal constants for each new binary.
      DTM = TIME - MAX(TEV0(I),TEV0(IG))
      IF (NAME(I).NE.NAMEI.OR.DTM.LE.DT0) THEN
          CALL QTIDES(I,IG,IM,SEMI,ECC)
          NAMEI = NAME(I)
          DEDT = 0.0
          DTPREV = 1.0
          NST = NSTAB(SEMI,SEMI1,ECC,ECC1,ZI,BODYI(1),BODYI(2),ZM3)
          IF (NST.EQ.0) THEN
              PCRIT = 0.98*QPERI*(1.0 - PERT)
              PCR = stability(BODYI(1),BODYI(2),ZM3,ECC,ECC1,ZI)
              PCR = PCR*SEMI
*       Specify reduced peri if old criterion < PCRIT/2 (avoids switching).
              IF (PCR.LT.0.5*PCRIT) THEN
                  PCRIT = 0.75*PCRIT
              END IF
              IF (PCRIT.LT.PCR.AND.PERT.LT.0.01.AND.
     &            ITIME.LT.20) THEN
                  ALPH = 360.0*ZI/TWOPI
                  FAIL = QPERI*(1-PERT) - PCR
                  WRITE (6,42)  TTOT, ALPH, ECC, ECC1, QPERI, FAIL, PERT
   42             FORMAT (' NEWSTAB    T INC EI EO QP FAIL PERT ',
     &                                 F7.1,F7.2,2F8.4,1P,3E10.2)
              END IF
          ELSE
              PCRIT = 1.01*QPERI
          END IF
          PMIN1 = SEMI1*(1.0 - ECC1)
          IF (PMIN1.LT.PCRIT) THEN
              NK = 1 + 10.0*ECC1/(1.0 - ECC1)
              TCHECK = TIME + TOFF + NK*TK1
          ELSE
              TCHECK = 1.0D+10
          END IF
      END IF
*
*       Determine time-step for integration from several criteria.
      ge=30*ECC+45*ECC**3+3.75*ECC**5
*       Apply factor of 2 correction for agreement with classical result.
      ge = 0.5*ge
      slr=SEMI*(1.0 - ECC**2)
      zq=ge/(slr**5*SEMI*sqrt(SEMI))
*       Correct for wrong eccentricity dependence.
      zq = zq*(1.0 - ECC**2)
      EDQ = zq*(cq(1) + cq(2))
      TQ = ECC/EDQ
      zgr = cgr/(slr*SEMI*SQRT(SEMI))
*       Note cgr is angular velocity (hence multiply by 2*pi; 12/03).
      TGR = TWOPI*ECC/zgr
*
*       Delay first time if apsidal motion or GR are important (ECC < 0.9).
      IF (ECC.LT.0.90.AND.ITIME.EQ.0) THEN
          IF (MIN(TQ,TGR)*TSTAR.LT.TG) THEN
              IQ = 2
              GO TO 40
          END IF
      END IF
*
*       Include check on tidal dissipation for large eccentricity.
      IF (KSTARM(IM).GE.0.AND.ECC.GT.0.9) THEN
          TT = ECC/ABS(EDT)
          TQ = MIN(TT,TQ)
*       Note that EDT from HICIRC agrees well with actual value from RKINT.
      END IF
*
*       Adopt harmonic mean of inner period and growth/quadrupole time.
      DT = 0.5*SQRT(TK*MIN(TG/TSTAR,TQ))
      IF (ECC.GT.0.98) DT = 0.1*DT
      DT = MIN(DT,0.001/(ABS(EDAV) + 1.0D-20))
*       Ensure conservative step first time and limit increase to factor 2.
      IF (ITIME.EQ.0) THEN
          DT = MIN(DT,1.0D-03*TG/TSTAR)
      ELSE
          DT = MIN(DT,2.0D0*DTPREV)
      END IF
*
*       Include safety check to avoid problems.
      ETRY = ECC + EDAV*DT1
      IF (ECC.GT.0.96.AND.ETRY.GT.0.98) THEN
          ETRY = MIN(ETRY,EMAX)
          QPTRY = SEMI*(1.0 - ETRY)
          CALL HICIRC(QPTRY,ETRY,I,IG,BODYI,TG,TC,EC,EDT,W)
          IF (TC.LT.10.0) THEN
              DT1 = 0.5*DT1
              DT = 0.5*DT
          END IF
      END IF
*
*       Ensure frequent termination check inside 5*RM to avoid overshooting.
      IF (QPERI.LT.5.0*RM.AND.EDAV.GT.0.0.AND.KSTARM(IM).GE.0) THEN
*       Restrict HIMOD interval so E + EDAV*DT1 gives new peri outside 4*RM.
          DTT = (QPERI - 4.0*RM)/(SEMI*EDAV)
          DT1 = MIN(DTT,DT1)
          DT1 = MAX(DT,DT1)
      END IF 
*
*       Modify the inner orbit using averaging of de/dt & dh/dt.
      DTPREV = DT
      CALL HIMOD(EVEC,A1,ZM3,ZMB,EI0,A2,ICALL,DT1,DT,ECC,XREL,VREL)
*
      IF (e0.ge.1.0) THEN
          IQ = 1
          GO TO 40
      END IF
      SEMI = a0
      ECC = e0
      ITIME = ITIME + 1
      DTSUM = DTSUM + DT1
      IF (DTSUM.LT.DT0) GO TO 4
*
*       Update step counter and next time for checking.
      NEINT = NEINT + ITIME
      TMDIS(IM) = TIME + DTSUM
*
*       Transform back to KS variables after integration.
   40 IF (ITIME.GT.0) THEN
          CALL PHYSKS(XREL,VREL,UI,V)
          DO 45 K = 1,4
              UM(K,IM) = UI(K)
              UMDOT(K,IM) = 0.25*V(K)
              IF (K.LT.4) THEN
                  YREL(K,IM) = XREL(K)
                  ZREL(K,IM) = VREL(K)
              END IF
   45     CONTINUE
          HM(IM) = -0.5D0*BODY(I)/a0
          CALL HIRECT(IM)
          ZMU = CM(1,IM)*CM(2,IM)/BODY(I)
          DECORR = ZMU*(HM0 - HM(IM))
          EMERGE = EMERGE - DECORR
          ECOLL = ECOLL + DECORR
          EGRAV = EGRAV + DECORR
*       Define circularization if e0 < 0.002 (end of Kozai cycle).
          IF (e0.LT.0.002) THEN
              KM = KSTARM(IM)
              KSTARM(IM) = 10
              IF (KM.LT.0) NCIRC = NCIRC + 1
              TMDIS(IM) = 1.0D+10
              TEV(N+IPAIR) = MIN(TEV(I),TEV(IG)) - 2.0*STEPX
              TB = DAYS*SEMI*SQRT(SEMI/BODY(I))
              ALPH = ZI*360.0/TWOPI
              WRITE (6,48)  NAME(I), NAME(IG), KM, e0, a0*SU, TB, ALPH
   48         FORMAT (' END KOZAI    NM K* E A P INC ',
     &                               2I6,I4,F7.3,F6.1,2F7.1)
          END IF
*       Include rectification of active SPIRAL (suppressed).
          IC = 0
          DO 50 K = 1,NCHAOS
              IF (NAMEC(K).EQ.NZERO - NAMEM(IM)) IC = K
*             IF (NAMEC(K).EQ.NZERO + NAME(I)) IC = K
   50     CONTINUE
          IF (IC.EQ.-3) THEN
              RP(IC) = a0*(1.0 - ECC)
              ES(IC) = ECC
              TOSC(IC) = TIME
              DH = HM0 - HM(IM)
              WRITE (6,60)  TTOT, TMDIS(IM), ECC, RP(IC), DH
   60         FORMAT (' UPDATE:    T TMD E RP DH ',
     &                             2F10.2,F8.4,1P,2E10.2)
          END IF
      END IF
*
      RETURN
*
      END
