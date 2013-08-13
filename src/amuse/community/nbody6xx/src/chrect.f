      SUBROUTINE CHRECT(IPAIR,DMR)
*     
*     
*     Rectification of chaotic orbit.
*     -------------------------------
*     
      INCLUDE 'common6.h'
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &     BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &     RP(NTMAX),ES(NTMAX),CM(2,NTMAX),IOSC(NTMAX),
     &     NAMEC(NTMAX)
      REAL*8  WW(6),W(4),WG(2),QG(2),WSCALE(2),QSCALE(2)
      CHARACTER*8  WHICH1
      LOGICAL  SLEEP
      DATA  WW  /2.119,3.113,8.175,3.742,4.953,9.413/
*     
*     
*     Define c.m. & KS indices and search current names for chaos index.
      I = N + IPAIR
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      IC = 0
      DO 1 K = 1,NCHAOS
         IF (NAMEC(K).EQ.NAME(I)) IC = K
    1 CONTINUE
*     
*     Include case of chain chaos without identified NAMEC.
      IF (IC.EQ.0.AND.NCHAOS.GT.0) THEN
         if(rank.eq.0)
     &        WRITE (6,2)  NCHAOS, IPAIR, KSTAR(I), NAME(I1), NAME(I2),
     &        LIST(1,I1), STEP(I1), STEP(I)
 2       FORMAT (' WARNING!    CHRECT    NCH KS K* NAM NP DT1 DTI ',
     &        3I4,2I6,I4,1P,2E10.2)
*     See whether wrong component name (+ NZERO) saved as NAMEC in CHAOS2.
         NAM2 = 0
         DO 6 L = 1,2
            NAM1 = KSAVE(2*L) - NAME(I2)
            IF (KSAVE(2*L-1).LT.0.AND.NAM1.EQ.NAME(I1)) THEN
               NAM2 = NZERO + NAME(I2)
            END IF
 6       CONTINUE
         IC = NCHAOS
         NAMC = NAMEC(IC)
*     Check identification for correct value (two K*=-2 are possible).
         DO 8 K = 1,NCHAOS
            IF (NAMEC(K).EQ.NAM2) IC = K
 8       CONTINUE
         NAMEC(IC) = NAME(I)
         IF (rank.eq.0.and.(NAM2.EQ.NAMC.OR.IC.LT.NCHAOS)) THEN
            WRITE (6,9)  IC, NCHAOS, NAM2, NAMC, NAME(I)
 9          FORMAT (' CHRECT RESTORE    IC NCH NM2 NMC NMI ',2I4,3I8)
         END IF
      ELSE IF (NCHAOS.EQ.0) THEN
         if(rank.eq.0)WRITE (6,3)  NCHAOS, IPAIR, KSTAR(I), NAME(I)
 3       FORMAT (' CHRECT RESTORE    NCH KS K* NAM ',3I4,I6)
*     Restore case of former merger with KSTARM < 0 to chaos table.
         NCHAOS = 1
         IC = 1
         NAMEC(NCHAOS) = NAME(I)
      END IF
*     
*     Save variables for diagnostic output.
      TIME0 = TOSC(IC)
      ES0 = ES(IC)
      TC = 0.0
      IDIS = KSTAR(I)
*     
*     Obtain current values of KS variables in case of spiral.
      IF (KSTAR(I).EQ.-2) THEN
*     Skip on call from RESET/MDOT with arbitrary phase (only updating).
         IF (ABS(TDOT2(IPAIR)).GT.1.0D-12.AND.DMR.GE.0.0) THEN
*     Rectify KS variables in order to obtain correct pericentre.
            CALL KSRECT(IPAIR)
*     Reduce eccentric anomaly by pi for inward motion.
            IF (TDOT2(IPAIR).LT.0.0D0) THEN
               CALL KSAPO(IPAIR)
            END IF
*     Transform from outward motion to exact pericentre.
            CALL KSPERI(IPAIR)
         END IF
*     
*     Form current two-body elements.
         SEMI = -0.5*BODY(I)/H(IPAIR)
         ECC2 = (1.0 - R(IPAIR)/SEMI)**2 +
     &        TDOT2(IPAIR)**2/(BODY(I)*SEMI)
         ECC = SQRT(ECC2)
*     
*     Update periastron and eccentricity (ECC modulation or mass loss).
         QPERI = SEMI*(1.0D0 - ECC)
         RP(IC) = QPERI
         ES(IC) = ECC
*     
*     Update orbital parameters after merger, mass loss or radius change.
         IF (DMR.GE.0.0D0) THEN
            KSTAR(I) = -KSTAR(I)
            CALL SPIRAL(IPAIR)
            IF (IPHASE.LT.0.OR.KSTAR(I).GT.0) GO TO 30
         END IF
*     
*     Check circularization time for spiral after radius/orbit expansion.
         ICIRC = -1
         CALL TCIRC(QPERI,ECC,I1,I2,ICIRC,TC)
         TC = MAX(TC,0.01D0)
*     Restrict lookup time to TC/2 for exit from possible merger.
         TEV(I1) = MIN(TEV(I1),TIME + 0.5*TC/TSTAR)
*     TEV(I2) = TEV(I1)
      END IF
*     
*     Form (new) semi-major axis, eccentricity & periastron distance.
      SEMI = -0.5*BODY(I)/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      ECC = SQRT(ECC2)
      QPERI = SEMI*(1.0D0 - ECC)
*     
*     Set latest values of periastron & eccentricity and update time.
      RP(IC) = QPERI
      ES(IC) = ECC
      TOSC(IC) = TIME
*     
*     Include check for SLEEP after recent WD/NS formation.
      SLEEP = .FALSE.
      KM = MAX(KSTAR(I1),KSTAR(I2))
      IF (KM.GE.10) THEN
*     Determine index of most recent degenerate object formation.
         IF (KSTAR(I2).LT.10) THEN
            J1 = I1
         ELSE IF (KSTAR(I1).LT.10) THEN
            J1 = I2
         ELSE
            J1 = I1
            IF (EPOCH(I2).GT.EPOCH(I1)) J1 = I2
         END IF
         IF (TIME - TEV0(J1).LT.2.0*STEPX.AND.TC.GT.3000.0) THEN
            SLEEP = .TRUE.
         END IF
      END IF
*     
*     Check circularization time after chain regularization.
      IF (IPHASE.EQ.8) THEN
         ICIRC = -1
         CALL TCIRC(QPERI,ECC,I1,I2,ICIRC,TC)
         IF (TC.GT.3000.0) THEN
            if(rank.eq.0)WRITE (6,4)  NAME(I1), TC
 4          FORMAT (' CHAIN SLEEP:    NM TC ',I6,1P,E10.2)
            SLEEP = .TRUE.
         END IF
      END IF
*     
*     Re-initialize all chaos parameters after expansion or check SLEEP.
      IF (DMR.GT.0.01.OR.KSTAR(I).EQ.-1.OR.SLEEP) THEN
         IF (rank.eq.0.and.DMR.GT.0.01.AND.TC.GT.100.0.AND.KM.LT.5) THEN
            WRITE (6,5)  NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2),
     &           KSTAR(I), TPHYS, RADIUS(I1), RADIUS(I2),
     &           QPERI, SEMI, ECC, ES0, BODY(I)*ZMBAR, TC
 5          FORMAT (' CHRECT:    NAM K* TP R* QP A E E0 M TC ',
     &           2I6,3I4,F8.1,1P,4E10.2,0P,3F7.3,F7.1)
         END IF
*     
*     Reset spiral indicator for degenerate component and long t_{circ}.
         IF (SLEEP) THEN
            KSTAR(I) = 0
            NSLP = NSLP + 1
            TK = SEMI*SQRT(SEMI/BODY(I))
            TB = YRS*TK
            XP = (TIME - TIME0)/(TWOPI*TK)
            QPS = SEMI*(1.0 - ECC)/MAX(RADIUS(I1),RADIUS(I2))
            if(rank.eq.0)
     &           WRITE (6,10)  TTOT, NAME(I1), NAME(I2), KSTAR(I1), 
     &           KSTAR(I2), ECC, ES0, QPS, SEMI, TC, TB, XP
 10         FORMAT (' SLEEP SPIRAL    T NM K* E E0 QP/S A TC TB DTK ',
     &           F9.2,2I6,2I4,2F8.4,1P,5E9.1)
*     Update chaos variables at end of routine SPIRAL (argument < 0).
            II = -I
            CALL SPIRAL(II)
            GO TO 30
         END IF
*     
*     Save eccentricity, binding energy & J0 and initialize EDEC & IOSC.
         ZMU = BODY(I1)*BODY(I2)/BODY(I)
         CJ = ZMU*SQRT(BODY(I))
         EB0(IC) = ZMU*H(IPAIR)
         ZJ0(IC) = CJ*SQRT(QPERI*(1.0 + ECC))
         EDEC(IC) = 0.0
         IOSC(IC) = 1
         ZN = 0.0
         DO 12 K = 1,2
            IK = I1 + K - 1
            IF (KSTAR(IK).EQ.3.OR.KSTAR(IK).EQ.5.OR.
     &           KSTAR(IK).EQ.6.OR.KSTAR(IK).EQ.9) THEN
               CALL GIANT(IPAIR,IK,WG,QG,WSCALE,QSCALE,ZN,QL)
               W(K) = WG(1)
            ELSE
               IP = 3
               IF (KSTAR(IK).EQ.0) IP = 1
               W(K) = WW(IP)
            END IF
 12      CONTINUE
*     
*     Set new chaos boundary parameters (ECRIT, AR & BR).
         CALL CHAOS0(QPERI,ECC,EB0(IC),ZJ0(IC),BODY(I1),BODY(I2),
     &        RADIUS(I1),RADIUS(I2),W,ECRIT(IC),AR(IC),BR(IC),IDIS)
*     
         RCOLL = RADIUS(I1) + RADIUS(I2)
         IF (IDIS.EQ.-1.AND.KSTAR(I).EQ.-1) THEN
            IOSC(IC) = 2
            if(rank.eq.0)
     &           WRITE (6,15)  TTOT, IPAIR, NAME(I1), NAME(I2), 
     &           KSTAR(I1), KSTAR(I2), RADIUS(I1), RADIUS(I2), QPERI,
     &           SEMI, ES0, ECC
 15         FORMAT (' CHAOS => SPIRAL    T KS NAM K* R* QP A E0 E ',
     &           F9.2,I4,2I6,2I4,1P,4E10.2,0P,2F7.3)
*     Activate spiral indicator and save time, pericentre & eccentricity.
            KSTAR(I) = -2
            TOSC(IC) = TIME
            RP(IC) = QPERI
            ES(IC) = ECC
            NSP = NSP + 1
            IF (KZ(8).GT.3) THEN
               CALL BINEV(IPAIR)
            END IF
            GO TO 30
         END IF
*     
*     Combine the two stars inelastically in case of chaos disruption.
         IF (IDIS.GT.0.AND.QPERI.LT.RCOLL) THEN
            R1 = MAX(RADIUS(I1),RADIUS(I2))
            WHICH1 = ' CHAOS  '
            IF (KSTAR(I).EQ.-2) WHICH1 = ' SPIRAL '
            if(rank.eq.0)
     &           WRITE (6,20)  WHICH1, IPAIR, NAME(I1), NAME(I2),
     &           KSTAR(I1), KSTAR(I2), R1, R(IPAIR), QPERI,
     &           SEMI, ECC, ES0, ZN
 20         FORMAT (' DISRUPTED',A8,'  KS NM K* R* R QP A E E0 n ',
     &           I4,2I6,2I4,1P,4E10.2,0P,3F7.3)
*     Update chaos variables at end of routine SPIRAL (argument < 0).
            IQCOLL = 1
            IF (KSTAR(I).EQ.-2) IQCOLL = 2
            II = -I
            CALL SPIRAL(II)
            KSTAR(I) = 0
            CALL XVPRED(I,0)
            KSPAIR = IPAIR
            CALL CMBODY(R(IPAIR),2)
            DMR = -1.0
            GO TO 30
         END IF
      ELSE IF (KSTAR(I).EQ.-2) THEN
         I1 = 2*IPAIR - 1
         I2 = I1 + 1
         IF (QPERI.LT.2.0*MAX(RADIUS(I1),RADIUS(I2))) THEN
*     Determine indices for primary & secondary star (donor & accretor).
            J1 = I1
            IF (RADIUS(I2).GT.RADIUS(I1)) J1 = I2
*     Define mass ratio and evaluate Roche radius for the primary.
            Q0 = BODY(J1)/(BODY(I) - BODY(J1))
            Q1 = Q0**0.3333
            Q2 = Q1**2
            RL1 = 0.49*Q2/(0.6*Q2 + LOG(1.0D0 + Q1))*SEMI
*     Check Roche radius but skip RESET call (no second EMERGE correction).
            IF (RADIUS(J1).GT.RL1.AND.IPHASE.NE.7) THEN
               if(rank.eq.0)
     &              WRITE (6,25) NAME(I1), NAME(I2), KSTAR(I1),
     &              KSTAR(I2), ECC, ES0, RCOLL, RL1, QPERI, SEMI, TC
 25            FORMAT (' DISRUPTED SPIRAL    NM K* E E0 RC RL QP A',
     &              ' TC ',2I6,2I4,2F7.3,1P,5E10.2)
*     Obtain KS variables at pericentre and enforce collision or CE.
               CALL KSPERI(IPAIR)
               KSTAR(I) = 0
               KSPAIR = IPAIR
               IQCOLL = 2
               CALL CMBODY(QPERI,2)
               DMR = -1.0
*     Note that same result achieved by TERM SPIRAL in case IPHASE = 7.
            END IF
         END IF
      END IF
*     
 30   RETURN
*     
      END
