      SUBROUTINE SPIRAL(IPAIR)
*     
*     
*     Tidal circularization of binary orbit.
*     --------------------------------------
*     
*     Rational function approximations of solution to the Hut
*     evolution equation with spin. Ref: A & A 99, 126, eqn (A15).
*     Developed by Rosemary Mardling (31/1/97 & 25/5/00). See 12/red/36.
*     
      INCLUDE 'common6.h'
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &     BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &     RP(NTMAX),ES(NTMAX),CM(2,NTMAX),IOSC(NTMAX),
     &     NAMEC(NTMAX)
      COMMON/SLOW0/  RANGE,ISLOW(10)
      common/spins/angmom0,rg2(2),m21,r21,semi0,C1,C2,C3,C4,semi
      REAL*8  WW(3),QQ(3),W(2),Q(2),AT0(2),M21,WG(2),QG(2),WSCALE(2),
     &     QSCALE(2),A(2),B(2),C(6),meanmotion,RJ(2),ROL(2)
      REAL*8  M0,MC1,MC2,CORERD
      EXTERNAL CORERD
      DATA  WW  /2.119,3.113,8.175/
      DATA  QQ  /0.4909,0.4219,0.2372/
      DATA  A  /6.306505,-7.297806/
      DATA  B  /32.17211,13.01598/
      DATA  C  /5.101417,24.71539,-9.627739,1.733964,
     &     -2.314374,-4.127795/
      DATA  ECCM  /0.002/
      SAVE  IONE,IWARN
      DATA  IONE,IWARN /0,0/
*     
*     
*     Check for just table updating at KS termination.
      IREM = 0
      IF (IPAIR.LT.0) THEN
         I = -IPAIR
         DO 1 K = 1,NCHAOS
            IF (NAMEC(K).EQ.NAME(I)) IREM = K
 1       CONTINUE
         IONE = 0
         IF (IREM.GT.0) GO TO 30
*     Include case of removing a specified chaos index (I <= NCHAOS).
         IF (I.LE.NCHAOS) THEN
            IREM = I
            GO TO 30
         END IF
         GO TO 100
      END IF
*     
*     Define c.m. & KS indices and search current names for chaos index.
      I = N + IPAIR
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      NSPIR = NSPIR + 1
      IC = 0
      DO 2 K = 1,NCHAOS
         IF (NAMEC(K).EQ.NAME(I)) IC = K
    2 CONTINUE
*     
*     Set non-zero index on call from CHRECT (KSTAR = -KSTAR) for CE skip.
      ISTAR = 0
      IF (KSTAR(I).GT.0) THEN
         KSTAR(I) = -KSTAR(I)
         ISTAR = 1
      END IF
*     
*     Include case of chain chaos without identified NAMEC.
      IF (IC.EQ.0.AND.NCHAOS.GT.0) THEN
         if(rank.eq.0)
     &        WRITE (6,3)  NCHAOS, KSTAR(I), IPAIR, NAME(I1), NAME(I2),
     &        STEP(2*IPAIR-1), STEP(I)
 3       FORMAT (' WARNING!    SPIRAL    NCH K* KS NAME STEP1 STEPI',
     &        3I4,2I6,1P,2E10.2)
*     Search for component names (set in CHAOS2).
         DO 5 K = 1,NCHAOS
            IF (NAMEC(K).EQ.NAME(I1).OR.NAMEC(K).EQ.NAME(I2)) IC = K
 5       CONTINUE
*     
*     Include safety check just in case.
         IF (rank.eq.0.and.IC.EQ.0) THEN
            WRITE (6,6)  NAME(I), KSTAR(I), (NAMEC(K),K=1,NCHAOS)
 6          FORMAT (' DANGER!    SPIRAL    NAMI K* NAMEC ',15I6)
            STOP
         END IF
         NAMEC(IC) = NAME(I)
      END IF
*     
*     Copy spiral parameters and set semi-major axis & period.
      TIME0 = TOSC(IC)
      RP0 = RP(IC)
      ES0 = ES(IC)
      HI = H(IPAIR)
      SEMI0 = -0.5*BODY(I)/H(IPAIR)
      TK = SEMI0*SQRT(SEMI0/BODY(I))
*     
*     Use actual elements for perturbed binary.
      IF (LIST(1,I1).GT.0) THEN
         ECC2 = (1.0 - R(IPAIR)/SEMI0)**2 +
     &        TDOT2(IPAIR)**2/(BODY(I)*SEMI0)
         ECC0 = SQRT(ECC2)
         ES0 = ECC0
         RP0 = SEMI0*(1.0 - ECC0)
         RP(IC) = RP0
      END IF
*     
*     Specify index J1 as biggest radius to be used with AT0(1).
      IF (RADIUS(I1).GE.RADIUS(I2)) THEN
         J1 = I1
         J2 = I2
      ELSE
         J1 = I2
         J2 = I1
      END IF
*     
*     Define oscillation period (dimensionless time) and damping constants.
      ZN = 0.0
      DO 10 K = 1,2
         IK = I1 + K - 1
         IF (K.EQ.1) THEN
            IK = J1
         ELSE
            IK = J2
         END IF
*     Specify polytropic index for each star (n = 3, 2 or 3/2).
         IF (KSTAR(IK).EQ.3.OR.KSTAR(IK).EQ.5.OR.
     &        KSTAR(IK).EQ.6.OR.KSTAR(IK).EQ.9) THEN
            CALL GIANT(IPAIR,IK,WG,QG,WSCALE,QSCALE,ZN,QL)
            W(K) = WG(1)
            Q(K) = QG(1)
*     Note: rg2 should really be given by k2 in HRDIAG. 
            rg2(k)= 0.1*(1.0 - CM(K,IC)/BODY(IK))
         ELSE
            QL = 1.0E+04
            IP = 3
            IF (KSTAR(IK).GE.3) IP = 2
            IF (KSTAR(IK).EQ.4.OR.KSTAR(IK).EQ.7) IP = 3
            IF (KSTAR(IK).EQ.8) IP = 3
            IF (KSTAR(IK).EQ.0) IP = 1
            W(K) = WW(IP)
            Q(K) = QQ(IP)
            rg2(k)= 0.21
            IF (KSTAR(IK).LE.2.OR.KSTAR(IK).EQ.7) rg2(k) = 0.1
            IF (KSTAR(IK).EQ.4) rg2(k)= 0.1*(1.0 - CM(K,IC)/BODY(IK))
         END IF
         TL = TWOPI*RADIUS(IK)*SQRT(RADIUS(IK)/BODY(IK)/W(K))
         AT0(K) = 1.0/(QL*TL)
 10   CONTINUE
*     
*     Form mass, radius & pericentre ratio.
      IF (RADIUS(I1).GE.RADIUS(I2)) THEN
         M21 = BODY(I2)/BODY(I1)
         R21 = RADIUS(I2)/RADIUS(I1)
         RP1 = RP(IC)/RADIUS(I1)
         rad = radius(i1)
      ELSE
         M21 = BODY(I1)/BODY(I2)
         R21 = RADIUS(I1)/RADIUS(I2)
         RP1 = RP(IC)/RADIUS(I2)
         rad = radius(i2)
      END IF
*     
*     Define initial angular momentum from the scaled semi-major axis.
      semi0 = RP1/(1.0 - ES0)
*     
*     Form the initial mean motion in N-body units.
      meanmotion = sqrt((body(i1)+body(i2))/(rad*semi0)**3)
*     
*     Convert from angular momentum to omega (denoted spin1 & spin2).
      IF (KSTAR(J1).LE.2.OR.(KSTAR(J1).GE.7.AND.KSTAR(J1).NE.9)) THEN
         SPIN1 = SPIN(J1)/(rg2(1)*BODY(J1)*RADIUS(J1)**2)
      ELSE 
         KW = KSTAR(J1)
         M0 = BODY0(J1)*SMU
         MC1 = CM(1,IC)*SMU
         IF (MC1.LE.0.0D0.OR.MC1.GT.M0) THEN
            MC1 = 0.3 + 0.1*FLOAT(KW - 3)
            IF(KW.EQ.9) MC1 = MIN(0.3D0,0.95*M0)
            CM(1,IC) = MC1/ZMBAR
         END IF
         ZDUM = 2.0D0
         RC1 = CORERD(KW,MC1,M0,ZDUM)/SU
         SPIN1 = SPIN(J1)/(rg2(1)*BODY(J1)*RADIUS(J1)**2 +
     &        0.21*MC1/SMU*RC1**2)
      END IF
      IF (KSTAR(J2).LE.2.OR.(KSTAR(J2).GE.7.AND.KSTAR(J2).NE.9)) THEN
         spin2 = SPIN(J2)/(rg2(2)*BODY(J2)*RADIUS(J2)**2)
*     Produce diagnostic information for recent WD.
         IF (rank.eq.0.and.KSTAR(J2).GE.10.AND.KSTAR(J2).LE.12.AND.
     &        TIME.LT.TEV0(J2) + 0.01) THEN
            WRITE (95,780)  TTOT, NAME(J2), NAME(J1), KSTAR(J1),
     &           spin1, spin2, meanmotion, SPIN(J2)
 780        FORMAT (' WD SPIN    T NM K*1 ROT1 ROT2 <n> S2 ',
     &           F8.1,2I6,I4,1P,4E10.2)
            CALL FLUSH(95)
         END IF
      ELSE
         KW = KSTAR(J2)
         M0 = BODY0(J2)*SMU
         MC2 = CM(2,IC)*SMU
         IF (MC2.LE.1.0D-10.OR.MC2.GT.M0) THEN
            MC2 = 0.3 + 0.1*FLOAT(KW - 3)
            IF(KW.EQ.9) MC2 = MIN(0.3D0,0.95*M0)
            CM(2,IC) = MC2/ZMBAR
         END IF
         ZDUM = 2.0D0
         RC2 = CORERD(KW,MC2,M0,ZDUM)/SU
         spin2 = SPIN(J2)/(rg2(2)*BODY(J2)*RADIUS(J2)**2 +
     &        0.21*MC2/SMU*RC2**2)
      END IF
*     
      IF (IONE.EQ.0) THEN
         IONE = IONE + 1
         IF (ABS(SPIN1).GT.0.D0.AND.ABS(SPIN2).GT.0.D0) THEN
            TS = 1.0D+06*365.0*TWOPI*TSTAR
            if(rank.eq.0)
     &           WRITE (6,790)  TS/spin1, TS/spin2, TS/meanmotion,
     &           DAYS*TK, RADIUS(J1)*SU, RADIUS(J2)*SU
 790        FORMAT (' BEGIN    TROT P TK R* ',1P,2E9.1,0P,4F9.2)
         END IF
      END IF
*     
*     Scale the spins by mean motion and define angular momentum.
      spin10=spin1/meanmotion
      spin20=spin2/meanmotion
      angmom0=(m21/(1+m21))*semi0**2*sqrt(1-es0**2)+rg2(1)*spin10+
     &     m21*r21**2*rg2(2)*spin20
*     
*     Evaluate damping coefficients (Mardling & SJA, M.N. 321, 398, 2001).
      cf = 54.0*twopi/5.0
      C1 = cf*(AT0(1)*(Q(1)/W(1))**2*(1.0 + M21)*M21)/semi0**8
      C2 = cf*(AT0(2)*(Q(2)/W(2))**2*((1.0 + M21)/M21**2)*R21**8)/
     &     semi0**8
      C3 = (cf/9.0)*(AT0(1)*(Q(1)/W(1))**2*M21**2)/rg2(1)/semi0**6
      C4 = (cf/9.0)*(AT0(2)*(Q(2)/W(2))**2/M21**2)*R21**6/rg2(2)/
     &     semi0**6
*     
*     Adopt WD scaling for any NS or BH to avoid numerical problem.
      IF (KSTAR(I1).GE.13) THEN
         C1 = 1.0D-04*C1
      END IF
      IF (KSTAR(I2).GE.13) THEN
         C2 = 1.0D-04*C2
      END IF
*     
      IF (time - time0.LE.0.0D0) GO TO 100
*     
*     Obtain dominant terms of the derivatives (cf. routine HUT/DERIV2).
      fac=1-es0**2
      udot1=-(es0/fac**6.5)*(C1 + C2)
      udot2=(1.0/fac)**6*C3
      udot3=(1.0/fac)**6*C4
*     
*     Choose the step from smallest time-scale (! time0 < time possible).
      taux = min(abs(1.0/udot1),abs(1.0/udot2),abs(1.0/udot3))
      nstep = 1 + 100.0*sqrt(ABS(time - time0)/taux)
      nstep = min(nstep,100)
*     Include extra steps for long intervals and/or AGB type 5 or 6.
      IF (TIME - TIME0.GT.0.1) nstep = nstep + 20
      IF (KSTAR(J1).EQ.5.OR.KSTAR(J1).EQ.6) nstep = nstep + 10
*     Evaluate the equilibrium angular velocity.
      e2 = es0**2
      f2 = ((0.3125*e2 + 5.625)*e2 + 7.5)*e2 + 1.0
      f5 = (0.375*e2 + 3.0)*e2 + 1.0
      omeq = f2*meanmotion/(f5*fac**1.5)
*     Increase number of steps on slow primary rotation.
      IF (omeq - spin1.GT.0.2*omeq) nstep = nstep + 10
      dtau1 = ABS(time - time0)/float(nstep)
*     
*     Check circularization time after merger or large interval.
      IF (IPHASE.EQ.7.OR.TIME - TIME0.GT.1.0) THEN
         ICIRC = -1
         A0 = -0.5*BODY(I)/H(IPAIR)
         QPERI = A0*(1.0 - es0)
         CALL TCIRC(QPERI,es0,I1,I2,ICIRC,TC)
         DT = MIN(0.1D0*TC,1.0D0)
*     Ensure careful integration with DT=0.1*TC for rapid evolution.
         IF (TC.LT.10.0.AND.(KSTAR(J1).GT.1.OR.ECC.GT.0.1)) THEN
            nstep = 500*(1.0 + 100.0/TC)
            nstep = MIN(nstep,5000)
            if(rank.eq.0)
     &           WRITE (6,12)  NAME(J1), nstep, ES0, QPERI,
     &           TIME-TIME0, TC
 12         FORMAT (' SPIRAL RESET    NM # E QP T-T0 TC ',
     &           2I6,F7.3,1P,3E10.2)
         END IF
         dtau1 = ABS(DT)/float(nstep)
      END IF
*     
      IONE = IONE + 1
      IF (rank.eq.0.and.IONE.LE.2.OR.MOD(IONE,10000).EQ.0) THEN
         WRITE (6,13) LIST(1,I1), NAME(I1),es0, omeq, meanmotion,
     &        spin1, spin2
 13      FORMAT (' SPIRAL    NP NM es0 omeq <n> spin ',
     &        I4,I6,F9.5,1P,4E10.2)
      END IF
*     
*     Integrate equations for eccentricity and angular velocities.
      call hut(es0,spin10,spin20,ecc,spin1,spin2,nstep,dtau1)
*     
*     Ensure that no overshooting takes place.
      IF (spin10.LT.1.0.and.spin1.gt.1.0) THEN
         IF (rank.eq.0.and.MOD(IWARN,100).EQ.0) THEN
            WRITE (66,14) TTOT, NAME(J1), KSTAR(J1), nstep, ecc,
     &           RP0*SU/(1.0-es0), spin1-1.0
 14         FORMAT (' SPIN WARNING   T NM K* # E A S-1 ',
     &           F8.1,I7,I4,I6,F7.3,1P,2E10.1)
            CALL FLUSH(66)
         END IF
         spin1 = 1.0
         IWARN = IWARN + 1
      END IF
      IF (spin20.LT.1.0.and.spin2.gt.1.0) THEN
         IF (rank.eq.0.and.MOD(IWARN,100).EQ.0) THEN
            WRITE (66,14) TTOT, NAME(J2), KSTAR(J2), nstep, ecc,
     &           RP0*SU/(1.0-es0), spin2-1.0
            CALL FLUSH(66)
         END IF
         spin2 = 1.0
         IWARN = IWARN + 1
      END IF
*     
*     Re-scale the semi-major axis and angular velocities to N-body units.
      semi = rad*semi0*semi
      spin1 = meanmotion*spin1
      spin2 = meanmotion*spin2
      ecc = MAX(ecc,0.001d0)
      ecc = MIN(ecc,0.999d0)
*     
*     Convert back to angular momenta.
      IF (KSTAR(J1).LE.2.OR.(KSTAR(J1).GE.7.AND.KSTAR(J1).NE.9)) THEN
         SPIN(J1) = rg2(1)*BODY(J1)*RADIUS(J1)**2*spin1
      ELSE
         SPIN(J1) = (rg2(1)*BODY(J1)*RADIUS(J1)**2 +
     &        0.21*MC1/SMU*RC1**2)*spin1
      END IF
      IF (KSTAR(J2).LE.2.OR.(KSTAR(J2).GE.7.AND.KSTAR(J2).NE.9)) THEN
         SPIN(J2) = rg2(2)*BODY(J2)*RADIUS(J2)**2*spin2
      ELSE
         SPIN(J2) = (rg2(2)*BODY(J2)*RADIUS(J2)**2 +
     &        0.21*MC2/SMU*RC2**2)*spin2
      END IF
*     
*     Obtain the tidal contributions from integration.
      A0 = RP0/(1.0 - ES0)
      DH = -0.5*BODY(I)*(1.0/SEMI - 1.0/A0)
*     
*     Update energy and semi-major axis.
      H(IPAIR) = H(IPAIR) + DH
      SEMI = -0.5*BODY(I)/H(IPAIR)
*     
*     Set new pericentre from final elements and update reference values.
      IF (LIST(1,I1).EQ.0) THEN
         QPERI = semi*(1.0 - ecc)
      ELSE
         QPERI = RP0*(1.0 + ES0)/(1.0 + ECC)
      END IF
      RP(IC) = QPERI
      ES(IC) = ECC
      TOSC(IC) = TIME
*     
*     Define synchronous state at the end and adopt ECC = ECCM.
      IF (ECC.LT.ECCM) THEN
         IONE = 0
         KSTAR(I) = 10
         NCIRC = NCIRC + 1
         TB = DAYS*SEMI*SQRT(SEMI/BODY(I))
*     Ensure the dominant star has synchronous rotation.
         IF(KSTAR(J1).LE.2.OR.(KSTAR(J1).GE.7.AND.KSTAR(J1).NE.9))THEN
            SPIN(J1) = rg2(1)*BODY(J1)*RADIUS(J1)**2*meanmotion
         ELSE
            SPIN(J1) = (rg2(1)*BODY(J1)*RADIUS(J1)**2 +
     &           0.21*MC1/SMU*RC1**2)*meanmotion
         END IF
         DOM = (spin1 - meanmotion)/meanmotion
         if(rank.eq.0)
     &        WRITE (6,15)  TTOT, IPAIR, ES0, ECC, 
     &        RP0, R(IPAIR), H(IPAIR), TB, ZN, DOM
 15      FORMAT (' END SPIRAL    T KS E0 E RP0 R H P n DOM/OM ',
     &        F9.2,I4,2F8.4,1P,4E10.2,0P,F5.1,F7.3)
         ECC = ECCM
      END IF
*     
*     Ensure apocentre or arbitrary phase is replaced by pericentre.
      ERR = ABS(RP0 - R(IPAIR))/RP0
      IF (R(IPAIR).GT.SEMI.OR.ERR.GT.1.0E-05) THEN
*     Reduce eccentric anomaly by pi for inward motion.
         IF (TDOT2(IPAIR).LT.0.0D0) THEN
            CALL KSAPO(IPAIR)
         END IF
*     Transform from outward motion (anomaly < pi) to exact pericentre.
         IF (GAMMA(IPAIR).LT.1.0D-04) THEN
*     Note consistency problem for large GAMMA (HDOT effect not included).
            TT0 = TIME
            CALL KSPERI(IPAIR)
*     Restore TIME just in case (possible bug in scheduling of #I1).
            TIME = TT0
         END IF
      END IF
*     
*     Form KS coordinate scaling factor from pericentre ratio.
      C1 = SQRT(QPERI/RP0)
*     
*     Set KS velocity scaling from energy relation with RP0 instead of R.
      C2 = SQRT((BODY(I) + H(IPAIR)*QPERI)/(BODY(I) + HI*RP0))
*     
*     Scale KS variables to yield the prescribed elements (set TDOT2 >= 0).
      R(IPAIR) = 0.0D0
      TDOT2(IPAIR) = 0.0
      DO 18 K = 1,4
         U(K,IPAIR) = C1*U(K,IPAIR)
         UDOT(K,IPAIR) = C2*UDOT(K,IPAIR)
         U0(K,IPAIR) = U(K,IPAIR)
         R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
         TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0*U(K,IPAIR)*UDOT(K,IPAIR)
 18   CONTINUE
      TDOT2(IPAIR) = MAX(TDOT2(IPAIR),0.0D0)
*     
*     Perform energy correction to maintain conservation.
      ZMU = BODY(I1)*BODY(I2)/BODY(I)
      ECOLL = ECOLL + ZMU*(HI - H(IPAIR))
      EGRAV = EGRAV + ZMU*(HI - H(IPAIR))
*     
*     Rectify large eccentricity deviation from integrated value.
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      ECCF = SQRT(ECC2)
      IF (ABS(ECC - ECCF).GT.0.01) THEN
         CALL KSRECT(IPAIR)
      END IF
*     
*     Check Roche time and update RADIUS of primary frequently.
      HI = H(IPAIR)
*     Specify semi-major axis and binding energy for evaluating Roche time.
      AC = QPERI*(1.0 + ECC)
      H(IPAIR) = -0.5*BODY(I)/AC
      CALL TRFLOW(IPAIR,DTR)
      H(IPAIR) = HI
*     
*     Obtain Roche radius and stellar radius for each component.
      J = I1
      DO 20 K = 1,2
         Q1 = BODY(J)/(BODY(I) - BODY(J))
         ROL(K) = RL(Q1)
         RJ(K) = RADIUS(J)
         J = I2
 20   CONTINUE
*     
*     Determine indices for primary & secondary star (donor & accretor).
      IF (RJ(1)/ROL(1).GE.RJ(2)/ROL(2)) THEN
         J1 = I1
         J2 = I2
      ELSE
         J1 = I2
         J2 = I1
      END IF
*     
*     Form critical mass ratio using approximate core mass.
      IF(KSTAR(J1).EQ.2)THEN
         QC = 4.D0
      ELSE IF(KSTAR(J1).GE.3.AND.KSTAR(J1).LE.6)THEN
         ZPARS7 = 0.3
         ZMC = CM(1,IC)
         IF (J1.EQ.I2) ZMC = CM(2,IC)
         QC = (1.67D0-ZPARS7+2.D0*(ZMC/BODY(J1))**5)/2.13D0
*     
*     Consider using condition of Hjellming & Webbink, 1987, ApJ, 318, 794.
*     QC = 0.362D0 + 1.D0/(3.D0*(1.D0 - MASSC(1)/MASS(1)))
      ELSE IF(KSTAR(J1).EQ.8.OR.KSTAR(J1).EQ.9)THEN
         QC = 0.784D0
      ELSE
         QC = 1000.0
      ENDIF
*     
*     Adopt common envelope evolution for eccentric binaries with Q1 > QC.
      Q1 = BODY(J1)/BODY(J2)
      IF (DTR.LT.0.1/TSTAR.AND.Q1.GT.QC.AND.ISTAR.EQ.0) THEN
         KSPAIR = IPAIR
         TIME = TBLOCK
         IQCOLL = 1
         CALL EXPEL(J1,J2,ICASE)
*     Exit before updating TEV in the unlikely case of coalescence.
         IF (IPHASE.EQ.-1) GO TO 100
         TEV(I1) = TEV(I) + 2.0*STEPX
         TEV(I2) = TEV(I1)
         GO TO 100
      END IF
*     
*     Enforce termination at constant angular momentum on Roche condition.
      IF (DTR.LE.0.1/TSTAR) THEN
         AF = SEMI*(1.0 - ECC**2)/(1.0 - ECCM**2)
         H(IPAIR) = -0.5*BODY(I)/AF
         ECOLL = ECOLL + ZMU*(HI - H(IPAIR))
         EGRAV = EGRAV + ZMU*(HI - H(IPAIR))
         KSTAR(I) = 10
         if(rank.eq.0)WRITE (6,21)  ES0, ECC, AF/SEMI, DTR, SEMI, AF
 21      FORMAT (' WARNING!    SPIRAL TERM    E0 E A/A0 DTR A AF ',
     &        2F8.4,F6.2,1P,3E10.2)
*     Ensure the dominant star has synchronous rotation.
         meanmotion = SQRT(BODY(I)/AF**3)
         IF(KSTAR(J1).LE.2.OR.(KSTAR(J1).GE.7.AND.KSTAR(J1).NE.9))THEN
            SPIN(J1) = rg2(1)*BODY(J1)*RADIUS(J1)**2*meanmotion
         ELSE
            SPIN(J1) = (rg2(1)*BODY(J1)*RADIUS(J1)**2 +
     &           0.21*MC1/SMU*RC1**2)*meanmotion
         END IF
         DOM = (spin1 - meanmotion)/meanmotion
         if(rank.eq.0)WRITE (6,26)  NAME(J1), Q1, QC, DOM
 26      FORMAT (' ENFORCED ROCHE   NM(J1) Q1 QC DOM/OM',I8,2F8.2,F7.3)
*     Modify KS variables to yield circular velocity with energy H.
         CALL EXPAND(IPAIR,R(IPAIR))
         CALL RESOLV(IPAIR,1)
      END IF
*     
*     Initialize KS for perturbed motion (cf. backwards step in KSPERI).
      IF (LIST(1,I1).GT.0) THEN
*     Ensure unperturbed motion for small GAMMA.
         IF (GAMMA(IPAIR).LT.1.0D-10.AND.ES0.LT.0.95) THEN
            LIST(1,I1) = 0
         END IF
         IF (GAMMA(IPAIR).LT.5.0D-07.AND.ES0.LT.0.3) THEN
            LIST(1,I1) = 0
         END IF
         CALL RESOLV(IPAIR,1)
         IMOD = KSLOW(IPAIR)
         CALL KSPOLY(IPAIR,IMOD)
      END IF
*     
*     Rectify the orbit to yield consistent variables (only at end).
      IF (KSTAR(I).EQ.10) THEN
         CALL KSRECT(IPAIR)
*     
*     Re-evaluate eccentricity after rectification.
         SEMI = -0.5*BODY(I)/H(IPAIR)
         ECC2 = (1.0-R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
         ECC = SQRT(ECC2)
*     
*     Deform the orbit to small eccentricity (H = const).
         IF (ECC.GT.ECCM) THEN
            CALL DEFORM(IPAIR,ECC,ECCM)
*     
            if(rank.eq.0)WRITE (6,22)  ECC, ECCM, TIME-TIME0, SEMI,
     &           BODY(I1)*ZMBAR, BODY(I2)*ZMBAR
 22         FORMAT (' DEFORM SPIRAL    E EF T-TOSC SEMI M1 M2 ',
     &           2F8.4,F9.4,1P,E10.2,0P,2F6.2)
         END IF
*     
*     Determine indices for primary & secondary star (donor & accretor).
         IF (BODY(I1)/RADIUS(I1)**3.LT.BODY(I2)/RADIUS(I2)**3) THEN
            J1 = I1
            J2 = I2
         ELSE
            J2 = I1
            J1 = I2
         END IF
*     
*     Define mass ratio and evaluate Roche radius for the primary.
         Q0 = BODY(J1)/BODY(J2)
         Q1 = Q0**0.3333
         Q2 = Q1**2
         RL1 = 0.49*Q2/(0.6*Q2 + LOG(1.0D0 + Q1))*SEMI
*     
*     Update Roche look-up time (does not depend on choice of primary).
         CALL TRFLOW(IPAIR,DTR)
         TEV(I) = TIME + DTR
         TK = DAYS*SEMI*SQRT(SEMI/BODY(I))
*     
         if(rank.eq.0)
     &        WRITE (6,24)  IC, NAME(J1), NAME(J2), KSTAR(J1), 
     &        KSTAR(J2), SEMI*SU, RL1*SU, RADIUS(J1)*SU, RADIUS(J2)*SU,
     &        TPHYS, TK, DTR
 24      FORMAT (' ROCHE CHECK    IC NAM K* A RL R* TP TK DTR ',
     &        I4,2I6,2I4,4F7.1,F9.2,1P,2E9.1)
         ZM1 = BODY(J1)*SMU
         CALL TRDOT(J1,DTM,ZM1)
         TEV(J1) = TIME + DTM
*     IF (DTR.LT.0.1/TSTAR) THEN
*     TEV(I1) = TEV(I) + 2.0*STEPX
*     TEV(I2) = TEV(I1)
*     END IF
*     Include enforcement of Roche coalescence to prevent shrinkage stop.
         IF (DTR.EQ.0.0D0.AND.SEMI.LT.0.5*RADIUS(J1)) THEN
            KSPAIR = IPAIR
            IQCOLL = 1
            if(rank.eq.0)WRITE (6,23)
 23         FORMAT (' ENFORCED COAL')
            CALL CMBODY(R(IPAIR),2)
            IF (IPHASE.LT.0) GO TO 100
         END IF
      END IF
*     
*     Include occasional diagnostics of spiral evolution (every output).
      IF (ABS(TIME - TADJ).LT.STEP(I1)) THEN
         IF (RP1.LT.3.0.AND.(ECC.LT.0.1.OR.KSTAR(J1).GE.3)) THEN
            NP = LIST(1,I1)
            TK = DAYS*SEMI*SQRT(SEMI/BODY(I))
            if(rank.eq.0)
     &           WRITE (6,25)  NAME(I1), NAME(I2), IPAIR, NCHAOS,
     &           KSTAR(I1), KSTAR(I2), NP, ECC, RP1, TK
 25         FORMAT (' SPIRAL    NAM KS NCH K* NP E RP P ',
     &           2I6,I5,4I4,F8.4,F6.2,F8.2)
         END IF
      END IF
*     
*     Include enforced circularization for difficult cases.
      IF ((ECC.LT.0.6.AND.GAMMA(IPAIR).LT.2.0D-07).OR.
     &     (ECC.LT.0.8.AND.GAMMA(IPAIR).LT.1.0D-08).OR.
     &     (ECC.LT.0.95.AND.GAMMA(IPAIR).LT.1.0D-09)) THEN
         ICIRC = 0
         CALL TCIRC(QPERI,ECC,I1,I2,ICIRC,TC) 
         IF (TC.GT.200.0.OR.(TC.GT.10.AND.ECC.LT.0.7)) THEN
            CALL KSPERI(IPAIR)
            CALL KSAPO(IPAIR)
            CALL DEFORM(IPAIR,ECC,ECCM)
            KSTAR(I) = 10
            NCIRC = NCIRC + 1
            CALL RESOLV(IPAIR,1)
            CALL KSPOLY(IPAIR,1)
            NP = LIST(1,I1)
            if(rank.eq.0)
     &           WRITE (6,28)  IPAIR, NP, NAME(I1),
     &           ECC, TC, GAMMA(IPAIR)
 28         FORMAT (' ENFORCED CIRC   KS NP NM E TC G  ',
     &           2I4,I6,F7.3,1P,2E9.1)
         END IF
      END IF
*     
*     Set standard binary on large ECC and TC if small GAMMA (large omeq).
      IF (ECC.GT.0.95.AND.GAMMA(IPAIR).LT.1.0D-09) THEN
         ICIRC = 0
         CALL TCIRC(QPERI,ECC,I1,I2,ICIRC,TC) 
         IF (TC.GT.1.0D+04) THEN
            KSTAR(I) = 0
            IREM = IC
            if(rank.eq.0)WRITE (6,29)  NAME(I1), LIST(1,I1), ECC, TC
 29         FORMAT (' FROZEN CIRC    NM NP E TC ',I7,I4,F9.5,1P,E9.1)
            GO TO 30
         END IF
      END IF
*     
*     Check for terminated or escaped chaotic binaries at end of spiral.
 30   IF (KSTAR(I).EQ.10.OR.IREM.GT.0) THEN
         IF (IREM.GT.0) GO TO 50
         J = 1
*     See whether case #J indicates current or escaped/disrupted KS binary.
 40      DO 45 JPAIR = 1,NPAIRS
            IF (NAMEC(J).EQ.NAME(N+JPAIR)) THEN
*     Update #J if KSTAR > 0, otherwise consider next member.
               IF (KSTAR(N+JPAIR).GT.0) THEN
                  GO TO 50
               ELSE
                  GO TO 70
               END IF
            END IF
 45      CONTINUE
*     
*     Skip during multiple regularizations (KSTAR not visible).
         IF (NSUB.GT.0) THEN
            GO TO 70
         END IF
*     
*     Skip removal if chaos binary is member of single/double merger.
         IF (NMERGE.GT.0) THEN
            DO 48 JPAIR = 1,NPAIRS
               IF (NAME(N+JPAIR).LT.0) THEN
                  IF (NAMEC(J).EQ.NZERO + NAME(2*JPAIR-1).OR.
     &                 NAMEC(J).EQ.NZERO + NAME(2*JPAIR).OR.
     &                 NAMEC(J).LT.-2*NZERO.OR.
     &                 NAMEC(J).EQ.NAME(2*JPAIR)) THEN
                     if(rank.eq.0)
     &                    WRITE (71,46)  NCHAOS, NAMEC(J), 
     &                    NAME(N+JPAIR), NAME(2*JPAIR)
 46                  FORMAT (' MERGED SPIRAL    NCH NAM ',I4,3I7)
                     CALL FLUSH(71)
                     GO TO 70
                  END IF
               END IF 
 48         CONTINUE
         END IF
*     
*     Update chaos variables for #IC and any disrupted or escaped binaries.
 50      NCHAOS = NCHAOS - 1
*     Copy chaos index in case of KS termination.
         IF (IREM.GT.0) THEN
            J = IREM
         END IF
*     
         DO 60 L = J,NCHAOS
            L1 = L + 1
            DO 55 K = 1,4
               EOSC(K,L) = EOSC(K,L1)
 55         CONTINUE
            EB0(L) = EB0(L1)
            ZJ0(L) = ZJ0(L1)
            ECRIT(L) = ECRIT(L1)
            AR(L) = AR(L1)
            BR(L) = BR(L1)
            EDEC(L) = EDEC(L1)
            TOSC(L) = TOSC(L1)
            RP(L) = RP(L1)
            ES(L) = ES(L1)
            CM(1,L) = CM(1,L1)
            CM(2,L) = CM(2,L1)
            IOSC(L) = IOSC(L1)
            NAMEC(L) = NAMEC(L1)
 60      CONTINUE
*     Ensure next chaos location contains zero core masses.
         CM(1,NCHAOS+1) = 0.0
         CM(2,NCHAOS+1) = 0.0
*     Consider the same location again after each removal (J <= NCHAOS).
         J = J - 1
 70      J = J + 1
         IF (J.LE.NCHAOS.AND.IREM.EQ.0) GO TO 40
      END IF
*     
*     Check optional diagnostics on transition to ECC = 0 or SLEEP.
      IF (KZ(8).GT.3.AND.KSTAR(I).NE.-2) THEN
*     Note case IPAIR < 0 from CHRECT with I denoting c.m.
         IF (IPAIR.LT.0) IPAIR = I - N
         CALL BINEV(IPAIR)
      END IF
*     

 100  RETURN
*     
      END
