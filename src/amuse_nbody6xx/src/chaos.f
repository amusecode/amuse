      SUBROUTINE CHAOS(IPAIR,I1,I2,QPERI,ECC,IS,ZMU,RKS,SEMI1,ECC1,IDIS)
*     
*     
*     Chaotic tidal interactions.
*     ---------------------------
*     
*     Theory of Rosemary Mardling, Ap. J. XX, YYY, 1995.
*     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*     
      INCLUDE 'common6.h'
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &     BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &     RP(NTMAX),ES(NTMAX),CM(2,NTMAX),IOSC(NTMAX),
     &     NAMEC(NTMAX)
      REAL*8  DE2(2),DE3(2),WW(6),W(4),ALF(4),TL(2),AT(2),TDYN(2),
     &     WG(2),QG(2),WSCALE(2),QSCALE(2),A0(3),A2(3),EOSC0(2),
     &     QL(2),TD(2)
      REAL*8  RAN2
      INTEGER  IS(2),KG(2)
      CHARACTER*8  WHICH1
      SAVE  KICKS,NDEC,EOSC0,ZJOSC
      DATA  WW  /2.119,3.113,8.175,3.742,4.953,9.413/
      DATA  ECCM2  /0.00000399/
*     
*     
*     Define c.m. index and search current names for chaos index.
      IDIS = 0
      I = N + IPAIR
      IC = 0
      DO 1 K = 1,NCHAOS
         IF (NAMEC(K).EQ.NAME(I)) IC = K
    1 CONTINUE
*     
*     See whether new case could be an old SPIRAL exchanged in chain.
      IF (IC.GT.0) THEN
*     Remove tables and re-initialize on large epoch or standard binary.
         IF (TIME - TOSC(IC).GT.1.0.OR.KSTAR(I).EQ.0) THEN
            II = -I
            CALL SPIRAL(II)
            IC = 0
         END IF
      END IF
*     
*     Increase counter for new chaos and initialize index & variables.
      IF (IC.EQ.0) THEN
         NCHAOS = NCHAOS + 1
         IC = NCHAOS
         NAMEC(IC) = NAME(I)
         KICKS = 0
         NDEC = 0
         IOSC(IC) = 0
         KSTAR(I) = -1
         TOSC(IC) = TIME
*     Ensure next location contains zero core mass (avoids confusion).
         CM(1,NCHAOS+1) = 0.0
         CM(2,NCHAOS+1) = 0.0
         DO 5 K = 1,4
            EOSC(K,IC) = 0.0D0
 5       CONTINUE
         IF (KZ(8).GT.3.AND.H(IPAIR).LT.0.0) THEN
            CALL BINEV(IPAIR)
         END IF
         IF (NCHAOS.GT.NTMAX) THEN
            if(rank.eq.0)WRITE (6,6)  NAME(I1), NCHAOS, ECC, QPERI
 6          FORMAT (' FATAL ERROR!    CHAOS    NM NCH E QP ',
     &           I6,I4,F8.4,1P,E9.1)
            if(rank.eq.0)WRITE (6,7)  (NAMEC(K),K=1,NCHAOS)
 7          FORMAT (' NAMEC  ',12I6,(/,12I6))
            STOP
         END IF
      END IF
*     
*     Check termination for rare case of wide chaos (minimum 1000 calls).
      QP = QPERI/MAX(RADIUS(I1),RADIUS(I2))
      IF (QP.GT.7.0.AND.NDEC.GT.1000) THEN
*     Activate spiral indicator and save time, pericentre & eccentricity.
         KSTAR(I) = -2
         TOSC(IC) = TIME
         RP(IC) = QPERI
         ES(IC) = ECC
         IOSC(IC) = 2
         NSP = NSP + 1
         SEMI = -0.5*BODY(I)/H(IPAIR)
         if(rank.eq.0)WRITE (6,8)  NAME(I1), NAME(I2), NDEC, 
     &        LIST(1,I1), QP, ECC, SEMI
 8       FORMAT (' WIDE CHAOS    NAM NDEC NP QP/S E A ',
     &        3I6,I4,F5.1,F8.4,1P,E10.2)
         NDEC = 0
         KICKS = 0
         GO TO 34
      END IF
*     
*     Define oscillation period (dimensionless time) and damping constants.
      ZN = 0.0
      QD = 0.0
      DO 10 K = 1,2
         J = K + 2
         IK = I1 + K - 1
         TDYN(K) = RADIUS(IK)*SQRT(RADIUS(IK)/BODY(IK))
*     Specify polytropic index for each star (n = 2, 3 or 3/2).
         IF (KSTAR(IK).EQ.3.OR.KSTAR(IK).EQ.5.OR.
     &        KSTAR(IK).EQ.6.OR.KSTAR(IK).EQ.9) THEN
            CALL GIANT(IPAIR,IK,WG,QG,WSCALE,QSCALE,ZN,QD)
            W(K) = WG(1)
            W(J) = WG(2)
            QL(K) = QD
            KG(K) = 1
         ELSE
            QL(K) = 1.0D+04
            KG(K) = 0
            IP = 3
            IF (KSTAR(IK).GE.3) IP = 2
            IF (KSTAR(IK).EQ.4.OR.KSTAR(IK).EQ.7) IP = 3
            IF (KSTAR(IK).EQ.8) IP = 3
            IF (KSTAR(IK).EQ.0) IP = 1
            W(K) = WW(IP)
            W(J) = WW(IP+3)
         END IF
         ALF(K) = 2.0*TDYN(K)/SQRT(W(K))
         ALF(J) = 3.0*TDYN(K)/SQRT(W(J))
         TL(K) = TWOPI*TDYN(K)/SQRT(W(K))
 10   CONTINUE
*     
*     Save initial eccentricity, binding energy & J0 the first time.
      CJ = ZMU*SQRT(BODY(I))
      IF (IOSC(IC).EQ.0.OR.LIST(1,I1).GT.0) THEN
         SEMI = -0.5*BODY(I)/H(IPAIR)
         ECC0 = ECC
         EB0(IC) = ZMU*H(IPAIR)
         ZJ0(IC) = CJ*SQRT(QPERI*(1.0 + ECC0))
         EDEC(IC) = 0.0
         IOSC(IC) = 1
*     Reset diagnostic indicator for CHAOS0 (output only first time).
         IF (NDEC.GT.0) IDIS = KSTAR(I)
*     
*     Initialize chaos boundary parameters (ECRIT, AR & BR).
         CALL CHAOS0(QPERI,ECC,EB0(IC),ZJ0(IC),BODY(I1),BODY(I2),
     &        RADIUS(I1),RADIUS(I2),W,ECRIT(IC),AR(IC),BR(IC),IDIS)
*     
*     Begin spiralling stage if chaos boundary has been crossed.
         IF (IDIS.EQ.-1) THEN
*     Include safety check on skipping rare case of hyperbolic orbit.
            IF (ECC.GT.1.0) THEN
               NCHAOS = NCHAOS - 1
               KSTAR(I) = 0
               GO TO 80
            END IF
*     Restrict tidal circularization to short time-scales (< 100 Myr).
            ICIRC = -1
            CALL TCIRC(QPERI,ECC,I1,I2,ICIRC,TC)
            IF (TC.GT.100.0) THEN
               NCHAOS = NCHAOS - 1
               KSTAR(I) = 0
               GO TO 80
            END IF
*     Skip initialization for difficult hierarchy and small EMAX.
            IF (LIST(1,I1).EQ.1) THEN
               ICIRC = -1
               JCOMP = LIST(2,I1)
               CALL INDUCE(IPAIR,EMAX,EMIN,ICIRC,TC,ANGLE,TG,EDAV)
               IF (EMAX.LT.0.9) THEN
                  NCHAOS = NCHAOS - 1
                  KSTAR(I) = 0
                  GO TO 80
               END IF
            END IF
            IOSC(IC) = 2
            if(rank.eq.0)
     &           WRITE (6,12)  TTOT, NAME(I1), NAME(I2), KSTAR(I1),
     &           KSTAR(I2), LIST(1,I1), BODY(I1)*ZMBAR,
     &           BODY(I2)*ZMBAR, RADIUS(I1)*SU,
     &           RADIUS(I2)*SU, QPERI, SEMI, ECC, ZN, QD
 12         FORMAT (' NEW SPIRAL    T NM K* NP M1 M2 R* QP A E n Q ',
     &           F9.2,2I6,3I3,2F5.1,2F6.1,
     &           1P,2E10.2,0P,F7.3,F5.1,F7.1)
*     Activate spiral indicator and save time, pericentre & eccentricity.
            KSTAR(I) = -2
            TOSC(IC) = TIME
            RP(IC) = QPERI
            ES(IC) = ECC0
            NSP = NSP + 1
*     Rectify the KS solution on transition to standard circularization.
            CALL KSRECT(IPAIR)
*     Initialize perturbed KS (small STEP after integration in KSPERI).
            IF (LIST(1,I1).GT.0) THEN
               IMOD = 1
               CALL KSPOLY(IPAIR,IMOD)
               KSLOW(IPAIR) = 1
            END IF
            GO TO 34
         END IF
*     
*     Reduce chaos index on disruption if current case is last.
         IF (IDIS.GT.0.AND.IC.EQ.NCHAOS) THEN
            NCHAOS = NCHAOS - 1
            GO TO 80
         END IF
*     
*     Print NEW CHAOS/CAPTURE the first time.
         IF (NDEC.EQ.0) THEN
            NCHA = NCHA + 1
            WHICH1 = ' CHAOS  '
            IF (SEMI.LT.0.0) WHICH1 = ' CAPTURE'
            GA = GAMMA(IPAIR)*(RMIN/RKS)**3
            if(rank.eq.0)
     &           WRITE (6,15)  WHICH1, TTOT, NAME(I1), NAME(I2), 
     &           KSTAR(I1), KSTAR(I2), LIST(1,I1), BODY(I1)*ZMBAR,
     &           BODY(I2)*ZMBAR, RADIUS(I1)*SU,
     &           RADIUS(I2)*SU, QPERI, SEMI, ECC, GA, ZN, QD
 15         FORMAT (' NEW',A8,'  T NM K* NP M1 M2 R* QP A E G n Q ',
     &           F9.2,2I6,3I3,2F5.1,2F6.1,1P,2E10.2,
     &           0P,F9.5,F7.3,F5.1,F7.1)
         END IF
      END IF
*     
*     Obtain energy dissipation from separate modes.
      CALL TIDES2(QPERI,BODY(I1),BODY(I2),RADIUS(I1),RADIUS(I2),IS,ECC,
     &     KG,WSCALE,QSCALE,DE2,DE3)
*     
*     Employ an emergency procedure for zero energy change (TIDES2 bug).
      IF (DE2(1).EQ.0.0D0.OR.DE2(2).EQ.0.0D0) THEN
         EB = BODYM*H(IPAIR)
         DE2(1) = -0.001*EB*RAN2(IDUM)
         DE2(2) = -0.001*EB*RAN2(IDUM)
      END IF
*     
*     Evaluate time-scale for Kochanek-type damping (linear & non-linear).
      EOSC0(1) = EOSC(1,IC)**2 + EOSC(2,IC)**2
      EOSC0(2) = EOSC(3,IC)**2 + EOSC(4,IC)**2
      EOSC0(1) = SQRT(EOSC0(1))
      EOSC0(2) = SQRT(EOSC0(2))
*     Adopt non-linear dissipation time scale of Kumar & Goodman 1995.
      QNL1 = QL(1)/MAX(SQRT(EOSC0(1)),0.00001D0)
      QNL2 = QL(2)/MAX(SQRT(EOSC0(2)),0.00001D0)
      SEMI = -0.5*BODY(I)/H(IPAIR)
      SEMI = ABS(SEMI)
      TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
      TD(1) = (1.0/QL(1) + 1.0/QNL1)*TK
      TD(2) = (1.0/QL(2) + 1.0/QNL2)*TK
*     
      RM = MAX(RADIUS(I1),RADIUS(I2))
      IF (NDEC.GT.10000.AND.QPERI.GT.4.0*RM) THEN
         DO 18 K = 1,2
            DE2(K) = 1000.0*DE2(K)
            DE3(K) = 1000.0*DE3(K)
 18      CONTINUE
      END IF
*     Sum old and new oscillation energies for all modes.
      E20 = 0.0
      E30 = 0.0
      E2T = 0.0
      E3T = 0.0
      ZJOSC = 0.0
      DO 20 K = 1,2
         J = K + 2
         AT(K) = -TD(K)/TL(K)
         E20 = E20 + EOSC(K,IC)
         E30 = E30 + EOSC(J,IC)
         EOSC(K,IC) = EOSC(K,IC)*EXP(AT(K))
         EOSC(J,IC) = EOSC(J,IC)*EXP(AT(K))
         EDEC(IC) = EDEC(IC) - (EOSC(K,IC) + EOSC(J,IC))
         IF (IOSC(IC).EQ.-1) THEN
            DELTA = 0.5*TWOPI
         ELSE
            DELTA = TWOPI*RAN2(IDUM1)
         END IF
         EOSC(K,IC) = EOSC(K,IC) +
     &        2.0*SQRT(EOSC(K,IC)*DE2(K))*COS(DELTA) + DE2(K)
         IF (IOSC(IC).EQ.-1) THEN
            IF (K.EQ.2) IOSC(IC) = 1
         ELSE
            DELTA = TWOPI*RAN2(IDUM1)
         END IF
         EOSC(J,IC) = EOSC(J,IC) +
     &        2.0*SQRT(EOSC(J,IC)*DE3(K))*COS(DELTA) + DE3(K)
*     Ensure that oscillation energies are not negative.
         EOSC(K,IC) = MAX(EOSC(K,IC),0.0D0)
         EOSC(J,IC) = MAX(EOSC(J,IC),0.0D0)
         E2T = E2T + EOSC(K,IC)
         E3T = E3T + EOSC(J,IC)
         ZJOSC = ZJOSC + ALF(K)*EOSC(K,IC) + ALF(J)*EOSC(J,IC)
 20   CONTINUE
*     
*     Specify change in oscillation energy and sum decayed energy.
*     DET = (E2T - E20) + (E3T - E30)
      EDEC(IC) = EDEC(IC) + E20 + E30
      TOSC(IC) = TIME
*     
*     Set new binding energy & semi-major axis.
      HI = H(IPAIR)
      HNEW = (EB0(IC) - EDEC(IC) - (E2T + E3T))/ZMU
      SEMI1 = -0.5*BODY(I)/HNEW
*     
*     Calculate the new eccentricity.
      EFAC = (ZJ0(IC) - ZJOSC)/CJ
      ECC2 = 1.0 - EFAC**2/SEMI1
      ECC2 = MAX(ECC2,ECCM2)
      ECC1 = SQRT(ECC2)
      PERI1 = SEMI1*(1.0D0 - ECC1)
*     
*     Switch off chaos indicator on transition to hyperbolic orbit.
      IF (HNEW.GT.0.0) THEN
         KSTAR(I) = 0
*     Reduce index if current case is last (otherwise updated in KSTERM).
         IF (IC.EQ.NCHAOS) THEN
            NCHAOS = NCHAOS - 1
         END IF
         if(rank.eq.0)
     &        WRITE (6,25)  IPAIR, NDEC, KICKS, ECC, ECC1, SEMI1, QPERI,
     &        HNEW - HI
 25      FORMAT (' TERMINATED CHAOS    IPAIR NDEC KICK E E1 A1 QP DH ',
     &        3I4,2F9.5,1P,3E10.2)
         GO TO 80
      END IF
*     
*     Update total energy loss due to change in binding energy.
      NDEC = NDEC + 1
      H(IPAIR) = HNEW
      DEB = ZMU*(HI - H(IPAIR))
      ECOLL = ECOLL + DEB
      EGRAV = EGRAV + DEB
      E(10) = E(10) + DEB
*     
*     Check energy or eccentricity criterion for chaotic case.
      IF (KZ(27).EQ.2.AND.(IOSC(IC).EQ.1.OR.IOSC(IC).EQ.-1)) THEN
         IF (EDEC(IC).GT.-(ECRIT(IC) - EB0(IC))) THEN
            IOSC(IC) = 2
            if(rank.eq.0)
     &           WRITE (6,30)  TTOT, NAME(I1), NAME(I2), NDEC, KICKS,
     &           ECC1, SEMI1, ECRIT(IC), EDEC(IC)
 30         FORMAT (' END CHAOS    T NM NDEC KICKS E A ECRIT EDEC ',
     &           F9.2,4I6,F8.4,1P,3E10.2)
*     Activate spiral indicator and save time, pericentre & eccentricity.
            KSTAR(I) = -2
            TOSC(IC) = TIME
            RP(IC) = PERI1
            ES(IC) = ECC1
            NSP = NSP + 1
            NDEC = 0
            KICKS = 0
         ELSE
            EPS = (2.0*ECRIT(IC) - EB0(IC) + EDEC(IC))/(ZMU*BODY(I))
            ECCM = (1.0 + EPS*BR(IC))/(1.0 - EPS*AR(IC))
            IF (ECC1.LT.ECCM) THEN
               IOSC(IC) = -1
               KICKS = KICKS + 1
               IF (KICKS.LT.3) THEN
                  if(rank.eq.0)WRITE (6,32)  NDEC, EPS, ECCM, ECC1
 32               FORMAT (' NEW CHAOS KICK     NDEC EPS ECCM E ',
     &                 I5,1P,E10.2,0P,2F8.4)
               END IF
            END IF
         END IF
      END IF
*     
*     Check for hierarchical configuration on first call.
 34   IF (NDEC.LE.1.AND.SEMI.GT.0.0.AND.SEMI.LT.2.0*RMIN) THEN
         NP1 = LIST(1,I1) + 1
         DO 39 L = 2,NP1
            J = LIST(L,I1)
            RIJ2 = 0.0
            VIJ2 = 0.0
            RDOT = 0.0
            A12 = 0.0
            A22 = 0.0
            A1A2 = 0.0
            DO 35 K = 1,3
               RIJ2 = RIJ2 + (X(K,I) - X(K,J))**2
               VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,J))**2
               RDOT = RDOT + (X(K,I) - X(K,J))*(XDOT(K,I) -XDOT(K,J))
               K1 = K + 1
               IF (K1.GT.3) K1 = 1
               K2 = K1 + 1
               IF (K2.GT.3) K2 = 1
               A0(K) = (X(K1,I1)-X(K1,I2))*(XDOT(K2,I1)-XDOT(K2,I2))
     &              - (X(K2,I1)-X(K2,I2))*(XDOT(K1,I1)-XDOT(K1,I2))
               A2(K) = (X(K1,J) - X(K1,I))*(XDOT(K2,J) - XDOT(K2,I))
     &              - (X(K2,J) - X(K2,I))*(XDOT(K1,J) - XDOT(K1,I))
               A12 = A12 + A0(K)**2
               A22 = A22 + A2(K)**2
               A1A2 = A1A2 + A0(K)*A2(K)
 35         CONTINUE
            RIP = SQRT(RIJ2)
            A1 = 2.0/RIP - VIJ2/(BODY(I) + BODY(J))
            A1 = 1.0/A1
*     Include impact parameter test for recoil check.
            A4 = RDOT**2/(A1*(BODY(I) + BODY(J)))
            ECCP = SQRT((1.0D0 - RIP/A1)**2 + A4)
            PMIN = A1*(1.0D0 - ECCP)
            IF (PMIN.LT.3.0*SEMI) THEN
               TM = MIN(TEV(I1),TEV(I2)) - TIME
               if(rank.eq.0)
     &              WRITE (6,36)  IPAIR, NAME(J), ECCP, PMIN/SEMI,
     &              RDOT/RIP, SEMI, A1, RIP, TM
 36            FORMAT (' RECOIL:    KS NMJ E1 PM/A RD A0 A1 RP TM ',
     &              I4,I6,F7.3,F5.1,F6.1,1P,4E10.2)
            END IF
*     Accept semi-major axis ratio below 25.
            IF (1.0/A1.GT.0.04/SEMI.AND.IDIS.LE.0) THEN
               RA = SEMI*(1.0 + ECC)
               SR = PMIN/RA
               GA = 2.0*BODY(J)*(RA/PMIN)**3/BODY(I)
*     Determine inclination (8 bins of 22.5 degrees).
               FAC = A1A2/SQRT(A12*A22)
               FAC = ACOS(FAC)
               IN = 1 + FAC*360.0/(TWOPI*22.5)
               if(rank.eq.0)
     &              WRITE (6,38)  IPAIR, NAME(J), H(IPAIR), SEMI, A1,
     &              PMIN, GA, ECCP, SR, IN
 38            FORMAT (' HIERARCHY:    KS NMJ H A0 A1 RP GA E1 SR ',
     &              'IN',I5,I6,F7.0,1P,3E9.1,0P,2F6.2,F6.1,I3)
            END IF
 39      CONTINUE
      END IF
*     
*     Check for terminated or escaped chaotic binaries first time.
      IF (NDEC.EQ.1.AND.NCHAOS.GT.1) THEN
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
*     Skip during mergers or multiple regularizations (KSTAR not visible).
         IF (NMERGE.GT.0.OR.NSUB.GT.0) GO TO 70
*     
*     Update chaos variables for #IC and any disrupted or escaped binaries.
 50      NCHAOS = NCHAOS - 1
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
*     Consider the same location again after each removal (J <= NCHAOS).
         J = J - 1
 70      J = J + 1
         IF (J.LE.NCHAOS) GO TO 40
      END IF
*     
*     Check optional binary diagnostics on transition from chaotic state.
 80   IF (KZ(8).GT.3.AND.KSTAR(I).NE.-1) THEN
*     Skip output for unperturbed case (CALL KSTIDE from UNPERT).
         IF (LIST(1,I1).GT.0.AND.H(IPAIR).LT.0.0) THEN
            CALL BINEV(IPAIR)
         END IF
      END IF
*     
      RETURN
*     
      END
