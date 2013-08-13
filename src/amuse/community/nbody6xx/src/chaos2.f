      SUBROUTINE CHAOS2(I1,I2,ECC,HI,IS,BODYI,ZMU,RSTAR,SEMI1,ECC1,DH,
     &                                                      IDIS,KSTARI)
*
*
*       Chain chaotic treatment of tidal capture.
*       -----------------------------------------
*
*       Theory of Rosemary Mardling, Ap. J. XX, YYY, 1995.
*       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &            NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK,RSTAR(2)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
*       Note renaming of NAMEC to NAMEX in COMMON/CHREG/ to avoid conflict.
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEX(6),NSTEP1,KZ27,KZ30
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIJ(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &               BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &               RP(NTMAX),ES(NTMAX),CX(2,NTMAX),IOSC(NTMAX),
     &               NAMEC(NTMAX)
      REAL*8  DE2(2),DE3(2),WW(6),W(4),ALF(4),TL(2),AT(2),TDYN(2),WG(2),
     &        QG(2),WSCALE(2),QSCALE(2),EOSC0(2),QL(2),TD(2)
      REAL*8  RAN2
      INTEGER  IS(2),KG(2)
      CHARACTER*8  WHICH1
      SAVE  KICKS,NDEC
      DATA  WW  /2.119,3.113,8.175,3.742,4.953,9.413/
      DATA  ECCM2  /0.00000399/
*
*
*       Skip treatment for circular or hyperbolic stage (ISYNC = 1 or -1).
      IF (ISYNC.NE.0) GO TO 100
*
      IF (ECC.LT.0.0021) THEN
          KSTARI = 10
          ISYNC = 1
          GO TO 100
      END IF
*
*       See whether binary c.m. can be identified from the chaos table.
      IC = 0
      NAM1 = NZERO + NAMEX(I1)
      NAM2 = NZERO + NAMEX(I2)
      DO 1 K = 1,NCHAOS
          IF (NAMEC(K).EQ.NAM1.OR.NAMEC(K).EQ.NAM2) IC = K
    1 CONTINUE
*
*       Skip existing spiral (long t_circ for both original and new case).
      IF (IC.GT.0) THEN
          K = 0
          IF (NAMEX(I1) + NAMEX(I2).EQ.KSAVE(2)) K = 1
          IF (NAMEX(I1) + NAMEX(I2).EQ.KSAVE(4)) K = 3
          IF (K.GT.0) THEN
              KSTARI = KSAVE(K)
              IF (KSTARI.EQ.-2) THEN
                  ISYNC = 1
                  GO TO 100
              END IF
          END IF
      END IF
*
*       Increase counter for new chaos and initialize index & variables.
      IF (IC.EQ.0) THEN
          NCHAOS = NCHAOS + 1
          IC = NCHAOS
          NAMEC(IC) = NAM1
          NDEC = 0
          IOSC(IC) = 0
          KSTARI = -1
          TOSC(IC) = TIME + TIMEC
          DO 5 K = 1,4
              EOSC(K,IC) = 0.0D0
    5     CONTINUE
          IF (NCHAOS.GT.NTMAX) THEN
              if(rank.eq.0)
     &        WRITE (6,6)  NAME(I1), NCHAOS, ECC
    6         FORMAT (' FATAL ERROR!    CHAOS2    NM NCH E ',
     &                                            I6,I4,F8.4)
              STOP
          END IF
      END IF
*
*       Define oscillation period (dimensionless time) and damping constants.
      DO 10 K = 1,2
	  J = K + 2
          IK = I1
          IF (K.EQ.2) IK = I2
          TDYN(K) = SQRT(RSTAR(K)**3/M(IK))
*       Specify polytropic index for each star (n = 2, 3 or 3/2).
          IF (IS(K).EQ.3.OR.IS(K).EQ.5) THEN
*       Set fudged KS index = 0 to identify second particle in routine GIANT.
              CALL GIANT2(K,IK,WG,QG,WSCALE,QSCALE,ZN,QD)
              W(K) = WG(1)
              W(J) = WG(2)
              QL(K) = QD
              KG(K) = 1
          ELSE
              QL(K) = 1.0D+04
              KG(K) = 0
              IP = 3
              IF (IS(K).GE.3) IP = 2
              IF (IS(K).EQ.4.OR.IS(K).EQ.6) IP = 3
              IF (IS(K).EQ.0) IP = 1
              W(K) = WW(IP)
	      W(J) = WW(IP+3)
          END IF
          ALF(K) = 2.0*TDYN(K)/SQRT(W(K))
          ALF(J) = 3.0*TDYN(K)/SQRT(W(J))
          TL(K) = TWOPI*TDYN(K)/SQRT(W(K))
   10 CONTINUE
*
*       Save initial eccentricity, binding energy & J0.
      SEMI = -0.5*BODYI/HI
      CJ = ZMU*SQRT(BODYI)
      IF (IOSC(IC).EQ.0) THEN
          ECC0 = ECC
          EB0(IC) = ZMU*HI
          ZJ0(IC) = CJ*SQRT(QPERI*(1.0 + ECC0))
          EDEC(IC) = 0.0
          IOSC(IC) = 1
          TOSC(IC) = TIME + TIMEC
          RP(IC) = QPERI
          ES(IC) = ECC0
          KICKS = 0
*       Initialize chaos boundary parameters (ECRIT, AR & BR).
          CALL CHAOS0(QPERI,ECC,EB0(IC),ZJ0(IC),M(I1),M(I2),
     &                 RSTAR(1),RSTAR(2),W,ECRIT(IC),AR(IC),BR(IC),IDIS)
*
*       Check whether chaotic stage should be replaced by spiralling.
          IF (IDIS.EQ.-1) THEN
              IOSC(IC) = 2
              if(rank.eq.0)
     &        WRITE (6,12)  IS(1), IS(2), M(I1), M(I2), RSTAR(1),
     &                      RSTAR(2), QPERI, SEMI, ECC
   12         FORMAT (' CHAIN SPIRAL    K* M1 M2 R* QP A E ',
     &                                  2I4,1P,6E10.2,0P,F8.4)
*       Activate spiral indicator and save time, pericentre & eccentricity.
	      KSTARI = -2
              ISYNC = 1
              TOSC(IC) = TIME + TIMEC
              RP(IC) = QPERI
              ES(IC) = ECC0
              NAMC = 0
              DO 14 J = 1, NCHAOS
                  IF (NAMEC(J).EQ.NAM1.OR.NAMEC(J).EQ.NAM2) THEN
                      NAMC = NAMEC(J)
                  END IF
   14         CONTINUE
              IC = 0
*       Obtain diagnostic output after strong interaction.
              CALL RECOIL(3)
              GO TO 90
          END IF
*
          NCHA = NCHA + 1
          WHICH1 = ' CHAOS  '
          if(rank.eq.0)
     &    WRITE (6,15)  WHICH1, NAMEX(I1), NAMEX(I2), IS(1), IS(2),
     &                  M(I1)*ZMBAR, M(I2)*ZMBAR, RSTAR(1), RSTAR(2),
     &                  QPERI, SEMI, ECC
   15     FORMAT (' CHAIN',A8,'  NAM K* M1 M2 R* QP A E ',
     &                           2I6,2I4,2F5.1,1P,4E10.2,0P,F8.4)
*
*       Obtain diagnostic output after strong interaction.
          CALL RECOIL(3)
*       Exit on positive indicator (collision is too difficult here).
          IF (IDIS.GT.1) GO TO 90
      END IF
*
*       Obtain energy dissipation from separate modes.
      CALL TIDES2(QPERI,M(I1),M(I2),RSTAR(1),RSTAR(2),IS,ECC,KG,
     &                                            WSCALE,QSCALE,DE2,DE3)
*
*       Evaluate time-scale for Kochanek-type damping (linear & non-linear).
      EOSC0(1) = EOSC(1,IC)**2 + EOSC(2,IC)**2
      EOSC0(2) = EOSC(3,IC)**2 + EOSC(4,IC)**2
      EOSC0(1) = SQRT(EOSC0(1))
      EOSC0(2) = SQRT(EOSC0(2))
*       Adopt non-linear dissipation time scale of Kumar & Goodman 1995.
      QNL1 = QL(1)/MAX(SQRT(EOSC0(1)),0.00001D0)
      QNL2 = QL(2)/MAX(SQRT(EOSC0(2)),0.00001D0)
      SEMI = ABS(SEMI)
      TK = TWOPI*SEMI*SQRT(SEMI/BODYI)
      TD(1) = (1.0/QL(1) + 1.0/QNL1)*TK
      TD(2) = (1.0/QL(2) + 1.0/QNL2)*TK
*
*       Include check on wide interaction near the boundary.
      RM = MAX(RSTAR(1),RSTAR(2))
      IF (NDEC.GT.10000.AND.QPERI.GT.4.0*RM) THEN
          DO 18 K = 1,2
              DE2(K) = 1000.0*DE2(K)
              DE3(K) = 1000.0*DE3(K)
   18     CONTINUE
      END IF
*
*       Sum old and new oscillation energies for all modes.
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
     &                2.0*SQRT(EOSC(K,IC)*DE2(K))*COS(DELTA) + DE2(K)
          IF (IOSC(IC).EQ.-1) THEN
              IF (K.EQ.2) IOSC(IC) = 1
          ELSE
              DELTA = TWOPI*RAN2(IDUM1)
          END IF
          EOSC(J,IC) = EOSC(J,IC) +
     &                 2.0*SQRT(EOSC(J,IC)*DE3(K))*COS(DELTA) + DE3(K)
*       Ensure that oscillation energies are not negative.
          EOSC(K,IC) = MAX(EOSC(K,IC),0.0D0)
          EOSC(J,IC) = MAX(EOSC(J,IC),0.0D0)
          E2T = E2T + EOSC(K,IC)
          E3T = E3T + EOSC(J,IC)
	  ZJOSC = ZJOSC + ALF(K)*EOSC(K,IC) + ALF(J)*EOSC(J,IC)
   20 CONTINUE
*
*       Specify change in oscillation energy and sum decayed energy.
*     DET = (E2T - E20) + (E3T - E30)
      EDEC(IC) = EDEC(IC) + E20 + E30
*
*       Set new binding energy & semi-major axis.
      HNEW = (EB0(IC) - EDEC(IC) - (E2T + E3T))/ZMU
      DH = HNEW - HI
*       Note that the net change in EGRAV is taken care of at termination.
      SEMI1 = -0.5*BODYI/HNEW
*
*       Calculate the new eccentricity and specify pericentre.
      EFAC = (ZJ0(IC) - ZJOSC)/CJ
      ECC2 = 1.0 - EFAC**2/SEMI1
      ECC2 = MAX(ECC2,ECCM2)
      ECC1 = SQRT(ECC2)
      PERI1 = SEMI1*(1.0D0 - ECC1)
*
*       Switch off chaos indicator on transition to hyperbolic orbit.
      IF (HNEW.GT.0.0) THEN
          KSTARI = 0
          ISYNC = -1
*       Reduce index if current case is last (otherwise updated later).
          IF (IC.EQ.NCHAOS) THEN
              NCHAOS = NCHAOS - 1
          END IF
          IC = 0
          if(rank.eq.0)
     &    WRITE (6,25)  NAMEX(I1), NAMEX(I2), NDEC, KICKS, ECC, ECC1,
     &                  SEMI1, QPERI
   25     FORMAT (' TERMINATED CHAOS    NAM NDEC KICK E E1 A1 QP ',
     &                                  2I6,2I4,2F8.4,1P,2E10.2)
      ELSE IF (ECC.GT.1.0) THEN
          VINF = SQRT(2.0*HI)*VSTAR1
          P = DAYS*SEMI1*SQRT(SEMI1/BODYI)
          if(rank.eq.0)
     &    WRITE (6,26)  TIME+TOFF, NAMEX(I1), NAMEX(I2), ECC, ECC1,
     &                  VINF, SEMI1, P
   26     FORMAT (' CHAIN CAPTURE    T NAM E E1 VINF A1 P ',
     &                               F10.3,2I6,F9.5,F8.4,F5.1,1P,2E9.1)
      END IF
*
*       Check energy or eccentricity criterion for chaotic case.
      IF (IOSC(IC).EQ.1.OR.IOSC(IC).EQ.-1) THEN
          IF (EDEC(IC).GT.-(ECRIT(IC) - EB0(IC))) THEN
              IOSC(IC) = 2
              if(rank.eq.0)
     &        WRITE (6,30)  NAMEX(I1), NAMEX(I2), NDEC, KICKS, ECC1,
     &                      SEMI1, ECRIT(IC), EDEC(IC)
   30         FORMAT (' END CHAOS    NAM NDEC KICKS E A ECRIT EDEC ',
     &                               2I6,2I5,F8.4,1P,3E10.2)
*       Activate spiral indicator and save time, pericentre & eccentricity.
              KSTARI = -2
              ISYNC = 1
              TOSC(IC) = TIME + TIMEC
              RP(IC) = PERI1
              ES(IC) = ECC1
              IC = 0
          ELSE
              EPS = (2.0*ECRIT(IC) - EB0(IC) + EDEC(IC))/(ZMU*BODYI)
              ECCM = (1.0 + EPS*BR(IC))/(1.0 - EPS*AR(IC))
              IF (ECC1.LT.ECCM) THEN
                  IOSC(IC) = -1
                  KICKS = KICKS + 1
                  IF (KICKS.LT.3) THEN
                      if(rank.eq.0)
     &                WRITE (6,32)  NDEC, EPS, ECCM, ECC1
   32                 FORMAT (' NEW KICK    NDEC EPS ECCM E ',
     &                                      I5,1P,E10.2,0P,2F8.4)
                  END IF
              END IF
          END IF
      END IF
*
*       Reduce chaos index on disruption if current case is last.
*     IF (IDIS.GT.0.AND.IC.EQ.NCHAOS) THEN
*         NCHAOS = NCHAOS - 1
*     END IF
      NDEC = NDEC + 1
*
*       Skip energy change for small ECC1 and reduce NCHAOS if new SPIRAL.
      IF (ECC1.LT.0.0021) THEN
          KSTARI = 10
          IF (IC.EQ.0) NCHAOS = NCHAOS - 1
          IC = 0
      END IF
*
*       Update dissipation indicator in case of change.
   90 K = 0
      IF (NAMEX(I1) + NAMEX(I2).EQ.KSAVE(2)) K = 1
      IF (NAMEX(I1) + NAMEX(I2).EQ.KSAVE(4)) K = 3
      IF (K.GT.0) THEN
          KSAVE(K) = KSTARI
      ELSE IF (KSTARI.LT.0) THEN
*       Include possibility of exchange with KSTAR < 0.
          L = 1
          IF (NN.GT.3) L = 3
          KSAVE(L) = KSTARI
          KSAVE(L+1) = NAMEX(I1) + NAMEX(I2)
      END IF
*
  100 RETURN
*
      END
