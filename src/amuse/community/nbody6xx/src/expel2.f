      SUBROUTINE EXPEL2(J1,J2,ICASE)
*
*
*       Common envelope stage for chain system.
*       ---------------------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &            NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK,XCM(3),XREL(3),VCM(3),VREL(3)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIK(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/EBSAVE/  EBS
      COMMON/KSAVE/  K10,K20
      REAL*8 LUMS(10),TSCLS(20),GB(10),M01,M1,TM,TN,LUM1,LUM2,AJ1,R1,
     &       M02,M2,AJ2,R2,SEP,MI(2),MC1,MC2,RCC
      REAL*8 JSPIN1,JSPIN2,MENV,RENV,K2
      LOGICAL COALS
*
*
*       Define global indices such that body #I1 is giant-like.
      IF(KSTAR(J1).GE.2.AND.KSTAR(J1).LE.9.AND.KSTAR(J1).NE.7)THEN
          I1 = J1
          I2 = J2
      ELSE
          I1 = J2
          I2 = J1
      ENDIF
*
*       Save original total mass and reduced mass.
      ZMB0 = BODY(I1) + BODY(I2)
      ZMU0 = BODY(I1)*BODY(I2)/ZMB0
*
*       Accumulate c.m. variables for dominant bodies.
      RIJ2 = 0.d0
      VIJ2 = 0.d0
      VCM2 = 0.d0
      RDOT = 0.d0
      DO 5 K = 1,3
          XREL(K) = X(K,I1) - X(K,I2)
          VREL(K) = XDOT(K,I1) - XDOT(K,I2)
          RIJ2 = RIJ2 + XREL(K)**2
          VIJ2 = VIJ2 + VREL(K)**2
          RDOT = RDOT + XREL(K)*VREL(K)
          XCM(K) = (BODY(I1)*X(K,I1) + BODY(I2)*X(K,I2))/ZMB0
          VCM(K) = (BODY(I1)*XDOT(K,I1) + BODY(I2)*XDOT(K,I2))/ZMB0
          VCM2 = VCM2 + VCM(K)**2
    5 CONTINUE
*
*       Form binding energy per unit mass and original semi-major axis.
      RIJ0 = SQRT(RIJ2)
      DMINC = RIJ0
      HI = 0.5d0*VIJ2 - ZMB0/RIJ0
      SEMI0 = 2.d0/RIJ0 - VIJ2/ZMB0
      SEMI0 = 1.d0/SEMI0
      SEMI00 = SEMI0
*
      ECC2 = (1.D0 - RIJ0/SEMI0)**2 + RDOT**2/(SEMI0*ZMB0)
      ECC = SQRT(ECC2)
*
*       Obtain SEMI0 & HI from EBS because of some velocity bug (19/8/96).
      SEMI0 = -0.5d0*M(K10)*M(K20)/EBS
      HI = -0.5d0*(M(K10) + M(K20))/SEMI0
*
*       Define non-dominant members as perturbers (JLIST(4) = 0 if NCH = 3).
      I3 = JLIST(3)
      I4 = JLIST(4)
      JPERT(1) = I3
      NP = 1
      IF (NCH.EQ.4) THEN
          JPERT(2) = I4
          NP = 2
      END IF
*       Evaluate potential energy of #I3/I4 before mass loss & displacement.
      CALL NBPOT(2,NP,POT1)
*
*       Update the stars to previous latest time (only for original KS pair).
*     ID = 0
*       Check original identity (NCH=3 may have been reduced from B-B).
*     DO 2 L = 1,2
*         IF (NAME(I1) + NAME(I2).EQ.KSAVE(2*L)) THEN
*             ID = 1
*         END IF
*   2 CONTINUE
      TEV1 = MAX(TEV0(I1),TEV0(I2))
*
*       Specify basic parameters for both stars.
      M01 = BODY0(I1)*ZMBAR
      M1 = BODY(I1)*ZMBAR
      MC1 = 0.D0
      AJ1 = TEV1*TSTAR - EPOCH(I1)
      JSPIN1 = SPIN(I1)*SPNFAC
      KW1 = KSTAR(I1)
      M02 = BODY0(I2)*ZMBAR
      M2 = BODY(I2)*ZMBAR
      MC2 = 0.D0
      AJ2 = TEV1*TSTAR - EPOCH(I2)
      JSPIN2 = SPIN(I2)*SPNFAC
      KW2 = KSTAR(I2)
      SEP = SEMI0*SU
      ECC0 = ECC
      ITER = 0
*
*       Perform common envelope evolution (note: SEP in SU).
    1 CALL COMENV(M01,M1,MC1,AJ1,JSPIN1,KW1,
     &            M02,M2,MC2,AJ2,JSPIN2,KW2,ECC,SEP,COALS)
*
*       Obtain consistent radii for the stars (skip #I2 on coalescence).
      CALL star(KW1,M01,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
      CALL hrdiag(M01,AJ1,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &            R1,LUM1,KW1,MC1,RCC,MENV,RENV,K2)
      if(rank.eq.0.and.kw1.ne.kstar(i1))then
         write(38,*)' EXPEL2 TYPE CHANGE1 ',kstar(i1),kw1
         write(38,*)' EXPEL2 TYPE CHANGE1 ',i1,name(i1),time
      endif
      IF(COALS)THEN
         KW2 = KSTAR(I2)
         R2 = 0.d0
         LUM2 = 0.d0
         JSPIN2 = 0.d0
      ELSE
          CALL star(KW2,M02,M2,TM,TN,TSCLS,LUMS,GB,ZPARS)
          CALL hrdiag(M02,AJ2,M2,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                R2,LUM2,KW2,MC2,RCC,MENV,RENV,K2)
          if(rank.eq.0.and.kw2.ne.kstar(i2))then
             write(38,*)' EXPEL2 TYPE CHANGE2 ',kstar(i2),kw2
          endif
      ENDIF
*
      DMSUN = M1 + M2 - ZMB0*SMU
      IF (ABS(DMSUN).LT.1.0D-12.AND..NOT.COALS) THEN
*         if(rank.eq.0)
*    &    WRITE (78,7) TIME, ECC, DMSUN, SEMI0*SU-SEP
*   7     FORMAT (' FAIL - CHAIN   T E DMS DSEP ',F10.4,F8.4,1P,2E10.2)
*         CALL FLUSH(78)
          IPHASE = 0
          GO TO 100
      END IF
*
*       Check common envelope condition again after circularization (10/08).
      IF (.NOT.COALS) THEN
          IF (ECC0.GT.0.001D0.AND.ECC.LE.0.001D0) THEN
              Q1 = M1/M2
              Q2 = 1.D0/Q1
              RL1 = RL(Q1)
              RL2 = RL(Q2)
              IF (R1.GT.RL1*SEP.OR.R2.GT.RL2*SEP) THEN
                  ITER = ITER + 1
                  IF (ITER.EQ.1) THEN
                      if(rank.eq.0)
     &                WRITE (6,8)  RL1*SEP, RL2*SEP
    8                 FORMAT (' ENFORCED CHAIN CE    RL*SEP ',1P,2E10.2)
                      ECC0 = ECC
                      GO TO 1
                  END IF
              END IF
          END IF
      END IF
*
*       Copy new values of radius and luminosity.
      RADIUS(I1) = R1/SU
      RADIUS(I2) = R2/SU
      ZLMSTY(I1) = LUM1
      ZLMSTY(I2) = LUM2
      SPIN(I1) = JSPIN1/SPNFAC
      SPIN(I2) = JSPIN2/SPNFAC
*
*       Update initial & chain masses (M2 = 0 for coalescence).
      BODY0(I1) = M01/ZMBAR
      BODY0(I2) = M02/ZMBAR
      M(K10) = M1/ZMBAR
      M(K20) = M2/ZMBAR
      ECOLL = ECOLL + ZMU0*HI
      CHCOLL = CHCOLL + ZMU0*HI
      EGRAV = EGRAV + ZMU0*HI
*
*       Distinguish between coalescence and surviving binary.
      EPOCH(I1) = TEV1*TSTAR - AJ1
      IF(COALS)THEN
          MI(1) = M1
          MI(2) = M2
          BODY0(I2) = 0.D0
*       Check for TZ object formed by CE evolution.
          IF(KSTAR(I2).GE.13.AND.KW1.GE.13)THEN
              NTZ = NTZ + 1
              if(rank.eq.0)
     &        WRITE (6,10)  M1, M2
   10         FORMAT (' NEW TZ    M1 M2 ',2F7.2)
          ENDIF
          IPAIR = -1
          CALL COAL(IPAIR,KW1,KW2,MI)
*       Reverse case indicator to denote exit from collision routine.
          ICASE = -ICASE
      ELSE
*
*       Update evolution times and TMDOT.
          EPOCH(I2) = TEV1*TSTAR - AJ2
          TEV(I1) = TEV1
          TEV0(I1) = TEV(I1)
          TEV(I2) = TEV1
          TEV0(I2) = TEV(I2)
          TMDOT = MIN(TMDOT,TEV1)
*
*       Copy new semi-major axis & masses and specify mass loss.
          SEMI = SEP/SU
          BODY(I1) = M1/ZMBAR
          BODY(I2) = M2/ZMBAR
          ZMB = BODY(I1) + BODY(I2)
          ZMU = BODY(I1)*BODY(I2)/ZMB
          DM = ZMB0 - ZMB
          ZMASS = ZMASS - DM
*
*       Copy chain perturbers and obtain potential energy w.r.t. I1 & I2.
*         NP = 0
*         DO 15 L = 3,NCH
*             NP = NP + 1
*             JPERT(NP) = JLIST(L)
*  15     CONTINUE
*         CALL NBPOT(2,NP,POT1)
*
*       Distinguish between bound and hyperbolic orbit.
          IF (ECC.LT.1.0) THEN
*
*       Specify eccentricity factor for transformation to apocentre.
              EFAC = SQRT((1.0 - ECC)/(1.0 + ECC))
*
*       Set new coordinates and velocities for relative & absolute motion.
              DO 20 K = 1,3
                  XREL(K) = XREL(K)*SEMI/RIJ0
                  XREL(K) = (1.0 + ECC)*XREL(K)
                  X(K,I1) =  XCM(K) + BODY(I2)*XREL(K)/ZMB
                  X(K,I2) =  XCM(K) - BODY(I1)*XREL(K)/ZMB
                  VREL(K) = SQRT(ZMB/(SEMI*VIJ2))*VREL(K)
                  VREL(K) = EFAC*VREL(K)
                  XDOT(K,I1) = VCM(K) + BODY(I2)*VREL(K)/ZMB
                  XDOT(K,I2) = VCM(K) - BODY(I1)*VREL(K)/ZMB
   20         CONTINUE
*
          ELSE
*       Form relative velocity by ratio of new and old pericentre value.
              V2 = ZMB*(2.0/RIJ0 - 1.0/SEMI)
              V20 = VREL(1)**2 + VREL(2)**2 + VREL(3)**2
              DO 30 K = 1,3
                  VREL(K) = SQRT(V2/V20)*VREL(K)
                  XDOT(K,I1) = VCM(K) + BODY(I2)*VREL(K)/ZMB
                  XDOT(K,I2) = VCM(K) - BODY(I1)*VREL(K)/ZMB
                  X(K,I1) =  XCM(K) + BODY(I2)*XREL(K)/ZMB
                  X(K,I2) =  XCM(K) - BODY(I1)*XREL(K)/ZMB
   30         CONTINUE
          END IF
*
*       Evaluate new interaction energy and perform differential correction.
*         CALL NBPOT(2,NP,POT2)
*         ECOLL = ECOLL + (POT2 - POT1)
*
*       Specify new energies (note BODY(I2) & ZMU = 0 for coalescence).
          HF = -0.5d0*ZMB/SEMI
*
*       Update energy corrections (change in H and mass loss kinetic energy). 
          ECOLL = ECOLL - ZMU*HF
          CHCOLL = CHCOLL - ZMU*HF
          ECDOT = ECDOT + 0.5*DM*VCM2
          EGRAV = EGRAV - ZMU*HF
*
*       Include potential energy terms due to all non-chain members.
          POTJ = 0.d0
      RX = 10.
*       Note no initialization of neighbours done in the chain case.
          DO 50 J = IFIRST,NTOT
              IF (J.NE.I1.AND.J.NE.I2.AND.J.NE.I3.AND.J.NE.I4) THEN
                  RIJ2 = (X(1,J) - XCM(1))**2 + 
     &                   (X(2,J) - XCM(2))**2 +
     &                   (X(3,J) - XCM(3))**2
                  POTJ = POTJ - BODY(J)/SQRT(RIJ2)
      RX = MIN(RX,RIJ2)
              END IF
   50     CONTINUE
*
*       Obtain non-dominant chain contributions after mass loss and new XREL.
          CALL NBPOT(2,NP,POT2)
*
*       See whether tidal terms should be included.
          IF (KZ(14).EQ.1) THEN
              ECDOT = ECDOT - 0.5d0*DM*(TIDAL(1)*XCM(1)**2 +
     &                                  TIDAL(3)*XCM(3)**2)
          END IF
*
*       Add both potential energy contributions to yield final correction.
          ECDOT = ECDOT + DM*POTJ + (POT2 - POT1)
*
          if(rank.eq.0)then
          WRITE (6,55)  DM*POTJ, SQRT(RX), SEMI00*SU, SEMI0*SU, EBS
   55     FORMAT (' EXPEL2!    DM*POTJ RX A00 A0 EBS ',1P,6E10.2)
          WRITE (6,60)  (NAME(JLIST(K)),K=1,4), KSTAR(I1), KSTAR(I2),
     &                  KW1, KW2, M1, M2, DM*ZMBAR, ECC, R1, R2,
     &                  SEMI0*SU, SEMI*SU, EBS, POT2-POT1
   60     FORMAT (' CHAIN CE    NAM K0* K* M1 M2 DM E R1 R2 A0 A EB DP',
     &                          4I6,4I3,3F5.1,F8.4,2F6.1,2F8.1,1P,2E9.1)
          end if
*       Flag hyperbolic chain CE.
          IF (rank.eq.0.and.SEMI0.LT.0.0) THEN
              WRITE (6,61)  TIME+TOFF
   61         FORMAT (' HYPERB CHAIN CE    T ',F7.1)
          END IF
*
          KSTAR(I1) = KW1
          KSTAR(I2) = KW2
          IPHASE = 0
          IF (ECC.LE.0.001D0) NCE = NCE + 1
*       Include possible NS or BH kick.
          IF (KW1.GE.10.AND.KW1.LE.14) THEN
*       Save chain members to prevent possible over-writing in HIVEL.
              DO 62 L = 1,NCH
                  JPERT(L) = JLIST(L)
   62         CONTINUE
              CALL KICK(I1,1,KW1)
              DO 65 L = 1,NCH
                  JLIST(L) = JPERT(L)
   65         CONTINUE
          END IF
      END IF
*
*       Initialize neighbour force polynomials on significant mass loss.
      IF (DM*SMU.GT.0.1) THEN
*       Save reference body in temporary variables.
          DO 70 K = 1,3
              XREL(K) = X(K,ICH)
              VREL(K) = XDOT(K,ICH)
              XCM(K) = 0.0
              VCM(K) = 0.0
   70     CONTINUE
          SAVEB = BODY(ICH)
*
*       Form new c.m. body.
          ZM = 0.0
          DO 80 L = 1,NCH
              J = JLIST(L)
              IF (J.EQ.0) GO TO 80
              ZM = ZM + BODY(J)
              DO 75 K = 1,3
                  XCM(K) = XCM(K) + BODY(J)*X(K,J)
                  VCM(K) = VCM(K) + BODY(J)*XDOT(K,J)
   75         CONTINUE
   80     CONTINUE
          DO 85 K = 1,3
              X(K,ICH) = XCM(K)/ZM
              XDOT(K,ICH) = VCM(K)/ZM
   85     CONTINUE
          BODY(ICH) = ZM
*
*       Loop over all c.m. neighbours (lists contain old name).
          NNB1 = LIST(1,ICH) + 1
          DO 95 L = 2,NNB1
              J = LIST(L,ICH)
              IF (J.EQ.I1.OR.J.EQ.I3.OR.J.EQ.I4) GO TO 95
              CALL DTCHCK(TIME,STEP(J),DTK(40))
              DO 90 KK = 1,3
                  X0DOT(KK,J) = XDOT(KK,J)
   90         CONTINUE
              CALL FPOLY1(J,J,0)
              CALL FPOLY2(J,J,0)
   95     CONTINUE
*
*       Restore reference body.
          DO 99 K = 1,3
              X(K,ICH) = XREL(K)
              XDOT(K,ICH) = VREL(K)
   99     CONTINUE
          BODY(ICH) = SAVEB
      END IF
*
*       Re-define indices of colliding bodies with J1 as new c.m.
  100 J1 = I1
      J2 = I2
*
      RETURN
*
      END
