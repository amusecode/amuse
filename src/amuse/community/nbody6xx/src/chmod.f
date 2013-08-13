      SUBROUTINE CHMOD(ISUB,KCASE)
*
*
*       Modification of chain member(s).
*       --------------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
     &            NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK,XCM(3),VCM(3),KSCH,Y(NMX8)
      INTEGER  ISORT(NMX)
      LOGICAL  KSLOW2,KCOLL
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/SLOW1/   TK2(0:NMX),EJUMP,KSCH(NMX),KSLOW2,KCOLL
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIK(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/SLOW3/  GCRIT,KZ26
      COMMON/INCOND/  X4(3,NMX),XDOT4(3,NMX)
*
*
*       Identify the dominant perturber (skip if none or NN >= 6).
      ITRY = 0
      JCLOSE = 0
    1 NNB = LISTC(1)
      IF (NNB.EQ.0.OR.NN.GE.6) GO TO 10
      PMAX = 0.0
      DO 2 L = 2,NNB+1
          J = LISTC(L)
          RIJ2 = (X(1,J) - X(1,ICH))**2 + (X(2,J) - X(2,ICH))**2 +
     &                                    (X(3,J) - X(3,ICH))**2
          PIJ = BODY(J)/(RIJ2*SQRT(RIJ2))
          IF (PIJ.GT.PMAX) THEN
              PMAX = PIJ
              RJMIN2 = RIJ2
              JCLOSE = J
          END IF
    2 CONTINUE
*
*       Form the scalar product R*V for sign of radial velocity.
      RDOT = (X(1,JCLOSE) - X(1,ICH))*(XDOT(1,JCLOSE) - XDOT(1,ICH)) +
     &       (X(2,JCLOSE) - X(2,ICH))*(XDOT(2,JCLOSE) - XDOT(2,ICH)) +
     &       (X(3,JCLOSE) - X(3,ICH))*(XDOT(3,JCLOSE) - XDOT(3,ICH))
*
*       Check for rejection (RIJ > 3*MIN(RSUM,RMIN); RDOT > 0 & G < 0.05).
      RIJ = SQRT(RJMIN2)
*
*       Include test on fast escaper approaching a perturber having RDOT > 0.
      IF (GPERT.GT.0.05.AND.RDOT.GT.0) THEN
*       Bypass test and repeated diagnostics for large perturbation.
          IF (GPERT.GT.0.5) GO TO 5
          RDX = 0.0
          L = 0
*       Evaluate actual distance and relative radial velocity of end member.
    3     RJX2 = (XCH(L+1) - X(1,JCLOSE))**2 +
     &           (XCH(L+2) - X(2,JCLOSE))**2 +
     &           (XCH(L+3) - X(3,JCLOSE))**2
          RDI = (XCH(L+1) - X(1,JCLOSE))*(VCH(L+1) - XDOT(1,JCLOSE)) +
     &          (XCH(L+2) - X(2,JCLOSE))*(VCH(L+2) - XDOT(2,JCLOSE)) +  
     &          (XCH(L+3) - X(3,JCLOSE))*(VCH(L+3) - XDOT(3,JCLOSE))
*       See whether any chain member is approaching the intruder.
          IF (RDI.LT.0.0) THEN
              RJX = SQRT(RJX2)
              RDX = RDI/RJX
          END IF
*       Consider the last chain member similarly.
          IF (L.EQ.0) THEN
              L = 3*(NN - 1)
              GO TO 3
          END IF
*       Bypass RDOT < 0 test for approaching ejection candidate.
          IF (RDX.LT.0.0.AND.RJX.LT.2.0*RSUM/FLOAT(NN-1)) THEN
              if(rank.eq.0)
     &        WRITE (6,4)  NAME(JCLOSE), GPERT, RIJ, RJX, RDX
    4         FORMAT (' TRY ABSORB    NAM PERT RIJ RJX RDX ',
     &                                I6,F6.2,1P,4E9.1)
              GO TO 5
          END IF
      END IF
*
*       Include conditions for skipping (large RIJ & size or small GPERT).
      IF (RIJ.GT.3.0*MIN(RSUM,RMIN)) GO TO 10
      IF (RDOT.GT.0.0.AND.GPERT.LT.0.05) GO TO 10
      IF (RSUM.GT.5.0*RMIN.AND.GPERT.LT.1.0) GO TO 10
*       Allow triple hierarchy subject to maximum membership of 6.
    5 IF (NN.GT.3.AND.NAME(JCLOSE).LT.0) GO TO 10
      IF (NN.GT.4.AND.JCLOSE.GT.N) GO TO 10
*
*       Perform impact parameter test: a*(1 - e) < RSUM.
      VR2 = (XDOT(1,ICH) - XDOT(1,JCLOSE))**2 +
     &      (XDOT(2,ICH) - XDOT(2,JCLOSE))**2 +
     &      (XDOT(3,ICH) - XDOT(3,JCLOSE))**2
      AINV = 2.0/RIJ - VR2/(MASS + BODY(JCLOSE))
      ECC2 = (1.0 - RIJ*AINV)**2 + RDOT**2*AINV/(MASS + BODY(JCLOSE))
      ECC = SQRT(ECC2)
      SEMI1 = 1.0/AINV
      PMIN = SEMI1*(1.0 - ECC)
*       Widen the impact parameter test to be on safe side.
      IF (PMIN.GT.1.5*RSUM.AND.GPERT.LT.0.05) GO TO 10
*
*       Check option for increasing or decreasing regularization parameters.
      IF (KZ(16).GT.2) THEN
          IF (GPERT.GT.0.05.AND.NCH.LE.4.AND.KSMAG.LE.5) THEN
              KSMAG = KSMAG + 1
              RMIN = 1.2*RMIN
              RMIN2 = RMIN**2
              DTMIN = 1.3*DTMIN
              IF (KSMAG.GE.4) THEN
                  if(rank.eq.0)
     &            WRITE (77,700)  TIME+TOFF, KSMAG, GPERT, RMIN, RIJ
  700             FORMAT (' INCREASE    T KSMAG GPERT RMIN RIJ  ',
     &                                  F9.1,I4,1P,3E10.2)
                  CALL FLUSH(77)
              END IF
          ELSE IF (GPERT.LT.0.001.AND.KSMAG.GT.1) THEN
              KSMAX = KSMAG - 1
              RMIN = 0.83*RMIN
              RMIN2 = RMIN**2
              DTMIN = 0.77*DTMIN
              if(rank.eq.0)
     &        WRITE (77,705)  TIME+TOFF, KSMAG, GPERT, RMIN, RIJ
  705         FORMAT ('   REDUCE    T KSMAG GPERT RMIN RIJ  ',
     &                              F9.1,I4,1P,3E10.2)
              CALL FLUSH(77)
          END IF
      END IF
*
*       Delay accepting very small binary (suppressed; eccentricity effect).
*     IF (JCLOSE.GT.N.AND.GPERT.LT.0.25) THEN
*         IF (100.0*R(JCLOSE-N).LT.RIJ) GO TO 10
*     END IF
*
*       Include additional criteria (V^2 > VC^2/2; JCL > N & RIJ < RSUM).
      RDOT = RDOT/RIJ
      VC2 = (MASS + BODY(JCLOSE))/RIJ
      IF ((RSUM + RIJ.LT.RMIN).OR.
     &    (RIJ.LT.2.0*RMIN.AND.PMIN.LT.RMIN).OR.
     &    (RIJ.LT.RSUM.AND.RDOT**2.GT.0.5*VC2).OR.
     &    (JCLOSE.GT.N.AND.RIJ.LT.RSUM).OR.
     &    (GPERT.GT.0.2.AND.RDOT.LT.0.0)) THEN
*
*       Do not allow multiple absorptions (NAMES(NMX,5) = 0 each new chain).
          IF (NAME(JCLOSE).NE.NAMES(NMX,5).AND.NAMES(NMX,5).EQ.0) THEN
              NAMES(NMX,5) = NAME(JCLOSE)
*       Include possible return after ejection in bound orbit.
          ELSE IF (RDOT.GT.0.0) THEN
              GO TO 10
          ELSE IF (NN.GT.4.AND.JCLOSE.GT.N.AND.RSUM.GT.2.0*RMIN) THEN
              GO TO 10
          END IF
*
          RSUM = RSUM + RIJ
          IF (KZ(30).GT.1) THEN
              IF (JCLOSE.GT.N) THEN
                  SEMI = -0.5*BODY(JCLOSE)/H(JCLOSE-N)
              ELSE
                  SEMI = 0.0
              END IF
              if(rank.eq.0)
     &        WRITE (6,6)  NSTEP1, JCLOSE, NAME(JCLOSE), GPERT, RIJ,
     &                     RSUM, PMIN, SEMI
    6         FORMAT (' ABSORB:    # JCLOSE NMJ GPERT RIJ RSUM PMIN A ',
     &                             3I6,F6.2,1P,4E9.1)
          END IF
*
*       Switch off any slow-down before increasing the chain (cf. new step).
          IF (KSLOW2) THEN
              FAC = 1.0
              DO 7 K = 1,NCH-1
                  FAC = MAX(KSCH(K),FAC)
    7         CONTINUE
              CALL YCOPY(Y)
              GSAVE = GCRIT
              GCRIT = 0.0
              CALL SLOW
              GCRIT = GSAVE
          ELSE
              FAC = 1.0
          END IF
*
*       Treat extreme binary as inert on significant slow-down.
          IF (JCLOSE.GT.N.AND.FAC.GT.10.0) THEN
              SEMI = -0.5*BODY(JCLOSE)/H(JCLOSE-N)
              RM = 1.0
              DO 8 K = 1,NCH-1
                  RM = MIN(1.0/RINV(K),RM)
    8         CONTINUE
              SMALL = 0.01*(RSUM - RIJ)
              IF (rank.eq.0.and.RM.LT.SMALL.AND.SEMI.LT.SMALL) THEN
                  WRITE (6,9)  NAME(JCLOSE), GPERT, SEMI, RIJ,
     &                                       (1.0/RINV(K),K=1,NCH-1)
    9             FORMAT (' ABSORB INERT    NM G A RIJ R ',
     &                                      I6,F6.2,1P,6E10.2)
*       Switch off option #26 temporarily for routine SETSYS.
                  KZ(26) = 1
                  CALL ABSORB(ISUB)
                  KZ(26) = 2
                  KCASE = 1
                  GO TO 50
              END IF
          END IF
*
*       Absorb the perturber (single particle or binary).
          CALL ABSORB(ISUB)
*
*       Reduce block-time since new c.m. step may be very small.
          TBLOCK = MIN(TIME,TBLOCK)
*
*       Activate indicator for new chain treatment and try a second search.
          KCASE = 1
          ITRY = ITRY + 1
          IF (ITRY.EQ.1) GO TO 1
      END IF
*
*       Exit if any particles have been absorbed.
   10 IF (ITRY.GT.0) GO TO 50
*
*       Place index of the smallest INVERSE distance in ISORT(1).
      CALL HPSORT(NN-1,RINV,ISORT)
*
*       Determine index of escaper candidate (single star or binary).
      KCASE = 0
      JESC = 0
*       Distinguish two cases of each type (beginning & end of chain).
      IF (ISORT(1).EQ.1) THEN
          IESC = INAME(1)
          KCASE = 1
      ELSE IF (ISORT(1).EQ.NN-1) THEN
          IESC = INAME(NN)
          KCASE = 1
*       Check for possible binary escaper (NN = 3 implies single escaper).
      ELSE IF (ISORT(1).EQ.2) THEN
          IESC = INAME(1)
          JESC = INAME(2)
          IBIN = 1
*       Switch binary index in case last separation is smallest (NN = 4).
          IF (RINV(1).LT.RINV(NN-1)) IBIN = NN - 1
          KCASE = 2
      ELSE IF (ISORT(1).EQ.NN-2) THEN
          IESC = INAME(NN-1)
          JESC = INAME(NN)
          IBIN = NN - 1
          KCASE = 2
      END IF
*
*       Try escape check if middle distance is largest.
      IF (KCASE.EQ.0.AND.NN.GE.5.AND.
     &    1.0/RINV(ISORT(1)).GT.2.0*RMIN) THEN
          R1 = 1.0/RINV(1)
          R2 = 1.0/RINV(NN-1)
*       Set relevant indices for beginning or end of chain.
          IF (R1.GT.R2) THEN
              IESC = INAME(1)
              JX = INAME(2)
              IB = 1
              ISORT(1) = 1
              R3 = 1.0/RINV(2)
          ELSE
              IESC = INAME(NN)
              JX = INAME(NN-1)
              IB = NN - 1
              ISORT(1) = NN - 1
              R3 = 1.0/RINV(NN-2)
          END IF
*       Define binary indices for large second separation.
          IF (R3.GT.MAX(R1,R2)) THEN
              JESC = JX
              IBIN = IB
              KCASE = 2
          ELSE
              KCASE = 1
          END IF
      END IF
*
*       Include safety test for abnormal configuration.
      IF (KCASE.EQ.0) GO TO 60
*
      IF (rank.eq.0.and.KZ(30).GT.2) THEN
          WRITE (6,12)  IESC, JESC, NSTEP1, ISORT(1),
     &                  (1.0/RINV(K),K=1,NN-1)
   12     FORMAT (' CHMOD:    IESC JESC # ISORT1 R ',2I3,I6,I3,1P,5E9.1)
      END IF
*
*       Copy chain variables to standard form.
      LK = 0
      DO 20 L = 1,NCH
          DO 15 K = 1,3
              LK = LK + 1
              X4(K,L) = XCH(LK)
              XDOT4(K,L) = VCH(LK)
   15     CONTINUE
   20 CONTINUE
*
*       First check case of escaping binary (JESC > 0 & RB < RJB/4).
      IF (JESC.GT.0) THEN
          RB = 1.0/RINV(IBIN)
          JBIN = IBIN + 1
          IF (IBIN.EQ.NN - 1) JBIN = IBIN - 1
          RJB = 1.0/RINV(JBIN)
*       Consider removal of outermost particle instead if binary is wide.
          IF (RB.GT.0.25*RJB) THEN
*       Change pointer to end of chain and redefine IESC.
              IF (ISORT(1).EQ.2) THEN
                  ISORT(1) = 1
              ELSE
                  ISORT(1) = NN - 1
              END IF
              IESC = JESC
              GO TO 30
          END IF
      ELSE
          GO TO 30
      END IF
*
*       Form coordinates & velocities of escaping binary (local c.m. frame).
      BCM = BODYC(IESC) + BODYC(JESC)
      RI2 = 0.0
      RDOT = 0.0
      VREL2 = 0.0
      DO 25 K = 1,3
          XCM(K) = (BODYC(IESC)*X4(K,IESC) + BODYC(JESC)*X4(K,JESC))/BCM
          VCM(K) = (BODYC(IESC)*XDOT4(K,IESC) +
     &              BODYC(JESC)*XDOT4(K,JESC))/BCM
          RI2 = RI2 + XCM(K)**2
          RDOT = RDOT + XCM(K)*VCM(K)
          VREL2 = VREL2 + (XDOT4(K,IESC) - XDOT4(K,JESC))**2
   25 CONTINUE
*
*       Convert to relative distance & radial velocity w.r.t. inner part.
      FAC = BODY(ICH)/(BODY(ICH) - BCM)
      RI = SQRT(RI2)
      RDOT = FAC*RDOT/RI
      RI = FAC*RI
*       Adopt arithmetic mean of RSUM and RMAXS for delaying escape.
      DESC = 0.5*SQRT(RSUM*RMAXS(ISUB))
*       Reduce delay test for large length contrasts.
      IM = ISORT(1)
      RM = 1.0/RINV(IM)
      CX = RSUM/(RSUM - RM)
      IF (CX.GT.25.0) DESC = 25.0*RGRAV
      RESC = MAX(3.0*RGRAV,DESC)
      AINV = 2.0/RB - VREL2/BCM
      EB = -0.5*BODYC(IESC)*BODYC(JESC)*AINV
      HI = 0.5*RDOT**2 - BODY(ICH)/RI
*
*       Employ parabolic escape criterion (terminate if RI > RMIN & NCH < 5).
      IF ((RI.GT.RESC.AND.RDOT.GT.0.0.AND.RB*AINV.GT.0.99).OR.
     &    (HI.GT.0.0.AND.RI.GT.3.0*RMIN)) THEN
*       Decide between weakly bound and escape orbit (RMAXS may be large).
          IF (RDOT**2.LT.2.0*BODY(ICH)/RI) THEN
*       Define effective perturbation using remaining size.
              GB = 2.0*BCM/(BODY(ICH) - BCM)*((RSUM - RJB - RB)/RI)**3
*       Split into two KS solutions if binaries are well separated.
              IF (NCH.EQ.4.AND.GB.LT.0.001) THEN
                  KCASE = -1
                  GO TO 50
              END IF
*       Enforce termination if RI > MIN(RMAXS/2,RMIN) and NCH <= 4.
              IF (RI.GT.MIN(0.5*RMAXS(ISUB),RMIN).AND.NCH.LE.4) THEN
                  KCASE = -1
                  GO TO 50
*       Accept binary escape for small perturbation & V**2 > M/R (NCH > 4).
              ELSE IF (GB.LT.0.01.AND.NCH.GT.4.AND.
     &                 (RDOT**2.GT.BODY(ICH)/RI.OR.RI.GT.RMIN)) THEN
                  if(rank.eq.0)
     &            WRITE (6,28)  IESC, JESC, NAMEC(IESC), NAMEC(JESC),
     &                          RI, RDOT**2, 2.0*BODY(ICH)/RI, RB
                  CM(9) = CM(9) - EB
                  GO TO 40
              ELSE
*       Check splitting into two KS solutions (R1 + R3 < R2/5 & RDOT > VC).
                  IF (NCH.EQ.4) THEN
                      R13 = 1.0/RINV(1) + 1.0/RINV(3)
                      VC2 = RDOT**2*RI/BODY(ICH)
*       Ensure RSUM > 0.1*RMIN for reducing differential force corrections.
                      IF (R13.LT.0.2/RINV(2).AND.VC2.GT.1.0.AND.
     &                    RSUM.GT.0.1*RMIN) THEN
                          KCASE = -1
                          GO TO 50
                      END IF
                  END IF
                  KCASE = 0
                  GO TO 60
              END IF
          ELSE
              IF (HI.GT.0.0) THEN
                  VINF = SQRT(2.0*HI)*VSTAR
              ELSE
                  VINF = 0.0
              END IF
              IF (KZ(30).GT.1.OR.VINF.GT.1.0) THEN
                  if(rank.eq.0)
     &            WRITE (6,28)  IESC, JESC, NAMEC(IESC), NAMEC(JESC),
     &                          RI, RDOT**2, 2.0*BODY(ICH)/RI, RB, VINF
   28             FORMAT (' CHAIN ESCAPE:    IESC JESC NM RI RDOT2 ',
     &                                      '2*M/R RB VINF ',
     &                                       2I3,2I6,1P,4E9.1,0P,F6.1)
              END IF
*       Enforce termination (KCASE < 0) if NCH <= 4 (final membership <= 2).
              IF (NCH.LE.4) THEN
                  KCASE = -1
                  GO TO 50
              END IF
              CM(9) = CM(9) - EB
              GO TO 40
          END IF
      ELSE
          KCASE = 0
          GO TO 60
      END IF
*
*       Form relative distance and radial velocity for single particle.
   30 RI = SQRT(X4(1,IESC)**2 + X4(2,IESC)**2 + X4(3,IESC)**2)
      RDOT = X4(1,IESC)*XDOT4(1,IESC) + X4(2,IESC)*XDOT4(2,IESC) +
     &                                  X4(3,IESC)*XDOT4(3,IESC)
      FAC = BODY(ICH)/(BODY(ICH) - BODYC(IESC))
      RDOT = FAC*RDOT/RI
      RI = FAC*RI
      IM = ISORT(1)
*       Ensure that escaper is not close to another body.
      RM = MIN(1.0/RINV(IM),RI)
      RB = 0.0
*
*       Check approximate escape criterion outside 3*RGRAV and DESC.
      DESC = 0.5*SQRT(RSUM*RMAXS(ISUB))
*       Reduce delay test for large length contrasts.
      CX = RSUM/(RSUM - RM)
      IF (CX.GT.25.0) DESC = 25.0*RGRAV
*       Note that arithmetic mean tends to delay escape.
      RESC = MAX(5.0*RGRAV,DESC)
      RESC = MIN(RESC,3.0*RMIN)
      IF (NN.EQ.3.AND.RESC.LT.0.001*RMIN) RESC = 5.0*RESC
      IF (RM.GT.RESC.AND.RDOT.GT.0.0) THEN
          IF (RDOT**2.LT.2.0*BODY(ICH)/RI) THEN
              FAC2 = 2.0*BODY(IESC)/(BODY(ICH) - BODYC(IESC))
              GI = FAC2*((RSUM - 1.0/RINV(IM))/RI)**3
*       Do not permit large length contrast even for bound orbit (NN = 4).
              IF (NN.EQ.4.AND.GI.LT.0.01) THEN
                  IB = ISORT(NN-1)
*       Prevent NN = 4 termination near small pericentre (subject to safety).
                  IF (1.0/RINV(IB).GT.0.1*RGRAV.OR.RSUM.GT.RMIN) THEN
                      KCASE = -1
                      GO TO 50
                  ELSE
                      KCASE = 0
                      GO TO 60
                  END IF
              END IF
              ER = 0.5*RDOT**2 - BODY(ICH)/RI
              RX = -BODY(ICH)/ER
              IF ((ER.LT.0.0.AND.RX.LT.MIN(2.0*RSUM,2.0*RMIN)).OR.
     &            (NN.EQ.3.AND.GI.GT.0.1)) THEN
                  KCASE = 0
                  GO TO 60
              END IF
              IF (RI.GT.MIN(0.5*RMAXS(ISUB),RMIN)) THEN
*       Decide between termination or removal for sub-parabolic motion.
                  IF (NN.LT.4) THEN
*       Delay termination for eccentric binary inside R = -m1*m2/2*E < a.
                      IB = 3 - IM
                      J1 = INAME(IB)
                      J2 = INAME(IB+1)
*       Obtain semi-major axis from ENERGY = -(m1*m2 + (m1+m2)*m3)/RGRAV.
                      ZMU = BODYC(J1)*BODYC(J2)/(BODYC(J1) + BODYC(J2))
                      ZF = BODYC(IESC)/ZMU
*       Write E = E_b + E_1 with E_b = -m1*m2/2*a and E_1 = MU3*ER.
                      SEMI = 0.5*RGRAV/(1.0 + ZF*(1.0 + RGRAV*ER/MASS))
                      RY = 1.0/RINV(IB)
                      IF (RY.LT.0.9*SEMI.AND.RI.LT.2.0*RMIN) THEN
                          if(rank.eq.0)
     &                    WRITE (6,31)  NSTEP1, RY/SEMI, RI, RDOT**2,
     &                                  2.0*BODY(ICH)/RI, SEMI
                          KCASE = 0
                          GO TO 60
                      END IF
                      KCASE = -1
                      GO TO 60
                  ELSE
                      KCASE = 1
                      GO TO 32
                  END IF
              ELSE
                  KCASE = 0
                  GO TO 60
              END IF
          ELSE IF (NN.GE.4) THEN
*       Accept escape of hyperbolic body if first or last separation > RMIN.
              RBIG = 1.0/RINV(IM)
              IF (RBIG.GT.RMIN) THEN
                  IF (IM.EQ.1.OR.IM.EQ.NN-1) THEN
                      KCASE = 1
                      GO TO 32
                  END IF
              END IF
*       Include three-body stability test for distant fourth body.
              IF (RBIG.GT.0.75*RSUM.AND.NN.EQ.4.AND.
     &            RSUM.GT.0.5*RMIN) THEN
                  CALL CHSTAB(ITERM)
                  IF (ITERM.LT.0) THEN
                      KCASE = -1
                      GO TO 60
                  END IF
              END IF
*       Skip if largest separation is not dominant.
              IF (RBIG.LT.0.66*RSUM) THEN
                  KCASE = 0
                  GO TO 60
              END IF
          ELSE
*       Note possibility a << RGRAV/6 after strong interaction with NN = 3.
              IB = 3 - IM
              J1 = INAME(IB)
              J2 = INAME(IB+1)
*       Obtain semi-major axis from ENERGY = -(m1*m2 + (m1+m2)*m3)/RGRAV.
              ZMU = BODYC(J1)*BODYC(J2)/(BODYC(J1) + BODYC(J2))
              ZF = BODYC(IESC)/ZMU
              ER = 0.5*RDOT**2 - BODY(ICH)/RI
*       Write E = E_b + E_1 with E_b = -m1*m2/2*a and E_1 = MU3*ER.
              SEMI = 0.5*RGRAV/(1.0 + ZF*(1.0 + RGRAV*ER/MASS))
              RY = 1.0/RINV(IB)
*       Note RGRAV = MIN(0.5*RSUM,RGRAV) is OK and RGRAV adjusted on QPMOD.
              IF (RY.LT.0.9*SEMI.AND.RI.LT.2.0*RMIN.AND.
     &            RI.LT.50.0*SEMI) THEN
                  if(rank.eq.0)
     &            WRITE (6,31)  NSTEP1, RY/SEMI, RI, RDOT**2,
     &                          2.0*BODY(ICH)/RI, SEMI
   31             FORMAT (' CHAIN DELAY    # R/A RI RD2 VP2 A ',
     &                                     I5,F6.2,1P,4E9.1)
                  KCASE = 0
                  GO TO 60
              END IF
          END IF
*
*       Check that escaper is well separated (ratio > 3).
          RR = RSUM - 1.0/RINV(IM)
          RATIO = 1.0/(RINV(IM)*RR)
          IF (RATIO.LT.3.0.AND.RSUM.LT.RMIN) THEN
              KCASE = 0
              GO TO 60
          END IF
*
   32     HI = 0.5*RDOT**2 - BODY(ICH)/RI
          IF (HI.GT.0.0) THEN
              VINF = SQRT(2.0*HI)*VSTAR
          ELSE
              VINF = 0.0
          END IF
          IF (KZ(30).GT.1.OR.VINF.GT.2.0) THEN
              if(rank.eq.0)
     &        WRITE (6,35)  IESC, NAMEC(IESC), RI, RDOT**2,
     &                      2.0*BODY(ICH)/RI, VINF
   35         FORMAT (' CHAIN ESCAPE:    IESC NM RI RDOT2 2*M/R VF ',
     &                                   I3,I6,1P,3E9.1,0P,F6.1)
          END IF
*       Ensure single body is removed in case of wide binary.
          JESC = 0
      ELSE
          KCASE = 0
          GO TO 60
      END IF
*
*       Reduce chain membership (NCH > 3) or specify termination.
   40 IF (NCH.GT.3) THEN
*       Subtract largest chain distance from system size (also binary).
          IM = ISORT(1)
          RSUM = RSUM - 1.0/RINV(IM) - RB
          CALL REDUCE(IESC,JESC,ISUB)
      ELSE
          KCASE = -1
      END IF
*
*       Set phase indicator < 0 to ensure new time-step list in INTGRT.
   50 IPHASE = -1
*
   60 RETURN
*
      END
