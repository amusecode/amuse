      SUBROUTINE IMPACT(I)
*
*
*       Multiple collision or merger search.
*       ------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      CHARACTER*8  WHICH1
      REAL*8  XX(3,3),VV(3,3)
      INTEGER LISTQ(100)
      SAVE LISTQ,QCHECK
      DATA IZARE,LISTQ(1),QCHECK /0,0,0.0D0/
*
*
*       Set index of KS pair & both components of c.m. body #I.
      IPAIR = I - N
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      NTTRY = NTTRY + 1
      PERT1 = 0.0
      PERT2 = 0.0
      JCOMP = IFIRST
      NP = 0
      KS2 = 0
      RMAX2 = 1.0
      TTOT = TIME + TOFF
      RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                               (X(3,I) - RDENS(3))**2
*
*       Search c.m. neighbours if binary has at most two perturbers.
      J1 = I1
      IF (LIST(1,J1).LE.2) J1 = I
      NNB2 = LIST(1,J1) + 1
*
*       Find the dominant body (JCOMP) and nearest perturber (JMAX).
      DO 10 L = 2,NNB2
          J = LIST(L,J1)
          RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                  (X(3,I) - X(3,J))**2
          IF (J1.EQ.I.AND.RIJ2.GT.1.0D+04*RMIN2) GO TO 10
          NP = NP + 1
          JLIST(NP) = J
          PERT = BODY(J)/(RIJ2*SQRT(RIJ2))
          IF (PERT.GT.PERT2) THEN 
              IF (PERT.GT.PERT1) THEN
                  RJMIN2 = RIJ2
                  JMAX = JCOMP
                  JCOMP = J
                  PERT2 = PERT1
                  PERT1 = PERT
              ELSE
                  JMAX = J
                  PERT2 = PERT
                  RMAX2 = RIJ2
              END IF
          END IF
   10 CONTINUE
*
*       Include safety check on rare case of no perturbers inside 100*RMIN.
      IF (NP.EQ.0) THEN
          JCOMP = LIST(2,I)
          JMAX = LIST(3,I)
          RJMIN2 = 1.0
      ELSE IF (NP.EQ.1) THEN
          JMAX = N
      END IF
*
      RDOT = (X(1,I) - X(1,JCOMP))*(XDOT(1,I) - XDOT(1,JCOMP)) +
     &       (X(2,I) - X(2,JCOMP))*(XDOT(2,I) - XDOT(2,JCOMP)) +
     &       (X(3,I) - X(3,JCOMP))*(XDOT(3,I) - XDOT(3,JCOMP))
*
*       Specify larger perturbation for optional chain regularization.
      IF ((KZ(30).GT.0.OR.KZ(30).EQ.-1).AND.NCH.EQ.0) THEN
          GSTAR = 100.0*GMIN
          KCHAIN = 1
      ELSE
          GSTAR = GMIN
          KCHAIN = 0
*       Specify indicator -1 for allowing TRIPLE & QUAD but not CHAIN.
          IF (KZ(30).EQ.-2) KCHAIN = -1
      END IF
*
*       Only accept inward motion or small secondary perturbation.
      PERT3 = 2.0*R(IPAIR)**3*PERT2/BODY(I)
      IF (RDOT.GT.0.0.OR.PERT3.GT.100.0*GSTAR) GO TO 100
*
*       Include impact parameter test to distinguish different cases.
      A2 = (XDOT(1,I) - XDOT(1,JCOMP))**2 + 
     &     (XDOT(2,I) - XDOT(2,JCOMP))**2 +
     &     (XDOT(3,I) - XDOT(3,JCOMP))**2
      RIJ = SQRT(RJMIN2)
      A3 = 2.0/RIJ - A2/(BODY(I) + BODY(JCOMP))
      SEMI1 = 1.0/A3
      A4 = RDOT**2/(SEMI1*(BODY(I) + BODY(JCOMP)))
      ECC1 = SQRT((1.0D0 - RIJ/SEMI1)**2 + A4)
      PMIN = SEMI1*(1.0D0 - ECC1)
*
*       Set semi-major axis, eccentricity & apocentre of inner binary.
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      A0 = SEMI
      ECC2 = (1.0D0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      ECC = SQRT(ECC2)
      APO = ABS(SEMI)*(1.0 + ECC)
*
*       Quit on hyperbolic orbit with large impact parameter.
      IF (ECC1.GT.1.0.AND.PMIN.GT.50.0*SEMI) GO TO 100
*
*       Include rectification for non-circular binaries with KSTAR = 10 & 12.
      IF (KZ(27).EQ.2.AND.KSTAR(I).GE.10.AND.KSTAR(I).LE.18) THEN
          IF (ECC.GT.0.1.AND.MOD(KSTAR(I),2).EQ.0) THEN
              RM = MAX(RADIUS(I1),RADIUS(I2))
              ICIRC = -1
              CALL INDUCE(IPAIR,EMAX,EMIN,ICIRC,TC,ANGLE,TG,EDAV)
              if(rank.eq.0)
     &        WRITE (6,15)  NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2),
     &                      KSTAR(I), LIST(1,I1), ECC, SEMI, RM, PMIN,
     &                      GAMMA(IPAIR), TC, 360.0*ANGLE/TWOPI, EMAX
   15         FORMAT (' NON-CIRCULAR    NM K* NP E A R* PM G TC IN EX ',
     &                                  2I7,4I4,F7.3,1P,5E9.1,0P,2F7.2)
*       Circularize the orbit instantaneously for short TC.
              IF (TC.LT.-100.0) THEN
*       Set temporary unperturbed orbit (avoids new KSPOLY in DEFORM).
                  NP = LIST(1,I1)
                  LIST(1,I1) = 0
                  DT1 = STEP(I1)
                  TIME0 = TIME
                  CALL KSRECT(IPAIR)
                  QP = SEMI*(1.0 - ECC)
                  ERR = ABS(QP - R(IPAIR))/QP
*       Deform orbit to circular eccentricity after determining apocentre.
                  IF (R(IPAIR).GT.SEMI.OR.ERR.GT.1.0D-04) THEN
*       Reduce eccentric anomaly by pi for inward motion.
                      IF (TDOT2(IPAIR).LT.0.0D0) THEN
                          CALL KSAPO(IPAIR)
                      END IF
                  END IF
*       Predict to pericentre and transform by pi to exact apocentre.
                  CALL KSPERI(IPAIR)
                  CALL KSAPO(IPAIR)
                  ECCM = 0.002
                  CALL DEFORM(IPAIR,ECC,ECCM)
                  LIST(1,I1) = NP
*        Resolv X & XDOT and initialize KS polynomial at apocentre time.
                  IF (NP.EQ.0) DT1 = 0.0
                  TIME = TIME0 - DT1
                  CALL RESOLV(IPAIR,1)
                  CALL KSPOLY(IPAIR,1)
              ELSE
                  KSTAR(I) = 0
              END IF
          END IF
      END IF
*
*       Form binding energy of inner & outer binary.
      EB = BODY(I1)*BODY(I2)*H(IPAIR)/BODY(I)
      IF(ABS(EB).LT.1.0D-10) EB = -1.0D-10
      EB1 = -0.5*BODY(JCOMP)*BODY(I)/SEMI1
*
*       Obtain the total perturbing force acting on body #I & JCOMP.
      CALL FPERT(I,JCOMP,NP,PERT)
*
*       Choose maximum of dominant scalar & total vectorial perturbation.
      PERT = PERT*RJMIN2/(BODY(I) + BODY(JCOMP))
      PERT4 = 2.0*RJMIN2*RIJ*PERT2/(BODY(I) + BODY(JCOMP))
      PERTM = MAX(PERT4,PERT)
*
*       Use combined semi-major axis for binary-binary collision.
      IF (JCOMP.GT.N) THEN
          JPAIR = JCOMP - N
          SEMI2 = -0.5D0*BODY(JCOMP)/H(JPAIR)
          J1 = 2*JPAIR - 1
          EB2 = -0.5*BODY(J1)*BODY(J1+1)/SEMI2
*       Define SEMI0 as smallest binary in case IPAIR denotes widest pair.
          SEMI0 = MIN(ABS(SEMI),ABS(SEMI2))
          SEMIX = MAX(SEMI,SEMI2)
          APO = APO + MAX(ABS(SEMI2),R(JPAIR))
          SEMI = SEMI + SEMI2
*       Do not allow negative or soft cross section.
          IF (1.0/SEMI.LT.0.5/RMIN) GO TO 100
*       Consider merger for PMIN > SEMI and large semi-major axis ratio.
          IF (PMIN.GT.SEMI.AND.SEMI2.GT.20.0*SEMI0) GO TO 30
      END IF
*
*       Check separation in case of chain regularization.
      IF (KCHAIN.GT.0) THEN
*       Form effective gravitational radius (combine triple & quad).
          EBT = EB + EB1
          ZMM = BODY(I1)*BODY(I2) + BODY(I)*BODY(JCOMP)
*       Set length of chain for decision-making (also used at termination).
          RSUM = R(IPAIR) + RIJ
          RI = R(IPAIR)
          IF (JCOMP.GT.N) THEN
              EBT = EBT + EB2
              ZMM = ZMM + BODY(J1)*BODY(J1+1)
              RSUM = RSUM + R(JPAIR)
              RI = MAX(R(JPAIR),RI)
          END IF
          RGRAV = ZMM/ABS(EBT)
*       Employ initial size as upper limit in case of weakly bound system.
          RGRAV = MIN(RGRAV,RMIN)
*       Save initial energy in binaries for routine SETSYS.
          EBCH0 = EBT - EB1
*       Use RIJ instead of RSUM in 3*RGRAV test (increases initial RIJ).
          IF (RIJ.GT.MAX(3.0*RGRAV,RMIN).OR.RSUM.GT.2.0*RMIN) GO TO 30
          GI = 2.0*BODY(JCOMP)*(RI/RIJ)**3/BODY(I)
*       Enforce KS orbit using MERGE for high eccentricity if PMIN > 10*RI.
          IF (ECC1.GT.0.99.AND.PMIN.GT.10.0*RI.AND.
     &        PERTM.LT.GMAX) GO TO 40
          IF (KZ(27).GT.0.AND.JCOMP.GT.N) THEN
              IF (SEMI0.LT.SEMI2) J1 = I1
              RT = 4.0*MAX(RADIUS(J1),RADIUS(J1+1))
*       Do not allow large distance ratio for nearly synchronous binary.
              IF (SEMI0.GT.RT.AND.RI.GT.25.0*SEMI0) GO TO 30
              IF (MIN(SEMI0,SEMI2).LT.0.05*RIJ) THEN
              IF (MAX(SEMI0,SEMI2).LT.0.1*RIJ) GO TO 30
              END IF
          END IF
      END IF
*
*       Include special case of strong interraction and large ECC1.
      IF (ECC1.GT.0.9.AND.GAMMA(IPAIR).GT.0.01) THEN
          IF (APO.LT.0.01*RMIN.AND.PMIN.LT.2.5*APO) GO TO 16
      END IF
*
*       Adopt triple, quad or chain regularization for strong interactions.
*     IF ((APO.GT.0.01*RMIN.OR.JCOMP.GT.N).AND.PMIN.GT.APO) GO TO 30
      IF ((APO.GT.0.01*RMIN.OR.JCOMP.GT.N).AND.PMIN.GT.1.5*APO) GO TO 30
*     IF (APO.LE.0.01*RMIN.AND.PMIN.GT.2.0*APO) GO TO 30
      IF ((RIJ.GT.RMIN.AND.SEMI1.GT.0.0).OR.RIJ.GT.2.0*RMIN) GO TO 100
      IF (PERTM.GT.100.0*GSTAR) GO TO 30
   16 IF (JCOMP.GT.N.AND.PMIN.GT.0.1*RMIN) THEN
          IF (PMIN.GT.A0 + SEMI2) GO TO 30
      END IF
      IF (JCOMP.GT.N.AND.PMIN.GT.4.0*SEMIX.AND.
     &   (ECC1.GT.0.9.AND.ECC1.LT.1.0)) GO TO 30
*
*       Check almost stable triples (factor 1.2 is experimental).
      IF (JCOMP.LE.N.AND.PMIN.GT.2.5*SEMI) THEN
          CALL HISTAB(IPAIR,JCOMP,PMIN,RSTAB)
          RA = SEMI1*(1.0 + ECC1)
          IF (SEMI1.LT.0.0) RA = RIJ
          GI = PERT*(RA/RIJ)**3
*       Use estimated apocentre perturbation for decision-making.
          IF (PMIN.GT.1.2*RSTAB) THEN
              IF (GI.LT.0.05) GO TO 30
*       Choose chain for critical case of highly eccentric outer orbit.
              IF (ECC1.LT.0.95) GO TO 100
          ELSE IF (PMIN.GT.0.7*RSTAB) THEN
*       Treat marginally stable triple according to external perturbation.
              IF (GI.LT.0.05) GO TO 30
              IF (GI.LT.1.0.OR.ECC1.LT.0.9) GO TO 100
          END IF
          IF (PMIN.GT.0.6*RSTAB.AND.PMIN.LT.0.9*RSTAB) GO TO 100
*       Delay for large distance ratio outside 0.5*RMIN.
          IF (RIJ.GT.MAX(10.0*APO,0.5*RMIN)) GO TO 100
          IF (RIJ.GT.10.0*APO) GO TO 100
          IF (PMIN.GT.2.5*APO) GO TO 40
      END IF
      IF (PMIN.GT.3.0*SEMI.AND.JCOMP.LE.N) GO TO 40
*
      IF (JCOMP.GT.N) THEN
          IF (RIJ.GT.10.0*APO) GO TO 100
      END IF
*       Skip chain if merged binary or chain c.m. (denoted by NAME <= 0).
      IF (NAME(I).LE.0.OR.NAME(JCOMP).LE.0) GO TO 100
*
*       Compare with existing subsystem of same type (if any).
      IF (NSUB.GT.0.AND.KCHAIN.LE.0) THEN
          PERIM = R(IPAIR) + RIJ
          IF (JCOMP.GT.N) PERIM = PERIM + R(JPAIR)
          IGO = 0
          CALL PERMIT(PERIM,IGO)
          IF (IGO.GT.0) GO TO 100
      END IF
*
*       Skip all multiple regs on zero option (mergers done by #15 > 0).
      IF (KZ(30).EQ.0) GO TO 100
*       Allow CHAIN only (#30 = -1) or TRIPLE & QUAD only (#30 < -1).
      IF (KZ(30).EQ.-2.AND.KCHAIN.EQ.0) GO TO 100
*
      WHICH1 = ' TRIPLE '
      IF (JCOMP.GT.N) WHICH1 = ' QUAD   '
      IF (KCHAIN.GT.0) WHICH1 = ' CHAIN  '
*
      IF (rank.eq.0.and.H(IPAIR).GT.0.0) THEN
          WRITE (6,18)  I, JCOMP, ECC, ECC1, SEMI1, RIJ, GAMMA(IPAIR)
   18     FORMAT (' HYP CHAIN    I J E E1 A1 RIJ G  ',
     &                           2I6,2F7.3,1P,3E9.1)
      END IF
*
      IF (KZ(15).GT.1.OR.KZ(30).GT.1) THEN
          RI = SQRT((X(1,I) - RDENS(1))**2 +
     &              (X(2,I) - RDENS(2))**2 +
     &              (X(3,I) - RDENS(3))**2)
          VI = SQRT(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2)
          if(rank.eq.0)
     &    WRITE (6,20)  WHICH1, TTOT, I, JCOMP, IPAIR, NAME(I1),
     &        NAME(I2), NAME(JCOMP), NAME(I), H(IPAIR), R(IPAIR),
     &                  BODY(I), BODY(JCOMP), PERT4, RIJ, PMIN,
     &                  EB1/EB, RI, VI, LIST(1,I1)
   20     FORMAT (/,' NEW',A8,1P,' T =',E12.5,' I =',3I8,' N =',4I8,0P,
     &              ' H =',F10.2,' R =',1P,E12.5,' M =',2E12.5,
     &              ' G4 =',E12.5,' R1 =',E8.1,' P =',E8.1,
     &              ' E1 =',E8.1,' RI, VI =',2E12.5,' NP =',I4)
      END IF
*
*       Include any close single or c.m. perturber (cf. routine SETSYS).
      IF (JMAX.NE.JCOMP.AND.SQRT(RMAX2).LT.MIN(2.0D0*RSUM,RMIN).AND.
     &    NAME(JMAX).GT.0) THEN
          IF (JCOMP.GT.N.AND.JMAX.GT.N) THEN
              JCMAX = 0
          ELSE
              if(rank.eq.0)
     &        WRITE (6,21)  NAME(JCOMP), NAME(JMAX), RSUM, SQRT(RMAX2)
   21         FORMAT (' B+2 CHAIN    NAM RSUM RMX ',2I7,1P,2E10.2)
              CALL XVPRED(JMAX,-1)
              JCMAX = JMAX
          END IF
      ELSE
          JCMAX = 0
      END IF
*
*       Save global index of intruder for TRIPLE or CHAIN.
      JCLOSE = JCOMP
*
*       Check B-B interaction for switch of IPAIR & JPAIR or inert binary.
      IF (KCHAIN.GT.0.AND.JCOMP.GT.N) THEN
          K1 = 2*JPAIR - 1
          if(rank.eq.0)
     &    WRITE (6,22)  NAME(I1), NAME(I2), NAME(K1), NAME(K1+1),
     &                  KSTAR(I), KSTAR(JCOMP), ECC, ECC1, A0, SEMI2,
     &                  RIJ, SEMI1, PMIN
   22     FORMAT (' CHAIN B-B    NAME K* E0 E1 A0 A2 RIJ A1 PM ',
     &                           4I7,2I4,2F7.3,1P,5E10.2)
          RT = 4.0*MAX(RADIUS(I1),RADIUS(I2))
          IF (SEMI0.LT.4.0*RT.AND.LIST(1,J1).EQ.0.OR.
     &        MIN(SEMI0,SEMI2).LT.0.01*RIJ) THEN
*       Ensure that widest binary comes first (more similar to triple).
              IF (SEMI0.LT.SEMI2) THEN
                  KPAIR = JPAIR
                  JPAIR = IPAIR
                  IPAIR = KPAIR
                  JCLOSE = N + JPAIR
              END IF
*       Check reduction of c.m. index (JPAIR becomes JPAIR - 1 if > IPAIR).
              IF (JPAIR.GT.IPAIR) JCLOSE = JCLOSE - 1
              IF (KZ(26).LT.2) THEN
*       Replace unperturbed near-synchronous binary by inert body in CHAIN.
                  JCOMP = 0
                  if(rank.eq.0)
     &            WRITE (6,25)  SEMI0, RIJ, R(JPAIR), GAMMA(JPAIR)
   25             FORMAT (' INERT BINARY    A RIJ R G ',1P,4E10.2)
              END IF
          ELSE
              JCLOSE = 0
          END IF
      END IF
*
*       Set phase indicator for calling TRIPLE or QUAD from MAIN.
      IPHASE = 4
      KSPAIR = IPAIR
*
*       Include the case of two interacting KS pairs.
      IF (JCOMP.GT.N) THEN
          IPHASE = 5
*       Switch pair indices and rename JCOMP if JPAIR has smaller size.
          IF (STEP(J1).LT.STEP(I1).AND.LIST(1,I1).GT.0) THEN
              KSPAIR = JPAIR
              JCOMP = I
              KS2 = IPAIR
          ELSE
              KS2 = JPAIR
          END IF
          IF (KZ(27).LE.0.AND.JPAIR.GT.IPAIR) THEN
              IF (JCLOSE.GT.0) JCLOSE = JCLOSE - 1
          END IF
*       Terminate smallest pair first and reduce second index if higher.
*         CALL KSTERM
          IF (KS2.GT.KSPAIR) KS2 = KS2 - 1
      END IF
*
*       See whether chain regularization indicator should be switched on.
      IF (KCHAIN.GT.0) THEN
          IPHASE = 8
      END IF
*
*       Save KS indices and delay initialization until end of block step.
      CALL DELAY(KCHAIN,KS2)
*
*       Prepare procedure for chain between hierarchy and single body (9/99).
      IF (NAME(I).LT.0.AND.NAME(I).GE.-NZERO.AND.JCOMP.LE.N) THEN
*       Indentify merged ghost particle JG.
          CALL FINDJ(I1,JG,IM)
          if(rank.eq.0)
     &    WRITE (6,28)  NAME(I), NAME(JCOMP), NAME(JG),ECC1, PMIN, RIJ
   28     FORMAT (' HI CHAIN    NAM E1 PM RIJ ',I7,2I6,F7.3,1P,2E10.2)
          JJ = JCOMP
*       Terminate the merger in the usual way.
          KSPAIR = IPAIR
          IPHASE = 7
          CALL RESET
          ZMU = BODY(2*NPAIRS-1)*BODY(2*NPAIRS)/BODY(NTOT)
          EBCH0 = EBCH0 + ZMU*H(NPAIRS)
*       Specify chain indicator and define the two single particles.
          IPHASE = 8
          JCMAX = JG
          JCLOSE = JJ
          KSPAIR = NPAIRS
*       Set relevant variables in DELAY before terminating inner binary.
          CALL DELAY(KCHAIN,KS2)
          CALL DELAY(IPHASE,-1)
*       Initialize new chain of the 4 members JMAX, JCLOSE & KS components.
          ISUB = 0
          CALL CHAIN(ISUB)
*       Note that IPHASE = -1 now and INTGRT goes back to the beginning.
      ELSE IF (NAME(I).LT.-NZERO.OR.NAME(JCOMP).LT.0.OR.
     &        (NAME(I).LT.0.AND.JCOMP.GT.N)) THEN
*       Continue until KS termination on MERGE2 or merger with JCOMP > N.
          IPHASE = 0
      END IF
*
      GO TO 100
*
*       Begin check for merger of stable hierarchical configuration.
   30 RA = SEMI1*(1.0 + ECC1)
      IF (SEMI1.LT.0.0) RA = RIJ
*
*       Identify formation of wide quadruple before merger is accepted.
      IF (JCOMP.GT.N.AND.ECC1.LT.1.0.AND.SEMI1.LT.0.1*RSCALE) THEN
          NNB = LISTQ(1) - 1
          K = 0
*       See whether current list contains first inner/outer component.
          NAM1 = NAME(2*JPAIR-1)
          DO 32 L = 2,NNB+2
              IF (NAM1.EQ.LISTQ(L)) K = K + 1
   32     CONTINUE
*       Generate diagnostics of first five outer orbits every half period.
          IF (K.LE.5.AND.TIME.GT.QCHECK.AND.KZ(37).GT.0) THEN
              ZMB = BODY(I) + BODY(JCOMP)
              RI = SQRT(RI2)
              TK = SEMI1*SQRT(SEMI1/ZMB)
              QCHECK = TIME + MIN(0.5*TWOPI*TK,0.1*TCR)
              TK = DAYS*TK
*       Check the stability criterion.
              PCR = stability(BODY(I1),BODY(I2),BODY(JCOMP),ECC,ECC1,
     &                                                  0.0D0)*SEMIX
              if(rank.eq.0)
     &        WRITE (89,33)  TTOT, NAME(2*IPAIR-1), NAM1, K, RI,
     &                       ECC1, EB, EB2, EB1, TK, PMIN, PCR
   33         FORMAT (' QUAD    T NAM LQ RI E1 EB EB2 EB1 P1 PM PC ',
     &                          F8.1,2I6,I4,F6.2,F8.4,1P,3E12.3,3E9.1)
              CALL FLUSH(89)
*       Remove two oldest members if list is too big.
              IF (NNB.GT.96) THEN
                  DO 34 K = 2,NNB
                      LISTQ(K) = LISTQ(K+2)
   34             CONTINUE
                  NNB = NNB - 2
              END IF
*       Add current names (inner & outer) at end and update membership.
              LISTQ(NNB+3) = NAME(2*IPAIR-1)
              LISTQ(NNB+4) = NAME(2*JPAIR-1)
              LISTQ(1) = NNB + 3
          END IF
      END IF
*
*       Allow temporary merger of inner part of extremely eccentric orbit.
      RFAC = 10.0*RMIN
      IF (ECC1.GT.0.99.AND.RA.GT.RFAC) THEN
          IF (RIJ.LT.0.1*SEMI1) RFAC = RA
      END IF
*
*       Increase apocentre tolerance to local scale factor for EB1 < EBS.
      EBS = 0.25*EBH/SQRT(1.0 + SQRT(RI2)/RSCALE)
      IF (EB1.LT.EBS) THEN
          H2 = (RC**2 + RI2)/FLOAT(NC+10)**0.66667
          RH = 6.0*SQRT(H2/CMSEP2)
          RFAC = MAX(RFAC,RH)
*       Extend maximum apocentre for massive systems (less perturbers).
          IF (BODY(I) + BODY(JCOMP).GT.10.0*BODYM) RFAC = 2.0*RFAC
      END IF
 
*       Skip merger for hyperbolic & soft binding energy or large apocentre.
*     ZF = 1.0
*     IF (BODY(I)*SMU.LT.0.4) ZF = 1.5
      IF (MIN(BODY(I1),BODY(I2)).LT.0.05*BODYM) THEN
*       Use orbital velocity condition instead of binding energy for planets.
          IF (SEMI1.GT.2.0*RMIN) GO TO 100
      ELSE IF (EB.GT.EBH.OR.EB1.GT.EBS.OR.RA.GT.RFAC) THEN
          GO TO 100
      END IF
*
*       Estimate the relative apocentre perturbations on body #I & JCOMP.
      IF (ECC1.LT.0.95) THEN
          PERT = PERT*(RA/RIJ)**3
      ELSE
          PERT = PERT*(ABS(SEMI1)/RIJ)**3
      END IF
      PERTA = PERT4*(RA/RIJ)**3
*
*       Check tidal capture option (synchronous or evolving binary orbit).
      IF (KZ(27).GT.0) THEN
*       Skip merger if outer component would suffer tidal dissipation.
***       IF (SEMI1*(1.0 - ECC1).LT.4.0*RADIUS(JCOMP)) GO TO 100
*       Do not allow merger if Roche overflow or mass loss during next orbit.
          TK = TWOPI*SEMI1*SQRT(SEMI1/(BODY(I) + BODY(JCOMP)))
          TM = MIN(TEV(I1),TEV(I2),TEV(JCOMP),TEV(I))
          IF (KZ(34).GT.0.AND.TM - TIME.LT.STEPX) THEN
              RT = 5.0*MAX(RADIUS(I1),RADIUS(I2))
              IF (A0.LT.RT.OR.KSTAR(I).GT.0) THEN
                  CALL TRFLOW(IPAIR,DTR)
                  IF ((MOD(KSTAR(I),2).EQ.1.AND.DTR.LT.STEPX).OR.
     &                 DTR.LT.TK) GO TO 100
              END IF
              IF (JCOMP.GT.N.AND.KSTAR(JCOMP).GT.0) THEN
                  CALL TRFLOW(JPAIR,DTR)
                  IF(MOD(KSTAR(JCOMP),2).EQ.1.OR.DTR.LT.STEPX) GOTO 100
              END IF
          END IF
*
*       Ensure SLEEP for circularizing binary with TCIRC > 1.0E+06.
          QPERI = A0*(1.0 - ECC)
          RM = MAX(RADIUS(I1),RADIUS(I2))
          IF (KSTAR(I).EQ.-2.AND.QPERI.GT.10.0*RM) THEN
              ICIRC = -1
              CALL INDUCE(IPAIR,EMAX,EMIN,ICIRC,TC,ANGLE,TG,EDAV)
              IF (TC.GT.1.0D+06) THEN
                  if(rank.eq.0)
     &            WRITE (6,35)  NAME(I1), KSTAR(I1), KSTAR(I2), ECC,
     &                          EMAX, QPERI/RM, EDAV, a0, PMIN, TC
   35             FORMAT (' IMPACT SLEEP    NM K* E EX QP/R ED A PM TC',
     &                                      I7,2I4,2F8.4,F6.1,1P,4E10.2)
                  KSTAR(I) = 0
                  NSLP = NSLP + 1
                  II = -I
                  CALL SPIRAL(II)
              END IF
          END IF
*
*       Delay merger for recently updated standard binary and short TCIRC.
          DT = MIN(TEV(I1),TEV(I2)) - TIME
          IF (KSTAR(I).EQ.0.AND.NAME(I).GT.0.AND.DT.LT.TK) THEN
              ICIRC = -1
              CALL TCIRC(QPERI,ECC,I1,I2,ICIRC,TC)
*             if(rank.eq.0)
*    &        WRITE (6,36) NAME(I1), ECC, TTOT, RADIUS(I1)*SU, QPERI, TC
*  36         FORMAT (' TCIRC    NAM E T R* QP TC ',
*    &                           I6,F7.3,F8.3,F7.1,1P,2E10.2)
*       Beware possible termination by routine HMDOT using QPERI < 3*RADIUS.
              IF (TC.LT.2000.0.AND.ECC.GT.0.002) GO TO 100
          END IF
          IF (KZ(19).GE.3) THEN
              IF (MIN(TEV(I1),TEV(I2)).LT.TIME + TK) GO TO 100
          END IF
*       Skip chaotic binary (KSTAR = -2 is treated below).
          IF (KSTAR(I).EQ.-1.OR.KSTAR(JCOMP).EQ.-1) GO TO 100
      END IF
*
*       Ensure consistency of estimated perturbations with termination.
      PERT = PERT + PERTA
      IF (NP.LE.3) THEN
          NP = LIST(1,I)
*       Copy neighbour list (routine FPERT skips body #JCOMP).
          DO 38 L = 1,NP
              JLIST(L) = LIST(L+1,I)
   38     CONTINUE
      END IF
*
*       Allow highly eccentric outer orbit (estimated PERT may be excessive).
      IF (ECC1.GT.0.98.AND.RIJ.LT.0.1*SEMI1) THEN
          PERT = PERT4
          GO TO 40
      END IF
*
*       Evaluate the actual perturbation.
      CALL FPERT(I,JCOMP,NP,PERT2)
      PERT2 = PERT2*RJMIN2/(BODY(I) + BODY(JCOMP))
      GI = PERT2*(RA/RIJ)**3
*       Switch to direct integration for planetary systems if GI > 1D-04.
      IF (MIN(BODY(I1),BODY(I2)).LT.0.05*BODYM) THEN
          IF (GI.GT.1.0D-04) GO TO 100
      END IF
      IF (PERT4.GT.GMAX.OR.GI.GT.0.05) GO TO 100
 
*       Skip merger if an outer binary is fairly perturbed or not hard.
      IF (JCOMP.GT.N) THEN
          IF (GAMMA(JPAIR).GT.1.0E-03.OR.EB2.GT.EBH) GO TO 100
      END IF
*
*       Ensure the inner semi-major axis is used for subsequent tests.
   40 SEMI = -0.5*BODY(I)/H(IPAIR)
      IF (NMERGE.EQ.MMAX - 1) THEN
          IF (NWARN.LT.1000) THEN
              NWARN = NWARN + 1
              if(rank.eq.0)
     &        WRITE (6,41)  NMERGE
   41         FORMAT (5X,'WARNING!    MERGER LIMIT    NMERGE =',I4)
          END IF
          GO TO 100
      END IF
*
*       Do not allow merger in the inner region of perturbed eccentric orbit.
      IF (RIJ.LT.SEMI1.AND.LIST(1,I1).GT.0) THEN
*       Note: moved down from label 30 with 0.1*SEMI1 replacing 2*PMIN.
          IF (ECC1.GT.0.95.AND.RIJ.LT.0.1*SEMI1) THEN
              GO TO 100
          END IF
      END IF
*
*     -----------------------------------------------------------------------
*       Form coefficients for stability test (Valtonen, Vistas Ast 32, 1988).
*     AM = (2.65 + ECC)*(1.0 + BODY(JCOMP)/BODY(I))**0.3333
*     FM = (2.0*BODY(JCOMP) - BODY(I))/(3.0*BODY(I))
*       Note that routine NEWTEV in MERGE2 replaces suppressed part.
*     IF (KZ(19).GE.3) THEN
*         TM = MIN(TEV(I1),TEV(I2),TEV(JCOMP),TEV(I))
*         IF (MIN(NAME(I),NAME(JCOMP)).LT.0.AND.TM-TIME.LT.0.2) THEN
*             GO TO 100
*         END IF
*     END IF
*
*       Expand natural logarithm for small arguments.
*     IF (ABS(FM).LT.0.67) THEN
*         BM = FM*(1.0 - (0.5 - ONE3*FM)*FM)
*     ELSE
*         BM = LOG(1.0D0 + FM)
*     END IF
*
*       Adopt mass dependent criterion of Harrington (A.J. 82, 753) & Bailyn.
*     PCRIT = AM*(1.0 + 0.7*BM)*SEMI
*     -----------------------------------------------------------------------
*
*       Form hierarchical stability ratio (Eggleton & Kiseleva 1995).
*     QL = BODY(I)/BODY(JCOMP)
*     Q1 = MAX(BODY(I2)/BODY(I1),BODY(I1)/BODY(I2))
*     Q3 = QL**0.33333
*     Q13 = Q1**0.33333
*     AR = 1.0 + 3.7/Q3 - 2.2/(1.0 + Q3) + 1.4/Q13*(Q3 - 1.0)/(Q3 + 1.0)
*     EK = AR*SEMI*(1.0D0 + ECC)
*
*       Choose the most dominant triple in case of two binaries.
      IF (JCOMP.GT.N) THEN
*       Adopt 10% fudge factor with linear dependence on smallest ratio.
          YFAC = 1.0 + 0.1*MIN(SEMI2/SEMI,SEMI/SEMI2)
      ELSE
          YFAC = 1.0
      END IF
*
*       Prepare inclination evaluation for triple or widest inner binary.
      IF (JCOMP.GT.N) THEN
*       Ensure widest inner binary (swap is OK for termination or ZARE).
          IF (SEMI.LT.SEMI2) THEN
              ECC2 = (1.0 - R(JPAIR)/SEMI2)**2 +
     &                             TDOT2(JPAIR)**2/(BODY(JCOMP)*SEMI2)
              ECC = SQRT(ECC2)
              KPAIR = IPAIR
              IPAIR = JPAIR
              JPAIR = KPAIR
              I1 = 2*IPAIR - 1
              I2 = I1 + 1
              JJ = I
              I = JCOMP
              JCOMP = JJ
              SEMIZ = SEMI2
              SEMI2 = SEMI
              SEMI = SEMIZ
          END IF
      END IF
*
*       Resolve binary (even perturbed KS not always done).
      IF (ECC1.LT.1.0) THEN
          CALL RESOLV(IPAIR,1)
      END IF
*
*       Copy coordinates and velocities to local variables.
      DO 42 K = 1,3
          XX(K,1) = X(K,I1)
          XX(K,2) = X(K,I2)
          XX(K,3) = X(K,JCOMP)
          VV(K,1) = XDOT(K,I1)
          VV(K,2) = XDOT(K,I2)
          VV(K,3) = XDOT(K,JCOMP)
  42  CONTINUE
*
*       Determine the inclination (in radians).
      CALL INCLIN(XX,VV,X(1,I),XDOT(1,I),ANGLE)
*
*       Employ the basic stability criterion for fast check (ECC1 < 1).
*     IF (ECC1.LT.1.0) THEN
*         Q = BODY(JCOMP)/BODY(I)
*         XFAC = (1.0 + Q)*(1.0 + ECC1)/SQRT(1.0 - ECC1)
*         PCRIT = 2.8*XFAC**0.4*SEMI*(1.0 - 0.6*ANGLE/TWOPI)
*         IF (PCRIT.GT.2.0*PMIN) GO TO 100
*     END IF
*
*       Evaluate the general stability function (Mardling MNRAS, 2008).
      IF (ECC1.LT.1.0.AND.YFAC.LT.1.02) THEN
          BJ = BODY(JCOMP)
          EOUT = ECC1
*       Increase tolerance near sensitive stability boundary (RM 10/2008).
          IF (EOUT.GT.0.8) THEN
              DE = 0.5*(1.0 - EOUT)
              DE = MIN(DE,0.01D0)
*       Allow extra tolerance after 1000 tries.
              IF (NMTRY.GE.1000) DE = MIN(1.0D0 - EOUT,0.02D0)
              EOUT = EOUT - DE
              PMIN = SEMI1*(1.0 - EOUT)
          END IF
          NST = NSTAB(SEMI,SEMI1,ECC,EOUT,ANGLE,BODY(I1),BODY(I2),BJ)
          IF (NST.EQ.0) THEN
              PCRIT = 0.98*PMIN*(1.0 - PERT)
              PCR = stability(BODY(I1),BODY(I2),BODY(JCOMP),ECC,EOUT,
     &                                                          ANGLE)
              PCR = PCR*SEMI
*       Specify reduced peri if old criterion < PCRIT/2 (avoids switching).
              IF (PCR.LT.0.5*PCRIT) THEN
                  PCRIT = 0.75*PCRIT
              END IF
              IF (PCRIT.LT.YFAC*PCR.AND.PERT.LT.0.02.AND.
     &            NMTRY.LT.10) THEN
                  ALPH = 360.0*ANGLE/TWOPI
                  FAIL = PMIN*(1-PERT) - YFAC*PCR
                  if(rank.eq.0)
     &            WRITE (6,43)  TTOT, ALPH, ECC, ECC1, PMIN, FAIL, PERT
   43             FORMAT (' NEWSTAB    T INC EI EO PMIN FAIL PERT ',
     &                                 F7.1,F7.2,2F8.4,1P,3E10.2)
              END IF
*     if(rank.eq.0)then
*     WRITE (57,444)  BODY(I1),(X(K,I1),K=1,3),(XDOT(K,I1),K=1,3)
*     WRITE (57,444)  BODY(I2),(X(K,I2),K=1,3),(XDOT(K,I2),K=1,3)
*     WRITE (57,444)BODY(JCOMP),(X(K,JCOMP),K=1,3),(XDOT(K,JCOMP),K=1,3)
* 444 FORMAT (' ',1P,E14.6,1P,3E18.10,3E14.6)
*     end if
              IF (NMTRY.GE.1000) THEN
*                 if(rank.eq.0)
*    &            WRITE (6,44)  TTOT, NAME(I1), ECC1, EOUT, PCRIT, PMIN
*  44             FORMAT (' MARGINAL STAB    T NM E1 EOUT PCR PMIN ',
*    &                                       F7.1,I7,2F8.4,1P,2E10.2)
                  NMTRY = 0
              END IF
          ELSE
              NMTRY = NMTRY + 1
              PCRIT = 1.01*PMIN
          END IF
      ELSE
          PCR = stability(BODY(I1),BODY(I2),BODY(JCOMP),ECC,ECC1,ANGLE)
          PCRIT = PCR*SEMI
      END IF
*
*       Check whether the main perturber dominates the outer component.
      IF (JMAX.NE.JCOMP) THEN
          RIJ2 = (X(1,JMAX) - X(1,JCOMP))**2 +
     &           (X(2,JMAX) - X(2,JCOMP))**2 +
     &           (X(3,JMAX) - X(3,JCOMP))**2
          FMAX = (BODY(JMAX) + BODY(JCOMP))/RIJ2
          IF (FMAX.GT.(BODY(I) + BODY(JCOMP))/RJMIN2) GO TO 100
      END IF
*
*       Determine time-scale for stability (absolute or approximate).
      PM1 = PMIN*(1.0 - 2.0*PERT)
      CALL TSTAB(I,ECC1,SEMI1,PM1,YFAC,ITERM)
      IF (ITERM.GT.0) GO TO 100
*
*       Check perturbed stability condition.
      IF (PMIN*(1.0 - PERT).LT.YFAC*PCRIT) GO TO 100
*
*       Extend active Roche case up to end of look-up time (bug fix 01/09).
      IF (KSTAR(I).GE.11.AND.MOD(KSTAR(I),2).NE.0) THEN
          DT = TEV(I) - TIME
          IF (DT.LT.10.0*STEPX) GO TO 100
          TMDIS(NMERGE+1) = MIN(TIME + DT, TMDIS(NMERGE+1))
      END IF
*
*       Obtain growth time and inclination for significant eccentricity.
      IF (SEMI1.GT.0.0.AND.ECC.GT.0.1) THEN
          ICIRC = -1
          CALL INDUCE(IPAIR,EMAX,EMIN,ICIRC,TC,ANGLE,TG,EDAV)
*       Prevent termination for TC < 2000 in HMDOT but allow small EMAX & DE.
          IF (KZ(27).EQ.2.AND.TC.LT.2000.0) THEN
              DE = ABS(EMAX - ECC)
              DT = MIN(TEV(I1),TEV(I2),TEV(I)) - TIME
*       Enforce an update next block-step for small remaining times.
              IF (DT.LT.0.1.AND.MOD(KSTAR(I),2).EQ.0) THEN
                  TEV(I1) = TIME + 0.1
                  TEV(I2) = TIME + 0.1
                  if(rank.eq.0)
     &            WRITE (6,46)  TIME+TOFF, NAME(I1), KSTAR(I), ECC, TC
   46             FORMAT (' ENFORCED TEV UPDATE    T NM K* E TC ',
     &                                             F9.2,I8,I4,F8.4,F7.1)
              END IF
*       Accept KSTAR = 10, coasting Roche or small eccentricity.
              IF (MOD(KSTAR(I),2).EQ.0) THEN
                  CALL TRFLOW(IPAIR,DTR)
                  IF (DTR.LT.STEP(I)) THEN
                      TEV(I) = TIME + DTR
                  END IF
              ELSE
                  IF (EMAX.GT.0.8.OR.DE.GT.0.2.OR.DT.LT.0.1) GO TO 100
              END IF
          END IF
      END IF
*
*       Check circularization and dissipation time (exclude Roche stages).
      IF (KZ(27).EQ.2.AND.KSTAR(I).LT.10.AND.ECC1.LT.1.0.AND.
     &    NAME(I).GT.0) THEN
          ECC0 = ECC
          CALL DECIDE(IPAIR,SEMI,ECC,EMAX,EMIN,TC,TG,EDAV,IQ)
          IF (IQ.GT.0.OR.IPHASE.LT.0) GO TO 100
          TK1 = TWOPI*SEMI1*SQRT(SEMI1/(BODY(I) + BODY(JCOMP)))
*         IF (TMDIS(NMERGE+1) - TIME.LT.TK1) GO TO 100
*         EK = EK*(1.0 + ECC)/(1.0 + ECC0)
          PCRIT = PCRIT*(1.0 + ECC)/(1.0 + ECC0)
      END IF
*
*       Perform safety check on radii for case of no circularization.
      IF (KZ(27).EQ.0) THEN
          IF (SEMI*(1.0 - ECC).LT.2.0*MAX(RADIUS(I1),RADIUS(I2))) THEN
              IF (KZ(19).EQ.0) THEN
                  GO TO 100
              END IF
          END IF
      END IF
*
*       Include rare case of circularizing outer binary.
      IF (KZ(27).EQ.2.AND.JCOMP.GT.N.AND.KSTAR(JCOMP).EQ.-2) THEN
          ECC2 = (1.0 - R(JPAIR)/SEMI2)**2 +
     &                               TDOT2(JPAIR)**2/(BODY(JCOMP)*SEMI2)
          ECCJ = SQRT(ECC2)
*       See whether to reduce the look-up time TMDIS (no skip here).
          CALL DECIDE(JPAIR,SEMI2,ECCJ,EMAXJ,EMIN,TC,TG,EDAV,IQ)
          IF (IPHASE.LT.0) GO TO 100
      END IF
*
*       Check Zare exchange stability criterion and create diagnostics.
      IF (SEMI1.GT.0.0) THEN
          CALL ZARE(I1,I2,SP)
          Q = BODY(JCOMP)/BODY(I)
*       Note inclination is determined by routine INCLIN for ECC < 0.1.
          IF (SP.LT.1.0.AND.ANGLE.LT.0.17) THEN
              IZARE = IZARE + 1
              IF (IZARE.LT.200) THEN
              if(rank.eq.0)then
              WRITE (7,48)  TTOT, Q, ECC, ECC1, SEMI, PMIN, PCRIT,
     &                      YFAC, SP
              WRITE (7,47) I,JCOMP,N,I1,I2,RIJ,SEMI1
   47         FORMAT (' I JCOMP N I1 I2 RIJ A1   ',5I6,1P,2E10.2)
              CALL FLUSH(7)
              WRITE (6,48)  TTOT, Q, ECC, ECC1, SEMI, PMIN, PCRIT,
     &                      YFAC, SP
   48         FORMAT (' ZARE TEST    T Q E E1 A PM PCR YF SP ',
     &                               F8.2,F5.1,2F7.3,1P,3E9.1,0P,2F6.2)
              end if
              END IF
          END IF
          if(rank.eq.0)
     &    WRITE (73,49)  TTOT, Q, ECC, ECC1, SEMI, PMIN, PCRIT,
     &                   TG, SP, 360.0*ANGLE/TWOPI, KSTAR(I)
   49     FORMAT (' STAB    T Q E E1 A PM PCR TG SP IN K* ',
     &                      F8.2,F5.1,2F7.3,1P,4E9.1,0P,F6.2,F7.1,I4)
          CALL FLUSH(73)
*         IF (KSTAR(I1).GE.10) TEV(I1) = 1.0E+10
*         IF (KSTAR(I2).GE.10) TEV(I2) = 1.0E+10
      END IF
*
*       Specify the final critical pericentre using the fudge factor.
      PCRIT = YFAC*PCRIT
*
      IF (NAME(I).GT.N.AND.NAME(JCOMP).GT.N.AND.ECC1.GT.1.0) GO TO 100
*       Skip if #JCOMP is a chain c.m. but allow bound double hierarchy.
      IF (NAME(JCOMP).EQ.0) GO TO 100
      IF (ECC1.GT.1.0.AND.MIN(NAME(I),NAME(JCOMP)).LT.0) GO TO 100
      DO 55 ISUB = 1,NSUB
          IF (NAME(JCOMP).EQ.NAMES(1,ISUB)) GO TO 100
   55 CONTINUE
*       Do not allow the formation of a SEPTUPLET.
*     IF ((NAME(I).LT.-2*NZERO.AND.JCOMP.GT.N).OR.
*    &     NAME(JCOMP).LT.-2*NZERO) GO TO 100
*
*       Include diagnostics for double hierarchy or optional standard case.
      IF (NAME(I).LT.0.OR.NAME(JCOMP).LT.0) THEN
          IF (KZ(15).GT.1) THEN
              WHICH1 = ' MERGE2 '
          RI = SQRT((X(1,I) - RDENS(1))**2 +
     &              (X(2,I) - RDENS(2))**2 +
     &              (X(3,I) - RDENS(3))**2)
          VI = SQRT(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2)
          if(rank.eq.0)
     &    WRITE (6,20)  WHICH1, TTOT, I, JCOMP, IPAIR, NAME(I1),
     &        NAME(I2), NAME(JCOMP), NAME(I), H(IPAIR), R(IPAIR),
     &                  BODY(I), BODY(JCOMP), PERT4, RIJ, PMIN,
     &                  EB1/EB, RI, VI, LIST(1,I1)
          END IF
*       Note rare case of two hierarchies merging and identify ghost names.
          IF (NAME(I).LT.0.AND.NAME(JCOMP).LT.0) THEN
              CALL FINDJ(I1,JI,IM)
              J1 = 2*JPAIR - 1
              CALL FINDJ(J1,JJ,JM)
              if(rank.eq.0)
     &        WRITE (6,60)  NAME(I1), NAME(JI), NAME(I1+1), NAME(J1),
     &                      NAME(JJ), NAME(J1+1), ECC, ECC1, SEMI,
     &                      SEMI1, PMIN, PCRIT
   60         FORMAT (' HI MERGE    NAM E E1 A A1 PM PC ',
     &                              6I7,2F7.3,1P,4E10.2)
          END IF
      ELSE IF (KZ(15).GT.1) THEN
          WHICH1 = ' MERGER '
          IF (JCOMP.GT.N) WHICH1 = ' QUAD   '
          RI = SQRT((X(1,I) - RDENS(1))**2 +
     &              (X(2,I) - RDENS(2))**2 +
     &              (X(3,I) - RDENS(3))**2)
          VI = SQRT(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2)
          if(rank.eq.0)
     &    WRITE (6,20)  WHICH1, TTOT, I, JCOMP, IPAIR, NAME(I1),
     &           NAME(I2), NAME(JCOMP), NAME(I), H(IPAIR), R(IPAIR),
     &                  BODY(I), BODY(JCOMP), PERT4, RIJ, PMIN,
     &                  EB1/EB, RI, VI, LIST(1,I1)
      END IF
*
*       Check for diagnostic output of quadruples.
      IF (SEMI1.GT.0.0.AND.JCOMP.GT.N.AND.KZ(37).GT.0) THEN
          ZMB = BODY(I) + BODY(JCOMP)
          TK = DAYS*SEMI1*SQRT(SEMI1/ZMB)
          if(rank.eq.0)
     &    WRITE (89,65)  TTOT, NAME(2*IPAIR-1), NAME(2*JPAIR-1),
     &                   LISTQ(1), SQRT(RI2), ECC1, EB, EB2, EB1,
     &                   TK, PMIN, PCRIT
   65     FORMAT (' QUAD#   T NAM LQ RI E1 EB EB2 EB1 P1 PM PC ',
     &                      F8.1,2I6,I4,F6.2,F8.4,1P,3E12.3,3E9.1)
      END IF
*
*       Generate a diagnostic file of stable hierarchies (suppressed).
      IF (ECC1.LT.-1.0) THEN
          RI = SQRT(RI2)/RC
          if(rank.eq.0)
     &    WRITE (80,70)  TPHYS, RI, NAME(JCOMP), QL, Q1, ECC, ECC1,
     &                   SEMI, SEMI1, PCRIT/PMIN, 360.0*ANGLE/TWOPI,EMAX
   70     FORMAT (F8.1,F5.1,I6,2F6.2,2F6.3,1P,2E10.2,0P,F5.2,F6.1,F6.3)
          CALL FLUSH(80)
      END IF
*
*       Copy pair index and set indicator for calling MERGE from MAIN.
      KSPAIR = IPAIR
      IPHASE = 6
*
*       Save KS indices and delay merger until end of block step.
      CALL DELAY(KS2,KS2)
*
  100 IF (IPHASE.NE.8) JCMAX = 0
*
      RETURN
*
      END
