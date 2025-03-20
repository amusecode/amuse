      SUBROUTINE DERQP4(P,Q,DP,DQ)
*
*
*       Equations of motion for chain regularization.
*       ---------------------------------------------
*
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8  M,MIJ,XNR(9),FNR(9),WI(9),AT(3),AQ(3),AP(9),FQ(9),WP(9),
     &        P(12),Q(12),DP(12),DQ(13),UPR(12),MB
      LOGICAL  SWITCH,GTYPE,GTYPE0
      COMMON/CREG/  M(4),X(12),XD(12),PP(12),QQ(12),TIME4,ENERGY,EPSR2,
     &              XR(9),W(9),R(6),TA(6),MIJ(6),CM(10),RMAX4,TMAX,
     &              DS,TSTEP,EPS,NSTEP4,NAME4(4),KZ15,KZ27,NREG,NFN
      COMMON/TPR/   SWITCH,GTYPE,GTYPE0
      COMMON/ICONF/  I1,I2,I3,I4
      COMMON/CLOSE/  RIJ(4,4),RCOLL,QPERI,SIZE(4),ECOLL3,IP(4)
      COMMON/CCOLL/  QK(12),PK(12),ICALL,ICOLL,NDISS4
      COMMON/BSSAVE/  EP(4),DSC,FACM,TFAC,ITFAC,JC
      COMMON/KSAVE/  K10,K20
      EQUIVALENCE  (T11,TA(1)),(T22,TA(2)),(T33,TA(3)),(T12,TA(4)),
     &             (T23,TA(5))
*     EQUIVALENCE  (M12,MIJ(1)),(M23,MIJ(2)),(M34,MIJ(3)),
*    &             (M13,MIJ(4)),(M24,MIJ(5)),(M14,MIJ(6))
*
*
      NFN = NFN + 1
      K = 0
      KP = 0
      DO 1 L = 1,3
          K1 = K + 1
          K2 = K + 2
          K3 = K + 3
          K4 = K + 4
          KP1 = KP + 1
          KP2 = KP + 2
          KP3 = KP + 3
*
*       Form physical vector & scalar distance.
          XR(KP1) = Q(K1)**2 - Q(K2)**2 - Q(K3)**2 + Q(K4)**2
          XR(KP2) = 2.D0*(Q(K1)*Q(K2) - Q(K3)*Q(K4))
          XR(KP3) = 2.D0*(Q(K1)*Q(K3) + Q(K2)*Q(K4))
          R(L) = Q(K1)**2 + Q(K2)**2 + Q(K3)**2 + Q(K4)**2
*       Set physical momentum times distance.
          WI(KP1) = Q(K1)*P(K1) - Q(K2)*P(K2) - Q(K3)*P(K3) +
     &                                          Q(K4)*P(K4)
          WI(KP2) = Q(K2)*P(K1) + Q(K1)*P(K2) - Q(K4)*P(K3) -
     &                                          Q(K3)*P(K4)
          WI(KP3) = Q(K3)*P(K1) + Q(K4)*P(K2) + Q(K1)*P(K3) +
     &                                          Q(K2)*P(K4)
          AT(L) = TA(L)*(P(K1)**2 + P(K2)**2 + P(K3)**2 +
     &                                         P(K4)**2) - MIJ(L)
          K = K4
          KP = KP3
    1 CONTINUE
*
*       Obtain irregular vectors & distances.
      R(4) = 0.0D0
      R(5) = 0.0D0
      R(6) = 0.0D0
      DO 2 K = 1,3
          XNR(K  ) = XR(K  ) + XR(K+3)
          XNR(K+3) = XR(K+3) + XR(K+6)
          XNR(K+6) = XNR(K ) + XR(K+6)
          R(4) = R(4) + XNR(K  )**2
          R(5) = R(5) + XNR(K+3)**2
          R(6) = R(6) + XNR(K+6)**2
    2 CONTINUE
*
*       Evaluate irregular physical forces.
      UNR = 0.0D0
      DO 5 I = 4,6
          DIST = SQRT(R(I))
          FC = MIJ(I)/DIST
          UNR = UNR + FC
          FC = FC/R(I)
          R(I) = DIST
          IK = 3*(I - 4)
          DO 3 K = 1,3
              FNR(IK+K) = FC*XNR(IK+K)
    3     CONTINUE
    5 CONTINUE
      W1W2 = WI(1)*WI(4) + WI(2)*WI(5) + WI(3)*WI(6)
      W2W3 = WI(4)*WI(7) + WI(5)*WI(8) + WI(6)*WI(9)
      UNRE = UNR + ENERGY
      R12 = R(1)*R(2)
      R13 = R(1)*R(3)
      R23 = R(2)*R(3)
      AQ(1) = AT(2)*R(3) + AT(3)*R(2) + T23*W2W3-R23*UNRE
      AQ(2) = AT(1)*R(3) + AT(3)*R(1) - R13*UNRE
      AQ(3) = AT(2)*R(1) + AT(1)*R(2) + T12*W1W2-R12*UNRE
      R123 = R12*R(3)
      AP(1) = 2.0D0*T11*R23
      AP(2) = 2.0D0*T22*R13
      AP(3) = 2.0D0*T33*R12
      WK12 = T12*R(3)
      WK23 = T23*R(1)
      DO 6 K = 1,3
          WP(K  ) = WK12*WI(K+3)
          WP(K+3) = WK12*WI(K  ) + WK23*WI(K+6)
          WP(K+6) = WK23*WI(K+3)
          FQ(K  ) = FNR(K  ) + FNR(K+6)
          FQ(K+3) = FNR(K  ) + FNR(K+3) + FNR(K+6)
          FQ(K+6) = FNR(K+3) + FNR(K+6)
    6 CONTINUE
*
*       Form regularized derivatives.
      KQ = 0
      L = 0
      DO 10 I = 1,3
          K1 = KQ + 1
          K2 = KQ + 2
          K3 = KQ + 3
          K4 = KQ + 4
          L1 = L + 1
          L2 = L + 2
          L3 = L + 3
*
          F1 = +FQ(L1)*Q(K1) + FQ(L2)*Q(K2) + FQ(L3)*Q(K3)
          F2 = -FQ(L1)*Q(K2) + FQ(L2)*Q(K1) + FQ(L3)*Q(K4)
          F3 = -FQ(L1)*Q(K3) - FQ(L2)*Q(K4) + FQ(L3)*Q(K1)
          F4 = +FQ(L1)*Q(K4) - FQ(L2)*Q(K3) + FQ(L3)*Q(K2)
*
          G1 = +WP(L1)*P(K1) + WP(L2)*P(K2) + WP(L3)*P(K3)
          G2 = -WP(L1)*P(K2) + WP(L2)*P(K1) + WP(L3)*P(K4)
          G3 = -WP(L1)*P(K3) - WP(L2)*P(K4) + WP(L3)*P(K1)
          G4 = +WP(L1)*P(K4) - WP(L2)*P(K3) + WP(L3)*P(K2)
*
          RRR = R123 + R123
          DP(K1) = -(2.D0*AQ(I)*Q(K1) + G1+RRR*F1)
          DP(K2) = -(2.D0*AQ(I)*Q(K2) + G2+RRR*F2)
          DP(K3) = -(2.D0*AQ(I)*Q(K3) + G3+RRR*F3)
          DP(K4) = -(2.D0*AQ(I)*Q(K4) + G4+RRR*F4)
*
          G1 = +WP(L1)*Q(K1) + WP(L2)*Q(K2) + WP(L3)*Q(K3)
          G2 = -WP(L1)*Q(K2) + WP(L2)*Q(K1) + WP(L3)*Q(K4)
          G3 = -WP(L1)*Q(K3) - WP(L2)*Q(K4) + WP(L3)*Q(K1)
          G4 = +WP(L1)*Q(K4) - WP(L2)*Q(K3) + WP(L3)*Q(K2)
*
          DQ(K1) = AP(I)*P(K1) + G1
          DQ(K2) = AP(I)*P(K2) + G2
          DQ(K3) = AP(I)*P(K3) + G3
          DQ(K4) = AP(I)*P(K4) + G4
*
          KQ = KQ + 4
          L = L + 3
   10 CONTINUE
*
*       Set the time derivative and check tolerance scaling (ITFAC > 1).
      DQ(13) = R123
      IF (ITFAC.GT.0) THEN
          TFAC = FACM*(R(1)*R(2) + R(1)*R(3) + R(2)*R(3))
          ITFAC = -1
      END IF
*
*       Procedure for obtaining GAMMA (only after first call).
*     IF (IDIAG.EQ.1) THEN
*         GAMMA  =  AT(1)*R23 + AT(2)*R13 + AT(3)*R12 +
*    &              (R(3)*T12*W1W2 + T23*R(1)*W2W3) - R123*UNRE
*         IDIAG  =  0
*       NB! IDIAG must be in COMMON and set > 0 by RCHAIN.
*     END IF
*
      IF (GTYPE) THEN
*       Use modified time transformation for critical case.
          GAMMA = AT(1)*R23 + AT(2)*R13 + AT(3)*R12 +
     &            R(3)*T12*W1W2 + T23*R(1)*W2W3 - R123*UNRE
          SIGIN = 1.D0/(R12 + R23 + R13)
          GSIGIN = 2.D0*GAMMA*SIGIN
          SUMR = R(1) + R(2) + R(3)
          DO 15 I = 1,3
              SI = (SUMR - R(I))*GSIGIN
              KQ = 4*(I - 1)
              DO 12 K = KQ+1,KQ+4
                  DQ(K) = SIGIN*DQ(K)
                  DP(K) = SIGIN*(DP(K) + SI*Q(K))
   12         CONTINUE
   15     CONTINUE
          DQ(13) = DQ(13)*SIGIN
          IF (ITFAC.NE.0) THEN
              TFAC = TFAC*SIGIN
              ITFAC = 0
          END IF
      END IF
*
*       Check osculating pericentre of closest pair (first call only).
      IF (ICALL.EQ.0) GO TO 50
*
*       Examine all close pair indices K1 & K2 in chain vector /ICONF/.
      QPMIN = 100.0
      IM0 = 1
      DO 20 IM = 1,3
          IF (IM.EQ.1) THEN
              K1 = I1
              K2 = I2
          ELSE IF (IM.EQ.2) THEN
              K1 = I2
              K2 = I3
          ELSE
              K1 = I3
              K2 = I4
          END IF
*
*       Set distance of nearest perturber.
          IF (IM.EQ.1.OR.IM.EQ.3) THEN
              RP = R(2)
          ELSE
              RP = MIN(R(1),R(3))
          END IF
*
*       Obtain pericentre for small perturbations (ignore mass effect).
          GI = (R(IM)/RP)**3
          IF (GI.LT.0.005) THEN
              K = 4*(IM - 1) + 1
              CALL PERI(Q(K),DQ(K),DQ(13),M(K1),M(K2),QPERI)
          ELSE
              QPERI = R(IM)
          END IF
*
*       Compare pericentre with previous mutual distances.
          RIJ(K1,K2) = MIN(RIJ(K1,K2),QPERI)
          RIJ(K2,K1) = MIN(RIJ(K2,K1),QPERI)
*
*       Save indices for smallest pericentre and switch off indicator.
          IF (IM.EQ.1.OR.R(IM).LE.R(IM0)) THEN
              QPMIN = QPERI
              RP0 = RP
              K10 = K1
              K20 = K2
              IM0 = IM
              ICALL = 0
          END IF
   20 CONTINUE
*
*       Check smallest pericentre and reset corresponding indices.
      RCOLL = MIN(RCOLL,QPMIN)
      QPERI = QPMIN
      IM = IM0
      K = 4*(IM - 1) + 1
      K1 = K10
      K2 = K20
*
*       Save minimum configuration (for EREL4 & TRANS4).
      DO 25 J = 1,12
          QK(J) = Q(J)
          PK(J) = P(J)
   25 CONTINUE
*
*       Check for tidal two-body interaction or stellar collision.
      IF (QPMIN.LT.4.0*MAX(SIZE(K1),SIZE(K2))) THEN
*
*       Exit if collision distance is not the smallest.
          IF (R(IM).GT.RP0) GO TO 50
*
*       Convert Q' to standard KS with T' = R and form radial velocity R'.
          RPR = 0.0D0
          DO 30 J = K,K+3
              UPR(J) = DQ(J)*R(IM)/DQ(13)
              RPR = RPR + 2.0D0*Q(J)*UPR(J)
   30     CONTINUE
*
*       Determine small semi-major axis from non-singular expressions.
          CALL EREL4(IM,EB,SEMI)
*
*       Obtain pericentre time interval from elements & Kepler's equation.
          MB = M(K1) + M(K2)
          CALL TPERI(SEMI,Q(K),UPR(K),MB,TP)
*
*       Activate collision indicator & B-S step selector (first time).
          IF (ICOLL.EQ.0) THEN
              ICOLL = -1
              JC = 1
          ELSE
*
*       Check convergence: radial velocity < 1.0E-06 parabolic velocity.
              IF (ABS(RPR).LT.1.0E-06*SQRT(2.0D0*MB*R(IM)).OR.
     &           (ABS(DSC).LT.1.0E-06)) THEN
*       Reset B-S step selector and set collision indicator to CHAIN index.
                  JC = 0
                  ICOLL = IM
              END IF
          END IF
*
*       Set regularized step for DIFSY4 using accelerated iteration rate.
          IF (R(IM).GT.2.0*QPMIN) THEN
              DSC = 2.0*ABS(TP)/DQ(13)
          ELSE
              DSC = ABS(TP)/DQ(13)
          END IF       
*
*       Ensure negative step beyond pericentre (restore step at end).
          IF (JC.GT.0.AND.RPR.GT.0.0D0) THEN
              DSC = -DSC
          END IF
          IF (JC.EQ.0) DSC = 1.0
*
*       Begin iteration if R(IM) is smallest or wide pair is past pericentre.
          IF (R(IM).LT.RP0.OR.RPR.GT.0.0D0) THEN
              DSFAC = 0.2*ABS(DS/DSC)
              IF (R(IM).LT.RP0) DSFAC = 1.0
              DS = DSFAC*DSC
          ELSE
*       Switch off indicators and carry on normally until next check.
              ICOLL = 0
              JC = 0
          END IF
      END IF
*
   50 RETURN
*
      END
