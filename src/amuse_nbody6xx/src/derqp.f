      SUBROUTINE DERQP(Q,CMX,ENERGY,P,CMV,CHTIME,DQ,DX,DE,DP,DV,DT)
*
*
*       Derivatives of chain variables.
*       -------------------------------
*
      INCLUDE 'commonc.h'
      LOGICAL  KSLOW,KCOLL
      REAL*8  KSCH
      COMMON/SLOW1/   TK2(0:NMX),EJUMP,KSCH(NMX),KSLOW,KCOLL
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
      COMMON/CALLS/  TPR,TKPR,STEP,IDER,ICALL,NFN,NREG,ITER,IMCIRC
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIJ(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/BSSAVE/  EP(4),DSC,FACM,TFAC,ITFAC,JC
      COMMON/EBSAVE/  EBS
      COMMON/KSAVE/  K10,K20
      REAL*8  Q(NMX4),CMX(3),P(NMX4),CMV(3)
      REAL*8  DQ(NMX4),DX(3),DP(NMX4),DV(3),FCM(3),UPR(NMX4)
      REAL*8  W(NMX4),AK(NMX4),DK(NMX),FNC(NMX3),FXTNL(NMX3)
      REAL*8  TP(NMX4),TQ(NMX4),UQ(NMX4),FAUX(4),AUX(-2:2),XAUX(3)
      INCLUDE "mpif.h"
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      COMMON/MPIDAT/group,rank,ierr,isize,status
*
*
*       Mikkola & Aarseth 1990 eqs. (77) -> (80).
      UC=0.0
      RSUM=0.0
*       Eq. (77).
      DO I=1,N-1
      KS = 4*(I-1)
      K1 = KS + 1
      K2 = KS + 2
      K3 = KS + 3
      K4 = KS + 4
*       Obtain W = L^T(Q)P ; L^T = transpose of L times P.
      W(K1)=(Q(K1)*P(K1)-Q(K2)*P(K2)-Q(K3)*P(K3)+Q(K4)*P(K4))
      W(K2)=(Q(K2)*P(K1)+Q(K1)*P(K2)-Q(K4)*P(K3)-Q(K3)*P(K4))
      W(K3)=(Q(K3)*P(K1)+Q(K4)*P(K2)+Q(K1)*P(K3)+Q(K2)*P(K4))
      W(K4)=(Q(K4)*P(K1)-Q(K3)*P(K2)+Q(K2)*P(K3)-Q(K1)*P(K4))
      RIJL=Q(K1)**2+Q(K2)**2+Q(K3)**2+Q(K4)**2
*       Form RSUM for decision-making.
      RSUM=RSUM+RIJL
      RINV(I)=1./RIJL
      A=.5D0*RINV(I)
      UC=UC+MKK(I)*RINV(I)
      DO K=1,4
      W(KS+K)=A*W(KS+K)
      END DO
      END DO
      LRI=N-1
*       Eq. (78).
      TKIN=0.0
      DO I=1,N-1
      AUX(-2)=0.5D0*TK2(I-1)
      AUX(-1)=0.5D0*TK1(I)
      AUX( 0)=TKK(I)
      AUX(+1)=0.5D0*TK1(I+1)
      AUX(+2)=0.5D0*TK2(I+1)
      L=4*(I-1)
      DK(I)=0.0
      DO K=1,4
      AA=0.0
      DO J=-2,2
      LJ=L+4*J
      if (lj.ge.0.and.lj.le.4*n-8) then
      AA=AA+AUX(J)*W(LJ+K)
      end if
      END DO
      AK(L+K)=AA
*       Eq. (79).
      DK(I)=DK(I)+AA*W(L+K)
      END DO
*       Eq. (80).
      TKIN=TKIN+DK(I)
      END DO
*
*       Obtain physical coordinates.
      DO K=1,3
      XI(K)=0.0
      END DO
      DO I=1,N-1
      L=3*(I-1)
      KS=4*(I-1)
      XC(L+1)=Q(KS+1)**2-Q(KS+2)**2-Q(KS+3)**2+Q(KS+4)**2
      XC(L+2)=2.D0*(Q(KS+1)*Q(KS+2)-Q(KS+3)*Q(KS+4))
      XC(L+3)=2.D0*(Q(KS+1)*Q(KS+3)+Q(KS+2)*Q(KS+4))
      DO K=1,3
      XI(L+3+K)=XI(L+K)+XC(L+K)
      END DO
      END DO
*
*       Evaluate external force (FXTNL = W', not the 'usual' force!).
      IF (IPERT.GT.0) THEN
      ISKIP=0
      CALL XTF(FXTNL,FCM,CMX,CHTIME)
*       Note that perturbation limit must be consistent with chlist.f.
      IF (GPERT.LT.1.0D-06) THEN
      IPERT=0
      ELSE
      IPERT=1
      END IF
      ELSE
      ISKIP=1
      DO I=1,3*(N-1)
      FXTNL(I)=0.0D0
      END DO
      DO I=1,3
      FCM(I)=0.0D0
      END DO
      END IF
*
*	Form non-chained contributions.
      UNC=0.0
      DO I=1,3*(N-1)
      FNC(I)=FXTNL(I)
      END DO
*
      DO I=1,N-2
      LI=3*(I-1)
      DO J=I+2,N
      LJ=3*(J-1)
      RIJ2=0.0
      IF(J.GT.I+2)THEN
      DO K=1,3
      XAUX(K)=XI(LJ+K)-XI(LI+K)
      RIJ2=RIJ2+XAUX(K)**2
      END DO
      ELSE
      DO K=1,3
      XAUX(K)=XC(LI+K)+XC(LI+K+3)
      RIJ2=RIJ2+XAUX(K)**2
      END DO
      END IF
      RIJ2IN=1./RIJ2
*       Introduce the inverse distances.
      LRI=LRI+1
      RINV(LRI)=SQRT(RIJ2IN)
*
      FM=MIJ(I,J)*RINV(LRI)
      UNC=UNC+FM
      FM=FM*RIJ2IN
*       Fij attraction.
      DO K=1,3
      FAUX(K)=-FM*XAUX(K)
      END DO
*       Add the contribution to interactions depending on Rij.
      DO IK=I,J-1
      L=3*(IK-1)
      DO K=1,3
      FNC(L+K)=FNC(L+K)+FAUX(K)
      END DO
      END DO
      END DO
      END DO
*
*	Evaluate UQ & TP.
      DO I=1,N-1
      L1=3*(I-1)+1
      KS=4*(I-1)
      KS1=KS+1
      CALL QFORCE(Q(KS1),FNC(L1),UQ(KS1))
      CALL VECTOR(Q(KS1),AK(KS1),TP(KS1))
*       The * operation of eq. (84).
      AK(KS+4)=-AK(KS+4)
      CALL VECTOR(P(KS1),AK(KS1),TQ(KS1))
*
      DO K=1,4
      UQ(KS+K)=UQ(KS+K)-2.0D0*MKK(I)*Q(KS+K)*RINV(I)**2
      TQ(KS+K)=TQ(KS+K)-4.D0*DK(I)*Q(KS+K)
      END DO
      END DO
*	NOTE: The division by R above (in TP & TQ) is delayed.
*
*	Proceed to final evaluation of derivatives (90)->(94).
      UPOT=UC+UNC
      G=1./(TKIN+UPOT)
      H=TKIN-UPOT
*       Reset EJUMP after significant slow-down change (Seppo's procedure).
*     IF (KJUMP) THEN
*         EJUMP = H - ENERGY
*         KJUMP = .false.
*     END IF
      GAMMA=(H-(ENERGY+EJUMP))*G
*
      GT= (1.-GAMMA)*G
      GU=-(1.+GAMMA)*G
*
      DO I=1,N-1
      KS=4*(I-1)
*       Apply the division by R here (to TP & TQ).
      GToverR=GT*RINV(I)
*       NOTE: TP & TQ never get 'correct' values (thus TP = R*TPtrue).
      DO K=1,4
      DQ(KS+K)=GToverR*TP(KS+K)
      DP(KS+K)=-GToverR*TQ(KS+K)-GU*UQ(KS+K)
      END DO
      END DO
      DT=G
      TPR=G
      DO K=1,3
      DX(K)=CMV(K)*G
      DV(K)=FCM(K)*G
      END DO
*	Evaluate E'.
      DE=0.0
      IF (ISKIP.EQ.0) THEN
      DO I=1,N-1
      L1=3*(I-1)+1
      KS=4*(I-1)
      KS1=KS+1
      CALL QFORCE(Q(KS1),FXTNL(L1),FAUX)
      DO K=1,4
      DE=DE+DQ(KS+K)*FAUX(K)
      END DO
      END DO
      END IF
*
*       Copy the time derivative for step control and reset indicator.
*     IF (IDER.GT.0) THEN
*         TPR = G
*         IDER = 0
*     END IF
*
*       Check osculating pericentre of closest pair (first call only).
      IF (ICALL.EQ.0) GO TO 50
      IF (JC.GT.0) GO TO 10
*
*       Perform a fast pericentre search (saves unnecessary complications).
      RM = 0.0
      DO 5 I = 1,N-1
          IF (RINV(I).GT.RM) THEN
              RM = RINV(I)
              IM = I
          END IF
    5 CONTINUE
*
      K1 = INAME(IM)
      K2 = INAME(IM+1)
      K = 4*(IM - 1) + 1
      CALL PERI(Q(K),DQ(K),TPR/KSCH(IM),M(K1),M(K2),QPERI)
*
*       Switch off indicator and exit for large pericentre.
      IF (QPERI.GT.4.0*MAX(SIZE(K1),SIZE(K2))) THEN
          ICALL = 0
          GO TO 50
      END IF
*
*       Examine all close pair indices K1 & K2 in chain vector INAME.
   10 QPMIN = 100.0
      RPMIN = 100.0
      IM0 = 1
      DO 20 IM = 1,N-1
          K1 = INAME(IM)
          K2 = INAME(IM+1)
*
*       Set distance of nearest perturber (note extra test for N > 4).
          IF (IM.EQ.1) THEN
              RP = 1.0/RINV(2)
          ELSE IF (IM.EQ.2.AND.N.GT.3) THEN
              RP = MIN(1.0/RINV(1),1.0/RINV(3))
          ELSE IF (IM.EQ.3.AND.N.GT.4) THEN
              RP = MIN(1.0/RINV(2),1.0/RINV(4))
          ELSE
              RP = 1.0/RINV(IM-1)
          END IF
          RPMIN = MIN(RP,RPMIN)
*
*       Determine pericentre for small perturbations by Mikkola's algorithm.
          GI = (1.0/(RINV(IM)*RP))**3
          IF (GI.LT.0.005) THEN
              K = 4*(IM - 1) + 1
              CALL PERI(Q(K),DQ(K),TPR/KSCH(IM),M(K1),M(K2),QPERI)
          ELSE
              QPERI = 1.0/RINV(IM)
          END IF
*
*       Compare pericentre with previous mutual distances (note symmetry).
          RIJ(K1,K2) = MIN(RIJ(K1,K2),QPERI)
          RIJ(K2,K1) = MIN(RIJ(K2,K1),QPERI)
*
*       Save indices for smallest pericentre and switch off indicator.
          IF (QPERI.LT.QPMIN.AND.IM.NE.IMCIRC) THEN
              G0 = GI
              QPMIN = QPERI
              RP0 = RP
              K10 = K1
              K20 = K2
              IM0 = IM
              ICALL = 0
          END IF
   20 CONTINUE
*
*       Check smallest pericentre of current chain and copy latest value.
      RCOLL = MIN(RCOLL,QPMIN)
      QPERI = QPMIN
*
*       Specify KS index and closest pair indices.
      IM = IM0
      KS1 = 4*(IM - 1) + 1
      K1 = K10
      K2 = K20
*
*       Save TPR and current configuration (for EREL, TRANSK & CHAIN).
      TKPR = TPR
      DO 25 I = 1,N-1
          KS = 4*(I - 1)
          DO 24 J = 1,4
              QK(KS+J) = Q(KS+J)
              PK(KS+J) = P(KS+J)
   24     CONTINUE
   25 CONTINUE
*
*       Check for tidal two-body interaction or stellar collision.
      IF (QPMIN.LT.4.0*MAX(SIZE(K1),SIZE(K2)).OR.ITER.GT.0) THEN
*
*       Terminate after 25 iterations (convergence problem).
          ITER = ITER + 1
          IF (ITER.GE.25) THEN
              if(rank.eq.0)WRITE (6,26)  ITER, IM, G0, QPMIN, ECC
   26         FORMAT (' WARNING!    NO CONVERGENCE    # IM GI QP ECC ',
     &                                                I5,I4,1P,3E10.2)
              JC = 0
              KCOLL = .false.
              ICOLL = -1
              ITER = 0
              IMCIRC = 0
              GO TO 50
          END IF
*
*       Exit if collision candidate distance is not the smallest.
          RB = 1.0/RINV(IM)
*         IF (RB.GT.1.001*RPMIN) THEN
*             JC = 0
*             DSC = 1.0
*             GO TO 50
*         END IF
*
*       Convert Q' to standard KS with T' = R and form radial velocity R'.
          RPR = 0.0D0
          DO 30 J = KS1,KS1+3
              UPR(J) = DQ(J)*RB*KSCH(IM)/TPR
              RPR = RPR + 2.0D0*Q(J)*UPR(J)
   30     CONTINUE
*
*       Determine small semi-major axis from non-singular expressions.
          CALL EREL(IM,EB,SEMI)
*
*       Exclude circularized binaries and consider second smallest peri.
          ECC = 1.0 - QPMIN/SEMI
*       Skip exclusion in GR case (denoted by large VSTAR1).
          IF (ECC.LT.0.002.AND.IMCIRC.EQ.0.AND.VSTAR1.LT.100.0) THEN
              ITER = 0
              IMCIRC = IM
              GO TO 10
          END IF
*
*       Temporary exit because SEMI < 0 does not converge (tested OK 3/99).
*         IF (SEMI.LT.0.0) THEN
*             JC = 0
*             GO TO 50
*         END IF
*
*       Obtain pericentre time interval from elements & Kepler's equation.
          MB = M(K1) + M(K2)
          CALL TPERI(SEMI,Q(KS1),UPR(KS1),MB,TPER)
*
*       Activate collision indicator & B-S step selector (first time).
          IF (ICOLL.EQ.0) THEN
              ICOLL = -1
              JC = 1
              KCOLL = .true.
          ELSE
*
*       Check convergence: radial velocity < 1.0E-09 parabolic velocity.
              IF (ABS(RPR).LT.1.0E-09*SQRT(2.0D0*MB*RB)) THEN
*       Reset B-S step selector and copy chain index to collision indicator.
                  JC = 0
                  ICOLL = IM
                  EBS = EB
                  ITER = 0
                  IMCIRC = 0
              END IF
          END IF
*
*       Set regularized step for DIFSY1 using accelerated iteration rate.
*         IF (RB.GT.2.0*QPMIN) THEN
*             DSC = 2.0*ABS(TPER)/TPR
*         ELSE
*             DSC = ABS(TPER)/TPR
*         END IF       
*
*       Evaluate regularized pericentre time (Stiefel & Scheifele, p. 85).
          HI = -0.5*MB/SEMI
*       Note use of Seppo's sign convention (TPER > 0 after peri).
          DSC = 2.0D0*(HI*TPER - 0.5D0*RPR)/MB
*       Scale from KS to chain time derivative.
          DSC = DSC*RB/TPR
*
*       Ensure negative step beyond pericentre (case RPR > 0 works OK).
          IF (JC.GT.0.AND.RPR.LT.0.0D0) THEN
              DSC = ABS(DSC)
          END IF
*       Restore step to dummy value at the end (not used).
          IF (JC.EQ.0) DSC = 1.0
*
*       Switch off iteration on large perturbation after five tries.
          IF (G0.GT.0.005.AND.ITER.GT.5) THEN
              JC = 0
              KCOLL = .false.
              ICOLL = -1
              ITER = 0
              IMCIRC = 0
*       Avoid apocentre region of secondary binary (algorithmic confusion).
          ELSE IF (RB.GT.SEMI.AND.IMCIRC.GT.0) THEN
              JC = 0
              ITER = 0
              IMCIRC = 0
              ICOLL = 0
              KCOLL = .false.
          END IF
      ELSE
          JC = 0
*       Enforce restart from saved variables on failed iteration.
          KCOLL = .false.
          ICOLL = -1
          ICALL = 0
          ITER = 0
      END IF
*
*       Increase function call counter.
   50 NFN = NFN + 1
*
      RETURN
*
      END
