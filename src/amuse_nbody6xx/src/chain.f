      SUBROUTINE CHAIN(ISUB)
*
*
*       Perturbed chain regularization. 
*       -------------------------------
*
*       Method of Mikkola & Aarseth, Celestial Mechanics 47, 375.
*       .........................................................
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      REAL*8   G0(3),Y(NMX8),R2(NMX,NMX),KSCH
      INCLUDE "mpif.h"
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      COMMON/MPIDAT/group,rank,ierr,isize,status
      INTEGER  IJ(NMX),IOLD(NMX)
      LOGICAL  CHECK,KSLOW,KCOLL,stopB,ICASE
      COMMON/SLOW1/   TK2(0:NMX),EJUMP,KSCH(NMX),KSLOW,KCOLL
      COMMON/SLOW2/   stepl,stopB
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NMX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NMX,5),ISYS(5)
      COMMON/INTFAC/  LX,LE,LP,LV,LT,J10,NHALF2
      COMMON/ECHAIN/  ECH
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
      COMMON/CALLS/  TPR,TKPR,STEP,IDER,ICALL,NFN,NREG,ITER,IMCIRC
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIJ(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/BSSAVE/  EP(4),DSC,FACM,TFAC,ITFAC,JC
      COMMON/EBSAVE/  EBS
      COMMON/KSAVE/  K1,K2
      COMMON/SLOW3/  GCRIT,KZ26
      COMMON/SWCALL/ NCALL
      EXTERNAL CHMOD
      SAVE

*
*
*       Main variables for chain regularization
*       ***************************************
*
*       -------------------------------------------------------------------
*       CHTIME  Local physical time (from CHAIN code).
*       CM(1-7) Coordinates, velocities & total mass of system.
*       CM(8)   Total energy of N-body system (copied in routine CHINIT).
*       CM(9)   Binary energy change (copied to CHCOLL in CHTERM).
*       ECH     Total energy of perturbed system (N-body interface).
*       ENERGY  Total energy.
*       EPS     Tolerance for DIFSY1 (1.0E-10 is recommended).
*       ICALL   Pericentre indicator (activated if R(IMIN) < EPSR2**0.5).
*       ICOLL   Collision indicator (activated in routine DERQP).
*       IDER    Indicator for setting time derivativative on first call).
*       IPERT   Perturbation indicator (=0: once per call; =1: every call).
*       I1-I4   Indices of final configuration (I1 & I2 is closest binary).
*       KSCH    Slow-down factor (= 1 when inactive).
*       KSLOW   Slow-down indicator (= .true. when active).
*       KZ27    Tidal dissipation option.
*       KZ30    Diagnostic option (copied in CHINIT; full output if > 1).
*       M       Particle mass (CM(7) = sum M(I), I = 1,N).
*       NAMEC   Particle identity (initialized to global name).
*       NFN     Number of function calls.
*       NPERT   Number of perturbers (for diagnostic output).
*       NREG    Number of regularization switches.
*       NSTEP1  Number of DIFSY1 calls.
*       P       Regularized momenta.
*       Q       Regularized coordinates.
*       RINV    Inverse chain distances.
*       RCOLL   Minimum two-body separation.
*       RGRAV   Gravitational radius ((sum M(I)*M(J))/ABS(ENERGY)).
*       RMAXS   Maximum size of unperturbed configuration.
*       RMAXC   Maximum of unperturbed & initial size (+ 20 per cent).
*       RSUM    Sum of all chain distances.
*       STEPS   Current integration interval (set by routine INTGRT).
*       TCR     Local crossing time ((sum M(I))**2.5/ABS(2*ENERGY)**1.5).
*       TIMEC   Local physical time in scaled units (= CHTIME).
*       TMAX    Maximum integration time (based on c.m. step).
*       TS      Global physical time at each return to N-body integration.
*       X       Particle coordinates (X(1), X(2), X(3) is X, Y, Z).
*       V       Velocity components.
*       -------------------------------------------------------------------
*
*
*       Save termination indicator and check for restart.
      ITERM = ISUB
      IF (ISUB.GT.0) THEN
*       Choose small step for termination (routine PERMIT & INTGRT).
          IF (STEPS(ISUB).EQ.0.0D0) THEN
              STEP = 1.0D-06*STEP
              GO TO 20
          END IF
*       Update maximum prediction interval at start of every call.
          CALL TCHAIN(ISUB,TSMIN)
          STEPS(ISUB) = TSMIN
*       Synchronize next time interval with subsystem step.
          TMAX = TIMEC + STEPS(ISUB)
          GO TO 20
      END IF
*
*       Copy initial conditions from N-body COMMON and prepare chain.
      TIMEC = 0.0D0
      CALL CHINIT(ISUB)
*
*       Initialize diagnostic & COMMON variables.
      RCOLL = 100.0
      QPERI = 100.0
      ZMI1 = 100.0
      ZMI2 = 100.0
      ZMI3 = 100.0
      ECOLL1 = 0.0D0
      DSC = 1.0
      DO 2 J = 1,10
          DO 1 K = 1,10
              RIJ(J,K) = 100.0
    1     CONTINUE
    2 CONTINUE
      ICALL = 0
      ICOLL = 0
      ISYNC = 0
      NCIRC = 0
      NSTEPX = 0
      IDER = 0
      IPERT = 1
      ITER = 0
      IMCIRC = 0
      NEXT = 0
      NDISS1 = 0
      JC = 0
      NSTEP1 = 0
      NREG = 0
      NFN = 0
      NCALL = 0
      KCASE = 0
      J10 = 10
      IPREV = 0
      ICASE = .FALSE.
      NAMES(NMX,5) = 0
*
*       Specify the tolerance for DIFSY1.
*     EPS = 1.0E-10
      EPS = 1.0E-12
*
*       Initialize subsystem time and set dummy variable for DIFSY.
      CHTIME = 0.0D0
      STIME = 0.0D0
*
*       Initialize chain regularization (new case or modified chain).
   10 NEQ = 8*N
      IF (KCASE.GT.0) THEN
          TMAX = TIMEC + STEPS(ISUB)
          CHTIME = TIMEC
          NEXT = NSTEP1 + 1
          ZMI1 = 100.0
          ZMI2 = 100.0
          ZMI3 = 100.0
          RCOLL = 100.0
      END IF
      Y(NEQ) = CHTIME
*
*       Define indices for DIFSY1 calls to derivative routine.
      NC = N - 1
      LX = 4*NC + 1
      LE = 4*NC + 4
      LP = 4*NC + 5
      LV = 8*NC + 5
      LT = 8*NC + 8
      NHALF2 = (NEQ/2)*2
*
*       Ensure routine XTPERT uses correct subsystem index (ISYS(5) is free).
      ISYS(5) = ISUB
*
*       Evaluate constants of the motion and save TPR.
      CALL CONST(X,V,M,N,ENER0,G0,ALAG)
      TPR = 1.0/ALAG
      TPR0 = TPR + 1.0
*
*       Select chain indices and transform from chain to KS variables.
      CALL SELECT
      CALL TRANSQ
*
*       Evaluate inverse distances (for routine INVERT).
      DO 11 I = 1,N-1
          KS1 = 4*(I - 1) + 1
          RK = Q(KS1)**2 + Q(KS1+1)**2 + Q(KS1+2)**2 + Q(KS1+3)**2
          RINV(I) = 1.0/RK
   11 CONTINUE
*
*       Define total energy (including any c.m. motion & external potential).
      ENERGY = ENER0 - 0.5D0*MASS*(CMV(1)**2 + CMV(2)**2 + CMV(3)**2)
*
*       Copy whole input array (for DIFSY call).
      CALL YCOPY(Y)
*
*       Find sum of mass products and set current gravitational radius.
      SUM = 0.0D0
      DO 14 L = 1,N-1
          DO 12 K = L+1,N
              SUM = SUM + MIJ(L,K)
   12     CONTINUE
   14 CONTINUE
      RGRAV = SUM/ABS(ENERGY)
*
*       Impose initial size as upper limit in case of weakly bound systems.
      RGRAV = MIN(RGRAV,0.5*RSUM)
      RSLOW = 0.5*RGRAV
*       Check slow-down option for chain (KZ26 copied in CHINIT).
      IF (KZ26.LT.2) RSLOW = 0.0
*
*       Set current crossing time and specify initial step using T' = 1/L.
      TCR = MASS**2.5/ABS(2.0D0*ENERGY)**1.5
*
*       Determine the smallest two-body time-scale from parabolic orbit.
      CALL R2SORT(IJ,R2)
      I1 = IJ(1)
      I2 = IJ(2)
      RM = SQRT(R2(I1,I2))
      VP2 = 2.0*(M(I1) + M(I2))/RM
      TP = RM/SQRT(VP2)
      TSTAR = MIN(TP,TCR)
      STEPIN = EPS**0.2*TSTAR*ALAG
      SMALL = 1.0D-04*STEPIN
*       Save step in case of backwards last iteration or collision.
      SAVEIT = STEPIN
*     RUNAV = 0.0001*RGRAV
*
*       (Re-)Initialize slow-down variables and minimum stellar radius.
      RSTAR = 0.0
      do i = 1,n-1
          ksch(i) = 1.0d0
          RSTAR = MAX(SIZE(i),RSTAR)
      end do
      RSTAR = MAX(SIZE(n),RSTAR)
*       Set pericentre search distance (twice 4*r_max; discrete intervals).
      RSTAR = 8.0*RSTAR  ! no longer used.
      KSLOW = .false.
      KCOLL = .false.
      EJUMP = 0.0D0
      stopB = .false.
*       Specify slow-down criterion with conservative value for N > 4.
      GCRIT = 1.0D-05
      IF (N.EQ.4) GCRIT = 5.0D-06
      IF (N.GT.4) GCRIT = 1.0D-06
*
      IF (KZ30.GT.1) THEN
          if(rank.eq.0)
     &    WRITE (6,15)  N, NPERT, ENERGY, RSUM, RGRAV, TCR, RMAXS(ISUB),
     &                  (NAMEC(I),I=1,N)
   15     FORMAT (' NEW CHAIN   N NP E RSUM RGRAV TCR RMAXS NAM ',
     &                          2I4,F10.5,1P,4E9.1,0P,6I6)
      END IF
      IF (RSUM.LT.1.0E-05) EPS = 1.0D-12
*
*       Evaluate binary energies (initially or after membership change).
      CALL RECOIL(0)
*
*       Assign integration step (initial case or modified membership).
      IF (KCASE.EQ.0) THEN
          STEP = STEPIN
      ELSE
          IPERT = 1
          KCASE = 0
          STEP = 0.01*STEPIN
          CALL TCHAIN(ISUB,TSMIN)
*       Copy current value in case STEPS = 0 which leads to termination.
          STEPS(ISUB) = TSMIN
*       Ensure termination if only two particles left or positive energy.
          IF (N.EQ.2.OR.ENERGY.GT.0.0) GO TO 70
      END IF
*
*       Perform regularized integration until termination or modification.
   20 RMAXC = 3.0*RGRAV
*       Replace perturbed boundary radius with MIN(3*RGRAV,2*RSUM).
      RMAXC = MIN(RMAXC,2.0*RSUM)
   21 STEP0 = STEP
      TIME0 = CHTIME
*
*       Modify DIFSY integration order according to relative perturbation.
*     JDIF = 10 - (1.0D+04*GPERT)**0.4
*     IF (JDIF.LT.J10) THEN
*         J10 = MAX(4,JDIF)
*     ELSE IF (JDIF.GT.J10) THEN
*         J10 = MIN(10,J10+1)
*     END IF
*
*       Obtain maximum step from two-body iteration based on c-functions.
      IF (ICOLL.EQ.0.AND.ICALL.EQ.0.AND.STEPS(ISUB).GT.0.0D0) THEN
*       Set slightly large interval to avoid small steps and check STEPS.
          DT = 1.01*TMAX - CHTIME
          IF (DT.LT.0.0) DT = ABS(DT)
          DT = MIN(1.01*STEPS(ISUB),DT)
*       Perform iteration of DTAU based on dominant two-body motion.
          CALL INVERT(DT,DTAU)
*       Increase step by 10% during expansion and check convergence value.
          IF (TPR.GT.TPR0) DTAU = 1.1*DTAU
*       Avoid choosing a negative step from previous backward integration.
          IF (STEP.GT.0.0) THEN
              STEP = MIN(DTAU,STEP)
          ELSE
              STEP = MIN(DTAU,SAVEIT)
          END IF
      END IF
*
*       Copy restart step after tidal dissipation or collision/failure.
      IF (ICASE.AND.IPREV.GT.0) THEN
          STEP = SAVEIT
          STEP0 = STEP
          IPREV = 0
          ICASE = .FALSE.
      END IF
*
*       Check time-step rejuvenation during standard integration.
      IF (MOD(NSTEP1,1000).EQ.0.AND.ICALL.EQ.0.AND.ICOLL.EQ.0) THEN
          CALL CONST(X,V,M,N,ENER0,G0,ALAG)
          SMALL = 1.0D-04*EPS**0.2*TSTAR*ALAG
*       Increase by modest factor to avoid danger of zero step.
          IF (STEP.LT.SMALL) STEP = SMALL
*         if(rank.eq.0)
*    &    WRITE (6,18)  NSTEP1, SMALL, STEP, (1.0/RINV(K),K=1,N-1)
*  18     FORMAT (' STEP INCREASE    # SM S R ',I8,1P,2E9.1,2X,5E10.2)
      END IF
*
*       Advance the solution one step.
      CALL DIFSY1(NEQ,EPS,STEP,STIME,Y)
*
*       Copy current physical time and save COMMON variables.
      CHTIME = Y(NEQ)
      TIMEC = CHTIME
      CALL YSAVE(Y)
      CALL TRANSX
      NSTEP1 = NSTEP1 + 1
      TPR0 = TPR
*
*       Save new step during standard integration for subsequent restart.
   22 IF (ICALL.EQ.0.AND.ICOLL.EQ.0) THEN
          SAVEIT = STEP
          ICASE = .TRUE.
*       Reset collision indicator and activate IPREV on failed iteration.
      ELSE IF (ICOLL.LT.0.AND..NOT.KCOLL) THEN
          ICOLL = 0
          IPREV = 1
          ICASE = .TRUE.
          KCOLL = .FALSE.
          NEXT = NSTEP1 + 1
          STEP = SAVEIT
*       Restore the original configuration and copy to input array.
          IF (IMCIRC.GT.0) ISYNC = 1
          TPR = TKPR
          DO 24 I = 1,N-1
              KS = 4*(I - 1)
              DO 23 J = 1,4
                  Q(KS+J) = QK(KS+J)
                  P(KS+J) = PK(KS+J)
   23         CONTINUE
   24     CONTINUE
          CALL YCOPY(Y)
          GO TO 21
      END IF
*
*       Check slow-down procedure for small minimum distance.
      IF (KZ26.GE.2.AND.(.not.stopB.or.KSLOW)) then
          RM = 0.0
          DO 25 I = 1,N-1
              IF (RM.LT.RINV(I)) THEN
                  RM = RINV(I)
                  IM = I
              END IF
   25     CONTINUE
          RM = 1.0/RM
          IF (RM.LT.RSLOW.AND.ICOLL.EQ.0) THEN
              CALL SLOW
          END IF
      END IF
*
*       Include possible emergency stop of slow-down (DERQP; inactive).
      if (stopB) then
          CALL YCOPY(Y)
          stopB = .false.
          IT = 0
   28     CALL DIFSY1(NEQ,EPS,stepl,STIME,Y)
          if (.not.stopB) then
              CALL YSAVE(Y)
          else if (IT.LT.5) then
              IT = IT + 1
              go to 28
          end if
          CALL SLOW
          go to 20
      end if
*
      IF (STEP.GT.10.0*STEP0.AND.ICOLL.EQ.0) THEN
          STEP = 10.0*ABS(STEP0)
      ELSE IF (STEP.EQ.0.0D0) THEN
          if(rank.eq.0) WRITE (6,*) ' Stepsize = 0!', char(7)
          STOP
      END IF
*
      IF (rank.eq.0.and.KZ30.GT.2) THEN
          WRITE (6,30)  STEP, TMAX-CHTIME, GPERT, (1.0/RINV(K),K=1,N-1)
   30     FORMAT (' CHAIN:   STEP TM-CHT G R  ',1P,8E9.1)
      END IF
*
*       Determine two-body distances for stability test and collision search.
      IF (CHTIME.GT.TIME0.AND.JC.EQ.0) THEN
          RM = 0.0
          RC2 = 0.0
          DO 35 K = 1,N-1
*       Find minimum separation for stability test and save chain index.
              IF (RINV(K).GT.RM) THEN
                  RM = RINV(K)
                  IX = K
              END IF
              RC2 = RC2 + RINV(K)**2
   35     CONTINUE
*       Update sum of 1/R^2 during forward integration (R^2 before 12/99).
          RM = 1.0/RM
          ZMI = RC2
          ZMI1 = ZMI2
          ZMI2 = ZMI3
          ZMI3 = ZMI
*         RUNAV = 0.9*RUNAV + 0.1*RM
*       Set search distance for closest separation.
          I1 = INAME(IX)
          I2 = INAME(IX+1)
*       Replace sum of radii by maximum value (08/08).
*         SX = SIZE(I1) + SIZE(I2)
          SX = MAX(SIZE(I1),SIZE(I2))
*       Turn off circular orbit indicators for overlapping stars.
          IF (IMCIRC.GT.0) THEN
*        Prevent repeated switching (old test RM < 2*MAX(SIZE); 07/08).
              IF (RM.LT.SX) THEN
                  SX = 0.5*RM
                  NCIRC = 0
                  ISYNC = 0
                  IMCIRC = 0
              ELSE
                  SX = 0.0
              END IF
          END IF
*       See whether another star is a more likely collider (compare SIZE/R).
          RY = RM
          DO 36 K = 1,N-1
              IF (K.NE.IX) THEN
                  J1 = INAME(K)
                  J2 = INAME(K+1)
                  SY = MAX(SIZE(J1),SIZE(J2))
                  IF (SY*RINV(K).GT.SX/RY) THEN
                      SX = SY
                      RY = 1.0/RINV(K)
                  END IF
              END IF
   36     CONTINUE
*       Include factor of 3 for criterion QPMIN < 4*MAX(SIZE(K1),SIZE(K2)).
          SX = 3.0*SX
          GX = 0.0
*       Compare sum of 1/R**3 neighbouring perturbations with dominant term.
          DO 38 K = 1,N-1
              IF (IABS(K-IX).EQ.1) THEN
                  GX = GX + RINV(K)**3
              END IF
   38     CONTINUE
*       Adopt safety ratio 0.01 for initializing close interaction search.
          PMAX = 0.01*RINV(IX)**3
      ELSE
          GX = 0.0
          PMAX = 1.0
      END IF
*
*       Switch on search indicator during iteration or just after pericentre.
      IF (ICOLL.LT.0) ICALL = 1
      IF (RY.LT.SX.AND.NSTEP1.GT.NEXT) THEN
          IF (ZMI3.LT.ZMI2.AND.ZMI2.GT.ZMI1) THEN
*       Skip pericentre determination on large perturbation (avoids looping).
              IF (ISYNC.EQ.0.AND.NCIRC.LT.25.AND.GX.LT.PMAX) THEN
                  ICALL = 1
                  IF (KZ26.GE.2.AND.KSLOW) THEN
                      GSAVE = GCRIT
                      GCRIT = 0.0
                      NCIRC = NCIRC + 1
                      CALL SLOW
                      GCRIT = GSAVE
                  END IF
              END IF
          END IF
      END IF
*
*       Restore indicators for further attempts after 100 steps.
      IF (ISYNC.GT.0.AND.NSTEP1.GT.NSTEPX + 100) THEN
          NCIRC = 0
          ISYNC = 0
          IMCIRC = 0
          NSTEPX = NSTEP1
      END IF
*
*       Check for collision, capture or circularization during last step.
      IF (ICOLL.LE.0) GO TO 50
      IPREV = ICOLL
*
*       Restore the minimum configuration from DERQP.
      TPR = TKPR
      DO 45 I = 1,N-1
          KS = 4*(I - 1)
          DO 40 J = 1,4
              Q(KS+J) = QK(KS+J)
              P(KS+J) = PK(KS+J)
   40     CONTINUE
   45 CONTINUE
*
*       Delay next search by one step to avoid the same pericentre.
      NEXT = NSTEP1 + 1
      ICALL = 0
*
*       Transform to physical variables.
      CALL TRANSX
*
*       Distinguish between tidal energy loss and collision (ICOLL holds IM).
      IM = ICOLL
      ICOLL = 0
*
*       Check for capture, circularization or collision.
      IF (QPERI.LT.4.0*MAX(SIZE(K1),SIZE(K2))) THEN
*
*       Distinguish between tidal energy loss and collision.
          J1 = K1
          J2 = K2
          IF (SIZE(K2).GT.SIZE(K1)) THEN
              J1 = K2
              J2 = K1
          END IF
*
*       Avoid tidal dissipation for circularized orbit.
          ICIRC = 0
          IF (KZ27.GT.0.AND.EBS.LT.0.0) THEN
              SEMI = -0.5*M(K1)*M(K2)/EBS
              ECC = 1.0 - QPERI/SEMI
              IF (ECC.GT.0.0022) ICIRC = 1
          END IF
*
*       Adopt collision criterion of Kochanek (Ap.J. 385, 684, 1992)
          IF (KZ27.LE.2) THEN
              FAC = 0.5*(M(J1) + M(J2))/M(J1)
              RCR = 1.7*FAC**0.3333*SIZE(J1)
          ELSE
              CLIGHT = 3.0D+05/VSTAR1
              RCR = 6.0*(M(J1) + M(J2))/CLIGHT**2
          END IF
*
          IF (QPERI.LT.RCR) THEN
*       Obtain global coordinates & velocities (ITERM < 0 denotes collision).
              ITERM = -1
              ISUB = -ISUB
              CALL CHTERM(ISUB)
*
*       Combine internal energy and external c.m. kinetic energy.
              HC = ENERGY + 0.5D0*MASS*(CM(4)**2 + CM(5)**2 + CM(6)**2)
*
*       Determine dominant two-body energy from non-singular expressions.
              CALL EREL(IM,EBS,SEMI)
*
*       Form composite body and terminate with new KS or continue chain.
              CALL CMBODY(HC,N)
*
*       Include chain restart (denoted by N < 0) for CE with > 3 members.
              IF (N.LT.0) THEN
                  N = -N
                  CALL CHINIT(ISUB)
                  GO TO 10
              END IF
*       Terminate for one single body after mass-less collision/coalescence.
              IF (N.EQ.0) THEN
                  ITERM = -1
                  GO TO 100
              END IF
*       Re-initialize with reduced membership NCH.
              N = N - 1
              CALL CHINIT(ISUB)
*       Activate indicator for new chain (terminates at once with N = 2).
              KCASE = 1
              GO TO 10
          ELSE IF (ICIRC.GT.0.AND.KZ27.GT.0.AND.ISYNC.EQ.0) THEN
*       Modify regularized variables due to tidal dissipation.
              CALL QPMOD(IM,ITERM)
*
*       Check reduction (ITERM = -2) or termination (ITERM = -1; N <= 4).
              IF (ITERM.EQ.-2) THEN
                  IF (N.EQ.2) GO TO 70
                  GO TO 10
              END IF
*       Re-determine the gravitational radius (cf. escape delay criterion).
              RGRAV = SUM/ABS(ENERGY)
              CALL YCOPY(Y)
              IF (ITERM.EQ.-1.AND.N.LE.4) THEN
                  CHTIME = Y(NEQ)
                  CALL YSAVE(Y)
                  CALL TRANSX
                  TIMEC = CHTIME
                  ECH = ENERGY
                  GO TO 70
              END IF
              GO TO 21
          ELSE
*       Allow for case of too long circularization time.
              ISYNC = 1
              ICOLL = -1
              KCOLL = .FALSE.
              GO TO 22
          END IF
      ELSE
*       Include case of large pericentre (first or second distance).
          IF (IMCIRC.GT.0) THEN
              ISYNC = 1
              ICOLL = -1
              GO TO 22
          END IF
          CALL YCOPY(Y)
          GO TO 21
      END IF
*
*       Check switching condition (Note: RINV now updated after switch).
   50 ISW = 0
      CALL SWCOND(CHECK)
      IF (CHECK) THEN
          DO 52 I = 1,N
              IOLD(I) = INAME(I)
   52     CONTINUE
          CALL SWITCH(Y)
          ISW = 1
*       See whether the chain name list is unchanged (even if CHECK =.TRUE.).
          DO 54 I = 1,N
              IF (IOLD(I).NE.INAME(I)) ISW = ISW + 1
   54     CONTINUE
          IF (ISW.GT.1) NREG = NREG + 1
*       Update slow-down if relevant (avoids perturbation jump).
          IF (KSLOW.AND.ISW.GT.1) THEN
              CALL SLOW
          END IF
      END IF
*
*       Check termination or strong perturbation (T > TMAX or GPERT > 0.01).
      IF (CHTIME.GT.TMAX.OR.GPERT.GT.0.01) THEN
          CALL YSAVE(Y)
          CALL TRANSX
          ECH = ENERGY
          TIMEC = CHTIME
          IF (rank.eq.0.and.KZ30.GT.2) THEN
              WRITE (6,55)  NSTEP1, T0S(ISUB)+TIMEC, TMAX-TIMEC,
     &                      (1.0/RINV(K),K=1,N-1)
   55         FORMAT (' CHAIN:  # T DTR R ',I5,F10.4,1P,6E9.1)
          END IF
*       Avoid checking after switch (just in case).
          IF (ISW.LE.1) THEN
              CALL CHMOD(ISUB,KCASE)
              IF (KCASE.GT.0) THEN
                  CALL RECOIL(1)
                  GO TO 10
              END IF
              IF (KCASE.LT.0) GO TO 60
          END IF
          IF (ISW.EQ.0.AND.N.EQ.3) THEN
              IF (RSUM.GT.4.0*RM) THEN
                  CALL CHSTAB(ITERM)
                  IF (ITERM.LT.0) GO TO 70
              END IF
          END IF
          IF (CHTIME.LT.TMAX.AND.STEPS(ISUB).GT.0.0D0) GO TO 21
          GO TO 60
      ELSE
*       Exit chain integration if time limit exceeded or STEPS = 0.
          IF (CHTIME.GT.TMAX) GO TO 60
*       Terminate on large step number and tiny STEP > 0.
          IF (NSTEP1.GT.100000.AND.STEP.LT.SMALL.AND.STEP.GT.0.0) THEN
              if(rank.eq.0)
     &        WRITE (6,58)  NSTEP1, N, STEP, (1.0/RINV(K),K=1,N-1)
   58         FORMAT (' ENFORCED CHAIN    # N STEP R',I8,I4,1P,9E9.1)
              GO TO 70
          END IF
          IF (STEPS(ISUB).GT.0.0D0) GO TO 21
      END IF
*
*       Check hierarchical stability condition for triple or quad.
   60 IF (N.EQ.3) THEN
          IF (RSUM.GT.4.0*RM) THEN
              CALL CHSTAB(ITERM)
              IF (ITERM.LT.0) GO TO 70
          END IF
      ELSE IF (N.EQ.4) THEN
          IF (RM.LT.0.1*RSUM) THEN
*       Find largest separation to distinguish triple or quad case.
              RX = 1.0
              DO 65 K = 1,N-1
                  RX = MIN(RX,RINV(K))
   65         CONTINUE
              RX = 1.0/RX
*       Check case of two binaries or degenerate triple and small binary.
              IF (RX.GT.0.7*RSUM) THEN
                  CALL CSTAB4(ITERM)
                  IF (ITERM.EQ.0.AND.IX.NE.2) THEN
                      CALL CSTAB2(ITERM)
                  END IF
*       Skip small middle distance (done by CSTAB3 called from CSTAB4).
              ELSE IF (RM.LT.0.01*RSUM.AND.IX.NE.2) THEN
                  CALL CSTAB2(ITERM)
              END IF
              IF (ITERM.LT.0) GO TO 70
          END IF
*       Reduce five/six-body system to triple if biggest binary < 0.04*RSUM.
      ELSE IF (N.GE.5) THEN
          CALL CSTAB5(ITERM)
          IF (ITERM.LT.0) GO TO 70
      END IF
*
*       See whether temporary or actual termination (continue if N > 5).
      IF (N.GT.5.OR.(KCASE.EQ.0.AND.STEPS(ISUB).GT.0.0D0)) GO TO 100
      IF (KCASE.LT.0) GO TO 70
      IF (TIMEC.GT.TMAX.AND.RSUM.LT.RMAXC) THEN
         IF (STEPS(ISUB).GT.0.0D0.OR.N.GT.4) GO TO 100
      END IF
*
*       Check for dominant binary.
   70 CALL RECOIL(2)
*
*       Set zero step to define termination (just in case).
      STEPS(ISUB) = 0.0D0
*
*       Copy Q & P for TRANSK via CHTERM & EREL (KZ26 > 2: rectification).
      IF (KZ26.GT.2) THEN
          DO 75 I = 1,N-1
              KS = 4*(I - 1)
              DO 72 J = 1,4
                  QK(KS+J) = Q(KS+J)
                  PK(KS+J) = P(KS+J)
   72         CONTINUE
              RIK = Q(KS+1)**2 + Q(KS+2)**2 + Q(KS+3)**2 + Q(KS+4)**2
              RINV(I) = 1.0/RIK
   75     CONTINUE
      END IF
*
*       Transform to global variables and begin new KS (I1 & I2), I3 & I4.
      CALL CHTERM(ISUB)
*
*       Activate termination index for routine INTGRT.
      ITERM = -1
*
      IF (rank.eq.0.and.KZ30.GT.1.AND.QPERI.LT.1.0) THEN
          WRITE (6,80)  RIJ(1,2), RIJ(1,3), RIJ(2,3), RCOLL, QPERI
   80     FORMAT (' END CHAIN   RIJ RCOLL QPERI ',1P,5E10.1)
      END IF
*
*       Update current time unless termination and set subsystem index.
  100 IF (ITERM.GE.0) TS(ISUB) = T0S(ISUB) + TIMEC
      ISUB = ITERM
*
      RETURN
*
      END
