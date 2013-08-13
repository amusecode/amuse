      SUBROUTINE QUAD(ISUB)
*
*
*       Four-body regularization. 
*       -------------------------
*
*       Method of Mikkola & Aarseth, Celestial Mechanics 47, 375.
*       .........................................................
*
*
      IMPLICIT REAL*8  (A-H,O-Z)
      PARAMETER  (NCMAX=10)
      REAL*8  M,MIJ,MB,MB1,SAVEX(12),SAVEXD(12)
      LOGICAL  SWITCH,GTYPE,GTYPE0
      COMMON/CREG/  M(4),X(12),XD(12),P(12),Q(12),TIME4,ENERGY,EPSR2,
     &              XR(9),W(9),R(6),TA(6),MIJ(6),CM(10),RMAX4,TMAX,
     &              DS,TSTEP,EPS,NSTEP4,NAME4(4),KZ15,KZ27,NREG,NFN
      COMMON/TPR/   SWITCH,GTYPE,GTYPE0
      COMMON/CONFIG/  R2(4,4),I1,I2,I3,I4
      COMMON/CLOSE/  RIJ(4,4),RCOLL,QPERI,SIZE(4),ECOLL3,IP(4)
      COMMON/CCOLL/  QK(12),PK(12),ICALL,ICOLL,NDISS4
      COMMON/BSSAVE/  EP(4),DSC,FACM,TFAC,ITFAC,JC
      COMMON/KSAVE/  K1,K2
      COMMON/EBSAVE/  EBS
      SAVE
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)

*
*
*       Main variables for four-body regularization
*       *******************************************
*
*       -------------------------------------------------------------------
*       CM(1-7) Coordinates, velocities & total mass of four-body system.
*       CM(8)   Total energy of N-body system (copied in routine START4).
*       CM(9)   Binary energy change (copied to BBCOL in START4).
*       DS      Regularized time-step (rescaled after switching & new TPR).
*       ENERGY  Total energy.
*       EPS     Tolerance for DIFSY4 (1.0E-08 is recommended).
*       EPSR2   Distance criterion for stabilized time transformation.
*       ICALL   Pericentre indicator (activated if R(IMIN) < EPSR2**0.5).
*       ICOLL   Collision indicator (activated in routine DERQP4).
*       IND     Sorting index (R(IND(6)) is largest distance).
*       I1-I4   Indices of final configuration (I1 & I2 is closest binary).
*       KZ15    Diagnostic option (copied in START4; full output if > 1).
*       M       Particle mass (CM(7) = M(1) + M(2) + M(3) + M(4)).
*       NAME4   Particle identity (initialized to global name).
*       NFN     Number of function calls.
*       NREG    Number of regularization switches.
*       NSTEP4  Number of DIFSY4 calls.
*       P       Regularized momenta.
*       Q       Regularized coordinates.
*       QPERI   Minimum two-body separation.
*       R       Mutual distances (R(1), R(2), R(3) are regularized).
*       RCOLL   Minimum two-body separation.
*       RGRAV   Gravitational radius ((sum M(I)*M(J))/ABS(ENERGY)).
*       RIJ     Minimum pairwise separations.
*       RMAXS   Maximum size of unperturbed configuration.
*       RMAX4   Maximum of unperturbed & initial size (+ 20 per cent).
*       R2      Square mutual distances (R2(I1,I2) is closest pair).
*       TCR     Local crossing time ((sum M(I))**2.5/ABS(2*ENERGY)**1.5).
*       TIME4   Local physical time in scaled units.
*       TMAX    Maximum integration time (based on c.m. step).
*       TSTEP   Initial time-step for regularization (physical units).
*       X       Particle coordinates (X(1), X(2), X(3) is X, Y, Z).
*       XD      Velocity components.
*       -------------------------------------------------------------------
*
*
*       Save termination indicator and check for restart.
      ITERM = ISUB
      IF (ISUB.GT.0) THEN
*       Synchronize next time interval with c.m. step unless termination.
          TMAX = TIME4 + STEPS(ISUB)
          IF (STEPS(ISUB).LE.0.0D0) DS = 0.01*DS
          GO TO 30
      END IF
*
*       Obtain initial conditions from N-body COMMON.
      CALL START4(ISUB)
*
*       Specify tolerance & DSC for DIFSY4 (relative energy error < TOL).
      EPS = 1.0D-08
      EPS = 1.0D-12
      DSC = 1.0
*
*       Initialize diagnostic & COMMON variables.
      TIME4 = 0.0D0
      RCOLL = 100.0
      DO 10 J = 1,4
          DO 5 K = 1,4
              RIJ(J,K) = 100.0
    5     CONTINUE
   10 CONTINUE
      ICALL = 0
      ICOLL = 0
      NDISS4 = 0
      JC = 0
      ITFAC = 0
      NSTEP4 = 0
      IQ = 0
      ECOLL3 = 0.0D0
      NREG = 0
      NFN = 0
      ITRY = 0
      ISAVE = 0
      RATIO = 0.0
      SWITCH = .FALSE.
*
*       Evaluate the initial total energy (superseded; done in START4).
*     CALL NEWSYS(X,XD,M,4,ENERGY,GAM)
*
*       Initialize regularized coordinates, momenta & particle indices.
      CALL NEWREG
*
*       Define gravitational radius & average mass factor for routine DERQP4.
      RGRAV = (MIJ(1) + MIJ(2) + MIJ(3) + MIJ(4) + MIJ(5) + MIJ(6))/
     &                                                      ABS(ENERGY)
      FACM = RGRAV*ABS(ENERGY)/6.0
*       Set small distance criterion for stabilized time transformation.
      EPSR2 = (0.2*RGRAV)**2
*
*       Ensure that unperturbed boundary exceeds system size (+ 20 per cent).
      RMAX = MAX(SQRT(R2(I1,I3)),SQRT(R2(I2,I4)))
      RMAX4 = MAX(RMAXS(ISUB),1.2*RMAX)
*
*       Specify local crossing time.
      TCR = (M(1) + M(2) + M(3) + M(4))**2.5/ABS(2.0D0*ENERGY)**1.5
*
*       Define a nominal crossing time in near-parabolic cases.
      TSTAR = RMAX*SQRT(RMAX/CM(7))
*
*       Determine the smallest two-body time-scale from parabolic orbit.
      RM = SQRT(R2(I1,I2))
      VP2 = 2.0*(M(I1) + M(I2))/RM
      TP = RM/SQRT(VP2)
      TSTAR = MIN(TP,TSTAR)
*
*       Set initial time-step in physical units (also used by RCHAIN).
      TSTEP = MIN(TCR,TSTAR)*EPS**0.1
*
*       Convert physical step to regularized time.
      DS = TSTEP/(R(1)*R(2)*R(3))
*
*       Form initial binary energies (scaled by internal energy).
      VREL2 = 0.0D0
      VREL21 = 0.0D0
      DO 25 K = 1,3
          J1 = 3*(I1-1) + K
          J2 = 3*(I2-1) + K
          J3 = 3*(I3-1) + K
          J4 = 3*(I4-1) + K
          VREL2 = VREL2 + (XD(J1) - XD(J2))**2
          VREL21 = VREL21 + (XD(J3) - XD(J4))**2
   25 CONTINUE
      EB0 = 2.0D0/SQRT(R2(I1,I2)) - VREL2/(M(I1) + M(I2))
      EB10 = 2.0D0/SQRT(R2(I3,I4)) - VREL21/(M(I3) + M(I4))
      EB0 = -0.5D0*M(I1)*M(I2)*EB0/ENERGY
      EB10 = -0.5D0*M(I3)*M(I4)*EB10/ENERGY
*
*       Perform regularized integration until termination is reached.
   30 CALL RCHAIN(IQ)
*
*       See whether temporary or actual termination (otherwise R > RMAX4).
      IF (TIME4.GT.TMAX) THEN
         IF (STEPS(ISUB).GT.0.0D0) GO TO 100
      END IF
*
*       Transform to physical variables and order particle indices.
      SWITCH = .FALSE.
      CALL ENDREG
      CALL STATUS(X,J1,J2,J3,J4)
*
*       Save the current configuration.
      DO 32 J = 1,12
          SAVEX(J) = X(J)
          SAVEXD(J) = XD(J)
   32 CONTINUE
      TSAVE = TIME4
      DS0 = DS
*
*       Check for collision during last step.
      IF (ICOLL.LE.0) GO TO 40
*
*       Restore the minimum configuration.
      DO 35 K = 1,12
          Q(K) = QK(K)
          P(K) = PK(K)
   35 CONTINUE
*
*       Transform to physical variables.
      SWITCH = .FALSE.
      CALL ENDREG
*
*       Order particle indices of the collision configuration.
*     CALL STATUS(X,J1,J2,J3,J4)
*
*       Distinguish between tidal energy loss and collision (ICOLL holds IM).
      IM = ICOLL
      ICOLL = 0
      IF (QPERI.LT.4.0*MAX(SIZE(K1),SIZE(K2))) THEN
          IF (QPERI.LT.0.75*(SIZE(K1) + SIZE(K2))) THEN
*
*       Obtain global coordinates & velocities (ITERM < 0 denotes collision).
              ITERM = -1
              ISUB = -ISUB
              CALL START4(ISUB)
*
*       Combine internal energy and external c.m. kinetic energy.
              H4 = ENERGY + 0.5D0*CM(7)*(CM(4)**2 + CM(5)**2 + CM(6)**2)
*
*       Determine dominant two-body energy from non-singular expressions.
              CALL EREL4(IM,EBS,SEMI)
              DMINC = MIN(RCOLL,DMINC)
*
*       Form composite body and begin KS regularization of closest pair.
              CALL CMBODY(H4,4)
              GO TO 100
          ELSE
*
*       Modify regularized variables and check stability (ITERM < 0).
              CALL QPMOD4(IM,ITERM)
              SWITCH = .TRUE.
              CALL NEWREG
              ICALL = 0
              IQ = ITERM
              IF (ITERM.EQ.-1) GO TO 100
              GO TO 30
          END IF
      END IF
*
*       Form output diagnostics (only hard binary printed if KZ15 < 2).
   40 VREL2 = 0.0D0
      VREL21 = 0.0D0
      RDOT = 0.0D0
*
      DO 50 K = 1,3
          J1 = 3*(I1-1) + K
          J2 = 3*(I2-1) + K
          J3 = 3*(I3-1) + K
          J4 = 3*(I4-1) + K
          VREL2 = VREL2 + (XD(J1) - XD(J2))**2
          VREL21 = VREL21 + (XD(J3) - XD(J4))**2
          RDOT = RDOT + (X(J1) - X(J2))*(XD(J1) - XD(J2))
   50 CONTINUE
*
*       Evaluate orbital elements.
      RB = SQRT(R2(I1,I2))
      RB1 = SQRT(R2(I3,I4))
      R13 = SQRT(R2(I1,I3))
      R24 = SQRT(R2(I2,I4))
      MB = M(I1) + M(I2)
      SEMI = 2.0D0/RB - VREL2/MB
      SEMI = 1.0/SEMI
      ECC = SQRT((1.0D0 - RB/SEMI)**2 + RDOT**2/(SEMI*MB))
      EB = -M(I1)*M(I2)/(2.0D0*SEMI*ENERGY)
      ET = -M(I1)*M(I2)/(2.0D0*SEMI*CM(8))
      MB1 = M(I3) + M(I4)
      SEMI1 = 2.0D0/RB1 - VREL21/MB1
      SEMI1 = 1.0/SEMI1
      EB1 = -M(I3)*M(I4)/(2.0D0*SEMI1*ENERGY)
*       Scale binding energy of relative motion by initial value.
      E1 = (1.0 - EB - EB1)/EB0
*
*       Determine nearest and most distant binary perturber.
      IF (R13.LT.R24) THEN
          IM = I3
          RM = R13
          IMAX = I4
          RMAX = R24
      ELSE
          IM = I4
          RM = R24
          IMAX = I3
          RMAX = R13
      END IF
*
*       Estimate perturbations on smallest binary.
      GB = 2.0*MB1/MB*(RB/RM)**3
      IF (EB1.LT.0.0) GB = GB*M(IM)/MB1
      G4 = 2.0*M(IMAX)/(MB + M(IM))*(RM/RMAX)**3
      DB = (EB - EB0)/EB0
      IF (ITRY.LT.0) GO TO 70
*
*       Find the optimum configuration for two binaries or compact triple.
      ITRY = ITRY + 1
      IF (RM/SEMI.LT.RATIO) GO TO 60
      IF (RB.LT.SEMI.OR.RB1.LT.SEMI1) GO TO 60
*       Avoid new regularization of any binary near small pericentre.
      RATIO = RM/SEMI
      ISAVE = ITRY
   60 IF (RATIO.LT.5.0.AND.ITRY.LE.5.
     & OR.RATIO.LT.3.0.AND.ITRY.LE.10) THEN
          GO TO 30
      END IF
*
      IF (ISAVE.EQ.ITRY) GO TO 70
      ITRY = -1
*       Restore the maximum triple configuration and obtain indices.
      DO 65 J = 1,12
          X(J) = SAVEX(J)
          XD(J) = SAVEXD(J)
   65 CONTINUE
      TIME4 = TSAVE
      CALL NEWREG
      DS = 0.5*DS0
      GO TO 40
*
*       Postpone temination by one small step unless restart case.
   70 IF (STEPS(ISUB).GT.0.0D0) THEN
          STEPS(ISUB) = 0.0D0
          GO TO 100
      END IF
*
*       Print final binary if relative energy increase exceeds 0.1.
      IF (DB.GT.0.1) THEN 
          if(rank.eq.0)
     &    WRITE (6,80)  MB, SEMI, ECC, EB, GB, G4, RB1, EB1, E1, ET
   80     FORMAT (/,' QUAD BINARY','  MB =',F7.4,'  A =',1P,E8.1,
     &     '  E =',0P,F5.2,'  EB =',F7.4,'  GB =',1P,E8.1,'  G4 =',E8.1,
     &     '  RB1 =',E8.1,'  EB1 =',0P,F5.2,'  E1 =',F5.2,'  ET =',F6.3)
      END IF
*
*       Check optional print diagnostics of four-body regularization.
      IF (KZ15.GT.1) THEN
          TC = TIME4/TCR
          EC = ENERGY/CM(8)
          if(rank.eq.0)
     &    WRITE (6,90)  I1, I2, I3, I4, RB, R13, R24, RGRAV, TC, NSTEP4,
     &                  NREG, DB, EC
   90     FORMAT (/,' END QUAD   ',4I3,'  RB =',1P,E8.1,'  R13 =',E8.1,
     &              '  R24 =',E8.1,'  RG =',E8.1,'  TC =',0P,F5.1,'  #',
     &                 I5,I4,'  DB =',F5.2,'  EC =',F6.3)
      END IF
*
*       Form binary energy change (set in BBCOLL; routine START4).
      CM(9) = (EB - EB0 - EB10)*ENERGY
      IF (EB1.GT.0.0) CM(9) = CM(9) + EB1*ENERGY
*
*       Transform to global variables and begin new KS (I1 & I2), I3 & I4.
      CALL START4(ISUB)
*
*       Activate termination index for routine INTGRT.
      ITERM = -1
*
*       Update current time unless termination and set subsystem index.
  100 IF (ITERM.GE.0) TS(ISUB) = T0S(ISUB) + TIME4
      ISUB = ITERM
*
      RETURN
*
      END
