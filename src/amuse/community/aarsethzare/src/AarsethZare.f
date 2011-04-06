*
*             T R I P L E
*             ***********
*
*
*       Three-body regularization program.
*       ----------------------------------
*
*       Method of Aarseth & Zare, Celestial Mechanics 10, 185.
*       ......................................................
*
*       Developed by Sverre Aarseth, IOA, Cambridge.
*       ............................................
*
      FUNCTION AarsethZare(TIME,BODY,POS,VEL,TCRIT,Etot,NREG)
c      FUNCTION AarsethZare(TIME,M,X,XDOT)
c      PROGRAM TRIPLE
*
      IMPLICIT  REAL*8  (A-H,M,O-Z)
      REAL*8 TIME, body(3), pos(3,3), vel(3,3), Etot, TCRIT
      INTEGER NREG
      COMMON/AZREG/  Q(8),P(8),R,R1,R2,ENERGY,M(3),X(3,3),XDOT(3,3),
     &               RCOLL,ERROR,C11,C12,C19,C20,C24,C25,NSTEPS,NAME(3)
      COMMON/CLOSE/  RIJ(3,3),ICALL
      COMMON/BSSAVE/  EP(4),TFAC,ITFAC,JC,NHALF2
      COMMON/NCALL/  NFN
      REAL*8  Y(17),ORDER(16),RK(3)
      INTEGER AarsethZare

      write (*,*) "Input:", TIME,BODY,POS,VEL,TCRIT,Etot,NREG
*
      AarsethZare = 1
      Etot = 0
*
*       COMMON variables
*       ****************
*
*       ------------------------------------------------------------------
*       C11     Inverse mass factor for DERQP (also C12,C19,C20,C24,C25).
*       ENERGY  Twice the initial total energy.
*       ERROR   Relative energy error at the end (ignore with T' = 1/U).
*       ICALL   Indicator for pericentre check (first function call only).
*       M       Particle mass.
*       NAME    Particle identity (initialized to 1,2,3).
*       NSTEPS  Number of DIFSY calls.
*       P       Regularized momenta.
*       Q       Regularized coordinates.
*       R       Distance between M(1) and M(2) (not regularized).
*       R1      Distance between M(1) and M(3).
*       R2      Distance between M(2) and M(3).
*       RCOLL   Minimum two-body separation (osculating pericentre).
*       RIJ     Minimum pairwise separations (osculating pericentre).
*       X       Particle coordinates (X(3,I) is Z-component).
*       XDOT    Velocity components (XDOT(3,I) is Z-component).
*       ------------------------------------------------------------------
*
*
      NEXP = 0
*       Read tolerance and termination time.
      TOL0 = 1.e-12
cccc  TCRIT = 100
c      READ (5,*)  TOL0, TCRIT
*
*       Generate next initial condition (place most active body in M(3)).
c      DO 2 I = 1,3
c          M(I) = body(i)
c          DO 3 K = 1,3
c             X(K,I) = pos(k,i)
c             XDOT(K,I) = vel(k,i)
c    3     CONTINUE
c    2 CONTINUE

      CALL DATA(body, pos, vel)
*
*       Save the tolerance.
      EPS = TOL0
*
*       Initialize diagnostic variables & counters.
      R12MIN = 100.0
      RMIN = 100.0
      RCOLL = 100.0
      DO 10 J = 1,3
          DO 5 K = 1,3
              RIJ(J,K) = 100.0
    5     CONTINUE
   10 CONTINUE
      NSTEPS = 0
      NREG = 0
      ICALL = 0
      NFN = 0
      JC = -1
*
*       Initialize local time & regularized time.
      TIME = 0.0D0
      TAU = 0.0D0
      Y(17) = 0.0D0
*       Specify the number of first-order equations for the integrator.
      N = 17
*
*       Obtain initial energy and transform to regularized variables.
      CALL TRANSF(1)
*
*       Define gravitational radius and pericentre check distance.
      RGRAV = (M(1)*M(2) + M(1)*M(3) + M(2)*M(3))/(0.5D0*ABS(ENERGY))
      RSTAR = 0.5*RGRAV
*
*       Form the two smallest distances (assume sensible reference body).
      R1 = Q(1)**2 + Q(2)**2 + Q(3)**2 + Q(4)**2
      R2 = Q(5)**2 + Q(6)**2 + Q(7)**2 + Q(8)**2
*
*       Set termination distance to maximum of separation and 10*RGRAV.
      RMAX0 = 1.1*MAX(R1,R2,10.0*RGRAV)
*
*       Specify the crossing time (also meaningful if ENERGY > 0).
      TCR = (M(1) + M(2) + M(3))**2.5/ABS(ENERGY)**1.5
      RAP = 10.0*RMAX0
      TNEXT = 0.0
*
*       Define a nominal crossing time for nearly parabolic systems. 
      RMAX = MAX(R1,R2)
      TSTAR = RMAX*SQRT(RMAX/(M(1) + M(2) + M(3)))
*
*       Determine the smallest two-body time-scale from parabolic orbit.
      IM = 1
      RM = R1
      IF (R2.LT.R1) THEN
          IM = 2
          RM = R2
      END IF
      VP2 = 2.0*(M(IM) + M(3))/RM
      TP = RM/SQRT(VP2)
      TSTAR = MIN(TP,TSTAR)
c      write (*,*) "TIMES=", TP, TSTAR
      TORB = TP
*     TCRIT = TORB*TCRIT
c     Reinstated: integrate for 100 orbits
c     Assure that TCRIT is properly defined if not given as input (Sept.2010).
      if (TCRIT.le.0) then
         TCRIT = TORB*100.0
      endif

*       Set TPR and initial step (ignore M1*M2/R if DT/DTAU = 1/U).
*     TPR = R1*R2/(M(1)*M(3)*R2 + M(2)*M(3)*R1)
      TPR = R1*R2/SQRT(R1 + R2)
      DTAU = MIN(TCR,TSTAR)*EPS**0.1/TPR
*
*       Try alternative expression for initial step (Seppo Mikkola).
*     rmxt=max(r1,r2)
*     dtau=.1*(EPS/1.0E-12)**0.1*rmxt*sqrt(rmxt/(m(1)+m(2)+m(3)) )/tpr

*       Initialize time constant & input array for the integrator.
      CONST = 0.0D0
      DO 20 K = 1,8
          CONST = CONST + Q(K)*P(K)
          Y(K) = Q(K)
          Y(K+8) = P(K)
   20 CONTINUE
*
*       Produce initial output.
      CALL TRANSF(2)
*
*       Specify reference values for each vector (Mikkola's method).
      DO 35 K = 1,4
          K1 = 4*(K - 1)
          SN1 = 0.25*(ABS(Y(K1+1)) + ABS(Y(K1+2)) + ABS(Y(K1+3)) +
     &                                              ABS(Y(K1+4)))
*       Initialize the reference vector.
          IF (NSTEPS.EQ.0) THEN
              ORDER(K1+1) = SN1
          END IF
*       Do not permit any small values because of relative tolerance.
          SN1 = 0.1*SN1 + 0.9*ORDER(K1+1)
          DO 34 L = 1,4
              ORDER(K1+L) = SN1
   34     CONTINUE
   35 CONTINUE
*
*       Advance the equations of motion by Bulirsch-Stoer integrator.
*     CALL DIFSY1(N,EPS,ORDER,DTAU,TAU,Y)
   30 CALL DIFSY1(N,EPS,DTAU,TAU,Y)
*
*       Copy regularized coordinates & momenta and obtain physical time.
c      SUMQP = 0.0D0
      DO 40 K = 1,8
          Q(K) = Y(K)
          P(K) = Y(K+8)
c          SUMQP = SUMQP + Q(K)*P(K)
*       Note that the momentum includes factor of 2 in AZ formulation.
   40 CONTINUE
*
*       Set explicit time (Baumgarte & Stiefel, 1974 & Aarseth, 1976).
c      TIME = (0.5D0*(SUMQP - CONST) - TAU)/ENERGY
       TIME = Y(17)
      write (*,*) "TIMES=", TP, TSTAR, TIME, X

*
*       Update relative distances (NB! not quite latest value).
      R1 = Q(1)**2 + Q(2)**2 + Q(3)**2 + Q(4)**2
      R2 = Q(5)**2 + Q(6)**2 + Q(7)**2 + Q(8)**2
*
*       Check minimum two-body separations and increase step counter.
      RMIN = MIN(RMIN,R)
      RM = MIN(R1,R2)
      R12MIN = MIN(R12MIN,RM)
      RMAX = MAX(R1,R2,R)
      NSTEPS = NSTEPS + 1
*
*       Check minimum two-body separations.
      RK(1) = R1
      RK(2) = R2
      RK(3) = R
*       Consider pairs 1-2, 1-3 & 2-3 with identified names.
      DO 44 K = 1,3
          DO 42 L = K+1,3
              I = NAME(K)
              J = NAME(L)
*       Use cyclic loop index (3,1,2) for distances R, R1 & R2.
              KK = K - 1
              IF (KK.EQ.0) KK = 3
              RIJ(I,J) = MIN(RIJ(I,J),RK(KK))
              RIJ(J,I) = MIN(RIJ(J,I),RK(KK))
   42     CONTINUE
   44 CONTINUE
*
*       Switch on search indicator inside RSTAR (reset in DERQP).
      IF (RM.LT.RSTAR) THEN
          ICALL = 1
      END IF
*
*       Use actual two-body force to decide on switching.
      F12 = (M(1) + M(2))/R**2
      F13 = (M(1) + M(3))/R1**2
      F23 = (M(2) + M(3))/R2**2
      IF (F12.LT.MAX(F13,F23)) GO TO 70
*
      IMIN = 1
*       Use a simple distance test to determine new reference body IMIN.
      IF (R2.LT.1.00001*R1) IMIN = 2
*
*       Transform to physical variables and rename the exchanged particles.
      CALL TRANSF(3)
*
      DO 50 K = 1,3
          TEMP1 = X(K,3)
          TEMP2 = XDOT(K,3)
          X(K,3) = X(K,IMIN)
          XDOT(K,3) = XDOT(K,IMIN)
          X(K,IMIN) = TEMP1
          XDOT(K,IMIN) = TEMP2
   50 CONTINUE
*
      TEMP1 = M(3)
      M(3) = M(IMIN)
      M(IMIN) = TEMP1
      NAME3 = NAME(3)
      NAME(3) = NAME(IMIN)
      NAME(IMIN) = NAME3
*
*       Transform back to regularized variables and initialize input array.
      CALL TRANSF(4)
      DO 60 K = 1,8
          Y(K) = Q(K)
          Y(K+8) = P(K)
   60 CONTINUE
*
*       Update regularization counter at the end of switching procedure.
      NREG = NREG + 1
*
*       Check termination criteria (TIME > TCRIT or RMAX > RMAX0).
   70 IF (TIME.GT.TCRIT) GO TO 90
*     IF (MAX(R1,R2).LT.RAP) GO TO 30
*     IF (TIME.LT.TNEXT) GO TO 30
      TNEXT = TNEXT + TORB
*
*       Obtain final output after transforming to physical variables.
      CALL TRANSF(2)
*
*       See whether the final configuration is well resolved.
      RATIO = R1/R2
*      IF (RATIO.GT.0.1.AND.RATIO.LT.10.0) GO TO 90
      IF (RATIO.GT.0.1.AND.RATIO.LT.10.0) GO TO 30
*
*       Set index of second binary component & escaper.
      IMIN = 1
      IF (R2.LT.R1) IMIN = 2
      IESC = 3 - IMIN
*
*       Evaluate orbital elements.
      RDOT = 0.0D0
      VREL2 = 0.0D0
      VESC = 0.0D0
      DO 75 K = 1,3
          RDOT = RDOT + (X(K,3) - X(K,IMIN))*(XDOT(K,3) - XDOT(K,IMIN))
          VREL2 = VREL2 + (XDOT(K,3) - XDOT(K,IMIN))**2
          VESC = VESC + (XDOT(K,3) - XDOT(K,IESC))**2
   75 CONTINUE
*
      RB = RM
      MB = M(3) + M(IMIN)
      SEMI = 2.0D0/RB - VREL2/MB
      VESC = SQRT(VESC)/SQRT(MB*ABS(SEMI))
*       Velocity of escaper w.r. to M(3) (scaled by binary velocity).
      SEMI = 1.0/SEMI
      E = SQRT((1.0 - RB/SEMI)**2 + RDOT**2/(SEMI*MB))
      EB = -0.5D0*M(3)*M(IMIN)/SEMI
      IF (TIME.LT.TORB) SEMI0 = SEMI
*
      IF (EB.LT.0.5*ENERGY) THEN
          EBE = EB/(0.5*ENERGY)
*       Binding energy scaled by total energy.
*
          WRITE (6,80)  TIME, NAME(IMIN), NAME(3), R1, R2, 
     &                  SEMI, E, NREG, ERROR, Etot
   80     FORMAT  ('TIME= ', E12.4, ' BINARY ',2I3,'  R1 =',1P,E9.1,
     &                  '  R2 =',E9.1,
     &                  '  A =',E12.4,'  E =',0P,F7.3,'  NREG =',I3,
     &                  '  DE/E =',E9.1,'  DEtot/Etot =',E9.1)
          CALL FLUSH(6)
c
c         copy the final parameters
          DO 85 L = 1,3
             I = NAME(L)
             BODY(I) = M(L)
             DO 84 K = 1,3
                POS(K,I) = X(K,L)
                VEL(K,I) = XDOT(K,L)
   84        CONTINUE
   85     CONTINUE

c         Last two lines added
          AarsethZare = 0
          return
          GO TO 90
      END IF
*     IF (TIME.LT.TCRIT.AND.ABS(SEMI-SEMI0).LT.0.5*SEMI0) GO TO 30
      IF (TIME.LT.TCRIT) GO TO 30
*
   90 NEXP = NEXP + 1
      CALL TRANSF(2)
*
c
c
*       Identify particles and copy  mass, coordinates and velocities.
      DO 95 L = 1,3
          I = NAME(L)
          BODY(I) = M(L)
          DO 94 K = 1,3
              POS(K,I) = X(K,L)
              VEL(K,I) = XDOT(K,L)
   94     CONTINUE
   95 CONTINUE
c
c
      RCOLL = MIN(RCOLL,R12MIN)
      Etot = Etot + ERROR
      WRITE (*,*)  "TORB+", TIME, TORB, RIJ(1,2), RIJ(1,3), RIJ(2,3)
      WRITE (6,100)  TIME/TORB, TIME, RMIN, R12MIN, ERROR, NSTEPS
  100 FORMAT  (/,'  T/TK =',F7.1, ' TIME=', F7.1, ' MIN(R) =',1PE8.1,
     &        '  MIN(R1,R2) =',E8.1,
     &                 '  DE/E =',E9.1,'  NSTEPS =',I8)
c      WRITE (6,105)  RIJ(1,2), RIJ(1,3), RIJ(2,3), RCOLL
c  105 FORMAT (' RIJ:  1-2 1-3 2-3 RCOLL ',1P4E10.2)
c      WRITE (6,110)  NSTEPS, NFN
c  110 FORMAT ('  NSTEPS =',I9,'  NFN =',I11)
*
*       Continue with the next experiment.
c      GO TO 1
      AarsethZare = 1
      return

*
      END
      include 'difsy1.f'
      include 'data.f'
      include 'peri.f'
      include 'block.f'
      include 'derqp.f'
      include 'transf.f'

