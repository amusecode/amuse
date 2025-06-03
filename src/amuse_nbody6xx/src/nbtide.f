      SUBROUTINE NBTIDE(I,J,RJMIN2)
*
*
*       Close two-body interaction.
*       ---------------------------
*
      INCLUDE 'common6.h'
      REAL*8  XI(3),XJ(3),VI(3),VJ(3),VCM(3),VR(3)
*
*
*       Copy coordinates of close bodies and predict velocity for body #J.
      RDOT = 0.0
      VREL2 = 0.0
      DT = TIME - T0(J)
      DO 10 K = 1,3
          XI(K) = X(K,I)
          XJ(K) = X(K,J)
          VI(K) = X0DOT(K,I)
          VJ(K) = (3.0*FDOT(K,J)*DT + 2.0*F(K,J))*DT + X0DOT(K,J)
          RDOT = RDOT + (XI(K) - XJ(K))*(VI(K) - VJ(K))
          VREL2 = VREL2 + (VI(K) - VJ(K))**2
   10 CONTINUE
*
*       Only consider approaching bodies.
      IF (RDOT.LE.0.0) GO TO 100
*
*       Predict coordinates & velocities at beginning of step for body #I.
      RDOT0 = 0.0
      DTI = TIME - T0(I) - STEP(I)
      DTJ = TIME - T0(J) - STEP(I)
      DO 20 K = 1,3
          XI(K) = ((FDOT(K,I)*DTI + F(K,I))*DTI + X0DOT(K,I))*DTI +
     &                                                           X0(K,I)
          XJ(K) = ((FDOT(K,J)*DTJ + F(K,J))*DTJ + X0DOT(K,J))*DTJ +
     &                                                           X0(K,J)
          VI(K) = (3.0*FDOT(K,I)*DTI + 2.0*F(K,I))*DTI + X0DOT(K,I)
          VJ(K) = (3.0*FDOT(K,J)*DTJ + 2.0*F(K,J))*DTJ + X0DOT(K,J)
          RDOT0 = RDOT0 + (XI(K) - XJ(K))*(VI(K) - VJ(K))
   20 CONTINUE
*
*       Check pericentre condition (old radial velocity < 0).
      IF (RDOT0.GT.0.0) GO TO 100
*
      ZM = BODY(I) + BODY(J)
      EREL = 0.5*VREL2 - ZM/SQRT(RJMIN2)
*
*       Specify the energy loss (experimental).
      DH = 0.1*VREL2/2.0
      HNEW = EREL - DH
      SEMI0 = -0.5*ZM/EREL
      SEMI = -0.5*ZM/HNEW
      if(rank.eq.0)
     &WRITE (6,50)  NAME(I), NAME(J), SEMI0, SEMI, SQRT(RJMIN2), DH
   50 FORMAT (5X,'NBTIDE:  NAMES A0 A RIJ DH  ',2I5,2F10.5,F8.4,F8.3)
*
*       Skip energy loss treatment unless final orbit is significantly bound.
      IF (SEMI.LT.0.0.OR.SEMI.GT.4.0*RMIN) GO TO 100
*
*       Predict coordinates & velocity for body #J to highest order.
      CALL XVPRED(J,0)
*
*       Set current velocities and form c.m. velocity.
      DO 30 K = 1,3
          VI(K) = XDOT(K,I)
          VJ(K) = XDOT(K,J)
          VCM(K) = (BODY(I)*XDOT(K,I) + BODY(J)*XDOT(K,J))/ZM
   30 CONTINUE
*
*       Introduce velocity change for each component and initialize X0DOT.
      FAC = SQRT((VREL2 - 2.0D0*DH)/VREL2)
      DO 40 K = 1,3
          VR(K) = FAC*(VI(K) - VJ(K))
          XDOT(K,I) = VCM(K) + BODY(J)*VR(K)/ZM
          XDOT(K,J) = VCM(K) - BODY(I)*VR(K)/ZM
          X0DOT(K,I) = VI(K)
          X0DOT(K,J) = VJ(K)
   40 CONTINUE
*
      DE = BODY(I)*BODY(J)*DH/ZM
*       Increase event counter and update total energy loss.
      NDISS = NDISS + 1
      ECOLL = ECOLL + DE
*
*       Set components and phase indicator for new KS regularization.
      ICOMP = MIN(I,J)
      JCOMP = MAX(I,J)
      IPHASE = 1
*
      if(rank.eq.0)
     &WRITE (6,60)  NAME(ICOMP), NAME(JCOMP), SEMI, BE(3) + DE , DE
   60 FORMAT (' TIDAL CAPTURE   NM A E DE ',2I5,F8.4,2F11.6)
*
  100 RETURN
*
      END
