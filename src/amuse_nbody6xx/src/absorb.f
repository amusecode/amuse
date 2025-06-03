      SUBROUTINE ABSORB(ISUB)
*
*
*       Absorption of chain member(s).
*       -----------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX3=3*NMX,NMX4=4*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK,XCM(3),VCM(3)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NCMAX,5),ISYS(5)
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIK(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/INCOND/  X4(3,NMX),XDOT4(3,NMX)
*
*
*       Define discrete time for new polynomial (DT2 < 0 is OK).
      DT2 = T0S(ISUB) + TIMEC - TPREV
      DT8 = (TBLOCK - TPREV)/8.0D0
*       Avoid zero since TBLOCK = TPREV during the first block-step.
      IF (DT8.EQ.0.0D0) DT8 = STEP(ICH)/8.0D0
*
*       Adopt the nearest truncated step (at most 8 subdivisions).
      IF (DT2.GT.0.0D0) THEN
          CALL STEPK(DT2,DTN2)
          DTN = NINT(DTN2/DT8)*DT8
      ELSE
*       Choose negative step if pericentre time < TPREV.
          DT2 = -DT2
          CALL STEPK(DT2,DTN2)
          DTN = -NINT(DTN2/DT8)*DT8
      END IF
*
*       Update time for new polynomial initializations (but check T - T0).
      TIME = TPREV + DTN
*
*       Avoid prediction skip by XVPRED in case TIME - T0 = 0.0.
      IF (TIME - T0(ICH).EQ.0.0D0) TIME = TIME + DT8/16.0D0
*
*       Re-define initial epoch for consistency (ignore phase error).
      T0S(ISUB) = TIME - TIMEC
*
      IF (KZ(30).GT.2) THEN
          WRITE (6,1)  TIME+TOFF, DT2, DT8
    1     FORMAT (' ABSORB:    TIME DT2 DT8 ',F12.6,1P,2E10.2)
      END IF
*
*       Increase membership of chain (JCLOSE: single body or KS pair).
      NCH0 = NCH
      CALL SETSYS
*
*       Improve coordinates & velocities of c.m. body to order F3DOT.
      CALL XVPRED(ICH,-1)
*
      SUM = 0.0
      DO 5 K = 1,3
          XCM(K) = 0.0
          VCM(K) = 0.0
    5 CONTINUE
*
*       Accumulate mass-weighted moments of absorbed particle(s).
      DO 15 L = NCH0+1,NCH
          J = JLIST(L)
          SUM = SUM + BODY(J)
          DO 10 K = 1,3
              XCM(K) = XCM(K) + BODY(J)*X(K,J)
              VCM(K) = VCM(K) + BODY(J)*XDOT(K,J)
   10     CONTINUE
   15 CONTINUE
*
*       Form combined c.m. of old chain and new perturber(s).
      DO 20 K = 1,3
          XCM(K) = (BODY(ICH)*X(K,ICH) + XCM(K))/(BODY(ICH) + SUM)
          VCM(K) = (BODY(ICH)*XDOT(K,ICH) + VCM(K))/(BODY(ICH) + SUM)
   20 CONTINUE
*
*       Define new relative coordinates & velocities and add to chain.
      LK = 3*NCH0
      DO 30 L = NCH0+1,NCH
          J = JLIST(L)
          SIZE(L) = RADIUS(J)
          ISTAR(L) = KSTAR(J)
          RIJ2 = 0.0
          DO 25 K = 1,3
              LK = LK + 1
              X4(K,L) = X(K,J) - XCM(K)
              XDOT4(K,L) = XDOT(K,J) - VCM(K)
              XCH(LK) = X4(K,L)
              VCH(LK) = XDOT4(K,L)
              RIJ2 = RIJ2 + (X4(K,L) - X4(K,L-1))**2
   25     CONTINUE
*       Initialize new inverse distance(s) (some value needed in chpert.f).
          RINV(L-1) = 1.0/SQRT(RIJ2)
   30 CONTINUE
*
*       Re-define old chain variables with respect to new c.m.
      LK = 0
      DO 40 L = 1,NCH0
          DO 35 K = 1,3
              LK = LK + 1
              XCH(LK) = XCH(LK) - (XCM(K) - X(K,ICH))
              VCH(LK) = VCH(LK) - (VCM(K) - XDOT(K,ICH))
   35     CONTINUE
   40 CONTINUE
*
*       Create ghost particle(s) and remove from neighbour lists.
      DO 50 L = NCH0+1,NCH
          J = JLIST(L)
          CALL GHOST(J)
   50 CONTINUE
*
*       Update total mass and initialize new c.m. body variables.
      BODY(ICH) = BODY(ICH) + SUM
      CM(7) = BODY(ICH)
      T0(ICH) = TIME
      DO 60 K = 1,3
          X(K,ICH) = XCM(K)
          X0(K,ICH) = XCM(K)
          XDOT(K,ICH) = VCM(K)
          X0DOT(K,ICH) = VCM(K)
   60 CONTINUE
*
*       Remove ghost particle(s) from neighbour list of #ICH.
      JPERT(1) = ICH
      JLAST = NTOT + 1
      DO 70 L = NCH0+1,NCH
          JLIST(1) = JLIST(L)
          CALL NBREM(JLAST,1,1)
   70 CONTINUE
*
*       Perform re-initialization of c.m. polynomials & perturber list.
      CALL REINIT(ISUB)
*
*       Check optional output for centre of mass condition.
      IF (KZ(30).GT.2) THEN
          DO 80 K = 1,6
              CM(K) = 0.0
   80     CONTINUE
*
          LK = 0
          DO 90 L = 1,NCH
              DO 85 K = 1,3
                  LK = LK + 1
                  CM(K) = CM(K) + BODYC(L)*XCH(LK)
                  CM(K+3) = CM(K+3) + BODYC(L)*VCH(LK)
   85         CONTINUE
   90     CONTINUE
*
          DO 95 K = 1,6
              CM(K) = CM(K)/CM(7)
   95     CONTINUE
*
          WRITE (6,99)  (CM(K),K=1,6)
   99     FORMAT (' ABSORB:   CM ',1P,6E9.1)
      END IF
*
      RETURN
*
      END
