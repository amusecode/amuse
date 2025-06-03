      SUBROUTINE STABLC(IM,ITERM,SEMI)
*
*
*       Chain stability test.
*       ---------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CLUMP/   BODYS(NMX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
     &                NAMES(NMX,5),ISYS(5)
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIJ(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/INCOND/  X4(3,NMX),XDOT4(3,NMX)
      REAL*8  XCM(3),VCM(3)
*
*
*       Check hierarchical stability of triple using simple criterion.
      IF (N.EQ.3) THEN
          CALL CHSTAB(ITERM)
*
*       Examine general case of N > 3.
      ELSE IF (N.GE.4) THEN
          IESC = 0
          JESC = 0
*       Distinguish between synchronous binary and large distance ratio.
          IF (ISYNC.NE.0) THEN
*       Specify removal of binary or single body (IM = 1, N-1; 1 < IM < N-1).
              IF (IM.EQ.1) THEN
                  IESC = INAME(1)
                  JESC = INAME(2)
              ELSE IF (IM.EQ.N - 1) THEN
                  IESC = INAME(N-1)
                  JESC = INAME(N)
              ELSE
                  IF (RINV(1).LT.RINV(N-1)) THEN
                      IESC = INAME(1)
                  ELSE
                      IESC = INAME(N)
                  END IF
              END IF
          ELSE
*       Check for distant binary or single body.
              IF (IM.EQ.1.OR.IM.EQ.N-1) THEN
                  JM = MIN(IM+1,N-2)
*       Reduce chain if close binary is hierarchical and at beginning or end.
                  IF (1.0/RINV(JM).GT.50.0*SEMI) THEN
                      IESC = INAME(IM)
                      JESC = INAME(IM+1)
                  END IF
              ELSE IF (IM.EQ.2.OR.IM.EQ.3) THEN
                  JM = 1
*       Identify index of the maximum end separation.
                  IF (1.0/RINV(JM).LT.1.0/RINV(N-1)) JM = N - 1
                  IF (1.0/RINV(JM).GT.50.0*SEMI) THEN
                      IESC = INAME(JM)
                  END IF
              END IF   
          END IF
*
*       Copy chain variables to standard form for escape test.
          IF (IESC.GT.0) THEN
              LK = 0
              DO 10 L = 1,N
                  DO 5 K = 1,3
                      LK = LK + 1
                      X4(K,L) = X(LK)
                      XDOT4(K,L) = V(LK)
    5             CONTINUE
   10         CONTINUE
          END IF
*
*       Check binary escape.
          IF (IESC.GT.0.AND.JESC.GT.0) THEN
*       Form coordinates & velocities of escaping binary (local c.m. frame).
              BCM = M(IESC) + M(JESC)
              RI2 = 0.0
              RDOT = 0.0
              DO 20 K = 1,3
                  XCM(K) = (M(IESC)*X4(K,IESC) + M(JESC)*X4(K,JESC))/BCM
                  VCM(K) = (M(IESC)*XDOT4(K,IESC) +
     &                      M(JESC)*XDOT4(K,JESC))/BCM
                  RI2 = RI2 + XCM(K)**2
                  RDOT = RDOT + XCM(K)*VCM(K)
   20         CONTINUE
*
*       Convert to relative distance & radial velocity w.r. to inner part.
              FAC = MASS/(MASS - BCM)
              RI = SQRT(RI2)
              RDOT = FAC*RDOT/RI
              RI = FAC*RI
              RESC = 0.5*RSUM
*
*       Employ parabolic escape criterion (reduce if RI > RSUM/2 & RDOT > 0).
              IF (RI.GT.RESC.AND.RDOT.GT.0.0) THEN
                  IF (RDOT**2.LT.2.0*MASS/RI) THEN
                      IESC = 0
                  END IF
              ELSE
                  IESC = 0
              END IF
*       Check single body escape.
          ELSE IF (IESC.GT.0.AND.JESC.EQ.0) THEN
*       Form relative distance and radial velocity for single particle.
              RI = SQRT(X4(1,IESC)**2 + X4(2,IESC)**2 + X4(3,IESC)**2)
              RDOT = X4(1,IESC)*XDOT4(1,IESC) + X4(2,IESC)*XDOT4(2,IESC)
     &                                        + X4(3,IESC)*XDOT4(3,IESC)
              FAC = MASS/(MASS - M(IESC))
              RDOT = FAC*RDOT/RI
              RI = FAC*RI
              RESC = 0.5*RSUM
*       Ensure that escaper is not close to another body.
              RM = MIN(1.0/RINV(IM),RI)
*
*       Check approximate escape criterion outside RESC (allow RX > 2*RESC).
              IF (RM.GT.RESC.AND.RDOT.GT.0.0) THEN
                  IF (RDOT**2.LT.2.0*MASS/RI) THEN
                      ER = 0.5*RDOT**2 - MASS/RI
                      RX = -MASS/ER
                      IF (ER.LT.0.0.AND.RX.LT.2.0*RESC) THEN
                          IESC = 0
                      END IF
                  END IF
              ELSE
                  IESC = 0
              END IF
          END IF
*
*       Remove identified binary or single escaper.
          IF (IESC.GT.0) THEN
              WRITE (6,40)  IESC, JESC, NAMEC(IESC), RI, RDOT**2,
     &                      2.0*MASS/RI
   40         FORMAT (' STABLC:    IESC JESC NAM RI RD2 2*M/R ',
     &                             2I4,I6,1P,3E9.1)
              ISUB = ISYS(5)
              CALL REDUCE(IESC,JESC,ISUB)
              ITERM = -2
          END IF
      END IF
*
      RETURN
*
      END
