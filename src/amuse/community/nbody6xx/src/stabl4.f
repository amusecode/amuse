      SUBROUTINE STABL4(ITERM)
*
*
*       Stability test of four-body system.
*       -----------------------------------
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      LOGICAL  SWITCH,GTYPE,GTYPE0,AUX
      COMMON/CREG/  M(4),X(12),XD(12),P(12),Q(12),TIME4,ENERGY,EPSR2,
     &              XR(9),W(9),R(6),TA(6),MIJ(6),CM(10),RMAX4,TMAX,
     &              DS,TSTEP,EPS,NSTEP4,NAME4(4),KZ15,KZ27,NREG,NFN
      COMMON/TPR/   SWITCH,GTYPE,GTYPE0
      COMMON/KSAVE/  K1,K2
      REAL*8  RC(3),VC(3),RC0(3),VC0(3)
*
*
*       Transform to physical variables (retain value of SWITCH).
      AUX = SWITCH
      SWITCH = .FALSE.
      CALL ENDREG
      SWITCH = AUX
*
*       Specify indices of two least dominant bodies (denoted K3 & K4).
      K3 = 0
      DO 1 L = 1,4
          IF (L.EQ.K1.OR.L.EQ.K2) GO TO 1
          IF (K3.EQ.0) THEN
              K3 = L
          ELSE
              K4 = L
          END IF
    1 CONTINUE
*
*       Initialize scalars for orbital elements.
      RB0 = 0.0D0
      RB = 0.0D0
      RB1 = 0.0D0
      RB2 = 0.0D0
      RDOT = 0.0D0
      RDOT1 = 0.0D0
      RDOT2 = 0.0D0
      VREL2 = 0.0D0
      VREL20 = 0.0D0
      VREL21 = 0.0D0
      VREL22 = 0.0D0
*
*       Define binary masses of smallest & widest pair (K1 & K2 and K3 & K4).
      MB0 = M(K1) + M(K2)
      MB = M(K3) + M(K4)
*
*       Form separations & velocities of MB0, MB and their relative orbit.
      DO 10 K = 1,3
          J1 = 3*(K1-1) + K
          J2 = 3*(K2-1) + K
          J3 = 3*(K3-1) + K
          J4 = 3*(K4-1) + K
          RC0(K) = (M(K1)*X(J1) + M(K2)*X(J2))/MB0
          RC(K) = (M(K3)*X(J3) + M(K4)*X(J4))/MB
          VC0(K) = (M(K1)*XD(J1) + M(K2)*XD(J2))/MB0
          VC(K) = (M(K3)*XD(J3) + M(K4)*XD(J4))/MB
          RB0 = RB0 + (X(J1) - X(J2))**2
          RB = RB + (X(J3) - X(J4))**2
          RB1 = RB1 + (RC(K) - RC0(K))**2
          RB2 = RB2 + (RC0(K) - X(J3))**2
          RDOT = RDOT + (X(J3) - X(J4))*(XD(J3) - XD(J4))
          RDOT1 = RDOT1 + (RC(K) - RC0(K))*(VC(K) - VC0(K))
          RDOT2 = RDOT2 + (RC0(K) - X(J3))*(VC0(K) - XD(J3))
          VREL2 = VREL2 + (XD(J3) - XD(J4))**2
          VREL20 = VREL20 + (XD(J1) - XD(J2))**2
          VREL21 = VREL21 + (VC(K) - VC0(K))**2
          VREL22 = VREL22 + (VC0(K) - XD(J3))**2
   10 CONTINUE
*
*       Determine semi-major axis of inner binary.
      RB0 = SQRT(RB0)
      SEMI0 = 2.0/RB0 - VREL20/MB0
      SEMI0 = 1.0/SEMI0
      E0 = 1.0 - RB0/SEMI0
*
*       Form semi-major axis & eccentricity of outer pair.
      RB = SQRT(RB)
      SEMI = 2.0D0/RB - VREL2/MB
      SEMI = 1.0/SEMI
*     E = SQRT((1.0D0 - RB/SEMI)**2 + RDOT**2/(SEMI*MB))
*
*       Evaluate orbital elements of relative c.m. motion.
      RB1 = SQRT(RB1)
      SEMI1 = 2.0D0/RB1 - VREL21/CM(7)
      SEMI1 = 1.0/SEMI1
      E1 = SQRT((1.0D0 - RB1/SEMI1)**2 + RDOT1**2/(SEMI1*CM(7)))
*
*       Consider the inner triple.
      MB2 = MB0 + M(K3)
      RB2 = SQRT(RB2)
      SEMI2 = 2.0D0/RB2 - VREL22/MB2
      SEMI2 = 1.0/SEMI2
      E2 = SQRT((1.0D0 - RB2/SEMI2)**2 + RDOT2**2/(SEMI2*MB2))
*
*       Obtain standard stability ratio (outer pericentre / inner apocentre).
      RATIO = SEMI1*(1.0D0 - E1)/(SEMI0*(1.0D0 + E0))
*
*       Form coefficients for stability test (Valtonen, Vistas Ast 32, 1988).
*     AM = (2.65 + E0)*(1.0 + MB0/MB)**0.3333
*     FM = (2.0*MB0 - MB)/(3.0*MB)
*
*       Expand natural logarithm for small arguments.
*     IF (ABS(FM).LT.0.67) THEN
*         BM = FM*(1.0 - (0.5 - 0.3333*FM)*FM)
*     ELSE
*         BM = LOG(1.0D0 + FM)
*     END IF
*
*       Adopt mass dependent criterion of Harrington (A.J. 80) & Bailyn.
*     PCRIT = AM*(1.0 + 0.7*BM)*SEMI0
*
*       Form hierarchical stability ratio (Kiseleva & Eggleton 1995).
*     Q0 = MB/MB0
*     Q1 = MAX(M(K2)/M(K1),M(K1)/M(K2))
*     Q3 = Q0**0.33333
*     Q13 = Q1**0.33333
*     AR = 1.0 + 3.7/Q3 - 2.2/(1.0 + Q3) + 1.4/Q13*(Q3 - 1.0)/(Q3 + 1.0)
*     PCRIT = AR*SEMI0*(1.0D0 + E0)
*
*       Check stability (AM 1997; inner triple or well separated quadruple).
      ITERM = 0
      IF (RB1.GT.5.0*RB2.AND.E2.LT.1.0) THEN
          Q1 = M(K3)/MB0
          XFAC = (1.0 + Q1)*(1.0 + E2)/SQRT(1.0 - E2)
          PCRIT = 2.8*XFAC**0.4*SEMI0
          PMIN = SEMI2*(1.0 - E2)
          IF (PCRIT.LT.PMIN) THEN
              ITERM = -1
              RATIO = SEMI2*(1.0D0 - E2)/(SEMI0*(1.0D0 + E0))
              if(rank.eq.0)
     &        WRITE (6,15)  SEMI0, SEMI2, E0, E2, RATIO, RB0, RB2,
     &                      PCRIT, PMIN
   15         FORMAT ('  STABT:    A0 A2 E0 E2 RATIO R0 R2 PCR PM ',
     &                             1P,2E10.2,0P,2F7.3,F6.2,1P,4E9.1)
          END IF
      ELSE IF (RB1.GT.5.0*MAX(RB0,RB).AND.E1.LT.1.0.AND.
     &         MIN(SEMI0,SEMI).GT.0.0) THEN
*       Choose smallest binary as third body and ignore fudge factor.
          IF (SEMI.GT.SEMI0) THEN
              Q1 = MB0/MB
              AIN = SEMI
          ELSE
              Q1 = MB/MB0
              AIN = SEMI0
          END IF
          XFAC = (1.0 + Q1)*(1.0 + E1)/SQRT(1.0 - E1)
          PCRIT = 2.8*XFAC**0.4*AIN
          PMIN = SEMI1*(1.0 - E1)
          IF (PCRIT.LT.PMIN) THEN
              ITERM = -1
              if(rank.eq.0)
     &        WRITE (6,20)  AIN, SEMI1, E0, E1, RATIO, RB0, RB1,
     &                      PCRIT, PMIN
   20         FORMAT ('  STABQ:    AIN A1 E0 E1 RATIO R0 R1 PCR PM ',
     &                             1P,2E10.2,0P,2F7.3,F6.2,1P,4E9.1)
          END IF
      END IF
*
      RETURN
*
      END
