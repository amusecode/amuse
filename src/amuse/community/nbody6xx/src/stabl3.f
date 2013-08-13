      SUBROUTINE STABL3(ITERM)
*
*
*       Stability test of three-body system.
*       ------------------------------------
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      COMMON/AZREG/  TIME3,TMAX,Q(8),P(8),R1,R2,R3,ENERGY,M(3),X3(3,3),
     &               XDOT3(3,3),CM(10),C11,C12,C19,C20,C24,C25,
     &               NSTEP3,NAME3(3),KZ15,KZ27
      COMMON/AZOUT/  MB,RB,RI,SEMI,SEMI1,E,E1
*
*
*       Transform to physical variables.
      CALL TRANS3(3)
*
*       Identify index of second binary component and the third body.
      IM = 1
      IF (R2.LT.R1) IM = 2
      I = 3 - IM
*
*       Form scalar product of distance & velocity for body #I.
      RIDOT = X3(1,I)*XDOT3(1,I) + X3(2,I)*XDOT3(2,I) +
     &                             X3(3,I)*XDOT3(3,I)
*
*       Set distance & radial velocity of body #I with respect to binary.
      MB = M(IM) + M(3)
      FAC = CM(7)/MB
      RI = SQRT(X3(1,I)**2 + X3(2,I)**2 + X3(3,I)**2)
      RIDOT = FAC*RIDOT/RI
      RI = FAC*RI
*
*       Evaluate orbital elements.
      RDOT = 0.0D0
      VREL2 = 0.0D0
      VI2 = 0.0D0
      DO 10 K = 1,3
          RDOT = RDOT + 
     &               (X3(K,3) - X3(K,IM))*(XDOT3(K,3) - XDOT3(K,IM))
          VREL2 = VREL2 + (XDOT3(K,3) - XDOT3(K,IM))**2
          VI2 = VI2 + XDOT3(K,I)**2
   10 CONTINUE
*
*       Determine semi-major axis & eccentricity of inner binary.
      RB = MIN(R1,R2)
      SEMI = 2.0D0/RB - VREL2/MB
      SEMI = 1.0/SEMI
      E = SQRT((1.0D0 - RB/SEMI)**2 + RDOT**2/(SEMI*MB))
*
*       Form outer semi-major axis & eccentricity.
      VI2 = VI2*FAC**2
      SEMI1 = 2.0D0/RI - VI2/CM(7)
      SEMI1 = 1.0/SEMI1
      E1 = SQRT((1.0D0 - RI/SEMI1)**2 + (RI*RIDOT)**2/(SEMI1*CM(7)))
*
*       Obtain standard stability ratio (outer pericentre / inner apocentre).
      RATIO = SEMI1*(1.0D0 - E1)/(SEMI*(1.0D0 + E))
*
*       Form coefficients for stability test (Valtonen, Vistas Ast 32, 1988).
*     AM = (2.65 + E)*(1.0 + M(I)/MB)**0.3333
*     FM = (2.0*M(I) - MB)/(3.0*MB)
*
*       Expand natural logarithm for small arguments.
*     IF (ABS(FM).LT.0.67) THEN
*         BM = FM*(1.0 - (0.5 - ONE3*FM)*FM)
*     ELSE
*         BM = LOG(1.0D0 + FM)
*     END IF
*       Define mass dependent criterion of Harrington (A.J. 80) & Bailyn.
*     PCRIT = AM*(1.0 + 0.7*BM)*SEMI
*
*       Form hierarchical stability ratio (Kiseleva & Eggleton 1995).
*     Q0 = MB/M(I)
*     Q1 = MAX(M(3)/M(IM),M(IM)/M(3))
*     Q3 = Q0**0.33333
*     Q13 = Q1**0.33333
*     AR = 1.0 + 3.7/Q3 - 2.2/(1.0 + Q3) + 1.4/Q13*(Q3 - 1.0)/(Q3 + 1.0)
*     PCRIT = AR*SEMI*(1.0D0 + E)
*
*       Adopt the semi-analytical stability criterion (MA 1997).
      Q1 = M(I)/MB
      IF (E1.LT.1.0) THEN
          XFAC = (1.0 + Q1)*(1.0 + E1)/SQRT(1.0 - E1)
      ELSE
          XFAC = 40.0*(1.0 + Q1)
      END IF
      PCRIT = 2.8*XFAC**0.4*SEMI
      PMIN = SEMI1*(1.0D0 - E1)
*
*       Set negative termination index if system is stable and RB > SEMI.
      IF (PCRIT.LT.PMIN.AND.E1.LT.1.0.AND.1./RB.LT.1./SEMI) THEN
          ITERM = -1
*       Obtain Zare's stability parameter (valid for small inclinations).
          M1 = M(IM)
          M2 = M(3)
          M3 = M(I)
          CALL STABLZ(M1,M2,M3,SP)
          if(rank.eq.0)
     &    WRITE (6,20)  SEMI, SEMI1, E, E1, RI, RATIO, SP, PCRIT, PMIN
   20     FORMAT ('  STABT:    A A1 E E1 RI RATIO SP PCR PM ',
     &                           1P,2E10.2,0P,2F6.2,F9.5,2F6.2,1P,2E9.1)
*       Terminate if escaper is outside 3*SEMI (> 0).
      ELSE IF (E1.GT.1.0.AND.1./RI.LT.1./(3.0*SEMI)) THEN
          ITERM = -1
      ELSE
          ITERM = 0
      END IF
*
      RETURN
*
      END
