      SUBROUTINE STABLZ(M1,M2,M3,SP)
*
*
*       Zare stability parameter.
*       -------------------------
*
*       Computes the stability parameter (SP) of a 3-body hierarchical
*       binary system with (M1,M2) constituting the inner binary.
*       The system is stable if SP > 1 is returned.
*       Reference:- K. Zare, Celestial Mechanics 16, 35 (1977).
*       Developed by Murray Alexander 1984.
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      COMMON/AZREG/  TIME3,TMAX,Q(8),P(8),R1,R2,R3,ENERGY,M(3),X3(3,3),
     &               XDOT3(3,3),CM(10),C11,C12,C19,C20,C24,C25,
     &               NSTEP3,NAME3(3),KZ15,KZ27
      REAL*8  M1,M2,M3,AM(3),A(6)
      INTEGER  IC(3,2)
      DATA  IC/2,3,1,3,1,2/
      DATA  ERR/1.0D-6/
*
*
*       Evaluate the total angular momentum.
      DO 10 IA = 1,3
          AM(IA) = 0.0D0
   10 CONTINUE
      DO 30 IA = 1,3
          JA = IC(IA,1)
          KA = IC(IA,2)
          DO 20 I = 1,3
              AM(IA) = AM(IA) + M(I)*(X3(JA,I)*XDOT3(KA,I) -
     &                                X3(KA,I)*XDOT3(JA,I))
   20     CONTINUE
   30 CONTINUE
*
*       Set the bifurcation parameter h*c**2 (note ENERGY = 2*h).
      C1 = 0.5D0*ENERGY*(AM(1)**2 + AM(2)**2 + AM(3)**2)
      A(1) = M3 + M1
      A(2) = 3.0D0*M3 + 2.0D0*M1
      A(3) = 3.0D0*M3 + M1
      A(4) = - A(3)
      A(5) = - (3.0D0*M2 + 2.0D0*M1)
      A(6) = - (M2 + M1)
*
*       Solve quintic equation for the central collinear configuration.
      ICOUNT = 0
      X = 1.0D0
   50 F1 = A(1)
      FP1 = F1*5.0D0
      DO 60 I = 2,5
          F1 = F1*X + A(I)
          FP1 = FP1*X + (6-I)*A(I)
   60 CONTINUE
*
      F1 = F1*X + A(6)
      DX = -F1/FP1
      IF (ABS(DX).GT.ERR) THEN
          X = X + DX
          ICOUNT = ICOUNT + 1
          IF (ICOUNT.LT.20) GO TO 50
          if(rank.eq.0)
     &    WRITE (6,70)  X
   70     FORMAT (5X,'WARNING!   NO CONVERGENCE IN STABLZ   X =',F8.4)
      END IF
*
*       Compute critical value CCRIT of C1.
      Y = 1.0D0 + X
      FX = M1*M3 + M2*(M3/Y + M1/X)
      GX = M1*M3 + M2*(M3*Y**2 + M1*X**2)
      CCRIT = -0.5D0*FX**2*GX/(M1 + M2 + M3)
*
*       Form stability parameter SP.
      SP = C1/CCRIT
*
      RETURN
*
      END
