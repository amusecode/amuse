      SUBROUTINE ZARE(I1,I2,SP)
*
*
*       Zare stability parameter.
*       -------------------------
*
*       Computes the stability parameter (SP) of a 3-body hierarchical
*       binary system with inner binary (I1,I2) and outer body JCOMP.
*       The system is stable to exchange if SP > 1 is returned.
*       Reference:- K. Zare, Celestial Mechanics 16, 35 (1977).
*       Developed by Murray Alexander 1984.
*
      INCLUDE 'common6.h'
      REAL*8  M1,M2,M3,M(3),AM(3),A(6),XCM(3),VCM(3),X3(3,3),XDOT3(3,3)
      INTEGER  IC(3,2),INAME(3)
      DATA  IC/2,3,1,3,1,2/
      DATA  ERR/1.0D-6/
*
*
*       Copy masses and global indices I1, I2 & JCOMP (set in IMPACT).
      M1 = BODY(I1)
      M2 = BODY(I2)
      M3 = BODY(JCOMP)
      INAME(1) = I1
      INAME(2) = I2
      INAME(3) = JCOMP
*
*       Transform to local centre of mass frame.
      DO 1 K = 1,3
          XCM(K) = 0.0
          VCM(K) = 0.0
          I = INAME(K)
          M(K) = BODY(I)
    1 CONTINUE
      DO 3 L = 1,3
          I = INAME(L)
          DO 2 K = 1,3
              XCM(K) = XCM(K) + M(L)*X(K,I)
              VCM(K) = VCM(K) + M(L)*XDOT(K,I)
    2     CONTINUE
    3 CONTINUE
      DO 5 L = 1,3
          I = INAME(L)
          DO 4 K = 1,3
              X3(K,L) = X(K,I) - XCM(K)/(M1 + M2 + M3)
              XDOT3(K,L) = XDOT(K,I) - VCM(K)/(M1 + M2 + M3)
    4     CONTINUE
    5 CONTINUE
*
*       Evaluate the total angular momentum, kinetic & potential energy.
      ZK3 = 0.0
      DO 10 IA = 1,3
          AM(IA) = 0.0D0
          ZK3 = ZK3 + M(IA)*(XDOT3(1,IA)**2 + XDOT3(2,IA)**2
     &                                      + XDOT3(3,IA)**2)
   10 CONTINUE
      POT3 = 0.0
      DO 30 IA = 1,3
          JA = IC(IA,1)
          KA = IC(IA,2)
          RIJ2 = 0.0
          DO 20 I = 1,3
              AM(IA) = AM(IA) + M(I)*(X3(JA,I)*XDOT3(KA,I) -
     &                                X3(KA,I)*XDOT3(JA,I))
              RIJ2 = RIJ2 + (X3(I,JA) - X3(I,KA))**2
   20     CONTINUE
          POT3 = POT3 - M(JA)*M(KA)/SQRT(RIJ2)
   30 CONTINUE
*
*       Set the bifurcation parameter h*c**2.
      ENERGY = 0.5*ZK3 + POT3
      C1 = ENERGY*(AM(1)**2 + AM(2)**2 + AM(3)**2)
      A(1) = M3 + M1
      A(2) = 3.0D0*M3 + 2.0D0*M1
      A(3) = 3.0D0*M3 + M1
*     A(4) = - A(3)
*       Note bug fix by Douglas Heggie 1/11/2000.
      A(4) = - (3.0D0*M2 + M1)
      A(5) = - (3.0D0*M2 + 2.0D0*M1)
      A(6) = - (M2 + M1)
*
*       Solve quintic equation for the central collinear configuration.
      ICOUNT = 0
      S = 1.0D0
*  50 F1 = A(1)
*     FP1 = F1*5.0D0
*     DO 60 I = 2,5
*         F1 = F1*S + A(I)
*         FP1 = FP1*S + (6-I)*A(I)
*  60 CONTINUE
*
*       Replace by iteration of f(s)/s**2 = 0 for faster convergence (DCH).
   50 F1 = ((A(1)*S + A(2))*S + A(3))*S + A(4) + A(5)/S + A(6)/S**2
      FP1 = (3.0*A(1)*S + 2.0*A(2))*S + A(3) - (2.0*A(6)/S + A(5))/S**2
*
*     F1 = F1*S + A(6)
      DX = -F1/FP1
      IF (ABS(DX).GT.ERR) THEN
          S = S + DX
          ICOUNT = ICOUNT + 1
          IF (ICOUNT.LT.20) GO TO 50
          if(rank.eq.0)
     &    WRITE (6,70)  S
   70     FORMAT (' WARNING!    ZARE    NO CONVERGENCE    S =',F8.4)
      END IF
*
*       Compute critical value CCRIT of C1.
      Y = 1.0D0 + S
      FX = M1*M3 + M2*(M3/Y + M1/S)
      GX = M1*M3 + M2*(M3*Y**2 + M1*S**2)
      CCRIT = -0.5D0*FX**2*GX/(M1 + M2 + M3)
*
*       Form stability parameter SP.
      SP = C1/CCRIT
*
      RETURN
*
      END
