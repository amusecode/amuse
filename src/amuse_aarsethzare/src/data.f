      SUBROUTINE DATA(BODY, POS, VEL)
*
*
*       Initial conditions.
*       -------------------
*
      IMPLICIT  REAL*8  (A-H,M,O-Z)
      REAL*8 body(3), pos(3,3), vel(3,3)
      COMMON/AZREG/  Q(8),P(8),R,R1,R2,ENERGY,M(3),X(3,3),XDOT(3,3),
     &               RCOLL,ERROR,C11,C12,C19,C20,C24,C25,NSTEPS,NAME(3)
      REAL*8  SUM(7)
*
*
      DO 1 K = 1,7
          SUM(K) = 0.0D0
    1 CONTINUE
*
*       Read initial conditions (or replace by automatic procedure).
c      DO 2 I = 1,3
c          READ (5,*)  M(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3)
c    2 CONTINUE
*
*       Copy initial conditions (or replace by automatic procedure).
      DO 2 I = 1,3
         M(I) = BODY(I)
         DO 9 K = 1,3
             X(K,I) = POS(K,I)
             XDOT(K,I) = VEL(K,I)
    9     CONTINUE
    2 CONTINUE
*
*       Carry on until M(3) < 0.
      IF (M(3).LT.0.0) STOP
*
*       Form total mass & centre of mass terms.
      DO 4 I = 1,3
          SUM(7) = SUM(7) + M(I)
          DO 3 K = 1,3
              SUM(K) = SUM(K) + M(I)*X(K,I)
              SUM(K+3) = SUM(K+3) + M(I)*XDOT(K,I)
    3     CONTINUE
    4 CONTINUE
*
*       Initialize NAME and define X & XDOT with respect to local c.m.
      DO 6 I = 1,3
          NAME(I) = I
          DO 5 K = 1,3
              X(K,I) = X(K,I) - SUM(K)/SUM(7)
              XDOT(K,I) = XDOT(K,I) - SUM(K+3)/SUM(7)
    5     CONTINUE
    6 CONTINUE
*
      RETURN
*
      END
