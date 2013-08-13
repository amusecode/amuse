      REAL*8 FUNCTION STEPI(F,FDOT,F2DOT,F3DOT,ETA)
*
*
*       Fast time-step expression.
*       --------------------------
*
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8  F(3),FDOT(3),F2DOT(3),F3DOT(3)
*
*
*       Set square FDOT & F2DOT and sum of all absolute values.
      FD2 = FDOT(1)**2 + FDOT(2)**2 + FDOT(3)**2
      F2D2 = F2DOT(1)**2 + F2DOT(2)**2 + F2DOT(3)**2
      FI = ABS(F(1)) + ABS(F(2)) + ABS(F(3))
      FD = ABS(FDOT(1)) + ABS(FDOT(2)) + ABS(FDOT(3))
      F2D = ABS(F2DOT(1)) + ABS(F2DOT(2)) + ABS(F2DOT(3))
      F3D = ABS(F3DOT(1)) + ABS(F3DOT(2)) + ABS(F3DOT(3))
*
*       Obtain time-step by simplified relative criterion.           
      STEPI = ETA*(2.0*FI*F2D + FD2)/(2.0*FD*F3D + F2D2)
      STEPI = SQRT(STEPI)
*
      RETURN
*
      END
