      SUBROUTINE IBLOCK
*
*
*       Initialization of block steps.
*       -----------------------------
*
      INCLUDE 'common6.h'
*
*
*       Form discrete steps in powers of 2.
      DTK(1) = 1.0
      DO 1 K = 2,64
         DTK(K) = 0.5D0*DTK(K-1)
    1 CONTINUE
*
*       Initialize previous and current block time.
      TPREV = 0.0D0
      TBLOCK = 0.0D0
*
      RETURN
*
      END
