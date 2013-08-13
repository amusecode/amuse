      SUBROUTINE XTRNLP(Q1,Q3,FP)
*
*
*       External KS perturbation.
*       -------------------------
*
      INCLUDE 'common6.h'
      REAL*8  FP(3)
*
*
*       See whether to include the galactic tidal force.
      IF (TIDAL(1).GT.0.0) THEN
          FP(1) = FP(1) + TIDAL(1)*Q1
          FP(3) = FP(3) + TIDAL(3)*Q3
*       Omit Coriolis terms which do not affect the binding energy.
      END IF
*
*       NB! Second call from routine KSPOLY uses velocity argument.
*
      RETURN
*
      END
