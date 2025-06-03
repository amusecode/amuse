      SUBROUTINE XTRNLU(X,M,NB,UG)
*
*
*       Harmonic test potential.
*       ------------------------
*
      IMPLICIT REAL*8 (A-H,M,O-Z)
      REAL*8 X(1),M(1)
*
*
      UG=0.0
      L=0
      DO I=1,NB
      DO K=1,3
      L=L+1
      UG=UG-5.0*X(L)**2*M(I)
      END DO
      END DO
*
      RETURN
      END
