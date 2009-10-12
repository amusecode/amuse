!     Helper function: calculate degeneracy parameter xi given lnf
      double precision function calculate_degeneracy_parameter(lnf)
      implicit none
      double precision, intent(in) :: lnf
      double precision :: uf, wf, xi

!      if (auf < -100.0) then
!         xi = auf
!      else if (auf > 100.0) then
!         uf = dexp(auf)
!         xi = 2.0d0*sqrt(uf)
!      else
         uf = dexp(lnf)
         wf = dsqrt(1.0d0 + uf)
         xi = dlog(uf) + 2.0d0*(wf - dlog(1.0d0 + wf))
!      end if

      calculate_degeneracy_parameter = xi
      end function

!     Find the parameter "log f" from the positron degeneracy parameter XI
      FUNCTION FIND_POSITRON_F(XI)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: XI
      DOUBLE PRECISION :: FIND_POSITRON_F
      DOUBLE PRECISION, EXTERNAL :: CALCULATE_DEGENERACY_PARAMETER
      DOUBLE PRECISION, EXTERNAL :: SOLVE_F_INTERVAL

!     Asymptotic expansions
      IF (XI > 15.0) THEN
         FIND_POSITRON_F = log(0.25*XI**2 + 1)
         return
      ELSE IF (XI < -5) THEN
         FIND_POSITRON_F = XI - 2.0 + log(4.0)
         return
      END IF
!     If logf is below -100.0 we really don't care enough to determine it
!     more accurately; all Fermi functions go to 0 as exp(xi) in this region
!     anyway.
      FIND_POSITRON_F = SOLVE_F_INTERVAL(CALCULATE_DEGENERACY_PARAMETER,
     &                     XI,   -5.d0, 15.0) 
      END FUNCTION

!     Calculate Fermi-Dirac integrals for positrons, passing in the appropriate
!     values of F and G for Positrons
      SUBROUTINE POSITRON_FDIRAC ( F, G )
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT (IN) :: F, G
!     Because the original FDIRAC uses a COMMON block and we're essentially
!     the same function, we do the same...
      DOUBLE PRECISION :: DE(16), DP(16), DE_STORE(16)
      COMMON /STATFD/ DE
      COMMON /STATFP/ DP

!     If the degeneracy parameter for positrons is very small, we don't care
!     because all contributions are going to be essentially 0 anyway.
      IF (F < 1.0D-300) THEN
         DP(:) = 0.0D0
         RETURN
      END IF

!     Store electron data so that we can restore it after the FD functions
!     have been calculated
      DE_STORE(:) = DE(:)

!     Compute positron FD functions
      CALL FDIRAC(F, G)

!     Set positron data, then restore backed-up electron data
      DP(:) = DE(:)
      DE(:) = DE_STORE(:)
      END SUBROUTINE


!     Solve function: find the point x in the interval (x1, x2) where the
!     function value equals y, using bisection.
!     Returns x, or x1 - 1 if the function did not converge.
!     Code based on the Numerical Recipes routine RTSEC
      FUNCTION SOLVE_F_INTERVAL(FUNC,Y,X1,X2)
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-8
      INTEGER, PARAMETER :: MAXITER = 10
      DOUBLE PRECISION, INTENT(IN) :: Y, X1, X2
      DOUBLE PRECISION, EXTERNAL :: FUNC
      DOUBLE PRECISION :: SOLVE_F_INTERVAL
      DOUBLE PRECISION :: DX, FH, FL, XL, XH
      INTEGER :: I

      FL = FUNC(X1) - Y
      FH = FUNC(X2) - Y
      IF (ABS(FL) < ABS(FH)) THEN
!        Swap highest/lowest variables
         XH = X1
         XL = X2
!        Swap function values
         DX = FL
         FL = FH
         FH = DX
      ELSE
         XL = X1
         XH = X2
      ENDIF
!     Loop until converged to the required accuracy
      DO I = 1, MAXITER
!        Calculate correction
         DX = (XL - XH) * FH / (FH - FL)
!        On each iteration, flip boundaries
         XL = XH
         FL = FH
         XH = XH + DX
         FH = FUNC(XH) - Y
!        Converged?
         IF (DABS(DX) < EPS .OR. FH == 0.0D0) EXIT
      END DO
      SOLVE_F_INTERVAL = XH
      END FUNCTION

