      SUBROUTINE TIDES3(QPERI,M1,M2,VSTAR,HI,ECC,DE)
*
*
*       GR tidal energy loss for interacting stars.
*       -------------------------------------------
*
      IMPLICIT REAL*8 (A-H,M,O-Z)
      REAL*8  DE(2)
*
*
*       Specify constants in N-body units.
      GM = 1.0
      C = 3.0D+05/VSTAR
      FAC = (GM/C**2)**3.5
*
*     FE = 1.0 + 73.0/24.0*ECC2 + 37.0/96.0*ECC**4
*     FE2 = SQRT(ECC**2 - 1.0D0)*(301.0/ ....
*       Adopt simplified eccentricity factor near ECC = 1.
      GE = 425.0*3.1415/(32.0*SQRT(2.0D0))
*
*       Set mass and pericentre expression in N-body units.
      FM = SQRT(M1 + M2)*M1**2*M2**2/QPERI**3.5
*       Form energy change in N-body units (total in DE(1)).
      DE(1) = 8.0/15.0*FAC*FM*C**2*GE
      DE(2) = 0.0
*
*       Obtain specific energy change for diagnostic purposes.
      MU = M1*M2/(M1 + M2)
      DH = DE(1)/MU
*
*       Print diagnostics for capture of hyperbolic orbit.
      IF (rank.eq.0.and.HI.GT.0.0.AND.DH.GT.HI) THEN
          WRITE (6,5)  M1, M2, HI, DH, DE(1), QPERI
    5     FORMAT (' CAPTURE CHECK!    M1 M2 H DH DE QP ',
     &                                1P,6E9.1)
      END IF
*
      RETURN
*
      END
