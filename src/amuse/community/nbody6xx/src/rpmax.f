      REAL*8 FUNCTION RPMAX(M1,M2,VSTAR,DV,PERI)
*
*
*       Maximum periastron distance for GR capture.
*       -------------------------------------------
*       Quinlan & Shapiro, Ap.J. 343, 725, 1989.
*
*
      IMPLICIT REAL*8 (A-H,M,O-Z)
      DATA FAC /2.6790/
*
*
*       Note M1, M2 and DV are in N-body units.
      MB = M1 + M2
      MFAC = (M1*M2)**(2.0/7.0)*MB**(3.0/7.0)
      CLIGHT = 3.0D+05/VSTAR
      VFAC = CLIGHT**(10.0/7.0)*DV**(4.0/7.0)
      RPMAX = FAC*MFAC/VFAC
*
*       Include diagnostics for testing (suppressed).
*     WRITE (6,5)  M1, M2, DV, RPMAX, PERI
*   5 FORMAT (' RPMAX    M1 M2 DV RPX PERI ',1P,5E10.2)
*     CALL FLUSH(6)
      RETURN
      END
