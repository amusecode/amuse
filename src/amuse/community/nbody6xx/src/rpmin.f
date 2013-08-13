      REAL*8 FUNCTION RPMIN(M1,M2,VSTAR,H,PERI)
*
*       Minimum pericentre distance for energy change (H < 0).
*       ------------------------------------------------------
*
      IMPLICIT REAL*8 (A-H,M,O-Z)
      DATA FUDGE /1.0D-02/
*
*       Adopt simplified eccentricity factor valid for ECC = 1.
      GE = 425.0*3.1415/(32.0*SQRT(2.0D0))
      CLIGHT = 3.0D+05/VSTAR
      FAC = (1.0/CLIGHT)**3.5
*
*       Set mass expression for energy change in N-body units.
      FM = SQRT(M1 + M2)*M1**2*M2**2
*       Form scaled energy (without r_min).
      DE = 8.0/15.0*FAC*FM*CLIGHT**2*GE
*
*       Adopt a fraction of the energy (-MU*H) for determination of RPMIN.
      MU = M1*M2/(M1 + M2)
      DEMAX = -FUDGE*MU*H
      RPMIN = (DE/DEMAX)**(2.0/7.0)
*
*       Include some diagnostics (suppressed).
*     SEMI = -(M1 + M2)/H
*     WRITE (6,5)  SEMI, RPMIN, PERI
*   5 FORMAT (' RPMIN    A RPM PERI ',1P,3E10.2)
*
      RETURN
      END
