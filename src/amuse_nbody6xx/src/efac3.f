      FUNCTION EFAC3(ZETA,ICASE)
*
*
*       Tidal capture efficiency factor (l = 3).
*       ----------------------------------------
*
*       Fitting functions for the third spherical harmonic index based on
*       Ray et. al. (A & A 184, 164) and Lee & Ostriker (Ap. J. 310, 176).
*       Developed at IOA by S. Portegies Zwart & T. Meinen (A & A 280, 174).
*
      REAL*8  ZETA,COEFF(6),EFAC3
*
*
*       Select coefficients for a given polytropic index (ICASE = 1, 2, 3).
      IF (ICASE.EQ.1) THEN
*       Polytropic index 1.5.
          COEFF(1) = -0.909
          COEFF(2) =  1.574
          COEFF(3) = 12.37
          COEFF(4) =-57.40
          COEFF(5) = 80.10
          COEFF(6) =-46.43
      ELSE
          IF (ICASE.EQ.2) THEN
*       Polytropic index 2.
              COEFF(1) = -1.040
              COEFF(2) = -1.354
              COEFF(3) = 37.64
              COEFF(4) =-139.9
              COEFF(5) = 168.2
              COEFF(6) = -66.53
          ELSE
*       Polytropic index 3.
              COEFF(1) = -1.703
              COEFF(2) = 2.653
              COEFF(3) =-14.34
              COEFF(4) = 12.85
              COEFF(5) = -0.492
              COEFF(6) = -3.600
          END IF
      END IF
*
*       Obtain the tidal energy factor from fourth-order polynomial fit.
      Y = LOG10(ZETA)
      EFAC3 = ((((COEFF(6)*Y + COEFF(5))*Y + COEFF(4))*Y + COEFF(3))*Y
     &                                     + COEFF(2))*Y + COEFF(1)
      EFAC3 = 10.0**EFAC3
*
      RETURN
*
      END
