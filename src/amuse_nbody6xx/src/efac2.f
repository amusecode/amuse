      FUNCTION EFAC2(ZETA,ICASE)
*
*
*       Tidal capture efficiency factor (l = 2).
*       ----------------------------------------
*
*       Fitting functions for the second spherical harmonic index based on
*       Ray et. al. (A & A 184, 164) and Lee & Ostriker (Ap. J. 310, 176).
*       Developed at IOA by S. Portegies Zwart & T. Meinen (A & A 280, 174).
*
      REAL*8  ZETA,COEFF(6),EFAC2
*
*
*       Select coefficients for a given polytropic index (ICASE = 1, 2, 3).
      IF (ICASE.EQ.1) THEN
*       Polytropic index 1.5.
          COEFF(1) =-0.397
          COEFF(2) = 1.678
          COEFF(3) = 1.277
          COEFF(4) =-12.42
          COEFF(5) = 9.446
          COEFF(6) =-5.550
      ELSE
          IF (ICASE.EQ.2) THEN
*       Polytropic index 2.
              COEFF(1) =-0.517
              COEFF(2) =-0.906
              COEFF(3) = 23.88
              COEFF(4) =-93.49
              COEFF(5) = 112.3
              COEFF(6) =-44.15
          ELSE
*       Polytropic index 3.
              COEFF(1) =-1.124
              COEFF(2) = 0.877
              COEFF(3) =-13.37 
              COEFF(4) = 21.55
              COEFF(5) =-16.48
              COEFF(6) = 4.124
          END IF
      END IF
*
*       Obtain the tidal energy factor from fifth-order polynomial fit.
      Y = LOG10(ZETA)
      EFAC2 = ((((COEFF(6)*Y + COEFF(5))*Y + COEFF(4))*Y + COEFF(3))*Y
     &                                     + COEFF(2))*Y + COEFF(1)
      EFAC2 = 10.0**EFAC2
*
      RETURN
*
      END
