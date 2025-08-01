      SUBROUTINE PERI(U,UPR,TPR,M1,M2,Q)
*
*
*       Pericentre determination.
*       -------------------------
*
*       Routine for estimating the closest two-body approach based on the
*       KS variable U and its derivative UPR. The derivative of time, TPR,
*       is also used to convert to physical variables. This routine can be
*       used in any two-body KS, AZ and chain regularization which provides
*       a KS transformation U of the relative vector and its derivative
*       with respect to an arbitrary time variable.
*       NB!! Only reliable for the first function call with Bulirsch-Stoer.
*       Developed by Seppo Mikkola, May 1990.
*
*       U = KS coordinates of the relative position vector of bodies M1, M2.
*       UPR = derivatives of U used by equation of motion.
*       TPR = derivative of physical time used by equation of motion.
*       M1,M2 = masses of the bodies.
*       Q = osculating pericentre distance.
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      REAL*8  U(4),UPR(4)
*
*
*       Form the physical distance and square KS velocity.
      R = U(1)**2 + U(2)**2 + U(3)**2 + U(4)**2
      RD = UPR(1)**2 + UPR(2)**2 + UPR(3)**2 + UPR(4)**2
*
*       Evaluate the square of the KS cross-product.
      CROSS = 0.0D0
      DO 2 I = 1,3
          DO 1 J = I+1,4
              CROSS = CROSS + (U(I)*UPR(J) - U(J)*UPR(I))**2
    1     CONTINUE
    2 CONTINUE
*
*       Transform to physical velocity and include the mass factor.
      COEF = (2.0D0*R/TPR)**2/(M1 + M2)
*
*       Form semi-latus rectum (note that P contains factor(s) of R).
      P = CROSS*COEF
*
*       Set inverse semi-major axis (division by R is safe because of P).
      OA = (2.0D0 - COEF*RD)/R
*
*       Obtain eccentricity & pericentre distance (avoid argument < 0).
      ECC2 = MAX(1.0D0 - P*OA,0.0D0)
      ECC = SQRT(ECC2)
      Q = P/(1.0D0 + ECC)
*
      RETURN
*
      END
