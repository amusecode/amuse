      SUBROUTINE EXPAND(IPAIR,SEMI0)
*
*
*       Expansion (contraction) of KS orbit.
*       ------------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Evaluate square regularized velocity.
      V20 = 0.0
      DO 10 K = 1,4
          V20 = V20 + UDOT(K,IPAIR)**2
   10 CONTINUE
*
*       Form KS coordinate & velocity scaling factors (general point is OK).
      I = N + IPAIR
      IF (SEMI0.GE.0.0D0) THEN
          SEMI = -0.5D0*BODY(I)/H(IPAIR)
      ELSE
          SEMI = SEMI0
      END IF
      C2 = SQRT(SEMI/SEMI0)
      V2 = 0.5*(BODY(I) + H(IPAIR)*R(IPAIR)*(SEMI/SEMI0))
      C1 = SQRT(V2/V20)
*
*       Re-scale KS variables to new energy (H < 0: constant eccentricity).
      R(IPAIR) = 0.0D0
*       Retain sign of radial velocity for unperturbed KS (apo or peri).
*     TDOT2(IPAIR) = 0.0D0
      DO 20 K = 1,4
          U(K,IPAIR) = C2*U(K,IPAIR)
          UDOT(K,IPAIR) = C1*UDOT(K,IPAIR)
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
*         TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0*U(K,IPAIR)*UDOT(K,IPAIR)
   20 CONTINUE
*
*       Re-initialize KS polynomials for perturbed or hyperbolic motion.
      T0(2*IPAIR-1) = TIME
      IF (LIST(1,2*IPAIR-1).GT.0.OR.H(IPAIR).GT.0.0) THEN
          CALL RESOLV(IPAIR,1)
          CALL KSPOLY(IPAIR,1)
      ELSE
*       Determine new interval for unperturbed motion (>= TK).
          TK = TWOPI*SEMI*SQRT(ABS(SEMI)/BODY(I))
          CALL TPERT(IPAIR,GMIN,DT)
*       Adopt c.m. step instead if integer argument exceeds 10**9.
          IF (DT.LT.2.0E+09*TK) THEN
              K = 1 + INT(0.5D0*DT/TK)
*       Restrict Kepler period to c.m. step (case of very wide orbit).
              STEP(2*IPAIR-1) = MIN(FLOAT(K)*TK,STEP(I))
          ELSE
              STEP(2*IPAIR-1) = STEP(I)
          END IF
      END IF
*
      RETURN
*
      END
