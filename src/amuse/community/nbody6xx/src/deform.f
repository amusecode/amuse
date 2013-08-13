      SUBROUTINE DEFORM(IPAIR,ECC0,ECC)
*
*
*       Deformation of elliptic orbit.
*       ------------------------------
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
*       Form KS coordinate & velocity scaling factors at peri or apocentre.
      I = N + IPAIR
      SEMI = -0.5*BODY(I)/H(IPAIR)
      IF (R(IPAIR).LT.SEMI) THEN
          EFAC = (1.0 - ECC)/(1.0 - ECC0)
          RNEW = SEMI*(1.0 - ECC)
      ELSE
          EFAC = (1.0 + ECC)/(1.0 + ECC0)
          RNEW = SEMI*(1.0 + ECC)
          RNEW = EFAC*R(IPAIR)
      END IF
      C1 = SQRT(EFAC)
      V2 = 0.5*(BODY(I) + H(IPAIR)*RNEW)
      IF(V2.LE.0.D0)THEN
         C2 = 1.0D-06
      ELSE
         C2 = SQRT(V2/V20)
      ENDIF
*
*       Re-scale KS variables at constant energy to new eccentricity (ECC).
      R(IPAIR) = 0.0D0
*       Retain sign of radial velocity for unperturbed KS (apo or peri).
      TDOT2(IPAIR) = 0.0D0
      DO 20 K = 1,4
          U(K,IPAIR) = C1*U(K,IPAIR)
          UDOT(K,IPAIR) = C2*UDOT(K,IPAIR)
          U0(K,IPAIR) = U(K,IPAIR)
          R(IPAIR) = R(IPAIR) + U(K,IPAIR)**2
          TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0*U(K,IPAIR)*UDOT(K,IPAIR)
   20 CONTINUE
*
*       Re-initialize KS polynomials for perturbed motion.
      T0(2*IPAIR-1) = TIME
      IF (LIST(1,2*IPAIR-1).GT.0) THEN
          CALL RESOLV(IPAIR,1)
          CALL KSPOLY(IPAIR,1)
      ELSE
*       Determine new interval for unperturbed motion (>= TK).
          TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
          CALL TPERT(IPAIR,GMIN,DT)
*       Adopt c.m. step instead if integer argument exceeds 10**9.
          IF (DT.LT.2.0E+09*TK) THEN
              K = 1 + INT(0.5D0*DT/TK)
*       Restrict Kepler period to c.m. step (case of very wide orbit).
              STEP(2*IPAIR-1) = FLOAT(K)*MIN(TK,STEP(I))
          ELSE
              STEP(2*IPAIR-1) = STEP(I)
          END IF
      END IF
*
      RETURN
*
      END
