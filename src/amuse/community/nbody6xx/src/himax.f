      SUBROUTINE HIMAX(I,IM,ECC,SEMI,EMAX,EMIN,ZI,TG,EDAV)
*
*
*       Maximum eccentricity of hierarchical binary.
*       --------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),YREL(3,MMAX),ZREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      REAL*8  A1(3),A2(3),XREL(3),VREL(3),EI(3),HI(3),HO(3),BHAT(3),
     &        UI(4),V(4)
*
*
*       Copy KS variables to local scalars.
      DO 1 K = 1,4
          UI(K) = UM(K,IM)
          V(K) = UMDOT(K,IM)
    1 CONTINUE
*
*       Transform to physical variables and multiply by 4 (momentum formula).
      CALL KSPHYS(UI,V,XREL,VREL)
      DO 2 K = 1,3
          VREL(K) = 4.0*VREL(K)
    2 CONTINUE
*
*       Specify index for outer body (second KS component) & pair index.
      JCOMP = I + 1
      IPAIR = KVEC(I)
*
*       Ensure that unperturbed binary is resolved.
      IF (LIST(1,2*IPAIR-1).EQ.0.OR.X(1,JCOMP).EQ.X(1,I)) THEN
          CALL RESOLV(IPAIR,1)
      END IF
*
*       Obtain inner & outer angular momentum by cyclic notation.
      A12 = 0.0
      A22 = 0.0
      A1A2 = 0.0
      RI2 = 0.0
      VI2 = 0.0
      RVI = 0.0
      DO 5 K = 1,3
          K1 = K + 1
          IF (K1.GT.3) K1 = 1
          K2 = K1 + 1
          IF (K2.GT.3) K2 = 1
          A1(K) = XREL(K1)*VREL(K2) - XREL(K2)*VREL(K1)
          A2(K) = (X(K1,JCOMP) - X(K1,I))*(XDOT(K2,JCOMP) - XDOT(K2,I))
     &          - (X(K2,JCOMP) - X(K2,I))*(XDOT(K1,JCOMP) - XDOT(K1,I))
          A12 = A12 + A1(K)**2
          A22 = A22 + A2(K)**2
          A1A2 = A1A2 + A1(K)*A2(K)
          RI2 = RI2 + XREL(K)**2
          VI2 = VI2 + VREL(K)**2
          RVI = RVI + XREL(K)*VREL(K)
    5 CONTINUE
*
*       Evaluate orbital parameters for outer orbit from KS elements.
      ZMB = BODY(I) + BODY(JCOMP)
      SEMI1 = -0.5*ZMB/H(IPAIR)
      ECC2 = (1.0 - R(IPAIR)/SEMI1)**2 + TDOT2(IPAIR)**2/(SEMI1*ZMB)
      ECC1 = SQRT(ECC2)
*       Determine inclination in radians.
      FAC = A1A2/SQRT(A12*A22)
      ZI = ACOS(FAC)
*
*       Construct the Runge-Lenz vector (Heggie & Rasio 1995, IAU174, Eq.5).
      EI2 = 0.0
      DO 10 K = 1,3
          EI(K) = (VI2*XREL(K) - RVI*VREL(K))/BODY(I) -
     &                                                 XREL(K)/SQRT(RI2)
          EI2 = EI2 + EI(K)**2
   10 CONTINUE
*
*       Define unit vectors for inner eccentricity and angular momenta.
      COSJ = 0.0
      SJSG = 0.0
      DO 15 K = 1,3
          EI(K) = EI(K)/SQRT(EI2)
          HI(K) = A1(K)/SQRT(A12)
          HO(K) = A2(K)/SQRT(A22)
          COSJ = COSJ + HI(K)*HO(K)
          SJSG = SJSG + EI(K)*HO(K)
   15 CONTINUE
*
*       Form unit vector BHAT and scalars AH & BH (Douglas Heggie, 10/9/96).
      AH = 0.0
      BH = 0.0
      DO 16 K = 1,3
          K1 = K + 1
          IF (K1.GT.3) K1 = 1
          K2 = K1 + 1
          IF (K2.GT.3) K2 = 1
          BHAT(K) = HI(K1)*EI(K2) - HI(K2)*EI(K1)
          AH = AH + EI(K)*HO(K)
          BH = BH + BHAT(K)*HO(K)
   16 CONTINUE
*
*       Evaluate the expressions A & Z.
      A = COSJ*SQRT(1.0 - EI2)
      Z = (1.0 - EI2)*(2.0 - COSJ**2) + 5.0*EI2*SJSG**2
*
*       Obtain maximum inner eccentricity (Douglas Heggie, Sept. 1995).
      Z2 = Z**2 + 25.0 + 16.0*A**4 - 10.0*Z - 20.0*A**2 - 8.0*A**2*Z
      EMAX = ONE6*(Z + 1.0 - 4.0*A**2 + SQRT(Z2))
      EMAX = SQRT(EMAX)
*
*       Form minimum eccentricity (Douglas Heggie, Sept. 1996).
      AZ = A**2 + Z - 2.0
      IF (AZ.GE.0.0) THEN
          AZ1 = 1.0 + Z - 4.0*A**2
          EMIN2 = ONE6*(AZ1 - SQRT(AZ1**2 - 12.0*AZ))
      ELSE
          EMIN2 = 1.0 - 0.5*(A**2 + Z)
      END IF
      EMIN2 = MAX(EMIN2,0.0D0)
      EMIN = SQRT(EMIN2)
*
*       Estimate eccentricity growth time-scale.
      TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
      TK1 = TWOPI*ABS(SEMI1)*SQRT(ABS(SEMI1)/ZMB)
      TG = TK1**2*ZMB*(1.0 - ECC1**2)**1.5/(BODY(JCOMP)*TK)
*
*       Evaluate numerical precession factor (involves elliptic integral).
      CONST = PFAC(A,Z)
      CONST = CONST*4.0/(1.5*TWOPI*SQRT(6.0))
*
*       Convert growth time to units of 10**6 yrs.
      TG = CONST*TG*TSTAR
*
*       Form doubly averaged eccentricity derivative (Douglas Heggie 9/96).
      YFAC = 15.0*BODY(JCOMP)/(4.0*ZMB)*TWOPI*TK/TK1**2
      YFAC = YFAC*ECC*SQRT(1.0 - ECC**2)/(1.0 - ECC1**2)**(1.5)
      EDAV = YFAC*AH*BH
*
      RETURN
*
      END
