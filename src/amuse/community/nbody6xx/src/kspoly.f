      SUBROUTINE KSPOLY(IPAIR,IMOD)
*
*
*       Initialization of KS polynomials.
*       ---------------------------------
*
      INCLUDE 'common6.h'
      COMMON/SLOW0/  RANGE,ISLOW(10)
      REAL*8  A1(3,4),A(8),FP(4),FPDOT(3),UI(4),UIDOT(4)
      REAL*8  B2(3,4),U2(4),Q2(4)
*
*
*       Specify indices of c.m. & components and perturber membership (+ 1).
      I = N + IPAIR
      I2 = 2*IPAIR
      I1 = I2 - 1
      NNB1 = LIST(1,I1) + 1
      SEMI = -0.5*BODY(I)/H(IPAIR)
*
*       Predict current coordinates & velocities of the perturbers.
      DO 5 L = 2,NNB1
          J = LIST(L,I1)
          CALL XVPRED(J,0)
    5 CONTINUE
*
*       Ensure current values of existing KS solution.
      IF (IPHASE.EQ.0) THEN
          CALL RESOLV(IPAIR,1)
      END IF
*
*       Initialize variables for accumulating contributions.
      DO 10 K = 1,3
          FP(K) = 0.0D0
          FPDOT(K) = 0.0D0
   10 CONTINUE
*       Set dummy index for summation of c.m. or resolved components.
      JDUM = 0
*
*       Obtain the perturbation & first derivative.
      DO 40 L = 2,NNB1
          J = LIST(L,I1)
          IF (J.GT.N) THEN
*       See whether c.m. approximation applies.
              RIJ2 = (X(1,J) - X(1,I))**2 + (X(2,J) - X(2,I))**2 +
     &                                      (X(3,J) - X(3,I))**2
              K = J - N
              IF (RIJ2.LT.CMSEP2*R(K)**2.AND.LIST(1,2*K-1).GT.0) THEN
                  JDUM = 2*K - 1
                  J = JDUM
              END IF
          END IF
*
*       Sum over both components (reverse sign for second component).
   20     II = I1
          DO 30 KCOMP = 1,2
              DO 25 K = 1,3
                  A(K) = X(K,J) - X(K,II)
                  A(K+3) = XDOT(K,J) - XDOT(K,II)
   25         CONTINUE
*       Current velocities are predicted in routines INTGRT, KSMOD, etc.
              RIJ2 = A(1)*A(1) + A(2)*A(2) + A(3)*A(3)
              A8 = BODY(J)/(RIJ2*SQRT(RIJ2))
              A9 = 3.0D0*(A(1)*A(4) + A(2)*A(5) + A(3)*A(6))/RIJ2
              IF (KCOMP.EQ.2) A8 = -A8
*
              DO 28 K = 1,3
                  FP(K) = FP(K) + A(K)*A8
                  FPDOT(K) = FPDOT(K) + (A(K+3) - A(K)*A9)*A8
   28         CONTINUE
              II = I2
   30     CONTINUE
*
          IF (J.EQ.JDUM) THEN
              J = J + 1
              GO TO 20
          END IF
   40 CONTINUE
*
*       Check for external perturbations.
      IF (KZ(14).GT.0.AND.KZ(14).LT.3) THEN
          Q1 = X(1,I1) - X(1,I2)
          Q3 = X(3,I1) - X(3,I2)
          CALL XTRNLP(Q1,Q3,FP)
*
*       Use same formalism for the first derivative (omit Coriolis force).
          VX = XDOT(1,I1) - XDOT(1,I2)
          VZ = XDOT(3,I1) - XDOT(3,I2)
          CALL XTRNLP(VX,VZ,FPDOT)
      END IF
*
*       Transform to regularized force derivative using T' = R.
      DO 45 K = 1,3
          FPDOT(K) = R(IPAIR)*FPDOT(K)
   45 CONTINUE
*
*       Save the relative perturbation.
      FP(4) = SQRT(FP(1)**2 + FP(2)**2 + FP(3)**2)
      GAMMA(IPAIR) = FP(4)*R(IPAIR)**2/BODY(I)
*
*       Copy U & UDOT to scalars and set current transformation matrix.
      DO 48 K = 1,4
          UI(K) = U(K,IPAIR)
          UIDOT(K) = UDOT(K,IPAIR)
   48 CONTINUE
      CALL MATRIX(UI,A1)
*
*       Construct regularized polynomials from explicit derivatives.
      TDOT2(IPAIR) = 0.0D0
      TDOT3(IPAIR) = 0.0D0
      TDOT4 = 0.0D0
      TDOT5 = 0.0D0
      TDOT6 = 0.0D0
      HDOT(IPAIR) = 0.0D0
      HDOT2(IPAIR) = 0.0D0
      HDOT3(IPAIR) = 0.0D0
      HDOT4(IPAIR) = 0.0D0
*
*       Scale perturbing force & first derivative by modification factor.
      IF (IMOD.GT.1) THEN
          ZMOD = FLOAT(ISLOW(IMOD))
          DO 50 K = 1,3
              FP(K) = ZMOD*FP(K)
              FPDOT(K) = ZMOD*FPDOT(K)
   50     CONTINUE
      END IF
*
*       Form regularized force & two derivatives of time & binding energy.
      DO 60 K = 1,4
          A(K) = A1(1,K)*FP(1) + A1(2,K)*FP(2) + A1(3,K)*FP(3)
          A(K+4) = A1(1,K)*FPDOT(1) + A1(2,K)*FPDOT(2) +
     &                                A1(3,K)*FPDOT(3)
          FU(K,IPAIR) = 0.5D0*H(IPAIR)*U(K,IPAIR) + 0.5D0*R(IPAIR)*A(K)
          TDOT2(IPAIR) = TDOT2(IPAIR) + 2.0D0*U(K,IPAIR)*UDOT(K,IPAIR)
          TDOT3(IPAIR) = TDOT3(IPAIR) + 2.0D0*UDOT(K,IPAIR)**2 +
     &                                      2.0D0*U(K,IPAIR)*FU(K,IPAIR)
          HDOT(IPAIR) = HDOT(IPAIR) + 2.0D0*UDOT(K,IPAIR)*A(K)
          HDOT2(IPAIR) = HDOT2(IPAIR) + 2.0D0*FU(K,IPAIR)*A(K)
          FP0(K,IPAIR) = 0.5D0*R(IPAIR)*A(K)
          U2(K) = FU(K,IPAIR)
   60 CONTINUE
*
*       Set regularized velocity matrix (Levi-Civita matrix not required).
      CALL MATRIX(UIDOT,A1)
*
*       Include the whole (L*F)' term in explicit derivatives of FU & H.
      DO 65 K = 1,4
          A(K+4) = A(K+4) + A1(1,K)*FP(1) + A1(2,K)*FP(2) +
     &                                      A1(3,K)*FP(3)
          HDOT2(IPAIR) = HDOT2(IPAIR) + 2.0D0*UDOT(K,IPAIR)*A(K+4)
          FUDOT(K,IPAIR) = 0.5D0*(H(IPAIR)*UDOT(K,IPAIR) + 
     &                            HDOT(IPAIR)*U(K,IPAIR) +
     &                            TDOT2(IPAIR)*A(K) + R(IPAIR)*A(K+4))
          HDOT3(IPAIR) = HDOT3(IPAIR) + 2.0*FUDOT(K,IPAIR)*A(K) +
     &                                            4.0*FU(K,IPAIR)*A(K+4)
          TDOT4 = TDOT4 + 2.0*FUDOT(K,IPAIR)*U(K,IPAIR) +
     &                                     6.0*FU(K,IPAIR)*UDOT(K,IPAIR)
          FD0(K,IPAIR) = 0.5D0*(TDOT2(IPAIR)*A(K) + R(IPAIR)*A(K+4) +
     &                                           HDOT(IPAIR)*U(K,IPAIR))
   65 CONTINUE
*
*       Form L(U'') and add two of the three terms in (L*F)'' (11/1/98).
      CALL MATRIX(U2,B2)
      DO 66 K = 1,4
          Q2(K) = B2(1,K)*FP(1) + B2(2,K)*FP(2) + B2(3,K)*FP(3) +
     &            2.0*(A1(1,K)*FPDOT(1) + A1(2,K)*FPDOT(2) +
     &                                    A1(3,K)*FPDOT(3))
          HDOT3(IPAIR) = HDOT3(IPAIR) + 2.0*UDOT(K,IPAIR)*Q2(K)
   66 CONTINUE
*
*       Form higher derivatives by boot-strapping (Q2 added 11/1/98).
      DO 70 K = 1,4
          FUDOT2(K,IPAIR) = 0.5*(H(IPAIR)*FU(K,IPAIR) +
     &                   HDOT2(IPAIR)*U(K,IPAIR) + TDOT3(IPAIR)*A(K)) +
     &                   TDOT2(IPAIR)*A(K+4) + HDOT(IPAIR)*UDOT(K,IPAIR)
          FUDOT3(K,IPAIR) = 0.5*(H(IPAIR)*FUDOT(K,IPAIR) + TDOT4*A(K) + 
     &                        HDOT3(IPAIR)*U(K,IPAIR)) +
     &                   1.5*(HDOT(IPAIR)*FU(K,IPAIR) +
     &                        HDOT2(IPAIR)*UDOT(K,IPAIR) +
     &                        TDOT3(IPAIR)*A(K+4))
          FUDOT2(K,IPAIR) = FUDOT2(K,IPAIR) + 0.5*R(IPAIR)*Q2(K)
          FUDOT3(K,IPAIR) = FUDOT3(K,IPAIR) + 1.5*TDOT2(IPAIR)*Q2(K)
          HDOT3(IPAIR) = HDOT3(IPAIR) + 2.0*UDOT(K,IPAIR)*Q2(K)
          HDOT4(IPAIR) = HDOT4(IPAIR) + 2.0*FUDOT2(K,IPAIR)*A(K) +
     &                   6.0*(FUDOT(K,IPAIR)*A(K+4) + FU(K,IPAIR)*Q2(K))
          TDOT5 = TDOT5 + 2.0*FUDOT2(K,IPAIR)*U(K,IPAIR) +
     &             8.0*FUDOT(K,IPAIR)*UDOT(K,IPAIR) + 6.0*FU(K,IPAIR)**2
          TDOT6 = TDOT6 + 2.0*FUDOT3(K,IPAIR)*U(K,IPAIR) +
     &                    10.0*FUDOT2(K,IPAIR)*UDOT(K,IPAIR) +
     &                    20.0*FUDOT(K,IPAIR)*FU(K,IPAIR)
   70 CONTINUE
*
*       Check maximum square step (soft binaries & weak hyperbolic pairs).
      IF (ABS(H(IPAIR)).GT.ECLOSE) THEN
          A2 = 0.5/ABS(H(IPAIR))
      ELSE
          A2 = MIN(R(IPAIR)/BODY(I),0.5/ABS(H(IPAIR)))
      END IF
*
*       Assign a conservative value of the initial step (except IMOD > 1).
      FAC = 0.8
      IF (IMOD.GT.1.OR.H(IPAIR).GT.0.0) FAC = 1.0
      DTU = FAC*ETAU*SQRT(A2)/(1.0 + 1000.0*GAMMA(IPAIR))**0.333
*
*       Include convergence criterion DH = H'*DTU + H''*DTU**2/2 = 0.001*|H|.
      IF (IMOD.EQ.1.AND.HDOT2(IPAIR).NE.0.0D0) THEN
          DH = 1.0E-03*MAX(ABS(H(IPAIR)),ECLOSE)
          XF = 2.0*DH/ABS(HDOT2(IPAIR))
          YF = HDOT(IPAIR)/HDOT2(IPAIR)
          DTU1 = SQRT(XF + YF**2) - ABS(YF)
          DTU = MIN(DTU1,DTU)
      END IF
*
*       Initialize reference energy.
      H0(IPAIR) = H(IPAIR)
*
*       Skip Stumpff part and use ZMOD*TK for unperturbed binary.
      IF (LIST(1,I1).EQ.0.AND.SEMI.GT.0.0) THEN
          TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
          STEP(I1) = MIN(TK,STEP(I))
          GO TO 80
      END IF
*
*       Obtain Stumpff coefficients.
   75 Z = -0.5D0*H(IPAIR)*DTU**2
      CALL STUMPF(IPAIR,Z)
      Z5 = SF(6,IPAIR)
      Z6 = SF(7,IPAIR)
      S6 = Z6*TDOT6/720.0
*
*       Convert regularized step to physical units and initialize time T0.
      STEP(I1) = ((((S6 + TDOT5*DTU*Z5/120.0D0 + TDOT4/24.0D0)*DTU +
     &               TDOT3(IPAIR)*ONE6)*DTU + 0.5D0*TDOT2(IPAIR))*DTU +
     &                                                     R(IPAIR))*DTU
*
*       Ensure that c.m. step is not exceeded (H > 0 is OK).
      IF (STEP(I1).GT.STEP(I).AND.H(IPAIR).LT.0.0) THEN
          DTU = 0.5*DTU
          GO TO 75
      END IF
*
*       Initialize step & time and include possible slow-down factor.
   80 DTAU(IPAIR) = DTU
      T0(I1) = TIME
      IF (IMOD.GT.1) THEN
          STEP(I1) = ZMOD*STEP(I1)
      END IF
*
*       Include factorials in force and first derivative.
      DO 85 K = 1,4
          FU(K,IPAIR) = 0.5D0*FU(K,IPAIR)
          FUDOT(K,IPAIR) = ONE6*FUDOT(K,IPAIR)
   85 CONTINUE
*
      RETURN
*
      END
