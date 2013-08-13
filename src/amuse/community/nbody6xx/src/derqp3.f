      SUBROUTINE DERQP3(Q,P,QPR,PPR,TPR)
*
*
*       Equations of motion for AZ regularization.
*       ------------------------------------------
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      COMMON/AZREG/  TIME3,TMAX,V(8),W(8),R1,R2,R3,ENERGY,M(3),X3(3,3),
     &               XDOT3(3,3),CM(10),C11,C12,C19,C20,C24,C25,
     &               NSTEP3,NAME3(3),KZ15,KZ27
*       Note COMMON locations of Q & P are replaced by dummy variables.
      COMMON/CLOSE/  RIJ(4,4),RCOLL,QPERI,SIZE(4),ECOLL3,IP(4)
      COMMON/AZCOLL/  RK(3),QK(8),PK(8),ICALL,ICOLL,NDISS3
      COMMON/BSSAVE/  EP(4),DSC,FACM,TFAC,ITFAC,JC
      COMMON/AZCONF/  ICONF(3)
      REAL*8  Q(8),P(8),QPR(8),PPR(8),S2(8),S5(8),S8(8),UPR(8)
*
*
*       Form scalar distances & coefficients.
      R1 = Q(1)*Q(1) + Q(2)*Q(2) + Q(3)*Q(3) + Q(4)*Q(4)
      R2 = Q(5)*Q(5) + Q(6)*Q(6) + Q(7)*Q(7) + Q(8)*Q(8)
      C3 = Q(1)*P(1) - Q(2)*P(2) - Q(3)*P(3) + Q(4)*P(4)
      C4 = Q(5)*P(5) - Q(6)*P(6) - Q(7)*P(7) + Q(8)*P(8)
      C5 = Q(2)*P(1) + Q(1)*P(2) - Q(4)*P(3) - Q(3)*P(4)
      C6 = Q(6)*P(5) + Q(5)*P(6) - Q(8)*P(7) - Q(7)*P(8)
      C7 = Q(3)*P(1) + Q(4)*P(2) + Q(1)*P(3) + Q(2)*P(4)
      C8 = Q(7)*P(5) + Q(8)*P(6) + Q(5)*P(7) + Q(6)*P(8)
      C9 = P(1)*P(1) + P(2)*P(2) + P(3)*P(3) + P(4)*P(4)
      C10 = P(5)*P(5) + P(6)*P(6) + P(7)*P(7) + P(8)*P(8)
      C13 = C11*R2
      C14 = C12*R1
      C15 = C12*C10
      C16 = C11*C9
      C17 = R2*ENERGY
      C18 = R1*ENERGY
*       Note that twice the energy is saved in COMMON.
      C21 = Q(1)*Q(1) - Q(2)*Q(2) - Q(3)*Q(3) + Q(4)*Q(4)
     &    - Q(5)*Q(5) + Q(6)*Q(6) + Q(7)*Q(7) - Q(8)*Q(8)
      C22 = Q(1)*Q(2) - Q(3)*Q(4) - Q(5)*Q(6) + Q(7)*Q(8)
      C23 = Q(1)*Q(3) + Q(2)*Q(4) - Q(5)*Q(7) - Q(6)*Q(8)
      C22 = C22 + C22
      C23 = C23 + C23
      RR = C21*C21 + C22*C22 + C23*C23
      R3 = SQRT(RR)
      A = C25/R3
*
*       Set first derivative of the physical time (standard transformation).
      TPR = R1*R2
      B = A*TPR/RR
*
*       Combine vectorial components with relevant coefficients.
      S2(1) = Q(1)*C4 + Q(2)*C6 + Q(3)*C8
      S2(2) =-Q(2)*C4 + Q(1)*C6 + Q(4)*C8
      S2(3) =-Q(3)*C4 - Q(4)*C6 + Q(1)*C8
      S2(4) = Q(4)*C4 - Q(3)*C6 + Q(2)*C8
      S2(5) = Q(5)*C3 + Q(6)*C5 + Q(7)*C7
      S2(6) =-Q(6)*C3 + Q(5)*C5 + Q(8)*C7
      S2(7) =-Q(7)*C3 - Q(8)*C5 + Q(5)*C7
      S2(8) = Q(8)*C3 - Q(7)*C5 + Q(6)*C7
      S5(1) = P(1)*C4 + P(2)*C6 + P(3)*C8
      S5(2) =-P(2)*C4 + P(1)*C6 + P(4)*C8
      S5(3) =-P(3)*C4 - P(4)*C6 + P(1)*C8
      S5(4) = P(4)*C4 - P(3)*C6 + P(2)*C8
      S5(5) = P(5)*C3 + P(6)*C5 + P(7)*C7
      S5(6) =-P(6)*C3 + P(5)*C5 + P(8)*C7
      S5(7) =-P(7)*C3 - P(8)*C5 + P(5)*C7
      S5(8) = P(8)*C3 - P(7)*C5 + P(6)*C7
      S8(1) = Q(1)*C21 + Q(2)*C22 + Q(3)*C23
      S8(2) =-Q(2)*C21 + Q(1)*C22 + Q(4)*C23
      S8(3) =-Q(3)*C21 - Q(4)*C22 + Q(1)*C23
      S8(4) = Q(4)*C21 - Q(3)*C22 + Q(2)*C23
      S8(5) =-Q(5)*C21 - Q(6)*C22 - Q(7)*C23
      S8(6) = Q(6)*C21 - Q(5)*C22 - Q(8)*C23
      S8(7) = Q(7)*C21 + Q(8)*C22 - Q(5)*C23
      S8(8) =-Q(8)*C21 + Q(7)*C22 - Q(6)*C23
      C1 = C17 - C15 + C19 + A*R2
      C2 = C18 - C16 + C20 + A*R1
*
*       Form derivatives for standard equations of motion.
      DO 10 I = 1,4
          K = I + 4
          QPR(I) = C13*P(I) + C24*S2(I)
          QPR(K) = C14*P(K) + C24*S2(K)
          PPR(I) = C1*Q(I) - C24*S5(I) - B*S8(I)
          PPR(K) = C2*Q(K) - C24*S5(K) - B*S8(K)
   10 CONTINUE
*
*       Check the Bulirsch-Stoer tolerance scaling TFAC = TPR*U (ITFAC > 1).
      IF (ITFAC.GT.0) THEN
          TFAC = FACM*(R1 + R2 + TPR/R3)
          ITFAC = -1
      END IF
*
*       Evaluate the regularized Hamiltonian (denoted GAMMA in 1974 paper).
      GAMMA = 0.5D0*(-R1*C17 + R2*C11*C9 + R1*C12*C10 - C20*R2 - C19*R1
     &                            - A*TPR) + C24*(C3*C4 + C5*C6 + C7*C8)
      G2 = 1.0/SQRT(R1 + R2)
*       Time transformation of 1976 paper DT/DTAU = R1*R2/(R1 + R2)**0.5.
*
*       ------------------------------------------------------------------
*     U = M(1)*M(3)*R2 + M(2)*M(3)*R1 + 0.5D0*A*TPR
*     G2 = 1.0/U
*       Alternative transformation DT/DTAU = R1*R2/U allows explicit time.
*       ------------------------------------------------------------------
*
*       Set the explicit time derivative and check tolerance scaling.
      TPR = TPR*G2
      IF (ITFAC.NE.0) THEN
          TFAC = TFAC*G2
          ITFAC = 0
      END IF
*
*       Include stabilization terms due to modified time transformation.
      GAMMA = GAMMA*G2
*     S6 = C19 + A*R2
      DO 20 I = 1,8
          QPR(I) = G2*QPR(I)
          PPR(I) = G2*(PPR(I) + GAMMA*Q(I))
*         PPR(I) = G2*(PPR(I) + GAMMA*(S6*Q(I) - B*S8(I)))
*         IF (I.EQ.4) S6 = C20 + A*R1
   20 CONTINUE
*
*       Check osculating pericentre of closest pair (first call only).
      IF (ICALL.EQ.0) GO TO 50
*
*       Find index & distance of closest pair and set corresponding KS index.
      IF (R1.LT.R2) THEN
          IM = 1
          RM = R1
      ELSE
          IM = 2
          RM = R2
      END IF    
      K = 1 + 4*(IM - 1)
*
*       Obtain pericentre by Mikkola's algorithm (perturbation < 0.005).
      GI = 2.0*M(3-IM)*(RM/R3)**3/(M(IM) + M(3))
      IF (GI.LT.0.005) THEN
          CALL PERI(Q(K),QPR(K),TPR,M(IM),M(3),QPERI)
      ELSE
          QPERI = RM
      END IF
*
*       Identify the close bodies and check mutual distances.
      I = ICONF(IM)
      J = ICONF(3)
      RIJ(I,J) = MIN(RIJ(I,J),QPERI)
      RIJ(J,I) = MIN(RIJ(J,I),QPERI)
*
*       Ckeck minimum two-body distance and switch off indicator.
      RCOLL = MIN(RCOLL,QPERI)
      ICALL = 0
*
*       Check for tidal two-body interaction or stellar collision.
      IF (QPERI.LT.4.0*MAX(SIZE(IM),SIZE(3))) THEN
*
*       Convert QPR to standard KS with T' = R and form radial velocity R'.
          RPR = 0.0D0
          DO 25 J = K,K+3
              UPR(J) = QPR(J)*RM/TPR
              RPR = RPR + 2.0D0*Q(J)*UPR(J)
   25     CONTINUE
*
*       Save minimum configuration.
          DO 30 J = 1,8
              QK(J) = Q(J)
              PK(J) = P(J)
   30     CONTINUE
*
*       Determine small semi-major axis from non-singular expressions.
          CALL EREL3(IM,EB,SEMI)
*
*       Obtain pericentre time interval from elements & Kepler's equation.
          MB = M(IM) + M(3)
          CALL TPERI(SEMI,Q(K),UPR(K),MB,TP)
*
*       Activate collision indicator & B-S step selector (first time).
          IF (ICOLL.EQ.0) THEN
              ICOLL = -1
              JC = 1
          ELSE
*
*       Check convergence: radial velocity < 1.0E-06 parabolic velocity.
              IF (ABS(RPR).LT.1.0E-08*SQRT(2.0D0*MB*RM).OR.
     &           (ABS(DSC).LT.1.0E-08)) THEN
*       Reset B-S step selector and set collision indicator to AZ pair index.
                  JC = 0
                  ICOLL = IM
              END IF
          END IF
*
*       Set regularized step for DIFSY3 using accelerated iteration rate.
          IF (RM.GT.2.0*QPERI) THEN
              DSC = 2.0*ABS(TP)/TPR
          ELSE
              DSC = ABS(TP)/TPR
          END IF       
*
*       Ensure negative step beyond pericentre (restore step at end).
          IF (JC.GT.0.AND.RPR.GT.0.0D0) THEN
              DSC = -DSC
          END IF
          IF (JC.EQ.0) DSC = 1.0
      END IF
*
   50 RETURN
*
      END
