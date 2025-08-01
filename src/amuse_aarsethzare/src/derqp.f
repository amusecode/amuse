      SUBROUTINE DERQP(Q,P,QPR,PPR,TPR)
*
*
*       Equations of motion.
*       --------------------
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      COMMON/AZREG/  V(8),W(8),R,R1,R2,ENERGY,M(3),X(3,3),XDOT(3,3),
     &               RCOLL,ERROR,C11,C12,C19,C20,C24,C25,NSTEPS,NAME(3)
*       Note COMMON locations of Q & P are replaced by dummy variables.
      COMMON/CLOSE/  RIJ(3,3),ICALL
      COMMON/BSSAVE/  EP(4),TFAC,ITFAC,JC,NHALF2
      COMMON/NCALL/  NFN
      REAL*8  Q(8),P(8),QPR(8),PPR(8),S2(8),S5(8),S8(8)
*
*
      NFN = NFN + 1
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
      R = SQRT(RR)
      A = C25/R
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
*       Check tolerance scaling TFAC = TPR*U (ITFAC > 1).
      IF (ITFAC.GT.0) THEN
          TFAC = M(3)*(M(1)*R2 + M(2)*R1)
          ITFAC = -1
      END IF
*
*       -------------------------------------
      IF (K.GT.0) RETURN
*       Dummy test for skipping stabilization.
*       --------------------------------------
*
*       Adopt stabilization with DT/DTAU = 1/U (Cortina paper 1976).
      GAMMA = 0.5D0*(-R1*C17 + R2*C11*C9 + R1*C12*C10 - C20*R2 - C19*R1
     &                            - A*TPR) + C24*(C3*C4 + C5*C6 + C7*C8)
*     U = M(1)*M(3)*R2 + M(2)*M(3)*R1 + 0.5D0*A*TPR
*     G2 = 1.0/U
      G2 = 1.0/SQRT(R1 + R2)
      G2 = 1.0
      TPR = TPR*G2
      IF (ITFAC.NE.0) THEN
          TFAC = TFAC*G2
          ITFAC = 0
      END IF
*
*       Include stabilization terms due to modified time transformation.
*     GAMMA = G2*GAMMA
*     GAMMA = 0.0
*     S6 = C19 + A*R2
      DO 20 I = 1,8
          QPR(I) = G2*QPR(I)
          PPR(I) = G2*(PPR(I) + GAMMA*Q(I))
*         PPR(I) = G2*(PPR(I) + GAMMA*(S6*Q(I) - B*S8(I)))
*         IF (I.EQ.4) S6 = C20 + A*R1
   20 CONTINUE
*
*       Check minimum two-body separation (first call only!).
      IF (ICALL.EQ.0) RETURN
*
*       Set minimum configuration index & separation.
      IF (R1.LT.R2) THEN
          IM = 1
          RI = R1
      ELSE
          IM = 2
          RI = R2
      END IF
*
*       Form the explicit time derivative for routine PERI.
*     TPR = TPR*G2
*
*       Determine the osculating pericentre by Mikkola's algorithm.
      K = 1 + 4*(IM - 1)
      CALL PERI(Q(K),QPR(K),TPR,M(IM),M(3),QPERI)
*
*       Define relative perturbation on the close pair.
      GI = 2.0*M(3-IM)*(RI/R)**3/(M(IM) + M(3))
*
*       Identify the close bodies and check mutual distances.
      IF (GI.LT.0.005) THEN
          I = NAME(IM)
          J = NAME(3)
          RIJ(I,J) = MIN(RIJ(I,J),QPERI)
          RIJ(J,I) = MIN(RIJ(J,I),QPERI)
      END IF
*
*       Ckeck minimum two-body distance and switch off indicator.
      IF (GI.LT.0.005) RCOLL = MIN(RCOLL,QPERI)
      ICALL = 0
*
      RETURN
*
      END
