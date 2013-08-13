      SUBROUTINE NEWREG
*
*
*       Initialization of chain regularization.
*       ---------------------------------------
*
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8  M,MIJ,M12,M23,M34,M13,M24,M14
      LOGICAL  SWITCH,GTYPE,GTYPE0
      COMMON/CREG/  M(4),X(12),XD(12),P(12),Q(12),TIME4,ENERGY,EPSR2,
     &              XR(9),W(9),R(6),TA(6),MIJ(6),CM(10),RMAX4,TMAX,
     &              DS,TSTEP,EPS,NSTEP4,NAME4(4),KZ15,KZ27,NREG,NFN
      COMMON/TPR/   SWITCH,GTYPE,GTYPE0
      COMMON/SAVEP/  PI(12)
      COMMON/ICONF/  I1,I2,I3,I4
      EQUIVALENCE  (T11,TA(1)),(T22,TA(2)),(T33,TA(3)),(T12,TA(4)),
     &             (T23,TA(5))
      EQUIVALENCE  (M12,MIJ(1)),(M23,MIJ(2)),(M34,MIJ(3)),
     &             (M13,MIJ(4)),(M24,MIJ(5)),(M14,MIJ(6))
*
*
*       Set physical momenta unless switching case (with PI in COMMON).
      IF (.NOT.SWITCH) THEN
          DO 4 I = 1,4
              DO 2 K = 1,3
                  IK = 3*(I - 1) + K
                  PI(IK) = M(I)*XD(IK)
    2         CONTINUE
    4     CONTINUE
      END IF
*
*       Determine vector chain for regularization.
      CALL STATUS(X,I1,I2,I3,I4)
      IP1 = 3*(I1 - 1)
      IP2 = 3*(I2 - 1)
      IP3 = 3*(I3 - 1)
      IP4 = 3*(I4 - 1)
*
*       Form appropriate momenta & relative coordinates. 
      DO 5 K = 1,3
          W(K  ) = -PI(IP1+K)
          W(K+3) = -PI(IP1+K) - PI(IP2+K)
          W(K+6) = +PI(IP4+K)
          XR(K  ) = X(IP2+K) - X(IP1+K)
          XR(K+3) = X(IP3+K) - X(IP2+K)
          XR(K+6) = X(IP4+K) - X(IP3+K)
    5 CONTINUE
*
*       Perform KS transformations.
      DO 7 L = 1,3
          L1 = 3*(L - 1) + 1
          L2 = L1 + 1
          L3 = L2 + 1
          R2 = XR(L1)**2 + XR(L2)**2 + XR(L3)**2
          R(L) = SQRT(R2)
          LQ1 = 4*(L - 1) + 1
          LQ2 = LQ1 + 1
          LQ3 = LQ2 + 1
          LQ4 = LQ3 + 1
          Q(LQ4) = 0.0D0
          A = R(L) + ABS(XR(L1))
          Q(LQ1) = SQRT(.5D0*A)
          B = 1./(2.0D0*Q(LQ1))
          Q(LQ2) = XR(L2)*B
          Q(LQ3) = XR(L3)*B
          IF (XR(L1).LT.0.0D0) THEN
              U1 = Q(LQ1)
              Q(LQ1) = Q(LQ2)
              Q(LQ2) = U1
              U3 = Q(LQ3)
              Q(LQ3) = Q(LQ4)
              Q(LQ4) = U3
          END IF
          P(LQ1) = 2.D0*(+Q(LQ1)*W(L1) + Q(LQ2)*W(L2) + Q(LQ3)*W(L3))
          P(LQ2) = 2.D0*(-Q(LQ2)*W(L1) + Q(LQ1)*W(L2) + Q(LQ4)*W(L3))
          P(LQ3) = 2.D0*(-Q(LQ3)*W(L1) - Q(LQ4)*W(L2) + Q(LQ1)*W(L3))
          P(LQ4) = 2.D0*(+Q(LQ4)*W(L1) - Q(LQ3)*W(L2) + Q(LQ2)*W(L3))
    7 CONTINUE
*
*       Set mass factors (note scaling by 0.25 for DERQP4).
      T11 = 0.25D0*(.5D0/M(I1) + .5D0/M(I2))
      T22 = 0.25D0*(.5D0/M(I2) + .5D0/M(I3))
      T33 = 0.25D0*(.5D0/M(I3) + .5D0/M(I4))
      T12 = -0.25D0/M(I2)
      T23 = -0.25D0/M(I3)
      M12 = M(I1)*M(I2)
      M23 = M(I2)*M(I3)
      M34 = M(I3)*M(I4)
      M13 = M(I1)*M(I3)
      M24 = M(I2)*M(I4)
      M14 = M(I1)*M(I4)
*
      RETURN
*
      END
