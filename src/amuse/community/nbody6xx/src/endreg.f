      SUBROUTINE ENDREG
*
*
*       Transformation to physical variables.
*       -------------------------------------
*
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8  M,MIJ,CM(3)
      LOGICAL  SWITCH,GTYPE,GTYPE0
      COMMON/CREG/  M(4),X(12),XD(12),P(12),Q(12),TIME4,ENERGY,EPSR2,
     &              XR(9),W(9),R(6),TA(6),MIJ(6),CMX(10),RMAX4,TMAX,
     &              DS,TSTEP,EPS,NSTEP4,NAME4(4),KZ15,KZ27,NREG,NFN
      COMMON/TPR/   SWITCH,GTYPE,GTYPE0
      COMMON/SAVEP/  PI(12)
      COMMON/ICONF/  I1,I2,I3,I4
*
*
      DO 1 L = 1,3
          L1 = 3*(L - 1) + 1
          L2 = L1 + 1
          L3 = L2 + 1
*
          K1 = 4*(L - 1) + 1
          K2 = K1 + 1
          K3 = K2 + 1
          K4 = K3 + 1
*
          XR(L1) = Q(K1)**2 - Q(K2)**2 - Q(K3)**2 + Q(K4)**2
          XR(L2) = 2.D0*(Q(K1)*Q(K2) - Q(K3)*Q(K4))
          XR(L3) = 2.D0*(Q(K1)*Q(K3) + Q(K2)*Q(K4))
          R(L) = Q(K1)**2 + Q(K2)**2 + Q(K3)**2 + Q(K4)**2
          A = 0.5D0/R(L)
          W(L1) = (Q(K1)*P(K1) - Q(K2)*P(K2) - Q(K3)*P(K3) +
     &                                         Q(K4)*P(K4))*A
          W(L2) = (Q(K2)*P(K1) + Q(K1)*P(K2) - Q(K4)*P(K3) -
     &                                         Q(K3)*P(K4))*A
          W(L3) = (Q(K3)*P(K1) + Q(K4)*P(K2) + Q(K1)*P(K3) +
     &                                         Q(K2)*P(K4))*A
    1 CONTINUE
*
      IP1 = I1*3 - 3
      IP2 = I2*3 - 3
      IP3 = I3*3 - 3
      IP4 = I4*3 - 3
*
      DO 2 K = 1,3
          PI(IP1+K) = -W(K  )
          PI(IP2+K) =  W(K  ) - W(K+3)
          PI(IP3+K) =  W(K+3) - W(K+6)
          PI(IP4+K) =  W(K+6)
          X(IP1+K) = -XR(K)
          X(IP2+K) = 0.0D0
*       Note terms X(IP1+K) = 0 & X(IP2+K) = XR(K) in early formulation.
          X(IP3+K) = X(IP2+K) + XR(K+3)
          X(IP4+K) = X(IP3+K) + XR(K+6)
    2 CONTINUE
*
*       Skip reduction to centre of mass after switching.
      IF (SWITCH) RETURN
*
*       Obtain velocities from momenta and form c.m. coordinates.
      CM(1) = 0.0D0
      CM(2) = 0.0D0
      CM(3) = 0.0D0
      SMASS = M(1) + M(2) + M(3) + M(4)
      DO 4 I = 1,4
          LI = 3*(I - 1)
          DO 3 K = 1,3
              XD(LI+K) = PI(LI+K)/M(I)
              CM(K) = CM(K) + M(I)*X(LI+K)/SMASS
    3     CONTINUE
    4 CONTINUE
*
*       Express coordinates in the c.m. frame.
      L = 0
      DO 6 I = 1,4
          DO 5 K = 1,3
              L = L + 1
              X(L) = X(L) - CM(K)
    5     CONTINUE
    6 CONTINUE
*
      RETURN
*
      END
