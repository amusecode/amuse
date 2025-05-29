      SUBROUTINE TRANS4
*
*
*       Transformation to physical momenta & separations.
*       -------------------------------------------------
*
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8  M,MIJ,XNR(9)
      COMMON/CREG/  M(4),X(12),XD(12),P(12),Q(12),TIME4,ENERGY,EPSR2,
     &              XR(9),W(9),R(6),TA(6),MIJ(6),CMX(10),RMAX4,TMAX,
     &              DS,TSTEP,EPS,NSTEP4,NAME4(4),KZ15,KZ27,NREG,NFN
      COMMON/CCOLL/  QK(12),PK(12),ICALL,ICOLL,NDISS4
      COMMON/ICONF/  I1,I2,I3,I4
      COMMON/SAVEP/  PI(12)
*
*
*       Form KS coordinates & chain momenta from configuration QK & PK.
      DO 1 L = 1,3
          L1 = 3*(L - 1) + 1
          L2 = L1 + 1
          L3 = L2 + 1
*
          J1 = 4*(L - 1) + 1
          J2 = J1 + 1
          J3 = J2 + 1
          J4 = J3 + 1
*
          XR(L1) = QK(J1)**2 - QK(J2)**2 - QK(J3)**2 + QK(J4)**2
          XR(L2) = 2.D0*(QK(J1)*QK(J2) - QK(J3)*QK(J4))
          XR(L3) = 2.D0*(QK(J1)*QK(J3) + QK(J2)*QK(J4))
          R(L) = QK(J1)**2 + QK(J2)**2 + QK(J3)**2 + QK(J4)**2
*
          A = 0.5D0/R(L)
          W(L1) = (QK(J1)*PK(J1) - QK(J2)*PK(J2) - QK(J3)*PK(J3) +
     &                                             QK(J4)*PK(J4))*A
          W(L2) = (QK(J2)*PK(J1) + QK(J1)*PK(J2) - QK(J4)*PK(J3) -
     &                                             QK(J3)*PK(J4))*A
          W(L3) = (QK(J3)*PK(J1) + QK(J4)*PK(J2) + QK(J1)*PK(J3) +
     &                                             QK(J2)*PK(J4))*A
    1 CONTINUE
*
      J1 = 3*(I1 - 1)
      J2 = 3*(I2 - 1)
      J3 = 3*(I3 - 1)
      J4 = 3*(I4 - 1)
*
*       Obtain physical momenta of configuration I1, I2, I3, I4.
      DO 10 K = 1,3
          PI(J1+K) = -W(K  )
          PI(J2+K) =  W(K  ) - W(K+3)
          PI(J3+K) =  W(K+3) - W(K+6)
          PI(J4+K) =  W(K+6)
          R(K+3) = 0.0D0
   10 CONTINUE
*
*       Update irregular distances.
      DO 15 K = 1,3
          XNR(K  ) = XR(K  ) + XR(K+3)
          XNR(K+3) = XR(K+3) + XR(K+6)
          XNR(K+6) = XNR(K ) + XR(K+6)
          R(4) = R(4) + XNR(K  )**2
          R(5) = R(5) + XNR(K+3)**2
          R(6) = R(6) + XNR(K+6)**2
   15 CONTINUE
      DO 20 I = 4,6
          R(I) = SQRT(R(I))
   20 CONTINUE
*
      RETURN
*
      END
