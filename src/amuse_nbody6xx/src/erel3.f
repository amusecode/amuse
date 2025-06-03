      SUBROUTINE EREL3(IM,EB,SEMI)
*
*
*       Dominant two-body energy in three-body system.
*       ----------------------------------------------
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      COMMON/AZREG/  TIME3,TMAX,Q(8),P(8),R1,R2,R3,ENERGY,M(3),X3(3,3),
     &               XDOT3(3,3),CM(10),C11,C12,C19,C20,C24,C25,
     &               NSTEP3,NAME3(3),KZ15,KZ27
      COMMON/AZCOLL/  RK(3),QK(8),PK(8),ICALL,ICOLL,NDISS3
      REAL*8  P1(8)
*
*
*       Specify KS index for least dominant motion (M(3-IM) & M(3)).
      K = 4*(2 - IM)
*
*       Form physical momenta from collision configuration.
      P1(K+1) = QK(K+1)*PK(K+1) - QK(K+2)*PK(K+2) - QK(K+3)*PK(K+3) +
     &                                              QK(K+4)*PK(K+4)
      P1(K+2) = QK(K+2)*PK(K+1) + QK(K+1)*PK(K+2) - QK(K+4)*PK(K+3) -
     &                                              QK(K+3)*PK(K+4)
      P1(K+3) = QK(K+3)*PK(K+1) + QK(K+4)*PK(K+2) + QK(K+1)*PK(K+3) +
     &                                              QK(K+2)*PK(K+4)
*
*       Set distance & physical momentum (eqn (53)).
      RI = QK(K+1)**2 + QK(K+2)**2 + QK(K+3)**2 + QK(K+4)**2
      P1(K+1) = 0.5D0*P1(K+1)/RI
      P1(K+2) = 0.5D0*P1(K+2)/RI
      P1(K+3) = 0.5D0*P1(K+3)/RI
*
*       Form square momentum (c.m. motion of M(IM) & M(3) added below).
      P2 = 0.0D0
      DO 15 J = 1,3
          P2 = P2 + P1(K+J)**2
   15 CONTINUE
*
*       Obtain binding energy from total energy and perturbing function.
      MB = M(IM) + M(3)
      MU2 = M(3-IM)*MB/(M(3-IM) + MB)
      VP = 0.5D0*P2/MU2 - M(1)*M(2)/R3 - M(3-IM)*M(3)/MAX(R1,R2)
      EB = 0.5D0*ENERGY - VP
*
*       Set semi-major axis.
      SEMI = -M(IM)*M(3)/(2.0D0*EB)
*
      RETURN
*
      END
