      SUBROUTINE EREL(IM,EB,SEMI)
*
*
*       Dominant two-body energy in chain system.
*       -----------------------------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      COMMON/KSAVE/  K1,K2
*
*
*       Obtain physical momenta & separations of pericentre configuration.
      CALL TRANSK
*
*       Form kinetic energy terms of dominant c.m. (K1 + K2) and the rest.
      ZK = 0.0D0
      P0 = 0.0D0
      DO 10 K = 1,3
          PS = 0.0
*       Exclude dominant bodies using names of chain members set in DERQP.
          DO 5 I = 1,N
              IF (INAME(I).NE.K1.AND.INAME(I).NE.K2) THEN
                  J = 3*(I - 1)
                  PS = PS + PI(J+K)
                  ZK = ZK + PI(J+K)**2/MC(I)
              END IF
    5     CONTINUE
          P0 = P0 + PS**2
   10 CONTINUE
*
*       Evaluate potential energy due to non-dominant chain distances.
      POT = 0.0D0
      DO 20 I = 1,N-1
          IF (I.EQ.IM) GO TO 20
          L = I + 1
          POT = POT + MIJ(I,L)*RINV(I)
   20 CONTINUE
*
*       Add non-chained contributions to the potential energy.
      IJ = N - 1
      DO 30 I = 1,N-2
          DO 25 J = I+2,N
              IJ = IJ + 1
              POT = POT + MIJ(I,J)*RINV(IJ)
   25     CONTINUE
   30 CONTINUE
*
*       Obtain binding energy from total energy and perturbing function.
      MB = M(K1) + M(K2)
      VP = 0.5D0*(P0/MB + ZK) - POT
      EB = ENERGY - VP
*
*       Set semi-major axis.
      SEMI = -M(K1)*M(K2)/(2.0D0*EB)
*
      RETURN
*
      END
