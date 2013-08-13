      SUBROUTINE EREL4(IM,EB,SEMI)
*
*
*       Dominant two-body energy in four-body system.
*       ---------------------------------------------
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      COMMON/CREG/  M(4),X(12),XD(12),P(12),Q(12),TIME4,ENERGY,EPSR2,
     &              XR(9),W(9),R(6),TA(6),MIJ(6),CM(10),RMAX4,TMAX,
     &              DS,TSTEP,EPS,NSTEP4,NAME4(4),KZ15,KZ27,NREG,NFN
      COMMON/SAVEP/  PI(12)
      COMMON/KSAVE/  K1,K2
*
*
*       Obtain physical momenta of pericentre configuration.
      CALL TRANS4
*
*       Form kinetic energy terms of dominant c.m. (K1 + K2) and the rest.
      ZK = 0.0D0
      P0 = 0.0D0
      DO 10 K = 1,3
          PS = 0.0
          DO 8 I = 1,4
*       Note PI is associated with the bodies (PI sequential in CHAIN code).
              IF (I.NE.K1.AND.I.NE.K2) THEN
                  J = 3*(I - 1)
                  PS = PS + PI(J+K)
                  ZK = ZK + PI(J+K)**2/M(I)
              END IF
    8     CONTINUE
          P0 = P0 + PS**2
   10 CONTINUE
*
*       Evaluate potential energy due to non-singular terms.
      POT = 0.0D0
      DO 20 I = 1,6
          IF (I.EQ.IM) GO TO 20
          POT = POT + MIJ(I)/R(I)
   20 CONTINUE
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
