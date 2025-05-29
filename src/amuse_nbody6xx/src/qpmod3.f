      SUBROUTINE QPMOD3(IM,ITERM)
*
*
*       Modification of AZ variables for tidal dissipation.
*       ---------------------------------------------------
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      COMMON/AZREG/  TIME3,TMAX,Q(8),P(8),R1,R2,R3,ENERGY,M(3),X3(3,3),
     &               XDOT3(3,3),CM(10),C11,C12,C19,C20,C24,C25,
     &               NSTEP3,NAME3(3),KZ15,KZ27
      COMMON/CLOSE/  RIJ(4,4),RCOLL,QPERI,SIZE(4),ECOLL3,IP(4)
      COMMON/AZCOLL/  RK(3),QK(8),PK(8),ICALL,ICOLL,NDISS3
      REAL*8  DE(2),RADIUS(2),P1(8)
      INTEGER IS(2)
      INCLUDE "mpif.h"
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      COMMON/MPIDAT/group,rank,ierr,isize,status
*
*
*       Obtain tidal energy loss (bodies #IM & 3 with separation QPERI).
      RADIUS(1) = SIZE(IM)
      RADIUS(2) = SIZE(3)
      IS(1) = IP(IM)
      IS(2) = IP(3)
*
      CALL TIDES(QPERI,M(IM),M(3),RADIUS(1),RADIUS(2),IS,DE)
*
*       Determine physical momenta of M(IM) & M(3-IM).
      DO 10 KCOMP = 1,2
          K = 4*(KCOMP - 1)
*
*       Form product of half Levi-Civita matrix & regularized momentum.
          P1(K+1) = Q(K+1)*P(K+1) - Q(K+2)*P(K+2) - Q(K+3)*P(K+3) +
     &                                              Q(K+4)*P(K+4)
          P1(K+2) = Q(K+2)*P(K+1) + Q(K+1)*P(K+2) - Q(K+4)*P(K+3) -
     &                                              Q(K+3)*P(K+4)
          P1(K+3) = Q(K+3)*P(K+1) + Q(K+4)*P(K+2) + Q(K+1)*P(K+3) +
     &                                              Q(K+2)*P(K+4)
          P1(K+4) = 0.0D0
*
*       Divide by the distance to obtain physical momentum (eqn (53)).
          RI = Q(K+1)**2 + Q(K+2)**2 + Q(K+3)**2 + Q(K+4)**2
          P1(K+1) = 0.5D0*P1(K+1)/RI
          P1(K+2) = 0.5D0*P1(K+2)/RI
          P1(K+3) = 0.5D0*P1(K+3)/RI
   10 CONTINUE
*
*       Define KS index for closest bodies and second interaction.
      K1 = 4*(IM - 1)
      K2 = 4*(2 - IM)
*
*       Evaluate kinetic energy term of least dominant interaction.
      P2 = 0.0D0
      DO 15 K = 1,3
          P2 = P2 + P1(K2+K)**2
*       Form relative momentum of dominant motion from c.m. condition.
          P1(K1+K) = 2.0D0*P1(K1+K) + P1(K2+K)
   15 CONTINUE
*
*       Set consistent third distance for energy relation.
      A21 = Q(1)*Q(1) - Q(2)*Q(2) - Q(3)*Q(3) + Q(4)*Q(4)
     &    - Q(5)*Q(5) + Q(6)*Q(6) + Q(7)*Q(7) - Q(8)*Q(8)
      A22 = Q(1)*Q(2) - Q(3)*Q(4) - Q(5)*Q(6) + Q(7)*Q(8)
      A23 = Q(1)*Q(3) + Q(2)*Q(4) - Q(5)*Q(7) - Q(6)*Q(8)
      A22 = A22 + A22
      A23 = A23 + A23
      R3 = SQRT(A21*A21 + A22*A22 + A23*A23)
*
*       Define radius, mass & reduced mass for the dominant bodies.
      RM = MIN(R1,R2)
      MB = M(IM) + M(3)
      MU = M(IM)*M(3)/MB
*
*       Obtain binding energy from total energy and perturbing function.
      MU2 = M(3-IM)*MB/(M(3-IM) + MB)
      VP = 0.5D0*P2/MU2 - M(1)*M(2)/R3 - M(3-IM)*M(3)/MAX(R1,R2)
      H = (0.5D0*ENERGY - VP)/MU
*
*       Set semi-major axis, eccentricity & pericentre (assume R' = 0).
      SEMI = -0.5D0*MB/H
      ECC = 1.0 - RM/SEMI
      PERI = SEMI*(1.0D0 - ECC)
*
*       Determine new eccentricity from angular momentum conservation.
      DH = -(DE(1) + DE(2))/MU
      AM0 = SEMI*(1.0D0 - ECC**2)
      ECC2 = ECC**2 + 2.0D0*AM0*DH/MB
      IF (ECC2.GT.0.0D0) THEN
          ECC1 = SQRT(ECC2)
      ELSE
          ECC1 = 0.0
      END IF
*
*       Update binding energy and set new semi-major axis & pericentre.
      H1 = H + DH
      SEMI1 = -0.5D0*MB/H1
      PERI1 = SEMI1*(1.0D0 - ECC1)
*
*       Form KS coordinate scaling factor from pericentre ratio.
      C1 = SQRT(PERI1/PERI)
*       Specify KS velocity scaling from angular momentum conservation.
      C2 = 1.0/C1
*
*       Modify the KS coordinates & physical momentum of dominant motion.
      RI = 0.0D0
      DO 20 K = 1,4
          Q(K1+K) = C1*Q(K1+K)
          RI = RI + Q(K1+K)**2
          P1(K1+K) = C2**2*P1(K1+K)
          P1(K1+K) = 0.5D0*(P1(K1+K) - P1(K2+K))
   20 CONTINUE
*
*       Transform to new regularized momentum (eqn (50)).
      K = K1
      P(K+1) = 2.0D0*(+Q(K+1)*P1(K+1) + Q(K+2)*P1(K+2) +
     &                                  Q(K+3)*P1(K+3))
      P(K+2) = 2.0D0*(-Q(K+2)*P1(K+1) + Q(K+1)*P1(K+2) +
     &                                  Q(K+4)*P1(K+3))
      P(K+3) = 2.0D0*(-Q(K+3)*P1(K+1) - Q(K+4)*P1(K+2) +
     &                                  Q(K+1)*P1(K+3))
      P(K+4) = 2.0D0*(+Q(K+4)*P1(K+1) - Q(K+3)*P1(K+2) +
     &                                  Q(K+2)*P1(K+3))
*
*       Update the smallest distance (new pericentre has increased).
      IF (IM.EQ.1) THEN
          R1 = RI
      ELSE
          R2 = RI
      END IF
*
*       Transform to physical variables and obtain consistent energy.
      CALL TRANS3(3)
      CALL TRANS3(0)
*
*       Update diagnostic variables (note DE > 0).
      ECOLL3 = ECOLL3 + (DE(1) + DE(2))
      NDISS3 = NDISS3 + 1
*
*       Perform stability test (ITERM < 0 denotes termination).
      CALL STABL3(ITERM)
*
*       Print diagnostic if eccentricity > 0.99.
      IF (ECC.GT.0.99) THEN
          if(rank.eq.0)
     &    WRITE (6,30)  NAME3(IM), NAME3(3), SEMI1, ECC, ECC1, H, QPERI
   30     FORMAT (' NEW QPMOD3   NAM AF E0 EF H QP ',
     &                              2I5,1P,E10.2,0P,2F8.3,F9.1,1P,E10.2)
          CALL FLUSH(6)
      END IF
*
      RETURN
*
      END
