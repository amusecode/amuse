      SUBROUTINE QPMOD4(IM,ITERM)
*
*
*       Modification of chain variables for tidal dissipation.
*       ------------------------------------------------------
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      COMMON/CREG/  M(4),X(12),XD(12),P(12),Q(12),TIME4,ENERGY,EPSR2,
     &              XR(9),W(9),R(6),TA(6),MIJ(6),CM(10),RMAX4,TMAX,
     &              DS,TSTEP,EPS,NSTEP4,NAME4(4),KZ15,KZ27,NREG,NFN
      COMMON/ICONF/  I1,I2,I3,I4
      COMMON/CLOSE/  RIJ(4,4),RCOLL,QPERI,SIZE(4),ECOLL3,IP(4)
      COMMON/CCOLL/  QK(12),PK(12),ICALL,ICOLL,NDISS4
      COMMON/SAVEP/  PI(12)
      COMMON/KSAVE/  K1,K2
      REAL*8  DE(2),RADIUS(2)
      INTEGER IS(2)
      INCLUDE "mpif.h"
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      COMMON/MPIDAT/group,rank,ierr,isize,status
*
*
*       Obtain tidal energy loss (bodies #K1 & K2 with separation QPERI).
      RADIUS(1) = SIZE(K1)
      RADIUS(2) = SIZE(K2)
      IS(1) = IP(K1)
      IS(2) = IP(K2)
*
      CALL TIDES(QPERI,M(K1),M(K2),RADIUS(1),RADIUS(2),IS,DE)
*
*       Evaluate binding energy & semi-major axis from non-singular terms.
      CALL EREL4(IM,EB,SEMI)
*
*       Define mass & reduced mass for the dominant bodies.
      MB = M(K1) + M(K2)
      MU = M(K1)*M(K2)/MB
*
*       Set energy per unit mass, eccentricity & pericentre (assume R' = 0).
      H = EB/MU
      ECC = 1.0 - R(IM)/SEMI
      PERI = SEMI*(1.0D0 - ECC)
*
*       Determine new eccentricity from angular momentum conservation.
      DH = -(DE(1) + DE(2))/MU
      AM0 = SEMI*(1.0D0 - ECC**2)
      ECC2 = ECC**2 + 2.0D0*AM0*DH/MB
      IF (ECC2.GT.0.0D0) THEN
          ECC1 = SQRT(ECC2)
      ELSE
          ECC1 = 0.001
      END IF
*
*       Update binding energy and set new semi-major axis & pericentre.
      H1 = H + DH
      SEMI1 = -0.5D0*MB/H1
      PERI1 = SEMI1*(1.0D0 - ECC1)
*
*       Form KS coordinate scaling factor from pericentre ratio.
      C1 = SQRT(PERI1/PERI)
*
*       Specify velocity scaling from angular momentum conservation.
      C2 = 1.0/C1
*
*       See whether circular orbit condition applies.
      IF (ECC1.LE.0.001) THEN
          AM = SEMI1*(1.0D0 - ECC1**2)
          C2 = SQRT(AM/AM0)/C1
      END IF
*
*       Modify the KS coordinates for dominant bodies and update distance.
      KS = 4*(IM - 1)
      R(IM) = 0.0D0
      DO 20 K = 1,4
          Q(KS+K) = C1*Q(KS+K)
          R(IM) = R(IM) + Q(KS+K)**2
   20 CONTINUE
*
*       Set modified physical momenta for the two critical bodies (K1 & K2).
      J1 = 3*(K1 - 1)
      J2 = 3*(K2 - 1)
*       Note PI is associated with the bodies (hence J1 # 3*(IM - 1).
      DO 30 K = 1,3
          P1K = -MU*(1.0 - C2**2)*(PI(J1+K)/M(K1) - PI(J2+K)/M(K2))
          PI(J1+K) = PI(J1+K) + P1K
          PI(J2+K) = PI(J2+K) - P1K
   30 CONTINUE
*
*       Form physical chain momenta. 
      IP1 = 3*(I1 - 1)
      IP2 = 3*(I2 - 1)
      IP4 = 3*(I4 - 1)
      DO 40 K = 1,3
          W(K  ) = -PI(IP1+K)
          W(K+3) = -PI(IP1+K) - PI(IP2+K)
          W(K+6) = +PI(IP4+K)
   40 CONTINUE
*
*       Re-determine all regularized momenta by KS transformation.
      DO 50 L = 1,3
          L1 = 3*(L - 1) + 1
          L2 = L1 + 1
          L3 = L2 + 1
          LQ1 = 4*(L - 1) + 1
          LQ2 = LQ1 + 1
          LQ3 = LQ2 + 1
          LQ4 = LQ3 + 1
          P(LQ1) = 2.D0*(+Q(LQ1)*W(L1) + Q(LQ2)*W(L2) + Q(LQ3)*W(L3))
          P(LQ2) = 2.D0*(-Q(LQ2)*W(L1) + Q(LQ1)*W(L2) + Q(LQ4)*W(L3))
          P(LQ3) = 2.D0*(-Q(LQ3)*W(L1) - Q(LQ4)*W(L2) + Q(LQ1)*W(L3))
          P(LQ4) = 2.D0*(+Q(LQ4)*W(L1) - Q(LQ3)*W(L2) + Q(LQ2)*W(L3))
   50 CONTINUE
*
*       Obtain consistent value of the total energy (instead of correcting).
      CALL ENDREG
      E0 = ENERGY
      CALL NEWSYS(X,XD,M,4,ENERGY,GAM)
*
*       Update diagnostic variables (Note: DE(1) + DE(2) is not sufficient).
*     ECOLL3 = ECOLL3 + (DE(1) + DE(2))
      ECOLL3 = ECOLL3 + (E0 - ENERGY)
      NDISS4 = NDISS4 + 1
*
*       Perform stability test (ITERM < 0 denotes termination).
      CALL STABL4(ITERM)
*
*       Print diagnostic if eccentricity > 0.99.
      IF (ECC.GT.0.99) THEN
          if(rank.eq.0)
     &    WRITE (6,60)  NAME4(K1), NAME4(K2), SEMI1, ECC, ECC1, H, QPERI
   60     FORMAT (' NEW QPMOD4   NAM AF E0 EF H QP ',
     &                              2I5,1P,E10.2,0P,2F8.3,F9.1,1P,E10.2)
      END IF
*
      RETURN
*
      END
