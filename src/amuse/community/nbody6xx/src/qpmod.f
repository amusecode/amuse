      SUBROUTINE QPMOD(IM,ITERM)
*
*
*       Modification of chain variables for tidal dissipation.
*       ------------------------------------------------------
*
      INCLUDE 'commonc.h'
      INCLUDE 'common2.h'
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/CCOLL2/  QK(NMX4),PK(NMX4),RIJ(NMX,NMX),SIZE(NMX),VSTAR1,
     &                ECOLL1,RCOLL,QPERI,ISTAR(NMX),ICOLL,ISYNC,NDISS1
      COMMON/KSAVE/  K1,K2
      REAL*8  G0(3),RADIUS(2),DE(2)
      INTEGER IS(2)
      DATA  ECCM,ECCM2  /0.002,0.00000399/
      INCLUDE "mpif.h"
      INTEGER group,rank,ierr,isize,status(MPI_STATUS_SIZE)
      COMMON/MPIDAT/group,rank,ierr,isize,status
*
*
*       Skip for small dissipation (ISYNC > 0 delays further calls).
      IF (KZ27.EQ.1) THEN
          R1 = MAX(SIZE(K1),SIZE(K2))
          IF (ABS(QPERI - 4.0*R1).LT.0.01*QPERI) THEN
              ISYNC = 1
              ITERM = 0
              GO TO 90
          END IF
      END IF
*
*       Copy radius & stellar type for interacting bodies #K1 & K2.
      RADIUS(1) = SIZE(K1)
      RADIUS(2) = SIZE(K2)
      IS(1) = ISTAR(K1)
      IS(2) = ISTAR(K2)
      KSTARI = 0
*
*       Evaluate binding energy & semi-major axis from non-singular terms.
      CALL EREL(IM,EB,SEMI)
*
*       Define mass & reduced mass for the dominant bodies.
      MB = M(K1) + M(K2)
      MU = M(K1)*M(K2)/MB
*
*       Set energy per unit mass, eccentricity & pericentre (assume R' = 0).
      H = EB/MU
      ECC = 1.0 - 1.0/(RINV(IM)*SEMI)
      PERI = 1.0/RINV(IM)
*     PERI = SEMI*(1.0D0 - ECC)
      AM0 = SEMI*(1.0D0 - ECC**2)
      NDISS1 = NDISS1 + 1
*
*       Choose between chaos treatment and PT or GR formulation.
      IF (KZ27.EQ.2) THEN
          IDIS = 0
          CALL CHAOS2(K1,K2,ECC,H,IS,MB,MU,RADIUS,SEMI1,ECC1,DH,IDIS,
     &                                                           KSTARI)
*       Exit on reaching circular orbit or spiralling stage (also collision).
          IF (KSTARI.EQ.-1) GO TO 2
          GO TO 90
      END IF
*
*       Consider sequential circularization or GR case.
      IF (KZ27.EQ.1) THEN
          AM0 = SEMI*(1.0D0 - ECC**2)
          ECC2 = ECCM2
          ECC1 = SQRT(ECC2)
          ACIRC = AM0/(1.0 - ECC2)
*       Accept circularized orbit directly if ACIRC < 4*R.
          IF (ACIRC.LT.4.0*R1) THEN
              SEMI1 = ACIRC
          ELSE
*       Obtain E1 by (1 + E1) = AM0/(4*R1) and A1 by A1*(1 - E1) = 4*R1.
              ECC1 = 0.25*AM0/R1 - 1.0
              ECC1 = MAX(ECC1,ECCM)
*       Set semi-major axis from angular momentum conservation.
              SEMI1 = AM0/(1.0 - ECC1**2)
          END IF
*       Form the corresponding energy change.
          DH = 0.5*MB*(1.0/SEMI - 1.0/SEMI1)
          DE(1) = -MU*DH
          DE(2) = 0.0
      ELSE
*       Obtain the tidal energy change for GR.
          CALL TIDES3(QPERI,M(K1),M(K2),VSTAR1,H,ECC,DE)
*       Include safety check on energy loss to prevent new SEMI < R.
          DH = -(DE(1) + DE(2))/MU
          IF (H + DH.LT.-0.5*MB*RINV(IM)) THEN
              DH = -0.5*MB*RINV(IM) - H
          END IF
          SEMI1 = -0.5*MB/(H + DH)
          ECC1 = 1.0 - PERI/SEMI1
          ECC1 = MAX(ECC1,0.0D0)
*       Note: Do not impose minimum ECCM unless with 2nd alternative C2.
*         ECC1 = MAX(ECC1,ECCM)
      END IF
*
*       Determine new eccentricity from angular momentum conservation.
*     ECC2 = ECC**2 + 2.0D0*AM0*DH/MB
*     ECC2 = MAX(ECC2,ECCM2)
*     ECC1 = SQRT(ECC2)
*
*       Adopt instantaneous circularization instead of standard PT.
*     ECC2 = ECCM2
*     ECC1 = SQRT(ECC2)
*     SEMI1 = AM0/(1.0 - ECC2)
*     DH = 0.5*MB*(1.0/SEMI - 1.0/SEMI1)
*     DE(1) = -MU*DH
*     DE(2) = 0.0

*       Skip on final hyperbolic energy.
      IF (H + DH.GT.0.0) GO TO 90
*
*       Set new pericentre and update binding energy.
    2 PERI1 = SEMI1*(1.0D0 - ECC1)
      HI = H
      H = H + DH
*
*       Form KS coordinate scaling factor from pericentre ratio.
      C1 = SQRT(PERI1/PERI)
*
*       Specify KS velocity scaling for conserved angular momentum.
      IF (KZ27.EQ.1.OR.KSTARI.EQ.-2) THEN
          C2 = 1.0/C1**2
      ELSE
          C2 = SQRT((MB + H*PERI1)/(MB + HI*PERI))
*       Note: since PERI1 = PERI this is same as used for KZ27 = 2.
      END IF
*
*       See whether circular orbit condition applies.
      IF (ECC1.LE.ECCM.AND.KZ27.EQ.1) THEN
          AM = SEMI1*(1.0D0 - ECC1**2)
          C2 = (AM/AM0)/C1**2
      END IF
*
*       Include alternative formulation of Mikkola (which is equivalent).
      IF (KZ27.LT.-1) THEN
*       Note this derivation may contain assumption of J = const (10/6/99).
          V02 = MB*(2.0/QPERI - 1.0/SEMI)
          A1 = MB/(QPERI*V02)
          A2 = MB/(V02*SEMI1)
*       Adopt reduced energy loss in case of imaginary solution.
          IF (A1**2.GT.A2) THEN
              A3 = A1 + SQRT(A1**2 - A2)
          ELSE
              A3 = A1
              if(rank.eq.0)WRITE (6,5)  A1, A2, SQRT(1.0/A3), A3, C1
    5         FORMAT (' WARNING!    QPMOD    A1 A2 SQRT(1/A3) A3 C1  ',
     &                                       1P,2E10.2,0P,3F12.6)
          END IF
          C2 = A3
      END IF
*
*       Derive scale factor from velocity ratio for chaotic motion.
      IF (KZ27.EQ.2.AND.KSTARI.EQ.-1) THEN
          C2 = SQRT((H + MB/PERI1)/(HI + MB/PERI))
*       Note the use of physical velocity here vs UDOT in routine KSTIDE.
      END IF
*
*       Modify the KS coordinates for dominant bodies and update distance.
      KS = 4*(IM - 1)
      RM = 0.0D0
      DO 10 K = 1,4
          Q(KS+K) = C1*Q(KS+K)
          RM = RM + Q(KS+K)**2
   10 CONTINUE
      RINV(IM) = 1.0/RM
*
*       Change the corresponding physical momenta (IM & IM + 1 in chain).
      J1 = 3*(IM - 1)
      J2 = 3*IM
      DO 20 K = 1,3
*         P11 = A3*(M(K2)*PI(J1+K) - M(K1)*PI(J2+K))/MB
          VB = (PI(J1+K) + PI(J2+K))/MB
*         P1K = -MU*(1.0 - C2)*(PI(J1+K)/M(K1) - PI(J2+K)/M(K2))
          P1K = C2*(M(K2)*PI(J1+K) - M(K1)*PI(J2+K))/MB
          PI(J1+K) = M(K1)*VB + P1K
          PI(J2+K) = M(K2)*VB - P1K
*         PI(J1+K) = PI(J1+K) + P1K
*         PI(J2+K) = PI(J2+K) - P1K
   20 CONTINUE
*
*       Form physical chain momenta (first & last, then intermediate values).
      L = 3*(N - 2)
      DO 30 K = 1,3
          WC(K) = -PI(K)
          WC(L+K) = PI(L+K+3)
   30 CONTINUE
*
      DO 40 I = 2,N-2
          L = 3*(I - 1)
          DO 35 K = 1,3
              WC(L+K) = WC(L+K-3) - PI(L+K)
   35     CONTINUE
   40 CONTINUE
*
*       Re-determine all regularized momenta by KS transformation.
      DO 50 L = 1,N-1
          L1 = 3*(L - 1) + 1
          L2 = L1 + 1
          L3 = L2 + 1
          KS1 = 4*(L - 1) + 1
          KS2 = KS1 + 1
          KS3 = KS2 + 1
          KS4 = KS3 + 1
          P(KS1) = 2.D0*(+Q(KS1)*WC(L1) + Q(KS2)*WC(L2) + Q(KS3)*WC(L3))
          P(KS2) = 2.D0*(-Q(KS2)*WC(L1) + Q(KS1)*WC(L2) + Q(KS4)*WC(L3))
          P(KS3) = 2.D0*(-Q(KS3)*WC(L1) - Q(KS4)*WC(L2) + Q(KS1)*WC(L3))
          P(KS4) = 2.D0*(+Q(KS4)*WC(L1) - Q(KS3)*WC(L2) + Q(KS2)*WC(L3))
   50 CONTINUE
*
*       Evaluate potential energy due to non-chained distances (KS terms OK).
      IT = 0
   52 POT1 = 0.0
      IJ = N - 1
      DO 60 I = 1,N-2
          DO 55 J = I+2,N
              IJ = IJ + 1
              POT1 = POT1 + MIJ(I,J)*RINV(IJ)
   55     CONTINUE
   60 CONTINUE
*
*       Obtain the potential energy of the new system (KS terms unchanged).
      IF (IT.EQ.0) THEN
*       Save modified configuration in QK & PK (for EREL & TRANSK).
          DO 70 I = 1,N-1
              KS = 4*(I - 1)
              DO 65 J = 1,4
                  QK(KS+J) = Q(KS+J)
                  PK(KS+J) = P(KS+J)
   65         CONTINUE
   70     CONTINUE
*
*       Set new values of inverse distances (only non-chained values needed).
          CALL TRANSK
          POT0 = POT1
          IT = IT + 1
          GO TO 52
      END IF
*
*       Update collision energy & internal energy (MU*DH is not sufficient).
      ECOLL1 = ECOLL1 - MU*DH + (POT1 - POT0)
      ENERGY = ENERGY + MU*DH - (POT1 - POT0)
*
      IF (KZ30.GT.1) THEN
*       Obtain consistent value of the total energy (for checking purposes).
          CALL TRANSX
          CALL CONST(X,V,M,N,EN1,G0,ALAG)
          DPOT = POT1 - POT0
          if(rank.eq.0)WRITE (6,75)EN1,EN1-ENERGY,MU*DH,SEMI,SEMI1,DPOT
   75     FORMAT (/,' QPMOD:   EN1 ERR DE A A1 DP ',
     &                         F10.6,1P,5E10.2)
      END IF
*
*       Activate indicator for synchronous orbit.
      IF (ECC.GT.ECCM.AND.ECC1.LE.ECCM) THEN
          ISYNC = 1
      END IF
*
*       Print diagnostic if eccentricity > 0.99 or #30 > 1.
      IF (ECC.GT.0.99.OR.KZ30.GT.1) THEN
          if(rank.eq.0)
     &        WRITE (6,80)NAMEC(K1),NAMEC(K2),SEMI1,ECC,ECC1,H, QPERI
   80     FORMAT (' NEW QPMOD    NAM AF E0 EF H QP ',
     &                           2I6,1P,E10.2,0P,2F8.4,F9.1,1P,2E10.2)
      END IF
*
*       Perform stability test (ITERM = -1: termination; = -2: reduction).
      CALL STABLC(IM,ITERM,SEMI1)
*
   90 RETURN
*
      END
