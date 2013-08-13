      SUBROUTINE RENEW(I)
*
*
*       Re-activation of dormant KS binary.
*       -----------------------------------
*
      INCLUDE 'common6.h'
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
*
*
*       Restore masses of binary components and initialize time.
      IPAIR = I - N
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      BODY(I1) = BODYC(9)
      BODY(I2) = BODYC(10)
      T0(I1) = TIME
*
*       Specify semi-major axis and Kepler period.
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
*
*       Estimate consistent unperturbed interval.
      CALL TPERT(IPAIR,GMIN,DT)
*
*       Specify unperturbed motion (subject to DT > 0).
      K = 1 + INT(0.5D0*DT/TK)
      K = MAX(K,1)
      STEP(I1) = FLOAT(K)*TK
*
*       Obtain coordinates & velocities of unperturbed binary.
      CALL RESOLV(IPAIR,1)
*
*       Copy global indices of all chain members to JPERT for routine NBPOT.
      DO 1 L = 1,NCH
          JPERT(L) = JLIST(L)
    1 CONTINUE
*
*       Evaluate potential energy with respect to inner components.
      JLIST(1) = I1
      JLIST(2) = I2
      CALL NBPOT(2,NCH,POT1)
*
*       Include interaction of inner c.m. & other members to give net effect.
      JLIST(1) = I
      CALL NBPOT(1,NCH,POT2)
      DPHI = POT1 - POT2
*
*       Form new KS polynomials if nearest perturber is significant (> GMIN).
      IF (DT.LT.0.0) THEN
          CALL KSLIST(IPAIR)
          CALL KSPOLY(IPAIR,1)
      END IF
*
*       Subtract binding energy from temporary save in ECOLL.
      EB = BODY(I1)*BODY(I2)*H(IPAIR)/BODY(I)
      ECOLL = ECOLL - EB + DPHI
*
      IF (rank.eq.0.and.KZ(30).GT.1) THEN
          WRITE (6,5)  IPAIR, EB, R(IPAIR), GAMMA(IPAIR), DPHI, K
    5     FORMAT (' RENEW:    KS =',I3,'  EB =',F8.2,'  R =',1P,E8.1,
     &            '  G =',E8.1,'  DPHI =',E8.1,'  DT/TK =',0P,I4)
      END IF
*
      RETURN
*
      END
