      SUBROUTINE CHLIST(I)
*
*
*       Perturber list for chain regularization.
*       -----------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  M,MASS,MC,MIJ,MKK
      PARAMETER  (NMX=10,NMX3=3*NMX,NMXm=NMX*(NMX-1)/2)
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CPERT/  RGRAV,GPERT,IPERT,NPERT
*
*
*       Prepare perturber selection using maximum of RSUM & RGRAV.
      RPERT = MAX(RSUM,RGRAV)
*       Restrict effective size to RMIN (energy may be near zero).
      RPERT = MIN(RPERT,RMIN)
      RCRIT2 = 2.0*RPERT**2/BODY(I)
      RCRIT3 = RCRIT2*RPERT/GMIN
*       Base fast search on maximum binary mass (2*BODY1).
      RCRIT2 = 2.0*RCRIT2*BODY1*CMSEP2
      PMAX = 0.0
*
*       Select new perturbers from chain c.m. neighbour list.
      NNB1 = 1
      NNB2 = LIST(1,I) + 1
      DO 10 L = 2,NNB2
          J = LIST(L,I)
          W1 = X(1,J) - X(1,I)
          W2 = X(2,J) - X(2,I)
          W3 = X(3,J) - X(3,I)
          RSEP2 = W1*W1 + W2*W2 + W3*W3
*       Include any merged c.m. bodies in the fast test.
          IF (RSEP2.LT.RCRIT2.OR.NAME(J).LT.0) THEN
              RIJ3 = RSEP2*SQRT(RSEP2)
*       Estimate perturbed distance from tidal limit approximation.
              IF (RIJ3.LT.BODY(J)*RCRIT3) THEN
                  NNB1 = NNB1 + 1
                  LISTC(NNB1) = J
                  PMAX = MAX(BODY(J)/RIJ3,PMAX)
              END IF
          END IF
   10 CONTINUE
*
*       Save perturber membership.
      LISTC(1) = NNB1 - 1
*
*       Form current perturbation from nearest body (effective value RPERT).
      GPERT = 2.0*PMAX*RPERT**3/BODY(ICH)
*
*       Obtain actual perturbation (note GAMX zero first time).
      CALL CHPERT(GAMX)
*
*     if(rank.eq.0)
*    &WRITE (6,20)  ICH, NNB2-1, LISTC(1), RSUM, RPERT, GPERT, GAMX
*  20 FORMAT (' CHLIST:   ICH NB NP RSUM RPERT GPERT GX ',3I4,1P,4E9.1)
*
      IF (GAMX.GT.0.0D0) THEN
          GPERT = GAMX
      END IF
      NPERT = LISTC(1)

*       Ensure consistent perturbation indicator after re-determination.
      IF (GPERT.GT.1.0D-06.AND.NPERT.GT.0) THEN
          IPERT = 1
      ELSE
          IPERT = 0
      END IF
*
      RETURN
*
      END
