      SUBROUTINE KSLIST(IPAIR)
*
*
*       KS perturber selection.
*       -----------------------
*
      INCLUDE 'common6.h'
*
*
*       Set component & c.m. index and form semi-major axis & eccentricity.
      I1 = 2*IPAIR - 1
      I = N + IPAIR
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      EB = -0.5*BODY(I1)*BODY(I1+1)/SEMI
      ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      ECC = SQRT(ECC2)
      IT = 0
*
*       Use semi-major axis and/or RMIN for perturber selection.
      IF (EB.LT.EBH) THEN
          RAP = SEMI*(1.0 + ECC)
      ELSE
*       Include a tapered criterion depending on energy for soft binaries.
          IF (EB.LT.0.0) THEN
              ZFAC = 1.0 + ABS(EB - EBH)/ABS(EBH)
          ELSE
              ZFAC = 1.0
          END IF
*       Adopt actual apocentre for perturber selection if R > SEMI.
          IF (SEMI.GT.0.0D0.AND.R(IPAIR).GT.SEMI) THEN
              RAP = SEMI*(1.0 + ECC)
          ELSE
              RAP = MAX(ZFAC*RMIN,R(IPAIR))
          END IF
*       Ensure extra perturbers at new regularizations (R may be small).
          IF (IPHASE.GT.0.AND.SEMI.GT.0.0) THEN
              RAP = SEMI*(1.0 + ECC)
          END IF
      END IF
*
      RAP = MIN(RAP,20.0*RMIN)
      RCRIT2 = 2.0*RAP**2/BODY(I)
      RCRIT3 = RCRIT2*RAP/GMIN
*       Base fast search on maximum binary mass (2*BODY1).
      RCRIT2 = 2.0*RCRIT2*BODY1*CMSEP2
*
*       Select new perturbers from the neighbour list.
    6 NNB1 = 1
      NNB2 = LIST(1,I) + 1
      DO 10 L = 2,NNB2
          J = LIST(L,I)
          W1 = X(1,J) - X(1,I)
          W2 = X(2,J) - X(2,I)
          W3 = X(3,J) - X(3,I)
          RSEP2 = W1*W1 + W2*W2 + W3*W3
*       Include any merged c.m. or chain c.m. bodies in the fast test.
          IF (RSEP2.LT.RCRIT2.OR.NAME(J).LE.0) THEN
              RIJ3 = RSEP2*SQRT(RSEP2)
*       Estimate unperturbed distance from tidal limit approximation.
              IF (RIJ3.LT.BODY(J)*RCRIT3) THEN
                  NNB1 = NNB1 + 1
                  LIST(NNB1,I1) = J
              ELSE IF (J.GT.N) THEN
*       Employ a more generous criterion for possible wide binary.
                  RJ = BODY(J)/H(J-N)
                  IF (RSEP2.LT.CMSEP2*RJ**2) THEN
                      NNB1 = NNB1 + 1
                      LIST(NNB1,I1) = J
                  END IF
              END IF
          END IF
   10 CONTINUE
*
*       Ensure at least one perturber first time (max 5 tries except CHAOS).
      IF (NNB1.EQ.1.AND.IPHASE.GT.0.AND.NNB2.GT.1.AND.TIME.GT.0.0) THEN
          RCRIT2 = 2.0*RCRIT2
          RCRIT3 = 2.0*RCRIT3
          IT = IT + 1
*       Skip repeat for small size (next KSLIST requires many periods).
          IF ((SEMI*SU.GT.10.0.AND.IT.LE.5).OR.KSTAR(I).EQ.-1) GO TO 6
      END IF
*
*       Check case of no perturbers (dual purpose).
      IF (NNB1.EQ.1) THEN
*       Add distant perturber for hyperbolic orbit.
          IF (SEMI.LT.0.0) THEN
              NNB1 = 2
              LIST(2,I1) = LIST(2,I)
              GO TO 20
          END IF
*       Retain the previous perturber after partial reflection.
*         IF (TDOT2(IPAIR).GT.0.0D0.AND.KZ(25).GT.0) THEN
*             NNB1 = 2
*         ELSE IF (KZ(27).LE.0) THEN
          IF (KZ(27).LE.0) THEN
*       Specify one unperturbed period at apocentre (NB! check STEP(I)).
              STEP(I1) = TWOPI*SEMI*SQRT(SEMI/BODY(I))
              STEP(I1) = MIN(STEP(I1),STEP(I))
          ELSE
*       Maintain perturbed motion during Chaos event (not ROCHE/SPIRAL).
              IF (KSTAR(I).EQ.-1) THEN
                  IF (LIST(1,I1).GT.0) THEN
                      NNB1 = 2
                      LIST(2,I1) = N
                  END IF
              ELSE
                  STEP(I1) = TWOPI*SEMI*SQRT(SEMI/BODY(I))
                  STEP(I1) = MIN(STEP(I1),STEP(I))
              END IF
          END IF
      END IF
*
*       Save perturber membership.
   20 LIST(1,I1) = NNB1 - 1
*
      RETURN
*
      END
