      SUBROUTINE NBSORT(NBL,IBL,NNB,NBLIST)
*
*
*       Sorting of neighbour lists.
*       ---------------------------
*
      INCLUDE 'common6.h'
      INTEGER  IBL(LMAX),NBLIST(NMAX),LP(NMAX)
*
*
*       Initialize counter and list pointers.
      NNB = 0
      DO 1 L = 1,NBL
          LP(L) = 2
    1 CONTINUE
*
*       Reset terminator count and minimum particle index.
    5 IEND = 0
      JMIN = NTOT + 1
*
*       Compare current location of all neighbour lists.
      DO 10 LL = 1,NBL
          I = IBL(LL)
          L = LP(LL)
*       Check each list member or count completed searches.
          IF (L.LT.LIST(1,I) + 2) THEN
              J = LIST(L,I)
              IF (J.LT.JMIN) THEN
                  JMIN = J
                  LMIN = LL
*       Update pointer for all J <= JMIN (except latest JMIN).
              ELSE IF (J.LE.JMIN) THEN
                  LP(LL) = LP(LL) + 1
              END IF
          ELSE
              IEND = IEND + 1
          END IF
   10 CONTINUE
*
*       Increase pointer & membership and save new index in NBLIST.
      IF (IEND.LT.NBL) THEN
          LP(LMIN) = LP(LMIN) + 1
          NNB = NNB + 1
          NBLIST(NNB) = JMIN
          GO TO 5
      END IF
*
*       Exit on rare case of one member with zero neighbours (10/2008).
      IF (NNB.EQ.0) GO TO 70
*
*       Add any integration particles missing from joint neighbour lists.
      NEW = 0
      FAC = NNB/FLOAT(NTOT)
      L = 1
   20 I = IBL(L)
*       Skip search if body #I is outside the list range.
      IF (NBLIST(1).GT.I.OR.NBLIST(NNB).LT.I) GO TO 50
      IF (NBLIST(1).EQ.I.OR.NBLIST(NNB).EQ.I) GO TO 60
*
*       Search sequential list using incremental expectation values (> 0).
      IG = I*FAC
      IG = MAX(IG,1)
   30 IF (NBLIST(IG).GT.I) THEN
          INC = (NBLIST(IG) - I)*FAC
          INC = MIN(INC,3)
          IF (INC.LE.1) THEN
              IG = IG - 1
              IF (NBLIST(IG).LT.I) GO TO 50
          ELSE
              IG = IG - INC
              IG = MAX(IG,1)
              IF (NBLIST(IG).LT.I) THEN
                  IG = IG + INC - 1
                  IF (NBLIST(IG).GT.I) THEN
                      IG = IG - 1
                      IF (NBLIST(IG).GT.I) IG = IG - 1
                  END IF
              ELSE
                  IF (NBLIST(IG).GT.I) THEN
                      IG = IG - 1
                      IF (NBLIST(IG).GT.I) IG = IG - 1
                  END IF
              END IF
          END IF
          GO TO 30
      END IF
*
*       Treat case of NBLIST < I until boundary exceeded or NBLIST > I.
   40 IF (NBLIST(IG).LT.I) THEN
          INC = (I - NBLIST(IG))*FAC
          INC = MAX(INC,1)
          IG = IG + INC
          IF (IG.GE.NNB) THEN
              IF (NBLIST(NNB).EQ.I) GO TO 60
              GO TO 50
          END IF
          IF (INC.EQ.1) THEN
              IF (NBLIST(IG).GT.I) GO TO 50
          END IF
*       Reduce index by 1 to avoid resonance oscillations (NBLIST > I).
   45     IF (NBLIST(IG).GT.I) THEN
              IG = IG - 1
              IF (NBLIST(IG).LT.I) GO TO 50
              GO TO 45
          END IF
          GO TO 40
      END IF
*
*       Skip addition if body #I is already a member (NBLIST(IG) = I).
      GO TO 60
*
*       Increase counter and add body #I at end of array.
   50 NEW = NEW + 1
      NBLIST(NNB+NEW) = I
*
*       Continue until last body has been checked.
   60 L = L + 1
      IF (L.LE.NBL) GO TO 20
*
*       Specify the total membership (Note: new members not sequential).
      NNB = NNB + NEW
*
   70 RETURN
*
      END
