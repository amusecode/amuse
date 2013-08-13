      SUBROUTINE UPDATE(IPAIR)
*
*
*       List modifications after KS termination.
*       ----------------------------------------
*
      INCLUDE 'common6.h'
      DATA IW /0/
      SAVE IW
*
*
*       Adjust neighbour lists to new sequence (skip last or only pair).
      ICM = N + IPAIR
      IF (IPAIR.GT.NPAIRS) GO TO 60
*
*       Set largest index of any regularized member.
      ILAST = 2*NPAIRS + 2
*
*       Rename lists containing single regularized components.
      DO 50 J = 1,NTOT
*       Skip renaming if first neighbour exceeds regularized components.
          IF (LIST(2,J).GT.ILAST) GO TO 50
          NNB = LIST(1,J)
          IF (NNB.EQ.0) GO TO 50
          L = 2
*       Determine list index and new location of any regularized component.
          INEW = 0
   10     IF (LIST(L,J).EQ.2*IPAIR-1) THEN
              I = 2*IPAIR - 1
              INEW = 2*NPAIRS + 1
          ELSE IF (LIST(L,J).EQ.2*IPAIR) THEN
              I = 2*IPAIR
              INEW = 2*NPAIRS + 2
          ELSE IF (INEW.EQ.0) THEN
*       Define as zero for single body test below (at least one < ILAST).
              I = 0
*       Note that L = 2 may be single body and L = 3 current KS component.
          END IF
          L = L + 1
          IF (L.LE.NNB+1.AND.LIST(L,J).LT.ILAST) GO TO 10
*
*       Rename regularized components > 2*IPAIR and reduce by 2 if <= ILAST.
          DO 20 K = 2,L
*       Note that L determined above is list index of standard particle.
              IF (LIST(K,J).EQ.I) THEN
                  LIST(K,J) = INEW
              ELSE IF (LIST(K,J).LE.ILAST.AND.
     &                 LIST(K,J).GT.2*IPAIR) THEN
                  LIST(K,J) = LIST(K,J) - 2
              END IF
   20     CONTINUE
*
*       Check that list of single KS components is sequential (up to ILAST).
          L = 2
   30     IF (LIST(L+1,J).LT.LIST(L,J).AND.L.LE.NNB) THEN
              JL = LIST(L,J)
              LIST(L,J) = LIST(L+1,J)
              LIST(L+1,J) = JL
          END IF
          L = L + 1
          IF (LIST(L,J).LE.ILAST.AND.L.LE.NNB+1) GO TO 30
   50 CONTINUE
*
*       Replace c.m. by components and reduce subsequent members by one.
   60 DO 80 J = 1,NTOT
          NNB = LIST(1,J)
          L = NNB + 1
          MOVE = 0
*       Check for removal of current c.m. and reduce by one if > N + IPAIR.
   70     IF (L.EQ.1.OR.LIST(L,J).LT.ICM) GO TO 80
*       See whether ICM is on neighbour list (otherwise skip splitting).
          IF (LIST(L,J).EQ.ICM) MOVE = 1
          IF (LIST(L,J).NE.ICM) THEN
              LIST(L,J) = LIST(L,J) - 1
              L = L - 1
              GO TO 70
          END IF
*
*       Skip on zero ID or remove c.m. name by moving up subsequent members.
          IF (MOVE.EQ.0) GO TO 80
   72     IF (L.LE.NNB) THEN
              LIST(L,J) = LIST(L+1,J) 
              L = L + 1
              GO TO 72
          END IF
*
          NNB = NNB - 1
*       Expand the list to include both components since c.m. was deleted.
          KCASE = 2
*       Only move neighbours down by one if the list has too many members.
          IF (NNB.GT.LMAX-3) KCASE = 1
          IF (NNB.EQ.0) GO TO 76
*       In this special case L = 2 already.
          L = NNB + 1
*       Take special precaution if last neighbour is a regularized body.
          IF (LIST(L,J).LE.JCOMP) THEN
              L = L + 1
              GO TO 76
          END IF
*
*       Move members down by two (or one) to make room for components.
   74     LIST(L+KCASE,J) = LIST(L,J)
          IF (L.GT.2.AND.LIST(L-1,J).GE.ICOMP) THEN
              L = L - 1
              GO TO 74
          END IF
*
*       Rename deleted c.m. appropriately and increase membership by 2 or 1.
   76     LIST(L,J) = ICOMP
*       Do not over-write the list if NNB > NNBMAX after removal of c.m.
          IF (KCASE.EQ.2) LIST(L+1,J) = JCOMP
          LIST(1,J) = NNB + KCASE
          IF (KCASE.EQ.1.AND.IW.LT.10) THEN
              if(rank.eq.0)
     &        WRITE (6,78)  NNB, J, JCOMP
   78         FORMAT (5X,'WARNING!    UPDATE    NNB J JCOMP ',3I6)
              IW = IW + 1
          END IF
   80 CONTINUE
*
*       Modify the list of previously regularized binaries.
      NNB = LISTR(1) - 1
      L = 0
   91 L = L + 2
   92 IF (L.GT.NNB + 1) GO TO 96
      J = LISTR(L)
      K = LISTR(L+1)
*       First check the current two-body separation of any old pairs.
      RJK2 = (X(1,J) - X(1,K))**2 + (X(2,J) - X(2,K))**2 +
     &                              (X(3,J) - X(3,K))**2
*       Remove pair if RJK > 4*RMIN when special procedure not needed.
      IF (RJK2.LT.16.0*RMIN**2) GO TO 91
      DO 94 K = L,NNB
          LISTR(K) = LISTR(K+2)
   94 CONTINUE
      NNB = NNB - 2
      GO TO 92
*
*       Add ICOMP & JCOMP to LISTR (maximum of MLR/2 - 1 pairs).
   96 IF (NNB.GT.MLR - 4) THEN
*       Note that NNB is one less than the actual membership.
          DO 98 K = 2,NNB
              LISTR(K) = LISTR(K+2)
   98     CONTINUE
          NNB = NNB - 2
*       Removal of the oldest KS pair.
      END IF
      LISTR(NNB+3) = ICOMP
      LISTR(NNB+4) = JCOMP
      LISTR(1) = NNB + 3
*
*       Copy flag index of disrupted pair (set in KSTERM).
      IFLAG = JLIST(3)
*       Add primordial pairs to LISTD (skip new KS pairs or primordials).
      IF (IFLAG.EQ.0.OR.IABS(JLIST(1) - JLIST(2)).EQ.1) GO TO 110
*
*       Check list of disrupted component names.
      NNB = LISTD(1) - 1
      KCOMP = 0
      DO 100 K = 2,NNB+1,2
          IF (LISTD(K).EQ.JLIST(1).AND.LISTD(K+1).EQ.JLIST(2)) KCOMP = 1
  100 CONTINUE
*
*       Include both components unless already members.
      IF (KCOMP.EQ.0) THEN
          IF (NNB.GT.MLD - 4) THEN
              DO 102 K = 2,NNB
                 LISTD(K) = LISTD(K+2)
  102         CONTINUE
              NNB = NNB - 2
          END IF
*       Add most recent names at the end (limit is MLD/2 - 1 pairs).
          LISTD(NNB+3) = JLIST(1)
          LISTD(NNB+4) = JLIST(2)
          LISTD(1) = NNB + 3
      END IF
      IF (rank.eq.0.and.IFLAG.NE.-1) THEN
          WRITE (8,104)  IPAIR, IFLAG, JLIST(1), JLIST(2)
  104     FORMAT (' LISTD INCONSISTENCY!!  IPAIR IFLAG NAMES ',2I5,2I8)
      END IF
*
*       Update list of high velocity particles containing c.m. members.
  110 NNB = LISTV(1)
      DO 130 L = 2,NNB+1
          IF (LISTV(L).EQ.ICM) THEN
*       Remove old c.m. and reduce the membership.
              DO 125 K = L,NNB
                  LISTV(K) = LISTV(K+1)
  125         CONTINUE
              LISTV(1) = LISTV(1) - 1
          END IF
*       Reduce higher particle locations by one.
          IF (LISTV(L).GT.ICM) THEN
              LISTV(L) = LISTV(L) - 1
          END IF
  130 CONTINUE
*
      RETURN
*
      END
