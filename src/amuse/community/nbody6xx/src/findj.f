      SUBROUTINE FINDJ(I,J,IMERGE)
*
*
*       Find merger ghost.
*       ------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
*
*
*       Check whether the search index is a c.m. body or KS component.
      J = -1
      IMERGE = 0
      IF (I.GE.IFIRST) THEN
          ICM = I
          IF (I.LE.N.OR.NAME(ICM).GE.0) THEN
              if(rank.eq.0)
     &        WRITE (6,1)  I, NAME(I), IFIRST, N
    1         FORMAT (' WARNING!    FINDJ    I NAM I* N ',4I6)
              GO TO 30
          END IF
      ELSE
*       Set KS pair and c.m. index (exit if no merger).
          JP = KVEC(I)
          ICM = N + JP
          IF (NAME(ICM).GE.0) GO TO 30
      END IF
*
*       Locate current position in the merger table.
      IMERGE = 0
      DO 5 K = 1,NMERGE
          IF (NAMEM(K).EQ.NAME(ICM)) THEN
              IMERGE = K
          END IF
    5 CONTINUE
      IF (IMERGE.EQ.0) GO TO 30
*
*       Find index of ghost particle with specified name (search NTOT).
      DO 10 K = 1,NTOT
          IF (NAME(K).EQ.NAMEG(IMERGE)) THEN
              J = K
              GO TO 20
          END IF
   10 CONTINUE
*
   20 IF (J.GT.0.AND.BODY(J).GT.0.0D0) THEN
          if(rank.eq.0)
     &    WRITE (6,25)  I, JP, NAME(ICM), J, BODY(J)
   25     FORMAT (' WARNING!    FINDJ    I JP NAM J MJ',2I5,2I6,1P,E9.1)
          J = -2
      END IF
*
*       Assign nominal value of merger index if not identified.
   30 IF (IMERGE.EQ.0) IMERGE = 1
*
      RETURN
*
      END
