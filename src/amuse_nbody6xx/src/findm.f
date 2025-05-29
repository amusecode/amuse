      SUBROUTINE FINDM(I,ITERM,MG)
*
*
*       Find ghost mass.
*       ----------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      REAL*8 MG
*
*
*       Distinguish between KS component and single particle or c.m.
      ITERM = 0
      IF (I.LT.IFIRST) THEN
          IPAIR = KVEC(I)
          ICM = N + IPAIR
*       Consider standard binary, simple hierarchy or high-order system.
          IF (NAME(ICM).GT.0) THEN
              IM = 0
*       Obtain merger index for later and identify component of CM.
              DO 3 K = 1,NMERGE
                  IF (NAMEG(K).EQ.NAME(ICM)) IM = K
    3         CONTINUE
              IF (IM.EQ.0) THEN
                  ITERM = -1
                  GO TO 50
              END IF
              J1 = 2*IPAIR - 1
              IF (NAME(J1).EQ.NAME(I)) K1 = 3
              IF (NAME(J1+1).EQ.NAME(I)) K1 = 4
          ELSE IF (NAME(ICM).GE.-2*NZERO) THEN
              IM = 0
              DO 5 K = 1,NMERGE
                  IF (NAMEM(K).EQ.NAME(ICM)) IM = K
    5         CONTINUE
              IF (IM.EQ.0) THEN
                  ITERM = -1
                  GO TO 50
              END IF
              J1 = 2*IPAIR - 1
              K1 = 1
              IF (NAME(J1+1).EQ.NAME(I)) THEN
                  K1 = 2
              END IF
          ELSE IF (NAME(ICM).LT.-2*NZERO) THEN
              IM = 0
              DO 6 K = 1,NMERGE
                  IF (NAMEM(K).EQ.NAME(ICM)) IM = K
    6         CONTINUE
              IF (IM.EQ.0) THEN
                  ITERM = -1
                  GO TO 50
              END IF
              JH = 0
*       Search c.m. ghost name to get KS pair index.
              DO 8 J = N+1,NTOT
                  IF (NAME(J).EQ.NAMEG(IM)) JH = J
    8         CONTINUE
              IF (JH.EQ.0) THEN
                  ITERM = -1
                  GO TO 50
              END IF
              JPAIR = JH - N
              J1 = 2*JPAIR - 1
*       Compare component names in order to decide appropriate CM index.
              IF (NAME(J1).EQ.NAME(I)) THEN
                  K1 = 1
              ELSE IF (NAME(J1+1).EQ.NAME(I)) THEN
                  K1 = 2
              ELSE
                  K1 = 2
              END IF
          ELSE
              ITERM = -1
              GO TO 50
          END IF
      ELSE
*       Determine merger index of ghost particle.
          IM = 0
          DO 10 K = 1,NMERGE
              IF (NAMEG(K).EQ.NAME(I)) IM = K
   10     CONTINUE
          IF (IM.EQ.0) THEN
              ITERM = -1
              GO TO 50
          END IF
          IF (I.LE.N) THEN
              J1 = 0
*       Identify the location of the corresponding KS component.
              DO 15 J = 1,IFIRST
                  IF (NAME(J).EQ.-NAMEM(IM)) J1 = J
   15         CONTINUE
              IF (J1.EQ.0) THEN
                  ITERM = -1
                  GO TO 50
              END IF
              IPAIR = KVEC(J1)
*       Decide the relevant component by comparing look-up times.
              K1 = 1
              IF (TEV(2*IPAIR).LT.TEV(I)) K1 = 2
          ELSE
              ICM = 0
              DO 20 K = N+1,NTOT
                  IF (NAMEM(IM).EQ.NAME(K)) ICM = K
   20         CONTINUE
              IF (ICM.EQ.0) THEN
                  ITERM = -1
                  GO TO 50
              END IF
              IPAIR = ICM - N
              J = I
              K1 = 1
              IF (TEV(2*IPAIR).LT.TEV(J)) K1 = 2
          END IF
      END IF
*
*       Copy the desired mass from merger table.
      MG = CM(K1,IM)
*
   50 RETURN
*
      END
