      SUBROUTINE NEWTEV(NNB,IX)
*
*
*       New look-up time of hierarchy.
*       ------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      REAL*8  M0,M1,MG,TSCLS(20),LUMS(10),GB(10)
*
*
*       Initialize iteration counter for binary component.
      ITER = 0
*
*       Determine smallest look-up time for single stars.
    1 ITER = ITER + 1
      TM = 1.1D+10
      DO 5 L = 1,NNB
          I = JLIST(L)
*       Include safety check (just in case).
          IF (I.LE.0) GO TO 5
*       Replace any binary c.m. by primary component.
          IF (I.GT.N) THEN
              I = 2*(I - N) - 1
          END IF
          IF (TEV(I).LT.TM) THEN
              TM = TEV(I)
              IX = I
          END IF
    5 CONTINUE
*
*       Skip increase of smallest look-up time above 10 Myr.
      IF ((TEV(IX) - TIME)*TSTAR.GT.10.0) THEN
          GO TO 20
      END IF
*
*       Distinguish between ghost, single star or KS component.
      IF (BODY(IX).EQ.0.0D0) THEN
          CALL FINDM(IX,ITERM,MG)
          IF (ITERM.LT.0) GO TO 20
          M1 = MG*SMU
      ELSE IF (IX.GE.IFIRST) THEN
          M1 = BODY(IX)*SMU
      ELSE
          ITER = ITER + 1
          IPAIR = KVEC(IX)
          IF (NAME(N+IPAIR).GT.0) THEN
              M1 = BODY(IX)*SMU
          ELSE
*       Find the merger index by identifying original c.m. name.
              IM = 0
              DO 10 K = 1,NMERGE
                  IF (NAMEM(K).EQ.NAME(N+IPAIR)) IM = K
   10         CONTINUE
              IF (NAME(N+IPAIR).LT.-2*NZERO) THEN
                  WRITE (6,12)  NAME(N+IPAIR), NAME(I), IM
   12             FORMAT (' NEWTEV WARNING    NAMC NAMI IM ',2I7,I4)
              END IF
              IF (IM.EQ.0) GO TO 20
*       Choose relevant component for standard triple first.
              K = 1
              IF (IX.EQ.2*IPAIR) K = 2
              M1 = CM(K,IM)*SMU
*       Copy mass from K=3 or K=4 for merged quadruple system.
              IF (CM(3,IM).GT.0.0D0) THEN
                  K = K + 2
                  M1 = CM(K,IM)*SMU
*       Restore mass of second component if quad is excluded.
              ELSE IF (IX.EQ.2*IPAIR) THEN
                  M1 = BODY(IX)*SMU
              END IF
              IF (M1.EQ.0.0D0) THEN
                  M1 = CM(K,IM)*SMU
                  WRITE (6,15)  NAME(IX), K, M1
   15             FORMAT (' DANGER!    NEWTEV    NM K M1 ',2I6,F7.3)
              END IF
          END IF
      END IF
*
*       Obtain stellar evolution time-scales. 
      M0 = BODY0(IX)*SMU
      KW = KSTAR(IX)
      CALL star(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
*
*       Determine new time-scale for changes in radius & mass.
      AGE = TEV0(IX)*TSTAR - EPOCH(IX)
      CALL trdot2(KW,AGE,TM,TN,TSCLS,DTM,DTR)
      DT = DTM/TSTAR
*
*       Increase look-up time by at most a factor 2 if due in 10 Myr.
      DT = MIN(DT,10.0D0/TSTAR)
*       Impose limit to prevent multiple increases without updating.
      DTX = TEV0(IX) + 2.0*DT - TEV(IX)
      DT = MIN(DTX,DT)
      DT = MAX(DT,0.0D0)
      TEV(IX) = TEV(IX) + DT
*
*       Perform a second iteration for possible binary component.
      IF (ITER.LE.2) GO TO 1
*
   20 RETURN
*
      END
