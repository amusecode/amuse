      SUBROUTINE BINEV(IPAIR)
*
*
*       Binary evolution.
*       -----------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      REAL*8  M1,M2
      LOGICAL  FIRST
      SAVE  FIRST,IPREV,ISKIP
      DATA  FIRST,IPREV,ISKIP  /.TRUE.,0,0/
*
*
*       Open unit #17 the first time.
      IF (rank.eq.0.and.FIRST) THEN
          OPEN (UNIT=17,STATUS='UNKNOWN',FORM='FORMATTED',FILE='BINEV')
          FIRST = .FALSE.
          WRITE (17,5)
    5     FORMAT ('   TPHYS  NAM1  NAM2  K1  K2  KC   M1   M2     R1',
     &            '    R2   RI    ECC    SEMI   PERIOD  C',/)
      END IF
*
*       Define KS & c.m. indices and set component masses & radii.
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      I = N + IPAIR
      M1 = BODY(I1)*SMU
      M2 = BODY(I2)*SMU
      R1 = RADIUS(I1)*SU
      R2 = RADIUS(I2)*SU
*
*       Skip identical or alternating sequences (i.e. KSTAR(I) = -1 & 0).
      ISUM = KSTAR(I1) + KSTAR(I2) + KSTAR(I)
      IF (ISUM.EQ.IPREV.OR.ISUM.EQ.ISKIP) GO TO 20
      ISKIP = IPREV
      IPREV = ISUM
*
*       Form basic two-body elements (distinguish mergers & ghosts).
      IF (NAME(I).GT.0.AND.BODY(I).GT.0.0D0) THEN
          SEMI = -0.5*BODY(I)/H(IPAIR)
          BODYC = BODY(I)
          ECC2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODYC*SEMI)
      ELSE
          CALL FINDJ(I1,J,IM)
          IF (J.LE.0) RETURN
          BODYC = CM(1,IM) + CM(2,IM)
          SEMI = -0.5*BODYC/HM(IM)
          RJ = UM(1,IM)**2 + UM(2,IM)**2 + UM(3,IM)**2 + UM(4,IM)**2
          ECC2 = (1.0 - RJ/SEMI)**2
      END IF
*
      ECC = SQRT(ECC2)
      P = DAYS*SEMI*SQRT(ABS(SEMI)/BODYC)
      P = MIN(P,99999.9D0)
      P = MAX(P,-1.0D0)
      A0 = MIN(SEMI*SU,9999.9D0)
      A0 = MAX(A0,-1.0D0)
      RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                               (X(3,I) - RDENS(3))**2
      RI = MIN(SQRT(RI2),99.9D0)
      R1 = MIN(R1,999.9D0)
      R2 = MIN(R2,999.9D0)
*
*       Print one line for every new stage (type change of #I1, #I2 or c.m.).
      if(rank.eq.0)then
      WRITE (17,10)  TPHYS, NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2),
     &               KSTAR(I), M1, M2, R1, R2, RI, ECC, A0, P, IQCOLL
   10 FORMAT (F8.1,2I6,3I4,2F5.1,F7.1,F6.1,F5.1,F7.3,F8.1,F9.1,I3)
      CALL FLUSH(17)
      end if
*
   20 RETURN
*
      END
