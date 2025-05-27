      SUBROUTINE LEVELS
*
*
*       Diagnostics for block time-step levels.
*       ---------------------------------------
*       plus estimation of parallel speed-up (R.Sp.)
*
      INCLUDE 'common6.h'
      INTEGER JHIST,JHISTR,JHISTU
      COMMON/BLKLVL/JHIST(0:NMAX),JHISTR(0:NMAX),JHISTU(0:NMAX)
      INTEGER  IHIST(NMAX),IHISTR(NMAX),IHISTU(NMAX)
      INTEGER IPES,IPROC(9),IY(1024),IYR(1024)
      REAL*8 XSPEED(9),XSPEDR(9)
      DATA IPROC/4,8,16,32,64,128,256,512,1024/
*
*
*       Initialize histogram counters.
      DO 10 J = 1,NMAX
          IHIST(J) = 0
          IHISTR(J) = 0
          IHISTU(J) = 0
   10 CONTINUE
*
*       Loop over all single particles & c.m.
      JMAXI = 0
      JMAXR = 0
      FAC = 1.0/LOG(1.9999999)
      DO 20 I = IFIRST,NTOT
          IF (BODY(I).EQ.0.0D0) GO TO 20
          J = MAX(1,1 - INT(LOG(STEP(I))*FAC))
          IF(J.GT.NMAX)J = NMAX
          IHIST(J) = IHIST(J) + 1
          JMAXI = MAX(J,JMAXI)
          J = MAX(1,1 - INT(LOG(STEPR(I))*FAC))
          IF(J.GT.NMAX)J = NMAX
          IHISTR(J) = IHISTR(J) + 1
          JMAXR = MAX(J,JMAXR)
   20 CONTINUE
*
      JMAXU = 0
*       Loop over KS binaries
      DO 25 IPAIR = 1,NPAIRS
          I1 = 2*IPAIR - 1
          J = MAX(1,1 - INT(LOG(STEP(I1))*FAC))
          IF(J.GT.NMAX)J = NMAX
          IHISTU(J) = IHISTU(J) + 1
          JMAXU = MAX(J,JMAXU)
   25 CONTINUE
*
*       Print histograms of block-steps (STEPR with KZ(33) > 1).
      JMAXI=MIN(JMAXI,NMAX)
      JMAXR=MIN(JMAXR,NMAX)
      JMAXU=MIN(JMAXU,NMAX)
      if(rank.eq.0)then
	  WRITE (6,30)  (IHIST(J),J=1,JMAXI)
   30 FORMAT (' STEP I ',22I5,(/,24I5))
      IF (KZ(33).GT.1)WRITE (6,301)  (IHISTR(J),J=1,JMAXR)
  301 FORMAT (' STEP R ',50I5,(/,24I5))
      IF (KZ(8).GT.0 .OR. NBIN0.GT.0 )WRITE (6,3001)
     *     (IHISTU(J),J=1,JMAXU)
 3001 FORMAT (' STEP U ',50I5,(/,24I5))
      end if
*
*       IPROC contains a list of possible processor numbers,
*       presently   4 to 1024 in powers of 2.
*
      DO 35 NN = 1,9
      IPES = IPROC(NN)
      DO 36 I = 1,IPES
      IY(I) = 0
 36   IYR(I) = 0
*
*       JHIST contains accumulated frequencies from INTGRT for irr. steps. 
*       JHIST(80) = 5 e.g. means that it occurred 5 times since the last
*                   call of LEVELS that 80 force calculations had
*                   to be done in parallel
*       IY maps this information to a certain number of processors,
*       e.g. in the above case for 64 processors it would mean that
*       five times 64 processors could be fully used, and five times
*       only 16 processors (to get the remaining forces).
*       JHIST(80) = 5 results in incrementing IY(64) and IY(16) by 5.
*
*       (IYE, JHISTR do the same for the regular steps)
*
      DO 37 I = 1,NMAX
      IF(JHIST(I).EQ.0)GOTO 37
      J = MOD(I,IPES)
      K = I/IPES
      IF(J.GT.0)IY(J)=IY(J)+JHIST(I)
      IF(K.GT.0)IY(IPES)=IY(IPES)+K*JHIST(I)
 37   CONTINUE
*
      DO 371 I = 1,NMAX
      IF(JHISTR(I).EQ.0)GOTO 371
      J = MOD(I,IPES)
      K = I/IPES
      IF(J.GT.0)IYR(J)=IYR(J)+JHISTR(I)
      IF(K.GT.0)IYR(IPES)=IYR(IPES)+K*JHISTR(I)
 371  CONTINUE
*
      ISPEED = 0
      ISPEDR = 0
      ITOT = 0
      ITOTR = 0
      DO 38 I = 1,IPES
      ISPEED = ISPEED + I*IY(I)
      ISPEDR = ISPEDR + I*IYR(I)
      ITOT = ITOT + IY(I)
      ITOTR = ITOTR + IYR(I)
 38   CONTINUE
*
*       Estimate theoretical speedup by (sum I*IY) / (sum IY);
*       enumerator is the total number of force calculations necessary;
*       denominator is the number of force calculation cycles necessary
*       on a parallel machine. Separately done for regular/irregular
*       steps. Communication not yet included in any way.
*
      IF(ITOT.NE.0)THEN
      XSPEED(NN) = REAL(ISPEED)/REAL(ITOT)
      ELSE
      XSPEED(NN) = 0.D0
      END IF
      IF(ITOTR.NE.0)THEN
      XSPEDR(NN) = REAL(ISPEDR)/REAL(ITOTR)
      ELSE
      XSPEDR(NN) = 0.D0
      END IF
 35   CONTINUE
*
	  if(rank.eq.0)then
      WRITE (6,40)  (IPROC(J),XSPEED(J),J=1,9)
 40   FORMAT (' Max Speedup Irr: ',1P,9(I5,D9.2))
      WRITE (6,41)  (IPROC(J),XSPEDR(J),J=1,9)
 41   FORMAT (' Max Speedup Reg: ',1P,9(I5,D9.2))
      end if
*
c$$$      DO 50 J = 1,NMAX
c$$$          JHIST(J) = 0
c$$$          JHISTR(J) = 0
c$$$   50 CONTINUE
*
      RETURN
*
      END
