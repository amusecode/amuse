      SUBROUTINE IMF2(BODY10,BODYN)
*
*
*       Mass function for binaries & single stars.
*       ------------------------------------------
*       BODY10,BODYN are in M_sun.
*
      INCLUDE 'common6.h'
      INTEGER KCM(NMAX)
      REAL*8  RAN2,IMFBD
      REAL*8  LM,UM,BCM
      EXTERNAL IMFBD
      COMMON/WORK1/  BCM(NMAX)
*
*=========================  if KZ(20)=2:
* KTG3 MF with alpha1=1.3 :
      DATA  G1,G2,G3,G4  /0.19,1.55,0.050,0.6/
* KTG3 MF with alpha1=1.1 :
c      DATA  G1,G2,G3,G4  /0.28,1.14,0.010,0.1/
*=========================
*
*       Generate initial mass function (N-NBIN0 singles & 2*NBIN0 binaries).
      KDUM = IDUM1
      ZMASS = 0.0D0
      DO 10 I = 1,N+NBIN0
    5     XX = RAN2(KDUM)
*
*       Choose between Kroupa et al (M.N. 262, 545) & Eggleton (book).
          IF (KZ(20).EQ.2.OR.KZ(20).EQ.4) THEN
              ZM = 0.08 + (G1*XX**G2 + G3*XX**G4)/(1.0 - XX)**0.58
          ELSE IF (KZ(20).EQ.3.OR.KZ(20).EQ.5) THEN
              ZM = 0.3*XX/(1.0 - XX)**0.55
          ELSE IF (KZ(20).EQ.6) THEN
              LM = BODYN
              UM = BODY10
              ZM = IMFBD(XX,LM,UM)
          END IF
*       Include possibility of setting non-MS types. 
          KSTAR(I) = 1
*
*       See whether the mass falls within the specified range.
          IF (ZM.GE.BODYN.AND.ZM.LE.BODY10) THEN
              BODY(I) = ZM
              ZMASS = ZMASS + BODY(I)
          ELSE
              GO TO 5
          END IF
   10 CONTINUE
*
*       See whether to skip mass function for binaries.
      IF (NBIN0.EQ.0) GO TO 50
*
*       Merge binary components in temporary variable for sorting.
      DO 20 I = 1,NBIN0
          BCM(I) = BODY(2*I-1) + BODY(2*I)
          JLIST(I) = I
          KCM(2*I-1) = KSTAR(2*I-1)
          KCM(2*I) = KSTAR(2*I)
   20 CONTINUE
*
*       Sort total binary masses in increasing order.
      IF (NBIN0.GT.1) THEN
          CALL SORT1(NBIN0,BCM,JLIST)
      END IF
*
*       Save scaled binary masses in decreasing order for routine BINPOP.
      DO 30 I = 1,NBIN0
          JB = JLIST(NBIN0-I+1)
          BODY0(2*I-1) = MAX(BODY(2*JB-1),BODY(2*JB))/ZMASS
          BODY0(2*I) = MIN(BODY(2*JB-1),BODY(2*JB))/ZMASS
          IF (KZ(20).GT.3.AND.KZ(20).LE.5) THEN
*       Adopt correlation (m1/m2)' = (m1/m2)**0.4 & constant sum (Eggleton).
              ZMB = BODY0(2*I-1) + BODY0(2*I)
              RATIO = BODY0(2*I-1)/BODY0(2*I)
              BODY0(2*I) = ZMB/(1.0 + RATIO**0.4)
              BODY0(2*I-1) = ZMB - BODY0(2*I)
          END IF
          KSTAR(2*I-1) = KCM(2*JB-1)
          KSTAR(2*I) = KCM(2*JB)
   30 CONTINUE
*
*       Merge binary components into single stars for scaling purposes.
      ZMB = 0.0
      DO 40 I = 1,NBIN0
          BODY(NBIN0-I+1) = BCM(I)
          ZMB = ZMB + BCM(I)
   40 CONTINUE
*
      if(rank.eq.0)
     &WRITE (6,45)  NBIN0, BODY(1), BODY(NBIN0), ZMB/FLOAT(NBIN0)
   45 FORMAT (//,12X,'BINARY STAR IMF:    NB =',I6,
     &               '  RANGE =',1P,2E10.2,'  <MB> =',E9.2)
*
*       Move the single stars up to form compact array of N members.
   50 IF (N.LE.NBIN0) GO TO 90
      NS = 0
      DO 60 L = 1,N-NBIN0
          BODY(NBIN0+L) = BODY(2*NBIN0+L)
          NS = NS + 1
          BCM(NS) = BODY(NBIN0+L)
          KCM(NS) = KSTAR(2*NBIN0+L)
*         JLIST(NS) = NBIN0 + L
          JLIST(NS) = NS
   60 CONTINUE
*
*       Sort masses of single stars in increasing order.
      CALL SORT1(NS,BCM,JLIST)
*
*       Copy the masses of single stars to COMMON in decreasing order.
      ZMS = 0.0
      DO 70 I = 1,NS
          BODY(N-I+1) = BCM(I)
          ZMS = ZMS + BCM(I)
          KSTAR(N-I+1+NBIN0) = KCM(JLIST(I))
   70 CONTINUE
*
      if(rank.eq.0)
     &WRITE (6,80)  N-NBIN0, BODY(NBIN0+1), BODY(N), ZMS/FLOAT(N-NBIN0)
   80 FORMAT (/,12X,'SINGLE STAR IMF:    NS =',I6,'  RANGE =',1P,2E10.2,
     &                                            '  <MS> =',E9.2)
*
*       Replace input value by actual mean mass in solar units.
   90 ZMBAR = ZMASS/FLOAT(N)
*
      RETURN
*
      END
