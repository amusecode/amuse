      SUBROUTINE BINOUT
*
*
*       Binary analysis & output.
*       -------------------------
*
      INCLUDE 'common6.h'
      COMMON/ECHAIN/  ECH
      INTEGER  JD(15),JEB(15),JE(10)
*
*
*       Define semi-major axis & binding energy of hard binary (= KT).
      A0 = 1.0/FLOAT(NZERO)
      EB0 = -0.25/FLOAT(NZERO)
*       Initialize counters & variables.
      DO 10 J = 1,15
          JD(J) = 0
          JEB(J) = 0
   10 CONTINUE
      DO 20 J = 1,10
          JE(J) = 0
   20 CONTINUE
      EBMAX = 0.0
      DISP = 0.0
      EMAX = 0.0
      JOR = 0
      JEX = 0
      JC = 0
      JLAST = 0
      KLAST = 0
      E(1) = 0.0D0
      E(2) = 0.0D0
      NPOP(1) = 0
      NPOP(2) = 0
      NPOP(3) = N - 2*NPAIRS
      NEWB = 0
      K10 = 0
*
*       Obtain relevant distributions of KS binaries.
      DO 50 J = 1,NPAIRS
          I = N + J
          IF (H(J).GE.0.0.OR.BODY(I).EQ.0.0D0) GO TO 50
          J1 = 2*J - 1
          J2 = 2*J
*
*       Count original & exchanged binaries and core members.
          IF (IABS(NAME(J1) - NAME(J2)).EQ.1) THEN
              IF (NAME(J1).LE.2*NBIN0) JOR = JOR + 1
          ELSE
              IF (MIN(NAME(J1),NAME(J2)).LE.2*NBIN0) JEX = JEX + 1
          END IF
          RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                                   (X(3,I) - RDENS(3))**2
          IF (RI2.LT.RC2) JC = JC + 1
*
*       Adopt logarithmic distribution of semi-major axis (factor 2).
          SEMI = -0.5D0*BODY(I)/H(J)
          IF (SEMI.LT.A0) THEN
              K = 2 + LOG10(A0/SEMI)/LOG10(2.0)
              K = MIN(K,10)
              JLAST = MAX(JLAST,K)
              JD(K) = JD(K) + 1
          ELSE
              JD(1) = JD(1) + 1
          END IF
*
*       Form eccentricity dispersion & histogram.
          ECC2 = (1.0 - R(J)/SEMI)**2 + TDOT2(J)**2/(BODY(I)*SEMI)
          DISP = DISP + ECC2
          IF (ECC2.GT.EMAX) EMAX = ECC2
          K = 1 + 10.0*SQRT(ECC2)
          IF (K.LE.10) JE(K) = JE(K) + 1
*
*       Set up logarithmic distribution of binding energy (factor 2).
          EB = -0.5*BODY(J1)*BODY(J2)/SEMI
          IF (EB.LT.EB0) THEN
              K = 2 + LOG10(EB/EB0)/LOG10(2.0)
              K = MIN(K,14)
              KLAST = MAX(KLAST,K)
              JEB(K) = JEB(K) + 1
              EBMAX = MIN(EBMAX,EB)
          ELSE
              JEB(1) = JEB(1) + 1
          END IF
*
*       Define flag to distinguish primordial or non-primordial binary.
          IP = 1
          IF (LIST(2,J2).EQ.0) IP = 2
          IF (IP.EQ.2) NEWB = NEWB + 1
*
*       Sum the respective energies & populations.
          E(IP) = E(IP) + EB
          NPOP(IP) = NPOP(IP) + 1
*
*       Perform consistency check in case of over-writing.
          IF (LIST(2,J2).NE.-1.AND.LIST(2,J2).NE.0) THEN
              if(rank.eq.0)
     &        WRITE (6,35)  J, LIST(1,J1), LIST(1,J2), LIST(2,J2),
     &                      NAME(N+J)
   35         FORMAT (/,5X,'WARNING!   FLAG PROBLEM   PAIR NP NP2 FLAG',
     &                                                     ' NAME ',5I6)
          END IF
*
*       Produce special diagnostics for new binaries (unit 18).
          IF (IP.EQ.2.AND.H(J).LT.0.0) THEN
              VR = ((X(1,I) - RDENS(1))*XDOT(1,I) +
     &              (X(2,I) - RDENS(2))*XDOT(2,I) +
     &              (X(3,I) - RDENS(3))*XDOT(3,I))/SQRT(RI2)
              K = 0
              IF (NAME(I).LT.0) K = -1
              GX = GAMMA(J)*(SEMI*(1.0 + SQRT(ECC2))/R(J))**3
              if(rank.eq.0)
     &        WRITE (18,40)  TTOT, NAME(J1), NAME(J2), LIST(2,J2), K,
     &                       BODY(J1), BODY(J2), EB, SEMI, SQRT(ECC2),
     &                       GX, SQRT(RI2), VR
   40         FORMAT (' T =',F7.1,'  NAME = ',2I6,2I3,'  M =',2F9.4,
     &                       '  EB =',F10.5,'  A =',F8.5,'  E =',F5.2,
     &                      '  GX =',F6.3,'  RI =',F6.2,'  VR =',F4.1)
              CALL FLUSH(18)
          END IF
*       Reset c.m. Roche flag to standard type for non-zero eccentricity.
          IF (KSTAR(I).GE.10.AND.ECC2.GT.4.0D-06) KSTAR(I) = 0
          IF (KSTAR(I).GE.10) K10 = K10 + 1
   50 CONTINUE
*
*       Set fractional gain of binding energy (initial energy = -0.25).
      IF (TIME.LE.0.0D0) THEN
          JOR = NBIN0
          DB = 0.0
      ELSE
          E(9) = EMERGE
          DB = -4.0*(E(1) + E(2) + E(5) + E(7) + E(9) + ECOLL + EMDOT
     &                                                + EKICK - EBIN0)
      END IF
*
*       Print eccentricity & energy distributions.
      IF (NPAIRS.EQ.0) GO TO 70 
      if(rank.eq.0)
     &WRITE (6,60)  JOR, JEX, DB, SBCOLL, BBCOLL, CHCOLL, JC, NCHAOS,
     &              NEWB, K10, (JD(J),J=1,JLAST)
   60 FORMAT (/,' BINARIES',4X,'OR =',I5,'  EX =',I3,'  DB =',F7.3,
     &          '  SB =',F8.4,'  BB =',F8.4,'  CH =',F8.4,'  NC =',I3,
     &         '  NCH =',I4,'  NEWB =',I4,'  CIRC =',I4,'  N(A) =',10I4)
*
      IF (DISP.GT.0.0D0) DISP = SQRT(DISP/FLOAT(NPAIRS))
      EMAX = SQRT(EMAX)
      if(rank.eq.0)
     &WRITE (6,65)  DISP, EMAX, (NPOP(J),J=1,8), (JEB(K),K=1,KLAST)
   65 FORMAT (' <E> =',F5.2,'  EMAX =',F7.4,'  NPOP =',I5,I3,2I6,I4,3I3,
     &                                                 '  EB/KT =',14I4)
*
*       Form the basic internal energy (binaries & single particles).
   70 ETOT = 0.0D0
      E(3) = ZKIN - POT + ETIDE
      DO 80 J = 1,3
          ETOT = ETOT + E(J)
   80 CONTINUE
*
*       Sum all separate energy components to monitor conservation.
      ETOT = ETOT + ESUB + EMERGE + EMDOT + ECDOT + ECOLL
*
*       Include energy of chain if NCH > 0.
      IF (NCH.GT.0) THEN
          ETOT = ETOT + ECH
      END IF
*
*       Form the net energy gain in binary interactions (also in events.f).
      DEGRAV = EBIN + ESUB + EBESC + EMESC + EMERGE + EGRAV - EBIN0
      if(rank.eq.0) WRITE (6,90)  (E(J),J=1,10), ETOT, DETOT, DEGRAV
   90 FORMAT (' ENERGIES   ',10F10.5,'  ETOT =',F11.6,'  DETOT =',F10.6,
     &                                                '  DEGRAV =',F7.3)
*
      RETURN
*
      END
