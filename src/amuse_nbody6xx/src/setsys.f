      SUBROUTINE SETSYS
*
*
*       Selection of chain system.
*       --------------------------
*
      INCLUDE 'common6.h'
      PARAMETER  (NMX=10,NMX3=3*NMX,NMXm=NMX*(NMX-1)/2)
      REAL*8  M,MASS,MC,MIJ,MKK
      COMMON/CHAIN1/  XCH(NMX3),VCH(NMX3),M(NMX),
     &                ZZ(NMX3),WC(NMX3),MC(NMX),
     &                XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX),
     &                MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),NN
      COMMON/CHAINC/  XC(3,NCMAX),UC(3,NCMAX),BODYC(NCMAX),ICH,
     &                LISTC(LMAX)
      COMMON/CHREG/  TIMEC,TMAX,RMAXC,CM(10),NAMEC(6),NSTEP1,KZ27,KZ30
      COMMON/BINARY/  ZM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAG(MMAX)
      SAVE JSAVE
      INTEGER JSAVE(3)
*
*
*       Check whether new (or renewed) chain or addition of member(s).
      IF (NCH.GT.0) GO TO 10
*
*       Include treatment for common envelope stage (denoted by NCH < 0).
      IF (NCH.LT.0) THEN
          NCH = -NCH
*       Copy members from JLIST (set in CHTERM and modified in CMBODY).
          DO 1 L = 1,NCH
              J = JLIST(L)
              NAMEC(L) = NAME(J)
              BODYC(L) = BODY(J)
              M(L) = BODY(J)
    1     CONTINUE
          GO TO 50
      END IF
*
*       Initialize chain indices, names & masses for maximum membership.
      DO 5 L = 1,4
          JLIST(L) = 2*NPAIRS + L
          NAMEC(L) = NAME(2*NPAIRS+L)
          BODYC(L) = BODY(2*NPAIRS+L)
          M(L) = BODY(2*NPAIRS+L)
    5 CONTINUE
      CM(9) = EBCH0
*
*       Include treatment for near-synchronous binary as inert body (B-B).
      IF (JCLOSE.GT.N.AND.KZ(27).GT.0.AND.KZ(26).LT.2) THEN
          RSUM = RMIN
          NCH = 2
          GO TO 10
      END IF
*
*       Define chain membership for three-body or four-body case.
      IF (JCLOSE.LE.N.AND.JCLOSE.GT.0) THEN
          NCH = 3
          JLIST(3) = JCLOSE
          NAMEC(3) = NAME(JCLOSE)
          BODYC(3) = BODY(JCLOSE)
          M(3) = BODY(JCLOSE)
*       See whether a second single body should be added.
          IF (JCMAX.LE.N.AND.JCMAX.GE.IFIRST) THEN
              JLIST(4) = JCMAX
              NAMEC(4) = NAME(JCMAX)
              BODYC(4) = BODY(JCMAX)
              M(4) = BODY(JCMAX)
              NCH = 4
          END IF
      ELSE
          NCH = 4
      END IF
*
*       Check for addition of binary (NCH <= 4).
      IF (JCMAX.GT.N.AND.NCH.LE.4) THEN
          KSP2 = JCMAX - N
          IF (KSP2.GT.KSPAIR) KSP2 = KSP2 - 1
          KSPAIR = KSP2
          JCOMP = 0
*       Save current members to prevent over-writing in KSTERM.
          DO 6 L = 1,NCH
              JSAVE(L) = JLIST(L)
    6     CONTINUE
          CALL KSTERM
*       Note that second binary will now come first in N-body arrays.
          DO 7 L = 1,NCH
              JLIST(L) = JSAVE(L)
    7     CONTINUE
*       Add terminated KS components to chain arrays.
          DO 8 L = 1,2
              NCH = NCH + 1
              JLIST(NCH) = 2*NPAIRS + L
              NAMEC(NCH) = NAME(2*NPAIRS+L)
              BODYC(NCH) = BODY(2*NPAIRS+L)
              M(NCH) = BODY(2*NPAIRS+L)
    8     CONTINUE
      END IF
      GO TO 50
*
*       Improve coordinates & velocities of single perturber or c.m. body.
   10 CALL XVPRED(JCLOSE,-1)
*
*       Expand membership and save chain variables (single body or KS pair).
      IF (JCLOSE.LE.N) THEN
          NCH = NCH + 1
          JLIST(NCH) = JCLOSE
          NAMEC(NCH) = NAME(JCLOSE)
          BODYC(NCH) = BODY(JCLOSE)
          M(NCH) = BODY(JCLOSE)
      ELSE
          KSPAIR = JCLOSE - N
*       Check for synchronous tidal binary (save dormant energy in ECOLL).
          IF (KZ(27).GT.0.AND.KZ(26).LT.2) THEN
              SEMI = -0.5*BODY(JCLOSE)/H(KSPAIR)
              IF (SEMI.LT.0.01*RSUM) THEN
                  NCH = NCH + 1
                  JLIST(NCH) = JCLOSE
                  NAMEC(NCH) = NAME(JCLOSE)
                  BODYC(NCH) = BODY(JCLOSE)
                  M(NCH) = BODY(JCLOSE)
*       Define temporary KS ghost by saving masses and binding energy.
                  T0(2*KSPAIR-1) = 1.0E+06
                  LIST(1,2*KSPAIR-1) = 0
                  BODYC(9) = BODY(2*KSPAIR-1)
                  BODYC(10) = BODY(2*KSPAIR)
                  ZMU = BODY(2*KSPAIR-1)*BODY(2*KSPAIR)/BODY(JCLOSE)
                  ECOLL = ECOLL + ZMU*H(KSPAIR)
                  BODY(2*KSPAIR-1) = 0.0D0
                  BODY(2*KSPAIR) = 0.0D0
                  GO TO 50
              END IF
          END IF
*
*       Re-activate any merged binary before terminating as last pair.
          JG = 0
          IF (NAME(JCLOSE).LT.0) THEN
*       Identify merger index and ghost for addition to chain.
              DO 12 K = 1,NMERGE
                  IF (NAMEM(K).EQ.NAME(JCLOSE)) THEN
                      IM = K
                  END IF
   12         CONTINUE
*       Note ghost must be single for maximum chain membership of 6.
              DO 14 J = 1,N
                  IF (BODY(J).EQ.0.0D0.AND.NAME(J).EQ.NAMEG(IM)) THEN
                      JG = J
                  END IF
   14         CONTINUE
              if(rank.eq.0)
     &        WRITE (6,15)  NAME(JCLOSE), NAME(JG), RSUM,  R(JCLOSE-N)
   15         FORMAT (' SETSYS HIARCH    NM NMG RSUM RB ',
     &                                   I6,I5,1P,2E10.2)
              IPHASE = 7
              CALL RESET
              KSPAIR = NPAIRS
          END IF
*
*       Save global indices of existing members (KSTERM uses JLIST).
          DO 16 L = 1,NCH
              JPERT(L) = JLIST(L)
   16     CONTINUE
*
*       Add energy of absorbed binary to the current initial energy.
          EB = BODY(2*KSPAIR-1)*BODY(2*KSPAIR)*H(KSPAIR)/BODY(N+KSPAIR)
          CM(9) = CM(9) + EB
          EBCH0 = EBCH0 + EB
*
*       Check saving of c.m. type and names.
          IF (KSAVE(1).EQ.0.OR.KSTAR(N+KSPAIR).NE.0) THEN
              KSAVE(1) = KSTAR(N+KSPAIR)
              KSAVE(2) = NAME(2*KSPAIR-1) + NAME(2*KSPAIR)
          END IF
*
*       Terminate KS pair and copy components (JCOMP=0 excludes ghost).
          IPHASE = 8
          JCOMP = 0
          CALL KSTERM
*
*       Restore old members.
          DO 18 L = 1,NCH
              JLIST(L) = JPERT(L)
   18     CONTINUE
*
*       Add terminated KS components to chain arrays.
          DO 20 L = 1,2
              NCH = NCH + 1
              JLIST(NCH) = 2*NPAIRS + L
              NAMEC(NCH) = NAME(2*NPAIRS+L)
              BODYC(NCH) = BODY(2*NPAIRS+L)
              M(NCH) = BODY(2*NPAIRS+L)
   20     CONTINUE
*       See whether to include merger ghost.
          IF (JG.GT.0) THEN
              NCH = NCH + 1
              JLIST(NCH) = JG
              NAMEC(NCH) = NAME(JG)
              BODYC(NCH) = BODY(JG)
              M(NCH) = BODY(JG)
          END IF
          IF (NCH.GT.6) THEN
              if(rank.eq.0)
     &        WRITE (6,30)  NCH
   30         FORMAT (' DANGER!    NCH ',I4)
              STOP
          END IF
      END IF
*
*       Specify membership for chain COMMON.
   50 NN = NCH
*
      RETURN
*
      END
