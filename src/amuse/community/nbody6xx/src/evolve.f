      SUBROUTINE EVOLVE(IPAIR,LEVEL)
*
*
*       Binary diagnostics.
*       -------------------
*
      INCLUDE 'common6.h'
      COMMON/GAMDOT/  DGAM
      REAL*8  XCM(3),VCM(3),WORK(10)
      INTEGER  IFLAG(KMAX),IOUT(LMAX)
*
*
      ZKT = -BE(1)/FLOAT(NZERO)
      M = 0
      IF (LEVEL.EQ.0) GO TO 500
*     NBSTAT = NBSTAT + 1
*       Prepare diagnostic output for regularized binary.
*
      DO 10 IP = 1,NPAIRS
          IFLAG(IP) = 0
   10 CONTINUE
      JPAIR = IPAIR
*
   20 CALL RESOLV(JPAIR,1)
      ICM = N + JPAIR
      J2 = 2*JPAIR
      J1 = J2 - 1
      SEMI = -0.5*BODY(ICM)/H(JPAIR)
      ZN = SQRT(BODY(ICM)/ABS(SEMI)**3)
      RDOT = (X(1,J1) - X(1,J2))*(XDOT(1,J1) - XDOT(1,J2))
     &     + (X(2,J1) - X(2,J2))*(XDOT(2,J1) - XDOT(2,J2))
     &     + (X(3,J1) - X(3,J2))*(XDOT(3,J1) - XDOT(3,J2))
      ECC = (1.0 - R(JPAIR)/SEMI)**2 + RDOT**2/(SEMI*BODY(ICM))
      ECC = SQRT(ECC)
      APO = ABS(SEMI)*(1.0 + ECC)
      RI = SQRT((X(1,ICM) - RDENS(1))**2 + (X(2,ICM) - RDENS(2))**2 +
     &                                     (X(3,ICM) - RDENS(3))**2)
      ZK = 0.5*BODY(ICM)*(XDOT(1,ICM)**2+XDOT(2,ICM)**2+XDOT(3,ICM)**2)
      ZK = ZK/ZKT
      EB = BODY(J1)*BODY(J2)*H(JPAIR)/BODY(ICM)
      EB = -EB/ZKT
      TK = TWOPI/ZN
      NK = TCR0/TK
*
*       Print time interval since last binary output and set new time.
      DT = 1.0E+06
      K = KZ(4)
      DO 40 L = 1,K
          IF (TIME - TLASTB(L).LT.DT) DT = TIME - TLASTB(L)
   40 CONTINUE
      TLASTB(LEVEL) = TIME
*
      IF (LEVEL.LT.0) DGAM = 1.
      GS = SIGN(GAMMA(JPAIR),DGAM)
      IF (NAME(ICM).GT.0) THEN
          if(rank.eq.0)
     &    WRITE (4,50)  LEVEL, TIME, NAME(J1), NAME(J2), EB, SEMI, ECC,
     &                  NK, R(JPAIR), GS, LIST(1,J1), RI, ZK, DT, TK, M
      ELSE
          if(rank.eq.0)
     &    WRITE (4,50)  LEVEL, TIME, NZERO + NAME(J1), NAME(J2), EB,
     &                  SEMI, ECC, NK, R(JPAIR), GS, LIST(1,J1), RI, ZK,
     &                  DT, TK, M, NAME(J1), -(NAME(ICM) + NZERO)
      END IF
*
   50 FORMAT (I3,F8.2,I6,I5,F7.1,F9.5,F6.2,I7,F9.5,F6.2,
     &                                         I5,F7.2,F6.1,2F10.6,5I5)
      M = M + 1
      IFLAG(JPAIR) = -1
*
      NNB = LIST(1,J1)
      DO 60 L = 1,NNB
          IOUT(L) = 0
   60 CONTINUE
      RJMIN2 = 1.E+06
*
      NNB1 = NNB + 1
      DO 150 L = 2,NNB1
          J = LIST(L,J1)
          RIJ2 = (X(1,ICM) - X(1,J))**2 + (X(2,ICM) - X(2,J))**2
     &                                  + (X(3,ICM) - X(3,J))**2
          VIJ2 = (XDOT(1,ICM) - XDOT(1,J))**2 +
     &           (XDOT(2,ICM) - XDOT(2,J))**2 +
     &           (XDOT(3,ICM) - XDOT(3,J))**2
          BB = BODY(ICM) + BODY(J)
          BBI = 1.0/BB
          RIJ = SQRT(RIJ2)
          EB = 0.5*VIJ2 - BB/RIJ
          SEMIJ = -0.5*BB/EB
          RDOT = (X(1,ICM) - X(1,J))*(XDOT(1,ICM) - XDOT(1,J)) +
     &           (X(2,ICM) - X(2,J))*(XDOT(2,ICM) - XDOT(2,J)) +
     &           (X(3,ICM) - X(3,J))*(XDOT(3,ICM) - XDOT(3,J))
          ECCJ = (1.0 - RIJ/SEMIJ)**2 + RDOT**2/(SEMIJ*BB)
          ECCJ = SQRT(ECCJ)
          P = SEMIJ*(1.0 - ECCJ)
          RDOT = RDOT/RIJ
          EB = -BODY(J)*BODY(ICM)*EB*BBI/ZKT
          GEFF = 2.0*BODY(J)*(ABS(SEMI)/RIJ)**3/BODY(ICM)
          IF (GEFF.GT.99.9) GEFF = 99.9
*       Exclude apocentres > A with A = 2/N for hard binary.
          APO = MIN(APO,2.0D0/FLOAT(N))
*
*       Print small impact parameter & significant effective perturbation.
          IF (P.LT.3.*APO.AND.GEFF.GT.GPRINT(1)) THEN
              if(rank.eq.0)
     &        WRITE (4,100)  NAME(J), EB, SEMIJ, ECCJ, RIJ, GEFF, P,RDOT
  100         FORMAT (16X,I6,F7.1,F9.5,F6.2,F16.5,F6.2,F9.5,F7.2)
              IOUT(L-1) = 1
              IF (J.GT.N.AND.IFLAG(J-N).EQ.0) IFLAG(J-N) = 1
          END IF
*
*       Save parameters for the nearest perturber.
          IF (RIJ2.LT.RJMIN2) THEN
              RJMIN2 = RIJ2
              LMIN = L
              WORK(1) = EB
              WORK(2) = SEMIJ
              WORK(3) = ECCJ
              WORK(4) = RIJ
              WORK(5) = GEFF
              WORK(6) = P
              WORK(7) = RDOT
          END IF
  150 CONTINUE
*
      IF (IOUT(LMIN-1).EQ.0) THEN
          JMIN = LIST(LMIN,J1)
          if(rank.eq.0)
     &    WRITE (4,100)  NAME(JMIN), (WORK(K),K=1,7)
          IF (JMIN.GT.N.AND.IFLAG(JMIN-N).EQ.0) IFLAG(JMIN-N) = 1
      END IF
*
      DO 190 IP = 1,NPAIRS
          IF(IFLAG(IP).GT.0) THEN
              JPAIR = IP
*       Only include full output for active binaries.
              IF (GAMMA(JPAIR).GT.GMAX) GO TO 20
          END IF
  190 CONTINUE
*
*       Set output times to facilitate DTCRIT estimate in routine KSINT.
      KK = KZ(4)
      DO 200 K = 1,KK
          TLASTB(K) = TIME
  200 CONTINUE
      GO TO 900
*
*       Search for binaries of single particles & c.m. bodies.
  500 SIMAX = 0.01*TCR0
      IF (IPAIR.EQ.0) THEN
          ILOW = IFIRST
          IHIGH = NTOT
      ELSE
          ILOW = IPAIR
          IHIGH = IPAIR
      END IF
      NAM = IPAIR
      IF (NAM.GT.0) NAM = NAME(NAM)
*
      DO 600 I = ILOW,IHIGH
          IF (STEP(I).GT.SIMAX) GO TO 600
          J2 = 0
          RIJ = 1.0E+06
          NNB1 = LIST(1,I) + 1
*
          DO 250 L = 2,NNB1
              J = LIST(L,I)
              IF (STEP(J).GT.SIMAX) GO TO 250
              RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                      (X(3,I) - X(3,J))**2
              IF (RIJ2.GT.RIJ) GO TO 250
              RIJ = RIJ2
              J2 = J
  250     CONTINUE
*
          IF (J2.LT.I.AND.IPAIR.EQ.0) GO TO 600
          RIJ = SQRT(RIJ)
          VIJ2 = (XDOT(1,I) - XDOT(1,J2))**2 +
     &           (XDOT(2,I) - XDOT(2,J2))**2 +
     &           (XDOT(3,I) - XDOT(3,J2))**2
          BODYIJ = BODY(I) + BODY(J2)
          BBI = 1.0/BODYIJ
          EREL = 0.5*VIJ2 - BODYIJ/RIJ
          IF (EREL.GT.-0.1*ECLOSE) GO TO 600
          SEMI = -0.5*BODYIJ/EREL
          ZN = SQRT(BODYIJ/SEMI**3)
          RDOT = (X(1,I) - X(1,J2))*(XDOT(1,I) - XDOT(1,J2)) +
     &           (X(2,I) - X(2,J2))*(XDOT(2,I) - XDOT(2,J2)) +
     &           (X(3,I) - X(3,J2))*(XDOT(3,I) - XDOT(3,J2))
          ECC = (1.0 - RIJ/SEMI)**2 + RDOT**2/(SEMI*BODYIJ)
          ECC = SQRT(ECC)
*
          DO 520 K = 1,3
              XCM(K) = (BODY(I)*X(K,I) + BODY(J2)*X(K,J2))*BBI
              VCM(K) = (BODY(I)*XDOT(K,I) + BODY(J2)*XDOT(K,J2))*BBI
  520     CONTINUE
*
          RI = SQRT((XCM(1) - RDENS(1))**2 + (XCM(2) - RDENS(2))**2 +
     &                                       (XCM(3) - RDENS(3))**2)
          ZK = 0.5*BODYIJ*(VCM(1)**2 + VCM(2)**2 + VCM(3)**2)
          ZK = ZK/ZKT
          EB = BODY(I)*BODY(J2)*EREL*BBI
          EB = -EB/ZKT
          TK = TWOPI/ZN
          NK = TCR0/TK
          GB = 0.0
*
*       Estimate the relative perturbation at apocentre.
          DT = T0(I) - T0(J2)
          DO 540 K = 1,3
              GB = GB + (F(K,I) - (F(K,J2) + 3.0*FDOT(K,J2)*DT))**2
  540     CONTINUE
          GB = 2.0*SQRT(GB)*RIJ**2*BBI
          GB = GB*(SEMI*(1.0 + ECC)/RIJ)**3
          IF (GB.GT.99.9) GB = 99.9
*
*       Set time interval since last binary output.
          DT = 1.0E+06
          K = KZ(4)
          DO 550 L = 1,K
              IF (TIME - TLASTB(L).LT.DT) DT = TIME - TLASTB(L)
  550     CONTINUE
*
*       Only output soft binaries for dominant motion at apocentre.
          IF (EB.LT.1.0.AND.GB.GT.0.5) GO TO 600
*
          IF (I.LE.N.AND.J2.LE.N) THEN
              if(rank.eq.0)
     &        WRITE (4,50)  LEVEL, TIME, NAME(I), NAME(J2), EB, SEMI,
     &                      ECC, NK, RIJ, GB, LIST(1,I), RI, ZK, DT, TK,
     &                      NAM
          ELSE IF (I.LE.N.AND.J2.GT.N) THEN
              JJ2 = 2*(J2-N)
              JJ1 = JJ2-1
              if(rank.eq.0)
     &        WRITE (4,50)  LEVEL, TIME, NAME(I), NAME(J2), EB, SEMI,
     &                      ECC, NK, RIJ, GB, LIST(1,I), RI, ZK, DT, TK,
     &                      NAM, NAME(JJ1), NAME(JJ2)
          ELSE IF (I.GT.N.AND.J2.LE.N) THEN
              II2 = 2*(I-N)
              II1 = II2-1
              if(rank.eq.0)
     &        WRITE (4,50)  LEVEL, TIME, NAME(I), NAME(J2), EB, SEMI,
     &                      ECC, NK, RIJ, GB, LIST(1,I), RI, ZK, DT, TK,
     &                      NAM, NAME(II1), NAME(II2)
          ELSE
              II2 = 2*(I-N)
              II1 = II2-1
              JJ2 = 2*(J2-N)
              JJ1 = JJ2-1
              if(rank.eq.0)
     &        WRITE (4,50)  LEVEL, TIME, NAME(I), NAME(J2), EB, SEMI,
     &                      ECC, NK, RIJ, GB, LIST(1,I), RI, ZK, DT, TK,
     &                      NAM, NAME(II1), NAME(II2), NAME(JJ1),
     &                      NAME(JJ2)
          END IF
  600 CONTINUE
*
      IF (IPAIR.EQ.0) THEN
          TLASTS = TIME
      ELSE
          TLASTT = TIME
      END IF
*
  900 RETURN
*
      END
