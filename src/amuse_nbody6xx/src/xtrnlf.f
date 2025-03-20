      SUBROUTINE XTRNLF(XI,XIDOT,FIRR,FREG,FD,FDR,KCASE)
*
*
*       External force & first derivative.
*       ----------------------------------
*
      INCLUDE 'common6.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2
      REAL*8  XI(3),XIDOT(3),FIRR(3),FREG(3),FD(3),FDR(3),
     &        XG(3),XGDOT(3),FM(3),FMD(3),FS(3),FSD(3)
*
*
*       See whether to include a linearized galactic tidal force (two cases).
      IF (KZ(14).LE.2) THEN
          FIRR(1) = FIRR(1) + TIDAL(4)*XIDOT(2)
          FIRR(2) = FIRR(2) - TIDAL(4)*XIDOT(1)
          FD(1) = FD(1) + TIDAL(4)*(FIRR(2) + FREG(2))
          FD(2) = FD(2) - TIDAL(4)*(FIRR(1) + FREG(1))
*       Add smooth part to the regular components (KCASE = 1 in REGINT).
          IF (KCASE.GT.0) THEN
              FREG(1) = FREG(1) + TIDAL(1)*XI(1)
              FREG(3) = FREG(3) + TIDAL(3)*XI(3)
              FDR(1) = FDR(1) + TIDAL(1)*XIDOT(1)
              FDR(3) = FDR(3) + TIDAL(3)*XIDOT(3)
          END IF
      END IF
*
*       Consider point-mass, disk and/or logarithmic halo model.
      IF (KZ(14).EQ.3.AND.KCASE.GT.0) THEN
          DO 5 K = 1,3
              XG(K) = RG(K) + XI(K)
              XGDOT(K) = VG(K) + XIDOT(K)
    5     CONTINUE
*       Employ differential instead of linearized forms for better accuracy.
          IF (GMG.GT.0.0D0) THEN
              CALL FNUC(RG,VG,FS,FSD)
              CALL FNUC(XG,XGDOT,FM,FMD)
              DO 10 K = 1,3
                  FREG(K) = FREG(K) + (FM(K) - FS(K))
                  FDR(K) = FDR(K) + (FMD(K) - FSD(K))
   10         CONTINUE
          END IF
*
*       Include Miyamoto disk for positive disk mass.
          IF (DISK.GT.0.0D0) THEN
              CALL FDISK(RG,VG,FS,FSD)
              CALL FDISK(XG,XGDOT,FM,FMD)
              DO 20 K = 1,3
                  FREG(K) = FREG(K) + (FM(K) - FS(K))
                  FDR(K) = FDR(K) + (FMD(K) - FSD(K))
   20         CONTINUE
          END IF
*
*       Check addition of logarithmic halo potential to regular force.
          IF (V02.GT.0.0D0) THEN
              CALL FHALO(RG,VG,FS,FSD)
              CALL FHALO(XG,XGDOT,FM,FMD)
              DO 30 K = 1,3
                  FREG(K) = FREG(K) + (FM(K) - FS(K))
                  FDR(K) = FDR(K) + (FMD(K) - FSD(K))
   30         CONTINUE
          END IF
      END IF
*
*       Include optional Plummer potential in the regular force.
      IF ((KZ(14).EQ.4.OR.MP.GT.0.0D0).AND.KCASE.GT.0) THEN
          RI2 = AP2
          RRDOT = 0.0
          DO 40 K = 1,3
              RI2 = RI2 + XI(K)**2
              RRDOT = RRDOT + XI(K)*XIDOT(K)
   40     CONTINUE
          IF (TIME + TOFF.GT.TDELAY) THEN
              ZMDOT = -MP0*MPDOT/(1.0 + MPDOT*(TIME+TOFF - TDELAY))**2
          ELSE
              ZMDOT = 0.0
          END IF
          ZF = 1.0/RI2
          ZF2 = ZF**1.5
*       Absorb scaling factor in 1/R3 term as MP*ZF2 (cf. Heggie & Hut p.73).
          FMP = MP*ZF2
          DO 50 K = 1,3
              FREG(K) = FREG(K) - XI(K)*FMP
              FDR(K) = FDR(K) - (XIDOT(K) - 3.0*RRDOT*ZF*XI(K))*FMP
              FDR(K) = FDR(K) - ZMDOT*ZF2*XI(K)
   50     CONTINUE
      END IF
*
      RETURN
*
      END
