      SUBROUTINE DECIDE(IPAIR,SEMI,ECC,EMAX,EMIN,TC,TG,EDAV,IQ)
*
*
*       Merger decision.
*       ----------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
*
*
*       Initialize the merger skip indicator and define KS indices.
      IQ = 0
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
*       Prevent over-writing in NMERGE+1.
      IF(NMERGE.GE.MMAX-1)THEN
         IQ = 1
         GOTO 40
      ENDIF
*
*       Accept active spiral but ensure some updating to increase pericentre.
      DMR = 0.D0
      IF (KSTAR(N+IPAIR).EQ.-2) THEN
          PMIN = SEMI*(1.0 - ECC)
          ICIRC = -1
          CALL TCIRC(PMIN,ECC,I1,I2,ICIRC,TC1)
*       Restrict circularization time to account for stellar evolution.
          TC1 = MIN(TC1,500.0D0)
*       Adopt termination interval 1/2 of predicted t_{circ} (at const R*).
          DT1 = 0.5*(TC1 + 0.1)/TSTAR
          IF (ECC.LT.0.1.AND.MAX(KSTAR(I1),KSTAR(I2)).LE.1) THEN
              DT1 = 10.0*DT1
          END IF
*       Rectify SPIRAL just in case.
          IF(LIST(1,I1).GT.0)THEN
*            IQ = 1
*            GO TO 30
*       Accept perturbed spiral case for now!
             DMR = -1.D0
             CALL CHRECT(IPAIR,DMR)
          ELSE
             CALL CHRECT(IPAIR,DMR)
          ENDIF
          GO TO 30
      END IF
*
*       Obtain eccentricity derivative due to #JCOMP (KS resolved in INDUCE).
      CALL EDOT(IPAIR,JCOMP,SEMI,ECC,ECCDOT)
*
*       Define c.m. index and set maximum stellar radius and peri.
      I = N + IPAIR
      RM = MAX(RADIUS(I1),RADIUS(I2))
      PM = SEMI*(1.0 - ECC)/RM
      PM = MIN(PM,99.9D0)
*
*     WRITE (75,1)  KSTAR(I1), KSTAR(I2), KSTAR(I), LIST(1,I1), NAME(I),
*    &              TTOT, ECC, EMIN, EMAX, TG, EDAV, PM
*   1 FORMAT (' DECIDE:    K* NP NM T E EM EX TG ED PM ',
*    &                     4I4,I6,F10.3,F8.4,2F6.2,2F7.2,F6.1)
*     CALL FLUSH(75)
*
*       See whether to deform binary orbit (skip chaos case).
      IF (KSTAR(I).NE.-1.AND.TG.LT.-10.0) THEN
*       Change eccentricity slowly by a fraction 0.1*(1 - E) of decay time.
          ECC0 = ECC
          DT1 = 0.1*(1.0 - ECC)*TG/TSTAR
          ECC = ECC + EDAV*DT1
          ECC = MAX(ECC,EMIN)
*       Impose temporary limit of 0.99*EMAX to allow for apsidal motion.
          ECC = MIN(ECC,0.99*EMAX)
*
*       Evaluate t_{circ} and (de/dt)_{circ} except for near-collision.
          PMIN = SEMI*(1.0 - ECC)
          IF (PMIN.GT.RM) THEN
              ICIRC = -1
              CALL TCIRC(PMIN,ECC,I1,I2,ICIRC,TC1)
*             CALL EDOT(IPAIR,JCOMP,SEMI,ECC,ECCDOT)
              CALL ECIRC(PMIN,ECC,I1,I2,ICIRC,TG,TC2,ECC2,EDT)
*             WRITE (6,10)  ECCDOT, R(IPAIR), SEMI, TDOT2(IPAIR)
*  10         FORMAT (' ECIRC:   ECDOT R A TD2 ',1P,4E10.2)
          ELSE
              TC1 = TG
              TC2 = TC1
              EDT = 2.0*EDAV
          END IF
*
*       Compare eccentricity derivatives and check t_{circ} (standard case).
          IF (ABS(EDT).GT.ABS(EDAV).AND.TC1.LT.50.0) THEN
*       Enforce new spiral on short circularization time (IQ > 0: no merger).
              IF (KSTAR(I).GE.0) IQ = 1
              ECC = ECC0
              IF (KSTAR(I).EQ.-2) GO TO 40
              GO TO 30
          END IF
*
*       Rectify spiral orbit before changing eccentricity.
          IF (KSTAR(I).EQ.-2) THEN
              CALL CHRECT(IPAIR,DMR)
          END IF
*
*       Transform to exact apocentre unless done already.
          IF (ABS(TDOT2(IPAIR)).GT.1.0D-12) THEN
              IF (R(IPAIR).GT.SEMI.AND.TDOT2(IPAIR).LT.0.0) THEN
                  CALL KSAPO(IPAIR)
              END IF
              CALL KSPERI(IPAIR)
              CALL KSAPO(IPAIR)
*       Check spiral case already at peri.
          ELSE IF (KSTAR(I).EQ.-2.AND.R(IPAIR).LT.SEMI.AND.
     &             ABS(TDOT2(IPAIR)).LT.1.0D-12) THEN
              CALL KSAPO(IPAIR)
          END IF
*
*       Deform orbit (H = const at apo), rectify and transform to X & XDOT.
          CALL DEFORM(IPAIR,ECC0,ECC)
          CALL KSRECT(IPAIR)
          CALL RESOLV(IPAIR,1)
*
*       Rectify spiral orbit (terminate on collision).
          IF (KSTAR(I).EQ.-2) THEN
              CALL CHRECT(IPAIR,DMR)
              IF (IPHASE.LT.0) THEN
                  IQ = 2
                  GO TO 40
              END IF
          END IF
*
*       Re-initialize KS polynomials for perturbed motion (also in DEFORM).
          IF (LIST(1,I1).GT.0) THEN
              IF (R(IPAIR).LT.SEMI.AND.
     &            ABS(TDOT2(IPAIR)).LT.1.0D-12) THEN
                  CALL KSAPO(IPAIR)
              END IF
              T0(I1) = TIME
              IMOD = KSLOW(IPAIR)
              CALL KSPOLY(IPAIR,IMOD)
              TDOT2(IPAIR) = -1.0D-20
          END IF
*
*       Produce diagnostic output for large eccentricity.
          IF (ECC.GT.0.9) THEN
              WRITE (75,20)  NAME(I1), TTOT, ECC0, ECC, EMIN, EMAX,
     &                       ECCDOT, EDT, TG, TC1, EDAV, PM
   20         FORMAT (' DECIDE:    NM T E0 E1 EM EX ED EDT TG TC EDAV ',
     &                            'PM ', I5,F9.2,4F7.3,1P,6E9.1)
              CALL FLUSH(75)
          END IF
      ELSE
          DT1 = 10.0/TSTAR
      END IF
*
*       See whether growth time exceeds merger disruption time from TSTAB.
   30 TG1 = TIME + DT1
      IF (IQ.EQ.0) THEN
          TMDIS(NMERGE+1) = MIN(TG1,TMDIS(NMERGE+1))
*     WRITE (6,55) TIME, NAME(I1), NAME(JCOMP), ECC, EMAX,DT1, TG, SEMI
*  55 FORMAT (' WATCH!   T NM NMJ E EX DT1 TG A  ',
*    &                   F10.4,2I7,2F8.4,1P,3E10.2)
      END IF
*     WRITE (6,35)  LIST(1,I1), ECC, EDAV*DT1, EMAX, TG1+TOFF, GAMMA(IPAIR)
*  35 FORMAT (' DECIDE:    NP E DE EMAX TG1 G ',I4,3F7.3,F10.3,1P,E9.1)
*
   40 RETURN
*
      END
