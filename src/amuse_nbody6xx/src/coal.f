      SUBROUTINE COAL(IPAIR,KW1,KW2,MASS)
*
*
*       Coalescence of Roche/CE binary.
*       -------------------------------
*
      INCLUDE 'common6.h'
      CHARACTER*8  WHICH1
      REAL*8  CM(6)
      REAL*8  MASS(2)
      LOGICAL  FIRST
      SAVE  FIRST
      DATA  FIRST /.TRUE./
*
*
*       Distinguish between KS and chain regularization.
      IF (IPAIR.GT.0) THEN
*       Define discrete time for prediction & new polynomials (T <= TBLOCK).
          I = N + IPAIR
          DT = 0.1d0*STEP(I)
          IF (DT.GT.2.4E-11) THEN
              TIME2 = TIME - TPREV
              CALL STEPK(DT,DTN)
              TIME = TPREV + INT((TIME2 + DT)/DTN)*DTN
              TIME = MIN(TBLOCK,TIME)
          ELSE
              TIME = MIN(T0(I) + STEP(I),TBLOCK)
          END IF
*
*       Set zero energy (EB correction done in routines EXPEL & EXPEL2).
          EB = 0.d0
          RCOLL = R(IPAIR)
          DMIN2 = MIN(DMIN2,RCOLL)
          VINF = 0.0
*
*       Define indicator for different cases, including hyperbolic KS.
          IF ((KSTAR(I).LE.10.AND.IQCOLL.NE.0).OR.IQCOLL.EQ.-2) THEN
              WHICH1 = ' BINARY '
              IQCOLL = 0
              EB1 = BODY(2*IPAIR-1)*BODY(2*IPAIR)*H(IPAIR)/BODY(I)
          ELSE
              WHICH1 = '  ROCHE '
*       Save energy for correction (Roche COAL but not via CMBODY & EXPEL).
              IF (IQCOLL.EQ.0) THEN
                  EB = BODY(2*IPAIR-1)*BODY(2*IPAIR)*H(IPAIR)/BODY(I)
                  EB1 = EB
                  EGRAV = EGRAV + EB
              END IF
              IQCOLL = 3
          END IF
          IF (H(IPAIR).GT.0.0) THEN
              WHICH1 = ' HYPERB '
              NHYP = NHYP + 1
              IQCOLL = -1
              VINF = SQRT(2.0*H(IPAIR))*VSTAR
              EB1 = BODY(2*IPAIR-1)*BODY(2*IPAIR)*H(IPAIR)/BODY(I)
          END IF
*
*       Check optional diagnostics for binary evolution.
          IF (KZ(8).GT.3) THEN
              CALL BINEV(IPAIR)
          END IF
*
*       Remove any circularized KS binary from the CHAOS/SYNCH table.
          IF (KSTAR(I).GE.10.AND.NCHAOS.GT.0) THEN
              II = -I
              CALL SPIRAL(II)
          END IF
*
*       Terminate KS pair and set relevant indices for collision treatment.
          IPHASE = 9
          KSPAIR = IPAIR
          T0(2*IPAIR-1) = TIME
          CALL KSTERM
          I1 = 2*NPAIRS + 1
          I2 = I1 + 1
      ELSE
*       Copy dominant indices, two-body separation and binding energy.
          I1 = JLIST(1)
          I2 = JLIST(2)
          RCOLL = DMINC
          EB = 0.0
          VINF = 0.d0
          IQCOLL = 5
          WHICH1 = '  CHAIN '
*       Note that new chain TIME already quantized in routine CHTERM.
      END IF
*
*       Define global c.m. coordinates & velocities from body #I1 & I2.
      ICOMP = I1
      ZM = BODY(I1) + BODY(I2)
      DO 5 K = 1,3
          CM(K) = (BODY(I1)*X(K,I1) + BODY(I2)*X(K,I2))/ZM
          CM(K+3) = (BODY(I1)*XDOT(K,I1) + BODY(I2)*XDOT(K,I2))/ZM
    5 CONTINUE
*
*       Form central distance (scaled by RC), central distance and period.
      RI2 = 0.d0
      RIJ2 = 0.d0
      VIJ2 = 0.d0
      DO 10 K = 1,3
          RI2 = RI2 + (X(K,I1) - RDENS(K))**2
          RIJ2 = RIJ2 + (X(K,I1) - X(K,I2))**2
          VIJ2 = VIJ2 + (XDOT(K,I1) - XDOT(K,I2))**2
   10 CONTINUE
      RI = SQRT(RI2)
      RIJ = SQRT(RIJ2)
      SEMI = 2.d0/RIJ - VIJ2/ZM
      SEMI = 1.d0/SEMI
      TK = DAYS*SEMI*SQRT(ABS(SEMI)/ZM)
      ECC = 1.0 - RCOLL/SEMI
      IF (IQCOLL.EQ.5) THEN
          EB1 = -0.5*BODY(I1)*BODY(I2)/SEMI
      END IF
*
*        Choose appropriate reference body containing #ICH neighbour list.
      IF (IPAIR.GT.0) THEN
          ICM = I1
      ELSE
          ICM = JLIST(3)
          IF (LIST(1,I1).GT.LIST(1,ICM)) ICM = I1
          IF (LIST(1,I2).GT.LIST(1,ICM)) ICM = I2
*       Include possible 4th chain member just in case.
          IF (JLIST(4).GT.0) THEN
              I4 = JLIST(4)
              IF (LIST(1,I4).GT.LIST(1,ICM)) ICM = 4
          END IF
      END IF
*
*       Form perturber list for corrections and add chain #I1 to NB lists.
      NNB = LIST(1,ICM)
      L2 = 1
      JMIN = N + 1
      RIJ2 = 1.0d+10
*       Copy neighbour list of former c.m. (less any #I2) to JPERT.
      DO 15 L = 1,NNB
          JPERT(L) = LIST(L+1,ICM)
          IF (JPERT(L).EQ.I2) THEN      !unlikely condition but no harm.
              L2 = L
              JPERT(L) = N
              IF (I2.EQ.N) JPERT(L) = N - 1
          ELSE
*       Determine closest external member for possible KS.
              JJ = JPERT(L)
              RI2 = 0.d0
              DO 12 K = 1,3
                  RI2 = RI2 + (CM(K) - X(K,JJ))**2
   12         CONTINUE
              IF (RI2.LT.RIJ2) THEN
                  JMIN = JJ
                  RIJ2 = RI2
              ENDIF
          END IF
   15 CONTINUE
*
*       Restore at least chain collider mass in local neighbour lists.
      IF (IPAIR.LE.0) THEN
          NNB = NNB + 1
          JPERT(NNB) = ICM
*       Include the case of old c.m. coming from #I1 or #I2.
          IF (ICM.EQ.I1.OR.ICM.EQ.I2) JPERT(NNB) = JLIST(3)
          JLIST(1) = I1
          IF (ICM.EQ.I1) JLIST(1) = JLIST(3)
          CALL NBREST(ICM,1,NNB)
      END IF
*
*       Evaluate potential energy with respect to colliding bodies.
      JLIST(1) = I1
      JLIST(2) = I2
      CALL NBPOT(2,NNB,POT1)
*
*       Specify new mass from sum and initialize zero mass ghost in #I2.
      ZMNEW = (MASS(1) + MASS(2))/ZMBAR
      DM = ZM - ZMNEW
      IF (DM.LT.1.0D-10) DM = 0.d0
*       Delay inclusion of any mass loss until after energy correction.
      ZM1 = BODY(I1)*ZMBAR
      ZM2 = BODY(I2)*ZMBAR
      BODY(I1) = ZM
      BODY(I2) = 0.d0
      NAME1 = NAME(I1)
      NAME2 = NAME(I2)
      IF(BODY0(I1).LT.BODY0(I2))THEN
         BODY0(I1) = BODY0(I2)
         EPOCH(I1) = EPOCH(I2)
         TEV(I1) = TEV(I2)
         SPIN(I1) = SPIN(I2)
         RADIUS(I1) = RADIUS(I2)
         NAME2 = NAME(I2)
         NAME(I2) = NAME(I1)
         NAME(I1) = NAME2
      ENDIF
      T0(I1) = TIME
      T0(I2) = TADJ + DTADJ
      CALL DTCHCK(TIME,STEP(I2),DTK(40))
*
*       Start new star from current time unless ROCHE case with TEV0 > TIME.
      TEV0(I1) = MAX(TEV0(I1),TEV0(I2))
      IF(IQCOLL.NE.3) TEV(I1) = MAX(TIME,TEV0(I1))
      TEV(I2) = 1.0d+10
      VI = SQRT(XDOT(1,I2)**2 + XDOT(2,I2)**2 + XDOT(3,I2)**2)
*
*	Set T0 = TIME for any other chain members.
      IF (IPAIR.LT.0) THEN
         DO 18 L = 1,NCH
            J = JLIST(L)
            IF (J.NE.I1.AND.J.NE.I2) THEN
               T0(J) = TIME
            END IF
   18    CONTINUE
      END IF
*
*       Check that a mass-less primary has type 15 for kick velocity.
      IF(ZMNEW*ZMBAR.LT.0.001.AND.KW1.NE.15)THEN
         if(rank.eq.0)then
         WRITE(6,*)' ERROR COAL: mass1 = 0.0 and kw1 is not equal 15'
         WRITE(6,*)' I KW ',I1,KW1
         end if
         STOP
      END IF
*
      DO 20 K = 1,3
          X(K,I1) = CM(K)
          X0(K,I1) = CM(K)
          XDOT(K,I1) = CM(K+3)
          X0DOT(K,I1) = CM(K+3)
*       Ensure that ghost will escape next output (far from fast escapers).
          X0(K,I2) = MIN(1.0d+04 + (X(K,I2)-RDENS(K)),
     &                   1000.d0*RSCALE*(X(K,I2)-RDENS(K))/RI)
          X(K,I2) = X0(K,I2)
          X0DOT(K,I2) = SQRT(0.004d0*ZMASS/RSCALE)*XDOT(K,I2)/VI
          XDOT(K,I2) = X0DOT(K,I2)
          F(K,I2) = 0.d0
          FDOT(K,I2) = 0.d0
          D0(K,I2) = 0.0
          D1(K,I2) = 0.0
          D2(K,I2) = 0.d0
          D3(K,I2) = 0.d0
          D0R(K,I2) = 0.0
          D1R(K,I2) = 0.0
          D2R(K,I2) = 0.d0
          D3R(K,I2) = 0.d0
   20 CONTINUE
*
*       Obtain potential energy w.r.t. new c.m. and apply tidal correction.
      CALL NBPOT(1,NNB,POT2)
      DP = POT2 - POT1
      ECOLL = ECOLL + DP
*
      J2 = JPERT(L2)
      JPERT(L2) = I2
*
*       Remove the ghost particle #I2 from perturber lists containing #I1.
      JLIST(1) = I2
      CALL NBREM(I1,1,NNB)
      JLIST(1) = I1
      JPERT(L2) = J2
*
*       Include correction procedure in case of mass loss (cf routine MIX).
      IF (KZ(19).GE.3.AND.DM.GT.0.0) THEN
*
*       Reduce mass of composite body and update total mass (check SN mass).
          BODY(I1) = ZMNEW
          BODY(I1) = MAX(BODY(I1),0.d0)
          IF (ABS(BODY(I1)).LT.1.0d-10) TEV(I1) = 1.0d+10
          ZMASS = ZMASS - DM
*
*       Adopt ILIST procedure from NBODY4 with NNB available in FCORR.
          ILIST(1) = NNB
          DO 22 L = 1,NNB
              ILIST(L+1) = JPERT(L)
   22     CONTINUE
*
*       Delay velocity kick until routine MDOT on type 13/14/15 in ROCHE.
*       [Note: NS,BH should not be created here unless possibly for
*        AIC of WD or TZ object from COMENV.]
          KW = KW1
          IF (KW1.EQ.13.OR.KW1.EQ.14) THEN
             IF(KSTAR(I1).GE.13.OR.KSTAR(I2).GE.13) KW = 0
*            IF(KSTAR(I1).GE.10.OR.KSTAR(I2).GE.10) KW = 0
          END IF
          IF (KW1.GE.10.AND.KW1.LE.12) THEN
             IF(KSTAR(I1).GE.10.OR.KSTAR(I2).GE.10) KW = 0
          END IF
*
*       Perform total force & energy corrections (new polynomial set later).
          CALL FCORR(I1,DM,KW)
*
*       Specify commensurate time-step (not needed for BODY(I1) = 0).
          CALL DTCHCK(TIME,STEP(I1),DTK(40))
*
*       Set IPHASE = -3 to preserve ILIST.
          IPHASE = -3
*
*       Initialize new polynomials of neighbours & #I1 for DM > 0.1 DMSUN.
          IF (DM*ZMBAR.GT.0.1) THEN
*
*       Include body #I1 at the end (counting from location #2).
              NNB2 = NNB + 2
              ILIST(NNB2) = I1
*
*       Obtain new F & FDOT and time-steps.
              DO 30 L = 2,NNB2
                  J = ILIST(L)
                  IF (L.EQ.NNB2) THEN
                      J = I1
                  ELSE IF (T0(J).LT.TIME) THEN
                      CALL XVPRED(J,-2)
                      CALL DTCHCK(TIME,STEP(J),DTK(40))
                  END IF
                  DO 25 K = 1,3
                      X0DOT(K,J) = XDOT(K,J)
   25             CONTINUE
*       Create ghost for rare case of massless first component.
                  IF (BODY(J).EQ.0.0D0) THEN
                      DO 26 K = 1,3
                         X0(K,I1) = MIN(1.0d+04 + (X(K,I1)-RDENS(K)),
     &                             1000.d0*RSCALE*(X(K,I1)-RDENS(K))/RI)
                          X(K,I1) = X0(K,I1)
                          X0DOT(K,I1) = SQRT(0.004d0*ZMASS/RSCALE)*
     &                                               XDOT(K,I1)/VI
                          XDOT(K,I1) = X0DOT(K,I1)
                          F(K,I1) = 0.d0
                          FDOT(K,I1) = 0.d0
                          D2(K,I1) = 0.d0
                          D3(K,I1) = 0.d0
   26                 CONTINUE
                      T0(I1) = TADJ + DTADJ
                      if(rank.eq.0) WRITE (6,28)  NAME(I1), KW1
   28                 FORMAT (' MASSLESS PRIMARY!    NAM KW ',I8,I4)
                  ELSE
                      CALL FPOLY1(J,J,0)
                      CALL FPOLY2(J,J,0)
                  END IF
   30         CONTINUE
          END IF
*         TPREV = TIME - STEPX
      END IF
*
*       See whether closest neighbour forms a KS pair (skip chain).
      IF (IPAIR.GT.0.AND.BODY(I1).GT.0.0D0) THEN
          IF (JMIN.LE.N.AND.RIJ2.LT.RMIN2) THEN
              DO 35 K = 1,3
                  X0DOT(K,JMIN) = XDOT(K,JMIN)
   35         CONTINUE
              ICOMP = MIN(I1,JMIN)
              JCOMP = MAX(I1,JMIN)
              CALL KSREG
              if(rank.eq.0)then
              WRITE (6,36) NAME(ICOMP), NAME(JCOMP), LIST(1,2*NPAIRS-1),
     &                     R(NPAIRS), H(NPAIRS), STEP(NTOT)
   36         FORMAT (' COAL KS    NM NP R H DTCM  ',2I6,I4,1P,3E10.2)
              end if
              I2 = JMIN
*       Note that T0(I2) may not have a large value after #I2 is exchanged.
              T0(I2) = TADJ + DTADJ
          ELSE
*       Initialize force polynomial for new single body.
              CALL FPOLY1(ICOMP,ICOMP,0)
              CALL FPOLY2(ICOMP,ICOMP,0)
          END IF
      END IF
*
*       Update energy loss & collision counters (EB = 0 for CHAIN COAL).
      ECOLL = ECOLL + EB
      E(10) = E(10) + EB
      NPOP(9) = NPOP(9) + 1
      NCOAL = NCOAL + 1
      EB = EB1
*
*       Open unit #12 the first time.
      IF (rank.eq.0.and.FIRST) THEN
          OPEN (UNIT=12,STATUS='UNKNOWN',FORM='FORMATTED',FILE='COAL')
          FIRST = .FALSE.
*
*       Print cluster scaling parameters at start of the run.
          IF (NCOAL.EQ.1) THEN
              WRITE (12,40)  RBAR, BODYM*ZMBAR, BODY1*ZMBAR, TSCALE,
     &                       NBIN0, NZERO
   40         FORMAT (/,6X,'MODEL:    RBAR =',F5.1,'  <M> =',F6.2,
     &                     '  M1 =',F6.1,'  TSCALE =',F6.2,
     &                     '  NB =',I5,'  N0 =',I6,//)
              WRITE (12,45)
   45         FORMAT ('    TIME  NAME  NAME  K1  K2  IQ  M1   M2',
     &                '   DM    R1     R2    r/Rc   R     ECC      P',/)
          END IF
      END IF
*
      if(rank.eq.0)then
      WRITE (12,50)  TTOT, NAME1, NAME2, KSTAR(I1), KSTAR(I2), IQCOLL,
     &               ZM1, ZM2, DM*ZMBAR, RADIUS(I1)*SU, RADIUS(I2)*SU,
     &               RI/RC, RIJ*SU, ECC, TK
   50 FORMAT (1X,F7.1,2I6,3I4,3F5.1,2F7.2,F6.1,F7.2,F9.5,1P,E9.1)
      CALL FLUSH(12)
      end if
*
      if(rank.eq.0)then
      WRITE (6,55)  WHICH1, IQCOLL, NAME1, NAME2, KSTAR(I1), KSTAR(I2),
     &              KW1, ZMNEW*ZMBAR, RCOLL, EB, DP, DM*ZMBAR, VINF
   55 FORMAT (/,A8,'COAL    IQ =',I3,'  NAME =',2I6,'  K* =',3I3,
     &             '  M =',F6.2,'  RCOLL =',1P,E8.1,'  EB =',E9.1,
     &             '  DP =',E9.1,'  DM =',0P,F6.2,'  VINF =',F4.1)
      end if
*
      KSTAR(I1) = KW1
      KSTAR(I2) = 15
*       Specify IPHASE < 0 for new sorting.
      IPHASE = -1
      IQCOLL = 0
*
      RETURN
*
      END
