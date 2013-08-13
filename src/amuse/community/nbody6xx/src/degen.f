      SUBROUTINE DEGEN(I1,I2,ICASE)
*
*
*       Binary output for degenerate stars.
*       -----------------------------------
*
      INCLUDE 'common6.h'
      REAL*8 EB(KMAX),SEMI(KMAX),ECC(KMAX),RCM(KMAX),ECM(KMAX)
      CHARACTER*8  WHICH1
      LOGICAL FIRST,FIRST2,FIRST3
      SAVE FIRST,FIRST2,FIRST3
      DATA FIRST,FIRST2,FIRST3 /.TRUE.,.TRUE.,.TRUE./
*
*
*       Skip output for start and end of merger.
      IF (IPHASE.EQ.6.OR.IPHASE.EQ.7) GO TO 50
      IF (ICASE.EQ.7) GO TO 40
*
*       See whether KS binaries contain any degenerate stars.
      NB = 0
      DO 1 IPAIR = I1,I2
          J2 = 2*IPAIR
          J1 = J2 - 1
          IF (KSTAR(J1).GT.9.OR.KSTAR(J2).GT.9) THEN
              NB = NB + 1
          END IF
    1 CONTINUE
*
*       Open unit #4 the first time.
      IF (NB.GT.0.AND.FIRST) THEN
          OPEN (UNIT=4,STATUS='UNKNOWN',FORM='FORMATTED',FILE='DEGEN')
          FIRST = .FALSE.
*
*       Print cluster scaling parameters at start of the run.
          if(rank.eq.0)then
          WRITE (4,2)  RBAR, BODYM*ZMBAR, BODY1*ZMBAR, TSCALE,
     &                 NBIN0, NZERO
    2     FORMAT (/,6X,'MODEL:    RBAR =',F5.1,'  <M> =',F6.2,
     &                 '  M1 =',F6.1,'  TSCALE =',F6.2,
     &                 '  NB =',I4,'  N0 =',I6,//)
*
          WRITE (4,3)
    3     FORMAT (' #   TPHYS    A     E     Rp/R*      P    r',
     &            '     M1   M2    K*         NAME',/)
          end if
      END IF
*
*       Form binding energy and central distance for degenerate stars.
      JPAIR = 0
      IPRINT = 0
      DO 20 IPAIR = I1,I2
          J2 = 2*IPAIR
          J1 = J2 - 1
          IF (KSTAR(J1).LT.10.AND.KSTAR(J2).LT.10) GO TO 20
          JPAIR = JPAIR + 1
          ICM = N + IPAIR
*       Avoid division by zero for merged or synchronous ghost binary.
          IF (BODY(J1).GT.0.0) THEN
              EB(JPAIR) = BODY(J1)*BODY(J2)*H(IPAIR)/
     &                                             (BODY(J1) + BODY(J2))
              SEMI(JPAIR) = -0.5d0*BODY(ICM)/H(IPAIR)
              ECC2 = (1.d0 - R(IPAIR)/SEMI(JPAIR))**2 +
     &                       TDOT2(IPAIR)**2/(BODY(ICM)*SEMI(JPAIR))
              ECC(JPAIR) = SQRT(ECC2)
*       Set zero eccentricity after common envelope stage (still large R).
              IF (R(IPAIR).GT.2.0*SEMI(JPAIR)) THEN
                  ECC(JPAIR) = 0.d0
              END IF
              EB(JPAIR) = MAX(EB(JPAIR),-9.99999d0)
          ELSE
              EB(JPAIR) = 0.d0
              SEMI(JPAIR) = R(IPAIR)
              ECC(JPAIR) = 0.d0
          END IF
          RCM(JPAIR) = SQRT((X(1,ICM) - RDENS(1))**2 +
     &                      (X(2,ICM) - RDENS(2))**2 +
     &                      (X(3,ICM) - RDENS(3))**2)
*       Obtain binding energy (per unit mass) of c.m. motion.
          VJ2 = XDOT(1,ICM)**2 + XDOT(2,ICM)**2 + 
     &                           XDOT(3,ICM)**2
*         POTJ = 0.0
*         DO 5 J = IFIRST,NTOT
*             IF (J.EQ.ICM) GO TO 5
*             RIJ2 = (X(1,ICM) - X(1,J))**2 + (X(2,ICM) - X(2,J))**2 +
*    &                                        (X(3,ICM) - X(3,J))**2 
*             POTJ = POTJ + BODY(J)/SQRT(RIJ2)
*   5     CONTINUE
*         IF (TIME.EQ.TADJ) THEN
*             POTJ = -PHI(ICM)
*         ELSE
              POTJ = 0.d0
*         END IF
          ECM(JPAIR) = 0.5d0*VJ2 - POTJ
*       Check for external tidal field (NB! already included on GRAPE).
*         IF (KZ(14).GT.0) THEN
*             CALL XTRNLV(ICM,ICM)
*             ECM(JPAIR) = ECM(JPAIR) + HT/(BODY(ICM) + 1.0D-20)
*         END IF
          TPHYS = TTOT*TSTAR
          TK = SEMI(JPAIR)*SQRT(ABS(SEMI(JPAIR))/(BODY(ICM) + 1.0d-20))
          TK = DAYS*TK
          TK = MIN(TK,999999.9d0)
          RBIG = MAX(RADIUS(J1),RADIUS(J2))
          RATIO = SEMI(JPAIR)*(1.d0 - ECC(JPAIR))/RBIG
          RATIO = MIN(RATIO,99.9d0)
          IF (SEMI(JPAIR).LT.0.0.AND.RATIO.GT.5.0) GO TO 20
          SEMI(JPAIR) = SEMI(JPAIR)*RBAR*AU
          IF (IPRINT.EQ.0.AND.ICASE.EQ.0) THEN
              if(rank.eq.0) WRITE (4,*)
          END IF
          if(rank.eq.0)then
          WRITE (4,10)  ICASE, TPHYS, SEMI(JPAIR), ECC(JPAIR), RATIO,
     &                  TK, RCM(JPAIR), BODY(J1)*ZMBAR, BODY(J2)*ZMBAR,
     &                  KSTAR(J1), KSTAR(J2), KSTAR(ICM),
     &                  NAME(J1), NAME(J2)
   10     FORMAT (I2,F8.1,F8.2,F7.3,F6.1,F9.1,F6.2,2F5.1,3I4,2I6)
          end if
          IPRINT = IPRINT + 1
   20 CONTINUE
*
*       Close file at end of main output.
      IF (IPRINT.GT.0.AND.ICASE.EQ.0) THEN
          CALL FLUSH(4)
      END IF
*
*       Search for neutron stars at main output.
      IF (ICASE.EQ.0) THEN
          DO 30 J = 1,N
              IF (KSTAR(J).EQ.13) THEN
                  IF (FIRST2) THEN
                      OPEN (UNIT=33,STATUS='UNKNOWN',FORM='FORMATTED',
     &                                                 FILE='NS')
                      FIRST2 = .FALSE.
                  END IF
                  IF (J.LT.IFIRST) THEN
                      WHICH1 = ' BINARY '
                      I = KVEC(J) + N
                      VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + 
     &                                     XDOT(3,I)**2
                  ELSE
                      WHICH1 = ' SINGLE '
                      VI2 = XDOT(1,J)**2 + XDOT(2,J)**2 + 
     &                                     XDOT(3,J)**2
                  END IF
                  if(rank.eq.0)then
                  WRITE (33,25)  WHICH1, J, NAME(J), IFIRST, KSTAR(J),
     &                           TPHYS, SQRT(VI2)*VSTAR
   25             FORMAT (1X,A8,'NS','  J NAM I* K* TPH V ',
     &                                  2I6,I5,I4,2F7.1)
                  CALL FLUSH(33)
                  end if
              ELSE IF (KSTAR(J).GT.13) THEN
                  IF (J.LT.IFIRST) THEN
                      JCM = N + KVEC(J)
                      IF (NAME(JCM).LT.0) GO TO 30
                  END IF
                  IF (FIRST3) THEN
                      OPEN (UNIT=34,STATUS='UNKNOWN',FORM='FORMATTED',
     &                                                 FILE='BH')
                      FIRST3 = .FALSE.
                  END IF
                  VI2 = XDOT(1,J)**2 + XDOT(2,J)**2 + 
     &                                 XDOT(3,J)**2
                  VI = SQRT(VI2)
                  WHICH1 = ' AIC    '
                  IF (KSTAR(J).EQ.14) WHICH1 = ' BH     '
                  IF (KSTAR(J).EQ.15) THEN
*       Ensure that errant massless remnant will escape.
                      RI = SQRT(X(1,J)**2 + X(2,J)**2 + X(3,J)**2)
                      DO 26 L = 1,3
                          X0(L,J) = 1000.0*RSCALE*X(L,J)/RI
                          X(L,J) = X0(L,J)
                          X0DOT(L,J) = SQRT(0.004*ZMASS/RSCALE)*
     &                                 XDOT(L,J)/VI
                          XDOT(L,J) = X0DOT(L,J)
   26                 CONTINUE
                  ELSE
                      VI = VI*VSTAR
                      if(rank.eq.0)then
                      WRITE (34,28) WHICH1,J,NAME(J),KSTAR(J),TPHYS,VI
   28                 FORMAT (1X,A8,'J NAM K* TPH V ',2I6,I4,2F7.1)
                      CALL FLUSH(34)
                      end if
                  END IF
              END IF
   30     CONTINUE
      END IF
*
*       Print single neutron stars at escape time.
   40 IF (ICASE.EQ.7) THEN
          J = I1
          VI2 = XDOT(1,J)**2 + XDOT(2,J)**2 + XDOT(3,J)**2
          VI = SQRT(VI2)*VSTAR
          WHICH1 = ' ESCAPE '
          if(rank.eq.0)then
          WRITE (33,25)  WHICH1, J, NAME(J), IFIRST, KSTAR(J), TPHYS, VI
          end if
      END IF
*
*       Include counter for doubly degenerate binary.
      IF (ICASE.EQ.3.OR.ICASE.EQ.4) THEN
          IPAIR = I1
          J1 = 2*IPAIR - 1
          J2 = J1 + 1
          IF (KSTAR(J1).GE.10.AND.KSTAR(J2).GE.10) THEN
              NDD = NDD + 1
              A = -0.5d0*SU*BODY(N+IPAIR)/H(IPAIR)
              if(rank.eq.0)then
              WRITE (6,48)  IPAIR, NAME(J1), NAME(J2), KSTAR(J1),
     &                      KSTAR(J2), KSTAR(N+IPAIR), R(IPAIR), A, TK
   48         FORMAT (' NEW DD    KS NM K* R A P ',I4,2I6,3I4,1P,3E10.2)
              end if
          END IF
      END IF
*
   50 RETURN
*
      END
