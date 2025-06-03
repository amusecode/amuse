      SUBROUTINE HIARCH(IPAIR)
*
*
*       Hierarchical system diagnostics.
*       --------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  A1(3),A2(3),XREL(3),VREL(3),EI(3),HI(3),HO(3),
     &        TK0(100),TLAST(100),TTERM(100)
      INTEGER  LAST(100),IPRINT(100)
      LOGICAL  FIRST
      SAVE  FIRST,TLAST,TTERM,RAP,TK0,LAST,IPRINT,NL
      DATA  FIRST /.TRUE./
      DATA  LAST,IPRINT  /200*0/
*
*
*       Open unit #10 the first time.
      IF (FIRST) THEN
          OPEN (UNIT=10,STATUS='UNKNOWN',FORM='FORMATTED',FILE='HIARCH')
          FIRST = .FALSE.
*
*       Print cluster scaling parameters at start of the run.
          IF (rank.eq.0.and.NMERG.EQ.0) THEN
              WRITE (10,1)  RBAR, BODYM*ZMBAR, BODY1*ZMBAR, TSCALE,
     &                      NBIN0, NZERO
    1         FORMAT (/,6X,'MODEL:    RBAR =',F5.1,'  <M> =',F6.2,
     &                     '  M1 =',F6.1,'  TSCALE =',F6.2,
     &                     '  NB =',I4,'  N0 =',I6,//)
          END IF
*
*       Define the saved variables in case of COMMON dump during merger.
          TK0(1) = 1.0
          RAP = 1.0
          NL = 0
          TLAST(1) = 0.0
          TTERM(1) = 0.0
          if(rank.eq.0)then
          WRITE (10,2)
    2     FORMAT ('      TIME   A         A1        E1    PMIN',
     &            '    P1/P0  Q     PCR   MB    Q1     I      NAME',
     &            '        KC',/)
          WRITE (10,3)
    3     FORMAT ('             r/Rc                           ',
     &           '           Rc   GAM   IT        N0   N1   KS  NAME',/)
          end if
      END IF
*
*       Set c.m. and component indices.
      I = N + IPAIR
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      TTOT = TIME + TOFF
*
*       Skip printing termination at line #2 until after first new merger.
      IF (IPHASE.EQ.7.AND.NL.EQ.0) THEN
          GO TO 30
      END IF
*
*       Determine index of current merger (may have been set previously).
      IL = 0
      DO 4 L = 1,NL
          IF (NAME(I1).EQ.LAST(L)) THEN
              IL = L
          END IF
    4 CONTINUE
*
*       Increase membership and set identifier = 0 temporarily for new case.
      IF (IL.EQ.0) THEN
          NL = NL + 1
          IL = NL
          LAST(IL) = 0
      END IF
*
*       Decide whether new merger or termination should be printed.
      IF (IPHASE.EQ.6) THEN
          IF (NAME(I1).EQ.LAST(IL)) THEN
              TLAST(IL) = TTOT
*       Avoid repeated output of same system if terminated recently (< TCR).
              IF (TTOT - TTERM(IL).LT.TCR) THEN
                  IPRINT(IL) = 0
                  GO TO 30
              ELSE
                  IPRINT(IL) = 1
              END IF
          ELSE
*       Save identity and merge time and set print indicator for line #2.
              LAST(IL) = NAME(I1)
              TLAST(IL) = TTOT
              TTERM(IL) = 0.0D0
              IPRINT(IL) = 1
          END IF
      ELSE
          TTERM(IL) = TTOT
*       Skip output of line #2 if line #1 was not printed.
          IF (IPRINT(IL).EQ.0) THEN
              GO TO 30
          END IF
          IPRINT(IL) = 0
      END IF
*
*       Form semi-major axis and eccentricity (skip hyperbolic orbits).
      SEMI = -0.5*BODY(I)/H(IPAIR)
      IF (SEMI.LT.0.0) GO TO 30
      ECC2 = (1.0D0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      ECC = SQRT(ECC2)
*       Skip highly eccentric fly-by's which are not formally stable.
      IF (ECC.GT.0.99) GO TO 30
      TK = TWOPI*SEMI*SQRT(SEMI/BODY(I))
*
*       Decide between new and terminated merger.
      IF (IPHASE.EQ.6) THEN
*       Ensure that unperturbed binary is resolved.
          IF (LIST(1,I1).EQ.0.OR.X(1,I1).EQ.X(1,I2)) THEN
              CALL RESOLV(IPAIR,1)
          END IF
          Q = BODY(I)/BODY(JCOMP)
          Q1 = MAX(BODY(I2)/BODY(I1),BODY(I1)/BODY(I2))
          RAP = SEMI*(1.0D0 + ECC)
          RIJ2 = 0.0
          VIJ2 = 0.0
          RDOT = 0.0
          A12 = 0.0
          A22 = 0.0
          A1A2 = 0.0
          RI2 = 0.0
          VI2 = 0.0
          RVI = 0.0
          DO 5 K = 1,3
              RIJ2 = RIJ2 + (X(K,I) - X(K,JCOMP))**2
              VIJ2 = VIJ2 + (XDOT(K,I) - XDOT(K,JCOMP))**2
              RDOT = RDOT + (X(K,I) - X(K,JCOMP))*
     &                                       (XDOT(K,I) - XDOT(K,JCOMP))
              K1 = K + 1
              IF (K1.GT.3) K1 = 1
              K2 = K1 + 1
              IF (K2.GT.3) K2 = 1
              A1(K) = (X(K1,I1) - X(K1,I2))*(XDOT(K2,I1) - XDOT(K2,I2))
     &              - (X(K2,I1) - X(K2,I2))*(XDOT(K1,I1) - XDOT(K1,I2))
              A2(K) = (X(K1,JCOMP) - X(K1,I))*
     &                                     (XDOT(K2,JCOMP) - XDOT(K2,I))
     &              - (X(K2,JCOMP) - X(K2,I))*
     &                                     (XDOT(K1,JCOMP) - XDOT(K1,I))
              A12 = A12 + A1(K)**2
              A22 = A22 + A2(K)**2
              A1A2 = A1A2 + A1(K)*A2(K)
*       Form relative vectors and scalars for inner binary.
              XREL(K) = X(K,I1) - X(K,I2)
              VREL(K) = XDOT(K,I1) - XDOT(K,I2)
              RI2 = RI2 + XREL(K)**2
              VI2 = VI2 + VREL(K)**2
              RVI = RVI + XREL(K)*VREL(K)
    5     CONTINUE
*
*      Evaluate orbital parameters for outer orbit (skip hyperbolic case).
          RIJ = SQRT(RIJ2)
          ZMB = BODY(I) + BODY(JCOMP)
          SEMI1 = 2.0/RIJ - VIJ2/ZMB
          SEMI1 = 1.0/SEMI1
          IF (SEMI1.LT.0.0) GO TO 30
          ECC1 = SQRT((1.0 - RIJ/SEMI1)**2 + RDOT**2/(SEMI1*ZMB))
          PMIN = SEMI1*(1.0D0 - ECC1)
          TK1 = TWOPI*SEMI1*SQRT(SEMI1/ZMB)
          PMIN = MIN(PMIN/RAP,999.0D0)
          TK1 = MIN(TK1/TK,99999.0D0)
*       Determine inclination (8 bins of 22.5 degrees).
          FAC = A1A2/SQRT(A12*A22)
          FAC = ACOS(FAC)
*         II = 1 + FAC*360.0/(TWOPI*22.5)
          IDEG = 360.0*FAC/TWOPI
*
*       Construct Runge-Lenz vector (Heggie & Rasio, 1995, IAU174, Eq.(5)).
          EI2 = 0.0
          DO 6 K = 1,3
              EI(K) = (VI2*XREL(K) - RVI*VREL(K))/BODY(I) -
     &                                                 XREL(K)/SQRT(RI2)
              EI2 = EI2 + EI(K)**2
    6     CONTINUE
*
*       Define unit vectors for inner eccentricity and angular momenta.
          COSJ = 0.0
          SJSG = 0.0
          DO 8 K = 1,3
              EI(K) = EI(K)/SQRT(EI2)
              HI(K) = A1(K)/SQRT(A12)
              HO(K) = A2(K)/SQRT(A22)
              COSJ = COSJ + HI(K)*HO(K)
              SJSG = SJSG + EI(K)*HO(K)
    8     CONTINUE
*
*       Form the expressions A & Z.
          A = COSJ*SQRT(1.0 - EI2)
          Z = (1.0 - EI2)*(2.0 - COSJ**2) + 5.0*EI2*SJSG**2
*
*       Obtain maximum inner eccentricity (Douglas Heggie, Sept. 1995).
          Z2 = Z**2 + 25.0 + 16.0*A**4 - 10.0*Z - 20.0*A**2 - 8.0*A**2*Z
          EMAX = ONE6*(Z + 1.0 - 4.0*A**2 + SQRT(Z2))
          EMAX = SQRT(EMAX)
          KCM = KSTAR(I)
          IF (NAME(I).LT.0) KCM = -10
          NEWHI = NEWHI + 1
*
          if(rank.eq.0)then
          WRITE (10,10)  TTOT, SEMI, SEMI1, ECC1, PMIN, TK1, Q,
     &                   PCRIT/SEMI, BODY(I)/BODYM, Q1, IDEG,
     &                   NAME(I1), NAME(I2), NAME(JCOMP), KCM
   10     FORMAT (/,' #1',F8.1,1P,2E10.2,0P,F6.2,F7.2,F8.1,3F6.2,F5.1,
     &                    I5,2I5,I6,I4)
*
          WRITE (30,10)  TTOT, SEMI, SEMI1, ECC1, PMIN, TK1, Q,
     &                   PCRIT/SEMI, BODY(I)/BODYM, Q1, IDEG,
     &                   NAME(I1), NAME(I2), NAME(JCOMP), KCM
          end if
*
          RBIG = MAX(RADIUS(I1),RADIUS(I2))
          PMIN = SEMI*(1.0 - ECC)
          PMIN2 = SEMI*(1.0 - EMAX)
          IF (KZ(19).EQ.3) THEN
              R1 = PMIN/RBIG
              R2 = PMIN2/RBIG
          ELSE
              R1 = 0.0
              R2 = 0.0
          END IF
          ZI = FAC
          EM = SQRT(SIN(ZI)**2 + EI2*COS(ZI)**2)
          if(rank.eq.0)
     &    WRITE (30,11)  SQRT(EI2), EMAX, KSTAR(I1), KSTAR(I2),KSTAR(I),
     &                   SEMI, A, Z, PMIN2, R1, R2, EM
   11     FORMAT (' INNER:    E EMAX K* SEMI A Z PM2 R1 R2 EM ',
     &                        2F8.4,3I3,1P,4E10.2,0P,2F7.1,F7.3)
          CALL FLUSH(30)
*       Save parameters for termination diagnostics.
          TK0(IL) = TK
*
      ELSE IF (IPHASE.EQ.7) THEN
          PMIN = SEMI*(1.0D0 - ECC)
          PMIN = MIN(PMIN/RAP,999.0D0)
          IF (TK0(IL).LE.0.0D0) TK0(IL) = TK
          TK1 = MIN(TK/TK0(IL),99999.0D0)
          NK = (TTOT - TLAST(IL))/TK0(IL)
          NK1 = (TTOT - TLAST(IL))/TK
          NK = MIN(NK,99999999)
          NK1 = MIN(NK1,9999)
          NK = MAX(NK,0)
          RI = 0.0
          DO 15 K = 1,3
              RI = RI + (X(K,I) - RDENS(K))**2
   15     CONTINUE
          RR = RC/RSCALE
*
*       Define type index (0, 1 or KSTAR; note: #I1 is inner c.m.).
          IT = 0
          IF (TIME.GT.MIN(TEV(I1),TEV(I2),TEV(I))) IT = 1
          IF (KSTAR(I1).NE.0.AND.IT.EQ.0) IT = KSTAR(I1)
*
          if(rank.eq.0)then
          WRITE (10,20)  TTOT, SQRT(RI)/RC, SEMI, ECC, PMIN, TK1,
     &                   RR, GAMMA(IPAIR), IT, NK, NK1, NPAIRS, NAME(I2)
   20     FORMAT (' #2',F8.1,1P,2E10.2,0P,F6.2,F7.2,F8.1,2F6.2,I4,I10,
     &                  3I6)
          CALL FLUSH(10)
          WRITE (30,20)  TTOT, SQRT(RI)/RC, SEMI, ECC, PMIN, TK1,
     &                   RR, GAMMA(IPAIR), IT, NK, NK1, NPAIRS, NAME(I2)
          end if
      END IF
*
*       Restrict memory of previous cases by removing the oldest one.
   30 IF (NL.GE.100) THEN
          LT = 1
*       Skip any current mergers (defined by TTERM = 0.0).
   32     IF (TTERM(LT).LE.0.0D0) THEN
              LT = LT + 1
              GO TO 32
          END IF
          DO 35 L = LT,NL-1
              TLAST(L) = TLAST(L+1)
              TTERM(L) = TTERM(L+1)
              TK0(L) = TK0(L+1)
              LAST(L) = LAST(L+1)
              IPRINT(L) = IPRINT(L+1)
   35     CONTINUE
          NL = NL - 1
      END IF
*
      RETURN
*
      END
