      SUBROUTINE MLOSS
*
*
*       Mass loss from evolving stars.
*       ------------------------------
*
*       Original scheme of Elena Terlevich 1983.
*       ----------------------------------------
*
      INCLUDE 'common6.h'
      REAL*8  A(9)
*
*
*       Find the heaviest body (exclude any merged binary components).
   99 IMAX = 0
    1 BMAX = 0.0D0
      DO 2 I = 1,N
          IF (BODY(I).LT.BMAX.OR.I.LE.IMAX) GO TO 2
          BMAX = BODY(I)
          IMAX = I
    2 CONTINUE
*
      ZMSTAR = ZMBAR*BMAX
*       Obtain total evolution time for ZMSTAR in solar masses.
      TMS = (2.55D+03 + 6.69D2*ZMSTAR**2.5D0 + ZMSTAR**4.5D0)/
     &      (3.27D-02*ZMSTAR**1.5D0 + 3.46D-01*ZMSTAR**4.5D0)
      TG = 0.15*TMS
      THE = TMS*1.37*ZMSTAR**(-0.881D0)
      TEV1 = TMS + TG + THE
*       Scale the evolution time to model units.
      TMDOT = TEV1/TSCALE
*
      IF (TMDOT.GT.TIME) THEN
*       Update the maximum mass.
          BODY1 = BMAX
*       Set phase indicator = -1 for new time-step list in INTGRT.
          IPHASE = -1
          RETURN
      END IF
*
      IBODY = IMAX
*       Include special treatment for regularized components.
      IF (IMAX.GT.2*NPAIRS) GO TO 4
      IPAIR = KVEC(IMAX)
      IBODY = N + IPAIR
*       Continue search if merged binary component is identified.
      IF (NAME(IBODY).LT.0) GO TO 1
*
*       Obtain current coordinates & velocities (unperturbed KS pair OK).
      CALL RESOLV(IPAIR,2)
*
*       Copy X0DOT since routine RESOLV skipped in KSTERM (IPHASE = -1).
      DO 3 K = 1,3
          X0DOT(K,2*IPAIR-1) = XDOT(K,2*IPAIR-1)
          X0DOT(K,2*IPAIR) = XDOT(K,2*IPAIR)
    3 CONTINUE
*
*       Set variable white dwarf mass (Iben & Renzini, Ann. Rev. 21, 298).
    4 ZMSTAR = 0.38 + 0.15*ZMSTAR
*       Assume neutron star instead if mass > 6 solar masses.
      IF (ZMBAR*BMAX.GT.6.0) ZMSTAR = 1.5
      DM = BODY(IMAX) - ZMSTAR/ZMBAR
      ZMASS = ZMASS - DM
      ZMDOT = ZMDOT + DM*ZMBAR
      NMDOT = NMDOT + 1
      BODY(IMAX) = ZMSTAR/ZMBAR
      IF (TIDAL(1).NE.0.0D0) RTIDE = (ZMASS/TIDAL(1))**0.3333
*
*       Save old velocity for FDOT correction.
      VI2 = 0.0
      DO 5 K = 1,3
          A(K+6) = XDOT(K,IMAX)
          VI2 = VI2 + XDOT(K,IMAX)**2
    5 CONTINUE
*
      DE = 0.5*DM*(VI2 - TIDAL(1)*X(1,IMAX)**2 - TIDAL(3)*X(3,IMAX)**2)
      A0 = 1.0
      IF (ZMBAR*BMAX.LE.6.0) GO TO 8
*
*       Assign a high velocity to neutron star (at least 4*rms velocity).
      A0 = 4.0*SQRT(0.5D0*ZMASS/(RSCALE*VI2))
      IF (IMAX.GT.2*NPAIRS) GO TO 6
      A2 = 2.0*(BODY(2*IPAIR-1) + BODY(2*IPAIR))/R(IPAIR)
*       Add escape velocity from the regularized pair.
      A0 = SQRT(A0**2 + A2/VI2)
    6 DO 7 K = 1,3
          XDOT(K,IMAX) = A0*XDOT(K,IMAX)
          X0DOT(K,IMAX) = XDOT(K,IMAX)
    7 CONTINUE
*
*       Correct total energy, forces & first derivatives.
    8 POTI = 0.0D0
      DO 20 J = 1,NTOT
          IF (J.EQ.IMAX) GO TO 20
          RIJDOT = 0.0D0
          RDVDOT = 0.0D0
*
          DO 10 K = 1,3
              A(K) = X(K,IMAX) - X(K,J)
              A(K+3) = A(K+6) - XDOT(K,J)
              RIJDOT = RIJDOT + A(K)*A(K+3)
              RDVDOT = RDVDOT + A(K)*(XDOT(K,IMAX) - A(K+6))
   10     CONTINUE
*
          RIJ2 = A(1)**2 + A(2)**2 + A(3)**2
          IF (J.LE.N) POTI = POTI + BODY(J)/SQRT(RIJ2)
          IF (J.LE.2*NPAIRS.OR.J.EQ.IBODY) GO TO 20
          A3 = 1.0/(RIJ2*SQRT(RIJ2))
          A4 = BODY(IMAX)*A3
          A5 = DM*A3
          A6 = 3.0*RIJDOT/RIJ2
          A7 = 3.0*RDVDOT/RIJ2
*
          DO 11 K = 1,3
              A(K+3) = (A(K+3) - A(K)*A6)*A5
              IF (A0.GT.1.0) THEN
*       Include FDOT corrections due to increased velocity.
                  A(K+3) = A(K+3) + (XDOT(K,IMAX) - A(K+6))*A4
                  A(K+3) = A(K+3) - A7*A(K)*A4
              END IF
   11     CONTINUE
*
*       Use neighbour list to distinguish irregular & regular terms.
          NNB = LIST(1,J) + 1
          DO 12 L = 2,NNB
              IF (LIST(L,J).EQ.IBODY)  GO TO 16
              IF (LIST(L,J).GT.IBODY)  GO TO 13
   12     CONTINUE
*
   13     DO 14 K = 1,3
              F(K,J) = F(K,J) - 0.5*A(K)*A5
              FR(K,J) = FR(K,J) - A(K)*A5
              FDOT(K,J) = FDOT(K,J) - ONE6*A(K+3)
              D1R(K,J) = D1R(K,J) - A(K+3)
   14     CONTINUE
          GO TO 20
*
   16     DO 18 K = 1,3
              F(K,J) = F(K,J) - 0.5*A(K)*A5
              FI(K,J) = FI(K,J) - A(K)*A5
              FDOT(K,J) = FDOT(K,J) - ONE6*A(K+3)
              D1(K,J) = D1(K,J) - A(K+3)
   18     CONTINUE
   20 CONTINUE
*
*       Correct total energy for mass loss effect.
      DE = DE - DM*POTI
      BE(3) = BE(3) - DE
      RI = SQRT((X(1,IMAX) - RDENS(1))**2 + (X(2,IMAX) - RDENS(2))**2 +
     &                                      (X(3,IMAX) - RDENS(3))**2)
      if(rank.eq.0)
     &WRITE (6,30)  NAME(IMAX), BMAX, BMAX*ZMBAR, ZMDOT, DE, BE(3),
     &              TMDOT, RI, LIST(1,IMAX)
   30 FORMAT (/,'   MASS LOSS   NAME =',I5,'  MI =',F8.4,'  M* =',F5.1,
     &                       '  ZM* =',F6.1,'  DE =',F9.5,'  E =',F10.6,
     &                      '  TMDOT =',F6.1,'  RI =',F5.2,'  NNB =',I3)
*
*       Set option = 2 to skip next energy check (reduced by 1 in ADJUST).
      KZ(19) = 2
*
      IF (ZMBAR*BMAX.LE.6.0) GO TO 50
*       Include correction to total energy from the increased velocity.
      DE = 0.5*BODY(IMAX)*VI2*(A0**2 - 1.0)
      BE(3) = BE(3) + DE
      if(rank.eq.0)
     &WRITE (6,45)  IMAX, DE, TEV1, A0, SQRT(VI2), BE(3)
   45 FORMAT ('   RECOIL   I =',I5,'  DE =',F9.5,'  T* =',F5.1,
     &                      '  V/V0 =',F6.2,'  V0 =',F5.2,'  E =',F10.6)
*
   50 IF (IMAX.LE.2*NPAIRS) THEN
*       Save pair index in KSPAIR and terminate regularization.
          KSPAIR = IPAIR
          IPHASE = -1
*       Indicator for skipping routine RESOLV to preserve velocity X0DOT.
          CALL KSTERM
      ELSE
*       Reduce the steps by velocity factor.
          STEP(IMAX) = MAX(STEP(IMAX)/A0,TIME - T0(IMAX))
          STEPR(IMAX) = MAX(STEPR(IMAX)/A0,STEP(IMAX))
      END IF
*
*       Set next mass loss time & current maximum mass before returning.
      GO TO 99
*
      END
