      SUBROUTINE CMCORR
*
*
*       Center of mass & total force corrections.
*       -----------------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Initialize centre of mass variables.
      DO 10 K = 1,3
          CMR(K) = 0.0D0
          CMRDOT(K) = 0.0D0
   10 CONTINUE
*
*       Form c.m. coordinate & velocity displacements.
      DO 20 I = IFIRST,NTOT
          DO 15 K = 1,3
              CMR(K) = CMR(K) + BODY(I)*X(K,I)
              CMRDOT(K) = CMRDOT(K) + BODY(I)*XDOT(K,I)
   15     CONTINUE
   20 CONTINUE
*
      DO 30 K = 1,3
          CMR(K) = CMR(K)/ZMASS
          CMRDOT(K) = CMRDOT(K)/ZMASS
   30 CONTINUE
*
*       Apply c.m. corrections to X & XDOT and accumulate energy changes.
      ERRX = 0.0D0
      ERRV = 0.0D0
      DO 40 I = IFIRST,NTOT
          DO 35 K = 1,3
              XI2 = X(K,I)**2
              VI2 = XDOT(K,I)**2
              X(K,I) = X(K,I) - CMR(K)
              XDOT(K,I) = XDOT(K,I) - CMRDOT(K)
              ERRX = ERRX - TIDAL(K)*BODY(I)*(X(K,I)**2 - XI2)
              ERRV = ERRV + BODY(I)*(XDOT(K,I)**2 - VI2)
   35     CONTINUE
   40 CONTINUE
*
*       Adjust the total energy to new kinetic energy & tidal potential.
      BE(3) = BE(3) + 0.5*(ERRX + ERRV)
      E(11) = E(11) - 0.5*(ERRX + ERRV)
*
*       Perform a consistent shift of the density centre.
      DO 50 K = 1,3
          RDENS(K) = RDENS(K) - CMR(K)
   50 CONTINUE
*
*       Subtract tidal corrections from total force & first derivative.
      IF (KZ(14).GT.0.AND.KZ(14).LE.2) THEN
          DO 60 I = IFIRST,NTOT
*       Skip ghosts to avoid spurious prediction inside 1.0E+10.
              IF (BODY(I).EQ.0.0D0) GO TO 60
              DO 55 K = 1,3
                  DF = TIDAL(K)*CMR(K)
                  DD = TIDAL(K)*CMRDOT(K)
                  FR(K,I) = FR(K,I) - DF
                  F(K,I) = F(K,I) - 0.5*DF
                  D1R(K,I) = D1R(K,I) - DD
                  FRDOT(K,I) = FRDOT(K,I) - DD
                  FDOT(K,I) = FDOT(K,I) - ONE6*DD
                  IF (K.EQ.1) THEN
                      FI(1,I) = FI(1,I) - TIDAL(4)*CMRDOT(2)
                  ELSE IF (K.EQ.2) THEN
                      FI(2,I) = FI(2,I) + TIDAL(4)*CMRDOT(1)
                  END IF
   55         CONTINUE
   60     CONTINUE
      END IF
*
*       Re-determine X0 & X0DOT consistently with current corrected X & XDOT.
      DO 70 I = IFIRST,NTOT 
          IF (BODY(I).EQ.0.0D0) GO TO 70
          DT = TIME - T0(I)
          DTR = TIME - T0R(I)
          A1 = 0.2*DT
          A2 = DT/24.0
*
          DO 65 K = 1,3
              F2DOTK = D3R(K,I)*DTR + D2R(K,I) + D3(K,I)*DT + D2(K,I)
              F3DOTK = D3R(K,I) + D3(K,I)
              DV0 = (((F3DOTK*A2 + ONE6*F2DOTK)*DT +
     &                                3.0*FDOT(K,I))*DT + 2.0*F(K,I))*DT
              X0DOT(K,I) = XDOT(K,I) - DV0
              DX0 = ((((F3DOTK*A1 + F2DOTK)*A2 + FDOT(K,I))*DT +
     &                                       F(K,I))*DT + X0DOT(K,I))*DT
              X0(K,I) = X(K,I) - DX0
              X(K,I) = X0(K,I)
              XDOT(K,I) = X0DOT(K,I)
   65     CONTINUE
   70 CONTINUE
*
*       Ensure consistent coordinates & velocities for binary components.
      DO 80 IPAIR = 1,NPAIRS
          IF (BODY(N+IPAIR).GT.0.0D0) THEN
              CALL RESOLV(IPAIR,1)
          END IF
   80 CONTINUE
*
      RETURN
*
      END
