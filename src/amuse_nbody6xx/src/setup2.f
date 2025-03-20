      SUBROUTINE SETUP2
*
*
*       Initial conditions in astrophysical units (#22 = -1).
*       -----------------------------------------------------
*
      INCLUDE 'common6.h'
*
*
*       Define GM & PC in cgs units and AU in pc (if preferred).
      GM = 6.67D-08*1.989D+33
      PC = 3.0857D+18
      AU = 2.0627E+05
*       Define RBAR (use appropriate value if N-body units skipped in SCALE).
      RBAR = 1.0
*
*       Read initial conditions from unit #10.
      ZMASS = 0.0
      DO 10 I = 1,N
          READ (10,*)  BODY(I), (X(K,I),K=1,3), (XDOT(K,I),K=1,3)
          ZMASS = ZMASS + BODY(I)
   10 CONTINUE
*
*       Form velocity scale in km/sec from mass in Msun & X,Y,Z in pc.
      VSTAR = 1.0D-05*SQRT(GM/PC)
      VSTAR = VSTAR*SQRT(ZMASS/RBAR)
*
*       Convert to N-body units before optional virial theorem scaling.
      DO 20 I = 1,N
          BODY(I) = BODY(I)/ZMASS
          DO 15 K = 1,3
              X(K,I) = X(K,I)/RBAR
              XDOT(K,I) = XDOT(K,I)/VSTAR
   15     CONTINUE
   20 CONTINUE
*       Note: setting KZ(5) = 3 skips scaling to N-body units (ensure RBAR).
*
      RETURN
*
      END
