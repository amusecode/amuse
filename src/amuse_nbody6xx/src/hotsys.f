      SUBROUTINE HOTSYS
*
*
*       Hot initial system.
*       -------------------
*
      INCLUDE 'common6.h'
*
*
*       Determine the rms velocity from current kinetic energy.
      VRMS = 0.0
      DO 10 I = 1,N
          DO 5 K = 1,3
              VRMS = VRMS + BODY(I)*XDOT(K,I)**2
    5     CONTINUE
   10 CONTINUE
      VRMS = SQRT(VRMS/ZMASS)
*
*       Define GM & PC in cgs units and form velocity scale in km/sec.
      GM = 6.67E-08*1.989E+33
      PC = 3.0857E+18
      VSTAR = 1.0E-05*SQRT(GM/PC)
*
*       Ensure ZMBAR & RBAR > 0 (=0: assume <M>/Sun = 1, RBAR = 1 pc).
      IF (ZMBAR.LE.0.0D0) ZMBAR = FLOAT(N)/ZMASS
      IF (RBAR.LE.0.0D0) RBAR = 1.0
*
*       Scale to working units of RBAR in pc & ZMBAR in solar masses.
      VSTAR = VSTAR*SQRT(ZMASS*ZMBAR/RBAR)
*
*       Read central velocity dispersion and form scaling factor.
      READ (5,*)  SIGMA0
      VSCALE = SIGMA0/(VSTAR*VRMS)
*
*       Scale the velocities to central velocity dispersion of SIGMA0 km/sec.
      DO 20 I = 1,N
          DO 15 K = 1,3
              XDOT(K,I) = VSCALE*XDOT(K,I)
              X0DOT(K,I) = XDOT(K,I)
   15     CONTINUE
   20 CONTINUE
*
*       Rescale crossing time, output times & termination time.
      RATIO = SIGMA0/VSTAR
      TCR = TCR/RATIO
      TCR0 = TCR
      DTADJ = DTADJ/RATIO
      DELTAT = DELTAT/RATIO
      TCRIT = TCRIT/RATIO
      VC = VSCALE*VC
*
      WRITE (6,30)  SIGMA0, VRMS, VSCALE
   30 FORMAT (/,12X,'HOT SYSTEM:   SIGMA0 =',F5.1,'  <V> =',F6.3,
     &                                              '  VSCALE =',F6.3,/)
*
      RETURN
*
      END
