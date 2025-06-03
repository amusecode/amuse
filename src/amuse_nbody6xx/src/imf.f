      SUBROUTINE IMF(BODY10,BODYN)
*
*
*       Initial mass function.
*       ----------------------
*
      INCLUDE 'common6.h'
*
*
*       Generate Scalo IMF (Eggleton, Fitchett & Tout, Ap.J. 347, 998).
      ITER = 1
*       Set upper limit of mass spectrum & corresponding mass fraction.
      BODYM = BODY10
      XM = 0.998
*       Find initial value of XM for upper & lower mass by iteration.
    1 Y0 = (1.0 - XM)**0.75 + 0.032*(1.0 - XM)**0.25
      Y2 = (0.75 + 0.008/SQRT(1.0 - XM))/
     &                              ((1.0 - XM) + 0.032*SQRT(1.0 - XM))
      Y1 = BODYM - 0.19*XM/Y0
      YPR = -0.19*(1.0 + XM*Y2)/Y0
*       Try next guess of Newton-Raphson iteration.
      XM = XM - Y1/YPR
*       Check upper & lower limits.
      IF (XM.GT.1.0) XM = 0.99999
      IF (XM.LT.0.0) XM = 0.001
*       Set current mass and check the convergence.
      BODYS = 0.19*XM/Y0
      IF (ABS(BODYS - BODYM).GT.1.0D-06*BODYM) GO TO 1
*
*       Save upper value of XM and perform second iteration.
      IF (ITER.EQ.1) THEN
          X1 = XM
          BODYM = BODYN
          XM = 0.4
          ITER = 2
          GO TO 1
      END IF
*
*       Set incremental fraction and assign individual masses.
      DX = (X1 - XM)/FLOAT(N - 1)
      ZMASS = 0.0D0
      DO 10 I = 1,N
          XI = X1 - DX*FLOAT(I-1)
          ZM0 = (1.0 - XI)**0.75 + 0.032*(1.0 - XI)**0.25
          BODY(I) = 0.19*XI/ZM0
          ZMASS = ZMASS + BODY(I)
   10 CONTINUE
*
      if(rank.eq.0)
     &WRITE (6,20)  BODY(1), BODY(N), ZMASS/FLOAT(N)
   20 FORMAT (/,12X,'REALISTIC MASS FUNCTION:','   BODY(1) =',1PE9.2,
     &                                '  BODY(N) =',E9.2,'  <M> =',E9.2)
*
*       Replace input value by actual mean mass in solar units.
      ZMBAR = ZMASS/FLOAT(N)
*
      RETURN
*
      END
