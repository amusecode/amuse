      REAL*8 FUNCTION RPMAX2(R1,R2,M1,M2,KS1,KS2,VINF)
*
*
*       Maximum periastron factor for capture.
*       --------------------------------------
*
*       Developed by Rosemary Mardling (March 1996).
*       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*
      REAL*8  R1,R2,M1,M2,VINF,M21,A(8),B(2),P(2)
      INTEGER KS1,KS2,KS(2)
      DATA  A /-600.2227,3538.748,1703.075,77.31540,
     &         191.8260,-2.203482,1202.116,70.80281/
      DATA B /0.54,1.08/
*
*
*       Define mass and radius ratios and swap stellar type if necessary.
      IF (R1.GE.R2) THEN
          M21 = M2/M1
          R21 = R2/R1
          KS(1) = KS1
	  KS(2) = KS2
      ELSE
	  M21 = M1/M2
	  R21 = R1/R2
	  KS(1) = KS2
	  KS(2) = KS1
      END IF
*
*       Choose polytropic index for each star (adopt n=1.5 for giants).
      DO 1 I = 1,2
          IF(KS(I).EQ.0.OR.KS(I).EQ.3.OR.KS(I).EQ.5.OR.
     &       KS(I).EQ.6.OR.KS(I).EQ.9)THEN
              IP = 1
              PP = B(1)
          ELSE
              IP = 5
              PP = B(2)
          ENDIF
          P(I) = (((A(IP+3)*M21 + A(IP+2))*M21 + A(IP+1))*M21 
     &                                         + A(IP))/VINF**PP
    1 CONTINUE
*
*       Calculate maximum periastron separation for capture (note arg < 0).
      Z = P(1) + (R21**5/M21**2)*P(2)
      IF (Z.LE.0.0) THEN
*       Adopt factor 2 for collision test in case of small mass ratio.
          RPMAX2 = 2.0
      ELSE
          RPMAX2 = Z**(1.0/6.0)
*         RPMAX2 = (P(1) + (R21**5/M21**2)*P(2))**(1.0/6.0)
      END IF
*
      RETURN
*
      END
