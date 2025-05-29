      SUBROUTINE GIANT3(I,BODYI,W,Q,ZN,QL)
*
*
*       Structure constants of giant star.
*       ----------------------------------
*
*       Theory of Rosemary Mardling, Ap. J. XX, YYY, 1995.
*       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*
      INCLUDE 'common6.h'
      REAL*8  W(2),Q(2),SW(2),MC
*     DATA  WW  /2.119,3.113,8.175,3.742,4.953,9.413/
*     DATA  QQ  /0.4909,0.4219,0.2372,0.4677,0.3560,0.1519/
      DATA  A0,A1,A2,A3  /0.944525,-0.392030,6.01655D-02,-3.34790D-03/
      DATA  B0,B1,B2  /-0.3789,1.481,-0.1018/
      DATA  C0,C1,C2,C3  /1.452104,3.923872,-11.88722,13.46106/
      DATA  E0,E1,E2,E3  /1.934977,2.214222,-4.855796,4.025394/
*
*
*       Set typical core mass of 0.3/0.5 Msun.
      MC = 0.3 + 0.1*FLOAT(KSTAR(I) - 3)
      IF(KSTAR(I).EQ.9) MC = 0.5*BODYI*ZMBAR
*
*       Form ratio of core mass and mass.
      SIG = MC/(BODYI*ZMBAR)
*
*       Include safety check on mass ratio just in case.
      IF (SIG.GT.0.9) SIG = 0.9
*
*       Define mass, radius, envelope mass and luminosity in S.U.
      ZM = BODYI*ZMBAR
      RSI = RADIUS(I)*SU
      ZME = ZM - MC
*       Obtain L^{1/3} from giant relation L = 1.98D+05*M_c^6.
      ZL3 = 58.3*MC**2
*       Evaluate damping constant from Zahn's theory (R.M. 13/5/97).
*       FAC = (GM)^{1/2}*M_{env}^{1/3)/(R^{5/6}*L^{1/3}) = 8.48D+03 for S.U.
      QL = 8.48D+03*SQRT(ZM)*(ZME/ZL3)**0.33/RSI**0.833
*
*       Set effective frequencies, overlap integrals and structure constants.
      DO 10 K = 1,2
*         K1 = 3*K - 2
	  IF (K.EQ.1) THEN
              SW(K) = ((C3*SIG + C2)*SIG + C1)*SIG + C0
          ELSE
              SW(K) = ((E3*SIG + E2)*SIG + E1)*SIG + E0
	  END IF
	  W(K) = SW(K)**2
          Q(K) = ((A3*SW(K) + A2)*SW(K) + A1)*SW(K) + A0
   10 CONTINUE
*
*       Evaluate new polytropic index.
      ZN = (B2*SW(1) + B1)*SW(1) + B0
*     WRITE (24,20)  IC, IPAIR, KSTAR(J), MC/(BODYI*SMU), RSI, ZN, QL
*  20 FORMAT (' GIANT:    IC KS K* MC/M R* n Q ',3I4,3F6.2,F7.1)
*     CALL FLUSH(24)
*
*       Include warning if n > 5.
*     IF (ZN.GE.5.0) THEN
*         WRITE (6,30)  IC, IPAIR, KSTAR(J), CM(L,IC)/BODY(I),
*    &                  RADIUS(I)/RIN, ZN, QL
*  30     FORMAT (' GIANT:    WARNING!    IC KS K* MC/M R/R0 n QL ',
*    &                                    3I4,F6.2,F6.1,F6.2,F7.1)
*     END IF

      RETURN
*
      END
