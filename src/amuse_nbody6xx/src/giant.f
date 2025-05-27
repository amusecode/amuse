      SUBROUTINE GIANT(IPAIR,I,W,Q,WSCALE,QSCALE,ZN,QL)
*
*
*       Structure constants of giant star.
*       ----------------------------------
*
*       Theory of Rosemary Mardling, Ap. J. XX, YYY, 1995.
*       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*
      INCLUDE 'common6.h'
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &               BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &               RP(NTMAX),ES(NTMAX),CM(2,NTMAX),IOSC(NTMAX),
     &               NAMEC(NTMAX)
      REAL*8  WW(6),QQ(6),W(2),Q(2),WSCALE(2),QSCALE(2),SW(2)
      DATA  WW  /2.119,3.113,8.175,3.742,4.953,9.413/
      DATA  QQ  /0.4909,0.4219,0.2372,0.4677,0.3560,0.1519/
      DATA  A0,A1,A2,A3  /0.944525,-0.392030,6.01655D-02,-3.34790D-03/
      DATA  B0,B1,B2  /-0.3789,1.481,-0.1018/
      DATA  C0,C1,C2,C3  /1.452104,3.923872,-11.88722,13.46106/
      DATA  E0,E1,E2,E3  /1.934977,2.214222,-4.855796,4.025394/
*
*
*       Determine the appropriate chaos index (adopt NCHAOS+1 if not found).
      IC = 0
      J = N + IPAIR
      DO 1 K = 1,NCHAOS
          IF (NAMEC(K).EQ.NAME(J)) IC = K
    1 CONTINUE
*
*       Decide between first & second binary component for saving core mass.
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      L = 1
*       Note possibility I = I2 on first call.
      IF (RADIUS(I1).GE.RADIUS(I2)) THEN
          IF (I.EQ.I2) L = 2
      END IF
*
*       Set typical core mass of 0.3/0.5 Msun for general binary in next bin.
      IF (IC.EQ.0) THEN
          IC = NCHAOS + 1
          CM(L,IC) = (0.3 + 0.1*FLOAT(KSTAR(I) - 3))/ZMBAR
*       Include rare case of hyperbolic encounter as the first event.
      ELSE IF (CM(L,IC).LE.0.0D0) THEN
          CM(L,IC) = (0.3 + 0.1*FLOAT(KSTAR(I) - 3))/ZMBAR
      END IF
      IF(KSTAR(I).EQ.9) CM(L,IC) = 0.5*BODY(I)
*
*       Form ratio of core mass and mass.
      SIG = CM(L,IC)/BODY(I)
*
*       Include safety check on mass ratio just in case.
      IF (SIG.GT.0.9.OR.SIG.LT.0.0) THEN
*         WRITE (6,5)  IC, KSTAR(I), NAME(I), BODY0(I)*ZMBAR,
*    &                 BODY(I)*ZMBAR, CM(L,IC)*ZMBAR
*   5     FORMAT (' WARNING!    GIANT    IC K* NM M0 M MC ',
*    &                                   2I4,I6,3F7.2)
          SIG = 0.9
          CM(L,IC) = 0.9*BODY(I)
      END IF
*
*       Define mass, core mass, radius, envelope mass and luminosity in S.U.
      ZM = BODY(I)*ZMBAR
      ZMC = CM(L,IC)*ZMBAR
      RSI = RADIUS(I)*SU
      ZME = ZM - ZMC
*       Obtain L^{1/3} from GB L(Mc) relation for Pop II M = 0.8Msun.
      ZL = 1.98D+05*ZMC**6
*       Evaluate damping constant from Zahn's theory (R.M. 13/5/97).
*       FAC = (GM)^{1/2}*M_{env}^{1/3)/(R^{5/6}*L^{1/3}) = 8.48D+03 for S.U.
      QL = 8.48D+03*SQRT(ZM)*(ZME/ZL)**0.33/RSI**0.833
*
*       Set effective frequencies, overlap integrals and structure constants.
      DO 10 K = 1,2
          K1 = 3*K - 2
	  IF (K.EQ.1) THEN
              SW(K) = ((C3*SIG + C2)*SIG + C1)*SIG + C0
          ELSE
              SW(K) = ((E3*SIG + E2)*SIG + E1)*SIG + E0
	  END IF
	  W(K) = SW(K)**2
          Q(K) = ((A3*SW(K) + A2)*SW(K) + A1)*SW(K) + A0
          WSCALE(K) = SQRT(W(K)/WW(K1))
          QSCALE(K) = (Q(K)/QQ(K1)/WSCALE(K))**2
          QSCALE(K) = MAX(QSCALE(K),0.0001D0)
   10 CONTINUE
*
*       Evaluate new polytropic index.
      ZN = (B2*SW(1) + B1)*SW(1) + B0
*     WRITE (24,20)  IC, IPAIR, KSTAR(J), CM(L,IC)/BODY(I), RSI, ZN, QL
*  20 FORMAT (' GIANT:    IC KS K* MC/M R* n Q ',3I4,3F6.2,F7.1)
*     CALL FLUSH(24)
*
*       Include warning if n > 5 (note limit QSCALE >= 0.0001).
*     QSM = MIN(QSCALE(1),QSCALE(2))
*     IF (ZN.GE.5.0.OR.QSM.LT.0.00011) THEN
*         WRITE (6,30)  IC, IPAIR, KSTAR(J), CM(L,IC)/BODY(I),
*    &                  RADIUS(I)*SU, QSM, ZN, QL
*  30     FORMAT (' GIANT:    WARNING!    IC KS K* MC/M R/R0 QSM n QL ',
*    &                                    3I4,F6.2,F7.1,F8.4,F6.2,F7.1)
*     END IF

      RETURN
*
      END
