***
      SUBROUTINE MIX(J1,J2,DM)
*
*     Author : J. R. Hurley
*     Date :   7th July 1998
*
*       Evolution parameters for mixed star.
*       ------------------------------------
*
      INCLUDE 'common6.h'
*
      REAL*8 TSCLS(20),LUMS(10),GB(10),TMS1,TMS2,TMS3,TN
      REAL*8 M01,M02,M03,M1,M2,M3,AGE1,AGE2,AGE3,MC3,LUM,RM,RCC
      REAL*8 MENV,RENV,K2E
      REAL*8 MCH,MXNS
      PARAMETER(MCH=1.44D0,MXNS = 3.0d0)
      LOGICAL  FIRST
      SAVE  FIRST
      DATA  FIRST /.TRUE./
      DATA  rg2 /0.1d0/
*
*
*       Define global indices with body #I1 being most evolved.
      IF(KSTAR(J1).GE.KSTAR(J2))THEN
          I1 = J1
          I2 = J2
          IF(KSTAR(J1).LE.1.AND.BODY(J1).LT.BODY(J2))THEN
             I1 = J2
             I2 = J1
          ENDIF
      ELSE
          I1 = J2
          I2 = J1
      END IF
*
*       Specify case index for collision treatment.
      K1 = KSTAR(I1)
      K2 = KSTAR(I2)
      ICASE = KTYPE(K1,K2)
      IF(rank.eq.0.and.ICASE.GT.100)THEN
         WRITE(38,*)' MIX ERROR ICASE>100 ',ICASE,K1,K2
      ENDIF
*
*       Record name and stellar type of the colliding stars.
      JC = NCOLL + 1
      ITYPE(1) = NAME(I1)
      ITYPE(2) = NAME(I2)
      ITYPE(3) = KSTAR(I1)
      ITYPE(4) = KSTAR(I2)
*
*       Set physical time and initialize mass loss & time increment.
      TTOT = TIME + TOFF
      TPHYS = TTOT*TSTAR
      DM = 0.D0
      DT = 0.D0
      RS1 = RADIUS(I1)*SU
      RS2 = RADIUS(I2)*SU
*
*       Evolve the stars to the current time unless they have been
*       evolved further by recent Roche interaction.
*     TEV1 = MAX(TIME,TEV0(I1))
*       Update the stars to the latest previous time.
      TEV1 = MAX(TEV0(I1),TEV0(I2))
*
*       Determine evolution time scales for first star.
      M01 = BODY0(I1)*ZMBAR
      M1 = BODY(I1)*ZMBAR
      AGE1 = TEV1*TSTAR - EPOCH(I1)
      CALL star(K1,M01,M1,TMS1,TN,TSCLS,LUMS,GB,ZPARS)
*
*       Obtain time scales for second star.
      M02 = BODY0(I2)*ZMBAR
      M2 = BODY(I2)*ZMBAR
      AGE2 = TEV1*TSTAR - EPOCH(I2)
      CALL star(K2,M02,M2,TMS2,TN,TSCLS,LUMS,GB,ZPARS)
*
*       Check for planetary systems - defined as HeWDs - and low-mass WDs!
      IF(K1.EQ.10.AND.M1.LT.0.05)THEN
         ICASE = K2
         IF(K2.LE.1)THEN
            ICASE = 1
            AGE1 = 0.D0
         ENDIF
      ELSEIF(K1.GE.11.AND.M1.LT.0.15.AND.ICASE.EQ.6)THEN
         ICASE = 1
      ENDIF
      IF(K2.EQ.10.AND.M2.LT.0.05)THEN
         ICASE = K1
         IF(K1.LE.1)THEN
            ICASE = 1
            AGE2 = 0.D0
         ENDIF
      ENDIF
*
      if(rank.eq.0)WRITE(38,67)K1,M01,M1
 67   FORMAT(' MIX OLD *1:',I4,2F10.6)
      if(rank.eq.0)WRITE(38,68)K2,M02,M2
 68   FORMAT(' MIX OLD *2:',I4,2F10.6)
*
*       Specify total mass.
      M3 = M1 + M2
      M03 = M01 + M02
      MC3 = 0.d0
      KW = ICASE
      AGE3 = 0.d0
      TMS3 = 0.d0
*
*       Check energy budget for partial disruption.
      ZMB = BODY(I1) + BODY(I2)
      RIJ2 = 0.0
      VIJ2 = 0.0
      RDOT = 0.0
      DO 50 K = 1,3
          RIJ2 = RIJ2 + (X(K,I1) - X(K,I2))**2
          VIJ2 = VIJ2 + (XDOT(K,I1) - XDOT(K,I2))**2
          RDOT = RDOT + (X(K,I1) - X(K,I2))*(XDOT(K,I1) - XDOT(K,I2))
   50 CONTINUE
      RIJ = SQRT(RIJ2)
*      Form binding energy of secondary and orbital kinetic energy.
      EB = -0.5*6.68D-08*4.0D+66*M2**2/(RADIUS(I2)*SU*7.0D+10)
      ZMU = SMU*BODY(I1)*BODY(I2)/(BODY(I1) + BODY(I2))
      ZKE = 0.5*2.0D+33*ZMU*1.0D+10*VIJ2*VSTAR**2
      IF (ZKE + EB.GT.0.0.AND.RIJ.LT.0.5*(RADIUS(I1)+RADIUS(I2)).AND.
     &    MAX(K1,K2).LT.2) THEN
          ZM = ABS(ZKE/EB)
          FAC = 1.0/SQRT(ZM)
       if(rank.eq.0)WRITE (6,55)  FAC, RIJ*SU, SQRT(VIJ2)*VSTAR, EB, ZKE
   55     FORMAT (' ENFORCED MASS LOSS    FAC RIJ VIJ EB KE ',
     &                                    F7.3,F8.2,F8.2,1P,2E10.2)
*         DM = FAC*M2
          DM = 0.3D0*M2
          M3 = M3 - DM
          M03 = M3
          M02 = M02 - DM
      END IF
*
*       Evaluate apparent age and other parameters.
*
      IF(ICASE.EQ.1)THEN
*       Specify new age based on complete mixing.
         IF(K1.EQ.7) KW = 7
         CALL star(KW,M03,M3,TMS3,TN,TSCLS,LUMS,GB,ZPARS)
         AGE3 = 0.1d0*TMS3*(AGE1*M01/TMS1 + AGE2*M02/TMS2)/M03
      ELSEIF(ICASE.EQ.3.OR.ICASE.EQ.6.OR.ICASE.EQ.9)THEN
         MC3 = M1
         CALL gntage(MC3,M3,KW,ZPARS,M03,AGE3)
      ELSEIF(ICASE.EQ.4)THEN
         MC3 = M1
         AGE3 = AGE1/TMS1
         IF(AGE3.GT.1.D0)THEN
            if(rank.eq.0)then
              WRITE(6,*)' WARNING! MIX AGE WRONG FOR KW=4 ',age1,tms1
              WRITE(6,*)K1,M01,M1,TEV1,EPOCH(I1)
            end if
            AGE3 = 0.99D0
         ENDIF
         CALL gntage(MC3,M3,KW,ZPARS,M03,AGE3)
      ELSEIF(ICASE.EQ.7)THEN
         CALL star(KW,M03,M3,TMS3,TN,TSCLS,LUMS,GB,ZPARS)
         AGE3 = TMS3*(AGE2*M2/TMS2)/M3
      ELSEIF(ICASE.LE.12)THEN
*       Ensure that a new WD has the initial mass set correctly.
         M03 = M3
         DT = 1.0d+04
      ELSEIF(ICASE.EQ.13.OR.ICASE.EQ.14)THEN
*       Set unstable Thorne-Zytkow object with fast mass loss of envelope
*       unless the less evolved star is a WD, NS or BH.
         IF(K2.LT.10)THEN
            M03 = M1
            M3 = M1
            DM = M2
            NTZ = NTZ + 1
            if(rank.eq.0)WRITE (6,2)  K1, K2, NAME(I1), NAME(I2)
    2       FORMAT (' NEW TZ    K* NM ',2I4,2I6)
         ELSEIF(ICASE.EQ.13)THEN
            NGB = NGB + 1
            IF(M3.GT.MXNS) KW = 14
         ENDIF
         DT = 1.0D+04
      ELSEIF(ICASE.EQ.15)THEN
         DM = M3
         M3 = 0.D0
      ELSEIF(ICASE.GT.100)THEN
*       Common envelope case which should only be used after COMENV.
         KW = K1
         AGE3 = AGE1
         M3 = M1
         M03 = M01
         DM = M2
      ELSE
*       This should not be reached.
        if(rank.eq.0)WRITE(6,*)' ERROR MIX: ICASE NOT CAUGHT!!!'
        if(rank.eq.0)WRITE(6,*)' K1 K2 -> K3 ',K2,K2,KW
        STOP
      ENDIF
*
      if(rank.eq.0)WRITE(38,69)KW,M03,M3
 69   FORMAT(' MIX NEW *3:',I4,2F10.6)
*
*       Determine consistent stellar type and specify mass loss.
      IF(KW.LE.14)THEN
         KW0 = KW 
         CALL star(KW,M03,M3,TMS3,TN,TSCLS,LUMS,GB,ZPARS)
         CALL hrdiag(M03,AGE3,M3,TMS3,TN,TSCLS,LUMS,GB,ZPARS,
     &               RM,LUM,KW,MC3,RCC,MENV,RENV,K2E)
         IF(KW.NE.KW0)then
            DM = M1 + M2 - M3
            if(rank.eq.0)WRITE (6,5)  NAME(I1), NAME(I2), KW0, KW, DM
    5       FORMAT (' MIX TYPE CHANGE:  NAM K* DM ',2I6,2I4,F6.2)
         ENDIF
      ENDIF
      KSTAR(I1) = KW
      RADIUS(I1) = RM/SU
      ZLMSTY(I1) = LUM
*       Update initial mass, epoch & evolution times.
      BODY0(I1) = M03/ZMBAR
      EPOCH(I1) = TEV1*TSTAR - AGE3
      TEV(I1) = TEV1 + DT
      TEV0(I1) = TEV1
      DM = DM/ZMBAR
*
*       Record final type of mixed star.
      ITYPE(5) = KSTAR(I1)
*
*       Check for blue straggler formation (TM < TPHYS & KSTAR <= 1).
      IF(TMS3.LT.TPHYS.AND.KSTAR(I1).LE.1)THEN
          if(rank.eq.0)WRITE (6,8)  NAME(I1), M3, TMS3, TPHYS, AGE3
    8     FORMAT (' NEW BS (MIX):    NAM M3 TM TP AGE ',
     &                               I6,F6.2,3F8.1)
*       Prepare spin diagnostics for fort.91 (rg2=0.1 is assumed).
          SEMI = 2.0/RIJ - VIJ2/ZMB
          SEMI = 1.0/SEMI
          PD = DAYS*SEMI*SQRT(ABS(SEMI)/ZMB)
          PD = MAX(PD,0.0D0)
          ECC2 = (1.0 - RIJ/SEMI)**2 + RDOT**2/(SEMI*ZMB)
          ECC = SQRT(ECC2)
          SPIN1 = SPIN(I1)/(rg2*BODY(I1)*RADIUS(I1)**2)
          SPIN2 = SPIN(I2)/(rg2*BODY(I2)*RADIUS(I2)**2)
          TS = 1.0D+06*365.0*TWOPI*TSTAR
          NBS = NBS + 1
          if(rank.eq.0)
     &        WRITE (91,9)  TIME+TTOT, NAME(I1), NAME(I2), M3, ECC, PD,
     &                  TS/SPIN1, TS/SPIN2
    9     FORMAT (' NEW BS    T NAM M3 ECC P ROT ',
     &                        F8.1,2I7,F6.2,F8.4,1P,3E9.1)
          CALL FLUSH(91)
      ENDIF
*
      if(rank.eq.0)
     &WRITE (6,10)  IQCOLL, K1, K2, KSTAR(I1), M1, M2, M3, RS1, RS2,
     &              RADIUS(I1)*SU, DM*ZMBAR
   10 FORMAT (' NEW STAR:    IQ K1 K2 K* M1 M2 M3 R1 R2 R* DM ',
     &                       4I3,3F6.1,3F7.1,F5.1)
*
*       Open unit #13 the first time.
      IF(rank.eq.0.and.FIRST.AND.IQCOLL.NE.3)THEN
          OPEN (UNIT=13,STATUS='UNKNOWN',FORM='FORMATTED',FILE='COLL')
          FIRST = .FALSE.
*
*       Print cluster scaling parameters at start of the run.
          IF(NCOLL.EQ.0)THEN
              WRITE (13,20)  RBAR, BODYM*ZMBAR, BODY1*ZMBAR, TSCALE,
     &                       NBIN0, NZERO
   20         FORMAT (/,6X,'MODEL:    RBAR =',F5.1,'  <M> =',F6.2,
     &                     '  M1 =',F6.1,'  TSCALE =',F6.2,
     &                     '  NB =',I4,'  N0 =',I6,//)
              WRITE (13,25)
   25         FORMAT ('    TIME  NAME  NAME  K1  K2  KC  M1   M2   MC',
     &                '   DM    R1     R2    r/Rc   R     ECC      P',/)
          ENDIF
      ENDIF
*
*       Form central distance (scaled by RC) and period in days.
      RI2 = 0.d0
      RIJ2 = 0.d0
      VIJ2 = 0.d0
      DO 30 K = 1,3
          RI2 = RI2 + (X(K,I1) - RDENS(K))**2
          RIJ2 = RIJ2 + (X(K,I1) - X(K,I2))**2
          VIJ2 = VIJ2 + (XDOT(K,I1) - XDOT(K,I2))**2
   30 CONTINUE
      RI = SQRT(RI2)/RC
      RIJ = SQRT(RIJ2)
      SEMI = 2.d0/RIJ - VIJ2/ZMB
      SEMI = 1.d0/SEMI
      TK = DAYS*SEMI*SQRT(ABS(SEMI)/ZMB)
      ECC = 1.0 - RIJ/SEMI
*
*       Accumulate collision diagnostics on unit #13 (filename COLL).
      IF (rank.eq.0.and.IQCOLL.NE.3) THEN
          WRITE (13,35)  TTOT, (ITYPE(K),K=1,5), M1, M2, M3,
     &                   DM*ZMBAR, RS1, RS2, RI/RC, RIJ*SU, ECC, TK
   35     FORMAT (1X,F7.1,2I6,3I4,4F5.1,2F7.2,F6.2,F7.2,F9.5,1P,E9.1)
          CALL FLUSH(13)
      END IF
*
*       Re-define indices of colliding bodies with J1 as new c.m.
      J1 = I1
      J2 = I2
*
      IF(rank.eq.0.and.KSTAR(I1).GT.12)THEN
          WRITE (15,40)  K1, K2, KSTAR(I1), BODY(I1)*ZMBAR,
     &                   BODY(I2)*ZMBAR, M3
   40     FORMAT (' MIX:    K1 K2 K* M1 M2 M3 ',3I4,3F7.2)
          CALL FLUSH(15)
      ENDIF
*
      RETURN
      END
***
