***
      SUBROUTINE ROCHE(IPAIR)
*
*
*       Roche lobe overflow.
*       Updated 6/1/98 & 28/11/01 by J. Hurley.
*       ---------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/RCHE/ CMRCH(13,MMAX),NAMER(2,MMAX),KSTARR(3,MMAX)
      REAL*8 TSCLS(20),LUMS(10),GB(10),TM,TN,TPHYS0
      REAL*8 M01,M1,AGE,R1,LUM,MC,RCC,SEP,M1CE,M2CE,RCE,KM0,KM,M2,MCX
      REAL*8 MASS0(2),MASS(2),MASSC(2),RAD(2),RADX(2),ROL(2),
     &       RADC(2),Q(2),AJ(2),TKH(2),TMS(2),TBGB(2),
     &       JSPIN(2),OSPIN(2),DJSPIN(2),DTMI(2)
      REAL*8 MENV(2),RENV(2),K2STR(2)
      REAL*8 JORB,OORB,RLPERI,DJORB,DJMB,DJGR,DJT,RDMIN,RDISK
      REAL*8 JSPBRU,OSPBRU
      REAL*8 DMS(2),DMA(2),DMR(2),IVSQM,VWIND2,VORB2
      REAL*8 DM1,DM2,DM22,DSPIN,DSPIN0,DELET,DELET1,EQSPIN,MEW,MT2
      REAL*8 TC,TSYN,TTID,FAC,FAC0,ECC0,DTMMIN,DJTK(2),DSPINK(2)
      REAL*8 MCH,EPSNOV,EDDFAC,GAMM1
      PARAMETER(MCH=1.44D0,EPSNOV=0.D0,EDDFAC=100.D0,GAMM1=-1.D0)
      REAL*8 K2,K3,ACC1,ACC2,BETA,XI,AURSUN
      PARAMETER(K3=0.21D0,ACC1=3.920659D+08,ACC2=1.5D0)
      PARAMETER(BETA=0.125D0,XI=1.D0,AURSUN=214.95D0)
      REAL*8 MLWIND,RL
      EXTERNAL MLWIND,RL
      LOGICAL COALS,DISK,NOVAE,SUPEDD,ISAVE
      CHARACTER*5 CH5
      SAVE IWARN,IGR,IMB
      DATA IWARN,IGR,IMB /0,0,0/
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./
*
*
*       Set components & c.m. index and initialize counters & interval.
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      I = N + IPAIR
      ITER = 0
      IQCOLL = 0
      DT1 = 0.0
      ISAVE = .FALSE.
      JKICK = 0
      INEW = 0
*
*       Form semi-major axis and eccentricity.
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      ECC2 = (1.D0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(SEMI*BODY(I))
      ECC = SQRT(ECC2)
*
*       Set age, mass ratio, Roche radius & radius (SU) for each star.
      J = I1
      DO 5 K = 1,2
          Q(K) = BODY(J)/(BODY(I) - BODY(J))
          ROL(K) = RL(Q(K))*SEMI*SU
          RAD(K) = RADIUS(J)*SU
          DMS(K) = 0.D0
          J = I2
    5 CONTINUE
*
*       Determine indices for primary & secondary star (donor & accretor).
      IF (RAD(1)/ROL(1).GE.RAD(2)/ROL(2)) THEN
          J1 = I1
          J2 = I2
          RL1 = ROL(1)
      ELSE
          J1 = I2
          J2 = I1
          RL1 = ROL(2)
          RR = RAD(1)
          RAD(1) = RAD(2)
          RAD(2) = RR
          RR = ROL(1)
          ROL(1) = ROL(2)
          ROL(2) = RR
          RR = Q(1)
          Q(1) = Q(2)
          Q(2) = RR
      END IF
*
*       Exit if physical radius is smaller than Roche radius.
      IF (RAD(1).LE.RL1) THEN
          GO TO 200
      ELSE
*       Set the time to the latest evolution time of the components.
          TPHYS0 = MAX(TEV0(J1),TEV0(J2))*TSTAR
          TPHYS = TPHYS0 + TOFF*TSTAR
*       Save the original stellar types and initialize coalescence flag.
          KW1 = KSTAR(J1)
          KW2 = KSTAR(J2)
          COALS = .FALSE.
*         if(rank.eq.0)
*    &    write(6,*)' start roche ',name(j1),name(j2),kw1,kw2,
*    &                              rad(1)/rl1,rad(2)/rol(2)
*
*       Set separation & masses (solar units) and mass ratios.
          SEP = SEMI*SU
          SEMI0 = SEMI
          MASS(1) = BODY(J1)*ZMBAR
          MASS(2) = BODY(J2)*ZMBAR
          MASS0(1) = BODY0(J1)*ZMBAR
          MASS0(2) = BODY0(J2)*ZMBAR
*
*       Set the orbital spin and angular momentum.
          TB = YRS*SEMI*SQRT(SEMI/BODY(I))
          TK = DAYS*TB/YRS
          OORB = TWOPI/TB
          JORB = MASS(1)*MASS(2)/(MASS(1)+MASS(2))
     &           *SQRT(1.D0-ECC*ECC)*SEP*SEP*OORB
*
          DM1 = 0.D0
          DM2 = 0.D0
*
*       Initialize local arrays.
          J = J1
          DO 97 K = 1,2
             M01 = MASS0(K)
             M1 = MASS(K)
             MC = 0.D0
             KW = KSTAR(J)
             AGE = TEV0(J)*TSTAR - EPOCH(J)
             CALL star(KW,M01,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
             CALL hrdiag(M01,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                   R1,LUM,KW,MC,RCC,M1CE,RCE,K2)
             MASSC(K) = MC
             AJ(K) = AGE
             MENV(K) = M1CE
             RENV(K) = RCE
             K2STR(K) = K2
             RADX(K) = R1
             IF(K.EQ.1)THEN
                RADX(K) = MIN(RADX(K),ROL(K))
                RADX(K) = MAX(RCC,RADX(K))
             ENDIF
             JSPIN(K) = MAX(SPIN(J)*SPNFAC,1.0D-10)
             OSPIN(K) = JSPIN(K)/(K2*RADX(K)*RADX(K)*(M1-MC) + 
     &                  K3*RCC*RCC*MC)
*            IF(K.EQ.1)THEN
*               IF(OSPIN(1).GT.OORB)THEN
*       Force co-rotation of primary and orbit to ensure that the
*       tides do not lead to unstable Roche.
*                  OSPIN(1) = OORB
*                  JSPIN(1) = OSPIN(1)*(K2*R1*R1*(M1-MC)+K3*RCC*RCC*MC)
*               ENDIF
*            ENDIF
             RLPERI = ROL(K)*(1.D0 - ECC)
             DMR(K) = MLWIND(KW,LUM,R1,M1,MC,RLPERI,ZMET)
             DMA(K) = 0.D0
             TMS(K) = TM
             TBGB(K) = TSCLS(1)
             J = J2
   97     CONTINUE
*
*       Implement any mass loss between TEV0 and the Roche starting time.
*       This is a check for the case where one star is catching up on the
*       other and Roche occurs in the meantime!
          J = J1
          VORB2 = ACC1*(MASS(1)+MASS(2))/SEP
          IVSQM = 1.D0/SQRT(1.D0-ECC*ECC)
          DO 98 K = 1,2
             DT = 1.0D+06*(TPHYS0 - TEV0(J)*TSTAR)
             IF(DT.EQ.0.0) GOTO 99
*            VWIND2 = 2.D0*BETA*ACC1*MASS(3-K)/RAD(3-K)
*            OMV2 = (1.D0 + VORB2/VWIND2)**(3.D0/2.D0)
*            DMA(K) = IVSQM*ACC2*DMR(3-K)*((ACC1*MASS(K)/VWIND2)**2)/
*    &                (2.D0*SEP*SEP*OMV2)
*            DMA(K) = MIN(DMA(K),0.8D0*DMR(3-K))
             DMS(K) = (DMR(K) - DMA(K))*DT
             DJSPIN(K) = (2.D0/3.D0)*(DMR(K)*RADX(K)*RADX(K)*OSPIN(K) -
     &                    XI*DMA(K)*RADX(3-K)*RADX(3-K)*OSPIN(3-K))*DT
             KW = KSTAR(J)
             IF(KW.LT.10)THEN
                DMS(K) = MIN(DMS(K),MASS(K)-MASSC(K))
                CALL MAGBRK(KW,MASS(K),MENV(K),RAD(K),OSPIN(K),DJMB)
                DJSPIN(K) = DJSPIN(K) + DJMB*DT
             ENDIF
             MASS(K) = MASS(K) - DMS(K)
             JSPIN(K) = MAX(JSPIN(K) - DJSPIN(K),1.0D-10)
             IF(KW.LE.2.OR.KW.EQ.7)THEN
                M01 = MASS0(K)
                MASS0(K) = MASS(K)
                BODY0(J) = MASS(K)/ZMBAR
                CALL star(KW,MASS0(K),MASS(K),TM,TN,TSCLS,LUMS,GB,ZPARS)
                IF(KW.EQ.2)THEN
                   IF(GB(9).LT.MASSC(K).OR.M01.GT.ZPARS(3))THEN
                      MASS0(K) = M01
                   ELSE
                      EPOCH(J) = TM + (TSCLS(1)-TM)*(AJ(K)-TMS(K))/
     &                                (TBGB(K)-TMS(K))
                      EPOCH(J) = TEV0(J)*TSTAR - EPOCH(J)
                   ENDIF
                ELSE
                   EPOCH(J) = TEV0(J)*TSTAR - AJ(K)*TM/TMS(K)
                ENDIF
             ELSE
                CALL star(KW,MASS0(K),MASS(K),TM,TN,TSCLS,LUMS,GB,ZPARS)
             ENDIF
             AGE = TPHYS0 - EPOCH(J)
             IF((AGE-TN).GT.1.0D-02)THEN
                IF(KW.LE.6)THEN
                   AGE = MIN(AGE,0.9999*TSCLS(11))
                ELSE
                   AGE = MIN(AGE,0.9999*TSCLS(8))
                ENDIF
                EPOCH(J) = TPHYS0 - AGE
             ENDIF
   99        J = J2
   98     CONTINUE
*
          SEP = SEP*(MASS(1)+MASS(2)+DMS(1)+DMS(2))/(MASS(1)+MASS(2))
          SEMI = SEP/SU
          TB = YRS*SEMI*SQRT(SEMI*ZMBAR/(MASS(1)+MASS(2)))
          TK = DAYS*TB/YRS
          OORB = TWOPI/TB
          JORB = MASS(1)*MASS(2)/(MASS(1)+MASS(2))
     &           *SQRT(1.D0-ECC*ECC)*SEP*SEP*OORB
*
*       Define first/second Roche overflow stage and increase event counter.
          IF (KSTAR(I).LE.10.OR.
     &        (KSTAR(I).GT.10.AND.MOD(KSTAR(I),2).EQ.0)) THEN
              KSTAR(I) = KSTAR(I) + 1
              NRO = NRO + 1
              INEW = 1
              if(rank.eq.0)
     &        WRITE (6,8)  NAME(J1), NAME(J2), TPHYS, KW1, KW2,
     &                     KSTAR(I), MASS, SEP, TK, RAD, RL1
    8         FORMAT (' NEW ROCHE    NAM TP K* M A P R* RL ',
     &                               2I6,F9.2,3I4,2F8.4,F8.2,1P,4E10.2)
              IF(rank.eq.0.and.KSTAR(I).EQ.50)THEN
                 WRITE(6,9)NAME(J1),NAME(J2),KW1,KW2
    9            FORMAT(' WARNING: TOO MUCH ROCHE  NAM K* ',2I6,2I3)
              ENDIF
              IGR = 0
              IMB = 0
              IF (KZ(8).GT.3) THEN
                  CALL binev(IPAIR)
              END IF
*
              IF(rank.eq.0.and.FIRST)THEN
                 OPEN(UNIT=85,STATUS='UNKNOWN',FORM='FORMATTED',
     &                FILE='ROCHE')
                 FIRST = .FALSE.
                 WRITE (85,94)
   94            FORMAT (/,' NAM1  NAM2  K1 K2   TPHYS     AGE1     ',
     &                     ' AGE2     M01    M02    M1     M2     Z ',
     &                     '     e        P        JSPIN1      JSPIN2')
              ENDIF
              CH5 = ' NEW '
              if(rank.eq.0)then
              WRITE(85,95)NAME(J1),NAME(J2),KSTAR(J1),KSTAR(J2),
     &                    KSTAR(I),TPHYS,AJ(1),AJ(2),
     &                    MASS0(1),MASS0(2),MASS(1),MASS(2),
     &                    ZMET,ECC,TK,JSPIN(1),JSPIN(2),CH5
   95         FORMAT(2I7,3I3,3F10.3,4F7.3,F7.4,F6.3,1P,3E12.4,A5)
              CALL FLUSH(85)
              end if
          ENDIF
*
*       Evaluate 1/10 of nuclear time interval (yrs).
          M1 = MASS(1)
          CALL trdot(J1,DTM,M1)
          DTM = DTM*TSTAR
          DTMI(1) = DTM
          DTMI(2) = DTMI(1)
*
*       Include time-scale check for GR & MB braking for small separations.
*       (Note: use old orbital expressions)
          IF (SEP.LT.10.0) THEN
              DSEPG = -1.663D-09*MASS(1)*MASS(2)*(MASS(1) + MASS(2))/
     &                                          (SEP*SEP*SEP)*1.0D+06
              DTGR = -0.1D0*SEP/DSEPG
              DTM = MIN(DTGR,DTM)
              IF (KSTAR(J1).LE.1.AND.
     &            MASS(1).GT.0.35.AND.MASS(1).LT.1.25) THEN
                  DSEPM = -4.574D-07*(MASS(1) + MASS(2))**2*
     &                        RADX(1)**3/(MASS(2)*SEP**4)*1.0D+06
                  DTMB = -0.1D0*SEP/DSEPM
                  DTM = MIN(DTMB,DTM)
              END IF
          END IF
*
*       Adopt conservative time interval (expansion by factor 2 allowed).
          KM0 = DTM*1.0D+03/TB
          ZMU0 = BODY(J1)*BODY(J2)/BODY(I)
          HI = H(IPAIR)
          GO TO 50
      END IF
*
*       Specify period in yrs and N-body units.
   10 TK = SEMI*SQRT(SEMI/BODY(I))
      TB = YRS*TK
      TK = TWOPI*TK
*
*       Save energy, reduced mass & semi-major axis for corrections.
*
      HI = H(IPAIR)
      ZMU0 = BODY(J1)*BODY(J2)/BODY(I)
      SEMI0 = SEMI
      ECC2 = ECC*ECC
*
*       Obtain Eddington limit during one orbit.
      DME = 2.08D-03*EDDFAC*(1.D0/(1.D0 + ZPARS(11)))*RAD(2)*TB
      SUPEDD = .FALSE.
      NOVAE = .FALSE.
      DISK = .FALSE.
      RDMIN = 0.0425D0*SEP*(Q(2)*(1.D0+Q(2)))**(1.D0/4.D0)
      IF(RDMIN.GT.RAD(2)) DISK = .TRUE.
*
* Adopt Kelvin-Helmholz time from the modified classical expression.
*
      J = J1
      DO 15 K = 1,2
         TKH(K) = 1.0D+07*MASS(K)/(RAD(K)*ZLMSTY(J))
         IF(KSTAR(J).LE.1.OR.KSTAR(J).EQ.7.OR.KSTAR(J).GE.10)THEN
            TKH(K) = TKH(K)*MASS(K)
         ELSE
            TKH(K) = TKH(K)*(MASS(K) - MASSC(K))
         END IF
         J = J2
   15 CONTINUE
*
* Define dynamical timescale for the primary.
      TDYN = 5.05D-05*SQRT(RAD(1)**3/MASS(1))
* Identify special cases.
*
      IF(KSTAR(J1).EQ.2)THEN
         QC = 4.D0
      ELSEIF(KSTAR(J1).EQ.3.OR.KSTAR(J1).EQ.5.OR.KSTAR(J1).EQ.6)THEN
*        QC = (1.67D0-ZPARS(7)+2.D0*(MASSC(1)/MASS(1))**5)/2.13D0
*
* Alternatively use condition of Hjellming & Webbink, 1987, ApJ, 318, 794.
         QC = 0.362D0 + 1.D0/(3.D0*(1.D0 - MASSC(1)/MASS(1)))
      ELSEIF(KSTAR(J1).EQ.8.OR.KSTAR(J1).EQ.9)THEN
         QC = 0.784D0
      ELSE
         QC = 3.D0
      ENDIF
*
      IF(KSTAR(J1).EQ.0.AND.Q(1).GT.0.695)THEN
*
* This will be dynamical mass transfer of a similar nature to common-envelope
* evolution. The result is always a single star.
*
         TAUM = SQRT(TKH(1)*TDYN)
         DM1 = MASS(1)
         IF (KSTAR(J2).LE.1)THEN
*
* Restrict accretion to thermal timescale of secondary.
*
            DM2 = TAUM/TKH(2)*DM1
            MASS(2) = MASS(2) + DM2
*
* Rejuvenate the main sequence star.
*
            MASS0(2) = MASS(2)
            CALL star(KSTAR(J2),MASS0(2),MASS(2),TM,TN,
     &                                 TSCLS,LUMS,GB,ZPARS)
* If the star has no convective core then the effective age decreases,
* otherwise it will become younger still.
            IF(MASS(2).LT.0.35.OR.MASS(2).GT.1.25)THEN
               AJ(2) = TM/TMS(2)*AJ(2)*(MASS(2) - DM2)/MASS(2)
            ELSE
               AJ(2) = TM/TMS(2)*AJ(2)
            END IF
            EPOCH(J2) = TPHYS0 - AJ(2)
*
*       Check for blue straggler formation (TM < TPHYS).
            IF(rank.eq.0.and.TM.LT.TPHYS)THEN
               WRITE(6,16)NAME(J2),MASS(2),TM,TPHYS,AJ(2)
   16          FORMAT(' NEW BS (ROCHE):  NAM M3 TM TP AGE ',
     &                                   I6,F6.2,3F8.1)
            ENDIF
*
         ELSEIF(KSTAR(J2).LE.6)THEN
*
* Add all the material to the giant's envelope.
*
            DM2 = DM1
            MASS(2) = MASS(2) + DM2
            IF(KSTAR(J2).EQ.2)THEN
               MASS0(2) = MASS(2)
               CALL star(KSTAR(J2),MASS0(2),MASS(2),TM,TN,
     &                                    TSCLS,LUMS,GB,ZPARS)
               AJ(2) = TM + TSCLS(1)*(AJ(2)-TMS(2))/TBGB(2)
               EPOCH(J2) = TPHYS0 - AJ(2)
            ENDIF
         ELSEIF(KSTAR(J2).LE.12)THEN
*
* Form a new giant envelope.
*
            DM2 = DM1
*     --02/28/13 9:51-lwang-debug--------------------------------------*
***** Note:modife kst--------------------------------------------------**
            IF (MIN(KSTAR(J1),KSTAR(J2)).LT.0) THEN
               KST = 0
            ELSE
               KST = KTYPE(KSTAR(J1),KSTAR(J2))
            END IF
*     --02/28/13 9:52-lwang-end-debug----------------------------------*
            IF(KST.GT.100) KST = KST - 100
            IF(KST.EQ.4)THEN
               AJ(2) = AJ(2)/TMS(2)
               MASSC(2) = MASS(2)
            ENDIF
*
* Check for planets or low-mass WDs.
*
            IF((KSTAR(J2).EQ.10.AND.MASS(2).LT.0.05D0).OR.
     &         (KSTAR(J2).GE.11.AND.MASS(2).LT.0.5D0))THEN
               KST = KSTAR(J1)
               MASS(1) = MASS(2) + DM2
               MASS(2) = 0.D0
            ELSE
               MASS(2) = MASS(2) + DM2
               CALL gntage(MASSC(2),MASS(2),KST,ZPARS,MASS0(2),AJ(2))
               EPOCH(J2) = TPHYS0 - AJ(2)
            ENDIF
            KSTAR(J2) = KST
         ELSE
*
* The neutron star or black hole simply accretes at the Eddington rate.
*
            DM2 = MIN(DME*TAUM/TB,DM1)
            IF(DM2.LT.DM1) SUPEDD = .TRUE.
            MASS(2) = MASS(2) + DM2
         ENDIF
         COALS = .TRUE.
         IF(MASS(2).GT.0.D0) MASS(1) = 0.D0
         KW1 = KSTAR(J2)
         GOTO 60
      ELSEIF(((ABS(ABS(2*KSTAR(J1)-11)-3).EQ.2.OR.KSTAR(J1).EQ.9).
     &        AND.(Q(1).GT.QC.OR.RADX(1).LE.RADC(1))).OR.
     &        (KSTAR(J1).EQ.2.AND.Q(1).GT.QC).OR.
     &        (KSTAR(J1).EQ.4.AND.Q(1).GT.QC))THEN
*
* Common-envelope evolution for critical giants and subgiants with extreme q.
*
         M1CE = MASS(1)
         M2CE = MASS(2)
         if(rank.eq.0)
     &   WRITE (6,20)  MASS, KSTAR(J1), KSTAR(J2), RAD, ROL, SEP
   20    FORMAT (' NEW CE    M K R RL A ',2F7.3,2I3,4F9.3,F9.3)
         KW1 = KSTAR(J1)
         KW2 = KSTAR(J2)
         CALL comenv(MASS0(1),MASS(1),MASSC(1),AJ(1),JSPIN(1),KW1,
     &               MASS0(2),MASS(2),MASSC(2),AJ(2),JSPIN(2),KW2,
     &               ECC,SEP,COALS)
         if(rank.eq.0)
     &   WRITE (6,25)  MASS, KW1, KW2, SEP
   25    FORMAT (' END CE    M K A ',2F7.3,2I3,F9.3)
*
* Next step should be made without changing the time.
*
         DTM = 0.D0
         EPOCH(J1) = TPHYS0 - AJ(1)
         IF(COALS)THEN
            GOTO 60
         ELSE
            IF(KSTAR(J1).LT.13.AND.(KW1.EQ.13.OR.KW1.EQ.14))THEN
               JKICK = J1
               KWK = KW1
            ENDIF
            IF (KSTAR(J2).LT.13.AND.(KW2.EQ.13.OR.KW2.EQ.14))THEN
               JKICK = J2
               KWK = KW2
            ENDIF
            IF(KSTAR(J1).LT.10.AND.KW1.GE.10.AND.KW1.LE.12)THEN
               JKICK = J1
               KWK = KW1
            ENDIF
            IF(KSTAR(J2).LT.10.AND.KW2.GE.10.AND.KW2.LE.12)THEN
               JKICK = J2
               KWK = KW2
            ENDIF
            KSTAR(J1) = KW1
            KSTAR(J2) = KW2
            SEMI = SEP/SU
            DM1 = M1CE - MASS(1)
            DM2 = MASS(2) - M2CE
            DM22 = DM2
            DMS(1) = 0.D0
            DMS(2) = DM1 - DM2
*           jkick = 0
         ENDIF
         EPOCH(J2) = TPHYS0 - AJ(2)
*
      ELSEIF(KSTAR(J1).GE.10.AND.KSTAR(J1).LE.12.AND.Q(1).GE.0.628)THEN
*
* Dynamic transfer from a white dwarf.  Secondary will have KW > 9.
*
         TAUM = SQRT(TKH(1)*TDYN)
         DM1 = MASS(1)
         IF(EDDFAC.LT.10.D0)THEN
            DM2 = MIN(DME*TAUM/TB,DM1)
            IF(DM2.LT.DM1) SUPEDD = .TRUE.
         ELSE
            DM2 = DM1
         ENDIF
         MASS(2) = MASS(2) + DM2
         IF(KSTAR(J1).EQ.10.AND.KSTAR(J2).EQ.10)THEN
*
* Assume the energy released by ignition of the triple-alpha reaction
* is enough to destroy the star.
*
            KSTAR(J2) = 15
            MASS(2) = 0.D0
         ELSEIF(KSTAR(J1).EQ.10.OR.KSTAR(J2).EQ.10)THEN
*
* Should be helium overflowing onto a CO or ONe core in which case the
* helium swells up to form a giant envelope so a HeGB star is formed.
* Allowance for the rare case of CO or ONe flowing onto He is made.
*
            KSTAR(J2) = 9
            IF(KSTAR(J2).EQ.10) MASSC(J2) = DM2
            CALL gntage(MASSC(2),MASS(2),KSTAR(J2),ZPARS,MASS0(2),AJ(2))
            EPOCH(J2) = TPHYS0 - AJ(2)
         ELSEIF(KSTAR(J2).LE.12)THEN
            MASS0(2) = MASS(2)
*           IF(KSTAR(J1).EQ.12.AND.KSTAR(J2).EQ.11)THEN
*
* Mixture of ONe and CO
* or rapid accretion of CO on to CO 
* will result in an ONe product.
*
               KSTAR(J2) = 12
*           ENDIF
         ENDIF
         MASS(1) = 0.D0
         KW1 = KSTAR(J2)
*       See whether Chandrasekhar's mass is exceeded (SN).
         IF(KSTAR(J2).LE.11.AND.MASS(2).GT.MCH)THEN
            if(rank.eq.0)
     &      WRITE (6,27)  NAME(J1), NAME(J2), KSTAR(J2), MASS(2)
   27       FORMAT (' ROCHE SN    NAM K2* M2 ',2I6,I4,F6.2)
            MASS(2) = 0.D0
            NAS = NAS + 1
            KW1 = 15
         ENDIF
         COALS = .TRUE.
         GOTO 60
*
      ELSEIF(KSTAR(J1).EQ.13)THEN
*
* Gamma ray burster?
*
         DM1 = MASS(1)
         MASS(1) = 0.D0
         KW1 = 14
         DM2 = DM1
         MASS(2) = MASS(2) + DM2
         COALS = .TRUE.
         NGB = NGB + 1
         GOTO 60
      ELSEIF(KSTAR(J1).EQ.14)THEN
*
* Both stars are black holes.  Let them merge quietly.
*
         DM1 = MASS(1)
         MASS(1) = 0.D0
         KW1 = KSTAR(J2)
         DM2 = DM1
         MASS(2) = MASS(2) + DM2
         COALS = .TRUE.
         GOTO 60
      ELSE
*
*       Form mass transfer for one Kepler orbit.
*
         DM1 = 3.0D-06*TB*(LOG(RAD(1)/ROL(1))**3)*
     &         MIN(MASS(1),5.D0)**2
         IF(KSTAR(J1).EQ.2)THEN
            MEW = (MASS(1) - MASSC(1))/MASS(1)
            DM1 = MAX(MEW,0.01D0)*DM1
         ELSEIF(KSTAR(J1).GE.10)THEN
            DM1 = DM1*1.0D+03/MAX(RAD(1),1.0D-04)
         ENDIF
         KST = KSTAR(J2)
*
* Limit mass transfer to the thermal rate for remaining giant-like stars and
* to the dynamical rate for all others.
*
         IF(KSTAR(J1).GE.2.AND.KSTAR(J1).LE.9.AND.KSTAR(J1).NE.7)THEN
            DM1 = MIN(DM1,MASS(1)*TB/TKH(1))
         ELSEIF(RAD(1).GT.10.D0*ROL(1).OR.(KSTAR(J1).LE.1.AND.
     &          KSTAR(J2).LE.1.AND.Q(1).GT.QC))THEN
*
* Although zeta_ad is formally infinite we do not expect a star to overfill
* its Roche lobe by more than a factor of 10 before the stars merge.
*
            if(rank.eq.0)
     &      WRITE (6,28)  SEP, DM1, RAD(1), ROL(1)
   28       FORMAT (' OVERFILL    A DM1 RAD ROL ',F7.1,1P,3E10.2)
*
*       Form new star by combining the components inelastically.
            IQCOLL = 3
            CH5 = ' OVER'
            GOTO 160
         ELSE
            DM1 = MIN(DM1,MASS(1)*TB/TDYN)
         ENDIF
*
* Calculate wind mass loss from the stars during one orbit.
*
         VORB2 = ACC1*(MASS(1)+MASS(2))/SEP
         IVSQM = 1.D0/SQRT(1.D0-ECC*ECC)
         J = J1
         DO 30 K = 1,2
            RLPERI = ROL(K)*(1.D0-ECC)
            DMR(K) = MLWIND(KSTAR(J),ZLMSTY(J),RADX(K),
     &                      MASS(K),MASSC(K),RLPERI,ZMET)
            VWIND2 = 2.D0*BETA*ACC1*MASS(K)/RADX(K)
            OMV2 = (1.D0 + VORB2/VWIND2)**(3.D0/2.D0)
            DMA(3-K) = IVSQM*ACC2*DMR(K)*((ACC1*MASS(3-K)/VWIND2)**2)/
     &                 (2.D0*SEP*SEP*OMV2)
            DMA(3-K) = MIN(DMA(3-K),DMR(K))
            J = J2
 30      CONTINUE
*
         DO 32 K = 1,2
            DMS(K) = (DMR(K)-DMA(K))*TB
 32      CONTINUE
*
*       Increase time-scale to relative mass loss of 0.5% but < c.m. step.
         KM = 5.0D-03/MAX(ABS(DM1+DMS(1))/MASS(1),DMS(2)/MASS(2))
         KM = MIN(KM,2.D0*KM0)
         KM0 = KM
         DTM = KM*TB/1.0D+06
         DT1 = DTM/TSTAR
*
* Take the stellar evolution timestep into account but don't let it
* be overly restrictive for long lived phases.
*
         IF(ITER.LE.100) DTM = MIN(DTM,DTMI(1),DTMI(2))
         DTM = MIN(DTM,(TIME-TEV0(I))*TSTAR)
         DTM = MAX(DTM,1.0D-10)
         KM = DTM*1.0D+06/TB
*
*       Decide between accreted mass by secondary and/or system mass loss.
*
         TAUM = MASS(2)/DM1*TB
         IF(KSTAR(J2).LE.2.OR.KSTAR(J2).EQ.4)THEN
*
* Limit mass loss according to the thermal timescale of the secondary.
*
            DM2 = MIN(1.D0,10.D0*TAUM/TKH(2))*DM1
         ELSEIF(KSTAR(J2).GE.7.AND.KSTAR(J2).LE.9)THEN
*
* Naked helium star secondary swells up to a core helium burning star
* or SAGB star unless the primary is also a helium star.
*
            IF(KSTAR(J1).GE.7)THEN
               DM2 = MIN(1.D0,10.D0*TAUM/TKH(2))*DM1
            ELSE
               DM2 = DM1
               DMCHK = DM2 - 1.05D0*DMS(2)
               IF(DMCHK.GT.0.D0.AND.DM2/MASS(2).GT.1.0D-04)THEN
                  KST = MIN(6,2*KSTAR(J2)-10)
                  IF(KST.EQ.4)THEN
                     AJ(2) = AJ(2)/TMS(2)
                     MCX = MASS(2)
                  ELSE
                     MCX = MASSC(2)
                  ENDIF
                  M2 = MASS(2) + KM*(DM2 - DMS(2))
                  CALL gntage(MCX,M2,KST,ZPARS,MASS0(2),AJ(2))
                  EPOCH(J2) = TPHYS0 + DTM - AJ(2)
               ENDIF
            ENDIF
         ELSEIF(KSTAR(J1).LE.6.AND.
     &           (KSTAR(J2).GE.10.AND.KSTAR(J2).LE.12))THEN
*
* White dwarf secondary.
*
            IF(DM1/TB.LT.2.71D-07)THEN
               IF(DM1/TB.LT.1.03D-07)THEN
*
* Accrete until a nova explosion blows away most of the accreted material.
*
                  NOVAE = .TRUE.
                  DM2 = MIN(DM1,DME)
                  IF(DM2.LT.DM1) SUPEDD = .TRUE.
                  DM22 = EPSNOV*DM2
               ELSE
*
* Steady burning at the surface.
*
                  DM2 = DM1
               ENDIF
*
            ELSE
*
* Make a new giant envelope.
*
               DM2 = DM1
*
* Check for planets or low-mass WDs.
*
               IF((KSTAR(J2).EQ.10.AND.MASS(2).LT.0.05D0).OR.
     &            (KSTAR(J2).GE.11.AND.MASS(2).LT.0.5D0))THEN
                  KST = KSTAR(J2)
               ELSE
                  KST = MIN(6,3*KSTAR(J2)-27)
                  M2 = MASS(2) + KM*(DM2 - DMS(2))
                  CALL gntage(MASSC(2),M2,KST,ZPARS,MASS0(2),AJ(2))
                  EPOCH(J2) = TPHYS0 + DTM - AJ(2)
               ENDIF
            ENDIF
         ELSEIF(KSTAR(J2).GE.10)THEN
*
* Impose the Eddington limit.
*
            DM2 = MIN(DM1,DME)
            IF(DM2.LT.DM1) SUPEDD = .TRUE.
         ELSE
*
* We have a giant whose envelope can absorb any transferred material.
*
            DM2 = DM1
         ENDIF
         IF(.NOT.NOVAE) DM22 = DM2
*
         IF(KST.GE.10.AND.KST.LE.12)THEN
            MT2 = MASS(2) + KM*(DM22 - DMS(2))
            IF(KSTAR(J1).LE.10.AND.KST.EQ.10.AND.MT2.GE.0.7)THEN
*
* HeWD can only accrete helium-rich material up to a mass of 0.7 when
* it is destroyed in a possible Type 1a SN.
*
               MASS(1) = MASS(1) - KM*(DM1 + DMS(1))
               MASS(2) = 0.D0
               KW1 = KSTAR(J1)
               COALS = .TRUE.
               GOTO 60
            ELSEIF(KSTAR(J1).LE.10.AND.KST.GE.11)THEN
*
* CO and ONeWDs accrete helium-rich material until the accumulated
* material exceeds a mass of 0.15 when it ignites. For a COWD with
* mass less than 0.95 the system will be destroyed in a possible
* Type 1a SN. COWDs with mass greater than 0.95 and ONeWDs will survive
* with all the material converted to ONe.
*
** Now changed to an ELD for all COWDs when 0.15 accreted (JH 11/01/00).
*
               IF((MT2-MASS0(2)).GE.0.15)THEN
                  IF(KST.EQ.11)THEN
                     MASS(1) = MASS(1) - KM*(DM1 + DMS(1))
                     MASS(2) = 0.D0
                     KW1 = KSTAR(J1)
                     COALS = .TRUE.
                     GOTO 60
                  ENDIF
                  MASS0(2) = MT2
               ENDIF
            ELSE
               MASS0(2) = MT2
            ENDIF
*
* If the Chandrasekhar limit is exceeded for a white dwarf then destroy
* the white dwarf in a supernova. If the WD is ONe then a neutron star
* will survive the supernova and we let HRDIAG take care of this when
* the stars are next updated.
* However, first check if CO WD is accreting at a fast enough rate to 
* convert to ONe (C. Tout 21/11/08). 
*
            IF(KST.EQ.11.AND.DM22.GT.0.4D0*DME/EDDFAC) KST = 12
            IF((KST.EQ.10.OR.KST.EQ.11).AND.MT2.GE.MCH)THEN
               DM1 = MCH - MASS(2) + KM*DMS(2)
               MASS(1) = MASS(1) - DM1 - KM*DMS(1)
               MASS(2) = 0.D0
               KW1 = KSTAR(J1)
               COALS = .TRUE.
               GOTO 60
            ENDIF
         ENDIF
*
*       Modify time-step & mass loss terms by speed-up factor.
         DM1 = KM*DM1
         DM2 = KM*DM2
         DM22 = KM*DM22
         DME = KM*DME
*       Save last Roche step (yrs).
         DTM0 = DTM*1.0D+06
*
* Transform to angular momentum for calculation of orbital changes.
*
         OORB = TWOPI/TB
         JORB = MASS(1)*MASS(2)/(MASS(1)+MASS(2))
     &          *SQRT(1.D0-ECC*ECC)*SEP*SEP*OORB
*
* Calculate orbital angular momentum change due to system mass loss.
*
         DJORB = ((DMR(1)+Q(1)*DMA(1))*MASS(2)*MASS(2) +
     &            (DMR(2)+Q(2)*DMA(2))*MASS(1)*MASS(1))*
     &           SEP*SEP*OORB/(MASS(1)+MASS(2))**2
         DJORB = DJORB*DTM0
*
* For super-Eddington mass transfer rates, for gamm1 = -2.0,
* and for novae systems, assume that material is lost from
* the system as if a wind from the secondary.
* If gamm1 = -1.0 then assume the lost material carries with it
* the specific angular momentum of the primary and for all
* gamm1 > 0.0 assume that it takes away a fraction gamma of
* the orbital angular momentum.
*
*
         IF(SUPEDD.OR.NOVAE.OR.GAMM1.LT.-1.5D0)THEN
            DJORB = DJORB + (DM1 - DM22)*MASS(1)*MASS(1)*
     &              SEP*SEP*OORB/(MASS(1)+MASS(2))**2
         ELSEIF(GAMM1.GE.0.D0)THEN
            DJORB = DJORB + GAMM1*(DM1 - DM2)*SEP*SEP*OORB
         ELSE
            DJORB = DJORB + (DM1 - DM2)*MASS(2)*MASS(2)*
     &              SEP*SEP*OORB/(MASS(1)+MASS(2))**2
         ENDIF
*
* For very close systems include gravitational radiation.
*
         IF(SEP.LE.10.0)THEN
            CALL GRRAD(MASS(1),MASS(2),SEP,ECC,JORB,DJGR,DELET1)
            DJGR = DJGR*DTM0
            DJORB = DJORB + DJGR
            IGR = IGR + 1
*           IF(IGR.LT.10.OR.ABS(DJGR)/JORB.GT.0.001)THEN
            IF(SEP.LT.1.05*(RAD(1) + RAD(2)))THEN
            IF (rank.eq.0.and.MOD(IGR,10).EQ.0) 
     &      WRITE (6,45)  MASS, SEP, DJGR, DTM
   45       FORMAT (' GR BRAKE    M1 M2 SEP DJ DTM ',2F7.3,1P,3E10.2)
            ENDIF
         ENDIF
*
         DMS(1) = KM*DMS(1)
         IF(KSTAR(J1).LT.10) DMS(1) = MIN(DMS(1),MASS(1) - MASSC(1))
         DMS(2) = KM*DMS(2)
         IF(KSTAR(J2).LT.10) DMS(2) = MIN(DMS(2),MASS(2) - MASSC(2))
*
* Calculate spin changes owing to stellar winds and magnetic braking.
*
         DJSPIN(1) = (2.D0/3.D0)*(DMR(1)*RADX(1)*RADX(1)*OSPIN(1) -
     &                XI*DMA(1)*RADX(2)*RADX(2)*OSPIN(2))*DTM0
         DJSPIN(2) = (2.D0/3.D0)*(DMR(2)*RADX(2)*RADX(2)*OSPIN(2) -
     &                XI*DMA(2)*RADX(1)*RADX(1)*OSPIN(1))*DTM0
         IF(KSTAR(J1).LT.10.AND.MASS(1).GE.0.35)THEN
            CALL MAGBRK(KSTAR(J1),MASS(1),MENV(1),RAD(1),OSPIN(1),DJMB)
            DJSPIN(1) = DJSPIN(1) + DJMB*DTM0
            IMB = IMB + 1
            IF(JSPIN(1).GT.0.D0)THEN
               IF(IMB.LT.2.OR.ABS(DJMB)/JSPIN(1).GT.0.03)THEN
                  if(rank.eq.0)
     &            WRITE (6,40)  MASS, SEP, DJMB, DTM
   40             FORMAT (' MB BRAKE    M1 M2 SEP DJSPN DTM ',
     &                                  2F7.3,1P,3E10.2)
               ENDIF
            ENDIF
         ENDIF
         IF(KSTAR(J2).LT.10.AND.MASS(2).GE.0.35)THEN
            CALL MAGBRK(KSTAR(J2),MASS(2),MENV(2),RAD(2),OSPIN(2),DJMB)
            DJSPIN(2) = DJSPIN(2) + DJMB*DTM0
            IMB = IMB + 1
            IF(JSPIN(2).GT.0.D0)THEN
               IF(IMB.LT.10.OR.ABS(DJMB)/JSPIN(2).GT.0.001)THEN
                  if(rank.eq.0)
     &            WRITE (6,40)  MASS, SEP, DJMB, DTM
               ENDIF
            ENDIF
         ENDIF
*
* Adjust the spin angular momentum of each star owing to mass transfer
* and conserve total angular momentum.
*
         DJT = DM1*RADX(1)*RADX(1)*OSPIN(1)
         DJSPIN(1) = DJSPIN(1) + DJT
         DJORB = DJORB - DJT
         IF(DISK)THEN
*
* Alter spin of the degenerate secondary by assuming that material
* falls onto the star from the inner edge of a Keplerian accretion
* disk and that the system is in a steady state.
*
            DJT = DM2*TWOPI*AURSUN*SQRT(AURSUN*MASS(2)*RADX(2))
            DJSPIN(2) = DJSPIN(2) - DJT
            DJORB = DJORB + DJT
*
         ELSE
*
* No accretion disk.
* Calculate the angular momentum of the transferred material by
* using the radius of the disk (see Ulrich & Burger) that would
* have formed if allowed.
*
            RDISK = 1.7D0*RDMIN
            DJT = DM2*TWOPI*AURSUN*SQRT(AURSUN*MASS(2)*RDISK)
            DJSPIN(2) = DJSPIN(2) - DJT
            DJORB = DJORB + DJT
*
         ENDIF
*
* Adjust the secondary spin if a nova eruption has occurred.
*
         IF(NOVAE)THEN
            DJT = (DM2 - DM22)*RADX(2)*RADX(2)*OSPIN(2)
            DJSPIN(2) = DJSPIN(2) + DJT
         ENDIF
*
* Calculate any changes owing to tidal evolution (includes eccentricity).
*
         J = J1
         DELET = 0.D0
         DO 500 K = 1,2
            IF((KSTAR(J).LE.9.AND.RAD(K).GE.0.01D0*ROL(K)).OR.
     &         (KSTAR(J).GE.10.AND.J.EQ.J1))THEN
*
               CALL BSETID(kstar(j),mass(k),massc(k),menv(k),radx(k),
     &                     radc(k),renv(k),zlmsty(j),ospin(k),k2str(k),
     &                     q(3-k),sep,ecc,oorb,delet1,dspin,eqspin,djt)
*
               DELET = DELET + DELET1*DTM0
               IF(DTM0.GT.0.D0.AND.ABS(DSPIN).GE.TINY)THEN
                  DSPIN0 = DSPIN
                  IF(DSPIN.GE.0.D0)THEN
                     DSPIN = MIN(DTM0*DSPIN,EQSPIN-OSPIN(K))/DTM0
                  ELSE
                     DSPIN = MAX(DTM0*DSPIN,EQSPIN-OSPIN(K))/DTM0
                  ENDIF
                  DJT = DJT*DSPIN/DSPIN0
               ELSE
                  DJT = 0.D0
               ENDIF
               DJORB = DJORB + DJT*DTM0
               DJSPIN(K) = DJSPIN(K) - DJT*DTM0
*
            ENDIF
*
            J = J2
 500     CONTINUE
*
* Update the initial and current masses.
*
         KSTAR(J2) = KST
         MASS(1) = MASS(1) - DM1 - DMS(1)
         IF(KSTAR(J1).LE.1.OR.KSTAR(J1).EQ.7) MASS0(1) = MASS(1)
         MASS(2) = MASS(2) + DM22 - DMS(2)
         IF(KSTAR(J2).LE.1.OR.KSTAR(J2).EQ.7) MASS0(2) = MASS(2)
*
* Update spins and ensure that the star does not spin up beyond break-up.
* For a HG star check if the initial mass can be reduced.
*
         J = J1
         DO 502 K = 1,2
            JSPIN(K) = MAX(JSPIN(K) - DJSPIN(K),1.0D-10)
            OSPBRU = TWOPI*SQRT(MASS(K)*AURSUN**3/RADX(K)**3)
            JSPBRU = (K2STR(K)*(MASS(K)-MASSC(K))*RADX(K)*RADX(K) +
     &                K3*MASSC(K)*RADC(K)*RADC(K))*OSPBRU
            IF(JSPIN(K).GT.JSPBRU)THEN
               DJORB = DJORB - (JSPIN(K) - JSPBRU)
               JSPIN(K) = JSPBRU
            ENDIF
            IF(KSTAR(J).EQ.2.AND.MASS0(K).LE.ZPARS(3))THEN
               M01 = MASS0(K)
               MASS0(K) = MASS(K)
               CALL star(KSTAR(J),MASS0(K),MASS(K),TMSNEW,TN,TSCLS,
     &                   LUMS,GB,ZPARS)
               IF(GB(9).LT.MASSC(K)) MASS0(K) = M01
            ENDIF
            J = J2
 502     CONTINUE
*
* Combine the mass lost during transfer with the secondary's wind.
*
         DMS(2) = DMS(2) + DM1 - DM22
*
* Adjust the orbital angular momentum and reset the orbital separation.
*
         ECC = MAX(ECC - DELET,0.001D0)
         JORB = MAX(JORB - DJORB,0.D0)
         SEP = (MASS(1) + MASS(2))*JORB*JORB/
     &         ((MASS(1)*MASS(2)*TWOPI)**2*AURSUN**3*(1.D0-ECC*ECC))
         SEP = MAX(SEP,RAD(2))
         SEMI = SEP/SU
*
*       Produce diagnostics for different types of CV and X-ray binaries.
*        if(rank.eq.0)
*    &   WRITE (19,38)  NAME(J1), NAME(J2), KSTAR(J1), KSTAR(J2),
*    &                   TPHYS, MASS, SEP, RAD(1)/ROL(1),
*    &                   DM1/DTM/1.0E06, DM2/DTM/1.0E06, DTM
*  38    FORMAT (' ROCHE    NAM K* T M A R/RL MD12 DT ',
*    &                       2I6,2I4,F9.2,2F6.1,F7.2,F6.2,1P,3E9.1)
*        CALL flush(19)
*
*       Update KS variables to SEMI1 while maintaining low eccentricity.
*        SEMI1 = SEMI
*        CALL decay(IPAIR,SEMI1)
*
      ENDIF
*
* Rejuvenate secondary and age primary if they are on the main sequence.
*
      IF(KSTAR(J1).LE.2.OR.KSTAR(J1).EQ.7)THEN
         CALL star(KSTAR(J1),MASS0(1),MASS(1),TM,TN,
     &             TSCLS,LUMS,GB,ZPARS)
         IF(KSTAR(J1).EQ.2)THEN
            AJ(1) = TM + (TSCLS(1)-TM)*(AJ(1)-TMS(1))/(TBGB(1)-TMS(1))
         ELSE
            AJ(1) = TM/TMS(1)*AJ(1)
         ENDIF
         EPOCH(J1) = TPHYS0 - AJ(1)
      ENDIF
*
      IF(KSTAR(J2).LE.2.OR.KSTAR(J2).EQ.7)THEN
         CALL star(KSTAR(J2),MASS0(2),MASS(2),TM,TN,
     &             TSCLS,LUMS,GB,ZPARS)
         IF(KSTAR(J2).EQ.2)THEN
            AJ(2) = TM + (TSCLS(1)-TM)*(AJ(2)-TMS(2))/(TBGB(2)-TMS(2))
         ELSEIF((MASS(2).LT.0.35.OR.MASS(2).GT.1.25).AND.
     &      KSTAR(J2).NE.7)THEN
            AJ(2) = TM/TMS(2)*AJ(2)*(MASS(2) - DM22)/MASS(2)
         ELSE
            AJ(2) = TM/TMS(2)*AJ(2)
         END IF
         EPOCH(J2) = TPHYS0 - AJ(2)
*
*       Check for blue straggler formation (TM < TPHYS).
         IF(KSTAR(J2).LE.1)THEN
            IF(TMS(2).GT.TPHYS.AND.TM.LT.TPHYS)THEN
               NBR = NBR + 1
               if(rank.eq.0)
     &         WRITE(6,48)NAME(J2),MASS(2),DM2,SEMI*SU,TM,TPHYS,
     &                    KSTAR(J1)
   48          FORMAT (' NEW BR    NAM M2 DM2 A TM TP KW1 ',
     &                             I6,2F7.3,3F8.1,I4)
            ENDIF
         ENDIF
*
      ENDIF
      TPHYS0 = TPHYS0 + DTM
      TPHYS = TPHYS + DTM
*
*       Accumulate the net system mass loss.
   50 ZMSY = ZMSY + DMS(1) + DMS(2)
*
*       Update masses.
      BODY(J1) = MASS(1)/ZMBAR
      BODY(J2) = MASS(2)/ZMBAR
      BODY0(J1) = MASS0(1)/ZMBAR
      BODY0(J2) = MASS0(2)/ZMBAR
*
*       Update mass loss & binary mass and re-define binding energy.
      DMT = DMS(1) + DMS(2)
      DMT = DMT/ZMBAR
      BODY(I) = BODY(I) - DMT
      H(IPAIR) = -0.5D0*BODY(I)/SEMI
*
*       Include mass loss from system.
      DT2 = 1.D0
      IF (DMT.GT.0.0) THEN
          DT2 = MIN(DT2,0.02D0*BODY(J1)*DTM/(DMT*TSTAR))
*       Modify KS variables at constant eccentricity (expand or contract).
*         CALL expand(IPAIR,SEMI1)
          CALL expand(IPAIR,SEMI0)
*       Resolve KS components.
          CALL resolv(IPAIR,1)
*       Copy neighbour list for force corrections (no need!).
*         IPHASE = 0
*         NNB = LIST(1,I)
*         DO 55 L = 1,NNB+1
*             ILIST(L) = LIST(L,I)
*  55     CONTINUE
*       Correct neighbour forces and obtain the corresponding energy loss.
          IF (DMT*SMU.LT.0.05) THEN
              CALL ficorr(I,DMT)
          ELSE
              CALL fcorr(I,DMT,0)
          END IF
          ZMASS = ZMASS - DMT
      ELSE IF (ABS(SEMI - SEMI0).GT.1.0D-10*SEMI) THEN
*       Modify KS variables at constant eccentricity to treat change in SEMI.
          CALL expand(IPAIR,SEMI0)
      END IF
*
*       Subtract energy correction terms due to change in relative motion.
      ZMU = BODY(J1)*BODY(J2)/BODY(I)
      EMDOT = EMDOT + ZMU0*HI - ZMU*H(IPAIR)
      EGRAV = EGRAV + ZMU0*HI - ZMU*H(IPAIR)
      DE = ZMU0*HI - ZMU*H(IPAIR)
*
*       Increase event counter and update mass-loss.
      NROCHE = NROCHE + 1
      ZMRO = ZMRO + DM1
*
      TB = YRS*SEMI*SQRT(SEMI/BODY(I))
      TK = DAYS*TB/YRS
*
*       Check termination by coalescence.
   60 IF(COALS)THEN
          TPHYS = (TIME+TOFF)*TSTAR
          BODY0(J1) = MASS0(1)/ZMBAR
          BODY0(J2) = MASS0(2)/ZMBAR
          SPIN(J1) = JSPIN(1)/SPNFAC
          SPIN(J2) = JSPIN(2)/SPNFAC
          KW2 = 15
*       Ensure zero BODY0 for secondary to avoid confusion in COAL
*       unless both stars heve zero MASS.
          II = J2
          K = 2
          IF(MASS(1).EQ.0.0.AND.MASS(2).NE.0.0)THEN
             II = J1
             K = 1
          ENDIF
          BODY0(II) = 0.D0
          TEV(J1) = TPHYS0/TSTAR
          CH5 = ' COAL'
          TK = 0.0
          if(rank.eq.0)
     &    WRITE(85,95)NAME(J1),NAME(J2),KW1,KW2,KSTAR(I),
     &                TPHYS,AJ(K),TK,MASS0(K),MASS0(3-K),
     &                MASS(K),TK,ZMET,ECC,TK,JSPIN(K),TK,CH5
*
          CALL coal(IPAIR,KW1,KW2,MASS)
          GO TO 200
      ENDIF
*
* Obtain the stellar parameters for the next step (first J1, then J2).
*
      J = J1
      IQ = 0
      DO 70 K = 1,2
         AGE = TPHYS0 - EPOCH(J)
         M01 = MASS0(K)
         M1 = MASS(K)
         MC = MASSC(K)
         KW = KSTAR(J)
         CALL star(KW,M01,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
         CALL hrdiag(M01,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               R1,LUM,KW,MC,RCC,M1CE,RCE,K2)
         CALL trdot2(KW,AGE,TM,TN,TSCLS,DTMI(K),DTR)
         IF(K.EQ.1) DTMAX = DTR/TSTAR
*
*       Update evolution times.
*
         TEV(J) = TPHYS0/TSTAR
         TEV0(J) = TEV(J)
*        if(rank.eq.0)
*    &   write(6,*)' roche star ',j,name(j),kstar(j),kw,log10(lum)
*
*       Quit Roche on any transition that involves a supernova as the mass
*       will change suddenly.
*
         IF((KSTAR(J).LT.13.AND.KW.GE.13).OR.
     &      (KSTAR(J).LT.10.AND.KW.GE.10.AND.KW.LE.12))THEN
            IQ = J
            BODY0(J) = MASS0(K)/ZMBAR
            BODY(J) = MASS(K)/ZMBAR
            SPIN(J) = JSPIN(K)/SPNFAC
            KWK = KW
            JK = J
            GO TO 68
         ENDIF
*
* Save relevant solar quantities.
*
         AJ(K) = AGE
         IF(KW.NE.KSTAR(J))THEN
            EPOCH(J) = TPHYS0 - AGE
            KSTAR(J) = KW
            MASS0(K) = M01
            BODY0(J) = MASS0(K)/ZMBAR
            MASS(K) = M1
            BODY(J) = MASS(K)/ZMBAR
         ENDIF
         RADIUS(J) = R1/SU
         RAD(K) = R1
         RADC(K) = RCC
         MASSC(K) = MC
         RENV(K) = RCE
         MENV(K) = M1CE
         K2STR(K) = K2
         ZLMSTY(J) = LUM
         TMS(K) = TM
         TBGB(K) = TSCLS(1)
         Q(K) = MASS(K)/MASS(3-K)
         ROL(K) = RL(Q(K))*SEP
         RADX(K) = R1
         IF(J.EQ.J1)THEN
            RADX(K) = MIN(RADX(K),ROL(K))
            RADX(K) = MAX(RCC,RADX(K))
         ENDIF
         OSPIN(K) = JSPIN(K)/(K2*RADX(K)*RADX(K)*(M1-MC) +
     &                        K3*RCC*RCC*MC)
         SPIN(J) = JSPIN(K)/SPNFAC
   68    J = J2
   70 CONTINUE
*
*       Check Roche termination due to sudden mass loss.
      TEV0(I) = TEV0(J1)
      IF(JKICK.GT.0)THEN
*       Implement velocity kick without mass loss for NS & BH.
         KSTAR(I) = KSTAR(I) + 1
         CH5 = ' KICK'
         if(rank.eq.0)
     &   WRITE(85,95)NAME(J1),NAME(J2),KSTAR(J1),KSTAR(J2),KSTAR(I),
     &               TPHYS,AJ(1),AJ(2),MASS0(1),MASS0(2),
     &               MASS(1),MASS(2),ZMET,ECC,TK,JSPIN(1),JSPIN(2),CH5
         CALL KICK2(JKICK)
         GOTO 160
      ELSE IF(IQ.GT.0)THEN
         TEV(I) = 1.0D+10
         TEV(J1) = TPHYS0/TSTAR
         TEV(J2) = TEV(J1)
         KSTAR(I) = KSTAR(I) + 1
         KSTAR(JK) = KWK
         CH5 = ' KICK'
         if(rank.eq.0)
     &   WRITE(85,95)NAME(J1),NAME(J2),KSTAR(J1),KSTAR(J2),KSTAR(I),
     &               TPHYS,AJ(1),AJ(2),MASS0(1),MASS0(2),
     &               MASS(1),MASS(2),ZMET,ECC,TK,JSPIN(1),JSPIN(2),CH5
         WRITE (6,71)  NAME(IQ), NAME(J1), NAME(J2), KSTAR(J1),
     &                 KSTAR(J2), KW1, MASS(1), MASS(2), RAD(1), ROL(1)
   71    FORMAT (' ROCHE TERM    NAM K* KW M R1 RL1 ',
     &                           3I6,3I4,2F7.3,2F8.2)
         GOTO 140
      ENDIF
*
***
      IF(ITER.EQ.0.AND.INEW.EQ.1)THEN
*
*       Circularize the orbit if necessary.
*       Force co-rotation of primary with orbit if necessary.
*
         TC = 1.0D+04
         TSYN = 1.0D+20
         ITERB = 0
         ECC0 = ECC
         FAC0 = OSPIN(1)/OORB
         SEMI0 = SEMI
*
  80     CONTINUE
*
         FAC = OSPIN(1)/OORB
         IF(FAC.GT.1.D0) FAC = 1.D0/FAC
*
         IF((ECC-0.001D0).GT.TINY.OR.FAC.LT.0.98D0)THEN
* Calculate changes owing to tidal evolution.
            DELET = 0.D0
            J = J1
            DO K = 1,2
               DSPINK(K) = 0.D0
               DJTK(K) = 0.D0
               IF((KSTAR(J).LE.9.AND.RAD(K).GE.0.01D0*ROL(K)).OR.
     &            (KSTAR(J).GE.10.AND.J.EQ.J1))THEN
*
                  CALL BSETID(kstar(j),mass(k),massc(k),menv(k),radx(k),
     &                      radc(k),renv(k),zlmsty(j),ospin(k),k2str(k),
     &                      q(3-k),sep,ecc,oorb,delet1,dspin,eqspin,djt)
*
                  DELET = DELET + DELET1
                  DSPINK(K) = DSPIN
                  DJTK(K) = DJT
               ENDIF
               J = J2
            ENDDO
* Include gravitational radiation.
            CALL GRRAD(MASS(1),MASS(2),SEP,ECC,JORB,DJGR,DELET1)
            DELET = DELET + DELET1
* Include magnetic braking.
            CALL MAGBRK(KSTAR(J1),MASS(1),MENV(1),RAD(1),OSPIN(1),DJMB)
            DJSPIN(1) = DJTK(1) + DJMB
*
            TTID = 1.0D+20
            DTM0 = 1.0D+20
            DTMMIN = 0.D0
            IF(ECC.GT.0.001D0)THEN
               IF(ABS(DELET).GT.1.0D-14)THEN
* Allow eccentricity to change by 10% at most with a minimum of 0.01.
                  DTM0 = MAX(0.1D0*ECC,0.0005D0)/ABS(DELET)
                  DTMMIN = MIN(0.01D0,ECC-0.001D0)/ABS(DELET)
                  TC = ECC/ABS(DELET)
                  TC = TC/1.0D+06
                  TTID = TC
               ELSE
* Perform instant circularization conserving orbital angular momentum. 
                  ECC = 0.001D0
                  SEP = JORB*(MASS(1)+MASS(2))/(MASS(1)*MASS(2)*OORB)
                  SEP = SQRT(SEP)
                  TB = (SEP/AURSUN)*SQRT(SEP/(AURSUN*(MASS(1)+MASS(2))))
                  TK = DAYS*TB/YRS
                  OORB = TWOPI/TB
               ENDIF
            ELSE
               DELET = 0.D0
            ENDIF
            IF(ABS(DJSPIN(1)).GT.1.0D-14)THEN
               IF(ABS(DSPINK(1)).GT.1.0D-20)THEN
                  TSYN = ABS(OORB-OSPIN(1))/ABS(DSPINK(1))
                  DTM0 = MIN(DTM0,0.3D0*TSYN)
                  TSYN = TSYN/1.0D+06
                  TTID = MIN(TTID,TSYN)
               ELSE
* Allow jspin to change by 50% at most.
                  DTM0 = MIN(DTM0,0.5D0*JSPIN(1)/ABS(DJSPIN(1)))
               ENDIF
            ENDIF
*
            IF(ITERB.EQ.0.AND.NWARN.LT.50)THEN
               if(rank.eq.0)
     &         WRITE(6,710)NAME(J1),KSTAR(J1),KSTAR(J2),ECC0,FAC0,
     &                     TC,TSYN
  710          FORMAT(' CIRC & SYNCH NM K* ECC0 SPIN1/OORB TC TSYN ',
     &                                I7,2I4,F7.4,F9.3,1P,2E10.2)
               NWARN = NWARN + 1
            ENDIF
*
            IF(TTID.LT.MIN(TPHYS,10.D0))THEN
*     --02/28/13 9:54-lwang-debug--------------------------------------*
***** Note:debug-------------------------------------------------------**
               ITRO = 0
               DTM0 = MAX(DTM0,DTMMIN)
               DELET = DELET*DTM0
 81            DJORB = DJGR*DTM0
*     --02/28/13 9:54-lwang-end-debug----------------------------------*
* Include magnetic braking for both stars.
               DJSPIN(1) = DJMB*DTM0
               CALL MAGBRK(KSTAR(J2),MASS(2),MENV(2),RAD(2),OSPIN(2),
     &                     DJMB)
               DJSPIN(2) = DJMB*DTM0
* Check equilibrium spin and add tidal changes.
               DO K = 1,2
                  IF(DSPINK(K).GT.1.0D-10)THEN
                     DSPIN0 = MIN(DTM0*DSPINK(K),EQSPIN-OSPIN(K))/DTM0
                     DJTK(K) = DSPIN0*DJTK(K)/DSPINK(K)
                  ELSEIF(DSPINK(K).LT.-1.0D-10)THEN
                     DSPIN0 = MAX(DTM0*DSPINK(K),EQSPIN-OSPIN(K))/DTM0
                     DJTK(K) = DSPIN0*DJTK(K)/DSPINK(K)
                  ELSE
                     DSPINK(K) = 0.D0
                  ENDIF
                  DJORB = DJORB + DJTK(K)*DTM0
                  DJSPIN(K) = DJSPIN(K) - DJTK(K)*DTM0
               ENDDO
* Update spins and ensure that the star does not spin up beyond break-up.
               DO K = 1,2
                  JSPIN(K) = MAX(JSPIN(K) - DJSPIN(K),1.0D-10)
                  OSPBRU = TWOPI*SQRT(MASS(K)*AURSUN**3/RADX(K)**3)
                  JSPBRU = (K2STR(K)*(MASS(K)-MASSC(K))*RADX(K)*RADX(K)
     &                      + K3*MASSC(K)*RADC(K)*RADC(K))*OSPBRU
                  IF(JSPIN(K).GT.JSPBRU)THEN
                     DJORB = DJORB - (JSPIN(K) - JSPBRU)
                     JSPIN(K) = JSPBRU
                  ENDIF
                  OSPIN(K) = JSPIN(K)/
     &                      (K2STR(K)*RADX(K)*RADX(K)*(MASS(K)-MASSC(K))
     &                       + K3*RADC(K)*RADC(K)*MASSC(K))
               ENDDO
*     --02/28/13 9:58-lwang-add----------------------------------------*
***** Note:add --------------------------------------------------------**
*     Impose 2.5 % limit on change of JORB (Chris Tout 1/2013).
               IF (ABS(DJORB).GT.0.1*JORB.AND.ITRO.LE.5) THEN
                   DTM0 = 0.1*JORB*DTM0/ABS(DJORB)
                   DELET = 0.1*JORB*DELET/ABS(DJORB)
                   ITRO = ITRO + 1
                   IF (JORB.GT.0.0D0) WRITE (6,83)  DJORB/JORB, JORB
   83              FORMAT (' ROCHE LIMIT    DJO/JORB JORB ',1P,2E10.2)
                   GO TO 81
               ENDIF
*     --02/28/13 9:58-lwang-end-add------------------------------------*
* Adjust the eccentricity and orbital angular momentum.
               ECC = MAX(ECC - DELET,0.001D0)
               JORB = MAX(JORB - DJORB,0.D0)
* Reset the orbital separation and Roche-lobe radii.
               SEP = (MASS(1) + MASS(2))*JORB*JORB/
     &             ((MASS(1)*MASS(2)*TWOPI)**2*AURSUN**3*(1.D0-ECC*ECC))
               SEP = MAX(SEP,RAD(2))
*     --02/28/13 9:58-lwang-add----------------------------------------*
***** Note:add---------------------------------------------------------**
* Perform coalescence check after shrinkage.
               IF (RAD(1)+RAD(2).LT.SEP)THEN
                  COALS = .TRUE.
                  GO TO 60
               ENDIF
*     --02/28/13 9:59-lwang-end-add------------------------------------*
               DO K = 1,2
                  ROL(K) = RL(Q(K))*SEP
               ENDDO
               RADX(1) = MIN(RAD(1),ROL(1))
               RADX(1) = MAX(RADC(1),RADX(1))
*
               TB = (SEP/AURSUN)*SQRT(SEP/(AURSUN*(MASS(1)+MASS(2))))
               TK = DAYS*TB/YRS
               OORB = TWOPI/TB
               ITERB = ITERB + 1
               IF(ITERB.LE.100) GOTO 80
            ELSE
               IF(rank.eq.0.and.ECC.GT.0.02D0)THEN
                  WRITE(6,*)' WARNING: will not circularize'
               ENDIF
               IF(FAC.LT.0.98D0.AND.KSTAR(J1).LT.10)THEN
                  if(rank.eq.0)
     &            WRITE(6,*)' WARNING: will not synchronize'
               ENDIF
               ITERB = ITERB + 1
            ENDIF
         ENDIF
         IF(ITERB.GT.0)THEN
            if(rank.eq.0)
     &      WRITE(6,711)ECC,OSPIN(1)/OORB,SEP,TK,RAD(1),ROL(1)
  711       FORMAT(' NEW PARAMS    ECC SPIN1/OORB A P R1 RL1 ',
     &                             F7.4,F9.3,F8.2,1P,3E10.2)
            SPIN(J1) = JSPIN(1)/SPNFAC
            SPIN(J2) = JSPIN(2)/SPNFAC
*       Set new orbital elements (ECC, SEMI & JORB have changed).
            HI = H(IPAIR)
            SEMI = SEP/SU
*       Specify consistent orbital energy from new semi-major axis.
            H(IPAIR) = -0.5d0*BODY(I)/SEMI
*       Change the KS variables at constant eccentricity.
            CALL EXPAND(IPAIR,SEMI0)
*       Transform to pericentre (note TIME may change).
            CALL KSPERI(IPAIR)
*       Modify the eccentricity at constant energy.
            CALL DEFORM(IPAIR,ECC0,ECC)
            ZMU = BODY(J1)*BODY(J2)/BODY(I)
            EMDOT = EMDOT + ZMU*(HI - H(IPAIR))
            EGRAV = EGRAV + ZMU*(HI - H(IPAIR))
         ENDIF
*
      ENDIF
*
*       Check if components should be saved for output, cf. HRPLOT
      J = J1
      IF(TEV(J).GE.TPLOT.AND..NOT.ISAVE)THEN
         ISAVE = .TRUE.
*       Check if components are already saved (temporary)
         II = 0
         DO 700 K = 1,NRSAVE
            IF(NAME(J1).EQ.NAMER(1,K)) II = NAME(J1)
            IF(NAME(J1).EQ.NAMER(2,K)) II = NAME(J1)
  700    CONTINUE
         IF(II.EQ.0)THEN
            NRSAVE = NRSAVE + 1
            IF(NRSAVE.GT.MMAX)THEN
               IWARN = IWARN + 1
               IF (rank.eq.0.and.IWARN.LT.20) THEN
                   WRITE(6,*)' ROCHE WARNING: NRSAVE > MMAX'
               END IF
               NRSAVE = MMAX
            ENDIF
            NAMER(1,NRSAVE) = NAME(J1)
            NAMER(2,NRSAVE) = NAME(J2)
            KSTARR(1,NRSAVE) = KSTAR(I)
            KSTARR(2,NRSAVE) = KSTAR(J1)
            KSTARR(3,NRSAVE) = KSTAR(J2)
            CMRCH(1,NRSAVE) = SEMI
            DO 701 K = 1,3
               CMRCH(K+1,NRSAVE) = X(K,I)
               CMRCH(K+4,NRSAVE) = XDOT(K,I)
  701       CONTINUE
            CMRCH(8,NRSAVE) = BODY0(J1)
            CMRCH(9,NRSAVE) = BODY(J1)
            CMRCH(10,NRSAVE) = EPOCH(J1)
            CMRCH(11,NRSAVE) = BODY0(J2)
            CMRCH(12,NRSAVE) = BODY(J2)
            CMRCH(13,NRSAVE) = EPOCH(J2)
         ENDIF
      ENDIF
*
* See whether the primary still fills its Roche lobe.
*
      IF(RAD(1).GT.ROL(1))THEN
*
* Test for a contact system.
*
         IF (RAD(2).GT.ROL(2).AND.ITER.GT.0) THEN
            NCONT = NCONT + 1
            if(rank.eq.0)
     &      WRITE (6,72)  MASS(1), MASS(2), KSTAR(J1), KSTAR(J2), SEP
   72       FORMAT (' CONTACT    M1 M2 KW1 KW2 A ',2F7.3,2I3,F8.3)
*
*       Form new star by combining the components inelastically.
            IQCOLL = 3
            CH5 = ' CONT'
            GO TO 160
         END IF
*
         IF (IWARN.LT.100.AND.RAD(1).GT.2.0*ROL(1)) THEN
            IWARN = IWARN + 1
            if(rank.eq.0)
     &      WRITE (6,73)  KSTAR(J1), KSTAR(J2), MASS(1), MASS(2),
     &                    RAD(1), ROL(1), DM1, DM2, DTM
   73       FORMAT (' BIG ROCHE    K12 M12 R1 RL1 DM12 DT ',
     &                             2I3,2F7.3,2F8.3,1P,3E10.2)
          END IF
*
*       Produce diagnostics for cases involving degenerate objects.
          IF (ITER.EQ.1.AND.INEW.GT.0.AND.KSTAR(J2).GE.10) THEN
              IF (ABS(DTM).GT.1.0D-20) THEN
                  PD = DAYS*SEMI*SQRT(SEMI/BODY(I))
                  if(rank.eq.0)
     &            WRITE (22,58)  NAME(J1), NAME(J2),
     &                           KSTAR(J1), KSTAR(J2),
     &                           MASS(1), MASS(2), TPHYS, SEMI*SU, PD,
     &                           DM1/DTM/1.0E+06, DM2/DTM/1.0E+06
   58             FORMAT (' DEGEN ROCHE    NAM K* M TP A P M1D M2D ',
     &                             2I6,2I4,2F6.2,F8.1,2F7.1,1P,2E9.1)
                  CALL flush(22)
              END IF
          END IF
*
          ITER = ITER + 1
          IF(TEV0(I).LT.TIME) GOTO 10
          IF (ITER.EQ.1.AND.GAMMA(IPAIR).LT.GMIN) GO TO 10
      ELSE
          if(rank.eq.0)
     &    WRITE (6,76)  KSTAR(J1), KSTAR(J2), MASS(1),
     &                  MASS(2),RAD(1), ROL(1), DM1, DM2, DTM
   76     FORMAT(' END ROCHE    K12 M12 R1 RL1 DM DT ',
     &                          2I3,2F7.3,2F8.2,2F8.3,1P,E10.2)
*       Check optional diagnostics for degenerate objects.
          IF(MAX(KSTAR(J1),KSTAR(J2)).GE.10)THEN
             IF(KZ(8).GT.3)THEN
                CALL degen(IPAIR,IPAIR,4)
             ELSE
                NDD = NDD + 1
             ENDIF
          ENDIF
          IF (KSTAR(J1).GE.10.AND.KSTAR(J1).LE.12) NWD = NWD + 1
*       Define end of Roche stage.
          KSTAR(I) = KSTAR(I) + 1
          CALL trflow(IPAIR,DTR)
          TEV(I) = TIME + DTR
          IF (KZ(8).GT.3) THEN
              CALL binev(IPAIR)
          END IF
          CH5 = ' END '
          if(rank.eq.0)
     &    WRITE(85,95)NAME(J1),NAME(J2),KSTAR(J1),KSTAR(J2),KSTAR(I),
     &                TPHYS,AJ(1),AJ(2),MASS0(1),MASS0(2),
     &                MASS(1),MASS(2),ZMET,ECC,TK,JSPIN(1),JSPIN(2),CH5
          CALL FLUSH(85)
      END IF
*
*       Check standard exit condition.
      IF(ITER.EQ.0) GOTO 150
*
  140 CONTINUE
*       Set new c.m. time only during active stage (otherwise by TRFLOW).
      IF (MOD(KSTAR(I),2).EQ.1.AND.IQ.EQ.0) THEN
          TEV(I) = TPHYS0/TSTAR
*       Advance update time (extra for critical configurations).
          IF(LIST(1,I1).GT.0)THEN
             DTMAX = MIN(DTMAX,50.D0*DT1)
          ELSE
             DTMAX = MIN(DTMAX,10.D0*DT1)
          ENDIF
          DTMAX = MIN(DTMAX,DT2)
*         DTMAX = MIN(DTMAX,STEPX)
          TEV(I) = TEV(I) + MAX(DTMAX,1.0D-07)
          TEV(J1) = TEV(I) + 2.0*STEPX + STEP(I1)
          TEV(J2) = TEV(J1)
          TMDOT = MIN(TMDOT,TEV(I))
          Z = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
          ECC = SQRT(Z)
          IF (ECC.GT.0.01) THEN
              QP = SEMI*(1.0 - ECC)
              ICIRC = -1
              CALL tcirc(QP,ECC,I1,I2,ICIRC,TC)
              if(rank.eq.0)
     &        WRITE (6,142)  TIME+TOFF, NAME(J1), KSTAR(J1), KSTAR(J2),
     &                       LIST(1,I1), MASS(1), MASS(2), ECC,
     &                       SEMI*SU, TC
  142         FORMAT (' ECCENTRIC ROCHE    T NM K* NP M E A TC ',
     &                                     F8.1,I6,3I4,3F6.2,2F7.1)
              KSTAR(I) = 0
              TEV(I) = 1.0D+10
          END IF
      ELSEIF(IQ.EQ.0)THEN
*       Update look-up time.
          CALL trdot(J1,DTM,MASS(1))
          TEV(J1) = TEV(J1) + DTM
          CALL trdot(J2,DTM,MASS(2))
          TEV(J2) = MIN(TEV(J1),TEV(J2) + DTM)
*
*       Ensure that close degenerate system gets GR check in mdot. 
          KW = MIN(KSTAR(J1),KSTAR(J2))
          IF(KW.GT.6) TEV(J2) = MIN(TEV(J2),TIME+0.1D0)
*
          TEV(J1) = TEV(J2)
      END IF
*
*       Restore times and step and update mass on GRAPE.
      TPHYS = (TIME+TOFF)*TSTAR
*
*       Include consistency check on mass difference.
      DM = (BODY(J1) + BODY(J2)) - BODY(I)
      IF(ABS(DM).GT.0.0001*BODY(I))THEN
          if(rank.eq.0)
     &    WRITE (6,145)  TIME+TOFF, BODY(J1), BODY(J2), BODY(I)
  145     FORMAT (' DANGER!    ROCHE    T M ',F10.3,1P,3E12.4)
          BODY(I) = BODY(I) + DM
      ENDIF
*
      GO TO 200
*
  150 IF(IPHASE.GE.0) GOTO 140
  160 TPHYS = (TIME+TOFF)*TSTAR
*
*       Treat case of OVERFILL or CONTACT binary in the same way.
      IF(IQCOLL.EQ.3)THEN
*       Predict c.m. to highest order and save KS index in COMMON.
          if(rank.eq.0)
     &    WRITE(85,95)NAME(J1),NAME(J2),KSTAR(J1),KSTAR(J2),KSTAR(I),
     &                TPHYS,AJ(1),AJ(2),MASS0(1),MASS0(2),
     &                MASS(1),MASS(2),ZMET,ECC,TK,JSPIN(1),JSPIN(2),CH5
          IF(TIME-T0(I).LE.STEP(I))THEN
              CALL xvpred(I,0)
          ENDIF
          KSPAIR = IPAIR
          CALL cmbody(R(IPAIR),2)
          IQCOLL = 0
      ENDIF
*
  200 RETURN
*
      END
***
