      SUBROUTINE MDOT
*
*
*       Mass loss from evolving stars.
*       Updated 6/1/98 by J. Hurley
*       ------------------------------
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  BM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &               BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &               RP(NTMAX),ES(NTMAX),CM(2,NTMAX),IOSC(NTMAX),
     &               NAMEC(NTMAX)
      INTEGER JX(2),KACC
      REAL*8 MASS(2),MASSC(2),RAD(2),RADC(2),LUMIN(2)
      REAL*8 AGE0(2),TM0(2),TBGB(2),MENV(2),RENV(2),K2STR(2)
      REAL*8 TSCLS(20),LUMS(10),GB(10),TM,TN
      REAL*8 M0,M1,RM,AGE,LUM,MC,RCC,ME,RE,RNEW,LNEW,MCH
      PARAMETER (MCH=1.44d0)
      REAL*8 M10,RXL1,RXL2,Q,ROL(2),RLPERI,RM0
      REAL*8 EPS,ALPHA2,TOL
      PARAMETER(EPS=1.0d-06,ALPHA2=0.09d0,TOL=1.0d-10)
      REAL*8 DTX(2),DMS(2),DMA(2),DMX(2),EPCH0(2),TEVK,DIFF
      REAL*8 DT,DSEP,DJMB,DJGR,DJORB,DJT,DJTI,DTGR,DSPIN,JSPBRU,OSPBRU
      REAL*8 RJ,HJ,TJ2,VORB2,VWIND2,IVSQM,OMV2,DTXMIN,DMXMAX,DELET
      REAL*8 JSPIN(2),OSPIN(2),DJSPIN(2),JORB,OORB,AURSUN
      PARAMETER(AURSUN=214.95d0)
      REAL*8 K2,K3,ACC1,ACC2,BETA,XI
      PARAMETER(K3=0.21d0,ACC1=3.920659d+08,ACC2=1.5d0)
      PARAMETER(BETA=0.125d0,XI=1.d0,MAXBLK=40)
      REAL*8 MLWIND,RL
      EXTERNAL MLWIND,RL
      CHARACTER*8  WHICH1
      LOGICAL ICORR,IKICK
      SAVE IWARN
      DATA IWARN /0/
*
*       Define current time (unit of million years) and update counter.
      TTOT = TIME + TOFF
      TPHYS = TTOT*TSTAR
      IPOLY = 0
      TIME0 = TIME
      IBLUE = 0
*
*       Find global index of next star to be checked (NB! reset KS & IKS).
    1 KS = 0
      KSPAIR = 0
      IKS = 0
      ITRY = 0
      IGHOST = 0
      I0 = 0
      IHDOT = 0
      IF (IPHASE.LT.0) IPOLY = -1
      IPHASE = 0
      IQCOLL = 0
*       Restore synchronized TIME in case of RESET, CHRECT & KSPERI.
      TIME = TIME0
*
*       Determine smallest look-up time (< TIME!).
      DO 5 J = 1,NTOT
          IF (TEV(J).LE.TIME) THEN
              I = J
              GO TO 10
          END IF
    5 CONTINUE
*
*       Determine new evolution time (current TMDOT may be escaped star).
      KW = 0
      I = 1
      GO TO 70
*
*       Include braking or merger termination at time for Roche overflow.
   10 IF (I.GT.N.AND.NAME(I).LT.0.AND.KZ(34).GT.0) THEN
          KSPAIR = I - N
          CALL BRAKE2(KSPAIR,ITERM)
          IF(ITERM.GT.0)THEN
             NWARN = NWARN + 1
             IF (rank.eq.0.and.NWARN.LT.1000) THEN
                 WRITE(14,900) KSPAIR, NAME(I), TTOT, R(KSPAIR)
 900             FORMAT(' MDOT MERGER TERM:    KS NAM T R ',
     &                                         I4,I6,F10.4,1P,E10.2)
                 CALL FLUSH(14)
             END IF
             JCOMP = 0
             IPHASE = 7
             CALL RESET
          ENDIF
          GO TO 1
      END IF
*
*       Check for hierarchy containing NS/BH (same TEV for inner components).
      IF (I.LT.IFIRST.AND.KZ(28).GT.0) THEN
          IP = KVEC(I)
          IF (NAME(N+IP).LT.0.AND.KSTAR(I).GE.13) THEN
*       Distinguish between inner binary components and outer component(s).
              IF (I.EQ.2*IP-1) THEN
                  CALL BRAKE3(IP,ITERM)
*       Update single outer KS member or terminate rare case of binary.
              ELSE IF (NAME(2*IP).LE.NZERO) THEN
                  TEV0(I) = TIME
                  TEV(I) = TIME + 100.0
                  ITERM = 0
              ELSE
                  ITERM = 1
              END IF
              IF (ITERM.GT.0) THEN
                  IPHASE = 7
                  KSPAIR = IP
                  CALL RESET
              END IF
              GO TO 1
          END IF
      END IF
*
*       Check possible Roche overflow condition (skip any ghost c.m.).
      IF (I.GT.N.AND.KZ(34).GT.0.AND.BODY(I).GT.0.0D0.AND.
     &    KSTAR(I).GE.10) THEN
          IPAIR = I - N
*       Consider odd types (11, 13, etc.) for Roche updating or treatment.
          IF (MOD(KSTAR(I),2).NE.0) THEN
              CALL ROCHE(IPAIR)
              IF (IQCOLL.NE.0.OR.IPHASE.LT.0) GO TO 1
*       Ensure that look-up time is increased in case of no mass transfer.
              IF (TEV(I).LE.TIME) THEN
                 IF (MOD(KSTAR(I),2).NE.0) KSTAR(I) = KSTAR(I) + 1
                 CALL TRFLOW(IPAIR,DTR)
                 TEV(I) = TIME + DTR
                 if(rank.eq.0)
     &           WRITE(6,902)I,NAME(I),IPAIR,KSTAR(I),TEV(I),TTOT
 902             FORMAT(' WARNING   MDOT NO ROCHE ',3i6,i4,1x,1p,2e13.4)
              ENDIF
          ELSE
*       Include other cases (even types outside 10 R_sun).
              CALL TRFLOW(IPAIR,DTR)
              IF (DTR.LT.STEP(I)) THEN
                  TEV(I) = TIME + DTR
                  CALL ROCHE(IPAIR)
                  IF (TEV(I).LE.TIME) TEV(I) = 1.000002d0*TIME
              ELSE
                  IF(KZ(34).EQ.1)THEN
                      CALL SYNCH(IPAIR)
                  ELSE
                      TEV(I) = TIME + DTR
                  ENDIF
              END IF
          END IF
          GO TO 1
      END IF
*
*       Check for special treatment of possible hierarchical system (IHDOT).
      IF (I.LT.IFIRST) THEN
          KSPAIR = KVEC(I)
          KSX = MAX(KSTAR(2*KSPAIR-1),KSTAR(2*KSPAIR))
          ICM = N + KSPAIR
          I2 = 2*KSPAIR
          IF(I2.EQ.I) I2 = I2 - 1
*
*       Select merger if #I is first or the second member is a c.m. body.
          IF (NAME(ICM).LT.0.AND.
     &       (I.EQ.2*KSPAIR - 1.OR.NAME(2*KSPAIR).GT.NZERO)) THEN
*       Activate indicator (distinguish inner and outer binary).
              IHDOT = 1
*       Exclude inner components and outer binary of double hierarchy.
              IF (NAME(ICM).LT.-2*NZERO) THEN
                  IF (I.LT.2*KSPAIR.OR.NAME(I).GT.NZERO) THEN
                      IPHASE = 7
                      CALL RESET
                      GO TO 1
                  END IF
              END IF
          END IF
*
*       Include Roche consistency check (single components not permitted).
          IF (KSTAR(ICM).GT.10.AND.MOD(KSTAR(ICM),2).EQ.1) THEN
              IF(BODY(I).EQ.0.0.OR.BODY(I2).EQ.0.0)THEN
                 TEV(I) = TIME + MIN(0.1D0,TEV(I)-TEV0(I))
                 TEV(I2) = TEV(I)
                 GOTO 1
              ENDIF
              IWARN = IWARN + 1
              IF (rank.eq.0.and.MOD(IWARN,100).EQ.0) THEN
                 WRITE(6,904) I, KSTAR(ICM), TEV(I) - TEV(ICM)
 904             FORMAT(' WARNING!    MDOT    I K* D(TEV) ',2I5,1P,E9.1)
              END IF
              CALL TRFLOW(KSPAIR,DTR)
              TEV(ICM) = TIME + DTR
              IF (DTR.GT.TSTAR) THEN
                 KSTAR(ICM) = KSTAR(ICM) + 1
              ELSE
                 TEV(I) = 1.000002d0*TIME
                 GO TO 1
              END IF
          END IF
      END IF
*
*       Determine relevant indices in case of merger (zero mass done below).
      IF (IHDOT.GT.0.AND.BODY(I).GT.0.0D0) THEN
*
*       Identify ghost address and merger index from c.m. of KSPAIR.
          CALL FINDJ(I,IGHOST,IMERGE)
*
          IF(rank.eq.0.and.IGHOST.LT.0)THEN
             WRITE(6,*)' MDOT GHOST NOT FOUND ',I,NAME(I),IMERGE
             STOP
          ENDIF
*
*       Determine index of mass-losing star for inner or outer binary.
          IF (IGHOST.LE.N) THEN
              IF (TEV(IGHOST).LT.TEV(I)) I = IGHOST
              IF (rank.eq.0.and.I.EQ.IGHOST) WRITE(6,906) I, NAME(I),
     &                                       KSTAR(I), TEV(I), BODY(I)
 906          FORMAT(' GHOST SWITCH:    I NAM K* TEV BODY ',
     &                                   3I5,F9.2,1P,E10.2)
          ELSE
*       Save KS component for updating look-up time.
              I0 = I
*       Compare look-up times of #I and IGHOST <= N on first KS component.
              IF (I.LT.2*KSPAIR) THEN
                  IF (TEV(IGHOST).LT.TEV(I)) I = IGHOST
              ELSE
*       Form ghost pair index and select KS component by comparing TEV.
                  JPAIR = IGHOST - N
                  I = 2*JPAIR - 1
                  IF (TEV(2*JPAIR).LT.TEV(I)) I = 2*JPAIR
                  IHDOT = 2
                  if(rank.eq.0)
     &            WRITE(6,908) IGHOST, I0, I, NAME(I), KSTAR(I), TEV(I)
 908              FORMAT(' OUTER GHOST:    IG I0 I NM K* TEV',5I6,F9.2)
                  IF(KSTAR(N+KSPAIR).GE.10) THEN
                     CALL RESET
                     GOTO 1
                  ENDIF
              END IF
          END IF
      ELSE IF (BODY(I).LE.0.0D0.OR.NAME(I).EQ.0) THEN
*       Distinguish between merger ghost and member of CHAIN/TRIPLE/QUAD.
          IMERGE = 0
          ICM = 0
          J = I
*       Replace #I by corresponding ghost c.m. for zero-mass KS component.
          IF (I.LT.IFIRST) J = N + KVEC(I)
*
*       Identify merger index and determine corresponding c.m. body.
          DO 15 K = 1,NMERGE
              IF (NAMEG(K).EQ.NAME(J)) IMERGE = K
   15     CONTINUE
          IF (IMERGE.GT.0) THEN
              DO 20 K = N+1,NTOT
                  IF (NAMEM(IMERGE).EQ.NAME(K)) ICM = K
   20         CONTINUE
              KSPAIR = ICM - N
              if(icm.eq.0) kspair = j - n
*       Terminate if ghost belongs to double hierarchy (bug fix 03/2001).
              IF (NAMEM(IMERGE).LT.-2*NZERO) THEN
                  IPHASE = 7
                  CALL RESET
                  GO TO 1
              END IF
          END IF
*
*       Define KS index and set indicator on successful identification.
          IF (ICM.GT.N) THEN
              IHDOT = 1
              IGHOST = J
              IF (IGHOST.GT.N) THEN
                  IHDOT = 2
                  JPAIR = IGHOST - N
                  I = 2*JPAIR - 1
                  IF (TEV(2*JPAIR).LT.TEV(I)) I = 2*JPAIR
                  I0 = IGHOST
              END IF
          ELSE
*       Extend TEV for other single ghosts or c.m. of compact subsystem.
              TEV(I) = TEV(I) + 0.01d0
*       Allow an extra amount for MS stars (depending on TEV - TEV0).
              IF(KSTAR(I).LE.1.AND.(TEV(I)-TEV0(I))*TSTAR.LT.10.0)THEN
                  TEV(I) = TEV(I) + 0.5*(TEV(I) - TEV0(I))
              END IF
              GO TO 1
          END IF
      END IF
*
*       Skip any c.m. particles after increasing look-up time.
      IF (I.GT.N) THEN
          I1 = 2*(I - N) - 1
          I2 = I1 + 1
          TM = MIN(TEV(I1),TEV(I2))
*         IF(KSTAR(I).GE.10)THEN
*            IWARN = IWARN + 1
*            IF (rank.eq.0.and.MOD(IWARN,100).EQ.0) THEN
*               WRITE(6,910) I, NAME(I), KSTAR(I1), KSTAR(I2), KSTAR(I),
*    &                       TTOT, TM - TEV(I)
*910            FORMAT(' WARNING!    MDOT:    I NAM K* T D(TEV) ',
*    &                               2I6,3I4,F9.2,1P,E10.2)
*            END IF
*         ENDIF
          IF (KSTAR(I).LT.10) THEN
              TEV(I) = 1.0d+10
          ELSE IF (KZ(34).GT.0) THEN
              IPAIR = I - N
              CALL TRFLOW(IPAIR,DTR)
              DTR = MAX(DTR,0.0D0)
              IF (DTR.LT.STEP(I)) THEN
                  TEV(I) = TIME + DTR
                  CALL ROCHE(IPAIR)
              END IF
              TEV(I) = MAX(TIME + DTR,1.000002d0*TIME)
          ELSE
              TEV(I) = 1.000002d0*TIME
          END IF
          GO TO 1
      END IF
*
*       Copy relevant mass (standard case or merger ghost member).
*       Also determine if the star is in a standard binary or is the inner
*       component of a hierachy and if so set the mass accretion flag.
      KACC = 0
      JX(1) = I
      ECC2 = 0.D0
      RJ = 0.D0
      HJ = 0.D0
      TJ2 = 0.D0
      IF(IGHOST.EQ.0.AND.IHDOT.EQ.0)THEN
         MASS(1) = BODY(I)*ZMBAR
         IF(I.LT.IFIRST)THEN
            IF(NAME(ICM).GT.0)THEN
               KACC = 1
               JX(2) = I2
               MASS(2) = BODY(I2)*ZMBAR
               RJ = R(KSPAIR)
               HJ = H(KSPAIR)
               TJ2 = TDOT2(KSPAIR)
            ENDIF
         ENDIF
      ELSE
         IF(IHDOT.EQ.1)THEN
            K = 1
            IF (I.GE.IFIRST) K = 2
            IF(I.NE.2*KSPAIR)THEN
               KACC = 2
               J = 3 - K
               JX(2) = IGHOST
               IF(K.EQ.2) JX(2) = 2*KSPAIR - 1
               MASS(2) = BM(J,IMERGE)*ZMBAR
               DO 25 II = 1,4
                  RJ = RJ + UM(II,IMERGE)*UM(II,IMERGE)
                  TJ2 = TJ2 + 2.d0*UM(II,IMERGE)*UMDOT(II,IMERGE)
  25           CONTINUE
               HJ = HM(IMERGE)
            ENDIF
         ELSE
            K = 3
            IF (I.EQ.2*JPAIR) K = 4
         END IF
         MASS(1) = BM(K,IMERGE)*ZMBAR
      END IF
      IF(KACC.GT.0)THEN
         J = N + KSPAIR
*       Check for zero c.m. mass (i.e. NAME < 0 & > -NZERO in [[B,S],[B,S]]).
         IF (BODY(J).EQ.0.0D0) THEN
             IPHASE = 7
             CALL RESET
             GO TO 1
         END IF
         SEMI = -0.5d0*BODY(J)/HJ
         SEP = SEMI*SU
         ECC2 = (1.d0 - RJ/SEMI)**2 + TJ2**2/(BODY(J)*SEMI)
         IF(ECC2.LT.1.0)THEN
            KACC = 2
            Q = MASS(1)/MASS(2)
            ROL(1) = RL(Q)*SEP
            RXL1 = RADIUS(JX(1))/ROL(1)
            Q = 1.d0/Q
            ROL(2) = RL(Q)*SEP
            RXL2 = RADIUS(JX(2))/ROL(2)
            ISWAP = 0
            IF(KSTAR(JX(1)).LE.1.AND.
     &         (KSTAR(JX(2)).EQ.5.OR.KSTAR(JX(2)).EQ.6))THEN
               ISWAP = 1
            ELSEIF(RXL2.GT.RXL1)THEN
               ISWAP = 1
            ENDIF
            IF(ISWAP.GT.0)THEN
               J1 = JX(1)
               JX(1) = JX(2)
               JX(2) = J1
               M10 = MASS(1)
               MASS(1) = MASS(2)
               MASS(2) = M10
               RXL1 = ROL(1)
               ROL(1) = ROL(2)
               ROL(2) = RXL1
            ENDIF
            VORB2 = ACC1*(MASS(1)+MASS(2))/SEP
            IVSQM = 1.D0/SQRT(1.D0-ECC2)
            DIFF = ABS(TEV0(JX(2)) - TEV0(JX(1)))
            IF(DIFF.GT.TOL)THEN
               TEV(JX(1)) = MAX(TEV0(JX(1)),TEV0(JX(2)))
            ELSE
               TEV0(JX(2)) = TEV0(JX(1))
               IF(TEV(JX(2)).NE.TEV(JX(1)))THEN
                  TEV(JX(1)) = MIN(TEV(JX(1)),TEV(JX(2)))
               ENDIF
            ENDIF
            TEV(JX(2)) = TEV(JX(1))
            DIFF = ABS(TEV0(JX(2)) - TEV0(JX(1)))
            ECC = SQRT(ECC2)
            OORB = TWOPI*SQRT((MASS(1)+MASS(2))/(SEP/AURSUN)**3)
            JORB = MASS(1)*MASS(2)/(MASS(1)+MASS(2))
     &             *SQRT(1.D0-ECC2)*SEP*SEP*OORB
         ELSE
            if(rank.eq.0)
     &      WRITE(6,912)JX,NAME(N+KSPAIR),SQRT(ECC2),TTOT
 912        FORMAT(' MDOT ECC2 > 1.0 I1 I2 N(ICM) e T ', 3I6,F6.2,E9.2)
            KACC = 1
         ENDIF
      ELSE
         KACC = 1
      ENDIF
*
      NMDOT = NMDOT + 1
*
      DO 200 K = 1,KACC
*
         I = JX(K)
*       Set interval since last mass update.
         DTX(K) = 1.0D+06*(TEV(I) - TEV0(I))*TSTAR
*       Set the initial mass and current type.
         M0 = BODY0(I)*ZMBAR
         M1 = MASS(K)
         MC = 0.D0
         KW = KSTAR(I)
*
*       Obtain stellar parameters at previous epoch.
         AGE = TEV0(I)*TSTAR - EPOCH(I)
         AGE = MAX(AGE,0.D0)
         AGE0(K) = AGE
         EPCH0(K) = EPOCH(I)
         CALL star(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
         CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               RM0,LUM,KW,MC,RCC,ME,RE,K2)
         RAD(K) = RM0
         LUMIN(K) = LUM
         TM0(K) = TM
         TBGB(K) = TSCLS(1)
         MASSC(K) = MC
         RADC(K) = RCC
         MENV(K) = ME
         RENV(K) = RE
         K2STR(K) = K2
*
*       Ensure that type change occurs at time TEV.
         IF(KW.NE.KSTAR(I).OR.DABS(M1-M0)/M0.GT.0.05D0)THEN
            IF (rank.eq.0.and.KW.GE.13) THEN
                WRITE (6,190)  I, NAME(I), KSTAR(I), KW, MASS(K)
  190           FORMAT (' NS/BH FORMATION    ',2I7,2I4,F7.2)
            END IF
            KWOLD = KW
            KW = KSTAR(I)
            M1 = MASS(K)
           
*        Start new full sev output
            if(rank.eq.0)then
            RI = SQRT(X(1,I)**2 + X(2,I)**2 + X(3,I)**2)
            VI = SQRT(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2)
            write(111,191) 0,0,name(i),ttot,m0,radius(i),
     *          kw,kwold,ri,vi,age0(k),mass(k),massc(k),menv(k),
     *          rad(k),radc(k),renv(k),lumin(k),tm0(k)
 191        format(1x,2i5,i8,1p,3e12.5,2i4,11e12.5)
            end if
         ENDIF
*
*       Evaluate mass loss due to stellar wind.
         RLPERI = 0.D0
         IF(KACC.GT.1) RLPERI = ROL(K)*(1.D0 - ECC)
         DMX(K) = MLWIND(KW,LUM,RM0,M1,MC,RLPERI,ZMET)
         DMA(K) = 0.D0
*
         IF(KACC.GT.1)THEN
*
*       Calculate how much of wind mass loss will be accreted by the
*       companion (Boffin & Jorissen, A&A 1988, 205, 155).
            VWIND2 = 2.D0*BETA*ACC1*M1/RM0
            OMV2 = (1.d0 + VORB2/VWIND2)**(3.D0/2.D0)
            DMA(3-K) = IVSQM*ACC2*DMX(K)*((ACC1*MASS(3-K)/VWIND2)**2)
     &                  /(2.D0*SEP*SEP*OMV2)
            DMA(3-K) = MIN(DMA(3-K),0.8d0*DMX(K))
*
         ENDIF
*
*       Convert the spin angular momentum of the star to physical units.
         JSPIN(K) = MAX(SPIN(I)*SPNFAC,1.0D-10)
*
*       Evaluate the spin of the star.
         OSPIN(K) = JSPIN(K)/(K2*RM0*RM0*(M1-MC)+K3*RCC*RCC*MC)
*
 200  CONTINUE
*
      DT = DTX(1)
*
      DO 210 K = 1,KACC
*
         I = JX(K)
         KW = KSTAR(I)
*
*       Check for unequal evolution times in a recently formed binary.
         IF(K.EQ.2)THEN
            IF(DTX(1).EQ.0.D0.OR.DTX(2).EQ.0.D0) DT = DTX(2)
         ENDIF
*
*       Calculate the change in spin angular momentum due to mass loss and
*       also due to mass accretion and/or tides if the star is in a binary.
         DJSPIN(K) = (2.D0/3.D0)*DMX(K)*RAD(K)*RAD(K)*OSPIN(K)
*       Restrict the time-step for extreme conditions (JH 05/09).
         IF(DJSPIN(K).GT.0.0D0)THEN
            DT = MIN(DT,0.3D0*JSPIN(K)/DJSPIN(K))
         ENDIF
*
         IF(KACC.GT.1)THEN
*
            IF(DMA(K).GT.0.D0)THEN
               DJSPIN(K) = DJSPIN(K) - (2.d0/3.d0)*XI*DMA(K)*
     &                     RAD(3-K)*RAD(3-K)*OSPIN(3-K)
*
*       Produce diagnostics for symbiotic stars.
               IF(rank.eq.0.and.DMA(K).GE.3.16E-9.AND.KW.GE.10)THEN
                  WRITE(20,913)NAME(I),NAME(JX(3-K)),KSTAR(I),
     &                     KSTAR(JX(3-K)),TPHYS,MASS,SEP,DMX(3-K),DMA(K)
 913              FORMAT(' WINDAC    NAM K* T M A MLOSS MACC ',
     &                               2I6,2I4,F9.2,2F6.2,F8.1,1P,2E9.1)
                  CALL FLUSH(20)
               ENDIF
            ENDIF
*
*       Tidal synchronisation (only operates on circular orbits).
            Q = MASS(3-K)/MASS(K)
            IF(K.EQ.1) DJT = 0.D0
            IF(KZ(34).EQ.2.AND.DT.GT.0.D0.AND.ECC.LE.0.01D0)THEN
               IF((KW.LE.9.AND.RAD(K).GE.0.01D0*ROL(K)).OR.
     &            (KW.GE.10.AND.Q.GE.1.D0))THEN
*
                  CALL BSETID(kw,mass(k),massc(k),menv(k),rad(k),
     &                       radc(k),renv(k),lumin(k),ospin(k),k2str(k),
     &                       q,sep,ecc,oorb,delet,dspin,eqspin,djti)
*
                  IF(ABS(DSPIN).GT.TINY)THEN
                     DSPIN0 = DSPIN
                     IF(DSPIN.GE.0.D0)THEN
                        DSPIN = MIN(DT*DSPIN,OORB-OSPIN(K))/DT
                     ELSE
                        DSPIN = MAX(DT*DSPIN,OORB-OSPIN(K))/DT
                     ENDIF
                     DJTI = DJTI*DSPIN/DSPIN0
                  ELSE
                     DJTI = 0.D0
                  ENDIF
                  DJSPIN(K) = DJSPIN(K) - DJTI
                  DJT = DJT + DJTI
*
               ENDIF
            ENDIF
         ENDIF
*
*       Include magnetic braking for stars with convective envelopes.
         CALL MAGBRK(KW,MASS(K),MENV(K),RAD(K),OSPIN(K),DJMB)
*       Limit to a 3% angular momentum change for the star owing to MB.
         IF(DJMB.GT.0.0D0)THEN
            DT = MIN(DT,0.03D0*JSPIN(K)/DJMB)
            DJSPIN(K) = DJSPIN(K) + DJMB
         ENDIF
*
*       Evaluate the mass loss from the star in the interval DT.
         DMS(K) = (DMX(K) - DMA(K))*DT
         DMR = ABS(DMS(K)/(MASS(K) + 1.0d-10))
*
*       Restrict accumulated mass loss to maximum of 2%.
         IF(DMR.GT.0.02)THEN
            DT = DT*0.02/DMR
            DMS(K) = 0.02*MASS(K)
         ENDIF
*
*       Check that mass loss does not exceed the envelope mass.
         IF(KSTAR(I).LT.10)THEN
            DML = MAX(MASS(K) - MASSC(K),1.0D-07)
            IF(DML.LT.DMS(K))THEN
               DT = (DML/DMS(K))*DT
               DMS(K) = DML
            ENDIF
         ENDIF
         DTX(K) = DT
*
 210  CONTINUE
*
*       Include angular momentum loss owing to mass loss and/or
*       gravitational radiation for close binary systems.
      DSEP = 0.D0
      DJORB = 0.D0
      DTGR = 2.0D+10
      IF(KACC.EQ.2)THEN
*
         DTXMIN = MIN(DTX(1),DTX(2))
         CALL GRRAD(MASS(1),MASS(2),SEP,ECC,JORB,DJGR,DELET)
         DJORB = DJT + DJGR
*        if(rank.eq.0)
*    &   write(*,*)' grrad ',sep,ecc,djgr,delet,dtxmin
*
*       Limit orbital angular momentum change to 2%.
*
         IF(ABS(DJORB).GT.TINY.AND.(DTXMIN.GT.0.D0.OR.DIFF.EQ.0.D0))THEN
            DTGR = 0.02D0*JORB/ABS(DJORB)
            IF(DELET.GT.TINY.AND.ECC.GT.0.0011D0)THEN
               DTGR = MIN(DTGR,0.05D0*ECC/DELET)
            ENDIF
            DTGR = MAX(DTGR,100.D0)
*           DTGR = MAX(DTGR,1.0D-04)
            DTXMIN = MIN(DTXMIN,DTGR)
            DTGR = DTGR/1.0D+06
            DJORB = DJORB*DTXMIN
            JORB = JORB - DJORB
            JORB = MAX(JORB,1.D0)
            ECC0 = ECC
            ECC = MAX(ECC0 - DELET*DTXMIN,0.001D0)
            ECC2 = ECC*ECC
            SEP1 = (MASS(1) + MASS(2))*JORB*JORB/
     &             ((MASS(1)*MASS(2)*TWOPI)**2*AURSUN**3*(1.D0-ECC2))
            DSEP = SEP - SEP1
            IF(DSEP.GT.0.D0)THEN
               Q = MASS(1)/MASS(2)
               RXL1 = 0.9D0*RAD(1)/RL(Q)
               SEP1 = MAX(SEP1,RXL1)
               DSEP = SEP - SEP1
            ENDIF
            DTX(1) = DTXMIN
            DTX(2) = DTX(1)
         ELSE
            DJORB = 0.D0
         ENDIF
*
*       Orbital changes owing to mass loss dealt with in hcorr for now.
         DMXMAX = MAX(DMX(1),DMX(2))
*        IF(DMXMAX.GT.0.D0)THEN
*           Q = MASS(1)/MASS(2)
*           DJORB = DJORB + (DMX(1)+Q*DMA(1))*MASS(2)*MASS(2)*DTX(1)
*           Q = 1.D0/Q
*           DJORB = DJORB + (DMX(2)+Q*DMA(2))*MASS(1)*MASS(1)*DTX(2)
*           DJORB = DJORB*SEP*SEP*OORB/(MASS(1)+MASS(2))**2
*        ENDIF
*
      ENDIF
*
      DO 220 K = 1,KACC
*       Set the initial mass and current type.
         I = JX(K)
         DT = DTX(K)
         M0 = BODY0(I)*ZMBAR
         M10 = M0
         M1 = MASS(K)
         MC = MASSC(K)
         KW = KSTAR(I)
*       Set indicator for mass loss correction.
         ICORR = .FALSE.
         DMS(K) = (DMX(K) - DMA(K))*DT
         DMR = ABS(DMS(K)/(MASS(K) + 1.0d-10))
         IF(DMR.GT.TINY)THEN
            ICORR = .TRUE.
            M1 = M1 - DMS(K)
*       Check rejuvenation of MS, HG or HE star.
            IF(KW.LE.2.OR.KW.EQ.7)THEN
               M0 = M1
               CALL STAR(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
               IF(KW.EQ.2)THEN
                  IF(GB(9).LT.MASSC(K).OR.M10.GT.ZPARS(3))THEN
                     M0 = M10
                  ELSE
                     EPCH0(K) = TM + (TSCLS(1)-TM)*(AGE0(K)-TM0(K))/
     &                               (TBGB(K)-TM0(K))
                     EPCH0(K) = TEV0(I)*TSTAR - EPCH0(K)
                  ENDIF
               ELSE
                  EPCH0(K) = TEV0(I)*TSTAR - AGE0(K)*TM/TM0(K)
               ENDIF
            ENDIF
         ENDIF
         DJSPIN(K) = DJSPIN(K)*DT
*
*       Check for blue straggler formation (TM < TPHYS).
         IF(DMR.GT.TINY.AND.KW.LE.1.AND.NAME(I).NE.IBLUE)THEN
            TPHYS2 = (TEV0(I)+TOFF)*TSTAR - EPOCH0
            IF(TM0(K).GT.TPHYS2.AND.TM.LT.0.98D0*TPHYS2)THEN
               J = JX(1)
               IF(I.EQ.J.AND.KACC.EQ.2) J = JX(2)
               if(rank.eq.0)
     &         WRITE(6,914)NAME(I),M1,TM,TPHYS2,EPCH0(K)+TOFF*TSTAR,
     &                     KSTAR(J)
 914           FORMAT(' NEW BS (MDOT):   NAM M TM TP EP KW1 ',
     &                                   I6,F6.2,3F8.1,I4)
               IBLUE = NAME(I)
               NBS = NBS + 1
            ENDIF
         ENDIF
*
*       Set current age to actual look-up value (allows looping).
         TEVK = TEV0(I) + DT/(1.0D+06*TSTAR)
*       Set indicator for mass loss correction.
         AGE = TEVK*TSTAR - EPCH0(K)
         AGE0(K) = AGE
*
*       Determine stellar evolution time scales and luminosity.
         CALL star(KW,M0,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
*
*       Include temporary error check.
         IF((AGE-TN).GT.1.0d-02)THEN
*           IF(I.GE.IFIRST.OR.(I.LT.IFIRST.AND.NAME(ICM).GT.0))THEN
               if(rank.eq.0)
     &         WRITE(6,994)NAME(I), KW, DMS(K), AGE, TN
*           ENDIF
 994        FORMAT(' MDOT WARNING! AGE > TN   NM KW DMS AGE TN ',
     &                                   I6,I4,F7.3,1P,2E9.1)
            IF(KW.LE.6)THEN
               AGE = MIN(AGE,0.9999D0*TSCLS(11))
            ELSE
               AGE = MIN(AGE,0.9999D0*TSCLS(5))
               AGE = MIN(AGE,1.0001D0*TN)
            ENDIF
         ENDIF
*##
*       Obtain stellar parameters at current epoch.
         CALL hrdiag(M0,AGE,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &               RM,LUM,KW,MC,RCC,ME,RE,K2)
*
         IF(rank.eq.0.and.KW.EQ.14)PRINT*,' MDOT: T,M0,MC,AGE=',
     &                TTOT,M0,MC,AGE
*
*       Check for change of CO WD to ONe on rapid accretion (21/11/08).
         IF(DMR.GT.TINY.AND.KW.EQ.11)THEN
            DME = 2.08D-03*(1.D0/(1.D0 + ZPARS(11)))*RM
            IF(ABS(DMA(K)).GT.0.4d0*DME) KW = 12
         ENDIF
*
*       Check various aspects related to change of type.
         IF(KW.NE.KSTAR(I))THEN
            IF(KW.EQ.15)THEN
               if(rank.eq.0)
     &         WRITE(6,915)I,NAME(I),M0,KSTAR(I),TEVK*TSTAR-EPCH0(K)
 915           FORMAT (' MDOT WARNING! SN KW=15   I NM M0 KSTAR AGE ',
     &                                            2I6,F7.3,I4,E9.1)
               IF (NAME(N+KSPAIR).LT.0) THEN
                  IPHASE = 7
                  CALL RESET
                  GO TO 1
               ENDIF
            ELSEIF(KW.GE.13)THEN
*
*       Force new NS or BH to have a period of one second.
               OSPIN(K) = 2.0D+08
               JSPIN(K) = K3*RCC*RCC*MC*OSPIN(K)
            ENDIF
*       Check inclusion of mass loss on type change.
*       Note DM should be zero if ICORR not already true and no SN.
            DMS(K) = MASS(K) - M1
            DMR = ABS(DMS(K)/(MASS(K) + 1.0d-10))
            ICORR = .TRUE.
         ENDIF
*
*       Set mass loss and new radius in N-body units.
         DMSUN = DMS(K)
         DM = DMS(K)/ZMBAR
         RNEW = RM/SU
         LNEW = LUM
         KW0 = KSTAR(I)
*
*       Check mass-loss treatment for inner binary components of merger.
         IF (IHDOT.GT.0) THEN
            IF((KW.EQ.KSTAR(I).AND.AGE.GE.0.99*TN).OR.
     &         (KSTAR(I).LE.6.AND.KW.GT.6).OR.
     &         (KSTAR(I).LE.9.AND.KW.GT.9)) THEN
               ITERM = 1
            ELSE
               IF (IHDOT.EQ.1) THEN
                  CALL HMDOT(I,IMERGE,M1,KW,MC,DMSUN,RNEW,ITERM)
               ELSE
                  CALL HMDOT2(I,IGHOST,IMERGE,M1,KW,MC,DMSUN,RNEW,ITERM)
               END IF
            ENDIF
*       Terminate on non-zero indicator.
            IF (ITERM.GT.0) THEN
               IPHASE = 7
               CALL RESET
               GO TO 1
            END IF
         END IF
*
*       Define KS index & Roche indicator and update core mass (chaos only).
         IF (I.LT.IFIRST.AND.IHDOT.EQ.0) THEN
            KSPAIR = KVEC(I)
            IF (.NOT.ICORR) THEN
               IF (NAME(N+KSPAIR).GT.0) ITRY = 2
            END IF
            IF (KZ(27).GT.1) THEN
               IC = 0
*       Search existing chaotic binaries (assign next index if none found).
               DO 30 KK = 1,NCHAOS
                  IF (NAMEC(KK).EQ.NAME(N+KSPAIR)) IC = KK
   30          CONTINUE
               IF (IC.EQ.0) IC = NCHAOS + 1
*       Save core mass in chaos variable for the appropriate KS component.
               KK = 1
               IF (I.EQ.2*KSPAIR) KK = 2
*       Note that KW <= 1 ensures initialization by MC = 0.
               CM(KK,IC) = MC/ZMBAR
            END IF
         END IF
*
*       Include special procedures for KS components.
         IF (I.LT.IFIRST.AND.ICORR.AND.IHDOT.EQ.0) THEN
*       Distinguish between KS pair and merger configuration.
            SEMI = -0.5d0*BODY(N+KSPAIR)/H(KSPAIR)
            IF (NAME(N+KSPAIR).GT.0) THEN
*       Set random phase for neutron star or BH formation (negative index).
               IF (KSTAR(I).LT.10.AND.KW.GE.10) THEN
                  JPAIR = -KSPAIR
                  RADIUS(I) = RNEW
                  KSTAR(I) = -KSTAR(I)
                  CALL KSAPO(JPAIR)
               END IF
*       Terminate for large mass loss or soft binary and re-determine index.
               IF ((DMR.GT.0.2.AND.R(KSPAIR).GT.RMIN).OR.
     &            (DM.GT.0.0.AND.H(KSPAIR) + DM/SEMI.GT.-ECLOSE.AND.
     &            KSTAR(N+KSPAIR).GE.0).OR.
*    &            (KW.NE.KSTAR(I).AND.KW.GE.13)) THEN
     &            (KW.NE.KSTAR(I).AND.
     &            (KW.GE.13.OR.(KW.GE.11.AND.KZ(25).GE.1).OR.
     &            (KW.EQ.10.AND.KZ(25).GT.1)))) THEN
                  I = I + 2*(NPAIRS - KSPAIR)
                  JX(K) = I
                  IF(KACC.EQ.2)THEN
                     JX(3-K) = JX(3-K) + 2*(NPAIRS - KSPAIR)
                  ENDIF
*       Predict current KS variables and save at end of routine RESOLV.
                  CALL RESOLV(KSPAIR,3)
                  IPHASE = 2
                  JCOMP = 0
                  CALL KSTERM
                  KS = 1
               ELSE IF (DM.NE.0.D0) THEN
*       Implement mass loss and expand KS orbit at constant eccentricity.
                  CALL HCORR(I,DM,RNEW)
                  ITRY = 1
               END IF
            ELSE
*       Adopt KS treatment for single outer component or terminate merger.
               IF (I.EQ.2*KSPAIR.AND.NAME(I).LE.NZERO.AND.
     &             NAME(N+KSPAIR).GT.-2*NZERO.AND.DM.GT.0.0D0.AND.
*    &             H(KSPAIR) + DM/SEMI.LT.-ECLOSE.AND.KW.LT.13) THEN
     &             H(KSPAIR) + DM/SEMI.LT.-ECLOSE.AND.
     &             (KW.LT.10.OR.(KW.LT.11.AND.KZ(25).LT.2).OR.
     &             (KW.LT.13.AND.KZ(25).LT.1))) THEN
                  CALL HCORR(I,DM,RNEW)
               ELSEIF(DM.GT.0.D0)THEN
                  IPHASE = 7
                  CALL RESET
                  GO TO 1
               END IF
            END IF
         END IF
*
*       Check for end of blue straggler evolution (allow 5% extra).
         IF (M1.GT.1.05*TURN.AND.KSTAR(I).EQ.1.AND.KW.EQ.2) THEN
            if(rank.eq.0)
     &      WRITE(6,920) I, NAME(I), KW, TPHYS, AGE, M1, M0
 920        FORMAT(' END BS:    I NAM KW TP AGE M1 M0 ',
     &                          2I6,I4,2F8.1,2F6.1)
         END IF
*
*       Perform neighbour force corrections on significant mass loss (no WD).
         IF (ICORR) THEN
*
*       Include optional diagnostics for mass loss orbit (filename MDOT).
*           IF (KZ(21).GT.2) THEN
*              CALL MTRACE(I,DM)
*           END IF
*
*       Update the mass (single body, standard KS or merger).
            IF(IGHOST.EQ.0)THEN
               BODY(I) = M1/ZMBAR
            ELSE
               J = 2*KSPAIR - 2 + IHDOT
               BODY(J) = BODY(J) - DM
            ENDIF
            BODY0(I) = M0/ZMBAR
            IF(DMSUN.LT.TOL) GOTO 250
*
*       Accumulate total mass loss (solar units) and reduce cluster mass.
            ZMDOT = ZMDOT + DMSUN
            ZMASS = ZMASS - DM
*       Update the maximum single body mass but skip compact subsystems.
            IF(MASS(K)/ZMBAR.GE.0.99*BODY1.AND.NSUB.EQ.0)THEN
               BODY1 = 0.d0
               DO 35 J = 1,N
                  BODY1 = MAX(BODY1,BODY(J))
  35           CONTINUE
            ENDIF
*
*       Update the mass loss counters for types > 2.
            IF(KW.EQ.3)THEN
               ZMRG = ZMRG + DMSUN
            ELSEIF(KW.EQ.4)THEN
               ZMHE = ZMHE + DMSUN
            ELSEIF(KW.EQ.5.OR.KW.EQ.6)THEN
               ZMRS = ZMRS + DMSUN
            ELSEIF(KW.GE.7.AND.KW.LE.9)THEN
               ZMNH = ZMNH + DMSUN
            ELSEIF(KW.GE.10.AND.KW.LE.12)THEN
               ZMWD = ZMWD + DMSUN
            ELSEIF(KW.EQ.13.OR.KW.EQ.15)THEN
               ZMSN = ZMSN + DMSUN
            ELSEIF(KW.EQ.14)THEN
               ZMBH = ZMBH + DMSUN
            ENDIF
*
*       Check optional diagnostics.
            IF (rank.eq.0.and.KZ(19).GT.3) THEN
               WRITE(6,36) I, NAME(I), KW, KSTAR(I), BODY(I)*ZMBAR,
     &                     DMSUN, ZMDOT, TPHYS
  36           FORMAT(' MDOT:    I NM KW K* MS DMS ZMDOT T6 ',
     &                           4I5,F6.1,F7.2,F7.1,F8.1)
            END IF
*
*       Replace any KS components by corresponding c.m. for main procedures.
            IF (I.LT.IFIRST) THEN
               IKS = I
               I = N + KSPAIR
               I1 = 2*KSPAIR - 1
*       Predict coordinates & velocities of any unperturbed KS components.
               IF (LIST(1,I1).EQ.0) THEN
                  CALL RESOLV(KSPAIR,1)
               END IF
            END IF
*
*       Switch to the associated c.m. during force correction (save #I).
            IF (IGHOST.GT.0) THEN
               II = I
               I = N + KSPAIR
            END IF
*
            NNB = LIST(1,I)
            DO 38 L = 1,NNB+1
                ILIST(L) = LIST(L,I)
   38       CONTINUE
*       Ensure at least one neighbour.
            IF (NNB.EQ.0) THEN
                ILIST(2) = N
                IF (I.EQ.N) ILIST(2) = N - 1
                LIST(2,I) = ILIST(2)
                LIST(1,I) = 1
                NNB = 1
            END IF
*       Include body #I at the end (counting from location #2).
            NNB2 = NNB + 2
            ILIST(NNB2) = I
*
*       Check prediction of neighbours and body #I to current time.
            DO 40 L = 2,NNB2
               J = ILIST(L)
               IF (T0(J).LT.TIME) THEN
                  CALL XVPRED(J,-2)
               END IF
   40       CONTINUE
*
*       Define logical variable to control optional WD kicks (avoids NaN).
            IKICK = .FALSE.
            IF (KZ(25).EQ.1.AND.(KW.EQ.10.OR.KW.EQ.11)) IKICK = .TRUE.
            IF (KZ(25).EQ.2.AND.KW.EQ.12) IKICK = .TRUE.
*
       IF(rank.eq.0.and.KW.GE.10)PRINT*,
     *   'MDOT2:KW,N,IKICK,M0,MC,AGE,TN,DM,DMSUN=',
     *   KW,NAME(I),IKICK,M0,MC,AGE,TN,DM,DMSUN

*       Perform total force & energy corrections (delay dF if DMSUN > 0.1).
            IF (DMSUN.LT.0.05.AND.(KW.LT.10.OR..NOT.IKICK)) THEN
                CALL FICORR(I,DM)
            ELSE
                CALL FCORR(I,DM,KW)
            END IF
*
*       Set flag to ensure new sorting after CALL FPOLY (need IPHASE < 0).
            IF (I.GT.N) IPOLY = -1
*
*       Initialize new polynomials of neighbours & #I for DMSUN > 0.1.
            IF (ABS(DMSUN).GT.0.1.OR.KW.GE.13.OR.IKICK) THEN
*
*       Obtain new F & FDOT and time-steps (no gain skipping WDs).
               DO 50 L = 2,NNB2
                  J = ILIST(L)
                  IF (L.EQ.NNB2) J = I
*                 CALL DTCHCK(TIME,STEP(J),DTK(MAXBLK)) ! no effect (08/10).
                  DO 45 KK = 1,3
                     X0DOT(KK,J) = XDOT(KK,J)
   45             CONTINUE
                  CALL FPOLY1(J,J,0)
                  CALL FPOLY2(J,J,0)
                  IPOLY = -1
   50          CONTINUE
*       Check optional c.m. correction.
               IF (KZ(14).GE.3.AND.KZ(31).GT.0) THEN
                  CALL XVPRED(IFIRST,NTOT)
                  CALL CMCORR
               END IF
            END IF
            TPREV = TIME - STEPX
*       Set indicator < 0 for new sorting.
            IF (IGHOST.GT.0) I = II
            IF(IPHASE.EQ.-3) IPHASE = 0
         END IF
*
*       Ensure that massless supernova remnant will escape next output.
         IF (KW.EQ.15) THEN
            T0(I) = TADJ + DTADJ
            STEP(I) = 1.0D+06
            STEPR(I) = 1.0D+06
            RI = SQRT(X(1,I)**2 + X(2,I)**2 + X(3,I)**2)
            VI = SQRT(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2)
            DO 55 L = 1,3
*       Ensure that ghost will escape next output (far from fast escapers).
               X0(L,I) = MIN(1.0d+04+X(L,I),1000.0*RSCALE*X(L,I)/RI)
               X(L,I) = X0(L,I)
               X0DOT(L,I) = SQRT(0.004*ZMASS/RSCALE)*XDOT(L,I)/VI
               XDOT(L,I) = X0DOT(L,I)
               F(L,I) = 0.0D0
               FDOT(L,I) = 0.0D0
               D0(L,I) = 0.0
               D1(L,I) = 0.0
               D2(L,I) = 0.0D0
               D3(L,I) = 0.0D0
               D1R(L,I) = 0.0
               D2R(L,I) = 0.0D0
               D3R(L,I) = 0.0D0
   55       CONTINUE
         END IF
*
 250     CONTINUE
*
*       Restore index in case of KS component (used as c.m. above).
         IF (IKS.GT.0) THEN
            I = IKS
            IKS = 0
*       Re-initialize KS polynomials for perturbed motion.
            IF (LIST(1,I1).GT.0) THEN
               CALL RESOLV(KSPAIR,1)
               CALL KSPOLY(KSPAIR,1)
            END IF
         END IF
*
*       Update event counters for types > 2.
         IF(KW.NE.KW0)THEN
            IF(KW.EQ.3)THEN
               NRG = NRG + 1
            ELSEIF(KW.EQ.4)THEN
               NHE = NHE + 1
            ELSEIF(KW.EQ.5)THEN
               NRS = NRS + 1
            ELSEIF(KW.GE.7.AND.KW.LE.9.AND.KW0.LE.6)THEN
               NNH = NNH + 1
            ELSEIF(KW.GE.10.AND.KW.LE.12.AND.KW0.LE.9)THEN
               NWD = NWD + 1
            ELSEIF((KW.EQ.13.OR.KW.EQ.15).AND.KW0.LE.12)THEN
               NSN = NSN + 1
            ELSEIF(KW.EQ.14.AND.KW0.LE.13)THEN
               NBH = NBH + 1
            ENDIF
         ENDIF
*
*       Include consistency warnings.
         IF (RNEW - RADIUS(I).GT.0.5*RADIUS(I)) THEN
            if(rank.eq.0)
     &      WRITE(43,924) I, NAME(I), TPHYS, DT/1.0d+06,
     &                    KSTAR(I), KW, M0, M1, RADIUS(I)*SU, RNEW*SU
 924        FORMAT(' EXPAND!    I NM TP DTP K* KW M0 M R RN ',
     &                          2I6,F7.1,F7.3,2I4,2F7.1,2F7.1)
            CALL FLUSH(43)
         END IF
*
*       Update R, L, classification type & spin angular momentum.
         RADIUS(I) = RNEW
         ZLMSTY(I) = LNEW
         KSTAR(I) = KW
         JSPIN(K) = MAX(JSPIN(K) - DJSPIN(K),1.0D-10)
*
*       Ensure that the star does not spin up beyond break-up.
         OSPBRU = TWOPI*SQRT(MASS(K)*AURSUN**3/RAD(K)**3)
         JSPBRU = (K2STR(K)*(MASS(K)-MASSC(K))*RAD(K)*RAD(K) +
     &             K3*MASSC(K)*RADC(K)*RADC(K))*OSPBRU
         IF(JSPIN(K).GT.JSPBRU) JSPIN(K) = JSPBRU
         SPIN(I) = JSPIN(K)/SPNFAC
*
*       Update epoch and check binary diagnostics for transition to new type.
         IF (KW.NE.KW0) THEN
*           EPOCH(I) = AGE0(K) + EPCH0(K) - AGE
*           IF(KW.GT.6)THEN
            IF(rank.eq.0.and.KW.LT.0)THEN
               WRITE(6,925)I,NAME(I),KW0,KW,M10,M1,DMS(K),RM,AGE
 925           FORMAT(' MDOT CHANGE: I NM K* M0 M1 DM R AGE',
     &                               2I6,2I4,2F6.1,F7.3,F7.2,F8.1)
            ENDIF
            IF (I.LT.IFIRST.AND.KZ(8).GT.3) THEN
               CALL BINEV(KSPAIR)
            END IF
         END IF
*
*       Include optional diagnostics (new type or significant mass loss).
         IF (KZ(19).GT.3.AND.(KW0.NE.KW.OR.ICORR)) THEN
            IF (KW0.NE.KW) THEN
               WHICH1 = ' TYPE   '
            ELSE
               WHICH1 = ' MASS   '
            END IF
            if(rank.eq.0)
     &      WRITE(6,926)WHICH1, TPHYS, I, NAME(I), DMR, KW0, KW, M0, M1,
     &                  RADIUS(I)*SU, EMDOT
 926        FORMAT(' NEW',A8,' TPHYS I NAM DM/M KW0 KW M0 M R EMD ',
     &                         F7.1,2I5,F6.2,2I3,2F6.1,F7.1,F10.5)
         END IF
*
*       Base new time scale for changes in radius & mass on stellar type.
         EPOCH(I) = AGE0(K) + EPCH0(K) - AGE
*        TEV(I) = (AGE0(K) + EPCH0(K))/TSTAR
         TEV(I) = TEVK
         TEV0(I) = TEV(I)
*        IF(I.EQ.IGHOST.OR.BODY(I).LE.0.0) RM0 = M1
         CALL TRDOT(I,DTM,M1)
         TEV(I) = TEV(I) + DTM
         IF(IHDOT.EQ.2)THEN
            TEV(I0) = TEV(I)
            IF(NAME(IGHOST).LT.0.OR.NAME(IGHOST).GT.NZERO)THEN
                TEV(IGHOST) = TEV(I)
            ENDIF
         ENDIF
*
*       See if former KS pair can be regularized again.
         IF (KS.GT.0) THEN
            ICOMP = IFIRST
            JCOMP = IFIRST + 1
            RIJ2 = (X(1,ICOMP) - X(1,JCOMP))**2 +
     &             (X(2,ICOMP) - X(2,JCOMP))**2 +
     &             (X(3,ICOMP) - X(3,JCOMP))**2
            IF (RIJ2.LT.RMIN22) THEN
*       Set temporary IPHASE > 0 for KSREG.
               IPHASE = 1
               CALL KSREG
               KSPAIR = NPAIRS
*       Restore current time to prevent small subsequent steps.
               TIME = TBLOCK
               IF ((KW.EQ.13.OR.KW.EQ.14).AND.H(NPAIRS).LT.0.0) THEN
*              IF (KW.GE.10.AND.KW.LE.14.AND.H(NPAIRS).LT.0.0) THEN
                  J = NTOT
                  SEMI = -0.5d0*BODY(J)/H(NPAIRS)
                  RA = R(NPAIRS)/SEMI
                  ECC2 = (1.0 - RA)**2 + TDOT2(NPAIRS)**2/(BODY(J)*SEMI)
                  TK = DAYS*SEMI*SQRT(SEMI/BODY(J))
                  if(rank.eq.0)
     &            WRITE(6,928) KW, SQRT(ECC2), RA, SEMI*SU, TK, STEP(J),
     &                         BODY(J)*ZMBAR, (XDOT(KK,J)*VSTAR,KK=1,3)
 928              FORMAT(' WD/NS BINARY    KW E R/A A P DT M V ',
     &                                  I4,2F6.2,1P,3E10.2,0P,4F7.1)
               END IF
               IF(KW.GE.10.AND.H(NPAIRS).LT.0.0)THEN
                  IF(KZ(8).GT.3)THEN
                     CALL DEGEN(NPAIRS,NPAIRS,3)
                  ELSE
                     J1 = 2*NPAIRS - 1
                     J2 = J1 + 1
                     IF(KSTAR(J1).GE.10.AND.KSTAR(J2).GE.10)THEN
                        NDD = NDD + 1
                     ENDIF
                  ENDIF
               ENDIF
            ELSE
               KSPAIR = 0
            END IF
            KS = 0
*       Ensure IPHASE < 0 on RETURN for new sorting after KSTERM.
            IPOLY = -1
         END IF
*
*       Note any formation of black holes or TZ object.
         IF(KW.EQ.14.AND.KSTAR(I).LT.14.AND.I.LE.N)THEN
            if(rank.eq.0)
     &      WRITE(6,930) I, NAME(I), KW0, KW, KSTAR(I), M0, M1, DMR
 930        FORMAT(' NEW BH/TZ    I NM K0 KW K* M0 M1 DM/M ',
     &                            2I6,3I4,3F7.2)
         END IF
*
*       Perform consistency check on massive WD (skip HMDOT treatment).
         IF(I.LE.N)THEN
            IF((KSTAR(I).GE.10.AND.KSTAR(I).LE.12).AND.IHDOT.EQ.0)THEN
               IF(rank.eq.0.and.BODY(I)*ZMBAR.GT.MCH)THEN
               WRITE(6,932)I,KW,BODY0(I)*ZMBAR,BODY(I)*ZMBAR,RADIUS(I)
 932           FORMAT(' DANGER!  MDOT   I K* M0 M R ',2I5,2F7.2,1P,E9.1)
               WRITE(6,934) NAME(I), TTOT, TEV0(I), TEV(I)
 934           FORMAT(' NAM T TEV0 TEV ',I6,3F10.2)
               STOP
               ENDIF
            ENDIF
         ENDIF
*
*       Include stellar radius check in case of NaN.
         IF (RADIUS(I).GE.0.0.AND.RADIUS(I).LT.1.0) GO TO 105
         if(rank.eq.0)
     &   WRITE(6,936) I, KSTAR(I), M1, RADIUS(I)
 936     FORMAT(' DANGER!    MDOT    I K* M1 R* ',I5,I4,F7.2,1P,E10.1)
         STOP
 105     CONTINUE
*
*       Check for chaotic binary.
         IF(ITRY.GT.0.AND.IGHOST.EQ.0)THEN
            II = I
            I = N + KSPAIR
            IPAIR = KSPAIR
            IF(KSTAR(I).LT.0)THEN
               IF(KW.GE.10.AND.DMR.GT.0.0)THEN
                  if(rank.eq.0)
     &            WRITE(6,940) NAME(II), KSTAR(II), KW, DMR
 940              FORMAT(' DEGEN SPIRAL    NAM K* DMR ',I5,2I4,1P,E9.1)
               ENDIF
               CALL CHRECT(IPAIR,DMR)
               IF (IPHASE.LE.0) IPOLY = -1
               IF(DMR.LT.0.0) GOTO 1
            ENDIF
            IF(KW.GE.10.AND.ICORR)THEN
               IF(KZ(8).GT.3)THEN
                  CALL DEGEN(IPAIR,IPAIR,3)
               ELSE
                  J1 = 2*IPAIR - 1
                  J2 = J1 + 1
                  IF(KSTAR(J1).GE.10.AND.KSTAR(J2).GE.10) NDD = NDD + 1
               ENDIF
            ENDIF
            ITRY = 0
         ENDIF
*
 220  CONTINUE
*
      IF(KSPAIR.EQ.0) KACC = 1
      IF(KACC.GT.1)THEN
         TEV(JX(1)) = MIN(TEV(JX(1)),TEV(JX(2)))
         IF(DTGR.LT.1.0D+10)THEN
            TEV(JX(1)) = MIN(TEV(JX(1)),TEV0(JX(1))+DTGR/TSTAR)
         ENDIF
         TEV(JX(2)) = TEV(JX(1))
      ENDIF
*
*       Check for Roche overflow condition (DTR < STEP) and implement any
*       orbit change owing to gravitational radiation, mass loss or tides.
      IF(KSPAIR.GT.0.AND.IGHOST.EQ.0)THEN
         I = N + KSPAIR
         IPAIR = KSPAIR
*        IF(H(IPAIR).LT.0.D0.AND.NAME(ICM).GT.0)THEN
         IF(H(IPAIR).LT.0.D0.AND.NAME(I).GT.0.AND.KSTAR(I).GE.0)THEN
            CALL BRAKE(IPAIR,DSEP,ECC)
            IF(IQCOLL.NE.0.OR.IPHASE.LT.0) GO TO 1
*
*       Include optional look-up time control for compact object binaries.
            IF (KSX.GE.13.AND.KZ(28).GT.2) THEN
              if(rank.eq.0)
     &         WRITE (6,944)  TTOT, NAME(2*IPAIR-1),KSTAR(2*IPAIR),
     &                        DTGR/TSTAR
  944          FORMAT (' GR CHECK   T NAM K* DTGR ',F8.2,I6,I4,1P,E9.1)
*              IF (KSX.GE.13.AND.KZ(28).GT.0) THEN
*              GE = (1.0 - ECC2)**3.5/(1.0 + 3.0*ECC2)
*              ZMX = MAX(BODY(2*IPAIR-1),BODY(2*IPAIR))
*              RATIO = MIN(BODY(2*IPAIR-1),BODY(2*IPAIR))/ZMX
*       Replace physical time-scale by N-body units (cf. Book, Eq.12.31).
*              TZ = GE*SEMI**4/(RATIO*(1.0 + RATIO)*ZMX**3)
*              TZ = 5.0/64.0*CLIGHT**5*TZ
*              IF (rank.eq.0.and.KZ(28).GT.2) THEN
*                 WRITE (6,944)  TTOT, NAME(2*IPAIR-1),KSTAR(2*IPAIR),TZ
*              END IF
*       Limit the time-scale in case of shrinkage from non-GR orbit.
*              TZ = MIN(TZ,100.0D0)
*       Specify new look-up as 1% of current GR time-scale (ignore old TEV).
*              TEV(2*IPAIR-1) = TIME + 0.01*TZ
*              TEV(2*IPAIR) = TEV(2*IPAIR-1)
*              TMDOT = MIN(TMDOT,TEV(2*IPAIR))
            END IF
*
            IF(KSTAR(I).GT.0.AND.KZ(34).GT.0)THEN
               JPAIR = -IPAIR
               CALL TRFLOW(JPAIR,DTR)
               IF(DTR.LT.STEP(I))THEN
                  CALL ROCHE(IPAIR)
                  IF(IQCOLL.NE.0.OR.IPHASE.LT.0) GO TO 1
               ELSE
                  TEV(I) = TIME + DTR
                  IF(KSTAR(I).GT.10.AND.MOD(KSTAR(I),2).EQ.1)THEN
                     J1 = 2*IPAIR - 1
                     J2 = J1 + 1
                     TEV(J1) = MAX(1.000002d0*TEV(I),TEV(J1))
                     TEV(J2) = MAX(1.000002d0*TEV(I),TEV(J2))
                  ENDIF
               ENDIF
               IF (KZ(34).EQ.1.AND.KSTAR(I).GT.0.AND.
     &            MOD(KSTAR(I),2).EQ.0)THEN
                  CALL SYNCH(IPAIR)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
*
*       Determine the time for next stellar evolution check.
   70 TMDOT = 1.0d+10
      DO 80 J = 1,NTOT
         IF(TEV(J).LE.TMDOT)THEN
            TMDOT = TEV(J)
         ENDIF
   80 CONTINUE
*
*       See whether any other stars should be considered.
      IF (TMDOT.LT.TIME) GO TO 1
*
*       Update any merged circularizing binary using TMDIS.
      IMERGE = 1
  104 IF (NMERGE.GT.0) THEN
          IF (KSTARM(IMERGE).LT.0.AND.TMDIS(IMERGE).LT.TIME) THEN
              KSPAIR = 0
*       Determine the KS index from c.m. identification.
              DO 106 I = N+1,NTOT
                  IF (NAME(I).EQ.NAMEM(IMERGE)) THEN
                      KSPAIR = I - N
                  END IF
  106         CONTINUE
              IF (KSPAIR.GT.0) THEN
                  IC = 0
                  DO 107 K = 1,NCHAOS
                      IF (NAMEC(K).EQ.NAME(N+KSPAIR)) IC = K
  107             CONTINUE
                  IC = MAX(IC,1)
                  SEMI = -0.5*(BM(1,IMERGE) + BM(2,IMERGE))/HM(IMERGE)
                  I1 = 2*KSPAIR - 1
                  if(rank.eq.0)
     &            WRITE (6,108) NAMEG(IMERGE), NAME(N+KSPAIR),KSTAR(I1),
     &                          ES(IC), TMDIS(IMERGE), SEMI
  108             FORMAT (' TMDIS TERM    NAMG NAMCM K* E TMD A ',
     &                                    2I7,I4,F8.4,F8.2,1P,E10.2)
*       Terminate hierarchy (RESET calls SPIRAL via CHRECT).
                  IPHASE = 7
                  CALL RESET
*       Consider the same merger index again.
                  IMERGE = IMERGE - 1
              END IF
          END IF
          IMERGE = IMERGE + 1
      END IF
      IF (IMERGE.LE.NMERGE) GO TO 104
*
*       Update the maximum single body mass but skip compact subsystems.
      IF(NSUB.EQ.0)THEN
         BODY1 = 0.d0
         DO 110 J = 1,N
            BODY1 = MAX(BODY1,BODY(J))
  110    CONTINUE
      ENDIF
*
      DO 122 J = N+1,NTOT
         IF(KSTAR(J).EQ.0.AND.NAME(J).GT.0.AND.TEV(J).LT.9.9E+09.AND.
     &      BODY(J).GT.0.0.AND.KZ(28).EQ.0)THEN
            if(rank.eq.0)
     &      WRITE(6,556)J,NAME(J),TEV(J)
            TEV(J) = 1.0d+10
 556        FORMAT(' MDOT TEV SMALL ',2I8,1P,E10.2)
         ENDIF
 122  CONTINUE
*
*       Ensure re-initialization of ICPERT (INTGRT) and KBLIST (SUBINT).
      IF (IPOLY.LT.0) THEN
          IPHASE = -1
          NBPREV = 0
      END IF
*
      RETURN
*
      END
***
