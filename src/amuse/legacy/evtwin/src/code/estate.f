! FXP: Fudged exponential function, used to avoid too large or too small numbers
! in the ionisation state of the Saha equation. Usesprecomputed values of
! the limiting exponentials.
      FUNCTION FXP(X)
      IMPLICIT NONE
      DOUBLE PRECISION :: FXP
      DOUBLE PRECISION, INTENT(IN) :: X
      DOUBLE PRECISION, PARAMETER :: FXP_LOW = -50.0D0
      DOUBLE PRECISION, PARAMETER :: FXP_HIGH = 50.0D0
      DOUBLE PRECISION, PARAMETER :: LOW = 1.928749847963917820563444D-22
      DOUBLE PRECISION, PARAMETER :: HIGH = 5.184705528587072045056000D+21
      IF (X>FXP_HIGH) THEN
         FXP = HIGH
         RETURN
      ELSEIF (X>FXP_LOW) THEN
         FXP = EXP(X)
         RETURN
      END IF
      FXP = LOW
      END FUNCTION



! STATEF: calculate equation of state variables
! AF is the logarithm of the variable f that parametrizes the el. degeneracy
! AT is the logarithm of the temperature
! Composition variables are passed using a common block.
      SUBROUTINE STATEF ( AF, AT )
      USE CONSTANTS
      USE SETTINGS
      USE OPACITY_CO
      IMPLICIT REAL*8 (A-H, L-Z)
c......../........./........./........./........./........./........./..
      COMMON /STAT1 / CSX(10), CS(90,127,10), CNU(60,41,2), W(4000), JZ
      COMMON /STAT2 / AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA, 
     :                SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, WR(22),
     &                WMU, DELTA, PHI, EXT, FKT, FKR, PRANDTL
      COMMON /STATFD/ RE, PE, QE, RET, PET, QET, REF, PEF, QEF, FDD(7)
      COMMON /STATFP/ RP, PP, QP, RPT, PPT, QPT, RPF, PPF, QPF, PFDD(7)
      COMMON /ABUND / XA(9), NA(9), NEO, NIO, NZZ, AVM, NE
      COMMON /ATDATA/ CH2, C1, C2, C3, CHI(26,9), COM(27), CAN(9), 
     :                CBN(9), KZN(9)
      DIMENSION HA(26), VA(26)
      DOUBLE PRECISION :: DKZN(9)
      PARAMETER (EPS = 1.0D-30)
      ! Inline (statement) functions.
      ! Fudged exponential function, use to prevent too large or too small numbers
      !FXP(VX) = DEXP(DMAX1(-50.0D0, DMIN1(50.0D0, VX)))

      ! We need the charges KZN in floating point operations a few times
      ! Do the int->float conversion once on startup to safe a few cycles
      DKZN(:) = KZN(:)
 
      ! Calculate F and G parameters, needed to calculate Fermi-Dirac integrals
      F = DEXP(AF)
      T = DEXP(AT)
      !T = MAX(1000.0, T)         ! Don't trust T < 1000K
      UF = F/(1.0D0 + F)
      WF = DSQRT(1.0D0 + F)
      PSI = AF + 2.0D0*(WF - DLOG(1.0D0 + WF))
      G = CTE*T*WF
* Evaluate Fermi-Dirac integrals according to Eggleton, Faulkner &
* Flannery (1973): RE, PE, SE and UE correspond to rho*, P*, S* and U*. 
* PSI is the usual electron degeneracy parameter
      CALL FDIRAC ( F, G )
      PE = G*PE
      PET = PET + 1.0D0
      PEF = PEF + 0.5D0*UF
      QE = QE/(RE*WF)
      SE = QE + 2.0D0*WF - PSI
      SEF = QE*(QEF - REF - 0.5D0*UF) - 1.0D0/WF
      SET = QE*(QET - RET)
      UE = SE + PSI - PE/(RE*CTE*T)
! Contributions from positrons
      RP = 0.0d0
      PP = 0.0d0
      SP = 0.0d0
      UP = 0.0d0
      RPT = 0.0d0
      PPT = 0.0d0
      SPT = 0.0d0
      RPF = 0.0d0
      PPF = 0.0d0
      SPF = 0.0d0
      XIP = -PSI - 2.0 / (CTE*T)
      IF (EOS_INCLUDE_PAIRPRODUCTION .AND. (XIP > -15.0)) THEN
!        Calculate effective degeneracy parameter for positrons
!        See also Cox & Giuli ¤24.9
!        NB: PF and PG are later returned as dP/dlogf and Pgas
         PAF = FIND_POSITRON_F(XIP)
         PF = DEXP(PAF)
         PUF = PF/(1.0D0 + PF)
         PWF = DSQRT(1.0D0 + PF)
         PG = CTE*T*PWF
         CALL POSITRON_FDIRAC(PF, PG)
!        Derivatives are now known with respect to the *positron* degeneracy
!        parameter pf, but we need to know them with respect to the *electron*
!        degeneracy parameter f. The conversion factor is easy, however,
!        FF = dlog pf / dlog f = - sqrt( (1+f)/(1+pf) )
!        (because dxi/dxip = -1)
         IF (RP > 0.0d0) then
         FF = -WF/PWF

         RPF = RPF * FF
         PPF = PPF * FF
         QPF = QPF * FF

!        Derivatives with respect to temperature get an extra term because
!        the effective positron degeneracy parameter depends on the temperature
!        as well as on the electron degeneracy parameter.
         RPT = RPT - 2.0D0 * RPF / ( CTE*T*G )
         PPT = PPT - 2.0D0 * PPF / ( CTE*T*G )
         QPT = QPT - 2.0D0 * QPF / ( CTE*T*G )

         PP = PG*PP
         PPT = PPT + 1.0D0  + PUF / PG
         PPF = PPF + 0.5D0 * PUF * FF
         QP = QP/(RP*PWF)
         SP = QP + 2.0D0*PWF - XIP
         SPF = QP*(QPF - RPF - 0.5D0*PUF) + 1.0D0/WF
         SPT = QP*(QPT - RPT)
         UP = SP + XIP - PP/(RP*CTE*T)
         end if
      end if
* Evaluate some quantities that do not depend on the state of ionization:
* the NA are the element number densities (per baryon); AVM is the average
* mass per baryon (in amu); NEO and NIO are the numbers of electrons and
* ions assuming complete ionization.
      NA(1:9) = XA(1:9)/CBN(1:9)
      NIO = SUM(NA(1:9))
      NEO = DOT_PRODUCT(DKZN(1:9), NA(1:9))
      NZZ = DOT_PRODUCT(DKZN(1:9)**2, NA(1:9))
      AVM = DOT_PRODUCT(CAN(1:9), NA(1:9))
* PRESSI gives a crude model for pressure ionization and 
* a model for Coulomb interactions, returning corrections to the electron
* chemical potential, pressure, entropy and internal energy.
* TI is 1eV/kT, DE is 1 amu * number of electrons/cm3
      TI = CEVB/T
      DE = RE*CD
      CALL PRESSI ( 1, TI, PSI, DE, REF, RET, F, DC, DVT, DVF, DPA, 
     &     DPAT, DPAF, DSA, DSAT, DSAF, DUA )
      DE = MAX(RE-RP, 1.0d-50)*CD
      DV = DC - PSI
      DVF = DVF - WF
* Contributions of the completely ionized species
      NE = 0.0D0
      NEF = 0.0D0
      NET = 0.0D0
      SI = 0.0D0
      SIF = 0.0D0
      SIT = 0.0D0
      UI = 0.0D0
      DO I = KION+1, 9
         NE = NE + DKZN(I)*NA(I)
         VM = CAN(I)*DSQRT(CAN(I))
         SI = SI - NA(I)*DLOG(NA(I)/VM + EPS)
      END DO
* Calculate ionization of the first KION elements.
      DVT_OVER_DVF = DVT / DVF
      DO I = KION, 1, -1
* compute potentials VA and number ratios HA of ionization state J
* relative to the ground state 1
         VA(1) = -CHI(1,I)*TI + DV
         HA(1) = FXP(VA(1))*COM(KZN(I))/COM(KZN(I)+1)
         SHA = 1.0D0 + HA(1)
         SJHA = HA(1)
         SCHA = CHI(1,I)*TI*HA(1)
         DO J = 2, KZN(I)
            VA(J) = -CHI(J,I)*TI + J*DV
            HA(J) = HA(J-1)*FXP(VA(J) - VA(J-1))*COM(KZN(I)+1-J)
     &           /COM(KZN(I)+2-J)
            SHA = SHA + HA(J)
            SJHA = SJHA + J*HA(J)
            SCHA = SCHA + CHI(J,I)*TI*HA(J)
         END DO
         SCHA_OVER_SHA = SCHA / SHA
         SJHA_OVER_SHA = SJHA / SHA
         VM = CAN(I)*DSQRT(CAN(I))
         SI = SI + NA(I)*DLOG(VM)
         IF ( I > 1 ) THEN
* contributions to electron number density NE, entropy SI and
* internal energy UI for Helium and heavier
            VX = NA(I)/SHA
            SI = SI - VX*DLOG(VX/COM(KZN(I)+1) + EPS)
            DO J = 1, KZN(I)
               NX = HA(J)*VX
               NXF = NX*DVF*(J - SJHA_OVER_SHA)
               NXT = NXF*DVT_OVER_DVF + NX*(CHI(J,I)*TI - SCHA_OVER_SHA)
               NE = NE + J*NX
               NEF = NEF + J*NXF
               NET = NET + J*NXT
               SIF = SIF - VA(J)*NXF
               SIT = SIT - VA(J)*NXT
               SI = SI - NX*DLOG(MAX(NX/COM(KZN(I) + 1 - J), EPS))
               UI = UI + CHI(J,I)*NX
            END DO
         END IF
      END DO
* Ionization and molecular dissciation of Hydrogen.
* partition function for H2 from Vardya (1960), Webbink (1975)
      DH2TI = CH2*TI
      EXPDH2TI = DEXP(-MIN(DH2TI, 100.0D0))
      D1TI = C1*TI
      D2TI = (C2*TI)**2
      D3TI = (C3*TI)**3
      ZET = 1.0D0 - (1.0D0 + DH2TI)*EXPDH2TI
      DZH2T = -DH2TI**2*EXPDH2TI/ZET
      DZH2TT = (DH2TI - 2.0D0 - DZH2T)*DZH2T
      ZH2 = 6608.8D0*ZET*DH2TI**(-2.5D0)*DEXP(-D1TI - D2TI - D3TI)
      ZH2 = MIN(ZH2, 1.0D100)       ! Avoid overflows later
      ZH2T = 2.5D0 + D1TI + 2.0D0*D2TI + 3.0D0*D3TI + DZH2T
      ZH2TT = -D1TI - 4.0D0*D2TI - 9.0D0*D3TI + DZH2TT
      ZH2S = ZH2*DSQRT(8.0D0)/VM
      H2A = CEN*(ZH2S/4.0D0)*DE/(T*DSQRT(T))/EXPDH2TI
      H2BT = DH2TI + 1.5D0 - ZH2T
      H2AT = RET - H2BT
* solve for densities of H+, H, and H2
      QA = 2.0D0*H2A + HA(1)*(1.0D0 + HA(1))
      QB = NE + HA(1)*(NE - NA(1))
      QC = NA(1)*NE
      HG = 2.0D0*QC/(DSQRT(QB*QB + 4.0D0*QA*QC) + QB)
      HI = HA(1)*HG
      NE = NE + HI                     ! Ionisation electrons / baryon
      NP = RP / dabs(RE - RP) * NE     ! Positrons / baryon
      EN = 1.0D0/NE
      H2 = H2A*HG*HG*EN
      NI = NIO - H2
* derivatives w.r.t. F and T
      QA = NE + 4.0D0*HG*H2A
      QB = HA(1)*(NE - 2.0D0*H2)
      QD = 1.0D0/(QA + QB)
      QC = 2.0D0*H2*QD
      QAF = (NEF - NE*REF )*QC
      QAT = (NET - NE*H2AT)*QC
      QF = HG*QD
      QBF = DVF*QF
      QBT = (CHI(1,1)*TI + DVT)*QF
      HGF = QAF - QB*QBF
      HGT = QAT - QB*QBT
      HIF = HA(1)*(QAF + QA*QBF)
      HIT = HA(1)*(QAT + QA*QBT)
      NEF = NEF + HIF
      NET = NET + HIT
      H2F = H2*REF  + EN*(2.0D0*H2A*HG*HGF - H2*NEF)
      H2T = H2*H2AT + EN*(2.0D0*H2A*HG*HGT - H2*NET)
* hydrogen contribution to entropy, internal energy
      SIF = SIF - VA(1)*HIF - H2BT*H2F
      SIT = SIT - VA(1)*HIT - H2BT*H2T + H2*(ZH2T + ZH2TT)
! Avoid overflow problems when HI, HG, H2 -> 0
      SI = SI - HI*DLOG(DABS(HI/COM(1)) + EPS)
     &        - HG*DLOG(DABS(HG/COM(2)) + EPS)
     &        - H2*(DLOG(DABS(H2/ZH2S) + EPS) - ZH2T)
      UI = UI + CHI(1,1)*HI + 0.5*CH2*(HI + HG)
* DB is 1 amu * number of baryons/cm3; RHO is the mass density in g/cm3
      DB = EN*DE
      DL = DLOG(DB)
      RHO = DB*AVM
      ARHO = DLOG(RHO)
      RT = (RE*RET + RP*RPT)/(RE + RP) - EN*NET
      RF = (RE*REF + RP*RPF)/(RE + RP) - EN*NEF
* second call to PRESSI compensates for spurious pressure and entropy terms
      DE = DB*NEO
      WMU = 1.0D0/(NEO+NI+H2)
      CALL PRESSI ( 0, TI, PSI, DE, RF, RT, F, DC, DVT, DVF, DPB, DPBT, 
     :            DPBF, DSB, DSBT, DSBF, DUB )
* pressure terms
      PE = CB*PE              ! Partial pressure due to electrons
      PP = CB*PP              ! Partial pressure due to positrons
      TCR = T*CR
      P0 = TCR*DB
      PI = NI*P0
      T2 = T*T
      T4 = T2*T2
      PR = CA*T4/3.0D0
      B = 4.0D0*PR/P0
! When positrons are present, ignore pressure ionisation - everything is
! fully ionised anyway.
      IF (PP > 0.0D0) THEN
         DPA = 0.0
         DPB = 0.0
         DPAF = 0.0
         DPBF = 0.0
         DPAT = 0.0
         DPBT = 0.0
         DSAF = 0.0
         DSBF = 0.0
         DSAT = 0.0
         DSBT = 0.0
         DSA = 0.0
         DSB = 0.0
      END IF
      PG = PE + PP + PI + TCR*(DPA - DPB)
      P = MAX(PG + PR, 1.0D-200)
      PF = (PE*PEF + PP*PPF      + PI*RF - H2F*P0 + TCR*(DPAF - DPBF))/P
      PT = (PE*PET + PP*PPT + PI + PI*RT - H2T*P0 + TCR*(DPAT - DPBT)+PR*4.D0)/P
      AP = DLOG(P)
* entropy, in erg/g/K
      DSF = NEF*DSA + NE*DSAF - NEO*DSBF - RF*B
      DST = NET*DSA + NE*DSAT - NEO*DSBT - (RT - 3.0D0)*B
      SF = CR*(-NI*RF           + NEF*SE + NE*SEF + NP*SPF + SIF + DSF)/AVM
      ST = CR*( NI*(1.5D0 - RT) + NET*SE + NE*SET + NP*SPT + SIT + DST)/AVM
      S = CR*(SE*NE + SP*NP + DSA*NE - DSB*NEO + B + (1.5D0*AT - DL + 2.5D0
     :        - DLOG(CEN))*NI + SI)/AVM
* internal energy, in erg/g
      U = TCR*(UE*NE + UP*NP + DUA*NE - DUB*NEO + 1.5D0*NI + ZH2T*H2 
     :       + 0.75D0*B + TI*UI)/AVM
* other thermodynamic quantities
      Q = MIN(PT*SF - PF*ST, 1.0d300)
      SCP = -Q/PF             ! Specific heat at constant pressure
      SCV = ST - SF*RT/RF     ! Specific heat at constant volume
      GRADA = SF/Q
      GAMMA = Q/(RT*SF-RF*ST)
      ZT = DSQRT(DABS((NE*REF/WF + NZZ)/NI))

! FIXME: the expression for phi is valid only for ideal gas+radiation!
      DELTA = RF*PT/PF - RT      ! Coefficient of thermal expansion (K&W (6.6))
      PHI   = 1.0                ! Coefficient of chemical expansion (K&W (6.6))
      
* TST ought to be zero, if all the above programming is correct
      TST = SF/CR - P*(RT*PF-RF*PT)/(TCR*RHO)
*** END OF THERMODYNAMIC CALCULATION. BEGINNING OF TABLE-BASED CALCULATIONS
      FRHO = ARHO/CLN
      TF = AT/CLN

! Cap values to allowed range for opacity tables
      FRHO = MIN(10.0D0, MAX(FRHO,-12.0D0))
      TF = MIN(9.3D0, MAX(TF,  3.0D0))
      FK = GET_OPACITY(FRHO, TF)
      XHI = 4.0D0*CL*PR/(FK*RHO*RHO*SCP*T)   ! cm**2/s
! Plasma and photon viscosity
! The plasma viscosity is from Spitzer (1962)
! Photon viscosity from Thomas (1930) and Weinberg (1971) with the
! correction factor (10/9) of van den Horn & van Weert (1981).
      !LLAMBDA = LOG(3.0D0/(2.0D0*ECHAR**3) * SQRT((BOLTZM*T)**3/(CPI*DB*ZT**5)))
      LLAMBDA = LOG(3.0D0/(2.0D0*ECHAR**3*ZT**2) * SQRT((BOLTZM*T)**3/(CPI*NE)))
      MU_PLASMA = 0.406*SQRT( AVM/(NI+NEO)*(BOLTZM*T)**5 )/( (ZT*ECHAR)**4*LLAMBDA)
      MU_RAD = 8.0D0*PR/(3.0D0*CL*FK*RHO)    ! gram/cm s
      NU = ( MU_PLASMA + MU_RAD ) / RHO
! Prandtl number
      PRANDTL = NU/XHI
! Derivatives of opacity: (dlogk/dlogT)_rho and (dlogk/dlogrho)_T
! FIXME: these are the results for Kramer's opacity, should take
! derivative of spline instead to get value from tables.
      FKR =  1.0d0
      FKT = -3.5d0
* Neutrino loss rates from Itoh et al (1983-1992)
      EN = 0.0D0
      IF ( TF >= 7.0D0 .AND. FRHO >= 0.0D0 ) THEN
         TT = 20.0D0*(TF - 6.95D0)
         IT = MAX0(1, MIN0(59, INT(TT)))
         TT = TT - IT
         TU = 1.0D0 - TT
         RN = 4.0D0*(FRHO + 0.25D0)
         IR = MAX0(1, MIN0(40, INT(RN)))
         RN = RN - IR
         RU = 1.0D0 - RN
         ENP = TT*(RN*CNU(IT + 1, IR + 1, 1) + RU*CNU(IT + 1, IR, 1))
     &           + TU*(RN*CNU(IT, IR + 1, 1) + RU*CNU(IT, IR, 1))
         ENB = TT*(RN*CNU(IT + 1, IR + 1, 2) + RU*CNU(IT + 1, IR, 2))
     &           + TU*(RN*CNU(IT, IR + 1, 2) + RU*CNU(IT, IR, 2))
         EN = -10.0D0**ENP - NZZ*10.0D0**ENB
      END IF
! Neutrino rate, from Itoh et al. (1996)
!      CALL GET_NEUTRINO_RATE(T, RHO, NIO, NEO, EN)
      RETURN
      END

! Calculate the Fermi-Dirac integrals, as parameterized by functions f and g
! Returns its results via a common block.
      SUBROUTINE FDIRAC ( F, G )
      IMPLICIT REAL*8 (A-H, L-Z)
      DIMENSION FF(4), GG(4), C(4,4,3), VW(4,4), VX(4,4)
      DATA C/ 2.315472D0,  7.128660D0,  7.504998D0,  2.665350D0, 
     :        7.837752D0, 23.507934D0, 23.311317D0,  7.987465D0, 
     :        9.215560D0, 26.834068D0, 25.082745D0,  8.020509D0, 
     :        3.693280D0, 10.333176D0,  9.168960D0,  2.668248D0, 
     :        2.315472D0,  6.748104D0,  6.564912D0,  2.132280D0, 
     :        7.837752D0, 21.439740D0, 19.080088D0,  5.478100D0, 
     :        9.215560D0, 23.551504D0, 19.015888D0,  4.679944D0, 
     :        3.693280D0,  8.859868D0,  6.500712D0,  1.334124D0, 
     :        1.157736D0,  3.770676D0,  4.015224D0,  1.402284D0, 
     :        8.283420D0, 26.184486D0, 28.211372D0, 10.310306D0, 
     :       14.755480D0, 45.031658D0, 46.909420D0, 16.633242D0, 
     :        7.386560D0, 22.159680D0, 22.438048D0,  7.664928D0/
      COMMON /STATFD/ D(3,3), XTT, XFT, XFF, XTTT, XFTT, XFFT, XFFF 
c Evaluate Fermi-Dirac integrals (Eggleton, Faulkner & Flannery 1973).
c Matrix D contains rho*, P* and Q* and 1st logarithmic derivatives w.r.t. 
c T, f. XTT etc are 2nd and 3rd log derivs of rho*
      VF = 1.0D0/(1.0D0 + F)
      VG = 1.0D0/(1.0D0 + G)
      UF = F*VF
      UG = G*VG
      FDF = G*G + G
      FDF = UF*FDF*DSQRT(FDF)
      FF(1) = 1.0D0
      GG(1) = 1.0D0
      DO I = 2, 4
         FF(I) = F*FF(I - 1)
         GG(I) = G*GG(I - 1)
         FDF = FDF*VF*VG
      END DO
      IND = 4
      DO I = 1, 3
         VX(1:IND, 1:IND) = 0.0d0
         DO IJ = 1, 4
            DO IK = 1, 4
               VW(1, 1) = C(IK, IJ, I)*GG(IJ)*FF(IK)
               DO IL = 1, IND - 1
                  VW(IL + 1, 1) = (IJ - 1)*VW(IL, 1)
                  DO  IM = 1, IND - IL
                     VW(IL, IM + 1) = (IK - 1)*VW(IL, IM)
                  end do
               end do
               VX(1:IND, 1:IND) = VX(1:IND, 1:IND) + VW(1:IND, 1:IND)
            end do
         end do
         
         WV = 1.0D0/VX(1, 1)
         VX(1, 2:IND) = VX(1, 2:IND)*WV
         VX(2:IND, 1:IND-1) = VX(2:IND, 1:IND-1)*WV
         D(I, 1) = FDF*VX(1, 1)
         D(I, 2) = VX(2, 1) + 1.5D0 - 1.5D0*UG
         WW = 0.5D0*D(I, 2) - 4.0D0
         D(I, 3) = VX(1, 2) + 1.0D0 + WW*UF
         IF ( I == 1 ) THEN
* second and third order derivatives of density, needed in PRESSI
            XTT = VX(3, 1) - VX(2, 1)**2 - 1.5D0*UG*VG
            XFT = VX(2, 2) - VX(2, 1)*VX(1, 2) + 0.5D0*UF*XTT
            XFF = VX(1, 3) - VX(1, 2)**2 + UF*(XFT + WW*VF - 0.25D0*XTT*UF)
            XTTT = VX(4, 1) + VX(2, 1)*(2.0D0*VX(2, 1)**2 - 3.0D0*VX(3, 1))
     &              - 1.5D0*(1.0-G)*UG*VG**2
            XFTT = VX(3, 2) - 2.0D0*VX(2, 1)*VX(2, 2) - VX(1, 2)*(VX(3, 1) 
     &             - 2.0D0*VX(2, 1)**2) + 0.5D0*UF*XTTT
            XFFT = VX(2, 3) - 2.0D0*VX(2, 2)*VX(1, 2) - VX(2, 1)*(VX(1, 3) 
     &            - 2.0D0*VX(1,2)**2) + UF*(XFTT + 0.5D0*VF*XTT - 0.25D0*UF*XTTT)
            XFFF = VX(1, 4) + VX(1, 2)*(2.0D0*VX(1, 2)**2 - 3.0D0*VX(1, 3)) 
     &         + UF*(1.5D0*(XFFT + VF*XFT) - UF*(0.75D0*(XFTT + VF*XTT) - 
     &         0.125D0*UF*XTTT) + WW*(1.0D0 - F)*VF**2)
         END IF
         IND = 2
      END DO
      RETURN
      END

! --------------------------------------------------------------------------
! PRESSI
! --------------------------------------------------------------------------
! Non-ideal corrections to the equation-of-state: pressure ionisation and
! Coulomb interactions.
! As explained in Pols&al. (1995), the effect of pressure ionisation is
! calculated twice: once for the actual number of electrons Ne and once for
! the total number of electrons Ne0, to make sure the corrections to
! pressure and entropy vanish in the limit of complete ionisation.
!
! Input parameters:
!  IPR      - 1 or 0, if 1 also calculates Coloumb interactions
!  TI       - 1eV/kT, recipocal temperature
!  PSI      - electron degeneracy parameter
!  XI       - electron number density in funny units, 1 amu * electrons/cm3
!  XF       - logarithmic derivative of electron density to f, dln rhoe/dln f
!  XT       - logarithmic derivative of electron density to T, dln rhoe/dln T
!  F        - Eggleton F parameter (alternative degeneracy parameter)
! Output parameters:
!  DC       - Correction for electron chemical potential (sign flipped)
!  DCT      - Derivative of above with respect to log T
!  DCF      - Derivative of above with respect to log f 
!  DP       - Correction for the pressure
!  DPT      - Derivative of above with respect to log T
!  DPF      - Derivative of above with respect to log f
!  DS       - Correction for the entropy
!  DST      - Derivative of above with respect to log T
!  DSF      - Derivative of above with respect to log f
!  DU       - Correction for the internal energy
! The output parameters are calculated by taking derivatives of the free
! energy, which requires high order derivatives of the electron density.
! These are calculated by FDIRAC and passed in the STATFD common block.
! --------------------------------------------------------------------------
      SUBROUTINE PRESSI ( IPR, TI, PSI, XI, XF, XT, F, 
     :                     DC, DCT, DCF, DP, DPT, DPF, DS, DST, DSF, DU )
      USE CONSTANTS;
* Computes effect (`pressure ionisation') on EOS of a free energy contribution 
* delta(F) = -R.T.Ne.G(X,Y), where X = XI = Ne/V = ne, and Y = YI = XIH/(R.T).
* Also does Coulomb interaction.
* delta(el. chem. pot.) = -DC = 1/(R.T) dF/dNe
* delta(P) = R.T.DP = -dF/dV 
* delta(S) = R.Ne.DS = -dF/dT.
* Works firstly with indep. vbles X, Y (called XI, YI), which are
* effectively Ne, V and T, but then has to get derivatives of DC, DP, DS
* with respect to independent variables f and T.
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON /STATFD/ FD(9), XTT, XFT, XFF, XTTT, XFTT, XFFT, XFFF
      COMMON /ABUND / XA(9), NA(9), NEO, NIO, NZZ, AVM, NE1
      DOUBLE PRECISION, PARAMETER :: CA1 = 0.89752D0
      DOUBLE PRECISION, PARAMETER :: CA2 = 0.768D0
      DOUBLE PRECISION, PARAMETER :: CA3 = 0.208D0
      DOUBLE PRECISION, PARAMETER :: CP1 = 3.0D0
      DOUBLE PRECISION, PARAMETER :: CP2 = 0.35D0
      DOUBLE PRECISION, PARAMETER :: CP3 = 2.0D0
      DOUBLE PRECISION, PARAMETER :: CP4 = 3.0D-2
      CBRT(VX) = DABS(VX)**C3RD
* Pressure ionization
      YI = 13.595D0*TI
      EE = 1.0D0/(1.0D0 + XI/CP4)
      WX = MIN((CP1/XI)**CP2, 300.0D0)
      CC = DEXP(-WX)
      BB = YI + PSI - CP3*DLOG(EE)
c extra free energy is -R.T.Ne.GI(X,Y)
c GI = exp{-WX(X)}{Y + psi(X,Y) + 2ln(1+X/const)}
c It OUGHT to have 3psi/5 for psi, but works worse!
      GI = CC*BB
      FF1 = 1.0D0/(F + 1.0D0)
c AA = dlnf/dpsi; THeta = dlnne/dpsi; also needed in Coulomb corrections
c f, T derivs, then X, Y.
      AA = DSQRT(FF1)
      TH = AA*XF
      FF2 = -0.5D0*F*FF1
      AF = AA*FF2
      THF = AF*XF + AA*XFF
      THT = AA*XFT
      RXF = 1.0D0/XF
      THX = THF*RXF
      THY = THX*XT - THT
c first and second derivatives dPSI/dlnX ... d2PSI/dlnY2
      PSX = 1.0D0/TH
      PSY = XT*PSX
      W2 = PSX*PSX
      PSXX = -THX*W2
      PSXY = -THY*W2
      PSYY = PSXY*XT + (XFT*XT*RXF - XTT)*PSX
c derivative -dlnEE/dlnX; -d2lnEE/dlnX2 = -EE*dlnEE/dlnX
      EEX = CP3*(XI/CP4)*EE
c derivatives of BB
      BX = PSX + EEX
      BY = YI + PSY
      BXX = PSXX + EE*EEX
      BXY = PSXY
      BYY = YI + PSYY
c derivatives of CC
      CX = CP2*WX*CC
      CXX = CP2*CX*(WX - 1.0D0)
C derivatives of GI
      DGDX = CC*BX + CX*BB
      DGDY = CC*BY
      DGDXX = CC*BXX + 2.0D0*CX*BX + CXX*BB
      DGDXY = CC*BXY + CX*BY
      DGDYY = CC*BYY
      IF ( IPR == 1 ) THEN
* Coulomb interaction.
c further derivatives of AA, THeta
         AFF = 3.0D0*AF*FF2 + AF
         THFF = AFF*XF + 2.0D0*AF*XFF + AA*XFFF
         THFT = AF*XFT + AA*XFFT 
         THTT = AA*XFTT 
c d2TH/dlnX2, etc
         TOF = XT*RXF
         WXX = THFF - THX*XFF
         WXY = THX*XFT - THFT
         THXX = WXX*RXF*RXF
         THXY = WXY*RXF + THXX*XT
         THYY = TOF*(TOF*WXX + 2.0D0*WXY) + THTT - THX*XTT
* GAM is the plasma interaction parameter. Note that THC = ZT**2*NIO/NEO
         THC = TH + DABS(NZZ/NEO)
         WW = 1.0D0/CA2
         GAM = CBRT(XI*(CPL*NEO/NIO)**2*C3RD)*TI/CEVB*THC
c new BB and EE, and their derivates
         BB = (CA1*DSQRT(3.0D0/GAM))**WW
         EE = (GAM/(GAM + CA3))**WW
         RBE = 1.0D0/(EE + BB)
c further addition GE to free en; adds to previous GI.
c GE = GE(GAMma), GAM = const. * X**0.33 * Y * (THeta(X,Y) + const)
         GE = (NIO/NEO)*CA1*GAM*RBE**CA2
         EEG = CA3/(GAM + CA3)
         BEG = (EEG*EE - 0.5D0*BB)*RBE
         BEGG = (EEG*EEG*(1.0D0 - GAM*CA2/CA3)*EE + 0.25D0*BB)*RBE
         GEG = 1.0D0 - BEG
         DGDG = GE*GEG
         DGDGG = GE*(GEG*GEG + (BEG*BEG - BEGG)*WW)
         RTHC = 1.0D0/THC
         WX = THX*RTHC
         WY = THY*RTHC
c dlnGAM/dlnX, etc
         GAMX = C3RD + WX
         GAMY = 1.0D0 + WY
         GAMXX = THXX*RTHC - WX*WX
         GAMXY = THXY*RTHC - WX*WY
         GAMYY = THYY*RTHC - WY*WY
         GI = GI + GE
c derivatives w.r.t. X, Y; in effect w.r.t. Ne, V, and T, since X = Ne/V
         DGDX = DGDX + DGDG*GAMX
         DGDY = DGDY + DGDG*GAMY
         DGDXX = DGDXX + DGDGG*GAMX**2 + DGDG*GAMXX
         DGDXY = DGDXY + DGDGG*GAMX*GAMY + DGDG*GAMXY
         DGDYY = DGDYY + DGDGG*GAMY**2 + DGDG*GAMYY
      END IF
* evaluate changes to el. chem. potential (-DC), pressure (R.T.DP), and
* entropy (R.Ne.DS), and their derivatives w.r.t. log(f) and log(T)
      DC = DGDX + GI
      WW = DGDXX + DGDX
      WT = WW*XT
      DCF = WW*XF
      DCT = WT - DGDXY - DGDY
      DP = -XI*DGDX
      DPF = -XI*DCF
      DPT = XI*(DGDXY - WT) + DP
      DS = GI - DGDY
      DST = DGDX - DGDXY
      DSF = DST*XF
      DST = DST*XT - DGDY + DGDYY
      DU = -DGDY
      RETURN
      END

      SUBROUTINE NUCRAT ( TL )
      USE CONSTANTS;
* Compute rates of (at present) 20 nuclear rections, and the corresponding
* energy and neutrino release
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON /STAT1 / CAT(119230), CRT(200,20), KCSX
      COMMON /STAT2 / W1(4), RHO, W2(4), ZT, W3(10), RRT(1), RPP, R33,  
     :     R34, RBE, RBP, RPC, RPN, RPO, R3A, RAC, RAN, RAO, RANE, 
     :     RCC, RCO, ROO, RGNE, RGMG, RCCG, RPNG, EX, ENX, WMU, DELTA,
     :     PHI, EXT, FKT, FKR, PRANDTL
      COMMON /ABUND / XA(9), N1, N4, N12, N14, N16, N20, N24, N28, N56,
     &                NE, NI, NZZ, AVM, NE1
c......../........./........./........./........./........./........./..
      COMMON /NCDATA/ QRT(20), QNT(20), CZA(92), CZB(92), CZC(92),
     &                CZD(92), VZ(9)
      DOUBLE PRECISION :: RTT(20)
      DOUBLE PRECISION, PARAMETER :: CSA = 0.624
      DOUBLE PRECISION, PARAMETER :: CSB = 0.316
      DOUBLE PRECISION, PARAMETER :: CSC = 0.460
      DOUBLE PRECISION, PARAMETER :: CSD = 0.38
      DOUBLE PRECISION, PARAMETER :: CXD = 0.86
      CBRT(VX) = DABS(VX)**C3RD
* RHB is 'baryon density': 1 amu * number of baryons per cm3
      RHB = RHO/AVM
* Electron screening theory from Graboske, DeWitt, Grossman & Cooper (1973),
* for strong (ZA, ZB, ZC) are intermediate screening (ZD). The reaction
* dependent charge parameters are stored in CZA ... CZD.
      WC = DOT_PRODUCT(XA(10:18), VZ(1:9))
      WC = WC/NI
      WB = NE/NI
      WA = ZT*ZT/(WB*WB)
      XB = CBRT(DABS(WB))
      VL = CPL*DSQRT(DABS(NI)*RHB*DEXP(-3.0D0*TL))
      ZA = CSA*XB*VL**(2.0/3.0)
      ZB = CSB*XB*ZA
      ZC = CSC/(XB*XB)
      ZD = CSD*WC*WA*(VL/(WA*ZT))**CXD
* Reaction rates interpolated in T, mostly from Caughlan & Fowler (1988)
      TF = TL/CLN
      RRT(2:20) = 0.0d0
      RTT(2:20) = 0.0d0
      TT = 50.0D0*(TF - 6.0D0) + 1.0D0
      IF ( TT >= 1.0D0 ) THEN
         IT = MAX0(1, MIN0(199, INT(TT)))
         TT = TT - IT
         TU = 1.0D0 - TT
         DO IJ = 1, 19
            RR = TU*CRT(IT, IJ) + TT*CRT(IT + 1, IJ)
            IF ( RR >= -50.0D0 ) THEN
               SCRN = ZD*CZD(IJ)
               STRN = ZA*CZA(IJ) + ZB*CZB(IJ)
               DSTR = ZC*CZC(IJ)
               IF (DSTR < 0.29D0*STRN) SCRN = MIN(SCRN, STRN - DSTR)
               RRT(IJ + 1) = EXP(CLN*RR + SCRN)
! Derivative (dR/dlogT)_rho, by differentiation of the interpolation function
               RTT(IJ+1) = 50.0/CLN * (CRT(IT + 1, IJ) - CRT(IT, IJ))
            END IF
         END DO
      END IF
* Multiply with density and abundances to get rates per baryon per second,
* note that abundances of He3 and Be7 are not needed in equilibrium
      RPP = RHB*N1*N1*RPP/2.0D0
      R33 = RHB*R33/2.0D0
      R34 = RHB*N4*R34
      RBE = RHB*NE*RBE
      RBP = RHB*N1*RBP
      RPC = RHB*N1*N12*RPC
      RPN = RHB*N1*N14*RPN
      RPO = RHB*N1*N16*RPO
      R3A = RHB*RHB*N4*N4*N4*R3A/6.0D0
      RAC = RHB*N4*N12*RAC
      RAN = RHB*N4*N14*RAN
      RAO = RHB*N4*N16*RAO
      RANE = RHB*N4*N20*RANE
      RCC = RHB*N12*N12*RCC/2.0D0
      RCO = RHB*N12*N16*RCO
      ROO = RHB*N16*N16*ROO/2.0D0
      RGNE = N20*RGNE
      RGMG = N24*RGMG
* Branching of pN and CC reactions
      FPNG = 8.0D-4
      RPNA = (1.0D0 - FPNG)*RPN
      RPNG = FPNG*RPN
      RPN = RPNA
      FCCG = RCCG
      RCCA = (1.0D0 - FCCG)*RCC
      RCCG = FCCG*RCC
      RCC = RCCA
* PP chain in equilibrium, RPP becomes effective rate of 2 H1 -> 0.5 He4
      F34 = 0.0D0
      RPP = MIN(RPP, 1.0D303)
      R33 = MIN(R33, 1.0D303)
      R34 = MIN(R34, 1.0D303)
      IF (R34 > 1.0D-20)
     :     F34 = 2.0D0/(1.0D0 + DSQRT(1.0D0 + 8.0D0*RPP*R33/(R34*R34)))
      RPP = RPP*(1.0D0 + F34)
      PP2 = 1.0D0
      IF (RBE + RBP > 1.0D-20) PP2 = RBE/(RBE + RBP)
      PP3 = 1.0D0 - PP2
      QPP = QRT(1) + 0.5D0*QRT(2)
      QNPP = (QNT(1) + F34*(QNT(4)*PP2 + QNT(5)*PP3))/(1.0D0 + F34)
* calculate energy release and neutrino loss, in erg/gram/sec
      EX = QPP*RPP
      ENX = QNPP*RPP
      EXT = QPP*RPP*RTT(2)
      DO I = 6, 20
         EX = EX + QRT(I)*RRT(I + 1)
         ENX = ENX + QNT(I)*RRT(I + 1)
         EXT = EXT + QRT(I)*RRT(I + 1)*RTT(I + 1)
      END DO
      EXT = EXT / MAX(EX, 1.0D-50)
      ! Guard against overflows
      EX = MIN(CME*EX/AVM, 1.0D303)
      ENX = -MIN(CME*ENX/AVM, 1.0D303)
      RETURN
      END
      
! DIFFUSION_COEFFICIENTS:
! Compute atomic diffusion coefficients, according to Paquette et al. 1986
! Input:
!  * RHO:      Density [g/cm^3]
!  * T:        Temperature
!  * Abundances from the common block ABUND
! Output:
!  * Ddiff(9): Diffusion coefficients [cm^2/s] for the 9 isotopes H, He, C,
!              N, O, Ne, Mg, Si and Fe.
! This function only calculates the first order coefficients at the moment,
! but calculating the second order corrections (as well as the coefficient
! for thermodiffusion) is not difficult.
      SUBROUTINE DIFFUSION_COEFFICIENTS(RHO, T, Ddiff, NREF)
      USE CONSTANTS
      IMPLICIT NONE
! Input and output variables
      DOUBLE PRECISION, INTENT(IN) :: RHO, T
      DOUBLE PRECISION, INTENT(OUT) :: Ddiff(9)
      INTEGER, INTENT(IN) :: NREF            ! Reference abundance (H)
! Local variables
! Lengthscales for the plasma
      DOUBLE PRECISION :: LAMBDA_I2, LAMBDA_D2, LAMBDA2
! Baryon density
      DOUBLE PRECISION :: NB
! Convenience, to precalculate as much as possible
      DOUBLE PRECISION :: KT, GAMMA2
! Collision integrals Omega (eqn. (18) in Paquette et al)
      DOUBLE PRECISION :: OMEGA1(3), OMEGA22
! Dimensionless collision integrals (eqn. (65) in Paquette et al)
      DOUBLE PRECISION :: F1(3), F22
! Local variables
      DOUBLE PRECISION :: PSI_ST, EPS_ST, GAMMA_ST2, A, E_PSI_ST
      DOUBLE PRECISION :: DPSI_N1, DPSI_N
      INTEGER :: I, J, N
! COMMON block variables
!     /ABUND /
      DOUBLE PRECISION :: XA(9), NA(9), NEO, NIO, NZZ, AVM, NE
      COMMON /ABUND / XA, NA, NEO, NIO, NZZ, AVM, NE
!     /ATDATA/ 
      DOUBLE PRECISION :: CH2(4), CHI(26,9), COM(27), CAN(9), CBN(9)
      INTEGER :: KZN(9)
      COMMON /ATDATA/ CH2, CHI, COM, CAN, CBN, KZN
!     /COLINT/
      DOUBLE PRECISION :: DC(4, 50, 3), DD(4, 50)
      COMMON /COLINT/ DC, DD

      NB = RHO/(AVM*AMU)
      KT = BOLTZM*T

!     Typical distance between ions (squared)
      LAMBDA_I2 = (3.0/(CPI4*NB*NIO))**(2.0D0*C3RD)
!     Debye length (squared) 
      LAMBDA_D2 = (KT/(CPI4 * ECHAR**2 * NB * (NZZ+NE)))
!     Use max of lambda_i and lambda_D as typical distance
      LAMBDA2 = MAX(LAMBDA_D2, LAMBDA_I2)
!     Precompute some factors
      GAMMA2 = (4.0*KT/(ECHAR**2*KZN(NREF)))**2*LAMBDA2
      DO I=1, 9
         GAMMA_ST2 = GAMMA2/KZN(I)**2
!        Reduced mass of the two particles
         A = (CAN(I)*CAN(NREF)/(CAN(I)+CAN(NREF)))
!        Units of Omega: length^2 * velocity
         EPS_ST = CPI4*LAMBDA2/GAMMA_ST2 * SQRT(2.0D0*KT/(CPI4 * AMU*A))
!        Find interpolation interval for dimensionless collision integral
         E_PSI_ST = LOG(1.0D0+GAMMA_ST2)
         PSI_ST = LOG(E_PSI_ST)
!        Evaluate the collision integral, for repulsive coefficients
         IF (PSI_ST < 3.0D0) THEN
!           Use spline interpolation to evaluate the collision integrals
!           Determine interval in the table
            N = MIN(50, 1+FLOOR((PSI_ST +7d0)/0.2d0))
            DPSI_N1 = (-7.0D0 + 0.2D0*(N+1)) - PSI_ST
            DPSI_N = PSI_ST - (-7.0D0 + 0.2D0*N)
            DO J=1,3
               F1(J) = DEXP(DC(1,N,J)*DPSI_N1**3 + DC(2,N,J)*DPSI_N**3
     &                    + DC(3,N,J)*DPSI_N1    + DC(4,N,J)*DPSI_N)
            END DO
            F22 = DEXP(DD(1,N)*DPSI_N1**3 + DD(2,N)*DPSI_N**3
     &               + DD(3,N)*DPSI_N1    + DD(4,N)*DPSI_N)
         ELSE
            F1(1) = 1.00141D0*E_PSI_ST - 3.18209D0
            F1(2) = 0.99559D0*E_PSI_ST - 1.29553D0
            F1(3) = 1.99814D0*E_PSI_ST - 0.64413D0
            F22   = 1.99016D0*E_PSI_ST - 4.56958D0
         END IF
         OMEGA1(:) = EPS_ST*F1(:)
         OMEGA22   = EPS_ST*F22

!        Diffusion coefficient
         Ddiff(I) = 3.0*KT/(16.0*RHO*(NA(NREF)+NA(I))*A*OMEGA1(1))
      END DO

      END SUBROUTINE

c Read nuclear reaction (QRT) and neutrino (QNT) Q values, in MeV; constants 
c for electron screening (CZA, CZB, CZC, CZD, VZ); atomic parameters (CBN, KZN),
c with masses (CAN) consistent with Q-values; ionization potentials (CHI) and 
c statistical weights (COM); molecular hydrogen parameters (CH2)
      SUBROUTINE LOAD_ATOMIC_DATA(FILE)
      IMPLICIT REAL*8 (A-H, L-Z)
      INTEGER, INTENT(IN) :: FILE
      COMMON /ATDATA/ CH2(4), CHI(26,9), COM(27), CAN(9), CBN(9), KZN(9)
      COMMON /NCDATA/ QRT(20), QNT(20), CZA(92), CZB(92), CZC(92),
     &                CZD(92), VZ(9)
      INTEGER Z1(92), Z2(92), I, J
      DOUBLE PRECISION :: CXA,CXB,CXC,CXD
      DATA Z1 /1,2,2,4,4,6,7,8,4,6,7,8,10,6,8,8,10,12,0,0,1,1,1,0,1,1,1,
     &    2,2,2,2,1,1,2,1,1,2,2,1,1,2,1,2,2,1,2,2,0,0,1,1,1,1,1,4,1,1,1,
     &    4,4,4,1,1,1,1,4,4,4,0,0,0,1,0,0,0,1,0,0,0,1,1,1,4,4,4,4,1,1,1,
     &    1,1,1/
      DATA Z2 /1,2,2,0,1,1,1,1,2,2,2,2, 2,6,6,8, 0, 0,0,0,1,3,4,4,6,7,8,
     &     6,7,8,3,5,6,6,8,8,8,8,9,9,9,10,10,10,10,10,10,11,11,11,11,11,
     &     11,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,13,13,13,13,
     &     13,13,13,13,13,13,13,13,13,13,13,11,11,11,14,14,14,7,6,10/
  992 FORMAT (1P, 12(10D12.4,/), 32(9D12.4,/), 4D12.4,/, 9I12)

      READ (FILE,992) 
     & QRT, QNT, CZA(1:20), CZB(1:20), CZC(1:20), CZD(1:20),
     & VZ, CBN, CAN, COM, CHI, CH2, KZN 
      CLOSE (FILE)

!     Compute screening factors: mainly because not all of those that are
!     needed for the nucleosynthesis code are in the data file.
      CXA = 5.0/3.0
      CXC = CXA - 1.0
      CXB = 2.0*CXC
      CXD = 1.86
      DO J = 1, 92
         CZA(J) = (Z1(J)+Z2(J))**CXA - Z1(J)**CXA - Z2(J)**CXA
         CZB(J) = (Z1(J)+Z2(J))**CXB - Z1(J)**CXB - Z2(J)**CXB
         CZC(J) = -((Z1(J)+Z2(J))**CXC - Z1(J)**CXC - Z2(J)**CXC)
         CZD(J) = (Z1(J)+Z2(J))**CXD - Z1(J)**CXD - Z2(J)**CXD
      END DO
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Opacity tables and opacity related functions !
! Should really be a separate file...          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate a bicubic spline interpolation fit for the temperature
! and density opacity fit.  This is based on routines from Princeton
! Splines, D. McCune found by L. Dray.
      SUBROUTINE OPSPLN
      IMPLICIT REAL*8 (A-H, L-Z)
      PARAMETER (KT = 127,KR = 90)
      REAL*8 MAT(4,KT)
      COMMON /STAT1 / CSX(10), CS(90,127,10), CHAT(8920), JCSX
      COMMON /STAT3/ FSPL(4,4,127,90,10), TFM(127), FRM(90), WW(12)
      DO JT = 1, KT
         TFM(JT) = 2.95D0 + 0.05D0*DFLOAT(JT)
      END DO
      DO JR = 1, KR
         FRM(JR) = 0.25D0*DFLOAT(JR) - 12.25D0
      END DO
      DO JX = 1, 10
         DO JQ = 1, KT
             DO IQ = 1, KR
                FSPL(1, 1, JQ, IQ, JX) = CS(IQ, JQ, JX)
            END DO
         END DO
* Construct splines in the T direction.
         DO IR = 1, KR
            DO IT = 1, KT
               MAT(1, IT) = FSPL(1, 1, IT, IR, JX)
            END DO
            CALL SPLINE ( KT, TFM, MAT )
            DO IT = 1, KT - 1
               FSPL(2, 1, IT, IR, JX) = MAT(2, IT)
               FSPL(3, 1, IT, IR, JX) = MAT(3, IT)
               FSPL(4, 1, IT, IR, JX) = MAT(4, IT)
            END DO
         END DO
* Construct splines in the rho direction.
         DO IT = 1, KT - 1
* Construct splines for each T coeff
            DO IC = 1, 4
               DO IR = 1, KR
                  MAT(1, IR) = FSPL(IC, 1, IT, IR, JX)
               END DO
               MAT(2, 1) = 0.0D0
               MAT(3, 1) = 0.0D0
               MAT(2, KR) = 0.0D0
               MAT(3, KR) = 0.0D0
               CALL SPLINE ( KR, FRM, MAT )
               DO IR = 1, KR - 1
                  FSPL(IC, 2, IT, IR, JX) = MAT(2, IR)
                  FSPL(IC, 3, IT, IR, JX) = MAT(3, IR)
                  FSPL(IC, 4, IT, IR, JX) = MAT(4, IR)
               END DO
            END DO
         END DO
      END DO
      RETURN
      END



! As OPSPLN for the CO-enhanced tables SdM
! Calculate a bicubic spline interpolation fit for the temperature
! and density opacity fit.  This is based on routines from Princeton
! Splines, D. McCune found by L. Dray.
      SUBROUTINE OPSPLN_CO
      USE SETTINGS
      USE OPACITY_CO
      IMPLICIT NONE
      DOUBLE PRECISION ::  MAT(4,CO_MT)
      INTEGER :: I, J, K, IC, IC1, IC2

* define T and R intervals in opac tables
* NB: R = rho/T_6^4 and not rho as in the ld tables
      DO I=1,CO_MT
         opT(I) = 3.0D0 + 5.0D-2*DFLOAT(I-1)
      ENDDO
      DO J=1,CO_MR
         opR(J) = -8.0D0 + 5.0D-1*DFLOAT(J-1)
      ENDDO         
* Setup CO spline tables
      DO K=1,305
         DO J=1,CO_MR
            DO I=1,CO_MT
               MAT(1,I) = SPLINE_OPAC_CO(1,1,I,J,K)
            ENDDO
* Construct splines in the T direction
            CALL SPLINE(CO_MT, opT, MAT)
            DO I=1,CO_MT-1
               DO IC1=2,4
                  SPLINE_OPAC_CO(IC1,1,I,J,K) = MAT(IC1,I)
               ENDDO
            ENDDO               
         ENDDO         
         
* Construct splines in the R (= rho/T_6^3) direction
         DO I=1,CO_MT-1
            DO IC = 1,4
               DO J=1,CO_MR
                  MAT(1,J) = SPLINE_OPAC_CO(IC,1,I,J,K)
               ENDDO
               MAT(2,1)  = 0.0D0
               MAT(3,1)  = 0.0D0
               MAT(2,CO_MR) = 0.0D0
               MAT(3,CO_MR) = 0.0D0   
               
               CALL SPLINE(CO_MR, opR, MAT)            
               DO J=1,CO_MR-1
                  DO IC2 = 2,4                       
                     SPLINE_OPAC_CO(IC,IC2,I,J,K) = MAT(IC2,J)
                  ENDDO         !IC2                     
               ENDDO            !J           
            ENDDO               !IC 
         ENDDO                  !I
      ENDDO                     !K
      END



! Calculate the coefficients of a 1-D cubic spline:
! Forsythe, Malcolm, Moler, Computer Methods for Mathematical
! Computations, Prentice-Hall, 1977, p.76
      SUBROUTINE SPLINE ( K, X, F )
      IMPLICIT REAL*8 (A-H, L-Z)
      DIMENSION X(*),F(4,*)
      F(2:4,K) = 0.0D0
* Set up a tridiagonal system for A*y=B where y(i) are the second
* derivatives at the knots.
* f(2,i) are the diagonal elements of A
* f(4,i) are the off-diagonal elements of A
* f(3,i) are the B elements/3, and will become c/3 upon solution
      F(4,1) = X(2)-X(1)
      F(3,2) = (F(1,2) - F(1,1))/F(4,1)
      DO I = 2, K - 1
         F(4,I) = X(I+1) - X(I)
         F(2,I) = 2.0D0*(F(4,I-1) + F(4,I))
         F(3,I+1) = (F(1,I+1) - F(1,I))/F(4,I)
         F(3,I) = F(3,I+1) - F(3,I)
      END DO
* Boundaries.
      F(2,2) = F(4,1) + 2.0D0*F(4,2)
      F(3,2) = F(3,2)*F(4,2)/(F(4,1) + F(4,2))
      F(2,K-1) = 2.0D0*F(4,K-2) + F(4,K-1)
      F(3,K-1) = F(3,K-1)*F(4,K-2)/(F(4,K-1) + F(4,K-2))
* Forward elimination.
      T = F(4,2)/F(2,2)
      F(2,3) = F(2,3) - T*(F(4,2) - F(4,1))
      F(3,3) = F(3,3) - T*F(3,2)
      DO I = 4, K - 2
         T = F(4,I-1)/F(2,I-1)
         F(2,I) = F(2,I)-T*F(4,I-1)
         F(3,I) = F(3,I)-T*F(3,I-1)
      END DO
      T = (F(4,K-2) - F(4,K-1))/F(2,K-2)
      F(2,K-1) = F(2,K-1) - T*F(4,K-2)
      F(3,K-1) = F(3,K-1) - T*F(3,K-2)
* Back substitution.
      F(3,K-1) = F(3,K-1)/F(2,K-1)
      DO IB = 1, K - 4
         I = K - 1 - IB
         F(3,I) = (F(3,I) - F(4,I)*F(3,I+1))/F(2,I)
      END DO
      F(3,2) = (F(3,2) - (F(4,2) - F(4,1))*F(3,3))/F(2,2)
* Reset d array to step size.
      F(4,1) = X(2) - X(1)
      F(4,K-1) = X(K) - X(K-1)
* Set f(3,1) for not-a-knot.
      F(3,1) = (F(3,2)*(F(4,1) + F(4,2)) - F(3,3)*F(4,1))/F(4,2)
      F(3,K) = F(3,K-1) + (F(3,K-1) - F(3,K-2))*F(4,K-1)/F(4,K-2)
* Compute the polynomial coefficients.
      DO I = 1, K - 1
         F(2,I) = (F(1,I+1) - F(1,I))/F(4,I) - F(4,I)*(F(3,I+1) 
     :            + 2.0D0*F(3,I))
         F(4,I) = (F(3,I+1) - F(3,I))/F(4,I)
         F(3,I) = 3.0D0*F(3,I)
         F(4,I) = F(4,I)
      END DO
      RETURN
      END



! Read opacity tables
      SUBROUTINE LOAD_OPACITY(FILE)
      USE SETTINGS
      USE CONSTANTS
      IMPLICIT REAL*8 (A-H, L-Z)
      INTEGER, INTENT(IN) :: FILE
      COMMON /STAT1 / CSX(10), CS(90, 127, 10), CHAT(8920), KCSX
  990 FORMAT (1X, 10F7.3)
C Read opacity, nuclear reaction and neutrino loss rate data
      READ (FILE,*) KCSX, CZS, CH_OPAC
c Don't override CH value from init.dat (CH>0)!
      IF (CH<0.0D0) CH = CH_OPAC
      READ (FILE,990) CSX
      DO I = 1, KCSX
         READ (FILE,990) ((CS(JR, JT, I), JR = 1, 90), JT = 1, 127)
      END DO
! log (Z/Z_sun), used now and then, eg. Vink mass loss rate
      CLOGZ = LOG10(MAX(CZS/CZSN, 1.0D-40))
      REWIND (FILE)
c Set up coefficients for cubic spline interpolation in opacity
      CALL OPSPLN
      END



! Read CO enhanced opacity tables
      SUBROUTINE LOAD_OPACITY_CO(FILE)
      USE SETTINGS
      USE CONSTANTS
      USE OPACITY_CO
      IMPLICIT NONE
      INTEGER :: I,J,K
      DOUBLE PRECISION :: temp
      INTEGER, INTENT(IN) :: FILE
      
      IF (.NOT. ALLOCATED(SPLINE_OPAC_CO)) ALLOCATE(SPLINE_OPAC_CO(4,4,CO_MT,CO_MR,305))

 9901 FORMAT (E10.2)
 9902 FORMAT (F5.2, 31F7.3)
      READ (FILE,9901) CZS

      DO K=1,305
         READ (FILE,*)   
         DO I=1,CO_MT
            READ (FILE,9902) temp, (SPLINE_OPAC_CO(1,1,I,J,K), J=1,CO_MR)
         ENDDO
      ENDDO
      REWIND (FILE)
c Calculate initial H and basic CO abundance as scaled from solar
c Don't override CH value from init.dat (CH>0)!
      IF (CH<0.0D0) CH = 0.76D0 - 3.0D0*CZS
      CBASE = CZS*CC         
      OBASE = CZS*CO 
! log (Z/Z_sun), used now and then, eg. Vink mass loss rate
      CLOGZ = LOG10(MAX(CZS/CZSN, 1.0D-40))
c Set up coefficients for cubic spline interpolation in opacity
      CALL OPSPLN_CO
      END



      SUBROUTINE LOAD_REACTION_AND_NEUTRINO_RATES(FILE)
      IMPLICIT REAL*8 (A-H, L-Z)
      INTEGER, INTENT(IN) :: FILE
      COMMON /STAT1 / CSX(10), CS(90, 127, 10), CHAT(8920), KCSX
 991  FORMAT (1X, 10F7.3)
c Read nuclear reaction and neutrino loss rate data
      READ (FILE, 991) CHAT
      REWIND (FILE)
      END


 
* Calculate a bicubic spline interpolation fit for the temperature
* and density opacity fit.  Do not stop if the input lies outside the
* array but rather use the value at the nearest edge-point.
* Uses the older opacity tables
      SUBROUTINE OPACTY ( JX, TF, FRHO, FKL, FKH )
      IMPLICIT REAL*8 (A-H, L-Z)
      INTEGER MT, MR
      PARAMETER (MT = 127, MR = 90)
      COMMON /STAT1 / CSX(10), CS(90,127,10), CNU(60,41,2), W(4000),
     :                KCSX
      COMMON /STAT3 / F(4,4,127,90,10), TFM(127), FRM(90),
     :                FKLM(6), FKHM(6)
      XTF = DMAX1(DMIN1(TF, TFM(MT)), TFM(1))
      XFR = DMAX1(DMIN1(FRHO, FRM(MR)), FRM(1))
* Find interval in which target point lies.
         I1 = 1 + (MT - 1)*(XTF - TFM(1))/(TFM(MT) - TFM(1))
         I2 = 1 + (MR - 1)*(XFR - FRM(1))/(FRM(MR) - FRM(1))
         DT = TF - TFM(I1)
         DR = FRHO - FRM(I2)
* Evaluate the splines.
         FKL = F(1, 1, I1, I2, JX) + DR*(F(1, 2, I1, I2, JX)
     &   + DR*(F(1, 3, I1, I2, JX) + DR*F(1, 4, I1, I2, JX)))
     &   + DT*(F(2, 1, I1, I2, JX) + DR*(F(2, 2, I1, I2, JX)
     &   + DR*(F(2, 3, I1, I2, JX) + DR*F(2, 4, I1, I2, JX)))
     &   + DT*(F(3, 1, I1, I2, JX) + DR*(F(3, 2, I1, I2, JX)
     &   + DR*(F(3, 3, I1, I2, JX) + DR*F(3, 4, I1, I2, JX)))
     &   + DT*(F(4, 1, I1, I2, JX) + DR*(F(4, 2, I1, I2, JX)
     &   + DR*(F(4, 3, I1, I2, JX) + DR*F(4, 4, I1, I2, JX))))))
         FKH = F(1, 1, I1, I2, JX + 1) + DR*(F(1, 2, I1, I2, JX + 1)
     &   + DR*(F(1, 3, I1, I2, JX + 1) + DR*F(1, 4, I1, I2, JX + 1)))
     &   + DT*(F(2, 1, I1, I2, JX + 1) + DR*(F(2, 2, I1, I2, JX + 1)
     &   + DR*(F(2, 3, I1, I2, JX + 1) + DR*F(2, 4, I1, I2, JX + 1)))
     &   + DT*(F(3, 1, I1, I2, JX + 1) + DR*(F(3, 2, I1, I2, JX + 1)
     &   + DR*(F(3, 3, I1, I2, JX + 1) + DR*F(3, 4, I1, I2, JX + 1)))
     &   + DT*(F(4, 1, I1, I2, JX + 1) + DR*(F(4, 2, I1, I2, JX + 1)
     &   + DR*(F(4, 3, I1, I2, JX + 1) + DR*F(4, 4, I1, I2, JX + 1))))))
      RETURN      ! We normally don't want any output here...
      IF ( (XTF /= TF) .OR. (XFR /= FRHO) ) WRITE (10,100) TF, FRHO
      IF(FKL >  10. .OR. FKH >  10.) WRITE(10,*) "OPACITY TOO HIGH", FKL, FKH
      IF(FKL < -15. .OR. FKH < -15.) WRITE(10,*) "OPACITY TOO LOW", FKL, FKH
  100 FORMAT ('OPACITY OUT OF RANGE',' TF FRHO ', 2F9.4)
      RETURN
      END


      SUBROUTINE OPACTY_CO ( JJX, FT, FR, FKL, FKH )
c     -----------------------------------------------------------------------
c     Computes the opacity for a given composition by interpolating in
c     Temperature and density/R recently changed by SdM
c     -----------------------------------------------------------------------
c     INPUT variables - JX: index referring to a certain composition
c                     - FT: log10 T 
c                     - FR: log10 R = log (T6/rho^3)
c
c     OUTPUT Variables
c       - FKL: log opacity for nearest lower composition value in table
c       - FKH: log opacity for nearest higher composition value in table
c     The actual opacity is found by doing a linear interpolation between
c     FKL and FKH
c     
c     Uses:
c      + MODULE OPACITY_CO
c        - SPLINE_OPAC_CO: spline coefficients for CO enhanced opacity tables
c        - opT(1:CO_MT): log10 T at which CO opacities are defined
c        - opR(1:CO_MR): log10 R at which CO opacities are defined 
c      + COMMON /STAT3/
c        - F: spline coefficients for Pols opacity tables
c        - FTM(127): log10 T at which Pols opacities are defined
c        - FRM(90): log10 Rho at which Pols opacities are defined
c     -----------------------------------------------------------------------
      USE OPACITY_CO
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JJX
      DOUBLE PRECISION, INTENT(IN) :: FT, FR
      DOUBLE PRECISION, INTENT(OUT) :: FKL, FKH

      DOUBLE PRECISION :: XFT, XFR, DT, DR
      INTEGER :: JX, NEXT_JX, I1, I2

* Calculate a bicubic spline interpolation fit for the temperature
* and density opacity fit.  Do not stop if the input lies outside the
* array but rather use the value at the nearest edge-point.
         XFT = DMAX1(DMIN1(FT, opT(CO_MT)), opT(1))
         XFR = DMAX1(DMIN1(FR, opR(CO_MR)), opR(1))
* Find interval in which target point lies. Index I1 is index of
* nearest lower(!) entry in opT/opR
        
         I1 = 1 + int((CO_MT - 1)*(XFT - opT(1))/(opT(CO_MT) - opT(1)))
         I2 = 1 + int((CO_MR - 1)*(XFR - opR(1))/(opR(CO_MR) - opR(1)))

         DT = XFT - opT(I1)
         DR = XFR - opR(I2)         

* Evaluate the splines.      
         NEXT_JX = 61
         JX = JJX
         IF (JX+NEXT_JX>305) THEN
            JX = JX-NEXT_JX
         END IF

         FKL = SPLINE_OPAC_CO(1, 1, I1, I2, JX) + DR*(SPLINE_OPAC_CO(1, 2, I1, I2, JX)
     &   + DR*(SPLINE_OPAC_CO(1, 3, I1, I2, JX) + DR* SPLINE_OPAC_CO(1, 4, I1, I2, JX)))
     &   + DT*(SPLINE_OPAC_CO(2, 1, I1, I2, JX) + DR*(SPLINE_OPAC_CO(2, 2, I1, I2, JX)
     &   + DR*(SPLINE_OPAC_CO(2, 3, I1, I2, JX) + DR* SPLINE_OPAC_CO(2, 4, I1, I2, JX)))
     &   + DT*(SPLINE_OPAC_CO(3, 1, I1, I2, JX) + DR*(SPLINE_OPAC_CO(3, 2, I1, I2, JX)
     &   + DR*(SPLINE_OPAC_CO(3, 3, I1, I2, JX) + DR* SPLINE_OPAC_CO(3, 4, I1, I2, JX)))
     &   + DT*(SPLINE_OPAC_CO(4, 1, I1, I2, JX) + DR*(SPLINE_OPAC_CO(4, 2, I1, I2, JX)
     &   + DR*(SPLINE_OPAC_CO(4, 3, I1, I2, JX) + DR* SPLINE_OPAC_CO(4, 4, I1, I2, JX))))))

         FKH = SPLINE_OPAC_CO(1, 1, I1, I2, JX + NEXT_JX) + DR*(SPLINE_OPAC_CO(1, 2, I1, I2, JX + NEXT_JX)
     &   + DR*(SPLINE_OPAC_CO(1, 3, I1, I2, JX + NEXT_JX) + DR* SPLINE_OPAC_CO(1, 4, I1, I2, JX + NEXT_JX)))
     &   + DT*(SPLINE_OPAC_CO(2, 1, I1, I2, JX + NEXT_JX) + DR*(SPLINE_OPAC_CO(2, 2, I1, I2, JX + NEXT_JX)
     &   + DR*(SPLINE_OPAC_CO(2, 3, I1, I2, JX + NEXT_JX) + DR* SPLINE_OPAC_CO(2, 4, I1, I2, JX + NEXT_JX)))
     &   + DT*(SPLINE_OPAC_CO(3, 1, I1, I2, JX + NEXT_JX) + DR*(SPLINE_OPAC_CO(3, 2, I1, I2, JX + NEXT_JX)
     &   + DR*(SPLINE_OPAC_CO(3, 3, I1, I2, JX + NEXT_JX) + DR* SPLINE_OPAC_CO(3, 4, I1, I2, JX + NEXT_JX)))
     &   + DT*(SPLINE_OPAC_CO(4, 1, I1, I2, JX + NEXT_JX) + DR*(SPLINE_OPAC_CO(4, 2, I1, I2, JX + NEXT_JX)
     &   + DR*(SPLINE_OPAC_CO(4, 3, I1, I2, JX + NEXT_JX) + DR* SPLINE_OPAC_CO(4, 4, I1, I2, JX + NEXT_JX))))))

      RETURN      ! We normally don't want any output here...
      IF ( (XFT /= FT) .OR. (XFR /= FR) ) WRITE (10,100) FT, FR
   
 100  FORMAT ('OPACITY OUT OF RANGE',' FT FR ', 2F9.4)
      RETURN
      END



      FUNCTION GET_OPACITY(FRHO, FT)
c     ------------------------------------------------------------------
c     Calculates interpolated opacities as function of density and
c     temperature for the composition which is currently stored in
c     /ABUND/:  NA - Abundances, by number of all elements
c              AVM - Average mass [AMU] per baryon
c     INPUT:
c        FRHO - log10 rho, density
c        FT   - log10 T, temperature
c     RETURNS:
c        The total opacity, in cm^2/g
c
c     ------------------------------------------------------------------
c     Control variable:
c     KOP = 1: for opacities as implemented by Pols & al.
c           2: wrong! opacities assuming a never changing composition 
c           3: CO enhanced opacities, after Eldridge & Tout
c           4: as 3, while every element heavier than oxygen contributes
c              to the opacity as oxygen does
c     ------------------------------------------------------------------
      USE CONSTANTS
      USE SETTINGS
      USE OPACITY_CO
      IMPLICIT REAL*8 (A-H, L-Z)
      DOUBLE PRECISION, INTENT(IN)  :: FRHO, FT ! density, temperature
      COMMON /STAT1 / CSX(10), CS(90,127,10), CNU(60,41,2), W(4000),
     :                KCSX
      COMMON /ABUND / XA(9), NA(9), NEO, NIO, NZZ, AVM, NE
      COMMON /ATDATA/ CH2, C1, C2, C3, CHI(26,9), COM(27), CAN(9), 
     :                CBN(9), KZN(9)
      DOUBLE PRECISION :: TKAPPA (2,2)
      INTEGER :: NC, NO, IIC, IIO, JX=1

      XH = NA(1)*CAN(1)/AVM      ! mass fraction of Hydrogen
      XHE = NA(2)*CAN(2)/AVM     ! mass fraction of Helium

      IF (KOP == 1) THEN
* Opacity tables from Alexander & Ferguson (1994; molecular), Itoh (1983;
* electron conduction) and Iglesias & Rogers (1992; the rest)
c Find XF interval in opacity table.
         XF  = XH + XH + XHE + 1.0D0
         ! FIXME: should use a more efficient binary search here...
         DO JX=1,KCSX-1
            IF( CSX(JX+1)<XF ) EXIT
         ENDDO  
         XT =  (XF-CSX(JX))/(CSX(JX+1)-CSX(JX))
         XU = 1.0D0 - XT
         ! XT = MIN(1.0D0,(MAX(0.0D0,XF))) ! make sure XT between 0 and 1.
         ! XU = MIN(1.0D0,(MAX(0.0D0,XU))) ! make sure XU between 0 and 1.

         ! bicubic spline interpolation      
         CALL OPACTY ( JX, FT, FRHO, FKL, FKH)
         FK = XT*10.0D0**MIN(FKH, 3.0D2) + XU*10.0D0**MIN(FKL, 3.0D2)        
      ELSEIF (KOP > 1 .OR. KOP < 0) THEN     ! KOP<0 for debug purposes
c ------------------------------------------------------------------------
c Use new opacity tables as implemented after ELdridge and Tout MNRAS 2003
c ------------------------------------------------------------------------
c OPAL tables (Iglesias & Rogers 1996) for different H, C, and O
c fractions extended to lower temperatures with Alexander & Ferguson
c 1994 (! or more up to date?) (only for different H fractions) and
c extended to higher temperatures with Buchler & Yuehc 1976. Finally the
c opacity owing to electron conduction is added according to Hubbard &
c Lamp (1969) and Canuto 1970

c The new opacity tables use the variable FR rather than the density
         FR = FRHO - 3.0D0*(FT -6.0D0) ! FR = log10 (Density/T_6^3)
c Determine in which interval of XH we are and prepare for linear interpolation
c FIXME: no good reason not to use a binary search here...
         DO JX=1,CO_MH-1
            IF ( opH(JX+1)>XH ) EXIT
         ENDDO
c JX will be the index of the element in opH with nearest lower H
c abundance on which the opacity table is defined. Note that opH has a
c descending H abundance
         JX = MIN(JX, CO_MH-1)
         XT =  (XH-opH(JX))/(opH(JX+1)-opH(JX))
         XT = MIN(1.0D0,(MAX(0.0D0,XT))) ! make sure XT between 0 and 1.
! ?? Could we extrapolate for X > 0.7, at low Z
         XU = 1.0D0 - XT
         XU = MIN(1.0D0,(MAX(0.0D0,XU))) ! make sure XU between 0 and 1.

c Composition parameters used for interpolation: (FC=0: Carbon 
c abundance as solar but scaled with metallicity, FC=1 max enhancement)
         XXC = NA(3)*CAN(3)/AVM  ! mass fraction of carbon
         XXO = NA(5)*CAN(5)/AVM  ! mass fraction of oxygen
         FC = (XXC - CBASE)/(1.0D0-CZS-XH)
         FO = (XXO - OBASE)/(1.0D0-CZS-XH)         
         FC = Min(1.0D0, Max(0.0D0, FC))
         FO = Min(1.0D0, Max(0.0D0, FO))

c Include every element heavier than oxygen in that abundance
! slightly different from JJE's code, but it seems right to me=SdM
         
         IF (KOP >= 4.AND.XH<1.0D-10) THEN  
            XXO = 1.0D0 - XH- XHE -XXC            
            FO = (XXO - OBASE)/(1.0D0-CZS-XH)         
         ENDIF

c For KOP=2, or if the enhancement of C and O is very small the opacity
c is taken for a mixture where C and O are never enhanced and only CBASE
c and OBASE scale with Z
         IF(KOP == 2 .OR. (FC<1.0D-6 .AND. FO<1.0D-6)) THEN   
            CALL OPACTY_CO(1+(JX-1)*61, FT, FR, FKL_CO, FKH_CO)
            FK_CO = XT*10.0D0**MIN(FKH_CO,3.0D2) + 
     &              XU*10.0D0**MIN(FKL_CO,3.0D2)
            FK = FK_CO
         ELSE                   !Do full 5d interpolation for KOP>=3
c for very small  C and O, neglect them
            IF (FC < 1.0D-20) FC = 0.0D0
            IF (FO < 1.0D-20) FO = 0.0D0
c Find relevant CO fraction interval in opacity table
            IIC = 1
            IIO = 1
            DO K2=1,7
               IF(COcompos(K2) < FC)  IIC = K2                          
               IF(COcompos(K2) < FO)  IIO = K2                  
            ENDDO          

            IIC = MAX(MIN(7,IIC), 1)
            IIO = MAX(MIN(7,IIO),1)
c For each FC there are tables with 8 different FO's except for 
c FC = 0.6 (7 tables) and FC = 1.0 (6 tables)
c For construction of the table see explanation JJE, 10 October 2003
c "implementation of Enhanced CO Opacity Table", website
c IIO and IIC are swapped compared to JJE's implementation
            IF(IIC<=7) THEN
               JJX = (JX-1)*61 + (IIC-1)*8 + IIO
            ELSE               
               JJX = (JX-1)*61 +       6*8 + 7 + IIO
            ENDIF

c Calculate opacity at all four lower and upper levels for C and O abundance
c Varying O direction +/- 1, Varying in C direction +/- 8 (or 7)
c IIO and IIC are swapped compared to JJE's implementation
            DO NC=0,1
               DO NO = 0,1
                  IF(IIC <= 6) THEN
                     JK = JJX + 8*NC + NO
                  ELSE
                     JK = JJX + 7*NC + NO
                  ENDIF            
    
                  CALL OPACTY_CO(JK, fT,fR,FKL_CO,FKH_CO)
                  FK_CO = XT*10.0D0**MIN(FKH_CO,3.0D2) + 
     &                    XU*10.0D0**MIN(FKL_CO,3.0D2)
c Store opacities of nearest point in the table for later interpolation
                  tkappa(NC+1,NO+1) = FK_CO
               ENDDO
            ENDDO
c Final 2D lineair interpolation for C and O            
            eC = ( fC-COcompos(IIC) )/( COcompos(IIC+1)-COcompos(IIC) )
            eO = ( fO-COcompos(IIO) )/( COcompos(IIO+1)-COcompos(IIO) )
           
            eC = Min(Max(eC, 1.0D-50), 1.0D0)
            eO = Min(Max(eO, 1.0D-50), 1.0D0)

            FK_CO = (1.0D0-eC) * (1.0D0-eO) * tkappa(1,1) +
     &                     eC  * (1.0D0-eO) * tkappa(2,1) +
     &              (1.0D0-eC) *        eO  * tkappa(1,2) +
     &                     eC  *        eO  * tkappa(2,2)

            FK = FK_CO
         ENDIF
      ENDIF
      GET_OPACITY = FK
      END
