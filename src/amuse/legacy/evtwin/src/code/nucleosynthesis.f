!     Subroutines for nucleosynthesis post-procesing
!     originally by Richard Stancliffe (Stancliffe & al. 2005 MNRAS).
      MODULE NUCLEOSYNTHESIS
      USE MESH

      INTEGER, PARAMETER :: NVAR_NUC = 50
      INTEGER, PARAMETER :: NUC_NUM_RATES = 92
      INTEGER, SAVE :: KH_NUC = 0

!     Global switch to determine whether nucleosynthesis has been switched
!     on or not
      LOGICAL, SAVE :: NUCLEOSYNTHESIS_ENABLED = .FALSE.

!     Input quantities needed to calculate the reaction rates for the
!     nucleosynthesis code
      DOUBLE PRECISION, SAVE :: HT(2,18,NM)

!     Nucleosynthesis data structures. These are the analogues of the H and
!     DH arrays.
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: HNUCPR(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: HNUC(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: DHNUC(:,:,:)

!     Nucleosynthesis reaction rates
      DOUBLE PRECISION, SAVE :: NUCSYN_CRT(200,NUC_NUM_RATES)
      DOUBLE PRECISION, SAVE :: NUCSYN_NRT(200,45)

!     Misc. options
!     Short-circuit decay of meta-stable Al26. Shouldn't matter but messes
!     up convergence some times.
      LOGICAL, PARAMETER :: INSTANTLY_DECAY_AL26M = .TRUE.

      END MODULE

      MODULE ZAMS_NUCLEOSYNTHESIS_ABUNDANCES
!     Default ZAMS abundance mixture, as a fraction of Z.
!     Based on Anders & Grevesse (1989)
      DOUBLE PRECISION, SAVE :: CXD     = 2.406e-03
      DOUBLE PRECISION, SAVE :: CXHE3   = 1.468e-03
      DOUBLE PRECISION, SAVE :: CXLI7   = 4.687e-07
      DOUBLE PRECISION, SAVE :: CXBE7   = 0.000e+00
      DOUBLE PRECISION, SAVE :: CXB11   = 2.368e-07
      DOUBLE PRECISION, SAVE :: CXC12   = 1.764e-01
      DOUBLE PRECISION, SAVE :: CXC13   = 1.829e-03
      DOUBLE PRECISION, SAVE :: CXC14   = 0.000e+00
      DOUBLE PRECISION, SAVE :: CXN14   = 5.212e-02
      DOUBLE PRECISION, SAVE :: CXN15   = 2.186e-04
      DOUBLE PRECISION, SAVE :: CXO16   = 5.031e-01
      DOUBLE PRECISION, SAVE :: CXO18   = 1.086e-03
      DOUBLE PRECISION, SAVE :: CXO17   = 1.948e-04
      DOUBLE PRECISION, SAVE :: CXF19   = 2.030e-05
      DOUBLE PRECISION, SAVE :: CXNE21  = 2.068e-04
      DOUBLE PRECISION, SAVE :: CXNE20  = 9.221e-02
      DOUBLE PRECISION, SAVE :: CXNE22  = 6.525e-03
      DOUBLE PRECISION, SAVE :: CXNA22  = 0.000e+00
      DOUBLE PRECISION, SAVE :: CXNA23  = 1.673e-03
      DOUBLE PRECISION, SAVE :: CXMG24  = 2.580e-02
      DOUBLE PRECISION, SAVE :: CXMG25  = 3.391e-03
      DOUBLE PRECISION, SAVE :: CXMG26  = 3.889e-03
      DOUBLE PRECISION, SAVE :: CXAL26M = 0.000e+00
      DOUBLE PRECISION, SAVE :: CXAL27  = 2.906e-03
      DOUBLE PRECISION, SAVE :: CXAL26G = 0.000e+00
      DOUBLE PRECISION, SAVE :: CXSI28  = 3.272e-02
      DOUBLE PRECISION, SAVE :: CXSI30  = 1.179e-03
      DOUBLE PRECISION, SAVE :: CXSI29  = 1.717e-03
      DOUBLE PRECISION, SAVE :: CXP31   = 4.087e-03
      DOUBLE PRECISION, SAVE :: CXS33   = 1.615e-04
      DOUBLE PRECISION, SAVE :: CXS32   = 1.984e-02
      DOUBLE PRECISION, SAVE :: CXS34   = 9.351e-04
      DOUBLE PRECISION, SAVE :: CXFE57  = 1.431e-03
      DOUBLE PRECISION, SAVE :: CXFE60  = 0.000e+00
      DOUBLE PRECISION, SAVE :: CXFE56  = 5.858e-02
      DOUBLE PRECISION, SAVE :: CXFE58  = 1.853e-04
      DOUBLE PRECISION, SAVE :: CXFE59  = 0.000e+00
      DOUBLE PRECISION, SAVE :: CXCO59  = 1.683e-04
      DOUBLE PRECISION, SAVE :: CXNI61  = 4.307e-05
      DOUBLE PRECISION, SAVE :: CXNI59  = 0.000e+00
      DOUBLE PRECISION, SAVE :: CXNI58  = 2.478e-03
      DOUBLE PRECISION, SAVE :: CXNI60  = 9.812e-04
      END MODULE

      SUBROUTINE ALLOCATE_NUCLEOSYNTHESIS_DATA(KH)
      USE NUCLEOSYNTHESIS
      USE MESH
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KH
      
      IF (.NOT. NUCLEOSYNTHESIS_ENABLED) RETURN

      IF (NVAR_NUC > NVAR) THEN
         WRITE(*, *) 'Stellar structure/solver code does not hold enough data to solve nucleosynthesis'
         WRITE(*, *) 'Increase the value of NVAR in mesh.f'
         STOP
      END IF
      IF (ALLOCATED(HNUCPR)) DEALLOCATE(HNUCPR)
      IF (ALLOCATED(HNUC)) DEALLOCATE(HNUC)
      IF (ALLOCATED(DHNUC)) DEALLOCATE(DHNUC)
      ALLOCATE(HNUCPR(2, NVAR_NUC, KH))
      ALLOCATE(HNUC(2, NVAR_NUC, KH))
      ALLOCATE(DHNUC(2, NVAR_NUC, KH))

      KH_NUC = KH
      END SUBROUTINE

      SUBROUTINE SET_INITIAL_NUCLEOSYNTHESIS_ABUNDANCES
      USE ZAMS_NUCLEOSYNTHESIS_ABUNDANCES
      USE NUCLEOSYNTHESIS
      USE CONSTANTS
      USE SETTINGS
      IMPLICIT NONE
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EP(3)
      INTEGER :: KH, KTW, KW(260)
      COMMON H, DH, EP, KH, KTW, KW
!     /ABUND /
      DOUBLE PRECISION :: XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, YA(9), WNE0, WA(2), AVM, WNE
      COMMON /ABUND / XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, YA, WNE0, WA, AVM, WNE
!     /ATDATA /
      DOUBLE PRECISION :: CH2(4), CHI(26,9), COM(27), CAN(9), CBN(9)
      INTEGER ::  KZN(9)
      COMMON /ATDATA/ CH2, CHI, COM, CAN, CBN, KZN

      INTEGER :: JK, JSTAR, I
      DOUBLE PRECISION :: XD, XHE3, XLI7, XBE7, XB11, XC12, XC13, XC14,
     &                  XN14, XN15, XO16,XO17, XO18, XF19, XNE20, XNE21, XNE22, 
     &                  XNA22, XNA23,XMG24, XMG25, XMG26, XAL26M, XAL26G, 
     &                  XAL27, XSI28,XSI29, XSI30, XP31, XS32, XS33, XS34, 
     &                  XFE56, XFE57,XFE58, XFE59, XFE60, XCO59, XNI58, XNI59,
     &                  XNI60, XNI61 
      DOUBLE PRECISION :: AF, AT

      IF (.NOT. NUCLEOSYNTHESIS_ENABLED) RETURN

C Beginning of minor variables, in mass order (except gallinoes)
      HNUC(1:2, 1 , 1:KH) = 0d0       ! gallinoes
      HNUC(1:2, 2 , 1:KH) = 0d0       ! neutrons
      HNUC(1:2, 3 , 1:KH) = CXD * CZS
      HNUC(1:2, 4 , 1:KH) = CXHE3 * CZS
      HNUC(1:2, 5 , 1:KH) = CXLI7 * CZS
      HNUC(1:2, 6 , 1:KH) = CXBE7 * CZS
      HNUC(1:2, 7 , 1:KH) = CXB11 * CZS
      HNUC(1:2, 8 , 1:KH) = CXC13 * CZS
      HNUC(1:2, 9 , 1:KH) = CXC14 * CZS
      HNUC(1:2, 10, 1:KH) = CXN15 * CZS
      HNUC(1:2, 11, 1:KH) = CXO17 * CZS
      HNUC(1:2, 12, 1:KH) = CXO18 * CZS
      HNUC(1:2, 13, 1:KH) = CXF19 * CZS
      HNUC(1:2, 14, 1:KH) = CXNE21 * CZS
      HNUC(1:2, 15, 1:KH) = CXNE22 * CZS
      HNUC(1:2, 16, 1:KH) = CXNA22 * CZS
      HNUC(1:2, 17, 1:KH) = CXNA23 * CZS
      HNUC(1:2, 18, 1:KH) = CXMG24 * CZS
      HNUC(1:2, 19, 1:KH) = CXMG25 * CZS
      HNUC(1:2, 20, 1:KH) = CXMG26 * CZS
      HNUC(1:2, 21, 1:KH) = CXAL26M * CZS
      HNUC(1:2, 22, 1:KH) = CXAL26G * CZS
      HNUC(1:2, 23, 1:KH) = CXAL27 * CZS
      HNUC(1:2, 24, 1:KH) = CXSI28 * CZS
      HNUC(1:2, 25, 1:KH) = CXSI29 * CZS
      HNUC(1:2, 26, 1:KH) = CXSI30 * CZS
      HNUC(1:2, 27, 1:KH) = CXP31 * CZS
      HNUC(1:2, 28, 1:KH) = CXS32 * CZS
      HNUC(1:2, 29, 1:KH) = CXS33 * CZS
      HNUC(1:2, 30, 1:KH) = CXS34 * CZS
      HNUC(1:2, 31, 1:KH) = CXFE56 * CZS
      HNUC(1:2, 32, 1:KH) = CXFE57 * CZS
      HNUC(1:2, 33, 1:KH) = CXFE58 * CZS
      HNUC(1:2, 34, 1:KH) = CXFE59 * CZS
      HNUC(1:2, 35, 1:KH) = CXFE60 * CZS
      HNUC(1:2, 36, 1:KH) = CXCO59 * CZS
      HNUC(1:2, 37, 1:KH) = CXNI58 * CZS
      HNUC(1:2, 38, 1:KH) = CXNI59 * CZS
      HNUC(1:2, 39, 1:KH) = CXNI60 * CZS
      HNUC(1:2, 40, 1:KH) = CXNI61 * CZS
! Remesh using inputs - not necessarily consistent with above
      DO JSTAR=1, 2
         DO JK = 1, KH
            XH = H(24*(JSTAR-1)+5, JK)
            XHE = H(24*(JSTAR-1)+9, JK)
            XC = H(24*(JSTAR-1)+10, JK)
            XN = H(24*(JSTAR-1)+16, JK)
            XO = H(24*(JSTAR-1)+3, JK)
            XNE = H(24*(JSTAR-1)+11, JK)
            AF = H(24*(JSTAR-1)+1, JK)
            AT = H(24*(JSTAR-1)+2, JK)
            CALL STATEF ( AF, AT )
            DO I=1,6
               HNUC(JSTAR, 40+I, JK) = CAN(I)*YA(I)/AVM
            END DO
         END DO
! For high core temperatures, log10 T > 6.7, set D, Li7 to 0 (should all
! have burned on the pre-main sequence)
         IF (AT > 6.5*CLN) THEN
            HNUC(JSTAR, 3 , 1:KH) = 0.0d0
            HNUC(JSTAR, 4 , 1:KH) = 0.0d0
            HNUC(JSTAR, 5 , 1:KH) = 0.0d0
            HNUC(JSTAR, 6 , 1:KH) = 0.0d0
            HNUC(JSTAR, 7 , 1:KH) = 0.0d0
         END IF
      END DO

      END SUBROUTINE

      SUBROUTINE LOAD_NUCLEOSYNTHESIS_RATES ( RATES_IN, NRATES_IN )
      USE NUCLEOSYNTHESIS
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: RATES_IN, NRATES_IN
      DOUBLE PRECISION :: CAT, CRT
      INTEGER :: KCSX
      COMMON /STAT1 / CAT(119230), CRT(200,20), KCSX
  990 FORMAT (1X, 10F7.3)

      READ (RATES_IN, 990) NUCSYN_CRT
      READ (NRATES_IN, 990) NUCSYN_NRT
      REWIND (RATES_IN)
      REWIND (NRATES_IN)
!     Make sure the first 20 reaction rates (structure and nucleosynthesis)
!     agree
      NUCSYN_CRT(:,1:20) = CRT(:, 1:20)
      END SUBROUTINE


      SUBROUTINE FUNCS2 ( K, K1, K2 )
      USE MESH
      USE CONSTANTS
      USE SETTINGS
      USE NUCLEOSYNTHESIS
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON /TVBLES/ DT, ZQ(243), JHOLD, JM2, JM1
      COMMON /STAT2 / AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA,
     & SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, STAT2_RATES(20), EX, ENX,
     & WMU, DELTA, PHIE, EXT, FKT, FKR, PRANDTL
      COMMON /INF   / W(NVAR),DW(NVAR),WJUNK(NFUNC),WJUNK2(NVAR,NFUNC)
! RJS common block to get nucleosynthesis variables
      COMMON /INF2  / WN(50),DWN(50),YDOT(50),Y(50),SG,DMT,WN2(8),
     :     YLIN(50,110)
      COMMON /ABUND / XA1(9), NA(9), VNE0, WA(2), AVM, VNE
      COMMON /ABUND2/ XA2(50), WW2(50)
C Rates from nucrat2 - format is Rabc
C a=light reactant,b=heavy reactant,c=light product, blank if photon
      COMMON /RATES / RPP, R33, R34, RBE, RBP, RPC, RPN, RPO, R3A, RAC, !10
     :     RAN, RAO, RANE, RCC, RCO, ROO, RGNE, RGMG, RCCG, RPNG,       !20
     :     RPD,RPBE,RPLI,REBE,RPC13,RPN15A,RPO17a,                      !27
     :     RAC13n, RAN15g, RAO17n, RLIA, RPB11, RPC14, RAC14,           !34
     :     RPO18,RPO18a,RAO18,RAO18n,RPF19,RPF19A,RAF19p,RPNE21,        !42
     :     RANE21,RANE21n,RPNE22,RANE22,RANE22n,RnNA22p,RnNA22a,RPNA22, !50
     :     RPNA23,RPNA23n,RPNA23a,RPMG24,RAMG24,RPMG25T,RPMG25M,RPMG25G,!58
     :     RAMG25,RAMG25n,RAMG25p,RPMG26,RPMG26Tn,RPMG26Mn,RPMG26Gn,    !65
     :     RAMG26,RAMG26n,RAMG26p,Rg26Tp,Rn26Tp,Rn26Ta,Rp26T,Rg26Mp,    !73
     :     Rn26Mp,Rn26Ma,Rp26M,Rg26Gp,Rn26Gp,Rn26Ga,Rp26G,RPAL27,       !81
     :     RPAL27a,RAAL27n,RANA23nT,RANA23nM,RANA23nG,RPSI28,RPSI29,    !88
     :     RPSI30,RPN15,RPO17,RPNE20                                    !92
      COMMON /NRATES/ RnFe56,RnFe57,RnFe58,RnCo59,RnNi58,RnNi59,RnNi60,
     :     Rnp,RnHe3,RnLi7,RnC12,RnC13,RnC14,RnN14,RnN15,RnO16,RnO18,
     :     RnF19,RnNe20,RnNe21,RnNe22,RnNa23,RnMg24,RnMg25,RnMg26,
     :     RnAl27,RnSi28,RnSi29,RnSi30,RnP31,RnS32,RnS33,RnS34,Rn26G,
     :     RnS33a,RnN14p,RnNi59p,RnNi59a,RnO17a,Rn26Gad,RnS32a,RnFe59,
     :     RnFe60,RnS34s,RnNi61s
      COMMON /DECAY / RDAL26G,RDNA22,RD26M,RDFe59,RDFe60,RDNi59,RDn,RDC14
      DIMENSION BARYN(50)
      DATA BARYN/1.0, 1.0, 2.0, 3.0, 7.0, 7.0, 11.0, 13.0, 14.0,  15.0,
     : 17.0, 18.0, 19.0, 21.0, 22.0, 22.0, 23.0, 24.0, 25.0, 26.0,
     : 26.0,26.0,27.0,28.0,29.0,30.0,31.0,32.0,33.0,34.0,56.0,57.0,58.0,
     : 59.0,60.0,59.0,58.0,59.0,60.0,61.0,1.0,4.0,12.0,14.0,16.0,20.0,
     : 1.0,1.0,1.0,1.0/
      COMMON /NEUTRO/ FRAC(50,NM)
C Set up composition for reaction calculation
      XA1(1) = W(5)  !X1
      XA1(2) = W(9)  !X4
      XA1(3) = W(10) !X12
      XA1(4) = W(16) !X14
      XA1(5) = W(3)  !X16
      XA1(6) = W(11) !X20
      XA1(7) = CMG*CZS
      XA1(8) = CSI*CZS
      XA1(9) = CFE*CZS
C Temperature
      AT = W(2)
C Number abundance stuff from FUNCS1
      VNE0 = HT(1,5,K)
      WA(1) = HT(1,6,K)
      WA(2) = HT(1,7,K)
      AVM = HT(1,8,K)
      VNE = HT(1,9,K)
      NA(1:9) = HT(1, 10:18, K)
C Blank any -ve abundances. There shouldn't be any anyway...
      DO I = 1, 9
         IF (XA1(I).LT.0d0) THEN
            XA1(I) = 0d0
            NA(I) = 0d0
         END IF
      END DO
      XA2(:) = WN(:)
      XA2(2) = 0d0 ! neutrons are annoying - we treat them differently
C Proton fudging: if structure says there are none, there are none
      IF (XA1(1).LT.1d-40) XA2(41) = 0d0
C Work out abundances for nucrat2 - RJS
      DO I = 1, 50
         WW2(I) = XA2(I)/BARYN(I)
      END DO
      IF ( K.LE.K1 ) THEN
         YLIN(1:50, 1:110) = 0.0d0
      END IF
      DAM = DW(4)
      RHO = HT(1,1, K)
      SG = HT(1,2, K)
      DM = W(4)
      DMT = DAM/DT
      DMK = HT(1,4, K)
      CALL NUCRAT2 ( AT )
C  EVALUATE COMPOSITION CHANGE
C 1/9/03 Minor element reactions from Clayton 1963
      YDOT(1) = 0d0
      YDOT(2) =  RDn - RAC13n - RAO17n - RANE21n - RAO18n  - RANE22n  +
     :     RnNA22p + RnNA22a - RPNA23n - RAMG25n - RPMG26Mn - RPMG26Gn
     :     - RAMG26n + Rn26Mp + Rn26Ma + Rn26G
     :     + Rn26Gp + Rn26Ga - RAAL27n - RANA23nM - RANA23nG
     :     + RnFe56 + RnFe57 + RnFe58 + RnCo59 + RnNi58 + RnNi59
     :     + RnNi60 + Rnp + RnHe3 + RnLi7 + RnC12 + RnC13 + RnC14
     :     + RnN14 + RnN15 + RnO16 + RnO18 + RnF19 + RnNe20 + RnNe21
     :     + RnNe22 + RnNa23 + RnMg24 + RnMg25 + RnMg26 + RnAl27
     :     + RnSi28 + RnSi29 + RnSi30 + RnP31 + RnS32 + RnS33 + RnS34
     :     + RnS33a + RnN14p + RnNi59p + RnNi59a + RnO17a + RnS32a
      YDOT(3) = RPD - RPP - Rnp
      YDOT(4) = - RPD + 2.0*R33 + R34 + RnHe3
      YDOT(5) =  RPLI - REBE + RLIA + RnLi7
      YDOT(6) = -R34 + RPBE + REBE
      YDOT(7) = RPB11 - RLIA
      YDOT(8) = RPC13 - RPC + RAC13n + RnC13 - RnC12
      YDOT(9) = RPC14 + RAC14 + RnC14 - RnN14p - RnO17a - RnC13 + RDC14
      YDOT(10)= RPN15A - RPN + RAN15g + RnN15 -  RnC14 - RnN14 - RPO18a
     :     + RPN15 - RPC14
      YDOT(11)= - RPO + RPO17a + RAO17n - RnO16 + RnO17a + RPO17
      YDOT(12)= - RAN - RAC14 + RPO18 + RPO18a + RAO18 + RAO18n + RnO18
     :     - RPO17
      YDOT(13)= - RAN15g + RPF19 + RPF19A + RAF19p - RnNA22a - RnO18
     :     + RnF19 - RPO18
      YDOT(14)= - RAO18n + RPNE21 + RANE21 + RANE21n + RnNe21 - RnNe20
     :     - RPNE20
      YDOT(15)= - RAO18 - RAF19p + RPNE22 + RANE22 + RANE22n - RnNA22p
     :     - RDNA22 - RnNe21 + RnNe22
      YDOT(16)= RnNA22p + RnNA22a + RPNA22 - RPNE21 + RDNA22
      YDOT(17)= RPNA23 + RPNA23n + RPNA23a - RPNE22 - Rn26Ga
     :     - Rn26Ma + RANA23nM + RANA23nG - RnNe22 + RnNa23
      YDOT(18)= RPMG24 + RAMG24 - RCC - RANE - RANE21n - RPNA23
     :     - RPAL27a - RnNa23 + RnMg24
      YDOT(19)= RPMG25M + RPMG25G + RAMG25 + RAMG25n + RAMG25p
     :     - RANE21 - RANE22n - Rg26Mp - Rg26Gp + RnMg25 - RnMg24
      YDOT(20)= RPMG26 + RPMG26Mn + RPMG26Gn
     :     + RAMG26 + RAMG26n + RAMG26p - RANE22 - Rn26Mp
     :     - Rn26Gp - RDAL26G - RD26M + RnMg26 - RnMg25
      YDOT(21)= Rg26Mp + Rn26Mp + Rn26Ma + Rp26M - RPMG25M - RPMG26Mn
     :     - RANA23nM + RD26M
      YDOT(22)= Rg26Gp + Rn26Gp + Rn26Ga + Rp26G - RPMG25G - RPMG26Gn
     :     - RANA23nG + RDAL26G + Rn26G
      YDOT(23)= RPAL27 + RPAL27a + RAAL27n - RPMG26 - RnMg26 + RnAl27
     :     - Rn26G
C assume decay of Si28(p,g)P29 -> Si29
      YDOT(24)= RPSI28 - RCO - RAMG24 - RAMG25n - RPAL27 - RnAl27 +
     :     RnSi28
C assume decay of Si29(p,g)P30 -> Si30, Mg26(a,p)Al29 -> Si29
      YDOT(25)= RPSI29 - RAMG25 - RAMG26n - RPSI28 - RnSi28 + RnSi29
     :          - RnS32a - RAMG26p
      YDOT(26)= RPSI30 - RAMG26 - RPSI29 - RnSi29 + RnSi30 - RnS33a
      YDOT(27)= - RPSI30 - RnSi30 + RnP31
      YDOT(28)= - RnP31 + RnS32 + RnS32a
      YDOT(29)= - RnS32 + RnS33 + RnS33a
      YDOT(30)= - RnS33 + RnS34 + RnS34s
      YDOT(31)= RnFe56 - RnNi59a
      YDOT(32)= RnFe57 - RnFe56
      YDOT(33)= RnFe58 - RnFe57
      YDOT(34)= - RnFe58 + RnFe59 + RDFe59
      YDOT(35)= - RnFe59 + RnFe60 + RDFe60
      YDOT(36)= RnCo59 - RDFe59 - RnNi59p - RDNi59
      YDOT(37)= RnNi58
      YDOT(38)= RnNi59 - RnNi58 + RnNi59a + RnNi59p + RDNi59
      YDOT(39)= RnNi60 - RnNi59
      YDOT(40)= RnNi61s - RnNi60
C Protons duplicate
      YDOT(41)= 2.0*RPP + RPD + Rnp + RPLI + RPBE + RPB11 - 2.0*R33
     :     + RPC + RPC13 + RPC14 + RPN + RPN15A + RPO18a + RPO
     :     + RPO17a - RDn + RPO18 + RPF19 + RPF19A - RAF19p + RPNE21
     :     + RPNE22 - RnNa22p + RPNA22 + RPNA23 + RPNA23n + RPNA23a
     :     + RPMG24 + RPAL27a + RPMG25M + RPMG25G + RPMG26 + RPN15
     :     + RPMG26Mn + RPMG26Gn - RAMG26p - RAMG25p
     :     + Rp26M + Rp26G + RPAL27 + RPSI28 + RPSI29 + RPSI30
     :     + RPO17 + RPNE20
C Helium duplicate - RPO17a is O17(p,a)N14
      YDOT(42)= - R33 + R34 - 2.0*RPLI - 2.0*RPBE + RLIA + RAC13n + RAN
     :     - RPN15A - RPO18a + RAC14 + RAN15g + RAO17n + RAO + 3.0*R3A
     :     + RAC + RAO18 + RAO18n - RPO17a - RPF19A + RAF19p - RCC
     :     + RANE21 + RANE21n + RANE22 + RANE22n - RnNa22a - RPNA23a
     :     - Rn26Ga - Rn26Ma + RANA23nM + RANA23nG + RANE + RAMG24
     :     - RPAL27a + RAMG25 + RAMG25n + RAMG25p + RAMG26 + RAMG26n
     :     + RAMG26p + RAAL27n - RnS32a - RnS33a - RnNi59a
C Carbon duplicate - RCC to alpha and Ne20 - should be ~50-50 with
C this and Na23 + p
      YDOT(43)= - RPB11 + RPC - RPN15A - R3A + RAC + 2.0*RCC + RnC12 + RCO
C Nitrogen duplicate
      YDOT(44)= - RPC13 + RnN14p + RPN + RAN - RPO17a
C Oxygen duplicate
      YDOT(45)= - RAC13n + RPO + RAO - RAC + RnO16 - RPF19A + RCO - RPN15
C Neon duplicate
      YDOT(46)= - RAO - RPF19 - RnF19 - RCC + RnNe20 + RANE - RPNA23a + RPNE20
C Make sure we don't run into problems when calculating the derivative for
C rates for which the abundance is 0 (the rate will be zero, but we
C calculate the derivative by rate/X, so X should be non-zero).
      DO I=1, 50
         IF (XA2(I).LE.0d0) XA2(I) = -9.9d0
      END DO
C  EVALUATE DERIVATIVES
C YLIN(derived by,equation)
C RJS 1/9/03 - Added species for PP-I,II chain from Clayton 1963
C 1st eq gallinoes

C 2nd eq neutrons
      YLIN(2, 2) = (RDn + RnNA22p + RnNA22a + Rn26Mp + Rn26Ma +
     :              Rn26Gp + Rn26Ga + RnFe56 + RnFe57 + RnFe58 +
     :              RnCo59 + RnNi58 + RnNi59 + RnNi60 + Rnp +
     :              RnHe3 + RnLi7 + RnC12 + RnC13 + RnC14 + RnN14 +
     :              RnN15 + RnO16 + RnO18 + RnF19 + RnNe20 + RnNe21 +
     :              RnNe22 + RnNa23 + RnMg24 + RnMg25 + RnMg26 +
     :              RnAl27 + RnSi28 + RnSi29 + RnSi30 + RnP31 + RnS32
     :              + RnS33 + RnS34 + Rn26G+ RnS33a + RnN14p + RnNi59p
     :              + RnNi59a + RnO17a + RnS32a)/XA2(2)
      YLIN(4, 2) = RnHe3/XA2(4)
      YLIN(5, 2) = RnLi7/XA2(5)
      YLIN(8, 2) = (RnC13- RAC13n)/XA2(8)
      YLIN(9, 2) = RnC14/XA2(9)
      YLIN(10,2) = RnN15/XA2(10)
      YLIN(11,2) = (RnO17a - RAO17n)/XA2(11)
      YLIN(12,2) = (RnO18 - RAO18n)/XA2(12)
      YLIN(13,2) = RnF19/XA2(13)
      YLIN(14,2) = (RnNe21 - RANE21n)/XA2(14)
      YLIN(15,2) = (RnNe22 - RANE22n)/XA2(15)
      YLIN(16,2) = (RnNA22p + RnNA22a)/XA2(16)
      YLIN(17,2) = (RnNa23 - RPNA23n - RANA23nM - RANA23nG)/XA2(17)
      YLIN(18,2) = RnMg24/XA2(18)
      YLIN(19,2) = (RnMg25 - RAMG25n)/XA2(19)
      YLIN(20,2) = (RnMg26 - RAMG26n - RPMG26Mn - RPMG26Gn)/XA2(20)
      YLIN(21,2) = (Rn26Mp + Rn26Ma)/XA2(21)
      YLIN(22,2) = (Rn26Gp + Rn26Ga + Rn26G)/XA2(22)
      YLIN(23,2) = (- RAAL27n + RnAl27)/XA2(23)
      YLIN(24,2) = RnSi28/XA2(24)
      YLIN(25,2) = RnSi29/XA2(25)
      YLIN(26,2) = RnSi30/XA2(26)
      YLIN(27,2) = RnP31/XA2(27)
      YLIN(28,2) = (RnS32 + RnS32a)/XA2(28)
      YLIN(29,2) = (RnS33 + RnS33a)/XA2(29)
      YLIN(30,2) = RnS34/XA2(30)
      YLIN(31,2) = RnFe56/XA2(31)
      YLIN(32,2) = RnFe57/XA2(32)
      YLIN(33,2) = RnFe58/XA2(33)
      YLIN(36,2) = RnCo59/XA2(36)
      YLIN(37,2) = RnNi58/XA2(37)
      YLIN(38,2) = (RnNi59 + RnNi59p + RnNi59a)/XA2(38)
      YLIN(39,2) = RnNi60/XA2(39)
      YLIN(41,2) = (Rnp - RPNA23n - RPMG26Mn - RPMG26Gn)/XA2(41)
      YLIN(42,2) = - (RAC13n + RAO17n + RAO18n + RANE21n + RANE22n
     :     + RANA23nM + RANA23nG + RAMG25n + RAMG26n + RAAL27n)/XA2(42)
      YLIN(43,2) = RnC12/XA2(43)
      YLIN(45,2) = RnO16/XA2(45)
      YLIN(46,2) = RnNe20/XA2(46)
C 3rd eq D
      YLIN(2, 3) = - Rnp/XA2(2)
      YLIN(3, 3) = RPD/XA2(3)
      YLIN(41,3) = (RPD - Rnp - 2.0*RPP)/XA2(41)
C 4th eq He-3
      YLIN(2, 4) = RnHe3/XA2(2)
      YLIN(3, 4) = -RPD/XA2(3)
      YLIN(4, 4) = (4.0*R33+R34+RnHe3)/XA2(4)
      YLIN(41,4) = - RPD/XA2(41)
      YLIN(42,4) = R34/XA2(42)
C 5th eq Li-7
      YLIN(2, 5) = RnLi7/XA2(2)
      YLIN(5, 5) = (RPLI+RLIA+RnLi7)/XA2(5)
      YLIN(6, 5) = -REBE/XA2(6)
      YLIN(41,5) = RPLI/XA2(41)
      YLIN(42,5) = RLIA/XA2(42)
C 6th eq Be-7
      YLIN(4, 6) = -R34/XA2(4)
      YLIN(6, 6) = (RPBE + REBE)/XA2(6)
      YLIN(41,6) = RPBE/XA2(41)
      YLIN(42,6) = - R34/XA2(42)
C 7th eq B-11
      YLIN(5, 7) = - RLIA/XA2(5)
      YLIN(7, 7) = RPB11/XA2(7)
      YLIN(41,7) = RPB11/XA2(41)
      YLIN(42,7) = - RLIA/XA2(42)
C CNO elements. Will assume beta-decay is instantaneous (which is ok on all but the very shortest timescales)
C 8th eq C-13
      YLIN(2, 8) = (RnC13 - RnC12)/XA2(2)
      YLIN(8, 8) = (RPC13+RAC13n+RnC13)/XA2(8)
      YLIN(41,8) = (RPC13 - RPC)/XA2(41)
      YLIN(42,8) = RAC13n/XA2(42)
      YLIN(43,8) = (- RPC - RnC12)/XA2(43)
C 9th eq C-14
      YLIN(2, 9) = (RnC14 - RnC13 - RnN14p - RnO17a)/XA2(2)
      YLIN(8, 9) = - RnC13/XA2(8)
      YLIN(9, 9) = (RPC14 + RAC14 + RnC14 + RDC14)/XA2(9)
      YLIN(11,9) = - RnO17a/XA2(11)
      YLIN(41,9) = RPC14/XA2(41)
      YLIN(42,9) = RAC14/XA2(42)
      YLIN(44,9) = - RnN14p/XA2(44)
C 10th eq N-15
C Assume instantaneous C15 -> N15 beta decay
      YLIN(2 ,10) = (RnN15 - RnC14 - RnN14)/XA2(2)
      YLIN(9 ,10) = -(RnC14 + RPC14)/XA2(9)
      YLIN(10,10) = (RPN15A+RAN15g+RnN15+RPN15)/XA2(10)
      YLIN(12,10) = - RPO18a/XA2(12)
      YLIN(41,10) = (RPN15A - RPN - RPC14 + RPO18a+RPN15)/XA2(41)
      YLIN(42,10) = RAN15g/XA2(42)
      YLIN(44,10) = (- RnN14 - RPN)/XA2(44)
C 11th eq O-17
      YLIN(2 ,11) = (RnO17a - RnO16)/XA2(2)
      YLIN(11,11) = (RPO17a+RAO17n+RnO17a+RPO17)/XA2(11)
      YLIN(41,11) = (RPO17a+RPO17)/XA2(41)
      YLIN(42,11) = RAO17n/XA2(42)
      YLIN(45,11) = - RnO16/XA2(45)
C 12th eq O18
      YLIN(2, 12) = RnO18/XA2(2)
      YLIN(9 ,12) = - RAC14/XA2(9)
      YLIN(11,12) = - RPO17/XA2(11)
      YLIN(12,12) = (RPO18 + RPO18a + RAO18 + RAO18n + RnO18)/XA2(12)
      YLIN(41,12) = (RPO18 + RPO18a - RPO17)/XA2(41)
      YLIN(42,12) = (RAO18 + RAO18n - RAN)/XA2(42)
      YLIN(44,12) = - RAN/XA2(44)
C 13th eq F19
      YLIN(2 ,13) = (RnF19 - RnNA22a - RnO18)/XA2(2)
      YLIN(10,13) = - RAN15g/XA2(10)
      YLIN(12,13) = - (RnO18 + RPO18)/XA2(12)
      YLIN(13,13) = (RPF19 + RPF19A + RAF19p + RnF19)/XA2(13)
      YLIN(16,13) = - RnNA22a/XA2(16)
      YLIN(41,13) = (RPF19 + RPF19A - RPO18)/XA2(41)
      YLIN(42,13) = (RAF19p - RAN15g)/XA2(42)
C 14th eq Ne21
      YLIN(2, 14) = (RnNe21-RnNe20)/XA2(2)
      YLIN(12,14) = - RAO18n/XA2(12)
      YLIN(14,14) = (RPNE21 + RANE21 + RANE21n + RnNe21)/XA2(14)
      YLIN(41,14) = (RPNE21 - RPNE20)/XA2(41)
      YLIN(42,14) = (RANE21 + RANE21n - RAO18n)/XA2(42)
      YLIN(46,14) = -RPNE20/XA2(46)
C 15th eq Ne22
      YLIN(2 ,15) = (RnNe22 - RnNe21 - RnNA22p)/XA2(2)
      YLIN(12,15) = - RAO18/XA2(12)
      YLIN(13,15) = - RAF19p/XA2(13)
      YLIN(14,15) = -RnNe21/XA2(14)
      YLIN(15,15) = (RPNE22 + RANE22 + RANE22n + RnNe22)/XA2(15)
      YLIN(16,15) = - (RnNA22p + RDNA22)/XA2(16)
      YLIN(41,15) = RPNE22/XA2(41)
      YLIN(42,15) = (RANE22 + RANE22n - RAO18 - RAF19p)/XA2(42)
C 16th eq Na22
      YLIN(2 ,16) = (RnNA22p + RnNA22a)/XA2(2)
      YLIN(14,16) = - RPNE21/XA2(14)
      YLIN(16,16) = (RnNA22p + RnNA22a + RPNA22 + RDNA22)/XA2(16)
      YLIN(41,16) = (RPNA22 - RPNE21)/XA2(41)
C 17th eq Na23
      YLIN(2 ,17) = (RnNa23 - Rn26Ma - Rn26Ga - RnNe22)/XA2(2)
      YLIN(15,17) = - (RPNE22+RnNe22)/XA2(15)
      YLIN(17,17) = (RPNA23 + RPNA23n + RPNA23a + RANA23nM
     :     + RANA23nG + RnNa23)/XA2(17)
      YLIN(21,17) = - Rn26Ma/XA2(21)
      YLIN(22,17) = - Rn26Ga/XA2(22)
      YLIN(41,17) = (RPNA23 + RPNA23n + RPNA23a - RPNE22)/XA2(41)
      YLIN(42,17) = (RANA23nM + RANA23nG)/XA2(42)
C 18th eq Mg24
      YLIN(2 ,18) = (RnMg24 - RnNa23)/XA2(2)
      YLIN(14,18) = - RANE21n/XA2(14)
      YLIN(17,18) = (- RnNa23 - RPNA23)/XA2(17)
      YLIN(18,18) = (RPMG24 + RAMG24 + RnMg24)/XA2(18)
      YLIN(23,18) = - RPAL27a/XA2(23)
      YLIN(41,18) = (RPMG24 - RPAL27a - RPNA23)/XA2(41)
      YLIN(42,18) = (RAMG24 - RANE - RANE21n)/XA2(42)
      YLIN(46,18) = - RANE/XA2(46)
C 19th eq Mg25
      YLIN(2 ,19) = (RnMg25 - RnMg24)/XA2(2)
      YLIN(14,19) = - RANE21/XA2(14)
      YLIN(15,19) = - RANE22n/XA2(15)
      YLIN(18,19) = - RnMg24/XA2(18)
      YLIN(19,19) = (RPMG25M + RPMG25G + RAMG25 + RAMG25n
     :      + RAMG25p + RnMg25)/XA2(19)
      YLIN(21,19) = - Rg26Mp/XA2(21)
      YLIN(22,19) = - Rg26Gp/XA2(22)
      YLIN(41,19) = (RPMG25M + RPMG25G)/XA2(41)
      YLIN(42,19) = (RAMG25 + RAMG25n + RAMG25p)/XA2(42)
C 20th eq Mg26
      YLIN(2 ,20) = (RnMg26 - RnMg25 - Rn26Mp - Rn26Gp)/XA2(2)
      YLIN(15,20) = - RANE22/XA2(15)
      YLIN(19,20) = - RnMg25/XA2(19)
      YLIN(20,20) = (RPMG26 + RPMG26Mn + RPMG26Gn
     :     + RAMG26 + RAMG26n + RAMG26p + RnMg26)/XA2(20)
      YLIN(21,20) = - (Rn26Mp + RD26M)/XA2(21)
      YLIN(22,20) = - (Rn26Gp + RDAL26G)/XA2(22)
      YLIN(41,20) = (RPMG26 + RPMG26Mn + RPMG26Gn)/XA2(41)
      YLIN(42,20) = (RAMG26 + RAMG26n + RAMG26p)/XA2(42)
C 21st eq Al26M
      YLIN(2 ,21) = (Rn26Mp + Rn26Ma)/XA2(2)
      YLIN(17,21) = - RANA23nM/XA2(17)
      YLIN(19,21) = - RPMG25M/XA2(19)
      YLIN(20,21) = - RPMG26Mn/XA2(20)
      YLIN(21,21) = (Rg26Mp + Rn26Mp + Rn26Ma + Rp26M + RD26M)/XA2(21)
      YLIN(41,21) = (Rp26M - RPMG25M - RPMG26Mn)/XA2(41)
      YLIN(42,21) = - RANA23nM/XA2(42)
C 22nd eq Al26G
      YLIN(2 ,22) = (Rn26Gp + Rn26Ga + Rn26G)/XA2(2)
      YLIN(17,22) = - RANA23nG/XA2(17)
      YLIN(19,22) = - RPMG25G/XA2(19)
      YLIN(20,22) = - RPMG26Gn/XA2(20)
      YLIN(22,22) = (Rg26Gp + Rn26Gp + Rn26Ga + Rp26G + RDAL26G
     :     + Rn26G)/XA2(22)
      YLIN(41,22) = (Rp26G - RPMG25G - RPMG26Gn)/XA2(41)
      YLIN(42,22) = - RANA23nG/XA2(42)
C 23rd eq Al27
      YLIN(2 ,23) = (RnAl27 - RnMg26 - Rn26G)/XA2(2)
      YLIN(20,23) = (-RnMg26 - RPMG26)/XA2(20)
      YLIN(22,23) = -Rn26G/XA2(22)
      YLIN(23,23) = (RPAL27 + RPAL27a + RAAL27n + RnAl27)/XA2(23)
      YLIN(41,23) = (RPAL27 + RPAL27a - RPMG26 - Rp26G - Rp26M)/XA2(41)
      YLIN(42,23) = RAAL27n/XA2(42)
C 24th eq Si28
      YLIN(2, 24) = (RnSi28 - RnAl27)/XA2(2)
      YLIN(18,24) = - RAMG24/XA2(18)
      YLIN(19,24) = - RAMG25n/XA2(19)
      YLIN(23,24) = (-RnAl27 - RPAL27)/XA2(23)
      YLIN(24,24) = (RnSi28 + RPSI28)/XA2(24)
      YLIN(41,24) = (RPSI28 - RPAL27)/XA2(41)
      YLIN(42,24) = (- RAMG24 - RAMG25n)/XA2(42)
      YLIN(43,24) = - RCO/XA2(43)
C 25th eq Si29
      YLIN(2 ,25) = (RnSi29 - RnSi28 - RnS32a)/XA2(2)
      YLIN(19,25) = - RAMG25/XA2(19)
      YLIN(20,25) = - (RAMG26n + RAMG26p)/XA2(20)
      YLIN(24,25) = - (RnSi28 + RPSI28)/XA2(24)
      YLIN(25,25) = (RPSI29 + RnSi29)/XA2(25)
      YLIN(28,25) = - RnS32a/XA2(28)
      YLIN(41,25) = (RPSI29 - RPSI28)/XA2(41)
      YLIN(42,25) = - (RAMG25 + RAMG26n + RAMG26p)/XA2(42)
C 26th eq Si30
      YLIN(2 ,26) = (RnSi30 - RnSi29 - RnS33a)/XA2(2)
      YLIN(20,26) = - RAMG26/XA2(20)
      YLIN(25,26) = - (RPSI29 + RnSi29)/XA2(25)
      YLIN(26,26) = (RPSI30+RnSi30)/XA2(26)
      YLIN(29,26) = - RnS33a/XA2(29)
      YLIN(41,26) = (RPSI30 - RPSI29)/XA2(41)
      YLIN(42,26) = - RAMG26/XA2(42)
C 27th eq P31
      YLIN(2 ,27) = (RnP31 - RnSi30)/XA2(2)
      YLIN(26,27) = - (RPSI30 + RnSi30)/XA2(26)
      YLIN(27,27) = RnP31/XA2(27)
      YLIN(41,27) = - RPSI30/XA2(41)
C 28th eq S32
      YLIN(2 ,28) = (RnS32 - RnP31 + RnS32a)/XA2(2)
      YLIN(27,28) = - RnP31/XA2(27)
      YLIN(28,28) = (RnS32 + RnS32a)/XA2(28)
C 29th eq S33
      YLIN(2 ,29) = (RnS33 - RnS32 + RnS33a)/XA2(2)
      YLIN(28,29) = - RnS32/XA2(28)
      YLIN(29,29) = (RnS33 + RnS33a)/XA2(29)
C 30th eq S34
      YLIN(2 ,30) = (RnS34 - RnS33 + RnS34s)/XA2(2)
      YLIN(29,30) = - RnS33/XA2(29)
      YLIN(30,30) = (RnS34 + RnS34s)/XA2(30)
C 31st eq Fe56
      YLIN(2 ,31) = (RnFe56 - RnNi59a)/XA2(2)
      YLIN(31,31) = RnFe56/XA2(31)
      YLIN(38,31) = -RnNi59a/XA2(38)
C 32nd eq Fe57
      YLIN(2 ,32) = (RnFe57 - RnFe56)/XA2(2)
      YLIN(31,32) = - RnFe56/XA2(31)
      YLIN(32,32) = RnFe57/XA2(32)
C 33rd eq Fe58
      YLIN(2 ,33) = (RnFe58-RnFe57)/XA2(2)
      YLIN(32,33) = - RnFe57/XA2(32)
      YLIN(33,33) = RnFe58/XA2(33)
C 34th eq Fe59
      YLIN(2 ,34) = (RnFe59 - RnFe58)/XA2(2)
      YLIN(33,34) = - RnFe58/XA2(33)
      YLIN(34,34) = (RnFe59 + RDFe59)/XA2(34)
C 35th eq Fe60
      YLIN(2 ,35) = (RnFe60 - RnFe59)/XA2(2)
      YLIN(34,35) = - RnFe59/XA2(34)
      YLIN(35,35) = (RnFe60 + RDFe60)/XA2(35)
C 36th eq Co59
      YLIN(2 ,36) = (RnCo59 - RnNi59p)/XA2(2)
      YLIN(34,36) = - RDFe59/XA2(34)
      YLIN(36,36) = RnCo59/XA2(36)
      YLIN(38,36) = (-RnNi59p - RDNi59)/XA2(38)
C 37th eq Ni58
      YLIN(2 ,37) = RnNi58/XA2(2)
      YLIN(37,37) = RnNi58/XA2(37)
C 38th eq Ni59
      YLIN(2 ,38) = (RnNi59 - RnNi58 + RnNi59p + RnNi59a)/XA2(2)
      YLIN(37,38) = - RnNi58/XA2(37)
      YLIN(38,38) = (RnNi59+ RnNi59p + RnNi59a + RDNi59)/XA2(38)
C 39th eq Ni60
      YLIN(2 ,39) = (RnNi60 - RnNi59)/XA2(2)
      YLIN(38,39) = - RnNi59/XA2(38)
      YLIN(39,39) = RnNi60/XA2(39)
C 40th eq Ni61
      YLIN(2 ,40) = (RnNi61s - RnNi60)/XA2(2)
      YLIN(39,40) = - RnNi60/XA2(39)
      YLIN(40,40) = RnNi61s/XA2(40)
C 41st eq H - duplicate
      YLIN(2 ,41) = (Rnp - RnN14p - RDn - Rn26Mp - Rn26Gp - RnNi59p)/XA2(2)
      YLIN(3 ,41) = RPD/XA2(3)
      YLIN(4 ,41) = - 4.0*R33/XA2(4)
      YLIN(5 ,41) = RPLI/XA2(5)
      YLIN(6 ,41) = RPBE/XA2(6)
      YLIN(7 ,41) = RPB11/XA2(7)
      YLIN(8 ,41) = RPC13/XA2(8)
      YLIN(9 ,41) = RPC14/XA2(9)
      YLIN(10,41) = (RPN15A+RPN15)/XA2(10)
      YLIN(11,41) = (RPO17a+RPO17)/XA2(11)
      YLIN(12,41) = RPO18a/XA2(12)
      YLIN(13,41) = (RPF19 + RPF19A  - RAF19p)/XA2(13)
      YLIN(14,41) = RPNE21/XA2(14)
      YLIN(15,41) = RPNE22/XA2(15)
      YLIN(16,41) = RPNA22/XA2(16)
      YLIN(17,41) = (RPNA23 + RPNA23n + RPNA23a)/XA2(17)
      YLIN(18,41) = RPMG24/XA2(18)
      YLIN(19,41) = (RPMG25M + RPMG25G - RAMG25p)/XA2(19)
      YLIN(20,41) = (RPMG26 + RPMG26Mn + RPMG26Gn - RAMG26p)/XA2(20)
      YLIN(21,41) = (Rp26M - Rn26Mp - Rg26Mp)/XA2(21)
      YLIN(22,41) = (Rp26G - Rn26Gp - Rg26Gp)/XA2(22)
      YLIN(23,41) = (RPAL27 + RPAL27a)/XA2(23)
      YLIN(24,41) = RPSI28/XA2(24)
      YLIN(25,41) = RPSI29/XA2(25)
      YLIN(26,41) = RPSI30/XA2(26)
      YLIN(38,41) = - RnNi59p/XA2(38)
      YLIN(41,41) = (4.0*RPP + RPD + Rnp + RPLI + RPBE + RPB11 - RAF19p
     :     + RPC13 + RPC14 + RPN + RPN15A + RPO18a + RPO + RPO17 + RPF19
     :     + RPF19A + RPNE21 + RPNE22 + RPNA22 + RPNA23 + RPNA23n
     :     + RPNA23a + RPMG24 + RPAL27a + RPMG25M + RPMG25G + RPMG26
     :     + RPMG26Mn + RPMG26Gn + Rp26M + Rp26G + RPAL27 + RPSI28
     :     + RPSI29 + RPSI30 + RPN15 + RPO17 + RPNE20)/XA2(41)
      YLIN(42,41) = - (RAF19p + RAMG25p + RAMG26p)/XA2(42)
      YLIN(43,41) = RPC/XA2(43)
      YLIN(44,41) = (RPN - RnN14p)/XA2(44)
      YLIN(45,41) = RPO/XA2(45)
      YLIN(46,41) = RPNE20/XA2(46)
C 42nd eq He4 - duplicate
      YLIN(2 ,42) = (- Rn26Ga - Rn26Ma - RnNa22a - RnS32a - RnS33a
     :     - RnNi59a)/XA2(2)
      YLIN(4 ,42) = (- 4.0*R33 + R34)/XA2(4)
      YLIN(5 ,42) = (RLIA - 2.0*RPLI)/XA2(5)
      YLIN(6 ,42) = -2.0*RPBE/XA2(6)
      YLIN(8 ,42) = RAC13n/XA2(8)
      YLIN(9 ,42) = RAC14/XA2(9)
      YLIN(10,42) = (RAN15g - RPN15A)/XA2(10)
      YLIN(11,42) = (RAO17n - RPO17a)/XA2(11)
      YLIN(12,42) = - RPO18a/XA2(12)
      YLIN(13,42) = (RAF19p - RPF19A)/XA2(13)
      YLIN(14,42) = (RANE21 + RANE21n)/XA2(14)
      YLIN(15,42) = (RANE22 + RANE22n)/XA2(15)
      YLIN(16,42) = - RnNa22a/XA2(16)
      YLIN(17,42) = (RANA23nM + RANA23nG)/XA2(17)
      YLIN(18,42) = RAMG24/XA2(18)
      YLIN(19,42) = (RAMG25 + RAMG25n + RAMG25p)/XA2(19)
      YLIN(20,42) = (RAMG26 + RAMG26n + RAMG26p)/XA2(20)
      YLIN(21,42) = - Rn26Ma/XA2(21)
      YLIN(22,42) = - Rn26Ga/XA2(22)
      YLIN(23,42) = (RAAL27n - RPAL27a)/XA2(23)
      YLIN(28,42) = - RnS32a/XA2(28)
      YLIN(29,42) = - RnS33a/XA2(29)
      YLIN(38,42) = - RnNi59a/XA2(38)
      YLIN(41,42) = (-2.0*RPLI - 2.0*RPBE - RPN15A - RPO18a - RPO17a
     :     - RPF19A - RPNA23a - RPAL27a)/XA2(41)
      YLIN(42,42) = (R34+RLIA+RAC13n+RAN + RAC14 + RAN15g + RAO17n
     :     + RAO + 9.0*R3A + RAC + RAF19p + RANE21 + RANE21n + RANE22
     :     + RANE22n + RANA23nM + RANA23nG + RANE + RAMG24 + RAMG25
     :     + RAMG25n + RAMG25p + RAMG26 + RAMG26n + RAMG26p + RAAL27n)
     :     /XA2(42)
      YLIN(43,42) = (RAC - 2.0*RCC)/XA2(43)
      YLIN(44,42) = RAN/XA2(44)
      YLIN(45,42) = RAO/XA2(45)
      YLIN(46,42) = RANE/XA2(46)
C 43rd eq C12 - duplicate
      YLIN(7 ,43) = - RPB11/XA2(7)
      YLIN(10,43) = - RPN15A/XA2(10)
      YLIN(41,43) = (RPC - RPB11 - RPN15A)/XA2(41)
      YLIN(42,43) = (RAC - 9.0*R3A)/XA2(42)
      YLIN(43,43) = (RPC + RAC + 4.0*RCC + RCO)/XA2(43)
      YLIN(45,43) = RCO/XA2(45)
C 44th eq N14 - duplicate
      YLIN(2 ,44) = RnN14p/XA2(2)
      YLIN(8 ,44) = - RPC13/XA2(8)
      YLIN(11,44) = - RPO17a/XA2(11)
      YLIN(41,44) = (RPN - RPC13 - RPO17a)/XA2(41)
      YLIN(44,44) = (RPN + RAN + RnN14p)/XA2(44)
C 45th eq O16 - duplicate
      YLIN(2, 45) = RnO16/XA2(2)
      YLIN(8 ,45) = - RAC13n/XA2(8)
      YLIN(10,45) = - RPN15/XA2(10)
      YLIN(13,45) = - RPF19A/XA2(13)
      YLIN(41,45) = (RPO - RPF19A - RPN15)/XA2(41)
      YLIN(42,45) = (RAO - RAC13n - RAC)/XA2(42)
      YLIN(43,45) = (RCO - RAC)/XA2(43)
      YLIN(45,45) = (RPO + RAO + RnO16 + RCO)/XA2(45)
C 46th eq Ne20 - duplicate
      YLIN(2 ,46) = (RnNe20 - RnF19)/XA2(2)
      YLIN(13,46) = (-RPF19 - RnF19)/XA2(13)
      YLIN(17,46) = - RPNA23a/XA2(17)
      YLIN(41,46) = (RPNE20 - RPF19 - RPNA23a)/XA2(41)
      YLIN(42,46) = (- RAO + RANE)/XA2(42)
      YLIN(43,46) = - 2.0*RCC/XA2(43)
      YLIN(45,46) = - RAO/XA2(45)
      YLIN(46,46) = (RnNe20 + RANE + RPNE20)/XA2(46)

!     Shortcut the cycling of Al26M -> Mg26
      IF (INSTANTLY_DECAY_AL26M) THEN
         YDOT(20) = YDOT(20) - YDOT(21)
         YLIN(:, 20) = YLIN(:,20) - YLIN(:, 21)
         YDOT(21) = 0.0d0
         YLIN(:, 21) = 0.0d0
      END IF
C Remove any NaN issues, by blanking rates if abundance = 0
      DO I=1, 50
         IF (XA2(I).LE. -9.9d0) THEN
            YLIN(I,1:50) = 0d0
         END IF
      END DO
      Y(1:50) = WN(1:50)
      DO I = 1, 50
         YLIN(1:50, I) = YLIN(1:50, I)*BARYN(I)*DMK
         YLIN(I, I) = YLIN(I, I) + DMK/DT
         YDOT(I) = (YDOT(I)*BARYN(I)+DWN(I)/DT)*DMK
      END DO
! Calculate neutron out rates for use later
! FIXME: we only use these when the nucleosynthesis has actually converged,
! so we might as well only bother to calculate these at that point. Should
! speed things up a bit.
      XA2(2) = WN(2)
      FRAC(:,:) = 0.0d0
      IF (XA2(2) == 0.0d0) RETURN      ! Don't bother if there are no neutrons
      DO I = 1, 50
         WW2(I) = XA2(I)/BARYN(I)
      END DO
      CALL NUCRAT2(AT)
      RATETOT = RDn + RnNA22p + RnNA22a + Rn26Mp + Rn26Ma + Rn26G
     :     + Rn26Gp + Rn26Ga
     :     + RnFe56 + RnFe57 + RnFe58 + RnCo59 + RnNi58 + RnNi59
     :     + RnNi60 + Rnp + RnHe3 + RnLi7 + RnC12 + RnC13 + RnC14
     :     + RnN14 + RnN15 + RnO16 + RnO18 + RnF19 + RnNe20 + RnNe21
     :     + RnNe22 + RnNa23 + RnMg24 + RnMg25 + RnMg26 + RnAl27
     :     + RnSi28 + RnSi29 + RnSi30 + RnP31 + RnS32 + RnS33 + RnS34
     :     + RnS33a + RnN14p + RnNi59p + RnNi59a + RnO17a + RnS32a
     :     + RnS34s + RnFe59 + RnFe60 + RnNi61s
      IF (RATETOT == 0d0) RATETOT = 1d0
      FRAC(1 ,K) = - RnNi61s
      FRAC(2 ,K) = 0d0
      FRAC(3 ,K) = - Rnp
      FRAC(4 ,K) = RnHe3
      FRAC(5 ,K) = RnLi7
      FRAC(6 ,K) = 0d0
      FRAC(7 ,K) = 0d0
      FRAC(8 ,K) = RnC13 - RnC12
      FRAC(9 ,K) = RnC14 - RnN14p - RnO17a - RnC13
      FRAC(10,K) = RnN15 -  RnC14 - RnN14
      FRAC(11,K) = - RnO16 + RnO17a
      FRAC(12,K) = RnO18
      FRAC(13,K) = - RnNA22a - RnO18 + RnF19
      FRAC(14,K) = RnNe21 - RnNe20
      FRAC(15,K) = - RnNe21 + RnNe22 - RnNA22p
      FRAC(16,K) = RnNA22p + RnNA22a
      FRAC(17,K) = - Rn26Ga - Rn26Ma - RnNe22 + RnNa23
      FRAC(18,K) = - RnNa23 + RnMg24
      FRAC(19,K) = RnMg25 - RnMg24
      FRAC(20,K) = - Rn26Mp - Rn26Gp + RnMg26 - RnMg25
      FRAC(21,K) = Rn26Mp + Rn26Ma
      FRAC(22,K) = Rn26Gp + Rn26Ga + Rn26G
      FRAC(23,K) = - RnMg26 + RnAl27 - Rn26G
      FRAC(24,K) = - RnAl27 + RnSi28
      FRAC(25,K) = - RnSi28 + RnSi29 - RnS32a
      FRAC(26,K) = - RnSi29 + RnSi30 - RnS33a
      FRAC(27,K) = - RnSi30 + RnP31
      FRAC(28,K) = - RnP31 + RnS32 + RnS32a
      FRAC(29,K) = - RnS32 + RnS33 + RnS33a
      FRAC(30,K) = - RnS33 + RnS34 + RnS34s
      FRAC(31,K) = RnFe56 - RnNi59a
      FRAC(32,K) = RnFe57 - RnFe56
      FRAC(33,K) = RnFe58 - RnFe57
      FRAC(34,K) = RnFe59 - RnFe58
      FRAC(35,K) = - RnFe59 + RnFe60
      FRAC(36,K) = RnCo59 - RnNi59p
      FRAC(37,K) = RnNi58
      FRAC(38,K) = RnNi59 - RnNi58 + RnNi59a + RnNi59p
      FRAC(39,K) = RnNi60 - RnNi59
      FRAC(40,K) = RnNi61s - RnNi60
      FRAC(41,K) = Rnp - RDn - Rn26Mp - Rn26Gp - RnNi59p - RnN14p
      FRAC(42,K) = - Rn26Ga - Rn26Ma - RnS32a - RnNa22a - RnS33a
     :                - RnNi59a! NEEDS RnHe4
      FRAC(43,K) = RnC12
      FRAC(44,K) = RnN14p
      FRAC(45,K) = RnO16
      FRAC(46,K) = - RnF19
      FRAC(1:48,K) = FRAC(1:48,K)/RATETOT
      RETURN
      END SUBROUTINE

      SUBROUTINE EQUNS2 ( JK, K1, K2, KL, KQ, KEQ )
      USE MESH
      IMPLICIT REAL*8 (A-H, L-Z)
!     The COMMON block INE contains
!     FN2(3,NFUNC), DFN2(3,NVAR,NFUNC), EQU(NEQ), DEQU(NVAR,3,NEQ)
!     The first block contains a copy of the elements in the COMMON block
!     INF2 from FUNCS2, specifically (in order)
!     YDOT(50), Y(50), SG, DMT, PADDING(SIZE)
!     The amount of padding to fill up the first element is 3*NFUNC - 102
      COMMON /INE   / XT(3, 50), X(3, 50), SG(3), DMT(3), WW(3*(NFUNC-102)),
     :  DXT(3, NVAR, NFUNC), EQU(NEQ), DEQU(NVAR, 3, NEQ)
C Next-to-surface boundary conditions, consistent with equns1
C FIXME: needs to account for thermohaline mixing and accretion
      !IF (JK /= K1) THEN
      DEQU(1:KEQ,1:3,1:KEQ) = 0.0d0
      EQU(1:KEQ) = 0.0d0
      !IF ( JK.LE.K1 + 1 ) THEN
      IF ( JK + KL == K1 + 1 ) THEN
         DO I1 = 1, KEQ
            DEQU(I1, 2, I1) = DEQU(I1, 2, I1) - 1.0
            DEQU(I1, 3, I1) = 1.0
            EQU(I1) = X(3,I1) - X(2,I1)
         END DO
      ELSE IF ( JK > K2 ) THEN
      !ELSE IF ( JK + KL == K2 + 1 ) THEN
C Next-to-central boundary conditions, consistent with equns1
         SG23 = 0.5D0*(SG(3) + SG(2))
         DO I1 = 1, KEQ
            DEQU(1:KEQ, 3, I1) = DXT(3 - KL, 1:KEQ, I1)
            DEQU(I1, 2, I1) = - SG23
            DEQU(I1, 3, I1) = DEQU(I1, 3, I1) + SG23
            EQU(I1) = SG23*(X(3, I1) - X(2, I1)) + XT(3 - KL, I1)
         END DO
      ELSE
C RJS: Interior points
C Needs checking, may not be conservative - ask Pierre!
         SG12 = 0.5D0*(SG(2) + SG(1)) - PSTV(KQ*DMT(2), 0.0D0)
         SG23 = 0.5D0*(SG(3) + SG(2)) - PSTV(-KQ*DMT(2), 0.0D0)
         DO I1 = 1, KEQ
            DO I2 = 1, KEQ
               DEQU(I2, 2, I1) = -DXT(2, I2, I1)
            END DO
            DEQU(I1, 1, I1) = SG12
            DEQU(I1, 2, I1) = DEQU(I1, 2, I1) - SG12 - SG23
            DEQU(I1, 3, I1) = SG23 
            EQU(I1) = SG23*(X(3, I1) - X(2, I1))
     :                - SG12*(X(2, I1) - X(1, I1)) - XT(2, I1)
         END DO
         !DO I1 = 1, KEQ
         !   DO I2 = 1, KEQ
         !      DO I3 = 1, 3
         !         if (dabs(dequ(i1,i3,i2))>1.0d0) print *, I1, I3, I2, DEQU(I1, I3, I2)
         !      END DO
         !   END DO
         !END DO
      END IF
      !END IF
      !print *, 'DEQU', JK
      !print *, DEQU(:,:,:)
      RETURN
      END SUBROUTINE


! ------------------------------------------------------------------------------
!  NEUTRON
!   Perform neutron captures
!   These are not done during the main part of the nucleosynthesis
!   calculation because the reactions are so fast compared to the other
!   reactions occurring at the same time that they are almost instantaneous
!   anyway. Trying to put them in the same matrix would cause
!   non-convergence.
! ------------------------------------------------------------------------------
      SUBROUTINE NEUTRON
      USE MESH
      USE NUCLEOSYNTHESIS
C RJS - deal with neutrons, basically assuming they're in equilibrium
      IMPLICIT REAL*8(A-H, O-Z)
      COMMON /NEUTRO/ FRAC(50,NM)
      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, KW(260)
      DIMENSION BARYN(50)
      DATA BARYN/62.0, 1.0, 2.0, 3.0, 7.0, 7.0, 11.0, 13.0, 14.0,  15.0,
     :     17.0, 18.0, 19.0, 21.0, 22.0, 22.0, 23.0, 24.0, 25.0, 26.0,
     :     26.0,26.0,27.0,28.0,29.0,30.0,31.0,32.0,33.0,34.0,56.0,57.0,
     :     58.0,59.0,60.0,59.0,58.0,59.0,60.0,61.0,1.0,4.0,12.0,14.0,16.0,
     :     20.0,1.0,1.0,1.0,1.0/ 

      DO K=1,KH
         DO I=1,50
            DHNUC(1,I,K) = DHNUC(1,I,K) - FRAC(I,K)*BARYN(I)*DHNUC(1,2,K)
         END DO
C reset neutron abundance to zero
      END DO
      DHNUC(1,2,1:KH) = 0d0
      END SUBROUTINE

! ------------------------------------------------------------------------------
!  NUCRAT2
!   Compute thermonuclear reaction rates for minor isotopes, for
!   nucleosynthesis (funcs2).
! ------------------------------------------------------------------------------
!  Input:
!     TL          - log temperature (in Kelvin)
!     /STAT2/     - Miscelaneous EoS quantities
!     /ABUND/     - Abundances of major isotopes used to compute STAT
!     /ABUND2/    - Abundances of minor isotopes
!  Output:
!     /RATES/     - Reaction rates
!     /NRATES/    - Neutron reaction rates
!     /DECAY/     - Decay rates
! ------------------------------------------------------------------------------
      SUBROUTINE NUCRAT2(TL)
      USE NUCLEOSYNTHESIS
      USE CONSTANTS
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON /STAT2 / AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA, 
     : SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, RPP, R33, R34, RBE, RBP,  
     : RPC, RPN, RPO, R3A, RAC, RAN, RAO, RANE, RCC, RCO, ROO, RGNE, 
     : RGMG, RCCG, RPNG, EX, ENX, WMU, DELTA, PHIE, EXT, FKT, FKR, PRANDTL
      COMMON /ABUND / XA(9), N(9), NE, NI, NZZ, AVM, NE1
      COMMON /NCDATA/ QRT(20), QNT(20), CZA(92), CZB(90), CZC(90),
     &                CZD(92), VZ(9)
      COMMON /ABUND2/ XA2(50),NNG, Nn, NN2,NN3,NNL7,NNB7,NN11,NN13,NN14,
     : NN15,NN17,NN18,NN19,NNE21,NNE22,NNA22,NNA23,NMG24,NMG25,NMG26,
     : N26M,N26G,NAL27,NSI28,NSI29,NSI30,NP31,NS32,NS33,NS34,
     : NFE56,NFE57,NFE58,NFE59,NFE60,NCO59,NNI58,NNI59,NNI60,NNI61,NN1,
     : NN4,NN12,NNN14,NN16,NN20,W(4)
      COMMON /RATES / RRT2(92)
      COMMON /NRATES/ NRATE(45)
      COMMON /DECAY / RDAL26G,RDNA22,RD26M,RDFe59,RDFe60,RDNi59,RDn,RDC14
      DATA CSA, CSB, CSC, CSD, CXD /0.624, 0.316, 0.460, 0.38, 0.86/
      CBRT(VX) = DEXP(DLOG(VX)/3.0D0)
* RHB is 'baryon density': 1 amu * number of baryons per cm3
      RHB = RHO/AVM
* Electron screening theory from Graboske, DeWitt, Grossman & Cooper (1973),
* for strong (ZA, ZB, ZC) are intermediate screening (ZD). The reaction
* dependent charge parameters are stored in CZA ... CZD.
      WC = DOT_PRODUCT(N(1:9), VZ(1:9))
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
      DO J = 1, NUC_NUM_RATES
         RN = 0.0D0
         TT = 50.0*(TF - 6.0) + 1.0
         IF ( TT .GE. 1.0D0 ) THEN
            IT = MAX(1, MIN(199, INT(TT)))
            TT = TT - IT
            TU = 1.0D0 - TT
            RR = TU*NUCSYN_CRT(IT, J) + TT*NUCSYN_CRT(IT+1, J)
            IF ( RR .GE. -50.0D0 ) THEN
               SCRN = ZD*CZD(J)
               STRN = ZA*CZA(J) + ZB*CZB(J)
               DSTR = ZC*CZC(J)
               IF (DSTR .LT. 0.29*STRN) SCRN = MIN(SCRN, STRN - DSTR)
               RN = EXP(CLN*RR + SCRN)
            END IF
         END IF
         RRT2(J) = RN
      END DO
C Sort out rates
      RPP = RRT2(1)
      R33 = RRT2(2)
      R34 = RRT2(3)
      RBE = RRT2(4)
      RBP = RRT2(5)
      RPC = RRT2(6)
      RPN = RRT2(7)
      RPO = RRT2(8)
      R3A = RRT2(9)
      RAC = RRT2(10)
      RAN = RRT2(11)
      RAO = RRT2(12)
      RANE = RRT2(13)
      RCC = RRT2(14)
      RCO = RRT2(15)
      ROO = RRT2(16)
      RGNE = RRT2(17)
      RGMG = RRT2(18)
      RCCG = RRT2(19)
      RPNG = RRT2(20)
* Multiply with density and abundances to get rates per baryon per second,
* note that abundances of He3 and Be7 are not needed in equilibrium
      RPP = RHB*NN1*NN1*RPP/2.0
      R33 = RHB*NN3*NN3*R33/2.0
      R34 = RHB*NN3*NN4*R34
      RBE = RHB*NE*RBE
      RBP = RHB*NN1*RBP
      RPC = RHB*NN1*NN12*RPC
      RPN = RHB*NN1*NNN14*RPN
      RPO = RHB*NN1*NN16*RPO
      R3A = RHB*RHB*NN4*NN4*NN4*R3A/6.0
      RAC = RHB*NN4*NN12*RAC
      RAN = RHB*NN4*NNN14*RAN
      RAO = RHB*NN4*NN16*RAO
      RANE = RHB*NN4*NN20*RANE
      RCC = RHB*NN12*NN12*RCC/2.0
      RCO = RHB*NN12*NN16*RCO
      ROO = RHB*NN16*NN16*ROO/2.0
      RGNE = NN20*RGNE
      RGMG = NMG24*RGMG
* Branching of pN and CC reactions
      FPNG = 8.0D-4
      RPNA = (1.0 - FPNG)*RPN
      RPNG = FPNG*RPN
      RPN = RPNA
      FCCG = RCCG
      RCCA = (1.0 - FCCG)*RCC
      RCCG = FCCG*RCC
      RCC = RCCA
C Put rates back to RRT
      RRT2(1) = RPP
      RRT2(2) = R33
      RRT2(3) = R34
      RRT2(4) = RBE
      RRT2(5) = RBP
      RRT2(6) = RPC
      RRT2(7) = RPN
      RRT2(8) = RPO
      RRT2(9) = R3A
      RRT2(10) = RAC
      RRT2(11) = RAN
      RRT2(12) = RAO
      RRT2(13) = RANE
      RRT2(14) = RCC
      RRT2(15) = RCO
      RRT2(16) = ROO
      RRT2(17) = RGNE
      RRT2(18) = RGMG
      RRT2(19) = RCCG
      RRT2(20) = RPNG
C Minor variable reaction rates - 3/9/03 RJS
      RRT2(21) = RHB*NN1*NN2*RRT2(21)
      RRT2(22) = RHB*NN1*NNB7*RRT2(22)
      RRT2(23) = RHB*NN1*NNL7*RRT2(23)
      RRT2(24) = RHB*NE*NNB7*RRT2(24)
      RRT2(25) = RHB*NN1*NN13*RRT2(25)
      RRT2(26) = RHB*NN1*NN15*RRT2(26)
      RRT2(27) = RHB*NN1*NN17*RRT2(27) 
      RRT2(28) = RHB*NN4*NN13*RRT2(28)
      RRT2(29) = RHB*NN4*NN15*RRT2(29)
      RRT2(30) = RHB*NN4*NN17*RRT2(30)
      RRT2(31) = RHB*NN4*NNL7*RRT2(31)
      RRT2(32) = RHB*NN1*NN11*RRT2(32)
      RRT2(91) = RHB*NN1*NN17*RRT2(91) 
      RRT2(92) = RHB*NN1*NN20*RRT2(92) 
C C14 reactions
      RRT2(33) = RHB*NN1*NN14*RRT2(33)
      RRT2(34) = RHB*NN4*NN14*RRT2(34)
C O18 reactions
      RRT2(35) = RHB*NN1*NN18*RRT2(35)
      RRT2(36) = RHB*NN1*NN18*RRT2(36)
      RRT2(37) = RHB*NN4*NN18*RRT2(37)
      RRT2(38) = RHB*NN4*NN18*RRT2(38)
C F19 reactions
      RRT2(39) = RHB*NN1*NN19*RRT2(39)
      RRT2(40) = RHB*NN1*NN19*RRT2(40)
      RRT2(41) = RHB*NN4*NN19*RRT2(41)
C Ne21
      RRT2(42) = RHB*NN1*NNE21*RRT2(42)
      RRT2(43) = RHB*NN4*NNE21*RRT2(43)
      RRT2(44) = RHB*NN4*NNE21*RRT2(44)
C Ne22
      RRT2(45) = RHB*NN1*NNE22*RRT2(45)
      RRT2(46) = RHB*NN4*NNE22*RRT2(46)
      RRT2(47) = RHB*NN4*NNE22*RRT2(47)
C Na22
      RRT2(48) = RHB*Nn*NNA22*RRT2(48)
      RRT2(49) = RHB*Nn*NNA22*RRT2(49)
      RRT2(50) = RHB*NN1*NNA22*RRT2(50)
C Na23
      RRT2(51) = RHB*NN1*NNA23*RRT2(51)
      RRT2(52) = RHB*NN1*NNA23*RRT2(52)
      RRT2(53) = RHB*NN1*NNA23*RRT2(53)
C Mg24
      RRT2(54) = RHB*NN1*NMG24*RRT2(54)
      RRT2(55) = RHB*NN4*NMG24*RRT2(55)
C Mg25
      RRT2(56) = RHB*NN1*NMG25*RRT2(56)
      RRT2(57) = RHB*NN1*NMG25*RRT2(57)
      RRT2(58) = RHB*NN1*NMG25*RRT2(58)
      RRT2(59) = RHB*NN4*NMG25*RRT2(59)
      RRT2(60) = RHB*NN4*NMG25*RRT2(60)
      RRT2(61) = RHB*NN4*NMG25*RRT2(61)
C Mg26
      RRT2(62) = RHB*NN1*NMG26*RRT2(62)
      RRT2(63) = RHB*NN1*NMG26*RRT2(63)
      RRT2(64) = RHB*NN1*NMG26*RRT2(64)
      RRT2(65) = RHB*NN1*NMG26*RRT2(65)
      RRT2(66) = RHB*NN4*NMG26*RRT2(66)
      RRT2(67) = RHB*NN4*NMG26*RRT2(67)
      RRT2(68) = RHB*NN4*NMG26*RRT2(68)
C Blank beause I put in a reaction that didn't exist!
      RRT2(69) = 0d0
      RRT2(70) = 0d0
      RRT2(71) = 0d0
      RRT2(72) = 0d0
C Al26M
      RRT2(73) = RHB*N26M*RRT2(73)
      RRT2(74) = RHB*Nn*N26M*RRT2(74)
      RRT2(75) = RHB*Nn*N26M*RRT2(75)
      RRT2(76) = RHB*NN1*N26M*RRT2(76)
C Al26G
      RRT2(77) = RHB*N26G*RRT2(77)
      RRT2(78) = RHB*Nn*N26G*RRT2(78)
      RRT2(79) = RHB*Nn*N26G*RRT2(79)
      RRT2(80) = RHB*NN1*N26G*RRT2(80)
C Al27
      RRT2(81) = RHB*NN1*NAL27*RRT2(81)
      RRT2(82) = RHB*NN1*NAL27*RRT2(82)
      RRT2(83) = RHB*NN4*NAL27*RRT2(83)
C Na23(a,n)Al26TGM
      RRT2(84) = RHB*NN4*NNA23*RRT2(84)
      RRT2(85) = RHB*NN4*NNA23*RRT2(85)
      RRT2(86) = RHB*NN4*NNA23*RRT2(86)
C Si reactions
      RRT2(87) = RHB*NN1*NSI28*RRT2(87)
      RRT2(88) = RHB*NN1*NSI29*RRT2(88)
      RRT2(89) = RHB*NN1*NSI30*RRT2(89)
C N15
      RRT2(90) = RHB*NN1*NN15*RRT2(90)

C Unstable particle decays
C Al26G t_0.5 = 0.72 Myr
      CLN2 = 0.69314718
      RDAL26G = (RHB/26.0)*N26G*CLN2/2.27d13
C C14 t_0.5 = 5730 yr
      RDC14 = (RHB/14.0)*NN14*CLN2/1.81d11
C Na22 t_0.5 = 2.6 yr
      RDNA22 = (RHB/22.0)*NNA22*CLN2/8.199d7
C Al26M t_0.5 = 6 s Should just shortcut the network...
      RD26M = (RHB/26.0)*N26M*CLN2/6.0
C Fe59 t_0.5 = 44.6 d
      RDFe59 = (RHB/59.0)*NFE59*CLN2/3.85d6
C Fe60 t_0.5 = 1.5 Myr
      RDFe60 = (RHB/60.0)*NFE60*CLN2/4.73d13
C Ni59 t_0.5 = 0.075 Myr
      RDNi59 = (RHB/59.0)*NNi59*CLN2/2.365d12
C Free n t_0.5 = 10.3 min
      RDn = (RHB/1.0)*Nn*CLN2/6.18d2
C (n,g) reactions
      IF (Nn == 0.0D0) THEN
         NRATE(:) = 0.0d0
         RETURN
      END IF
      DO J = 1, 45
         RN = 0.0D0
         TT = 50.0*(TF - 6.0) + 1.0
         IF ( TT .GE. 1.0D0 ) THEN
            IT = MAX(1, MIN(199, INT(TT)))
            TT = TT - IT
            TU = 1.0D0 - TT
            RR = TU*NUCSYN_NRT(IT, J) + TT*NUCSYN_NRT(IT+1, J)
            IF (RR.GE.-50.0d0) THEN
               RN = EXP(CLN*(RR + 20.0D0))*1.0D-20
            END IF
         END IF
         NRATE(J) = RN
      END DO
      NRATE(1) = RHB*Nn*NFE56*NRATE(1)
      NRATE(2) = RHB*Nn*NFE57*NRATE(2)
      NRATE(3) = RHB*Nn*NFE58*NRATE(3)
      NRATE(4) = RHB*Nn*NCO59*NRATE(4)
      NRATE(5) = RHB*Nn*NNI58*NRATE(5)
      NRATE(6) = RHB*Nn*NNI59*NRATE(6)
      NRATE(7) = RHB*Nn*NNI60*NRATE(7)
      NRATE(8) = RHB*Nn*NN1*NRATE(8)
      NRATE(9) = RHB*Nn*NN3*NRATE(9)
      NRATE(10) = RHB*Nn*NNL7*NRATE(10)
      NRATE(11) = RHB*Nn*NN12*NRATE(11)
      NRATE(12) = RHB*Nn*NN13*NRATE(12)
      NRATE(13) = RHB*Nn*NN14*NRATE(13)
      NRATE(14) = RHB*Nn*NNN14*NRATE(14)
      NRATE(15) = RHB*Nn*NN15*NRATE(15)
      NRATE(16) = RHB*Nn*NN16*NRATE(16)
      NRATE(17) = RHB*Nn*NN18*NRATE(17)
      NRATE(18) = RHB*Nn*NN19*NRATE(18)
      NRATE(19) = RHB*Nn*NN20*NRATE(19)
      NRATE(20) = RHB*Nn*NNE21*NRATE(20)
      NRATE(21) = RHB*Nn*NNE22*NRATE(21)
      NRATE(22) = RHB*Nn*NNA23*NRATE(22)
      NRATE(23) = RHB*Nn*NMG24*NRATE(23)
      NRATE(24) = RHB*Nn*NMG25*NRATE(24)
      NRATE(25) = RHB*Nn*NMG26*NRATE(25)
      NRATE(26) = RHB*Nn*NAL27*NRATE(26)
      NRATE(27) = RHB*Nn*NSI28*NRATE(27)
      NRATE(28) = RHB*Nn*NSI29*NRATE(28)
      NRATE(29) = RHB*Nn*NSI30*NRATE(29)
      NRATE(30) = RHB*Nn*NP31*NRATE(30)
      NRATE(31) = RHB*Nn*NS32*NRATE(31)
      NRATE(32) = RHB*Nn*NS33*NRATE(32)
      NRATE(33) = RHB*Nn*NS34*NRATE(33)
      NRATE(34) = RHB*Nn*N26G*NRATE(34)
      NRATE(35) = RHB*Nn*NS33*NRATE(35)
      NRATE(36) = RHB*Nn*NNN14*NRATE(36)
      NRATE(37) = RHB*Nn*NNI59*NRATE(37)
      NRATE(38) = RHB*Nn*NNI59*NRATE(38)
      NRATE(39) = RHB*Nn*NN17*NRATE(39)
      NRATE(40) = RHB*Nn*N26G*NRATE(40)
      NRATE(41) = RHB*Nn*NS32*NRATE(41)
      NRATE(42) = RHB*Nn*NFE59*NRATE(42)
      NRATE(43) = RHB*Nn*NFE60*NRATE(43)
      NRATE(44) = RHB*Nn*NS34*NRATE(44)
      NRATE(45) = RHB*Nn*NNI61*NRATE(45)
      RETURN
      END SUBROUTINE

! ------------------------------------------------------------------------------
!  CHECKS2
!   Check nucleosynthesis variables for small an negative values
!   Modifies HNUC and DHNUC
! ------------------------------------------------------------------------------
      SUBROUTINE CHECKS2(JSTAR)
      USE NUCLEOSYNTHESIS
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JSTAR
      DOUBLE PRECISION :: H, DH, EPS, DEL, DH0
      INTEGER :: KH, KTW, KW
      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, KW(260)
      INTEGER :: IK, IJ
! Composition variables - should not be negative!
      DO IK = 1, KH_NUC
         DO IJ = 1, NVAR_NUC
            IF ( HNUC(JSTAR, IJ, IK) + DHNUC(JSTAR, IJ, IK) < 0.0 ) THEN
               DHNUC(JSTAR, IJ, IK) = -HNUC(JSTAR, IJ, IK)
               HNUC(JSTAR, IJ, IK) = 0.0D0
            END IF
            !IF ( HNUC(JSTAR, IJ, IK) < 1.0d-20) HNUC(JSTAR, IJ, IK) = 0.0d0
         END DO
! Set H abundance in nucleosynthesis code to 0 if structure says it is 0
         !IF (H(5 + (JSTAR-1)*24, IK) == 0.0d0) HNUC(JSTAR, 41, IK) = 0.0d0
      END DO
      END SUBROUTINE

! ------------------------------------------------------------------------------
!  UPDATE2
!   Update nucleosynthesis variables
! ------------------------------------------------------------------------------
      SUBROUTINE UPDATE2
      USE NUCLEOSYNTHESIS
      IMPLICIT NONE
      IF (ALLOCATED(HNUC)) THEN
         HNUCPR(:,:,:) = HNUC(:,:,:)
         HNUC(:,:,:) = HNUC(:,:,:) + DHNUC(:,:,:)
      END IF
      END SUBROUTINE


! ------------------------------------------------------------------------------
!  BACKUP2
!   Backup nucleosynthesis variables
! ------------------------------------------------------------------------------
      SUBROUTINE BACKUP2
      USE NUCLEOSYNTHESIS
      IMPLICIT NONE
      IF (ALLOCATED(HNUC)) THEN
         HNUC(:,:,:) = HNUCPR(:,:,:)
         DHNUC(:,:,:) = 0.0d0
      END IF
      END SUBROUTINE

