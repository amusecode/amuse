! FUNCS1:
!  Evaluate the functions of independent variables needed to set up the
!  difference equations.
!  JK is the number of the meshpoint to calculate values for
!  JI determines from where the function is called, or which variable was changed
!      >0: The number of the variable that has been varied for this call
!       0: No variable changed
!      -1: PRINTB
!      -2: REMESH
!  The last two options affect what is calculated by the code, the others
!   are used for caching the output from the equation-of-state.
!  This routine first copies independent variables from X in INF to NAMEIN
!  The calculated quantities are in the common block NAMEOUT and then
!   transfered to Y in INF.
      SUBROUTINE FUNCS1 ( JK, JI )
      USE MESH
      USE MESH_ENC
      USE EXTRA_ELEMENTS
      USE FUDGE_CONTROL
      USE CONSTANTS
      USE SETTINGS
      USE MASSLOSS      
      USE PLOTVARIABLES
      USE EXPLICIT_FUNCTIONS
      USE NUCLEOSYNTHESIS
!      IMPLICIT REAL*8 (A-H, L-Z)  
      IMPLICIT NONE
C DECLARING SUBROUTINE VARIABLES
      INTEGER, INTENT(IN) :: JK, JI
C DECLARING COMMON BLCOK VARIABLES EXPLICITLY
!    /     /
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, KW(260)
!     /QUERY /
      DOUBLE PRECISION ::  ML, QL, XL, UC(21)
      INTEGER :: JMOD, JB, JNN, JTER, JOC, JKH
!     /TVBLES/
      DOUBLE PRECISION :: DT, ZQ, HSPN(2), SRLF(2), ZET(2), XIT(2), AGE,
     $     BM,MC(2), OM, PER, SM, ENC, TC(2), TFR, T0, M0, MTA, OM0,
     $     OMTA, A0,ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS, SE, TN(2),
     $     WMH, WMHE,MH, MHE, MCO, VMG, BE(2), LH, LHE, LC, LNU, LTH,
     $     MCB(8),MSB(6), RCB(8), TCT(8), HR(81), PPR(81)
      INTEGER :: JHOLD, JM2, JM1
!     /VBLES /
      DOUBLE PRECISION :: LOLEDD, DG, EG, GRADT, ETH, EGR, R, QQ, QMU, 
     : WL, WCV, HP, TW, PHIMU, GMR, SEP, M3,PX(NPX), SX(NPX,NM+1),QA(NM)
!     /STORE /
      DOUBLE PRECISION ::  HPR(NVAR,NM), MS(60025)
!     /INF   /
      DOUBLE PRECISION :: X(NVAR), DX(NVAR), Y(NFUNC), Z(NVAR*NFUNC)
!     /NAMEIN/  SIZE(NAMEIN) = 3*NFUNC
!     This is an abused commonblock, containing different things in
!     different subroutines: in [FUNCS1] it contains the INdepent
!     variables for the current star, the binary orbit, and their
!     changes, after executing CALL NAMES1(JSTAR,1). In [EQUNS1] it
!     contains the Dependent variables or 'functions', at the current,
!     previous and anteprevious meshpoint instead.  We need to fill up
!     the unused remaining entries to make sure that the size of the
!     commonblock is the same everywhere. -SdM
      DOUBLE PRECISION :: AF, AT, VX16, M, VX1, VQK, AR, L, VX4, VX12,
     $     VX20,AI, VOMEGA, PHI, PHIS, VX14,          !16 stellar vars
     $     OA, ECC, XI, MB, VX24, VMENC, VMEA, VMET,  ! 8 binary vars
     $     DAF, DAT, DX16, DM, DX1, DVQK, DAR, DL, DX4, DX12,
     $     DX20,DAI, DOMEGA, DPHI, DPHIS, DX14,       !16 stellar var. changes
     $     DOA, DE, DXI, DMB,DX24, DMENC, DMEA, DMET, ! 8 bin var. changes
     $     FILLUP_COMMONBLOCK_NAMEIN(3*NFUNC-(NVSTAR+NVBIN)),  ! Padding
     $     STARVAR, BINVAR, DSTARVAR, DBINVAR                  ! Padding
!     /NAMOUT/   ! SIZE(NAMEOUT) = NFUNC
      DOUBLE PRECISION ::  BCP, BCT, VP, VPK, VR, VRK, VT, VTK, VL, LK,
     $     LQ, MT, VM, VMK, SG, WT, X1, X1T, X16, X16T, X4, X4T, X12,
     $     X12T, X20, X20T, BCM, VI, VIK, VPHI, PHIK, BCF, BCS, BCPH,
     $     X14, X14T, AVMU, SGTH, OMEGA, OMEGAT, SGAM, SI, BCA, BCE, XIM
     $     , XIK, DLRK, BCMB, X24, X24T, MENC, MENCT, MEA, MEAT, MET,
     $     METT, FILLUP_COMMONBLOCK_NAMOUT(2 + NSFSTAR + NXFUNC  )
!     give nameout same size as in other subroutines:   -SdM
!     NFUNC = NSFSTAR1 (42) +   ! nicely filled
!             NSFBIN   (16) +   ! of which only 14 are really used,  needs fillup
!             NSFSTAR2 (42) +   ! not used but filled with  FILLUP...SFSTAR2
!             NXFUNC   (60) +   ! not used but filled with  FILLUP...XFUNC
!     /XIN1  /
      DOUBLE PRECISION :: XVS(NSVAR+1          :NSVAR +  NXVSTAR), 
     :                    XVB(NSVAR+2*NXVSTAR+1:NSVAR +2*NXVSTAR+NXVBIN),
     :                   DXVS(NSVAR+1          :NSVAR +  NXVSTAR), 
     :                   DXVB(NSVAR+2*NXVSTAR+1:NSVAR +2*NXVSTAR+NXVBIN)
!     /XOUT1 /
      DOUBLE PRECISION :: XFS(NSFUNC+1          :NSFUNC +  NXFSTAR), 
     :                    XFB(NSFUNC+2*NXFSTAR+1:NSFUNC +2*NXFSTAR+NXFBIN)
!     /STAT2 
      DOUBLE PRECISION ::  AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA
     $     ,SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, RPP, R33, R34, RBE,
     $     RBP,RPC, RPN, RPO, R3A, RAC, RAN, RAO, RANE, RCC, RCO, ROO,
     $     RGNE,RGMG, RCCG, RPNG, EX, ENX, WMU, DELTA, PHIE, EXT, FKT,
     &     FKR, PRANDTL
!     /ABUND /
      DOUBLE PRECISION :: XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, 
     &     NA(9), NEO, NIO, NZZ, AVMA, NE
!     /ACCRET /
      DOUBLE PRECISION :: X1AC, X4AC, X12AC, X14AC, X16AC, X20AC, X24AC
      DOUBLE PRECISION :: XAC(7, 2)
!     /ATDATA/ 
      DOUBLE PRECISION :: CH2(4), CHI(26,9), COM(27), CAN(9), CBN(9)
      INTEGER :: KZN(9)
!      /LSMOOTH/
      DOUBLE PRECISION :: LK_PREV(NM), LQ_PREV(NM), ENUC_PREV(NM)
! Various stuff, needed for gravitational settling
      INTEGER :: NREF = 1                    ! Reference abundance (H)
      DOUBLE PRECISION :: D12, D0, SGAT, SGATF, GI, APR, WPI(9), Ddiff(9), XA(9)
      DOUBLE PRECISION :: LAMBDA_D2, BETA_D2, WPP, ZZ, KF_NUM, KF_DENOM, KT, A12

C DECLARING COMMON BLOCKS
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW
      COMMON /QUERY / ML, QL, XL, UC, JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /TVBLES/ DT, ZQ, HSPN, SRLF, ZET, XIT, AGE, BM, 
     : MC, OM, PER, SM, ENC, TC, TFR, T0, M0, MTA, OM0, OMTA, A0, 
     : ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS, SE, TN, WMH, WMHE, 
     : MH, MHE, MCO, VMG, BE, LH, LHE, LC, LNU, LTH, MCB, 
     : MSB, RCB, TCT, HR, PPR, JHOLD, JM2, JM1
      COMMON /VBLES / LOLEDD, DG, EG, GRADT, ETH, EGR, R, QQ, QMU, 
     : WL, WCV, HP, TW, PHIMU, GMR, SEP, M3,PX, SX,QA
      COMMON /STORE / HPR, MS
c The first 40 vbles of /INF/ are the guessed quantities at each meshpt,
c the next 40 are their changes since previous t-step. Only 11 or 13 or 19
c vbles are used currently (34 if KTW = 2). The next 100 vbles are 
c functions of the first 80, to be used in constructing the diff. equns.
c In NAMES1, X, DX of /INF/ are equivalenced to /NAMEIN/, Y to /NAMOUT/
! NAME1 first lists the independent variables and then the changes in the
!  last time step.
      COMMON /INF   / X, DX, Y, Z
      COMMON /NAMEIN/ AF, AT, VX16, M, VX1, VQK, AR, L, VX4, VX12, VX20, 
     : AI, VOMEGA, PHI, PHIS, VX14, VX24,          ! 17 entries for *1
     : STARVAR(NVSTAR - 17),                       ! Padding
     : OA, ECC, XI, MB,                            ! 4 entries for bin
     : BINVAR(NVBIN - 4),                          ! Padding
     : DAF, DAT, DX16, DM, DX1, DVQK, DAR, DL, DX4, DX12, DX20, 
     : DAI, DOMEGA, DPHI, DPHIS, DX14, DX24,       ! 17 entries for *1
     : DSTARVAR(NVSTAR - 17),                      ! Padding
     : DOA, DE, DXI, DMB,                          ! 4 entries for bin
     : DBINVAR(NVBIN - 4),                         ! Padding
     : FILLUP_COMMONBLOCK_NAMEIN                   ! fill up remaining space
      COMMON /NAMOUT/ BCP, BCT, VP, VPK, VR, VRK, VT,
     : VTK, VL, LK, LQ, MT, VM, VMK, SG, WT, X1, X1T, X16, X16T, X4, 
     : X4T, X12, X12T, X20, X20T, BCM, VI, VIK, VPHI, PHIK, BCF, BCS, 
     : BCPH, X14, X14T, AVMU, SGTH, OMEGA, OMEGAT, SGAM, SI,  ! Star 1 (42 values)
     : BCA, BCE, XIM, XIK, DLRK, BCMB, X24, X24T, MENC, MENCT, MEA, 
     : MEAT, MET, METT, FILLUP_COMMONBLOCK_NAMOUT
      ! Extra variables: for star (S) and binary (B), with changes
      COMMON /XIN1  / XVS, XVB, DXVS, DXVB
      ! Extra output variables: for star (S) and binary (B)
      COMMON /XOUT1 / XFS, XFB
      COMMON /STAT2 / AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA, 
     : SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, RPP, R33, R34, RBE, RBP,  
     : RPC, RPN, RPO, R3A, RAC, RAN, RAO, RANE, RCC, RCO, ROO, RGNE, 
     : RGMG, RCCG, RPNG, EX, ENX, WMU, DELTA, PHIE, EXT, FKT, FKR, PRANDTL
      COMMON /ABUND / XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, NA, NEO, NIO, NZZ, AVMA, NE
      COMMON /ACCRET/ X1AC, X4AC, X12AC, X14AC, X16AC, X20AC, X24AC, XAC
      COMMON /ATDATA/ CH2, CHI, COM, CAN, CBN, KZN
      COMMON /LSMOOTH/ LK_PREV, LQ_PREV, ENUC_PREV

C DECLARING LOCAL FUNCTIONS
      DOUBLE PRECISION :: CBRT, RLOBE
      DOUBLE PRECISION :: VX

C DECLARING EXTERNAL FUNCTIONS
      DOUBLE PRECISION :: STEP, PSTV, CALC_MDOT_JVINK

C DECLARING LOCAL VARIABLES
      INTEGER :: JSTAR, IKK, J, JSTAR_ACCRETOR

      DOUBLE PRECISION :: AVMU_PREV, AVMU_NEW, R2, R3, WTH, ENUC,
     $     DENC, TARGET_DAF, TARGET_DAT, DS, MOR3, R0, LOR3,  R2MU, FP,
     $     FT, GOR, APMU, GRADR, S2, WC1, MU, LOM, WC2, WC3, WC4, GTA,
     $     ATMU, VPP, VTT, XMC, VMF, VMMU, VRR, MUK, MK, B, WOS,
     $     CON, FAC, F_WL, VML, DF0, SGF, DEG, UG, W, BETA, FAC1, FAC2,
     $     V2, SGTH_NUMER, SGTH_DENOM, DTH, SGTH2, REC, VES, VMU, APP,
     $     MMK, ABSL, AVM, XH_SURFACE, XHE_SURFACE, Z_EFF, RHT, RAD, W2,
     $     TET, ROT, MDTDDW, ALFAC, MDTL, MDTSW, MDOT_JVINK, L_PREV, T_PREV,
     $     MDOT_KUDRITZKI, MDOT_WR, MDOT_CMI, MDOT, S_X, S_HE, S_C, S_O,
     &     F_PHOTON_TIRING, FDT, M2, E2, EF, WE25, RL, RGR, ETGR, OATGR,
     &     RLF, XI1, WW, GAM, WE5, WE65, WE2, WE3, WEA, WEB, DMT0, CIF, ALPHA

      DOUBLE PRECISION :: DG_LEDOUX, DSC, GRADMU, TKH, NU
      DOUBLE PRECISION :: APER, DAP

C For determining mass transfer rates
      DOUBLE PRECISION :: RR, RPR, GMRR, GMRP(2)

      DOUBLE PRECISION :: RA2S(2), AIS(2), OATMB(2)
      DOUBLE PRECISION :: ET(2), MP(2), PA(2), BR(2), ZEP(2), W1P(2),
     $     W2P(2), W3P(2),W4P(2), SPP(2), DSP(2), SUP(2), DMP(2),BRP(2),
     $     ENP(2), MKP(2), PHIP(2), PHISP(2), PAP(2), SDP(2)
      DOUBLE PRECISION :: MDOTEDD(2), LEDD(2), FMDOTLOST(2), PHIRATIO(2)

C     related to spinup
      DOUBLE PRECISION :: Rmin_stream, R_spinup
      DOUBLE PRECISION :: WOWCRIT, WOWCRITP(1:2) 

C  DEFINING INLINE STATEMENT FUNCTIONS
      CBRT(VX) = DABS(VX)**C3RD
      RLOBE(VX) = 0.49D0*VX*VX/(0.6D0*VX*VX + DLOG(1.0D0 + VX))
      JSTAR = 1
! NAMES1 is a lengthy equivalence statement, between INF, NAMEIN and NAMOUT
! After this call, the variables in INF are stored in NAMEIN
 1    CALL NAMES1 ( JSTAR, 1 )
c M is the mass variable
      MT = DM/DT
      DMP(JSTAR) = DM

!     Compute the rotational period from the rotational frequency
      OMEGA = VOMEGA
      APER = 2.0*CPI/(DABS(OMEGA) * CSDAY)                        ! CGS
      DAP = -APER/DABS(OMEGA) * DOMEGA

c Calculate the mean-molecular weight (for thermohaline mixing)
c Use the values from the previous timestep instead of the current one; this
c seems to be more numerically stable
c FIXME: now that we have the value of the molecular weight gradient (for the
c purpose of semiconvection) we should rework the TH mixing prescription.
      XH =  H(5 + (JSTAR-1)*24, JK)
      XHE = H(9 + (JSTAR-1)*24, JK)
      XC =  H(10 + (JSTAR-1)*24, JK)
      XN =  H(16 + (JSTAR-1)*24, JK)
      XO =  H(3 + (JSTAR-1)*24, JK)
      XNE = H(11 + (JSTAR-1)*24, JK)
      XMG = MAX(0.0D0, 1.0D0 - XH - XHE - XC - XN - XO - XNE - XSI - XFE)
      AVMU_PREV = 1.0/(0.5 + 1.5*XH + (42.0*XHE + 14.0*XC + 12.0*XN + 10.5*XO +
     &     8.4*XNE + 7.0*XMG + 6.0*XSI - 3.0*XFE)/168.0)
C set up the composition variables (sanely)
      XH =  MAX(0.0D0, MIN(VX1,  1.0D0))
      XHE = MAX(0.0D0, MIN(VX4,  1.0D0))
      XC =  MAX(0.0D0, MIN(VX12, 1.0D0))
      XN =  MAX(0.0D0, MIN(VX14, 1.0D0))
      XO =  MAX(0.0D0, MIN(VX16, 1.0D0))
      XNE = MAX(0.0D0, MIN(VX20, 1.0D0))

      IF (USE_MG24_EQN) THEN
         XMG = MAX(0.0D0, MIN(VX24, 1.0D0))
      ELSE
         XMG = 1.0D0 - XH - XHE - XC - XN - XO - XNE - XSI - XFE
         IF ( XMG < 0.0D0 ) XMG = 0.0D0
      END IF

      AVMU_NEW = 1.0/(0.5 + 1.5*XH + (42.0*XHE + 14.0*XC + 12.0*XN + 10.5*XO +
     &     8.4*XNE + 7.0*XMG + 6.0*XSI - 3.0*XFE)/168.0)

      AVMU = AVMU_SMOOTH*AVMU_NEW + (1.0-AVMU_SMOOTH)*AVMU_PREV
      
c Equation of state
      CALL STATEL ( JI, AF, AT, JSTAR )
! Store composition variables in an array as well, for convenience
      XA(1) = XH
      XA(2) = XHE
      XA(3) = XC
      XA(4) = XN
      XA(5) = XO
      XA(6) = XNE
      XA(7) = XMG
      XA(8) = XSI
      XA(9) = XFE
      ENP(JSTAR) = U + P/RHO
      R2 = DEXP(2.0D0*AR) - CT(8)
      R = DSQRT(DABS(R2))
      R3 = R*ABS(R2)
      WTH = 2.0D0*KTH/(KTH**4 + 1.0D0) * LUMINOSITY_FUDGE
      IF ( JMOD == 0 ) WTH = 0.0D0
      ENUC = LLUMI_SMOOTH*EX + (1.0D0-LLUMI_SMOOTH)*ENUC_PREV(JK)
      DENC = 0.0D0
      IF (USEMENC .AND. M<TH(4,1)) THEN
         !LOM = ENUC + EN + ENX + MENC*MEA - WTH*T*(SF*DAF + ST*DAT)/DT
         TARGET_DAF = GET_IPH(M, 1)-AF
         TARGET_DAT = GET_IPH(M, 2)-AT
         DS = SF*TARGET_DAF + ST*TARGET_DAT
         DENC = IMPOSE_ENTROPY_FACTOR*(T*DS/DT)
         LOM = ENUC + EN + ENX + DENC * LUMINOSITY_FUDGE !LOM = dL/dM
      ELSE
         LOM = ENUC + EN + ENX + ENC - WTH*T*(SF*DAF + ST*DAT)/DT
      ENDIF
      LOM = LOM + ENC_PARACHUTE

      IF ( JK == KH ) THEN
! Approximations for the centre (M=0); the normal expressions become singular
         MOR3 = CPI4*RHO*C3RD
         LOR3 = LOM*MOR3
         R0 = 1.0D-10
         M3 = CBRT(MOR3 * R0**3)
         MP(JSTAR) = M3**3
      ELSE
! Values away from the centre
         MOR3 = ABS(M)/R3
         LOR3 = L/R3 !note the difference with LOM=dL/dm and not L/M
         M3 = CBRT(ABS(M))
         R0 = R
         MP(JSTAR) = ABS(M)
      END IF
      R2MU = 3.0D0*MOR3**C3RD/(CPI4*RHO)                          ! CGS *and* codeunits
! Correction factors FP and FT to the pressure and temperature gradients for
! rotating stars. Could (should?) use the more sophisticated scheme from
! Endal&Sofia 1976 for this. See code by Alexander Heger. This can be adapted
! for doing implicit calculations easily, but will require an extra equation
! to solve for R0...
! Limit FP to 0.1, smaller values seem to give convergence problems -SdM
      FP = max( 0.1D0, 1.0D0 - 2.0D0*C3RD/( MOR3*(CG2*APER)**2 ))
      FT = 1.0D0
C pressure gradient equation; gravity corrected for rotation.
      GOR = CG*MOR3 * FP                                          ! CGS
      PHIMU = 5.0D21*GOR*R2MU                                     ! Codeunits
      APMU = - RHO*PHIMU/P                                        ! Codeunits
      HP = DMIN1(P/(GOR*1.0D11*R0*RHO), DSQRT(P/(CG*RHO*RHO)))    ! CGS
C temperature gradient equation
      GRADR = 0.25D0*FK*P*LOR3/(CPI4*CL*GOR*PR) * FT
      GRADMU = EXPL_VAR(JK, EXPLV_GRADMU, JSTAR)
      DG = GRADR - GRADA
      DG_LEDOUX = DG - CONVECTION_LEDOUX * PHIE*GRADMU/DELTA
C mixing-length theory; the cubic equation for wl/chi has one real solution
      S2 = GRADA*SCP*T
      WC1 = MIN(0.5D0*S2*(CALP*CALP*HP/(9.0D0*XHI))**2, 1.0D302)
      WC2 = MIN(546.75D0*WC1*DMAX1(0.0D0, DG_LEDOUX) + 73.0D0, 1.0D302)
      WC3 = MAX(1.0D-200, CBRT(WC2 + DSQRT(WC2*WC2 + 12167.0D0))) ! Overflow guard
      WC4 = DMAX1(WC3 - 23.0D0/WC3 - 2.0D0,0.0D0)
      WCV = WC4*XHI/(3.0D0*CALP*HP)
      WL = CALP*HP*WCV
C if GRADR <= GRADA, i.e. DG <= 0, we get WCV = 0 and GRADT = GRADR
      GRADT = GRADR - 4.0D0*HP*WCV**3/(CALP*S2*XHI)
      GTA = GRADT - GRADA
      ATMU = GRADT*APMU
C mesh spacing equation, with modified pressure, temperature gradient equations
      VP = CT(4)*AP + CT(5)*DLOG(P + CT(9))
      VPP = CT(4) + CT(5)*P/(P + CT(9))
      VT = CT(7)*DLOG(T/(T + CT(10)))
      VTT = CT(7)*CT(10)/(T + CT(10))
C VR and VM must go to zero at centre like r**2, m**(2/3)
      XMC = CT(6)*CBRT(MC(JSTAR)*MC(JSTAR))
      VMF = XMC + M3*M3
      IF ( JK == KH ) VMF = XMC
      VM = DLOG(DABS(XMC/VMF))
      VMMU = - 1.0D0/VMF
      VR = - CT(3)*DLOG(DABS(R2/CT(8) + 1.0D0))
      VRR = - CT(3)/(R2 + CT(8))
C QQ is the quantity that the meshpoints should be at equal intervals of
      QQ = VP + VT + VM + VR
      QMU = VPP*APMU + VTT*ATMU + VMMU + VRR*R2MU
c L dependence for mesh function
      QMU = QMU + CT(2)*EX*1.5D-8*M3
      MUK = VQK/QMU                                                  ! Codeuints
      MK = 1.5D0*M3*MUK
      MKP(JSTAR) = MK
c Convert d/dm to d/dk = (d/dm)(dm/dk)
      VMK = VMMU*MUK
      VPK = VPP*APMU*MUK
      VTK = VTT*ATMU*MUK
      VRK = VRR*R2MU*MUK
! Calculate ratio of timestep to local radiative diffusion timescale.
! This is used to change the differencing scheme in the energy and diffusion
!  equations. See Sugimoto 1970.
      WT = R2*XHI*DT/(VRK*VRK)
! Above is probably wrong - should it be something like this in stead?
      !WT = XHI*DT/(1.0D22 * R2MU*MUK)
      WT = OFF_CENTRE_WEIGHT * XHI*DT/(R2*1.0D22)
      !IF (XH < 1.0d-9 .and. L>1.0d6) print *, XHI*DT/(R2*1.0d22)
      WT = SIGN(MIN(ABS(WT), 1.0D302), WT)
      TW = WT 
c potential equation, moment of inertia, angular frequency
      VPHI = PHI
      PHIP(JSTAR) = PHI
      PHISP(JSTAR) = PHIS
      PHIK = PHIMU*MUK
      VI = AI
      VIK = R2*MK/1.5D0                                           ! Codeunits
      SI = VIK/MK                                                 ! Codeunits
c Total angular momentum, code units
      XFS(FX_AM) = XVS(NTAM)                                      ! Codeunits
      XFS(FX_AMK) = VIK/APER

C composition equations
C hydrogen equation with ZAMS baryon correction when KCN = 1
      X1 = VX1
      X1T = (2.0*KX*((1-KCN)*RPP + RPC + RPNG + (1-2*KCN)*(RPN + RPO))
     :     + DX1/DT)*MK
C helium equation
      X4 = VX4
      X4T = (4.0*(-KX*(0.5*RPP + RPN + RPO)*(1-KCN)
     :     + KY*(3.0*R3A + RAC + 1.5*RAN + RAO + RANE)
     :     - KZ*(RCC + RCO + 2.0*ROO + RGNE + RGMG)) + DX4/DT)*MK
C carbon equation 
      X12 = VX12
      X12T = (12.0*(KX*(RPC - RPN) - KY*(R3A - RAC)
     :            + KZ*(2.0*(RCC + RCCG) + RCO)) + DX12/DT)*MK
C nitrogen equation
      X14 = VX14
      X14T = (14.0*(KX*(RPN + RPNG - RPC - RPO) + KY*RAN) + DX14/DT)*MK
C oxygen equation
      X16 = VX16
      X16T = (16.0*(KX*(RPO - RPNG) - KY*(RAC - RAO)
     :            + KZ*(RCO + 2.0*ROO - RGNE)) + DX16/DT)*MK
C neon equation
      X20 = VX20
      X20T = (20.0*(KY*(RANE - RAN - RAO) + KZ*(RGNE - RGMG - RCC))
     :            + DX20/DT)*MK
c Magnesium equation
      X24 = VX24
      X24T = (24.0*(KZ*(RGMG - RANE - RCCG - RCO - ROO)) + DX24/DT)*MK

c Artificial composition adjustment
      IF(ADJ_COMP) THEN
         CIF = MIXING_FUDGE*IMPOSE_COMPOSITION_FACTOR
         X1T  = X1T  + CIF*(X1  - GET_IPH(M,  5)) / DT*MK;
         X4T  = X4T  + CIF*(X4  - GET_IPH(M,  9)) / DT*MK;
         X12T = X12T + CIF*(X12 - GET_IPH(M, 10)) / DT*MK;
         X14T = X14T + CIF*(X14 - GET_IPH(M, 16)) / DT*MK;
         X16T = X16T + CIF*(X16 - GET_IPH(M,  3)) / DT*MK;
         X20T = X20T + CIF*(X20 - GET_IPH(M, 11)) / DT*MK;
      ENDIF

C Mixing coefficient, energy equation, mass-transfer equation, central BCs.
      B = PR/PG
      IF ( JK == KH ) THEN
         WOS = 1.0D10
         CON = 0.0D0
         VL = L
         LK = LOM*MUK*DSQRT(XMC)
         LQ = 0.0D0
         FAC = 0.0D0
         XIK = 0.0D0
         VM = M
         VR = R2
      ELSE
         WOS =(2.5D0+B*(2.0D1+1.6D1*B))*(1.5D0*CU*M3/DABS(APMU*M)+1.0D0)
         !CON = 6.0D-22*R2*M3/(R2MU*R2MU*MUK)*WC1**C3RD*XHI  
! The mixing coefficient is calculated approximately and possibly reduced
!  by underrelaxation (using the cube root of WCV is convenient because it
!  has about the right behavior near the edge of the convective zone; physically
!  it is nonsense of course).
         F_WL = MIXING_BOOST*WC1 + (1.0D0 - MIXING_BOOST)*WCV
         CON = 6.0D-22*R2*M3/(R2MU*R2MU*MUK)*CBRT(F_WL)*XHI         
         VML = DSQRT(VMF)/M3
         VL = L*VML                   
         LK = VML*MK*(LOM - C3RD*XMC/VMF*LOR3/MOR3)
         LQ = VML*WTH*SCP*T*APMU*MUK*GTA
         DF0 = CDF*1.0D22*CG*M/R
!         FAC = SMF(1, PHI, DF0)
         FAC = STEP(PHI, DF0)
!         XIK = CMT*DSQRT(SMF(2, 2.0D0*PHIS, DF0))/R*FAC*MK
         XIK = CMT*DSQRT(PSTV(2.0D0*PHIS, DF0))/R*FAC*MK
      END IF

! Smoothen energy generation rate
      !LK = EGEN_SMOOTH*LK + (1.0D0 - EGEN_SMOOTH)*LK_PREV(JK)
      !LQ = EGEN_SMOOTH*LQ + (1.0D0 - EGEN_SMOOTH)*LQ_PREV(JK)

! Conversion factor for diffusion coefficients (converting to code units)
      SGF = (CPI4*R2*RHO)**2/MK
C Fudged convective diffusion coefficient: COS > 0 gives overshooting
C FIXME: when Ledoux convection is used, currently does overshooting from
C the Ledoux boundary, which doesn't make muxh sense - it should do
C overshooting from the Schwarzschild boundary.
      IF ( XH >= 1.0D-7 ) DEG = COS/WOS              
      IF ( XH < 1.0D-7 ) DEG = CPS/WOS              
      EG = DG_LEDOUX + DEG
      UG = PSTV(EG, 0.0D0)       ! used to be  UG = PS(EG)

      SELECT CASE (CONVECTION_SCHEME)
         CASE (1)       ! Eggleton 1972
            ! This is where the convective mixing coefficient is fudged by 
            !  taking the square rather than the cube root, as per the 
            !  mixing-length model.
            SG = CRD*CON*UG*UG

         CASE (2)       ! Pols&Tout 2001, Stancliffe, Tout & Pols 2004
            W = DG / GRADR
            ! Fudge factor beta - beta = inf means no fudge.
            ! Use different values of beta in the exhausted core and in the
            ! envelope; should be set through init.dat, but hardcoded for
            ! now.
            BETA = 1.0D0                     ! Envelope
            IF (XH < 1.0D-7) BETA = 5.0D-5   ! Core; Stancliffe&al 2004
            SG = CON*CBRT(UG) * BETA*W / ( 1.0D0 - (1.0D0-BETA)*W )
         CASE DEFAULT   ! Die
            SG = 0
      END SELECT

      ! Semi-convection, after Langer, Sugimoto & Fricke 1983 (A&A)
      ! CSMC = 0.04 is the effciency of semiconvection, alpha in (10) of LSF83.
      ! Stability condition: DG > 0 (Schwarzschild unstable), DG_LEDOUX<0
      ! (Ledoux unstable)
      DSC = 0.0D0
      IF (DG >= 0.0D0 .AND. DG_LEDOUX < 0.0D0) THEN
         DSC = -CSMC/6.0D0 * XHI * DG / DG_LEDOUX
      END IF
      ! Convert to code units and add to the overall diffusion coefficient
      ! SGF in 1e11 g/cm**2, DSC in g/cm**2, so SGF*DSC is in 1e-22 [1e33 g/s]
      SG = SG + 1.0d-22*DSC*SGF
      
      ! Artificial extra mixing: ARTMIX is the diffusion coefficient in cm^2/s
      ! 0 by default, but can be set from init.dat.
      SG = SG + ARTMIX*SGF
      
      ! For ZAMS runs (KCN = 1) always weakly mix inner and outer meshpoints
      IF (KCN == 1) THEN
         IF ( JK < 0.075D0*KH.OR.JK > 0.8D0*KH ) SG = MIN(SG, 1.0D-4*CRD*CON)
      END IF
      !IF ( JK < 0.075D0*KH ) SG = MIN(SG, 1.0D-2*CRD*CON)

      SG = SG * MIXING_FUDGE
      
      WL = SG
      XIM = XI
      DLRK = 0.0D0
      FAC1 = 0.0D0
      FAC2 = 0.0D0
      V2 = 0.0D0
! Mixing coefficient for thermohaline mixing (except for mu gradient)
      IF (CTH > 0.0D0) THEN
         ! radius is ok ~ 10^10 cm
         ! hp is >~ radius (seems ok)
         ! ca is ok (7.566e-15 in cgs)
         ! cl is ok (3e10 in cgs)
         ! T ~ 10^7 K near the centre
         ! gta ~ 10^-3 (small)
         ! scp ~ 10^8
         ! fk ~ 1
         ! rho ~ 10-100
         ! avmu ~ 0.588
         ! mk ~ -1.6e-3
         SGTH_NUMER = 16.0 * CPI * R2 * HP * CA * CL * T**3
         SGTH_DENOM = -GTA * SCP * FK * RHO * AVMU * ABS(MK)
         ! DTH is in 1e-11 cm^2/s
         DTH = 0.0d0
         IF (SGTH_DENOM /= 0.0d0) DTH = SGTH_NUMER/SGTH_DENOM

         ! correct for radiation pressure, see Kippenhahn, Ruschenplatt & 
         ! Thomas, A&A 91, 1980 equations (34) and (A6).
         ! NB: B = Prad/Pgas, not beta = Pgas/P
         ! Correction factor is phi/delta = beta / (4-3beta) = 1/(1+4*B)
         ! FIXME: should use actual phi/delta from EoS
         DTH = DTH * 1.0/(1.0+4.0*B)

         ! Convert to code mass units (1e33g)
         DTH = 1.0D-33 * DTH

         ! SGF is in 1e11 g/cm**2, so SGTH is in [1e33 g]/s
         SGTH2 = DTH*SGF
         
         ! >0 is 'unphysical'
         SGTH = MIN(CTH * SGTH2, 0.0D0) * MIXING_FUDGE
      ELSE
         SGTH = 0.0D0
      ENDIF

c Angular momentum transport, specific angular momentum
!  ( (R2+CT(8))*DAR - C3RD*R2MU/M3*DM )/R0**2 can be CGS or code units
!  Final units: 1/s^2 * [SI]*[M] = 10^55 g cm^2/s^2
      OMEGAT = ( DOMEGA + 2*OMEGA/R0**2*
     &            ( (R2+CT(8))*DAR - C3RD*R2MU/M3*DM ) )/DT * SI * MK
c Innermost meshpoint is always at constant M&R (or should be!)
      IF (JK == KH) THEN
         OMEGAT = DOMEGA/DT * SI * MK
      ENDIF

      ! The transport coefficient for angular momentum transport is the sum of
      !  different coefficients, but most of these depend on the local shear,
      !  which has to be found from equns1.
      ! In terms of units, the dimension of the mixing coefficients should be
      !  in code units.
      SGAM = CON*UG**2 * MIXING_FUDGE

      ! When relaxing the AM diffusion term, always have some (weak!) mixing,
      !  this helps convergence.
      SGAM = SGAM - (1.0-AMADV_SMOOTH)*1.0D-16
      !SGAM = SGAM - 1.0D-16

      ! Enforce rigid rotation on the first model
      !IF (JMOD==1) SGAM = -1.0D21
      IF (JMOD==1) SGAM = CON
      
      ! Always make outer few meshpoints rotate as a solid body
      ! This is designed to aid numerical stability
      !IF ( JK < 0.075D0*KH ) SGAM = CON
      !IF ( JK > 0.95D0*KH ) SGAM = CON

! Dynamical shear instability. For this, we need both the Richardson number
!  and the size of the unstable region. For the former we need the gradient
!  of omega, for the latter we need the size of the unstable region, here
!  estimated to be one pressure scaleheight (compare notes from Yoon).
! We omit the composition gradient from the Richardson number for simplicity.

! All quantities in the Richardson number should be in CGS units.
!  RHO/P * PSTV(-GTA, 0.0D0)*GOR**2 is in CGS
!  R2MU is in CGS (or code units)
!  MUK/R0 is in 10^22 g^(2/3)/10^11 cm = 10^11 g^(2/3)/cm
!  So R2MU*MUK/R0 is in 10^11 cm
      XFS(FX_RICH) = RHO/P * PSTV(-GTA, 0.0D0)*GOR**2 * (1.0d11*R2MU*MUK/R0)**2

! The mixing coefficient should be in [1e33g]/s. 
!  HP**2/t_dyn is in cm**2/s
!  SGF is in 1e11 g/cm**2
! So SGF*HP**2/t_dyn is in 1e11 g/s = 1e-22 [1e33g]/s
      XFS(FX_DDSI) = AMADV_SMOOTH*CDSI * 1.0D-22 * HP**2/SQRT(CG*MOR3) * SGF * MIXING_FUDGE

! Secular shear instability, apart from domega/dk
      XFS(FX_DSSI) = 0.0D0
      IF (GTA<0.0D0) THEN
         REC = 2500.0D0    ! Critical Reynolds number
         XFS(FX_DSSI) = CSSI*C3RD * XHI * 1.0D-33*R**3 * HP / ( GOR*PSTV(-GTA, 0.0D0)*(MUK*R2MU)**2 )
         XFS(FX_DSSI) = XFS(FX_DSSI) * SGF * MIXING_FUDGE
         XFS(FX_SSSI) = PRANDTL*REC * 1.0D-18*RHO/P * PSTV(-GTA, 0.0D0) * (1.0D11*GOR*R2MU*MUK/R0)**2/8.0D0
      END IF

! Eddington-Sweet circulation; size of the unstable region should be
!  limited to a velocity scale-height; the velocity itself should be
!  reduced by an adverse mu gradient.
! To estimate the size of the unstable region, we use the pressure scale-height
!  instead.
      VES = 0.0D0
      VMU = 0.0D0
      IF (GTA<0.0D0) THEN
         ! Circulation velocity VES is in 1e-11 cm/s
         ! We use the Zahn (1992) expression (3.40), which is equivalent to 
         ! Heger (2000) except in the last term (after some thermodynamics).
         ! Following Maeder&Zahn (1998) we include the gravitational
         ! term in the energy generation rate. This means we can use LOM
         ! instead of ENUC+EN+ENX. NB: LOM != L/M.
         ! This expression now differs from Maeder&Zahn in not including
         ! the stabilising effect of the mu gradient (which is
         ! incorporated here by using "mu currents", wrong according to
         ! Maeder&Zahn) and different factors that arise due to
         ! differential rotation.
         VES = -GRADA/GTA * LOR3 * OMEGA**2/(CG*MOR3 * GOR) * 
     &           (2.0*LOM/LOR3 - 2.0/MOR3 + OMEGA**2/(CG*MOR3*CPI*RHO))/R0**2
         ! Stabilising current driven by molecular weight gradient.
         ! Based on Kippenhahn & al. (1980).
         ! FIXME: add ther "phi/delta grad_mu" term to -GTA in both VES and VMU
         ! First, calculate the local thermal timescale. We need the
         ! typical size of a fluid blob, expression (29) of Kippenhahn & al.
         ! We omit a factor Hp because this gets divided out in the
         ! expression for v_mu anyway.
         ! NU is the dimensionless plasma viscosity. All quantities here
         ! are in physical units, giving VMU in cm/s = 1e11 (1e-11 cm)/s.
         NU = XHI*PRANDTL * 2.5d0 *CL*FK*RHO**2 / PR
         TKH = 24.0d0*DSQRT(0.6*NU*GRADA/(-GTA)) * (B/(1.0+B)) *
     &            SCP / (16.0*CL*CA*T**3)
         VMU = 1.0D11*PHI*GRADMU/(DELTA * GTA * TKH)

!         ! Get this from D_{th} = d * v_{mu} ~ r * v_{mu}, see Heger PhD thesis
!         ! DTH is now in 1e22 cm^2/s, R0 is in 1e11 cm, so DTH/R0 is in 
!         !  1e11 cm^2/s = 1e22 (1e-11 cm^2/s)
!         VMU = 1.0D22*DTH/R0
      END IF
! Diffusion coefficient, according to Heger (originally due to Mestel)
! Export variables to EQUNS1, to complete calculation there.
! VES * HP * SGF should be in 1e33 g/s
!  VES (and VMU) should be in 1e-11 cm/s
!  HP is in cm
!  SGF is in 1e11 g/cm**2
! -> 1e-33 * VES*HP*SGF is in 1e33 g/s
      XFS(FX_DES) = 1.0d-33*CESC*MAX(0.0d0, DABS(VES) - DABS(VMU))*HP*SGF
      !XFS(FX_VES) = AMADV_SMOOTH * VES * MIXING_FUDGE
      !XFS(FX_VMU) = VMU
      !XFS(FX_HP) = HP
      !XFS(FX_SGF) = SGF
! Diffusion coefficient, according to Maeder & Zahn (1998)
! Since this uses R0 rather than HP, the conversion factor is 1e-22
! This doesn't include the stabilising effect of the mu gradient, which
! needs to be added in equns1
      !XFS(FX_DES) = AMADV_SMOOTH*CESC * VES*HP*1.0D-11 * MIXING_FUDGE * SGF
      XFS(FX_DES) = AMADV_SMOOTH*MIXING_FUDGE * CESC * 1d-22*DABS(VES*R0) * SGF
      XFS(FX_DGMU) = 1.0/(GTA * DELTA * MUK * APMU * AVMU)

! Gravitational settling
! Use (3.19) from Montmerle & Michaud, 1976 ApJS 31, ignoring the atomic
! diffusion term (which should be negligible) and including only the
! gravitational settling term.
! Reference species is either H or He, depending on which is more abundant
! (so this will not work if neither is present anymore).
      NREF = 1
      IF (NA(2)>NA(1)) NREF = 2
      IF (CGRS > 0.0D0 .AND. JMOD > 2 .AND. XA(NREF)>0.0D0) THEN
         ! Compute atomic diffusion coefficients relative to the reference
         ! species for all elements [cm^2/s].
         ! Paquette & al. (1986) equation (40)
         CALL DIFFUSION_COEFFICIENTS(RHO, T, Ddiff, NREF)

         ! Compute advection velocities [1e11 cm/s] in the restframe of the
         ! reference particles Montmerle & Michaud (1976) equation (3.19)
         APR = 2.0 * R * 1.0D-11 * APMU/R2MU     ! dlogP/dr, cgs
         DO J=1,9
            GI = NA(J)/(NA(NREF)+NA(J))
            KF_NUM = (1+CAN(J)*KZN(J)*GI) * (1.0+0.5*(KZN(J)+1)*GI)
            KF_DENOM = (1+CAN(J)*GI) * (1.0+0.5*(KZN(J)+1)*KZN(J)*GI)
            WPI(J) = Ddiff(J)*(1.0+GI)*(CBN(J)-1 + KF_NUM/KF_DENOM * (CAN(J)-KZN(J))) / (1.0+CAN(J)*GI) * APR
            ! Convert to code units
            WPI(J) = 1.0d-11 * WPI(J)
         END DO
         ! The velocity of the reference particles in the rest frame of the
         ! star is found from the requirement that the total mass flux must
         ! vanish.
         ! NB: Fe and Si are inert, so don't contribute to the mass flux.
         ! They *do* contribute to the sum of all mass fractions and hence
         ! the normalisation of the sum from which the hydrogen velocity is
         ! calculated.
         WPI(NREF) = 0
         WPP = 0
         DO J=1,7
            WPP = WPP - WPI(J) * XA(J)
         END DO
         ! Correct the normalisation for the presence of Fe and Si.
         WPP = WPP / (1.0-(XFE+XSI))
         ! Now transform all velocities to the rest frame of the star.
         WPI(1:7) = WPI(1:7) + WPP

         ! Export fluxes to equns1
         DO J=1,7
            XFS(FX_FH+J-1) = MIXING_FUDGE * CGRS*WPI(J) * RHO * XA(J) * CPI4*R2
         END DO
      ELSE
         XFS(FX_FH:FX_FH+6) = 0.0D0
      ENDIF

      IF ( JSTAR /= 1 ) THEN
c Put here the *internal* quantities that depend on *TWIN* vbles. Quantities 
c XX needed from both stars should have been stored earlier as XXP(JSTAR) = XX.
c Approx. to binary separation, and potential difference between L1, L2
         APP = OA*OA*(MP(1) + MP(2))/(MP(1)*MP(2)*CG1)**2
         CALL POTENT ( MP(1)/MP(2), DPHI )
         DF0 = CDF*1.0D22*CG*(MP(1) + MP(2))*DPHI/APP 
!         DFL = 1.0D22*CG*(MP(1) + MP(2))*DPHI/APP 
c Differential mass flux XI, semidetached or in contact: XIK is dXI/dK
c Constants CMT, CLT are read in as data (fort.22). SMF is a smoothing
c function to prevent a jump discontinuity at L1 surface.
!         FAC1 = SMF(1, PHIP(1), DF0)
!         FAC2 = SMF(1, PHIP(2), DF0)
         FAC1 = STEP(PHIP(1), DF0)
         FAC2 = STEP(PHIP(2), DF0)
         MMK = 0.5D0*(MKP(1) + MKP(2))
!         V2 = (SMF(2, PHISP(1), DF0) - SMF(2, PHISP(2), DF0))*
!     :      (FAC1 + FAC2) 
         V2 = (PSTV(PHISP(1), DF0) - PSTV(PHISP(2), DF0))*(FAC1 + FAC2)
! Mass transfer: if XIK>0, transfer from *2->*1, otherwise from *1->*2
         XIK = CMT*V2/(DSQRT(DABS(V2) + 1.0D11)*APP)*MMK * MDOT_SMOOTH
c Heat transfer due to differential rotation, CLT rad/sec
         DLRK = - CLT*(ENP(1) - ENP(2))*FAC1*FAC2*MMK
      END IF
      IF ( JI > 0 ) GO TO 4   
      IKK = KH + 2 - JK
      IF ( JSTAR == KTW ) THEN 
C sundry numbers for PRINTB, FUNCS2; SX(40 to 45,JK) for TWIN vbles  
         SX(45, IKK) = FAC1
         SX(44, IKK) = FAC2
         SX(43, IKK) = V2
         SX(42, IKK) = XIK
         SX(41, IKK) = ENP(1) - ENP(2)
         SX(40, IKK) = - DLRK
      ENDIF
c do. for non-TWIN vbles
      ETH = - T*(SF*DAF + ST*DAT)/DT + SCP*T*APMU*GTA*MT/(1.5D0*M3) 
      EGR = 1.0D22*GOR*R2 - U

c FIXME: this does NOT work properly in TWIN mode
! Reaction rates
      SX(50, IKK) = RPP
      SX(51, IKK) = RPC
      SX(52, IKK) = RPNG 
      SX(53, IKK) = RPN
      SX(54, IKK) = RPO
      SX(55, IKK) = RAN

! Cp ds/dlog p
      SX(56, IKK) = SCP * ST/PT
      
! Values of LK and LQ
      SX(57, IKK) = LK
      SX(58, IKK) = LQ

! Thermohaline mixing and convection stuff
      SX(31, IKK) = AVMU
      SX(23, IKK) = SGTH
      SX(30, IKK) = SG
      SX(7, IKK) = GRADT

! Rotation rate and mixing coefficients
      SX(59, IKK) = OMEGA
      SX(60, IKK) = XFS(FX_RICH)
      SX(61, IKK) = XFS(FX_DDSI)
      SX(62, IKK) = XFS(FX_DSSI)
      SX(63, IKK) = XFS(FX_VES)
      SX(64, IKK) = XFS(FX_VMU)
      SX(65, IKK) = SGF

      SX(66, IKK) = XFS(FX_SSSI) ! needed in printb SdM

      SX(76, IKK) = XFS(FX_DES)
      SX(78, IKK) = XFS(FX_DGMU)
c values needed for nucleosynthesis calculations in FUNCS2, EQUNS2
      HT(JSTAR, 1, JK) = RHO
      HT(JSTAR, 2, JK) = SG
      HT(JSTAR, 3, JK) = ZT
      HT(JSTAR, 4, JK) = MK
      HT(JSTAR, 5, JK) = NEO
      HT(JSTAR, 6, JK) = NIO
      HT(JSTAR, 7, JK) = NZZ
      HT(JSTAR, 8, JK) = AVMA
      HT(JSTAR, 9, JK) = NE
      HT(JSTAR, 10:18, JK) = NA(1:9)
 4    IF ( JK /= 1 ) GO TO 2  
C pressure, temperature surface boundary conditions
      ABSL = DABS(L)
      BCP = DLOG(FK/(R*1.0D11*GOR)*(1.5D0*PG + 0.75D0*PR))
      BCT = DLOG(1.0D11*ABSL/(0.75D0*CPI4*CL*R2*PR))

! Eddington factor, used in (some) mass loss prescriptions below
      LOLEDD = FK*LOR3/(CPI4*CL*CG*MOR3)
      LEDD(JSTAR) = L / LOLEDD

! Eddington accretion rate, optionally used to limit accretion rate
      MDOTEDD(JSTAR) = CPI4*R*CL/FK * 1.0D-22   ! in 1e33g/s

! Store accreted abundances for companion star
! For a single star with accretion composition switched on (CCAC=1.0), the
! accreted abundance is read from init.dat
      JSTAR_ACCRETOR = 3-JSTAR
      XAC(1:7, JSTAR_ACCRETOR) = XA(1:7)

! ----------------------------------------------!
! Determine mass loss rate due to stellar winds !
! ----------------------------------------------!
c Metalicity dependent mass loss scaling
c get mass-fraction from baryon-fraction (X1, X4) !SdM  
c     AVM = sum ( mass in amu * baryonfract / baryon nr )
      AVM = CAN(1)*X1/CBN(1)  + !independent
     :      CAN(2)*X4/CBN(2)  +
     :      CAN(3)*X12/CBN(3) +
     :      CAN(4)*X14/CBN(4) +
     :      CAN(5)*X16/CBN(5) +
     :      CAN(6)*X20/CBN(6) +
     :      CAN(7)*XMG/CBN(7) + !sometimes independent
     :      CAN(8)*XSI/CBN(8) + !constant
     :      CAN(9)*XFE/CBN(9)   !constant          
c Determine effective surface metallicity i.o.t. scale Mdot with Z_EFF !SdM
      XH_SURFACE = X1/CBN(1)*CAN(1)/AVM   ! massfraction H surface      
      XHE_SURFACE = X4/CBN(2)*CAN(2)/AVM  ! massfraction He surface 
      Z_EFF = 1.D0 - XH_SURFACE - XHE_SURFACE
      Z_EFF = MIN(1.D0, MAX(Z_EFF, 0.D0))
c mass loss rate for dynamo-driven wind (MDTDDW); Alfven radius squared (RA2)
      RHT = (0.755D0*ABSL**0.47D0 + 0.05D0*L**0.8D0)/M**0.31D0
      RAD = R/RHT
      W2 = R*ABSL/M
      TET = 120.0D0*R/W2**C3RD*RAD**2.7D0
      RO = 1.67D0*DABS(APER/TET)
      ROT = 1.0D0/(RO**2 + 1.0D0)
      MDTDDW = 1.54D-17*W2*RAD**2*ROT**3.67D0
      BP = 5.4D-3*DSQRT(MOR3)*W2**C3RD*RAD**3.4D0*ROT**1.21D0
      IF (BP<1.0D-99) BP = 0.0D0
      ALFAC = 6.2D-08*(R2/MOR3*(BP*BP/MDTDDW)**2)**C3RD
      RA2 = 0.0D0
      IF (CHL>0.0D0) RA2 = CHL*(R*ALFAC)**2

c Eggleton's Reimers-like wind
c cool superwind (Mdot prop. to Lum/(Env. binding energy))
      MDTSW = 0.0D0
      IF ( (BE(JSTAR) /= 0.0D0) .AND. (TN(JSTAR) /= 0.0D0) )
     : MDTSW = M*DMIN1(1.3D-5*L/DABS(BE(JSTAR)), 1.0D1/(TN(JSTAR)*CSY))

! Use smart mass loss routine that determines the mass loss rate based
! on recipes appropriate for the stellar parameters, falling back to the
! de Jager rate when no other applicable rate can be found.
      MDOT = 0.0
      IF (SMART_MASS_LOSS>0.0) THEN
         S_X = CAN(1)*X1/CBN(1)
         S_HE = CAN(2)*X4/CBN(2) 
         S_C = CAN(3)*X12/CBN(3)
         S_O = CAN(5)*X16/CBN(5)
         MDOT = CALCULATE_MDOT(
     &       T, L/CLSN, L/LOLEDD/CLSN, M/CMSN, R/CRSN, S_X, S_HE, S_O, S_C, CZS
     &      )
         MDOT = MDOT * CMSN/CSY * SMART_MASS_LOSS * MDOT_SMOOTH
      ELSE
c Mass-loss rate for luminous stars from de Jager et al (1988)
         MDTL = 0.D0
         IF(CMJ > 0.D0) CALL DJNVDH ( DLOG10(T), DLOG10(L/CLSN), MDTL )
c Mass loss rate for massive stars, Vink et al. (1999, 2000, 2001)
         MDOT_JVINK = 0.0D0
!         T_PREV = EXP(H(2,1))
!         L_PREV = H(8,1)
!         IF (CMV > 0.0D0) MDOT_JVINK = 
!     &      CMV*CALC_MDOT_JVINK(M/CMSN, L_PREV/CLSN, T_PREV)*CMSN/CSY
         IF (CMV > 0.0D0) MDOT_JVINK =
     &      CMV*CALC_MDOT_JVINK(M/CMSN, L/CLSN, T)
c Mass loss rate for Wolf-Rayet stars, Nugis&Lamers 2000
         MDOT_WR = 0.0D0
         IF (CMNL > 0.0D0) THEN
            S_X = CAN(1)*X1/CBN(1)
            S_HE = CAN(2)*X4/CBN(2) 
            S_C = CAN(3)*X12/CBN(3)
            S_O = CAN(5)*X16/CBN(5)
            MDOT_WR = CMNL*
     &         CALC_MDOT_WR_NL(LOG10(L/CLSN),LOG10(T),S_X,S_HE,S_C,S_O,CLOGZ)
     &         * MDOT_SMOOTH
         END IF

c Scale with the effective metallicity: (as Vink et al)
! The Vink rate and the Nugis & Lamers rate already depend on the
! metallicity
         IF (ZSCALING_MDOT /= 0.D0) THEN
            MDTL = MDTL*(Z_EFF/CZSN)**ZSCALING_MDOT
         ENDIF
! The expected mass loss rate is the maximum of the different recipes
! Convert from Msun/yr to code units (1e33g/s)
         MDOT = MAX(CMJ*MDTL, MDOT_JVINK, MDOT_WR) * CMSN/CSY
      END IF
c Total wind mass loss rate, in 1e33 g/s
      ZEP(JSTAR) = CML*MDTDDW + DMAX1(MDOT, CMR*MDTSW)
      IF ( MUTATE .OR. JMOD == 0 ) ZEP(JSTAR) = 0.0D0
      ZEP(JSTAR) = ZEP(JSTAR) * MDOT_SMOOTH

!     Rotation rate over the critical rotation rate, simply defined by
!     setting gravitational acceleration equal to the centrifugal one,
!     taking into account the eddington factor (e.g. Heger, Langer,
!     Woosley 2000, eq 50)
      WOWCRIT = OMEGA/SQRT(CG*MOR3*(1.0D0-LOLEDD))
      WOWCRITP(JSTAR) = WOWCRIT

!     Rotational enhanced mass loss, according to options:
!      1. Heger, Langer & Langer (2000) [they cite Friend & Abbott 1986, but
!         the expression does not appear in that paper]
!      2. Maeder & Meynet (2000) (Paper VI from their series)
      IF (CMDOTROT_HLW>0.0d0) THEN
         ZEP(JSTAR) = ZEP(JSTAR) / (1.0d0 - WOWCRIT)**0.43
      END IF
!     Rotational enhanced mass loss from Maeder & Meynet (2000)
!     Empirical force multipliers, from Lamers et al. 1995
      IF (CMDOTROT_MM>0.0d0) THEN
         ALPHA = 0.52
         IF (AT < 4.30*CLN) THEN
            ALPHA = -0.71 + 0.22*AT/CLN
         ELSEIF (AT > 4.30*CLN .AND. AT < 4.35*CLN) THEN
            ALPHA = 0.236+0.284*(AT/CLN - 4.30)/0.05
         ELSEIF (AT > 5.0*CLN) THEN
            ALPHA = 1.0
         END IF
         ZEP(JSTAR) = ZEP(JSTAR) * 
     &   ( (1.0-LOLEDD) / DABS(1.0 - (2.0*WOWCRIT/3.0)**2 - LOLEDD) )**
     &         (1.0D0/ALPHA-1.0D0)
      END IF
! Artificial mass loss/gain - rate depends on configuration
      MDOT_CMI = 0.0D0
      IF ( JMOD /= 0 ) THEN
         SELECT CASE (CMI_MODE)
            CASE (1) ! Exponential
               MDOT_CMI = MDOT_SMOOTH*CMI*MP(1)
            CASE (2) ! Linear
               MDOT_CMI = MDOT_SMOOTH*CMI*CMSN
         END SELECT
         IF ((M+MDOT_CMI*DT) > UC(13)*CMSN) MDOT_CMI = (M-UC(13)*CMSN)/DT
      END IF
!     Adjust luminosity for the kinitic energy of the wind (photon tiring)
      F_PHOTON_TIRING = CPHOTONTIRE*MIN(1.0D22*GOR/LOR3 * ZEP(JSTAR)/R, 1.0D0)
      L = (1.0D0 - F_PHOTON_TIRING) * L
      VL = L
C pressure, temperature surface boundary conditions
      ABSL = DABS(L)
      BCP = DLOG(FK/(R*1.0D11*GOR)*(1.5D0*PG + 0.75D0*PR))
      BCT = DLOG(1.0D11*ABSL/(0.75D0*CPI4*CL*R2*PR))
!MvdS - transport wind mass loss in Mo/yr to printb
      WML = -ZEP(JSTAR)*CSY/CMSN
c Determination of mass of current component (M), and other cpt (OM).
      IF ( KTW == 1 .AND. JB == 1 ) THEN
         OM = MB - M                     ! M2; M = M1
      END IF
      IF ( KTW == 1 .AND. JB == 2 ) THEN
         ZEP(JSTAR) = 0.0D0
         RA2 = 0.0D0
         BP = 0.0D0
         FDT = AGE*CSY + DT - T0
         IF ( FDT.LE.0.0D0 ) FDT = DT
         M2 = OM0 + FDT*OMTA 
         OA = A0 + FDT*ATA 
         ECC = E0 + FDT*ETA 
         OM = M0 + FDT*MTA                ! M1; M = M2
         MB = M2 + OM
      END IF
      IF ( KTW == 2 ) OM = MB - M
c end of M, OM determination
      BM = MB
      HORB = OA
      SE = ECC
c quantities depending on the orbit
      MU = M*OM/MB
      E2 = ECC*ECC
      EF = 1.0D0/DSQRT(1.0D0 - E2)
      WE25 = EF**5
      SEP = (OA*EF/(CG1*MU))**2/MB
      PER = DSQRT(SEP**3/MB)/CG2
      RL = SEP*RLOBE(CBRT(M/OM))
c quantities for spin up/down due to mass transfer -SdM
C impact parameter for a mass transfer stream: following Ulrichj&Bulger76   
      Rmin_stream = sep * 4.25D-2*( (M/OM)*(1.0D0+(M/OM)) )**0.25 !if current star is accretor
      
      if (Rmin_stream>R) then      
!     Stream misses star -> accretion disk formed Assume AM is accreted
!     rom the inner kepler orbit of the disk, which has radius equal to
!     stellar radiusr
         R_spinup = R
      else
!     Stream hits the star, accretion by impact.  According to U&B76 a
!     disk would have formed with radius of 1.7*rmin_stream, if the
!     accreting star was not in its way.
         R_spinup = 1.7*Rmin_stream
      endif

C Wind accretion, Bipolar reemission,  AM transfer
C -------------------------------------------------
C IMPROVE THIS
      IF ( KTW == 1) THEN
         PAP(2) = CPA*(RLOBE(CBRT(OM/M)))**2
         BRP(2) = CBR   
         SUP(2) = CSU*OA  
         SDP(1) = CSD*OA
      ELSE ! If doing TWIN
C PAP: Partial Accretion of the wind of the companion:
C [TODO We should replace this by Bondi-Hoyle accretion or something like that]
         PAP(JSTAR) = CPA*(R/SEP)**2
C BRP: Bipolar Reemision: i.e. material accreted through RL overflow can
C be reemitted again along the poles. This is the way in TWIN to do
C non-conservative mass transfer. It assumed (later on in the code) that
C the material lost has the specific angular momentum of the orbit of
C the accreting star. For comparison: beta (as for example in de mink
C etal 2008) = Maccr / Mtransfer = (1-BRP). 
C [Todo: make this a parameter in init.dat]
         BRP(JSTAR) = 0.0D0
C SUP: Spin Up: specific angular momentum of the material accreted,
C which is added to the spin of the accreting star.  Warning:
C Switching this of CSU=0 can give strange effects in wider
C binaries: accreting mass but conserving spin leads to slower
C rotation, opposite to what is realistic. Only current solution is
C setting the efficency of spin orbit coupling very high to prevent
C this. Switching it on also has unwanted (but physical) effects: the accreting
C star will reach critical rotation, especially in wider system where the tides are ineffective
         SUP(JSTAR) = CSU*CG1*DSQRT(M*R_spinup)
C SDP: Spin Down: specific angular momentum lost from the donor star.
C The specific angular momentum lost due to RLOF is that of a disk with the
C radius of the star, *not* the specific angular momentum of a sphere with the
C radius of the star.
         SDP(JSTAR) = CSD*(CG1*R**2)/(CG2*APER) 
      END IF

C Radius for determination of RLOF: average the star's current radius and
C its radius on the previous timestep to smooth out kinks in the mass
C transfer rate.
      RPR = DSQRT(DEXP(2.0D0*H(7, 1)) - CT(8))
      RR = 0.5d0 * (R + RPR)
C effective gravity (centrifugal and gravity) at the surface - efective
C gravity at the Rochelobe surface (assuming synchronous rotation of
C stars with the rotation of orbit)
      GMR  = 1.0D22*CG*(MB*C3RD/SEP**3*(R2 - RL*RL) + M*(RL - R)/(R*RL))
      GMRR = 1.0D22*CG*(MB*C3RD/SEP**3*(RR*RR - RL*RL) + M*(RL - RR)/(RR*RL))
c GR terms for circularisation and ang. mom. loss
      RGR = CGRT*MU/(MB*PER)*(SEP/PER)**5*WE25
      ETGR = RGR/9.6D1*(3.04D2 + 1.21D2*E2)*ECC
      OATGR = CGW*RGR*(1.0D0 + 0.875D0*E2)*OA
! MvdS Magnetic braking/1d50 according to Rappaport, Verbunt & Joss, 1983, CMB is the switch
      IF (CMB > 0.0D0) OATMB(JSTAR) = CMB*3.8D-3*M*R**4*(CPI/(APER*4.32D4))**3
      IF(TEFF >= 5.D3) THEN  !Use full MB strength for Teff<5000K
         IF(QCNV < 0.02.AND.QCNV > 0.D0) OATMB(JSTAR) =
     &        OATMB(JSTAR)*DEXP(1.D0 - 2.D-2/QCNV)
         IF(QCNV.LE.0.D0) OATMB(JSTAR) = 0.D0
      END IF
c BC for surface potential eigenvalue
      !BCPH = PHIS + GMRR
      BCPH = PHIS + GMR
c BC for the potential on L1 surface
      BCF = PHI + GMR
! Ratio of difference in potential with RLOF surface and potential at surface
      GMRP(JSTAR) = GMR
      PHIRATIO(JSTAR) = 1.0D-22*GMR / (CG*MOR3*R2 + CG*MB*R2/SEP**3)
c The cheaper model for RLOF, used if CMS > 0 
      RLF = AR - DLOG(RL)
      IF ( CMS > 0.0D0 ) XI = CMS*(PSTV(RLF, 0.0D0))**3
      IF ( CMS > 0.0D0 .AND. JSTAR == 1 ) XI1 = XI
      IF ( CMS > 0.0D0 .AND. JSTAR == 2 ) XI = XI1
      MTR = -XI*CSY/CMSN  !MvdS - transport to printb
c equilibrium-tide model for tidal friction
      WW = 0.4D0*R2*M/(DABS(AI) + 1.0D-10)
      GAM = CTF*(2.0D0/(WW + WW**3.2D0) + RAD**8)
!MvdS Use CSO as spin-orbit coupling (tidal friction) switch
      TFR = CSO*9.0D0*GAM/(3.0D22*M*R2/L)**C3RD*(R/SEP)**8*OM*MB/(M*M)
      WE5 = WE25*WE25
      WE65 = EF**1.3D1
      WE2 = (1.0D0 + E2*(7.5D0 + E2*(5.625D0 + E2*0.3125D0)))*WE65
      WE3 = (1.0D0 + E2*(3.0D0 + E2*0.375D0))*WE5
      WEA = (5.5D0 + E2*(8.25D0 + E2*0.6875D0))*WE5
      WEB = (9.0D0 + E2*(33.75D0 + E2*(16.875D0 + E2*0.703D0)))*WE65
      W1P(JSTAR) = OM/(MB*M)
      W2P(JSTAR) = TFR*WE2
      W3P(JSTAR) = MU*SEP**2/(AI*EF)*TFR*WE3
      W4P(JSTAR) = (RA2 + 2.0D0*C3RD*R2)/AI
! Total spin angular momentum and change in spin angular momentum.
! These expressions are different depending on whether differential rotation
! is taken into account or not.
      IF (RIGID_ROTATION) THEN
! dlog(Iw) = dI/I - dP/P, angular momentum change for solid body rotation.
         DSP(JSTAR) = DAI/AI - DAP/APER
! CG1/CG2 = 2pi/day * 1e5, 2pi/(APER*day) = omega, in s^-1.
         SPP(JSTAR) = CG1*AI/(CG2*APER)
      ELSE
! FIXME: We're not actually interested in the total angular momentum here,
! what we're interested in is that the angular momentum transport through
! the surface works out. In other words, the rate of change of the surface 
! angular momentum should enter the BC. Angular momentum transport in the
! interior then has to make sure that angular momentum is transported back
! to the surface.
! Global conservation of angular momentum is build into the expressions by
! construction (also in the case of advection?) and so *should* be ok.
! Caveat: the code should not crash in the situation where the spin
! direction of the surface is reversed (APER becomes negative). Note the
! way in which APER becomes negative: first APER->+inf, then APER->-inf and
! finally -inf<APER<0. For this reason, it may be better to switch to
! 1/APER as an independent variable.
         DSP(JSTAR) = DXVS(NTAM) / ( XVS(NTAM)+1.0D-16)    ! For differential rotation
         SPP(JSTAR) = CG1 * DXVS(NTAM) / CG2
      END IF
      ET(JSTAR) = ECC*TFR*(WEB - PER/APER*WEA )
      IF ( JI == 0 ) THEN
         SRLF(JSTAR) = RLF
         HSPN(JSTAR) = SPP(JSTAR)
         ZQ = OATGR*DT
      END IF
      XFS(FX_MACC) = 0.0D0
      IF ( KTW == 1 ) THEN ! Non twin case
c BCs for masses, spins, orbital AM and eccentricity in non-TWIN case
         IF ( JB == 1 ) THEN
C     ZET(1): wind mass loss that is not accreted by companion
            ZET(1) = (1.0D0 - PAP(2))*ZEP(1)
C     XIT(1): mass transferred to companion (*2) by RLOF + Wind accretion ?
            XIT(1) = XI + PAP(2)*ZEP(1)
C     ZET(2): mass emitted from pole of companion (after accretion)  
            ZET(2) = BRP(2)*XIT(1)
C     BCM: boundary condition for total mass of *1
            BCM = DMP(1)/DT + ZET(1) + XIT(1) - MDOT_CMI
            ! The accreted mass flux for the star, used to accrete
            ! material of a given composition.
            ! CCAC is a switch that determines whether this composition
            ! change is taken into account or neglected.
            IF (JMOD>0) XFS(FX_MACC) = CCAC*MAX(MDOT_CMI - ZET(1) - XIT(1), 0.0D0)
            BCS = SPP(1)*(DSP(1)/DT + W3P(1) + W4P(1)*ZET(1)) - W2P(1)*OA 
     :               + OATMB(1)
            ! For differential rotation
            !BCS = ( 0.0*DXVS(NTAM) / DT - SI/APER**2*DAP/DT + SI/APER*ZET(1) )
            BCA = DOA/DT + (W1P(1)*ZET(1) + W2P(1))*OA - W3P(1)*SPP(1)
     :               + OATGR + M/(OM*MB)*ZET(2)*OA
            BCE = DE/DT + ET(1) + ETGR
            BCMB = DMB/DT + ZET(1) + ZET(2)
!MvdS: Transport different dHorb/dt's to plotb - case not TWIN
!  SP:spin (wind+mag.br.), MB:mag.br., SO:spin-orbit, ML:system mass loss, GW:grav.waves
           DHSP(1) = W4P(1)*SPP(1)*ZEP(1) ! + OATMB(1) ! ?? Onno
           DHMB(1) = RA2*SPP(1)*ZEP(1)/AI + OATMB(1)
!           DHMT(1) = M/(OM*MB)*ZET(2)*OA
           DHSO(1) = (W2P(1)*OA - W3P(1)*SPP(1))
           DHML    = W1P(1)*ZEP(1)*OA
           DHGW    = OATGR
         ELSE
           BCM = M - M2
         END IF
      END IF
      !MvdS: Total orbital angular momentum loss
      DHDT = DOA/DT
      !MvdS: Need RA2 and AI for both stars to calculate mag.br.
      RA2S(JSTAR) = RA2
      AIS(JSTAR)  = AI
!     Boundary conditions depending on both components
      IF ( JSTAR == 2 ) THEN
         IF(COMPUTE_TWINBC_THEOLDWAY) THEN
c Put here the BCs that depend on *both* components (TWIN case, KTW = 2).
c For PASW and BPRE, determine ZET, XIT from ZEP, XI, PAP, BRP
            DMT0 = 0.001D0*MP(1)/(TN(1) + 1.0D0) + MP(2)/(TN(2) + 1.0D0)
            DO J = 1, 2
c: PA: Partial accretion of the wind of the companion (if that wind is
cstronger than its own wind)
               PA(J) = PAP(J)*STEP(ZEP(J) - ZEP(3 - J), 0.0D0) 
            END DO
c: XI: mass accreted from RLOF + partial wind accretion
            DO J = 1, 2
               XIT(J) = PSTV((3 - 2*J)*XI, 0.0D0) + PA(3 - J)*ZEP(J)
            END DO
! ----------------------------------------------------------------------------
!        Non conservative mass transfer
! ----------------------------------------------------------------------------
!        There are several ways to treat non-conservative mass transfer:
!
!         1. The mass accretion rate can be limited by the Eddington rate
!            of the accretor (caveat: it is not clear whether this applies
!            when accreting from a disk rather than through direct impact of
!            the tidal stream; naively not, but it depends on whether the
!            accretion disk becomes thick or not).
!         2. Angular momentum limit: quench mass accretion rate if the star
!            is close to critical rotation. This really only works if the star
!            is spun up by the accretion.
!         3. Radiative feedback from the hot-spot, in case of a direct impact.
!         4. Just set a constant fraction of material to be accreted  
!
!        For now, we implement (1), (2) and (4)
! ----------------------------------------------------------------------------

!        Eddington-limited accretion (based on Beer, Dray, King & Wynn 2007)
!        Fraction of transferred material lost (confusingly, this is sometimes
!        called "beta" and sometimes "1-beta", although Beer & al. call it psi).
            DO J = 1, 2
               FMDOTLOST(J) = 0.0D0
!               if (ji == 0) print '(1X,1P,I3,3E24.16)', J,
!     &            XIT(3-J)*CSY/CMSN, MDOTEDD(J)*CSY/CMSN,
!     &            (LEDD(3-J)/(1.0d22*PHISP(J)*DABS(XIT(3-J)))+1.0d0)*PHIRATIO(J)
               IF (XIT(3-J) >= MDOTEDD(J) .AND. CMTEL > 0.0d0) THEN
!        Cannot accrete all, lose a bit
!                  FMDOTLOST(J) = 1.0d0 - MDOTEDD(J)/XIT(3-J)
!        This is expression (7) in Beer & al, but with factors rearranged and
!        using the potential at the Roche lobe surface in lieu of the potential
!        in L1 (these are similar and ought to be the same anyway).
                 FMDOTLOST(J) = (1.0d0-LEDD(J)/(DABS(GMRP(J)*XIT(3-J))))
                 FMDOTLOST(J) = FMDOTLOST(J) * PHIRATIO(J)
                 FMDOTLOST(J) = MIN(1.0d0, MAX(FMDOTLOST(J), 0.0d0))
               END IF
            ENDDO

!        Combine all contributions to non-conservative mass transfer to get the
!        total fraction of material ejected.
!        Get total fraction of mass accreted by multiplying all terms that may
!        contribute to non-conservative mass loss (alternative: use the most
!        restrictive).
!        CBR stands for "Bipolar reemission" and is a constant in init.dat
            DO J = 1, 2
               BRP(J) = (1.0d0 - CBR)
               BRP(J) = BRP(J)*(1.0d0 - CMTEL * FMDOTLOST(J))
               BRP(J) = BRP(J)*(1.0d0 - CMTWL * WOWCRITP(JSTAR))
               BRP(J) = 1.0d0 - BRP(J)
            END DO

c: Calculate fraction of material ejected from the system
            DO J = 1, 2
               BR(J) = BRP(J)*STEP(XIT(3 - J) - XIT(J), 0.0D0)
            END DO
c: Mass ejected from the system: the fraction of the wind that is not accreted
c  by the companion and the contribution from non-conservative Roche lobe
c  overflow
            DO J = 1, 2
               ZET(J) = (1.0D0 - PA(3 - J))*ZEP(J) + BR(J)*XIT(3 - J)
            END DO
! The net mass accretion flux for each of the two components
            IF (JMOD>0) THEN
               Y(FX_MACC)   = CCAC*MAX(MDOT_CMI - ZET(1) + XIT(2) - XIT(1), 0.0D0)
               XFS(FX_MACC) = CCAC*MAX(MDOT_CMI - ZET(2) - XIT(2) + XIT(1), 0.0D0)
            END IF
! Mass boundary condition for each of the two components...
            Y(27) = DMP(1)/DT + ZET(1) + XIT(1) - XIT(2) ! *1
            BCM   = DMP(2)/DT + ZET(2) - XIT(1) + XIT(2) ! *2
! ... and for the binary
            BCMB  = DMB/DT + ZET(1) + ZET(2) ! binary
! AM boundary condition for each of the two components ...
            Y(33) = SPP(1)*(DSP(1)/DT + W3P(1) + W4P(1)*ZET(1)) - W2P(1)*OA 
     :           - SUP(1)*XIT(2) + SDP(1)*XIT(1) + OATMB(1) ! *1
            BCS   = SPP(2)*(DSP(2)/DT + W3P(2) + W4P(2)*ZET(2))-W2P(2)*OA 
     :           - SUP(2)*XIT(1) + SDP(2)*XIT(2) + OATMB(2) ! *2
! ... and for the binary
            BCA   = DOA/DT + (W1P(1)*ZET(1) + W2P(1))*OA - W3P(1)*SPP(1)
     :           + (W1P(2)*ZET(2) + W2P(2))*OA - W3P(2)*SPP(2) + OATGR
     :           + (SUP(1) - SDP(2))*XIT(2) + (SUP(2) - SDP(1))*XIT(1) ! orbit
!Eccentricity boudary condition
            BCE   = DE/DT + ET(1) + ET(2) + ETGR
!MvdS: Export different contributions to dHorb/dt's to PLOTB
!  SP:spin (wind+mag.br.), MB:mag.br., SO:spin-orbit, ML:system mass loss, GR:grav.waves
            DHSP(1) = W4P(1)*SPP(1)*ZET(1)
            DHSP(2) = W4P(2)*SPP(2)*ZET(2)
            DHMB(1) = RA2S(1)*SPP(1)*ZET(1)/AIS(2) + OATMB(1)
            DHMB(2) = RA2S(2)*SPP(2)*ZET(2)/AIS(2) + OATMB(2)
            DHSO(1) = (W2P(1)*OA - W3P(1)*SPP(1))
            DHSO(2) = (W2P(2)*OA - W3P(2)*SPP(2))
            DHML = (W1P(1)*ZET(1) + W1P(2)*ZET(2))*OA
            DHGW = OATGR

         ELSE !New twin boundary equations (needed for spinup)

! Reformulation of boundary conditions will be implemented here, .....-SdM

         ENDIF

      END IF
      ! End of surface boundary condition block
 2    CONTINUE
      ! Copy output from NAMEOUT to INF
      CALL NAMES1 ( JSTAR, 2 )
      ! Store variables for explicit calculations
      if (ji <= 0) then
         expl_var(JK, EXPLV_AVMU, JSTAR) = AVMU_NEW
         expl_var(JK, EXPLV_LOGP, JSTAR) = AP
      end if

      IF ( JSTAR == KTW .OR. JI < 0 ) RETURN
      JSTAR = 2
      GO TO 1                            
      END

      FUNCTION STEP (F, D)
c Smoothed step function (=0 for F<=-D, =1 for F>>D)
c D is the smoothing length
      IMPLICIT REAL*8 (A-H, L-Z)
      STEP = SMF (1, F, D)
      RETURN
      END

      FUNCTION PSTV (F, D)
c Smoothed triangle function (=0 for F<=-D, =F for F>>D)
c D is the smoothing length
      IMPLICIT REAL*8 (A-H, L-Z)
      PSTV = SMF (2, F, D)
      RETURN
      END

      FUNCTION SMF (IHT, F, D)
c Smoothed step function (IHT=1) or triangle function (IHT=2)
      IMPLICIT REAL*8 (A-H, L-Z)
      X = F + D  
      IF ( X.LE.1.0D-150 ) THEN  ! Make sure that X**2>0 numerically
         SMF = 0.0D0
      ELSE
         IF(X < 1.D130) THEN    ! safety net for big numbers !SdM
            Y = X*X             ! moved to location where it is needed !SdM
            SMF = Y/(D*D + Y)
         ELSE
            SMF = 1.D0
         ENDIF
      END IF
      IF ( IHT == 2 ) SMF = X*SMF           
      RETURN
      END

      SUBROUTINE EQUNS1 ( JK, KL, KQ )
      USE MESH
      USE MESH_ENC
      USE EXTRA_ELEMENTS
      USE CONSTANTS
      USE SETTINGS
      USE FUDGE_CONTROL
      IMPLICIT REAL*8 (A-H, L-Z)    
      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, KW(260)
c in NAMES2, X of /INE/ is equivalenced to /NAMEIN/, Y to /NAMOUT/
      COMMON /NAMEIN/ BCP(3), BCT(3), VP(3), VPK(3), VR(3), VRK(3), 
     : VT(3), VTK(3), L(3), LK(3), LQ(3), MT(3), VM(3), VMK(3), SG(3), 
     : WT(3), X1(3), X1T(3), X16(3), X16T(3), X4(3), X4T(3), X12(3), 
     : X12T(3), X20(3), X20T(3), BCM(3), VI(3), VIK(3), PHI(3), PHIK(3), 
     : BCF(3), BCS(3), BCPH(3), X14(3), X14T(3), AVMU(3), SGTH(3), OMEGA(3),
     : OMEGAT(3), SGAM(3), SI(3),
     : BCA(3), BCE(3), XI(3), XIK(3), DLRK(3), BCMB(3), X24(3), X24T(3), 
     : MENC(3), MENCT(3), MEA(3), MEAT(3), MET(3), METT(3), 
     : FILLUP_COMMONBLOCK_NAMEIN_EQUNS1(3*(NFUNC-56)) 
      COMMON /NAMOUT/ EQU(NFUNC)
      ! Extra variables
      COMMON /XIN2  / XVARS(3, NSFUNC+1:NSFUNC+NXFSTAR), XVARB(3, NSFUNC+2*NXFSTAR+1:NSFUNC+2*NXFSTAR+NXFBIN)
      ! Extra equations
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /XOUT2 / XEQUS(NSEQ+1:NSEQ+NXESTAR), XEQUB(NSEQ+2*NXESTAR+1:NSEQ+2*NXESTAR+NXVBIN)
      COMMON /ABUND / XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, XW(14)
C VAR(3),(2),(1): values at current, previous and anteprevious meshpoints.
C (3) is nearest centre if KL=0, KQ=1; nearest surface if KL=1, KQ=-1
      DOUBLE PRECISION :: X1AC, X4AC, X12AC, X14AC, X16AC, X20AC, X24AC
      DOUBLE PRECISION :: XAC(7, 2)
      COMMON /ACCRET/ X1AC, X4AC, X12AC, X14AC, X16AC, X20AC, X24AC, XAC
      JSTAR = 1
 1    CALL NAMES2 ( JSTAR, 1 )
C second-order difference equations at interior points
      IF ( 3.LE.JK + KL .AND. JK.LE.KH ) THEN
         
         ! Molecular weight gradient, for thermohaline mixing and angular momentum transport
         DMU12 = KQ * (AVMU(1)-AVMU(2))
         DMU23 = KQ * (AVMU(2)-AVMU(3))

         ! Gradient in rotational velocity, for rotational shear
         DW12 = (OMEGA(2) - OMEGA(1))**2
         DW23 = (OMEGA(3) - OMEGA(2))**2

         ! Richardson number (minus omega gradient), for shear instability
         RICH12 = 0.5*(XVARS(1, FX_RICH)+XVARS(2, FX_RICH))
         RICH23 = 0.5*(XVARS(2, FX_RICH)+XVARS(3, FX_RICH))
         
         ! Mixing coefficient for dynamical shear instability
         DDSI12 = 0.0D0
         DDSI23 = 0.0D0
         IF ( DW12 > 0.0D0 .AND. RICH12>0.0D0 .AND. RICH12/DW12 < 1.0D0 ) THEN
            TON = (1.0D0 - RICH12/DW12)**2   ! Turn on factor
            DDSI12 = 0.5D0*( XVARS(1, FX_DDSI)+XVARS(2, FX_DDSI) )*TON
         END IF
         IF ( DW23 > 0.0D0 .AND. RICH23>0.0D0 .AND. RICH23/DW23 < 1.0D0 ) THEN
            TON = (1.0D0 - RICH23/DW23)**2   ! Turn on factor
            DDSI23 = 0.5D0*( XVARS(3, FX_DDSI)+XVARS(2, FX_DDSI) )*TON
         END IF

         ! Mixing coefficients secular shear instability
         RIS12 = 0.5D0*(XVARS(1, FX_SSSI)+XVARS(2, FX_SSSI))
         RIS23 = 0.5D0*(XVARS(2, FX_SSSI)+XVARS(3, FX_SSSI))
         DSSI12 = 0.0D0
         DSSI23 = 0.0D0
         IF ( DW12 > 0.0D0 .AND. RIS12>0.0D0 .AND. RIS12/DW12 < 1.0D0 ) THEN
            DSSI12 = 0.5D0*(XVARS(1, FX_DSSI)+XVARS(2, FX_DSSI))*DW12
         END IF
         IF ( DW23 > 0.0D0 .AND. RIS23>0.0D0 .AND. RIS23/DW23 < 1.0D0 ) THEN
            DSSI23 = 0.5D0*(XVARS(2, FX_DSSI)+XVARS(3, FX_DSSI))*DW23
         END IF

         !HP12 = 0.5*(XVARS(1, FX_HP)+XVARS(2, FX_HP))
         !VES12 = 0.5*(XVARS(1, FX_VES)+XVARS(2, FX_VES))
         !VMU12 = 0.5*(XVARS(1, FX_VMU)+XVARS(2, FX_VMU)) * CFMU*DMU12
         !VES12 = PSTV(DABS(VES12) - DABS(VMU12), 0.0D0)
         !print *, JK, 1.0D-11*HP12*VES12, 1.0d22*DDSI12/XVARS(1, FX_SGF)

!         ! Eddington-Sweet circulation, compute net circulation velocity
!         HP12 = 0.5*(XVARS(1, FX_HP)*XVARS(1, FX_SGF)+XVARS(2, FX_HP)*XVARS(2, FX_SGF))
!         HP23 = 0.5*(XVARS(2, FX_HP)*XVARS(2, FX_SGF)+XVARS(3, FX_HP)*XVARS(3, FX_SGF))
!         VMU12 = 0.5*(XVARS(1, FX_VMU)+XVARS(2, FX_VMU)) * CFMU*DMU12
!         VMU23 = 0.5*(XVARS(2, FX_VMU)+XVARS(3, FX_VMU)) * CFMU*DMU23
!         VES12 = 0.5*(XVARS(1, FX_VES)+XVARS(2, FX_VES))
!         VES23 = 0.5*(XVARS(2, FX_VES)+XVARS(3, FX_VES))
!         VES12 = PSTV(DABS(VES12) - DABS(VMU12), 0.0D0)
!         VES23 = PSTV(DABS(VES23) - DABS(VMU23), 0.0D0)

         ! Eddington-Sweet circulation, diffusion approximation following Heger.
         DES12 = 0.5*(XVARS(1, FX_DES) + XVARS(2, FX_DES))
         DES23 = 0.5*(XVARS(2, FX_DES) + XVARS(3, FX_DES))

         ! Alternatively, the expression from Maeder & Zahn (1998)
         ! Factor by which to multiply VES to include stabilising effect of
         ! mu gradient
!         VESMU121 = 1.0d0/max(1.0d0, 1.0d0-XVARS(1, FX_DGMU)*DMU12)
!         VESMU122 = 1.0d0/max(1.0d0, 1.0d0-XVARS(2, FX_DGMU)*DMU12)
!         VESMU231 = 1.0d0/max(1.0d0, 1.0d0-XVARS(3, FX_DGMU)*DMU23)
!         VESMU232 = 1.0d0/max(1.0d0, 1.0d0-XVARS(2, FX_DGMU)*DMU23)
!         VESMU12 = 0.5d0*(VESMU121+VESMU122)
!         VESMU23 = 0.5d0*(VESMU231+VESMU232)
!
!         DES12 = 0.5*(VESMU121*XVARS(1, FX_DES) + VESMU122*XVARS(2, FX_DES))
!         DES23 = 0.5*(VESMU231*XVARS(3, FX_DES) + VESMU232*XVARS(2, FX_DES))

         ! Combined diffusion coefficients for chemical mixing
         ! Convection and thermohaline mixing
         MT_SMOOTH = 0.0D0
         S12 = 0.5D0*(SG(1) + SG(2)) - MIXING_FUDGE*PSTV(KQ*MT(2), MT_SMOOTH)
     &       + 0.5D0*(SGTH(1)+SGTH(2))*PSTV(DMU12, 0.0D0)
         S23 = 0.5D0*(SG(2) + SG(3)) - MIXING_FUDGE*PSTV(-KQ*MT(2), MT_SMOOTH)
     &       + 0.5D0*(SGTH(2)+SGTH(3))*PSTV(DMU23, 0.0D0)

         ! Dynamical shear instability
         S12 = S12 + CFC*DDSI12
         S23 = S23 + CFC*DDSI23

         ! Secular shear instability
         S12 = S12 + CFC*DSSI12
         S23 = S23 + CFC*DSSI23

         ! Eddington-Sweet circulation
         S12 = S12 + CFC*DES12
         S23 = S23 + CFC*DES23

         EQU(1)     = S23*(X1(3)  - X1(2))  - S12*(X1(2) - X1(1))  - X1T(2)
         EQU(2)     = S23*(X16(3) - X16(2)) - S12*(X16(2) -X16(1)) - X16T(2)
         EQU(3)     = S23*(X4(3)  - X4(2))  - S12*(X4(2) - X4(1))  - X4T(2)
         EQU(4)     = S23*(X12(3) - X12(2)) - S12*(X12(2) -X12(1)) - X12T(2)
         EQU(5)     = S23*(X20(3) - X20(2)) - S12*(X20(2) -X20(1)) - X20T(2)
         EQU(EN14)  = S23*(X14(3) - X14(2)) - S12*(X14(2) -X14(1)) - X14T(2)
         EQU(EMG24) = S23*(X24(3) - X24(2)) - S12*(X24(2) -X24(1)) - X24T(2)

         ! Add advection terms (contributions from gravitational settling)
         IF (CGRS > 0.0D0) THEN
            EQU(1)     = EQU(1)    - KQ*(XVARS(3, FX_FH)  - XVARS(2, FX_FH))
            EQU(2)     = EQU(2)    - KQ*(XVARS(3, FX_FO) - XVARS(2, FX_FO))
            EQU(3)     = EQU(3)    - KQ*(XVARS(3, FX_FHe)  - XVARS(2, FX_FHe))
            EQU(4)     = EQU(4)    - KQ*(XVARS(3, FX_FC) - XVARS(2, FX_FC))
            EQU(5)     = EQU(5)    - KQ*(XVARS(3, FX_FNe) - XVARS(2, FX_FNe))
            EQU(EN14)  = EQU(EN14) - KQ*(XVARS(3, FX_FN) - XVARS(2, FX_FN))
            EQU(EMG24)  = EQU(EMG24) - KQ*(XVARS(3, FX_FMG) - XVARS(2, FX_FMG))
         END IF
         EQU(ESUMX) = EQU(1)+EQU(2)+EQU(3)+EQU(4)+EQU(5)+EQU(EN14)+EQU(EMG24);

         ! Angular momentum transport
         ! The mixing coefficients here need to be weighed with the specific
         !  moment of inertia, so we need to recalculate some quantities
         MT_SMOOTH = 1.0D-22
         MT_SMOOTH = 1.0D-3*DABS(MT(2))
         MT_SMOOTH = 0.0D0
         S12 = 0.5D0*(SGAM(1)*SI(1) + SGAM(2)*SI(2)) - 1.0D0*PSTV(KQ*MT(2), MT_SMOOTH)*SI(2)
         S23 = 0.5D0*(SGAM(2)*SI(2) + SGAM(3)*SI(3)) - 1.0D0*PSTV(-KQ*MT(2), MT_SMOOTH)*SI(2)
         !S12 = 0.5D0*(SGAM(1)*SI(1) + SGAM(2)*SI(2)) - MT(2)*SI(2)
         !S23 = 0.5D0*(SGAM(2)*SI(2) + SGAM(3)*SI(3)) - MT(3)*SI(3)
         !S12 = ( 0.5D0*(SGAM(1) + SGAM(2)) - PSTV(KQ*MT(2), 0.0D0) )*SI(2)
         !S23 = ( 0.5D0*(SGAM(2) + SGAM(3)) - PSTV(-KQ*MT(2), 0.0D0) )*SI(2)
         ! Dynamical shear instability
         IF ( DW12 > 0.0D0 .AND. RICH12>0.0D0 .AND. RICH12/DW12 < 1.0D0 ) THEN
            TON = (1.0D0 - RICH12/DW12)**2   ! Turn on factor
            DDSI12 = 0.5D0*( XVARS(1, FX_DDSI)*SI(1)+XVARS(2, FX_DDSI)*SI(2) )*TON
            !DDSI12 = 0.5D0*( XVARS(1, FX_DDSI)+XVARS(2, FX_DDSI) )*SI(2)*TON
         END IF
         IF ( DW23 > 0.0D0 .AND. RICH23>0.0D0 .AND. RICH23/DW23 < 1.0D0 ) THEN
            TON = (1.0D0 - RICH23/DW23)**2   ! Turn on factor
            DDSI23 = 0.5D0*( XVARS(3, FX_DDSI)*SI(3)+XVARS(2, FX_DDSI)*SI(2) )*TON
            !DDSI23 = 0.5D0*( XVARS(3, FX_DDSI)+XVARS(2, FX_DDSI) )*SI(2)*TON
         END IF
         S12 = S12 + DDSI12
         S23 = S23 + DDSI23
         
         ! Secular shear instability
         DSSI12 = 0.0D0
         DSSI23 = 0.0D0
         IF ( DW12 > 0.0D0 .AND. RIS12>0.0D0 .AND. RIS12/DW12 < 1.0D0 ) THEN
            DSSI12 = 0.5D0*(XVARS(1, FX_DSSI)*SI(1)+XVARS(2, FX_DSSI)*SI(2))*DW12
         END IF
         IF ( DW23 > 0.0D0 .AND. RIS23>0.0D0 .AND. RIS23/DW23 < 1.0D0 ) THEN
            DSSI23 = 0.5D0*(XVARS(2, FX_DSSI)*SI(2)+XVARS(3, FX_DSSI)*SI(3))*DW23
         END IF
         !print *, JK, DSSI12, DSSI23
         S12 = S12 + DSSI12
         S23 = S23 + DSSI23

!         ! Eddington-Sweet circulation
!         ! FIXME: use the diffusion coefficient exported from funcs1, as
!         ! for chemical transport - or use the "proper" advection treatment.
!         HP12 = 0.5*(XVARS(1, FX_HP)*XVARS(1, FX_SGF)*SI(1)+XVARS(2, FX_HP)*XVARS(2, FX_SGF)*SI(2))
!         HP23 = 0.5*(XVARS(2, FX_HP)*XVARS(2, FX_SGF)*SI(2)+XVARS(3, FX_HP)*XVARS(3, FX_SGF)*SI(3))
!         S12 = S12 + 1.d-33*CESC*VES12 * HP12
!         S23 = S23 + 1.d-33*CESC*VES23 * HP23

         ! Eddington-Sweet circulation, diffusion approximation following Heger.
         DES12 = 0.5*(XVARS(1, FX_DES)*SI(1) + XVARS(2, FX_DES)*SI(2))
         DES23 = 0.5*(XVARS(2, FX_DES)*SI(2) + XVARS(3, FX_DES)*SI(3))
         S12 = S12 + DES12
         S23 = S23 + DES23

         
         ! Goldreich-Schubert-Fricke instability
         ! FIXME
         ! The mixing coefficient is similar to that of the Eddington-Sweet
         !  circulation, and has a similar order of magnitude. For now, we take
         !  the easy way out and just apply the Eddington-Sweet term again
         ! Stability condition: j=i omega increases outward
         VGSF12 = VES12 * (OMEGA(2) - OMEGA(1))/(VR(2) - VR(1))
         VGSF23 = VES23 * (OMEGA(3) - OMEGA(2))/(VR(3) - VR(2))
         IF ( (SI(2+KL)*OMEGA(2+KL) - SI(3-KL)*OMEGA(3-KL)) < 0) THEN
            S23 = S23 + 1.d-33*CGSF*VGSF23 * HP23
         END IF
         IF ( (SI(1+KL)*OMEGA(1+KL) - SI(2-KL)*OMEGA(2-KL)) < 0) THEN
            S12 = S12 + 1.d-33*CGSF*VGSF12 * HP12
         END IF

         ! Diffusion equation for angular momentum
         XEQUS(EAMT) = S23*(OMEGA(3) - OMEGA(2)) - S12*(OMEGA(2) - OMEGA(1)) - OMEGAT(2)
      END IF
C next-to-surface boundary conditions for second-order equations
      IF ( JK + KL == 2 ) THEN
c Modified equations for accreting of matter with different composition
         MTA = PSTV(MT(2), 0.0D0)
         SG2 = 0.5D0*(SG(2)+SG(3))
     &       + 0.5D0*CTH*(SGTH(2)+SGTH(3))*PSTV(KQ * (AVMU(2)-AVMU(3)), 0.0D0)

         ! Thermohaline mixing of surface layer
         ! Accretion abundance set in FUNCS1
         !X1AC =    0.598490994028748
         !X4AC =    0.381392784203137
         !X16AC =   1.009965586641970E-002
         !X12AC =   3.539808668161780E-003
         !X20AC =   1.851094969072195E-003
         !X14AC =   1.045628497798272E-003
         !X24AC =   6.838266745312392E-004
         
         IF (KL == 0) THEN    ! (3) is nearest to centre, so we have the grid points | 3 | 2 | AC |
            S12 = 0.5D0*(SG(1) + SG(2)) - MIXING_FUDGE*PSTV(KQ*MT(2), 0.0D0)
            S23 = 0.5D0*(SG(2) + SG(3)) - MIXING_FUDGE*PSTV(-KQ*MT(2), 0.0D0)
            EQU(1)     = S23*(X1(3)  - X1(2))  - S12*(X1(2) - XAC(1,JSTAR))   - X1T(2)
            EQU(2)     = S23*(X16(3) - X16(2)) - S12*(X16(2) -XAC(5,JSTAR))  - X16T(2)
            EQU(3)     = S23*(X4(3)  - X4(2))  - S12*(X4(2) - XAC(2,JSTAR))   - X4T(2)
            EQU(4)     = S23*(X12(3) - X12(2)) - S12*(X12(2) -XAC(3,JSTAR))  - X12T(2)
            EQU(5)     = S23*(X20(3) - X20(2)) - S12*(X20(2) -XAC(6,JSTAR))  - X20T(2)
            EQU(EN14)  = S23*(X14(3) - X14(2)) - S12*(X14(2) -XAC(4,JSTAR))  - X14T(2)
            EQU(EMG24) = S23*(X24(3) - X24(2)) - S12*(X24(2) -XAC(7,JSTAR))  - X24T(2)
         ELSE                 ! (3) is nearest to surface, we have the grid points | 1 | 2 | 3 | AC |
            S12 = 0.5D0*(SG(2) + SG(3)) - PSTV(KQ*MT(3), 0.0D0)
     &          + 0.5D0*(SGTH(2)+SGTH(3))*PSTV(KQ * (AVMU(2)-AVMU(3)), 0.0D0)
            S23 = -XVARS(3, FX_MACC) 

            EQU(1)     = S23*(XAC(1,JSTAR)   - X1(3))  - S12*(X1(3) - X1(2))  - X1T(3)
            EQU(2)     = S23*(XAC(5,JSTAR)  - X16(3)) - S12*(X16(3) -X16(2)) - X16T(3)
            EQU(3)     = S23*(XAC(2,JSTAR)   - X4(3))  - S12*(X4(3) - X4(2))  - X4T(3)
            EQU(4)     = S23*(XAC(3,JSTAR)  - X12(3)) - S12*(X12(3) -X12(2)) - X12T(3)
            EQU(5)     = S23*(XAC(6,JSTAR)  - X20(3)) - S12*(X20(3) -X20(2)) - X20T(3)
            EQU(EN14)  = S23*(XAC(4,JSTAR)  - X14(3)) - S12*(X14(3) -X14(2)) - X14T(3)
            EQU(EMG24) = S23*(XAC(7,JSTAR)  - X24(3)) - S12*(X24(3) -X24(2)) - X24T(3)
         END IF

         ! Advection terms (from gravitational settling)
         IF (CGRS > 0.0D0) THEN
            EQU(1)     = EQU(1)    - KQ*(XVARS(3, FX_FH)  - XVARS(2, FX_FH))
            EQU(2)     = EQU(2)    - KQ*(XVARS(3, FX_FO) - XVARS(2, FX_FO))
            EQU(3)     = EQU(3)    - KQ*(XVARS(3, FX_FHe)  - XVARS(2, FX_FHe))
            EQU(4)     = EQU(4)    - KQ*(XVARS(3, FX_FC) - XVARS(2, FX_FC))
            EQU(5)     = EQU(5)    - KQ*(XVARS(3, FX_FNe) - XVARS(2, FX_FNe))
            EQU(EN14)  = EQU(EN14) - KQ*(XVARS(3, FX_FN) - XVARS(2, FX_FN))
            EQU(EMG24)  = EQU(EMG24) - KQ*(XVARS(3, FX_FMG) - XVARS(2, FX_FMG))
         END IF

         EQU(ESUMX) = EQU(1)+EQU(2)+EQU(3)+EQU(4)+EQU(5)+EQU(EN14)+EQU(EMG24);

         ! Angular momentum transport boundary condition (use normal BC ?)
         S12 = 0.5D0*(SGAM(1)*SI(1) + SGAM(2)*SI(2)) - PSTV(KQ*MT(2), 0.0D0)*SI(2)
         S23 = 0.5D0*(SGAM(2)*SI(2) + SGAM(3)*SI(3)) - PSTV(-KQ*MT(2), 0.0D0)*SI(2)
         ! Dynamical shear instability
         ! FIXME: add dynamical shear instability for surface layer
         ! Always need at least a little bit of mixing to correlate meshpoints
         XEQUS(EAMT) = S23*(OMEGA(3) - OMEGA(2)) - S12*(OMEGA(2) - OMEGA(1)) - OMEGAT(2)
         !XEQUS(EAMT) = -S23*(OMEGA(3) - OMEGA(2)) - OMEGAT(3)
         !XEQUS(EAMT) = OMEGAT(2)
         !print *, XEQUS(EAMT), S12, S23
         !print *, OMEGAT(3), OMEGA(3), OMEGA(2)
         !print *, OMEGAT(3), SI(3)*OMEGA(3)
         !print *, MT(1), MT(2), MT(3)
         XEQUS(EAMT) = OMEGA(3) - OMEGA(2)
      END IF
C next-to-central boundary conditions for second-order equations
      IF ( JK + KL == KH + 1 ) THEN
         S23 = KQ*0.5D0*(SG(2)+SG(3)) 
     &       + KQ*0.5D0*(SGTH(2)+SGTH(3))*PSTV(KQ * (AVMU(2)-AVMU(3)), 0.0D0)
         EQU(1)     = S23*(X1(3)  - X1(2))  + X1T(3 - KL)
         EQU(2)     = S23*(X16(3) - X16(2)) + X16T(3 - KL)
         EQU(3)     = S23*(X4(3)  - X4(2))  + X4T(3 - KL)
         EQU(4)     = S23*(X12(3) - X12(2)) + X12T(3 - KL)
         EQU(5)     = S23*(X20(3) - X20(2)) + X20T(3 - KL)
         EQU(EN14)  = S23*(X14(3) - X14(2)) + X14T(3 - KL)
         EQU(EMG24) = S23*(X24(3) - X24(2)) + X24T(3 - KL)
         EQU(ESUMX) = EQU(1)+EQU(2)+EQU(3)+EQU(4)+EQU(5)+EQU(EN14)+EQU(EMG24);

         ! Angular momentum transport
         S23 = KQ*0.5D0*(SGAM(2)*SI(2) + SGAM(3)*SI(3))
         ! FIXME: add other mixing coefficients for central boundary condition
         XEQUS(EAMT) = S23*(OMEGA(3) - OMEGA(2)) + OMEGAT(3 - KL)
      END IF
      
C first-order difference equations at interior points.
c KL=0, KQ=1 means that XX(3) is deeper in star than XX(2); KL=1, KQ=-1
c means nearer surface. XXK is sort-of dXX/dK, except that its sign is
c *independent* of KL, KQ, so a factor of KQ is put in.
      IF ( 2.LE.JK .AND. JK.LE.KH ) THEN
         WTA = 0.5D0*KQ
         WTB = 0.5D0*(WT(2) + WT(3))
         WTC = WTA*WTB/(1.0D0 + WTB)
         !IF (WTC < 0.3) WTC = 0.0d0
         !IF (WTC > 0.7) WTC = 1.0d0
         WTD = KQ - WTC
         EQU(6) = VP(3) - VP(2) - WTC*VPK(3 - KL) - WTD*VPK(2 + KL)
         EQU(7) = VR(3) - VR(2) - WTD*VRK(3 - KL) - WTC*VRK(2 + KL)
         EQU(8) = VT(3) - VT(2) - WTC*VTK(3 - KL) - WTD*VTK(2 + KL)
c Luminosity equn: MT gives advective term, DLRK heat transfer in contact
         !WTC = (0.5 + 0.5D0*(1.0D0 - EGEN_SMOOTH)) * KQ
         !WTD = KQ - WTC
         II = 2*JSTAR - 3
         EQU(9) = L(3) - L(2) - WTD*LK(3 - KL) - WTC*LK(2 + KL) 
     :        - LQ(2)*PSTV(KQ*MT(2),0.0D0) + LQ(3)*PSTV(-KQ*MT(3),0.0D0)
     :        + II*KQ*DLRK(2 + KL)
         EQU(10) = VM(3) - VM(2) - WTA*(VMK(3) + VMK(2))
         EQU(11) = VI(3) - VI(2) - WTA*(VIK(3) + VIK(2))
         EQU(12) = PHI(3) - PHI(2) - WTA*(PHIK(3) + PHIK(2))
         EQU(19) = XI(3) - XI(2) - WTA*(XIK(3) + XIK(2))
         XEQUS(ETAM) = XVARS(3, FX_AM) - XVARS(2, FX_AM) - WTA*(XVARS(3, FX_AMK) + XVARS(2, FX_AMK))
      END IF
C surface boundary conditions for first-order equations and `eigenvalues'
      IF ( JK == 1 ) THEN
         EQU(6) = BCM(3)
         EQU(7) = BCP(3)
         EQU(8) = BCT(3)
         EQU(9) = BCF(3)
         EQU(10) = BCS(3)
         EQU(11) = BCPH(3)
         EQU(17) = BCA(3)
         EQU(18) = BCE(3)
         EQU(20) = BCMB(3)
      END IF
C central boundary conditions for first-order equations
      IF ( JK == KH + 1 ) THEN
         EQU(6) = VM(3)    
         EQU(7) = L(3)     
         EQU(8) = VR(3)   
         EQU(9) = VI(3) 
         EQU(19) = XI(3)
         XEQUS(ETAM) = XVARS(3, FX_AM)
      END IF
      CALL NAMES2 ( JSTAR, 2 )
      IF ( KTW == JSTAR ) RETURN
      JSTAR = 2
      GO TO 1
      END

! JBE = 1: Copy variables from X in INF to NAMEIN
! JBE = 2: Copy variables from NAMEOUT to Y in INF
      SUBROUTINE NAMES1 ( JSTAR, JBE )
      USE MESH
      IMPLICIT REAL*8 (A-H, L-Z)
      INTEGER :: I, II
      COMMON /INF   / X(NVAR), DX(NVAR), Y(NFUNC), Z(NVAR,NFUNC)
      COMMON /NAMEIN/ VAR(NVSTAR), BINVAR(NVBIN),
     &               DVAR(NVSTAR), DBINVAR(NVBIN),
     &               V(3*NFUNC-2*(NVSTAR+NVBIN))
      COMMON /NAMOUT/ W(NFUNC)
      COMMON /XIN1  / XVS(NXVSTAR), XVB(NXVBIN), DXVS(NXVSTAR), DXVB(NXVBIN)
      COMMON /XOUT1 / XWS(NXFSTAR), XWB(NXFBIN)
      DIMENSION IM(4), IN(4), IP(4), IQ(4)
c nv1=16 each, nv2=8 both; nf1=42 each, nf2=16 both: 16/8/16, 42/16/42
      DATA IM /24, 16, 42, 42/ 
      DATA IN / 0, 24,  0, 58/ 
      DATA IP /24, 24, 43, 43/ 
      DATA IQ / 0, 24, 58, 58/ 
      II = 2*JBE + JSTAR - 2
      IF ( JBE == 1 ) THEN          ! X, DX->NAMEIN
         ! Copy normal variables (up to NSFUNC)
         !V( 1 : IM(II) ) = X( IN(II)+1 : IN(II)+IM(II) )
         I = (JSTAR-1)*24
         VAR(1:16) = X(1+I:16+I)
         BINVAR(1:8) = X(17:24)
         !V( IP(II)+1 : IP(II)+IM(II) ) = DX( IQ(II)+1 : IQ(II)+IM(II) )
         DVAR(1:16) = DX(1+I:16+I)
         DBINVAR(1:8) = DX(17:24)
         
         ! Copy extra variables, star 1 or star 2
         I = NSVAR+(JSTAR-1)*NXVSTAR
         VAR(17:NVSTAR) = X(I+1:I+NXVSTAR)
         DVAR(17:NVSTAR) = DX(I+1:I+NXVSTAR)

         ! Extra binary variables
         I = NSVAR+2*NXVSTAR
         BINVAR(9:NVBIN) = X(I+1:I+NXVBIN)
         DBINVAR(9:NVBIN) = DX(I+1:I+NXVBIN)

         ! Copy extra variables, star 1 or star 2
         I = NSVAR+(JSTAR-1)*NXVSTAR
         XVS(1:NXVSTAR) = X(I+1:I+NXVSTAR)
         DXVS(1:NXVSTAR) = DX(I+1:I+NXVSTAR)

         ! Extra binary variables
         I = NSVAR+2*NXVSTAR
         XVB(1:NXVBIN) = X(I+1:I+NXVBIN)
         DXVB(1:NXVBIN) = DX(I+1:I+NXVBIN)
      ELSE IF ( JBE == 2 ) THEN     ! NAMEOUT->Y
         Y(IN(II)+1 : IN(II)+IM(II) ) = W( 1:IM(II) )
         Y(IP(II) : IQ(II) ) = W( IP(II) : IQ(II) )
         
         ! Copy extra functions, star 1 or star 2
         ! The blocklayout is different here: the star blocks come first.
         !  followed by the binary parameters block.
         I = NSFUNC+(JSTAR-1)*NXFSTAR
         Y(I+1:I+NXFSTAR) = XWS(1:NXFSTAR)
         I = NSFUNC+2*NXFSTAR
         Y(I+1:I+NXFBIN) = XWB(1:NXFBIN)
      END IF
      RETURN
      END
      
! JBE = 1: Copy variables from X in INE to NAMEIN
! JBE = 2: Copy variables from NAMEOUT (EQUations) to Y in INE
      SUBROUTINE NAMES2 ( JSTAR, JBE )
      USE MESH
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON /INE   / X(3*NFUNC), DX(3*NVAR*NFUNC), Y(NEQ), Z(3*NVAR*NVAR)
      COMMON /NAMEIN/ V(3*NFUNC)
      COMMON /NAMOUT/ W(NFUNC)
      COMMON /XIN2  / XVS(3, NXFSTAR), XVB(3, NXFBIN)
      COMMON /XOUT2 / XWS(NXVSTAR), XWB(NXVBIN)
      DIMENSION IM(4), IN(4), IP(4), IQ(4)
c nf1=42*3each, nf2=8*3 both; ne1=16 each, ne2=8 both: 126/48/126, 16/8/16
      DATA IM, IN /126, 126, 16, 16,   0, 174,  0, 24/
      DATA IP, IQ /127, 127, 17, 17, 174, 174, 24, 24/ 
      II = 2*JBE + JSTAR - 2
      IF ( JBE == 1 ) THEN          ! X->NAMEIN         
         V(1:IM(II)) = X(IN(II)+1 : IN(II)+IM(II))   ! Stars
         V(IP(II):IQ(II)) = X(IP(II):IQ(II))       ! Binary
         
         ! Extra functions
         ! Stars
         DO I=1, NXFSTAR
            XVS(1,I) = X(3*NSFUNC + (JSTAR-1)*3*NXFSTAR + (I-1)*3 + 1)
            XVS(2,I) = X(3*NSFUNC + (JSTAR-1)*3*NXFSTAR + (I-1)*3 + 2)
            XVS(3,I) = X(3*NSFUNC + (JSTAR-1)*3*NXFSTAR + (I-1)*3 + 3)
         END DO
         ! Binary
         DO I=1, NXFBIN
            XVB(1,I) = X(3*NSFUNC + 3*2*NXFSTAR + (I-1)*3 + 1)
            XVB(2,I) = X(3*NSFUNC + 3*2*NXFSTAR + (I-1)*3 + 2)
            XVB(3,I) = X(3*NSFUNC + 3*2*NXFSTAR + (I-1)*3 + 3)
         END DO
      ELSE  IF ( JBE == 2 ) THEN    ! NAMEOUT->Y
         Y(IN(II)+1:IN(II)+IM(II)) = W(1:IM(II))
         Y(IP(II):IQ(II)) = W(IP(II):IQ(II))
         
         ! Extra equations
         I = NSEQ+(JSTAR-1)*NXESTAR
         Y(I+1:I+NXESTAR) = XWS(1:NXESTAR)
         I = NSEQ+2*NXESTAR
         Y(I+1:I+NXEBIN) = XWB(1:NXEBIN)
      END IF
      RETURN
      END

   
! Remesh:
!  Adjust the mesh spacing of a model
!  KH2 - new number of meshpoints
!        KH in the common block holds the current number of meshpoints
!   OA - New orbital angular momentum
!   P1 - New rotation period
!  ECC - New eccentricity
      SUBROUTINE REMESH ( KH2, JCH, BM, TM, P1, ECC, OA, JSTAR, JF )
c JCH = 1, 2: initialises some new variables, depending on JF. JCH=3: also
c constructs new mesh spacing. JCH=4: also initialises compos. to uniformity
      USE MESH
      USE MESH_ENC;
      USE EXTRA_ELEMENTS;
      USE CONSTANTS
      USE SETTINGS
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON H(NVAR,NM), DH(NVAR,NM), EP(3), KH, KTW, ID(260)
      COMMON /STORE / HPR(NVAR,NM), MS(60025)
      COMMON /TVBLES/ ZQ(12), MC(2), OM(230), JHOLD(3)
      COMMON /VBLES / LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QMU, WL, 
     : WCV, HP, WT, PHIMU, GMR, SEP, M3, PX(NPX*(NM+2)), QA(NM)
      COMMON /STAT2 / WU(3), P, RHO, WS(36), EX, ENX, WMU, DELTA, PHIE,
     & EXT, FKT, FKR, PRANDTL
      COMMON /INF   / VAR(NVAR), DVAR(NVAR), 
     : FILLUP_COMMONBLOCK_INF(NFUNC+NVAR*NFUNC)
      COMMON /ATDATA/ CH2(4), CHI(26,9), COM(27), CAN(9), CBN(9), KZN(9)
      COMMON /ABUND / XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, XW(14)
      DOUBLE PRECISION :: X1AC, X4AC, X12AC, X14AC, X16AC, X20AC, X24AC
      DOUBLE PRECISION :: XAC(7, 2)
      COMMON /ACCRET/ X1AC, X4AC, X12AC, X14AC, X16AC, X20AC, X24AC, XAC
      DOUBLE PRECISION :: NEW_H(NVAR,NM), NEW_DH(NVAR,NM)
      DOUBLE PRECISION :: EX_MAX

* Set initial composition.
* The composition variables are NOT the actual mass fractions if we 
* use non-integer atomic masses, so we have to compute what they are
* The composition variables used in the code are baryon number fractions
* We calculate this even if we don't want to do a ZAMS run because we need to
* know the baryon number densities of Fe, Si and Mg.
      CHE = 1.0D0 - CH - CZS
      CN = 1.0D0 - CC - CO - CNE - CMG - CSI - CFE
      XH = CH*CBN(1)/CAN(1)
      XHE = CHE*CBN(2)/CAN(2)
      XC = CC*CZS*CBN(3)/CAN(3)
      XN = CN*CZS*CBN(4)/CAN(4)
      XO = CO*CZS*CBN(5)/CAN(5)
      XNE = CNE*CZS*CBN(6)/CAN(6)
      XMG = CMG*CZS
      XSI = CSI*MAX(CZS, 1.0D-4)
      XFE = CFE*MAX(CZS, 1.0D-4)
      VMA = XH + XHE + XC + XN + XO + XNE + XMG + XSI + XFE
      XFE = XFE / VMA
      XSI = XSI / VMA
      XMG = XMG / VMA

      ! Mg24 is never stored in the input model!
      XMG = MAX(1.0 - (XH+XHE+XC+XN+XO+XNE)/VMA-(XSI+XFE), 0.0D0);
      IF (USE_MG24_EQN) THEN
         !H(21,1:KH) = XMG
         DO IK=1, KH
            H(NMG24,IK) = MAX(0.0D0, 1.0 - (H(5,IK) + H(9, IK) + H(10, IK) + H(3, IK) + H(11, IK) + H(16, IK) + XFE + XSI))
         END DO
      END IF
      IF ( JCH >= 4 ) THEN
         H( 5,1:KH) = XH/VMA
         H( 3,1:KH) = XO/VMA
         H( 9,1:KH) = XHE/VMA
         H(10,1:KH) = XC/VMA
         H(11,1:KH) = XNE/VMA
         H(16,1:KH) = XN/VMA
      END IF

! Now we must also convert the abundances of the accreted material. If not
! set from init.dat, set from initial abundances.
      IF (X1AC < 0.0) X1AC = CH*CBN(1)/CAN(1)
      IF (X4AC < 0.0) X4AC = CHE*CBN(2)/CAN(2)
      IF (X12AC < 0.0) X12AC = CC*CZS*CBN(3)/CAN(3)
      IF (X14AC < 0.0) X14AC = CN*CZS*CBN(4)/CAN(4)
      IF (X16AC < 0.0) X16AC = CO*CZS*CBN(5)/CAN(5)
      IF (X20AC < 0.0) X20AC = CNE*CZS*CBN(6)/CAN(6)
      IF (X24AC < 0.0) X24AC = CMG*CZS*CBN(7)/CAN(7)
! make sure XH is 1-everything else and abundancies sum to 1
      X1AC = MAX(0.0D0, 1.0D0 -(X4AC+X12AC+X14AC+X16AC+X20AC+X24AC+CSI*CZS+CFE*CZS))

      XH = X1AC*CBN(1)/CAN(1)
      XHE = X4AC*CBN(2)/CAN(2)
      XC = X12AC*CBN(3)/CAN(3)
      XN = X14AC*CBN(4)/CAN(4)
      XO = X16AC*CBN(5)/CAN(5)
      XNE = X20AC*CBN(6)/CAN(6)
      XMG = X24AC*CBN(7)/CAN(7)
      !XSI = CSI*CZS
      !XFE = CFE*CZS
      
      VMA = XH + XHE + XC + XN + XO + XNE + XMG + XSI + XFE

      X1AC  = XH / VMA
      X4AC  = XHE / VMA
      X12AC = XC / VMA
      X14AC = XN / VMA
      X16AC = XO / VMA
      X20AC = XNE / VMA
      X24AC = XMG / VMA

      ! Initialise accretion abundances for both stars
      XAC(1, 1:2) = X1AC  
      XAC(2, 1:2) = X4AC  
      XAC(3, 1:2) = X12AC 
      XAC(4, 1:2) = X14AC 
      XAC(5, 1:2) = X16AC 
      XAC(6, 1:2) = X20AC 
      XAC(7, 1:2) = X24AC 

      QE = 0.0D0              ! Modified mesh spacing function
      MC(JSTAR) = TM          ! New mass after remesh
      VAR(:) = H(:, KH)
      DVAR(:) = 0.0D0
      CALL FUNCS1 ( KH, -2)
      HPC = DSQRT(P/(CG*RHO*RHO))
      MC(JSTAR) = 3.5D-33*RHO*HPC**3
      VD = TM/H(4, 1)
      H(4, 1:KH) = VD*H(4, 1:KH)                ! Scale mass
      H(6, 1:KH) = 1.0D0                        ! Gradient of mesh spacing
      IF ( JF == 0 ) H(12, 1:KH) = 1.0D0        ! Moment of Inertia (later)
      !H(13, 1:KH) = P1                          ! New rotation period
      H(17, 1:KH) = OA                          ! New orbital angular momentum
      H(18, 1:KH) = ECC                         ! New eccentricity
      H(20, 1:KH) = BM                          ! New total binary mass

!     Rotation rate/rotational frequency
      H(13, 1:KH) = 2.0*CPI/(P1 * CSDAY)
      DH(13, 1:KH) = -2.0*CPI/(P1**2 * CSDAY) * DH(13, 1:KH)

      EX_MAX = EX
      IF (EX_MAX < 1.0D-3) ENC_PARACHUTE = H(8, 1) / H(4, 1)
      DO IK = 1, KH
         VAR(:) = H(:, IK)
         CALL FUNCS1 ( IK, -2 )                 ! Calculate stuff
         QA(IK) = QQ                            ! Store mesh spacing function
C Preliminary integration for M.I. and potential
         IF ( .NOT. ( IK == 1 .OR. JF == 2 ) ) THEN
            H(12, IK) = H(12, IK - 1) + R*R*M3/(DABS(QMU))
            H(14, IK) = H(14, IK - 1) + PHIMU/DABS(QMU)
         END IF
         EX_MAX = MAX(EX, 0.0D0)
c Integration for L-dependent mesh spacing.
         !QE = QE + CT(2)*EX*HT(4, IK)*1.0D-8
         !QA(IK) = QA(IK) + QE                   ! Update meshspacing
      END DO
! Find values of mesh spacing function at the external points
! Needed to calculate the new mesh, where the meshspacing gradient is
! constant.
      Q1 = QA(1)
      Q2 = (KH2 - 1.0D0)/(QA(KH) - QA(1))
      IF ( JCH >= 3 ) THEN
c If required, redistribute the mesh-points, either by changing
c their number or by changing the mesh-distribution function
         DO IK = 1, KH
            VX = H(5, IK) + 1.0D-10             ! Fudge hydrogen fraction
            H(5, IK) = DLOG(VX)
            QA(IK) = (QA(IK) - Q1)*Q2 + 1.0D0   ! Adjust meshspacing
         END DO
         IH = 1
         DO IK = 1, KH2
            DK = 0.0D0
            IF ( IK == KH2 ) IH = KH
            IF ( IK /= 1 .AND. IK /= KH2 ) THEN
               ! Find the proper value for the meshspacing function at
               ! this meshpoint
               DO I = 1, 50
                  ! Sanity check: abort if we're running out of the mesh
                  ! boundary
                  IF ( IH+1 > KH) THEN
                     WRITE (0, *) 
     &                'REMESH running outside mesh boundary, aborting'
                     STOP
                  END IF
                  IF ( IK >= QA(IH + 1) ) IH = IH + 1
                  IF ( IK < QA(IH + 1) ) EXIT   ! Break loop
               END DO
               DK = (IK - QA(IH))/(QA(IH + 1) - QA(IH))
            END IF
            ! Linear interpolation for new H and DH
            NEW_H(:, IK) = H(:, IH) + DK*(H(:, IH + 1) - H(:, IH))
            NEW_DH(:, IK) = DH(:, IH) + DK*(DH(:, IH + 1) - DH(:, IH))
         END DO
         Q2 = Q2 * (KH - 1.0) / (KH2 - 1.0)
         KH = KH2
         DO IK = 1, KH
            ! Un-fudge the hydrogen abundance
            NEW_H(5, IK) = DEXP(NEW_H(5, IK)) - 1.0D-10
            IF ( NEW_H(5, IK) < 1.0D-5 ) NEW_H(5, IK) = 0.0D0
            H(:, IK) = NEW_H(:, IK)
            DH(:, IK) = NEW_DH(:, IK)
         END DO
      END IF
      QK = 1.0D0/Q2                             ! Gradient of mesh spacing
      SI = H(12, KH)                            ! Potential at the surface
C Some new variables that may not have been read in from stored model
      DO IK = 1, KH
         H(6, IK) = QK                          ! Gradient of mesh spacing
         IF (JF /= 2) THEN
            H(12, IK) = (SI - H(12, IK))*DABS(QK)  ! Moment of inertia of interior
            H(14, IK) = - GMR - H(14, IK)*DABS(QK) ! Gravitational potential
            H(15, IK) = - GMR                      ! Potential at the stellar surface
         END IF
      END DO
      IF (EX_MAX < 1.0D-3) ENC_PARACHUTE = 1.1*H(8, 1) / H(4, 1)
      H(NTAM, 1:KH) = H(12, 1:KH)*H(13, 1);
      RETURN
      END SUBROUTINE


      SUBROUTINE POTENT ( QQ, DPHI )
c Solve for dimensionless L1, L2 potentials
      IMPLICIT REAL*8 (A-H, L-Z)
      Q = QQ
      IF ( Q < 1.0D0 ) Q = 1.0D0/QQ
      Q5 = 1.0D0 + Q
      XL1 = (Q/(3.0D0*Q5))**0.333333333D0
      XL2 = 2.0D0 - Q*XL1/Q5
      Q3 = 3.0D0 + Q
      Q2 = Q + Q
      Q4 = 3.0D0 + Q2
      Q1 = 2.0D0 + Q
      Q6 = 1.0D0/Q5
      DO 1 IJ = 1, 4
c Newton-Raphson iteration for L1, L2 points
      XL1 = XL1 + 
     : (Q - XL1*(Q2 - XL1*(Q - XL1*(Q3 - XL1*(Q4 - XL1*Q5)))))
     : /(Q2 - XL1*(Q2 - XL1*(3.0D0*Q3 - XL1*(4.0D0*Q4 - XL1*5.0D0*Q5))))
 1    XL2 = XL2 + (Q - XL2*(Q2 - XL2*(Q1 - XL2*(Q3 - XL2*(Q4 - XL2*Q5)))
     :   ))/(Q2 - XL2*(2.0D0*Q1 - XL2*(3.0D0*Q3 - XL2*(4.0D0*Q4 
     :   - XL2*5.0D0*Q5))))
      DPHI = Q*Q6/XL1 + Q6/(1.0D0 - XL1) + 0.5D0*(XL1 - Q6)**2
     :    - (Q*Q6/XL2 + Q6/(XL2 - 1.0D0) + 0.5D0*(XL2 - Q6)**2)
      RETURN
      END

      SUBROUTINE DJNVDH ( TT, LL, MDOT )
c Mass loss rate for luminous star, de Jager et al (1988)
c Input parameters:
c     TT: log10 effective temperature/K
c     LL: luminosity
c Output parameter:
c     MDOT: Mass loss rate in solar masses per year
c ------------------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H, L-Z)
      DIMENSION A(6,6)
      DATA A/ 
     :  6.34916D0,3.41678D0,-1.08683D0, 0.13095D0,0.22427D0,0.11968D0, !J=1
     : -5.04240D0,0.15629D0, 0.41952D0,-0.09825D0,0.46591D0,0.0D0,     !J=2
     : -0.83426D0,2.96244D0,-1.37272D0, 0.13025D0,2*0.0D0,             !J=3
     : -1.13925D0,0.33659D0,-1.07493D0, 3*0.0D0,                       !J=4
     : -0.12202D0,0.57576D0, 10*0.0D0/                                 !J=5,6
      CHEB(X, I) = DCOS(I*DACOS(MIN(MAX(X, -1.0D0), 1.0D0)))
      MDOT = 0.0D0
      IF ( TT < 3.301D0 .OR. LL < 2.501D0 ) RETURN
      IF ( TT > 4.799D0 .OR. LL > 6.699D0 ) RETURN
      DO 2 I = 0, 5
      DO 2 J = 0, 5 - I
 2    MDOT = MDOT + A(I + 1, J + 1)*CHEB((TT - 4.05D0)/0.75D0, I)
     :                             *CHEB((LL - 4.6D0)/2.1D0, J)
      MDOT = 10.0D0**(-MDOT)
      RETURN
      END

      FUNCTION CALC_MDOT_JVINK(Mass, Lum, Teff)
c ------------------------------------------------------------------------------
c J. Vinks Mass loss recipe for massive stars including a metallicity dependence
c Based on Vink et al (99, 00, 01) see http://www.astro.keele.ac.uk/~jsv/
c Implemented in STARS: SdM July 2006
c ------------------------------------------------------------------------------
c Input parameters:
c     CZS: metallicity (0.02 for solar)m taken from settings
c     Teff: effective temperature
c     logL: luminosity
c     Mass: stellar mass
c Output parameter:
c     mdot_vink: mass loss predicted by J Vink, in solar masses per year
c ------------------------------------------------------------------------------
      USE CONSTANTS
      USE SETTINGS
      IMPLICIT NONE  
      DOUBLE PRECISION :: CALC_MDOT_JVINK
! constants   
      DOUBLE PRECISION :: SMALL = 1.0D-40   
      DOUBLE PRECISION :: SIGMAE = 0.325
! input/output variabels
      DOUBLE PRECISION, INTENT(IN) :: LUM, MASS, TEFF ![Lsun],[Msun],[[K]
! Coefficients below and above the first bistability jump
      DOUBLE PRECISION :: COEFF_BELOW(7) = (/-6.688, 2.210, -1.339, -1.601, 1.07, 0., 0.85/)
      DOUBLE PRECISION :: COEFF_ABOVE(7) = (/-6.697, 2.194, -1.313, -1.226, 0.933, -10.92, 0.85/)
      DOUBLE PRECISION :: TEFF_SMOOTH = 1000.0
! parameters
      DOUBLE PRECISION :: COEFF(7), RATIO !=vesc/vinf 
! derived variables
      DOUBLE PRECISION :: GAMMA
      DOUBLE PRECISION :: CHAR_RHO, TEFF_JUMP1, TEFF_JUMP2
      DOUBLE PRECISION :: LOG_TEFF, lOG_MDOT, X, MDOTJ
! Use second bistability jump?
      LOGICAL :: USE_SECOND_JUMP = .TRUE.

! TODO: SIGMAE should be calculated as in Lamers&Leitherer 1993 instead of being
! constant, although this is probably a minor correction.
      GAMMA = 7.66D-5 *SIGMAE * LUM/MASS
     
c Determine postions of bistability jumps
      CHAR_RHO   =  -14.94 + ( 3.1875 * Gamma   ) + (0.85 * CLOGZ) !Eq 23, Vink (2001)
      TEFF_JUMP1 = ( 61.2  +   2.59  * CHAR_RHO ) *  1.0D3         !Eq 15, Vink (2001)
      TEFF_JUMP2 = ( 1.0D2 +   6.0D0 * CHAR_RHO ) *  1.0D3         !Eq  6, Vink (2001)
!      IF(Teffjump2 > Teffjump1) PRINT *, 
!     :     "Warning: invalid Bistability jumps, stellar par. unrealis."

c Determine postion with respect to temperature jumps 
      IF(TEFF < TEFF_JUMP1) THEN
         RATIO = 1.3
         COEFF(:) = COEFF_BELOW(:)
         IF (TEFF < TEFF_JUMP2) THEN ! below both jumps
            !RATIO = 0.7
            !IF (USE_SECOND_JUMP) COEFF(1) = -5.990 

! Smooth out the size of the (second) bistability jump over ~ 1000 K
            IF ( abs(TEFF - TEFF_JUMP2)<TEFF_SMOOTH ) THEN
               X = ( TEFF_JUMP2+TEFF_SMOOTH - TEFF ) / (2.0*TEFF_SMOOTH)
               !X = 0.5*(1.0D0 - COS(X*CPI))   ! Cosine interpolation, is smoother
               !RATIO = X*0.7 + (1.0D0 - X)*1.3
               IF (USE_SECOND_JUMP) COEFF(1) = (1.0d0-X)*COEFF_BELOW(1) + X*(-5.990)
            END IF
         ENDIF
      ELSE                          ! above both jumps
         ratio = 2.6
         COEFF(:) = COEFF_ABOVE(:)
      ENDIF

! Smooth out the size of the (first) bistability jump over ~ 1000 K
      IF ( abs(TEFF - TEFF_JUMP1)<TEFF_SMOOTH ) THEN
         X = ( TEFF_JUMP1+TEFF_SMOOTH - TEFF ) / (2.0*TEFF_SMOOTH)
         ! X = 0.5*(1.0D0 - COS(X*CPI))   ! Cosine interpolation, is smoother
         RATIO = X*1.3 + (1.0D0 - X)*2.6
         COEFF(:) = X*COEFF_BELOW(:) + (1.0D0 - X)*COEFF_ABOVE(:)
      END IF

c get mass loss
      LOG_TEFF = LOG10( MAX(SMALL, TEFF  / 4.0D4) )
      LOG_MDOT = COEFF(1)
     &        + COEFF(2) *  LOG10( MAX(SMALL, LUM   / 1.0D5) )
     &        + COEFF(3) *  LOG10( MAX(SMALL, MASS  / 3.0D1) )
     &        + COEFF(4) *  LOG10( MAX(SMALL, RATIO / 2.0D0) )
     &        + COEFF(5) *  LOG_TEFF
     &        + COEFF(6) * LOG_TEFF**2
     &        + COEFF(7) *  CLOGZ
    
      CALC_MDOT_JVINK = 10**lOG_MDOT
      IF (TEFF < TEFF_JUMP2) THEN
         CALL DJNVDH ( DLOG10(TEFF), DLOG10(LUM), MDOTJ )
         CALC_MDOT_JVINK = MAX(CALC_MDOT_JVINK, MDOTJ)
      END IF
      END
c --------------------------------------------------------------END JVINK MDOT

      SUBROUTINE STATEL ( JI, AF, AT, JSTAR )
c Either compute new EoS values, or recover old values from store
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON /STAT2 / VF(50)
      COMMON /STAT4 / SF(50, 2)
      DIMENSION IFN(43)
      DATA IFN/2, 2, 0, 3*2, 1, 2, 3*1, 3*2, 4*1, 2, 8*1, 
     :                  3*2, 1, 2, 3*1, 3*2, 4*1, 2 /
      ! IFN contains entries for the independent variables: a 2 is listed
      !  if a variation in the corresponding variable requires that the
      !  equation of state is reevaluated, a 1 means that this is not
      !  nescessary and the cached result can be used.
      ! JI=0,-1,-2 are `dummy variables'. 0 means that the EOS is calculated
      !  and stored, -1 and -2 indicate the request came from
      !  printb or remesh respectively.
      !  Don't compute eos for *1 if *2 has changed and vice versa.
      IF (JI<=0 .or. (JSTAR==1 .and. JI<=16).or.(JSTAR==2 .and. JI>=24)) then
         I = IFN(JI + 3)
         IF ( I /= 1 ) THEN
            CALL STATEF ( AF, AT )
            CALL NUCRAT ( AT )
            IF ( I > 0 ) RETURN
            SF(:, JSTAR) = VF(:)
            RETURN
         END IF
      ENDIF
      VF(:) = SF(:, JSTAR)
      RETURN
      END

      SUBROUTINE CHECKS
c Set some very small, or negative, values to zero
c In particular, this deals with composition variables
      USE MESH
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, KW(260)
      DIMENSION IX(20)
      DATA IX/3, 5, 9, 10, 11, 16, 18, 27, 29, 33, 34, 35, 40, 7*3/
! Composition variables - should not be negative!
      DO IK = 1, KH
         DO I = 1, 6
            IJ = IX(I)
            IF ( H(IJ, IK) + DH(IJ, IK) < 0.0 ) THEN
               DH(IJ, IK) = -H(IJ, IK)
               H(IJ, IK) = 0.0D0
            END IF
         END DO
      END DO
! Other variables
      DO IK = 1, KH
         DO I = 7, 13
            IJ = IX(I)
            IF ( H(IJ, IK) + DH(IJ, IK) <= 1.0D-12 ) THEN
               H(IJ, IK) = 0.0D0
               DH(IJ, IK) = 0.0D0
            END IF
         END DO
      END DO
      H(4, KH) = 0.0D0
      H(8, KH) = 0.0D0
      DH(4, KH) = 0.0D0
      DH(8, KH) = 0.0D0
      RETURN
      END

! ------------------------------------------------------------------------------
!  SCALE_EQN
!   Calculate the expected scaling for equations and (surface) boundary
!   conditions, based on the current values in the H() array and the
!   current timestep.
!   Calls FUNCS1() to compute quantities that enter the equations
! ------------------------------------------------------------------------------
!  Input:
!   H(:)          -  Current values of independent variables (taken from
!                    the nameless COMMON block)
!   DT            -  Current timestep, in s (taken from COMMON TVBLES)
!  Output:
!   EQN_SCALE(:)  -  Expected scaling for equations
!   SBC_SCALE(:)  -  Expected scaling for surface boundary conditions
! ------------------------------------------------------------------------------
      SUBROUTINE SCALE_EQN ( EQN_SCALE, SBC_SCALE )
      USE MESH
      USE EXTRA_ELEMENTS
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT) :: EQN_SCALE(1:NVAR), SBC_SCALE(1:NVAR)
      
!     Unnamed COMMON block: current values of variables
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, KW(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW

!     COMMON block INF: pass input to FUNCS1
      DOUBLE PRECISION :: Q(NVAR), QD(NVAR), DUMMY_INF(NFUNC+NVAR*NFUNC)
      COMMON /INF   / Q, QD, DUMMY_INF

!     COMMON block TVBLES: get the current timestep DT
      DOUBLE PRECISION :: DT, ZQ(243)
      INTEGER :: JHOLD, JM2, JM1
      COMMON /TVBLES/ DT, ZQ, JHOLD, JM2, JM1

!     COMMON block NAMOUT: get return values from FUNCS1
      DOUBLE PRECISION ::  BCP, BCT, VP, VPK, VR, VRK, VT, VTK, VL, LK,
     &     LQ, MT, VM, VMK, SG, WT, X1, X1T, X16, X16T, X4, X4T, X12,
     &     X12T, X20, X20T, BCM, VI, VIK, VPHI, PHIK, BCF, BCS, BCPH,
     &     X14, X14T, AVMU, SGTH, OMEGA, OMEGAT, SGAM, SI, BCA, BCE, XIM,
     &     XIK, DLRK, BCMB, X24, X24T, MENC, MENCT, MEA, MEAT, MET,
     &     METT, FILLUP_COMMONBLOCK_NAMOUT(2 + NSFSTAR + NXFUNC  )
      COMMON /NAMOUT/ BCP, BCT, VP, VPK, VR, VRK, VT, VTK, VL, LK,
     &     LQ, MT, VM, VMK, SG, WT, X1, X1T, X16, X16T, X4, X4T, X12,
     &     X12T, X20, X20T, BCM, VI, VIK, VPHI, PHIK, BCF, BCS, BCPH,
     &     X14, X14T, AVMU, SGTH, OMEGA, OMEGAT, SGAM, SI, BCA, BCE, XIM, 
     &     XIK, DLRK, BCMB, X24, X24T, MENC, MENCT, MEA, MEAT, MET,
     &     METT, FILLUP_COMMONBLOCK_NAMOUT

!     Local variables
      INTEGER :: IK, I
      INTEGER, PARAMETER :: JSTAR = 1  ! Only do single star for now

!     Initialise with default settings: all equations have a natural value of 1
      EQN_SCALE(:) = 0.0d0
      SBC_SCALE(:) = 1.0d0

      DO IK=KH, 1, -1
!        Copy current value of variables: single star
         QD(1:16) = DH(1 + 24*(JSTAR-1):16 + 24*(JSTAR-1), IK)
         Q(1:16) = H(1 + 24*(JSTAR-1):16 + 24*(JSTAR-1), IK) + QD(1:16)
!        Binary parameters
         QD(17:24) = DH(17:24, IK)
         Q(17:24) = H(17:24, IK) + QD(17:24)
!        Compute function values
         CALL FUNCS1 ( IK, -1 )

!        Composition equations
         EQN_SCALE(1) = MAX(EQN_SCALE(1), DABS(X1T)) 
         EQN_SCALE(2) = MAX(EQN_SCALE(2), DABS(X16T)) 
         EQN_SCALE(3) = MAX(EQN_SCALE(3), DABS(X4T)) 
         EQN_SCALE(4) = MAX(EQN_SCALE(4), DABS(X12T)) 
         EQN_SCALE(5) = MAX(EQN_SCALE(5), DABS(X20T)) 
         EQN_SCALE(EN14) = MAX(EQN_SCALE(EN14), DABS(X14T))
         EQN_SCALE(EMG24) = MAX(EQN_SCALE(EMG24), DABS(X24T))

!        Structure equations
         EQN_SCALE(6) = MAX(EQN_SCALE(6), DABS(VPK))
         EQN_SCALE(7) = MAX(EQN_SCALE(7), DABS(VRK))
         EQN_SCALE(8) = MAX(EQN_SCALE(8), DABS(VTK))
         EQN_SCALE(9) = MAX(EQN_SCALE(9), DABS(LK - LQ*MT))
         EQN_SCALE(10) = MAX(EQN_SCALE(10), DABS(VMK))
         EQN_SCALE(11) = MAX(EQN_SCALE(11), DABS(VIK))
         EQN_SCALE(12) = MAX(EQN_SCALE(12), DABS(PHIK))
      ENDDO
      EQN_SCALE(1) = 1.0d0
      EQN_SCALE(2) = 1.0d0
      EQN_SCALE(3) = 1.0d0
      EQN_SCALE(4) = 1.0d0
      EQN_SCALE(5) = 1.0d0
      EQN_SCALE(EN14) = 1.0d0
      EQN_SCALE(EMG24) = 1.0d0
      EQN_SCALE(ESUMX) = MAX(EQN_SCALE(1), EQN_SCALE(2), EQN_SCALE(3),
     &                       EQN_SCALE(4), EQN_SCALE(5), EQN_SCALE(EN14),
     &                       EQN_SCALE(EMG24))

      DO I = 1, NVAR
         IF (EQN_SCALE(I) == 0.0d0) EQN_SCALE(I) = 1.0d0
      END DO
      RETURN
      END SUBROUTINE
