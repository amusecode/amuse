!     For some variables it is inconvenient to calculate them implicitly, so
!     we calculate them explicitly, between timesteps/iterations
      module explicit_functions
      use mesh

      integer, parameter :: EXPLV_GRADMU = 1
      integer, parameter :: EXPLV_AVMU   = 2
      integer, parameter :: EXPLV_LOGP   = 3
      integer, parameter :: NUM_EXPLV = 4

      double precision :: expl_var(NM, NUM_EXPLV, 2)

      end module

      MODULE PRINTB_GLOBAL_VARIABLES
C     Module for exchanging variables between subroutines within PRINTB only

      DOUBLE PRECISION ::  DMTR, PERC, RCZ, DRCZ, DMT, TET,
     &     DMSW, OMS, FR, FL, SDC, SDM, SDS, STC, STM, STS, PSURF, PCNTR, VK2,
     &     RAF, DPHI, EPHI, F1, DF, HTOT, AJ, DL, TKH, MEX(12), ECC


C     Meshpoint were maximum temperature is reached
      INTEGER :: IM_TMAX
      END MODULE

      SUBROUTINE CHECK_STOP_CONDITIONS ( JSTAR, JO, JCM )
      USE MESH
      USE SETTINGS
      USE CONSTANTS
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JSTAR
      INTEGER, INTENT(INOUT) :: JO, JCM
!     / /
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, KW(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW
!     /QUERY /
      DOUBLE PRECISION :: ML, QL, XL, UC(21)
      INTEGER ::  JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /QUERY / ML, QL, XL, UC, JMOD, JB, JNN, JTER, JOC, JKH
!     /VBLES/
      DOUBLE PRECISION :: LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,WL,
     &     WCV, HP, WT, PHIM, GMR, SEP, M3, PX(NPX), SX(NPX,NM+1),QA(NM)
      COMMON /VBLES / LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,
     & WL, WCV, HP, WT, PHIM, GMR, SEP, M3, PX, SX, QA
!     /TVBLES/
      DOUBLE PRECISION :: DT, ZQ(1), HSPN(2), RLF(2), ZET(2),
     &     XIT(2), AGE,BM, MC(2), OM, PER, SM, ENC, TC(2), TFR, T0, M0,
     &     MTA, OM0, OMTA,A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS,
     &     SE, TN(2), WMH,WMHE, MH, MHE, MCO, VMG, BE(2), LH, LHE, LC,
     &     LNU, LTH, MCB(8),MSB(6), RCB(8), TCT(8), QR(81), PPR(81)
      INTEGER :: JHOLD, JM2, JM1 
      COMMON /TVBLES/ DT, ZQ, HSPN, RLF, ZET, XIT, AGE, 
     & BM, MC, OM, PER, SM, ENC, TC, TFR, T0, M0, MTA, OM0, OMTA, 
     & A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS, SE, TN, WMH, 
     & WMHE, MH, MHE, MCO, VMG, BE, LH, LHE, LC, LNU, LTH, MCB, 
     & MSB, RCB, TCT, QR, PPR, JHOLD, JM2, JM1
!     /INF /
      DOUBLE PRECISION :: Q(NVAR), QD(NVAR), 
     &     FILLUP_COMMONBLOCK_INF(NFUNC+NVAR*NFUNC)
      COMMON /INF   / Q, QD, FILLUP_COMMONBLOCK_INF

      DOUBLE PRECISION :: SDC, DMT, TKH
      SDC = DLOG10(SX(3,2))      ! log10 central density
      DMT = SX(34, KH+1)         ! dM/dt
      TKH = 1.0D22*CG*Q(4)*Q(4)/(R*Q(8)*CSY)

! Conditions for terminating *1 or *2; the UC are consts from fort.23
      IF ( JB == 1 .AND. RLF(1) > UC(1) ) JO = 4                  !0.1D0
      IF ( AGE >= UC(2) ) JO = 5                  !2.0D10
      IF ( LC > UC(3) ) JO = 6                                    !1.0D2
      IF ( JB == 2 .AND. RLF(JSTAR) > UC(4) ) JO = 7              !0.0D0 
      IF ( LHE > UC(5) .AND. SDC > UC(6) .AND. MHE == 0.0D0 ) JO = 8
      IF ( MHE > UC(7) .AND. SDC > UC(8) ) JO = 9
      IF ( DABS(DMT) > UC(9)*SM/TKH ) JO = 10    !30.0D0
      JCM = 1
      IF ( SX(10, 2) == 0.0D0 ) JCM = 2
      IF ( SX(11, 2) == 0.0D0 ) JCM = 3
      ! Reduce accuracy for He depleted cores (don't want this/not needed)
      !IF ( SX(11,2) < UC(10) ) EPS = UC(11)     !0.15D0, 1.0D-4 !! fudge
      ! restore old value when XHe < 1e-6
      !IF ( SX(11,2) < 1.0D-6 .AND. EPS > 1.0D-6 ) EPS = 1.0D-6     
      IF ( SM+1.0D-6 >= UC(13)) CMI = 0.0D0
      IF ( SX(10,2) < UC(15) ) JO = 51                      ! Stop at TAMS
      IF ( PX(17) > UC(16) .AND. UC(16) > 0.0D0) JO = 52    ! Radius > limit
      IF ( MH > UC(14)) JO = 13        ! End of ZAHB construction
      IF ( MH < MHE .OR. MHE < MCO ) JO = 14
      END SUBROUTINE



! Adjust the weights in the MSF according to rules-of-thumb
! FIXME: two stars in a binary may want different MSFs, this is not
! currently taken into consideration.
      SUBROUTINE ADJUST_MESH_SPACING(  )
      USE MESH
      USE SETTINGS
      IMPLICIT NONE
!     /TVBLES/
      DOUBLE PRECISION :: DT, ZQ(1), HSPN(2), RLF(2), ZET(2),
     &     XIT(2), AGE,BM, MC(2), OM, PER, SM, ENC, TC(2), TFR, T0, M0,
     &     MTA, OM0, OMTA,A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS,
     &     SE, TN(2), WMH,WMHE, MH, MHE, MCO, VMG, BE(2), LH, LHE, LC,
     &     LNU, LTH, MCB(8),MSB(6), RCB(8), TCT(8), QR(81), PPR(81)
      INTEGER :: JHOLD, JM2, JM1 
      COMMON /TVBLES/ DT, ZQ, HSPN, RLF, ZET, XIT, AGE, 
     & BM, MC, OM, PER, SM, ENC, TC, TFR, T0, M0, MTA, OM0, OMTA, 
     & A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS, SE, TN, WMH, 
     & WMHE, MH, MHE, MCO, VMG, BE, LH, LHE, LC, LNU, LTH, MCB, 
     & MSB, RCB, TCT, QR, PPR, JHOLD, JM2, JM1

! Set defaults
      !CT(:) = INITIAL_CT(:)
! PPN phase, increase resolution in the outer regions to reduce numerical 
! diffusion
      !IF (MHE > 0.99*MH) CT(3) = CT(3) * (1.0 + 50.0*(MHE/MH - 0.99)) 
      !IF (MH > 0.95*SM .and. CT(3) < 4.0*INITIAL_CT(3)) CT(3) = CT(3) * 1.05
      END SUBROUTINE



! Adjust the timestep control parameter according to rules-of-thumb
      SUBROUTINE UPDATE_TIMESTEP_PARAMETERS( JSTAR )
      USE MESH
      USE SETTINGS
      USE PLOTVARIABLES
      IMPLICIT NONE      
      INTEGER, INTENT(IN) :: JSTAR
!     / /
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, KW(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW
!     /VBLES/
      DOUBLE PRECISION :: LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,WL,
     &     WCV, HP, WT, PHIM, GMR, SEP, M3, PX(NPX), SX(NPX,NM+1),QA(NM)
      COMMON /VBLES / LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,
     & WL, WCV, HP, WT, PHIM, GMR, SEP, M3, PX, SX, QA
!     /TVBLES/
      DOUBLE PRECISION :: DT, ZQ(1), HSPN(2), RLF(2), ZET(2),
     &     XIT(2), AGE,BM, MC(2), OM, PER, SM, ENC, TC(2), TFR, T0, M0,
     &     MTA, OM0, OMTA,A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS,
     &     SE, TN(2), WMH,WMHE, MH, MHE, MCO, VMG, BE(2), LH, LHE, LC,
     &     LNU, LTH, MCB(8),MSB(6), RCB(8), TCT(8), QR(81), PPR(81)
      INTEGER :: JHOLD, JM2, JM1 
      COMMON /TVBLES/ DT, ZQ, HSPN, RLF, ZET, XIT, AGE, 
     & BM, MC, OM, PER, SM, ENC, TC, TFR, T0, M0, MTA, OM0, OMTA, 
     & A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS, SE, TN, WMH, 
     & WMHE, MH, MHE, MCO, VMG, BE, LH, LHE, LC, LNU, LTH, MCB, 
     & MSB, RCB, TCT, QR, PPR, JHOLD, JM2, JM1
      DOUBLE PRECISION, SAVE :: RLF_PREV(2) = (/0.0, 0.0/)
      DOUBLE PRECISION, SAVE :: QCNV_PREV(2) = (/0.0, 0.0/)

      IF ( SX(10, KH+1) < CXB ) MH = SX(9, KH+1) 
      CDD = CDC(1)
! Allow larger or smaller time-increments, according to rules-of-thumb
      IF ( MHE <= MH ) THEN
! End of main sequence, reduce timestep to resolve hook
         IF (SX(10,2) < 5.0D-2 .AND. SX(10,2) > 1.0D-5) CDD = CDC(1)*CDC(6)
! Hertzsprung gap/end of core hydrogen burning
         IF (SX(10,2) < 1.0D-5 .AND. SX(11,2) > 1.0D-5
     &       .AND. MH < 0.22 .AND. SX(20,2) < 1.0) CDD = CDC(1)*CDC(7)
! First dredge up; can probably be extended t any dredgeup episode
         IF (MH > 0.12 .AND. MH < 0.4 .AND. QCNV>0.01 .AND. QCNV > QCNV_PREV(JSTAR)) CDD = CDC(1)*CDC(8)
! Core helium burning
         IF (SX(10,2) < 1.0D-5.AND.SX(11,2) < 0.95D0) CDD = CDC(1)*CDC(2)
! End of core helium burning
         IF ( SX(11,2) < 0.1D0 ) CDD = CDC(1)*CDC(3)
! Double shell burning/shells become close
         IF ( MHE > 0.75D0*MH ) CDD = CDC(1)*CDC(4)
! Reduce timestep when approaching Roche lobe overflow -SdM
         IF (RLF(JSTAR) > -2.0D-2 .AND. RLF(JSTAR) > RLF_PREV(JSTAR))
     &      CDD = MIN(CDD, CDC(1)*CDC(9))
! Keep timestep small if star is close to RLOF but shrinking
         IF (RLF(JSTAR) > -2.0D-2 .AND. RLF(JSTAR) < RLF_PREV(JSTAR))
     &      CDD = MIN(CDD, CDC(1)*CDC(10))
! Switch off Reimers-like wind for He stars
         IF ( MHE > 0.96D0*SX(9, KH+1) ) CMR = 0.0D0
      END IF
      RLF_PREV(JSTAR) = RLF(JSTAR)
      QCNV_PREV(JSTAR) = QCNV

! Dirty hacks to help get through He-flash; not all work well and none
! work well enough
      IF ( LHE > 1.0D-2 .AND. DLOG10(SX(3, 2)) > 5.3 ) CDD = CDD*0.1 
      IF ( LHE > 1.0 .AND. DLOG10(SX(3, 2)) > 5.3 ) CDD = CDD*0.5
      !IF ( LHE > 1.0D5 ) CDD = CDD*0.5
      !IF ( LHE > 1.0D7 ) EPS = 1.0D-4
      !IF ( LHE < 1.0D7 .AND. EPS>1.0D-6) EPS = 1.0D-6
      !IF ( LHE > 1.0D6 ) EPS = 1.0D-2
      !IF ( LHE < 1.0D6 .AND. EPS>1.0D-4) EPS = 1.0D-4
      !IF ( LHE > 1.0D0 ) WANTED_EPS = 1.0D-12
      !IF ( LHE > 2.0D6 ) LLUMI_SMOOTH = 0.0D0
      !IF ( LHE < 2.0D6 ) LLUMI_SMOOTH = 1.0D0
      END SUBROUTINE



      SUBROUTINE UPDATE_EXPLICIT_QUANTITIES( JO, JSTAR )
      USE MESH
      USE MESH_ENC
      USE FUDGE_CONTROL
      USE CONSTANTS
      USE SETTINGS
      USE EXPLICIT_FUNCTIONS
      USE FGB2HB_COMPOSITION
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: JO
      INTEGER, INTENT(IN) :: JSTAR
!     / /
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, KW(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW
!     /QUERY /
      DOUBLE PRECISION :: ML, QL, XL, UC(21)
      INTEGER ::  JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /QUERY / ML, QL, XL, UC, JMOD, JB, JNN, JTER, JOC, JKH
!     /TVBLES/
      DOUBLE PRECISION :: DT, ZQ(1), HSPN(2), RLF(2), ZET(2),
     &     XIT(2), AGE,BM, MC(2), OM, PER, SM, ENC, TC(2), TFR, T0, M0,
     &     MTA, OM0, OMTA,A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS,
     &     SE, TN(2), WMH,WMHE, MH, MHE, MCO, VMG, BE(2), LH, LHE, LC,
     &     LNU, LTH, MCB(8),MSB(6), RCB(8), TCT(8), QR(81), PPR(81)
      INTEGER :: JHOLD, JM2, JM1 
      COMMON /TVBLES/ DT, ZQ, HSPN, RLF, ZET, XIT, AGE, 
     : BM, MC, OM, PER, SM, ENC, TC, TFR, T0, M0, MTA, OM0, OMTA, 
     : A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS, SE, TN, WMH, 
     : WMHE, MH, MHE, MCO, VMG, BE, LH, LHE, LC, LNU, LTH, MCB, 
     : MSB, RCB, TCT, QR, PPR, JHOLD, JM2, JM1
      DOUBLE PRECISION :: W1, DLOGMU, DLOGP
      INTEGER :: ik

! Update the composition of accreted material, for instance during ZAHB
! construction. This does not matter for the core, but it may matter for
! the envelope or the hydrogen burning shell (according to Haili Hu).
      call update_accretion_abundance

! update constant extra energy source
      W1 = DT*CET/CSY
      ENC = ENC*(1.0D0 + W1)*CEA
      IF (CEA /= 0.0D0) ENC = ENC/(CEA + W1*ENC) 

      IF (JNN > 1) ENC_PARACHUTE = ENC_PARACHUTE*0.1
      IF (ENC_PARACHUTE < 1.0d-4) ENC_PARACHUTE = 0.0d0

      IF (ADJ_MEA) THEN
         ! Monitor convergence to target model.
         CALL CHECK_CONVERGENCE_TO_TARGET_STRUCTURE()
         IF ( CURR_DIFFSQR < BEST_DIFFSQR ) THEN
            BEST_DIFFSQR = CURR_DIFFSQR
            BEST_MOD = JMOD+1
            MUTANT_H(1:24, 1:KH) = H(1:24, 1:KH)
         END IF
         IF ( BEST_DIFFSQR<1.0D-4 ) JO = 53
      ENDIF
      
      IF ( EPS>1.0D-4 ) THEN
         EPS = MAX(EPS*0.8, 1.0D-4)
         print *, 'Set EPS to ', EPS
      END IF

! Adjust the timescale for the artificial energy term. This causes it to be
! turned on smoothly.
      IF (USEMENC .AND. IMPOSE_ENTROPY_FACTOR<1.0D0) THEN
         IMPOSE_ENTROPY_FACTOR = MIN(2.0D0*IMPOSE_ENTROPY_FACTOR, 1.0D0)
         !PRINT *, 'Set eart factor to ', IMPOSE_ENTROPY_FACTOR
      ENDIF
      
! Adjust the fudge-factor for the artificial composition adjustment. This causes
! the composition profile to be changed more gradually by turning the factor
! on smoothly.
      IF (ADJ_COMP .AND. IMPOSE_COMPOSITION_FACTOR<1.0D0) THEN
         IMPOSE_COMPOSITION_FACTOR = MIN(2.0D0*IMPOSE_COMPOSITION_FACTOR, 1.0D0);
         !print *, 'Set composition factor to ', IMPOSE_COMPOSITION_FACTOR
      END IF

! Explicit functions otherwise used in funcs1
! Calculate molecular weight gradient.
! Don't bother with boundary points (can they even matter?)
      expl_var(1, EXPLV_GRADMU, jstar) = 0.0d0
      expl_var(KH, EXPLV_GRADMU, jstar) = 0.0d0
      do ik=2, KH-1
         dlogmu = 0.5d0 *
     &      (expl_var(ik-1,EXPLV_AVMU,jstar) - expl_var(ik+1,EXPLV_AVMU,jstar))
         dlogmu = dlogmu / expl_var(ik, EXPLV_AVMU, jstar)
         dlogp = 0.5d0 *
     &      (expl_var(ik-1,EXPLV_LOGP,jstar) - expl_var(ik+1,EXPLV_LOGP,jstar))
         if (dlogp /= 0.0d0) expl_var(ik, EXPLV_GRADMU, jstar) = dlogmu / dlogp
      end do
      END SUBROUTINE



! FIND_STELLAR_TYPE:
!
! Determine the stellar type, more or less as in HPT00 (some of these
! will never be returned obviously).
! Return value:
!  0 - Main sequence star, Qenv > 0.1
!  1 - Main sequence star, other
!  2 - Hertzsprung gap
!  3 - RGB
!  4 - Core helium burning/horizontal branch
!  5 - Early AGB
!  6 - Thermal pulsing AGB
!  7 - Helium main sequence star / WR star
!  8 - Helium Hertzsprung gap star
!  9 - Helium giant
! 10 - Helium white dwarf
! 11 - CO white dwarf
! 12 - O/Ne white dwarf
! 13 - neutron star
! 14 - black hole
! 15 - massless remnant
      FUNCTION FIND_STELLAR_TYPE ( )
      USE MESH
      USE SETTINGS
      USE PLOTVARIABLES
      IMPLICIT NONE
      INTEGER :: FIND_STELLAR_TYPE
!     / /
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, KW(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW
!     /VBLES/
      DOUBLE PRECISION :: LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,WL,
     &     WCV, HP, WT, PHIM, GMR, SEP, M3, PX(NPX), SX(NPX,NM+1),QA(NM)
      COMMON /VBLES / LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,
     & WL, WCV, HP, WT, PHIM, GMR, SEP, M3, PX, SX, QA
!     /TVBLES/
      DOUBLE PRECISION :: DT, ZQ(1), HSPN(2), RLF(2), ZET(2),
     &     XIT(2), AGE,BM, MC(2), OM, PER, SM, ENC, TC(2), TFR, T0, M0,
     &     MTA, OM0, OMTA,A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS,
     &     SE, TN(2), WMH,WMHE, MH, MHE, MCO, VMG, BE(2), LH, LHE, LC,
     &     LNU, LTH, MCB(8),MSB(6), RCB(8), TCT(8), QR(81), PPR(81)
      INTEGER :: JHOLD, JM2, JM1 
      COMMON /TVBLES/ DT, ZQ, HSPN, RLF, ZET, XIT, AGE, 
     & BM, MC, OM, PER, SM, ENC, TC, TFR, T0, M0, MTA, OM0, OMTA, 
     & A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS, SE, TN, WMH, 
     & WMHE, MH, MHE, MCO, VMG, BE, LH, LHE, LC, LNU, LTH, MCB, 
     & MSB, RCB, TCT, QR, PPR, JHOLD, JM2, JM1


      ! Main sequence stars: hydrogen in the core
      IF ( SX(10, 2) > 1.0D-5 ) THEN
         FIND_STELLAR_TYPE = 0
         IF (QCNV < 0.1D0) THEN
            FIND_STELLAR_TYPE = 1
         END IF
         RETURN
      END IF

      ! Hertzsprung gap/end of core hydrogen burning
      IF (SX(10,2) < 1.0D-5 .AND. SX(11,2) > 1.0D-5
     &       .AND. MH < 0.22 .AND. SX(20,2) < 1.0) THEN
         FIND_STELLAR_TYPE = 2
         RETURN
      END IF

      ! Check for evolved stars without nuclear burning
      IF (LHE < 1.0 .AND. LH < 1.0 .AND. LC < 1.0) THEN
         IF (SX(11,2) > 1.0D-5) THEN      ! He white dwarf
            FIND_STELLAR_TYPE = 10
            RETURN
         END IF
         IF (SX(12,2) > 1.0D-5) THEN      ! C/O white dwarf
            FIND_STELLAR_TYPE = 11
            RETURN
         END IF
         ! If none of the above, call it an O/Ne white dwarf
         ! It's not likely to be a neutron star
         FIND_STELLAR_TYPE = 12
         RETURN
      END IF

      ! Red giant branch: helium core, but no helium burning
      IF ( LHE < 3.0 .AND. SX(11, 2) > 1.0D-5 ) THEN
         FIND_STELLAR_TYPE = 3
         RETURN
      END IF

      ! AGB: inert C/O core
      IF ( SX(11, 2) < 1.0D-5 .AND. LC < 1.0 ) THEN
         ! Early AGB: shell sources are well seperated
         IF ( MHE < 0.75*MH ) THEN     ! 0.75? Should be closer?
            FIND_STELLAR_TYPE = 5
         ELSE                          ! T-P AGB
            FIND_STELLAR_TYPE = 6
         END IF
         RETURN
      END IF

      ! Some kind of helium star (hopefully, otherwise it's weird)
      IF ( SX(11, 2) > 1.0D-5 .AND. SX(10, KH+1) < 0.1D0 ) THEN
         FIND_STELLAR_TYPE = 7
         RETURN
      END IF

      ! Some kind of helium star (hopefully, otherwise it's weird)
      IF ( SX(11, 2) < 1.0D-5 .AND. SX(10, KH+1) < 0.1D0 ) THEN
         IF (LC < 1.0) THEN
            FIND_STELLAR_TYPE = 8
         ELSE
            FIND_STELLAR_TYPE = 9
         END IF
         RETURN
      END IF

      ! We should probably ever get here... what ever this is, call it a
      ! Hertzsprung gap star? Presumably it's out of thermal equilibrium...
      FIND_STELLAR_TYPE = 2
      END FUNCTION



      SUBROUTINE COMPUTE_OUTPUT_QUANTITIES ( JSTAR )
      USE MESH
      USE MESH_ENC
      USE FUDGE_CONTROL
      USE CONSTANTS
      USE SETTINGS
      USE PLOTVARIABLES
      USE PRINTB_GLOBAL_VARIABLES
      IMPLICIT NONE
!     Subroutine Arguments
      INTEGER, INTENT (IN) :: JSTAR
!     / /
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, KW(260)
!     /STORE /
      DOUBLE PRECISION :: HPR(NVAR,NM), HT(4,NM), MS(60025)
!     /QUERY /
      DOUBLE PRECISION :: ML, QL, XL, UC(21)
      INTEGER ::  JMOD, JB, JNN, JTER, JOC, JKH
!     /TVBLES/
      DOUBLE PRECISION :: DT, ZQ(1), HSPN(2), RLF(2), ZET(2),
     $     XIT(2), AGE,BM, MC(2), OM, PER, SM, ENC, TC(2), TFR, T0, M0,
     $     MTA, OM0, OMTA,A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS,
     $     SE, TN(2), WMH,WMHE, MH, MHE, MCO, VMG, BE(2),  LH, LHE, LC,
     $     LNU, LTH, MCB(8),MSB(6), RCB(8), TCT(8), QR(81), PPR(81)
      INTEGER :: JHOLD, JM2, JM1 
!     /VBLS/
      DOUBLE PRECISION :: LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,WL,
     $     WCV, HP, WT, PHIM, GMR, SEP, M3, PX(NPX), SX(NPX,NM+1),QA(NM)
!     /INF /
      DOUBLE PRECISION :: Q(NVAR), QD(NVAR), 
     $     FILLUP_COMMONBLOCK_INF(NFUNC+NVAR*NFUNC)
!     /ABUND /
      DOUBLE PRECISION :: WA1(9), YA(9), WNE0, WA(2), AVM, WNE
!     /STAT2 /
      DOUBLE PRECISION :: PL, RL, U, P, RHO, FK, T, SF, ST, ZT, GRADA,
     $     SCP,RF, RT, XHI, S, PR, PG, PF, PT, EN, RPP, R33, R34, RBE,
     $     RBP,RPC, RPNA, RPO, R3A, RAC, RAN, RAO, RANE, RCCA, RCO ,ROO,
     $     RGNE, RGMG, RCCG, RPNG, EX, ENX, WMU, DELTA, PHI, EXT, FKT,
     &     FKR, PRANDTL
!     /ATDATA /
      DOUBLE PRECISION :: CH2(4), CHI(26,9), COM(27), CAN(9), CBN(9)
      INTEGER ::  KZN(9)
!     /NCDATA /
      DOUBLE PRECISION :: QPP, Q33, Q34, QBE, QBP, QPC, QPNA, QPO, Q3A,
     $     QAC,QAN, QAO, QANE, QCCA, QCO, QOO, QGNE, QGMG, QCCG, QPNG,
     $     CNUC(109)
!     /LSMOOTH /
      DOUBLE PRECISION :: LK_PREV(NM), LQ_PREV(NM), ENUC_PREV(NM)

      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW
      COMMON /STORE / HPR, HT, MS
      COMMON /QUERY / ML, QL, XL, UC, JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /TVBLES/ DT, ZQ, HSPN, RLF, ZET, XIT, AGE, 
     : BM, MC, OM, PER, SM, ENC, TC, TFR, T0, M0, MTA, OM0, OMTA, 
     : A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS, SE, TN, WMH, 
     : WMHE, MH, MHE, MCO, VMG, BE, LH, LHE, LC, LNU, LTH, MCB, 
     : MSB, RCB, TCT, QR, PPR, JHOLD, JM2, JM1
      COMMON /VBLES / LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,
     : WL, WCV, HP, WT, PHIM, GMR, SEP, M3, PX, SX, QA
! DG is Schwarzschild value, EG is ov value, of gradr-grada
      COMMON /INF   / Q, QD, FILLUP_COMMONBLOCK_INF
      COMMON /ABUND / WA1, YA, WNE0, WA, AVM, WNE
      COMMON /STAT2 / PL, RL, U, P, RHO, FK, T, SF, ST, ZT, GRADA, SCP, 
     : RF, RT, XHI, S, PR, PG, PF, PT, EN, RPP, R33, R34, RBE, RBP, 
     : RPC, RPNA, RPO, R3A, RAC, RAN, RAO, RANE, RCCA, RCO,
     : ROO, RGNE, RGMG, RCCG, RPNG, EX, ENX, WMU, DELTA, PHI, EXT, FKT,
     & FKR, PRANDTL
      COMMON /ATDATA/ CH2, CHI, COM, CAN, CBN, KZN
      COMMON /NCDATA/ QPP, Q33, Q34, QBE, QBP, QPC, QPNA, QPO, Q3A, QAC,
     : QAN, QAO, QANE, QCCA, QCO, QOO, QGNE, QGMG, QCCG, QPNG, CNUC
      COMMON /LSMOOTH/ LK_PREV, LQ_PREV, ENUC_PREV

      INTEGER :: ICNV, ICNV1

      INTEGER :: IC, IS, IO, JE, IG, I, IKK, IK, J
     $     , IJ, IP, IMAX,  JJ,  KK

      DOUBLE PRECISION, EXTERNAL :: PSTV !!!is this a good way?? SdM

      DOUBLE PRECISION :: MCNV, MCNV1
      DOUBLE PRECISION :: FTOUT(NM), FPOUT(NM), XR(NM), XM(NM)
      DOUBLE PRECISION :: EXLIM, TMAX, HPC, WF, PDG,
     $     PEG, RICH, DOMEGA2
      DOUBLE PRECISION :: CGS_SGF, DMULEFT, DMURIGHT, DMU, DWLEFT,
     $     DWRIGHT, DW, Dcgs_CON, Dcgs_THL, Dcgs_DSI, Dcgs_SSI, Dcgs_SHI
     $     , Dcgs_GSF, Dcgs_ESC, TON, RIS, ES_HP, ES_VMU, ES_VES
      DOUBLE PRECISION :: SGF, DRMAX, DDR, DTY

      IC = 1
      IS = 1
      IO = 1
      IM_TMAX = 1
! Set lower limit for nuclear burning to 10*L/M, for printing purposes only (Onno)
      EXLIM = 10.0*H(8 + 24*(JSTAR-1),1)/H(4 + 24*(JSTAR-1),1)
      JE = 1
      MEX(1:12) = 0.0D0
      TMAX = 0.0D0
      PX(1:45) = 0.0D0
C      ZQ(38:80) = 0.0D0    ! replaced by explicit statements below. -SdM
      WMH= 0.0D0                !ZQ(38)
      WMHE= 0.0D0
      MH= 0.0D0 
      MHE= 0.0D0 
      MCO= 0.0D0 
      VMG= 0.0D0
      BE(1:2) = 0.0D0
      LH = 0.0D0
      LHE = 0.0D0 
      LC = 0.0D0
      LNU = 0.0D0
      LTH = 0.0D0               !ZQ(50)
      MCB(1:8) = 0.0D0            !ZQ(51:58)
      MSB(1:6) = 0.0D0            !ZQ(59:64)
      RCB(1:8) = 0.0D0            !ZQ(65:72)
      TCT(1:8) = 0.0D0            !ZQ(73:80)

      DL = 0.0D0
      MCNV = 0.D0
      MCNV1 = 0.D0
      ICNV = 0
      ICNV1 = 0

      ! Rotational periods for the centre and surface
      PSURF = 0.0D0
      PCNTR = 0.0D0
      DTY = DT/CSY
      DO IKK = 1, KH
         IK = KH + 1 - IKK
         DO I = 1, 16
            QD(I) = DH(I + 24*JSTAR - 24, IK)
            Q(I) = H(I + 24*JSTAR - 24, IK) + QD(I)
         END DO
         QD(17:24) = DH(17:24, IK)
         Q(17:24) = H(17:24, IK) + QD(17:24)
         CALL FUNCS1 ( IK, -1 )
         if (ik==1 .and. jstar == 1) then
            !print *, rlf(jstar), mtr
         end if
         IF ( IKK == 1 ) THEN
            HPC = DSQRT(P/(CG*RHO*RHO))
            MC(JSTAR) = 3.5D-33*RHO*HPC**3
            !PCNTR = exp(Q(13))
            !PCNTR = Q(13)
            PCNTR = 2.0*CPI/(DABS(Q(13)) * CSDAY)
         END IF
         IF ( IKK == KH ) THEN
            !PSURF = exp(Q(13))
            !PSURF = Q(13)
            PSURF = 2.0*CPI/(DABS(Q(13)) * CSDAY)
            !print *, JMOD, Q(42), Q(4)
            !print *, PCNTR, PSURF
         END IF
! Evaluate the functions to be printed
         WF = DSQRT(1.0D0 + DEXP(Q(1)))
         PX(1) = DMIN1(99.0D0, 2.0D0*(WF - DLOG(WF + 1.0D0)) + Q(1))
         PX(2) = P
         PX(3) = RHO
         PX(4) = T
         TMAX = MAX(TMAX, T)
         IF ( T == TMAX ) IM_TMAX = IKK + 1
         PX(5) = FK
         PX(6) = GRADA
         PX(7) = GRAD
         PX(8) = DG    
         PDG = PX(8)
         IF ( IKK == KH .AND. PDG > 0.0D0 ) PDG = 0.0D0
         PX(39) = EG    
         PEG = PX(39)
         IF ( IKK == KH .AND. PEG > 0.0D0 ) PEG = 0.0D0
         PX(9) = Q(4)/CMSN
! Correct abundances for actual atomic masses
         DO I = 1, 7
            PX(I + 9) = CAN(I)*YA(I)/AVM
            IF (PX(I + 9)<1.0D-99) PX(I + 9) = 0.0D0
         END DO
         PX(17) = R/CRSN
         PX(18) = Q(8)/CLSN
         PX(19) = ETH
         PX(20) = EX
         PX(21) = -EN 
         PX(27) = U
         PX(28) = S
         PX(29) = LOLEDD
         PX(30) = WL
         PX(31) = CR*RHO*T/PG       ! Mean molecular weight
         PX(32) = WT
         PX(33) = WNE
         PX(34) = QD(4)/(CMSN*DTY)  ! dM/dt
!         PX(34) = WNE0
         PX(35) = WCV
         PX(36) = Q(12)
         PX(37) = Q(14)
         PX(38) = Q(19)
         DO I = 40, 45
            PX(I) = SX(I, IKK + 1)
         END DO
! Reaction rates
         DO I = 50, 54
            PX(I) = SX(I, IKK + 1)
         END DO
! Thermohaline mixing
         PX(23) = SX(23, IKK + 1)
         PX(31) = SX(31, IKK + 1)
         PX(30) = SX(30, IKK + 1)
         PX(7) = SX(7, IKK + 1)
! Differential rotation
         PX(59) = SX(59, IKK + 1)
         PX(60) = SX(60, IKK + 1)
         
! LK and LQ
         LK_PREV(IK) = SX(57, IKK + 1)
         LQ_PREV(IK) = SX(58, IKK + 1)
         ENUC_PREV(IK) = EX
         
         DL = DL + 0.5D0*(PX(40) + SX(40, IKK))/CLSN
! locate convective/radiative radius boundaries (DG = 0)
         IF ( IKK > 1 ) THEN

         IF ( IC < 8 ) THEN
            IF ( SX(8, IKK)*PDG <= 0.0D0 ) THEN
               IC = IC + 1  !Start a new convective layer
               RCB(IC) = (SX(17, IKK)*PDG-PX(17)*SX(8, IKK))/(PDG-SX(8, IKK))
               IF ( PDG > 0.0D0 ) 
     :            TCT(IC + 1) = (PX(17) - RCB(IC))/WCV
               IF ( SX(8, IKK) > 0.0D0 ) TCT(IC)=
     :            TCT(IC) + (RCB(IC) - SX(17, IKK))/SX(35, IKK)
            ELSE
! Integrate dr/(conv. vel.), for convective turnover time
               IF ( WCV > 0.0D0 ) TCT(IC + 1) = TCT(IC + 1) +
     :             (PX(17) - SX(17,IKK)) * (WCV + SX(35,IKK))/
     :                ( WCV*WCV + SX(35,IKK)*(SX(35,IKK) + WCV) )
            END IF
         END IF
! locate convective/semiconvective mass boundaries (EG = CGR = 0.01)

         IF ( (PX(39)-CGR)*(SX(39,IKK)-CGR) <= 0.0D0 .AND. IS <= 6 ) THEN
            MSB(IS) = (SX(9,IKK)*(PX(39) - CGR) - PX(9)*(SX(39,IKK) - CGR))
     :              /ABS(PX(39) - SX(39,IKK))
            IS = IS + 1
         END IF

! MvdS  Switch the counting of MCNV on and off
         IF(PEG > 0.D0.AND.PX(4) < 1.D7) MCNV1 = MCNV1 + PX(9) - SX(9,IKK)
         IF(PEG > 0.D0.AND.PX(4) < 1.D7 .AND. SX(39,IKK)*PEG < 0.0D0) THEN
            MCNV = MCNV + SX(9,IKK) -
     &             (SX(9,IKK)*PEG - PX(9)*SX(39,IKK))/(PEG - SX(39,IKK))
            ICNV = 1
            ICNV1 = 1
         ENDIF
         IF(PEG < 0.D0.AND.PX(4) < 1.D7.AND.
     &                                     SX(39,IKK)*PEG < 0.0D0) THEN
            MCNV = MCNV + (SX(9,IKK)*PEG - PX(9)*SX(39,IKK))/
     &                                  (PEG - SX(39,IKK)) - SX(9,IKK)
            IF(ICNV == 0.AND.PX(4) < 1.D5) MCNV = MCNV + SX(9,IKK)
            ICNV = 0
         ENDIF
         IF (ICNV == 1) MCNV = MCNV + PX(9) - SX(9,IKK)
         IF (IKK == KH.AND.ICNV1 == 0.AND.MCNV == 0.D0.AND.MCNV1 > 0.D0)
     &      MCNV = MCNV1 + SX(9,2)
! locate overshoot mass boundaries (EG = 0)
         IF ( IO < 8 .AND. SX(39, IKK)*PEG <= 0.0D0 ) THEN
            MCB(IO) = (SX(9,IKK)*PEG-PX(9)*SX(39,IKK))/ABS(PEG-SX(39,IKK))
            IO = IO + 1
         END IF
! locate burning shell boundaries (XH = CXB, XHE = CXB, XC = 0.05)
         IF ( PX(10) > CXB .AND. SX(10,IKK) < CXB ) MH = (PX(9)*
     :    (CXB-SX(10,IKK))+SX(9,IKK)*(PX(10)-CXB))/(PX(10)-SX(10,IKK))
         IF ( PX(11) > CXB .AND. SX(11,IKK) < CXB ) MHE = (PX(9)*
     :    (CXB-SX(11,IKK))+SX(9,IKK)*(PX(11)-CXB))/(PX(11)-SX(11,IKK))  
         IF ( PX(12) > 0.05D0 .AND. SX(12,IKK) < 0.05D0 ) MCO = 
     :    (PX(9)*(0.05D0 - SX(12, IKK)) + SX(9, IKK)*(PX(12) - 0.05D0))
     :    /(PX(12) - SX(12, IKK))
! Locate boundaries of nuclear burning regions, EX > EXLIM (Onno)
         IF((PX(20)-EXLIM)*(SX(20,IKK)-EXLIM) <= 0.0D0.AND.JE < 12.AND.
     &        IKK >= 2) THEN
            MEX(JE)= (SX(9,IKK)*(PX(20)-EXLIM)-PX(9)*(SX(20,IKK)-EXLIM))
     :           /(PX(20)-SX(20,IKK))
            JE = JE  + 1
         END IF
! some homology invariants
         WF = 1.0D0/DLOG10(PX(2)/SX(2,IKK))
         PX(24) = WF*DLOG10(PX(3)/SX(3,IKK))
         PX(25) = -WF*DLOG10(PX(17)/(DABS(SX(17,IKK))+1.0D-10))
         PX(26) = -WF*DLOG10(PX(9)/(DABS(SX(9,IKK))+1.0D-10))
! Some integrated quantities, to be printed in the short summary
         END IF   ! IKK>1
         PX(22) = PX(9) - SX(9,IKK)
         IF ( IKK == 1 ) PX(22) = PX(9)
         WMH = WMH + PX(22)*PX(10)
         WMHE = WMHE + PX(22)*PX(11)
         IF ( .NOT. (IKK == 1 .OR. MH == 0.0D0))  THEN
            IF ( MH > SX(9, IKK) ) THEN
               BE(JSTAR) = (PX(9) - MH)*EGR
            ELSE
               BE(JSTAR) = BE(JSTAR) + PX(22)*EGR
            END IF
         END IF
         WF = PX(22)*CMSN/CLSN 
         LH = LH + WF*CME*((QPP + 0.5*Q33)*RPP + QPC*RPC + QPO*RPO
     :             + QPNA*RPNA + QPNG*RPNG)
         LHE = LHE + WF*CME*(Q3A*R3A + QAC*RAC + QAN*RAN
     :             + QAO*RAO + QANE*RANE)
         LC = LC + WF*CME*(QCCA*RCCA + QCCG*RCCG + QCO*RCO +QOO*ROO
     :             + QGNE*RGNE + QGMG*RGMG)
         LNU = LNU - WF*EN
         LTH = LTH + WF*ETH
! save the last 5 places for TWIN variables
         SX(1:39, IKK + 1) = PX(1:39)
      END DO
! Compute FT and FP correction factors
      !XR(1:KH) = SX(17, 2:KH+1)
      !XM(1:KH) = SX(9, 2:KH+1)
      !CALL POTENTIAL(KH, XM, XR, SX(3, 2:NM+1), SX(59, 2:NM+1), FTOUT, FPOUT)
      !DO IK=1, KH
      !   print *, IK, FPOUT(IK), FTOUT(IK)/FPOUT(IK)
      !ENDDO
! Richardson number, angular momentum diffusion coefficients
      DO IKK = 1, KH-1
         RICH = SX(60, IKK+1)
         DOMEGA2 = (SX(59, IKK+1) - SX(59, IKK+2))**2
         IF (RICH>0.0D0 .AND. DOMEGA2>0.0D0 .and. RICH/DOMEGA2 < 1.0D0) THEN
            ! Unstable
            SX(60, IKK+1) = RICH/(1.0D-32 + DOMEGA2);
            SX(61, IKK+1) = SX(61, IKK+1)/SX(65, IKK+1)
         ELSE
            ! Stable
            IF (RICH<0.0D0 .OR. DOMEGA2<0.0D0) THEN
               SX(60, IKK+1) = 1.0D32
            ELSE
               SX(60, IKK+1) = RICH/(1.0D-32 + DOMEGA2);
            END IF
            SX(61, IKK+1) = 0
         END IF
!      SX(59, IKK) = OMEGA
      END DO


!------------------------------------------------------------------ 
! DIFFUSION COEFFICIENTS [cm^2/s] 
!     The diffussion coefficients are partially calculated in FUNCS1
!     (available here through SX) except for the gradients which are
!     calculated in EQUNS1. In this section the gradients are calculated
!     again for the inner meshpoints and stored in SX to able to print
!     them in the mdl file for plotting purposes. -SdM
!------------------------------------------------------------------ 
!     Memory help: Variables available from calculations in FUNCS1
!      SX(23, IKK) = SGTH
!      SX(59, IKK) = OMEGA
!      SX(60, IKK) = XFS(FX_RICH)
!      SX(61, IKK) = XFS(FX_DDSI)
!      SX(62, IKK) = XFS(FX_DSSI)
!      SX(63, IKK) = XFS(FX_VES)
!      SX(64, IKK) = XFS(FX_VMU)
!      SX(65, IKK) =  XFS(FX_SGF)
!      SX(66, IKK) =  XFS(FX_SSSI)
!------------------------------------------------------------------ 
!     We loop from the center to surface when IKK runs from 2 to KH-2
      DO IKK = 3, KH-1 
! Conversion factor: SGF = (4 pi r^2 rho)^2/(dm/dk) [1e11 g/cm**2]
         cgs_SGF =  SX(65, IKK)*1e22
! Molecular weight gradient (DMU)
         DMUleft  =  ( SX(31, IKK) - SX(31, IKK-1)  )
         DMUright =  ( SX(31,IKK+1 )  - SX(31, IKK) )
         DMU      = 0.5*( DMUleft  +DMUright )
! Omega gradient**2 (DW)
         DWleft  = ( SX(59, IKK-1) - SX(59, IKK)   )
         DWright = ( SX(59,  IKK)  - SX(59, IKK+1) )
         DW      = 0.5*( DWleft  +DWright )
         DW      = DW**2 
!------------------------------------------------------------------ 
! [Dcgs_CON] Convective mixing  
!     Either eggletons 72 scheme or pols&Tout01 set by convection_scheme
!     parameter + possibly artificial extra mixing + weakly mixing of
!     inner and outer boundary points 
!     = SG/SGF,
         Dcgs_CON = SX(30, IKK) / cgs_SGF
!     =====================================
! [Dcgs_THL] Thermohaline mixing 
!     = 0.0 or SGTH * DMU / SGF
         Dcgs_THL = SX(23, IKK)* PSTV(DMU, 0.0D0) / cgs_SGF
!     =====================================
! [Dcgs_DSI] Dynamical Shear instability   
!      =  0.0 or  XFS(FX_DDSI)* (1.0D0 - RICH/DW)**2 /cgs_SGF
         RICH = SX(60, IKK)
         Dcgs_DSI = 0.0D0
         IF( DW > 0.0D0 .AND. RICH > 0.0D0 .AND. RICH/DW < 1.0D0 ) THEN
            TON = (1.0D0 - RICH/DW)**2 ! Turn on factor
            Dcgs_DSI= TON*SX(61, IKK)/cgs_SGF
         ENDIF
!     =====================================
! [Dcgs_SSI] Secular Shear Instability   
!     = 0.0 or XFS(FX_DSSI)*DW
         RIS = SX(66, IKK)
         Dcgs_SSI = 0.0D0
         IF ( DW > 0.0D0 .AND. RIS>0.0D0 .AND. RIS/DW < 1.0D0 ) THEN
            Dcgs_SSI =  SX(62, IKK)*DW / cgs_SGF        
         ENDIF
!     =====================================
! [Dcgs_ESC] Eddington-Sweet circulation,
!     compute net circulation velocity
!         HP = SX(67, IKK)
!         ES_HP  = HP * SX(65, IKK) !!!should this be the non converted ////sgf????
!         ES_VMU = SX(64, IKK) * CFMU*DMU
!         ES_VES = SX(63, IKK) 
!         ES_VES = PSTV( DABS(ES_VES) - DABS(ES_VMU), 0.0D0)
!         Dcgs_ESC = 1d-33*ES_VES * ES_HP /  cgs_SGF  
!     SX(76, ) is the diffusion coefficient from Maeder&Zahn (1998),
!     apart from a factor to account for the molecular weight gradient.
         Dcgs_ESC = SX(76, IKK)/(max(1.0d0, 1.0d0-SX(78, IKK)*DMU) * cgs_SGF)  
!     =====================================
! [Dcgs_GSF] Goldreich-Schubert-Fricke instability
! Not yet effective for chemical transport
         Dcgs_GSF = 0.0D0 ! not yet properly implemented
!     =====================================
! [Dcgs_SHI] Solberg-Hoiland Instability
         Dcgs_SHI = 0.0D0 !not yet implemented
!------------------------------------------------------------------ 
!     store gradients in SX:
         SX(69, IKK) = DMU
         SX(70, IKK) = DW
!    store diffusion coefficients in SX:
         SX(71, IKK) =   Dcgs_CON ! Convection + artificial mixing
         SX(72, IKK) =   Dcgs_THL ! Thermohaline
         SX(73, IKK) =   Dcgs_SHI ! Solberg-Hoiland 
         SX(74, IKK) =   Dcgs_DSI ! Dynamical Shear
         SX(75, IKK) =   Dcgs_SSI ! Secular Shear
         SX(76, IKK) =   Dcgs_ESC ! Eddington-Sweet
         SX(77, IKK) =   Dcgs_GSF ! Goldberg-Schubert-Fricke
      END DO
!------------------------------------------------------------------ 
!     set (next-to) boundary values to zero
      DO IKK = 1,2
        SX(69, IKK) =   0.0D0
        SX(70, IKK) =   0.0D0 
        SX(71, IKK) =   0.0D0 ! Convection + artificial mixing
        SX(72, IKK) =   0.0D0 ! Thermohaline
        SX(73, IKK) =   0.0D0 ! Solberg Hoiland 
        SX(74, IKK) =   0.0D0 ! Dynamical Shear
        SX(75, IKK) =   0.0D0 ! Secular Shear
        SX(76, IKK) =   0.0D0 ! Eddington Sweet
        SX(77, IKK) =   0.0D0 ! Goldberg Schubert Frick
      ENDDO
      DO IKK = KH, KH+1
        SX(69, IKK) =   0.0D0
        SX(70, IKK) =   0.0D0 
        SX(71, IKK) =   0.0D0 ! Convection + artificial mixing
        SX(72, IKK) =   0.0D0 ! Thermohaline
        SX(73, IKK) =   0.0D0 ! Solberg Hoiland 
        SX(74, IKK) =   0.0D0 ! Dynamical Shear
        SX(75, IKK) =   0.0D0 ! Secular Shear
        SX(76, IKK) =   0.0D0 ! Eddington Sweet
        SX(77, IKK) =   0.0D0 ! Goldberg Schubert Frick
      ENDDO

!------------------------------------------ END DIFFUSION COEFFICIENTS

        !     SX(60, KH+1) = 1.0D0;!SX(60, KH+1)/MAX((SX(59, KH-1) - SX(59, KH))**2, 1.0d-64)
! variables needed for MB in FUNCS1
      SM = SX(9, KH+1)
      TEFF = PX(4)
      QCNV = MCNV/SM

! convective envelope turnover time
      DRMAX = 0.0D0
      IMAX = 1
      DO I = 1, IC
         IF ( I == 1 ) DDR = RCB(1)
         IF ( I > 1 ) DDR = RCB(I) - RCB(I - 1)
         IF ( DDR >= DRMAX .AND. TCT(I) /= 0.0D0)  IMAX = I
         IF ( DDR >= DRMAX .AND. TCT(I) /= 0.0D0)  DRMAX = DDR
      END DO
      TET = 0.4D0*TCT(IMAX)*CRSN/5.76D-7+1.0D-2
      RCZ = RCB(IMAX)/PX(17)
      DRCZ = DRMAX/PX(17)

      TKH = 1.0D22*CG*Q(4)*Q(4)/(R*Q(8)*CSY)
      TN(JSTAR) = 4.0D10*SM/Q(8)
      DMT = SX(34, KH+1)  ! dM/dt
!      DMT = QD(4)/(CMSN*DTY)  
      DMTR = XIT(JSTAR)*CSY/CMSN  ! dMtrans/dt 
      DMSW = ZET(JSTAR)*CSY/CMSN  ! dMwind/dt
      OMS = OM/CMSN
      FR = DLOG10(SX(17, KH+1))
      FL = DLOG10(SX(18, KH+1))
      SDC = DLOG10(SX(3,2))
      SDM = DLOG10(SX(3,IM_TMAX))
      SDS = DLOG10(SX(3, KH+1))
      STC = DLOG10(SX(4,2))
      STM = DLOG10(SX(4,IM_TMAX))
      STS = DLOG10(SX(4, KH+1))
      PERC = 0.5D0*DLOG10(PX(17)**3/(8.157D0*SM))
! IF M.I. isn't being solved for, VK2 = gyr. rad. squared is nonsense.
      VK2 = DMIN1(SX(36,KH)/(Q(4)*R*R), 1.0D0)
      RAF = DSQRT(DABS(RA2))/R
      CALL POTENT ( Q(4)/OM, DPHI )
      EPHI = 1.0D22*CG*(Q(4) + OM)*DPHI/SEP
      F1 = Q(14)/EPHI
      DF = 0.0
      IF ( KTW == 2 ) DF = (H(62 - 24*JSTAR, 1) + 
     :                     DH(62 - 24*JSTAR, 1))/EPHI - F1
      HORB = Q(17)            ! This seems to not have been updated before?
      ECC  = Q(18)
      HTOT = HORB + HSPN(1) + HSPN(2)

      AJ = AGE + DTY
      IF ( JNN == 0 .AND. KTW == JSTAR ) AGE = AJ

      END SUBROUTINE



      SUBROUTINE PRINTB ( DTY, JO, JCM, RZ, JSTAR, IFT )
      USE MESH
      USE MESH_ENC
      USE FUDGE_CONTROL
      USE CONSTANTS
      USE SETTINGS
      IMPLICIT NONE
!     Subroutine Arguments
      DOUBLE PRECISION, INTENT(INOUT) :: DTY
      DOUBLE PRECISION, INTENT(OUT) :: RZ
      INTEGER, INTENT (IN) :: JSTAR, IFT
      INTEGER, INTENT (INOUT):: JO, JCM
!     /QUERY /
      DOUBLE PRECISION :: ML, QL, XL, UC(21)
      INTEGER ::  JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /QUERY / ML, QL, XL, UC, JMOD, JB, JNN, JTER, JOC, JKH
!     /TVBLES/
      DOUBLE PRECISION :: DT, ZQ(1), HSPN(2), RLF(2), ZET(2),
     &     XIT(2), AGE,BM, MC(2), OM, PER, SM, ENC, TC(2), TFR, T0, M0,
     &     MTA, OM0, OMTA,A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS,
     &     SE, TN(2), WMH,WMHE, MH, MHE, MCO, VMG, BE(2), LH, LHE, LC,
     &     LNU, LTH, MCB(8),MSB(6), RCB(8), TCT(8), QR(81), PPR(81)
      INTEGER :: JHOLD, JM2, JM1 
      COMMON /TVBLES/ DT, ZQ, HSPN, RLF, ZET, XIT, AGE, 
     & BM, MC, OM, PER, SM, ENC, TC, TFR, T0, M0, MTA, OM0, OMTA, 
     & A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS, SE, TN, WMH, 
     & WMHE, MH, MHE, MCO, VMG, BE, LH, LHE, LC, LNU, LTH, MCB, 
     & MSB, RCB, TCT, QR, PPR, JHOLD, JM2, JM1
      INTEGER :: IG, JMAD
      INTEGER, SAVE :: LMOD(2)

! ----------------------------------------------------------
! Compute quantities
! ----------------------------------------------------------
      CALL COMPUTE_OUTPUT_QUANTITIES ( JSTAR )
      RZ = RLF(JSTAR)

! -------------------------------------------------------------------
! Update quantities that have to be calculated explicitly
! -------------------------------------------------------------------
      CALL UPDATE_EXPLICIT_QUANTITIES( JO, JSTAR )

! -------------------------------------------------------------------
! Update the control parameters for the next timestep
! -------------------------------------------------------------------
      CALL UPDATE_TIMESTEP_PARAMETERS( JSTAR )

! -------------------------------------------------------------------
! Adjust the weights of the mesh spacing function
! -------------------------------------------------------------------
      CALL ADJUST_MESH_SPACING

! -------------------------------------------------------------------
! Check termination conditions
! -------------------------------------------------------------------
      CALL CHECK_STOP_CONDITIONS ( JSTAR, JO, JCM )
! Don't look at age during ZAHB construction
      IF ( IFT == 24 .AND. JO == 5) JO = 0

! -------------------------------------------------------------------
! Output
! -------------------------------------------------------------------
      IG = MAX0(JSTAR, JB)
      IF ( IFT == 24 ) IG = 25

      JMAD = JMOD + 1                  ! Model number for output, count from 1
      IF ( JNN == 0 ) JMAD = -JMOD     ! Backup, print negative model number

! Decide whether to print the interior or not, 1,2: file.out1,2
      IF ( MOD(JNN, KT1) == 1 .OR. JMOD == 0 ) CALL WRITE_INTERNAL_DETAILS(IG)

! Print the short summary, for every KT4'th model, 1,2: file.out1,2
      IF ( KT4>0 .AND. MOD(JMOD, KT4) == 0 ) CALL WRITE_SUMMARY ( JSTAR, JMAD, IG )

! Funny composition profile - flag this as a problem
      IF ( MH < MHE .OR. MHE < MCO ) WRITE (IG,*) 'ccc', MH, MHE, MCO
      CALL FLUSH ( IG )

! Write some quantities of each model for plotting purposes, 31,32: file.plt1,2
! Print a line for every KT4'th model
! Don't output models during ZAHB construction
      IF (IG <= 2 .AND. MOD(JMAD, KT4) == 0 .AND. JMAD > 0) THEN
! Overwrite last line if model didn't converge  MvdS
         IF (KT4>0) THEN
            DO WHILE (JMAD <= LMOD(IG)) 
               BACKSPACE(30+IG)
               LMOD(IG) = LMOD(IG) - KT4
            END DO
         END IF
         LMOD(IG) = JMAD
         CALL WRITE_PLT_FILE(JSTAR, JMAD, 30+IG, KT4)
      END IF

! Write detailed model data for plotting purposes, 33,34: file.mdl1,2
! Write one model every KT1'th model
! Don't output models during ZAHB construction
! FIXME: in case of convergence failure, go back and delete previous
! model, as done for the PLT file above
      IF (IG <= 2 .AND. MOD(JMAD, KT1) == 0) CALL WRITE_MDL_FILE(JSTAR, JMAD, 32+IG)
      RETURN
      END



      SUBROUTINE WRITE_INTERNAL_DETAILS ( IG )
      USE SETTINGS
      USE CONSTANTS
      USE MESH_ENC
      USE PLOTVARIABLES
      USE PRINTB_GLOBAL_VARIABLES
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IG
!     / /
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, KW(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW
!     /QUERY /
      DOUBLE PRECISION :: ML, QL, XL, UC(21)
      INTEGER ::  JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /QUERY / ML, QL, XL, UC, JMOD, JB, JNN, JTER, JOC, JKH
!     /VBLES/
      DOUBLE PRECISION :: LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,WL,
     &     WCV, HP, WT, PHIM, GMR, SEP, M3, PX(NPX), SX(NPX,NM+1),QA(NM)
      COMMON /VBLES / LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,
     & WL, WCV, HP, WT, PHIM, GMR, SEP, M3, PX, SX, QA
!     /TVBLES/
      DOUBLE PRECISION :: DT, ZQ(1), HSPN(2), RLF(2), ZET(2),
     &     XIT(2), AGE,BM, MC(2), OM, PER, SM, ENC, TC(2), TFR, T0, M0,
     &     MTA, OM0, OMTA,A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS,
     &     SE, TN(2), WMH,WMHE, MH, MHE, MCO, VMG, BE(2), LH, LHE, LC,
     &     LNU, LTH, MCB(8),MSB(6), RCB(8), TCT(8), QR(81), PPR(81)
      INTEGER :: JHOLD, JM2, JM1 
      COMMON /TVBLES/ DT, ZQ, HSPN, RLF, ZET, XIT, AGE, 
     & BM, MC, OM, PER, SM, ENC, TC, TFR, T0, M0, MTA, OM0, OMTA, 
     & A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS, SE, TN, WMH, 
     & WMHE, MH, MHE, MCO, VMG, BE, LH, LHE, LC, LNU, LTH, MCB, 
     & MSB, RCB, TCT, QR, PPR, JHOLD, JM2, JM1
      INTEGER :: I, J, IJ, IKK, IP, IK
!     realigned, to have 10 on a row, SdM
!     extended to include diffussion coefficients 71-77, SdM
      CHARACTER*5 CHAR(0:NPX)
      DATA CHAR/'     ',' psi ','  P  ',' rho ','  T  ','kappa','grada',' grad',
     :'gr-ga','  m  ','  H1 ',' He4 ',' C12 ',' N14 ',' O16 ',' Ne20',
     :' Mg24','  r  ','  L  ',' Eth ',' Enuc',' Eneu','  dm ','SGTH ',
     :'n/n+1','lr/lp','lm/lp',' U   ',' S   ','L/Edd',' w.l ','  mu ',
     :'  wt ',' Ne  ','dm/dt','  wcv',' M.I.','phi  ','  xi ',' DGOS',
     :'dDLdk','Denth',' xik ','v**2 ',' F2  ',' F1  ','  46 ','  47 ',
     :'  48 ','  49 ',' RPP ',' RPC ',' RPNG',' RPN ',' RPO ',' RAN ',
     :'dS/dP','  LK ','  LQ ','Omega','RichN','DDSI ','  62 ','  63 ',
     :'  64 ','  65 ','  66 ','  67 ','  68 ','  69 ','  70 ', ' 71 ',
     :' 72 ', ' 73 ','  74 ','  75 ','  76 ','  77 ','  78 ','  79 ',
     :'  80 ' /


99001 FORMAT (/, '  K', 15(4X, A5, 1X),/)
99002 FORMAT (I4, 1P, 15D10.3)
99005 FORMAT (/, ' K ', 12(5X, A5, 2X),/)
99006 FORMAT (I4, 1P, 12D12.5)
      IF(KT3 < 1) RETURN
      
      WRITE (IG, 99001) (CHAR(KSX(I)), I = 1, 15)
!     Print the interior details on first `page', if required
      DO IKK = 1, KH
         IK = KH + 1 - IKK
         IF ( MOD(IK-1,KT2) == 0 ) WRITE(IG,99002) IK,
     :        (SX(KSX(J), IKK + 1), J=1,15)
         IF ( MOD(IK,10*KT2) == 0 ) WRITE(IG, 99001) 
     :        (CHAR(KSX(J)), J = 1, 15)
      END DO
!     Write further `pages' for each detailed model, if required
      DO I = 2, KT3
         IJ = 15*I - 15
         WRITE(IG, 99005) (CHAR(KSX(IP + IJ)), IP=1,12)
         DO IKK = 1, KH
            IK = KH + 1 - IKK
!     Have to choose the SX's that are wanted for printing
            IF ( MOD(IK-1, KT2) == 0 ) WRITE(IG, 99006) IK,
     :           (SX(KSX(IP + IJ), IKK + 1), IP = 1, 12)
            IF ( MOD(IK, 10*KT2) == 0 ) WRITE(IG, 99005) 
     :           (CHAR(KSX(IP + IJ)), IP = 1, 12)
         END DO
      END DO
      END SUBROUTINE



      SUBROUTINE WRITE_PLT_FILE ( JSTAR, JMAD, IG, KT4 )
      USE CONSTANTS
      USE MESH_ENC
      USE PLOTVARIABLES
      USE PRINTB_GLOBAL_VARIABLES
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JSTAR, JMAD, IG, KT4
!     / /
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, KW(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW
!     /QUERY /
      DOUBLE PRECISION :: ML, QL, XL, UC(21)
      INTEGER ::  JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /QUERY / ML, QL, XL, UC, JMOD, JB, JNN, JTER, JOC, JKH
!     /VBLES/
      DOUBLE PRECISION :: LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,WL,
     &     WCV, HP, WT, PHIM, GMR, SEP, M3, PX(NPX), SX(NPX,NM+1),QA(NM)
      COMMON /VBLES / LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,
     & WL, WCV, HP, WT, PHIM, GMR, SEP, M3, PX, SX, QA
!     /TVBLES/
      DOUBLE PRECISION :: DT, ZQ(1), HSPN(2), RLF(2), ZET(2),
     &     XIT(2), AGE,BM, MC(2), OM, PER, SM, ENC, TC(2), TFR, T0, M0,
     &     MTA, OM0, OMTA,A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS,
     &     SE, TN(2), WMH,WMHE, MH, MHE, MCO, VMG, BE(2), LH, LHE, LC,
     &     LNU, LTH, MCB(8),MSB(6), RCB(8), TCT(8), QR(81), PPR(81)
      INTEGER :: JHOLD, JM2, JM1 
      COMMON /TVBLES/ DT, ZQ, HSPN, RLF, ZET, XIT, AGE, 
     & BM, MC, OM, PER, SM, ENC, TC, TFR, T0, M0, MTA, OM0, OMTA, 
     & A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS, SE, TN, WMH, 
     & WMHE, MH, MHE, MCO, VMG, BE, LH, LHE, LC, LNU, LTH, MCB, 
     & MSB, RCB, TCT, QR, PPR, JHOLD, JM2, JM1
      LOGICAL :: FIRST_TIME(2) = (/.TRUE., .TRUE./)
      INTEGER :: I, J
      DOUBLE PRECISION :: DTY
      DTY = DT/CSY

! Write number of cols to the first line of the star.plt[12] files
      IF(FIRST_TIME(JSTAR)) THEN
         REWIND(IG)
         WRITE(IG,'(I4)') 83
         FIRST_TIME(JSTAR) = .FALSE.
      ENDIF
      WRITE (IG,12100) JMAD,AJ,DTY,SM,MH,MHE,MCO,
     &      LOG10(PX(17)),LOG10(PX(18)),LOG10(PX(4)),LOG10(SX(4,2)),STM,
     &      LOG10(SX(3,2)),SDM,BE(JSTAR),LH,LHE,LC,LNU,LTH,
     &      PSURF,VK2,RCZ,DRCZ,TET,RAF,BP,
     &      PER,RLF(JSTAR),F1,DMT,DMSW,DMTR,HORB,DHDT,DHGW,DHML,DHSO(JSTAR),
     &      DHMT(JSTAR),OMS,ECC,
     &      (PX(J),J=10,16),(SX(I,IM_TMAX),I=10,16),(SX(J,2),J=10,16),
     &      (MCB(J),J=1,6),(MSB(J),J=1,6),(MEX(J),J=1,6),QCNV,SX(2,2), PCNTR
      CALL FLUSH(IG)
12100 FORMAT (I6,ES17.9,ES14.6,12ES13.5,7ES12.4,3ES13.5,17ES12.4,
     &  39ES13.5,ES14.6,ES13.5,ES14.6)
      END SUBROUTINE


!     Write structure output (the values of some physical parameters on
!     each meshpoint) to a file usually called star.mdl1 or star.mdl2
!
!     INPUT VARIABLES
!     IG = indicates which starr: 1 for primary, 2 for secondary component
!     JMAD = model number
      SUBROUTINE WRITE_MDL_FILE ( JSTAR, JMAD, IG )
      USE MESH
      USE SETTINGS
      USE PLOTVARIABLES
      USE PRINTB_GLOBAL_VARIABLES
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JSTAR, IG, JMAD
!     / /
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, KW(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW
!     /VBLES/
      DOUBLE PRECISION :: LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,WL,
     &     WCV, HP, WT, PHIM, GMR, SEP, M3, PX(NPX), SX(NPX,NM+1),QA(NM)
      COMMON /VBLES / LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,
     & WL, WCV, HP, WT, PHIM, GMR, SEP, M3, PX, SX, QA
  
      LOGICAL :: FIRST_TIME(2) =(/ .true., .true./)       
      INTEGER :: KK, J
      INTEGER, PARAMETER :: NWRT5 = 42
      INTEGER :: IPX(NWRT5)   
!     IPX: List of Indices of the entries in the array SX containing
!     physcial pars as calculated during latest call cfuncs. This
!     determines the order in which these pars are written to the output
!     file .mdl[12]. 
      DATA IPX/9, 17,  2,  3,  4,    5,  6 , 8, 10, 11,   
     &        12, 13, 14, 15, 16,   18, 19, 20, 21, 28,   
     &        27, 50, 51, 52, 53,   54, 55, 56, 31, 23,
     &        30,  7, 59, 69, 70,   71, 72, 73, 74, 75,
     &        76, 77/! 67, 68, 60, 60, 61,   66, 63, 64, 65/!, 24, 25, 26/ 

!     If this function is called for the first time, write a header
!     (expected by onno's yorick routines)
      IF ( FIRST_TIME(JSTAR)) THEN
         REWIND(IG)
         WRITE (IG,'(2I6,F7.3)') KH, NWRT5, COS      
         FIRST_TIME(JSTAR) = .FALSE.
      ENDIF
! Every time this fuction is called write a block header, and a datablock
      WRITE (IG,12201) JMAD, AJ
      DO KK = 1, KH
         ! Limit range of output exponent so it fits in the output format
         WRITE (IG,12202) 
     &     (dsign(MIN(dabs(SX(IPX(J),KK+1)),1.0D99),
     &             SX(IPX(J),KK+1) ), J=1, NWRT5)  
      END DO
      CALL FLUSH(IG)
12201 FORMAT (I6,1P,E17.9,0P)
12202 FORMAT (1P,E13.6,4E11.4,80E11.3,0P)
      END SUBROUTINE


      SUBROUTINE WRITE_SUMMARY ( JSTAR, JMAD, IG )
      USE CONSTANTS
      USE MESH_ENC
      USE PLOTVARIABLES
      USE PRINTB_GLOBAL_VARIABLES
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JSTAR, JMAD, IG
!     / /
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, KW(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW
!     /QUERY /
      DOUBLE PRECISION :: ML, QL, XL, UC(21)
      INTEGER ::  JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /QUERY / ML, QL, XL, UC, JMOD, JB, JNN, JTER, JOC, JKH
!     /VBLES/
      DOUBLE PRECISION :: LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,WL,
     &     WCV, HP, WT, PHIM, GMR, SEP, M3, PX(NPX), SX(NPX,NM+1),QA(NM)
      COMMON /VBLES / LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,
     & WL, WCV, HP, WT, PHIM, GMR, SEP, M3, PX, SX, QA
!     /TVBLES/
      DOUBLE PRECISION :: DT, ZQ(1), HSPN(2), RLF(2), ZET(2),
     &     XIT(2), AGE,BM, MC(2), OM, PER, SM, ENC, TC(2), TFR, T0, M0,
     &     MTA, OM0, OMTA,A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS,
     &     SE, TN(2), WMH,WMHE, MH, MHE, MCO, VMG, BE(2), LH, LHE, LC,
     &     LNU, LTH, MCB(8),MSB(6), RCB(8), TCT(8), QR(81), PPR(81)
      INTEGER :: JHOLD, JM2, JM1 
      COMMON /TVBLES/ DT, ZQ, HSPN, RLF, ZET, XIT, AGE, 
     & BM, MC, OM, PER, SM, ENC, TC, TFR, T0, M0, MTA, OM0, OMTA, 
     & A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS, SE, TN, WMH, 
     & WMHE, MH, MHE, MCO, VMG, BE, LH, LHE, LC, LNU, LTH, MCB, 
     & MSB, RCB, TCT, QR, PPR, JHOLD, JM2, JM1
11993 FORMAT (
     :I6, ' M/dty/    Po/P*/e/   xi/zet/   tN/tKH/  LH/LHE/LC Pcr/Rcz ', 
     :'  MH/Me/MC/  XH     XHe     XC      XN      XO      XNe     XMg',
     :'      psi   logrho  logT ', /,
     :'  age/cm/Horb     rf1/Fn   mt/rf2/dF tET/DL/BE /Lnu/Lth  ',
     :'DRcz/RA/Bp  MH/Me convective/semiconvective boundaries         ',
     :'     k**2           logR    logL ',  /,
     : 3(1P, D16.9, 5D10.3, 0P, F8.3, 7F8.5, 2F8.4, F9.5, 2X, A4, /), 
     :  (1P, D16.9, 5D10.3, 0P, F8.3, 7F8.3, 2F8.4, F9.5, 2X, A4, /),  
     :   1P, D16.9, 5D10.3, 0P, 8F8.3, 24X, I6, I2, /)

      INTEGER :: JJ, I
      DOUBLE PRECISION :: DTY
      JJ = MAX0(JB, JSTAR)
      DTY = DT/CSY
! KTW=2: suffix JSTAR for things set in PRINTB, 1 for things set in FUNCS1
      WRITE (IG, 11993) 
     &    JMAD, SM, PER, DMTR, TN(JSTAR), LH, PERC, MH,
     &    (SX(I, 2), I = 10, 16), SX(1,2), SDC, STC, 'cntr', 
     &    DTY, PSURF, DMSW, TKH, LHE, RCZ, MHE, (SX(I, IM_TMAX), I = 10, 16),
     &    SX(1, IM_TMAX), SDM, STM,'Tmax', 
     &    AJ, SE, DMT, TET, LC, DRCZ, MCO, (SX(I, KH), I = 10, 16), PX(1), SDS, 
     &    STS, 'srfc', 
     &    OMS, RLF, DL, LNU, RAF, WMH, MCB, FR, FL, '    ', 
     &    HTOT, F1, DF, BE(JSTAR), LTH, BP, WMHE, MSB, VK2, JMAD, JJ
      IF (ADJ_MEA) THEN
         WRITE (IG, '(1P, 3D16.9, I6)'), CURR_MAXDIFFSQR, CURR_DIFFSQR, BEST_DIFFSQR, BEST_MOD
      END IF
      CALL FLUSH ( IG )
      END SUBROUTINE
