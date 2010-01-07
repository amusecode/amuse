! Read first model from JIP and write to JOP
! Read init.dat from IT
      SUBROUTINE STAR12 ( JO, JCM, JOP, JIP, KSV, IT )
      USE MESH
      USE MESH_ENC
      USE FUDGE_CONTROL
      USE SETTINGS
      USE EXTRAPOLATE_DH
      USE EXPLICIT_FUNCTIONS
      USE NUCLEOSYNTHESIS
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, ID(130), IE(130) 
      COMMON /QUERY / ML(3), UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /TVBLES/ DT, ZQ(81), PR(81), PPR(81), IHOLD, JM2(2)

      DOUBLE PRECISION :: CDD_1, CDD_2
C Read an initial model
      JNN = 0
      EXPL_VAR(:, :, :) = 0.0D0
      CALL BEGINN ( JOP, JIP, DTY, KP, KR1, KR2, KSV, KT5, IT, JO, JF )
      IF ( JO == -2 ) GO TO 3
      CALL PRINTB ( DTY, JO, JCM, RZ, 1, IT)
      IF ( JO.EQ.1 ) GO TO 3         
      IF ( KTW.EQ.2 ) THEN
      !average change wanted by *1
         CDD_1 = ZQ(30)
         CALL PRINTB ( DTY, JO, JCM, RZ, 2, IT )
      !average change wanted by *2
         CDD_2 = ZQ(30)
      !set average change to minimum of the two components
         ZQ(30) = min(CDD_1, CDD_2)
      ENDIF
      IF ( KP.EQ.0 ) GO TO 3
      CALL NEXTDT ( DTY, JO, IT )
      IF ( JO.EQ.3 ) GO TO 3
      ITER = KR1
C Begin evolutionary loop of KP time steps
      JNN = 1
      JTER = 0
      DO
         JO = 0
C Solve for structure, mesh, and major composition variables
         JOC = 1

         CALL SMART_SOLVER ( ITER, ID, KT5, JO )

C If no convergence, restart from 2 steps back, DT decreased substantially
    5    IF (JO /= 0) THEN
            CALL BACKUP ( DTY, JO )
            IF (USE_QUADRATIC_PREDICTIONS) THEN
               CALL RESTORE_PREVIOUS_MODEL()
               CALL PREDICT_DH(DTY, DH)
            END IF
            IF ( JO.EQ.2 ) EXIT     ! Abort if timestep below limit
            GO TO 4
         END IF
C If required, solve for minor composition vbles; mesh, structure fixed.
         JOC = 2
         CALL SOLVER ( KR_NUCSYN, IE, KT5, JO )
         IF (JO /= 0) THEN
            JO = 0
            DHNUC(1, :, :) = 0.0d0
            CALL SOLVER ( KR_NUCSYN, IE, KT5, JO )
         END IF
         IF (JO /= 0) THEN
            JO = 0
            DO IK = 1, KH_NUC
               DO IJ = 1, NVAR_NUC
                  IF ( HNUC(JSTAR, IJ, IK) < 1.0d-20) HNUC(JSTAR, IJ, IK) = 0.0d0
               END DO
               IF (H(5 + (JSTAR-1)*24, IK) == 0.0d0) HNUC(JSTAR, 41, IK) = 0.0d0
            END DO
            CALL SOLVER ( KR_NUCSYN, IE, KT5, JO )
         END IF
         IF (JO /= 0) JO = 17
C If no convergence, restart from 2 steps back, DT decreased substantially
         IF (JO /= 0) THEN
            CALL BACKUP ( DTY, JO )
            IF (USE_QUADRATIC_PREDICTIONS) THEN
               CALL RESTORE_PREVIOUS_MODEL()
               CALL PREDICT_DH(DTY, DH)
            END IF
            IF ( JO.EQ.2 ) EXIT     ! Abort if timestep below limit
            GO TO 4
         END IF
!Mvds: If mass < 0, probably NaN, so exit with code 99
! (works for Lahey F95, but probably compiler-dependent!)
         IF(H(4,1).LT.0.) JO = 99
C If model didn't converge, give up
         IF ( JO >= 1 ) EXIT

         PRZ = RZ
         CALL PRINTB ( DTY, JO, JCM, RZ, 1, IT )
         IF ( KTW.EQ.2 ) THEN
            !average change wanted by *1
            CDD_1 = ZQ(30)
            CALL PRINTB ( DTY, JO, JCM, RZ, 2, IT )
            !average change wanted by *2
            CDD_2 = ZQ(30)
            !set average change to minimum of the two components
            ZQ(30) = min(CDD_1, CDD_2)
         ENDIF

         IF ( JO.EQ.14 ) GO TO 5

! Check for manual termination
         READ (11, *) JSTOP
         IF ( JSTOP.NE.0 ) JO = 15
         CLOSE (11)
         WRITE (11, *) 0
         CLOSE (11)

         IF ( JO >= 2 ) EXIT
         CALL UPDATE ( DTY )
         CALL UPDATE2
         IF (USE_QUADRATIC_PREDICTIONS) CALL STORE_CURRENT_MODEL()
!        IF (MOD(JMOD,KSV).EQ.0 .OR. JNN.EQ.1)   !MvdSnew: should save every model and model 1000 iso 1001, etc. JNN=1 saves the first model
         IF ( MOD(JMOD, KSV) == 1 .OR. KSV == 1 )   !MvdS: added 2nd option to save every model, doesn't work with mod
     &       CALL OUTPUT ( KP, 15, JO, JF )
    4    CALL NEXTDT ( DTY, JO, IT )
         IF (USE_QUADRATIC_PREDICTIONS) CALL PREDICT_DH(DTY, DH)
         IF ( JO == 3 ) EXIT
         IF ( JO == 5 ) EXIT        ! Reached age limit
         IF ( JNN >= 4 ) ITER = KR2
         JNN = JNN + 1
C Don't break the run just as RLOF starts!
         KR = DABS(RZ)/(DABS(RZ - PRZ) + 1.0D-6)
         IF ( KP == 1 ) EXIT
         
         IF (.NOT.( (JNN<=KP .OR. KP<0) .OR. IHOLD<3 .OR. JTER>2 .OR. KR<=10 )) EXIT
      ENDDO
 3    CONTINUE
C Output the last converged model, unless it's at the He flash
      CALL OUTPUT ( KP, JOP, JO, JF )
! Save last model to file.mod as well (suggested by Selma)
      IF ( .NOT. (MOD(JMOD, KSV).EQ.1 .OR. KSV.EQ.1) )
     &   CALL OUTPUT ( KP, 15, JO, JF )
      IF ( JO.NE.8.AND.JO.NE.13 ) WRITE (JB, *) 'switch to *', 3 - JB
      IF ( JO.EQ.8 ) WRITE (JB, *) 'He flash'
      IF ( JO.EQ.13 ) WRITE (JB,*) 'start ZAHB'
      CALL FLUSH ( JB )
      RETURN
      END

! Smart solver iterator:
! Tries changing a `fudge factor' from 0 to 1.
! 0 means fully fudged, 1 means not fudged.
! This is called with various factors by the top level solver to attempt
!  relaxation of the solution. It works because FORTRAN passes arguments by
!  reference, so changing the value in this function will change the global.
      SUBROUTINE FUDGED_SOLVER ( ITER, IG, KT5, JO, FACTOR, FH, DFH )
      USE MESH
      USE MESH_ENC
      USE SETTINGS
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ITER, KT5, IG(130);
      INTEGER, INTENT(INOUT) :: JO
      DOUBLE PRECISION, INTENT(INOUT)  :: FH(NVAR,NM), DFH(NVAR,NM)
      DOUBLE PRECISION, INTENT(INOUT) :: FACTOR;
      DOUBLE PRECISION, PARAMETER :: START_FACTOR = 1.0D-10
      INTEGER, PARAMETER :: VERBOSITY = 0;      ! Increase for more output
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, ID(130), IE(130)
      DOUBLE PRECISION :: LH(NVAR,NM), DLH(NVAR,NM)
      INTEGER :: KT6;
      INTEGER :: IFTRY;
      DOUBLE PRECISION :: PREV_FACTOR, OFACTOR;
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, ID 

      IF (JO /= 0) THEN
         KT6 = ITER+30+KT5+1;    ! Supress all status output from SOLVER
         IF (VERBOSITY > 1) KT6 = KT5

! Restore H and DH from before the unsuccesful call to SOLVER
         H(1:NVAR, 1:KH) = FH(1:NVAR, 1:KH)
         DH(1:NVAR, 1:KH) = DFH(1:NVAR, 1:KH)
         
         OFACTOR = FACTOR
         FACTOR = 0.0
         JO = 0
         CALL SOLVER ( ITER, IG, KT6, JO )

         FACTOR = START_FACTOR
         IFTRY = 0
         LH(1:NVAR, 1:KH) = H(1:NVAR, 1:KH)
         DLH(1:NVAR, 1:KH) = DH(1:NVAR, 1:KH)
         PREV_FACTOR = FACTOR
         IF (VERBOSITY>0) WRITE(1, *) 'Starting relaxation'
         DO WHILE (FACTOR<1.0D0 .AND. JO == 0)
            IF (FACTOR < 0.8D0) THEN
               FACTOR = SQRT(FACTOR);
            ELSE
               FACTOR = MIN(1.0D0, 1.1D0*FACTOR)
            ENDIF
            IF (VERBOSITY>0) WRITE(1, *) 'Factor now at', FACTOR

! Restore H from before the last call to SOLVER; keep DH
            H(1:NVAR, 1:KH) = FH(1:NVAR, 1:KH)
            CALL SOLVER ( ITER, IG, KT6, JO )
            IF (JO /= 0) THEN
! Not converged, backtrack and try with a reduced factor
               IFTRY = IFTRY+1
               IF (IFTRY < 5) JO = 0
               DH(1:NVAR, 1:KH) = DLH(1:NVAR, 1:KH)
               FACTOR = PREV_FACTOR*FACTOR
            ELSE
               IFTRY=0
               DLH(1:NVAR, 1:KH) = DH(1:NVAR, 1:KH)
               PREV_FACTOR = FACTOR
            ENDIF
         END DO
         IF (JO == 0) THEN
            IF (VERBOSITY>0) WRITE(1, *) 'Succes! :)'
         ELSE
            IF (VERBOSITY>0) WRITE(1, *) 'Failed :('
         ENDIF
! Set the factor back to 1.0, which means to use normal physics
         FACTOR = OFACTOR
      END IF
      END SUBROUTINE

! SOLVER wrapper that tries a few clever things to recover in case the normal
! SOLVER fails to converge.
      SUBROUTINE SMART_SOLVER ( ITER, IG, KT5, JO )
      USE MESH
      USE MESH_ENC
      USE SETTINGS
      USE FUDGE_CONTROL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ITER, KT5, IG(130);
      INTEGER, INTENT(OUT) :: JO
      DOUBLE PRECISION :: FH(NVAR,NM), DFH(NVAR,NM), LH(NVAR,NM), DLH(NVAR,NM)
      DOUBLE PRECISION :: BCKDH0, LASTEGEN_SMOOTH, PREV_FUDGE, BCKAVMU_SMOOTH
      INTEGER :: IFTRY, KT6
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, ID(130), IE(130)
      INTEGER, PARAMETER :: VERBOSITY = 0;
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, ID 
      
! Store H and DH for possible reuse
      FH(1:NVAR, 1:KH) = H(1:NVAR, 1:KH)
      DFH(1:NVAR, 1:KH) = DH(1:NVAR, 1:KH)
      CALL SOLVER ( ITER, IG, KT5, JO )

! Detect circumstances that would force a smaller timestep; try to work around
!  that first.

! The first is trivial: if we are about to converge, then stopping now is silly
!      IF ( ALLOW_EXTENSION .AND. (JO == 1) .AND. (ERR < ERRPPR)
!     &    .AND. ( (EPS - ERRPPR)/(ERR - ERRPPR) < 0.5D0*ITER ) ) THEN
!         write (1, *) 'Expecting convergence soon - extending iterations'
!         JO = 0
!         H(1:NVAR, 1:KH) = FH(1:NVAR, 1:KH)
!         DH(1:NVAR, 1:KH) = DFH(1:NVAR, 1:KH)
!         CALL SOLVER ( (3*ITER+1)/2, IG, KT5, JO )
!         IF (JO == 0) THEN
!            write(1, *) 'Succes! :)'
!         ELSE
!            write(1, *) 'Failed :('
!         ENDIF
!      END IF

! Try reducing DH0
      IF (JO /= 0) THEN
         KT6 = 3*ITER+KT5;
         IF (VERBOSITY>0) write(1, *) 'Attempting to reduce DH0'
         BCKDH0 = DH0
         DH0 = DH0 / 16.0D0
! Restore H and DH from before the last unsuccesful call to SOLVER
         JO = 0
         H(1:NVAR, 1:KH) = FH(1:NVAR, 1:KH)
         DH(1:NVAR, 1:KH) = DFH(1:NVAR, 1:KH)
         CALL SOLVER ( ITER, IG, KT6, JO )
         IF (JO == 0) THEN
            IF (VERBOSITY>0) write(1, *) 'Succes! :)'
         ELSE
            IF (VERBOSITY>0) write(1, *) 'Failed :('
         ENDIF
         DH0 = BCKDH0
      ENDIF

      IF (JO /= 0) THEN
         KT6 = 3*ITER+KT5;
         IF (VERBOSITY>0) write(1, *) 'Attempting to reduce DH0 more'
         BCKDH0 = DH0
         DH0 = DH0 / 128.0D0
! Restore H and DH from before the last unsuccesful call to SOLVER
         JO = 0
         H(1:NVAR, 1:KH) = FH(1:NVAR, 1:KH)
         DH(1:NVAR, 1:KH) = DFH(1:NVAR, 1:KH)
         CALL SOLVER ( ITER, IG, KT6, JO )
         IF (JO == 0) THEN
            IF (VERBOSITY>0) write(1, *) 'Succes! :)'
         ELSE
            IF (VERBOSITY>0) write(1, *) 'Failed :('
         ENDIF
         DH0 = BCKDH0
      ENDIF

! Try setting DH to 0; this shouldn't usually help, but sometimes reducing
!  the timestep seems to help convergence more than it seems it should,
!  and it is actually clearing DH that helps.
      IF (JO /= 0) THEN
         KT6 = 3*ITER+KT5;
         IF (VERBOSITY>0) write(1, *) 'Trying to clear DH'
! Restore H from before the last unsuccesful call to SOLVER
         JO = 0
         H(1:NVAR, 1:KH) = FH(1:NVAR, 1:KH)
         DH(1:NVAR, 1:KH) = 0.0D0
         CALL SOLVER ( ITER, IG, KT6, JO )
         IF (JO == 0) THEN
            IF (VERBOSITY>0) write(1, *) 'Succes! :)'
         ELSE
            IF (VERBOSITY>0) write(1, *) 'Failed :('
         ENDIF
      ENDIF

      IF (ALLOW_AMADVRELAXATION .AND. JO /= 0) THEN
         IF (JO /= 0 .AND. VERBOSITY>0) WRITE (1, *) 'Attempting to underrelax angular momentum transport'
         CALL FUDGED_SOLVER ( ITER, IG, KT5, JO, AMADV_SMOOTH, FH, DFH);
         IF (JO == 0) THEN
            WRITE (1, *) 'Succes! [Angular momentum transport]'
         ELSE
            WRITE (1, *) 'Failure! [Angular momentum transport]'
         END IF
      END IF

      IF (ALLOW_OVERRELAXATION) THEN
! Try boosting the convective mixing (!) using some weird scheme that seems to
!  work well in some cases.
         IF (JO /= 0 .AND. VERBOSITY>0) WRITE (1, *) 'Attempting to overrelax convective mixing'
         CALL FUDGED_SOLVER ( ITER, IG, KT5, JO, MIXING_BOOST, FH, DFH);
      END IF

      IF (ALLOW_UNDERRELAXATION) THEN
! Try reducing the convective mixing
         IF (JO /= 0 .AND. VERBOSITY>0) WRITE (1, *) 'Attempting to underrelax convective mixing'
         CALL FUDGED_SOLVER ( ITER, IG, KT5, JO, MIXING_FUDGE, FH, DFH);

! Try smoothing out the thermal energy contribution to the luminosity equation
         IF (JO /= 0 .AND. VERBOSITY>0) WRITE (1, *) 'Attempting to underrelax thermal energy release'
         CALL FUDGED_SOLVER ( ITER, IG, KT5, JO, LUMINOSITY_FUDGE, FH, DFH);
      END IF
      
      IF (ALLOW_EGENRELAXATION) THEN
         IF (JO /= 0 .AND. VERBOSITY>0) WRITE (1, *) 'Attempting to underrelax energy generation'
         CALL FUDGED_SOLVER ( ITER, IG, KT5, JO, EGEN_SMOOTH, FH, DFH);

         IF (JO /= 0 .AND. VERBOSITY>0) WRITE (1, *) 'Attempting to underrelax upstream/downstream egen terms'
         CALL FUDGED_SOLVER ( ITER, IG, KT5, JO, LLUMI_SMOOTH, FH, DFH);
      END IF
      
      IF (ALLOW_MDOTRELAXATION) THEN
         IF (JO /= 0 .AND. VERBOSITY>0) WRITE (1, *) 'Attempting to underrelax CMI term'
         CALL FUDGED_SOLVER ( ITER, IG, KT5, JO, MDOT_SMOOTH, FH, DFH);
      END IF

      IF (ALLOW_AVMURELAXATION .AND. .NOT. USE_PREVIOUS_MU) THEN
         IF (JO /= 0 .AND. VERBOSITY>0) WRITE (1, *) 'Attempting to underrelax mean molecular weight'
         CALL FUDGED_SOLVER ( ITER, IG, KT5, JO, AVMU_SMOOTH, FH, DFH);
      END IF

! Various ways to improve the solution that we found, eg, things that were
! computed implicitly from values at the current timestep (for stbility
! reasons).
! Store previous values of some of these settings (to be restored when the 
! solution has been polished; these options need to stay swtiched on/off
! during the solution polishing step)
      BCKAVMU_SMOOTH = AVMU_SMOOTH

! If converged with a value of mu from the previous timestep, try with the
! value of mu from the current timestep now
      IF (ALLOW_AVMURELAXATION .AND. USE_PREVIOUS_MU .AND. JO==0) THEN
         IF (VERBOSITY>0) WRITE (1, *) 'Switching to current molecular weight'
         AVMU_SMOOTH = 1.0D0
         CALL SOLVER ( ITER, IG, KT5/2, JO )
      END IF

! Restore previous values of some settings
      AVMU_SMOOTH = BCKAVMU_SMOOTH

      END SUBROUTINE

      ! Function for reading init.dat (fort.22) in either `classical'
      ! numbers-only format or more modern `with names' format      
      FUNCTION READ_INIT_DAT(IT, KH2, KR1, KR2, KSV, KT5, JCH)
      USE ZAMS_NUCLEOSYNTHESIS_ABUNDANCES
      USE NUCLEOSYNTHESIS
      USE MESH
      USE CONSTANTS
      USE SETTINGS
      USE FUDGE_CONTROL
      USE EXTRAPOLATE_DH
      USE MASSLOSS
      IMPLICIT NONE
      LOGICAL :: READ_INIT_DAT
      INTEGER, INTENT(IN) :: IT
      INTEGER, INTENT(OUT) :: KH2, KR1, KR2, KSV, KT5, JCH
      INTEGER :: IOERROR
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      DOUBLE PRECISION :: CH_OPAC
      INTEGER :: KH, KTW, 
     & KE1, KE2, KE3, KBC, KEV, KFN, KL, JH(3),KP_VAR(40), KP_EQN(40), KP_BC(40), 
     & KE1_2, KE2_2, KE3_2, KBC_2, KEV_2, KFN_2, KL_2, JH_2(3),KP_VAR_2(40), KP_EQN_2(40), KP_BC_2(40)
      DOUBLE PRECISION :: CDC_MS, CDC_HEC, CDC_HES, CDC_DBLSH, CDC5, 
     & CDC_EMS, CDC_HG, CDC_1DUP, CDC_RLOF, CDC_RLOF_REDUCE,
     &      UNUSED1, UNUSED(17)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, 
     & KE1, KE2, KE3, KBC, KEV, KFN, KL, JH, KP_VAR, KP_EQN, KP_BC, 
     & KE1_2, KE2_2, KE3_2, KBC_2, KEV_2, KFN_2, KL_2, JH_2, KP_VAR_2, KP_EQN_2, KP_BC_2
      DOUBLE PRECISION :: X1AC, X4AC, X12AC, X14AC, X16AC, X20AC, X24AC
      DOUBLE PRECISION :: XAC(7, 2)
      LOGICAL :: ACCRET_COMPOSITION    ! Obsolete switch for composition accretion
      COMMON /ACCRET/ X1AC, X4AC, X12AC, X14AC, X16AC, X20AC, X24AC, XAC
      ! Namelist for the new init.dat I/O
      NAMELIST /INIT_DAT/ KH2, KR1, KR2, JCH, KTH, KX, KY, KZ, 
     & KCL, KION, KAM, KOP, KCC, KNUC, KCN,
     & KT1, KT2, KT3, KT4, KT5, KSV,
     & EPS, WANTED_EPS, DEL, DH0, CDC_MS, CDC_HEC, CDC_HES, CDC_DBLSH, CDC5, 
     & CDC_EMS, CDC_HG, CDC_1DUP, CDC_RLOF, CDC_RLOF_REDUCE,
     & KE1, KE2, KE3, KBC, KEV, KFN, KL, JH, KP_VAR, KP_EQN, KP_BC, 
     & KE1_2, KE2_2, KE3_2, KBC_2, KEV_2, KFN_2, KL_2, JH_2, KP_VAR_2, KP_EQN_2, KP_BC_2,
     & KSX, KN, KJN, 
     & CT1, CT2, CT3, CT, 
     & CH, CC, CN, CO, CNE, CMG, CSI, CFE, 
     & CALP, CU, COS, CPS, CRD, CXB, CGR, CEA, CET,
     & CMT, CMS, CMI, CMR, CMJ, CMV, CMK, CMNL, CMRR, CMVW, CMSC, CMW,
     & CMAL, SMART_MASS_LOSS, CML, CHL, CTF, CLT, CMDOTROT_HLW, CMDOTROT_MM,
     & CMTEL, CMTWL,
     & CPA, CBR, CSU, CSD, CDF, CGW, CSO, CMB, CPHOTONTIRE, UNUSED1,
     & CTH, UNUSED, USE_FUDGE_CONTROL, ALLOW_EXTENSION, ALLOW_UNDERRELAXATION,
     & ALLOW_OVERRELAXATION, CONVECTION_SCHEME, ALLOW_EGENRELAXATION,
     & USE_PREVIOUS_MU, OFF_CENTRE_WEIGHT, ALLOW_MDOTRELAXATION,
     & USE_SMOOTH_REMESHER, RELAX_LOADED_MODEL, CONVECTION_LEDOUX,
     & ALLOW_AVMURELAXATION, USE_QUADRATIC_PREDICTIONS, CLIMIT,
     & CMI_MODE, ZSCALING_MDOT, CDSI, CSHI, CSSI,
     & CESC, CGSF, CFMU, CFC, ARTMIX, CGRS, CCAC, CSMC
    
      ! Namelist for the accretion of matter
      NAMELIST /ACCRET/ X1AC, X4AC, X12AC, X14AC, X16AC, X20AC, X24AC, ACCRET_COMPOSITION

      ! Namelist for initial (ZAMS) nucleosynthesis abundances
      NAMELIST /ABUND/ CXD, CXHE3, CXLI7, CXBE7, CXB11, CXC12, CXC13,
     &CXC14, CXN14, CXN15, CXO16, CXO18, CXO17, CXF19, CXNE21, CXNE20,
     &CXNE22, CXNA22, CXNA23, CXMG24, CXMG25, CXMG26, CXAL26M, CXAL27,
     &CXAL26G, CXSI28, CXSI30, CXSI29, CXP31, CXS33, CXS32, CXS34, CXFE57,
     &CXFE60, CXFE56, CXFE58, CXFE59, CXCO59, CXNI61, CXNI59, CXNI58,
     &CXNI60

      ! FORMAT specifier for the old init.dat format
  993 FORMAT (8I5, /, 7I5, /, 6I5, /, 1P, 8D8.1, 0P, /, 2(10I4, /,
     & 6(20I3, /)), 3(15I3, /), I3, /, 2(20I3, /), 10F5.2, 1P, 3D8.1, /,
     & 0P, 7F6.3, /, 1P, 5(9D9.2, /), 0P)

      ! Default values
      CH_OPAC = CH
      CPHOTONTIRE = 0.0D0  ! Do not keep track of photon tiring (do not take kinetic energy of the wind from the luminosity)
      CMV = 0.0D0    ! Vink mass loss, disabled by default
      CMK = 0.0D0    ! Kudritzki 2002 mass loss, disabled by default
      CMNL = 0.0D0   ! Nugis&Lamers, for WR stars, disabled by default
      CH = -1.0D0    ! Use default value from opacity table 
      ARTMIX = 0.0D0 ! No artificial mixing
      CGRS  = 0.0D0  ! No gravitational settling/diffusion
      CCAC = 0.0D0   ! Default (single star): don't accrete composition
      CDC(:) = 1.0D0 ! Default timestep control
      CDC_EMS = 1.0d0! Timestep parameter at end of MS
      CDC_HG = 1.0d0 ! Timestep parameter in HG
      CDC_1DUP = 1.0d0! Timestep parameter at 1DUP
      CDC_RLOF = 1.0d0! Timestep parameter if star is close to RLOF
      CDC_RLOF_REDUCE = 1.0d0! Timestep if star is close to RLOF but contracting
      CSMC = 0.04d0  ! Efficiency of semi-convection
      CONVECTION_LEDOUX = 0.0d0  ! Default: Schwarzschild convection
      CMDOTROT_HLW = 0.0d0 ! Enhanced mass loss from Heger, Langer & Woosely
      CMDOTROT_MM = 0.0d0  ! Enhanced mass loss from Maeder & Meynet
      CMTEL = 0.0d0  ! No Eddington-limited accretion
      CMTWL= 0.0d0   ! No angular momentum limited accretion

! Improved mass loss rates, can be used in a general mass loss recipe
! instead of the current de Jager rate. These individual mass loss
! recipes can be switched on and off there.
      SMART_MASS_LOSS = 0.0
      CMRR = multiplier_reimers
      CMVW = multiplier_vasiliadis_wood
      CMSC = multiplier_schroeder
      CMW = multiplier_wachter
      CMAL = multiplier_achmad

      IF (KTW == 2) CCAC = 1.0D0 ! Default (binary/TWIN): do accrete composition
      X1AC = -1.0
      X4AC = -1.0
      X12AC = -1.0
      X14AC = -1.0
      X16AC = -1.0
      X20AC = -1.0
      X24AC = -1.0
      CMI_MODE = 1
! First, try to read the NAMELIST formatted init.dat
      READ (IT, NML=INIT_DAT, IOSTAT=IOERROR)
      KFN = NFUNC    ! Just copy all functions; common error to forget. 
      KFN_2 = NFUNC  ! Just copy all functions; common error to forget. 
      IF (KE2_2 /= 0) NUCLEOSYNTHESIS_ENABLED = .TRUE.
      IF (CH < 0.0) CH = CH_OPAC    ! Use CH from opacity table
      IF (USE_PREVIOUS_MU) AVMU_SMOOTH = 0.0 ! Use mu from prev. timestep
      IF (WANTED_EPS > EPS) WANTED_EPS = EPS
! Turn on mass loss recipes in smart mass loss routine, pass on the
! exponent for metallicity scaling that is to be used.
      multiplier_dejager = CMJ
      multiplier_schroeder = CMSC
      multiplier_reimers = CMRR
      multiplier_vasiliadis_wood = CMVW
      multiplier_wachter = CMW
      multiplier_achmad = CMAL
      multiplier_vink = CMV
      multiplier_kudritzki = CMK
      multiplier_nl = CMNL
      METALLICITY_SCALING_EXPONENT = ZSCALING_MDOT
! Store a copy of the mesh spacing function coefficients, we want to
! dynamically change them during some evolutionary stages and we need to
! be able to restore the original settings.
      INITIAL_CT(:) = CT(:)
      IF (IOERROR == 0) THEN
         READ_INIT_DAT = .TRUE.
         
! Allocate memory for storage of previous models in parabola extrapolation
         IF (USE_QUADRATIC_PREDICTIONS) CALL INITIALISE_PARABOLA_STORAGE_SPACE(KH2, KE1+KE2+KEV)

! Export some local names to their global counter parts
         CDC(1) = CDC_MS 
         CDC(2) = CDC_HEC
         CDC(3) = CDC_HES
         CDC(4) = CDC_DBLSH
         CDC(5) = CDC5
         CDC(6) = CDC_EMS
         CDC(7) = CDC_HG
         CDC(8) = CDC_1DUP
         CDC(9) = CDC_RLOF
         CDC(10) = CDC_RLOF_REDUCE
         
! Set some leftover junk that's never actually used
         UNUSED1 = 0.0
         UNUSED(:) = 0.0D0
         CQ1 = UNUSED1
         CQ(:) = UNUSED(:)
         ! Try to get accretion information from the same file
         ! This is a bit ugly, but do we want to use another fort.nnn file for
         ! this?
         READ (IT, NML=ACCRET, IOSTAT=IOERROR)

         ! Finally, read the nucleosynthesis abundances
         READ (IT, NML=ABUND, IOSTAT=IOERROR)
         
         RETURN
      END IF
      
! Fall back to the `old' bare numbers init.dat; close file to seek from start
      REWIND (IT)
      READ  (IT, 993, IOSTAT=IOERROR) KH2, KR1, KR2, JCH, KTH, KX, KY, KZ,
     & KCL, KION, KAM, KOP, KCC, KNUC, KCN, KT1, KT2, KT3, KT4, KT5, KSV, 
     & EPS, DEL, DH0, CDC_MS, CDC_HEC, CDC_HES, CDC_DBLSH, CDC5,
     & KE1, KE2, KE3, KBC, KEV, KFN, KL, JH, KP_VAR, KP_EQN, KP_BC, 
     & KE1_2, KE2_2, KE3_2, KBC_2, KEV_2, KFN_2, KL_2, JH_2, KP_VAR_2, KP_EQN_2, KP_BC_2,
     & KSX, KN, KJN, CT1, CT2, CT3, CT, 
     & CC, CN, CO, CNE, CMG, CSI, CFE, 
     & CALP, CU, COS, CPS, CRD, CXB, CGR, CEA, CET, 
     & CMT, CMS, CMI, CMR, CMJ, CML, CHL, CTF, CLT, 
     & CPA, CBR, CSU, CSD, CDF, CGW, CSO, CMB, UNUSED1,
     & CTH, UNUSED
      KFN = NFUNC    ! Just copy all functions; common error to forget. 

      INITIAL_CT(:) = CT(:)
! Export some local names to their global counter parts
      CQ1 = UNUSED1
      CQ(:) = UNUSED(:)
      CDC(1) = CDC_MS 
      CDC(2) = CDC_HEC
      CDC(3) = CDC_HES
      CDC(4) = CDC_DBLSH
      CDC(5) = CDC5
      WANTED_EPS = EPS

      IF (IOERROR == 0) THEN
         IF (USE_QUADRATIC_PREDICTIONS) CALL INITIALISE_PARABOLA_STORAGE_SPACE(KH2, KE1+KE2+KEV)
         READ_INIT_DAT = .TRUE.
      ELSE
         ! Error during read operation
         READ_INIT_DAT = .FALSE.
      END IF
      ! Rewind the file, in case we want to read it again later
      REWIND (IT)
      END FUNCTION

      SUBROUTINE BEGINN ( JOP, JIP, DTY, KP, KR1, KR2, KSV, KT5, IT, JO, JF )
      USE MESH
      USE MESH_ENC
      USE EXTRA_ELEMENTS
      USE CONSTANTS
      USE SETTINGS
      USE MODEL_INITIALISER
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON H(NVAR,NM), DH(NVAR,NM), EP(3), KH, KTW, KW(260)
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /TVBLES/ DT, ZQ(9), AGE, BM, MC(2), OM, PER, SM, ENC,
     : TC(64), PR(81), PPR(81), JHOLD, JM2, JM1
      COMMON /STORE / HPR(NVAR,NM), MS(9999), ST(9999),
     &  SDT(9999), SCM(9999), SANG(9999), SE(9999), WW(16), WX(15)
      DOUBLE PRECISION :: X1AC, X4AC, X12AC, X14AC, X16AC, X20AC, X24AC
      DOUBLE PRECISION :: XAC(7, 2)
      COMMON /ACCRET/ X1AC, X4AC, X12AC, X14AC, X16AC, X20AC, X24AC, XAC
      LOGICAL :: READ_INIT_DAT
      LOGICAL :: STATUS
      INTEGER :: IOERROR, KEQ
      CHARACTER*500 :: JIP_NAME
  993 FORMAT (8I5, /, 7I5, /, 6I5, /, 1P, 8D8.1, 0P, /, 2(10I4, /,
     & 6(20I3, /)), 3(15I3, /), I3, /, 2(20I3, /), 10F5.2, 1P, 3D8.1, /,
     & 0P, 7F6.3, /, 1P, 5(9D9.2, /), 0P)
  995 FORMAT (1X, 1P, 2D14.6, D17.9, 5D14.6, 0P, 6I5)
      H(:, 1:KH) = 0.0D0
      DH(:, 1:KH) = 0.0D0
C      ZQ(1:81) = 0.0D0  ! May generate index out-of-bound warning - replaced by explicit statements below - SdM
      ZQ(1:9) = 0.0D0      
      AGE = 0.0D0 
      BM = 0.0D0 
      MC(1:2) = 0.0D0 
      OM = 0.0D0 
      PER = 0.0D0 
      SM = 0.0D0 
      ENC = 0.0D0 
      TC(1:64) = 0.0D0  !was ZQ(18:91)

      ST(1:9999) = 1.0D20
      KB = 1
C Read miscellaneous data, usually unchanged during one evol run      
      ! Try to read in init.dat; abort in case this fails
      STATUS = READ_INIT_DAT(IT, KH2, KR1, KR2, KSV, KT5, JCH)
      IF (STATUS .EQV. .FALSE.) THEN
         INQUIRE(UNIT=IT, NAME=JIP_NAME)
         WRITE(0, *) 'Error reading init.dat data from "', TRIM(JIP_NAME),'"'
         STOP
      END IF

!     Just constructed a ZAHB model; DON'T homogenise the model,
!     whatever the input .dat file says.
!     Without this flag, the post-he flash model will have its
!     composition reset to ZAMS composition if we started from a
!     homogeneous model - which is not what we want.
!     We should also force the accretion mode to exponential.
      IF (JO == 13) THEN
         JCH = MIN(3, JCH)
         JO = -1
         CMI_MODE = 1
      END IF

      IF (MUTATE .AND. KH>0) THEN
         ! Don't rezone the model
         KH2 = KH
      END IF
      
! Read opacity data and construct splines 
! KOP (read from init.dat) sets which type of opacity tables
      IF (KOP<=1) THEN
         ! Iglesias & Rogers (1992), as implemented by Pols & al. 1995
         CALL LOAD_OPACITY(20)
      ELSE
         ! Iglesias & Rogers (1996), as implemented by Eldridge & Tout 2003
         CALL LOAD_OPACITY_CO(41)
      ENDIF

      ! Abort if the requested number of meshpoints is larger than the
      ! size of the array
      IF (KH2 > NM) THEN
         WRITE (0, *) 'Cannot rezone to ', KH2, 'meshpoints. Maximum size is ', NM, 'meshpoints.'
         JO = -2
         RETURN
      END IF

      ! Autodetect if we should solve for Mg24 or not by checking if the
      ! corresponding equation is in the list of equations
      USE_MG24_EQN = .FALSE.
      DO II = 51, 100
         IF (KW(II) == EMG24 .OR. KW(II) == ESUMX) THEN
            USE_MG24_EQN = .TRUE.
            EXIT
         ENDIF
         IF (KW(II) == 0) EXIT   ! Break loop if end of list found
      ENDDO

      ! Detect whether rotation is treated as solid body rotation or
      ! whether we consider differential rotation. The switch is on whether
      ! or not the rotational period (var. 13) is listed as an eigenvalue.
      RIGID_ROTATION = .TRUE.
      KEQ = KW(1)+KW(2)    ! Total number of first+second order equations
      DO II = 11, 10 + KEQ
         IF (KW(II) == 13) THEN  ! Rotational pertiod not an EV
            RIGID_ROTATION = .FALSE.
            EXIT
         END IF
      END DO
      
C Read data for initial model (often last model of previous run)
C e.g. SM = stellar mass, solar units; DTY = next timestep, years	
      II = 1
 11   READ  (JIP, *, IOSTAT=IOERROR)
     : SM, DTY, AGE, PER, BMS, ECC, P1, ENC, KH, KP, JMOD, JB, JIN, JF
      IF ( IOERROR /= 0 )THEN
         INQUIRE(UNIT=JIP, NAME=JIP_NAME)
         WRITE(0, *) 'Error reading "', TRIM(JIP_NAME),'"'
         STOP
      END IF
      IF ( JMOD == 0 .AND. KP>0 ) KP = KP + 1
      IF ( JIP.EQ.13 .OR. JIP.EQ.14 ) DTY = CT3*DTY
      IF ( II.EQ.2 ) JB = 1
      IF ( JMOD.LT.10 .AND. II.EQ.1 ) WRITE (JB, 993) 
     & KH2, KR1, KR2, JCH, KTH, KX, KY, KZ, KCL, KION, KAM, KOP, KCC, KNUC,
     & KCN, KT1, KT2, KT3, KT4, KT5, KSV, EP, CDC(1:5), KW, KSX, 
     & KN, KJN, CT1, CT2, CT3, CT,
     & CC, CN, CO, CNE, CMG, CSI, CFE,
     & CALP, CU, COS, CPS, CRD, CXB, CGR, 
     & CEA, CET, CMT, CMS, CMI, CMR, CMJ, 
     & CML, CHL, CTF, CLT, CPA, CBR, CSU, CSD, CDF, CGW, CSO, CMB, 
     & CQ1, CTH, CQ
      WRITE (JB, 995) SM, DTY, AGE, PER, BMS, ECC, P1, ENC, KH,
     : KP, JMOD, JB, JIN, JF
C Read the initial model 
      DO IK = 1, KH
         READ (JIP, *, IOSTAT=IOERROR) (H(IJ,IK), IJ = 1, JIN)
         IF ( IOERROR /= 0 )THEN
            INQUIRE(UNIT=JIP, NAME=JIP_NAME)
            WRITE(0, *) 'Error reading "', TRIM(JIP_NAME),'"'
            STOP
         END IF
      END DO
! Read DH(:) if the JF flag says it's available (bit 2 is set)
      IF (IAND(JF, 4)==4) THEN
         JF = JF - 4            ! Unset the flag
         DO IK = 1, KH
            READ (JIP, *, IOSTAT=IOERROR) (DH(IJ,IK), IJ = 1, JIN)
            IF ( IOERROR /= 0 )THEN
               INQUIRE(UNIT=JIP, NAME=JIP_NAME)
               WRITE(0, *) 'Error reading "', TRIM(JIP_NAME),'"'
               STOP
            END IF
         END DO
      ENDIF
      TM = CMSN*SM
      ! Decrease timestep for mutation runs
      IF ( MUTATE .AND. JMOD>0 ) DTY = MIN(DTY, 1.0D3)
      ! Make sure the timestep is small when starting ZAHB construction
      IF (IT == 24) THEN
         AGE = 0.0D0
         DTY = MIN(DTY, 1.0D3)
      END IF
C Convert some things to `cgs' units: 10**11 cm, 10**33 gm, 10**33 erg/s
      IF ( II.EQ.2 ) GO TO 8
      DT = CSY*DTY 
      CMI = CMI/CSY
      CMS = CMS/CSY
      CMT = CMT*1.0D-11
      IF ( CEA.LE.0.0D0 ) CEA = 1.0D-10
c Initialise orbital angular momentum
      BM = CMSN*BMS
      OM = BM - TM
      OA = CG1*TM*OM*(CG2*PER/BM)**C3RD*DSQRT(1.0D0 - ECC*ECC)
      IF ( JB.EQ.1 ) GO TO 8
c For *2 of a binary, read in the mass-loss history of *1
      DO IK = 1, 9999
         READ (3, 102, END=8) MS(IK), WW, SDT(IK), WW, ST(IK), SE(IK),  
     :     WX, SCM(IK), WW,  SANG(IK), WW, IR
         IF ( IK.GT.1 .AND. ST(IK).GT.1.0D19 ) GO TO 8
      END DO
  102 FORMAT (3(1P, D16.9, 5D10.3, 0P, F8.3, 7F8.5, 3F8.4, /),
     :         (1P, D16.9, 5D10.3, 0P, F8.3, 7F8.3, 3F8.4, /),
     :          1P, D16.9, 5D10.3, 0P, 11F8.3, I6)
 8    AGE = AGE - DTY
      CALL NEXTDT ( DTY, JO, IT )
C optionally rezone the model, e.g. for different no. of meshpoints.
c also initialise some variables that were not in input model.
      IF (USE_SMOOTH_REMESHER) THEN
      !   CALL REMESH ( KH, JCH, BM, H(4,1), H(13,1), ECC, OA, II, JF )
         CALL REMESHER ( KH2, JCH, BM, TM, P1, ECC, OA, II, JF )
      ELSE
         CALL REMESH ( KH2, JCH, BM, TM, P1, ECC, OA, II, JF )
      END IF
!     Load the nucleosynthesis data structures
      CALL ALLOCATE_NUCLEOSYNTHESIS_DATA(KH)
      CALL SET_INITIAL_NUCLEOSYNTHESIS_ABUNDANCES

! ZAHB construction, we need to take special care here since the ZAHB
! starting model will not have the correct envelope abundance. This is
! corrected by accreting material of the right composition, but we
! cannot do this immediately because we need to H burning shell to stay
! active and dumping pure He on top of it will make it less numerically
! stable.
! Composition accretion is switched on in printb when the stellar mass
! reaches the desired core mass.
      IF (IT == 24) CCAC = 0.0
      IF ( II.EQ.KTW ) GO TO 12
! Backup the model in case we need to go back a timestep
! Backup the extra variables
      IF (NXVAR>0) HPR(41:NVAR, 1:KH) = H(41:NVAR, 1:KH)
      HPR(1:16, 1:KH) = H(1:16, 1:KH)

      II = 2
      GO TO 11
 12   IF ( II.EQ.1 ) GO TO 14
      H(25:40, 1:KH) = H(1:16, 1:KH)
      HPR(25:40, 1:KH) = H(1:16, 1:KH)
      H(1:16, 1:KH) = HPR(1:16, 1:KH)   

c store some numbers for possible restart with BACKUP
 14   JHOLD = 2
C     PR(1:81) = ZQ(1:81) ! May generate index out-of-bound warning
C     Replaced by explicit statement belwo -SdM
      PR(1:81) = (/ ZQ(1:9), AGE, BM, MC(1:2), OM, PER, SM, ENC, TC(1:64) /)      
      PPR(1:81) = PR(1:81)
      JM1 = JMOD
c Determine whether I and phi are computed or not, for OUTPUT
      JF = 0
      DO I = 11, 50
         IF ( KW(I).EQ.12 .OR. KW(I).EQ.14 ) JF = JF + 1
      ENDDO
      IF ( JF.EQ.1 ) JF = 0
    4 RETURN
      END

      SUBROUTINE FGB2HB ( JOP, JJB, JO3 )
      USE MESH
      USE MESH_ENC
      USE FUDGE_CONTROL
      USE FGB2HB_COMPOSITION
c Evolve a standard ZAHB model (stored on fort.12, with init.dat in fort.24)
c to the required ZAHB model: same total mass (SM) and core mass (VMH)
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON H(NVAR,NM), DH(NVAR,NM), EP(3), KH, KTW, ID(260)
      COMMON /QUERY / ML(3), UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /TVBLES/ ZQ(16), SM, ENC(24), VMH, VME(202), JHOLD(3)
      COMMON /ABUND / XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, XW(14)
      INTEGER :: CMI_MODE_BCK
      UC(13) = SM
! FIXME: the core mass must be >= the mass of the ZAHB construction model
! Probably doen't matter much in practice.
      UC(14) = MAX(0.40, VMH)
      REWIND (JOP)
      CMI_MODE_BCK = CMI_MODE
      CMI_MODE = 1
      CALL STORE_PRE_FLASH_COMPOSITION
      CALL STAR12 ( JO3, JC1, JOP, 12, KSV, 24 ) 
      CALL CLEANUP_PRE_FLASH_COMPOSITION
      CMI_MODE = CMI_MODE_BCK
      UC(13) = 1.0D3
      UC(14) = 1.0D3
      REWIND (22)
      REWIND (12)
      REWIND (24)
      REWIND (JOP)
      RETURN
      END

      SUBROUTINE MS2BSS ( JOP, JO3 )
c Convert a normal main sequence star to a different type of star, say
c a collision remnant.
      USE MESH
      USE MESH_ENC
      USE FUDGE_CONTROL
      USE CONSTANTS
      IMPLICIT NONE
      INTEGER :: JOP, JO3, ISB, IP1, IM1, IP2, IM2, KPT, KP
      DOUBLE PRECISION :: ML(3), UC(21), DT, ZQ1(9), AGE, ZQ2(4), SM, PER, 
     & ENC(24), VMH, VME(202)
      INTEGER :: JMOD, JB, JNN, JTER, JOC, JKH, JHOLD(3), JC1, KSV, JIP, JOP2
      INTEGER :: KTH(101), KH, KTW, KW(260)
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EP(3)
      DOUBLE PRECISION :: SM0, DTY0, AGE0, PER0, BMS0, ECC0, P0, ENC0
      DOUBLE PRECISION :: SM1, DTY1, AGE1, PER1, BMS1, ECC1, P1, ENC1
      INTEGER :: S_KPT
      COMMON H, DH, EP, KH, KTW, KW
      COMMON /QUERY / ML, UC, JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /TVBLES/ DT, ZQ1, AGE, ZQ2, PER, SM, ENC, VMH, VME, JHOLD
      COMMON /TN1/ SM0, DTY0, AGE0, PER0, BMS0, ECC0, P0, ENC0
      COMMON /T0 / SM1, DTY1, AGE1, PER1, BMS1, ECC1, P1, ENC1
      COMMON /SAVEINIT/ S_KPT

      ! First stage: Evolve a normal main sequence star to the right mass and
      !  core hydrogen content. Read settings from the normal init.dat
      CALL READ_TARGET_MODEL()
      UC(15) = GET_IPH(0.0D0, 5)
      REWIND (JOP)
      JIP = JOP
      CALL STAR12 ( JO3, JC1, JOP, JIP, KSV, 22 )
      REWIND (22)
      REWIND (JOP)
      UC(15) = -1
      IF ( JO3 == 51 ) WRITE (JB, *) 'Begin mutation'
      
      ! Termination code should be 51 (end of MS)
      IF (JO3 /= 51) RETURN

      ! Second stage: mutate the model by tweaking the composition profiles and
      !  puting an artificial contribution to luminosity in the energy equation
      ! Stop when we have matched the target radius
      
      ! Enable some flags that control various parts of the mutation process
      ADJ_MEA = .TRUE.
      ADJ_COMP = .TRUE.
      USEMENC = .TRUE.
      CURR_DIFFSQR = 1.0D3;
      JIP = 27 - JIP
      JO3 = -1
      BEST_DIFFSQR = CURR_DIFFSQR
      MUTANT_H(1:24, 1:KH) = H(1:24, 1:KH)
      REWIND (JIP)
      CALL STAR12 ( JO3, JC1, JOP, JIP, KSV, MUTATE_DAT )
      REWIND (MUTATE_DAT)
      REWIND (JOP)
      REWIND (JIP)
      
      IF (JO3 == 0 .OR. JO3 == 2 .OR. JO3 == 5 .OR. JO3 == 12 .OR. JO3 == 15) JO3 = 53
      
      IF (JO3 == 53) THEN
         H(1:24, 1:KH) = MUTANT_H(1:24, 1:KH)
         DT = 1.0D1*CSY
      END IF
      
      H(NMENC, 1:KH) = 0.0D0
      H(NMEA, 1:KH) = 0.0D0
      H(NMET, 1:KH) = 0.0D0
      
      ! Re-read some values from init.dat
      !READ (23, *) ISB, KTW, IP1, IM1, IP2, IM2, KPT, KP
      !CLOSE (23)
      AGE = 0
      
      ! Rotational period; set it as when it was read in from the target model
      ! Optionally replace with the value from init.run
      IF (P1 >= 0.0) THEN
         H(13, 1:KH) = P1
      ELSE
         H(13, 1:KH) = P0
      END IF
      CALL OUTPUT ( S_KPT, 14, 0, 0 )
      CALL OUTPUT ( S_KPT, 13, 0, 0 )
      REWIND (14)
      REWIND (13)
      CALL OUTPUT ( S_KPT, 65, 0, 0 )
      CLOSE (65)

      ! Termination code should be 53 (Found convergence minumum)
      IF (JO3 /= 53) RETURN
      
      WRITE (JB, *) 'Mutation done', BEST_MOD, BEST_DIFFSQR
      ! Third stage: evolve the remnant
      ! Let the normal code handle this
      JO3 = -1
      MUTATE = .FALSE.
      ADJ_MEA = .FALSE.
      ADJ_COMP = .FALSE.
      USEMENC = .FALSE.
      RETURN
      END

      SUBROUTINE NEXTDT ( DTY, JO, IT )
      USE MESH
      USE COMPARE_FLOATS;
      USE CONSTANTS
      USE SETTINGS
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON H(NVAR,NM), DH(NVAR,NM), EP(3), KH, KTW, KW(260)
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /TVBLES/ DT, ZQ(9), AGE, BM(10), T0, M0, MTA, OM0, OMTA, 
     : A0, ATA, E0, ETA, CDD, BP(213), JHOLD, JM1(2)
      COMMON /STORE / HPR(NVAR,NM), MS(9999), ST(9999),
     &  SDT(9999), SCM(9999), SANG(9999), SE(9999), WW(16), WX(15)
c Find change from last timestep
      SUM = 0.0D0
      DO IK = 1, KH
         DO IJ = 1, KN
            SUM = SUM + DABS(DH(KJN(IJ), IK))
         END DO
      END DO
      IF ( feq(SUM, 0.0D0) ) FAC = 1.0D0
      IF ( SUM.GT.0.0D0 ) FAC = SUM/(CDD*KH)
!      IF ( IT.NE.12 ) WRITE (JB, *) FAC    !MvdS: Commented disturbing line in output
      IF ( feq(SUM, 0.0D0) ) GO TO 9
!      WRITE (JB, *) SUM/KH 
c Specify next DTY; simple, for *1
      IF ( JHOLD.GT.3 ) DTY = DMAX1( CT1, DMIN1(CT2, 1.0D0/FAC))*DTY
! Limit DTY in such a way that we don't exceed the age specified in the
! input file (ignore for ZAHB construction).
! In this way, we get a model of *exactly* a specific age.
! Control the timestep so that the approach to the final age is a bit smooth
      IF (IT /= 24) THEN
      DTY_MIN = UC(12)/CSY
      DTY_MAX = UC(2) - AGE
      IF ( AGE+2*DTY < UC(2) .AND. AGE+3*DTY >= UC(2)) THEN
! We expect three more timesteps, start constraining the timestep
         DTY = 0.4*MAX(DTY_MAX, DTY_MIN )
      END IF
      IF ( AGE+DTY < UC(2) .AND. AGE+2*DTY >= UC(2)) THEN
! We expect maybe two more timesteps, constrain
         DTY = 0.6*MAX(DTY_MAX, DTY_MIN )
      END IF
      IF ( AGE+DTY >= UC(2) .AND. AGE>UC(12)) THEN
! This is our final timestep
         DTY = MIN(DTY_MAX, DTY )
      END IF
      IF ( DTY_MAX <= DTY_MIN ) JO = 5
      END IF

    9 IF ( JB.EQ.1 ) GO TO 39
      KB = 1
      IF ( AGE.LE.ST(KB) ) GO TO 81
c For *2, find KB such that ST(KB) < AGE <= ST(KB + 1)
      DO WHILE ( .NOT. (AGE.GT.ST(KB) .AND. AGE.LE.ST(KB + 1)) ) 
         KB = KB + 1
      END DO
 81   KA = KB + 1
      IF ( AGE + DTY.LE.ST(KB) ) GO TO 85
c For *2, find KA such that ST(KA - 1) < AGE + DTY <= ST(KA)
 84   IF ( AGE + DTY.GT.ST(KA - 1) .AND. AGE + DTY.LE.ST(KA) ) GO TO 85
      KA = KA + 1
      GO TO 84
 85   JD = 2.0D0 + 1.0D0/(FAC + 0.1D0)
      IF ( KA.LE.KB + JD ) GO TO 87
c If KA >> KB, pull back a bit
      KA = KB + JD
      DTY = 0.5D0*(ST(KA) + ST(KA - 1)) - AGE
      DH(:, 1:KH) = 0.0D0
 87   IF ( ST(KA).GT.1.0D19 ) JO = 3
      WRITE (2, 962) KB, KA, JD, AGE, DTY, AGE + DTY, ST(KB), ST(KA)
 962  FORMAT (3I5, 1P, 10D16.8)
   39 DT = CSY*DTY
      IF ( feq(CT1, 1.0D0) .AND. feq(CT2, 1.0D0) ) GO TO 6
      IF ( JHOLD.GT.3 .AND. FAC.LT.2.0D0 .AND. JNN.GT.2 ) GO TO 6
c clear DH in some circumstances
      DH(:, 1:KH) = 0.0D0
 6    IF ( JB.EQ.1 .OR. JO.EQ.3 ) RETURN
c For *2, some factors relating to accretion from *1. Ignored if this is *1
      IF ( KA.EQ.1 ) KA = 2
      IPS = KA - 1
      T0 = CSY*ST(IPS)
      M0 = CMSN*MS(IPS)
      FDT = ST(KA) - ST(IPS)
      IF ( FDT.LE.0.0D0 ) JO = 11
      TDF = 0.0D0
      IF ( FDT.GT.0.0D0 ) TDF = 1.0D0/(FDT*CSY)
      MTA = CMSN*(MS(KA) - MS(IPS))*TDF
      OM0 = CMSN*SCM(IPS)
      OMTA = CMSN*(SCM(KA) - SCM(IPS))*TDF
      A0 = SANG(IPS)
      ATA = (SANG(KA) - SANG(IPS))*TDF
      E0 = SE(IPS)
      ETA = (SE(KA) - SE(IPS))*TDF
      RETURN
      END

      SUBROUTINE UPDATE ( DTY )
      USE MESH
      USE MESH_ENC
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON H(NVAR,NM), DH(NVAR,NM), EP(3), KH, KTW, KW(260)
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /TVBLES/ DT, ZQ(9), AGE, BM(71), PR(81), PPR(81), JHOLD,
     : JM2, JM1
      COMMON /STORE / HPR(NVAR,NM), MS(9999), ST(50026)
c Store certain current and previous values, for possible emergency backup
      AGE = AGE + DTY
      JMOD = JMOD + 1
      PPR(1:81) = PR(1:81)
!      PR(1:81) = ZQ(1:81)  ! May generate index out-of-bound warning
      PR(1:81) =  (/ZQ(1:9), AGE, BM(1:71)/)
      JM2 = JM1
      JM1 = JMOD
c Update the stored models: previous and anteprevious
      HPR(:, 1:KH) = H(:, 1:KH)
      H(:, 1:KH) = H(:, 1:KH) + DH(:, 1:KH)
! Update nucleosynthesis
      JHOLD = JHOLD + 1 
      RETURN
      END

      SUBROUTINE BACKUP ( DTY, JO )
      USE MESH
      USE MESH_ENC
      USE CONSTANTS
      USE SETTINGS
      USE NUCLEOSYNTHESIS
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON H(NVAR,NM), DH(NVAR,NM), EP(3), KH, KTW, KW(260)
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /TVBLES/ DT, ZQ(81), PR(81), PPR(81), JHOLD, JM2, JM1
      COMMON /STORE / HPR(NVAR,NM), MS(60025)
      ZQ(1:81) = PPR(1:81)
c *don't* backup DT!!
      JMOD = JM2
      DT = CT3*DT
      DTY = DT/CSY
      IF ( DT.LT.UC(12) .OR. JMOD.LE.2 .OR. CT3.GT.0.9999D0 ) JO = 2
      DH(:, 1:KH) = 0.0D0 
      H(:, 1:KH) = HPR(:, 1:KH)
      CALL BACKUP2
      JHOLD = -1
      RETURN
      END

      SUBROUTINE OUTPUT ( KP, JOP, JO, JF )
c Write an intermediate or final model to disc
      USE MESH
      USE MESH_ENC
      USE CONSTANTS
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON H(NVAR,NM), DH(NVAR,NM), EP(3), KH, KTW, KW(260)
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /TVBLES/ DT, ZQ(9), AGE, BM, MC(2), OM, PER, SM, ENC,
     : TC(18), SE, TN(207), JHOLD, JM2, JM1
! Save some parameters so that they can be restored after the He flash
      DOUBLE PRECISION, SAVE :: AJ
      INTEGER, SAVE :: JMD, KP_BCK
      INTEGER*4 I1(2), I2(2)
      DATA I1, I2 / 0, 24, 24, 16/
   92 FORMAT (1X, 1P, 40D23.15, 0P)
   95 FORMAT (1X, 1P, 8D23.15, 0P, 6I6)
   99 FORMAT (1X, 1P, I6, 4D23.15, 0P)
      IF ( JO.EQ.8 ) THEN        ! He flash
         ! Backup age, model number and number of models, to be restored after the flash
         AJ = AGE
         JMD = JMOD
         KP_BCK = KP
         RETURN   ! Comment this line to save post-He-flash model.
      ELSEIF ( JO.EQ.13 ) THEN
         ! Restore age and model number, increase the timestep
         AGE = AJ
         JMOD = JMD
         DT = 5.0D4*CSY
         ! The number of remaining models to calculate. Do at least a few.
         KP = MAX(KP_BCK - JMOD, 10)
         IF (KP_BCK < 0) KP = KP_BCK
      END IF
      IF ( MUTATE .AND. JMOD>0 ) DT = MIN(DT, 1.0D3*CSY)
      JP = JOP
      IF ( ( JOP.EQ.13.OR.JOP.EQ.14 ).AND.JO.NE.13 ) JP = 27 - JOP
      DO II = 1, KTW
         SMOUT = H(4 + I1(II), 1)/CMSN
         APER = 2.0*CPI/(DABS(H(13 + I1(II), 1)) * CSDAY)
         WRITE (JP, 95) SMOUT, DT/CSY, AGE, PER, BM/CMSN, SE,
     :         APER, ENC, KH, KP, JMOD, JB, I2(II), JF
         DO IK = 1, KH
            APER = 2.0*CPI/(DABS(H(13 + I1(II), IK)) * CSDAY)
            WRITE (JP, 92) (H(IJ + I1(II), IK), IJ = 1, 12), APER,
     &                     (H(IJ + I1(II), IK), IJ = 14, I2(II))
         END DO
      END DO
      CALL FLUSH(JP) !MvdS

      RETURN
      END
