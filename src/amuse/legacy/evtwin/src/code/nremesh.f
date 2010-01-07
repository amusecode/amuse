      module model_initialiser
      
      contains

! remesher:
! Function based on the origianl REMESH function
! Input (in COMMON blocks):
!  H(,) - array of independent variables (in unnamed COMMON)
!  DH(,) - list of changes in independent variables (in unnamed COMMON)
!  KH - current number of meshpoints in H(,)
! Input options:
!  KH2 - new number of meshpoints for this model
!  JCH - switch to determine whether to construct new mesh spacing function
!        and initialise composition or not.
!  BM  - Total mass in the binary (EV)
!  TM  - New mass of the star
!  P1  - New rotational period of the star (solid body)
!  ECC - New eccentricity of the binary
!  OA  - New orbital angular momentum
!  JSTAR - Labels which if the stars to initislise variables for
!  JF - Switch to decide which variables to recompute
! TODO: this could be split up into different subroutines for each of the
! individual tasks.
      subroutine remesher( KH2, JCH, BM, TM, P1, ECC, OA, JSTAR, JF )
      use constants
      use mesh
      use init_dat
      use settings
      use fudge_control
      use extra_elements
      implicit none
      integer, intent(in) :: KH2, JCH, JSTAR, JF
      double precision, intent(in) :: BM, TM, P1, ECC, OA
      integer :: nm_current, nm_next, nm_target
      integer :: ik, ih, i
      integer :: JO
      integer :: KR1, KR2, KSV, KT5
      type(init_dat_settings) :: initdat
      DOUBLE PRECISION :: NH(NVAR,NM), NDH(NVAR,NM), NNDH(NVAR,NM)
      DOUBLE PRECISION :: Q1, Q2, DK, DTY, MF, VD, DTB, AGEB, Pcrit, VMA, HPC
      DOUBLE PRECISION :: SI
      LOGICAL :: EQUILIBRIUM
      LOGICAL :: NEWMESH
! Unnamed COMMON block
      DOUBLE PRECISION :: H(NVAR,NM),DH(NVAR,NM),EPS,DEL,DH0
      INTEGER :: KH,KTW,ID(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, ID
! COMMON block INF
      DOUBLE PRECISION :: VAR(NVAR), DVAR(NVAR), 
     :     FILLUP_COMMONBLOCK_INF(NFUNC+NVAR*NFUNC)
      COMMON /INF   / VAR, DVAR,   FILLUP_COMMONBLOCK_INF
! COMMON block VBLES
      DOUBLE PRECISION :: LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QMU, WL, 
     : WCV, HP, WT, PHIMU, GMR, SEP, M3, PX(NPX*(NM+2)), QA(NM)
      COMMON /VBLES / LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QMU, WL, 
     : WCV, HP, WT, PHIMU, GMR, SEP, M3, PX, QA
! COMMON block TVBLES
      DOUBLE PRECISION :: DT, ZQ(1), HSPN(2), RLF(2), ZET(2),
     $     XIT(2), AGE,BBM, MC(2), OM, PER, SM, ENC, TC(2), TFR, T0, M0,
     $     MTA, OM0, OMTA,A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS,
     $     SE, TN(2), WMH,WMHE, MH, MHE, MCO, VMG, BE(2), LH, LHE, LC,
     $     LNU, LTH, MCB(8),MSB(6), RCB(8), TCT(8), QR(81), PPR(81)
      INTEGER :: JHOLD, JM2, JM1
      COMMON /TVBLES/ DT, ZQ, HSPN, RLF, ZET, XIT, AGE, 
     : BBM, MC, OM, PER, SM, ENC, TC, TFR, T0, M0, MTA, OM0, OMTA, 
     : A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS, SE, TN, WMH, 
     : WMHE, MH, MHE, MCO, VMG, BE, LH, LHE, LC, LNU, LTH, MCB, 
     : MSB, RCB, TCT, QR, PPR, JHOLD, JM2, JM1
! COMMON block STORE
      DOUBLE PRECISION :: HPR(NVAR,NM), MS(9999), ST(9999),
     &  SDT(9999), SCM(9999), SANG(9999), SSE(9999), WWW(16), WX(15)
      COMMON /STORE / HPR, MS, ST, SDT, SCM, SANG, SSE, WWW, WX
! COMMON block QUERY
      DOUBLE PRECISION :: ML, QL, XL, UC(21)
      INTEGER :: JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /QUERY / ML, QL, XL, UC, JMOD, JB, JNN, JTER, JOC, JKH
! COMMON block ATDATA
      DOUBLE PRECISION :: CH2(4), CHI(26,9), COM(27), CAN(9), CBN(9), KZN(9)
      COMMON /ATDATA/ CH2, CHI, COM, CAN, CBN, KZN
! COMMON block STAT2
      DOUBLE PRECISION :: WU(3), P, RHO, WS(36), EX, ENX, WMU, DELTA,
     & PHI, EXT, FKT, FKR, PRANDTL
      COMMON /STAT2 / WU, P, RHO, WS, EX, ENX, WMU, DELTA, PHI, EXT, FKT,
     & FKR, PRANDTL
! COMMON bloack ABUND
      DOUBLE PRECISION :: XH0, XHE0, XC0, XN0, XO0, XNE0, XMG0, XSI0, XFE0
      DOUBLE PRECISION :: CHE
      DOUBLE PRECISION :: XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, XW(14)
      COMMON /ABUND / XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, XW
! COMMON block ACCRET
      DOUBLE PRECISION :: X1AC, X4AC, X12AC, X14AC, X16AC, X20AC, X24AC
      DOUBLE PRECISION :: XAC(7, 2)
      COMMON /ACCRET/ X1AC, X4AC, X12AC, X14AC, X16AC, X20AC, X24AC, XAC

! Backup current settings so we can restore them when we're done
      DTB = DT
      AGEB = AGE
      call push_init_dat(initdat, KH2, KR1, KR2, KSV, KT5, JCH)
! Change settings to a reasonable set of defaults
      call load_basic_init_dat(ik, KR1, KR2, KSV, KT5, ik)
      KOP = initdat%KOP
      KX = 0; KY = 0; KZ = 0; KTH = 0
      CMI = 0.0
      CRD = 0.0D0
      KT1 = 100
      KT2 = 0
      KT3 = 0
      KT4 = 100
      KT5 = 100
      KSV = 100
      JOC = 1
      JTER = 0

! Set initial composition.
! The composition variables are NOT the actual mass fractions if we 
! use non-integer atomic masses, so we have to compute what they are
! The composition variables used in the code are baryon number fractions
! We calculate this even if we don't want to do a ZAMS run because we need to
! know the baryon number densities of Fe, Si and Mg.
! Note that the actual abundances of Si and Fe are forced to by non-zero.
! This is because the EoS becomes poorly defined if there are not enough free
! electrons. It has no impact on the opacity and only a small impact on the
! mean molecular weight and the mass loss rate.
      CHE = 1.0D0 - CH - CZS
      CN = 1.0D0 - CC - CO - CNE - CMG - CSI - CFE

      XH0 = CH*CBN(1)/CAN(1)
      XHE0 = CHE*CBN(2)/CAN(2)
      XC0 = CC*CZS*CBN(3)/CAN(3)
      XN0 = CN*CZS*CBN(4)/CAN(4)
      XO0 = CO*CZS*CBN(5)/CAN(5)
      XNE0 = CNE*CZS*CBN(6)/CAN(6)
      XMG0 = CMG*CZS*CBN(7)/CAN(7)
      XSI0 = CSI*CZS*CBN(8)/CAN(8)
      XFE0 = CFE*CZS*CBN(9)/CAN(9)

      XH = XH0
      XHE = XHE0
      XC = XC0
      XN = XN0
      XO = XO0
      XNE = XNE0
      XMG = XMG0
      XSI = MAX(XSI0, CSI*1.0D-4)
      XFE = MAX(XFE0, CFE*1.0D-4)

      VMA = XH + XHE + XC + XN + XO + XNE + XMG + XSI + XFE
      XH = XH / VMA
      XHE = XHE / VMA
      XC = XC / VMA
      XN = XN / VMA
      XO = XO / VMA
      XNE = XNE / VMA
      XFE = XFE / VMA
      XSI = XSI / VMA
      XMG = XMG / VMA

! Initialise composition variables
      IF ( JCH >= 4 ) THEN
         H( 5,1:KH) = XH
         H( 3,1:KH) = XO
         H( 9,1:KH) = XHE
         H(10,1:KH) = XC
         H(11,1:KH) = XNE
         H(16,1:KH) = XN
      END IF
! We should always do this for Mg24, since that's never stored
      IF (USE_MG24_EQN) THEN
         DO IK=1, KH
            H(21,IK) = MAX(0.0D0, 1.0 - (H(5,IK) + H(9, IK) + H(10, IK) + H(3, IK) + H(11, IK) + H(16, IK) + XFE + XSI))
         END DO
      END IF

! Now we must also convert the abundances of the accreted material. If not
! set from init.dat, set from initial abundances.
! FIXME: this will cause problems if we ever need to call REMESH twice in
! the same run
      X1AC = X1AC*CBN(1)/CAN(1)
      X4AC = X4AC*CBN(2)/CAN(2)
      X12AC = X12AC*CBN(3)/CAN(3)
      X14AC = X14AC*CBN(4)/CAN(4)
      X16AC = X16AC*CBN(5)/CAN(5)
      X20AC = X20AC*CBN(6)/CAN(6)
      X24AC = X24AC*CBN(7)/CAN(7)
      IF (X1AC < 0.0) X1AC = XH0
      IF (X4AC < 0.0) X4AC = XHE0
      IF (X12AC < 0.0) X12AC = XC0
      IF (X14AC < 0.0) X14AC = XN0
      IF (X16AC < 0.0) X16AC = XO0
      IF (X20AC < 0.0) X20AC = XNE0
      IF (X24AC < 0.0) X24AC = XMG0
      VMA = X1AC + X4AC + X12AC + X14AC + X16AC + X20AC + X24AC + XFE + XSI

      X1AC  = X1AC  / VMA
      X4AC  = X4AC  / VMA
      X12AC = X12AC / VMA
      X14AC = X14AC / VMA
      X16AC = X16AC / VMA
      X20AC = X20AC / VMA
      X24AC = X24AC / VMA
! make sure XH is 1-everything else and abundancies sum to 1
      X1AC = MAX(0.0D0, 1.0D0-(X4AC+X12AC+X14AC+X16AC+X20AC+X24AC+XFE0+XSI0))

      ! Initialise accretion abundances for both stars
      XAC(1, 1:2) = X1AC  
      XAC(2, 1:2) = X4AC  
      XAC(3, 1:2) = X12AC 
      XAC(4, 1:2) = X14AC 
      XAC(5, 1:2) = X16AC 
      XAC(6, 1:2) = X20AC 
      XAC(7, 1:2) = X24AC 

! Set initial values of some other variables
! Typical mass-scale for the interior (needed for mesh spacing function)
      MC(JSTAR) = TM          ! New mass after remesh
      VAR(:) = H(:, KH)
      DVAR(:) = 0.0D0
      CALL FUNCS1 ( KH, -2)
      HPC = DSQRT(P/(CG*RHO*RHO))
      MC(JSTAR) = 3.5D-33*RHO*HPC**3

! Initialise binary (orbital) parameters
      H(17, 1:KH) = OA                          ! New orbital angular momentum
      H(18, 1:KH) = ECC                         ! New eccentricity
      H(20, 1:KH) = BM                          ! New total binary mass

! First: change the number of meshpoints or the mesh spacing function
! Determine if we need to calculate a new mesh (independent of the number
! of meshpoints)
      NEWMESH = .FALSE.
      if (JCH>3) NEWMESH = .TRUE.
! Set new number of meshpoints
      nm_current = KH
      nm_target = KH2
      nm_next = nm_target
      print *, 'NREMESH from', nm_current, 'to', nm_target
      do while(NEWMESH .or. nm_current /= nm_next)
         NEWMESH = .FALSE.
         print *, 'Trying ', nm_next
! Store old model, so we can go back if needed
         NH(:,1:nm_current) = H(:,1:nm_current)
         NDH(:,1:nm_current) = DH(:,1:nm_current)
! Find values of mesh spacing function
         do ik=1, nm_current
            VAR(:) = H(:, IK)
            CALL FUNCS1 ( IK, -2 )                 ! Calculate stuff
            QA(IK) = QQ                            ! Store mesh spacing function
         end do
! Interpolate model onto new mesh
! Find values of mesh spacing function at the external points
! Needed to calculate the new mesh, where the meshspacing gradient is
! constant.
         Q1 = QA(1)
         Q2 = (nm_next - 1.0D0)/(QA(nm_current) - QA(1))
         do ik = 1, nm_current
            QA(IK) = (QA(IK) - Q1)*Q2 + 1.0D0   ! Adjust meshspacing
         end do
         ih = 1
         do ik = 1, nm_next
            DK = 0.0D0
            if ( ik == nm_next ) ih = nm_current
            if ( ik /= 1 .and. ik /= nm_next ) then
               ! Find the proper value for the meshspacing function at
               ! this meshpoint
               do i = 1, 50
                  ! Sanity check: abort if we're running out of the mesh
                  ! boundary
                  if ( ih+1 > nm_current) then
                     write (0, *) 
     &                'REMESH running outside mesh boundary, aborting'
                     stop
                  end if
                  if ( ik >= qa(ih + 1) ) ih = ih + 1
                  if ( ik < qa(ih + 1) ) exit   ! Break loop
               end do
               DK = (ik - QA(ih))/(QA(ih + 1) - QA(ih))
            END IF
            ! Linear interpolation for new H and DH
            H(:, ik) = NH(:, ih) + DK*(NH(:, ih + 1) - NH(:, ih))
            NNDH(:, ik) = NDH(:, ih) + DK*(NDH(:, ih + 1) - NDH(:, ih))
         END DO
         !H(6, 1:KH) = 1.0/Q2              ! Gradient of mesh spacing

! Now see if the model will converge properly if we let the code iterate
         DTY = DT/CSY
         JO = 0
         JNN = 0
         KH = nm_next
         CALL PRINTB ( DTY, JO, IK, DK, 1, 22 )
         AGE = AGE - DTY
         CALL NEXTDT ( DTY, JO, 22 )
         JNN = 1
         DH(:, 1:nm_next) = 0.0D0
         CALL SOLVER(20, ID, KT5, JO)
         IF (JO == 0) THEN
            print *, 'Converged ok'
! If yes, pick next number of meshpoints
            nm_current = nm_next
            nm_next = nm_target
            H(:,1:nm_current) = H(:,1:nm_current) + DH(:,1:nm_current)
         ELSE
! If no, pick a smaller number of meshpoints in between the current value
! and the target and try again.
            nm_next = (nm_current+nm_next)/2
            print *, 'Cannot converge, reduce to ', nm_next
! Restore backup copies of H and DH
            H(:,1:nm_current) = NH(:,1:nm_current)
            DH(:,1:nm_current) = NDH(:,1:nm_current)
         END IF
      end do
      print *, 'NREMESH finished with ', nm_current, '(wanted ', nm_target,')'
      if (nm_current < nm_target) then
         print *, '*** NREMESH failed ***'
         stop
      end if
      DH(:, 1:nm_current) = NNDH(:,1:nm_current)
      DH(:, 1:nm_current) = 0.0
      KH = nm_current

! Second: scale the mass
      VD = TM/H(4, 1)
      print *, 'Scaling mass by factor', VD
      DO WHILE (DABS(1.0D0 - VD) > 0.1)
         VD = MAX(0.9,MIN(VD, 1.1))
         H(4, 1:KH) = VD*H(4, 1:KH)       ! Scale mass

         KTH = 1
         JHOLD = 4
         EQUILIBRIUM = equilibrate_model(KT5)
         IF (.NOT. EQUILIBRIUM) THEN !H(:,1:nm_current) = NH(:,1:nm_current)
            print *, '*** Failed ***'
            stop
         end if

         VD = TM/H(4, 1)
      END DO
      H(4, 1:KH) = VD*H(4, 1:KH)          ! Scale mass

! Third: scale the surface rotational period
! Make the whole star rotate with the surface rate, if desired
      IF (START_WITH_RIGID_ROTATION) H(13,2:KH) = H(13,1)
! Now scale the surface rate, similar to the way the mass is scaled
! If we forced the star to rigid rotation before then this will set the
! rotation profile throughout the entire star
      VD = P1/H(13, 1)
      Pcrit = 2.0*CPI/sqrt(CG*TM/exp(3*H(7, 1)))/CSDAY
      print *, 'Scaling rotational period by factor', VD
! For rotation rates within 1/4 of critical be a bit more careful. Here we
! need to approach the desired rate smoothly, adjusting the structure of
! the star at each step.
! FIXME: This should be more akin to the bit that updates the code on the
! new mesh, ie, not necessarily bring the star into equilibrium.
      IF (P1/Pcrit < 4.0) THEN
         print *, P1, 'close to critical rate', Pcrit
         VD = 4.0*Pcrit/H(13,1)
         H(13, 1:KH) = VD*H(13, 1:KH)     ! Scale rotational period
         EQUILIBRIUM = equilibrate_model(KT5)
         DO WHILE (EQUILIBRIUM .AND. DABS(P1 - H(13,1)) > 1.0d-6)
            VD = MAX(0.9,MIN(VD, 1.1))
            H(13, 1:KH) = VD*H(13, 1:KH)
            KTH = 1
            JHOLD = 4
            EQUILIBRIUM = equilibrate_model(KT5)
            VD = P1/H(13, 1)
            IF (.NOT. EQUILIBRIUM) THEN
               print *, '*** Failed ***'
               stop
            end if
         END DO
      END IF
      H(13, 1:KH) = VD*H(13, 1:KH)        ! Scale rotational period

! Compute moment of inertia and surface potential
      IF (JF /= 2) THEN
         Q2 = H(6,1)
         H(6, 1:KH) = 1.0D0
         H(12, 1:KH) = 1.0D0              ! Moment of Inertia
         H(14, 1:KH) = 0.0D0              ! Gravitational potential
         DO IK = 2, KH
            VAR(:) = H(:, IK)
            CALL FUNCS1 ( IK, -2 )
            H(12, IK) = H(12, IK - 1) + R*R*M3/(DABS(QMU))
            H(14, IK) = H(14, IK - 1) + PHIMU/DABS(QMU)
         END DO
         H(6, 1:KH) = Q2
         H(15, 1:KH) = - GMR              ! Potential at the stellar surface
         SI = H(12, KH)                   ! Total moment of inertia

         DO IK = 1, KH
            H(12, IK) = (SI - H(12, IK))*DABS(Q2)  ! Moment of inertia of interior
            H(14, IK) = - GMR - H(14, IK)*DABS(Q2) ! Gravitational potential
         END DO
      END IF
! Total angular momentum integration
      H(NTAM, KH) = 0.0D0
      DO IK = 1, KH-1
         H(NTAM, IK) = H(NTAM, IK+1) + (H(12, IK) - H(12, IK+1))/H(13, IK)
      END DO

      IF (RELAX_LOADED_MODEL) THEN
         KTH = 1
         JHOLD = 4
         print *, 'Equilibrating...'
         EQUILIBRIUM = equilibrate_model(KT5)
         IF (.NOT. EQUILIBRIUM) THEN
            print *, '*** Failed ***'
            stop
         end if
         print *, 'Done'
         DH(:, 1:KH) = 0.0
      END IF
      
! Restore old init.dat
      call pop_init_dat(initdat, ik, KR1, KR2, KSV, KT5, ik)
      DT = DTB
      AGE = AGEB
      end subroutine

      function equilibrate_model(KT5)
      use mesh
      use constants
      implicit none
      logical :: equilibrate_model
      integer, intent(in) :: KT5
      double precision :: DTY
      integer :: I, JO, IDUMMY
      double precision :: ddummy
! Unnamed COMMON block
      DOUBLE PRECISION :: H(NVAR,NM),DH(NVAR,NM),EPS,DEL,DH0
      INTEGER :: KH,KTW,ID(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, ID
! COMMON block TVBLES
      DOUBLE PRECISION :: DT, ZQ(1), HSPN(2), RLF(2), ZET(2),
     $     XIT(2), AGE,BBM, MC(2), OM, PER, SM, ENC, TC(2), TFR, T0, M0,
     $     MTA, OM0, OMTA,A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS,
     $     SE, TN(2), WMH,WMHE, MH, MHE, MCO, VMG, BE(2), LH, LHE, LC,
     $     LNU, LTH, MCB(8),MSB(6), RCB(8), TCT(8), QR(81), PPR(81)
      INTEGER :: JHOLD, JM2, JM1
      COMMON /TVBLES/ DT, ZQ, HSPN, RLF, ZET, XIT, AGE, 
     : BBM, MC, OM, PER, SM, ENC, TC, TFR, T0, M0, MTA, OM0, OMTA, 
     : A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS, SE, TN, WMH, 
     : WMHE, MH, MHE, MCO, VMG, BE, LH, LHE, LC, LNU, LTH, MCB, 
     : MSB, RCB, TCT, QR, PPR, JHOLD, JM2, JM1

      JO = 0
      DTY = DT/CSY
      DO I=1, 40
         CALL SOLVER(20, ID, KT5, JO)
         IF (JO /= 0) EXIT
         AGE = 0.0D0
         CALL PRINTB ( DTY, JO, idummy, ddummy, 1, 22 )
         H(:,1:KH) = H(:,1:KH) + DH(:,1:KH)
         CALL NEXTDT ( DTY, JO, 22 )
         IF (DABS(LTH) < 1.0D-8 .or. LTH < 0.0D0) EXIT
      END DO
      equilibrate_model = .false.
      IF (LTH < 1.0D-6 .and. JO == 0) equilibrate_model = .true.

      end function
      
      end module

