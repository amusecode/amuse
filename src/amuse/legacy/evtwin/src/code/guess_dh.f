
!     Estimate what DH should be
      IF (JO /= 0) THEN
         KT6 = 3*ITER+KT5;
         JO = 0
         IF (.true. .or. VERBOSITY>0) write(1, *) 'Trying to guess DH'
!        Restore H from before the last unsuccesful call to SOLVER
         H(1:NVAR, 1:KH) = FH(1:NVAR, 1:KH)
         DH(1:NVAR, 1:KH) = 0.0d0
         !CALL GUESS_DH
         CALL GUESS_DH2
         CALL SOLVER ( ITER, IG, 0, JO )
         IF (JO == 0) THEN
            IF (.true. .or. VERBOSITY>0) write(1, *) 'Succes! :)'
         ELSE
            IF (.true. .or. VERBOSITY>0) write(1, *) 'Failed :('
         ENDIF
      END IF





! GUESS_DH
!  From the current values in H(:,:), estimate what the changes DH(:,:) are
!  likely to be and fill these in.
!  This is useful when (re)starting a model
      SUBROUTINE GUESS_DH
      USE MESH
      USE SETTINGS
      USE NUCLEOSYNTHESIS
      IMPLICIT NONE
!     COMMON blocks
!     /      /
      DOUBLE PRECISION :: H, DH, EP
      INTEGER :: KH, KTW, ID
      COMMON H(NVAR,NM), DH(NVAR,NM), EP(3), KH, KTW, ID(260)
!     /VBLS/
      DOUBLE PRECISION :: LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,WL,
     &     WCV, HP, WT, PHIM, GMR, SEP, M3, PX(NPX), SX(NPX,NM+1),QA(NM)
      COMMON /VBLES / LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,
     & WL, WCV, HP, WT, PHIM, GMR, SEP, M3, PX, SX, QA
!     /TVBLES/
      DOUBLE PRECISION :: DT, ZQ(1), HSPN(2), RLF(2), ZET(2),
     &     XIT(2), AGE,BM, MC(2), OM, PER, SM, ENC, TC(2), TFR, T0, M0,
     &     MTA, OM0, OMTA,A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS,
     &     SE, TN(2), WMH,WMHE, MH, MHE, MCO, VMG, BE(2),  LH, LHE, LC,
     &     LNU, LTH, MCB(8),MSB(6), RCB(8), TCT(8), QR(81), PPR(81)
      INTEGER :: JHOLD, JM2, JM1 
      COMMON /TVBLES/ DT, ZQ, HSPN, RLF, ZET, XIT, AGE, 
     & BM, MC, OM, PER, SM, ENC, TC, TFR, T0, M0, MTA, OM0, OMTA, 
     & A0, ATA, E0, ETA, CDD, BP, HORB, RO, RA2, RS, SE, TN, WMH, 
     & WMHE, MH, MHE, MCO, VMG, BE, LH, LHE, LC, LNU, LTH, MCB, 
     & MSB, RCB, TCT, QR, PPR, JHOLD, JM2, JM1
!     /STAT2 /
      DOUBLE PRECISION :: PL, RL, U, P, RHO, FK, T, SF, ST, ZT, GRADA,
     &     SCP,RF, RT, XHI, S, PR, PG, PF, PT, EN, RPP, R33, R34, RBE,
     &     RBP,RPC, RPNA, RPO, R3A, RAC, RAN, RAO, RANE, RCCA, RCO ,ROO,
     &     RGNE, RGMG, RCCG, RPNG, EX, ENX, WMU, DELTA, PHI, EXT, FKT,
     &     FKR, PRANDTL
      COMMON /STAT2 / PL, RL, U, P, RHO, FK, T, SF, ST, ZT, GRADA, SCP, 
     & RF, RT, XHI, S, PR, PG, PF, PT, EN, RPP, R33, R34, RBE, RBP, 
     & RPC, RPNA, RPO, R3A, RAC, RAN, RAO, RANE, RCCA, RCO,
     & ROO, RGNE, RGMG, RCCG, RPNG, EX, ENX, WMU, DELTA, PHI, EXT, FKT,
     & FKR, PRANDTL
!     /INF /
      DOUBLE PRECISION :: Q(NVAR), QD(NVAR), INFPADDING(NFUNC+NVAR*NFUNC)
      COMMON /INF   / Q, QD, INFPADDING
!     /NAMOUT/   ! SIZE(NAMEOUT) = NFUNC
      DOUBLE PRECISION ::  BCP, BCT, VP, VPK, VR, VRK, VT, VTK, VL, LK,
     &     LQ, MT, VM, VMK, SG, WT1, X1, X1T, X16, X16T, X4, X4T, X12,
     &     X12T, X20, X20T, BCM, VI, VIK, VPHI, PHIK, BCF, BCS, BCPH,
     &     X14, X14T, AVMU, SGTH, OMEGA, OMEGAT, SGAM, SI, BCA, BCE, XIM
     &     , XIK, DLRK, BCMB, X24, X24T, MENC, MENCT, MEA, MEAT, MET,
     &     METT, FILLUP_COMMONBLOCK_NAMOUT(2 + NSFSTAR + NXFUNC  )
      COMMON /NAMOUT/ BCP, BCT, VP, VPK, VR, VRK, VT,
     & VTK, VL, LK, LQ, MT, VM, VMK, SG, WT1, X1, X1T, X16, X16T, X4, 
     & X4T, X12, X12T, X20, X20T, BCM, VI, VIK, VPHI, PHIK, BCF, BCS, 
     & BCPH, X14, X14T, AVMU, SGTH, OMEGA, OMEGAT, SGAM, SI,
     & BCA, BCE, XIM, XIK, DLRK, BCMB, X24, X24T, MENC, MENCT, MEA, 
     & MEAT, MET, METT, FILLUP_COMMONBLOCK_NAMOUT
!     Local variables
      DOUBLE PRECISION :: DLR(KH), DMR(KH)
      DOUBLE PRECISION :: DLL(KH), DML(KH)
      DOUBLE PRECISION :: LOM(KH)
      DOUBLE PRECISION :: M(2,2), invM(2,2), detM, x(2), B(2)
      DOUBLE PRECISION :: GRADF, GRADT, DELT, DELF, VPP, APK, MK
      INTEGER :: IK, I

!     Clear DH array
      DH(:,1:KH) = 0.0d0

!     Determine the actual luminosity gradient by finite differencing
      !DL(2:KH-1) = H(8, 1:KH-2) - H(8, 3:KH)
      !DL(1) = H(8,1) - H(8, 2)
      !DL(KH) = H(8,KH-1) - H(8, KH)
      !DM(2:KH-1) = H(4, 1:KH-2) - H(4, 3:KH)
      !DM(1) = H(4,1) - H(4, 2)
      !DM(KH) = H(4,KH-1) - H(4, KH)
      DLR(1:KH-1) = H(8, 1:KH-1) - H(8, 2:KH)
      DLR(KH) = H(8,KH) - H(8, KH-1)
      DLL(2:KH) = H(8, 2:KH) - H(8, 1:KH-1)
      DLL(1) = H(8,1) - H(8, 2)

      DMR(1:KH-1) = H(4, 1:KH-1) - H(4, 2:KH)
      DMR(KH) = H(4,KH) - H(4, KH-1)
      DML(2:KH) = H(4, 2:KH) - H(4, 1:KH-1)
      DML(1) = H(4,1) - H(4, 2)

      LOM(1:KH) = 0.5d0*(DLL(1:KH)/DML(1:KH) + DLR(1:KH)/DMR(1:KH))

!     Estimate changes in variables
      DO IK=1, KH
         Q(17:24) = H(17:24, IK)
         QD(:) = 0.0d0

!        Numerical derivatives of GRAD
         Q(:) = H(:, IK)
         Q(1) = Q(1) + 1.0d-7*Q(1)
         CALL FUNCS1 ( IK, -1 )
         GRADF = GRAD
         DELF = Q(1) - H(1, IK)

         Q(:) = H(:, IK)
         Q(2) = Q(2) + 1.0d-7*Q(2)
         CALL FUNCS1 ( IK, -1 )
         GRADT = GRAD
         DELT = Q(2) - H(2, IK)

!        Actual un-perturbed function values
         Q(:) = H(:, IK)
         CALL FUNCS1 ( IK, -1 )

         VPP = CT(4) + CT(5)*P/(P + CT(9))
         APK = VPK/VPP
         GRADT = (GRADT-GRAD)/DELT
         GRADF = (GRADF-GRAD)/DELF
         MK = HT(1, 4, IK)

!        Setup coefficients for the matrix equation
         M(1,1) = 1 - GRAD*PT - GRADT*APK
         M(1,2) = -GRAD*PF - GRADF*APK
         M(2,1) = ST
         M(2,2) = SF
         B(1) = 0.0d0
         B(2) = DT * (EX+EN+ENX - LOM(IK)) / T

!        Invert the matrix
         detM = M(1,1)*M(2,2) - M(1,2)*M(2,1)
         invM(1,1) = M(2,2); invM(1,2) = -M(1,2)
         invM(2,2) = M(1,1); invM(2,1) = -M(2,1)
         invM(:,:) = invM(:,:) / detM

!        Compute estimates for dlogT, dlogf
         x(:) = matmul(invM(:,:), B(:))
         DH(1, IK) = X(2) * 1.0d-7
         DH(2, IK) = X(1) * 1.0d-7

         !print '(I6, 6E16.8)', IK, x
         !print '(I6, 2E16.8)', -IK, DH(2,IK), DH(1, IK)
         !print '(I6, 6E16.8)', -IK, B(:), GRAD, GRADA
         !print '(I6, 6E16.8)', -IK, GRADT, GRADF, APK
         !print '(I6, 6E16.8)', -IK, LOM(IK), EX+EN+ENX
         !print '(I6, 6E16.8)', -IK, M(1,1:2), M(2,1:2), detM
         !print '(I6, 6E16.8)', -IK, invM(1,1:2), invM(2,1:2)
         !print *, ''
      END DO
      END SUBROUTINE
      
! Try guessing changes DH based on the values in the current timestep
      SUBROUTINE GUESS_DH2
      USE MESH
      IMPLICIT NONE
      INTEGER :: K
      INTEGER :: KH, KTW, ID(130), IE(130)
      INTEGER :: IHOLD, JM2(2)
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      DOUBLE PRECISION :: DT, ZQ(81), PREV(81), PPR(81)
      DOUBLE PRECISION :: XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, XW(14)
      DOUBLE PRECISION :: AT, AF, AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA, 
     & SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, RPP, R33, R34, RBE, RBP,  
     & RPC, RPNA, RPO, R3A, RAC, RAN, RAO, RANE, RCCA, RCO, ROO, RGNE, 
     & RGMG, RCCG, RPNG, EX, ENX, WMU, DELTA, PHI, EXT, FKT, FKR, PRANDTL
      DOUBLE PRECISION :: DM(NM), DL(NM), L(NM), ETH
      DOUBLE PRECISION :: F1, F2
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, ID, IE
      COMMON /TVBLES/ DT, ZQ, PREV, PPR, IHOLD, JM2
      COMMON /ABUND / XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, XW
      COMMON /STAT2 / AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA, SCP, 
     & RF, RT, XHI, S, PR, PG, PF, PT, EN, RPP, R33, R34, RBE, RBP, 
     & RPC, RPNA, RPO, R3A, RAC, RAN, RAO, RANE, RCCA, RCO,
     & ROO, RGNE, RGMG, RCCG, RPNG, EX, ENX, WMU, DELTA, PHI, EXT, FKT,
     & FKR, PRANDTL

! First: reset DH to 0
      DH(1:NVAR, 1:KH) = 0.0D0

! Now make eucated guesses about what the change in each variable should look
!  like approximately. In a way this is assigning the changes according to an
!  explicit integration scheme and then improving on it in the implicit scheme.
! 1) Guess DAT and DAF (change in ln T and ln f) from the energy generation
!    rate: knowing dL/dm and knowing the nuclear energy generation rate, we
!    can calculate what dS/dt is. We need one extra assumption to convert this
!    to a guess for dlnT/dt and dlnf/dt; here use dP/dt=0 for simplicity; in
!    reality, this should be the equation of hydrostatic equilibrium
      
! Compute dm and dL, second order
      DM(1) = H(4,1) - H(4, 2)
      DL(1) = H(8,1) - H(8, 2)
      DO K=2, KH-1
         DM(K) = H(4, K-1) - H(4, K+1)
         DL(K) = H(8, K-1) - H(8, K+1)
      END DO
      DM(KH) = H(4,KH) - H(4, KH-1)
      DL(KH) = H(8,KH) - H(8, KH-1)

      DO K=1, KH
         XH = H(5, K)
         XHE = H(9, K)
         XC = H(10, K)
         XN = H(16, K)
         XO = H(3, K)
         XNE = H(11, K)
         AT = H(2, K)
         AF = H(1, K)
         CALL STATEF ( AF, AT )
         CALL NUCRAT ( AT )
         
         ETH = DL(K)/DM(K) - EX - EN - ENX

         IF (K == 1) THEN
            L(K) = 0
         ELSE
            L(K) = L(K-1) + (EX + EN + ENX)*DM(K)
         END IF
         
         !print *, K, ETH
         F1 = ST*PF - SF*PT
         F2 = ( RT*PF - RF*PT ) / ( RT*SF - RF*ST )

         DH(2, K) = -ETH/T  * ( PF - SF * F2)/F1 * DT * 1.0D-4
         DH(1, K) = ETH/T  * ( PT + ST * F2)/F1 * DT * 1.0D-5
         DH(7, K) = -SIGN(1.0D-1*H(7,K), -ETH)
         
         DH(8, K) = ( L(K) - H(8, K) ) * 1.0D-3

         DH(5, K) = -2.0*(RPP + RPC + RPNG + RPNA + RPO)*DT
         DH(9, K) = -4.0*(-(0.5*RPP + RPNA + RPO)
     &     + (3.0*R3A + RAC + 1.5*RAN + RAO + RANE)
     &     - (RCCA + RCO + 2.0*ROO + RGNE + RGMG))*DT
!         DH(10, K) = -12.0*((RPC - RPN) - (R3A - RAC)
!     :            + (2.0*(RCC + RCCG) + RCO))*DT
!         DH(16, K) = -14.0*(RPN + RPNG - RPC - RPO + RAN)*DT
!         DH(3, K) = -16.0*((RPO - RPNG) - (RAC - RAO)
!     :            + (RCO + 2.0*ROO - RGNE))*DT
!         DH(11, K) = -20.0*((RANE - RAN - RAO) + (RGNE - RGMG - RCC))*DT
         !IF (H(2, K) > 16.1D0) THEN
         !   DH(2, K) = (16.1D0 - H(2, K))*2.
         !END IF
         !print *, K, DH(2, K)
      END DO

      END SUBROUTINE

