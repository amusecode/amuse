      MODULE SOLVER_GLOBAL_VARIABLES
      USE MESH
      IMPLICIT NONE

      INTEGER :: BLOCK_VARS = 0
      INTEGER :: BLOCK_NMESH = 0

      DOUBLE PRECISION, ALLOCATABLE :: S(:, :)
      DOUBLE PRECISION, ALLOCATABLE :: C(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: XD(:)
      DOUBLE PRECISION, ALLOCATABLE :: ER(:)       ! Typical value of variables
      DOUBLE PRECISION, ALLOCATABLE :: EQN_SCALE(:)      ! Typical value of eqns
      DOUBLE PRECISION, ALLOCATABLE :: SBC_SCALE(:)    ! Typical value of sbcs
      DOUBLE PRECISION, ALLOCATABLE :: DDH(:,:)
!      DOUBLE PRECISION, ALLOCATABLE :: GRADF(:)

      ! Quantities needed for "linesearch" algorithm, see equations (9.7.7)
      ! through (9.7.9) in Numerical Recipes (2nd edition)
      DOUBLE PRECISION :: RESIDUE         ! Called g in NR
      DOUBLE PRECISION :: PREV_RESIDUE    ! Value of RESIDUE on previous iter
      DOUBLE PRECISION :: RESIDUE_DERIV   ! Called g' in NR
      
      INTEGER :: KEQ, KVB, KVC
      INTEGER :: KQ, KEE
      INTEGER :: KJ2, KJ5, KJ6, KJ10, KJ12, KI4
      INTEGER :: KE1, KE2, KE3, KE4, KBC, KEV, KFN, KL, JH1, JH2, JH3
      INTEGER :: KI1, KI2, KI3, KJ1, KJ11, KJ3, KJ4, KJ7, KJ8, KJ9
      INTEGER :: IK_FIRST, IK_LAST
      INTEGER :: KD(3*40)     ! 3 sets of 40 integers, for 40 variables

      END MODULE

      SUBROUTINE SOLVER ( ITER, IG, KT5, JO )
      USE MESH
      USE SETTINGS
      USE SOLVER_GLOBAL_VARIABLES
c     Explicitly type all variables (Steve, 5/08):
      IMPLICIT NONE

c     Arguments:
      INTEGER :: iter, ig, kt5, jo

c     External functions
      DOUBLE PRECISION, EXTERNAL :: LINESEARCH

c     Local variables:
      INTEGER :: i, ii, ij, ik, ikm, jk, k, kk
      DOUBLE PRECISION :: err, errpr, ert, es, fac, frr, vx, vz

      DIMENSION IKM(NVAR), ERT(NVAR), ES(NEQ), IG(130)
      LOGICAL :: CONVERGED

c     Common variables:
      INTEGER :: kh, ktw, kw
      INTEGER :: jmod, jb, jnn, jter, joc, jkh

      DOUBLE PRECISION :: h, dh, eps, del, dh0
      DOUBLE PRECISION :: ml, ql, xl, uc
      DOUBLE PRECISION :: errs, errppr, besterr, besth, bestdh

      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, KW(260)
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /SOLMON/ ERRS, ERRPPR, BESTERR, BESTH(NVAR,NM), BESTDH(NVAR,NM)

      CHARACTER*7 VARNAMES(1:NVAR)

      DATA VARNAMES/'   ln f','   ln T','    X16','      M','     X1',
     &              '      C','   ln r','      L','     X4','    X12',
     &              '    X20','      I','   Prot','    phi','  phi_s',
     &              '    X14','  H_orb','      e','     xi','    M_B',
     &              ' *1:X24','  22   ','  23   ','  24   ','*2:ln f',
     &              '*2:ln T','*2: X16','*2:   M','*2:  X1','*2:   C',
     &              '*2:ln r','*2:   L','*2:  X4','*2: X12','*2: X20',
     &              '*2:   I','*2:Prot','*2: phi','*2:phis','*2: X14',
     &              '  41   ','  42   ','  43   ','  44   ','  45   ',
     &              '  46   ','  47   ','  48   ','  49   ','  50   ',
     &              '  51   ','  52   ','  53   ','  54   ','  55   ',
     &              '  56   ','  57   ','  58   ','  59   ','  60   ',
     &              '  61   ','  62   ','  63   ','  64   ','  65   ',
     &              '  66   ','  67   ','  68   ','  69   ','  70   ',
     &              '  71   ','  72   ','  73   ','  74   ','  75   ',
     &              '  76   ','  77   ','  78   ','  79   ','  80   ',
     &              '  81   ','  82   ','  83   ','  84   ','  85   ',
     &              '  86   ','  87   ','  88   ','  89   ','  90   ',
     &              '  91   ','  92   ','  93   ','  94   ','  95   ',
     &              '  96   ','  97   ','  98   ','  99   ',' 100   ',
     &              ' 101   ',' 102   ',' 103   ',' 104   ',' 105   ',
     &              ' 106   ',' 107   ',' 108   ',' 109   ',' 110   ',
     &              ' 111   ',' 112   ',' 113   ',' 114   ',' 115   ',
     &              ' 116   ',' 117   ',' 118   ',' 119   ',' 120   '/

c-----------------------------------------------------------------------
c     Explicitly initialize local variables (Steve, 5/08).

      IKM(:) = 0
      ES(:) = 0
      ERRPR = 0
      VX = 0
      VZ = 0
c-----------------------------------------------------------------------

      ! Set number of equations/variables.
      ! Used to be a block copy, but does not work reliably when using
      ! aggressive optimisations in the compiler.
      KE1 = IG(1)
      KE2 = IG(2)
      KE3 = IG(3)                   ! Third order equations, always 0
      KBC = IG(4)
      KEV = IG(5)
      KFN = IG(6)
      KL  = IG(7)
      JH1 = IG(8)
      JH2 = IG(9)
      JH3 = IG(10)
      KD(1:120) = IG(11:130)        ! Permutations of variables/equations
      KQ = 1 - 2*KL
      KEQ = KE1 + KE2               ! Total number of equations in block matrix
      IF ( KEQ == 0 ) RETURN        ! No equns to solve
      KVB = KEQ + KEV               ! Total number of variables (incl. EVs)
      KVC = KVB + KVB               ! ?
      KE4 = KE3 + KE2               ! Total number of equations
C LOCATE VERTICAL (KJ) AND HORIZONTAL (KI) PARTITIONS IN MATRIX
      KJ1 = KBC + 1
      KJ2 = KEQ + 1
      KJ3 = KEQ + KBC
      KJ4 = KJ3 + 1
      KJ5 = KEQ + KEQ
      KJ6 = KJ5 + 1
      KJ7 = KJ5 + KBC
      KJ8 = KJ7 + 1
      KJ9 = KJ5 + KEQ
      KJ10 = KJ9 + 1
      KJ11 = KJ9 + KEV
      KJ12 = KJ10 + KEV
      KI1 = KEQ - KBC
      KI2 = KI1 + 1
      KI3 = KI1 + KEV
      KI4 = KI3 + 1

!     If the number of variables/meshpoints is larger than the size of the
!     arrays (eg, we loaded a new model), we need to reallocate them.
      IF ( ALLOCATED(S) .AND. (KH>BLOCK_NMESH .OR. KVB>BLOCK_VARS) ) THEN
         DEALLOCATE( S )
         DEALLOCATE( C )
         DEALLOCATE( XD )
         DEALLOCATE( DDH )
         DEALLOCATE( ER )
         DEALLOCATE( EQN_SCALE )
         DEALLOCATE( SBC_SCALE )
!         DEALLOCATE( GRADF )
      END IF
!     Allocate memory for storage of the Jacobian
      IF (.NOT. ALLOCATED(S)) THEN
         BLOCK_VARS = KVB
         BLOCK_NMESH = KH
         ALLOCATE( S(3*KVB+1, KVB) )
         ALLOCATE( C(KH+1, KVB, KVB+1) )
         ALLOCATE( XD(KH) )
         ALLOCATE( DDH(KH, KVB) )
         ALLOCATE( ER(NVAR) )
         ALLOCATE( EQN_SCALE(NVAR) )
         ALLOCATE( SBC_SCALE(NVAR) )
!         ALLOCATE( GRADF(KEQ*KH + KEV) )
      END IF
C FIRST MESHPOINT IS IK_FIRST, LAST IS IK_LAST
      IK_FIRST = 1 + KL*(KH - 1)
      IK_LAST = KH + 1 - IK_FIRST
C DETERMINE 'TYPICAL', USUALLY MAXIMAL, VALUES FOR EACH VARIABLE
      DO IJ = 1, KVB
         I = KD(IJ)
         IF ( JNN > 1 ) ES(I) = ER(I) 
         ER(I) = 1.0D0
         DO K = 1, KH
            ER(I) = DMAX1(ER(I), DABS(H(I, K)))
         END DO
         IF ( ER(I) == 0.0D0 ) ER(I) = 1.0D0
         IF ( JNN > 1 .AND. ES(I) /= 1.0D0 ) ER(I) = 0.5D0*(ER(I) + ES(I))
      END DO

!     Now determine typical values for each equation, based on current
!     values in H(,) and the current timestep.
!     Unfortunately we need to go through FUNCS1() to do this (perhaps this
!     could better be collected in printb and then passed in here).
      CALL SCALE_EQN ( EQN_SCALE, SBC_SCALE )

!     KEE is the left-most KEQxKVB sub-block we need to solve. If there are
!     no higher order equations to solve, the left-most block (1, derivatives
!     wrt variables on the previous meshpoint) is always filled with 0s, so
!     we skip it in that case.
      KEE = 1
      IF (KE4 == 0) KEE = 2
C BEGIN ITERATIVE LOOP
      FAC = 1.0D0
      BESTERR = 1.0D0
      CONVERGED = .FALSE.
      C(:,:,:) = 0.0D0
      S(:, :) = 0.0D0
      RESIDUE = 0.0D0
      PREV_RESIDUE = 0.0D0
      DO JTER = 1, ITER
      ERRPPR = ERRPR
      ERRPR = ERR
      ERR = 0.0D0
      JKH = 0
      IF ( JNN == JH1 .AND. JTER == JH2 ) JKH = 1
C Solve system of linear equations by Gaussian elimination
      CALL INVERT_MATRIX( JO )
      IF (JO == 1) EXIT
C ESTIMATE ACCURACY OF ITERATION
      DO IJ = 1, KVB
         VX = 0.0D0
         DO IK = 1, KH
            DDH(IK, IJ) = MIN(CLIMIT, C(IK, IJ, 1))*ER(KD(IJ))
            VZ = DABS(C(IK, IJ, 1))
            IF ( VZ >= VX ) THEN
               VX = VZ
               IKM(IJ) = IK
            END IF
            ERR = ERR + VZ
         END DO
         ERT(IJ) = VX
      END DO
      ERR = ERR/(KVB*KH)

C Alternative to using a line search algorithm: decrease FAC, more or less
C ad-hoc
      IF (JTER > 1 .AND. ERR > ERRPR) THEN
         FAC = 0.8D0*FAC;
      ELSE
         FAC = MIN(DEL/DMAX1(ERR, DEL), 1.1D0*FAC)
      ENDIF
      !IF (JTER > 1) FAC = LINESEARCH( JO )
      !FAC = 1.0d0
      !IF (JTER > 1 .AND. RESIDUE > PREV_RESIDUE) FAC = LINESEARCH( JO )
      FRR = DLOG10(ERR)
      ERT(1:KVB) = DLOG10(ERT(1:KVB) + 1.3D-10)

C Write convergence information to summary file
      IF ( JTER == KT5+1 )
     & WRITE (JB, '(A4,2X,A3,1X,A4,2X,A3,1X, 17(1X,A7),/,3(20X,17(1X,A7),/))') 
     &        'Iter', 'Err', 'Ferr', 'Fac', (VARNAMES(KD(IJ)), IJ=1,KVB)
      IF ( JTER > KT5 )
     : WRITE (JB, 992) JTER, FRR, LOG10(RESIDUE), FAC, (IKM(IJ), ERT(IJ), IJ = 1, KVB) 
  992 FORMAT (1X, I2, 3F6.2, 17(I4, F4.1),/, 3(21X, 17(I4, F4.1),/))
      CALL FLUSH ( JB )
C APPLY CORRECTIONS, SCALED DOWN IF TOO LARGE
      DO I = 1, KVB
         IJ = KD(I)
         DH(IJ, 1:KH) = DH(IJ, 1:KH) - FAC*DDH(1:KH, I)
      END DO
      CALL CHECKS
      IF ( ERR < BESTERR ) THEN
         BESTH(1:NVAR,1:KH) = H(1:NVAR,1:KH)
         BESTDH(1:NVAR,1:KH) = DH(1:NVAR,1:KH)
         BESTERR = ERR
      END IF
      ERRS = ERR
      IF ( ERR > 1.0D1*DEL ) JO = 1
      DO I = 1, KVB
         IF (ERT(I) > 0.0d0) JO = 1
      END DO
      IF (ERR < EPS) CONVERGED = .TRUE.
      !IF (ERR < EPS .OR. RESIDUE < EPS*EPS) CONVERGED = .TRUE.
      IF (RESIDUE < EPS*EPS) CONVERGED = .TRUE.
      IF ( CONVERGED .OR. JO == 1 ) THEN
         RETURN
      END IF
C CONTINUE ITERATING IF NOT YET ACCURATE ENOUGH
      END DO
      JO = 1
      RETURN
      END SUBROUTINE

      SUBROUTINE INVERT_MATRIX(JO)
      USE MESH
      USE SOLVER_GLOBAL_VARIABLES
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: JO

C Blank common block
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, KW(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW
C Local variables
      INTEGER :: JK, KK, JKLAST, II, IJ

      PREV_RESIDUE = RESIDUE
      RESIDUE = 0.0D0
      RESIDUE_DERIV = 0.0D0
!      GRADF(:) = 0.0D0
      JK = IK_FIRST + KL
C EVALUATE FUNCTIONS AT FIRST MESHPOINT
      XD(:) = 0.0D0
      CALL DIFFERENCES ( JK, KI2, KEQ, 80 + KEV ) 
      CALL DIVIDE_ROWS ( KI2, KJ6, KEQ, JK, JO ) 
      IF ( JO == 1 ) RETURN
      JK = JK + KQ
      IF ( IK_FIRST == IK_LAST ) RETURN
C DITTO SECOND, ELIMINATING SOME UNKNOWNS
      CALL DIFFERENCES ( JK, 1, KEQ, 40 )
      CALL ELIMINATE   ( 1, KJ2, KEQ, KJ3, KI2, KJ4, KJ5, JK - KQ,   1 ) 
      CALL DIVIDE_ROWS ( 1, KJ4, KEQ, JK, JO) 
      IF ( JO == 1 ) RETURN
      DO JK = IK_FIRST + KL + 2*KQ, IK_LAST + KL, KQ
C DITTO REMAINING MESHPOINTS EXCEPT LAST
         CALL DIFFERENCES ( JK, 1, KEQ, 40 )
         CALL ELIMINATE   ( 1,   1, KE4, KBC, KI2, KJ1, KEQ, JK - 2*KQ, 1 )
         CALL ELIMINATE   ( 1, KJ1, KE4, KEQ,   1, KJ4, KJ5, JK - KQ,   2 )
         CALL ELIMINATE   ( 1, KJ2, KEQ, KJ3, KI2, KJ4, KJ5, JK - KQ,   3 )
         CALL DIVIDE_ROWS ( 1, KJ4, KEQ, JK, JO)
         IF ( JO == 1 ) RETURN
      END DO
      JK = IK_LAST + KL + KQ
C DITTO LAST MESHPOINT
      CALL DIFFERENCES ( JK, 1, KI3, 80 ) 
      CALL ELIMINATE   ( 1, KJ2, KE4, KJ3, KI2, KJ4, KJ5, JK - 2*KQ, 1 )
      CALL ELIMINATE   ( 1, KJ4, KE4, KJ5,   1, KJ8, KJ9, JK - KQ,   2 )
      CALL ELIMINATE   ( 1, KJ6, KI3, KJ7, KI2, KJ8, KJ9, JK - KQ,   3 )
C SOLVE FOR CORRECTIONS AT LAST MESHPOINT
      CALL DIVIDE_ROWS ( 1, KJ8, KI3, JK, JO)
      IF ( JO == 1 ) RETURN
C BY BACK-SUBSTITUTION FIND CORRECTIONS THROUGHOUT
      II = 1
      JKLAST = IK_FIRST + KL
      IF (KBC == 0) JKLAST = JKLAST-KQ
      DO JK = IK_LAST + KL, JKLAST, -KQ
         IF ( JK == IK_FIRST + KL ) II = KI2
         KK = JK + KQ   
         DO IJ = 1, KI2-1
            C(JK, II:KEQ, KI4) = C(JK, II:KEQ, KI4) - C(JK, II:KEQ, IJ)*C(KK, IJ, KI4)
         END DO
         KK = IK_LAST + 1 - KL
         DO IJ = KI2, KI3
            C(JK, II:KEQ, KI4) = C(JK, II:KEQ, KI4) - C(JK, II:KEQ, IJ)*C(KK, IJ, KI4)
         END DO
         !CALL PRINTC ( JK, JKH, II, KEQ, 1, KI4 )
      END DO
      IF ( KEV /= 0 ) THEN
         DO JK = 1, KH
            C(JK, KJ2:KVB, 1) = C(IK_LAST + 1 - KL, KJ2 - KBC:KVB - KBC, KI4)
         END DO
      END IF
      IF ( KBC /= 0 ) THEN
         C(1:KH, 1:KBC, 1) = C(KL+1:KH+KL, 1 + KI1:KBC + KI1, KI4)
      END IF
      C(1:KH, 1 + KBC:KI1 + KBC, 1) = C(2 - KL:KH + 1 - KL, 1:KI1, KI4)
      !CALL PRINTC ( JK, JKH, 1, KVB, 1, 1 )
      END SUBROUTINE

C Line search algoritm, as described in Numerical Recipes, 2nd edition, ¤9.7
C Slightly different from the implemantation there.
C This function returns the value "lambda" (or "FAC") by which the current
C set of corrections are to be multiplied.
C FIXME: in order for this to work properly, equations need to be scaled
C with their typical values, otherwise checking whether the residuals are
C "large" doesn't make sense.
      FUNCTION LINESEARCH(JO)
      USE MESH
      USE SOLVER_GLOBAL_VARIABLES
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: JO
      DOUBLE PRECISION :: LINESEARCH
C Blank common block
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, KW(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW
C Local variables
      DOUBLE PRECISION, PARAMETER :: ALPHA = 1.0D-6
      DOUBLE PRECISION :: A, B, DISC, RHS1, RHS2, SLOPE
      DOUBLE PRECISION :: LAMBDA, LAMBDA_MIN, NEW_LAMBDA, PREV_LAMBDA
      DOUBLE PRECISION :: INIT_RESIDUE
      DOUBLE PRECISION :: DH_OLD(KVB, KH)
      INTEGER :: IK, IJ, I

C     Backup DH
      DO IJ = 1, KVB
         I = KD(IJ)
         DH_OLD(IJ, 1:KH) = DH(I, 1:KH)
      END DO

C     Determine downward slope, grad f.dx
      SLOPE = 0.0D0
C     Contribution from normal variables
!      DO IK=1, KH
!         DO IJ = 1, KEQ
!            SLOPE = SLOPE + GRADF((IK-1)*KEQ+IJ) * DDH(IK, IJ)
!         END DO
!      END DO
C     Contribution from eigenvalues
!      DO IJ = KEQ+1, KVB
!         SLOPE = SLOPE + GRADF((KH-1)*KEQ+IJ) * DDH(KH, IJ)
!      END DO

C     Find allowed range of values for "lambda"
      LAMBDA_MIN = 0.0D0
      DO IJ = 1, KVB
         I = KD(IJ)
         LAMBDA_MIN = MAX( LAMBDA_MIN, DABS(DDH(IK, IJ))/ER(I) )
      END DO
      LAMBDA_MIN = 1.0D-4 / LAMBDA_MIN
      LAMBDA = 1.0D0
      LINESEARCH = 1.0d0
      SLOPE = -2.0D0*RESIDUE

C     Store residual at original point
      INIT_RESIDUE = RESIDUE

      DO WHILE (.TRUE.)
C        New value of DH
         DO IJ = 1, KVB
            I = KD(IJ)
            DH(I, 1:KH) = DH_OLD(IJ, 1:KH) + LAMBDA*DDH(1:KH, IJ)
         END DO
C        Evaluate functions for new parameters
         !print *, 'Line search:', LAMBDA, RESIDUE
         CALL INVERT_MATRIX( JO )
         IF ( JO == 1 ) EXIT
C        Check for convergence on lambda (parameters)
         IF (LAMBDA < LAMBDA_MIN) EXIT
C        Check for convergence on function value
         IF (RESIDUE <= INIT_RESIDUE + ALPHA*LAMBDA*SLOPE) EXIT
C        Determine new value of lambda
         IF (LAMBDA == 1.0d0) THEN     ! FIXME: float comparison
            NEW_LAMBDA = -SLOPE/(2.0d0*(RESIDUE - INIT_RESIDUE - SLOPE))
         ELSE
            RHS1 = RESIDUE - INIT_RESIDUE - LAMBDA*SLOPE
            RHS2 = PREV_RESIDUE - INIT_RESIDUE - PREV_LAMBDA*SLOPE
            A = (RHS1/LAMBDA**2 - RHS2/PREV_LAMBDA**2) / (LAMBDA - PREV_LAMBDA)
            B = -(PREV_LAMBDA*RHS1/LAMBDA**2 - LAMBDA*RHS2/PREV_LAMBDA**2) /
     &            (LAMBDA - PREV_LAMBDA)
            IF (A == 0.0d0) THEN    ! This will NEVER happen in practice...
               NEW_LAMBDA = -SLOPE/(2.0d0*B)
            ELSE
               DISC = B**2 - 3.0d0*A*SLOPE
               IF (DISC < 0.0d0) THEN
                  NEW_LAMBDA = 0.5*LAMBDA
               ELSE IF (B<= 0.0D0) THEN
                  NEW_LAMBDA = (-B + DSQRT(DISC))/(3.0D0*A)
               ELSE
                  NEW_LAMBDA = -SLOPE/(B+DSQRT(DISC))
               END IF
               NEW_LAMBDA = MIN(0.5d0*LAMBDA, NEW_LAMBDA)
            END IF
         END IF
         PREV_LAMBDA = LAMBDA
         LAMBDA = MAX(NEW_LAMBDA, 0.1d0 * LAMBDA)
      END DO

C     Restore DH and return
      DO IJ = 1, KVB
         I = KD(IJ)
         DH(I, 1:KH) = DH_OLD(IJ, 1:KH)
      END DO
      LINESEARCH = LAMBDA
      END FUNCTION

      ! Compute FUNCS and derivatives
      SUBROUTINE COMPFUNCS ( K )
      USE MESH
      USE SOLVER_GLOBAL_VARIABLES
      IMPLICIT NONE
c     Parameter:
      INTEGER :: k

c     Local variables:
      INTEGER :: i, ji
      DOUBLE PRECISION :: dx, dvx, vx, sign_vx

c     Common variables:
      INTEGER :: kh, ktw, kw
      INTEGER :: jmod, jb, jnn, jter, joc, jkh

      DOUBLE PRECISION :: h, dh, eps, del, dh0
      DOUBLE PRECISION :: ml, ql, xl, uc
      DOUBLE PRECISION :: var, dvar, fn1, dfn1
      DOUBLE PRECISION :: fn2, dfn2, equ, dequ

      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, KW(260)
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /INF   / VAR(NVAR), DVAR(NVAR), FN1(NFUNC), DFN1(NVAR,NFUNC)
      COMMON /INE   / FN2(3,NFUNC), DFN2(3,NVAR,NFUNC), EQU(NEQ), DEQU(NVAR,3,NEQ)
      DOUBLE PRECISION :: FF(NFUNC)

c-----------------------------------------------------------------------
c     Explicitly initialize all local variables (Steve, 5/08).

      i = 0
      ji = 0
      dx = 0
      dvx = 0
      vx = 0
c-----------------------------------------------------------------------

      DVAR(1:NVAR) = DH(1:NVAR, K - KL)
       VAR(1:NVAR) =  H(1:NVAR, K - KL) + DVAR(1:NVAR)
      IF ( JOC == 1 ) THEN
         XD(1:KVC) = XD(1 + KVB:KVC + KVB)
C EVALUATE AND STORE THE REQUIRED FUNCS OF THE INDEP. VBLES
         CALL FUNCS1 ( K - KL, 0 )
         FF(:) = FN1(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The effect of the above code is the following: After the call to FUNCS1, the
!  array FN1() contains the values for independent variables and boundary cond.
! This is mapped onto the array FN2 for three subsequent meshpoints. The end
!  result is that the order of arrays in NAMEIN in EQUNS1 corresponds to the
!  order of variables in NAMEOUT in FUNCS1.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C VARYING INDEP. VBLES IN TURN, EVALUATE NUMERIC DERIVATIVES OF FUNCS
         DO I = 1, KVB
            JI = KD(I)
            VX = VAR(JI)
            DVX = DVAR(JI)
            SIGN_VX = 1
            IF ( DVX < 0.0D0 ) SIGN_VX = -1
            DX = SIGN_VX*DH0*MAX(ABS(VX), 1.0D0)  
            XD(I + KVC) = ER(JI)/DX
            VAR(JI) = VX + DX
            DVAR(JI) = DVX + DX
            CALL FUNCS1 ( K - KL, JI )
            DFN1(I, 1:NFUNC) = FN1(1:NFUNC)
            DVAR(JI) = DVX
            VAR(JI) = VX
         END DO
         FN1(:) = FF(:)
      ELSE
C ALTERNATIVE TO BE USED IF DERIVS COMPUTED ANALYTICALLY
         CALL FUNCS2 ( K - KL, 1, KH )
      END IF
      END SUBROUTINE COMPFUNCS

      ! JZ is the `start' index in the list of permutations
      SUBROUTINE DIFFERENCES ( K, JX, JY, JZ )
      USE MESH
      USE SOLVER_GLOBAL_VARIABLES
      IMPLICIT NONE

c     Parameters:
      INTEGER, INTENT(IN) :: k, jx, jy, jz

c     Local variables:
      INTEGER :: i, iee, ieq, ii, ivb, j, jeq, jj, jvb3, jvb4
      DOUBLE PRECISION :: ds

c     Common variables:
      INTEGER :: jmod, jb, jnn, jter, joc, jkh

      DOUBLE PRECISION :: ml, ql, xl, uc
      DOUBLE PRECISION :: var, dvar, fn1, dfn1
      DOUBLE PRECISION :: fn2, dfn2, equ, dequ

      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /INF   / VAR(NVAR), DVAR(NVAR), FN1(NFUNC), DFN1(NVAR,NFUNC)
      COMMON /INE   / FN2(3,NFUNC), DFN2(3,NVAR,NFUNC), EQU(NEQ), DEQU(NVAR,3,NEQ)

c-----------------------------------------------------------------------
c     Explicitly initialize all local variables (Steve, 5/08).

      i = 0
      iee = 0
      ieq = 0
      ii = 0
      ivb = 0
      j = 0
      jeq = 0
      jj = 0
      jvb3 = 0
      jvb4 = 0
      ds = 0
c-----------------------------------------------------------------------

      IF ( JKH >= 1 .AND. ( K < 4.OR.K > BLOCK_NMESH-2 )) THEN
         JKH = JH3 + 1
      ELSE
         JKH = 1
      END IF
C REDEFINE PREVIOUS CURRENT MESHPOINT AS CURRENT PREVIOUS MESHPOINT
      IF ( (KL == 0.AND.K <= BLOCK_NMESH) .OR. (KL == 1.AND.K >= 2) ) THEN
          FN2(KEE:2, 1:NFUNC)        =  FN2(KEE+1:3, 1:NFUNC)
         DFN2(KEE:2, 1:KVB, 1:NFUNC) = DFN2(KEE+1:3, 1:KVB, 1:NFUNC)

C Call function to perform the following tasks:
C * evaluate argument list of vbles and increments at current meshpoint
C * evaluate and store the required funcs of the indep. vbles
C * varying indep. vbles in turn, evaluate numeric derivatives of funcs
C * Or: Store result of derivs computed analytically
C These values can be recomputed (typically) or stored in a cache
         CALL COMPFUNCS( K )
         FN2(3, 1:NFUNC) = FN1(1:NFUNC)
         DFN2(3, 1:KVB, 1:NFUNC) = DFN1(1:KVB, 1:NFUNC)
      END IF
      IF ( JOC == 1 ) THEN
C EVALUATE AND STORE THE DIFFERENCE EQUNS WHICH ARE TO BE SATISFIED
         CALL EQUNS1 ( K, KL, KQ )
         DO IEQ = JX, JY
            JEQ = KD(IEQ + JZ)
            S(KJ12, IEQ) = EQU(JEQ)
            RESIDUE = RESIDUE + 0.5D0*(EQU(JEQ)/EQN_SCALE(JEQ))**2
         END DO
C VARYING INDEP. VBLES IN TURN, EVALUATE NUMERIC DERIVATIVES OF EQUNS
         JVB3 = KEE*KVB - KVB                   ! First column in this block
         S(KJ5+KEQ+1:KJ5+KVB, JX:JY) = 0.0D0    ! Clear derivatives to EVs
         DO IEE = KEE, 3
            ! Loop over all normal variables
            DO IVB = 1, KEQ
               JVB3 = JVB3 + 1
               JVB4 = KEQ*(IEE - 1) + IVB
               FN1(1:NFUNC) = FN2(IEE, 1:NFUNC)
               FN2(IEE, 1:NFUNC) = DFN2(IEE, IVB, 1:NFUNC)
               CALL EQUNS1 ( K, KL, KQ )
               FN2(IEE, 1:NFUNC) = FN1(1:NFUNC)
               DO IEQ = JX, JY
                  JEQ = KD(IEQ + JZ)
                  DS = (EQU(JEQ) - S(KJ12, IEQ))*XD(JVB3)
                  S(JVB4, IEQ) = DS
                  RESIDUE_DERIV = RESIDUE_DERIV - S(KJ12, IEQ)**2
                  ! Store EQU.J vector
!                  I = BLOCK_NMESH*KEQ - (K + (2 - IEE))*KEQ + IVB
!                  GRADF(I) = GRADF(I) + S(KJ12,IEQ)*DS
               END DO
            END DO
            ! Loop over all eigen values; these derivatives have their own
            ! column off to the right of the matrix.
            DO IVB = KEQ+1, KVB
               JVB3 = JVB3 + 1
               JVB4 = KJ5 + IVB                 ! Eigenvalues in column
               FN1(1:NFUNC) = FN2(IEE, 1:NFUNC)
               FN2(IEE, 1:NFUNC) = DFN2(IEE, IVB, 1:NFUNC)
               CALL EQUNS1 ( K, KL, KQ )
               FN2(IEE, 1:NFUNC) = FN1(1:NFUNC)
               DO IEQ = JX, JY
                  JEQ = KD(IEQ + JZ)
                  DS = (EQU(JEQ) - S(KJ12, IEQ))*XD(JVB3)
                  S(JVB4, IEQ) = S(JVB4, IEQ) + DS
                  RESIDUE_DERIV = RESIDUE_DERIV - S(KJ12, IEQ)**2/BLOCK_NMESH
                  ! Store EQU.J vector (eigenvalues go at the end)
!                  I = BLOCK_NMESH*KEQ + IVB - KEQ
!                  GRADF(I) = GRADF(I) + S(KJ12, IEQ)*DS
               END DO
            END DO
         END DO
      ELSE
C ALTERNATIVE TO BE USED IF DERIVS COMPUTED ANALYTICALLY
         CALL EQUNS2 ( K, 1, BLOCK_NMESH, KEQ )
         DO J = JX, JY
            JJ = KD(J + JZ)
            DO I = 1, KEQ
               II = KD(I)
               S(I, J) = DEQU(I, 1, JJ)*ER(II)
               S(I + KEQ, J) = DEQU(I, 2, JJ)*ER(II)
               S(I + KJ5, J) = DEQU(I, 3, JJ)*ER(II)
            END DO
            IF ( KEV /= 0 ) THEN
               DO I = KJ2, KVB
                  S(I + KJ5, J) = (DEQU(I, 1, JJ) + DEQU(I, 2, JJ) 
     :                         + DEQU(I, 3, JJ)) * ER(KD(I))
               END DO
            END IF
            S(KJ12, J) = EQU(JJ)
         END DO
      END IF
      RETURN
      END

      SUBROUTINE DIVIDE_ROWS ( ID1, JD1, ID2, K, JO )
      USE MESH
      USE SOLVER_GLOBAL_VARIABLES
      IMPLICIT NONE
c     Row flags
      INTEGER, PARAMETER :: FORCED_PIVOT = 1
      INTEGER, PARAMETER :: CURRENT_PIVOT = 2
      INTEGER, PARAMETER :: LATER_PIVOT = 3
      INTEGER, PARAMETER :: USED_PIVOT = 4
      INTEGER, PARAMETER :: POSSIBLE_PIVOT = 5
      INTEGER, PARAMETER :: SINGULAR = 6

c     Value to promote a "future" pivot to a current candidate pivot
      INTEGER, PARAMETER :: NEXT_PIVOT = 2

c     Parameters:
      INTEGER :: id1, jd1, id2, k, jo

c     Local variables:
      INTEGER :: i, ii, im, j, jc, jd2, jd3, jj, jl, jmi
      DOUBLE PRECISION :: vm, vs, vt, vx
      LOGICAL :: ALL_DONE

c     Common variables:
      INTEGER :: jmod, jb, jnn, jter, joc, jkh

      DOUBLE PRECISION :: ml, ql, xl, uc

      INTEGER :: ITYPE(NVAR), JMAX(NVAR+1)
      LOGICAL :: IDONE(NVAR)
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH

c-----------------------------------------------------------------------
c     Explicitly initialize all local variables (Steve, 5/08).

      i = 0
      ii = 0
      im = 0
      j = 0
      jc = 0
      jd2 = 0
      jd3 = 0
      jj = 0
      jl = 0
      jmi = 0
      vm = 0
      vs = 0
      vt = 0
      vx = 0
c-----------------------------------------------------------------------

      JC = 1
      JD2 = JD1 + ID2 - ID1
      IF ( ID1 > ID2 .OR. JD2 + 1 > KJ12 ) RETURN
      !IF ( JKH >= 2 ) CALL PRINTS ( K, JKH, ID1, ID2, 0 )
      IDONE(ID1:ID2) = .FALSE.
      ITYPE(ID1:ID2) = POSSIBLE_PIVOT
      ALL_DONE = .FALSE.
      DO WHILE (JC < 5000 .AND. .NOT. ALL_DONE)
      JC = JC + 1
      DO I = ID1, ID2
         IF ( ITYPE(I) < LATER_PIVOT ) ITYPE(I) = ITYPE(I) + NEXT_PIVOT
         IF ( ITYPE(I) < POSSIBLE_PIVOT ) IDONE(JMAX(I)) = .TRUE.
      END DO
      VT = 0.0
c Locate the most significant (remaining) row. `Significance' of a row is 
c the ratio of largest |element| = VM to sum of remaining |elements| = VS 
      DO I = ID1, ID2
         IF ( ITYPE(I) >= POSSIBLE_PIVOT ) THEN
            VM = 0.0
            DO J = JD1, JD2
               JJ = J - JD1 + ID1
               IF ( .NOT. IDONE(JJ) ) THEN
                  VX = DABS(S(J, I))
                  IF ( VX >= VM ) THEN
                     VM = VX
                     JL = JJ
                  END IF
               END IF
            END DO
            IF ( JL < 1 .OR. JL > BLOCK_VARS+1 ) GO TO 8
            JMAX(I) = JL
            VS = 0.0D0
            DO J = JD1, JD2
               IF ( J - JD1 + ID1 /= JL ) VS = VS + DABS(S(J, I))
            END DO
            IF ( VM ==  0.0D0 ) THEN
               ITYPE(I) = SINGULAR
               VX = 0.0D0
            END IF
            IF ( VS == 0.0D0 ) THEN    ! Only one non-zero element in row
               ITYPE(I) = FORCED_PIVOT
               IF ( VM > 0.0D0 ) THEN
                  IDONE(JL) = .TRUE.
                  VX = 2.0D0
               END IF
            ELSE
               VX = VM/(VM + VS)
            END IF
            IF ( VX >= VT ) THEN
               VT = VX
               IM = I
            END IF
         END IF
      END DO
      IF ( IM < 1 .OR. IM > NVAR ) GO TO 8
      IF ( ITYPE(IM) == POSSIBLE_PIVOT ) ITYPE(IM) = CURRENT_PIVOT
c Largest element moduluswise of most significant row is the leading 
c pivot; eliminate elements above and below it
      DO I = ID1, ID2
         IM = ITYPE(I)
         IF ( IM < LATER_PIVOT ) THEN
            JMI = JMAX(I) + JD1 - ID1
            JD3 = JD1
            IF ( IM == FORCED_PIVOT ) JD3 = JD2 + 1
            VX = 1.0D0/S(JMI, I)
            S(JD3:KJ12, I) = VX*S(JD3:KJ12, I)
            S(JMI, I) = 1.0D0
            DO II = ID1, ID2
               IF ( ITYPE(II) > LATER_PIVOT ) THEN
                  VX = S(JMI, II)
                  S(JD3:KJ12, II) = S(JD3:KJ12, II) - VX*S(JD3:KJ12, I)
                  S(JMI, II) = 0.0D0
               END IF
            END DO
         END IF
         IDONE(I) = .FALSE.
      END DO
c Are we done now?
      ALL_DONE = .TRUE.
      DO I = ID1, ID2
         IF ( ITYPE(I) == POSSIBLE_PIVOT .OR. ITYPE(I) <= CURRENT_PIVOT ) THEN
            ALL_DONE = .FALSE.
            EXIT
         END IF
      END DO
      END DO
      DO I = ID1, ID2
         DO J = JD2 + 1, KJ12
            C(K, JMAX(I), J - KJ12 + KI4) = S(J, I)
         END DO
      END DO
      IF ( JKH >= 2 ) GO TO 9
      IF ( JC < 5000 ) RETURN
   8  JO = 1
   9  CONTINUE
c Some emergency printout
!      WRITE (10, 991) K, JNN, JTER, JKH, JL, IM, JC, ITYPE, JMAX, IDONE
!  991 FORMAT (7I5, /, 40I3, /, 41I3, /, 40I3)
      !CALL PRINTC ( K, JKH, ID1, ID2, JD2 + 1 - KJ12 + KI4, KI4 )
      RETURN
      END

      SUBROUTINE ELIMINATE ( IS1, JS1, IS2, JS2, IS3, JS4, JS5, K, JCL )
      USE MESH
      USE SOLVER_GLOBAL_VARIABLES
      IMPLICIT NONE

c     Parameters:
      INTEGER :: is1, js1, is2, js2, is3, js4, js5, k, jcl

c     Local variables:
      INTEGER :: i, jj, jjj
      DOUBLE PRECISION :: vx

c     Common variables:
      INTEGER :: jmod, jb, jnn, jter, joc, jkh

      DOUBLE PRECISION :: ml, ql, xl, uc

      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH

c-----------------------------------------------------------------------
c     Explicitly initialize all local variables (Steve, 5/08).

      i = 0
      jj = 0
      jjj = 0
      vx = 0
c-----------------------------------------------------------------------

      !IF ( JKH >= 2 ) CALL PRINTS ( K, JKH, IS1, IS2, JCL )
      IF ( JS1 > JS2 .OR. 1 > IS2 ) RETURN
      IF ( JS4 <= JS5 ) THEN
         DO JJ = JS1, JS2
            JJJ = JJ - JS1 + IS3
            DO I = IS1, IS2
               VX = S(JJ, I)
               S(JS4:JS5, I) = S(JS4:JS5, I) - VX*C(K, JJJ, 1:JS5 - JS4 + 1)
            END DO
         END DO
      END IF
      DO JJ = JS1, JS2
         JJJ = JJ - JS1 + IS3
         DO I = IS1, IS2
            VX = S(JJ, I)
            S(KJ10:KJ12, I) = S(KJ10:KJ12, I) - VX*C(K, JJJ, KJ10 - KJ12 + KI4:KI4)
         END DO
      END DO
      RETURN
      END

      SUBROUTINE PRINTS ( JK, JKH, IPS1, IPS2, JCL )
      USE MESH
      USE SOLVER_GLOBAL_VARIABLES
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, KW(260)
      WRITE (10,*) 'printS', JCL, JK, JKH, KEE, KEQ, KJ2, KJ5, KJ6, KJ12
      IF ( KEE /= 2 ) THEN
         DO I = IPS1, IPS2
            WRITE (10, 991) JK, (S(J, I), J = 1, KEQ)
         END DO
         WRITE (10, 991) JK
      END IF
      DO I = IPS1, IPS2
         WRITE (10, 991) JK, (S(J, I), J = KJ2, KJ5)
      END DO
      WRITE (10, 991) JK
      DO I = IPS1, IPS2
         WRITE (10, 991) JK, (S(J, I), J = KJ6, KJ12)
      END DO
      WRITE (10, 991) JK
  991 FORMAT (I5, 1P, 17D9.1, /, 17D9.1)
      RETURN
      END

      SUBROUTINE PRINTC ( JK, JKH, IPC1, IPC2, JPC1, JPC2 )
      USE MESH
      USE SOLVER_GLOBAL_VARIABLES
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, KW(260)
      IF ( JKH < 2 .OR. (JK >= 4 .AND. JK <= KH - 2) ) RETURN
      WRITE (10,*) 'printC', JK, JKH, JPC1, JPC2
      DO I = IPC1, IPC2
         WRITE (10, 99002) (C(JK, I, J), J = JPC1, JPC2)
      END DO
99002 FORMAT (1P, 15D10.2)
      RETURN
      END
