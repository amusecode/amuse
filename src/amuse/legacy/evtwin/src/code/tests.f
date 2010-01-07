!     Compare the output of NUCRAT and NUCRAT2 for the reactions that exist
!     in both
      SUBROUTINE COMPARE_NUCRAT
      USE MESH
      USE NUCLEOSYNTHESIS
      USE CONSTANTS
      IMPLICIT NONE

      INTEGER :: K, I
c     Common variables:
      INTEGER :: kh, ktw, kw
      INTEGER :: jmod, jb, jnn, jter, joc, jkh

      DOUBLE PRECISION :: h, dh, eps, del, dh0
      DOUBLE PRECISION :: ml, ql, xl, uc
      DOUBLE PRECISION :: var, dvar, fn1, dfn1
      DOUBLE PRECISION :: fn2, dfn2, equ, dequ
      DOUBLE PRECISION :: VARNUC, DVARNUC, FN1NUC, DFN1NUC, RRT1(20)
!     /STAT2 
      DOUBLE PRECISION ::  AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA
     $     ,SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, RPP, R33, R34, RBE,
     $     RBP,RPC, RPN, RPO, R3A, RAC, RAN, RAO, RANE, RCC, RCO, ROO,
     $     RGNE,RGMG, RCCG, RPNG, EX, ENX, WMU, DELTA, PHIE, EXT, FKT,
     &     FKR, PRANDTL
      DOUBLE PRECISION :: RRT
!     /RATES /
      DOUBLE PRECISION :: RRT2
!     /ABUND /
      DOUBLE PRECISION :: XA, NA, NE, NI, NZZ, AVM, NE1
      DOUBLE PRECISION :: XA2, NA2
      DOUBLE PRECISION :: CH2(4), CHI(26,9), COM(27), CAN(9), CBN(9)
      INTEGER ::  KZN(9)

      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, KW(260)
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /INF   / VAR(NVAR), DVAR(NVAR), FN1(NFUNC), DFN1(NVAR,NFUNC)
      COMMON /INF2  / VARNUC(50), DVARNUC(50), FN1NUC(110), DFN1NUC(50, 110)
      COMMON /INE   / FN2(3,NFUNC), DFN2(3,NVAR,NFUNC), EQU(NEQ), DEQU(NVAR,3,NEQ)
      COMMON /STAT2 / AP, ARHO, U, P, RHO, FK, T, SF, ST, ZT, GRADA, 
     : SCP, RF, RT, XHI, S, PR, PG, PF, PT, EN, RRT(20), EX, ENX, WMU, DELTA, PHIE, EXT, FKT, FKR, PRANDTL
      COMMON /RATES / RRT2(90)
      COMMON /ABUND / XA(9), NA(9), NE, NI, NZZ, AVM, NE1
      COMMON /ABUND2/ XA2(50), NA2(50)
      COMMON /ATDATA/ CH2, CHI, COM, CAN, CBN, KZN

      DO K=1, KH
!        Setup COMMON blocks: structure
         DVAR(1:NVAR) = DH(1:NVAR, K)
          VAR(1:NVAR) =  H(1:NVAR, K) + DVAR(1:NVAR)

         VAR(2) = VAR(2) + LOG(12.0)

!        Setup COMMON blocks: nucleosynthesis
         DVARNUC(1:50) = DHNUC(1, 1:50,K)
         VARNUC(1:50) = HNUC(1, 1:50,K) + DVARNUC(1:50)

!        Calculate the structure
         CALL FUNCS1 ( K, 0 )
         RRT1(:) = RRT(:)

!        Calculate the nucleosynthesis
         CALL FUNCS2 ( K, 1, KH )

!        Compare the output
   95    FORMAT (1X, 1P, I6, 6D15.8)
         I = 3
         write (*, 95) K, H(2, K)/CLN, H(5, K)*CAN(1)/CBN(1)/AVM, HNUC(1,41, K), RRT1(I), RRT2(I)
      END DO

      END SUBROUTINE
