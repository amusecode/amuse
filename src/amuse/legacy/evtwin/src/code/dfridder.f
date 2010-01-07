      subroutine dfridders(IVAR, IEE, JX, JY, JZ, K)
      use mesh
      use solver_global_variables
      implicit none
      integer, intent(in) :: IVAR, IEE, K, JX, JY, JZ
      integer, parameter :: NTAB = 5
      double precision, parameter :: CON = 2.0d0**0.5
      double precision, parameter :: CON2 = CON**2
      double precision, parameter :: BIG = 1.0d30
      double precision, parameter :: SAFE = 2.0d0

!     Local variables
      integer :: iloop, jloop, ji, jeq, ieq, JK
      double precision :: vx, dvx, dx, invDX
      double precision :: err(NEQ), errt(NEQ), fac
      double precision :: a(NTAB, NTAB, NEQ)
      double precision :: DE(NEQ)

c     Common variables:
      double precision :: H, DH, EPS, DEL, DH0
      integer :: KH, KTW, KW
      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, KW(260)

      integer :: JMOD, JB, JNN, JTER, JOC, JKH
      double precision :: ML, QL, XL, UC
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH

      double precision :: VAR, DVAR, FN1, DFN1
      COMMON /INF   / VAR(NVAR), DVAR(NVAR), FN1(NFUNC), DFN1(NVAR,NFUNC)

      double precision :: FN2, DFN2, EQU, DEQU
      COMMON /INE   / FN2(3,NFUNC), DFN2(3,NVAR,NFUNC), EQU(NEQ), DEQU(NVAR,3,NEQ)
      double precision :: old_FN2(NFUNC)!, old_DFN2(3,NVAR,NFUNC), old_EQU(NEQ)
      double precision :: dummy

      old_FN2(:) = FN2(IEE,:)
      !old_DFN2(:,:,:) = DFN2(:,:,:)
      JK = K - KL
      JK = JK + (3 - IEE)
      
      if (JK > KH) JK = KH
      if (JK < 1) JK = 1

!     foreach variable
      !if (K > KH) then
         DVAR(1:NVAR) = DH(1:NVAR, JK)
          VAR(1:NVAR) =  H(1:NVAR, JK) + DVAR(1:NVAR)
      !else
      !   DVAR(1:NVAR) = DH(1:NVAR, K - KL - 3+KEE)
      !    VAR(1:NVAR) =  H(1:NVAR, K - KL - 3+KEE) + DVAR(1:NVAR)
      !end if
      CALL FUNCS1 ( JK, 0 )

      JI = KD(IVAR)
      VX = VAR(JI)
      DVX = DVAR(JI)
      DX = SIGN(DH0*MAX(ABS(VX), 1.0D0), DVX)
      !dummy = VX + DX
      !call void(dummy)
      !DX = dummy - VX

!     compute funcs(x-dx) and equns(x-dx)
      !invDX = 0.5d0 * ER(JI)/DX
      !VAR(JI) = VX - DX
      !DVAR(JI) = DVX - DX
      invDX = ER(JI)/DX
      VAR(JI) = VX
      DVAR(JI) = DVX
      CALL FUNCS1 ( JK, JI )
      FN2(IEE, 1:NFUNC) = FN1(1:NFUNC)
      CALL EQUNS1 ( K, KL, KQ )
      do ieq = JX, JY
         jeq = KD(ieq + JZ)
         DE(jeq) = EQU(jeq)
      end do

!     compute funcs(x+dx) and equns(x+dx)
      VAR(JI) = VX + DX
      DVAR(JI) = DVX + DX
      CALL FUNCS1 ( JK, JI )
      FN2(IEE, 1:NFUNC) = FN1(1:NFUNC)
      CALL EQUNS1 ( K, KL, KQ )

!     store all derivatives
      do ieq = JX, JY
         jeq = KD(ieq + JZ)
         a(1,1, jeq) = -(DE(jeq) - EQU(jeq)) * invDX
         DEQU(IVAR, IEE, jeq) = a(1, 1, jeq)
      end do
      !return
      err(:)=BIG
      do iloop=2,NTAB
         DX = DX/CON
         invDX = invDx * CON
!     compute funcs(x-dx) and equns(x-dx)
         !VAR(JI) = VX - DX
         !DVAR(JI) = DVX - DX
         !CALL FUNCS1 ( K - KL, JI )
         !FN2(IEE, 1:NFUNC) = FN1(1:NFUNC)
         !CALL EQUNS1 ( K, KL, KQ )
         !do ieq = JX, JY
         !   jeq = KD(ieq + JZ)
         !   DE(jeq) = EQU(jeq)
         !end do

!     compute funcs(x+dx) and equns(x+dx)
         VAR(JI) = VX + DX
         DVAR(JI) = DVX + DX
         CALL FUNCS1 ( K - KL, JI )
         FN2(IEE, 1:NFUNC) = FN1(1:NFUNC)
         CALL EQUNS1 ( K, KL, KQ )

!     store all derivatives
         do ieq = JX, JY
            jeq = KD(ieq + JZ)
            a(1,iloop,jeq) = -(DE(jeq) - EQU(jeq)) * invDX
         end do

         fac=CON2
         do jloop=2,iloop
            do ieq = JX, JY
               jeq = KD(ieq + JZ)
               a(jloop,iloop,jeq)=(a(jloop-1,iloop,jeq)*fac - a(jloop-1,iloop-1,jeq))/(fac-1.0d0)
            end do
            fac = fac * CON2
            do ieq = JX, JY
               jeq = KD(ieq + JZ)
               errt(jeq) = max(abs(a(jloop,iloop,jeq) - a(jloop-1,iloop, jeq)),
     &                         abs(a(jloop,iloop,jeq) - a(jloop-1,iloop-1,jeq)))
               if (errt(jeq) <= err(jeq)) then
                  err(jeq) = errt(jeq)
                  ! Store result
                  DEQU(IVAR, IEE, jeq) = a(jloop, iloop, jeq)
               endif
            end do
         end do
         !if(abs(a(iloop,iloop)-a(iloop-1,iloop-1)).ge.SAFE*err)return
      end do
!     end foreach variable

      FN2(IEE,:) = old_FN2(:)
      !DFN2(:,:,:) = old_DFN2(:,:,:)
      return
      end subroutine

      subroutine void(x)
      implicit none
      double precision, intent(inout) :: x
      x = x*1.0d0
      end subroutine
