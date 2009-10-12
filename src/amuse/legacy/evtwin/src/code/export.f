! Export a stellar model.
! This module is disgned to work along with the MUSE library and exports
! the data for a stellar model, to be used alongside a tool like
! MakeMeAStar.
! Order of variables stored at each meshpoint:
!  Mass coordinate [Msun], Radius [Rsun], log density [cgs],
!  log pressure [cgs], XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE
! Mespoints will be passed out from *surface* to *centre*

      function get_number_of_meshpoints()
      use mesh
      implicit none
      integer :: get_number_of_meshpoints
      double precision :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      integer :: KH, KTW, KW(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW
      
      get_number_of_meshpoints = KH
      end function



      subroutine export_stellar_model(jstar, model)
      use mesh
      implicit none
      double precision :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      integer :: KH, KTW, KW(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW

      double precision :: LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,WL,
     &     WCV, HP, WT, PHIM, GMR, SEP, M3, PX(NPX), SX(NPX,NM+1),QA(NM)
      COMMON /VBLES / LOLEDD, DG, EG, GRAD, ETH, EGR, R, QQ, QM,
     & WL, WCV, HP, WT, PHIM, GMR, SEP, M3, PX, SX, QA

      double precision :: CH2(4), CHI(26,9), COM(27), CAN(9), CBN(9), KZN(9)
      COMMON /ATDATA/ CH2, CHI, COM, CAN, CBN, KZN

      double precision :: XA(9), NA(9)
      double precision :: NEO, NIO, NZZ, AVM, NE
      COMMON /ABUND / XA, NA, NEO, NIO, NZZ, AVM, NE

      double precision :: HPR(NVAR,NM), HT(60025+4*NM)
      COMMON /STORE / HPR, HT

      integer, intent(in) :: jstar     ! Which component of a binary
      double precision, intent(out) :: model(13,KH)
      double precision :: HH(NVAR,NM)

      integer :: ik, ikk, n

      ! Should we export the *current* or *previous* model?
      ! FIXME: make this an option! At the moment we output the previous
      ! model.
      HH(:,:) = H(:,:)
      H(:,:) = HPR(:,:)

      CALL COMPUTE_OUTPUT_QUANTITIES ( JSTAR )

      do ik=1, KH
         ikk = 1 + (KH+1-ik)
         model(1, ik) = SX(9, ikk)
         model(2, ik) = SX(17, ikk)
         model(3, ik) = log(SX(3, ikk))
         model(4, ik) = log(SX(2, ikk))
         ! Convert *all* abundances to mass fractions
         XA(1) = H(5, ik)
         XA(2) = H(9, ik)
         XA(3) = H(10, ik)
         XA(4) = H(16, ik)
         XA(5) = H(3, ik)
         XA(6) = H(11, ik)
         XA(7) = 1.0D0 - SUM(XA(1:6)) - SUM(XA(8:9))
         do n=1, 9
            NA(n) = XA(n) * CAN(n)/CBN(n)
         end do
         AVM = SUM(NA(1:9))
         NA(1:9) = NA(1:9) / AVM
         model(5:13, ik) = NA(1:9)
         !model(4:10, ik) = SX(10:16, ikk)
         !model(11:12, ik) = XA(8:9)    ! Fe and Si: constant
      end do
      H(:,:) = HH(:,:)

      end subroutine
