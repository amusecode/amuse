!     Code for making sure the composition profile of the post He-flash
!     model corresponds to the pre-flash model.
      module fgb2hb_composition

      integer, save :: interpolation_table_size = 0;
      double precision, save, allocatable :: pre_flash_mass(:)
      double precision, save, allocatable :: pre_flash_h(:,:)
      double precision, save, allocatable :: pre_flash_he(:,:)
      double precision, save, allocatable :: pre_flash_c(:,:)
      double precision, save, allocatable :: pre_flash_n(:,:)
      double precision, save, allocatable :: pre_flash_o(:,:)
      double precision, save, allocatable :: pre_flash_ne(:,:)
      double precision, save, allocatable :: pre_flash_mg(:,:)

      contains

!     Save the pre-flash composition profile
      subroutine store_pre_flash_composition
      use mesh_enc
      implicit none
!     Unnamed COMMON block, for getting the composition profile before the
!     flash
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, KW(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW
!     Common block for elemental abundances
      DOUBLE PRECISION :: XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, XW(14)
      COMMON /ABUND / XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, XW

      interpolation_table_size = kh;
      allocate( pre_flash_mass(kh) )
      allocate( pre_flash_h(4, kh) )
      allocate( pre_flash_he(4, kh) )
      allocate( pre_flash_c(4, kh) )
      allocate( pre_flash_n(4, kh) )
      allocate( pre_flash_o(4, kh) )
      allocate( pre_flash_ne(4, kh) )
      allocate( pre_flash_mg(4, kh) )

!     Construct tables for interpolation
      pre_flash_mass(1:KH) = H(4, 1:KH)
      pre_flash_h(1, 1:KH) = H(5, 1:KH)
      pre_flash_he(1, 1:KH) = H(9, 1:KH)
      pre_flash_c(1, 1:KH) = H(10, 1:KH)
      pre_flash_n(1, 1:KH) = H(16, 1:KH)
      pre_flash_o(1, 1:KH) = H(3, 1:KH)
      pre_flash_ne(1, 1:KH) = H(11, 1:KH)
      pre_flash_mg(1, 1:KH) = XMG!H(24, 1:KH)

!     Accretion abundances
      call iptable_init (interpolation_table_size, pre_flash_mass(:),
     &      pre_flash_h(1, :), pre_flash_h(2, :),
     &      pre_flash_h(3, :), pre_flash_h(4, :))

      call iptable_init (interpolation_table_size, pre_flash_mass(:),
     &      pre_flash_he(1, :), pre_flash_he(2, :),
     &      pre_flash_he(3, :), pre_flash_he(4, :))

      call iptable_init (interpolation_table_size, pre_flash_mass(:),
     &      pre_flash_c(1, :), pre_flash_c(2, :),
     &      pre_flash_c(3, :), pre_flash_c(4, :))

      call iptable_init (interpolation_table_size, pre_flash_mass(:),
     &      pre_flash_n(1, :), pre_flash_n(2, :),
     &      pre_flash_n(3, :), pre_flash_n(4, :))

      call iptable_init (interpolation_table_size, pre_flash_mass(:),
     &      pre_flash_o(1, :), pre_flash_o(2, :),
     &      pre_flash_o(3, :), pre_flash_o(4, :))

      call iptable_init (interpolation_table_size, pre_flash_mass(:),
     &      pre_flash_ne(1, :), pre_flash_ne(2, :),
     &      pre_flash_ne(3, :), pre_flash_ne(4, :))

      call iptable_init (interpolation_table_size, pre_flash_mass(:),
     &      pre_flash_mg(1, :), pre_flash_mg(2, :),
     &      pre_flash_mg(3, :), pre_flash_mg(4, :))

!     Abundances for regular composition adjustment during the run
!     We set the above to accrete material of the right composition, we
!     set these to counter the effect of numerical diffusion.
      TH(4, 1:KH) = H(4, 1:KH)
      TH(5, 1:KH) = H(5, 1:KH)
      TH(9, 1:KH) = H(9, 1:KH)
      TH(10, 1:KH) = H(10, 1:KH)
      TH(16, 1:KH) = H(16, 1:KH)
      TH(3, 1:KH) = H(3, 1:KH)
      TH(11, 1:KH) = H(11, 1:KH)
      call iptable_init (kh, TH(4,:), TH(5,:), THb(5,:), THc(5,:), THd(5,:))
      call iptable_init (kh, TH(4,:), TH(9,:), THb(9,:), THc(9,:), THd(9,:))
      call iptable_init (kh, TH(4,:), TH(10,:), THb(10,:), THc(10,:), THd(10,:))
      call iptable_init (kh, TH(4,:), TH(16,:), THb(16,:), THc(16,:), THd(16,:))
      call iptable_init (kh, TH(4,:), TH(3,:), THb(3,:), THc(3,:), THd(3,:))
      call iptable_init (kh, TH(4,:), TH(11,:), THb(11,:), THc(11,:), THd(11,:))

      end subroutine

!     Set the accretion abundance as appropriate for the current mass of
!     the star, based on the pre-flash model
      subroutine update_accretion_abundance
      use mesh_enc
      implicit none
!     Unnamed COMMON block, for getting the composition profile before the
!     flash
      DOUBLE PRECISION :: H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0
      INTEGER :: KH, KTW, KW(260)
      COMMON H, DH, EPS, DEL, DH0, KH, KTW, KW
!     Common block for accretion abundances
      DOUBLE PRECISION :: X1AC, X4AC, X12AC, X14AC, X16AC, X20AC, X24AC
      DOUBLE PRECISION :: XAC(7, 2)
      COMMON /ACCRET/ X1AC, X4AC, X12AC, X14AC, X16AC, X20AC, X24AC, XAC

      double precision :: m

!     Break out if we're not supposed to do anything
      if (interpolation_table_size == 0) return

!     Old-style accretion: set only the surface composition
!     This is good enough for the envelope composition, maybe not for
!     seismology. It seems to be a little more numerically stable,
!     however.
      XAC(1, 1) = pre_flash_h(1, 1)
      XAC(2, 1) = pre_flash_he(1, 1)
      XAC(3, 1) = pre_flash_c(1, 1)
      XAC(4, 1) = pre_flash_n(1, 1)
      XAC(5, 1) = pre_flash_o(1, 1)
      XAC(6, 1) = pre_flash_ne(1, 1)
      XAC(7, 1) = pre_flash_mg(1, 1)
      return

      m = H(4, 1)
      if (m >= pre_flash_mass(1)) return;

      XAC(1, 1) = iptable_eval(interpolation_table_size, m, pre_flash_mass(:),
     &      pre_flash_h(1, :), pre_flash_h(2, :),
     &      pre_flash_h(3, :), pre_flash_h(4, :))

      XAC(2, 1) = iptable_eval(interpolation_table_size, m, pre_flash_mass(:),
     &      pre_flash_he(1, :), pre_flash_he(2, :),
     &      pre_flash_he(3, :), pre_flash_he(4, :))

      XAC(3, 1) = iptable_eval(interpolation_table_size, m, pre_flash_mass(:),
     &      pre_flash_c(1, :), pre_flash_c(2, :),
     &      pre_flash_c(3, :), pre_flash_c(4, :))

      XAC(4, 1) = iptable_eval(interpolation_table_size, m, pre_flash_mass(:),
     &      pre_flash_n(1, :), pre_flash_n(2, :),
     &      pre_flash_n(3, :), pre_flash_n(4, :))

      XAC(5, 1) = iptable_eval(interpolation_table_size, m, pre_flash_mass(:),
     &      pre_flash_o(1, :), pre_flash_o(2, :),
     &      pre_flash_o(3, :), pre_flash_o(4, :))

      XAC(6, 1) = iptable_eval(interpolation_table_size, m, pre_flash_mass(:),
     &      pre_flash_ne(1, :), pre_flash_ne(2, :),
     &      pre_flash_ne(3, :), pre_flash_ne(4, :))

      XAC(7, 1) = iptable_eval(interpolation_table_size, m, pre_flash_mass(:),
     &      pre_flash_mg(1, :), pre_flash_mg(2, :),
     &      pre_flash_mg(3, :), pre_flash_mg(4, :))

      end subroutine

!     Clean up allocated memory
      subroutine cleanup_pre_flash_composition
      implicit none

      interpolation_table_size = 0

      deallocate( pre_flash_mass )
      deallocate( pre_flash_h )
      deallocate( pre_flash_he )
      deallocate( pre_flash_c )
      deallocate( pre_flash_n )
      deallocate( pre_flash_o )
      deallocate( pre_flash_ne )
      deallocate( pre_flash_mg )
      end subroutine

      end module
