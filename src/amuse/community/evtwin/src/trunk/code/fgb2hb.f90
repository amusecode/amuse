module resolve_helium_flash
  use real_kind

  implicit none
  integer, save, private :: interpolation_table_size = 0;
  real(double), save, private, allocatable :: pre_flash_mass(:)
  real(double), save, private, allocatable :: pre_flash_h(:,:)
  real(double), save, private, allocatable :: pre_flash_he(:,:)
  real(double), save, private, allocatable :: pre_flash_c(:,:)
  real(double), save, private, allocatable :: pre_flash_n(:,:)
  real(double), save, private, allocatable :: pre_flash_o(:,:)
  real(double), save, private, allocatable :: pre_flash_ne(:,:)
  real(double), save, private, allocatable :: pre_flash_mg(:,:)


contains

   !     Save the pre-flash composition profile
   subroutine store_pre_flash_composition
      use real_kind
      use interpolate
      use indices
      use mesh
      use indices
      implicit none

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
      pre_flash_mass(1:KH)  = H(VAR_MASS, 1:KH)
      pre_flash_h(1, 1:KH)  = H(VAR_H1,   1:KH)
      pre_flash_he(1, 1:KH) = H(VAR_HE4,  1:KH)
      pre_flash_c(1, 1:KH)  = H(VAR_C12,  1:KH)
      pre_flash_n(1, 1:KH)  = H(VAR_N14,  1:KH)
      pre_flash_o(1, 1:KH)  = H(VAR_O16,  1:KH)
      pre_flash_ne(1, 1:KH) = H(VAR_NE20, 1:KH)
      pre_flash_mg(1, 1:KH) = H(VAR_MG24, 1:KH)

      !     Accretion abundances
      call iptable_init (interpolation_table_size, pre_flash_mass(:),  &
           pre_flash_h(1, :), pre_flash_h(2, :),  &
           pre_flash_h(3, :), pre_flash_h(4, :))

      call iptable_init (interpolation_table_size, pre_flash_mass(:),  &
           pre_flash_he(1, :), pre_flash_he(2, :),  &
           pre_flash_he(3, :), pre_flash_he(4, :))

      call iptable_init (interpolation_table_size, pre_flash_mass(:),  &
           pre_flash_c(1, :), pre_flash_c(2, :),  &
           pre_flash_c(3, :), pre_flash_c(4, :))

      call iptable_init (interpolation_table_size, pre_flash_mass(:),  &
           pre_flash_n(1, :), pre_flash_n(2, :),  &
           pre_flash_n(3, :), pre_flash_n(4, :))

      call iptable_init (interpolation_table_size, pre_flash_mass(:),  &
           pre_flash_o(1, :), pre_flash_o(2, :),  &
           pre_flash_o(3, :), pre_flash_o(4, :))

      call iptable_init (interpolation_table_size, pre_flash_mass(:),  &
           pre_flash_ne(1, :), pre_flash_ne(2, :),  &
           pre_flash_ne(3, :), pre_flash_ne(4, :))

      call iptable_init (interpolation_table_size, pre_flash_mass(:),  &
           pre_flash_mg(1, :), pre_flash_mg(2, :),  &
           pre_flash_mg(3, :), pre_flash_mg(4, :))

   end subroutine store_pre_flash_composition



   !>     Set the accretion abundance as appropriate for the current mass of
   !!     the star, based on the pre-flash model
   !<
   subroutine update_accretion_abundance
      use real_kind
      use interpolate
      use accretion_abundances
      use mesh
      use indices

      implicit none
      real(double) :: m

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

      m = H(VAR_MASS, 1)
      if (m >= pre_flash_mass(1)) return;

      xac(1, 1) = iptable_eval(interpolation_table_size, m, pre_flash_mass(:),  &
           pre_flash_h(1, :), pre_flash_h(2, :),  &
           pre_flash_h(3, :), pre_flash_h(4, :))

      xac(2, 1) = iptable_eval(interpolation_table_size, m, pre_flash_mass(:),  &
           pre_flash_he(1, :), pre_flash_he(2, :),  &
           pre_flash_he(3, :), pre_flash_he(4, :))

      xac(3, 1) = iptable_eval(interpolation_table_size, m, pre_flash_mass(:),  &
           pre_flash_c(1, :), pre_flash_c(2, :),  &
           pre_flash_c(3, :), pre_flash_c(4, :))

      xac(4, 1) = iptable_eval(interpolation_table_size, m, pre_flash_mass(:),  &
           pre_flash_n(1, :), pre_flash_n(2, :),  &
           pre_flash_n(3, :), pre_flash_n(4, :))

      xac(5, 1) = iptable_eval(interpolation_table_size, m, pre_flash_mass(:),  &
           pre_flash_o(1, :), pre_flash_o(2, :),  &
           pre_flash_o(3, :), pre_flash_o(4, :))

      xac(6, 1) = iptable_eval(interpolation_table_size, m, pre_flash_mass(:),  &
           pre_flash_ne(1, :), pre_flash_ne(2, :),  &
           pre_flash_ne(3, :), pre_flash_ne(4, :))

      xac(7, 1) = iptable_eval(interpolation_table_size, m, pre_flash_mass(:),  &
           pre_flash_mg(1, :), pre_flash_mg(2, :),  &
           pre_flash_mg(3, :), pre_flash_mg(4, :))

     end subroutine update_accretion_abundance



     !     Clean up allocated memory
     subroutine cleanup_pre_flash_composition
      use real_kind
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
   end subroutine cleanup_pre_flash_composition

   !> ------------------------------------------------------------------------------
   !!  FGB2HB
   !!   Construct a horizontal branch (HB) model based on a first giant branch
   !!   model (FGB) that is close to the helium flash.
   !!   Preserves the total mass and composition of the star.
   !!
   !!  \todo  FIXME: only really works for single stars.
   !!
   !! ------------------------------------------------------------------------------
   !!  Input:
   !!   JOP - FORTRAN file handle where output model is stored (for STAR12)
   !!  Output:
   !!   JO3 - Return code from STAR12, to indicate convergence success/failure
   !!         Returns 13 to indicate that the ZAHB construction was successful
   !!   COMMON H(:,:) - the post-flash core He burning model
   !! ------------------------------------------------------------------------------
   !<

   subroutine fgb2hb ( jop, jo3 )
      use real_kind
      use mesh
      use mesh_enc
      use control
      use test_variables
      use current_model_properties
      use constants
      use init_dat
      use stopping_conditions
      
      implicit none
      integer :: jop, jo3, wanted_kh, kt5, jch
      integer :: cmi_mode_bck,ksv
      real(double) :: bm_fgb, bper_fgb
      
      if(debug_level.ge.2) write(6,'(/,A20,2I12)')'fgb2hb: ',jop,jo3
      
      ! Evolve a standard ZAHB model (stored on fort.12, with init.dat in fort.24)
      ! to the required ZAHB model: same total mass (SM) and core mass (VMH)
      uc(13) = sm

      !> \todo FIXME: the core mass must be >= the mass of the ZAHB construction model
      !! Probably doen't matter much in practice.
      !< 
      uc(14) = max(0.40d0, mh)
      rewind (jop)
      
      ! Save some variables to use on HB:
      cmi_mode_bck = cmi_mode
      cmi_mode = 1
      bm_fgb = bm
      bper_fgb = bper

      if (read_init_dat(24, wanted_kh, ksv, kt5, jch) .eqv. .false.) then
         jo3 = -3
         return
      end if
      rewind (24)

      call store_pre_flash_composition
      call star12 ( jo3, wanted_kh, jch, jop, 12, ksv, kt5, 24 )
      call cleanup_pre_flash_composition

      if (read_init_dat(22, wanted_kh, ksv, kt5, jch) .eqv. .false.) then
         jo3 = -3
         return
      end if
      rewind (22)

      ! Restore variables to FGB values:
      cmi_mode = cmi_mode_bck
      bm = bm_fgb
      bper = bper_fgb
      
      uc(13) = 1.0d3
      uc(14) = 1.0d3

      rewind (22)
      rewind (12)
      rewind (24)
      rewind (jop)

   end subroutine fgb2hb

end module resolve_helium_flash
