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

   subroutine make_post_flash_model(jo)
      use real_kind
      use mesh
      use mesh_enc
      use control
      use settings
      use test_variables
      use current_model_properties
      use constants
      use init_dat
      use stopping_conditions
      use indices
      use nucleosynthesis, only: hnuc
      use polytrope
      use interpolate
      
      implicit none
      integer, intent(inout) :: jo
      integer :: jo3, wanted_kh, kt5, jch, kh2
      integer :: cmi_mode_bck,ksv
      real(double) :: bm_fgb, bper_fgb, ccmi, xa(9)
      real(double) :: sm1,dty1,age1,per1,bms1,ecc1,p1,enc1, tm, oa
      integer      :: kh1,kp1,jmod1,jb1,jn1,jf1
      integer      :: it, jk
      real(double) :: hn1(50, nm)
      type(interpolate_t) :: pf_h, pf_he, pf_c, pf_n, pf_o, pf_ne, pf_mg, pf_si, pf_fe
      real(double) :: mcore, mstar
      real(double) :: age_backup
      integer      :: jmod_backup, jnn_backup

      ! Evolve a standard ZAHB model (stored on fort.12, with init.dat in fort.24)
      ! to the required ZAHB model: same total mass (SM) and core mass (VMH)
      uc(13) = sm

      !> \todo FIXME: the core mass must be >= the mass of the ZAHB construction model
      !! Probably doen't matter much in practice.
      !< 
      mcore = max(0.40d0, mh)
      mstar = sm

      ! Backup age, model number and number of models, to be restored after the flash
      age_backup  = age
      jmod_backup = jmod
      jnn_backup  = jnn

      ! Save some variables to use on HB:
      cmi_mode_bck = cmi_mode
      cmi_mode = 1
      bm_fgb = bm
      bper_fgb = bper
      kh2 = kh

      if (read_init_dat(24, wanted_kh, ksv, kt5, jch) .eqv. .false.) then
         jo3 = -3
         return
      end if
      rewind (24)
      kr1 = 200
      kr2 = 50
      kt5 = max(kr1, kr2)
      !kt4 = 10000
      ccmi = cmi
      cmi = 0.0d0
      ky = 0
      !artmix = 1.0d-8

      call store_pre_flash_composition
      call make_interpolation_table(kh, H(VAR_MASS, 1:KH), H(VAR_H1,   1:KH), pf_h)
      call make_interpolation_table(kh, H(VAR_MASS, 1:KH), H(VAR_HE4,  1:KH), pf_he)
      call make_interpolation_table(kh, H(VAR_MASS, 1:KH), H(VAR_C12,  1:KH), pf_c)
      call make_interpolation_table(kh, H(VAR_MASS, 1:KH), H(VAR_N14,  1:KH), pf_n)
      call make_interpolation_table(kh, H(VAR_MASS, 1:KH), H(VAR_O16,  1:KH), pf_o)
      call make_interpolation_table(kh, H(VAR_MASS, 1:KH), H(VAR_NE20, 1:KH), pf_ne)
      call make_interpolation_table(kh, H(VAR_MASS, 1:KH), H(VAR_MG24, 1:KH), pf_mg)
      call make_interpolation_table(kh, H(VAR_MASS, 1:KH), H(VAR_SI28, 1:KH), pf_si)
      call make_interpolation_table(kh, H(VAR_MASS, 1:KH), H(VAR_FE56, 1:KH), pf_fe)

      do jk = 1, kh
         if (H(VAR_MASS, jk) < mcore*CMSN) cycle;
         H(VAR_H1,   jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_h)
         H(VAR_HE4,  jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_he)
         H(VAR_C12,  jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_c)
         H(VAR_N14,  jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_n)
         H(VAR_O16,  jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_o)
         H(VAR_NE20, jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_ne)
         H(VAR_MG24, jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_mg)
         H(VAR_SI28, jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_si)
         H(VAR_FE56, jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_fe)
      end do

      ! Construct a suitable core
      oa    = h(VAR_HORB,  1)
      xa(1) = H(VAR_H1,   1)
      xa(2) = H(VAR_HE4,  1)
      xa(3) = H(VAR_C12,  1)
      xa(4) = H(VAR_N14,  1)
      xa(5) = H(VAR_O16,  1)
      xa(6) = H(VAR_NE20, 1)
      xa(7) = H(VAR_MG24, 1)
      xa(8) = H(VAR_SI28, 1)
      xa(9) = H(VAR_FE56, 1)
      tm = 2.5d0
      uc(14) = 0.40
      call generate_starting_model(2.5d0, h, dh, hn1, sm1,dty1,age1,per1,bms1,ecc1,p1,enc1,kh1,kp1,jmod1,jb1,jn1,jf1)
      call remesh ( kh2, 4, tm*CMSN+0.1, tm*CMSN, p1, 0.0d0, 10.*oa, 1, 2 )
      hpr = h
      call star12_loop ( jo, ksv, kt5, kp1, 24, dty1 )

      age = 0.0d0
      dty1 = 1.0d3
      cmi = -1.0d-6
      artmix = 0.0d0
      jmod = 0
      jo = 0
      call remesh ( kh2, 3, tm*CMSN+0.1, tm*CMSN, p1, 0.0d0, 10.*oa, 1, 2 )

      it = 24
      call printb ( jo, 1, it)
      call nextdt ( dty1, jo, it )

      ! Begin evolutionary loop of KP time steps:
      jnn = 1
      uc(13) = 1.0d3
      cmi = -1.0d-6 / CSY
      uc(14) = mcore
      kx = 1
      do_kp: do
         jo = 0
         ! Solve for structure, mesh, and major composition variables
         joc = 1

         if (h(VAR_MASS, 1)/CMSN < 0.25d0 * (3.d0*mcore + mstar)) then
            cmi = ccmi
            uc(13) = mstar
            uc(14) = mcore
            ccac = 1.0d0
         end if
         if (abs(cmi) <= 1.0d-100) then
            h(VAR_MASS, 1:kh) = h(VAR_MASS, 1:kh) / h(VAR_MASS, 1) * mstar * CMSN
         end if
         call smart_solver ( kr2, id, kt5, jo )

         ! If no convergence, restart from 2 steps back, DT decreased substantially
         if (jo /= 0) then
            call backup ( dty1, jo )
            if ( jo.eq.2 ) exit do_kp    ! Abort if timestep below limit
            goto 4
         end if

         ! If model didn't converge, give up
         if ( jo >= 1 ) exit do_kp

         call printb ( jo, 1, it )

         if ( jo >= 2 ) exit do_kp
         call update ( dty1 )

   4     continue
         call nextdt ( dty1, jo, it )
         if ( jo == 3 ) exit do_kp
         jnn = jnn + 1
      end do do_kp

      dty1 = 1.0d3
      do jk = 1, kh
         if (H(VAR_MASS, jk) < mcore / CMSN) then
         H(VAR_H1,   jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_h)
         H(VAR_HE4,  jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_he)
         H(VAR_C12,  jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_c)
         H(VAR_N14,  jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_n)
         H(VAR_O16,  jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_o)
         H(VAR_NE20, jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_ne)
         H(VAR_MG24, jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_mg)
         H(VAR_SI28, jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_si)
         H(VAR_FE56, jk) = evaluate_interpolation_table(H(VAR_MASS, jk), pf_fe)
         end if
      end do
      hpr = h
      dh = 0
      !call star12_loop ( jo, ksv, kt5, kp1, 24, dty1 )

      call destroy_interpolation_table(pf_h)
      call destroy_interpolation_table(pf_he)
      call destroy_interpolation_table(pf_c)
      call destroy_interpolation_table(pf_n)
      call destroy_interpolation_table(pf_o)
      call destroy_interpolation_table(pf_ne)
      call destroy_interpolation_table(pf_mg)
      call destroy_interpolation_table(pf_si)
      call destroy_interpolation_table(pf_fe)

      if (read_init_dat(22, wanted_kh, ksv, kt5, jch) .eqv. .false.) then
         jo3 = -3
         return
      end if
      rewind (22)
      uc(13) = 1.0d3
      uc(14) = 1.0d3
      jo = 0

      ! Restore age, model number and number of models
      age  = age_backup
      jmod = jmod_backup
      jnn  = jnn_backup

      return
   end subroutine make_post_flash_model

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
