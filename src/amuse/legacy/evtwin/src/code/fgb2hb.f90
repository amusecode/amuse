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
   use fgb2hb_composition
   use test_variables
   use current_model_properties
   use constants
   
   implicit none
   integer :: jop, jo3
   integer :: cmi_mode_bck,jc1,ksv
   real(double) :: bm_fgb, bper_fgb
   
   if(debug.ge.2) write(6,'(/,A20,2I12)')'fgb2hb: ',jop,jo3
   
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

   call store_pre_flash_composition
   call star12 ( jo3, jc1, jop, 12, ksv, 24 )
   call cleanup_pre_flash_composition

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

