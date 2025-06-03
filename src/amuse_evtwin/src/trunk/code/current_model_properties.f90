!> \file current_model_properties.f90  Store properties of the current model
!!
!! This module used to be the common block query and stores:
!! - Properties of the current loop iteration (primary mass, mass ratio, orbital period)
!! - Model numbers (j...)
!!
!! \todo Sort out which model number (j...) is which


module current_model_properties
   use real_kind
   implicit none

   integer :: jmod    ! Current model number
   integer :: jb      ! Single star (1) or binary (2). Compare ktw
   integer :: jnn     ! Number of calculated models. Different from jmod because of back-tracking
   integer :: jter    ! Current iteration number in the solver (FIXME: why is this even a global???)
   integer :: joc     ! Flag to decide whether we're solving for structure (1), or nucleosynthesis (2 or 3, for *1,2)
   integer :: jkh     ! Something about debugging solver with printS/printC; related to jh3 (whatever that does). 
                      ! So no, not sure what it's for

   ! Variable number holding the maximum error and the value of the maximum error for the last iteration
   integer :: var_max_err
   integer :: point_max_err
   real(double) :: max_err

   ! Values used for time-step control
   real(double), save :: rlf_prev(2) = (/0.0, 0.0/)   ! Previous value of Roche-lobe filling factor
   real(double), save :: qcnv_prev(2) = (/0.0, 0.0/)  ! Previous value of mass fraction of convective envelope
   real(double), save :: lnuc_prev(2) = (/0.0, 0.0/)  ! Previous value of nuclear burning luminosity
   real(double), save :: lhe_prev(2) = (/0.0, 0.0/)   ! Previous value of He burning luminosity
   real(double), save :: lh_prev(2) = (/0.0, 0.0/)    ! Previous value of H burning luminosity
 end module current_model_properties

