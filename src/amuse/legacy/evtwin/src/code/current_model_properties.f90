!> current_model_properties:
!!
!! This module used to be the common block query and stores:
!! 
!! * Properties of the current loop iteration (primary mass, mass ratio, orbital period)
!! 
!! * Termination conditions (UC)
!! 
!! * Model numbers (j...)
!!
!! \todo Sort out which model number (j...) is which
!<

module current_model_properties
   
   use real_kind
   
   implicit none
   
   !Loop iteration properties:
   real(double) :: ml       !> Log of the primary mass in the current iteration of the mass loop (from init.run)
   real(double) :: ql       !> Log of the mass ratio in the current iteration of the mass-ratio loop (from init.run)
   real(double) :: xl       !> Log of the orbital period in the current iteration of the Porb loop (from init.run)
   
   !Terminal conditions:
   real(double) :: uc(21)   !> Termination and special-mode (e.g. He-flash evation) conditions (from init.run)
   
   integer :: jmod    !> Current model number
   integer :: jb      !> Single star (1) or binary (2). Compare ktw
   integer :: jnn     !> Number of calculated models. Different from jmod because of back-tracking
   integer :: jter    !> Current iteration number in the solver (FIXME: why is this even a global???)
   integer :: joc     !> Flag to decide whether we're solving for structure (1), or nucleosynthesis (2 or 3, for *1,2)
   integer :: jkh     !> Something about debugging solver with printS/printC; related to jh3 (whatever that does). 
                      !! So no, not sure what it's for

   ! Variable number holding the maximum error and the value of the maximum error for the last iteration
   integer :: var_max_err
   integer :: point_max_err
   real(double) :: max_err
   
end module current_model_properties

