!> 
!! The variables in this module used to be in the common block vbles
!< 

module structure_variables
   use real_kind
   use mesh
   
   implicit none
   !These variables are for the current mesh point (?)

   real(double), allocatable, save :: sx(:,:)  ! sx(:,k) stores a large number of variables for the mesh points k=1:NM
   
end module structure_variables

