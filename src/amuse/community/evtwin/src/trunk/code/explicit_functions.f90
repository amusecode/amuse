!> \brief For some variables it is inconvenient to calculate them implicitly, so
!! we calculate them explicitly, between timesteps/iterations
module explicit_functions
   use real_kind
   use mesh

   implicit none
   ! FIXME: these are really redundant with the semi-implicit functions, and should probably be integrated with those
   integer, parameter :: explv_gradmu = 1
   integer, parameter :: explv_avmu   = 2
   integer, parameter :: explv_logp   = 3
   integer, parameter :: num_explv = 4

   real(double), allocatable, save :: expl_var(:,:,:)
   real(double), allocatable, save :: radacc(:,:,:)
end module explicit_functions


