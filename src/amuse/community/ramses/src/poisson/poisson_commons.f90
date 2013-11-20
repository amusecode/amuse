module poisson_commons 
  use amr_commons
  use poisson_parameters

  real(dp),allocatable,dimension(:)  ::phi,phi_old       ! Potential
  real(dp),allocatable,dimension(:)  ::rho               ! Density
  real(dp),allocatable,dimension(:,:)::f                 ! 3-force

  real(dp),allocatable,dimension(:)  ::rho_top   ! Density at last CIC level                                 

  ! Multigrid lookup table for amr -> mg index mapping
  integer, allocatable, dimension(:) :: lookup_mg   ! Lookup table

  ! Communicator arrays for multigrid levels
  type(communicator), allocatable, dimension(:,:) :: active_mg
  type(communicator), allocatable, dimension(:,:) :: emission_mg

  ! Minimum MG level
  integer :: levelmin_mg

  ! Multigrid safety switch
  logical, allocatable, dimension(:) :: safe_mode

  ! Multipole coefficients
  real(dp),dimension(1:ndim+1)::multipole

end module poisson_commons
