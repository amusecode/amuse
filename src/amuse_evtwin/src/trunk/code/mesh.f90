module mesh
  use real_kind
  use indices, only: nvar, nfunc, neq

  implicit none
  ! Define the maximum number of meshpoints
  !integer, parameter :: NM = 1000
  integer, save :: max_nm = 200
  integer :: NM

  ! Number of variables in the PX and SX arrays (output variables)
  integer, parameter :: npx = 90

  ! Various arrays for storing the structure of the stars:
  !   H(nvar, nm) - Array of stellar structure variables
  !  DH(nvar, nm) - Changes since last timestep
  ! HPR(nvar, nm) - the previous converged model
  real(double), allocatable, save :: h(:,:)
  real(double), allocatable, save :: dh(:,:)
  real(double), allocatable, save :: hpr(:,:)

  ! Various control parameters:
  ! KH    - Total number of grid points in use
  ! KTW   - TWIN mode (2) or not (1)
  ! ISB   - Single star (1) or binary (2)
  ! ID(:) - List of variables and equations to solve (input for solver)
  ! IE(:) - Ditto, for nucleosynthesis part of the run
  integer, save :: kh, ktw, isb, id(130), ie(130)
end module mesh

