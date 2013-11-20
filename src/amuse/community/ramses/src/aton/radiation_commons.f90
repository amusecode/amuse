module radiation_commons
  use amr_commons
  use timing
  implicit none

  ! Photon number density on the AMR grid
  real(dp),allocatable,dimension(:)  ::Erad

  ! Photon emission source field on the AMR grid
  real(dp),allocatable,dimension(:)  ::Srad

  ! Number of cells in each dimension.
  integer::grid_size_x,grid_size_y,grid_size_z

  ! Number of boundary cells on each side of each dimension.
  integer::boundary_size=0

  ! Radiation feedback activated
  logical::radiation_feedback=.true.

  ! Factor used by cooling in CUDATON.
  real(kind=8)::fudgecool=0.1

  ! Speed of light to use for radiation. (We use a lower value than the actual
  ! physical value.) Units: m/s
  real(kind=8)::c_light=2.99792458e8


  ! Radiation namelist parameters:

  logical::allow_gpu_overload=.false.  ! See aton_fortran.h.

  ! ATON version to use. The options are 'gpu' or 'cpu'.
  ! 'cpu' is experimental.
  character(len=10)::rad_aton_version='gpu'

  ! Maximum time step in user units.
  ! -1 means there is no maximum.
  ! This only needs to be set for some test cases where the hydro time
  ! step is too large.
  ! TODO(tstranex): Calculate this automatically.
  real(kind=8)::rad_max_time_step=-1

  real(kind=8)::rad_light_speed_factor=1.0

  real(kind=8)::rad_escape_fraction=0.01

  ! Number of sources.
  integer::rad_num_sources=0
  real(kind=8)::rad_source_x=0.0
  real(kind=8)::rad_source_y=0.0
  real(kind=8)::rad_source_z=0.0
  real(kind=8)::rad_source_rate=0.0  ! [photons / second]

  ! Boundary conditions:
  ! 0 (default) means transmissive (zero-gradient) boundary conditions
  ! 1 means periodic boundary conditions
  integer::rad_boundary_condition=0

  ! Options for imposing a flux at the x min boundary:
  logical::rad_flux_at_x_min_boundary=.false.
  real(kind=8)::rad_flux_x=0.0  ! initial photon flux [photons / s / m^2]
  real(kind=8)::rad_flux_y=0.0
  real(kind=8)::rad_flux_z=0.0
  real(kind=8)::rad_density=0.0  ! inital photon number density [photons / m^3]

  
  ! MPI variables:
  ! TODO(tstranex): Comment on these.

  integer::tag=102

  integer::num_cpu_x=1,num_cpu_y=1,num_cpu_z=1
  integer::my_i,my_j,my_k

  real(dp)::my_xmin,my_ymin,my_zmin,my_xmax,my_ymax,my_zmax

  integer,dimension(:),allocatable::sendbuf,recvbuf,icount
  integer,dimension(:,:,:),allocatable::cpu_rad
  integer::info

  integer,dimension(:),allocatable::receive_request,send_request
  integer::n_receive_request,n_send_request
  integer,dimension(:,:),allocatable::request_status
  integer,dimension(:,:),allocatable::boundary_request_status

  real(kind=8),dimension(:,:),allocatable::boundary_send,boundary_recv
  integer::boundary_memory_size

  type(timer_state)::ramses_timer
  type(timer_state)::total_timer,full_memory_timer,boundary_memory_timer,boundary_timer,aton_timer,mpi_timer

end module radiation_commons
