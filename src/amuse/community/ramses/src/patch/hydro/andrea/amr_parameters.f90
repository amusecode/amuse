module amr_parameters

  ! Define real types
  integer,parameter::sp=kind(1.0E0)
#ifndef NPRE
  integer,parameter::dp=kind(1.0E0) ! default
#else
#if NPRE==4
  integer,parameter::dp=kind(1.0E0) ! real*4
#else
  integer,parameter::dp=kind(1.0D0) ! real*8
#endif
#endif
  
  ! Number of dimensions
#ifndef NDIM
  integer,parameter::ndim=1
#else
  integer,parameter::ndim=NDIM
#endif
  integer,parameter::twotondim=2**ndim
  integer,parameter::threetondim=3**ndim
  integer,parameter::twondim=2*ndim

  ! Vectorization parameter
#ifndef NVECTOR
  integer,parameter::nvector=500  ! Size of vector sweeps
#else
  integer,parameter::nvector=NVECTOR
#endif
  ! Run control
  logical::verbose =.false.   ! Write everything
  logical::hydro   =.false.   ! Hydro activated
  logical::pic     =.false.   ! Particle In Cell activated
  logical::poisson =.false.   ! Poisson solver activated
  logical::cosmo   =.false.   ! Cosmology activated
  logical::star    =.false.   ! Star formation activated
  logical::debug   =.false.   ! Debug mode activated
  logical::static  =.false.   ! Static mode activated
  logical::tracer  =.false.   ! Tracer particles activated

  ! Mesh parameters
  integer::geom=1             ! 1: cartesian, 2: cylindrical, 3: spherical
  integer::nx=1,ny=1,nz=1     ! Number of coarse cells in each dimension
  integer::levelmin=1         ! Full refinement up to levelmin
  integer::nlevelmax=1        ! Maximum number of level
  integer::ngridmax=0         ! Maximum number of grids
  integer::nexpand=1          ! Number of mesh expansion
  integer::nexpand_bound=1    ! Number of mesh expansion for virtual boundaries
  real(dp)::boxlen=1.0D0      ! Box length along x direction
  character(len=128)::ordering='hilbert'

  ! Step parameters
  integer::nrestart=0         ! New run or backup file number
  integer::nstepmax=1000000   ! Maximum number of time steps
  integer::ncontrol=1         ! Write control variables
  integer::fbackup=1000000    ! Backup data to disk
  integer::nremap=0           ! Load balancing frequency (0: never)

  ! Output parameters
  integer::iout=1             ! Increment for output times
  integer::ifout=1            ! Increment for output files
  integer::iback=1            ! Increment for backup files
  integer::noutput=1          ! Total number of outputs
  integer::foutput=1000000    ! Frequency of outputs
  integer::output_mode=0      ! Output mode (for hires runs)

  ! Cosmology
  real(dp)::boxlen_ini        ! Box size in h-1 Mpc
  real(dp)::omega_b=0.0D0     ! Omega Baryon
  real(dp)::omega_m=1.0D0     ! Omega Matter
  real(dp)::omega_l=0.0D0     ! Omega Lambda
  real(dp)::omega_k=0.0D0     ! Omega Curvature
  real(dp)::h0     =1.0D0     ! Hubble constant in km/s/Mpc
  real(dp)::aexp   =1.0D0     ! Current expansion factor
  real(dp)::hexp   =0.0D0     ! Current Hubble parameter
  real(dp)::n_star =0.1D0     ! Star formation density threshold in H/cc
  real(dp)::t_star =0.0D0     ! Star formation time scale in Gyr
  real(dp)::T2_star=0.0D0     ! Typical ISM temperature
  real(dp)::g_star =1.6D0     ! Typical ISM polytropic index
  real(dp)::del_star=2.D2     ! Minimum overdensity to define ISM
  real(dp)::eta_sn =0.0D0     ! Supernova mass fraction
  real(dp)::yield  =0.0D0     ! Supernova yield
  real(dp)::f_w    =1.0D1     ! Wind loading factor
  real(dp)::z_ave  =0.0D0     ! Average metal abundance
  real(dp)::T2_start          ! Starting gas temperature
  real(dp)::t_delay           ! Feedback time delay
  integer ::ndebris=100       ! Number of debris particles
  real(dp)::J21    =0.0D0     ! UV flux at threshold in 10^21 units
  real(dp)::a_spec =1.0D0     ! Slope of the UV spectrum
  logical ::pressure_fix=.false.
  logical ::cooling=.false.
  logical ::metal=.false.
  logical ::haardt_madau=.false.

  ! Halo parameters
  real(dp)::Mvir=1.0D0        ! Halo mass in solar units
  real(dp)::overdensity=1.0d0 ! Overdensity=101 for z=0.
  real(dp)::Mdot_BH=0.0d0     ! Black hole accretion rate

  ! Jet parameters
  real(dp)::t_jet=4D1         ! Jet life-time in Myr
  real(dp)::r_jet=5D0         ! Jet radius in fine cell units
  real(dp)::h_jet=4D0         ! Jet height in fine cell units
  real(dp)::e_jet=1D2         ! Mass loading factor

  ! Output times
  integer,parameter::MAXOUT=1000
  real(dp),dimension(1:MAXOUT)::aout=1.1       ! Output expansion factors
  real(dp),dimension(1:MAXOUT)::tout=0.0       ! Output times

  ! Refinement parameters for each level
  integer,parameter::MAXLEVEL=100
  real(dp),dimension(1:MAXLEVEL)::m_refine   =-1.0 ! Lagrangian threshold
  real(dp),dimension(1:MAXLEVEL)::r_refine   =-1.0 ! Radius of refinement region
  real(dp),dimension(1:MAXLEVEL)::x_refine   = 0.0 ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::y_refine   = 0.0 ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::z_refine   = 0.0 ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::exp_refine = 2.0 ! Exponent for distance
  real(dp),dimension(1:MAXLEVEL)::a_refine   = 1.0 ! Ellipticity (Y/X)
  real(dp),dimension(1:MAXLEVEL)::b_refine   = 1.0 ! Ellipticity (Z/X)

  ! Initial condition files for each level
  logical::multiple=.false.
  character(LEN=80),dimension(1:MAXLEVEL)::initfile=' '
  character(LEN=20)::filetype='ascii'

  ! Initial condition regions parameters
  integer,parameter::MAXREGION=100
  integer                           ::nregion=0
  character(LEN=10),dimension(1:MAXREGION)::region_type='square'
  real(dp),dimension(1:MAXREGION)   ::x_center=0.
  real(dp),dimension(1:MAXREGION)   ::y_center=0.
  real(dp),dimension(1:MAXREGION)   ::z_center=0.
  real(dp),dimension(1:MAXREGION)   ::length_x=1.E10
  real(dp),dimension(1:MAXREGION)   ::length_y=1.E10
  real(dp),dimension(1:MAXREGION)   ::length_z=1.E10
  real(dp),dimension(1:MAXREGION)   ::exp_region=2.0

  ! Boundary conditions parameters
  integer,parameter::MAXBOUND=100
  logical                           ::simple_boundary=.false.
  integer                           ::nboundary=0
  integer                           ::icoarse_min=0
  integer                           ::icoarse_max=0
  integer                           ::jcoarse_min=0
  integer                           ::jcoarse_max=0
  integer                           ::kcoarse_min=0
  integer                           ::kcoarse_max=0
  integer ,dimension(1:MAXBOUND)    ::boundary_type=0
  integer ,dimension(1:MAXBOUND)    ::ibound_min=0
  integer ,dimension(1:MAXBOUND)    ::ibound_max=0
  integer ,dimension(1:MAXBOUND)    ::jbound_min=0
  integer ,dimension(1:MAXBOUND)    ::jbound_max=0
  integer ,dimension(1:MAXBOUND)    ::kbound_min=0
  integer ,dimension(1:MAXBOUND)    ::kbound_max=0

end module amr_parameters
