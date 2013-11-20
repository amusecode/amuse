module amr_commons
  use amr_parameters

  logical::output_done=.false.                  ! Output just performed
  logical::init=.false.                         ! Set up or run
  logical::balance=.false.                      ! Load balance or run
  logical::shrink=.false.                       ! Shrink mesh or run
  integer::nstep=0                              ! Time step
  integer::nstep_coarse=0                       ! Coarse step
  integer::nstep_coarse_old=0                   ! Old coarse step
  integer::nflag,ncreate,nkill                  ! Refinements
  integer::ncoarse                              ! nx.ny.nz
  integer::ngrid_current                        ! Actual number of octs

  real(dp)::emag_tot=0.0D0                      ! Total magnetic energy
  real(dp)::ekin_tot=0.0D0                      ! Total kinetic energy
  real(dp)::eint_tot=0.0D0                      ! Total internal energy
  real(dp)::epot_tot=0.0D0                      ! Total potential energy
  real(dp)::epot_tot_old=0.0D0                  ! Old potential energy
  real(dp)::epot_tot_int=0.0D0                  ! Time integrated potential
  real(dp)::const=0.0D0                         ! Energy conservation
  real(dp)::aexp_old=1.0D0                      ! Old expansion factor
  real(dp)::rho_tot=0.0D0                       ! Mean density in the box
  real(dp)::t=0.0D0                             ! Time variable

  ! MPI variables
  integer::ncpu,ndomain,myid,overload=1

  ! Friedman model variables
  integer::n_frw
  real(dp),allocatable,dimension(:)::aexp_frw,hexp_frw,tau_frw,t_frw

  ! Initial conditions parameters from grafic
  integer                  ::nlevelmax_part
  real(dp)                 ::aexp_ini=10.
  real(dp),dimension(1:MAXLEVEL)::dfact=1.0d0,astart
  real(dp),dimension(1:MAXLEVEL)::vfact
  real(dp),dimension(1:MAXLEVEL)::xoff1,xoff2,xoff3,dxini
  integer ,dimension(1:MAXLEVEL)::n1,n2,n3

  ! Level related arrays
  real(dp),dimension(1:MAXLEVEL)::dtold,dtnew ! Time step at each level
  real(dp),dimension(1:MAXLEVEL)::rho_max     ! Maximum density at each level
  integer ,dimension(1:MAXLEVEL)::nsubcycle=2 ! Subcycling at each level

  ! Pointers for each level linked list
  integer,allocatable,dimension(:,:)::headl
  integer,allocatable,dimension(:,:)::taill
  integer,allocatable,dimension(:,:)::numbl
  integer,allocatable,dimension(:,:)::numbtot

  ! Pointers for each level boundary linked list
  integer,allocatable,dimension(:,:)::headb
  integer,allocatable,dimension(:,:)::tailb
  integer,allocatable,dimension(:,:)::numbb

  ! Pointers for free memory grid linked list
  integer::headf,tailf,numbf,used_mem,used_mem_tot

  ! Tree arrays
  real(dp),allocatable,dimension(:,:)::xg      ! grids position
  integer ,allocatable,dimension(:,:)::nbor    ! neighboring father cells
  integer ,allocatable,dimension(:)  ::father  ! father cell
  integer ,allocatable,dimension(:)  ::next    ! next grid in list
  integer ,allocatable,dimension(:)  ::prev    ! previous grid in list
  integer ,allocatable,dimension(:)  ::son     ! sons grids
  integer ,allocatable,dimension(:)  ::flag1   ! flag for refine
  integer ,allocatable,dimension(:)  ::flag2   ! flag for expansion

  ! Global indexing
  integer ,allocatable,dimension(:)  ::cpu_map  ! domain decomposition
  integer ,allocatable,dimension(:)  ::cpu_map2 ! new domain decomposition

  ! Hilbert key
  real(qdp),allocatable,dimension(:)::hilbert_key
  real(qdp),allocatable,dimension(:)::bound_key,bound_key2
  real(qdp)                         ::order_all_min,order_all_max

  ! Recursive bisection                                                                               
  real(dp),allocatable,dimension(:)    ::bisec_wall         ! bisection wall positions                
  integer ,allocatable,dimension(:,:)  ::bisec_next         ! next 2 child cells in bisection         
  integer::bisec_root                                       ! root of bisection tree                  

  integer,allocatable,dimension(:)     ::bisec_indx         ! map from leaf cell id to cpu id         
  real(dp),allocatable,dimension(:,:)  ::bisec_cpubox_min   ! cpu domains boxes                       
  real(dp),allocatable,dimension(:,:)  ::bisec_cpubox_max
  real(dp),allocatable,dimension(:,:)  ::bisec_cpubox_min2  ! cpu domains boxes for new decomp        
  real(dp),allocatable,dimension(:,:)  ::bisec_cpubox_max2

  integer,allocatable,dimension(:)     ::bisec_cpu_load     ! CPU loads (for stats)                   
  integer,allocatable,dimension(:,:)   ::bisec_hist         ! histograms for load computation         
  integer,allocatable,dimension(:)     ::bisec_hist_bounds  ! histogram splitting boundaries          
  integer,allocatable,dimension(:)     ::new_hist_bounds
  integer,allocatable,dimension(:)     ::bisec_ind_cell     ! histo swap id -> cell id map (big)      
  integer,allocatable,dimension(:)     ::cell_level         ! store the level of the cells (big)      

  real(dp)::bisec_res                                       ! resolution parameters                   
  integer ::bisec_nres

  ! Communication structure
  type communicator
     integer                            ::ngrid
     integer                            ::npart
     integer     ,dimension(:)  ,pointer::igrid
     integer     ,dimension(:,:),pointer::f
     real(kind=8),dimension(:,:),pointer::u
     integer     ,dimension(:,:),pointer::fp
     real(kind=8),dimension(:,:),pointer::up
#ifdef ATON
     real(kind=8),dimension(:,:),pointer::u_radiation
#endif
  end type communicator

  ! Active grid, emission and reception communicators
  type(communicator),allocatable,dimension(:)  ::active
  type(communicator),allocatable,dimension(:,:)::boundary
  type(communicator),allocatable,dimension(:,:)::emission
  type(communicator),allocatable,dimension(:,:)::reception

  ! Types for physical boundary conditions
  CHARACTER(LEN=20)::type_hydro  ='hydro'
  CHARACTER(LEN=20)::type_accel  ='accel'
  CHARACTER(LEN=20)::type_flag   ='flag'

  ! Units specified by the user in the UNITS_PARAMS namelist for non-cosmo runs.
  ! These values shouldn't be used directly. Instead call units() in amr/units.f90.
  real(dp)::units_density=1.0  ! [g/cm^3]
  real(dp)::units_time=1.0     ! [seconds]
  real(dp)::units_length=1.0   ! [cm]

end module amr_commons

