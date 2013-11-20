module pm_commons
  use amr_parameters
  use pm_parameters
  use random
  ! Sink particle related arrays
  real(dp),allocatable,dimension(:)::msink,r2sink,v2sink,c2sink,oksink_new,oksink_all
  real(dp),allocatable,dimension(:)::tsink,tsink_new,tsink_all
  real(dp),allocatable,dimension(:)::msink_new,msink_all,r2k,v2sink_new,c2sink_new
  real(dp),allocatable,dimension(:)::v2sink_all,c2sink_all
  real(dp),allocatable,dimension(:)::xmsink
  real(dp),allocatable,dimension(:)::dMsink_overdt
  real(dp),allocatable,dimension(:)::delta_mass,delta_mass_new,delta_mass_all
  real(dp),allocatable,dimension(:)::vol_gas_agn,mass_gas_agn
  real(dp),allocatable,dimension(:)::mass_blast_agn,vol_blast_agn,p_agn
  real(dp),allocatable,dimension(:)::vol_gas_agn_all,mass_gas_agn_all
  real(dp),allocatable,dimension(:)::wden,weth,wvol,wden_new,weth_new,wvol_new
  real(dp),allocatable,dimension(:,:)::divsink,divsink_new
  real(dp),allocatable,dimension(:)::total_volume
  real(dp),allocatable,dimension(:,:)::wmom,wmom_new
  real(dp),allocatable,dimension(:,:)::vsink,vsink_new,vsink_all
  real(dp),allocatable,dimension(:,:)::fsink,fsink_new,fsink_all
  real(dp),allocatable,dimension(:,:,:)::vsnew,vsold
  real(dp),allocatable,dimension(:,:,:)::fsink_partial,sink_jump
  real(dp),allocatable,dimension(:,:)::lsink,lsink_new,lsink_all !sink angular momentum
  real(dp),allocatable,dimension(:,:)::xsink,xsink_new,xsink_all
  real(dp),allocatable,dimension(:)::acc_rate,acc_lum !sink accretion rate and luminosity
  real(dp),allocatable,dimension(:,:)::weighted_density,weighted_volume,weighted_ethermal
  real(dp),allocatable,dimension(:,:,:)::weighted_momentum
  integer,allocatable,dimension(:)::idsink,idsink_new,idsink_old,idsink_all
  integer,allocatable,dimension(:)::level_sink,level_sink_new,level_sink_all
  logical,allocatable,dimension(:)::ok_blast_agn,ok_blast_agn_all
  integer,allocatable,dimension(:)::idsink_sort,ind_blast_agn,new_born,new_born_all
  integer::ncloud_sink
  integer::nindsink=0

  ! Particles related arrays
  real(dp),allocatable,dimension(:,:)::xp       ! Positions
  real(dp),allocatable,dimension(:,:)::vp       ! Velocities
  real(dp),allocatable,dimension(:)  ::mp       ! Masses
#ifdef OUTPUT_PARTICLE_POTENTIAL
  real(dp),allocatable,dimension(:)  ::ptcl_phi ! Potential of particle added by AP for output purposes 
#endif
  real(dp),allocatable,dimension(:)  ::tp       ! Birth epoch
  real(dp),allocatable,dimension(:)  ::zp       ! Birth metallicity
  integer ,allocatable,dimension(:)  ::nextp    ! Next particle in list
  integer ,allocatable,dimension(:)  ::prevp    ! Previous particle in list
  integer ,allocatable,dimension(:)  ::levelp   ! Current level of particle
  integer ,allocatable,dimension(:)  ::idp      ! Identity of particle
  ! Tree related arrays
  integer ,allocatable,dimension(:)  ::headp    ! Head particle in grid
  integer ,allocatable,dimension(:)  ::tailp    ! Tail particle in grid
  integer ,allocatable,dimension(:)  ::numbp    ! Number of particles in grid
  ! Global particle linked lists
  integer::headp_free,tailp_free,numbp_free=0,numbp_free_tot=0
  ! Local and current seed for random number generator
  integer,dimension(IRandNumSize) :: localseed=-1
end module pm_commons
