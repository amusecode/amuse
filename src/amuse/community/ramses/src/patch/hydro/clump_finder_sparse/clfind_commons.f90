module clfind_commons
  use amr_commons, ONLY: qdp,dp

  integer::nparts,nparts_tot,npeaks,npeaks_tot
  real(dp)::tot_mass
  real(dp)::relevance_threshold=1.5
  real(dp)::density_threshold=0.
  real(dp)::mass_threshold=0.
  logical::merge_unbound=.false.,clinfo=.false.

  ! Linked list for saddlepoint values
  integer::nmatmax=500000
  integer::last_imat_free
  integer,allocatable,dimension(:)::icurrent
  integer,allocatable,dimension(:)::imat_next,imat_prev
  real(dp),allocatable,dimension(:)::saddle_dens,saddle_dens_tot
  
  ! Peak patch properties
  real(dp),allocatable,dimension(:,:)::clump_size_tot,center_of_mass_tot,clump_momentum_tot,clump_force_tot
  real(dp),allocatable,dimension(:,:,:)::second_moments,second_moments_tot,Icl_d_3by3_tot,Icl_3by3_tot
  real(dp),allocatable,dimension(:)::min_dens_tot,av_dens_tot,phi_min_tot
  real(dp),allocatable,dimension(:)::max_dens_tot,e_kin_int_tot,e_bind_tot,e_thermal_tot
  !real(dp),allocatable,dimension(:)::e_kin_int_tot4,e_bind_tot4,e_thermal_tot4
  real(dp),allocatable,dimension(:)::clump_mass_tot,clump_vol_tot,clump_mass_tot4
  real(dp),allocatable,dimension(:,:)::peak_pos_tot,bulk_momentum_tot
  real(dp),allocatable,dimension(:)::saddle_max_tot
  real(dp),allocatable,dimension(:)::relevance_tot
  real(dp),allocatable,dimension(:)::phi_ref, phi_ref_tot
  real(dp),allocatable,dimension(:)::Psurf,Psurf_tot,v_therm_tot,v_rms_tot,m4_tot
  real(dp),allocatable,dimension(:)::e_bind_iso_tot,e_therm_iso_tot,e_kin_iso_tot,grav_term_tot,clump_virial_tot
  real(dp),allocatable,dimension(:)::peak_check,ball4_check,isodens_check,clump_check
  real(dp),allocatable,dimension(:)::Icl_tot,Icl_d_tot,Icl_dd_tot
  logical,allocatable,dimension(:)::contracting

  ! Test particles properties
  real(dp),allocatable,dimension(:)::denp ! Density of the cell containing a test particle. Davide: used by the clump finder.
  integer,allocatable,dimension(:)::iglobalp,icellp,levp,testp_sort ! Used to sort test particles by density  
  integer,allocatable,dimension(:)::n_cells_tot,minmatch_tot,new_peak
  integer,allocatable,dimension(:)::sort_index
  integer,allocatable,dimension(:)::occupied,occupied_all ! Tells whether there is already a sink in a clump.
  integer,allocatable,dimension(:)::form,form_all ! Tells whether a sink has to be formed within a clump.

end module clfind_commons
