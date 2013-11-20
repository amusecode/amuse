subroutine compute_clump_properties(ntest)
  use amr_commons
  use hydro_commons, ONLY:uold
  use clfind_commons
  use poisson_commons, ONLY:phi,f
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ntest

  !----------------------------------------------------------------------------
  ! this subroutine performs one loop over all "test-particles" and collects the 
  ! relevant information from the cells where the particles sit. After a lot
  ! of mpi-communication, all necessary peak-patch properties can be computed
  !----------------------------------------------------------------------------


  integer::ipart,grid,info,i,j,peak_nr,ilevel
  real(dp)::zero=0.
  !variables needed temporarily store cell properties
  real(dp)::d,vol
  real(dp),dimension(1:3)::vd
  ! variables related to the size of a cell on a given level
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp),dimension(1:nlevelmax)::volume
  real(dp),dimension(1:3)::skip_loc,xcell
  real(dp),dimension(1:twotondim,1:3)::xc
  integer::nx_loc,ind,ix,iy,iz
  !peak-patch related arrays before sharing information with other cpus
  real(dp),dimension(1:npeaks_tot)::max_dens,phi_min
  real(dp),dimension(1:npeaks_tot)::min_dens,clump_mass,clump_vol
  real(dp),dimension(1:npeaks_tot,1:3)::center_of_mass,clump_momentum,peak_pos,clump_force
  integer,dimension(1:npeaks_tot)::n_cells


  min_dens=huge(zero);  max_dens=0.d0; n_cells=0; phi_min=huge(zero); second_moments=0.
  clump_mass=0.d0; clump_vol=0.d0; clump_momentum=0.d0; center_of_mass=0.d0; clump_force=0.d0 
  peak_pos=0.

  !------------------------------------------
  ! compute volume of a cell in a given level
  !------------------------------------------
  do ilevel=1,nlevelmax
     ! Mesh spacing in that level
     dx=0.5D0**ilevel 
     nx_loc=(icoarse_max-icoarse_min+1)
     scale=boxlen/dble(nx_loc)
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     volume(ilevel)=vol_loc
  end do

  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)
     xc(ind,2)=(dble(iy)-0.5D0)
     xc(ind,3)=(dble(iz)-0.5D0)
  end do


  !---------------------------------------------------------------------------
  ! loop over all test particles to collect information from the cells
  !---------------------------------------------------------------------------
  do ipart=1,ntest     
     ! peak number after merge
     peak_nr=flag2(icellp(ipart)) 

     if (peak_nr /=0 ) then
        
        ! Cell coordinates
        ind=(icellp(ipart)-ncoarse-1)/ngridmax+1 ! cell position
        grid=icellp(ipart)-ncoarse-(ind-1)*ngridmax ! grid index
        dx=0.5D0**levp(ipart)
        xcell(1:ndim)=(xg(grid,1:ndim)+xc(ind,1:ndim)*dx-skip_loc(1:ndim))*scale
        
        ! gas density and momentum density
        d=uold(icellp(ipart),1)
        do i=1,ndim
           vd(i)=uold(icellp(ipart),i+1)
        end do

        ! Cell volume
        vol=volume(levp(ipart))

        ! number of leaf cells per clump
        n_cells(peak_nr)=n_cells(peak_nr)+1
        
        ! find min density and potential
        min_dens(peak_nr)=min(d,min_dens(peak_nr))
        phi_min(peak_nr)=min(phi(icellp(ipart)),phi_min(peak_nr))

        ! find max density and peak location
        if(d>=max_dens(peak_nr))then
           max_dens(peak_nr)=d
           peak_pos(peak_nr,1:ndim)=xcell(1:ndim)
        end if

        ! find clump mass
        clump_mass(peak_nr)=clump_mass(peak_nr)+vol*d

        ! clump volume
        clump_vol(peak_nr)=clump_vol(peak_nr)+vol

        do i=1,ndim

           !ceter of mass location
           center_of_mass(peak_nr,i)=center_of_mass(peak_nr,i)+xcell(i)*vol*d

           !center of mass velocity
           clump_momentum(peak_nr,i)=clump_momentum(peak_nr,i)+vd(i)*vol

           ! center of mass force
           clump_force(peak_nr,i)=clump_force(peak_nr,i)+f(icellp(ipart),i)*vol*d
           
           ! compute second order moments
           do j=1,ndim
              second_moments(peak_nr,i,j)=second_moments(peak_nr,i,j)+xcell(i)*xcell(j)*vol*d
           end do
        end do
     end if
  end do

  ! a lot of MPI communication to collect the results from the different cpu's
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(n_cells,n_cells_tot,npeaks_tot,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(min_dens,min_dens_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(phi_min,phi_min_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(max_dens,max_dens_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(clump_mass,clump_mass_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(clump_force,clump_force_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(clump_vol,clump_vol_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(clump_momentum,clump_momentum_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(center_of_mass,center_of_mass_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(second_moments,second_moments_tot,9*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI      
  n_cells_tot=n_cells
  min_dens_tot=min_dens
  phi_min_tot=phi_min
  max_dens_tot=max_dens
  clump_mass_tot=clump_mass
  clump_vol_tot=clump_vol
  clump_momentum_tot=clump_momentum
  clump_force_tot=clump_force
  center_of_mass_tot=center_of_mass
  second_moments_tot=second_moments
#endif
  !clean some wannabe peaks (due to MPI)
  do j=1,npeaks_tot
     do i=1,ndim
        if (max_dens(j)<max_dens_tot(j))peak_pos(j,i)=0.d0
     end do
  end do
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(peak_pos,peak_pos_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  peak_pos_tot=peak_pos
#endif
  
  !calculate total mass above threshold
  tot_mass=sum(clump_mass_tot)
  ! compute further properties of the clumps
  av_dens_tot(1:npeaks_tot)=clump_mass_tot(1:npeaks_tot)/(clump_vol_tot(1:npeaks_tot)+tiny(0.d0))
  do i=1,ndim
     center_of_mass_tot(1:npeaks_tot,i)=center_of_mass_tot(1:npeaks_tot,i)/(clump_mass_tot(1:npeaks_tot)+tiny(0.d0))
     clump_force_tot(1:npeaks_tot,i)=clump_force_tot(1:npeaks_tot,i)/(clump_mass_tot(1:npeaks_tot)+tiny(0.d0))
  end do

end subroutine compute_clump_properties
!################################################################
!################################################################
!################################################################
!################################################################
subroutine compute_clump_properties_round2(ntest,all_bound)
  use amr_commons
  use hydro_commons, ONLY:uold,gamma
  use poisson_commons, ONLY:phi,f
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::all_bound
  integer::ntest

  !----------------------------------------------------------------------------
  ! This subroutine performs another loop over all particles and collects 
  ! more information like binding energies, etc, that can not be created by
  ! just summing up cell properties.
  !----------------------------------------------------------------------------

  integer::ipart,ilevel,info,i,peak_nr,j,ii,jj
  integer::grid,nx_loc,ix,iy,iz,ind
  real(dp)::d,vol,M,ekk,phi_rel,de,c_sound,d0,v_bulk2,p,v_rms
  real(dp)::t_larson1,cont_speed=0.
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,cty
  real(dp)::dx,dx_loc,scale,vol_loc,abs_err,A1=0.,A2=0.,A3=0.
  real(dp),dimension(1:nlevelmax)::volume
  real(dp),dimension(1:3)::vd,xcell,xpeak,v_cl,rrel,vrel,frel,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:npeaks_tot)::e_kin_int,e_bind,e_thermal,v_therm,m4
  real(dp),dimension(1:npeaks_tot)::clump_mass4,e_kin_iso,e_bind_iso,e_therm_iso,grav_term,Icl_d,Icl
  real(dp),dimension(1:npeaks_tot,1:3)::clump_size,bulk_momentum,contractions
  real(dp),dimension(1:npeaks_tot,1:3,1:3)::Icl_d_3by3,Icl_3by3
  real(dp),dimension(1:3,1:3)::eigenv,a

  !  first, get minimum potential on saddle surface
  call get_phi_ref(ntest)
  
  !initialize arrays
  e_kin_int=0.d0; clump_size=0.d0; e_bind=0.d0; e_thermal=0.d0
  clump_mass4=0.d0
  v_therm=0.; bulk_momentum=0.; m4=0.
  e_kin_iso=0.; e_bind_iso=0.; e_therm_iso=0.
  grav_term=0.d0; Icl_d=0.d0; Icl=0.; Icl_dd_tot=0.
  Icl_3by3=0.;  Icl_d_3by3=0.
  contracting=.false.

  ! Conversion factor from user units to cgs units                                             
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  d0 = density_threshold/scale_nH;
  if(cosmo)d0=d0/aexp**3
  !lifetime of first larson core in code units                                                  
  cty=scale_t/(365.25*24.*3600.)
  t_larson1=larson_lifetime/cty
  if (merge_stars)cont_speed=-1./t_larson1



  !------------------------------------------
  ! compute volume of a cell in a given level
  !------------------------------------------
  do ilevel=1,nlevelmax
     ! Mesh spacing in that level
     dx=0.5D0**ilevel 
     nx_loc=(icoarse_max-icoarse_min+1)
     scale=boxlen/dble(nx_loc)
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     volume(ilevel)=vol_loc
  end do

  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
 
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)
     xc(ind,2)=(dble(iy)-0.5D0)
     xc(ind,3)=(dble(iz)-0.5D0)
  end do

  !---------------------------------------------------------------------------
  ! loop over all test particles to collect information from the cells
  !---------------------------------------------------------------------------
  do ipart=1,ntest     
     ! peak number after merge
     peak_nr=flag2(icellp(ipart)) 

     if (peak_nr /=0 ) then
        
        ! Cell coordinates
        ind=(icellp(ipart)-ncoarse-1)/ngridmax+1 ! cell position
        grid=icellp(ipart)-ncoarse-(ind-1)*ngridmax ! grid index
        dx=0.5D0**levp(ipart)
        xcell(1:ndim)=(xg(grid,1:ndim)+xc(ind,1:ndim)*dx-skip_loc(1:ndim))*scale

        ! gas density and energydensity
        d=uold(icellp(ipart),1)
        de=uold(icellp(ipart),ndim+2)
        do i=1,ndim
           vd(i)=uold(icellp(ipart),i+1)
           xpeak(i)=peak_pos_tot(peak_nr,i)
        end do
        
        vol=volume(levp(ipart))                  


        M=clump_mass_tot(peak_nr)
        v_cl(1:ndim)=clump_momentum_tot(peak_nr,1:ndim)/M
        
        !properties of the cell relative to center of mass
        rrel=xcell(1:3)-center_of_mass_tot(peak_nr,1:3)
        vrel=vd(1:3)/d-v_cl(1:3)
        frel=f(icellp(ipart),1:3)-clump_force_tot(peak_nr,1:3)

        

        do i=1,ndim
           ! size relative to center of mass
           clump_size(peak_nr,i)=clump_size(peak_nr,i)+rrel(i)**2 * vol

           ! internal kinetic energy
           e_kin_int(peak_nr)=e_kin_int(peak_nr)+vrel(i)**2*d*vol*0.5
        end do

        ! potential energy using the acutal phi W= 0.5*int phi_rel*rho
        phi_rel=(phi(icellp(ipart))-phi_ref_tot(peak_nr))*scale
        e_bind(peak_nr)=e_bind(peak_nr)-phi_rel*d*vol*5.d-1
                
        ! thermal energy
        ekk=0.
        do i=1,3 
           ekk=ekk+0.5*vd(i)**2/d                          
        end do
        p=(de-ekk)*(gamma-1)
        e_thermal(peak_nr)=e_thermal(peak_nr)+1.5*vol*p

        ! sound speed
        c_sound=(de-ekk)/d*gamma/(gamma-1)

        !Mass weighted thermal Velocity
        v_therm(peak_nr)=v_therm(peak_nr)+c_sound*d*vol/M        
                
        !properties for regions close to peak (4 cells away)
        if (((xpeak(1)-xcell(1))**2.+(xpeak(2)-xcell(2))**2.+(xpeak(3)-xcell(3))**2.) .LE. 16.*volume(nlevelmax)**(2./3.))then
           do i=1,3
              bulk_momentum(peak_nr,i)=bulk_momentum(peak_nr,i)+(vd(i)/d-v_cl(i))*vol*(d-d0)
           end do
           m4(peak_nr)=m4(peak_nr)+(d-d0)*vol
           clump_mass4(peak_nr)=clump_mass4(peak_nr)+d*vol           
        end if

        !properties for region enclosed by isopotential surface 
        if (phi_rel<0.)then
           do i=1,3
              !not strictly correct since v_cl is av. vel of WHOLE clump
              e_kin_iso(peak_nr)=e_kin_iso(peak_nr)+(vd(i)/d-v_cl(i))**2*d*vol*0.5
           end do
           e_bind_iso(peak_nr)=e_bind_iso(peak_nr)-phi_rel*d*vol*0.5
           e_therm_iso(peak_nr)=e_therm_iso(peak_nr)+1.5*p*vol
        endif

        !terms for virial theorem analysis
        do i=1,3
           grav_term(peak_nr) = grav_term(peak_nr) + frel(i) * rrel(i) * vol*d
           Icl_d(peak_nr)     = Icl_d(peak_nr)     + vrel(i) * rrel(i) * vol*d
           Icl(peak_nr)       = Icl(peak_nr)       + rrel(i) * rrel(i) * vol*d
           do j=1,3
              Icl_d_3by3(peak_nr,i,j)=  Icl_d_3by3(peak_nr,i,j)   + ( vrel(j) * rrel(i)  +  vrel(i) * rrel(j) )   * vol*d
              Icl_3by3(peak_nr,i,j)  =  Icl_3by3(peak_nr,i,j)     +   rrel(j) * rrel(i)                           * vol*d
           end do
        end do
     end if
  end do
  !---------------------------------------------------------------------------
  ! a lot of MPI communication to collect the results from the different cpu's
  !---------------------------------------------------------------------------
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(e_kin_int,e_kin_int_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(e_bind,e_bind_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(clump_size,clump_size_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(m4,m4_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(bulk_momentum,bulk_momentum_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(e_thermal,e_thermal_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(v_therm,v_therm_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(e_therm_iso,e_therm_iso_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(e_kin_iso,e_kin_iso_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(e_bind_iso,e_bind_iso_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(grav_term,grav_term_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(Icl_d,Icl_d_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(Icl,Icl_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(Icl_d_3by3,Icl_d_3by3_tot,9*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(Icl_3by3,Icl_3by3_tot,9*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(clump_mass4,clump_mass_tot4,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  e_kin_int_tot=e_kin_int
  e_bind_tot=e_bind
  clump_size_tot=clump_size
  m4_tot=m4
  bulk_momentum_tot=bulk_momentum
  e_thermal_tot=e_thermal
  v_therm_tot=v_therm
  e_bind_iso_tot=e_bind_iso
  e_kin_iso_tot=e_kin_iso
  e_therm_iso_tot=e_therm_iso
  grav_term_tot=grav_term
  Icl_d_tot=Icl_d
  Icl_tot=Icl
  Icl_d_3by3_tot=Icl_d_3by3
  Icl_3by3_tot=Icl_3by3
  clump_mass_tot4=clump_mass4
#endif

  !second time derivative of I
  Icl_dd_tot(1:npeaks_tot)=2.*(grav_term_tot(1:npeaks_tot)-Psurf_tot(1:npeaks_tot)+2*e_kin_int_tot(1:npeaks_tot)+2*e_thermal_tot(1:npeaks_tot))

  all_bound=.true.

  do j=npeaks_tot,1,-1
     if (relevance_tot(j)>0.)then
        !compute eigenvalues and eigenvectors of Icl_d_3by3_tot
        a=Icl_3by3_tot(j,1:3,1:3)
        abs_err=1.d-8*Icl_tot(j)**2+1.d-40
        call jacobi(a,eigenv,abs_err)
        A1=a(1,1); A2=a(2,2); A3=a(3,3)

        !compute the contractions along the eigenvectors of Icl
        contractions(j,1:3)=0.
        do ii=1,3
           do jj=1,3
              contractions(j,1)=contractions(j,1)+Icl_d_3by3_tot(j,ii,jj)*eigenv(1,ii)*eigenv(1,jj)
              contractions(j,2)=contractions(j,2)+Icl_d_3by3_tot(j,ii,jj)*eigenv(2,ii)*eigenv(2,jj)
              contractions(j,3)=contractions(j,3)+Icl_d_3by3_tot(j,ii,jj)*eigenv(3,ii)*eigenv(3,jj)
           end do
        end do

        !Check wether clump is contracting fast enough along all axis
        if (Icl_tot(j)>0)then 
           contracting(j)=.true.
           contracting(j)=contracting(j) .and. contractions(j,1)/(A1+tiny(0.d0)) < cont_speed 
           contracting(j)=contracting(j) .and. contractions(j,2)/(A2+tiny(0.d0)) < cont_speed 
           contracting(j)=contracting(j) .and. contractions(j,3)/(A3+tiny(0.d0)) < cont_speed 
        end if
        
        !compute peak check for smbh sink formation
        v_rms=2.*e_kin_int_tot(j)/clump_mass_tot(j)
        v_bulk2=(bulk_momentum_tot(j,1)**2+bulk_momentum_tot(j,2)**2&
             +bulk_momentum_tot(j,3)**2)/(m4_tot(j)**2+tiny(0.d0))     
        peak_check(j)=scale*(phi_ref_tot(j)-phi_min_tot(j))/((v_therm_tot(j)**2+v_rms+v_bulk2)*0.5+tiny(0.d0))

        !compute other checks (currently not needed for sink formation)
        isodens_check(j)=scale*e_bind_iso_tot(j)/(tiny(0.d0)+2*e_kin_iso_tot(j)+2*e_therm_iso_tot(j))
        clump_check(j)=(-1.*grav_term_tot(j)+Psurf_tot(j))/(tiny(0.d0)+2*e_kin_int_tot(j)+2*e_thermal_tot(j))
        
        !update the all_bound property
        all_bound=all_bound.and.(isodens_check(j)>1.)

     endif
  end do

  !write to the log file some information that could be of interest for debugging etc.
  if(myid==1 .and. clinfo .and. .not. smbh .and. sink)then 
     write(*,'(135A)')'==========================================================================================='
     write(*,'(135A)')'Cl_N     t1[y]      t2[y]      t3[y] |I_d|/I_dd[y] tidal_Fg   Psurf      e_kin      e_therm'
     write(*,'(135A)')'==========================================================================================='
     do j=npeaks_tot,1,-1
        if (relevance_tot(j)>0.)then
           write(*,'(I4,2X,8(E8.2E2,3X))'),j&
                ,A1/(contractions(j,1)+tiny(0.d0))*cty,A2/(contractions(j,2)+tiny(0.d0))*cty,A3/(contractions(j,3)+tiny(0.d0))*cty&
                ,abs(Icl_d_tot(j))/Icl_dd_tot(j)*cty&
                ,grav_term_tot(j),-1.*Psurf_tot(j)&
                ,e_kin_int_tot(j),e_thermal_tot(j)
        end if
     end do
     write(*,'(135A)')'==========================================================================================='
  end if
     
end subroutine compute_clump_properties_round2
!################################################################
!################################################################
!################################################################
!################################################################
subroutine write_clump_properties(to_file)
  use amr_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::to_file

  !---------------------------------------------------------------------------
  ! this routine writes the clump properties to screen and to file
  !---------------------------------------------------------------------------

  integer::i,j,jj,ilun,n_rel,info,nx_loc
  real(dp)::rel_mass,scale
  character(LEN=5)::nchar
  real(dp),dimension(1:npeaks_tot)::peakd
  integer,dimension(1:npeaks_tot)::ind_sort
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  
  !sort clumps by peak density in ascending order
  do i=1,npeaks_tot
     peakd(i)=max_dens_tot(i)
     ind_sort(i)=i
  end do
  call quick_sort_dp(peakd,ind_sort,npeaks_tot) 

  If(to_file)then
     ilun=20
  else 
     ilun=6
  end if

  !print results in descending order to screen/file
  if(myid==1)then
     rel_mass=0.
     n_rel=0
     if (to_file .eqv. .true.) then
        call title(ifout-1,nchar)
        open(unit=20,file=TRIM('output_'//TRIM(nchar)//'/clump_info.txt'),form='formatted')
        open(unit=21,file=TRIM('output_'//TRIM(nchar)//'/clump_masses.txt'),form='formatted')
     end if
     if(smbh)then
        if(verbose)write(ilun,'(135A)')'Cl_N #leaf-cells  peak_x [uu] peak_y [uu] peak_z [uu] size_x [cm] size_y [cm] size_z [cc]'//&
             ' |v|_CM [u.u.] rho- [H/cc] rho+ [H/cc] rho_av [H/cc] M_cl [M_sol] V_cl [AU^3] rel.  peak_check   isodens_check   clump_check '
        do j=npeaks_tot,1,-1
           jj=ind_sort(j)
           if (relevance_tot(jj) > 0)then
              if(verbose)write(ilun,'(I6,X,I10,16(1X,1PE14.7))')jj&          
                   ,n_cells_tot(jj)&
                   ,peak_pos_tot(jj,1),peak_pos_tot(jj,2),peak_pos_tot(jj,3)&
                   ,(5.*clump_size_tot(jj,1)/clump_vol_tot(jj))**0.5*scale_l &
                   ,(5.*clump_size_tot(jj,2)/clump_vol_tot(jj))**0.5*scale_l &
                   ,(5.*clump_size_tot(jj,3)/clump_vol_tot(jj))**0.5*scale_l &
                   ,(clump_momentum_tot(jj,1)**2+clump_momentum_tot(jj,2)**2+ &
                   clump_momentum_tot(jj,3)**2)**0.5/clump_mass_tot(jj)*scale_l/scale_t&
                   ,min_dens_tot(jj)*scale_nH,max_dens_tot(jj)*scale_nH&
                   ,clump_mass_tot(jj)/clump_vol_tot(jj)*scale_nH&
                   ,clump_mass_tot(jj)*scale_d*dble(scale_l)**3/1.98892d33&
                   ,clump_vol_tot(jj)*(scale_l)**3&
                   ,relevance_tot(jj)&
                   ,peak_check(jj)&
!                   ,ball4_check(jj)&
                   ,isodens_check(jj)&
                   ,clump_check(jj)
              
              rel_mass=rel_mass+clump_mass_tot(jj)*scale_d*dble(scale_l)**3/1.98892d33
              n_rel=n_rel+1
           end if
        end do
        
     else
        if(clinfo .and. (to_file .eqv. .false.))write(ilun,'(135A)')' Cl_N #leaf-cells peak_x [uu] peak_y [uu] peak_z [uu] size_x [AU]'//&
             ' size_y [AU] size_z [AU]  |v|_CM [u.u.]  rho- [H/cc]  rho+ [H/cc] rho_av[H/cc] M_cl[M_sol] V_cl [AU^3]    rel.   '//&
             ' peak_check  clump_check '
        do j=npeaks_tot,1,-1
           jj=ind_sort(j)
           
           if (relevance_tot(jj) > 0)then
              if(clinfo .and. (to_file .eqv. .false.))then
                 write(ilun,'(I6,X,I10,3(X,F11.5),3(X,F11.5),X,F13.5,3(X,E12.3E2),5(X,E11.2E2))')&
                      jj&
                      ,n_cells_tot(jj)&
                      ,peak_pos_tot(jj,1)&
                      ,peak_pos_tot(jj,2)&
                      ,peak_pos_tot(jj,3)&
                      ,(5.*clump_size_tot(jj,1)/clump_vol_tot(jj))**0.5*(scale_l/1.496d13)&
                      ,(5.*clump_size_tot(jj,2)/clump_vol_tot(jj))**0.5*(scale_l/1.496d13)&
                      ,(5.*clump_size_tot(jj,3)/clump_vol_tot(jj))**0.5*scale_l/1.496d13&
                      ,(clump_momentum_tot(jj,1)**2+clump_momentum_tot(jj,2)**2+clump_momentum_tot(jj,3)**2)**0.5/clump_mass_tot(jj)&
                      ,min_dens_tot(jj)*scale_nH&
                      ,max_dens_tot(jj)*scale_nH&
                      ,clump_mass_tot(jj)/clump_vol_tot(jj)*scale_nH&
                      ,clump_mass_tot(jj)*scale_d*dble(scale_l)**3/1.98892d33&
                      ,clump_vol_tot(jj)*(scale_l/1.496d13)**3&
                      ,relevance_tot(jj)&
                      ,isodens_check(jj)&
                      ,clump_check(jj)
                 
                 rel_mass=rel_mass+clump_mass_tot(jj)*scale_d*dble(scale_l)**3/1.98892d33
                 n_rel=n_rel+1
              end if
           end if
        end do
     end if
     if(to_file)then
        write(21,*)n_rel
        do j=npeaks_tot,1,-1
           jj=ind_sort(j)
           if (relevance_tot(jj)>0)write(21,*)clump_mass_tot(jj)*scale_d*dble(scale_l)**3/1.98892d33
        end do
     else
        if(clinfo .and. (to_file .eqv. .false.))write(ilun,'(A,1PE12.5)')'total mass above threshold =',tot_mass*scale_d*dble(scale_l)**3/1.98892d33
        if(clinfo .and. (to_file .eqv. .false.))write(ilun,'(A,I6,A,1PE12.5)')'total mass in',n_rel,' listed clumps =',rel_mass
     endif
     if (to_file)then
        close(20)
        close(21)
     end if
  end if

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

end subroutine write_clump_properties
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine saddlepoint_search(ntest)
  use amr_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ntest

  !---------------------------------------------------------------------------
  ! subroutine which creates a npeaks_tot**2 sized array of saddlepoint densities
  ! by looping over all testparticles and passing them to neighborcheck
  ! with case 4, which means that saddlecheck will be called for each neighboring
  ! leaf cell. There it is checked, whether the two cells (original cell and 
  ! neighboring cell) are connected by a new densest saddle.
  !---------------------------------------------------------------------------

  integer::ipart,ip,ilevel,next_level
  integer::i,j,info,dummyint
  integer,dimension(1:nvector)::ind_cell
  real(dp),allocatable,dimension(:)::temp,temp_tot

  ! saddle array for 1 cpu
  allocate(saddle_dens(1:npeaks_tot,1:npeaks_tot))
  saddle_dens=0.

  ! loop 'testparts', pass the information of nvector parts to neighborsearch 
  ip=0
  do ipart=1,ntest
     ip=ip+1
     ilevel=levp(testp_sort(ipart)) ! level
     next_level=0
     if(ipart<ntest)next_level=levp(testp_sort(ipart+1)) !level of next particle
     ind_cell(ip)=icellp(testp_sort(ipart))
     if (flag2(ind_cell(ip))==0)print*,'alert neighborsearch',ind_cell(ip),myid,ilevel
     if(ip==nvector .or. next_level /= ilevel)then
        call neighborsearch(ind_cell,ip,dummyint,ilevel,4)
        ip=0
     endif
  end do
  if (ip>0)call neighborsearch(ind_cell,ip,dummyint,ilevel,4)

  ! share the results among MPI domains (communicate line by line in case of big arrays) 
#ifndef WITHOUTMPI
  allocate(temp(1:npeaks_tot),temp_tot(1:npeaks_tot))
    do i=1,npeaks_tot
     temp(1:npeaks_tot)=saddle_dens(1:npeaks_tot,i)
     temp_tot=0.d0
     call MPI_ALLREDUCE(temp,temp_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info) 
     saddle_dens_tot(1:npeaks_tot,i)=temp_tot(1:npeaks_tot)
  end do
  deallocate(temp,temp_tot)
#endif
#ifdef WITHOUTMPI
  saddle_dens_tot=saddle_dens
#endif

  ! check symmetry 
  if (debug)then
     do i=1,npeaks_tot
        do j=1,i
           if(saddle_dens_tot(i,j)/=saddle_dens_tot(j,i).and.myid==1)then 
              write(*,*),'Alert! asymmetric saddle point array!',i,j,saddle_dens_tot(i,j),saddle_dens_tot(j,i)
           endif
        end do
     end do
  end if

  ! compute saddle_max value and relevance
  saddle_max_tot=maxval(saddle_dens_tot,dim=1)
  do i=1,npeaks_tot
     if (saddle_max_tot(i)>0.)then
        relevance_tot(i)=max_dens_tot(i)/saddle_max_tot(i)
     else
        relevance_tot(i)=max_dens_tot(i)/min_dens_tot(i)
     end if
  end do
  
  !from here only saddle_dens_tot is used
  deallocate(saddle_dens)

end subroutine saddlepoint_search
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine merge_clumps(ntest)
  use amr_commons
  use clfind_commons
  implicit none
  integer::ntest

  !---------------------------------------------------------------------------
  ! This routine merges the irrelevant clumps 
  ! -clumps are sorted by ascending max density
  ! -irrelevent clumps are merged to most relevant neighbor
  !---------------------------------------------------------------------------

  integer::j,i,ii,merge_count,final_peak,merge_to,ipart
  integer::peak,next_peak
  real(dp)::max_val
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,d0
  integer,dimension(1:npeaks_tot)::old_peak,ind_sort
  real(dp),dimension(1:npeaks_tot)::peakd

  if (verbose)write(*,*)'Now merging clumps'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  d0 = density_threshold/scale_nH;
  if(cosmo)d0=d0/aexp**3


  ! Sort clumps by peak density in ascending order
  do i=1,npeaks_tot
     peakd(i)=max_dens_tot(i)
     ind_sort(i)=i
  end do
  call quick_sort_dp(peakd,ind_sort,npeaks_tot) 


  if (smbh .eqv. .false.) then
     do i=1,npeaks_tot
        ii=ind_sort(i)
        new_peak(ii)=ii
        
        ! If the relevance is below the threshold -> merge
        if (relevance_tot(ii)<relevance_threshold.and.relevance_tot(ii)>.5) then
           
           ! Go through the ii-th line in the saddle point array to find the neighbor to merge to
           merge_to=0; max_val=0.
           do j=1,npeaks_tot
              if (saddle_dens_tot(ii,j)>max_val)then
                 merge_to=j
                 max_val=saddle_dens_tot(ii,j)
              end if
           end do
           
           ! Store new peak index
           new_peak(ii)=merge_to
           if(clinfo .and. myid==1)then
!              if(merge_to>0)then
              write(*,*)'clump ',ii,'merged to ',merge_to
 !             endif
           endif
           
           ! Update clump properties
           if (merge_to>0)then
              do j=1,ndim
                 clump_momentum_tot(merge_to,j)=&
                      clump_momentum_tot(merge_to,j)+clump_momentum_tot(ii,j)
                 clump_momentum_tot(ii,j)=0.
                 center_of_mass_tot(merge_to,j)=&
                      (clump_mass_tot(merge_to)*center_of_mass_tot(merge_to,j) &
                      +clump_mass_tot(ii)*center_of_mass_tot(ii,j)) &
                      /(clump_mass_tot(merge_to)+clump_mass_tot(ii))
                 center_of_mass_tot(ii,j)=0.
                 clump_force_tot(merge_to,j)=&
                      (clump_mass_tot(merge_to)*clump_force_tot(merge_to,j) &
                      +clump_mass_tot(ii)*clump_force_tot(ii,j)) &
                      /(clump_mass_tot(merge_to)+clump_mass_tot(ii))
                 clump_force_tot(ii,j)=0.
              end do
              n_cells_tot(merge_to)=n_cells_tot(merge_to)+n_cells_tot(ii)
              clump_vol_tot(merge_to)=clump_vol_tot(ii)+clump_vol_tot(merge_to)
              max_dens_tot(merge_to)=max(max_dens_tot(merge_to),max_dens_tot(ii))
              min_dens_tot(merge_to)=min(min_dens_tot(merge_to),min_dens_tot(ii))
              clump_mass_tot(merge_to)=clump_mass_tot(merge_to)+clump_mass_tot(ii)
              phi_min_tot(merge_to)=min(phi_min_tot(merge_to),phi_min_tot(ii))
           end if
           n_cells_tot(ii)=0
           clump_vol_tot(ii)=0.
           max_dens_tot(ii)=0.
           min_dens_tot(ii)=0.
           clump_mass_tot(ii)=0.
           
           ! Update saddle point array
           do j=1,npeaks_tot
              if (merge_to>0)then
                 if(saddle_dens_tot(ii,j)>saddle_dens_tot(merge_to,j))then
                    saddle_dens_tot(merge_to,j)=saddle_dens_tot(ii,j)
                    saddle_dens_tot(j,merge_to)=saddle_dens_tot(ii,j)
                 end if
                 saddle_dens_tot(merge_to,merge_to)=0.
              end if
              saddle_dens_tot(ii,j)=0.        
              saddle_dens_tot(j,ii)=0.
           end do
           
           ! Update saddle_max value
           if (merge_to>0)then
              saddle_max_tot(merge_to)=0
              do j=1,npeaks_tot
                 if (saddle_dens_tot(merge_to,j)>saddle_max_tot(merge_to))then
                    saddle_max_tot(merge_to)=saddle_dens_tot(merge_to,j)
                 end if
              end do
           end if
           
           ! Update relevance of clumps
           if (merge_to>0)then
              if (saddle_max_tot(merge_to)>1.d-40)then
                 relevance_tot(merge_to)=max_dens_tot(merge_to)/saddle_max_tot(merge_to)
              else 
                 relevance_tot(merge_to)=max_dens_tot(merge_to)/min_dens_tot(merge_to)
              end if
           end if
           relevance_tot(ii)=0.
        end if
     end do
  else
     ! If SMBH merge all peaks
     do i=1,npeaks_tot
        ii=ind_sort(i)
        new_peak(ii)=ii
        
        ! Go through the ii-th line in the saddle point array to find the neighbor to merge to
        merge_to=0; max_val=0.
        do j=1,npeaks_tot
           if (saddle_dens_tot(ii,j)>max_val)then
              merge_to=j
              max_val=saddle_dens_tot(ii,j)
           end if
        end do
        
        ! Store new peak index
        if(merge_to>0)new_peak(ii)=merge_to
        if(verbose .and. myid==1)then
           if(merge_to>0)then
              write(*,*)'clump ',ii,'merged to ',merge_to
           endif
        endif
        
        ! Update clump properties
        if (merge_to>0)then
           do j=1,ndim
              clump_momentum_tot(merge_to,j)=&
                   clump_momentum_tot(merge_to,j)+clump_momentum_tot(ii,j)
              clump_momentum_tot(ii,j)=0.
              center_of_mass_tot(merge_to,j)=&
                   (clump_mass_tot(merge_to)*center_of_mass_tot(merge_to,j) &
                   +clump_mass_tot(ii)*center_of_mass_tot(ii,j)) &
                   /(clump_mass_tot(merge_to)+clump_mass_tot(ii))
              center_of_mass_tot(ii,j)=0.
           end do
           n_cells_tot(merge_to)=n_cells_tot(merge_to)+n_cells_tot(ii)
           clump_vol_tot(merge_to)=clump_vol_tot(ii)+clump_vol_tot(merge_to)
           max_dens_tot(merge_to)=max(max_dens_tot(merge_to),max_dens_tot(ii))
           min_dens_tot(merge_to)=min(min_dens_tot(merge_to),min_dens_tot(ii))
           clump_mass_tot(merge_to)=clump_mass_tot(merge_to)+clump_mass_tot(ii)
           phi_min_tot(merge_to)=min(phi_min_tot(merge_to),phi_min_tot(ii))
           n_cells_tot(ii)=0
           clump_vol_tot(ii)=0.
           max_dens_tot(ii)=0.
           min_dens_tot(ii)=0.
           clump_mass_tot(ii)=0.
        end if
        
        ! Update saddle point array
        do j=1,npeaks_tot
           if (merge_to>0)then
              if(saddle_dens_tot(ii,j)>saddle_dens_tot(merge_to,j))then
                 saddle_dens_tot(merge_to,j)=saddle_dens_tot(ii,j)
                 saddle_dens_tot(j,merge_to)=saddle_dens_tot(ii,j)
              end if
              saddle_dens_tot(merge_to,merge_to)=0.
              saddle_dens_tot(ii,j)=0.
              saddle_dens_tot(j,ii)=0.
           end if
        end do
        
        ! Update saddle_max value
        if (merge_to>0)then
           saddle_max_tot(merge_to)=0
           do j=1,npeaks_tot
              if (saddle_dens_tot(merge_to,j)>saddle_max_tot(merge_to))then
                 saddle_max_tot(merge_to)=saddle_dens_tot(merge_to,j)
              end if
           end do
        end if
        
        ! Update relevance of clumps
        if (merge_to>0)then
           if (saddle_max_tot(merge_to)>1.d-40)then
              relevance_tot(merge_to)=max_dens_tot(merge_to)/saddle_max_tot(merge_to)
           else
              relevance_tot(merge_to)=max_dens_tot(merge_to)/min_dens_tot(merge_to)
           end if
           relevance_tot(ii)=0.
        end if
     end do
  end if

  if (verbose)write(*,*)'Done merging clumps 0'

! Change new_peak so that it points to the end point of the merging 
! history and not only to the clump it has been merged to in first place 
  do i=1,npeaks_tot
     peak=i
     merge_count=0
     next_peak=new_peak(peak)
     do while(peak.NE.next_peak.AND.next_peak>0)
        !construct old_peak to walk the tree from trunk to leafs
        old_peak(next_peak)=peak
        !go to next peak
        merge_count=merge_count+1
        peak=next_peak
        next_peak=new_peak(peak)
     end do
     final_peak=next_peak
     !now we walk the other way, so next_peak will always be the previous one
     next_peak=peak
     do j=merge_count,1,-1
        next_peak=old_peak(next_peak)
        new_peak(next_peak)=final_peak
     end do
  end do

  if (verbose)write(*,*)'Done merging clumps I'

! Remove peaks using HOP-based criterion
  if (smbh)then
     do i=1,npeaks_tot
        ii=ind_sort(i)
        if( max_dens_tot(new_peak(ii)) < 3.0*d0 )then
           new_peak(ii)=0
           n_cells_tot(ii)=0
           clump_vol_tot(ii)=0.
           max_dens_tot(ii)=0.
           min_dens_tot(ii)=0.
           clump_mass_tot(ii)=0.
           ! Update saddle point array
           do j=1,npeaks_tot
              saddle_dens_tot(ii,j)=0.
              saddle_dens_tot(j,ii)=0.
           end do
           relevance_tot(ii)=0.        
        end if
     end do
     if (verbose)write(*,*)'Done merging clumps II'  
  endif

  !update flag 2
  do ipart=1,ntest
     if (flag2(icellp(ipart))>0)flag2(icellp(ipart))=new_peak(flag2(icellp(ipart)))
  end do


end subroutine merge_clumps
!################################################################                 
!################################################################ 
!################################################################                 
!################################################################     
subroutine allocate_peak_patch_arrays
  use amr_commons, ONLY:ndim,dp
  use clfind_commons
  implicit none
  real(dp)::zero=0.

  ! Allocate peak-patch_properties
  allocate(n_cells_tot(1:npeaks_tot))
  allocate(clump_size_tot(1:npeaks_tot,1:ndim))
  allocate(peak_pos_tot(1:npeaks_tot,1:ndim))
  allocate(center_of_mass_tot(1:npeaks_tot,1:ndim))
  allocate(clump_force_tot(1:npeaks_tot,1:ndim))
  allocate(second_moments(1:npeaks_tot,1:ndim,1:ndim)) 
  allocate(second_moments_tot(1:npeaks_tot,1:ndim,1:ndim))
  allocate(min_dens_tot(1:npeaks_tot))
  allocate(av_dens_tot(1:npeaks_tot))
  allocate(max_dens_tot(1:npeaks_tot))
  allocate(clump_mass_tot(1:npeaks_tot))
  allocate(clump_mass_tot4(1:npeaks_tot))
  allocate(clump_vol_tot(1:npeaks_tot))
  allocate(saddle_max_tot(1:npeaks_tot))
  allocate(relevance_tot(1:npeaks_tot))
  allocate(saddle_dens_tot(1:npeaks_tot,1:npeaks_tot))
  allocate(clump_momentum_tot(1:npeaks_tot,1:ndim))
  allocate(e_kin_int_tot(npeaks_tot))
  allocate(e_bind_tot(npeaks_tot))
  allocate(e_thermal_tot(npeaks_tot))
  allocate(phi_min_tot(npeaks_tot))
  allocate(minmatch_tot(npeaks_tot))
  allocate(new_peak(npeaks_tot))
  allocate(phi_ref(npeaks_tot))
  allocate(phi_ref_tot(npeaks_tot))
  allocate(Psurf(npeaks_tot))
  allocate(Psurf_tot(npeaks_tot))
  allocate(v_therm_tot(npeaks_tot))
  allocate(m4_tot(npeaks_tot))
  allocate(bulk_momentum_tot(1:npeaks_tot,1:ndim))
  allocate(e_kin_iso_tot(npeaks_tot))
  allocate(e_bind_iso_tot(npeaks_tot))
  allocate(e_therm_iso_tot(npeaks_tot))
  allocate(peak_check(npeaks_tot))
  allocate(isodens_check(npeaks_tot))
  allocate(clump_check(npeaks_tot))
  allocate(grav_term_tot(npeaks_tot))
  allocate(contracting(npeaks_tot))
  allocate(Icl_tot(npeaks_tot))
  allocate(Icl_d_tot(npeaks_tot))
  allocate(Icl_dd_tot(npeaks_tot))
  allocate(Icl_d_3by3_tot(npeaks_tot,1:3,1:3))  
  allocate(Icl_3by3_tot(npeaks_tot,1:3,1:3))  

  !initialize all peak based arrays
  n_cells_tot=0
  saddle_max_tot=0.
  relevance_tot=1.
  clump_size_tot=0.
  min_dens_tot=huge(zero)
  max_dens_tot=0.
  av_dens_tot=0.
  clump_mass_tot=0.
  clump_mass_tot4=0.
  clump_vol_tot=0.
  peak_pos_tot=0.
  center_of_mass_tot=0.
  clump_force_tot=0.
  second_moments=0.; second_moments_tot=0.
  saddle_dens_tot=0.
  clump_momentum_tot=0.
  e_kin_int_tot=0.
  e_bind_tot=0.
  e_thermal_tot=0.
  phi_min_tot=0.
  minmatch_tot=1
  new_peak=0
  phi_ref=huge(zero)
  Psurf=0.;Psurf_tot=0.
  grav_term_tot=0.d0
  isodens_check=-1.;  clump_check=-1.; peak_check=-1.; 
  contracting=.false.
  Icl_tot=0.; Icl_d_tot=0.; Icl_dd_tot=0.; Icl_d_3by3_tot=0.; Icl_3by3_tot=0.

end subroutine allocate_peak_patch_arrays
!################################################################                 
!################################################################ 
!################################################################                 
!################################################################     
subroutine deallocate_all
  use clfind_commons
  use amr_commons, only:smbh
  implicit none

  deallocate(n_cells_tot)
  deallocate(clump_size_tot)
  deallocate(peak_pos_tot)
  deallocate(center_of_mass_tot,clump_force_tot)
  deallocate(second_moments)
  deallocate(second_moments_tot)
  deallocate(min_dens_tot)
  deallocate(av_dens_tot)
  deallocate(max_dens_tot)
  deallocate(clump_mass_tot)
  deallocate(clump_vol_tot)
  deallocate(saddle_max_tot)
  deallocate(relevance_tot)
  deallocate(saddle_dens_tot)
  deallocate(clump_momentum_tot)
  deallocate(e_kin_int_tot)
  deallocate(e_bind_tot,grav_term_tot)
  deallocate(e_thermal_tot)
  deallocate(phi_min_tot)
  deallocate(minmatch_tot)
  deallocate(new_peak)
  deallocate(phi_ref,phi_ref_tot)
  deallocate(Psurf,Psurf_tot)
  deallocate(v_therm_tot)
  deallocate(m4_tot,bulk_momentum_tot)
  deallocate(e_kin_iso_tot,e_bind_iso_tot,e_therm_iso_tot)
  deallocate(peak_check,isodens_check,clump_check)
  deallocate(contracting)
  deallocate(Icl_dd_tot,Icl_d_tot,Icl_tot,Icl_d_3by3_tot,Icl_3by3_tot)
  if (.not. smbh)deallocate(clump_mass_tot4)

end subroutine deallocate_all
!################################################################                 
!################################################################ 
!################################################################                 
!################################################################     
subroutine get_phi_ref(ntest)
  use amr_commons
  use hydro_commons
  use pm_commons
  use clfind_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ntest

  !---------------------------------------------------------------
  ! This subroutine finds the minimum potential on the saddle 
  ! surface of the peak patch by looping over all "test-particles"
  ! This loop is also used to compute the surface pressure term.
  !---------------------------------------------------------------

  integer::info   
  integer::ipart,ip,ilevel,next_level
  integer,dimension(1:nvector)::ind_cell

  Psurf=0.
  ! loop 'testparts', pass the information of nvector parts to neighborsearch 
  ip=0
  do ipart=1,ntest
     ip=ip+1
     ilevel=levp(testp_sort(ipart)) ! level
     if (verbose.and.ilevel/=nlevelmax)print*,'not all particles in max level',ilevel
     next_level=0
     if(ipart<ntest)next_level=levp(testp_sort(ipart+1)) !level of next particle
     ind_cell(ip)=icellp(testp_sort(ipart))
     if(ip==nvector .or. next_level /= ilevel)then
        call neighborsearch(ind_cell,ip,0,ilevel,5)
        call surface_int(ind_cell,ip,ilevel)
        ip=0
     endif
  end do
  if (ip>0)then 
     call neighborsearch(ind_cell,ip,0,ilevel,5)
     call surface_int(ind_cell,ip,ilevel)
  endif
   
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(phi_ref,phi_ref_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(Psurf,Psurf_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  phi_ref_tot=phi_ref
  Psurf_tot=Psurf
#endif

end subroutine get_phi_ref
!################################################################                 
!################################################################ 
!################################################################                 
!################################################################     
subroutine trim_clumps(ntest)
  use amr_commons
  use clfind_commons
  use pm_commons, only:ir_cloud
  implicit none
  integer::ntest

  !---------------------------------------------------------------------------
  ! this routine trims the clumps down to the intersection of the clump with 
  ! the accretion zone of the sink. Cells that are too far away from the peak
  ! are removed from the clump by setting flag2 to 0.
  !---------------------------------------------------------------------------

  integer::ipart,nx_loc,ind
  real(dp)::dx,scale,dx_loc,r2
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  integer ::ix,iy,iz,grid,peak_nr

  real(dp),dimension(1:3)::skip_loc,xcell
  real(dp),dimension(1:twotondim,1:3)::xc

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  
  ! Mesh spacing in max level
  dx=0.5D0**nlevelmax
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)
     xc(ind,2)=(dble(iy)-0.5D0)
     xc(ind,3)=(dble(iz)-0.5D0)
  end do

  !update flag 2
  do ipart=1,ntest
     peak_nr=flag2(icellp(ipart))
     if (peak_nr /=0 ) then
        ind=(icellp(ipart)-ncoarse-1)/ngridmax+1 ! cell position
        grid=icellp(ipart)-ncoarse-(ind-1)*ngridmax ! grid index
        dx=0.5D0**levp(ipart)
        xcell(1:ndim)=(xg(grid,1:ndim)+xc(ind,1:ndim)*dx-skip_loc(1:ndim))*scale
        r2=(peak_pos_tot(peak_nr,1)-xcell(1))**2&
             +(peak_pos_tot(peak_nr,2)-xcell(2))**2&
             +(peak_pos_tot(peak_nr,3)-xcell(3))**2.
        if (r2 > (ir_cloud*dx_loc)**2.)then        
           !remove cell from clump
           flag2(icellp(ipart))=0
        end if
     end if
  end do

end subroutine trim_clumps
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine jacobi(A,x,err2)
  use amr_commons, only:myid,dp
  implicit none
  real(dp)::err2
  real(dp),dimension(3,3)::A,x

  !---------------------------------------------------------------------------
  ! Compute eigenvalues and eigenvectors using the jacobi-Method 
  ! as for example described in Numerical Recipes. 
  ! Returns eigenvalues as diagonal elements of A
  !---------------------------------------------------------------------------

  integer::n
  integer::i,j,k
  real(dp)::b2, bar
  real(dp)::beta, coeff, c, s, cs, sc
  
  n=3
  ! x is identity matrix initially
  x = 0.0
  do i=1,n
     x(i,i) = 1.0
  end do

  ! sum all squared off-diagonal elements 
  b2 = 0.0
  do i=1,n
     do j=1,n
        if (i.ne.j) b2 = b2 + A(i,j)**2
     end do
  end do

  if (b2 <= err2) then
     if (myid==1)write(*,*), 'returning. maybe err2 too small? ',err2
     return
  endif

  ! average for off-diagonal elements /2
  bar = 0.5*b2/9.

  do while (b2 > err2)
     do i=1,n-1
        do j=i+1,n
           if (A(j,i)**2 <= bar) cycle  ! do not touch small elements
           b2 = b2 - 2.0*A(j,i)**2
           bar = 0.5*b2/9.
           ! calculate coefficient c and s for Givens matrix
           beta = (A(j,j)-A(i,i))/(2.0*A(j,i))
           coeff = 0.5*beta*(1.0+beta**2)**(-0.5)
           s = (max(0.5+coeff,0.0))**0.5
           c = (max(0.5-coeff,0.0))**0.5
           ! update rows i and j
           do k=1,n
              cs =  c*A(i,k)+s*A(j,k)
              sc = -s*A(i,k)+c*A(j,k)
              A(i,k) = cs
              A(j,k) = sc
           end do
           ! find new matrix A_{k+1} 
           do k=1,n
              cs =  c*A(k,i)+s*A(k,j)
              sc = -s*A(k,i)+c*A(k,j)
              A(k,i) = cs
              A(k,j) = sc
              cs =  c*x(k,i)+s*x(k,j)
              sc = -s*x(k,i)+c*x(k,j)
              x(k,i) = cs
              x(k,j) = sc
           end do
        end do
     end do
  end do
end subroutine jacobi
!################################################################
!################################################################
!################################################################
!################################################################
subroutine write_clump_map(ntest)
  use amr_commons
  use clfind_commons
  implicit none
  integer::ntest

  !---------------------------------------------------------------------------
  ! This routine writes a csv-file of cell center coordinates and clump number
  ! for each cell which is in a clump. Makes only sense to be called when the 
  ! clump finder is called at output-writing and not for sink-formation.
  !---------------------------------------------------------------------------

  integer::ind,grid,ix,iy,iz,ipart,nx_loc,peak_nr
  real(dp)::scale,dx
  real(dp),dimension(1:3)::xcell,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  character(LEN=5)::myidstring,nchar 

  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)

  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)
     xc(ind,2)=(dble(iy)-0.5D0)
     xc(ind,3)=(dble(iz)-0.5D0)
  end do

  !prepare file output for peak map
  call title(ifout-1,nchar)
  call title(myid,myidstring)
  open(unit=20,file=TRIM('output_'//TRIM(nchar)//'/clump_map.csv'//myidstring),form='formatted')

  !loop parts
  do ipart=1,ntest     
     peak_nr=flag2(icellp(ipart)) 
     if (peak_nr /=0 ) then
        ! Cell coordinates
        ind=(icellp(ipart)-ncoarse-1)/ngridmax+1 ! cell position
        grid=icellp(ipart)-ncoarse-(ind-1)*ngridmax ! grid index
        dx=0.5D0**levp(ipart)
        xcell(1:ndim)=(xg(grid,1:ndim)+xc(ind,1:ndim)*dx-skip_loc(1:ndim))*scale
        !peak_map
        write(20,'(F11.8,A,F11.8,A,F11.8,A,I8)')xcell(1),',',xcell(2),',',xcell(3),',',peak_nr
     end if
  end do
  close(20)
end subroutine write_clump_map
