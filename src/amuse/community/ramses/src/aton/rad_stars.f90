! Adapted from pm/feedback.f90.

subroutine compute_Srad()
  use pm_commons
  use amr_commons
  use hydro_commons
  use radiation_commons
  use observe_commons
  implicit none
  integer::ilevel

  observe_total_star_source=0d0
  observe_num_stars=0

  if (.not.pic) then
     if (myid.eq.1) then
        write(*,*)"Not computing Srad for star sources because pic=false."
     end if
     return
  end if

  if (myid.eq.1) then
     write(*,*)"Computing Srad."
  end if

  Srad = 0.0
  call compute_Srad_level(levelmin)

end subroutine compute_Srad

subroutine compute_Srad_level(ilevel)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::t0,scale,dx_min,vsn,rdebris,ethermal
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc
  real(dp),dimension(1:3)::skip_loc
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Gather star particles only.

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count star particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).gt.0.and.tp(ipart).ne.0)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather star particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only star particles
              if(idp(ipart).gt.0.and.tp(ipart).ne.0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                 call process_particle(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip>0)call process_particle(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

111 format('   Entering compute_Srad for level ',I2)

end subroutine compute_Srad_level


subroutine process_particle(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use hydro_parameters, ONLY: imetal
  use observe_commons
  use radiation_commons
  implicit none
  real(kind=8),parameter ::mH      = 1.6600000d-24  ! from cooling_module.f90
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine feedback. Each debris particle
  ! dumps mass, momentum and energy in the nearest grid cell using array
  ! uold.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc
  real(dp)::xxx,mmm,t0,ESN,mejecta,zloss
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp)::age,ageprime,fraction
  real(dp)::star_source,mass,emission_rate
  real(dp)::time_step_begin,time_step_end
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  integer,dimension(1:nvector),save::list1
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mloss,mzloss,ethermal,ekinetic,dteff
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Massive star lifetime from Myr to code units
  t0=10.*1d6*(365.*24.*3600.)/scale_t

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<=0.0D0.or.x(j,idim)>=6.0D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in sn2'
     write(*,*)ilevel,ng,np
     stop
  end if

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=x(j,idim)
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igd(j,idim)=id(j,idim)/2
     end do
  end do
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
  do j=1,np
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do

  ! Compute parent cell adresses
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)

        ! Move up to levelmin.
        do i = 1,ilevel-levelmin
           indp(j) = father(indp(j))
        end do

     end if
  end do


  ! FIXME(tstranex): is the dt calculation ok since we are working on the coarse
  ! level only?

  ! Compute individual time steps                                               
  do j=1,np
     if(ok(j))then
        dteff(j)=dtnew(levelmin)
     endif
  end do

  ! Update the radiation source field based on the stars.

  ! UV photon emission rate
  ! emission_rate = 4800 [photons / Myr / baryon]
  emission_rate = 4800d0 / (1d6*365*24*60*60) / mH  ! [photons / s / g]

  do j=1,np
     if(ok(j))then
        age = t - tp(ind_part(j)) ! Age at the end of the step.
        ageprime = age - dteff(j) ! Age at the start of the step.
        fraction = 0.0
        if (ageprime.lt.0) then
           ! The star was formed mid step.
           fraction = age/dteff(j)
        else if (ageprime.lt.t0) then
           ! The star was formed before the current step and either dies mid
           ! step or later.
           fraction = min((t0-ageprime)/dteff(j), 1.0)
        end if

        mass = mp(ind_part(j)) * scale_d * scale_l**3  ! [g]
        star_source = fraction * &
             & mass*emission_rate*rad_escape_fraction  ! [photon/s]

        Srad(indp(j)) = Srad(indp(j)) + star_source

        observe_total_star_source = observe_total_star_source + star_source
        observe_num_stars = observe_num_stars + 1

     endif
  end do

end subroutine process_particle
