subroutine courant_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  include 'mpif.h'
  integer::ilevel
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::info ,nleaf,ngrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

  real(dp)::dt_lev,dx,vol
  real(kind=8)::mass_loc,ekin_loc,eint_loc,dt_loc
  real(kind=8)::mass_all,ekin_all,eint_all,dt_all
  real(dp),dimension(1:nvector,1:nvar),save::uu
  real(dp),dimension(1:nvector,1:ndim),save::gg

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0
  dt_all=dtnew(ilevel); dt_loc=dt_all

  ! Mesh spacing at that level
  dx=0.5D0**ilevel
  vol=dx**ndim

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim        
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=ind_grid(i)+iskip
        end do
        
        ! Gather leaf cells
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
           end if
        end do

        ! Gather hydro variables
        do ivar=1,nvar
           do i=1,nleaf
              uu(i,ivar)=uold(ind_leaf(i),ivar)
           end do
        end do
        
        ! Gather gravitational acceleration
        gg=0.0d0
        if(poisson)then
           do idim=1,ndim
              do i=1,nleaf
                 gg(i,idim)=f(ind_leaf(i),idim)
              end do
           end do
        end if
        
        ! Compute total mass
        do i=1,nleaf
           mass_loc=mass_loc+uu(i,1)*vol
        end do
        
        ! Compute total energy
        do i=1,nleaf
           ekin_loc=ekin_loc+uu(i,ndim+2)*vol
        end do
        
        ! Compute total internal energy
        do i=1,nleaf
           eint_loc=eint_loc+uu(i,ndim+2)*vol
        end do
        do ivar=1,ndim
           do i=1,nleaf
              eint_loc=eint_loc-0.5d0*uu(i,1+ivar)**2/uu(i,1)*vol
           end do
        end do
        
        ! Compute CFL time-step
        if(nleaf>0)then
           call cmpdt(uu,gg,dx,dt_lev,nleaf)
           dt_loc=min(dt_loc,dt_lev)
        end if
        
     end do
     ! End loop over cells
     
  end do
  ! End loop over grids

  ! Compute global quantities
  call MPI_ALLREDUCE(mass_loc,mass_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ekin_loc,ekin_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(eint_loc,eint_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(  dt_loc,  dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)

  mass_tot=mass_tot+mass_all
  ekin_tot=ekin_tot+ekin_all
  eint_tot=eint_tot+eint_all
  dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)

  ! Set accretor variables to vacuum (Bondi runs)
  if(ilevel==nlevelmax)call bondi_fine(ilevel)

111 format('   Entering courant_fine for level ',I2)

end subroutine courant_fine


subroutine bondi_fine(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_parameters
  implicit none
  integer::ilevel
  
  integer::igrid,ncache,i,ind,iskip,ngrid
  integer::ix,iy,iz,idim,ivar
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::xmass,ymass,zmass,emass,dx,scale
  real(dp),dimension(1:twotondim,1:ndim)::xc
  real(dp),dimension(1:nvector),save::rr
 
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh size
  dx=0.5d0**ilevel

  ! Rescaling factor
  scale=dble(icoarse_max-icoarse_min+1)/boxlen
  
  ! Scale black hole parameters
  emass=gravity_params(2)*scale
  xmass=gravity_params(3)*scale+dble(icoarse_min)
  ymass=gravity_params(4)*scale+dble(jcoarse_min)
  zmass=gravity_params(5)*scale+dble(kcoarse_min)

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx-xmass
!     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx-ymass
!     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx-zmass
  end do
  
  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        
        do i=1,ngrid
           rr(i)=0.0d0
        end do
        do i=1,ngrid
           rr(i)=rr(i)+(xg(ind_grid(i),1)+xc(ind,1))**2
        end do
#if NDIM>1
        do i=1,ngrid
           rr(i)=rr(i)+(xg(ind_grid(i),2)+xc(ind,2))**2
        end do
#endif
#if NDIM>2
        do i=1,ngrid
           rr(i)=rr(i)+(xg(ind_grid(i),3)+xc(ind,3))**2
        end do
#endif
        do i=1,ngrid
           rr(i)=sqrt(rr(i))
        end do

        do i=1,ngrid
           if(rr(i)<emass)then
              uold(ind_cell(i),1)=100.*smallr
           end if
        end do
        do i=1,ngrid
           if(rr(i)<emass)then
              uold(ind_cell(i),ndim+2)=100.*smallr*smallc**2/gamma/(gamma-1.0)
           end if
        end do
        do idim=1,ndim
           do i=1,ngrid
              if(rr(i)<emass)then
                 uold(ind_cell(i),idim+1)=0.0d0
              end if
           end do
        end do
     end do

  end do

  ! Update boundaries
  do ivar=1,nvar
     call make_virtual_fine_dp(uold(1,ivar),ilevel)
  end do

111 format('   Entering bondi_fine for level ',I2)

end subroutine bondi_fine
