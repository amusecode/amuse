subroutine courant_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::info,nleaf,ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

  real(dp)::dt_lev,dx,vol,scale
  real(kind=8)::mass_loc,ekin_loc,eint_loc,dt_loc
  real(kind=8)::mass_all,ekin_all,eint_all,dt_all
  real(kind=8),dimension(3)::comm_buffin,comm_buffout
  real(dp),dimension(1:nvector,1:nvar),save::uu
  real(dp),dimension(1:nvector,1:ndim),save::gg

  if(ilevel==nlevelmax)call bondi_fine(ilevel)

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0
  dt_all=dtnew(ilevel); dt_loc=dt_all

  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
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
#ifndef WITHOUTMPI
  comm_buffin(1)=mass_loc
  comm_buffin(2)=ekin_loc
  comm_buffin(3)=eint_loc
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,3,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dt_loc     ,dt_all      ,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
  mass_all=comm_buffout(1)
  ekin_all=comm_buffout(2)
  eint_all=comm_buffout(3)
#endif
#ifdef WITHOUTMPI
  mass_all=mass_loc
  ekin_all=ekin_loc
  eint_all=eint_loc
  dt_all=dt_loc
#endif
  mass_tot=mass_tot+mass_all
  ekin_tot=ekin_tot+ekin_all
  eint_tot=eint_tot+eint_all
  dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)

111 format('   Entering courant_fine for level ',I2)

end subroutine courant_fine

subroutine bondi_fine(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  
  integer::igrid,ncache,i,ind,iskip,ngrid,info
  integer::ix,iy,iz,idim,ivar,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::xjet,yjet,zjet,rjet,rhoSjet,vjet,tjet,rhojet,rjetcgs
  real(dp)::dtoverdx,dx,scale,mdot,dx_loc,vol_loc,f0_jet,hjet
  real(kind=8)::Mcylinder,Vcylinder,Ecylinder,Mcylinder_all,Vcylinder_all,Ecylinder_all
  real(dp)::rho_cylinder,cs2_cylinder
  real(dp)::Mvirphu,rvir,cnfw,eta0,gammaKS,d_init,p_init,p_fudge_factor
  real(dp)::delta_rho
  real(dp),dimension(1:twotondim,1:ndim)::xc
  real(dp),dimension(1:nvector),save::rr,zz,ekk,ett
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:3)::skip_loc

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh size
  dx=0.5d0**ilevel

  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**3
  
  ! Jet parameters in cgs
  rjet=r_jet*dx_loc               ! Jet radius
  hjet=h_jet*dx_loc               ! Jet height
  tjet=t_jet*(1d6*365.*24.*3600.) ! Jet life time
  xjet=boxlen/2.
  yjet=boxlen/2.
  zjet=boxlen/2.

  ! Jet ejection rate
  f0_jet=Mdot_BH/(6.28*rjet**2*hjet)

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  Mcylinder=0.0
  Ecylinder=0.0
  Vcylinder=0.0

  ! Loop over active grids by vector sweeps
  if(t.gt.0.0.and.t.lt.tjet)then
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
        
        ! Gather cell centre positions
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
           end do
        end do
        ! Rescale position from code units to user units
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
           end do
        end do

        do i=1,ngrid
           rr(i)=sqrt((xx(i,1)-xjet)**2+(xx(i,2)-yjet)**2)
        end do
        do i=1,ngrid
           zz(i)=(xx(i,3)-zjet)
        end do

        ! Compute volume internal energy
        do i=1,ngrid
           ett(i)=uold(ind_cell(i),ndim+2)
        end do
        do i=1,ngrid
           ekk(i)=0.0d0
        end do
        do idim=1,ndim
           do i=1,ngrid
              ekk(i)=ekk(i)+0.5*uold(ind_cell(i),idim+1)**2/uold(ind_cell(i),1)
           end do
        end do
        do i=1,ngrid
           ett(i)=ett(i)-ekk(i)
        end do

        ! Cylinder
        do i=1,ngrid
           if(rr(i)<rjet.and.abs(zz(i))<hjet)then
              Vcylinder=Vcylinder+vol_loc
              Mcylinder=Mcylinder+uold(ind_cell(i),1)*vol_loc
              Ecylinder=Ecylinder+ett(i)*vol_loc
              delta_rho = f0_jet*exp(-5.0*(rr(i)/rjet)**2)*2.*zz(i)/hjet*dtold(ilevel)
              uold(ind_cell(i),1)=uold(ind_cell(i),1)+abs(delta_rho)*e_jet
              uold(ind_cell(i),4)=uold(ind_cell(i),4)+delta_rho*sqrt(0.2)*3d10
              uold(ind_cell(i),5)=uold(ind_cell(i),5)+abs(delta_rho)*0.1*(3d10)**2
           end if
        end do

     end do

  end do
  endif

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(Vcylinder,Vcylinder_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(Mcylinder,Mcylinder_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(Ecylinder,Ecylinder_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  Vcylinder=Vcylinder_all
  Mcylinder=Mcylinder_all
  Ecylinder=Ecylinder_all
#endif

  ! Update Black Hole accretion rate
  if(t.gt.0.0.and.t.lt.tjet.and.myid==1)then
     rho_cylinder=Mcylinder/Vcylinder
     cs2_cylinder=gamma*(gamma-1.0)*Ecylinder/Mcylinder
     Mdot_BH=2d72*rho_cylinder/cs2_cylinder**1.5     
     if(myid==1)write(*,*)'Mdot_BH (Msol/yr)=',Mdot_BH/2d33*(3600.*24.*365.)
  endif

  ! Update boundaries
  do ivar=1,nvar
     call make_virtual_fine_dp(uold(1,ivar),ilevel)
  end do

111 format('   Entering bondi_fine for level ',I2)

end subroutine bondi_fine
