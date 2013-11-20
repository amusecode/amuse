!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_boundary_hydro(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_parameters
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------
  integer::ibound,boundary_dir,idim,inbor
  integer::i,ncache,ivar,igrid,ngrid,ind,iperp1,iperp2,iplane,icell
  integer::iskip,iskip_ref,gdim,nx_loc,ix,iy,iz
  integer,dimension(1:8)::ind_ref,alt
  integer,dimension(1:2,1:4)::ind0
  integer,dimension(1:nvector),save::ind_grid,ind_grid_ref
  integer,dimension(1:nvector),save::ind_cell,ind_cell_ref

  real(dp)::switch,dx,dx_loc,scale,emag
  real(dp),dimension(1:3)::gs,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar+3),save::uu
  real(dp),dimension(1:nvector,1:twotondim),save::uu_ref

  if(.not. simple_boundary)return
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  ! Loop over boundaries
  do ibound=1,nboundary

     ! Compute direction of reference neighbors
     boundary_dir=boundary_type(ibound)-10*(boundary_type(ibound)/10)
     if(boundary_dir==1)inbor=2
     if(boundary_dir==2)inbor=1
     if(boundary_dir==3)inbor=4
     if(boundary_dir==4)inbor=3
     if(boundary_dir==5)inbor=6
     if(boundary_dir==6)inbor=5

     ! Compute index of reference cells
     ! Reflexive boundary
     if(boundary_type(ibound)== 1)ind_ref(1:8)=(/2,1,4,3,6,5,8,7/)
     if(boundary_type(ibound)== 2)ind_ref(1:8)=(/2,1,4,3,6,5,8,7/)
     if(boundary_type(ibound)== 3)ind_ref(1:8)=(/3,4,1,2,7,8,5,6/)
     if(boundary_type(ibound)== 4)ind_ref(1:8)=(/3,4,1,2,7,8,5,6/)
     if(boundary_type(ibound)== 5)ind_ref(1:8)=(/5,6,7,8,1,2,3,4/)
     if(boundary_type(ibound)== 6)ind_ref(1:8)=(/5,6,7,8,1,2,3,4/)
     ! Free boundary
     if(boundary_type(ibound)==11)ind_ref(1:8)=(/1,1,3,3,5,5,7,7/)
     if(boundary_type(ibound)==12)ind_ref(1:8)=(/2,2,4,4,6,6,8,8/)
     if(boundary_type(ibound)==13)ind_ref(1:8)=(/1,2,1,2,5,6,5,6/)
     if(boundary_type(ibound)==14)ind_ref(1:8)=(/3,4,3,4,7,8,7,8/)
     if(boundary_type(ibound)==15)ind_ref(1:8)=(/1,2,3,4,1,2,3,4/)
     if(boundary_type(ibound)==16)ind_ref(1:8)=(/5,6,7,8,5,6,7,8/)
     ! Imposed boundary (used only for flag1 and B perp)
     if(boundary_type(ibound)==21)ind_ref(1:8)=(/1,1,3,3,5,5,7,7/)
     if(boundary_type(ibound)==22)ind_ref(1:8)=(/2,2,4,4,6,6,8,8/)
     if(boundary_type(ibound)==23)ind_ref(1:8)=(/1,2,1,2,5,6,5,6/)
     if(boundary_type(ibound)==24)ind_ref(1:8)=(/3,4,3,4,7,8,7,8/)
     if(boundary_type(ibound)==25)ind_ref(1:8)=(/1,2,3,4,1,2,3,4/)
     if(boundary_type(ibound)==26)ind_ref(1:8)=(/5,6,7,8,5,6,7,8/)
     ! For magnetic field, we have only free boundary amd imposed boundary
     if(boundary_dir==1)alt(1:8)=-(/2,1,2,1,2,1,2,1/)
     if(boundary_dir==2)alt(1:8)=+(/1,2,1,2,1,2,1,2/)
     if(boundary_dir==3)alt(1:8)=-(/1,1,2,2,1,1,2,2/)
     if(boundary_dir==4)alt(1:8)=+(/2,2,1,1,2,2,1,1/)
     if(boundary_dir==5)alt(1:8)=-(/1,1,1,1,2,2,2,2/)
     if(boundary_dir==6)alt(1:8)=+(/2,2,2,2,1,1,1,1/)

     ! Velocity sign switch for reflexive boundary conditions
     gs=(/1,1,1/)
     if(boundary_type(ibound)==1.or.boundary_type(ibound)==2)gs(1)=-1
     if(boundary_type(ibound)==3.or.boundary_type(ibound)==4)gs(2)=-1
     if(boundary_type(ibound)==5.or.boundary_type(ibound)==6)gs(3)=-1
     
     ! Direction of gravity vector for hydrostatic equilibrium
     if(boundary_dir==1.or.boundary_dir==2)gdim=1
     if(boundary_dir==3.or.boundary_dir==4)gdim=2
     if(boundary_dir==5.or.boundary_dir==6)gdim=3

     if(boundary_dir==1)ind0(1:2,1:4)=RESHAPE((/2,4,6,8,1,3,5,7/),SHAPE=(/2, 4/))
     if(boundary_dir==2)ind0(1:2,1:4)=RESHAPE((/1,3,5,7,2,4,6,8/),SHAPE=(/2, 4/))
     if(boundary_dir==3)ind0(1:2,1:4)=RESHAPE((/3,4,7,8,1,2,5,6/),SHAPE=(/2, 4/))
     if(boundary_dir==4)ind0(1:2,1:4)=RESHAPE((/1,2,5,6,3,4,7,8/),SHAPE=(/2, 4/))
     if(boundary_dir==5)ind0(1:2,1:4)=RESHAPE((/5,6,7,8,1,2,3,4/),SHAPE=(/2, 4/))
     if(boundary_dir==6)ind0(1:2,1:4)=RESHAPE((/1,2,3,4,5,6,7,8/),SHAPE=(/2, 4/))

     if(boundary_dir==1)then
        iperp1=6; iperp2=nvar+1
     endif
     if(boundary_dir==2)then
        iperp1=nvar+1; iperp2=6
     endif
     if(boundary_dir==3)then
        iperp1=7; iperp2=nvar+2
     endif
     if(boundary_dir==4)then
        iperp1=nvar+2; iperp2=7
     endif
     if(boundary_dir==5)then
        iperp1=8; iperp2=nvar+3
     endif
     if(boundary_dir==6)then
        iperp1=nvar+3; iperp2=8
     endif

     ! Loop over grids by vector sweeps
     ncache=boundary(ibound,ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=boundary(ibound,ilevel)%igrid(igrid+i-1)
        end do
        
        ! Gather neighboring reference grid
        do i=1,ngrid
           ind_grid_ref(i)=son(nbor(ind_grid(i),inbor))
        end do

        ! Wall and free boundary conditions
        if((boundary_type(ibound)/10).ne.2)then
           
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              
              ! Gather neighboring reference cell
              iskip_ref=ncoarse+(ind_ref(ind)-1)*ngridmax
              do i=1,ngrid
                 ind_cell_ref(i)=iskip_ref+ind_grid_ref(i)
              end do
              
              ! Gather reference hydro variables
              do ivar=1,nvar+3
                 do i=1,ngrid
                    uu(i,ivar)=uold(ind_cell_ref(i),ivar)
                 end do
              end do
              ! Remove magnetic energy
              do i=1,ngrid
                 emag=0.125d0*((uu(i,6)+uu(i,nvar+1))**2+ &
                      &        (uu(i,7)+uu(i,nvar+2))**2+ &
                      &        (uu(i,8)+uu(i,nvar+3))**2)
                 uu(i,5)=uu(i,5)-emag
              end do

              ! Scatter to boundary region
              do ivar=1,nvar+3
                 switch=1
                 if(ivar>1.and.ivar<=4)switch=gs(ivar-1)
                 if(ivar.ne.(5+gdim).and.ivar.ne.(nvar+gdim))then
                    do i=1,ngrid
                       uold(ind_cell(i),ivar)=uu(i,ivar)*switch
                    end do
                 endif
                 if(ivar==(5+gdim))then
                    do i=1,ngrid
                       uold(ind_cell(i),5+gdim)=uu(i,5+gdim)+(uu(i,nvar+gdim)-uu(i,5+gdim))*alt(ind)
                    end do
                 endif
                 if(ivar==(nvar+gdim))then
                    do i=1,ngrid
                       uold(ind_cell(i),nvar+gdim)=uu(i,nvar+gdim)+(uu(i,nvar+gdim)-uu(i,5+gdim))*alt(ind)
                    end do
                 endif
              end do
              ! Add magnetic energy
              do i=1,ngrid
                 emag=0.125d0*((uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))**2+ &
                      &        (uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))**2+ &
                      &        (uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))**2)
                 uold(ind_cell(i),5)=uold(ind_cell(i),5)+emag
              end do

           end do
           ! End loop over cells

        ! Imposed boundary conditions
        else
              
           ! Loop over cells
           do ind=1,twotondim

              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              
              ! Compute cell center in code units
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
              
              call boundana(xx,uu,dx_loc,ibound,ngrid)
              
              ! Remove magnetic energy
              do i=1,ngrid
                 emag=0.125d0*((uu(i,6)+uu(i,nvar+1))**2+ &
                      &        (uu(i,7)+uu(i,nvar+2))**2+ &
                      &        (uu(i,8)+uu(i,nvar+3))**2)
                 uu(i,5)=uu(i,5)-emag
              end do

              ! Scatter variables
              do ivar=1,nvar+3
                 do i=1,ngrid
                    uold(ind_cell(i),ivar)=uu(i,ivar)
                 end do
              end do

           end do
           ! End loop over cells

           ! Correct normal magnetic field component

           ! Loop over planes
           do iplane=1,2

              ! Loop over cells
              do icell=1,twotondim/2

                 ind=ind0(iplane,icell)
                 
                 iskip=ncoarse+(ind-1)*ngridmax
                 do i=1,ngrid
                    ind_cell(i)=iskip+ind_grid(i)
                 end do
                 
                 ! Gather neighboring reference cell
                 if(iplane==1)then
                    iskip_ref=ncoarse+(ind_ref(ind)-1)*ngridmax
                    do i=1,ngrid
                       ind_cell_ref(i)=iskip_ref+ind_grid_ref(i)
                    end do
                 else
                    iskip_ref=ncoarse+(ind0(1,icell)-1)*ngridmax
                    do i=1,ngrid
                       ind_cell_ref(i)=iskip_ref+ind_grid(i)
                    end do
                 endif
                    
                 ! Compute correction
                 do i=1,ngrid
                    uu_ref(i,icell)=uold(ind_cell_ref(i),iperp1)-uold(ind_cell(i),iperp2)
                 end do
                 
                 do i=1,ngrid
                    uold(ind_cell(i),5+gdim)=uold(ind_cell(i),5+gdim)+uu_ref(i,icell)
                    uold(ind_cell(i),nvar+gdim)=uold(ind_cell(i),nvar+gdim)+uu_ref(i,icell)
                 end do
                 
              end do
              ! End loop over cells

           end do
           ! End loop over planes
                 
           ! Add magnetic energy

           ! Loop over cells
           do ind=1,twotondim

              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              
              do i=1,ngrid
                 emag=0.125d0*((uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))**2+ &
                      &        (uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))**2+ &
                      &        (uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))**2)
                 uold(ind_cell(i),5)=uold(ind_cell(i),5)+emag
              end do
              
           end do
           ! End loop over cells
           
        end if

     end do
     ! End loop over grids

  end do
  ! End loop over boundaries

111 format('   Entering make_boundary_hydro for level ',I2)

end subroutine make_boundary_hydro
