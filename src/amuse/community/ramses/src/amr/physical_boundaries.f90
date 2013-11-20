!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_boundary_coarse
  use amr_commons
  implicit none
  !----------------------------------------------------------------------
  ! This routine makes one cell width cubic buffer around inner cells 
  ! at the coarse level by performing a dilatation on array flag2.
  ! Array flag1 is used as temorary work space.
  ! This routine is called by subroutine authorize_coarse.
  !----------------------------------------------------------------------
  integer::nxny,nxm1,nym1,nzm1
  integer::ismooth,ind,i,j,k,ibound,i_nbor
  integer::ind1,ind2,i1,i2,j1,j2,k1,k2
  integer,dimension(1:3)::n_nbor

  if(verbose)write(*,*)'  Entering init_virtual_coarse'

  ! Constants
  nxny=nx*ny; nxm1=nx-1; nym1=ny-1; nzm1=nz-1

  ! Initialize flag2 to zero for physical boundaries
  do ibound=1,nboundary
     do k=kbound_min(ibound),kbound_max(ibound)
     do j=jbound_min(ibound),jbound_max(ibound)
     do i=ibound_min(ibound),ibound_max(ibound)
        ind=1+i+j*nx+k*nxny
        flag2(ind)=0
     end do
     end do
     end do
  end do

  ! Loop over steps
  n_nbor(1:3)=(/1,2,3/)
  flag2(0)=0
  do ismooth=1,ndim  
     ! Initialize flag1 to zero
     do ibound=1,nboundary
        do k=kbound_min(ibound),kbound_max(ibound)
        do j=jbound_min(ibound),jbound_max(ibound)
        do i=ibound_min(ibound),ibound_max(ibound)
           ind=1+i+j*nx+k*nxny
           flag1(ind)=0
        end do
        end do
        end do
     end do
     ! Count neighbors with flag2=1 and flag1 cells accordingly.
     do ibound=1,nboundary
        do k=kbound_min(ibound),kbound_max(ibound)
        do j=jbound_min(ibound),jbound_max(ibound)
        do i=ibound_min(ibound),ibound_max(ibound)
           ind=1+i+j*nx+k*nxny
           i_nbor=0
           if(ndim>0)then
              i1=i-1; i2=i+1
              if(i==0)i1=nxm1
              if(i==nxm1)i2=0
              ind1=1+i1+j*nx+k*nxny
              ind2=1+i2+j*nx+k*nxny
              i_nbor=i_nbor+flag2(ind1)+flag2(ind2)
           end if
           if(ndim>1)then
              j1=j-1; j2=j+1
              if(j==0)j1=nym1
              if(j==nym1)j2=0
              ind1=1+i+j1*nx+k*nxny
              ind2=1+i+j2*nx+k*nxny
              i_nbor=i_nbor+flag2(ind1)+flag2(ind2)
           end if
           if(ndim>2)then
              k1=k-1; k2=k+1
              if(k==0)k1=nzm1
              if(k==nzm1)k2=0
              ind1=1+i+j*nx+k1*nxny
              ind2=1+i+j*nx+k2*nxny
              i_nbor=i_nbor+flag2(ind1)+flag2(ind2)
           end if
           if(i_nbor>=n_nbor(ismooth))flag1(ind)=1
        end do
        end do
        end do
     end do
     ! Set flag2 to 1 for cells with flag1=1
     do ibound=1,nboundary
        do k=kbound_min(ibound),kbound_max(ibound)
        do j=jbound_min(ibound),jbound_max(ibound)
        do i=ibound_min(ibound),ibound_max(ibound)
           ind=1+i+j*nx+k*nxny
           if(flag1(ind)==1)flag2(ind)=1
        end do
        end do
        end do
     end do
  end do
  ! End loop over steps

end subroutine init_boundary_coarse
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_boundary_fine(ilevel)
  use amr_commons
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine makes one cell width cubic buffer around inner cells 
  ! at level ilevel by performing a dilatation on array flag2.
  ! This routine is called by subroutine authorize_fine.
  ! Array flag1 is used as temporary work space.
  ! -------------------------------------------------------------------
  integer::ismooth,ibound
  integer,dimension(1:3)::n_nbor
  integer::i,ncache,iskip
  integer::igrid,ind,ngrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector,0:twondim),save::igridn

  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Initialize flag2 to 0 for physical boundaries
  do ibound=1,nboundary
     ncache=boundary(ibound,ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=boundary(ibound,ilevel)%igrid(igrid+i-1)
        end do
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ngrid
              flag2(ind_cell(i))=0
           end do
        end do
     end do
  end do
  flag2(0)=0

  ! Loop over steps
  n_nbor(1:3)=(/1,2,3/)
  do ismooth=1,ndim
     ! Initialize flag1 to 0 for physical boundaries
     do ibound=1,nboundary
        ncache=boundary(ibound,ilevel)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
              ind_grid(i)=boundary(ibound,ilevel)%igrid(igrid+i-1)
           end do
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              do i=1,ngrid
                 flag1(ind_cell(i))=0
              end do
           end do
        end do
     end do
     ! Count neighbors and set flag1 accordingly
     do ibound=1,nboundary
        ncache=boundary(ibound,ilevel)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
              ind_grid(i)=boundary(ibound,ilevel)%igrid(igrid+i-1)
           end do
           call getnborgrids_check(ind_grid,igridn,ngrid)
           do ind=1,twotondim
              call count_nbors2(igridn,ind,n_nbor(ismooth),ngrid)
           end do
        end do
     end do
     ! Set flag2=1 for cells with flag1=1
     do ibound=1,nboundary
        ncache=boundary(ibound,ilevel)%ngrid
        do igrid=1,ncache,nvector
           ngrid=MIN(nvector,ncache-igrid+1)
           do i=1,ngrid
              ind_grid(i)=boundary(ibound,ilevel)%igrid(igrid+i-1)
           end do
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=iskip+ind_grid(i)
              end do
              do i=1,ngrid
                 if(flag1(ind_cell(i))==1)flag2(ind_cell(i))=1
              end do
           end do
        end do
     end do
  end do
  ! End loop over steps

111 format('   Entering init_boundary_fine for level ',I2)

end subroutine init_boundary_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_boundary_coarse
  use amr_commons
  use hydro_commons
  implicit none
  !---------------------------------------------------------------------
  ! This routine set up boundary conditions for flag at the coarse level
  !---------------------------------------------------------------------
  integer::nxny,boundary_dir
  integer::i,j,k,ibound,iskip,ind,ind_ref

  if(verbose)write(*,*)'  Entering make_boundary_coarse'

  ! Constants
  nxny=nx*ny

  ! Loop over boundaries
  do ibound=1,nboundary

     ! Compute offset depending on boundary direction
     boundary_dir=boundary_type(ibound)-10*(boundary_type(ibound)/10)
     if(boundary_dir==1)iskip=+1
     if(boundary_dir==2)iskip=-1
     if(boundary_dir==3)iskip=+nx
     if(boundary_dir==4)iskip=-nx
     if(boundary_dir==5)iskip=+nxny
     if(boundary_dir==6)iskip=-nxny

     ! flag1 boundary cell as reference cell
     do k=kbound_min(ibound),kbound_max(ibound)
        do j=jbound_min(ibound),jbound_max(ibound)
           do i=ibound_min(ibound),ibound_max(ibound)
              ind=1+i+j*nx+k*nxny
              ind_ref=1+i+j*nx+k*nxny+iskip
              flag1(ind)=flag1(ind_ref)
           end do
        end do
     end do

  end do
  ! End loop over boundaries

end subroutine make_boundary_coarse
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_boundary_flag(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_parameters
  implicit none
  integer::ilevel
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------
  integer::ibound,boundary_dir,idim,inbor
  integer::i,ncache,ivar,igrid,ngrid,ind
  integer::iskip,iskip_ref,gdim,nx_loc,ix,iy,iz
  integer,dimension(1:8)::ind_ref,alt
  integer,dimension(1:nvector),save::ind_grid,ind_grid_ref
  integer,dimension(1:nvector),save::ind_cell,ind_cell_ref
  integer,dimension(1:nvector),save::fff

  if(.not. simple_boundary)return
  if(verbose)write(*,111)ilevel

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
     ! Imposed boundary (used only for flag1)
     if(boundary_type(ibound)==21)ind_ref(1:8)=(/1,1,3,3,5,5,7,7/)
     if(boundary_type(ibound)==22)ind_ref(1:8)=(/2,2,4,4,6,6,8,8/)
     if(boundary_type(ibound)==23)ind_ref(1:8)=(/1,2,1,2,5,6,5,6/)
     if(boundary_type(ibound)==24)ind_ref(1:8)=(/3,4,3,4,7,8,7,8/)
     if(boundary_type(ibound)==25)ind_ref(1:8)=(/1,2,3,4,1,2,3,4/)
     if(boundary_type(ibound)==26)ind_ref(1:8)=(/5,6,7,8,5,6,7,8/)

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

           do i=1,ngrid
              fff(i)=flag1(ind_cell_ref(i))
           end do
           do i=1,ngrid
              flag1(ind_cell(i))=fff(i)
           end do
           
        end do
        ! End loop over cells

     end do
     ! End loop over grids

  end do
  ! End loop over boundaries

111 format('   Entering make_boundary_flag for level ',I2)

end subroutine make_boundary_flag


