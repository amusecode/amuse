subroutine refine
  use amr_commons
  implicit none

  integer::ilevel

  if(verbose)write(*,*)'Entering refine' 

  call refine_coarse
  call build_comm(1)
  call make_virtual_fine_int(cpu_map(1),1)
  do ilevel=1,nlevelmax-1
     call refine_fine(ilevel)
     call build_comm(ilevel+1)
     call make_virtual_fine_int(cpu_map(1),ilevel+1)
  end do

  if(verbose)write(*,*)'Complete refine'

end subroutine refine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine refine_coarse
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::nxny,i,j,k
  integer::ind,info,ibound
  logical::boundary_region
  logical::ok_free,ok_all
  integer,dimension(1:nvector),save::ind_cell_tmp
  
  if(verbose)write(*,*)'  Entering refine_coarse'
  
  ! Constants
  nxny=nx*ny
  
  ! Compute cell authorization map
  call authorize_coarse

  ! Gather coarse cells for refinement
  ncreate=0
  do k=0,nz-1
  do j=0,ny-1
  do i=0,nx-1
     ind=1+i+j*nx+k*nxny
     if(flag2(ind)==1.and.flag1(ind)==1.and.son(ind)==0)then
        ncreate=ncreate+1
     end if
  end do
  end do
  end do

  ! Check for free memory
  ok_free=(numbf-ncreate)>0
  if(.not. ok_free)then
     write(*,*)'No more free memory'
     write(*,*)'Increase ngridmax'
#ifndef WITHOUTMPI
     call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
     stop
#endif
  end if

  ! Refine marked cells
  ibound=0; boundary_region=.false.
  do k=kcoarse_min,kcoarse_max
  do j=jcoarse_min,jcoarse_max
  do i=icoarse_min,icoarse_max
     ind=1+i+j*nx+k*nxny
     if(flag2(ind)==1.and.flag1(ind)==1.and.son(ind)==0)then
        call make_grid_coarse(ind,ibound,boundary_region)
     end if
  end do
  end do
  end do
  boundary_region=.true.
  do ibound=1,nboundary
     do k=kbound_min(ibound),kbound_max(ibound)
     do j=jbound_min(ibound),jbound_max(ibound)
     do i=ibound_min(ibound),ibound_max(ibound)
        ind=1+i+j*nx+k*nxny
        if(flag2(ind)==1.and.flag1(ind)==1.and.son(ind)==0)then
           call make_grid_coarse(ind,ibound,boundary_region)
        end if
     end do
     end do
     end do
  end do

  if(verbose)write(*,112)ncreate

  ! De-refine coarse cells marked for de-refinement
  nkill=0
  ibound=0; boundary_region=.false.
  do k=kcoarse_min,kcoarse_max
  do j=jcoarse_min,jcoarse_max
  do i=icoarse_min,icoarse_max
     ind=1+i+j*nx+k*nxny
     if(flag1(ind)==0.and.son(ind)>0)then
        nkill=nkill+1
        ind_cell_tmp(1)=ind
        call kill_grid(ind_cell_tmp,1,1,ibound,boundary_region)
     end if
  end do
  end do
  end do
  boundary_region=.true.
  do ibound=1,nboundary
     do k=kbound_min(ibound),kbound_max(ibound)
     do j=jbound_min(ibound),jbound_max(ibound)
     do i=ibound_min(ibound),ibound_max(ibound)
        ind=1+i+j*nx+k*nxny
        if(flag1(ind)==0.and.son(ind)>0)then
           nkill=nkill+1
           ind_cell_tmp(1)=ind
           call kill_grid(ind_cell_tmp,1,1,ibound,boundary_region)
        end if
     end do
     end do
     end do
  end do

  if(verbose)write(*,113)nkill

  ! Compute grid number statistics at level 1
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbl(myid,1),numbtot(1,1),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(numbl(myid,1),numbtot(2,1),1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(numbl(myid,1),numbtot(3,1),1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(used_mem,used_mem_tot,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbtot(1,1)=numbl(myid,1)
  numbtot(2,1)=numbl(myid,1)
  numbtot(3,1)=numbl(myid,1)
  used_mem_tot=used_mem
#endif
  numbtot(4,1)=numbtot(1,1)/ncpu

112 format('   ==> Make ',i6,' sub-grids')
113 format('   ==> Kill ',i6,' sub-grids')

end subroutine refine_coarse
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine make_grid_coarse(ind_cell,ibound,boundary_region)
  use amr_commons
  implicit none
  integer::ind_cell,ibound
  logical::boundary_region
  !----------------------------------------------------------
  ! This routine create one new grid at level 1, contained
  ! in father cells where ind_cell=1+ix+iy*nx+iz*nxny
  ! is the actual cell number of the father coarse cell.
  !----------------------------------------------------------
  integer::j,igrid,nxny,iskip,icpu,nx_loc
  integer::ix,iy,iz,ind_grid_son,idim
  integer,dimension(1:twondim)::ixn,iyn,izn

  real(dp)::dx_loc,scale
  real(dp),dimension(1:3)::xc,skip_loc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  integer ,dimension(1:nvector),save::cc

  ! Local constants
  nxny=nx*ny
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=scale

  ! Get one new grid from free memory
  igrid=headf
  ind_grid_son=igrid
  headf=next(headf)
  numbf=numbf-1
  used_mem=ngridmax-numbf

  ! Compute grid center
  iz=(ind_cell-1)/nxny
  iy=(ind_cell-1-iz*nxny)/nx
  ix=(ind_cell-1-iz*nxny-iy*nx)
  if(ndim>0)xg(ind_grid_son,1)=(dble(ix)+0.5D0)
  if(ndim>1)xg(ind_grid_son,2)=(dble(iy)+0.5D0)
  if(ndim>2)xg(ind_grid_son,3)=(dble(iz)+0.5D0)

  ! Connect to father cell
  son(ind_cell)=ind_grid_son
  father(ind_grid_son)=ind_cell
  
  ! Connect to neighboring father cells
  do j=1,twondim
     ixn(j)=ix
     iyn(j)=iy
     izn(j)=iz
  end do

  ! With periodic boundary conditions
  if(ndim>0)then
     if(ix>0)then
        ixn(1)=ix-1
     else
        ixn(1)=nx-1
     end if
     if(ix<nx-1)then
        ixn(2)=ix+1
     else
        ixn(2)=0
     end if
  end if
#if NDIM>1
  if(ndim>1)then
     if(iy>0)then
        iyn(3)=iy-1
     else
        iyn(3)=ny-1
     end if
     if(iy<ny-1)then
        iyn(4)=iy+1
     else
        iyn(4)=0
     end if
  end if
#endif
#if NDIM>2
  if(ndim>2)then
     if(iz>0)then
        izn(5)=iz-1
     else
        izn(5)=nz-1
     end if
     if(iz<nz-1)then
        izn(6)=iz+1
     else
        izn(6)=0
     end if
  end if
#endif
  do j=1,twondim
     nbor(ind_grid_son,j)=1+ixn(j)+iyn(j)*nx+izn(j)*nxny
  end do

  ! Update cpu map
  if(boundary_region)then
     do j=1,twotondim
        iskip=ncoarse+(j-1)*ngridmax
        cpu_map(iskip+ind_grid_son)=0
     end do
  else
     do j=1,twotondim
        iskip=ncoarse+(j-1)*ngridmax
        iz=(j-1)/4
        iy=(j-1-4*iz)/2
        ix=(j-1-2*iy-4*iz)
        if(ndim>0)xc(1)=(dble(ix)-0.5D0)/2.0d0
        if(ndim>1)xc(2)=(dble(iy)-0.5D0)/2.0d0
        if(ndim>2)xc(3)=(dble(iz)-0.5D0)/2.0d0
        ! Compute new cell position
        do idim=1,ndim
           xx(1,idim)=xg(ind_grid_son,idim)+xc(idim)
        end do
        ! Rescale position from code units to user units
        do idim=1,ndim
           xx(1,idim)=(xx(1,idim)-skip_loc(idim))*scale
        end do
        call cmp_cpumap(xx,cc,1)
        cpu_map(iskip+ind_grid_son)=cc(1)
     end do
  end if

  ! Connect new grid to level 1 grids linked list
  if(boundary_region)then
     igrid=ind_grid_son
     if(numbb(ibound,1)>0)then
        next(igrid)=0
        prev(igrid)=tailb(ibound,1)
        next(tailb(ibound,1))=igrid
        tailb(ibound,1)=igrid
        numbb(ibound,1)=numbb(ibound,1)+1
     else
        next(igrid)=0
        prev(igrid)=0
        headb(ibound,1)=igrid
        tailb(ibound,1)=igrid
        numbb(ibound,1)=1
     end if
  else
     igrid=ind_grid_son
     icpu=cpu_map(ind_cell)
     if(numbl(icpu,1)>0)then
        next(igrid)=0
        prev(igrid)=taill(icpu,1)
        next(taill(icpu,1))=igrid
        taill(icpu,1)=igrid
        numbl(icpu,1)=numbl(icpu,1)+1
     else
        next(igrid)=0
        prev(igrid)=0
        headl(icpu,1)=igrid
        taill(icpu,1)=igrid
        numbl(icpu,1)=1
     end if
  end if

end subroutine make_grid_coarse
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine refine_fine(ilevel)
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !---------------------------------------------------------
  ! This routine refines cells at level ilevel if cells
  ! are flagged for refinement and are not already refined.
  ! This routine destroys refinements for cells that are 
  ! not flagged for refinement and are already refined.
  ! For single time-stepping, numerical rules are 
  ! automatically satisfied. For adaptive time-stepping,
  ! numerical rules are checked before refining any cell.
  !---------------------------------------------------------
  integer::ncache,ngrid
  integer::igrid,icell,i
  integer::ind,iskip,info,icpu,ibound
  integer::ncreate_tmp,nkill_tmp
  logical::boundary_region
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector),save::ind_grid_tmp,ind_cell_tmp
  logical,dimension(1:nvector),save::ok
  logical::ok_free,ok_all

  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  !--------------------------
  ! Compute authorization map
  !--------------------------
  call authorize_fine(ilevel)

  if(.not. shrink)then
  !---------------------------------------------------
  ! Step 1: if cell is flagged for refinement and
  ! if it is not already refined, create a son grid.
  !---------------------------------------------------
        
  !------------------------------------
  ! Refine cells marked for refinement
  !------------------------------------
  ncreate=0
  do icpu=1,ncpu+nboundary  ! Loop over cpus and boundaries
     if(icpu==myid)then
        ibound=0
        boundary_region=.false.
        ncache=active(ilevel)%ngrid
     else if(icpu<=ncpu)then
        ibound=0
        boundary_region=.false.
        ncache=reception(icpu,ilevel)%ngrid
     else
        ibound=icpu-ncpu
        boundary_region=.true.
        ncache=boundary(ibound,ilevel)%ngrid
     end if
     do igrid=1,ncache,nvector  ! Loop over grids
        ngrid=MIN(nvector,ncache-igrid+1)
        if(myid==icpu)then
           do i=1,ngrid
              ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
           end do
        else if(icpu<=ncpu)then
           do i=1,ngrid
              ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
           end do
        else
           do i=1,ngrid
              ind_grid(i)=boundary(ibound,ilevel)%igrid(igrid+i-1)           
           end do
        end if
        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           ! Gather flagged, unrefined and authorized cells
           do i=1,ngrid
              ok(i)= flag2(ind_cell(i))==1 .and. &
                   & flag1(ind_cell(i))==1 .and. &
                   & son  (ind_cell(i))==0
           end do
           ! Count cells for refinement
           ncreate_tmp=0
           do i=1,ngrid
              if(ok(i))ncreate_tmp=ncreate_tmp+1
           end do
           ncreate=ncreate+ncreate_tmp

           ! Check for free memory
           if(ncreate_tmp>=numbf) then
              write(*,*)'No more free memory'
              write(*,*)'Increase ngridmax'
#ifndef WITHOUTMPI
              call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
              stop
#endif
           end if

           ! Refine selected cells
           if(ncreate_tmp>0)then
              icell=0
              do i=1,ngrid
                 if(ok(i))then
                    icell=icell+1
                    ind_grid_tmp(icell)=ind_grid(i)
                    ind_cell_tmp(icell)=ind_cell(i)
                 end if
              end do
              call make_grid_fine(ind_grid_tmp,ind_cell_tmp,ind, &
                   & ilevel+1,ncreate_tmp,ibound,boundary_region)
           end if
        end do
     end do
  end do
  if(verbose)write(*,112)ncreate
  endif

  !-----------------------------------------------------
  ! Case 2: if cell is not flagged for refinement,but
  ! it is refined, then destroy its child grid.
  !-----------------------------------------------------
  nkill=0  
  do icpu=1,ncpu+nboundary  ! Loop over cpus and boundaries
     if(icpu==myid)then
        ibound=0
        boundary_region=.false.
        ncache=active(ilevel)%ngrid
     else if(icpu<=ncpu)then
        ibound=0
        boundary_region=.false.
        ncache=reception(icpu,ilevel)%ngrid
     else
        ibound=icpu-ncpu
        boundary_region=.true.
        ncache=boundary(ibound,ilevel)%ngrid
     end if
     do igrid=1,ncache,nvector  ! Loop over grids
        ngrid=MIN(nvector,ncache-igrid+1)
        if(myid==icpu)then
           do i=1,ngrid
              ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
           end do
        else if(icpu<=ncpu)then
           do i=1,ngrid
              ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
           end do
        else
           do i=1,ngrid
              ind_grid(i)=boundary(ibound,ilevel)%igrid(igrid+i-1)           
           end do
        end if
        do ind=1,twotondim     ! Loop over cells
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           if(shrink)then
              ! Gather unauthorized and refined cells
              do i=1,ngrid
                 ok(i)= flag2(ind_cell(i))==0 .and. &
                      & son  (ind_cell(i))>0
              end do
           else
              ! Gather unflagged and refined cells
              do i=1,ngrid
                 ok(i)= flag1(ind_cell(i))==0 .and. &
                      & son  (ind_cell(i))>0
              end do
           endif
           ! Count cells for de-refinement
           nkill_tmp=0
           do i=1,ngrid
              if(ok(i))then
                 nkill_tmp=nkill_tmp+1
              end if
           end do
           nkill=nkill+nkill_tmp
           ! De-refine selected cells
           if(nkill_tmp>0)then
              icell=0
              do i=1,ngrid
                 if(ok(i))then
                    icell=icell+1
                    ind_cell_tmp(icell)=ind_cell(i)
                 end if
              end do
              call kill_grid(ind_cell_tmp,ilevel+1,nkill_tmp,ibound,boundary_region)
           end if           
        end do  ! End loop over cells
     end do
  end do
  if(verbose)write(*,113)nkill

  ! Compute grid number statistics at level ilevel+1
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbl(myid,ilevel+1),numbtot(1,ilevel+1),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(numbl(myid,ilevel+1),numbtot(2,ilevel+1),1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(numbl(myid,ilevel+1),numbtot(3,ilevel+1),1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(used_mem,used_mem_tot,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbtot(1,ilevel+1)=numbl(myid,ilevel+1)
  numbtot(2,ilevel+1)=numbl(myid,ilevel+1)
  numbtot(3,ilevel+1)=numbl(myid,ilevel+1)
  used_mem_tot=used_mem
#endif
  numbtot(4,ilevel+1)=numbtot(1,ilevel+1)/ncpu

111 format('   Entering refine_fine for level ',I2)
112 format('   ==> Make ',i6,' sub-grids')
113 format('   ==> Kill ',i6,' sub-grids')

end subroutine refine_fine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine make_grid_fine(ind_grid,ind_cell,ind,ilevel,nn,ibound,boundary_region)
  use amr_commons
  use hydro_commons
  use poisson_commons, ONLY:f, phi,phi_old
#ifdef RT
  use rt_hydro_commons
#endif
#ifdef ATON
  use radiation_commons, ONLY:Erad
#endif
  implicit none
  integer::nn,ind,ilevel,ibound
  logical::boundary_region
  integer,dimension(1:nvector)::ind_grid,ind_cell
  !--------------------------------------------------------------
  ! This routine create new grids at level ilevel (ilevel >= 2)
  ! contained in father cell ind_cell(:).
  ! ind_grid(:) is the number of the grid that contains 
  ! the father cell and ind = 1, 2, 3, 4, 5, 6, 7, or 8.
  ! The actual father cell number is:
  ! ind_cell = ncoarse + (ind-1)*ngridmax + ind_grid
  ! WARNING: USE THIS ROUTINE WITH CARE, SINCE IT ASSUMES THAT
  ! ALL FATHER CELL'S NEIGHBORS DO EXIST !!!
  !--------------------------------------------------------------
  integer ::idim,igrid,iskip,icpu
  integer ::i,j,ix,iy,iz,ivar,nx_loc
  integer ,dimension(1:nvector)          ,save::ind_grid_son
  integer ,dimension(1:nvector,0:twondim),save::ind_fathers
  integer ,dimension(1:nvector,0:twondim),save::igridn
  integer ,dimension(1:nvector,1:twondim),save::indn

  real(dp)::dx,dx_loc,scale
  real(dp),dimension(1:3)::xc,skip_loc
#ifdef SOLVERmhd
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3),save::u2
  integer ,dimension(1:nvector,0:twondim),save::ind1
#else
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar),save::u2
#endif
#ifdef RT
  real(dp),dimension(1:nvector,0:twondim  ,1:nrtvar),save::urt1
  real(dp),dimension(1:nvector,1:twotondim,1:nrtvar),save::urt2
#endif  
  real(dp),dimension(1:nvector,0:twondim  ,1:ndim),save::g1=0.0
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::g2=0.0

  real(dp),dimension(1:nvector,1:ndim),save::xx
  integer ,dimension(1:nvector),save::cc

  logical::error

  ! Mesh spacing in father level
  dx=0.5D0**(ilevel-1)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Get nn new grids from free memory
  do i=1,nn
     igrid=headf
     ind_grid_son(i)=igrid
     headf=next(headf)
     numbf=numbf-1
     used_mem=ngridmax-numbf
  end do
  
  ! Set new grids position
  iz=(ind-1)/4
  iy=(ind-1-4*iz)/2
  ix=(ind-1-2*iy-4*iz)
  if(ndim>0)xc(1)=(dble(ix)-0.5D0)*dx
  if(ndim>1)xc(2)=(dble(iy)-0.5D0)*dx
  if(ndim>2)xc(3)=(dble(iz)-0.5D0)*dx
  do idim=1,ndim
     do i=1,nn
        xg(ind_grid_son(i),idim)=xg(ind_grid(i),idim)+xc(idim)
     end do
  end do

  ! Connect new grids to father cells
  do i=1,nn
     son(ind_cell(i))=ind_grid_son(i)
  end do
  do i=1,nn
     father(ind_grid_son(i))=ind_cell(i)
  end do

  ! Connect news grids to neighboring father cells
  call getnborgrids(ind_grid,igridn,nn)
  call getnborcells(igridn,ind,indn,nn)
  error=.false.
  do j=1,twondim
     do i=1,nn
        nbor(ind_grid_son(i),j)=indn(i,j)
        if(indn(i,j)==0)then
           error=.true.
        end if
     end do
  end do
  if(error)then
     do j=1,twondim
        do i=1,nn
           if(indn(i,j)==0)then
              if(cpu_map(ind_cell(i))==myid)then
                 write(*,*)'Fatal error in make_grid_fine'
                 write(*,*)myid,cpu_map(ind_cell(i))
                 write(*,*)ilevel,j,ibound,boundary_region
                 stop
              endif
           end if
        end do
     end do
  end if
  
  ! Update cpu map
  if(boundary_region)then
     do j=1,twotondim
        iskip=ncoarse+(j-1)*ngridmax
        do i=1,nn
           cpu_map(iskip+ind_grid_son(i))=0
        end do
     end do
  else
     do j=1,twotondim
        iz=(j-1)/4
        iy=(j-1-4*iz)/2
        ix=(j-1-2*iy-4*iz)
        if(ndim>0)xc(1)=(dble(ix)-0.5D0)*dx/2.0d0
        if(ndim>1)xc(2)=(dble(iy)-0.5D0)*dx/2.0d0
        if(ndim>2)xc(3)=(dble(iz)-0.5D0)*dx/2.0d0
        ! Compute cell coordinates
        do idim=1,ndim
           do i=1,nn
              xx(i,idim)=xg(ind_grid_son(i),idim)+xc(idim)
           end do
        end do
        ! Rescale position from code units to user units
        do idim=1,ndim
           do i=1,nn
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
           end do
        end do
        call cmp_cpumap(xx,cc,nn)
        iskip=ncoarse+(j-1)*ngridmax
        do i=1,nn
           cpu_map(iskip+ind_grid_son(i))=cc(i)
        end do
     end do
  end if

  ! Connect news grids to level ilevel linked list
  if(boundary_region)then
     do i=1,nn
        igrid=ind_grid_son(i)
        if(numbb(ibound,ilevel)>0)then
           next(igrid)=0
           prev(igrid)=tailb(ibound,ilevel)
           next(tailb(ibound,ilevel))=igrid
           tailb(ibound,ilevel)=igrid
           numbb(ibound,ilevel)=numbb(ibound,ilevel)+1
        else
           next(igrid)=0
           prev(igrid)=0
           headb(ibound,ilevel)=igrid
           tailb(ibound,ilevel)=igrid
           numbb(ibound,ilevel)=1
        end if
     end do
  else
     do i=1,nn
        igrid=ind_grid_son(i)
        icpu=cpu_map(ind_cell(i))
        if(numbl(icpu,ilevel)>0)then
           next(igrid)=0
           prev(igrid)=taill(icpu,ilevel)
           next(taill(icpu,ilevel))=igrid
           taill(icpu,ilevel)=igrid
           numbl(icpu,ilevel)=numbl(icpu,ilevel)+1
        else
           next(igrid)=0
           prev(igrid)=0
           headl(icpu,ilevel)=igrid
           taill(icpu,ilevel)=igrid
           numbl(icpu,ilevel)=1
        end if
     end do
  end if

  ! Interpolate parent variables to get new children ones
  if(.not.init .and. .not.balance)then
     ! Get neighboring father cells
     do i=1,nn
        ind_fathers(i,0)=father(ind_grid_son(i))
     end do
     do j=1,twondim
        do i=1,nn
           ind_fathers(i,j)=nbor(ind_grid_son(i),j)
        end do
     end do
     !============================
     ! Interpolate hydro variables
     !============================
     if(hydro)then
        do j=0,twondim
           ! Gather hydro variables
#ifdef SOLVERmhd
           do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                 do i=1,nn
                    u1(i,j,ivar)=uold(ind_fathers(i,j),ivar)
                 end do
#ifdef SOLVERmhd
              end do
#else
           end do
#endif
#ifdef SOLVERmhd
           ! Gather son index
           do i=1,nn
              ind1(i,j)=son(ind_fathers(i,j))
           end do
#endif
        end do
        ! Interpolate
#ifdef SOLVERmhd
        call interpol_hydro(u1,ind1,u2,nn)
#else
        call interpol_hydro(u1,u2,nn)
#endif
        ! Scatter to children cells
        do j=1,twotondim
           iskip=ncoarse+(j-1)*ngridmax
#ifdef SOLVERmhd
           do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                 do i=1,nn
                    uold(iskip+ind_grid_son(i),ivar)=u2(i,j,ivar)
                 end do
#ifdef SOLVERmhd
              end do
#else
           end do
#endif
        enddo
     end if
#ifdef RT
     !============================
     ! Interpolate RT variables
     !============================
     if(rt)then
        do j=0,twondim
           ! Gather hydro variables
           do ivar=1,nrtvar
              do i=1,nn
                 urt1(i,j,ivar)=rtuold(ind_fathers(i,j),ivar)
              end do
           end do
        end do
        ! Interpolate
        call rt_interpol_hydro(urt1,urt2,nn)
        ! Scatter to children cells
        do j=1,twotondim
           iskip=ncoarse+(j-1)*ngridmax
           do ivar=1,nrtvar
              do i=1,nn
                 rtuold(iskip+ind_grid_son(i),ivar)=urt2(i,j,ivar)
              end do
           end do
        enddo
     end if
#endif
     !==============================
     ! Interpolate gravity variables
     !==============================
     if(poisson)then
        ! Scatter to children cells
        do j=1,twotondim
           iskip=ncoarse+(j-1)*ngridmax
           do idim=1,ndim
              do i=1,nn
                 f(iskip+ind_grid_son(i),idim)=f(ind_fathers(i,0),idim)
              end do
           end do
           do i=1,nn
              phi(iskip+ind_grid_son(i))=phi(ind_fathers(i,0))
              phi_old(iskip+ind_grid_son(i))=phi_old(ind_fathers(i,0))
           end do
        end do
     end if
     !===========================
     ! Interpolate ATON variables
     !===========================
#ifdef ATON
     if(aton)then
        do j=1,twotondim
           iskip=ncoarse+(j-1)*ngridmax
           do i=1,nn
              Erad(iskip+ind_grid_son(i))=Erad(ind_fathers(i,0))
           end do
        enddo
     end if
#endif
  endif

end subroutine make_grid_fine
!###############################################################
!###############################################################
!###############################################################
!###############################################################
subroutine kill_grid(ind_cell,ilevel,nn,ibound,boundary_region)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
#ifdef RT
  use rt_hydro_commons
  use rt_parameters
#endif
#ifdef ATON
  use radiation_commons, ONLY:Erad
#endif
  implicit none
  integer::nn,ilevel,ibound
  logical::boundary_region
  integer,dimension(1:nvector)::ind_cell
  !----------------------------------------------------
  ! This routine destroy the grids at level ilevel
  ! contained in father cell ind_cell(:)
  !----------------------------------------------------
  integer::igrid,iskip,icpu
  integer::i,j,idim,ind,ivar
  integer,dimension(1:nvector),save::ind_grid_son,ind_cell_son
#ifdef RT
  real(dp),dimension(nIons)::xion
#endif

#ifdef RT
  if(upload_equilibrium_x) then                                       
     ! Enforce equilibrium on ionization states when merging, to      
     ! prevent unnatural values (e.g when merging hot and cold cells).
     do i=1,nn                                                        
        call calc_equilibrium_xion(uold(ind_cell(i),1:nvar) &
             , rtuold(ind_cell(i),1:nrtvar), xion)    
        uold(ind_cell(i),iIons:iIons+nIons-1)=xion*uold(ind_cell(i),1)
     enddo                                                            
  endif                                                               
#endif

  ! Gather son grids
  do i=1,nn
     ind_grid_son(i)=son(ind_cell(i))
  end do
  
  ! Disconnect son grids from father cells
  do i=1,nn
     son(ind_cell(i))=0
  end do

  ! Disconnect son grids from level ilevel linked list
  if(boundary_region)then
     do i=1,nn
        igrid=ind_grid_son(i)
        if(prev(igrid).ne.0) then
           if(next(igrid).ne.0)then
              next(prev(igrid))=next(igrid)
              prev(next(igrid))=prev(igrid)
           else
              next(prev(igrid))=0
              tailb(ibound,ilevel)=prev(igrid)
           end if
        else
           if(next(igrid).ne.0)then
              prev(next(igrid))=0
              headb(ibound,ilevel)=next(igrid)
           else
              headb(ibound,ilevel)=0
              tailb(ibound,ilevel)=0
           end if
        end if
        numbb(ibound,ilevel)=numbb(ibound,ilevel)-1 
     end do
  else
     do i=1,nn
        igrid=ind_grid_son(i)
        icpu=cpu_map(ind_cell(i))
        if(prev(igrid).ne.0) then
           if(next(igrid).ne.0)then
              next(prev(igrid))=next(igrid)
              prev(next(igrid))=prev(igrid)
           else
              next(prev(igrid))=0
              taill(icpu,ilevel)=prev(igrid)
           end if
        else
           if(next(igrid).ne.0)then
              prev(next(igrid))=0
              headl(icpu,ilevel)=next(igrid)
           else
              headl(icpu,ilevel)=0
              taill(icpu,ilevel)=0
           end if
        end if
        numbl(icpu,ilevel)=numbl(icpu,ilevel)-1 
     end do
  end if

  ! Reset grid variables
  do idim=1,ndim
     do i=1,nn
        xg(ind_grid_son(i),idim)=0.0D0
     end do
  end do
  do i=1,nn
     father(ind_grid_son(i))=0
  end do
  do j=1,twondim
     do i=1,nn
        nbor(ind_grid_son(i),j)=0
     end do
  end do
  if(pic)then
     do i=1,nn
        headp(ind_grid_son(i))=0
        tailp(ind_grid_son(i))=0
        numbp(ind_grid_son(i))=0
     end do
  end if

  ! Reset cell variables
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,nn
        ind_cell_son(i)=iskip+ind_grid_son(i)
     end do
     ! Tree variables
     do i=1,nn
        son     (ind_cell_son(i))=0
        flag1   (ind_cell_son(i))=0
        flag2   (ind_cell_son(i))=0
        cpu_map (ind_cell_son(i))=0
        cpu_map2(ind_cell_son(i))=0
     end do
     ! Gravity variables
     if(poisson)then
        do i=1,nn
           rho(ind_cell_son(i))=0.0D0
           phi(ind_cell_son(i))=0.0D0
           phi_old(ind_cell_son(i))=0.0D0
        end do
        do idim=1,ndim
           do i=1,nn
              f(ind_cell_son(i),idim)=0.0D0
           end do
        end do
     end if
     ! Hydro variables
     if(hydro)then
#ifdef SOLVERmhd
        do ivar=1,nvar+3
#else
        do ivar=1,nvar
#endif
           do i=1,nn
              uold(ind_cell_son(i),ivar)=0.0D0
              unew(ind_cell_son(i),ivar)=0.0D0
           end do
#ifdef SOLVERmhd
        end do
#else
        end do
#endif
     end if
#ifdef RT
     ! RT variables
     if(rt)then
        do ivar=1,nrtvar
           do i=1,nn
              rtuold(ind_cell_son(i),ivar)=0.0D0
              rtunew(ind_cell_son(i),ivar)=0.0D0
           end do
        end do
     end if
#endif
#ifdef ATON
     if(aton)then
        do i=1,nn
           Erad(ind_cell_son(i))=0.0D0
        end do
     end if
#endif
  end do

  ! Put son grids at the tail of the free memory linked list
  do i=1,nn
     igrid=ind_grid_son(i)
     next(tailf)=igrid
     prev(igrid)=tailf
     next(igrid)=0
     tailf=igrid
     numbf=numbf+1
  end do
  
end subroutine kill_grid
