!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_tree
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------
  ! This subroutine build the particle linked list at the 
  ! coarse level for ALL the particles in the box.
  ! This routine should be used only as initial set up for
  ! the particle tree.
  !------------------------------------------------------
  integer::ipart,idim,i,nxny,ilevel
  integer::npart1,info,icpu,nx_loc
  logical::error
  real(dp),dimension(1:3)::xbound
  integer,dimension(1:nvector),save::ix,iy,iz
  integer,dimension(1:nvector),save::ind_grid,ind_part
  logical,dimension(1:nvector),save::ok=.true.
  real(dp),dimension(1:3)::skip_loc
  real(dp)::scale

  if(verbose)write(*,*)'  Entering init_tree'

  ! Local constants
  nxny=nx*ny
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  !----------------------------------
  ! Initialize particle linked list
  !----------------------------------
  prevp(1)=0; nextp(1)=2
  do ipart=2,npartmax-1
     prevp(ipart)=ipart-1
     nextp(ipart)=ipart+1
  end do
  prevp(npartmax)=npartmax-1; nextp(npartmax)=0
  ! Free memory linked list
  headp_free=npart+1
  tailp_free=npartmax
  numbp_free=tailp_free-headp_free+1
  if(numbp_free>0)then
     prevp(headp_free)=0
  end if
  nextp(tailp_free)=0
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,&
       & MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbp_free_tot=numbp_free
#endif

  !--------------
  ! Periodic box
  !--------------
  do idim=1,ndim
     do ipart=1,npart
        if(xp(ipart,idim)/scale+skip_loc(idim)<0.0d0) &
             & xp(ipart,idim)=xp(ipart,idim)+(xbound(idim)-skip_loc(idim))*scale
        if(xp(ipart,idim)/scale+skip_loc(idim)>=xbound(idim)) &
             & xp(ipart,idim)=xp(ipart,idim)-(xbound(idim)-skip_loc(idim))*scale
     end do
     if(sink)then
        do ipart=1,nsink
           if(xsink(ipart,idim)/scale+skip_loc(idim)<0.0d0) &
                & xsink(ipart,idim)=xsink(ipart,idim)+(xbound(idim)-skip_loc(idim))*scale
           if(xsink(ipart,idim)/scale+skip_loc(idim)>=xbound(idim)) &
                & xsink(ipart,idim)=xsink(ipart,idim)-(xbound(idim)-skip_loc(idim))*scale
        end do
     endif
  end do
 
  !----------------------------------
  ! Reset all linked lists at level 1
  !----------------------------------
  do i=1,active(1)%ngrid
     headp(active(1)%igrid(i))=0
     tailp(active(1)%igrid(i))=0
     numbp(active(1)%igrid(i))=0
  end do
  do icpu=1,ncpu
     do i=1,reception(icpu,1)%ngrid
        headp(reception(icpu,1)%igrid(i))=0
        tailp(reception(icpu,1)%igrid(i))=0
        numbp(reception(icpu,1)%igrid(i))=0
     end do
  end do

  !------------------------------------------------
  ! Build linked list at level 1 by vector sweeps
  !------------------------------------------------
  do ipart=1,npart,nvector
     npart1=min(nvector,npart-ipart+1)
     ! Gather particles
     do i=1,npart1
        ind_part(i)=ipart+i-1
     end do
     ! Compute coarse cell
#if NDIM>0
     do i=1,npart1
        ix(i)=int(xp(ind_part(i),1)/scale+skip_loc(1))
     end do
#endif
#if NDIM>1
     do i=1,npart1
        iy(i)=int(xp(ind_part(i),2)/scale+skip_loc(2))
     end do
#endif
#if NDIM>2
     do i=1,npart1
        iz(i)=int(xp(ind_part(i),3)/scale+skip_loc(3))
     end do
#endif
     ! Compute level 1 grid index
     error=.false.
     do i=1,npart1
        ind_grid(i)=son(1+ix(i)+nx*iy(i)+nxny*iz(i))
        if(ind_grid(i)==0)error=.true.
     end do
     if(error)then
        write(*,*)'Error in init_tree'
        write(*,*)'Particles appear in unrefined regions'
        call clean_stop
     end if
     ! Add particle to level 1 linked list
     call add_list(ind_part,ind_grid,ok,npart1)
  end do

  if (sink.and.nsink>0.and.nrestart==0)call create_cloud(1)

  ! Sort particles down to levelmin
  do ilevel=1,levelmin-1
     call make_tree_fine(ilevel)
     call kill_tree_fine(ilevel)
     ! Update boundary conditions for remaining particles
     call virtual_tree_fine(ilevel)
  end do

  !show particle tree at startup
  ! do ilevel=1,nlevelmax
  !    call count_parts(ilevel)
  ! end do
  ! if(myid==1)print*,'+++++++++++++++++++++++++++++++++++++++++doneinit'
  

end subroutine init_tree
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_tree_fine(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
  !-----------------------------------------------------------------------
  ! This subroutine checks if particles have moved from their parent grid
  ! to one of the 3**ndim neighboring sister grids. The particle is then 
  ! disconnected from the parent grid linked list, and connected to the
  ! corresponding sister grid linked list. If the sister grid does
  ! not exist, the particle is left to its original parent grid.
  ! Particles must not move to a distance greater than direct neighbors
  ! boundaries. Otherwise an error message is issued and the code stops.
  !-----------------------------------------------------------------------
  integer::idim,nx_loc
  real(dp)::dx,scale
  real(dp),dimension(1:3)::xbound
  real(dp),dimension(1:3)::skip_loc
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,icpu
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx=0.5D0**ilevel    
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle  <--- Very important !!!
              next_part=nextp(ipart)
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig
              ! Gather nvector particles
              if(ip==nvector)then
                 call check_tree(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
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
     if(ip>0)call check_tree(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do
  ! End loop over cpus

  ! Periodic boundaries
  if(sink)then
     do idim=1,ndim
        do ipart=1,nsink
           if(xsink(ipart,idim)/scale+skip_loc(idim)<0.0d0) &
                & xsink(ipart,idim)=xsink(ipart,idim)+(xbound(idim)-skip_loc(idim))*scale
           if(xsink(ipart,idim)/scale+skip_loc(idim)>=xbound(idim)) &
                & xsink(ipart,idim)=xsink(ipart,idim)-(xbound(idim)-skip_loc(idim))*scale
        end do
     end do
  endif

111 format('   Entering make_tree_fine for level ',I2)

end subroutine make_tree_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine check_tree(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by make_tree_fine.
  !-----------------------------------------------------------------------
  logical::error
  integer::i,j,idim,nx_loc
  real(dp)::dx,xxx,scale
  real(dp),dimension(1:3)::xbound
  ! Grid-based arrays
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_father
  ! Particle-based arrays
  integer,dimension(1:nvector),save::ind_son,igrid_son
  integer,dimension(1:nvector),save::list1,list2
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:3)::skip_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel    
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_father(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_father,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Compute particle position in 3-cube
  error=.false.
  ind_son(1:np)=1
  ok(1:np)=.false.
  do idim=1,ndim
     do j=1,np
        i=int((xp(ind_part(j),idim)/scale+skip_loc(idim)-x0(ind_grid_part(j),idim))/dx/2.0D0)
        if(i<0.or.i>2)error=.true.
        ind_son(j)=ind_son(j)+i*3**(idim-1)
        ! Check if particle has escaped from its parent grid
        ok(j)=ok(j).or.i.ne.1 
     end do
  end do
  if(error)then
     write(*,*)'Problem in check_tree at level ',ilevel
     write(*,*)'A particle has moved outside allowed boundaries'
     do idim=1,ndim
        do j=1,np
           i=int((xp(ind_part(j),idim)/scale+skip_loc(idim)-x0(ind_grid_part(j),idim))/dx/2.0D0)
           if(i<0.or.i>2)then
              write(*,*)xp(ind_part(j),idim),x0(ind_grid_part(j),idim)
           endif
     end do
  end do
     stop
  end if

  ! Compute neighboring grid index
  do j=1,np
     igrid_son(j)=son(nbors_father_cells(ind_grid_part(j),ind_son(j)))
  end do

  ! If escaped particle sits in unrefined cell, leave it to its parent grid.
  ! For ilevel=levelmin, this should never happen.
  do j=1,np
     if(igrid_son(j)==0)ok(j)=.false.
  end do

  ! Periodic box
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           xxx=xp(ind_part(j),idim)/scale+skip_loc(idim)-xg(igrid_son(j),idim)
           if(xxx> xbound(idim)/2.0)then
              xp(ind_part(j),idim)=xp(ind_part(j),idim)-(xbound(idim)-skip_loc(idim))*scale
           endif
           if(xxx<-xbound(idim)/2.0)then
              xp(ind_part(j),idim)=xp(ind_part(j),idim)+(xbound(idim)-skip_loc(idim))*scale
           endif
        endif
     enddo
  enddo

  ! Switch particles linked list
  do j=1,np
     if(ok(j))then
        list1(j)=ind_grid(ind_grid_part(j))
        list2(j)=igrid_son(j)
     end if
  end do
  call remove_list(ind_part,list1,ok,np)
  call add_list(ind_part,list2,ok,np)

end subroutine check_tree
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kill_tree_fine(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine sorts particle between ilevel grids and their 
  ! ilevel+1 children grids. Particles are disconnected from their parent 
  ! grid linked list and connected to their corresponding child grid linked 
  ! list. If the  child grid does not exist, the particle is left to its 
  ! original parent grid. 
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::i,ig,ip,npart1,icpu
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel+1)==0)return
  if(verbose)write(*,111)ilevel

  ! Reset all linked lists at level ilevel+1
  do i=1,active(ilevel+1)%ngrid
     headp(active(ilevel+1)%igrid(i))=0
     tailp(active(ilevel+1)%igrid(i))=0
     numbp(active(ilevel+1)%igrid(i))=0
  end do
  do icpu=1,ncpu
     do i=1,reception(icpu,ilevel+1)%ngrid
        headp(reception(icpu,ilevel+1)%igrid(i))=0
        tailp(reception(icpu,ilevel+1)%igrid(i))=0
        numbp(reception(icpu,ilevel+1)%igrid(i))=0
     end do
  end do

  ! Sort particles between ilevel and ilevel+1

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig   
              if(ip==nvector)then
                 call kill_tree(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
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
     if(ip>0)call kill_tree(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

111 format('   Entering kill_tree_fine for level ',I2)

end subroutine kill_tree_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kill_tree(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine kill_tree_fine.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc
  real(dp)::dx,xxx,scale
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  integer,dimension(1:nvector),save::list1,list2
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:3)::skip_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel   
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Compute lower left corner of grid
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-dx
     end do
  end do

  ! Select only particles within grid boundaries
  ok(1:np)=.true.
  do idim=1,ndim
     do j=1,np
        xxx=(xp(ind_part(j),idim)/scale+skip_loc(idim)-x0(ind_grid_part(j),idim))/dx
        ok(j)=ok(j) .and. (xxx >= 0.d0 .and. xxx < 2.0d0)
     end do
  end do
  
  ! Determines in which son particles sit
  ind_son(1:np)=0
  do idim=1,ndim
     do j=1,np
        i=int((xp(ind_part(j),idim)/scale+skip_loc(idim)-x0(ind_grid_part(j),idim))/dx)
        ind_son(j)=ind_son(j)+i*2**(idim-1)
     end do
  end do 
  do j=1,np
     ind_son(j)=ncoarse+ind_son(j)*ngridmax+ind_grid(ind_grid_part(j))
  end do

  ! Determine which son cell is refined
  igrid_son(1:np)=0
  do j=1,np
     if(ok(j))igrid_son(j)=son(ind_son(j))
  end do
  do j=1,np
     ok(j)=igrid_son(j)>0
  end do

  ! Compute particle linked list
  do j=1,np
     if(ok(j))then
        list1(j)=ind_grid(ind_grid_part(j))
        list2(j)=igrid_son(j)
     end if
  end do

  ! Remove particles from their original linked lists
  call remove_list(ind_part,list1,ok,np)
  ! Add particles to their new linked lists
  call add_list(ind_part,list2,ok,np)

end subroutine kill_tree
!################################################################
!################################################################
!################################################################
!################################################################
subroutine merge_tree_fine(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine disconnects all particles contained in children grids
  ! and connects them to their parent grid linked list.
  !---------------------------------------------------------------
  integer::igrid,iskip,icpu
  integer::i,ind,ncache,ngrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_grid_son
  logical,dimension(1:nvector),save::ok

  if(numbtot(1,ilevel)==0)return
  if(ilevel==nlevelmax)return
  if(verbose)write(*,111)ilevel

  ! Loop over cpus
  do icpu=1,ncpu
     if(icpu==myid)then
        ncache=active(ilevel)%ngrid
     else
        ncache=reception(icpu,ilevel)%ngrid
     end if
     ! Loop over grids by vector sweeps
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        if(icpu==myid)then
           do i=1,ngrid
              ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
           end do
        else
           do i=1,ngrid
              ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
           end do
        end if
        ! Loop over children grids
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ngrid
              ind_grid_son(i)=son(ind_cell(i))
           end do
           do i=1,ngrid
              ok(i)=ind_grid_son(i)>0
           end do
           do i=1,ngrid
           if(ok(i))then
           if(numbp(ind_grid_son(i))>0)then
              if(numbp(ind_grid(i))>0)then
                 ! Connect son linked list at the tail of father linked list
                 nextp(tailp(ind_grid(i)))=headp(ind_grid_son(i))
                 prevp(headp(ind_grid_son(i)))=tailp(ind_grid(i))
                 numbp(ind_grid(i))=numbp(ind_grid(i))+numbp(ind_grid_son(i))
                 tailp(ind_grid(i))=tailp(ind_grid_son(i))
              else
                 ! Initialize father linked list
                 headp(ind_grid(i))=headp(ind_grid_son(i))
                 tailp(ind_grid(i))=tailp(ind_grid_son(i))
                 numbp(ind_grid(i))=numbp(ind_grid_son(i))
              end if

           end if
           end if
           end do
        end do
        ! End loop over children
     end do
     ! End loop over grids
  end do
  ! End loop over cpus

111 format('   Entering merge_tree_fine for level ',I2)

end subroutine merge_tree_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine virtual_tree_fine(ilevel)
  use pm_commons
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-----------------------------------------------------------------------
  ! This subroutine move particles across processors boundaries.
  !-----------------------------------------------------------------------
  integer::igrid,ipart,jpart,ncache_tot,next_part
  integer::ip,ipcom,npart1,icpu,ncache
  integer::info,buf_count,tag=101,tagf=102,tagu=102
  integer::countsend,countrecv
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE,2*ncpu)::statuses
  integer,dimension(2*ncpu)::reqsend,reqrecv
  integer,dimension(ncpu)::sendbuf,recvbuf
#endif
  integer,dimension(1:nvector),save::ind_part,ind_list,ind_com
  logical::ok_free,ok_all
  integer::particle_data_width

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

#ifdef WITHOUTMPI
  return
#endif

#ifndef WITHOUTMPI
  ! Count particle sitting in virtual boundaries
  do icpu=1,ncpu
     reception(icpu,ilevel)%npart=0
     do igrid=1,reception(icpu,ilevel)%ngrid
        reception(icpu,ilevel)%npart=reception(icpu,ilevel)%npart+&
             & numbp(reception(icpu,ilevel)%igrid(igrid))
     end do
     sendbuf(icpu)=reception(icpu,ilevel)%npart
  end do

  ! Calculate how many particle properties are being transferred
  particle_data_width = twondim+1
  if(star.or.sink) then
     if(metal) then
        particle_data_width=twondim+3
     else
        particle_data_width=twondim+2
     endif
  endif

#ifdef OUTPUT_PARTICLE_POTENTIAL
  particle_data_width=particle_data_width+1
#endif
  
  ! Allocate communication buffer in emission
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%npart
     if(ncache>0)then
        ! Allocate reception buffer
        allocate(reception(icpu,ilevel)%fp(1:ncache,1:3))
        allocate(reception(icpu,ilevel)%up(1:ncache,1:particle_data_width))
     end if
  end do

  ! Gather particle in communication buffer
  do icpu=1,ncpu
     if(reception(icpu,ilevel)%npart>0)then
     ! Gather particles by vector sweeps
     ipcom=0
     ip=0
     do igrid=1,reception(icpu,ilevel)%ngrid
        npart1=numbp(reception(icpu,ilevel)%igrid(igrid))
        ipart =headp(reception(icpu,ilevel)%igrid(igrid))
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle  <--- Very important !!!
           next_part=nextp(ipart)
           ip=ip+1
           ipcom=ipcom+1
           ind_com (ip)=ipcom
           ind_part(ip)=ipart
           ind_list(ip)=reception(icpu,ilevel)%igrid(igrid)
           reception(icpu,ilevel)%fp(ipcom,1)=igrid
           if(ip==nvector)then
              call fill_comm(ind_part,ind_com,ind_list,ip,ilevel,icpu)
              ip=0
           end if
           ipart=next_part  ! Go to next particle
        end do
     end do
     if(ip>0)call fill_comm(ind_part,ind_com,ind_list,ip,ilevel,icpu)
     end if
  end do

  ! Communicate virtual particle number to parent cpu
  call MPI_ALLTOALL(sendbuf,1,MPI_INTEGER,recvbuf,1,MPI_INTEGER,MPI_COMM_WORLD,info)

  ! Allocate communication buffer in reception
  do icpu=1,ncpu
     emission(icpu,ilevel)%npart=recvbuf(icpu)
     ncache=emission(icpu,ilevel)%npart
     if(ncache>0)then
        ! Allocate reception buffer
        allocate(emission(icpu,ilevel)%fp(1:ncache,1:3))
        allocate(emission(icpu,ilevel)%up(1:ncache,1:particle_data_width))
     end if
  end do

  ! Receive particles
  countrecv=0
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%npart
     if(ncache>0)then
        buf_count=ncache*3
        countrecv=countrecv+1
        call MPI_IRECV(emission(icpu,ilevel)%fp,buf_count, &
             & MPI_INTEGER,icpu-1,&
             & tagf,MPI_COMM_WORLD,reqrecv(countrecv),info)
        buf_count=ncache*particle_data_width
        countrecv=countrecv+1
        call MPI_IRECV(emission(icpu,ilevel)%up,buf_count, &
             & MPI_DOUBLE_PRECISION,icpu-1,&
             & tagu,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Send particles
  countsend=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%npart
     if(ncache>0)then
        buf_count=ncache*3
        countsend=countsend+1
        call MPI_ISEND(reception(icpu,ilevel)%fp,buf_count, &
             & MPI_INTEGER,icpu-1,&
             & tagf,MPI_COMM_WORLD,reqsend(countsend),info)
        buf_count=ncache*particle_data_width
        countsend=countsend+1
        call MPI_ISEND(reception(icpu,ilevel)%up,buf_count, &
             & MPI_DOUBLE_PRECISION,icpu-1,&
             & tagu,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Compute total number of newly created particles
  ncache_tot=0
  do icpu=1,ncpu
     ncache_tot=ncache_tot+emission(icpu,ilevel)%npart
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,&
       & MPI_COMM_WORLD,info)
  ok_free=(numbp_free-ncache_tot)>=0
  if(.not. ok_free)then
     write(*,*)'No more free memory for particles'
     write(*,*)'Increase npartmax'
     write(*,*)numbp_free,ncache_tot
     write(*,*)myid
     write(*,*)emission(1:ncpu,ilevel)%npart
     write(*,*)'============================'
     write(*,*)reception(1:ncpu,ilevel)%npart
     call MPI_ABORT(MPI_COMM_WORLD,1,info)
  end if

  ! Scatter new particles from communication buffer
  do icpu=1,ncpu
     ! Loop over particles by vector sweeps
     ncache=emission(icpu,ilevel)%npart
     do ipart=1,ncache,nvector
        npart1=min(nvector,ncache-ipart+1)
        do ip=1,npart1
           ind_com(ip)=ipart+ip-1
        end do
        call empty_comm(ind_com,npart1,ilevel,icpu)
     end do
  end do

  ! Deallocate temporary communication buffers
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%npart
     if(ncache>0)then
        deallocate(emission(icpu,ilevel)%fp)
        deallocate(emission(icpu,ilevel)%up)
     end if
     ncache=reception(icpu,ilevel)%npart
     if(ncache>0)then
        deallocate(reception(icpu,ilevel)%fp)
        deallocate(reception(icpu,ilevel)%up)
     end if
  end do
#endif

111 format('   Entering virtual_tree_fine for level ',I2)
end subroutine virtual_tree_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine fill_comm(ind_part,ind_com,ind_list,np,ilevel,icpu)
  use pm_commons
  use amr_commons
  implicit none
  integer::np,ilevel,icpu
  integer,dimension(1:nvector)::ind_part,ind_com,ind_list
  integer::current_property
  integer::i,idim
  logical,dimension(1:nvector),save::ok=.true.

  ! Gather particle level and identity
  do i=1,np
     reception(icpu,ilevel)%fp(ind_com(i),2)=levelp(ind_part(i))
     reception(icpu,ilevel)%fp(ind_com(i),3)=idp   (ind_part(i))
  end do
  
  ! Gather particle position and velocity
  do idim=1,ndim
     do i=1,np
        reception(icpu,ilevel)%up(ind_com(i),idim     )=xp(ind_part(i),idim)
        reception(icpu,ilevel)%up(ind_com(i),idim+ndim)=vp(ind_part(i),idim)
     end do
  end do
  
  current_property = twondim+1
  ! Gather particle mass
  do i=1,np
     reception(icpu,ilevel)%up(ind_com(i),current_property)=mp(ind_part(i))
  end do
  current_property = current_property+1

#ifdef OUTPUT_PARTICLE_POTENTIAL
  ! Gather particle potential
  do i=1,np
     reception(icpu,ilevel)%up(ind_com(i),current_property)=ptcl_phi(ind_part(i))
  end do
  current_property = current_property+1
#endif
  
  ! Gather particle birth epoch
  if(star.or.sink)then
     do i=1,np
        reception(icpu,ilevel)%up(ind_com(i),current_property)=tp(ind_part(i))
     end do
     if(metal)then
        do i=1,np
           reception(icpu,ilevel)%up(ind_com(i),current_property+1)=zp(ind_part(i))
        end do
     end if
  end if

  ! following line is not strictly necessary, but in case one adds extra data later
  current_property = current_property + 2 
  
  ! Remove particles from parent linked list
  call remove_list(ind_part,ind_list,ok,np)
  call add_free(ind_part,np)
  
end subroutine fill_comm
!################################################################
!################################################################
!################################################################
!################################################################
subroutine empty_comm(ind_com,np,ilevel,icpu)
  use pm_commons
  use amr_commons
  implicit none
  integer::np,icpu,ilevel
  integer,dimension(1:nvector)::ind_com
  
  integer::i,idim,igrid
  integer,dimension(1:nvector),save::ind_list,ind_part
  logical,dimension(1:nvector),save::ok=.true.
  integer::current_property

  ! Compute parent grid index
  do i=1,np
     igrid=emission(icpu,ilevel)%fp(ind_com(i),1)
     ind_list(i)=emission(icpu,ilevel)%igrid(igrid)
  end do

  ! Add particle to parent linked list
  call remove_free(ind_part,np)
  call add_list(ind_part,ind_list,ok,np)

  ! Scatter particle level and identity
  do i=1,np
     levelp(ind_part(i))=emission(icpu,ilevel)%fp(ind_com(i),2)
     idp   (ind_part(i))=emission(icpu,ilevel)%fp(ind_com(i),3)
  end do

  ! Scatter particle position and velocity
  do idim=1,ndim
  do i=1,np
     xp(ind_part(i),idim)=emission(icpu,ilevel)%up(ind_com(i),idim     )
     vp(ind_part(i),idim)=emission(icpu,ilevel)%up(ind_com(i),idim+ndim)
  end do
  end do

  current_property = twondim+1

  ! Scatter particle mass
  do i=1,np
     mp(ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),current_property)
  end do
  current_property = current_property+1

#ifdef OUTPUT_PARTICLE_POTENTIAL
  ! Scatter particle phi
  do i=1,np
     ptcl_phi(ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),current_property)
  end do
  current_property = current_property+1
#endif

  ! Scatter particle birth eopch
  if(star.or.sink)then
     do i=1,np
        tp(ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),current_property)
     end do
     if(metal)then
        do i=1,np
           zp(ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),current_property+1)
        end do
     end if
  end if

  ! As with the gather routine, we leave this in case extra properties are
  ! added later:
  current_property = current_property+2

end subroutine empty_comm
!################################################################
!################################################################
!################################################################
!################################################################
