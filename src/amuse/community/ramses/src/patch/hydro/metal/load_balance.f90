!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine load_balance
  use amr_commons
  use pm_commons
  use hydro_commons, ONLY: nvar, uold
  use poisson_commons, ONLY: f
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------
  ! This routine performs parallel load balancing.
  !------------------------------------------------
  integer::igrid,ncache,ilevel,i,ind,jlevel,info
  integer::idim,ivar,icpu,jcpu,kcpu
  integer::nxny,ix,iy,iz,iskip

  if(ncpu==1)return

#ifndef WITHOUTMPI
  if(myid==1)write(*,*)'Load balancing AMR grid...'
  
  ! Put all particle in main tree trunk
  if(pic)then
     call make_tree_fine(levelmin)
     do ilevel=levelmin-1,1,-1
        call merge_tree_fine(ilevel)
     end do
  endif

  balance=.true.

  if(verbose)then
     write(*,*)'Input mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if

  !-------------------------------------------
  ! Compute new cpu map using choosen ordering
  !-------------------------------------------
  call cmp_new_cpu_map

  !------------------------------------------------------
  ! Expand boundaries to account for new mesh partition
  !------------------------------------------------------
  call flag_coarse
  call refine_coarse
  call build_comm(1)
  call make_virtual_fine_int(cpu_map (1),1)
  call make_virtual_fine_int(cpu_map2(1),1)
  do i=1,nlevelmax-1
     call flag_fine(i,2)
     call refine_fine(i)
     call build_comm(i+1)
     call make_virtual_fine_int(cpu_map (1),i+1)
     call make_virtual_fine_int(cpu_map2(1),i+1)
  end do

  !--------------------------------------
  ! Update physical boundary conditions
  !--------------------------------------
  do ilevel=nlevelmax,1,-1
     if(hydro)then
        do ivar=1,nvar
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
        end do
        if(poisson)then
           do idim=1,ndim
              call make_virtual_fine_dp(f(1,idim),ilevel)
           end do
        end if
     end if
     if(simple_boundary)then
        if(hydro)call make_boundary_hydro(ilevel)
     end if
  end do

  !--------------------------------------
  ! Rearrange oct's between cpu's
  !--------------------------------------
  do ilevel=1,nlevelmax
     do icpu=1,ncpu
        if(icpu==myid)then
           ncache=active(ilevel)%ngrid
        else
           ncache=reception(icpu,ilevel)%ngrid
        end if
        ! Disconnect from old linked list
        do i=1,ncache
           if(icpu==myid)then
              igrid=active(ilevel)%igrid(i)
           else
              igrid=reception(icpu,ilevel)%igrid(i)
           end if
           kcpu=cpu_map (father(igrid))
           jcpu=cpu_map2(father(igrid))
           if(kcpu.ne.jcpu)then
              if(prev(igrid).ne.0) then
                 if(next(igrid).ne.0)then
                    next(prev(igrid))=next(igrid)
                    prev(next(igrid))=prev(igrid)
                 else
                    next(prev(igrid))=0
                    taill(kcpu,ilevel)=prev(igrid)
                 end if
              else
                 if(next(igrid).ne.0)then
                    prev(next(igrid))=0
                    headl(kcpu,ilevel)=next(igrid)
                 else
                    headl(kcpu,ilevel)=0
                    taill(kcpu,ilevel)=0
                 end if
              end if
              numbl(kcpu,ilevel)=numbl(kcpu,ilevel)-1 
           end if
        end do        
        ! Connect to new linked list
        do i=1,ncache
           if(icpu==myid)then
              igrid=active(ilevel)%igrid(i)
           else
              igrid=reception(icpu,ilevel)%igrid(i)
           end if
           kcpu=cpu_map (father(igrid))
           jcpu=cpu_map2(father(igrid))
           if(kcpu.ne.jcpu)then
              if(numbl(jcpu,ilevel)>0)then
                 next(igrid)=0
                 prev(igrid)=taill(jcpu,ilevel)
                 next(taill(jcpu,ilevel))=igrid
                 taill(jcpu,ilevel)=igrid
                 numbl(jcpu,ilevel)=numbl(jcpu,ilevel)+1
              else
                 next(igrid)=0
                 prev(igrid)=0
                 headl(jcpu,ilevel)=igrid
                 taill(jcpu,ilevel)=igrid
                 numbl(jcpu,ilevel)=1
              end if
           end if
        end do
     end do
  end do
  !--------------------------------------
  ! Compute new grid number statistics
  !--------------------------------------
  do ilevel=1,nlevelmax
  call MPI_ALLREDUCE(numbl(myid,ilevel),numbtot(1,ilevel),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(numbl(myid,ilevel),numbtot(2,ilevel),1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(numbl(myid,ilevel),numbtot(3,ilevel),1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(used_mem          ,used_mem_tot     ,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  numbtot(4,ilevel)=numbtot(1,ilevel)/ncpu
  end do

  !--------------------------------------
  ! Set old cpu map to new cpu map
  !--------------------------------------
  bound_key=bound_key2
  nxny=nx*ny
  do iz=kcoarse_min,kcoarse_max
  do iy=jcoarse_min,jcoarse_max
  do ix=icoarse_min,icoarse_max
     ind=1+ix+iy*nx+iz*nxny
     cpu_map(ind)=cpu_map2(ind)
  end do
  end do
  end do
  do ilevel=1,nlevelmax
     ! Build new communicators
     call build_comm(ilevel)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           active(ilevel)%igrid(i)=active(ilevel)%igrid(i)+iskip
        end do
        do i=1,active(ilevel)%ngrid
           cpu_map(active(ilevel)%igrid(i))=cpu_map2(active(ilevel)%igrid(i))
        end do
        do i=1,active(ilevel)%ngrid
           active(ilevel)%igrid(i)=active(ilevel)%igrid(i)-iskip
        end do
     end do
     call make_virtual_fine_int(cpu_map(1),ilevel)
  end do

  if(pic)then
     ! Sort particles down to nlevelmax
     do ilevel=1,nlevelmax-1
        call kill_tree_fine(ilevel)
        call virtual_tree_fine(ilevel)
     end do
     call virtual_tree_fine(nlevelmax)
     do ilevel=nlevelmax-1,levelmin,-1
        call merge_tree_fine(ilevel)
     end do
  end if

  !--------------------------------------------
  ! Shrink boundaries around new mesh partition
  !--------------------------------------------
  shrink=.true.
  do i=nlevelmax-1,1,-1
     call flag_fine(i,2)
     call refine_fine(i)
     call build_comm(i+1)
  end do  
  call flag_coarse
  call refine_coarse
  call build_comm(1)
  shrink=.false.

  balance=.false.

  if(verbose)then
     write(*,*)'Output mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if
#endif

999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine load_balance
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_new_cpu_map
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !---------------------------------------------------
  ! This routine computes the new cpu map using 
  ! the choosen ordering to balance load across cpus.
  !---------------------------------------------------
  integer::igrid,ncell,ncell_loc,ncache,ngrid
  integer::ncode,bit_length,ilevel,i,ind,idim
  integer::nx_loc,ny_loc,nz_loc,nfar
  integer::info,icpu,jcpu
  integer::nxny,ix,iy,iz,iskip
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::dx,scale,weight
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  integer(kind=8)::incost_tot,ind_long
  integer(kind=8),dimension(0:ncpu)::incost_new,incost_old
  integer,dimension(1:ncpu)::cost_loc,cost_old,cost_new
  real(kind=8),dimension(0:ncpu)::bound_key_loc
  real(kind=8),dimension(1:nvector),save::order_min,order_max

  ! Local constants
  nxny=nx*ny
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  if(cosmo)scale=scale/boxlen

  !----------------------------------------
  ! Compute cell ordering and cost
  ! for leaf cells with cpu map = myid.
  ! Store cost in flag1 and MAXIMUM  
  ! ordering key in hilbert_key of kind=8.
  !----------------------------------------
  ncell=0
  ncell_loc=1
  dx=1.0*scale
  do iz=0,nz-1
  do iy=0,ny-1
  do ix=0,nx-1
     ind=1+ix+iy*nx+iz*nxny
     if(cpu_map(ind)==myid.and.son(ind)==0)then
        xx(1,1)=(dble(ix)+0.5d0-dble(icoarse_min))*scale
#if NDIM>1
        xx(1,2)=(dble(iy)+0.5d0-dble(jcoarse_min))*scale
#endif
#if NDIM>2
        xx(1,3)=(dble(iz)+0.5d0-dble(kcoarse_min))*scale
#endif
        call cmp_minmaxorder(xx,order_min,order_max,dx,ncell_loc)
        ncell=ncell+1
        flag1(ncell)=1
        hilbert_key(ncell)=order_max(1)
     end if
  end do
  end do
  end do
  ! Loop over levels
  do ilevel=1,nlevelmax
     icost=1
     if(ilevel>levelmin)icost=2^(n_subcycle(ilevel-levelmin))
     ! Cell size and cell center offset
     dx=0.5d0**ilevel
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5d0)*dx-dble(icoarse_min)
#if NDIM>1
        xc(ind,2)=(dble(iy)-0.5d0)*dx-dble(jcoarse_min)
#endif
#if NDIM>2
        xc(ind,3)=(dble(iz)-0.5d0)*dx-dble(kcoarse_min)
#endif
     end do
     ! Loop over cpus
     do icpu=1,ncpu
        if(icpu==myid)then
           ncache=active(ilevel)%ngrid
        else
           ncache=reception(icpu,ilevel)%ngrid
        end if
        ! Loop over grids by vector sweeps
        do igrid=1,ncache,nvector
           ! Gather nvector grids
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
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ngrid
                 ind_cell(i)=ind_grid(i)+iskip
              end do
              do idim=1,ndim
              ncell_loc=0
              do i=1,ngrid
              if(cpu_map(ind_cell(i))==myid.and.son(ind_cell(i))==0)then
                 ncell_loc=ncell_loc+1
                 xx(ncell_loc,idim)=(xg(ind_grid(i),idim)+xc(ind,idim))*scale
              end if
              end do
              end do
              if(ncell_loc>0)call cmp_minmaxorder(xx,order_min,order_max,dx*scale,ncell_loc)
              ncell_loc=0
              do i=1,ngrid
              if(cpu_map(ind_cell(i))==myid.and.son(ind_cell(i))==0)then
                 ncell    =ncell    +1
                 ncell_loc=ncell_loc+1
                 flag1(ncell)=icost+int(dble(numbp(father(ind_grid(i))))/8.d0)
                 hilbert_key(ncell)=order_max(ncell_loc)
              end if
              end do
           end do
           ! End loop over cells
        end do
        ! End loop over grids
     end do
     ! End loop over cpus
  end do
  ! End loop over levels

  !------------------------------------------------
  ! Sort ordering key and store new index in flag2
  !------------------------------------------------
  if(ncell>0)call quick_sort(hilbert_key(1),flag2(1),ncell)

  !-----------------------------
  ! Balance cost across cpus
  !-----------------------------
  cost_loc = 0 ! Compute local and global cost
  cost_loc(myid) = ncell
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(cost_loc,cost_old,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
  incost_tot = int(0,kind=8)
  incost_old(0) = int(0,kind=8)*
  do icpu = 1,ncpu
     incost_tot = incost_tot + int(cost_old(icpu),kind=8)
     incost_old(icpu) = incost_tot
  end do
  incost_new(0) = int(0,kind=8)
  do icpu = 1,ncpu
     cost_new(icpu) = incost_tot/int(ncpu,kind=8) ! Exact load balancing
     if (icpu <= mod(incost_tot,int(ncpu,kind=8))) then
        cost_new(icpu) = cost_new(icpu)+1
     endif
     incost_new(icpu) = incost_new(icpu-1) + int(cost_new(icpu),kind=8)
  end do

  !-----------------------------
  ! Compute new cpu boundaries
  !-----------------------------
  bound_key_loc=0.0d0; bound_key2=0.0d0
  if(ncell>0)then
     ! First cpu on local domain
     icpu=0
     do while(incost_new(icpu)<incost_old(myid-1))
        icpu=icpu+1
     end do
     ! Last cpu on local domain
     jcpu=ncpu
     do while(incost_new(jcpu)>incost_old(myid  ))
        jcpu=jcpu-1
     end do
     ! Compute Hilbert key at boundaries
     do i=icpu,jcpu
        ind_long=incost_new(i)-incost_old(myid-1)
        if(ind_long>0.and.ind_long<=ncell)then
           bound_key_loc(i)=hilbert_key(ind_long)
        endif
     end do
  end if
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(bound_key_loc,bound_key2,ncpu+1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
  bound_key2(0)   =order_all_min
  bound_key2(ncpu)=order_all_max

  if(verbose)then
     write(*,*)bound_key
     write(*,*)bound_key2
  endif

  !----------------------------------------
  ! Compute new cpu map
  !----------------------------------------
  cpu_map2=0
  ncell_loc=1
  do iz=0,nz-1
  do iy=0,ny-1
  do ix=0,nx-1
     ind=1+ix+iy*nx+iz*nxny
     xx(1,1)=(dble(ix)+0.5d0-dble(icoarse_min))*scale
#if NDIM>1
     xx(1,2)=(dble(iy)+0.5d0-dble(jcoarse_min))*scale
#endif
#if NDIM>2
     xx(1,3)=(dble(iz)+0.5d0-dble(kcoarse_min))*scale
#endif
     call cmp_ordering(xx,order_max,ncell_loc)
     cpu_map2(ind)=ncpu ! default value
     do icpu=1,ncpu
        if( order_max(1).ge.bound_key2(icpu-1).and. &
          & order_max(1).lt.bound_key2(icpu  ))then
           cpu_map2(ind)=icpu
        endif
     end do
  end do
  end do
  end do
  ! Loop over levels
  do ilevel=1,nlevelmax
     ! Cell size and cell center offset
     dx=0.5d0**ilevel
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5d0)*dx-dble(icoarse_min)
#if NDIM>1
        xc(ind,2)=(dble(iy)-0.5d0)*dx-dble(jcoarse_min)
#endif
#if NDIM>2
        xc(ind,3)=(dble(iz)-0.5d0)*dx-dble(kcoarse_min)
#endif
     end do     
     ncache=active(ilevel)%ngrid
     ! Loop over grids by vector sweeps
     do igrid=1,ncache,nvector
        ! Gather nvector grids
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
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=(xg(ind_grid(i),idim)+xc(ind,idim))*scale
              end do
           end do
           if(ngrid>0)call cmp_ordering(xx,order_max,ngrid)
           do i=1,ngrid
              cpu_map2(ind_cell(i))=ncpu ! default value
              do icpu=1,ncpu
                 if( order_max(i).ge.bound_key2(icpu-1).and. &
                   & order_max(i).lt.bound_key2(icpu  ))then
                    cpu_map2(ind_cell(i))=icpu
                 endif
              end do
           end do
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

  ! Update virtual boundaries for new cpu map
  call make_virtual_coarse_int(cpu_map2(1))
  do ilevel=1,nlevelmax
     call make_virtual_fine_int(cpu_map2(1),ilevel)
  end do 

end subroutine cmp_new_cpu_map
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_cpumap(x,c,nn)
  use amr_parameters
  use amr_commons
  implicit none
  integer ::nn
  integer ,dimension(1:nvector)::c
  real(dp),dimension(1:nvector,1:ndim)::x

  integer::i,icpu
  real(kind=8),dimension(1:nvector),save::order

  call cmp_ordering(x,order,nn)
  do i=1,nn
     c(i)=ncpu ! default value
     do icpu=1,ncpu
        if(    order(i).ge.bound_key(icpu-1).and. &
             & order(i).lt.bound_key(icpu  ))then
           c(i)=icpu
        endif
     end do
  end do

end subroutine cmp_cpumap
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_ordering(x,order,nn)
  use amr_parameters
  use amr_commons
  implicit none
  integer ::nn
  real(dp),dimension(1:nvector,1:ndim)::x
  real(kind=8),dimension(1:nvector)::order
  !--------------------------------------------------------
  ! This routine computes the index key of the input cell
  ! according to its position in space and for the chosen
  ! ordering. Position x are in user units.
  !-----------------------------------------------------
  integer,dimension(1:nvector),save::ix,iy,iz
  integer::i,ncode,bit_length,nx_loc
  real(kind=8)::scale,bscale

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  if(cosmo)scale=scale/boxlen

  if(ordering=='planar')then
     ! Planar domain decomposition
     do i=1,nn
        order(i)=x(i,3)
     end do
  end if

  if(ordering=='angular')then
#if NDIM>1
     ! Angular domain decomposition
     do i=1,nn
        order(i)=atan(x(i,2)/(x(i,1)+1d-10))
     end do
#endif
  end if

  if(ordering=='hilbert')then
     ! Hilbert curve domain decomposition
     bscale=2**(nlevelmax+1)
     ncode=nx_loc*int(bscale)
     bit_length=int(log(dble(ncode))/log(2.))+1
     bscale=bscale/scale
     
     do i=1,nn
        ix(i)=int(x(i,1)*bscale)
#if NDIM>1           
        iy(i)=int(x(i,2)*bscale)
#endif
#if NDIM>2
        iz(i)=int(x(i,3)*bscale)
#endif
     end do

     if(ndim==1)then
        call hilbert1d(ix,order,nn)
     else if(ndim==2)then
        call hilbert2d(ix,iy,order,bit_length,nn)
     else if (ndim==3)then
        call hilbert3d(ix,iy,iz,order,bit_length,nn)
     end if

  end if

end subroutine cmp_ordering
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_minmaxorder(x,order_min,order_max,dx,nn)
  use amr_parameters
  use amr_commons
  implicit none
  integer ::nn
  real(dp)::dx
  real(dp),dimension(1:nvector,1:ndim)::x
  real(kind=8),dimension(1:nvector)::order_min,order_max
  !-----------------------------------------------------
  ! This routine computes the minimum and maximum index
  ! key contained in the input cell and for the chosen 
  ! ordering.
  !-----------------------------------------------------
  integer,dimension(1:nvector),save::ix,iy,iz
  integer::i,ncode,bit_length,nxny,nx_loc

  real(dp)::theta1,theta2,theta3,theta4,dxx,dxmin  
  real(kind=8)::dkey,scale,bscaleloc,bscale

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  if(cosmo)scale=scale/boxlen

  dxmin=scale/dble(2**nlevelmax)

  if(ordering=='planar')then
     ! Planar domain decomposition
     dxx=0.5d0*dx
     do i=1,nn
        order_min(i)=x(i,3)-dxx
        order_max(i)=x(i,3)+dxx
     end do
  end if

  if(ordering=='angular')then
     ! Angular domain decomposition
     dxx=0.5d0*dx
#if NDIM>1
     do i=1,nn
        theta1=atan((x(i,2)-dxx)/(x(i,1)-dxx+1d-10))
        theta2=atan((x(i,2)-dxx)/(x(i,1)+dxx+1d-10))
        theta3=atan((x(i,2)+dxx)/(x(i,1)+dxx+1d-10))
        theta4=atan((x(i,2)+dxx)/(x(i,1)-dxx+1d-10))
        order_min(i)=min(theta1,theta2,theta3,theta4)
        order_max(i)=max(theta1,theta2,theta3,theta4)
     end do
#endif
  end if

  if(ordering=='hilbert')then
     ! Hilbert curve domain decomposition
     bscale=2**(nlevelmax+1)
     bscaleloc=2**nlevelmax*dxmin/dx
     ncode=nx_loc*int(bscaleloc)
     bit_length=int(log(dble(ncode))/log(2.0d0))
     bscaleloc=bscaleloc/scale
     bscale   =bscale   /scale
     

     do i=1,nn
        ix(i)=int(x(i,1)*bscaleloc)
#if NDIM>1           
        iy(i)=int(x(i,2)*bscaleloc)
#endif
#if NDIM>2
        iz(i)=int(x(i,3)*bscaleloc)
#endif
     end do

     if(ndim==1)then
        call hilbert1d(ix,order_min,nn)
     else if(ndim==2)then
        call hilbert2d(ix,iy,order_min,bit_length,nn)
     else if (ndim==3)then
        call hilbert3d(ix,iy,iz,order_min,bit_length,nn)
     end if

     dkey=(dble(bscale)/dble(bscaleloc))**ndim
     do i=1,nn
        order_max(i)=(order_min(i)+1.0d0)*dkey
        order_min(i)=(order_min(i))*dkey
     end do

  end if

end subroutine cmp_minmaxorder
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
