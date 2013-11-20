!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine load_balance
  use amr_commons
  use pm_commons
  use hydro_commons, ONLY: nvar, uold
#ifdef RT
  use rt_hydro_commons, ONLY: nrtvar, rtuold
#endif
  use poisson_commons, ONLY: phi, f
  use bisection
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
  integer,dimension(nlevelmax,3)::comm_buffin,comm_buffout

  if(ncpu==1)return

#ifndef WITHOUTMPI
  if(myid==1)write(*,*)'Load balancing AMR grid...'
  
  ! Put all particle in main tree trunk
  if(pic.and.(.not.init))then
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
  ! Compute new cpu map using chosen ordering
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
#ifdef SOLVERmhd
        do ivar=1,nvar+3
#else
        do ivar=1,nvar
#endif
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
        end do
#else
        end do
#endif
        if(simple_boundary)then
           call make_boundary_hydro(ilevel)
        end if
     end if
#ifdef RT
     if(rt)then
        do ivar=1,nrtvar
           call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
        end do
        if(simple_boundary)then
           call rt_make_boundary_hydro(ilevel)
        end if
     endif
#endif
     if(poisson)then
        call make_virtual_fine_dp(phi(1),ilevel)
        do idim=1,ndim
           call make_virtual_fine_dp(f(1,idim),ilevel)
        end do
     end if
  end do

  !--------------------------------------
  ! Rearrange octs between cpus
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
     comm_buffin(ilevel,1)=numbl(myid,ilevel)
     comm_buffin(ilevel,2)=numbl(myid,ilevel)
     comm_buffin(ilevel,3)=numbl(myid,ilevel)
  end do
  call MPI_ALLREDUCE(comm_buffin(1,1),comm_buffout(1,1),nlevelmax,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(comm_buffin(1,2),comm_buffout(1,2),nlevelmax,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(comm_buffin(1,3),comm_buffout(1,3),nlevelmax,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(used_mem        ,used_mem_tot     ,1        ,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  do ilevel=1,nlevelmax
     numbtot(1,ilevel)=comm_buffout(ilevel,1)
     numbtot(2,ilevel)=comm_buffout(ilevel,2)
     numbtot(3,ilevel)=comm_buffout(ilevel,3)
     numbtot(4,ilevel)=numbtot(1,ilevel)/ncpu
  end do

  !--------------------------------------
  ! Set old cpu map to new cpu map
  !--------------------------------------
  if(ordering/='bisection') then
     bound_key=bound_key2
  else
     bisec_cpubox_min=bisec_cpubox_min2
     bisec_cpubox_max=bisec_cpubox_max2
  end if

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
           cpu_map(active(ilevel)%igrid(i)+iskip)=cpu_map2(active(ilevel)%igrid(i)+iskip)
        end do
     end do
     call make_virtual_fine_int(cpu_map(1),ilevel)
  end do

  if(pic.and.(.not.init))then
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
  use pm_commons
  use bisection
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
  integer::info,icpu,jcpu,isub,idom,jdom
  integer::nxny,ix,iy,iz,iskip
  integer::ind_long
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::dx,scale,weight
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(kind=8)::incost_tot,local_cost,cell_cost
  real(kind=8),dimension(0:ndomain)::incost_new,incost_old
  integer(kind=8),dimension(1:overload)::npart_sub
  integer(kind=8)::wflag
  integer,dimension(1:overload)::ncell_sub
  real(kind=8),dimension(1:ndomain)::cost_loc,cost_old,cost_new
  real(qdp),dimension(0:ndomain)::bound_key_loc
  real(kind=8),dimension(0:ndomain)::bigdbl,bigtmp
  integer,dimension(1:nvector),save::dom
  real(qdp),dimension(1:nvector),save::order_min,order_max
  integer,dimension(1:MAXLEVEL),save::niter_cost

  real(dp),dimension(1:1,1:ndim),save :: xx_tmp
  integer,dimension(1:1),save :: c_tmp

  ! Local constants
  nxny=nx*ny
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  
  ! Compute time step related cost
  if(cost_weighting)then
     niter_cost(levelmin)=1
     if (nlevelmax - levelmin - 1 > 31) write(*,*) 'Warning load_balance: niter_cost may need to become a kind=8 integer'
     do ilevel=levelmin+1,nlevelmax
        niter_cost(ilevel)=nsubcycle(ilevel-1)*niter_cost(ilevel-1)
     end do
  else
     niter_cost(levelmin:nlevelmax)=1
  endif

  if(verbose) print *,"Entering cmp_new_cpu_map"

  if(ordering/='bisection') then      ! begin if not bisection

  !----------------------------------------
  ! Compute cell ordering and cost
  ! for leaf cells with cpu map = myid.
  ! Store cost in flag1 and MAXIMUM  
  ! ordering key in hilbert_key of kind=16
  !----------------------------------------
  ncell=0
  npart_sub=0
  ncell_sub=0
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
        call cmp_dommap(xx,dom,ncell_loc)
        ncell=ncell+1
        isub=(dom(1)-1)/ncpu+1
        ncell_sub(isub)=ncell_sub(isub)+1
        flag1(ncell)=0
        hilbert_key(ncell)=order_max(1)
     end if
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
              if(ncell_loc>0)then
                 call cmp_minmaxorder(xx,order_min,order_max,dx*scale,ncell_loc)
                 call cmp_dommap(xx,dom,ncell_loc)
              end if
              ncell_loc=0
              do i=1,ngrid
                 if(cpu_map(ind_cell(i))==myid.and.son(ind_cell(i))==0)then
                    ncell    =ncell    +1
                    ncell_loc=ncell_loc+1
                    isub=(dom(ncell_loc)-1)/ncpu+1
                    ncell_sub(isub)=ncell_sub(isub)+1
                    flag1(ncell)=8*10 ! Magic number
                    if(pic)then
                       flag1(ncell)=flag1(ncell)+numbp(ind_grid(i))
                    endif
                    wflag = flag1(ncell)*niter_cost(ilevel)
                    if (wflag > 2147483647) then 
                       write(*,*) ' wrong type for flag1 --> change to integer kind=8: ',wflag
                       stop
                    endif
                    flag1(ncell)=flag1(ncell)*niter_cost(ilevel)
                    npart_sub(isub)=npart_sub(isub)+flag1(ncell)
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
  if (ncell>0) call quick_sort(hilbert_key(1),flag2(1),ncell)

  !-----------------------------
  ! Balance cost across cpus
  !-----------------------------
  cost_loc = 0 ! Compute local and global cost
  do isub=1,overload
     cost_loc(myid+(isub-1)*ncpu) = dble(npart_sub(isub))
  end do
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(cost_loc,cost_old,ndomain,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
  incost_tot = 0D0
  incost_old(0) = 0D0
  do idom = 1,ndomain
     incost_tot = incost_tot + cost_old(idom)
     incost_old(idom) = incost_tot
  end do
  incost_new(0) = 0D0
  do idom = 1,ndomain
     cost_new(idom) = incost_tot/dble(ndomain) ! Exact load balancing
     incost_new(idom) = incost_new(idom-1) + cost_new(idom)
  end do

  !-----------------------------
  ! Compute new cpu boundaries
  !-----------------------------
  bound_key_loc=0.0d0; bound_key2=0.0d0
  ncell_loc=0
  do isub=1,overload
     if(ncell_sub(isub)>0)then
        ! First cpu on local domain
        idom=0
        do while(incost_new(idom)<incost_old(myid-1+(isub-1)*ncpu))
           idom=idom+1
           if (idom > ndomain) exit 
        end do
        ! Compute Hilbert key at boundaries
        i=idom
        local_cost=incost_old(myid-1+(isub-1)*ncpu)
        do ind_long=1,ncell_sub(isub)
           cell_cost=dble(flag1(flag2(ind_long+ncell_loc)))
           local_cost=local_cost+cell_cost
           if (i > ndomain) exit
           if(incost_new(i)<local_cost)then
              bound_key_loc(i)=hilbert_key(ind_long+ncell_loc)
              i=i+1
           endif
        end do
     end if
     ncell_loc=ncell_loc+ncell_sub(isub)
  end do
#ifndef WITHOUTMPI
#ifdef QUADHILBERT
  bigdbl= real(bound_key_loc,kind=8)
  bigtmp= 0.0d0
  call MPI_ALLREDUCE(bigdbl,bigtmp,ndomain+1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  ! if call to mpi_sum with mpi_type=mpi_real16 is supported by mpi_allreduce we can do: 
  !call MPI_ALLREDUCE(bound_key_loc,bound_key2,ndomain+1,MPI_REAL16,MPI_SUM,MPI_COMM_WORLD,info)
  bound_key2         = real(bigtmp,kind=qdp)
#else
  call MPI_ALLREDUCE(bound_key_loc,bound_key2,ndomain+1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#endif
  bound_key2(0)      =order_all_min
  bound_key2(ndomain)=order_all_max

  else     ! doing bisection
     ! update the bisection                                                                             
     call build_bisection(update=.true.)
  end if   ! end if not bisection

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
     cpu_map2(ind)=ncpu ! default value                                                               

     if(ordering/='bisection') then
        call cmp_ordering(xx,order_max,ncell_loc)
        cpu_map2(ind)=ncpu ! default value
        do idom=1,ndomain
           if( order_max(1).ge.bound_key2(idom-1).and. &
                & order_max(1).lt.bound_key2(idom))then
              cpu_map2(ind)=mod(idom-1,ncpu)+1
           endif
        end do
     else
        xx_tmp(1,:) = xx(1,:)
        call cmp_bisection_cpumap(xx_tmp,c_tmp,1)
        cpu_map2(ind) = c_tmp(1)
     end if
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
           if(ordering/='bisection') then                                                                     
              if(ngrid>0)call cmp_ordering(xx,order_max,ngrid)
              do i=1,ngrid
                 cpu_map2(ind_cell(i))=ncpu ! default value
                 do idom=1,ndomain
                    if( order_max(i).ge.bound_key2(idom-1).and. &
                         & order_max(i).lt.bound_key2(idom))then
                       cpu_map2(ind_cell(i))=mod(idom-1,ncpu)+1
                    endif
                 end do
              end do
           else
              do i=1,ngrid
                 ! compute cpu_map2 using bisection                                                      
                 xx_tmp(1,:) = xx(i,:)
                 call cmp_bisection_cpumap(xx_tmp,c_tmp,1)
                 cpu_map2(ind_cell(i)) = c_tmp(1)
              end do
           endif
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
  use bisection
  implicit none
  integer ::nn
  integer ,dimension(1:nvector)::c
  real(dp),dimension(1:nvector,1:ndim)::x

  integer::i,idom
  real(qdp),dimension(1:nvector),save::order

  if(ordering /= 'bisection') then
     call cmp_ordering(x,order,nn)
     do i=1,nn
        c(i)=ndomain ! default value
        do idom=1,ndomain
           if(    order(i).ge.bound_key(idom-1).and. &
                & order(i).lt.bound_key(idom  ))then
              c(i)=idom
           endif
        end do
     end do
     do i=1,nn
        c(i)=MOD(c(i)-1,ncpu)+1
!        c(i)=c(i)-((c(i)-1)/ncpu)*ncpu
     end do
  else
     call cmp_bisection_cpumap(x,c,nn)
  end if

end subroutine cmp_cpumap
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_dommap(x,c,nn)
  use amr_parameters
  use amr_commons
  use bisection
  implicit none
  integer ::nn
  integer ,dimension(1:nvector)::c
  real(dp),dimension(1:nvector,1:ndim)::x

  integer::i,idom
  real(qdp),dimension(1:nvector),save::order

  call cmp_ordering(x,order,nn)
  do i=1,nn
     c(i)=ndomain ! default value
     do idom=1,ndomain
        if(    order(i).ge.bound_key(idom-1).and. &
             & order(i).lt.bound_key(idom  ))then
           c(i)=idom
        endif
     end do
  end do
  
end subroutine cmp_dommap
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_ordering(x,order,nn)
  use amr_parameters
  use amr_commons
  implicit none
  integer ::nn
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(dp),dimension(1:nvector,1:ndim)::x
  real(qdp),dimension(1:nvector)::order
  !--------------------------------------------------------
  ! This routine computes the index key of the input cell
  ! according to its position in space and for the chosen
  ! ordering. Position x are in user units.
  !-----------------------------------------------------
  integer,dimension(1:nvector),save::ix,iy,iz
  integer::i,ncode,bit_length,nx_loc
  integer::temp,info
  real(kind=8)::scale,bscale,xx,yy,zz,xc,yc,zc

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)

  if(ordering=='planar')then
     ! Planar domain decomposition
     do i=1,nn
        order(i)=x(i,1)
     end do
  end if

#if NDIM>1
  if(ordering=='angular')then
     ! Angular domain decomposition
     xc=boxlen/2.
     yc=boxlen/2.
     zc=boxlen/2.
     do i=1,nn
        xx=x(i,1)-xc+1d-10
        yy=x(i,2)-yc
#if NDIM>2
        zz=x(i,3)
#endif
        if(xx>0.)then
           order(i)=atan(yy/xx)+acos(-1.)/2.
        else
           order(i)=atan(yy/xx)+acos(-1.)*3./2.
        endif
#if NDIM>2
        if(zz.gt.zc)order(i)=order(i)+2.*acos(-1.)
#endif
     end do
  end if
#endif

  if(ordering=='hilbert')then
     ! Hilbert curve domain decomposition
     bscale=2**(nlevelmax+1)
     ncode=nx_loc*int(bscale)
     bscale=bscale/scale
     
     temp=ncode
     do bit_length=1,32
        ncode=ncode/2
        if(ncode<=1) exit
     end do
     if(bit_length==32) then
        write(*,*)'Error in cmp_minmaxorder'
#ifndef WITHOUTMPI
        call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
        stop
#endif
     end if

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
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer ::nn
  integer ::temp,info
  real(dp)::dx
  real(dp),dimension(1:nvector,1:ndim)::x
  real(qdp),dimension(1:nvector)::order_min,order_max
  !-----------------------------------------------------
  ! This routine computes the minimum and maximum index
  ! key contained in the input cell and for the chosen 
  ! ordering.
  !-----------------------------------------------------
  integer,dimension(1:nvector),save::ix,iy,iz
  integer::i,ncode,bit_length,nxny,nx_loc

  real(dp)::theta1,theta2,theta3,theta4,dxx,dxmin  
  real(kind=8)::scale,bscaleloc,bscale,xx,yy,zz,xc,yc,zc
  real(qdp)::dkey,oneqdp=1.0

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dxmin=scale/dble(2**nlevelmax)

  if(ordering=='planar')then
     ! Planar domain decomposition
     dxx=0.5d0*dx
     do i=1,nn
        order_min(i)=x(i,1)-dxx
        order_max(i)=x(i,1)+dxx
     end do
  end if

#if NDIM>1
  if(ordering=='angular')then
     ! Angular domain decomposition
     dxx=0.5d0*dx
     xc=boxlen/2.
     yc=boxlen/2.
     zc=boxlen/2.
     do i=1,nn
        if(dx==boxlen)then
           order_min(i)=0.
           order_max(i)=4.*acos(-1.)
        else
           ! x- y-
           yy=x(i,2)-yc-dxx
           xx=x(i,1)-xc-dxx
           if(xx.ge.0.)then
              xx=xx+1d-10
              theta1=atan(yy/xx)+acos(-1.)/2.
           else
              xx=xx-1d-10
              theta1=atan(yy/xx)+acos(-1.)*3./2.
           endif
           ! x+ y-
           xx=x(i,1)-xc+dxx
           if(xx.gt.0.)then
              xx=xx+1d-10
              theta2=atan(yy/xx)+acos(-1.)/2.
           else
              xx=xx-1d-10
              theta2=atan(yy/xx)+acos(-1.)*3./2.
           endif
           
           ! x+ y+
           yy=x(i,2)-yc+dxx
           if(xx.gt.0.)then
              xx=xx+1d-10
              theta3=atan(yy/xx)+acos(-1.)/2.
           else
              xx=xx-1d-10
              theta3=atan(yy/xx)+acos(-1.)*3./2.
           endif
           ! x- y+
           xx=x(i,1)-xc-dxx
           if(xx.ge.0.)then
              xx=xx+1d-10
              theta4=atan(yy/xx)+acos(-1.)/2.
           else
              xx=xx-1d-10
              theta4=atan(yy/xx)+acos(-1.)*3./2.
           endif
           order_min(i)=min(theta1,theta2,theta3,theta4)
           order_max(i)=max(theta1,theta2,theta3,theta4)
#if NDIM>2
           zz=x(i,3)
           if(zz.gt.zc)then
              order_min(i)=order_min(i)+2.*acos(-1.)
              order_max(i)=order_max(i)+2.*acos(-1.)
           endif
#endif
        endif
     end do
  end if
#endif

  if(ordering=='hilbert')then
     ! Hilbert curve domain decomposition
     bscale=2**(nlevelmax+1)
     bscaleloc=2**nlevelmax*dxmin/dx
     ncode=nx_loc*int(bscaleloc)
     bscaleloc=bscaleloc/scale
     bscale   =bscale   /scale
     
     temp=ncode
     do bit_length=1,32
        ncode=ncode/2
        if(ncode<=1) exit
     end do
     if(bit_length==32) then
        write(*,*)'Error in cmp_minmaxorder'
#ifndef WITHOUTMPI
        call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
        stop
#endif
     end if

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

     dkey=(real(bscale,kind=qdp)/real(bscaleloc,kind=qdp))**ndim
     do i=1,nn
        order_max(i)=(order_min(i)+oneqdp)*dkey
        order_min(i)=(order_min(i))*dkey
     end do

  end if

end subroutine cmp_minmaxorder
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine defrag
  use amr_commons
  use pm_commons
  use poisson_commons
  use hydro_commons
#ifdef RT
  use rt_hydro_commons
#endif
  implicit none

  integer::ncache,ngrid2,igridmax,i,igrid,ibound,ilevel
  integer::iskip1,iskip2,igrid1,igrid2,ind1,ind2,icell1,icell2
  integer::ind,idim,ivar,istart

  if(verbose)write(*,*)'Defragmenting main memory...'

  ngrid2=0
  igridmax=0
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              cpu_map2(igrid)=ngrid2+i
              igridmax=max(igridmax,igrid)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  
  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              icell1=father(igrid)
              if(icell1>ncoarse)then
                 ind1=(icell1-ncoarse-1)/ngridmax+1
                 iskip1=ncoarse+(ind1-1)*ngridmax
                 igrid1=(icell1-iskip1)
                 igrid2=cpu_map2(igrid1)
                 iskip2=ncoarse+(ind1-1)*ngridmax
                 icell2=iskip2+igrid2
              else
                 icell2=icell1
              end if
              flag2(ngrid2+i)=icell2
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     father(igrid)=flag2(igrid)
  end do

  do ind=1,twondim
  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              icell1=nbor(igrid,ind)
              if(icell1>ncoarse)then
                 ind1=(icell1-ncoarse-1)/ngridmax+1
                 iskip1=ncoarse+(ind1-1)*ngridmax
                 igrid1=(icell1-iskip1)
                 igrid2=cpu_map2(igrid1)
                 iskip2=ncoarse+(ind1-1)*ngridmax
                 icell2=iskip2+igrid2
              else
                 icell2=icell1
              end if
              flag2(ngrid2+i)=icell2
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     nbor(igrid,ind)=flag2(igrid)
  end do
  end do

  do idim=1,ndim
  ngrid2=0
  do igrid=1,igridmax
     hilbert_key(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              hilbert_key(ngrid2+i)=real(xg(igrid,idim),kind=qdp)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     xg(igrid,idim)=real(hilbert_key(igrid),kind=8)
  end do
  end do

  if(pic)then

  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              flag2(ngrid2+i)=headp(igrid)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     headp(igrid)=flag2(igrid)
  end do

  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              flag2(ngrid2+i)=tailp(igrid)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     tailp(igrid)=flag2(igrid)
  end do

  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              flag2(ngrid2+i)=numbp(igrid)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     numbp(igrid)=flag2(igrid)
  end do

  endif

  do ind=1,twotondim
  iskip2=ncoarse+(ind-1)*ngridmax
  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              igrid1=son(iskip2+igrid)
              if(igrid1>0)then
                 igrid2=cpu_map2(igrid1)
              else
                 igrid2=0
              end if
              flag2(ngrid2+i)=igrid2
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     son(iskip2+igrid)=flag2(igrid)
  end do
  end do

  do ind=1,twotondim
  iskip2=ncoarse+(ind-1)*ngridmax
  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              flag2(ngrid2+i)=cpu_map(iskip2+igrid)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     cpu_map(iskip2+igrid)=flag2(igrid)
  end do
  end do

  do ind=1,twotondim
  iskip2=ncoarse+(ind-1)*ngridmax
  ngrid2=0
  do igrid=1,igridmax
     flag2(igrid)=0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              flag2(ngrid2+i)=flag1(iskip2+igrid)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     flag1(iskip2+igrid)=flag2(igrid)
  end do
  end do

  if(hydro)then

#ifdef SOLVERmhd
  do ivar=1,nvar+3
#else
  do ivar=1,nvar
#endif
  do ind=1,twotondim
  iskip2=ncoarse+(ind-1)*ngridmax
  ngrid2=0
  do igrid=1,igridmax
     hilbert_key(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              hilbert_key(ngrid2+i)=real(uold(iskip2+igrid,ivar),kind=qdp)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     uold(iskip2+igrid,ivar)=real(hilbert_key(igrid),kind=8)
  end do
  end do
  end do

  end if

#ifdef RT
  if(rt)then

  do ivar=1,nrtvar
  do ind=1,twotondim
  iskip2=ncoarse+(ind-1)*ngridmax
  ngrid2=0
  do igrid=1,igridmax
     hilbert_key(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              hilbert_key(ngrid2+i)=real(rtuold(iskip2+igrid,ivar),kind=qdp)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     rtuold(iskip2+igrid,ivar)=real(hilbert_key(igrid),kind=8)
  end do
  end do
  end do

  end if
#endif

  if(poisson)then

  do ind=1,twotondim
  iskip2=ncoarse+(ind-1)*ngridmax
  ngrid2=0
  do igrid=1,igridmax
     hilbert_key(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              hilbert_key(ngrid2+i)=real(phi(iskip2+igrid),kind=qdp)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     phi(iskip2+igrid)=real(hilbert_key(igrid),kind=8)
  end do
  end do

  do idim=1,ndim
  do ind=1,twotondim
  iskip2=ncoarse+(ind-1)*ngridmax
  ngrid2=0
  do igrid=1,igridmax
     hilbert_key(igrid)=0.0D0
  end do
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           igrid=istart
           do i=1,ncache
              hilbert_key(ngrid2+i)=real(f(iskip2+igrid,idim),kind=qdp)
              igrid=next(igrid)
           end do
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  do igrid=1,igridmax
     f(iskip2+igrid,idim)=real(hilbert_key(igrid),kind=8)
  end do
  end do
  end do

  end if

  ngrid2=0
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
           if(ibound<=ncpu)then
              headl(ibound,ilevel)=ngrid2+1
              taill(ibound,ilevel)=ngrid2+ncache
           else
              headb(ibound-ncpu,ilevel)=ngrid2+1
              tailb(ibound-ncpu,ilevel)=ngrid2+ncache
           end if
           prev(ngrid2+1)=0
           do i=2,ncache
              prev(ngrid2+i)=ngrid2+i-1
           end do
           do i=1,ncache-1
              next(ngrid2+i)=ngrid2+i+1
           end do
           next(ngrid2+ncache)=0
           ngrid2=ngrid2+ncache
        end if
     end do
  end do
  headf=ngrid2+1
  tailf=ngridmax
  numbf=ngridmax-ngrid2
  prev(headf)=0
  next(tailf)=0
  do i=ngrid2+2,ngridmax
     prev(i)=i-1
  end do
  do i=ngrid2+1,ngridmax-1
     next(i)=i+1
  end do

  do i=1,nlevelmax
     call build_comm(i)
  end do

  ngrid_current=ngrid2
 
end subroutine defrag
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
