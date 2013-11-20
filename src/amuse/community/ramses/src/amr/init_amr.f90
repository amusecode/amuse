subroutine init_amr
  use amr_commons
  use hydro_commons
  use pm_commons  
  use poisson_commons
  use bisection
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'  
#endif
  integer::i,idim,ncell,iskip,ind,ncache,ilevel,ibound,nboundary2
  integer::ncpu2,ndim2,nx2,ny2,nz2,ngridmax2,nlevelmax2
  integer::noutput2,iout2,ifout2,ilun,info
  integer::ix,iy,iz,ix_max,iy_max,iz_max,nxny,nx_loc
  real(dp)::mass_sph2 
  integer,dimension(:),allocatable::ind_grid,iig,pos,grid
  real(dp),dimension(1:MAXOUT)::aout2=1.1d0 
  real(dp),dimension(1:MAXOUT)::tout2=0.0d0 
  real(dp),dimension(:),allocatable::xxg
  integer ,dimension(1:nvector)::c
  real(dp),dimension(1:nvector,1:ndim)::x
  real(qdp),dimension(1:nvector)::order_min,order_max
  logical::ok
  real(dp)::dx_loc,scale
  character(LEN=128)::ordering2
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp),allocatable,dimension(:)::bxmin,bxmax

  if(verbose.and.myid==1)write(*,*)'Entering init_amr'

  ! Constants
  ncoarse=nx*ny*nz
  ncell=ncoarse+twotondim*ngridmax
  nxny=nx*ny
  ix_max=0; iy_max=0; iz_max=0
  if(ndim>0)ix_max=1 
  if(ndim>1)iy_max=1
  if(ndim>2)iz_max=1
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)

  ! Initial time step for each level
  dtold=0.0D0
  dtnew=0.0D0

  ! Allocate AMR cell-based arrays
  allocate(flag1(0:ncell)) ! Note: starting from 0
  allocate(flag2(0:ncell)) ! Note: starting from 0
  allocate(son  (1:ncell)) ! Son index
  flag1=0; flag2=0; son=0

  ! Allocate MPI cell-based arrays
  allocate(cpu_map    (1:ncell)) ! Cpu map
  allocate(cpu_map2   (1:ncell)) ! New cpu map for load balance
  allocate(hilbert_key(1:ncell)) ! Ordering key
  cpu_map=0; cpu_map2=0; hilbert_key=0.0d0

  ! Bisection ordering: compute array boundaries and
  ! allocate arrays if needed
  nbilevelmax=ceiling(log(dble(ncpu))/log(2.0))
  nbinodes=2**(nbilevelmax+1)-1
  nbileafnodes=2**nbilevelmax
  bisec_nres=2**(nlevelmax+1)
  bisec_res=scale/dble(bisec_nres)

  if(ordering=='bisection') then
    ! allocate bisection tree structure
    allocate(bisec_wall(1:nbinodes))
    allocate(bisec_next(1:nbinodes,1:2))
    allocate(bisec_indx(1:nbinodes))
    bisec_wall=0.0d0; bisec_next=0; bisec_indx=0; bisec_root=0
    ! allocate some other bisection stuff
    allocate(bisec_cpubox_min (1:ncpu,1:ndim))
    allocate(bisec_cpubox_max (1:ncpu,1:ndim))
    allocate(bisec_cpubox_min2(1:ncpu,1:ndim))
    allocate(bisec_cpubox_max2(1:ncpu,1:ndim))
    allocate(bisec_cpu_load(1:ncpu))
    bisec_cpubox_min=0;  bisec_cpubox_max=0;
    bisec_cpubox_min2=0; bisec_cpubox_max2=0;
    bisec_cpu_load=0;
    ! allocate histograms
    allocate(bisec_hist(1:nbileafnodes,1:bisec_nres))
    allocate(bisec_hist_bounds(1:(nbileafnodes+1)))
    allocate(new_hist_bounds  (1:(nbileafnodes+1)))
    allocate(bisec_ind_cell(1:ncell))    ! big array
    allocate(cell_level    (1:ncell))    ! big array
    bisec_hist=0
    bisec_hist_bounds=0; new_hist_bounds=0
    bisec_ind_cell=0; cell_level=0
  end if

 bisection_or_ordering:if(ordering /= 'bisection') then ! use usual ordering machinery

    ! Cpu boundaries in chosen ordering
    ndomain=ncpu*overload
    allocate(bound_key (0:ndomain))
    allocate(bound_key2(0:ndomain))

    ! Compute minimum and maximum ordering key
    dx_loc=scale
    x(1,1)=0.5*scale
#if NDIM>1
    x(1,2)=0.5*scale
#endif
#if NDIM>2
    x(1,3)=0.5*scale
#endif
    call cmp_minmaxorder(x,order_min,order_max,dx_loc,1)
    order_all_min=order_min(1)
    order_all_max=order_max(1)
    do iz=kcoarse_min,kcoarse_max
       do iy=jcoarse_min,jcoarse_max
          do ix=icoarse_min,icoarse_max
             ind=1+ix+iy*nx+iz*nxny
             x(1,1)=(dble(ix)+0.5d0-dble(icoarse_min))*scale
#if NDIM>1
             x(1,2)=(dble(iy)+0.5d0-dble(jcoarse_min))*scale
#endif
#if NDIM>2
             x(1,3)=(dble(iz)+0.5d0-dble(kcoarse_min))*scale
#endif
             call cmp_minmaxorder(x,order_min,order_max,dx_loc,1)
             order_all_min=min(order_all_min,order_min(1))
             order_all_max=max(order_all_max,order_max(1))
          end do
       end do
    end do

    ! Set initial cpu boundaries
    do i=0,ndomain-1
#ifdef QUADHILBERT
       bound_key(i)=order_all_min+real(i,16)/real(ndomain,16)* &
            & (order_all_max-order_all_min)
#else
       bound_key(i)=order_all_min+real(i,8)/real(ndomain,8)* &
            & (order_all_max-order_all_min)
#endif
    end do
    bound_key(ndomain)=order_all_max

      else ! Init bisection balancing

       call build_bisection(update=.false.)

  end if bisection_or_ordering

  ! Compute coarse cpu map
  do iz=kcoarse_min,kcoarse_max
  do iy=jcoarse_min,jcoarse_max
  do ix=icoarse_min,icoarse_max
     ind=1+ix+iy*nx+iz*nxny
     x(1,1)=(dble(ix)+0.5d0-dble(icoarse_min))*scale
#if NDIM>1
     x(1,2)=(dble(iy)+0.5d0-dble(jcoarse_min))*scale
#endif
#if NDIM>2
     x(1,3)=(dble(iz)+0.5d0-dble(kcoarse_min))*scale
#endif
     call cmp_cpumap(x,c,1)
     cpu_map(ind)=c(1)
  end do
  end do
  end do

  ! Allocate linked list for each level
  allocate(headl(1:ncpu,1:nlevelmax))
  allocate(taill(1:ncpu,1:nlevelmax))
  allocate(numbl(1:ncpu,1:nlevelmax))
  allocate(numbtot(1:10,1:nlevelmax))
  headl=0    ! Head grid in the level
  taill=0    ! Tail grid in the level
  numbl=0    ! Number of grids in the level
  numbtot=0  ! Total number of grids in the level

  ! Allocate communicators
  allocate(active(1:nlevelmax))
  allocate(emission(1:ncpu,1:nlevelmax))
  allocate(reception(1:ncpu,1:nlevelmax))
  do ilevel=1,nlevelmax
     active(ilevel)%ngrid=0
     do i=1,ncpu
        emission (i,ilevel)%ngrid=0
        emission (i,ilevel)%npart=0
        reception(i,ilevel)%ngrid=0
        reception(i,ilevel)%npart=0
     end do
  end do
  ! Allocate lookup array for multigrid fine
  if(poisson) allocate(lookup_mg(1:ngridmax))

  ! Allocate physical boundary for each level
  allocate(headb   (1:MAXBOUND,1:nlevelmax))
  allocate(tailb   (1:MAXBOUND,1:nlevelmax))
  allocate(numbb   (1:MAXBOUND,1:nlevelmax))
  allocate(boundary(1:MAXBOUND,1:nlevelmax))
  do i=1,MAXBOUND
     do ilevel=1,nlevelmax
        headb   (i,ilevel)=0       ! Head grid in boundary
        tailb   (i,ilevel)=0       ! Tail grid in boundary
        numbb   (i,ilevel)=0       ! Number of grids in boundary
        boundary(i,ilevel)%ngrid=0 ! Communicators
     end do
  end do

  ! Allocate grid center coordinates
  allocate(xg(1:ngridmax,1:ndim))
  xg=0.0D0

  ! Allocate tree arrays
  allocate(father(1:ngridmax))
  allocate(nbor  (1:ngridmax,1:twondim))
  allocate(next  (1:ngridmax))
  allocate(prev  (1:ngridmax))
  father=0; nbor=0; next=0; prev=0

  ! Allocate pointer to particles linked lists
  if(pic)then
     allocate(headp(1:ngridmax))
     allocate(tailp(1:ngridmax))
     allocate(numbp(1:ngridmax))
     headp=0; tailp=0; numbp=0
  endif

  ! Initialize AMR grid linked list
  do i=1,ngridmax-1
     next(i)=i+1
  end do
  do i=2,ngridmax
     prev(i)=i-1
  end do
  headf=1          ! Pointer to first grid in free memory
  tailf=ngridmax   ! Pointer to last grid in free memory
  prev(headf)=0; next(tailf)=0
  numbf=ngridmax   ! Number of grids in free memory
  used_mem=ngridmax-numbf

  !----------------------------
  ! Read amr file for a restart
  !----------------------------
  if(nrestart>0)then
     ilun=myid+10
     call title(nrestart,nchar)
     fileloc='output_'//TRIM(nchar)//'/amr_'//TRIM(nchar)//'.out'
     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)
     inquire(file=fileloc, exist=ok)
     if(.not. ok)then
        write(*,*)'Restart failed:'
        write(*,*)'File '//TRIM(fileloc)//' not found'
        call clean_stop
     end if
     if(debug)write(*,*)'amr.tmp opened for processor ',myid
     open(unit=ilun,file=fileloc,form='unformatted')
     ! Read grid variables
     read(ilun)ncpu2
     read(ilun)ndim2
     read(ilun)nx2,ny2,nz2
     read(ilun)nlevelmax2
     read(ilun)ngridmax2
     read(ilun)nboundary2
     read(ilun)ngrid_current
     read(ilun)boxlen
     if(ncpu2.ne.ncpu)then
        write(*,*)'Number of processes not compatible'
        write(*,*)'ncpu should be set equal to',ncpu2
        call clean_stop
     end if
     ! Read time variables
     read(ilun)noutput2,iout2,ifout2
     if(noutput2>MAXOUT)then
       write(*,*) 'Error: noutput>MAXOUT'
       call clean_stop
     end if
     read(ilun)tout2(1:noutput2)
     read(ilun)aout2(1:noutput2)
     ! Check compatibility with current parameters
     if((ndim2.ne.ndim).or.(nx2.ne.nx).or.(ny2.ne.ny).or.(nz2.ne.nz).or.&
          & (nboundary2.ne.nboundary).or.(nlevelmax2>nlevelmax).or.&
          & (ngrid_current>ngridmax).or.(noutput2>noutput) )then
        write(*,*)'File amr.tmp is not compatible with namelist'
        write(*,*)'         ndim   nx   ny   nz nlevelmax noutput   ngridmax nboundary'
        write(*,'("amr.tmp  =",4(I4,1x),5x,I4,4x,I4,3x,I8)')&
             & ndim2,nx2,ny2,nz2,nlevelmax2,noutput2,ngrid_current,nboundary2
        write(*,'("namelist =",4(I4,1x),5x,I4,4x,I4,3x,I8)')&
             & ndim ,nx ,ny ,nz ,nlevelmax ,noutput, ngridmax     ,nboundary
        if(myid==1)write(*,*)'Restart failed'
        call clean_stop 
     end if
     ! Old output times
     tout(1:noutput2)=tout2(1:noutput2)
     aout(1:noutput2)=aout2(1:noutput2)
     iout=iout2
     ifout=ifout2
     read(ilun)t
     read(ilun)dtold(1:nlevelmax2)
     read(ilun)dtnew(1:nlevelmax2)
     read(ilun)nstep,nstep_coarse
     nstep_coarse_old=nstep_coarse
     read(ilun)const,mass_tot_0,rho_tot
     read(ilun)omega_m,omega_l,omega_k,omega_b,h0,aexp_ini,boxlen_ini
     read(ilun)aexp,hexp,aexp_old,epot_tot_int,epot_tot_old
     if(cosmo)then
        read(ilun)mass_sph
     else
        read(ilun)mass_sph2
     endif
     if(myid==1)write(*,*)'Restarting at t=',t,' nstep_coarse=',nstep_coarse

     ! Compute movie frame number if applicable
     if(imovout>0) then
        do i=2,imovout
           if(aendmov>0)then
              if(aexp>amovout(i-1).and.aexp<amovout(i)) then
                 imov=i
              endif
           else
              if(t>tmovout(i-1).and.t<tmovout(i)) then
                 imov=i
              endif
           endif
        enddo
        if(aendmov>0)then
           if(myid==1)write(*,*) "Frame number, aexp ",imov, amovout(imov)
        else
           if(myid==1)write(*,*) "Frame number, t ",imov, tmovout(imov)
        endif
     endif

     ! Read levels variables
     read(ilun)headl(1:ncpu,1:nlevelmax2)
     read(ilun)taill(1:ncpu,1:nlevelmax2)
     read(ilun)numbl(1:ncpu,1:nlevelmax2)
     read(ilun)numbtot(1:10,1:nlevelmax2)
     ! Read boundary linked list
     if(simple_boundary)then
        read(ilun)headb(1:nboundary,1:nlevelmax2)
        read(ilun)tailb(1:nboundary,1:nlevelmax2)
        read(ilun)numbb(1:nboundary,1:nlevelmax2)
     end if
     ! Read free memory
     read(ilun)headf,tailf,numbf,used_mem,used_mem_tot
     headf=ngrid_current+1
     tailf=ngridmax
     numbf=ngridmax-ngrid_current
     prev(headf)=0
     next(tailf)=0
     ! Read cpu boundaries
     read(ilun)ordering2
     if(ordering2.ne.ordering)then
        if(myid==1)write(*,*)'Ordering is uncompatible'
        call clean_stop
     endif
     if(ordering=='bisection') then
        read(ilun)bisec_wall(1:nbinodes)
        read(ilun)bisec_next(1:nbinodes,1:2)
        read(ilun)bisec_indx(1:nbinodes)
        read(ilun)bisec_cpubox_min(1:ncpu,1:ndim)
        read(ilun)bisec_cpubox_max(1:ncpu,1:ndim)
     else
        read(ilun)bound_key(0:ndomain)
     endif
     ! Read coarse level
     read(ilun)son(1:ncoarse)
     read(ilun)flag1(1:ncoarse)
     read(ilun)cpu_map(1:ncoarse)
     ! Read fine levels
     do ilevel=1,nlevelmax2
        do ibound=1,nboundary+ncpu
           if(ibound<=ncpu)then
              ncache=numbl(ibound,ilevel)
           else
              ncache=numbb(ibound-ncpu,ilevel)
           end if
           if(ncache>0)then
              allocate(ind_grid(1:ncache))
              allocate(xxg(1:ncache))
              allocate(iig(1:ncache))
              allocate(pos(1:ncache))
              allocate(grid(1:ncache))
              ! Read grid index
              read(ilun)ind_grid
              ! Read next index
              read(ilun)iig
              do i=1,ncache
                 next(ind_grid(i))=iig(i)
              end do
              ! Read prev index
              read(ilun)iig
              do i=1,ncache
                 prev(ind_grid(i))=iig(i)
              end do              
              ! Read grid center
              do idim=1,ndim
                 read(ilun)xxg
                 do i=1,ncache
                    xg(ind_grid(i),idim)=xxg(i)
                 end do
              end do
              ! Read father index
              read(ilun)iig
              if(ngridmax.ne.ngridmax2.and.ilevel>1)then
                 do i=1,ncache
                    pos(i)=(iig(i)-ncoarse-1)/ngridmax2
                 end do
                 do i=1,ncache
                    grid(i)=iig(i)-ncoarse-pos(i)*ngridmax2
                 end do
                 do i=1,ncache
                    iig(i)=ncoarse+pos(i)*ngridmax+grid(i)
                 end do
              end if
              do i=1,ncache
                 father(ind_grid(i))=iig(i)
              end do
              ! Read nbor index
              do ind=1,twondim
                 read(ilun)iig
                 if(ngridmax.ne.ngridmax2.and.ilevel>1)then
                    do i=1,ncache
                       pos(i)=(iig(i)-ncoarse-1)/ngridmax2
                    end do
                    do i=1,ncache
                       grid(i)=iig(i)-ncoarse-pos(i)*ngridmax2
                    end do
                    do i=1,ncache
                       iig(i)=ncoarse+pos(i)*ngridmax+grid(i)
                    end do
                 end if
                 do i=1,ncache
                    nbor(ind_grid(i),ind)=iig(i)
                 end do
              end do
              ! Read son index
              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax
                 read(ilun)iig
                 do i=1,ncache
                    son(ind_grid(i)+iskip)=iig(i)
                 end do
              end do
              ! Read cpu map
              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax
                 read(ilun)iig
                 do i=1,ncache
                    cpu_map(ind_grid(i)+iskip)=iig(i)
                 end do
              end do
              ! Read refinement map
              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax
                 read(ilun)iig
                 do i=1,ncache
                    flag1(ind_grid(i)+iskip)=iig(i)
                 end do
              end do
              deallocate(xxg,iig,pos,grid,ind_grid)
           end if
        end do
     end do
     close(ilun)

#ifndef WITHOUTMPI
     if(debug)write(*,*)'amr.tmp read for processor ',myid
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(verbose)write(*,*)'AMR backup files read completed'

     ! Build communicators
     do ilevel=1,nlevelmax
        call build_comm(ilevel)
     end do

  end if
  
end subroutine init_amr


