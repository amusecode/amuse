! Copyright (c) 2006-2011 IDRIS/CNRS
! Author: Philippe Wautelet (IDRIS/CNRS), wautelet@idris.fr
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.

!TODO: use pack/unpack!

#ifdef IOPROC

#ifndef WITHOUTMPI

subroutine backup_amr_send
  use amr_commons
  use hydro_commons
  use io_parameters
  use pm_commons
  use timer
  use mpi
  implicit none

  integer::idx,ierr
  integer,parameter::tag=TAG_BAK_AMR
  integer::ilevel,ibound,ncache,istart,i,igrid,idim,ind,iskip
  integer,allocatable,dimension(:)::ind_grid,iig
  real(dp),allocatable,dimension(:)::xdp

  if(verbose)write(*,*)'Entering backup_amr_send'

  !-----------------------------------
  ! Output restart dump file = amr.bak
  !-----------------------------------

  call start_timer()

  allocate(iig(19+(3*(ncpu+nboundary)+10)*nlevelmax+3*ncoarse))
  allocate(xdp(19+2*noutput+ncpu))

  iig(1)=ncpu
  iig(2)=ndim
  iig(3)=nx;iig(4)=ny;iig(5)=nz
  iig(6)=nlevelmax
  iig(7)=ngridmax
  iig(11)=nboundary
  iig(14)=ngrid_current
  iig(8)=noutput
  iig(9)=iout
  iig(10)=ifout
  iig(12)=nstep;iig(13)=nstep_coarse
  iig(15)=headf;iig(16)=tailf;iig(17)=numbf;iig(18)=used_mem;iig(19)=used_mem_tot
  iskip=0
  do i=1,nlevelmax
     iig(iskip+20:iskip+19+ncpu)=headl(:,i)
     iig(iskip+20+  ncpu*nlevelmax:iskip+19+ncpu*(nlevelmax+1))=taill(:,i)
     iig(iskip+20+2*ncpu*nlevelmax:iskip+19+2*ncpu*nlevelmax+ncpu)=numbl(:,i)
     iig((i-1)*10+20+3*ncpu*nlevelmax:i*10+19+3*ncpu*nlevelmax)=numbtot(:,i)
     iskip=iskip+ncpu
  end do
  idx=20+(3*ncpu+10)*nlevelmax
  if(simple_boundary)then
     iskip=0
     do i=1,nlevelmax
        idx=20+(3*ncpu+10)*nlevelmax
        iig(iskip+idx:iskip+idx-1+nboundary)=headb(:,i);idx=idx+nboundary*nlevelmax
        iig(iskip+idx:iskip+idx-1+nboundary)=tailb(:,i);idx=idx+nboundary*nlevelmax
        iig(iskip+idx:iskip+idx-1+nboundary)=numbb(:,i);idx=idx+nboundary*nlevelmax
        iskip=iskip+nboundary
     end do
  end if
  iig(idx:idx-1+ncoarse)=son(1:ncoarse)    ;idx=idx+ncoarse
  iig(idx:idx-1+ncoarse)=flag1(1:ncoarse)  ;idx=idx+ncoarse
  iig(idx:idx-1+ncoarse)=cpu_map(1:ncoarse)

  xdp(1)=boxlen
  xdp(2)=t
  xdp(3)=const;xdp(4)=mass_tot_0;xdp(5)=rho_tot
  xdp(6)=omega_m;xdp(7)=omega_l;xdp(8)=omega_k;xdp(9)=omega_b
  xdp(10)=h0;xdp(11)=aexp_ini;xdp(12)=boxlen_ini
  xdp(13)=aexp;xdp(14)=hexp;xdp(15)=aexp_old;xdp(16)=epot_tot_int;xdp(17)=epot_tot_old
  xdp(18)=mass_sph
  xdp(19:18+noutput)=tout(1:noutput)
  xdp(19+noutput:18+2*noutput)=aout(1:noutput)
  if(ordering=='bisection') then

  else
    xdp(19+2*noutput:19+2*noutput+ncpu)=bound_key(0:ncpu)
  endif

  idx = 19+(3*(ncpu+nboundary)+10)*nlevelmax+3*ncoarse


  call MPI_SEND(iig,size(iig),MPI_INTEGER,0,tag,MPI_COMM_IOGROUP,ierr)
  call MPI_SEND(xdp,size(xdp),MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)

  deallocate(iig,xdp)

  call MPI_SEND(dtold,nlevelmax,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
  call MPI_SEND(dtnew,nlevelmax,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)

  ! Write cpu boundaries
  if(ordering=='bisection') then
    call MPI_SEND(bisec_wall,nbinodes,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
    call MPI_SEND(bisec_next,nbinodes,MPI_INTEGER,         0,tag,MPI_COMM_IOGROUP,ierr)
    call MPI_SEND(bisec_indx,2*nbinodes,MPI_INTEGER,       0,tag,MPI_COMM_IOGROUP,ierr)
    call MPI_SEND(bisec_cpubox_min,ncpu*ndim,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
    call MPI_SEND(bisec_cpubox_max,ncpu*ndim,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
  else
     !Send already done (bound_key)
  endif

  ! Write fine levels
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
           allocate(ind_grid(1:ncache),xdp(1:ncache),iig(1:ncache))
           ! Write grid index
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           call MPI_SEND(ind_grid,ncache,MPI_INTEGER,0,tag,MPI_COMM_IOGROUP,ierr)
           ! Write next index
           do i=1,ncache
              iig(i)=next(ind_grid(i))
           end do
           call MPI_SEND(iig,ncache,MPI_INTEGER,0,tag,MPI_COMM_IOGROUP,ierr)
           ! Write prev index
           do i=1,ncache
              iig(i)=prev(ind_grid(i))
           end do
           call MPI_SEND(iig,ncache,MPI_INTEGER,0,tag,MPI_COMM_IOGROUP,ierr)
           ! Write grid center
           do idim=1,ndim
              do i=1,ncache
                 xdp(i)=xg(ind_grid(i),idim)
              end do
              call MPI_SEND(xdp,ncache,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
           end do
           ! Write father index
           do i=1,ncache
              iig(i)=father(ind_grid(i))
           end do
           call MPI_SEND(iig,ncache,MPI_INTEGER,0,tag,MPI_COMM_IOGROUP,ierr)
           ! Write nbor index
           do ind=1,twondim
              do i=1,ncache
                 iig(i)=nbor(ind_grid(i),ind)
              end do
              call MPI_SEND(iig,ncache,MPI_INTEGER,0,tag,MPI_COMM_IOGROUP,ierr)
           end do
           ! Write son index
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ncache
                 iig(i)=son(ind_grid(i)+iskip)
              end do
              call MPI_SEND(iig,ncache,MPI_INTEGER,0,tag,MPI_COMM_IOGROUP,ierr)
           end do
           ! Write cpu map
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ncache
                 iig(i)=cpu_map(ind_grid(i)+iskip)
              end do
              call MPI_SEND(iig,ncache,MPI_INTEGER,0,tag,MPI_COMM_IOGROUP,ierr)
           end do
           ! Write refinement map
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ncache
                 iig(i)=flag1(ind_grid(i)+iskip)
              end do
              call MPI_SEND(iig,ncache,MPI_INTEGER,0,tag,MPI_COMM_IOGROUP,ierr)
           end do
           deallocate(xdp,iig,ind_grid)
        end if
     end do
  end do

  call stop_timer('AMR I/O processes backup',writing=.true.)

end subroutine backup_amr_send

subroutine backup_amr_recv
  use amr_commons
  use hydro_commons
  use io_commons
  use pm_commons
  use mpi
  implicit none

  integer::count,i,idx,ierr,src
  integer,parameter::tag=TAG_BAK_AMR
  integer,dimension(MPI_STATUS_SIZE)::status
  integer::ilun
  integer::ilevel,ibound,ncache,idim,ind,iskip
  integer,allocatable,dimension(:)::bisec_int,ind_grid,iig
  logical,allocatable,dimension(:)::list_recv
  real(dp),allocatable,dimension(:)::bisec_dp,xdp
  character(LEN=5)::cpuchar,iochar,nchar
  character(LEN=MAXLINE)::filename

  if(verbose)write(*,*)'Entering backup_amr_recv'
  allocate(list_recv(ncpu_iogroup-1))
  list_recv(:)=.false.

  count=0
  do while(count<ncpu_iogroup-1)
     ! Allocate receive buffers
     idx=19+(3*(ncpu+nboundary)+10)*nlevelmax+3*ncoarse
     allocate(iig(idx))
     allocate(xdp(19+2*noutput+ncpu))

     ! Select a source
     call MPI_RECV(iig,idx,MPI_INTEGER,MPI_ANY_SOURCE,tag,MPI_COMM_IOGROUP,status,ierr)
     src=status(MPI_SOURCE)
     if(list_recv(src).EQV..false.)then
        list_recv(src)=.true.
     else
        print *,'Error: unexpected message received by ',myid_world
        call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
     end if

     call MPI_RECV(xdp,19+2*noutput+ncpu,MPI_DOUBLE_PRECISION,src,tag, &
                   MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)

     ! Generate filename
     ilun=myid_world+10
     call title(count_bak,nchar)
     call title(iogroup2comp(src),cpuchar)
     call title(myid_io,iochar)
     filename='ionode_'//TRIM(iochar)//'/process_'//TRIM(cpuchar)//'/'
     filename=TRIM(filename)//'amr_'//TRIM(nchar)//'.out'
     filename=TRIM(filename)//TRIM(cpuchar)
     nbfiles=nbfiles+1
     filelist(nbfiles)=trim(filename)
     open(unit=ilun,file=trim(scratchdir)//trim(filename),status="replace",form="unformatted",action="write",iostat=ierr)
     if(ierr/=0)then
        print *,'Error: open file failed in backup_amr_recv'
        call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
     end if

     !-----------------------------------
     ! Output restart dump file = amr.bak
     !-----------------------------------

     ! Write grid variables
     write(ilun)iig(1) !ncpu
     write(ilun)iig(2) !ndim
     write(ilun)iig(3),iig(4),iig(5) !nx,ny,nz
     write(ilun)iig(6) !nlevelmax
     write(ilun)iig(7) !ngridmax
     write(ilun)iig(11) !nboundary
     write(ilun)iig(14) !ngrid_current
     write(ilun)xdp(1) !boxlen
     ! Write time variables
     write(ilun)iig(8),iig(9),iig(10) !noutput,iout,ifout
     write(ilun)xdp(19:18+noutput)   !tout
     write(ilun)xdp(19+noutput:18+2*noutput) !aout
     write(ilun)xdp(2) !t
     call MPI_RECV(dtold,iig(6),MPI_DOUBLE_PRECISION,src,tag, &
                   MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
     write(ilun)dtold(1:nlevelmax)
     call MPI_RECV(dtnew,iig(6),MPI_DOUBLE_PRECISION,src,tag, &
                   MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
     write(ilun)dtnew(1:nlevelmax)
     write(ilun)iig(12),iig(13) !nstep,nstep_coarse
     write(ilun)xdp(3),xdp(4),xdp(5) !const,mass_tot_0,rho_tot
       !omega_m,omega_l,omega_k,omega_b,h0,aexp_ini,boxlen_ini
     write(ilun)xdp(6),xdp(7),xdp(8),xdp(9),xdp(10),xdp(11),xdp(12)
     write(ilun)xdp(13),xdp(14),xdp(15),xdp(16),xdp(17) !aexp,hexp,aexp_old,epot_tot_int,epot_tot_old
     write(ilun)xdp(18) !mass_sph
     ! Write levels variables
     write(ilun)iig(20:19+ncpu*nlevelmax) !headl
     write(ilun)iig(20+  ncpu*nlevelmax:19+2*ncpu*nlevelmax) !taill
     write(ilun)iig(20+2*ncpu*nlevelmax:19+3*ncpu*nlevelmax) !numbl
     iskip=0
     do i=1,nlevelmax
        idx=20+2*ncpu*nlevelmax
        numblio(1:ncpu,i,src)=iig(iskip+idx:iskip+idx-1+ncpu)
        iskip=iskip+ncpu
     end do
     write(ilun)iig(20+3*ncpu*nlevelmax:19+(3*ncpu+10)*nlevelmax) !numbtot
     idx=20+(3*ncpu+10)*nlevelmax
     ! Write boundary linked list
     if(simple_boundary)then
        write(ilun)iig(idx:idx-1+nboundary*nlevelmax);idx=idx+nboundary*nlevelmax !headb
        write(ilun)iig(idx:idx-1+nboundary*nlevelmax);idx=idx+nboundary*nlevelmax !tailb
        write(ilun)iig(idx:idx-1+nboundary*nlevelmax)                             !numbb
        iskip=0
        do i=1,nlevelmax
           numbb(1:nboundary,i)=iig(iskip+idx:iskip+idx-1+nboundary)
           iskip=iskip+nboundary
        end do
     end if

     ! Write free memory
     write(ilun)iig(15),iig(16),iig(17),iig(18),iig(19) !headf,tailf,numbf,used_mem,used_mem_tot

    ! Write cpu boundaries
    write(ilun)ordering
    if(ordering=='bisection') then
     allocate(bisec_int(nbinodes*2),bisec_dp(max(nbinodes,ncpu*ndim)))
     call MPI_RECV(bisec_dp,nbinodes,MPI_DOUBLE_PRECISION,src,tag, &
                   MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
     write(ilun)bisec_dp(1:nbinodes)  !bisec_wall(1:nbinodes)
     call MPI_RECV(bisec_int,2*nbinodes,MPI_INTEGER,src,tag, &
                   MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
     write(ilun)bisec_int(1:2*nbinodes) !bisec_next(1:nbinodes,1:2)
     call MPI_RECV(bisec_int,nbinodes,MPI_INTEGER,src,tag, &
                   MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
     write(ilun)bisec_int(1:nbinodes) !bisec_indx(1:nbinodes)
     call MPI_RECV(bisec_dp,ncpu*ndim,MPI_DOUBLE_PRECISION,src,tag, &
                   MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
     write(ilun)bisec_dp(1:ncpu*ndim)  !bisec_cpubox_min(1:ncpu,1:ndim)
     call MPI_RECV(bisec_dp,ncpu*ndim,MPI_DOUBLE_PRECISION,src,tag, &
                   MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
     write(ilun)bisec_dp(1:ncpu*ndim)  !bisec_cpubox_max(1:ncpu,1:ndim)
     deallocate(bisec_int,bisec_dp)
    else
      write(ilun)xdp(19+2*noutput:19+2*noutput+ncpu) !bound_key
    endif

     ! Write coarse level
     write(ilun)iig(idx:idx-1+ncoarse);idx=idx+ncoarse !son
     write(ilun)iig(idx:idx-1+ncoarse);idx=idx+ncoarse !flag1
     write(ilun)iig(idx:idx-1+ncoarse)                 !cpu_map

     deallocate(iig,xdp)

    ! Write fine levels
    do ilevel=1,nlevelmax
      do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
          ncache=numblio(ibound,ilevel,src)
        else
          ncache=numbb(ibound-ncpu,ilevel)
        end if
        if(ncache>0)then
          allocate(ind_grid(1:ncache),xdp(1:ncache),iig(1:ncache))
          ! Write grid index
          call MPI_RECV(ind_grid,ncache,MPI_INTEGER,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
           write(ilun)ind_grid
           ! Write next index
           call MPI_RECV(iig,ncache,MPI_INTEGER,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
           write(ilun)iig
           ! Write prev index
           call MPI_RECV(iig,ncache,MPI_INTEGER,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
           write(ilun)iig
           ! Write grid center
           do idim=1,ndim
             call MPI_RECV(xdp,ncache,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
             write(ilun)xdp
           end do
           ! Write father index
           call MPI_RECV(iig,ncache,MPI_INTEGER,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
           write(ilun)iig
           ! Write nbor index
           do ind=1,twondim
             call MPI_RECV(iig,ncache,MPI_INTEGER,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
             write(ilun)iig
           end do
           ! Write son index
           do ind=1,twotondim
             call MPI_RECV(iig,ncache,MPI_INTEGER,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
             write(ilun)iig
           end do
           ! Write cpu map
           do ind=1,twotondim
             call MPI_RECV(iig,ncache,MPI_INTEGER,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
             write(ilun)iig
           end do
           ! Write refinement map
           do ind=1,twotondim
             call MPI_RECV(iig,ncache,MPI_INTEGER,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
             write(ilun)iig
           end do
           deallocate(xdp,iig,ind_grid)
         end if
       end do
     end do

     close(ilun)

     count=count+1
  end do

  deallocate(list_recv)

end subroutine backup_amr_recv

subroutine output_info_send
  use amr_commons
  use hydro_commons
  use io_parameters
  use pm_commons
  use mpi
  implicit none

  integer,parameter::tag=TAG_OUT_INF
  integer::nx_loc,ny_loc,nz_loc,ierr
  real(dp)::scale
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8),dimension(13)::msg

  if(verbose)write(*,*)'Entering output_info_send'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Local constants
  nx_loc=nx; ny_loc=ny; nz_loc=nz
  if(ndim>0)nx_loc=(icoarse_max-icoarse_min+1)
  if(ndim>1)ny_loc=(jcoarse_max-jcoarse_min+1)
  if(ndim>2)nz_loc=(kcoarse_max-kcoarse_min+1)
  scale=boxlen/dble(nx_loc)

  msg(1)  = scale
  msg(2)  = t
  msg(3)  = aexp
  msg(4)  = omega_m
  msg(5)  = omega_l
  msg(6)  = omega_k
  msg(7)  = omega_b
  msg(8)  = scale_l
  msg(9)  = scale_d
  msg(10) = scale_t
  msg(11) = h0
  msg(12) = nstep_coarse+0.1 !Not very clean but useful (one less message)
  msg(13) = ndomain

  call MPI_SEND(msg,size(msg),MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)

  if(ordering=='bisection') then
    call MPI_SEND(bisec_cpubox_min,ncpu*ndim,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
    call MPI_SEND(bisec_cpubox_max,ncpu*ndim,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
    call MPI_SEND(bisec_cpu_load,ncpu,       MPI_INTEGER,         0,tag,MPI_COMM_IOGROUP,ierr)
  else
    call MPI_SEND(bound_key(0:ncpu),ncpu+1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
  endif
end subroutine output_info_send

subroutine output_info_recv
  use amr_commons
  use hydro_commons
  use io_commons
  use pm_commons
  use mpi
  implicit none

  integer::ilun,icpu,ierr,idom
  integer,parameter::tag=TAG_OUT_INF
  character(LEN=MAXLINE)::filename
  character(LEN=5)::nchar,iochar
  real(kind=8),dimension(13)::msg

  if(verbose)write(*,*)'Entering output_info_recv'

  call MPI_RECV(msg,size(msg),MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
  ndomain = msg(13)

  ! Generate filename
  ilun=100
  call title(count_bak,nchar)
  call title(myid_io,iochar)
  filename='ionode_'//TRIM(iochar)//'/process_00001/'
  filename=TRIM(filename)//'info_'//TRIM(nchar)//'.txt'
  nbfiles=nbfiles+1
  filelist(nbfiles)=trim(filename)
  open(unit=ilun,file=trim(scratchdir)//trim(filename),status="replace",form="formatted",action="write",iostat=ierr)
  if(ierr/=0)then
     print *,'Error: open file failed in output_info_recv'
print *,'filename=',trim(filename)
     call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
  end if

  ! Write run parameters
  write(ilun,'("ncpu        =",I11)')ncpu
  write(ilun,'("ndim        =",I11)')ndim
  write(ilun,'("levelmin    =",I11)')levelmin
  write(ilun,'("levelmax    =",I11)')nlevelmax
  write(ilun,'("ngridmax    =",I11)')ngridmax
  write(ilun,'("nstep_coarse=",I11)')int(msg(12)) !nstep_coarse
  write(ilun,*)

  ! Write physical parameters
  write(ilun,'("boxlen      =",E23.15)')msg(1)  !scale
  write(ilun,'("time        =",E23.15)')msg(2)  !t
  write(ilun,'("aexp        =",E23.15)')msg(3)  !aexp
  write(ilun,'("H0          =",E23.15)')msg(11) !h0
  write(ilun,'("omega_m     =",E23.15)')msg(4)  !omega_m
  write(ilun,'("omega_l     =",E23.15)')msg(5)  !omega_l
  write(ilun,'("omega_k     =",E23.15)')msg(6)  !omega_k
  write(ilun,'("omega_b     =",E23.15)')msg(7)  !omega_b
  write(ilun,'("unit_l      =",E23.15)')msg(8)  !scale_l
  write(ilun,'("unit_d      =",E23.15)')msg(9)  !scale_d
  write(ilun,'("unit_t      =",E23.15)')msg(10) !scale_t
  write(ilun,*)

  ! Write ordering information
  write(ilun,'("ordering type=",A80)')ordering
  if(ordering=='bisection') then
    allocate(bisec_cpubox_min(ncpu,ndim))
    allocate(bisec_cpubox_max(ncpu,ndim))
    allocate(bisec_cpu_load(ncpu))
    call MPI_RECV(bisec_cpubox_min,ncpu*ndim,MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(bisec_cpubox_max,ncpu*ndim,MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(bisec_cpu_load,ncpu,MPI_INTEGER,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
    do icpu=1,ncpu
      ! write 2*ndim floats for cpu bound box
      write(ilun,'(E23.15)')bisec_cpubox_min(icpu,:),bisec_cpubox_max(icpu,:)
      ! write 1 float for cpu load
      write(ilun,'(E23.15)')dble(bisec_cpu_load(icpu))
    end do
    deallocate(bisec_cpubox_min,bisec_cpubox_max,bisec_cpu_load)
  else
    allocate(bound_key(0:ncpu))
    call MPI_RECV(bound_key(0:ncpu),ncpu+1,MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
    write(ilun,'("   DOMAIN   ind_min                 ind_max")')
    do idom=1,ndomain
      write(ilun,'(I8,1X,E23.15,1X,E23.15)')idom,bound_key(idom-1),bound_key(idom)
    end do
    deallocate(bound_key)
  endif

 close(ilun)

end subroutine output_info_recv

#endif

#endif
