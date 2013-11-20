! Copyright (c) 2006-2011 IDRIS/CNRS
! Author: Philippe Wautelet (IDRIS/CNRS), wautelet@idris.fr
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.

#ifdef IOPROC

#ifndef WITHOUTMPI

subroutine backup_part_send
  use amr_commons
  use pm_commons
  use io_parameters
  use mpi
  use timer
  implicit none

  integer,parameter::tag=TAG_BAK_PAR
  integer::i,idim,ipart,ierr
  real(dp),allocatable,dimension(:)::xdp
  integer ,allocatable,dimension(:)::ii
  integer ,allocatable,dimension(:)::ll
  integer ,allocatable,dimension(:)::ok

  if(verbose)write(*,*)'Entering backup_part_send'

  call start_timer()

  allocate(ii(3+IRandNumSize),xdp(2))
  ii(1:IRandNumSize)=localseed(:)
  ii(IRandNumSize+1)=npart
  ii(IRandNumSize+2)=nstar_tot
  ii(IRandNumSize+3)=nsink
  call MPI_SEND(ii,size(ii),MPI_INTEGER,0,tag,MPI_COMM_IOGROUP,ierr)
  xdp(1)=mstar_tot
  xdp(2)=mstar_lost
  call MPI_SEND(xdp,2,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
  deallocate(ii,xdp)

  allocate(xdp(1:npart))
  allocate(ok (1:npart))
  ! Write position
  do idim=1,ndim
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           xdp(ipart)=xp(i,idim)
           ok(ipart)=i
        end if
     end do
     call MPI_SEND(xdp,npart,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
  end do
  ! Write velocity
  do idim=1,ndim
     ipart=0
     do i=1,npart
        ipart=ipart+1
        xdp(ipart)=vp(ok(i),idim)
     end do
     call MPI_SEND(xdp,npart,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
  end do
  ! Write mass
  ipart=0
  do i=1,npart
     ipart=ipart+1
     xdp(ipart)=mp(ok(i))
  end do
  call MPI_SEND(xdp,npart,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
  ! Write identity
  allocate(ii(1:npart))
  ipart=0
  do i=1,npart
     ipart=ipart+1
     ii(ipart)=idp(ok(i))
  end do
  call MPI_SEND(ii,npart,MPI_INTEGER,0,tag,MPI_COMM_IOGROUP,ierr)
  deallocate(ii)
  ! Write level
  allocate(ll(1:npart))
  ipart=0
  do i=1,npart
     ipart=ipart+1
     ll(ipart)=levelp(ok(i))
  end do
  call MPI_SEND(ll,npart,MPI_INTEGER,0,tag,MPI_COMM_IOGROUP,ierr)
  deallocate(ll)
  if(star.or.sink)then
     ! Write birth epoch
     ipart=0
     do i=1,npart
        ipart=ipart+1
        xdp(ipart)=tp(ok(i))
     end do
     call MPI_SEND(xdp,npart,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
     if(metal)then
        ! Write metallicity
        ipart=0
        do i=1,npart
           ipart=ipart+1
           xdp(ipart)=zp(ok(i))
        end do
     call MPI_SEND(xdp,npart,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
     endif
  end if
  deallocate(xdp)
  deallocate(ok)
  ! Write sink properties
  if(sink.and.nsink>0)then
    call MPI_SEND(msink(1:nsink),nsink,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
    call MPI_SEND(xsink(1:nsink,1:ndim),nsink*ndim,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
  endif

  call stop_timer('part I/O processes backup',writing=.true.)
end subroutine backup_part_send

subroutine backup_part_recv
  use amr_commons
  use io_commons
  use pm_commons
  use mpi
  implicit none

  integer,parameter::tag=TAG_BAK_PAR
  integer,dimension(MPI_STATUS_SIZE)::status
  logical,allocatable,dimension(:)::list_recv
  integer::idim,ilun,count,src,ierr
  real(dp),allocatable,dimension(:)::xdp
  integer ,allocatable,dimension(:)::ii
  integer ,allocatable,dimension(:)::ll
  character(LEN=MAXLINE)::filename
  character(LEN=5)::cpuchar,iochar,nchar

  if(verbose)write(*,*)'Entering backup_part_recv'

  allocate(list_recv(ncpu_iogroup-1))
  list_recv(:)=.false.

  count=0
  do while(count<ncpu_iogroup-1)
     ! Select a source
     allocate(ii(3+IRandNumSize),xdp(2))
     call MPI_RECV(ii,size(ii),MPI_INTEGER,MPI_ANY_SOURCE,tag,MPI_COMM_IOGROUP,status,ierr)

     src=status(MPI_SOURCE)
     if(list_recv(src).EQV..false.)then
        list_recv(src)=.true.
     else
        print *,'Error: unexpected message received by ',myid_world
        call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
     end if

     ! Generate filename
     ilun=myid_world+10
     call title(count_bak,nchar)
     call title(iogroup2comp(src),cpuchar)
     call title(myid_io,iochar)
     filename='ionode_'//TRIM(iochar)//'/process_'//TRIM(cpuchar)//'/'
     filename=TRIM(filename)//'part_'//TRIM(nchar)//'.out'
     filename=TRIM(filename)//TRIM(cpuchar)
     nbfiles=nbfiles+1
     filelist(nbfiles)=trim(filename)
     open(unit=ilun,file=trim(scratchdir)//trim(filename),form="unformatted",status="replace",action="write",iostat=ierr)
     if(ierr/=0)then
        print *,'Error: open file failed in backup_part_recv'
        call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
     end if

     call MPI_RECV(xdp,2,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)

     ! Write header
     write(ilun)ncpu
     write(ilun)ndim
     write(ilun)ii(IRandNumSize+1)  !npart
     npart=ii(IRandNumSize+1)
     write(ilun)ii(1:IRandNumSize) !localseed
     write(ilun)ii(IRandNumSize+2) !nstar_tot
     write(ilun)xdp(1) !mstar_tot
     write(ilun)xdp(2) !mstar_lost
     nsink=ii(IRandNumSize+3)
     write(ilun)nsink
     deallocate(ii,xdp)

     allocate(xdp(1:npart))
     ! Write position
     do idim=1,ndim
        call MPI_RECV(xdp,npart,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
        write(ilun)xdp
     end do
     ! Write velocity
     do idim=1,ndim
        call MPI_RECV(xdp,npart,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
        write(ilun)xdp
     end do
     ! Write mass
     call MPI_RECV(xdp,npart,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
     write(ilun)xdp
     ! Write identity
     allocate(ii(1:npart))
     call MPI_RECV(ii,npart,MPI_INTEGER,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
     write(ilun)ii
     deallocate(ii)
     ! Write level
     allocate(ll(1:npart))
     call MPI_RECV(ll,npart,MPI_INTEGER,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
     write(ilun)ll
     deallocate(ll)
     if(star.or.sink)then
        ! Write birth epoch
        call MPI_RECV(xdp,npart,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
       write(ilun)xdp
        if(metal)then
           ! Write metallicity
           call MPI_RECV(xdp,npart,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
           write(ilun)xdp
        endif
     end if
     deallocate(xdp)

    ! Write sink properties
    if(sink.and.nsink>0)then
      allocate(xdp(ndim*nsink))
      call MPI_RECV(xdp,nsink,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
      write(ilun)xdp(1:nsink)
      call MPI_RECV(xdp,nsink*ndim,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
      write(ilun)xdp
      deallocate(xdp)
    endif

     close(ilun)

     count = count+1
  end do

  deallocate(list_recv)

end subroutine backup_part_recv
#endif

#endif
