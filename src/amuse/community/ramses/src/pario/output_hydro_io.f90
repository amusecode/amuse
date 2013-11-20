! Copyright (c) 2006-2011 IDRIS/CNRS
! Author: Philippe Wautelet (IDRIS/CNRS), wautelet@idris.fr
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.

#ifdef IOPROC

#ifndef WITHOUTMPI

subroutine backup_hydro_send
  use amr_commons
  use hydro_commons
  use io_parameters
  use mpi
  use timer
  implicit none

  integer,parameter::tag=TAG_BAK_HYD
  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,istart,ibound,ierr
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp

  if(verbose)write(*,*)'Entering backup_hydro_send'

  call start_timer()

  call MPI_SEND(nboundary,1,MPI_INTEGER,0,tag,MPI_COMM_IOGROUP,ierr)
  call MPI_SEND(gamma,1,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)

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
           allocate(ind_grid(1:ncache),xdp(1:ncache))
           ! Loop over level grids
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do ivar=1,nvar
                 if(ivar==1)then ! Write density
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,1)
                    end do
                 else if(ivar>=2.and.ivar<=ndim+1)then ! Write velocity field
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,ivar)/uold(ind_grid(i)+iskip,1)
                    end do
                 else if(ivar==ndim+2)then ! Write pressure
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,ndim+2)
                       xdp(i)=xdp(i)-0.5d0*uold(ind_grid(i)+iskip,2)**2/uold(ind_grid(i)+iskip,1)
#if NDIM>1
                       xdp(i)=xdp(i)-0.5d0*uold(ind_grid(i)+iskip,3)**2/uold(ind_grid(i)+iskip,1)
#endif
#if NDIM>2
                       xdp(i)=xdp(i)-0.5d0*uold(ind_grid(i)+iskip,4)**2/uold(ind_grid(i)+iskip,1)
#endif
                       xdp(i)=(gamma-1d0)*xdp(i)
                    end do
                 else ! Write passive scalars if any
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,ivar)/uold(ind_grid(i)+iskip,1)
                    end do
                 endif
                 call MPI_SEND(xdp,ncache,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
              end do
           end do
           deallocate(ind_grid, xdp)
        end if
     end do
  end do

  call stop_timer('hydro I/O processes backup',writing=.true.)

end subroutine backup_hydro_send

subroutine backup_hydro_recv
  use amr_commons
  use hydro_commons
  use io_commons
  use mpi
  implicit none

  integer::count,ierr,src
  integer,parameter::tag=TAG_BAK_HYD
  integer,dimension(MPI_STATUS_SIZE)::status
  integer::ivar,ncache,ind,ilevel,ilun,ibound
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::cpuchar,iochar,nchar
  character(LEN=MAXLINE)::filename
  logical,allocatable,dimension(:)::list_recv

  if(verbose)write(*,*)'Entering backup_hydro_recv'
  allocate(list_recv(ncpu_iogroup-1))
  list_recv(:)=.false.

  count=0
  do while(count<ncpu_iogroup-1)
     ! Select a source
     ! mpi_probe is used because it is too early to receive
     ! the message (the filename is not yet determined by example)
     call MPI_PROBE(MPI_ANY_SOURCE,tag,MPI_COMM_IOGROUP,status,ierr)

     src=status(MPI_SOURCE)
     if(list_recv(src).EQV..false.)then
        list_recv(src)=.true.
     else
        print *,'Error: unexpected message received by ',myid_world
        call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
     end if

     ! Generate filename
     ilun=ncpu_world+myid_world+10
     call title(count_bak,nchar)
     call title(iogroup2comp(src),cpuchar)
     call title(myid_io,iochar)
     filename='ionode_'//TRIM(iochar)//'/process_'//TRIM(cpuchar)//'/'
     filename=TRIM(filename)//'hydro_'//TRIM(nchar)//'.out'
     filename=TRIM(filename)//TRIM(cpuchar)
     nbfiles=nbfiles+1
     filelist(nbfiles)=trim(filename)
     open(unit=ilun,file=trim(scratchdir)//trim(filename),status="replace",form="unformatted",action="write",iostat=ierr)
     if(ierr/=0)then
        print *,'Error: open file failed in backup_hydro_recv'
        call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
     end if

     write(ilun)ncpu
     write(ilun)nvar
     write(ilun)ndim
     write(ilun)nlevelmax
     call MPI_RECV(nboundary,1,MPI_INTEGER,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
     write(ilun)nboundary
     call MPI_RECV(gamma,1,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
     write(ilun)gamma

    do ilevel=1,nlevelmax
      do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
          ncache=numblio(ibound,ilevel,src)
        else
          ncache=numbb(ibound-ncpu,ilevel)
        end if
        write(ilun)ilevel
        write(ilun)ncache
        if(ncache>0)then
          allocate(xdp(1:ncache))
          ! Loop over cells
          do ind=1,twotondim
            do ivar=1,nvar
              call MPI_RECV(xdp,ncache,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
              write(ilun)xdp
            end do
          end do
          deallocate(xdp)
        end if
      end do
    end do

     close(ilun)

     count = count+1
  end do

  deallocate(list_recv)

end subroutine backup_hydro_recv
#endif

#endif
