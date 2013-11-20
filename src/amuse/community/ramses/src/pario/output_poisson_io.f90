! Copyright (c) 2006-2011 IDRIS/CNRS
! Author: Philippe Wautelet (IDRIS/CNRS), wautelet@idris.fr
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.

#ifdef IOPROC

#ifndef WITHOUTMPI

subroutine backup_poisson_send
  use amr_commons
  use poisson_commons
  use io_parameters
  use timer
  use mpi
  implicit none

  integer,parameter::tag=TAG_BAK_POI
  integer::i,ierr
  integer::ivar,ncache,ind,ilevel,igrid,iskip,istart,ibound
  integer,allocatable,dimension(:)::ind_grid
  real(dp),allocatable,dimension(:)::xdp

  if(verbose)write(*,*)'Entering backup_poisson_send'

  call start_timer()

  call MPI_SEND(nboundary,1,MPI_INTEGER,0,tag,MPI_COMM_IOGROUP,ierr)

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
              do ivar=1,ndim
                 do i=1,ncache
                    xdp(i)=f(ind_grid(i)+iskip,ivar)
                 end do
                 call MPI_SEND(xdp,ncache,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_IOGROUP,ierr)
              end do
           end do
           deallocate(ind_grid, xdp)
        end if
     end do
  end do

  call stop_timer('poisson I/O processes backup',writing=.true.)
end subroutine backup_poisson_send

subroutine backup_poisson_recv
  use amr_commons
  use io_commons
  use pm_commons
  use mpi
  implicit none

  integer,parameter::tag=TAG_BAK_POI
  integer,dimension(MPI_STATUS_SIZE)::status
  logical,allocatable,dimension(:)::list_recv
  integer::ibound,ilevel,ilun,ind,ivar,ncache,count,src,ierr
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=MAXLINE)::filename
  character(LEN=5)::cpuchar,iochar,nchar

  if(verbose)write(*,*)'Entering backup_poisson_recv'

  allocate(list_recv(ncpu_iogroup-1))
  list_recv(:)=.false.

  count=0
  do while(count<ncpu_iogroup-1)
     ! Select a source
     call MPI_RECV(nboundary,1,MPI_INTEGER,MPI_ANY_SOURCE,tag,MPI_COMM_IOGROUP,status,ierr)

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
     filename=TRIM(filename)//'grav_'//TRIM(nchar)//'.out'
     filename=TRIM(filename)//TRIM(cpuchar)
     nbfiles=nbfiles+1
     filelist(nbfiles)=trim(filename)
     open(unit=ilun,file=trim(scratchdir)//trim(filename),form="unformatted",status="replace",action="write",iostat=ierr)
     if(ierr/=0)then
        print *,'Error: open file failed in backup_part_recv'
        call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
     end if

    write(ilun)ncpu
    write(ilun)ndim
    write(ilun)nlevelmax
    write(ilun)nboundary
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
             do ivar=1,ndim
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

end subroutine backup_poisson_recv
#endif

#endif
