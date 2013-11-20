! Copyright (c) 2006-2011 IDRIS/CNRS
! Author: Philippe Wautelet (IDRIS/CNRS), wautelet@idris.fr
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.

#ifdef IOPROC

#ifndef WITHOUTMPI

!subroutine output_cool_send => see in cooling_module

subroutine output_cool_recv
  use amr_commons
  use cooling_module,only:if_species_abundances
  use io_commons
  use mpi
  implicit none
  character(LEN=MAXLINE)::filename
  character(LEN=5)::nchar,iochar
  integer nn,nt,ierr,ilun
  integer,parameter::tag=TAG_OUT_COO
  real(kind=8),dimension(:),allocatable::data

  ! Generate filename
  ilun=myid_world+10
  call title(count_bak,nchar)
  call title(myid_io,iochar)
  filename='ionode_'//TRIM(iochar)//'/process_00001/'
  filename=TRIM(filename)//'cooling_'//TRIM(nchar)//'.out'
  nbfiles=nbfiles+1
  filelist(nbfiles)=trim(filename)
  open(unit=ilun,file=trim(scratchdir)//trim(filename),status="replace",form="unformatted",action="write",iostat=ierr)
  if(ierr/=0)then
     print *,'Error: open file failed in output_cool_recv'
     call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
  end if

  call MPI_RECV(nn,1,MPI_INTEGER,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
  call MPI_RECV(nt,1,MPI_INTEGER,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
  write(ilun)nn,nt
  allocate(data(nn*nt*6))

  call MPI_RECV(data,nn,   MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
  write(ilun)data(1:nn)
  call MPI_RECV(data,nt,   MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
  write(ilun)data(1:nt)
  call MPI_RECV(data,nn*nt,MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
  write(ilun)data(1:nn*nt)
  call MPI_RECV(data,nn*nt,MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
  write(ilun)data(1:nn*nt)
  call MPI_RECV(data,nn*nt,MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
  write(ilun)data(1:nn*nt)
  call MPI_RECV(data,nn*nt,MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
  write(ilun)data(1:nn*nt)
  call MPI_RECV(data,nn*nt,MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
  write(ilun)data(1:nn*nt)
  call MPI_RECV(data,nn*nt,MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
  write(ilun)data(1:nn*nt)
  call MPI_RECV(data,nn*nt,MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
  write(ilun)data(1:nn*nt)
  call MPI_RECV(data,nn*nt,MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
  write(ilun)data(1:nn*nt)
  call MPI_RECV(data,nn*nt,MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
  write(ilun)data(1:nn*nt)
  call MPI_RECV(data,nn*nt,MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
  write(ilun)data(1:nn*nt)
  call MPI_RECV(data,nn*nt,MPI_DOUBLE_PRECISION,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
  write(ilun)data(1:nn*nt)
  if (if_species_abundances) then
     call MPI_RECV(data,nn*nt*6,MPI_DOUBLE_PRECISION,&
                   1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
     write(ilun)data(1:6*nn*nt)
  end if

  close(ilun)

  deallocate(data)

end subroutine output_cool_recv

#endif

#endif
