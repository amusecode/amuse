! Copyright (c) 2007-2011 IDRIS/CNRS
! Author: Philippe Wautelet (IDRIS/CNRS), wautelet@idris.fr
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.

#ifndef WITHOUTMPI

subroutine comp_offsets_sizes_by_level()
  use amr_commons,only:mpi_comm_comp,myid,ncpu,nboundary,nlevelmax,numbb,numbl
  use io_commons
  use mpi
  implicit none

  ! The multiplier is used to increaze the size of ncache_max by this factor to determine the buffer size

  integer::i,j,status,output
  integer(kind=mpi_offset_kind)::size_loc,offset,sz
  integer(kind=mpi_offset_kind),dimension(nlevelmax,ncpu)::sizes_loc_by_level

  if(multiplier<1.d0) print*,'Warning: multiplier should be greater or equal to 1.0d0'

  !Compute indices and variable sizes
  size_loc=0
  size_loc_by_level(:)=0
  ncache_max=0
  do i=1,nlevelmax
     do j=1,ncpu
        size_loc_by_level(i) = size_loc_by_level(i)+numbl(j,i)
     end do
     do j=1,nboundary
        size_loc_by_level(i) = size_loc_by_level(i)+numbb(j,i)
     end do
     size_loc = size_loc + size_loc_by_level(i)
     if(size_loc_by_level(i)>ncache_max) ncache_max = size_loc_by_level(i)
  end do
  call MPI_ALLREDUCE(size_loc,size_max,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_COMP,status)
  call MPI_ALLGATHER(size_loc_by_level, nlevelmax,MPI_INTEGER8,&
                     sizes_loc_by_level,nlevelmax,MPI_INTEGER8,MPI_COMM_COMP,status)

  ! Compute offset of each level for each process
  offset=0
  do i=1,nlevelmax
     level_offsets(i)=offset
     do j=1,ncpu
        if(j<myid) level_offsets(i)=level_offsets(i)+sizes_loc_by_level(i,j)
        offset=offset+sizes_loc_by_level(i,j)
     end do
  end do

  ! Group the different outputs to fill the buffer_size
  ! This allows to reduce the number of I/O and to increaze
  ! the average size of the I/O accesses
  sz=0
  output=1
  outputs=.false.
  buffer_size=max(1.0d0,multiplier)*ncache_max
  do i=1,nlevelmax-1
    if(sz+size_loc_by_level(i+1) > buffer_size)then
      sz=0
      outputs(i)=.true.
    end if
    !TODO: do not do a allreduce but compute from the sizes_loc_by_level
    call MPI_ALLREDUCE(MPI_IN_PLACE,outputs(i),1,MPI_LOGICAL,MPI_LOR,MPI_COMM_COMP,status)
    sz=sz+size_loc_by_level(i+1)
  end do
  outputs(nlevelmax)=.true.

end subroutine comp_offsets_sizes_by_level

#endif

#ifndef WITHOUTMPI
!Version with MPI
module timer
  use amr_commons,only:myid,MPI_COMM_COMP
  use mpi

  implicit none

  integer,private     ::ierr
  real(kind=8),private::te,tb,ti
  real(kind=8),private::tmin,tmax

  contains

  subroutine start_timer
    call MPI_BARRIER(MPI_COMM_COMP,ierr)
    tb = MPI_WTIME()
    ti = tb
  end subroutine start_timer

  subroutine stop_timer(filetype,writing)
    character(len=*),intent(in)::filetype
    logical,intent(in)::writing

    character(len=28)::line

    te = MPI_WTIME()

    call MPI_REDUCE(te-tb,tmin,1,MPI_DOUBlE_PRECISION,MPI_MIN,0,MPI_COMM_COMP,ierr)
    call MPI_REDUCE(te-tb,tmax,1,MPI_DOUBlE_PRECISION,MPI_MAX,0,MPI_COMM_COMP,ierr)

    if(myid==1) then
      line=trim(filetype)//':'
      if(writing) then
        write(*,'(" Time to write ",A28," max=",F8.3,"s min=",F8.3,"s")')line,tmax,tmin
      else
        write(*,'(" Time to read ",A28," max=",F8.3,"s min=",F8.3,"s")')line,tmax,tmin
      end if
    end if
  end subroutine stop_timer

  subroutine interm_timer(text)
    character(len=*),intent(in)::text

    te = MPI_WTIME()

    if(myid==1) print *,'  ',text,': ',te-ti,'s'
    ti=te
  end subroutine interm_timer
end module
#else
!Version without MPI
module timer
  implicit none

  integer,private :: ticks_beg, ticks_end, ticks_interm
  integer,private :: ticks, ticks_max, ticks_by_sec
  real,private    :: elapsed

  contains

  subroutine start_timer
    call SYSTEM_CLOCK(COUNT_RATE=ticks_by_sec,COUNT_MAX=ticks_max)

    call SYSTEM_CLOCK(COUNT=ticks_beg)
  end subroutine start_timer

  subroutine stop_timer(filetype,writing)
    character(len=*),intent(in)::filetype
    logical,intent(in)::writing

    character(len=28)::line

    call SYSTEM_CLOCK(COUNT=ticks_end)
    ticks = ticks_end-ticks_beg
    if(ticks_end<ticks_beg) ticks = ticks + ticks_max
    elapsed = real(ticks)/ticks_by_sec

    line=trim(filetype)//':'
    if(writing) then
      write(*,'(" Time to write ",A28,F8.3,"s")')line,elapsed
    else
      write(*,'(" Time to read ",A28,F8.3,"s")')line,elapsed
    end if
  end subroutine stop_timer

  subroutine interm_timer(text)
   character(len=*),intent(in)::text

    call SYSTEM_CLOCK(COUNT=ticks_end)
    ticks = ticks_end-ticks_interm
    if(ticks_end<ticks_interm) ticks = ticks + ticks_max
    elapsed = real(ticks)/ticks_by_sec

    print *,'  ',text,': ',elapsed,'s'

    ticks_interm = ticks_end
  end subroutine interm_timer
end module
#endif
