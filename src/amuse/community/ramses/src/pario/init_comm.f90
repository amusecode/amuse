! Copyright (c) 2006-2011 IDRIS/CNRS
! Author: Philippe Wautelet (IDRIS/CNRS), wautelet@idris.fr
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.

#ifndef WITHOUTMPI

! This subroutine is executed by every process
subroutine init_comm
  use amr_commons
  use mpi
  implicit none

  integer::color_iogroup,color_comp,i,ierr
  integer::nbgbl,szbgbl,szsmbl

  ! Check if ncpus are coherent
  if(ncpu_io==0)then
     if(ncpu==0)then
        ncpu = ncpu_world
     else
        ncpu_io = ncpu_world-ncpu
     end if
  else
     if(ncpu==0)then
        ncpu = ncpu_world-ncpu_io
     end if
  end if

#ifndef IOPROC
  if(ncpu_io/=0) then
     if(myid_world==1)write(*,*)'Error:'
     if(myid_world==1)write(*,*)'Compiled without specialized I/O processes support'
     call clean_stop
  end if
#endif
  if(ncpu_io/=0.AND.(fileformat_in/='std'.OR.fileformat_out/='std')) then
    if(myid_world==1)write(*,*)'Error:'
    if(myid_world==1)write(*,*)'Specialized I/O processes only with std fileformats'
    call clean_stop
  end if

  if(ncpu+ncpu_io/=ncpu_world.OR.ncpu<0.OR.ncpu_io<0)then
     if(myid_world==1)write(*,*)'Error in the namelist:'
     if(myid_world==1)write(*,*)'ncpu/ncpu_io distribution is incoherent'
     call clean_stop
  end if

  if(ncpu_io>ncpu)then
     if(myid_world==1)write(*,*)'Error in the namelist:'
     if(myid_world==1)write(*,*)'ncpu should be at least ncpu_io'
     call clean_stop
  end if

  ! Distribute CPUs
  if(ncpu_io/=0)then
     if(mod(ncpu_world,ncpu_io)==0)then
        nbgbl  = ncpu_io
        szbgbl = ncpu_world/ncpu_io
        szsmbl = 0
     else
        nbgbl  = mod(ncpu_world,ncpu_io)
        szbgbl = ncpu_world/ncpu_io+1
        szsmbl = ncpu_world/ncpu_io
     end if
     if(myid_world<=nbgbl*szbgbl)then
        color_iogroup = (myid_world-1)/szbgbl
        if(mod(myid_world,szbgbl)==1)then
           color_comp = 0
        else
           color_comp = 1
        end if
     else
        color_iogroup = (myid_world-1-nbgbl*szbgbl)/szsmbl+nbgbl
        if(mod(myid_world-nbgbl*szbgbl,szsmbl)==1)then
           color_comp = 0
        else
           color_comp = 1
        end if
     end if
  else
     ! The consequence of this is that MPI_COMM_COMP=MPI_COMM_IOGROUP=MPI_COMM_WORLD
     color_iogroup = 0
     color_comp    = 1
  end if

  ! Initialize communicators
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,color_iogroup,myid_world,MPI_COMM_IOGROUP,ierr)
  call MPI_COMM_RANK(MPI_COMM_IOGROUP,myid_iogroup,ierr)
  call MPI_COMM_SIZE(MPI_COMM_IOGROUP,ncpu_iogroup,ierr)
  if(ncpu_io>0)then
     myid_iogroup = myid_iogroup+1
  else
     ! Necessary to prevent problems later (a wrong branch can be taken by a process)
     myid_iogroup=MPI_UNDEFINED
  end if

  call MPI_COMM_SPLIT(MPI_COMM_WORLD,color_comp,myid_world,MPI_COMM_COMP,ierr)
  call MPI_COMM_RANK(MPI_COMM_COMP,myid,ierr)
  if(color_comp==1)then
     myid    = myid+1
     myid_io = MPI_UNDEFINED
  else
     myid_io = myid+1
     myid    = MPI_UNDEFINED
     allocate(iogroup2world(0:ncpu_iogroup-1),iogroup2comp(0:ncpu_iogroup-1))
     do i=0,ncpu_iogroup-1
        iogroup2world(i)=myid_world-1+i
     end do
     iogroup2comp(0)=MPI_UNDEFINED
     do i=1,ncpu_iogroup-1
        iogroup2comp(i)=myid_world-myid_io+i
     end do
  end if

  ncoarse = nx*ny*nz

  my_iogroup=color_iogroup+1

  ncpu = ncpu
end subroutine init_comm

#endif
