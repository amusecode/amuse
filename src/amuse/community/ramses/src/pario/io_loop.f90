! Copyright (c) 2006-2011 IDRIS/CNRS
! Author: Philippe Wautelet (IDRIS/CNRS), wautelet@idris.fr
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.

#ifdef IOPROC
subroutine io_loop
  use amr_commons
  use io_commons
  use mpi
  implicit none

  integer::i,ierr,ilun=ILUN_IO,size,type_oper
  real(kind=8)::freespace_perm,freespace_scratch,min_freespace_scratch_glob
  real(kind=8)::speed_comp2scratch,speed_comp2scratch_min,speed_scratch2perm,speed_scratch2perm_min
  real(kind=8),dimension(3)::buf1,buf2
  integer(kind=8)::KB,MB,GB,TB
  integer,parameter::tag=100
  real(kind=8)::total_size
  logical::exist,flag,error
  character(len=5)::iochar
  character(len=MAXLINE)::src,dest,hostname,filename,status
  real(kind=8)::t1,t2,t3,t4
  integer,dimension(MPI_STATUS_SIZE)::stat

  KB=1024.d0
  MB=KB*1024.d0
  GB=MB*1024.d0
  TB=GB*1024.d0

  error=.false.

  allocate(numblio(1:MAX(ncpu_iogroup-1,ncpu),1:nlevelmax,1:ncpu_iogroup-1))
  allocate(numbb(1:MAXBOUND,1:nlevelmax))
  allocate(filelist(MAXFILES*(ncpu_iogroup-1)))

  if(index(scratchdir,'/',back=.true.)/=len_trim(scratchdir)) scratchdir=trim(scratchdir)//'/'
  if(index(permdir   ,'/',back=.true.)/=len_trim(permdir))    permdir   =trim(permdir)//'/'

  hostname=' ' ! Necessary to have a valid Fortran string
  call gethname(hostname)

  call getfreespace(trim(scratchdir)//char(0),freespace_scratch)
  call getfreespace(trim(permdir)   //char(0),freespace_perm)
  call MPI_ALLREDUCE(freespace_scratch,min_freespace_scratch_glob,1,MPI_INTEGER8,MPI_MIN,MPI_COMM_COMP,ierr)

  call CreateDirs(trim(scratchdir)//char(0),trim(permdir)//char(0),&
                 my_iogroup,ncpu_iogroup,iogroup2comp(1))

  ! Receive number of the current output and backup
  if(nrestart>0)then
     call MPI_RECV(count_out,1,MPI_INTEGER,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
     call MPI_RECV(count_bak,1,MPI_INTEGER,1,tag,MPI_COMM_IOGROUP,MPI_STATUS_IGNORE,ierr)
     count_out=count_out-1
     count_bak=count_bak-1
  end if

  ! Create the IO stat file
  call title(my_iogroup,iochar)
  filename=trim(permdir)//'ionode_'//iochar//'/IOstats_'//iochar//'.txt'
  if (nrestart>0) then
     inquire(file=filename,exist=exist)
     if (exist) then
        status='old'
     else
        status='new'
     end if
     open(unit=ilun,file=filename,status=status,form='formatted',position='append',action='write',iostat=ierr)
     if(ierr/=0)then
        print *,'Error: open file failed in io_loop'
     end if
  else
     open(unit=ilun,file=filename,status='replace',form='formatted',position='rewind',action='write',iostat=ierr)
     if(ierr/=0)then
        print *,'Error: open file failed in io_loop'
     end if
  endif

  ! Write basic infos to the IO stat file
  if(nrestart>0)then
     write(ilun,*) '------------------------------------------------------------'
     write(ilun,*) 'RESTART #',nrestart
  end if
  write(ilun,*) '------------------------------------------------------------'
  write(ilun,*) 'Ionode: ',my_iogroup
  write(ilun,*) 'Computational processes managed:'
  write(ilun,'(t2,9I6)') iogroup2comp(1:ncpu_iogroup-1)
  write(ilun,*) 'Hostname: ',trim(hostname)
  write(ilun,*) 'Scratch dir: ',trim(scratchdir)
  write(ilun,*) 'Perm dir:    ',trim(permdir)
  write(ilun,*) 'Free space:'

  if(min_freespace_scratch_glob<KB)then
     write(ilun,'(A,F0.0,A)') '  Smallest scratch dir: ',min_freespace_scratch_glob,'B'
  else if(min_freespace_scratch_glob<MB)then
     write(ilun,'(A,F7.2,A)') '  Smallest scratch dir: ',min_freespace_scratch_glob/KB,'KB'
  else if(min_freespace_scratch_glob<GB)then
     write(ilun,'(A,F7.2,A)') '  Smallest scratch dir: ',min_freespace_scratch_glob/MB,'MB'
  else if(min_freespace_scratch_glob<TB)then
     write(ilun,'(A,F7.2,A)') '  Smallest scratch dir: ',min_freespace_scratch_glob/GB,'GB'
  else
     write(ilun,'(A,F7.2,A)') '  Smallest scratch dir: ',min_freespace_scratch_glob/TB,'TB'
  endif
  if(freespace_scratch<KB)then
     write(ilun,'(A,F0.0,A)') '  Scratch dir: ',freespace_scratch,'B'
  else if(freespace_scratch<MB)then
     write(ilun,'(A,F7.2,A)') '  Scratch dir: ',freespace_scratch/KB,'KB'
  else if(freespace_scratch<GB)then
     write(ilun,'(A,F7.2,A)') '  Scratch dir: ',freespace_scratch/MB,'MB'
  else if(freespace_scratch<TB)then
     write(ilun,'(A,F7.2,A)') '  Scratch dir: ',freespace_scratch/GB,'GB'
  else
     write(ilun,'(A,F7.2,A)') '  Scratch dir: ',freespace_scratch/TB,'TB'
  endif
  if(freespace_perm<KB)then
     write(ilun,'(A,F7.0,A)') '  Perm dir:    ',freespace_perm,'B'
  else if(freespace_perm<MB)then
     write(ilun,'(A,F7.2,A)') '  Perm dir:    ',freespace_perm/KB,'KB'
  else if(freespace_perm<GB)then
     write(ilun,'(A,F7.2,A)') '  Perm dir:    ',freespace_perm/MB,'MB'
  else if(freespace_perm<TB)then
     write(ilun,'(A,F7.2,A)') '  Perm dir:    ',freespace_perm/GB,'GB'
  else
     write(ilun,'(A,F7.2,A)') '  Perm dir:    ',freespace_perm/TB,'TB'
  endif
  write(ilun,*) '------------------------------------------------------------'
  ! Flush the stats (not very clean but it works...)
  close(ilun)
  open(unit=ilun,file=filename,status='old',form='formatted',position='append',action='write',iostat=ierr)
  if(ierr/=0)then
     print *,'Error: re-open file failed in io_loop'
  end if

  do
     nbfiles=0

     t1 = MPI_WTIME()
     ! The following busywait loop is necessary to prevent timeouts
     ! Other solution: play with MP_TIMEOUT (on IBM) and just do a MPI_RECV
     do
        ! Do not use mpi_status_ignore because there is a Blue Gene/P bug
        ! that will crash the application if used (20090717)
        call MPI_IPROBE(1,tag,MPI_COMM_IOGROUP,flag,stat,ierr)
        if(flag.eqv..true.)exit
     end do
     t2 = MPI_WTIME()

     ! 0=output, 1=backup, 2=end of computation, 3=output+backup
     call MPI_RECV(type_oper,1,MPI_INTEGER,1,tag,MPI_COMM_IOGROUP,status,ierr)

     select case(type_oper)
       case(0)
          call local_output()

       case(1)
          call local_backup()

       case(2)
          close(ilun)
          call clean_stop

       case(3)
          call local_output()
          call local_backup()
     end select

     call getfreespace(trim(scratchdir)//char(0),freespace_scratch)

     t3 = MPI_WTIME()
     total_size=0
     if(scratchdir/=permdir)then
        do i=1,nbfiles
           src= trim(scratchdir)//trim(filelist(i))//char(0)
           dest=trim(permdir)   //trim(filelist(i))//char(0)
           call TransferFile(src,dest,size)
           if(size==-1)then
              error=.true.
              write(ilun,*) 'ERROR: File ',src,' has not been transfered to permdir'
          else
              total_size=total_size+size
           end if
        end do
     end if
     t4 = MPI_WTIME()

     speed_comp2scratch = total_size/(MB*(t3-t2))
     speed_scratch2perm = total_size/(MB*(t4-t3))

     call getfreespace(trim(permdir)//char(0),freespace_perm)
     buf1(1)=freespace_scratch
     buf1(2)=speed_comp2scratch
     buf1(3)=speed_scratch2perm
     call MPI_ALLREDUCE(buf1,buf2,3,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_COMP,ierr)
     min_freespace_scratch_glob=buf2(1)
     speed_comp2scratch_min=buf2(2)
     speed_scratch2perm_min=buf2(3)

     ! Write stats
     select case(type_oper)
        case(0)
           write(ilun,*) 'Output #',count_out
        case(1)
           write(ilun,*) 'Backup #',count_bak
        case(3)
           write(ilun,*) 'Output #',count_out,' and Backup #',count_bak
     end select
     if(error)then
        write(ilun,*) 'ERROR(S) DURING TRANSFER (see error output for more info)'
        write(ilun,*) 'File(s) have been kept on scratchdir'
     end if
     write(ilun,*) 'Number of files: ',nbfiles
     if(total_size<KB)then
        write(ilun,'(A,F7.0,A)') ' Data volume: ',total_size,'B'
     else if(total_size<MB)then
       write(ilun,'(A,F7.2,A)') ' Data volume: ',total_size/KB,'KB'
     else if(total_size<GB)then
        write(ilun,'(A,F7.2,A)') ' Data volume: ',total_size/MB,'MB'
     else if(total_size<TB)then
        write(ilun,'(A,F7.2,A)') ' Data volume: ',total_size/GB,'GB'
     else
        write(ilun,'(A,F7.2,A)') ' Data volume: ',total_size/TB,'TB'
     end if
     write(ilun,'(A,F8.2,A)') ' Wait for data:             ',t2-t1,'s'
     write(ilun,'(A,F8.2,A,F0.2,A,F0.2,A)') &
          ' Receives and local writes: ',t3-t2,'s (',speed_comp2scratch,'MB/s, slowest ',speed_comp2scratch_min,'MB/s)'
     write(ilun,'(A,F8.2,A,F0.2,A,F0.2,A)') &
          ' Scratch to perm transfer:  ',t4-t3,'s (',speed_scratch2perm,'MB/s, slowest ',speed_scratch2perm_min,'MB/s)'
     write(ilun,*) 'Free space (after writes, before transfers for scratch):'
     if(min_freespace_scratch_glob<KB)then
        write(ilun,'(A,F0.0,A)') '  Smallest scratch dir: ',min_freespace_scratch_glob,'B'
     else if(min_freespace_scratch_glob<MB)then
        write(ilun,'(A,F7.2,A)') '  Smallest scratch dir: ',min_freespace_scratch_glob/KB,'KB'
     else if(min_freespace_scratch_glob<GB)then
        write(ilun,'(A,F7.2,A)') '  Smallest scratch dir: ',min_freespace_scratch_glob/MB,'MB'
     else if(min_freespace_scratch_glob<TB)then
        write(ilun,'(A,F7.2,A)') '  Smallest scratch dir: ',min_freespace_scratch_glob/GB,'GB'
     else
        write(ilun,'(A,F7.2,A)') '  Smallest scratch dir: ',min_freespace_scratch_glob/TB,'TB'
     endif
     if(freespace_scratch<KB)then
        write(ilun,'(A,F0.0,A)') '  Scratch dir: ',freespace_scratch,'B'
     else if(freespace_scratch<MB)then
        write(ilun,'(A,F7.2,A)') '  Scratch dir: ',freespace_scratch/KB,'KB'
     else if(freespace_scratch<GB)then
        write(ilun,'(A,F7.2,A)') '  Scratch dir: ',freespace_scratch/MB,'MB'
     else if(freespace_scratch<TB)then
        write(ilun,'(A,F7.2,A)') '  Scratch dir: ',freespace_scratch/GB,'GB'
     else
        write(ilun,'(A,F7.2,A)') '  Scratch dir: ',freespace_scratch/TB,'TB'
     endif
     if(freespace_perm<KB)then
        write(ilun,'(A,F7.0,A)') '  Perm dir:    ',freespace_perm,'B'
     else if(freespace_perm<MB)then
        write(ilun,'(A,F7.2,A)') '  Perm dir:    ',freespace_perm/KB,'KB'
     else if(freespace_perm<GB)then
        write(ilun,'(A,F7.2,A)') '  Perm dir:    ',freespace_perm/MB,'MB'
     else if(freespace_perm<TB)then
        write(ilun,'(A,F7.2,A)') '  Perm dir:    ',freespace_perm/GB,'GB'
     else
        write(ilun,'(A,F7.2,A)') '  Perm dir:    ',freespace_perm/TB,'TB'
     endif
     write(ilun,*) '------------------------------------------------------------'
     ! Flush the stats (not very clean but it works...)
     close(ilun)
     open(unit=ilun,file=filename,status='old',form='formatted',position='append',action='write',iostat=ierr)
     if(ierr/=0)then
        print *,'Error: re-open file failed in io_loop'
     end if     
  end do

end subroutine io_loop

subroutine local_output
  use amr_commons
  use io_commons
  use mpi
  implicit none

  integer ierr

  print *,'Error: local_output does not exist anymore!'
  call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
end subroutine local_output

subroutine local_backup
  use amr_commons
  use io_commons
  implicit none

  count_bak = count_bak+1

  if(myid_world==1) call output_info_recv
  if(cooling.and.myid_world==1) call output_cool_recv
  call backup_amr_recv
  if(hydro) call backup_hydro_recv
  if(pic)   call backup_part_recv
  if(poisson) call backup_poisson_recv(.true.)

end subroutine local_backup
#endif
