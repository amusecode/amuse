!> \file source_input.F90

!> \brief Handles readin of sourcefiles
!<

module source_input_mod
use myf03_mod
use particle_system_mod, only: particle_system_type
use particle_system_mod, only: source_type
use global_mod, only: psys, GV, PLAN
implicit none

!> source header type 
!-------------------------
type source_header_type 
   integer(i8b) :: NsrcSnap  !< number of sources in snapshot
   integer(i8b) :: NsrcFile  !< number of sources in file
   integer(i8b) :: NumFiles  !< number of files per snapshot
   real(r8b) :: TotalRays    !< sum of all rays to be cast this snap
   real(r8b) :: Lunit        !< Luminosities are photons/s * Lunit
end type source_header_type

contains

!> reads in a source header, prints the header to the screen if verbose
!========================================================================
subroutine read_source_header(snapfile, shead, lun, closefile)

  character(*), intent(in) :: snapfile  !< file containing source header
  type(source_header_type), intent(inout) :: shead  !< source header to read
  integer(i4b), intent(out) :: lun  !< output lun assigned to snapfile
  logical, optional, intent(in) :: closefile      !< close file when done?

  logical :: closef

  if (.not. present(closefile) ) then
     closef = .true.
  else
     closef = closefile
  end if

  call open_formatted_file_r(snapfile,lun)

  read(lun,*) shead%NsrcSnap
  read(lun,*) shead%NsrcFile
  read(lun,*) shead%NumFiles
  read(lun,*) shead%TotalRays
  read(lun,*) shead%Lunit
  read(lun,*)
  read(lun,*)
  if (closef) close(lun)
  
end subroutine read_source_header


subroutine get_planning_data_sources()

  character(clen), parameter :: myname = 'get_planning_data_sources'
  logical, parameter :: crash = .true.
  integer, parameter :: verb = 2

  integer(i4b) :: loglun
  type(source_header_type) :: shead
  character(clen) :: snapfile ! snapshot file name
  character(clen) :: logfile

  integer(i4b) :: iSnap, fSnap    ! initial and final snapshot numbers
  integer(i4b) :: sfiles          ! files/snap for sources   
  integer(i4b) :: i,j,lun


  ! these global variables are read from the config file
  !======================================================
  iSnap = GV%StartSnapNum
  fSnap = GV%EndSnapNum
    
  sfiles = GV%SourceFilesPerSnap


  ! open up the planning data log file
  !======================================================
  logfile = trim(GV%OutputDir) // "/" // "source_headers.log"
  call open_formatted_file_w(logfile,loglun)

  ! read all source headers and write to log file
  !===================================================
  write(loglun,*) 
  write(loglun,'(A)') "reading all source header(s) ... "
  do i = iSnap,fSnap
     do j = 1,sfiles
        call form_snapshot_file_name(GV%SourcePath,GV%SourceFileBase,i,j,snapfile)
        write(loglun,'(I3,"  ",A)') i,trim(snapfile)
        call mywrite('   srcfile = '//trim(snapfile) , verb)

        call read_source_header(snapfile,shead,lun,closefile=.true.)
        
        PLAN%snap(i)%RaysFromSrcHeader = shead%TotalRays
        GV%Lunit = shead%Lunit
        
        write(loglun,'(A,"  ",ES15.5)') '  Total Rays:', shead%TotalRays
        write(loglun,'(A,"  ",ES15.5)') '  L Unit:', shead%Lunit
        
     end do
  end do

  ! close headers log file
  !========================
  close(loglun)

  
end subroutine get_planning_data_sources


!> reads a source snapshot into arc array (src is allocated in this routine)
!=============================================================================
subroutine read_src_snapshot()

  character(clen), parameter :: myname="read_src_snapshot"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2
  character(clen) :: str, fmt

  character(clen) :: snapfile  
  type(source_header_type) :: shead
  integer(i4b) :: lun, err
  integer(i4b) :: i, N, Nall
  logical :: closefile
  integer(i4b) :: fn
  real(r4b) :: vel(3)
  real(r8b) :: MB

  fn=1
  call form_snapshot_file_name(GV%SourcePath,GV%SourceFileBase,GV%CurSnapNum,fn,snapfile)
      
  closefile = .false.
  call read_source_header(snapfile,shead,lun,closefile)

  Nall = shead%NsrcSnap        
  N    = shead%NsrcFile

  MB = GV%bytespersrc * real(Nall) / 2**20
  GV%MB = GV%MB + MB

  write(str,"(T4,A,F10.4,A,I10,A)") "allocating ", MB, " MB for ", Nall, " sources"
  call mywrite(str,verb)

  if (allocated(psys%src)) deallocate(psys%src)
  allocate (psys%src(Nall),stat=err)
  if (err /= 0) then
     call myerr("cant allocate sources in particle system", myname, crash)
  end if

  write(str,"(T4,A,A)") "reading source snapshot file: ", trim(snapfile) 
  call mywrite(str,verb)
     
  do i = 1,N 
     read(lun,*) psys%src(i)%pos(1), &
                 psys%src(i)%pos(2), &
                 psys%src(i)%pos(3), &
                 vel(1), &
                 vel(2), &
                 vel(3), &
                 psys%src(i)%L, &
                 psys%src(i)%SpcType, &
                 psys%src(i)%EmisPrf

#ifdef incVel
     psys%src(i)%vel = vel 
#endif

  end do

  close(lun)
 
  
end subroutine read_src_snapshot


!> Reorders the sources according to the array order. 
!! The array order is not preserved.
!========================================================
subroutine order_sources(srcs,order) 
  type(source_type), intent(inout) :: srcs(:)  !< sources to order
  integer(i8b), intent(inout) :: order(:)         !< desired order
  
  type(source_type) :: dumsrc
  integer(i8b) :: i,goal,nsrc
  
  if (size(srcs) /= size(order)) stop "size(srcs) /= size(order)"
  nsrc = size(srcs)

  do i = 1,nsrc
     dumsrc = srcs(i)
     goal = order(i)
     do while(goal < i)
        goal=order(goal)
        order(i)=goal
     end do
     srcs(i) = srcs(goal)
     srcs(goal)=dumsrc
  end do
  
end subroutine order_sources



!> Orders the sources from brightest to dimmest.  Also
!! sets the cumulative luminosity function entry for each source
!================================================================
subroutine order_sources_lum(src) 
  use m_mrgrnk, only: mrgrnk
  
  type(source_type), intent(inout) :: src(:) !< source array
  
  type(source_type), allocatable :: srclist(:)
  real(r8b), allocatable :: lumlist(:)
  integer(i8b), allocatable :: order(:)
  integer(i8b) :: i,N
  integer(i8b) :: err
  real(r8b) :: Ltot
  
  !      write(*,'(A)') "sorting the sources from brightest to dimmest ... "
  
  Ltot = 0.0
  N = size(src)

  allocate( srclist(N), stat=err )
  if(err/=0) stop "  order_sources_lum> cant allocate srclist"
      
  allocate( lumlist(N) , stat=err)
  if(err/=0) stop "  order_sources_lum> cant allocate lumlist"
    
  allocate( order(N) , stat=err)
  if(err/=0) stop "  order_sources_lum> cant allocate order"

  srclist = src(1:N)
  lumlist = src(1:N)%L
  call mrgrnk(-lumlist(1:N), order)
  
  ! sorting the negative of the luminosities puts the brightest first

  Ltot = 0.0d0
  do i = 1,N
     src(i) = srclist(order(i))
     Ltot = Ltot + src(i)%L
  end do
 
  lumlist = lumlist / Ltot
  
  do i = 1,N
     src(i)%Lcdf = lumlist(order(i))
  end do
  
  do i = 2,N
     src(i)%Lcdf = src(i-1)%Lcdf + src(i)%Lcdf 
  end do
  
  deallocate( srclist, lumlist, order )

end subroutine order_sources_lum


!> forms a snapshot name from a path, a file base, a snapshot number
!! and a file number.
!===================================================================
subroutine form_snapshot_file_name(Path,FileBase,SnapNum,FileNum,SnapFile)
  character(*), intent(in) :: Path       !< path to snapshot dir
  character(*), intent(in) :: FileBase   !< file base names
  integer(i4b), intent(in) :: SnapNum      !< snapshot number
  integer(i4b), intent(in) :: FileNum      !< file number in snapshot
  character(*), intent(out) :: SnapFile  !< file name to return
  
  character(10) :: FileNumChar
  
  write(FileNumChar,"(I6)") FileNum
  
  100 format(A,"/",A,"_",I3.3,".",A)
  write(SnapFile,100) trim(Path), trim(FileBase), SnapNum, &
                      trim(adjustl(FileNumChar))
  
end subroutine form_snapshot_file_name


!> outputs source snapshot header to the screen
!-----------------------------------------------
subroutine source_header_to_screen(shead)
  type(source_header_type), intent(in) :: shead !< source header to print

97 format(T2,"=",T58,"=")
98 format(T2,"= ",A,T58,"=")
99 format(T2,57("="))
100 format(T2,"=",T4,A,T30,I6,T58,"=")
200 format(T2,"=",T4,A,T30,ES12.5,T58,"=")
201 format(T2,"=",T4,A,T58,"=")
  
  write(*,99) 
  write(*,98) "source header data"
  write(*,99) 
  write(*,97) 
  write(*,100) "sources in snapshot:", shead%NsrcSnap
  write(*,100) "sources in file:", shead%NsrcFile
  write(*,100) "number of files in snap:", shead%NumFiles
  write(*,200) "rays for this snap:*", shead%TotalRays
  write(*,201) "*(only if RayScheme='header' in config file)"
  write(*,200) "luminosity unit (photons/s):", shead%Lunit
  write(*,97) 
  write(*,99) 
  
end subroutine source_header_to_screen


!> outputs source snapshot header to file
!-----------------------------------------------
subroutine source_header_to_file(shead,lun)
  type(source_header_type), intent(in) :: shead !< source header to print
  integer(i4b), intent(in) :: lun

97 format(T2,"=",T58,"=")
98 format(T2,"= ",A,T58,"=")
99 format(T2,57("="))
100 format(T2,"=",T4,A,T30,I6,T58,"=")
200 format(T2,"=",T4,A,T30,ES12.5,T58,"=")
201 format(T2,"=",T4,A,T58,"=")
  
  write(lun,99) 
  write(lun,98) "source header data"
  write(lun,99) 
  write(lun,97) 
  write(lun,100) "sources in snapshot:", shead%NsrcSnap
  write(lun,100) "sources in file:", shead%NsrcFile
  write(lun,100) "number of files in snap:", shead%NumFiles
  write(lun,200) "rays for this snap:*", shead%TotalRays
  write(lun,201) "*(only if RayScheme='header' in config file)"
  write(lun,200) "luminosity unit (photons/s):", shead%Lunit
  write(lun,97) 
  write(lun,99) 
  
end subroutine source_header_to_file

!> outputs the source data currently loaded in the particle system to screen
!! -------------------------------------------------------------------------
subroutine source_info_to_screen(psys,str,lun)

  type(particle_system_type), intent(in) :: psys     !< particle system
  character(*), optional, intent(in) :: str          !< arbitrary string
  integer(i4b), optional, intent(in) :: lun          !< if present goes to file
  integer(i4b) :: outlun

  integer(i8b), parameter :: srclimit = 20     !< max number of sources to screen
  integer(i8b) :: i
  
  outlun=stdout
  if (present(lun)) outlun=lun

  100 format(72("-"))
  110 format(T2,A,3ES10.3)
  111 format(T2,A,ES10.3)
  
  write(outlun,100)
  if (present(str)) write(outlun,"(A)") trim(str)
  write(outlun,"(A,I4,A,I13,A)") "source data for first ", srclimit, &
                                 " of ", size(psys%src), "  sources"
  write(outlun,*)   
  write(outlun,*) 
  120 format(T1,A, T15,A,     T37,A,     T51,A,     T65,A)
  121 format(T1,I2,T5,3ES10.3,T37,ES10.3,T51,ES10.3,T65,I3)
  write(outlun,120) "Src","Position","Luminosity","Spectrum","Emis Prf"
!  write(outlun,100)
  
  do i = 1, size(psys%src)
     write(outlun,121) i, psys%src(i)%pos, psys%src(i)%L, psys%src(i)%SpcType, &
          psys%src(i)%EmisPrf
     if (i==srclimit) exit
  end do
  
  write(outlun,100) 
  write(outlun,*)
  
  
end subroutine source_info_to_screen

end module source_input_mod
