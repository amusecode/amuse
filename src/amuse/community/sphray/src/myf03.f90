!> \file myf03.f90

!> \brief My Fortran 2003/8 module
!!
!! Makes use of the Fortran 2003 and 2008 standards.  The features I've 
!! included here are supported by version 4.5 of the GNU Fortran compiler. 
!! The module contains several useful variables and subroutines. 
!! 1. preconnected logical unit numbers (lun) for standard in, out and error.
!! 2. portable kind type parameters for 1/2/4/8 byte integers 
!! 3. portable kind type parameters for 4/8/16 byte reals. 
!! 4. default character length
!! 5. a verbosity variable to be used with write wrapper
!! 6. a command line type to standardize recieving arguments
!! portable precision definitions and subroutines for reading and 
!! parsing files 
!<

module myf03_mod
#if (__GNUC__ > 3 && __GNUC_MINOR__ > 7) || __INTEL_COMPILER > 1200
use iso_fortran_env    
#endif
implicit none
private

public :: stderr, stdin, stdout
public :: i1b, i2b, i4b, i8b
public :: r4b, r8b
public :: clen 
public :: command_line_type

public :: myf03_verbosity  
public :: myf03_initialize_command_line
public :: get_free_lun
public :: get_env

public :: open_formatted_file_r
public :: open_unformatted_file_r
public :: open_unformatted_file_sr

public :: open_formatted_file_w
public :: open_unformatted_file_w
public :: open_unformatted_file_sw

public :: scanfile
public :: mywrite
public :: myerr

  interface scanfile
     module procedure scanfile_r4b, scanfile_r8b, scanfile_i4b, &
                      scanfile_i8b, scanfile_chr, scanfile_log
  end interface

! pre connected logical unit numbers
!------------------------------------------------------------------------
#if __GNUC__ > 3 && __GNUC_MINOR__ > 7 
  integer, parameter :: stderr = error_unit  !< preconnected std error lun 
  integer, parameter :: stdin  = input_unit  !< preconnected std in lun
  integer, parameter :: stdout = output_unit !< preconnected std out lun
#else
  integer, parameter :: stderr = 0  !< preconnected std error lun 
  integer, parameter :: stdin  = 5  !< preconnected std in lun
  integer, parameter :: stdout = 6  !< preconnected std out lun
#endif
! selected_int_kind(r)     -10^r < n < 10^r
!------------------------------------------------------------------------

  integer, parameter :: i1b = selected_int_kind(2)  !< 1 byte integer type
  integer, parameter :: i2b = selected_int_kind(4)  !< 2 byte integer type
  integer, parameter :: i4b = selected_int_kind(9)  !< 4 byte integer type
  integer, parameter :: i8b = selected_int_kind(18) !< 8 byte integer type

! selected_real_kind(p,r)  p is decimal precision, r is exponent range
!------------------------------------------------------------------------
  integer, parameter :: r4b  = selected_real_kind(p=6,r=37)    !< 4 byte real type
  integer, parameter :: r8b  = selected_real_kind(p=15,r=307)  !< 8 byte real type





  integer, parameter :: clen  = 512 !< default character variable length

  integer :: myf03_verbosity=3        !< global verbosity threshold

!> command line type. 
!----------------------------
  type command_line_type
     integer :: nargs       !< number of arguments excluding executable
     integer :: len         !< number of chars in whole command line 
     character(clen) :: str !< the whole command line  
     character(clen), allocatable :: args(:)  !< individual command line args
     integer, allocatable :: arglens(:)       !< number of chars in each arg
  end type command_line_type


contains



!> fills in the command line variable
!-----------------------------------------------
  function myf03_initialize_command_line(verb) result(cmnd)
    type(command_line_type) :: cmnd
    integer :: verb 
    integer :: stat
    integer :: i
    character(clen) :: str

    character(clen), parameter :: myname="myf03_initialize_command_line"
    logical, parameter :: crash = .true.


#if (__GNUC__ > 3 && __GNUC_MINOR__ > 7) || __INTEL_COMPILER > 1200
    cmnd%nargs = command_argument_count()
    
    call get_command( command=cmnd%str, length=cmnd%len )
#else
    cmnd%nargs = 0
    cmnd%str = ' '
    
#endif
    call mywrite('', verb) 
    call mywrite('command string: ' // trim(cmnd%str), verb)
    write(str,'(A,I5)') 'command nargs :', cmnd%nargs
    call mywrite( str, verb) 
    
    allocate( cmnd%args( 0:cmnd%nargs ), cmnd%arglens( 0:cmnd%nargs ) )

    do i = 0,cmnd%nargs
       call get_command_argument( i, &
            value=cmnd%args(i), &
            length=cmnd%arglens(i) )

       if (i==0) then
          write(str,'(A,A)') 'exe: ', trim(cmnd%args(i))
       else
          write(str,'(A,I2,A,A)') 'arg ', i, ': ', trim(cmnd%args(i))
       endif
       call mywrite(str, verb)

       call mywrite('', verb)
    end do
    
  end function myf03_initialize_command_line



!> reports an error message and stops execution
!-----------------------------------------------
   subroutine myerr(str,routine,crash)
     character(len=*), intent(in) :: str
     character(len=*), intent(in) :: routine
     logical, intent(in) :: crash

     write(*,*) "*****************************"
     write(*,*) "error detected"
     write(*,*) "routine: ", trim(routine)
     write(*,*) trim(str)
     write(*,*) "*****************************"
     if (crash) stop

   end subroutine myerr



!> verbosity dependent write
!-----------------------------------------------
   subroutine mywrite( str, verb, lun, adv )
     character(len=*), intent(in) :: str
     integer, intent(in) :: verb
     integer, optional, intent(in) :: lun
     logical, optional, intent(in) :: adv
     character(3) :: sadv
     character(10) :: fmt

     character(clen), parameter :: myname='mywrite'
     logical, parameter :: crash=.true.

     if (present(adv)) then
        if (adv) then
           sadv = "yes"
        else
           sadv = "no "
        end if
     else
        sadv = "yes"
     end if

     fmt='(A)'

     if (present(lun)) then
        if (verb <= myf03_verbosity) write(lun,fmt,advance=sadv) trim(str) 
     else
        if (verb <= myf03_verbosity) write(stdout,fmt,advance=sadv) trim(str)
     endif

   end subroutine mywrite



!> returns a free logical unit number
!------------------------------------
   function get_free_lun() result(lun)

    integer(i4b) :: lun                        !< free lun
    integer(i4b) :: i                          !< loop counter
    integer(i4b), parameter :: minlun = 110    !< min lun to check
    integer(i4b), parameter :: maxlun = 1000   !< max lun to check
    logical :: badlun                          !< true if lun already connected

    character(clen), parameter :: myname="get_free_lun"
 
    do i = minlun,maxlun
      inquire(unit=i, opened=badlun)
      if (.not. badlun) then
         lun = i
         return
      end if
    end do

    write(*,*) '  checked luns from, ', minlun, ' to ', maxlun

    call myerr( str=" no free logical unit numbers", &
                routine=myname, &
                crash=.true. )

  end function get_free_lun



!> returns an environment variable
!----------------------------------------
  function get_env(env) result(var)
    character(len=*) :: env !< name of environment variable requested
    character(clen) :: var  !< value of environment variable requested

    call get_environment_variable(env,  value=var,  trim_name=.true.)

  end function get_env



!> opens a formatted file for reading and returns the lun
!--------------------------------------------------------
  subroutine open_formatted_file_r(filename,lun)
  
    character(len=*), intent(in) :: filename !< name of file to open
    integer(i4b), intent(out) :: lun     !< lun 

    lun = get_free_lun() 
    open(unit=lun, file=filename, action="read")

  end subroutine open_formatted_file_r



!> opens an unformatted file for reading and returns the lun
!-------------------------------------------------------------
  subroutine open_unformatted_file_r(filename,lun)
  
    character(len=*), intent(in) :: filename  !< name of file to open
    integer(i4b), intent(out)  :: lun     !< lun

    lun = get_free_lun()
    open(unit=lun, file=filename, form='unformatted', action="read")

  end subroutine open_unformatted_file_r

!> opens an unformatted file for stream reading and returns the lun
!--------------------------------------------------------------------
  subroutine open_unformatted_file_sr(filename,lun)
  
    character(len=*), intent(in) :: filename  !< name of file to open
    integer(i4b), intent(out)  :: lun     !< lun

    lun = get_free_lun()
    open(unit=lun, file=filename, form="unformatted", &
         action="read", access="stream", position="rewind")

  end subroutine open_unformatted_file_sr



!> opens a formatted file for writing and returns the lun
!--------------------------------------------------------
  subroutine open_formatted_file_w(filename,lun)
  
    character(len=*), intent(in) :: filename !< name of file to open
    integer(i4b), intent(out) :: lun     !< lun 

    lun = get_free_lun()
    open(unit=lun, file=filename, action="write")

  end subroutine open_formatted_file_w


!> opens an unformatted file for writing and returns the lun
!------------------------------------------------------------
  subroutine open_unformatted_file_w(filename,lun)
  
    character(len=*), intent(in) :: filename  !< name of file to open
    integer(i4b), intent(out)  :: lun             !< lun

    lun = get_free_lun()
    open(unit=lun, file=filename, form='unformatted', &
         action="write")
   
  end subroutine open_unformatted_file_w


!> opens an unformatted file for stream writing and returns the lun
!-------------------------------------------------------------------
  subroutine open_unformatted_file_sw(filename,lun)
  
    character(len=*), intent(in) :: filename  !< name of file to open
    integer(i4b), intent(out)  :: lun             !< lun

    lun = get_free_lun()
    open(unit=lun, file=filename, form='unformatted', &
         action="write", access="stream", position="rewind")

  end subroutine open_unformatted_file_sw



  ! scan file routines, one for each type
  !==================================================================
  subroutine scanfile_r4b(filename,keyword,var) 
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: keyword
    real(r4b), intent(out) :: var

    integer :: keylen
    integer :: strlen
    integer :: keypos
    character(clen) :: string
    integer :: ierr
    integer :: lun
    integer :: count

    ierr=0
    count=0
    keylen = len(adjustl(trim(keyword)))
    call open_formatted_file_r(filename,lun)
    do while( ierr==0 ) 
       read(lun,'(A)',iostat=ierr) string
       strlen = len(adjustl(trim(string)))
       keypos = index( string, adjustl(trim(keyword)) )
       if (keypos == 1 .and. ierr==0) then
          count = count + 1
          read(string(keypos+keylen:strlen),*) var
       endif
    end do
    close(lun)

    if (count/=1) then
       write(*,'(A,A)') ' keyword ', trim(keyword)
       write(*,'(A,A)') ' file:   ', trim(filename)
       if (count==0) then
          write(*,*) ' *** keyword not found in file'
          stop
       else if (count>1) then
          write(*,*) ' *** keyword found ', count, ' times'
          stop
       endif
    endif

  end subroutine scanfile_r4b


  subroutine scanfile_r8b(filename,keyword,var) 
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: keyword
    real(r8b), intent(out) :: var

    integer :: keylen
    integer :: strlen
    integer :: keypos
    character(clen) :: string
    integer :: ierr
    integer :: lun
    integer :: count

    ierr=0
    count=0
    keylen = len(adjustl(trim(keyword)))
    call open_formatted_file_r(filename,lun)
    do while( ierr==0 ) 
       read(lun,'(A)',iostat=ierr) string
       strlen = len(adjustl(trim(string)))
       keypos = index( string, adjustl(trim(keyword)) )
       if (keypos == 1 .and. ierr==0) then
          count = count + 1
          read(string(keypos+keylen:strlen),*) var
       endif
    end do
    close(lun)

    if (count/=1) then
       write(*,'(A,A)') ' keyword ', trim(keyword)
       write(*,'(A,A)') ' file:   ', trim(filename)
       if (count==0) then
          write(*,*) ' *** keyword not found in file'
          stop
       else if (count>1) then
          write(*,*) ' *** keyword found ', count, ' times'
          stop
       endif
    endif

  end subroutine scanfile_r8b


  subroutine scanfile_i4b(filename,keyword,var) 
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: keyword
    integer(i4b), intent(out) :: var

    integer :: keylen
    integer :: strlen
    integer :: keypos
    character(clen) :: string
    integer :: ierr
    integer :: lun
    integer :: count

    ierr=0
    count=0
    keylen = len(adjustl(trim(keyword)))
    call open_formatted_file_r(filename,lun)
    do while( ierr==0 ) 
       read(lun,'(A)',iostat=ierr) string
       strlen = len(adjustl(trim(string)))
       keypos = index( string, adjustl(trim(keyword)) )
       if (keypos == 1 .and. ierr==0) then
          count = count + 1
          read(string(keypos+keylen:strlen),*) var
       endif
    end do
    close(lun)

    if (count/=1) then
       write(*,'(A,A)') ' keyword ', trim(keyword)
       write(*,'(A,A)') ' file:   ', trim(filename)
       if (count==0) then
          write(*,*) ' *** keyword not found in file'
          stop
       else if (count>1) then
          write(*,*) ' *** keyword found ', count, ' times'
          stop
       endif
    endif

  end subroutine scanfile_i4b



  subroutine scanfile_i8b(filename,keyword,var) 
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: keyword
    integer(i8b), intent(out) :: var

    integer :: keylen
    integer :: strlen
    integer :: keypos
    character(clen) :: string
    integer :: ierr
    integer :: lun
    integer :: count

    ierr=0
    count=0
    keylen = len(adjustl(trim(keyword)))
    call open_formatted_file_r(filename,lun)
    do while( ierr==0 ) 
       read(lun,'(A)',iostat=ierr) string
       strlen = len(adjustl(trim(string)))
       keypos = index( string, adjustl(trim(keyword)) )
       if (keypos == 1 .and. ierr==0) then
          count = count + 1
          read(string(keypos+keylen:strlen),*) var
       endif
    end do
    close(lun)

    if (count/=1) then
       write(*,'(A,A)') ' keyword ', trim(keyword)
       write(*,'(A,A)') ' file:   ', trim(filename)
       if (count==0) then
          write(*,*) ' *** keyword not found in file'
          stop
       else if (count>1) then
          write(*,*) ' *** keyword found ', count, ' times'
          stop
       endif
    endif

  end subroutine scanfile_i8b



  subroutine scanfile_chr(filename,keyword,var)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: keyword
    character(clen), intent(out) :: var

    integer :: keylen
    integer :: strlen
    integer :: keypos
    character(clen) :: string
    integer :: ierr
    integer :: lun
    integer :: count

    ierr=0
    count=0
    keylen = len(adjustl(trim(keyword)))
    call open_formatted_file_r(filename,lun)
    do while( ierr==0 ) 
       read(lun,'(A)',iostat=ierr) string
       strlen = len(adjustl(trim(string)))
       keypos = index( string, adjustl(trim(keyword)) )
       if (keypos == 1 .and. ierr==0) then
          count = count + 1
          read(string(keypos+keylen:strlen),'(A)') var
          var = adjustl(trim(var))
       endif
    end do
    close(lun)

    if (count/=1) then
       write(*,'(A,A)') ' keyword ', trim(keyword)
       write(*,'(A,A)') ' file:   ', trim(filename)
       if (count==0) then
          write(*,*) ' *** keyword not found in file'
          stop
       else if (count>1) then
          write(*,*) ' *** keyword found ', count, ' times'
          stop
       endif
    endif
  end subroutine scanfile_chr



  subroutine scanfile_log(filename,keyword,var)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: keyword
    logical, intent(out) :: var

    integer :: keylen
    integer :: strlen
    integer :: keypos
    character(clen) :: string
    integer :: ierr
    integer :: lun
    integer :: count

    ierr=0
    count=0
    keylen = len(adjustl(trim(keyword)))
    call open_formatted_file_r(filename,lun)
    do while( ierr==0 ) 
       read(lun,'(A)',iostat=ierr) string
       strlen = len(adjustl(trim(string)))
       keypos = index( string, adjustl(trim(keyword)) )
       if (keypos == 1 .and. ierr==0) then
          count = count + 1
          read(string(keypos+keylen:strlen),*) var
       endif
    end do
    close(lun)

    if (count/=1) then
       write(*,'(A,A)') ' keyword ', trim(keyword)
       write(*,'(A,A)') ' file:   ', trim(filename)
       if (count==0) then
          write(*,*) ' *** keyword not found in file'
          stop
       else if (count>1) then
          write(*,*) ' *** keyword found ', count, ' times'
          stop
       endif
    endif

  end subroutine scanfile_log





end module myf03_mod


