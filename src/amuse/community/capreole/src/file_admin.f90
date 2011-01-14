!>
!! \brief This module contains parameters ans routines having to do
!!  with file in/ouput
!!
!! Module for C2-Ray/Capreole (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2008-05-27
!!
!! This module is also accepted by the F compiler (Dec 9, 2003)
!<
module file_admin

  implicit none

  private

  !> directory where files with results are placed
  character(len=50),public,parameter :: results_dir="./"

  !> set to not 5 if input via file (not implemented in all codes, check main program!)  
  integer,public,parameter :: stdinput=5 
  integer,public,parameter :: log_unit=30 !< unit number of log file(s)
  integer,public,parameter :: ah3=40  !< unit number of ah3 files (if applicable)
  integer,public,parameter :: iterdump=50  !< unit number of iterdump files (if applicable)
  logical,public :: file_input=.false. !< input from file?

  public :: flag_for_file_input

contains
  
  !> Set the file_input logical to distinguish between input
  !! through file from input through standard input 
  subroutine flag_for_file_input(flag)

    logical,intent(in) :: flag !< true if input from file

    ! Sets flag for inpit via file (and not via std input)
    file_input=flag
    
  end subroutine flag_for_file_input

end module file_admin
