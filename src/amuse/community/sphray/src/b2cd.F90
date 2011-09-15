!> \file b2cd.F90

!> \brief Handles the impact parameter to column depth table

module b2cd_mod
use myf03_mod
implicit none
private

public :: b2cdfac
public :: read_b2cd_file

integer(i4b) :: Nentries  !< table entries (from b=0 to b=1)
real(r8b) :: db           !< spacing in impact parameter
real(r8b), allocatable :: b2cd_table(:)  !< stores line integrations
                                         !< 0:Nentries-1


contains


!>  read in the impact parameter -> column depth file
!----------------------------------------------------
subroutine read_b2cd_file(b2cd_file)

  character(clen), intent(in) :: b2cd_file  ! impact parameter to CD table

  integer, parameter :: verb = 1
  character(clen) :: str
  integer(i4b) :: lun, err, i
  logical :: fthere

  character(clen), parameter :: myname = 'read_b2cd_file'
  logical, parameter :: crash = .true.

  write(str,'(A,T27,A)') 'using b2cd file: ', trim(b2cd_file)
  call mywrite(str, verb) 

  inquire(file=b2cd_file, exist=fthere)
  if (.not. fthere) then
     call myerr('cant find b2cd file: ' // trim(b2cd_file), myname, crash)
  end if

  call open_formatted_file_r(b2cd_file,lun)
  read(lun,*) Nentries
  write(str,'(A,I5)') '  b2cd table entries:    ', Nentries
  call mywrite('', verb+1)
  call mywrite(str, verb+1) 

  allocate ( b2cd_table(0:Nentries-1), stat=err )
  if(err/=0) call myerr('cant allocate b2cd_table', myname, crash)
     
  str = '  reading in SPH impact parameter -> column depth table'
  call mywrite(str, verb+1)

  do i = 0,Nentries-1
     read(lun,*) b2cd_table(i)
  end do
  db = 1.0d0 / (Nentries - 1)

  close(lun)
  call sum_b2cd_table()


end subroutine read_b2cd_file

!> Interpolates from the line integral table
!---------------------------------------------
function interpolate_b2cd(i,xfrac) result(b2cd)
  character(clen), parameter :: myname="interpolate_b2cd"
  logical, parameter :: crash=.true.
  integer(i4b) :: i   !< entry just below evaluation point
  real(r8b) :: xfrac  !< fraction of the way to the next entry
  real(r8b) :: b2cd   !< interpolated value
  real(r8b) :: dy     !< difference between adjacent table entries
  
  if (i >= Nentries-1) then
     write(*,*) "i = ", i
     call myerr(" out of bounds in interpolate b2cd",myname,crash)
  end if

  dy = b2cd_table(i+1) - b2cd_table(i)
  b2cd = b2cd_table(i) + xfrac * dy

end function interpolate_b2cd


!> checks normalization of impact parameter to column depth table
!----------------------------------------------------------------
subroutine sum_b2cd_table()
use physical_constants_mod, only: pi

  real(r8b) :: fac,sum,a,b,ab2,xfrac
  integer(i4b) :: i

  integer, parameter :: verb=2
  character(clen) :: str

  sum = 0.0d0
  xfrac = 0.5d0
  
  do i = 0,Nentries-2
     a = i * db
     b = (i+1) * db
     ab2 = (b+a)/2.0d0

     fac = db/6.0d0 * ( b2cd_table(i) + &
                        4.0d0 * interpolate_b2cd(i,xfrac) + &
                        b2cd_table(i+1) )

     sum = sum + 2.0d0 * pi * ab2 * fac 
  end do

  write(str,'(A,ES10.4)') &
       "  normalization of b2cd table (should be close to 1.0) = ", sum
  
  call mywrite(str, verb) 
  call mywrite('', verb)
  
end subroutine sum_b2cd_table

!> scales the table entry to the column depth factor using the 
!! smoothing length, the code cgs length, and the impact parameter
!! normalized to the smoothing length
function b2cdfac(b_norm,hsml,cgs_len) result(cdfac)

  real(r8b) :: b_norm  !< normalized impact parameter
  real(r8b) :: hsml    !< smoothing length
  real(r8b) :: cgs_len !< cgs length
  real(r8b) :: cdfac   !< output column depth factor

  integer :: b_indx
  real(r8b) :: xfrac
  real(r8b) :: b2cd

  b_indx = floor( (Nentries-1) * b_norm)  ! same as floor( b_norm / db )
  xfrac = b_norm / db - b_indx

  if (b_indx >= Nentries-1) then
     b_indx = Nentries-2
     xfrac = 1.0d0
  end if

  b2cd = interpolate_b2cd( b_indx, xfrac )
  
  cdfac = b2cd / ( hsml*hsml*cgs_len*cgs_len )

end function b2cdfac


end module b2cd_mod
