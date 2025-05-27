!> \file gadget_public_header_class.f90

!> \brief Handles GADGET 2.0 Public version style headers.  
!!
!! Contains the means to read/write headers from a file or open lun.  
!! In addition, provides types to handle units and constants in a 
!! standard way.  Note that where 6 element arrays are needed that
!! correspond to particle types, we have indexed them from 0-5 as 
!! opposed to the default Fortran 1-6. 
!< 

module gadget_public_header_class
use myf03_mod
use gadget_general_class

#ifdef useMPI
  use mpi
#endif

implicit none
private

public :: gadget_public_header_type
public :: gadget_public_header_read_lun
public :: gadget_public_header_read_file
public :: gadget_public_header_write_lun
public :: gadget_public_header_print_lun
public :: gadget_public_header_return_gyr
public :: gadget_public_header_broadcast



!> Header type
!-----------------------------------------
type gadget_public_header_type
   integer(i4b) :: npar_file(0:5)   !< number of particles in snapshot file
   real(r8b)    :: mass(0:5)        !< mass of each particle type if constant
   real(r8b)    :: a                !< scale factor or time
   real(r8b)    :: z                !< redshift
   integer(i4b) :: flag_sfr         !< flag for star formation
   integer(i4b) :: flag_feedback    !< flag for feedback
   integer(i4b) :: npar_all(0:5)    !< number of particles in whole snapshot
   integer(i4b) :: flag_cooling     !< flag for radiative cooling
   integer(i4b) :: nfiles           !< number of files in a this snapshot
   real(r8b)    :: boxlen           !< box length
   real(r8b)    :: OmegaM           !< omega matter
   real(r8b)    :: OmegaL           !< omega lambda
   real(r8b)    :: h                !< little hubble
   integer(i4b) :: flag_age         !< flag for stellar age
   integer(i4b) :: flag_metals      !< flag for metallicity
   integer(i4b) :: npar_hw(0:5)     !< 64 bit part of npar 
   integer(i4b) :: flag_entr_ics    !< ICs contain entropy instead of energy
   integer(i4b) :: unused(15)       !< padding to make 256 bytes
end type gadget_public_header_type



contains


!> uses the scale factor and OmegaM+OmegaL to 
!! calculate the time in Gyr since the Big Bang
!----------------------------------------------------
function gadget_public_header_return_gyr(this) result(t_gyr)
  type(gadget_public_header_type) :: this
  type(gadget_constants_type) :: const
  real(r8b) :: t_gyr
  real(r8b) :: aeq, pre, arg, H0, km_per_Mpc, sec_per_Gyr

  if (this%OmegaL == 0.0) then
     t_gyr = 1.0d0
     return
  endif

  km_per_Mpc = const%CM_PER_MPC / 1.0d5 
  sec_per_Gyr = const%SEC_PER_MEGAYEAR * 1.0d3 

  H0 = this%h * 100.0d0 / km_per_Mpc  ! 1/s
  aeq = ( this%OmegaM/this%OmegaL )**(1.0d0/3.0d0)
  pre = 2.0d0 / ( 3.0d0 * sqrt(this%OmegaL) )
  arg = (this%a/aeq)**(3.0d0/2.0d0) + sqrt(1 + (this%a/aeq)**3 )
  t_gyr = pre * log(arg) / H0  ! in s
  t_gyr = t_gyr / sec_per_Gyr 

end function gadget_public_header_return_gyr


!> reads a gadget header from an open raw binary file 
!--------------------------------------------------------------
subroutine gadget_public_header_read_lun(this, lun)
  type(gadget_public_header_type) :: this
  integer(i4b), intent(in) :: lun

  read(lun) this%npar_file(:), this%mass(:), this%a, this%z, &              
       this%flag_sfr, this%flag_feedback, this%npar_all(:), &    
       this%flag_cooling, this%nfiles, this%boxlen, this%OmegaM, &         
       this%OmegaL, this%h, this%flag_age, this%flag_metals, &    
       this%npar_hw(:), this%flag_entr_ics, this%unused(:)      

end subroutine gadget_public_header_read_lun



!> reads a gadget header from snapfile and closes it afterwards
!--------------------------------------------------------------
subroutine gadget_public_header_read_file(this, snapfile)
  type(gadget_public_header_type) :: this
  character(*), intent(in) :: snapfile
  integer(i4b) :: lun    

  call open_unformatted_file_r( snapfile, lun )
  call gadget_public_header_read_lun(this, lun)
  close(lun)

end subroutine gadget_public_header_read_file


!> writes a gadget header to an open file
!--------------------------------------------------------------
subroutine gadget_public_header_write_lun(this, lun)
  type(gadget_public_header_type) :: this
  integer(i4b), intent(in) :: lun

  write(lun) this%npar_file(:), this%mass(:), this%a, this%z, &              
       this%flag_sfr, this%flag_feedback, this%npar_all(:), &    
       this%flag_cooling, this%nfiles, this%boxlen, this%OmegaM, &         
       this%OmegaL, this%h, this%flag_age, this%flag_metals, &    
       this%npar_hw(:), this%flag_entr_ics, this%unused(:)      

end subroutine gadget_public_header_write_lun




!> formatted print of header to lun (including standard out)
!---------------------------------------------------------------
subroutine gadget_public_header_print_lun(this, lun)
  type(gadget_public_header_type), intent(in) :: this 
  integer(i4b), intent(in) :: lun
  integer(i8b) :: i

  character(clen) :: n1,n2,n3
  character(clen) :: clmn_fmt
  character(clen) :: type_fmt
  character(clen) :: totl_fmt
  character(clen) :: scal_fmt
  character(clen) :: flag_fmt
  character(clen) :: star_fmt
  character(clen) :: line_fmt

  clmn_fmt = "(T17,A,T32,A,T47,A)"
  type_fmt = "(A,I1,A,A,A,T17,A,T32,A,T47,A)"
  totl_fmt = "(A,T17,A,T47,A)"
  scal_fmt = "(A,A,T25,A,A,T45,A,A)"
  flag_fmt = "(A,6(A,I1))"
  star_fmt = "(78('='))"
  line_fmt = "(78('-'))"

  write(lun,*)
  write(lun,star_fmt)
  write(lun,clmn_fmt) "n_all", "mass", "n_file"
  write(lun,line_fmt)
  do i = 0,5
     write(n1,'(I20)') this%npar_all(i)
     write(n2,'(E12.5)') this%mass(i)
     write(n3,'(I20)') this%npar_file(i)
     write(lun,type_fmt) "type",i,"(",gadget_ptype_names(i),")",&
          trim(adjustl(n1)), trim(adjustl(n2)), trim(adjustl(n3))
  end do
  write(lun,line_fmt)
  write(n1,'(I20)') sum(this%npar_all)
  write(n2,'(I20)') sum(this%npar_file)
  write(lun,totl_fmt) "total:", trim(adjustl(n1)), trim(adjustl(n2)) 
  write(lun,*)

  write(n1,'(F12.5)') this%a
  write(n2,'(F12.5)') this%z
  write(n3,'(F12.5)') this%h

  write(lun,scal_fmt) "a = ", trim(adjustl(n1)), &
       "z = ", trim(adjustl(n2)), &
       "h = ", trim(adjustl(n3))

  write(n1,'(F12.5)') this%OmegaM
  write(n2,'(F12.5)') this%OmegaL
  write(n3,'(F12.5)') this%boxlen

  write(lun,scal_fmt) "OmegaM = ", trim(adjustl(n1)), &
       "OmegaL = ", trim(adjustl(n2)), &
       "BoxSize = ", trim(adjustl(n3))
  write(lun,*)

  write(n1,'(I20)') this%nfiles
  write(lun,'(A,A)') "nfiles = ", trim(adjustl(n1))
  write(lun,flag_fmt) "/flags/ ",&
       "  sfr=",this%flag_sfr, &
       ", feedback=",this%flag_feedback, &
       ", cooling=",this%flag_cooling, &
       ", age=",this%flag_age, &
       ", metals=",this%flag_metals, &
       ", entr_ics=", this%flag_entr_ics
  write(lun,*) 
  write(lun,star_fmt)

end subroutine gadget_public_header_print_lun






!> broadcasts header
!-----------------------------------------
subroutine gadget_public_header_broadcast( this )
  type(gadget_public_header_type) :: this
  integer :: count
  integer :: root
  integer :: ierr

  root = 0

#ifdef useMPI
  count = 6
  call mpi_bcast( this%npar_file, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( this%npar_all, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( this%npar_hw, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( this%mass, count, mpi_double_precision, root, mpi_comm_world, ierr )

  count = 1
  call mpi_bcast( this%a, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( this%z, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( this%flag_sfr, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( this%flag_feedback, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( this%flag_cooling, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( this%nfiles, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( this%boxlen, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( this%OmegaM, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( this%OmegaL, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( this%h, count, mpi_double_precision, root, mpi_comm_world, ierr )
  call mpi_bcast( this%flag_age, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( this%flag_metals, count, mpi_integer, root, mpi_comm_world, ierr )
  call mpi_bcast( this%flag_entr_ics, count, mpi_integer, root, mpi_comm_world, ierr )
#endif

end subroutine gadget_public_header_broadcast




end module gadget_public_header_class
