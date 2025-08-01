!> \file gadget_owls_header_class.F90

!> \brief Handles HDF5 GADGET 3.0 OWLS/GIMIC style headers.  
!!
!! Contains the means to read/write headers from a file.  
!< 

module gadget_owls_header_class
use myf03_mod
use gadget_general_class
use gadget_public_header_class

#ifdef useHDF5
use hdf5_wrapper
#endif

implicit none
private

public :: gadget_owls_header_type

public :: gadget_owls_header_read_lun
public :: gadget_owls_header_read_file
public :: gadget_owls_header_write_lun
public :: gadget_owls_header_print_lun


!> Header type
!-----------------------------------------
type gadget_owls_header_type
   character(200) :: run_label      !< name of OWLS/GIMIC run
   integer(i4b)   :: npar_file(0:5) !< number of particles in snapshot file
   integer(i4b)   :: npar_all(0:5)  !< number of particles in whole snapshot
   integer(i4b)   :: npar_hw(0:5)   !< 64 bit part of npar 
   real(r8b)      :: mass(0:5)      !< mass of each particle type if constant
   real(r8b)      :: a              !< scale factor or time
   real(r8b)      :: time_gyr       !< time in Gyr
   real(r8b)      :: z              !< redshift
   real(r8b)      :: boxlen         !< box length
   integer(i4b)   :: nfiles         !< number of files in a this snapshot
   real(r8b)      :: OmegaM         !< omega matter
   real(r8b)      :: OmegaB         !< omega baryon
   real(r8b)      :: OmegaL         !< omega lambda
   real(r8b)      :: h              !< little hubble
   integer(i4b)   :: flag_sfr       !< flag for star formation
   integer(i4b)   :: flag_cooling   !< flag for radiative cooling
   integer(i4b)   :: flag_age       !< flag for stellar age
   integer(i4b)   :: flag_metals    !< flag for metallicity
   integer(i4b)   :: flag_feedback  !< flag for feedback
end type gadget_owls_header_type



contains
  



!> reads an OWLS/GIMIC gadget header from an open hdf5 file
!--------------------------------------------------------------
subroutine gadget_owls_header_read_lun(this, fh)
  type(gadget_owls_header_type) :: this
  integer, intent(in) :: fh

#ifdef useHDF5
  call hdf5_read_attribute(fh,'Header/RunLabel',this%run_label)
  call hdf5_read_attribute(fh,'Header/NumPart_ThisFile',this%npar_file)
  call hdf5_read_attribute(fh,'Header/NumPart_Total',this%npar_all)
  call hdf5_read_attribute(fh,'Header/NumPart_Total_HighWord',this%npar_hw)
  call hdf5_read_attribute(fh,'Header/MassTable',this%mass)
  call hdf5_read_attribute(fh,'Header/ExpansionFactor',this%a)

  call hdf5_read_attribute(fh,'Header/Time_GYR',this%time_gyr)
  call hdf5_read_attribute(fh,'Header/Redshift',this%z)
  call hdf5_read_attribute(fh,'Header/BoxSize',this%boxlen)
  call hdf5_read_attribute(fh,'Header/NumFilesPerSnapshot',this%nfiles)
  call hdf5_read_attribute(fh,'Header/Omega0',this%OmegaM)
  call hdf5_read_attribute(fh,'Header/OmegaBaryon',this%OmegaB)

  call hdf5_read_attribute(fh,'Header/OmegaLambda',this%OmegaL)
  call hdf5_read_attribute(fh,'Header/HubbleParam',this%h)
  call hdf5_read_attribute(fh,'Header/Flag_Sfr',this%flag_sfr)
  call hdf5_read_attribute(fh,'Header/Flag_Cooling',this%flag_cooling)
  call hdf5_read_attribute(fh,'Header/Flag_StellarAge',this%flag_age)
  call hdf5_read_attribute(fh,'Header/Flag_Metals',this%flag_metals)
  call hdf5_read_attribute(fh,'Header/Flag_Feedback',this%flag_feedback)
#else 
  write(*,*) 'useHDF5 macro not defined in Makefile' 
  stop
#endif

end subroutine gadget_owls_header_read_lun



!> reads an OWLS/GIMIC gadget header from an hdf5 file
!--------------------------------------------------------------
subroutine gadget_owls_header_read_file(this, snapfile)
  type(gadget_owls_header_type) :: this
  character(*), intent(in) :: snapfile
  integer :: fh

#ifdef useHDF5
  call hdf5_open_file( fh, snapfile, readonly=.true. )
  call gadget_owls_header_read_lun(this, fh)
  call hdf5_close_file( fh )
#else 
  write(*,*) 'useHDF5 macro not defined in Makefile' 
  stop
#endif

end subroutine gadget_owls_header_read_file


!> writes a gadget header from an open hdf5 file
!--------------------------------------------------------------
subroutine gadget_owls_header_write_lun(this, fh)
  type(gadget_owls_header_type), intent(in) :: this
  integer, intent(in) :: fh

#ifdef useHDF5
  call hdf5_write_attribute(fh,'Header/RunLabel',this%run_label)
  call hdf5_write_attribute(fh,'Header/NumPart_ThisFile',this%npar_file)
  call hdf5_write_attribute(fh,'Header/NumPart_Total',this%npar_all)
  call hdf5_write_attribute(fh,'Header/NumPart_Total_HW',this%npar_hw)
  call hdf5_write_attribute(fh,'Header/MassTable',this%mass)
  call hdf5_write_attribute(fh,'Header/ExpansionFactor',this%a)

  call hdf5_write_attribute(fh,'Header/Time_GYR',this%time_gyr)
  call hdf5_write_attribute(fh,'Header/Redshift',this%z)
  call hdf5_write_attribute(fh,'Header/BoxSize',this%boxlen)
  call hdf5_write_attribute(fh,'Header/NumFilesPerSnapshot',this%nfiles)
  call hdf5_write_attribute(fh,'Header/Omega0',this%OmegaM)
  call hdf5_write_attribute(fh,'Header/OmegaBaryon',this%OmegaB)

  call hdf5_write_attribute(fh,'Header/OmegaLambda',this%OmegaL)
  call hdf5_write_attribute(fh,'Header/HubbleParam',this%h)
  call hdf5_write_attribute(fh,'Header/Flag_Sfr',this%flag_sfr)
  call hdf5_write_attribute(fh,'Header/Flag_Cooling',this%flag_cooling)
  call hdf5_write_attribute(fh,'Header/Flag_StellarAge',this%flag_age)
  call hdf5_write_attribute(fh,'Header/Flag_Metals',this%flag_metals)
  call hdf5_write_attribute(fh,'Header/Flag_Feedback',this%flag_feedback)
#else 
  write(*,*) 'useHDF5 macro not defined in Makefile' 
  stop
#endif

end subroutine gadget_owls_header_write_lun




!> formatted print of header to lun (including standard out)
!---------------------------------------------------------------
subroutine gadget_owls_header_print_lun(this, lun)
  type(gadget_owls_header_type), intent(in) :: this 
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
  flag_fmt = "(A,5(A,I1))"
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
       ", metals=",this%flag_metals
  write(lun,*) 
  write(lun,star_fmt)

end subroutine gadget_owls_header_print_lun





end module gadget_owls_header_class
