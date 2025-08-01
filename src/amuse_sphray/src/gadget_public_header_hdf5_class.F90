!> \file gadget_public_header_hdf5_class.F90

!> \brief Handles HDF5 GADGET 2.0 Public version style headers.  
!!
!! Contains the means to read/write headers from a file.  
!< 

module gadget_public_header_hdf5_class
use myf03_mod
use gadget_general_class
use gadget_public_header_class

#ifdef useHDF5
use hdf5_wrapper
#endif

implicit none
private
 
public :: gadget_public_header_hdf5_read_lun
public :: gadget_public_header_hdf5_read_file
public :: gadget_public_header_hdf5_write_lun


contains


!> reads a gadget header from an hdf5 file
!--------------------------------------------------------------
subroutine gadget_public_header_hdf5_read_lun(this, fh)
  type(gadget_public_header_type) :: this
  integer, intent(in) :: fh

#ifdef useHDF5
  call hdf5_read_attribute(fh,'Header/NumPart_ThisFile',this%npar_file)
  call hdf5_read_attribute(fh,'Header/NumPart_Total',this%npar_all)
  call hdf5_read_attribute(fh,'Header/NumPart_Total_HW',this%npar_hw)
  call hdf5_read_attribute(fh,'Header/MassTable',this%mass)
  call hdf5_read_attribute(fh,'Header/Time',this%a)

  call hdf5_read_attribute(fh,'Header/Redshift',this%z)
  call hdf5_read_attribute(fh,'Header/BoxSize',this%boxlen)
  call hdf5_read_attribute(fh,'Header/NumFilesPerSnapshot',this%nfiles)
  call hdf5_read_attribute(fh,'Header/Omega0',this%OmegaM)
  call hdf5_read_attribute(fh,'Header/OmegaLambda',this%OmegaL)
  call hdf5_read_attribute(fh,'Header/HubbleParam',this%h)

  call hdf5_read_attribute(fh,'Header/Flag_Sfr',this%flag_sfr)
  call hdf5_read_attribute(fh,'Header/Flag_Cooling',this%flag_cooling)
  call hdf5_read_attribute(fh,'Header/Flag_StellarAge',this%flag_age)
  call hdf5_read_attribute(fh,'Header/Flag_Metals',this%flag_metals)
  call hdf5_read_attribute(fh,'Header/Flag_Feedback',this%flag_feedback)
  call hdf5_read_attribute(fh,'Header/Flag_Entropy_ICs',this%flag_entr_ics)
#else 
  write(*,*) 'useHDF5 macro not defined in Makefile' 
  stop
#endif

end subroutine gadget_public_header_hdf5_read_lun


!> reads a gadget header from an hdf5 file
!--------------------------------------------------------------
subroutine gadget_public_header_hdf5_read_file(this, snapfile)
  type(gadget_public_header_type) :: this
  character(*), intent(in) :: snapfile
  integer :: fh

#ifdef useHDF5
  call hdf5_open_file( fh, snapfile, readonly=.true. )
  call gadget_public_header_hdf5_read_lun(this, fh)
  call hdf5_close_file( fh )
#else 
  write(*,*) 'useHDF5 macro not defined in Makefile' 
  stop
#endif

end subroutine gadget_public_header_hdf5_read_file



!> writes a gadget header to an hdf5 file
!--------------------------------------------------------------
subroutine gadget_public_header_hdf5_write_lun(this, fh)
  type(gadget_public_header_type), intent(in) :: this
  integer, intent(in) :: fh

#ifdef useHDF5
  call hdf5_write_attribute(fh,'Header/NumPart_ThisFile',this%npar_file)
  call hdf5_write_attribute(fh,'Header/NumPart_Total',this%npar_all)
  call hdf5_write_attribute(fh,'Header/NumPart_Total_HW',this%npar_hw)
  call hdf5_write_attribute(fh,'Header/MassTable',this%mass)
  call hdf5_write_attribute(fh,'Header/Time',this%a)

  call hdf5_write_attribute(fh,'Header/Redshift',this%z)
  call hdf5_write_attribute(fh,'Header/BoxSize',this%boxlen)
  call hdf5_write_attribute(fh,'Header/NumFilesPerSnapshot',this%nfiles)
  call hdf5_write_attribute(fh,'Header/Omega0',this%OmegaM)
  call hdf5_write_attribute(fh,'Header/OmegaLambda',this%OmegaL)
  call hdf5_write_attribute(fh,'Header/HubbleParam',this%h)

  call hdf5_write_attribute(fh,'Header/Flag_Sfr',this%flag_sfr)
  call hdf5_write_attribute(fh,'Header/Flag_Cooling',this%flag_cooling)
  call hdf5_write_attribute(fh,'Header/Flag_StellarAge',this%flag_age)
  call hdf5_write_attribute(fh,'Header/Flag_Metals',this%flag_metals)
  call hdf5_write_attribute(fh,'Header/Flag_Feedback',this%flag_feedback)
  call hdf5_write_attribute(fh,'Header/Flag_Entropy_ICs',this%flag_entr_ics)
#else 
  write(*,*) 'useHDF5 macro not defined in Makefile' 
  stop
#endif

end subroutine gadget_public_header_hdf5_write_lun



end module gadget_public_header_hdf5_class
