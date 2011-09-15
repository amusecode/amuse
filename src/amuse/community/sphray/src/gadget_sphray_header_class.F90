!> \file gadget_sphray_header_class.F90

!> \brief Handles Gadget style headers with extra SPHRAY info.  
!!
!!
!< 

module gadget_sphray_header_class
use myf03_mod
use gadget_general_class
use gadget_public_header_class
use gadget_owls_header_class

#ifdef useHDF5
use hdf5_wrapper
#endif

implicit none
private

public :: gadget_sphray_header_type
public :: gadget_sphray_header_copy_public
public :: gadget_sphray_header_copy_owls
public :: gadget_sphray_header_write_lun
public :: gadget_sphray_header_hdf5_write_lun



!> Gadget header type
!-----------------------------------------
type gadget_sphray_header_type
   integer(i4b) :: npar_file(0:5)  !< number of particles in snapshot file
   real(r8b)    :: mass(0:5)       !< mass of each particle type if constant
   real(r8b)    :: a               !< scale factor or time
   real(r8b)    :: z               !< redshift
   integer(i4b) :: flag_sfr        !< flag for star formation
   integer(i4b) :: flag_feedback   !< flag for feedback
   integer(i4b) :: npar_all(0:5)   !< number of particles in whole snapshot
   integer(i4b) :: flag_cooling    !< flag for radiative cooling
   integer(i4b) :: nfiles          !< number of files in a this snapshot
   real(r8b)    :: boxlen          !< box length
   real(r8b)    :: OmegaM          !< omega matter
   real(r8b)    :: OmegaL          !< omega lambda
   real(r8b)    :: h               !< little hubble
   integer(i4b) :: flag_age        !< flag for stellar age
   integer(i4b) :: flag_metals     !< flag for metallicity
   integer(i4b) :: npar_hw(0:5)    !< 64 bit part of npar 
   integer(i4b) :: flag_entr_ics   !< ICs contain entropy instead of energy

   real(r8b)    :: OmegaB          !< omega baryon 
   integer(i8b) :: rays_traced     !< number of rays traced
   integer(i4b) :: flag_Hmf        !< output has Hydorgen mass fraction? 
   integer(i4b) :: flag_Hemf       !< output has Helium mass fraction? 
   integer(i4b) :: flag_helium     !< output has xHeI and xHeII?
   integer(i4b) :: flag_gammaHI    !< output has HI photoionization rate? 
   integer(i4b) :: flag_cloudy     !< output has Cloudy xHIeq? 
   integer(i4b) :: flag_eos        !< output has EOS info?
   integer(i4b) :: flag_incsfr     !< output has SFR info?
   real(r8b)    :: time_gyr        !< time since BB from cosmo variables [Gyr]
   integer(i4b) :: unused(2)       !< spacer
end type gadget_sphray_header_type


contains


!> writes a gadget header to an open file
!--------------------------------------------------------------
subroutine gadget_sphray_header_write_lun(this, lun)
  type(gadget_sphray_header_type), intent(in) :: this
  integer(i4b), intent(in) :: lun

  write(lun) this%npar_file(:), this%mass(:), this%a, this%z, &              
       this%flag_sfr, this%flag_feedback, this%npar_all(:), &    
       this%flag_cooling, this%nfiles, this%boxlen, this%OmegaM, &         
       this%OmegaL, this%h, this%flag_age, this%flag_metals, &    
       this%npar_hw(:), this%flag_entr_ics, this%OmegaB, &
       this%rays_traced, this%flag_Hmf, this%flag_Hemf, this%flag_helium, &
       this%flag_gammaHI, this%flag_cloudy, this%flag_eos, this%flag_incsfr, &
       this%time_gyr, this%unused(:)      

end subroutine gadget_sphray_header_write_lun


!> writes a gadget header to an open HDF5 file
!--------------------------------------------------------------
subroutine gadget_sphray_header_hdf5_write_lun(this, fh)
  type(gadget_sphray_header_type), intent(in) :: this
  integer(i4b), intent(in) :: fh

#ifdef useHDF5
  call hdf5_write_attribute(fh,'Header/NumPart_ThisFile',this%npar_file)
  call hdf5_write_attribute(fh,'Header/NumPart_Total',this%npar_all)
  call hdf5_write_attribute(fh,'Header/NumPart_Total_HighWord',this%npar_hw)
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
  call hdf5_write_attribute(fh,'Header/Flag_Entropy_In_ICs',this%flag_entr_ics)
  call hdf5_write_attribute(fh,'Header/RaysTraced',this%rays_traced)

  call hdf5_write_attribute(fh,'Header/Flag_Hmf',this%flag_Hmf)
  call hdf5_write_attribute(fh,'Header/Flag_Hemf',this%flag_Hemf)
  call hdf5_write_attribute(fh,'Header/Flag_Helium',this%flag_helium)
  call hdf5_write_attribute(fh,'Header/Flag_GammaHI',this%flag_gammaHI)
  call hdf5_write_attribute(fh,'Header/Flag_Cloudy',this%flag_cloudy)
  call hdf5_write_attribute(fh,'Header/Flag_EoS',this%flag_eos)
  call hdf5_write_attribute(fh,'Header/Flag_IncSFR',this%flag_incsfr)
#endif


end subroutine gadget_sphray_header_hdf5_write_lun


!> copies an OWLS/GIMIC header into the sphray style header
!--------------------------------------------------------------
subroutine gadget_sphray_header_copy_owls( this, owlshead )
  type(gadget_sphray_header_type) :: this
  type(gadget_owls_header_type) :: owlshead

  this%npar_file(0:5)   = owlshead%npar_file(0:5)   
  this%mass(0:5)        = owlshead%mass(0:5)        
  this%a                = owlshead%a                
  this%z                = owlshead%z                
  this%flag_sfr         = owlshead%flag_sfr         
  this%flag_feedback    = owlshead%flag_feedback    
  this%npar_all(0:5)    = owlshead%npar_all(0:5)    
  this%flag_cooling     = owlshead%flag_cooling     
  this%nfiles           = owlshead%nfiles           
  this%boxlen           = owlshead%boxlen           
  this%OmegaM           = owlshead%OmegaM           
  this%OmegaL           = owlshead%OmegaL           
  this%h                = owlshead%h                
  this%flag_age         = owlshead%flag_age         
  this%flag_metals      = owlshead%flag_metals      
  this%npar_hw(0:5)     = owlshead%npar_hw(0:5)     

  this%OmegaB           = owlshead%OmegaB
  this%time_gyr         = owlshead%time_gyr

end subroutine gadget_sphray_header_copy_owls


!> copies a Public header into the sphray style header
!--------------------------------------------------------------
subroutine gadget_sphray_header_copy_public( this, pubhead )
  type(gadget_sphray_header_type) :: this
  type(gadget_public_header_type) :: pubhead

  this%npar_file(0:5)   = pubhead%npar_file(0:5)   
  this%mass(0:5)        = pubhead%mass(0:5)        
  this%a                = pubhead%a                
  this%z                = pubhead%z                
  this%flag_sfr         = pubhead%flag_sfr         
  this%flag_feedback    = pubhead%flag_feedback    
  this%npar_all(0:5)    = pubhead%npar_all(0:5)    
  this%flag_cooling     = pubhead%flag_cooling     
  this%nfiles           = pubhead%nfiles           
  this%boxlen           = pubhead%boxlen           
  this%OmegaM           = pubhead%OmegaM           
  this%OmegaL           = pubhead%OmegaL           
  this%h                = pubhead%h                
  this%flag_age         = pubhead%flag_age         
  this%flag_metals      = pubhead%flag_metals      
  this%npar_hw(0:5)     = pubhead%npar_hw(0:5)     
  this%flag_entr_ics    = pubhead%flag_entr_ics    

end subroutine gadget_sphray_header_copy_public


end module gadget_sphray_header_class
