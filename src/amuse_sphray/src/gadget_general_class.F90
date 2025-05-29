!> \file gadget_general_class.f90

!> \brief Handles universal GADGET stuff.  
!!
!! Provides types to handle units and constants in a 
!! standard way.  Note that where 6 element arrays are needed that
!! correspond to particle types, we have indexed them from 0-5 as 
!! opposed to the default Fortran 1-6. 
!< 

module gadget_general_class
use myf03_mod

#ifdef useMPI
  use mpi
#endif

#ifdef useHDF5
  use hdf5_wrapper
#endif

implicit none
private


public :: gadget_ptype_names

public :: gadget_units_type
public :: gadget_units_set
public :: gadget_units_print_lun
public :: gadget_units_read_lun
public :: gadget_units_read_file
public :: gadget_units_write_lun
public :: gadget_units_broadcast

public :: gadget_constants_type
public :: gadget_constants_read_lun
public :: gadget_constants_read_file
public :: gadget_constants_write_lun

public :: gadget_data_attributes_type
public :: gadget_data_attributes_set
public :: gadget_data_attributes_read_lun
public :: gadget_data_attributes_write_lun

public :: form_gadget_snapshot_file_name




!> Particle type names
!-----------------------------------------
character(5), parameter :: gadget_ptype_names(0:5) = &
     (/"gas  ","halo ","disk ","bulge","stars","bndry"/)


!> Units type (default= kpc/h, 1.0e10 Msolar/h, km/s)
!-----------------------------------------------------
type gadget_units_type
   real(r8b) :: cgs_length   = 3.0856780d21  !<  [cm h^-1]
   real(r8b) :: cgs_mass     = 1.989d43      !<  [g h^-1]
   real(r8b) :: cgs_velocity = 1.0d5         !<  [cm s^-1]
   real(r8b) :: cgs_time     = 3.085678d16   !<  [s h^-1]
   real(r8b) :: cgs_density  = 6.769911d-22  !<  [g cm^-3 h^2]
   real(r8b) :: cgs_pressure = 6.769911d-12  !<  [ba = g cm^-1 s^-2 h^2]
   real(r8b) :: cgs_energy   = 1.989d53      !<  [erg = g cm^2 s^-2 h^-1] 
end type gadget_units_type


!> Constants type ( 23 doubles )
!-----------------------------------------
type gadget_constants_type
   real(r8b) :: PI = 3.1415927d0                 !< [pure]                  
   real(r8b) :: GAMMA = 1.6666667d0              !< [pure]
   real(r8b) :: GRAVITY = 6.6720000d-08          !< [cm^3 g^-1s^-2]    
   real(r8b) :: SOLAR_MASS = 1.9890000d+33       !< [g]                        
   real(r8b) :: SOLAR_LUM = 3.8260000d+33        !< [erg/s]
   real(r8b) :: RAD_CONST = 7.5650000d-15        !< [erg cm^-3 K^-4]  
   real(r8b) :: AVOGADRO = 6.0222000d+23         !< [pure] 
   real(r8b) :: BOLTZMANN = 1.3806000d-16        !< [cm^2 g s^-2 K^-1] 
   real(r8b) :: GAS_CONST = 83142500.d0          !< [erg K^-1 mol^-1]  
   real(r8b) :: C = 2.9979000d+10                !< [cm/s]   
   real(r8b) :: PLANCK = 6.6262000d-27           !< [cm^2 g s^-1]  
   real(r8b) :: CM_PER_MPC = 3.0856780d+24       !< [pure]    
   real(r8b) :: PROTONMASS = 1.6726000d-24       !< [g]
   real(r8b) :: ELECTRONMASS = 9.1095300d-28     !< [g] 
   real(r8b) :: ELECTRONCHARGE = 4.8032000d-10   !< [esu]
   real(r8b) :: HUBBLE = 3.2407789d-18           !< [s^-1 h]
   real(r8b) :: T_CMB0 = 2.7280000d0             !< [K]
   real(r8b) :: SEC_PER_MEGAYEAR = 3.1550000d+13 !< [pure]
   real(r8b) :: SEC_PER_YEAR = 31550000.d0       !< [pure]
   real(r8b) :: STEFAN = 7.5657000d-15     !< = rad. const. [erg cm^-3 K^-4]  
   real(r8b) :: THOMPSON = 6.6524587d-25   !< [cm^2]
   real(r8b) :: EV_TO_ERG = 1.6021765d-12  !< [pure]
   real(r8b) :: Z_SOLAR = 0.012663729d0    !< [Mass Fraction] 
end type gadget_constants_type


!> Attributes for data in HDF5 files.
!-----------------------------------------------------
type gadget_data_attributes_type
   real(r8b) :: CGSConversionFactor 
   real(r4b) :: h_scale_exponent
   real(r4b) :: aexp_scale_exponent
   character(200) :: VarDescription
end type gadget_data_attributes_type



type(gadget_data_attributes_type), public :: pos_attrs
type(gadget_data_attributes_type), public :: id_attrs
type(gadget_data_attributes_type), public :: mass_attrs
type(gadget_data_attributes_type), public :: T_attrs
type(gadget_data_attributes_type), public :: rho_attrs
type(gadget_data_attributes_type), public :: ye_attrs
type(gadget_data_attributes_type), public :: xHI_attrs
type(gadget_data_attributes_type), public :: hsml_attrs
type(gadget_data_attributes_type), public :: lasthit_attrs

type(gadget_data_attributes_type), public :: vel_attrs
type(gadget_data_attributes_type), public :: xHI_cloudy_attrs
type(gadget_data_attributes_type), public :: Hmf_attrs
type(gadget_data_attributes_type), public :: xHeI_attrs
type(gadget_data_attributes_type), public :: xHeII_attrs
type(gadget_data_attributes_type), public :: Hemf_attrs
type(gadget_data_attributes_type), public :: gammaHI_attrs
type(gadget_data_attributes_type), public :: eos_attrs
type(gadget_data_attributes_type), public :: sfr_attrs


contains


!!============================================
!!
!!    UNITS    
!!
!!============================================



!> sets user supplied units
!--------------------------------------------------------------
subroutine gadget_units_set(this, cgs_length, cgs_mass, cgs_velocity)
  type(gadget_units_type) :: this
  real(r8b) :: cgs_length
  real(r8b) :: cgs_mass
  real(r8b) :: cgs_velocity

  this%cgs_length   = cgs_length
  this%cgs_mass     = cgs_mass
  this%cgs_velocity = cgs_velocity

  this%cgs_density  = this%cgs_mass / this%cgs_length**3
  this%cgs_energy   = this%cgs_mass * this%cgs_velocity**2
  this%cgs_time     = this%cgs_length / this%cgs_velocity
  this%cgs_pressure = this%cgs_mass / &
       (this%cgs_length**3 / this%cgs_velocity**2)

end subroutine gadget_units_set


!> formatted print of units to lun (including standard out)
!---------------------------------------------------------------
subroutine gadget_units_print_lun(this, lun, h)
  type(gadget_units_type), intent(in) :: this 
  integer(i4b), intent(in) :: lun
  real(r8b), intent(in) :: h

  type(gadget_constants_type) :: gconst

  character(clen) :: n1
  character(clen) :: star_fmt
  character(clen) :: unit_fmt
  character(clen) :: line_fmt

  real(r8b), parameter :: cm_per_km = 1.0d5

  ! binding energy of 1.0d8 [Msolar/h] halo @ z=0 
  ! Delta_c = 18 pi^2, units = erg/h
  ! http://arxiv.org/pdf/astro-ph/0010468v3
  real(r8b), parameter :: E8 = 5.45d53 * 1.0d-1

  real(r8b) :: rho_crit_0

  rho_crit_0 = 3.0d0 * gconst%hubble**2 / (8.0d0 * gconst%pi * gconst%gravity)

  star_fmt = "(78('='))"
  line_fmt = "(78('-'))"
  unit_fmt = "(A,T30,A)"

  write(lun,*)
  write(lun,star_fmt)
  write(lun,*) "   Units"
  write(lun,line_fmt)
  write(lun,*) "  h = ", h, " = H0[km/s/Mpc] / 100 "
  write(n1,'(ES20.6)') this%cgs_length
  write(lun,unit_fmt) "length [cm/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_mass
  write(lun,unit_fmt) "mass [g/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_velocity
  write(lun,unit_fmt) "velocity [cm/s]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_density
  write(lun,unit_fmt) "density [g/cm^3 h^2]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_energy
  write(lun,unit_fmt) "energy [erg/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_time
  write(lun,unit_fmt) "time [s/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_pressure
  write(lun,unit_fmt) "pressure [ba h^2]:", trim(adjustl(n1)) 

  write(lun,*) 

  write(n1,'(ES20.6)') rho_crit_0
  write(lun,unit_fmt) "rho_crit_0 [g/cm^3 h^2]:", trim(adjustl(n1))
  write(n1,'(ES20.6)') E8
  write(lun,unit_fmt) "E8 [erg/h]:", trim(adjustl(n1))

  write(lun,*)

  write(n1,'(ES20.6)') this%cgs_length / gconst%cm_per_mpc
  write(lun,unit_fmt) "length [Mpc/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_mass / gconst%solar_mass
  write(lun,unit_fmt) "mass [Msolar/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_velocity / cm_per_km
  write(lun,unit_fmt) "velocity [km/s]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_density / rho_crit_0
  write(lun,unit_fmt) "density [rho_crit_0 h^2]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_energy / E8
  write(lun,unit_fmt) "energy [E8]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_time / gconst%sec_per_megayear
  write(lun,unit_fmt) "time [Myr/h]:", trim(adjustl(n1)) 

  write(n1,'(ES20.6)') this%cgs_pressure / gconst%boltzmann 
  write(lun,unit_fmt) "pressure/k_b [cm^-3 K h^2]:", trim(adjustl(n1)) 

  write(lun,star_fmt)

end subroutine gadget_units_print_lun


!> reads units from an open hdf5 file
!--------------------------------------------------------------
subroutine gadget_units_read_lun(this, fh)
  type(gadget_units_type), intent(out) :: this
  integer, intent(in) :: fh
  
#ifdef useHDF5
  call hdf5_read_attribute(fh, 'Units/UnitLength_in_cm', this%cgs_length)
  call hdf5_read_attribute(fh, 'Units/UnitMass_in_g', this%cgs_mass)
  call hdf5_read_attribute(fh, 'Units/UnitVelocity_in_cm_per_s', this%cgs_velocity)
  call hdf5_read_attribute(fh, 'Units/UnitDensity_in_cgs', this%cgs_density)
  call hdf5_read_attribute(fh, 'Units/UnitEnergy_in_cgs', this%cgs_energy)
  call hdf5_read_attribute(fh, 'Units/UnitPressure_in_cgs', this%cgs_pressure)
  call hdf5_read_attribute(fh, 'Units/UnitTime_in_s', this%cgs_time)
#else
  write(*,*) 'useHDF5 macro not defined in Makefile'
  stop
#endif
  
end subroutine gadget_units_read_lun


!> reads units from an hdf5 file
!--------------------------------------------------------------
subroutine gadget_units_read_file(this, snapfile)
  type(gadget_units_type), intent(out) :: this
  character(*), intent(in) :: snapfile
  integer :: fh
  
#ifdef useHDF5
  call hdf5_open_file( fh, snapfile, readonly=.true. )
  call gadget_units_read_lun(this, fh)
  call hdf5_close_file( fh )
#else
  write(*,*) 'useHDF5 macro not defined in Makefile'
  stop
#endif

end subroutine gadget_units_read_file


!> writes units to an open hdf5 file
!--------------------------------------------------------------
subroutine gadget_units_write_lun(this, fh)
  type(gadget_units_type), intent(in) :: this
  integer, intent(in) :: fh

#ifdef useHDF5
  call hdf5_write_attribute(fh, 'Units/UnitLength_in_cm', this%cgs_length)
  call hdf5_write_attribute(fh, 'Units/UnitMass_in_g', this%cgs_mass)
  call hdf5_write_attribute(fh, 'Units/UnitVelocity_in_cm_per_s', this%cgs_velocity)
  call hdf5_write_attribute(fh, 'Units/UnitDensity_in_cgs', this%cgs_density)
  call hdf5_write_attribute(fh, 'Units/UnitEnergy_in_cgs', this%cgs_energy)
  call hdf5_write_attribute(fh, 'Units/UnitPressure_in_cgs', this%cgs_pressure)
  call hdf5_write_attribute(fh, 'Units/UnitTime_in_s', this%cgs_time)
#endif

end subroutine gadget_units_write_lun


!> broadcasts units
!-----------------------------------------
subroutine gadget_units_broadcast( units )
  type(gadget_units_type) :: units
  integer :: count
  integer :: root
  integer :: ierr

#ifdef useMPI
  count = 7
  root = 0
  call mpi_bcast( units, count, mpi_double_precision, root, mpi_comm_world, ierr )
  if (ierr /= 0) then
     call mpi_barrier(mpi_comm_world, ierr)
     call mpi_finalize(ierr)
  endif
#else 
  write(*,*) 'useMPI macro not defined in Makefile'
  stop
#endif

end subroutine gadget_units_broadcast




!!============================================
!!
!!    CONSTANTS
!!
!!============================================


!> reads constants from an open hdf5 file
!--------------------------------------------------------------
subroutine gadget_constants_read_lun(this, fh)
  type(gadget_constants_type), intent(out) :: this
  integer, intent(in) :: fh
  
#ifdef useHDF5
  call hdf5_read_attribute(fh,'Constants/PI',this%pi)
  call hdf5_read_attribute(fh,'Constants/GAMMA',this%gamma)
  call hdf5_read_attribute(fh,'Constants/GRAVITY',this%gravity)
  call hdf5_read_attribute(fh,'Constants/SOLAR_MASS',this%solar_mass)
  call hdf5_read_attribute(fh,'Constants/SOLAR_LUM',this%solar_lum)
  call hdf5_read_attribute(fh,'Constants/RAD_CONST',this%rad_const)
  call hdf5_read_attribute(fh,'Constants/AVOGADRO',this%avogadro)
  call hdf5_read_attribute(fh,'Constants/BOLTZMANN',this%boltzmann)
  call hdf5_read_attribute(fh,'Constants/GAS_CONST',this%gas_const)
  call hdf5_read_attribute(fh,'Constants/C',this%c)
  call hdf5_read_attribute(fh,'Constants/PLANCK',this%planck)
  call hdf5_read_attribute(fh,'Constants/CM_PER_MPC',this%cm_per_mpc)
  call hdf5_read_attribute(fh,'Constants/PROTONMASS',this%protonmass)
  call hdf5_read_attribute(fh,'Constants/ELECTRONMASS',this%electronmass)
  call hdf5_read_attribute(fh,'Constants/ELECTRONCHARGE',this%electroncharge)
  call hdf5_read_attribute(fh,'Constants/HUBBLE',this%hubble)
  call hdf5_read_attribute(fh,'Constants/T_CMB0',this%t_cmb0)
  call hdf5_read_attribute(fh,'Constants/SEC_PER_MEGAYEAR',this%sec_per_megayear)
  call hdf5_read_attribute(fh,'Constants/SEC_PER_YEAR',this%sec_per_year)
  call hdf5_read_attribute(fh,'Constants/STEFAN',this%stefan)
  call hdf5_read_attribute(fh,'Constants/THOMPSON',this%thompson)
  call hdf5_read_attribute(fh,'Constants/EV_TO_ERG',this%ev_to_erg)
  call hdf5_read_attribute(fh,'Constants/Z_Solar',this%z_solar)
#else 
  write(*,*) 'useHDF5 macro not defined in Makefile' 
  stop
#endif  

end subroutine gadget_constants_read_lun

  
!> reads constants from an hdf5 file
!--------------------------------------------------------------
subroutine gadget_constants_read_file(this, snapfile)
  type(gadget_constants_type), intent(out) :: this
  character(*), intent(in) :: snapfile
  integer :: fh

#ifdef useHDF5
  call hdf5_open_file( fh, snapfile, readonly=.true. )
  call gadget_constants_read_lun(this, fh)
  call hdf5_close_file( fh )
#else 
  write(*,*) 'useHDF5 macro not defined in Makefile' 
  stop
#endif  

end subroutine gadget_constants_read_file


!> writes constants to an open hdf5 file
!--------------------------------------------------------------
subroutine gadget_constants_write_lun(this, fh)
  type(gadget_constants_type), intent(in) :: this
  integer, intent(in) :: fh

#ifdef useHDF5
  call hdf5_write_attribute(fh,'Constants/PI',this%pi)
  call hdf5_write_attribute(fh,'Constants/GAMMA',this%gamma)
  call hdf5_write_attribute(fh,'Constants/GRAVITY',this%gravity)
  call hdf5_write_attribute(fh,'Constants/SOLAR_MASS',this%solar_mass)
  call hdf5_write_attribute(fh,'Constants/SOLAR_LUM',this%solar_lum)
  call hdf5_write_attribute(fh,'Constants/RAD_CONST',this%rad_const)
  call hdf5_write_attribute(fh,'Constants/AVOGADRO',this%avogadro)
  call hdf5_write_attribute(fh,'Constants/BOLTZMANN',this%boltzmann)
  call hdf5_write_attribute(fh,'Constants/GAS_CONST',this%gas_const)
  call hdf5_write_attribute(fh,'Constants/C',this%c)
  call hdf5_write_attribute(fh,'Constants/PLANCK',this%planck)
  call hdf5_write_attribute(fh,'Constants/CM_PER_MPC',this%cm_per_mpc)
  call hdf5_write_attribute(fh,'Constants/PROTONMASS',this%protonmass)
  call hdf5_write_attribute(fh,'Constants/ELECTRONMASS',this%electronmass)
  call hdf5_write_attribute(fh,'Constants/ELECTRONCHARGE',this%electroncharge)
  call hdf5_write_attribute(fh,'Constants/HUBBLE',this%hubble)
  call hdf5_write_attribute(fh,'Constants/T_CMB0',this%t_cmb0)
  call hdf5_write_attribute(fh,'Constants/SEC_PER_MEGAYEAR',this%sec_per_megayear)
  call hdf5_write_attribute(fh,'Constants/SEC_PER_YEAR',this%sec_per_year)
  call hdf5_write_attribute(fh,'Constants/STEFAN',this%stefan)
  call hdf5_write_attribute(fh,'Constants/THOMPSON',this%thompson)
  call hdf5_write_attribute(fh,'Constants/EV_TO_ERG',this%ev_to_erg)
  call hdf5_write_attribute(fh,'Constants/Z_Solar',this%z_solar)
#else 
  write(*,*) 'useHDF5 macro not defined in Makefile' 
  stop
#endif

end subroutine gadget_constants_write_lun





!!============================================
!!
!!    DATA ATTRIBUTES
!!
!!============================================



!> writes data attributes to an hdf5 file  
!------------------------------------------------    
subroutine gadget_data_attributes_write_lun( attr, fh, grp_tag, dat_tag  )
  type(gadget_data_attributes_type) :: attr
  integer(i4b) :: fh
  character(*) :: grp_tag
  character(*) :: dat_tag

  character(clen) :: tag

#ifdef useHDF5

  tag = trim(grp_tag) // trim(dat_tag) // '/CGSConversionFactor'
  call hdf5_write_attribute( fh, tag, attr%CGSConversionFactor )

  tag = trim(grp_tag) // trim(dat_tag) // '/h-scale-exponent'
  call hdf5_write_attribute( fh, tag, attr%h_scale_exponent )

  tag = trim(grp_tag) // trim(dat_tag) // '/aexp-scale-exponent'
  call hdf5_write_attribute( fh, tag, attr%aexp_scale_exponent )

  tag = trim(grp_tag) // trim(dat_tag) // '/VarDescription'
  call hdf5_write_attribute( fh, tag, attr%VarDescription )

#else

  stop " *** called write_data_attr_lun w/o defining useHDF5 *** "

#endif

end subroutine gadget_data_attributes_write_lun



!> reads data attributes from an hdf5 file  
!------------------------------------------------    
subroutine gadget_data_attributes_read_lun( attr, fh, grp_tag, dat_tag  )
  type(gadget_data_attributes_type) :: attr
  integer(i4b) :: fh
  character(*) :: grp_tag
  character(*) :: dat_tag

  character(clen) :: tag

#ifdef useHDF5

  tag = trim(grp_tag) // trim(dat_tag) // '/CGSConversionFactor'
  call hdf5_read_attribute( fh, tag, attr%CGSConversionFactor )

  tag = trim(grp_tag) // trim(dat_tag) // '/h-scale-exponent'
  call hdf5_read_attribute( fh, tag, attr%h_scale_exponent )

  tag = trim(grp_tag) // trim(dat_tag) // '/aexp-scale-exponent'
  call hdf5_read_attribute( fh, tag, attr%aexp_scale_exponent )

  tag = trim(grp_tag) // trim(dat_tag) // '/VarDescription'
  call hdf5_read_attribute( fh, tag, attr%VarDescription )

#else

  stop " *** called read_data_attr_lun w/o defining useHDF5 *** "

#endif

end subroutine gadget_data_attributes_read_lun





!> sets data attributes using supplied units
!------------------------------------------------    
subroutine gadget_data_attributes_set(gunits)
  type(gadget_units_type), intent(in) :: gunits

  pos_attrs%CGSConversionFactor = gunits%cgs_length
  vel_attrs%CGSConversionFactor = gunits%cgs_velocity
  id_attrs%CGSConversionFactor = 1.0d0
  mass_attrs%CGSConversionFactor = gunits%cgs_mass
  T_attrs%CGSConversionFactor = 1.0d0
  rho_attrs%CGSConversionFactor = gunits%cgs_density
  ye_attrs%CGSConversionFactor = 1.0d0
  xHI_attrs%CGSConversionFactor = 1.0d0
  hsml_attrs%CGSConversionFactor = gunits%cgs_length
  lasthit_attrs%CGSConversionFactor = 1.0d0

  xHI_cloudy_attrs%CGSConversionFactor = 1.0d0
  Hmf_attrs%CGSConversionFactor = 1.0d0
  xHeI_attrs%CGSConversionFactor = 1.0d0
  xHeII_attrs%CGSConversionFactor = 1.0d0
  Hemf_attrs%CGSConversionFactor = 1.0d0
  gammaHI_attrs%CGSConversionFactor = 1.0d0
  eos_attrs%CGSConversionFactor = 1.0d0
  sfr_attrs%CGSConversionFactor = 1.0d0



  pos_attrs%VarDescription = "Co-moving coordinates. Physical: r = a x = Coordinates h^-1 a U_L [cm]"
  vel_attrs%VarDescription = "Co-moving velocities. Physical v_p = a dx/dt  = Velocities a^1/2 U_V [cm/s]"
  id_attrs%VarDescription = "Unique particle identifier"
  mass_attrs%VarDescription = "Particle mass. Physical m = Mass h^-1 U_M [g]"
  T_attrs%VarDescription = "Temperature [K]"
  rho_attrs%VarDescription = "Co-moving mass densities. Physical rho = Densities h^2 a^-3 U_M U_L^-3 [g/cm^3]"
  ye_attrs%VarDescription = "Electron abundance = n_e / n_H"
  xHI_attrs%VarDescription = "Hydrogen neutral fraction = n_HI / n_H"
  hsml_attrs%VarDescription = "Co-moving smoothing length. Physical h = SmoothingLength h^-1 a U_L [cm]"
  lasthit_attrs%VarDescription = "Index of last ray to cross this particle"

  xHI_cloudy_attrs%VarDescription = "Hydrogen neutral fraction from Cloudy Tables xHI = n_HI / n_H"
  Hmf_attrs%VarDescription = "Hydrogen mass fraction = mass in atomic Hydrogen / particle mass"
  xHeI_attrs%VarDescription = "Helium I fraction = n_HeI / n_He"
  xHeII_attrs%VarDescription = "Helium II fraction = n_HeII / n_He"
  Hemf_attrs%VarDescription = "Helium mass fraction = mass in atomic Helium / particle mass"
  gammaHI_attrs%VarDescription = "Hydrogen photoionization rate [1/s]"
  eos_attrs%VarDescription = "Flag for effective equation os state. 1 if currently on EoS" // &
       ", 0 if has never been on EoS, -ExpansionFactor if left the EoS at ExpansionFactor"
  sfr_attrs%VarDescription = "Star formation rate. Physical sfr = StarformationRate SOLAR_MASS SEC_PER_YEAR^-1 [g/s]"



  pos_attrs%h_scale_exponent = -1.0d0
  vel_attrs%h_scale_exponent = 0.0d0
  id_attrs%h_scale_exponent = 0.0d0
  mass_attrs%h_scale_exponent = -1.0d0
  T_attrs%h_scale_exponent = 0.0d0
  rho_attrs%h_scale_exponent = 2.0d0
  ye_attrs%h_scale_exponent = 0.0d0
  xHI_attrs%h_scale_exponent = 0.0d0
  hsml_attrs%h_scale_exponent = -1.0d0
  lasthit_attrs%h_scale_exponent = 0.0d0

  xHI_cloudy_attrs%h_scale_exponent = 0.0d0
  Hmf_attrs%h_scale_exponent = 0.0d0
  xHeI_attrs%h_scale_exponent = 0.0d0
  xHeII_attrs%h_scale_exponent = 0.0d0
  Hemf_attrs%h_scale_exponent = 0.0d0
  gammaHI_attrs%h_scale_exponent = 0.0d0
  eos_attrs%h_scale_exponent = 0.0d0
  sfr_attrs%h_scale_exponent = 0.0d0



  pos_attrs%aexp_scale_exponent = 1.0d0
  vel_attrs%aexp_scale_exponent = 0.5d0
  id_attrs%aexp_scale_exponent = 0.0d0
  mass_attrs%aexp_scale_exponent = 0.0d0
  T_attrs%aexp_scale_exponent = 0.0d0
  rho_attrs%aexp_scale_exponent = -3.0d0
  ye_attrs%aexp_scale_exponent = 0.0d0
  xHI_attrs%aexp_scale_exponent = 0.0d0
  hsml_attrs%aexp_scale_exponent = 1.0d0
  lasthit_attrs%aexp_scale_exponent = 0.0d0

  xHI_cloudy_attrs%aexp_scale_exponent = 0.0d0
  Hmf_attrs%aexp_scale_exponent = 0.0d0
  xHeI_attrs%aexp_scale_exponent = 0.0d0
  xHeII_attrs%aexp_scale_exponent = 0.0d0
  Hemf_attrs%aexp_scale_exponent = 0.0d0
  gammaHI_attrs%aexp_scale_exponent = 0.0d0
  eos_attrs%aexp_scale_exponent = 0.0d0
  sfr_attrs%aexp_scale_exponent = 0.0d0


end subroutine gadget_data_attributes_set











!> forms a snapshot file name from, path, base, SnapNum, and FileNum.
!! For example, path = "/home/galtay/data/snapshots", base = "snap",
!! SnapNum = 2 and FileNum = 15 would return in SnapFile, 
!! "/home/galtay/data/snapshots/snap_002.15".  Setting the hdf5bool to 
!! true appends ".hdf5" to the file. 
!------------------------------------------------------------------------
subroutine form_gadget_snapshot_file_name(path, base, SnapNum, FileNum, &
     SnapFile, hdf5bool)

  character(*), intent(in)     :: path        !< path to snapshot dir
  character(*), intent(in)     :: base        !< base snapshot name
  integer(i4b), intent(in)     :: SnapNum     !< snapshot number
  integer(i4b), intent(in)     :: FileNum     !< file number of snapshot
  character(clen), intent(out) :: SnapFile    !< snapshot filename
  logical, intent(in)          :: hdf5bool    !< hdf5 file? 

  character(clen) :: SnapFileTmp
  character(10) :: FileNumChar
  character(clen) :: fmt
  logical :: Fthere

  write(FileNumChar,"(I6)") FileNum
  fmt = "(A,'/',A,'_',I3.3)"

  ! first write a file with no extension
  !--------------------------------------
  write(SnapFileTmp,fmt) trim(path), trim(base), SnapNum
  SnapFile = trim(SnapFileTmp)
  if (hdf5bool) SnapFile = trim(SnapFile) // ".hdf5"
  inquire( file=SnapFile, exist=Fthere )


  ! if the file number is 0 and a file with no extension exists then return
  ! otherwise append the FileNum to the file name
  !-------------------------------------------------------------------------
  if (FileNum == 0 .and. Fthere) then
     return
  else
     SnapFile = trim(SnapFileTmp) // "." // trim(adjustl(FileNumChar))
     if (hdf5bool) SnapFile = trim(SnapFile) // ".hdf5"
  end if


end subroutine form_gadget_snapshot_file_name










end module gadget_general_class


