!> \file output.F90

!> \brief The module that handles output
!<

module output_mod
use myf03_mod
use gadget_general_class
use gadget_public_header_class
use gadget_owls_header_class
use gadget_public_input_hdf5_mod
use gadget_sphray_header_class
use particle_system_mod
use oct_tree_mod, only: oct_tree_type
use physical_constants_mod
use global_mod, only: PLAN, GV
use global_mod, only: saved_gheads
use config_mod, only: write_config_hdf5_lun

#ifdef useHDF5
use hdf5_wrapper
#endif

implicit none

contains
 
!=======================================================================
!=======================================================================
!
! OUTPUT ROUTINES
!
!=======================================================================
!=======================================================================


!> put the config choices in a sphray header for output.
! ==========================================================
 
  subroutine config_to_ghead(fnum, ghead)


    integer, intent(in) :: fnum                    !< which snapshot file
    type(gadget_sphray_header_type), intent(out) :: ghead !< gadget header

    ghead%npar_all = 0
    ghead%npar_all(0) = saved_gheads(GV%CurSnapNum,fnum)%npar_all(0)     

    ghead%npar_file = 0
    ghead%npar_file(0) = saved_gheads(GV%CurSnapNum,fnum)%npar_file(0)     

    ghead%mass = 0.0

    ghead%nfiles    = saved_gheads(GV%CurSnapNum,fnum)%nfiles

    ghead%a = PLAN%snap(GV%CurSnapNum)%ScalefacAt
    ghead%z = saved_gheads(GV%CurSnapNum,fnum)%z
    ghead%time_gyr = GV%time_elapsed_myr / 1.0d3

    ghead%boxlen   = saved_gheads(GV%CurSnapNum,fnum)%boxlen
    ghead%OmegaM   = saved_gheads(GV%CurSnapNum,fnum)%OmegaM
    ghead%OmegaL   = saved_gheads(GV%CurSnapNum,fnum)%OmegaL
    ghead%OmegaB   = saved_gheads(GV%CurSnapNum,fnum)%OmegaB
    ghead%h        = saved_gheads(GV%CurSnapNum,fnum)%h

    ghead%npar_hw = 0
    ghead%npar_hw(1) = saved_gheads(GV%CurSnapNum,fnum)%npar_hw(1)
       
    ghead%flag_sfr      = saved_gheads(GV%CurSnapNum,fnum)%flag_sfr
    ghead%flag_feedback = saved_gheads(GV%CurSnapNum,fnum)%flag_feedback
    ghead%flag_cooling  = saved_gheads(GV%CurSnapNum,fnum)%flag_cooling
    ghead%flag_age      = saved_gheads(GV%CurSnapNum,fnum)%flag_age
    ghead%flag_metals   = saved_gheads(GV%CurSnapNum,fnum)%flag_metals
    ghead%flag_entr_ics = saved_gheads(GV%CurSnapNum,fnum)%flag_entr_ics
    
    ghead%rays_traced = GV%src_rayn

#ifdef incHmf
    ghead%flag_Hmf = 1
#else
    ghead%flag_Hmf = 0
#endif
!-------------------------
#ifdef incHemf
    ghead%flag_Hemf = 1
#else
    ghead%flag_Hemf = 0
#endif
!-------------------------
#ifdef incHe
    ghead%flag_helium = 1
#else
    ghead%flag_helium = 0
#endif
!-------------------------
#ifdef outGammaHI
    ghead%flag_gammaHI = 1
#else
    ghead%flag_gammaHI = 0
#endif
!-------------------------
#ifdef incCloudy
    ghead%flag_cloudy = 1
#else
    ghead%flag_cloudy = 0
#endif
!-------------------------
#ifdef incEOS
    ghead%flag_eos = 1
#else
    ghead%flag_eos = 0
#endif
!-------------------------
#ifdef incSFR
    ghead%flag_incsfr = 1
#else
    ghead%flag_incsfr = 0
#endif
!-------------------------

    ghead%unused = 0

  end subroutine config_to_ghead

!> outputs the whole shebang to a file
!=======================================

  subroutine output_total_snap(psys)

     character(clen), parameter :: myname="output_total_snap"
     logical, parameter :: crash = .true.
     integer, parameter :: verb = 2

     type(particle_system_type), intent(inout) :: psys
     type(gadget_sphray_header_type) :: ghead
     type(gadget_constants_type) :: gconst
     type(gadget_units_type) :: gunits

     character(3) :: label
     character(4) :: ext
     character(clen) :: filename
     integer(i4b) :: lun
     integer(i8b) :: Nread, Nfile
     integer(i4b) :: ipar, ifile

     integer(i8b) :: i,j

     real(r8b) :: scale
     real(r8b) :: hub
     real(r8b) :: nHe_over_nH
     real(r8b) :: Hmf, Hemf
     real(r4b) :: mu
     real(r4b), allocatable :: uint(:)
     real(r4b), allocatable :: rblock3(:,:)

     character(clen) :: tag
     character(clen) :: group_name

     ! set defaults
     !==============
     nHe_over_nH = 0.0
     ipar=0
     group_name = 'PartType0/'

     scale = PLAN%snap(GV%CurSnapNum)%ScalefacAt
     hub   = GV%LittleH

     write(*,*) 'output number: ', GV%OutputIndx
     write(*,*) 'output type:   ', GV%OutputType
     write(*,*) "writing total state of system"
     write(*,*) "time (elapsed code) ", GV%time_elapsed_code
     write(*,*) "time (elapsed myr)  ", GV%time_elapsed_myr

     if (GV%Comoving) then
        call particle_system_scale_physical_to_comoving(psys, scale, hub)
     endif

     100 format(I3.3)
     write(label,100) GV%OutputIndx


     call particle_system_set_ye( psys, GV%H_mf, GV%He_mf, GV%NeBackground )


     !================================================================
     ! GADGET standard formatted output
     !================================================================
     !================================================================
     if (GV%OutputType == 1) then
        Nread = 0

        do ifile = 0, GV%ParFilesPerSnap-1

           call config_to_ghead( ifile, ghead )

           !    form file name
           write(ext,'(I3)') ifile

           filename = trim(GV%OutputDir) // "/" // trim(GV%OutputFileBase) // "_" // label 
           if (GV%ParFilesPerSnap > 1) then
              filename = trim(filename) // "." // trim(adjustl(ext))
           end if
           write(*,*) "writing snapshot state to ", trim(filename)
           write(*,*) 
           call open_unformatted_file_w(filename,lun)

           Nfile = ghead%npar_file(0)

           call gadget_sphray_header_write_lun(ghead,lun)

           allocate( rblock3(3,Nfile) )
           rblock3(1,:) = psys%par(Nread+1:Nread+Nfile)%pos(1) 
           rblock3(2,:) = psys%par(Nread+1:Nread+Nfile)%pos(2) 
           rblock3(3,:) = psys%par(Nread+1:Nread+Nfile)%pos(3) 
           write(lun) rblock3

#ifdef incVel
           rblock3(1,:) = psys%par(Nread+1:Nread+Nfile)%vel(1) 
           rblock3(2,:) = psys%par(Nread+1:Nread+Nfile)%vel(2) 
           rblock3(3,:) = psys%par(Nread+1:Nread+Nfile)%vel(3) 
#else
           rblock3 = 0.0e0
#endif

           write(lun) rblock3
           deallocate( rblock3 )

           write(lun) psys%par(Nread+1:Nread+Nfile)%id 
           write(lun) psys%par(Nread+1:Nread+Nfile)%mass


           ! do some conversions
           !------------------------
           allocate( uint(Nfile) )

           j=1
           do i = Nread+1, Nread+Nfile
!------------------------------
#ifdef incHmf
              Hmf=psys%par(i)%Hmf
#else
              Hmf=GV%H_mf
#endif     
!-------------------------------
              mu      = 4.0d0 / ( 3.0d0 * Hmf + 1.0d0 + 4.0d0 * Hmf * psys%par(i)%ye )
              uint(j) = ( gconst%BOLTZMANN * psys%par(i)%T ) / &
                        ( (gconst%GAMMA - 1.0d0) * mu * gconst%PROTONMASS  )
              uint(j) = uint(j) * GV%cgs_mass / GV%cgs_enrg

              j = j+1
           end do


           write(lun) uint(1:Nfile)
           deallocate( uint )

           write(lun) psys%par(Nread+1:Nread+Nfile)%rho         
           write(lun) psys%par(Nread+1:Nread+Nfile)%ye
           write(lun) psys%par(Nread+1:Nread+Nfile)%xHI
           write(lun) psys%par(Nread+1:Nread+Nfile)%hsml 

           write(lun) psys%par(Nread+1:Nread+Nfile)%T


#ifdef incHmf
           write(lun) psys%par(Nread+1:Nread+Nfile)%Hmf
#endif


#ifdef incHemf
           write(lun) psys%par(Nread+1:Nread+Nfile)%Hemf
#endif


#ifdef incHe
           write(lun) psys%par(Nread+1:Nread+Nfile)%xHeI
           write(lun) psys%par(Nread+1:Nread+Nfile)%xHeII
#endif


#ifdef outGammaHI
           do ipar = Nread+1, Nread+Nfile
              if (psys%par(ipar)%time > 0.0) then
                 psys%par(ipar)%gammaHI = psys%par(ipar)%gammaHI / psys%par(ipar)%time
              else
                 psys%par(ipar)%gammaHI = 0.0
              end if
           end do
           write(lun) psys%par(Nread+1:Nread+Nfile)%gammaHI

           psys%par(Nread+1:Nread+Nfile)%gammaHI = 0.0
           psys%par(Nread+1:Nread+Nfile)%time = 0.0
#endif

 
#ifdef incCloudy
           write(lun) psys%par(Nread+1:Nread+Nfile)%xHI_cloudy
#endif

#ifdef incEOS
           write(lun) psys%par(Nread+1:Nread+Nfile)%eos
#endif

#ifdef incSFR
           write(lun) psys%par(Nread+1:Nread+Nfile)%sfr
#endif

           write(lun) psys%par(Nread+1:Nread+Nfile)%lasthit


           close(lun)
           Nread = Nread + Nfile

        end do

     end if
 

     !================================================================     
     ! GADGET HDF5 formatted output
     !================================================================
     !================================================================
     if (GV%OutputType == 2) then

#ifdef useHDF5

        Nread = 0

        do ifile = 0, GV%ParFilesPerSnap-1

           call config_to_ghead( ifile, ghead )

           !    form file name
           write(ext,'(I3)') ifile

           filename = trim(GV%OutputDir) // "/" // trim(GV%OutputFileBase) // "_" // label 
           if (GV%ParFilesPerSnap > 1) then
              filename = trim(filename) // "." // trim(adjustl(ext))
           end if
           filename = trim(filename) // ".hdf5"
           call mywrite("writing snapshot state to "//trim(filename), verb)
           call mywrite('',verb) 

           Nfile = ghead%npar_file(0)

           call hdf5_create_file( lun, filename )


           !================================================================
           ! write Header
           !================================================================
           call hdf5_create_group( lun, 'Header/' )
           call gadget_sphray_header_write_lun( ghead, lun )
!           call ghead%write_Gsphray_header_hdf5_lun(lun) 

           !================================================================
           ! write Units
           !================================================================
           call hdf5_create_group( lun, 'Units/' )
           gunits%cgs_length   = GV%cgs_len
           gunits%cgs_mass     = GV%cgs_mass
           gunits%cgs_velocity = GV%cgs_vel
           gunits%cgs_density  = GV%cgs_rho
           gunits%cgs_energy   = GV%cgs_enrg
           gunits%cgs_pressure = GV%cgs_prs
           gunits%cgs_time     = GV%cgs_time
           call gadget_units_write_lun( gunits, lun )          
!           call gunits%write_Gowls_units_lun(lun) 

           call gadget_data_attributes_set(gunits)
                     
           !================================================================
           ! write config parameters
           !================================================================
           call hdf5_create_group( lun, 'Config/' )
           call write_config_hdf5_lun(lun)

           if (Nfile == 0) then
              call hdf5_close_file(lun)
              cycle
           endif
           
           !================================================================
           ! write data
           !================================================================
           allocate( rblock3(3,Nfile) )
           rblock3(1,:) = psys%par(Nread+1:Nread+Nfile)%pos(1) 
           rblock3(2,:) = psys%par(Nread+1:Nread+Nfile)%pos(2) 
           rblock3(3,:) = psys%par(Nread+1:Nread+Nfile)%pos(3) 
           tag = trim(group_name)//'Coordinates'
           call hdf5_write_data( lun, tag, rblock3 )
           call gadget_data_attributes_write_lun( pos_attrs, lun, group_name, 'Coordinates' )
!           call pos_attrs%write_lun( lun, group_name, 'Coordinates' )


#ifdef incVel
           rblock3(1,:) = psys%par(Nread+1:Nread+Nfile)%vel(1) 
           rblock3(2,:) = psys%par(Nread+1:Nread+Nfile)%vel(2) 
           rblock3(3,:) = psys%par(Nread+1:Nread+Nfile)%vel(3) 
#else
           rblock3 = 0.0e0
#endif
           tag = trim(group_name)//'Velocity'
           call hdf5_write_data( lun, tag, rblock3 )
           call gadget_data_attributes_write_lun( vel_attrs, lun, group_name, 'Coordinates' )
!           call vel_attrs%write_lun( lun, group_name, 'Velocity' )
           deallocate( rblock3 )



           tag = trim(group_name)//'ParticleIDs'
           call hdf5_write_data( lun, tag, psys%par(Nread+1:Nread+Nfile)%id )
           call gadget_data_attributes_write_lun( id_attrs, lun, group_name, 'Coordinates' )
!           call id_attrs%write_lun( lun, group_name, 'ParticleIDs' )

           tag = trim(group_name)//'Mass'
           call hdf5_write_data( lun, tag, psys%par(Nread+1:Nread+Nfile)%mass )
           call gadget_data_attributes_write_lun( mass_attrs, lun, group_name, 'Coordinates' )
!           call mass_attrs%write_lun( lun, group_name, 'Mass' )


           tag = trim(group_name)//'Density'
           call hdf5_write_data( lun, tag, psys%par(Nread+1:Nread+Nfile)%rho )
           call gadget_data_attributes_write_lun( rho_attrs, lun, group_name, 'Coordinates' )
!           call rho_attrs%write_lun( lun, group_name, 'Density' )


           tag = trim(group_name)//'ElectronFraction'
           call hdf5_write_data( lun, tag, psys%par(Nread+1:Nread+Nfile)%ye )
           call gadget_data_attributes_write_lun( ye_attrs, lun, group_name, 'Coordinates' )
!           call ye_attrs%write_lun( lun, group_name, 'ElectronFraction' )

 
           tag = trim(group_name)//'HydrogenOneFraction'
           call hdf5_write_data( lun, tag, psys%par(Nread+1:Nread+Nfile)%xHI )
           call gadget_data_attributes_write_lun( xHI_attrs, lun, group_name, 'Coordinates' )
!           call xHI_attrs%write_lun( lun, group_name, 'HydrogenOneFraction' )


           tag = trim(group_name)//'SmoothingLength'
           call hdf5_write_data( lun, tag, psys%par(Nread+1:Nread+Nfile)%hsml )
           call gadget_data_attributes_write_lun( hsml_attrs, lun, group_name, 'Coordinates' ) 
!           call hsml_attrs%write_lun( lun, group_name, 'SmoothingLength' )


           tag = trim(group_name)//'Temperature'
           call hdf5_write_data( lun, tag, psys%par(Nread+1:Nread+Nfile)%T )
           call gadget_data_attributes_write_lun( T_attrs, lun, group_name, 'Coordinates' )
!           call T_attrs%write_lun( lun, group_name, 'Temperature' )


          
#ifdef incHmf
           tag = trim(group_name)//'ElementAbundance/Hydrogen'
           call hdf5_write_data( lun, tag, psys%par(Nread+1:Nread+Nfile)%Hmf )
           call gadget_data_attributes_write_lun( Hmf_attrs, lun, group_name, 'Coordinates' )
!           call Hmf_attrs%write_lun( lun, group_name, 'ElementAbundance/Hydrogen' )
#endif


#ifdef incHemf
           tag = trim(group_name)//'ElementAbundance/Helium'
           call hdf5_write_data( lun, tag, psys%par(Nread+1:Nread+Nfile)%Hemf )
           call gadget_data_attributes_write_lun( Hemf_attrs, lun, group_name, 'Coordinates' )
!           call Hemf_attrs%write_lun( lun, group_name, 'ElementAbundance/Helium' )
#endif


#ifdef incHe
           tag = trim(group_name)//'HeliumOneFraction'
           call hdf5_write_data( lun, tag, psys%par(Nread+1:Nread+Nfile)%xHeI )
           call gadget_data_attributes_write_lun( xHeI_attrs, lun, group_name, 'Coordinates' )
!           call xHeI_attrs%write_lun( lun, group_name, 'HeliumOneFraction' )

           tag = trim(group_name)//'HeliumTwoFraction'
           call hdf5_write_data( lun, tag, psys%par(Nread+1:Nread+Nfile)%xHeII )
           call gadget_data_attributes_write_lun( xHeII_attrs, lun, group_name, 'Coordinates' )
!           call xHeII_attrs%write_lun( lun, group_name, 'HeliumTwoFraction' )
#endif



#ifdef outGammaHI
           tag = trim(group_name)//'HydrogenOneGamma'
           call hdf5_write_data( lun, tag, psys%par(Nread+1:Nread+Nfile)%gammaHI )
           call gadget_data_attributes_write_lun( gammaHI_attrs, lun, group_name, 'Coordinates' )
!           call gammaHI_attrs%write_lun( lun, group_name, 'HydrogenOneGamma' )
#endif


 
#ifdef incCloudy
           tag = trim(group_name)//'HydrogenOneFraction_Cloudy'
           call hdf5_write_data( lun, tag, psys%par(Nread+1:Nread+Nfile)%xHI_cloudy )
           call gadget_data_attributes_write_lun( xHI_cloudy_attrs, lun, group_name, 'Coordinates' )
!           call xHI_cloudy_attrs%write_lun( lun, group_name, 'HydrogenOneFraction_Cloudy' )
#endif



#ifdef incEOS
           tag = trim(group_name)//'OnEquationOfState'
           call hdf5_write_data( lun, tag, psys%par(Nread+1:Nread+Nfile)%eos )
           call gadget_data_attributes_write_lun( eos_attrs, lun, group_name, 'Coordinates' ) 
!           call eos_attrs%write_lun( lun, group_name, 'OnEquationOfState' )
#endif


#ifdef incSFR
           tag = trim(group_name)//'StarFormationRate'
           call hdf5_write_data( lun, tag, psys%par(Nread+1:Nread+Nfile)%sfr )
           call gadget_data_attributes_write_lun( sfr_attrs, lun, group_name, 'Coordinates' )
!           call sfr_attrs%write_lun( lun, group_name, 'StarFormationRate' )
#endif


#ifdef incH2
           tag = trim(group_name)//'MolecularHydrogenMassFraction'
           call hdf5_write_data( lun, tag, psys%par(Nread+1:Nread+Nfile)%fH2 )
           call gadget_data_attributes_write_lun( fh2_attrs, lun, group_name, 'Coordinates' )
!           call fh2_attrs%write_lun( lun, group_name, 'MolecularHydrogenMassFraction' )
#endif


           call hdf5_close_file(lun)
           Nread = Nread + Nfile

        end do

#else

        stop " *** want hdf5 output w/o defining useHDF5 macro *** "

#endif
        
    

     end if
 
     !================================================================
     !================================================================

     if (GV%Comoving) then
        call particle_system_scale_comoving_to_physical(psys, scale, hub)
     endif


  end subroutine output_total_snap



  
!> writes the global ionization state of the system to a file
!=============================================================
  subroutine ion_frac_out(psys,tree)
  use iliev_comparison_project_mod

     type(particle_system_type), intent(in) :: psys    !< particle system
     type(oct_tree_type), intent(in) :: tree

     real :: Nionfrac, Mionfrac, Vionfrac
     real :: CallsPerCross
     real :: minGHI, maxGHI, GHI
     integer :: i
     integer :: nothit

     100 format(A,T24,A,T44,A,T64,A)
     101 format(A,T44,A,T64,A)
     105 format(T21,ES12.3,T43,ES12.3,T63,ES12.3)
     106 format(T43,ES12.3,T63,ES12.3)
     110 format(T21,I15,T43,I15,T63,ES11.5)

     write(*,100) "time:", "code units", "Myrs", "seconds"
     write(*,105) GV%start_time_code + GV%itime * GV%dt_code, &
                  GV%start_time_myr  + GV%itime * GV%dt_myr, &
                  GV%start_time_s    + GV%itime * GV%dt_s
     write(*,*) 

     write(*,100) "time elapsed:", "code units", "Myrs", "seconds"
     write(*,105) GV%time_elapsed_code, &
                  GV%time_elapsed_myr, &
                  GV%time_elapsed_s
     write(*,*) 


!    calculate the ionized number, volume, and mass fractions
     Nionfrac = particle_system_mean_xHII_number_weight(psys)
     Mionfrac = particle_system_mean_xHII_mass_weight(psys, DfltH_mf=GV%H_mf) 
     Vionfrac = particle_system_mean_xHII_volume_weight(psys)

     GV%nwionfrac = Nionfrac
     GV%mwionfrac = Mionfrac
     GV%vwionfrac = Vionfrac


     write(*,100) "neutral fraction:", "number weighted", "mass weighted", "volume weighted"
     write(*,105) 1.0d0-Nionfrac, 1.d0-Mionfrac, 1.0d0-Vionfrac
     write(*,*) 

     150 format (6ES15.5)
     write(GV%ionlun,150) GV%start_time_myr, &
                          GV%time_elapsed_myr, &
                          1.0d0-Nionfrac, 1.0d0-Mionfrac, 1.0d0-Vionfrac

     write(*,100) "rays cast:", "source", "diffuse", "diffuse/source"
     write(*,105) GV%TotalSourceRaysCast, GV%TotalDiffuseRaysCast, &
                  GV%TotalDiffuseRaysCast/GV%TotalSourceRaysCast 
     write(*,*) 
     
     if (GV%TotalPhotonsCast > 0.) then
        write(*,100) "photons:", "total cast", "per second", "% leaving box"
        write(*,105) GV%TotalPhotonsCast, GV%IonizingPhotonsPerSec, &
                     100.0 * GV%PhotonsLeavingBox / GV%TotalPhotonsCast
     else
        write(*,*) "zero luminosity source test"
     end if
     write(*,*)

     if (GV%ParticleCrossings .GT. 0) then
        CallsPerCross = GV%TotalDerivativeCalls / GV%ParticleCrossings
     else
        CallsPerCross = 0.0
     end if
     write(*,101) "solver(since last screen output):", &
                  "par crossings", "evals/crossings"
     write(*,106) GV%ParticleCrossings, CallsPerCross
     GV%TotalDerivativeCalls = 0.0
     GV%ParticleCrossings = 0.0

                  
     write(*,*)
     write(*,*) "Peak Updates = ", GV%PeakUpdates
     write(*,*) 



     160 format(T1,A,T20,ES12.5,T34,ES12.5)
     161 format(T1,A,T20,ES12.5)
     162 format(T1,A,T20,I12,T34,I12)
     write(*,160) "Min/Max xHI      = ", minval(psys%par(:)%xHI), &
                                        maxval(psys%par(:)%xHI)

     write(*,160) "Min/Max xHII     = ", minval(psys%par(:)%xHII), &
                                        maxval(psys%par(:)%xHII)

#ifdef incHe
     write(*,160) "Min/Max xHeI     = ", minval(psys%par(:)%xHeI), &
                                        maxval(psys%par(:)%xHeI)
     write(*,160) "Min/Max xHeII    = ", minval(psys%par(:)%xHeII), &
                                        maxval(psys%par(:)%xHeII)
     write(*,160) "Min/Max xHeIII   = ", minval(psys%par(:)%xHeIII), &
                                        maxval(psys%par(:)%xHeIII)
#endif

     write(*,160) "Min/Max T        = ", minval(psys%par(:)%T), &
                                        maxval(psys%par(:)%T)

#ifdef outGammaHI
     minGHI = huge(1.0)
     maxGHI = tiny(1.0)
     do i = 1, size(psys%par(:))
        if (psys%par(i)%time > 0.0) then
           GHI = psys%par(i)%gammaHI / psys%par(i)%time
           if (GHI < minGHI) minGHI = GHI
           if (GHI > maxGHI) maxGHI = GHI
        endif
     enddo
     write(*,160) "Min/Max GHI      = ", minGHI, maxGHI
#endif

     write(*,162) "Min/Max LastHit  = ", minval(psys%par(:)%lasthit), &
                                        maxval(psys%par(:)%lasthit)


     nothit = 0
     do i = 1, size(psys%par)
        if (psys%par(i)%lasthit == 0) then
           nothit = nothit + 1
        endif
     end do
     write(*,161) "Fraction Not Hit = ", real(nothit) / size(psys%par)


! do test specific outputs
     if (GV%DoTestScenario) then
        if (GV%TestScenario=="iliev_test1") then
           call iliev_test1_screen_out(psys,tree,GV)
        else if (GV%TestScenario=="iliev_test2") then
           call iliev_test2_screen_out(psys,tree,GV)
        else if (GV%TestScenario=="iliev_test3") then
           call iliev_test3_screen_out(psys,GV)
        else if (GV%TestScenario=="iliev_test4") then
           call iliev_test4_screen_out(psys,GV)
        end if
     end if

     write(*,*) 

  end subroutine ion_frac_out




end module output_mod
