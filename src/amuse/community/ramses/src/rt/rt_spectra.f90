! These are modules for reading, integrating and interpolating RT-relevant
! values from spectral tables, specifically SED (spectral energy 
! distrbution) tables for stellar particle sources, as functions of age 
! and metallicity, and UV-spectrum tables, as functions of redshift.
! The modules are:
! spectrum_integrator_module
!    For dealing with the integration of spectra in general
! SED_module
!    For dealing with the SED data
! UV_module
!    For dealing with the UV data
!_________________________________________________________________________

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Module for integrating a wavelength-dependent spectrum
!_________________________________________________________________________
!
MODULE spectrum_integrator_module
!_________________________________________________________________________
  use amr_parameters,only:dp
  implicit none

  PUBLIC integrateSpectrum, f1, fLambda, fdivLambda, fSig, fSigLambda,   &
         fSigdivLambda, trapz1

  PRIVATE   ! default

CONTAINS

!*************************************************************************
FUNCTION integrateSpectrum(X, Y, N, e0, e1, species, func, doPrint)

! Integrate spectral weighted function in energy interval [e0,e1]
! X      => Wavelengths [angstrom]
! Y      => Spectral luminosity per angstrom at wavelenghts [XX A-1] 
! N      => Length of X and Y
! e0,e1  => Integrated interval [ev]
! species=> ion species, used as an argument in fx
! func   => Function which is integrated (of X, Y, species)
!-------------------------------------------------------------------------
  use rt_parameters,only:c_cgs,eV_to_erg, hp
  use amr_commons,only:myid
  real(kind=8):: integrateSpectrum, X(N), Y(N), e0, e1
  integer :: N, species
  interface
     real(kind=8) function func(wavelength,intensity,species)
       use amr_parameters,only:dp
       real(kind=8)::wavelength,intensity
       integer::species
     end function func
  end interface!----------------------------------------------------------
  real(kind=8),dimension(:),allocatable:: xx, yy, f
  real(dp):: la0, la1
  integer :: i
  logical,optional::doPrint
!-------------------------------------------------------------------------
  integrateSpectrum=0.
  if(N .le. 2) RETURN
  ! Convert energy interval to wavelength interval
  la0 = X(1) ; la1 = X(N)
  if(e1.gt.0) la0 = max(la0, 1.d8 * hp * c_cgs / e1 / eV_to_erg)                
  if(e0.gt.0) la1 = min(la1, 1.d8 * hp * c_cgs / e0 / eV_to_erg)
  if(la0 .ge. la1) RETURN         
  ! If we get here, the [la0, la1] inverval is completely within X
  allocate(xx(N)) ; allocate(yy(N)) ; allocate(f(N))
  xx =  la0   ;   yy =  0.   ;   f = 0.
  i=2
  do while ( i.lt.N .and. X(i).le.la0 )
     i = i+1                      !              Below wavelength interval
  enddo                           !   X(i) is now the first entry .gt. la0 
  ! Interpolate to value at la0
  yy(i-1) = Y(i-1) + (xx(i-1)-X(i-1))*(Y(i)-Y(i-1))/(X(i)-X(i-1))
  f(i-1)  = func(xx(i-1), yy(i-1), species)
  do while ( i.lt.N .and. X(i).le.la1 )              ! Now within interval
     xx(i) = X(i) ; yy(i) = Y(i) ; f(i) = func(xx(i),yy(i),species)
     i = i+1
  enddo                          ! i=N or X(i) is the first entry .gt. la1
  xx(i:) = la1                   !             Interpolate to value at la1
  yy(i) = Y(i-1) + (xx(i)-X(i-1))*(Y(i)-Y(i-1))/(X(i)-X(i-1))
  f(i)  = func(xx(i),yy(i),species)

  !if(present(doPrint)) then
  !   if(doprint) then
  !      write(*,*) e0,e1,la0,la1,N
  !      write(*,*) '***************'
  !      do i=1,N
  !         write(*,*) xx(i),f(i),yy(i)
  !      end do
  !      stop
  !   endif
  !endif

  integrateSpectrum = trapz1(xx,f,i)
  deallocate(xx) ; deallocate(yy) ; deallocate(f)

END FUNCTION integrateSpectrum

!*************************************************************************
! FUNCTIONS FOR USE WITH integrateSpectrum:
! lambda  => wavelengths in Angstrom
! f       => function of wavelength (a spectrum in some units)
! species => 1=HI, 2=HeI or 3=HeII
!_________________________________________________________________________
FUNCTION f1(lambda, f, species)
  real(kind=8):: f1, lambda, f
  integer :: species
  f1 = f
END FUNCTION f1

FUNCTION fLambda(lambda, f, species)
  real(kind=8):: fLambda, lambda, f
  integer :: species
  fLambda = f * lambda
END FUNCTION fLambda

FUNCTION fdivLambda(lambda, f, species)
  real(kind=8):: fdivlambda, lambda, f
  integer :: species
  fdivLambda = f / lambda
END FUNCTION fdivLambda

FUNCTION fSig(lambda, f, species)
  real(kind=8):: fSig, lambda, f
  integer :: species
  fSig = f * getCrosssection_Hui(lambda,species)
END FUNCTION fSig

FUNCTION fSigLambda(lambda, f, species)
  real(kind=8):: fSigLambda, lambda, f
  integer :: species
  fSigLambda = f * lambda * getCrosssection_Hui(lambda,species)
END FUNCTION fSigLambda

FUNCTION fSigdivLambda(lambda, f, species)
  real(kind=8):: fSigdivLambda, lambda, f
  integer :: species
  fSigdivLambda = f / lambda * getCrosssection_Hui(lambda,species)
END FUNCTION fSigdivLambda
!_________________________________________________________________________

!*************************************************************************
FUNCTION trapz1(X,Y,N,cum)

! Integrates function Y(X) along the whole interval 1..N, using a very 
! simple staircase method and returns the result.
! Optionally, the culumative integral is returned in the cum argument.
!-------------------------------------------------------------------------
  integer :: N,i
  real(kind=8):: trapz1
  real(kind=8):: X(N),Y(N)
  real(kind=8),optional::cum(N)
  real(kind=8),allocatable::cumInt(:)
!-------------------------------------------------------------------------
  allocate(cumInt(N))
  cumInt(:)=0.d0
  if (N.le.1) RETURN
  do i=2,N
     cumInt(i)= cumInt(i-1) + abs(X(i)-X(i-1)) * (Y(i)+Y(i-1)) / 2.d0
  end do
  trapz1 = cumInt(N)
  if(present(cum)) cum=cumInt
  deallocate(cumInt)
END FUNCTION trapz1

!*************************************************************************
FUNCTION getCrosssection_Hui(lambda, species)

! Gives an atom-photon cross-section of given species at given wavelength, 
! as given by Hui and Gnedin (1997).
! lambda  => Wavelength in angstrom
! species => 1=HI, 2=HeI or 3=HeII
! returns :  photoionization cross-section in cm^2
!------------------------------------------------------------------------
  use rt_parameters,only:c_cgs, eV_to_erg, ionEvs, hp
  real(kind=8)      :: lambda, getCrosssection_Hui
  integer           :: species
  real(kind=8)      :: E, E0, cs0, P, ya, yw, y0, y1, x, y
!------------------------------------------------------------------------
  E = hp * c_cgs/(lambda*1.d-8) / ev_to_erg         ! photon energy in ev
  if ( E .lt. ionEvs(species) ) then            ! below ionization energy
     getCrosssection_Hui=0.
     RETURN
  endif
  select case (species)
  case (1) ! HI
     E0 = 4.298d-1 ; cs0 = 5.475d-14    ; P  = 2.963 
     ya = 32.88    ; yw  = 0            ; y0 = 0         ; y1 = 0
  case (2) ! HeI
     E0 = 1.361d1    ; cs0 = 9.492d-16  ; P  = 3.188 
     ya = 1.469      ; yw  = 2.039      ; y0 = 0.4434    ; y1 = 2.136
  case (3) ! HeII
     E0 = 1.720      ; cs0 = 1.369d-14  ; P  = 2.963 
     ya = 32.88      ; yw  = 0          ; y0 = 0         ; y1 = 0
  end select

  x = E/E0 - y0
  y = sqrt(x**2+y1**2)

  getCrosssection_Hui = &
       cs0 * ((x-1.)**2 + yw**2) * y**(0.5*P-5.5)/(1.+sqrt(y/ya))**P

END FUNCTION getCrosssection_Hui


END MODULE spectrum_integrator_module

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Module for Stellar Energy Distribution table. 
!_________________________________________________________________________
!
MODULE SED_module
!_________________________________________________________________________
  use amr_parameters,only:dp
  use rt_parameters,only:nGroups
  implicit none

  PUBLIC nSEDgroups                                                      &
      , init_SED_table, inp_SED_table, update_SED_group_props            &
      , update_star_RT_feedback, star_RT_feedback

  PRIVATE   ! default

  ! Light properties for different spectral energy distributions----------
  integer::nSEDgroups=nGroups         ! Default: All groups are SED groups
  integer::SED_nA, SED_nZ=8           ! Number of age bins and Z bins
  ! Age and z logarithmic intervals and lowest values:
  real(dp),parameter::SED_dlgA=0.02d0
  real(dp)::SED_dlgZ
  real(dp)::SED_lgA0, SED_lgZ0
  real(dp),allocatable,dimension(:)::SED_ages,SED_zeds![Gyr],[m_met/m_gas]
  ! SED_table: iAges, imetallicities, igroups, properties 
  !                                         (Lum, Lum-acc, egy, csn, cse).
  ! Lum is photons per sec per solar mass (eV per sec per solar mass in 
  ! the case of SED_isEgy=true). Lum-acc is accumulated lum.
  real(dp),allocatable,dimension(:,:,:,:)::SED_table
  ! ----------------------------------------------------------------------

CONTAINS

!*************************************************************************
SUBROUTINE init_SED_table()

! Initiate SED properties table, which gives photon luminosities, 
! integrated luminosities, average photon cross sections and energies of 
! each photon group as a function of stellar population age and
! metallicity.  The SED is read from a directory specified by sed_dir.
!-------------------------------------------------------------------------
  use amr_commons,only:myid
  use rt_parameters
  use spectrum_integrator_module
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  ! Temporary SSP/SED parameters (read from SED files):
  integer:: nAges, nzs, nLs              ! # of bins of age, z, wavelength
  real(kind=8),allocatable::ages(:), Zs(:), Ls(:), rebAges(:)
  real(kind=8),allocatable::SEDs(:,:,:)           ! SEDs f(lambda,age,met)
  real(kind=8),allocatable::tbl(:,:,:), tbl2(:,:,:), reb_tbl(:,:,:)
  integer::i,ia,iz,ip,ii,dum
  character(len=128)::fZs, fAges, fSEDs                        ! Filenames
  logical::ok,okAge,okZ
  real(kind=8)::dlgA, pL0, pL1, tmp, minAge
  integer::locid,ncpu,ierr
  integer::nv=3+2*nIons  ! # vars in SED table: L,Lacc,egy,nions*(csn,egy)
!-------------------------------------------------------------------------
  if(myid==1) &
        write(*,*) 'Stars are photon emitting, so initializing SED table'
  if(SED_DIR.eq.'')call get_environment_variable('RAMSES_SED_DIR',SED_DIR)
  inquire(FILE=TRIM(sed_dir)//'/all_seds.dat', exist=ok)
  if(.not. ok)then
     if(myid.eq.1) then 
        write(*,*)'Cannot access SED directory ',TRIM(sed_dir)
        write(*,*)'Directory '//TRIM(sed_dir)//' not found'
        write(*,*)'You need to set the RAMSES_SED_DIR envvar' // &
                  ' to the correct path, or use the namelist.'
     endif
     call clean_stop
  end if
  write(fZs,'(a,a)')   trim(sed_dir),"/metallicity_bins.dat"
  write(fAges,'(a,a)') trim(sed_dir),"/age_bins.dat"
  write(fSEDs,'(a,a)') trim(sed_dir),"/all_seds.dat"
  inquire(file=fZs, exist=okZ)
  inquire(file=fAges, exist=okAge)
  inquire(file=fSEDs, exist=ok)
  if(.not. ok .or. .not. okAge .or. .not. okZ) then
     if(myid.eq.1) then 
        write(*,*) 'Cannot read SED files...'
        write(*,*) 'Check if SED-directory contains the files ',  &
                   'metallicity_bins.dat, age_bins.dat, all_seds.dat'
     endif
     call clean_stop
  end if

  ! READ METALLICITY BINS-------------------------------------------------
  open(unit=10,file=fZs,status='old',form='formatted')
  read(10,'(i8)') nzs
  allocate(zs(nzs))
  do i = 1, nzs
     read(10,'(e14.6)') zs(i)
  end do
  close(10)
  ! READ AGE BINS---------------------------------------------------------
  open(unit=10,file=fAges,status='old',form='formatted')
  read(10,'(i8)') nAges
  allocate(ages(nAges))
  do i = 1, nAges
     read(10,'(e14.6)') ages(i)
  end do
  close(10)
  ages = ages*1.e-9                       !         Convert from yr to Gyr
  if(ages(1) .ne. 0.) ages(1) = 0.
  ! READ SEDS-------------------------------------------------------------
  open(unit=10,file=fSEDs,status='old',form='unformatted')
  read(10) nLs, dum
  allocate(Ls(nLs))
  read(10) Ls(:)
  allocate(SEDs(nLs,nAges,nzs))
  do iz = 1, nzs
     do ia = 1, nAges
        read(10) SEDs(:,ia,iz)
     end do
  end do
  close(10)
  ! If MPI then share the SED integration between the cpus:
#ifndef WITHOUTMPI
  call MPI_COMM_RANK(MPI_COMM_WORLD,locid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
#endif
#ifdef WITHOUTMPI
  locid=0
  ncpu=1
#endif

  ! Perform SED integration of luminosity, csn and egy per (age,Z) bin----
  allocate(tbl(nAges,nZs,nv))
  do ip = 1,nSEDgroups                                ! Loop photon groups
     tbl=0.
     pL0 = groupL0(ip) ; pL1 = groupL1(ip)! eV interval of photon group ip
     do iz = 1, nzs                                     ! Loop metallicity
     do ia = locid+1,nAges,ncpu                                 ! Loop age
        tbl(ia,iz,1) = getSEDLuminosity(Ls,SEDs(:,ia,iz),nLs,pL0,pL1)
        tbl(ia,iz,3) = getSEDEgy(Ls,SEDs(:,ia,iz),nLs,pL0,pL1)
        do ii = 1,nIons                                     ! Loop species
           tbl(ia,iz,2+ii*2) = getSEDcsn(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,ii)
           tbl(ia,iz,3+ii*2) = getSEDcse(Ls,SEDs(:,ia,iz),nLs,pL0,pL1,ii)
        end do ! End species loop
     end do ! End age loop
     end do ! End Z loop

#ifndef WITHOUTMPI
     allocate(tbl2(nAges,nzs,nv))
     call MPI_ALLREDUCE(tbl,tbl2,nAges*nzs*nv,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     tbl = tbl2
     deallocate(tbl2)
#endif     

     ! Integrate the cumulative luminosities:
     do iz = 1, nzs ! Loop metallicity
        tmp = trapz1( ages, tbl(:,iz,1), nAges, tbl(:,iz,2) )
     end do
     tbl(:,:,2) = tbl(:,:,2) * Gyr2sec ! Convert from sec-1 to Myr-1

     ! Now the SED properties are in tbl...just need to rebin it, to get !
     ! even log-intervals between bins for fast interpolation.           !
     dlgA = SED_dlgA ; SED_dlgZ = -SED_nz
     ages(1)=minAge     ! Bc zero age impossible for logarithmic rebinning
     call rebin_log(dlgA, dble(SED_dlgZ)                                 &
          , tbl(2:nAges,:,:), nAges-1, nZs, ages(2:nAges), zs, nv        &
          , reb_tbl, SED_nA, SED_nZ, rebAges, SED_Zeds)
     SED_nA=SED_nA+1                              ! Make room for zero age
     if(ip .eq. 1 ) allocate(SED_table(SED_nA, SED_nZ, nSEDgroups, nv))
     SED_table(1, :,ip,:) = reb_tbl(1,:,:)            ! Zero age properties
     SED_table(1, :,ip,2) = 0.                        !  Lacc=0 at zero age
     SED_table(2:,:,ip,:) = reb_tbl
     deallocate(reb_tbl)

  end do ! End photon group loop

  SED_lgZ0 = log10(SED_Zeds(1))                  ! Interpolation intervals
  SED_lgA0 = log10(rebAges(1))
  allocate(SED_ages(SED_nA))
  SED_ages(1)=0.d0 ; SED_ages(2:)=rebAges ;    ! Must have zero initial age

  deallocate(SEDs) ; deallocate(tbl)
  deallocate(ages) ; deallocate(rebAges) 
  deallocate(zs)
  deallocate(Ls)
 
  if (myid==1) call write_SEDtable

END SUBROUTINE init_SED_table

!*************************************************************************
SUBROUTINE update_SED_group_props()

! Compute and assign to all SED photon groups an average of the  
! quantities group_csn and group_egy from all the star particles in the
! simulation, weighted by the luminosity of the particles.
! If there are no stars, we assign the SED properties of zero-age,
! zero-metallicity stellar populations.
!-------------------------------------------------------------------------
  use amr_commons
  use pm_parameters
  use pm_commons
  use rt_parameters
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer :: i, ip, ii, nstars, info
  real(dp),save,allocatable,dimension(:)::  L_star
  real(dp),save,allocatable,dimension(:,:)::csn_star, cse_star
  real(dp),save,allocatable,dimension(:)::  egy_star
  real(dp),save,allocatable,dimension(:)::  sum_L_cpu,sum_L_all
  real(dp),save,allocatable,dimension(:,:)::sum_csn_cpu,sum_csn_all
  real(dp),save,allocatable,dimension(:,:)::sum_cse_cpu,sum_cse_all
  real(dp),save,allocatable,dimension(:)::sum_egy_cpu,sum_egy_all
  real(dp):: mass, age, Z
!-------------------------------------------------------------------------
  if(.not. allocated(L_star)) then
     allocate(L_star(nSEDgroups)) 
     allocate(egy_star(nSEDgroups)) 
     allocate(csn_star(nSEDgroups,nIons)) 
     allocate(cse_star(nSEDgroups,nIons)) 
     allocate(sum_L_cpu(nSEDgroups))    
     allocate(sum_L_all(nSEDgroups)) 
     allocate(sum_egy_cpu(nSEDgroups)) 
     allocate(sum_egy_all(nSEDgroups)) 
     allocate(sum_csn_cpu(nSEDgroups,nIons)) 
     allocate(sum_csn_all(nSEDgroups,nIons)) 
     allocate(sum_cse_cpu(nSEDgroups,nIons)) 
     allocate(sum_cse_all(nSEDgroups,nIons)) 
  endif
  sum_L_cpu   = 0d0 ! Accumulated luminosity, avg cross sections and 
  sum_egy_cpu = 0d0 ! photon energies for all stars belonging to  
  sum_csn_cpu = 0d0 ! 'this' cpu
  sum_cse_cpu = 0d0
  do i=1,npartmax
     if(levelp(i).le.0 .or. idp(i).eq.0 .or. tp(i).eq.0.)                &
        cycle ! not a star
     ! particle exists and is a star
     nstars = nstars + 1
     mass = mp(i)
     call getAgeGyr(tp(i), age)                         !     age = [Gyrs]
     if(metal) then
        Z = max(zp(i), 10.d-5)                          ! [m_metals/m_tot]
     else
        Z = max(z_ave*0.02, 10.d-5)                     ! [m_metals/m_tot]
     endif
     call inp_SED_table(age, Z, 1, .false., L_star)     !  [# s-1 m_sun-1]
     call inp_SED_table(age, Z, 3, .true., egy_star(:)) !             [ev]
     do ii=1,nIons
        call inp_SED_table(age, Z, 2+2*ii, .true., csn_star(:,ii))! [cm^2]
        call inp_SED_table(age, Z, 3+2*ii, .true., cse_star(:,ii))! [cm^2]
     end do

     do ip=1,nSEDgroups
        L_star(ip) = L_star(ip) * mass             !       [# photons s-1]       
        sum_L_cpu(ip)    =   sum_L_cpu(ip)   + L_star(ip)
        sum_egy_cpu(ip) =  sum_egy_cpu(ip)   + L_star(ip) * egy_star(ip)
        sum_csn_cpu(ip,:)= sum_csn_cpu(ip,:) + L_star(ip) * csn_star(ip,:)
        sum_cse_cpu(ip,:)= sum_cse_cpu(ip,:) + L_star(ip) * cse_star(ip,:)
     end do

  end do

  ! Sum up for all cpus
#ifdef WITHOUTMPI
  sum_L_all   = sum_L_cpu
  sum_egy_all = sum_egy_cpu
  sum_csn_all = sum_csn_cpu
  sum_cse_all = sum_cse_cpu
#else
  call MPI_ALLREDUCE(sum_L_cpu,   sum_L_all,   nSEDgroups,               &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(sum_egy_cpu, sum_egy_all, nSEDgroups,               &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(sum_csn_cpu, sum_csn_all, nSEDgroups*nIons,         &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(sum_cse_cpu, sum_cse_all, nSEDgroups*nIons,         &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#endif
  
  ! ...and take averages weighted by luminosities
  do ip=1,nSEDgroups
     ! No update for non-SED groups (L0>L1):
     if(groupL0(ip).ne. 0d0 .and. groupL1(ip) .ne. 0d0 .and. &
          (groupL0(ip) .ge. groupL1(ip)) ) cycle
     if(sum_L_all(ip) .gt. 0.) then
        group_egy(ip)   = sum_egy_all(ip)   / sum_L_all(ip)
        group_csn(ip,:) = sum_csn_all(ip,:) / sum_L_all(ip)
        group_cse(ip,:) = sum_cse_all(ip,:) / sum_L_all(ip)
     else ! no stars -> assign zero-age zero-metallicity props
        group_egy(ip)       = SED_table(1,1,ip,3)
        do ii=1,nIons
           group_csn(ip,ii) = SED_table(1,1,ip,2+2*ii)
           group_cse(ip,ii) = SED_table(1,1,ip,3+2*ii)
        enddo
     endif
  end do
  call updateRTgroups_coolConstants
  if(myid==1) write(*,*) &
                    'SED Photon groups updated through stellar polling'
  !call write_group_props(.true.,6)

END SUBROUTINE update_SED_group_props

!*************************************************************************
SUBROUTINE update_star_RT_feedback(ilevel)

! Turn on RT advection if needed.
! Update photon group properties from stellar populations.
!-------------------------------------------------------------------------
  use amr_parameters
  use amr_commons
  use rt_parameters
  use pm_commons
  integer::ilevel
  logical,save::groupProps_init=.false.
!-------------------------------------------------------------------------
  if(rt_star.and. nstar_tot .gt. 0) then
     if(.not.rt_advect) then ! Turn on RT advection due to newborn stars:
        if(myid==1) write(*,*) '*****************************************'
        if(myid==1) write(*,*) 'Stellar RT turned on at a=',aexp      
        if(myid==1) write(*,*) '*****************************************'
        rt_advect=.true.     
     endif
     ! Set group props from stellar populations:                 
     if(sedprops_update .gt. 0 .and. .not.groupProps_init) then
        call update_SED_group_props                                         
        groupProps_init=.true.        
     else if(sedprops_update .gt. 0 .and. ilevel==levelmin &
          .and. mod(nstep_coarse,sedprops_update)==0) then               
        call update_SED_group_props                         
     endif
  endif
END SUBROUTINE update_star_RT_feedback

!*************************************************************************
SUBROUTINE star_RT_feedback(ilevel, dt)

! This routine adds photons from radiating stars to appropriate cells in 
! the  hydro grid.
! ilegel =>  grid level in which to perform the feedback
! ti     =>  initial time for the timestep (code units)
! dt     =>  real timestep length in code units
!-------------------------------------------------------------------------
  use pm_commons
  use amr_commons
  use rt_parameters
  integer:: ilevel
  real(dp):: dt
  integer:: igrid, jgrid, ipart, jpart, next_part
  integer:: i, ig, ip, npart1, npart2, icpu
  integer,dimension(1:nvector),save:: ind_grid, ind_part, ind_grid_part
!-------------------------------------------------------------------------
  if(.not.rt_advect)RETURN
  if(nstar_tot .le. 0 ) return
  if(numbtot(1,ilevel)==0)return ! number of grids in the level
  if(verbose)write(*,111)ilevel
  ! Gather star particles only.
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel) ! grid index
     ig=0                     ! index of grid with stars (in ind_grid)
     ip=0                     ! index of star particle   (in ind_part)
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)   ! Number of particles in the grid
        npart2=0              ! number of selected (i.e. star) particles
        
        ! Count star particles in the grid
        if(npart1 > 0)then
          ipart = headp(igrid)        ! particle index
           ! Loop over particles
           do jpart = 1, npart1
              ! Save next particle       <--- Very important !!!
              next_part = nextp(ipart)
              !if(idp(ipart) .gt. 0 .and. tp(ipart) .ne. 0.d0) then 
              if(idp(ipart) .ne. 0 .and. tp(ipart) .ne. 0.d0) then 
                 npart2 = npart2+1     ! only stars
              endif
              ipart = next_part        ! Go to next particle
           end do
        endif

        ! Gather star particles within the grid
        if(npart2 > 0)then        
           ig = ig+1
           ind_grid(ig) = igrid
           ipart = headp(igrid)
           ! Loop over particles
           do jpart = 1, npart1
              ! Save next particle      <--- Very important !!!
              next_part = nextp(ipart)
              ! Select only star particles
              if(idp(ipart) .ne. 0 .and. tp(ipart) .ne. 0.d0) then
                 if(ig==0)then
                    ig=1      
                    ind_grid(ig)=igrid
                 end if
                 ip = ip+1
                 ind_part(ip) = ipart
                 ind_grid_part(ip) = ig ! points to grid a star is in  
              endif
              if(ip == nvector)then
                 if(static) then
                    call star_RT_vsweep_pp( &
                         ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel)
                 else
                    call star_RT_vsweep( &
                         ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel)
                 endif
                 ip = 0
                 ig = 0
              end if
              ipart = next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid = next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip > 0) then
        if(static) then
           call star_RT_vsweep_pp( &
                ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel)
        else
           call star_RT_vsweep( &
                ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel)
        endif
     endif
  end do 
  ! End loop over cpus

111 format('   Entering star_rt_feedback for level ',I2)

END SUBROUTINE star_RT_feedback

!*************************************************************************
!*************************************************************************
! START PRIVATE SUBROUTINES AND FUNCTIONS*********************************

!*************************************************************************
FUNCTION getSEDLuminosity(X, Y, N, e0, e1)

! Compute and return luminosity in energy interval (e0,e1) [eV]   
! in SED Y(X). Assumes X is in Angstroms and Y in Lo/Angstroms/Msun. 
! (Lo=[Lo_sun], Lo_sun=[erg s-1]. total solar luminosity is 
! Lo_sun=10^33.58 erg/s)
! returns: Photon luminosity in, [# s-1 Msun-1], 
!                             or [eV s-1 Msun-1] if SED_isEgy=true
!-------------------------------------------------------------------------
  use rt_parameters,only:c_cgs, hp, ev_to_erg, SED_isEgy
  use spectrum_integrator_module
  real(kind=8):: getSEDLuminosity, X(n), Y(n), e0, e1
  integer :: N, species
  real(kind=8),parameter :: const=1.0e-8/hp/c_cgs
  real(kind=8),parameter :: Lsun=3.8256d33 ! Solar luminosity [erg/s]
  ! const is a div by ph energy => ph count.  1e-8 is a conversion into 
  ! cgs, since wly=[angstrom] h=[erg s-1], c=[cm s-1]
!-------------------------------------------------------------------------
  species          = 1                   ! irrelevant but must be included
  if(.not. SED_isEgy) then               !  Photon number per sec per Msun
     getSEDLuminosity = const &
          * integrateSpectrum(X, Y, N, e0, e1, species, fLambda)
     getSEDLuminosity = getSEDLuminosity*Lsun  ! Scale by solar luminosity
  else                             ! SED_isEgy=true -> eV per sec per Msun
     getSEDLuminosity = integrateSpectrum(X, Y, N, e0, e1, species, f1)
     ! Scale by solar lum and convert to eV (bc group energies are in eV)
     getSEDLuminosity = getSEDLuminosity/eV_to_erg*Lsun 
  endif
END FUNCTION getSEDLuminosity

!*************************************************************************
FUNCTION getSEDEgy(X, Y, N, e0, e1)

! Compute average energy, in eV, in energy interval (e0,e1) [eV] in SED
! Y(X). Assumes X is in Angstroms and Y is energy weight per angstrom
! (not photon count). 
!-------------------------------------------------------------------------
  use rt_parameters,only:c_cgs, eV_to_erg, hp
  use spectrum_integrator_module
  real(dp):: getSEDEgy, X(N), Y(N), e0, e1, norm
  integer :: N,species
  real(dp),parameter :: const=1.d8*hp*c_cgs/eV_to_erg! energy conversion
!-------------------------------------------------------------------------
  species      = 1                       ! irrelevant but must be included
  norm         = integrateSpectrum(X, Y, N, e0, e1, species, fLambda)
  getSEDEgy    = const * &
                 integrateSpectrum(X, Y, N, e0, e1, species, f1) / norm
END FUNCTION getSEDEgy

!*************************************************************************
FUNCTION getSEDcsn(X, Y, N, e0, e1, species)

! Compute and return average photoionization 
! cross-section, in cm^2, for a given energy interval (e0,e1) [eV] in 
! SED Y(X). Assumes X is in Angstroms and that Y is energy weight per 
! angstrom (not photon #).
! Species is a code for the ion in question: 1=HI, 2=HeI, 3=HeIII
!-------------------------------------------------------------------------
  use spectrum_integrator_module
  use rt_parameters,only:ionEVs
  real(kind=8):: getSEDcsn, X(N), Y(N), e0, e1, norm
  integer :: N, species
!-------------------------------------------------------------------------
  if(e1 .gt. 0. .and. e1 .le. ionEvs(species)) then
     getSEDcsn=0. ; RETURN    ! [e0,e1] below ionization energy of species
  endif
  norm     = integrateSpectrum(X, Y, N, e0, e1, species, fLambda)
  getSEDcsn= integrateSpectrum(X, Y, N, e0, e1, species, fSigLambda)/norm
END FUNCTION getSEDcsn

!************************************************************************
FUNCTION getSEDcse(X, Y, N, e0, e1, species)

! Compute and return average energy weighted photoionization 
! cross-section, in cm^2, for a given energy interval (e0,e1) [eV] in 
! SED Y(X). Assumes X is in Angstroms and that Y is energy weight per 
! angstrom (not photon #).
! Species is a code for the ion in question: 1=HI, 2=HeI, 3=HeIII
!-------------------------------------------------------------------------
  use spectrum_integrator_module
  use rt_parameters,only:ionEVs
  real(dp):: getSEDcse, X(N), Y(N), e0, e1, norm
  integer :: N, species
!-------------------------------------------------------------------------
  if(e1 .gt. 0. .and. e1 .le. ionEvs(species)) then
     getSEDcse=0. ; RETURN    ! [e0,e1] below ionization energy of species
  endif
  norm      = integrateSpectrum(X, Y, N, e0, e1, species, f1)
  getSEDcse = integrateSpectrum(X, Y, N, e0, e1, species, fSig) / norm
END FUNCTION getSEDcse

!*************************************************************************
SUBROUTINE rebin_log(xint_log, yint_log,                                 &
               data,       nx,       ny,     x,     y,     nz,           & 
               new_data, new_nx, new_ny, new_x, new_y      )

! Rebin the given 2d data into constant logarithmic intervals, using 
! linear interpolation.
! xint_log,  => x and y intervals in the rebinned data. If negative,
! yint_log      these values represent the new number of bins.
! data       => The 2d data to be rebinned
! nx,ny      => Number of points in x and y in the original data
! nz         => Number of values in the data
! x,y        => x and y values for the data
! new_data  <=  The rebinned 2d data
! new_nx,ny <=  Number of points in x and y in the rebinned data
! new_x,y   <=  x and y point values for the rebinned data
!-------------------------------------------------------------------------
  real(kind=8):: xint_log, yint_log
  integer,intent(in):: nx, ny, nz
  integer::new_nx, new_ny
  real(kind=8),intent(in):: x(nx),y(ny)
  real(kind=8),intent(in):: data(nx,ny,nz)
  real(kind=8),dimension(:,:,:),allocatable:: new_data
  real(dp),dimension(:),allocatable:: new_x, new_lgx, new_y, new_lgy
  real(kind=8):: dx0, dx1, dy0, dy1, x_step, y_step
  real(kind=8):: x0lg, x1lg, y0lg, y1lg
  integer :: i, j, ix, iy, ix1, iy1
!-------------------------------------------------------------------------
  if(allocated(new_x)) deallocate(new_x)
  if(allocated(new_y)) deallocate(new_y)
  if(allocated(new_data)) deallocate(new_data)

  ! Find dimensions of the new_data and initialize it
  x0lg = log10(x(1));   x1lg = log10(x(nx)) 
  y0lg = log10(y(1));   y1lg = log10(y(ny))
  
  if(xint_log .lt. 0 .and. nx .gt. 1) then
     new_nx=-xint_log                             ! xint represents wanted
     xint_log = (x1lg-x0lg)/(new_nx-1)            !     number of new bins
  else
     new_nx = (x1lg-x0lg)/xint_log + 1
  endif
  allocate(new_x(new_nx)) ; allocate(new_lgx(new_nx))
  do i = 0, new_nx-1                              !  initialize the x-axis
     new_lgx(i+1) = x0lg + i*xint_log 
  end do
  new_x=10.d0**new_lgx

  if(yint_log .lt. 0 .and. ny .gt. 1) then        ! yint represents wanted
     new_ny=-yint_log                             !     number of new bins
     yint_log = (y1lg-y0lg)/(new_ny-1)
  else
     new_ny = (y1lg-y0lg)/yint_log + 1
  endif
  allocate(new_y(new_ny)) ; allocate(new_lgy(new_ny))
  do j = 0, new_ny-1                              !      ...and the y-axis
     new_lgy(j+1) = y0lg + j*yint_log
  end do
  new_y=10.d0**new_lgy

  ! Initialize new_data and find values for each point in it
  allocate(new_data(new_nx, new_ny, nz))
  do j = 1, new_ny  
     call locate(y, ny, new_y(j), iy) 
     ! y(iy) <= new_y(j) <= y(iy+1)
     ! iy is lower bound, so it can be zero but not larger than ny
     if(iy < 1) iy=1
     if (iy < ny) then
        iy1  = iy + 1
        y_step = y(iy1) - y(iy)
        dy0  = max(new_y(j) - y(iy),    0.0d0)  / y_step
        dy1  = min(y(iy1)   - new_y(j), y_step) / y_step
     else
        iy1  = iy
        dy0  = 0.0d0 ;  dy1  = 1.0d0
     end if

     do i = 1, new_nx
        call locate(x, nx, new_x(i), ix)
        if(ix < 1) ix=1
        if (ix < nx) then
           ix1  = ix+1
           x_step = x(ix1)-x(ix)
           dx0  = max(new_x(i) - x(ix),    0.0d0)  / x_step
           dx1  = min(x(ix1)   - new_x(i), x_step) / x_step
        else
           ix1  = ix
           dx0  = 0.0d0 ;  dx1  = 1.0d0
        end if

        if (abs(dx0+dx1-1.0d0) .gt. 1.0d-6 .or.                          &
                                         abs(dy0+dy1-1.0d0) > 1.0d-6) then
           write(*,*) 'Screwed up the rebin interpolation ... '
           write(*,*) dx0+dx1,dy0+dy1
           call clean_stop
        end if

        new_data(i,j,:) =                                                &
             dx0 * dy0 * data(ix1,iy1,:) + dx1 * dy0 * data(ix, iy1,:) + &
             dx0 * dy1 * data(ix1,iy, :) + dx1 * dy1 * data(ix, iy, :)
     end do
  end do

  deallocate(new_lgx)
  deallocate(new_lgy)

END SUBROUTINE rebin_log

!*************************************************************************
SUBROUTINE write_SEDtable()

! Write the SED properties to a file (this is just in debugging, to check
! if the SEDs are being read correctly).
!-------------------------------------------------------------------------
  character(len=128)::filename
  integer::ip, i, j
!-------------------------------------------------------------------------
  do ip=1,nSEDgroups
     write(filename,'(A, I1, A)') 'SEDtable', ip, '.list'
     open(10, file=filename, status='unknown')
     write(10,*) SED_nA, SED_nZ
     
     do j = 1,SED_nz
        do i = 1,SED_nA
           write(10,900)                                                 &
                SED_ages(i)        ,    SED_zeds(j)        ,             &
                SED_table(i,j,ip,1),    SED_table(i,j,ip,2),             &
                SED_table(i,j,ip,3),    SED_table(i,j,ip,4),             &
                SED_table(i,j,ip,5),    SED_table(i,j,ip,6),             &
                SED_table(i,j,ip,7),    SED_table(i,j,ip,8),             &
                SED_table(i,j,ip,9)
        end do
     end do
     close(10)
  end do
900 format (ES15.4, ES15.4, ES15.4, ES15.4, f15.4                        &
           ,ES15.4, ES15.4, ES15.4, ES15.4, ES15.4, ES15.4)

END SUBROUTINE write_SEDtable

!*************************************************************************
SUBROUTINE inp_SED_table(age, Z, nProp, same, ret)

! Compute SED property by interpolation from table.
! input/output:
! age   => Star population age [Gyrs]
! Z     => Star population metallicity [m_metals/m_tot] 
! nprop => Number of property to fetch 
!          1=log(photon # intensity [# Msun-1 s-1]),
!          2=log(cumulative photon # intensity [# Msun-1]),
!          3=avg_egy, 2+2*iIon=avg_csn, 3+2*iIon=avg_cse 
! same  => If true then assume same age and Z as used in last call. 
!          In this case the interpolation indexes can be recycled.
! ret   => The interpolated values of the sed property for every photon 
!          group
!-------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  real(dp), intent(in):: age, Z
  real(dp):: lgAge, lgZ
  integer:: nProp
  logical:: same
  real(dp),dimension(:):: ret
  integer,save:: ia, iz
  real(dp),save:: da, da0, da1, dz, dz0, dz1
!-------------------------------------------------------------------------
  ! ia, iz: lower indexes: 0<ia<sed_nA etc.
  ! da0, da1, dz0, dz1: proportional distances from edges: 
  ! 0<=da0<=1, 0<=da1<=1 etc.
  if(.not. same) then
     if(age.le.0d0) then
        lgAge=-4d0 
     else
        lgAge = log10(age)
     endif
     lgZ=log10(Z)
     ia = min(max(floor((lgAge-SED_lgA0)/SED_dlgA ) + 2, 1  ),  SED_nA-1 )
     da = SED_ages(ia+1)-SED_ages(ia)
     da0= min( max(   (age-SED_ages(ia)) /da,       0. ), 1.          )
     da1= min( max(  (SED_ages(ia+1)-age)/da,       0. ), 1.          )
     
     iz = min(max(floor((lgZ-SED_lgZ0)/SED_dlgZ ) + 1,   1  ),  SED_nZ-1 )
     dz = sed_Zeds(iz+1)-SED_Zeds(iz)
     dz0= min( max(   (Z-SED_zeds(iz)) /dz,         0. ),  1.         )
     dz1= min( max(  (SED_Zeds(iz+1)-Z)/dz,         0. ),  1.         )

     if (abs(da0+da1-1.0d0) > 1.0d-5 .or. abs(dz0+dz1-1.0d0) > 1.0d-5) then
        write(*,*) 'Screwed up the sed interpolation ... '
        write(*,*) da0+da1,dz0+dz1
        call clean_stop
     end if
  endif

  ret = da0 * dz0 * SED_table(ia+1, iz+1, :, nProp) + &
        da1 * dz0 * SED_table(ia,   iz+1, :, nProp) + &
        da0 * dz1 * SED_table(ia+1, iz,   :, nProp) + &
        da1 * dz1 * SED_table(ia,   iz,   :, nProp)

END SUBROUTINE inp_SED_table

!*************************************************************************
SUBROUTINE getNPhotonsEmitted(age1_Gyr, dt_Gyr, Z, ret)

! Compute number of photons emitted by a stellar particle per solar mass 
! over a timestep.
! input/output:
! age1_Gyr => Star population age [Gyrs] at timestep end
! dt_gyr   => Timestep length in Gyr
! Z        => Star population metallicity [m_metals/m_tot] 
! ret      => # of photons emitted per solar mass over timestep
!-------------------------------------------------------------------------
  use rt_parameters
  real(dp),intent(in):: age1_Gyr, dt_Gyr, Z
  real(dp),dimension(nSEDgroups):: ret, Lc0, Lc1, SEDegy
!-------------------------------------------------------------------------
  ! Lc0 = cumulative emitted photons at the start of the timestep
  call inp_SED_table(age1_Gyr-dt_Gyr, Z, 2, .false., Lc0)
  ! Lc1 = cumulative emitted photons at the end of the timestep
  call inp_SED_table(age1_Gyr, Z, 2, .false., Lc1)
  ret = max(Lc1-Lc0,0.)
  if(SED_isEgy) then ! Integrate correct energy rather than # of photons
     ! Divide emitted energy by group energy -> Photon count
     ret = ret / group_egy(1:nSEDgroups)
  endif
END SUBROUTINE getNPhotonsEmitted

!*************************************************************************
SUBROUTINE star_RT_vsweep(ind_grid,ind_part,ind_grid_part,ng,np,dt,ilevel)

! This routine is called by subroutine star_rt_feedback.
! Each star particle dumps a number of photons into the nearest grid cell
! using array rtunew.
! Radiation is injected into cells at level ilevel, but it is important
! to know that ilevel-1 cells may also get some radiation. This is due
! to star particles that have just crossed to a coarser level. 
!
! ind_grid      =>  grid indexes in amr_commons (1 to ng)
! ind_part      =>  star indexes in pm_commons(1 to np)
! ind_grid_part =>  points from star to grid (ind_grid) it resides in
! ng            =>  number of grids
! np            =>  number of stars
! dt            =>  timestep length in code units
! ilevel        =>  amr level at which we're adding radiation
!-------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  use rt_hydro_commons
  use rt_parameters
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  real(dp)::dt
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,ip
  real(dp)::dx,dx_loc,scale,vol_loc
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,ngroups),save::part_NpInp
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,3),save::id=0,igd=0,icd=0
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  ! units and temporary quantities
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, scale_Np  & 
            , scale_Fp, age, z, scale_inp, scale_Nphot, dt_Gyr           &
            , dt_loc_Gyr, scale_msun
  real(dp),parameter::vol_factor=2**ndim   ! Vol factor for ilevel-1 cells
!-------------------------------------------------------------------------
  if(.not. metal) z = log10(max(z_ave*0.02, 10.d-5))![log(m_metals/m_tot)]
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_Np, scale_Fp)
  dt_Gyr = dt*scale_t*sec2Gyr
  ! Mesh spacing in ilevel
  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max - icoarse_min + 1)
  skip_loc = (/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1) = dble(icoarse_min)
  if(ndim>1)skip_loc(2) = dble(jcoarse_min)
  if(ndim>2)skip_loc(3) = dble(kcoarse_min)
  scale = boxlen/dble(nx_loc) ! usually scale == 1
  dx_loc = dx*scale
  vol_loc = dx_loc**ndim
  scale_inp = rt_esc_frac * scale_d / scale_np / vol_loc / m_sun    
  scale_nPhot = vol_loc * scale_np * scale_l**ndim / 1.d50
  scale_msun = scale_d * scale_l**ndim / m_sun    

  ! Lower left corners of 3x3x3 grid-cubes (with given grid in center)
  do idim = 1, ndim
     do i = 1, ng
        x0(i,idim) = xg(ind_grid(i),idim) - 3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i) = father(ind_grid(i))
  end do
  call get3cubefather(&
          ind_cell, nbors_father_cells, nbors_father_grids, ng, ilevel)
  ! now nbors_father cells are a cube of 27 cells with ind_cell in the 
  ! middle and nbors_father_grids are the 8 grids that contain these 27
  ! cells (though only one of those is fully included in the cube)

  ! Rescale position of stars to positions within 3x3x3 cell supercube
  do idim = 1, ndim
     do j = 1, np
        x(j,idim) = xp(ind_part(j),idim)/scale + skip_loc(idim)
        x(j,idim) = x(j,idim) - x0(ind_grid_part(j),idim)
        x(j,idim) = x(j,idim)/dx 
        ! so 0<x<2 is bottom cell, ...., 4<x<6 is top cell 
     end do
  end do

   ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim) = x(j,idim) ! So id=0-5 is the cell (in the 
     end do                    ! 3x3x3 supercube) containing the star
  end do

  ! Compute parent grids
  do idim = 1, ndim
     do j = 1, np
        igd(j,idim) = id(j,idim)/2 ! must be 0, 1 or 2
     end do
  end do
  do j = 1, np
     kg(j) = 1 + igd(j,1) + 3*igd(j,2) + 9*igd(j,3) ! 1 to 27
  end do
  do j = 1, np
     igrid(j) = son(nbors_father_cells(ind_grid_part(j),kg(j))) 
     ! grid (not cell) containing the star 
  end do

  ! Check if particles are entirely in level ilevel.
  ! This should always be the case in post-processing
  ok(1:np) = .true.
  do j = 1, np
     ok(j) = ok(j) .and. igrid(j) > 0
  end do ! if ok(j) is true then particle j's cell contains a grid. 
  ! Otherwise it is a leaf, and we need to fill it with radiation.

  ! Compute parent cell position within it's grid
  do idim = 1, ndim
     do j = 1, np
        if( ok(j) ) then
           icd(j,idim) = id(j,idim) - 2*igd(j,idim) ! 0 or 1
        end if
     end do
  end do
  do j = 1, np
     if( ok(j) ) then
        icell(j) = 1 + icd(j,1) + 2*icd(j,2) + 4*icd(j,3) ! 1 to 8
     end if
  end do

  ! Compute parent cell adress and particle radiation contribution
  do j = 1, np
     if(metal)z= max(zp(ind_part(j)), 10.d-5)      !      [m_metals/m_tot]
     call getAgeGyr(tp(ind_part(j)), age)          !  End-of-dt age [Gyrs]
     ! Possibilities:     Born i) before dt, ii) within dt, iii) after dt:
     dt_loc_Gyr = max(min(dt_Gyr, age), 0.)
     call getNPhotonsEmitted(age,dt_loc_Gyr,z,part_NpInp(j,1:nSEDgroups))
     part_NpInp(j,:) = part_NpInp(j,:)*mp(ind_part(j))*scale_inp !#photons

     if(showSEDstats .and. nSEDgroups .gt. 0) then
        step_nPhot = step_nPhot+part_NpInp(j,1)*scale_nPhot
        step_nStar = step_nStar+dt_loc_Gyr*Gyr2sec/scale_t
        step_mStar = step_mStar+mp(ind_part(j)) * scale_msun             &
                                          * dt_loc_Gyr * Gyr2sec / scale_t
     endif

     if( ok(j) )then
        indp(j) = ncoarse + (icell(j)-1)*ngridmax + igrid(j)
     else
        indp(j) = nbors_father_cells(ind_grid_part(j),kg(j))
     end if
  end do
  ! Increase photon density in cell due to stellar radiation
  do j=1,np
     if( ok(j) ) then                                      !   ilevel cell
        do ip=1,nSEDgroups
           rtunew(indp(j),iGroups(ip)) &
                = rtunew(indp(j),iGroups(ip))+part_NpInp(j,ip)
        end do
     else                                                  ! ilevel-1 cell
        do ip=1,nSEDgroups
           rtunew(indp(j),iGroups(ip)) = rtunew(indp(j),iGroups(ip))     &
                + part_NpInp(j,ip) / vol_factor
        end do
     endif
     !begin debug
     !     if(rtunew(indp(j),iGroups(1))*scale_np*rt_c_cgs .gt.1.d11) then
     !        write(*,777),ilevel, myid, rtunew(indp(j),                   &
     !           iGroups(1))*scale_np*rt_c_cgs, &
     !             irad(j,1), irad(j,2),                                 &
     !             mp(ind_part(j)),                                      &
     !             mp(ind_part(j))*scale_d*scale_l**3/m_sun,             &
     !             vol_loc*scale_l**3,                                   &
     !             dt, scale_t,                                          &
     !             ok(j)
     !        print*,'*******************************'
     !     endif
     !end debug
     !777 format('HERE:',I2, ' ', I2, 8(1pe11.3), L, ' AHA')
  end do

END SUBROUTINE star_RT_vsweep

!*************************************************************************
SUBROUTINE star_RT_vsweep_pp(ind_grid, ind_part, ind_grid_part, ng, np,  &
                             dt, ilevel)
! This routine is called by subroutine star_rt_feedback,
! in the case of rt-post-processing. The difference from doing coupled 
! RT is that in the post-processing case, the pointer from a grid to a 
! particle is always true, so we can trust completely that 
! ind_grid(ind_grid_part(j)) refers to the grid that particle ind_part(j)
! resides in, and we don't need to bother searching for the particle's 
! grid.  Each star particle dumps a number of photons into the nearest 
! grid cell using array rtunew.
!
! ind_grid      =>  grid indexes in amr_commons (1 to ng)
! ind_part      =>  star indexes in pm_commons(1 to np)
! ind_grid_part =>  points from star to grid (ind_grid) it resides in
! ng            =>  number of grids
! np            =>  number of stars
! dt            =>  timestep length in code units
! ilevel        =>  amr level at which we're adding radiation
!-------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  use rt_hydro_commons
  use rt_parameters
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  real(dp)::dt
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,ip
  real(dp)::dx,dx_loc,scale,vol_loc
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  ! Particle based arrays
  real(dp),dimension(1:nvector,nGroups),save::part_NpInp
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,3),save::id=0
  !integer ,dimension(1:nvector,1:ndim),save::id
  integer ,dimension(1:nvector),save::icell,indp
  real(dp),dimension(1:3)::skip_loc
  ! units and temporary quantities
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, scale_Np  & 
            , scale_Fp, age, z, scale_inp, scale_Nphot, dt_Gyr           &
            , dt_loc_Gyr, scale_msun
!-------------------------------------------------------------------------
  if(.not. metal) z = log10(max(z_ave*0.02, 10.d-5))![log(m_metals/m_tot)]
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_np, scale_Fp)
  dt_Gyr = dt*scale_t*sec2Gyr
  ! Mesh spacing in ilevel
  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max - icoarse_min + 1)
  skip_loc = (/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1) = dble(icoarse_min)
  if(ndim>1)skip_loc(2) = dble(jcoarse_min)
  if(ndim>2)skip_loc(3) = dble(kcoarse_min)
  scale = boxlen/dble(nx_loc) ! usually scale == 1
  dx_loc = dx*scale
  vol_loc = dx_loc**ndim
  scale_inp = rt_esc_frac * scale_d / scale_np / vol_loc / m_sun    
  scale_nPhot = vol_loc * scale_np * scale_l**ndim / 1.d50
  scale_msun = scale_d * scale_l**ndim / m_sun    

  ! Lower left corners of parent grid
  do idim = 1, ndim
     do i = 1, ng
        x0(i,idim) = xg(ind_grid(i),idim) - dx
     end do
  end do

  ! Rescale position of star to position within 3x3x3 cell supercube
  do idim = 1, ndim
     do j = 1, np
        x(j,idim) = xp(ind_part(j),idim)/scale + skip_loc(idim)
        x(j,idim) = x(j,idim) - x0(ind_grid_part(j),idim)
        x(j,idim) = x(j,idim)/dx 
        ! so 0<x<1 is bottom cell, ...., 1<x<2 is top cell
     end do
  end do

  ! Compute position of parent cell within parent grid
  do idim=1,ndim
     do j=1,np
        id(j,idim) = x(j,idim) ! id=0-1 is cell in parent containing star
     end do
  end do
  do j = 1, np
     icell(j) = 1 + id(j,1) + 2*id(j,2) + 4*id(j,3) ! 1 to 8
  end do

  ! Compute parent cell adress and particle effective mass
  do j = 1, np
     if(metal) z = max(zp(ind_part(j)), 10.d-5)   !       [m_metals/m_tot]
     call getAgeGyr(tp(ind_part(j)), age)         !   End-of-dt age [Gyrs]
     !Possibilities:      Born i) before dt, ii) within dt, iii) after dt:
     dt_loc_Gyr = max(min(dt_Gyr, age), 0.)
     call getNPhotonsEmitted(age,dt_loc_Gyr, z,part_NpInp(j,1:nSEDGroups))
     part_NpInp(j,:) = part_NpInp(j,:)*mp(ind_part(j))*scale_inp !#photons

     if(showSEDstats .and. nSEDgroups .gt. 0) then
        step_nPhot = step_nPhot+part_NpInp(j,1)*scale_nPhot
        step_nStar = step_nStar+dt_loc_Gyr*Gyr2sec/scale_t
        step_mStar = step_mStar+mp(ind_part(j)) * scale_msun             &
                                          * dt_loc_Gyr * Gyr2sec / scale_t
     endif

     indp(j)= ncoarse + (icell(j)-1)*ngridmax + ind_grid(ind_grid_part(j))
  end do

  do j=1,np                       ! Update hydro variables due to feedback
     do ip=1,nSEDgroups
        rtunew(indp(j),iGroups(ip)) &
             = rtunew(indp(j),iGroups(ip)) + part_NpInp(j,ip)
                
     end do
  end do

END SUBROUTINE star_RT_vsweep_pp

END MODULE SED_module


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Module for UV table of redshift dependent photon fluxes, cross sections
! and photon energies, per photon group. This is mainly useful for 
! calculating heating rates due to the homogeneous UV background in the 
! cooling-module, and for  radiative transfer of the UV background, where 
! we then assign intensities, average cross sections and energies to UV 
! photon groups.
! 
! There are two tables of UV properties:
! UV_rates_table: Integrated ionization and heating rates for a 
!                 homogeneous UV background.
! UV_groups_table: UV properties integrated over each photon group 
!                 wavelength interval: photon flux, average energy, 
!                 average cross section, energy weighted cross section.
!_________________________________________________________________________
!
MODULE UV_module
!_________________________________________________________________________
  use amr_parameters,only:dp
  use spectrum_integrator_module
  use rt_parameters,only:c_cgs, eV_to_erg, hp, nIons, ionEVs, nGroups

  implicit none

  PUBLIC nUVgroups, iUVvars_cool, UV_Nphot_cgs, init_UV_background       &
       , inp_UV_rates_table, inp_UV_groups_table, UV_minz, UV_maxz       &
       , update_UVsrc

  PRIVATE   ! default

  ! UV properties for different redshifts---------------------------------
  integer::UV_nz
  real(dp),allocatable,dimension(:)::UV_zeds
  ! Redshift interval
  real(dp)::UV_minz, UV_maxz  
  ! Table of integrated ionization and heating rates per ion species
  ! (UV_nz, nIons, 2):(z, nIon, ionization rate [# s-1]: Hrate [erg s-1]) 
  real(dp),allocatable,dimension(:,:,:)::UV_rates_table
  ! rt_n_UVsrc vectors of redshift dependent fluxes for UV background
  real(dp),allocatable::UV_fluxes_cgs(:)    !                    [#/cm2/s]
  real(dp),allocatable::UV_Nphot_cgs(:)     !        Photon density [cm-3]
  integer::nUVgroups=0                      !        # of UV photon groups 
  ! UV group indexes among nGroup groups, 
  ! UV Np indexes among solve_cooling vars:
  integer,allocatable::iUVgroups(:),iUVvars_cool(:) 
  ! Table of photon group props (UV_nz, nUVgroups, 1+2*nIons): 
  !(z, group, flux [#/cm2/s]: nIons*csn [cm-2]: nIons*egy [ev])
  real(dp),allocatable,dimension(:,:,:)::UV_groups_table
  ! The tables do not have constant redshift intervals, so we need to
  ! first locate the correct interval when doing interpolation.
  real(dp),parameter::fourPi=12.566371
  ! ----------------------------------------------------------------------

CONTAINS

!*************************************************************************
SUBROUTINE init_UV_background()
! Initiate UV props table, which gives photon fluxes and average photon 
! cross sections and energies of each photon group as a function of z.
! The data is read from a file specified by uv_file containing: 
! a) number of redshifts
! b) redshifts (increasing order)
! c) number of wavelengths
! d) wavelengths (increasing order) [Angstrom]
! e) fluxes per (redshift,wavelength) [photons cm-2 s-1 A-1 sr-1]
!-------------------------------------------------------------------------
  use amr_commons,only:myid
  use rt_parameters
  use SED_module
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer:: nLs                                 ! # of bins of wavelength
  real(kind=8),allocatable  :: Ls(:)            ! Wavelengths
  real(kind=8),allocatable  :: UV(:,:)          ! UV f(lambda,z)
  real(kind=8),allocatable  :: tbl(:,:), tbl2(:,:)
  integer::i,ia,iz,ip,ii,dum,locid,ncpu,ierr
  logical::ok
  real(kind=8)::da, dz, pL0,pL1
!-------------------------------------------------------------------------
  ! First check if there is any need for UV setup:
  if(rt_UVsrc_nHmax .le. 0d0 .and. .not. rt_UV_hom) return

  if(myid==1) print*,'Initializing UV background'

  ! Read UV spectra from file---------------------------------------------
  if(UV_FILE .eq. '') &
       call get_environment_variable('RAMSES_UV_FILE', UV_FILE)
  inquire(file=TRIM(uv_file), exist=ok)
  if(.not. ok)then
     if(myid.eq.1) then 
        write(*,*)'Cannot access UV file ',TRIM(uv_file)
        write(*,*)'File '//TRIM(uv_file)//' not found'
        write(*,*)'You need to set the RAMSES_UV_FILE envvar' // &
                  ' to the correct path, or use the namelist var uv_file'
     endif
     call clean_stop
  end if
  ! Read redshifts, wavelengths and spectra
  open(unit=10,file=TRIM(uv_file),status='old',form='unformatted')
  read(10) UV_nz ; allocate(UV_zeds(UV_nz)) ; read(10) UV_zeds(:)
  read(10) nLs   ; allocate(Ls(nLs))        ; read(10) Ls(:)
  allocate(UV(nLs,UV_nz))                   ; read(10) UV(:,:)
  close(10)
  ! If mpi then share the UV integration between the cpus:
#ifndef WITHOUTMPI
  call MPI_COMM_RANK(MPI_COMM_WORLD,locid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
#endif
#ifdef WITHOUTMPI
  locid=0 ; ncpu=1
#endif

  ! Shift the highest z in the table (10) to reionization epoch,
  ! so that we start injecting at z_reion
  if(z_reion .gt. UV_zeds(UV_nz-1))  UV_zeds(UV_nz) = z_reion
  UV_minz = UV_zeds(1) ; UV_maxz=UV_zeds(UV_nz)
  
  ! Non-propagated UV background -----------------------------------------
  if(rt_UV_hom) then
     if(myid==1) print*,'The UV background is homogeneous'
     allocate(UV_rates_table(UV_nz, nIons, 2))
     allocate(tbl(UV_nz, 2))
     do ii = 1, nIons
        tbl=0.
        do iz = locid+1,UV_nz,ncpu
           tbl(iz,1)= getUV_Irate(Ls,UV(:,iz),nLs,ii)
           tbl(iz,2)= getUV_Hrate(Ls,UV(:,iz),nLs,ii)
        end do
#ifndef WITHOUTMPI
        allocate(tbl2(UV_nz,2))
        call MPI_ALLREDUCE(tbl, tbl2, UV_nz*2,  &
             MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        tbl = tbl2
        deallocate(tbl2)
#endif     
        UV_rates_table(:,ii,:)=tbl 
     end do
     deallocate(tbl)
     if (myid==1) call write_UVrates_table
  endif

  ! Note: Need to take special care of the highest redshift. Some UV     !
  !       spectra (read: Faucher-Giguere) have zero fluxes at the highest!
  !       redshift, and integrations of energies and cross sections of   !
  !       these fluxes results in a division by zero and NaNs.           !
  !       The approach taken is then to keep zero flux at the highest    !
  !       redshift but otherwise use the same values of energies and     !
  !       cross sections as the next-highest redshift.                   !

  ! Propagated UV background----------------------------------------------
  if(rt_UVsrc_nHmax .gt. 0.d0) then ! UV propagation from diffuse cells--
     if(myid==1) print*,'The UV background is propagated'
     if(myid==1 .and. rt_UV_hom) then
          print*,'ATT: UV background is BOTH homogeneous and propagated'
          print*,'  You likely don''t want this duplicated background...'
       endif
     rt_isDiffuseUVsrc=.true.
     if(nUVgroups .eq. 0) nUVgroups = nGroups ! Default: All groups are UV
     nSEDgroups=nGroups-nUVgroups
     ! SED groups are the first nSEDgroups, UV groups the last nUVgroups
     allocate(iUVgroups(nUVgroups)) ; allocate(iUVvars_cool(nUVgroups))
     do i=1,nUVgroups                 !      Initialize UV group indexes
        iUVgroups(i) = nGroups-nUVgroups+i        ! UV groups among groups
        iUVvars_cool(i) = 4+iUVgroups(i) !UV Np's in vars in solve_cooling
     end do
     if(nUVgroups .gt. 0) then
        allocate(UV_fluxes_cgs(nUVgroups))             ; UV_fluxes_cgs=0.
        allocate(UV_Nphot_cgs(nUVgroups))              ; UV_Nphot_cgs=0.
     endif

     ! Initialize photon groups table-------------------------------------
     allocate(UV_groups_table(UV_nz, nUVgroups, 2+2*nIons))                   
     allocate(tbl(UV_nz, 2+2*nIons))                                      
     do ip = 1,nUVgroups           !                    Loop photon groups
        tbl=0.
        pL0 = groupL0(nSEDgroups+ip) !  Energy interval of photon group ip
        pL1 = groupL1(nSEDgroups+ip) !
        do iz = locid+1,UV_nz,ncpu
           tbl(iz,1) =        getUVFlux(Ls,UV(:,iz),nLs,pL0,pL1)
           if(tbl(iz,1) .eq. 0.d0) cycle     ! Can't integrate zero fluxes
           tbl(iz,2) =        getUVEgy(Ls,UV(:,iz),nLs,pL0,pL1)
           do ii = 1,nIons
              tbl(iz,1+ii*2)= getUVcsn( Ls,UV(:,iz),nLs,pL0,pL1,ii)
              tbl(iz,2+ii*2)= getUVcse( Ls,UV(:,iz),nLs,pL0,pL1,ii)
           end do
        end do
#ifndef WITHOUTMPI
        allocate(tbl2(UV_nz,1+2*nIons))
        call MPI_ALLREDUCE(tbl,tbl2,UV_nz*(1+2*nIons),&
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        tbl = tbl2
        deallocate(tbl2)
#endif     
        if(tbl(UV_nz,1) .eq. 0.d0) &            !                Zero flux
             tbl(UV_nz,2:)=tbl(UV_nz-1,2:)
        UV_groups_table(:,ip,:)=tbl 
     end do
     deallocate(tbl) ; deallocate(Ls) ; deallocate(UV)

     call update_UVsrc
     if (myid==1) call write_UVgroups_tables
  endif ! End propagated UV background

END SUBROUTINE init_UV_background

!*************************************************************************
SUBROUTINE inp_UV_rates_table(z, ret)
! Compute UV heating and ionization rates by interpolation from table.
! z     => Redshift
! ret  <=  [nIons,2] interpolated values of all UV rates.
!          1=ionization rate [# s-1], 2=heating rate [erg s-1]
!-------------------------------------------------------------------------
  real(dp), intent(in):: z
  real(dp):: ret(nIons,2), dz0, dz1
  integer:: iz0, iz1
!-------------------------------------------------------------------------
  ret=0. ; if(z .gt. UV_maxz) RETURN
  call inp_1d(UV_zeds, UV_nz, z, iz0, iz1, dz0, dz1)
  ret = dz0*UV_rates_table(iz1, :, :) + dz1*UV_rates_table(iz0, :, :)
END SUBROUTINE inp_UV_rates_table

!*************************************************************************
SUBROUTINE inp_UV_groups_table(z, ret)
! Compute UV properties by interpolation from table.
! z    => Redshift
! ret <=  [nGroups,1+2*nIons) interpolated values of all UV properties.
!         1=ph. flux [#/cm2/s], 2*i=group_csn[cm-2], 1+2*i=group_egy [ev]
!-------------------------------------------------------------------------
  real(dp), intent(in):: z
  real(dp):: ret(nUVgroups,2+2*nIons), dz0, dz1
  integer:: iz0, iz1
!-------------------------------------------------------------------------
  ret=0. ; if(z .gt. UV_maxz) RETURN
  call inp_1d(UV_zeds, UV_nz, z, iz0, iz1, dz0, dz1)
  ret = dz0 * UV_groups_table(iz1, :, :) + dz1 * UV_groups_table(iz0, :,:)
END SUBROUTINE inp_UV_groups_table

!*************************************************************************
SUBROUTINE update_UVsrc

! Update UV background source properties. So as not to do too much of
! this, this should only be done every coarse timestep.
!-------------------------------------------------------------------------
  use rt_parameters
  use SED_module
  use amr_commons,only:t,levelmin,myid,aexp
  implicit none
  integer::ilevel,i
  real(dp),allocatable,save::UVprops(:,:) !Each group: flux, egy, csn, cse
  real(dp)::scale_Np, scale_Fp, redshift
!-------------------------------------------------------------------------
  if(.not.rt_isDiffuseUVsrc) return
  if(nUVgroups.le.0) then
     if(myid==1) write(*,*) 'No groups dedicated to the UV background!'
     RETURN
  endif
  if(.not. allocated(UVprops)) allocate(UVprops(nUVgroups,2+2*nIons))
  call rt_units(scale_Np, scale_Fp)

  redshift=1./aexp-1.

  ! Turn on RT after z=UV_maxz:
  if(redshift.le.UV_maxz .and. .not. rt_advect) then
     if(myid==1) then
        write(*,*) '*****************************************************'
        write(*,*) 'Turned on RT advection and the UV background'
        write(*,*) '*****************************************************'
     endif
     rt_advect=.true.
  endif

  if(redshift .gt. UV_maxz) return ! UV background not turned on yet

  call inp_UV_groups_table(redshift, UVprops)
  UV_fluxes_cgs(:)      = UVprops(:,1)
  UV_Nphot_cgs          = UV_fluxes_CGS/rt_c_cgs
  group_egy(iUVgroups)  = UVprops(:,2)
  do i=1,nIons
     group_csn(iUVgroups,i)  = UVprops(:,1+2*i)
     group_cse(iUVgroups,i)  = UVprops(:,2+2*i)
  enddo

  call updateRTgroups_CoolConstants     

  if(myid==1) then
     write(*,*) 'Updated UV fluxes [# cm-2 s-1] to'
     write(*,900) UV_fluxes_cgs
     call write_group_props(.true.,6)
  endif

900 format (20f16.6)
901 format (20(1pe16.6))
END SUBROUTINE update_UVsrc

!*************************************************************************
! START PRIVATE SUBROUTINES AND FUNCTIONS*********************************

!*************************************************************************
FUNCTION getUV_Irate(X, Y, N, species)
! Compute and return photoionization rate per ion, given the UV background
! Y(X). Assumes X is in Angstroms and Y in [# cm-2 s-1 sr-1 A-1].
! returns: Photoionization rate in # s-1
!-------------------------------------------------------------------------
  real(kind=8):: getUV_Irate, X(N), Y(N)
  integer :: N, species
!-------------------------------------------------------------------------
  getUV_Irate = fourPi *  & 
           integrateSpectrum(X,Y,N,dble(ionEvs(species)),0d0,species,fsig)
END FUNCTION getUV_Irate

!*************************************************************************
FUNCTION getUV_Hrate(X, Y, N, species)
! Compute and return heating rate per ion, given the UV background   
! Y(X). Assumes X is in Angstroms and Y in [# cm-2 s-1 sr-1 A-1].
! returns: Heating rate in erg s-1
!-------------------------------------------------------------------------
  real(kind=8):: getUV_Hrate, X(N), Y(N), e0
  integer :: N, species
  real(kind=8),parameter :: const1=fourPi*1.d8*hp*c_cgs
  real(kind=8),parameter :: const2=fourPi*eV_to_erg
!-------------------------------------------------------------------------
  e0=ionEvs(species)
  getUV_Hrate = &
        const1*integrateSpectrum(X,Y,N, e0, 0.d0, species, fsigDivLambda)&
       -const2*ionEvs(species) *                                         &
               integrateSpectrum(X,Y,N, e0, 0.d0, species, fsig) 
END FUNCTION getUV_Hrate

!*************************************************************************
FUNCTION getUVFlux(X, Y, N, e0, e1)
! Compute and return UV photon flux in energy interval (e0,e1) [eV]   
! in UV spectrum Y(X). Assumes X is in [A] and Y in # cm-2 s-1 sr-1 A-1.
! returns: Photon flux in # cm-2 s-1
!-------------------------------------------------------------------------
  real(kind=8):: getUVflux, X(N), Y(N), e0, e1
  integer :: N, species
!-------------------------------------------------------------------------
  species          = 1                   ! irrelevant but must be included
  getUVflux = fourPi*integrateSpectrum(X, Y, N, e0, e1, species, f1)
END FUNCTION getUVflux

!*************************************************************************
FUNCTION getUVEgy(X, Y, N, e0, e1)
! Compute average photon energy, in eV, in energy interval (e0,e1) [eV] in
! UV spectrum Y(X). Assumes X is in [A] and Y is # cm-2 s-1 sr-1 A-1.
!-------------------------------------------------------------------------
  real(dp):: getUVEgy, X(N), Y(N), e0, e1, norm
  integer :: N,species
  real(dp),parameter :: const=1.d8*hp*c_cgs/eV_to_erg    ! unit conversion
!-------------------------------------------------------------------------
  species      = 1                       ! irrelevant but must be included
  norm         = integrateSpectrum(X, Y, N, e0, e1, species, f1)
  getUVEgy  = const * &
            integrateSpectrum(X, Y, N, e0, e1, species, fdivLambda) / norm
END FUNCTION getUVEgy

!*************************************************************************
FUNCTION getUVcsn(X, Y, N, e0, e1, species)
! Compute and return average photoionization cross-section [cm2] for given 
! energy interval (e0,e1) [eV] in UV spectrum Y. Assumes X is in Angstroms
! and and Y is # cm-2 s-1 sr-1 A-1.
! Species is a code for the ion in question: 1=HI, 2=HeI, 3=HeIII
!-------------------------------------------------------------------------
  real(kind=8):: getUVcsn, X(N), Y(N), e0, e1, norm
  integer :: N, species
!-------------------------------------------------------------------------
  if(e1 .gt. 0. .and. e1 .le. ionEvs(species)) then
     getUVcsn=0. ; RETURN    ! [e0,e1] below ionization energy of species
  endif
  norm     = integrateSpectrum(X, Y, N, e0, e1, species, f1)
  getUVcsn = integrateSpectrum(X, Y, N, e0, e1, species, fSig)/norm
END FUNCTION getUVcsn

!************************************************************************
FUNCTION getUVcse(X, Y, N, e0, e1, species)
! Compute average energy weighted photoionization cross-section [cm2] for
! given energy interval (e0,e1) [eV] in UV spectrum Y. Assumes X is in 
! Angstroms and that Y is energy intensity per angstrom.
! Species is a code for the ion in question: 1=HI, 2=HeI, 3=HeIII
!-------------------------------------------------------------------------
  real(dp):: getUVcse, X(N), Y(N), e0, e1, norm
  integer :: N, species
!-------------------------------------------------------------------------
  if(e1 .gt. 0. .and. e1 .le. ionEvs(species)) then
     getUVcse=0. ; RETURN    ! [e0,e1] below ionization energy of species
  endif
  norm     = integrateSpectrum(X, Y, N, e0, e1, species, fdivLambda)
  getUVcse = integrateSpectrum(X, Y, N, e0, e1, species, fSigdivLambda)  &
           / norm
END FUNCTION getUVcse

!*************************************************************************
SUBROUTINE write_UVrates_table()
! Write the UV rates to a file (this is just in 
! debugging, to check if the UV spectra are being read correctly).
!-------------------------------------------------------------------------
  character(len=128)::filename
  integer::i, j
!-------------------------------------------------------------------------
  write(filename,'(A, I1, A)') 'UVrates.list'
  open(10, file=filename, status='unknown')
  write(10,*) UV_nz
  
  do i = 1,UV_nz
     write(10,900) UV_zeds(i), UV_rates_table(i,:,:)
  end do
  close(10)
900 format (f21.6, 20(1pe21.6))
END SUBROUTINE write_UVrates_table

!*************************************************************************
SUBROUTINE write_UVgroups_tables()

! Write the UV photon group properties to files (this is just in 
! debugging, to check if the UV spectra are being read correctly).
!-------------------------------------------------------------------------
  character(len=128)::filename
  integer::ip, i, j
!-------------------------------------------------------------------------
  do ip=1,nUVgroups
     write(filename,'(A, I1, A)') 'UVtable', ip, '.list'
     open(10, file=filename, status='unknown')
     write(10,*) UV_nz
     
     do i = 1,UV_nz
        write(10,901)                                                    &
                UV_zeds(i)           ,                                   &
                UV_groups_table(i,ip,1), UV_groups_table(i,ip,2),        &
                UV_groups_table(i,ip,3), UV_groups_table(i,ip,4),        &
                UV_groups_table(i,ip,5), UV_groups_table(i,ip,6),        &
                UV_groups_table(i,ip,7), UV_groups_table(i,ip,8)
     end do
     close(10)
  end do
901 format (f21.6,   f21.6,   f21.6,   1pe21.6, &
          & 1pe21.6, 1pe21.6, 1pe21.6, 1pe21.6, 1pe21.6  )
END SUBROUTINE write_UVgroups_tables

END MODULE UV_module

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!*************************************************************************
SUBROUTINE locate(xx,n,x,j)
! Locates position j of a value x in an ordered array xx of n elements
! After: xx(j) <= x <= xx(j+1) (assuming increasing order)
! j is lower bound, so it can be zero but not larger than n
!-------------------------------------------------------------------------
  use amr_commons,only:dp
  integer ::  n,j,jl,ju,jm
  real(dp)::  xx(n),x
!-------------------------------------------------------------------------
  jl = 0
  ju = n+1
  do while (ju-jl > 1) 
     jm = (ju+jl)/2
     if ((xx(n) > xx(1)) .eqv. (x > xx(jm))) then
        jl = jm
     else
        ju = jm
     endif
  enddo
  j = jl
END SUBROUTINE locate

!*************************************************************************
SUBROUTINE inp_1d(xax,nx,x,ix0,ix1,dx0,dx1)
! Compute variables by interpolation from table with non-equal intervals.
! xax      => Axis of x-values in table
! nx       => Length of x-axis
! x        => x-value to interpolate to
! ix0,ix1 <=  Lower and upper boundaries of x in xax
! dx0,dx1 <=  Weights of ix0 and ix1 indexes
!-------------------------------------------------------------------------
  use amr_commons,only:dp
  real(dp), intent(in)::xax(nx), x
  real(dp):: x_step, dx0, dx1
  integer:: nx, ix0, ix1
!-------------------------------------------------------------------------
  call locate(xax, nx, x, ix0)
  if(ix0 < 1) ix0=1
  if (ix0 < nx) then
     ix1  = ix0+1
     x_step = xax(ix1) - xax(ix0)
     dx0  = max(           x - xax(ix0), 0.0d0 ) / x_step
     dx1  = min(xax(ix1) - x           , x_step) / x_step
  else
     ix1  = ix0
     dx0  = 0.0d0 ;  dx1  = 1.0d0
  end if

  if (abs(dx0+dx1-1.0d0) .gt. 1.0d-5) then
     write(*,*) 'Screwed up the 1d interpolation ... '
     write(*,*) dx0+dx1
     call clean_stop
  end if
  !ret = dx0 * table(ix1, :) + dx1 * table(ix0, :)

END SUBROUTINE inp_1d

