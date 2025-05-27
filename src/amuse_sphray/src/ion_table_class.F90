!> \file ion_table_class.F90

!> \brief Handles the CLOUDY ionization rate table used in some 
!! versions of Gadget simulations.
!<

module ion_table_class
use myf03_mod

#ifdef useHDF5
use hdf5_wrapper
#endif

#ifdef useMPI
use mpi
#endif

implicit none
private

real(r8b), parameter :: zero= 0.0d0
real(r8b), parameter :: half= 0.5d0
real(r8b), parameter :: one= 1.0d0
real(r8b), parameter :: two = 2.0d0
real(r8b), parameter :: four = 4.0d0

real(r8b), parameter :: pi = 3.14159265358979323846264338327950288
real(r8b), parameter :: four_pi= four * pi
real(r8b), parameter :: HIth_nu= 3.2899d15    ! [Hz]
real(r8b), parameter :: PLANCK= 6.626068d-27  ! [erg s]

 
public :: ion_table_type
public :: ion_header_type
public :: ion_spectrum_type
public :: mini_spec_type
public :: read_ion_table_file
public :: set_mini_spectrum
public :: return_gammaHI_at_z
public :: return_logflux_at_z
public :: integrate_flux_at_z
public :: interpolate_ion_table
public :: gammaHI_from_mini_spec_thin
public :: gammaHI_from_mini_spec_shield
public :: broadcast_ion_table


type ion_spectrum_type
   integer :: s_gammaHI
   integer :: s_z
   integer :: s_logryd
   integer :: s_logflux
   real(r8b), allocatable :: gammaHI(:)
   real(r8b), allocatable :: z(:)
   real(r8b), allocatable :: logryd(:)    ! [Rydbergs]
   real(r8b), allocatable :: logflux(:,:) ! [ergs/s/Hz/cm^2/sr]
   character(clen) :: model_name
end type ion_spectrum_type

type ion_header_type
   character(clen) :: cloudy_version
   type(ion_spectrum_type) :: ispec
end type ion_header_type


type ion_table_type
   type(ion_header_type) :: ihead
   integer :: s_ibal
   integer :: s_z
   integer :: s_logt
   integer :: s_logd
   real(r8b) :: logd_max
   real(r8b) :: logt_max
   real(r8b), allocatable :: ibal(:,:,:)  ! ibal( iz, it, id )
   real(r8b), allocatable :: z(:)
   real(r8b), allocatable :: logt(:)
   real(r8b), allocatable :: logd(:)
end type ion_table_type


!> stores a piece of the spectrum at a specific (interpolated) z
!-----------------------------------------------------------------
type mini_spec_type
   real(r8b) :: min_logryd                  !< lowest requested energy
   real(r8b) :: max_logryd                  !< highest requested energy
   real(r8b) :: z                           !< redshift
   real(r8b) :: gammaHI                     !< interpolated gammaHI from table
   real(r8b), allocatable :: logryd(:)      !< log energy samples
   real(r8b), allocatable :: ryd(:)         !< energy samples
   real(r8b), allocatable :: logflux(:)     !< interpolated log flux from nu1-nu2
   real(r8b), allocatable :: flux(:)        !< interpolated flux from nu1-nu2
   real(r8b), allocatable :: sigmaHI(:)     !< HI photoionization cs @ logryd
   real(r8b), allocatable :: sigma_ratio(:) !< sigma(nu) / sigma(nu_th)
end type mini_spec_type



contains


  !> sets up a mini spectrum
  !------------------------------------------------------------------------
  subroutine set_mini_spectrum( min_logryd, max_logryd, z, itab, mini_spec )
    real(r8b), intent(in) :: min_logryd
    real(r8b), intent(in) :: max_logryd
    real(r8b), intent(in) :: z
    type(ion_table_type), intent(in) :: itab
    type(mini_spec_type), intent(out) :: mini_spec

    real(r8b), allocatable :: logflux_all(:) !< interp logflux @ all energy
                                             
    integer(i4b) :: nall, n
    integer(i4b) :: ilwr, iupr
    integer(i4b) :: i

    mini_spec%min_logryd = min_logryd
    mini_spec%max_logryd = max_logryd
    mini_spec%z          = z

    mini_spec%gammaHI = return_gammaHI_at_z( itab, z )

    ! first get interplated (in redshift) logflux at all energies
    !-------------------------------------------------------------
    nall = size(itab%ihead%ispec%logryd)
    allocate( logflux_all( nall ) )
    logflux_all = return_logflux_at_z( itab, z )

    ! now find indices that bracket requested energies
    !---------------------------------------------------
    do i = 1,nall
       if ( itab%ihead%ispec%logryd(i) > min_logryd ) then
          ilwr = i-1
          exit
       endif
    end do

    do i = 1,nall
       if ( itab%ihead%ispec%logryd(i) > max_logryd ) then
          iupr = i
          exit
       endif
    end do

    n = iupr - ilwr + 1
    allocate( mini_spec%logryd(n), mini_spec%ryd(n), mini_spec%logflux(n),      &
              mini_spec%flux(n), mini_spec%sigmaHI(n), mini_spec%sigma_ratio(n) )

    mini_spec%logryd  = itab%ihead%ispec%logryd(ilwr:iupr)
    mini_spec%logflux = logflux_all(ilwr:iupr)
    mini_spec%flux    = 10**mini_spec%logflux
    mini_spec%ryd     = 10**mini_spec%logryd
    do i = 1,n
       mini_spec%sigmaHI(i)     = HI_photo_cs_verner( mini_spec%ryd(i) )
       mini_spec%sigma_ratio(i) = HI_photo_cs_verner( mini_spec%ryd(i) ) / HI_photo_cs_verner( 1.0d0 ) 
    end do




  end subroutine set_mini_spectrum


  !> calculates gammaHI from mini spectrum in optically thin limit
  !-------------------------------------------------------------------------
  function gammaHI_from_mini_spec_thin( mini_spec ) result( gammaHI )
    type(mini_spec_type), intent(in) :: mini_spec
    real(r8b) :: gammaHI

    gammaHI = int_tabulated_trap_rule( mini_spec%logryd, &
         four_pi * log(10.d0) * mini_spec%flux * mini_spec%sigmaHI / PLANCK )   

  end function gammaHI_from_mini_spec_thin



  !> calculates gammaHI from mini spectrum in optically w/ shielding
  !-------------------------------------------------------------------------
  function gammaHI_from_mini_spec_shield( mini_spec, tauHI) result( gammaHI )
    type(mini_spec_type), intent(in) :: mini_spec
    real(r8b), intent(in) :: tauHI
    real(r8b) :: gammaHI

    gammaHI = int_tabulated_trap_rule( mini_spec%logryd, &
         four_pi * log(10.d0) * mini_spec%flux * mini_spec%sigmaHI / PLANCK * &
         exp(-tauHI * mini_spec%sigma_ratio) )   

  end function gammaHI_from_mini_spec_shield





  !> HI photo ionization x-section (Analytic) [cm^2]
  !-------------------------------------------------------------------------
  function HI_photo_cs_analytic( Ry ) result( sigma )
    real(r8b), intent(in) :: Ry  !< energy [Rydbergs]
    real(r8b) :: sigma           !< cross section [cm^2]
    real(r8b), parameter :: A0 = 6.30d-18


    real(r8b) :: eps 
    real(r8b) :: ieps
    real(r8b) :: efac
    real(r8b) :: num
    real(r8b) :: den

    eps = sqrt( Ry - one )
    ieps = one/eps
    efac = four - ( four * atan(eps) ) * ieps
    num = exp(efac)
    den = one - exp(-two * pi * ieps)
  
    sigma = A0 * (Ry)**(-4) * num / den


  end function HI_photo_cs_analytic



  !> HI photo ionization x-section at 1 Rydberg (Verner) [cm^2]
  !-------------------------------------------------------------------------
  function HIth_photo_cs_verner() result( sigma )
    real(r8b) :: sigma
    sigma = HI_photo_cs_verner( one )
  end function HIth_photo_cs_verner

 
  !> HI photo ionization x-section (Verner) [cm^2]
  !-------------------------------------------------------------------------
  function HI_photo_cs_verner( Ry ) result( sigma )
    real(r8b), intent(in) :: Ry  !< energy [Rydbergs]
    real(r8b) :: sigma           !< cross section [cm^2]
    real(r8b), parameter :: Eth = 13.6d0
    real(r8b), parameter :: Emax = 5.0d4
    real(r8b), parameter :: E0 = 4.298d-1
    real(r8b), parameter :: sig0 = 5.475d4
    real(r8b), parameter :: ya = 3.288d1
    real(r8b), parameter :: P = 2.963d0

    real(r8b) :: eV
    real(r8b) :: x
    real(r8b) :: y

    eV = Ry * 13.6d0
    x = eV / E0
    y = x
  
    sigma = sig0 * (x-1)**2 * y**(0.5d0 * P - 5.5d0) * (1 + sqrt(y/ya))**(-P)
    sigma = sigma * 1.0d-18

  end function HI_photo_cs_verner


  !> HI photo ionization x-section (Osterbrok) [cm^2]
  !----------------------------------------------------------------------------
  function HI_photo_cs_osterbrok( Ry ) result( sigma )    
    real(r8b), intent(in) :: Ry    !< energy [Rydbergs]
    real(r8b) :: sigma             !< cross section 
    real(r8b), parameter :: sigma0 = 6.3d-18

    if (Ry < one) then
       sigma = zero
    else
       sigma = sigma0 * ( Ry )**(-3)
    endif

  end function HI_photo_cs_osterbrok




  ! integrates a tabulated function using the trap. rule
  !-------------------------------------------------------------------------
  function int_tabulated_trap_rule( xarr, yarr ) result(sum)
    real(r8b), intent(in) :: xarr(:)
    real(r8b), intent(in) :: yarr(:)
    real(r8b) :: sum
    real(r8b) :: h
    integer :: n
    integer :: i

    n = size(xarr)
    if ( size(yarr) /= n ) then
       write(*,*) 'size mismatch in int_tabulated_trap_rule'
       stop
    endif

    sum = zero    
    do i = 1, n-1
       h = xarr(i+1) - xarr(i)
       sum = sum + h * half * (yarr(i+1) + yarr(i))
    end do

  end function int_tabulated_trap_rule



  ! integrate flux
  !---------------------------------------------
  function integrate_flux_at_z( itab, z ) result(sum)
    type(ion_table_type) :: itab
    real(r8b) :: z
    real(r8b) :: sum

    real(r8b), allocatable :: logfluxall(:) 
    real(r8b), allocatable :: ryd(:)
    real(r8b), allocatable :: sigma(:)
    real(r8b), allocatable :: nu(:)
    real(r8b), allocatable :: logflux(:)
    real(r8b), allocatable :: xarr(:)
    real(r8b), allocatable :: yarr(:)

    integer :: n
    integer :: i
    integer :: icnt
    integer :: nabove
    real(r8b) :: GHI


    n = size(itab%ihead%ispec%logryd)
    allocate( logfluxall(n) )
    logfluxall = return_logflux_at_z( itab, z )
    

    icnt = 0
    do i = 1,n
       if ( itab%ihead%ispec%logryd(i) > zero ) then
          icnt = i-1
          exit
       end if
    end do

    nabove = size( itab%ihead%ispec%logryd(icnt:icnt+49) )


    allocate( ryd(nabove), sigma(nabove), nu(nabove), logflux(nabove) )
    allocate( xarr(nabove), yarr(nabove) )

    xarr = itab%ihead%ispec%logryd(icnt:icnt+49)
    logflux = logfluxall(icnt:icnt+49)
    ryd = 10**xarr
    ryd(1) = 1.0d0
    do i = 1,size(ryd) 
       sigma(i) = HI_photo_cs_verner( ryd(i) )
    end do
    nu = ryd * HIth_nu
        
    yarr = four_pi * log(10.d0) * 10**logflux * sigma / PLANCK
    
    sum = int_tabulated_trap_rule( xarr, yarr )

    GHI = return_gammaHI_at_z( itab, z )


  end function integrate_flux_at_z





  ! returns the HI photoionization rate at z
  !---------------------------------------------
  function return_gammaHI_at_z( itab, z ) result(GHI)
    type(ion_table_type) :: itab
    real(r8b) :: z
    real(r8b) :: GHI

    real(r8b) :: zsearch, zlow, zhigh
    real(r8b) :: dz1, dz2, z1, z2
    real(r8b) :: frac1, frac2, dGHI
    integer :: n, i, iz1, iz2

    zsearch = z
    n = size(itab%ihead%ispec%z)

    zlow  = itab%ihead%ispec%z(1)
    zhigh = itab%ihead%ispec%z(n)
    
    if (z <= zlow) then
       GHI = itab%ihead%ispec%gammaHI(1)
       return
    endif
    if (z >= zhigh) then
       GHI = itab%ihead%ispec%gammaHI(n)
       return
    endif

    ! bracket the redshift
    do i = 2, n
       if (itab%ihead%ispec%z(i) > z) then
          iz2 = i
          iz1 = iz2-1
          exit
       end if
    end do

    z1 = itab%ihead%ispec%z(iz1)
    z2 = itab%ihead%ispec%z(iz2)

    dz1 = z  - z1
    dz2 = z2 - z

    frac1 = dz1 / (z2-z1)
    frac2 = dz2 / (z2-z1)

    dGHI = itab%ihead%ispec%gammaHI(iz2) - itab%ihead%ispec%gammaHI(iz1)
    GHI = itab%ihead%ispec%gammaHI(iz1) + frac1 * dGHI

  end function return_gammaHI_at_z



  ! returns logryd and logflux at z
  !---------------------------------------------
  function return_logflux_at_z( itab, z ) result(logflux)
    type(ion_table_type) :: itab
    real(r8b) :: z
    real(r8b) :: logflux(size(itab%ihead%ispec%logryd))

    real(r8b) :: zsearch, zlow, zhigh
    real(r8b) :: dz1, dz2, z1, z2
    real(r8b) :: frac1, frac2, dflux
    integer :: n, i, iz1, iz2

    zsearch = z
    n = size(itab%ihead%ispec%z)

    zlow  = itab%ihead%ispec%z(1)
    zhigh = itab%ihead%ispec%z(n)
    
    if (z <= zlow) then
       logflux = itab%ihead%ispec%logflux(1,:)
       return
    endif
    if (z >= zhigh) then
       logflux = itab%ihead%ispec%logflux(n,:)
       return
    endif

    ! bracket the redshift
    do i = 2, n
       if (itab%ihead%ispec%z(i) > z) then
          iz2 = i
          iz1 = iz2-1
          exit
       end if
    end do

    z1 = itab%ihead%ispec%z(iz1)
    z2 = itab%ihead%ispec%z(iz2)

    dz1 = z  - z1
    dz2 = z2 - z

    frac1 = dz1 / (z2-z1)
    frac2 = dz2 / (z2-z1)

    do i = 1, size(itab%ihead%ispec%logryd)
       dflux = itab%ihead%ispec%logflux(iz2,i) - itab%ihead%ispec%logflux(iz1,i)
       logflux(i) = itab%ihead%ispec%logflux(iz1,i) + frac1 * dflux
    enddo



  end function return_logflux_at_z








  subroutine read_ion_table_file( file, itab )
    character(clen), parameter :: myname="read_ion_table_file" 
    logical, parameter :: crash=.true.
    integer, parameter :: verb=2
    character(clen) :: str

    character(*) :: file
    type(ion_table_type) :: itab
    integer :: fh

#ifdef useHDF5    
    integer :: rank
    integer :: dims(3)
    character(clen) :: fmt

    integer :: nz, nt, nd

    fmt = "(T7, A, 4I8)"
   
    call hdf5_open_file( fh, file, readonly=.true. )

    call hdf5_read_attribute(fh, 'header/cloudy_version', itab%ihead%cloudy_version)
    call hdf5_read_attribute(fh, 'header/spectrum/model_name', itab%ihead%ispec%model_name)

!    call mywrite('',verb) 
    call mywrite("      cloudy version = "//trim(itab%ihead%cloudy_version), verb)
    call mywrite("      model name     = "//trim(itab%ihead%ispec%model_name), verb)
    call mywrite('',verb) 

    dims=0
    call hdf5_get_dimensions(fh, 'header/spectrum/gammahi', rank, dims)
    allocate( itab%ihead%ispec%gammaHI( dims(1) ) )
    call hdf5_read_data(fh, 'header/spectrum/gammahi', itab%ihead%ispec%gammaHI )
    write(str,fmt) "rank and dims of gammaHI:       ", rank, dims 
    call mywrite(str,verb)
    itab%ihead%ispec%s_gammaHI = dims(1) 

    dims=0
    call hdf5_get_dimensions(fh, 'header/spectrum/logenergy_ryd', rank, dims)
    allocate( itab%ihead%ispec%logryd( dims(1) ) )
    call hdf5_read_data(fh, 'header/spectrum/logenergy_ryd', itab%ihead%ispec%logryd )
    write(str,fmt) "rank and dims of logenergy_ryd: ", rank, dims
    call mywrite(str,verb) 
    itab%ihead%ispec%s_logryd = dims(1)

    dims=0
    call hdf5_get_dimensions(fh, 'header/spectrum/logflux', rank, dims)
    allocate( itab%ihead%ispec%logflux( dims(1), dims(2) ) )
    call hdf5_read_data(fh, 'header/spectrum/logflux', itab%ihead%ispec%logflux )
    write(str,fmt) "rank and dims of logflux:       ", rank, dims
    call mywrite(str,verb) 
    itab%ihead%ispec%s_logflux = product( dims(1:2) )

    dims=0
    call hdf5_get_dimensions(fh, 'header/spectrum/redshift', rank, dims)
    allocate( itab%ihead%ispec%z( dims(1) ) )
    call hdf5_read_data(fh, 'header/spectrum/redshift', itab%ihead%ispec%z )
    write(str,fmt) "rank and dims of redshift:      ", rank, dims 
    call mywrite(str,verb) 
    itab%ihead%ispec%s_z = dims(1)

    dims=0
    call hdf5_get_dimensions(fh, 'ionbal', rank, dims)
    allocate( itab%ibal( dims(1), dims(2), dims(3) ) )
    call hdf5_read_data(fh, 'ionbal', itab%ibal )
    write(str,fmt) "rank and dims of ionbal:        ", rank, dims
    call mywrite(str,verb) 
    itab%s_ibal = product( dims(1:3) )

    dims=0
    call hdf5_get_dimensions(fh, 'logd', rank, dims)
    allocate( itab%logd( dims(1) ) )
    call hdf5_read_data(fh, 'logd', itab%logd )
    write(str,fmt) "rank and dims of logd:          ", rank, dims
    call mywrite(str,verb) 
    nd = dims(1)
    itab%s_logd = dims(1) 

    dims=0
    call hdf5_get_dimensions(fh, 'logt', rank, dims)
    allocate( itab%logt( dims(1) ) )
    call hdf5_read_data(fh, 'logt', itab%logt )
    write(str,fmt) "rank and dims of logt:          ", rank, dims 
    call mywrite(str,verb) 
    nt = dims(1)
    itab%s_logt = dims(1) 

    dims=0
    call hdf5_get_dimensions(fh, 'redshift', rank, dims)
    allocate( itab%z( dims(1) ) )
    call hdf5_read_data(fh, 'redshift', itab%z )
    write(str,fmt) "rank and dims of z:             ", rank, dims
    call mywrite(str,verb) 
    nz = dims(1)
    itab%s_z = dims(1) 

    call hdf5_close_file( fh )

    itab%ibal = log10( itab%ibal )

    itab%logd_max = itab%logd(nd)
    itab%logt_max = itab%logt(nt)


    fmt = "(T7, A, 2ES15.5)"

    call mywrite('',verb) 
    call mywrite("     spectrum:",verb)

    write(str,fmt) "min/max G_HI    = ", minval(itab%ihead%ispec%gammaHI), maxval(itab%ihead%ispec%gammaHI)
    call mywrite(str,verb)

    write(str,fmt) "min/max logryd  = ", minval(itab%ihead%ispec%logryd), maxval(itab%ihead%ispec%logryd)
    call mywrite(str,verb)

    write(str,fmt) "min/max logflux = ", minval(itab%ihead%ispec%logflux), maxval(itab%ihead%ispec%logflux)
    call mywrite(str,verb)

    write(str,fmt) "min/max z       = ", minval(itab%ihead%ispec%z), maxval(itab%ihead%ispec%z)
    call mywrite(str,verb)

    call mywrite('',verb)
    call mywrite("     ionization balance:",verb)
    write(str,fmt) "min/max log10(ionbal)  = ", minval(itab%ibal), maxval(itab%ibal)
    call mywrite(str,verb)
 
    write(str,fmt) "min/max logd           = ", minval(itab%logd), maxval(itab%logd)
    call mywrite(str,verb)

    write(str,fmt) "min/max logt           = ", minval(itab%logt), maxval(itab%logt)
    call mywrite(str,verb)

    write(str,fmt) "min/max z              = ", minval(itab%z), maxval(itab%z)
    call mywrite(str,verb)
    call mywrite('',verb) 

#else

    stop "trying to read hdf5 ion table, but haven't defined hdf5 macro in Makefile"

#endif

  end subroutine read_ion_table_file




  function interpolate_ion_table( itab, z, logt, logd ) result (ioneq)
    type(ion_table_type) :: itab
    real(r8b) :: z      !< input redshift
    real(r8b) :: logt   !< input log temperature
    real(r8b) :: logd   !< input log density
    real(r8b) :: ioneq  !< returned ionization equilibrium

    integer :: nz, nt, nd !< size of interpolation grid
    integer :: iz, it, id !< lower bracketing index
    integer :: fz, ft, fd !< upper bracketing index
    real(r8b)    :: rz, rt, rd !< real index (between bracketing index)

    real(r8b) :: wiz, wfz, wit, wft, wid, wfd  !< fractional distances between bracketing indices
    real(r8b) :: w111, w211, w121, w221, w112, w212, w122, w222  !< weights
    real(r8b) :: dz, dt, dd                    !< difference between bracketing indices

    integer :: i !< general loop index

    ! get sizes of interpolation grid
    !----------------------------------
    nz = size(itab%z)
    nt = size(itab%logt)
    nd = size(itab%logd)

    ! do some checks
    !-------------------    
    if ( z < itab%z(1) )  then
       write(*,*) "warning z below limit in table"
       stop
    endif

    if ( z > itab%z(nz) ) then
       write(*,*) "warning z above limit in table"
       stop
    endif

    if ( logt < itab%logt(1) )  then
!       write(*,*) "warning logt below limit in table"
       logt = itab%logt(1)
    endif

    if ( logt > itab%logt(nt) ) then
!       write(*,*) "warning logt above limit in table"
       logt = itab%logt(nt)
    endif

    if ( logd < itab%logd(1) )  then
!       write(*,*) "warning logd below limit in table"
       logd = itab%logd(1)
    endif

    if ( logd > itab%logd(nd) ) then
!       write(*,*) "warning logd above limit in table"
       logd = itab%logd(nd)
    endif

    
    

    ! find bracketing indices
    !--------------------------
    do i = 1, nz
       if (itab%z(i) > z) then
          fz = i
          iz = fz-1
          exit
       end if
    end do

    do i = 1, nt
       if (itab%logt(i) > logt) then
          ft = i
          it = ft-1
          exit
       end if
    end do

    do i = 1, nd
       if (itab%logd(i) > logd) then
          fd = i
          id = fd-1
          exit
       end if
    end do


    if (z==itab%z(1)) then
       iz = 1
       fz = 2
    else if (z==itab%z(nz)) then
       fz = nz
       iz = nz-1
    endif

    if (logt==itab%logt(1)) then
       it = 1
       ft = 2
    else if (logt==itab%logt(nt)) then
       ft = nt
       it = nt-1
    endif

    if (logd==itab%logd(1)) then
       id = 1
       fd = 2
    else if (logd==itab%logd(nd)) then
       fd = nd
       id = nd-1
    endif


    dz = itab%z(fz) - itab%z(iz)
    wiz = ( z - itab%z(iz) ) / dz
    wfz = ( itab%z(fz) - z ) / dz

    dt = itab%logt(ft) - itab%logt(it)
    wit = ( logt - itab%logt(it) ) / dt
    wft = ( itab%logt(ft) - logt ) / dt

    dd = itab%logd(fd) - itab%logd(id)
    wid = ( logd - itab%logd(id) ) / dd
    wfd = ( itab%logd(fd) - logd ) / dd

    rz = iz + wiz
    rt = it + wit
    rd = id + wid

 
    w111 = wfz * wft * wfd
    w211 = wiz * wft * wfd
    w121 = wfz * wit * wfd
    w221 = wiz * wit * wfd
    w112 = wfz * wft * wid
    w212 = wiz * wft * wid
    w122 = wfz * wit * wid
    w222 = wiz * wit * wid


    ioneq = &
         w111 * itab%ibal( iz, it, id ) + &
         w211 * itab%ibal( fz, it, id ) + &
         w121 * itab%ibal( iz, ft, id ) + &
         w221 * itab%ibal( fz, ft, id ) + &
         w112 * itab%ibal( iz, it, fd ) + &
         w212 * itab%ibal( fz, it, fd ) + &
         w122 * itab%ibal( iz, ft, fd ) + &
         w222 * itab%ibal( fz, ft, fd ) 

    
    ioneq = 10**ioneq



   
  end function interpolate_ion_table



  !> broadcasts the ionization table to all processing elements
  !--------------------------------------------------------------
  subroutine broadcast_ion_table( MyPE, itab )
    integer, intent(in) :: MyPE
    type(ion_table_type) :: itab
    integer :: count
    integer :: root
    integer :: ierr
    
    root = 0

#ifdef useMPI

    count = 1
    call mpi_bcast( itab%logd_max, count, mpi_double_precision, root, mpi_comm_world, ierr )
    call mpi_bcast( itab%logt_max, count, mpi_double_precision, root, mpi_comm_world, ierr )

    call mpi_bcast( itab%s_logd,    count, mpi_integer, root, mpi_comm_world, ierr )
    call mpi_bcast( itab%s_logt,    count, mpi_integer, root, mpi_comm_world, ierr )
    call mpi_bcast( itab%s_z,       count, mpi_integer, root, mpi_comm_world, ierr )
    call mpi_bcast( itab%s_ibal,    count, mpi_integer, root, mpi_comm_world, ierr )

    call mpi_bcast( itab%ihead%ispec%s_logflux, count, mpi_integer, root, mpi_comm_world, ierr )
    call mpi_bcast( itab%ihead%ispec%s_logryd,  count, mpi_integer, root, mpi_comm_world, ierr )
    call mpi_bcast( itab%ihead%ispec%s_z,       count, mpi_integer, root, mpi_comm_world, ierr )
    call mpi_bcast( itab%ihead%ispec%s_gammaHI, count, mpi_integer, root, mpi_comm_world, ierr )


    call mpi_barrier(mpi_comm_world, ierr)
    if (MyPE /= root) then
       allocate( itab%logd( itab%s_logd ) )
       allocate( itab%logt( itab%s_logt ) )
       allocate( itab%z   ( itab%s_z    ) )
       allocate( itab%ibal( itab%s_z, itab%s_logt, itab%s_logd ) )

       allocate( itab%ihead%ispec%logflux ( itab%ihead%ispec%s_z, itab%ihead%ispec%s_logryd  ) )
       allocate( itab%ihead%ispec%logryd  ( itab%ihead%ispec%s_logryd  ) )
       allocate( itab%ihead%ispec%z       ( itab%ihead%ispec%s_z       ) )
       allocate( itab%ihead%ispec%gammaHI ( itab%ihead%ispec%s_gammaHI ) )
    endif
    call mpi_barrier(mpi_comm_world, ierr)


    count = itab%s_logd  
    call mpi_bcast( itab%logd, count, mpi_double_precision, root, mpi_comm_world, ierr )

    count = itab%s_logt  
    call mpi_bcast( itab%logt, count, mpi_double_precision, root, mpi_comm_world, ierr )

    count = itab%s_z 
    call mpi_bcast( itab%z, count, mpi_double_precision, root, mpi_comm_world, ierr )

    count = itab%s_ibal 
    call mpi_bcast( itab%ibal, count, mpi_double_precision, root, mpi_comm_world, ierr )


    count = clen
    call mpi_bcast( itab%ihead%cloudy_version,   count, mpi_character, root, mpi_comm_world, ierr )
    call mpi_bcast( itab%ihead%ispec%model_name, count, mpi_character, root, mpi_comm_world, ierr )

    count = itab%ihead%ispec%s_logflux
    call mpi_bcast( itab%ihead%ispec%logflux(1,1), count, mpi_double_precision, root, mpi_comm_world, ierr )

    count = itab%ihead%ispec%s_logryd 
    call mpi_bcast( itab%ihead%ispec%logryd, count, mpi_double_precision, root, mpi_comm_world, ierr )

    count = itab%ihead%ispec%s_z 
    call mpi_bcast( itab%ihead%ispec%z, count, mpi_double_precision, root, mpi_comm_world, ierr )

    count = itab%ihead%ispec%s_gammaHI 
    call mpi_bcast( itab%ihead%ispec%gammaHI, count, mpi_double_precision, root, mpi_comm_world, ierr )



#endif


  end subroutine broadcast_ion_table


end module ion_table_class
