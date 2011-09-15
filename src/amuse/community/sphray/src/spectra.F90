!> \file spectra.F90

!> \brief spectra handling
!<

module spectra_mod
use myf03_mod
implicit none
private

  public :: rn2freq, read_spectra_file

  integer(i8b), parameter :: maxfreqbins = 1000  !< maximum frequency bins

  integer(i8b) :: Nspectra      !< number of user defined spectra
  integer(i8b) :: Nfreqs        !< number of entries

  real(r8b), allocatable :: nus(:,:)  !< frequency in HI ionizing units (SpecNum,freq)
  real(r8b), allocatable :: cdf(:,:)  !< cumulative distribution function (SpecNum,cdf)

contains


!>   reads all user defined spectra 
!---------------------------------------
  subroutine read_spectra_file(spectra_file)

    character(clen), intent(in) :: spectra_file  !< file containing spectra

    character(clen), parameter :: myname = 'read_spectra_file'
    logical, parameter :: crash = .true.
    integer, parameter :: verb = 1
    character(clen) :: str

    real(r8b) :: Numin, Numax
    integer(i4b) :: lun, err, i, j
    logical :: fthere


    write(str,'(A,A)') 'using spectra file:       ', trim(spectra_file) 
    call mywrite(str, verb) 

    inquire(file=spectra_file,exist=fthere)
    if (.not. fthere) then
       call myerr('cant find spectra file: ' // trim(spectra_file), myname, crash)
    end if

    call open_formatted_file_r(spectra_file,lun)
    read(lun,*) Nspectra
    if (Nspectra<=0) Nspectra=1  ! just to make sure spectra is allocated

    write(str,'(A,I5)') '  number of user defined spectra = ', Nspectra 
    call mywrite('', verb+1)
    call mywrite(str, verb+1) 

    do i = 1,Nspectra

       read(lun,*) Nfreqs, Numin, Numax

       write(str,'(A,I5,A,I5,A)') '  user defined spectrum ',i,' has ', &
                                  Nfreqs,' entries'
       call mywrite(str, verb+1) 

       if (i == 1) then
          allocate ( nus(Nfreqs,Nspectra), cdf(Nfreqs,Nspectra), stat=err )
          if(err/=0) then
             call myerr('cant allocate nus and cdf', myname, crash) 
          end if
       end if

       do j = 1,Nfreqs
          read(lun,*) nus(j,i), cdf(j,i)
       end do

    end do

    close(lun)

    call mywrite('', verb+1) 


  end subroutine read_spectra_file


  !> returns a frequency given a spectral type
  !-------------------------------------------
  function rn2freq(spectype) result(freq)
  use mt19937_mod, only: genrand_real1

    character(clen), parameter :: myname = 'rn2freq'
    logical, parameter :: crash = .true.

    real(r4b) :: spectype  !< spectral user defined type
    real(r4b) :: freq      !< returned frequency

    real(r8b) :: rn
    integer(i8b) :: specbin
    integer(i8b) :: specnum

    if (spectype <= 0.0) then
       freq = abs(spectype)
       return
    else
       specnum = int(spectype)
       rn = genrand_real1()
       do specbin = 1,Nfreqs
          if ( rn <= cdf(specbin,specnum) ) then
             freq = nus(specbin,specnum)
             return
          end if
       end do

       call myerr("went through loop w/o finding frequency", myname, crash)

    end if

   end function rn2freq

end module spectra_mod
