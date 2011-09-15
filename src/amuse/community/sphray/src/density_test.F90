!> \file density_test.F90

!> \brief Module contains calls to the density test (in development)
!<

! this program illustrates the use of the density routines..
program densitytest
  use particle_system_mod, only: particle_system_type
  use snapshot_mod, only: particle_header_type
  use snapshot_mod, only: read_particle_header
  use input_mod, only: read_initial_snapshot
  use nblist_mod, only: nblist_type, prepare_neighbor_search
  use kernel_density_mod, only: initdensity, hsmdens, density
  use global_mod, only: global_variables_type
  use config_mod, only: config_variables_type
  use sphpar_mod, only: sphpar_type
  use physical_constants_mod
  use myf90, only: clen
  implicit none

  
  type(particle_system_type) :: psys
  type(particle_header_type) :: phead
  type(nblist_type) :: nblist
  type(sphpar_type) :: lpart
  type(global_variables_type) :: gvars
  type(config_variables_type) :: cvars

  real :: small
  real :: cellsize(3)
  integer :: cells1d(3)
  real :: zonesize(3),halfzone(3),quarzone(3)  

  character(clen) :: snapfile,rhofile,ionfile,tempfile
  character(clen) :: buffer
  integer :: i,j,k

  real :: lpos(3)
  real :: DfltHsml, DfltMass, DfltTemp
  real :: zcoord
  real, allocatable :: rho(:,:), xHII(:,:), T(:,:)

  real :: rhomin, rhomax
  real :: xHIImin, xHIImax
  real :: maxmass
  logical :: verbose,closefile
  integer :: lun
  real :: MBalloc
  real :: hh

    call getarg(1,snapfile)
    call getarg(2,rhofile)
    call getarg(3,ionfile)
    call getarg(4,tempfile)
    call getarg(5,buffer)
    read(buffer,*) zonesize(1)
    call getarg(5,buffer)
    read(buffer,*) zonesize(2)
    call getarg(5,buffer)
    read(buffer,*) zonesize(3)

    cells1D(1) = 256
    cells1D(2) = 256
    cells1D(3) = 256

!    zonesize(3) = zonesize(3)/2
!    cells1D(3) = cells1D(3)/2

    allocate( rho(cells1d(1),cells1d(3)) )
    allocate( xHII(cells1d(1),cells1d(3)) )
    allocate( T(cells1d(1),cells1d(3)) )

    halfzone = zonesize / 2.0
    quarzone = zonesize / 4.0
    cellsize = zonesize / cells1d
    write(*,*) "cellsize = ", cellsize
    write(*,*) "zonesize = ", zonesize
    
    verbose = .false. 
    closefile = .true.
    call read_particle_header(snapfile,verbose,closefile,phead,lun)
    call set_necessary_globals()
    call read_initial_snapshot(gvars,snapfile,verbose,psys,MBalloc)

    rhomin = minval(psys%par(:)%rho)
    rhomax = maxval(psys%par(:)%rho)
    
    xHIImin = minval(psys%par(:)%xHII)
    xHIImax = maxval(psys%par(:)%xHII)

  

! everything up to this point just reads in a particle system, assigns
! global variables, and changes the boundry conditions to periodic.

! intiialize kernel density routines  
#ifdef incmass
  maxmass = maxval(psys%par(:)%mass)
  call initdensity(maxmass)
#else
  call initdensity(gvars%iso_mass)
#endif

! initialize nblist search 
  call prepare_neighbor_search(psys,nblist,DfltMass,DfltTemp)
  write(*,*) 

! this calculates the density at the particle positions
! (and outputs) note that the actual density calculation uses a 
! 'local' particle type, which is a particle with some special 
!  members exclusive for density calculation 
! (this also allows the density to be probed at any point, because lpart
! need not be initialized from real particles!!)





  write(*,*) "default smoothing length = ", DfltHsml
  write(*,*) "npar = ", psys%npar
  write(*,*) "min/max readin rho", rhomin, rhomax
  write(*,*) "min/max readin xHII", xHIImin, xHIImax
#ifdef incT
  write(*,*) "min/max T", minval(psys%par(:)%T), maxval(psys%par(:)%T)
#endif
  write(*,*) "zonesize = ", zonesize

  write(*,*) "tops = ", psys%box%top
  write(*,*) "bots = ", psys%box%bot


10 format (A,3I4,3ES10.2)
100 format(A,I5,A,3F8.2)

  hh=DfltHsml
  do i = 1,cells1d(1)

     if (mod(i,64)==0) write(*,'(A,F10.2)') "% done ", real(i)/cells1d(1) 

     do j = 1,cells1d(3)
!        do k = 1,cells1d

        lpart%pos(1) = -halfzone(1) + (i-0.5) * cellsize(1) 
        lpart%pos(2) = -halfzone(3) + (j-0.5) * cellsize(3)
        lpart%pos(3) = 0.0

 
        ! if we are outside the computational volume just set to 0.0
        
        small = 0.0        
        if (lpart%pos(1).lt.psys%box%bot(1)+small .or. &
            lpart%pos(2).lt.psys%box%bot(2)+small .or. & 
            lpart%pos(3).lt.psys%box%bot(3)+small .or. &
            lpart%pos(1).gt.psys%box%top(1)-small .or. &
            lpart%pos(2).gt.psys%box%top(2)-small .or. & 
            lpart%pos(3).gt.psys%box%top(3)-small ) then
            rho(i,j) = 0.0
            xHII(i,j) = 0.0
            T(i,j) = 0.0
            cycle
         end if

        lpart%rho=1.0
        lpart%xHII=0.0
        lpart%T =0.0
        lpart%nnb=0
        lpart%gradrho=0.
        lpart%drhodh=0.
        lpart%hsmooth=hh
        lpart%fi=0.
        lpart%dfi=0.

        call hsmdens(lpart,psys,nblist)
        call density(lpart,psys,nblist)

!        if (lpart%xHII==1.00) lpart%xHII = xHIImax

        ! stop execution if averaged value is below particle extremes
        if (lpart%rho.le.0.0) stop "rho .le. 0.0"

!        if (lpart%rho.lt.rhomin-small) then
!           write(*,*) "***************************"
!           write(*,*) "rho below min", lpart%rho
!           write(*,*) "pos = ", lpart%pos
!           write(*,*) "nnb = ", lpart%nnb
!           write(*,*) "***************************"
!           stop
!           lpart%rho = rhomin
!        end if

        if (lpart%rho.gt.rhomax+small) then
           write(*,*) "***************************"
           write(*,*) "rho above max", lpart%rho
           write(*,*) "pos = ", lpart%pos
           write(*,*) "nnb = ", lpart%nnb
           write(*,*) "***************************"
           stop
           lpart%rho = rhomax
        end if
        

        if (lpart%xHII.lt. xHIImin) then
           write(*,*) "***************************"
           write(*,*) "xHII below min", lpart%xHII
           write(*,*) "pos = ", lpart%pos
           write(*,*) "nnb = ", lpart%nnb
           write(*,*) "***************************"
           stop
           lpart%xHII = xHIImin
        else if (lpart%xHII.gt.xHIImax) then
           write(*,*) "***************************"
           write(*,*) "xHII above max", lpart%xHII
           write(*,*) "pos = ", lpart%pos
           write(*,*) "nnb = ", lpart%nnb
           write(*,*) "***************************"
           stop
           lpart%xHII = xHIImax
        end if

        rho(i,j)  = lpart%rho
        xHII(i,j) = lpart%xHII
        T(i,j)    = lpart%T
        hh= lpart%hsmooth

!        write(*,*) "rho = ", rho(i,j)
!        write(*,*) "xHII = ", xHII(i,j)
!        write(*,*) "T    = ", T(i,j)

     end do
  end do

  write(*,*) "min/max rho", minval( rho ), maxval( rho )
  write(*,*) "min/max xHII", minval( xHII ), maxval( xHII )
  write(*,*) "min/max T", minval(T), maxval(T)
  write(*,*) 
  write(*,*) "==========================================================="

  open(unit=10,file=rhofile,status="UNKNOWN",form="unformatted")  
  open(unit=20,file=ionfile,status="UNKNOWN",form="unformatted")  
  open(unit=30,file=tempfile,status="UNKNOWN",form="unformatted")  
  write(10) rho
  write(20) xHII
  write(30) T
  close(10)
  close(20)
  close(30)

contains

 subroutine set_necessary_globals()
 
   integer :: NsphNnb, Npar
   real :: BoxVol, VperPar

   DfltMass = phead%real(2)%var
   DfltTemp = phead%real(3)%var
   NsphNnb = phead%int(4)%var

   Npar = phead%int(1)%var
   gvars%iso_mass = DfltMass
   cvars%IsoMass = DfltMass
   cvars%NsphNnb = NsphNnb

   psys%box%bbound=0
   psys%box%tbound=0
   psys%box%bot = minval(phead%real(5:7)%var)
   psys%box%top = maxval(phead%real(8:10)%var)
 
   DfltHsml = 0.0
   
 end subroutine set_necessary_globals

end program densitytest
