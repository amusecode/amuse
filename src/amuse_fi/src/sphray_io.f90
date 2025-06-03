module sphray_io_mod
implicit none

public
public :: read_sphray_snapshot,write_sphray_snapshot,sphray_snap_type

integer, parameter, private :: i4b = selected_int_kind(9)       
integer, parameter, private :: i8b = selected_int_kind(18)      
integer, parameter, private :: r4b  = selected_real_kind(6,37)  
integer, parameter, private :: r8b  = selected_real_kind(15,307)
integer, parameter, private :: taglen = 24
integer, parameter, private :: NdataBlocks = 30
integer, parameter, private :: NheadEntries = 25

  ! names of the header entries
  character(taglen), parameter :: Htags(NheadEntries) = &  
       (/ "NparSnap                ", "NparFile                ", &  ! 1,2
          "NumFiles                ", "NsphNeighbors           ", &  ! 3,4
          "RayBoundaryConditions   ", "OnTheSpot               ", &  ! 5,6
          "TotalRaysTraced         ", "Time                    ", &  ! 7,1
          "IsoMass                 ", "IsoTemperature          ", &  ! 2,3
          "ScaleFactor             ", "BoxLowX                 ", &  ! 4,5
          "BoxLowY                 ", "BoxLowZ                 ", &  ! 6,7
          "BoxHighX                ", "BoxHighY                ", &  ! 8,9
          "BoxHighZ                ", "CodeLenINcm             ", &  ! 10,11
          "CodeLumINergs           ", "CodeTimeINs             ", &  ! 12,13 
          "OmegaMatter             ", "OmegaBaryon             ", &  ! 14,15
          "OmegaLambda             ", "LittleHubble            ", &  ! 16,17
          "RecombRayTolerance      " /)                              ! 18

  ! names of the possible snapshot data fields
  character(taglen), parameter :: Dtags(NdataBlocks) = &  
       (/ "id                      ", "xpos                    ", &  ! 1,2
          "ypos                    ", "zpos                    ", &  ! 3,4
          "xvel                    ", "yvel                    ", &  ! 5,6
          "zvel                    ", "hsml                    ", &  ! 7,8
          "rho                     ", "mass                    ", &  ! 9,10
          "temperature             ", "xHII                    ", &  ! 11,12
          "xHeII                   ", "xHeIII                  ", &  ! 13,14
          "xHIIrc                  ", "xHeIIrc                 ", &  ! 15,16
          "xHeIIIrc                ", "lasthit                 ", &  ! 17,18
          "undefined               ", "undefined               ", &  ! 19,20
          "undefined               ", "undefined               ", &  ! 21,22
          "undefined               ", "undefined               ", &  ! 23,24
          "undefined               ", "undefined               ", &  ! 25,26
          "undefined               ", "undefined               ", &  ! 27,28
          "undefined               ", "undefined               " /)  ! 29,30


  
type sphray_snap_type
! header entries
 integer(i8b) :: NparSnap
 integer(i8b) :: NparFile
 integer(i8b) :: Nfiles
 integer(i8b) :: NsphNbrs
 integer(i8b) :: RayBCs
 integer(i8b) :: OnTheSpot
 integer(i8b) :: NraysTraced

 real(r8b) :: Time
 real(r8b) :: IsoMass
 real(r8b) :: IsoTemp
 real(r8b) :: ScaleFac
 real(r8b) :: BoxLow(3)
 real(r8b) :: BoxHigh(3)
 real(r8b) :: CGSlen
 real(r8b) :: CGSlum
 real(r8b) :: CGStime
 real(r8b) :: OmegaM
 real(r8b) :: OmegaB
 real(r8b) :: OmegaL
 real(r8b) :: LittleH
 real(r8b) :: RecRayTol

 integer(i4b) :: DataPresent(30)
 
! data blocks
 integer(i8b), allocatable :: id(:)    !< particle id
 real(r4b), allocatable    :: pos(:,:) !< x,y,z coordinates
 real(r4b), allocatable    :: vel(:,:) !< x,y,z velocities
 real(r4b), allocatable    :: hsml(:)  !< smoothing length
 real(r4b), allocatable    :: rho(:)   !< density = mass * NsphNnb / hsml^3 
 real(r4b), allocatable    :: mass(:)  !< particle mass
 real(r4b), allocatable    :: T(:)     !< temperature in K       
 real(r4b), allocatable    :: xHII(:)  !< HII ionization fraction
 real(r4b), allocatable    :: xHeII(:) !< HeII ionization fraction
 real(r4b), allocatable    :: xHeIII(:)   !< HeIII ionization fraction
 real(r4b), allocatable    :: xHIIrc(:)   !< HII recombination fraction
 real(r4b), allocatable    :: xHeIIrc(:)  !< HeII recombination fraction
 real(r4b), allocatable    :: xHeIIIrc(:) !< HeIII recombination fraction
 integer(i8b), allocatable :: lasthit(:)  !< last ray to cross this particle

end type

contains

subroutine read_sphray_snapshot(s,file)
  type(sphray_snap_type) :: s
  integer, parameter :: lun = 28
  character(200), optional :: file
  character(200) :: snapfile
  integer :: i
  integer(i8b) :: Npar
  character(taglen) :: tag

  if(present(file)) then
   snapfile=file
  else 
   write(*,*) "enter snapshot file to read:"
   read(*,"(A)") snapfile
  endif
  
  write(*,*) "reading ", trim(snapfile)
  open(unit=lun,file=snapfile,status="old",form="unformatted")
  read(lun) tag, &
            tag, s%NparSnap,    tag, s%NparFile,   tag, s%Nfiles, &
            tag, s%NsphNbrs,    tag, s%RayBCs,     tag, s%OnTheSpot,  &
            tag, s%NraysTraced, tag, s%Time,       tag, s%IsoMass,    &
            tag, s%IsoTemp,     tag, s%ScaleFac,   tag, s%BoxLow(1),  &
            tag, s%BoxLow(2),   tag, s%BoxLow(3),  tag, s%BoxHigh(1), &
            tag, s%BoxHigh(2),  tag, s%BoxHigh(3), tag, s%CGSlen,     &
            tag, s%CGSlum,      tag, s%CGStime,    tag, s%OmegaM,     &
            tag, s%OmegaB,      tag, s%OmegaL,     tag, s%LittleH,    &
            tag, s%RecRayTol, s%DataPresent
  
  
  do i = 1,18
     if (s%DataPresent(i) .EQ. 1) then
        write(*,*) Dtags(i), " present"
     else if (s%DataPresent(i) .EQ. 0) then
        write(*,*) Dtags(i), " absent"
     end if
  end do
  Npar = s%NparSnap

  write(*,*) "Npar = ", Npar
 
  if (s%DataPresent(1)) then
     allocate(s%id(Npar))
     read(lun) tag,s%id
     write(*,*) "min/max id = ", minval(s%id), maxval(s%id)
  end if
  
  if (s%DataPresent(2)) then
     allocate(s%pos(Npar,3))
     read(lun) tag, s%pos(:,1)
     read(lun) tag, s%pos(:,2)
     read(lun) tag, s%pos(:,3)
     write(*,*) "min/max pos = ", minval(s%pos), maxval(s%pos)
  end if
  
  if (s%DataPresent(5)) then
     allocate(s%vel(Npar,3))
     read(lun) tag, s%vel(:,1)
     read(lun) tag, s%vel(:,2)
     read(lun) tag, s%vel(:,3)
     write(*,*) "min/max vel = ", minval(s%vel), maxval(s%vel)
  end if
  
  if (s%DataPresent(8)) then
     allocate(s%hsml(Npar))
     read(lun) tag, s%hsml
     write(*,*) "min/max hsml = ", minval(s%hsml), maxval(s%hsml)
  end if
  
  if (s%DataPresent(9)) then
     allocate(s%rho(Npar))
     read(lun) tag, s%rho
     write(*,*) "min/max rho = ", minval(s%rho), maxval(s%rho)
  end if
  
  if (s%DataPresent(10)) then
     allocate(s%mass(Npar))
     read(lun) tag, s%mass
     write(*,*) "min/max mass = ", minval(s%mass), maxval(s%mass)
  end if
  
  if (s%DataPresent(11)) then
     allocate(s%T(Npar))
     read(lun) tag, s%T
     write(*,*) "min/max T = ", minval(s%T), maxval(s%T)
  end if
  
  if (s%DataPresent(12)) then
     allocate(s%xHII(Npar))
     read(lun) tag, s%xHII
     write(*,*) "min/max xHII = ", minval(s%xHII), maxval(s%xHII)
  end if
  
  if (s%DataPresent(13)) then
     allocate(s%xHeII(Npar))
     read(lun) tag, s%xHeII
     write(*,*) "min/max xHeII = ", minval(s%xHeII), maxval(s%xHeII)
  end if
  
  if (s%DataPresent(14)) then
     allocate(s%xHeIII(Npar))
     read(lun) tag, s%xHeIII
     write(*,*) "min/max xHeIII = ", minval(s%xHeIII), maxval(s%xHeIII)
  end if
  
  if (s%DataPresent(15)) then
     allocate(s%xHIIrc(NPar))
     read(lun) tag, s%xHIIrc
     write(*,*) "min/max xHIIrc = ", minval(s%xHIIrc), maxval(s%xHIIrc)
  end if
  
  if (s%DataPresent(16)) then
     allocate(s%xHeIIrc(Npar))
     read(lun) tag, s%xHeIIrc
     write(*,*) "min/max xHeIIrc = ", minval(s%xHeIIrc), maxval(s%xHeIIrc)
  end if
  
  if (s%DataPresent(17)) then
     allocate(s%xHeIIIrc(Npar))
     read(lun) tag, s%xHeIIIrc
     write(*,*) "min/max xHeIIIrc = ", minval(s%xHeIIIrc), maxval(s%xHeIIIrc)
  end if
  
  if (s%DataPresent(18)) then
     allocate(s%lasthit(Npar))
     read(lun) tag, s%lasthit
     write(*,*) "min/max lasthit = ", minval(s%lasthit), maxval(s%lasthit) 
  end if
  
end subroutine read_sphray_snapshot
 

! writes whatever is loaded into the variables in this module into a 
! snapshot. 
subroutine write_sphray_snapshot(s,file)
 type(sphray_snap_type) :: s
 character(200),optional :: file
 character(200) :: newfile 
 integer(i8b) :: i
 character(taglen) :: tag

   if(present(file)) then
    newfile=file
   else 
    newfile = "this_is_a_new_sphray_snapshot_change_my_name.unf"
   endif
   
   open(unit=150,file=newfile,form="unformatted")

   tag = "header"
   write(150) tag, &
        Htags(1),  s%NparSnap,    Htags(2),  s%NparFile,   Htags(3),  s%Nfiles, &
        Htags(4),  s%NsphNbrs,    Htags(5),  s%RayBCs,     Htags(6),  s%OnTheSpot,  &
        Htags(7),  s%NraysTraced, Htags(8),  s%Time,       Htags(9),  s%IsoMass,    &
        Htags(10), s%IsoTemp,     Htags(11), s%ScaleFac,   Htags(12), s%BoxLow(1),  &
        Htags(13), s%BoxLow(2),   Htags(14), s%BoxLow(3),  Htags(15), s%BoxHigh(1), &
        Htags(16), s%BoxHigh(2),  Htags(17), s%BoxHigh(3), Htags(18), s%CGSlen,     &
        Htags(19), s%CGSlum,      Htags(20), s%CGStime,    Htags(21), s%OmegaM,     &
        Htags(22), s%OmegaB,      Htags(23), s%OmegaL,     Htags(24), s%LittleH,    &
        Htags(25), s%RecRayTol, s%DataPresent

   if (s%DataPresent(1)) write(150) Dtags(1), s%id 
   if (s%DataPresent(2)) write(150) Dtags(2), s%pos(:,1) 
   if (s%DataPresent(3)) write(150) Dtags(3), s%pos(:,2) 
   if (s%DataPresent(4)) write(150) Dtags(4), s%pos(:,3) 
   if (s%DataPresent(5)) write(150) Dtags(5), s%vel(:,1) 
   if (s%DataPresent(6)) write(150) Dtags(6), s%vel(:,2) 
   if (s%DataPresent(7)) write(150) Dtags(7), s%vel(:,3)
   if (s%DataPresent(8)) write(150) Dtags(8), s%hsml
   if (s%DataPresent(9)) write(150) Dtags(9), s%rho
   if (s%DataPresent(10)) write(150) Dtags(10), s%mass
   if (s%DataPresent(11)) write(150) Dtags(11), s%T
   if (s%DataPresent(12)) write(150) Dtags(12), s%xHII
   if (s%DataPresent(13)) write(150) Dtags(13), s%xHeII
   if (s%DataPresent(14)) write(150) Dtags(14), s%xHeIII
   if (s%DataPresent(15)) write(150) Dtags(15), s%xHIIrc
   if (s%DataPresent(16)) write(150) Dtags(16), s%xHeIIrc
   if (s%DataPresent(17)) write(150) Dtags(17), s%xHeIIIrc
   if (s%DataPresent(18)) write(150) Dtags(18), s%lasthit
   close(150)

end subroutine write_sphray_snapshot

subroutine alloc_sphray_snapshot(s)
 type(sphray_snap_type) :: s
 integer(i8b) :: Npar
 
  Npar = s%NparSnap

  if(Npar.LE.0) call error('set number of particles first')
   
  if (s%DataPresent(1)) then
     allocate(s%id(Npar))
  end if
  
  if (s%DataPresent(2)) then
     allocate(s%pos(Npar,3))
  end if
  
  if (s%DataPresent(5)) then
     allocate(s%vel(Npar,3))
  end if
  
  if (s%DataPresent(8)) then
     allocate(s%hsml(Npar))
  end if
  
  if (s%DataPresent(9)) then
     allocate(s%rho(Npar))
  end if
  
  if (s%DataPresent(10)) then
     allocate(s%mass(Npar))
  end if
  
  if (s%DataPresent(11)) then
     allocate(s%T(Npar))
  end if
  
  if (s%DataPresent(12)) then
     allocate(s%xHII(Npar))
  end if
  
  if (s%DataPresent(13)) then
     allocate(s%xHeII(Npar))
  end if
  
  if (s%DataPresent(14)) then
     allocate(s%xHeIII(Npar))
  end if
  
  if (s%DataPresent(15)) then
     allocate(s%xHIIrc(NPar))
  end if
  
  if (s%DataPresent(16)) then
     allocate(s%xHeIIrc(Npar))
  end if
  
  if (s%DataPresent(17)) then
     allocate(s%xHeIIIrc(Npar))
  end if
  
  if (s%DataPresent(18)) then
     allocate(s%lasthit(Npar))
  end if

end subroutine
  
end module sphray_io_mod


