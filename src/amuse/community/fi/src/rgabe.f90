program read_snapshot

   include 'globals.h'
   integer :: NFILES        

   character*200 fileroot,filename
   character*6 nostring

!    this is an array of the possible data fields that can appear in an
!    output/restart snapshot.
  integer, parameter :: ntagl = 24

  ! names for the integer part of the header
  integer*4, parameter :: Nint_header_entries = 7
  character(ntagl), parameter :: i_header_variable_names(Nint_header_entries) = &
       (/ "NparSnap                ", &    ! 1
          "NparFile                ", &    ! 2
          "NumFiles                ", &    ! 3
          "NsphNeighbors           ", &    ! 4
          "RayBoundaryConditions   ", &    ! 5
          "OnTheSpot               ", &    ! 6
          "TotalRaysTraced         " /)    ! 7

  ! names for the real part of the header
  integer*4, parameter :: Nreal_header_entries = 18
  character(ntagl), parameter :: r_header_variable_names(Nreal_header_entries) = &
       (/ "Time                    ", &  ! 1
          "IsoMass                 ", &  ! 2
          "IsoTemperature          ", &  ! 3
          "ScaleFactor             ", &  ! 4
          "BoxLowX                 ", &  ! 5
          "BoxLowY                 ", &  ! 6
          "BoxLowZ                 ", &  ! 7
          "BoxHighX                ", &  ! 8
          "BoxHighY                ", &  ! 9
          "BoxHighZ                ", &  ! 10
          "CodeLenINcm             ", &  ! 11
          "CodeLumINergs_per_s     ", &  ! 12
          "CodeTimeINs             ", &  ! 13
          "OmegaMatter             ", &  ! 14
          "OmegaBaryon             ", &  ! 15
          "OmegaLambda             ", &  ! 16
          "LittleHubble            ", &  ! 17
          "RecombRayTolerance      " /)  ! 18


  ! names of the possible snapshot data fields
  integer*4, parameter :: ndatablocks = 30
  character(ntagl), parameter :: data_field_names(ndatablocks) = &  
       (/ "id                      ",&  ! 1
          "xpos                    ",&  ! 2
          "ypos                    ",&  ! 3
          "zpos                    ",&  ! 4
          "xvel                    ",&  ! 5
          "yvel                    ",&  ! 6
          "zvel                    ",&  ! 7
          "hsml                    ",&  ! 8
          "rho                     ",&  ! 9
          "mass                    ",&  ! 10
          "temperature             ",&  ! 11
          "xHII                    ",&  ! 12
          "xHeII                   ",&  ! 13
          "xHeIII                  ",&  ! 14
          "xHIIrc                  ",&  ! 15
          "xHeIIrc                 ",&  ! 16
          "xHeIIIrc                ",&  ! 17
          "lasthit                 ",&  ! 18
          "undefined               ",&
          "undefined               ",&
          "undefined               ",&
          "undefined               ",&
          "undefined               ",&
          "undefined               ",&
          "undefined               ",&
          "undefined               ",&
          "undefined               ",&
          "undefined               ",&
          "undefined               ",&
          "undefined               " /)

  type integer_header_entry
     character(ntagl) :: name
     integer*8 :: var
  end type integer_header_entry

  type real_header_entry
     character(ntagl) :: name
     real*8 :: var
  end type real_header_entry

  type particle_header_type
     type(integer_header_entry) :: int(Nint_header_entries)
     type(real_header_entry) :: real(Nreal_header_entries)
     integer*4 :: DataPresent(ndatablocks)
  end type particle_header_type

   real(kind=4),save :: rblock(nsphmax)
   integer(kind=4),save :: iblock(nsphmax)
   integer(kind=8),save :: iblock2(nsphmax)
   type(particle_header_type) :: phead
   real :: codemass,codevel,mscale,vscale,tscale,lscale
   integer no,k,i
      character(ntagl) :: DataFieldLabel

   call initmem(nbodsmax,nsphmax,ncells)

   unitm_in_msun=1.e9
   unitl_in_kpc=1.
   timescale=1.488e19*unitl_in_kpc*sqrt(unitl_in_kpc/unitm_in_msun)
   lengthscale=3.086e21*unitl_in_kpc
   velscale=lengthscale/timescale

   print*,'gabe file conversion'
   print*,'fileroot?'
   read*,fileroot
   
   filename=trim(fileroot)//'.unf'

   k=len_trim(fileroot)
   read(fileroot(k-2:k),*) no
   print *,'opening...  '//trim(filename)
   print*,' detected No.', no

! now, read in the header
   open (1, file=filename, form='unformatted')
   read(1) DataFieldLabel, phead

   print*,DataFieldLabel
   do i = 1,Nint_header_entries
    if (phead%int(i)%name /= i_header_variable_names(i)) then
      print*,phead%int(i)%name,i_header_variable_names(i)
      write(*,*) "integer part of snapshot header corrupted"
      stop
    end if
   end do

   do i = 1,Nreal_header_entries
    if (phead%real(i)%name /= r_header_variable_names(i)) then
      write(*,*) "real part of header corrupted"
      stop
    end if
   end do

        do i = 1,3
           if (phead%real(4+i)%var > phead%real(7+i)%var) then
              write(*,*) "Box Low > Box High ... exiting"
              stop
            end if
         end do


   codemass=phead%real(12)%var * phead%real(13)%var**3 / phead%real(11)%var**2
   codevel=phead%real(11)%var/phead%real(13)%var

   if(phead%int(3)%var.NE.1) then
    print*,'numfiles>1 not yet..'
    stop
   endif 
   
    print*,' assumed that all is gas'
   if(phead%int(1)%var.GT.nsphmax) then
    print*,'max number of particle exceed'
    stop
   endif 
   nsph=phead%int(1)%var
   nbodies=nsph

   if(phead%DataPresent(1).EQ.1) then 
    read(1) DataFieldLabel, iblock2(1:nsph)
   endif 
   
   if(phead%DataPresent(2).EQ.1) then
    read(1) DataFieldLabel, rblock(1:nsph) 
    pos(1:nsph,1)=rblock(1:nsph)
   endif 
   if(phead%DataPresent(3).EQ.1) then 
    read(1) DataFieldLabel,rblock(1:nsph)
    pos(1:nsph,2)=rblock(1:nsph)    
   endif
   if(phead%DataPresent(4).EQ.1) then 
    read(1) DataFieldLabel,rblock(1:nsph)
    pos(1:nsph,3)=rblock(1:nsph)
   endif 
  
   if(phead%DataPresent(5).EQ.1) then 
    read(1) DataFieldLabel,rblock(1:nsph);vel(1:nsph,1)=rblock(1:nsph)
   endif
   if(phead%DataPresent(6).EQ.1) then 
    read(1) DataFieldLabel,rblock(1:nsph); vel(1:nsph,2)=rblock(1:nsph)
   endif
   if(phead%DataPresent(7).EQ.1) then 
    read(1) DataFieldLabel,rblock(1:nsph); vel(1:nsph,3)=rblock(1:nsph)
   endif
   if(phead%DataPresent(8).EQ.1) then 
    read(1) DataFieldLabel,rblock(1:nsph); hsmooth(1:nsph)=rblock(1:nsph)
   endif
   if(phead%DataPresent(9).EQ.1) then 
    read(1) DataFieldLabel,rblock(1:nsph); rho(1:nsph)=rblock(1:nsph)
   endif
   if(phead%DataPresent(10).EQ.1) then 
    read(1) DataFieldLabel,rblock(1:nsph); mass(1:nsph)=rblock(1:nsph)
   endif
   if(phead%DataPresent(11).EQ.1) then 
    read(1) DataFieldLabel,rblock(1:nsph); temperat(1:nsph)=rblock(1:nsph)
   endif
   if(phead%DataPresent(12).EQ.1) then 
    read(1) DataFieldLabel,rblock(1:nsph); elecfrac(1:nsph)=rblock(1:nsph)
   endif
   if(phead%DataPresent(13).EQ.1) then 
    read(1) DataFieldLabel,rblock(1:nsph); elecfrac(1:nsph)=elecfrac(1:nsph)+0.1*rblock(1:nsph)
   endif
   if(phead%DataPresent(14).EQ.1) then 
    read(1) DataFieldLabel,rblock(1:nsph); elecfrac(1:nsph)=elecfrac(1:nsph)+0.1*2*rblock(1:nsph)
   endif
   if(phead%DataPresent(15).EQ.1) then 
    read(1) DataFieldLabel,rblock(1:nsph)
   endif
   if(phead%DataPresent(16).EQ.1) then 
    read(1) DataFieldLabel,rblock(1:nsph)
   endif
   if(phead%DataPresent(17).EQ.1) then 
    read(1) DataFieldLabel,rblock(1:nsph)
   endif
   if(phead%DataPresent(18).EQ.1) then 
    read(1) DataFieldLabel,iblock2(1:nsph)
   endif


   close(1)
   print *,'Done with reading.'

   mass(1:nsph)=mass(1:nsph)/unitm_in_msun/solarmass*codemass
   pos(1:nsph,1:3)=pos(1:nsph,1:3)/unitl_in_kpc*phead%real(11)%var/kpc
   hsmooth(1:nsph)=hsmooth(1:nsph)/unitl_in_kpc*phead%real(11)%var/kpc
   vel(1:nsph,1:3)=vel(1:nsph,1:3)/velscale*codevel
   
   totptag=nbodies
   tnow=0
!   tnow=phead%time/timescale*phead%CodeTime_s
   massres=0.

   if(phead%DataPresent(11).NE.1) then 
    temperat(1:nsph)=8000
   endif


   firstsnap=0
   uentropy=.FALSE.
   output=0
   output(1:3)=1
   output(13)=1
   output(17:18)=1
   
   write(nostring,'(I6.6)') no   
   call outbods(trim(fileroot(1:k-4)//'.'//nostring))
end program




