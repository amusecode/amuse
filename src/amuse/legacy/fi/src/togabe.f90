subroutine basicinput
 include 'globals.h'
 character*24 filenaam
 integer i,mapnr,nfeedback
 real,allocatable :: weights(:)
 
 call set_parameters(0)

 print*,'filenaam?'
 read*,filenaam

 call readbods(filenaam)

 call heattabel
        
 call initpars

 if(periodic) call initbc

 call initnewstar

 if(usepm) then
  if(verbosity.GT.0) print*,' ...initPM...'
  call initpm
 endif

 call postprocessread

end subroutine

program snapreader
 use StarsMod
 use ionMod
 include 'globals.h'
 integer i,p,nsources
 real timefac,starage,ltot,lout,l
 real ttemp,temp,xelec,nh,xh,xhe1,xhe2
 real,save ::xhii(nsphmax),xheii(nsphmax),xheiii(nsphmax)
 real(kind=4),save :: rblock(nsphmax)
 integer(kind=4),save :: iblock(nsphmax)
 integer(kind=8),save :: iblock2(nsphmax)
 real codemass,codevel
 integer*8,save :: ids(nsphmax)

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

 type(particle_header_type) phead
 character(ntagl) :: DataFieldLabel
 real hboxl,maxx,minx

 call initmem(nbodsmax,nsphmax,ncells)
 call basicinput

 do p=1,nsph 
  ttemp=temperat(p)
  nh=densconst*rho(p)
  call Ionh(ttemp,temp,xelec,xh,xhe1,xhe2,nh,crionrate)
  xhii(p)=xh
  xheii(p)=xhe1
  xheiii(p)=xhe2
!  print*, temperat(p),temp,elecfrac(p),xelec
 enddo

 timefac=timescale/year
 ltot=0
 lout=0
 nsources=0
 print*,'writing GABE sources'
 open(unit=1,file='sources.raw',status='unknown')
 do p=nbodies-nstar+1,nbodies 
  starage=(tnow-tform(p))*timefac
  l=Nlya(starage)*unitm_in_msun*mass(p)
  if(starage.LT.3.e7) then 
   nsources=nsources+1
   write(1,'( 4g16.6 )') pos(p,1:3), l
   lout=lout+l
  endif
  ltot=ltot+l 
 enddo
 close(1)
 write(*,'("wrote ",i," sources with ",g8.2," of total luminosity")'),nsources,lout/ltot

    print*, 'convert to GABE unit'

    phead%real(11)%var=3.08568e21 
    phead%real(12)%var=3.826e33
    phead%real(13)%var=3.1556926e14

   codemass=phead%real(12)%var * phead%real(13)%var**3 / phead%real(11)%var**2
   codevel=phead%real(11)%var/phead%real(13)%var
    print*, 'massunit (msol):', codemass/solarmass 
    print*, 'velunit (km/s):', codevel/1.e5 
    print*, 'length unit (kpc):',phead%real(11)%var/kpc

    mass(1:nsph)=mass(1:nsph)*unitm_in_msun*solarmass/codemass
    pos(1:nsph,1:3)=pos(1:nsph,1:3)*unitl_in_kpc/phead%real(11)%var*kpc
    hsmooth(1:nsph)=hsmooth(1:nsph)*unitl_in_kpc/phead%real(11)%var*kpc
    rho(1:nsph)=rho(1:nsph)*unitm_in_msun*solarmass/codemass &
                 /unitl_in_kpc**3*phead%real(11)%var**3/kpc**3
    vel(1:nsph,1:3)=vel(1:nsph,1:3)*velscale/codevel
    ids(1:nsph)=nbexist(1:nsph)
    tnow=tnow*timescale/phead%real(13)%var

    maxx=maxval(pos(1:nsph,1:3))
    minx=minval(pos(1:nsph,1:3))

    hboxl=max(maxx,abs(minx))

    phead%int(1:7)%name=i_header_variable_names(1:7)
    phead%real(1:21)%name=r_header_variable_names(1:21)

    phead%int(1)%var=nsph
    phead%int(2)%var=nsph
    phead%int(3)%var=1
    phead%int(4)%var=nsmooth
    phead%real(1)%var=0
    phead%real(2)%var=0
    phead%real(3)%var=0
    phead%real(4)%var=1.
    phead%real(8:10)%var=hboxl
    phead%real(5:7)%var=-hboxl
    phead%int(5)%var=0
    if(periodic) phead%int(5)%var=1
!    phead%MeanHsml=sum(hsmooth(1:nsph))/nsph
    phead%real(14)%var=0.3
    phead%real(15)%var=0.04
    phead%real(16)%var=0.7
    phead%real(17)%var=0.7
    phead%real(18)%var=0
    phead%real(19)%var=0
    phead%real(20)%var=0
    phead%real(21)%var=0
    phead%int(6)%var=0
    phead%int(7)%var=0
    phead%Datapresent=0.

    
    phead%DataPresent(1:9)=1
    phead%DataPresent(10:11)=1
    phead%DataPresent(12:14)=1
    phead%DataPresent(18)=1
    
! convert to gadget smoothing lengths
    hsmooth(1:nsph)=2*hsmooth(1:nsph)    
    
    print*, 'writing GABE file'
    open(unit=1,file='gas_000.1',status='unknown',form='UNFORMATTED')
    DataFieldLabel="header";write(1) DataFieldLabel, phead
    DataFieldLabel=data_field_names(1);iblock2(1:nsph)= ids;write(1) DataFieldLabel,iblock2(1:nsph)
    DataFieldLabel=data_field_names(2);rblock(1:nsph)=pos(1:nsph,1);write(1) DataFieldLabel,rblock(1:nsph)
    DataFieldLabel=data_field_names(3);rblock(1:nsph)=pos(1:nsph,2);write(1) DataFieldLabel,rblock(1:nsph)
    DataFieldLabel=data_field_names(4);rblock(1:nsph)=pos(1:nsph,3);write(1) DataFieldLabel,rblock(1:nsph)
    DataFieldLabel=data_field_names(5);rblock(1:nsph)=vel(1:nsph,1);write(1) DataFieldLabel,rblock(1:nsph)
    DataFieldLabel=data_field_names(6);rblock(1:nsph)=vel(1:nsph,2);write(1) DataFieldLabel,rblock(1:nsph)
    DataFieldLabel=data_field_names(7);rblock(1:nsph)=vel(1:nsph,3);write(1) DataFieldLabel,rblock(1:nsph)
    DataFieldLabel=data_field_names(8);rblock(1:nsph)= hsmooth(1:nsph);write(1) DataFieldLabel,rblock(1:nsph)
    DataFieldLabel=data_field_names(9);rblock(1:nsph)= rho(1:nsph);write(1) DataFieldLabel,rblock(1:nsph)
    DataFieldLabel=data_field_names(10);rblock(1:nsph)= mass(1:nsph);write(1) DataFieldLabel,rblock(1:nsph)
    DataFieldLabel=data_field_names(11);rblock(1:nsph)= temperat(1:nsph);write(1) DataFieldLabel,rblock(1:nsph)
    DataFieldLabel=data_field_names(12);rblock(1:nsph)= xhii(1:nsph);write(1) DataFieldLabel,rblock(1:nsph)
    DataFieldLabel=data_field_names(13);rblock(1:nsph)= xheii(1:nsph);write(1) DataFieldLabel,rblock(1:nsph)
    DataFieldLabel=data_field_names(14);rblock(1:nsph)= xheiii(1:nsph);write(1) DataFieldLabel,rblock(1:nsph)


!    rblock(1:nsph)= elecfrac(1:nsph);write(1) rblock(1:nsph)

    DataFieldLabel=data_field_names(18);iblock2(1:nsph)= 0   ;write(1) DataFieldLabel,iblock2(1:nsph)

    close(1)
    print*, 'done'
end program
