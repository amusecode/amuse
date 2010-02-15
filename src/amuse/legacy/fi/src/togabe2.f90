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
 use sphray_io_mod
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
 type(sphray_snap_type) :: s
 character*200 :: outfilename

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

    s%CGSlen=3.08568e21 
    s%CGSlum=3.826e33
    s%CGStime=3.1556926e14
    codemass=s%CGSlum * s%CGStime**3/ s%CGSlen**2
    codevel=s%CGSlen/s%CGStime
    print*, 'massunit (msol):', codemass/solarmass 
    print*, 'velunit (km/s):', codevel/1.e5 
    print*, 'length unit (kpc):',s%CGSlen/kpc

    mass(1:nsph)=mass(1:nsph)*unitm_in_msun*solarmass/codemass
    pos(1:nsph,1:3)=pos(1:nsph,1:3)*unitl_in_kpc/s%CGSlen*kpc
    hsmooth(1:nsph)=hsmooth(1:nsph)*unitl_in_kpc/s%CGSlen*kpc
    rho(1:nsph)=rho(1:nsph)*unitm_in_msun*solarmass/codemass &
                 /unitl_in_kpc**3*s%CGSlen**3/kpc**3
    vel(1:nsph,1:3)=vel(1:nsph,1:3)*velscale/codevel
    ids(1:nsph)=nbexist(1:nsph)
    tnow=tnow*timescale/s%CGStime

    maxx=maxval(pos(1:nsph,1:3))
    minx=minval(pos(1:nsph,1:3))

    hboxl=max(maxx,abs(minx))

!    phead%int(1:7)%name=i_header_variable_names(1:7)
!    phead%real(1:21)%name=r_header_variable_names(1:21)

 s%NparSnap=nsph
 s%NparFile=nsph
 s%Nfiles=1
 s%NsphNbrs=nsmooth
 s%RayBCs=0
 if(periodic) s%RayBCS=1
 s%OnTheSpot=0
 s%NraysTraced=0

 s%Time=0
 s%IsoMass=0
 s%IsoTemp=0
 s%ScaleFac=1.
 s%BoxLow(1:3)=-hboxl
 s%BoxHigh(1:3)=hboxl
 s%OmegaM=0.3
 s%OmegaB=0.04
 s%OmegaL=0.7
 s%LittleH=0.7
 s%RecRayTol=0

 s%DataPresent=0
 s%DataPresent(1:9)=1
 s%DataPresent(10:11)=1
 s%DataPresent(12:14)=1
 s%DataPresent(18)=1
    
 call alloc_sphray_snapshot(s)
    
! convert to gadget smoothing lengths
    hsmooth(1:nsph)=2*hsmooth(1:nsph)    
    
    s%id(1:nsph)= ids(1:nsph)
    s%pos(1:nsph,1)=pos(1:nsph,1)
    s%pos(1:nsph,2)=pos(1:nsph,2)
    s%pos(1:nsph,3)=pos(1:nsph,3)
    s%vel(1:nsph,1)=vel(1:nsph,1)
    s%vel(1:nsph,2)=vel(1:nsph,2)
    s%vel(1:nsph,3)=vel(1:nsph,3)
    s%hsml(1:nsph)= hsmooth(1:nsph)
    s%rho(1:nsph)= rho(1:nsph)
    s%mass(1:nsph)= mass(1:nsph)
    s%T(1:nsph)= temperat(1:nsph)
    s%xHII(1:nsph)= xhii(1:nsph)
    s%xHeII(1:nsph)= xheii(1:nsph)
    s%xHeIII(1:nsph)= xheiii(1:nsph)
    s%lasthit(1:nsph)=0

    print*, 'writing GABE file'
    outfilename='gas_000.1'
    call write_sphray_snapshot(s,outfilename)
    
    print*, 'done'

end program
