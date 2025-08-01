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
 integer i,p
 real timefac,starage
 real ttemp,temp,xelec,xhii,xheii,xheiii,nh

 call initmem(nbodsmax,nsphmax,ncells)
 call basicinput

      print*,densconst,crionrate

! insert code after this
 timefac=timescale/year
 open(unit=1,file='sources.gabe',status='unknown')
 do p=nbodies-nstar+1,nbodies 
  starage=(tnow-tform(p))*timefac
  if(starage.LT.3.e7) write(1,'( 4g16.6 )') pos(p,1:3), Nlya(starage)*unitm_in_msun*mass(p)
 enddo
 close(1)

 ttemp=10000
 temp=10000
 nh=1.
 call Ionh(ttemp,temp,xelec,xhii,xheii,xheiii,nh,crionrate)
 print*,temp,xelec,xhii,xheii,xheiii

end program
