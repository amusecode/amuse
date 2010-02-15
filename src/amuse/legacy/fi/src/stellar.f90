subroutine starevolv
  use StarsMod
  include 'globals.h'
     
  integer p
  real efuvnew,starage,timefac,fluxfac
      
  if(cosmo) call terror('no cosmo starevolv (yet)')
 
  timefac=timescale/year
  fluxfac=heatconst*unitm_in_msun*10**flxscale
 
  do p=nbodies-nstar+1,nbodies-nbh
    starage=(tnow-tform(p))*timefac
    efuvnew=fluxfac*dFUV(starage)
    starfuv(p)=efuvnew*mass(p)
    if(starfuv(p).LT.0.) starfuv(p)=0.0
  enddo
 
  do p=nbodies-nbh+1,nbodies
    starfuv(p)=0.
  enddo
   
end subroutine


subroutine initnewstar
  use StarsMod
  use ElementsMod
  include 'globals.h'

  if(.not.radiate) return
  call InitElements(metallicity)
  if(verbosity.GT.0)print*,' ...InitStars...'
  call InitStars(datadir,zQ())
  
end subroutine


subroutine newstar
include 'globals.h'	 
  real mjeans,tcsound2,ranfb,newstarmass,dvx,dvy,dvz
  logical formstar,force
  integer nform,nstable,p,i,ib,nforced,ioerror
  character*4 option         
  integer, allocatable :: forced(:)
  real fthresh
  parameter(fthresh=0.5)

  allocate(forced(nsph))
   newstarmass=sfeff*8.*massres/nsmooth  

! 3 modes:
!  gerritsen: sf mass fraction of gas mass, constant delay 
!  nieuw: sfmass constant, delay with time inversely prop. to gas mass
! ( in order to keep physical star formation rate independent of particle mass)
!  molecular: sf mass fraction of gas mass, only if enough h2frac!
	  
  nform=0 
  nstable=0
  nforced=0
              
  if(sfmode.EQ.'gerritsen') then
    do p=1,nsph
      if(rho(p).EQ.0) cycle
      force=.false.
      tcsound2=csound(p)**2
      mjeans=rho(p)*pi/6.*(tcsound2*pi/rho(p))**1.5     
      if(mjeans.GT.masscrit) then
        tcollaps(p)=tnow        
        nstable=nstable+1
      endif
      formstar=(tnow-tcollaps(p))*SQRT(4*pi*rho(p)).GT.tcollfac
      if(rho(p).GT.rhomax.AND.mass(p).GT.0) then
        formstar=.TRUE.
        force=.TRUE.
      endif
      formstar=formstar.AND.(tnow-tform(p)).GT.tsnbeg.AND.mass(p).GT.0    
      if(formstar) then
        nform=nform+1
        bodlist(nform)=p       
        forced(nform)=0
        if(force) forced(nform)=1        
      endif
    enddo
  endif

  if(sfmode.EQ.'nieuw') then
    do p=1,nsph
      if(rho(p).EQ.0) cycle
      force=.false.
      tcsound2=csound(p)**2
      mjeans=rho(p)*pi/6.*(tcsound2*pi/rho(p))**1.5   
      if(mjeans.GT.masscrit) then
        tcollaps(p)=tnow       
        nstable=nstable+1
      endif
      formstar=mjeans.LE.masscrit.AND.  &
        exp(-(tnow-tstarform)*SQRT(4*pi*rho(p))* &
        mass(p)*sfeff/(tcollfac*newstarmass)).LT.ranfb(rnseed)
      if(rho(p).GT.rhomax.AND.mass(p).GT.0) then
        formstar=.TRUE.
        force=.TRUE.
      endif
      formstar=formstar.AND.(tnow-tform(p)).GT.tsnbeg.AND.mass(p).GT.0    
      if(formstar) then
        nform=nform+1
        bodlist(nform)=p       
        forced(nform)=0
        if(force) forced(nform)=1        
      endif
    enddo
  endif

  if(sfmode.EQ.'molecular') then
    do p=1,nsph
      if(rho(p).EQ.0) cycle
      force=.false.
      tcsound2=csound(p)**2
      mjeans=rho(p)*pi/6.*(tcsound2*pi/rho(p))**1.5   
      if(mjeans.GT.masscrit) then
        tcollaps(p)=tnow       
        nstable=nstable+1
      endif
      formstar=mjeans.LE.masscrit.AND.h2frac(p).GT.sfeff+fthresh.AND.  &
        exp(-(tnow-tstarform)*SQRT(4*pi*rho(p))/tcollfac).LT.ranfb(rnseed)
      if(rho(p).GT.rhomax.AND.mass(p).GT.0) then
        formstar=.TRUE.
        force=.TRUE.
      endif
      formstar=formstar.AND.(tnow-tform(p)).GT.tsnbeg.AND.mass(p).GT.0    
      if(formstar) then
        nform=nform+1
        bodlist(nform)=p       
        forced(nform)=0
        if(force) forced(nform)=1        
      endif
    enddo
  endif

  if(sfmode.EQ.'molecular2') then
    do p=1,nsph          
      if(rho(p).EQ.0) cycle
      force=.false.
      tcsound2=csound(p)**2
      mjeans=rho(p)*pi/6.*(tcsound2*pi/rho(p))**1.5	  
      if(mjeans.GT.masscrit) then
        tcollaps(p)=tnow	  
        nstable=nstable+1
      endif
! note that the condition on sfeff is a numerical req. and not 
! part of physical requirement (in which case I think it should
! be part of resetting clock)
      formstar=h2frac(p).GT.sfeff.AND. &   
        (tnow-tcollaps(p))*SQRT(4*pi*rho(p)).GT.tcollfac
      if(rho(p).GT.rhomax.AND.mass(p).GT.0) then
        formstar=.TRUE.
        force=.TRUE.
      endif
      formstar=formstar.AND.(tnow-tform(p)).GT.tsnbeg.AND.mass(p).GT.0	   
      if(formstar) then
        nform=nform+1
        bodlist(nform)=p	  
        forced(nform)=0
        if(force) forced(nform)=1	    
      endif
    enddo
  endif
	 
  if(nbodies+nform.GT.nbodsmax) CALL terror('too many stars')

! shift bh to end
  if(nform.GT.0) then
    do i=nbodies+nform,nbodies+nform-nbh+1,-1
      p=i-nform
      mass(i)=mass(p)
      pos(i,1)=pos(p,1)
      pos(i,2)=pos(p,2)
      pos(i,3)=pos(p,3)
      vel(i,1)=vel(p,1)
      vel(i,2)=vel(p,2)
      vel(i,3)=vel(p,3)
      acc(i,1)=acc(p,1)
      acc(i,2)=acc(p,2)
      acc(i,3)=acc(p,3)
      acc(i,4)=acc(p,4)
      phi(i)=phi(p)
      phiext(i)=phiext(p)
      epsgrav(i)=epsgrav(p)
      itimestp(i)=itimestp(p)
      otimestp(i)=otimestp(p)
      tform(i)=tform(p)
      starfuv(i)=starfuv(p)
      tfeedb(i)=tfeedb(p)
      tvel(i)=tvel(p)
      snentropy(i)=snentropy(p)
      hsmooth(i)=hsmooth(p)
      nbexist(i)=nbexist(p)
    enddo

    open(unit=sfrfilenr,file=TRIM(outputfile)//'.sfh', &
           status='unknown',POSITION='APPEND',IOSTAT=ioerror)

    do i=1,nform
      p=nbodies+i-nbh
      ib=bodlist(i)
      tform(p)=tnow
      mass(p)=newstarmass
      if(sfmode.EQ.'gerritsen'.OR.sfmode.EQ.'molecular') &
        mass(p)=sfeff*mass(ib)
      if(sfmode.EQ.'molecular2') mass(p)=h2frac(ib)*mass(ib)
      if(mass(ib)-mass(p).lt.removgas.OR.forced(i).EQ.1) &
        mass(p)=mass(ib)
      pos(p,1)=pos(ib,1)
      pos(p,2)=pos(ib,2)
      pos(p,3)=pos(ib,3)
      vel(p,1)=vel(ib,1)
      vel(p,2)=vel(ib,2)
      vel(p,3)=vel(ib,3)
      acc(p,1)=acc(ib,1)
      acc(p,2)=acc(ib,2)
      acc(p,3)=acc(ib,3)
      acc(p,4)=acc(ib,4)
      phi(p)=phi(ib)
      phiext(p)=phiext(ib)
      tfeedb(p)=-2
      itimestp(p)=itimestp(ib)
      otimestp(p)=otimestp(ib)
      tvel(p)=tvel(ib)
      epsgrav(p)=eps
      if(adaptive_eps.OR.eps.EQ.0) epsgrav(p)=epsgrav(ib)       
!           epsgrav(p)=epsgrav(ib)
! check
      snentropy(p)=0.
      hsmooth(p)=hsmooth(ib)
      nbexist(p)=totptag+i

      if(ioerror.EQ.0) then
        write(sfrfilenr,*) ib,p,forced(i)
        write(sfrfilenr,*) tnow,tcollaps(ib)
        write(sfrfilenr,*) pos(p,1),pos(p,2),pos(p,3)
        write(sfrfilenr,*) vel(p,1),vel(p,2),vel(p,3)
        write(sfrfilenr,*) rho(ib)
        write(sfrfilenr,*) temperat(ib)          
        write(sfrfilenr,*) elecfrac(ib)
        write(sfrfilenr,*) fuvheat(ib)
        write(sfrfilenr,*) h2frac(ib)
        write(sfrfilenr,*) mass(p),mass(ib)
      endif
      tform(ib)=tnow
      tcollaps(ib)=tnow
      mass(ib)=MAX(mass(ib)-mass(p),0.)       
      if(sfmode.EQ.'molecular'.AND.sfeff.LT.1.) then
        h2frac(ib)=max((h2frac(ib)-sfeff)/(1.-sfeff),0.) 
      endif
      if(sfmode.EQ.'molecular2'.AND.sfeff.LT.1.) then
        h2frac(ib)=max((h2frac(ib)-sfeff)/(1.-sfeff),0.) 
      endif
      if(forced(i).EQ.1) nforced=nforced+1
    
      if(mass(ib)+mass(p).NE.0) then
        dvx=csound(ib)*(2*ranfb(rnseed)-1)
        dvy=csound(ib)*(2*ranfb(rnseed)-1)
        dvz=csound(ib)*(2*ranfb(rnseed)-1)
        vel(p,1)=vel(p,1)+dvx*mass(ib)/(mass(ib)+mass(p))
        vel(ib,1)=vel(ib,1)-dvx*mass(p)/(mass(ib)+mass(p))
        vel(p,2)=vel(p,2)+dvy*mass(ib)/(mass(ib)+mass(p))
        vel(ib,2)=vel(ib,2)-dvy*mass(p)/(mass(ib)+mass(p))
        vel(p,3)=vel(p,3)+dvz*mass(ib)/(mass(ib)+mass(p))
        vel(ib,3)=vel(ib,3)-dvz*mass(p)/(mass(ib)+mass(p))
      endif
      if(mass(p).LE.0) call terror('danger starmass')
    enddo

    if(ioerror.EQ.0) close(sfrfilenr)
  endif

  if(verbosity.GT.0) print*,'<starform> unstable, forming,forced:', &
    nsph-nstable,nform,nforced

  nstar=nstar+nform
  nbodies=nbodies+nform
  totptag=totptag+nform
  tstarform=tnow
  deallocate(forced)
      
end subroutine
