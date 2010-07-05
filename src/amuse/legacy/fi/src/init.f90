subroutine initdump(n)
  use fixedhalomod
  use ElementsMod
  use CoolingMod
  use IonMod
  use MoleculeMod
  use H2CoolMod
  use StarsMod
  include 'globals.h' 
  integer n

  call readdump(n)
  call set_parameters(1)
  call initfixedhalo(poshalo,trim(datadir)//halofile)
  if(masshalo.NE.fixedhalomass) call terror('halo reinit error')

  call heattabel

  call InitStars(datadir,zQ())

  if(usepm) then
    if(verbosity.GT.0) print*,' ...initPM...'
    call initpm 
  endif
  call outstate(n)      
end subroutine

subroutine initpm
  use pmgravMod
  include 'globals.h' 
  integer i
  real r,bmin(3),bmax(3)
   
  bmin=-hboxsize
  bmax=hboxsize 
  call pmgravInit(npm, bmin, bmax)
   
  tpm=tnow
  rcut=rcutQ()
  rcut2=rcut**2
  pmsoft=softQ()
  pmdr=ninterp/rcut2      
  do i=0,ninterp
    r=1/pmsoft*sqrt(real(i)/pmdr)
    pmpot(i)=ficorrect(r)
    pmacc(i)=forcecorrect(r)
  enddo
  pmpot(ninterp+1)=0 
  pmacc(ninterp+1)=0 
  if(verbosity.GT.0) write(*,'("  > pm box, rcut, rsoft:", 3f8.2)') &
    pboxsize,sqrt(rcut2),pmsoft
end subroutine

subroutine inithalo
  use fixedhalomod                 
  include 'globals.h'
  
  if(halofile.EQ.'none'.OR.halofile.EQ.'NONE'.OR. &
    halofile.EQ.'None') then
    fixthalo=.FALSE.
    return
  endif

  print*,' ...dark halo...'
  poshalo=0.
  call initfixedhalo(poshalo,trim(datadir)//halofile)
  masshalo=fixedhalomass
end subroutine

subroutine initeps
  include 'globals.h'
  integer i,j,low
  real epsavg
  if(adaptive_eps.OR.eps.EQ.0) then
    if(usesph) then
      call seteps('coll')
    else
      call seteps('all ')
    endif

    if(adaptive_eps) then
      return
    else 
      epsavg=0.
      low=1
      if(usesph.AND.nbodies.GT.nsph) low=1+nsph
      do i=low,nbodies
        epsavg=epsavg+epsgrav(i)
      enddo
      epsavg=epsavg/(nbodies-low+1)
    endif
  else
    epsavg=eps
  endif
  
  if(eps.EQ.0) eps=epsavg  ! eps in snapshot      
  epsgrav(nsph+1:nbodies)=epsavg

  if(.not.usesph) epsgrav(1:nsph)=epsavg

end subroutine

subroutine initrnd
  include 'globals.h'   
  rnseed=initialseed
  call setRND
end subroutine

subroutine initpars
  use ElementsMod
  include 'globals.h'
  integer p,i
  real deldr2,xw,xw2,deldrg,xw3,xw4

  if(verbosity.GT.0) print*,' ...initpars...'

  nsnap=firstsnap

  call initrnd

  radiate=radiate.AND.usesph
  starform=starform.AND.usesph

  tiny=1.e-20

  gamma1=gamma-1.
  if(gamma.EQ.1.) gamma1=1.

  heat_par1=heat_par1*3.46e14*unitl_in_kpc**2.5/unitm_in_msun**1.5
  heat_par2=heat_par2*2.34e-17/SQRT(unitl_in_kpc*unitm_in_msun)
  cool_par=cool_par*2.34e-17/SQRT(unitl_in_kpc*unitm_in_msun)
  mhboltz=1.211e-8
  mhboltz=mhboltz*4.30e4*unitm_in_msun/unitl_in_kpc
!  mhboltz=mhydrogen/kboltz*4.30e4*unitm_in_msun/unitl_in_kpc
  mumhkbol=meanmwt*mhboltz
  mumhkgam=meanmwt*mhboltz/gamma
  mumhkg1=meanmwt*mhboltz*gamma1
  if(radiate) graineff=graineff*zQ()/solarzQ()     
  
  heatconst=2.0693e14*unitl_in_kpc**2.5/unitm_in_msun**1.5
!  nsconst: rho-> H nr dens  (strictly speaking fhydrogn = f_H/meanmwt)
  densconst=4.078e-8*fhydrogn*unitm_in_msun/unitl_in_kpc**3
  timescale=1.488e19*unitl_in_kpc*sqrt(unitl_in_kpc/unitm_in_msun)
  lengthscale=3.086e21*unitl_in_kpc
  velscale=lengthscale/timescale
  flxscale=-41.282-2.*alog10(unitl_in_kpc)
  massscale=1.989e33*unitm_in_msun

  tsnbeg=tsnbeg*year/timescale
  tbubble=tbubble*year/timescale
!  snenergy=nsupernova*esupernova*unitl_in_kpc/(unitm_in_msun*unitm_in_msun)
  snenergy=unitm_in_msun*nsupernova* &
    (esupernova/1.e51)*1.1695e13*unitl_in_kpc/(unitm_in_msun*unitm_in_msun)
 
  tcollfac=tcollfac*sfeff
  masscrit=masscrit/unitm_in_msun
     
  hboxsize=0.5*pboxsize

  if(periodic) then
    call terror('periodic tb fixed ')
  endif

  deldr2=4./ninterp
  deldr2i=1./deldr2

  do i=0,1+ninterp
    xw2=i*deldr2
    xw=SQRT(xw2)
    if(xw2.LE.1.) then
      wsmooth(i)=1.-1.5*xw2+0.75*xw*xw2
      dwsmooth(i)=-3.+2.25*xw
    else
      wsmooth(i)=0.25*(2.-xw)**3
      dwsmooth(i)=-0.75*(2.-xw)**2/(xw+tiny)
    endif
    if(xw2.GE.4) then
      wsmooth(i)=0.
      dwsmooth(i)=0.
    endif
  enddo

  deldrg=2./ninterp

  if(radiate) call rdintpol

end subroutine

subroutine initpos
  include 'globals.h'
  call corrpos(itimestp,'desync')
end subroutine

subroutine smoothall
  include 'globals.h'
  real temp(nsphmax)
  integer i,j

  call terror(' TBD smoothing')

  do i=1,3
    do j=1,nsph
      temp(j)=vel(j,i)
    enddo
!    call omp_smooth(temp)
    do j=1,nsph
      veltpos(j,i)=vel(j,i)
      vel(j,i)=temp(j)
    enddo
  enddo

  do j=1,nsph
    temp(j)=ethermal(j)
  enddo
!  call omp_smooth(temp)
  do j=1,nsph
    ethold(j)=ethermal(j)
    ethermal(j)=temp(j) 
  enddo  
    
end subroutine

subroutine initstep
  include 'globals.h'
  integer i

  do i=1,nbodies
    itimestp(i)=1
    otimestp(i)=1
  enddo

!  upbin=max_tbin
  upbin=0
  stime=0.0
  endstep=.TRUE.
  npactive=nbodies
  tsteppos=0.

  do i=1,npactive
    pactive(i)=i
  enddo

  if(usesph) then
    nsphact=nsph
  else
    nsphact=0
  endif

end subroutine

subroutine initsys
  include 'globals.h'
  integer ij
  real rtime

  call set_parameters(0)

  call startout

  call heattabel

  call initpars

  call inithalo

  if(periodic) call initbc

  if(radiate) call initnewstar

  call readbods(inputfile)

  if(usepm) then
    if(verbosity.GT.0) PRINT*,' ...initPM...'
    call initpm 
  endif
          
  call postprocessread

  call initpos
                          
  call outstate(0)

end subroutine

subroutine ethdotinit   
  include'globals.h'
   
  if(sph_visc.EQ.'bulk') call terror(' bulk visc. not implemented ')
           
  if(.NOT.uentropy) then          
    call alldethacc
  else
    call alldentacc
  endif
          
end subroutine

subroutine initepsgas
  include 'globals.h' 
  integer i
  real hsmthavg

  if(epsgas.GT.0) then
    do i=1,nsph
      epsgrav(i)=epsgas
    enddo
    return
  endif
  if(eps.GT.0) then
    do i=1,nsph
      epsgrav(i)=eps
    enddo
    return
  endif
  call seteps('sph ')
  hsmthavg=0.
  do i=1,nsph
    hsmthavg=hsmthavg+epsgrav(i)
  enddo
  hsmthavg=hsmthavg/nsph
  do i=1,nsph
    epsgrav(i)=hsmthavg
  enddo    
end subroutine

subroutine rdintpol 
  include 'globals.h'
  integer i,ioerror
  real dummy
!  waardes van de spline soft. (als functie van (r/eps)^2, voor 0+1/ninter2-4)

!  lees uvpot in: spline soft tabel voor uv flux
!  nb: uv interpolatie gaat tot grotere r!

  open(UNIT=uvfilenr,file=trim(datadir)//uvfile, status='OLD',iostat=ioerror)
  if(ioerror.NE.0) then
    print*,' stop -- data file error:',uvfile       
    print*,' please check directory:',datadir
  endif   
  do i=1,uvinterp
    read(uvfilenr,*)dummy,smoothuv(i)
  enddo        
  close(uvfilenr)

  smoothuv(0)=1
  do i=0,uvinterp
    smoothuv(i)=smoothuv(i)*3       
  enddo
end subroutine

subroutine nbexistinit
  include 'globals.h'
  integer i
   
  totptag=nbodies         
  do i=1,nbodies
    nbexist(i)=i
  enddo  
end subroutine

subroutine initbc
  call terror('not implemented')
end subroutine
