! fromgabe
  function Hui_HII_recombB(T)
    real :: T
    real :: Hui_HII_recombB, lambda
    real, parameter :: T_HI=1.57807d5

    lambda = 2.d0 * T_HI / T
    Hui_HII_recombB = 2.753d-14 * lambda**1.500d0 / &
         ( 1.d0 + (lambda/2.740d0)**0.407d0 )**2.242d0

  end function Hui_HII_recombB


 subroutine projectparticles(mapmode)
 use MakeMapMod
 use StarsMod
 use ElementsMod
 use CoolingMod
 use fccoMod
  include 'globals.h'
  character(len=5) :: mapmode
  integer :: np,test,inear,i,p,bandnr,C
  real :: band(20),taufactor,clustersize,clustertime,rhoconst
  real :: hsearch,getlocalscale
  real :: pe,qc,fcoc,g0,zmetal,ethtoent
  real, external :: pRND
  real :: fha=0.3
  real, allocatable :: himass(:),gpos(:,:),hsm(:),opacs(:)
  integer, allocatable :: pindex(:)
  real, external :: Hui_HII_recombB

! set some constants 
 clustersize=0.04/unitl_in_kpc
 clustertime=10**7*year
 taufactor=6.25e-8*unitm_in_msun/unitl_in_kpc**2*fhydrogn*zQ()/solarzQ()
 rhoconst=6.7678e-32*unitm_in_msun/unitl_in_kpc**3
 zmetal=zQ()/solarzQ()
 C=elementnrQ("C   ")

 select case(mapmode)
 
 case('XXXXX')
  print*,'you have found the projectparticle debug mode!'
  print*,'N?'
  read*,np
 case('gas')
  np=nsph   
 case('dark')
  np=nbodies-nstar-nsph
 case('all')
  np=nbodies 
 case('HI')
  np=nsph
 case('H2')
  np=nsph
 case('Ha')
  np=nstar+nsph
 case('Ha2')
  np=nstar+nsph
 case('FUV')
  np=nstar+nsph
 case('stars')
  np=nstar
 case('new')
  np=nstar
 case('unst')
  np=nsph
 case('C')
  np=nsph
  call fccoInit
 case('CO')
  np=nsph   
  call fccoInit
 case('LC')
  np=nsph
  call fccoInit
 case('LCO')
  np=nsph   
  call fccoInit
 case('Ccool')
  np=nsph
 case default
  call startree
  bandnr=IACHAR(mapmode(1:1))-IACHAR('0')
  if(bandnr.LT.1.OR.bandnr.GT.nbandsQ()) &
    call terror('mapmode error')
  np=nsph+nstar
 end select

 allocate(gpos(np,1:3),hsm(np),himass(np),opacs(np),pindex(np),stat=test)
 if(test.NE.0) call terror('project mem. alloc. fails')

 select case(mapmode)
 
 case('XXXXX')
  print*,'give particle info: x,y,z,size,mass,opac'
  do i=1,np
   read*,gpos(i,1),gpos(i,2),gpos(i,3),hsm(i),himass(i),opacs(i) 
   pindex(i)=i
  enddo
 case('gas')
  np=nsph   
  do i=1,nsph
   himass(i)=unitm_in_msun*mass(i)
   gpos(i,1:3)=pos(i,1:3)
   hsm(i)=epsgrav(i)
   opacs(i)=0.
   pindex(i)=nbexist(i)
  enddo
 case('dark')
  np=nbodies-nstar-nsph
  do i=1,np
   p=nsph+i
   himass(i)=unitm_in_msun*mass(p)
   gpos(i,1:3)=pos(p,1:3)
   hsm(i)=epsgrav(p)
   opacs(i)=0.
   pindex(i)=nbexist(p)
  enddo
 case('all')
  np=nbodies 
  do i=1,nbodies
   himass(i)=unitm_in_msun*mass(i)
   gpos(i,1:3)=pos(i,1:3)
   hsm(i)=epsgrav(i)
   opacs(i)=0.
   pindex(i)=nbexist(i)
  enddo
 case('stars')
  np=nstar
  do i=1,nstar
   p=nbodies-nstar+i
   himass(i)=unitm_in_msun*mass(p)
   gpos(i,1:3)=pos(p,1:3)
   hsm(i)=epsgrav(p)
   pindex(i)=nbexist(p)
  enddo
 case('new')
  np=0
  do i=1,nstar-nbh
   p=nbodies-nstar+i
   if((tnow-tform(p))*timescale.LE.5.e6*year) then
    np=np+1
    himass(np)=unitm_in_msun*mass(p)
    gpos(np,1:3)=pos(p,1:3)
    hsm(np)=epsgrav(p)
    pindex(np)=nbexist(p)
   endif
  enddo
 case('HI')
  np=nsph
  do i=1,nsph
   himass(i)=unitm_in_msun*mass(i)*max(0.,1-elecfrac(i))*max(0.,1-h2frac(i))
   gpos(i,1:3)=pos(i,1:3)
   hsm(i)=hsmooth(i)
   opacs(i)=0.
   pindex(i)=nbexist(i)
  enddo
 case('H2')
   np=nsph
   do i=1,nsph
    himass(i)=unitm_in_msun*mass(i)*max(0.,1-elecfrac(i))*max(0.,h2frac(i))
    gpos(i,1:3)=pos(i,1:3)
    hsm(i)=hsmooth(i)
    opacs(i)=0.
    pindex(i)=nbexist(i)
   enddo

  case('Ha')
   np=nstar+nsph
   do i=1,nsph
    gpos(i,1:3)=pos(i,1:3)
    hsm(i)=hsmooth(i)
    himass(i)=0
    opacs(i)=0.75*taufactor* &
      mass(i)*max(0.,1-elecfrac(i)) !   *max(0.,1-h2frac(i))
    pindex(i)=nbexist(i)
   enddo
   do i=1,nstar
    p=nbodies-nstar+i
    himass(i+nsph)=unitm_in_msun*mass(p)*ha((tnow-tform(p))*timescale/year)
    hsm(i+nsph)=epsgrav(p)
    gpos(i+nsph,1:3)=pos(p,1:3)
    opacs(i+nsph)=0.
    pindex(i+nsph)=nbexist(p)
   enddo  

  case('Ha2')
   np=nstar+nsph
   do i=1,nsph
    temperat(i)=10000
!    print*,mass(i),elecfrac(i),rho(i),temperat(i),fha
    gpos(i,1:3)=pos(i,1:3)
    hsm(i)=hsmooth(i)
    himass(i)=unitm_in_msun*mass(i)*min(elecfrac(i),1.)* &
  &    elecfrac(i)*densconst*rho(i)*fha*Hui_HII_recombB(temperat(i))
    opacs(i)=0.75*taufactor* &
      mass(i)*max(0.,1-elecfrac(i)) ! *max(0.,1-h2frac(i))
    pindex(i)=nbexist(i)
   enddo
   do i=1,nstar
    p=nbodies-nstar+i
    himass(i+nsph)=0
!    himass(i+nsph)=unitm_in_msun*mass(p)*ha((tnow-tform(p))*timescale/year)
    hsm(i+nsph)=epsgrav(p)
    gpos(i+nsph,1:3)=pos(p,1:3)
    opacs(i+nsph)=0.
    pindex(i+nsph)=nbexist(p)
   enddo  

  case('FUV')
   np=nstar+nsph
   do i=1,nsph
    gpos(i,1:3)=pos(i,1:3)
    hsm(i)=hsmooth(i)
    himass(i)=0
    opacs(i)=3.*taufactor*mass(i)*max(0.,1-elecfrac(i)) ! *max(0.,1-h2frac(i))
    pindex(i)=nbexist(i)
   enddo
   do i=1,nstar
    p=nbodies-nstar+i
    himass(i+nsph)=unitm_in_msun*mass(p)*dFUV((tnow-tform(p))*timescale/year)
    gpos(i+nsph,1:3)=pos(p,1:3)
    hsm(i+nsph)=epsgrav(p)
    opacs(i+nsph)=0.
    pindex(i+nsph)=nbexist(p)
  enddo

  case('unst')
  np=0
  do i=1,nsph
   if(rho(i)*pi/6.*(csound(i)**2*pi/rho(i))**1.5.GT.masscrit) then
    np=np+1
    gpos(np,1:3)=pos(i,1:3)
    hsm(np)=hsmooth(i)
    himass(np)=unitm_in_msun*mass(i)*max(0.,1-elecfrac(i))*max(0.,1-h2frac(i))
    pindex(np)=nbexist(i)
   endif
  enddo
  
  case('C')
   np=nsph
   do i=1,nsph
    gpos(i,1:3)=pos(i,1:3)
    hsm(i)=hsmooth(i)
    if(uentropy) then 
      ethtoent=gamma1/rho(p)**gamma1
    else 
      ethtoent=1
    endif
    pe=(gamma1*ethermal(i)/ethtoent+vdisp(i)**2/3.)*velscale**2*rho(i)*rhoconst/kboltz
    g0=fuvheat(i)/heatconst/2/1.71  ! /2: half is blocked, /1.71: Habing -> Draine
    fcoc=fcco(rho(i)*densconst,G0,Pe,zmetal)
    qc=q1c(temperat(i))
    himass(i)=unitm_in_msun*mass(i)*max(0.,1-elecfrac(i))* & 
      min(fcoc,h2frac(i))
    pindex(i)=nbexist(i)
   enddo

  case('CO')
   np=nsph
   do i=1,nsph
    gpos(i,1:3)=pos(i,1:3)
    hsm(i)=hsmooth(i)    
    if(uentropy) then 
      ethtoent=gamma1/rho(p)**gamma1
    else 
      ethtoent=1
    endif
    pe=(gamma1*ethermal(i)/ethtoent+vdisp(i)**2/3.)*velscale**2*rho(i)*rhoconst/kboltz
    g0=fuvheat(i)/heatconst/2/1.71  ! /2: half is blocked, /1.71: Habing -> Draine
    fcoc=fcco(rho(i)*densconst,G0,Pe,zmetal)
    qc= q1co(temperat(i))
    himass(i)=unitm_in_msun*mass(i)*max(0.,1-elecfrac(i))* &
      min(fcoc,h2frac(i))
    pindex(i)=nbexist(i)
   enddo
   
  case('LC')
   np=nsph
   do i=1,nsph
    gpos(i,1:3)=pos(i,1:3)
    hsm(i)=hsmooth(i)
    if(uentropy) then 
      ethtoent=gamma1/rho(p)**gamma1
    else 
      ethtoent=1
    endif
    pe=(gamma1*ethermal(i)/ethtoent+vdisp(i)**2/3.)*velscale**2*rho(i)*rhoconst/kboltz
    g0=fuvheat(i)/heatconst/2/1.71  ! /2: half is blocked, /1.71: Habing -> Draine
    fcoc=fcco(rho(i)*densconst,G0,Pe,zmetal)
    qc=q1c(temperat(i))
    himass(i)=unitm_in_msun*mass(i)*max(0.,1-elecfrac(i))* & 
      min(fcoc,h2frac(i))*qc*zmetal
    pindex(i)=nbexist(i)
   enddo

  case('LCO')
   np=nsph
   do i=1,nsph
    gpos(i,1:3)=pos(i,1:3)
    hsm(i)=hsmooth(i)    
    if(uentropy) then 
      ethtoent=gamma1/rho(p)**gamma1
    else 
      ethtoent=1
    endif
    pe=(gamma1*ethermal(i)/ethtoent+vdisp(i)**2/3.)*velscale**2*rho(i)*rhoconst/kboltz
    g0=fuvheat(i)/heatconst/2/1.71  ! /2: half is blocked, /1.71: Habing -> Draine
    fcoc=fcco(rho(i)*densconst,G0,Pe,zmetal)
    qc= q1co(temperat(i))
    himass(i)=unitm_in_msun*mass(i)*max(0.,1-elecfrac(i))* &
      min(fcoc,h2frac(i))*qc*zmetal
    pindex(i)=nbexist(i)
   enddo
   
  case('Ccool')
   np=nsph
   do i=1,nsph
    gpos(i,1:3)=pos(i,1:3)
    hsm(i)=hsmooth(i)    
    himass(i)=unitm_in_msun*mass(i)*cool_par*rho(i)* &
      ElementCoolFunc(elecfrac(i),temperat(i),C)
    pindex(i)=nbexist(i)   
   enddo
   print*,'lum:',sum(himass(1:nsph))
    
  case default
   np=nsph+nstar
   do i=1,nsph
    gpos(i,1:3)=pos(i,1:3)
    hsm(i)=hsmooth(i)
    himass(i)=0
    opacs(i)=axavQ(bandnr)*taufactor* &
      mass(i)*max(0.,1-elecfrac(i)) ! *max(0.,1-h2frac(i))   
    pindex(i)=nbexist(i)   
   enddo
   do i=1,nstar
    p=nbodies-nstar+i
    band=mbands((tnow-tform(p))*timescale/year)
    gpos(i+nsph,1:3)=pos(p,1:3)
    himass(i+nsph)=unitm_in_msun*mass(p)*band(bandnr)
    hsearch=getlocalscale(pos(p,1:3))
    call nnscalefix(pos(p,1:3),hsearch,inear,targetnn,nn_tol,root)
!    hsearch=0.05
!    inear=0
!    call search(root,2*hsearch,pos(p,1:3),inear,bodlist)
!    print*,p,hsearch,inear
!    hsm(nsph+i)= &
!      MIN(hsearch,clustersize*(1+((tnow-tform(p))*timescale/clustertime)))   
!    hsearch=eps
!    hsm(nsph+i)=0.05-hsearch*log( MAX(pRND(nbexist(p)+4),0.00001) )
    hsm(nsph+i)=-hsearch*log( MAX(pRND(nbexist(p)+4),0.00001) )/2. &
                -hsearch*log( MAX(pRND(3*nbexist(p)+4),0.00001) )/2.
!    hsm(nsph+i)=hsearch*2*pRND(nbexist(p)+4)
    pindex(nsph+i)=nbexist(p)
    opacs(nsph+i)=0.
   enddo  
  end select

  if(np.GT.0) call project(np,gpos,hsm,himass,opacs,pindex)

  deallocate(gpos,hsm,himass,opacs,pindex)

 end subroutine



subroutine maphimap(nmap)
 use MakeMapMod
 use StarsMod
 use ElementsMod
 include 'globals.h'
 integer nmap
 integer,parameter :: ndigits=6,maxmap=100
 character(len=ndigits) :: nstring
 character(len=80) :: filenaam,opnaam
 character(len=8) :: identify(maxmap)
 character(len=5) :: mapmode(maxmap)
 integer :: imsize(maxmap,2),proj(maxmap),ext(maxmap)
 integer, save :: maps=0 
 real  :: focus(maxmap,3), pointing(maxmap,3), &
           upvector(maxmap,3),width(maxmap),zmin(maxmap)
 integer :: ioerror,i,j,skipfactor,oldseed
 
! reset rndtable to fixed seed (to get consistency for random stars) 
 oldseed=rnseed
 rnseed=initialseed
 call setRND
 rnseed=oldseed

 skipfactor=1  
 call itos(nmap,ndigits,nstring)
  
! read image info (if present)
 open(unit=upars, file='image', status='OLD', iostat=ioerror)
 if(ioerror.NE.0) RETURN

 read(upars,*,iostat=ioerror) maps, skipfactor
   if(skipfactor.LE.0) skipfactor=1
   if(maps.GT.maxmap) maps=maxmap
   if(maps.GT.0) then
    do i=1,maps
    read(upars,*,iostat=ioerror) mapmode(i)
    read(upars,*,iostat=ioerror) identify(i)
    read(upars,*,iostat=ioerror) imsize(i,1:2)
    read(upars,*,iostat=ioerror) focus(i,1:3)
    read(upars,*,iostat=ioerror) pointing(i,1:3)
    read(upars,*,iostat=ioerror) upvector(i,1:3)
    read(upars,*,iostat=ioerror) width(i)
    read(upars,*,iostat=ioerror) proj(i)
    read(upars,*,iostat=ioerror) ext(i)
    read(upars,*,iostat=ioerror) zmin(i)
    if(width(i).eq.0.and.proj(i).eq.0) width(i)=rsize
    if(width(i).eq.0.and.proj(i).eq.1) width(i)=45.
    enddo
   endif
   if(ioerror.NE.0) then
    print*,' stop -- error reading image info'
    stop
   endif
  close(upars)
  
 if(mod(nmap,skipfactor).NE.0) return 

 do j=1,maps
  call Initmap(imsize(j,1:2),width(j),focus(j,1:3),pointing(j,1:3),upvector(j,1:3),proj(j),ext(j),zmin(j))

  if(verbosity.GT.0) print*,'<map> ',identify(j),imsize(j,1:2)

  call projectparticles(mapmode(j))

  filenaam=trim(outputfile)//'-'//nstring//'_'//trim(identify(j))//'.fits'
  opnaam='XXXXX'
  if(ext(j).EQ.1.AND.j.EQ.maps) opnaam='opacity'
  call EndMap(filenaam,opnaam)
 enddo
end subroutine






