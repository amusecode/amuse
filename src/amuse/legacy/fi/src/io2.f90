subroutine outbods(fname)
 include 'globals.h'
 integer,parameter :: ndigits=6
 character(len=ndigits) :: nstring
 character(len=24) :: filenaam
 character(len=*) :: fname
 if(fname(1:6).EQ.'XXXXXX') then
  call itos(nsnap,ndigits,nstring)
  filenaam=TRIM(outputfile)//'.'//nstring
  nsnap=nsnap+1 
 else
  filenaam=fname
 endif 
  call writebods(filenaam)
end subroutine outbods


subroutine writebods(filenaam)
 include 'globals.h'
 character(len=*) :: filenaam 
 integer :: nbods,nkeep,p
 integer(kind=4), parameter :: nihead=32,nrhead=32,nphead=64
 integer(kind=4) :: ioversion
 integer(kind=4) :: header(8),ihead(nihead)
 real(kind=8) :: rhead(nrhead)
 integer(kind=4) :: phead(nphead)
 
 call prepareoutput
 
 ioversion=1
  
 header(1)=nihead
 header(2)=nrhead
 header(3)=nphead
 header(4)=0
 header(5)=0
 header(6)=0
 header(7)=0
 header(8)=0
 
 ihead(1:nihead)=0
 rhead(1:nrhead)=0
 phead(1:nphead)=0
  
 ihead(1)=nbodies
 ihead(2)=nsph
 ihead(3)=nstar
 ihead(4)=nbodies-nstar-nsph
 ihead(5)=totptag
 ihead(6)=nbh

 ihead(8)=syncflag
 if(syncflag.NE.0) then
  print*, ' ** warning ** async output??'
  stop 
 endif

 
 rhead(1)=tnow
 rhead(2)=mtot
 rhead(3)=ektot
 rhead(4)=eptot
 rhead(5)=snheat
 rhead(6)=eradiate
 rhead(7)=esofttot
 rhead(8)=enerror
 rhead(9)=estar
 rhead(10)=efuvheat
 rhead(11)=eradcool
 rhead(12)=tpos
 rhead(13)=teth
 rhead(14)=trad
 rhead(15)=meanmwt
 rhead(16)=timescale
 rhead(17)=unitm_in_msun
 rhead(18)=lengthscale
 rhead(19)=densconst
 rhead(20)=massres
 rhead(21)=fhydrogn
 rhead(22)=tstarform
 rhead(23)=tsnfeedback
 rhead(24)=tbh
 rhead(25)=eps

! rhead(31)=redshift
! rhead(32)=redshift0
 
 phead(1:nphead)=output(1:nphead) 
 if(nsph.eq.0) then
  phead(10:25)=0
 endif
 if(nstar.eq.0) then
  phead(33:35)=0
 endif
  
 open(unit=ubodsout,FILE=filenaam,status='UNKNOWN',form='UNFORMATTED')
 
 WRITE(ubodsout) ioversion
 WRITE(ubodsout) (header(p),p=1,8)
 WRITE(ubodsout) (ihead(p),p=1,nihead)
 WRITE(ubodsout) (rhead(p),p=1,nrhead)
 WRITE(ubodsout) (phead(p),p=1,nphead)

 IF(phead(1).EQ.1) WRITE(ubodsout) (mass(p),p=1,nbodies)
 IF(phead(2).EQ.1) WRITE(ubodsout) (pos(p,1),p=1,nbodies)
 IF(phead(2).EQ.1) WRITE(ubodsout) (pos(p,2),p=1,nbodies)
 IF(phead(2).EQ.1) WRITE(ubodsout) (pos(p,3),p=1,nbodies)
 IF(phead(3).EQ.1) WRITE(ubodsout) (vel(p,1),p=1,nbodies)
 IF(phead(3).EQ.1) WRITE(ubodsout) (vel(p,2),p=1,nbodies)
 IF(phead(3).EQ.1) WRITE(ubodsout) (vel(p,3),p=1,nbodies)
 IF(phead(4).EQ.1) WRITE(ubodsout) (epsgrav(p),p=1,nbodies)
 IF(phead(5).EQ.1) WRITE(ubodsout) (tform(p),p=1,nbodies)
 IF(phead(6).EQ.1) WRITE(ubodsout) (acc(p,1),p=1,nbodies)
 IF(phead(6).EQ.1) WRITE(ubodsout) (acc(p,2),p=1,nbodies)
 IF(phead(6).EQ.1) WRITE(ubodsout) (acc(p,3),p=1,nbodies)

IF(nsph.GT.0) THEN
 IF(phead(10).EQ.1) WRITE(ubodsout) (rho(p),p=1,nsph)
 IF(phead(11).EQ.1) THEN
  IF(uentropy) THEN
   WRITE(ubodsout) (tempvect(p),p=1,nsph)
  ELSE
   WRITE(ubodsout) (ethermal(p),p=1,nsph)  
  ENDIF
 ENDIF
 IF(phead(12).EQ.1) THEN
  IF(uentropy) THEN
   WRITE(ubodsout) (entropy(p),p=1,nsph)
  ELSE
   WRITE(ubodsout) (tempvect(p),p=1,nsph)  
  ENDIF
 ENDIF
 IF(phead(13).EQ.1) WRITE(ubodsout) (hsmooth(p),p=1,nsph) 
 IF(phead(14).EQ.1) WRITE(ubodsout) (fuvheat(p),p=1,nsph)
 IF(phead(15).EQ.1) WRITE(ubodsout) (esnthdt(p),p=1,nsph)
 IF(phead(16).EQ.1) WRITE(ubodsout) (tcollaps(p),p=1,nsph)
 IF(phead(17).EQ.1) WRITE(ubodsout) (temperat(p),p=1,nsph)
 IF(phead(18).EQ.1) WRITE(ubodsout) (elecfrac(p),p=1,nsph)
 IF(phead(19).EQ.1) WRITE(ubodsout) (csound(p),p=1,nsph)
 IF(phead(20).EQ.1) WRITE(ubodsout) &
  (csound(p)**2*rho(p)/gamma,p=1,nsph)
 IF(phead(21).EQ.1) WRITE(ubodsout) (hsmdivv(p),p=1,nsph)
 IF(phead(22).EQ.1) WRITE(ubodsout) (mumaxdvh(p),p=1,nsph) 
 IF(phead(23).EQ.1) WRITE(ubodsout) (hsmcurlv(p),p=1,nsph)
 IF(phead(24).EQ.1) WRITE(ubodsout) (vdisp(p),p=1,nsph)
 IF(phead(25).EQ.1) WRITE(ubodsout) (h2frac(p),p=1,nsph)
ENDIF

IF(nstar.GT.0) THEN 
 IF(phead(34).EQ.1) WRITE(ubodsout) (starfuv(p),p=nbodies-nstar+1,nbodies)
 IF(phead(35).EQ.1) WRITE(ubodsout) (snentropy(p),p=nbodies-nstar+1,nbodies)
ENDIF

 IF(phead(40).EQ.1) WRITE(ubodsout) (phi(p),p=1,nbodies)
 IF(phead(41).EQ.1) WRITE(ubodsout) (phiext(p),p=1,nbodies)
 IF(phead(42).EQ.1) WRITE(ubodsout) (nbexist(p),p=1,nbodies)
 IF(phead(43).EQ.1) WRITE(ubodsout) (itimestp(p),p=1,nbodies)

 close(ubodsout)

end subroutine writebods

subroutine readbods(filenaam)
 include 'globals.h'
 character(len=*) :: filenaam 
 integer :: nbods,nkeep,p,ioerr
 integer(kind=4),parameter ::  nihead=32,nrhead=32,nphead=64
 integer(kind=4) :: ioversion
 integer(kind=4) :: header(8),ihead(nihead)
 real(kind=8) :: rhead(nrhead)
 integer(kind=4) :: phead(nphead)
 
 print*,' ...reading data... ',filenaam
 open(unit=ubodsin,FILE=filenaam,status='OLD',form='UNFORMATTED',IOSTAT=ioerr)
 
 IF(ioerr.NE.0) THEN
  print*, ' stop -- error reading file ', filenaam
  stop         
 ENDIF
  
 READ(ubodsin) ioversion
 READ(ubodsin) (header(p),p=1,8)
 
 ihead(1:nihead)=0
 rhead(1:nrhead)=0
 phead(1:nphead)=0
 

 if(header(1).GT.nihead.OR.header(2).GT.nrhead.OR. &
                 header(3).GT.nphead.OR.ioversion.NE.1) then
  print*,' warning: file input from newer file format?'
 endif
 
 READ(ubodsin) (ihead(p),p=1,MIN(nihead,header(1)))
 READ(ubodsin) (rhead(p),p=1,MIN(nrhead,header(2)))
 READ(ubodsin) (phead(p),p=1,MIN(nphead,header(3)))

 nbodies=ihead(1)
 nsph=ihead(2)
 nstar=ihead(3)
 IF(nbodies-nstar-nsph.NE.ihead(4)) THEN
  print*,' inconsistent body count'
  STOP
 ENDIF  
 IF(nbodies.GT.nbodsmax.OR.nsph.gt.nsphmax) THEN
  PRINT*,nbodies,nbodsmax,nsph,nsphmax
  PRINT*,' particle overflow '  
  PRINT*,'(increase compiled limits)'
  STOP
 ENDIF
 IF(nbodies-nstar-nsph.GT.0.AND.fixthalo) THEN
  PRINT*,'** warning ** rigid halo used and DM particles found' 
 ENDIF
 totptag=ihead(5)
 nbh=ihead(6)
 
 syncflag=ihead(8)
 if(syncflag.NE.0) print*,' ** warning ** async input'
  
 tnow=rhead(1)
 mtot=rhead(2)
 ektot=rhead(3)
 eptot=rhead(4)
 snheat=rhead(5)
 eradiate=rhead(6)
 esofttot=rhead(7)
 enerror=rhead(8)
 estar=rhead(9)
 efuvheat=rhead(10)
 eradcool=rhead(11)
 tpos=rhead(12)
 teth=rhead(13)
 trad=rhead(14)
 if(meanmwt.LE.0) meanmwt=rhead(15)
 massres=rhead(20)
 ! tstarform=rhead(22)
 tstarform=tnow
 tsnfeedback=tnow
 tbh=tnow
 if(rhead(25).NE.0) eps=rhead(25) 

 input(1:nphead)=phead(1:nphead)

 pordercount=pordercount+1
 ppropcount=ppropcount+1

 IF(phead(1).EQ.1) READ(ubodsin) (mass(p),p=1,nbodies)
 IF(phead(2).EQ.1) READ(ubodsin) (pos(p,1),p=1,nbodies)
 IF(phead(2).EQ.1) READ(ubodsin) (pos(p,2),p=1,nbodies)
 IF(phead(2).EQ.1) READ(ubodsin) (pos(p,3),p=1,nbodies)
 IF(phead(3).EQ.1) READ(ubodsin) (vel(p,1),p=1,nbodies)
 IF(phead(3).EQ.1) READ(ubodsin) (vel(p,2),p=1,nbodies)
 IF(phead(3).EQ.1) READ(ubodsin) (vel(p,3),p=1,nbodies)
 IF(phead(4).EQ.1) READ(ubodsin) (epsgrav(p),p=1,nbodies)
 IF(phead(5).EQ.1) READ(ubodsin) (tform(p),p=1,nbodies)
 IF(phead(6).EQ.1) READ(ubodsin) (acc(p,1),p=1,nbodies)
 IF(phead(6).EQ.1) READ(ubodsin) (acc(p,2),p=1,nbodies)
 IF(phead(6).EQ.1) READ(ubodsin) (acc(p,3),p=1,nbodies)


IF(nsph.GT.0) THEN
 IF(phead(10).EQ.1) READ(ubodsin) (rho(p),p=1,nsph)
 IF(phead(11).EQ.1) THEN
  IF(uentropy) THEN
   READ(ubodsin) (tempvect(p),p=1,nsph)
  ELSE
   READ(ubodsin) (ethermal(p),p=1,nsph)  
  ENDIF
 ENDIF
 IF(phead(12).EQ.1) THEN
  IF(uentropy) THEN
   READ(ubodsin) (entropy(p),p=1,nsph)
  ELSE
   READ(ubodsin) (tempvect(p),p=1,nsph)  
  ENDIF
 ENDIF
 IF(phead(13).EQ.1) READ(ubodsin) (hsmooth(p),p=1,nsph) 
 IF(phead(14).EQ.1) READ(ubodsin) (fuvheat(p),p=1,nsph)
 IF(phead(15).EQ.1) READ(ubodsin) (esnthdt(p),p=1,nsph)
 IF(phead(16).EQ.1) READ(ubodsin) (tcollaps(p),p=1,nsph)
 IF(phead(17).EQ.1) READ(ubodsin) (temperat(p),p=1,nsph)
 IF(phead(18).EQ.1) READ(ubodsin) (elecfrac(p),p=1,nsph)
 IF(phead(19).EQ.1) READ(ubodsin) (csound(p),p=1,nsph)
 IF(phead(20).EQ.1) READ(ubodsin) (csound(p),p=1,nsph)
 IF(phead(21).EQ.1) READ(ubodsin) (hsmdivv(p),p=1,nsph)
 IF(phead(22).EQ.1) READ(ubodsin) (mumaxdvh(p),p=1,nsph)
 IF(phead(23).EQ.1) READ(ubodsin) (hsmcurlv(p),p=1,nsph)
 IF(phead(24).EQ.1) READ(ubodsin) (vdisp(p),p=1,nsph)
 IF(phead(25).EQ.1) READ(ubodsin) (h2frac(p),p=1,nsph)
 ENDIF

IF(nstar.GT.0) THEN
 IF(phead(34).EQ.1) READ(ubodsin) (starfuv(p),p=nbodies-nstar+1,nbodies)
 IF(phead(35).EQ.1) READ(ubodsin) (snentropy(p),p=nbodies-nstar+1,nbodies)
ENDIF
 
 IF(phead(40).EQ.1) READ(ubodsin) (phi(p),p=1,nbodies)
 IF(phead(41).EQ.1) READ(ubodsin) (phiext(p),p=1,nbodies)
 IF(phead(42).EQ.1) READ(ubodsin) (nbexist(p),p=1,nbodies)
 IF(phead(43).EQ.1) READ(ubodsin) (itimestp(p),p=1,nbodies)

 close(ubodsin)
 
end subroutine readbods


subroutine prepareoutput
 include 'globals.h'
 integer :: nkeep,nbods
 if(uentropy) then
  tempvect(1:nsph)=entropy(1:nsph)/gamma1*rho(1:nsph)**gamma1
 else
  if(.NOT.isotherm) then
   tempvect(1:nsph)=ethermal(1:nsph)*gamma1/rho(1:nsph)**gamma1 
  else
   ethermal(1:nsph)=csound(1:nsph)**2
   tempvect(1:nsph)=csound(1:nsph)**2
   temperat(1:nsph)=csound(1:nsph)**2*meanmwt*mhboltz
  endif
 endif 
end subroutine

subroutine itos(i,n,s)
 integer,intent(in) :: i,n
 character(len=n), intent(inout) :: s
 character(len=11) :: nstring
 data nstring/'0123456789X'/
 integer :: j,k,l

 if(i.LT.0.OR.i.GE.10**n) then
  do k=1,n
  s(k:k)=nstring(11:11)
  enddo
  return
 endif 
 j=1
 do k=1,n
 l=1+mod(i,10*j)/j
 s(n-k+1:n-k+1)=nstring(l:l)
 j=j*10
 enddo

end subroutine

! insert autmatically generated code below

subroutine writedump(n) 
 include 'globals.h' 
 integer i,n 
 real dummy 
 logical ldummy
 open(unit=uboddump,file=dumpfile,status='UNKNOWN',form='UNFORMATTED') 

write(uboddump) n 

      write(uboddump)  output,firstsnap,stepout,steplog,verbosity,       &
   datadir,inputfile,outputfile,halofile
      write(uboddump)  pboxsize,usesph,radiate,starform,cosmo
      write(uboddump)  unitl_in_kpc,unitm_in_msun
      write(uboddump)  dtime,tstepcrit,tstpcr2,freev,freea,freevexp,   &
   freeaexp,nsteps,max_tbin,minppbin,sqrttstp,acc_tstp, &
   freetstp
      write(uboddump)  bh_tol,eps,gdgtol,nn_tol,targetnn,usequad,         &
   directsum,selfgrav,fixthalo,adaptive_eps,gdgop
      write(uboddump)  epsgas,gamma,alpha,beta,epssph,courant,removgas, &
   consthsm,nsmtol,nsmooth,smoothinput,consph,sphinit,uentropy, &
   isotherm,eps_is_h,hupdatemethod,sph_visc
      write(uboddump)  graineff,crionrate,heat_par1,heat_par2,cool_par,    &
   optdepth
      write(uboddump)  comove
      write(uboddump)  tcollfac,masscrit,sfeff,tbubble,sne_eff,tsnbeg,   &
   rhomax,sfmode,feedback
      write(uboddump)  rndtable,rnseed
      write(uboddump)  nbodies,nsnap
      write(uboddump)  rmin,rsize,incells,incellsg
      write(uboddump)  root
      write(uboddump)  nttot,ntmin,ntmax,ntavg,nttotfuv,ntminfuv,        &
   ntmaxfuv,ntavgfuv
      write(uboddump)  nstot,nsmin,nsmax,nsavg
      write(uboddump)  tnow,tpos
      write(uboddump)  tiny
      write(uboddump)  mtot,etot,ektot,eptot,mstar,mgas,snheat,esofttot, &
   enerror,amvec,cmpos,cmvel
      write(uboddump)  input
      write(uboddump)  npactive,nsphact
      write(uboddump)  etol,dummy,tsteppos,active_bin
      write(uboddump)  syncflag,entropyflag
      write(uboddump)  teth,ethtot
      write(uboddump)  symmetry
      write(uboddump)  dummy,dummy,gamma1,nstar,nsph
      write(uboddump)  massres
      write(uboddump)  deldr2i
      write(uboddump)  wsmooth,dwsmooth
      write(uboddump)  nntot,nnmin,nnmax,nnavg
      write(uboddump)  hboxsize
      write(uboddump)  poshalo,masshalo
      write(uboddump)  eradiate,trad,meanmwt,     &
   mhboltz,fhydrogn,mumhkbol,mumhkgam, &
   mumhkg1,efuvheat,eradcool,mcold,mluke,mwarm,mhot,estar, &
   heatconst,densconst,timescale,lengthscale,velscale,flxscale, &
   massscale
      write(uboddump)  tstarform,tsnfeedback,snenergy,totptag
      write(uboddump)  h2time
      write(uboddump)  smoothuv
      write(uboddump)  tbh,nbh
      write(uboddump)  tpm,rcut,rcut2,pmsoft,pmdr,pmpot,       &
   pmacc
      write(uboddump)  searchpos,searchh,searchdelta,searchacc4,     &
   searchn,searchreuse,reuseflag,ncalls,nsearches
 write(uboddump) ( itimestp(i),i=1,nbodies) 
 write(uboddump) ( nbexist(i),i=1,nbodies) 
 write(uboddump) ( order_bodlist(i),i=1,nbodies) 
 write(uboddump) ( otimestp(i),i=1,nbodies) 
 write(uboddump) ( pactive(i),i=1,nbodies) 
 write(uboddump) ( acc(i,1:ndim+1),i=1,nbodies) 
 write(uboddump) ( csound(i),i=1,nsph) 
 write(uboddump) ( derad(i),i=1,nsph) 
 write(uboddump) ( dethdt(i),i=1,nsph) 
 write(uboddump) ( dethold(i),i=1,nsph) 
 write(uboddump) ( drhodh(i),i=1,nsph) 
 write(uboddump) ( elecfrac(i),i=1,nsph) 
 write(uboddump) ( epsgrav(i),i=1,nbodies) 
 write(uboddump) ( esnthdt(i),i=1,nsph) 
 write(uboddump) ( ethermal(i),i=1,nsph) 
 write(uboddump) ( ethold(i),i=1,nsph) 
 write(uboddump) ( fuvheat(i),i=1,nbodies) 
 write(uboddump) ( h2frac(i),i=1,nsph) 
 write(uboddump) ( hsmcurlv(i),i=1,nsph) 
 write(uboddump) ( hsmdivv(i),i=1,nsph) 
 write(uboddump) ( hsmooth(i),i=1,nbodies) 
 write(uboddump) ( mass(i),i=1,nbodies) 
 write(uboddump) ( mumaxdvh(i),i=1,nsph) 
 write(uboddump) ( oldderad(i),i=1,nsph) 
 write(uboddump) ( phi(i),i=1,nbodies) 
 write(uboddump) ( phiext(i),i=1,nbodies) 
 write(uboddump) ( pos(i,1:ndim),i=1,nbodies) 
 write(uboddump) ( rho(i),i=1,nsph) 
 write(uboddump) ( snentropy(i),i=1,nbodies) 
 write(uboddump) ( starfuv(i),i=1,nbodies) 
 write(uboddump) ( tcollaps(i),i=1,nsph) 
 write(uboddump) ( temperat(i),i=1,nsph) 
 write(uboddump) ( tfeedb(i),i=1,nbodies) 
 write(uboddump) ( tform(i),i=1,nbodies) 
 write(uboddump) ( tvel(i),i=1,nbodies) 
 write(uboddump) ( vdisp(i),i=1,nsph) 
 write(uboddump) ( vel(i,1:ndim),i=1,nbodies) 
 write(uboddump) ( veltpos(i,1:ndim),i=1,nsph) 
 close(uboddump) 
end subroutine

subroutine readdump(n) 
 include 'globals.h' 
 integer ioerror,n,i 
 real dummy 
 logical ldummy
 open(unit=uboddump,file=dumpfile,status='OLD',form='UNFORMATTED',iostat=ioerror) 
 if(ioerror.NE.0) then 
 print*,' dump file not present' 
 stop 
 endif 

read(uboddump) n 

      read(uboddump)  output,firstsnap,stepout,steplog,verbosity,       &
   datadir,inputfile,outputfile,halofile
      read(uboddump)  pboxsize,usesph,radiate,starform,cosmo
      read(uboddump)  unitl_in_kpc,unitm_in_msun
      read(uboddump)  dtime,tstepcrit,tstpcr2,freev,freea,freevexp,   &
   freeaexp,nsteps,max_tbin,minppbin,sqrttstp,acc_tstp, &
   freetstp
      read(uboddump)  bh_tol,eps,gdgtol,nn_tol,targetnn,usequad,         &
   directsum,selfgrav,fixthalo,adaptive_eps,gdgop
      read(uboddump)  epsgas,gamma,alpha,beta,epssph,courant,removgas, &
   consthsm,nsmtol,nsmooth,smoothinput,consph,sphinit,uentropy, &
   isotherm,eps_is_h,hupdatemethod,sph_visc
      read(uboddump)  graineff,crionrate,heat_par1,heat_par2,cool_par,    &
   optdepth
      read(uboddump)  comove
      read(uboddump)  tcollfac,masscrit,sfeff,tbubble,sne_eff,tsnbeg,   &
   rhomax,sfmode,feedback
      read(uboddump)  rndtable,rnseed
      read(uboddump)  nbodies,nsnap
      read(uboddump)  rmin,rsize,incells,incellsg
      read(uboddump)  root
      read(uboddump)  nttot,ntmin,ntmax,ntavg,nttotfuv,ntminfuv,        &
   ntmaxfuv,ntavgfuv
      read(uboddump)  nstot,nsmin,nsmax,nsavg
      read(uboddump)  tnow,tpos
      read(uboddump)  tiny
      read(uboddump)  mtot,etot,ektot,eptot,mstar,mgas,snheat,esofttot, &
   enerror,amvec,cmpos,cmvel
      read(uboddump)  input
      read(uboddump)  npactive,nsphact
      read(uboddump)  etol,dummy,tsteppos,active_bin
      read(uboddump)  syncflag,entropyflag
      read(uboddump)  teth,ethtot
      read(uboddump)  symmetry
      read(uboddump)  dummy,dummy,gamma1,nstar,nsph
      read(uboddump)  massres
      read(uboddump)  deldr2i
      read(uboddump)  wsmooth,dwsmooth
      read(uboddump)  nntot,nnmin,nnmax,nnavg
      read(uboddump)  hboxsize
      read(uboddump)  poshalo,masshalo
      read(uboddump)  eradiate,trad,meanmwt,     &
   mhboltz,fhydrogn,mumhkbol,mumhkgam, &
   mumhkg1,efuvheat,eradcool,mcold,mluke,mwarm,mhot,estar, &
   heatconst,densconst,timescale,lengthscale,velscale,flxscale, &
   massscale
      read(uboddump)  tstarform,tsnfeedback,snenergy,totptag
      read(uboddump)  h2time
      read(uboddump)  smoothuv
      read(uboddump)  tbh,nbh
      read(uboddump)  tpm,rcut,rcut2,pmsoft,pmdr,pmpot,       &
   pmacc
      read(uboddump)  searchpos,searchh,searchdelta,searchacc4,     &
   searchn,searchreuse,reuseflag,ncalls,nsearches
 read(uboddump) ( itimestp(i),i=1,nbodies) 
 read(uboddump) ( nbexist(i),i=1,nbodies) 
 read(uboddump) ( order_bodlist(i),i=1,nbodies) 
 read(uboddump) ( otimestp(i),i=1,nbodies) 
 read(uboddump) ( pactive(i),i=1,nbodies) 
 read(uboddump) ( acc(i,1:ndim+1),i=1,nbodies) 
 read(uboddump) ( csound(i),i=1,nsph) 
 read(uboddump) ( derad(i),i=1,nsph) 
 read(uboddump) ( dethdt(i),i=1,nsph) 
 read(uboddump) ( dethold(i),i=1,nsph) 
 read(uboddump) ( drhodh(i),i=1,nsph) 
 read(uboddump) ( elecfrac(i),i=1,nsph) 
 read(uboddump) ( epsgrav(i),i=1,nbodies) 
 read(uboddump) ( esnthdt(i),i=1,nsph) 
 read(uboddump) ( ethermal(i),i=1,nsph) 
 read(uboddump) ( ethold(i),i=1,nsph) 
 read(uboddump) ( fuvheat(i),i=1,nbodies) 
 read(uboddump) ( h2frac(i),i=1,nsph) 
 read(uboddump) ( hsmcurlv(i),i=1,nsph) 
 read(uboddump) ( hsmdivv(i),i=1,nsph) 
 read(uboddump) ( hsmooth(i),i=1,nbodies) 
 read(uboddump) ( mass(i),i=1,nbodies) 
 read(uboddump) ( mumaxdvh(i),i=1,nsph) 
 read(uboddump) ( oldderad(i),i=1,nsph) 
 read(uboddump) ( phi(i),i=1,nbodies) 
 read(uboddump) ( phiext(i),i=1,nbodies) 
 read(uboddump) ( pos(i,1:ndim),i=1,nbodies) 
 read(uboddump) ( rho(i),i=1,nsph) 
 read(uboddump) ( snentropy(i),i=1,nbodies) 
 read(uboddump) ( starfuv(i),i=1,nbodies) 
 read(uboddump) ( tcollaps(i),i=1,nsph) 
 read(uboddump) ( temperat(i),i=1,nsph) 
 read(uboddump) ( tfeedb(i),i=1,nbodies) 
 read(uboddump) ( tform(i),i=1,nbodies) 
 read(uboddump) ( tvel(i),i=1,nbodies) 
 read(uboddump) ( vdisp(i),i=1,nsph) 
 read(uboddump) ( vel(i,1:ndim),i=1,nbodies) 
 read(uboddump) ( veltpos(i,1:ndim),i=1,nsph) 
 close(uboddump) 
end subroutine
