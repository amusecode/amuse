subroutine clean
  include 'globals.h'
  integer :: p,j,ndel
  real :: tm,cm(3),cmv(3)
  logical del
  ndel=0
  do p=1,nbodies
    del=.false.
    do j=1,ndim
      if(pos(p,j).GT.pboxsize/2.OR.pos(p,j).LT.-pboxsize/2) then
        del=.true.
      endif
    enddo
    if(del)then 
      ndel=ndel+1
      templist(ndel)=p
    endif
  enddo

  if(periodic.AND.ndel.GT.0) then
    print*,ndel
    call terror(' particles outside box?') 
  endif

  tm=0;cm=0;cmv=0 
  do j=1,ndel
    p=templist(j)
    tm=tm+mass(p)
    cm(1:ndim)=cm(1:ndim)+mass(p)*pos(p,1:ndim)
    cmv(1:ndim)=cmv(1:ndim)+mass(p)*vel(p,1:ndim)
    mass(p)=0.
  enddo

  if(ndel.GT.0) then
    if(verbosity.GT.0) print*,' mass erased:',tm
    open(unit=erasenr,file=TRIM(outputfile)//'.del',status='unknown',position='APPEND')
    write(erasenr,*) tnow,tm,cm/tm,cmv/tm
    close(erasenr)
  endif 

end subroutine

subroutine copysph(p,i) ! copy sph part p to spot i
  include 'globals.h'
  integer i,p

  if(i.GT.nsphmax.OR.p.GT.nsphmax) then
    print*, p,i
    call terror(' sph part copy error')
  endif
  call copypart(p,i)
  rho(i)=rho(p)
  drhodh(i)=drhodh(p)
  csound(i)=csound(p)
  derad(i)=derad(p)
  hsmdivv(i)=hsmdivv(p)
  mumaxdvh(i)=mumaxdvh(p)
  hsmcurlv(i)=hsmcurlv(p)
  veltpos(i,1)=veltpos(p,1)
  veltpos(i,2)=veltpos(p,2)
  veltpos(i,3)=veltpos(p,3)
  fuvheat(i)=fuvheat(p)
  tcollaps(i)=tcollaps(p)
  esnthdt(i)=esnthdt(p)
  temperat(i)=temperat(p)
  elecfrac(i)=elecfrac(p)
  ethold(i)=ethold(p)
  dethold(i)=dethold(p)
  h2frac(i)=h2frac(p)
  vdisp(i)=vdisp(p)
  entropy(i)=entropy(p) ! equiv ethermal
  dentdt(i)=dentdt(p)   ! equiv dethdt
end subroutine

subroutine copypart(p,i) ! copy part p to spot i
  include 'globals.h'
  integer i,p
  if(i.GT.nbodsmax.OR.p.GT.nbodsmax) then
    print*, p,i
    call terror(' sph part copy error')
  endif
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
end subroutine

subroutine partremoval
  include 'globals.h'
  integer i,lastsph,lastdm,laststar,lastbh 
  integer onsph,onstar,onbh,onbodies
  i=1

  pordercount=pordercount+1
  
  lastsph=nsph
  lastdm=nbodies-nstar
  laststar=nbodies-nbh
  lastbh=nbodies
  do while(i.LE.lastsph) 
    if(mass(i).LT.tiny) then
! AVE, AMUSE
! keep a record of the removed ids
      removedidssph(nremovals+1) = nbexist(i)
      nremovals=nremovals+1
! ---
      if(mass(lastsph).GE.tiny) then
        call copysph(lastsph,i)
        mass(lastsph)=0
      endif
      lastsph=lastsph-1
    else
      i=i+1
    endif 
  enddo

  do while(i.LE.lastdm) 
    if(mass(i).LT.tiny) then
      if(mass(lastdm).GE.tiny) then
        call copypart(lastdm,i)
        mass(lastdm)=0
      endif
      lastdm=lastdm-1
    else
      i=i+1
    endif 
  enddo

  do while(i.LE.laststar) 
    if(mass(i).LT.tiny) then
      if(mass(laststar).GE.tiny) then
        call copypart(laststar,i)
        mass(laststar)=0
      endif
      laststar=laststar-1
    else
      i=i+1
    endif 
  enddo

  do while(i.LE.lastbh) 
    if(mass(i).LT.tiny) then
      if(mass(lastbh).GE.tiny) then
        call copypart(lastbh,i)
        mass(lastbh)=0
      endif
        lastbh=lastbh-1
    else
      i=i+1
    endif 
  enddo
   
  nsph=lastsph
  nstar=lastbh-lastdm ! bh is also a star
  nbh=lastbh-laststar
  nbodies=lastbh 
end subroutine 

subroutine zeropart(p,n)
  include 'globals.h'
  integer p,n

  mass(p:p+n-1)=0
  pos(p:p+n-1,1)=0
  pos(p:p+n-1,2)=0
  pos(p:p+n-1,3)=0
  vel(p:p+n-1,1)=0
  vel(p:p+n-1,2)=0
  vel(p:p+n-1,3)=0
  acc(p:p+n-1,1)=0
  acc(p:p+n-1,2)=0
  acc(p:p+n-1,3)=0
  acc(p:p+n-1,4)=0
  phi(p:p+n-1)=0
  phiext(p:p+n-1)=0
  epsgrav(p:p+n-1)=0
  itimestp(p:p+n-1)=0
  otimestp(p:p+n-1)=0
  tform(p:p+n-1)=0
  starfuv(p:p+n-1)=0
  tfeedb(p:p+n-1)=0
  tvel(p:p+n-1)=0
  snentropy(p:p+n-1)=0
  hsmooth(p:p+n-1)=0
  nbexist(p:p+n-1)=0
end subroutine

subroutine zerosph(p,n)
  include 'globals.h'
  integer p,n

  call zeropart(p,n)

  rho(p:p+n-1)=0
  drhodh(p:p+n-1)=0
  csound(p:p+n-1)=0
  derad(p:p+n-1)=0
  hsmdivv(p:p+n-1)=0
  mumaxdvh(p:p+n-1)=0
  hsmcurlv(p:p+n-1)=0
  veltpos(p:p+n-1,1)=0
  veltpos(p:p+n-1,2)=0
  veltpos(p:p+n-1,3)=0
  fuvheat(p:p+n-1)=0
  tcollaps(p:p+n-1)=0
  esnthdt(p:p+n-1)=0
  temperat(p:p+n-1)=0
  elecfrac(p:p+n-1)=0
  ethold(p:p+n-1)=0
  dethold(p:p+n-1)=0
  h2frac(p:p+n-1)=0
  vdisp(p:p+n-1)=0
  entropy(p:p+n-1)=0
  dentdt(p:p+n-1)=0
end subroutine

function addspace_bh(n) result(p)
  include 'globals.h'
  integer p,n

  if(nbodies+n.GT.nbodsmax) call terror('add_bh: exceed nbodsmax')

  pordercount=pordercount+1

  p=nbodies+1
  call zeropart(p,n)
  nbodies=nbodies+n
  nbh=nbh+n
  nstar=nstar+n
 
end function

function addspace_star(n) result(p)
  include 'globals.h'
  integer i,p,n

  if(nbodies+n.GT.nbodsmax) call terror('add_star: exceed nbodsmax')

  pordercount=pordercount+1

  do i=1,min(n,nbh)
    call copypart(nbodies-nbh+i,nbodies+n-min(n,nbh)+i)
  enddo

  p=nbodies-nbh+1
  call zeropart(p,n)
  nbodies=nbodies+n
  nstar=nstar+n
end function

function addspace_dm(n) result(p)
  include 'globals.h'
  integer i,p,n

  if(nbodies+n.GT.nbodsmax) call terror('add_dm: exceed nbodsmax')

  pordercount=pordercount+1

  do i=1,min(n,nbh)
    call copypart(nbodies-nbh+i,nbodies+n-min(n,nbh)+i)
  enddo

  do i=1,min(n,nstar-nbh)
    call copypart(nbodies-nstar+i,nbodies-nbh+n-min(n,nstar-nbh)+i)
  enddo

  p=nbodies-nstar+1
  call zeropart(p,n)
  nbodies=nbodies+n 
end function

function addspace_gas(n) result(p)
  include 'globals.h'
  integer i,p,n

  if(nbodies+n.GT.nbodsmax.OR.nsph+n.GT.nsphmax) &
    call terror('add_gas: exceed nbodsmax')

  pordercount=pordercount+1

  do i=1,min(n,nbh)
    call copypart(nbodies-nbh+i,nbodies+n-min(n,nbh)+i)
  enddo

  do i=1,min(n,nstar-nbh)
    call copypart(nbodies-nstar+i,nbodies-nbh+n-min(n,nstar-nbh)+i)
  enddo

  do i=1,min(n,nbodies-nsph-nstar)
    call copypart(nsph+i,nbodies-nstar+n-min(n,nbodies-nsph-nstar)+i)
  enddo

  p=nsph+1
  call zerosph(p,n)
  nbodies=nbodies+n
  nsph=nsph+n 
end function

