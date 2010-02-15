subroutine mortonsort
 include 'globals.h'
 integer i,jg,js,jh,jb,k
  
 call maketree

 ttree=tnow-1. ! because tree is not valid after sorting
    
 jg=0
 js=0
 jh=0
 jb=0
 do i=1,nbodies
  k=order_bodlist(i)
  if(k.LE.nsph) then
   jg=jg+1
   bodlist(jg)=k
   templist(jg)=i
  endif
  if(k.GT.nsph.AND.k.LE.nbodies-nstar) then
   jh=jh+1
   bodlist(jh+nsph)=k
   templist(jh+nsph)=i
  endif
  if(k.GT.nbodies-nstar.AND.k.LE.nbodies-nbh) then
   js=js+1
   bodlist(js+nbodies-nstar)=k
   templist(js+nbodies-nstar)=i
  endif
  if(k.GT.nbodies-nbh) then
   jb=jb+1
   bodlist(jb+nbodies-nbh)=k
   templist(jb+nbodies-nbh)=i
  endif

 enddo
! bh=star 
 if(jg.NE.nsph.OR. &
    js.NE.nstar-nbh.OR. & 
    jh.NE.nbodies-nstar-nsph.OR. &
    jb.NE.nbh) call terror('ordering error') 

 do i=1,nbodies
  order_bodlist(templist(i))=i
 enddo

 tempvect(1:nbodies)=mass(bodlist(1:nbodies))
 mass(1:nbodies)=tempvect(1:nbodies)

 do i=1,ndim
  tempvect(1:nbodies)=pos(bodlist(1:nbodies),i)
  pos(1:nbodies,i)=tempvect(1:nbodies)
 enddo

 do i=1,ndim
  tempvect(1:nbodies)=vel(bodlist(1:nbodies),i)
  vel(1:nbodies,i)=tempvect(1:nbodies)
 enddo

 do i=1,ndim+1
  tempvect(1:nbodies)=acc(bodlist(1:nbodies),i)
  acc(1:nbodies,i)=tempvect(1:nbodies)
 enddo

 tempvect(1:nbodies)=phi(bodlist(1:nbodies))
 phi(1:nbodies)=tempvect(1:nbodies)

 tempvect(1:nbodies)=phiext(bodlist(1:nbodies))
 phiext(1:nbodies)=tempvect(1:nbodies)

 tempvect(1:nbodies)=epsgrav(bodlist(1:nbodies))
 epsgrav(1:nbodies)=tempvect(1:nbodies)

 templist(1:nbodies)=itimestp(bodlist(1:nbodies))
 itimestp(1:nbodies)=templist(1:nbodies)

 templist(1:nbodies)=otimestp(bodlist(1:nbodies))
 otimestp(1:nbodies)=templist(1:nbodies)

 tempvect(1:nbodies)=tform(bodlist(1:nbodies))
 tform(1:nbodies)=tempvect(1:nbodies)

 tempvect(1:nbodies)=starfuv(bodlist(1:nbodies))
 starfuv(1:nbodies)=tempvect(1:nbodies)

 tempvect(1:nbodies)=tfeedb(bodlist(1:nbodies))
 tfeedb(1:nbodies)=tempvect(1:nbodies)

 tempvect(1:nbodies)=tvel(bodlist(1:nbodies))
 tvel(1:nbodies)=tempvect(1:nbodies)

 tempvect(1:nbodies)=snentropy(bodlist(1:nbodies))
 snentropy(1:nbodies)=tempvect(1:nbodies)

 tempvect(1:nbodies)=hsmooth(bodlist(1:nbodies))
 hsmooth(1:nbodies)=tempvect(1:nbodies)

 templist(1:nbodies)=nbexist(bodlist(1:nbodies))
 nbexist(1:nbodies)=templist(1:nbodies)

 tempvect(1:nsph)=rho(bodlist(1:nsph))
 rho(1:nsph)=tempvect(1:nsph)

 tempvect(1:nsph)=drhodh(bodlist(1:nsph))
 drhodh(1:nsph)=tempvect(1:nsph)

 tempvect(1:nsph)=csound(bodlist(1:nsph))
 csound(1:nsph)=tempvect(1:nsph)

 tempvect(1:nsph)=derad(bodlist(1:nsph))
 derad(1:nsph)=tempvect(1:nsph)

 tempvect(1:nsph)=hsmdivv(bodlist(1:nsph))
 hsmdivv(1:nsph)=tempvect(1:nsph)

 tempvect(1:nsph)=mumaxdvh(bodlist(1:nsph))
 mumaxdvh(1:nsph)=tempvect(1:nsph)

 tempvect(1:nsph)=hsmcurlv(bodlist(1:nsph))
 hsmcurlv(1:nsph)=tempvect(1:nsph)

 do i=1,ndim
  tempvect(1:nsph)=veltpos(bodlist(1:nsph),i)
  veltpos(1:nsph,i)=tempvect(1:nsph)
 enddo

 tempvect(1:nsph)=fuvheat(bodlist(1:nsph))
 fuvheat(1:nsph)=tempvect(1:nsph)

 tempvect(1:nsph)=tcollaps(bodlist(1:nsph))
 tcollaps(1:nsph)=tempvect(1:nsph)

 tempvect(1:nsph)=esnthdt(bodlist(1:nsph))
 esnthdt(1:nsph)=tempvect(1:nsph)

 tempvect(1:nsph)=temperat(bodlist(1:nsph))
 temperat(1:nsph)=tempvect(1:nsph)

 tempvect(1:nsph)=elecfrac(bodlist(1:nsph))
 elecfrac(1:nsph)=tempvect(1:nsph)

 tempvect(1:nsph)=ethold(bodlist(1:nsph))
 ethold(1:nsph)=tempvect(1:nsph)

 tempvect(1:nsph)=dethold(bodlist(1:nsph))
 dethold(1:nsph)=tempvect(1:nsph)

 tempvect(1:nsph)=h2frac(bodlist(1:nsph))
 h2frac(1:nsph)=tempvect(1:nsph)

 tempvect(1:nsph)=vdisp(bodlist(1:nsph))
 vdisp(1:nsph)=tempvect(1:nsph)

 tempvect(1:nsph)=ethermal(bodlist(1:nsph))
 ethermal(1:nsph)=tempvect(1:nsph)

 tempvect(1:nsph)=dethdt(bodlist(1:nsph))
 dethdt(1:nsph)=tempvect(1:nsph)
 
end subroutine
