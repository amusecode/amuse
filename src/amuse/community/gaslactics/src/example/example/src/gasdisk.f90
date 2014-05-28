function rnd()
 real rnd
 call random_number(rnd)
end function

function nrnd()
 real,save :: x,r=-1    
 real :: nrnd,r1,y
 if(r.ge.0) then
  nrnd=x
  r=-1
 endif
 r=1
 do while(r.ge.1)
  x=2*rnd()-1
  y=2*rnd()-1
  r=sqrt(x*x+y*y)
 enddo
  if(r.ne.0) then
   r1=2*sqrt(-log(r))
   x=r1/r*x
   nrnd=r1/r*y
  else
   nrnd=0.
  endif
end function

function sech2rnd()
 real sech2rnd
 sech2rnd=atanh(1-2*rnd())
end function 
 
program gas
  use hydrostaticdiskMod
 real, parameter :: pi=3.141592654
 real, allocatable :: rad(:),crad(:),csigtab(:)
 real nrnd,mass,drs,vrot,fipart,rpart,rrpart,
  rr(20001),vc(20001),dum1,dum2,cm,dr
 integer ngas,nrs,bin,bin1,imax,i,j,jm,nfis,nr
 real m,x(3),v(3)
 real*8 m8,x8(3),v8(3)

 print*,'ngas, nr?'
 read*, ngas,nrs
 print*,' velocity dispersion?'
 read*, vdisp
 print*,' csound?'
 read*, vdispt
 
 call readgasdisk
 nr=nrQ()
 if(mod(ngas,nrs).NE.0) print*, 'warning, ngas/nrs not int'
 nfis=ngas/nrs
 
 vdisp=0 
 vdispt=csQ()  
 print*,vdispt
 allocate(rad(0:nr),csigtab(0:nr),crad(0:nrs))
 rad=radQ() 
 csigtab=csigtabQ()  

 imax=0
 do i=1,nr
  if(csigtab(i).GT.csigtab(imax)) imax=i
 enddo
 mass=csigtab(imax)

 dr=rad(imax)/imax
 drs=rad(imax)/nrs
 print*,' total mass of gasdisk:',mass,nrs,drs
 crad(0)=0
 do i=1,nrs-1
  cm=i*mass/nrs
  do j=0,imax
   if(csigtab(j).lt.cm) jm=j
  enddo
  crad(i)=dr*jm+dr*(cm-csigtab(jm))/(csigtab(jm+1)-csigtab(jm))
 enddo
 crad(nrs)=rad(imax)

 open(1,file='newfreqdbh.dat',status='old',err=999)
 read(1,*);read(1,*)
 nreq=0
77 read(1,*,end=99) rr(nreq+1),vc(nreq+1),dum1,dum2
    nreq=nreq+1
   goto 77
99 close(1)

 open(unit=1,file='gas',status='unknown',form='unformatted')
 write(1) ngas,ngas,0

m=mass/ngas
print*,'pmass:',m  
 do i=1,ngas

  bin=(i-1)/nfis+1
  rpart=sqrt(rnd()*(crad(bin)**2-crad(bin-1)**2)+crad(bin-1)**2)
  fipart=2*pi*(mod(i-1,nfis)+rnd())/nfis
 
  if(rpart.le.rr(nreq)) then
   bin=1
   bin1=nreq 
   do while((bin1-bin).GT.1)
    ibin=(bin+bin1)/2
    if(rpart.GT.rr(ibin)) then
      bin=ibin  
    else
      bin1=ibin
    endif
   enddo
   db=(rpart-rr(bin))/(rr(bin1)-rr(bin))
   vrot=vc(bin)*(1-db)+db*vc(bin1)
  else
   vrot=vc(nreq)*sqrt(rr(nreq)/rpart) 
  endif
  x(1)=rpart*cos(fipart)
  x(2)=rpart*sin(fipart)
  x(3)=getz(rpart,rnd())
  if(rnd().GT.0.5) x(3)=-x(3)
  
! drift correction: checkcheck
  vrot=sqrt(max(0.,vrot**2-driftcorfac(rpart)*vdispt**2))
!  print*,rpart,vrot,driftcorfac(rpart)
  v(1)=-vrot*sin(fipart)
  v(2)=vrot*cos(fipart)
 
  v(1)=v(1)+vdisp*nrnd()
  v(2)=v(2)+vdisp*nrnd()
  v(3)=v(3)+vdisp*nrnd()
 
  m8=m
  x8=x
  v8=v
  write(1) m8,x8,v8
 enddo
 
 close(1)
 print*,' finished making gasdisk, output written to file: gas'

 stop
999 print*,' oops - error reading freqs'
end program
