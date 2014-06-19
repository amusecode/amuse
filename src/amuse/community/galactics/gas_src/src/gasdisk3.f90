function rnd()
 integer dummmy
 real*4 rnd,ran1
 external ran1
 rnd=ran1(dummy)
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
 
subroutine generate_gas(ngas,vdisp,iseed,m,x,y,z,vx,vy,vz,cs)
  use hydrostaticdiskMod
 real, parameter :: pi=3.141592654
 real, allocatable :: rad(:),crad(:),csigtab(:)
 real nrnd,mass,drs,vrot,fipart,rpart,rrpart,zpart,cm,dr
 integer ngas,nrs,bin,bin1,imax,i,j,jm,nfis,nr,ibuf(100)
 real m(ngas),x(ngas),y(ngas),z(ngas),vx(ngas),vy(ngas),vz(ngas),cs(ngas)
 real fr,fz,phi
 character*60 filename
 integer*8 iseed,itmp
 integer*8,parameter :: skip=100

 nrs=1000
 
 filename='dbh.dat'
 call readharmfile(filename,ibuf)
 
 call readgasdisk
 nr=nrQ()
 
 vdisp=0 
 vdispt=csQ()  
 allocate(rad(0:nr),csigtab(0:nr),crad(0:nrs))
 rad=radQ() 
 csigtab=csigtabQ()  

 imax=0
 do i=1,nr
  if(csigtab(i).GT.csigtab(imax)) imax=i
!  print*,rad(i),1-csigtab(i)/csigtab(nr)
 enddo
 mass=csigtab(imax)

 dr=rad(imax)/imax
 drs=rad(imax)/nrs
 crad(0)=0
 do i=1,nrs-1
  cm=i*mass/nrs
  do j=0,imax
   if(csigtab(j).lt.cm) jm=j
  enddo
  crad(i)=dr*jm+dr*(cm-csigtab(jm))/(csigtab(jm+1)-csigtab(jm))
 enddo
 crad(nrs)=rad(imax)

 write(0,*) 'cdisp:',vdispt
 write(0,*) 'total mass of gasdisk:',mass,nrs,drs
 write(0,*) 'rmax, mass:',rad(imax),mass/ngas

m=mass/ngas
cs=vdispt
 do i=1,ngas
   itmp=(iseed+i)*skip
   call ran_seed(itmp)
!  bin=(i-1)/nfis+1
!  rpart=sqrt(rnd()*(crad(bin)**2-crad(bin-1)**2)+crad(bin-1)**2)
  rpart=getrad(rnd())
!  fipart=2*pi*(mod(i-1,nfis)+rnd())/nfis
  fipart=2*pi*rnd()
  zpart=getz(rpart,rnd())

  x(i)=rpart*cos(fipart)
  y(i)=rpart*sin(fipart)
  z(i)=zpart
  if(rnd().GT.0.5) z(i)=-z(i)
  
  call force(rpart,zpart,fr,fz,phi)
  vrot=sqrt(max(0.,-rpart*fr))
!  vrot=sqrt(max(0.,rpart*fr))
! drift correction: checkcheck
!  print*,'zz',rpart,zpart,fr,phi
  vrot=sqrt(max(0.,vrot**2-driftcorfac(rpart)*vdispt**2))
  vx(i)=-vrot*sin(fipart)
  vy(i)=vrot*cos(fipart)
  vz(i)=0
 
  vx(i)=vx(i)+vdisp*nrnd()
  vy(i)=vy(i)+vdisp*nrnd()
  vz(i)=vz(i)+vdisp*nrnd()

 enddo
end subroutine
