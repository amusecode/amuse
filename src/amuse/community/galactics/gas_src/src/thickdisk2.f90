module thickdiskMOD
 private
 public :: initthickdisk,fi,sets,fint,dfidr,dfidz,dfint,ddfint, &
  firstcall,readdisk,writedisk

 real, parameter :: pi=3.141592654
 real, parameter :: dx=.05
 character(*), parameter :: datafile='gaspot.dat'

! fi approximately ok to rmax/dx/2 ..
! dr,dslice not super important 
 real :: rmax,zmax,kmax,dk,dr,dslice,router,totalmass
 integer :: nr,nk,nslices
 real, allocatable :: s(:,:)
 real, allocatable :: slicez(:)

 real, parameter :: fifac=1.2
 integer :: nfi 
 real :: dfi
 real, allocatable :: fitable(:,:),fir(:),dfidrtable(:,:),dfidztable(:,:), &
                       dfitable(:),ddfitable(:)

 logical, save :: firstcall=.true.

 contains

subroutine readdisk
 if (.not. firstcall) return
 firstcall=.false.
 open(unit=1,file=datafile,status='OLD',form='UNFORMATTED')
 read(1) rmax,zmax,kmax,dk,dr,dslice,router,totalmass,dfi
 read(1) nr,nk,nslices,nfi
 allocate(s(0:nslices,0:nk),slicez(0:nslices))
 allocate(fitable(0:nfi,0:nfi),dfidrtable(0:nfi,0:nfi),dfidztable(0:nfi,0:nfi), &
       dfitable(0:nfi),ddfitable(0:nfi),fir(0:nfi)) 
 read(1) s
 read(1) slicez
 read(1) fitable
 read(1) fir
 read(1) dfidrtable
 read(1) dfidztable
 read(1) dfitable
 read(1) ddfitable
 close(1)
 write(0,*) ' thickdisk read:', datafile
end subroutine

subroutine writedisk
 print*,' thickdisk writing:', datafile
 open(unit=1,file=datafile,status='UNKNOWN',form='UNFORMATTED')
 write(1) rmax,zmax,kmax,dk,dr,dslice,router,totalmass,dfi
 write(1) nr,nk,nslices,nfi
 write(1) s
 write(1) slicez
 write(1) fitable
 write(1) fir
 write(1) dfidrtable
 write(1) dfidztable
 write(1) dfitable
 write(1) ddfitable
 close(1)
end subroutine

subroutine initthickdisk(rm,zm,lnr,nz,tmass)
 real, external :: dens
 real, optional :: rm,zm,tmass
 integer, optional :: lnr,nz
 real :: mm
 integer :: iz

 firstcall=.false.
 print*,' init thickdiskmod'

 if(present(rm)) then
  rmax=rm
 else
  print*,'rmax?'
  read*,rmax
 endif  

 if(present(zm)) then
  zmax=zm
 else
  print*,'zmax?'
  read*,zmax
 endif  

 if(present(lnr)) then
  nr=lnr
 else
  print*,'# of radii?'
  read*,nr
 endif  

 if(present(nz)) then
  nslices=nz
 else
  print*,'# of slices?'
  read*,nslices
 endif  

 if(present(tmass)) then
  totalmass=tmass
 else
  print*,'target mass?'
  read*,totalmass
 endif  

 dr=rmax/nr
 dslice=zmax/nslices
 nk=nr/dx
 kmax=1/dr 
 dk=dx/rmax
 router=rmax/dx

 allocate(s(0:nslices,0:nk),slicez(0:nslices))

 s=0
 do iz=0,nslices
  slicez(iz)=dslice*iz
 enddo
 print*,'total mass:',totalmass

 nfi=nr
 allocate(fitable(0:nfi,0:nfi),dfitable(0:nfi),ddfitable(0:nfi),fir(0:nfi)) 
 allocate(dfidrtable(0:nfi,0:nfi),dfidztable(0:nfi,0:nfi)) 
 fir(0)=0
 do i=0,nfi
  fir(i)=router*fifac**(-nfi+i)
 enddo
 print*,'rmin:', fir(1)
 fitable=0
 dfidrtable=0
 dfidztable=0

end subroutine

 function fint(r,z) result(fi)
  real fi,r,z,du,dt,lz
  integer u,t  
  lz=abs(z)
  if(r.ge.router.or.lz.ge.router) then
   fi=totalmass/SQRT(r*r+lz*lz)
   return
  else
  u=log(r/router)/log(fifac)+nfi
  t=log(lz/router)/log(fifac)+nfi
  if(u.lt.0) u=0
  if(t.lt.0) t=0
  du= (r-fir(u))/(fir(u+1)-fir(u))
  dt= (lz-fir(t))/(fir(t+1)-fir(t))
   fi= du*dt*fitable(u+1,t+1)+(1-du)*dt*fitable(u,t+1)+ &
       du*(1-dt)*fitable(u+1,t)+(1-dt)*(1-du)*fitable(u,t)
  endif
 end function

 function dfidr(r,z) result(dfi)
  real dfi,r,z,du,dt,lz
  integer u,t  
  lz=abs(z)
  if(r.ge.router.or.lz.ge.router) then
   dfi=totalmass*r/(r*r+lz*lz)**(1.5)
   return
  else
  u=log(r/router)/log(fifac)+nfi
  t=log(lz/router)/log(fifac)+nfi
  if(u.lt.0) u=0
  if(t.lt.0) t=0
  du= (r-fir(u))/(fir(u+1)-fir(u))
  dt= (lz-fir(t))/(fir(t+1)-fir(t))
   dfi= du*dt*dfidrtable(u+1,t+1)+(1-du)*dt*dfidrtable(u,t+1)+ &
       du*(1-dt)*dfidrtable(u+1,t)+(1-dt)*(1-du)*dfidrtable(u,t)
  endif
 end function

 function dfidz(r,z) result(dfi)
  real dfi,r,z,du,dt,lz,s
  integer u,t  
  s=1; if(z.lt.0) s=-1
  lz=abs(z)
  if(r.ge.router.or.lz.ge.router) then
   dfi=-s*totalmass*lz/(r*r+lz*lz)**(1.5)
   return
  else
  u=log(r/router)/log(fifac)+nfi
  t=log(lz/router)/log(fifac)+nfi
  if(u.lt.0) u=0
  if(t.lt.0) t=0
  du= (r-fir(u))/(fir(u+1)-fir(u))
  dt= (lz-fir(t))/(fir(t+1)-fir(t))
   dfi= du*dt*dfidztable(u+1,t+1)+(1-du)*dt*dfidztable(u,t+1)+ &
       du*(1-dt)*dfidztable(u+1,t)+(1-dt)*(1-du)*dfidztable(u,t)
  endif
  dfi=dfi*s
 end function


 function dfint(r) result(fi)
  real fi,r,du
  integer u  
  if(r.ge.router) then
   fi=totalmass/(r*r)
   return
  else
   u=log(r/router)/log(fifac)+nfi
   if(u.lt.0) u=0
   du= (r-fir(u))/(fir(u+1)-fir(u))
   fi=(1-du)*dfitable(u)+du*dfitable(u+1)
  endif  
 end function

 function ddfint(r) result(fi)
  real fi,r,du
  integer u  
  if(r.ge.router) then
   fi=-2*totalmass/(r*r*r)
   return
  else
  u=log(r/router)/log(fifac)+nfi
  if(u.lt.0) u=0
  du= (r-fir(u))/(fir(u+1)-fir(u))
  fi=(1-du)*ddfitable(u)+du*ddfitable(u+1)
  endif  
 end function
 
 subroutine fillfi
  integer i,j
  
  do i=0,nfi
   do j=0,nfi
    fitable(i,j)=fi(fir(i),fir(j))
   enddo
  enddo

  do j=0,nfi
   dfidrtable(0,j)=0.
   do i=1,nfi
    dfidrtable(i,j)=(fi(1.05*fir(i),fir(j))-fi(0.95*fir(i),fir(j)))/(0.1*fir(i))
   enddo
  enddo
 ! dfidrtable(0,0)=(fi(0.01*fir(1),0.)-fi(0.,0.))/(0.01*fir(1))

  do i=0,nfi
   dfidztable(i,0)=0.
   do j=1,nfi
    dfidztable(i,j)=(fi(fir(i),1.05*fir(j))-fi(fir(i),0.95*fir(j)))/(0.1*fir(j))
   enddo
  enddo
 ! dfidztable(0,0)=(fi(0.,0.01*fir(1))-fi(0.,0.))/(0.01*fir(1))

 dfitable(0)=(fi(0.1*fir(1),0.)-fi(0.,0.))/(0.1*fir(1))
 do i=1,nfi
  dfitable(i)=(fi(1.05*fir(i),0.)-fi(0.95*fir(i),0.))/(0.1*fir(i))
 enddo

 ddfitable(0)=(dfitable(1)-dfitable(0))/(fir(1))
 ddfitable(nfi)=(dfitable(nfi)-dfitable(nfi-1))/(fir(nfi)-fir(nfi-1))
 do i=1,nfi-1
  ddfitable(i)=(dfitable(i+1)-dfitable(i-1))/(fir(nfi+1)-fir(nfi-1))
 enddo

 end subroutine

subroutine sets(dens)
  integer :: i,j,iz
  real :: k,r,mm
  real, external :: dens
 do iz=0,nslices
 do i=0,nk
  k=i*dk
  s(iz,i)=0.5*bessj0(k*rmax)*dslice*dens(rmax,slicez(iz))*rmax
  do j=1,nr-1
   r=j*dr
   s(iz,i)=s(iz,i)+bessj0(r*k)*dslice*dens(r,slicez(iz))*r
  enddo
  s(iz,i)=-2*pi*s(iz,i)*dr
 enddo
 enddo
 mm=-(s(0,0)+s(nslices,0)+2*sum(s(1:nslices-1,0)))
! print*,'mass:', mm
 if(mm.EQ.0) call error('zero mass??')
 s=s/mm*totalmass
 call fillfi
 print*,'filled fi'
end subroutine

function ftilde(k,z) result(f)
  real k,z,f,az,dz
 az=abs(z)
 dz=dslice
 if(k.EQ.0.) then
  f=1.
  return
 endif  
 if(az.GT.dz/2.) then
  f=exp(-k*az)/(dz*k)*(exp(k*dz/2.)-exp(-k*dz/2.))
 else
  f=1/(dz*k)*(2.-exp(-k*dz/2.)*(exp(-k*az)+exp(k*az)))
 endif
end function

function fi(r,z) 
  real :: r,z,fi,k,zs
  integer :: i,j,iz
 fi=0
 do i=-nslices,nslices
  zs=i*dslice
  iz=i
  if(iz.LT.0) iz=-iz
  fi=fi+0.5*(s(iz,0)*ftilde(0.,z-zs)*bessj0(0)+ &
            s(iz,nk)*ftilde(kmax,z-zs)*bessj0(r*kmax))  
  do j=1,nk-1
   k=j*dk
   fi=fi+s(iz,j)*ftilde(k,z-zs)*bessj0(r*k)
  enddo
 enddo
 fi=fi*dk
end function

 subroutine error(string,i)
  character(*) :: string
  integer, optional :: i
  
  print*,'error detected in thickdiskmod'

  if(present(i)) then
   print*,string,i
  else
   print*,string
  endif
  stop
end subroutine

end module

module hydrostaticdiskMod

 private
 public :: initdisk,dens,zEquilibrium,getz,totalmassQ,driftcorfac,radQ, &
                 setsimplez,readgasdisk,writegasdisk,csigtabQ,firstcall, &
		 nrQ,csQ,getrad
		 


 real, parameter :: pi=3.141592654
 character(*), parameter :: datafile='gasdens.dat'
 integer :: nr,nz
 real :: rmax,zmax,dr,dz
 
 real, allocatable :: rad(:),zet(:)
 real, allocatable :: sigtab(:)
 real, allocatable :: csigtab(:)
 real, allocatable :: denstab(:,:)
 real, allocatable :: cdenstab(:,:)
 real, allocatable :: ddenstab(:)
 real :: cs,gamma,totalmass

 logical, save :: firstcall=.true.
 
 contains

function csQ()
 real :: csQ
 csQ=cs
end function



function radQ()
 real :: radQ(0:nr)
 radQ=rad
end function

function nrQ()
 integer :: nrQ
 nrQ=nr
end function

function csigtabQ()
 real :: csigtabQ(0:nr)
 csigtabQ=csigtab
end function

subroutine initDisk(sigma,rextent,zextent,nrad,nzet,csnd,gam,tmass)
 real, external :: sigma
 real, optional :: rextent,zextent,tmass,csnd,gam
 integer, optional :: nrad,nzet
 print*,'init hydrostatic disk'

 firstcall=.false.
 if(present(rextent)) then
  rmax=rextent
 else
  print*,'rmax?'
  read*,rmax
 endif  

 if(present(zextent)) then
  zmax=zextent
 else
  print*,'zmax?'
  read*,zmax
 endif  

 if(present(nrad)) then
  nr=nrad
 else
  print*,'nrad?'
  read*,nr
 endif  

 if(present(nzet)) then
  nz=nzet
 else
  print*,'nz?'
  read*,nz
 endif  

 if(present(csnd)) then
  cs=csnd
 else
  print*,'csound?'
  read*,cs
 endif  

 if(present(gam)) then
  gamma=gam
 else
  print*,'gamma?'
  read*,gamma
 endif  

 dr=rmax/nr;dz=zmax/nz
 allocate(rad(0:nr), sigtab(0:nr),csigtab(0:nr), &
             denstab(0:nr,0:nz),cdenstab(0:nr,0:nz),zet(0:nz)) 
 totalmass=0

 print*,dr,dz,nr,nz
 do i=0,nr
  r=i*dr
  rad(i)=r
  sigtab(i)=sigma(r)
 enddo
 csigtab(0)=0.
 do i=1,nr
  csigtab(i)=csigtab(i-1)+2*Pi*dr*(rad(i)*sigtab(i)+rad(i-1)*sigtab(i-1))/2.
 enddo
 totalmass=csigtab(nr)
 print*,'hydrostaticdisk mass:',totalmass
 if(present(tmass)) then
  print*,' mass corrected to:',tmass
  sigtab=sigtab*tmass/totalmass  
  csigtab=csigtab*tmass/totalmass  
 endif 
 do j=0,nz
  zet(j)=j*dz
 enddo
 denstab=0
 cdenstab=0
 allocate(ddenstab(0:nr))
 call setdd
end subroutine

function totalmassQ()
 real totalmassQ
 totalmassQ=totalmass
end function 

function dens(r,z)
  real dens,r,z,du,dt,lz
  integer u,t  
  lz=abs(z)
  if(r.gt.rmax.or.lz.gt.zmax) then
   dens=0
   return
  else
   u=r/dr
   t=lz/dz
   du=r/dr-u
   dt=lz/dz-t
   dens= du*dt*denstab(u+1,t+1)+(1-du)*dt*denstab(u,t+1)+ &
       du*(1-dt)*denstab(u+1,t)+(1-dt)*(1-du)*denstab(u,t)
  endif  
end function

subroutine zEquilibrium(pot)
 real,external :: pot
 real :: r,z,fi,fi0,signorm
 integer :: i,j
 do i=0,nr
  r=i*dr
  fi0=pot(r,0.)
  do j=0,nz
   z=j*dz
   fi=pot(r,z)
   denstab(i,j)=exp(-gamma*(fi-fi0)/cs**2)   
  enddo   
  signorm=(denstab(i,0)+denstab(i,nz)+2*sum(denstab(i,1:nz-1)))*dz
  denstab(i,0:nz)=denstab(i,0:nz)*sigtab(i)/signorm
  cdenstab(i,0)=0.
  do j=1,nz
   cdenstab(i,j)=cdenstab(i,j-1)+dz*(denstab(i,j)+denstab(i,j-1))/2.
  enddo
 enddo
 call setdd
end subroutine

subroutine setsimplez(zdisk)
 real :: r,z,fi,fi0,signorm,zdisk,con
 integer :: i,j
 do i=0,nr
  r=i*dr
  fi0=pot(r,0.)
  do j=0,nz
   z=j*dz
   con = exp(-abs(z/zdisk))
   con = (2.0*con/(1.0 + con*con))**2
   denstab(i,j)=con   
  enddo   
  signorm=(denstab(i,0)+denstab(i,nz)+2*sum(denstab(i,1:nz-1)))*dz
  denstab(i,0:nz)=denstab(i,0:nz)*sigtab(i)/signorm
  cdenstab(i,0)=0.
  do j=1,nz
   cdenstab(i,j)=cdenstab(i,j-1)+dz*(denstab(i,j)+denstab(i,j-1))/2.
  enddo
 enddo
 call setdd
end subroutine

subroutine setdens(n,m,indens)
 integer :: n,m,i,j
 real :: indens(0:n,0:m),r,signorm
 if(nr.NE.n.OR.nz.NE.m) call error('array mismatch')
 denstab=indens
 do i=0,nr
  r=i*dr
  signorm=(denstab(i,0)+denstab(i,nz)+2*sum(denstab(i,1:nz-1)))*dz
  denstab(i,0:nz)=denstab(i,0:nz)*sigtab(i)/signorm
  cdenstab(i,0)=0.
  do j=1,nz
   cdenstab(i,j)=cdenstab(i,j-1)+dz*(denstab(i,j)+denstab(i,j-1))/2.
  enddo
 enddo
 call setdd
end subroutine

subroutine readgasdisk
 if (.not. firstcall) return
 firstcall=.false.
 open(unit=1,file=datafile,status='old',form='UNFORMATTED')
 read(1) rmax,zmax,nr,nz,cs,gamma,totalmass,dr,dz 
 allocate(rad(0:nr), sigtab(0:nr),csigtab(0:nr), &
             denstab(0:nr,0:nz),cdenstab(0:nr,0:nz),zet(0:nz))  
 read(1) rad,zet
 read(1) sigtab 
 read(1) csigtab
 read(1) denstab
 read(1) cdenstab
 close(1)
 allocate(ddenstab(0:nr))
 call setdd
 write(0,*) ' hydrostaticdisk read:', datafile
end subroutine

subroutine writegasdisk
 print*,' hydrostaticdisk writing:', datafile
 open(unit=1,file=datafile,status='unknown',form='UNFORMATTED')
 write(1) rmax,zmax,nr,nz,cs,gamma,totalmass,dr,dz 
 write(1) rad,zet
 write(1) sigtab 
 write(1) csigtab
 write(1) denstab
 write(1) cdenstab
 close(1)
end subroutine

function driftcorfac(r) result(x)
 real :: x,r,du,d,dd
 integer :: u

  d=dens(r,0.)
  if(r.gt.rmax.OR.d.EQ.0) then
   x=0
   return
  else
   u=r/dr
   du=r/dr-u
   dd=du*ddenstab(u+1)+(1-du)*ddenstab(u)
   x=-r*dd/d
  endif  
end function

subroutine setdd
 ddenstab(0)=(denstab(1,0)-denstab(0,0))/(rad(1)-rad(0))
 ddenstab(nr)=(denstab(nr,0)-denstab(nr-1,0))/(rad(nr)-rad(nr-1))
 do i=1,nr-1
 ddenstab(i)=(denstab(i+1,0)-denstab(i-1,0))/(rad(i+1)-rad(i-1))
 enddo
end subroutine

function getz(r,rnd)
 real getz,rnd,z1,z2,du,y1,y2
 integer :: u
 getz=0
 u=r/dr
 du=r/dr-u
 if(u.GE.nr.OR.u.LT.0) call error('r out of range',u)
 y1=cdenstab(u,nz)*rnd
 y2=cdenstab(u+1,nz)*rnd
 z1=invertcumul(nz,zet(0:nz),cdenstab(u,0:nz),y1)
 z2=invertcumul(nz,zet(0:nz),cdenstab(u+1,0:nz),y2)
 getz=du*z2+(1-du)*z1
end function 

function getrad(rnd)
 real getrad,rnd,y
 integer :: u
 getr=0
 y=csigtab(nr)*rnd
 getrad=invertcumul(nr,rad(0:nr),csigtab(0:nr),y)
end function 

function invertcumul(n,xlist,ylist,y) result(x)
 integer :: n,bin,up,low
 real :: xlist(n),ylist(n),x,y,u,s  
 if(ylist(1).LT.ylist(n)) then
  s=1
 else
  s=-1
 endif   
 if(s*y.LE.s*ylist(1)) then
  x=xlist(1)
  return
 endif
 if(s*y.GE.s*ylist(n)) then
  x=xlist(n)
  return
 endif
 up=n;low=1
 do while((up-low).GT.1)
   bin=(low+up)/2
   if(s*y.LT.s*ylist(bin)) then
    up=bin
   else
    low=bin
   endif   
 enddo 
 bin=up
 u=(y-ylist(bin-1))/(ylist(bin)-ylist(bin-1))
 x=(1-u)*xlist(bin-1)+u*xlist(bin)  
 return
end function

 subroutine error(string,i)
  character(*) :: string
  integer, optional :: i
  
  print*,'error detected in hydrostat'

  if(present(i)) then
   print*,string,i
  else
   print*,string
  endif
  stop
 end subroutine

end module

function sigma(r) result(x)
 real gasm,gasr,gastr,gasdt,csnd,gamma,gassig
 integer gasnr
 common /gasstuf/ gasm,gasr,gastr,gasdt,csnd,gamma,gasnr,gassig,gasz
 real r,x,tr,dtr
  dtr=(r-gastr)/gasdt
  tr=1
  if(dtr.gt.-20) then
   if(dtr.lt.10) then
    tr=1/(1.+exp(dtr))   
   else
    tr=0
  endif
  endif
   x=tr*gassig/(1+r/gasr)  
 end function

function thickdiskdens(r,z) result(ddens)
 use hydrostaticdiskMod
 real ddens,r,z
 if(firstcall) call readgasdisk
 ddens=dens(r,z)
end function

function thickdiskpot(r,z) result(pot)
  use thickdiskMOD
 real pot,r,z
 if(firstcall) call readdisk
 pot=fint(r,z)
end function 

function thickdiskforce(r) result(x) ! -fi'
  use thickdiskMOD
 real x,r,z
 if(firstcall) call readdisk
 x=-dfint(r)
end function 

function thickdiskrforce(r,z) result(x) ! -dfi/dr
  use thickdiskMOD
 real x,r,z
 if(firstcall) call readdisk
 x=-dfidr(r,z)
end function 

function thickdiskzforce(r,z) result(x) ! -dfi/dz
  use thickdiskMOD
 real x,r,z
 if(firstcall) call readdisk
 x=-dfidz(r,z)
end function 


function thickdiskf2(r) result(x) ! fi''
  use thickdiskMOD
 real x,r,z
 if(firstcall) call readdisk
 x=ddfint(r)
end function 

subroutine setgasdisk()
 use hydrostaticdiskMod
 use thickdiskMod
 
 real gasm,gasr,gastr,gasdt,csnd,gamma,gassig
 integer gasnr
 common /gasstuf/ gasm,gasr,gastr,gasdt,csnd,gamma,gasnr,gassig,gasz
 real rextent,zextent
 real,external :: sigma
 integer nr,nz,nslices
 parameter(pi=3.1415926535)

 gassig=gasm/(2*pi*gasr**2*(gastr/gasr-log(1+gastr/gasr)))
 rextent=gastr+10*gasdt
 zextent=gasz
 nr=gasnr
 nz=nr/2
 nslices=nr/5
 print*,'rextent,zextent:',rextent,zextent
 call initDisk(sigma,rextent,zextent,nr,nz,csnd,gamma,gasm)
 print*,'gas sigma0:',sigma(0.)
 call initthickdisk(rextent,zextent,nr,nslices,gasm)
 
end subroutine 

subroutine firstgasdens(zdisk)
 use hydrostaticdiskMod
 real zdisk
 call setsimplez(zdisk)
end subroutine 

subroutine solveZ
 use hydrostaticdiskMod
 real, external :: pot
 call zEquilibrium(pot)
end subroutine

subroutine solvepot
  use thickdiskMOD
 real, external :: thickdiskdens
 call sets(thickdiskdens)
end subroutine

subroutine writegas
 use hydrostaticdiskMod
 use ThickdiskMod
 call writegasdisk
 call writedisk
 print*,' gas dens+pot files written'
end subroutine 

