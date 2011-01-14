module fixedhalomod
 implicit none

private
public :: initfixedhalo,endfixedhalo,fixedhaloacc,fixedhalopot,fixedhalotidalfield
public :: fixedhalopos, fixedhalomass

  integer, parameter :: iounit=19 
  real,parameter :: pi=3.14159265358979

  real :: fixedhalopos(3),fixedhalomass

  integer :: nhalo
  real, allocatable :: radhalo(:),fourpirho(:),pothalo(:),mhalo(:)
  
 contains

subroutine initfixedhalo(pos,halofile)
 character(len=*) :: halofile
 real :: pos(3)
 integer :: i,ioerror
  
 if(allocated(radhalo)) then
  print*,' ** warning: fixedhalomod already initialized **'
  return
 endif 

 fixedhalopos=pos
 
 open(iounit,FILE=halofile, STATUS='OLD',iostat=ioerror)
 if(ioerror.NE.0) call Error(' stop -- error opening halofile: '//halofile)
 
 read(iounit,*) nhalo
 if(nhalo.LE.10.OR.nhalo.GT.10000000) call Error('nhalo of file: '//halofile)

 allocate(radhalo(nhalo),fourpirho(nhalo),mhalo(nhalo), &
         pothalo(nhalo))
 do i=1,nhalo
   read(iounit,*,iostat=ioerror) radhalo(i),fourpirho(i),mhalo(i),pothalo(i)
   if(ioerror.NE.0) call Error(' error reading halofile: '//halofile)
   if(i.GE.2) then
      if(radhalo(i).LE.radhalo(i-1)) &
         call Error(' ordering in halofile: '//halofile,i)
      if(mhalo(i).LT.mhalo(i-1)) mhalo(i)=mhalo(i-1)   
   endif
 enddo
 close(iounit)

 fixedhalomass=mhalo(nhalo)
 pothalo(1:nhalo)=pothalo(1:nhalo)-(fixedhalomass/radhalo(nhalo)+pothalo(nhalo))
 fourpirho=4*pi*radhalo**2*fourpirho 	
         
end subroutine

function fixedhalopot(pos) result(pot)
 real pos(3),pot,r,r2
 r2=sum((pos-fixedhalopos)**2)
 r=sqrt(r2)
 if(r.GE.radhalo(nhalo)) then
  pot=-fixedhalomass/r
 else
   if(r.LE.radhalo(2)) then
      pot=pothalo(2)-mhalo(2)/2/radhalo(2)* &
            (1-r2/radhalo(2)**2)
   else
      pot=linearlookup(nhalo,radhalo,pothalo,r) 
   endif
 endif
end function

function fixedhaloacc(pos) result(acc)
 real pos(3),acc(3),r,r2,m

 r2=sum((pos-fixedhalopos)**2)
 r=sqrt(r2)
 if(r.GE.radhalo(nhalo)) then
   acc=-(pos-fixedhalopos)*fixedhalomass/r2/r
 else
   if(r.LT.radhalo(2)) then
      r=radhalo(2)
      r2=r**2
   endif
   m=hermitelookup(nhalo,radhalo,mhalo,fourpirho,r)   
   acc=-(pos-fixedhalopos)*m/r2/r
 endif
end function

function fixedhalotidalfield(pos) result(tide)
 real pos(3),tide(6),r,r2,m,dm,dr3i,dr5i,dx,dy,dz

 dx=pos(1)-fixedhalopos(1)
 dy=pos(2)-fixedhalopos(2)
 dz=pos(3)-fixedhalopos(3)
 r2=dx**2+dy**2+dz**2
 r=sqrt(r2)
 if(r.GE.radhalo(nhalo)) then
   dr3i=fixedhalomass/(r2*r)
   dr5i=3*fixedhalomass/(r2*r2*r)
 else  
   if(r.LT.radhalo(2)) then
      dr3i=mhalo(2)/radhalo(2)**3
      dr5i=0
   else
      call hermitelookup_deriv(nhalo,radhalo,mhalo,fourpirho,r,m,dm)
      dr3i=m/(r2*r)
      dr5i=(3*m-dm*r)/(r2*r2*r)
   endif
 endif
 tide(1)=(-dr3i+dx**2*dr5i)
 tide(4)=(-dr3i+dy**2*dr5i)
 tide(6)=(-dr3i+dz**2*dr5i)
 tide(2)=(dx*dy*dr5i)
 tide(3)=(dx*dz*dr5i)
 tide(5)=(dy*dz*dr5i)
end function

subroutine endfixedhalo

 deallocate(radhalo,fourpirho,mhalo,pothalo)

end subroutine

! -------------------utility functions below

 function findbin(n,ylist,y) result (bin)
  integer :: n,bin
  real :: ylist(n),y,s

  if(ylist(1).LT.ylist(n)) then
   s=1
  else
   s=-1
  endif  
  
 bin=findbin_const(n,ylist,y)
 if( (bin.LE.1.AND.s*ylist(1).GE.s*y).OR. &
      (bin.GE.n+1.AND.s*ylist(n).LE.s*y).OR. &
      (s*ylist(bin-1).LE.s*y.AND.s*ylist(bin).GE.s*y)) return
 bin=findbin_bisect(n,ylist,y)
 if( (bin.LE.1.AND.s*ylist(1).GE.s*y).OR. &
      (bin.GE.n+1.AND.s*ylist(n).LE.s*y).OR. &
      (s*ylist(bin-1).LE.s*y.AND.s*ylist(bin).GE.s*y)) return
 call Error('findbin failure') 

 end function

 function findbin_const(n,ylist,y) result (bin)
  integer :: n,bin
  real :: ylist(n),y,s,dr

  if(ylist(1).LT.ylist(n)) then
   s=1
  else
   s=-1
  endif  
  
  if(s*y.LE.s*ylist(1)) then
   bin=0
   return
  endif
  if(s*y.GE.s*ylist(n)) then
   bin=n+1  
   return
  endif

  dr=(ylist(n)-ylist(1))/(n-1)
  bin=INT((y-ylist(1))/dr)+2 ! 2 because: top index+fortran start at 1

 end function
 
 function findbin_bisect(n,ylist,y) result(bin)
  integer :: n,bin,up,low
  real :: ylist(n),y,s
  
  if(ylist(1).LT.ylist(n)) then
   s=1
  else
   s=-1
  endif  
  
  if(s*y.LE.s*ylist(1)) then
   bin=0
   return
  endif
  if(s*y.GE.s*ylist(n)) then
   bin=n+1  
   return
  endif
  up=n
  low=1
  do while((up-low).GT.1)
   bin=(low+up)/2
   if(s*y.LT.s*ylist(bin)) then
    up=bin
   else
    low=bin
   endif   
  enddo 
  bin=up
  end function

 function interpolatelinear(bin,n,xlist,ylist,x) result(y)
  integer :: n,bin,up,low
  real :: xlist(n),ylist(n),x,y,u,s
  if(bin.LE.1) then
   y=ylist(1)
   return
  endif 
  if(bin.GT.n) then
   y=ylist(n)
   return
  endif 
  if(xlist(bin).EQ.xlist(bin-1)) then
   y=(ylist(bin-1)+ylist(bin))/2
   return
  endif
  u=(x-xlist(bin-1))/(xlist(bin)-xlist(bin-1))
  y=(1-u)*ylist(bin-1)+u*ylist(bin)  
  return
 end function

 subroutine interpolatecubic(bin,n,xlist,ylist,yderiv,x,y,d)
  integer :: bin,n
  real :: xlist(n),ylist(n),yderiv(n),x,y,d
  real :: u,dy,dx,y1,yd1,yd2  
  if(bin.LE.1) then
   y=ylist(1)
   d=0.
   return
  endif 
  if(bin.GT.n) then
   y=ylist(n)
   d=0.
   return
  endif 
  dx=xlist(bin)-xlist(bin-1)
  dy=ylist(bin)-ylist(bin-1)
  if(dx.EQ.0) then
   y=(ylist(bin-1)+ylist(bin))/2
   d=(yderiv(bin-1)+yderiv(bin))/2
   return
  endif
  y1=ylist(bin-1)
  yd2=yderiv(bin)
  yd1=yderiv(bin-1)
  u=(x-xlist(bin-1))/(xlist(bin)-xlist(bin-1))
  y=u**3*(-2*dy+dx*(yd1+yd2))+u**2*(3*dy-dx*(2*yd1+yd2))+dx*yd1*u+y1
  d=(3*u**2*(-2*dy+dx*(yd1+yd2))+2*u*(3*dy-dx*(2*yd1+yd2))+dx*yd1)/ &
      (xlist(bin)-xlist(bin-1))
 end subroutine

 function hermitelookup(n,xlist,ylist,yderiv,x) result(y)
  integer :: n,bin
  real :: x,y,xlist(n),ylist(n),yderiv(n),d
  bin=findbin(n,xlist,x)
  call interpolatecubic(bin,n,xlist,ylist,yderiv,x,y,d)
 end function

 subroutine hermitelookup_deriv(n,xlist,ylist,yderiv,x,y,d)
  integer :: n,bin
  real :: x,y,xlist(n),ylist(n),yderiv(n),d
  bin=findbin(n,xlist,x)
  call interpolatecubic(bin,n,xlist,ylist,yderiv,x,y,d)
 end subroutine

 function linearlookup(n,xlist,ylist,x) result(y)
  integer :: n,bin
  real :: x,y,xlist(n),ylist(n)
  bin=findbin(n,xlist,x)
  y=interpolatelinear(bin,n,xlist,ylist,x)
 end function

subroutine Error(string,value)
 character(len=*) :: string
 integer, optional :: value
 
 print*,'--fixedHalo Mod Error--'
 print*,string
 if(present(value)) print*,value
 stop

end subroutine  

end module

subroutine accpothalo(ppos,pacc,pphi,option)
 use fixedhalomod
 include 'globals.h'
 character*4 :: option
 real :: ppos(3),pacc(3),pphi
 if(option.NE.'pot ') pacc=pacc+fixedhaloacc(ppos)
 if(option.NE.'acc ') pphi=pphi+fixedhalopot(ppos)
end subroutine

subroutine halotidalfield(ppos,ptide)
  use fixedhalomod
 real :: ppos(3),ptide(6)
 ptide=ptide+fixedhalotidalfield(ppos)
end subroutine

subroutine testhalo
 use fixedhalomod
 real pos(3),pot,acc(3),tide(6),dens,r
 integer i
 
 call initfixedhalo((/0.,0.,0./),"/home/inti/code/stars/data/hk_96_2.halo")

10 read*,pos
pot=fixedhalopot(pos)
acc=fixedhaloacc(pos)
tide=fixedhalotidalfield(pos)
dens=-(tide(1)+tide(4)+tide(6))/4/3.14159265358979
 r=sqrt(sum(pos**2))
 write(*,'(6g16.8)'), r,pot,acc,dens

 call endfixedhalo
end subroutine
