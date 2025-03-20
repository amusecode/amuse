module HICubeMod
 implicit none
 private
 public:: InitHICube, projectparticle,EndCube

 real,save :: xsize,ysize,vsize
 integer,save :: nx,ny,nv
 real,save :: dx,dy,dv
 real,save :: centerv
 real,save,dimension(3) :: lookat,xax,yax,direction

 real, allocatable,save :: cube(:,:,:),scratch(:,:,:)
 
 real,parameter :: maxx=2,maxy=2,maxv=3
 real,parameter :: Pi=3.141592654

 contains

 subroutine InitHICube
  real :: up(3)
  integer :: test
  
  print*,' image pixel size (WxH), velocity channels?'
  read*,nx,ny,nv
  print*,' image width, velocity range?'
  read*,xsize,vsize
  print*,' look at?'
  read*, lookat
  print*,' direction?'
  read*,direction
  print*,' upvector'
  read*,up
  print*,' system velocity?'
  read*,centerv

  allocate(cube(nx,ny,nv),scratch(nx,ny,nv),STAT=test)  
  if(test.NE.0) call cubeError(1)
   cube=0.

  
  dx=(xsize)/(nx-1)
  dy=dx
  ysize=dx*(ny-1)
  dv=(vsize)/(nv-1)
  
  print*,dx,dy,dv
  call normalize(direction)
  call outer(xax,direction,up)
  call normalize(xax)
  call outer(yax,xax,direction)
  call normalize(yax)
  print*,xax,yax

 end subroutine

 subroutine endCube

 call writefits

 deallocate(cube,scratch)
 
 end subroutine

 subroutine cubeError(i)
  integer i
  
  print*,i
  stop
 end subroutine

 subroutine projectparticle(pos,hsmooth,vel,dopplerwidth,mass)
  real :: pos(3),hsmooth,vel(3),dopplerwidth,mass,offset(3)
  integer i,j,k,isize,jsize,ksize
!  real, allocatable :: pimage(:,:,:)
  
  isize=2*maxx*hsmooth/dx+1
  jsize=2*maxy*hsmooth/dy+1
  ksize=2*maxv*dopplerwidth/dv+1
  if(isize.le.nx.and.jsize.le.ny.and.ksize.le.nv) then
  
  call positionp(i,j,k,offset,pos,hsmooth,vel,dopplerwidth)

!  allocate(pimage(isize,jsize,ksize))

  call makepimage(scratch,offset,isize,jsize,ksize,hsmooth,dopplerwidth)
  
  call addpimage(i,j,k,isize,jsize,ksize,scratch,mass)

  endif

!  deallocate(pimage)

 end subroutine


 subroutine addpimage(i,j,k,isize,jsize,ksize,pimage,mass)
  integer i,j,k,isize,jsize,ksize,n,m,l
  real mass,pimage(isize,jsize,ksize)
  
!  do n=1,isize
!   do m=1,jsize
!    do l=1,ksize
  do l=1,ksize
   do m=1,jsize
    do n=1,isize
    if((n+i.GE.1.AND.n+i.LE.nx).AND. &
       (m+j.GE.1.AND.m+j.LE.ny).AND. &
       (k+l.GE.1.AND.k+l.LE.nv)) then
       cube(n+i,m+j,l+k)=cube(n+i,m+j,l+k)+pimage(n,m,l)*mass
    endif    
    enddo
   enddo
  enddo  
  
 end subroutine

 subroutine makepimage(pimage,offset,isize,jsize,ksize,hsmooth,dw)
  integer :: isize,jsize,ksize
  real :: offset(3),hsmooth,dw,pimage(isize,jsize,ksize)
  real :: r,di,dj,dk,dens,norm,hinv,dwinv
  integer :: i,j,k
  
  hinv=1/hsmooth
  dwinv=1/dw
  norm=0.
  if(isize.eq.1.or.jsize.eq.1) then
   pimage(1:isize,1:jsize,1:ksize)=0.
   do k=1,ksize
    dk=-maxv*dw-offset(3)+(k-1)*dv
    dens=exp(-(dk*dwinv)**2/2)
    norm=norm+dens
    pimage((isize+1)/2,(jsize+1)/2,k)=dens
   enddo  
  else 
!  do i=1,isize
!   do j=1,jsize
!    do k=1,ksize
 do k=1,ksize
  do j=1,jsize
    do i=1,isize
      di=-maxx*hsmooth-offset(1)+(i-1)*dx    
      dj=-maxy*hsmooth-offset(2)+(j-1)*dy
      dk=-maxv*dw-offset(3)+(k-1)*dv
      r=SQRT(di**2+dj**2)*hinv
      dens=spline(r)*exp(-(dk*dwinv)**2/2)
      norm=norm+dens
      pimage(i,j,k)=dens
    enddo
   enddo
  enddo  
  endif
  
  if(norm.ne.0) then 
!  do i=1,isize
!   do j=1,jsize
!    do k=1,ksize
    do k=1,ksize
   do j=1,jsize
  do i=1,isize
      pimage(i,j,k)=pimage(i,j,k)/norm
    enddo
   enddo
  enddo  
  endif
 
 end subroutine

 subroutine positionp(i,j,k,offset,pos,hsmooth,vel,dw)
  integer :: i,j,k
  real :: offset(3),pos(3),hsmooth,vel(3),dw
  real x,y,v
  
  x=dot_product(xax,pos-lookat)
  y=dot_product(yax,pos-lookat)
  v=dot_product(direction,vel)-centerv
  
  i=(x-maxx*hsmooth+xsize/2)/dx
  j=(y-maxy*hsmooth+ysize/2)/dy
  k=(v-maxv*dw+vsize/2)/dv
  
  offset(1)=(x-maxx*hsmooth+xsize/2)-i*dx
  offset(2)=(y-maxy*hsmooth+ysize/2)-j*dy
  offset(3)=(v-maxv*dw+vsize/2)-k*dv

 end subroutine

 function spline(r) result(w)
  real :: r,w
  w=0
  if(r.GT.2) return
  if(r.LT.1) then
   w=1-3/2.*r**2+3/4.*r**3
  else
   w=1/4.*(2-r)**3
  endif
   w=w*10/(7*Pi)
 end function
 
 subroutine normalize(n)
  real n(3),v
   v=SQRT(n(1)**2+n(2)**2+n(3)**2)
   n(1)=n(1)/v
   n(2)=n(2)/v
   n(3)=n(3)/v
 end subroutine
	
 subroutine outer(o,v,w)
  real v(3),w(3),o(3)
  o(1)=v(2)*w(3)-w(2)*v(3)
  o(2)=v(3)*w(1)-w(3)*v(1)
  o(3)=v(1)*w(2)-w(1)*v(2)
 end subroutine
 
 subroutine writefits
      integer status,unit,blocksize,bitpix,naxis,naxes(3)
      integer i,j,k,group,fpixel,nelements
      character filename*80
      logical simple,extend

      status=0

      filename='HICube.fits'
      
      call deletefile(filename,status)
      call ftgiou(unit,status)

      blocksize=1
      call ftinit(unit,filename,blocksize,status)

      simple=.true.
      bitpix=-64
      naxis=3
      naxes(1)=nx
      naxes(2)=ny
      naxes(3)=nv
      extend=.true.

      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

      group=1
      fpixel=1
      nelements=naxes(1)*naxes(2)*naxes(3)
      call ftpprd(unit,group,fpixel,nelements,cube,status)

!      call ftpkyj(unit,'EXPOSURE',1500,'Total Exposure Time',status)

      call ftclos(unit, status)
      call ftfiou(unit, status)
      
      if (status .gt. 0)call printerror(status)
      end subroutine
          
      subroutine printerror(status)


      integer status
      character errtext*30,errmessage*80

      if (status .le. 0)return
      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
      enddo
      end subroutine
      
      
      subroutine deletefile(filename,status)


      integer status,unit,blocksize
      character*(*) filename

      if (status .gt. 0)return

      call ftgiou(unit,status)

      call ftopen(unit,filename,1,blocksize,status)

      if (status .eq. 0)then
          call ftdelt(unit,status)
      else if (status .eq. 103)then
          status=0
          call ftcmsg
      else
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if

      call ftfiou(unit, status)
      end subroutine
      
end module


program makecube
 use HICubeMod
 include 'globals.h'
 real gpos(3),hsm,gvel(3),dw,himass
 integer i,mode
 character*24 filenaam
 call initmem(nbodsmax,nsphmax,ncells) 
 
 CALL set_parameters(0)

 print*,'filenaam?'
 read*,filenaam
 CALL readbods(filenaam)

 call heattabel        
 CALL initpars

 IF(periodic) CALL initbc

 CALL initnewstar

 CALL postprocessread
 
 call InitHICube

 print*,' map mode( 0=HI,1=H2)?'
 read*, mode

if(mode.eq.0) then
 do i=1,nsph
  gpos=pos(i,1:3)
  gvel=vel(i,1:3)*velscale/1.e5
  hsm=hsmooth(i)
  dw=SQRT(kboltz/amu*temperat(i))/1.e5
  himass=unitm_in_msun*mass(i)*max(0.,1-elecfrac(i))*max(0.,1.-h2frac(i))
  if( himass.gt.0) call projectparticle(gpos,hsm,gvel,dw,himass)
 enddo
else
 do i=1,nsph
  gpos=pos(i,1:3)
  gvel=vel(i,1:3)*velscale/1.e5
  hsm=hsmooth(i)
  dw=SQRT(kboltz/amu/2*temperat(i))/1.e5
  himass=unitm_in_msun*mass(i)*max(0.,1-elecfrac(i))*h2frac(i)
  if( himass.gt.0) call projectparticle(gpos,hsm,gvel,dw,himass)
 enddo
endif

 
 call EndCube
 
end program
