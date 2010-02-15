module makemapMod
 implicit none
 private
 public:: InitMap, project,EndMap

 real,save :: xsize,ysize
 integer,save :: nx,ny
 integer,save :: maptype
 real,save :: dx,zmin,camdis
 real,save,dimension(3) :: lookat,xax,yax,zax,direction,from

 real, allocatable,save :: opdepth(:,:),pic(:,:)
 
 real, parameter :: sqrt2i=0.70710678118654752440
 real,parameter :: sqrt2=1.4142135623730950488
 real,parameter :: maxx=sqrt2
 real,parameter :: Pi=3.141592654

 integer,save :: projection, extinction

 integer,parameter :: Nstarmap=1000
 real :: mi(0:Nstarmap)

 type pparticle_type
  integer :: magic_number
  real :: pos(2)
  real :: size
  real :: mass
  real :: opac
 end type

 contains

 subroutine InitMap(imsize,width,focus,pointing,upvector,proj,ext,zm)
  integer, optional :: imsize(2),proj,ext
  real, optional :: width,focus(3),pointing(3),upvector(3),zm
  real :: up(3)
  integer :: test

  call InitStarMap
         
  if(present(proj)) then
   projection=proj
  else
   print*,' projection? (parallel=0, perspective=1)'
   read*, projection
  endif

  if(present(ext)) then
   extinction=ext
  else
   print*,'extinction (0=no,1=yes)?'
   read*,extinction
  endif

  zmin=-1.e30
  if(projection.ge.1) then
   if(present(zm)) then
    zmin=zm
   else
    print*,' zmin?'
    read*,zmin
   endif
  endif 
  
  if(present(imsize)) then
   nx=imsize(1)
   ny=imsize(2)
  else
   print*,' image pixel size (WxH)?'
   read*,nx,ny
  endif
  
  if(present(width)) then
   xsize=width
  else
   if(projection.eq.0) print*,' image width?'
   if(projection.eq.1) print*,' image angle?'    
   read*,xsize
  endif
  
  if(present(focus)) then
   lookat=focus
  else
   print*,' look at?'
   read*, lookat
  endif
  
  if(present(pointing)) then
   direction=pointing
  else
   if(projection.eq.0) print*,' direction?'
   if(projection.eq.1) print*,' from?'    
   read*,direction
  endif
  
  if(present(upvector)) then
   up=upvector
  else 
   print*,' upvector'
   read*,up
  endif
  
  allocate(pic(nx,ny),STAT=test)  
  if(test.NE.0) call Error(1)
   pic=0.
  
  allocate(opdepth(nx,ny),STAT=test)  
  if(test.NE.0) call Error(1)
   opdepth=1.
   
  if(projection.eq.0) from=lookat-direction  
  if(projection.eq.1) then
   from=direction
   direction=lookat-from
  endif
  camdis=norm(direction)
  if(projection.eq.1) xsize=2*camdis*tan(xsize/2./180.*Pi)
  call normalize(direction)
  zax=direction
  call outer(xax,direction,up)
  call normalize(xax)
  call outer(yax,xax,direction)
  call normalize(yax)

  dx=xsize/(nx-1)
  ysize=dx*(ny-1)
 end subroutine

 subroutine endMap(filename,tauname)
  character, optional :: filename*80,tauname*80
  character :: file*80
  integer :: naxes(3)
 
  if(present(filename)) then
   file=filename
  else
   print*,' uit filenaam?'
   read*, file
  endif
 
  naxes(1)=nx
  naxes(2)=ny
  if(file(1:5).NE.'XXXXX') call writefits(file,2,naxes,pic)
 
  if(extinction.eq.1) then
   
  if(present(tauname)) then
   file=tauname
  else
   print*,' tau filenaam?'
   read*, file
  endif
 
  naxes(1)=nx
  naxes(2)=ny
  if(file(1:5).NE.'XXXXX') call writefits(file,2,naxes,opdepth)

  endif  
 
  deallocate(pic,opdepth)

 end subroutine

 subroutine Error(i)
  integer i
  
  print*,' error making map',i
  stop
 end subroutine

 subroutine project(n,pos,sizes,mass,opac,pindex)
  integer,intent(in) :: n  ! arrsize necessary because pos() stored wrong way?
  real,intent(in) :: pos(:,:),sizes(:),mass(:)
  real,intent(in),optional :: opac(:)
  integer, intent(in), optional :: pindex(:)
  real ::ppos(3),psize,pmass,popac
  real, allocatable :: scratch(:,:),zdist(:)
  integer, allocatable :: order(:)
  integer i,j,test,xmin,xmax,ymin,ymax
  integer mythread,totalthread,low(2),up(2)
!  integer mythread,totalthread,low,up
!$  integer omp_get_thread_num,omp_get_num_threads     
  integer :: magic_number   ! used to get consistent 'random' numbers 


 if((extinction.EQ.1).AND.(.not.present(opac))) then
  call error(99)
 endif 

 allocate(zdist(n),order(n),stat=test)
 if(test.ne.0) call Error(5)

 if(extinction.EQ.1) then
  do i=1,n
   call camcoor(pos(i,1:3),ppos)
   zdist(i)=ppos(3)
  enddo
  call mrgrnk(n,zdist,order)
 else
  do i=1,n
   order(i)=i
  enddo
 endif

!$omp parallel private(scratch,i,j,up,low,test,mythread,totalthread,ppos, &
!$omp pmass,popac,psize,xmin,xmax,ymin,ymax,magic_number)   
 allocate(scratch(min(nx,ny),min(nx,ny)),STAT=test)
 if(test.ne.0) call Error(5)
 
 mythread=0
 totalthread=1
!$ mythread=omp_get_thread_num()
!$ totalthread=omp_get_num_threads()
! low=(mythread*nx)/totalthread+1
! up=((mythread+1)*nx)/totalthread
 call verdeel(mythread,totalthread,nx,ny,low,up) 
 do j=1,n
  i=order(j)
  magic_number=i
  if(present(pindex)) magic_number=pindex(i)
  call camcoor(pos(i,1:3),ppos)
  pmass=0
  popac=0
  if(ppos(3).ge.zmin) then
   pmass=mass(i)
   psize=sizes(i)
   if(present(opac)) popac=opac(i)/dx**2
   if(projection.eq.1) then   
    psize=psize*camdis/ppos(3)
    ppos(1)=ppos(1)*camdis/ppos(3)
    ppos(2)=ppos(2)*camdis/ppos(3)  
    pmass=pmass*camdis**2/ppos(3)**2
    popac=popac*camdis**2/ppos(3)**2
   endif
   xmin=FLOOR((ppos(1)-maxx*psize+xsize/2.)/dx)+1
   xmax=CEILING((ppos(1)+maxx*psize+xsize/2.)/dx)+1
   ymin=FLOOR((ppos(2)-maxx*psize+ysize/2.)/dx)+1
   ymax=CEILING((ppos(2)+maxx*psize+ysize/2.)/dx)+1
   if(xmin.GT.up(1).OR.xmax.LT.low(1).OR. &
      ymin.GT.up(2).OR.ymax.LT.low(2)) then
!   if(xmin.GT.up.OR.xmax.LT.low.OR. &
!      ymin.GT.ny.OR.ymax.LT.1) then
      pmass=0   
      popac=0
   endif   
  endif
  if(extinction.EQ.1) then
   if(pmass.GT.0.OR.popac.gt.0) call projectopacity(low,up,pic,opdepth, &
                                  scratch,ppos,psize,pmass,popac,magic_number)
  else 
   if(pmass.gt.0) call projectparticle(low,up,pic,scratch,ppos,psize,pmass)
  endif
 enddo

 deallocate(scratch) 
!$omp end parallel 
 deallocate(zdist,order)
 
 end subroutine

 subroutine camcoor(ipos,opos)
  real, intent(in) :: ipos(3)
  real, intent(out) :: opos(3)
  
  opos(1)=dot_product(xax,ipos-from)
  opos(2)=dot_product(yax,ipos-from)
  opos(3)=dot_product(zax,ipos-from)

 end subroutine camcoor

 subroutine projectparticle(low,up,img,temp,pos,hsmooth,mass)
  real :: pos(3),hsmooth,mass,offset(2),img(:,:),temp(:,:)
  integer low(2),up(2),i,j,size
!  integer low,up,i,j,size
    
  size=CEILING(2*maxx*hsmooth/dx)+2
  if(2*maxx*hsmooth.LE.sqrt2*dx) size=3
    
  if(size.le.nx.AND.size.le.ny) then

   call positionp(i,j,offset,pos,hsmooth)

   call makepimage(temp,offset,size,hsmooth)
  
   call addpimage(low,up,img,i,j,size,temp,mass)
   
  endif

 end subroutine

 subroutine projectopacity(low,up,img,opdepth,temp,pos,hsmooth,mass,opac,magic_number)
   real :: pos(3),hsmooth,mass,offset(2),img(:,:),temp(:,:),opdepth(:,:),opac
  integer low(2),up(2),i,j,size
!  integer low,up,i,j,size
  integer :: magic_number 


  size=CEILING(2*maxx*hsmooth/dx)+2
  if(2*maxx*hsmooth.LE.sqrt2*dx) size=3
!  if(2*maxx*hsmooth/dx.LT.1.5) opac=opac*(2*maxx*hsmooth/dx/1.5)**2
    
  if(size.le.nx.AND.size.le.ny) then

   call positionp(i,j,offset,pos,hsmooth)

   if(opac.EQ.0) then
    call makesimage(temp,offset,size,hsmooth,magic_number)
   else   
    call makepimage(temp,offset,size,hsmooth)
   endif
   
   call addpopacimage(low,up,img,temp,opdepth,i,j,size,mass,opac)
   
  endif
   
 end subroutine

 subroutine addpimage(low,up,img,i,j,size,pimage,mass)
  integer low(2),up(2),i,j,size,n1,n2,m1,m2
!  integer low,up,i,j,size,n1,n2,m1,m2
  real mass,pimage(:,:),img(:,:)
  
 n1=max(1,low(1)-i); n2=min(size,up(1)-i)
 m1=max(1,low(2)-j); m2=min(size,up(2)-j)
! n1=max(1,low-i); n2=min(size,up-i)
! m1=max(1,1-j); m2=min(size,ny-j)
 img(i+n1:i+n2,j+m1:j+m2)=img(i+n1:i+n2,j+m1:j+m2)+mass*pimage(n1:n2,m1:m2)
  
 end subroutine

 subroutine addpopacimage(low,up,img,pimage,tau,i,j,size,mass,opac)
  integer low(2),up(2),i,j,size,n1,n2,m1,m2
!  integer low,up,i,j,size,n1,n2,m1,m2
  real mass,opac,pimage(:,:),img(:,:),tau(:,:)
  
 n1=max(1,low(1)-i); n2=min(size,up(1)-i)
 m1=max(1,low(2)-j); m2=min(size,up(2)-j)
! n1=max(1,low-i); n2=min(size,up-i)
! m1=max(1,1-j); m2=min(size,ny-j)

! if(opac.gt.0) tau(i+n1:i+n2,j+m1:j+m2)=tau(i+n1:i+n2,j+m1:j+m2)*exp(-opac/2*pimage(n1:n2,m1:m2))
! if(opac.gt.0) tau(i+n1:i+n2,j+m1:j+m2)=tau(i+n1:i+n2,j+m1:j+m2)*exp(-opac/2*pimage(n1:n2,m1:m2))

  if(mass.gt.0.and.opac.gt.0) then
   tau(i+n1:i+n2,j+m1:j+m2)=tau(i+n1:i+n2,j+m1:j+m2)*(1.-min(1.,opac/2*pimage(n1:n2,m1:m2)))
   img(i+n1:i+n2,j+m1:j+m2)=img(i+n1:i+n2,j+m1:j+m2)+mass*pimage(n1:n2,m1:m2)*tau(i+n1:i+n2,j+m1:j+m2)
   tau(i+n1:i+n2,j+m1:j+m2)=tau(i+n1:i+n2,j+m1:j+m2)*(1.-min(1.,opac/2*pimage(n1:n2,m1:m2)))
   return
  endif
  if(mass.gt.0) img(i+n1:i+n2,j+m1:j+m2)=img(i+n1:i+n2,j+m1:j+m2)+mass*pimage(n1:n2,m1:m2)*tau(i+n1:i+n2,j+m1:j+m2)
  if(opac.gt.0) tau(i+n1:i+n2,j+m1:j+m2)=tau(i+n1:i+n2,j+m1:j+m2)*(1.-min(1.,opac*pimage(n1:n2,m1:m2)))
  
 end subroutine


 subroutine makepimage(pimage,offset,size,hsmooth)
  integer :: size
  real :: offset(2),hsmooth,pimage(:,:)
  real :: r,di,dj,dens,norm,hinv
  integer :: i,j
  
  if(size.eq.3) then
   i=2;j=2   
   if(offset(1)+maxx*hsmooth.LT.dx/2) i=1
   if(offset(2)+maxx*hsmooth.LT.dx/2) j=1
   if(offset(1)+maxx*hsmooth.GT.3*dx/2.) i=3
   if(offset(2)+maxx*hsmooth.GT.3*dx/2.) j=3
   pimage(1:3,1:3)=0
   pimage(i,j)=1
   return
  endif 

  norm=0.
  hinv=1/hsmooth
  do j=1,size
   dj=(-maxx*hsmooth-offset(2)+(j-1)*dx)**2
   do i=1,size
      di=-maxx*hsmooth-offset(1)+(i-1)*dx    
      r=SQRT(di**2+dj)*hinv
      dens=spline(r)
      norm=norm+dens
      pimage(i,j)=dens
   enddo
  enddo  
  if(norm.EQ.0)  call Error(6)
  pimage(1:size,1:size)=pimage(1:size,1:size)/norm

 end subroutine

 subroutine positionp(i,j,offset,pos,hsmooth)
  integer :: i,j
  real :: offset(2),pos(3),hsmooth
  real x,y
  x=pos(1) ; y=pos(2)
    
  i=FLOOR((x-maxx*hsmooth+xsize/2.)/dx)-1
  j=FLOOR((y-maxx*hsmooth+ysize/2.)/dx)-1
  offset(1)=(x-maxx*hsmooth+xsize/2.)-(i+1)*dx
  offset(2)=(y-maxx*hsmooth+ysize/2.)-(j+1)*dx
 end subroutine

subroutine InitStarMap
 real :: dr,df,du,m(0:Nstarmap)
 integer :: n,i,j
 
 n=Nstarmap
 dr=sqrt(2.)/n
 df=1./n
 
 m(0)=0
 do i=1,n
  m(i)=m(i-1)+.5*dr*(dr*(i-1)*spline(dr*(i-1))+dr*i*spline(dr*i))
 enddo

 mi(0)=0.
 mi(n)=dr*n
 do i=1,n-1
  j=0
  do while(m(j)/m(n).lt.df*i)
   j=j+1
  enddo
  du=(df*i*m(n)-m(j-1))/(m(j)-m(j-1))
  mi(i)=(1-du)*dr*(j-1)+du*dr*j
 enddo 

end subroutine

subroutine makesimage(pimage,offset,size,hsmooth,magic_number)
  integer :: size
  real :: offset(2),hsmooth,pimage(:,:)
  integer :: k,i,j,ix,iy,im
  real :: rr,x(2),dr,mir,bright,norm,pRND
  integer, parameter :: ndepth=4
  real, parameter :: brightfac=1.5
  integer magic_number

  if(size.eq.3) then
   i=2;j=2   
   if(offset(1)+maxx*hsmooth.LT.dx/2) i=1
   if(offset(2)+maxx*hsmooth.LT.dx/2) j=1
   if(offset(1)+maxx*hsmooth.GT.3*dx/2.) i=3
   if(offset(2)+maxx*hsmooth.GT.3*dx/2.) j=3
   pimage(1:3,1:3)=0
   pimage(i,j)=1
   return
  endif 

  im=0
  norm=0.
  pimage(1:size,1:size)=0.
  do k=0,ndepth
   bright=1./brightfac**k
   do i=1,2**k
    rr=2
    do while(rr.gt.1)
      x(1)=pRND(magic_number+2*im)
      x(2)=pRND(3*magic_number+2*im+1)
      x=x*2-1
      rr=sum(x**2)
      im=im+1
    enddo
    j=rr*Nstarmap
    dr=rr*Nstarmap-j
    mir=(1-dr)*mi(j)+dr*mi(j+1)
    x=hsmooth*x*mir/sqrt(rr)
    x=x+offset+maxx*hsmooth
    ix=NINT(x(1)/dx)+1
    iy=NINT(x(2)/dx)+1
    if(ix.lt.1.or.ix.gt.size.or.iy.lt.1.or.iy.gt.size) then
!$omp critical
      print*,'project out of range'
      print*, x,ix,iy
      stop
!$omp end critical      
    endif 
    pimage(ix,iy)=pimage(ix,iy)+bright
    norm=norm+bright
   enddo
  enddo
  if(norm.le.0) then
   print*,norm
   stop
  endif  
  pimage(1:size,1:size)=pimage(1:size,1:size)/norm
end subroutine




 function spline(r) result(w)
  real :: r,w
  w=0
!  if(r.GT.2) return
  if(r.GT.sqrt2) return
  w=1-3/2.*r**2+sqrt2i*r**3

!  if(r.LT.1) then
!   w=1-3/2.*r**2+3/4.*r**3
!  else
!   w=1/4.*(2-r)**3
!  endif
!   w=w*10/(7*Pi)
 end function
 
 real function norm(n)
  real n(3)
   norm=SQRT(sum(n**2))
 end function

 subroutine normalize(n)
  real n(3),v
   v=norm(n)
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


subroutine verdeel(i,n,nx,ny,low,up)
 integer:: n,deler1,deler2,nx,ny,i,ix,iy
 integer :: low(2),up(2)
 
 deler1=n**(0.5)
 do while(mod(n,deler1).NE.0.AND.deler1.GT.1)
  deler1=deler1-1
 enddo
 deler2=n/deler1

  iy=mod(i,deler1)
  ix=i/deler1

  low(1)=1+(nx*ix)/deler2
  up(1)=(nx*(ix+1))/deler2
  low(2)=1+(ny*iy)/deler1
  up(2)=(ny*(iy+1))/deler1

end subroutine


 
 subroutine writefits(filename,naxis,naxes,pic)
      integer status,unit,blocksize,bitpix,naxis,naxes(3)
      integer i,j,k,group,fpixel,nelements
      real pic(:,:)
      character filename*80
      logical simple,extend

      status=0

!      filename='HImap.fits'
      
      call deletefile(filename,status)
      call ftgiou(unit,status)

      blocksize=1
      call ftinit(unit,filename,blocksize,status)

      simple=.true.
      bitpix=-64
!      naxis=2
!      naxes(1)=nx
!      naxes(2)=ny
      extend=.true.

      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

      group=1
      fpixel=1
      nelements=1
      do i=1,naxis
       nelements=naxes(i)*nelements
      enddo
      call ftpprd(unit,group,fpixel,nelements,pic,status)

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






