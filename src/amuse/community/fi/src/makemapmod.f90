module makemapMod
 implicit none
 private
 public:: InitMap, project,EndMap,EraseMap
 public:: nx,ny,pic,opdepth,extinction

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

subroutine EraseMap()

 pic=0.
 opdepth=1.

end subroutine

 subroutine EndMap()
 
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
  real ::ppos(3),psize,pmass,popac,d2
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
    d2=sum(ppos**2)   
    psize=psize*camdis/ppos(3)
    ppos(1)=ppos(1)*camdis/ppos(3)
    ppos(2)=ppos(2)*camdis/ppos(3)  
    pmass=pmass*camdis**2/d2
    popac=popac*camdis**2/d2
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


end module
