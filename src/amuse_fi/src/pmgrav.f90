module pmgravMod
 implicit none
 private
 public :: pmgravInit,pmgravEnd,pmgrav,pmgravaccpot,dmeshQ,softQ,rcutQ, & 
            ficorrect,forcecorrect
 real,parameter ::  PI=3.14159265358979323846
 integer :: ngrid            ! treePM grid size
 real  :: rsoft=.75          ! treePM softening parameter (cells)
 real :: rcut=4.5            ! cutoff radius (for tree) (*rsoft)
 real :: dmesh
 real :: boxmin(3)=0.,boxmax(3)=0.
 
 real, allocatable :: pot(:,:,:)
 real, allocatable :: accx(:,:,:),accy(:,:,:),accz(:,:,:)
 complex, allocatable :: fpot(:,:,:),fkernel(:,:,:)

 integer*8 :: plan,iplan

 logical, parameter :: periodic=.false.
   
 contains
 
 subroutine pmgravInit(n,bmin,bmax,rs)
 include 'fftw3.f'
  integer :: n,test,nthreads,i(3),i1(3)
  real :: bmin(3),bmax(3),bc(3),bs,dx(3)
  real, optional :: rs
!$  integer,external :: omp_get_max_threads
  
  ngrid=n
  bs=maxval(bmax-bmin)
  if(bs.LE.0) call Error("zero box size")
  
  if(periodic) dmesh=bs/ngrid
  if(.not.periodic) dmesh=bs/(ngrid/2-6)
  bc=0.5*(boxmax+boxmin)
  if(periodic) then
   boxmin=bc-dmesh*ngrid/2.
   boxmax=bc+dmesh*ngrid/2.
  else
   boxmin=bc-dmesh*ngrid/4.
   boxmax=bc+dmesh*3*ngrid/4.
  endif
  if(present(rs)) rsoft=rs  


  allocate(pot(ngrid,ngrid,ngrid), & 
	      fkernel(ngrid/2+1,ngrid,ngrid), &
	      fpot(ngrid/2+1,ngrid,ngrid),stat=test)
  if(test.NE.0) call Error("memory allocation failure 1") 

  if(.not.periodic) &
  allocate(accx(ngrid/2-2,ngrid/2-2,ngrid/2-2), & 
	   accy(ngrid/2-2,ngrid/2-2,ngrid/2-2), &
	   accz(ngrid/2-2,ngrid/2-2,ngrid/2-2), stat=test)
  if(periodic) &
  allocate(accx(ngrid,ngrid,ngrid), & 
	   accy(ngrid,ngrid,ngrid), &
	   accz(ngrid,ngrid,ngrid), stat=test)	   
  if(test.NE.0) call Error("memory allocation failure 2") 

! fftw init  
 nthreads=1
!$ nthreads=omp_get_max_threads()
!$ call dfftw_init_threads
!$ call dfftw_plan_with_nthreads(nthreads)

 call dfftw_plan_dft_r2c_3d(plan,n,n,n,pot,fpot,FFTW_ESTIMATE)
 call dfftw_plan_dft_c2r_3d(iplan,n,n,n,fpot,pot,FFTW_ESTIMATE)

 call makekernel

 end subroutine 
 
 subroutine pmgravEnd
 
  if(allocated(pot)) deallocate(pot,accx,accy,accz,fpot,fkernel)
  
! fftw stuff  
 call dfftw_destroy_plan(plan) 
 call dfftw_destroy_plan(iplan) 
!$ call dfftw_cleanup_threads
 
 end subroutine
 
 subroutine makeKernel
  integer :: i,j,k
  real :: x,y,z,fx,fy,fz,fac,kx,ky,kz,r,u,k2

  pot=0.
  fkernel=0.
  
!$omp parallel do  private(i,j,k,x,y,z,fx,fy,fz,fac,kx,ky,kz,r,u,k2) &
!$omp shared(pot,rsoft,ngrid)
  do k=1,ngrid
   do j=1,ngrid
    do i=1,ngrid
      x=(i-1.)/ngrid
      if(x.GE.0.5) x=x-1.
      y=(j-1.)/ngrid
      if(y.GE.0.5) y=y-1.
      z=(k-1.)/ngrid
      if(z.GE.0.5) z=z-1.
      r=sqrt(x**2+y**2+z**2)
      u=0.5*r/rsoft*ngrid
      fac=1-erfc(u)
      if(r.GT.0) then 
       pot(i,j,k)=-fac/r
      else
       pot(i,j,k)=-1/sqrt(PI)/rsoft*ngrid
      endif
    enddo
   enddo
  enddo  
  call dfftw_execute(plan)

!$omp parallel do  private(i,j,k,x,y,z,fx,fy,fz,fac,kx,ky,kz,r,u,k2) &
!$omp shared(fpot,fkernel,rsoft,ngrid)
 do k=1,ngrid
  do j=1,ngrid
   do i=1,ngrid/2+1
    kx=(i-1)
    if(kx.GT.ngrid/2) kx=(i-1)-ngrid
    ky=(j-1)
    if(ky.GT.ngrid/2) ky=(j-1)-ngrid
    kz=(k-1)
    if(kz.GT.ngrid/2) kz=(k-1)-ngrid
    k2=kx**2+ky**2+kz**2
    if(k2.GT.0) then
      fx=1;fy=1;fz=1
      if(kx.NE.0) fx=sin(PI*kx/ngrid)*ngrid/PI/kx
      if(ky.NE.0) fy=sin(PI*ky/ngrid)*ngrid/PI/ky
      if(kz.NE.0) fz=sin(PI*kz/ngrid)*ngrid/PI/kz
      fac=(fx*fy*fz)**(-4)
      fkernel(i,j,k)=fpot(i,j,k)*fac
    endif
   enddo
  enddo
 enddo   
 fkernel(1,1,1)=fpot(1,1,1) 
 end subroutine
 
 
 subroutine pmgrav(n,mass,pos)
  integer, intent(in) :: n
  integer :: k
  real,intent(in) :: mass(:), pos(:,:)

 if(.not.allocated(pot)) call Error(" pmgravInit not called")

 pot=0.
 call cic(n,mass,pos) 
 call dfftw_execute(plan)

!$omp parallel do shared(fpot,fkernel,dmesh,ngrid)
 do k=1,ngrid
 fpot(1:ngrid/2+1,1:ngrid,k)=fpot(1:ngrid/2+1,1:ngrid,k)* &
      fkernel(1:ngrid/2+1,1:ngrid,k)/REAL(ngrid)**4/dmesh
 enddo
 call dfftw_execute(iplan)
 
 call setacc
 
 end subroutine
 
 subroutine setacc
  integer :: i,j,k
!$omp parallel do private( i,j,k) &
!$omp shared(pot,accx,accy,accz)
 do k=3,ngrid/2-2
  do j=3,ngrid/2-2
   do i=3,ngrid/2-2
 accx(i,j,k)=-1/dmesh* &
  (pot(i-2,j,k)/12-2*pot(i-1,j,k)/3+2*pot(i+1,j,k)/3-pot(i+2,j,k)/12)
 accy(i,j,k)=-1/dmesh* &
  (pot(i,j-2,k)/12-2*pot(i,j-1,k)/3+2*pot(i,j+1,k)/3-pot(i,j+2,k)/12)
 accz(i,j,k)=-1/dmesh* &
  (pot(i,j,k-2)/12-2*pot(i,j,k-1)/3+2*pot(i,j,k+1)/3-pot(i,j,k+2)/12)
    enddo
   enddo
  enddo  
 end subroutine

 subroutine setacc2
  integer :: i,j,k,i1,i2,i22,i11,j1,j2,j11,j22,k1,k2,k11,k22
!$omp parallel do private( i,j,k,i1,i2,i22,i11,j1,j2,j11,j22,k1,k2,k11,k22) &
!$omp shared(pot,accx,accy,accz)
 do k=1,ngrid
  do j=1,ngrid
   do i=1,ngrid
 i2=MODULO(i-3,ngrid)+1
 i1=MODULO(i-2,ngrid)+1
 i11=MODULO(i,ngrid)+1
 i22=MODULO(i+1,ngrid)+1
 j2=MODULO(j-3,ngrid)+1
 j1=MODULO(j-2,ngrid)+1
 j11=MODULO(j,ngrid)+1
 j22=MODULO(j+1,ngrid)+1
 k2=MODULO(k-3,ngrid)+1
 k1=MODULO(k-2,ngrid)+1
 k11=MODULO(k,ngrid)+1
 k22=MODULO(k+1,ngrid)+1
 accx(i,j,k)=-1/dmesh* &
  (pot(i2,j,k)/12-2*pot(i1,j,k)/3+2*pot(i11,j,k)/3-pot(i22,j,k)/12)
 accy(i,j,k)=-1/dmesh* &
  (pot(i,j2,k)/12-2*pot(i,j1,k)/3+2*pot(i,j11,k)/3-pot(i,j22,k)/12)
 accz(i,j,k)=-1/dmesh* &
  (pot(i,j,k2)/12-2*pot(i,j,k1)/3+2*pot(i,j,k11)/3-pot(i,j,k22)/12)
    enddo
   enddo
  enddo  
 end subroutine

 
 subroutine pmgravaccpot(ppos,pacc,ppot)
  real :: ppos(3)
  real, optional :: pacc(3),ppot
  real :: dx(3)
  integer :: i(3),i1(3)

  call getindices(ppos,i,i1,dx)
  if(i1(1).GT.ngrid) i1(1)=1
  if(i1(2).GT.ngrid) i1(2)=1
  if(i1(3).GT.ngrid) i1(3)=1
  if(present(pacc)) then
   pacc(1)=pacc(1)+interpolate(i,i1,dx,accx)
   pacc(2)=pacc(2)+interpolate(i,i1,dx,accy)
   pacc(3)=pacc(3)+interpolate(i,i1,dx,accz) 
  endif
  if(present(ppot)) ppot=ppot+interpolate(i,i1,dx,pot)
 end subroutine
  
 subroutine cic(n,mass,pos)
  integer,intent(in) :: n
  real,intent(in) :: mass(:),pos(:,:)
  integer :: p,i(3),j,k,i1(3),ii,iii,jj,kk
  real :: dx(3),ppos(3),pmass    
  real, allocatable :: localpot(:,:,:) 
  integer :: mythread,totalthread
  integer :: imin(3),imax(3),pmin,pmax,test
!$  integer, external :: omp_get_thread_num,omp_get_num_threads
!$omp parallel &
!$omp default(private) &
!$omp shared(pot,ngrid,pos,mass,n)
  mythread=0
  totalthread=1
!$ mythread=omp_get_thread_num()
!$ totalthread=omp_get_num_threads()
  pmin=(mythread*n)/totalthread+1
  pmax=MIN(((mythread+1)*n)/totalthread,n)
  if(pmax-pmin+1.GT.0) then
! find local extent 
  imin=ngrid
  imax=1
  do p=pmin,pmax
   ppos=pos(p,1:3)
   pmass=mass(p)
   if(pmass.GT.0) then
    call getindices(ppos,i,i1,dx) 
    imin=min(imin,i)   
    imax=max(imax,i1)  
   endif  
  enddo
  allocate(localpot(imin(1):imax(1),imin(2):imax(2),imin(3):imax(3)), &
        stat=test)  
  if(test.NE.0) call Error("local mem. alloc fails")  
  localpot=0.
  do p=pmin,pmax
   ppos=pos(p,1:3)
   pmass=mass(p)
   if(pmass.GT.0) then
    call getindices(ppos,i,i1,dx)
    localpot(i1(1),i1(2),i1(3))=localpot(i1(1),i1(2),i1(3))+pmass*dx(1)*dx(2)*dx(3)
    localpot(i(1), i1(2),i1(3))=localpot(i(1), i1(2),i1(3))+pmass*(1-dx(1))*dx(2)*dx(3)
    localpot(i1(1),i(2), i1(3))=localpot(i1(1),i(2), i1(3))+pmass*dx(1)*(1-dx(2))*dx(3)
    localpot(i(1), i(2), i1(3))=localpot(i(1), i(2), i1(3))+pmass*(1-dx(1))*(1-dx(2))*dx(3)
    localpot(i1(1),i1(2),i(3)) =localpot(i1(1),i1(2),i(3)) +pmass*dx(1)*dx(2)*(1-dx(3))
    localpot(i(1), i1(2),i(3)) =localpot(i(1), i1(2),i(3)) +pmass*(1-dx(1))*dx(2)*(1-dx(3))
    localpot(i1(1),i(2), i(3)) =localpot(i1(1),i(2), i(3)) +pmass*dx(1)*(1-dx(2))*(1-dx(3))
    localpot(i(1), i(2), i(3)) =localpot(i(1), i(2), i(3)) +pmass*(1-dx(1))*(1-dx(2))*(1-dx(3))
   endif
  enddo  
!$omp critical
  do k=imin(3),imax(3)
   do j=imin(2),imax(2)
    do ii=imin(1),imax(1)
    iii=ii;jj=j;kk=k
    if(ii.GT.ngrid) iii=1
    if(j.GT.ngrid) jj=1
    if(k.GT.ngrid) kk=1
    pot(iii,jj,kk)=pot(iii,jj,kk)+localpot(ii,j,k)
    enddo
   enddo
  enddo  
 !$omp end critical  
  deallocate(localpot)
  endif  
!$omp end parallel  
 end subroutine 

 subroutine getindices(ppos,i,i1,dx)
  real ppos(3),dx(3)
  integer :: i(3),i1(3),k
   do k=1,3
    dx(k)=(ppos(k)-boxmin(k))/dmesh
    if(dx(k).lt.0.OR.dx(k).GT.ngrid) call Error("part. outside box")   
    i(k)=FLOOR(dx(k))+1
    i1(k)=i(k)+1
    dx(k)=dx(k)-(i(k)-1)
   enddo   
 end subroutine

 subroutine spread(mass,i,i1,dx,arr)
  real :: dx(3),arr(:,:,:),mass
  integer :: i(3),i1(3)
  arr(i1(1),i1(2),i1(3))=arr(i1(1),i1(2),i1(3))+mass*dx(1)*dx(2)*dx(3)
  arr(i(1), i1(2),i1(3))=arr(i(1), i1(2),i1(3))+mass*(1-dx(1))*dx(2)*dx(3)
  arr(i1(1),i(2), i1(3))=arr(i1(1),i(2), i1(3))+mass*dx(1)*(1-dx(2))*dx(3)
  arr(i(1), i(2), i1(3))=arr(i(1), i(2), i1(3))+mass*(1-dx(1))*(1-dx(2))*dx(3)
  arr(i1(1),i1(2),i(3)) =arr(i1(1),i1(2),i(3)) +mass*dx(1)*dx(2)*(1-dx(3))
  arr(i(1), i1(2),i(3)) =arr(i(1), i1(2),i(3)) +mass*(1-dx(1))*dx(2)*(1-dx(3))
  arr(i1(1),i(2), i(3)) =arr(i1(1),i(2), i(3)) +mass*dx(1)*(1-dx(2))*(1-dx(3))
  arr(i(1), i(2), i(3)) =arr(i(1), i(2), i(3)) +mass*(1-dx(1))*(1-dx(2))*(1-dx(3))
 end subroutine

 function interpolate(i,i1,dx,arr) result(y)
  real :: dx(3),arr(:,:,:),y
  integer :: i(3),i1(3)
 y= arr(i1(1),i1(2),i1(3))*dx(1)*dx(2)*dx(3)+ &
    arr(i(1),i1(2),i1(3))  *(1-dx(1))*dx(2)*dx(3)+ &
    arr(i1(1),i(2),i1(3))  *dx(1)*(1-dx(2))*dx(3)+ &
    arr(i(1),i(2),i1(3)) *(1-dx(1))*(1-dx(2))*dx(3)+ &
    arr(i1(1),i1(2),i(3))  *dx(1)*dx(2)*(1-dx(3))+ &
    arr(i(1),i1(2),i(3)) *(1-dx(1))*dx(2)*(1-dx(3))+ &
    arr(i1(1),i(2),i(3)) *dx(1)*(1-dx(2))*(1-dx(3))+ &
    arr(i(1),i(2),i(3))*(1-dx(1))*(1-dx(2))*(1-dx(3))
 end function

 function ficorrect(r) result(y)
  real :: y,r
  y=erfc(r/2)
 end function

 function forcecorrect(r) result(y)
  real :: y,r
  y=erfc(r/2)+r/sqrt(PI)*exp(-r**2/4)
 end function
 
 real function dmeshQ()
  dmeshQ=dmesh
 end function

 real function softQ()
  softQ=rsoft*dmesh
 end function

 real function rcutQ()
  rcutQ=rcut*rsoft*dmesh
 end function

  
  subroutine Error(string,i)
  character(*) :: string
  integer, optional :: i

  print*,' pmgrav Error detected:'

  if(present(i)) then
   print*,string,i
  else
   print*,string
  endif

  stop
 end subroutine

end module 



subroutine testpm
 use pmgravMod
 real mass(1),pos(1,3),bmin(3),bmax(3),ppos(3),pot,pacc(3)
 integer i
 mass(1)=1
 pos=0.
 bmin=-1
 bmax=1

 call pmgravInit(256,bmin,bmax)
 call pmgrav(1,mass,pos)

ppos=0.
 do i=0,1000
  ppos=i/1000.
  pot=0.
  pacc=0.
  call pmgravaccpot(ppos,ppot=pot,pacc=pacc)
  write(*,*) sqrt(sum((ppos)**2)),sqrt(sum(pacc**2))
 enddo

end
