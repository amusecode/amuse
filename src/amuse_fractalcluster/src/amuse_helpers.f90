module fractalMod

 integer :: nstar=0, iseed=-1, n_generated=0
 real :: fdim=1.6
 
 real, allocatable :: r(:,:), v(:,:)

 
contains

function amuse_get_state(id,x,y,z,vx,vy,vz,n) result(err)
  integer :: n,i,id(n),err
  real :: x(n),y(n),z(n)
  real :: vx(n),vy(n),vz(n)
  err=0
  x=r(1,id(1:n))
  y=r(2,id(1:n))
  z=r(3,id(1:n))
  vx=v(1,id(1:n))
  vy=v(2,id(1:n))
  vz=v(3,id(1:n))    
end function

function amuse_generator() result(err)
  integer :: err, Count
  err=0
  if(nstar.LE.0) then 
    err=-1
    return
  endif  
  if(allocated(r)) deallocate(r,v)
  allocate(r(3,nstar),v(3,nstar))
  if (iseed.LE.0) then
     call system_clock(Count)
     iseed = Count
  endif
  n_generated = nstar
  call makefractal(n_generated,r,v,fdim,iseed)
end function

end module
