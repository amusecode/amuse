function get_state(id,x,y,z,vx,vy,vz,n) result(err)
  use fractalMod
  integer :: n,id(n),err
  real :: x(n),y(n),z(n)
  real :: vx(n),vy(n),vz(n)
  err=amuse_get_state(id,x,y,z,vx,vy,vz,n)
end function

function generator() result(err)
  use fractalMod
  integer :: err
  err=amuse_generator()
end function  

function set_fractal_dimension(fd) result(err)
  use fractalMod
  integer :: err
  real :: fd
  err=0
  fdim=fd
end function

function get_fractal_dimension(fd) result(err)
  use fractalMod
  integer :: err
  real :: fd
  err=0
  fd=fdim
end function

function set_random_seed(i) result(err)
  use fractalMod
  integer i,err
  err=0
  iseed=i
end function

function get_random_seed(i) result(err)
  use fractalMod
  integer i,err
  err=0
  i=iseed
end function

function get_nstar(i) result(err)
  use fractalMod
  integer i,err
  err=0
  i=nstar
end function

function set_nstar(i) result(err)
  use fractalMod
  integer i,err
  err=0
  nstar=i
end function
