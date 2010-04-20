function setup_module() result(ret)
  use amuse_helpers
  integer :: ret
  
  ret=amuse_init()
end function

function cleanup_module() result(ret)
  use amuse_helpers
  integer :: ret
  
  ret=free_grid()
end function  
  
function setup_mesh(mx,my,mz,xlen,ylen,zlen) result(ret)
  use amuse_helpers
  integer :: ret,mx,my,mz
  real*8 :: xlen,ylen,zlen
  character(len=10), parameter :: default="reflective"
  
  ret=amuse_init_mesh(mx,my,mz)
  if(ret.NE.0) return
  ret=amuse_init_coords(xlen,ylen,zlen)
  if(ret.NE.0) return
  ret=amuse_init_hydro()  
  if(ret.NE.0) return
  ret=amuse_set_boundary(default,default,default,default,default,default)
  
end function
  
function initialize_grid(t0) result(ret)
  use amuse_helpers
  integer :: ret
  real*8 :: t0
  
  ret=amuse_initialize_grid(t0)
  
end function  

function set_boundary(lowx,highx,lowy,highy,lowz,highz) result(ret)
  use amuse_helpers
  integer :: ret
  character(len=10) :: lowx,highx,lowy,highy,lowz,highz

  ret=amuse_set_boundary(lowx,highx,lowy,highy,lowz,highz)
  
end function  

! this is somewhat inefficient: blocksize=1
function fill_grid_state(i,j,k,rho_in,rhvx_in,rhvy_in,rhvz_in,en_in) result(ret)
  use amuse_helpers
  integer :: ret,i,j,k
  real*8 :: rho_in,rhvx_in,rhvy_in,rhvz_in,en_in
  real*8 :: lstate(neq)
#ifdef MPI
  integer retsum,ierr
#endif
 
  lstate(RHO)=rho_in
  lstate(RHVX)=rhvx_in
  lstate(RHVY)=rhvy_in
  lstate(RHVZ)=rhvz_in
  lstate(EN)=en_in
  ret=fill_grid(i,j,k,lstate)
#ifdef MPI
  call MPI_REDUCE(ret,retsum,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
  if(retsum.NE.1) ret=-1
#endif
  if(ret.NE.1) then
    ret=-1
    return
  endif  
  ret=0 
end function

function get_grid_state(i,j,k,rho_out,rhvx_out,rhvy_out,rhvz_out,en_out) result(ret)
  use amuse_helpers
  integer :: ret,i,j,k
  real*8 :: rho_out,rhvx_out,rhvy_out,rhvz_out,en_out
  real*8 :: lstate(neq)
#ifdef MPI
  integer retsum,ierr
  real*8 :: tmp(neq)
#endif
  
  ret=retrieve_grid(i,j,k,lstate)
#ifdef MPI
  call MPI_REDUCE(ret,retsum,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
  if(retsum.NE.1) ret=-1
#endif
  if(ret.NE.1) then
    return
  endif  
#ifdef MPI
  tmp=lstate
  call MPI_REDUCE(tmp,lstate,neq,MPI_DOUBLE_PRECISION,MPI_SUM,0.,MPI_COMM_NEW,ierr)
#endif
  rho_out=lstate(RHO)
  rhvx_out=lstate(RHVX)
  rhvy_out=lstate(RHVY)
  rhvz_out=lstate(RHVZ)
  en_out=lstate(EN)
  ret=0
end function

function set_gravity_field(i,j,k,fx,fy,fz) result(ret)
  use amuse_helpers
  integer :: ret,i,j,k
  real*8 :: fx,fy,fz
  real*8 :: force(3)
#ifdef MPI
  integer retsum,ierr
#endif
 
  force(1)=fx
  force(2)=fy
  force(3)=fz
  ret=fill_gforce(i,j,k,force)
#ifdef MPI
  call MPI_REDUCE(ret,retsum,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
  if(retsum.NE.1) ret=-1
#endif
  if(ret.NE.1) then
    ret=-1
    return
  endif  
  ret=0 
end function

function get_gravity_field(i,j,k,fx,fy,fz) result(ret)
  use amuse_helpers
  integer :: ret,i,j,k
  real*8 :: fx,fy,fz
  real*8 :: force(3)
#ifdef MPI
  integer retsum,ierr
  real*8 :: tmp(3)
#endif
  
  ret=retrieve_gforce(i,j,k,force)
#ifdef MPI
  call MPI_REDUCE(ret,retsum,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
  if(retsum.NE.1) ret=-1
#endif
  if(ret.NE.1) then
    return
  endif  
#ifdef MPI
  tmp=force
  call MPI_REDUCE(tmp,force,3,MPI_DOUBLE_PRECISION,MPI_SUM,0.,MPI_COMM_NEW,ierr)
#endif
  fx=force(1)
  fy=force(2)
  fz=force(3)
  ret=0
end function




function get_position_of_index(i,j,k,xout,yout,zout) result(ret)
  use amuse_helpers
  integer :: ret,i,j,k
  real*8 :: xout,yout,zout
#ifdef MPI
  integer retsum,ierr
#endif

  ret=amuse_get_pos_of_index(i,j,k,xout,yout,zout)
#ifdef MPI
  call MPI_REDUCE(ret,retsum,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
  if(retsum.NE.1) ret=-1
#endif
end function

function get_index_of_position(xin,yin,zin,i,j,k) result(ret)
  use amuse_helpers
  integer :: ret,i,j,k
  real*8 :: xin,yin,zin

  ret=amuse_get_index_of_pos(xin,yin,zin,i,j,k)
end function

function evolve(tend) result(ret)
  use amuse_helpers
  integer :: ret
  real*8 :: tend

  ret=amuse_evolve(tend)
  state=>stold
end function

function get_time(tnow) result(ret)
  use amuse_helpers
  integer :: ret
  real*8 :: tnow

  ret=amuse_get_time(tnow)
end function

