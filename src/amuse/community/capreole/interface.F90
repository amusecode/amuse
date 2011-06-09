function initialize_code() result(ret)
  use amuse_helpers
  integer :: ret
  character(len=10), parameter :: default="reflective"
  ret=amuse_set_boundary(default,default,default,default,default,default)
  
end function

function cleanup_code() result(ret)
  use amuse_helpers
  integer :: ret
  
  ret=free_grid()
end function  

function commit_parameters() result(ret)
  use amuse_helpers
  integer :: ret
  ret = amuse_commit_parameters()
end function  

function recommit_parameters() result(ret)
  use amuse_helpers
  integer :: ret
  
  ret=-1
end function  
  
function setup_mesh(mx,my,mz,xlen,ylen,zlen) result(ret)
  use amuse_helpers
  integer :: ret,mx,my,mz
  real*8 :: xlen,ylen,zlen
  
  ret=amuse_init()
  if(ret.NE.0) return
  ret=amuse_init_mesh(mx,my,mz)
  if(ret.NE.0) return
  ret=amuse_init_coords(xlen,ylen,zlen)
  if(ret.NE.0) return
  
end function


function get_mesh_size(mx,my,mz) result(ret)
  use amuse_helpers
  integer :: ret,mx,my,mz
  
  ret = amuse_get_mesh(mx,my,mz)
  
end function
  
function initialize_grid() result(ret)
  use amuse_helpers
  integer :: ret
  real*8 :: t0
  
  t0  = 0.
  ret = amuse_initialize_grid(t0)
  
end function  

function set_boundary_innerxstate(rho_in,rhvx_in,rhvy_in,rhvz_in,en_in) result(ret)
  use amuse_helpers
  integer :: ret
  real*8 :: tmp(neq),rho_in,rhvx_in,rhvy_in,rhvz_in,en_in    
  tmp(RHO)=rho_in
  tmp(RHVX)=rhvx_in
  tmp(RHVY)=rhvy_in
  tmp(RHVZ)=rhvz_in
  tmp(EN)=en_in
  ret=amuse_set_boundary_innerxstate(tmp)
end function  

function set_boundary_innerystate(rho_in,rhvx_in,rhvy_in,rhvz_in,en_in) result(ret)
  use amuse_helpers
  integer :: ret
  real*8 :: tmp(neq),rho_in,rhvx_in,rhvy_in,rhvz_in,en_in    
  tmp(RHO)=rho_in
  tmp(RHVX)=rhvx_in
  tmp(RHVY)=rhvy_in
  tmp(RHVZ)=rhvz_in
  tmp(EN)=en_in
  ret=amuse_set_boundary_innerystate(tmp)
end function  

function set_boundary_innerzstate(rho_in,rhvx_in,rhvy_in,rhvz_in,en_in) result(ret)
  use amuse_helpers
  integer :: ret
  real*8 :: tmp(neq),rho_in,rhvx_in,rhvy_in,rhvz_in,en_in    
  tmp(RHO)=rho_in
  tmp(RHVX)=rhvx_in
  tmp(RHVY)=rhvy_in
  tmp(RHVZ)=rhvz_in
  tmp(EN)=en_in
  ret=amuse_set_boundary_innerzstate(tmp)
end function  

function set_boundary_outerxstate(rho_in,rhvx_in,rhvy_in,rhvz_in,en_in) result(ret)
  use amuse_helpers
  integer :: ret
  real*8 :: tmp(neq),rho_in,rhvx_in,rhvy_in,rhvz_in,en_in    
  tmp(RHO)=rho_in
  tmp(RHVX)=rhvx_in
  tmp(RHVY)=rhvy_in
  tmp(RHVZ)=rhvz_in
  tmp(EN)=en_in
  ret=amuse_set_boundary_outerxstate(tmp)
end function  

function set_boundary_outerystate(rho_in,rhvx_in,rhvy_in,rhvz_in,en_in) result(ret)
  use amuse_helpers
  integer :: ret
  real*8 :: tmp(neq),rho_in,rhvx_in,rhvy_in,rhvz_in,en_in    
  tmp(RHO)=rho_in
  tmp(RHVX)=rhvx_in
  tmp(RHVY)=rhvy_in
  tmp(RHVZ)=rhvz_in
  tmp(EN)=en_in
  ret=amuse_set_boundary_outerystate(tmp)
end function  

function set_boundary_outerzstate(rho_in,rhvx_in,rhvy_in,rhvz_in,en_in) result(ret)
  use amuse_helpers
  integer :: ret
  real*8 :: tmp(neq),rho_in,rhvx_in,rhvy_in,rhvz_in,en_in    
  tmp(RHO)=rho_in
  tmp(RHVX)=rhvx_in
  tmp(RHVY)=rhvy_in
  tmp(RHVZ)=rhvz_in
  tmp(EN)=en_in
  ret=amuse_set_boundary_outerzstate(tmp)
end function  


function set_boundary(lowx,highx,lowy,highy,lowz,highz) result(ret)
  use amuse_helpers
  integer :: ret
  character(len=10) :: lowx,highx,lowy,highy,lowz,highz

  ret=amuse_set_boundary(lowx,highx,lowy,highy,lowz,highz)
  
end function  

! this is somewhat inefficient: blocksize=1
function set_grid_state(i,j,k,rho_in,rhvx_in,rhvy_in,rhvz_in,en_in) result(ret)
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
  call MPI_ALLREDUCE(ret,retsum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_NEW,ierr)
  ret=1
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
  call MPI_ALLREDUCE(ret,retsum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_NEW,ierr)
  ret=1
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



function get_grid_momentum_density(i,j,k,rhvx_out,rhvy_out,rhvz_out) result(ret)
  use amuse_helpers
  integer :: ret,i,j,k
  real*8 :: rhvx_out,rhvy_out,rhvz_out
  real*8 :: lstate(neq)
#ifdef MPI
  integer retsum,ierr
  real*8 :: tmp(neq)
#endif
  
  ret=retrieve_grid(i,j,k,lstate)
#ifdef MPI
  call MPI_ALLREDUCE(ret,retsum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_NEW,ierr)
  ret=1
  if(retsum.NE.1) ret=-1
#endif
  if(ret.NE.1) then
    return
  endif  
#ifdef MPI
  tmp=lstate
  call MPI_REDUCE(tmp,lstate,neq,MPI_DOUBLE_PRECISION,MPI_SUM,0.,MPI_COMM_NEW,ierr)
#endif
  rhvx_out=lstate(RHVX)
  rhvy_out=lstate(RHVY)
  rhvz_out=lstate(RHVZ)
  ret=0
end function


function get_grid_energy_density(i,j,k,en_out) result(ret)
  use amuse_helpers
  integer :: ret,i,j,k
  real*8 :: en_out
  real*8 :: lstate(neq)
#ifdef MPI
  integer retsum,ierr
  real*8 :: tmp(neq)
#endif
  
  ret=retrieve_grid(i,j,k,lstate)
#ifdef MPI
  call MPI_ALLREDUCE(ret,retsum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_NEW,ierr)
  ret=1
  if(retsum.NE.1) ret=-1
#endif
  if(ret.NE.1) then
    return
  endif  
#ifdef MPI
  tmp=lstate
  call MPI_REDUCE(tmp,lstate,neq,MPI_DOUBLE_PRECISION,MPI_SUM,0.,MPI_COMM_NEW,ierr)
#endif
  en_out=lstate(EN)
  ret=0
end function


function get_grid_density(i,j,k,rho_out) result(ret)
  use amuse_helpers
  integer :: ret,i,j,k
  real*8 :: rho_out
  real*8 :: lstate(neq)
#ifdef MPI
  integer retsum,ierr
  real*8 :: tmp(neq)
#endif
  
  ret=retrieve_grid(i,j,k,lstate)
#ifdef MPI
  call MPI_ALLREDUCE(ret,retsum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_NEW,ierr)
  ret=1
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
  call MPI_ALLREDUCE(ret,retsum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_NEW,ierr)
  ret=1
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
  call MPI_ALLREDUCE(ret,retsum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_NEW,ierr)
  ret=1
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
  real*8 posout(3),tmp(3)
#endif

  ret=amuse_get_pos_of_index(i,j,k,xout,yout,zout)
#ifdef MPI
  call MPI_ALLREDUCE(ret,retsum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_NEW,ierr)
  ret=1
  if(retsum.NE.1) ret=-1
#endif
  if(ret.NE.1) then
    return
  endif  
#ifdef MPI
  tmp(1)=xout
  tmp(2)=yout
  tmp(3)=zout
  call MPI_REDUCE(tmp,posout,3,MPI_DOUBLE_PRECISION,MPI_SUM,0.,MPI_COMM_NEW,ierr)
  xout=posout(1)
  yout=posout(2)
  zout=posout(3)
#endif
  ret=0
end function

function get_index_of_position(xin,yin,zin,i,j,k) result(ret)
  use amuse_helpers
  integer :: ret,i,j,k
  real*8 :: xin,yin,zin

  ret=amuse_get_index_of_pos(xin,yin,zin,i,j,k)
end function

function evolve_model(tend) result(ret)
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

