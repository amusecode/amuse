integer function initialize_code()
!  use amuse_helpers
    implicit none
!  character(len=10), parameter :: default="reflective"
!  ret=amuse_set_boundary(default,default,default,default,default,default)
    initialize_code = 0
!  ret = amuse_initialize_code()
end function

function cleanup_code() result(ret)
    implicit none
    integer :: ret
    ret = 0
end function  

function commit_parameters() result(ret)
    implicit none
    integer :: ret
    ret = 0
end function  

function recommit_parameters() result(ret)
    implicit none
    integer :: ret
    ret = 0
end function  

function setup_mesh(mx,my,mz,xlen,ylen,zlen) result(ret)
    implicit none
    integer :: ret,mx,my,mz
    real*8 :: xlen,ylen,zlen
    ret = 0
end function


function get_mesh_size(mx,my,mz) result(ret)
    implicit none
    integer :: ret,mx,my,mz
    ret = 0
end function
  
function initialize_grid() result(ret)
    implicit none
    integer :: ret
    real*8 :: t0
    ret = 0
end function  

function get_boundary_index_range_inclusive(index_of_boundary, minx, maxx, miny, maxy, minz, maxz) result(ret)
    implicit none
    integer, intent(in) :: index_of_boundary
    integer, intent(out) :: minx, maxx, miny, maxy, minz, maxz
    integer :: ret
    ret = 0
end function


function set_boundary_innerxstate(rho_in,rhvx_in,rhvy_in,rhvz_in,en_in) result(ret)
    implicit none
    integer :: ret
    real*8 :: rho_in,rhvx_in,rhvy_in,rhvz_in,en_in    
    ret = 0
end function  

function set_boundary_innerystate(rho_in,rhvx_in,rhvy_in,rhvz_in,en_in) result(ret)
    implicit none
    integer :: ret
    real*8 :: rho_in,rhvx_in,rhvy_in,rhvz_in,en_in    
    ret = 0
end function  

function set_boundary_innerzstate(rho_in,rhvx_in,rhvy_in,rhvz_in,en_in) result(ret)
    implicit none
    integer :: ret
    real*8 :: rho_in,rhvx_in,rhvy_in,rhvz_in,en_in    
    ret = 0
end function  

function set_boundary_outerxstate(rho_in,rhvx_in,rhvy_in,rhvz_in,en_in) result(ret)
    implicit none
    integer :: ret
    real*8 :: rho_in,rhvx_in,rhvy_in,rhvz_in,en_in    
    ret = 0
end function  

function set_boundary_outerystate(rho_in,rhvx_in,rhvy_in,rhvz_in,en_in) result(ret)
    implicit none
    integer :: ret
    real*8 :: rho_in,rhvx_in,rhvy_in,rhvz_in,en_in    
    ret = 0
end function  

function set_boundary_outerzstate(rho_in,rhvx_in,rhvy_in,rhvz_in,en_in) result(ret)
    implicit none
    integer :: ret
    real*8 :: rho_in,rhvx_in,rhvy_in,rhvz_in,en_in    
    ret = 0
end function  


function set_boundary(lowx,highx,lowy,highy,lowz,highz) result(ret)
    implicit none
    integer :: ret
    character(len=10) :: lowx,highx,lowy,highy,lowz,highz
    ret = 0
end function  

function set_grid_state(i,j,k,rho_in,rhvx_in,rhvy_in,rhvz_in,en_in,n) result(ret)
    implicit none
    integer :: n
    integer :: ret,ii,i(n),j(n),k(n)
    real*8 :: rho_in(n),rhvx_in(n),rhvy_in(n),rhvz_in(n),en_in(n)
    integer,allocatable :: retsum(:)
    ret = 0
end function

function get_grid_state(i,j,k,rho_out,rhvx_out,rhvy_out,rhvz_out,en_out,n) result(ret)
    implicit none
    integer :: n
    integer :: ret,ii,i(n),j(n),k(n)
    real*8 :: rho_out(n),rhvx_out(n),rhvy_out(n),rhvz_out(n),en_out(n)
    integer,allocatable :: retsum(:)
    ret = 0
end function

function get_grid_momentum_density(i,j,k,rhvx_out,rhvy_out,rhvz_out,n) result(ret)
    implicit none
    integer :: n
    integer :: ret,ii,i(n),j(n),k(n)
    real*8 :: rhvx_out(n),rhvy_out(n),rhvz_out(n)
    integer,allocatable :: retsum(:)
    ret = 0
end function

function get_grid_energy_density(i,j,k,en_out,n) result(ret)
    implicit none
    integer :: n
    integer :: ret,ii,i(n),j(n),k(n)
    real*8 :: en_out(n)
    integer,allocatable :: retsum(:)
    ret = 0
end function


function get_grid_density(i,j,k,rho_out,n) result(ret)
    implicit none
    integer :: n
    integer :: ret,ii,i(n),j(n),k(n)
    real*8 :: rho_out(n),rhvx_out(n),rhvy_out(n),rhvz_out(n),en_out(n)
    integer,allocatable :: retsum(:)
    ret = 0
end function

function set_gravity_field(i,j,k,fx,fy,fz) result(ret)
    implicit none
    integer :: ret,i,j,k
    real*8 :: fx,fy,fz
    real*8 :: force(3)
    ret = 0
end function

function get_gravity_field(i,j,k,fx,fy,fz) result(ret)
    implicit none
    integer :: ret,i,j,k
    real*8 :: fx,fy,fz
    real*8 :: force(3)
    ret = 0
end function

function get_position_of_index(i,j,k,xout,yout,zout,n) result(ret)
    implicit none
    integer :: n
    integer :: ret,ii,i(n),j(n),k(n)
    real*8 :: xout(n),yout(n),zout(n)
    integer,allocatable :: retsum(:)
    ret = 0
end function

function get_index_of_position(xin,yin,zin,i,j,k,n) result(ret)
    implicit none
    integer :: n
    integer :: ret,ii,i(n),j(n),k(n)
    real*8 :: xin(n),yin(n),zin(n)
    integer,allocatable :: retsum(:)
    ret = 0
end function

function evolve_model(tend) result(ret)
    implicit none
    integer :: ret
    real*8 :: tend
    ret = 0
end function

function get_time(tnow) result(ret)
    implicit none
    integer :: ret
    real*8 :: tnow
    ret = 0
end function

function get_boundary_state(i,j,k,index_of_boundary,rho_out,rhvx_out,rhvy_out,rhvz_out,en_out,n) result(ret)
    implicit none
    integer :: n, error
    integer :: ret,ii, i0, j0, k0
    integer, intent(in) :: i(n),j(n),k(n), index_of_boundary(n)
    real*8, intent(out) :: rho_out(n),rhvx_out(n),rhvy_out(n),rhvz_out(n),en_out(n)
    real*8, pointer, dimension(:,:,:,:) :: boundary_grid
    integer,allocatable :: retsum(:)
    ret = 0
end function

function set_boundary_state(i,j,k,rho_in,rhvx_in,rhvy_in,rhvz_in,en_in,index_of_boundary,n) result(ret)
    implicit none
    integer :: n, error
    integer :: ret,ii, i0, j0, k0
    integer,intent(in) :: i(n),j(n),k(n), index_of_boundary(n)
    real*8, intent(in) :: rho_in(n),rhvx_in(n),rhvy_in(n),rhvz_in(n),en_in(n)
    integer,allocatable :: retsum(:)
    real*8, pointer, dimension(:,:,:,:) :: boundary_grid
    ret = 0
end function

function get_boundary_position_of_index(i,j,k,index_of_boundary,xout,yout,zout,n) result(ret)
    implicit none
    integer, intent(in) :: n
    integer :: ret,ii, i0, j0, k0
    integer, intent(in) :: i(n),j(n),k(n), index_of_boundary(n)
    real*8, intent(out) :: xout(n),yout(n),zout(n)
    ret = 0
end function

function set_parallel_decomposition(nx, ny, nz) result(ret)
    implicit none
    integer :: ret, ntotal
    integer, intent(in) :: nx, ny, nz
    ret = 0
end function set_parallel_decomposition

function get_parallel_decomposition(nx, ny, nz) result(ret)
    implicit none
    integer :: ret, ntotal
    integer, intent(out) :: nx, ny, nz
end function get_parallel_decomposition

function get_timestep(outputvalue) result(ret)
    implicit none
    integer :: ret
    integer, intent(out) :: outputvalue
    ret = 0
end function get_timestep

function set_timestep(inputvalue) result(ret)
    implicit none
    integer :: ret
    integer, intent(in) :: inputvalue
    ret = 0
end function set_timestep

function set_gamma(inputvalue) result(ret)
    implicit none
    integer :: ret
    double precision :: inputvalue
    ret = 0
end function

function get_gamma(outputvalue) result(ret)
    implicit none
    integer :: ret
    double precision :: outputvalue
    ret = 0
end function

function get_hydro_state_at_point(x1, x2, x3, vx, vy, vz, &
     rho_out, rhovx_out, rhovy_out, rhovz_out, rhoen_out, npoints) result(ret)
    implicit none
    integer :: ret, ii, retsum
    integer, intent(in) :: npoints
    real*8, intent(in) :: x1(npoints), x2(npoints), x3(npoints)
    real*8, intent(in) :: vx(npoints), vy(npoints), vz(npoints)
    real*8, intent(out):: rho_out(npoints), rhovx_out(npoints), rhovy_out(npoints)
    real*8, intent(out):: rhovz_out(npoints), rhoen_out(npoints)
    ret = 0
end function
