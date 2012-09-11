function initialize_code() result(ret)
  use amuse_helpers
  use my_mpi, only: mpi_basic
  implicit none
  include "stopcond.inc"
  integer :: set_support_for_condition
  integer :: ret
  character(len=10), parameter :: default="reflective"
  ret=amuse_set_boundary(default,default,default,default,default,default)
  
  ret = set_support_for_condition(TIMEOUT_DETECTION)
  ret = set_support_for_condition(NUMBER_OF_STEPS_DETECTION)
  call mpi_basic()
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

function get_boundary_index_range_inclusive(index_of_boundary, minx, maxx, miny, maxy, minz, maxz) result(ret)
    use sizes, only: mbc
    use mesh, only: sx,ex,sy,ey,sz,ez,meshx,meshy,meshz
    use problem
    implicit none
    integer, intent(in) :: index_of_boundary
    integer, intent(out) :: minx, maxx, miny, maxy, minz, maxz
    integer :: ret
    minx = 1
    miny = 1
    minz = 1
    
    maxx = 1
    maxy = 1 
    maxz = 1
    
    if (index_of_boundary < 1 .OR. index_of_boundary > 6) then
        ret = -1
        return
    end if
    
    if (domainboundaryconditions((index_of_boundary + 1)/ 2,MODULO(index_of_boundary + 1, 2) + 1) .EQ. PROBLEM_DEF)  then
        select case (index_of_boundary)
            case (1)
                maxx = mbc
                maxy = meshy
                maxz = meshz
            case (2)
                maxx = mbc
                maxy = meshy
                maxz = meshz
            case (3)
                maxx = mbc * 2 + meshx
                maxy = mbc
                maxz = meshz
            case (4)
                maxx = mbc * 2 + meshx
                maxy = mbc
                maxz = meshz
            case (5)
                maxx = mbc * 2 + meshx
                maxy = mbc * 2 + meshy
                maxz = mbc
            case (6)
                maxx = mbc * 2 + meshx
                maxy = mbc * 2 + meshy
                maxz = mbc
        end select
    end if
    ret = 0
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

function set_grid_state(i,j,k,rho_in,rhvx_in,rhvy_in,rhvz_in,en_in,n) result(ret)
  use amuse_helpers
  integer :: n
  integer :: ret,ii,i(n),j(n),k(n)
  real*8 :: rho_in(n),rhvx_in(n),rhvy_in(n),rhvz_in(n),en_in(n)
  real*8 :: lstate(neq)
  integer,allocatable :: retsum(:)
#ifdef MPI
  integer ierr
#endif
  allocate(retsum(n))
 
  do ii=1,n
    lstate(RHO)=rho_in(ii)
    lstate(RHVX)=rhvx_in(ii)
    lstate(RHVY)=rhvy_in(ii)
    lstate(RHVZ)=rhvz_in(ii)
    lstate(EN)=en_in(ii)
    retsum(ii)=fill_grid(i(ii),j(ii),k(ii),lstate)
  enddo
#ifdef MPI
  if(rank.NE.0) then
    call MPI_REDUCE(retsum,0,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
  else
    call MPI_REDUCE(MPI_IN_PLACE,retsum,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
  endif
#endif  
  ret=0
  if(any(retsum.NE.1)) ret=-1
  deallocate(retsum)
end function

function get_grid_state(i,j,k,rho_out,rhvx_out,rhvy_out,rhvz_out,en_out,n) result(ret)
  use amuse_helpers
  integer :: n
  integer :: ret,ii,i(n),j(n),k(n)
  real*8 :: rho_out(n),rhvx_out(n),rhvy_out(n),rhvz_out(n),en_out(n)
  real*8 :: lstate(neq)
  integer,allocatable :: retsum(:)
#ifdef MPI
  integer ierr
#endif
  allocate(retsum(n))
  
  do ii=1,n
    retsum(ii)=retrieve_grid(i(ii),j(ii),k(ii),lstate)
    rho_out(ii)=lstate(RHO)
    rhvx_out(ii)=lstate(RHVX)
    rhvy_out(ii)=lstate(RHVY)
    rhvz_out(ii)=lstate(RHVZ)
    en_out(ii)=lstate(EN)
  enddo

#ifdef MPI
  if(rank.NE.0) then
    call MPI_REDUCE(retsum,0,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(rho_out,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(rhvx_out,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(rhvy_out,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(rhvz_out,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(en_out,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
  else
    call MPI_REDUCE(MPI_IN_PLACE,retsum,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,rho_out,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,rhvx_out,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,rhvy_out,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,rhvz_out,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,en_out,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
  endif
#endif  
  ret=0
  if(any(retsum.NE.1)) ret=-1
  deallocate(retsum)
end function



function get_grid_momentum_density(i,j,k,rhvx_out,rhvy_out,rhvz_out,n) result(ret)
  use amuse_helpers
  integer :: n
  integer :: ret,ii,i(n),j(n),k(n)
  real*8 :: rhvx_out(n),rhvy_out(n),rhvz_out(n)
  real*8 :: lstate(neq)
  integer,allocatable :: retsum(:)
#ifdef MPI
  integer ierr
#endif
  allocate(retsum(n))
  
  do ii=1,n
    retsum(ii)=retrieve_grid(i(ii),j(ii),k(ii),lstate)
    rhvx_out(ii)=lstate(RHVX)
    rhvy_out(ii)=lstate(RHVY)
    rhvz_out(ii)=lstate(RHVZ)
  enddo

#ifdef MPI
  if(rank.NE.0) then
    call MPI_REDUCE(retsum,0,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(rhvx_out,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(rhvy_out,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(rhvz_out,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
  else
    call MPI_REDUCE(MPI_IN_PLACE,retsum,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,rhvx_out,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,rhvy_out,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,rhvz_out,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
  endif
#endif  
  ret=0
  if(any(retsum.NE.1)) ret=-1
  deallocate(retsum)
end function


function get_grid_energy_density(i,j,k,en_out,n) result(ret)
  use amuse_helpers
  integer :: n
  integer :: ret,ii,i(n),j(n),k(n)
  real*8 :: en_out(n)
  real*8 :: lstate(neq)
  integer,allocatable :: retsum(:)
#ifdef MPI
  integer ierr
#endif
  allocate(retsum(n))
  
  do ii=1,n
    retsum(ii)=retrieve_grid(i(ii),j(ii),k(ii),lstate)
    en_out(ii)=lstate(EN)
  enddo

#ifdef MPI
  if(rank.NE.0) then
    call MPI_REDUCE(retsum,0,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(en_out,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
  else
    call MPI_REDUCE(MPI_IN_PLACE,retsum,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,en_out,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
  endif
#endif  
  ret=0
  if(any(retsum.NE.1)) ret=-1
  deallocate(retsum)
end function


function get_grid_density(i,j,k,rho_out,n) result(ret)
  use amuse_helpers
  integer :: n
  integer :: ret,ii,i(n),j(n),k(n)
  real*8 :: rho_out(n),rhvx_out(n),rhvy_out(n),rhvz_out(n),en_out(n)
  real*8 :: lstate(neq)
  integer,allocatable :: retsum(:)
#ifdef MPI
  integer ierr
#endif
  allocate(retsum(n))
  
  do ii=1,n
    retsum(ii)=retrieve_grid(i(ii),j(ii),k(ii),lstate)
    rho_out(ii)=lstate(RHO)
  enddo

#ifdef MPI
  if(rank.NE.0) then
    call MPI_REDUCE(retsum,0,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(rho_out,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
  else
    call MPI_REDUCE(MPI_IN_PLACE,retsum,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,rho_out,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
  endif
#endif  
  ret=0
  if(any(retsum.NE.1)) ret=-1
  deallocate(retsum)
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


function get_position_of_index(i,j,k,xout,yout,zout,n) result(ret)
  use amuse_helpers
  integer :: n
  integer :: ret,ii,i(n),j(n),k(n)
  real*8 :: xout(n),yout(n),zout(n)
  integer,allocatable :: retsum(:)
#ifdef MPI
  integer ierr
#endif
  allocate(retsum(n))
  
  do ii=1,n
    retsum(ii)=amuse_get_pos_of_index(i(ii),j(ii),k(ii),xout(ii),yout(ii),zout(ii))
  enddo

#ifdef MPI
  if(rank.NE.0) then
    call MPI_REDUCE(retsum,0,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(xout,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(yout,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(zout,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
  else
    call MPI_REDUCE(MPI_IN_PLACE,retsum,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,xout,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,yout,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,zout,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
  endif
#endif  
  ret=0
  if(any(retsum.NE.1)) ret=-1
  deallocate(retsum)
end function

function get_index_of_position(xin,yin,zin,i,j,k,n) result(ret)
  use amuse_helpers
  integer :: n
  integer :: ret,ii,i(n),j(n),k(n)
  real*8 :: xin(n),yin(n),zin(n)
  integer,allocatable :: retsum(:)
#ifdef MPI
  integer ierr
#endif
  allocate(retsum(n))
  
  do ii=1,n
    retsum(ii)=amuse_get_index_of_pos(xin(ii),yin(ii),zin(ii),i(ii),j(ii),k(ii))
  enddo
  
  ret=0
  if(any(retsum.NE.1)) ret=-1
  deallocate(retsum)
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


function get_boundary_state(i,j,k,index_of_boundary,rho_out,rhvx_out,rhvy_out,rhvz_out,en_out,n) result(ret)
  use amuse_helpers
  implicit none
  
  integer :: n, error
  integer :: ret,ii, i0, j0, k0
  integer, intent(in) :: i(n),j(n),k(n), index_of_boundary(n)
  real*8,  intent(out) :: rho_out(n),rhvx_out(n),rhvy_out(n),rhvz_out(n),en_out(n)
  real*8 :: lstate(neq)
  real(kind=dp),pointer,dimension(:,:,:,:) :: boundary_grid
  integer,allocatable :: retsum(:)
#ifdef MPI
  integer ierr
#endif
  allocate(retsum(n))
  
  do ii=1,n
    i0 = i(ii)
    j0 = j(ii)
    k0 = k(ii)
    retsum(ii) = 0
    if (is_on_boundary_grid(index_of_boundary(ii), i0, j0, k0, .FALSE.)) then
        retsum(ii) = 1
        boundary_grid => get_boundary_grid_pointer(index_of_boundary(ii))
        error = on_boundary_grid(index_of_boundary(ii), i0, j0, k0)
    
        rho_out(ii)=boundary_grid(i0, j0, k0, RHO)
        rhvx_out(ii)=boundary_grid(i0, j0, k0, RHVX)
        rhvy_out(ii)=boundary_grid(i0, j0, k0, RHVY)
        rhvz_out(ii)=boundary_grid(i0, j0, k0, RHVZ)
        en_out(ii)=boundary_grid(i0, j0, k0, EN)
    else
        rho_out(ii)=0
        rhvx_out(ii)=0
        rhvy_out(ii)=0
        rhvz_out(ii)=0
        en_out(ii)=0
    end if
  enddo

#ifdef MPI
  if(rank.NE.0) then
    call MPI_REDUCE(retsum,0,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(rho_out,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(rhvx_out,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(rhvy_out,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(rhvz_out,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(en_out,0,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
  else
    call MPI_REDUCE(MPI_IN_PLACE,retsum,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,rho_out,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,rhvx_out,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,rhvy_out,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,rhvz_out,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
    call MPI_REDUCE(MPI_IN_PLACE,en_out,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_NEW,ierr)
  endif
#endif  
  ret=0
  if(any(retsum.NE.1)) ret=-1
  deallocate(retsum)
end function

    

function set_boundary_state(i,j,k,rho_in,rhvx_in,rhvy_in,rhvz_in,en_in,index_of_boundary,n) result(ret)
  use amuse_helpers
  implicit none
  
  integer :: n, error
  integer :: ret,ii, i0, j0, k0
  integer,intent(in) :: i(n),j(n),k(n), index_of_boundary(n)
  real*8, intent(in) :: rho_in(n),rhvx_in(n),rhvy_in(n),rhvz_in(n),en_in(n)
  integer,allocatable :: retsum(:)
  real(kind=dp),pointer,dimension(:,:,:,:) :: boundary_grid
#ifdef MPI
  integer ierr
#endif
  allocate(retsum(n))
 
  do ii=1,n
    retsum(ii) = 0
    i0 = i(ii)
    j0 = j(ii)
    k0 = k(ii)
    if (is_on_boundary_grid(index_of_boundary(ii), i0, j0, k0, .TRUE.)) then
        if (is_on_boundary_grid(index_of_boundary(ii), i0, j0, k0, .FALSE.)) then
            retsum(ii) = 1
        end if
        error = on_boundary_grid(index_of_boundary(ii), i0, j0, k0)
        boundary_grid => get_boundary_grid_pointer(index_of_boundary(ii))
        boundary_grid(i0, j0, k0, RHO)=rho_in(ii)
        boundary_grid(i0, j0, k0, RHVX)=rhvx_in(ii)
        boundary_grid(i0, j0, k0, RHVY)=rhvy_in(ii)
        boundary_grid(i0, j0, k0, RHVZ)=rhvz_in(ii)
        boundary_grid(i0, j0, k0, EN)=en_in(ii)
    end if
  enddo
#ifdef MPI
  if(rank.NE.0) then
    call MPI_REDUCE(retsum,0,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
  else
    call MPI_REDUCE(MPI_IN_PLACE,retsum,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_NEW,ierr)
  endif
#endif  
  ret=0
  if(any(retsum.NE.1)) ret=-1
  deallocate(retsum)
end function




function get_boundary_position_of_index(i,j,k,index_of_boundary,xout,yout,zout,n) result(ret)
  use amuse_helpers
  implicit none
  
  integer, intent(in) :: n
  integer :: ret,ii, i0, j0, k0
  integer, intent(in) :: i(n),j(n),k(n), index_of_boundary(n)
  real*8, intent(out) :: xout(n),yout(n),zout(n)
  
   do ii=1,n
    i0 = i(ii)
    j0 = j(ii)
    k0 = k(ii)
    select case (index_of_boundary(ii))
        case(1)
            xout(ii) = dx*(0.5_dp-real(i0,dp))
            yout(ii) = dy*(0.5_dp+real(j0-1,dp))
            zout(ii) = dz*(0.5_dp+real(k0-1,dp))
        case(2)
            xout(ii) = dx*(0.5_dp+real(ex + i0 - 1,dp))
            yout(ii) = dy*(0.5_dp+real(j0-1,dp))
            zout(ii) = dz*(0.5_dp+real(k0-1,dp))
        case(3)
            xout(ii) = dx*(0.5_dp+real(i0-mbc-1,dp))
            yout(ii) = dy*(0.5_dp-real(j0,dp))
            zout(ii) = dz*(0.5_dp+real(k0-1,dp))
        case(4)
            xout(ii) = dx*(0.5_dp+real(i0-mbc-1,dp))
            yout(ii) = dy*(0.5_dp+real(ey + j0 - 1,dp))
            zout(ii) = dz*(0.5_dp+real(k0-1,dp))
        case(5)
            xout(ii) = dx*(0.5_dp+real(i0-mbc-1,dp))
            yout(ii) = dy*(0.5_dp+real(j0-mbc-1,dp))
            zout(ii) = dz*(0.5_dp-real(k0,dp))
        case(6)
            xout(ii) = dx*(0.5_dp+real(i0-mbc-1,dp))
            yout(ii) = dy*(0.5_dp+real(j0-mbc-1,dp))
            zout(ii) = dz*(0.5_dp+real(ez + k0 - 1,dp))
    end select
  enddo

  ret=0
end function


function set_parallel_decomposition(nx, ny, nz) result(ret)
  use amuse_helpers
  implicit none
  integer :: ret, ntotal
  integer, intent(in) :: nx, ny, nz
  integer :: localdims(NPDIM)
  localdims = 0
  ntotal = 1
  ret = 0
  if (nx.GT.0) then
    ntotal = ntotal * nx
  end if
  if (ny.GT.0) then
    ntotal = ntotal * ny
  end if
  if (nz.GT.0) then
    ntotal = ntotal * nz
  end if
  if(ntotal .GT. npr) then
    call set_dims(dims)
    ret = -1
    return
  end if
  localdims(1) = nx
  localdims(2) = ny
  localdims(3) = nz
  call set_dims(localdims)
end function set_parallel_decomposition

function get_parallel_decomposition(nx, ny, nz) result(ret)
  use amuse_helpers
  implicit none
  integer :: ret, ntotal
  integer, intent(out) :: nx, ny, nz
  integer :: localdims(NPDIM)
  localdims = 0
  ret = 0
  call get_dims(localdims)
  
  nx = localdims(1)
  ny = localdims(2)
  nz = localdims(3)
end function get_parallel_decomposition



function get_timestep(outputvalue) result(ret)
  use amuse_helpers
  implicit none
  integer :: ret
  integer, intent(out) :: outputvalue
  
  outputvalue = dt
end function get_timestep


function set_timestep(inputvalue) result(ret)
  use amuse_helpers
  implicit none
  integer :: ret
  integer, intent(in) :: inputvalue
  
end function set_timestep
