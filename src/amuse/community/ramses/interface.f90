module amuse_globals
  integer::ilevel,idim,ivar,info
  real(kind=8)::tt1,tt2
  real(kind=4)::real_mem,real_mem_tot
end module amuse_globals

integer function initialize_code()
    implicit none
    call read_params
    initialize_code = 0
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
    use amuse_globals
!  use amr_commons
!  use hydro_commons
!  use pm_commons
  use poisson_commons
  use cooling_module
#ifdef RT
  use rt_hydro_commons
#endif
    implicit none
    integer :: ret
    real*8 :: t0
#ifndef WITHOUTMPI
  tt1=MPI_WTIME(info)
#endif

  call init_amr                      ! Initialize AMR variables
  call init_time                     ! Initialize time variables
  if(hydro)call init_hydro           ! Initialize hydro variables
#ifdef RT
  if(rt.or.neq_chem) &
       call rt_init_hydro            ! Initialize radiation variables
#endif
  if(poisson)call init_poisson       ! Initialize poisson variables
#ifdef ATON
  if(aton)call init_radiation        ! Initialize radiation variables
#endif
  if(nrestart==0)call init_refine    ! Build initial AMR grid
  if(cooling.and..not.neq_chem) &
       call set_table(dble(aexp))    ! Initialize cooling look up table
  if(pic)call init_part              ! Initialize particle variables
  if(pic)call init_tree              ! Initialize particle tree
  if(nrestart==0)call init_refine_2  ! Build initial AMR grid again

#ifndef WITHOUTMPI
  tt2=MPI_WTIME(info)
  if(myid==1)write(*,*)'Time elapsed since startup:',tt2-tt1
#endif

  if(myid==1)then
     write(*,*)'Initial mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if

  nstep_coarse_old=nstep_coarse
    ret = 0
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')
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
    use amuse_globals
    use amuse_error_code
!  use amr_commons
  use hydro_commons
!  use pm_commons
  use poisson_commons
!  use cooling_module
#ifdef RT
  use rt_hydro_commons
#endif
    implicit none
    integer :: ret
    real*8 :: tend
    
    tout(noutput) = tend
    
  if(myid==1)write(*,*)'Starting time integration' 

  do ! Main time loop

#ifndef WITHOUTMPI
     tt1=MPI_WTIME(info)
#endif

     if(verbose)write(*,*)'Entering amr_step_coarse'

     epot_tot=0.0D0  ! Reset total potential energy
     ekin_tot=0.0D0  ! Reset total kinetic energy
     mass_tot=0.0D0  ! Reset total mass
     eint_tot=0.0D0  ! Reset total internal energy
#ifdef SOLVERmhd
     emag_tot=0.0D0  ! Reset total magnetic energy
#endif

     ! Make new refinements
     if(levelmin.lt.nlevelmax .and..not.static)then
        call refine_coarse
        do ilevel=1,levelmin
           call build_comm(ilevel)
           call make_virtual_fine_int(cpu_map(1),ilevel)
           if(hydro)then
#ifdef SOLVERmhd
              do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
              end do
#else
              end do
#endif
              if(simple_boundary)call make_boundary_hydro(ilevel)
           endif
#ifdef RT
           if(rt)then
              do ivar=1,nrtvar
                 call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
              end do
              if(simple_boundary)call rt_make_boundary_hydro(ilevel)
           endif
#endif
           if(poisson)then
              call make_virtual_fine_dp(phi(1),ilevel)
              do idim=1,ndim
                 call make_virtual_fine_dp(f(1,idim),ilevel)
              end do
           end if
           if(ilevel<levelmin)call refine_fine(ilevel)
        end do
     endif

     ! Call base level
     call amr_step(levelmin,1)

     if(levelmin.lt.nlevelmax .and..not. static)then
        do ilevel=levelmin-1,1,-1
           ! Hydro book-keeping
           if(hydro)then
              call upload_fine(ilevel)
#ifdef SOLVERmhd
              do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
              end do
#else
              end do
#endif
              if(simple_boundary)call make_boundary_hydro(ilevel)
           end if
#ifdef RT
           ! Radiation book-keeping
           if(rt)then
              call rt_upload_fine(ilevel)
              do ivar=1,nrtvar
                 call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
              end do
              if(simple_boundary)call rt_make_boundary_hydro(ilevel)
           end if
#endif
           ! Gravity book-keeping
           if(poisson)then
              call make_virtual_fine_dp(phi(1),ilevel)
              do idim=1,ndim
                 call make_virtual_fine_dp(f(1,idim),ilevel)
              end do
           end if
        end do
        
        ! Build refinement map
        do ilevel=levelmin-1,1,-1
           call flag_fine(ilevel,2)
        end do
        call flag_coarse
     endif

     ! New coarse time-step
     nstep_coarse=nstep_coarse+1

#ifndef WITHOUTMPI
     tt2=MPI_WTIME(info)
     if(mod(nstep_coarse,ncontrol)==0)then
        call getmem(real_mem)
        call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
        if(myid==1)then
           write(*,*)'Time elapsed since last coarse step:',tt2-tt1
           call writemem(real_mem_tot)
        endif
     endif
#endif

     write(*,*) "Error code:", error_code
     if (error_code.eq.1) then
        error_code = 0
        if (t>=tout(noutput)) then
            ret = 0
        else
            ret = -1
        endif
        exit
     endif
  end do
end function

function get_time(tnow) result(ret)
    use amr_commons
    implicit none
    integer :: ret
    real*8 :: tnow
    tnow = t
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
