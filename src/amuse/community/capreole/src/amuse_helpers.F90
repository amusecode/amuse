module amuse_helpers
  use my_mpi
  use evolution
  use mesh
  use grid
  use times
  use hydro
  use scaling, only: SCDENS,SCMOME,SCENER,SCTIME,SCLENG
  use boundary, only: boundaries,REFLECTIVE,REFLECTIVE_SHIFT, &
    OUTFLOW,PROBLEM_DEF,PERIODIC
  use sizes, only: neq,neuler,mbc,nrOfDim,RHO,RHVX,RHVY,RHVZ,EN  
  use problem
  use atomic, only: gamma, gamma1
  
  integer, private, save :: nstep
  
  contains

  function amuse_init() result(ret)
    integer :: ret
!  call setup_clocks ()
    call mpi_setup()
    nstep = 0
    ret=0    
  end function

  function amuse_evolve(tend) result(ret)
    include 'stopcond.inc'
    integer :: ret,nf
    real*8 :: tend
    integer :: is_any_condition_set
    integer :: is_stopping_condition_enabled
    integer :: is_number_of_steps_detection_enabled
    integer :: is_timeout_detection_enabled
    integer :: get_stopping_condition_number_of_steps_parameter 
    integer :: get_stopping_condition_timeout_parameter 
    integer :: clock_init, clock_current, count_rate, count_max
    integer :: max_number_of_steps
    integer :: timeout
    integer :: number_of_steps_innerloop
    integer :: stopping_index
    integer :: next_index_for_stopping_condition
    integer :: set_stopping_condition_info
    integer :: reset_stopping_conditions, error
    real(kind=dp)    :: nexttime   ! timer for output
    integer :: nframe              ! integers for output
    
    error = reset_stopping_conditions()
    error = is_stopping_condition_enabled(NUMBER_OF_STEPS_DETECTION, is_number_of_steps_detection_enabled)
    error = is_stopping_condition_enabled(TIMEOUT_DETECTION, is_timeout_detection_enabled)
    error = get_stopping_condition_number_of_steps_parameter(max_number_of_steps)
    error = get_stopping_condition_timeout_parameter(timeout)
    if(error /= 0) then
        ret = -2
        return
    end if
    call SYSTEM_CLOCK(clock_init, count_rate, count_max)
    number_of_steps_innerloop = 0
    if(tend-time.GT.0) then
      lastframe=0
      nf=1
      do while(nf.GE.lastframe)
        lastframe=lastframe+1
        frametime=tend/lastframe
        nf=nint(time/frametime)
      enddo
      
      ! Set time for next output
      nframe=nint(time/frametime)
      nexttime=real(nframe+1,dp)*frametime
      nframe=nframe+1    ! update output counter
    
      do
        error = evolve_step(nstep, nexttime)
        if (time >= nexttime) then
          state => stold
          nstep = 0
          nframe=nframe+1
          nexttime=nexttime+frametime
        endif

        if (nframe > lastframe .or. error /= 0) exit ! end the integration loop
        
        if (is_number_of_steps_detection_enabled.GT.0) then
            number_of_steps_innerloop = number_of_steps_innerloop +1
            if (number_of_steps_innerloop.GE.max_number_of_steps) then
              stopping_index = next_index_for_stopping_condition()
              error = set_stopping_condition_info(stopping_index, NUMBER_OF_STEPS_DETECTION)
              exit
            endif
        endif
        if (is_timeout_detection_enabled.GT.0) then
            call SYSTEM_CLOCK(clock_current, count_rate, count_max)
            if ((clock_current-clock_init).GE.timeout) then
              stopping_index = next_index_for_stopping_condition()
              error = set_stopping_condition_info(stopping_index, TIMEOUT_DETECTION)
              exit
            endif
        endif
        
      end do
    endif
    ret=0
    if(abs(time-tend).LT.1.e-15) time=tend
    if(error /= 0) then
        ret = -2
        return
    end if
    !if(time.NE.tend) then
    !  print*,time.LE.tend,time-tend,time,tend,nf,lastframe,frametime
    !  ret=-1
    !endif 
  end function

  function amuse_endrun() result(ret)
    integer :: ret
    call mpi_end()
    ret=0
  end function

  function amuse_init_mesh(nx,ny,nz) result(ret)
    integer :: ret,nx,ny,nz

    meshx=nx
    meshy=ny
    meshz=nz

    
    ret=0
  end function
  
  
  function amuse_get_mesh(nx,ny,nz) result(ret)
    integer :: ret,nx,ny,nz
    
    nx = meshx
    ny = meshy
    nz = meshz
    
    ret=0
  end function

  function amuse_init_coords(xlen,ylen,zlen) result(ret)
    integer :: ret,i,j,k
    real*8 :: xlen,ylen,zlen
    character(len=10),parameter :: str_length_unit="default"

    xlength=xlen
    ylength=ylen
    zlength=zlen
    
    ret=0

  end function

  function amuse_commit_parameters() result(ret)
    implicit none
    integer :: i,j,k,ret
    character(len=10),parameter :: str_length_unit="default"

    ret = 0
    write(unit=log_unit,fmt="(2/,A,/)") "----- Mesh -----"
    write(unit=log_unit,fmt="(A,3I5)") "1) Number of mesh points: ", &
      meshx,meshy,meshz
    
    if( meshx.eq.0 .or. meshy.eq.0 .or. meshz.eq.0) then
        ret = -1
        return
    end if
    call fnd3ddecomp()

    write(unit=log_unit,fmt=*) "Local mesh: ",sx,ex,sy,ey,sz,ez
    write(unit=log_unit,fmt="(a,3(e10.3),a)") & 
      "2) Size of grid box : ", &
      xlength,ylength,zlength,str_length_unit

    call init_grid()
    ret = init_boundary()
    if(RET.NE.0) return
    
    
    ret=amuse_init_hydro()  
    if(RET.NE.0) return
    
    xlength=xlength/scleng
    ylength=ylength/scleng
    zlength=zlength/scleng

    dx=xlength/real(max(1,meshx),dp)
    dy=ylength/real(max(1,meshy),dp)
    dz=zlength/real(max(1,meshz),dp)

    do i=sx-mbc,ex+mbc            
      x(i)=dx*(real(i-1,dp)+0.5_dp)
    enddo
    do j=sy-mbc,ey+mbc
      y(j)=dy*(real(j-1,dp)+0.5_dp)
    enddo
    do k=sz-mbc,ez+mbc
      z(k)=dz*(real(k-1,dp)+0.5_dp)
    enddo

    xedge(1:2)=1.0
    yedge(1:2)=1.0
    zedge(1:2)=1.0
    
    if (domainboundaryconditions(1,1).eq.PROBLEM_DEF) then
        innerxpressure=amuse_get_pressure(innerxstate)
    endif
    
    if (domainboundaryconditions(1,2).eq.PROBLEM_DEF) then
        outerxpressure=amuse_get_pressure(outerxstate)
    endif
    
    if (domainboundaryconditions(2,1).eq.PROBLEM_DEF) then
        innerypressure=amuse_get_pressure(innerystate)
    endif
    if (domainboundaryconditions(2,2).eq.PROBLEM_DEF) then
        outerypressure=amuse_get_pressure(outerystate)
    endif
    if (domainboundaryconditions(3,1).eq.PROBLEM_DEF) then
        innerzpressure=amuse_get_pressure(innerzstate)
    endif
    if (domainboundaryconditions(3,2).eq.PROBLEM_DEF) then
        outerzpressure=amuse_get_pressure(outerzstate)
    endif
    
  end function
  
  function amuse_init_hydro() result(ret)
    integer :: ret
    call init_hydro()
    ret=0
  end function

  function amuse_get_pressure(istate) result(ret)
    real*8 :: ret
    real*8 :: istate(neq)    
    ret = gamma1*(istate(EN)- &
      (istate(RHVX)**2+istate(RHVY)**2+istate(RHVZ)**2)/istate(RHO))
      
  end function  
  
  function amuse_set_boundary_innerxstate(istate) result(ret)
    integer :: ret
    real*8 :: istate(neq)    
    ret=0
    innerxstate=istate 
    if(amuse_get_pressure(istate).LE.0) ret=-1
  end function  

  function amuse_set_boundary_innerystate(istate) result(ret)
    integer :: ret
    real*8 :: istate(neq)    
    ret=0
    innerystate=istate  
    innerypressure=amuse_get_pressure(innerystate)
    if(amuse_get_pressure(istate).LE.0) ret=-1
  end function  

  function amuse_set_boundary_innerzstate(istate) result(ret)
    integer :: ret
    real*8 :: istate(neq)    
    ret=0
    innerzstate=istate  
    innerzpressure=amuse_get_pressure(innerzstate)
    if(amuse_get_pressure(istate).LE.0) ret=-1
  end function  

  function amuse_set_boundary_outerxstate(istate) result(ret)
    integer :: ret
    real*8 :: istate(neq)    
    ret=0
    outerxstate=istate  
    ret=0
    if(amuse_get_pressure(istate).LE.0) ret=-1
  end function  

  function amuse_set_boundary_outerystate(istate) result(ret)
    integer :: ret
    real*8 :: istate(neq)    
    outerystate=istate  
    outerypressure=amuse_get_pressure(outerystate)
    ret=0
    if(amuse_get_pressure(istate).LE.0) ret=-1
  end function  

  function amuse_set_boundary_outerzstate(istate) result(ret)
    integer :: ret
    real*8 :: istate(neq)    
    ret=0
    outerzstate=istate  
    outerzpressure=amuse_get_pressure(outerzstate)
    if(amuse_get_pressure(istate).LE.0) ret=-1
  end function  

  function amuse_set_boundary(lowx,highx,lowy,highy,lowz,highz) result(ret)
    integer :: ret
    character(len=*) :: lowx,highx,lowy,highy,lowz,highz
    logical :: periods(3)=.FALSE.


    select case (lowx)
    case("reflective")
      domainboundaryconditions(1,1)=REFLECTIVE    
    case("ref_shift")
      domainboundaryconditions(1,1)=REFLECTIVE_SHIFT    
    case("periodic")
      domainboundaryconditions(1,1)=PERIODIC
      periods(1)=.TRUE.
    case("outflow")
      domainboundaryconditions(1,1)=OUTFLOW          
    case("interface")
      domainboundaryconditions(1,1)=PROBLEM_DEF          
    case default
      ret=-1
      return
    end select

    select case (highx)
    case("reflective")
      domainboundaryconditions(1,2)=REFLECTIVE    
    case("ref_shift")
      domainboundaryconditions(1,2)=REFLECTIVE_SHIFT    
    case("periodic")
      domainboundaryconditions(1,2)=PERIODIC
      periods(1)=.TRUE.
    case("outflow")
      domainboundaryconditions(1,2)=OUTFLOW          
    case("interface")
      domainboundaryconditions(1,2)=PROBLEM_DEF          
    case default
      ret=-1
      return
    end select

    select case (lowy)
    case("reflective")
      domainboundaryconditions(2,1)=REFLECTIVE    
    case("ref_shift")
      domainboundaryconditions(2,1)=REFLECTIVE_SHIFT    
    case("periodic")
      domainboundaryconditions(2,1)=PERIODIC
      periods(2)=.TRUE.
    case("outflow")
      domainboundaryconditions(2,1)=OUTFLOW          
    case("interface")
      domainboundaryconditions(2,1)=PROBLEM_DEF          
    case default
      ret=-1
      return
    end select

    select case (highy)
    case("reflective")
      domainboundaryconditions(2,2)=REFLECTIVE    
    case("ref_shift")
      domainboundaryconditions(2,2)=REFLECTIVE_SHIFT    
    case("periodic")
      domainboundaryconditions(2,2)=PERIODIC
      periods(2)=.TRUE.
    case("outflow")
      domainboundaryconditions(2,2)=OUTFLOW          
    case("interface")
      domainboundaryconditions(2,2)=PROBLEM_DEF          
    case default
      ret=-1
      return
    end select

    select case (lowz)
    case("reflective")
      domainboundaryconditions(3,1)=REFLECTIVE    
    case("ref_shift")
      domainboundaryconditions(3,1)=REFLECTIVE_SHIFT    
    case("periodic")
      domainboundaryconditions(3,1)=PERIODIC
      periods(3)=.TRUE.
    case("outflow")
      domainboundaryconditions(3,1)=OUTFLOW          
    case("interface")
      domainboundaryconditions(3,1)=PROBLEM_DEF          
    case default
      ret=-1
      return
    end select

    select case (highz)
    case("reflective")
      domainboundaryconditions(3,2)=REFLECTIVE    
    case("ref_shift")
      domainboundaryconditions(3,2)=REFLECTIVE_SHIFT    
    case("periodic")
      domainboundaryconditions(3,2)=PERIODIC
      periods(3)=.TRUE.
    case("outflow")
      domainboundaryconditions(3,2)=OUTFLOW          
    case("interface")
      domainboundaryconditions(3,2)=PROBLEM_DEF          
    case default
      ret=-1
      return
    end select


#ifdef MPI  
    call set_periods(periods)
#endif


    ret=check_boundaries()
  end function

  function fill_grid(ix,iy,iz,istate) result(ret)
    integer :: ret,ix,iy,iz
    real*8 :: istate(neq)
    ret=0
    if( in_grid(ix,iy,iz) )then
      state(ix,iy,iz,1:neq)=istate(1:neq)
      ret=ret+1
    endif
  end function

  function retrieve_grid(ix,iy,iz,ostate) result(ret)
    integer :: ret,ix,iy,iz
    real*8 :: ostate(neq)
    ret=0
    ostate(1:neq)=0.
    if( in_grid(ix,iy,iz) ) then
      ostate(1:neq)=state(ix,iy,iz,1:neq)
      ret=ret+1
    endif
  end function

  function fill_gforce(ix,iy,iz,force) result(ret)
    integer :: ret,ix,iy,iz
    real*8 :: force(nrOfDim)
    ret=0
    if( in_grid(ix,iy,iz) )then
      gforce(ix,iy,iz,1:nrOfDim)=force(1:nrOfDim)
      ret=ret+1
    endif
  end function

  function retrieve_gforce(ix,iy,iz,force) result(ret)
    integer :: ret,ix,iy,iz
    real*8 :: force(nrOfDim)
    ret=0
    force(1:nrOfDim)=0.
    if( in_grid(ix,iy,iz) ) then
      force(1:nrOfDim)=gforce(ix,iy,iz,1:nrOfDim)
      ret=ret+1
    endif
  end function


  function amuse_initialize_grid(t0) result(ret)
    integer :: ret
    real*8 :: t0

    state(:,:,:,RHO)=state(:,:,:,RHO)/scdens
    state(:,:,:,RHVX)=state(:,:,:,RHVX)/scmome
    state(:,:,:,RHVY)=state(:,:,:,RHVY)/scmome
    state(:,:,:,RHVZ)=state(:,:,:,RHVZ)/scmome
    state(:,:,:,EN)=state(:,:,:,EN)/scener

  ! Fill boundary conditions
    call boundaries(OLD,domainboundaryconditions,problemboundary) 

  ! Initialize the ionic concentrations
    ret=amuse_init_ionic()
    if(ret.NE.0) return

    ret=amuse_init_time(t0)

  end function

  function amuse_init_ionic() result(ret)
    integer :: ret
    ret=0
  end function  

  function amuse_init_time(t0) result(ret)
    integer :: ret
    real*8 :: t0
    time=t0
    dt=0.0d0
    frametime=t0
    LastFrame=0
    ret=0
  end function

  function amuse_get_pos_of_index(ix,iy,iz,xout,yout,zout) result(ret)
    integer :: ret,ix,iy,iz
    real*8 :: xout,yout,zout
    ret=0
    xout=0
    yout=0
    zout=0
    if( in_grid(ix,iy,iz) ) then
        xout=x(ix)
        yout=y(iy)
        zout=z(iz)
        ret=ret+1
    endif  
  end function

  function amuse_get_index_of_pos(xin,yin,zin,ix,iy,iz) result(ret)
    integer :: ret,ix,iy,iz
    real*8 :: xin,yin,zin
    ix=floor(xin/dx)+1
    iy=floor(yin/dy)+1
    iz=floor(zin/dz)+1
    ret=0    
  end function

  function in_grid(i,j,k)
    logical :: in_grid
    integer :: i,j,k
    in_grid=.FALSE.
    if(i.LT.sx.OR.i.GT.ex.OR. &
       j.LT.sy.OR.j.GT.ey.OR. &
       k.LT.sz.OR.k.GT.ez) return
    in_grid=.TRUE.
  end function

  function free_grid() result(ret)
    integer :: ret
    
    if (ALLOCATED(state1)) then
        deallocate(state1,state2,pressr,gforce)
    end if
    if (ALLOCATED(x)) then
        deallocate(x,y,z)
    end if
    
    ret=0
  end function

  function amuse_get_time(tnow) result(ret)
    integer :: ret
    real*8 tnow
    tnow=time
    ret=0
  end function

  function check_boundaries() result(ret)
    integer :: ret,i
    logical :: periods(nrOfDim)

#ifdef MPI  
    call get_periods(periods)
#endif
  
    ret=0
    do i=1,nrOfDim
#ifdef MPI  
      if(periods(i)) then
        if(domainboundaryconditions(i,1).NE.PERIODIC.OR. &
             domainboundaryconditions(i,2).NE.PERIODIC) then
           ret=-1
           return
        endif         
      else
        if(domainboundaryconditions(i,1).EQ.PERIODIC.OR. &
             domainboundaryconditions(i,2).EQ.PERIODIC) then
           ret=-2
           return
        endif    
      endif
#else
      if(domainboundaryconditions(i,1).EQ.PERIODIC.AND. &
             domainboundaryconditions(i,2).NE.PERIODIC) then
         ret=-1
         return
      endif         
      if(domainboundaryconditions(i,1).NE.PERIODIC.AND. &
             domainboundaryconditions(i,2).EQ.PERIODIC) then
         ret=-2
         return
      endif         
#endif
    enddo
    
  end function  
  
end module





