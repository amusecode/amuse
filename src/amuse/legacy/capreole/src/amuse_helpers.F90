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
  
  contains

  function amuse_init() result(ret)
    integer :: ret
!  call setup_clocks ()
    call mpi_setup()
    ret=0    
  end function

  function amuse_evolve(tend) result(ret)
    integer :: ret,nf
    real*8 :: tend

    if(tend-time.GT.0) then
      lastframe=0
      nf=1
      do while(nf.GE.lastframe)
        lastframe=lastframe+1
        frametime=tend/lastframe
        nf=nint(time/frametime)
      enddo
      call evolve()
    endif
    ret=0
    if(time.NE.tend) ret=-1
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

    write(unit=log_unit,fmt="(2/,A,/)") "----- Mesh -----"
    write(unit=log_unit,fmt="(A,3I5)") "1) Number of mesh points: ", &
      meshx,meshy,meshz

    call fnd3ddecomp()

    write(unit=log_unit,fmt=*) "Local mesh: ",sx,ex,sy,ey,sz,ez
    ret=0
  end function

  function amuse_init_coords(xlen,ylen,zlen) result(ret)
    integer :: ret,i,j,k
    real*8 :: xlen,ylen,zlen
    character(len=10),parameter :: str_length_unit="default"

    xlength=xlen
    ylength=ylen
    zlength=zlen

    write (unit=log_unit,fmt="(a,3(e10.3),a)") & 
      "2) Size of grid box : ", &
      xlength,ylength,zlength,str_length_unit

    call init_grid()

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

    ret=0

  end function

  function amuse_init_hydro() result(ret)
    integer :: ret
    call init_hydro()
    ret=0
  end function

  function amuse_set_boundary_innerxstate(istate) result(ret)
    integer :: ret
    real*8 :: istate(neq)    
    ret=0
    innerxstate=istate  
    innerxpressure=gamma1*(istate(EN)- &
      (istate(RHVX)**2+istate(RHVY)**2+istate(RHVZ)**2)/istate(RHO))
    if(innerxpressure.LE.0) ret=-1
    ret=0
  end function  

  function amuse_set_boundary_innerystate(istate) result(ret)
    integer :: ret
    real*8 :: istate(neq)    
    ret=0
    innerystate=istate  
    innerypressure=gamma1*(istate(EN)- &
      (istate(RHVX)**2+istate(RHVY)**2+istate(RHVZ)**2)/istate(RHO))
    if(innerypressure.LE.0) ret=-1
    ret=0
  end function  

  function amuse_set_boundary_innerzstate(istate) result(ret)
    integer :: ret
    real*8 :: istate(neq)    
    ret=0
    innerzstate=istate  
    innerzpressure=gamma1*(istate(EN)- &
      (istate(RHVX)**2+istate(RHVY)**2+istate(RHVZ)**2)/istate(RHO))
    if(innerzpressure.LE.0) ret=-1
    ret=0
  end function  

  function amuse_set_boundary_outerxstate(istate) result(ret)
    integer :: ret
    real*8 :: istate(neq)    
    ret=0
    outerxstate=istate  
    outerxpressure=gamma1*(istate(EN)- &
      (istate(RHVX)**2+istate(RHVY)**2+istate(RHVZ)**2)/istate(RHO))
    if(outerxpressure.LE.0) ret=-1
    ret=0
  end function  

  function amuse_set_boundary_outerystate(istate) result(ret)
    integer :: ret
    real*8 :: istate(neq)    
    ret=0
    outerystate=istate  
    outerypressure=gamma1*(istate(EN)- &
      (istate(RHVX)**2+istate(RHVY)**2+istate(RHVZ)**2)/istate(RHO))
    if(outerypressure.LE.0) ret=-1
    ret=0
  end function  

  function amuse_set_boundary_outerzstate(istate) result(ret)
    integer :: ret
    real*8 :: istate(neq)    
    ret=0
    outerzstate=istate  
    outerzpressure=gamma1*(istate(EN)- &
      (istate(RHVX)**2+istate(RHVY)**2+istate(RHVZ)**2)/istate(RHO))
    if(outerzpressure.LE.0) ret=-1
    ret=0
  end function  

  function amuse_set_boundary(lowx,highx,lowy,highy,lowz,highz) result(ret)
    integer :: ret
    character(len=10) :: lowx,highx,lowy,highy,lowz,highz
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
    
    deallocate(state1,state2,pressr,gforce)
    deallocate(x,y,z)
    
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





