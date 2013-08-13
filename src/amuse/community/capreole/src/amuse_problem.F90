module problem
  use sizes, only: mbc,neq,RHO,RHVX,RHVY,RHVZ,EN,nrofDim
  use precision, only: dp
  use atomic, only : gamma1
  use hydro, only:  state,pressr,set_state_pointer,NEW,OLD,restart_state,gforce
  use mesh, only: sx,ex,sy,ey,sz,ez,meshx,meshy,meshz
  use boundary, only: boundaries,REFLECTIVE,OUTFLOW,PROBLEM_DEF,X_IN,X_OUT,Y_IN, &
       Y_OUT,Z_IN,Z_OUT

  real(kind=dp) :: innerxstate(neq),outerxstate(neq), &
                   innerystate(neq),outerystate(neq), &
                   innerzstate(neq),outerzstate(neq)
  real(kind=dp) :: innerxpressure,outerxpressure, &
                   innerypressure,outerypressure, &
                   innerzpressure,outerzpressure
  
  real(kind=dp),dimension(:,:,:,:),allocatable,target,public :: boundary_left_x1
  real(kind=dp),dimension(:,:,:,:),allocatable,target,public :: boundary_right_x1
  real(kind=dp),dimension(:,:,:,:),allocatable,target,public :: boundary_left_x2
  real(kind=dp),dimension(:,:,:,:),allocatable,target,public :: boundary_right_x2
  real(kind=dp),dimension(:,:,:,:),allocatable,target,public :: boundary_left_x3
  real(kind=dp),dimension(:,:,:,:),allocatable,target,public :: boundary_right_x3
  
  integer, dimension(nrofDim,2) :: domainboundaryconditions

  contains

  function get_boundary_grid_pointer(index_of_boundary) result(ret)
    integer,intent(in) :: index_of_boundary
    real(kind=dp),pointer,dimension(:,:,:,:) :: ret

    ! Point state to appropriate array
    select case (index_of_boundary)
        case(1)
            ret => boundary_left_x1
        case(2)
            ret => boundary_right_x1
        case(3)
            ret => boundary_left_x2
        case(4)
            ret => boundary_right_x2
        case(5)
            ret => boundary_left_x3
        case(6)
            ret => boundary_right_x3
    end select
  end function get_boundary_grid_pointer
  
  
  function is_on_boundary_grid(index_of_boundary, i, j, k, is_set) result(ret)
    implicit none
    integer,intent(in) :: i, j, k
    integer,intent(in) :: index_of_boundary
    logical :: ret
    logical,intent(in) :: is_set
    integer :: il, iu
    integer :: jl, ju
    
    if(is_set .OR. sx .EQ. 1) then
        il = sx - mbc
    else
        il = sx
    endif
        
    if(is_set .OR. ex .EQ. meshx) then
        iu = ex + mbc
    else
        iu = ex
    endif
    
    if(is_set .OR. sy .EQ. 1) then
        jl = sy - mbc
    else
        jl = sy
    endif
        
    if(is_set .OR. ey .EQ. meshy) then
        ju = ey + mbc
    else
        ju = ey
    endif
    
    ! Point state to appropriate array
    ret = .TRUE.
    ! print *, "I:", sx, ex, i, meshx, il, iu, is_set, (i-mbc).LT.il , (i-mbc).GT.iu
    ! print *, "J:", sy, ey, j
    ! print *, "K:", sz, ez, k
    ! print *, "index_of_boundary:",index_of_boundary, j.LT.1 .OR. j.GT.mbc
    
    select case (index_of_boundary)
        case(1)
            if(sx .NE. 1) then
                ret = .FALSE.
            else if(i.LT.1 .OR. i.GT.mbc) then
                ret = .FALSE.
            else if (j.LT.sy.OR.j.GT.ey) then
                ret = .FALSE.
            else if (k.LT.sz.OR.k.GT.ez) then
                ret = .FALSE.
            end if
        case(2)
            if(ex .NE. meshx) then
                ret = .FALSE.
            else if(i.LT.1 .OR. i.GT.mbc) then
                ret = .FALSE.
            else if (j.LT.sy.OR.j.GT.ey) then
                ret = .FALSE.
            else if (k.LT.sz.OR.k.GT.ez) then
                ret = .FALSE.
            end if
        case(3)
            if(sy .NE. 1) then
                ret = .FALSE.
            else if((i-mbc).LT.il .OR. (i-mbc).GT.iu) then
                ret = .FALSE.
            else if(j.LT.1 .OR. j.GT.mbc) then
                ret = .FALSE.
            else if (k.LT.sz.OR.k.GT.ez) then
                ret = .FALSE.
            end if
        case(4)
            if(ey .NE. meshy) then
                ret = .FALSE.
            else if((i-mbc).LT.il .OR. (i-mbc).GT.iu) then
                ret = .FALSE.
            else if(j.LT.1 .OR. j.GT.mbc) then
                ret = .FALSE.
            else if (k.LT.sz.OR.k.GT.ez) then
                ret = .FALSE.
            end if
        case(5)
            if(sz .NE. 1) then
                ret = .FALSE.
            else if((i-mbc).LT.il .OR. (i-mbc).GT.iu) then
                ret = .FALSE.
            else if((j-mbc).LT.jl .OR. (j-mbc).GT.ju) then
                ret = .FALSE.
            else if(k.LT.1 .OR. k.GT.mbc) then
                ret = .FALSE.
            end if
        case(6)
            if(ez .NE. meshz) then
                ret = .FALSE.
            else if((i-mbc).LT.il .OR. (i-mbc).GT.iu) then
                ret = .FALSE.
            else if((j-mbc).LT.jl .OR. (j-mbc).GT.ju) then
                ret = .FALSE.
            else if(k.LT.1 .OR. k.GT.mbc) then
                ret = .FALSE.
            end if
    end select
    ! print *, "RET:", ret
  end function is_on_boundary_grid
  
  
  function on_boundary_grid(index_of_boundary, i, j, k) result(ret)
    implicit none
    integer,intent(inout) :: i, j, k
    integer,intent(in) :: index_of_boundary
    integer :: ret
    ret = 0
    
    select case (index_of_boundary)
        case(3)
            i = i - mbc
        case(4)
            i = i - mbc
        case(5)
            i = i - mbc
            j = j - mbc
        case(6)
            i = i - mbc
            j = j - mbc
    end select
  end function on_boundary_grid
  
  function init_boundary()
    integer :: init_boundary
    integer :: status
    allocate(boundary_left_x1(1:mbc,sy:ey,sz:ez,neq), STAT=status)
    allocate(boundary_right_x1(1:mbc,sy:ey,sz:ez,neq), STAT=status)
    allocate(boundary_left_x2(sx-mbc:ex+mbc,1:mbc,sz:ez,neq), STAT=status)
    allocate(boundary_right_x2(sx-mbc:ex+mbc,1:mbc,sz:ez,neq), STAT=status)
    allocate(boundary_left_x3(sx-mbc:ex+mbc,sy-mbc:ey+mbc,1:mbc,neq), STAT=status)
    allocate(boundary_right_x3(sx-mbc:ex+mbc,sy-mbc:ey+mbc,1:mbc,neq), STAT=status)
    if(status.NE.0) then
      init_boundary = -1
    else
      init_boundary = 0
    end if
  end function
  
  function get_pressure(istate) result(ret)
    real*8 :: ret
    real*8 :: istate(neq)    
    ret = gamma1*(istate(EN)- &
      0.5*(istate(RHVX)**2+istate(RHVY)**2+istate(RHVZ)**2)/istate(RHO))
      
  end function  
  
  subroutine problemboundary(boundary_id,newold)
    integer,intent(in) :: boundary_id
    integer,intent(in) :: newold
    
    integer :: i,j,k
    real(kind=dp) :: pres
    
    ! Point state to appropriate array
    state => set_state_pointer(newold)

    select case (boundary_id)
    case (X_IN)
       do k=sz,ez
          do j=sy,ey
             do i=1, mbc
                state(sx-i,j,k,RHO)=boundary_left_x1(mbc - i + 1, j, k, RHO)
                state(sx-i,j,k,RHVX)=boundary_left_x1(mbc - i + 1, j, k, RHVX)
                state(sx-i,j,k,RHVY)=boundary_left_x1(mbc - i + 1, j, k, RHVY)
                state(sx-i,j,k,RHVZ)=boundary_left_x1(mbc - i + 1, j, k, RHVZ)
                state(sx-i,j,k,EN)=boundary_left_x1(mbc - i + 1, j, k, EN)
                pressr(sx-i,j,k)=get_pressure(boundary_left_x1(mbc - i + 1, j, k,:))
             enddo
          enddo
       enddo
    case (X_OUT)
       do k=sz, ez
          do j=sy, ey
             do i=1,mbc
                state(ex+i,j,k,RHO)=boundary_right_x1(i, j, k, RHO)
                state(ex+i,j,k,RHVX)=boundary_right_x1(i, j, k, RHVX)
                state(ex+i,j,k,RHVY)=boundary_right_x1(i, j, k, RHVY)
                state(ex+i,j,k,RHVZ)=boundary_right_x1(i, j, k, RHVZ)
                state(ex+i,j,k,EN)=boundary_right_x1(i, j, k, EN)
                pressr(ex+i,j,k)=get_pressure(boundary_right_x1(i, j, k,:))
             enddo
          enddo
       enddo
    case (Y_IN)
       do k=sz, ez
          do j=1, mbc
             do i=sx-mbc,ex+mbc
                state(i,sy-j,k,RHO)=boundary_left_x2(i, mbc - j + 1, k, RHO)
                state(i,sy-j,k,RHVX)=boundary_left_x2(i, mbc - j + 1, k, RHVX)
                state(i,sy-j,k,RHVY)=boundary_left_x2(i, mbc - j + 1, k, RHVY)
                state(i,sy-j,k,RHVZ)=boundary_left_x2(i, mbc - j + 1, k, RHVZ)
                state(i,sy-j,k,EN)=boundary_left_x2(i, mbc - j + 1, k, EN)
                pressr(i,sy-j,k)=get_pressure(boundary_left_x2(i, mbc - j + 1, k,:))
             enddo
          enddo
       enddo
    case (Y_OUT)
       do k=sz, ez
          do j=1, mbc
             do i=sx-mbc,ex+mbc
                state(i,ey+j,k,RHO)=boundary_right_x2(i, j, k, RHO)
                state(i,ey+j,k,RHVX)=boundary_right_x2(i, j, k, RHVX)
                state(i,ey+j,k,RHVY)=boundary_right_x2(i, j, k, RHVY)
                state(i,ey+j,k,RHVZ)=boundary_right_x2(i, j, k, RHVZ)
                state(i,ey+j,k,EN)=boundary_right_x2(i, j, k, EN)
                pressr(i,ey+j,k)=get_pressure(boundary_right_x2(i, j, k,:))
             enddo
          enddo
       enddo
    case (Z_IN)
       do k=1, mbc
          do j=sy-mbc,ey+mbc
             do i=sx-mbc,ex+mbc
                state(i,j,sz-k,RHO)=boundary_left_x3(i, j, mbc - k + 1, RHO)
                state(i,j,sz-k,RHVX)=boundary_left_x3(i, j, mbc - k + 1, RHVX)
                state(i,j,sz-k,RHVY)=boundary_left_x3(i, j, mbc - k + 1, RHVY)
                state(i,j,sz-k,RHVZ)=boundary_left_x3(i, j, mbc - k + 1, RHVZ)
                state(i,j,sz-k,EN)=boundary_left_x3(i, j, mbc - k + 1, EN)
                pressr(i,j,sz-k)=get_pressure(boundary_left_x3(i, j, mbc - k + 1,:))
             enddo
          enddo
       enddo
    case (Z_OUT)
       do k=1, mbc
          do j=sy-mbc,ey+mbc
             do i=sx-mbc,ex+mbc
                state(i,j,ez+k,RHO)=boundary_right_x3(i, j, k, RHO)
                state(i,j,ez+k,RHVX)=boundary_right_x3(i, j, k, RHVX)
                state(i,j,ez+k,RHVY)=boundary_right_x3(i, j, k, RHVY)
                state(i,j,ez+k,RHVZ)=boundary_right_x3(i, j, k, RHVZ)
                state(i,j,ez+k,EN)=boundary_right_x3(i, j, k, EN)
                pressr(i,j,ez+k)=get_pressure(boundary_right_x3(i, j, k,:))
             enddo
          enddo
       enddo
    end select
        
  end subroutine

  subroutine apply_grav_force(dt,newold)
                  
    real(kind=dp),intent(in) :: dt
    integer,intent(in) :: newold
    real(kind=dp) :: u


    integer :: i,j,k

    ! Point state to appropriate array
    state => set_state_pointer(newold)


    do k=sz-mbc,ez+mbc
       do j=sy-mbc,ey+mbc
          do i=sx-mbc,ex+mbc
             u=state(i,j,k,EN)-(state(i,j,k,RHVX)**2+ &
                                state(i,j,k,RHVY)**2+ &
                                state(i,j,k,RHVZ)**2)/state(i,j,k,RHO)
!             state(i,j,k,EN)=state(i,j,k,EN)+ &
!                  dt*(state(i,j,k,RHVX)*gforce(i,j,k,1)+ &
!                  state(i,j,k,RHVY)*gforce(i,j,k,2)+ &
!                  state(i,j,k,RHVZ)*gforce(i,j,k,3))
             state(i,j,k,RHVX)=state(i,j,k,RHVX)+ &
                  dt*state(i,j,k,RHO)*gforce(i,j,k,1)
             state(i,j,k,RHVY)=state(i,j,k,RHVY)+ &
                  dt*state(i,j,k,RHO)*gforce(i,j,k,2)
             state(i,j,k,RHVZ)=state(i,j,k,RHVZ)+ &
                  dt*state(i,j,k,RHO)*gforce(i,j,k,3)
             state(i,j,k,EN)=u+(state(i,j,k,RHVX)**2+ &
                                state(i,j,k,RHVY)**2+ &
                                state(i,j,k,RHVZ)**2)/state(i,j,k,RHO)    
          enddo
       enddo
    enddo

  end subroutine apply_grav_force

end module
