module problem
  use sizes, only: mbc,neq,RHO,RHVX,RHVY,RHVZ,EN,nrofDim
  use precision, only: dp
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
  
  integer, dimension(nrofDim,2) :: domainboundaryconditions

  contains

  subroutine problemboundary(boundary_id,newold)
    integer,intent(in) :: boundary_id
    integer,intent(in) :: newold
    
    integer :: i,j,k
    real(kind=dp) :: pres
    
    ! Point state to appropriate array
    state => set_state_pointer(newold)

    select case (boundary_id)
    case (X_IN)
       do k=sz-mbc,ez+mbc
          do j=sy-mbc,ey+mbc
             do i=sx-mbc,sx-1
                state(i,j,k,RHO)=innerxstate(RHO)
                state(i,j,k,RHVX)=innerxstate(RHVX)
                state(i,j,k,RHVY)=innerxstate(RHVY)
                state(i,j,k,RHVZ)=innerxstate(RHVZ)
                state(i,j,k,EN)=innerxstate(EN)
                pressr(i,j,k)=innerxpressure
             enddo
          enddo
       enddo
    case (X_OUT)
       do k=sz-mbc,ez+mbc
          do j=sy-mbc,ey+mbc
             do i=ex+1,ex+mbc
                state(i,j,k,RHO)=outerxstate(RHO)
                state(i,j,k,RHVX)=outerxstate(RHVX)
                state(i,j,k,RHVY)=outerxstate(RHVY)
                state(i,j,k,RHVZ)=outerxstate(RHVZ)
                state(i,j,k,EN)=outerxstate(EN)
                pressr(i,j,k)=outerxpressure
             enddo
          enddo
       enddo
    case (Y_IN)
       do k=sz-mbc,ez+mbc
          do j=sy-mbc,sy-1
             do i=sx-mbc,ex+mbc
                state(i,j,k,RHO)=innerystate(RHO)
                state(i,j,k,RHVX)=innerystate(RHVX)
                state(i,j,k,RHVY)=innerystate(RHVY)
                state(i,j,k,RHVZ)=innerystate(RHVZ)
                state(i,j,k,EN)=innerystate(EN)
                pressr(i,j,k)=innerypressure
             enddo
          enddo
       enddo
    case (Y_OUT)
       do k=sz-mbc,ez+mbc
          do j=ey+1,ey+mbc
             do i=sx-mbc,ex+mbc
                state(i,j,k,RHO)=outerystate(RHO)
                state(i,j,k,RHVX)=outerystate(RHVX)
                state(i,j,k,RHVY)=outerystate(RHVY)
                state(i,j,k,RHVZ)=outerystate(RHVZ)
                state(i,j,k,EN)=outerystate(EN)
                pressr(i,j,k)=outerypressure
             enddo
          enddo
       enddo
    case (Z_IN)
       do k=sz-mbc,sz-1
          do j=sy-mbc,ey+mbc
             do i=sx-mbc,ex+mbc
                state(i,j,k,RHO)=innerzstate(RHO)
                state(i,j,k,RHVX)=innerzstate(RHVX)
                state(i,j,k,RHVY)=innerzstate(RHVY)
                state(i,j,k,RHVZ)=innerzstate(RHVZ)
                state(i,j,k,EN)=innerzstate(EN)
                pressr(i,j,k)=innerzpressure
             enddo
          enddo
       enddo
    case (Z_OUT)
       do k=ez+1,ez+mbc
          do j=sy-mbc,ey+mbc
             do i=sx-mbc,ex+mbc
                state(i,j,k,RHO)=outerzstate(RHO)
                state(i,j,k,RHVX)=outerzstate(RHVX)
                state(i,j,k,RHVY)=outerzstate(RHVY)
                state(i,j,k,RHVZ)=outerzstate(RHVZ)
                state(i,j,k,EN)=outerzstate(EN)
                pressr(i,j,k)=outerzpressure
             enddo
          enddo
       enddo
    end select
        
  end subroutine

  subroutine apply_grav_force(dt,newold)
                  
    real(kind=dp),intent(in) :: dt
    integer,intent(in) :: newold

    integer :: i,j,k

    ! Point state to appropriate array
    state => set_state_pointer(newold)

    do k=sz-mbc,ez+mbc
       do j=sy-mbc,ey+mbc
          do i=sx-mbc,ex+mbc
             state(i,j,k,RHVX)=state(i,j,k,RHVX)+ &
                  dt*state(i,j,k,RHO)*gforce(i,j,k,1)
             state(i,j,k,RHVY)=state(i,j,k,RHVY)+ &
                  dt*state(i,j,k,RHO)*gforce(i,j,k,2)
             state(i,j,k,RHVZ)=state(i,j,k,RHVZ)+ &
                  dt*state(i,j,k,RHO)*gforce(i,j,k,3)
             state(i,j,k,EN)=state(i,j,k,EN)+ &
                  dt*(state(i,j,k,RHVX)*gforce(i,j,k,1)+ &
                  state(i,j,k,RHVY)*gforce(i,j,k,2)+ &
                  state(i,j,k,RHVZ)*gforce(i,j,k,3))
          enddo
       enddo
    enddo

  end subroutine apply_grav_force

end module
