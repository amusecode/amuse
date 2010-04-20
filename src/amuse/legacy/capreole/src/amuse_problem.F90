module problem
  use sizes, only: mbc,neq,RHO,RHVX,RHVY,RHVZ,EN,nrofDim
  use precision, only: dp
  use hydro, only:  state,pressr,set_state_pointer,NEW,OLD,restart_state,gforce
  use mesh, only: sx,ex,sy,ey,sz,ez,meshx,meshy,meshz
  
  integer, dimension(nrofDim,2) :: domainboundaryconditions

  contains

  subroutine problemboundary(boundary_id,newold)
    integer,intent(in) :: boundary_id
    integer,intent(in) :: newold
! dummy
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
