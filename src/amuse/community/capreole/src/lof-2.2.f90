module oddeven

  ! Module for Capreole (2D)
  ! Author: Garrelt Mellema
  ! Date: 2004-05-12 (previous 2003-08-26)
  ! This module is also accepted by the F compiler
  !
  ! This module contains the routines related to the correction
  ! of the odd-even decoupling / carbuncle phenomenon
  !
  ! This is achieved through a local oscillation filter. The
  ! state is searched for oscillations, after which diffusion
  ! is applied.
  !
  ! This version used allocation in order to limit the use of 
  ! stack memory.
 
  use precision
  use scaling
  use sizes
  use mesh
  use grid
  use atomic
  use geometry
  use hydro
  use times

  private
  
  ! Shock detection parameter
  real(kind=dp),parameter,private :: shockDetectFactor = 0.5_dp
  ! The diffusion coefficient used
  real(kind=dp),parameter,private :: eta=0.1_dp

  integer,private :: i,j,k,ieq,m
  real(kind=dp),private :: deltaPressure

  ! Counters
  integer :: nrOfXCorrections,nrOfYCorrections,nrOfZCorrections

  ! Flag for shock detection
  ! ( exported via state(neq) )
  integer,dimension(:,:,:),allocatable,save :: flag

  ! Diffusive fluxes
  real(kind=dp),dimension(:,:,:,:),allocatable,save :: fdiff
  real(kind=dp),dimension(:,:,:,:),allocatable,save :: gdiff
  real(kind=dp),dimension(:,:,:,:),allocatable,save :: hdiff

  public :: odd_even

contains
  !========================================================================
  
  subroutine odd_even (newold,action)
    
    ! This routine handles the odd-even fix interface with the
    ! integration routine. The actual work is done in the set and
    ! apply routines for the different coordinate directions.
    
    integer,intent(in) :: action ! what to do:
    ! 1: fix x-direction
    ! 2: fix y-direction
    ! 0: fix all directions
    ! other values: do nothing
    integer,intent(in) :: newold

    ! Point state to appropriate array
    state => set_state_pointer(newold)
 
    ! allocate flag and flux arrays if not allocated before
    if (.not.(allocated(flag))) then
       allocate(flag(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc))
       allocate(fdiff(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc,neq))
       !allocate(gdiff(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc,neq))
       !allocate(hdiff(sx-mbc:ex+mbc,sy-mbc:ey+mbc,sz-mbc:ez+mbc,neq))
    endif

    ! Figure out what to do.
    select case (action)
    case (0)
       call set_odd_even_x ()
       !state(sx:ex,sy:ey,sz:ez,neq)=state(sx:ex,sy:ey,sz:ez,neq)+ &
       !     real(flag(sx:ex,sy:ey,sz:ez),dp)
       call apply_odd_even_x ()
       call set_odd_even_y ()
       !state(sx:ex,sy:ey,sz:ez,neq)=state(sx:ex,sy:ey,sz:ez,neq)+ &
       !     real(flag(sx:ex,sy:ey,sz:ez),dp)
       call apply_odd_even_y ()
       call set_odd_even_z ()
       !state(sx:ex,sy:ey,sz:ez,neq)=state(sx:ex,sy:ey,sz:ez,neq)+ &
       !     real(flag(sx:ex,sy:ey,sz:ez),dp)
       call apply_odd_even_z ()
    case (1)
       call set_odd_even_x ()
       call apply_odd_even_x ()
       !state(sx:ex,sy:ey,sz:ez,neq)=state(sx:ex,sy:ey,sz:ez,neq)+ &
       !     real(flag(sx:ex,sy:ey,sz:ez),dp)
    case (2)
       call set_odd_even_y ()
       call apply_odd_even_y ()
       !state(sx:ex,sy:ey,sz:ez,neq)=state(sx:ex,sy:ey,sz:ez,neq)+ &
       !     real(flag(sx:ex,sy:ey,sz:ez),dp)
    case (3)
       call set_odd_even_z ()
       call apply_odd_even_z ()
       !state(sx:ex,sy:ey,sz:ez,neq)=state(sx:ex,sy:ey,sz:ez,neq)+ &
       !     real(flag(sx:ex,sy:ey,sz:ez),dp)
    case default
       ! Do nothing
    end select

  end subroutine odd_even

  !---------------------------------------------------------------------------

  subroutine set_odd_even_x
      
    ! This routine find shocks in the perpendicular direction,
    ! then looks for density oscillations along those shocks,
    ! and calculates a diffusive flux if found.
    
    integer :: nrOfFlags,cellCount
    integer :: imin,iplus
    
    ! Reset flags
    flag(:,:,:) = 0
    cellCount=0
    nrOfXCorrections=0
    
    ! Reset diffusive flux
    fdiff(:,:,:,:)=0.0d0
    
    ! Find the shocks in y and z direction, set flag to 1.
    nrOfFlags=detect_shock_y (1)
    nrOfFlags=detect_shock_z (1)
    
    ! Test for odd-even pattern.
    ! Use 4-point up and down pattern for this
    ! Only for flagged cells with shock in y/z-direction.
    ! Set flag to 2 in those zones.
    do k=sz,ez
       do j=sy,ey
          do i=sx,ex
             cellCount = 0
             do m=-2,1
                if (flag(i+m,j,k) == 1 .or. flag(i+m,j,k) == 2)  &
                     cellCount = cellCount + 1
             enddo
             if(cellCount == 4) then
                if(state(i-2,j,k,RHO) > state(i-1,j,k,RHO).and. &
                     state(i-1,j,k,RHO) < state(i,j,k  ,RHO).and. &
                     state(i,j,k  ,RHO) > state(i+1,j,k,RHO)) &
                     flag(i,j,k) = 2
                if(state(i-2,j,k,RHO) < state(i-1,j,k,RHO).and. &
                     state(i-1,j,k,RHO) > state(i,j,k,RHO).and.  & 
                     state(i,j,k  ,RHO) < state(i+1,j,k,RHO)) &
                     flag(i,j,k) = 2
!               if(cellCount == 5) then
!                  if(state(i-2,j,k,RHO) > state(i-1,j,k,RHO).and. &
!                       state(i-1,j,k,RHO) < state(i,j,k  ,RHO).and. &
!                       state(i,j,k  ,RHO) > state(i+1,j,k,RHO).and. &
!                       state(i+1,j,k,RHO) < state(i+2,j,k,RHO)) &
!                       flag(i,j,k) = 2
!                  if(state(i-2,j,k,RHO) < state(i-1,j,k,RHO).and. &
!                       state(i-1,j,k,RHO) > state(i,j,k,RHO).and.  & 
!                       state(i,j,k  ,RHO) < state(i+1,j,k,RHO).and.  &
!                       state(i+1,j,k,RHO) > state(i+2,j,k,RHO)) &
!                       flag(i,j,k) = 2
             endif
          enddo
       enddo
    enddo
    
    ! Calculate diffusive flux for the cells flagged with a 2
    do k=sz,ez
       do j=sy,ey
          do i=sx,ex
             if ( flag(i,j,k) == 2 ) then
                nrOfXCorrections = nrOfXCorrections + 1
                imin=max(i-1,1)      ! diffusive flux at grid edge is zero
                iplus=min(i+1,meshx) ! diffusive flux at grid edge is zero
                do ieq=1,neq
                   fdiff(i,j,k,ieq)=state(imin,j,k,ieq)-state(i,j,k,ieq)
                   fdiff(i+1,j,k,ieq)=state(i,j,k,ieq)-state(iplus,j,k,ieq)
                end do
             end if
          end do
       end do
    end do
    
  end subroutine set_odd_even_x

!------------------------------------------------------------------------------

  subroutine apply_odd_even_x
    
    ! Apply the diffusive flux
    
    if (nrOfXCorrections > 0) then
       do ieq=1,neq
          do k=sz,ez
             do j=sy,ey
                do i=sx,ex
                   state(i,j,k,ieq)=state(i,j,k,ieq)+eta* &
                        (fdiff(i,j,k,ieq)-fdiff(i+1,j,k,ieq))
                enddo
             enddo
          enddo
       enddo
       
    end if
    
  end subroutine apply_odd_even_x
  
  !--------------------------------------------------------------------------
  
  subroutine set_odd_even_y 
    
    ! This routine find shocks in the perpendicular direction,
    ! then looks for density oscillations along those shocks,
    ! and calculates a diffusive flux if found.
    
    integer :: nrOfFlags,cellCount
    integer :: jmin,jplus
    
    ! Reset flags
    flag(:,:,:) = 0
    cellCount=0
    nrOfYCorrections=0
    
    ! Reset diffusive flux
    fdiff(:,:,:,:)=0.0d0
    
    ! Find the shocks in x and z direction, set flag to 3
    nrOfFlags=detect_shock_x (3)
    nrOfFlags=detect_shock_z (3)
    
    ! Test for odd-even pattern in the y-direction.
    ! Use 4-point up and down pattern for this
    ! Only for flagged cells with shock in x/z-direction.
    ! Set flag to 4 in those zones.
    do k=sz,ez
       do j=sy,ey
          do i=sx,ex
             cellCount = 0
             do m=-2,1
                if (flag(i,j+m,k) == 3 .or. flag(i,j+m,k) == 4) &
                     cellCount = cellCount + 1
             enddo
             if(cellCount == 4) then
                if(state(i,j-2,k,RHO) > state(i,j-1,k,RHO).and. &
                     state(i,j-1,k,RHO) < state(i,j,k  ,RHO).and. &
                     state(i,j,k  ,RHO) > state(i,j+1,k,RHO)) &
                     flag(i,j,k) = 4
                if(state(i,j-2,k,RHO) < state(i,j-1,k,RHO).and. &
                     state(i,j-1,k,RHO) > state(i,j,k,RHO).and.  & 
                     state(i,j,k  ,RHO) < state(i,j+1,k,RHO)) &
                     flag(i,j,k) = 4
             endif
!                  if(state(i,j-1,k,RHO) < state(i,j,k  ,RHO).and. &
!                       state(i,j,k  ,RHO) > state(i,j+1,k,RHO).and. &
!                       state(i,j+1,k,RHO) < state(i,j+2,k,RHO)) &
!                       flag(i,j,k) = 4
!                  if(state(i,j-1,k,RHO) > state(i,j,k,RHO).and.  & 
!                       state(i,j,k  ,RHO) < state(i,j+1,k,RHO).and.  &
!                       state(i,j+1,k,RHO) > state(i,j+2,k,RHO)) &
!                       flag(i,j,k) = 4

!               if(cellCount == 5) then
!                  if(state(i,j-2,k,RHO) > state(i,j-1,k,RHO).and. &
!                       state(i,j-1,k,RHO) < state(i,j,k  ,RHO).and. &
!                       state(i,j,k  ,RHO) > state(i,j+1,k,RHO).and. &
!                       state(i,j+1,k,RHO) < state(i,j+2,k,RHO)) &
!                       flag(i,j,k) = 4
!                  if(state(i,j-2,k,RHO) < state(i,j-1,k,RHO).and. &
!                       state(i,j-1,k,RHO) > state(i,j,k,RHO).and.  & 
!                       state(i,j,k  ,RHO) < state(i,j+1,k,RHO).and. &
!                       state(i,j+1,k,RHO) > state(i,j+2,k,RHO)) &
!                       flag(i,j,k) = 4
!               endif

          enddo
       enddo
    enddo
    
    ! Calculate diffusive flux for the cells flagged with a 4
    do k=sz,ez
       do j=sy,ey
          do i=sx,ex
             if (flag(i,j,k) == 4) then
                nrOfYCorrections = nrOfYCorrections + 1
                jmin=max(j-1,1)      ! diffusive flux at grid edge is zero
                jplus=min(j+1,meshy) ! diffusive flux at grid edge is zero
                do ieq=1,neq
                   fdiff(i,j,k,ieq)=state(i,jmin,k,ieq)-state(i,j,k,ieq)
                   fdiff(i,j+1,k,ieq)=state(i,j,k,ieq)-state(i,jplus,k,ieq)
                end do
             end if
          end do
       end do
    end do
    
  end subroutine set_odd_even_y
  
  !--------------------------------------------------------------------------
  
  subroutine apply_odd_even_y ()
    
    ! Apply the diffusive flux
    
    if (nrOfYCorrections > 0) then
       do ieq=1,neq
          do k=sz,ez
             do j=sy,ey
                do i=sx,ex
                   state(i,j,k,ieq)=state(i,j,k,ieq)+eta* &
                        (fdiff(i,j,k,ieq)-fdiff(i,j+1,k,ieq))
                enddo
             enddo
          enddo
       enddo
       
    end if
    
  end subroutine apply_odd_even_y
  
  !--------------------------------------------------------------------------
  
  subroutine set_odd_even_z 
    
    ! This routine find shocks in the perpendicular direction,
    ! then looks for density oscillations along those shocks,
    ! and calculates a diffusive flux if found.
    
    integer :: nrOfFlags,cellCount
    integer :: kmin,kplus
    
    ! Reset flags
    flag(:,:,:) = 0
    cellCount=0
    nrOfZCorrections=0
    
    ! Reset diffusive flux
    fdiff(:,:,:,:)=0.0d0
    
    ! Find the shocks in x and y direction, set flag to 5
    nrOfFlags=detect_shock_x (5)
    nrOfFlags=detect_shock_y (5)
    
    ! Test for odd-even pattern.
    ! Use 4-point up and down pattern for this
    ! Check only for flagged cells with shock in x/y-direction.
    ! Set flag to 6 in those zones.
    do k=sz,ez
       do j=sy,ey
          do i=sx,ex
             cellCount = 0
             do m=-2,1
                if (flag(i,j,k+m) == 5 .or. flag(i,j,k+m) == 6) &
                     cellCount = cellCount + 1
             enddo
             if(cellCount == 4) then
                if(state(i,j,k-2,RHO) > state(i,j,k-1,RHO).and. &
                     state(i,j,k-1,RHO) < state(i,j,k  ,RHO).and. &
                     state(i,j,k  ,RHO) > state(i,j,k+1,RHO)) &
                     flag(i,j,k) = 6
                if(state(i,j,k-2,RHO) < state(i,j,k-1,RHO).and. &
                     state(i,j,k-1,RHO) > state(i,j,k,RHO).and.  & 
                     state(i,j,k  ,RHO) < state(i,j,k+1,RHO)) &
                     flag(i,j,k) = 6
             endif
             
!                  if(state(i,j,k-1,RHO) < state(i,j,k  ,RHO).and. &
!                       state(i,j,k  ,RHO) > state(i,j,k+1,RHO).and. &
!                       state(i,j,k+1,RHO) < state(i,j,k+2,RHO)) &
!                       flag(i,j,k) = 6
!                  if(state(i,j,k-1,RHO) > state(i,j,k,RHO).and.  & 
!                       state(i,j,k  ,RHO) < state(i,j,k+1,RHO).and.  &
!                       state(i,j,k+1,RHO) > state(i,j,k+2,RHO)) &
!                       flag(i,j,k) = 6

!               if(cellCount == 5) then
!                  if(state(i,j,k-2,RHO) > state(i,j,k-1,RHO).and. &
!                       state(i,j,k-1,RHO) < state(i,j,k  ,RHO).and. &
!                       state(i,j,k  ,RHO) > state(i,j,k+1,RHO).and. &
!                       state(i,j,k+1,RHO) < state(i,j,k+2,RHO)) &
!                       flag(i,j,k) = 6
!                  if(state(i,j,k-2,RHO) < state(i,j,k-1,RHO).and. &
!                       state(i,j,k-1,RHO) > state(i,j,k,RHO).and.  & 
!                       state(i,j,k  ,RHO) < state(i,j,k+1,RHO).and. &
!                       state(i,j,k+1,RHO) > state(i,j,k+2,RHO)) &
!                       flag(i,j,k) = 6
!               endif

          enddo
       enddo
    enddo
    
    ! Calculate diffusive flux for the cells flagged with a 6
    do k=sz,ez
       do j=sy,ey
          do i=sx,ex
             if (flag(i,j,k) == 6) then
                nrOfZCorrections = nrOfZCorrections + 1
                kmin=max(k-1,1)      ! diffusive flux at grid edge is zero
                kplus=min(k+1,meshz) ! diffusive flux at grid edge is zero
                do ieq=1,neq
                   fdiff(i,j,k,ieq)=state(i,j,kmin,ieq)-state(i,j,k,ieq)
                   fdiff(i,j,k+1,ieq)=state(i,j,k,ieq)-state(i,j,kplus,ieq)
                end do
             end if
          end do
       end do
    end do
    
  end subroutine set_odd_even_z
  
  !--------------------------------------------------------------------------
  
  subroutine apply_odd_even_z ()
    
    ! Apply the diffusive flux
    
    if (nrOfZCorrections > 0) then
       do ieq=1,neq
          do k=sz,ez
             do j=sy,ey
                do i=sx,ex
                   state(i,j,k,ieq)=state(i,j,k,ieq)+eta* &
                        (fdiff(i,j,k,ieq)-fdiff(i,j,k+1,ieq))
                enddo
             enddo
          enddo
       enddo
       
    end if
    
  end subroutine apply_odd_even_z
  
  !--------------------------------------------------------------------------
  
  function detect_shock_x (marker)
    
    ! Detect shocks in the x-direction and flag cells with marker
    ! This is done by looking at the pressure jumps.
    
    ! Output: number of cells flagged
    
    integer :: detect_shock_x
    
    integer,intent(in) :: marker
    
    detect_shock_x=0
    
    do k=sz,ez
       do j=sy,ey
          do i=sx,ex
             deltaPressure = (pressr(i+1,j,k)-pressr(i,j,k)) / &
                  max(pressr(i+1,j,k),pressr(i,j,k))
             if(-(deltaPressure) > shockDetectFactor) then
                flag(i-4,j,k) = marker
                flag(i-3,j,k) = marker
                flag(i-2,j,k) = marker
                flag(i-1,j,k) = marker
                flag(i  ,j,k) = marker
                detect_shock_x = detect_shock_x + 5
             endif
             if((deltaPressure) > shockDetectFactor) then
                flag(i+5,j,k) = marker
                flag(i+4,j,k) = marker
                flag(i+3,j,k) = marker
                flag(i+2,j,k) = marker
                flag(i+1,j,k) = marker
                detect_shock_x = detect_shock_x + 5
             endif
          enddo
       enddo
    enddo
    
  end function detect_shock_x
  
  !--------------------------------------------------------------------------
  
  function detect_shock_y (marker)
    
    ! Detect shocks in y-direction and flag cells with marker
    
    ! Output: number of cells flagged
    
    integer :: detect_shock_y
    
    integer,intent(in) :: marker
    
    detect_shock_y=0
    do k=sz,ez
       do j=sy,ey
          do i=sx,ex
             deltaPressure = (pressr(i,j+1,k)-pressr(i,j,k)) / &
                  max(pressr(i,j+1,k),pressr(i,j,k))
             if(-(deltaPressure) > shockDetectFactor) then
                flag(i,j-4,k) = marker
                flag(i,j-3,k) = marker
                flag(i,j-2,k) = marker
                flag(i,j-1,k) = marker
                flag(i  ,j,k) = marker
                detect_shock_y = detect_shock_y + 5
             endif
             if((deltaPressure) > shockDetectFactor) then
                flag(i,j+5,k) = marker
                flag(i,j+4,k) = marker
                flag(i,j+3,k) = marker
                flag(i,j+2,k) = marker
                flag(i,j+1,k) = marker
                detect_shock_y = detect_shock_y + 5
             endif
          enddo
       enddo
    enddo
    
  end function detect_shock_y
  
  !--------------------------------------------------------------------------
  
  function detect_shock_z (marker)
    
    ! detect shocks in the z-direction and flag cells with marker
    
    ! Output: number of cells flagged
    
    integer :: detect_shock_z
      
    integer,intent(in) :: marker
    
    detect_shock_z=0
    do k=sz,ez
       do j=sy,ey
          do i=sx,ex
             deltaPressure = (pressr(i,j,k+1)-pressr(i,j,k)) / &
                  max(pressr(i,j,k+1),pressr(i,j,k))
             if(-(deltaPressure) > shockDetectFactor) then
                flag(i,j,k-4) = marker
                flag(i,j,k-3) = marker
                flag(i,j,k-2) = marker
                flag(i,j,k-1) = marker
                flag(i  ,j,k) = marker
                detect_shock_z = detect_shock_z + 5
             endif
             if((deltaPressure) > shockDetectFactor) then
                flag(i,j,k+5) = marker
                flag(i,j,k+4) = marker
                flag(i,j,k+3) = marker
                flag(i,j,k+2) = marker
                flag(i,j,k+1) = marker
                detect_shock_z = detect_shock_z + 5
             endif
          enddo
       enddo
    enddo
    
  end function detect_shock_z

end module oddeven
