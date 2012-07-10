module protection

  ! Module for Capreole 
  ! Author: Garrelt Mellema
  ! Date: 2006-12-30
  !
  ! This module contains the routines related to the negative
  ! pressure / density / total energy correction

  ! Memory saving version

  use precision, only:dp
  use sizes, only: mbc,neq, neuler, RHO, EN
  use mesh, only: meshx, meshy, meshz, sx,ex,sy,ey,sz,ez
  use hydro, only: state, pressr, set_state_pointer
  use times, only: time
  use abundances, only: mu
  use atomic, only: boltzm, gamma1
  use geometry, only: presfunc
#ifdef __IBMC__
  USE XLFUTILITY, only: flush => flush_
#endif

  implicit none

  private

  public:: presprot

contains
  !========================================================================

  subroutine presprot(inegative,icall,newold)

    ! This routine protects the pressure, density, and energy density 
    ! from becoming negative.
    ! We do this in two stages, if 1) fails, do 2).
    ! 1) diffuse  
    ! 2) setting a minimal temperature tmin
    
    ! the diffusion coefficient used in stage 1)
    real(kind=dp),parameter :: eta=0.05d0 

    ! the minimum temperature applied in stage 2)
    real(kind=dp),parameter :: tmin=1.0d1

    integer,intent(out) :: inegative ! control integer; /= 0 if fixes failed
    integer,intent(in)  :: icall ! call id (to distinguish between different
                                 !  calls to this subroutine)
    integer,intent(in) :: newold

    integer :: i,j,k,ieq
    real(kind=dp)    :: pnew
    integer :: imin,jmin,iplus,jplus,kmin,kplus
    integer :: ierror
    integer :: problem_counter, nproblem

    ! A structure to contain the position and diffusive
    ! fluxes for a problem
    type problem
       integer,dimension(3) ::position
       real(kind=dp),dimension(2,3,neq) :: dflux
    end type problem

    ! Not more than 1000 problems allowed
    type(problem),dimension(1000) :: problem_list

    ! Point state to appropriate array
    state => set_state_pointer(newold)

    inegative=0

    problem_counter=0
    do k=sz-1,ez+1
       do j=sy-1,ey+1
          do i=sx-1,ex+1
             ! Check for negative pressure/density/energy
             if (pressr(i,j,k) <= 0.0.or.state(i,j,k,RHO) <= 0.0.or. &
                  state(i,j,k,EN) <= 0.0) then
                ! Report to log file
                write(30,"(A,1PE10.3,2(2X,E10.3))") &
                     "Negative pressure/density/energy: ", &
                     pressr(i,j,k),state(i,j,k,RHO),state(i,j,k,EN)
                write(30,"(A,3(I4,X),A,1PE10.3)") " at ", &
                     i,j,k," time = ",time
                write(30,*) "call ",icall
                call flush(30)
                ! Set a control variable
                problem_counter=problem_counter+1

                if (problem_counter < 1001) then
                   problem_list(problem_counter)%position=(/ i, j, k /)
                   ! Pressure fix 1: diffuse with four neighbours
                   ! diffusion parameter is eta (used below)
                   ! To keep things conservative, express everything
                   ! as fluxes
                   imin=max(i-1,1) ! diffusive flux at edge 1 is zero
                   iplus=min(i+1,meshx) ! diffusive flux at edge 1 is zero
                   jmin=max(j-1,1)
                   jplus=min(j+1,meshy)
                   kmin=max(k-1,1)
                   kplus=min(k+1,meshz)
                   ! Only set diffuse flux if the neighbouring cell has
                   ! enough positive pressure itself, or if we are
                   ! correcting a negative density or energy.
                   if (-(pressr(imin,j,k)/pressr(i,j,k)) > eta .or. &
                        -(state(imin,j,k,RHO)/state(i,j,k,RHO)) > eta) then
                      !.or. &
                      !  pressr(i,j,k) > 0.0d0 ) then
                      problem_list(problem_counter)%dflux(1,1,:)= &
                           (state(imin,j,k,:)-state(i,j,k,:))
                   else
                      problem_list(problem_counter)%dflux(1,1,:)=0.0
                   endif
                   if (-(pressr(iplus,j,k)/pressr(i,j,k)) > eta .or. &
                        -(state(iplus,j,k,RHO)/state(i,j,k,RHO)) > eta) then
                      !.or. &
                      !  pressr(i,j,k) > 0.0d0 ) then
                      problem_list(problem_counter)%dflux(2,1,:)= &
                           (-state(iplus,j,k,:)+state(i,j,k,:))
                   else
                      problem_list(problem_counter)%dflux(2,1,:)=0.0
                   endif
                   if (-(pressr(i,jmin,k)/pressr(i,j,k)) > eta .or. &
                        -(state(i,jmin,k,RHO)/state(i,j,k,RHO)) > eta) then
                      !.or. &
                      !  pressr(i,j,k) > 0.0d0 ) then
                      problem_list(problem_counter)%dflux(1,2,:)= &
                           (state(i,jmin,k,:)-state(i,j,k,:))
                   else
                      problem_list(problem_counter)%dflux(1,2,:)=0.0
                   endif
                   if (-(pressr(i,jplus,k)/pressr(i,j,k)) > eta .or. &
                        -(state(i,jplus,k,RHO)/state(i,j,k,RHO)) > eta) then
                      !.or. &
                      !  pressr(i,j,k) > 0.0d0 ) then
                      problem_list(problem_counter)%dflux(2,2,:)= &
                           (-state(i,jplus,k,:)+state(i,j,k,:))
                   else
                      problem_list(problem_counter)%dflux(2,2,:)=0.0
                   endif
                   if (-(pressr(i,j,kmin)/pressr(i,j,k)) > eta .or.  &
                        -(state(i,j,kmin,RHO)/state(i,j,k,RHO)) > eta) then
                      !.or. &
                      !  pressr(i,j,k) > 0.0d0 ) then
                      problem_list(problem_counter)%dflux(1,3,:)= &
                           (state(i,j,kmin,:)-state(i,j,k,:))
                   else
                      problem_list(problem_counter)%dflux(1,3,:)=0.0
                   endif
                   if (-(pressr(i,j,kplus)/pressr(i,j,k)) > eta .or. &
                        -(state(i,j,kplus,RHO)/state(i,j,k,RHO)) > eta) then
                      !.or. &
                      !  pressr(i,j,k) > 0.0d0 ) then
                      problem_list(problem_counter)%dflux(2,3,:)= &
                           (-state(i,j,kplus,:)+state(i,j,k,:))
                   else
                      problem_list(problem_counter)%dflux(2,3,:)=0.0
                   endif
                else
                   write(30,*) "ERROR: Too many problems!"
                endif
             endif
          enddo
       enddo
    enddo
    
    do nproblem=1,problem_counter
       !write(30,*) 
       ! Apply fluxes
       i=problem_list(nproblem)%position(1)
       j=problem_list(nproblem)%position(2)
       k=problem_list(nproblem)%position(3)
       !write(30,"(A,5(1pe10.3))") "Old: ",state(i,j,k,1:neuler)
       !write(30,"(A,5(1pe10.3))") " ",state(i-1,j,k,1:neuler)
       !write(30,"(A,5(1pe10.3))") " ",state(i+1,j,k,1:neuler)
       !write(30,"(A,5(1pe10.3))") " ",state(i,j-1,k,1:neuler)
       !write(30,"(A,5(1pe10.3))") " ",state(i,j+1,k,1:neuler)
       !write(30,"(A,5(1pe10.3))") " ",state(i,j,k-1,1:neuler)
       !write(30,"(A,5(1pe10.3))") " ",state(i,j,k+1,1:neuler)
       state(i,j,k,:)=state(i,j,k,:)+eta* &
            (problem_list(nproblem)%dflux(1,1,:)- &
            problem_list(nproblem)%dflux(2,1,:) + &
            problem_list(nproblem)%dflux(1,2,:)- &
            problem_list(nproblem)%dflux(2,2,:) + &
            problem_list(nproblem)%dflux(1,3,:)- &
            problem_list(nproblem)%dflux(2,3,:))
       if (i > 1) state(i-1,j,k,:)=state(i-1,j,k,:)-eta* &
            problem_list(nproblem)%dflux(1,1,:)
       if (i < meshx) state(i+1,j,k,:)=state(i+1,j,k,:)+eta* &
            problem_list(nproblem)%dflux(2,1,:)
       if (j > 1) state(i,j-1,k,:)=state(i,j-1,k,:)-eta* &
            problem_list(nproblem)%dflux(1,2,:)
       if (j < meshy) state(i,j+1,k,:)=state(i,j+1,k,:)+eta* &
            problem_list(nproblem)%dflux(2,2,:)
       if (k > 1) state(i,j,k-1,:)=state(i,j,k-1,:)-eta* &
            problem_list(nproblem)%dflux(1,3,:)
       if (k < meshz) state(i,j,k+1,:)=state(i,j,k+1,:)+eta* &
            problem_list(nproblem)%dflux(2,3,:)
       !write(30,"(A,5(1pe10.3))") "New: ",state(i,j,k,1:neuler)
       !write(30,"(A,5(1pe10.3))") " ",state(i-1,j,k,1:neuler)
       !write(30,"(A,5(1pe10.3))") " ",state(i+1,j,k,1:neuler)
       !write(30,"(A,5(1pe10.3))") " ",state(i,j-1,k,1:neuler)
       !write(30,"(A,5(1pe10.3))") " ",state(i,j+1,k,1:neuler)
       !write(30,"(A,5(1pe10.3))") " ",state(i,j,k-1,1:neuler)
       !write(30,"(A,5(1pe10.3))") " ",state(i,j,k+1,1:neuler)
    enddo

    if (problem_counter > 0) then
       ! Recalculate the pressure after the diffusion
       ! the presfunc routine needs to be supplied
       call presfunc(sx,ex,sy,ey,sz,ez,newold,ierror)
       
       do k=sz,ez
          do j=sy,ey
             do i=sx,ex
                ! Check if the pressure is still negative
                if (pressr(i,j,k) <= 0.0) then ! check result of fix 1
                   write(30,*) "Still Negative pressure: ", &
                        pressr(i,j,k)," at ",i,j,k
                   call flush(30)
                   
                   ! Pressure fix 2: set temperature to minimum value
                   pnew=state(i,j,k,RHO)*boltzm*tmin/mu
                   !pnew=temper2pressr(tmin,rho2n(state(i,j,k,RHO)), &
                   !     electrondens(rho2n(state(i,j,k,RHO)), &
                   !     state(i,j,k,XHI:XHII))
                   state(i,j,k,EN)=state(i,j,k,EN)+(pnew-pressr(i,j,k))/gamma1
                   pressr(i,j,k)=pnew
                endif
                
                ! Check for negative densities and energies
                ! These are fatal. inegative is used to communicate
                ! this to the calling program.
                if (state(i,j,k,RHO) <= 0.0) then
                   write(30,*) "Still negative density: ", &
                        state(i,j,k,RHO)," at ",i,j,k
                   inegative=1
                endif
                if (state(i,j,k,EN) <= 0.0) then
                   write(30,*) "Still negative energy: ", &
                        state(i,j,k,EN)," at ",i,j,k
                   inegative=2
                endif
             enddo
          enddo
       enddo
    endif
    
  end subroutine presprot
  
end module protection
