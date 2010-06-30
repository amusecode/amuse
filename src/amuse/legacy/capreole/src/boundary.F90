module boundary

  ! Module for Capreole
  ! Author: Garrelt Mellema
  ! Date: 2010-06-14
  !
  ! This module contains the routines related with handling the grid
  ! boundaries. This is mostly the internal (MPI) boundaries, but also
  ! the outside boundaries.
  !
  ! History:
  ! 2004-05-11: first dated version
  ! 2007-10-05: clean-up, added only's to use statements.
  ! 2010-01-25: added option for problem specific boundary conditions here
  ! 2010-06-14: added periodic boundary conditions.

  use precision, only: dp
  use my_mpi
  use sizes, only: neq, mbc, RHO, RHVX, RHVY, RHVZ, EN, nrofDim
  use mesh, only: sx,ex,sy,ey,sz,ez,meshx,meshy,meshz
  use grid, only: vol,volx
  use geometry, only: presfunc
  use hydro, only: state, NEW, OLD, set_state_pointer

  implicit none
  private

  ! Define outflow types
  integer,parameter,public :: REFLECTIVE=-1
  integer,parameter,public :: REFLECTIVE_SHIFT=-2
  integer,parameter,public :: OUTFLOW=1
  integer,parameter,public :: PERIODIC=2
  integer,parameter,public :: PROBLEM_DEF=3
  integer,parameter,public :: X_IN=1
  integer,parameter,public :: X_OUT=2
  integer,parameter,public :: Y_IN=3
  integer,parameter,public :: Y_OUT=4
  integer,parameter,public :: Z_IN=5
  integer,parameter,public :: Z_OUT=6

  public :: boundaries

contains

  !=======================================================================

  subroutine boundaries (newold, &
       domainboundaryconditions, &
       problem_boundary_routine)
    
    ! Deals with the boundary conditions: internal boundaries (exchanges 
    ! boundary cells between neighbours) and external boundaries.
    ! extend of boundary: mbc

    integer,intent(in) :: newold
    integer,dimension(nrOfDim,2),intent(in) :: domainboundaryconditions
    interface
       subroutine problem_boundary_routine (boundary_id, state_id)
         integer,intent(in) :: boundary_id
         integer,intent(in) :: state_id
       end subroutine problem_boundary_routine
    end interface
  
#ifdef MPI
    real(kind=dp),dimension(mbc,1-mbc:ey-sy+1+mbc,1-mbc:ez-sz+1+mbc,neq) :: &
         xplane1
    real(kind=dp),dimension(mbc,1-mbc:ey-sy+1+mbc,1-mbc:ez-sz+1+mbc,neq) :: &
         xplane2
    real(kind=dp),dimension(mbc,1-mbc:ey-sy+1+mbc,1-mbc:ez-sz+1+mbc,neq) :: &
         xplane3
    real(kind=dp),dimension(mbc,1-mbc:ey-sy+1+mbc,1-mbc:ez-sz+1+mbc,neq) :: &
         xplane4

    real(kind=dp),dimension(1-mbc:ex-sx+1+mbc,mbc,1-mbc:ez-sz+1+mbc,neq) :: &
         yplane1
    real(kind=dp),dimension(1-mbc:ex-sx+1+mbc,mbc,1-mbc:ez-sz+1+mbc,neq) :: &
         yplane2
    real(kind=dp),dimension(1-mbc:ex-sx+1+mbc,mbc,1-mbc:ez-sz+1+mbc,neq) :: &
         yplane3
    real(kind=dp),dimension(1-mbc:ex-sx+1+mbc,mbc,1-mbc:ez-sz+1+mbc,neq) :: &
         yplane4

    real(kind=dp),dimension(1-mbc:ex-sx+1+mbc,1-mbc:ey-sy+1+mbc,mbc,neq) :: &
         zplane1
    real(kind=dp),dimension(1-mbc:ex-sx+1+mbc,1-mbc:ey-sy+1+mbc,mbc,neq) :: &
         zplane2
    real(kind=dp),dimension(1-mbc:ex-sx+1+mbc,1-mbc:ey-sy+1+mbc,mbc,neq) :: &
         zplane3
    real(kind=dp),dimension(1-mbc:ex-sx+1+mbc,1-mbc:ey-sy+1+mbc,mbc,neq) :: &
         zplane4

    integer :: status(MPI_STATUS_SIZE)
    integer :: i,j,k,ieq
    
    integer :: sizex,sizey,sizez

    integer,parameter :: xfromright=101,xfromleft=102  ! tags
    integer ::  request1,request2
    integer,parameter ::  yfromdown=201,yfromup=202  ! tags
    integer ::  request3,request4
    integer,parameter ::  zfrombelow=301,zfromabove=302  ! tags
    integer ::  request5,request6
#endif

    integer :: ierror
 
    ! Point state to appropriate array
    state => set_state_pointer(newold)

    ! Account for non-existing neighbours, these are real boundaries
    ! the [inner,outer][x,y]bound routines need to be supplied

    if (nbrleft == MPI_PROC_NULL) &
         call innerxbound(newold,domainboundaryconditions(1,1),problem_boundary_routine)
    if (nbrright == MPI_PROC_NULL) &
         call outerxbound(newold,domainboundaryconditions(1,2),problem_boundary_routine)
    if (nbrdown == MPI_PROC_NULL) &
         call innerybound(newold,domainboundaryconditions(2,1),problem_boundary_routine)
    if (nbrup == MPI_PROC_NULL)   &
         call outerybound(newold,domainboundaryconditions(2,2),problem_boundary_routine)
    if (nbrbelow == MPI_PROC_NULL) &
         call innerzbound(newold,domainboundaryconditions(3,1),problem_boundary_routine)
    if (nbrabove == MPI_PROC_NULL) &
         call outerzbound(newold,domainboundaryconditions(3,2),problem_boundary_routine)

#ifdef MPI
    ! Sizes in x and y direction

    sizex=ex-sx+1
    sizey=ey-sy+1
    sizez=ez-sz+1
    
    ! Exchange mbc wide yrows with left and right neighbours

    do ieq=1,neq   ! put planes to be sent in array
       do k=sz-mbc,ez+mbc
          do j=sy-mbc,ey+mbc
             do i=sx,sx+mbc-1
                xplane1(i-sx+1,j-sy+1,k-sz+1,ieq)=state(i,j,k,ieq)
             enddo
             do i=ex-mbc+1,ex
                xplane2(i-ex+mbc,j-sy+1,k-sz+1,ieq)=state(i,j,k,ieq)
             enddo
          enddo
       enddo
    enddo

    ! SendReceive planes to left and right neighbours 

    call MPI_SENDRECV( &
         xplane1,mbc*(sizey+2*mbc)*(sizez+2*mbc)*neq, &
         MPI_DOUBLE_PRECISION,nbrleft,xfromright, &
         xplane4,mbc*(sizey+2*mbc)*(sizez+2*mbc)*neq, &
         MPI_DOUBLE_PRECISION,nbrright,xfromright, &
         MPI_COMM_NEW,status,ierror)

    call MPI_SENDRECV( &
         xplane2,mbc*(sizey+2*mbc)*(sizez+2*mbc)*neq, &
         MPI_DOUBLE_PRECISION,nbrright,xfromleft, &
         xplane3,mbc*(sizey+2*mbc)*(sizez+2*mbc)*neq, & 
         MPI_DOUBLE_PRECISION,nbrleft,xfromleft, &
         MPI_COMM_NEW,status,ierror)
      
    ! Fill the boundaries with the received planes

    if (nbrleft /= MPI_PROC_NULL) then
       do ieq=1,neq
          do k=sz-mbc,ez+mbc
             do j=sy-mbc,ey+mbc
                do i=sx-mbc,sx-1
                   state(i,j,k,ieq)=xplane3(i-sx+mbc+1,j-sy+1,k-sz+1,ieq)
                enddo
             enddo
          enddo
       enddo
    endif
    if (nbrright /= MPI_PROC_NULL) then
       do ieq=1,neq
          do k=sz-mbc,ez+mbc
             do j=sy-mbc,ey+mbc
                do i=ex+1,ex+mbc
                   state(i,j,k,ieq)=xplane4(i-ex,j-sy+1,k-sz+1,ieq)
                enddo
             enddo
          enddo
       enddo
    endif
    
    ! Wait for the sends to be completed
    call MPI_BARRIER(MPI_COMM_NEW,ierror)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! Y 
    ! Exchange mbc thick y-planes with the neighbours

    do ieq=1,neq    ! put planes to be sent in array
       do k=sz-mbc,ez+mbc
          do i=sx-mbc,ex+mbc
             do j=sy,sy+mbc-1
                yplane1(i-sx+1,j-sy+1,k-sz+1,ieq)=state(i,j,k,ieq)
             enddo
             do j=ey-mbc+1,ey
                yplane2(i-sx+1,j-ey+mbc,k-sz+1,ieq)=state(i,j,k,ieq)
             enddo
          enddo
       enddo
    enddo
    
    ! SendReceive planes to down and up neighbours

    call MPI_SENDRECV( &
         yplane1,(sizex+2*mbc)*(sizez+2*mbc)*mbc*neq, &
         MPI_DOUBLE_PRECISION,nbrdown,yfromup,  &
         yplane4,(sizex+2*mbc)*(sizez+2*mbc)*mbc*neq, &
         MPI_DOUBLE_PRECISION,nbrup,yfromup, &
         MPI_COMM_NEW,status,ierror)
          
    call MPI_SENDRECV( &
         yplane2,(sizex+2*mbc)*(sizez+2*mbc)*mbc*neq, &
         MPI_DOUBLE_PRECISION,nbrup,yfromdown,   &
         yplane3,(sizex+2*mbc)*(sizez+2*mbc)*mbc*neq, &
         MPI_DOUBLE_PRECISION,nbrdown,yfromdown, &
         MPI_COMM_NEW,status,ierror)
      
    ! Fill the boundaries with the received planes

    if (nbrdown /= MPI_PROC_NULL) then
       do ieq=1,neq
          do k=sz-mbc,ez+mbc
             do j=sy-mbc,sy-1
                do i=sx-mbc,ex+mbc
                   state(i,j,k,ieq)=yplane3(i-sx+1,j-sy+mbc+1,k-sz+1,ieq)
                enddo
             enddo
          enddo
       enddo
    endif
    if (nbrup /= MPI_PROC_NULL) then
       do ieq=1,neq
          do k=sz-mbc,ez+mbc
             do j=ey+1,ey+mbc
                do i=sx-mbc,ex+mbc
                   state(i,j,k,ieq)=yplane4(i-sx+1,j-ey,k-sz+1,ieq)
                enddo
             enddo
          enddo
       enddo
    endif
    
    ! Wait for the sends to be completed
    
    call MPI_BARRIER(MPI_COMM_NEW,ierror)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! Z
    ! Exchange z-planes with the neighbours

    do ieq=1,neq    ! put planes to be sent in array
       do j=sy-mbc,ey+mbc
          do i=sx-mbc,ex+mbc
             do k=sz,sz+mbc-1
                zplane1(i-sx+1,j-sy+1,k-sz+1,ieq)=state(i,j,k,ieq)
             enddo
             do k=ez-mbc+1,ez
                zplane2(i-sx+1,j-sy+1,k-ez+mbc,ieq)=state(i,j,k,ieq)
             enddo
          enddo
       enddo
    enddo
    
    ! SendReceive planes to down and above neighbours

    call MPI_SENDRECV( &
         zplane1,(sizex+2*mbc)*(sizey+2*mbc)*mbc*neq, &
         MPI_DOUBLE_PRECISION,nbrbelow,zfromabove,  &
         zplane3,(sizex+2*mbc)*(sizey+2*mbc)*mbc*neq, &
         MPI_DOUBLE_PRECISION,nbrabove,zfromabove, &
         MPI_COMM_NEW,status,ierror)

    call MPI_SENDRECV( &
         zplane2,(sizex+2*mbc)*(sizey+2*mbc)*mbc*neq, &
         MPI_DOUBLE_PRECISION,nbrabove,zfrombelow,   &
         zplane4,(sizex+2*mbc)*(sizey+2*mbc)*mbc*neq, &
         MPI_DOUBLE_PRECISION,nbrbelow,zfrombelow, &
         MPI_COMM_NEW,status,ierror)

    ! Fill the boundaries with the received planes

    if (nbrbelow /= MPI_PROC_NULL) then
       do ieq=1,neq
          do k=sz-mbc,sz-1
             do j=sy-mbc,ey+mbc
                do i=sx-mbc,ex+mbc
                   state(i,j,k,ieq)=zplane3(i-sx+1,j-sy+1,k-sz+mbc+1,ieq)
                enddo
             enddo
          enddo
       enddo
    endif
    if (nbrabove /= MPI_PROC_NULL) then
       do ieq=1,neq
          do k=ez+1,ez+mbc
             do j=sy-mbc,ey+mbc
                do i=sx-mbc,ex+mbc
                   state(i,j,k,ieq)=zplane4(i-sx+1,j-sy+1,k-ez,ieq)
                enddo
             enddo
          enddo
       enddo
    endif
    
    ! Wait for the sends to be completed
    
    call MPI_BARRIER(MPI_COMM_NEW,ierror)
#endif

    ! Calculate the pressure in the newly set boundary cells
    ! the presfunc routine has to be supplied

    call presfunc(sx-mbc,ex+mbc,sy-mbc,ey+mbc,sz-mbc,ez+mbc,newold,ierror)
    !! call presfunc(ex+1,ex+mbc,sy-1,ey+1,0)
    !! call presfunc(sx-1,ex+1,sy-mbc,sy-1,0)
    !! call presfunc(sx-1,ex+1,ey+1,ey+mbc,0)

  end subroutine boundaries

  !==========================================================================

  subroutine innerxbound (newold, boundarycondition, problem_boundary)
    
    ! This routine resets the inner x boundary
    
    integer :: i,j,k,ieq
    
    integer,intent(in) :: newold
    integer,intent(in) :: boundarycondition

    interface
       subroutine problem_boundary (boundary_id, state_id)
         integer,intent(in) :: boundary_id
         integer,intent(in) :: state_id
       end subroutine problem_boundary
    end interface
  
    state => set_state_pointer(newold)

    select case (boundarycondition)
    case (PERIODIC)
       if (sx /= 1 .or. ex /= meshx) then
          write(*,*) "Error, applying periodic boundary conditions in MPI case"
          write(*,*) "Fix mpi.F90"
       endif
       do ieq=1,neq
          do k=sz-mbc,ez+mbc
             do j=sy-mbc,ey+mbc
                do i=sx-mbc,sx-1
                   state(i,j,k,ieq)=state(ex+i,j,k,ieq)
                enddo
             enddo
          enddo
       enddo
    case (REFLECTIVE)
       do k=sz-mbc,ez+mbc
          do j=sy-mbc,ey+mbc
             do i=sx-mbc,sx-1
                state(i,j,k,RHO)=state(2*sx-1-i,j,k,RHO)
                state(i,j,k,RHVX)=-state(2*sx-1-i,j,k,RHVX)
                state(i,j,k,RHVY:neq)=state(2*sx-1-i,j,k,RHVY:neq)
             enddo
          enddo
       enddo
    case (REFLECTIVE_SHIFT)
       do k=sz-mbc,ez+mbc
          do j=sy-mbc,ey+mbc
             do i=sx-mbc,sx-1
                state(i,j,k,RHO)=state(2*sx-i,j,k,RHO)
                state(i,j,k,RHVX)=-state(2*sx-i,j,k,RHVX)
                state(i,j,k,RHVY:neq)=state(2*sx-i,j,k,RHVY:neq)
             enddo
          enddo
       enddo
    case(OUTFLOW)
       do ieq=1,neq
          do k=sz-mbc,ez+mbc
             do j=sy-mbc,ey+mbc
                do i=sx-mbc,sx-1
                   state(i,j,k,ieq)=state(sx,j,k,ieq)
                enddo
             enddo
          enddo
       enddo
    case(PROBLEM_DEF)
       call problem_boundary(X_IN,newold)
    end select

  end subroutine innerxbound

  !==========================================================================

  subroutine outerxbound (newold, boundarycondition, problem_boundary)

    ! This routine resets the outer x boundary

    integer :: i,j,k,ieq
    
    integer,intent(in) :: newold
    integer,intent(in) :: boundarycondition

    interface
       subroutine problem_boundary (boundary_id, state_id)
         integer,intent(in) :: boundary_id
         integer,intent(in) :: state_id
       end subroutine problem_boundary
    end interface
  
    state => set_state_pointer(newold)

    select case (boundarycondition)
    case (PERIODIC)
       if (sx /= 1 .or. ex /= meshx) then
          write(*,*) "Error, applying periodic boundary conditions in MPI case"
          write(*,*) "Fix mpi.F90"
       endif
       do ieq=1,neq
          do k=sz-mbc,ez+mbc
             do j=sy-mbc,ey+mbc
                do i=ex+1,ex+mbc
                   state(i,j,k,ieq)=state(i-ex,j,k,ieq)
                enddo
             enddo
          enddo
       enddo
    case (REFLECTIVE)
       do k=sz-mbc,ez+mbc
          do j=sy-mbc,ey+mbc
             do i=ex+1,ex+mbc
                state(i,j,k,RHO)=state(2*ex+1-i,j,k,RHO)
                state(i,j,k,RHVX)=-state(2*ex+1-i,j,k,RHVX)
                state(i,j,k,RHVY:neq)=state(2*ex+1-i,j,k,RHVY:neq)
             enddo
          enddo
       enddo
    case (REFLECTIVE_SHIFT)
       do k=sz-mbc,ez+mbc
          do j=sy-mbc,ey+mbc
             do i=ex+1,ex+mbc
                state(i,j,k,RHO)=state(2*ex-i,j,k,RHO)
                state(i,j,k,RHVX)=-state(2*ex-i,j,k,RHVX)
                state(i,j,k,RHVY:neq)=state(2*ex-i,j,k,RHVY:neq)
             enddo
          enddo
       enddo
    case(OUTFLOW)
       do ieq=1,neq
          do k=sz-mbc,ez+mbc
             do j=sy-mbc,ey+mbc
                do i=ex+1,ex+mbc
                   state(i,j,k,ieq)=state(ex,j,k,ieq)
                enddo
             enddo
          enddo
       enddo
    case(PROBLEM_DEF)
       call problem_boundary(X_OUT,newold)
    end select

  end subroutine outerxbound

  !==========================================================================

  subroutine innerybound (newold, boundarycondition, problem_boundary)

    ! This routine resets the inner y boundary

    integer :: i,j,k,ieq

    integer,intent(in) :: newold
    integer,intent(in) :: boundarycondition

    interface
       subroutine problem_boundary (boundary_id, state_id)
         integer,intent(in) :: boundary_id
         integer,intent(in) :: state_id
       end subroutine problem_boundary
    end interface
  
    state => set_state_pointer(newold)

    select case (boundarycondition)
    case (PERIODIC)
       if (sy /= 1 .or. ey /= meshy) then
          write(*,*) "Error, applying periodic boundary conditions in MPI case"
          write(*,*) "Fix mpi.F90"
       endif
       do ieq=1,neq
          do k=sz-mbc,ez+mbc
             do j=sy-mbc,sy-1
                do i=sx-mbc,ex+mbc
                   state(i,j,k,ieq)=state(i,ey+j,k,ieq)
                enddo
             enddo
          enddo
       enddo
    case (REFLECTIVE)
       do k=sz-mbc,ez+mbc
          do j=sy-mbc,sy-1
             do i=sx-mbc,ex+mbc
                state(i,j,k,RHO)=state(i,2*sy-1-j,k,RHO)
                state(i,j,k,RHVX)=state(i,2*sy-1-j,k,RHVX)
                state(i,j,k,RHVY)=-state(i,2*sy-1-j,k,RHVY)
                state(i,j,k,RHVZ:neq)=state(i,2*sy-1-j,k,RHVZ:neq)
             enddo
          enddo
       enddo
    case (REFLECTIVE_SHIFT)
       do k=sz-mbc,ez+mbc
          do j=sy-mbc,sy-1
             do i=sx-mbc,ex+mbc
                state(i,j,k,RHO)=state(i,2*sy-j,k,RHO)
                state(i,j,k,RHVX)=state(i,2*sy-j,k,RHVX)
                state(i,j,k,RHVY)=-state(i,2*sy-j,k,RHVY)
                state(i,j,k,RHVZ:neq)=state(i,2*sy-j,k,RHVZ:neq)
             enddo
          enddo
       enddo
    case(OUTFLOW)
       do ieq=1,neq
          do k=sz-mbc,ez+mbc
             do j=sy-mbc,sy-1
                do i=sx-mbc,ex+mbc
                   state(i,j,k,ieq)=state(i,sy,k,ieq)
                enddo
             enddo
          enddo
       enddo
    case(PROBLEM_DEF)
       call problem_boundary(Y_IN,newold)
    end select

  end subroutine innerybound

  !==========================================================================

  subroutine outerybound (newold, boundarycondition, problem_boundary)

    ! This routine resets the outer y boundary

    integer :: i,j,k,ieq

    integer,intent(in) :: newold
    integer,intent(in) :: boundarycondition

    interface
       subroutine problem_boundary (boundary_id, state_id)
         integer,intent(in) :: boundary_id
         integer,intent(in) :: state_id
       end subroutine problem_boundary
    end interface
  
    state => set_state_pointer(newold)

    select case (boundarycondition)
    case (PERIODIC)
       if (sy /= 1 .or. ey /= meshy) then
          write(*,*) "Error, applying periodic boundary conditions in MPI case"
          write(*,*) "Fix mpi.F90"
       endif
       do ieq=1,neq
          do k=sz-mbc,ez+mbc
             do j=ey+1,ey+mbc
                do i=sx-mbc,ex+mbc
                   state(i,j,k,ieq)=state(i,j-ey,k,ieq)
                enddo
             enddo
          enddo
       enddo
    case (REFLECTIVE)
       do k=sz-mbc,ez+mbc
          do j=ey+1,ey+mbc
             do i=sx-mbc,ex+mbc
                state(i,j,k,RHO)=state(i,2*ey+1-j,k,RHO)
                state(i,j,k,RHVX)=state(i,2*ey+1-j,k,RHVX)
                state(i,j,k,RHVY)=-state(i,2*ey+1-j,k,RHVY)
                state(i,j,k,RHVZ:neq)=state(i,2*ey+1-j,k,RHVZ:neq)
             enddo
          enddo
       enddo
    case (REFLECTIVE_SHIFT)
       do k=sz-mbc,ez+mbc
          do j=ey+1,ey+mbc
             do i=sx-mbc,ex+mbc
                state(i,j,k,RHO)=state(i,2*ey-j,k,RHO)
                state(i,j,k,RHVX)=state(i,2*ey-j,k,RHVX)
                state(i,j,k,RHVY)=-state(i,2*ey-j,k,RHVY)
                state(i,j,k,RHVZ:neq)=state(i,2*ey-j,k,RHVZ:neq)
             enddo
          enddo
       enddo
    case(OUTFLOW)
       do ieq=1,neq
          do k=sz-mbc,ez+mbc
             do j=ey+1,ey+mbc
                do i=sx-mbc,ex+mbc
                   state(i,j,k,ieq)=state(i,ey,k,ieq)
                enddo
             enddo
          enddo
       enddo
    case(PROBLEM_DEF)
       call problem_boundary(Y_OUT,newold)
    end select
    
  end subroutine outerybound

  !==========================================================================

  subroutine innerzbound (newold, boundarycondition, problem_boundary)

    ! This routine resets the inner z boundary

    integer,intent(in) :: newold
    integer,intent(in) :: boundarycondition
    
    interface
       subroutine problem_boundary (boundary_id, state_id)
         integer,intent(in) :: boundary_id
         integer,intent(in) :: state_id
       end subroutine problem_boundary
    end interface
  
    integer :: i,j,k,ieq

    state => set_state_pointer(newold)

    select case (boundarycondition)
    case (PERIODIC)
       if (sz /= 1 .or. ez /= meshz) then
          write(*,*) "Error, applying periodic boundary conditions in MPI case"
          write(*,*) "Fix mpi.F90"
       endif
       do ieq=1,neq
          do k=sz-mbc,sz-1
             do j=sy-mbc,ey+mbc
                do i=sx-mbc,ex+mbc
                   state(i,j,k,ieq)=state(i,j,ez+k,ieq)
                enddo
             enddo
          enddo
       enddo
    case (REFLECTIVE)
       do k=sz-mbc,sz-1
          do j=sy-mbc,ey+mbc
             do i=sx-mbc,ex+mbc
                state(i,j,k,RHO)=state(i,j,2*sz-1-k,RHO)
                state(i,j,k,RHVX)=state(i,j,2*sz-1-k,RHVX)
                state(i,j,k,RHVY)=state(i,j,2*sz-1-k,RHVY)
                state(i,j,k,RHVZ)=-state(i,j,2*sz-1-k,RHVZ)
                state(i,j,k,EN:neq)=state(i,j,2*sz-1-k,EN:neq)
             enddo
          enddo
       enddo
    case (REFLECTIVE_SHIFT)
       do k=sz-mbc,sz-1
          do j=sy-mbc,ey+mbc
             do i=sx-mbc,ex+mbc
                state(i,j,k,RHO)=state(i,j,2*sz-k,RHO)
                state(i,j,k,RHVX)=state(i,j,2*sz-k,RHVX)
                state(i,j,k,RHVY)=state(i,j,2*sz-k,RHVY)
                state(i,j,k,RHVZ)=-state(i,j,2*sz-k,RHVZ)
                state(i,j,k,EN:neq)=state(i,j,2*sz-k,EN:neq)
             enddo
          enddo
       enddo
    case(OUTFLOW)
       do ieq=1,neq
          do k=sz-mbc,sz-1
             do j=sy-mbc,ey+mbc
                do i=sx-mbc,ex+mbc
                   state(i,j,k,ieq)=state(i,j,sz,ieq)
                enddo
             enddo
          enddo
       enddo
    case(PROBLEM_DEF)
       call problem_boundary(Z_IN,newold)
    end select

  end subroutine innerzbound

  !==========================================================================

  subroutine outerzbound (newold, boundarycondition, problem_boundary)

    ! This routine resets the outer y boundary

    integer,intent(in) :: newold
    integer,intent(in) :: boundarycondition

    interface
       subroutine problem_boundary (boundary_id, state_id)
         integer,intent(in) :: boundary_id
         integer,intent(in) :: state_id
       end subroutine problem_boundary
    end interface
  
    integer :: i,j,k,ieq

    state => set_state_pointer(newold)

    select case (boundarycondition)
    case (PERIODIC)
       if (sz /= 1 .or. ez /= meshx) then
          write(*,*) "Error, applying periodic boundary conditions in MPI case"
          write(*,*) "Fix mpi.F90"
       endif
       do ieq=1,neq
          do k=ez+1,ez+mbc
             do j=sy-mbc,ey+mbc
                do i=sx-mbc,ex+mbc
                   state(i,j,k,ieq)=state(i,j,k-ez,ieq)
                enddo
             enddo
          enddo
       enddo
    case (REFLECTIVE)
       do k=ez+1,ez+mbc
          do j=sy-mbc,ey+mbc
             do i=sx-mbc,ex+mbc
                state(i,j,k,RHO)=state(i,j,2*ez+1-k,RHO)
                state(i,j,k,RHVX)=state(i,j,2*ez+1-k,RHVX)
                state(i,j,k,RHVY)=state(i,j,2*ez+1-k,RHVY)
                state(i,j,k,RHVZ)=-state(i,j,2*ez+1-k,RHVZ)
                state(i,j,k,EN:neq)=state(i,j,2*ez+1-k,EN:neq)
             enddo
          enddo
       enddo
    case (REFLECTIVE_SHIFT)
       do k=ez+1,ez+mbc
          do j=sy-mbc,ey+mbc
             do i=sx-mbc,ex+mbc
                state(i,j,k,RHO)=state(i,j,2*ez-k,RHO)
                state(i,j,k,RHVX)=state(i,j,2*ez-k,RHVX)
                state(i,j,k,RHVY)=state(i,j,2*ez-k,RHVY)
                state(i,j,k,RHVZ)=-state(i,j,2*ez-k,RHVZ)
                state(i,j,k,EN:neq)=state(i,j,2*ez-k,EN:neq)
             enddo
          enddo
       enddo
    case(OUTFLOW)
       do ieq=1,neq
          do k=ez+1,ez+mbc
             do j=sy-mbc,ey+mbc
                do i=sx-mbc,ex+mbc
                   state(i,j,k,ieq)=state(i,j,ez,ieq)
                enddo
             enddo
          enddo
       enddo
    case(PROBLEM_DEF)
       call problem_boundary(Z_OUT,newold)
    end select

  end subroutine outerzbound


end module boundary
