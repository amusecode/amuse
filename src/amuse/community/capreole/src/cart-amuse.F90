module problem

  ! Module for Capreole3D (F90)
  ! Author: Garrelt Mellema
  ! Date: 2009-10-08

  ! This is the problem module. It contains the routines which define the
  ! problem being solved:

  ! Contents:
  ! init_problem - sets up the hydro variables according to the specified
  !                  problem (restart or fresh start)
  ! fresh_start_state - sets up the hydro variables for a fresh start
  ! inflow - sets the inflow boundary conditions (dummy here)
  ! apply_grav_force - apply gravity force (dummy here)

  ! This version: interface for AMUSE (cartesian coordinates). This is
  !  now a 3D simplified Riemann problem: density & pressure are set in
  !  two regions, an inner cube, half the size of the domain, and the
  !  region outside of that.

  use file_admin, only: stdinput, log_unit, file_input
  use precision, only: dp
  use cgsconstants, only: m_p
  use my_mpi
  use sizes, only: mbc,neq,RHO,RHVX,RHVY,RHVZ,EN,nrofDim
  use scaling, only: SCDENS, SCMOME, SCENER
  use atomic, only: gamma1
  use abundances, only: mu
  use mesh, only: sx,ex,sy,ey,sz,ez,meshx,meshy,meshz
  use grid, only: xlength,ylength,zlength,x,y,z
  use hydro, only:  state,pressr,set_state_pointer,NEW,OLD,restart_state,gforce
  use boundary, only: boundaries,REFLECTIVE,OUTFLOW,PROBLEM_DEF,X_IN,X_OUT,Y_IN, &
       Y_OUT,Z_IN,Z_OUT
  use ionic, only: init_ionic

  implicit none

  integer, dimension(nrofDim,2) :: domainboundaryconditions

  real(kind=dp) :: density_inner
  real(kind=dp) :: density_outer
  real(kind=dp) :: pressure_inner
  real(kind=dp) :: pressure_outer

contains

  subroutine init_problem (restart,restartfile)
    
    ! This routine initializes all hydro variables
    
    ! This may be a fresh start or a restart of a saved run
    
    !> tells you whether it's a new run or a restart
    logical,intent(in) :: restart 
    character(len=19),intent(in) :: restartfile !< file from which to restart

    ! Local variables
    real(kind=dp) :: r_interface !< dummy needed for calling init_ionic
    integer :: ierror !< error flag

    ! Set domain boundary conditions
    domainboundaryconditions(2:3,:)=REFLECTIVE    
    domainboundaryconditions(1,2)=REFLECTIVE
    domainboundaryconditions(1,1)=PROBLEM_DEF

    ! Fill the state variable
    if (.not.restart) then 

       ! Fresh start
       call fresh_start_state( )
       
    else

       ! Read state from restartfile
       call restart_state(restartfile,ierror)

       ! Scale to code scaling
       state(:,:,:,RHO)=state(:,:,:,RHO)/scdens
       state(:,:,:,RHVX)=state(:,:,:,RHVX)/scmome
       state(:,:,:,RHVY)=state(:,:,:,RHVY)/scmome
       state(:,:,:,RHVZ)=state(:,:,:,RHVZ)/scmome
       state(:,:,:,EN)=state(:,:,:,EN)/scener

       call boundaries(OLD,domainboundaryconditions,problemboundary) ! Fill boundary conditions

    endif
       
    ! Initialize the ionic concentrations
    call init_ionic(restart,r_interface)

  end subroutine init_problem

  !==========================================================================

  subroutine fresh_start_state ( )
    
    ! This routine initializes all hydro variables for a fresh start
    
    ! Case: amuse
    
    integer :: mx1,my1,mz1,temperature_setup
    character(len=512) :: densityfield_file
    real,allocatable,dimension(:,:,:) :: tempdens
    real(kind=dp) :: maxdens !< maximum density
    integer :: i,j,k,ieq

#ifdef MPI       
    integer :: status(MPI_STATUS_SIZE)
    integer :: request,nextproc
    integer,parameter :: outputcircle=601
    integer :: ierror
#endif

    ! Ask for the input if you are processor 0.
    if (rank == 0) then
       if (.not.file_input) write (*,"(A,$)") "1) Inner density: "
       read (stdinput,*) density_inner
       if (.not.file_input) write (*,"(A,$)") "2) Inner pressure: "
       read (stdinput,*) pressure_inner
       if (.not.file_input) write (*,"(A,$)") "3) Outer density: "
       read (stdinput,*) density_outer
       if (.not.file_input) write (*,"(A,$)") "4) Outer pressure: "
       read (stdinput,*) pressure_outer
    endif

    ! report input parameters
    if (rank == 0) then
       write(log_unit,"(A)") & 
            "Problem: cart-amuse (cartesian: AMUSE)"
       write (log_unit,"(A,1PE10.3)") "1) Inner density: ",density_inner
       write (log_unit,"(A,1PE10.3)") "2) Inner pressure: ",pressure_inner
       write (log_unit,"(A,1PE10.3)") "3) Outer density: ",density_outer
       write (log_unit,"(A,1PE10.3)") "4) Outer pressure: ",pressure_outer
    endif
    
#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    ! (densityfield_file is distributed in the MPI ring below)
    call MPI_BCAST(density_inner,1,MPI_DOUBLE_PRECISION,0, & 
         MPI_COMM_NEW,ierror)
    call MPI_BCAST(pressure_inner,1,MPI_DOUBLE_PRECISION,0, & 
         MPI_COMM_NEW,ierror)
    call MPI_BCAST(density_outer,1,MPI_DOUBLE_PRECISION,0, & 
         MPI_COMM_NEW,ierror)
    call MPI_BCAST(pressure_outer,1,MPI_DOUBLE_PRECISION,0, & 
         MPI_COMM_NEW,ierror)
#endif
    
    ! Set density in state array and pressure
    do k=sz-mbc,ez+mbc
       do j=sy-mbc,ey+mbc
          do i=sx-mbc,ex+mbc
             ! inner densities in the inner cube
             if (abs(x(i)-0.5*xlength) < 0.25*xlength .and. &
                  abs(y(j)-0.5*ylength) < 0.25*ylength .and. &
                  abs(z(k)-0.5*zlength) < 0.25*zlength) then
                state(i,j,k,RHO)=density_inner
                pressr(i,j,k)=pressure_inner
             else
                state(i,j,k,RHO)=density_outer
                pressr(i,j,k)=pressure_outer
             endif
          enddo
       enddo
    enddo
    
    ! Set the initial momenta to zero and calculate the energy density
    do k=sz-mbc,ez+mbc
       do j=sy-mbc,ey+mbc
          do i=sx-mbc,ex+mbc
             state(i,j,k,RHVX)=0.0d0
             state(i,j,k,RHVY)=0.0d0
             state(i,j,k,RHVZ)=0.0d0
             state(i,j,k,EN)=pressr(i,j,k)/gamma1
          enddo
       enddo
    enddo

  end subroutine fresh_start_state
  
  !==========================================================================

  subroutine problemboundary (boundary_id,newold)
    
    ! This routine resets the inner boundary to the inflow condition

    ! Version: dummy routine

    integer,intent(in) :: boundary_id
    integer,intent(in) :: newold
    
    integer :: i,j,k

    ! Point state to appropriate array
    state => set_state_pointer(newold)

    select case (boundary_id)
    case (X_IN)
       do k=sz-1,ez+1
          do j=sy-1,ey+1
             do i=1-mbc,1
                state(i,j,k,RHO)=density_outer
                state(i,j,k,RHVX)=0.0
                state(i,j,k,RHVY)=0.0
                state(i,j,k,RHVZ)=0.0
                pressr(i,j,k)=pressure_outer
                state(i,j,k,EN)=pressr(i,j,k)/gamma1
             enddo
          enddo
       enddo
    case (X_OUT)
    case (Y_IN)
    case (Y_OUT)
    case (Z_IN)
    case (Z_OUT)
    end select
    
  end subroutine problemboundary
  
  !==========================================================================


  subroutine apply_grav_force(dt,newold)
                
    ! Dummy routine
             
    real(kind=dp),intent(in) :: dt
    integer,intent(in) :: newold
    
  end subroutine apply_grav_force
    

end module problem
