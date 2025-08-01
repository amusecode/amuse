module output
  
  ! Module for Capreole (3D)
  ! Author: Garrelt Mellema
  ! Date: 2007-10-05 (previous: unknown)
  !
  ! This module contains the routines related with handling the data
  ! output
  !
  ! History:
  ! 2007-10-05: clean-up, added only's to use statements, added this
  !             introduction.
  !
  ! Version: Output using the ah3 format.
  !
  ! Contents:
  ! - init_output: prepares the generics for the output
  ! - make_output: does the output of one snap shot

  use file_admin, only: stdinput,ah3,log_unit,file_input
  use precision, only: dp
  use scaling, only: SCDENS,SCMOME,SCENER,SCTIME,SCLENG
  use sizes, only: nrOfDim,neq,neuler,RHO,RHVX,RHVY,RHVZ,EN
  use my_mpi
  use mesh, only: sx,ex,sy,ey,sz,ez
  use grid, only: x,y,z,dx,dy,dz
  use atomic, only: gamma
  use hydro, only: state
  use times, only: time

  implicit none

  private

  character(len=1),private :: runid
  character(len=8),private :: revdate

  public :: init_output, make_output

#ifdef MPI
    integer :: mpi_ierror
#endif

contains
  
  !========================================================================
  subroutine init_output(restart,restartfile,ierror)
    
    ! This routine prepares the data output stream.
    
    logical,intent(in) :: restart
    character(len=19),intent(in) :: restartfile
    integer,intent(out) :: ierror

    character(len=19) :: filename ! name of output file
    character(len=8) :: date
    logical :: test_exist

    ! Run ID
    if (rank == 0) then
       if (.not.restart) then
          ! Construct date-identifier: d-m-y
          call date_and_time(date)
          revdate(1:2)=date(7:8)
          revdate(3:4)=date(5:6)
          revdate(5:8)=date(1:4)
          if (.not.file_input) write(*,"(A,$)") "Run ID (one letter): "
          read(stdinput,*) runid
          ! Inquire if it exists
          do
             filename=revdate//"_"//runid//"_0000.ah3"
             inquire(file=filename,exist=test_exist)
             if (.not.test_exist) exit
             if (runid < "z") then
                runid=achar(iachar(runid)+1)
                write(log_unit,*) "RunID already in use, trying ",runid
             else
                write(log_unit,*) "All runids are in use?"
                ierror=1
             endif
          enddo
       else
          revdate=restartfile(1:8)
          runid=restartfile(10:10)
       endif
    endif

#ifdef MPI
    ! Distribute runid over nodes
    call MPI_BCAST(runid,1,MPI_CHARACTER,0,MPI_COMM_NEW,mpi_ierror)
#endif     

  end subroutine init_output

  !========================================================================
  subroutine make_output(nframe)
    
    ! This routine takes care of the output
    ! It deals with a parallel run by using an MPI-ring
    ! structure.
    
    integer,intent(in) :: nframe
    
    ! AH3D Output variables
    character(len=80),parameter :: banner="Capreole (F90) 3D Hydrodynamics"
    integer,parameter :: refinementFactor=1
    integer,parameter :: level=1

    integer   :: ierror,i,j,k,ieq ! counters and error flags
    character(len=19) :: filename ! name of output file
    character(len=4)  :: string_frame

#ifdef MPI
    integer :: status(MPI_STATUS_SIZE)
    integer :: request,nextproc
    integer,parameter :: outputcircle=601
#endif

    ierror=0
    write(log_unit,*) "Writing Frame: ",nframe
      
    ! Construct file name and open the file
    write(string_frame,"(i4.4)") nframe
    filename=revdate//"_"//runid//"_"//string_frame//".ah3"
    
    ! AH3D output format (24-10-2002)

    ! File name: 24102002_a_0000.ah3

    ! Header:
    ! string 80 bytes
    ! int    nrOfDim
    ! int    nrOfVars
    ! int    nrOfGrids
    ! int    refinementFactor
    ! int    frameNr
    ! double gamma
    ! double time

    ! Grids:
    ! int    nrOfCells1 [, nrOfCells2, nrOfCells3]
    ! double corner1 [, corner2, corner3]
    ! double cellSize1 [,cellSize2, cellSize3]
    ! int    level

    ! Cells;
    ! double rho
    ! double rho*v1 [, rho*v1, rho*v3]
    ! double rho*e
    ! [double var(nrOfVars-(2+nrOfDim))]

    if (rank.eq.0) &
         open(unit=ah3,file=filename,form="UNFORMATTED",status="NEW", &
         IOSTAT=ierror)
#ifdef MPI
    ! Distribute runid over nodes
    call MPI_BCAST(ierror,1,MPI_INTEGER,0,MPI_COMM_NEW,mpi_ierror)
#endif     
    errortest: if (ierror == 0) then
       ranktest: if (rank == 0) then
          ! Header
          write(ah3) banner
          write(ah3) nrOfDim
          write(ah3) neq
          write(ah3) npr
          write(ah3) refinementFactor
          write(ah3) nframe
          write(ah3) gamma
          write(ah3) time*sctime
          close(ah3)
          
          ! First ring: Grid 
          open(unit=ah3,file=filename,form="UNFORMATTED",status="OLD", &
               position="APPEND")
          write(ah3) ex-sx+1,ey-sy+1,ez-sz+1
          write(ah3) x(sx)*scleng,y(sy)*scleng,z(sz)*scleng
          write(ah3) dx*scleng,dy*scleng,dz*scleng
          write(ah3) level
          close(ah3)
          
#ifdef MPI
          ! Send filename to next processor
          if (npr > 1) then
             call MPI_ISSEND(filename,19,MPI_CHARACTER,rank+1,outputcircle, &
                  MPI_COMM_NEW,request,mpi_ierror)
             write(log_unit,*) rank,"sent filename to ",rank+1
             ! Wait for the circle to complete
             call MPI_RECV(filename,19,MPI_CHARACTER,npr-1,outputcircle, &
                  MPI_COMM_NEW,status,mpi_ierror)
             write(log_unit,*) rank,"received filename from ",npr-1
          endif
#endif
          
          ! Second ring: Cell data
          open(unit=ah3,file=filename,form="UNFORMATTED",status="OLD", &
               position="APPEND")
          write(ah3) (((state(i,j,k,RHO)*scdens,i=sx,ex),j=sy,ey),k=sz,ez)
          write(ah3) (((state(i,j,k,RHVX)*scmome,i=sx,ex),j=sy,ey),k=sz,ez)
          write(ah3) (((state(i,j,k,RHVY)*scmome,i=sx,ex),j=sy,ey),k=sz,ez)
          write(ah3) (((state(i,j,k,RHVZ)*scmome,i=sx,ex),j=sy,ey),k=sz,ez)
          write(ah3) (((state(i,j,k,EN)*scener,i=sx,ex),j=sy,ey),k=sz,ez)
          if (neq > neuler) then
             write(ah3) ((((state(i,j,k,ieq),i=sx,ex),j=sy,ey),k=sz,ez), &
                  ieq=neuler+1,neq)
          endif
          close(ah3)
          ! Send filename to next processor
#ifdef MPI
          if (npr > 1) then
             call MPI_ISSEND(filename,19,MPI_CHARACTER,rank+1,outputcircle, &
                  MPI_COMM_NEW,request,mpi_ierror)
             write(log_unit,*) rank,"sent filename to ",rank+1
             ! Wait for the circle to complete
             call MPI_RECV(filename,19,MPI_CHARACTER,npr-1,outputcircle, &
                  MPI_COMM_NEW,status,mpi_ierror)
             write(log_unit,*) rank,"received filename from ",npr-1," end of loop."
          endif
#endif
#ifdef MPI
       else
          
          ! First ring: Grids
          ! Receive filename from previous processor
          call MPI_RECV(filename,19,MPI_CHARACTER,rank-1,outputcircle, &
               MPI_COMM_NEW,status,mpi_ierror)
          write(log_unit,*) rank,"received filename from ",rank-1
          
          if (mpi_ierror == 0) then ! if ok
             open(unit=ah3,file=filename,form="UNFORMATTED",status="OLD", &
                  position="APPEND")
             ! Grid 
             write(ah3) ex-sx+1,ey-sy+1,ez-sz+1
             write(ah3) x(sx)*scleng,y(sy)*scleng,z(sz)*scleng
             write(ah3) dx*scleng,dy*scleng,dz*scleng
             write(ah3) level
             close(ah3)
             
             ! Determine next processor (npr -> 0)
             nextproc=mod(rank+1,npr)
             ! Send filename along
             call MPI_ISSEND(filename,19,MPI_CHARACTER,nextproc, &
                  outputcircle,MPI_COMM_NEW,request,mpi_ierror)
             write(log_unit,*) rank,"sent filename to ",nextproc
          endif
          
          ! Second ring: Cells
          ! Receive filename from previous processor
          call MPI_RECV(filename,19,MPI_CHARACTER,rank-1,outputcircle, &
               MPI_COMM_NEW,status,mpi_ierror)
          write(log_unit,*) rank,"received filename from ",rank-1
          
          if (mpi_ierror == 0) then ! if ok
             open(unit=ah3,file=filename,form="UNFORMATTED",status="OLD", &
                  position="APPEND")
             ! Cells
             write(ah3) (((state(i,j,k,RHO)*scdens,i=sx,ex),j=sy,ey),k=sz,ez)
             write(ah3) (((state(i,j,k,RHVX)*scmome,i=sx,ex),j=sy,ey),k=sz,ez)
             write(ah3) (((state(i,j,k,RHVY)*scmome,i=sx,ex),j=sy,ey),k=sz,ez)
             write(ah3) (((state(i,j,k,RHVZ)*scmome,i=sx,ex),j=sy,ey),k=sz,ez)
             write(ah3) (((state(i,j,k,EN)*scener,i=sx,ex),j=sy,ey),k=sz,ez)
             if (neq > neuler) then
                write(ah3) ((((state(i,j,k,ieq),i=sx,ex),j=sy,ey),k=sz,ez), &
                     ieq=neuler+1,neq)
             endif
             close(ah3)
             
             ! Determine next processor (npr -> 0)
             nextproc=mod(rank+1,npr)
             ! Send filename along
             call MPI_ISSEND(filename,19,MPI_CHARACTER,nextproc, &
                  outputcircle,MPI_COMM_NEW,request,mpi_ierror)
             write(log_unit,*) rank,"sent filename to ",nextproc
          else
             write(log_unit,*) "MPI error in output routine"
             ierror=ierror+1
          endif
#endif
       endif ranktest
    else
       write(log_unit,*) "Error opening output file, continuing without writing."
    endif errortest
    
  end subroutine make_output
  
end module output
