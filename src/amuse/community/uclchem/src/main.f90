! 2022 UCLCHEM v3.0
! The canonical main file where the core code is written is actually found in wrap.f90
! main.f90 just provides a simple fortran interface to the core code so that a binary can
! be built and used directly from the command line.
PROGRAM uclchem

USE uclchemwrap, only: cloud,hot_core,cshock
USE io, only: inputId
USE constants, only: dp
IMPLICIT NONE
    CHARACTER (LEN=100):: modelType
    CHARACTER (LEN=100):: paramFile
    CHARACTER(:), ALLOCATABLE :: paramDict
    REAL(dp) :: abundances(500),dissipationResult
    INTEGER :: success,fileLength
    !Any subset of the parameters can be passed in a file on program start
    !see example.inp
    CALL GET_COMMAND_ARGUMENT(1, modelType)  
    CALL GET_COMMAND_ARGUMENT(2, paramFile)

    OPEN(unit=inputId, file=paramFile, action="read", &
    form="unformatted", access="stream")
    INQUIRE(unit=inputId, size=fileLength)

    ALLOCATE(character(fileLength) :: paramDict)
    READ(inputId) paramDict
    CLOSE(inputId)

    SELECT CASE(modelType)
    CASE("CLOUD")
        CALL cloud(paramDict,"",abundances,success)
    CASE("HOTCORE")
        !call hot_core (temp_indx,maxTemp,....)
        CALL hot_core(3,3.0d2,paramDict,"",abundances,success)
    CASE("CSHOCK")
        !call cshock(vs,timestep_factor,minimum_temp,....)
        CALL cshock(20.0d0,0.01d0,0.0d0,paramDict,"",abundances,dissipationResult,success)
    CASE default
        write(*,*) 'Model type not recognised'
        WRITE(*,*) 'Supported models are: CLOUD, CSHOCK, HOTCORE'
        STOP
    END SELECT
END PROGRAM uclchem