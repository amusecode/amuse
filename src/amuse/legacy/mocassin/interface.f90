MODULE mocassin_interface
    use common_mod
    use constants_mod
    use dust_mod
    use grid_mod
    use iteration_mod
    use output_mod
    use set_input_mod
    use xSec_mod
    
    implicit none
    
    logical :: areStarsDefined = .false.
    
    type(grid_type) :: grid3D(maxGrids)       ! the 3D Cartesian  grid
    
CONTAINS

    SUBROUTINE set_default_values()
        IMPLICIT NONE
        ! set default values and set non oprional values to 0 or 0. or "zero"
        ! copy pasted from set_input_mod.f90

        lgIsotropic   = .false.       
        lgRecombination = .false.
        lgAutoPackets = .false.
        lgMultiDustChemistry = .false.
        lgMultiChemistry = .false.
        lgTalk        = .false.
        lgHdenConstant= .false.
        lgDfile       = .false.
        lgDebug       = .false.
        lgDlaw        = .false.
        lgDust        = .false.
        lgDustConstant= .false.
        lgGas         = .true.
        lgMdMg        = .false.
        lgMdMh        = .false.
        lgNeInput     = .false.
        lgNeutral     = .true.
        lgOutput      = .false. 
        lgPlaneIonization = .false.
        lg1D          = .false.
        lgDustScattering = .true.
        lgSymmetricXYZ= .false.
        lgEcho        = .false.
        lgNosource    = .false.

        nPhotonsDiffuse = 0        
        nStars        = 0
        nxIn(:)       = 0
        nyIn(:)       = 0
        nzIn(:)       = 0
        nxIn(1)       = 30
        nyIn(1)       = 30
        nzIn(1)       = 30        
        maxIterateMC  = 30
        maxPhotons    = 0
        nbins         = 600        
        MdMgValue     = 0.
        MdMgFile      = "none"               
        nAngleBins    = 0        
        nGrids        = 1
        NdustValue     = 0.
        NdustFile      = "none"
        nstages        = 7
        
        ! qheat
        Tmax           =700.
        nTbins         =300
        lgWritePss     =.false.
        minaQHeat     = 1.e-3
        minConvQHeat  = 99.

        fillingFactor = 1.
        contCube      = -1.
        convIncPercent= 0.
        convWriteGrid = 0.
        nPhotIncrease = 0.

        densityLaw    = 0.

        !multiPhotoSources = "none"
        densityFile   = "none"
        dustFile      = "none"
        gridList      = "none"

        nu0           = 0.  
        Ldiffuse      = 0. 
        Tdiffuse      = 0.
        shapeDiffuse  = "none"
        Hdensity      = 0.
        H0Start       = 3.e-5
        LPhot         = 0.
        minConvergence= 95.
        NeStart       = 0.
        nuMax         = 15.
        nuMin         = 1.001e-5
        nuStepSize    = 0.075
        resLinesTransfer = 101.
        Rnx           = -1.
        Rny           = -1.
        Rnz           = -1.
        R_in          = -1.
        R_out         = 0.
        TeStart       = 10000.
        XHILimit      = 0.05
        
        lg2D = .false.
        
    END SUBROUTINE
    
    FUNCTION cleanup_code()
        IMPLICIT NONE
        include 'mpif.h'
        INTEGER cleanup_code
        INTEGER :: iGrid
        
        ! free all space allocated to the 3D grid
        do iGrid=1, nGrids
            call freeGrid(grid3D(iGrid))
        end do
    
        cleanup_code=0
    END FUNCTION

    FUNCTION recommit_parameters()
        IMPLICIT NONE
        include 'mpif.h'
        INTEGER recommit_parameters
        recommit_parameters=-1
    END FUNCTION
    
    FUNCTION iterate()
        IMPLICIT NONE
        include 'mpif.h'
        INTEGER iterate
        
        call MCIterationDriver(grid3D(1:nGrids))
        iterate=-1
    END FUNCTION
    
    FUNCTION initialize_code()
        IMPLICIT NONE
        include 'mpif.h'
        INTEGER initialize_code
        
        call set_default_values()
        
        call mpi_comm_rank(MPI_COMM_WORLD, taskid, ierr)
        call mpi_comm_size(MPI_COMM_WORLD, numtasks, ierr)
     
        
        allocate(abundanceFile(0:100))
        
        initialize_code=0
    END FUNCTION
    
    function redirect_outputs_to(stdoutfile, stderrfile)
        implicit none
        
        include 'mpif.h'
        
        character(LEN=*) , INTENT(IN) :: stdoutfile, stderrfile
        character(1024) :: fullname
        integer :: rank, redirect_outputs_to
        
        CLOSE(UNIT=5)
        
        call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
        
        CLOSE(UNIT=6)
        
        write( fullname, '(A,".",I3.3)' )  stdoutfile, rank
        OPEN(UNIT=6, FILE=TRIM(fullname), ACCESS="APPEND")
        CLOSE(UNIT=0)
        
        write( fullname, '(A,".",I3.3)' )  stderrfile, rank
        OPEN(UNIT=0, FILE=TRIM(fullname), ACCESS="APPEND")
        redirect_outputs_to = 0
    end function
    
    
    FUNCTION commit_parameters()
        IMPLICIT NONE
        include 'mpif.h'
        INTEGER commit_parameters, iGrid
        integer             :: err              ! allocation error status
        
        print *, "nstages ::", nstages
        allocate(lgDataAvailable(3:nElements, 1:nstages), stat=err)
        if (err /= 0) then
           print*, '! commit_parameters: allocation error for lgDataAvailable pointer'
           commit_parameters = -1
           return
        end if
        
        do iGrid = 1, nGrids
            call initCartesianGrid(grid3D(iGrid), nxIn(iGrid), nyIn(iGrid), nzIn(iGrid))
        end do
        
        commit_parameters=0

    END FUNCTION
    
    
    
    FUNCTION commit_particles()
        IMPLICIT NONE
        INTEGER commit_particles
        INTEGER i
        
        do i = 1, nStars
           nPhotons(i) = NphotonsTot/nStars
           deltaE(i) = Lstar(i)/nPhotons(i)
        end do  
        
    END FUNCTION

    FUNCTION commit_grid()
        IMPLICIT NONE
        INTEGER commit_grid
        INTEGER i, iGrid
        DOUBLE PRECISION test
        Integer :: nxA,nyA,nzA
        
        
        
        ! initialize opacities x sections array
        call initXSecArray()
        
        ! set the ionzing continuum according to the contShape variable
        call setContinuum()
        
        call fillGrid(grid3D(1:nGrids))
        
        do i = 1, nStars
            nxA = size(grid3D(1)%xAxis)
            nyA = size(grid3D(1)%yAxis)
            nzA = size(grid3D(1)%zAxis)
            starPosition(i)%x = starPosition(i)%x/grid3D(1)%xAxis(nxA)
            starPosition(i)%y = starPosition(i)%y/grid3D(1)%yAxis(nyA)
            starPosition(i)%z = starPosition(i)%z/grid3D(1)%zAxis(nzA)
        end do
        
        
        call setStarPosition(grid3D(1)%xAxis,grid3D(1)%yAxis,grid3D(1)%zAxis, grid3D(1:nGrids))
        
        if (taskid==0) then
           do iGrid = 1, nGrids
              print*, 'Grid : ', iGrid 
              print*, 'active cells: ', grid3D(iGrid)%nCells
           end do
        end if

        ! prepare atomica data stuff
        if (lgGas) call makeElements()
        print*, 'active elements: ', nElementsUsed
        
        ! if grains are included, calculate the dust opacity     
        if (lgDust) then
           if (taskid==0) print*, '! mocassin: calling dustDriver'
           do iGrid = 1, nGrids
              call dustDriver(grid3D(iGrid))
           end do
           if (taskid==0) print*, '! mocassin: dustDriver done'
        end if
        
        if (Ldiffuse>0.) then
           if (taskid==0) print*, '! mocassin: calling setLdiffuse'
           call setLdiffuse(grid3D(1:nGrids))
           test=0.
           do iGrid = 1, nGrids
              do i = 1, grid3D(igrid)%ncells
                 test = test+grid3d(igrid)%LdiffuseLoc(i)
              end do
           end do
           if (taskid==0) then
              print*, '! mocassin: setLdiffuse done, total Ldiffuse: ', test
           end if
        end if
        
        
    END FUNCTION

    FUNCTION setup_mesh( ni, nj, nk, xlength, ylength, zlength, grid_index)
        IMPLICIT NONE
        INTEGER setup_mesh, ni, nj, nk, grid_index
        DOUBLE PRECISION xlength, ylength, zlength
        nxIn(grid_index) = ni
        nyIn(grid_index) = nj
        nzIn(grid_index) = nk  
        IF (grid_index .EQ. 1) THEN
            Rnx = xlength
            Rny = ylength
            Rnz = zlength
        END IF
        setup_mesh = 0
    END FUNCTION
    
    
    FUNCTION setup_auto_convergence(convergence_level_increase, number_of_photons_increase , maximum_number_of_photons)
        IMPLICIT NONE
        DOUBLE PRECISION , INTENT(IN) :: convergence_level_increase, number_of_photons_increase
        INTEGER , INTENT(IN) :: maximum_number_of_photons
        INTEGER :: setup_auto_convergence
        lgAutoPackets = .true.
        convIncPercent = convergence_level_increase
        nPhotIncrease  = number_of_photons_increase
        maxPhotons     = maximum_number_of_photons
        setup_auto_convergence = 0
    END FUNCTION

    FUNCTION has_auto_convergence(value)
        IMPLICIT NONE
        LOGICAL , INTENT(OUT) :: value
        INTEGER :: has_auto_convergence
        value = lgAutoPackets
        has_auto_convergence = 0
    END FUNCTION
    
    FUNCTION uset_auto_convergence()
        IMPLICIT NONE
        INTEGER :: uset_auto_convergence
        lgAutoPackets = .false.
        uset_auto_convergence = 0
    END FUNCTION

    FUNCTION get_position_of_index( &
            i, j, k, index_of_grid, &
            x, y, z)
        IMPLICIT NONE
        INTEGER :: i, j, k, index_of_grid, get_position_of_index
        DOUBLE PRECISION :: x, y, z
        PRINT *, "ngrids", nGrids
        IF(index_of_grid .GT. nGrids) THEN
            get_position_of_index = -1
            return
        END IF
        PRINT *, index_of_grid, grid3D(index_of_grid)%nx
        IF(i .GT. grid3D(index_of_grid)%nx .OR. i .LT. 1) THEN
            get_position_of_index = -2
            return
        END IF
        IF(j .GT. grid3D(index_of_grid)%ny .OR. j .LT. 1) THEN
            get_position_of_index = -3
            return
        END IF
        IF(k .GT. grid3D(index_of_grid)%nz .OR. k .LT. 1) THEN
            get_position_of_index = -4
            return
        END IF
        
        x = grid3D(index_of_grid)%xAxis(i)
        y = grid3D(index_of_grid)%xAxis(j)
        z = grid3D(index_of_grid)%xAxis(k)
        
        get_position_of_index=0
    END FUNCTION
    
    FUNCTION set_abundancies_filename(index, filename)
        IMPLICIT NONE
        CHARACTER(LEN=*) filename
        INTEGER index, set_abundancies_filename
        
        
        !abundanceFile(index) = TRIM(filename)
        abundanceFile = TRIM(filename)
        set_abundancies_filename = 0
    END FUNCTION

    FUNCTION get_abundancies_filename(index, filename)
        IMPLICIT NONE
        CHARACTER(LEN=*) filename
        INTEGER index, get_abundancies_filename
        
        filename = abundanceFile(index)
        get_abundancies_filename = 0
    END FUNCTION
    
    
    FUNCTION set_input_directory(path)
        IMPLICIT NONE
        CHARACTER(LEN=*) path
        INTEGER set_input_directory
        home = TRIM(path)
        set_input_directory = 0
    END FUNCTION
    

    FUNCTION get_input_directory(path)
        IMPLICIT NONE
        CHARACTER(LEN=*) path
        INTEGER get_input_directory
        
        path = home
        get_input_directory = 0
    END FUNCTION
    
    FUNCTION set_constant_hydrogen_density(value)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: value
        INTEGER set_constant_hydrogen_density
        lgHdenConstant = .true.
        Hdensity = value
        set_constant_hydrogen_density = 0
    END FUNCTION

    FUNCTION get_constant_hydrogen_density(value)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(OUT) :: value
        INTEGER get_constant_hydrogen_density
        value = Hdensity
        get_constant_hydrogen_density = 0
    END FUNCTION

    FUNCTION has_constant_hydrogen_density(value)
        IMPLICIT NONE
        LOGICAL, INTENT(OUT) :: value
        INTEGER has_constant_hydrogen_density
        
        value = lgHdenConstant
        has_constant_hydrogen_density = 0
    END FUNCTION

    FUNCTION define_stars(x, y, z, temperature, luminocity, n)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        DOUBLE PRECISION, INTENT(IN), DIMENSION(N) :: temperature, luminocity, x, y, z
        !INTEGER, INTENT(IN), DIMENSION(N) :: id
        INTEGER :: define_stars, i
        
        
        if (areStarsDefined) then
            deallocate(deltaE)
            deallocate(spID)        
            deallocate(tStep)        
            deallocate(TStellar)        
            deallocate(Lstar)
            deallocate(ContShape)        
            deallocate(ContShapeIn)        
            deallocate(nPhotons)        
            deallocate(starPosition)       
        end if
        
        nStars = n
        areStarsDefined = .true.
        
        allocate(deltaE(0:nStars))
        allocate(spID(0:nStars))        
        allocate(tStep(nStars))        
        allocate(TStellar(0:nStars))        
        allocate(Lstar(nStars))
        allocate(ContShape(0:nStars))        
        allocate(ContShapeIn(0:nStars))        
        allocate(nPhotons(nStars))        
        allocate(starPosition(nStars))       

        TStellar=0.
        Lstar=0.
        ContShape='none'
        ContShapeIn='none'
        spID='none'
        tStep=0.
        nPhotons=0.
        deltaE=0.
        print *, "DELTA E::", deltaE

        do i = 1, nStars
            starPosition(i)%x = x(i)
            starPosition(i)%y = y(i)
            starPosition(i)%z = z(i)
            ContShape(i)='blackbody'
            spID='blackbody'
            ContShapeIn(i)=ContShape(i)
            TStellar(i) = temperature(i)
            Lstar(i) = luminocity(i)
        end do
        
        define_stars = 0
    END FUNCTION
    
    FUNCTION set_total_number_of_photons(value)
        IMPLICIT NONE
        integer, INTENT(IN) :: value
        INTEGER set_total_number_of_photons
        nPhotonsTot = value
        set_total_number_of_photons = 0
    END FUNCTION

    FUNCTION get_total_number_of_photons(value)
        IMPLICIT NONE
        integer, INTENT(OUT) :: value
        INTEGER get_total_number_of_photons
        value = nPhotonsTot
        get_total_number_of_photons = 0
    END FUNCTION

    FUNCTION set_total_number_of_points_in_frequency_mesh(value)
        IMPLICIT NONE
        integer, INTENT(IN) :: value
        INTEGER set_total_number_of_points_in_frequency_mesh
        nbins = value
        set_total_number_of_points_in_frequency_mesh = 0
    END FUNCTION

    FUNCTION get_total_number_of_points_in_frequency_mesh(value)
        IMPLICIT NONE
        integer, INTENT(OUT) :: value
        INTEGER get_total_number_of_points_in_frequency_mesh
        value = nbins
        get_total_number_of_points_in_frequency_mesh = 0
    END FUNCTION
    
    FUNCTION set_initial_nebular_temperature(value)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: value
        INTEGER set_initial_nebular_temperature
        TeStart = value
        set_initial_nebular_temperature = 0
    END FUNCTION

    FUNCTION get_initial_nebular_temperature(value)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(OUT) :: value
        INTEGER get_initial_nebular_temperature
        value = TeStart
        get_initial_nebular_temperature = 0
    END FUNCTION
    
    FUNCTION set_symmetricXYZ(value)
        IMPLICIT NONE
        logical, INTENT(IN) :: value
        INTEGER set_symmetricXYZ
        lgSymmetricXYZ = value
        set_symmetricXYZ = 0
    END FUNCTION

    FUNCTION get_symmetricXYZ(value)
        IMPLICIT NONE
        logical, INTENT(OUT) :: value
        INTEGER get_symmetricXYZ
        value = lgSymmetricXYZ
        get_symmetricXYZ = 0
    END FUNCTION
    
    FUNCTION set_maximum_number_of_monte_carlo_iterations(value)
        IMPLICIT NONE
        integer, INTENT(IN) :: value
        INTEGER set_maximum_number_of_monte_carlo_iterations
        maxIterateMC = value
        set_maximum_number_of_monte_carlo_iterations = 0
    END FUNCTION

    FUNCTION get_maximum_number_of_monte_carlo_iterations(value)
        IMPLICIT NONE
        integer, INTENT(OUT) :: value
        INTEGER get_maximum_number_of_monte_carlo_iterations
        value = maxIterateMC
        get_maximum_number_of_monte_carlo_iterations = 0
    END FUNCTION
    
    FUNCTION set_minimum_convergence_level(value)
        IMPLICIT NONE
        double precision, INTENT(IN) :: value
        INTEGER set_minimum_convergence_level
        minConvergence = value
        set_minimum_convergence_level = 0
    END FUNCTION

    FUNCTION get_minimum_convergence_level(value)
        IMPLICIT NONE
        double precision, INTENT(OUT) :: value
        INTEGER get_minimum_convergence_level
        value = minConvergence
        get_minimum_convergence_level = 0
    END FUNCTION
    
    FUNCTION set_high_limit_of_the_frequency_mesh(value)
        IMPLICIT NONE
        double precision, INTENT(IN) :: value
        INTEGER set_high_limit_of_the_frequency_mesh
        nuMax = value
        set_high_limit_of_the_frequency_mesh = 0
    END FUNCTION

    FUNCTION get_high_limit_of_the_frequency_mesh(value)
        IMPLICIT NONE
        double precision, INTENT(OUT) :: value
        INTEGER get_high_limit_of_the_frequency_mesh
        value = nuMax
        get_high_limit_of_the_frequency_mesh = 0
    END FUNCTION


    FUNCTION set_low_limit_of_the_frequency_mesh(value)
        IMPLICIT NONE
        double precision, INTENT(IN) :: value
        INTEGER set_low_limit_of_the_frequency_mesh
        nuMin = value
        set_low_limit_of_the_frequency_mesh = 0
    END FUNCTION

    FUNCTION get_low_limit_of_the_frequency_mesh(value)
        IMPLICIT NONE
        double precision, INTENT(OUT) :: value
        INTEGER get_low_limit_of_the_frequency_mesh
        value = nuMin
        get_low_limit_of_the_frequency_mesh = 0
    END FUNCTION
    
    FUNCTION set_inner_radius_of_the_ionised_region(value)
        IMPLICIT NONE
        double precision, INTENT(IN) :: value
        INTEGER set_inner_radius_of_the_ionised_region
        R_in = value
        set_inner_radius_of_the_ionised_region = 0
    END FUNCTION

    FUNCTION get_inner_radius_of_the_ionised_region(value)
        IMPLICIT NONE
        double precision, INTENT(OUT) :: value
        INTEGER get_inner_radius_of_the_ionised_region
        value = R_in
        get_inner_radius_of_the_ionised_region = 0
    END FUNCTION

    FUNCTION set_outer_radius_of_the_ionised_region(value)
        IMPLICIT NONE
        double precision, INTENT(IN) :: value
        INTEGER set_outer_radius_of_the_ionised_region
        R_out = value
        set_outer_radius_of_the_ionised_region = 0
    END FUNCTION

    FUNCTION get_outer_radius_of_the_ionised_region(value)
        IMPLICIT NONE
        double precision, INTENT(OUT) :: value
        INTEGER get_outer_radius_of_the_ionised_region
        value = R_out
        get_outer_radius_of_the_ionised_region = 0
    END FUNCTION

    FUNCTION set_convergence_limit(value)
        IMPLICIT NONE
        double precision, INTENT(IN) :: value
        INTEGER set_convergence_limit
        XHILimit = value
        set_convergence_limit = 0
    END FUNCTION

    FUNCTION get_convergence_limit(value)
        IMPLICIT NONE
        double precision, INTENT(OUT) :: value
        INTEGER get_convergence_limit
        value = XHILimit
        get_convergence_limit = 0
    END FUNCTION
    
    FUNCTION set_number_of_ionisation_stages(value)
        IMPLICIT NONE
        integer, INTENT(IN) :: value
        INTEGER set_number_of_ionisation_stages
        nstages = value
        set_number_of_ionisation_stages = 0
    END FUNCTION

    FUNCTION get_number_of_ionisation_stages(value)
        IMPLICIT NONE
        integer, INTENT(OUT) :: value
        INTEGER get_number_of_ionisation_stages
        value = nstages
        get_number_of_ionisation_stages = 0
    END FUNCTION

    FUNCTION set_write_snapshot_every_iteration(value)
        IMPLICIT NONE
        logical, INTENT(IN) :: value
        INTEGER set_write_snapshot_every_iteration
        lgOutput = value
        set_write_snapshot_every_iteration = 0
    END FUNCTION

    FUNCTION get_write_snapshot_every_iteration(value)
        IMPLICIT NONE
        logical, INTENT(OUT) :: value
        INTEGER get_write_snapshot_every_iteration
        value = lgOutput
        get_write_snapshot_every_iteration = 0
    END FUNCTION
        

    FUNCTION get_number_of_elements_used(value)
        IMPLICIT NONE
        integer, INTENT(OUT) :: value
        INTEGER get_number_of_elements_used
        value = nElementsUsed
        get_number_of_elements_used = 0
    END FUNCTION
    
    FUNCTION get_grid_electron_temperature(i,j,k,index_of_grid,electron_temperature, n)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(IN), DIMENSION(N) :: i,j,k, index_of_grid
        
        double precision, INTENT(OUT), DIMENSION(N) :: electron_temperature
        
        INTEGER get_grid_electron_temperature
        INTEGER :: index, active_cell_index


        
        DO index = 1,n
            
            IF (index_of_grid(index).GT.nGrids) THEN
                get_grid_electron_temperature = -1
                return
            END IF
            active_cell_index = grid3D(index_of_grid(index))%active(i(index), j(index), k(index))
            if (active_cell_index .eq. 0) then
                
                    electron_temperature(index) = 0
                
            else
                
                    electron_temperature(index) = grid3D(index_of_grid(index))%Te(active_cell_index)
                
            end if
        END DO
        get_grid_electron_temperature = 0
    END FUNCTION

    FUNCTION get_grid_electron_density(i,j,k,index_of_grid,electron_density, n)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(IN), DIMENSION(N) :: i,j,k, index_of_grid
        
        double precision, INTENT(OUT), DIMENSION(N) :: electron_density
        
        INTEGER get_grid_electron_density
        INTEGER :: index, active_cell_index


        
        DO index = 1,n
            
            IF (index_of_grid(index).GT.nGrids) THEN
                get_grid_electron_density = -1
                return
            END IF
            active_cell_index = grid3D(index_of_grid(index))%active(i(index), j(index), k(index))
            if (active_cell_index .eq. 0) then
                
                    electron_density(index) = 0
                
            else
                
                    electron_density(index) = grid3D(index_of_grid(index))%Ne(active_cell_index)
                
            end if
        END DO
        get_grid_electron_density = 0
    END FUNCTION

    FUNCTION get_grid_active(i,j,k,index_of_grid,is_active, n)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(IN), DIMENSION(N) :: i,j,k, index_of_grid
        
        LOGICAL, INTENT(OUT), DIMENSION(N) :: is_active
        
        INTEGER get_grid_active
        INTEGER :: index, active_cell_index


        
        DO index = 1,n
            
            IF (index_of_grid(index).GT.nGrids) THEN
                is_active(index) = .FALSE.
                CONTINUE
            END IF
            active_cell_index = grid3D(index_of_grid(index))%active(i(index), j(index), k(index))
            if (active_cell_index .eq. 0) then
                is_active(index) = .FALSE.
            else
                is_active(index) = .TRUE.
            end if
        END DO
        get_grid_active = 0
    END FUNCTION
    
    FUNCTION set_emit_rate_of_photons(value)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: value
        INTEGER set_emit_rate_of_photons
        LPhot = value
        set_emit_rate_of_photons = 0
    END FUNCTION

    FUNCTION get_emit_rate_of_photons(value)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(OUT) :: value
        INTEGER get_emit_rate_of_photons
        value = LPhot
        get_emit_rate_of_photons = 0
    END FUNCTION
END MODULE


