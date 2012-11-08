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
    logical :: is_grid_committed = .false.
    real, pointer :: hydrogen_density_input(:,:,:)
    integer :: totPercentOld
    integer :: total_number_of_photons = 1000000
    
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
        lginputGasMass = .false.
        lgDoSeedRandom = .true.

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
        
        totPercentOld = 0.
        
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
        integer :: error
        
        ! call MCIterationDriver(grid3D(1:nGrids))
        
        
        
        do while (.true.)
           error = step()
           

           if (nIterateMC > 1 .and. totPercent < 95. .and. lgAutoPackets & 
                & .and. nPhotonsTot < maxPhotons .and. totPercentOld > 0.) then

              if ( (totPercent-totPercentOld)/totPercentOld <= convIncPercent ) then
                 nPhotons = nPhotons*nPhotIncrease
                 deltaE   = deltaE/nPhotIncrease

                 if (taskid==0) &
                      & print*, "! iterateMC: [talk] Total number of energy packets &
                      &increased to ", nPhotons

                 if (Ldiffuse>0.) then
                    nPhotonsDiffuseLoc = nPhotonsDiffuseLoc*nPhotIncrease
                    if (taskid==0) &
                         & print*, "! iterateMC: [talk] number of diffuse energy packets &
                         &per cell increased to ", nPhotonsDiffuseLoc
                 end if              

              end if              
           end if
           
        
           totPercentOld = totPercent

           if ( totPercent >= minConvergence ) then
              
              if (lgOutput .and. taskid==0) then
                 print*, "! iterateMC: [talk] convergence reached after ", &
                      &                           nIterateMC, " iterations. & Finishing up ... "
                 
                 ! output results at this itearation stage (every 3 iterations)

                 call writeGrid(grid3D(1:nGrids))
                 if (lgGas) call outputGas(grid3D(1:nGrids)) 
              end if
              
              call mpi_barrier(mpi_comm_world, ierr)

              exit
           else if (nIterateMC >= maxIterateMC) then
              
              if (lgOutput .and. taskid==0) then

                 print*, " ! iterateMC: maximum number of iterations reached. Finishing up ... ",&
                      & maxIterateMC

                 call writeGrid(grid3D(1:nGrids))
                 if (lgGas) call outputGas(grid3D(1:nGrids))
              end if

              call mpi_barrier(mpi_comm_world, ierr)
                
              exit
              
           else
              
              if (lgOutput .and. taskid == 0 ) then
                 ! output results at this itearation stage (every ? iterations)
                 if ( mod(nIterateMC, 1) == 0 ) then
                    if (lgGas) call outputGas(grid3D(1:nGrids)) 
                 end if
              end if
              
              call mpi_barrier(mpi_comm_world, ierr)
                 

              
           end if
           
        end do
        
        iterate=0
    END FUNCTION
    
    FUNCTION step()
        IMPLICIT NONE
        include 'mpif.h'
        INTEGER step
        integer :: error, iStar
        
        !reinitialize the number of photons per star
        !if a change in the total number of photons occured
        if(total_number_of_photons.NE.nPhotonsTot) then
            nPhotonsTot = total_number_of_photons
            do iStar=1, nStars              
                nPhotons(iStar) = nPhotonsTot / nStars
                deltaE(iStar) = Lstar(iStar)/nPhotons(iStar)
            end do
            
            ! deal with round-off in number of photons by summing the counts
            nPhotonsTot = 0.0
            do iStar=1, nStars           
                nPhotonsTot = nPhotonsTot+nPhotons(iStar) 
            end do
            total_number_of_photons = nPhotonsTot
        end if
        
        
        error = step_monte_carlo(grid3D(1:nGrids))
        print *, "step error", error
        if (error.NE.0) then
            step = error
            return
        end if
        
       
        error = calculate_convergence(grid3D(1:nGrids))
        
        nPhotonsTot = 0.0
        do iStar=1, nStars             
            nPhotonsTot = nPhotonsTot+nPhotons(iStar)
        end do
        total_number_of_photons = nPhotonsTot
        
        totPercentOld = totPercent
        
        
        ! step up MC iterations counter
        nIterateMC = nIterateMC + 1
              
        
        step=error
    END FUNCTION
    
    
    FUNCTION initialize_code()
        IMPLICIT NONE
        include 'mpif.h'
        INTEGER initialize_code
                
        call set_default_values()
        
        call mpi_comm_rank(MPI_COMM_WORLD, taskid, ierr)
        call mpi_comm_size(MPI_COMM_WORLD, numtasks, ierr)
     
        
        allocate(abundanceFile(0:100))
        abundanceFile = ''
        
        initialize_code=0
    END FUNCTION
    
    function set_random_seed(value)
        integer, intent(in) :: value
        integer set_random_seed
        integer, allocatable :: seed(:)
        integer :: seedSize, ierr
        
        call random_seed(seedSize) 
        set_random_seed = 0
        allocate(seed(1:seedSize), stat= ierr)
        if (ierr /= 0) then
            print*, "can't allocate array memory: seed"
            set_random_seed = -1
            return
        end if

        seed = value
        lgDoSeedRandom = .false. ! no random seeding anywhere else
        
        call random_seed(put = seed)
        
        deallocate(seed)
        
    end function
    
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
        INTEGER :: commit_parameters, iGrid      
        integer :: err              ! allocation error status
        
        print *, "nstages ::", nstages
        allocate(lgDataAvailable(3:nElements, 1:nstages), stat=err)
        if (err /= 0) then
           print*, '! commit_parameters: allocation error for lgDataAvailable pointer'
           commit_parameters = -1
           return
        end if
        
        do iGrid = 1, nGrids
            print *, nxIn(iGrid), nyIn(iGrid), nzIn(iGrid)
            call initCartesianGrid(grid3D(iGrid), nxIn(iGrid), nyIn(iGrid), nzIn(iGrid))
        
            
        end do
        
        commit_parameters = setup_grid(grid3D(1:nGrids))
        
        allocate(hydrogen_density_input(1:nxIn(1), 1:nyIn(1), 1:nzIn(1)), stat = err)
        if (err /= 0) then
            print *, "! can't allocate hydrogen input grid memory"
            commit_parameters = -2
            return
        end if
        
        hydrogen_density_input = 0.0
        
        is_grid_committed = .false. 
           
        commit_parameters=0

    END FUNCTION
    
    FUNCTION setup_grid(grid)
        
        type(grid_type), dimension(:),intent(inout) :: grid
        integer setup_grid
        integer i, iG, err
        
        
        if (lgSymmetricXYZ ) then

            ! create the grid axes, forcing it to have grid points at the
            ! centre
            do i = 1, grid(1)%nx
                grid(1)%xAxis(i) = real(i-1)/real(grid(1)%nx-1) 
                grid(1)%xAxis(i) = grid(1)%xAxis(i) * Rnx
            end do

            do i = 1, grid(1)%ny
                grid(1)%yAxis(i) = real(i-1)/real(grid(1)%ny-1) 
                grid(1)%yAxis(i) = grid(1)%yAxis(i) * Rny
            end do

            do i = 1, grid(1)%nz
                grid(1)%zAxis(i) = real(i-1)/real(grid(1)%nz-1) 
                grid(1)%zAxis(i) = grid(1)%zAxis(i) * Rnz
            end do

        else ! not lgSymmetricXYZ           

            if (mod(grid(1)%nx,2) == 0.) then
                print*, "! fillGrid: the automatic grid option &
                  & requires odd integer nx if not symmetric"
                setup_grid = -1
                return
            end if

            if (mod(grid(1)%ny,2) == 0.) then
                print*, "! fillGrid: the automatic grid option &
                  & requires odd integer ny if not symmetric"
                setup_grid = -1
                return
            end if

            if (mod(grid(1)%nz,2) == 0.) then
                print*, "! fillGrid: the automatic grid option &
                  & requires odd integer nz if not symmetric"
                setup_grid = -1
                return
            end if
           
            ! create the grid axes, forcing it to have grid points at the
            ! centre
            do i = 1, grid(1)%nx
                grid(1)%xAxis(i) = 2.*real(i-1)/real(grid(1)%nx-1) - 1.
                grid(1)%xAxis(i) = grid(1)%xAxis(i) * Rnx
            end do


            do i = 1, grid(1)%ny
                grid(1)%yAxis(i) = 2.*real(i-1)/real(grid(1)%ny-1) - 1.
                grid(1)%yAxis(i) = grid(1)%yAxis(i) * Rny
            end do

            do i = 1, grid(1)%nz
                grid(1)%zAxis(i) = 2.*real(i-1)/real(grid(1)%nz-1) - 1.
                grid(1)%zAxis(i) = grid(1)%zAxis(i) * Rnz
            end do
        
        end if
        
        call locate(grid(1)%xAxis, 0., iOrigin)
        call locate(grid(1)%yAxis, 0., jOrigin)
        call locate(grid(1)%zAxis, 0., kOrigin)

        if (taskid == 0) print*, 'Origin at mother grid cell:  ' , iOrigin, jOrigin, kOrigin
        
        allocate(dl(nGrids), stat = err)
        if (err /= 0) then
           print*, "! fillGrid: can't allocate dl memory"
           setup_grid = -1
           return
        end if
        dl = 0.
        
        do iG = 1, nGrids
           ! find geometric corrections
           grid(iG)%geoCorrX = (grid(iG)%xAxis(grid(iG)%nx) - grid(iG)%xAxis(grid(iG)%nx-1))/2.
           if (.not. lg1D) then
              grid(iG)%geoCorrY = (grid(iG)%yAxis(grid(iG)%ny) - grid(iG)%yAxis(grid(iG)%ny-1))/2.
              grid(iG)%geoCorrZ = (grid(iG)%zAxis(grid(iG)%nz) - grid(iG)%zAxis(grid(iG)%nz-1))/2.
           else
              grid(iG)%geoCorrY = 0.
              grid(iG)%geoCorrZ = 0.
           end if
 
           if (taskid==0) print*, "Geometric grid corrections for grid ", &
                & iG, ' : ', grid(iG)%geoCorrX, grid(iG)%geoCorrY, grid(iG)%geoCorrZ
           
           ! find linear increment
           dl(iG) =  abs(grid(iG)%xAxis(2) - grid(iG)%xAxis(1))
           do i = 2, grid(iG)%nx-1
              dl(iG) = min(dl(iG), abs(grid(iG)%xAxis(i+1)-grid(iG)%xAxis(i)) )
           end do
           do i = 1, grid(iG)%ny-1
              dl(iG) = min(dl(iG), abs(grid(iG)%yAxis(i+1)-grid(iG)%yAxis(i)) )
           end do
           do i = 1, grid(iG)%nz-1
              dl(iG) = min(dl(iG), abs(grid(iG)%zAxis(i+1)-grid(iG)%zAxis(i)) )
           end do
           dl(iG) = dl(iG)/50.                                                                                                 
        end do
        
        setup_grid = 0
    END FUNCTION
    
    FUNCTION commit_particles()
        IMPLICIT NONE
        INTEGER commit_particles
        INTEGER i
        Integer :: nxA,nyA,nzA
        
        nPhotonsTot = total_number_of_photons
        
        do i = 1, nStars
           nPhotons(i) = nPhotonsTot/nStars
           deltaE(i) = Lstar(i)/nPhotons(i)
        end do
        
        ! initialize opacities x sections array
        call initXSecArray()
        
        ! set the ionzing continuum according to the contShape variable
        call setContinuum()
        
        
        do i = 1, nStars
            nxA = size(grid3D(1)%xAxis)
            nyA = size(grid3D(1)%yAxis)
            nzA = size(grid3D(1)%zAxis)
            starPosition(i)%x = starPosition(i)%x/grid3D(1)%xAxis(nxA)
            starPosition(i)%y = starPosition(i)%y/grid3D(1)%yAxis(nyA)
            starPosition(i)%z = starPosition(i)%z/grid3D(1)%zAxis(nzA)
        end do
        
        
        call setStarPosition(grid3D(1)%xAxis,grid3D(1)%yAxis,grid3D(1)%zAxis, grid3D(1:nGrids))
        
        commit_particles=0
        
    END FUNCTION
    
    
    FUNCTION recommit_particles()
        IMPLICIT NONE
        INTEGER recommit_particles
        INTEGER i
        
        nPhotonsTot = total_number_of_photons
        
        do i = 1, nStars
           nPhotons(i) = nPhotonsTot/nStars
           deltaE(i) = Lstar(i)/nPhotons(i)
        end do

        recommit_particles=0
        
    END FUNCTION
    
     FUNCTION setup_subgrid_references(grid)
        
        type(grid_type), dimension(:),intent(inout) :: grid
        integer :: setup_subgrid_references
        integer i, j, k, jG, iG
        
        setup_subgrid_references = 0
        
        if (nGrids>1) then

           do iG = 1, nGrids
              
              do i = 1, grid(iG)%nx
                 do j = 1, grid(iG)%ny
                    do k = 1, grid(iG)%nz
                 
                       do jG = 2, nGrids
                          if (jG /= iG) then
                             if (lgSymmetricXYZ) then                               

                                if ( ( &
                                     & (grid(iG)%xAxis(i) > grid(jG)%xAxis(1)) .or.&
                                     & (grid(iG)%xAxis(i) >= grid(jG)%xAxis(1) .and. grid(jG)%xAxis(1)==0.) & 
                                     & ) .and. &
                                     & grid(iG)%xAxis(i)<grid(jG)%xAxis(grid(jG)%nx) .and.&
                                     & ( & 
                                     & (grid(iG)%yAxis(j) > grid(jG)%yAxis(1)) .or. & 
                                     & (grid(iG)%yAxis(j) >= grid(jG)%yAxis(1) .and. grid(jG)%yAxis(1)==0. ) & 
                                     & ) .and. &
                                     & grid(iG)%yAxis(j)<grid(jG)%yAxis(grid(jG)%ny) .and. &
                                     & ( & 
                                     & (grid(iG)%zAxis(k) > grid(jG)%zAxis(1)) .or.& 
                                     & (grid(iG)%zAxis(k) >= grid(jG)%zAxis(1) .and. grid(jG)%zAxis(1)==0. ) & 
                                     & ) .and. &
                                     & grid(iG)%zAxis(k)<grid(jG)%zAxis(grid(jG)%nz) ) then
                                   grid(iG)%active(i,j,k) = -jG

                                end if



                          else
                                if ( (grid(iG)%xAxis(i) > grid(jG)%xAxis(1) .and. &
                                     &grid(iG)%xAxis(i)<grid(jG)%xAxis(grid(jG)%nx)) .and.&
                                     & (grid(iG)%yAxis(j)>grid(jG)%yAxis(1) .and. & 
                                     & grid(iG)%yAxis(j)<grid(jG)%yAxis(grid(jG)%ny)) .and. &
                                     & (grid(iG)%zAxis(k)>grid(jG)%zAxis(1) .and. & 
                                     & grid(iG)%zAxis(k)<grid(jG)%zAxis(grid(jG)%nz)) ) then
                                   grid(iG)%active(i,j,k) = -jG
                                end if

                             end if

                          end if
                       end do

                    end do
                 end do
              end do

           end do

        end if
        setup_subgrid_references = 0
        
    END FUNCTION

    FUNCTION commit_grid()
        IMPLICIT NONE
        INTEGER commit_grid
        INTEGER i, iGrid
        DOUBLE PRECISION test
        
        if (lgHdenConstant) then
            call fillGrid(grid3D(1:nGrids))
        else
            commit_grid = setup_mother_grid(grid3D(1))
            if (commit_grid.NE.0) then
                return
            end if
            
            if (grid3D(1)%nCells.EQ.0) then
                commit_grid = -1
                return
            end if 
            
            nPhotonsDiffuseLoc = nPhotonsDiffuse/grid3D(1)%nCells

            commit_grid = setup_subgrid_references(grid3D(1:nGrids))
            
            if (commit_grid.NE.0) then
                return
            end if
        end if
        
        !if(associated(hydrogen_density_input)) deallocate(hydrogen_density_input)
        
        is_grid_committed = .true.
        
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
    
    FUNCTION set_abundancies_filename(filename, index)
        IMPLICIT NONE
        CHARACTER(LEN=*) filename
        INTEGER index, set_abundancies_filename
        
        if (index > 1) then
            abundanceFile(index) = TRIM(filename)
        else
            abundanceFile = TRIM(filename)
        end if
        
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
    
    
    FUNCTION set_output_directory(path)
        IMPLICIT NONE
        CHARACTER(LEN=*) path
        INTEGER set_output_directory
        outputDir = TRIM(path)
        set_output_directory = 0
    END FUNCTION
    

    FUNCTION get_output_directory(path)
        IMPLICIT NONE
        CHARACTER(LEN=*) path
        INTEGER get_output_directory
        
        path = outputDir
        get_output_directory = 0
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


    FUNCTION get_percentage_converged(value)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(OUT) :: value
        INTEGER get_percentage_converged
        
        value = totPercent
        
        get_percentage_converged = 0
    END FUNCTION
    
    FUNCTION has_constant_hydrogen_density(value)
        IMPLICIT NONE
        LOGICAL, INTENT(OUT) :: value
        INTEGER has_constant_hydrogen_density
        
        value = lgHdenConstant
        has_constant_hydrogen_density = 0
    END FUNCTION

    FUNCTION set_has_constant_hydrogen_density(value)
        IMPLICIT NONE
        LOGICAL, INTENT(IN) :: value
        INTEGER set_has_constant_hydrogen_density
        
        lgHdenConstant = value
        set_has_constant_hydrogen_density = 0
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
        INTEGER :: set_total_number_of_photons
        
        total_number_of_photons = value
       
        set_total_number_of_photons = 0
    END FUNCTION

    FUNCTION get_total_number_of_photons(value)
        IMPLICIT NONE
        integer, INTENT(OUT) :: value
        INTEGER get_total_number_of_photons
        value = total_number_of_photons
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
    
    FUNCTION get_grid_hydrogen_density(i,j,k,index_of_grid,hydrogen_density, n)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(IN), DIMENSION(N) :: i,j,k, index_of_grid
        
        double precision, INTENT(OUT), DIMENSION(N) :: hydrogen_density
        
        INTEGER get_grid_hydrogen_density
        INTEGER :: index, active_cell_index

        if (.not.is_grid_committed) then
            DO index = 1,n
                hydrogen_density(index) = hydrogen_density_input(i(index), j(index), k(index))
            END DO
        else
        
            DO index = 1,n
                
                IF (index_of_grid(index).GT.nGrids) THEN
                    get_grid_hydrogen_density = -1
                    return
                END IF
                active_cell_index = grid3D(index_of_grid(index))%active(i(index), j(index), k(index))
                if (active_cell_index .eq. 0) then
                    
                        hydrogen_density(index) = 0
                    
                else
                    
                        hydrogen_density(index) = grid3D(index_of_grid(index))%Hden(active_cell_index)
                    
                end if
            END DO
        end if
        get_grid_hydrogen_density = 0
    END FUNCTION

    FUNCTION set_grid_hydrogen_density(i,j,k,hydrogen_density,index_of_grid, n)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(IN), DIMENSION(N) :: i,j,k, index_of_grid
        
        double precision, INTENT(IN), DIMENSION(N) :: hydrogen_density
        
        INTEGER set_grid_hydrogen_density
        INTEGER :: index, active_cell_index


        if (.not.is_grid_committed) then
            DO index = 1,n
                hydrogen_density_input(i(index), j(index), k(index)) = hydrogen_density(index)
            END DO
        else
            DO index = 1,n
                
                IF (index_of_grid(index).GT.nGrids) THEN
                    set_grid_hydrogen_density = -1
                    return
                END IF
                active_cell_index = grid3D(index_of_grid(index))%active(i(index), j(index), k(index))
                if (active_cell_index .ne. 0) then
                    
                    grid3D(index_of_grid(index))%Hden(active_cell_index) = hydrogen_density(index)
                    
                end if
            END DO
        end if
        
        set_grid_hydrogen_density = 0
    END FUNCTION



    FUNCTION get_grid_electron_density(i,j,k,index_of_grid,density, n)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(IN), DIMENSION(N) :: i,j,k, index_of_grid
        
        double precision, INTENT(OUT), DIMENSION(N) :: density
        
        INTEGER get_grid_electron_density
        INTEGER :: index, active_cell_index

        if (.not.is_grid_committed) then
            DO index = 1,n
                density(index) = 0.0
            END DO
            get_grid_electron_density = -1
            return
        else
        
            DO index = 1,n
                
                IF (index_of_grid(index).GT.nGrids) THEN
                    get_grid_electron_density = -1
                    return
                END IF
                active_cell_index = grid3D(index_of_grid(index))%active(i(index), j(index), k(index))
                if (active_cell_index .eq. 0) then
                    
                        density(index) = 0.0
                    
                else
                    
                        density(index) = grid3D(index_of_grid(index))%Ne(active_cell_index)
                    
                end if
            END DO
        end if
        get_grid_electron_density = 0
    END FUNCTION

    FUNCTION set_grid_electron_density(i,j,k,density,index_of_grid, n)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(IN), DIMENSION(N) :: i,j,k, index_of_grid
        
        double precision, INTENT(IN), DIMENSION(N) :: density
        
        INTEGER set_grid_electron_density
        INTEGER :: index, active_cell_index


        if (is_grid_committed) then
            DO index = 1,n
                
                IF (index_of_grid(index).GT.nGrids) THEN
                    set_grid_electron_density = -1
                    return
                END IF
                active_cell_index = grid3D(index_of_grid(index))%active(i(index), j(index), k(index))
                if (active_cell_index .ne. 0) then
                    
                    grid3D(index_of_grid(index))%Ne(active_cell_index) = density(index)
                    
                end if
            END DO
        else
            set_grid_electron_density = -1
            return
        end if
        
        set_grid_electron_density = 0
    END FUNCTION

    
    

    FUNCTION get_grid_ion_density(i,j,k,element,bin, index_of_grid, density, n)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(IN), DIMENSION(N) :: i,j,k,element,bin, index_of_grid
        
        double precision, INTENT(OUT), DIMENSION(N) :: density
        
        INTEGER get_grid_ion_density
        INTEGER :: index, active_cell_index

        if (.not.is_grid_committed) then
            DO index = 1,n
                density(index) = 0.0
            END DO
            get_grid_ion_density = -1
            return
        else
        
            DO index = 1,n
                
                IF (index_of_grid(index).GT.nGrids) THEN
                    get_grid_ion_density = -1
                    return
                END IF
                active_cell_index = grid3D(index_of_grid(index))%active(i(index), j(index), k(index))
                if (active_cell_index .eq. 0) then
                    
                        density(index) = 0.0
                    
                else
                        
                        density(index) = grid3D(index_of_grid(index))%ionDen(active_cell_index, element(index), bin(index))
                    
                end if
            END DO
        end if
        get_grid_ion_density = 0
    END FUNCTION

    FUNCTION set_grid_ion_density(i,j,k,element, bin,density,index_of_grid, n)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(IN), DIMENSION(N) :: i,j,k,element, bin, index_of_grid
        
        double precision, INTENT(IN), DIMENSION(N) :: density
        
        INTEGER set_grid_ion_density
        INTEGER :: index, active_cell_index


        if (is_grid_committed) then
            DO index = 1,n
                
                IF (index_of_grid(index).GT.nGrids) THEN
                    set_grid_ion_density = -1
                    return
                END IF
                active_cell_index = grid3D(index_of_grid(index))%active(i(index), j(index), k(index))
                if (active_cell_index .ne. 0) then
                    
                    grid3D(index_of_grid(index))%ionDen(active_cell_index, element(index), bin(index)) = density(index)
                    
                end if
            END DO
        else
            set_grid_ion_density = -1
            return
        end if
        
        set_grid_ion_density = 0
    END FUNCTION
    
    
    
    FUNCTION get_max_indices(index_of_grid,ni,nj,nk)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: index_of_grid
        INTEGER, INTENT(OUT) :: ni,nj,nk
        
        INTEGER :: get_max_indices

        IF (index_of_grid .GT. nGrids) THEN
            ni = 0
            nj = 0
            nk = 0
            get_max_indices = -1
            CONTINUE
        END IF
        
        ni = size(grid3D(index_of_grid)%xAxis)
        nj = size(grid3D(index_of_grid)%yAxis)
        nk = size(grid3D(index_of_grid)%zAxis)
        
        get_max_indices = 0
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
    
    
    
    function setup_mother_grid(grid)
        implicit none

        type(grid_type), intent(inout) :: grid      ! the grid

        ! local variables
        real                           :: denominator   ! denominator
        real                           :: dV           ! volume element
        real                           :: gasCell      ! mass of gas in current cell
        real                           :: H0in         ! estimated H0 at the inner radius for regionI
        real, pointer                  :: MdMg(:,:,:)  ! Md/Mg
        real                           :: MhMg         ! mass oh hydrogen over mass of gas
        real                           :: norm, scale  ! normalisation and scaling for meanField
        real                           :: radius       ! distance from the origin
        real                           :: totalMass    ! total ionized mass 
        real                           :: totalVolume  ! total active volume
        real                      :: echoVolume, vol   ! just echo volume

        real, dimension(nElements) :: aWeight
        real, parameter :: amu = 1.66053e-24 ! [g]

        real, pointer                  :: NdustTemp(:,:,:) ! temporary dust number density arra
        real, pointer                  :: dustAbunIndexTemp(:,:,:) ! temporary dust abundance index array
        real, pointer                  :: twoDscaleJTemp(:)
        


        integer                        :: i,j,k        ! counters
        integer                        :: index        ! general index
        integer                        :: ios, err     ! I/O and allocation error status
        integer                        :: elem, ion    ! counters
        integer                        :: nspec, ai    ! counters
        integer                        :: nsp          ! pointer
        integer                        :: nu0P         ! 

        integer                        :: yTop, xPmap  ! 2D indeces
        integer                        :: setup_mother_grid
                                                       ! with one of the axes
        character(len=40)              :: keyword      ! character string readers
        

        print*, 'in setup_mother_grid'

        ! this is the mother grid
        grid%motherP = 0
        setup_mother_grid = -1
        
        if (lg2D) then
           yTop = 1
        else
           yTop = grid%ny
        end if                


        grid%active = 1
        
          ! set up dust data
          if (lgDust) then
              allocate(NdustTemp(1:grid%nx,1:grid%ny,1:grid%nz), stat = err)
              if (err /= 0) then
                 print*, "! setMotherGrid: can't allocate NdustTemp memory"
                 return
              end if

              NdustTemp = 0.              

              if (lgMultiDustChemistry) then
                 allocate(dustAbunIndexTemp(1:grid%nx,1:grid%ny,1:grid%nz), stat = err)
                 if (err /= 0) then
                    print*, "! setMotherGrid: can't allocate dustAbunIndexTemp memory"
                    return
                 end if
                 dustAbunIndexTemp = 0.              
              end if

              ! set grains mass density [g/cm^3]
!              allocate(dustComPoint(nDustComponents))
!              dustComPoint = 0
!              dustComPoint(1) = 1

!              nSpecies = 0
!              nSpeciesMax = 0
!print*, 'heer', ndustcomponents
!              do icomp = 1, nDustComponents
!print*, icomp
!                 close(13)
!                 open(file =   dustSpeciesFile(icomp), action="read",unit=13, &
!                      &position="rewind",status="old", iostat = ios)
!                 if (ios /= 0 ) then
!                    print*, "! setMotherGrid: can't open file ", dustSpeciesFile(icomp)
!                    stop
!                 end if
!                 read(13, *) nSpeciesPart(icomp)
!print*, nspeciespart(icomp)
!                 close(13)
!                 nSpecies = nSpecies+nSpeciesPart(icomp)
!print*, nspecies
!                 if (nSpeciesMax < nSpeciesPart(icomp)) nSpeciesMax = nSpeciesPart(icomp)
!print*, nspeciesmax
!              end do
!
!              allocate(rho(1:nSpecies), stat = err)
!              if (err /= 0) then
!                 print*, "! setMotherGrid: can't allocate rho memory"
!                 stop
!              end if
!              rho=0.
!              allocate(grainVn(1:nSpecies), stat = err)
!              if (err /= 0) then
!                 print*, "! setMotherGrid: can't allocate grainVn memory"
!                 stop
!              end if
!              grainVn=0.
!              allocate(MsurfAtom(1:nSpecies), stat = err)
!              if (err /= 0) then
!                 print*, "! setMotherGrid: can't allocate surfAtom memory"
!                 stop
!              end if
!              MsurfAtom=0
!              
!              do icomp = 1, nDustComponents
!                 if (icomp > 1) dustComPoint(icomp) = dustComPoint(icomp-1)+nSpeciesPart(icomp)
!                 close(13)
!                 open(file =   dustSpeciesFile(icomp), action="read",unit=13, &
!                      & position="rewind",status="old", iostat = ios)
!                 if (ios /= 0 ) then
!                    print*, "! setMotherGrid: can't open file ", dustSpeciesFile(icomp)
!                    stop
!                 end if
!                 read(13, *) nSpeciesPart(icomp)              
!
!                 do i = 1, nSpeciesPart(icomp)
!                    read(13,*) extFile
!                    close(14)
!                    open(file=extFile,unit=14,  action="read", position="rewind",status="old", iostat = ios)
!                    if (ios /= 0 ) then
!                       print*, "! setMotherGrid: can't open file ", extFile
!                       stop
!                    end if
!                    read(14,*) readChar
!                    read(14,*) readChar, readReal, rho(dustComPoint(icomp)+i-1), grainVn(dustComPoint(icomp)+i-1), &
!                         &MsurfAtom(dustComPoint(icomp)+i-1)
!                    close(14)
!                 end do
!                 close(13)
!              end do
!
              if (lgMdMg .or. lgMdMh) then

                 allocate(MdMg(1:grid%nx,1:grid%ny,1:grid%nz), stat = err)
                 if (err /= 0) then
                    print*, "! setMotherGrid: can't allocate MdMg memory"
                    return
                 end if
                 MdMg = 0.
                 
                 if (lgDustConstant) then
                    MdMg = MdMgValue
                 else
                    close(20)
                    open(unit=20, file=MdMgFile,  action="read", position="rewind",status="old", iostat = ios)
                    if (ios /= 0 ) then
                       print*, "! setMotherGrid: can't open MdMgFile file ", MdMgFile
                       return
                    end if
                    read(20,*)keyword
                    if (keyword .ne. "#") backspace 20
                    
                    if (lgMultiDustChemistry) then
                       do i = 1, grid%nx
                          do j = 1, yTop
                             do k = 1, grid%nz
                                read(20, *) index, index, index, MdMg(i,j,k), dustAbunIndexTemp(i,j,k)
                             end do
                          end do
                       end do
                    else                    
                       do i = 1, grid%nx
                          do j = 1, yTop
                             do k = 1, grid%nz
                                read(20, *) index, index, index, MdMg(i,j,k)
                             end do
                          end do
                       end do
                    end if
                    close(20)
                 end if

              else

                 ! Ndust was directly defined by  the user
                 if (lgDustConstant) then
                    NdustTemp = NdustValue
                 else
                    close(20)
                    open(unit=20, file=NdustFile,  action="read", position="rewind",status="old", iostat = ios)
                    if (ios /= 0 ) then
                       print*, "! setMotherGrid: can't open NdustFile file ", NdustFile
                       return
                    end if
                    read(20,*)keyword
                    if (keyword .ne. "#") backspace 20

                    if (lgMultiDustChemistry) then
                       do i = 1, grid%nx
                          do j = 1, yTop
                             do k = 1, grid%nz
                                read(20, *) grid%xAxis(i), grid%yAxis(j), &
                                     & grid%zAxis(k), NdustTemp(i,j,k), dustAbunIndexTemp(i,j,k)
                             end do
                          end do
                       end do
                    else
                       print*,grid%nx,yTop,grid%nz
                       do i = 1, grid%nx
                          do j = 1, yTop
                             do k = 1, grid%nz
                                read(20, *) grid%xAxis(i), grid%yAxis(j), grid%zAxis(k), NdustTemp(i,j,k)
!                                write(6,'(3i4,4es11.3)')i,j,k,grid%xAxis(i), grid%yAxis(j), grid%zAxis(k), NdustTemp(i,j,k)
                             end do
                          end do
                       end do
                    end if
                    close(20)
                 end if

              end if              

           end if
           
        if (lg2D) grid%yAxis = grid%xAxis             


        ! set active cells pointers
        grid%nCells = 0
        do i = 1, grid%nx
            do j = 1, yTop
                do k = 1, grid%nz

                   if (lgDust .and. lgGas) then
                      if (hydrogen_density_input(i,j,k) > 0. .or. NdustTemp(i,j,k)>0.) then
                         grid%nCells = grid%nCells + 1
                         grid%active(i,j,k) = grid%nCells
                      else
                         grid%active(i,j,k) = 0
                         hydrogen_density_input(i,j,k) = 0.
                         NdustTemp(i,j,k) = 0.
                         if (lgMultiDustChemistry) dustAbunIndexTemp(i,j,k) = 0.
                      end if
                   else if ( lgDust .and. (.not.lgGas) ) then
                      if (NdustTemp(i,j,k)>0.) then
                         grid%nCells = grid%nCells + 1
                         grid%active(i,j,k) = grid%nCells
                      else
                         grid%active(i,j,k) = 0
                         NdustTemp(i,j,k) = 0.
                         if (lgMultiDustChemistry) dustAbunIndexTemp(i,j,k) = 0.
                      end if
                   else if ( (.not.lgDust) .and. lgGas) then 
                      if (hydrogen_density_input(i,j,k) > 0.) then
                         grid%nCells = grid%nCells + 1
                         grid%active(i,j,k) = grid%nCells
                      else
                         grid%active(i,j,k) = 0
                         hydrogen_density_input(i,j,k) = 0.
                      end if
                   else

                      print*, '! setup_mother_grid: no gas and no dust? The grid is empty.'
                      setup_mother_grid = -1
                      return
                   end if
                end do
            end do
        end do


        allocate(TwoDscaleJtemp(grid%nCells))
        TwoDscaleJtemp = 1.


        if (lg2D) then
             do i = 1, grid%nx
                do j = 2, grid%ny
                   do k = 1, grid%nz
                      radius = 1.e10*sqrt( (grid%xAxis(i)/1.e10)*&
                           &(grid%xAxis(i)/1.e10) + &
                           &(grid%yAxis(j)/1.e10)*(grid%yAxis(j)/1.e10) ) 

                      call locate(grid%xAxis, radius, xPmap)
                      if (xPmap < grid%nx) then
                         if (radius >= (grid%xAxis(xPmap)+grid%xAxis(xPmap+1))/2.) &
                              & xPmap = xPmap+1
                      end if
                      grid%active(i,j,k) = grid%active(xPmap, 1, k)
                      
                      if (grid%active(xPmap,1,k)>0) &
                           & TwoDScaleJtemp(grid%active(xPmap,1,k)) = &
                           & TwoDScaleJtemp(grid%active(xPmap,1,k))+1.
                      
                   end do
                end do
             end do

             grid%nCells = 0
             do i = 1,  grid%nx
                do k = 1,  grid%nz
                   if (grid%active(i,1,k) > 0) grid%nCells = grid%nCells +1                   

                end do
             end do

             allocate(TwoDscaleJ(grid%nCells))
             do i = 1, grid%nCells
                TwoDscaleJ(i) = TwoDscaleJtemp(i)
             end do
             deallocate(TwoDscaleJtemp)

          end if



        print*, '! setup_mother_grid: active cells :', grid%nCells



        ! allocate grid arrays
        if (lgGas .and. lgDust .and. lgPhotoelectric) then
            allocate(grid%JPEots(1:grid%nCells, 1:nbins), stat = err)
            if (err /= 0) then
                print*, "! setup_mother_grid: can't allocate JPEots memory"
                setup_mother_grid = -1
                return
            end if
            grid%JPEots=0
        end if

        if (lgGas) then

            allocate(grid%Hden(0:grid%nCells), stat = err)
            if (err /= 0) then
                print*, "! setMotherGrid: can't allocate grid memory"
                
                setup_mother_grid = -1
                return
            
            end if

            allocate(grid%recPDF(0:grid%nCells, 1:nbins), stat = err)
            if (err /= 0) then
                print*, "Can't allocate grid memory, 8"
                
                setup_mother_grid = -1
                return
            end if

            allocate(grid%totalLines(0:grid%nCells), stat = err)
            if (err /= 0) then
                print*, "Can't allocate grid memory, 10"
                
                setup_mother_grid = -1
                return
            end if

            allocate(grid%ionDen(0:grid%nCells, 1:nElementsUsed, 1:nstages), stat = err)
            if (err /= 0) then
                print*, "! setMotherGrid: can't allocate grid memory,ionDen"  
                
                setup_mother_grid = -1
                return
            end if
            allocate(ionDenUsed(1:nElementsUsed, 1:nstages), stat = err)
            if (err /= 0) then
                print*, "! setMotherGrid: can't allocate grid memory,ionDen"  
                
                setup_mother_grid = -1
                return
            end if
            allocate(grid%Ne(0:grid%nCells), stat = err)
            if (err /= 0) then
                print*, "! setMotherGrid: can't allocate grid memory"
                
                setup_mother_grid = -1
                return
            end if
            allocate(grid%Te(0:grid%nCells), stat = err)
            if (err /= 0) then
                print*, "! setMotherGrid:can't allocate grid memory"
                
                setup_mother_grid = -1
                return
            end if

            grid%Hden = 0.
            grid%Ne = 0.
            grid%Te = 0.        
            grid%ionDen = 0.
            grid%recPDF = 0.
            grid%totalLines = 0.
            
        end if

          if (Ldiffuse>0.) then
             allocate(grid%LdiffuseLoc(0:grid%nCells), stat = err)
             if (err /= 0) then
                print*, "! setMotherGrid: can't allocate grid memory : LdiffuseLoc "
                
                setup_mother_grid = -1
                return
             end if
             grid%LdiffuseLoc=0.
          end if


          allocate(grid%opacity(0:grid%nCells, 1:nbins), stat = err)
          if (err /= 0) then
             print*, "! setMotherGrid: can't allocate grid memory : opacity "
             
                setup_mother_grid = -1
                return
          end if
          allocate(grid%Jste(0:grid%nCells, 1:nbins), stat = err)
          if (err /= 0) then
             print*, "! setMotherGrid: can't allocate grid memory : Jste"
             
                setup_mother_grid = -1
                return
          end if

          allocate(grid%escapedPackets(0:grid%nCells, 0:nbins,0:nAngleBins), stat = err)
          if (err /= 0) then
             print*, "! setMotherGrid: can't allocate grid memory : Jste"
             
                setup_mother_grid = -1
                return
          end if

          if (lgEquivalentTau) then

             allocate(SEDnoExt(1:nbins), stat = err)
             if (err /= 0) then
                print*, "! setMotherGrid: can't allocate grid memory : SEDnoExt"
                
                setup_mother_grid = -1
                return
             end if
             allocate(equivalentTau(1:nbins), stat = err)
             if (err /= 0) then
                print*, "! setMotherGrid: can't allocate grid memory : equivalentTau"
                
                setup_mother_grid = -1
                return
             end if

          end if
          if (lgDust) then

             allocate(grid%Ndust(0:grid%nCells), stat = err)
             if (err /= 0) then
                print*, "! grid: can't allocate Ndust memory"
                
                setup_mother_grid = -1
                return
             end if
             grid%Ndust=0.
             ! be 2.02.44
             allocate(grid%dustAbunIndex(0:grid%nCells), stat = err)
             if (err /= 0) then
                print*, "! grid: can't allocate dustAbunIndex memory"
                
                setup_mother_grid = -1
                return
             end if
             grid%dustAbunIndex=0.
             ! be 2.02.44 end
             
             if (.not.lgGas) then
                allocate(grid%dustPDF(0:grid%nCells, 1:nbins), stat = err)
                if (err /= 0) then
                   print*, "! grid: can't allocate dustPDF memory"
                   
                setup_mother_grid = -1
                return
                end if
                grid%dustPDF = 0.
             end if

          end if
!          print* ,'a'

!BS10          if (lgDebug) then
             allocate(grid%Jdif(0:grid%nCells, 1:nbins), stat = err)
             if (err /= 0) then
                print*, "! grid: can't allocate grid memory"
                
                setup_mother_grid = -1
                return
             end if

             grid%Jdif = 0. 

             ! allocate pointers depending on nLines
             if (nLines > 0) then

                allocate(grid%linePackets(0:grid%nCells, 1:nLines), stat = err)
                if (err /= 0) then
                   print*, "! initCartesianGrid: can't allocate grid%linePackets memory"
                   
                setup_mother_grid = -1
                return
                end if

                allocate(grid%linePDF(0:grid%nCells, 1:nLines), stat = err)
                if (err /= 0) then
                   print*, "! initCartesianGrid: can't allocate grid%linePDF memory"
                   
                setup_mother_grid = -1
                return
                end if
                
                grid%linePackets = 0.
                grid%linePDF     = 0.
                
             end if
             
!BS10          end if

          allocate(grid%lgConverged(0:grid%nCells), stat = err)
          if (err /= 0) then
             print*, "Can't allocate memory to lgConverged array"
             
                setup_mother_grid = -1
                return
          end if

          allocate(grid%lgBlack(0:grid%nCells), stat = err)
          if (err /= 0) then
             print*, "Can't allocate memory to lgBlack array"
             
                setup_mother_grid = -1
                return
          end if


          if (lgNeInput) then
             allocate(grid%NeInput(0:grid%nCells), stat = err)
             if (err /= 0) then
                print*, "! setMotherGrid: can't allocate grid memory, grid%NeInput"
                 setup_mother_grid = -1
                return
             end if
             grid%NeInput = 0.
          end if


          grid%opacity = 0.        
          grid%Jste = 0.        
          grid%lgConverged = 0
          grid%lgBlack = 0           

          do i = 1, grid%nx
             do j = 1, yTop
                do k = 1, grid%nz
                   if (grid%active(i,j,k)>0) then

                      if (lgGas) then
                         grid%Hden(grid%active(i,j,k)) = hydrogen_density_input(i,j,k)
                         grid%Te(grid%active(i,j,k)) = TeStart
                      end if
                      if (lgDust) then
                         grid%Ndust(grid%active(i,j,k)) = NdustTemp(i,j,k)
                         if (lgMultiDustChemistry) &
                              & grid%dustAbunIndex(grid%active(i,j,k)) = &
                              &dustAbunIndexTemp(i,j,k)
                      end if
                   end if
                end do
             end do
          end do
!          print*, 'b'
          if (lgNeInput) then 
             grid%NeInput = grid%Hden
             ! 1.11 is just a starting guess for the ionization 
             ! factor
             grid%Hden = grid%Hden/1.11
          end if

          if (lgGas) then

             H0in  = 1.e-5

             do i = 1, grid%nx
                do j = 1, yTop
                   do k = 1, grid%nz
                      
                      if (grid%active(i,j,k)>0 ) then

                         ! calculate ionDen for H0
                         if (lgElementOn(1)) then
                            grid%ionDen(grid%active(i,j,k),elementXref(1),1) = H0in                                              
                            grid%ionDen(grid%active(i,j,k),elementXref(1),2) = &
                                 & 1.-grid%ionDen(grid%active(i,j,k),elementXref(1),1)
                         end if
                         if (lgElementOn(2)) then
                            grid%ionDen(grid%active(i,j,k),elementXref(2),1) = &
                                 & grid%ionDen(grid%active(i,j,k),elementXref(1),1)
                            grid%ionDen(grid%active(i,j,k),elementXref(2),2) = &
                                 & (1.-grid%ionDen(grid%active(i,j,k),elementXref(2),1))
                            grid%ionDen(grid%active(i,j,k),elementXref(2),3) = 0.
                         end if

                         ! initialize Ne
                         grid%Ne(grid%active(i,j,k)) =  grid%Hden(grid%active(i,j,k))
                         
                         ! initialize all heavy ions (careful that the sum over all ionisation 
                         ! stages of a given atom doesn't exceed 1.)
                         do elem = 3, nElements
                            do ion = 1, min(elem+1,nstages)
                               if (lgElementOn(elem)) then
                                  if (ion == 1) then
                                     grid%ionDen(grid%active(i,j,k),elementXref(elem),ion) &
                                          & = grid%ionDen(grid%active(i,j,k),1,1)
                                  else
                                     grid%ionDen(grid%active(i,j,k),elementXref(elem),ion) = 0.
                                  end if
                               end if
                            end do
                         end do
                      end if     ! active condition
                   
                   end do
                end do
             end do
             

          end if ! lgGas

           if (lgDust) then
              if(associated(NdustTemp)) deallocate(NdustTemp)
              if(lgMultiDustChemistry .and. associated(dustAbunIndexTemp)) deallocate(dustAbunIndexTemp)
           end if

           ! set up atomic weight array
           aWeight = (/1.0080, 4.0026, 6.941, 9.0122, 10.811, 12.0111, 14.0067, 15.9994, &
                & 18.9984, 20.179, 22.9898, 24.305, 26.9815, 28.086, 30.9738, 32.06, 35.453, &
                & 39.948, 39.102, 40.08, 44.956, 47.90, 50.9414, 51.996, 54.9380, 55.847, 58.9332, &
                & 58.71, 63.546, 65.37 /)

           totalDustMass = 0.
           totalMass = 0.
           totalVolume = 0.
           echoVolume = 0.

           do i = 1, grid%nx
              do j = 1, yTop
                 do k = 1, grid%nz
                    grid%echoVol(i,j,k)=0. ! initialize
                    if (grid%active(i,j,k)>0) then

                       dV = getVolume(grid,i,j,k)

                       if (lgGas) then
                          gasCell = 0.
                          do elem = 1, nElements
                             gasCell = gasCell + grid%elemAbun(grid%abFileIndex(i,j,k),elem)*&
                                  & aWeight(elem)*amu
                             totalMass = totalMass + &
                                  & grid%Hden(grid%active(i,j,k))*dV*grid%elemAbun(grid%abFileIndex(i,j,k),elem)*&
                                  & aWeight(elem)*amu
                          end do
                       end if

                       totalVolume = totalVolume + dV
!
! echoes only
!
                       if (lgEcho) then
                          grid%echoVol(i,j,k)=vEcho(grid,i,j,k,echot1,echot2,vol)
                          echoVolume = echoVolume + vol 
                       endif

                       if (lgDust .and. (lgMdMg.or.lgMdMh) ) then

                          
                          if (lgMdMh) then
                             MhMg=0.
                             do elem = 1, nElements
                                ! transform to MdMg
                                MhMg = MhMg+grid%elemAbun(grid%abFileIndex(i,j,k),elem)*&
                                     & aWeight(elem)
                             end do
                             MhMg = 1./MhMg                             
                             MdMg(i,j,k) = MdMg(i,j,k)*MhMg
                          end if

                          if (.not.lgGas) then
                             print*, '! setMotherGrid: Mass to dust ratio (MdMg) cannot be used in a pure dust (noGas)&
                                  & simulation. Ndust must be used instead.'
                             stop
                          end if

                          grid%Ndust(grid%active(i,j,k)) = gasCell*MdMg(i,j,k)*grid%Hden(grid%active(i,j,k))

                          denominator = 0.

                          if (lgMultiDustChemistry) then
                             nsp = grid%dustAbunIndex(grid%active(i,j,k))
                          else
                             nsp = 1
                          end if
!print*, nsp, 'here!'
                          do nspec = 1, nSpeciesPart(nsp)
                             do ai = 1, nSizes
                                denominator = denominator + &
                                     & (1.3333*Pi*( (grainRadius(ai)*1.e-4)**3)*&
                                     & rho(dustComPoint(nsp)+nspec-1)*grainWeight(ai)*&
                                     & grainAbun(nsp, nspec))
                             end do
                          end do
                          grid%Ndust(grid%active(i,j,k)) = grid%Ndust(grid%active(i,j,k))/&
                               & (denominator)
                       end if

                       ! calculate total dust mass
                       if (lgDust) then
                          if (lgMultiDustChemistry) then
                             nsp = grid%dustAbunIndex(grid%active(i,j,k))
                          else
                             nsp = 1
                          end if

                          do ai = 1, nsizes
                             do nspec = 1, nspeciesPart(nsp)
                                totalDustMass = totalDustMass + &
                                     &(1.3333*Pi*((grainRadius(ai)*1.e-4)**3)*&
                                     & rho(dustComPoint(nsp)-1+nspec)*grainWeight(ai)*&
                                     & grainAbun(nsp,nspec))*grid%Ndust(grid%active(i,j,k))*dV

                             end do
                          end do

                       end if

                    end if
                 end do
              end do
           end do



           if(associated(MdMg)) deallocate(MdMg)

           if (taskid == 0) then

              print*, 'Mothergrid :'
              if (lgGas) then
                 print*, 'Total gas mass of ionized region by mass [1.e45 g]: ', totalMass
              end if
              if (lgDust) then
                 print*, 'Total dust mass of ionized region by mass [1.e45 g]: ', totalDustMass
              end if
              print*, 'Total volume of the active region [e45 cm^3]: ', totalVolume
              if (lgEcho) then 
                 print*, 'Total volume of the echo [e45 cm^3]: ', echoVolume*846732407.," or ",echoVolume," ly^3"
                 open(unit=99, status='unknown', position='rewind', file='output/echo.out', action="write",iostat=ios)
                 write(99,*)'Total volume of the active region [e45 cm^3]: ', totalVolume
                 write(99,*)'Total volume of the echo [e45 cm^3]: ', echoVolume*846732407.," or ",echoVolume," ly^3"
                 if (echoVolume .eq. 0.) then
                    print*,'No dust in echo region. Stopping'
                    write(99,*)'No dust in echo region. Stopping'
                    stop
                 endif
                 close(99)
              endif

           end if
              
           ! if we are using a plane parallel ionization then we must find the luminosity 
           ! of the ionizing plane from the input meanField
           if (lgPlaneIonization) then

              print*, 'Flux above ', nu0, ' is ', meanField

              if (nu0 > 0.) then
                 call locate(nuArray, nu0, nu0P) 
                 if (nu0P >= nbins .or. nu0P <1) then
                    print*, "! setMotherGrid: insanity in nu0P", nu0P, nuArray(i), nuArray(nbins)
                    stop
                 end if
                 norm = 0.
                 do i = nu0P, nbins
                    norm = norm+inSpectrumPhot(i)*widflx(i)
                 end do
                 scale  = meanField/norm
                 norm = 0.
                 do i = 1, nbins
                    norm = norm+inSpectrumErg(i)*widFlx(i)
                 end do
                 meanField = norm*scale
              end if

              print*, 'Flux bolometric is ', meanField              

              Lstar(1) = (meanField/1.e36)*grid%xAxis(grid%nx)*grid%zAxis(grid%nz)
              deltaE(1) = Lstar(1)/nPhotons(1)
              
              ! AVE (17-nov-2010) moved this code into the if statement
              ! Lstar is set just before this statement
              if (taskid == 0) then
                print*, 'Total ionizing flux :', Lstar(1)
                print*, 'deltaE :', deltaE(1)
              end if
           end if
                 
            setup_mother_grid = 0
           print*, 'out setup_mother_grid'

    end function setup_mother_grid
    
    


    
    function step_monte_carlo(grid)
        implicit none
    
        include 'mpif.h'

        type(grid_type), intent(inout) :: grid(*)
        integer :: step_monte_carlo
        integer :: i,j,k
        ! local variables
        real, pointer :: budgetTemp(:,:)      ! temporary dust heating budget array
        real, pointer :: dustPDFTemp(:,:)     ! temporary dust emission PDF array
        real, pointer :: escapedPacketsTemp(:,:,:)!temporary escaped packets array
        real, pointer :: fEscapeResPhotonsTemp(:,:) ! temporary escape res line phot 
        real, pointer :: JDifTemp(:,:)        ! temporary diffuse field array
        real, pointer :: JSteTemp(:,:)        ! temporary stellar field array
        real, pointer :: linePacketsTemp(:,:) ! temporary line packets array
        real, pointer :: opacityTemp(:,:)     ! temporary opacities array
        real, pointer :: recPDFTemp(:,:)      ! temporary rec prob distribution function
        real, pointer :: linePDFTemp(:,:)     ! temporary line prob distribution function
        real, pointer :: totalLinesTemp(:)    ! temporary fraction of non-ionizing line phots

        integer, pointer       :: planeIonDistributionTemp(:,:) 
        integer, pointer       :: resLinePacketsTemp(:) ! temporary array for extra packets
        integer                :: err             ! allocation error status
        integer                :: freq,nS! counters
        integer                :: ifreq, ian      ! counters
        integer                :: ios,iG          ! I/O error status           
        integer                :: load,rest       ! 
        integer                :: size            ! size for mpi
        integer                :: iCell           ! cell index including non-active
        integer                :: iStar           ! star index
        integer                :: ai              ! grain size counter
        integer                :: imu             ! direction cosine            
        integer                :: cellLoc(3)      ! local cell counters
        integer                :: gpLoc           ! local grid counter
        
        integer                :: ii,jj,kk        ! counters
        integer                :: ngridloc
        integer                :: yTop

        integer                :: dcp, nsp
         
        step_monte_carlo = 0

        ! re-initialize MC estimators
        do iG = 1, nGrids
          grid(iG)%lgConverged(0:grid(iG)%nCells)    = 0
          grid(iG)%lgBlack(0:grid(iG)%nCells)        = 0
          if (lgGas) then
             ! zero out PDF arrays
             grid(iG)%recPDF(0:grid(iG)%nCells, 1:nbins) = 0.
             if (lgDebug) grid(iG)%linePDF(0:grid(iG)%nCells, 1:nLines)    = 0.
             grid(iG)%totalLines(0:grid(iG)%nCells) = 0.  
             
             ! zero out Balmer jump
             BjumpTemp = 0.
             Bjump     = 0.
          end if

          if (lgDust .and. .not.lgGas) then
             ! zero out dust PDF arrays
             grid(iG)%dustPDF(0:grid(iG)%nCells, 1:nbins)     = 0.
          end if
        end do

!*****************************************************************************

        iCell = 0

        do iG = 1, nGrids

          if (ig>1 .or. (.not.lg2D)) then
             yTop = grid(iG)%ny
          else if (iG ==1 .and. lg2D) then
             yTop = 1
          end if


          grid(iG)%opacity(0:grid(iG)%nCells, 1:nbins) = 0.

          if (lgGas) then
             if (taskid==0) print*, '! iterateMC: ionizationDriver in', iG
             ! calculate the opacities at every grid cell
             do i = 1, grid(iG)%nx
                do j = 1, yTop
                   do k = 1, grid(iG)%nz
                      iCell = iCell+1
                      if (mod(iCell-(taskid+1),numtasks)==0) &
                           & call ionizationDriver(grid(iG),i,j,k)
                   end do
                end do
             end do
             if (taskid==0) print*, '! iterateMC: ionizationDriver out', iG

             allocate(opacityTemp(0:grid(iG)%nCells, nbins), stat = err)
             if (err /= 0) then
                print*, "! iterateMC: can't allocate array memory: opacityTemp ", iG
                stop
             end if
             
             opacityTemp = 0.
             
             size = (grid(iG)%nCells+1)*nbins
             
             call mpi_allreduce(grid(iG)%opacity, opacityTemp, size, &
                  & mpi_real, mpi_sum, mpi_comm_world, ierr)
             
             
             do i = 1, grid(iG)%nx
                do j = 1, yTop
                   do k = 1, grid(iG)%nz
                      do freq = 1, nbins
                         
                         if (grid(iG)%active(i,j,k)>0) then
                            grid(iG)%opacity(grid(iG)%active(i,j,k),freq) = opacityTemp(grid(iG)%active(i,j,k),freq)
                         end if

                      end do
                   end do
                end do
             end do
             if ( associated(opacityTemp) ) deallocate(opacityTemp)

          end if

          ! add dust contribution to total opacity
          if (taskid==0) print*, '! iterateMC: adding dust contribution to total opacity ',iG           
          if (lgDust) then
             grid(iG)%scaOpac(0:grid(iG)%nCells, 1:nbins) = 0.
             grid(iG)%absOpac(0:grid(iG)%nCells, 1:nbins) = 0.

             if((nIterateMC==1 .and. lgEquivalentTau)) then
                grid(iG)%scaOpac(0:grid(iG)%nCells, 1:nbins) = 0.
                grid(iG)%absOpac(0:grid(iG)%nCells, 1:nbins) = 0.
             else


                if (ig>1 .or. (.not. lg2D)) then
                   yTop = grid(iG)%ny
                else if (iG ==1 .and. lg2D) then
                   yTop = 1
                end if

                do i = 1, grid(iG)%nx
                   do j = 1, yTop
                      do k = 1, grid(iG)%nz
                         if (grid(iG)%active(i,j,k)>0) then
                            if (lgMultiDustChemistry) then 
                               dcp = dustComPoint(grid(iG)%dustAbunIndex(grid(iG)%active(i,j,k)))
                               nsp = grid(iG)%dustAbunIndex(grid(iG)%active(i,j,k))
                            else
                               dcp = dustComPoint(1)
                               nsp = 1
                            end if

                            do nS = 1, nSpeciesPart(nsp)

                               do ai = 1, nSizes                                
                                  if (grid(iG)%Tdust(nS,ai,grid(iG)%active(i,j,k))<TdustSublime(dcp-1+nS)) then
                                     do freq = 1, nbins 
                                        grid(iG)%scaOpac(grid(iG)%active(i,j,k),freq) = &
                                             & grid(iG)%scaOpac(grid(iG)%active(i,j,k),freq) + & 
                                             & grainAbun(nsp,nS)*&
                                             & grainWeight(ai)*grid(iG)%Ndust(grid(iG)%active(i,j,k))*&
                                             & xSecArray(dustScaXsecP(nS+dcp-1,ai)+freq-1)
                                        grid(iG)%absOpac(grid(iG)%active(i,j,k),freq) = &
                                             & grid(iG)%absOpac(grid(iG)%active(i,j,k),freq) + & 
                                             & grainAbun(nsp,nS)&
                                             & *grainWeight(ai)*grid(iG)%Ndust(grid(iG)%active(i,j,k))*&
                                             & xSecArray(dustAbsXsecP(nS+dcp-1,ai)+freq-1)
                                     end do
                                  end if
                               end do
                            end do
                            
                            do freq = 1, nbins
                               grid(iG)%opacity(grid(iG)%active(i,j,k),freq) = &
                                    &grid(iG)%opacity(grid(iG)%active(i,j,k),freq) + &
                                    & (grid(iG)%scaOpac(grid(iG)%active(i,j,k),freq) + &
                                    &grid(iG)%absOpac(grid(iG)%active(i,j,k),freq))
                            end do
                         end if
                         
                      end do
                   end do
                end do
             end if
             
          end if
          if (taskid==0) print*, '! iterateMC: dust contribution to total opacity added ',iG                            
          
        end do ! ngrids

        call mpi_barrier(mpi_comm_world, ierr)


!*****************************************************************************

           

           
            if (taskid==0 .and. lgWritePss) then
            open(unit=89, file=trim(outputDir)//'qHeatPss.out',  action="write",status='unknown', position='rewind', iostat=ios)
            if (ios /= 0) then
             print*, "! iterationMC: can't open file for writing, output/qHeatPss.out"
             stop
            end if
            end if

            ! check if min convergence was reached to carry out resonant line transfer
            if (lgDust .and. lgGas .and. convPercent>=resLinesTransfer .and. lgResLinesFirst &
            &.and. (.not.nIterateMC==1)) then

            call initResLines(grid(1:nGrids))

            do iG = 1, nGrids
             allocate(fEscapeResPhotonsTemp(0:grid(iG)%nCells, 1:nResLines), stat &
                  &= err)
             if (err /= 0) then
                print*, "! iterateMC: can't allocate array&
                     & memory:fEscapeResPhotonsTemp"
                stop
             end if
             fEscapeResPhotonsTemp = 0.

             allocate(resLinePacketsTemp(0:grid(iG)%nCells), stat &
                  &= err)
             if (err /= 0) then
                print*, "! iterateMC: can't allocate array&
                     & memory:resLinePacketsTemp"
                stop
             end if
             resLinePacketsTemp = 0.

            end do
            end if           

            ! set the diffuse PDFs at every grid cell
            do iG = 1, nGrids

            if (ig>1 .or. (.not.lg2D) ) then
             yTop = grid(iG)%ny
            else if (iG ==1 .and. lg2D) then
             yTop = 1
            end if


            if (taskid==0) print*, '! iterateMC: emissionDriver in',iG           
            iCell = 0
            do i = 1, grid(iG)%nx
             do j = 1, yTop
                do k = 1, grid(iG)%nz
                   iCell = iCell+1
                   if (mod(iCell-(taskid+1),numtasks)==0) &
                        & call emissionDriver(grid,i,j,k,iG)

                end do
             end do
            end do

            if (taskid==0) print*, '! iterateMC: emissionDriver out', iG           
            if (taskid==0 .and. lgWritePss) close(89)

            if (lgDust .and. .not.lgGas) then 
             allocate(dustPDFTemp(0:grid(iG)%nCells, 1:nbins),&
                  & stat = err)
             if (err /= 0) then
                print*, "! iterateMC: can't allocate array memory:&
                     & dustPDFTemp ", iG
                stop
             end if
             dustPDFTemp = 0.
            end if

            if (lgGas) then
             call mpi_allreduce(BjumpTemp, Bjump, 1, mpi_real&
                  &, mpi_sum, mpi_comm_world, ierr)
             
             ! Balmer jump is in [erg/s/A] since factor of e-40 from gas emissivity
             ! and factor of e45 from volume calculations, so multiply by e5
             if (lgSymmetricXYZ) Bjump = Bjump*8.
             Bjump = Bjump*1.e5
             
             print*, "Balmer Jump: [erg/s/A] ", Bjump
             
             allocate(recPDFTemp(0:grid(iG)%nCells, 1:nbins),&
                  & stat = err)
             if (err /= 0) then
                print*, "! iterateMC: can't allocate array memory:&
                     & recPDFTemp "
                stop
             end if

             allocate(totalLinesTemp(0:grid(iG)%nCells), stat =&
                  & err)
             if (err /= 0) then
                print*, "! iterateMC: can't allocate array memory:&
                     & opacityTemp "
                stop
             end if
             
             recPDFTemp = 0.
             totalLinesTemp = 0.
             
            end if

            if (lgDebug) then
             allocate(linePDFTemp(0:grid(iG)%nCells, nLines),&
                  & stat = err)
             if (err /= 0) then
                print*, "! iterateMC: can't allocate array memory:&
                     & opacityTemp "
                stop
             end if
             linePDFTemp = 0.
            end if

            size =  (grid(iG)%nCells+1)*nbins



            if (lgGas) then
             call mpi_allreduce(grid(iG)%recPDF, recPDFTemp, size, mpi_real&
                  &, mpi_sum, mpi_comm_world, ierr)

             if (lgDebug) then
                size =  (grid(iG)%nCells+1)*nLines

                call mpi_allreduce(grid(iG)%linePDF, linePDFTemp, size,&
                     & mpi_real, mpi_sum, mpi_comm_world, ierr)
             end if
            end if

            if (lgDust .and. .not.lgGas) then

             call mpi_allreduce(grid(iG)%dustPDF, dustPDFTemp, size, mpi_real&
                  &, mpi_sum, mpi_comm_world, ierr)

             do i = 0, grid(iG)%nCells
                do freq = 1, nbins
                   grid(iG)%dustPDF(i,freq) = dustPDFTemp(i,freq)
            !if (nIterateMC>1 .and. ig == 181 .and. i==9) then
            !print*,i, grid(iG)%dustPDF(i,freq)
            !end if
                end do
             end do

             if (associated(dustPDFTemp)) deallocate(dustPDFTemp)
            end if

            call mpi_barrier(mpi_comm_world, ierr)

            size = (grid(iG)%nCells+1)

            if (lgGas) then
             call mpi_allreduce(grid(iG)%totalLines, totalLinesTemp, size,&
                  & mpi_real, mpi_sum, mpi_comm_world, ierr)
             
             do i = 1, grid(iG)%nx
                do j = 1, yTop
                   do k = 1, grid(iG)%nz
                      if (grid(iG)%active(i,j,k)>0) then
                         grid(iG)%totalLines(grid(iG)%active(i,j,k)) = totalLinesTemp(grid(iG)%active(i,j,k))
                         do freq = 1, nbins
                            grid(iG)%recPDF(grid(iG)%active(i,j,k),freq) = recPDFTemp(grid(iG)%active(i,j,k) ,freq)
                         end do
                         if (lgDebug) then
                            do freq = 1, nLines
                               grid(iG)%linePDF(grid(iG)%active(i,j,k),freq) = linePDFTemp(grid(iG)%active(i,j,k),freq)
                            end do
                         end if
                      end if
                   end do
                end do
             end do

             call mpi_barrier(mpi_comm_world, ierr)
             
             if ( associated(totalLinesTemp) )&
                  & deallocate(totalLinesTemp)
             if (lgDebug) then
                if ( associated(linePDFTemp) ) deallocate(linePDFTemp)
             end if
             if ( associated(recPDFTemp) ) deallocate(recPDFTemp) 
            end if

            ! check if min convergence was reached to carry out resonant line transfer
            if (lgDust .and. lgGas .and. convPercent>=resLinesTransfer .and. lgResLinesFirst &
               &.and. (.not.nIterateMC==1)) then

             if (iG==nGrids) lgResLinesFirst = .false.
             
             size = (grid(iG)%nCells+1)*nResLines

             call mpi_allreduce(grid(iG)%fEscapeResPhotons, fEscapeResPhotonsTemp, size, &
                  & mpi_real, mpi_sum, mpi_comm_world, ierr)
             
             grid(iG)%fEscapeResPhotons(0:grid(iG)%nCells, 1:nResLines) = fEscapeResPhotonsTemp
             
             size = grid(iG)%nCells+1

             if (associated(fEscapeResPhotonsTemp)) deallocate(fEscapeResPhotonsTemp)

             call mpi_allreduce(grid(iG)%resLinePackets, resLinePacketsTemp, size, &
                  & mpi_real, mpi_sum, mpi_comm_world, ierr)

             grid(iG)%resLinePackets(0:grid(iG)%nCells) = resLinePacketsTemp

             if (associated(resLinePacketsTemp)) deallocate(resLinePacketsTemp)
             
             call mpi_barrier(mpi_comm_world, ierr)

            end if


            end do ! nGrids loop

            !**********************************************************************************************

            do iG = 1, nGrids

            grid(iG)%Jste(0:grid(iG)%nCells, 1:nbins)    = 0.
            if (lgDebug) then
             grid(iG)%Jdif(0:grid(iG)%nCells, 1:nbins) = 0.
             grid(iG)%linePackets(0:grid(iG)%nCells, 1:nLines) = 0.
            end if
            end do


            totalEscaped = 0.

            do iG = 1, nGrids
            grid(iG)%escapedPackets(0:grid(iG)%nCells, 0:nbins,0:nAngleBins) = 0.
            end do

            do iStar = 1, nStars
            if(taskid==0) print*, 'iterateMC: Starting transfer for ionising source ', iStar

            load = int(nPhotons(iStar)/numtasks)
            rest = mod(nPhotons(iStar), numtasks)           


            if (lgPlaneIonization) then
             planeIonDistribution = 0
            end if


            ! send the photons through and evaluate the MC 
            ! estimators of Jste and Jdif at every grid cell
            if (taskid < rest) then
             load = load+1
             call energyPacketDriver(iStar,load, grid(1:nGrids))
            else
             call energyPacketDriver(iStar,load, grid(1:nGrids))
            end if

            call mpi_barrier(mpi_comm_world, ierr)
            end do

            if (Ldiffuse>0.) then

            if (emittingGrid>0) then                 
             ngridloc = emittingGrid
            else
             ngridloc = ngrids
            end if

            do gpLoc = 1, nGridloc

             if (gploc>1 .or. (.not. lg2D)) then
                yTop = grid(gploc)%ny
             else if (gploc ==1 .and. lg2D) then
                yTop = 1
             end if


             if(taskid==0) print*, 'iterateMC: Starting transfer for diffuse source grid: ', gpLoc
             do ii = 1,grid(gpLoc)%nx
                do jj = 1,yTop
                   do kk = 1,grid(gpLoc)%nz
             
                      if (grid(gpLoc)%active(ii,jj,kk)>0) then
                      
                         cellLoc(1)  = ii
                         cellLoc(2)  = jj
                         cellLoc(3)  = kk

                         load = int(nPhotonsDiffuseLoc/numtasks)
                         rest = mod(nPhotonsDiffuseLoc, numtasks)           


                         ! send the photons through and evaluate the MC 
                         ! estimators of Jste and Jdif at every grid cell
                         if (taskid < rest) then
                            load = load+1
                            call energyPacketDriver(iStar=0,n=load, grid=grid(1:nGrids), &
                                 & gpLoc=gpLoc, cellLoc=cellLoc)
                         else
                            call energyPacketDriver(iStar=0,n=load, grid=grid(1:nGrids), &
                                 & gpLoc=gpLoc, cellLoc=cellLoc)
                         end if

                         call mpi_barrier(mpi_comm_world, ierr)

                      end if

                   end do
                end do
             end do
            end do

            end if


            if (lgPlaneIonization) then

            allocate(planeIonDistributionTemp(grid(1)%nx, grid(1)%nz), stat = err)
            if (err /= 0) then
             print*, "! iterateMC: can't allocate grid memory, planeIonDistributionTemp"
             stop
            end if
            planeIonDistributionTemp = 0

            size = grid(1)%nx*grid(1)%nz

            call mpi_allreduce(planeIonDistribution, planeIonDistributionTemp, size, &
               & mpi_integer, mpi_sum, mpi_comm_world, ierr)

            planeIonDistribution = planeIonDistributionTemp

            if (taskid ==0) then
             open(file=trim(outputDir)//"planeIonDistribution.out",  action="write",unit=18, status="unknown")
             do i = 1, grid(1)%nx
                do k = 1, grid(1)%nz
                   write(18,*) i,k,planeIonDistribution(i,k)
                end do
             end do                 
             close(18)
            end if

            if (associated(planeIonDistributionTemp)) deallocate(planeIonDistributionTemp)

            end if

            do iG = 1, nGrids

            if (ig>1 .or. (.not. lg2D)) then
             yTop = grid(iG)%ny
            else if (iG ==1 .and. lg2D) then
             yTop = 1
            end if



            allocate(escapedPacketsTemp(0:grid(iG)%nCells, 0:nbins, 0:nAngleBins), stat = err)
            if (err /= 0) then
             print*, "! iterateMC: can't allocate grid memory: escapedPacketsTemp", iG
             stop
            end if

            escapedPacketsTemp  = 0.

            allocate(JSteTemp(0:grid(iG)%nCells, nbins), stat = err)
            if (err /= 0) then
             print*, "! iterateMC: can't allocate array memory: JsteTemp ", iG
             stop
            end if
            JSteTemp           = 0.

            if (lgDebug) then
             allocate(JDifTemp(0:grid(iG)%nCells, nbins), stat = err)
             if (err /= 0) then
                print*, "! iterateMC: can't allocate array memory: JdifTemp "
                stop
             end if
             allocate(linePacketsTemp(0:grid(iG)%nCells, nLines), stat = err)
             if (err /= 0) then
                print*, "! iterateMC: can't allocate array memory: linePacketsTemp "
                stop
             end if

             JDifTemp           = 0.              
             linePacketsTemp    = 0.
            end if


            size =  (grid(iG)%nCells+1)*(1+nbins)*(nAngleBins+1)

            call mpi_allreduce(grid(iG)%escapedPackets, escapedPacketsTemp, size, &
               & mpi_real, mpi_sum, mpi_comm_world, ierr)

            do i = 0, grid(iG)%nCells
             do freq = 0, nbins                       
                do imu = 0, nAngleBins                       
                   grid(iG)%escapedPackets(i, freq,imu) = escapedPacketsTemp(i, freq, imu)

                end do
             end do
            end do
             

            if ( associated(escapedPacketsTemp) ) deallocate(escapedPacketsTemp)

            !              if (taskid==0) call writeSED(grid)
            !              if (taskid==0 .and. contCube(1)>0. .and. contCube(2)>0. ) &
            !                   & call writeContCube(grid, contCube(1),contCube(2))    

            size =  (grid(iG)%nCells+1)*nbins

            if (lgDebug) then
             call mpi_allreduce(grid(iG)%JDif, JDifTemp, size, &
                  & mpi_real, mpi_sum, mpi_comm_world, ierr)
            end if

            call mpi_allreduce(grid(iG)%JSte, JSteTemp, size, &
               & mpi_real, mpi_sum, mpi_comm_world, ierr)

            size =  (grid(iG)%nCells+1)*nLines

            if (lgDebug) then
             call mpi_allreduce(grid(iG)%linePackets, linePacketsTemp, size, &
                  & mpi_real, mpi_sum, mpi_comm_world, ierr)
            end if


            do i = 1, grid(iG)%nx
             do j = 1, yTop
                do k = 1, grid(iG)%nz
                   if (grid(iG)%active(i,j,k)>0) then
                      do freq = 1, nbins
                         if (lgDebug) then
                            grid(iG)%JDif(grid(iG)%active(i,j,k),freq) = &
                                 & JDifTemp(grid(iG)%active(i,j,k),freq)
                         end if
                         grid(iG)%JSte(grid(iG)%active(i,j,k),freq) = &
                              & JSteTemp(grid(iG)%active(i,j,k),freq)
                         if (lg2D) grid(iG)%JSte(grid(iG)%active(i,j,k),freq) = &
                              & grid(iG)%JSte(grid(iG)%active(i,j,k),freq)/&
                              & TwoDscaleJ(grid(iG)%active(i,j,k))
                         
                      end do
                      if (lgDebug) then
                         do freq = 1, nLines
                            grid(iG)%linePackets(grid(iG)%active(i,j,k),freq) &
                                 & = linePacketsTemp(grid(iG)%active(i,j,k),freq)
                            if (lg2D) grid(iG)%JDif(grid(iG)%active(i,j,k),freq) = &
                                 & grid(iG)%JDif(grid(iG)%active(i,j,k),freq)/&
                                 & TwoDscaleJ(grid(iG)%active(i,j,k))

                         end do
                      end if
                   end if
                end do
             end do
            end do

            call mpi_barrier(mpi_comm_world, ierr)


              if (lgDebug) then
                 if ( associated(linePacketsTemp) )    deallocate(linePacketsTemp)
                 if ( associated(JDifTemp) )           deallocate(JDifTemp)
              end if
              if ( associated(JSteTemp) )           deallocate(JSteTemp)           


              do i = 0, grid(iG)%nCells
                 grid(iG)%Jste(i,:) = grid(iG)%Jste(i,:) * 1.e-9
                 
                 if (lgDebug) grid(iG)%Jdif(i,:) = grid(iG)%Jdif(i,:) * 1.e-9
                 do ifreq = 1, nbins
                    totalEscaped = totalEscaped+&
                         & grid(iG)%escapedPackets(i,ifreq, 0)
                    do ian = 0, nAngleBins
                       grid(iG)%escapedPackets(i,ifreq,ian) = grid(iG)%escapedPackets(i,ifreq,ian)                           
                    end do
                 end do
                 
                 if (lgSymmetricXYZ) then
                    grid(iG)%Jste(i,:) = grid(iG)%Jste(i,:)/8.
                    grid(iG)%escapedPackets(i,:,:) = grid(iG)%escapedPackets(i,:,:)/8.
                    totalEscaped = totalEscaped/8.
                    if (lgDebug) grid(iG)%Jdif(i,:) = grid(iG)%Jdif(i,:)/8.
                 end if
                 
              end do

           end do


           if (lgDust .and. (nIterateMC>1 .or. .not.lgEquivalentTau)) then
              print*, "! iterateMC: [Interactions] : total -- abs -- sca: "
              print*, "! iterateMC: [Interactions] ", absInt+scaInt, " -- ", &
                   &  absInt*100./(absInt+scaInt),"% -- ", &
                   &  scaInt*100./(scaInt+absInt),"%"
           end if
           
           print*, " total Escaped Packets :",  totalEscaped              

           if (taskid==0) call writeSED(grid)                          
           if (taskid==0 .and. contCube(1)>0. .and. contCube(2)>0. ) & 
                & call writeContCube(grid, contCube(1),contCube(2))    

!*******************************************************************

           do iG = 1, nGrids
              
              if (ig>1 .or. (.not.lg2D)) then
                 yTop = grid(iG)%ny
              else if (iG ==1 .and. lg2D) then
                 yTop = 1
              end if


              allocate(lgConvergedTemp(0:grid(iG)%nCells), stat &
                   &= err)
              if (err /= 0) then
                 print*, "! iterateMC: can't allocate array&
                      & memory:lgConvergedTemp"
                 stop
              end if

              allocate(lgBlackTemp(0:grid(iG)%nCells), stat &
                   &= err)
              if (err /= 0) then
                 print*, "! iterateMC: can't allocate array&
                      & memory:lgBlackTemp  "
                 stop
              end if

              if (lgGas) then

                 allocate(NeTemp(0:grid(iG)%nCells), stat = err)
                 if (err /= 0) then
                    print*, "! iterateMC: can't allocate array memory:&
                         & NeTemp ", iG
                    stop
                 end if
                 NeTemp           = 0.
                 allocate(TeTemp(0:grid(iG)%nCells), stat = err)
                 if (err /= 0) then
                    print*, "! iterateMC: can't allocate array memory:&
                         & TeTemp ", iG
                    stop
                 end if
                 TeTemp           = 0.
                 allocate(ionDenTemp(0:grid(iG)%nCells, &
                      & nElementsUsed, nStages), stat = err)
                 if (err /= 0) then
                    print*, "! iterateMC: can't allocate array memory: &
                         &ionDenTemp ", iG
                    stop
                 end if
                 ionDenTemp       = 0.
              
              end if
             
              if (lgDust) then
                 allocate(TdustTemp(0:nSpeciesMax,0:nSizes,0:grid(iG)%nCells), stat = err)
                 if (err /= 0) then
                    print*, "! iterateMC: can't allocate array memory:&
                         &TdustTemp ", iG
                    stop
                 end if
                 TdustTemp = 0.
              end if


              grid(iG)%lgConverged(0:grid(iG)%nCells) = 0
              grid(iG)%lgBlack(0:grid(iG)%nCells)     = 0
              grid(iG)%noHit       = 0.
              grid(iG)%noIonBal    = 0.
              grid(iG)%noTeBal     = 0.         
              lgConvergedTemp      = 0
              lgBlackTemp          = 0 

               if (lgTraceHeating.and.taskid==0) then
                 open(file=trim(outputDir)//"thermalBalance.out",  action="write",unit=57, status="unknown", iostat=ios)
                 if (ios /= 0) then
                    print*, "! iterationMC: can't open file for writing, output/thermalBalance.out"
                    stop
                 end if
              end if

              if(taskid==0) print*, 'iterateMC: updateCell in', iG
              iCell = 0
              do i = 1, grid(iG)%nx
                 do j = 1, yTop
                    do k = 1, grid(iG)%nz
                       iCell = iCell+1
                       if (mod(iCell-(taskid+1),numtasks)==0) &
                            & call updateCell(grid(iG),i,j,k)

                    end do
                 end do
              end do
              if(taskid==0) print*, 'iterateMC: updateCell out', iG

              if (lgTraceHeating.and.taskid==0) then
                 close (57)
              end if


              
             

              size =  (grid(iG)%nCells+1)

              call mpi_allreduce(grid(iG)%lgConverged, lgConvergedTemp, size, &
                   & mpi_integer, mpi_sum, mpi_comm_world, ierr)
              call mpi_allreduce(grid(iG)%lgBlack, lgBlackTemp, size, &
                   & mpi_integer, mpi_sum, mpi_comm_world, ierr)


              if (lgGas) then
                 call mpi_allreduce(NeTemp, grid(iG)%Ne, size, &
                      & mpi_real, mpi_sum, mpi_comm_world, ierr)
                 call mpi_allreduce(TeTemp, grid(iG)%Te, size, &
                      & mpi_real, mpi_sum, mpi_comm_world, ierr)
           
           
                 size = (grid(iG)%nCells+1)*nElementsUsed*nStages
                 
                 call mpi_allreduce(ionDenTemp, grid(iG)%ionDen, size, &
                      & mpi_real, mpi_sum, mpi_comm_world, ierr)

              end if

              if (lgDust) then
                 size =  (grid(iG)%nCells+1)*(nspeciesMax+1)*(nsizes+1)

                 call mpi_allreduce(TdustTemp,grid(iG)%Tdust,size, &
                      & mpi_real, mpi_sum, mpi_comm_world, ierr)

                 if (lgGas .and. convPercent>=resLinesTransfer .and. (.not.lgResLinesFirst) &
                      & .and. .not.(nIterateMC==1)) then

                    allocate(budgetTemp(0:nAbComponents,0:nResLines), stat=err)
                    if (err /= 0) then
                       print*, "! iterateMC: can't allocate array memory:&
                            &budgetTemp "
                       stop
                    end if
                    budgetTemp=0.

                    size = (nAbcomponents+1)*(nResLines+1)
                    
                    call mpi_allreduce(dustHeatingBudget,budgetTemp,size, &
                         & mpi_real, mpi_sum, mpi_comm_world, ierr)

                    dustHeatingBudget = budgetTemp


                    deallocate(budgetTemp)
                    
                 end if

              end if

              do i = 1, grid(iG)%nx
                 do j = 1, yTop
                    do k = 1, grid(iG)%nz
                       if (grid(iG)%active(i,j,k)>0) then
                          grid(iG)%lgConverged(grid(iG)%active(i,j,k)) = &
                               & lgConvergedTemp(grid(iG)%active(i,j,k))           
                          grid(iG)%lgBlack(grid(iG)%active(i,j,k)) = &
                               & lgBlackTemp(grid(iG)%active(i,j,k))
                       end if
                    end do
                 end do
              end do


              call mpi_barrier(mpi_comm_world, ierr)
           
              if ( associated(lgConvergedTemp) )  deallocate(lgConvergedTemp)
              if ( associated(lgBlackTemp) )  deallocate(lgBlackTemp)
              if (lgGas) then
                 if ( associated(NeTemp) )           deallocate(NeTemp)
                 if ( associated(TeTemp) )           deallocate(TeTemp)
                 if ( associated(ionDenTemp) )       deallocate(ionDenTemp)
              end if
              if (lgDust) then
                 if ( associated(TdustTemp)) deallocate(TdustTemp)
              end if
              
           end do
           
!******************************************************************************

           if (nGrids>1) then
!              print*, " ! iterateMC: integratePathTau stuff still not implemented for multiple grids.... skipping"
           else
              call writeTauNu(grid)
           end if


           

    end function step_monte_carlo
    
    function calculate_convergence(grid)
        implicit none
        
        include 'mpif.h'
        
        type(grid_type), intent(inout) :: grid(*)
        integer :: calculate_convergence
        integer :: iG,  icomp, icontrib ! counters
        integer :: i,j,k ! counters
        integer :: ios, ierr ! counters
        integer :: totCells
        integer :: size
        integer :: yTop
        real    :: totheatdust 
        real, pointer          :: noHitPercent(:)    ! percentage of no Hit cells
        real, pointer          :: noIonBalPercent(:) ! percentage of cell where ion bal not conv
        real, pointer          :: noTeBalPercent(:)  ! percentage of cell where Te bal not conv
    
        calculate_convergence = 0
        
        allocate(noHitPercent(nGrids))
        allocate(noIonBalPercent(nGrids))
        allocate(noTeBalPercent(nGrids))  
        
        noHitPercent    = 0.
        noIonBalPercent  = 0.
        noTeBalPercent  = 0.
        
         
        do iG = 1, nGrids
            size = 1
        
            call mpi_allreduce(grid(iG)%noHit, noHitPercent(iG), size, &
                   & mpi_real, mpi_sum,mpi_comm_world, ierr)           

            if (lgGas) then
                call mpi_allreduce(grid(iG)%noIonBal, noIonBalPercent(iG), size, &
                  & mpi_real, mpi_sum, mpi_comm_world, ierr)           

                call mpi_allreduce(grid(iG)%noTeBal, noTeBalPercent(iG), size, &
                  & mpi_real, mpi_sum, mpi_comm_world, ierr)                                 
            end if 
        end do
        
    ! reinitialize convPercent and totCells           
           totPercent  = 0.

           ! calculate the percentage of converged cells
           do iG = 1, nGrids

              totCells    = 0
              convPercent = 0.

              if (ig>1 .or. (.not.lg2D)) then
                 yTop = grid(iG)%ny
              else if (iG ==1 .and. lg2D) then
                 yTop = 1
              end if

              do i = 1, grid(iG)%nx
                 do j = 1, yTop
                    do k = 1, grid(iG)%nz
                       if (grid(iG)%active(i,j,k)>0)  then
                          if (.not.lgEcho) then 
                             convPercent = convPercent + grid(iG)%lgConverged(grid(iG)%active(i,j,k))
                             totCells    = totCells + 1
                          else ! if light echo then only count echo cells!
                             if (grid(iG)%echoVol(i,j,k).gt.0.0) then
                                convPercent = convPercent + grid(iG)%lgConverged(grid(iG)%active(i,j,k))
                                totCells    = totCells + 1
                             end if
                          endif
                       end if
                    end do
                 end do
              end do

              convPercent               = 100.*convPercent / totCells
              noHitPercent              = 100.*noHitPercent(iG) / totCells
              grid(iG)%noIonBal         = 100.*noIonBalPercent(iG) / totCells
              grid(iG)%noTeBal          = 100.*noTeBalPercent(iG) / totCells 
              totPercent = totPercent + convPercent*grid(iG)%nCells/100.

              if (taskid == 0) then
                 if (nIterateMC == 1) then
                    close(21)
                    open(unit=21, status='unknown', position='rewind', file=trim(outputDir)//'summary.out', iostat=ios)
                    if (ios /= 0) then
                       print*, "! iterationMC: can't open file for writing, summary.out -1"
                       stop
                    end if
                 else
                    close(21)
                    open(unit=21, status='unknown', position='append', file=trim(outputDir)//'summary.out', iostat=ios)
                    if (ios /= 0) then
                       print*, "! iterationMC: can't open file for writing, summary.out -2"
                       stop
                    end if
                 end if                
                 
                 print*, "! iterateMC:  Summary] Iteration ",nIterateMC,'; ', &
                      &int(convPercent),"% converged cells in grid ", iG
                 print*, "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                      &int(noHitPercent(iG)),"% no hit cells in grid ", iG
                 print*, "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                      &int(grid(iG)%noIonBal),"% Ion balance not reached in grid ", iG
                 print*, "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                      &int(grid(iG)%noTeBal),"% Te balance not reached in grid", iG
                 write(21,*) "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                      &int(convPercent),"% converged cells in grid ", iG
                 write(21,*) "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                      &int(noHitPercent(iG)),"% no hit cells in grid ", iG
                 write(21,*) "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                      &int(grid(iG)%noIonBal),"% Ion balance not reached in grid ", iG
                 write(21,*) "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                      &int(grid(iG)%noTeBal),"% Te balance not reached in grid", iG
                 if (lgGas .and. convPercent>=resLinesTransfer .and. .not.lgResLinesFirst .and. &
                      & (.not. nIterateMC==1) ) then
                    write(21,*) "! iterateMC: [Summary] Iteration ",nIterateMC,'; '
                    write(21,*) "Dust Budgets: "
                    do icomp = 0, nAbComponents
                       write(21,*) " Component ", icomp
                       totheatdust = 0.
                       do icontrib =0, nResLines
                          totheatdust = totheatdust+dustHeatingBudget(icomp,icontrib)
                       end do
                       do icontrib =0, nResLines
                          write(21,*) " Contribution", icontrib, dustHeatingBudget(icomp,icontrib)/&
                               & totheatdust
                       end do
                    end do
                 end if
              end if
              
           end do
           
           if (lgDust .and. convPercent>=resLinesTransfer .and. lgGas .and. &
                & (.not. nIterateMC==1) .and. (.not.lgResLinesFirst)) dustHeatingBudget = 0.

           totCells = 0
           do iG =1,nGrids
              totCells = totCells+grid(iG)%nCells
           end do

           totPercent      = 100.*totPercent / totCells           

           if (taskid==0) then
              print*, "! iterateMC: [Summary] Iteration ",nIterateMC,'; Total:  ', &
                   &totPercent,"% converged cells over all grids"
              print*, "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                   &nPhotons, " energy packets used"
              write(21, *) "! iterateMC: [Summary] Iteration ",nIterateMC,'; Total:  ', &
                   &totPercent,"% converged cells over all grids"
              write(21,*) "! iterateMC: [Summary] Iteration ",nIterateMC,'; ', &
                   &nPhotons, " energy packets used"
              close(21)
           
              ! write grid to files for warm start
              
              if ( totPercent >= convWriteGrid) then

                 ! write old grid to disk
                 call writeGrid(grid(1:nGrids))

              end if
           
           end if
                
        deallocate(noHitPercent)
        deallocate(noIonBalPercent)
        deallocate(noTeBalPercent)  
        
           


!           if (Ldiffuse>0. .and. nIterateMC > 1 .and. totPercent < 95. .and. lgAutoPackets & 
!                &  .and. totPercentOld > 0.) then
!
!              if ( (totPercent-totPercentOld)/totPercentOld <= convIncPercent ) then
!                 nPhotonsDiffuseLoc = nPhotonsDiffuseLoc*nPhotIncrease
!
!                 if (taskid==0) &
!                      & print*, "! iterateMC: [talk] number of diffuse energy packets &
!                      &per cell increased to ", nPhotonsDiffuseLoc
!              end if
!              
!           end if
    
    end function
    
END MODULE


