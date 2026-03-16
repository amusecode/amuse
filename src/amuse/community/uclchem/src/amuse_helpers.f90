MODULE uclchemhelper
    USE physicscore
    USE chemistry
    USE io
    USE constants
    use network
    IMPLICIT NONE
    type particle_type
    !Particles are defined here
        integer :: id
        double precision :: density !unit = cm^-3
        double precision :: temperature !unit = K
        double precision :: ionrate !unit = galactic cr ionization rate 
        double precision :: uvrad !unit = Habing
        double precision :: abundances(nSpec+1)
    end type
    integer :: nmols
    type(particle_type), allocatable :: particles(:)
  
    double precision :: tcurrent  ! time unit = yr
  
    integer :: nparticle
    integer :: tot_id

    character(len=500) :: out_species

    integer, parameter :: NMAX=1000000

    logical :: particles_searcheable=.FALSE.
    integer, allocatable :: pid(:)

CONTAINS
    function chem_initialize()  result(ret)
        integer :: ret
        tcurrent=0.
        nparticle=0
        tot_id=0
        nmols=nspec
        if(.not.allocated(particles)) allocate(particles(NMAX))
        particles(:)%density=0.
        
        ret=0
    end function

    function chem_commit_parameters() result(ret)
        integer :: ret
        ret=0
    end function

    function chem_commit_particles() result (ret)
        integer :: ret
        integer :: n
        integer :: i
     
        particles_searcheable=.FALSE.
        nparticle=clean_particles(particles)
        ret=0
      end function

    function chem_end() result(ret)
        integer :: ret 
        if(allocated(particles)) deallocate(particles)
        ret=0
     end function

    function set_particle_state(id,density,temperature,ionrate,uvrad) result(ret)
        integer :: ret
        integer :: id,index
        double precision :: density, temperature, ionrate, uvrad
        index=find_particle(id)
        if(index.LT.0) then
        ret=index
        return
        endif

        particles(index)%density=density
        particles(index)%temperature=temperature
        particles(index)%ionrate=ionrate
        particles(index)%uvrad=uvrad
        ret=0

    end function

    function get_particle_state(id,density,temperature,ionrate,uvrad) result(ret)
        integer :: ret
        integer :: id,index
        double precision :: density, temperature, ionrate, uvrad
        index=find_particle(id)
        if(index.LT.0) then
        ret=index
        return
        endif

        density=particles(index)%density
        temperature=particles(index)%temperature
        ionrate=particles(index)%ionrate
        uvrad=particles(index)%uvrad
        ret=0

    end function

    function get_particle_abundance(id, aid, abundance) result(ret)
        integer :: ret
        integer :: id,index,aid
        double precision :: abundance
        index=find_particle(id)
        if(index.LT.0) then
          ret=index
          return
        endif
        if(aid.LT.1.OR.aid.GT.500) then
          ret=-1
          return
        endif
        abundance=particles(index)%abundances(aid)
        ret=0
    end function

    function set_particle_abundance(id, aid, abundance) result(ret)
        integer :: ret
        integer :: id,index,aid
        double precision :: abundance
        index=find_particle(id)
        if(index.LT.0) then
          ret=index
          return
        endif
        if(aid.LT.1.OR.aid.GT.500) then
          ret=-1
          return
        endif
        particles(index)%abundances(aid)=abundance
        ret=0
      end function

    function add_particle(id,density,temperature,ionrate,uvrad) result(ret)
        use network
        integer :: ret
        integer :: i,id
        double precision :: density, temperature, ionrate, uvrad
        double precision :: x(500)
        particles_searcheable=.FALSE.
        id=new_id()  
        i=nparticle+1
      
        if(i.GT.NMAX) then
          ret=-1
          return
        endif
        particles(i)%id=id
        particles(i)%density=density
        particles(i)%temperature=temperature
        particles(i)%ionrate=ionrate
        particles(i)%uvrad=uvrad
        particles(i)%abundances=1.d-40
        
        if(density.GT.0) then
          !values taken from the UCLchem default parameters
          particles(i)%abundances(nh)  = 0.5 !H
          particles(i)%abundances(nh2) = 0.25   !H2
          particles(i)%abundances(nhe) = 0.1 !He
          particles(i)%abundances(ncx) = 1.77d-04 !C+ (C is fully ionised by default)
          particles(i)%abundances(nc) = 1.d-10 !C
          particles(i)%abundances(no) = 3.34d-04 !O
          particles(i)%abundances(nn) = 6.18d-05 !N
          particles(i)%abundances(nmg) = 2.256d-06 !Mg
          particles(i)%abundances(np) = 7.78d-08 !P
          particles(i)%abundances(nf) = 3.6d-08 !
          particles(i)%abundances(nsx) = 3.51d-6 !S
          particles(i)%abundances(nsix) = 1.78d-06 !Si
          particles(i)%abundances(nclx) = 3.39d-08 !Cl
          particles(i)%abundances(nelec) = 1.77d-04 + 3.51d-6 + 1.78d-06 + 3.39d-08 !electrons; sum of ions
        endif
      
        nparticle=nparticle+1
        ret=0
    end function

    function remove_particle(id) result(ret)
        integer :: ret
        integer :: i,id
        i=find_particle(id)
        if(i.LE.0) then
          ret=i
          return
        endif
        if(particles(i)%density.LT.0) then
          ret=-4
          return
        endif
        particles(i)%density=-1.
        ret=0
    end function

    function new_id()
        integer new_id
        tot_id=tot_id+1
        new_id=tot_id
    end function

    function find_particle(id_) result(index)
        !Function to find a particle with a given id. Taken from the Krome interface.
        use hashMod
        type(hash_type),save ::  hash
        integer id_,index
        integer, save :: nbod=0
        
        if(.NOT.particles_searcheable) then
          nbod=nparticle
          if(allocated(pid)) deallocate(pid)
          allocate(pid(nbod))
          pid(1:nbod)=particles(1:nbod)%id
          call initHash(nbod/2+1,nbod, pid,hash)
          particles_searcheable=.TRUE.
        endif
        
      
        index=find(id_,pid,hash)
      
        if(index.LE.0) then   
          index=-1
          return
        endif
        if(index.GT.nbod) then
          index=-2
          return
        endif
        if(pid(index).NE.id_) then   
          index=-3
          return
        endif
        
    end function

    function get_current_time(time) result(ret)
        integer :: ret
        double precision :: time 
        time = tcurrent 
        ret = 0
    end function
      
    function clean_particles(par) result(np)
        integer :: left,right,np
        type(particle_type), allocatable :: par(:)
        type(particle_type) :: tmp
        left=1
        if(.NOT.allocated(par)) then
          np = 0 
          return  
        endif
        right=size(par)
        if(right.EQ.0) then
          np=0
          return
        endif 
        do while(.TRUE.)
          do while(par(left)%density.GT.0.AND.left.LT.right)
            left=left+1
          enddo
          do while(par(right)%density.LE.0.AND.left.LT.right)
            right=right-1  
          enddo
          if(left.LT.right) then
            tmp=par(left)
            par(left)=par(right)
            par(right)=tmp
          else
            exit
          endif
        enddo
        if(par(left)%density.GT.0) left=left+1
        np=left-1
    end function  

    function simple_evolution(dictionary, outSpeciesIn) result(ret)
        !Evolves each particle seperately
        integer :: ret
        integer :: i, iret
        CHARACTER(LEN=*) :: dictionary(:)
        CHARACTER(LEN=*) :: outSpeciesIn

        ret = 0
        INCLUDE 'defaultparameters.f90'
        print *, outSpeciesIn
        do i=1,nparticle
            CALL dictionaryParser(dictionary(i), outSpeciesIn,ret)
            iret = evolve_1_particle(particles(i))
            ret = min(iret,ret)
        enddo
        tcurrent = timeInYears
        return

    end function

    function evolve_1_particle(part) result(ret)
        use cloud_mod
        type(particle_type) :: part
        integer :: ret        

        dstep=1


        call coreInitializePhysics(ret)
        CALL coreInitializePhysics(ret)

        CALL initializeChemistry(readabunds=.FALSE.)

        dstep=1
        !Make sure the saved abundances from previous steps are used
        abund(:,1) = part%abundances
        DO WHILE (((endAtFinalDensity) .and. (density(1) < finalDens)) .or. &
            &((.not. endAtFinalDensity) .and. (timeInYears < finalTime)))
            currentTimeold=currentTime
            !In standard UCLchem, each physics module has its own timestep scheme.
            CALL updateTargetTime
            !loop over parcels, counting from centre out to edge of cloud
            DO dstep=1,points
                !reset time if this isn't first depth point

                currentTime=currentTimeold
                !update chemistry from currentTime to targetTime
                CALL updateChemistry(ret)
                IF (ret .lt. 0) THEN
                    write(*,*) 'Error updating chemistry'
                    RETURN
                END IF

                !get time in years for output, currentTime is now equal to targetTime
                timeInYears= currentTime/SECONDS_PER_YEAR

                !Update physics so it's correct for new currentTime and start of next time step
                Call coreUpdatePhysics
                !Sublimation checks if Sublimation should happen this time step and does it
                CALL sublimation(abund)

                !write this depth step now time, chemistry and physics are consistent
                CALL output
            END DO
        END DO
        !Save the computed abundances to the particle set
        part%abundances=abund(:,1)
        currentTime = 0.0
    end function

    SUBROUTINE get_rates(dictionary,abundancesIn,speciesIndx,rateIndxs,&
        &speciesRates,successFlag,transfer,swap,bulk_layers)
        !Given a species of interest, some parameters and abundances, this subroutine
        !return the rate of all reactions that include that species plus some extra variables
        !to allow for the calculation of the rate of bulk/surface ice transfer.
        USE cloud_mod
        USE network, only : nspec
        CHARACTER(LEN=*):: dictionary
        DOUBLE PRECISION :: abundancesIn(500),speciesRates(500)
        DOUBLE PRECISION :: transfer,swap,bulk_layers
        INTEGER:: rateIndxs(500),speciesIndx, successFlag
        DOUBLE PRECISION :: ydot(nspec+1)
        INTEGER :: speci,bulk_version,surface_version
        !f2py intent(in) dictionary,abundancesIn,speciesIndx,rateIndxs
        !f2py intent(out) speciesRates,successFlag,transfer,swap,bulk_layers
        INCLUDE 'defaultparameters.f90'

        CALL dictionaryParser(dictionary, "",successFlag)
        IF (successFlag .lt. 0) THEN
            WRITE(*,*) 'Error reading parameter dictionary'
            RETURN
        END IF
        CALL coreInitializePhysics(successFlag)
        CALL initializePhysics(successFlag)
        IF (successFlag .lt. 0) then
            WRITE(*,*) 'Error initializing physics'
            RETURN
        END IF

        CALL initializeChemistry(readAbunds)
        dstep=1
        successFlag=1
        abund(:nspec,dstep)=abundancesIn(:nspec)
        abund(neq,dstep)=initialDens
        currentTime=0.0
        timeInYears=0.0

        targetTime=1.0d-7
        CALL updateChemistry(successFlag)

        CALL F(NEQ,currentTime,abund(:,dstep),ydot)

        speciesRates=rate(rateIndxs)

        IF ((specname(speciesIndx)(1:1) .eq. "#") .or.&
        & (specname(speciesIndx)(1:1) .eq. "@")) THEN
            DO speci=1,nSpec
                IF (specname(speci) .eq. "@"//specname(speciesIndx)(2:)) bulk_version=speci
                IF (specname(speci) .eq. "#"//specname(speciesIndx)(2:)) surface_version=speci
            END DO
            IF (YDOT(nsurface) .lt. 0) THEN
                transfer=YDOT(nsurface)*surfaceCoverage*abund(bulk_version,1)/safeBulk
            ELSE
                transfer=YDOT(nsurface)*surfaceCoverage*abund(surface_version,1)
            END If
            swap=totalSwap
            bulk_layers=bulkLayersReciprocal
        ELSE
            swap=0.0
            transfer=0.0
            bulk_layers=0.0
        END IF

    END SUBROUTINE get_rates

    SUBROUTINE dictionaryParser(dictionary, outSpeciesIn,successFlag)
        !Reads the input parameters from a string containing a python dictionary/JSON format
        !set of parameter names and values.
        !dictionary - lowercase keys matching the names of the parameters in the select case below
        !OutSpeciesIn - string containing the species to output
        !successFlag - integer flag to indicate success
        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        INTEGER, INTENT(OUT) :: successFlag
        INTEGER, ALLOCATABLE, DIMENSION(:) :: locations
        LOGICAL :: ChemicalDuplicateCheck
        INTEGER :: posStart, posEnd, whileInteger,inputindx
        CHARACTER(LEN=100) :: inputParameter, inputValue
        close(10)
        close(11)
        close(7)

        !always deallocate these so that if user didn't specify them,
        ! they don't remain from previous run
        IF (ALLOCATED(outIndx)) DEALLOCATE(outIndx)
        IF (ALLOCATED(outSpecies)) DEALLOCATE(outSpecies)

        !All reads use IOSTAT which will change successFlag from 0 if an error occurs
        !so set zero and check for non-zero each loop.
        successFlag=0

        whileInteger = 0

        posStart = scan(dictionary, '{')
        DO WHILE (whileInteger .NE. 1)
            posEnd = scan(dictionary, ':')
            inputParameter = dictionary(posStart+2:posEnd-2)
            dictionary = dictionary(posEnd:)
            posStart = scan(dictionary, ' ')
            IF (scan(dictionary, ',') .EQ. 0) THEN
                posEnd = scan(dictionary, '}')
                whileInteger = 1
            ELSE
                posEnd = scan(dictionary, ',')
            END IF
            inputValue = dictionary(posStart+1:posEnd-1)

            SELECT CASE (inputParameter)
                CASE('alpha')
                    !To provide alphas, set keyword alpha in inputdictionary with a dictionary value
                    !that dictionary should be index:value pairs for the alpha array    
                    posStart=scan(dictionary,'{')
                    posEnd=scan(dictionary,'}')
                    CALL coefficientParser(dictionary(posStart+1:posEnd),alpha)
                CASE('beta')
                    !To provide alphas, set keyword alpha in inputdictionary with a dictionary value
                    !that dictionary should be index:value pairs for the alpha array    
                    posStart=scan(dictionary,'{')
                    posEnd=scan(dictionary,'}')
                    CALL coefficientParser(dictionary(posStart+1:posEnd),beta)
                CASE('gamma')
                    !To provide alphas, set keyword alpha in inputdictionary with a dictionary value
                    !that dictionary should be index:value pairs for the alpha array    
                    posStart=scan(dictionary,'{')
                    posEnd=scan(dictionary,'}')
                    CALL coefficientParser(dictionary(posStart+1:posEnd),gama)
                CASE('initialtemp')
                    READ(inputValue,*,iostat=successFlag) initialTemp
                CASE('initialdens')
                    READ(inputValue,*,iostat=successFlag) initialDens
                CASE('finaldens')
                    READ(inputValue,*,iostat=successFlag) finalDens
                CASE('currenttime')
                    READ(inputValue,*,iostat=successFlag) currentTime
                CASE('finaltime')
                    READ(inputValue,*,iostat=successFlag) finalTime
                CASE('radfield')
                    READ(inputValue,*,iostat=successFlag) radfield
                CASE('zeta')
                    READ(inputValue,*,iostat=successFlag) zeta
                CASE('freezefactor')
                    READ(inputValue,*,iostat=successFlag) freezeFactor
                CASE('rout')
                    READ(inputValue,*,iostat=successFlag) rout
                CASE('rin')
                    READ(inputValue,*,iostat=successFlag) rin
                CASE('baseav')
                    READ(inputValue,*,iostat=successFlag) baseAv
                CASE('points')
                    READ(inputValue,*,iostat=successFlag) points
                CASE('endatfinaldensity')
                    Read(inputValue,*,iostat=successFlag) endAtFinalDensity
                CASE('freefall')
                    READ(inputValue,*,iostat=successFlag) freefall
                CASE('freefallfactor')
                    READ(inputValue,*,iostat=successFlag) freefallFactor
                CASE('desorb')
                    READ(inputValue,*,iostat=successFlag) desorb
                CASE('h2desorb')
                    READ(inputValue,*,iostat=successFlag) h2desorb
                CASE('crdesorb')
                    READ(inputValue,*,iostat=successFlag) crdesorb
                CASE('uvdesorb')
                    READ(inputValue,*,iostat=successFlag) uvdesorb
                CASE('thermdesorb')
                    READ(inputValue,*,iostat=successFlag) uvdesorb
                CASE('instantsublimation')
                    READ(inputValue,*,iostat=successFlag) instantSublimation
                CASE('cosmicrayattenuation')
                    READ(inputValue,*,iostat=successFlag) cosmicRayAttenuation
                CASE('ionmodel')
                    READ(inputValue,*,iostat=successFlag) ionModel
                CASE('improvedh2crpdissociation')
                    READ(inputValue,*,iostat=successFlag) improvedH2CRPDissociation
                CASE('ion')
                    READ(inputValue,*,iostat=successFlag) ion
                CASE('fhe')
                    READ(inputValue,*,iostat=successFlag) fhe
                CASE('fc')
                    READ(inputValue,*,iostat=successFlag) fc
                CASE('fo')
                    READ(inputValue,*,iostat=successFlag) fo
                CASE('fn')
                    READ(inputValue,*,iostat=successFlag) fn
                CASE('fs')
                    READ(inputValue,*,iostat=successFlag) fs
                CASE('fmg')
                    READ(inputValue,*,iostat=successFlag) fmg
                CASE('fsi')
                    READ(inputValue,*,iostat=successFlag) fsi
                CASE('fcl')
                    READ(inputValue,*,iostat=successFlag) fcl
                CASE('fp')
                    READ(inputValue,*,iostat=successFlag) fp
                CASE('ff')
                    READ(inputValue,*,iostat=successFlag) ff
                CASE('outspecies')
                    READ(inputValue,*,iostat=successFlag) nout
                    ALLOCATE(outIndx(nout))
                    ALLOCATE(outSpecies(nout))
                    IF (outSpeciesIn .eq. "") THEN
                        write(*,*) "Outspecies parameter set but no outspecies string given"
                        write(*,*) "general(parameter_dict,outSpeciesIn) requires a delimited string of species names"
                        write(*,*) "if outSpecies or columnFlag is set in the parameter dictionary"
                        successFlag=-1
                        RETURN
                    ELSE
                        READ(outSpeciesIn,*, END=22) outSpecies
                        IF (outSpeciesIn .eq. "") THEN
22                              write(*,*) "mismatch between outSpeciesIn and number given in dictionary"
                            write(*,*) "Number:",nout
                            write(*,*) "Species list:",outSpeciesIn
                            successFlag=-1
                            RETURN
                        END IF
                    END IF
                    !assign array indices for important species to the integers used to store them.
                    DO i=1,nspec
                        DO j=1,nout
                            IF (specname(i).eq.outSpecies(j)) outIndx(j)=i
                        END DO
                    END DO
                CASE('writestep')
                    READ(inputValue,*,iostat=successFlag) writeStep
                CASE('ebmaxh2')
                    READ(inputValue,*,iostat=successFlag) ebmaxh2
                CASE('epsilon')
                    READ(inputValue,*,iostat=successFlag) epsilon
                CASE('uvcreff')
                    READ(inputValue,*,iostat=successFlag) uvcreff
                CASE('ebmaxcr')
                    READ(inputValue,*,iostat=successFlag) ebmaxcr
                CASE('phi')
                    READ(inputValue,*,iostat=successFlag) phi
                CASE('ebmaxuvcr')
                    READ(inputValue,*,iostat=successFlag) ebmaxuvcr
                CASE('uv_yield')
                    READ(inputValue,*,iostat=successFlag) uv_yield
                CASE('metallicity')
                    READ(inputValue,*,iostat=successFlag) metallicity
                CASE('omega')
                    READ(inputValue,*,iostat=successFlag) omega
                CASE('reltol')
                    READ(inputValue,*,iostat=successFlag) reltol
                CASE('abstol_factor')
                    READ(inputValue,*,iostat=successFlag) abstol_factor
                CASE('abstol_min')
                    READ(inputValue,*,iostat=successFlag) abstol_min
                ! CASE('jacobian')
                !     READ(inputValue,*) jacobian
                CASE('abundsavefile')
                    READ(inputValue,*,iostat=successFlag) abundSaveFile
                    abundSaveFile = TRIM(abundSaveFile)
                    open(abundSaveID,file=abundSaveFile,status="unknown")
                CASE('abundloadfile')
                    READ(inputValue,*,iostat=successFlag) abundLoadFile
                    abundLoadFile = TRIM(abundLoadFile)
                    open(abundLoadID,file=abundLoadFile,status='old')
                CASE('outputfile')
                    READ(inputValue,*,iostat=successFlag) outFile
                    outputFile = trim(outFile)
                    fullOutput=.True.
                    open(outputId,file=outputFile,status='unknown',iostat=successFlag)
                    IF (successFlag .ne. 0) THEN
                        write(*,*) "An error occured when opening the output file!"//&
                                        & NEW_LINE('A')//&
                                    &" The failed file was ",outputFile&
                                    &, NEW_LINE('A')//"A common error is that the directory doesn't exist"&
                                    &//NEW_LINE('A')//"************************"
                        successFlag=-1
                        RETURN
                    END IF
                CASE('columnfile')
                    IF (trim(outSpeciesIn) .NE. '' ) THEN
                        columnOutput=.True.
                        READ(inputValue,*,iostat=successFlag) columnFile
                        columnFile = trim(columnFile)
                        open(columnId,file=columnFile,status='unknown')
                    ELSE
                        WRITE(*,*) "Error in output species. No species were given but a column file was given."
                        WRITE(*,*) "columnated output requires output species to be chosen."
                        successFlag=-1
                        RETURN
                    END IF

                CASE DEFAULT
                    WRITE(*,*) "Problem with given parameter: '", trim(inputParameter), "'."
                    WRITE(*,*) "This is either not supported yet, or invalid."
                    successFlag=-1
                    RETURN
            END SELECT
            dictionary = dictionary(posEnd:)
            IF (SCAN(dictionary,',') .eq. 0) whileInteger=1

            !check for failure
            IF (successFlag .ne. 0) THEN
                WRITE(*,*) "Error reading ",inputParameter
                write(*,*) "This is usually due to wrong type."
                successFlag=PARAMETER_READ_ERROR
                RETURN
            END IF 
        END DO

    END SUBROUTINE dictionaryParser

    SUBROUTINE coefficientParser(coeffDictString,coeffArray)
        !Similar to dictionaryParser, it reads a python dictionary
        !however, it's intended to read pairs of reaction indices and coefficient values
        !for the alpha, beta, and gama arrays.
        ! No return value, just modifies the coeffArray
        CHARACTER(LEN=*) :: coeffDictString
        REAL(dp), INTENT(INOUT) :: coeffArray(*)
        INTEGER :: inputIndx,posStart,posEnd
        CHARACTER(LEN=100) :: inputValue
        LOGICAL :: continue_flag
        
        continue_flag=.True.
        DO WHILE (continue_flag)
            !substring containing integer key
            posStart=1
            posEnd=SCAN(coeffDictString,':')
            !read it into index integer
            READ(coeffDictString(posStart:posEnd-1),*) inputindx

            !substring including alpha value for the index.
            posStart=posEnd+1
            posEnd=SCAN(coeffDictString,',')
            !last value will have a } instead of , so grab index and tell loop to finish
            IF (posEnd .eq. 0) THEN
                posEnd=SCAN(coeffDictString,"}")
                continue_flag=.False.
            END IF

            !read that substring
            inputValue=coeffDictString(posStart:posEnd-1)
            READ(inputValue,*) coeffArray(inputIndx)
            !update string to remove this entry
            coeffDictString=coeffDictString(posEnd+1:)
        END DO
    END SUBROUTINE coefficientParser

END MODULE uclchemhelper