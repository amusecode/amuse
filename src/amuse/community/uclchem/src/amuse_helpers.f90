MODULE uclchemhelper
    USE physicscore
    USE chemistry
    USE io
    USE constants
    IMPLICIT NONE
CONTAINS

    SUBROUTINE cloud(dictionary, outSpeciesIn,abundance_out,successFlag)
        !Subroutine to call a cloud model, used to interface with python
        ! Loads cloud specific subroutines and send to solveAbundances
        !
        !Args:
        ! dictionary - python parameter dictionary
        ! outSpeciesIn - list of species to output as a space separated string
        !Returns:
        ! abundance_out - list of abundances of species in outSpeciesIn
        ! successFlag - integer flag indicating success or fail

        USE cloud_mod

        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        DOUBLE PRECISION :: abundance_out(500)
        INTEGER :: successFlag
        !f2py intent(in) dictionary,outSpeciesIn
        !f2py intent(out) abundance_out,successFlag
        successFlag=1

        CALL solveAbundances(dictionary, outSpeciesIn,successFlag,initializePhysics,updatePhysics,updateTargetTime,sublimation)

        IF ((ALLOCATED(outIndx)) .and. (successFlag .ge. 0)) THEN 
            abundance_out(1:SIZE(outIndx))=abund(outIndx,1)
        END IF 
    END SUBROUTINE cloud

    SUBROUTINE solveAbundances(dictionary,outSpeciesIn,successFlag,&
        &modelInitializePhysics,modelUpdatePhysics,updateTargetTime,&
        &sublimation)
        ! Core UCLCHEM routine. Solves the chemical equations for a given set of parameters through time 
        ! for a specified physical model.
        ! Change behaviour of physics by sending different subroutine arguments - hence the need for model subroutines above
        ! dictionary - the parameter dictionary string reprenting a python dictionary
        ! outSpeciesIn - the species to output
        ! successFlag - Integer to indicate whether code completed successfully
        ! modelInitializePhysics - subroutine to initialize physics from a physics module
        ! modelUpdatePhysics - subroutine to update physics from a physics module
        ! updateTargetTime - subroutine to update the target time from a physics module
        ! sublimation - subroutine allowing physics module to directly modify abundances once per time step.
        CHARACTER(LEN=*) :: dictionary, outSpeciesIn
        EXTERNAL modelInitializePhysics,updateTargetTime,modelUpdatePhysics,sublimation

        INTEGER, INTENT(OUT) :: successFlag
        successFlag=1

        ! Set variables to default values
        INCLUDE 'defaultparameters.f90'
        !Read input parameters from the dictionary
        CALL dictionaryParser(dictionary, outSpeciesIn,successFlag)
        IF (successFlag .lt. 0) THEN
            successFlag=PARAMETER_READ_ERROR
            WRITE(*,*) 'Error reading parameter dictionary'
            RETURN
        END IF
        
        dstep=1
        currentTime=0.0
        timeInYears=0.0

        CALL fileSetup

        !Initialize core physics first then model specific
        !This allows model to overrule changes made by core
        call coreInitializePhysics(successFlag)
        CALL modelInitializePhysics(successFlag)

        IF (successFlag .lt. 0) then
            successFlag=PHYSICS_INIT_ERROR
            WRITE(*,*) 'Error initializing physics'
            RETURN
        END IF

        CALL initializeChemistry(readAbunds)
        CALL readInputAbunds !this won't do anything if no abundLoadFile was in input
        !CALL simpleDebug("Initialized")

    
        dstep=1

        call output

        !loop until the end condition of the model is reached
        DO WHILE (((endAtFinalDensity) .and. (density(1) < finalDens)) .or. &
            &((.not. endAtFinalDensity) .and. (timeInYears < finalTime)))
          
            currentTimeold=currentTime

            !Each physics module has a subroutine to set the target time from the current time
            CALL updateTargetTime
            !loop over parcels, counting from centre out to edge of cloud
            DO dstep=1,points
                !reset time if this isn't first depth point
                currentTime=currentTimeold
                !update chemistry from currentTime to targetTime
                CALL updateChemistry(successFlag)
                IF (successFlag .lt. 0) THEN
                    write(*,*) 'Error updating chemistry'
                    RETURN
                END IF

                !get time in years for output, currentTime is now equal to targetTime
                timeInYears= currentTime/SECONDS_PER_YEAR

                !Update physics so it's correct for new currentTime and start of next time step
                Call coreUpdatePhysics
                CALL modelUpdatePhysics
                !Sublimation checks if Sublimation should happen this time step and does it
                CALL sublimation(abund)

                !write this depth step now time, chemistry and physics are consistent
                CALL output
            END DO
        END DO
        CALL finalOutput
        CALL closeFiles
    END SUBROUTINE solveAbundances

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
END MODULE uclchemwrap