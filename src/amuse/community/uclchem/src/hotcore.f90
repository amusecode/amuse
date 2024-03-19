! hotcore simulates a hot core/corino by following the increase in temperature as a function of time
! It can also reproduce the episodic thermal sublimation seen in laboratory experiments for two phase models
! It is based on Viti et al. 2004 and Collings et al. 2004
MODULE hotcore
    USE physicscore
    USE network
    USE constants
    IMPLICIT NONE
    !Flags let physics module control when evap takes place.flag=0/1/2 corresponding to not yet/evaporate/done
    INTEGER :: solidflag,volcflag,coflag
    
    !Arrays for phase 2 temp profiles. parameters for equation chosen by index
    !arrays go [1Msun,5, 10, 15, 25,60]
    INTEGER, PARAMETER :: nMasses= 6 
    INTEGER :: tempIndx
    REAL(dp),PARAMETER :: tempa(nMasses)=(/1.927d-1,4.8560d-2,7.8470d-3,9.6966d-4,1.706d-4,4.74d-7/)
    REAL(dp),PARAMETER :: tempb(nMasses)=(/0.5339,0.6255,0.8395,1.085,1.289,1.98/)
    REAL(dp),PARAMETER :: solidtemp(nMasses)=(/20.0,19.6,19.45,19.3,19.5,20.35/)
    REAL(dp),PARAMETER :: volctemp(nMasses)=(/84.0,86.3,88.2,89.5,90.4,92.2/)
    REAL(dp),PARAMETER :: codestemp(nMasses)=(/95.0,97.5,99.4,100.8,101.6,103.4/)
    REAL(dp), allocatable :: monoFracCopy(:)
    REAL(dp) :: maxTemp
contains

    SUBROUTINE initializePhysics(successFlag)
        INTEGER, INTENT(OUT) :: successFlag
        successFlag=0

        ! Modules not restarted in python wraps so best to reset everything manually.
        IF (ALLOCATED(monoFracCopy)) DEALLOCATE(monoFracCopy)
        ALLOCATE(monoFracCopy(size(monoFractions)))
        coFlag=0 !reset sublimation
        solidFlag=0
        volcFlag=0
        monoFracCopy=monoFractions !reset monofractions

        IF (freefall) density=1.001*initialDens

        IF (tempindx .gt. nMasses) THEN
            write(*,*) "tempindx was ",tempindx
            write(*,*) "tempindx must be less than",nMasses
            write(*,*) "1=1Msol, 2=5M, 3=10M, 4=15M, 5=25M, 6=60M"
            successFlag=-1
            RETURN
        END IF 
    END SUBROUTINE

    !Called every time loop in main.f90. Sets the timestep for the next output from   
    !UCLCHEM. This is also given to the integrator as the targetTime in chemistry.f90 
    !but the integrator itself chooses an integration timestep.                       
    SUBROUTINE updateTargetTime
        IF (timeInYears .gt. 1.0d6) THEN !code in years for readability, targetTime in s
            targetTime=(timeInYears+1.0d5)*SECONDS_PER_YEAR
        ELSE  IF (timeInYears .gt. 1.0d5) THEN
            targetTime=(timeInYears+1.0d4)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 1.0d4) THEN
            targetTime=(timeInYears+1000.0)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 1000) THEN
            targetTime=(timeInYears+100.0)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 0.0) THEN
            targetTime=(timeInYears*10.0)*SECONDS_PER_YEAR
        ELSE
            targetTime=SECONDS_PER_YEAR*1.0d-7
        ENDIF
    END SUBROUTINE updateTargetTime

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !This is called every time/depth step from main.f90                               !
    !Update the density, temperature and av to their values at currentTime            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE updatePhysics
         IF (gasTemp(dstep) .lt. maxTemp) THEN
        !Below we include temperature profiles for hot cores, selected using tempindx
        !They are taken from Viti et al. 2004 with an additional distance dependence from Nomura and Millar 2004.
        !It takes the form T=A(t^B)*[(d/R)^-0.5], where A and B are given below for various stellar masses
            gasTemp(dstep)=(cloudSize/(rout*pc))*(real(dstep)/real(points))
            gasTemp(dstep)=gasTemp(dstep)**(-0.5)
            gasTemp(dstep)=initialTemp + ((tempa(tempindx)*(currentTime/SECONDS_PER_YEAR)**tempb(tempindx))*gasTemp(dstep))
            if (gasTemp(dstep) .gt. maxTemp) gasTemp(dstep)=maxTemp
        END IF
        dustTemp=gasTemp
    END SUBROUTINE updatePhysics

    SUBROUTINE sublimation(abund)
    ! This subroutine mimics episodic thermal desorption if the network is two pahse
        REAL(dp) :: abund(nspec+1,points)
        INTENT(INOUT) :: abund
        IF (.not. THREE_PHASE) THEN
            IF (instantSublimation) THEN
                instantSublimation=.False.
                CALL totalSublimation(abund)
            ELSE IF (coflag .ne. 2) THEN
                IF (gasTemp(dstep) .gt. solidtemp(tempindx) .and. solidflag .ne. 2) solidflag=1
                IF (gasTemp(dstep) .gt. volctemp(tempindx) .and. volcflag .ne. 2) volcflag=1
                IF (gasTemp(dstep) .gt. codestemp(tempindx)) coflag=1
                CALL thermalEvaporation(abund)
            END IF
        END IF
    END SUBROUTINE sublimation

    SUBROUTINE thermalEvaporation(abund)
        !Evaporation is based on Viti et al. 2004. A proportion of the frozen species is released into the gas phase
        !in specific events. These events are activated by flags (eg solidflag) which can be set in physics module.
        !The species evaporated are in lists, created by Makerates and based on groupings. see the viti 2004 paper.
        REAL(dp) :: abund(nspec+1,points)
        INTENT(INOUT) :: abund
       
            IF (sum(abund(iceList,dstep)) .gt. 1d-30) THEN
                !Solid Evap
                IF (solidflag .eq. 1) THEN
                    CALL partialSublimation(solidFractions,abund)
                    solidflag=2
                ENDIF
    
                !monotonic evaporation at binding energy of species
                CALL bindingEnergyEvap(abund)
    
                !Volcanic evap
                IF (volcflag .eq. 1) THEN
                    CALL partialSublimation(volcanicFractions,abund)
                    volcflag=2 !Set flag to 2 to stop it being recalled
                ENDIF
    
                !Co-desorption
                IF (coflag .eq. 1) THEN
                    CALL totalSublimation(abund)
                    coflag=2
                ENDIF
            ENDIF
    END SUBROUTINE thermalEvaporation

    SUBROUTINE partialSublimation(fractions, abund)
        REAL(dp) :: abund(nspec+1,points)
        REAL(dp) :: fractions(:)

        abund(gasiceList,dstep)=abund(gasiceList,dstep)+fractions*abund(iceList,dstep)
        abund(iceList,dstep)=(1.0-fractions)*abund(iceList,dstep)

    END SUBROUTINE partialSublimation

    SUBROUTINE totalSublimation(abund)
        REAL(dp) :: abund(nspec+1,points)
        abund(gasiceList,dstep)=abund(gasiceList,dstep)+abund(iceList,dstep)
        abund(iceList,dstep)=1d-30
    END SUBROUTINE totalSublimation

    SUBROUTINE bindingEnergyEvap(abund)
        REAL(dp) :: abund(nspec+1,points)
        REAL(dp), parameter :: SURFACE_SITE_DENSITY = 1.5d15
        INTENT(INOUT) :: abund
        INTEGER :: i
        !Subroutine to handle mono-evaporation. See viti 2004
        REAL(dp) en,newm,expdust,freq,kevap
        integer speci
        !mono evaporation at the binding energy of each species
        DO i=lbound(iceList,1),ubound(iceList,1)
            speci=iceList(i)
            en=bindingEnergy(i)*K_BOLTZ_SI
            expdust=bindingEnergy(i)/gasTemp(dstep)
            newm = mass(speci)*1.66053e-27
            freq = dsqrt((2*(SURFACE_SITE_DENSITY)*en)/((pi**2)*newm))
            kevap=freq*exp(-expdust)
            IF (kevap .ge. 0.99) THEN
                abund(gasiceList(i),dstep)=abund(gasiceList(i),dstep)+(monoFracCopy(i)*abund(speci,dstep))
                abund(speci,dstep)=(1.0-monoFracCopy(i))*abund(speci,dstep)
                monoFracCopy(i)=0.0
            END IF 
        END DO
    END SUBROUTINE bindingEnergyEvap
END MODULE hotcore 
