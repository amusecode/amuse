
!J-shock paramterization
!Based on James et al. 2019 A&A 634
!https://ui.adsabs.harvard.edu/abs/2020A%26A...634A..17J/abstract
MODULE jshock_mod
    USE physicscore
    USE network
    USE constants
    USE sputtering
    IMPLICIT NONE
 
    REAL(dp) :: tstart,maxTemp,vMin,mfp,tCool,tShock,d,dMax,maxDens
    REAL(dp) :: t_lambda, n_lambda

    REAL(dp) :: z2,vs,v0,at
    REAL(dp), allocatable :: tn(:),ti(:),tgc(:),tgr(:),tg(:)

    !*******************************************************************

CONTAINS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Checks inputs make sense and then calculates a few constants and!
    ! sets up variables for the shock paramterization that follows    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE initializePhysics(successFlag)
        INTEGER, INTENT(OUT) :: successFlag
        successFlag=1
        !Reset variables for python wrap.
        
        cloudSize=(rout-rin)*pc

        if (freefall) THEN
            write(*,*) "Cannot have freefall on during jshock"
            Write(*,*) "setting freefall=0 and continuing"
            freefall=.False.
        ENDIF
        IF (points .gt. 1) THEN
            WRITE(*,*) "Cannot have more than one point in shock"
            successFlag=-1
            RETURN
        END IF

        density=initialDens

        ! Determine the maximum temperature
        maxTemp = (5e3)*(vs/10)**2
        currentTimeOld=0.0


        ! Determine minimum velocity
        vMin = ((-2.058e-07*(vs**4) + 3.844e-05*(vs**3) - 0.002478*(vs**2) + 0.06183*(vs) - 0.4254)**2)**0.5

        ! Determine the shock width (of the order of the mean free path)
        mfp = ((SQRT(2.0)*(1e3)*(pi*(2.4e-8)**2))**(-1))/1d4
        tShock = mfp/(vs*1d5)
        ! Determine shock width
        tCool = (1/initialDens)*1d6*(60*60*24*365)
        ! Determine the maximum density attained
        maxDens = vs*initialDens*(1d2)
        ! Determine the rate constants
        t_lambda = LOG(maxTemp/initialTemp)
        n_lambda = LOG(maxDens/initialDens)


        if (allocated(tn)) deallocate(tn,ti,tgc,tgr,tg)
        allocate(tn(points),ti(points),tgc(points),tgr(points),tg(points))

        currentTimeOld=0.0
        CALL sputteringSetup       
    END SUBROUTINE initializePhysics

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Called every time loop in main.f90. Sets the timestep for the next output from   !
    !UCLCHEM. This is also given to the integrator as the targetTime in chemistry.f90 !
    !but the integrator itself chooses an integration timestep.                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE updateTargetTime
        IF (timeInYears .gt. 1e6) THEN
            targetTime=(timeInYears+1e5)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 1.0d4) THEN
            targetTime=(timeInYears+1000)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 1.0d3) THEN
            targetTime=(timeInYears+100.)*SECONDS_PER_YEAR
        ELSE IF (timeInYears*SECONDS_PER_YEAR .lt. tShock) THEN
            targetTime=currentTime+0.05*tShock
        ELSE
            targetTime=1.1*currentTime
        END IF
    END SUBROUTINE updateTargetTime

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate shock properties for current time and set density, temperature and Av  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE updatePhysics

        ! Determine the shock velocity at the current time
        v0 = vs*(DEXP(LOG(vMin/vs)*(currentTime/(finalTime*60*60*24*365))))
        IF (v0 .lt. vMin) THEN
            v0 = vMin
        END IF

        ! Determine whether shock is still increasing the temperature
        ! Or whether it is in the post-shock cooling phase
        ! Or whether the temperature is now constant
        IF (currentTime .le. tShock) THEN
            tn(dstep) = ((currentTime/tShock)**2)*(maxTemp) + initialTemp
            density = (((currentTime/tShock)**3)*(4*initialDens))
            WHERE (density .lt. initialDens) density = initialDens
        ELSE IF (currentTime .gt. tShock .AND. currentTime .le. tCool) THEN
            ! Otherwise we're in the cooling phase
            tn(dstep) = maxTemp*DEXP(-t_lambda*(currentTime/(tCool)))
            density = (4*initialDens)*DEXP(n_lambda*(currentTime/(tCool)))

            ! Ensure the gas does not cool below around 10 K
            IF (tn(dstep) .le. 10) THEN
                tn(dstep) = 10
            END IF

            where(density .gt. maxDens) density = maxDens
        ELSE
            tn(dstep) = 10
            density = maxDens
        END IF
        gasTemp(dstep)=tn(dstep)
        dustTemp(dstep)=gasTemp(dstep)
    END SUBROUTINE updatePhysics


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This subroutine must be in every physics module.                                !
    ! It receives the abundance array and performs any sublimation related activity   !
    ! In hot core that means following thermalEvaporation subroutine.                 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE sublimation(abund)
        REAL(dp),INTENT(INOUT) :: abund(nspec+1,points)
        REAL(dp) :: timeDelta
        timeDelta=(currentTime-currentTimeOld)

        IF ((sum(abund(iceList,dstep)) .gt. 1d-25) .AND. (v0 .gt. 0))&
        & CALL sputterIces(abund(:,dstep),v0,gasTemp(dstep),density(dstep),timeDelta)
        WHERE(abund.lt. 1.0d-50) abund=0.0d-50        
    END SUBROUTINE sublimation

END MODULE jshock_mod
