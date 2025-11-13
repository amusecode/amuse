! Models the chemistry in a collapsing prestellar core
! F. D. Priestley et al 2018 AJ 156 51 (https://ui.adsabs.harvard.edu/abs/2018AJ....156...51P/abstract)
! Uses the following parameterizations of MHD models.
! collapse = 2: Bonnor-Ebert sphere, overdensity factor 1.1 (Aikawa+2005)
! collapse = 3: Bonnor-Ebert sphere, overdensity factor 4 (Aikawa+2005)
! collapse = 4: magnetised filament, initially unstable to collapse (Nakamura+1995)
! collapse = 5: magnetised cloud, initially stable, collapse due to ambipolar diffusion (Fiedler+1993)
MODULE collapse_mod
   USE physicscore
   USE network
   USE constants
   IMPLICIT NONE
   
   INTEGER :: collapse_mode
   REAL(dp) :: maxTime
   REAL(dp), allocatable :: massInRadius(:),parcelRadius(:)
   REAL(dp) :: dt,drad
   CHARACTER (LEN=100) :: collapseFile
   LOGICAL :: writePhysics
CONTAINS
   
    SUBROUTINE initializePhysics(successFlag)
        INTEGER, INTENT(OUT) :: successFlag
       
        IF (ALLOCATED(parcelRadius)) DEALLOCATE(parcelRadius,massInRadius)
        ALLOCATE(parcelRadius(points),massInRadius(points))
      
         SELECT CASE(collapse_mode)
            CASE(1)
                maxTime=1.175d6
                finalTime=0.97*maxTime
            CASE(2) 
                maxTime=1.855d5
                finalTime=0.97*maxTime
            CASE(3)
            CASE(4)
            CASE DEFAULT
                write(*,*) "unacceptable collapse mode"
                successFlag=-1
                RETURN
         END SELECT

         DO dstep=1,points
               parcelRadius(dstep)=dstep*rout/float(points)
         END DO
         
         IF (writePhysics) OPEN(unit=66,file=collapseFile,status='unknown',err=99)
         IF (successFlag .lt. 0) THEN
            99 write(*,*) "could not open physics output file",collapseFile
            successFlag=-1
            RETURN
         END IF
         density=rhofit(rin,rho0fit(timeInYears),r0fit(timeInYears),afit(timeInYears))
         IF (collapse_mode .le. 2) CALL findmassInRadius
    END SUBROUTINE initializePhysics

    SUBROUTINE updateTargetTime
        IF (timeInYears .gt. 10000) THEN
            targetTime=(timeInYears+1000.0)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 1000) THEN
            targetTime=(timeInYears+100.0)*SECONDS_PER_YEAR
        ELSE IF (timeInYears .gt. 0.0) THEN
            targetTime=(timeInYears*10)*SECONDS_PER_YEAR
        ELSE
            targetTime=3.16d7*10.d-8
        ENDIF

       !IF (targetTime .gt. finalTime*SECONDS_PER_YEAR) targetTime=finalTime*SECONDS_PER_YEAR
    END SUBROUTINE updateTargetTime

    !This routine is formed for every parcel at every time step.
    !update any physics here. For example, set density
    SUBROUTINE updatePhysics
        !calculate column density. Remember dstep counts from core to edge
        !and coldens should be amount of gas from edge to parcel.
        call findcoldens(coldens(dstep),rin,rho0fit(timeInYears),r0fit(timeInYears),afit(timeInYears),rout)
        !calculate the Av using an assumed extinction outside of core (baseAv), depth of point and density
        av(dstep)= baseAv +coldens(dstep)/1.6d21
        !If collapse is one of the parameterized modes, find new density and radius
        
        IF ((collapse_mode .le. 2)) THEN
            !I changed rin to rout
            CALL findNewRadius(massInRadius(dstep),rout,rho0fit(timeInYears),&
                &r0fit(timeInYears),afit(timeInYears),parcelRadius(dstep))
        ELSE
            dt = targetTime - currentTime
            drad = vrfit(parcelRadius(dstep),rminfit(timeInYears),vminfit(timeInYears),avfit(timeInYears))*dt/pc
            parcelRadius(dstep) = parcelRadius(dstep) + drad
            IF (writePhysics) THEN
               write(66,*) timeInYears,parcelRadius(dstep),rhofit(parcelRadius(dstep),&
                        &rho0fit(timeInYears),r0fit(timeInYears),afit(timeInYears)),&
                        &vrfit(parcelRadius(dstep),rminfit(timeInYears),&
                        &vminfit(timeInYears),avfit(timeInYears))
            END IF
        END IF
        density(dstep)=rhofit(parcelRadius(dstep),rho0fit(timeInYears),r0fit(timeInYears),afit(timeInYears))        
    END SUBROUTINE updatePhysics

    !This module is isothermal and as such, no sublimation occurs.
    !This is a dummy subroutine.
    SUBROUTINE sublimation(abund)
        REAL(dp) :: abund(nspec+1,points)
        INTENT(IN) :: abund

    END SUBROUTINE sublimation


    ! finds initial mass within starting radius, assuming spherical symmetry
    SUBROUTINE findMassInRadius
      REAL(dp) :: rho0,r0,a
      INTEGER :: i,np,dstep
      REAL(dp) :: dr,drho

        rho0=rho0fit(timeInYears)
        r0=r0fit(timeInYears)
        a=afit(timeInYears)
      DO dstep=1,points
        np = 1000
        dr = parcelRadius(dstep)/np
        massInRadius(dstep) = 0.0d0

        DO i=1,np
           drho = 0.5d0*(rhofit(i*dr,rho0,r0,a)+rhofit((i-1)*dr,rho0,r0,a))
           massInRadius(dstep) = massInRadius(dstep) + drho*dr*(i*dr)**2
        END DO
      END DO
    END SUBROUTINE findMassInRadius

! finds radius enclosing a mass of massInRadius
    SUBROUTINE findNewRadius(massInRadius,r,rho0,r0,a,newRadius)
      REAL(dp),intent(in) :: massInRadius,r,rho0,r0,a
      REAL(dp),intent(out) :: newRadius
      INTEGER :: i
      REAL(dp) :: dr,drho,m1

      i=1
      dr = r/1.0d4
      m1 = 0.0d0
      DO WHILE (m1 .lt. massInRadius)
         drho = 0.5d0*(rhofit(i*dr,rho0,r0,a)+rhofit((i-1)*dr,rho0,r0,a))
         m1 = m1 + drho*dr*(i*dr)**2
         newRadius = i*dr
         i=i+1
      END DO
      IF (writePhysics) write(66,*) timeInYears,newRadius,rhofit(newRadius,rho0,r0,a),m1

    END SUBROUTINE findNewRadius

! finds column density to edge of cloud based on density profile
    SUBROUTINE findcoldens(coldens,rin,rho0,r0,a,rout)
      REAL(dp),intent(in) :: rin,rout,rho0,r0,a
      REAL(dp),intent(out) :: coldens
      INTEGER :: i,np
      REAL(dp) :: dr,drho,size,r1,r2

      np = 10000
      size = rout-rin
      dr = size/np
      coldens = 0.0d0
      IF (size .le. 0.0d0) return

      DO i=1,np
         r1 = rin + (i-1)*dr
         r2 = rin + i*dr
         drho = 0.5d0*(rhofit(r2,rho0,r0,a)+rhofit(r1,rho0,r0,a))
         coldens = coldens + drho*dr*pc
      END DO

    END SUBROUTINE findcoldens

! fit to density profile of hydrodynamic simulations
    REAL(dp) FUNCTION rhofit(r,rho0,r0,a)
      REAL(dp) :: r,rho0,r0,a
      REAL(dp) :: rau,unitrho,unitr,r75

      IF (collapse_mode .eq. 1) then
         rau = r*au
         rhofit = rho0/(1 + (rau/r0)**a)
      ELSE IF (collapse_mode .eq. 2) then
         rau = r*au
         rhofit = rho0/(1 + (rau/r0)**a)
      ELSE IF (collapse_mode .eq. 3) then
         unitrho = 2.2d4
         unitr = sqrt(1.38d-16*10/2/mh)*(2*pi*6.67d-8*2.2d4*mh)**(-0.5) ! distance unit equal to c_s * (2 pi G rho0)**-1/2
         unitr = unitr/pc
         rhofit = unitrho*rho0/(1+(r/unitr/r0)**2)**a
      ELSE IF (collapse_mode .eq. 4) then
         r75 = r/7.5d-1
         rhofit = rho0/(1 + (r75/r0)**a)
      END IF

    END FUNCTION rhofit

! fit to time evolution of central density
    REAL(dp) FUNCTION rho0fit(t)
      REAL(dp) :: t,logrho0,unitt
      IF (collapse_mode .eq. 1) then
         logrho0 = 61.8*(maxTime-t)**(-0.01) - 49.4
         rho0fit = 10**logrho0
      ELSE IF (collapse_mode .eq. 2) then
         logrho0 = 68.4*(maxTime-t)**(-0.01) - 55.7
         rho0fit = 10**logrho0
      ELSE IF (collapse_mode .eq. 3) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5) ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         logrho0 = 3.54*(5.47-t/unitt)**(-0.15) - 2.73
         rho0fit = 10**logrho0
      ELSE IF (collapse_mode .eq. 4) then
         IF (t .le. 6.0d0) then
            rho0fit = 2.0d3 + 1.7d3*(t/6.0 - 1.0)
         ELSE
            logrho0 = 5.3*(16.138-1d-6*t)**(-0.1) - 1.0
            rho0fit = 10**logrho0
         END IF
      END IF

    END FUNCTION rho0fit

! fit to time evolution of radius parameter
    REAL(dp) FUNCTION r0fit(t)
      REAL(dp) :: t,logr0,unitt

      IF (collapse_mode .eq. 1) then
         logr0 = -28.5*(maxTime-t)**(-0.01) + 28.93
         r0fit = 10**logr0
      ELSE IF (collapse_mode .eq. 2) then
         logr0 = -39.0*(maxTime-t)**(-0.01) + 38.7
         r0fit = 10**logr0
      ELSE IF (collapse_mode .eq. 3) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5) ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         logr0 = -1.34*(5.47-t/unitt)**(-0.15) + 1.47
         r0fit = 10**logr0
      ELSE IF (collapse_mode .eq. 4) then
         logr0 = -2.57*(16.138-1d-6*t)**(-0.1) + 1.85
         r0fit = 10**logr0
      END IF

    END FUNCTION r0fit

! fit to time evolution of density slope parameter
    REAL(dp) FUNCTION afit(t)
      REAL(dp) :: t,unitt

      IF (collapse_mode .eq. 1) then
         afit = 2.4d0
      ELSE IF (collapse_mode .eq. 2) then
         afit = 1.9 + 0.5*exp(-t/1e5)
      ELSE IF (collapse_mode .eq. 3) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5) ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         afit = 2.0 - 0.5*(t/unitt/5.47)**9
      ELSE IF (collapse_mode .eq. 4) then
         afit = 2.4 - 0.2*(1d-6*t/16.138)**40
      END IF

    END FUNCTION afit

! fit to radial velocity of hydrodynamical simulation
    REAL(dp) FUNCTION vrfit(r,rmin,vmin,a)
      REAL(dp) :: r,rmin,vmin,a
      REAL(dp) :: unitr,newRadius,rmid,r75

      IF (collapse_mode .eq. 3) then
         unitr = sqrt(1.38d-16*10/2/mh)*(2*pi*6.67d-8*2.2d4*mh)**(-0.5) ! distance unit equal to c_s * (2 pi G rho0)**-1/2
         unitr = unitr/pc
         newRadius = r/unitr - rmin
         IF (newRadius .lt. 0.0d0) then
            vrfit = vmin*((newRadius/rmin)**2 -1)
         ELSE
            vrfit = vmin*(exp(-2.0d0*a*newRadius) - 2*exp(-a*newRadius))
         END IF
         vrfit = sqrt(1.38d-16*10/2/mh)*vrfit ! convert to cm s-1 using c_s
      ELSE IF (collapse_mode .eq. 4) then
         rmid = 0.5
         r75 = r/7.5d-1
         newRadius = r75 - rmin
         IF (r75 .lt. rmin) then
            vrfit = vmin*((newRadius/rmin)**2 - 1)
         ELSE IF (r75 .le. rmid) then
            vrfit = (vmin-a)*(newRadius/(rmid-rmin))**0.3 - vmin
         ELSE
            vrfit = a/(1.0-rmid)*(r75-rmid) - a
         END IF
         vrfit = 1d3*vrfit ! convert to cm s-1 from 1e-2 km s-1
      END IF

    END FUNCTION vrfit

! fit to time evolution of radius of minimum velocity
    REAL(dp) FUNCTION rminfit(t)
      REAL(dp) :: t,unitt,tnew
      REAL(dp) :: t6

      IF (collapse_mode .eq. 3) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5) ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         tnew = t/unitt
         IF (tnew .eq. 0.0d0) then
            rminfit = 7.2d0
         ELSE IF (log(tnew) .lt. 1.6d0) then
            rminfit = -1.149*tnew + 7.2
         ELSE IF (log(tnew) .lt. 1.674d0) then
            rminfit = -9.2*log(tnew) + 16.25
         ELSE
            rminfit = -22.0*log(tnew) + 37.65
         END IF
      ELSE IF (collapse_mode .eq. 4) then
         t6 = 1d-6*t
         IF (t6 .le. 10.2) then
            rminfit = -0.0039*t6 + 0.49
         ELSE IF (t6 .le. 15.1) then
            rminfit = -0.0306*(t6-10.2) + 0.45
         ELSE
            rminfit = -0.282*(t6-15.1) + 0.3
         END IF
      END IF

    END FUNCTION rminfit

! fit to time evolution of minimum velocity
    REAL(dp) FUNCTION vminfit(t)
      REAL(dp) :: t,unitt,tnew
      REAL(dp) :: t6

      IF (collapse_mode .eq. 3) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5) ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         tnew = t/unitt
         IF (tnew .eq. 0.0d0) then
            vminfit = 0.0d0
         ELSE IF (log(tnew) .lt. 1.6d0) then
            vminfit = 0.0891*tnew
         ELSE IF (log(tnew) .lt. 1.674d0) then
            vminfit = 5.5*log(tnew) - 8.37
         ELSE
            vminfit = 18.9*log(tnew) - 30.8
         END IF
      ELSE IF (collapse_mode .eq. 4) then
         t6 = 1d-6*t
         vminfit = 3.44*(16.138-t6)**(-0.35) - 0.7
      END IF

    END FUNCTION vminfit

! fit to time evolution of velocity a-parameter (collapse 4) or velocity at r=0.5 (collapse 5)
    REAL(dp) FUNCTION avfit(t)
      REAL(dp) :: t,unitt,tnew
      REAL(dp) :: t6

      IF (collapse_mode .eq. 3) then
         unitt = (2*pi*6.67d-8*2.2d4*mh)**(-0.5) ! time unit equal to (2 pi G rho0)**-1/2
         unitt = unitt/SECONDS_PER_YEAR
         tnew = t/unitt
         IF (tnew .eq. 0.0d0) then
            avfit = 0.4d0
         ELSE IF (log(tnew) .lt. 1.6d0) then
            avfit = 0.0101*tnew + 0.4
         ELSE IF (log(tnew) .lt. 1.674d0) then
            avfit = 0.695*log(tnew) - 0.663
         ELSE
            avfit = 2.69*log(tnew) - 4.0
         END IF
      ELSE IF (collapse_mode .eq. 4) then
         t6 = 1d-6*t
         IF (t6 .le. 10.2) then
            avfit = 0.143*t6
         ELSE
            avfit = 0.217*(t6-10.2) + 1.46
         END IF
      END IF

    END FUNCTION avfit
END MODULE collapse_mod

!REQUIRED PHYSICS ENDS HERE, ANY ADDITIONAL PHYSICS CAN BE ADDED BELOW.
