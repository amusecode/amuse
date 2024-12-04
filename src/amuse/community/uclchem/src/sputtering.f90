!Sputtering module for the shock parameterizations
!Takes a velocity from the shock module and uses the shock sputtering
!treatment described in  Jimenez-Serra et al. 2008 A&A 482
!http://adsabs.harvard.edu/abs/2008A&A...482..549J
! to calulate the amount of material removed from dust grains
!
! Additionally uses Guillet et al. 2011 (http://www.aanda.org/10.1051/0004-6361/201015973)
! result that for shocks >19km/s you get vaporization which sends dust grain material
! into the gas phase. We make use of this by setting 19km/s as limit above which 
! refratory dust grain material is sputtered.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE sputtering
      USE network
      USE constants
      USE SurfaceReactions, only: GAS_DUST_DENSITY_RATIO,GRAIN_RADIUS
      USE physicscore, only: timeInYears
      IMPLICIT NONE

      INTEGER :: projectiles(6)
      REAL(dp) :: sConst,eta,epso
      !Speed at which refractory species are also removed from dust grains during sputtering. 19.0 km/s taken from Guillet et al. 2011. (see above)
      REAL(dp), PARAMETER :: VAPORIZE_SPEED=19.0 
      INTEGER, ALLOCATABLE :: sputters(:),gasSputters(:)
CONTAINS
    SUBROUTINE sputteringSetup
        INTEGER :: i,j,k,new_size
        LOGICAL :: found

        !Need abundance of major projectiles to get sputtering rate
        projectiles=(/nh2,nhe,nc,no,nsi,nco/)


        !check for refractory species and create a sublist of ice species
        !that only includes species not in refractory list
        !This allows us to sputter only volatile species when shockvelocity is low
        IF (ALLOCATED(sputters)) DEALLOCATE(sputters,gasSputters)

        IF (refractoryList(1) .gt. 0) THEN
            new_size=size(iceList)-size(refractoryList)
            ALLOCATE(sputters(new_size),gasSputters(new_size))
            k=1
            DO i=1,size(iceList)
                found=.False.
                DO j=1,size(refractoryList)
                    IF (iceList(i) .eq. refractoryList(j)) found=.true.
                END DO
                IF (.not. found) THEN
                    sputters(k)=iceList(i)
                    gasSputters(k)=gasIceList(i)
                    k=k+1
                END IF
            
            END DO
        ELSE
            ALLOCATE(sputters(size(iceList)))
            ALLOCATE(gasSputters(size(iceList)))
            gasSputters=gasIceList
            sputters=iceList
        END IF
    END SUBROUTINE sputteringSetup

    SUBROUTINE sputterIces(abund,shockVel,gasTemp,density,timeDelta)
      ! Sputter ices following Jimenez-Serra 2008 treatment
      ! Args:
      !     abund: abundances of all species
      !     shockVel: relative velocity of dust and colliding gas (shockVel in c-shock)
      REAL(dp), INTENT(INOUT) :: abund(nspec+1)
      REAL(dp), INTENT(IN) :: shockVel,gasTemp,density,timeDelta
      REAL(dp) :: sputterRate,abundChangeFrac,totalMantle, grainNumberDensity
      INTEGER :: iSpec

      !Constant relating mass and speed of projectile to energy
      sConst=(shockVel*shockVel*km*km)/(2.0*gasTemp*K_BOLTZ)
      sConst=sqrt(sConst)

      !loop over projectile species and get rates of change of mantle for each, summing them
      sputterRate=0.0
      DO iSpec=1,SIZE(projectiles) !!!! Make projectiles array in initialize
          sputterRate=sputterRate+iceYieldRate(mass(projectiles(iSpec))*MH,density*abund(projectiles(iSpec)),gasTemp)
      END DO

      grainNumberDensity=density/GAS_DUST_DENSITY_RATIO
      !Total rate/cm3 (ie released particles /cm3/s) is sputterRate (per grain) multiplied by grain number density
      sputterRate=sputterRate*grainNumberDensity

      !integrate that forward from currentTimeOld to currentTime. to get total number of particles released
      abundChangeFrac=sputterRate*(timeDelta)!/density
      !I think that commented out dens is required for units. However, sputtering doesn't happen if it is uncommented
      !and sputtering matches Jimenez-Serra et al. 2008 curves when it's commented out.


      !if M particles are released and there are N particles on the grain total
      !then a species with X particles on the grain will release M*(X/N)
      !this is M/N and we'll multiply by X below
      totalMantle=sum(abund(iceList))
      abundChangeFrac=abundChangeFrac/totalMantle
      if (abundChangeFrac .gt. 1.0d0) abundChangeFrac=1.0d0
      if (abundChangeFrac .lt. 0.00d0) abundChangeFrac=0.0d0

    write(87,*) timeInYears,shockVel,abundChangeFrac,timeDelta/SECONDS_PER_YEAR
      !multiply M/N by x and add to gas phase
      if (shockVel .ge. VAPORIZE_SPEED) THEN
        abund(gasIceList)=abund(gasIceList)+abundChangeFrac*abund(iceList)
        abund(iceList)=abund(iceList)-abundChangeFrac*abund(iceList)
      ELSE
        abund(gasSputters)=abund(gasSputters)+abundChangeFrac*abund(sputters)
        abund(sputters)=abund(sputters)-abundChangeFrac*abund(sputters)
      END IF
  END SUBROUTINE

  !Function calculates rate of change of ice mantle abundance of a species!
  !due to the impact of molecules of a given mass. actual rate is         !
  !proportional to projectile abundance                                   !
  FUNCTION iceYieldRate(projectileMass,projectileDensity,gasTemp)
      REAL(dp) :: iceYieldRate
      REAL(dp) :: projectileMass,projectileDensity,gasTemp
      REAL(dp) :: lowerLimit,upperLimit,s

      REAL(dp), PARAMETER :: iceBindingEnergy=0.53*1.6d-12,targetMass=18.0*MH   
      REAL(dp), PARAMETER :: iceYieldEfficiency=0.8 !

      !eta is effectively reduced mass of the collision
      eta=4.*iceYieldEfficiency*projectileMass*targetMass*((projectileMass+targetMass)**(-2.0))
      epso=max(1.,4.*eta)
      s=sConst*sqrt(projectileMass)

      !Lower limit is xth in Jimenez-Serra et al. 2008
      lowerLimit=sqrt(epso*iceBindingEnergy/(eta*K_BOLTZ*gasTemp))

      !Upper limit is just where the integrand goes to zero
      upperLimit=iceYieldIntegralLimit(lowerLimit,projectileMass,gasTemp)
      !calculate eq B.1 from Jimenez-Serra et al. 2008
      IF ((upperlimit-lowerLimit) .gt. 1d-4) THEN
          !first get integral from Eq B.1 including 1/s factor
          iceYieldRate=trapezoidIntegrate(lowerLimit,upperLimit,projectileMass,gasTemp)/s
          !multiply through by constants
          iceYieldRate=iceYieldRate*GRAIN_RADIUS*GRAIN_RADIUS*sqrt(8.0*K_BOLTZ*gasTemp*pi/projectileMass)
          !need projectile number density
          iceYieldRate=iceYieldRate*projectileDensity
      ELSE
          iceYieldRate=0.0
      ENDIF
  END FUNCTION


  !Function calculates integrand from Eq B.1 of Jimenez-Serra et al. 2008 !
  !                                                                       !
  !Inputs are mass of projectile and x. Returns value of integrand at x   !
  !allowing trapezium rule to integrate from xth to infinity              !
  FUNCTION iceYieldIntegrand(x,projectileMass,gasTemp)
      REAL(dp) :: iceYieldIntegrand,x,projectileMass,gasTemp
      REAL(dp) :: yield,s,eps

      REAL(dp), PARAMETER :: yieldConst=8.3d-4
      REAL(dp), PARAMETER :: iceBindingEnergy=0.53*1.6d-12

      !this is s from exp(x+s) in eq B.1, varies only with mass so constant precalculated in initialize
      s=sConst*sqrt(projectileMass)

      !epsilon is calculated from inmpact energy (Ep)
      eps=(x**2)*K_BOLTZ*gasTemp
      !and some other factors
      eps=eta*eps/iceBindingEnergy
      !this yield is for ice averaged over all angles. There's a different one for cores (Appendix B Jimenez-Serra 2008)
      !it's 2 times the normal incidence yield, but there's a factor of 0.5 in integrand so we drop both
      yield=yieldConst*((eps-epso)**2)/(1.+((eps/30.)**(1.3333)))
      iceYieldIntegrand=yield*(x**2)*(DEXP(-((x-s)**2))-DEXP(-((x+s)**2)))
  END FUNCTION iceYieldIntegrand

  !Function to calculate the upper limit beyond which there's no point   
  !evaluating the ice yield integrand. Ie trapezoids from upper limit to 
  !upperlimit+dx will have an area~0                                     
  FUNCTION iceYieldIntegralLimit(xth,projectileMass,gasTemp)
      REAL(dp) :: iceYieldIntegralLimit
      REAL(dp) :: xth,projectileMass,gasTemp
      INTEGER :: i
      i=1
      !Take upperlimit to be half way between lowerLimit and 1000.
      iceYieldIntegralLimit=xth+(1d3-xth)*(0.5**i)
      !decrease upper limit for as long as f(upperlimit) is <1.0e-20 and 
      !difference between lower and upper limit is not zero.
      DO WHILE (iceYieldIntegrand(iceYieldIntegralLimit,projectileMass,gasTemp) .lt. 1d-200 .and.&
          & (iceYieldIntegralLimit-xth) .gt. 1.0d-3)
              i=i+1
              iceYieldIntegralLimit=xth+(1d3-xth)*(0.5**i)
      END DO
  END FUNCTION iceYieldIntegralLimit



  Function trapezoidIntegrate(lowerLimit,upperlimit,projectileMass,gasTemp)
     !Subroutine that calculates an integral using the trapezoidal method. It repeatedly
    !tries smaller intervals until the area is small enough to be accurate.
    ! It used to take a function to integrate as an argument but I removed it
    ! since we just want to integrate iceYieldIntegrand anyway.
      REAL(dp) :: trapezoidIntegrate
      INTEGER,PARAMETER :: JMAX=25
      REAL(dp) lowerLimit,upperlimit,projectileMass,gasTemp
      REAL(dp), PARAMETER :: tolerance=1.e-3
      INTEGER j
      REAL(dp) olds
      olds=-1.e30
      DO j=1,JMAX
          call trapzd(lowerLimit,upperlimit,trapezoidIntegrate,j,projectileMass,gasTemp)
          if (abs(trapezoidIntegrate-olds).le.tolerance*abs(olds)) RETURN
          olds=trapezoidIntegrate
      END DO
  END Function trapezoidIntegrate

  SUBROUTINE trapzd(a,b,s,n,projectileMass,gasTemp)
    ! Subroutine to integrate with trapazoidal rule using n intervals.
      INTEGER n
      REAL(dp) a,b,s,projectileMass,gasTemp
      INTEGER it,j
      REAL(dp) del,sum,tnm,x
      IF(n.eq.1) THEN
          s=0.5*(b-a)*(iceYieldIntegrand(a,projectileMass,gasTemp)+&
          &iceYieldIntegrand(b,projectileMass,gasTemp))
      ELSE
          it=2**(n-2)
          tnm=it
          del=(b-a)/tnm
          x=a+0.5*del
          sum=0.
          DO j=1,it
              sum=sum+iceYieldIntegrand(x,projectileMass,gasTemp)
              x=x+del
          END DO
          s=0.5*(s+(b-a)*sum/tnm)
      ENDIF
      RETURN
  END SUBROUTINE trapzd
END MODULE sputtering
