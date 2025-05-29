!> \file cosmology.F90

!> \brief Module contains simple analytic cosmological functions
!<
module cosmology_mod
use myf03_mod
implicit none

  contains

!-------------------------------------------------------------------
!> returns the time since the big bang in a flat LCDM Universe
!! the units of time are seconds.  
  function tsinceBB(a,OmegaM,h) result (t)
  use physical_constants_mod, only: km2Mpc
    real(r8b) :: a      !< scale factor
    real(r8b) :: OmegaM !< Omega Matter
    real(r8b) :: h      !< little hubble
    real(r8b) :: t      !< time since big bang in seconds
    real(r8b) :: OmegaL,aeq
    real(r8b) :: pre,arg,H0

    H0 = h * 100.0d0 * km2Mpc ! Hubble parameter in inverse seconds
    OmegaL = 1.0 - OmegaM
    aeq = (OmegaM/OmegaL)**(1./3.)

    pre = 2./( 3. * sqrt(OmegaL) )
    arg = (a/aeq)**(3./2.) + sqrt(1 + (a/aeq)**3)
    t = pre * log(arg) / H0

  end function tsinceBB


!-------------------------------------------------------------------
!> returns the physical time between two scale factors (seconds)
!! in a flat LCDM Universe.
  subroutine da2dt(a1,a2,OmegaM,h,dt) 
  use physical_constants_mod, only: km2Mpc

    real(r8b), intent(in) :: a1     !< initial scale factor
    real(r8b), intent(in) :: a2     !< final scale factor
    real(r8b), intent(in) :: OmegaM !< Omega Matter
    real(r8b), intent(in) :: h      !< little Hubble
    real(r8b), intent(out) :: dt     !< 
    
    real(r8b) :: t1,t2
    
    if (OmegaM <= 0.0) stop "cosmology.f90>> OmegaMatter must be >= 0.0"
    if (a1 <= 0.0 .or. a1 > 1.0) then
       write(*,*) "cosmology.f90>> da2dt - a1 out of range "
       write(*,*) "a1 = ", a1
       stop
    end if

    if (a2 <= 0.0 .or. a2 > 1.0) then
       write(*,*) "cosmology.f90>> da2dt - a2 out of range "
       write(*,*) "a2 = ", a2
       stop
    end if

    if (a1 <= 0.0 .or. a1 > 1.0) then
       write(*,*) "cosmology.f90>> da2dt - a1 out of range "
       write(*,*) "a1 = ", a1
       stop
    end if

    if (a2 < a1) then
       write(*,*) "cosmology.f90>> da2dt -  a2 must be > a1"
       write(*,*) "a1,a2 = ", a1, a2
       stop
    end if

    t1 = tsinceBB(a1,OmegaM,h)
    t2 = tsinceBB(a2,OmegaM,h)
    dt = t2-t1
 
  end subroutine da2dt



end module cosmology_mod


