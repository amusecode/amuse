MODULE CONSTANTS
   use, intrinsic :: iso_fortran_env, dp=>real64 !define the size of our double precision numbers
   REAL(dp), parameter :: C  = 2.99792458D+10 !Speed of light in cgs
   REAL(dp), PARAMETER :: K_BOLTZ = 1.38065040D-16 ! Boltzmann constant cgs
   REAL(dp), PARAMETER :: HP = 6.62606896D-27 !Planck constant in cgs
   REAL(dp), PARAMETER :: REDUCED_PLANCK=1.054571628d-27
   REAL(dp), PARAMETER :: MH = 1.67262164D-24 !H nucleus mass in cgs
   REAL(dp), PARAMETER :: AMU=1.66053892d-24 !atomic mass unit in cgs
   REAL(dp), PARAMETER :: PI = 3.141592654
   REAL(dp), PARAMETER :: K_BOLTZ_SI=1.38d-23 !Boltzmann constant SI
   REAL(dp), PARAMETER :: PC=3.086d18 !parsec in cgs
   REAL(dp), PARAMETER :: au=2.063d5 !1 AU in cgs
   REAL(dp), PARAMETER :: KM=1.d5 !kilometre in cgs
   REAL(dp), PARAMETER :: SECONDS_PER_YEAR=3.16d7
   REAL(dp), PARAMETER :: T_CMB=2.73
   REAL(dp), PARAMETER :: EV = 1.60217646D-12 ! electron volt in erg
   REAL(dp), PARAMETER :: GRAV_G = 6.674d-8 !gravitational constant in cgs
   REAL(dp), PARAMETER :: SB_CONST=5.6704d-5 !Stefan Boltzmann constant in cgs

   !Error codes for python wrap
   INTEGER, PARAMETER :: PARAMETER_READ_ERROR=-1
   INTEGER, PARAMETER :: PHYSICS_INIT_ERROR=-2
   INTEGER, PARAMETER :: CHEM_INIT_ERROR=-3
   INTEGER, PARAMETER :: INT_UNRECOVERABLE_ERROR=-4
   INTEGER, PARAMETER :: INT_TOO_MANY_FAILS_ERROR=-5
CONTAINS
   !Hold over from heating branch
    SUBROUTINE pair_insertion_sort(array)
    REAL(dp), INTENT(inout) :: array(:)
    INTEGER :: i,j,last
    REAL(dp) :: t1,t2

    last=size(array)
    DO i=2,last-1,2
       t1=min(array(i),array(i+1))
       t2=max(array(i),array(i+1))
       j=i-1
       DO while((j.ge.1).and.(array(j).gt.t2))
          array(j+2)=array(j)
          j=j-1
       ENDDO
       array(j+2)=t2
       DO while((j.ge.1).and.(array(j).gt.t1))
          array(j+1)=array(j)
          j=j-1
       ENDDO
       array(j+1)=t1
    END DO

    IF(mod(last,2).eq.0)then
       t1=array(last)
       DO j=last-1,1,-1
          IF (array(j).le.t1) exit
          array(j+1)=array(j)
       END DO
       array(j+1)=t1
    ENDIF

  END SUBROUTINE pair_insertion_sort

END MODULE CONSTANTS