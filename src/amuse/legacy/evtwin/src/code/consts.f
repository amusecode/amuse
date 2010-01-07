! This module contains all mathematical and physical constants
! This does not include constants which are configurable through input files,
      module constants
      
      ! Mathematical constants: pi, 4pi, log(10) and 1/3
      double precision, parameter :: CPI = 3.1415926535897932384626433832795029D0
      double precision, parameter :: CPI4 = 4.0D0 * CPI
      double precision, parameter :: CLN = 2.3025850929940456840179914546843642D0
      double precision, parameter :: C3RD = 1.0D0/3.0D0

      ! Physical constants
      double precision, parameter :: CL = 2.99792458D10     ! Speed of light
      double precision, parameter :: PLANCK = 6.6260755D-27 
      double precision, parameter :: CG = 6.672D-8          ! Newton's constant
      double precision, parameter :: BOLTZM = 1.380658D-16 
      double precision, parameter :: ECHAR = 4.8032068D-10  ! Electrostatic charge

      double precision, parameter :: AMU = 1.6605402D-24    ! Atomic mass unit
      double precision, parameter :: AME = 9.1093897D-28    ! Electron mass

      ! Astrophysical parameters
      double precision, parameter :: CMSN = 1.9891D0        ! Solar mass
      double precision, parameter :: CLSN = 3.844D0         ! Solar luminosity
      double precision, parameter :: CRSN = 0.69598D0       ! Solar radius
      double precision, parameter :: CASN = 4.57D9          ! Solar age
      double precision, parameter :: CZSN = 0.02            ! solar metallicity 

      ! Units
      double precision, parameter :: CSY = 3.155692597D7    ! Seconds per year
      double precision, parameter :: CSDAY = 86400D0        ! Seconds per day
      double precision, parameter :: EVOLT = ECHAR/CL*1.0D8 ! erg/electron volt

      ! Derived constants
      double precision, parameter :: CA = 8.0D0*CPI**5*BOLTZM**4 / (15.0D0 * (CL*PLANCK)**3)
      double precision, parameter :: CME = 1.0D6*EVOLT/AMU  ! (erg/MeV) / (gram/amu)
      double precision, parameter :: CEVB = EVOLT/BOLTZM 
      double precision, parameter :: CR = BOLTZM/AMU        ! gas constant

      double precision, parameter :: LAMC = PLANCK/(AME*CL) ! Compton length of electron
      double precision, parameter :: CRHO = 8.0D0*CPI/LAMC**3

      double precision, parameter :: CB = CRHO*AME*CL**2 
      double precision, parameter :: CD = CRHO*AMU 
      
      double precision, parameter :: CTE = BOLTZM/(AME*CL**2) 
      double precision, parameter :: CGRT = 6.4D0*(CPI/4.32D4)**6*(1.0D11/CL)**5 

      ! Constants that cannot be computed at compile time
      double precision :: CEN
      double precision :: CPL
      double precision :: CG1
      double precision :: CG2
      
      contains

      ! Initialise collision integrals for atomic diffusion
      subroutine initialise_collision_integrals
      implicit none
      double precision :: dc(50*4*3), dd(50*4)
      common /colint/ dc, dd

      DC = (/
     &   +1.19599D-02,-2.39198D-02,-3.02547D+01,-2.94860D+01,
     &   -2.39198D-02,-1.48010D-02,-2.94860D+01,-2.87231D+01,
     &   -1.48010D-02,-1.77390D-02,-2.87231D+01,-2.79637D+01,
     &   -1.77390D-02,-1.74423D-02,-2.79637D+01,-2.72086D+01,
     &   -1.74423D-02,-1.80040D-02,-2.72086D+01,-2.64576D+01,
     &   -1.80040D-02,-1.83218D-02,-2.64576D+01,-2.57110D+01,
     &   -1.83218D-02,-1.86847D-02,-2.57110D+01,-2.49688D+01,
     &   -1.86847D-02,-1.90073D-02,-2.49688D+01,-2.42310D+01,
     &   -1.90073D-02,-1.93026D-02,-2.42310D+01,-2.34978D+01,
     &   -1.93026D-02,-1.95555D-02,-2.34978D+01,-2.27693D+01,
     &   -1.95555D-02,-1.97557D-02,-2.27693D+01,-2.20454D+01,
     &   -1.97557D-02,-1.98886D-02,-2.20454D+01,-2.13263D+01,
     &   -1.98886D-02,-1.99373D-02,-2.13263D+01,-2.06120D+01,
     &   -1.99373D-02,-1.98810D-02,-2.06120D+01,-1.99024D+01,
     &   -1.98810D-02,-1.96948D-02,-1.99024D+01,-1.91976D+01,
     &   -1.96948D-02,-1.93486D-02,-1.91976D+01,-1.84975D+01,
     &   -1.93486D-02,-1.88059D-02,-1.84975D+01,-1.78021D+01,
     &   -1.88059D-02,-1.80227D-02,-1.78021D+01,-1.71112D+01,
     &   -1.80227D-02,-1.69459D-02,-1.71112D+01,-1.64246D+01,
     &   -1.69459D-02,-1.55109D-02,-1.64246D+01,-1.57421D+01,
     &   -1.55109D-02,-1.36394D-02,-1.57421D+01,-1.50633D+01,
     &   -1.36394D-02,-1.12361D-02,-1.50633D+01,-1.43878D+01,
     &   -1.12361D-02,-8.18466D-03,-1.43878D+01,-1.37150D+01,
     &   -8.18466D-03,-4.34258D-03,-1.37150D+01,-1.30441D+01,
     &   -4.34258D-03,+4.65253D-04,-1.30441D+01,-1.23743D+01,
     &   +4.65253D-04,+6.45493D-03,-1.23743D+01,-1.17044D+01,
     &   +6.45493D-03,+1.38941D-02,-1.17044D+01,-1.10329D+01,
     &   +1.38941D-02,+2.31151D-02,-1.10329D+01,-1.03581D+01,
     &   +2.31151D-02,+3.45317D-02,-1.03581D+01,-9.67777D+00,
     &   +3.45317D-02,+4.86585D-02,-9.67777D+00,-8.98913D+00,
     &   +4.86585D-02,+6.61321D-02,-8.98913D+00,-8.28881D+00,
     &   +6.61321D-02,+8.77309D-02,-8.28881D+00,-7.57261D+00,
     &   +8.77309D-02,+1.14383D-01,-7.57261D+00,-6.83537D+00,
     &   +1.14383D-01,+1.47142D-01,-6.83537D+00,-6.07066D+00,
     &   +1.47142D-01,+1.87092D-01,-6.07066D+00,-5.27065D+00,
     &   +1.87092D-01,+2.35096D-01,-5.27065D+00,-4.42573D+00,
     &   +2.35096D-01,+2.91268D-01,-4.42573D+00,-3.52439D+00,
     &   +2.91268D-01,+3.53977D-01,-3.52439D+00,-2.55315D+00,
     &   +3.53977D-01,+4.18217D-01,-2.55315D+00,-1.49695D+00,
     &   +4.18217D-01,+4.73499D-01,-1.49695D+00,-3.40379D-01,
     &   +4.73499D-01,+5.02343D-01,-3.40379D-01,+9.29832D-01,
     &   +5.02343D-01,+4.82140D-01,+9.29832D-01,+2.32060D+00,
     &   +4.82140D-01,+3.92303D-01,+2.32060D+00,+3.82709D+00,
     &   +3.92303D-01,+2.20401D-01,+3.82709D+00,+5.42773D+00,
     &   +2.20401D-01,-5.31156D-02,+5.42773D+00,+7.08127D+00,
     &   -5.31156D-02,-3.94063D-01,+7.08127D+00,+8.72205D+00,
     &   -3.94063D-01,-5.99574D-01,+8.72205D+00,+1.02683D+01,
     &   -5.99574D-01,-4.71033D-01,+1.02683D+01,+1.16706D+01,
     &   -4.71033D-01,-4.68969D-01,+1.16706D+01,+1.29598D+01,
     &   -4.68969D-01,+2.34484D-01,+1.29598D+01,+1.41366D+01,
     &
     &   +1.34102D-02,-2.68205D-02,-2.55941D+01,-2.48408D+01,
     &   -2.68205D-02,-1.66309D-02,-2.48408D+01,-2.40939D+01,
     &   -1.66309D-02,-1.99547D-02,-2.40939D+01,-2.33511D+01,
     &   -1.99547D-02,-1.96575D-02,-2.33511D+01,-2.26130D+01,
     &   -1.96575D-02,-2.03264D-02,-2.26130D+01,-2.18796D+01,
     &   -2.03264D-02,-2.07272D-02,-2.18796D+01,-2.11511D+01,
     &   -2.07272D-02,-2.11843D-02,-2.11511D+01,-2.04276D+01,
     &   -2.11843D-02,-2.16034D-02,-2.04276D+01,-1.97091D+01,
     &   -2.16034D-02,-2.20005D-02,-1.97091D+01,-1.89959D+01,
     &   -2.20005D-02,-2.23602D-02,-1.89959D+01,-1.82879D+01,
     &   -2.23602D-02,-2.26732D-02,-1.82879D+01,-1.75853D+01,
     &   -2.26732D-02,-2.29256D-02,-1.75853D+01,-1.68881D+01,
     &   -2.29256D-02,-2.31017D-02,-1.68881D+01,-1.61965D+01,
     &   -2.31017D-02,-2.31822D-02,-1.61965D+01,-1.55103D+01,
     &   -2.31822D-02,-2.31439D-02,-1.55103D+01,-1.48298D+01,
     &   -2.31439D-02,-2.29591D-02,-1.48298D+01,-1.41548D+01,
     &   -2.29591D-02,-2.25945D-02,-1.41548D+01,-1.34853D+01,
     &   -2.25945D-02,-2.20100D-02,-1.34853D+01,-1.28212D+01,
     &   -2.20100D-02,-2.11574D-02,-1.28212D+01,-1.21624D+01,
     &   -2.11574D-02,-1.99787D-02,-1.21624D+01,-1.15087D+01,
     &   -1.99787D-02,-1.84040D-02,-1.15087D+01,-1.08598D+01,
     &   -1.84040D-02,-1.63488D-02,-1.08598D+01,-1.02153D+01,
     &   -1.63488D-02,-1.37108D-02,-1.02153D+01,-9.57474D+00,
     &   -1.37108D-02,-1.03658D-02,-9.57474D+00,-8.93745D+00,
     &   -1.03658D-02,-6.16232D-03,-8.93745D+00,-8.30266D+00,
     &   -6.16232D-03,-9.15489D-04,-8.30266D+00,-7.66934D+00,
     &   -9.15489D-04,+5.60137D-03,-7.66934D+00,-7.03625D+00,
     &   +5.60137D-03,+1.36664D-02,-7.03625D+00,-6.40181D+00,
     &   +1.36664D-02,+2.36209D-02,-6.40181D+00,-5.76409D+00,
     &   +2.36209D-02,+3.58820D-02,-5.76409D+00,-5.12070D+00,
     &   +3.58820D-02,+5.09546D-02,-5.12070D+00,-4.46870D+00,
     &   +5.09546D-02,+6.94376D-02,-4.46870D+00,-3.80447D+00,
     &   +6.94376D-02,+9.20141D-02,-3.80447D+00,-3.12357D+00,
     &   +9.20141D-02,+1.19403D-01,-3.12357D+00,-2.42059D+00,
     &   +1.19403D-01,+1.52233D-01,-2.42059D+00,-1.68896D+00,
     &   +1.52233D-01,+1.90763D-01,-1.68896D+00,-9.20788D-01,
     &   +1.90763D-01,+2.34333D-01,-9.20788D-01,-1.06834D-01,
     &   +2.34333D-01,+2.80392D-01,-1.06834D-01,+7.63360D-01,
     &   +2.80392D-01,+3.23003D-01,+7.63360D-01,+1.70085D+00,
     &   +3.23003D-01,+3.51091D-01,+1.70085D+00,+2.71586D+00,
     &   +3.51091D-01,+3.47669D-01,+2.71586D+00,+3.81513D+00,
     &   +3.47669D-01,+2.93155D-01,+3.81513D+00,+4.99784D+00,
     &   +2.93155D-01,+1.77064D-01,+4.99784D+00,+6.25091D+00,
     &   +1.77064D-01,+1.68817D-02,+6.25091D+00,+7.54647D+00,
     &   +1.68817D-02,-1.37836D-01,+7.54647D+00,+8.84609D+00,
     &   -1.37836D-01,-2.36519D-01,+8.84609D+00,+1.01126D+01,
     &   -2.36519D-01,-2.49035D-01,+1.01126D+01,+1.13224D+01,
     &   -2.49035D-01,-1.95933D-01,+1.13224D+01,+1.24724D+01,
     &   -1.95933D-01,-1.60453D-01,+1.24724D+01,+1.35754D+01,
     &   -1.60453D-01,+8.02267D-02,+1.35754D+01,+1.46398D+01,
     &
     &   +1.46130D-02,-2.92259D-02,-1.93212D+01,-1.85804D+01,
     &   -2.92259D-02,-1.81472D-02,-1.85804D+01,-1.78467D+01,
     &   -1.81472D-02,-2.17897D-02,-1.78467D+01,-1.71173D+01,
     &   -2.17897D-02,-2.14905D-02,-1.71173D+01,-1.63931D+01,
     &   -2.14905D-02,-2.22463D-02,-1.63931D+01,-1.56741D+01,
     &   -2.22463D-02,-2.27134D-02,-1.56741D+01,-1.49604D+01,
     &   -2.27134D-02,-2.32457D-02,-1.49604D+01,-1.42522D+01,
     &   -2.32457D-02,-2.37412D-02,-1.42522D+01,-1.35496D+01,
     &   -2.37412D-02,-2.42183D-02,-1.35496D+01,-1.28526D+01,
     &   -2.42183D-02,-2.46613D-02,-1.28526D+01,-1.21615D+01,
     &   -2.46613D-02,-2.50615D-02,-1.21615D+01,-1.14763D+01,
     &   -2.50615D-02,-2.54057D-02,-1.14763D+01,-1.07971D+01,
     &   -2.54057D-02,-2.56787D-02,-1.07971D+01,-1.01240D+01,
     &   -2.56787D-02,-2.58623D-02,-1.01240D+01,-9.45709D+00,
     &   -2.58623D-02,-2.59347D-02,-9.45709D+00,-8.79637D+00,
     &   -2.59347D-02,-2.58697D-02,-8.79637D+00,-8.14187D+00,
     &   -2.58697D-02,-2.56362D-02,-8.14187D+00,-7.49358D+00,
     &   -2.56362D-02,-2.51970D-02,-7.49358D+00,-6.85144D+00,
     &   -2.51970D-02,-2.45075D-02,-6.85144D+00,-6.21534D+00,
     &   -2.45075D-02,-2.35146D-02,-6.21534D+00,-5.58513D+00,
     &   -2.35146D-02,-2.21544D-02,-5.58513D+00,-4.96057D+00,
     &   -2.21544D-02,-2.03503D-02,-4.96057D+00,-4.34132D+00,
     &   -2.03503D-02,-1.80104D-02,-4.34132D+00,-3.72695D+00,
     &   -1.80104D-02,-1.50236D-02,-3.72695D+00,-3.11691D+00,
     &   -1.50236D-02,-1.12556D-02,-3.11691D+00,-2.51047D+00,
     &   -1.12556D-02,-6.54378D-03,-2.51047D+00,-1.90674D+00,
     &   -6.54378D-03,-6.90518D-04,-1.90674D+00,-1.30457D+00,
     &   -6.90518D-04,+6.54438D-03,-1.30457D+00,-7.02574D-01,
     &   +6.54438D-03,+1.54523D-02,-7.02574D-01,-9.90038D-02,
     &   +1.54523D-02,+2.63848D-02,-9.90038D-02,+5.08275D-01,
     &   +2.63848D-02,+3.97599D-02,+5.08275D-01,+1.12189D+00,
     &   +3.97599D-02,+5.60625D-02,+1.12189D+00,+1.74504D+00,
     &   +5.60625D-02,+7.58259D-02,+1.74504D+00,+2.38165D+00,
     &   +7.58259D-02,+9.95773D-02,+2.38165D+00,+3.03645D+00,
     &   +9.95773D-02,+1.27708D-01,+3.03645D+00,+3.71516D+00,
     &   +1.27708D-01,+1.60207D-01,+3.71516D+00,+4.42451D+00,
     &   +1.60207D-01,+1.96165D-01,+4.42451D+00,+5.17232D+00,
     &   +1.96165D-01,+2.32939D-01,+5.17232D+00,+5.96720D+00,
     &   +2.32939D-01,+2.64985D-01,+5.96720D+00,+6.81799D+00,
     &   +2.64985D-01,+2.82691D-01,+6.81799D+00,+7.73238D+00,
     &   +2.82691D-01,+2.72481D-01,+7.73238D+00,+8.71461D+00,
     &   +2.72481D-01,+2.20588D-01,+8.71461D+00,+9.76224D+00,
     &   +2.20588D-01,+1.23271D-01,+9.76224D+00,+1.08628D+01,
     &   +1.23271D-01,+5.64105D-05,+1.08628D+01,+1.19930D+01,
     &   +5.64105D-05,-1.01669D-01,+1.19930D+01,+1.31231D+01,
     &   -1.01669D-01,-1.45351D-01,+1.31231D+01,+1.42289D+01,
     &   -1.45351D-01,-1.12670D-01,+1.42289D+01,+1.52998D+01,
     &   -1.12670D-01,-1.25510D-01,+1.52998D+01,+1.63436D+01,
     &   -1.25510D-01,+4.14471D-02,+1.63436D+01,+1.73573D+01,
     &   +4.14471D-02,-2.07236D-02,+1.73573D+01,+1.83810D+01
     &/)

      DD = (/
     &   +1.18229D-02,-2.36458D-02,-2.55112D+01,-2.47319D+01,
     &   -2.36458D-02,-1.46794D-02,-2.47319D+01,-2.39583D+01,
     &   -1.46794D-02,-1.76226D-02,-2.39583D+01,-2.31882D+01,
     &   -1.76226D-02,-1.73748D-02,-2.31882D+01,-2.24223D+01,
     &   -1.73748D-02,-1.79780D-02,-2.24223D+01,-2.16606D+01,
     &   -1.79780D-02,-1.83439D-02,-2.16606D+01,-2.09032D+01,
     &   -1.83439D-02,-1.87580D-02,-2.09032D+01,-2.01502D+01,
     &   -1.87580D-02,-1.91359D-02,-2.01502D+01,-1.94017D+01,
     &   -1.91359D-02,-1.94907D-02,-1.94017D+01,-1.86579D+01,
     &   -1.94907D-02,-1.98074D-02,-1.86579D+01,-1.79186D+01,
     &   -1.98074D-02,-2.00758D-02,-1.79186D+01,-1.71842D+01,
     &   -2.00758D-02,-2.02817D-02,-1.71842D+01,-1.64545D+01,
     &   -2.02817D-02,-2.04080D-02,-1.64545D+01,-1.57298D+01,
     &   -2.04080D-02,-2.04344D-02,-1.57298D+01,-1.50099D+01,
     &   -2.04344D-02,-2.03360D-02,-1.50099D+01,-1.42949D+01,
     &   -2.03360D-02,-2.00828D-02,-1.42949D+01,-1.35848D+01,
     &   -2.00828D-02,-1.96388D-02,-1.35848D+01,-1.28796D+01,
     &   -1.96388D-02,-1.89603D-02,-1.28796D+01,-1.21790D+01,
     &   -1.89603D-02,-1.79946D-02,-1.21790D+01,-1.14830D+01,
     &   -1.79946D-02,-1.66781D-02,-1.14830D+01,-1.07913D+01,
     &   -1.66781D-02,-1.49334D-02,-1.07913D+01,-1.01036D+01,
     &   -1.49334D-02,-1.26670D-02,-1.01036D+01,-9.41952D+00,
     &   -1.26670D-02,-9.76488D-03,-9.41952D+00,-8.73845D+00,
     &   -9.76488D-03,-6.08817D-03,-8.73845D+00,-8.05974D+00,
     &   -6.08817D-03,-1.46688D-03,-8.05974D+00,-7.38248D+00,
     &   -1.46688D-03,+4.30758D-03,-7.38248D+00,-6.70557D+00,
     &   +4.30758D-03,+1.14919D-02,-6.70557D+00,-6.02763D+00,
     &   +1.14919D-02,+2.04026D-02,-6.02763D+00,-5.34694D+00,
     &   +2.04026D-02,+3.14302D-02,-5.34694D+00,-4.66134D+00,
     &   +3.14302D-02,+4.50547D-02,-4.66134D+00,-3.96821D+00,
     &   +4.50547D-02,+6.18616D-02,-3.96821D+00,-3.26426D+00,
     &   +6.18616D-02,+8.25526D-02,-3.26426D+00,-2.54546D+00,
     &   +8.25526D-02,+1.07939D-01,-2.54546D+00,-1.80685D+00,
     &   +1.07939D-01,+1.38897D-01,-1.80685D+00,-1.04234D+00,
     &   +1.38897D-01,+1.76234D-01,-1.04234D+00,-2.44487D-01,
     &   +1.76234D-01,+2.20387D-01,-2.44487D-01,+5.95659D-01,
     &   +2.20387D-01,+2.70814D-01,+5.95659D-01,+1.48870D+00,
     &   +2.70814D-01,+3.24880D-01,+1.48870D+00,+2.44673D+00,
     &   +3.24880D-01,+3.76096D-01,+2.44673D+00,+3.48274D+00,
     &   +3.76096D-01,+4.11895D-01,+3.48274D+00,+4.60901D+00,
     &   +4.11895D-01,+4.12194D-01,+4.60901D+00,+5.83413D+00,
     &   +4.12194D-01,+3.52176D-01,+5.83413D+00,+7.15818D+00,
     &   +3.52176D-01,+2.14860D-01,+7.15818D+00,+8.56675D+00,
     &   +2.14860D-01,+1.56932D-02,+8.56675D+00,+1.00269D+01,
     &   +1.56932D-02,-1.84529D-01,+1.00269D+01,+1.14908D+01,
     &   -1.84529D-01,-3.17740D-01,+1.14908D+01,+1.29104D+01,
     &   -3.17740D-01,-3.51796D-01,+1.29104D+01,+1.42538D+01,
     &   -3.51796D-01,-2.86379D-01,+1.42538D+01,+1.55127D+01,
     &   -2.86379D-01,-2.68214D-01,+1.55127D+01,+1.67029D+01,
     &   -2.68214D-01,+1.34107D-01,+1.67029D+01,+1.78287D+01
     &/)

      end subroutine
      
      subroutine initialise_constants
      implicit none
      call initialise_collision_integrals
      CEN = (PLANCK/(2.0D0*CPI*AMU*BOLTZM)**0.5D0)**3 / AMU 
      CPL = (CPI4/(AMU*BOLTZM**3))**0.5D0 * ECHAR**3 
      CG1 = 1.0D5*CG**0.5D0 
      CG2 = CG1*(CSDAY*1.0D-5) / (2.0D0 * CPI) 
      end subroutine
      
      end module

! This module defines all run parameters and options
      module settings 

      double precision :: CZS
      double precision :: CH = -1.0
      double precision :: CDC(10)
      double precision :: CT1, CT2, CT3
      double precision :: INITIAL_CT(10)
      double precision :: CT(10)
      double precision :: CC, CN, CO, CNE, CMG, CSI, CFE
     
      double precision :: CALP, CU, COS, CPS, CRD, CTH
      double precision :: CXB, CGR
      double precision :: CEA, CET
      double precision :: CMT, CMS, CMI, CML, CHL, CTF, CLT
      double precision :: CPA, CBR, CSU, CSD, CDF, CGW, CSO, CMB
      
      double precision :: CLIMIT = 1.0D-1    ! Limit changes in variables during iterations

      ! Desired accuracy. The solver will aim for an accuracy between
      !   EPS (in the unnamed COMMON block) and WANTED_EPS
      ! No effect if WANTED_EPS <= EPS
      double precision :: wanted_eps = 1.0D-8

      double precision :: CPHOTONTIRE = 0.0  ! Switch to include photon tiring

      ! Emergency energy generation term, normally set to 0.
      ! This cannot be set from the input file. It will be set by REMESH
      ! if there is NO nuclear energy generation in the initial model AT
      ! ALL. In that case, the first iteration(s) will return LOM = 0.0
      ! throughout the star because the thermal energy term is initially
      ! 0 as well. This is a numerical fudge to remove the resulting
      ! singularity. This term will be set to L/M (constant energy
      ! generation throughout the star) and will be reduced to 0 by printb.
      double precision :: ENC_PARACHUTE = 0.0
      
      ! Should the equation-of-state include the effects of pair production?
      ! This is only important in the very late burning stages of very massive
      ! stars. Positrons are only calculated if their degeneracy parameter
      ! >= -15.0 - otherwise they are negligible anyway.
      logical :: EOS_INCLUDE_PAIRPRODUCTION = .false.
      
      ! Turn on smart mass loss routine, that picks an appropriate
      ! recipe depending on the stellar parameters. This is an
      ! alternative for the de Jager rate and replaces it when
      ! SMART_MASS_LOSS is switched on. Off by default.
      double precision :: SMART_MASS_LOSS = 0.0

      ! Individual mass loss recipe switches.
      ! These also turn on recipes when SMART_MASS_LOSS is used,
      ! although that does store its own set of mass loss options (to
      ! keep it more modular).
      double precision :: CMR       ! Switch for Reimers-like mass loss rate
      double precision :: CMJ       ! Switch for de Jager mass loss rate
      double precision :: CMV       ! Switch for Vink mass loss rate
      double precision :: CMK       ! Switch for Kudritzki 2002 mass loss rate
      double precision :: CMNL      ! Switch for Nugis&Lamers mass loss rate (WR stars)
      double precision :: CMRR      ! Switch for Real Reimers mass loss rate
      double precision :: CMVW      ! Switch for Vasiliadis&Wood (AGB) rate
      double precision :: CMSC      ! Switch for Schroeder&Cuntz mass loss
      double precision :: CMW       ! Switch for Wachter&al (AGB) mass loss
      double precision :: CMAL      ! Switch for Achmad&Lamers (A supergiants)

      ! Rotationally enhanced mass loss rates, two options: Heger & al,
      ! Maeder & Meynet. Set one of these!
      double precision :: CMDOTROT_HLW ! Heger, Langer & Woosely
      double precision :: CMDOTROT_MM  ! Maeder & Meynet

      ! Non-conservative mass transfer options (depending on stellar parameters)
      double precision :: CMTEL     ! Eddington-limited accretion (0 or 1)
      double precision :: CMTWL     ! Angular momentum limited accretion

      ! Scaling with metallicity applied to de Jager mass loss rate in funcs1
      double precision :: ZSCALING_MDOT = 0.8

      double precision :: ARTMIX = 0.0D0  ! Artificial mixing coefficient [cm^2/s]

      double precision :: CCAC = 0.0D0 ! Switch for composition accretion

      double precision :: CGRS = 0.0D0 ! Switch for gravitational settling 

      double precision :: CSMC = 0.04D0! Semi-convection efficiency, after Langer (1991)

 
      double precision :: CDSI = 1.0D0 ! Switch for the Dynamical Shear Instability
      double precision :: CSHI = 1.0D0 ! Switch for the Solberg-Hoiland Instability (not implmented)
      double precision :: CSSI = 1.0D0 ! Switch for the Secular Shear Instability
      double precision :: CESC = 1.0D0 ! Switch for the Eddington-Sweet Circulation
      double precision :: CGSF = 1.0D0 ! Switch for the Goldreich-Schubert-Fricke instability

      double precision :: CFMU = 0.05D0 ! weight of mu gradient in rotational instabilities [see Hegers thesis page 36 and pinsonneault]
      double precision :: CFC = 1.0D0/30.0D0 ! ratio of tubulent viscosity over the diffusion coefficient [see Hegers thesis page 35]

     
      integer ::KTH, KX, KY, KZ
      integer ::KT1, KT2, KT3, KT4
      integer :: KCL, KION, KAM, KOP, KCC, KNUC, KCN
      integer :: KSX(45)
      integer :: KN, KJN(40)

! Number of allowed iterations for the nucleosynthesis code
      integer :: KR_NUCSYN = 60

c Variables derived from Settings and never changed
      double precision :: CLOGZ
      logical :: rigid_rotation = .true.     ! Rigid rotation or differential rotation?

! Switches for the new "smooth" remesher.
! The first two are set in init.dat, the last one in init.run
      logical :: use_smooth_remesher = .false.
      logical :: relax_loaded_model = .true.
      logical :: start_with_rigid_rotation = .true.

      ! Unused, but the code relies on this being defined:
      double precision :: CQ1, CQ(17) 
      
      end module

      MODULE OPACITY_CO
c     -----------------------------------------------------------------
c     This module contains all variables relating to the CO opacity
c     tables. Adapted from Eldridge implemented by SdM
c     -----------------------------------------------------------------
c     CO_MT:          opacity table size in log10 T direction
c     CO_MR:          "" in log R = rho^3/T6 direction
c     CO_MH:          "" in Hyrdrogen abundance direction
c     CO_CO:          "" in Carbon and oxygen direction 
c     SPLINE_OPAC_CO: Array for spline coefficients
c     CBASE, OBASE:   C and O abundance because of metallicity only
c     opT, opR:       log T and log R on which opacity table is defined
c     opH, COcompos:  Compositions on which opacity table is defined
c     -----------------------------------------------------------------

      INTEGER, PARAMETER :: CO_MH = 5  
      INTEGER, PARAMETER :: CO_MT = 141
      INTEGER, PARAMETER :: CO_MR = 31
      INTEGER, PARAMETER :: CO_CO = 8

      DOUBLE PRECISION :: CBASE,OBASE

      DOUBLE PRECISION :: opT(CO_MT), opR(CO_MR)
      DOUBLE PRECISION :: opH(CO_MH) = 
     &     (/0., 0.03, 0.1, 0.35, 0.7/)
      DOUBLE PRECISION :: COcompos(CO_CO) =  
     &     (/0.0d0, 0.01d0, 0.03d0, 0.1d0, 0.2d0, 0.4d0, 0.6d0, 1.0d0/)

      DOUBLE PRECISION, ALLOCATABLE :: SPLINE_OPAC_CO(:,:,:,:,:)

      END MODULE 
