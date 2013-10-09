!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MERCURY.INC    (ErikSoft   4 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Parameters that you may want to alter at some point:
!
! NMAX  = maximum number of bodies
! CMAX  = maximum number of close-encounter minima monitored simultaneously
! NMESS = maximum number of messages in message.in
! HUGE  = an implausibly large number
! NFILES = maximum number of files that can be open at the same time
!
      integer NMAX, CMAX, NMESS, NFILES
      real*8 HUGE
!
      parameter (NMAX = 2000)
      parameter (CMAX = 50)
      parameter (NMESS = 200)
      parameter (HUGE = 9.9d29)
      parameter (NFILES = 50)
!
!------------------------------------------------------------------------------
!
! Constants:
!
! DR = conversion factor from degrees to radians
! K2 = Gaussian gravitational constant squared
! AU = astronomical unit in cm
! MSUN = mass of the Sun in g
!
      real*8 PI,TWOPI,PIBY2,DR,K2,AU,MSUN
!
      parameter (PI = 3.141592653589793d0)
      parameter (TWOPI = PI * 2.d0)
      parameter (PIBY2 = PI * .5d0)
      parameter (DR = PI / 180.d0)
      parameter (K2 = 2.959872379828867d-04)  ! was: 2.959122082855911d-4
      parameter (AU = 1.4959787e13)
      parameter (MSUN = 1.988920e+33) ! was: 1.9891e33
