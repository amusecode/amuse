c*************************************************************************
c                        SWIFT.INC
c*************************************************************************
C Include file for SWIFT
c
c Author:  Hal Levison  
c Date:    2/2/93
c Last revision: 3/7/93

        implicit NONE    ! you got it baby

c...   Maximum array size
       integer  NPLMAX, NTPMAX
       parameter  (NPLMAX = 202)   ! max number of planets, including the Sun 
       parameter  (NTPMAX = 2000) ! max number of test particles

c...   Size of the test particle status flag
	integer NSTAT
        parameter  (NSTAT = 3)  

c...   convergence criteria for danby
	real*8 DANBYAC , DANBYB
	parameter (DANBYAC= 1.0d-14, DANBYB = 1.0d-13)

c...    loop limits in the Laguerre attempts
	integer NLAG1, NLAG2
	parameter(NLAG1 = 50, NLAG2 = 400)

c...    A small number
	real*8 TINY
        PARAMETER(TINY=4.D-15)

c...    trig stuff
	real*8 PI,TWOPI,PIBY2,DEGRAD
	parameter (PI = 3.141592653589793D0)
	parameter (TWOPI = 2.0D0 * PI)
	parameter (PIBY2 = PI/2.0D0)
	parameter (DEGRAD = 180.0D0 / PI)

c-------------------------------------------------------------------------
