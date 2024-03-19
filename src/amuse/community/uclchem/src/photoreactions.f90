MODULE photoreactions
USE constants
IMPLICIT NONE
    REAL(dp) :: UV_FAC=3.02

    !Below are arrays for self-shielding of CO and H2
    LOGICAL :: start=.True.
    INTEGER :: NUM_LAMBDA=30
    REAL(dp), DIMENSION(30) :: LAMBDA_GRID=(/ &
      &  910.0D0, 950.0D0,1000.0D0,1050.0D0,1110.0D0, &
      & 1180.0D0,1250.0D0,1390.0D0,1490.0D0,1600.0D0, &
      & 1700.0D0,1800.0D0,1900.0D0,2000.0D0,2100.0D0, &
      & 2190.0D0,2300.0D0,2400.0D0,2500.0D0,2740.0D0, &
      & 3440.0D0,4000.0D0,4400.0D0,5500.0D0,7000.0D0, &
      & 9000.0D0,12500.0D0,22000.0D0,34000.0D0,1.0D9/)
    REAL(dp), DIMENSION(30) :: XLAMBDA_GRID=(/ &
      & 5.76D0,5.18D0,4.65D0,4.16D0,3.73D0, &
      & 3.40D0,3.11D0,2.74D0,2.63D0,2.62D0, &
      & 2.54D0,2.50D0,2.58D0,2.78D0,3.01D0, &
      & 3.12D0,2.86D0,2.58D0,2.35D0,2.00D0, &
      & 1.58D0,1.42D0,1.32D0,1.00D0,0.75D0, &
      & 0.48D0,0.28D0,0.12D0,0.05D0,0.00D0/)
    REAL(dp), DIMENSION(30) :: XLAMBDA_DERIV

    logical :: startr=.True.

    !  12CO line shielding data from van Dishoeck & Black (1988, ApJ, 334, 771, Table 5)
    INTEGER, PARAMETER ::  DIMCO=7, DIMH2=6
    REAL(KIND=DP), DIMENSION(8) :: NCO_GRID=(/12.0D0,13.0D0,14.0D0,15.0D0,16.0D0,17.0D0,18.0D0,19.0D0/)
    REAL(KIND=DP), DIMENSION(6) :: NH2_GRID=(/18.0D0,19.0D0,20.0D0,21.0D0,22.0D0,23.0D0/)
    REAL(KIND=DP), DIMENSION(8,6) :: SCO_GRID=RESHAPE((/ &
      &  0.000D+00,-1.408D-02,-1.099D-01,-4.400D-01,-1.154D+00,-1.888D+00,-2.760D+00,-4.001D+00, &
      & -8.539D-02,-1.015D-01,-2.104D-01,-5.608D-01,-1.272D+00,-1.973D+00,-2.818D+00,-4.055D+00, &
      & -1.451D-01,-1.612D-01,-2.708D-01,-6.273D-01,-1.355D+00,-2.057D+00,-2.902D+00,-4.122D+00, &
      & -4.559D-01,-4.666D-01,-5.432D-01,-8.665D-01,-1.602D+00,-2.303D+00,-3.146D+00,-4.421D+00, &
      & -1.303D+00,-1.312D+00,-1.367D+00,-1.676D+00,-2.305D+00,-3.034D+00,-3.758D+00,-5.077D+00, &
      & -3.883D+00,-3.888D+00,-3.936D+00,-4.197D+00,-4.739D+00,-5.165D+00,-5.441D+00,-6.446D+00/), (/8,6/))
    REAL(KIND=DP), DIMENSION(8,6) :: SCO_DERIV
CONTAINS

!photodissociation rate of H2  accounting for self-shielding
REAL(dp) FUNCTION H2PhotoDissRate(NH2,radField,av,turbVel)
    !H2 Column density to surface, UV at that surface, visual extinction to surface, and turbulent velocity of cloud
    REAL(dp), INTENT(IN) :: NH2 ,radField,av ,turbVel
    !unshielded reaction rate, characteristic wavelendth of radiation, radiative linewidth
    REAL(dp), PARAMETER :: baseRate=5.18d-11,xl=1000.0,radWidth=8.0d7 
    REAL(dp) :: dopplerWidth

    dopplerWidth=turbVel/(xl*1.0d-8)
    H2PhotoDissRate = baseRate * (radField/1.7) * scatter(xl,av) * H2SelfShielding(NH2,dopplerWidth,radWidth)
END FUNCTION H2PhotoDissRate 

REAL(dp) FUNCTION COPhotoDissRate(NH2,NCO,radField,av)
    REAL(dp), INTENT(IN) :: NH2,NCO,radField,av!column densities of H2 and CO
    REAL(dp) :: ssf,lba,sca
    !calculate photodissociation rates for co (species # nco; reaction
    !# nrco) according to van dishoeck and black (apj 334, p771 (1988))
    !cocol is the co column density (in cm-2); the scaling of the pdrate
    !by a factor of 1.8 is due to the change from draine's is uv field
    !to the habing field
    ssf = COSelfShielding(NH2,NCO)
    lba = lbar(NCO,NH2)
    sca = scatter(lba,av)

    !The reason why rad is divided by 1.7 is that the alphas are for Draine and the rad is in 
    !Habing units
    COPhotoDissRate = (2.d-10) * (radfield/1.7) * ssf * sca
END FUNCTION COPhotoDissRate

FUNCTION cIonizationRate(alpha,gamma,gasTemp,NC,NH2,av,radfield) RESULT(RATE)
   REAL(DP) :: RATE
   REAL(DP), INTENT(IN) :: alpha,gamma,gasTemp,NC,NH2,av,radField
   REAL(DP) :: TAUC

!  Calculate the optical depth in the CI absorption band, accounting
!  for grain extinction and shielding by CI and overlapping H2 lines
   TAUC=gamma*av+1.1D-17*NC+(0.9D0*gasTemp**0.27D0*(NH2/1.59D21)**0.45D0)
!  Calculate the CI photoionization rate
   RATE=alpha*(radfield/1.7)*EXP(-TAUC)
   RETURN
END FUNCTION cIonizationRate
 
!-----------------------------------------------------------------------
!  H2 line self-shielding, adopting the treatment of
!  Federman, Glassgold & Kwan (1979, ApJ, 227, 466)
!-----------------------------------------------------------------------        
REAL(dp) FUNCTION H2SelfShielding(NH2,dopplerWidth,radWidth)
    REAL(dp), INTENT(IN) :: NH2,dopplerWidth,radWidth
    REAL(dp) ::  r, sj, sr, t, u, taud
    REAL(dp), PARAMETER :: FPARA=0.5,FOSC  = 1.0d-2
    !--------------------------------------------------------------
    !taud = opt. depth at line centre (assum. ortho:parah2=1)
    !pi**0.5 * e2 / (m(electr) * c) = 1.5e-2 cm2/s

    taud  = FPARA * NH2 * 1.5e-2 * FOSC / dopplerWidth
     
    !calculate doppler contribution of self shielding function sj
    IF (taud .eq. 0.0d0) THEN
       sj = 1.0d0
    ELSE IF (taud .lt. 2.0d0) THEN
       sj = exp(-0.6666667d0*taud)
    ELSE IF (taud .lt. 10.0d0) THEN
       sj = 0.638d0*taud**(-1.25d0)
    ELSE IF (taud .lt. 100.0d0) THEN
       sj = 0.505d0*taud**(-1.15d0)
    ELSE
       sj = 0.344d0*taud**(-1.0667d0)
    END IF

    !calculate wing contribution of self shielding function sr
    !IF (taud.lt.0.0d0)  taud=0.0d0
    IF (radWidth .eq. 0.0d0) then
       sr = 0.0d0
    ELSE
       r  = radWidth/(1.7724539d0*dopplerWidth)
       t  = 3.02d0 * ((r*1.0d+03)**(-0.064d0))
       u  = SQRT(taud*r)/t
       sr = r/(t*SQRT(0.78539816d0+u**2.0))
    ENDIF
    
    !calculate total self shielding function fgk
    H2SelfShielding = sj + sr
END FUNCTION H2SelfShielding
     


!calculate the influence of dust extinction (g=0.8, omega=0.3)
!wagenblast&hartquist, mnras237, 1019 (1989)     
REAL(dp) FUNCTION scatter(x1,av)


!---------------------------------------------------------------------
!         i/o variables type declaration
!       scat   : factor describing the influence of grain scattering
!                on the uv flux dependent on the total h number
!                density and wavelength of the uv radiation
!       x1      : wavelength (in angstrom)
!       cdntot : total h number density (in cm-2)
!   
!         program variables
!       av     : visual extinction in magnitudes (cf. savage et al.,
!                 1977 apj 216, p.291)
!        expo   : exponent
!        i      : loop index
!        tl     : tau(lambda)
!        tv     : tau(visual=5500a)
!        xlambda : function which calculates tl/tv
!        c(0)   : c(0) * exp(-k(0)*tau) : (rel.) intensity
!                 decrease for 0<=tau<=1 caused by grain
!                 scattering with g=0.8, omega=0.3
!                 (approximation)
!        c(i)   : sum (c(i) * exp(-k(i)*tau)) i=1,5  (rel.)
!                 intensity decrease for 1<=tau<=oo caused by
!                 grain scattering with g=0.8, omega=0.3.
!                 (cf. flannery, roberge, and rybicki 1980,
!                 apj 236, p.598).
!        k(0)   : see c0
!        k(i)   : see ci
!---------------------------------------------------------------------

!   i/o variables type declaration
    REAL(DP), INTENT(IN) :: x1,av
    !   
    !program variables type declaration
    REAL(dp), dimension(6) :: c=(/1.0d0,2.006d0,-1.438d0,7.364d-01,-5.076d-01,-5.920d-02/)
    REAL(dp), dimension(6) ::  k1=(/7.514d-01,8.490d-01,1.013d0,1.282d0,2.005d0,5.832d0/)
    REAL(dp)  :: expo, tl, tv
    INTEGER :: i

    !optical depth in the visual
    tv = av/ 1.086d0
      
    !make correction for the wavelength considered
    tl = tv * xlambda(x1)
       
    !calculate attuentuation  due to dust scattering
    scatter = 0.0d0
    IF (tl.lt.1.0d0) THEN
        expo = k1(1)*tl
        IF (expo.lt.100.0d0) THEN
            scatter = c(1) * EXP(-expo)
        ENDIF
    ELSE
        DO i=2,6
            expo = k1(i)*tl
            IF (expo.lt.100.0d0) THEN
            scatter = scatter + c(i)*EXP(-expo)
            ENDIF
        END DO
    ENDIF

END FUNCTION scatter

!=======================================================================
!
!  Determine the ratio of the optical depth at a given wavelength to
!  that at visual wavelength (λ=5500Å) using the extinction curve of
!  Savage & Mathis (1979, ARA&A, 17, 73, Table 2)
!
!-----------------------------------------------------------------------
!
!  Input parameters:
!  LAMBDA  = wavelength (in Å)
!
!  Program variables:
!  XLAMBDA       = Ratio of tau(λ)/tau(V) at the desired wavelength
!                  by 1D spline interpolation over the grid values
!  XLAMBDA_GRID  = tau(λ)/tau(V) values, determined by dividing the
!                  Aλ/E(B-V) values from Savage & Mathis (1979) by
!                  an assumed reddening of R=AV/E(B-V)=3.1
!  XLAMBDA_DERIV = 2nd derivative of XLAMBDA_GRID values from SPLINE
!  LAMBDA_GRID   = wavelengths (in Å) listed in Table 2
!  NUM_LAMBDA    = number of wavelengths
!  START         = .TRUE. when XLAMBDA is first called
!
!  Functions called:
!  SPLINE = second derivative of the supplied 1D function (in spline.f90)
!  SPLINT = 1-dimensional cubic spline interpolated value (in spline.f90)
!
!-----------------------------------------------------------------------
REAL(dp) FUNCTION xlambda(LAMBDA)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: LAMBDA

    REAL(dp) :: LAMBDA_VALUE

    IF (START) THEN
      CALL SPLINE(LAMBDA_GRID,XLAMBDA_GRID,NUM_LAMBDA,1.0D30,1.0D30,XLAMBDA_DERIV)
    ENDIF

    LAMBDA_VALUE=LAMBDA
    IF(LAMBDA.LT.LAMBDA_GRID(1)) LAMBDA_VALUE=LAMBDA_GRID(1)
    IF(LAMBDA.GT.LAMBDA_GRID(NUM_LAMBDA)) LAMBDA_VALUE=LAMBDA_GRID(NUM_LAMBDA)

    CALL SPLINT(LAMBDA_GRID,XLAMBDA_GRID,XLAMBDA_DERIV,NUM_LAMBDA,LAMBDA_VALUE,xlambda)
    IF(XLAMBDA.LT.0.0D0) XLAMBDA=0.0D0

    RETURN
END FUNCTION xlambda

!-----------------------------------------------------------------------
!self-shielding of CO due to 12CO self-shielding and H2 screening
!Use Van Dishoeck & Black APJ 334, 1988 for value 
!-----------------------------------------------------------------------


REAL(dp) FUNCTION COSelfShielding(NH2,NCO)
    REAL(dp), INTENT(IN) :: NH2, NCO
    REAL(dp) :: lognco,lognh2

    if (startr)  THEN
        call splie2(NCO_GRID,NH2_GRID,SCO_GRID,DIMCO,DIMH2,SCO_DERIV)
        startr = .false.
    END IF
   
    lognco = dlog10(NCO+1.0)
    lognh2 = dlog10(NH2+1.0)
       
    if (lognco.lt.NCO_GRID(1))      lognco = NCO_GRID(1)
    if (lognh2.lt.NH2_GRID(1))      lognh2 = NH2_GRID(1)
    if (lognco.gt.NCO_GRID(DIMCO))  lognco = NCO_GRID(DIMCO)
    if (lognh2.gt.NH2_GRID(DIMH2))  lognh2 = NH2_GRID(DIMH2)
       
    call splin2(NCO_GRID,NH2_GRID,SCO_GRID,SCO_DERIV,DIMCO,DIMH2,lognco,&
    &               lognh2,COSelfShielding)
    COSelfShielding = 10.0d0**COSelfShielding
END FUNCTION COSelfShielding
   
   
   
REAL(dp) FUNCTION lbar(u,w)
!calculate lambda bar (in a) according to equ. 4 of van dishoeck
!and black, apj 334, p771 (1988)
! --------------------------------------------------------------
!       i/o parameter
!       u : co column density in (cm-2)
!       w : h2 column density in (cm-2)
   
!        program variables
!        lu : log10(co column density in cm-2)
!        lw : log10(h2 column density in cm-2)
   
!--------------------------------------------------------------
    !i/o parameter type declaration
    REAL(dp)  u, w, lu, lw

    lu = dlog10(dabs(u)+1.0d0)
    lw = dlog10(dabs(w)+1.0d0)
    
    lbar = (5675.0d0 - 200.6d0*lw) - (571.6d0 - 24.09d0*lw)*lu +&
    &(18.22d0 - 0.7664d0*lw)*lu**2
       
    !lbar represents the mean of the wavelengths of the 33
    !dissociating bands weighted by their fractional contribution
    !to the total rate of each depth. lbar cannot be larger than
    !the wavelength of band 33 (1076.1a) and not be smaller than
    !the wavelength of band 1 (913.6a).
    if (lbar.gt.1076.1d0)  lbar = 1076.1d0
    if (lbar.lt. 913.6d0)  lbar =  913.6d0
END FUNCTION lbar

SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
!given an m by n tabulated function ya, and tabulated indepen-
!dent variables x1a (m values) and x2a (n values), this routine
!constructs one-dimensional natural cubic splines of the rows
!of ya and returns the second-derivatives in the array y2a.
!(copied from numerical recipes)

!--------------------------------------------------------------
!i/o parameter and program variables
          integer           nn
          parameter         (nn=100)
          integer           m, n, j, k
          REAL(dp)  x1a(m), x2a(n), ya(m,n), y2a(m,n), ytmp(nn),&
     &                      y2tmp(nn)
!--------------------------------------------------------------
    DO j=1,m
        DO k=1,n
            ytmp(k) = ya(j,k)
        END DO
        !values 1.0d30 signal a natural spline.
        call spline(x2a,ytmp,n,1.0d30,1.0d30,y2tmp)
        DO k=1,n
            y2a(j,k) = y2tmp(k)
        END DO
    END DO
return
!==============================================================
END SUBROUTINE splie2
  
SUBROUTINE spline(x,y,n,yp1,ypn,y2)

!calculate cubic spline for a set of points (x,y)
 
!(cf. "numerical recipes" 3.3 : routine spline)
!given arrays x and y of length n containing a tabulated
!function, i.e. y(i) = f(x(i)), with x(1) < x(2) < ... < x(n),
!and given values yp1 and ypn for the first derivative of the
!interpolating function at points 1 and n, respectively, this
!routine returns an array y2 of length n which contains the
!second derivatives of the interpolating function at the
!tabulated points x(i). if yp1 and/or ypn are equal to 1.0e+30
!or larger, the routine is signalled to set the corresponding
!boundary condition for a natural spline, with zero second
!derivative on that boundary.

!--------------------------------------------------------------
!i/o parameter
!x   : vector for independent variable x; dimension x(n)
!y   : vector for x-dependent variable y; dimension y(n)
!n   : dimension of vectors containing the tabulated function
!yp1 : 1. derivative of the interpolating function at point 1
!ypn : 1. derivative of the interpolating function at point n
!y2  : 2. derivative of the interpolating function
!--------------------------------------------------------------
   
!i/o parameter type declaration
    integer           n
    REAL(dp)  x(n), y(n), yp1, ypn, y2(n)
   
!program variables type declaration
    integer           i, k
    REAL(dp)  p, qn, sig, u(100), un
!--------------------------------------------------------------

    IF (yp1 .ge. 1.0d30) THEN
    !the lower boundary condition is set either to be
    !"natural"
        y2(1) =  0.0d0
        u(1)  =  0.0d0
    ELSE
    !or else to have a specified first derivative.
        y2(1) = -0.5d0
        u(1)  = (3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    ENDIF
   
    !this is the decomposition loop of the tridiagonal algorithm.
    !y2 and u are used for temporary storage of decomposed factors.
    DO  i=2,n-1
        sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p     = sig*y2(i-1) + 2.0d0
        y2(i) = (sig-1.0d0)/p
        u(i)  = (6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))&
        &                /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    END DO
       
    IF (ypn .ge. 1.0d30) THEN
    !     the upper boundary condition is set either to be
    !     "natural"
        qn = 0.0d0
        un = 0.0d0
    ELSE
    !     or else to have a specified first derivative.
        qn = 0.5d0
        un = (3.0d0/(x(n)-x(n-1))) *&
    &              (ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    ENDIF
       
    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
       
    !this is the backsubstitution loop of the tridiagonal algorithm
    DO k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
    END DO
END SUBROUTINE SPLINE
 
SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
!given x1a, x2a, ya, m, n as described in splie2 and y2a as
!produced by that routine; and given a desired interpolating
!point x1, x2; this routine returns an interpolated function
!value y by bicubic spline interpolation.
  
!--------------------------------------------------------------
!i/o parameter and program variables type declaration
integer           nn
parameter         (nn=100)
integer           m, n, j, k
REAL(dp)  x1a(m), x2a(n), ya(m,n), y2a(m,n), ytmp(nn),&
&                      y2tmp(nn), yytmp(nn), x1, x2, y
!--------------------------------------------------------------

! perform m evaluations of the row splines constructed by splie2
!using the one-dimensional spline evaluator splint.
    DO j=1,m
        DO k=1,n
            ytmp(k)  = ya(j,k)
            y2tmp(k) = y2a(j,k)
        END DO
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
    END DO
!construct the one-dimensional column spline and evaluate it.
    call spline(x1a,yytmp,m,1.0d30,1.0d30,y2tmp)
    call splint(x1a,yytmp,y2tmp,m,x1,y)
END SUBROUTINE splin2   
     
SUBROUTINE splint(xa,ya,y2a,n,x,y)
SAVE
!cubic spline interpolation
   
!(cf. "numerical recipes" 3.3 : routine splint, and 3.4.
!routine hunt)
!given the arrays xa and ya of length n, which tabulate a
!function (with the xa(i)'s in order), and given the array y2a,
!which is the output of routine cubspl, and given a value x,
!this routine returns a cubic-spline interpolated value y.
   
!--------------------------------------------------------------
!-i/o parameters
!-xa  : vector for independent variable x; dimension xa(n)
!-ya  : vector for x-dependent variable y; dimension ya(n)
!-y2a : 2. derivative of the interpolating function; dim. y2a(n)
!-n   : dimension of input vectors
!-x   : x value for which y is to be interpolated
!-y   : result of interpolation
!--------------------------------------------------------------
   
!--------------------------------------------------------------
    !i/o parameter type declaration
    integer           n,nstore
    REAL(dp)  x, xa(n), y, ya(n), y2a(n)
       
    !program variables type declaration
    integer           inc, jhi, jlo, jm
    REAL(dp)  h, a, b
    logical           ascnd
   
    !find interval xa(jlo) <= x <= xa(jlo+1) = xa(jhi)
    !ascnd is true if ascending order of table, false otherwise
    ascnd = xa(n).gt.xa(1)
    if (jlo.le.0 .or. jlo.gt.n) then
    !input guess not useful. go immediately to bisection.
        jlo = 0
        jhi = n+1
    ELSE           
        !set the hunting increment
        inc = 1
        
        !hunt up if ascending array or down if descending.
        IF (x.ge.xa(jlo) .eqv. ascnd) THEN
            !hunt up:
            jhi=jlo+inc
            IF (jhi .gt. n) THEN
                !done hunting since off end of table
                jhi=n+1
            ELSE
                !nstore is a work around for old 'go to' logic, if jhi exceeds n, that is fine
                !but the do while loop will break so jhi equals n temporarily and nstore holds
                !real value until we exit loop.
                nstore=1
                DO WHILE (((x.ge.xa(jhi)) .eqv. ascnd) .and. (jhi .lt. n))
                    !not done hunting
                    jlo=jhi
                    !so double increment
                    inc=inc+inc
                    !try again
                    jhi=jlo+inc
                    IF (jhi .gt. n) THEN
                        jhi=n
                        nstore=n+1
                    ENDIF
                END DO
                IF (nstore .eq. n+1) jhi=nstore
            !done hunting, value bracketed.
            END IF      
        ELSE
            jhi = jlo
            !hunt down:
            jlo = jhi-inc
            IF (jlo .lt. 1) THEN
            jlo=0
            ELSE
                nstore=1
                DO WHILE (((x.lt.xa(jlo)) .eqv. ascnd) .and. (jlo .gt. 1))
                    !not done hunting,
                    jhi = jlo
                    !so double the increment
                    inc = inc+inc
                    !and try again.
                    jlo = jhi-inc
                    IF (jlo .lt. 1) THEN
                        jlo=1
                        nstore=0
                    END IF
                END DO
                IF (nstore .eq. 0) jlo=nstore
            END IF
            !done hunting, since off end of table.
        ENDIF
    END IF   
    DO WHILE (jhi-jlo.ne.1)  
    !hunt is done, so begin final bisection phase:
        jm = (jhi+jlo)/2
        IF (x.gt.xa(jm) .eqv. ascnd) THEN
           jlo = jm
        ELSE
           jhi = jm
        ENDIF    
    END DO
       
    IF (jlo.eq.0)  THEN
        jlo = 1
        jhi = 2
    ENDIF
       
    !jlo and jhi now bracket the input value of x.
    !cubic spline polynomial is now evaluated.
    IF (jlo.eq.n)  THEN
        jlo = n-1
        jhi = n
    ENDIF
    h = xa(jhi) - xa(jlo)
    a = (xa(jhi) - x) / h
    b = (x - xa(jlo)) / h
    y = a*ya(jlo) + b*ya(jhi) +&
    &  ((a**3-a)*y2a(jlo) + (b**3-b)*y2a(jhi)) * (h**2)/6.0d0
END SUBROUTINE splint

END MODULE photoreactions


