! --------------------------------------------------------------------------
! get_neutrino_rate
! --------------------------------------------------------------------------
!  Get energy loss rate from neutrinos, in erg/g/s
!  Fits based on Itoh & al. (1996), based on a program by F. Timmes
!  available at http://cococubed.asu.edu/codes/nuloss/public_sneut4.f
! Input:
!  T     -  Temperature [K]
!  RHO   -  Density (g/cm^3)
!  NIO   -  Number of ions per baryon (from EoS)
!  NEO   -  Average charge per baryon (=number of electrons, from EoS)
! Output:
!  EN    -  Energy loss rate due to neutrinos, in erg/g/s
! --------------------------------------------------------------------------
      subroutine get_neutrino_rate(t, rho, nio, neo, en)
      implicit none
      double precision, intent(in) :: t, rho, nio, neo
      double precision, intent(out) :: en
      double precision :: abar, zbar, pair, plas, phot, brem, recomb

!     Compute average mass and charge per nucleus
      abar = 1.0d0 / nio
      zbar = neo * abar

!     Get neutrino rates
      call sneut4(t, rho, zbar, abar,
     &            pair, plas, phot, brem, recomb, en)
      return
      end subroutine

      subroutine sneut4(temp,den,zbar,abar,
     1                  pair,plas,phot,brem,recomb,total)
      implicit none 
      save  

c..this routine computes neutrino losses from the analytic fits of
c..itoh et al. apjs 102, 411, 1996

c..input:
c..temp = temperature (in K)
c..den  = density in gr/cm**3)
c..zbar = mean isotopic charge
c..abar = mean nucleon number 

c..output in erg/g/s:
c..pair   = pair neutrino contributions 
c..plas   = plasma neutrino contributions
c..phot   = photoneutrino contributions
c..brem   = bremstrahlung neutrino contributions
c..recomb = recombination neutrino contributions
c..sneut  = total neutrino loss rate 

c..declare
      double precision temp,den,zbar,abar,
     1                 pair,plas,phot,brem,recomb,total

      double precision xmue,t9,xl,xlp5,xl2,xl3,xl4,xl5,xl6,xl8,xl9,
     1                 xlm1,xlm2,xlm3,rm,xnum,xden,fpair,fphoto,gl,
     2                 zeta,zeta2,zeta3,qpair,gl2,gl12,gl32,gl72,gl6,
     3                 ft,fl,cc,xlnt,fxy,qv,tau,qphoto,den6,tfermi,t8,
     4                 t832,t83,t86,t8m2,t8m5,eta,etam1,etam2,fbrem,
     5                 gbrem,a0,a1,a2,a3,b1,b2,b3,c00,
     6                 c01,c02,c03,c04,c05,c06,c10,c11,c12,c13,c14,c15,
     7                 c16,c20,c21,c22,c23,c24,c25,c26,dd01,dd02,dd03,
     8                 dd04,dd05,dd11,dd12,dd13,dd14,dd15,dd21,dd22,
     9                 dd23,dd24,dd25,u,gamma,gm1,gm13,gm23,v,w,fb,
     &                 gb,gt,fliq,gliq,ifermi12,nu,nu2,nu3,b,c,d,f1,
     1                 f2,f3,z,bigj,cos1,cos2,cos3,cos4,cos5,sin1,sin2,
     2                 sin3,sin4,sin5,ye,last,deni

      double precision pi,fac1,fac2,fac3,third,twoth,con1,sixth
      parameter        (pi      = 3.1415926535897932384d0,
     1                  fac1    = 5.0d0 * pi / 3.0d0,
     2                  fac2    = 10.0d0 * pi,
     3                  fac3    = pi / 5.0d0,
     4                  third   = 1.0d0/3.0d0,
     5                  twoth   = 2.0d0/3.0d0,
     6                  con1    = 1.0d0/5.9302d0,
     7                  sixth   = 1.0d0/6.0d0)


c..theta is sin**2(theta_weinberg) = 0.2319 plus/minus 0.00005 (1996)
c..xnufam is the number of neutrino flavors = 3.02 plus/minus 0.005 (1998)
c..change theta and xnufam if need be, and the changes will automatically
c..propagate through the routine. cv and ca are the vector and axial currents.

      double precision theta,xnufam,cv,ca,cvp,cap,tfac1,tfac2,tfac3,
     1                 tfac4,tfac5,tfac6
      parameter        (theta  = 0.2319d0,
     1                  xnufam = 3.0d0,
     2                  cv     = 0.5d0 + 2.0d0 * theta,
     3                  cvp    = 1.0d0 - cv,
     4                  ca     = 0.5d0,
     5                  cap    = 1.0d0 - ca,
     6                  tfac1  = cv*cv + ca*ca + 
     7                           (xnufam-1.0d0) * (cvp*cvp+cap*cap),
     8                  tfac2  = cv*cv - ca*ca + 
     9                           (xnufam-1.0d0) * (cvp*cvp - cap*cap),
     &                  tfac3  = tfac2/tfac1,
     1                  tfac4  = 0.5d0 * tfac1,
     2                  tfac5  = 0.5d0 * tfac2,
     3                  tfac6  = cv*cv + 1.5d0*ca*ca + (xnufam - 1.0d0)*
     4                           (cvp*cvp + 1.5d0*cap*cap))


c..initialize and bail if its too cold
      pair   = 0.0d0
      plas   = 0.0d0
      phot   = 0.0d0
      brem   = 0.0d0
      recomb = 0.0d0
      total  = 0.0d0
      if (temp .lt. 1.0e7) return

c..some frequent factors
      xmue  = abar/zbar
      ye    = zbar/abar
      t9    = temp * 1.0d-9

      xl    = t9 * con1
      xlp5  = sqrt(xl)
      xl2   = xl*xl
      xl3   = xl2*xl
      xl4   = xl2*xl2
      xl5   = xl*xl4
      xl6   = xl2*xl4
      xl8   = xl4*xl4
      xl9   = xl8 * xl
      xlm1  = 1.0d0/xl
      xlm2  = xlm1*xlm1
      xlm3  = xlm1*xlm2

      rm    = den/xmue
      zeta  = (rm * 1.0d-9)**third * xlm1
      zeta2 = zeta * zeta
      zeta3 = zeta2 * zeta



c..pair neutrino section
c..for reactions like e+ + e- => nu_e + nubar_e 
c..equation 2.8 
      gl = 1.0d0 - 13.04d0*xl2 + 133.5d0*xl4 +1534.0d0*xl6 + 918.6d0*xl8

c..equation 2.7
      if (t9 .lt. 10.0) then
       xnum = (6.002d19 + 2.084d20*zeta + 1.872d21*zeta2)
     1         * exp(-5.5924d0*zeta)
       xden = zeta3 + 9.383d-1*xlm1 - 4.141d-1*xlm2 + 5.829d-2*xlm3
      else
       xnum = (6.002d19 + 2.084d20*zeta + 1.872d21*zeta2)
     1         * exp(-4.9924d0*zeta)
       xden = zeta3 + 1.2383d0*xlm1 - 8.141d-1*xlm2 
      end if
      fpair = xnum/xden

c..equation 2.6
      xnum   = 1.0d0 / (10.7480d0*xl2 + 0.3967d0*xlp5 + 1.005d0)
      xden   = (1.0d0 + rm/(7.692d7*xl3 + 9.715d6*xlp5) )**(-0.3d0)
      qpair  = xnum*xden

c..equation 2.5
      pair     = tfac4 * (1.0d0 + tfac3 * qpair) 
     1           * gl * exp(-2.0d0*xlm1) * fpair





c..plasma neutrino section 
c..for collective reactions like gamma_plasmon => nu_e + nubar_e
c..equation 4.6
      gl2  = 1.1095d11*rm/(temp*temp * sqrt(1.0d0+(1.019d-6*rm)**twoth))
      gl   = sqrt(gl2)
      gl12 = sqrt(gl)
      gl32 = gl * gl12
      gl72 = gl2 * gl32
      gl6  = gl2 * gl2 * gl2

c..equation 4.7
      ft   = 2.4d0 + 0.6d0*gl12 + 0.51d0*gl + 1.25d0*gl32

c..equation 4.8
      fl   = (8.6d0*gl2 + 1.35d0*gl72)/(225.0d0 - 17.0d0*gl + gl2)

c..equation 4.9 and 4.10
      cc   = log10(2.0d0*rm)
      xlnt = log10(temp)
      xnum = sixth * ( 17.5d0 + cc - 3.0d0 * xlnt)
      xden = sixth * (-24.5d0 + cc + 3.0d0 * xlnt)

c..equation 4.11
      if (abs(xnum) .gt. 0.7d0  .or.  xden .lt. 0.0d0) then
       fxy = 1.0d0
      else 
       fxy = 1.05d0 + (0.39d0 - 1.25d0*xnum - 0.35d0*sin(4.5d0*xnum)
     1                 - 0.3d0 * exp(-1.0d0*(4.5d0*xnum + 0.9d0)**2))
     2       * exp(-1.0d0* ( min(0.0d0, xden - 1.6d0 + 1.25d0*xnum)
     3                 / (0.57d0 - 0.25d0*xnum) )**2  )
      end if

c..equation 4.5
      qv = 3.0d21 * xl9 * gl6 * exp(-gl) * (ft + fl) * fxy

c..equation 4.1
      plas = 0.93153d0 * qv
      




c..photoneutrino process section 
c..for reactions like e- + gamma => e- + nu_e + nubar_e
c..                   e+ + gamma => e+ + nu_e + nubar_e
c..equation 3.8 for tau, equation 3.6 for cc,
c..and table 2 written out for speed
      if (temp .ge. 1.0d7  .and. temp .lt. 1.0d8) then
       tau  =  log10(temp * 1.0d-7)
       cc   =  0.5654d0 + tau
       c00  =  1.008d11
       c01  =  0.0d0
       c02  =  0.0d0
       c03  =  0.0d0
       c04  =  0.0d0
       c05  =  0.0d0
       c06  =  0.0d0
       c10  =  8.156d10
       c11  =  9.728d8
       c12  = -3.806d9
       c13  = -4.384d9
       c14  = -5.774d9
       c15  = -5.249d9
       c16  = -5.153d9
       c20  =  1.067d11
       c21  = -9.782d9 
       c22  = -7.193d9
       c23  = -6.936d9
       c24  = -6.893d9
       c25  = -7.041d9
       c26  = -7.193d9
       dd01 =  0.0d0
       dd02 =  0.0d0
       dd03 =  0.0d0
       dd04 =  0.0d0
       dd05 =  0.0d0
       dd11 = -1.879d10
       dd12 = -9.667d9
       dd13 = -5.602d9
       dd14 = -3.370d9
       dd15 = -1.825d9
       dd21 = -2.919d10
       dd22 = -1.185d10
       dd23 = -7.270d9
       dd24 = -4.222d9
       dd25 = -1.560d9

      else if (temp .ge. 1.0d8  .and. temp .lt. 1.0d9) then
       tau  =  log10(temp * 1.0d-8)
       cc   =  1.5654d0
       c00  =  9.889d10 
       c01  = -4.524d8
       c02  = -6.088d6 
       c03  =  4.269d7 
       c04  =  5.172d7 
       c05  =  4.910d7 
       c06  =  4.388d7
       c10  =  1.813d11
       c11  = -7.556d9 
       c12  = -3.304d9  
       c13  = -1.031d9
       c14  = -1.764d9  
       c15  = -1.851d9
       c16  = -1.928d9
       c20  =  9.750d10
       c21  =  3.484d10
       c22  =  5.199d9  
       c23  = -1.695d9  
       c24  = -2.865d9  
       c25  = -3.395d9  
       c26  = -3.418d9
       dd01 = -1.135d8   
       dd02 =  1.256d8   
       dd03 =  5.149d7   
       dd04 =  3.436d7   
       dd05 =  1.005d7
       dd11 =  1.652d9  
       dd12 = -3.119d9  
       dd13 = -1.839d9  
       dd14 = -1.458d9  
       dd15 = -8.956d8
       dd21 = -1.549d10  
       dd22 = -9.338d9  
       dd23 = -5.899d9  
       dd24 = -3.035d9  
       dd25 = -1.598d9

      else if (temp .ge. 1.0d9) then
       tau  =  log10(t9)
       cc   =  1.5654d0
       c00  =  9.581d10
       c01  =  4.107d8
       c02  =  2.305d8   
       c03  =  2.236d8   
       c04  =  1.580d8   
       c05  =  2.165d8   
       c06  =  1.721d8
       c10  =  1.459d12
       c11  =  1.314d11
       c12  = -1.169d11  
       c13  = -1.765d11  
       c14  = -1.867d11  
       c15  = -1.983d11  
       c16  = -1.896d11
       c20  =  2.424d11
       c21  = -3.669d9
       c22  = -8.691d9  
       c23  = -7.967d9  
       c24  = -7.932d9  
       c25  = -7.987d9  
       c26  = -8.333d9
       dd01 =  4.724d8
       dd02 =  2.976d8   
       dd03 =  2.242d8   
       dd04 =  7.937d7   
       dd05 =  4.859d7
       dd11 = -7.094d11
       dd12 = -3.697d11
       dd13 = -2.189d11  
       dd14 = -1.273d11  
       dd15 = -5.705d10
       dd21 = -2.254d10
       dd22 = -1.551d10
       dd23 = -7.793d9
       dd24 = -4.489d9
       dd25 = -2.185d9
      end if

c..equation 3.7, compute the expensive trig functions only one time
      cos1 = cos(fac1*tau)
      cos2 = cos(fac1*2.0d0*tau)
      cos3 = cos(fac1*3.0d0*tau)
      cos4 = cos(fac1*4.0d0*tau)
      cos5 = cos(fac1*5.0d0*tau)
      last = cos(fac2*tau)

      sin1 = sin(fac1*tau)
      sin2 = sin(fac1*2.0d0*tau)
      sin3 = sin(fac1*3.0d0*tau)
      sin4 = sin(fac1*4.0d0*tau)
      sin5 = sin(fac1*5.0d0*tau)

      a0 = 0.5d0*c00 
     1     + c01*cos1 + dd01*sin1 + c02*cos2 + dd02*sin2
     2     + c03*cos3 + dd03*sin3 + c04*cos4 + dd04*sin4
     3     + c05*cos5 + dd05*sin5 + 0.5d0*c06*last
      a1 = 0.5d0*c10 
     1     + c11*cos1 + dd11*sin1 + c12*cos2 + dd12*sin2
     2     + c13*cos3 + dd13*sin3 + c14*cos4 + dd14*sin4
     3     + c15*cos5 + dd15*sin5 + 0.5d0*c16*last
      a2 = 0.5d0*c20 
     1     + c21*cos1 + dd21*sin1 + c22*cos2 + dd22*sin2
     2     + c23*cos3 + dd23*sin3 + c24*cos4 + dd24*sin4
     3     + c25*cos5 + dd25*sin5 + 0.5d0*c26*last

c..equation 3.4
      xnum   = (a0 + a1*zeta + a2*zeta2)*exp(-cc*zeta)
      xden   = zeta3 + 6.290d-3*xlm1 + 7.483d-3*xlm2 + 3.061d-4*xlm3
      fphoto = xnum/xden

c..equation 3.3
      xnum   = 0.666d0*((1.0d0 + 2.045d0 * xl)**(-2.066d0))
      xden   = 1.0d0 + rm / (1.875d8*xl + 1.653d8*xl2 
     1                     + 8.449d8*xl3 - 1.604d8*xl4)
      qphoto = xnum/xden

c..equation 3.2
      phot   =  tfac4 * (1.0d0 - tfac3 * qphoto) * rm * xl5 * fphoto
      phot   = max(phot,0.0d0)





c..bremsstrahlung neutrino section 
c..for reactions like e- + (z,a) => e- + (z,a) + nu + nubar
c..                   n  + n     => n + n + nu + nubar
c..                   n  + p     => n + p + nu + nubar
c..equation 4.3
      den6   = den * 1.0d-6
      tfermi = 5.9302d9*(sqrt(1.0d0+1.018d0*(den6/xmue)**twoth)-1.0d0)

c.."weak" degenerate electrons only
      if (temp .gt. 0.3d0 * tfermi) then
       t8     = temp * 1.0d-8
       t832   = t8 * sqrt(t8)
       t83    = t8 * t8 * t8
       t86    = t83 * t83
       t8m2   = 1.0d0/(t8 * t8)
       t8m5   = t8m2 * t8m2/t8

c..equation 5.3
       eta   = rm/(7.05d6 * t832 + 5.12d4 * t83)
       etam1 = 1.0d0/eta
       etam2 = etam1 * etam1

c..equation 5.2
       xnum  = 1.0d0/(23.5d0 + 6.83d4*t8m2 + 7.81d8*t8m5)
       xden  = 1.26d0*(1.0d0+etam1)/(1.0d0+1.47d0*etam1+3.29d-2*etam2)
       fbrem = xnum + xden

c..equation 5.9
       xnum  = 1.0d0/( (1.0d0 + rm*1d-9)
     1                * (230.0d0 + 6.7d5*t8m2 + 7.66d9*t8m5) )
       c00   = 7.75d5*t832 + 247.0d0*t8**(3.85d0)
       c01   = 4.07d0 + 0.0240d0 * t8**(1.4d0)
       c02   = 4.59d-5 * t8**(-0.110d0)
       xden  = 1.0d0/(c00/rm + c01 + c02 * den**(0.656d0))
       gbrem = xnum + xden

c..equation 5.1
       brem  = 0.5738d0*zbar*ye*t86*den * (tfac4*fbrem - tfac5*gbrem)


c..liquid metal with c12 parameters (not too different for other elements)
c..equation 5.18 and 5.16
      else
       t8    = temp * 1.0d-8
       t86   = t8 * t8 * t8 * t8 * t8 * t8
       u     = fac3 * (log10(den) - 3.0d0)
       gamma = 2.275d-1 * zbar * zbar/t8 * (den6/abar)**third
       gm1   = 1.0d0/gamma
       gm13  = gm1**third
       gm23  = gm13 * gm13

c..equation 5.25 and 5.26
       v = -0.05483d0 - 0.01946d0*gm13 + 1.86310d0*gm23 - 0.78873d0*gm1
       w = -0.06711d0 + 0.06859d0*gm13 + 1.74360d0*gm23 - 0.74498d0*gm1

c..compute the expensive trig functions of equation 5.21 only once
       cos1 = cos(u)
       cos2 = cos(2.0d0*u)
       cos3 = cos(3.0d0*u)
       cos4 = cos(4.0d0*u)
       cos5 = cos(5.0d0*u)

       sin1 = sin(u)
       sin2 = sin(2.0d0*u)
       sin3 = sin(3.0d0*u)
       sin4 = sin(4.0d0*u)

c..equation 5.21
       fb =  0.5d0 * 0.17946d0  + 0.00945d0*u + 0.34529d0   
     1       - 0.05821d0*cos1 - 0.04969d0*sin1
     2       - 0.01089d0*cos2 - 0.01584d0*sin2
     3       - 0.01147d0*cos3 - 0.00504d0*sin3
     4       - 0.00656d0*cos4 - 0.00281d0*sin4
     5       - 0.00519d0*cos5 
      
c..equation 5.22
       ft =  0.5d0 * 0.06781d0 - 0.02342d0*u + 0.24819d0
     1       - 0.00944d0*cos1 - 0.02213d0*sin1
     2       - 0.01289d0*cos2 - 0.01136d0*sin2
     3       - 0.00589d0*cos3 - 0.00467d0*sin3
     4       - 0.00404d0*cos4 - 0.00131d0*sin4
     5       - 0.00330d0*cos5 

c..equation 5.23
       gb =  0.5d0 * 0.00766d0 - 0.01259d0*u + 0.07917d0
     1       - 0.00710d0*cos1 + 0.02300d0*sin1
     2       - 0.00028d0*cos2 - 0.01078d0*sin2
     3       + 0.00232d0*cos3 + 0.00118d0*sin3
     4       + 0.00044d0*cos4 - 0.00089d0*sin4
     5       + 0.00158d0*cos5

c..equation 5.24
       gt =  -0.5d0 * 0.00769d0  - 0.00829d0*u + 0.05211d0
     1       + 0.00356d0*cos1 + 0.01052d0*sin1
     2       - 0.00184d0*cos2 - 0.00354d0*sin2
     3       + 0.00146d0*cos3 - 0.00014d0*sin3
     4       + 0.00031d0*cos4 - 0.00018d0*sin4
     5       + 0.00069d0*cos5 

c..equation 5.19 and 5.20
       fliq = v*fb + (1.0d0 - v)*ft
       gliq = w*gb + (1.0d0 - w)*gt

c..equation 5.17
       brem = 0.5738d0*zbar*ye*t86*den * (tfac4*fliq - tfac5*gliq)
      end if




c..recombination neutrino section
c..for reactions like e- (continuum) => e- (bound) + nu_e + nubar_e
c..equation 6.11 solved for nu
      xnum = 1.10520d8 * den * ye /(temp*sqrt(temp))
      nu   = ifermi12(xnum)
      nu2  = nu * nu
      nu3  = nu2 * nu

c..equation 6.7 and table 12
      zeta = 1.579d5 * zbar * zbar/temp
      if (nu .ge. -20.0  .and. nu .lt. 0.0) then
       a1 = 1.51d-2
       a2 = 2.42d-1
       a3 = 1.21d0
       b  = 3.71d-2
       c  = 9.06e-1
       d  = 9.28d-1
       f1 = 0.0d0
       f2 = 0.0d0
       f3 = 0.0d0
      else if (nu .ge. 0.0  .and. nu .le. 10.0) then
       a1 = 1.23d-2
       a2 = 2.66d-1
       a3 = 1.30d0
       b  = 1.17d-1
       c  = 8.97e-1
       d  = 1.77d-1
       f1 = -1.20d-2
       f2 = 2.29d-2
       f3 = -1.04d-3
      end if

c..equation 6.13  and 6.14
      if (nu .ge. -20.0  .and.  nu .le. 10.0) then
       z    = zeta/(1.0d0 + f1*nu + f2*nu2 + f3*nu3)  
       bigj = (a1/z + a2*z**(-2.25) + a3*z**(-4.55)) * exp(nu)
     1        / (1.0d0 + b*exp(c*nu)*(1.0d0 + d*z))

c..equation 6.5
       recomb = tfac6 * 2.649d-18 * ye* zbar**13 * den * bigj
     1          / (exp(zeta + nu) + 1.0d0) 
      end if 




c..convert from erg/cm^3/s to erg/g/s 
c..comment these out to duplicate the itoh et al plots
      deni   = 1.0d0/den
      pair   = pair*deni
      plas   = plas*deni
      phot   = phot*deni
      brem   = brem*deni
      recomb = recomb*deni


c..the total neutrino loss rate
      total = plas + pair + phot + brem + recomb
      return
      end

      double precision function ifermi12(f)
      implicit none
      save

c..this routine applies a rational function expansion to get the inverse
c..fermi-dirac integral of order 1/2 when it is equal to f.
c..maximum error is 4.19d-9.   reference: antia apjs 84,101 1993

c..declare
      integer          i,m1,k1,m2,k2
      double precision f,an,a1(12),b1(12),a2(12),b2(12),rn,den,ff


c..load the coefficients of the expansion
      data  an,m1,k1,m2,k2 /0.5d0, 4, 3, 6, 5/
      data  (a1(i),i=1,5)/ 1.999266880833d4,   5.702479099336d3,
     1     6.610132843877d2,   3.818838129486d1,
     2     1.0d0/
      data  (b1(i),i=1,4)/ 1.771804140488d4,  -2.014785161019d3,
     1     9.130355392717d1,  -1.670718177489d0/
      data  (a2(i),i=1,7)/-1.277060388085d-2,  7.187946804945d-2, 
     1                    -4.262314235106d-1,  4.997559426872d-1,
     2                    -1.285579118012d0,  -3.930805454272d-1,
     3     1.0d0/
      data  (b2(i),i=1,6)/-9.745794806288d-3,  5.485432756838d-2,
     1                    -3.299466243260d-1,  4.077841975923d-1,
     2                    -1.145531476975d0,  -6.067091689181d-2/


      if (f .lt. 4.0d0) then
       rn = f + a1(m1)
       do i=m1-1,1,-1
        rn = rn*f + a1(i)
       enddo
       den = b1(k1+1)
       do i=k1,1,-1
        den = den*f + b1(i)
       enddo
       ifermi12 = log(f * rn/den)

      else
       ff = 1.0d0/f**(1.0d0/(1.0d0 + an))
       rn = ff + a2(m2)
       do i=m2-1,1,-1
        rn = rn*ff + a2(i)
       enddo
       den = b2(k2+1)
       do i=k2,1,-1
        den = den*ff + b2(i)
       enddo
       ifermi12 = rn/(den*ff)
      end if

      return
      end
