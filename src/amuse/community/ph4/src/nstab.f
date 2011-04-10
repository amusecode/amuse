c	general three-body stability algorithm
c
c	system is unstable if nstab=1 is returned
c	system is stable if nstab=0 is returned
c
c	Rosemary Mardling
c	School of Mathematical Sciences, Monash University
c
c	version as of 16-11-07
c	email mardling@sci.monash.edu.au to be added to updates email list
c	preprint on astro-ph by New Year :-)
c
c	sigma=period ratio (outer/inner) **should be > 1**
c	ei0=initial inner eccentricity
c	eo=outer eccentricity
c	relinc=relative inclination (radians)
c	m1, m2, m3=masses (any units; m3=outer body)
c
c       valid for all inclinations
c
c	MASS RATIO CONDITIONS

c	valid for systems with at least one of m_2/m_1>0.05 OR
c	m_3/m_1>0.05 (so that one could have, for example, m_2/m_1=0 and
c	m_3/m_1=0.1) OR BOTH m2/m1>0.01 AND m3/m1>0.01

c	**future version will include other resonances to cover smaller
c	mass ratios
c
c	assumes resonance angle phi=0 because resonance overlap
c	criterion doesn't recognize instability outside separatrix.
c
	integer function nstab(sigma,ei0,eo,relinc,m1,m2,m3)
	
	implicit real*8 (a-h,m,o-z)
	common/params/mm1,mm2,mm3

	pi=4.d0*datan(1.0d0)
	
	mm1=m1
	mm2=m2
	mm3=m3
	
	m12=m1+m2
	m123=m12+m3
	
	Mi2=m3/m123
	Mo2=(m1*m2/m12**2)*(m12/m123)**(2./3.)
	Mi3=(m3/m12)*(m12/m123)**(4./3.)*(m1-m2)/m12
	Mo3=(m1*m2/m12**2)*(m12/m123)*(m1-m2)/m12
	
	c22=3./8.
	c20=0.25
	c31=sqrt(3.)/4.
	c33=-sqrt(5.)/4.
	
	e=eo
	
c	inclination coefficients

	win=0
	   
	A=sqrt(1-ei0**2)*cos(relinc)
	Z=(1-ei0**2)*(1+sin(relinc)**2)
     +		+ 5*ei0**2*(sin(win)*sin(relinc))**2
	Del=z**2+25+16*A**4-10*Z-20*A**2-8*A**2*Z
	   
	eK=sqrt(abs((Z+1-4*A**2+sqrt(Del))/6.))
	cosIK=A/sqrt(1-eK**2)
	sinIK=sqrt(1-cosIK**2)
	
	gam222=0.25*(1+cosIK)**2	   
	gam22m2=0.25*(1-cosIK)**2	   
	gam220=0.5*sqrt(1.5)*sinIK**2 
     	gam200=0.5*(3*cosIK**2-1)
	
c	induced inner eccentricity
	ei=ein_induced(sigma,ei0,e,relinc)
	
c	octopole emax	   	   
	if(m1.ne.m2)then
	   eoctmax=eoct(sigma,ei0,e)
	   ei=max(eoctmax,ei)	   
	endif
		   
	ei=max(eK,ei)
	ei=min(ei,1.)
	
	n=sigma
	nstab=0
		
c	[n:1](222) resonance
	s221=-3*ei+(13./8.)*ei**3+(5./192.)*ei**5         

	f22n=flmn(2,2,n,e)/(1-e)**3

	An=abs(6*c22*s221*f22n*(Mi2+Mo2*sigma**0.666)*gam222)
	phi=0
	En=0.5*(sigma-n)**2-An*(1+cos(phi))
	   
c	[n+1:1](222) resonance	   
	f22n=flmn(2,2,n+1,e)/(1-e)**3

	An=abs(6*c22*s221*f22n*(Mi2+Mo2*sigma**0.666)*gam222)
	
	Enp1=0.5*(sigma-(n+1))**2-An*(1+cos(phi))
	if(En.lt.0.and.Enp1.lt.0)nstab=1
	
c	[n:1](22-2) resonance
	s22m1=-(ei**3*(4480 + 1880*ei**2 + 1091*ei**4))/15360.        

	f22n=flmn(2,2,n,e)/(1-e)**3

	An=abs(6*c22*s22m1*f22n*(Mi2+Mo2*sigma**0.666)*gam22m2)

	phi=0
	En=0.5*(sigma-n)**2-An*(1+cos(phi))
	   
c	[n+1:1](22-2) resonance	   
	f22n=flmn(2,2,n+1,e)/(1-e)**3

	An=abs(6*c22*s22m1*f22n*(Mi2+Mo2*sigma**0.666)*gam22m2)
	
	Enp1=0.5*(sigma-(n+1))**2-An*(1+cos(phi))
	if(En.lt.0.and.Enp1.lt.0)nstab=1

c	[n:1](202) resonance
	s201=(ei*(-9216 + 1152*ei**2 - 48*ei**4 + ei**6))/9216.	  

	f22n=flmn(2,2,n,e)/(1-e)**3

	An=abs(6*sqrt(c20*c22)*s201*f22n*(Mi2+Mo2*sigma**0.666)*gam220)
		
	phi=0
	En=0.5*(sigma-n)**2-An*(1+cos(phi))
	   
c	[n+1:1](202) resonance	   
	f22n=flmn(2,2,n+1,e)/(1-e)**3

	An=abs(6*sqrt(c20*c22)*s201*f22n*(Mi2+Mo2*sigma**0.666)*gam220)
	
	Enp1=0.5*(sigma-(n+1))**2-An*(1+cos(phi))
	if(En.lt.0.and.Enp1.lt.0)nstab=1

c	[n:1](002) resonance
	s201=(ei*(-9216 + 1152*ei**2 - 48*ei**4 + ei**6))/9216.	  

	f20n=flmn(2,0,n,e)/(1-e)**3

	An=abs(3*c20*s201*f20n*(Mi2+Mo2*sigma**0.666)*gam200)
		
	phi=0
	En=0.5*(sigma-n)**2-An*(1+cos(phi))
	   
c	[n+1:1](002) resonance	   
	f20n=flmn(2,0,n+1,e)/(1-e)**3

	An=abs(3*c20*s201*f20n*(Mi2+Mo2*sigma**0.666)*gam200)
	
	Enp1=0.5*(sigma-(n+1))**2-An*(1+cos(phi))
	if(En.lt.0.and.Enp1.lt.0)nstab=1

8	continue
	
	end

c---------------------------------------------------------------------------	
c	Asymptotic expression for f^(lm)_n(e) for all e<1 and n.
c
      real*8 function flmn(l,m,n,e)

      implicit real*8 (a-h,o-z)

	if(e.lt.5.e-3)then
	   if(m.eq.n)then
	      flmn=1
	      else
	      flmn=0
	   endif
	   return
	endif
	
	pi=4.d0*datan(1.0d0)
      
      rho=n*(1-e)**1.5
		
      xi=(acosh(1/e)-sqrt(1-e**2))/(1-e)**1.5         
      
      flmn=(1/(2*pi*n))*2.0**m*(sqrt(2*pi)/facfac(l,m))*
     .        ((1+e)**(real(3*m-l-1)/4.)/e**m)*
     .        (rho**(real(l+m+1)/2.))*
     .        exp(-rho*xi)

      end


c---------------------------------------------------------------------------
	real*8 function ein_induced(sigma,ei0,e,relinc)

	implicit real*8 (a-h,m,o-z)	
	common/params/m1,m2,m3

	pi=4.d0*datan(1.0d0)
	
	m123=m1+m2+m3	
	n=sigma

     	gam222=0.25*(1+cos(relinc))**2
     	gam220=0.5*sqrt(1.5)*sin(relinc)**2
     	gam200=0.5*(3*cos(relinc)**2-1)
	
     	f22n=flmn(2,2,n,e)/(1-e)**3
     	f20n=flmn(2,0,n,e)/(1-e)**3
     	
     	prod222=f22n*gam222
     	prod220=f22n*gam220
     	prod200=f20n*gam200
     	
     	prod=max(prod222,prod220,prod200)
		    	
     	a=4.5*(m3/m123)*(2*pi*n)*prod/sigma**2
     	
     	ein_induced=sqrt(ei0**2+a**2)
     	     	
     	end

c------------------------------------------------------------------------------
c	eoct.f
c	
c	calculates maximum eccentricity for arbitrary coplanar system
c	using Mardling (2007) MNRAS in press
c
	real*8 function eoct(sigma,ei0,eo)
	implicit real*8 (a-h,m,o-z)
	parameter (ne=100)
	real*8 func0,funcpi
	real*4 ei(0:ne),f0(0:ne),fpi(0:ne)
	common/params/m1,m2,m3

	pi=4.d0*datan(1.0d0)
		
	m12=m1+m2
	m123=m12+m3
	aoai=((m123/m12)*sigma**2)**0.3333
	al=1/aoai
	
	epso=sqrt(1-eo**2)
	
	eeq=1.25*al*eo/epso**2/abs(1-sqrt(al)*(m2/m3)/epso)

	AA=abs(1-ei0/eeq)  
		
	if(AA.lt.1)then
	   eoct=(1+AA)*eeq
	else
	   eoct=ei0+2*eeq
	endif
		
	end

c---------------------------------------------------------------------------	
	real*8 function acosh(x)
	real*8 x

	acosh=dlog(x+dsqrt(x**2-1.d0))

	end	

c---------------------------------------------------------------------------
        real*8 function sgn(x)
        real*8 x

        if(x.lt.0)then
           sgn=-1
        else
           sgn=1
        endif

        end       	

c---------------------------------------------------------------------------	
      real*8 function facfac(l,m)
      implicit real*8 (a-h,o-z)

      prod=1

      n=l+m-1

      do i=1,n,2
         prod=prod*i
      enddo

      facfac=prod

      end
	

