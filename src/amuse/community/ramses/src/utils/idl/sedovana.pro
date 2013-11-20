;+
; NAME:
;       SEDOVANA
;
; PURPOSE:
;       This procedure computes the analytical solution for the Sedov
;       blast wave in 1, 2 or 3 dimensions (planar, cylindrical and
;       spherical case).
;
; CATEGORY:
;       Mathematical support routines.
;
; CALLING SEQUENCE:
;       SEDOVANA,r,d,u,p,gamma=gamma,dim=dim
;
; OPTIONAL INPUTS:
;       gamma:     if set, the adiabatic exponent of the fluid.
;       Default: 1.4
;
;       dim:       if set, the dimensionality of the problem.
;       Default: 1.
;
; OUTPUTS:
;       The routine returns 4 arrays of unknown size (determined
;       automatically by the routine). The output variables are
;       the following dimensionless quantities: r (position from the
;       point like explosion, d (density), u (velocity) and p
;       (pressure). To recover the true value, you have to rescale
;       these dimensionless values to the true values, defining first
;       the total energy E_0, the initial mass density rho_0 and the
;       time t you consider and finally computing the true values
;       using  the following scaling laws:
;
;       r = r * (E_0/rho_0)^(1./(dim+2.)) * t^(2./(dim+2.))
;       d = d * rho_0
;       u = u * (E_0/rho_0)^(1./(dim+2.)) * t^(-dim/(dim+2.))
;       p = p * (E_0/rho_0)^(2./(dim+2.)) * t^(-2.*dim/(dim+2.)) * rho_0
;
; EXAMPLE:
;       To compute the analytical solution of the planar Sedov blast
;       wave for a gamma=1.6667 fluid, type:
;
;               SEDOVANA,r,d,u,p,gamma=1.6667,dim=1
;
; MODIFICATION HISTORY:
;       Written by:     Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;       Fevrier, 2001:  Comments and header added by Romain Teyssier.
;-
;###################################################
;###################################################
;###################################################
pro sedovana,r,d,u,p,gamma=gamma,dim=dim

IF N_PARAMS() NE 4 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'sedovana'
    RETURN
ENDIF

if not keyword_set(gamma) then gamma=1.4d0
if not keyword_set(dim) then dim=1.0d0

n=long(dim)

print,'gamma=',gamma
g=double(gamma)
n=double(n)
n1=1000
n2=1000

vmax=4.d0/(n+2.)/(g+1.)
vmin=2.d0/(n+2.)/g
;v=vmin+(vmax-vmin)*(DINDGEN(n1)+1.0d0)/double(n1)
v=vmin+10.d0^(-10.d0*(1.0d0-(DINDGEN(n1)+1)/double(n1)))*(vmax-vmin)
a2 = (1.-g)/(2.*(g-1.)+n)
a1 = (n+2.)*g/(2.+n*(g-1.)) * ( 2.*n*(2.-g)/g/(n+2.)^2 - a2 )
a3 = n/(2.*(g-1)+n)
a4 = a1*(n+2.)/(2.-g)
a5 = 2./(g-2.)
a6 = g/(2.*(g-1.)+n)
a7 = a1*(2.+n*(g-1.))/(n*(2.-g))

r1 = ((n+2.)*(g+1.)/4.*v)^(-2./(2.+n))* $
    ((g+1.)/(g-1.)*( (n+2.)*g/2.*v-1.) )^(-a2)* $
    ( (n+2.)*(g+1.) / ( (n+2)*(g+1)-2.*(2.+n*(g-1.)) ) * (1.-(2.+n*(g-1.))/2.*v) )^(-a1)

u1 = (n+2.)*(g+1.)/4.*v*r1

d1 = ((g+1.)/(g-1.)*((n+2.)*g/2.*v-1.))^(a3)* $
    ((g+1.)/(g-1.)*(1.-(n+2.)/2.*v  ))^(a5)* $
    ((n+2.)*(g+1.)/( (n+2)*(g+1)-2.*(2.+n*(g-1.))) *(1.-(2.+n*(g-1.))/2.*v) )^(a4)    
     
p1 = ((n+2.)*(g+1.)/4.*v)^(2.*n/(2.+n))* $
    ((g+1.)/(g-1.)*(1.-(n+2.)/2.*v  ))^(a5+1.)* $
    ((n+2.)*(g+1.)/( (n+2)*(g+1)-2.*(2.+n*(g-1.))) *(1.-(2.+n*(g-1.))/2.*v) )^(a4-2.*a1)    
   
r2=r1(0)*(FINDGEN(n2)+0.5d0)/double(n2) ;10.d0^(-4.d0+4.d0*FINDGEN(100)/100.d0)
u2=u1(0)*r2/r1(0)
d2=d1(0)*(r2/r1(0))^(n/(g-1.0d0))
p2=p1(0)*(r2/r2)

r=fltarr(n_elements(r1)+n_elements(r2)+2L)
r(0:n_elements(r2)-1)=r2
r(n_elements(r2):n_elements(r1)+n_elements(r2)-1)=r1
r(n_elements(r1)+n_elements(r2))=max(r1)
r(n_elements(r1)+n_elements(r2)+1)=max(r1)+1000.
d=r
d(0:n_elements(r2)-1)=d2
d(n_elements(r2):n_elements(r1)+n_elements(r2)-1)=d1
d(n_elements(r1)+n_elements(r2))=1./((g+1.)/(g-1.))
d(n_elements(r1)+n_elements(r2)+1)=1./((g+1.)/(g-1.))
u=r
u(0:n_elements(r2)-1)=u2
u(n_elements(r2):n_elements(r1)+n_elements(r2)-1)=u1
u(n_elements(r1)+n_elements(r2))=0.
u(n_elements(r1)+n_elements(r2)+1)=0.
p=r
p(0:n_elements(r2)-1)=p2
p(n_elements(r2):n_elements(r1)+n_elements(r2)-1)=p1
p(n_elements(r1)+n_elements(r2))=0.
p(n_elements(r1)+n_elements(r2)+1)=0.

d=d*(g+1.)/(g-1.)
u=u*4./(n+2.)/(g+1.)
p=p*8./(n+2.)^2./(g+1.)

nn=n_elements(r)
vol=r
for i=1,nn-1 do vol(i)=r(i)^n-r(i-1)^n
vol(0)=r(0)^n
const=1.
if n eq 1. then const=2.0d0
if n eq 2. then const=!DPI
if n eq 3. then const=4.*!DPI/3.

vol=vol*const
int1=(d*u*u/2.d0)*vol
int2=p/(g-1.d0)*vol
sum1=total(int1,/double)
sum2=total(int2,/double)
sum=sum1+sum2
print,'chi0=',sum^(-1./(2.+n))
chi0=sum^(-1./(2.+n))
r=r*chi0
u=u*chi0
p=p*chi0^2
end


