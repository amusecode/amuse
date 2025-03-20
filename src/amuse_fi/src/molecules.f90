module MoleculeMod
  implicit none
 private
 public :: InitMolecule, h2fraction,polyQ 

 real, parameter :: xifuv=3             ! 2-4
 real, parameter :: mu=17.5              ! formation rate const. mu, min: 3.5
 integer, parameter :: poly=1           ! polytropic index
 real, parameter :: nzero=1500          ! cm^-3
 real, parameter :: pzero=1.e4          ! K cm^-3
 real, parameter :: pcut=1000.		! K cm^-3 (transition pressure) 
 real, parameter :: k0=4.e-11		! dissociation rate s^-1
 real, parameter :: kboltz=1.38e-23     ! SI
 real, parameter :: amu=1.67e-27        ! SI
 real, parameter :: sqrtpi=1.772453851
 real, parameter :: year=3.15576e7
 real, parameter :: sigmav=4.9e-22	! cm^2
 real, parameter :: pc=3.086e18		! cm

 real, parameter :: Tform=1000,Tdes=3000,T4=10000	! temperature intervals 

 real, save :: metal,meanmwt,meann,fmin,n0con

 logical, parameter :: debug=.FALSE.

 contains

subroutine InitMolecule(metallicity)
  use ElementsMod
 character(len=*), optional :: metallicity

 
 print*,' using molecule formation routines!'
 print*,'  > formation rate constant:', mu
 if(poly.EQ.0) print*,'  > constant density clouds'
 if(poly.EQ.1) print*,'  > logotropic clouds' 
     
 if(present(metallicity)) then
  call InitElements(metallicity)
 else
  call InitElements()
 endif

 metal=zQ()/solarzQ()
 meanmwt=meanmwtQ() 
 meann=1./fractionQ(elementnrQ("H   "))
 n0con=xifuv*sigmav*pc
 
 if(poly.NE.0.AND.poly.NE.1) call molError(1)
 
 if(poly.eq.0) fmin=0.
 if(poly.eq.1) fmin=1.e-6
 
end subroutine 


! dt in year, n in cm^-3, T in K , g0 in units of Draine field, disp in km/s
subroutine h2fraction(fh2,dt,n,T,g0,disp,z,xe)
  real,intent(inout) :: fh2
  real,intent(in) :: n,T,g0,dt,disp
  real,intent(in), optional :: z,xe
  real ptot,lz,afuv,amean,ficonstant,dtau,snzero,ipar,Rf,acr,lxe
  integer :: i,j
 
 if(.not.present(z)) then
  lz=metal
 else
  lz=z
 endif
 
 if(.not.present(xe)) then
  lxe=0.0
 else
  lxe=xe
 endif  
 
 if(debug) print*,'*1*',fh2,dt,n,T,g0,disp,lz
 
 fh2=max(fmin,fh2)

 if(dt.LE.0) return

 if(n.GT.1.e10.AND.T.LE.Tform) then
  fh2=1.
  return
 endif    
 if(T.GT.T4) then
  fh2=fmin
  return
 endif 
 if(T.GT.Tdes) then 
  dtau=dt/taudiss(T,n,lxe)
  call fh2des(dtau,fh2) 
 endif
 if(debug) print*,'*2*',fh2,dtau
! ficonstant=6.6e-6*sqrtpi*(lz*xifuv)**.5
 ficonstant=3.39e-6*sqrtpi*(lz*xifuv)**.5
 ptot=(meann+lxe)*n*T+1.e6*meanmwt*meann*n*disp**2*amu/kboltz/3.
 amean=meanAv(ptot,lz)
 afuv=xifuv*atrans(fh2,amean)
 if(debug) print*,'*3*', ficonstant,ptot,amean,afuv
 if(T.GT.Tform) then
  dtau=2*g0*k0*ficonstant*dt*year
!  afuv=log(exp(afuv)+dtau) 
  afuv=afuv+log(1+dtau*exp(-afuv))
 endif
 if(debug) print*,'*4*', dtau,afuv
 if(T.LE.Tform) then
  Rf=rform(T,lz)
  dtau=2*n*Rf*dt*year
  ipar=g0*k0*ficonstant/n/Rf
  if(debug) print*,'*5*',Rf,dtau,ipar
  if(poly.EQ.0) then
   call solveh2_explicit(ipar,dtau,afuv)   
!   j=INT(10*dtau)+1
!   do i=1,j
!    call solveh2(ipar,dtau/j,afuv)
!   enddo 
   if(debug) print*,'*6*',j,afuv
  endif
  if(poly.EQ.1) then
   dtau=dtau*2./3.
   snzero=2./3.*nzero*n0con*(ptot/pzero)**.5
   call solveh2log_explicit(ipar,dtau,snzero,afuv)   
!   j=INT(10*dtau)+1
!   do i=1,j
!    call solveh2log(ipar,dtau/j,snzero,afuv)
!   enddo 
   if(debug) print*,'*6*',j,afuv,snzero
  endif
 endif
 acr=afuv/xifuv
 fh2=fmol(acr,amean)
 fh2=max(fmin,fh2)
end subroutine


subroutine molError(i)
  integer i
 
 if(i.eq.1) print*,' poly not equal to 0 or 1'
 
 stop
  
end subroutine


function rform(T,z) ! formation rate in cm^3 s^-1
 real :: rform,z,T
 real, parameter :: Tref=102
 rform=3.5e-18*z*T**.5/(1+T/Tref)**2*mu
 if(T.gt.Tform) rform=0
end function 

function fmol(acr,av) 
  real :: fmol,acr,av
 if(poly.EQ.0) then
  if(acr.gt.3./4.*av) then 
   fmol=0
  else
   fmol=(1.-4./3.*acr/av)**3
  endif
 endif
 if(poly.EQ.1) then
  fmol=exp(-4.*acr/av)
 endif
end function

function atrans(fm,av)
  real :: fm, atrans,av
 if(poly.eq.0) then
  atrans=3./4.*av*(1-fm**(1./3.))
 endif
 if(poly.eq.1) then
  atrans=-av/4.*log(fm)
 endif
end function

function meanAv(pe,z)
  real :: pe,z,meanAv
 if(pe.gt.pcut) then
  meanAv=.22*z*(nzero/100)*sqrt(pe/pzero)
 else
  meanAv=.22*z*(nzero/100)*pe/pzero*sqrt(pzero/pcut)
 endif 
end function

subroutine solveh2(alphaG,dbeta,sn)
 real, intent(in) :: alphaG, dbeta
 real, intent(inout) :: sn
 real, parameter :: precision=0.0001
 real :: fi,dfi,sn1, sn2,esn
 integer :: i
 integer,parameter :: imax=100
  
 i=0 
 sn2=sn 
10 continue
 i=i+1
 sn1=sn2
 esn=exp(-sn1)
 fi=sn1-sn-dbeta*(alphaG*esn-sn1)
 dfi=1+dbeta*(alphaG*esn+1)
 sn2=sn1-fi/dfi 
 if(i.gt.imax) return
 if(ABS(sn2-sn1).GT.precision*sn1) goto 10
 sn=sn2
end subroutine

subroutine solveh2log(alphaG,dbeta,snzero,sn)
 real, intent(in) :: alphaG, dbeta,snzero
 real, intent(inout) :: sn
 real, parameter :: precision=0.0001
 real :: fi,dfi,sn1, sn2,esn,esn1
 integer :: i
 integer,parameter :: imax=100
  
 i=0 
 if(sn.gt.50*snzero) sn=50*snzero
 sn2=sn     
10 continue
 i=i+1
 if(sn2.gt.50*snzero) sn2=50*snzero
 sn1=sn2
 esn=exp(-sn1)
 esn1=exp(sn1/snzero)
 fi=sn1-sn-dbeta*(alphaG*esn-snzero*(esn1-1))
 dfi=1+dbeta*(alphaG*esn+esn1)
 sn2=sn1-fi/dfi 
 if(i.gt.imax) return
 if(ABS(sn2-sn1).GT.precision*sn1) goto 10
 sn=sn2
end subroutine

subroutine solveh2log_explicit(alphaG,dbeta,snzero,sn)
 real, intent(in) :: alphaG, dbeta,snzero
 real, intent(inout) :: sn
 real, parameter :: fac=0.001
 real :: fi,dfi,sn1, sn2,esn,esn1,dtau,tau

 sn1=sn
 tau=0
 do while(tau.LT.dbeta)
   esn=exp(-sn1)
   esn1=exp(sn1/snzero)
   dfi=(alphaG*esn-snzero*(esn1-1))
   dtau=abs(fac*sn1/dfi)
   if(sn1.EQ.0) dtau=abs(1.e-8/dfi)
   if((dbeta-tau).LT.dtau) then
     dtau=dbeta-tau
     tau=dbeta
   else
     tau=tau+dtau
   endif
   if(dfi.EQ.0) then 
     exit
   endif
   sn1=sn1+dtau*dfi
 enddo
 sn=sn1
end subroutine

subroutine solveh2_explicit(alphaG,dbeta,sn)
 real, intent(in) :: alphaG, dbeta
 real, intent(inout) :: sn
 real, parameter :: fac=0.001
 real :: fi,dfi,sn1, sn2,esn,esn1,dtau,tau

 sn1=sn
 tau=0
 do while(tau.LT.dbeta)
   esn=exp(-sn1)
   dfi=(alphaG*esn-sn1)
   dtau=abs(fac*sn1/dfi)
   if(sn1.EQ.0) dtau=abs(1.e-8/dfi)
   if((dbeta-tau).LT.dtau) then
     dtau=dbeta-tau
     tau=dbeta
   else
     tau=tau+dtau
   endif
   if(dfi.EQ.0) then 
     exit
   endif  
   sn1=sn1+dtau*dfi
 enddo
 sn=sn1
end subroutine

subroutine fh2des(beta,fh2) ! collissional h2 destruction (gamma2=gamma1/20)
 real :: beta,fh2
 fh2=20*fh2*exp(-beta)/(19*fh2*exp(-beta)+(20-19*fh2))
end subroutine

 function taudiss(T,n,xe) result(tau)
  real :: tau
  real, intent(in) :: T,n,xe
  real :: tau1,tau2,lxe
  
  lxe=MIN(1.,xe)

  tau1=taudiss1(T,n)
  tau2=taudiss2(T,n)
  tau=1./((1-lxe)/tau1+lxe/tau2)
!  print*,tau,tau1,tau2
 end function 

 function taudiss1(T,n) result(tau)
  real :: tau
  real, intent(in) :: T,n
  real, parameter :: Tref=3162,tauref=3.2e10,Tmin=3000,Tmax=10000
!  if(T.LT.Tmin) then
!   tau=tauref*(Tref/Tmin)**12/n 
!   return
!  endif
  if(T.GT.Tmax) then
   tau=tauref*(Tref/Tmax)**12/n
   return
  endif  
  tau=tauref*(Tref/T)**12/n
 end function

 function taudiss2(T,n) result(tau)
  real :: tau
  real, intent(in) :: T,n
  real, parameter :: Tref=3162,tauref=3.2e7,Tmin=3000,Tmax=10000
!  if(T.LT.Tmin) then
!   tau=tauref*(Tref/Tmin)**12/n 
!   return
!  endif
  if(T.GT.Tmax) then
   tau=tauref*(Tref/Tmax)**6/n
   return
  endif  
  tau=tauref*(Tref/T)**6/n
 end function

 function polyQ()
  integer polyQ
  polyQ=poly
 end function 

end module

subroutine h2test
 use MoleculeMod
 real :: fh2,dt,n,T,g0,disp,p
 real :: tform,feq,tdiss,tdes,xe
 integer :: i,j,k,l
 
 print*,'molecule formation routine test'

 call InitMolecule

print*,' fh2,dt,n,T,disp,g0,xe'
10 read*,fh2,dt,n,T,disp,g0,xe
 if(n.eq.0) stop 
  call h2fraction(fh2,dt,n,T,g0,disp,xe=xe)
 print*,dt,fh2
 goto 10
end 
