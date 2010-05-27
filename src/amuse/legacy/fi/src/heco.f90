subroutine heattabel
  use ElementsMod
  use CoolingMod
  use IonMod
  use MoleculeMod
  use H2CoolMod
  use starsMod
  include "globals.h"
  integer hy

  if(verbosity.GT.0) print*,' ...InitElements...'
  call InitElements(metallicity)
  print*,'  > compiled metallicity:', zQ()/solarzQ()
  if(nstar.GT.0) then
    if(verbosity.GT.0)print*,' ...InitStars...'
    call InitStars(datadir,zQ())
  endif
  if(radiate) then
    if(verbosity.GT.0) print*,' ...InitCooling...'
    call InitCooling(datadir)
    call InitH2(datadir)
    if(verbosity.GT.0) print*,' ...InitIon...'
    call InitIon(unitm_in_msun,unitl_in_kpc)
    if(verbosity.GT.0) print*,' ...InitMolecules...'
    call InitMolecule
  endif

! solarz scaling with elements present 	
  hy=elementnrQ("H   ")
!  meanmwt=meanmwtQ()/fractionQ(hy)
  meanmwt=meanmwtQ()
  fhydrogn=wfractionQ(hy)
end subroutine

function heateff( t, fuv, ne)
  real fuv,t,ne,heateff,fuvtne
!  wolfire 1995 formula (2) / 0.05 (<- reference value for graineff)
  fuvtne=fuv*t**.5/ne
  heateff= 0.98 /(1+(fuvtne/1925)**.73) &
                  + 0.74*(t/1.e4)**.7/(1+(fuvtne/5000.))
  heateff=heateff*exp(-(t/2.e4)**2)
end function

function heateff2( t, fuv, ne)
  real fuv,t,ne,heateff2,fuvtne
! common sense grain heating efficiency? 
  heateff2=exp(-(t/2.e4)**2)
end function

function ffcooling(temp,p)
  include'globals.h'
  integer p
  real ffcooling,temp
  real gff,ffconst,xe
  parameter(gff=1.3,ffconst=1.42e-6)
  ffcooling=0
  if(temp.lt.1e5) return 
  ffcooling=gff*ffconst*cool_par*temp**.5*rho(p)*fhydrogn**2
end function

subroutine heco(heating,cooling,eth,p)
  use CoolingMod
  use IonMod
  use HeatingMod
  use H2coolMod
  include'globals.h'
  real cooling,heating,eth,nh
  real ttemp,temp,xelec,grainc,h2f,xh
  real heateff,heateff2
  integer p
	 
  nh=densconst*rho(p)
  h2f=h2frac(p)
  ttemp=eth*mumhkg1
  if(ionCool) then
    call Ion(ttemp,temp,xelec,nh,crionrate)
!    call Ionh(ttemp,temp,xelec,xhii,xheii,xheiii,nh,crionrate)
  else
    temp=ttemp
    xelec=0.1
  endif
	 
  if(simpleeff) then
    grainc=graineff*heateff2(temp,fuvheat(p)/heatconst,xelec*nh)
  else 
    grainc=graineff*heateff(temp,fuvheat(p)/heatconst,xelec*nh)
  endif
	 
  heating=fhydrogn*(fuvheat(p)*grainc+esnthdt(p)+ &
            grainc/graineff*heat_par1) + heat_par2*rho(p)*fhydrogn**2

  if(cosmrayheat) heating=heating+ &
                    fhydrogn*crionrate*1.6e-5*heatconst*crheat(xelec)

! constant (cosmological) heating term
!  heating=heating+fhydrogn*heat_par1*(1-xh)
         
!  cooling=cool_par*H2CoolFunc(xelec,temp,h2f)*rho(p)

  if(.NOT.H2Cooling) then
    cooling=cool_par*CoolFunc(xelec,temp)*rho(p)
  else
    cooling=cool_par*rho(p)*(CoolFunc(xelec,temp)+h2f*H2CoolFunc(temp,nh))
  endif

end subroutine
	 	 	
subroutine temperature
  use IonMod
  include 'globals.h'
  real ttemp,temp,xe,n
  integer p
	 
  if(isotherm) return 
!$omp parallel do private(ttemp,temp,xe,n)
  do p=1,nsph
    if(rho(p).EQ.0) cycle
!    if(uentropy) then
!      ttemp=mumhkg1*ethermal(p)/gamma1*rho(p)**gamma1
!    else
!      ttemp=mumhkg1*ethermal(p)
!    endif
    n=rho(p)*densconst
    ttemp=mumhkg1*csound(p)**2/gamma/gamma1	 
    if(radiate) then
      call Ion(ttemp,temp,xe,n,crionrate)
    else
      temp=ttemp
      xe=0.00001
    endif  
    temperat(p)=temp
    elecfrac(p)=xe
  enddo

end subroutine

subroutine molecules
  use MoleculeMod
  include 'globals.h' 
  real fh2,dt,n,T,g0,disp,lxe
  integer p

  dt=(tnow-h2time)*timescale/year
!$omp parallel do private(p,fh2,n,T,g0,disp,lxe) shared(dt)	
  do p=1,nsph
    if(rho(p).EQ.0) cycle
    fh2=h2frac(p)
    n=rho(p)*densconst
    T=temperat(p)
! 2pi/ 4pi and Habing -> Draine
    g0=fuvheat(p)/heatconst/2./1.71 
    disp=vdisp(p)*velscale/1.e5
    lxe=elecfrac(p)
    call h2fraction(fh2,dt,n,T,g0,disp,xe=lxe)	
    h2frac(p)=fh2  
  enddo

  h2time=tnow
end subroutine
        
