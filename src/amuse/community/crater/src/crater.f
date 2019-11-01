      function crater(rhoproj, L, v,
     &   theta, rhotarget, g, targtype, cratertype, Dyield,
     &   Dgault, Tform, Dfinal, Lyield,
     &   Lgault)
c      program crater
c
c     Short program to evaluate the scaling equations to determine
c     the diameter of a transient crater given details on the nature
c     of the projectile, conditions of impact, and state of the
c     target.  The diameter is evaluated by three independent methods,
c     yield scaling, pi-scaling and Gault's semi-empirical relations
c     supplemented by rules on how crater size depends on gravity and
c     angle of impact.
c
c     Updated Nov. 1997 to compute projectile size from a given 
c     transient crater diameter.  Projectile and crater diameter 
c     computation functions merged into a single program April 1998.
c     See Melosh, Impact Cratering, chapter 7 for more details
c
c     Corrected March 2016 in accord with scaling analysis by
c     Johnson et al. (2016) Icarus 271, 350-359 who argued that
c     simple/transient crater ratio (gama) is 1.25, not 1.56 and 
c     that the power in ratio of final crater to s/c transition (eta)
c     crater is 0.13, not 0.18. 
c
c     Copyright 1996, 1997, 1998 and 2016 by H. J. Melosh
c
      implicit none
      real*8 Cd(3),beta(3),Ct,Dt,Dstd,L,rhoproj,v,theta,rhotarget,g,W
      real*8 pi,third,pitwo,pifac,gearth,gmoon,m,dscale,anglefac,densfac
      real*8 rhomoon,Dstarmoon,Dstar,Dsimple,Dpr,Dprmoon,Dfinal
      real*8 Dpiscale,Dyield,Dgault,gsmall,Tform
      real*8 Lpiscale,Lyield,Lgault
      real*8 gama,eta
c      character ans,cratertype*14
      integer targtype,comptype, cratertype, crater
      pi=3.14159265
      third=1./3.
c
c     Initialize constants in AMUSE
      comptype = 1
      Dt = 0.0
      crater = 1

      write (*,*) "all parameters:", rhoproj, L, v,
     &   theta, rhotarget, g, targtype, cratertype, Dyield,
     &   Dgault, Tform, Dfinal, Lyield,
     &   Lgault      
c    constants for the Schmidt-Holsapple pi scaling,
c    gravity conversion factors and simple/complex
c    transition (see Johnson et al 2016)
c
      data   Cd/1.88,1.54 ,1.6 /
      data beta/0.22,0.165,0.22/
      data gearth,gmoon/9.8,1.67/
      data rhomoon,Dstarmoon,Dprmoon/2700.,1.8e4,1.4e5/
      data gama, eta /1.25, 0.13/
c
c***********************************************************************
c
c     get input conditions from the user's keyboard
c
c***********************************************************************
  10  continue
c$$$      write(*,999)
c$$$      read(*,*) comptype
c$$$	  if((comptype.ne.1).and.(comptype.ne.2)) then
c$$$	    write(*,1010)
c$$$		goto 10
c$$$      endif
c$$$      if(comptype.eq.2) then
c$$$        write(*,1400)
c$$$        read(*,*) Dt
c$$$		if(Dt.eq.0.0) then
c$$$          write(*,1500)
c$$$          read(*,*) Dfinal
c$$$		endif
c$$$      endif
c$$$      write(*,1000)
c$$$      read(*,*) rhoproj
c$$$      if(comptype.eq.1) then
c$$$        write(*,2000)
c$$$        read(*,*) L
c$$$      endif
c$$$      write(*,3000)
c$$$      read(*,*) v
c$$$      write(*,4000)
c$$$      read(*,*) theta
c$$$      write(*,5000)
c$$$      read(*,*) rhotarget
c$$$      write(*,6000)
c$$$      read(*,*) g
c$$$  20  write(*,7000)
c$$$      read(*,*) targtype
c$$$      if((targtype.lt.1).or.(targtype.gt.3)) then
c$$$        write(*,7010)
c$$$        goto 20
c$$$      endif
c
c     convert units to SI and compute some auxiliary quantites
c
      v=1000.*v                         !km/sec to m/sec
      Dt=1000.*Dt                       !km to m
      Dfinal=1000.*Dfinal               !km to m
      theta=theta*(pi/180.)             !degrees to radians
      anglefac=(sin(theta))**third      !impact angle factor
      densfac=(rhoproj**0.16667)/sqrt(rhotarget)
      pifac=(1.61*g)/v**2               !inverse froude length factor
      Ct=0.80                           !coefficient for formation time
      if(targtype.eq.1) Ct=1.3
      Dstar=(gmoon*rhomoon*Dstarmoon)/(g*rhotarget) !transition crater diameter
      Dpr  =(gmoon*rhomoon*Dprmoon  )/(g*rhotarget) !peak-ring crater diameter
c
c***********************************************************************
c
c          computation for specified projectile diameter
c
c***********************************************************************
      if(comptype.eq.1) then
        m=(pi/6.)*rhoproj*L**3          !projectile mass
        W=0.5*m*v**2                    !projectile kinetic energy
        pitwo=pifac*L                   !inverse froude number
        dscale=(m/rhotarget)**third     !scale for crater diameter
c
c     Pi Scaling (Schmidt and Holsapple 1987)
c
        Dpiscale=dscale*Cd(targtype)*pitwo**(-beta(targtype))
        Dpiscale=Dpiscale*anglefac
c
c     Yield Scaling (Nordyke 1962) with small correction for depth
c     of projectile penetration
c
        Dyield=0.0133*W**(1/3.4)+1.51*sqrt(rhoproj/rhotarget)*L
        Dyield=Dyield*anglefac*(gearth/g)**0.165
c
c     Gault (1974) Semi-Empirical scaling
c
        gsmall=0.25*densfac*(W**0.29)*anglefac
        if(targtype.eq.3) gsmall=0.015*densfac*(W**0.37)*anglefac**2
        if(gsmall.lt.100.) then
          Dgault=gsmall
        else
          Dgault=0.27*densfac*(W**0.28)*anglefac
       endif
       Dgault=Dgault*(gmoon/g)**0.165
c
c    Compute crater formation time from Schmidt and Housen
c
        Tform=(Ct*L/v)*pitwo**(-0.61)
c
c     Compute final crater type and diameter from pi-scaled transient dia.
c
        Dsimple=gama*Dpiscale
		if (Dsimple.lt.Dstar) then
		  Dfinal=Dsimple
		  cratertype= 0!'Simple'
		else
		  Dfinal=Dsimple*(Dsimple/Dstar)**eta
		  cratertype=1!'Complex'
		endif
		if((Dsimple.lt.Dstar*1.4).and.(Dsimple.gt.Dstar*0.71)) 
     &             cratertype=2!'Simple/Complex'
	        if(Dfinal.gt.Dpr) cratertype=3!'Peak-ring'
c
c    Print out results
c
c$$$        write(*,8000) Dyield,Dpiscale,Dgault,Tform
c$$$        write(*,8100) Dfinal/1000.,cratertype
c
c***********************************************************************
c
c        computation for specified crater diameter
c
c***********************************************************************
      elseif(comptype.eq.2) then
c
c     convert input crater rim-to-rim diameter to transient crater dia.
c
        if(Dt.eq.0.) then
          if(Dfinal.lt.Dstar) then
		    Dt=Dfinal/gama
		  else
		    Dt=(1./gama)*(Dfinal*Dstar**eta)**(1./(1.+eta))
		  endif
		endif
        dscale=((6.*rhotarget)/(pi*rhoproj))**third
c
c     Pi Scaling (Schmidt and Holsapple 1987)
c
        Dstd=Dt/anglefac
        Lpiscale=(Dstd*dscale*pifac**beta(targtype))/Cd(targtype)
        Lpiscale=Lpiscale**(1./(1.-beta(targtype)))
c
c     Yield Scaling (Nordyke 1962) without correction for projectile
c     penetration depth.
c
        Dstd=(Dt*(g/gearth)**0.165)/anglefac
        W=(Dstd/0.0133)**3.4
        Lyield=((12.*W)/(pi*rhoproj*v**2))**third
c
c     Gault (1974) Semi-Empirical scaling
c
        Dstd=Dt*(g/gmoon)**0.165
        if((Dstd.le.10.).and.(targtype.eq.3)) then
          W=((Dstd/0.015)/(densfac*anglefac**2))**2.70
        elseif(Dstd.lt.300.) then
          W=((Dstd/0.25)/(densfac*anglefac))**3.45
        else
          W=((Dstd/0.27)/(densfac*anglefac))**3.57
        endif
        Lgault=((12.*W)/(pi*rhoproj*v**2))**third
c
c    Compute crater formation time for Pi-scaled diameter
c
        Tform=(Ct*Lpiscale/v)*(pifac*Lpiscale)**(-0.61)
c
c    Print out results
c
c$$$        write(*,8500) Lyield,Lpiscale,Lgault,Tform
c$$$      else
c$$$        write(*,9500) comptype     !user should never get here
      endif
c
c***********************************************************************
c
c     Check to see if user wants to quit
c
c***********************************************************************
c$$$      write(*,9000) 
c$$$      read(*,*) ans
c$$$      if((ans.eq.'Y').or.(ans.eq.'y')) goto 10
c$$$      stop
c
  999 format(/
     &'            *** IMPACT SIZE ***'//
     &'This is a program to estimate the size of a gravity-'/
     &'dominated impact crater or the projectile that made it.'/
     &'Three different estimates are presented, but the '/
     &'pi-scaling method is currently considered the best!'//
     &'   enter the type of computation desired (1 or 2):'/
     &'      type 1, crater size'/
     &'      type 2, projectile size....................... ',$)
 1010 format(/
     &'   entry error! computation type must be 1 or 2.  Try again!')
 1400 format(/
     &'Crater descriptor:'//
     &'  enter the transient crater diameter in km'/
     &'  (if the final, not the transient crater diameter'/
     &'   is known, enter zero (0.0) here)................. ',$)
 1500 format(/
     &'  enter the final crater diameter in km............. ',$)
 1000 format(/
     &'Projectile descriptors:'//
     &'  enter the projectile density in kg/m^3............ ',$)
 2000 format(
     &'  enter the projectile diameter in m................ ',$)
 3000 format(/
     &'Impact conditions:'//
     &'  enter the impact velocity in km/sec............... ',$)
 4000 format(
     &'  enter the impact angle in degrees................. ',$)
 5000 format(/
     &'Target descriptors:'//
     &'  enter the target density in kg/m^3................ ',$)
 6000 format(
     &'  enter the acceleration of gravity in m/sec^2...... ',$)
 7000 format(
     &'  enter the target type, (1-3):'/
     &'    type 1 = liquid water'/
     &'    type 2 = loose sand  '/
     &'    type 3 = competent rock or saturated soil....... ',$)
 7010 format(//'target type must be 1, 2 or 3.  try again!'//)
 8000 format(//'Three scaling laws yield the following *transient*',
     &' crater diameters:'/
     &'  (note that diameters are measured at the pre-impact surface.'/
     &'   Rim-to-rim diameters are about 1.25X larger!)'//,
     &' Yield Scaling                       Dyield   = ',1pe10.2,' m'/
     &' Pi Scaling    (Preferred method!)   Dpiscale = ',1pe10.2,' m'/
     &' Gault Scaling                       Dgault   = ',1pe10.2,' m'//
     &' Crater Formation Time               Tform    = ',1pe10.2,' sec'
     & //)
 8100 format(
     &'Using the Pi-scaled transient crater, the *final* crater has'/
     &'rim-to-rim diameter=',1pe10.2,' km, and is of type: ',a14//)
 8500 format(//'Three scaling laws yield the following projectile',
     &' diameters:'/
     &'  (note that diameters assume a spherical projectile)'//
     &' Yield Scaling                       Lyield   = ',1pe10.2,' m'/
     &' Pi Scaling    (Preferred method!)   Lpiscale = ',1pe10.2,' m'/
     &' Gault Scaling                       Lgault   = ',1pe10.2,' m'//
     &' Crater Formation Time               Tform    = ',1pe10.2,' sec'
     & //)
 9000 format('Do you want to make another calculation (y/n)?......',$)
 9500 format('unrecognized computation type = ',i5,
     &' (should be either 1 or 2)')
      end
