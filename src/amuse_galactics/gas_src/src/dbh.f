      program dbh

      integer nrmaxsize,nr,lmax,nrmax
      
      real apot,fr,dr,potcor,plcon
      real psi0,psi2,psi4,dpsi,ds1s2,dapot
      real a,b,c,v0,q,psi00,psiout,rho1,sigbulge2,rmdisk,rdisk,
     & zdisk,outdisk,drtrunc,potcor1,v02,v03,rdisk2,diskconst,bulgea
      integer idiskflag,ibulgeflag,ihaloflag
            
      parameter(nrsize=20000)

      common /potconstants/ apot(20,0:nrsize), fr(20,0:nrsize), 
     +     dr, nr, lmax, potcor
      common /legendre/ plcon(0:40)
      common /gparameters/  a, b, c, v0, q, psi00, 
     +                      psiout, rho1, sigbulge2, 
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea
      common /flags/ idiskflag, ibulgeflag, ihaloflag

      real densnew,densold
      real gasm,gasr,gastr,gasdt,csnd,gamma,gassig
      integer gasnr
      common /gasstuf/ gasm,gasr,gastr,gasdt,csnd,gamma,gasnr,gassig,gasz

      parameter(pi=3.1415926535)
      real adens(20,0:nrsize),s1(0:nrsize),s2(0:nrsize)
      real rref(10),zref(10),pref(10),oldpref(10)
      real fr2(20,0:nrsize)
      real rscl,rv,rhs
      real bulgerscl,bulgerv,bulgerhs      
      external alldensnohalo,alldensnobulge,zerofunc
      
c  constants for legendre polynomials
      do l=0,40
         plcon(l) = sqrt((2*l+1)/(4.0*pi))
      enddo
c  set flags
      idiskflag = 0;ibulgeflag = 0;ihaloflag  = 0
c set default parameters
c disk: 
      rmdisk = 0.0;rdisk  = 1.0;zdisk  = .15
      outdisk = 5.0;drtrunc = .3
c bulge:
      rho1 = 0.;psiout = -2.0;sigbulge = .4
c halo:
      psi00 = 0.;q = 1.0
      v0 = 1.0;rking = 1

c Enter halo parameters
      write(*,*) 'Include halo? (0=no,1=dynamic,2=fixed,3=spherical df)'
      read(*,*) ihaloflag

      if( ihaloflag.eq.2) then
       print*, 'using MODEL.MASS file'
       q=1
      endif

      if( ihaloflag.eq.1) then
         write(*,*) 'central potential, v0, q, coreparam (=Rc/rK)^2, halo Ra?'
         read(*,*) psi00,v0,q,coreparam,ra
         v02 = v0**2
         v03 = v0**3
         rhochalo=3*v02/(4*pi*ra**2)
         write(*,*) 'rhochalo = ', rhochalo
         r1=v0/sqrt(4*pi*rhochalo)*exp(-psi00/v02)
c  fix Evans's B constant such that the core radius R_c is prop. to the King radius of the model.
         a=(2/pi)**2.5*(1-q**2)/q**2/v03/r1**4
         b=0
         c=(2*q**2-1)/(4*pi**2.5*q**2*v0)/r1**2
      endif
      
      if (ihaloflag.eq.3) then
        print*, 'halo parameters? (rs,rvir,rhos)'
        read*, rscl,rv,rhs	
	call sethalodf(rscl,rv,rhs)
      endif
            
c Enter disk parameters
      write(*,*) 'Include disk? (0=no,1=dynamic,2=fixed,3=fixed+dynamic)'
      read(*,*) idiskflag
      if( idiskflag.eq.1.OR.idiskflag.EQ.3) then
         write(*,*) 'Disk mass, scale length, radius, scale height and trunc. width?'
         read(*,*) rmdisk, rdisk, outdisk, zdisk, drtrunc
         rdisk2 = rdisk*rdisk
         diskconst = rmdisk/(4.0*pi*rdisk2*zdisk)
         do iref=1,4
           rref(iref)=iref*rdisk*2
           zref(iref)=0
           enddo
         do iref=5,7
           rref(iref)=(iref-5)*rdisk*4
           zref(iref)=zdisk/2.
           enddo
         do iref=8,10
           rref(iref)=rref(iref-3)
           zref(iref)=zdisk*2.
           enddo
         open(18,file='refpoints.dat',status='unknown')
      endif
      if( idiskflag.EQ.2.OR.idiskflag.EQ.3) then
        write(*,*) 'gasdisk mass, scale length, radius and trunc with?'
	read(*,*) gasm,gasr,gastr,gasdt
        write(*,*) 'csound,gamma, gas nr,maxz?'
	read(*,*) csnd,gamma,gasnr,gasz
      endif
            
c Enter bulge parameters
      write(*,*) 'Include bulge? (0=no,3=sphericaldf)'
      read(*,*) ibulgeflag
      if( ibulgeflag .eq. 1) then
         print*,'check bulge=1'
         write(*,*) 'Bulge central density (ish), cutoff potential, velocity dispersion?'
         read(*,*) rho1,psiout,sigbulge
         sigbulge2=sigbulge**2
         bulgea=rho1*exp((psi00-psiout)/sigbulge2)
      endif

      if (ibulgeflag.eq.3) then
        print*, 'bulge parameters? (rs,rvir,rhos)'
        read*, bulgerscl,bulgerv,bulgerhs
	call setbulgedf(bulgerscl,bulgerv,bulgerhs,halormin(),halormax())	
      endif

      write(*,*) 'radial step size, no. of bins?'
      read(*,*) dr,nr
      write(*,*) 'max. azimuthal harmonic l?'
      read(*,*) lmaxx

      if( ihaloflag .eq. 1 ) then
         print*,'warning: check!' 
         apot(1,0)=psi00*sqrt(4.*pi)
         dens00=totdens(0.,0.)
         rking=sqrt(9*v0**2/(8*pi*dens00))
c  find Rc which gives desired Rc/RK by bisection
         if (coreparam.le.0.) then
            b=0
         else
            rclo=0
            rchi=rking*sqrt(coreparam)
            rcmid=rchi/2.
            do while (rchi-rclo.gt.1.e-5*rcmid)
               b=sqrt(2/pi**5)*rcmid**2/q**2/v0/r1**4
               dens00=totdens(0.,0.)
               rcnew=sqrt(coreparam)*sqrt(9*v02/(8*pi*dens00))
               if (rcnew.gt.rcmid) then
                  rclo=rcmid
               else
                  rchi=rcmid
               endif
               rcmid=(rclo+rchi)/2.
            enddo
            b=sqrt(2/pi**5)*rcmid**2/q**2/v0/r1**4
            dens00=totdens(0.,0.)
            rking=sqrt(9*v02/(8*pi*dens00))
         endif
c  report final King radius and central density         
         write(*,*) 'central density=',dens00
         write(*,*) 'King radius=',rking
      endif


      call pgbegin(0,'/XWINDOW',5,4)
      call pgask(.false.)
      call pgsch(2.)
      write(*,*)

! init gas disk
      if( idiskflag.EQ.2.OR.idiskflag.EQ.3) then
        call setgasdisk()
! calculate 
! set plausible starting density for gas disk:
! (constant? stellar z-height? zero thickness?)
        call firstgasdens(zdisk)        
        densold=thickdiskdens(gasr,0.)
      endif

! calculate first guess potential:
      if(ihaloflag.eq.3) then
        call inithalodf(alldensnohalo)
!       psi00=0
       psi00=halodfpsi00()
      endif

c  cold start
c  initial potential, a good initial guess is important!
      if( ihaloflag .eq. 0 ) then
         q = 1
         psi00 = 0.0
         rcmid = 1.0 
         v02 = 1.0
      endif
      z = 0.0
      do ir=0,nr
         r=ir*dr
         if(ihaloflag.EQ.1) 
     &	 apot(1,ir)=psi00*sqrt(4.*pi)*exp(-100*(r*r+z*z)/(nr*dr)**2)
         if(ihaloflag.EQ.3)
     &	 apot(1,ir)=halodffi(r)*sqrt(4.*pi)-
     &   (rmdisk/sqrt(rdisk**2+r**2)-rmdisk/rdisk)*sqrt(4.*pi)
! fixed external potentials are added in pot function !!
!         print*,r,halodffi(r)
	 if(ihaloflag.EQ.2) apot(1,ir)=0.                  
	 do l=2,lmaxx,2
            apot(l/2+1,ir)=0
         enddo
      enddo

      lmaxstep=2
      lmax=0
      lmaxold=-2
      iter=0

! -------------------begin gas loop

!        print*,thickdiskdens(gasr,0.)
1      if( idiskflag.EQ.2.OR.idiskflag.EQ.3) then
!  begin with guess gas disk density:
         call solveZ    
!        print*,thickdiskdens(gasr,0.) 
! then potential:
         call solvepot
! repeat (because gas disk is probably subdominant):
         call solveZ        
!	 print*,thickdiskdens(gasr,0.)
         call solvepot
       endif

! need to redo:
      if(ihaloflag.eq.3) then
        call inithalodf(alldensnohalo)
        psi00=halodfpsi00()
      endif
! halo bigger than bulge assumed!     
      print*,'psi00:',psi00

      if(ibulgeflag.eq.3) then
        call initbulgedf(alldensnobulge)
        print*,bulgedfpsi00()
!       psi00=sdfpsi00()
      endif
      
! reset apot to psi00 at center:
	 a00=apot(1,0)
	 epot=extpot(0.,0.)
         print*,'reset central potential'
	 print*,'     psi00,pot,expot:',psi00,a00/sqrt(4.*pi),epot
         if(ihaloflag.EQ.1.OR.ihaloflag.EQ.3.OR.ibulgeflag.EQ.3) then
         do ir=0,nr
            apot(1,ir)=apot(1,ir)+(psi00-epot)*sqrt(4.*pi)-a00
         enddo
	 else
!         do ir=0,nr
!            apot(1,ir)=apot(1,ir)+fixedpot(0.)*sqrt(4.*pi)-a00
!         enddo	    
	 endif

      ntheta=(lmax)*4+2
      ntheta=max(10,ntheta)

c  now iterate. number of iterations and change in lmax depends on initial conditions.
c  iterate on first cycle until tidal radius is stable, then once for every harmonic added,
c  and until convergence once l has reached the desired maximum lmax.
      
c  if only a disk potential
      if( ihaloflag .eq. 2 .and. ibulgeflag .eq. 0 ) lmax = lmaxx-2
      drtidal = 1.
      rtidalold = 1e10
      iteroutside = 0

! get gas potential: will be added as external potential! 

      do while(iter.LE.299)
         iter=iter+1
c  determine number of harmonics, and number of azimuthal bins, to use this iteration
         if( lmax .eq. 0 .or. lmax .eq. lmaxx ) then
            if( drtidal .lt. 0.0005) then
               lmax=lmax+lmaxstep
               ntheta=lmax*4+4
            endif
         else
            lmax=lmax+lmaxstep
            ntheta=lmax*4+4
         endif
         if( lmax .eq. lmaxx+lmaxstep ) then
	  lmax=lmaxx
	  goto 199
         endif
c  Now get the harmonics of the density in this potential --> adens
c  NB that dens only gives the density without an approximate sech**2
c  component --- that part of the potential is not represented by the 
c  harmonics.
c  The function dens(r,z) returns the density - high-frequency cpt
c  The full density is returned by totdens
         adens(1,0)=dens(0.,0.)*sqrt(4.*pi)
         do l=2,lmax,2
            adens(l/2+1,0)=0
         enddo
         do ir=1,nr
            adens(1,ir)=0
         enddo
c  nrmax will mark the outermost radial bin with non-zero density.
         nrmax=nr
         do l=0,lmax,2
c  integrate density * spherical harmonic function over quadrant
c  use cos(theta) as independent variable.
            do ir=1,nrmax
               r=ir*dr
               s=0
               dctheta=1.0/ntheta
               s=s+polardens(r,1.0,l)+polardens(r,0.0,l)
               do is=1,ntheta-1,2
                  ctheta=is*dctheta
                  s=s+4*polardens(r,ctheta,l)
               enddo
               do is=2,ntheta-2,2
                  ctheta=is*dctheta
                  s=s+2*polardens(r,ctheta,l)
               enddo
               s=s*dctheta/3.
               s=s*4*pi
	       adens(l/2+1,ir)=s
c  mark the first even radial bin on which the density has fallen to zero.
	       if (l.eq.0 .and. s.eq.0.) then
                  nrmax=nr
                  nrmax=nrmax+mod(nrmax,2)
                  goto 77
               endif
            enddo	    
 77      enddo
 
c  now get the potential harmonics of this new density. (BT 2-208)
c  Simpson's rule integration. 
         do l=0,lmax,2
            s1(0)=0
            r = 2*dr
            s1(2)=(r*dr/3.)*(4*adens(l/2+1,1)*(1.0-dr/r)**(l+2)+adens(l/2+1,2))
            rold = r 
c  doesn't matter but should be nonzero
c  to avoid round-off error
            do ir=4,nr,2
               r=ir*dr
               s1a = (r*dr/3.)*(adens(l/2+1,ir-2)*(1.0-2*dr/r)**(l+2)+
     &              4*adens(l/2+1,ir-1)*(1.0-dr/r)**(l+2)+adens(l/2+1,ir))
               s1(ir) = s1a + s1(ir-2)*(rold/r)**(l+1)
               rold = r
            enddo
            s2(nr)=0
            rold = nr*dr
            do ir=nr-2,2,-2
               r=ir*dr
               s2a = (r*dr/3.)*(adens(l/2+1,ir+2)*(1.0+2*dr/r)**(1-l)+
     &              4*adens(l/2+1,ir+1)*(1.0+dr/r)**(1-l)+adens(l/2+1,ir))
               s2(ir) = s2a + s2(ir+2)*(r/rold)**l
!	       print*,s2a,s2(ir),adens(l/2+1,ir),ir                  
	       rold = r
            enddo
!	    stop
c  replace the potential harmonics with a mean of the previous iteration (25%) 
c  and the current one (75%). This damps out oscillations that otherwise occur.
c  if this is the first time this harmonic is calculated, use the entire new 
c  value.
            if(l.eq.0) then
             epot=extpot(0.,0.)
!	     print*,apot(1,2)/sqrt(4*pi),apot(1,nr)/sqrt(4*pi)
!	     print*,-(s1(2)+s2(2))*sqrt(4*pi),-(s1(nr)+s2(nr))*sqrt(4*pi)
C do the correction to psi00 here!	  
             psi2=-(s1(2)+s2(2))*sqrt(4*pi)
	     psi4=-(s1(4)+s2(4))*sqrt(4*pi)  
 	     psi0=2*psi2-psi4
	     dpsi=psi00-epot-psi0
	     ds1s2=-dpsi/sqrt(4*pi)
!	     print*,'a',ds1s2,s1(2)+s2(2)
	    else
	     ds1s2=0. 
	    endif
            do ir=2,nr,2
               if (l.le.lmaxold) then
                  apot(l/2+1,ir)=0.5*apot(l/2+1,ir)-
     &		  0.5*4*pi/(2.*l+1.)*(s1(ir)+s2(ir)+ds1s2)
               else
                  apot(l/2+1,ir)=-4*pi/(2.*l+1.)*
     +                 (s1(ir)+s2(ir)+ds1s2)
               endif
            enddo

c  
c  Calculate the 1st and 2nd-order radial gradients
c  
            do ir=2,nr,2
               r = ir*dr
               fr(l/2+1,ir)=-4*pi/(2.*l+1.)*(-(l+1)*s1(ir) + l*s2(ir))/r
               fr2(l/2+1,ir)=-4*pi/(2.*l+1.)*
     +              ((l+1)*(l+2)*s1(ir)/r**2+ 
     +              l*(l-1)*s2(ir)/r**2 -(2*l+1)*adens(l/2+1,ir))
            enddo
         enddo
c  now interpolate the gaps 
c  first quadratically interpolate the monopole back to the origin.
c  the remaining multipoles are zero there.
!	 apot(1,0)=(4*apot(1,2)-apot(1,4))/3.
	 apot(1,0)=(2*apot(1,2)-apot(1,4))
!         print*,apot(1,0),apot(1,2)
!         apot(1,0)=apot(1,2)
	 fr(1,0)=0.0
         fr2(1,0)=2*fr2(1,2)-fr2(1,4)
         fr2(1,0)=fr2(1,2)
         do l=2,lmax,2
            apot(l/2+1,0)=0
            fr(l/2+1,0)=0
            fr2(l/2+1,0)=0
         enddo
c  then linearly interpolate other bins.
         do ir=1,nr-1,2
            do l=0,lmax,2
               apot(l/2+1,ir)=(apot(l/2+1,ir-1)+apot(l/2+1,ir+1))/2.
               fr(l/2+1,ir)=(fr(l/2+1,ir-1)+fr(l/2+1,ir+1))/2.
               fr2(l/2+1,ir)=(fr2(l/2+1,ir-1)+fr2(l/2+1,ir+1))/2.
            enddo
         enddo
c plot current harmonic functions

c if you are only a disk potential then no need to iterate exit the loop
c
	 if( ihaloflag .eq. 0.and.ibulgeflag.eq.0 ) goto 199
c  finally reset the potential at the origin to psi00-extpot
c  Note that the fake disk potential is zero at the origin.
	 a00=apot(1,0)
	 epot=extpot(0.,0.)
!	 print*,'     psi00,pot,expot:',psi00,a00/sqrt(4.*pi),epot	 
	 if(ihaloflag.EQ.1.OR.ihaloflag.EQ.3.OR.ibulgeflag.EQ.3) then
         dapot=(psi00-epot)*sqrt(4.*pi)-a00
!	 print*,dapot
	 do ir=0,nr
            apot(1,ir)=apot(1,ir)+dapot
         enddo
	 else
!         do ir=0,nr
!            apot(1,ir)=apot(1,ir)+fixedpot(0.)*sqrt(4.*pi)-a00
!         enddo	    
	 endif
!	 print*,'a00:',apot(1,0)/sqrt(4.*pi)

!         if (a00/sqrt(4.*pi).gt.psi00-epot) then
!	    print*,a00/sqrt(4*pi),psi00
!            write(*,'(''Iter'',i4,'': lmax='',i4,
!     +           '', tidal radius is infinite'')') iter,lmax
!            iteroutside = iteroutside + 1
!            if( iteroutside .gt. 20 ) then
!                write(*,'(''nr='',i4,'' is too small'',
!     +           '', try larger number of radial bins  - exiting program'')') nr
!                goto 12345
!            endif
!            drtidal = 2.0*dr
!         else
            potr=apot(1,0)/sqrt(4.*pi)
            if(ihaloflag.eq.2) then
	        rtidal=nr*dr
		drtidal=0.
		goto 9
	    else
	    do ir=1,nr
c  look for new tidal radius of the model
c  defined as the radius at the equator where the potential is zero
               
	       
	       potrm1=potr
               potr=pot(ir*dr,0.)
!	       print*,ir*dr,potr,potrm1
!	       ,apot,dr,nr,lmax)
               if (potrm1*potr.le.0.) then
                  dpot = potr - potrm1
                  if( dpot .eq. 0.0 ) then
                     rtidal = (ir-1)*dr
                  else
                     rtidal=(ir-1-potrm1/dpot)*dr
                  endif
                  write(*,'(''   Iter'',i4,'': lmax='',i4,
     +                 '', tidal radius is '',g15.6)') iter,lmax,rtidal
                  drtidal = abs(rtidal - rtidalold)/rtidal
                  rtidalold = rtidal
                  
		  goto 9
               endif	       
            enddo
!	    endif
            write(*,'(''Iter'',i4,'': lmax='',i4,
     +           '', tidal radius is outside grid'')') iter,lmax
            drtidal = 1.0
 9       endif
c write out the changes in the potential at the reference points 
c at this iteration.
         if( idiskflag .eq. 1.or.idiskflag.eq.3 ) then 
            do iref=1,10
               oldpref(iref)=pref(iref)
               pref(iref)=pot(rref(iref),zref(iref))
!	       ,apot,dr,nr,lmax)
               enddo
            write(18,'(2i3,10g12.4)') iter,lmax,
     +                  (pref(iref)-oldpref(iref),iref=1,10)
          endif
c  now repeat with this new potential!
         lmaxold=lmax
         call dbhplot(adens,0,nr,dr)
         call dbhplot(apot,0,nr,dr)
      enddo
c  end of the main loop
 199  continue
 
 
 ! goto begin gas loop if not happy
      if( idiskflag.EQ.2.OR.idiskflag.EQ.3) then       
       densnew=thickdiskdens(gasr,0.)
       if(abs(densnew-densold)/densnew.GT.0.005) then
!        print*,densnew,densold
        print*,'repeating'
        densold=densnew
	goto 1
       endif 
!        print*,densnew,densold
        print*,' done?'       
       call writegas
      endif	
	      
      call pgend
      if( idiskflag .eq. 1.or.idiskflag.eq.3 ) close(18)
c  write final density, potential and potential gradient harmonics
c
C wrong!
      totalmass = fr(1,nr)/sqrt(4*pi)*(dr*nr)**2
      write(*,*) 'psi00=',psi00,'Rt/RK=',rtidal/rking
      write(*,*) 'Total mass=', totalmass

c  
c  Calculate force and potential for halo only
c  
      lmax = lmaxx
      halomass = 0.0
      bulgemass = 0.0
      haloedge = 0
      bulgeedge = 0
      if( ihaloflag.NE.0) call halopotential(halomass, haloedge)
      
      if( ibulgeflag.NE.0) call bulgepotential(bulgemass, bulgeedge)

      diskmass = totalmass - halomass - bulgemass
      diskedge = outdisk + 2.0*drtrunc

      open(20,file='mr.dat',status='unknown')
      write(20,*) diskmass, diskedge
      write(20,*) bulgemass, bulgeedge
      write(20,*) halomass, haloedge
      close(20)

      open(11,file='dbh.dat',status='unknown')
      write(11,'('' # psi00,v0,q,a,b,c,dr,nr,lmax='')')
      if( ihaloflag.NE.3) then
       write(11,'('' #'',7g15.5,i6,i4)') psi00,v0,q,a,b,c,dr,nr,lmaxx
      else
       write(11,'('' #'',7g15.5,i6,i4)') psi00,rscl,rv,rhs,b,c,dr,nr,lmaxx
      endif
      write(11,'('' # bulge psic, rho1, sig:'')')
      if(ibulgeflag.NE.3) then
       write(11,'('' #'',3g15.5)') psiout,rho1,sigbulge
      else
       write(11,'('' #'',3g15.5)') bulgerscl,bulgerv,bulgerhs      
      endif
      write(11,'('' # Mdisk, rdisk, zdisk, outdisk, drtrunc'')')
      write(11,'('' #'',5g15.5)') rmdisk, rdisk, zdisk, outdisk, drtrunc
      write(11,'('' #'',3i5)') idiskflag, ibulgeflag, ihaloflag
      write(11,'('' #  OUTPUT FROM DBH8. TOTAL POTENTIAL.'')')

      do ir=0,nr
         write(11,'(8g16.8)') ir*dr,(adens(l/2+1,ir),l=0,lmaxx,2)
      enddo
      write(11,*)
      write(11,*)
      do ir=0,nr
         write(11,'(8g16.8)') ir*dr,(apot(l/2+1,ir),l=0,lmaxx,2)
      enddo
      write(11,*)
      write(11,*)
      do ir=0,nr
         write(11,'(8g16.8)') ir*dr,(fr(l/2+1,ir),l=0,lmaxx,2)
      enddo
      write(11,*)
      write(11,*)
      do ir=0,nr
         write(11,'(8g16.8)') ir*dr,(fr2(l/2+1,ir),l=0,lmaxx,2)
      enddo
      write(11,*)
      do ir=0,nr
         psi=pot(ir*dr,0.)
!	 ,apot,dr,nr,lmax)
         write(11,'(2g16.8)') ir*dr,diskdens(ir*dr,0.0,psi)
      enddo
      close(11)

      write(*,*) 'Final model written to file ''dbh.dat''.'
12345 continue
      
      end
