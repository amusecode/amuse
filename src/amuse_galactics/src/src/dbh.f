      include 'commonblocks'

      real xa(5),ya(5)

      character ans*1

      open(file='check.dat',unit=10,status='replace')

c  constants for legendre polynomials
      do l=0,40
         plcon(l) = sqrt((2*l+1)/(4.0*pi))
      enddo
c
c     Enter halo parameters
c
      write(*,*) 'Include halo?'
      read(*,'(a)') ans

      if(ans.eq.'y') then
         write(*,*) 'outer radius, v0, a, drtrunchalo?'
         read(*,*) chalo, v0, a, drtrunchalo, cusp, outerhaloslope
         ihaloflag = 1
         haloconst = (2.**(1.-cusp))*v0*v0/4./pi/a/a
      endif
      
c Enter disk parameters
      write(*,*) 'Include disk?'
      read(*,'(a)') ans
      if( ans .eq. 'y' ) then
         write(*,*)
     *        'Disk mass, scale length, rad, scale height,trunc width'
         read(*,*) rmdisk, rdisk, outdisk, zdisk, drtrunc
         idiskflag = 1
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
      else
         rmdisk = 0.
      endif

c Enter bulge parameters

      write(*,*) 'Include bulge?'
      read(*,'(a)') ans
      if( ans .eq. 'y' ) then
         write(*,*) 'sersic, ppp, v0_bulge, a_bulge?'
         read(*,*) nnn, ppp, v0bulge, abulge
         if(ppp.lt.0) then
            ppp = 1. - 0.6097/nnn + 0.05463/(nnn*nnn)
         endif
         ibulgeflag = 1
         call setsersicparameters
      else
         v0bulge = 0.
         abulge = 0.
         ibulgeflag = 0
      endif

c     Enter blackhole parameters
      write(*,*) 'Include blackhole?'
      read(*,'(a)') ans
      if(ans.eq.'y') then
         write(*,*) 'blackhole mass?'
         read(*,*) bhmass
         ibhflag = 1
      else
         bhmass = 0.
         ibhflag = 0
      endif

c     Enter grid parameters
      write(*,*) 'radial step size, no. of bins?'
      read(*,*) dr,nr
      write(*,*) 'max. azimuthal harmonic l?'
      read(*,*) lmaxx

      open(file='in.gendenspsi',unit=40,status='old')
      read(40,*) npsi,nint
      close(40)

      if(ihaloflag.eq.1) call halopotentialestimate
      if(idiskflag.eq.1) call diskpotentialestimate
      call gentableE
      if(ibulgeflag.eq.1) call gendfsersic
      if(ihaloflag.eq.1) call gendfnfw
      if(ihaloflag.eq.1) call gendenspsihalo
      if(ibulgeflag.eq.1) call gendenspsibulge

      open(file='dfhalo.table',unit=41,status='replace')
      open(file='denspsihalo.dat',unit=42,status='old')
      do i=1,npsi
         read(42,*) energy
         dfhaloval = dfhalo(energy)
         write(41,*) energy, dfhaloval
      enddo

      apot(1,0)=psi0*sqrt(4.*pi)

c     initial potential is from our potential ansatz

      z = 0.0
      do ir=1,nr
         r=ir*dr
         apot(1,ir) = sqrt(4.*pi)*gettotalpsi(r)
         do l=2,lmaxx,2
            apot(l/2+1,ir)=0
         enddo
      enddo

      niter=(2+lmaxx/2)*10
      lmaxstep=2
      lmax=0
      ntheta=lmax*10+2
      ntheta=max(10,ntheta)

c     now iterate. number of iterations and change in lmax depends on
c     initial conditions.  iterate on first cycle until tidal radius is
c     stable, then once for every harmonic added, and until convergence
c     once l has reached the desired maximum lmax.
      
      drtidal = 2*dr
      rtidalold = 1e10
      tidalcheck = rtidalold
      lmaxold=-2
      iteroutside = 0
      frac=0.75

      do iter=0,100
         if( lmax .eq. 0 .or. lmax .eq. lmaxx ) then
            if(drtidal .lt. dr .and. iter.gt.10) then
               lmax=lmax+lmaxstep
               ntheta=lmax*4+2
            endif
         else
            lmax=lmax+lmaxstep
            ntheta=lmax*4+2
         endif
         if( lmax .eq. lmaxx+lmaxstep ) goto 199

c     Now get the harmonics of the density in this potential --> adens
c     NB that dens only gives the density without an approximate sech**2
c     component --- that part of the potential is not represented by the
c     harmonics.  The function dens(r,z) returns the density -
c     high-frequency cpt The full density is returned by totdens

         eps = 0.00001
         adens(1,0)=dens(eps,0.)*sqrt(4.*pi)

         do l=2,lmax,2
            adens(l/2+1,0)=0
         enddo
         do ir=1,nr
            adens(1,ir)=0
         enddo

c     nrmx will mark the outermost radial bin with non-zero density.
         nrmx=nr
         do l=0,lmax,2

c     integrate density * spherical harmonic function over quadrant use
c     cos(theta) as independent variable.

            do ir=1,nrmx
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
c     mark the first even radial bin on which the density has fallen to
c     zero.
               if (l.eq.0 .and. s.eq.0.) then
                  nrmx=nr
                  nrmx=nrmx+mod(nrmx,2)
                  goto 77
               endif
            enddo
 77      enddo

c     now get the potential harmonics of this new density. (BT 2-208)
c     Simpson's rule integration.

         do l=0,lmax,2
            s1(0)=0
            r = 2*dr
            s1(2)=(r*dr/3.)*(4*adens(l/2+1,1)*(1.0-dr/r)**(l+2)
     +           +adens(l/2+1,2))
            rold = r 
            do ir=4,nr,2
               r=ir*dr
               s1a = (r*dr/3.)*(adens(l/2+1,ir-2)*(1.0-2*dr/r)**(l+2)+
     &              4*adens(l/2+1,ir-1)*(1.0-dr/r)**(l+2)+
     +              adens(l/2+1,ir))
               s1(ir) = s1a + s1(ir-2)*(rold/r)**(l+1)
               rold = r
            enddo
            s2(nr)=0
            rold = nr*dr
            do ir=nr-2,2,-2
               r=ir*dr
               s2a = (r*dr/3.)*(adens(l/2+1,ir+2)*(1.0+2*dr/r)**(1-l)+
     &              4*adens(l/2+1,ir+1)*(1.0+dr/r)**(1-l)+
     +              adens(l/2+1,ir))
               s2(ir) = s2a + s2(ir+2)*(r/rold)**l
               rold = r
            enddo

c     replace the potential harmonics with a mean of the previous
c     iteration (25%) and the current one (75%). This damps out
c     oscillations that otherwise occur.  if this is the first time this
c     harmonic is calculated, use the entire new value.

            do ir=2,nr,2
               if (l.le.lmaxold) then
                  apot(l/2+1,ir)=frac*apot(l/2+1,ir)+
     +                 (1.-frac)*4*pi/(2.*l+1.)*(s1(ir)+s2(ir))
               else
                  apot(l/2+1,ir)=4*pi/(2.*l+1.)*(s1(ir)+s2(ir))
               endif
            enddo

c     Calculate the 1st and 2nd-order radial gradients
  
            do ir=2,nr,2
               r = ir*dr
               fr(l/2+1,ir)=-4*pi/(2.*l+1.)*(-(l+1)*s1(ir) + l*s2(ir))/r
               fr2(l/2+1,ir)=-4*pi/(2.*l+1.)*
     +              ((l+1)*(l+2)*s1(ir)/r**2+ 
     +              l*(l-1)*s2(ir)/r**2 -(2*l+1)*adens(l/2+1,ir))
            enddo
         enddo

c     now interpolate the gaps first quadratically interpolate the
c     monopole back to the origin.  the remaining multipoles are zero
c     there.

c         npoly = 5
c         rzero = 0.
c         xa(1) = (2.)
c         xa(2) = (4.)
c         xa(3) = (6.)
c         xa(4) = (8.)
c         xa(5) = (10.)
c         ya(1) = apot(1,2)
c         ya(2) = apot(1,4)
c         ya(3) = apot(1,6)
c         ya(4) = apot(1,8)
c         ya(5) = apot(1,10)
c         call polint(xa,ya,npoly,rzero,apot(1,0),dy)
         apot(1,0)=3*(apot(1,2)-apot(1,4))+apot(1,6)

         fr2(1,0)=2*fr2(1,2)-fr2(1,4)
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

c  finally reset the potential at the origin to psi0
c  Note that the fake disk potential is zero at the origin.

         if (ihaloflag.eq.0.and.ibulgeflag.eq.0) goto 199

         a00=apot(1,0)
         do ir=0,nr
            apot(1,ir)=apot(1,ir)+psi0*sqrt(4.*pi)-a00
         enddo

         if (a00/sqrt(4.*pi)-psi0.gt.psic) then
            write(*,'(''Iter'',i4,'': lmax='',i4,
     +           '', tidal radius is infinite'')') iter, lmax 
            iteroutside = iteroutside + 1
            if( iteroutside .gt. 40 ) then
                write(*,'(''nr='',i4,'' is too small'',
     +           '', try larger number of radial bins  
     +              - exiting program'')') nr
                goto 12345
            endif
            drtidal = 2.0*dr
         else
            potr=psi0
            do ir=1,nr
c     look for new tidal radius of the model defined as the radius at
c     the equator where the potential is equal to psic
               potrm1=potr
               potr=pot(ir*dr,0.)
               aa = potrm1 - psic
               bb = potr - psic
               if (aa*bb.le.0.) then
                  dpot = potr - potrm1
                  if( dpot .eq. 0.0 ) then
                     rtidal = (ir-1)*dr
                  else
                     rtidal=(ir-1-aa/dpot)*dr
                  endif
                  drtidal = abs(rtidal - rtidalold)
                  tidalcheck = abs(rtidal - rtidalold)/rtidal
                  write(*,'(''Iter'',i4,'': lmax='',i4,
     +                 '', tidal radius is '',g15.6)') iter,lmax,rtidal
                  rtidalold = rtidal
                  goto 9
               endif
            enddo
            write(*,'(''Iter'',i4,'': lmax='',i4,
     +           '', tidal radius is outside grid'')') iter,lmax
            drtidal = 2.0*dr
 9       endif
c write out the changes in the potential at the reference points 
c at this iteration.
         if( idiskflag .eq. 1 ) then 
            do iref=1,10
               oldpref(iref)=pref(iref)
               pref(iref)=pot(rref(iref),zref(iref))
            enddo
            write(18,'(2i3,10g12.4)') iter,lmax,
     +           (pref(iref)-oldpref(iref),iref=1,10)
         endif
c     now repeat with this new potential!
         lmaxold=lmax

      enddo
c     end of the main loop
 199  continue

      open(file='rtidal.dat',unit=20,status='replace')
      if( idiskflag .eq. 1 ) close(18)
      totalmass = fr(1,nr)/sqrt(4*pi)*(dr*nr)**2
      write(*,*) 'Total mass=', totalmass
      write(20,*) rtidal
c  
c  Calculate force and potential for halo only
c  
      lmax = lmaxx

      halomass = 0.0
      bulgemass = 0.0
      haloedge = 0
      bulgeedge = 0
      if( ihaloflag .eq. 1 ) call halopotential(halomass, haloedge)
      if( ibulgeflag .eq. 1 ) call bulgepotential(bulgemass, bulgeedge)
      if( ibhflag. eq. 1) call genblackhole(bhmass)
      if( idiskflag.eq.1) then
         diskmass = totalmass - halomass - bulgemass
         diskedge = outdisk + 2.0*drtrunc
      else
         diskmass = 0.
         diskedge = 0.
      endif

      open(20,file='mr.dat',status='unknown')
      write(20,*) diskmass, diskedge
      write(20,*) bulgemass, bulgeedge
      write(20,*) halomass, haloedge
      close(20)

      open(11,file='dbh.dat',status='unknown')
      write(11,'('' # chalo,v0,a,nnn,v0bulge,abulge,dr,nr,lmax='')')
      write(11,'('' #'',7g15.5,i6,i4)') chalo,v0,a,nnn,v0bulge,abulge,
     +     dr,nr,lmaxx
      write(11,'('' # psi0, haloconst, bulgeconst:'')')
      write(11,'('' #'',3g15.8)') psi0,haloconst,bulgeconst
      write(11,'('' # Mdisk, rdisk, zdisk, outdisk, drtrunc'')')
      write(11,'('' #'',5g15.5)') rmdisk, rdisk, zdisk, outdisk, drtrunc
      write(11,'('' # psic, psi0-psid, bhmass'')')
      write(11,'('' #'',3g15.5)') psic,psi0-psid,bhmass
      write(11,'('' #'',4i5)') idiskflag, ibulgeflag, ihaloflag, ibhflag
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
         write(11,'(2g16.8)') ir*dr,diskdens(ir*dr,0.,psi)
      enddo
      close(11)

      write(*,*) 'Final model written to file ''dbh.dat''.'

      write(*,*) haloedge

      call check

12345 continue

      open(file='tidalr.out',status='replace',unit=14)
      write(14,*) rtidal
      
      end

      subroutine check

      include 'commonblocks'
      
      real psi

      open(file='check.out',unit=66,status='replace')

      z = 0.
      do ir=1,nr
         r = ir*dr
         psi = pot(r,z)
         bdens = bulgedenspsi(psi)
         ebdens = sersicdens(r)
         hdens = halodenspsi(psi)
         ehdens = halodensity(r)
         ddens = diskdens(r,z,psi)
         eddens = rmdisk/4./pi/rdisk/rdisk/zdisk*exp(-r/rdisk)
         tdens = totdens(r,z)
         write(66,88) 
     +        log10(r),
     +        log10(ddens),log10(eddens),
     +        log10(bdens),log10(ebdens),
     +        log10(hdens),log10(ehdens),psi
      enddo
 88   format(10f12.3)
      close(10)
      return
      end

      subroutine checkhalo(haloedge)

      include 'commonblocks'
      
      real psi

      open(file='check.out',unit=66,status='replace')

      z = 0.
      rlogmin = log10(dr)
      rlogmax = log10(haloedge)
      drlog = (rlogmax-rlogmin)/float(nr)
      do ir=1,nr
         r = ir*dr
         psi = pot(r,z)
         hdens = halodenspsi(psi)
         ehdens = halodensity(r)
         write(66,88) log10(r),log10(hdens),ehdens,psi
c         write(66,88) log10(r),log10(ebdens),log10(ebpot),
c     +         log10(ebdens1),log10(ebdens2)
      enddo
 88   format(8f12.3)
      close(10)
      return
      end

      subroutine genblackhole(bhmass)

      open(file='blackhole',unit=20,status='replace')
      n = 1
      t = 0.
      write(20,*) n,t
      write(20,*) bhmass,t,t,t,t,t,t
      close(20)
      
      return
      end

