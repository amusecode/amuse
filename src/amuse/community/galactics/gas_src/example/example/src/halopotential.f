! tbd: ftheta
      subroutine halopotential(halomass, haloedge) 

      common /potconstants/ apot(20,0:20000), fr(20,0:20000), 
     +     dr, nr, lmax, potcor
      common /gparameters/  a, b, c, v0, q, psi00, 
     +                      psiout, rho1, sigbulge2, 
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea
      common /flags/ idiskflag, ibulgeflag, ihaloflag

      real halopot(20,0:20000), halofr(20,0:20000)

      parameter(pi=3.1415926535)
      real adens(20,0:20000),s1(0:20000),s2(0:20000)

!      idiskflag = 0
!      ibulgeflag = 0
!      ihaloflag = 1
c Now get the harmonics of the density in this potential --> adens
      ntheta=lmax*4+2
      adens(1,0)=halodens(0.,0.)*sqrt(4.*pi)
      do l=2,lmax,2
         adens(l/2+1,0)=0
         enddo
      do ir=1,nr
         adens(1,ir)=0
         enddo
c nrmax will mark the outermost radial bin with non-zero density.
      nrmax=nr
      do l=0,lmax,2
c integrate density * spherical harmonic function over quadrant
         do ir=1,nrmax
            r=ir*dr
            s=0
            dctheta=1.0/ntheta
            s=s+polarhalodens(r,1.0,l)+polarhalodens(r,0.0,l)
            do is=1,ntheta-1,2
               ctheta=is*dctheta
               s=s+4*polarhalodens(r,ctheta,l)
            enddo
            do is=2,ntheta-2,2
               ctheta=is*dctheta
               s=s+2*polarhalodens(r,ctheta,l)
            enddo
            s=s*dctheta/3.
            s=s*4*pi
            adens(l/2+1,ir)=s
c            if (l.eq.0 .and. s.eq.0.) then
c               nrmax=nr
c               nrmax=nrmax+mod(nrmax,2)
c               goto 77
c               endif
c      write(*,*) 'a(',l,') =',s,'at r=',r
            enddo
 77      enddo
c now get the potential harmonics of this new density. Simpson's rule integration. (BT 2-208)
      do l=0,lmax,2
         s1(0)=0
         do ir=2,nr,2
            r=ir*dr
            s1(ir)=s1(ir-2)+(dr/3.)*
     +           (adens(l/2+1,ir-2)*(r-2*dr)**(l+2)+4*adens(l/2+1,ir-1)*(r-dr)**(l+2)+
     +           adens(l/2+1,ir)*r**(l+2))
            enddo
         s2(nr)=0
         do ir=nr-2,2,-2
            r=ir*dr
            s2(ir)=s2(ir+2)+(dr/3.)*
     +           (adens(l/2+1,ir+2)*(r+2*dr)**(1-l)+4*adens(l/2+1,ir+1)*(r+dr)**(1-l)+
     +           adens(l/2+1,ir)*r**(1-l))
            enddo
         do ir=2,nr,2
            halopot(l/2+1,ir)=-4*pi/(2.*l+1.)*
     +              (s1(ir)/(ir*dr)**(l+1)+s2(ir)*(ir*dr)**l)
            enddo
c
c  Calculate the radial gradients
c
         do ir=2,nr,2
            halofr(l/2+1,ir)=-4*pi/(2.*l+1.)*
     +              (-(l+1)*s1(ir)/(ir*dr)**(l+2)+ l*s2(ir)*(ir*dr)**(l-1))
            enddo
         enddo
c now interpolate the gaps 
c first quadratically interpolate the monopole back to the origin.
c the remaining multipoles are zero there.
      halopot(1,0)=(4*halopot(1,2)-halopot(1,4))/3.
      fr(1,0)=0.0
      do l=2,lmax,2
         halopot(l/2+1,0)=0
         halofr(l/2+1,0)=0
         enddo
      do ir=1,nr-1,2
         do l=0,lmax,2
            halopot(l/2+1,ir)=(halopot(l/2+1,ir-1)+halopot(l/2+1,ir+1))/2.
            halofr(l/2+1,ir)=(halofr(l/2+1,ir-1)+halofr(l/2+1,ir+1))/2.
            enddo
         enddo
c finally reset the potential at the origin to psi00
c      write(*,*) 'Central monopole:',apot(1,0),apot(1,1),apot(1,2),apot(1,3)
      a00=halopot(1,0)
      do ir=0,nr
         halopot(1,ir)=halopot(1,ir)-a00
         enddo
c write halo final potential.
c
c Reset the potentials so that phi is 0 at infinity
c
      redge = nr*dr
      do l=0,lmax,2
          c0 = halopot(l/2+1,nr) + halofr(l/2+1,nr)*redge/(l+1)
          do i=0,nr
              halopot(l/2+1,i) = halopot(l/2+1,i) - c0
          enddo
      enddo

      open(11,file='h.dat',status='unknown')

      write(11,'('' # psi00,v0,q,a,b,c,dr,nr,lmax='')')
      write(11,'('' #'',7g15.5,i6,i4)') psi00,v0,q,a,b,c,dr,nr,lmax
      write(11,'('' # bulge psiout, rho1, sig:'')')
      write(11,'('' #'',3g15.5)') psiout,rho1,sqrt(sigbulge2)
      write(11,'('' # Mdisk, rdisk, zdisk, outdisk, drtrunc'')')
      write(11,'('' #'',5g15.5)') rmdisk, rdisk, zdisk, outdisk, drtrunc
!      write(11,'('' #'',3i5)') idiskflag, ibulgeflag, ihaloflag
      write(11,'('' #'',3i5)') 0, 0, ihaloflag
      write(11,'('' #  OUTPUT FROM DBH8. HALO POTENTIAL ONLY.'')')
      do ir=0,nr
         write(11,'(8g16.8)') ir*dr,(adens(l/2+1,ir),l=0,lmax,2)
         enddo
      write(11,*)
      write(11,*)
      do ir=0,nr
         write(11,'(8g16.8)') ir*dr,(halopot(l/2+1,ir),l=0,lmax,2)
         enddo
      write(11,*)
      write(11,*)
      do ir=0,nr
         write(11,'(8g16.8)') ir*dr,(halofr(l/2+1,ir),l=0,lmax,2)
         enddo
      close(11)
      halomass = halofr(1,nr)/sqrt(4*pi)*(nr*dr)**2
      do ir=0, nr
          if( adens(1,ir) .eq. 0.0 ) then
              haloedge = ir*dr
              goto 999
          endif
      enddo
      haloedge = nr*dr
999   write(*,*) 'Halo mass =', halomass
      write(*,*) 'Halo edge radius =', haloedge
      write(*,*) 'Halo harmonics written to file ''h.dat''.'
      return
      end



