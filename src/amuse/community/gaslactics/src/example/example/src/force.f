c this file contains the following subroutines:
c
c readharmfile: reads a file with the harmonic expansion of the potential and density.
c all the relevant stuff ends up in a common block that only needs to be seen by
c the routines in this file (I think!). This routine must be called first.
c
c dens(r,z): density at point (r,z)
c
c densrpsi(r,psi): density as a function of r and potential psi.
c
c pot(r,z): potential at point (r,z)
c
c force(r,z,fr,fz): radial and z  force at (r,z)
c all the rest is (slightly nobbled in places) numerical recipes routines.
c


c
c calculate the force by interpolating the derived radial functions
c
c      subroutine force(lm,s,z,fs,fz,fs1,fz1)
      subroutine force(s,z,fs,fz,pot)
      parameter (pi=3.1415926535)
      parameter (pc0=0.282094792, pc2=0.630783131, pc4=0.846284375)
c constants cn are sqrt((2*l+1)/(4*pi))
      parameter (pc6=1.017107236, pc8=1.163106623)

      dimension pc(20), p(20), dp(20)
      common /potconstants/ apot(20,0:20000), frad(20,0:20000), 
     +     dr, nr, lmax, potcor
      common /legendre/ plcon(0:40)
      common /gparameters/ a,b,c,v0,q,psi00, 
     +                      psiouta, rho1a, sigbulge, 
     +                      grmdisk, grdisk, gzdisk, goutdisk, gdrtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea
      common /diskpars/ sigr0, disksr, nrdisk

      common /flags/ idiskflag, ibulgeflag, ihaloflag

c      write(*,*) 'Force entered with s,z=',s,z
c      write(*,*) 'lmax =',lmax
c      write(*,*) 'apot starts as:'
c      do ll=1,3
c         write(*,*) (apot(ll,i),i=1,10)
c      enddo
      r=sqrt(s*s+z*z)
      ihi=int(r/dr)+1
      if (ihi.lt.1) ihi=1
      if (ihi.gt.nr) ihi=nr
      r1=dr*(ihi-1)
      r2=dr*ihi
      redge = nr*dr
      t=(r-r1)/(r2-r1)
      ihim1 = ihi - 1
      tm1 = 1.0 - t
      if (r.eq.0.) then
         fs = 0.0
         fz = 0.0
      else
         costheta=z/r
         ct2 = costheta*costheta
         sintheta=s/r
         
         do l=0,lmax,2
            pc(l/2+1) = sqrt((2.0*l + 1)/(4.0*pi))
            p(l/2+1) = plgndr1(l,costheta) 
			if( costheta .eq. 1.0 ) then
               dp(l/2+1) = 0.0
            else
               st2 = 1.0 - costheta*costheta
               dp(l/2+1) = l*(plgndr1(l-1, costheta) - costheta*p(l/2+1))/st2
            endif
         enddo
         do i=1,lmax/2+1
            p(i) = p(i)*pc(i)
            dp(i) = dp(i)*pc(i)
         enddo

         if( r .le. redge ) then
             fr = 0.0
             do i=1,lmax/2+1
                fr = fr + p(i)*(t*frad(i,ihi) + tm1*frad(i,ihim1))
             enddo
             fth = 0.0
             do i=2,lmax/2+1
                fth = fth - sintheta*dp(i)*(t*apot(i,ihi) + tm1*apot(i,ihim1))
             enddo
             pot = 0.0
             do i=1,lmax/2+1
                pot = pot + p(i)*(t*apot(i,ihi) + tm1*apot(i,ihim1))
             enddo
         else
! expressions below for fr and pot not correct if apot shifted!!
             fr = 0.0
             do i=1,lmax/2+1
                 l = 2*(i-1)
                 fr = fr - (l+1)*p(i)*apot(i,nr)/redge*(redge/r)**(l+2)
             enddo
             fth = 0.0
             do i=2,lmax/2+1
                 l = 2*(i-1)
                 fth = fth - sintheta*dp(i)*apot(i,nr)*(redge/r)**(l+1)
	     enddo
             pot = 0.0
             do i=1,lmax/2+1
                l = 2*(i-1)
                pot = pot + p(i)*apot(i,nr)*(redge/r)**(l+1)
             enddo
         endif
         if( idiskflag .eq. 1.OR.idiskflag.EQ.3 ) then
             pot = pot + appdiskpot(s,z)
         endif
 
      if( idiskflag .eq. 2 .or. idiskflag .eq. 3) then
       pot= pot + thickdiskpot(s,z)      
!       fr=fr+thickdiskforce(r)
! fth (or rather fs and fz directly) tbd, see below
      endif

      if(ihaloflag.eq.2) then
       pot= pot + fixedpot(r)
       fr=fr-fixedforce(r)
! fixed force negative, fr positive       
! fixed halo is spherical so no fth       
      endif 

         fs = -(sintheta*fr + costheta/r*fth)
         fz = -(costheta*fr - sintheta/r*fth)
 
!         write(*,*) 'Force before disk correction ',fr,fth
         if( idiskflag .eq. 1 .OR. idiskflag.EQ.3) then
            call appdiskforce(s,z,fsad,fzad)
c            write(*,*) 'Force of disk correction ',fsad,fzad
            fs = fs + fsad
            fz = fz + fzad
         endif

         if( idiskflag .eq. 2 .OR. idiskflag.EQ.3) then
!            print*,fs,fz,thickdiskrforce(s,z)
            fs = fs + thickdiskrforce(s,z)
            fz = fz + thickdiskzforce(s,z)
         endif
      endif
      return
      end
