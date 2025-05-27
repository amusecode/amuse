      subroutine halopotentialestimate

      include 'commonblocks'
      common /haloestimate/ hdens(0:nmax),hpot(0:nmax),hfr(0:nmax)

      eps = 0.0001
     
      do ir=0,nr
         r=ir*dr
         if(r.eq.0.) r = eps
         hdens(ir)=halodensity(r)
 77   enddo

c     now get the potential harmonics of this new density. (BT 2-208)
c     Simpson's rule integration.

      s1(0)=0
      r = 2*dr
      s1(2)=(r*dr/3.)*(4*hdens(1)*(1.0-dr/r)**2+hdens(2))
      rold = r 
      do ir=4,nr,2
         r=ir*dr
         s1a = (r*dr/3.)*(hdens(ir-2)*(1.0-2*dr/r)**2+
     +        4*hdens(ir-1)*(1.0-dr/r)**2+hdens(ir))
         s1(ir) = s1a + s1(ir-2)*rold/r
         rold = r
      enddo
      s2(nr)=0
      rold = nr*dr
      do ir=nr-2,2,-2
         r=ir*dr
         s2a = (r*dr/3.)*(hdens(ir+2)*(1.0+2*dr/r)+
     &        4*hdens(ir+1)*(1.0+dr/r)+hdens(ir))
         s2(ir) = s2a + s2(ir+2)
         rold = r
      enddo

      do ir = 2,nr,2
         r = ir*dr
         hpot(ir)=(4.*pi)*(s1(ir)+s2(ir))
         hfr(ir)=-(4.*pi)*s1(ir)/r
      enddo

      hpot(0)=3*(hpot(2)-hpot(4))+hpot(6)
      hfr(0)=0.0
      
c     then linearly interpolate other bins.
      
      do ir=1,nr-1,2
         hpot(ir)=(hpot(ir-1)+hpot(ir+1))/2.
         hfr(ir)=(hfr(ir-1)+hfr(ir+1))/2.
      enddo

      return
      end

      function haloforce(r)

      include 'commonblocks'
      common /haloestimate/ hdens(0:nmax),hpot(0:nmax),hfr(0:nmax)

      ihi=int(r/dr)+1
      if (ihi.lt.1) ihi=1
      if (ihi.gt.nr) ihi=nr
      r1=dr*(ihi-1)
      r2=dr*ihi
      t=(r-r1)/(r2-r1)
      tm1 = 1.0 - t
      if (r.eq.0.) then
         lmaxx=0
         costheta=0
      else
         costheta=z/r
         lmaxx=lmax
      endif
      haloforce = (t*hfr(ihi)+ tm1*hfr(ihi-1))

      return
      end

