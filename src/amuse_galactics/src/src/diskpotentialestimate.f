      subroutine diskpotentialestimate

      include 'commonblocks'
      common /diskestimate/ ddens(0:nmax),dpot(0:nmax),dfr(0:nmax)

      eps = dr/10.
      ntheta = 100
     
      do ir=0,nr
         r=ir*dr
         if(r.eq.0.) r = eps
         cthetamax = min(1.,10.*zdisk/r)
         dctheta=cthetamax/float(ntheta)
         s=s+dpolardiskdens(r,cthetamax)+dpolardiskdens(r,0.0)
         do is=1,ntheta-1,2
            ctheta=is*dctheta
            s=s+4*dpolardiskdens(r,ctheta)
         enddo
         do is=2,ntheta-2,2
            ctheta=is*dctheta
            s=s+2*dpolardiskdens(r,ctheta)
         enddo
         s=s*dctheta/3.
         ddens(ir)=s
 77   enddo

c     now get the potential harmonics of this new density. (BT 2-208)
c     Simpson's rule integration.

      s1(0)=0
      r = 2*dr
      s1(2)=(r*dr/3.)*(4*ddens(1)*(1.0-dr/r)**2+ddens(2))
      rold = r 
      do ir=4,nr,2
         r=ir*dr
         s1a = (r*dr/3.)*(ddens(ir-2)*(1.0-2*dr/r)**2+
     +        4*ddens(ir-1)*(1.0-dr/r)**2+ddens(ir))
         s1(ir) = s1a + s1(ir-2)*rold/r
         rold = r
      enddo
      s2(nr)=0
      rold = nr*dr
      do ir=nr-2,2,-2
         r=ir*dr
         s2a = (r*dr/3.)*(ddens(ir+2)*(1.0+2*dr/r)+
     &        4*ddens(ir+1)*(1.0+dr/r)+ddens(ir))
         s2(ir) = s2a + s2(ir+2)
         rold = r
      enddo

      do ir = 2,nr,2
         r = ir*dr
         dpot(ir)=(4.*pi)*(s1(ir)+s2(ir))
         dfr(ir)=-(4.*pi)*s1(ir)/r
      enddo

      dpot(0)=3*(dpot(2)-dpot(4))+dpot(6)
      dfr(0)=0.0
      
c     then linearly interpolate other bins.
      
      do ir=1,nr-1,2
         dpot(ir)=(dpot(ir-1)+dpot(ir+1))/2.
         dfr(ir)=(dfr(ir-1)+dfr(ir+1))/2.
      enddo

      return
      end

      function diskforce(r)

      include 'commonblocks'
      common /diskestimate/ ddens(0:nmax),dpot(0:nmax),dfr(0:nmax)

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
      diskforce = (t*dfr(ihi)+ tm1*dfr(ihi-1))

      return
      end

      function diskdensity(r)

      include 'commonblocks'
      common /diskestimate/ ddens(0:nmax),dpot(0:nmax),dfr(0:nmax)

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
      diskdensity = (t*ddens(ihi)+ tm1*ddens(ihi-1))

      return
      end
