      function rcirc(am)

      include 'commonblocks'

      real om

      dimension rtab(nmax),amtab(nmax),rtab2(nmax)
     +     ,rtab2zero(nmax)
      data ifirst /0/

      save ifirst,rtab,amtab,rtab2,rtab2zero
      
c make a spline fit to [Rcirc/sqrt(Lcirc)] vs. Lcirc.
c this is ~ constant at the origin, and goes like Lcirc**1.5 
c at large radii (where the potential is Kepler).

      if (ifirst.eq.0) then

         do i=1,nr
            r=(i-1)*dr
            call omekap(r,om,t)
            amtab(i)=om*r*r
            rtab(i)=1./sqrt(om)
         enddo
         slopeinf=1.5*rtab(nrdisk)/amtab(nr)
         call splined(amtab,rtab,nr,1.e32,slopeinf,rtab2)
         do i=1,nr
            rtab2zero(i) = 0.
         enddo         
         ifirst=1
      endif
      aam=abs(am)
      if (aam.gt.amtab(nr)) then
         rcirc=rtab(nr)*(aam/amtab(nr))**2
      else
         call splintd(amtab,rtab,rtab2zero,nr,aam,rc)
         rcirc=rc*sqrt(aam)
         if(rc.lt.0.) then
            write(0,*) rc
            stop
         endif
      endif

      return
      end
