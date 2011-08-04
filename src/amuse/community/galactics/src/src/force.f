      subroutine force(s,z,fs,fz,pot)

      parameter (pc0=0.282094792, pc2=0.630783131, pc4=0.846284375)
      parameter (pc6=1.017107236, pc8=1.163106623)
      real pc(20), p(20), dp(20)
      real pot

      include 'commonblocks'

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
               dp(l/2+1) = l*(plgndr1(l-1, costheta) - 
     +              costheta*p(l/2+1))/st2
            endif
         enddo
         do i=1,lmax/2+1
            p(i) = p(i)*pc(i)
            dp(i) = dp(i)*pc(i)
         enddo

         if( r .le. redge ) then
             frr = 0.0
             do i=1,lmax/2+1
                frr = frr + p(i)*(t*fr(i,ihi) + tm1*fr(i,ihim1))
             enddo
             fth = 0.0
             do i=2,lmax/2+1
                fth = fth - sintheta*dp(i)*(t*apot(i,ihi) + 
     +               tm1*apot(i,ihim1))
             enddo
             pot = 0.0
             do i=1,lmax/2+1
                pot = pot + p(i)*(t*apot(i,ihi) + tm1*apot(i,ihim1))
             enddo
         else
             frr = 0.0
             do i=1,lmax/2+1
                 l = 2*(i-1)
                 frr = frr-(l+1)*p(i)*apot(i,nr)/redge*(redge/r)**(l+2)
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
         if( idiskflag .eq. 1 ) then
             pot = pot + appdiskpot(s,z)
         endif

         fs = -(sintheta*frr + costheta/r*fth)
         fz = -(costheta*frr - sintheta/r*fth)
         if( idiskflag .eq. 1 ) then
            call appdiskforce(s,z,fsad,fzad)
            fs = fs + fsad
            fz = fz + fzad
         endif
      endif
      return
      end

