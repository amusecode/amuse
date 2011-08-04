      subroutine gendenspsibulge

      include 'commonblocks'

      real psi

      open(file='denspsibulge.dat',unit=50,status='replace')

      fcut_bulge = dfbulge(psic)
      fcut_bulge = 0.

      do i=1,npsi
         psi = tableE(i)
         call getbulgedens(psi,denspsibulge(i))
         if(i.eq.npsi) denspsibulge(i) = 0.
         write(50,*) tableE(i),denspsibulge(i)
      enddo
      
      close(50)
      return
      
      end

      subroutine getbulgedens(psi,rho)

      include 'commonblocks'
      real psi

      m = 100
      v02 = v0*v0
      vmin = -10.
      vmax = log(sqrt(2.*(psi-psic)))
      dlogv = (vmax - vmin)/float(m-1)
      rho = 0.
      do j=1,4
         vlog = vmin + dlogv*float(j-1)
         v = exp(vlog)
         e = psi - v*v/2.
         df_lowered = dfbulge(e)-fcut_bulge
         if(df_lowered.gt.0.) 
     +         rho = rho + 4.*pi*dlogv*v*v*v*coef(j)*df_lowered
      enddo
      do j=5,m-4
         vlog = vmin + dlogv*float(j-1)
         v = exp(vlog)
         e = psi - v*v/2.
         df_lowered = dfbulge(e)-fcut_bulge
         if(df_lowered.gt.0.) 
     +         rho = rho + 4.*pi*dlogv*v*v*v*df_lowered
      enddo
      do jj=4,2,-1
         j = m-jj+1
         vlog = vmin + dlogv*float(j-1)
         v = exp(vlog)
         e = psi - v*v/2.
         df_lowered = dfbulge(e)-fcut_bulge
         if(df_lowered.gt.0.) 
     +         rho = rho + 4.*pi*dlogv*v*v*v*coef(jj)*df_lowered
      enddo
      return
      end

