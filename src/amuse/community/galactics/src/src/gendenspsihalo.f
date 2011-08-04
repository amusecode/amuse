      subroutine gendenspsihalo

      include 'commonblocks'

      real psi

      open(file='denspsihalo.dat',unit=50,status='replace')

c      fcut_halo = dfhalo(psic)
      fcut_halo = 0.

      do i=1,npsi
         psi = tableE(i)
         call gethalodens(psi,denspsihalo(i))
         if(i.eq.npsi) denspsihalo(i) = 0.
         write(50,*) tableE(i),denspsihalo(i)
      enddo
      
      close(50)
      return
      
      end

      subroutine gethalodens(psi,rho)

      include 'commonblocks'

      real psi

      m = 200
      vmin = -10.
      vmax = log(sqrt(2.*(psi-psic)))
      dlogv = (vmax - vmin)/float(m-1)
      rho = 0.
      do j=1,4
         vlog = vmin + dlogv*float(j-1)
         v = exp(vlog)
         e = psi - v*v/2.
         df_lowered = dfhalo(e)-fcut_halo
         if(df_lowered.gt.0.) 
     +         rho = rho + 4.*pi*dlogv*v*v*v*coef(j)*df_lowered
      enddo
      do j=5,m-4
         vlog = vmin + dlogv*float(j-1)
         v = exp(vlog)
         e = psi - v*v/2.
         df_lowered = dfhalo(e)-fcut_halo
         if(df_lowered.gt.0.) 
     +         rho = rho + 4.*pi*dlogv*v*v*v*df_lowered
      enddo
      do jj=4,2,-1
         j = m-jj+1
         vlog = vmin + dlogv*float(j-1)
         v = exp(vlog)
         e = psi - v*v/2.
         df_lowered = dfhalo(e)-fcut_halo
         if(df_lowered.gt.0.) 
     +         rho = rho + 4.*pi*dlogv*v*v*v*coef(jj)*df_lowered
      enddo
      return
      end

      function coef(j)
      
      if(j.eq.1) coef=17./48.
      if(j.eq.2) coef=59./48.
      if(j.eq.3) coef=43./48.
      if(j.eq.4) coef=49./48.

      return
      end
