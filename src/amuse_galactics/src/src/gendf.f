      subroutine gentableE

      include 'commonblocks'

      r = 0.
      rmin = dr/100.
      psi0 = gettotalpsi(r)
      psic = gettotalpsi(chalo + 5.*drtrunchalo)
      psid = gettotalpsi(rmin)
      tableE(1) = psid
      tableE(npsi) = psic
      do i=1,npsi
         xx = (psi0-psic)/(psi0-psid)
         yy = float(i-1)/float(npsi-1)
         tableE(i) = psi0 - (psi0-psid)*(xx**yy)
      enddo

      return
      end

      subroutine gendfnfw

      include 'commonblocks'

      open(file='dfnfw.dat',unit=50,status='replace')

      do i=1,npsi
         energy = tableE(i)
         dfnfw(i) = getdfnfw(energy)
         if(dfnfw(i).gt.0.) then 
            dfnfw(i) = log(dfnfw(i))
            dfnfwlast = dfnfw(i)
         else
            dfnfw(i) = dfnfwlast
         endif
         write(50,*) tableE(i),dfnfw(i)
      enddo
      close(50)

      return
     
      end

      function getdfnfw(energy)

      include 'commonblocks'
      real psi
      
      external gettotalpsi

      tmax = sqrt(energy-psic)
      dt = tmax/float(nint-1)
      psi = energy

      racc = 1e-8
      call findbrackets(gettotalpsi,psi,rmin,rmax)
      call rtbis(gettotalpsi,psi,rmin,rmax,racc,rpsi)
      d2rhodpsi2 =  getd2rhonfwdpsi2(rpsi)
      t = 0.
      sum = dt*d2rhodpsi2

      do i=1,nint-2
         t = dt*float(i)
         psi = energy - t*t
         call findbrackets(gettotalpsi,psi,rmin,rmax)
         call rtbis(gettotalpsi,psi,rmin,rmax,racc,rpsi)
         d2rhodpsi2 =  getd2rhonfwdpsi2(rpsi)
         sum = sum + 2.*dt*d2rhodpsi2
      enddo
      getdfnfw = sum/sqrt(8.)/pi/pi

      return
      end

      subroutine gendfsersic

      include 'commonblocks'

      open(file='dfsersic.dat',unit=50,status='replace')

      do i=1,npsi
         energy = tableE(i)
         dfsersic(i) = getdfsersic(energy)
         if(dfsersic(i).gt.0.) then
            dfsersic(i) = log(dfsersic(i))
            dfsersiclast = dfsersic(i)
         else
            dfsersic(i) = dfsersiclast
         endif
         write(50,*) tableE(i),dfsersic(i)
      enddo
      close(50)
      return
     
      end

      function getdfsersic(energy)

      include 'commonblocks'

      external gettotalpsi
      real psi

      tmax = sqrt(energy-psic)
      dt = tmax/float(nint-1)
      psi = energy

      racc = 1e-8
      call findbrackets(gettotalpsi,psi,rmin,rmax)
      call rtbis(gettotalpsi,psi,rmin,rmax,racc,rpsi)
      d2rhodpsi2 =  getd2rhosersicdpsi2(rpsi)
      sum = dt*d2rhodpsi2
      do i=1,nint-2
         t = dt*float(i)
         psi = energy - t*t
         call findbrackets(gettotalpsi,psi,rmin,rmax)
         call rtbis(gettotalpsi,psi,rmin,rmax,racc,rpsi)
         d2rhodpsi2 =  getd2rhosersicdpsi2(rpsi)
         sum = sum + 2.*dt*d2rhodpsi2
      enddo
      getdfsersic = sum/sqrt(8.)/pi/pi
      return
      end

      function getd2rhosersicdpsi2(rad)

      include 'commonblocks'

      real haloforce, halodens,L2

      external rtbiss,gammq,gammln

      Den = sersicdens(rad)
      Denp = sersicdensprime(rad)
      Denpp = sersicdens2prime(rad)
      Force = sersicforce(rad)

      totalden = Den
      if(ihaloflag.eq.1) then
         Force = Force + haloforce(rad)
         totalden = totalden + halodensity(rad)
      endif
      if(idiskflag.eq.1) then
         Force = Force + diskforce(rad)
         totalden = totalden + diskdensity(rad)
      endif

      bbb = 4.*pi*totalden*Denp/Force
      ccc = 2.*Denp/rad
      ddd = Denpp
      getd2rhosersicdpsi2 = (bbb + ccc + ddd)/(Force**2.)

      return
      END

      function getd2rhonfwdpsi2(rad)

      include 'commonblocks'

      real haloforce, halodensity, halodensprime, halodens2prime, L2
      
      external rtbiss,gammq,gammln

      s = rad/a

      Den = halodensity(rad)
      Denp = halodensprime(rad)
      Denpp = halodens2prime(rad)
      Force = haloforce(rad)

      totalden = Den

      if(ibulgeflag.eq.1) then

         Force = Force + sersicforce(rad)
         totalden = totalden + sersicdens(rad)

      endif

      if(idiskflag.eq.1) then
         Force = Force + diskforce(rad)
         totalden = totalden + diskdensity(rad)
      endif

      bbb = 4.*pi*totalden*Denp/Force
      ccc = 2.*Denp/rad
      ddd = Denpp
      getd2rhonfwdpsi2 = (bbb + ccc + ddd)/(Force**2.)

c      write(77,*) log10(rad),log10(Den),log10(-Denp),(Denpp),
c     +     getd2rhonfwdpsi2
         
      return
      END

