
      function densrpsi(r,psi)
      parameter (pi=3.1415926535, sq8=2.828427125)

      common /potconstants/ apot(20,0:20000), fr(20,0:20000), 
     +     dr, nr, lmax, potcor
      common /gparameters/  a, b, c, v0, q, psi00, 
     +                      psiout, rho1, sigbulge2, 
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea

c Lowered Evans density      
      if (psi .ge. 0.0) then
         densrpsi=0
      else
         t2=sqrt(-2.*psi)
         r2 = r**2
         tpsiv02 = -2.0*psi/v02
         if( tpsiv02 .gt. 35 ) then
             e2psi = 0
         else
             e2psi = exp(-2.*psi/v02)
         endif
         densrpsi=2*pi*(
     +         (a*r2*v02*.03125 + 0.125*b)*v03
     +        *erfcen(2.0*t2/v0)*e2psi*e2psi
     +        +(c-a*r2*v02/2-b)*v03/sq8
     +        *erfcen(2.0*sqrt(-psi)/v0)*e2psi
     +        +a*v02*r2*(0.375*v02-psi/3.)*t2+b*0.5*v02*t2
     +        +c*(4./3.*psi-v02)*t2)
c         if (densrpsi.lt.0.) write(*,*) 'negative densrpsi at r,psi=',r,psi,densrpsi
         endif
      return
      end
