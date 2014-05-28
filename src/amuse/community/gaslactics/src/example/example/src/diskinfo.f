      subroutine scales(r1,r2,s1,s2)
       real r1,r2,h1,h2
      character*60 filename
      parameter (pi=3.1415926535)
      common /fileroot/ filename
      common /potconstants/ apot(20,0:20000), frad(20,0:20000), 
     +     dr, nr, lmax, potcor
      common /legendre/ plcon(0:40)
      common /gparameters/  a, b, c, v0, q, psi00, 
     +                      psiout, rho1, sigbulge2, 
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea
      common /diskpars/ sigr0, disksr, nrdisk
      dimension ibuf(100)

      filename='dbh.dat'
      call readharmfile(filename,ibuf)
      
      r1=0.
      s1=rmdisk/(2*pi*rdisk**2)
      
      r2=outdisk
      s2=s1*exp(-outdisk/rdisk)
      
      end 
