      program haloinfo
     
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
      real r,z

      filename='dbh.dat'
      call readharmfile(filename,ibuf)
      
10    read*,r,z
      if(r.lt.0) stop
      print*,r, halodens(r,z)
      goto 10
      print*,' '
      
      end 
