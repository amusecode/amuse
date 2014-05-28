c to fix: interpolation in rcirc does not seem to work well. after 
c that, hope that things come together!

      subroutine readdiskdf(filename1,ibuf)
      parameter(nrmax=1000)
      parameter (pi=3.1415926535)
      integer*4 ibuf(15)
      character*60 filename, filename1

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
      common /splines/ rr(0:nrmax), fdrat(0:nrmax), drat2(0:nrmax), 
     +                 fszrat(0:nrmax), szrat2(0:nrmax), nrspl

      call readharmfile(filename1,ibuf)
      filename = filename1
c Read in the correction functions
c
      write(0,*) 'Reading disk DF correction functions from cor'//filename

      open(17,file='cor'//filename,status='old')
      read(17,'(2x, 2g17.7,x,i4)') sigr0,disksr,nrspl
      do i=0,nrspl
         read(17,*) rr(i),fdrat(i),fszrat(i)
c         write(0,*) rr(i),fdrat(i),fszrat(i)
      enddo
      close(17)
      call splined(rr(0),fdrat(0),nrspl+1,1.e32,1.e32,drat2(0))
      call splined(rr(0),fszrat(0),nrspl+1,1.e32,1.e32,szrat2(0))

      rfid = 2.5*rdisk
      call omekap(rfid,fom,fka)
      sigr = sqrt(sigr2(rfid))
      sigden = diskdensf(rfid,0.0)*2.0*zdisk

c
c a close estimate of the disk surface density
c
      sigrcrit = 3.36*sigden/fka
      qtoomre = sigr/sigrcrit
      write(0,*) 'Toomre Q = ',qtoomre, ' at R = 2.5 R_d'

      open(12,file='toomre.dat',status='unknown')
      do i=0, nrspl
          call omekap(rr(i),fom,fka)
          sigr = sqrt(sigr2(rr(i)))
          sigden = diskdens(rr(i),0.0)*2.0*zdisk
          sigrcrit = 3.36*sigden/fka
          qtoomre = sigr/sigrcrit
          write(12,*) rr(i), qtoomre
      enddo
      close(12)

      return

      end
