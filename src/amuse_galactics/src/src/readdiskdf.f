      subroutine readdiskdf(filename1,ibuf)

      include 'commonblocks'

      integer*4 ibuf(15),jbuf(15),kbuf(15)
      character*30 filename, filename1

      common /fileroot/ filename

      filename1 = 'dbh.dat'

      call readharmfile(filename1,ibuf,jbuf,kbuf)
      filename = filename1

c
c Read in the correction functions
c
      open(17,file='cordbh.dat',status='old')
      read(17,'(2x, 2g17.7,1x,i4)') sigr0,disksr,nrspl
      do i=0,nrspl
         read(17,*) rr(i),fdrat(i),fszrat(i)
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

C      open(12,file='toomre.dat',status='unknown')
      ntoomre = int((outdisk + 2*drtrunc)/dr)
      xlindbladmax = 0.
      x50 = 0.
      do i=0, ntoomre
         rtoomre = float(i)*dr
         call omekap(rtoomre,fom,fka)
         sigr = sqrt(sigr2(rtoomre))
         sigden = diskdens(rtoomre,0.0)*2.0*zdisk
         sigrcrit = 3.36*sigden/fka
         qtoomre = sigr/sigrcrit
         xparam = fka*fka*rtoomre/6.283/sigden
         rcrise = fka*fka/2./fom/fom - 1.
         xlindblad = fom - fka/2.
         if(xlindblad.gt.xlindbladmax.and.i.gt.1) then
            xlindbladmax = xlindblad
            rxlindblad = rtoomre
         endif
C         write(12,666) rtoomre, qtoomre, xparam, fom,fka,rcrise
      enddo
C      close(12)
 666  format(6f16.3)
C      open(file='stability2.out',unit=10,status='replace')
C      write(10,*) xlindbladmax,rxlindblad,x50
C      close(10)

      return

      end
