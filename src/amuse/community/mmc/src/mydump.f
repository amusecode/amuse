      subroutine mydump(i)
*
*
*       common save or read.
*       --------------------
*
      include 'params.h'
c      include 'zdata.h'
*
      integer i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i14,i15,i16,i17,
     &        i18,i20,i21,i22,i23,i24,i25,i26,i27,i28,i29,i30,i31,k1,
     &        i,nco,i32,i33,i34
*
      parameter (i1=32,i2=nlagra+10,i3=8*nmax,
     &           i4=8*nbmax3,i5=350*nbmax3,i6=57,i7=26,
     &           i8=18,i9=4,i10=120,i11=7*nmax+nbmax3+nzonma+2*nsupzo+3,
     &           i12=48+4*nmax,i14=2*nmax,i15=24*nmax+24,i16=3*nmax+3,
     &           i17=nmax+1,i18=3*nmax+3,i20=3,i21=2,i22=1,i23=1,
     &           i24=6,i25=3,i26=2,i27=2,i28=1,i29=34,i30=225,i31=5,
     &           i32=200,i33=200,i34=2*nmax,k1=5*nmax)
*
      real*8  y1,y2,y3,y4,y5,y6,y7,y8,y9,y14,y15,rtidkg,smto,
     &        ys1,ys2,ys3,ys4,ys5,ys6,ys7,ys8,ys9,ys10,ys11,z1,zini
*
      real*4 timeold
*      
      integer iy10,iy11,iy12,iy13,nto,iys1,iys2,iys3,iys4,iys5,
     &        iys6,iys7,iflagns,iflagbh,itime
*
*     monte carlo commons
*

      character*200 datadir

      common /param/  y1(i1)
      common /coefi/  y2(i2)
      common /body/   y3(i3)
      common /binar/  y4(i4)
      common /bininf/ y5(i5)
      common /system/ y6(i6)
      common /proba/  y8(i8)
      common /corep/  y9(i9),nco
      common /iparam/ iy10(i10)
      common /ibody/  iy11(i11)
      common /isyste/ iy12(i12)
      common /kingm/  rtidkg
      common /uptime/ y14(i14)
      common /oldpot/ smto,nto
      common /integral/ z1(k1)
      common /randx/  iy13(i29)
      common /zset/ zini
      common /fflags/ iflagns,iflagbh
      common /runtime/ timeold
      common /iruntime/ itime
*
*     stellar evolution commons
*
      common /value1/ ys1(i20)
      common /value2/ ys2(i21)
      common /value4/ ys3(i22),iys1(i23)
      common /value5/ ys4(i24)
      common /points/ ys5(i25)
      common /tstepc/ ys6(i26)
      common /params/ ys7(i27)
      common /stellar/ ys8(i15),iys6(i17)
      common /binary/ ys9(i16),iys7(i18)
      common /mscff/ ys10(i32)
      common /gbcff/ ys11(i33)
*
      common /value3/ iys2(i28)
      common /rand3/ iys3(i29)
      common /types/ iys4(i30,i30)
      common /flags/ iys5(i31)
      common /AMUSE/ datadir
*
*       open restart file   -  restart.fil 
*
c       open(1,file='restart.fil',status='unknown',form='formatted')
       open(1,file=trim(datadir)//'/restart.fil',status='unknown',
     &      form='unformatted')
*
*
*       read all common variables saved for restart
*
      if(i.eq.2) then
*
c         read (1,*)   (y1(k),k=1,i1)
c         read (1,*)   (y2(k),k=1,i2)
c         read (1,*)   (y3(k),k=1,i3)
c         read (1,*)   (y4(k),k=1,i4)
c         read (1,*)   (y5(k),k=1,i5)
c         read (1,*)   (y6(k),k=1,i6)
c         read (1,*)   (y7(k),k=1,i7)
c         read (1,*)   (y8(k),k=1,i8)
c         read (1,*)   (y9(k),k=1,i9)
c         read (1,*)   (iy10(k),k=1,i10),nco
c         read (1,*)   (iy11(k),k=1,i11)
c         read (1,*)   (iy12(k),k=1,i12)
c         read (1,*)   (iy13(k),k=1,34)
*
         read (1)   y1
         read (1)   y2
         read (1)   y3
         read (1)   y4
         read (1)   y5
         read (1)   y6
c         read (1)   y7
         read (1)   y8
         read (1)   y9
         read (1)   y14
c         read (1)   y15
         read (1)   ys1
         read (1)   ys2
         read (1)   ys3
         read (1)   ys4
         read (1)   ys5
         read (1)   ys6
         read (1)   ys7
         read (1)   ys8
         read (1)   ys9
         read (1)   ys10
         read (1)   ys11
         read (1)   rtidkg
         read (1)   smto,nto
         read (1)   z1
         read (1)   zini
         read (1)   timeold
         read (1)   iy10,nco
         read (1)   iy11
         read (1)   iy12
         read (1)   iy13
         read (1)   iys1
         read (1)   iys2
         read (1)   iys3
         read (1)   iys4
         read (1)   iys5
         read (1)   iys6
         read (1)   iys7
         read (1)   iflagns,iflagbh
         read (1)   itime
*
      endif
*
*       save all common variables needed for restart
      if(i.eq.1) then
*
c         write (1,*)   (y1(k),k=1,i1)
c         write (1,*)   (y2(k),k=1,i2)
c         write (1,*)   (y3(k),k=1,i3)
c         write (1,*)   (y4(k),k=1,i4)
c         write (1,*)   (y5(k),k=1,i5)
c         write (1,*)   (y6(k),k=1,i6)
c         write (1,*)   (y7(k),k=1,i7)
c         write (1,*)   (y8(k),k=1,i8)
c         write (1,*)   (y9(k),k=1,i9)
c         write (1,*)   (iy10(k),k=1,i10),nco
c         write (1,*)   (iy11(k),k=1,i11)
c         write (1,*)   (iy12(k),k=1,i12)
c         write (1,*)   (iy13(k),k=1,34)
*
         write (1)   y1
         write (1)   y2
         write (1)   y3
         write (1)   y4
         write (1)   y5
         write (1)   y6
c         write (1)   y7
         write (1)   y8
         write (1)   y9
         write (1)   y14
c         write (1)   y15
         write (1)   ys1
         write (1)   ys2
         write (1)   ys3
         write (1)   ys4
         write (1)   ys5
         write (1)   ys6
         write (1)   ys7
         write (1)   ys8
         write (1)   ys9
         write (1)   ys10
         write (1)   ys11
         write (1)   rtidkg
         write (1)   smto,nto
         write (1)   z1
         write (1)   zini
         write (1)   timeold
         write (1)   iy10,nco
         write (1)   iy11 
         write (1)   iy12
         write (1)   iy13
         write (1)   iys1
         write (1)   iys2
         write (1)   iys3
         write (1)   iys4
         write (1)   iys5
         write (1)   iys6
         write (1)   iys7
         write (1)   iflagns,iflagbh
         write (1)   itime
*
      endif
*
*
      close(1)
*
*
      return
*
      end
*
*
*
*
