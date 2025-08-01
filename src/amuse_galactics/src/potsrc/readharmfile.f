c     this file contains the following subroutines:
c
c     readharmfile: reads a file with the harmonic expansion of the
c     potential and density.  all the relevant stuff ends up in a common
c     block that only needs to be seen by the routines in this file (I
c     think!). This routine must be called first.
c
c     dens(r,z): density at point (r,z)
c
c     densrpsi(r,psi): density as a function of r and potential psi.
c
c     pot(r,z): potential at point (r,z)
c
c     all the rest is (slightly nobbled in places) numerical recipes
c     routines.
c
      subroutine readharmfile(filename)

      include 'commonblocks.h'
      include 'equivalence.h'

      character*30 filename

      real psi0mpsid

C      open(file='in.gendenspsi',unit=1,status='old')
C      read(1,*) npsi
C      close(1)


c  constants for legendre polynomials
      do l=0,40
         plcon(l) = sqrt((2*l+1)/(4.0*pi))
      enddo
c
      open(14,file=filename,status='old')
      read(14,*)
      read(14,'(2x,7g15.5,i6,i4)') chalo,v0,a,nnn,v0bulge,abulge,
     +     dr,nr,lmax
      read(14,*)
      read(14,'(2x,3g15.8)') psi0, haloconst, bulgeconst
      read(14,*)
      read(14,'(2x,5g15.5)') rmdisk, rdisk, zdisk, outdisk, drtrunc
      read(14,*)
      read(14,'(2x,3g15.5)') psic,psi0mpsid,bhmass
      psid = psi0 - psi0mpsid
      read(14,'(2x,4i5)') idiskflag, ibulgeflag, ihaloflag, ibhflag
      read(14,*)
      do ir=0,nr
         read(14,'(8g16.8)') rdummy,(adens(l/2+1,ir),l=0,lmax,2)
         enddo
      read(14,*)
      read(14,*)
      do ir=0,nr
         read(14,'(8g16.8)') rdummy,(apot(l/2+1,ir),l=0,lmax,2)
         enddo
      read(14,*)
      read(14,*)
      do ir=0,nr
         read(14,'(8g16.8)') rdummy,(fr(l/2+1,ir),l=0,lmax,2)
         enddo
      close(14)

c
c these are needed for the bulgeconstants common block
c
      sigbulge2=sigbulge*sigbulge
c
c additional disk constants for efficiency - common diskconstants
c

      nrdisk= int((outdisk + 2.0*drtrunc)/dr) + 10
      rdisk2 = rdisk*rdisk
      diskconst = rmdisk/(4.0*pi*rdisk2*zdisk)
c
c and some more constants...
c
      v02=v0*v0
      v03=v0*v02
c
c calculate potential correction
c
      redge = nr*dr
      potcor = 0
      do l=0,lmax,2
          potcor = potcor - apot(l/2+1,nr) + fr(l/2+1,nr)*redge/(l+1)
      enddo
c      potcor = potcor*plcon(0)
c      potcor1 = potcor
c
c transfer gparams buffer for output
c

c      do i=1,15
c         ibuf1(i) = ibuf(i)
c      enddo
c      do j=1,3
c         jbuf1(j) = jbuf(j)
c      enddo
c      do k=1,1
c         kbuf1(k) = kbuf(k)
c      enddo
      
      return
      end

	
      subroutine fixpot

      include 'commonblocks.h'
      integer*4 i

      do i=0,nr
          apot(1,i) = apot(1,i) + potcor
      enddo

      return
      end
