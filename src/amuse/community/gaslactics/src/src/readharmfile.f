c this file contains the following subroutines:
c
c readharmfile: reads a file with the harmonic expansion of the potential 
c               and density.
c all the relevant stuff ends up in a common block that only needs to be seen by
c the routines in this file (I think!). This routine must be called first.
c
c dens(r,z): density at point (r,z)
c
c densrpsi(r,psi): density as a function of r and potential psi.
c
c pot(r,z): potential at point (r,z)
c
c all the rest is (slightly nobbled in places) numerical recipes routines.
c

      subroutine readharmfile(filename,ibuf1)
      parameter(pi=3.1415926535)
      parameter(nrsize=20000)
      common /potconstants/ apot(20,0:nrsize), frad(20,0:nrsize), 
     +     dr, nr, lmax, potcor
      common /legendre/ plcon(0:40)
      common /gparameters/  a, b, c, v0, q, psi00, 
     +                      psiout, rho1, sigbulge2, 
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea
      common /diskpars/ sigr0, disksr, nrdisk
      common /flags/ idiskflag, ibulgeflag, ihaloflag

      real adens(20,0:nrsize)
      integer*4 ibuf(18),ibufb(3), ibuf1(18)
c this allows the gparameters to be passes as a C structure
      character filename*60
      equivalence (a, ibuf),(idiskflag,ibufb)
      external alldensnohalo,alldensnobulge,zerofunc

c  constants for legendre polynomials
      do l=0,40
         plcon(l) = sqrt((2*l+1)/(4.0*pi))
      enddo
c
      open(14,file=filename,status='old')
      read(14,*)
      read(14,'(2x,7g15.5,i6,i4)') psi00,v0,q,a,b,c,dr,nr,lmax
      read(14,*)
      read(14,'(2x,3g15.5)') psiout,rho1,sigbulge
      read(14,*)
      read(14,'(2x,5g15.5)') rmdisk, rdisk, zdisk, outdisk, drtrunc
      read(14,'(2x,3i5)') idiskflag, ibulgeflag, ihaloflag
      read(14,*)
!      write(0,*) ' this model had psi00,v0,q,a,b,c,dr,nr,lmax='
!      write(0,'(x,7g15.5,i6,i4)') psi00,v0,q,a,b,c,dr,nr,lmax

c these are needed for the bulgeconstants common block
c       
      if(ibulgeflag.NE.3) then
      sigbulge2=sigbulge*sigbulge
      bulgea = rho1*exp((psi00 - psiout)/sigbulge2)
      endif
c
c additional disk constants for efficiency - common diskconstants
c
      nrdisk= int((outdisk + 2.0*drtrunc)/dr) + 10
      rdisk2 = rdisk*rdisk
      diskconst = rmdisk/(4.0*pi*rdisk2*zdisk)
c
c and some more constants...
c
      if(ihaloflag.NE.3) then
      v02=v0*v0
      v03=v0*v02
      endif
! in case ihaloflag=3: v0=rscl,q=rv,a=rhs ...
! in case ibulgeflag=3: psiout=bulgerscl, rho1=bulgerv,sigbulge=bulgerhs....

      if(ihaloflag.eq.3) call sethalodf(v0,q,a)
 !     print*,'##',v0,q,a,ihaloflag
      if(ibulgeflag.eq.3)
     &      call setbulgedf(psiout,rho1,sigbulge,halormin(),halormax())
! ugly
!      print*,psiout,rho1,sigbulge,halormax()
!      print*,v0,q,a

      if(ihaloflag.eq.3) then
        call inithalodf(alldensnohalo)
      endif

      if(ibulgeflag.eq.3) then
        call initbulgedf(alldensnobulge)
      endif
      
      write(0,*) 'Reading harmonics file '//trim(filename)
      do ir=0,nr
         read(14,'(8g16.8)') rr,(adens(l/2+1,ir),l=0,lmax,2)
         enddo
      read(14,*)
      read(14,*)
      do ir=0,nr
         read(14,'(8g16.8)') rr,(apot(l/2+1,ir),l=0,lmax,2)
         enddo
      read(14,*)
      read(14,*)
      do ir=0,nr
         read(14,'(8g16.8)') rr,(frad(l/2+1,ir),l=0,lmax,2)
         enddo
      close(14)
c      write(0,*) 'have read density and potential'
c      write(*,*) 'lmax =',lmax
c      write(*,*) 'apot starts as:'
c      do ll=1,3
c         write(*,*) (apot(ll,i),i=1,10)
c      enddo
c
c
c
c calculate potential correction
c
!!!!
      redge = nr*dr
      potcor = 0
      do l=0,lmax,2
          potcor = potcor + apot(l/2+1,nr) + frad(l/2+1,nr)*redge/(l+1)
      enddo
      potcor = potcor*plcon(0)
      potcor1 = potcor
!!!
c
c transfer gparams buffer for output
c

      do i=1,15
          ibuf1(i) = ibuf(i)
      enddo
      do i=1,3
          ibuf1(i+15) = ibufb(i)
      enddo
      
      return
      end
