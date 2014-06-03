      program getfreqs

      parameter(pi=3.1415926535)
      common /potconstants/ apot(20,0:20000), fr(20,0:20000), 
     +     dr, nr, lmax, potcor
      common/legendre/ plcon(0:40)
      common /gparameters/  a, b, c, v0, q, psi00, 
     +                      psiout, rho1, sigbulge2, 
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1

      real adens(20,0:20000),potmaj(0:20000),
     +     potmajtot(0:20000),potup(0:20000),surden(0:20000),
     +     potmin(0:20000),vcmaj(0:20000),vctot(0:20000),vertfreq(0:20000),vcbulge(0:20000),
     +     vertfreqbulge(0:20000),psi2(0:20000)

c  constants for legendre polynomials
      do l=0,40
         plcon(l) = sqrt((2*l+1)/(4.0*pi))
      enddo

      open(14,file='dbh.dat',status='old')
      read(14,*) 
      read(14,'(2x,7g15.5,i6,i4)') psi00,v0,q,a,b,c,dr,nr,lmax
      read(14,*)
      read(14,'(2x,3g15.5)') psiout,rho1,sigbulge
      read(14,*)
      read(14,'(2x,5g15.5)') rmdisk, rdisk, zdisk, outdisk, drtrunc
      read(14,'(2x,3i5)') idiskflag, ibulgeflag, ihaloflag
      read(14,*)
      write(0,'(x,''Harmonics up to l='',i4)') lmax
c
c READ HARMONICS IN FILE FILE INTO APOT ARRAY
c THESE CAN THEN BE USED BY POT(R,Z) TO GIVE THE CORRESPONDING FUNCTION
c
      do ir=0,nr
         read(14,*) rr,(adens(l/2+1,ir),l=0,lmax,2)
      enddo
      read(14,*)
      read(14,*)
c
c READ IN HARMONICS OF THE TOTAL POTENTIAL, TABULATE POT ON MAJ & MIN AXES
c
      do ir=0,nr
         read(14,*) rr,(apot(l/2+1,ir),l=0,lmax,2)
      enddo
      read(14,*)
      read(14,*)

c      write(*,*) 'Start of total potential harmonics:'
c      do ir=0,5
c         write(*,*) ir*dr,(apot(l/2+1,ir),l=0,lmax,2)
c      enddo

      open(15,file='freqdbh.dat',status='unknown')
      do ir=0,nr
         r=ir*dr
         potmajtot(ir)=pot(r,0.)
         potmin(ir)=pot(0.,r)
      enddo
c
c READ IN HARMONICS OF RADIAL FORCE, TABULATE ROTATION CURVE
c
      do ir=0,nr
         read(14,*) rr,(apot(l/2+1,ir),l=0,lmax,2)
      enddo
      read(14,*)
      read(14,*)

c      write(*,*) 'Start of total radial force harmonics:'
c      do ir=0,5
c         write(*,*) ir*dr,(apot(l/2+1,ir),l=0,lmax,2)
c      enddo

      do ir=0,nr
         r=ir*dr
         vctot(ir)=sqrt(max(0.,r*pot(r,0.)))
      enddo
c
c READ IN HARMONICS OF 2ND DERIVATIVE OF PSI, TABULATE IN PSI2
c         
      do ir=0,nr
         read(14,*) rr,(apot(l/2+1,ir),l=0,lmax,2)
      enddo
      read(14,*)

c      write(*,*) 'Start of total potential 2nd deriv harmonics:'
c      do ir=0,5
c         write(*,*) ir*dr,(apot(l/2+1,ir),l=0,lmax,2)
c      enddo

      do ir=0,nr
         r=ir*dr
         psi2(ir)=pot(r,0.)
      enddo
c 
c READ IN SURFACE DENSITY OF THE DISK, TABULATE
c
      do ir=0,nr
         read(14,'(2g16.8)') rr,surden(ir)
         enddo
      close(14)


      drsave=dr
      nrsave=nr
c 
c NOW REPEAT FOR HALO-ONLY HARMONICS FILE
c

      open(14,file='h.dat',status='old')
      read(14,*) 
      read(14,'(2x,7g15.5,i6,i4)') psi00,v0,q,a,b,c,drsave,nrsave,lmax
      read(14,*)
      read(14,'(2x,3g15.5)') psiout,rho1,sigbulge
      read(14,*)
      read(14,'(2x,5g15.5)') rmdisk, rdisk, zdisk, outdisk, drtrunc
      read(14,'(2x,3i5)') idiskflag, ibulgeflag, ihaloflag
      read(14,*)
      write(0,'(x,''Harmonics up to l='',i4)') lmax
      if (nr.ne.nrsave .or. dr.ne.drsave) then
         write(0,*) 'radial bins do not match.'
         stop
      endif

      idiskflag=0

c
c READ HALO POTENTIAL, TABULATE IT ON MAJOR & INCLINED AXES
c

      do ir=0,nr
         read(14,*) rr,(adens(l/2+1,ir),l=0,lmax,2)
      enddo
      read(14,*)
      read(14,*)

      do ir=0,nr
         read(14,*) rr,(apot(l/2+1,ir),l=0,lmax,2)
      enddo
      read(14,*)
      read(14,*)

c      write(*,*) 'Start of halo potential harmonics:'
c      do ir=0,5
c         write(*,*) ir*dr,(apot(l/2+1,ir),l=0,lmax,2)
c      enddo

      do ir=0,nr
         r=ir*dr
         potmaj(ir)=pot(r,0.)
         potup(ir)=pot(r*cos(0.05),r*sin(0.05))
      enddo

c
c READ HALO RADIAL FORCE, TABULATE ITS ROTATION CURVE AND VERT FREQ
c
      do ir=0,nr
         read(14,*) rr,(apot(l/2+1,ir),l=0,lmax,2)
      enddo

c      write(*,*) 'Start of halo radial force harmonics:'
c      do ir=0,5
c         write(*,*) ir*dr,(apot(l/2+1,ir),l=0,lmax,2)
c      enddo

      do ir=0,nr
         r=ir*dr
         if (ir.eq.0) then
            vcmaj(ir)=0
         else
            vcmaj(ir)=sqrt(max(0.,r*pot(r,0.)))
         endif
         if (ir.eq.0) then
            vertfreq(ir)=0
         else
            vertfreq(ir)=sqrt(max(0.,vcmaj(ir)**2+2.*(potup(ir)-potmaj(ir))/0.05**2))/r
         endif
      enddo

      close(14)

c REPEAT FOR BULGE


      open(14,file='b.dat',status='old')
      read(14,*) 
      read(14,'(2x,7g15.5,i6,i4)') psi00,v0,q,a,b,c,drsave,nrsave,lmax
      read(14,*)
      read(14,'(2x,3g15.5)') psiout,rho1,sigbulge
      read(14,*)
      read(14,'(2x,5g15.5)') rmdisk, rdisk, zdisk, outdisk, drtrunc
      read(14,'(2x,3i5)') idiskflag, ibulgeflag, ihaloflag
      read(14,*)
      write(0,'(x,''Harmonics up to l='',i4)') lmax
      if (nr.ne.nrsave .or. dr.ne.drsave) then
         write(0,*) 'radial bins do not match.'
         stop
      endif

      idiskflag=0
c 
c READ BULGE POTENTIAL, TABULATE
c
      do ir=0,nr
         read(14,*) rr,(adens(l/2+1,ir),l=0,lmax,2)
      enddo
      read(14,*)
      read(14,*)
      do ir=0,nr
         read(14,*) rr,(apot(l/2+1,ir),l=0,lmax,2)
      enddo
      read(14,*)
      read(14,*)

c      write(*,*) 'Start of bulge potential harmonics:'
c      do ir=0,5
c         write(*,*) ir*dr,(apot(l/2+1,ir),l=0,lmax,2)
c      enddo

      do ir=0,nr
         r=ir*dr
         potmaj(ir)=pot(r,0.)
         potup(ir)=pot(r*cos(0.05),r*sin(0.05))
      enddo
C
C READ BULGE RADIAL FORCE, MAKE ROT CURVE AND VERT FREQS.
C

      do ir=0,nr
         read(14,*) rr,(apot(l/2+1,ir),l=0,lmax,2)
      enddo

c      write(*,*) 'Start of bulge radial force harmonics:'
c      do ir=0,5
c         write(*,*) ir*dr,(apot(l/2+1,ir),l=0,lmax,2)
c      enddo

      do ir=0,nr
         r=ir*dr
         if (ir.eq.0) then
            vcbulge(ir)=0
         else
            vcbulge(ir)=sqrt(max(0.,r*pot(r,0.)))
         endif
         if (ir.eq.0) then
            vertfreqbulge(ir)=0
         else
            vertfreqbulge(ir)=sqrt(vcbulge(ir)**2+2.*(potup(ir)-potmaj(ir))/0.05**2)/r
         endif
      enddo

      close(14)

      write(15,'(''#'',16(3x,a10,3x))') 'RADIUS','OMEGA_H','NU_H',
     :     'SIGMA_D','VC_TOT','VC_B','NU_B','PSIMAJ_TOT','PSI'''''
      write(15,'(''#'')')

      ir=0
      write(15,'(9g16.8)') 0.,vcmaj(1)/dr,(4*vertfreq(1)-vertfreq(2))/3.,
     +     surden(ir),
     +     vctot(ir),vcbulge(ir),(4*vertfreqbulge(1)-vertfreqbulge(2))/3.,
     +     potmajtot(ir),
     +     psi2(ir)

      do ir=1,nr  
         r=ir*dr
         write(15,'(9g16.8)') r,vcmaj(ir)/r,vertfreq(ir),surden(ir),
     +           vctot(ir),vcbulge(ir),vertfreqbulge(ir),potmajtot(ir),
     +           psi2(ir)
         enddo
 99   close(15)
      end
