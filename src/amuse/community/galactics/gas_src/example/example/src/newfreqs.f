      program getfreqs

      parameter(pi=3.1415926535)
      common /potconstants/ apot(20,0:20000), fr(20,0:20000), 
     +     dr, nr, lmax, potcor
      common/legendre/ plcon(0:40)
      common /gparameters/  a, b, c, v0, q, psi00, 
     +                      psiout, rho1, sigbulge2, 
     +                      rmdisk, rdisk, zdisk, outdisk, drtrunc,
     +                      potcor1
      common /flags/ idiskflag, ibulgeflag, ihaloflag
      common /moreconstants/ v02, v03, rdisk2, diskconst, bulgea
      
      real adens(20,0:20000),potmaj(0:20000),
     +     potmajtot(0:20000),potup(0:20000),surden(0:20000),
     +     potmin(0:20000),vcmaj(0:20000),vctot(0:20000),vertfreq(0:20000),vcbulge(0:20000),
     +     vertfreqbulge(0:20000),psi2(0:20000),fr2(20,0:20000)

c  constants for legendre polynomials
      do l=0,40
         plcon(l) = sqrt((2*l+1)/(4.0*pi))
      enddo

      open(14,file='dbh.dat',status='old')
      open(15,file='newfreqdbh.dat',status='unknown')

      read(14,*) 
      read(14,'(2x,7g15.5,i6,i4)') psi00,v0,q,a,b,c,dr,nr,lmax
      read(14,*)
      read(14,'(2x,3g15.5)') psiout,rho1,sigbulge
      read(14,*)
      read(14,'(2x,5g15.5)') rmdisk, rdisk, zdisk, outdisk, drtrunc
      read(14,'(2x,3i5)') idiskflag, ibulgeflag, ihaloflag
      read(14,*)
      write(0,'(x,''Harmonics up to l='',i4)') lmax

      rdisk2 = rdisk*rdisk
      diskconst = rmdisk/(4.0*pi*rdisk2*zdisk)


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

      do ir=0,nr
         r=ir*dr
         potmajtot(ir)=pot(r,0.)
         potmin(ir)=pot(0.,r)
      enddo

      do ir=0,nr
         read(14,*) rr,(fr(l/2+1,ir),l=0,lmax,2)
      enddo
      read(14,*)
      read(14,*)

      do ir=0,nr
         r=ir*dr
         vctot(ir)=cfharm(r,0.,fr)
      enddo

      do ir=0,nr
         read(14,*) rr,(fr2(l/2+1,ir),l=0,lmax,2)
      enddo
      read(14,*)

      do ir=0,nr
         r=ir*dr
         psi2(ir)=cfharm(r,0.,fr2)
      enddo

      do ir=0,nr
         read(14,'(2g16.8)') rr,surden(ir)
      enddo
      
      close(14)

      print*,'total non-fixed mass:',vctot(nr)*ir*ir*dr*dr

       do ir=0,nr
        r=ir*dr
        if(ihaloflag.eq.2) then
	 vctot(ir)=vctot(ir)-fixedforce(r)
         psi2(ir)=psi2(ir)+fixedf2(r)
        endif
	if(idiskflag.EQ.2.OR.idiskflag.EQ.3) then
         vctot(ir)=vctot(ir)-thickdiskforce(r)
	 psi2(ir)=psi2(ir)+thickdiskf2(r)
	endif     	
	vctot(ir)=sqrt(max(0.,r*vctot(ir)))
       enddo

      
      drsave=dr
      nrsave=nr

      write(15,'(''#'',16(3x,a10,3x))') 'RADIUS','VC_TOT','PSIMAJ_TOT','PSI2'
      write(15,'(''#'')')

      ir=0
      write(15,'(9g16.8)') 0.,vctot(ir),potmajtot(ir),psi2(ir)

      do ir=1,nr  
         r=ir*dr
      write(15,'(9g16.8)') r,vctot(ir),potmajtot(ir),psi2(ir)
      enddo
 99   close(15)
      end
