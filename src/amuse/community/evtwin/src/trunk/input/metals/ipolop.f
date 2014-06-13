      program averopac

      implicit real (a-h,l,m,o-z)

      parameter (nt=127, nr=90)
      real :: op1(10,nr,nt),op2(10,nr,nt),op3(10,nr,nt)
      real :: x1,x2,x3,z1,z2,z3,feh1,feh2,feh3
      real :: xf1(10),xf2(10),xf3(10),dat(892,10)
      real :: minz,maxz
      integer :: n1,n2,i
      character :: fname1*30,fname2*30,fname3*30

      data rmin, rmax, tmin, tmax /-12.0, 7.75, 3.3, 9.3/

      write (*,'(/,1x,a,$)') 'Filename 1 ? '
      read (*,'(a)') fname1
      open (21, file=fname1, status='old')
      write (*,'(1x,a,$)') 'Filename 2 ? '
      read (*,'(a)') fname2
      open (22, file=fname2, status='old')
      
      read(21,*)n1,z1,x1
      read(22,*)n2,z2,x2
      
      read (21,'(1x,10f7.3)')xf1
      read (22,'(1x,10f7.3)')xf2
      
      if(n1.ne.n2) then
        print*,'The files have an unequal number of data blocks'
	print*,fname1,n1,fname2,n2
        write(6,'(1x,10f7.3)')xf1+z1
        write(6,'(1x,10f7.3)')xf2+z2
	goto 9999
      endif
            
      do i=1,n1
        if(xf1(i).ne.1.0.and.xf2(i).ne.1.0) then
          if(abs(xf1(i)+z1-(xf2(i)+z2)).gt.0.00001) then
	    print*,' File block',i,'has unequal XF'
	    print*,n1,'blocks'
            write(6,'(1x,10f7.3)')xf1+z1
            write(6,'(1x,10f7.3)')xf2+z2
	    goto 9999
	  endif
	endif
      enddo
      
      minz = min(z1,z2)
      maxz = max(z1,z2)
      
      feh1 = log10(z1/0.02)
      feh2 = log10(z2/0.02)
      
      do i=1,n1
        read (21,'(1x,10f7.3)') ((op1(i,jr,jt), jr=1,nr), jt=1,nt)
        read (22,'(1x,10f7.3)') ((op2(i,jr,jt), jr=1,nr), jt=1,nt)
      enddo
      
      read (21,'(1x,10f7.3)')dat
      
      print*,''
      write(6,'(A14,F8.6,2(A6,F8.6),A11,F7.4)')' File 1:  X = ',x1,
     &'  Y = ',  1-x1-z1,'  Z = ',z1,'  [Fe/H] = ',feh1
      write(6,'(A14,F8.6,2(A6,F8.6),A11,F7.4)')' File 2:  X = ',x2,
     &'  Y = ',  1-x2-z2,'  Z = ',z2,'  [Fe/H] = ',feh2
100   print*,''
      write(6,'(A47,F8.6,A3,F8.6,A4,$)')' What value for Z do you want i
     &n the new file (',minz,' - ',maxz,') ? '
      read*,z3
      
      if(z3.lt.minz.or.z3.gt.maxz) then
        write(6,'(A32,F8.6,A5,F8.6)')' The value of Z MUST be between ',
     & minz,' and ',maxz
	goto 100
      endif
      
      wt1 = abs((z3-z2)/(z2-z1))
      wt2 = 1.0 - wt1
      
      x3 = 0.76 - 3.*z3
      feh3 = log10(z3/0.02)
      
      print*,' The relative weights are ',wt1,wt2
      
      do i=1,n1
	xf3(i) = xf1(i)+z1-z3
	if(xf1(i).eq.1.0.and.xf2(i).eq.1.0) xf3(i)=1.
        do jt=1,nt
          do jr=1,nr
            if (op1(i,jr,jt).le.-99.0 .or. op2(i,jr,jt).le.-99.0) then
              op3(i,jr,jt) = -99.999
            else
              op3(i,jr,jt) = log10(wt1*10.0**(op1(i,jr,jt)) + 
     :              wt2*10.0**(op2(i,jr,jt)))
            endif
          enddo
        enddo
      enddo !i
      
      print*,''
      write(6,'(A14,F8.6,2(A6,F8.6),A11,F7.4)')' File 3:  X = ',x3,
     &'  Y = ',  1-x3-z3,'  Z = ',z3,'  [Fe/H] = ',feh3
      write (*,'(1x,a,$)') 'Filename (out) ? '
      read (*,'(a)') fname3
      
      open (23,file=fname3, status='new')
      write(23,'(I4,2F9.5)')n1,z3,x3
      write(23,'(1x,10f7.3)')xf3
      do i=1,n1
        write (23,'(1x,10f7.3)') ((op3(i,jr,jt),jr=1,nr), jt=1,nt)
      enddo
      write(23,'(1x,10f7.3)')dat
      close(23)

 9999 print*,''
      end
