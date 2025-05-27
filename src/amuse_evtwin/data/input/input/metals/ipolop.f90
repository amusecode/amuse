!> Program averopac:
!!
!! Interpolates opacity tables in the files phys.<z>.
!! Each file has 1 (header) + 1 (xfs) + nb * 9*127 (op) + 892 (dat) lines
!!             = 1 + 1 + nb*1143 + 892 lines  
!! The first two lines get adjusted for Z, the op blocks get interpolated and dat gets copied unchanged
!! Xfs contains one value of Xf for each block of op, hence max 10 blocks.
!! Xf = 2X + Y + 1 (?)
!!
!! In the latest (svn) version of the tables, the dat() table no longer exists and this part of the code should be adapted.
!<

program averopac
   implicit none
   integer, parameter :: double = kind(0.d0)
   integer, parameter :: nt=127, nr=90
   integer :: jr,jt
   
   real(double) :: op1(10,nr,nt),op2(10,nr,nt),op3(10,nr,nt)
   real(double) :: x1,x2,x3,y1,y2,y3,z1,z2,z3,feh1,feh2,feh3
   real(double) :: xfs1(10),xfs2(10),xfs3(10),dat(892,10)
   real(double) :: wt1,wt2
   real(double) :: minz,maxz
   real(double) :: rmin, rmax, tmin, tmax
   integer :: nb1,nb2,nb3,i,io
   character :: fname1*30,fname2*30,fname3*30,arg*30
   logical :: nodat
   
   !Table boundaries:
   rmin = -12.0
   rmax = 7.75
   tmin = 3.3
   tmax = 9.3
   
   
   if(iargc().eq.4) then
      call getarg(1,arg)
      fname1 = arg
      call getarg(2,arg)
      fname2 = arg
      call getarg(3,arg)
      fname3 = arg
      call getarg(4,arg)
      read(arg,*)z3
   else
      write(0,'(A)')'  '
      write(0,'(A)')'  Syntax:'
      write(0,'(A)')'    ipolop <infile1> <infile2> <outfile> <outZ>'
      write(0,'(A)')'    '
      write(0,'(A)')'    ipolop interpolates the opacity tables in the files infile1 and infile2'
      write(0,'(A)')'    and gererates the file outfile for metallicity outZ.  The metallicity'
      write(0,'(A)')'    outZ must be between those of the two input files.'
      write(0,'(A)')'    '
      stop
   end if
   
   
   !Open the input files:
   open (21, file=fname1, status='old', action='read')
   open (22, file=fname2, status='old', action='read')
   
   !Read header line:
   read(21,*)nb1,z1,x1
   read(22,*)nb2,z2,x2
   
   minz = min(z1,z2)
   maxz = max(z1,z2)
   
   if(z3.le.minz .or. z3.ge.maxz) then
      write(0,'(A)')'  '
      write(0,'(A,3(F9.6,A))')'  The output metallicity cannot be',z3,', but must lie between',minz,' and',maxz,'.'
      write(0,'(A)')'  '
      stop
   end if
   
   
   !Read XFs:
   read(21,'(1x,10f7.3)')xfs1
   read(22,'(1x,10f7.3)')xfs2
   
   
   !Check whether the files contain the same number of blocks
   if(nb1.ne.nb2) then
      write(0,'(A)')'  '
      write(0,'(A)')'  The files have an unequal number of data blocks:'
      write(0,'(2(2x,A,I4,3x))')trim(fname1)//':',nb1,trim(fname2)//':',nb2
      write(0,'(A)')'  '
      write(0,'(A)')'  The values for Xf in each file are:'
      write(0,'(1x,10f7.3)')xfs1+z1
      write(0,'(1x,10f7.3)')xfs2+z2
      write(0,'(A)')'  '
      write(0,'(A)')'  Check the two rows of Xf values above - there should be one for each block of opacities.'
      write(0,'(A)')'  In many cases, the first value of Xf differs, i.e. the first opacity block in the file.'
      write(0,'(A)')'  In such a case, you should remove the first value of Xf (and add an extra 0.000 at the end of the row)'
      write(0,'(A)')'  and the first block of 1143 lines of opacities, i.e., line 3-1145.'
      write(0,'(A)')'  '
      stop
   end if
   nb3 = nb1
   
   
   !Verify that the values for Xf are equal:
   do i=1,nb3
      if(xfs1(i).ne.1.0.and.xfs2(i).ne.1.0) then
         if(abs(xfs1(i)+z1-(xfs2(i)+z2)).gt.0.00001) then
            write(0,'(A)')'  '
            write(0,'(A,I4,A)')'   File block',i,' has unequal values for Xf'
            write(0,'(I4,A)')nb3,' blocks'
            write(0,'(1x,10f17.13)')xfs1+z1
            write(0,'(1x,10f17.13)')xfs2+z2
            write(0,'(A)')'  '
            stop
         end if
      end if
   end do
   
   
   feh1 = log10(z1/0.02)
   feh2 = log10(z2/0.02)
   
   
   !Read nb3 blocks of op(1:10,1:90,1:127), i.e., nb3 * 9 * 127 lines:
   do i=1,nb3
      read(21,'(1x,10f7.3)') ((op1(i,jr,jt), jr=1,nr), jt=1,nt)
      read(22,'(1x,10f7.3)') ((op2(i,jr,jt), jr=1,nr), jt=1,nt)
   end do
   
   
   !Read dat(892,10), i.e. 892 lines:
   read(21,'(1x,10f7.3)', iostat=io)dat
   nodat = .false.
   if(io.ne.0) then
      write(6,'(/,A)')'  The data block dat() seems to be absent in '//trim(fname1)
      write(6,'(A,/)')'    and will not be written to '//trim(fname3)//'.'
      nodat = .true.
   end if
   
   
   
   write(6,*)
   y1 = 1.d0-x1-z1
   y2 = 1.d0-x2-z2
   write(6,'(A14,F8.6,2(A6,F8.6),A11,F7.4,5x,A)')' File 1:  X = ',x1,'  Y = ',y1,'  Z = ',z1,'  [Fe/H] = ',feh1,trim(fname1)
   write(6,'(A14,F8.6,2(A6,F8.6),A11,F7.4,5x,A)')' File 2:  X = ',x2,'  Y = ',y2,'  Z = ',z2,'  [Fe/H] = ',feh2,trim(fname2)
   write(6,*)
   
   wt1 = abs((z3-z2)/(z2-z1))
   wt2 = 1.0 - wt1
   
   x3 = 0.76 - 3.*z3
   feh3 = log10(z3/0.02)
   
   write(6,'(A,2F10.4)')' The relative weights are ',wt1,wt2
   
   xfs3 = 0.
   do i=1,nb3
      xfs3(i) = xfs1(i)+z1-z3
      if(xfs1(i).eq.1.0.and.xfs2(i).eq.1.0) xfs3(i) = 1.0
      do jt=1,nt
         do jr=1,nr
            if (op1(i,jr,jt).le.-99.0 .or. op2(i,jr,jt).le.-99.0) then
               op3(i,jr,jt) = -99.999
            else
               op3(i,jr,jt) = log10(wt1*10.0**(op1(i,jr,jt)) + wt2*10.0**(op2(i,jr,jt)))
            end if
         end do
      end do
   end do !i
   
   
   write(6,*)
   y3 = 1.d0-x3-z3
   write(6,'(A14,F8.6,2(A6,F8.6),A11,F7.4,5x,A)')' File 3:  X = ',x3,'  Y = ',y3,'  Z = ',z3,'  [Fe/H] = ',feh3,trim(fname3)
   
   
   !Write output file:
   open (23,file=fname3, status='new', action='write')
   write(23,'(I4,2F9.5)')nb3,z3,x3  !Header
   write(23,'(1x,10f7.3)')xfs3      !XFs
   do i=1,nb3
      write (23,'(1x,10f7.3)') ((op3(i,jr,jt),jr=1,nr), jt=1,nt)  !op
   end do
   if(.not.nodat) write(23,'(1x,10f7.3)')dat       !dat
   close(23)
   
   write(6,*)
end program averopac

