program amrdir
  !---------------------------------------------------------------------
  ! This code create empty output and backup directories for RAMSES runs.
  ! f90 amrdir.f90 -o ~/bin/amrdir
  !---------------------------------------------------------------------
  implicit none
  integer::n1,n2,n
  logical::ok
  character*5::nchar
  character*50::filedir,filecmd
  integer::iargc,ifout
  character(len=128)::arg

  n = iargc()
  if (n.NE.2) then
     print *, 'usage: amrdir n1 n2'
     stop
  end if
  call getarg(1,arg)
  read(arg,'(I8)')n1
  call getarg(2,arg)
  read(arg,'(I8)')n2

  do ifout=n1,n2
     call title(ifout,nchar)
     filedir='output_'//TRIM(nchar)
     filecmd='mkdir -p '//TRIM(filedir)
     call system(filecmd)
     filedir='backup_'//TRIM(nchar)
     filecmd='mkdir -p '//TRIM(filedir)
     call system(filecmd)
  end do

end program amrdir

!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character*5::nchar

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif


end subroutine title

