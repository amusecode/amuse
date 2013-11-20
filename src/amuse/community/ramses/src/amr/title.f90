!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character(LEN=5)::nchar

  character(LEN=1)::nchar1
  character(LEN=2)::nchar2
  character(LEN=3)::nchar3
  character(LEN=4)::nchar4
  character(LEN=5)::nchar5

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

