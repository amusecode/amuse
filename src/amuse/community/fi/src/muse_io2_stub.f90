subroutine outbods(fname)
 character(len=*) :: fname
end subroutine outbods


subroutine writebods(filenaam)
 character(len=*) :: filenaam 
end subroutine writebods

subroutine readbods(filenaam)
 character(len=*) :: filenaam 
end subroutine readbods


subroutine prepareoutput
end subroutine

subroutine itos(i,n,s)
 integer,intent(in) :: i,n
 character(len=n), intent(inout) :: s
 character(len=11) :: nstring
 data nstring/'0123456789X'/
 integer :: j,k,l

 if(i.LT.0.OR.i.GE.10**n) then
  do k=1,n
  s(k:k)=nstring(11:11)
  enddo
  return
 endif 
 j=1
 do k=1,n
 l=1+mod(i,10*j)/j
 s(n-k+1:n-k+1)=nstring(l:l)
 j=j*10
 enddo

end subroutine

! insert autmatically generated code below

subroutine writedump(n) 
 integer n 
end subroutine

subroutine readdump(n) 
 integer n
end subroutine
