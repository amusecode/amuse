program test

  use io_ramses

  integer::i,j,k
  integer,parameter::n=128*128
  real(kind=8),dimension(1:n)::x,y,z
  real(kind=8),dimension(1:n,1:6)::v
  real(kind=4),dimension(1:128,1:128)::map
  integer,dimension(1:n)::l
  character(LEN=128)::file

  do j=0,127
     do k=0,127
        i=1+j+128*k
        x(i)=0.50904+dble(j)/dble(128)*0.002-0.001
        y(i)=0.46650!+dble(k)/dble(128)*0.002-0.001
        z(i)=0.50152+dble(k)/dble(128)*0.002-0.001
     end do
  end do
  file='output_00699'
!  call getcell(x(1),y(1),z(1),v(1,1),l(1),n,6,file,verbose=.true.)
  call getcell(x,y,z,v,l,n,6,file,verbose=.true.)

  do j=0,127
     do k=0,127
        i=1+j+128*k
        map(j+1,k+1)=v(i,1)
     end do
  end do
  open(1,file='toto.dat',form='unformatted')
  i=128
  write(1)i,i
  write(1)map
  close(1)

end program test

