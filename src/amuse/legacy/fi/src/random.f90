! simple code for random number sequences
! from numerical recipes..
function ranfb(ix)
  real :: y,ranfb
  integer :: ix,ix1,ix2,ia1,ia2,ic,j1,j2,j3

  ia1=273
  ia2=25301
  ic=226908345
  
  if (ix.lt.0) ix=-ix
  if (ix.gt.2**30) ix=MOD(ix,2**30)

  ix1=ix/(2**15)
  ix2=ix-ix1*(2**15)
  j1=MOD(ix1*ia2,2**15)*2**15
  j2=MOD(ix2*ia1,2**15)*2**15
  j3=ix2*ia2
  ix=MOD(j1+j2,2**30)
  ix=MOD(ix+j3,2**30)
  ix=MOD(ix+ic,2**30)
  y=FLOAT(ix)

  ranfb=y/2**30

end function

subroutine setRND
  include 'globals.h'
  integer i
  real ranfb
  do i=1,nrndtable
    rndtable(i)=ranfb(rnseed)
  enddo
end subroutine

function pRND(i)
  include 'globals.h'
  real :: pRND
  integer :: i
  pRND=rndtable(mod(i,nrndtable)+1)
end function
