program testgas
       integer*4 ibuf(18),ibufb(3), ibuf1(18)
       character filename*60

      common /gparameters/  a, b, c, v0, q, psi00, & 
     &                      psiout, rho1, sigbulge2, &
     &                      rmdisk, rdisk, zdisk, outdisk, drtrunc, &
     &                      potcor1


   integer, parameter :: n=1000,m=50
   real :: dens(n),z(n)
   real :: r,cs,gamma,mz,norm,mz2,rs(m),zm(m),z2m(m)

   filename='dbh.dat'
   call readharmfile(filename,ibuf)

  print*,pot(0.,0.),zdisk

  do j=1,m
  r=j/(1.*m)*16
  cs=10./65.6
  gamma=5./3.
  gamma=1.
  psi0=pot(r,0.)
  norm=0.
  mz2=0.
  mz=0.
  do i=1,n
   z(i)=(i-1)/(1.*n)*20*zdisk
   dens(i)=exp(-gamma*(pot(r,z(i))-psi0)/cs**2)
   mz=mz+dens(i)*z(i)
   mz2=mz2+dens(i)*z(i)**2
   norm=norm+dens(i)
  enddo
  mz=mz/norm
  mz2=mz2/norm
  rs(j)=r
  zm(j)=mz
  z2m(j)=mz2
  enddo 
 
  open(unit=1,file='testgas.data',status='UNKNOWN')
  do i=1,m
   write(1,*) rs(i),zm(i),z2m(i)
  enddo
  close(1)
 

 end program
