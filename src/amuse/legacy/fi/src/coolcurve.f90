program test
 use ElementsMod 
 use H2CoolMod
 use CoolingMod
 use IonMod
 character(len=200) :: indatadir
 real :: temp,n,h2,le,lp,cr,rtemp,xe,h2f
 
 call InitElements
 indatadir="/home/pelupes/fish/stars/data/"
 call InitIon(1.e9,1.)
 call InitH2(indatadir)
 call InitCooling(indatadir)
!10 read*, temp,n
! if(temp.eq.0) stop
! h2=H2CoolFunc(temp)
! lp=CoolFunc(1.,temp)
! write(*,5),temp, n, h2,lp
!goto 10 
 
 do i=1,100
  temp=10**(1+7*(i-1.)/99.)
  n=1.
  cr=18
  h2f=1.
  call ion(temp,rtemp,xe,n,cr)
  le=xe*xeCool(rtemp)
  lp=max(0.,(1-xe))*HCool(rtemp)
  h2=H2Cool(rtemp,n)
  write(*,5),rtemp,xe,le,lp,h2
 enddo

5 format( e12.4,' ',e12.4,' ',e12.4,' ',e12.4,' ',e12.4)
 end program test

