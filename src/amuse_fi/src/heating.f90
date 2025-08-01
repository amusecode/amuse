module HeatingMod
 implicit none
 private
 public :: crheat
 
 real, parameter :: sece = 35

contains

 function crheat_f1(x) result(f1)
  real :: x,f1
  f1=.9971*(1-(1-x**.2663)**1.3163)
 end function
 
 function crheat_f2(e) result(f2)
  real :: e,f2
  f2=2*e-13.6
 end function 
 
 function crheat_f3(x) result(f3)
  real :: x,logx,f3
  logx=log10(x)
  f3=-.265*(-2.37+logx)*(4.69+logx)*(15.62+7.21*logx+logx**2)
 end function
 
 function crheat(x)
  real :: crheat,x,xe
  xe=x 
  if(x.lt.5.e-5) xe=5.e-5
  if(x.gt.1) xe=1
  crheat=crheat_f1(xe)+1./(crheat_f2(sece)/crheat_f3(xe)+1./(1-crheat_f1(xe))) 
  crheat=crheat*sece
 end function
 
end module 

!program test
!use heatingMod
!real xe

!10 read*, xe
!   if(xe.eq.0) stop
!   print*, crheat(xe)
!   goto 10
   
!end program   
