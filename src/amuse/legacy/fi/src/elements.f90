module ElementsMod
 implicit none
 private
 public :: elementnrQ,abundancesQ,weightsQ,fractionsQ,wfractionsQ,zQ ,&
           meanmwtQ,nelementsQ,elementnamesQ,wfractionQ,InitElements, &
           solarzQ,fractionQ,abundanceQ
 
  integer, parameter :: nelements=8
  character(len=4),dimension(nelements),parameter :: elementnames = &
             (/"H   ","He  ","C   ","N   ","O   ","Ne  ","Si  ","Fe  "/) 
  integer, parameter :: H=1,He=2,metals(6)=(/ 3,4,5,6,7,8 /)
  
  real :: abundances(nelements)=(/0., 0., 0., 0., 0., 0., 0., 0. /)

!  real, parameter :: abundances(nelements)= &
!                      (/ 1.,.09,3.75e-5,8.7e-6,4.4e-5,2.6e-6,3.2e-6,3.2e-6 /)
!  real, parameter :: abundances(nelements)= &
!                      (/ 1.,.1,3.75e-4,8.7e-5,4.4e-4,2.6e-5,3.2e-5,3.2e-5 /)

  real, parameter :: weights(nelements) = &
                      (/ 1.0079,4.0026,12.011,14.007,15.999,20.18,28.086,55.845/)
  real, save :: fractions(nelements)=0., wfractions(nelements)=0., meanmwt=0.
 
  real, parameter :: solarz=0.01545
 
contains
 
  function elementnrQ(string) result(x)
   character(len=4), intent(in) :: string
   integer :: i,x
   x=0
   do i=1,nelements
    if(elementnames(i).EQ.string) x=i
   enddo
  end function
  
  subroutine InitElements(metallicity)
   character(len=*), optional :: metallicity
   
   if(any(abundances.NE.0)) return
   
   if(present(metallicity)) then
    select case (metallicity)
    case('2.5-solar')
      abundances=(/ 1.,.085,8.28e-4,2.1e-4,1.67e-3,3.e-4,8.88e-5,8.e-5 /)
    case('solar')
      abundances=(/ 1.,.085,3.31e-4,8.3e-5,6.76e-4,1.2e-4,3.55e-5,3.2e-5 /)
    case('0.2-solar')
      abundances=(/ 1.,.085,6.62e-5,1.66e-5,1.33e-4,2.4e-5,7.1e-6,6.4e-6 /)
    case('0.1-solar')
      abundances=(/ 1.,.085,3.31e-5,8.3e-6,6.76e-5,1.2e-5,3.55e-6,3.2e-6 /)
    case('0.02-solar')
      abundances=(/ 1.,.085,6.62e-6,1.66e-6,1.33e-5,2.4e-6,7.1e-7,6.4e-7 /)
    case('0.01-solar')
      abundances=(/ 1.,.065,3.31e-6,8.3e-7,6.76e-6,1.2e-6,3.55e-7,3.2e-7 /)
    case('FeSi-depl')
      abundances=(/ 1.,.085,3.31e-4,8.3e-5,6.76e-4,1.2e-4,1.18e-6,1.07e-6 /)    
    case('C-depl')    
      abundances=(/ 1.,.085,1.1e-5,8.3e-5,6.76e-4,1.2e-4,3.55e-5,3.2e-5 /)
    case default
      print*,'unknown metallicity (compilation error)'         
      stop
    end select
   else
    abundances=(/ 1.,.085,3.31e-4,8.3e-5,6.76e-4,1.2e-4,3.55e-5,3.2e-5 /)
   endif
      
   fractions=abundances/sum(abundances)
   meanmwt=sum(fractions*weights)
   wfractions=weights*abundances/sum(weights*abundances) 
  end subroutine   

  function abundancesQ() result(x)
   real :: x(nelements)
   x=abundances
  end function 
  function abundanceQ(i) result(x)
   integer :: i
   real :: x
   x=abundances(i)
  end function 
  function weightsQ() result(x)
   real :: x(nelements)
   x=weights
  end function 
  function fractionsQ() result(x)
   real :: x(nelements)
   x=fractions
  end function 
  function wfractionsQ() result(x)
   real :: x(nelements)
   x=wfractions
  end function 
  function meanmwtQ() result(x)
   real :: x
   x=meanmwt
  end function 
  function nelementsQ() result(x)
   integer :: x
   x=nelements
  end function 
  subroutine elementnamesQ(x)
   character(len=4),dimension(nelements) :: x
   x=elementnames
  end subroutine  
  function wfractionQ(n) result(x)
   integer :: n
   real :: x
   x=wfractions(n)
  end function
  function fractionQ(n) result(x)
   integer :: n
   real :: x
   x=fractions(n)
  end function
  function zQ() result(x)
   real :: x
   x=sum(wfractions(metals))
  end function
  function solarzQ() result(x)
   real :: x
   x=solarz
  end function

end module ElementsMod


!program test
! use ElementsMod
! real weights(8)
! call InitElements

! print*,meanmwt/fractions(1),wfractions(1)

! weights=wfractionsQ()
 
! print*,weights
! print*,weights(1),weights(2),sum(weights(3:8))
! print*, zQ()
!end program 
