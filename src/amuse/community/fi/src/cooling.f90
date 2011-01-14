! module for cooling
! usage: 
!  - include module "use CoolingMod"
!  - initialization by calling InitCooling once with directory of data files
!  - call Coolfunc(ne,T) 
!  - call EndCooling when done
! (see test program at the end for example)
! uses module Elementsmod 

module CoolingMod
 implicit none
 private
 public :: InitCooling, Coolfunc,ElementCool,ElementCoolFunc,EndCooling
 
 public :: Cool_nelements

  logical, parameter :: dbg=.FALSE.
  integer, parameter :: coolTableLength=160
  real, parameter :: logTmin=1,logTmax=9,dlogT=0.05
  real, parameter :: amu=1.6605e-24

  integer, save :: nelements
  integer, save :: H,He,C
  integer, save :: Cool_nelements
  real, save :: meanmwt
  real, save :: coolUnit=1.
  
  real,  allocatable, save :: abundances(:)
  real,  allocatable, save :: weights(:)
  real,  allocatable, save :: fractions(:), wfractions(:)
  real,  allocatable, save :: coolTemp(:), coolTable(:,:), coolTable2(:,:)
  character(len=4),dimension(:), allocatable, save :: elementnames
  character(len=200),save :: datadir
  
 contains
 
  subroutine EndCooling
   integer :: test
   deallocate(abundances,weights,fractions, elementnames, &
              wfractions,coolTemp,coolTable,coolTable2,STAT=test)
   if(test.NE.0) call CoolError(99)
  end subroutine EndCooling
 
  subroutine InitCooling(indatadir)
  use ElementsMod
  integer :: i,test
  character(len=200) :: indatadir
 
   datadir=indatadir
   if(allocated(abundances)) call EndCooling

   call InitElements
   nelements=nelementsQ()
   Cool_nelements=nelements
   if(nelements.LE.0) call CoolError(0)
  
   allocate(abundances(nelements),weights(nelements),fractions(nelements), &
            wfractions(nelements),coolTemp(coolTableLength), &
            coolTable(nelements,coolTableLength),elementnames(nelements), &
            coolTable2(nelements,coolTableLength),STAT=test)
   if(test.NE.0) call CoolError(1)

   H=elementnrQ("H   ")
   He=elementnrQ("He  ")
   C=elementnrQ("C   ")
   if(H.EQ.0.OR.He.EQ.0.OR.C.EQ.0) call CoolError(2)
   
   abundances=abundancesQ()
   weights=weightsQ()
   fractions=fractionsQ()
   meanmwt=meanmwtQ()
   wfractions=wfractionsQ()   
   call elementnamesQ(elementnames)
   
   coolTable=-50.
   coolTable2=-50.
   do i=1,nelements
    call ReadCool(i)
   enddo

  coolTemp=10**coolTemp
  coolTable=10**coolTable
  coolTable2=10**coolTable2

 end subroutine InitCooling
  
  
  subroutine ReadCool(i)
   integer :: i,j,jmax
   character(len=215) :: filename
   
   filename=trim(datadir)//trim(elementNames(i))//'.cool'
   if(dbg)print*,filename
   open(unit=19,file=filename,status='OLD')
    read(19,*) jmax
    do j=1,jmax
      read(19,*) coolTemp(j),coolTable(i,j)
   if(dbg) print*,coolTemp(j),coolTable(i,j)
    enddo
   close(19)    
   
   filename=trim(datadir)//trim(elementNames(i))//'.Hcool'
   if(dbg) print*,filename
   open(unit=19,file=filename,status='OLD')
    read(19,*) jmax
    do j=1,jmax
      read(19,*) coolTemp(j),coolTable2(i,j)
   if(dbg) print*,coolTemp(j),coolTable2(i,j)
    enddo
   close(19)    
     
  end subroutine ReadCool
  
  function ElementCool(ne,T,elem)
   real :: ElementCool,HElementCool
   real,intent(in) ::ne,T
   integer,intent(in) :: elem
   real :: dt,logT
   integer :: ti,ti1,ti2
   
   if(dbg) then
    if(ne.LT.0) call coolError(1)
    if(T.LT.10**logTmin.OR.T.GT.10**logTmax) call coolError(1)
   endif
   
   logT=log10(T)
   ti=(logT-logTmin)/dlogT+1
   ti1=max(ti,1);ti1=min(coolTableLength-1,ti1)
   dt=(T-coolTemp(ti1))/(coolTemp(ti1+1)-coolTemp(ti1))
   ElementCool=(1.-dt)*coolTable(elem,ti1)+dt*coolTable(elem,ti1+1)
   HElementCool=(1.-dt)*coolTable2(elem,ti1)+dt*coolTable2(elem,ti1+1)
   ElementCool=ne*ElementCool    
   ElementCool=ElementCool+HElementCool

  end function ElementCool
  
  function ElementCoolFunc(ne,T,elem,chm)
   real :: ElementCoolFunc
   real,intent(in) ::ne,T
   integer,intent(in) :: elem
   real,intent(in), optional ::chm(1:nelements)
   real :: chem(1:nelements)
   
   if(.not.present(chm)) then 
    chem=abundances
   else
    chem=chm
   endif 
   
   ElementCoolFunc=ElementCool(ne,T,elem)*chem(elem)* &
    (fractions(H)/meanmwt/amu)**2*coolUnit
  end function ElementCoolFunc

  
  
  function xeCool(T,chm)
   real :: xeCool
   real,intent(in) ::T
   real,intent(in), optional ::chm(1:nelements)
   real :: dt,logT,coolinterp(1:nelements),chem(1:nelements)
   integer :: ti,ti1,ti2
   
   if(.not.present(chm)) then 
    chem=abundances
   else
    chem=chm
   endif 
   
   if(dbg) then
    if(T.LT.10**logTmin.OR.T.GT.10**logTmax) call coolError(1)
   endif
   
   logT=log10(T)
   ti=(logT-logTmin)/dlogT+1
   ti1=max(ti,1);ti1=min(coolTableLength-1,ti1)
   dt=(T-coolTemp(ti1))/(coolTemp(ti1+1)-coolTemp(ti1))
   coolinterp=(1.-dt)*coolTable(1:nelements,ti1)+dt*coolTable(1:nelements,ti1+1)
   xeCool=sum(chem*coolinterp)    
  end function xeCool
  
  function HCool(T,chm)
   real :: HCool
   real,intent(in) :: T
   real,intent(in), optional ::chm(1:nelements)
   real :: dt,logT,coolinterp(1:nelements),chem(1:nelements)
   integer :: ti,ti1,ti2
   
   if(.not.present(chm)) then 
    chem=abundances
   else
    chem=chm
   endif 
   
   if(dbg) then
    if(T.LT.10**logTmin.OR.T.GT.10**logTmax) call coolError(1)
   endif
   logT=log10(T)
   ti=(logT-logTmin)/dlogT+1
   ti1=max(ti,1);ti1=min(coolTableLength-1,ti1)
   dt=(T-coolTemp(ti1))/(coolTemp(ti1+1)-coolTemp(ti1))
   coolinterp=(1.-dt)*coolTable2(1:nelements,ti1)+dt*coolTable2(1:nelements,ti1+1)
   HCool=sum(chem*coolinterp)    
  end function HCool
  
  function CoolFunc(ne,T,chm)
   real :: CoolFunc
   real,intent(in) :: ne,T
   real,intent(in), optional ::chm(1:nelements)
   real :: Tl
   Tl=T 
   if(Tl.LT.10**logTmin) Tl=10**logTmin
   if(.not.present(chm)) then
    CoolFunc=ne*xeCool(Tl)+HCool(Tl)   ! not ok. but for low ne
    CoolFunc=(fractions(H)/meanmwt/amu)**2*CoolFunc*coolUnit 
   else
    CoolFunc=ne*xeCool(Tl,chm)+HCool(Tl,chm)   ! not ok. but for low ne
    CoolFunc=(fractions(H)/meanmwt/amu)**2*CoolFunc*coolUnit 
   endif
  
  end function CoolFunc
  
!  function H2Cool(T) result(x)
!   real :: x,T
!   real,parameter :: rop=1.5
!      
!   x=4.77e-26/(1+rop)*(exp(-510/T)+1.35*rop*exp(-845/T))   

!  end function
  
!  function H2CoolFunc(ne,T,h2frac) result(x)
!   real :: x
!   real,intent(in) :: ne,T,h2frac
!   real :: Tl
!   Tl=T 
!   if(Tl.LT.10**logTmin) Tl=10**logTmin
!   x=ne*xeCool(Tl)+ &
!      (1-h2frac)*HCool(Tl)+h2frac*H2Cool(Tl)    ! not ok. but for low ne
!   x=(fractions(H)/meanmwt/amu)**2*x*coolUnit 
!  end function

  
  subroutine coolError(i)
   integer i
   
   print*,'cool error:',i
   stop
   
  end subroutine 
   
 end module CoolingMod


! program test
! use ElementsMod 
! use CoolingMod
! character(len=200) :: indatadir
! real :: temp
! call InitElements
! indatadir="/home/pelupes/fish/stars/data/"
! call InitCooling(indatadir)

! print*,0.01*xecool(2000.)/elementcool(0.01,2000.,8)

!10 read*, temp
! if(temp.eq.0) stop
! print*,temp,Coolfunc(0.1,10**temp)
!goto 10 
! write(*,5),temp, xecool(temp),HCool(temp),H2cool(temp)
!goto 10 
!5 format( e12.4,' ',e12.4,' ',e12.4,' ',e12.4)
! end program test

