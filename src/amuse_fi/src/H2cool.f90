module H2coolMod
 implicit none
 private
 public :: InitH2, H2CoolFunc
 
  logical, parameter :: debug=.FALSE.
  integer, parameter :: coolTableTLength=61
  integer, parameter :: coolTablenLength=9
  real, parameter :: logTmin=1,logTmax=4,dlogT=0.05
  real, parameter :: lognmin=0.3010300,lognmax=4.3010300,dlogn=.5
  real, parameter :: amu=1.6605e-24

  integer, save :: nelements
  integer, save :: H,He,C
  real, save :: meanmwt
  real, save :: Tfac
  
  real, save :: coolUnit=1.
  real,  allocatable, save :: abundances(:)
  real,  allocatable, save :: weights(:)
  real,  allocatable, save :: fractions(:), wfractions(:)
  real,  allocatable, save :: coolTemp(:), cooln(:),coolTable(:,:)
  character(len=4),dimension(:), allocatable, save :: elementnames
  character(len=200),save :: datadir

 contains

 subroutine EndH2
 integer :: test
  deallocate(abundances,weights,fractions,wfractions,coolTemp,cooln, &
  coolTable,elementnames,STAT=test)
  if(test.NE.0) call H2coolError(99)
 end subroutine
 
 subroutine InitH2(indatadir)
   use ElementsMod
  integer :: i,test
  character(len=200) :: indatadir
 
  datadir=indatadir
  if(allocated(abundances)) call EndH2
  call InitElements
  nelements=nelementsQ()
  if(nelements.LE.0) call H2CoolError(0)
  
  allocate(abundances(nelements),weights(nelements),fractions(nelements), &
            wfractions(nelements),coolTemp(coolTableTLength),cooln(coolTablenLength), &
            coolTable(coolTableTLength,coolTablenLength),elementnames(nelements), &
            STAT=test)
   if(test.NE.0) call H2CoolError(1)

   H=elementnrQ("H   ")
   He=elementnrQ("He  ")
   C=elementnrQ("C   ")
   if(H.EQ.0.OR.He.EQ.0.OR.C.EQ.0) call H2CoolError(2)
   
   abundances=abundancesQ()
   weights=weightsQ()
   fractions=fractionsQ()
   meanmwt=meanmwtQ()
   wfractions=wfractionsQ()   
   call elementnamesQ(elementnames)
   
   coolTable=-50.
   call ReadCool

  coolTemp=10**coolTemp
  cooln=10**cooln
  coolTable=10**coolTable
 
  Tfac=0.54
 
 end subroutine

 subroutine ReadCool
  integer:: i,j,imax,jmax
  character(len=215) :: filename
  
  filename=trim(datadir)//'H2.HH2cool'
  open(unit=19,file=filename,status='OLD')
  read(19,*) imax,jmax
  if(imax.ne.coolTableTLength.OR.jmax.ne.coolTablenLength) call H2coolError(3) 
  do j=1,jmax 
   do i=1,imax
    read(19,*) coolTemp(i),cooln(j),coolTable(i,j)
   if(debug) print*,coolTemp(j),coolTable(i,j)
   enddo
  enddo
  close(19)
 end subroutine
 
 function H2Cool(T,n)
  real :: H2cool
  real, intent(in) :: T,n
  real :: ln,dt,logT,dn,logn,lT
  integer :: ti,ti1,ti2,ni,ni1,ni2
  
  ln=max(10**lognmin,n)
  ln=min(10**lognmax,ln)
  lT=max(10**logTmin,T)
  lT=min(10**logTmax,lT)
 
  logT=log10(lT)
  ti=(logT-logTmin)/dlogT+1
  ti1=max(ti,1);ti1=min(coolTableTLength-1,ti1)
  dt=(lT-coolTemp(ti1))/(coolTemp(ti1+1)-coolTemp(ti1))
  
  logn=log10(ln)
  ni=(logn-lognmin)/dlogn+1
  ni1=max(ni,1);ni1=min(coolTablenLength-1,ni1)
  dn=(ln-cooln(ni1))/(cooln(ni1+1)-cooln(ni1))
 
  H2cool=(1-dt)*(1-dn)*coolTable(ti1,ni1)+   &
             dt*(1-dn)*coolTable(ti1+1,ni1)+ &
             (1-dt)*dn*coolTable(ti1,ni1+1)+ &
                 dt*dn*coolTable(ti1+1,ni1+1)
 
 end function

 function H2CoolFunc(T,n) result(x)
  real :: x,ln
  real, intent(in) :: T
  real, intent(in), optional :: n
  
  if(.NOT.present(n)) then
   ln=10**lognmin
  else
   ln=n
  endif  
    
  x=H2Cool(T/Tfac,ln)
  x=(fractions(H)/meanmwt/amu)**2*x*coolUnit
  
 end function
 
 subroutine H2coolError(i)
  integer ::i
  
  print*,' H2 cool error:',i
  stop
 
 end subroutine

end module H2coolMod

