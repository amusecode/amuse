! ionization eq solved including cosmic rays + photoion.

module IonMod
 implicit none
 private
 public :: InitIon,Ion,EndIon,Ionh
 
  logical, parameter :: dbg=.FALSE.
  integer, parameter :: ntemptable=256,ncrtable=256
  real, parameter :: logTmin=1,logTmax=9,logcrmin=-6,logcrmax=8
  real, parameter :: dlogT=(logTmax-logTmin)/ntemptable
  real, parameter :: dlogcr=(logcrmax-logcrmin)/ncrtable
  real, parameter :: amu=1.6605e-24
  real, parameter :: crunit=1.e-17
  real, parameter :: Eel=40
  real, parameter :: tablerho=1.

 
  integer, save :: nelements
  integer, save :: H,He,C
  real,save :: meanmwt, minion,tablecr


  real, dimension(ntemptable), save :: logTin
  real, dimension(ncrtable), save :: logcrin
  real, dimension(ntemptable,ncrtable),save :: logTout,logeout,loghout, &
                                                loghe1out,loghe2out

  real,  allocatable, save :: abundances(:)
  real,  allocatable, save :: weights(:)
  real,  allocatable, save :: fractions(:), wfractions(:)
  character(len=4),dimension(:), allocatable, save :: elementnames
    
 contains
 
  subroutine ion(rtemp,temp,elecfrac,rho,cr)
   real,intent(in) :: rtemp
   real,intent(out):: temp,elecfrac
   real,intent(in),optional :: rho,cr
   real lrho,lcr
   
   if(.not.present(cr)) then 
    lcr=tablecr/crunit
   else
    lcr=cr
   endif
   
   if(.not.present(rho)) then
    lrho=tablerho 
   else
    lrho=rho
   endif

   call eTinterpol(rtemp,temp,elecfrac,lrho,lcr)
  end subroutine ion

  subroutine ionh(rtemp,temp,xe,xhii,xheii,xheiii,rho,cr)
   real,intent(in) :: rtemp
   real,intent(out):: temp
   real,intent(in),optional :: rho,cr
   real,intent(out) :: xe,xhii,xheii,xheiii
   real lrho,lcr
   
   if(.not.present(cr)) then 
    lcr=tablecr/crunit
   else
    lcr=cr
   endif
   
   if(.not.present(rho)) then
    lrho=tablerho 
   else
    lrho=rho
   endif

   call hTinterpol(rtemp,temp,xe,xhii,xheii,xheiii,lrho,lcr)
  end subroutine ionh

 
  subroutine Ioniter(rtemp,temp,x,rho,cr)
   real,intent(in) :: rtemp,rho,cr
   real,intent(out) :: temp,x
   real :: xnew,ttemp
   integer :: j
  
    call eTinterpol(rtemp,temp,x,rho,cr)
    call ionize(temp,cr,rho,xnew)
    j=0
    do while(ABS((x-xnew)/x).GT.0.001.AND.j.LT.200)
    j=j+1
    x=(xnew+x)/2.
    temp=tempIon(rtemp,x)
    call ionize(temp,cr,rho,xnew)
    enddo
    x=(xnew+x)/2.
  end subroutine Ioniter
 
  subroutine EndIon
   integer :: test
   deallocate(abundances,weights,fractions, elementnames,wfractions,STAT=test)
   if(test.NE.0) call IonError(99)
  end subroutine EndIon
  
  subroutine InitIon(unitm_in_msun,unitl_in_kpc,cr)
  use ElementsMod
   real, intent(in) :: unitm_in_msun,unitl_in_kpc
   real, intent(in), optional :: cr
   integer :: test
  
   if(present(cr)) then
    tablecr=cr*crunit
   else
    tablecr=1.8*crunit
   endif
  
   if(allocated(abundances)) call EndIon

   call InitElements
   nelements=nelementsQ()
   if(nelements.LE.0) call IonError(0)
  
   allocate(abundances(nelements),weights(nelements),fractions(nelements), &
            wfractions(nelements),elementnames(nelements),STAT=test)
   if(test.NE.0) call IonError(1)

   H=elementnrQ("H   ")
   He=elementnrQ("He  ")
   C=elementnrQ("C   ")
   if(H.EQ.0.OR.He.EQ.0.OR.C.EQ.0) call IonError(2)
   
   abundances=abundancesQ()
   weights=weightsQ()
   fractions=fractionsQ()
   meanmwt=meanmwtQ()
   wfractions=wfractionsQ()   
   call elementnamesQ(elementnames)
   
   minion=fractions(C)/fractions(H)
   
   if(dbg) print*,"InitIon:",minion
   
   call preparetemp
  end subroutine InitIon
  
  real function tempIon( t, elecfrac)
   real :: t,elecfrac 
   tempion=t/(1+fractions(H)*elecfrac)
  end function
  
  subroutine eTinterpol(tin,t,e,rho,crin)
   real :: tin,t,e,logT,dt,h
   real :: lcr,dc,cr,rho,crin
   integer :: ti,ti1,ci,ci1
   
    cr=crin/rho

   if(dbg) then
    if(cr.LT.10**logcrmin.OR.cr.GT.10**logcrmax) call IonError(3)
   endif

   if(dbg) then
    if(tin.LT.10**logTmin.OR.tin.GT.10**logTmax) call IonError(4)
   endif
   
   logT=log10(tin)
   ti=(logT-logTmin)/dlogT+1
   ti1=max(ti,1);ti1=min(ntemptable-1,ti1)
   dt=(logT-logTmin)/dlogT-ti+1

   lcr=log10(cr)
   ci=(lcr-logcrmin)/dlogcr+1
   ci1=max(ci,1);ci1=min(ncrtable-1,ci1)
   dc=(lcr-logcrmin)/dlogcr+1-ci 

   t=(1-dt)*(1-dc)*logTout(ti1,ci1)+dt*(1-dc)*logTout(ti1+1,ci1)+ &
     (1-dt)*dc*logTout(ti1,ci1+1)+dt*dc*logTout(ti1+1,ci1+1)
   t=10.**t
   e=(1-dt)*(1-dc)*logeout(ti1,ci1)+dt*(1-dc)*logeout(ti1+1,ci1)+ &
     (1-dt)*dc*logeout(ti1,ci1+1)+dt*dc*logeout(ti1+1,ci1+1)
   e=10.**e
  end subroutine 

  subroutine hTinterpol(tin,t,e,xhii,xheii,xheiii,rho,crin)
   real :: tin,t,e,logT,dt,xhii,xheii,xheiii
   real :: lcr,dc,cr,rho,crin
   integer :: ti,ti1,ci,ci1
   
    cr=crin/rho

   if(dbg) then
    if(cr.LT.10**logcrmin.OR.cr.GT.10**logcrmax) call IonError(3)
   endif

   if(dbg) then
    if(tin.LT.10**logTmin.OR.tin.GT.10**logTmax) call IonError(4)
   endif
   
   logT=log10(tin)
   ti=(logT-logTmin)/dlogT+1
   ti1=max(ti,1);ti1=min(ntemptable-1,ti1)
   dt=(logT-logTmin)/dlogT-ti+1

   lcr=log10(cr)
   ci=(lcr-logcrmin)/dlogcr+1
   ci1=max(ci,1);ci1=min(ncrtable-1,ci1)
   dc=(lcr-logcrmin)/dlogcr+1-ci 

   t=(1-dt)*(1-dc)*logTout(ti1,ci1)+dt*(1-dc)*logTout(ti1+1,ci1)+ &
     (1-dt)*dc*logTout(ti1,ci1+1)+dt*dc*logTout(ti1+1,ci1+1)
   t=10.**t
   e=(1-dt)*(1-dc)*logeout(ti1,ci1)+dt*(1-dc)*logeout(ti1+1,ci1)+ &
     (1-dt)*dc*logeout(ti1,ci1+1)+dt*dc*logeout(ti1+1,ci1+1)
   e=10.**e
   xhii=(1-dt)*(1-dc)*loghout(ti1,ci1)+dt*(1-dc)*loghout(ti1+1,ci1)+ &
     (1-dt)*dc*loghout(ti1,ci1+1)+dt*dc*loghout(ti1+1,ci1+1)
   xhii=10.**xhii
   xheii=(1-dt)*(1-dc)*loghe1out(ti1,ci1)+dt*(1-dc)*loghe1out(ti1+1,ci1)+ &
     (1-dt)*dc*loghe1out(ti1,ci1+1)+dt*dc*loghe1out(ti1+1,ci1+1)
   xheii=10.**xheii
   xheiii=(1-dt)*(1-dc)*loghe2out(ti1,ci1)+dt*(1-dc)*loghe2out(ti1+1,ci1)+ &
     (1-dt)*dc*loghe2out(ti1,ci1+1)+dt*dc*loghe2out(ti1+1,ci1+1)
   xheiii=10.**xheiii
  end subroutine 


  
  subroutine preparetemp
   integer :: i,j,k
   real :: temp,rtemp, rho, cr,x,xnew,xhii,xheii,xheiii
   rho=tablerho
   do j=1,ncrtable
    do i=1,ntemptable
      logTin(i)=logTmin+(i-1)*dlogT
      rtemp=10**logTin(i)
      logcrin(j)=logcrmin+(j-1)*dlogcr
      cr=10**logcrin(j)
      x=0.1
      xnew=0.
      k=0
      do while(ABS((x-xnew)/x).GT.0.00001.AND.k.LT.200)
       k=k+1
       x=(xnew+x)/2.
       temp=tempIon(rtemp,x)
       call ionize(temp,cr,rho,xnew,xhii,xheii,xheiii)
      enddo
      x=(xnew+x)/2.
      logTout(i,j)=log10(temp)
      logeout(i,j)=log10(x)
      loghout(i,j)=log10(xhii)   
      loghe1out(i,j)=log10(xheii)   
      loghe2out(i,j)=log10(xheiii)   
    enddo
   enddo
  end subroutine
    
  subroutine ionize(temp,crrate,dens,elecfrac,xhii,xheii,xheiii)
      real :: temp,crrate,dens,elecfrac
      real,optional :: xhii,xheii,xheiii
      real :: ndens
      real :: xh,xhe1,xhe2
      real :: crh
      real :: cih,cihe1,cihe2
      real :: arh,arhe1,arhe2
      real :: abund
      integer :: i     
      
      ndens=dens

      call cfit(1,1,temp,cih)
      call rrfit(1,1,temp,arh)
      call cfit(2,2,temp,cihe1)
      call rrfit(2,2,temp,arhe1)
      call cfit(2,1,temp,cihe2)
      call rrfit(2,1,temp,arhe2)
      cih=cih*ndens
      arh=ndens*arh
      cihe1=cihe1*ndens
      arhe1=ndens*arhe1      
      cihe2=cihe2*ndens
      arhe2=ndens*arhe2
      
      elecfrac=minion+cih/(cih+arh)
      
      do i=1,8
      
      crh=crunit*crrate*(1+fih(elecfrac)+fihe(elecfrac))/elecfrac
      
      xh=(cih+crh)/(cih+crh+arh)      
      xhe1=(cihe1+crh)/((cihe1+crh)*(1+(cihe2+crh)/arhe2)+arhe1) 
      xhe2=xhe1*(cihe2+crh)/arhe2
!      elecfrac=.75*( xh+fractions(He)/fractions(H)*(xhe1+2*xhe2)+minion)+ &
!               .25*elecfrac
      elecfrac=.75*( xh+abundances(He)*(xhe1+2*xhe2)+minion)+ &
               .25*elecfrac
      enddo
      if(present(xhii)) xhii=xh               
      if(present(xheii)) xheii=abundances(He)*xhe1               
      if(present(xheiii)) xheiii=abundances(He)*xhe2               
      end subroutine

  
  
 real function fih(x)
  real :: x
  fih=0;if(x.gt.1) return
  fih=(Eel/13.6-1.)*.39*(1-x**.41)**1.76*(Eel/1000.)**(2.3*x)
 end function
      
 real function fihe(x)
  real :: x
  fihe=0;if(x.gt.1) return
  fihe=(Eel/24.6-1.)*0.055*(1-x**.46)**1.67*(Eel/1000.)**(2.3*x)
 end function
  
        subroutine cfit(iz,in,t,c)
!******************************************************************************
!*** This subroutine calculates rates of direct collisional ionization 
!*** for all ionization stages of all elements from H to He (Z=2)
!*** by use of the fits from G. S. Voronov, 1997, ADNDT, 65, 1
!*** Input parameters:  iz - atomic number 
!***                    in - number of electrons from 1 to iz 
!***                    t  - temperature, K
!*** Output parameter:  c  - rate coefficient, cm^3 s^(-1)
!******************************************************************************
       integer :: iz,in,i
       real :: t,c,te,u 
       REAL CF(5,2,2)
      
      DATA(CF(I, 1, 1),I=1,5)/   13.6,0.,2.91E-08,0.2320,0.39/
      DATA(CF(I, 2, 2),I=1,5)/   24.6,0.,1.75E-08,0.1800,0.35/
      DATA(CF(I, 2, 1),I=1,5)/   54.4,1.,2.05E-09,0.2650,0.25/
      
      c=0.0
       te=t*8.617385e-05
       u=cf(1,iz,in)/te
       if(u.gt.80.0)return
      c=cf(3,iz,in)*(1.0+cf(2,iz,in)*sqrt(u))/(cf(4,iz,in)+u)*u**cf(5,iz,in)*exp(-u)
      return
      end subroutine

            subroutine rrfit(iz,in,t,r)
!*** adapted from D. A. Verner, verner@pa.uky.edu 
!******************************************************************************
!*** This subroutine calculates rates of radiative recombination for all ions
!*** of all elements from H and He, case b rcombination for H:
!*** Input parameters:  iz - atomic number 
!***                    in - number of electrons from 1 to iz 
!***                    t  - temperature, K
!*** Output parameter:  r  - rate coefficient, cm^3 s^(-1)
!******************************************************************************
       real :: t,r,tt
       integer :: iz,in,i
       REAL rnew(4,2,2)
             
      data(rnew(i, 1, 1),i=1,4)/7.782e-11,0.7070,2.781e+00,1.656e+05/
      data(rnew(i, 2, 1),i=1,4)/1.891e-10,0.7524,9.370e+00,2.774e+06/
      data(rnew(i, 2, 2),i=1,4)/9.356e-10,0.7892,4.266e-02,4.677e+06/

      r=0.0
       tt=sqrt(t/rnew(3,iz,in))
       r=rnew(1,iz,in)/(tt*(tt+1.0)**(1.0-rnew(2,iz,in))* &
             (1.0+sqrt(t/rnew(4,iz,in)))**(1.0+rnew(2,iz,in)))
      return 
      end subroutine

      subroutine IonError(i)
       integer i
       
       print*,'Ion error:',i
       stop
      end subroutine
  
 end module IonMod
 
!program testion
! use ElementsMod
! use IonMod
! use CoolingMod
! integer i,j
! real x,xh,rtemp,temp,e,rho,cr,lambda,fractions(8),meanmwt
! character(len=200),save :: datadir
! datadir='/home/pelupes/fish/stars/data/'

! call InitCooling(datadir) 
! call InitIon(1.e9,1.)
! fractions=fractionsQ()
! meanmwt=meanmwtQ()
! amu=1.6605e-24
! cr=1800
! rho=.1
!10 read*,rtemp
!if(rtemp.eq.0) stop
! call ion(rtemp,temp,x,rho,cr)
! call ionh(rtemp,temp,xh,rho,cr)
! print*,temp,x,xh
! goto 10
!  rho=.3
!  x=.1
!   do i=1,256
!    rtemp=10**(1+(i-1)/32.)
!    call ion(rtemp,temp,x,rho,cr)
!    lambda=CoolFunc(x,temp)
!    print*,alog10(temp),alog10((amu/fractions(1)*meanmwt)**2*lambda)
!   enddo
!   do j=1,nelementsQ()
!   do i=1,256
!    rtemp=10**(1+(i-1)/32.)
!    call ion(rtemp,temp,x,rho,cr)
!    lambda=ElementCool(x,max(10.,temp),j)
!    print*,alog10(temp),alog10(lambda)
!   enddo
!   enddo


!   print*, ' '
!  enddo
 
!end program testion  
