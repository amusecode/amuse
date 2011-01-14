subroutine pmgravsum(p,option)
 use pmgravMod
 include 'globals.h'
 character (LEN=4), intent(in):: option
 integer, intent(in):: p
 real :: peps,ppos(3),pacc(3),ppot

 peps=epsgrav(p)
 ppos=pos(p,1:3)
 pacc=0.
 ppot=0.
 select case (option)
 case('acc ')
  call pmgravaccpot(ppos,pacc=pacc)
 case('pot ')
  call pmgravaccpot(ppos,ppot=ppot)
 case default
  call pmgravaccpot(ppos,pacc=pacc,ppot=ppot)
 end select
 acc(p,1:3)=acc(p,1:3)+pacc   
 phi(p)=phi(p)+ppot
end subroutine pmgravsum

subroutine pgravsum(p,nterms,option,lesofttot)
 include 'globals.h'
 character (LEN=4), intent(in):: option
 integer, intent(in):: p,nterms
 real,intent(inout):: lesofttot
 real :: peps,ppos(3),pacc(3),pphi

 peps=epsgrav(p)
 ppos=pos(p,1:3)
 call gravsum(peps,ppos,pphi,pacc,nterms,option,lesofttot)
 lesofttot=lesofttot*mass(p)*0.5
 select case (option)
 case('acc ')
   acc(p,1:3)=acc(p,1:3)+pacc   
 case('pot ')
   phi(p)=phi(p)+pphi+1.4/epsgrav(p)*mass(p) ! correct for self interaction potential
 case default
   acc(p,1:3)=acc(p,1:3)+pacc   
   phi(p)=phi(p)+pphi+1.4/epsgrav(p)*mass(p)
 end select
end subroutine pgravsum

subroutine gravsum(peps,ppos,pphi,pacc,nterms,option,lesofttot)
 include 'globals.h'
 character (LEN=4), intent(in):: option
 integer, intent(in):: nterms
 real,intent(in) :: peps,ppos(3)
 real,intent(inout):: lesofttot,pphi,pacc(3)

 pacc=0.
 pphi=0.
 lesofttot=0.
 select case (option)
 case('acc ')
  call gravsuma(peps,ppos,pacc,nterms,lesofttot)
 case('pot ')
  call gravsump(peps,ppos,pphi,nterms,lesofttot)
 case default
  call gravsumpa(peps,ppos,pphi,pacc,nterms,lesofttot)
 end select
end subroutine

! epseff=peps**2+epsgrav**2 expected
subroutine gravsuma(peps,ppos,pacc,nterms,lesofttot)
 include 'globals.h'
 integer, intent(in):: nterms
 real,intent(in) :: peps,ppos(3)
 real,intent(inout):: lesofttot,pacc(3)
 real :: dx,dy,dz,dist2,dist,epseff,epseffi,rinveff,drdeldrg, &
  &  r2inveff,r3inveff,phiquad,qr5inv,acci,drsm,accfac,rhinv
 integer :: i,q,smindex
 accfac=1
 do i=1,nterms
   q=bodlist(i)
   dx=ppos(1)-pos(q,1)
   dy=ppos(2)-pos(q,2)
   dz=ppos(3)-pos(q,3)
   dist2=dx*dx+dy*dy+dz*dz
   if(usepm) then
    if(dist2.GT.rcut2) cycle
    drdeldrg=dist2*pmdr
    smindex=FLOOR(drdeldrg)
    drsm=drdeldrg-smindex
    accfac=pmacc(smindex)*(1-drsm)+drsm*pmacc(smindex+1)
   endif
   dist=SQRT(dist2)     
   epseff=peps+epsgrav(q)  
   if(dist.GE.epseff) then
    rinveff=1./dist
    r3inveff=rinveff**3
   else
    epseffi=1./epseff
    rhinv=dist*epseffi
    if(rhinv.LT.0.5) then
      r3inveff=4./3.+rhinv**2*(4.*rhinv-4.8)
      rinveff=1.4-rhinv**2*(8./3.+rhinv**2*(3.2*rhinv-4.8))
    else
      r3inveff=8./3.-6.*rhinv+4.8*rhinv**2- &
            4./3.*rhinv**3-1./120./rhinv**3
      rinveff=1.6-1/(30.*rhinv)-rhinv**2* &
           (16./3.+rhinv*(-8.+rhinv*(4.8-rhinv*16./15.)))
    endif   
    rinveff=rinveff*2*epseffi   
    r3inveff=r3inveff*8*epseffi**3
   endif
              
   acci=accfac*mass(q)*r3inveff    
   pacc(1)=pacc(1)-dx*acci
   pacc(2)=pacc(2)-dy*acci
   pacc(3)=pacc(3)-dz*acci
   
   if(.not.usepm) then
   if(usequad.AND.q.GE.nbods1) then
    r2inveff=rinveff*rinveff
    qr5inv=r3inveff*r2inveff
    phiquad=(-.5*((dx**2-dz**2)*quad(q,1)+(dy**2-dz**2)*quad(q,4)) &
            -(dx*dy*quad(q,2)+dx*dz*quad(q,3)+dy*dz*quad(q,5)))*qr5inv
    phiquad=5.*phiquad*r2inveff
    pacc(1)=pacc(1)+dx*phiquad+(dx*quad(q,1)+dy*quad(q,2)+dz*quad(q,3))*qr5inv
    pacc(2)=pacc(2)+dy*phiquad+(dy*quad(q,4)+dx*quad(q,2)+dz*quad(q,5))*qr5inv
    pacc(3)=pacc(3)+dz*phiquad &
            +(dz*(-quad(q,1)-quad(q,4))+dx*quad(q,3)+dy*quad(q,5))*qr5inv
   endif
   endif
 enddo
end subroutine gravsuma

subroutine gravsump(peps,ppos,pphi,nterms,lesofttot)
 include 'globals.h'
 integer, intent(in):: nterms
 real,intent(in) :: peps,ppos(3)
 real,intent(inout):: lesofttot,pphi
 real :: dx,dy,dz,dist2,dist,epseff,epseffi,drdeldrg,rinveff, &
  &  r2inveff,r3inveff,phiquad,qr5inv,acci,drsm,potfac,rhinv
 integer :: i,q,smindex
 potfac=1 
 do i=1,nterms
   q=bodlist(i)
   dx=ppos(1)-pos(q,1)
   dy=ppos(2)-pos(q,2)
   dz=ppos(3)-pos(q,3)
   dist2=dx*dx+dy*dy+dz*dz
   if(usepm) then
    if(dist2.GT.rcut2) cycle
    drdeldrg=dist2*pmdr
    smindex=FLOOR(drdeldrg)
    drsm=drdeldrg-smindex
    potfac=pmpot(smindex)*(1-drsm)+drsm*pmpot(smindex+1)
   endif
   dist=SQRT(dist2)     
   epseff=peps+epsgrav(q)  
   if(dist.GE.epseff) then
    rinveff=1./dist
    r3inveff=rinveff**3
   else
    epseffi=1./epseff
    rhinv=dist*epseffi
    if(rhinv.LT.0.5) then
      r3inveff=4./3.+rhinv**2*(4.*rhinv-4.8)
      rinveff=1.4-rhinv**2*(8./3.+rhinv**2*(3.2*rhinv-4.8))
    else
      r3inveff=8./3.-6.*rhinv+4.8*rhinv**2- &
            4./3.*rhinv**3-1./120./rhinv**3
      rinveff=1.6-1/(30.*rhinv)-rhinv**2* &
           (16./3.+rhinv*(-8.+rhinv*(4.8-rhinv*16./15.)))
    endif   
    rinveff=rinveff*2*epseffi   
    r3inveff=r3inveff*8*epseffi**3
   endif
   
   pphi=pphi-potfac*mass(q)*rinveff
   lesofttot=lesofttot+mass(q)*(rinveff-dist2*r3inveff)
   
   if(.not.usepm) then   
   if(usequad.AND.q.GE.nbods1) then   
    r2inveff=rinveff*rinveff
    qr5inv=r3inveff*r2inveff
    phiquad=(-.5*((dx**2-dz**2)*quad(q,1)+(dy**2-dz**2)*quad(q,4)) &
            -(dx*dy*quad(q,2)+dx*dz*quad(q,3)+dy*dz*quad(q,5)))*qr5inv
    pphi=pphi+phiquad
   endif
   endif
 enddo
         
end subroutine gravsump

subroutine gravsumpa(peps,ppos,pphi,pacc,nterms,lesofttot)
 include 'globals.h'
 integer, intent(in):: nterms
 real,intent(in) :: peps,ppos(3)
 real,intent(inout):: lesofttot,pphi,pacc(3)
 real :: dx,dy,dz,dist2,dist,epseff,epseffi,drdeldrg,rinveff, &
  &  r2inveff,r3inveff,phiquad,qr5inv,acci,drsm,accfac,potfac,rhinv
 integer :: i,q,smindex

 accfac=1
 potfac=1 
 do i=1,nterms
   q=bodlist(i)
   dx=ppos(1)-pos(q,1)
   dy=ppos(2)-pos(q,2)
   dz=ppos(3)-pos(q,3)
   dist2=dx*dx+dy*dy+dz*dz
   if(usepm) then
    if(dist2.GT.rcut2) cycle
    drdeldrg=dist2*pmdr
    smindex=FLOOR(drdeldrg)
    drsm=drdeldrg-smindex
    accfac=pmacc(smindex)*(1-drsm)+drsm*pmacc(smindex+1)
    potfac=pmpot(smindex)*(1-drsm)+drsm*pmpot(smindex+1)
   endif
   dist=SQRT(dist2)     
   epseff=peps+epsgrav(q)  
   if(dist.GE.epseff) then
    rinveff=1./dist
    r3inveff=rinveff**3
   else
    epseffi=1./epseff
    rhinv=dist*epseffi
    if(rhinv.LT.0.5) then
      r3inveff=4./3.+rhinv**2*(4.*rhinv-4.8)
      rinveff=1.4-rhinv**2*(8./3.+rhinv**2*(3.2*rhinv-4.8))
    else
      r3inveff=8./3.-6.*rhinv+4.8*rhinv**2- &
            4./3.*rhinv**3-1./120./rhinv**3
      rinveff=1.6-1/(30.*rhinv)-rhinv**2* &
           (16./3.+rhinv*(-8.+rhinv*(4.8-rhinv*16./15.)))
    endif   
    rinveff=rinveff*2*epseffi   
    r3inveff=r3inveff*8*epseffi**3
   endif
   
      pphi=pphi-potfac*mass(q)*rinveff
      lesofttot=lesofttot+mass(q)*(rinveff-dist2*r3inveff)
   
    acci=accfac*mass(q)*r3inveff
    pacc(1)=pacc(1)-dx*acci
    pacc(2)=pacc(2)-dy*acci
    pacc(3)=pacc(3)-dz*acci
   
   if(.not.usepm) then
   if(usequad.AND.q.GE.nbods1) then
   
    r2inveff=rinveff*rinveff
    qr5inv=r3inveff*r2inveff
    phiquad=(-.5*((dx**2-dz**2)*quad(q,1)+(dy**2-dz**2)*quad(q,4)) &
            -(dx*dy*quad(q,2)+dx*dz*quad(q,3)+dy*dz*quad(q,5)))*qr5inv
    pphi=pphi+phiquad
    phiquad=5.*phiquad*r2inveff
    pacc(1)=pacc(1)+dx*phiquad+(dx*quad(q,1)+dy*quad(q,2)+dz*quad(q,3))*qr5inv
    pacc(2)=pacc(2)+dy*phiquad+(dy*quad(q,4)+dx*quad(q,2)+dz*quad(q,5))*qr5inv
    pacc(3)=pacc(3)+dz*phiquad &
            +(dz*(-quad(q,1)-quad(q,4))+dx*quad(q,3)+dy*quad(q,5))*qr5inv
   endif
   endif
 enddo
         
end subroutine gravsumpa

subroutine pfuvsum(p,nterms)
 include 'globals.h'
 integer, intent(in):: p,nterms
 real :: peps,ppos(3),pheat
 peps=epsgrav(p)
 ppos=pos(p,1:3)
 pheat=0.
 call fuvsum(peps,ppos,pheat,nterms)
 fuvheat(p)=fuvheat(p)+pheat         
end subroutine pfuvsum

subroutine fuvsum(peps,ppos,pheat,nterms)
 include 'globals.h'
 integer, intent(in):: nterms
 real,intent(in) :: ppos(3),peps
 real, intent(inout) :: pheat
 real :: pmass,dx,dy,dz,drdotdr,epseff,epseffi,drdeldrg,distance,r2inveff,drsm
 integer :: i,q,smindex
 real, parameter :: uvfac=uvinterp/uvedge
 
 do i=1,nterms
   q=bodlist(i)
   dx=ppos(1)-pos(q,1)
   dy=ppos(2)-pos(q,2)
   dz=ppos(3)-pos(q,3)
   IF(dx.GE.hboxsize.AND.periodic) dx=dx-pboxsize
   IF(dx.LT.-hboxsize.AND.periodic) dx=dx+pboxsize
   IF(dy.GE.hboxsize.AND.periodic) dy=dy-pboxsize
   IF(dy.LT.-hboxsize.AND.periodic) dy=dy+pboxsize
   IF(dz.GE.hboxsize.AND.periodic) dz=dz-pboxsize
   IF(dz.LT.-hboxsize.AND.periodic) dz=dz+pboxsize
   drdotdr=dx*dx+dy*dy+dz*dz
   epseff=peps+epsgrav(q)
   epseffi=1./epseff
   drdeldrg=drdotdr*epseffi*epseffi
   if(drdeldrg.GE.uvedge) then
    r2inveff=1/drdotdr
   else
    drsm=uvfac*drdeldrg
    smindex=INT(drsm)
    drsm=drsm-smindex
    r2inveff=(1-drsm)*smoothuv(smindex)+drsm*smoothuv(1+smindex)
    r2inveff=r2inveff*4*epseffi**2
   endif
   
   if(optdepth.LT.tiny) then
    pheat=pheat+starfuv(q)*r2inveff
   else
    distance=sqrt( drdotdr)
    pheat=pheat+starfuv(q)*r2inveff*exp(-optdepth*distance)
   endif
   
 enddo
end subroutine fuvsum
