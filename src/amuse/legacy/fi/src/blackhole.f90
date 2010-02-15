! bh is star particle, except:
!  tform indicates time of last 'activity'
!  tfeedb is the accretion rate ( >= 0)
!
subroutine skinnydensity_bh(n,h,ppos,pvel,lrho,ldrhodh, &
                              lhdivv,lhcurlv,lvdisp,lvbulk,lcsound)
  include 'globals.h'
  integer i,iwsm,n,p
  real hsminv,wnorm,distnorm,dr2p,drw,wsm,wmass,dr2, &
       wmass1,dwnorm,vdotdr,dwsm,dwmass,dx,dy,dz, &
       dvx,dvy,dvz,dwmnbi,curlvx,curlvy,curlvz,mv(3),mv2
  real h,ppos(3),pvel(3),lrho,ldrhodh,lhdivv, &
      lhcurlv,lvdisp,lvbulk,lcsound
  lrho=0.
  ldrhodh=0.
  lhdivv=0.
  lhcurlv=0.
  lvdisp=0.
  curlvx=0.;curlvy=0.;curlvz=0.
  mv=0.;mv2=0.
  lcsound=0.
  hsminv=1./h
  wnorm=piinv*hsminv*hsminv*hsminv
  dwnorm=piinv*hsminv**2*hsminv**2*hsminv
  distnorm=hsminv**2*deldr2i
  do i=1,n
    p=srlist(i)
    dx=ppos(1)-pos(p,1)
    dy=ppos(2)-pos(p,2)
    dz=ppos(3)-pos(p,3)		   
    if(dx.GE.hboxsize.AND.periodic) dx=dx-pboxsize
    if(dx.LT.-hboxsize.AND.periodic) dx=dx+pboxsize
    if(dy.GE.hboxsize.AND.periodic) dy=dy-pboxsize
    if(dy.LT.-hboxsize.AND.periodic) dy=dy+pboxsize
    if(dz.GE.hboxsize.AND.periodic) dz=dz-pboxsize
    if(dz.LT.-hboxsize.AND.periodic) dz=dz+pboxsize
    dr2=dx*dx+dy*dy+dz*dz   
    dr2p=dr2*distnorm
    if(ninterp.GE.dr2p) then
      iwsm=INT(dr2p)		   
      drw=dr2p-iwsm
      wsm=(1.-drw)*wsmooth(iwsm)+drw*wsmooth(1+iwsm)
      wmass=wnorm*wsm
      lrho=lrho+mass(p)*wmass
      ldrhodh=ldrhodh-3*hsminv*mass(p)*wmass
      lcsound=lcsound+mass(p)*csound(p)*wmass

      mv(1)=mv(1)+veltpos(p,1)*mass(p)*wmass
      mv(2)=mv(2)+veltpos(p,2)*mass(p)*wmass
      mv(3)=mv(3)+veltpos(p,3)*mass(p)*wmass
      mv2=mv2+(veltpos(p,1)**2+veltpos(p,2)**2+ &
          veltpos(p,3)**2)*mass(p)*wmass

      dwsm=(1.-drw)*dwsmooth(iwsm)+drw*dwsmooth(1+iwsm)
      dwmass=dwnorm*dwsm
      dwmnbi=mass(p)*dwmass
      dvx=pvel(1)-veltpos(p,1)
      dvy=pvel(2)-veltpos(p,2)
      dvz=pvel(3)-veltpos(p,3)
      vdotdr=dvx*dx+dvy*dy+dvz*dz 	    
      lhdivv=lhdivv-dwmnbi*vdotdr
      ldrhodh=ldrhodh-dr2*hsminv*dwmnbi
    
    
      curlvx=curlvx+dwmnbi*(dz*dvy-dy*dvz)
      curlvy=curlvy+dwmnbi*(dx*dvz-dz*dvx)
      curlvz=curlvz+dwmnbi*(dy*dvx-dx*dvy)
    endif
  enddo   
  lhdivv=h*lhdivv/lrho
  lhdivv=lhdivv/(1+h*ldrhodh/3/(rhomin+lrho)) ! f_i correction 
  lhcurlv=h*SQRT(curlvx**2+curlvy**2+curlvz**2)/lrho
  lcsound=lcsound/lrho 
  mv=mv/lrho
  mv2=mv2/lrho
   
  lvdisp=MAX(0.,(mv2-mv(1)**2-mv(2)**2-mv(3)**2))
  lvdisp=SQRT(lvdisp)
  lvbulk=SQRT(sum((pvel-mv)**2))
end subroutine

function bhmdot(bhmass,ppos,pvel,h,nneigh)
  include 'globals.h'
  real :: h,bhmass,hneigh,ppos(3),pvel(3), &
          lrho,ldrhodh,lhdivv,lhcurlv,lvdisp,lvbulk,lcsound, &
          bhmdot,mdot,mdot_edd,mdot_bondi,mdot_divv,ts,hbondi, &
          drhodt,reff
  integer :: i,j,nneigh
! real,parameter :: salpetertime=450. ! (Myr) 
! real, parameter :: bh_rad_eff=0.1
! integer, parameter :: bh_method=2  
  if(nneigh.eq.0) return   
  call skinnydensity_bh(nneigh,h,ppos,pvel,lrho,ldrhodh,lhdivv, &
    lhcurlv,lvdisp,lvbulk,lcsound)

! old mdot recipe
  mdot_bondi=4*pi*bhmass**2*lrho/(lcsound**2+lvbulk**2)**(3/2.) 
! new mdot recipe
  drhodt=-lrho*lhdivv/h
  hbondi=bhmass/lcsound**2
  reff=min(h,hbondi)
  if(drhodt.GT.0) then
    mdot_divv=4/3.*pi*reff**3*drhodt
  else
    mdot_divv=0.
  endif 
! Eddington
  ts=salpetertime*bh_rad_eff*year*1.e6/timescale
  mdot_edd=bhmass/ts

  select case( bh_method)
  case(1)
    mdot=mdot_bondi
    bhmdot=mdot
  case(2)
    mdot=mdot_divv
    bhmdot=mdot
  case(3)
    mdot=mdot_bondi
    bhmdot=min(mdot_edd,mdot)
  case(4)
    mdot=mdot_divv
    bhmdot=min(mdot_edd,mdot)
  case default
    mdot=mdot_edd
    bhmdot=mdot
  end select

  write(bhfilenr,'(20g14.6)') tpos,reff,lrho,lcsound,lvbulk,drhodt, &
    lcsound**2*mumhkg1/gamma/gamma1,mdot_bondi,mdot_divv,mdot_edd,bhmdot, &
    bhmass,h,ppos(1),ppos(2),ppos(3),ldrhodh,pvel(1),pvel(2),pvel(3)

end function

subroutine blackholes
  include 'globals.h'
  integer :: i,j,p,nneigh,ii,jj,k
  real :: mdot,bhmass,ppos(3),pvel(3),h,hneigh,dmass,dp(3)
  real :: bhmdot,d1,d2
  real,allocatable :: workvect(:)
  allocate(workvect(nsph))
 
  if(nbh.LE.0) return
 
  workvect(1:nsph)=-1.
  bodlist(1:nsph)=0

  open(unit=bhfilenr,file=TRIM(outputfile)//'.bh', &
   status='unknown',POSITION='APPEND')

  do i=1,nbh
    p=nbodies-nbh+i
    write(bhfilenr,'(i9,"  ")',advance='NO') nbexist(p)
  
    bhmass=mass(p)
    ppos=pos(p,1:3)
    pvel=vel(p,1:3)
    h=hsmooth(p)

    ii=0;jj=0;
    hneigh=0.;nneigh=0
    call gethsm(ppos,h,ii,jj)
    call search_d2(root,2*h,ppos,nneigh,srlist,tempvect)
    mdot=bhmdot(bhmass,ppos,pvel,h,nneigh)   

    hsmooth(p)=h   
    bhmass=mass(p)+mdot*(tpos-tbh)
    tfeedb(p)=mdot
 
    do k=1,nneigh
      if(tempvect(k).LT.4*h**2) then
        if(workvect(srlist(k)).LT.0.OR. &
           workvect(srlist(k)).GT.tempvect(k).OR. &
           (workvect(srlist(k)).EQ.tempvect(k).AND.p.GT.bodlist(srlist(k))) ) then
          workvect(srlist(k))=tempvect(k)
          bodlist(srlist(k))=p
        endif   
      endif
    enddo  
  enddo

  close(bhfilenr)

  if(tbh.EQ.tpos) then
    do i=1,nbh
      p=nbodies-nbh+i
      if(tfeedb(p).GT.0) tform(p)=tbh
    enddo
    return
  endif
 
  do i=1,nbh
    p=nbodies-nbh+i
    h=hsmooth(p)
    ppos=pos(p,1:3)
    hneigh=0.;nneigh=0
    call search_d2(root,2*h,ppos,nneigh,srlist,tempvect)
    call mrgrnk(nneigh,tempvect,templist)
    dmass=tfeedb(p)*(tpos-tbh)
    j=1
    dp=mass(p)*vel(p,1:3)
    do while(dmass.GT.0.AND.tempvect(templist(j)).LT.4*h**2.AND.j.LE.nneigh)
      k=srlist(templist(j))
      if(bodlist(k).LE.nbodies-nbh.OR.bodlist(k).GT.nbodies) then
        print*,nbodies,nbh,nneigh
        print*, k,templist(j),bodlist(k),workvect(k)
        call terror('BH id error')
      endif
      if(bodlist(k).EQ.p) then
        if(dmass.LT.mass(k)) then 
          dp=dp+dmass*vel(k,1:3)
          mass(k)=mass(k)-dmass
          dmass=0.
        else
          dmass=dmass-mass(k)
          dp=dp+mass(k)*vel(k,1:3)
          mass(k)=0.
      endif 
    endif
    j=j+1
    enddo    
    if(dmass.LT.0) call terror(' BH accr error')
!    if(dmass.GT.0) call terror(' BH accr mismatch')
    if(dmass.GT.tfeedb(p)*(tpos-tbh)) call terror(' BH accr mismatch')
    if(dmass.GT.0) then 
      print*,' accretion warning!'
      tfeedb(p)=tfeedb(p)-dmass/(tpos-tbh)
    endif
!    mass(p)=mass(p)+tfeedb(p)*(tpos-tbh)-dmass  
    mass(p)=mass(p)+tfeedb(p)*(tpos-tbh)   
    vel(p,1:3)=dp/mass(p)
    if(tfeedb(p).GT.0) tform(p)=tbh
  enddo
  tbh=tpos
  deallocate(workvect)
end subroutine

subroutine bhmergers
  include 'globals.h'
  integer :: p,i,j,k,nneigh
  real :: h2,dr,ppos(3),dv2,pvel(3),pacc(4),mp,bht,bheps,mdot

  call bhtree
  
  do i=1,nbh
    p=nbodies-nbh+i
    ppos=pos(p,1:3)
    pvel=vel(p,1:3)
    pacc=acc(p,1:4)
    mp=mass(p)
    bht=tform(p)
    bheps=hsmooth(p)
    mdot=tfeedb(p)
    h2=2*hsmooth(p)

! do search 
    nneigh=0
    call search(root,h2,ppos,nneigh,srlist)
    if(nneigh.GT.0) then
!sort?
      do j=1,nneigh
        k=srlist(j)
        if(p.NE.k) then
          if((mp.GT.mass(k).OR.(mp.EQ.mass(k).AND.p.GT.k)).AND. &
                                                  mp.GT.tiny) then
            dr=sqrt(sum((ppos-pos(k,1:3))**2)) 
            dv2=sum((pvel-vel(k,1:3))**2)     
            if(2*(mp+mass(k))/dr.GT.dv2) then       
              ppos=(ppos*mp+pos(k,1:3)*mass(k))/(mp+mass(k))
              pvel=(pvel*mp+vel(k,1:3)*mass(k))/(mp+mass(k))
              pacc=(pacc*mp+acc(k,1:4)*mass(k))/(mp+mass(k))
              mp=mp+mass(k)
              bht=MAX(bht,tform(k))
              bheps=MAX(bheps,hsmooth(k))
              mdot=mdot+tfeedb(k)

              mass(k)=0.           
            endif
          endif
        endif
      enddo
      pos(p,1:3)=ppos
      vel(p,1:3)=pvel
      acc(p,1:4)=pacc
      mass(p)=mp
      tform(p)=bht
      hsmooth(p)=bheps
      tfeedb(p)=mdot
    endif  
  enddo
end subroutine
