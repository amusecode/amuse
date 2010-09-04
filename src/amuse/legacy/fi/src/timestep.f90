function courant_timestep(p) result(dt)
! (eq 2.37/2.38 hk)
 include 'globals.h'
 integer :: p
 real :: dt
 dt=dtime
 if(p.GT.nsph) return
 dt=(hsmooth(p)*courant)/ &
   (csound(p)+1.2*(alpha*csound(p)+beta*mumaxdvh(p))+ABS(hsmdivv(p)))
 if(dt.le.0) then
   print*,p,dt,hsmooth(p),csound(p),mumaxdvh(p),hsmdivv(p)
   call terror('timestep error: courant')	  
 endif	 
end function

function sqrtacc_timestep(p) result(dt)
 include 'globals.h'
 integer :: p
 real :: dt
 dt=dtime
 if(acc(p,4).LE.0) return
 if(adaptive_eps) then
   dt=tstepcrit*sqrt(epsgrav(p)/acc(p,4))
 else
   dt=tstepcrit*sqrt(MIN(epsgrav(p),eps)/acc(p,4))		
 endif
 if(dt.le.0.) then
   print*,dt,epsgrav(p),acc(p,4),eps
   call terror('timestep error: sqrttstp')	  
 endif	 
end function

function acc_timestep(p) result(dt)
 include 'globals.h'
 integer :: p
 real :: dt
 dt=dtime
 if(acc(p,4).LE.0) return
 dt=tstpcr2/acc(p,4)
 if(dt.le.0.) then
   print*,dt,epsgrav(p),acc(p,4)
   call terror('timestep error: acc_tstp')	  
 endif	 
end function

function freeform_timestep(p) result(dt)
 include 'globals.h'
 integer :: p
 real :: dt,abvel
 dt=dtime
 abvel=sqrt(vel(p,1)**2+vel(p,2)**2+vel(p,3)**2)
 if(acc(p,4).LE.0.OR.abvel.LE.0) return
 dt=(abvel/freev)**freevexp*(acc(p,4)/freea)**freeaexp
 if(dt.le.0.) then
   print*,dt,abvel,acc(p,4)
   call terror('timestep error: freetstp')	  
 endif	 
end function

subroutine promote_bins(minppbin,npactive,templist,bodlist)
 integer :: npactive,minppbin,templist(npactive),bodlist(npactive)
 integer :: i,p,startbod,intt,ninbin
 character(len=70) mess
 
 if(minppbin.LE.1) return
 i=0
 p=1
 startbod=0
 do while(p.LE.npactive)
   if(MOD(p-1,minppbin).EQ.0) then
      intt=templist(bodlist(p))
      ninbin=0
      startbod=p
   else
      if(intt.LT.templist(bodlist(p))) then
	print*,intt,templist(bodlist(p))
     	call terror(' timestep inconsistency')
      endif
      if(intt.GT.templist(bodlist(p))) i=i+1
      templist(bodlist(p))=intt
   endif
   ninbin=ninbin+1
   p=p+1
 enddo
 if(ninbin.LE.minppbin/2.AND.startbod.NE.0) then
   if(startbod.EQ.1) call terror(' timestep shuffle error')
   templist(bodlist(startbod:npactive))=templist(bodlist(startbod-1))
 endif
 write(mess,'("<timestep> checked, promoted:", i9,i9)') npactive,i 
 call report(1, mess)
end subroutine 

function log2(i) ! integer log base 2 rounded down
 integer log2,i
 if(i.LE.0) return
 log2=30
 do while(IAND(2**log2,i).EQ.0)
  log2=log2-1
 enddo
end function 

subroutine refresh_itimestp
 include 'globals.h'
 integer :: p,pbin,npcheck
 real :: dt
 real :: courant_timestep,sqrtacc_timestep,acc_timestep,freeform_timestep
 integer :: log2

 npcheck=0
 do p=1,nbodies
   if(itimestp(p).GT.active_bin) then
      npcheck=npcheck+1

      dt=dtime
      if(usesph)   dt=MIN(dt,courant_timestep(p))
      if(sqrttstp)  dt=MIN(dt,sqrtacc_timestep(p))
      if(acc_tstp)  dt=MIN(dt,acc_timestep(p))
      if(freetstp) dt=MIN(dt,freeform_timestep(p))
      
      if(dt.GE.dtime) then
        tempvect(npcheck)=dtime
	templist(npcheck)=1
      else
        tempvect(npcheck)=dt
	if(dtime/tempvect(npcheck).GT.2**30) &
                   call terror('timestep way too small') 
	templist(npcheck)=INT(dtime/tempvect(npcheck))+1
      endif
   endif
 enddo 

 call mrgrnk(npcheck,tempvect,bodlist)
 
 if(npcheck.GE.nsmooth+1) then 
   if(templist(bodlist(nsmooth+1)).GT.max_tbin) & 
                 call terror(' timestep too small')
 endif
 
 if(minppbin.GT.1) then
  call promote_bins(minppbin,npcheck,templist,bodlist) 
 else
  if(verbosity.GT.0) print*,'<timestep> checked:',npcheck
 endif
 
 npcheck=0
 do p=1,nbodies
   if(itimestp(p).GT.active_bin) then
      npcheck=npcheck+1
      pbin=active_bin+1   
      do while(2**(pbin-1).LT.templist(npcheck).AND.2**(pbin-1).LT.max_tbin)
        pbin=pbin+1
      enddo
      itimestp(p)=max(itimestp(p)-1,pbin)
   endif   
 enddo         

end subroutine

subroutine timestep(itime)
 include 'globals.h'
 integer :: itime
 integer :: i,p,pbin,largest_bin
 real :: dt
 integer :: maxbin
 real :: courant_timestep,sqrtacc_timestep,acc_timestep,freeform_timestep
 integer :: log2
 
 maxbin=log2(2*max_tbin)

 call refresh_itimestp()

 largest_bin=maxval(itimestp(1:nbodies))
 
 itime=itime+2**(maxbin-largest_bin)

 active_bin=largest_bin
 do while(IAND(itime,2**(maxbin-active_bin)).EQ.0)
  active_bin=active_bin-1
 enddo

 tsteppos=dtime/2**largest_bin

 npactive=0
 nsphact=0

 if(active_bin.EQ.0) then
   if(itime.NE.2**maxbin) then    
      print*,itime,maxbin,2**maxbin
      call terror('itime inconsistent')
   endif
  else
   do p=1,nbodies
    if(itimestp(p).EQ.active_bin) then
      npactive=npactive+1
      if(p.LE.nsph.and.usesph) nsphact=nsphact+1
      pactive(npactive)=p
    endif
   enddo
 endif

 if(verbosity.GT.0) &
   print*, '<timestep> smallest, current,npactive:', largest_bin,active_bin,npactive

 if(active_bin.EQ.0.AND.verbosity.GT.0) print*, '<timestep> endstep'

 if(npactive.GT.0) then 
   call corrpos(otimestp,'sync')
   call corrpos(itimestp,'desync')
   otimestp(pactive(1:npactive))=itimestp(pactive(1:npactive))
 endif
 
end subroutine

