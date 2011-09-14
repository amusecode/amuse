subroutine pcond_fuvwalk(node,p,nterms)
   include 'globals.h'
   real h,h2,ppos(3),peps
   integer nterms,p,node
   h=hsmooth(p)
   ppos=pos(p,1:3)
   call cond_fuvwalk(node,ppos,0.,nterms)    
end subroutine

subroutine cond_fuvwalk(node,ppos,delta,nterms)
  include 'globals.h'
  real peps,ppos(3),delta,dist
  integer nterms,node
  ncalls=ncalls+1
  if(searchreuse.EQ.0.OR.reuseflag.EQ.0) then
    searchpos=ppos
    searchn=0
    searchdelta=delta
  else
    dist=sqrt(sum((ppos-searchpos)**2))
    if(delta+dist.LE.searchdelta) then
      nterms=searchn
      return
    endif
    searchpos=ppos
    searchn=0
    searchdelta=delta
  endif
  nsearches=nsearches+1
  call fuvwalk(node,searchpos,searchdelta,searchn)
  nterms=searchn
  searchreuse=1
end subroutine

subroutine pcond_treewalk(node,p,nterms)
  include 'globals.h'
  real h,h2,ppos(3),peps,acc4
  integer nterms,p,node
  h=hsmooth(p)
  ppos=pos(p,1:3)
  peps=epsgrav(p)
  acc4=acc(p,4)   
  call cond_treewalk(node,ppos,peps,acc4,0.,nterms)    
end subroutine

subroutine cond_treewalk(node,ppos,peps,acc4,delta,nterms)
  include 'globals.h'
  real peps,ppos(3),acc4,delta,dist
  integer nterms,node
  ncalls=ncalls+1
  if(searchreuse.EQ.0.OR.reuseflag.EQ.0) then
    searchpos=ppos
    searchh=peps
    searchn=0
    searchdelta=delta
    searchacc4=acc4   
  else
    dist=sqrt(sum((ppos-searchpos)**2))
    if(peps.LE.searchh.AND. &
          delta+dist.LE.searchdelta.AND. &
                          acc4.GE.searchacc4) then
      nterms=searchn
      return
    endif
    searchpos=ppos
    searchh=peps
    searchn=0
    searchdelta=delta
    searchacc4=acc4   
  endif
  nsearches=nsearches+1
  call treewalk(node,searchpos,searchh,searchacc4,searchdelta,searchn)
  nterms=searchn
  searchreuse=1
end subroutine

subroutine ptreewalk(node,p,nterms)
  include 'globals.h'
  integer, intent(in):: node,p
  integer, intent(inout)::nterms
  real :: ppos(3),peps,acc4
  ppos=pos(p,1:3)
  peps=epsgrav(p)
  acc4=acc(p,4)
  call treewalk(node,ppos,peps,acc4,0.,nterms)
end subroutine

subroutine treewalk(node,ppos,peps,acc4,delta,nterms)
  include 'globals.h'
  integer, intent(in):: node
  integer, intent(inout)::nterms
  real,intent(in) :: ppos(3),peps,acc4,delta
  if(periodic) then
    call terror('periodic treewalk broken?')
  endif
  if(acc4.LE.0.OR.(.not.gdgop)) then 
    if(usepm) then
      call pmtreewalk(node,ppos,delta,nterms)
    else
      if(usequad) then
        call quadtreewalk(node,ppos,peps,delta,nterms)
      else   
        call monotreewalk(node,ppos,delta,nterms)
      endif
    endif
  else
    if(usepm) then
      call gdgpmtreewalk(node,ppos,acc4,delta,nterms)
    else
      if(usequad) then
        call gdgquadtreewalk(node,ppos,peps,acc4,delta,nterms)
      else
        call gdgmonotreewalk(node,ppos,acc4,delta,nterms)
      endif
    endif   
  endif
end subroutine

subroutine pfuvwalk(node,p,nterms)
  include 'globals.h'
  integer, intent(in):: node,p
  integer, intent(inout)::nterms
  real :: ppos(3)
  ppos=pos(p,1:3)
  call fuvwalk(node,ppos,0.,nterms)
end subroutine

subroutine fuvwalk(node,ppos,delta,nterms)
  include 'globals.h'
  integer, intent(in):: node
  integer, intent(inout)::nterms
  real, intent(in) :: ppos(3),delta
  if(periodic) then
    call terror('periodic treewalk broken!')
  endif
  call monotreewalk(node,ppos,delta,nterms)
end subroutine

subroutine monotreewalk(node,ppos,delta,nterms)
  include 'globals.h'
  integer, intent(in):: node
  real,intent(in) :: ppos(3),delta 
  integer, intent(inout)::nterms
  integer :: daughter,nstack,cstack(nsubcell*81)
  real :: bh_tol2,dist2,dx,dy,dz,dist2cell 
  bh_tol2=bh_tol*bh_tol
  cstack(1)=node
  nstack=1
  do while(nstack.GT.0.AND.nstack.LE.nsubcell*80)
    daughter=cstack(nstack)
    nstack=nstack-1
    if(daughter.GT.0) then 
      if(daughter.LT.nbods1) then
        nterms=nterms+1
        bodlist(nterms)=daughter
      else  
        dx=ppos(1)-pos(daughter,1)
        dy=ppos(2)-pos(daughter,2)
        dz=ppos(3)-pos(daughter,3)
        dist2=dx**2+dy**2+dz**2
        if(bh_tol2*dist2.GT.(cellsize(daughter)+delta*bh_tol)**2) then
          dist2=dist2cell(ppos,bottom(daughter,1:3),cellsize(daughter))
          if(dist2.GT.(delta+0.1*cellsize(daughter))**2) then
            nterms=nterms+1
            bodlist(nterms)=daughter
            cycle
          endif
        endif  
        cstack(nstack+1)=subp(daughter,8)
        cstack(nstack+2)=subp(daughter,7)
        cstack(nstack+3)=subp(daughter,6)
        cstack(nstack+4)=subp(daughter,5)
        cstack(nstack+5)=subp(daughter,4)
        cstack(nstack+6)=subp(daughter,3)
        cstack(nstack+7)=subp(daughter,2)
        cstack(nstack+8)=subp(daughter,1)
        nstack=nstack+8
      endif 
    endif
  enddo
  if(nstack.NE.0) then
!$omp critical
    call terror(" treewalk error")
!$omp end critical
  endif            
end subroutine

subroutine quadtreewalk(node,ppos,peps,delta,nterms)
  include 'globals.h'
  integer, intent(in):: node
  real,intent(in) :: ppos(3),peps,delta 
  integer, intent(inout)::nterms
  integer :: daughter,nstack,cstack(nsubcell*81)
  real :: bh_tol2,dist2,dx,dy,dz,dist2cell 
  cstack(1)=node
  nstack=1
  bh_tol2=bh_tol*bh_tol
  do while(nstack.GT.0.AND.nstack.LE.nsubcell*80)
  daughter=cstack(nstack)
  nstack=nstack-1
    if(daughter.LT.nbods1) then
      if(daughter.GT.0) then 
        nterms=nterms+1
        bodlist(nterms)=daughter
      endif
    else  
      dx=ppos(1)-pos(daughter,1)
      dy=ppos(2)-pos(daughter,2)
      dz=ppos(3)-pos(daughter,3)
      dist2=dx**2+dy**2+dz**2
      if(bh_tol2*dist2.GT.(bh_tol*delta+cellsize(daughter))**2) then
        if(dist2.GT.(delta+peps+epsgrav(daughter))**2) then
          dist2=dist2cell(ppos,bottom(daughter,1:3),cellsize(daughter))
          if(dist2.GT.(delta+0.1*cellsize(daughter))**2) then
            nterms=nterms+1
            bodlist(nterms)=daughter
            cycle
          endif
        endif
      endif  
      cstack(nstack+1)=subp(daughter,8)
      cstack(nstack+2)=subp(daughter,7)
      cstack(nstack+3)=subp(daughter,6)
      cstack(nstack+4)=subp(daughter,5)
      cstack(nstack+5)=subp(daughter,4)
      cstack(nstack+6)=subp(daughter,3)
      cstack(nstack+7)=subp(daughter,2)
      cstack(nstack+8)=subp(daughter,1)
      nstack=nstack+8
    endif 
  enddo
  if(nstack.NE.0) then
!$omp critical
    call terror(" treewalk error")
!$omp end critical
  endif            
end subroutine

subroutine gdgmonotreewalk(node,ppos,acc4,delta,nterms)
  include 'globals.h'
  integer, intent(in):: node
  real,intent(in) :: ppos(3),acc4,delta 
  integer, intent(inout)::nterms
  integer :: daughter,nstack,cstack(nsubcell*81)
  real :: dist2,dx,dy,dz,dist2cell 
  cstack(1)=node
  nstack=1
  do while(nstack.GT.0.AND.nstack.LE.nsubcell*80)
    daughter=cstack(nstack)
    nstack=nstack-1
    if(daughter.LT.nbods1) then
      if(daughter.GT.0) then 
        nterms=nterms+1
        bodlist(nterms)=daughter
      endif
    else  
      dx=ppos(1)-pos(daughter,1)
      dy=ppos(2)-pos(daughter,2)
      dz=ppos(3)-pos(daughter,3)
      dist2=dx**2+dy**2+dz**2
      dist2=MAX((sqrt(dist2)-delta),0.)**2
      if(mass(daughter)*cellsize(daughter)**2.LT.gdgtol*acc4*dist2**2) then
        dist2=dist2cell(ppos,bottom(daughter,1:3),cellsize(daughter))
        if(dist2.GT.(delta+0.1*cellsize(daughter))**2) then
          nterms=nterms+1
          bodlist(nterms)=daughter
          cycle
        endif
      endif  
      cstack(nstack+1)=subp(daughter,8)
      cstack(nstack+2)=subp(daughter,7)
      cstack(nstack+3)=subp(daughter,6)
      cstack(nstack+4)=subp(daughter,5)
      cstack(nstack+5)=subp(daughter,4)
      cstack(nstack+6)=subp(daughter,3)
      cstack(nstack+7)=subp(daughter,2)
      cstack(nstack+8)=subp(daughter,1)
      nstack=nstack+8
    endif 
  enddo
  if(nstack.NE.0) then
!$omp critical
    call terror(" treewalk error")
!$omp end critical
  endif            
end subroutine

subroutine gdgquadtreewalk(node,ppos,peps,acc4,delta,nterms)
  include 'globals.h'
  integer, intent(in):: node
  real,intent(in) :: ppos(3),peps,acc4 ,delta
  integer, intent(inout)::nterms
  integer :: daughter,nstack,cstack(nsubcell*81)
  real :: dist2,dx,dy,dz,dist2cell 
  cstack(1)=node
  nstack=1
  do while(nstack.GT.0.AND.nstack.LE.nsubcell*80)
    daughter=cstack(nstack)
    nstack=nstack-1
    if(daughter.LT.nbods1) then
      if(daughter.GT.0) then 
        nterms=nterms+1
        bodlist(nterms)=daughter
      endif 
    else  
      dx=ppos(1)-pos(daughter,1)
      dy=ppos(2)-pos(daughter,2)
      dz=ppos(3)-pos(daughter,3)
      dist2=dx**2+dy**2+dz**2
      dist2=MAX(sqrt(dist2)-delta,0.)
      if(mass(daughter)*cellsize(daughter)**4.LT.gdgtol*acc4*dist2**6) then
        if(dist2.GT.(peps+epsgrav(daughter))) then
          dist2=dist2cell(ppos,bottom(daughter,1:3),cellsize(daughter))
          if(dist2.GT.(delta+0.1*cellsize(daughter))**2) then
            nterms=nterms+1
            bodlist(nterms)=daughter
            cycle
          endif
        endif
      endif  
      cstack(nstack+1)=subp(daughter,8)
      cstack(nstack+2)=subp(daughter,7)
      cstack(nstack+3)=subp(daughter,6)
      cstack(nstack+4)=subp(daughter,5)
      cstack(nstack+5)=subp(daughter,4)
      cstack(nstack+6)=subp(daughter,3)
      cstack(nstack+7)=subp(daughter,2)
      cstack(nstack+8)=subp(daughter,1)
      nstack=nstack+8
    endif 
  enddo
  if(nstack.NE.0) then
!$omp critical
    call terror(" treewalk error")
!$omp end critical
  endif            
end subroutine


subroutine pmtreewalk(node,ppos,delta,nterms)
  include 'globals.h'
  integer, intent(in):: node
  real,intent(in) :: ppos(3) ,delta
  integer, intent(inout)::nterms
  integer :: daughter,nstack,cstack(nsubcell*81)
  real :: bh_tol2,dist2,dx,dy,dz,dist2cell 
  cstack(1)=node
  nstack=1
  bh_tol2=bh_tol*bh_tol
  do while(nstack.GT.0.AND.nstack.LE.nsubcell*80)
    daughter=cstack(nstack)
    nstack=nstack-1
    if(daughter.LT.nbods1) then
      if(daughter.GT.0) then   
        nterms=nterms+1
        bodlist(nterms)=daughter
      endif
    else  
      dist2=dist2cell(ppos,bottom(daughter,1:3),cellsize(daughter))       
      if(dist2.GT.(rcut+delta)**2) cycle
      if(dist2.GT.(delta+0.1*cellsize(daughter))**2) then
        dx=ppos(1)-pos(daughter,1)
        dy=ppos(2)-pos(daughter,2)
        dz=ppos(3)-pos(daughter,3)
        dist2=dx**2+dy**2+dz**2
        if((bh_tol2*dist2).GT.(cellsize(daughter)+bh_tol*delta)**2) then
          if(dist2.LT.(rcut+delta)**2) then
            nterms=nterms+1
            bodlist(nterms)=daughter
          endif
          cycle
        endif 
      endif 
      cstack(nstack+1)=subp(daughter,8)
      cstack(nstack+2)=subp(daughter,7)
      cstack(nstack+3)=subp(daughter,6)
      cstack(nstack+4)=subp(daughter,5)
      cstack(nstack+5)=subp(daughter,4)
      cstack(nstack+6)=subp(daughter,3)
      cstack(nstack+7)=subp(daughter,2)
      cstack(nstack+8)=subp(daughter,1)
      nstack=nstack+8
    endif 
  enddo
  if(nstack.NE.0) then
!$omp critical
    call terror(" treewalk error")
!$omp end critical
  endif            
end subroutine

subroutine gdgpmtreewalk(node,ppos,acc4,delta,nterms)
  include 'globals.h'
  integer, intent(in):: node
  real,intent(in) :: ppos(3),acc4,delta 
  integer, intent(inout)::nterms
  integer :: daughter,nstack,cstack(nsubcell*81)
  real :: dist2,dx,dy,dz,dist2cell 
  cstack(1)=node
  nstack=1
  do while(nstack.GT.0.AND.nstack.LE.nsubcell*80)
    daughter=cstack(nstack)
    nstack=nstack-1
    if(daughter.LT.nbods1) then
      if(daughter.GT.0) then 
        nterms=nterms+1
        bodlist(nterms)=daughter
      endif
    else  
      dist2=dist2cell(ppos,bottom(daughter,1:3),cellsize(daughter))       
      if(dist2.GT.(rcut+delta)**2) cycle
      if(dist2.GT.(delta+0.1*cellsize(daughter))**2) then
        dx=ppos(1)-pos(daughter,1)
        dy=ppos(2)-pos(daughter,2)
        dz=ppos(3)-pos(daughter,3)
        dist2=dx**2+dy**2+dz**2
        dist2=MAX(sqrt(dist2)-delta,0.)**2
        if(mass(daughter)*cellsize(daughter)**2.LT.gdgtol*acc4*dist2**2) then
          if(dist2.LT.rcut2) then
            nterms=nterms+1
            bodlist(nterms)=daughter
          endif
          cycle
        endif 
      endif 
      cstack(nstack+1)=subp(daughter,8)
      cstack(nstack+2)=subp(daughter,7)
      cstack(nstack+3)=subp(daughter,6)
      cstack(nstack+4)=subp(daughter,5)
      cstack(nstack+5)=subp(daughter,4)
      cstack(nstack+6)=subp(daughter,3)
      cstack(nstack+7)=subp(daughter,2)
      cstack(nstack+8)=subp(daughter,1)
      nstack=nstack+8
    endif 
  enddo
  if(nstack.NE.0) then
!$omp critical
    call terror(" treewalk error")
!$omp end critical
  endif            
end subroutine

function dist2cell(spos,bot,csize) result(dist2)
  real :: spos(3),bot(3),csize,dist2
  integer :: k
  dist2=0 
  do k=1,3
    if(spos(k).LT.bot(k)) then
      dist2=dist2+(bot(k)-spos(k))**2
    else
      if(spos(k).GT.bot(k)+csize) then
        dist2=dist2+(spos(k)-(bot(k)+csize))**2
      endif
    endif
  enddo
end function

subroutine pretreewalk(k,kmax,nbuf,buf,ntodo,todo)
  include 'globals.h'
  integer :: k,kmax,nbuf,buf(nbuf),todo(nbuf),ntodo, &
    ib,p,q,nterms 
  real :: dist2,spos(3),ppos(3),qpos(3)
  real :: pacc4,qacc4,sacc4,peps,qeps,seps,delta,scale2,scale
  real :: getlocalscale
  ntodo=0
  if(buf(1).EQ.0) then
    p=pactive(k)
    ntodo=ntodo+1; todo(ntodo)=p
    if(.not.directsum) then
      ppos=pos(p,1:3); pacc4=acc(p,4); peps=epsgrav(p)
      scale=getlocalscale(ppos)/bh_tol
      scale=MIN(eps*3,scale)
      scale2=scale**2
!   if(gdgop.AND.pacc4.GT.0.AND.mass(p).GT.0) then
!    if(usequad) then
!      scale2=(mass(p)*scale2**2/gdgtol/pacc4)**(1./3.)      
!    else
!      scale2=sqrt(mass(p)*scale**2/gdgtol/pacc4)
!    endif
!   endif
      if(scale.LE.0.OR.scale.GT.rsize) call terror('pretreewalk error')
      spos=ppos; sacc4=pacc4;seps=peps;delta=0.
      do ib=1,MIN(nbuf-1,kmax-k)
        if(buf(ib+1).NE.0) cycle
        q=pactive(k+ib)
        qpos=pos(q,1:3)
        dist2=sum((ppos-qpos)**2)
        if(dist2.LT.2*scale2) then
          buf(ib+1)=1        
          ntodo=ntodo+1; todo(ntodo)=q
          spos=((ntodo-1)*spos)/ntodo+qpos/ntodo             
        endif
      enddo
      do ib=1,ntodo
        q=todo(ib)
        qpos=pos(q,1:3); qacc4=acc(q,4); qeps=epsgrav(q)
        seps=MAX(seps,qeps)
        dist2=sum((spos-qpos)**2)         
        delta=MAX(delta,sqrt(dist2))
        sacc4=MIN(sacc4,qacc4)
      enddo
      searchreuse=0
      call cond_treewalk(root,spos,seps,sacc4,delta,nterms)    
    endif
  endif
  do ib=1,nbuf-1
    buf(ib)=buf(ib+1)
  enddo 
  buf(nbuf)=0
end subroutine

subroutine prefuvwalk(k,kmax,nbuf,buf,ntodo,todo)
  include 'globals.h'
  integer :: k,kmax,nbuf,buf(nbuf),todo(nbuf),ntodo, &
    ib,p,q,nterms 
  real :: dist2,spos(3),ppos(3),qpos(3)
  real :: delta,scale2,scale
  ntodo=0
  if(buf(1).EQ.0) then
    p=pactive(k)
    ntodo=ntodo+1; todo(ntodo)=p
    if(.not.directsum) then
      ppos=pos(p,1:3)
      scale=epsgrav(p) ! =hsmooth(p) !!
      scale2=scale**2
      if(scale.LE.0.OR.scale.GT.rsize) call terror('prefuvwalk error')
      spos=ppos; delta=0.
      do ib=1,MIN(nbuf-1,kmax-k)
        if(buf(ib+1).NE.0) cycle
        q=pactive(k+ib)
        qpos=pos(q,1:3)
        dist2=sum((ppos-qpos)**2)
        if(dist2.LT.scale2) then
          buf(ib+1)=1        
          ntodo=ntodo+1; todo(ntodo)=q
          spos=((ntodo-1)*spos)/ntodo+qpos/ntodo             
        endif
      enddo
      do ib=1,ntodo
        q=todo(ib)
        qpos=pos(q,1:3);
        dist2=sum((spos-qpos)**2)         
        delta=MAX(delta,sqrt(dist2))
      enddo
      searchreuse=0
      call cond_fuvwalk(root,spos,delta,nterms)    
    endif
  endif
  do ib=1,nbuf-1
    buf(ib)=buf(ib+1)
  enddo 
  buf(nbuf)=0
end subroutine

