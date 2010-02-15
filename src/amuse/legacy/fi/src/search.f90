!
! basic search routines for range search 
! (with and without distance list)
!

subroutine search_d2(node,h,ppos,npnear,nblist,d2list)
  include 'globals.h'
  integer, intent(in):: node
  real, intent(in) :: h,ppos(3)
  integer, intent(inout)::npnear,nblist(*)
  real, intent(inout) :: d2list(*) 
  integer :: i,daughter,nstack,cnode,cstack(nsubcell*81)
  real :: dist2,xnode,ynode,znode,dx,dy,dz,hsearch2

  if(periodic.AND.h.GT.rsize/4) call terror('search h>rsize/4')
  hsearch2=h**2
  cstack(1)=node
  nstack=1
  do while(nstack.GT.0.AND.nstack.LE.nsubcell*80)
    daughter=cstack(nstack)
    nstack=nstack-1
    if(daughter.LT.nbods1) then              
      if(daughter.GT.0) then 
        dx=ppos(1)-pos(daughter,1)
        dy=ppos(2)-pos(daughter,2)
        dz=ppos(3)-pos(daughter,3)
        if(dx.GE.hboxsize.AND.periodic) dx=dx-pboxsize
        if(dx.LT.-hboxsize.AND.periodic) dx=dx+pboxsize
        if(dy.GE.hboxsize.AND.periodic) dy=dy-pboxsize
        if(dy.LT.-hboxsize.AND.periodic) dy=dy+pboxsize
        if(dz.GE.hboxsize.AND.periodic) dz=dz-pboxsize
        if(dz.LT.-hboxsize.AND.periodic) dz=dz+pboxsize
        dist2=dx**2+dy**2+dz**2
        if(dist2.LT.hsearch2) then
          npnear=npnear+1
          nblist(npnear)=daughter
          d2list(npnear)=dist2
        endif
      endif 
    else
      xnode=bottom(daughter,1)+0.5*cellsize(daughter)
      ynode=bottom(daughter,2)+0.5*cellsize(daughter)
      znode=bottom(daughter,3)+0.5*cellsize(daughter)
      dx=ppos(1)-xnode
      dy=ppos(2)-ynode
      dz=ppos(3)-znode      
      if(dx.GE.hboxsize.AND.periodic) xnode=xnode+pboxsize
      if(dx.LT.-hboxsize.AND.periodic) xnode=xnode-pboxsize
      if(dy.GE.hboxsize.AND.periodic) ynode=ynode+pboxsize
      if(dy.LT.-hboxsize.AND.periodic) ynode=ynode-pboxsize
      if(dz.GE.hboxsize.AND.periodic) znode=znode+pboxsize
      if(dz.LT.-hboxsize.AND.periodic) znode=znode-pboxsize
      dist2=0
      if(ppos(1).LT.xnode-0.5*cellsize(daughter)) then
        dist2=dist2+(ppos(1)-xnode+0.5*cellsize(daughter))**2
      else
        if(ppos(1).GT.xnode+0.5*cellsize(daughter)) then
          dist2=dist2+(ppos(1)-xnode-0.5*cellsize(daughter))**2
        endif
      endif
      if(ppos(2).LT.ynode-0.5*cellsize(daughter)) then
        dist2=dist2+(ppos(2)-ynode+0.5*cellsize(daughter))**2
      else
        if(ppos(2).GT.ynode+0.5*cellsize(daughter)) then
          dist2=dist2+(ppos(2)-ynode-0.5*cellsize(daughter))**2
        endif
      endif
      if(ppos(3).LT.znode-0.5*cellsize(daughter)) then
        dist2=dist2+(ppos(3)-znode+0.5*cellsize(daughter))**2
      else
        if(ppos(3).GT.znode+0.5*cellsize(daughter)) then
          dist2=dist2+(ppos(3)-znode-0.5*cellsize(daughter))**2
        endif
      endif 
      if(dist2.LT.hsearch2) then
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
    call terror(" search_d2 error")
!$omp end critical
  endif            
end subroutine search_d2

subroutine search(node,h,ppos,npnear,nblist)
  include 'globals.h'
  integer, intent(in):: node
  real, intent(in) :: h,ppos(3)
  integer, intent(inout)::npnear,nblist(*)
  integer :: i,daughter,nstack,cnode,cstack(nsubcell*81)
  real :: dist2,xnode,ynode,znode,dx,dy,dz,hsearch2
  if(periodic.AND.h.GT.rsize/4) call terror('search h>rsize/4')
  hsearch2=h**2
  cstack(1)=node
  nstack=1
  do while(nstack.GT.0.AND.nstack.LE.nsubcell*80)
    daughter=cstack(nstack)
    nstack=nstack-1
    if(daughter.LT.nbods1) then              
      if(daughter.GT.0) then 
        dx=ppos(1)-pos(daughter,1)
        dy=ppos(2)-pos(daughter,2)
        dz=ppos(3)-pos(daughter,3)
        if(dx.GE.hboxsize.AND.periodic) dx=dx-pboxsize
        if(dx.LT.-hboxsize.AND.periodic) dx=dx+pboxsize
        if(dy.GE.hboxsize.AND.periodic) dy=dy-pboxsize
        if(dy.LT.-hboxsize.AND.periodic) dy=dy+pboxsize
        if(dz.GE.hboxsize.AND.periodic) dz=dz-pboxsize
        if(dz.LT.-hboxsize.AND.periodic) dz=dz+pboxsize
        dist2=dx**2+dy**2+dz**2
        if(dist2.LT.hsearch2) then
          npnear=npnear+1
          nblist(npnear)=daughter
        endif
      endif 
    else
      xnode=bottom(daughter,1)+0.5*cellsize(daughter)
      ynode=bottom(daughter,2)+0.5*cellsize(daughter)
      znode=bottom(daughter,3)+0.5*cellsize(daughter)
      dx=ppos(1)-xnode
      dy=ppos(2)-ynode
      dz=ppos(3)-znode      
      if(dx.GE.hboxsize.AND.periodic) xnode=xnode+pboxsize
      if(dx.LT.-hboxsize.AND.periodic) xnode=xnode-pboxsize
      if(dy.GE.hboxsize.AND.periodic) ynode=ynode+pboxsize
      if(dy.LT.-hboxsize.AND.periodic) ynode=ynode-pboxsize
      if(dz.GE.hboxsize.AND.periodic) znode=znode+pboxsize
      if(dz.LT.-hboxsize.AND.periodic) znode=znode-pboxsize
      dist2=0
      if(ppos(1).LT.xnode-0.5*cellsize(daughter)) then
        dist2=dist2+(ppos(1)-xnode+0.5*cellsize(daughter))**2
      else
        if(ppos(1).GT.xnode+0.5*cellsize(daughter)) then
          dist2=dist2+(ppos(1)-xnode-0.5*cellsize(daughter))**2
        endif
      endif
      if(ppos(2).LT.ynode-0.5*cellsize(daughter)) then
        dist2=dist2+(ppos(2)-ynode+0.5*cellsize(daughter))**2
      else
        if(ppos(2).GT.ynode+0.5*cellsize(daughter)) then
          dist2=dist2+(ppos(2)-ynode-0.5*cellsize(daughter))**2
        endif
      endif
      if(ppos(3).LT.znode-0.5*cellsize(daughter)) then
        dist2=dist2+(ppos(3)-znode+0.5*cellsize(daughter))**2
      else
        if(ppos(3).GT.znode+0.5*cellsize(daughter)) then
          dist2=dist2+(ppos(3)-znode-0.5*cellsize(daughter))**2
        endif
      endif 
      if(dist2.LT.hsearch2) then
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
    call terror(" search error")
!$omp end critical
  endif            
end subroutine search

!
! basic search routines for complementary range search 
! (with and without distance list)
!

subroutine comsearch_d2(node,h,delta,ppos,npnear,nblist,d2list)
  include 'globals.h'
  integer, intent(in):: node
  real, intent(in) :: h,ppos(3),delta
  integer, intent(inout)::npnear,nblist(*)
  real,intent(inout) ::  d2list(*)
  integer :: i,daughter,nstack,cnode,cstack(nsubcell*81)
  real :: dist2,xnode,ynode,znode,dx,dy,dz,hsearch2

  if(periodic.AND. &
      (2*h.GT.rsize/4.OR.(2*hsmooth(node)+delta).GT.rsize/4)) &
    call terror('comsearch h>rsize/4')
  hsearch2=h**2
  cstack(1)=node
  nstack=1
  do while(nstack.GT.0.AND.nstack.LE.nsubcell*80)
    daughter=cstack(nstack)
    nstack=nstack-1
    if(daughter.LT.nbods1) then              
      if(daughter.GT.0) then 
        dx=ppos(1)-pos(daughter,1)
        dy=ppos(2)-pos(daughter,2)
        dz=ppos(3)-pos(daughter,3)
        if(dx.GE.hboxsize.AND.periodic) dx=dx-pboxsize
        if(dx.LT.-hboxsize.AND.periodic) dx=dx+pboxsize
        if(dy.GE.hboxsize.AND.periodic) dy=dy-pboxsize
        if(dy.LT.-hboxsize.AND.periodic) dy=dy+pboxsize
        if(dz.GE.hboxsize.AND.periodic) dz=dz-pboxsize
        if(dz.LT.-hboxsize.AND.periodic) dz=dz+pboxsize
        dist2=dx**2+dy**2+dz**2
        if(dist2.LT.hsearch2.OR. &
                   dist2.LT.(2*hsmooth(daughter)+delta)**2) then        
          npnear=npnear+1
          nblist(npnear)=daughter
          d2list(npnear)=dist2
        endif
      endif  
    else
      xnode=bottom(daughter,1)+0.5*cellsize(daughter)
      ynode=bottom(daughter,2)+0.5*cellsize(daughter)
      znode=bottom(daughter,3)+0.5*cellsize(daughter)
      dx=ppos(1)-xnode
      dy=ppos(2)-ynode
      dz=ppos(3)-znode      
      if(dx.GE.hboxsize.AND.periodic) xnode=xnode+pboxsize
      if(dx.LT.-hboxsize.AND.periodic) xnode=xnode-pboxsize
      if(dy.GE.hboxsize.AND.periodic) ynode=ynode+pboxsize
      if(dy.LT.-hboxsize.AND.periodic) ynode=ynode-pboxsize
      if(dz.GE.hboxsize.AND.periodic) znode=znode+pboxsize
      if(dz.LT.-hboxsize.AND.periodic) znode=znode-pboxsize
      dist2=0
      if(ppos(1).LT.xnode-0.5*cellsize(daughter)) then
        dist2=dist2+(ppos(1)-xnode+0.5*cellsize(daughter))**2
      else
        if(ppos(1).GT.xnode+0.5*cellsize(daughter)) then
          dist2=dist2+(ppos(1)-xnode-0.5*cellsize(daughter))**2
        endif
      endif
      if(ppos(2).LT.ynode-0.5*cellsize(daughter)) then
        dist2=dist2+(ppos(2)-ynode+0.5*cellsize(daughter))**2
      else
        if(ppos(2).GT.ynode+0.5*cellsize(daughter)) then
          dist2=dist2+(ppos(2)-ynode-0.5*cellsize(daughter))**2
        endif
      endif
      if(ppos(3).LT.znode-0.5*cellsize(daughter)) then
        dist2=dist2+(ppos(3)-znode+0.5*cellsize(daughter))**2
      else
        if(ppos(3).GT.znode+0.5*cellsize(daughter)) then
          dist2=dist2+(ppos(3)-znode-0.5*cellsize(daughter))**2
        endif
      endif 
      if(dist2.LT.hsearch2.OR. &
         dist2.LT.(2*hsmooth(daughter)+delta)**2) then
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
    call terror(" comsearch_d2 error")
!$omp end critical
  endif            
end subroutine comsearch_d2

subroutine comsearch(node,h,delta,ppos,npnear,nblist)
  include 'globals.h'
  integer, intent(in):: node
  real, intent(in) :: h,ppos(3),delta
  integer, intent(inout)::npnear,nblist(*)
  integer :: i,daughter,nstack,cnode,cstack(nsubcell*81)
  real :: dist2,xnode,ynode,znode,dx,dy,dz,hsearch2
  
  if(periodic.AND. &
      (2*h.GT.rsize/4.OR.(2*hsmooth(node)+delta).GT.rsize/4)) &
    call terror('comsearch h>rsize/4')
  hsearch2=h**2
  cstack(1)=node
  nstack=1
  do while(nstack.GT.0.AND.nstack.LE.nsubcell*80)
    daughter=cstack(nstack)
    nstack=nstack-1
    if(daughter.LT.nbods1) then              
      if(daughter.GT.0) then 
        dx=ppos(1)-pos(daughter,1)
        dy=ppos(2)-pos(daughter,2)
        dz=ppos(3)-pos(daughter,3)
        if(dx.GE.hboxsize.AND.periodic) dx=dx-pboxsize
        if(dx.LT.-hboxsize.AND.periodic) dx=dx+pboxsize
        if(dy.GE.hboxsize.AND.periodic) dy=dy-pboxsize
        if(dy.LT.-hboxsize.AND.periodic) dy=dy+pboxsize
        if(dz.GE.hboxsize.AND.periodic) dz=dz-pboxsize
        if(dz.LT.-hboxsize.AND.periodic) dz=dz+pboxsize
        dist2=dx**2+dy**2+dz**2
        if(dist2.LT.hsearch2.OR. &
             dist2.LT.(2*hsmooth(daughter)+delta)**2) then        
          npnear=npnear+1
          nblist(npnear)=daughter
        endif
      endif  
    else

      xnode=bottom(daughter,1)+0.5*cellsize(daughter)
      ynode=bottom(daughter,2)+0.5*cellsize(daughter)
      znode=bottom(daughter,3)+0.5*cellsize(daughter)
      dx=ppos(1)-xnode
      dy=ppos(2)-ynode
      dz=ppos(3)-znode      
      if(dx.GE.hboxsize.AND.periodic) xnode=xnode+pboxsize
      if(dx.LT.-hboxsize.AND.periodic) xnode=xnode-pboxsize
      if(dy.GE.hboxsize.AND.periodic) ynode=ynode+pboxsize
      if(dy.LT.-hboxsize.AND.periodic) ynode=ynode-pboxsize
      if(dz.GE.hboxsize.AND.periodic) znode=znode+pboxsize
      if(dz.LT.-hboxsize.AND.periodic) znode=znode-pboxsize
      dist2=0
      if(ppos(1).LT.xnode-0.5*cellsize(daughter)) then
        dist2=dist2+(ppos(1)-xnode+0.5*cellsize(daughter))**2
      else
        if(ppos(1).GT.xnode+0.5*cellsize(daughter)) then
          dist2=dist2+(ppos(1)-xnode-0.5*cellsize(daughter))**2
        endif
      endif
      if(ppos(2).LT.ynode-0.5*cellsize(daughter)) then
        dist2=dist2+(ppos(2)-ynode+0.5*cellsize(daughter))**2
      else
        if(ppos(2).GT.ynode+0.5*cellsize(daughter)) then
          dist2=dist2+(ppos(2)-ynode-0.5*cellsize(daughter))**2
        endif
      endif
      if(ppos(3).LT.znode-0.5*cellsize(daughter)) then
        dist2=dist2+(ppos(3)-znode+0.5*cellsize(daughter))**2
      else
        if(ppos(3).GT.znode+0.5*cellsize(daughter)) then
          dist2=dist2+(ppos(3)-znode-0.5*cellsize(daughter))**2
        endif
      endif 
      if(dist2.LT.hsearch2.OR. &
           dist2.LT.(2*hsmooth(daughter)+delta)**2) then
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
    call terror(" comsearch error")
!$omp end critical
  endif            
end subroutine comsearch

!
! particle-index search routines for range search 
! (with and without distance list)
!
subroutine psearch(node,p,npnear,nblist)
  include 'globals.h'
  integer, intent(in):: node,p
  real :: h,ppos(3)
  integer, intent(inout)::npnear,nblist(*)
  h=2*hsmooth(p)
  ppos=pos(p,1:3)
  call search(node,h,ppos,npnear,nblist)
end subroutine

subroutine psearch_d2(node,p,npnear,nblist,d2list)
  include 'globals.h'
  integer, intent(in):: node,p
  real :: h,ppos(3)
  integer, intent(inout)::npnear,nblist(*)
  real  :: d2list(*)
  h=2*hsmooth(p)
  ppos=pos(p,1:3)
  call search_d2(node,h,ppos,npnear,nblist,d2list)     
end subroutine

!
! particle-index search routines for complementary range search 
! (with and without distance list)
!

subroutine pcomsearch(node,p,npnear,nblist)
  include 'globals.h'
  integer, intent(in):: node,p
  real :: h,ppos(3), delta
  integer, intent(inout)::npnear,nblist(*)
  h=2*hsmooth(p)
  ppos=pos(p,1:3)
  call comsearch(node,h,0.,ppos,npnear,nblist)
end subroutine

subroutine pcomsearch_d2(node,p,npnear,nblist,d2list)
  include 'globals.h'
  integer, intent(in):: node,p
  real :: h,ppos(3),delta
  integer, intent(inout)::npnear,nblist(*)
  real :: d2list(*)
  h=2*hsmooth(p)
  ppos=pos(p,1:3)
  delta=0.
  call comsearch_d2(node,h,delta,ppos,npnear,nblist,d2list)
end subroutine

!
! particle-index search routines for conditional (ie reusable) range search 
! (with and without distance list)
!

subroutine pcond_srch(node,p,nneigh,nblist)
  include 'globals.h'
  real h,h2,ppos(3)
  integer nneigh,p,node,nblist(*)
  h=hsmooth(p)
  ppos=pos(p,1:3)
  call cond_srch(ppos,h,node,nneigh,nblist)
end subroutine

subroutine pcond_srch_d2(node,p,nneigh,nblist,d2list)
  include 'globals.h'
  real h,h2,ppos(3)
  integer nneigh,p,node,nblist(*)
  real :: d2list(*)
  h=hsmooth(p)
  ppos=pos(p,1:3)
  call cond_srch_d2(ppos,h,node,nneigh,nblist,d2list)
end subroutine

!
! driver search routines for conditional (ie reusable) range search 
! (with and without distance list)
!

subroutine cond_srch(ppos,h,node,nneigh,nblist)
  include 'globals.h'
  real h,ppos(3),dist
  integer nneigh,node,nblist(*)
  ncalls=ncalls+1
  if(searchreuse.EQ.0.OR.reuseflag.EQ.0) then
    searchpos=ppos
    searchh=h
    searchn=0
  else
    dist=sqrt(sum((ppos-searchpos)**2))
    if(h+dist/2.LE.searchh.AND.h.GT.searchmagic*searchh) then
      nneigh=searchn
      return
    endif
    searchn=0
    searchpos=ppos
    searchh=h
  endif
  nsearches=nsearches+1
  call search(node,2*searchh,searchpos,searchn,nblist)
  nneigh=searchn
  searchreuse=1
end subroutine

subroutine cond_srch_d2(ppos,h,node,nneigh,nblist,d2list)
  include 'globals.h'
  real h,ppos(3),dist
  integer nneigh,node,nblist(*)
  real :: d2list(*)
  ncalls=ncalls+1
  if(searchreuse.EQ.0.OR.reuseflag.EQ.0) then
    searchpos=ppos
    searchh=h
    searchn=0
  else
    dist=sqrt(sum((ppos-searchpos)**2))
    if(h+dist/2.LE.searchh.AND.h.GT.searchmagic*searchh) then
      nneigh=searchn
      return
    endif
    searchn=0
    searchpos=ppos
    searchh=h
  endif
  nsearches=nsearches+1
  call search_d2(node,2*searchh,searchpos,searchn,nblist,d2list)
  nneigh=searchn
  searchreuse=1
end subroutine

!
! particle-index search routines for conditional (ie reusable) 
! and complementary range search 
! (with and without distance list)
!

subroutine pcond_comsrch(node,p,nneigh,nblist)
  include 'globals.h'
  real h,h2,ppos(3),delta
  integer nneigh,p,node,nblist(*)
  h=hsmooth(p)
  ppos=pos(p,1:3)
  delta=0.
  call cond_comsrch(ppos,h,delta,node,nneigh,nblist)
end subroutine

subroutine pcond_comsrch_d2(node,p,nneigh,nblist,d2list)
  include 'globals.h'
  real h,h2,ppos(3),delta
  integer nneigh,p,node,nblist(*)
  real :: d2list(*)
  h=hsmooth(p)
  ppos=pos(p,1:3)
  delta=0.
  call cond_comsrch_d2(ppos,h,delta,node,nneigh,nblist,d2list)
end subroutine

!
! driver search routines for conditional (ie reusable)
!  and complementary range search 
! (with and without distance list)
!
subroutine cond_comsrch(ppos,h,delta,node,nneigh,nblist)
  include 'globals.h'
  real h,ppos(3),dist,delta
  integer nneigh,node,nblist(*)
  ncalls=ncalls+1
  if(searchreuse.EQ.0.OR.reuseflag.EQ.0) then
    searchpos=ppos
    searchh=h
    searchn=0
    searchdelta=delta
  else
    dist=sqrt(sum((ppos-searchpos)**2))
    if(h+dist/2.LE.searchh.AND. &
        delta+dist.LE.searchdelta.AND. &
         h.GT.searchmagic*searchh) then
      nneigh=searchn
      return
    endif
    searchn=0
    searchpos=ppos
    searchh=h 
    searchdelta=delta
  endif
  nsearches=nsearches+1
  call comsearch(node,2*searchh,searchdelta,searchpos,searchn,nblist)
  nneigh=searchn
  searchreuse=1
end subroutine

subroutine cond_comsrch_d2(ppos,h,delta,node,nneigh,nblist,d2list)
  include 'globals.h'
  real h,ppos(3),dist,delta
  integer nneigh,node,nblist(*)
  real d2list(*)
  ncalls=ncalls+1
  if(searchreuse.EQ.0.OR.reuseflag.EQ.0) then
    searchpos=ppos
    searchh=h
    searchn=0
    searchdelta=delta
  else
    dist=sqrt(sum((ppos-searchpos)**2))
    if(h+dist/2.LE.searchh.AND. &
         delta+dist.LE.searchdelta.AND. &
         h.GT.searchmagic*searchh) then
      nneigh=searchn
      return
    endif
    searchn=0
    searchpos=ppos
    searchh=h 
    searchdelta=delta
  endif
  nsearches=nsearches+1
  call comsearch_d2(node,2*searchh,searchdelta, &
                       searchpos,searchn,nblist,d2list)
  nneigh=searchn
  searchreuse=1
end subroutine

!
! pre search routines
!  these try to determine a search such that multiple particles
!  can reuse the same search list

subroutine presearch(k,kmax,nbuf,buf,ntodo,todo)
  include 'globals.h'
  integer :: k,kmax,nbuf,buf(nbuf),todo(nbuf),ntodo,ib,p,q,nneigh 
  real :: dist2,spos(3),ppos(3),ph,qh,sh,qpos(3)
  ntodo=0
  if(buf(1).EQ.0) then
    p=pactive(k)
    ntodo=ntodo+1 
    todo(ntodo)=p
    if(hsmooth(p).GT.0) then
      ppos=pos(p,1:3); ph=hsmooth(p)
      spos=ppos; sh=ph
      do ib=1,MIN(nbuf-1,kmax-k)
        if(buf(ib+1).NE.0) cycle
        q=pactive(k+ib)
        qpos=pos(q,1:3); qh=hsmooth(q)
        dist2=sum((ppos-qpos)**2)
        if(dist2.LT.ph**2) then
          buf(ib+1)=1  
          ntodo=ntodo+1; todo(ntodo)=q
          spos=((ntodo-1)*spos)/ntodo+qpos/ntodo      
        endif
      enddo
      do ib=1,ntodo
        q=todo(ib)
        qpos=pos(q,1:3); qh=hsmooth(q)
        dist2=sum((spos-qpos)**2)         
        sh=MAX(sh,qh+sqrt(dist2)/2)
      enddo
      sh=sh*1.02
      searchreuse=0
      call cond_srch(spos,sh,root,nneigh,srlist)
    endif
  endif
  do ib=1,nbuf-1
    buf(ib)=buf(ib+1)
  enddo 
  buf(nbuf)=0
end subroutine

subroutine precomsearch(k,kmax,nbuf,buf,ntodo,todo)
  include 'globals.h'
  integer :: k,kmax,nbuf,buf(nbuf),todo(nbuf),ntodo,ib,p,q,nneigh 
  real :: dist2,spos(3),ppos(3),ph,qh,sh,qpos(3),delta
  ntodo=0
  if(buf(1).EQ.0) then
    p=pactive(k)
    ntodo=ntodo+1 
    todo(ntodo)=p
    if(hsmooth(p).GT.0) then
      ppos=pos(p,1:3); ph=hsmooth(p)
      spos=ppos; sh=ph; delta=0
      do ib=1,MIN(nbuf-1,kmax-k)
        if(buf(ib+1).NE.0) cycle
        q=pactive(k+ib)
        qpos=pos(q,1:3); qh=hsmooth(q)
        dist2=sum((ppos-qpos)**2)
        if(dist2.LT.ph**2) then
          buf(ib+1)=1  
          ntodo=ntodo+1; todo(ntodo)=q
          spos=((ntodo-1)*spos)/ntodo+qpos/ntodo      
        endif
      enddo
      do ib=1,ntodo
        q=todo(ib)
        qpos=pos(q,1:3); qh=hsmooth(q)
        dist2=sum((spos-qpos)**2)         
        sh=MAX(sh,qh+sqrt(dist2)/2)
        delta=MAX(delta,sqrt(dist2))
      enddo
      searchreuse=0
      call cond_comsrch(spos,sh,delta,root,nneigh,srlist)
    endif
  endif
  do ib=1,nbuf-1
    buf(ib)=buf(ib+1)
  enddo 
  buf(nbuf)=0
end subroutine


