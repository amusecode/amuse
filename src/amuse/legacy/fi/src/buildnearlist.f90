subroutine seteps(option)
  character(len=4) :: option
  real :: epsold

  call setbox(option)
  call loadtree5(option,0)
  call epscal(option)
end subroutine

subroutine epscal(option)
  include 'globals.h'
  character(len=4) :: option
  integer :: low,up,p,inear
  real :: position(3),hsearch,getlocalscale

  select case (option)
  case('all ')
    low=1;up=nbodies
  case('sph ')
    low=1;up=nsph
  case('star')
    low=nbodies-nstar+1;up=nbodies
  case('coll')
    low=nsph+1;up=nbodies
  case default
    low=1;up=nbodies
  end select

  nstot=0
  nsmin=nbodies
  nsmax=0

  if(low.GT.up) return
  if(up-low.LE.targetnn) then
    epsgrav(low:up)=rsize/2
    return
  endif    

!$omp parallel do private(p,position,hsearch,inear) &
!$omp reduction(+ : nstot) &
!$omp reduction(min : nsmin) &
!$omp reduction(max : nsmax)
  do p=low,up
    position=pos(p,1:3)
    if(.not.adaptive_eps.or.epsgrav(p).eq.0) then
      hsearch=getlocalscale(position)
    else
      hsearch=epsgrav(p)
    endif  
    call nnscalefix(position,hsearch,inear,targetnn,nn_tol,root)
    epsgrav(p)=hsearch
    nstot=nstot+inear
    nsmax=MAX(nsmax,inear)
    nsmin=MIN(nsmin,inear)
  enddo
  nsavg=nstot/(up-low+1)

  if(verbosity.GT.0) print*,'<epscal> < a > t:',nsmin,nsavg,nsmax,nstot

end subroutine

subroutine hsmcal
  include 'globals.h'
  character(len=4) :: option
  integer :: low,up,p,inear
  real :: position(3),hsearch,getlocalscale

  low=1;up=nsph

  nntot=0
  nnmin=nbodies
  nnmax=0

  if(low.GT.up) return

!$omp parallel do private(p,position,hsearch,inear) &
!$omp reduction(+ : nntot) &
!$omp reduction(min : nnmin) &
!$omp reduction(max : nnmax)
  do p=low,up
    position=pos(p,1:3)
    if(hupdatemethod.EQ.'none'.or.hsmooth(p).le.0) then
      hsearch=getlocalscale(position)
    else
      hsearch=hsmooth(p)
    endif  
    call nnscalefix(position,hsearch,inear,nsmooth,nsmtol,root)
    hsmooth(p)=hsearch
    nntot=nntot+inear
    nnmax=MAX(nnmax,inear)
    nnmin=MIN(nnmin,inear)
  enddo
  nnavg=nntot/(up-low+1)
  if(verbosity.GT.0) print*,'<hsmcal> < a > t:',nnmin,nnavg,nnmax,nntot

end subroutine


subroutine nnscalefix(position,hsearch,inear,ntarget,ntol,rt)
  include 'globals.h'
  real,dimension(3),intent(in) :: position
  real, intent(in) :: ntol
  real, intent(inout) :: hsearch
  integer, intent(in) :: ntarget,rt 
  integer,intent(inout) :: inear

  if(hsearch.le.0) call terror('nnscalefix error 1')
  inear=0
  do while(inear.LT.ntarget)
    if(hsearch.GT.rsize) call terror('nnscalefix error 2')
    inear=0
    call search_d2(rt,2*hsearch,position,inear,srlist,tempvect)
    if(inear.GE.(1.-ntol)*ntarget.AND.inear.LE.(1.+ntol)*ntarget) return
    hsearch=hsearch*1.3
  enddo
 
  call mrgrnk(inear,tempvect,srlist) ! srlist can be overwritten  
  hsearch=SQRT(tempvect(srlist(ntarget)))/2.
  inear=ntarget
 
end subroutine nnscalefix

function getlocalscale(spos) result(scale)
  include 'globals.h'
  real :: spos(3),lpos(3),upos(3),scale,dist2box
  integer :: daughter,subind,subcellindex2,id

  lpos=bottom(root,1:3)
  upos=bottom(root,1:3)+rsize
  scale=dist2box(spos,lpos,upos)
  if(scale.EQ.0) then
    daughter=root
    scale=rsize
    id=1
    do while(daughter.GT.nbodsmax) 
      scale=scale/2
      lpos=bottom(daughter,1:3)+cellsize(daughter)/2.
      subind=subcellindex2(lpos,spos)
      subind=peanovalue(subind,id)
      daughter=subp(daughter,subind)
      id=pshapefromp(subind,id)
    enddo
  else
    scale=sqrt(scale) 
  endif 
end function getlocalscale

function dist2box(spos,lpos,upos) result(dist2)
  real :: spos(3),lpos(3),upos(3),dist2
  integer :: k
 
  dist2=0
 
  do k=1,3
    if(spos(k).LT.lpos(k)) then
      dist2=dist2+(lpos(k)-spos(k))**2
    else
      if(spos(k).GT.upos(k)) then
        dist2=dist2+(upos(k)-spos(k))**2
      endif
    endif
  enddo
end function

