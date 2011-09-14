subroutine maketree
  include 'globals.h'
  integer i

  if(treestatecount.EQ.ppropcount+pordercount) return
  call setbox('all ')
  call loadtree5('all ',0)
  incellsg=incells
  treestatecount=ppropcount+pordercount
  sphtreecount=sphtreecount-1
end subroutine
	
subroutine makesphtree
  include 'globals.h'

  if(sphtreecount.EQ.ppropcount+pordercount) return
  if(nsph.EQ.nbodies.AND. &
     treestatecount.EQ.ppropcount+pordercount) then
     call tree_reduction(root,incells,'sph ')
  else
    call setbox('sph ')
    call loadtree5('sph ',0)
  endif
  sphtreecount=ppropcount+pordercount
  treestatecount=treestatecount-1
end subroutine

subroutine startree
  include 'globals.h'
  
  call setbox('star')
  call loadtree5('star',0)
  incellsg=incells
end subroutine

subroutine bhtree
  call setbox('bh  ')
  call loadtree5('bh  ',0)
end subroutine
        
subroutine setbox(option)
  include 'globals.h'
  character(len=4) :: option
  integer :: low,up,k
  real posmax,posmin,absmax
  
  if(periodic) then
    rsize=pboxsize
    rmin(1)=-0.5*pboxsize
    rmin(2)=-0.5*pboxsize
    rmin(3)=-0.5*pboxsize
    return
  endif 
 
  select case (option)
  case('all ')
   low=1;up=nbodies
  case('sph ')
   low=1;up=nsph
  case('star')
   low=nbodies-nstar+1;up=nbodies
  case('coll')
   low=nsph+1;up=nbodies
  case('bh  ')
   low=nbodies-nbh+1;up=nbodies
  case default
   low=1;up=nbodies
  end select
 
  posmax=MAXVAL(pos(low:up,1:3))
  posmin=MINVAL(pos(low:up,1:3))
  absmax=1.001*MAX(ABS(posmax),ABS(posmin))
  k=CEILING(log(absmax)/log(2.))+1
 
  rsize=2.**k 
  rmin=-0.5*rsize
end subroutine setbox

! cell sort2: make indexing in order of cellsize
!  start: starting cell, numcells=number of cells to order
!         (useful for ordering subtrees)
!  ind: indexing in order of cellsize

subroutine cell_sort2(start,numcells,ind)
  include 'globals.h'
  integer,intent(out) :: ind(*)
  integer,intent(in)  :: start,numcells
  integer :: facell, lacell, nacell, k, i
 
  ind(numcells) = start
  facell = 1;lacell = 1;nacell = 1
10  if(lacell.LT.numcells) then  
      do k = 1, nsubcell 
        do i = facell,lacell 
          if (subp(ind(numcells-i+1),k).GE.nbods1) then 
            nacell = nacell + 1 
            ind(numcells-nacell+1) = subp(ind(numcells-i+1),k) 
          endif 
        end do 
      enddo
      facell = lacell + 1 
      lacell = nacell 
      go to 10 
    endif
  if (nacell .ne. numcells) then
!$omp critical
    print*, "cellsort2",nacell,numcells
    call terror ('  cell_sort2: INCONSISTENT CELL COUNT')
!$omp end critical
  endif 
end subroutine cell_sort2

! subcellindex: function to determine value of subcell index 

integer function subcellindex2(pcell,ppar)
  real,intent(in),dimension(3) :: pcell,ppar 
  integer, parameter,dimension(3) :: nindex=(/ 4,2,1 /)
  integer :: k 
 
  subcellindex2 = 1
  do k = 1, 3 
    if (ppar(k).GE.pcell(k)) subcellindex2 = subcellindex2 + nindex(k) 
  enddo
end function subcellindex2


! loadtree5: loadtree with explicit ordering of particles 
subroutine loadtree5(option,startcell)
  include 'globals.h'
  character (len=4), intent(in) :: option
  integer,intent(in) :: startcell
  integer :: ilower,iupper,i,j,success 
  integer :: lincells,lcell,hcell,llcell,lhcell,low,high
  integer :: mythread,totalthread,maxthread,shift,maxcell,lnumbodies,restbodies
  integer :: checkredu,checkbody,checkcell,totcell,totredu,totbody,targetn 
  integer :: omp_get_max_threads,omp_get_thread_num ,omp_get_num_threads
  real, parameter :: cellfactor=float(ncells)/nbodsmax
  real :: time1,time2,mintime,maxtime,tottime,utime1,utime2
   
  select case (option)
  case('all ')
    ilower = 1 
    iupper = nbodies 
  case('sph ')
    ilower = 1 
    iupper = nsph
  case('star')
    ilower = nbodies - nstar + 1 
    iupper = nbodies
  case('bh  ')
    ilower = nbodies - nbh +1
    iupper = nbodies 
  case default
    ilower = nsph + 1 
    iupper = nbodies
  end select
  
  incells=1
  root=nbods1+startcell
  subp(root,1:nsubcell)=0
  cellsize(root) = rsize
  pshape(root)=1
  do j=1,ndim
    pos(root,j) = rmin(j) + 0.5*rsize 
    bottom(root,j) = rmin(j) 
  enddo
 
  cellstart(root)=1
  do i=ilower,iupper
    order_bodlist(i-ilower+1)=i
  enddo
  restbodies=iupper-ilower+1
  npercell(root)=restbodies

  if(restbodies.eq.0) return
  if(restbodies.eq.1) then
    subp(root,1)=order_bodlist(1)
    call tree_reduction(root,1,option)
    return
  endif 

  hcell=root
  lcell=root-1 
  maxthread=1
!$ maxthread=omp_get_max_threads()
  lnumbodies=restbodies/maxthread
  if(lnumbodies.LT.100) lnumbodies=1

  call parttree(lcell,hcell,startcell+incells,lincells, &
                  lnumbodies,nbodcell,success) 
  incells=incells+lincells
  if(success.ne.0) then
    print*,success
    call terror(' loadtree error 1')
  endif
  restbodies=sum(npercell(lcell+1:hcell))
  totcell=0
  totredu=0
  totbody=0
 
  mintime=1.e10
  maxtime=0.
  tottime=0

!$omp parallel private(mythread,totalthread,shift,lincells,i,checkcell,checkbody,checkredu &
!$omp ,llcell,lhcell,maxcell,low,high,success,lnumbodies,targetn,time1,time2) &
!$omp shared(startcell,incells,lcell,hcell,option,restbodies) &
!$omp reduction( MIN : mintime) & 
!$omp reduction( MAX : maxtime) &
!$omp reduction(+ : totredu,totcell,totbody,tottime) if(hcell-lcell.gt.maxthread)
  call cpu_time(time1)
  mythread=0
  totalthread=1
!$ mythread=omp_get_thread_num()    
!$ totalthread=omp_get_num_threads()

  maxcell=0
  lhcell=lcell
  i=0
  do while(i.LE.mythread)
    llcell=lhcell
    shift=maxcell
    lnumbodies=restbodies-maxcell
    targetn=lnumbodies/(totalthread-i)
    do while(maxcell-shift.LT.targetn.AND. &
              (maxcell-shift+npercell(lhcell+1)-targetn)*(totalthread-i-1).LT.(targetn-maxcell+shift))
      if(lhcell.GE.hcell) call terror(' loadtree distr')
      lhcell=lhcell+1
      maxcell=maxcell+npercell(lhcell) 
    enddo
    i=i+1
  enddo 

  lnumbodies=maxcell-shift
  shift=shift*cellfactor
  maxcell=maxcell*cellfactor

  shift=shift+incells+startcell
  maxcell=MIN(nbodcell,maxcell+nbodsmax+incells+startcell)
 
  checkbody=0
  checkcell=0
  checkredu=0
  do i=llcell+1,lhcell
    low=i-1;high=i
    call parttree(low,high,shift+checkcell,lincells,1, maxcell,success)
    if(success.ne.0.OR.high.NE.low) then
!$omp critical
      print*,shift,maxcell,lnumbodies
      print*,incells,restbodies
      print*,option,success,mythread,totalthread
      call terror(' loadtree error 2')
!$omp end critical    
    endif   
    call tree_reduction(i,lincells+1,option)
    checkbody=checkbody+npercell(i)
    checkredu=checkredu+lincells+1
    checkcell=checkcell+lincells
  enddo
  totcell=totcell+checkcell
  totbody=totbody+checkbody
  totredu=totredu+checkredu
  call cpu_time(time2)
  mintime=MIN(mintime,time2-time1)
  maxtime=MAX(maxtime,time2-time1)
  tottime=tottime+time2-time1
!$omp end parallel
  totredu=totredu+lcell+1-nbods1-startcell
  totcell=totcell+incells
  if(restbodies.NE.totbody.OR.totredu.NE.totcell) then
    print*,restbodies,totredu,totcell,totbody
    call terror('loadtree mismatch')
  endif
  call tree_reduction(root,lcell+1-nbods1-startcell,option)
  incells=totcell
  if(verbosity.GT.0) then
    write(*,'(" <loadtree> time:", 3f8.2)') mintime,maxtime,tottime
    print*,'<loadtree> ',option,incells
  endif
end subroutine loadtree5

subroutine parttree(lcell,hcell,shift,nnewcell,mxincell,maxnewcell,success)
  include 'globals.h'
  integer, intent(inout) :: lcell,hcell,nnewcell,success
  integer, intent(in)    :: mxincell,shift,maxnewcell
  integer :: i,j, k,cellp(0:nsubcell)
  integer :: subcellindex2
  integer :: maxnsub,subind,particle,nsub,newcell,firstnew
  real :: ppos(3),cpos(3),pm1(nsubcell,ndim)
  data pm1/ 4* - 1., 4*1., 2* - 1., 2*1., 2* - 1., 2*1., -1., 1., &
                                              -1.,1.,-1., 1., -1., 1./
  success=0
  nnewcell=0
  if(hcell.LE.lcell) return
  maxnsub=0
  do i=lcell+1,hcell
    if(maxnsub.LT.npercell(i)) maxnsub=npercell(i)
  enddo
  newcell=nbodsmax+shift
  firstnew=nbodsmax+shift+1
  do while(maxnsub.GT.mxincell)
    maxnsub=0
    do i=lcell+1,hcell
      if(npercell(i).LT.2) then
        success=1
        return
      endif
! end tree if nbodies <= 8
      if(npercell(i).LE.8) then   
        subp(i,1:npercell(i))= &
          order_bodlist(cellstart(i):cellstart(i)+npercell(i)-1)
        cycle
      endif   
      do j=cellstart(i),cellstart(i)+npercell(i)-1
        particle=order_bodlist(j)
        cpos=pos(i,1:3)
        ppos=pos(particle,1:3)
        subind=subcellindex2(cpos,ppos)    
        subind=peanovalue(subind,pshape(i))
        srlist(j-cellstart(i)+1)=subind
        subp(i,subind)=subp(i,subind)+1
      enddo   
      cellp(0)=0
      nsub=0 
      do j=1,nsubcell
        cellp(j)=cellp(j-1)+nsub
        nsub=subp(i,j)
        if(nsub.GT.maxnsub) maxnsub=nsub
        if(nsub.GT.1) then
          nnewcell=nnewcell+1
          newcell=nbodsmax+nnewcell+shift
          if(newcell.GT.maxnewcell) then
            success=2
            return
          endif
          npercell(newcell)=nsub
          cellstart(newcell)=cellp(j)+cellstart(i)
          cellsize(newcell)=cellsize(i)/2
          do k=1,ndim
            pos(newcell,k)=pos(i,k)+pm1(indexfromp(j,pshape(i)),k)*cellsize(newcell)/2
!            pos(newcell,k)=pos(i,k)+pm1(j,k)*cellsize(newcell)/2
            bottom(newcell,k)=pos(newcell,k)-cellsize(newcell)/2
          enddo
          subp(newcell,1:nsubcell)=0
          subp(i,j)=newcell             
          pshape(newcell)=pshapefromp(j,pshape(i))
        endif
      enddo
             
      do j=cellstart(i),cellstart(i)+npercell(i)-1
        particle=order_bodlist(j)
        subind=srlist(j-cellstart(i)+1)
        bodlist(1+cellp(subind))=particle
        cellp(subind)=cellp(subind)+1
        if(subp(i,subind).EQ.1) then
          subp(i,subind)=particle
        endif
      enddo 
   
      order_bodlist(cellstart(i):cellstart(i)+npercell(i)-1)= &
        bodlist(1:npercell(i))   
    enddo
    lcell=firstnew-1
    hcell=newcell
    firstnew=hcell+1
  enddo
end subroutine parttree

subroutine tree_reduction(start,numcells,option)
  integer, intent(in) :: start,numcells
  character (len=4), intent(in) :: option

  if(numcells.LE.0) return
 
  select case (option)
  case('all ')
    call tree_reduction_grav(start,numcells)
  case('sph ')
    call tree_reduction_sph(start,numcells)
  case('star')
    call tree_reduction_fuv(start,numcells)
  case('bh  ')
    return
  case default
    call tree_reduction_grav(start,numcells)
  end select

end subroutine tree_reduction

subroutine tree_reduction_grav(start,numcells)
  include 'globals.h'
  integer,intent(in) :: start,numcells
  integer p,q,i,j,k,l,m,n,mupper

  call cell_sort2(start,numcells,bodlist)

  do i=1,numcells
    p=bodlist(i)
    mass(p)=0.
    starfuv(p)=0.
    pos(p,1)=0.
    pos(p,2)=0.
    pos(p,3)=0.
    epsgrav(p)=0.
  enddo

  if(usequad) then
    do k=1,2*ndim-1
      do i=1, numcells
        p=bodlist(i)
        quad(p,k)=0.
      enddo
    enddo    
  endif

  do i=1,numcells
    p=bodlist(i)
    do j=1,nsubcell 
      q=subp(p,j)
      if(q.GT.0) then
        mass(p)=mass(p)+mass(q)
        pos(p,1:3)=pos(p,1:3)+mass(q)*pos(q,1:3)
        starfuv(p)=starfuv(p)+starfuv(q)
        epsgrav(p)=epsgrav(p)+mass(q)*epsgrav(q)
      endif
    enddo
    if(mass(p).lt.0) call terror('tree_reduction_grav: mass<0')
    if(mass(p).ne.0) then
      pos(p,1:3)=pos(p,1:3)/mass(p)
      epsgrav(p)=epsgrav(p)/mass(p) 
    else
      pos(p,1:3)=bottom(p,1:3)+cellsize(p)/2
      epsgrav(p)=eps
    endif
  enddo
   
  if(usequad) then
    if(ndim.GT.2) then
      mupper=2
    else
      mupper=ndim
    endif


    do i=1,numcells
      p=bodlist(i)
      do j=1,nsubcell
        q=subp(p,j)
        if(q.GT.0) then
          do m=1,mupper
            do n=m,ndim
              l=(m-1)*(ndim-1)+n
              quad(p,l)=quad(p,l)+ & 
                mass(q)*(3.*(pos(q,m)-pos(p,m))*(pos(q,n)-pos(p,n)))
              if(m.EQ.n) then
                do k=1,ndim
                  quad(p,l)=quad(p,l)-mass(q)*(pos(q,k)-pos(p,k))**2
                enddo
              endif
              if(q.ge.nbods1) quad(p,l)=quad(p,l)+quad(q,l)
            enddo
          enddo
        endif 
      enddo
    enddo 
  endif
end subroutine tree_reduction_grav

subroutine tree_reduction_fuv(start,numcells)
  include 'globals.h'
  integer,intent(in) :: start,numcells
  integer p,q,i,j,k,l,m,n,mupper

  call cell_sort2(start,numcells,bodlist)

  do i=1,numcells
    p=bodlist(i)
    mass(p)=1.  ! not zero..
    starfuv(p)=0.
    pos(p,1)=0.
    pos(p,2)=0.
    pos(p,3)=0.
    epsgrav(p)=0.
  enddo

  do i=1,numcells
    p=bodlist(i)
    do j=1,nsubcell 
      q=subp(p,j)
      if(q.GT.0) then
        pos(p,1:3)=pos(p,1:3)+starfuv(q)*pos(q,1:3)
        starfuv(p)=starfuv(p)+starfuv(q)
        epsgrav(p)=epsgrav(p)+starfuv(q)*epsgrav(q)
      endif
    enddo
    if(starfuv(p).NE.0) then
      pos(p,1:3)=pos(p,1:3)/starfuv(p)
      epsgrav(p)=epsgrav(p)/starfuv(p) 
    else
      pos(p,1:3)=bottom(p,1:3)+cellsize(p)/2
      epsgrav(p)=eps
    endif
  enddo 
end subroutine tree_reduction_fuv

subroutine tree_reduction_sph(start,numcells)
  include 'globals.h'
  integer,intent(in) :: start,numcells
  integer p,q,i,j,k,l,m,n,mupper

  call cell_sort2(start,numcells,bodlist)

  do i=1,numcells
    p=bodlist(i)
    hsmooth(p)=0.
  enddo
 
! note: instead of max of hsmooth it is probably better to
! determine the max range (and convert that to 'hsmooth'
! so:
!  hsmooth=MAX(hsmooth,bottom(p)-(bottom(q)-hsmooth(q)))
!  hsmooth=MAX(hsmooth,bottom(q)+cellsize(q)+hsmooth(q)-
!   (bottom(p)+cellsize(p)))
! and for particles:
!  hsmooth=MAX(hsmooth, bottom(p)-(pos(q)-hsmooth(q)))
!  hsmooth=MAX(hsmooth, pos(q)+hsmooth(q)-(bottom(p)+cellsize(p))) 
  do i=1,numcells
    p=bodlist(i)
    do j=1,nsubcell 
      q=subp(p,j)
!   if(q.GT.0) then
!    hsmooth(p)=MAX(hsmooth(p),hsmooth(q))
!   endif
      if(q.LT.nbods1) then
        if(q.GT.0) then
          do k=1,3
            hsmooth(p)=MAX(pos(q,k)+2*hsmooth(q)-bottom(p,k)-cellsize(p),2*hsmooth(p))/2
            hsmooth(p)=MAX(bottom(p,k)-pos(q,k)+2*hsmooth(q),2*hsmooth(p))/2
          enddo
        endif
      else
        do k=1,3
          hsmooth(p)=MAX(bottom(q,k)+cellsize(q)+2*hsmooth(q)-bottom(p,k)-cellsize(p),2*hsmooth(p))/2
          hsmooth(p)=MAX(bottom(p,k)-bottom(q,k)+2*hsmooth(q),2*hsmooth(p))/2
        enddo
      endif
    enddo
  enddo 
end subroutine tree_reduction_sph
