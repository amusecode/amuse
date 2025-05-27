!> \file octtree3p.F90

!> \brief The peano hilbert oct tree module 
!!
!<

module QuadTreeMod
 use particleMod
 implicit none
 
 private
 public :: quadtree_type,cell_type
 public :: setparticleorder,buildtree,dist2cell,getlocalscale

 integer, parameter :: nsubcell=2**ndim

 integer, parameter :: partInCell=12

 type cell_type
  integer :: start
  integer :: next
  integer :: daughter
  integer :: pshape
  real,dimension(ndim) :: bot
  real,dimension(ndim) :: top
 end type

 type quadtree_type
  integer :: np=0,ncells=0,maxcells,nplastcell
  integer, allocatable :: order(:)
  integer, allocatable :: cellorder(:)
  type(cell_type), allocatable :: cell(:) 
 end type

        INTEGER peanovalue(8,24)
        DATA peanovalue/                                                 &
     &    1,    8,    4,    5,    2,    7,    3,    6,                   &
     &    1,    2,    4,    3,    8,    7,    5,    6,                   &
     &    7,    8,    6,    5,    2,    1,    3,    4,                   &
     &    3,    6,    4,    5,    2,    7,    1,    8,                   &
     &    1,    4,    8,    5,    2,    3,    7,    6,                   &
     &    5,    8,    4,    1,    6,    7,    3,    2,                   &
     &    1,    2,    8,    7,    4,    3,    5,    6,                   &
     &    3,    2,    4,    1,    6,    7,    5,    8,                   &
     &    7,    2,    6,    3,    8,    1,    5,    4,                   &
     &    5,    6,    4,    3,    8,    7,    1,    2,                   &
     &    7,    8,    2,    1,    6,    5,    3,    4,                   &
     &    7,    6,    8,    5,    2,    3,    1,    4,                   &
     &    3,    4,    6,    5,    2,    1,    7,    8,                   &
     &    1,    4,    2,    3,    8,    5,    7,    6,                   &
     &    5,    4,    8,    1,    6,    3,    7,    2,                   &
     &    5,    8,    6,    7,    4,    1,    3,    2,                   &
     &    1,    8,    2,    7,    4,    5,    3,    6,                   &
     &    7,    2,    8,    1,    6,    3,    5,    4,                   &
     &    5,    6,    8,    7,    4,    3,    1,    2,                   &
     &    3,    2,    6,    7,    4,    1,    5,    8,                   &
     &    7,    6,    2,    3,    8,    5,    1,    4,                   &
     &    5,    4,    6,    3,    8,    1,    7,    2,                   &
     &    3,    4,    2,    1,    6,    5,    7,    8,                   &
     &    3,    6,    2,    7,    4,    5,    1,    8 /

        INTEGER pshapefromi(8,24)
        DATA pshapefromi/                                                &
     &     2,    3,    4,    4,    5,    6,    5,    6,                  &
     &     1,    7,    8,    7,    9,   10,    8,   10,                  &
     &    11,    1,   11,   12,   13,    9,   13,   12,                  &
     &    12,    8,    1,    1,   12,    8,   10,   13,                  &
     &    14,   13,   12,   13,    1,    1,   15,   15,                  &
     &    10,   16,   10,    8,    1,    1,   15,   15,                  &
     &    17,    2,   18,   19,   20,    2,   20,   19,                  &
     &    18,   18,    2,    6,    4,    4,    2,   20,                  &
     &    21,   20,   21,   20,    2,    3,   18,   18,                  &
     &     6,    2,    6,   19,   22,    2,    4,   19,                  &
     &     3,   17,   23,   18,    3,   21,   23,   21,                  &
     &    18,   18,    5,    3,    4,    4,   21,    3,                  &
     &     3,    5,   23,    5,    3,   22,   23,    4,                  &
     &     5,   23,   17,   17,   21,   23,   22,   22,                  &
     &    22,   22,   19,   23,    5,    6,    5,    6,                  &
     &    19,    6,   17,   17,   19,   20,   22,   22,                  &
     &     7,   11,   14,   16,   24,   24,   14,   16,                  &
     &    12,    8,    7,   11,   12,    8,    9,    9,                  &
     &    16,    7,   15,    7,   16,   10,   24,   10,                  &
     &     9,    9,   24,   24,    7,   16,    7,    8,                  &
     &     9,    9,   24,   24,   14,   11,   12,   11,                  &
     &    15,   15,   14,   16,   10,   13,   14,   16,                  &
     &    11,   14,   11,   15,   13,   14,   13,   24,                  &
     &    21,   20,   21,   20,   17,   17,   19,   23 /

        INTEGER pshapefromp(8,24)
        DATA pshapefromp/                                                &
     &     2,    5,    5,    4,    4,    6,    6,    3,                  &
     &     1,    7,    7,    8,    8,   10,   10,    9,                  &
     &     9,   13,   13,   12,   12,   11,   11,    1,                  &
     &    10,   12,   12,    1,    1,    8,    8,   13,                  &
     &    14,    1,    1,   13,   13,   15,   15,   12,                  &
     &     8,   15,   15,   10,   10,    1,    1,   16,                  &
     &    17,    2,    2,   20,   20,   19,   19,   18,                  &
     &     6,   18,   18,    2,    2,    4,    4,   20,                  &
     &     3,   20,   20,   18,   18,   21,   21,    2,                  &
     &     4,   19,   19,    6,    6,    2,    2,   22,                  &
     &    18,   23,   23,   21,   21,    3,    3,   17,                  &
     &    21,    4,    4,    3,    3,   18,   18,    5,                  &
     &    22,    3,    3,    5,    5,   23,   23,    4,                  &
     &     5,   17,   17,   23,   23,   22,   22,   21,                  &
     &    23,    6,    6,   22,   22,    5,    5,   19,                  &
     &    20,   22,   22,   19,   19,   17,   17,    6,                  &
     &     7,   14,   14,   24,   24,   16,   16,   11,                  &
     &    11,    8,    8,    9,    9,   12,   12,    7,                  &
     &    24,   10,   10,   16,   16,    7,    7,   15,                  &
     &    16,    9,    9,    7,    7,   24,   24,    8,                  &
     &    12,   24,   24,   11,   11,    9,    9,   14,                  &
     &    13,   16,   16,   15,   15,   14,   14,   10,                  &
     &    15,   11,   11,   14,   14,   13,   13,   24,                  &
     &    19,   21,   21,   17,   17,   20,   20,   23 /

         INTEGER indexfromp(8,24)
         DATA indexfromp/                                                &
     &    1,    5,    7,    3,    4,    8,    6,    2,                   &
     &    1,    2,    4,    3,    7,    8,    6,    5,                   &
     &    6,    5,    7,    8,    4,    3,    1,    2,                   &
     &    7,    5,    1,    3,    4,    2,    6,    8,                   &
     &    1,    5,    6,    2,    4,    8,    7,    3,                   &
     &    4,    8,    7,    3,    1,    5,    6,    2,                   &
     &    1,    2,    6,    5,    7,    8,    4,    3,                   &
     &    4,    2,    1,    3,    7,    5,    6,    8,                   &
     &    6,    2,    4,    8,    7,    3,    1,    5,                   &
     &    7,    8,    4,    3,    1,    2,    6,    5,                   &
     &    4,    3,    7,    8,    6,    5,    1,    2,                   &
     &    7,    5,    6,    8,    4,    2,    1,    3,                   &
     &    6,    5,    1,    2,    4,    3,    7,    8,                   &
     &    1,    3,    4,    2,    6,    8,    7,    5,                   &
     &    4,    8,    6,    2,    1,    5,    7,    3,                   &
     &    6,    8,    7,    5,    1,    3,    4,    2,                   &
     &    1,    3,    7,    5,    6,    8,    4,    2,                   &
     &    4,    2,    6,    8,    7,    5,    1,    3,                   &
     &    7,    8,    6,    5,    1,    2,    4,    3,                   &
     &    6,    2,    1,    5,    7,    3,    4,    8,                   &
     &    7,    3,    4,    8,    6,    2,    1,    5,                   &
     &    6,    8,    4,    2,    1,    3,    7,    5,                   &
     &    4,    3,    1,    2,    6,    5,    7,    8,                   &
     &    7,    3,    1,    5,    6,    2,    4,    8 / 





 contains

 subroutine setparticleorder(psys,tree)
 type(quadtree_type) :: tree
 type(particlesys_type) :: psys
 type(particle_type) :: part
 integer i,goal
  call check(psys,tree)
  call orderpsys(psys,tree%order) 
 end subroutine

 subroutine check(psys,tree)
 type(quadtree_type) :: tree
 type(particlesys_type) :: psys

  if(tree%np.NE.psys%npart) call error('tree, psys')
  
 end subroutine

 subroutine maketree(psys,tree)
 type(quadtree_type) :: tree
  type(particlesys_type) :: psys
 integer :: npart,ncells,err
  if(allocated(tree%order)) call killtree(tree)
  tree%maxcells=psys%npart/2+2*psys%npart/partInCell
  tree%np=psys%npart
  allocate(tree%order(1:tree%np), tree%cellorder(1:max(tree%np,tree%maxcells)),tree%cell(0:tree%maxcells),stat=err)
 if(err.ne.0) call error('allocation      ',err)
 end subroutine

 subroutine killtree(tree)
 type(quadtree_type) :: tree
 integer :: err
  deallocate(tree%order, tree%cellorder,tree%cell,stat=err)
  if(err.ne.0) call error('deallocation     ',err)
  tree%np=0
 end subroutine

 integer function subcellindex2(pcell,ppar)
 real, dimension(ndim) :: pcell,ppar
 integer, parameter,dimension(3) :: nindex=(/ 1,2,4 /)
 integer :: k 

 subcellindex2 = 1
 do k = 1, ndim
  if (ppar(k).GT.pcell(k)) subcellindex2 = subcellindex2 + nindex(k)
 enddo
 end function subcellindex2

 subroutine buildtree(psys,tree)
  type(quadtree_type) :: tree
  type(particlesys_type) :: psys

  call maketree(psys,tree)

  call inittree(psys,tree)  

  call parttree(psys,tree)

  call patchtree(tree)

  call makecellorder(tree)

  call shrinkcells2(psys,tree)

!  call adjustbox(psys,tree%cell(1)%bot,tree%cell(1)%top)
  
 end subroutine
  
 subroutine inittree(psys,tree)
  integer :: npart,i
  type(quadtree_type) :: tree
  type(particlesys_type) :: psys
  real :: bot(ndim),top(ndim),rsize
 
 if(psys%npart.LE.0) call error('no particles   ')
 if(psys%npart.NE.tree%np) call error('psys != tree   ')

 tree%cell(0)%start=psys%npart+1
 tree%cell(1)%start=1
 tree%cell(1)%next=0
 tree%cell(1)%pshape=1
 tree%ncells=1
 tree%nplastcell=psys%npart
 
 do i=1,ndim
  bot(i)=MINVAL(psys%part(1:psys%npart)%pos(i))
  top(i)=MAXVAL(psys%part(1:psys%npart)%pos(i))
 enddo
 
 rsize=MAXVAL(top-bot)
 tree%cell(1)%bot(1:ndim)=0.5*(top+bot)-0.5*rsize
 tree%cell(1)%top(1:ndim)=tree%cell(1)%bot(1:ndim)+rsize

! tree%cell(1)%bot=0
! tree%cell(1)%top=1
 
 do i=1,psys%npart
  tree%order(i)=i
 enddo

 end subroutine

 recursive subroutine parttree(psys,tree)
 type(quadtree_type) :: tree
 type(particlesys_type) :: psys
 integer, parameter :: offset_factor(3,8)=(/ 0,0,0,1,0,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,1,1,1 /)  
 real :: ccpos(ndim),csize(ndim)
 integer :: i,j,subp(0:nsubcell+1),subind
 integer :: ccell,nsub,nppcell,newcell,sister,pshape
  ccell=tree%ncells
  nppcell=tree%nplastcell
  if(nppcell.GT.partInCell) then
  ccpos=(tree%cell(ccell)%top+tree%cell(ccell)%bot)/2
  csize=tree%cell(ccell)%top-tree%cell(ccell)%bot
  pshape=tree%cell(ccell)%pshape
  subp=0
  do j=tree%cell(ccell)%start,tree%cell(ccell)%start+nppcell-1
   subind=subcellindex2(ccpos,psys%part(tree%order(j))%pos)    
   subind=peanovalue(subind,pshape)
   subp(subind)=subp(subind)+1
  enddo
  subp(0)=0
  subp(nsubcell+1)=nppcell
  do j=nsubcell,1,-1
   subp(j)=subp(j+1)-subp(j) 
  enddo
  if(subp(1).NE.0) call error('subp count  ')
  do j=tree%cell(ccell)%start,tree%cell(ccell)%start+nppcell-1
   subind=subcellindex2(ccpos,psys%part(tree%order(j))%pos)
   subind=peanovalue(subind,pshape)
   subp(subind)=subp(subind)+1
   tree%cellorder(subp(subind))=tree%order(j)
  enddo 
!  tree%order(tree%cell(ccell)%start:tree%cell(ccell)%start+nppcell-1)= &
!   tree%cellorder(1:nppcell)   
  do j=tree%cell(ccell)%start,tree%cell(ccell)%start+nppcell-1
   tree%order(j)=tree%cellorder(j-tree%cell(ccell)%start+1)
  enddo 
  sister=0
  do j=1,nsubcell
   if(subp(j)-subp(j-1).GT.0) then 
    tree%ncells=tree%ncells+1
    newcell=tree%ncells
    if(newcell.GT.tree%maxcells) call error('cell overflow  ',newcell)
    tree%nplastcell=subp(j)-subp(j-1)
    if(tree%cell(ccell)%daughter.EQ.0)tree%cell(ccell)%daughter=newcell
    tree%cell(newcell)%start=subp(j-1)+tree%cell(ccell)%start
    tree%cell(newcell)%next=ccell
    tree%cell(newcell)%daughter=0
!    tree%cell(newcell)%bot=tree%cell(ccell)%bot+offset_factor(1:ndim,j)*csize/2
    tree%cell(newcell)%bot=tree%cell(ccell)%bot+offset_factor(1:ndim,indexfromp(j,pshape))*csize/2    
    tree%cell(newcell)%top=tree%cell(newcell)%bot+csize/2
    tree%cell(newcell)%pshape=pshapefromp(j,pshape)
    tree%cell(sister)%next=newcell
    sister=newcell
    call parttree(psys,tree)
    endif
  enddo    
  endif
 end subroutine

 subroutine patchtree(tree)
 type(quadtree_type) :: tree
 integer ::i
  tree%cell(0)%next=0
  do i=1,tree%ncells
   if(tree%cell(i)%next.LT.i) tree%cell(i)%next=tree%cell(tree%cell(i)%next)%next
  enddo
 end subroutine 

 subroutine makecellorder(tree)
 type(quadtree_type) :: tree
 integer ccount,dcell,i,up,low,pcell
  ccount=tree%ncells
  tree%cellorder(ccount)=1
  up=ccount+1
  low=ccount
  do while(up.GT.low)
   do i=up-1,low,-1
    pcell=tree%cellorder(i)
    dcell=tree%cell(pcell)%daughter
    if(dcell.NE.0) then
      ccount=ccount-1
      tree%cellorder(ccount)=dcell
      dcell=tree%cell(dcell)%next    
      do while(dcell.NE.tree%cell(pcell)%next)
       ccount=ccount-1
       tree%cellorder(ccount)=dcell
       dcell=tree%cell(dcell)%next     
      enddo
    endif
   enddo
   up=low
   low=ccount
  enddo
  if(low.NE.1) call error('ordering    ',low) 
 end subroutine
 
 integer function ppcell(tree,i)
 type(quadtree_type) :: tree  
 integer ::i
  ppcell=tree%cell(tree%cell(i)%next)%start-tree%cell(i)%start
 end function

 subroutine shrinkcells(psys,tree)
 type(quadtree_type) :: tree
 type(particlesys_type) :: psys
  integer :: i0,i,j,k,n,next
  real x
  
  do i0=1,tree%ncells
   i=tree%cellorder(i0)
   do k=1,ndim
    x=tree%cell(i)%top(k)
    tree%cell(i)%top(k)=tree%cell(i)%bot(k)
    tree%cell(i)%bot(k)=x
   enddo
   next=tree%cell(i)%daughter
   if(next.EQ.0) then
    do j=tree%cell(i)%start,tree%cell(tree%cell(i)%next)%start-1
      do k=1,ndim
       x=psys%part(tree%order(j))%pos(k)
       if(tree%cell(i)%bot(k).gt.x)tree%cell(i)%bot(k)=x
       if(tree%cell(i)%top(k).lt.x)tree%cell(i)%top(k)=x
      enddo
    enddo
   else
    do while(next.NE.tree%cell(i)%next)
      do k=1,ndim
       x=tree%cell(next)%bot(k)
       if(tree%cell(i)%bot(k).gt.x)tree%cell(i)%bot(k)=x
       x=tree%cell(next)%top(k)
       if(tree%cell(i)%top(k).lt.x)tree%cell(i)%top(k)=x
      enddo
      next=tree%cell(next)%next
    enddo 
   endif   
  enddo
! print*,tree%cell(1)%top
! print*,tree%cell(1)%bot
 end subroutine

 subroutine shrinkcells2(psys,tree)
 type(quadtree_type) :: tree
 type(particlesys_type) :: psys
  integer :: i0,i,j,k,n,next
  real x,hsm
  
  do i0=1,tree%ncells
   i=tree%cellorder(i0)
   do k=1,ndim
    x=tree%cell(i)%top(k)
    tree%cell(i)%top(k)=tree%cell(i)%bot(k)
    tree%cell(i)%bot(k)=x
   enddo
   next=tree%cell(i)%daughter
   if(next.EQ.0) then
    do j=tree%cell(i)%start,tree%cell(tree%cell(i)%next)%start-1
      hsm=psys%part(tree%order(j))%hsm
      do k=1,ndim
       x=psys%part(tree%order(j))%pos(k)
       if(tree%cell(i)%bot(k).gt.x-hsm)tree%cell(i)%bot(k)=x-hsm
       if(tree%cell(i)%top(k).lt.x+hsm)tree%cell(i)%top(k)=x+hsm
      enddo
    enddo
   else
    do while(next.NE.tree%cell(i)%next)
      do k=1,ndim
       x=tree%cell(next)%bot(k)
       if(tree%cell(i)%bot(k).gt.x)tree%cell(i)%bot(k)=x
       x=tree%cell(next)%top(k)
       if(tree%cell(i)%top(k).lt.x)tree%cell(i)%top(k)=x
      enddo
      next=tree%cell(next)%next
    enddo 
   endif   
  enddo
! print*,tree%cell(1)%top
! print*,tree%cell(1)%bot
 end subroutine

 function getlocalscale(pos,tree) result(scale)
 real :: pos(ndim),scale,dist2
 type(quadtree_type) :: tree
 integer :: next,old
  
 next=1
 old=1
 do while(tree%cell(next)%daughter.NE.0)
  old=next
  next=tree%cell(next)%daughter
  do while(dist2cell(pos,tree%cell(next)).GT.0)
   next=tree%cell(next)%next
  enddo
 enddo 
 if(next.eq.0) call error('getlocalscale: noc')
 scale=sum(tree%cell(old)%top-tree%cell(old)%bot)/ndim
 end function
 
 function dist2cell(pos,cell) result(dist2)
  type(cell_type) :: cell
  real :: pos(ndim),lpos(ndim),upos(ndim),dist2
  integer :: k
  dist2=0
  lpos=cell%bot
  upos=cell%top
  do k=1,ndim
   if(pos(k).LT.lpos(k)) then
    dist2=dist2+(lpos(k)-pos(k))**2
   else
    if(pos(k).GT.upos(k)) then
      dist2=dist2+(upos(k)-pos(k))**2
    endif
   endif
  enddo
 end function


 subroutine error(string,i)
  character(*) :: string
  integer, optional :: i
  
  print*,'error detected'

  if(present(i)) then
   print*,string,i
  else
   print*,string
  endif
  stop
end subroutine


end module
