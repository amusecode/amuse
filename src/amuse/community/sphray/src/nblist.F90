!> \file nblist.F90

!> \brief The neighbor list module.
!! 
!<
module nblist_mod
use myf90_mod
use particle_system_mod, only: particle_system_type
use particle_system_mod, only: transformation_type
use sphpar_mod, only: sphpar_type, par2sphpar
use oct_tree_mod, only: oct_tree_type
use global_mod
implicit none

 private
 public :: nblist_type
! public :: search,fullsearch,prepare_neighbor_search,getnbscale,resetnblist
 public :: fullsearch,prepare_neighbor_search,getnbscale,resetnblist

 integer,parameter :: MAX_NBLIST_LENGTH=10000
 integer,parameter :: ndim=3

 type nblist_type
  integer :: nnb
  integer :: maxnnb
  integer(i8b) :: searchcell
  logical :: reuseable
  integer :: searchimage
  integer :: nsearchimages
  real :: searchrange
  real :: searchpos(ndim)
  real :: DfltMass
  real :: DfltTemp
  type(transformation_type) :: searchtrafos(3**ndim)
  type(sphpar_type), allocatable :: sphpar(:)
  type(oct_tree_type) :: stree 
 end type nblist_type
 
 contains

 subroutine prepare_neighbor_search(psys,nblist,DfltMass,DfltTemp)
 use oct_tree_mod, only: buildtree
 use oct_tree_mod, only: setparticleorder
 type(particle_system_type) psys
 type(nblist_type) nblist
 real, intent(in) :: DfltMass
 real, intent(in) :: DfltTemp
 real :: MB
  call makenblist(nblist,MAX_NBLIST_LENGTH)
  call buildtree(psys,nblist%stree,MB)
  call setparticleorder(psys,nblist%stree)  
  call setsearchimages(psys,nblist)
  nblist%DfltMass = DfltMass
  nblist%DfltTemp = DfltTemp
 end subroutine prepare_neighbor_search

 subroutine fullsearch(psys,nblist)
 type(particle_system_type) :: psys
 type(nblist_type) nblist
 if(nblist%searchimage.EQ.0) call nbError('nblist uninit.')
  do while(nblist%searchimage.LE.nblist%nsearchimages)
   call dosearchtree(psys,nblist)
   if(nblist%searchcell.NE.0) return 
   nblist%searchimage=nblist%searchimage+1
   nblist%searchcell=1
  enddo 
  nblist%searchcell=0
 end subroutine fullsearch

 subroutine dosearchtree(psys,nblist)
 use oct_tree_mod, only: dist2cell

 type(particle_system_type) :: psys
 type(nblist_type) :: nblist
 integer :: i, si
 integer(i8b) :: ip, this, next, daughter
 real :: dist2
 real :: dx(ndim), searchpos(ndim)

! integer :: k

  si = nblist%searchimage
  searchpos = nblist%searchpos * nblist%searchtrafos(si)%fac + & 
              nblist%searchtrafos(si)%shift
  next=nblist%searchcell
  do while(next.ne.0)
  this=next
  daughter=nblist%stree%cell(this)%daughter
  next=nblist%stree%cell(this)%next
  if(daughter.EQ.0) then
   if(nblist%nnb+nblist%stree%cell(next)%start- &
      nblist%stree%cell(this)%start .GT. nblist%maxnnb) then 
    nblist%reuseable=.FALSE.
    nblist%searchcell=this
    return
   endif
   do i=nblist%stree%cell(this)%start,nblist%stree%cell(next)%start-1
    ip=nblist%stree%partorder(i)
    dx=searchpos-psys%par(ip)%pos
    dist2=sum(dx**2)
    if(dist2.LE.nblist%searchrange**2) then
      nblist%nnb=nblist%nnb+1
      call par2sphpar(nblist%sphpar(nblist%nnb), &
               psys%par(ip), nblist%DfltMass, nblist%DfltTemp, nblist%searchtrafos(si))
    endif
   enddo
  else
   dist2=dist2cell(searchpos,nblist%stree%cell(this))
   if(dist2.LE.nblist%searchrange**2) next=daughter  
  endif
  enddo
  nblist%searchcell=next
 end subroutine dosearchtree

 function getnbscale(pos,nblist)
 use oct_tree_mod, only: getlocalscale
 real :: getnbscale
 type(nblist_type) :: nblist
 real :: pos(ndim)
  getnbscale=getlocalscale(pos,nblist%stree)
 end function getnbscale

 subroutine makenblist(nblist,maxnnb)
 type(nblist_type) :: nblist
 integer :: maxnnb
  nblist%nnb=0
  nblist%maxnnb=maxnnb
  nblist%searchpos=0
  nblist%searchrange=0
  nblist%searchcell=1
  nblist%reuseable=.FALSE.
  nblist%searchimage=1
  nblist%searchtrafos(1)%fac=1
  nblist%searchtrafos(1)%shift=0  
  allocate(nblist%sphpar(maxnnb))
 end subroutine makenblist

 subroutine resetnblist(nblist,range,pos)
   type(nblist_type) :: nblist
   real, optional :: range,pos(ndim)

!   integer ::maxnnb

     nblist%nnb=0
     nblist%searchcell=1
     nblist%searchimage=1
     nblist%reuseable=.FALSE.
     if(present(pos)) then
        nblist%searchpos=pos
     else 
     endif
     if(present(range)) then
        nblist%searchrange=range
     else 
     endif
 end subroutine resetnblist

 subroutine killnblist(nblist,maxnnb)
 type(nblist_type) :: nblist
 integer :: maxnnb
  nblist%nnb=0
  nblist%maxnnb=maxnnb
  nblist%searchpos=0
  nblist%searchrange=0
  nblist%searchcell=0
  deallocate(nblist%sphpar)
 end subroutine killnblist

 subroutine setsearchimages(psys,nblist)
 type(nblist_type) nblist
 type(particle_system_type) psys
 integer i,j,k(ndim)
 real top(ndim),bot(ndim)
 integer bbound(ndim),tbound(ndim),l,nsearchimages
 real lx(ndim),hx(ndim)
 logical boxexists

  top=psys%box%top
  bot=psys%box%bot
  bbound=psys%box%bbound
  tbound=psys%box%tbound
 
  nsearchimages=0

  do i=0,3**ndim-1
  l=i
  do j=1,ndim
   k(j)=mod(l,3)-1 
   l=l/3
  enddo 
  boxexists=.TRUE.
  do j=1,ndim
   lx(j)=bot(j)+k(j)*(top(j)-bot(j))
   hx(j)=top(j)+k(j)*(top(j)-bot(j))
   if(k(j).eq.-1.AND.bbound(j).EQ.0) boxexists=.FALSE.  
   if(k(j).eq.1.AND.tbound(j).EQ.0) boxexists=.FALSE.  
  enddo
  if(boxexists) then
   nsearchimages=nsearchimages+1
   do j=1,ndim    
    nblist%searchtrafos(nsearchimages)%fac(j)=1
    nblist%searchtrafos(nsearchimages)%shift(j)=0.
    if(k(j).eq.-1) then
    nblist%searchtrafos(nsearchimages)%fac(j)=bbound(j)
    nblist%searchtrafos(nsearchimages)%shift(j)=(-3*nblist%searchtrafos(nsearchimages)%fac(j)+1)*bot(j)/2+&
      (1+nblist%searchtrafos(nsearchimages)%fac(j))*top(j)/2
    endif  
    if(k(j).eq.1) then
    nblist%searchtrafos(nsearchimages)%fac(j)=tbound(j)
    nblist%searchtrafos(nsearchimages)%shift(j)=(-3*nblist%searchtrafos(nsearchimages)%fac(j)+1)*top(j)/2+&
      (1+nblist%searchtrafos(nsearchimages)%fac(j))*bot(j)/2
    endif  
   enddo
  endif
  enddo
  nblist%nsearchimages=nsearchimages
 end subroutine setsearchimages

 subroutine nbError(string,i)
  character*(*) :: string
  integer, optional :: i
  
  print*,'error detected'

  if(present(i)) then
   print*,string,i
  else
   print*,string
  endif
  stop
 end subroutine nbError

end module nblist_mod

!   dist2=0
!   do k=1,ndim
!    if(nblist%searchpos(k).LT.tree%cell(this)%bot(k)) then
!      dist2=dist2+(nblist%searchpos(k)-tree%cell(this)%bot(k))**2
!    else
!      if(nblist%searchpos(k).GT.tree%cell(this)%top(k)) then
!       dist2=dist2+(tree%cell(this)%top(k)-nblist%searchpos(k))**2
!      endif
!    endif
!   enddo



