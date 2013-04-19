module rndmod

 integer, parameter :: nrndtable=2048
 real :: rndtable(nrndtable)

contains

! simple code for random number sequences
! from numerical recipes..
function ranfb(ix)
  real :: y,ranfb
  integer :: ix,ix1,ix2,ia1,ia2,ic,j1,j2,j3

  ia1=273
  ia2=25301
  ic=226908345
  
  if (ix.lt.0) ix=-ix
  if (ix.gt.2**30) ix=MOD(ix,2**30)

  ix1=ix/(2**15)
  ix2=ix-ix1*(2**15)
  j1=MOD(ix1*ia2,2**15)*2**15
  j2=MOD(ix2*ia1,2**15)*2**15
  j3=ix2*ia2
  ix=MOD(j1+j2,2**30)
  ix=MOD(ix+j3,2**30)
  ix=MOD(ix+ic,2**30)
  y=FLOAT(ix)

  ranfb=y/2**30

end function

end module

subroutine setRND(rnseed)
  use rndmod
  integer i,rnseed
  do i=1,nrndtable
    rndtable(i)=ranfb(rnseed)
  enddo
end subroutine

function pRND(i)
  use rndmod
  real :: pRND
  integer :: i
  pRND=rndtable(mod(i,nrndtable)+1)
end function

module map_helpers  
  integer :: randomseed=45678910
  
  integer :: nmax=0

  real,allocatable :: gpos(:,:),radius(:),mass(:),opac(:)
  integer,allocatable :: pindex(:)
  
  integer :: nbod=0,firstfree=1
  
  logical :: searcheable=.FALSE.
  
  integer :: imsize(2)=(/640,480/)
  real :: width=1.
  real :: angle=45.
  real :: focus(3) = (/ 0.,0.,0. /)
  real :: viewpoint(3) = (/ 0.,1.,0. /)
  real :: direction(3)=(/ 0.,-1.,0. /)
  real :: upvector(3)=(/ 0.,0.,1. /)
  integer :: ext=0
  real :: zm=0.001
  character(len=11) :: projection_mode="parallel"
  
  
contains

subroutine reallocate_arrays(n)
  integer :: n,test
  real,allocatable :: gpos_(:,:),radius_(:),mass_(:),opac_(:)
  integer,allocatable :: pindex_(:)

  if(nbod.GT.0) then
    if(.not.allocated(gpos)) call terror("mem alloc impossible error 1")
    allocate(gpos_(nbod,3),radius_(nbod),mass_(nbod),opac_(nbod),pindex_(nbod),stat=test)
    if(test.NE.0) call terror('mem alloc error')
    gpos_(1:nbod,1:3)=gpos(1:nbod,1:3)
    radius_(1:nbod)=radius(1:nbod)
    mass_(1:nbod)=mass(1:nbod)
    opac_(1:nbod)=opac(1:nbod)
    pindex_(1:nbod)=pindex(1:nbod)
  endif
    if(allocated(gpos)) deallocate(gpos,radius,mass,opac,pindex)
    allocate(gpos(n,3),radius(n),mass(n),opac(n),pindex(n),stat=test)
    if(test.NE.0) call terror('mem alloc error')
  if(nbod.GT.0) then
    gpos(1:nbod,1:3)=gpos_(1:nbod,1:3)
    radius(1:nbod)=radius_(1:nbod)
    mass(1:nbod)=mass_(1:nbod)
    opac(1:nbod)=opac_(1:nbod)
    pindex(1:nbod)=pindex_(1:nbod)
    deallocate(gpos_,radius_,mass_,opac_,pindex_)
  endif  
  nmax=n
end subroutine

function map_init() result(ret)
  integer :: ret
  ret=-1
  if(projection_mode.EQ."parallel") then
    call map_init_parallel_projection()
    ret=0
  else
    if(projection_mode.EQ."perspective") call map_init_perspective_projection()
    ret=0
  endif
end function

subroutine map_init_parallel_projection()
  use makemapmod
  call InitMap(imsize,width,focus,direction,upvector,0,ext,zm)  
end subroutine

subroutine map_init_perspective_projection()
  use makemapmod
  call InitMap(imsize,angle,focus,viewpoint,upvector,1,ext,zm)  
end subroutine

subroutine map_generate_projection()
  use makemapmod
  call setRND(randomseed)
  call clean_array()
  call project(nbod,gpos,radius,mass,opac,pindex)
end subroutine

subroutine map_reset()
  use makemapmod
  call EndMap()
end subroutine

subroutine map_erase()
  use makemapmod
  call EraseMap()
end subroutine

subroutine clean_array()
  integer :: n,i
  n=0
  i=0
  do while(n.LT.nbod)
    i=i+1
    if(mass(i).GE.0) then
      n=n+1
      mass(n)=mass(i)
      gpos(n,1:3)=gpos(i,1:3)
      radius(n)=radius(i)
      opac(n)=opac(i)
      pindex(n)=pindex(i)
    endif
  enddo
  firstfree=nbod+1
  searcheable=.FALSE.
end subroutine

function add_particle(id,m,x,y,z,r,o,i) result(ret)
  integer ret
  integer id,i
  real :: m,x,y,z,r,o
   
  if(firstfree.GT.nmax) then
    call clean_array()
    if(firstfree.GT.nmax) then
        call reallocate_arrays(max(nmax*2,1024))
    endif
  endif
   
  gpos(firstfree,1)=x  
  gpos(firstfree,2)=y  
  gpos(firstfree,3)=z  
  mass(firstfree)=m
  radius(firstfree)=r
  opac(firstfree)=o
  pindex(firstfree)=i
  if(pindex(firstfree).LT.0) pindex(firstfree)=firstfree
  id=pindex(firstfree)

  firstfree=firstfree+1
  nbod=nbod+1
   
  searcheable=.FALSE.
   
  ret=0
 end function

function set_particle_state(id,m,x,y,z,r,o) result(ret)
  integer ret
  integer id,index
  real :: m,x,y,z,r,o

  index=map_find_particle(id,nbod,pindex)
  if(index.LT.0) then
    ret=index
    return
  endif  

  gpos(index,1)=x  
  gpos(index,2)=y  
  gpos(index,3)=z  
  mass(index)=m
  radius(index)=r
  opac(index)=o
   
  ret=0
 end function

function set_particle_weight(id,m) result(ret)
  integer ret
  integer id,index
  real :: m

  index=map_find_particle(id,nbod,pindex)
  if(index.LT.0) then
    ret=index
    return
  endif  

  mass(index)=m
   
  ret=0
end function

function get_particle_state(id,m,x,y,z,r,o) result(ret)
  integer ret
  integer id,index
  real :: m,x,y,z,r,o

  index=map_find_particle(id,nbod,pindex)
  if(index.LT.0) then
    ret=index
    return
  endif  

  x=gpos(index,1)
  y=gpos(index,2)  
  z=gpos(index,3)  
  m=mass(index)
  r=radius(index)
  o=opac(index)
   
  ret=0
end function

function map_remove_particle(id) result(ret)
  integer :: id
  integer :: index,ret

  index=map_find_particle(id,nbod,pindex)
  if(index.LT.0) then
    ret=index
    return
  endif  

  mass(index)=-1
  nbod=nbod-1

  ret=0
end function

function map_find_particle(id,n,ids) result(p)
  use hashMod
  integer :: id,n,p,ids(*)
  integer, save :: nsearch
  
  p=-1
  
  if(.not.searcheable) then
    call clean_array()
    nsearch=n
    if(n.GT.0) call initHash(nsearch/2+1,nsearch,ids) 
    searcheable=.TRUE.
  endif  
   
  if(n.GT.0) p=find(id,ids)

end function

subroutine map_set_projection_mode(x)
  character(len=15), intent(in) :: x
  projection_mode=x
end subroutine
subroutine map_get_projection_mode(x)
  character(len=15), intent(out) :: x
  x=projection_mode
end subroutine


function map_set_random_seed(rs) result(ret)
  integer :: rs,ret
  randomseed=rs
  ret=0
end function 

function map_get_random_seed(rs) result(ret)
  integer :: rs,ret
  rs=randomseed
  ret=0
end function 

subroutine map_set_zm(x)
  real :: x
  zm=x
end subroutine
subroutine map_get_zm(x)
  real :: x
  x=zm
end subroutine

subroutine map_set_ext(x)
  integer :: x
  ext=x
end subroutine
subroutine map_get_ext(x)
  integer :: x
  x=ext
end subroutine

subroutine map_set_angle(x)
  real :: x
  angle=x
end subroutine
subroutine map_get_angle(x)
  real :: x
  x=angle
end subroutine

subroutine map_set_width(x)
  real :: x
  width=x
end subroutine
subroutine map_get_width(x)
  real :: x
  x=width
end subroutine

subroutine map_set_imsize(nx,ny)
  integer :: nx,ny
  imsize(1)=nx
  imsize(2)=ny
end subroutine
subroutine map_get_imsize(nx,ny)
  integer :: nx,ny
  nx=imsize(1)
  ny=imsize(2)
end subroutine

subroutine map_set_focus(x,y,z)
  real :: x,y,z
  focus(1)=x
  focus(2)=y
  focus(3)=z
end subroutine
subroutine map_get_focus(x,y,z)
  real :: x,y,z
  x=focus(1)
  y=focus(2)
  z=focus(3)
end subroutine

subroutine map_set_viewpoint(x,y,z)
  real :: x,y,z
  viewpoint(1)=x
  viewpoint(2)=y
  viewpoint(3)=z
end subroutine
subroutine map_get_viewpoint(x,y,z)
  real :: x,y,z
  x=viewpoint(1)
  y=viewpoint(2)
  z=viewpoint(3)
end subroutine

subroutine map_set_direction(x,y,z)
  real :: x,y,z
  direction(1)=x
  direction(2)=y
  direction(3)=z
end subroutine
subroutine map_get_direction(x,y,z)
  real :: x,y,z
  x=direction(1)
  y=direction(2)
  z=direction(3)
end subroutine

subroutine map_set_upvector(x,y,z)
  real :: x,y,z
  upvector(1)=x
  upvector(2)=y
  upvector(3)=z
end subroutine
subroutine map_get_upvector(x,y,z)
  real :: x,y,z
  x=upvector(1)
  y=upvector(2)
  z=upvector(3)
end subroutine

function map_get_pic(i,j,pvalues,n) result(ret)
  use makemapMod
  integer :: ret,n,i(n),j(n),k
  real :: pvalues(n)
  
  if(any(i.LT.1).OR.any(i.GT.imsize(1)).OR.any(j.LT.1).OR.any(j.GT.imsize(2))) then
    ret=-1
    return
  endif  

  do k=1,n
    pvalues(k)=pic(i(k),j(k))
  enddo

  ret=0
  
end function

function map_get_opdepth(i,j,pvalues,n) result(ret)
  use makemapMod
  integer :: ret,n,i(n),j(n),k
  real :: pvalues(n)
    
  if(any(i.LT.1).OR.any(i.GT.imsize(1)).OR.any(j.LT.1).OR.any(j.GT.imsize(2))) then
    ret=-1
    return
  endif  

  do k=1,n
    pvalues(k)=opdepth(i(k),j(k))
  enddo

  ret=0
  
end function


end module
