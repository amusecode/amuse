!> \file octtree3.F90

!> \brief the oct tree module 
!!
!<

!> oct-tree construction and types
module oct_tree_mod
use myf03_mod
use particle_system_mod
implicit none
 
 integer, private, parameter :: nsubcell = 2**3  !< number of daughter cells 
 integer, private, parameter :: bytespercell = 3 * 4 + 4 * 12 !< bytes per cell
 integer, private, parameter :: PartInCell = 12 !< dflt min particles in a leaf
 real, private :: tree_storage_factor !< guess at storage space for tree

 private check, error

!> oct-tree cell type
 type cell_type
    integer(i4b) :: start  !< index of head particle in the cell 
    integer(i4b) :: next   !< index of next cell on current or higher lvl  
    integer(i4b) :: daughter !< index of first daughter cell
    real(r4b) :: bot(1:3) !< oct-tree lower boundaries (fit particle positions)
    real(r4b) :: top(1:3) !< oct-tree upper boundaries (fit particle positions)
    real(r4b) :: botrange(1:3) !< AABB lower boundaries (fit smoothing lengths)
    real(r4b) :: toprange(1:3) !< AABB upper boundaries (fit smoothing lengths)
 end type cell_type
 
 !> oct-tree type
 type oct_tree_type
    integer(i8b) :: ppc    !< particles per cell
    integer(i8b) :: np     !< number of particles (rezero if theres a problem)
    integer(i8b) :: ncells     !< number of cells rezero if theres a problem
    integer(i8b) :: maxcells   !< maximum number of cells allowed
    integer(i8b) :: nplastcell !< work variable used to construct tree 
    integer(i4b), allocatable    :: partorder(:) !< ordered list of the pars
    integer(i4b), allocatable    :: cellorder(:) !< ordered list of the cells
    type(cell_type), allocatable :: cell(:)      !< the array of tree cells
 end type oct_tree_type


 contains



!> puts the particles in the particle system into the 
!! order described by the tree
!-----------------------------------------------------------------------------
subroutine setparticleorder(psys,tree)
  type(oct_tree_type) :: tree     !< the oct-tree
  type(particle_system_type) :: psys !< the particle system
  call check(psys,tree)
  call particle_system_order_particles(psys,tree%partorder)
end subroutine setparticleorder


!> makes sure the tree and the particle system know about the same number of 
!> particles
!-----------------------------------------------------------------------------
subroutine check(psys,tree)
  type(oct_tree_type), intent(in) :: tree  !< the oct-tree
  type(particle_system_type), intent(in) :: psys !< the particle system
  character(clen) :: string
  
  string = "N tree particles not equal to N psys particles"
  if(tree%np.NE.size(psys%par)) call error(string)
  
end subroutine check


!> constructs an oct-tree from the particle system with the minimum number of 
!> particles per leaf as an optional input (ppc)
!-----------------------------------------------------------------------------
 subroutine buildtree(psys, tree, MBalloc, ppc)

   character(clen), parameter :: myname="buildtree"
   logical, parameter :: crash=.true.
   integer, parameter :: verb=2
   character(clen) :: str,fmt

   type(particle_system_type), intent(in) :: psys !< input particle system
   type(oct_tree_type), intent(inout) :: tree     !< oct-tree
   real(r8b), intent(out) :: MBalloc              !< MB of memory allocated for tree
   integer(i4b), intent(in), optional :: ppc      !< optional minimum particles per leaf
   integer(i4b) :: err

   if (present(ppc)) then
      tree%ppc = ppc
   else
      tree%ppc = PartInCell  ! default
   end if

   fmt="(A,I4,A)"
   write(str,fmt) "attempting to build oct-tree with ", &
                        tree%ppc, " particles per cell:"
   call mywrite(str,verb-1) 
   call mywrite("",verb-1)

   tree_storage_factor = 2.5/tree%ppc

   tree_allocation: do 
      tree%maxcells = tree_storage_factor* size(psys%par)
      call maketree(psys,tree,MBalloc)
      call inittree(psys,tree)  
      call parttree(psys,tree,err)
      if(err /= 0) then
         call killtree(tree)
         call mywrite("   woops, guessed wrong for tree memory, increasing allocation",verb)
         tree_storage_factor = 1.414 * tree_storage_factor
      else
         exit tree_allocation
      end if
   end do tree_allocation

   call patchtree(tree)
   call makecellorder(tree)
   call shrinkcells2(psys,tree)

   fmt="(A,F8.4)"
   write(str,fmt) "   allocated tree cells / actual tree cells =  ", &
                          real(tree%maxcells)/tree%ncells
   call mywrite(str,verb) 
   call mywrite("",verb)
   
!  call adjustbox(psys%box, tree%cell(1)%bot, tree%cell(1)%top)
  
 end subroutine buildtree

!-----------------------------------------------
!> allocates the various parts of the oct-tree
subroutine maketree(psys,tree,MBalloc)

  character(clen), parameter :: myname="maketree"
  logical, parameter :: crash=.true.
  integer, parameter :: verb=2
  character(clen) :: str,fmt
  
  type(particle_system_type), intent(in) :: psys !< particle system
  type(oct_tree_type), intent(inout)   :: tree !< oct-tree
  real(r8b), intent(out) :: MBalloc !< MB of memory allocated for tree
  
  integer(i4b) :: err
  real(r8b) :: MB


  MBalloc = 0.0d0
  
  if (allocated(tree%partorder)) call killtree(tree)
  if (tree%maxcells == 0) call myerr(" maxcells=0",myname,crash)
  tree%np = size(psys%par)
    
  MB = tree%np * 4 / 2**20
  str = "cant allocate oct-tree%partorder [npar]"
  allocate(tree%partorder(1:tree%np),stat=err)
  if(err.ne.0) call myerr(str,myname,crash)     
  MBalloc = MBalloc + MB

  MB = max(tree%np,tree%maxcells) * 4 / 2**20
  str = "cant allocate oct-tree%cellorder max[npar,maxcells]"
  allocate(tree%cellorder(1:max(tree%np,tree%maxcells)),stat=err)
  if(err.ne.0) call myerr(str,myname,crash)
  MBalloc = MBalloc + MB     
  
  MB = (tree%maxcells+1) * bytespercell / 2**20
  str = "cant allocate oct-tree%cell [maxcells] "
  allocate(tree%cell(0:tree%maxcells),stat=err)
  if(err.ne.0) call myerr(str,myname,crash)     
  MBalloc = MBalloc + MB     

  fmt = "(A,F12.4,A)"
  write(str,fmt) "   allocated ", MBalloc, " MB for tree"
  call mywrite(str,verb) 
  
end subroutine maketree

!-----------------------------------------------
!> deallocates the various parts of the oct-tree
 subroutine killtree(tree)
   type(oct_tree_type) :: tree !< oct-tree
   integer(i4b) :: err !< deallocation error number
   character(clen) :: string  !< error string
     deallocate(tree%partorder, tree%cellorder, tree%cell, stat=err)
     string = "deallocation      "
     if(err.ne.0) call error(string,err)
     tree%np=0
 end subroutine

!--------------------------------------------------------------
!> determines which octant the vector ppos - pcell points into
function get_octant(pcell,ppos) result(octant)
  integer, parameter,dimension(3) :: nindex=(/ 1,2,4 /) 
  integer :: octant    !< 1-8 labeling the octant (ppos-pcell) points into.
  real :: pcell(3)     !< position of the center of the cell.
  real(r4b) :: ppos(3) !< an arbitrary point in space.
  integer :: k 
  octant = 1
  do k = 1, 3
     if (ppos(k).GT.pcell(k)) octant = octant + nindex(k)
  enddo
end function get_octant

!---------------------------------------------
!> dummy checks the particles per leaf number
subroutine setPPC(tree,ppc)
  integer ::ppc !< particles/leaf to this number if the current one is crazy
  type(oct_tree_type) :: tree !< oct-tree
  if(ppc.LT.1.OR.ppc.GT.100000) call error("ppc error")
  tree%ppc=ppc
end subroutine setPPC


  
!----------------------------------------------------------------------
!> initializes the tree variables for the zeroth, first, and last cells
 subroutine inittree(psys,tree)

   type(particle_system_type), intent(in) :: psys !< particle system
   type(oct_tree_type), intent(inout) :: tree !< oct-tree
   integer :: npar,i
   real :: bot(3),top(3),rsize

     npar = size(psys%par)

     if (npar <= 0) call error("no particles in psys   ")
     if (npar /= tree%np) call error("npar in psys /equal to npar in tree")  

     tree%cell(0)%start = npar+1
     tree%cell(1)%start = 1
     tree%cell(1)%next  = 0
     tree%cell(1)%daughter = 0
     tree%ncells=1
     tree%nplastcell=npar
 
     do i=1,3
        bot(i)=minval(psys%par(1:npar)%pos(i))
        top(i)=maxval(psys%par(1:npar)%pos(i))
     enddo
 
     rsize=maxval(top-bot)
     tree%cell(1)%bot(1:3)=0.5*(top+bot)-0.5*rsize
     tree%cell(1)%top(1:3)=tree%cell(1)%bot(1:3)+rsize
 
     do i=1,npar
        tree%partorder(i)=i
     enddo
     
 end subroutine inittree


!-------------------------------------------------------------------
!> recursive subroutine that constructs the rest of the tree cells
 recursive subroutine parttree(psys,tree,err)
 
   type(oct_tree_type) :: tree          !< oct-tree being built
   type(particle_system_type) :: psys   !< particle system
   integer(i4b) :: err                       !< error number
   integer, parameter :: off(24) = &
       (/ 0,0,0, 1,0,0, 0,1,0, 1,1,0, 0,0,1, 1,0,1, 0,1,1, 1,1,1 /)
   integer :: offset_factor(3,8)
                
   real :: ccpos(3),csize(3)
   integer :: j,subp(0:nsubcell+1),subind
   integer :: ccell,nppcell,newcell,sister
   
     err=0   
     offset_factor = reshape( off, (/3,8/) )
     ccell=tree%ncells
     nppcell=tree%nplastcell
     if(nppcell.GT.tree%ppc) then
        ccpos=(tree%cell(ccell)%top+tree%cell(ccell)%bot)/2
        csize=tree%cell(ccell)%top-tree%cell(ccell)%bot
        subp=0
        do j=tree%cell(ccell)%start,tree%cell(ccell)%start+nppcell-1
           subind=get_octant(ccpos,psys%par(tree%partorder(j))%pos)    
           subp(subind)=subp(subind)+1
        enddo
        subp(0)=0
        subp(nsubcell+1)=nppcell

        do j=nsubcell,1,-1
           subp(j)=subp(j+1)-subp(j) 
        enddo
        if(subp(1).NE.0) call error("subp count ")

        do j = tree%cell(ccell)%start, tree%cell(ccell)%start + nppcell - 1
           subind=get_octant(ccpos,psys%par(tree%partorder(j))%pos)
           subp(subind)=subp(subind)+1
           tree%cellorder(subp(subind))=tree%partorder(j)
        enddo

        do j=tree%cell(ccell)%start,tree%cell(ccell)%start+nppcell-1
           tree%partorder(j)=tree%cellorder(j-tree%cell(ccell)%start+1)
        enddo
        sister=0
        do j=1,nsubcell
           if(subp(j)-subp(j-1).GT.0) then 
              tree%ncells=tree%ncells+1
              newcell=tree%ncells

              if(newcell.GT.tree%maxcells) then
                err=1
                return
              end if

              tree%nplastcell=subp(j)-subp(j-1)

              if(tree%cell(ccell)%daughter.EQ.0) then
                 tree%cell(ccell)%daughter=newcell
              end if

              tree%cell(newcell)%start=subp(j-1)+tree%cell(ccell)%start
              tree%cell(newcell)%next=ccell
              tree%cell(newcell)%daughter=0
              tree%cell(newcell)%bot=tree%cell(ccell)%bot + &
                                     offset_factor(1:3,j)*csize/2
              tree%cell(newcell)%top=tree%cell(newcell)%bot+csize/2
              tree%cell(sister)%next=newcell
              sister=newcell
              call parttree(psys,tree,err)
              if(err.ne.0) return
           endif
        enddo
     endif
 end subroutine parttree

!------------------------------------------------------------------------------
!> fixes the cell pointers using the information from the fully built tree
 subroutine patchtree(tree)
   type(oct_tree_type) :: tree !> oct-tree being built
   integer ::i
     tree%cell(0)%next=0
     do i=1,tree%ncells
        if(tree%cell(i)%next < i) &
             tree%cell(i)%next=tree%cell(tree%cell(i)%next)%next
     enddo
 end subroutine patchtree


!----------------------------------
!> constructs the cell order array
 subroutine makecellorder(tree)
   type(oct_tree_type) :: tree !< oct-tree being built
   integer ::  ccount,dcell,i,up,pcell
   integer(i8b) :: low
   integer(i4b) :: err

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

     if(low.NE.1) then
        err=low
        call error("ordering ",err)
     endif

 end subroutine makecellorder

 !---------------------------------------------
 !> returns the number of particles in cell i
 function ppcell(tree,i) result (parincell)
   integer :: parincell  !< the number of particles in cell i 
   type(oct_tree_type) :: tree   !< the oct-tree 
   integer ::i !< the index of the cell to search
   parincell=tree%cell(tree%cell(i)%next)%start-tree%cell(i)%start
 end function

 !-----------------------------------------------------------
 !> makes the cell boundaries tight wrt to particle positions
 subroutine shrinkcells(psys,tree)
   type(particle_system_type), intent(in) :: psys !< particle system
   type(oct_tree_type), intent(out) :: tree   !< oct-tree being built
   integer :: i0,i,j,next
   real ::  x(3)
  
     do i0=1,tree%ncells
        i=tree%cellorder(i0)
        x=tree%cell(i)%top
        tree%cell(i)%top = tree%cell(i)%bot
        tree%cell(i)%bot = x
        next=tree%cell(i)%daughter

        if(next.EQ.0) then
           do j=tree%cell(i)%start,tree%cell(tree%cell(i)%next)%start-1
              x = psys%par(tree%partorder(j))%pos
              where( tree%cell(i)%bot .gt. x ) tree%cell(i)%bot=x
              where( tree%cell(i)%top .lt. x ) tree%cell(i)%top=x
           enddo
        else
           do while(next.NE.tree%cell(i)%next)
              x = tree%cell(next)%bot
              where( tree%cell(i)%bot .gt. x) tree%cell(i)%bot=x
              x = tree%cell(next)%top
              where( tree%cell(i)%top .lt. x) tree%cell(i)%top=x
              next = tree%cell(next)%next
           enddo
        endif
     enddo
 end subroutine shrinkcells

 !-----------------------------------------------------------
 !> makes the cell boundaries tight wrt to particle positions
 !! and the AABB boundaries tight wrt the smoothing lengths
 subroutine shrinkcells2(psys,tree)
   type(oct_tree_type) :: tree !< the oct-tree being built
   type(particle_system_type) :: psys !< particle system
   integer :: i0,i,j,next
   real(r4b) :: x(3),hsml
   character(clen) :: string
   
     do i0 = 1,tree%ncells
        i = tree%cellorder(i0)
        x = tree%cell(i)%top
        tree%cell(i)%top=tree%cell(i)%bot
        tree%cell(i)%bot=x
        tree%cell(i)%botrange=tree%cell(i)%bot
        tree%cell(i)%toprange=tree%cell(i)%top
        next=tree%cell(i)%daughter

        if(next.EQ.0) then
           if(tree%cell(tree%cell(i)%next)%start-tree%cell(i)%start.LT.1) then
              string = "zero part. cell?    "
              call error(string)
           end if
           do j=tree%cell(i)%start,tree%cell(tree%cell(i)%next)%start-1
              hsml=psys%par(tree%partorder(j))%hsml
              x(1:3)=psys%par(tree%partorder(j))%pos(1:3)

!             do these work correctly for each element of the array? yes
              tree%cell(i)%bot=min(tree%cell(i)%bot,x)
              tree%cell(i)%top=max(tree%cell(i)%top,x)
              tree%cell(i)%botrange=min(tree%cell(i)%botrange,x-hsml)
              tree%cell(i)%toprange=max(tree%cell(i)%toprange,x+hsml)
           enddo
        else
           do while(next.NE.tree%cell(i)%next)
              tree%cell(i)%bot=min(tree%cell(i)%bot,tree%cell(next)%bot)
              tree%cell(i)%top=max(tree%cell(i)%top,tree%cell(next)%top)
              tree%cell(i)%botrange= &
                   min(tree%cell(i)%botrange,tree%cell(next)%botrange)
              tree%cell(i)%toprange= &
                   max(tree%cell(i)%toprange,tree%cell(next)%toprange)
              next=tree%cell(next)%next
           enddo
        endif
     enddo
 end subroutine shrinkcells2


!--------------------------------------------------------------
!> returns the linear extent of the tree leaf that contains pos
 function getlocalscale(pos,tree) result(scale)

   real :: pos(3),scale
   type(oct_tree_type) :: tree
   integer :: next,this,daughter,old,end

   scale=sqrt(dist2cell(pos,tree%cell(1)))
   if(scale.GT.0) return
   old=1
   this=1
   next=1
   end=0
   do while(next.NE.end)
    this=next
    if(dist2cell(pos,tree%cell(this)).GT.0) then
     next=tree%cell(this)%next
!     print*,'a',next,dist2cell(pos,tree%cell(this))
    else
     daughter=tree%cell(this)%daughter
     if(daughter.EQ.0) exit
     old=this
     end=tree%cell(this)%next
!     print*,'b',daughter,dist2cell(pos,tree%cell(this))
     next=daughter
    endif       
   enddo

   scale=sum(tree%cell(old)%top-tree%cell(old)%bot)/3/2
   if(scale.LE.0) then
    call error('getlocalscale error')
   endif

!   print*,scale
!   print*,daughter,next,old
!   print*,tree%cell(old)%top,tree%cell(old)%bot
!   stop
 end function getlocalscale
 
!> calculates the distance squared from a point to the nearest of the 
!! cell walls.  returns zero for points inside the cell.
 function dist2cell(pos,cell) result(dist2)
  type(cell_type) :: cell
  real :: pos(3),lpos(3),upos(3),dist2
  integer :: k
  dist2=0
  lpos=cell%bot
  upos=cell%top
  do k=1,3
     if(pos(k).LT.lpos(k)) then
        dist2=dist2+(lpos(k)-pos(k))**2
     else
        if(pos(k).GT.upos(k)) then
           dist2=dist2+(upos(k)-pos(k))**2
        endif
     endif
  enddo

 end function dist2cell

!> a simple error handeling routine
 subroutine error(string,i)
   character(*) :: string !< error string
   integer(i4b), optional :: i !< optional error number
  
    print*,'error detected'

    if(present(i)) then
       print*,string,i
    else
       print*,string
    endif
    stop

 end subroutine error


end module oct_tree_mod
