!> \file raylist.F90

!> \brief The raylist module 
!!
!<

module raylist_mod
use myf03_mod
use ray_mod
use particle_system_mod, only: particle_system_type
use particle_system_mod, only: particle_type
use particle_system_mod, only: transformation_type
use oct_tree_mod, only: oct_tree_type

implicit none

 private
 public :: intersection_type
 public :: trace_ray
 public :: raylist_type
 public :: prepare_raysearch
 public :: kill_raylist

 real, parameter :: zero = 0.0d0
 real, parameter :: half = 0.5d0

 integer,parameter :: MAX_RAYLIST_LENGTH = 1000000  !< default maximum
 integer,parameter :: ndim = 3                      !< number of dimensions
 integer,parameter :: nimages = 3**ndim             !< number of search images


!> holds a particle index, an impact parameter, and a distance along a ray
!! usefull for the ionization and temperature updating
! -----------------------------------------------------------------------
   type intersection_type
      integer(i8b) :: pindx   !< particle index
      real :: b          !< impact parameter
      real :: d          !< distance along ray
      real :: dl         !< path length
   end type intersection_type
 
!> grouping of all things ray + impacts
!--------------------------------------- 
   type raylist_type
      integer :: nnb              !< number of intersections
      integer :: maxnnb           !< maximum number of intersections
      integer :: lastnnb          !< intersection where we stopped
      integer(i8b) :: searchcell  !< index of cell being searched
      logical :: reuseable        !< is this ray reusable?
      integer :: searchimage      !< which transformation of the particles?
      integer :: nsearchimages    !< how many images to search
      type(src_ray_type) :: ray   !< ray 
      type(transformation_type) :: trafo(nimages)  !< transformations  
      type(intersection_type), allocatable :: intersection(:) !< ray/par 
   end type raylist_type



contains


 
!> set intersection values
!-----------------------------------------
 function set_intersection(ray, part, i) result(t)
   type(intersection_type) :: t  !< the intersection
   type(src_ray_type) :: ray         !< the ray
   type(particle_type) :: part   !< the particle
   integer(i8b)  :: i            !< the particle index
   real(r8b) :: d
   
     t%pindx = i
     t%b = src_ray_dist2pt(ray, part%pos, d)
     t%d = d

 end function set_intersection

 
!> initialize raylist variables and set search images
!------------------------------------------------------
 subroutine prepare_raysearch(psys, raylist)

   type(particle_system_type) psys !< particle system
   type(raylist_type) raylist      !< ray list
  
   call make_raylist(MAX_RAYLIST_LENGTH, raylist)
   call setsearchimages(psys,raylist)

 end subroutine prepare_raysearch


!> check all of the search images for intersection
!---------------------------------------------------
 subroutine fullsearch(psys, searchtree, raylist)

   type(particle_system_type) :: psys  !< particle system
   type(oct_tree_type) :: searchtree   !< oct-tree to search
   type(raylist_type) raylist          !< raylist

   if(raylist%searchimage == 0) call raylistError('raylist init.')
   do while (raylist%searchimage <= raylist%nsearchimages)   
      call raysearch(psys, searchtree, raylist)
      if (raylist%searchcell /= 0) return 
      raylist%searchimage = raylist%searchimage + 1
      raylist%searchcell = 1
   enddo
   raylist%searchcell = 0

 end subroutine fullsearch


!> checks a single search image for intersection
!-----------------------------------------------
 subroutine raysearch(psys, searchtree, raylist)

   type(particle_system_type), intent(in) :: psys         !< particle system
   type(oct_tree_type), intent(in), target :: searchtree  !< oct-tree to search
   type(raylist_type) :: raylist                          !< raylist


   type(oct_tree_type), pointer :: tree      !< pointer to tree
   type(src_ray_type) :: curay               !< transformed ray
   integer(i8b) :: this, daughter, next      !< octree indices
   integer(i8b) :: par_in_cell               !< number of particles in current search cell
   logical :: par_hit                        !< ray / particle intersection test
   logical :: cell_hit                       !< ray / cell intersection test
   logical :: long                           !< past max distance? 
   integer(i8b) :: i, si, orderindx

   tree => searchtree
   si = raylist%searchimage

   ! curay%start = ray%start * fac + shift
   call src_ray_transform( raylist%ray, curay, raylist%trafo(si) )
   next = raylist%searchcell

   do while (next /= 0)

      this     = next
      daughter = tree%cell(this)%daughter
      next     = tree%cell(this)%next
           
      ! if we've reached a leaf
      !----------------------------
      if (daughter == 0) then

         ! return if we go over max intersections
         !--------------------------------------------
         par_in_cell = tree%cell(next)%start - tree%cell(this)%start
         if (raylist%nnb + par_in_cell > raylist%maxnnb) then             
            write(*,*) ' *** reached max intersections *** '
            raylist%reuseable  = .false.
            raylist%searchcell = this
            return            
         endif
         
         ! add intersected particles to list
         !----------------------------------------------
         do i = tree%cell(this)%start, tree%cell(next)%start - 1
            orderindx = tree%partorder(i)
            par_hit = src_ray_part_intersection( curay, psys%par(orderindx) ) 
            if (par_hit) then
               raylist%nnb = raylist%nnb + 1
               raylist%intersection(raylist%nnb) = &
                    set_intersection(curay, psys%par(orderindx), orderindx)
            endif
         enddo


      ! if we need to descend further
      !--------------------------------
      else

         cell_hit = src_ray_cell_intersection( curay, tree%cell(this) )
         if ( cell_hit ) next = daughter  

      endif

   enddo
   
   raylist%searchcell = next
   
 end subroutine raysearch


!> initializes values in the raylist
!-------------------------------------------
 subroutine make_raylist(maxnnb, raylist, ray)

   integer, intent(in) :: maxnnb       !< maximum number of intersections
   type(raylist_type)      :: raylist  !< the raylist
   type(src_ray_type),optional :: ray  !< the ray
   real(r8b) :: start(ndim)
   real(r8b) :: dir(ndim)
  
     raylist%nnb            = 0
     raylist%maxnnb         = maxnnb
     raylist%searchcell     = 1
     raylist%reuseable      = .false.
     raylist%searchimage    = 1
     raylist%nsearchimages  = 1
     raylist%trafo(1)%fac   = 1
     raylist%trafo(1)%shift = zero 

     allocate(raylist%intersection(maxnnb))

     if (present(ray)) then
        raylist%ray = ray
     endif

 end subroutine make_raylist

!> reset an already initialized raylist
!---------------------------------------
 subroutine reset_raylist(raylist,ray)

   type(raylist_type) :: raylist       !< the raylist
   type(src_ray_type), optional :: ray !< the ray
  
     raylist%nnb         = 0
     raylist%searchcell  = 1
     raylist%searchimage = 1
     raylist%reuseable   = .false.
     if (present(ray)) then
        raylist%ray = ray
     endif

 end subroutine reset_raylist

!> kill a raylist
!---------------------------------
 subroutine kill_raylist(raylist)

   type(raylist_type) :: raylist !< the raylist to kill
   real(r8b) :: start(ndim)
   real(r8b) :: dir(ndim) 
  
     raylist%nnb=0
     raylist%maxnnb=MAX_RAYLIST_LENGTH
     raylist%searchcell=0
     deallocate(raylist%intersection)

 end subroutine kill_raylist


!> create the transformations for the searchimages
!-------------------------------------------------
 subroutine setsearchimages(psys,raylist)

   type(particle_system_type), intent(in) :: psys
   type(raylist_type) :: raylist
   
   real :: top(ndim)        !< top box corner
   real :: bot(ndim)        !< bottom box corner
   integer :: bbound(ndim)  !< bottom BCs
   integer :: tbound(ndim)  !< top BCs
   integer :: nsearchimages !< number of valid periodic images

   integer :: i_image       !< loop over images counter
   integer :: i_dim         !< loop over dimensions counter
  
   integer :: vec(ndim)     !< vector from origin of base box to origin of image box
   integer :: fac           !< tmp storage transformation factor for each dimension
   real    :: shift         !< tmp storage transformation shift for each dimension

   integer :: l
   real :: lx(ndim),hx(ndim)
   logical :: boxexists

   top    = psys%box%tops
   bot    = psys%box%bots
   bbound = psys%box%bbound
   tbound = psys%box%tbound
   
   nsearchimages = 0
   
   do i_image = 0, nimages - 1
      l = i_image
      
      do i_dim = 1, ndim
         vec(i_dim)=mod(l,3)-1 
         l=l/3
      enddo
      
      boxexists = .true.
      
      do i_dim = 1, ndim
         lx(i_dim) = bot(i_dim) + vec(i_dim) * ( top(i_dim) - bot(i_dim) )
         hx(i_dim) = top(i_dim) + vec(i_dim) * ( top(i_dim) - bot(i_dim) )
         if( vec(i_dim) == -1 .and. bbound(i_dim) == 0) boxexists = .false.  
         if( vec(i_dim) ==  1 .and. tbound(i_dim) == 0) boxexists = .false.  
      enddo
            
      if (boxexists) then
         nsearchimages = nsearchimages + 1
         
         do i_dim = 1, ndim
            
            fac   = 1
            shift = zero
            
            if( vec(i_dim) == -1 ) then
               fac   = bbound(i_dim)
               shift = (-3 * fac + 1) * half * bot(i_dim) + &
                       (     fac + 1) * half * top(i_dim)
            endif
            
            if( vec(i_dim) == 1 ) then
               fac   = tbound(i_dim)
               shift = (-3 * fac + 1) * half * top(i_dim) + &
                       (     fac + 1) * half * bot(i_dim)
            endif

            raylist%trafo(nsearchimages)%fac(i_dim)   = fac
            raylist%trafo(nsearchimages)%shift(i_dim) = shift

         enddo
         
      endif ! boxexists
      
      
   enddo
   
   raylist%nsearchimages = nsearchimages


 end subroutine setsearchimages


!> error handling
!---------------------------------
 subroutine raylistError(string,i)
   character(*) :: string  !< error message
   integer, optional :: i  !< error number

     print*,' Error detected:'
  
     if(present(i)) then
        print*,string,i
     else
        print*,string
     endif
  
     stop

 end subroutine raylistError

!> sorts a raylist according to the distance along the ray with the particles 
!! closest to the origin of the ray first in the list.   The corresponding 
!! changes are made to pindx and b
!---------------------------------------------------------------------------
 subroutine sort3_raylist(raylist)
 use m_mrgrnk, only: mrgrnk
 implicit none

   type(raylist_type) :: raylist     !< raylist to sort

   integer(i8b) :: N 
   integer(i8b) :: indexx(raylist%nnb)
   real(r8b) :: darr(raylist%nnb)

   N = raylist%nnb
   darr(1:N)=raylist%intersection(1:N)%d
   call mrgrnk(darr,indexx)

   raylist%intersection(1:N)%d=raylist%intersection(indexx(1:N))%d  
   raylist%intersection(1:N)%b=raylist%intersection(indexx(1:N))%b 
   raylist%intersection(1:N)%pindx=raylist%intersection(indexx(1:N))%pindx  

 end subroutine sort3_raylist





!> given a ray creates a raylist with intersections
!------------------------------------------------------
 subroutine trace_ray(ray, raylist, psys, searchtree, dosort) 
   type(src_ray_type), intent(in) :: ray           !< the ray to trace
   type(raylist_type), intent(inout) :: raylist    !< the returned raylist
   type(particle_system_type), intent(in) :: psys  !< the particle system
   type(oct_tree_type), intent(in) :: searchtree   !< the oct-tree to search
   logical, intent(in), optional :: dosort         !< default is true

   logical :: wantsort

   wantsort = .true.
   if (present(dosort)) then
      if (.not. dosort) wantsort = .false.
   endif

   if (raylist%maxnnb /= MAX_RAYLIST_LENGTH) then      
      call prepare_raysearch(psys, raylist)
   end if

   call reset_raylist(raylist, ray)
   call fullsearch(psys, searchtree, raylist)

   if (wantsort) call sort3_raylist(raylist)

 end subroutine trace_ray




end module raylist_mod
