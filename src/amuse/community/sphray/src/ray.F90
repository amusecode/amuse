!> \file ray.F90

!> \brief the ray module 
!!
!<

module ray_mod
use myf03_mod
use particle_system_mod
use oct_tree_mod
use mt19937_mod, only: genrand_real1
use spectra_mod, only: rn2freq
use physical_constants_mod, only: HI_th_erg, M_H, M_He
implicit none

private

public :: src_ray_type
public :: src_ray_make
public :: src_ray_transform

public :: src_ray_class_from_dir
public :: src_ray_dist2pt
public :: src_ray_part_intersection
public :: src_ray_cell_intersection
public :: src_ray_pluecker

public :: raystatbuffsize
public :: raystat_type



integer(i8b) :: raystatbuffsize 
real(r8b), parameter :: zero = 0.0d0
real(r8b), parameter :: one = 1.0d0

  

!> basic ray type
!-----------------------------------------------------------------------
type src_ray_type
   real(r8b) :: start(3)  !< starting position
   real(r8b) :: dir(3)    !< unit vector direction
   real(r8b) :: length    !< length (determines when to stop tracing)
   integer(i4b) :: class  !< based on direction signs (MMM, PMM, ...)
   
   real(r8b) :: freq      !< freq in HI ionizing units
   real(r8b) :: enrg      !< enrg of a single photon in ergs
   real(r8b) :: pcnt      !< photon count (changes as the ray is depleted)
   real(r8b) :: pini      !< initial photons
   real(r8b) :: dt_s      !< time step associated with ray [s]
end type src_ray_type


!> ray stats that can be output for each ray
!-----------------------------------------------------------------------
type raystat_type
   integer(i4b) :: srcn
   real(r4b) :: start(3)
   real(r4b) :: ryd
end type raystat_type




contains




!> creates a source ray 
!-----------------------------------------------------------------------  
  subroutine src_ray_make(ray, src, rayn, dtray_s, Lunit, box, length)

    type(src_ray_type), intent(out) :: ray    !< ray to make
    type(source_type), intent(inout) :: src   !< source
    integer(i8b), intent(in) :: rayn          !< ray indx
    real(r8b), intent(in) :: dtray_s          !< time between rays [s]   
    real(r8b), intent(in) :: Lunit            !< converts src%lum -> photons/s
    type(box_type), intent(in) :: box         !< simulation box

    real(r8b), intent(in), optional :: length !< optional length (default=huge)

    real(r8b) :: xx,yy,zz,r
    real(r8b) :: rn1, rn2
    real(r8b) :: prate
    integer :: i

  
!  set the direction of the ray from the emmission profile (src%EmisPrf)
!     0  = isotropic
!    -1  = towards +z
!    -2  = towards -z
!
!  note that for all point sources the luminosity is interpreted as a Flux 
!  [photons/s].  

    select case (src%EmisPrf)
       
    ! this makes rays go in -z direction  
    !-----------------------------------------------------------------------  
    case(-2)

       ray%dir(1) = 0.0
       ray%dir(2) = 0.0
       ray%dir(3) = -1.0       
       
    ! this makes rays go in +z direction 
    !-----------------------------------------------------------------------  
    case(-1)

       ray%dir(1) = 0.0
       ray%dir(2) = 0.0
       ray%dir(3) = 1.0
                 
    ! random direction on the unit sphere
    !-----------------------------------------------------------------------  
    case(0)

       r=2.0d0 
       do while ( r .GT. 1.0d0 .and. r .NE. 0.0d0 )
          xx=(2.0d0 * genrand_real1()-1.0d0)   
          yy=(2.0d0 * genrand_real1()-1.0d0)   
          zz=(2.0d0 * genrand_real1()-1.0d0)   
          r=xx*xx+yy*yy+zz*zz
       enddo
       r = sqrt(r)
       ray%dir(1) = xx/r  ! it is important that ray%dir be a unit vector
       ray%dir(2) = yy/r
       ray%dir(3) = zz/r          
       
    case default

       write(*,*) "emission profile not recognized"
       write(*,*) "profile = ", src%EmisPrf
       stop          
       
    end select
    !-----------------------------------------------------------------------  

    ray%start = src%pos

    if ( present(length) ) then
       ray%length = length
    else
       ray%length = huge(1.0d0) * 0.1d0
    endif

  
!   set the class of the ray (what octant is it going into)
    call src_ray_class_from_dir(ray)


!   set the frequency and energy / photon of the ray
    ray%freq = rn2freq(src%SpcType)
    ray%enrg = ray%freq * HI_th_erg  


!   set the number of photons in the ray
    if (ray%enrg > 0.) then
       prate = src%L * Lunit ! photons per sec
    else 
       prate = 0.
    end if
    if (rayn > src%lastemit) then
       ray%dt_s = dtray_s * (rayn - src%lastemit) 
       ray%pini = prate * ray%dt_s
    else
       write(*,*) "make_source_ray> rayn .LT. src%lastemit in ray.f90"
       stop
    end if
    ray%pcnt = ray%pini
    src%lastemit = rayn

  end subroutine src_ray_make


!> returns a transformed ray while preserving the initial ray
!--------------------------------------------------------------
  subroutine src_ray_transform(inray, outray, trans)
    type(src_ray_type) :: inray       !< input ray
    type(src_ray_type) :: outray       !< output ray
    type(transformation_type) :: trans !< transformation
    
    outray%start  = inray%start * trans%fac + trans%shift
    outray%dir    = inray%dir   
    outray%length = inray%length  
    outray%class  = inray%class

    outray%freq = inray%freq
    outray%enrg = inray%enrg
    outray%pcnt = inray%pcnt
    outray%pini = inray%pini
    outray%dt_s = inray%dt_s

  end subroutine src_ray_transform



! pre computes the class of the ray for the Pluecker test
! ray label    class
!   MMM          0
!   PMM          1
!   MPM          2
!   PPM          3
!   MMP          4
!   PMP          5
!   MPP          6
!   PPP          7
!-----------------------------------------------------------
subroutine src_ray_class_from_dir( src_ray ) 
  type(src_ray_type), intent(inout) :: src_ray
  integer(i4b) :: i

  src_ray%class = 0
  do i = 1, 3
     if ( src_ray%dir(i) >= zero ) src_ray%class = src_ray%class + 2**(i-1)
  end do
  
end subroutine src_ray_class_from_dir






!> returns the perpendicular distance between a point and a ray
!--------------------------------------------------------------  
function src_ray_dist2pt(src_ray, pt, proj) result(perp)
  type(src_ray_type), intent(in) :: src_ray
  real(r4b), intent(in) :: pt(3)     !< point
  real(r8b), optional :: proj        !< distance along ray 
  real(r8b) :: perp                  !< distance perp. to ray
  
  real(r8b) :: dotp
  real(r8b) :: diff(3)
  real(r8b) :: perp2
  real(r8b) :: vec1(3)
  real(r8b) :: vec2(3)
  
  vec1 = pt - src_ray%start
  vec2 = src_ray%dir
  
  dotp  = dot_product( vec1, vec2 )
  diff  = vec1 - dotp * src_ray%dir 
  perp2 = dot_product( diff, diff )
  perp = sqrt(perp2)
  
  if( present(proj) ) proj = dotp
  
end function src_ray_dist2pt



!> tests for ray / particle intersection.  
!----------------------------------------------  
function src_ray_part_intersection(src_ray, part) result(hit)
  type(src_ray_type), intent(in) :: src_ray
  type(particle_type), intent(in) :: part       !< particle
  logical :: hit                    !< true or false result
  
  real(r8b) :: start2cen       !< distance^2 from ray start to part position
  real(r8b) :: end2cen         !< distance^2 from ray end to part position
  real(r8b) :: perp            !< perpendicular distance to particle 
  real(r8b) :: proj            !< projected distance along ray
  real(r8b) :: diff(3)         !< vector from ray start to part center
  
  ! if the perpendicular distance to the point is larger than hsml exit
  !---------------------------------------------------------------------
  perp = src_ray_dist2pt( src_ray, part%pos, proj ) 
  if (perp >= part%hsml) then
     hit = .false.
     return
  endif
  
  ! now our only concern is the position of the particle
  ! along the ray.  


  ! first we reject particles that cannot possibly be intersections
  !
  if ( proj < -part%hsml ) then
     hit = .false. 
     return
  endif

  if ( proj > src_ray%length + part%hsml ) then
     hit = .false. 
     return
  endif


  ! if the particle smoothing volume could contain the origin of 
  ! the ray, test for that.
  !---------------------------------------------------------------------
  if ( abs(proj) < part%hsml ) then
     diff = part%pos - src_ray%start
     start2cen = sum( diff*diff ) 
     if (start2cen < part%hsml*part%hsml) then
        hit = .true.
        return
     endif
  endif

  ! if the particle smoothing volume could contain the terminus of 
  ! the ray, test for that.
  !---------------------------------------------------------------------
  if ( abs(proj - src_ray%length) < part%hsml ) then
     diff = part%pos - (src_ray%start + src_ray%length * src_ray%dir)
     end2cen = sum( diff*diff )
     if (end2cen < part%hsml*part%hsml) then
        hit = .true.
        return
     endif
  endif

  ! at this point we must have a hit
  hit = .true.
  return

  
end function src_ray_part_intersection


!> tests for ray - AABB intersection. 
!----------------------------------------------  
function src_ray_cell_intersection(src_ray,cell) result(hit)
  type(src_ray_type), intent(in) :: src_ray  !< ray
  type(cell_type), intent(in) :: cell        !< cell
  logical :: hit                 !< true or false result
  real(r8b) :: bot(3)
  real(r8b) :: top(3)
  bot = cell%botrange - src_ray%start
  top = cell%toprange - src_ray%start
  hit = src_ray_pluecker(src_ray, bot, top)
end function src_ray_cell_intersection



!> pluecker test for line segment / cell intersection
!-----------------------------------------------------    
function src_ray_pluecker(src_ray, s2b, s2t) result( hit )

  type(src_ray_type), intent(in) :: src_ray

  real(r8b) :: s2b(3)         !< vector from ray start to lower cell corner
  real(r8b) :: s2t(3)         !< vector from ray start to upper cell corner
  logical :: hit              !< true or false result
  
  real(r8b) :: dir(3)
  real(r8b) :: dist

  real(r8b) :: e2b(3)       !< vector from ray end to lower cell corner
  real(r8b) :: e2t(3)       !< vector from ray end to upper cell corner

  dir  = src_ray%dir  
  dist = src_ray%length

  e2b = s2b - dir * dist
  e2t = s2t - dir * dist

  hit = .false.

  ! branch on ray direction
  !---------------------------
  select case( src_ray%class )

     ! MMM
     !-----------
  case(0)

     if(s2b(1) > zero .or. s2b(2) > zero .or. s2b(3) > zero) return ! on negative part of ray 
     if(e2t(1) < zero .or. e2t(2) < zero .or. e2t(3) < zero) return ! past length of ray      

     if ( dir(1)*s2b(2) - dir(2)*s2t(1) < zero .or.  &
          dir(1)*s2t(2) - dir(2)*s2b(1) > zero .or.  &
          dir(1)*s2t(3) - dir(3)*s2b(1) > zero .or.  &
          dir(1)*s2b(3) - dir(3)*s2t(1) < zero .or.  &
          dir(2)*s2b(3) - dir(3)*s2t(2) < zero .or.  &
          dir(2)*s2t(3) - dir(3)*s2b(2) > zero       ) return
     
     ! PMM
     !-----------
  case(1)
     
     if(s2t(1) < zero .or. s2b(2) > zero .or. s2b(3) > zero) return ! on negative part of ray 
     if(e2b(1) > zero .or. e2t(2) < zero .or. e2t(3) < zero) return ! past length of ray      
     
     if ( dir(1)*s2t(2) - dir(2)*s2t(1) < zero .or.  &
          dir(1)*s2b(2) - dir(2)*s2b(1) > zero .or.  &
          dir(1)*s2b(3) - dir(3)*s2b(1) > zero .or.  &
          dir(1)*s2t(3) - dir(3)*s2t(1) < zero .or.  &
          dir(2)*s2b(3) - dir(3)*s2t(2) < zero .or.  &
          dir(2)*s2t(3) - dir(3)*s2b(2) > zero       ) return
     
     ! MPM
     !-----------
  case(2)
     
     if(s2b(1) > zero .or. s2t(2) < zero .or. s2b(3) > zero) return ! on negative part of ray 
     if(e2t(1) < zero .or. e2b(2) > zero .or. e2t(3) < zero) return ! past length of ray      
     
     if ( dir(1)*s2b(2) - dir(2)*s2b(1) < zero .or.  &
          dir(1)*s2t(2) - dir(2)*s2t(1) > zero .or.  &
          dir(1)*s2t(3) - dir(3)*s2b(1) > zero .or.  &
          dir(1)*s2b(3) - dir(3)*s2t(1) < zero .or.  &
          dir(2)*s2t(3) - dir(3)*s2t(2) < zero .or.  &
          dir(2)*s2b(3) - dir(3)*s2b(2) > zero       ) return
     
     ! PPM
     !-----------
  case(3)
     
     if(s2t(1) < zero .or. s2t(2) < zero .or. s2b(3) > zero) return ! on negative part of ray 
     if(e2b(1) > zero .or. e2b(2) > zero .or. e2t(3) < zero) return ! past length of ray      
     
     if ( dir(1)*s2t(2) - dir(2)*s2b(1) < zero .or.  &
          dir(1)*s2b(2) - dir(2)*s2t(1) > zero .or.  &
          dir(1)*s2b(3) - dir(3)*s2b(1) > zero .or.  &
          dir(1)*s2t(3) - dir(3)*s2t(1) < zero .or.  &
          dir(2)*s2t(3) - dir(3)*s2t(2) < zero .or.  &
          dir(2)*s2b(3) - dir(3)*s2b(2) > zero       ) return
     
     ! MMP
     !-----------
  case(4)
     
     if(s2b(1) > zero .or. s2b(2) > zero .or. s2t(3) < zero) return ! on negative part of ray 
     if(e2t(1) < zero .or. e2t(2) < zero .or. e2b(3) > zero) return ! past length of ray      
     
     if ( dir(1)*s2b(2) - dir(2)*s2t(1) < zero .or.  &
          dir(1)*s2t(2) - dir(2)*s2b(1) > zero .or.  &
          dir(1)*s2t(3) - dir(3)*s2t(1) > zero .or.  &
          dir(1)*s2b(3) - dir(3)*s2b(1) < zero .or.  &
          dir(2)*s2b(3) - dir(3)*s2b(2) < zero .or.  &
          dir(2)*s2t(3) - dir(3)*s2t(2) > zero       ) return
     
     
     ! PMP
     !-----------
  case(5)
     
     if(s2t(1) < zero .or. s2b(2) > zero .or. s2t(3) < zero) return ! on negative part of ray 
     if(e2b(1) > zero .or. e2t(2) < zero .or. e2b(3) > zero) return ! past length of ray      
     
     if ( dir(1)*s2t(2) - dir(2)*s2t(1) < zero .or.  &
          dir(1)*s2b(2) - dir(2)*s2b(1) > zero .or.  &
          dir(1)*s2b(3) - dir(3)*s2t(1) > zero .or.  &
          dir(1)*s2t(3) - dir(3)*s2b(1) < zero .or.  &
          dir(2)*s2b(3) - dir(3)*s2b(2) < zero .or.  &
          dir(2)*s2t(3) - dir(3)*s2t(2) > zero       ) return
     
     
     ! MPP
     !-----------
  case(6)
     
     if(s2b(1) > zero .or. s2t(2) < zero .or. s2t(3) < zero) return ! on negative part of ray 
     if(e2t(1) < zero .or. e2b(2) > zero .or. e2b(3) > zero) return ! past length of ray      
     
     if ( dir(1)*s2b(2) - dir(2)*s2b(1) < zero .or.  &
          dir(1)*s2t(2) - dir(2)*s2t(1) > zero .or.  &
          dir(1)*s2t(3) - dir(3)*s2t(1) > zero .or.  &
          dir(1)*s2b(3) - dir(3)*s2b(1) < zero .or.  &
          dir(2)*s2t(3) - dir(3)*s2b(2) < zero .or.  &
          dir(2)*s2b(3) - dir(3)*s2t(2) > zero       ) return
     
     ! PPP
     !-----------
  case(7)
     
     if(s2t(1) < zero .or. s2t(2) < zero .or. s2t(3) < zero) return ! on negative part of ray 
     if(e2b(1) > zero .or. e2b(2) > zero .or. e2b(3) > zero) return ! past length of ray      
     
     if ( dir(1)*s2t(2) - dir(2)*s2b(1) < zero .or.  &
          dir(1)*s2b(2) - dir(2)*s2t(1) > zero .or.  &
          dir(1)*s2b(3) - dir(3)*s2t(1) > zero .or.  &
          dir(1)*s2t(3) - dir(3)*s2b(1) < zero .or.  &
          dir(2)*s2t(3) - dir(3)*s2b(2) < zero .or.  &
          dir(2)*s2b(3) - dir(3)*s2t(2) > zero       ) return
     
  case default
     call rayError('ray class.')
     
  end select
  
  hit=.true.
  
end function src_ray_pluecker



!> error handling
!-----------------------------      
  subroutine rayError(string,i)
    character(*) :: string  !< error string
    integer, optional :: i  !< error number
    
    print*,' Error detected:'
    
    if(present(i)) then
       print*,string,i
    else
       print*,string
    endif
    
    stop
  end subroutine rayError

  
end module ray_mod
