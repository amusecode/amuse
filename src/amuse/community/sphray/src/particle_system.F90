!> \file particle_system.F90

!> \brief Particles, sources, and box, types and subroutines.
!! 
!<

module particle_system_mod
use myf03_mod
use atomic_rates_mod, only: calc_colion_eq_fits
use mt19937_mod, only: genrand_real1
use m_mrgrnk, only: mrgrnk
implicit none
private

public :: particle_type
public :: particle_copy
public :: particle_transform
public :: particle_set_ye
public :: particle_set_ci_eq

public :: box_type
public :: box_adjust

public :: source_type

public :: particle_system_type
public :: particle_system_scale_comoving_to_physical
public :: particle_system_scale_physical_to_comoving
public :: particle_system_create_particle_random_access_list
public :: particle_system_create_particle_density_access_list
public :: particle_system_order_particles
public :: particle_system_mean_xHII_number_weight
public :: particle_system_mean_xHII_mass_weight
public :: particle_system_mean_xHII_volume_weight
public :: particle_system_set_ye
public :: particle_system_set_ci_eq
public :: particle_system_enforce_x_and_T_minmax
public :: particle_system_print_lun

public :: transformation_type
public :: return_bytes_per_particle
public :: return_bytes_per_source


!> Particle type. 
!=========================
type particle_type

   real(r4b)    :: pos(3)     !< x,y,z coordinates
   integer(i4b) :: id         !< particle id
   real(r4b)    :: mass       !< particle mass
   real(r4b)    :: T          !< temperature in K       
   real(r4b)    :: rho        !< density 
   real(r4b)    :: ye         !< electron fraction
   real(r4b)    :: xHI        !< nHI/nH 
   real(r4b)    :: xHII       !< nHII/nH
   real(r4b)    :: hsml       !< smoothing length
   integer(i8b) :: lasthit    !< indx of last ray to cross this particle

#ifdef incVel
   real(r4b)    :: vel(3)     !< x,y,z velocities
#endif

#ifdef incCloudy
   real(r4b)    :: xHI_cloudy !< cloudy eq solutions
#endif

#ifdef incHmf   
   real(r4b)    :: Hmf        !< Hydrogen mass fraction
#endif

#ifdef incHe
   real(r4b)    :: xHeI       !< nHeI/nHe 
   real(r4b)    :: xHeII      !< nHeII/nHe
   real(r4b)    :: xHeIII     !< nHeIII/nHe
#endif

#ifdef incHemf
   real(r4b)    :: Hemf       !< Helium mass fraction
#endif

#ifdef outGammaHI
   real(r4b)    :: gammaHI    !< time averaged HI photoionization rate
   real(r4b)    :: time       !< elapsed time in seconds - reset at outputs
#endif

#ifdef incEOS
   real(r4b)    :: eos        !< equation of state variable 
#endif

#ifdef incSFR
   real(r4b)    :: sfr        !< star formation rate
#endif

end type particle_type



!> source type
!================
type source_type

   integer(i4b) :: id        !< id added for AMUSE
   
   real(r4b)    :: pos(3)    !< x,y,z coordinates

#ifdef incVel
   real(r4b)    :: vel(3)    !< x,y,z velocities
#endif

   real(r4b)    :: L         !< luminosity
   real(r4b)    :: SpcType   !< spectral type
   integer(i4b) :: EmisPrf   !< emission profile
   real(r4b)    :: Lcdf      !< relates this luminosity to the other sources
   integer(i8b) :: lastemit  !< last ray emitted from this source
end type source_type



!> simulation box and boundary conditions
!==========================================
type box_type
   real(r8b)    :: tops(1:3)    !< upper x,y,z coordinates [code]
   real(r8b)    :: bots(1:3)    !< lower x,y,z coordinates [code]
   real(r8b)    :: lens(1:3)    !< side lengths [code] (xf-xi,yf-yi,zf-zi)
   real(r8b)    :: lens_cm(1:3) !< side lengths [cm] (xf-xi,yf-yi,zf-zi)
   real(r8b)    :: vol          !< box volume [code] (xlen*ylen*zlen)
   real(r8b)    :: vol_cm       !< box volume [cm^3] (xlen*ylen*zlen) 

   integer(i8b) :: bbound(1:3)  !< BCs for upper faces (0:vac 1:per -1:ref) 
   integer(i8b) :: tbound(1:3)  !< BCs for lower faces (0:vac 1:per -1:ref) 
end type box_type



!> particles, sources, and box
!=========================================
  type particle_system_type
     type(box_type) :: box                          !< the simulation box     
     type(particle_type), allocatable :: par(:)     !< all particles
     type(source_type), allocatable :: src(:)       !< all sources
     integer(i4b), allocatable :: acc_list(:)       !< access list
  end type particle_system_type


!> transformation type
!=======================
type transformation_type
   integer(i8b) :: fac(1:3)  !< newpos = fac * (oldpos - shift)
   real(i8b) :: shift(1:3)   !< newpos = fac * (oldpos - shift)
end type transformation_type





contains



!!============================================
!!
!!    BOX
!!
!!============================================

!> resets the box limits where the BCs are vacuum
!==================================================================
subroutine box_adjust(box,bot,top)
  type(box_type), intent(inout) :: box !< input box
  real(r4b), intent(in) :: bot(3)       !< new bottoms
  real(r4b), intent(in) :: top(3)       !< new tops
  where (box%bbound==0) box%bots = bot
  where (box%tbound==0) box%tops = top
end subroutine box_adjust




!!============================================
!!
!!    PARTICLE
!!
!!============================================


!> creates a copy of a particle
!============================================================
function particle_copy(this) result(copy)
  type(particle_type), intent(in) :: this  !< input particle
  type(particle_type) :: copy               !< particle copy
  copy = this
end function particle_copy


!> transforms a  particle 
!============================================================
subroutine particle_transform(this, transform)
  type(particle_type), intent(inout) :: this        !< input particle
  type(transformation_type), intent(in) :: transform !< transformation
  this%pos = transform%fac * (this%pos - transform%shift)
end subroutine particle_transform


!> set electron fraction, ye=ne/nH from ionization fractions
!==============================================================
subroutine particle_set_ye(par, dfltH_mf, dfltHe_mf, ne_bckgnd)

  type(particle_type) :: par
  real(r8b), intent(in) :: dfltH_mf
  real(r8b), intent(in) :: dfltHe_mf
  real(r8b), intent(in) :: ne_bckgnd
  integer(i8b) :: i
  real(r8b) :: Hmf
  real(r8b) :: Hemf
  real(r8b) :: nHe_over_nH
  
  par%ye = par%xHII + ne_bckgnd
  
#ifdef incHe
  
#ifdef incHmf
  Hmf = par%Hmf
#else
  Hmf = dfltH_mf
#endif
  
#ifdef incHemf
  Hemf = par%Hemf
#else
  Hemf = dfltHe_mf
#endif
  
  nHe_over_nH = 0.25d0 * Hemf / Hmf
  par%ye = par%ye + ( par%xHeII + 2.0d0 * par%xHeIII ) * nHe_over_nH
  
#endif
  
  
end subroutine particle_set_ye


!> sets ionization fractions to their collisional equilibrium values
!======================================================================
subroutine particle_set_ci_eq(par, caseA, DoH, DoHe, fit)

  type(particle_type) :: par
  logical, intent(in) :: caseA(2)   !< 1st slot for H, 2nd for He  
  logical, intent(in) :: DoH        !< set Hydrogen?
  logical, intent(in) :: DoHe       !< set Helium?
  character(*), intent(in) :: fit   !< one of ['hui','cen']
  real(r8b) :: T                    !< 8 byte temperature
  real(r8b) :: xvec(5)              !< [xHI,xHII,xHeI,xHeII,xHeIII]

  T = par%T
  call calc_colion_eq_fits(fit, T, caseA, xvec)

  if (DoH) then
     par%xHI = xvec(1)
     par%xHII = xvec(2)
  endif
  if (DoHe) then
#ifdef incHe
     par%xHeI = xvec(3)
     par%xHeII = xvec(4)
     par%xHeIII = xvec(5)
#else
     write(*,*) 
     write(*,*) 'In particle_set_ci_eq'
     write(*,*) 'DoHe = .true. but incHe macro not defined in Makefile'
     stop
#endif
  endif

end subroutine particle_set_ci_eq








!> figures out how many bytes of RAM are needed per particle
!========================================================================
function return_bytes_per_particle() result(bpp)
  integer(i4b)  :: bpp   !< bytes per particle

  bpp = 12       ! positions
  bpp = bpp + 4  ! ID
  bpp = bpp + 4  ! mass
  bpp = bpp + 4  ! temperature
  bpp = bpp + 4  ! rho
  bpp = bpp + 4  ! ye
  bpp = bpp + 8  ! H ionization fractions
  bpp = bpp + 4  ! hsml
  bpp = bpp + 8  ! last hit index

#ifdef incVel
  bpp = bpp + 12 ! velocities  
#endif

#ifdef incCloudy
  bpp = bpp + 4  ! cloudy table xHI
#endif

#ifdef incHmf
  bpp = bpp + 4  ! H mass fraction
#endif

#ifdef incHe
  bpp = bpp + 12 ! He ionization fractions
#endif

#ifdef incHemf
  bpp = bpp + 4  ! He mass fraction
#endif

#ifdef outGammaHI
  bpp = bpp + 4  ! GammaHI tracking
  bpp = bpp + 4  ! time var
#endif

#ifdef incEOS
  bpp = bpp + 4  ! Equation of State variable
#endif

#ifdef incSFR
  bpp = bpp + 4  ! Star formation rate
#endif

end function return_bytes_per_particle



!> figures out how many bytes of RAM are needed per source
!========================================================================
function return_bytes_per_source() result(bps)
  integer(i4b)  :: bps   !< bytes per source

  bps = 12       ! positions
  bps = bps + 4  ! luminosity
  bps = bps + 4  ! spectral type
  bps = bps + 4  ! emission profile
  bps = bps + 4  ! luminosity cumulative distribution function
  bps = bps + 8  ! last emit index

#ifdef incVel
  bps = bps + 12 ! velocities  
#endif

end function return_bytes_per_source








!!============================================
!!
!!    PARTICLE SYSTEM
!!
!!============================================

!> scales particles, sources, and the box from comoving to physical values.
! velocity is taken from Gadget code value to peculiar. 
!==========================================================================
subroutine particle_system_scale_comoving_to_physical(this, a, h)

  character(clen), parameter :: myname="scale_comoving_to_physical"
  integer, parameter :: verb=2
  character(clen) :: str,fmt

  type(particle_system_type) :: this
  real(r8b), intent(in) :: a    !< scale factor
  real(r8b), intent(in) :: h    !< hubble parameter (little h)

  call mywrite("   scaling comoving to physical coordinates", verb)
  fmt = "(A,F12.5,T22,A,T25,F12.5)"
  write(str,fmt) "   a = ", a, "h = ", h
  call mywrite(str,verb)

  ! particles
  !------------------------------------------------
  this%par%pos(1) = this%par%pos(1) * a / h
  this%par%pos(2) = this%par%pos(2) * a / h
  this%par%pos(3) = this%par%pos(3) * a / h

#ifdef incVel
  this%par%vel(1) = this%par%vel(1) * sqrt(a) 
  this%par%vel(2) = this%par%vel(2) * sqrt(a) 
  this%par%vel(3) = this%par%vel(3) * sqrt(a) 
#endif

  this%par%mass = this%par%mass / h
  this%par%hsml = this%par%hsml * a / h
  this%par%rho  = ( this%par%rho / (a*a*a) ) * (h*h)

  ! sources
  !------------------------------------------------
  this%src%pos(1) = this%src%pos(1) * a / h
  this%src%pos(2) = this%src%pos(2) * a / h
  this%src%pos(3) = this%src%pos(3) * a / h

#ifdef incVel
  this%src%vel(1) = this%src%vel(1) * sqrt(a) 
  this%src%vel(2) = this%src%vel(2) * sqrt(a)
  this%src%vel(3) = this%src%vel(3) * sqrt(a)
#endif

  ! box
  !------------------------------------------------
  this%box%tops = this%box%tops * a / h
  this%box%bots = this%box%bots * a / h
  
  this%box%lens    = this%box%lens * a / h
  this%box%lens_cm = this%box%lens_cm * a / h
  
  this%box%vol    = product( this%box%lens )
  this%box%vol_cm = product( this%box%lens_cm )
  

end subroutine particle_system_scale_comoving_to_physical


!> scales particles, sources, and the box from physical to comoving values.
! velocity is taken from peculiar to Gadget code value. 
!==========================================================================
subroutine particle_system_scale_physical_to_comoving(this, a, h)

  character(clen), parameter :: myname="scale_physical_to_comoving"
  integer, parameter :: verb=2
  character(clen) :: str,fmt

  type(particle_system_type) :: this
  real(r8b), intent(in) :: a   !< scale factor 
  real(r8b), intent(in) :: h   !< hubble parameter (little h)

  call mywrite("   scaling physical to comoving coordinates", verb)
  fmt = "(A,F12.5,T22,A,T25,F12.5)"
  write(str,fmt) "   a = ", a, "h = ", h
  call mywrite(str,verb)

  ! particles
  !------------------------------------------------
  this%par%pos(1) = this%par%pos(1) / a * h
  this%par%pos(2) = this%par%pos(2) / a * h
  this%par%pos(3) = this%par%pos(3) / a * h

#ifdef incVel
  this%par%vel(1) = this%par%vel(1) / sqrt(a)
  this%par%vel(2) = this%par%vel(2) / sqrt(a)
  this%par%vel(3) = this%par%vel(3) / sqrt(a)
#endif

  this%par%mass = this%par%mass * h
  this%par%hsml = this%par%hsml / a * h
  this%par%rho = this%par%rho * (a*a*a) / (h*h)

  ! sources
  !------------------------------------------------
  this%src%pos(1) = this%src%pos(1) / a * h
  this%src%pos(2) = this%src%pos(2) / a * h
  this%src%pos(3) = this%src%pos(3) / a * h

#ifdef incVel
  this%src%vel(1) = this%src%vel(1) / sqrt(a)
  this%src%vel(2) = this%src%vel(2) / sqrt(a)
  this%src%vel(3) = this%src%vel(3) / sqrt(a)
#endif

  ! box
  !------------------------------------------------
  this%box%tops = this%box%tops / a * h
  this%box%bots = this%box%bots / a * h
  
  this%box%lens    = this%box%lens / a * h
  this%box%lens_cm = this%box%lens_cm / a * h
  
  this%box%vol    = product( this%box%lens )
  this%box%vol_cm = product( this%box%lens_cm )


end subroutine particle_system_scale_physical_to_comoving



!> allows for accessing the particles in a random order
!------------------------------------------------------
subroutine particle_system_create_particle_random_access_list( psys )
  type(particle_system_type) :: psys

  integer(i4b) :: i
  integer(i4b) :: n
  real(r4b), allocatable :: randoms(:)

  n = size( psys%par )
  if (.not. allocated(psys%acc_list) ) allocate( psys%acc_list(n) )
  allocate( randoms(n) )

  do i = 1, size(randoms)
     randoms(i) = genrand_real1()
  end do
  
  call mrgrnk( randoms, psys%acc_list )

  deallocate( randoms )


end subroutine particle_system_create_particle_random_access_list



!> allows for accessing the particles from least to most dense
!--------------------------------------------------------------
subroutine particle_system_create_particle_density_access_list( psys )
  type(particle_system_type) :: psys

  integer(i4b) :: i
  integer(i4b) :: n
  real(r4b), allocatable :: rhos(:)

  n = size( psys%par )
  if (.not. allocated(psys%acc_list) ) allocate( psys%acc_list(n) )
  allocate( rhos(n) )

  do i = 1,n     
     rhos(i) = psys%par(i)%rho
  end do
  call mrgrnk( rhos, psys%acc_list )
  
  deallocate( rhos )

end subroutine particle_system_create_particle_density_access_list


! this routine rearranges the particles in the particle system so that 
! they are stored in the sequence given by the array order.  For example
! order = [3,1,2] takes the third particle to the first position, the first 
! particle to the second position, and the second particle to the third 
! position.  the array order is not preserved during the routine

!> reorders the particles according to the array order
!===========================================================================
subroutine particle_system_order_particles(this, order)
  type(particle_system_type), intent(inout) :: this !< input particle system
  integer(i4b), intent(inout) :: order(:)  !< desired order

  type(particle_type) :: par
  integer(i8b) :: i
  integer(i8b) :: goal
  integer(i8b) :: npar
  
  if (size(this%par) /= size(order)) stop "size(this%par) /= size(order)"
  npar = size(this%par)
  
  do i=1,npar 
     par=this%par(i)
     goal=order(i)
     do while(goal < i)
        goal=order(goal)
        order(i)=goal
     enddo
     this%par(i)=this%par(goal)
     this%par(goal)=par 
  enddo
  do i=1,npar
     order(i)=i
  enddo

end subroutine particle_system_order_particles


!> calculates the number weighted mean value of xHII
!========================================================================
function particle_system_mean_xHII_number_weight(this) result(numionfrac)
  type(particle_system_type), intent(in) :: this !< input particle system  
  real(r8b) :: numionfrac !< number weighted global ionization fraction
  integer :: i
     
  numionfrac = 0.0d0
  do i = 1,size(this%par)
     numionfrac = numionfrac + this%par(i)%xHII
  end do
  numionfrac = numionfrac / size(this%par)
  
end function particle_system_mean_xHII_number_weight
   
   
!> calculates the mass weighted mean value of xHII
!========================================================
function particle_system_mean_xHII_mass_weight(this, DfltH_mf) result(massionfrac)
  type(particle_system_type), intent(in) :: this !< input particle system  
  real(r8b) :: DfltH_mf     !< H_mf if we dont have a value for each par     
  real(r8b) :: massionfrac  !< ion fraction m weighted
  real(r8b) :: masstot      !< total volume
  real(r8b) :: Hmf          !< Hydrogen mass fraction
  integer :: i
  
  massionfrac = 0.0d0
  masstot = 0.0d0
  do i = 1,size(this%par)
#ifdef incHmf
     Hmf = this%par(i)%Hmf
#else
     Hmf = dfltH_mf
#endif
     massionfrac = massionfrac + this%par(i)%mass * Hmf * this%par(i)%xHII
     masstot = masstot + this%par(i)%mass * Hmf
  end do
  massionfrac = massionfrac / masstot
  
end function particle_system_mean_xHII_mass_weight

  
!> calculates the volume weighted mean of xHII
!========================================================
function particle_system_mean_xHII_volume_weight(this) result(volionfrac)
  type(particle_system_type), intent(in) :: this !< input particle system  
  real(r8b) :: volionfrac                    !< ion fraction v weighted
  real(r8b) :: voltot                        !< total volume
  real(r8b) :: h3                            !< hsml^3
  integer :: i
  
  volionfrac = 0.0d0
  voltot = 0.0d0
  do i = 1,size(this%par)
     h3 = this%par(i)%hsml * this%par(i)%hsml * this%par(i)%hsml
     volionfrac = volionfrac + h3 * this%par(i)%xHII
     voltot = voltot + h3
  end do
  volionfrac = volionfrac / voltot
  
end function particle_system_mean_xHII_volume_weight


!> set electron fraction, ye=ne/nH from ionization fractions
!==============================================================
subroutine particle_system_set_ye(psys, dfltH_mf, dfltHe_mf, ne_bckgnd)

  type(particle_system_type) :: psys
  real(r8b), intent(in) :: dfltH_mf
  real(r8b), intent(in) :: dfltHe_mf
  real(r8b), intent(in) :: ne_bckgnd
  integer(i8b) :: i

  do i = 1,size(psys%par)
     call particle_set_ye( psys%par(i), dfltH_mf, dfltHe_mf, ne_bckgnd)
  end do

end subroutine particle_system_set_ye


!> sets ionization fractions to their collisional equilibrium values
!======================================================================
subroutine particle_system_set_ci_eq(psys, caseA, DoH, DoHe, fit)

  type(particle_system_type) :: psys
  logical, intent(in) :: caseA(2)   !<  1st slot for H, 2nd for He  
  logical, intent(in) :: DoH        !<  set Hydrogen?
  logical, intent(in) :: DoHe       !<  set Helium?
  character(*), intent(in) :: fit   !< one of ['hui','cen']
  real(r8b) :: T                    !< 8 byte temperature
  real(r8b) :: xvec(5)              !< [xHI,xHII,xHeI,xHeII,xHeIII]
  integer(i8b) :: i
  
  do i = 1, size(psys%par)
     call particle_set_ci_eq( psys%par(i), caseA, DoH, DoHe, fit )
  end do

end subroutine particle_system_set_ci_eq



!> enforces a minimum and maximum value on the ionization fractions 
!! and temperatures
!=============================================================================
subroutine particle_system_enforce_x_and_T_minmax(psys,xmin,xmax,tmin,tmax)

  type(particle_system_type), intent(inout) :: psys !< inout particle system
  real(r8b), intent(in) :: xmin, xmax, tmin, tmax
  integer(i8b) :: i

  do i = 1,size(psys%par)
     if (psys%par(i)%xHI < xmin) psys%par(i)%xHI = xmin
     if (psys%par(i)%xHI > xmax) psys%par(i)%xHI = xmax
     if (psys%par(i)%xHII < xmin) psys%par(i)%xHII = xmin
     if (psys%par(i)%xHII > xmax) psys%par(i)%xHII = xmax
#ifdef incHe
     if (psys%par(i)%xHeI < xmin) psys%par(i)%xHeI = xmin
     if (psys%par(i)%xHeII < xmin) psys%par(i)%xHeII = xmin
     if (psys%par(i)%xHeIII < xmin) psys%par(i)%xHeIII = xmin
     if (psys%par(i)%xHeI > xmax) psys%par(i)%xHeI = xmax
     if (psys%par(i)%xHeII > xmax) psys%par(i)%xHeII = xmax
     if (psys%par(i)%xHeIII > xmax) psys%par(i)%xHeIII = xmax
#endif

     if (psys%par(i)%T < tmin) psys%par(i)%T = tmin
     if (psys%par(i)%T > tmax) psys%par(i)%T = tmax
  end do

end subroutine particle_system_enforce_x_and_T_minmax





!> outputs currently loaded particle data to the screen
!=================================================================
subroutine particle_system_print_lun(psys,str,lun)

  type(particle_system_type), intent(in) :: psys    !< particle system
  character(*), optional, intent(in) :: str          !< arbitrary string
  integer(i4b), optional, intent(in) :: lun          !< if present goes to file
  integer(i4b) :: outlun
  
  
  outlun=stdout
  if (present(lun)) outlun=lun


  99  format(72("-"))  
  100 format(A,T10,3ES15.5)
  101 format(A,T10,2I15,ES15.5)
  102 format(A,T10,2I15)
  103 format(A,T10,3I15)


  write(outlun,99) 
  if (present(str)) write(outlun,"(A)") trim(str)
  write(outlun,"(A,I15,A)") "particle data for ", size(psys%par), "  particles"

  write(outlun,*) 


  write(outlun,100) "xpos", minval(psys%par%pos(1)), &
       maxval(psys%par%pos(1)), meanval_real(psys%par%pos(1))

  write(outlun,100) "ypos", minval(psys%par%pos(2)), &
       maxval(psys%par%pos(2)), meanval_real(psys%par%pos(2))
  
  write(outlun,100) "zpos", minval(psys%par%pos(3)), &
       maxval(psys%par%pos(3)), meanval_real(psys%par%pos(3))   
    
  write(outlun,102) "id",   minval(psys%par%id), &
       maxval(psys%par%id)
  
  write(outlun,100) "mass", minval(psys%par%mass), &
       maxval(psys%par%mass), meanval_real(psys%par%mass)
  
  write(outlun,100) "T",    minval(psys%par(:)%T), &
       maxval(psys%par(:)%T), meanval_real(psys%par%T)
  
  write(outlun,100) "rho",  minval(psys%par%rho), &
       maxval(psys%par%rho), meanval_real(psys%par%rho)
  
  write(outlun,100) "ye",   minval(psys%par%ye), &
       maxval(psys%par%ye), meanval_real(psys%par%ye)
  
  write(outlun,100) "xHI",  minval(psys%par%xHI), &
       maxval(psys%par%xHI), meanval_real(psys%par%xHI)

  write(outlun,100) "xHII",  minval(psys%par%xHII), &
       maxval(psys%par%xHII), meanval_real(psys%par%xHII)

  write(outlun,100) "hsml", minval(psys%par%hsml), &
       maxval(psys%par%hsml), meanval_real(psys%par%hsml)

  write(outlun,101) "lasthit", minval(psys%par%lasthit), &
       maxval(psys%par%lasthit)

#ifdef incVel
  write(outlun,100) "xvel", minval(psys%par%vel(1)), &
       maxval(psys%par%vel(1)), meanval_real(psys%par%vel(1))
  
  write(outlun,100) "yvel", minval(psys%par%vel(2)), &
       maxval(psys%par%vel(2)), meanval_real(psys%par%vel(2))
  
  write(outlun,100) "zvel", minval(psys%par%vel(3)), &
       maxval(psys%par%vel(3)), meanval_real(psys%par%vel(3))
#endif

#ifdef incCloudy
  write(outlun,100) "xHI_cld", minval(psys%par%xHI_cloudy), &
       maxval(psys%par%xHI_cloudy), meanval_real(psys%par%xHI_cloudy)
#endif
  
#ifdef incHmf  
  write(outlun,100) "Hmf",    minval(psys%par%Hmf), &
       maxval(psys%par%Hmf), meanval_real(psys%par%Hmf)
#endif

#ifdef incHe
  write(outlun,100) "xHeI",   minval(psys%par%xHeI), &
       maxval(psys%par%xHeI), meanval_real(psys%par%xHeI)
  
  write(outlun,100) "xHeII",  minval(psys%par%xHeII), &
       maxval(psys%par%xHeII), meanval_real(psys%par%xHeII)
  
  write(outlun,100) "xHeIII", minval(psys%par%xHeIII), &
       maxval(psys%par%xHeIII), meanval_real(psys%par%xHeIII)
#endif

#ifdef incHemf  
  write(outlun,100) "Hemf",    minval(psys%par%Hemf), &
       maxval(psys%par%Hemf), meanval_real(psys%par%Hemf)
#endif

#ifdef outGammaHI
  write(outlun,100) "gammaHI",   minval(psys%par%gammaHI), &
       maxval(psys%par%gammaHI), meanval_real(psys%par%gammaHI)
  
  write(outlun,100) "time(s)",   minval(psys%par%time), &
       maxval(psys%par%time), meanval_real(psys%par%time)
#endif

#ifdef incEOS
  write(outlun,100) "eos",   minval(psys%par%eos), &
       maxval(psys%par%eos), meanval_real(psys%par%eos)  
#endif

#ifdef incSFR
  write(outlun,100) "sfr",   minval(psys%par%sfr), &
       maxval(psys%par%sfr), meanval_real(psys%par%sfr)  
#endif

  write(outlun,*)
  
  write(outlun,100) "Box Uppers = ", psys%box%tops
  write(outlun,100) "Box Lowers = ", psys%box%bots
  write(outlun,103) "Upr BCs    = ", psys%box%tbound
  write(outlun,103) "Lwr BCs    = ", psys%box%bbound
  
  write(outlun,99) 
  
end subroutine particle_system_print_lun


!> calculates the mean value w/o using the intrinsics
!=================================================================
function meanval_real(arr) result (mean)
  real(r4b), dimension(:), intent(in) :: arr  !< array to average
  real(r8b) :: mean                           !< mean value to return
  integer(i8b) :: i
  mean = 0.0d0
  do i = 1,size(arr)
     mean = mean + arr(i)
  end do
  mean = mean / size(arr)
end function meanval_real





end module particle_system_mod
