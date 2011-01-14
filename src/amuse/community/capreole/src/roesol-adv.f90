module hydrosolver
  
  ! Module for Capreole (f90)
  ! Author: Garrelt Mellema
  ! Date: 2004-05-11

  ! This module contains the Roe solver for solving the Euler equations
  ! in two dimensions on a curvilinear grid, including advected quantities.
  ! 
  ! The implementation is in the `fluctuations' form of LeVeque.
  ! It includes the transverse waves for improved 2D performance as an option. 
  ! The Euler quantities use `flux reconstruction', with stationary 
  ! extrapolation for the geometric source terms.

  ! The main routine is split into many subroutines for readability
  ! and to aid the compiler.

  ! This version has routines for constructing and destructing the arrays 
  ! needed for the solver, since otherwise they would be placed on the
  ! stack, which can cause problems.
  
  use precision, only: dp
  use sizes, only: nrOfDim,neq,neuler,mbc,RHO,EN
  use atomic, only: gamma1
  !use geometry

  implicit none

  private
  
  real(kind=dp),parameter,private    :: HALF=0.5d0 ! 1/2
  real(kind=dp),parameter,private    :: ONE=1.0d0  ! 1
  real(kind=dp),parameter,private    :: ZERO=0.0d0 ! 0

  ! Number of waves
  integer,parameter,private :: nwaves=5    ! number of waves
      
  ! Limiters
  integer,parameter,private :: NO_LIMITER=0
  integer,parameter,private :: MON_CEN=1
  integer,parameter,private :: VAN_LEER=2
  integer,parameter,private :: SUPERBEE=3
  integer,parameter,private :: USE_LIMITER=SUPERBEE

  ! Transverse wave directions
  integer,parameter,private :: RIGHT=0
  integer,parameter,private :: LEFT=-1

  ! The following variables are public so that they can be filled
  ! or used outside of the solver. They correspond to the one
  ! spatial dimension
  real(kind=dp),dimension(:,:),allocatable,public :: state1d
  ! volume times pressure
  real(kind=dp),dimension(:),allocatable,public :: wp 
  ! dstate will contain the state changes:
  real(kind=dp),dimension(:,:),allocatable,public :: dstate
  !$OMP THREADPRIVATE(state1d,wp,dstate)  

  ! The following variables are private since they are only
  ! needed inside the solver.
  real(kind=dp),dimension(:,:),allocatable :: w ! parameter vector
  
  real(kind=dp),dimension(:,:),allocatable :: fluxc ! central flux at i
  real(kind=dp),dimension(:,:),allocatable :: fluxl ! left flux at i-1/2
  real(kind=dp),dimension(:,:),allocatable :: fluxr ! right flux at i-1/2
  real(kind=dp),dimension(:,:),allocatable :: fludif ! flux difference at i-1/2
  real(kind=dp),dimension(:,:),allocatable :: stadif ! state difference at i-1/2
  !$OMP THREADPRIVATE(w,fluxc,fluxl,fluxr,fludif,stadif)

  real(kind=dp),dimension(:),allocatable :: rrt  ! 1/rho average
  real(kind=dp),dimension(:,:),allocatable :: vt ! roe average velocity 2
  real(kind=dp),dimension(:),allocatable :: ht ! roe average enthalpy
  real(kind=dp),dimension(:),allocatable :: absvt  ! absolute velocity
  real(kind=dp),dimension(:),allocatable :: vst2    ! sound speed squared
  real(kind=dp),dimension(:),allocatable :: vst   ! sound speed
  !$OMP THREADPRIVATE(rrt,vt,ht,absvt,vst2,vst)
  
  real(kind=dp),dimension(:,:),allocatable :: eiglam ! eigen values/characteristics
  real(kind=dp),dimension(:,:),allocatable :: sgn ! sign of eigenvalues
  
  real(kind=dp),dimension(:,:,:),allocatable :: eigv
  
  real(kind=dp),dimension(:),allocatable :: uvdif ! inner product of v and 
  real(kind=dp),dimension(:,:),allocatable :: a ! projection coefficient
  
  real(kind=dp),dimension(:,:,:),allocatable :: wave
  real(kind=dp),dimension(:,:,:),allocatable :: adq
  real(kind=dp),dimension(:,:),allocatable :: flux  ! 2nd order correction
  !$OMP THREADPRIVATE(eiglam,sgn,eigv,uvdif,a,wave,adq,flux)

  public :: constr_solver, destr_solver, solver
  
contains

  !============================================================================
  subroutine constr_solver (mesh)
    
    integer,intent(in) :: mesh

    ! Initializes solver variables

    allocate(state1d(1-mbc:mesh+mbc,neq))
    allocate(wp(1-mbc:mesh+mbc))
    allocate(dstate(1-mbc:mesh+mbc,neq))

    allocate(w(1-mbc:mesh+mbc,neuler))
    allocate(fluxc(1-mbc:mesh+mbc,neuler))
    allocate(fluxl(2-mbc:mesh+mbc,neuler))
    allocate(fluxr(2-mbc:mesh+mbc,neuler))
    allocate(fludif(2-mbc:mesh+mbc,neuler))
    allocate(stadif(2-mbc:mesh+mbc,neq))

    allocate(rrt(2-mbc:mesh+mbc))
    allocate(vt(2-mbc:mesh+mbc,nrOfDim))
    allocate(ht(2-mbc:mesh+mbc))
    allocate(absvt(2-mbc:mesh+mbc))
    allocate(vst2(2-mbc:mesh+mbc))
    allocate(vst(2-mbc:mesh+mbc))
    
    allocate(eiglam(2-mbc:mesh+mbc,nwaves))
    allocate(sgn(2-mbc:mesh+mbc,nwaves))
    
    allocate(eigv(2-mbc:mesh+mbc,neq,nwaves))
    
    allocate(uvdif(2-mbc:mesh+mbc))
    allocate(a(2-mbc:mesh+mbc,nwaves))
    
    allocate(wave(2-mbc:mesh+mbc,neuler,nwaves))
    allocate(adq(2-mbc:mesh+mbc,neq,LEFT:RIGHT))
    allocate(flux(2-mbc:mesh+mbc,neq))
    
  end subroutine constr_solver

!-----------------------------------------------------------------------------

  subroutine destr_solver

    ! Destructs solver variables

    deallocate(state1d)
    deallocate(wp)
    deallocate(dstate)

    deallocate(w)
    deallocate(fluxc)
    deallocate(fluxl)
    deallocate(fluxr)
    deallocate(fludif)
    deallocate(stadif)

    deallocate(rrt)
    deallocate(vt)
    deallocate(ht)
    deallocate(absvt)
    deallocate(vst2)
    deallocate(vst)
    
    deallocate(eiglam)
    deallocate(sgn)
    
    deallocate(eigv)
    
    deallocate(uvdif)
    deallocate(a)
    
    deallocate(wave)
    deallocate(adq)
    deallocate(flux)
    
  end subroutine destr_solver

!-----------------------------------------------------------------------------

  subroutine solver (mesh,dt,dx,dy,dz,V1,V2,V3,ij,ik,ierror)
    
    ! The central Roe solver routine, supplies dstate (normal and transverse
    ! state changes)
    
    ! parameters
    ! small    - a small number
    ! nwaves   - number of waves
    
    integer,intent(in) :: mesh
    real(kind=dp),intent(in) :: dt,dx,dy,dz
    integer,intent(in) :: ij,ik
    integer,intent(in) :: V1,V2,V3 ! indices of the the two velocities:
                                        ! V1: integration direction
                                        ! V2,V3: perpendicular direction
    integer,intent(out) :: ierror ! control integer

    integer :: i,ieq,mw
    real(kind=dp) :: dtdx!,dtdy

    ! smooth transition through zero characteristics
    !real(kind=dp) :: fctr
    ! steepness of transition function
    !real(kind=dp),dimension(nwaves) :: beta= (/ 1d4, 1d4, 1d4, 1d4 /)

    !real(kind=dp) :: rhol,rhor,ul,ur
    !real(kind=dp),dimension(2-mbc:mesh+mbc)   :: utadv

    integer :: ipres_error
    
    !------------------------------------------------------------------------

    ierror=0 ! initialise the control variable to 0
    ipres_error=0

    call initialize ()

    ! ---------------------------------------------------------------------
    ! Part 1: Roe solver
    ! ---------------------------------------------------------------------

    call calculate_parametervector ()

    call calculate_fluxes ()

    call stationary_extrapolation ()

    call calculate_differences ()

    call calculate_roe_averages (ierror)

    ! Only continue if we have a sound speed
    error: if (ierror == 0) then

       call calculate_jacobian (V1-1,V2-1,V3-1)
       
       call calculate_projection_coeffs (V1-1,V2-1,V3-1,stadif)
       
       call split_waves (stadif,adq)
       
       call flux_limiting (USE_LIMITER)
       
       call calculate_state_change ()
    
       ipres_error=pressure_fix ()
       
       ! Find final change in state
       dstate(1:mesh,1:neq)=dtdx*dstate(1:mesh,1:neq)
       
       
    endif error
    
    !--------------------------------------------------------------------------

  contains
    
    !--------------------------------------------------------------------------
    subroutine initialize
      
      ! Note: ifc 7.1 has a bug and does not map the flux
      ! adq, amq, and dstate arrays properly to this internal
      ! subroutine. Result: run time error "Adress Error"
      ! Fix: 
      ! 1) compile with -auto (this triggers other errors)
      ! 2) compile with -g

      ! Initialize some arrays to zero
      
      dstate=ZERO
      adq=ZERO
      flux=ZERO

    end subroutine initialize

    !--------------------------------------------------------------------------

    subroutine calculate_parametervector 

      ! ----------------------------------------------------------------------
      ! find the parameter-vector w from the state. 
      ! The state should already contain the volume terms.
      !
      ! 1 <-> sqrt(vol)*sqrt(rho)
      ! 2 <-> sqrt(vol)*sqrt(rho)*u
      ! 3 <-> sqrt(vol)*sqrt(rho)*v
      ! 4 <-> sqrt(vol)*sqrt(rho)*enthalpy
      ! wp<-> vol*pressure (is imported)
      ! ----------------------------------------------------------------------
      w(:,RHO)=sqrt(state1d(:,RHO))
      w(:,V1)=state1d(:,V1)/w(:,RHO)
      w(:,V2)=state1d(:,V2)/w(:,RHO)
      w(:,V3)=state1d(:,V3)/w(:,RHO)
      w(:,EN)=(state1d(:,EN)+wp(:))/w(:,RHO)
      
    end subroutine calculate_parametervector

    !--------------------------------------------------------------------------

    subroutine calculate_fluxes

      ! ---------------------------------------------------------------------
      ! calculate the fluxes at the cell centre
      ! ---------------------------------------------------------------------
      do i=1-mbc,mesh+mbc
         fluxc(i,RHO)=w(i,V1)*w(i,RHO)
         fluxc(i,V1)=w(i,V1)*w(i,V1) + wp(i)
         fluxc(i,V2)=w(i,V1)*w(i,V2)
         fluxc(i,V3)=w(i,V1)*w(i,V3)
         fluxc(i,EN)=w(i,V1)*w(i,EN)
      enddo
      
    end subroutine calculate_fluxes

    !--------------------------------------------------------------------------

    subroutine stationary_extrapolation ()

      ! ---------------------------------------------------------------------
      ! calculate the fluxes at the cell walls, using 
      ! `stationary extrapolation'
      ! ---------------------------------------------------------------------
      do ieq=1,neuler
         do i=2-mbc,mesh+mbc
            fluxl(i,ieq)=fluxc(i-1,ieq)
            fluxr(i,ieq)=fluxc(i,ieq)
         enddo
      enddo

    end subroutine stationary_extrapolation
    
    !--------------------------------------------------------------------------

    subroutine calculate_differences ()

      ! calculate the state differences at the cell walls.
      ! index i corresponds to i-1/2 (left cell well of cell i).
      do ieq=1,neq
         do i=2-mbc,mesh+mbc
            stadif(i,ieq)=state1d(i,ieq)-state1d(i-1,ieq)
         enddo
      enddo
      ! ---------------------------------------------------------------------
      ! calculate the flux differences at the cell walls
      ! ---------------------------------------------------------------------
      !do ieq=1,neuler
      !   do i=2-mbc,mesh+mbc
      !      fludif(i,ieq)=fluxr(i,ieq)-fluxl(i,ieq)
      !   enddo
      !enddo
      
      ! ---------------------------------------------------------------------
      ! calculate the state differences at the cell walls
      ! (for the advected quantities)
      ! ---------------------------------------------------------------------
      !do ieq=neuler+1,neq
      !   do i=2-mbc,mesh+mbc
      !      stadif(i,ieq)=state1d(i,ieq)-state1d(i-1,ieq)
      !   enddo
      !enddo

    end subroutine calculate_differences

    !--------------------------------------------------------------------------

    subroutine calculate_roe_averages (ierror)
      
      integer,intent(out) :: ierror

      ! calculate Roe averages
      do i=2-mbc,mesh+mbc
         rrt(i)=ONE/(w(i-1,RHO)+w(i,RHO))
         vt(i,V1-1)=(w(i-1,V1)+w(i,V1))*rrt(i) ! velocity
         vt(i,V2-1)=(w(i-1,V2)+w(i,V2))*rrt(i) ! velocity
         vt(i,V3-1)=(w(i-1,V3)+w(i,V3))*rrt(i) ! velocity
         ht(i)=(w(i-1,EN)+w(i,EN))*rrt(i)       ! enthalpy
      enddo
      
      ! absolute velocity
      absvt(:)=HALF*(vt(:,1)*vt(:,1)+vt(:,2)*vt(:,2)+vt(:,3)*vt(:,3))
      
      ! sound speed
      vst2(:)=gamma1*(ht(:)-absvt(:))
      
      do i=2-mbc,mesh+mbc
         ! check for negative values of vst2
         if (vst2(i) > 0.0) then         
            vst(i)=sqrt(vst2(i))
         else
            vst(i)=sqrt(abs(vst2(i)))
            ierror=ierror+1
            write(30,*) 'VST2 neagtive: ',i,ij,ik,vst2(i),V1,V2,V3
            write(30,*) 'velocity: ',vt(i-1:i,1:3)
            write(30,*) 'pressure: ',wp(i-1:i)
            !stop                       ! negative vst2 is fatal!
         endif
      enddo
      
    end subroutine calculate_roe_averages
    
    !--------------------------------------------------------------------------

    subroutine calculate_jacobian (itr,ip1,ip2)
      
      ! Calculates the eigenvalues and eigenvectors of the Jacobian
      ! for the coordinate direction itr.
      ! The results are eiglam and eigv
      
      ! itr  - index of transverse direction we are considering
      ! ip1  - index of first parallel direction
      ! ip2  - index of second parallel direction

      integer :: itr,ip1,ip2
      
      ! calculate the eigenvalues and vectors of the jacobian
      
      eiglam(:,1)=vt(:,itr)-vst(:)
      eiglam(:,2)=vt(:,itr)
      eiglam(:,3)=vt(:,itr)+vst(:)
      eiglam(:,4)=vt(:,itr)
      eiglam(:,5)=vt(:,itr)
      
      eigv(:,RHO,1:3)=ONE    ! density
      eigv(:,RHO,4:5)=ZERO
      
      ! shift velocity indices to state indices
      eigv(:,itr+1,1)=vt(:,itr)-vst(:)
      eigv(:,itr+1,2)=vt(:,itr)
      eigv(:,itr+1,3)=vt(:,itr)+vst(:)
      eigv(:,itr+1,4:5)=ZERO
      
      eigv(:,ip1+1,1)=vt(:,ip1)              
      eigv(:,ip1+1,2)=vt(:,ip1)
      eigv(:,ip1+1,3)=vt(:,ip1)
      eigv(:,ip1+1,4)=ONE
      eigv(:,ip1+1,5)=ZERO
      
      eigv(:,ip2+1,1)=vt(:,ip2)              
      eigv(:,ip2+1,2)=vt(:,ip2)              
      eigv(:,ip2+1,3)=vt(:,ip2)              
      eigv(:,ip2+1,4)=ZERO
      eigv(:,ip2+1,5)=ONE
      
      eigv(:,EN,1)=ht(:)-vt(:,itr)*vst(:) ! energy
      eigv(:,EN,2)=absvt(:)
      eigv(:,EN,3)=ht(:)+vt(:,itr)*vst(:)
      eigv(:,EN,4)=vt(:,ip1)
      eigv(:,EN,5)=vt(:,ip2)
      
    end subroutine calculate_jacobian
    
    !-------------------------------------------------------------------------

    subroutine calculate_projection_coeffs (itr,ip1,ip2,q)
      
      ! Projects the quantity q onto the eigenvectors of the Jacobian
      
      ! itr  - index of transverse direction we are considering (x=1, y=2, z=3)
      ! ip1   - index of first parallel direction
      ! ip2   - index of second parallel direction
      ! q     - quantity to project
      
      integer :: itr,ip1,ip2
      real(kind=dp),dimension(2-mbc:mesh+mbc,neq) :: q

      real(kind=dp),dimension(2-mbc:mesh+mbc) ::  vdotq
      
      ! calculate the projection of q onto the eigenvectors
      
      vdotq(:)=vt(:,itr)*q(:,itr+1)+vt(:,ip1)*q(:,ip1+1)+vt(:,ip2)*q(:,ip2+1)
      
      do i=2-mbc,mesh+mbc
         a(i,1)=HALF*(gamma1*(absvt(i)*q(i,RHO)+q(i,EN)-vdotq(i))- &
              vst(i)*(q(i,itr+1)-vt(i,itr)*q(i,RHO)))/(vst2(i))
         a(i,2)=gamma1*((ht(i)-2.0d0*absvt(i))*q(i,RHO)+ & 
              vdotq(i)-q(i,EN))/(vst(i)*vst(i))
         a(i,3)=HALF*(gamma1*(absvt(i)*q(i,RHO)+q(i,EN)-vdotq(i))+ &
              vst(i)*(q(i,itr+1)-vt(i,itr)*q(i,RHO)))/(vst2(i))
         a(i,4)=q(i,ip1+1)-vt(i,ip1)*q(i,RHO)
         a(i,5)=q(i,ip2+1)-vt(i,ip2)*q(i,RHO)
      enddo

    end subroutine calculate_projection_coeffs

!------------------------------------------------------------------------------

    subroutine split_waves (q,splitq)

      ! ---------------------------------------------------------------------
      ! Calculate the waves and divide them into left and right going
      ! (first order solution).
      !
      ! There may be cases where a smoother transition than this step
      ! function is better, this avoids problems with expansion points
      ! (entropy fix) and perhaps odd-even decoupling. The commented
      ! out lines do a transition using an arctan function.
      ! ---------------------------------------------------------------------

      real(kind=dp),dimension(2-mbc:mesh+mbc,neq) :: q
      real(kind=dp),dimension(2-mbc:mesh+mbc,neq,LEFT:RIGHT) :: splitq

      ! put splitq to zero
      splitq(:,:,:)=ZERO

      ! For the Euler quantities
      do ieq=1,neuler
         do i=2-mbc,mesh+mbc
            do mw=1,nwaves
               wave(i,ieq,mw)=a(i,mw)*eigv(i,ieq,mw)
               splitq(i,ieq,LEFT)=splitq(i,ieq,LEFT)+ &
                    min(eiglam(i,mw),ZERO)*wave(i,ieq,mw)
               splitq(i,ieq,RIGHT)=splitq(i,ieq,RIGHT)+ &
                    max(eiglam(i,mw),ZERO)*wave(i,ieq,mw)
            enddo
         enddo
      enddo

      ! For the advected quantities
      do ieq=neuler+1,neq
         do i=2-mbc,mesh+mbc
            splitq(i,ieq,LEFT)=splitq(i,ieq,LEFT) + &
                 min(eiglam(i,4),ZERO)*q(i,ieq)
            splitq(i,ieq,RIGHT)=splitq(i,ieq,RIGHT) + &
                 max(eiglam(i,4),ZERO)*q(i,ieq)
         enddo
      enddo

    end subroutine split_waves
      
    
    !=========================================================================

    subroutine flux_limiting (limfunc)
      
      ! Calculates the second order fluxes, using a limiter
      
      real(kind=dp),parameter ::  small=1.0d-24 ! small number
      
      integer,intent(in) :: limfunc

      integer :: isb
      
      real(kind=dp) :: wnorm2,wip,r,philim
      
      ! ---------------------------------------------------------------------
      ! apply the fluxcorrection to make 2nd order correction fluxes
      !
      ! This comes from LeVeque: take the inner product of the wave and
      ! the upwind wave and compare this with the absolute value of the
      ! wave. Since they are vectors, this is a better way to compare
      ! them.
      ! ---------------------------------------------------------------------
      dtdx=dt/dx
      ! Set flux to zero
      flux(:,:)=ZERO
      
      ! For each wave, find the limiter using the ratio of the
      ! inner product of the waves
      do mw=1,nwaves
         do i=1,mesh+1

            ! Find out whether to look left or right (for individual waves)
            if (eiglam(i,mw).lt.ZERO) then
               isb=i+1
            else
               isb=i-1
            endif

            ! calculate inner products (more accurate)
            wip=dot_product(wave(isb,1:neuler,mw),wave(i,1:neuler,mw))
            wnorm2=dot_product(wave(i,1:neuler,mw),wave(i,1:neuler,mw))
            
            ! calculate the ratio
            if (wnorm2.lt.small) then
               r=ZERO
            else
               r=wip/wnorm2
            endif
            
            ! Limiter
            philim = limiter(r,limfunc)
            
            ! set the 2nd order correction term flux
            do ieq=1,neuler   
               flux(i,ieq)=flux(i,ieq)+ &
                    abs(eiglam(i,mw))*(ONE-dtdx*abs(eiglam(i,mw)))* &
                    philim*wave(i,ieq,mw)
            enddo
         enddo
      enddo
      
      ! Limit the advective fluxes:
      !   eiglam is the velocity (eiglam(4))
      !   for accuracy take the minimum of all philim values
      
      do i=1,mesh+1
         philim=2.0d0

         ! Find out whether to look left or right (for individual waves)
         if (eiglam(i,4).lt.ZERO) then
            isb=i+1
         else
            isb=i-1
         endif

         ! calculate the ratio
         do ieq=neuler+1,neq 
            if (abs(stadif(i,ieq)).lt.1e-25) then
               r=2.0d0 
            else
               r=stadif(isb,ieq)/stadif(i,ieq)
            endif
            
            ! Limiter
            philim=min(philim,limiter(r,limfunc))
         enddo
         
         ! set the 2nd order correction term flux
         do ieq=neuler+1,neq  
            flux(i,ieq)=abs(eiglam(i,4))* &
                 (ONE-dtdx*abs(eiglam(i,4)))*philim*stadif(i,ieq)
         enddo
      enddo
      
    end subroutine flux_limiting
    
    !--------------------------------------------------------------------------
    
    function limiter(ratio,limfunc)
      
      ! Limiter function used in flux limiting.
      ! It limits a ratio using the function limfunc.
      
      real(kind=dp) :: limiter
      
      ! superbee parameter (1 to 2)
      real(kind=dp),parameter :: sbpar=1.8
      
      real(kind=dp),intent(in) :: ratio
      integer,intent(in) :: limfunc
      
      real(kind=dp) :: c
      
      select case (limfunc)
      case (NO_LIMITER)
         limiter = ZERO
      case(MON_CEN)
         ! monotinized centered 
         c = (1.d0 + ratio)/2.d0
         limiter = max(0.d0, min(c, 2.d0, 2.d0*ratio))
      case(VAN_LEER)
         ! van Leer
         limiter = (ratio + abs(ratio)) / (1.d0 + abs(ratio))
      case(SUPERBEE)
         ! Superbee (sbpar)
         limiter = max(ZERO,min(ONE,sbpar*ratio),min(sbpar,ratio))
         !!philim = min(philim,max(ZERO, min(ONE, sbpar*r), min(sbpar, r)))
      end select
      
    end function limiter

    !--------------------------------------------------------------------------

    subroutine calculate_state_change
      
      ! calculate the change in state (1st and 2nd order added)
      do ieq=1,neq
         do i=1,mesh
            dstate(i,ieq)=-adq(i,ieq,RIGHT)-adq(i+1,ieq,LEFT) &
                 +HALF*(flux(i,ieq)-flux(i+1,ieq))
         enddo
      enddo
      
    end subroutine calculate_state_change
    
    !--------------------------------------------------------------------------

    function pressure_fix ()
      
      ! ---------------------------------------------------------------------
      ! This routines checks whether the new solution will have positive
      ! pressure, density and energy density. If not two fixes are tried:
      ! 1) drop all the second order terms (flux)
      ! 2) add some diffusion
      !
      ! There is an issue with problems spreading after having fixed a
      ! point. It is necessary to search again for bad points in the
      ! neighbourhood after having fixed one point.
      ! ---------------------------------------------------------------------

      integer :: pressure_fix

      ! Coefficient of the diffusion
      real(kind=dp),parameter :: eta=0.0001

      real(kind=dp),dimension(mesh) :: ptest      ! negative pressure tests
      real(kind=dp),dimension(mesh) :: ptest2
      real(kind=dp),dimension(mesh) :: ptest3
      real(kind=dp),dimension(1-mbc:mesh+mbc,neq) :: diff

      integer :: imin

      pressure_fix=0
      diff(:,:)=ZERO

      ! ---------------------------------------------------------------------
      ! calculate test variable for negative pressure check
      ! Note: this check is not full proof since it does not incorporate
      ! the transverse fluxes
      ! ---------------------------------------------------------------------
      do i=1,mesh
         ptest(i)=prestest(i)
      enddo
      
      ! ---------------------------------------------------------------------
      ! check for negative pressure/internal energy/density and set change
      ! in state to first order if detected
      ! ---------------------------------------------------------------------
      p1loop: do i=1,mesh
         prtest1: if (ptest(i) <= ZERO .or.                         &
              (dx*state1D(i,EN)+dt*dstate(i,EN)) <= ZERO .or.      &
              (dx*state1D(i,RHO)+dt*dstate(i,RHO)) <= ZERO) then

            !!write(30,*) 'Roe solver: ptest1 failed: ', i,ij,V1,V2
            
            ! Set 2nd order fluxes to zero at the interfaces of i
            ! Recalculate dstate for i-1,i,i+1
            do ieq=1,neq
               flux(i,ieq)=ZERO
               flux(i+1,ieq)=ZERO
               dstate(i-1,ieq)=-adq(i-1,ieq,RIGHT)-adq(i,ieq,LEFT)+ &
                    HALF*flux(i-1,ieq)
               dstate(i,ieq)=-adq(i,ieq,RIGHT)-adq(i+1,ieq,LEFT)
               dstate(i+1,ieq)=-adq(i+1,ieq,RIGHT)-adq(i+2,ieq,LEFT)- &
                    HALF*flux(i+2,ieq)
            enddo
            
            ! Recalculate ptest for new dstate of i+1
            ! Do not do this if we are at the boundary
            if (i < mesh) ptest(i+1)=prestest(i+1)
            
            ! Recalculate ptest for new dstate of i-1 and below and fix 
            ! if needed.
            ! Do not do this if the cell was already flagged
            ! and if we are at the boundary
            imin=i-1
            backward1: do
               okbefore1: if (imin >= 1 .and. ptest(imin) > zero .and.        &
                    (dx*state1D(imin,EN)+dt*dstate(imin,EN)) > zero .and.     &
                    (dx*state1D(imin,RHO)+dt*dstate(imin,RHO)) > zero) then
                  ptest(imin)=prestest(imin)
                  badnow1: if (ptest(imin) <= ZERO .or.                       &
                       (dx*state1D(imin,EN)+dt*dstate(imin,EN)) <= ZERO .or.  &
                       (dx*state1D(imin,RHO)+dt*dstate(imin,RHO)) <= ZERO) then
                     !write(30,*) 'Roe solver: ptest1 back correction: ', & 
                     !     imin,ij,V1,V2
                     do ieq=1,neq
                        flux(imin,ieq)=0.0
                        dstate(imin-1,ieq)=-adq(imin-1,ieq,RIGHT)- &
                             adq(imin,ieq,LEFT) + HALF*flux(imin-1,ieq)
                        dstate(imin,ieq)=-adq(imin-1,ieq,RIGHT)- &
                             adq(imin,ieq,LEFT)
                     enddo
                  else
                     exit
                  endif badnow1
               else
                  exit
               endif okbefore1
               imin=imin-1
            end do backward1
            
         end if prtest1
      end do p1loop
      
      ! Check if it helped, set control variable if not
      do i=1,mesh
         ptest2(i)=prestest(i)
      enddo

      ! ---------------------------------------------------------------------
      ! check for negative pressure/internal energy/density and add diffusion
      ! if detected
      ! ---------------------------------------------------------------------
      p2loop: do i=1,mesh
         prtest2: if (ptest2(i) <= ZERO .or.                               &
              (dx*state1D(i,EN)+dt*dstate(i,EN)) <= ZERO .or.      &
              (dx*state1D(i,RHO)+dt*dstate(i,RHO)) <= ZERO) then
            !!write(30,*) 'Roe solver: ptest2 failed: ', i,ij,V1,V2
            do ieq=1,neq
               diff(i,ieq)=eta*(state1D(i-1,ieq)-state1D(i,ieq))/dtdx
               diff(i+1,ieq)=eta*(state1D(i,ieq)-state1D(i+1,ieq))/dtdx
               dstate(i-1,ieq)=-adq(i-1,ieq,RIGHT)-adq(i,ieq,LEFT)+ &
                    HALF*(flux(i-1,ieq)-flux(i,ieq))+        &
                    diff(i-1,ieq)-diff(i,ieq)
               dstate(i,ieq)=-adq(i,ieq,RIGHT)-adq(i+1,ieq,LEFT)+ &
                    HALF*(flux(i,ieq)-flux(i+1,ieq))+      &
                    diff(i,ieq)-diff(i+1,ieq)
               dstate(i+1,ieq)=-adq(i+1,ieq,RIGHT)-adq(i+2,ieq,LEFT)+ &
                    HALF*(flux(i+1,ieq)-flux(i+2,ieq))+        &
                    diff(i+1,ieq)-diff(i+2,ieq)
            enddo

            ! Recalculate ptest for new dstate of i+1
            ! Do not do this if we are at the boundary
            if (i < mesh) ptest2(i+1)=prestest(i+1)
            
            ! Recalculate ptest for new dstate of i-1 and below and fix 
            ! if needed.
            ! Do not do this if the cell was already flagged
            ! and if we are at the boundary
            imin=i-1
            backward2: do
               okbefore2: if (imin >= 1 .and. ptest2(imin) > zero .and.       &
                    (dx*state1D(imin,EN)+dt*dstate(imin,EN)) > zero .and.     &
                    (dx*state1D(imin,RHO)+dt*dstate(imin,RHO)) > zero) then
                  ptest2(imin)=prestest(imin)
                  badnow2: if (ptest2(imin) <= ZERO .or.                     &
                       (dx*state1D(imin,EN)+dt*dstate(imin,EN)) <= ZERO .or. &
                       (dx*state1D(imin,RHO)+dt*dstate(imin,RHO)) <= ZERO) then
                     !write(30,*) 'Roe solver: ptest2 back correction: ', & 
                     !     imin,ij,V1,V2
                     do ieq=1,neq
                        diff(imin,ieq)=eta*(state1D(imin-1,ieq)-       &
                             state1D(imin,ieq))/dtdx
                        dstate(imin-1,ieq)=-adq(imin-1,ieq,RIGHT)- &
                             adq(imin,ieq,LEFT)-                     &
                             diff(imin,ieq) +                        &
                             HALF*(flux(imin-1,ieq)-flux(imin,ieq))
                        dstate(imin,ieq)=-adq(imin,ieq,RIGHT)-     &
                             adq(imin+1,ieq,LEFT)+                   &
                             diff(imin,ieq)+                         &
                             HALF*(flux(imin,ieq)-flux(imin+1,ieq))
                     enddo
                  else
                     exit
                  endif badnow2
               else
                  exit
               endif okbefore2
               imin=imin-1
            enddo backward2
         endif prtest2
      enddo p2loop
      
      ! Check if it helped, set control variable if not
      do i=1,mesh
         ptest3(i)=prestest(i)
      enddo
      do i=1,mesh
         if (ptest3(i) <= ZERO .or.                                &
              (dx*state1D(i,EN)+dt*dstate(i,EN)) <= ZERO .or.      &
              (dx*state1D(i,RHO)+dt*dstate(i,RHO)) <= ZERO) then
            !write(30,*) 'Roe solver: ptest3 failed: ', i,ij,V1,V2
            pressure_fix=pressure_fix+1
         endif
      enddo
      
    end function pressure_fix

    !--------------------------------------------------------------------------

    function prestest(i)
      
      ! Pressure test function
      ! Calculates the new pressure times a constant, to allow
      ! checking the sign of the new pressure in pressure_fix.
      
      real(kind=dp) :: prestest
      
      integer,intent(in) :: i
      
      prestest=2.0d0*(dx*state1D(i,EN)+dt*dstate(i,EN))*                  &
            (dx*state1D(i,RHO)+dt*dstate(i,RHO))-                         &
            (dx*state1D(i,V1)+dt*dstate(i,V1))**2-  &
            (dx*state1D(i,V2)+dt*dstate(i,V2))**2-  &
            (dx*state1D(i,V3)+dt*dstate(i,V3))**2

    end function prestest
    
  end subroutine solver
    
end module hydrosolver





