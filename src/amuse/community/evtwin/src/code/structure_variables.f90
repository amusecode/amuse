!> 
!! The variables in this module used to be in the common block vbles
!< 

module structure_variables
   use real_kind
   use mesh
   
   implicit none
   !These variables are for the current mesh point (?)
   real(double) :: LoLedd        ! Eddington factor: luminosity over Eddington luminosity
   real(double) :: dg            ! Convective grad_r - grad_a  for Schwarzschild criterion
   real(double) :: eg            ! Convective grad_r - grad_a  for overshooting
   real(double) :: grad          ! grad T (used to be called gradt in some routines, change back?)
   real(double) :: eth           ! Thermal energy-generation rate(?)
   real(double) :: egr           ! Gravtitational energy-generation rate(?)
   real(double) :: r             ! Radius coordinate 
   real(double) :: qq            ! Determines mesh-point metric: mesh-point interval
   real(double) :: qm            ! Determines mesh-point metric: ?  
   
   real(double) :: wl            ! Convective velocity x mixing length
   real(double) :: wcv           ! Convective velocity
   real(double) :: hp            ! Pressure scale height
   real(double) :: wt            !> \todo Something to do with diffusion(?)
   real(double) :: phim          ! Gravitational potential correction for ... rotation?, 
                                 ! used to be called phimu in some routines, change back?
   real(double) :: gmr           ! Effective gravity at the surface(?)
   real(double) :: sep           ! Orbital separation
   real(double) :: m3            ! m^(1/3), m is the mass coordinate
   
   real(double) :: px(npx)       ! px(:) stores the variables of sx(:,k) for the mesh point of the *previous*(?) iteration
   real(double) :: sx(npx,NM+1)  ! sx(:,k) stores a large number of variables for the mesh points k=1:NM
   real(double) :: qa(NM)        ! qa(:) contains H(:,k) for mesh point k
   
end module structure_variables

