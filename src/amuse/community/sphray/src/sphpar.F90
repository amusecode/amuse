!> \file sphpar.F90

!> \brief module for determining smoothing lengths of particles
!<

module sphpar_mod
use myf03_mod
use particle_system_mod, only: particle_type
use particle_system_mod, only: transformation_type

!---------------------------------
!> local particle type. 
  type sphpar_type

     real(r4b)    :: pos(3)     !< x,y,z coordinates

#ifdef incVel
     real(r4b)    :: vel
#endif

     integer(i4b) :: id         !< particle id
     real(r4b)    :: mass       !< particle mass
     real(r4b)    :: T          !< temperature in K       
     real(r4b)    :: rho        !< density 
     real(r4b)    :: ye         !< electron fraction
     real(r4b)    :: xHI        !< nHI/nH 
     real(r4b)    :: xHII       !< * nHII/nH
     real(r4b)    :: hsml       !< smoothing length
     
#ifdef incCloudy
     real(r4b)    :: xHI_cloudy !< cloudy eq solutions
#endif
     
#ifdef incHmf   
     real(r4b)    :: Hmf        !< Hydrogen mass fraction
#endif

#ifdef incHe
     real(r4b)    :: xHeI
     real(r4b)    :: xHeII
     real(r4b)    :: xHeIII
#endif
     
#ifdef incHemf
     real(r4b)    :: Hemf       !< Helium mass fraction
#endif
     
#ifdef outGammaHI
     real(r4b)    :: gammaHI    !< * time averaged HI photoionization rate
     real(r4b)    :: time       !< * elapsed time in seconds - reset at outputs
#endif
     
     integer(i8b) :: lasthit    !< * indx of last ray to cross this particle
     
     integer(i4b) :: nnb     !< number of SPH neighbors
     real(r4b)    :: gradrho(3) !< vector gradient of the density
     real(r4b)    :: drhodh     !< derivative of rho wrt to h
     real(r4b)    :: fi         !< used for finding smoothing length
     real(r4b)    :: dfi        !< used for finding smoothing length

 end type sphpar_type

contains

!-----------------------------------------------------------------------
!> creates a transormed sph particle from an input particle
  subroutine par2sphpar(sphpar,par,pt)
    type(sphpar_type) :: sphpar  !< output sphpar
    type(particle_type) :: par   !< input particle
    type(transformation_type),optional :: pt  !< possible transformation

    sphpar%pos = par%pos

#ifdef incVel
    sphpar%vel = par%vel
#endif

    sphpar%id = par%id
    sphpar%mass = par%mass
    sphpar%T = par%T 
    sphpar%rho = par%rho
    sphpar%ye = par%ye
    sphpar%xHI = par%xHI
    sphpar%xHII = par%xHII
    sphpar%hsml = par%hsml

#ifdef incCloudy
    sphpar%xHI_cloudy = par%xHI_cloudy
#endif

#ifdef incHmf
    sphpar%Hmf = par%Hmf
#endif

#ifdef incHe
    sphpar%xHeI = par%xHeI
    sphpar%xHeII = par%xHeII
    sphpar%xHeIII = par%xHeIII
#endif

#ifdef incHemf
    sphpar%Hemf = par%Hemf
#endif

#ifdef outGammaHI
    sphpar%gammaHI = par%gammaHI
    sphpar%time = par%time
#endif

    sphpar%lasthit = par%lasthit

    if(present(pt)) sphpar%pos=pt%fac*(par%pos-pt%shift)

  end subroutine par2sphpar


end module sphpar_mod
