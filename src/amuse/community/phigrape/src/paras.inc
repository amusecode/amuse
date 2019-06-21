!---  general parameters to scale the size of arrays

!---  particles
!          NMAX                maximum number of particles
!          NMAX_LOC            maximum number of particles per processor (usually size of GRAPE memory)
!          NDIM                number of dimensions (must be 3)
!          NBHMAX              maximum number of BH particles      
      INTEGER    NMAX,NDIM,NMAX_LOC,NBHMAX
      PARAMETER (NMAX=2**20,NDIM=3,NMAX_LOC=2**17,NBHMAX=10)

!---  minimum time step
      INTEGER    ndtmin          ! minimum time step is dt_min = 2**mdtmin
      PARAMETER (ndtmin=-30)

!---  maximum number of processors
      INTEGER    NPEMAX
      PARAMETER (NPEMAX=32)
      
