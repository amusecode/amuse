*       PARAMS.H
*       --------
*
*
*       MONT-CAR parameters.
*       --------------------
*
      INTEGER NMAX,NBMAX2,NBMAX3,NBMAX4,NBMAXP,NZONMA,NSUPZO,NLAGRA
*
      PARAMETER (NMAX=700000,NBMAX2=100,NBMAX3=80000,NBMAX4=10,
c      PARAMETER (NMAX=50000,NBMAX2=100,NBMAX3=8000,NBMAX4=10,
c      PARAMETER (NMAX=2500000,NBMAX2=100,NBMAX3=250000,NBMAX4=10,
c     &           NBMAXP=10,NZONMA=10100,NSUPZO=200,NLAGRA=16)
     &           NBMAXP=10,NZONMA=40100,NSUPZO=200,NLAGRA=16)
*
*
*
*       --------------------------------------------------------------
*       NMAX      Maximum number of stars (single stars + 2*binaries + 
*                 safety number)
*       NBMAX2    Maximum number of 2-body binaries.
*       NBMAX3    Maximum number of 3-body binaries.
*       NBMAX4    Maximum number of 4-body binaries.
*       NBMAXP    Maximum number of primorodial binaries.
*       NZONMA    Maximum number of zons.
*       NSUPZO    Maximum number of super zones.
*       NLAGRA    Number of Lagrangian radii
*       --------------------------------------------------------------
*
