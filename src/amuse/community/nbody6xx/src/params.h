*       NBODY6 parameters.
*       ------------------
*
*       Choose between small or large run.
*      PARAMETER  (NMAX=640200,KMAX=5000,LMAX=300,MMAX=10,
*      PARAMETER  (NMAX=262200,KMAX=100,LMAX=300,MMAX=10,
*      PARAMETER  (NMAX=135000,KMAX=100,LMAX=300,MMAX=10,
*     PARAMETER  (NMAX=65566,KMAX=2010,LMAX=500,MMAX=100,
      PARAMETER  (NMAX=40010,KMAX=7010,LMAX=500,MMAX=100,
*     PARAMETER  (NMAX=65566,KMAX=1010,LMAX=256,MMAX=10,
*     PARAMETER  (NMAX=4010,KMAX=1010,LMAX=100,MMAX=10,
     &            MLD=22,MLR=22,MLV=10,MCL=10,NCMAX=10,NTMAX=100)
      parameter (maxpe=1024)
*
*
*       ------------------------------------------------------
*       NMAX    Maximum number of single bodies + 3*NBIN + NHI.
*       KMAX    Maximum number of KS solutions.
*       LMAX    Maximum size of neighbour lists.
*       MMAX    Maximum number of merged binaries.
*       MLD     Maximum number of disrupted KS components.
*       MLR     Maximum number of recently regularized KS pairs.
*       MLV     Maximum number of high-velocity particles.
*       MCL     Maximum number of interstellar clouds.
*       NCMAX   Maximum number of chain members (do not change).
*       NTMAX   Maximum number of circularizing binaries.
*       ------------------------------------------------------
*
