      MODULE MESH
      !MvdS: Define the maximum number of meshpoints
      INTEGER, PARAMETER :: NM = 500

      ! Number of variables in the PX and SX arrays (output variables)
      INTEGER, PARAMETER :: NPX = 80

      ! Variables for star 1, star 2 and binary parameters
      INTEGER, PARAMETER :: NSVSTAR = 16
      INTEGER, PARAMETER :: NSVBIN = 8
      
      ! Standard number of variables (40)
      INTEGER, PARAMETER :: NSVAR = NSVSTAR + NSVBIN + NSVSTAR
      ! Extra variables
      INTEGER, PARAMETER :: NXVSTAR = 40
      INTEGER, PARAMETER :: NXVBIN = 0
      INTEGER, PARAMETER :: NXVAR = 2*NXVSTAR + NXVBIN

      ! Total number of variables      
      INTEGER, PARAMETER :: NVAR = NSVAR + NXVAR
      INTEGER, PARAMETER :: NVSTAR = NSVSTAR + NXVSTAR
      INTEGER, PARAMETER :: NVBIN = NSVBIN + NXVBIN
      
      ! The arrays holding independent variables can hold up to NVAR numbers.
      ! The first NSVAR of these are treated in the `old' way, which means that
      ! the blocks of 16/8/16 for *1/bin/*2 are retained and respected.
      ! The next NXVAR variables are extra and treated seperately. The block
      ! order in which these are stored is *1/*2/bin
      ! The first extra variable for *1 is at element NSVAR+1 and at
      !  NSVAR+NVXSTAR+1 for *2, while binary parameters start at 
      !  NSVAR+2*NVXSTAR+1
      
      ! Standard number of dependent variables (`functions') (100)
      ! These are 42 for each star and 16 for binary parameters, put down
      ! as *1/bin/*2
      INTEGER, PARAMETER :: NSFSTAR = 42
      INTEGER, PARAMETER :: NSFBIN = 16
      INTEGER, PARAMETER :: NSFUNC = NSFSTAR + NSFBIN + NSFSTAR
      ! Extra functions
      INTEGER, PARAMETER :: NXFSTAR = 30
      INTEGER, PARAMETER :: NXFBIN = 0
      INTEGER, PARAMETER :: NXFUNC = 2*NXFSTAR + NXFBIN

      ! Total number of functions
      INTEGER, PARAMETER :: NFUNC = NSFUNC+NXFUNC

      ! The arrays holding dependent variables can hold up to NFUNC numbers.
      ! The first NSFUNC of these are treated in the `old' way, which means that
      ! the blocks of 42/16/42 for *1/bin/*2 are retained and respected.
      ! The next NXFUNC variables are extra and treated separately. The block
      ! order in which these are stored is *1/*2/bin
      ! The first extra function for *1 is at element NSFUNC+1 and at
      !  NSFUNC+NXFSTAR+1 for *2, while binary parameters start at 
      !  NSFUNC+2*NXFSTAR+1

      ! Number of possible equations
      INTEGER, PARAMETER :: NSESTAR = NSVSTAR
      INTEGER, PARAMETER :: NSEBIN = NSVBIN
      
      ! Standard number of equations (40)
      INTEGER, PARAMETER :: NSEQ = NSESTAR + NSEBIN + NSESTAR
      ! Extra variables
      INTEGER, PARAMETER :: NXESTAR = NXVSTAR
      INTEGER, PARAMETER :: NXEBIN = NXVBIN
      INTEGER, PARAMETER :: NXEQ = 2*NXESTAR + NXEBIN

      ! Total number of equations
      INTEGER, PARAMETER :: NEQ = NSEQ + NXEQ
      
      END MODULE MESH
      
