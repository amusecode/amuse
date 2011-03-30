module mesh
  use real_kind
  
  implicit none
  ! Define the maximum number of meshpoints
  integer, parameter :: NM = 50000
  
  ! Number of variables in the PX and SX arrays (output variables)
  integer, parameter :: npx = 80
  
  ! Variables for star 1, star 2 and binary parameters
  integer, parameter :: nsvstar = 16
  integer, parameter :: nsvbin = 8
  
  ! Standard number of variables (40)
  integer, parameter :: nsvar = nsvstar + nsvbin + nsvstar
  ! Extra variables
  integer, parameter :: nxvstar = 40
  integer, parameter :: nxvbin = 0
  integer, parameter :: nxvar = 2*nxvstar + nxvbin
  
  ! Total number of variables      
  integer, parameter :: nvar = nsvar + nxvar
  integer, parameter :: nvstar = nsvstar + nxvstar
  integer, parameter :: nvbin = nsvbin + nxvbin
  
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
  integer, parameter :: nsfstar = 42
  integer, parameter :: nsfbin = 16
  integer, parameter :: nsfunc = nsfstar + nsfbin + nsfstar
  ! Extra functions
  integer, parameter :: nxfstar = 30
  integer, parameter :: nxfbin = 0
  integer, parameter :: nxfunc = 2*nxfstar + nxfbin
  
  ! Total number of functions
  integer, parameter :: nfunc = nsfunc + nxfunc
  
  ! The arrays holding dependent variables can hold up to NFUNC numbers.
  ! The first NSFUNC of these are treated in the `old' way, which means that
  ! the blocks of 42/16/42 for *1/bin/*2 are retained and respected.
  ! The next NXFUNC variables are extra and treated separately. The block
  ! order in which these are stored is *1/*2/bin
  ! The first extra function for *1 is at element NSFUNC+1 and at
  !  NSFUNC+NXFSTAR+1 for *2, while binary parameters start at 
  !  NSFUNC+2*NXFSTAR+1
  
  ! Number of possible equations
  integer, parameter :: nsestar = nsvstar
  integer, parameter :: nsebin = nsvbin
  
  ! Standard number of equations (40)
  integer, parameter :: nseq = nsestar + nsebin + nsestar
  ! Extra variables
  integer, parameter :: nxestar = nxvstar
  integer, parameter :: nxebin = nxvbin
  integer, parameter :: nxeq = 2*nxestar + nxebin
  
  ! Total number of equations
  integer, parameter :: neq = nseq + nxeq
  
  ! Global variables, formerly in the unnamed COMMON block
  ! TODO: move these to their own proper module
  real(double), save :: h(nvar,NM), dh(nvar,NM), eps, del, dh0
  integer, save :: kh, ktw, isb, id(130), ie(130)
end module mesh

